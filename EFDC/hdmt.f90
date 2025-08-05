  ! ----------------------------------------------------------------------
  !   This file is a part of EFDC+
  !   Website:  https://eemodelingsystem.com/
  !   Repository: https://github.com/dsi-llc/EFDC_Plus.git
  ! ----------------------------------------------------------------------
  ! Copyright 2021-2024 DSI, LLC
  ! Distributed under the GNU GPLv2 License.
  ! ----------------------------------------------------------------------
  SUBROUTINE HDMT

  ! *** SUBROUTINE HDMT EXECUTES THE FULL HYDRODYNAMIC AND MASS TRANSPORT
  ! *** TIME INTEGRATION
  !
  !----------------------------------------------------------------------C
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2014-08-12    D H CHUNG         SET EXPLICIT PRECISIONS OF INTEGER & REAL
  !    2012-09       Chung Dang        OMP
  !    05/01/2002    John Hamrick      Modified calls to calbal and budget subroutines
  !    09-22-2004    Paul M. Craig     Merged DS and TT versions with the 06-04-2004 TT code
  !
  !
  ! *** *******************************************************************C
  !
  use GLOBAL
  use DRIFTER  ,only:DRIFTER_CALC
  use WINDWAVE ,only:WINDWAVEINIT,WINDWAVETUR,WINDWAVECAL
  use HIFREQOUT
  use RESTART_MODULE
  use EFDCOUT
  use CALCSERMOD,only: CALCSER
  use FIELDS
  use WATERQUALITY,only:WQ3D
  use Variables_MPI
  use MPI_All_Reduce
  use OMP_LIB
  use MOD_NETCDF
  use MPI
  use Communicate_Ghost_Routines
  use Mod_Map_Write_EE_Binary
  use Mod_Map_Gather_Sort
  use Mod_Map_Write_NetCDF
  use Variables_Propwash
  use Mod_GOTM

  implicit none
  integer(IK4) :: NS
  integer :: L, ND, NTMP1, NTMP2, NTMP, K, IMAX, JMAX, KMAX
  integer :: IMIN, JMIN, KMIN, NMD, LS, LP
  integer :: ILOGC, NRAMP
  integer :: LN, LNW, LSE, LF, LL, LE, LW

  integer :: M1, M2, LVARY, LLVARY, NVARY, NTSSAVE
  real :: SALMIN, HPPTMP, WTM, WTMP, TMP, AVTMP, ABTMP
  real :: DELVOL, SALMAX, TAUB2, DELTD2, DZDDELT
  real :: TAUBC, TAUBC2, UTMP, VTMP, CURANG, USGZ, VSGZ

  real(RKD), external :: DSTIME

  ! *** 3TL Timestep adjustment variables
  integer :: LVOLMAX, LLVOLMAX, LVRISE, LLVRISE
  real :: CTIM, DEPTHMIN(NQSIJ), PEAKFLOW(NQSER)
  real :: TDIFF, WTM1, WTM2, DTHMIN, PRISEMAX, PRISE
  real :: SAVEDT, SAVETQ
  real(RKD) :: DAYNEXT, VARYNEXT, DDT, DTINC, VUPDATEINC, VQSERINC, VQSERNEXT
  real(RKD) :: TTDS, T1TMP, TIME_RESTART

#ifdef _WIN
  logical, external :: KEY_PRESSED
  logical, external :: ISEXIT
#endif
  real(8)           :: DEL !RESTTIME
  logical           :: RES

  integer :: ierr

  ! *** Perform calculation on a single master process
  if( process_id == master_id )then
    write(*,'(A)')'STARTING HDMT 3TL'
  endif

  T1TMP = DSTIME(0)
  FOURDPI = 4./PI
  RES = .TRUE.

  ! *** INITIALIZE TIME FOR INITIAL CONDITION COMPUTATIONS
  TIMESEC = DBLE(TCON)*DBLE(TBEGIN)
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8
  DAYNEXT = TIMEDAY + 1.0
  DTWQ = 0.0

  ! *** *******************************************************************C
  !
  ! *** INITIALIZE COURNT NUMBER DIAGNOSTICS
  !
  ! *** *******************************************************************!
  ! *** REINITIALIZE VARIABLES
  if( ISRESTI == 0 )then
    do L = 2,LA
      H2P(L) = HP(L)
      H1P(L) = HP(L)
      H1U(L) = HU(L)
      H1UI(L) = HUI(L)
      H1V(L) = HV(L)
      H1VI(L) = HVI(L)
      UHDY1E(L) = UHDYE(L)
      VHDX1E(L) = VHDXE(L)
    enddo

    do K = 1,KC
      do L = 2,LA
        HPK(L,K) = HP(L)*DZC(L,K)
        H2PK(L,K) = HPK(L,K)
        H1PK(L,K) = HPK(L,K)
        U1(L,K) = U(L,K)
        V1(L,K) = V(L,K)
        UCTR1(L,K) = UCTR(L,K)
        VCTR1(L,K) = VCTR(L,K)
        UHDYF1(L,K) = UHDYF(L,K)
        VHDXF1(L,K) = VHDXF(L,K)
        UHDY1(L,K) = UHDY(L,K)
        VHDX1(L,K) = VHDX(L,K)
      enddo
    enddo
  endif

  ! *** *******************************************************************!
  ! *** INITIALIZE COURANT NUMBER DIAGNOSTICS
  if( ISINWV == 1 )then
    do K = 1,KC
      do L = 2,LA
        CFLUUU(L,K) = 0.
        CFLVVV(L,K) = 0.
        CFLWWW(L,K) = 0.
        CFLCAC(L,K) = 0.
      enddo
    enddo
  endif
  ILOGC = 0

  ! *** *******************************************************************C
  ! *** CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
  ! *** CALCULATE VELOCITY GRADIENTS
  !----------------------------------------------------------------------C
  do L = 2,LA
    LN = LNC(L)
    LS = LSC(L)
    LE = LEC(L)
    LW = LWC(L)
    LNW = LNWC(L)
    LSE = LSEC(L)

    !UV(L)  = 0.25*( HP(LS) *(U(LSE,KSZV(L)) +U(LS,KSZV(L))  ) + HP(L) *(U(LE,KSZV(L)) +U(L,KSZV(L))  ) )*HVI(L)
    !U1V(L) = 0.25*( H1P(LS)*(U1(LSE,KSZV(L))+U1(LS,KSZV(L)) ) + H1P(L)*(U1(LE,KSZV(L))+U1(L,KSZV(L)) ) )*H1VI(L)
    !VU(L)  = 0.25*( HP(LW) *(V(LNW,KSZU(L)) +V(LW,KSZU(L))  ) + HP(L) *(V(LN,KSZU(L)) +V(L,KSZU(L))  ) )*HUI(L)
    !V1U(L) = 0.25*( H1P(LW)*(V1(LNW,KSZU(L))+V1(LW,KSZU(L)) ) + H1P(L)*(V1(LN,KSZU(L))+V1(L,KSZU(L)) ) )*H1UI(L)
    if( LKSZ(LS,KSZ(L)) )then
      UV(L) = 0.5*( U(LE,KSZ(L)) + U(L,KSZ(L)) )
      U1V(L) = 0.5*( U1(LE,KSZ(L)) + U1(L,KSZ(L)) )
    else
      UV(L)= 0.25*( HP(LS)*( U(LSE,KSZU(LSE)) + U(LS,KSZU(LS)) ) + HP(L)*( U(LE,KSZU(LE)) + U(L,KSZU(L)) ) )*HVI(L)
      U1V(L)= 0.25*( H1P(LS)*( U1(LSE,KSZU(LSE)) + U1(LS,KSZU(LS)) ) + H1P(L)*( U1(LE,KSZU(LE)) + U1(L,KSZU(L)) ) )*H1VI(L)
    endif
    if( LKSZ(LW,KSZ(L)) )then
      VU(L)= 0.5*( V(LN,KSZ(L)) + V(L,KSZ(L)) )
      V1U(L)= 0.5*( V1(LN,KSZ(L)) + V1(L,KSZ(L)) )
    else
      VU(L)= 0.25*( HP(LW)*( V(LNW,KSZV(LNW)) + V(LW,KSZV(LW)) ) + HP(L)*( V(LN,KSZV(LN)) + V(L,KSZV(L)) ) )*HUI(L)
      V1U(L)= 0.25*( H1P(LW)*( V1(LNW,KSZV(LNW)) + V1(LW,KSZV(LW)) ) + H1P(L)*( V1(LN,KSZV(LN)) + V1(L,KSZV(L)) ) )*H1UI(L)
    endif

  enddo

  ! *** MPI Communication of ghost cells
  call communicate_ghost_cells(UV,  'UV')
  call communicate_ghost_cells(VU,  'VU')

  ! *** *******************************************************************!
  ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8
  if( ISWAVE > 2 .and. LSEDZLJ )then
    call WINDWAVEINIT
    call WINDWAVECAL
  else
    if( ISWAVE == 1 ) CALL WAVEBL
    if( ISWAVE == 2 ) CALL WAVESXY
    if( ISWAVE >= 3 .and. NWSER > 0 )then
      call WINDWAVEINIT
      call WINDWAVETUR
    endif
  endif

  !> @todo make this work with MPI
  ! *** NETCDF INIT
  NC_DATESTAMP = ''
  if( NCDFOUT > 0 )then
    ! *** Call subroutine to gather and map variables in prep for writing to the NetCDF file(s)
    call set_nc_flags
    if( TIMEDAY >= TBEGNCDF ) Call Map_Write_NetCDF

    ! *** Only do NETCDF writing on the master
    if( process_id == master_id )then
      call READCORN
      if( TIMEDAY >= TBEGNCDF ) CALL nc_output()
    endif ! *** end on master

  endif

  if( HFREOUT == 1 )then
    call Gather_High_Frequency
    if( process_id == master_id )then
      do NS = 1,NSUBSET
        if( NCDFOUT > 0 ) call nc_output_hf(1,ns)
        call HFREHYOUT(1,NS)
        call HFREWCOUT(1,NS)
        if( ISTRAN(8) >= 1 )then
          call HFREWQOUT(1,NS)
          if( ISRPEM > 0 ) CALL HFRERPEMOUT(1,NS)
        endif
      enddo
    endif
  endif

  ! *** *******************************************************************!
  ! *** FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  ISTL = 3
  IS2TL = 0
  DTDYN = DT  ! *** PMC - FOR INITIALIZATION

  ! *** *******************************************************************!
  ! *** CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
  if( ISHDMF >= 1 ) CALL CALHDMF3

  ! *** *******************************************************************!
  ! *** CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
  !------------------------------------------------------------------------!
  N = -1
  call CALTSXY(1)

  if( ISRESTI == 0 .or. Restart_In_Ver < 1210 )then
    do L = 2,LA
      TBX1(L) = (AVCON1*H1UI(L) + STBX(L)*SQRT(V1U(L)*V1U(L) + U1(L,KSZ(L))*U1(L,KSZ(L)))) * U1(L,KSZ(L))
      TBY1(L) = (AVCON1*H1VI(L) + STBY(L)*SQRT(U1V(L)*U1V(L) + V1(L,KSZ(L))*V1(L,KSZ(L)))) * V1(L,KSZ(L))
      TSX1(L) = TSX(L)
      TSY1(L) = TSY(L)
    enddo
  endif

  ! *** *******************************************************************!
  ! *** SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  call CALTBXY

  if( ISRESTI == 0 .or. Restart_In_Ver < 1200 )then
    ! *** *******************************************************************C
    ! *** SET BOTTOM AND SURFACE STRESSES
    do L = 2,LA
      TBX1(L) = (AVCON1*H1UI(L) + STBX(L)*SQRT(V1U(L)*V1U(L) + U1(L,KSZ(L))*U1(L,KSZ(L)))) * U1(L,KSZ(L))
      TBY1(L) = (AVCON1*H1VI(L) + STBY(L)*SQRT(U1V(L)*U1V(L) + V1(L,KSZ(L))*V1(L,KSZ(L)))) * V1(L,KSZ(L))
    enddo

    do L = 2,LA
      USGZ = U(L,KSZU(L))
      VSGZ = V(L,KSZV(L))
      TBX(L) = ( STBX(L)*SQRT(VU(L)*VU(L) + USGZ*USGZ) )*USGZ
      TBY(L) = ( STBY(L)*SQRT(UV(L)*UV(L) + VSGZ*VSGZ) )*VSGZ
    enddo

    call CALTSXY(1)
  endif
  N = 0

  ! *** MPI Communication routines
  call communicate_ghost_cells(TBX, 'TBX')
  call communicate_ghost_cells(TBY, 'TBY')
  call communicate_ghost_cells(TSX, 'TSX')
  call communicate_ghost_cells(TSY, 'TSY')
  call communicate_ghost_cells(TBX1, 'TBX1')
  call communicate_ghost_cells(TBY1, 'TBY1')
  call communicate_ghost_cells(TSX1, 'TSX1')
  call communicate_ghost_cells(TSY1, 'TSY1')

  !----------------------------------------------------------------------!
  ! *** SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
  do L = 2,LA
    HDFUFX(L) = 1.
    HDFUFY(L) = 1.
    HDFUF(L) = 1.
  enddo

  ! *** SET CORNER CORRECTIONS
  do L = 2,LA
    WCOREST(L) = 1.
    WCORWST(L) = 1.
    WCORNTH(L) = 1.
    WCORSTH(L) = 1.
  enddo

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  if( (ISRESTI == 0 .or. Restart_In_Ver < 1200) .and. (ISWAVE == 0 .or. LSEDZLJ) )then
    do L = 2,LA
      TVAR3W(L) = TSX1(LEC(L))
      TVAR3S(L) = TSY1(LNC(L))
      TVAR3E(L) = TBX1(LEC(L))
      TVAR3N(L) = TBY1(LNC(L))
    enddo

    do L = 2,LA
      QQ1(L,0)  = 0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX1(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY1(L))**2)
      QQ1(L,KC) = 0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX1(L))**2 + (RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY1(L))**2)
    enddo

    do L = 2,LA
      TVAR3W(L) = TSX(LEC(L))
      TVAR3S(L) = TSY(LNC(L))
      TVAR3E(L) = TBX(LEC(L))
      TVAR3N(L) = TBY(LNC(L))
    enddo

    do L = 2,LA
      QQ(L,0)  = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L))**2  + (RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L))**2 )
      QQ(L,KC) = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L))**2  + (RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L))**2 )
      QQSQR(L,0) = SQRT(QQ(L,0))
    enddo
  endif    ! *** END OF ISWAVE = 0

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  !----------------------------------------------------------------------!
  if( (ISRESTI == 0 .or. Restart_In_Ver < 1200) .and. ISWAVE >= 1 .and. .not. LSEDZLJ )then
    do L = 2,LA
      TVAR3S(L) = TSY1(LNC(L))
      TVAR3W(L) = TSX1(LEC(L))
      TVAR3E(L) = TBX1(LEC(L))
      TVAR3N(L) = TBY1(LNC(L))
    enddo

    do L = 2,LA
      TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
      TAUBC = 0.5*SQRT(TAUBC2)
      CTAUC(L) = TAUBC
      UTMP = 0.5*STCUV(L)*( U1(LEC(L),KSZ(L)) + U1(L,KSZ(L)) )+1.E-12
      VTMP = 0.5*STCUV(L)*( V1(LNC(L),KSZ(L)) + V1(L,KSZ(L)) )
      CURANG = ATAN2(VTMP,UTMP)
      TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
      TAUB2 = MAX(TAUB2,0.)
      QQ1(L,0)  = CTURB2*SQRT(TAUB2)
      QQ1(L,KC) = 0.5*CTURB2*SQRT( (TVAR3W(L)+TSX1(L))**2 + (TVAR3S(L)+TSY1(L))**2 )
    enddo

    do L = 2,LA
      TVAR3S(L) = TSY(LNC(L))
      TVAR3W(L) = TSX(LEC(L))
      TVAR3E(L) = TBX(LEC(L))
      TVAR3N(L) = TBY(LNC(L))
    enddo

    do L = 2,LA
      TAUBC2 = 0.25*( (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2 )
      TAUBC = SQRT(TAUBC2)
      CTAUC(L) = TAUBC
      UTMP = 0.5*STCUV(L)*(U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)))+1.E-12
      VTMP = 0.5*STCUV(L)*(V(LN,KSZV(LNC(L)))     + V(L,KSZV(L)))
      CURANG = ATAN2(VTMP,UTMP)
      TAUB2 =  TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
      TAUB2 =  MAX(TAUB2,0.)
      QQ(L,0)  = CTURB2*SQRT(TAUB2)
      QQ(L,KC) = 0.5*CTURB2*SQRT( (TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2 )
      QQSQR(L,0) = SQRT(QQ(L,0))
    enddo
  endif

  ! *** SET GRAIN STRESS
  if( (ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1) .and. .not. LSEDZLJ )then
    do L = 2,LA
      TAUBSED(L) = QQ(L,0)/CTURB2
      TAUBSND(L) = QQ(L,0)/CTURB2
    enddo
  endif

  ! *** *******************************************************************C
  !
  ! ***  SET SWITCHES FOR THREE TIME LEVEL STEP
  !
  ISTL = 3
  DELT = DT2
  DELTD2 = DT
  ROLD = 0.
  RNEW = 1.

  ! *** *******************************************************************!
  ! *** BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION
  ! *** SET CYCLE COUNTER AND CALL TIMER
  ! *** *******************************************************************
  NTIMER = 0
  N = 0
  TIME_RESTART = TIMEDAY + TIMERST

  ! *** INITIALIZE & RECORD TIME
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8

  if( process_id == master_id )then
    call TIMELOG(N,TIMEDAY,OUTDIR,0._8)
  endif

  ! *** Setup variable DT parameters
  NINCRMT = 1
  VARYTHRESH = DTSSFAC   !0.025
  if( VARYTHRESH > 0.0 )then
    DEPTHMIN = 1e12
    VARYMAX = DT
    VARYMIN = DTSSDHDT
    VARYINC = DTSSDHDT
    DTINC = 0.1_8*VARYINC
    NTSSAVE = NTS
    NTS = NTS*DT/VARYMIN
    DTMIN = VARYMIN
    NRAMP = 0

    VUPDATEINC = 15.       ! *** Time increment to look forward for QSER inflows (min)
    VARYNEXT = TIMEDAY + VUPDATEINC/1440.       ! ***                            (days)
    VQSERINC = 5.          ! *** Check inflows every 5 minutes                   (min)

    ! *** Start simulation with the minimum DT then update DT each day based on the next 24 hours of inflows
    SAVEDT = DT
    DT = VARYMIN
    DTDYN = DT
  endif

  DDT = DT
  NTIMER = 1

  ! *** ************************************************************************
  ! *** ************************************************************************
  ! *** BEGINNING OF THE MAIN TIME ITERATION LOOP FOR TWO TIME LEVEL SOLUTION
  ! *** ************************************************************************
  ! *** ************************************************************************
  N = 0
  TIMESEC = DBLE(TCON)*DBLE(TBEGIN)

  do while( TIMEDAY < TIMEEND )

    ! *** Check inflows to adjust delta T, if needed
    if( VARYTHRESH > 0.0 )then
      if( TIMEDAY >= VARYNEXT .and. NCTBC == 1 )then
        LVARY = 0
        PEAKFLOW = 0.0

        ! *** Check inflows every VQSERINC days
        VQSERNEXT = TIMEDAY + VQSERINC/1440.
        do while( VQSERNEXT < TIMEDAY+VUPDATEINC/1440.+0.0001 )
          do NS = 1,NQSER
            ! *** Compute inflows @ VQSERNEXT
            M2 = MTSQLAST(NS)
            do while( VQSERNEXT > TSFL(NS).TIM(M2) )
              M2 = M2 + 1
              if( M2 > TSFL(NS).NREC )then
                M2 = TSFL(NS).NREC
                EXIT
              endif
            enddo
            if( M2 > TSFL(NS).NREC ) CYCLE

            M1 = M2-1
            TDIFF = TSFL(NS).TIM(M2)-TSFL(NS).TIM(M1)
            WTM1 = (TSFL(NS).TIM(M2)-VQSERNEXT)/TDIFF
            WTM2 = (VQSERNEXT-TSFL(NS).TIM(M1))/TDIFF
            DELVOL = 0.
            do K = 1,KC
              DELVOL = DELVOL + WTM1*TSFL(NS).VAL(M1,K) + WTM2*TSFL(NS).VAL(M2,K)
            enddo

            if( DELVOL > PEAKFLOW(NS) )then
              PEAKFLOW(NS) = DELVOL                                  ! *** m3/s
            endif
          enddo
          VQSERNEXT = VQSERNEXT + VQSERINC/1440.
        enddo

        ! *** Compute metrics for each cell to determine limiting cell and inflow
        PRISEMAX = 0.0
        do LL = 1,NQSIJ
          NS = BCFL(LL).NQSERQ
          if( NS < 1 ) CYCLE
          L = BCFL(LL).L
          DTHMIN = MAX(DEPTHMIN(LL),0.5*HDRY)

          PRISE = PEAKFLOW(NS)/DTHMIN/DXYP(L)         ! *** 1/s = (m3/s) / (m) / (m2)
          if( PRISE > PRISEMAX )then
            LVARY  = L                  ! *** Cell Index
            LLVARY = LL                 ! *** QSER Index
            PRISEMAX = PRISE            ! *** 1/s inflow loading rate
          endif
        enddo

        ! *** Adjust DT, if needed
        TDIFF = VARYTHRESH / PRISEMAX
        DTDYN = TDIFF
        DTDYN = MIN(DTDYN,VARYMAX)
        DTDYN = MAX(DTDYN,VARYMIN)

        if( ABS(DT-DTDYN) >= 0.999*DTINC )then
          ! *** Report the change in delta t
          !write(20,'(A,I10,2I5,f15.5,F10.3,I6,E12.4,2f8.2)') 'DTA ',N, ISTL, NCTBC, TIMEDAY, DEPTHMIN(LLVARY), LVARY, PRISEMAX, dt, dtdyn   ! *** DELME
          PRINT    '(A,I10,2I5,f15.5,F10.3,I6,E12.4,2f8.2)', 'DTA ',N, ISTL, NCTBC, TIMEDAY, DEPTHMIN(LLVARY), LVARY, PRISEMAX, dt, dtdyn

          open(9,FILE = OUTDIR//'TIME.LOG',POSITION = 'APPEND')
          write(9,'(A,I10,F14.5,A,F10.3,F10.4,F9.4,A,F10.3)') 'DYNAMIC TIMESETP ADJUSTMENT ',N,TIMEDAY,'     [DT,L,HP]   NEW: ', DTDYN, PRISEMAX, DEPTHMIN(LLVARY),'     OLD: ',DT
          close(9)
        endif

        ! *** SET INCREMENTAL INCREASE IN COUNTER
        NINCRMT = INT(DTDYN/VARYINC)

        VARYNEXT = TIMEDAY + VUPDATEINC/1440.
        DEPTHMIN = 1e12
      endif

      ! *** Transition DT if not equal to DTDYN
      if( NCTBC == 1 )then
        if( ABS(DT-DTDYN) >= 0.999*DTINC )then
          if( DDT > DTDYN )then
            DDT = DTDYN            ! *** Ramp down immediately
          else
            if( NRAMP > 20 )then
              DDT = DDT + DTINC    ! *** Slowly ramp up
              NRAMP = 0
            endif
          endif
        endif
        DT      = DDT
        DT2     = 2.*DT
        DELT    = DT2
        DELTD2  = DT
        NRAMP = NRAMP + 1
      endif

      ! *** Check minimum depths
      do LL = 1,NQSIJ
        NS = BCFL(LL).NQSERQ
        L = BCFL(LL).L
        DEPTHMIN(LL) = MIN(DEPTHMIN(LL), HP(L))
      enddo

    endif     ! *** End of variable delta T

    TIMESEC  = TIMESEC  + DDT
    TIMEDAY  = TIMESEC/86400._8

    N = N + NINCRMT
    NITER = N

    ILOGC = ILOGC+1

    ! ***************************************************
    ! *** STARTUP TERMS - MOMENTUM
    if( N <= NLTS )then
      SNLT = 0.
    elseif( N > NLTS .and. N <= NTTS )then
      NTMP1 = N-NLTS
      NTMP2 = NTTS-NLTS+1
      SNLT = FLOAT(NTMP1)/FLOAT(NTMP2)
    else
      SNLT = 1.
    endif

    ! *** STARTUP TERMS - GRAVITY/DENSITY EFFECTS
    if( N <= NTSVB )then
      GP = GPO*(FLOAT(N)/FLOAT(NTSVB))
    else
      GP = GPO
    endif

    ! 2018-11-22, NTL: UPDATE TIME VARIABLE TOPOGRAPHY
    if( BATHY.IFLAG > 0 ) CALL UPDATETOPO(TIMEDAY,BELV,HP,HDRY)
    if( ROUGH.IFLAG > 0 ) CALL UPDATEFIELD(ROUGH,TIMEDAY,1,ZBR)

    !----------------------------------------------------------------------C

    !----------------------------------------------------------------------C
    ! *** REENTER HERE FOR TWO TIME LEVEL CORRECTION
500 continue


    ! *** *******************************************************************!
    ! *** UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS
    call CALQVS

    ! *** *******************************************************************!
    ! *** CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
    if( KC > 1 )then
      TTDS = DSTIME(0)
      if( ISGOTM > 0 )then

        ! *** Update AV, AB computed from GOTM
        call Advance_GOTM(ISTL)

        ! *** Apply maximum value if required
        if( ISAVBMX >= 1 )then
          !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,KS,LLWET,LKWET,AVMX,ABMX,HPI,AV,AB) PRIVATE(ND,K,LP,L,AVTMP,ABTMP)
          do ND = 1,NDM

            do K = 1,KS
              do LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)
                AVTMP = AVMX*HPI(L)
                ABTMP = ABMX*HPI(L)
                AV(L,K) = MIN(AV(L,K),AVTMP)
                AB(L,K) = MIN(AB(L,K),ABTMP)
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO
        endif

        ! ****************************************************************************
        !call MPI_barrier(MPI_Comm_World, ierr)
        !call communicate_ghost_3d0(AV)
        ! ****************************************************************************

        ! *** Compute the inverse of the average AV at U and V interfaces
        if( IGRIDV == 0 )then
          !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,KS,LWC,LSC,LLWET,LKWET,HPI,AV,AVUI,AVVI) PRIVATE(ND,K,LP,L)
          do ND = 1,NDM
            do K = 1,KS
              do LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)
                AVUI(L,K) = 2.0/( AV(L,K) + AV(LWC(L),K) )
                AVVI(L,K) = 2.0/( AV(L,K) + AV(LSC(L),K) )
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO
        else
          !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,KS,LWC,LSC,LLWET,LKWET,SUB,SVB,SVB3D,SUB3D,HP,HU,HV,AV,AVUI,AVVI)    &
          !$OMP             PRIVATE(ND,K,LP,L)
          do ND = 1,NDM
            do K = 1,KS
              do LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)
                AVUI(L,K) = SUB(L)*(1. + SUB3D(L,K))/( AV(L,K)*HP(L) + SUB3D(L,K)*AV(LWC(L),K)*HP(LWC(L) ) )*HU(L)
                AVVI(L,K) = SVB(L)*(1. + SVB3D(L,K))/( AV(L,K)*HP(L) + SVB3D(L,K)*AV(LSC(L),K)*HP(LSC(L) ) )*HV(L)
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO
        endif
      else
        ! *** ORIGINAL EFDC MELLOR-YAMADA TURBULENCE CLOSURE
        call CALAVB
      endif
      TAVB = TAVB + (DSTIME(0)-TTDS)
    endif

    ! *** *******************************************************************!
    ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
    if( ISTL == 3 )then
      TTDS = DSTIME(0)
      if( ISWAVE > 2 .and. LSEDZLJ )then
        call WINDWAVECAL
      else
        if( ISWAVE == 1 ) CALL WAVEBL
        if( ISWAVE == 2 ) CALL WAVESXY
        if( ISWAVE >= 3 .and. NWSER > 0 ) CALL WINDWAVETUR
      endif
      TTBXY = TTBXY + (DSTIME(0)-TTDS)
    endif

    ! *** *******************************************************************!
    ! *** UPDATE TIME VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
    call CALCSER
    call CALVEGSER

    PSERT(0) = 0.
    if( NPSER >= 1 ) CALL CALPSER

    ! *** *******************************************************************!
    ! *** CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
    ! *** CALEXP   -  PRODUCTION VERSION WITH HORIZONTAL MOMENTUM SOURCE
    ! ***             AND 3D IMPLICIT VEGETATION DRAG
    TTDS = DSTIME(0)
    call CALEXP
    TCEXP = TCEXP + (DSTIME(0)-TTDS)

    ! *** *******************************************************************!
    ! *** SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
    TTDS = DSTIME(0)
    call CALPUV9C
    TPUV = TPUV + (DSTIME(0)-TTDS) - TMPITMP - TTWAIT
    DSITIMING(2) = DSITIMING(2)    + TMPITMP

    ! *** *******************************************************************C
    ! *** ADVANCE TIME SNAPSHOT INTERNAL VARIABLES FOR THREE TIME LEVEL STEP
    if( ISTL == 3 )then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,K,L)
      do ND = 1,NDM
        LF = 2 + (ND-1)*LDM
        LL = MIN(LF+LDM-1,LA)

        do K = 1,KC
          do L = LF,LL
            UHDYF2(L,K) = UHDYF1(L,K)   ! *** FLOW BY LAYER FOR THE ENTIRE DEPTH
            UHDYF1(L,K) = UHDYF(L,K)
            VHDXF2(L,K) = VHDXF1(L,K)
            VHDXF1(L,K) = VHDXF(L,K)

            UHDY2(L,K) = UHDY1(L,K)     ! *** FLOW BY LAYER
            UHDY1(L,K) = UHDY(L,K)
            VHDX2(L,K) = VHDX1(L,K)
            VHDX1(L,K) = VHDX(L,K)

            U2(L,K) = U1(L,K)           ! *** N-2 TIME, NOT INTERVAL AVERAGE
            V2(L,K) = V1(L,K)           ! *** N-2 TIME, NOT INTERVAL AVERAGE
            U1(L,K) = U(L,K)
            V1(L,K) = V(L,K)
            W2(L,K) = W1(L,K)           ! *** N-2 TIME, NOT INTERVAL AVERAGE
            W1(L,K) = W(L,K)
            UCTR2(L,K) = UCTR1(L,K)     ! *** N-2 TIME, NOT INTERVAL AVERAGE   DELME - WHAT IS THIS?
            VCTR2(L,K) = VCTR1(L,K)     ! *** N-2 TIME, NOT INTERVAL AVERAGE
            UCTR1(L,K) = UCTR(L,K)
            VCTR1(L,K) = VCTR(L,K)
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif

    ! *** *******************************************************************C
    ! *** ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
    ! *** MODE FORCING
    !----------------------------------------------------------------------C
    if( ISTL == 3 )then
      TTDS = DSTIME(0)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          TSX1(L) = TSX(L)
          TSY1(L) = TSY(L)
        enddo
      enddo
      !$OMP END PARALLEL DO

      call CALTSXY(0)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          DU(L,KS) = DU(L,KS) - CDZUU(L,KS)*TSX(L)
          DV(L,KS) = DV(L,KS) - CDZUV(L,KS)*TSY(L)
        enddo
      enddo
      !$OMP END PARALLEL DO

      TTBXY = TTBXY+(DSTIME(0)-TTDS)
    endif

    ! *** *******************************************************************C
    ! *** SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
    ! *** USING THE INTERNAL SHEARS STORED IN DU & DV
    TTDS = DSTIME(0)
    TMPITMP = 0.
    if( KC > 1 )then
      call CALUVW
    else
      !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LA,LDM,UHDYE,UHDYF,UHDY,HUI,DYIU,U,VHDXE,VHDXF,VHDX,HVI,DXIV,V,W)   &
      !$OMP                           PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM
        LF = 2 + (ND-1)*LDM
        LL = MIN(LF+LDM-1,LA)

        do L = LF,LL
          UHDYF(L,1) = UHDYE(L)
          UHDY(L,1)  = UHDYE(L)
          U(L,1)     = UHDYE(L)*HUI(L)*DYIU(L)
          VHDXF(L,1) = VHDXE(L)
          VHDX(L,1)  = VHDXE(L)
          V(L,1)     = VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)     = 0.
        enddo
      enddo
      !$OMP END PARALLEL DO
      call CALUVW
    endif

    ! *** Compute the cell center of velocity components
    if( ISGOTM > 0 )then
      do ND = 1,NDM
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            UCTR(L,K) = 0.5*(U(L,K) + U(LEC(L),K))
            VCTR(L,K) = 0.5*(V(L,K) + V(LNC(L),K))
          enddo
        enddo
      enddo
    endif

    TUVW = TUVW + (DSTIME(0)-TTDS) - TMPITMP
    DSITIMING(3) = DSITIMING(3) + TMPITMP

    ! *** *******************************************************************C
    ! *** CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
    ! *** AT TIME LEVEL (N+1)
    if( ISTRANACTIVE > 0 .or. ISGOTM > 0 )then
      call CALCONC
    else
      ! *** Call Propwash module if that option is selected
      if( propwash_on .and. ISTL == 3 )then
        call Propwash_Calc_Sequence(0)
      endif
    endif

    ! *** Communicate variables with MPI calls before writing to files and looping computations
    call MPI_barrier(MPI_Comm_World, ierr)
    TTDS = DSTIME(0)
    call Communicate_CON1(0)

    TMPITMP = DSTIME(0) - TTDS
    DSITIMING(6) = DSITIMING(6) + TMPITMP

    !----------------------------------------------------------------------C
    ! *** CHECK RANGE OF SALINITY AND DYE CONCENTRATION
    if( ISMMC == 1 .and. DEBUG )then

      SALMAX = -100000.
      SALMIN = 100000.
      do K = 1,KC
        do L = 2,LA
          if( SAL(L,K) > SALMAX )then
            SALMAX = SAL(L,K)
            IMAX = IL(L)
            JMAX = JL(L)
            KMAX = K
          endif
          if( SAL(L,K) < SALMIN )then
            SALMIN = SAL(L,K)
            IMIN = IL(L)
            JMIN = JL(L)
            KMIN = K
          endif
        enddo
      enddo

      write(6,6001)N
      write(6,6002)SALMAX,IMAX,JMAX,KMAX
      write(6,6003)SALMIN,IMIN,JMIN,KMIN

      SALMAX = -100000.
      SALMIN = 100000.
      do K = 1,KC
        do L = 2,LA
          if( DYE(L,K,1) > SALMAX )then
            SALMAX = DYE(L,K,1)
            IMAX = IL(L)
            JMAX = JL(L)
            KMAX = K
          endif
          if( DYE(L,K,1) < SALMIN )then
            SALMIN = DYE(L,K,1)
            IMIN = IL(L)
            JMIN = JL(L)
            KMIN = K
          endif
        enddo
      enddo

      write(6,6004)SALMAX,IMAX,JMAX,KMAX
      write(6,6005)SALMIN,IMIN,JMIN,KMIN

      SALMAX = -100000.
      SALMIN = 100000.
      do K = 1,KC
        do L = 2,LA
          if( SFL(L,K) > SALMAX )then
            SALMAX = SFL(L,K)
            IMAX = IL(L)
            JMAX = JL(L)
            KMAX = K
          endif
          if( SFL(L,K) < SALMIN )then
            SALMIN = SFL(L,K)
            IMIN = IL(L)
            JMIN = JL(L)
            KMIN = K
          endif
        enddo
      enddo

      write(6,6006)SALMAX,IMAX,JMAX,KMAX
      write(6,6007)SALMIN,IMIN,JMIN,KMIN

    endif

    if( ISMMC == 2 .and. DEBUG )then
      SALMAX = -100000.
      SALMIN = 100000.
      do K = 1,KC
        do L = 2,LA
          if( TEM(L,K) > SALMAX )then
            SALMAX = TEM(L,K)
            IMAX = IL(L)
            JMAX = JL(L)
            KMAX = K
          endif
          if( TEM(L,K) < SALMIN )then
            SALMIN = TEM(L,K)
            IMIN = IL(L)
            JMIN = JL(L)
            KMIN = K
          endif
        enddo
      enddo

      write(6,6001)N
      write(6,6008)SALMAX,IMAX,JMAX,KMAX
      write(6,6009)SALMIN,IMIN,JMIN,KMIN
    endif

6001 FORMAT('  N = ',I10)
6002 FORMAT('  SALMAX = ',F14.4,5X,'I,J,K = ',(3I10))
6003 FORMAT('  SALMIN = ',F14.4,5X,'I,J,K = ',(3I10))
6004 FORMAT('  DYEMAX = ',F14.4,5X,'I,J,K = ',(3I10))
6005 FORMAT('  DYEMIN = ',F14.4,5X,'I,J,K = ',(3I10))
6006 FORMAT('  SFLMAX = ',F14.4,5X,'I,J,K = ',(3I10))
6007 FORMAT('  SFLMIN = ',F14.4,5X,'I,J,K = ',(3I10))
6008 FORMAT('  TEMMAX = ',F14.4,5X,'I,J,K = ',(3I10))
6009 FORMAT('  TEMMIN = ',F14.4,5X,'I,J,K = ',(3I10))

    ! *** *******************************************************************C
    !
    ! *** CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
    ! *** CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOUBLE TIME
    ! *** STEP TRANSPORT FIELD
    !
    !----------------------------------------------------------------------C
    if( ISWQFLUX == 1 )then
      NTMP = MOD(NCTBC,2)
      if( ISTL == 3 ) DTWQ = DTWQ + DT
      if( NTMP == 0 .and. ISTL == 3 )then

        ! *** CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM
          LF = (ND-1)*LDMWET+1
          LL = MIN(LF+LDMWET-1,LAWET)
          do LP = LF,LL
            L = LWET(LP)
            HWQ(L) = 0.25*(H2P(L) + 2.*H1P(L) + HP(L))
          enddo
        enddo
        !$OMP END DO

        if( ISICM >= 1 )then
          !$OMP DO PRIVATE(ND,LF,LL,LP,L)
          do ND = 1,NDM
            LF = (ND-1)*LDMWET+1
            LL = MIN(LF+LDMWET-1,LAWET)
            do LP = LF,LL
              L = LWET(LP)
              TVAR3E(L) = UHDY2E(LEC(L))
              TVAR3N(L) = VHDX2E(LNC(L))
            enddo
          enddo
          !$OMP END DO

          !$OMP DO PRIVATE(ND,K,LP,L)
          do ND = 1,NDM
            do K = 1,KC
              do LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)
                TVAR2E(L,K) = UHDY2(LEC(L),K)
                TVAR2N(L,K) = VHDX2(LNC(L),K)
              enddo
            enddo
          enddo
          !$OMP END DO

          !$OMP DO PRIVATE(ND,LF,LL,LP,L,HPPTMP)
          do ND = 1,NDM
            LF = (ND-1)*LDMWET+1
            LL = MIN(LF+LDMWET-1,LAWET)
            do LP = LF,LL
              L = LWET(LP)
              HPPTMP = H2WQ(L) + DT2*DXYIP(L)*( QSUME(L) - (TVAR3E(L)-UHDY2E(L)+TVAR3N(L)-VHDX2E(L)) )
              HWQ(L) = SPB(L)*HPPTMP + (1.-SPB(L))*HWQ(L)
            enddo
          enddo
          !$OMP END DO
        endif
        !$OMP END PARALLEL

        ! *** ADD CHANNEL INTERACTIONS
        if( MDCHH >= 1 )then
          do NMD = 1,MDCHH
            if( MDCHTYP(NMD) == 1 )then
              HWQ(LMDCHH(NMD)) = HWQ(LMDCHH(NMD))  +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
              HWQ(LMDCHU(NMD)) = HWQ(LMDCHU(NMD))  -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
            endif
            if( MDCHTYP(NMD) == 2 )then
              HWQ(LMDCHH(NMD)) = HWQ(LMDCHH(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
              HWQ(LMDCHV(NMD)) = HWQ(LMDCHV(NMD)) - DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
            endif
            if( MDCHTYP(NMD) == 3 )then
              HWQ(LMDCHH(NMD)) = HWQ(LMDCHH(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
              HWQ(LMDCHU(NMD)) = HWQ(LMDCHU(NMD)) - DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
              HWQ(LMDCHV(NMD)) = HWQ(LMDCHV(NMD)) - DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
            endif
          enddo
        endif      ! *** END ADD CHANNEL INTERACTIONS

        ! *** CALL WATER QAULITY KINETICS, SEDIMENT AND RPEM PROCESSES
        if( ISTRAN(8) >= 1 ) CALL WQ3D
        if( ISTRAN(4) >= 1 ) CALL CALSFT

        if( ISICM >= 1 )then
          do L = 2,LA
            H2WQ(L) = HWQ(L)
          enddo
        endif

        ! *** Communicate variables with MPI calls before writing to files and looping computations
        if( ISTRAN(8) >= 1 )then
          call MPI_barrier(MPI_Comm_World, ierr)
          TTDS = DSTIME(0)
          call communicate_ghost_cells(WQV, NWQV)
          TMPITMP = DSTIME(0) - TTDS
          DSITIMING(6) = DSITIMING(6) + TMPITMP

          DTWQ = 0.0
        endif
      endif

    endif         ! *** END OF WQ SECTION

    ! *** *******************************************************************C
    ! *** UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING AN EQUATION OF STATE
    if( BSC > 1.E-6 )then
      if( ISTL == 3 )then
        call CALBUOY(.TRUE.)
      else
        call CALBUOY(.FALSE.)
      endif
    endif

    ! *** *******************************************************************C
    !
    ! *** CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,LN,LS,LNW,LSE,LE,LW)            &
    !$OMP             SHARED(NDM,LDMWET,LAWET,LWET,LEC,LNC,LSC,LWC,LNWC,LSEC,LSWC,LKSZ)   &
    !$OMP             SHARED(SUB,SVB,UV,VU,HP,U,V,U1V,V1U,HVI,HUI,KSZ,KSZU,KSZV,HPW,HPE,HPN,HPS)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)
      do LP = LF,LL
        L = LWET(LP)
        LN = LNC(L)
        LS = LSC(L)
        LE = LEC(L)
        LW = LWC(L)
        LNW = LNWC(L)
        LSE = LSEC(L)

        U1V(L) = UV(L)
        V1U(L) = VU(L)

        if( LKSZ(LS,KSZ(L)) )then
          UV(L) = 0.5*( U(LE,KSZ(L)) + U(L,KSZ(L)) )
        else
          UV(L)= 0.25*( HP(LS)*( U(LSE,KSZU(LSE)) + U(LS,KSZU(LS)) ) + HP(L)*( U(LE,KSZU(LE)) + U(L,KSZU(L)) ) )*HVI(L)
        endif
        if( LKSZ(LW,KSZ(L)) )then
          VU(L)= 0.5*( V(LN,KSZ(L)) + V(L,KSZ(L)) )
        else
          VU(L)= 0.25*( HP(LW)*( V(LNW,KSZV(LNW)) + V(LW,KSZV(LW)) ) + HP(L)*( V(LN,KSZV(LN)) + V(L,KSZV(L)) ) )*HUI(L)
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! *** *******************************************************************C
    ! *** CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
    ! *** AT TIME LEVEL (N)
    if( ISTL /= 2 .and. ISHDMF >= 1 )then
      TTDS = DSTIME(0)
      call CALHDMF3
      THMDF = THMDF+(DSTIME(0)-TTDS)
    endif

    ! *** *******************************************************************C
    !
    ! *** UPDATE BOTTOM STRESSES AND SURFACE AND BOTTOM TURBULENT
    ! *** INTENSITIES
    if( ISTL == 3 )then
      if( ISCDMA == 2 )then
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM
          LF = (ND-1)*LDMWET+1
          LL = MIN(LF+LDMWET-1,LAWET)
          do LP = LF,LL
            L = LWET(LP)
            TBX1(L)   = TBX(L)
            TBY1(L)   = TBY(L)
            QQ2(L,0)  = QQ(L,0)  + QQ1(L,0)
            QQ2(L,KC) = QQ(L,KC) + QQ1(L,KC)
            QQ1(L,0)  = QQ(L,0)
            QQ1(L,KC) = QQ(L,KC)
          enddo
        enddo
        !$OMP END PARALLEL DO

      else  ! *** if( ISCDMA < 2 )then

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM
          LF = (ND-1)*LDMWET+1
          LL = MIN(LF+LDMWET-1,LAWET)
          do LP = LF,LL
            L = LWET(LP)
            TBX1(L)   = TBX(L)
            TBY1(L)   = TBY(L)
            QQ2(L,0)  = QQ1(L,0)  + QQ1(L,0)
            QQ2(L,KC) = QQ1(L,KC) + QQ1(L,KC)
            QQ1(L,0)  = QQ(L,0)
            QQ1(L,KC) = QQ(L,KC)
          enddo
        enddo
        !$OMP END PARALLEL DO

      endif
    endif

    ! *** *******************************************************************C
    !
    ! *** CALCULATE BOTTOM STRESS AT LEVEL (N+1)

    TTDS = DSTIME(0)
    call CALTBXY

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,AVCON1,USGZ,VSGZ)       &
    !$OMP             SHARED(ICALTB,ISAVCOMP,NDM,LDMWET,LAWET,LWET,KSZ,KSZU,KSZV) &
    !$OMP             SHARED(HUI,HVI,TBX,STBX,VU,U,TBY,STBY,UV,V)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)

      if( ICALTB > 0 .and. ISAVCOMP == 0 )then
        do LP = LF,LL
          L = LWET(LP)
          USGZ   = U(L,KSZU(L))
          VSGZ   = V(L,KSZV(L))
          TBX(L) = AVCON1*HUI(L)*USGZ
          TBY(L) = AVCON1*HVI(L)*VSGZ
        enddo
      else
        do LP = LF,LL
          L = LWET(LP)
          USGZ   = U(L,KSZU(L))
          VSGZ   = V(L,KSZV(L))
          TBX(L) = ( STBX(L)*SQRT(VU(L)*VU(L) + USGZ*USGZ) )*USGZ
          TBY(L) = ( STBY(L)*SQRT(UV(L)*UV(L) + VSGZ*VSGZ) )*VSGZ
        enddo
      endif
    enddo
    !$OMP END PARALLEL DO

    TTBXY = TTBXY + (DSTIME(0)-TTDS)

    ! *** *******************************************************************C
    !
    ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
    !
    !----------------------------------------------------------------------C
    TTDS = DSTIME(0)
    if( ISWAVE == 0 .or. LSEDZLJ )then

      ! ***  STANDARD CALCULATIONS - NO CORNER CORRECTS
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,L,K,LF,LL,LP,TMP)           &
      !$OMP             SHARED(NDM,LDMWET,LAWET,LWET,KC,LEC,LNC)             &
      !$OMP             SHARED(TVAR3S,TVAR3W,TVAR3E,TVAR3N,TBX,TBY,TSX,TSY)  &
      !$OMP             SHARED(CTURB2,RSSBCE,RSSBCW,RSSBCN,RSSBCS,QQ)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)

        do LP = LF,LL
          L = LWET(LP)
          TVAR3W(L) = TSX(LEC(L))
          TVAR3S(L) = TSY(LNC(L))
          TVAR3E(L) = TBX(LEC(L))
          TVAR3N(L) = TBY(LNC(L))
        enddo

        do LP = LF,LL
          L = LWET(LP)
          TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2
          QQ(L,0)  = 0.5*CTURB2*SQRT(TMP)

          TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2
          QQ(L,KC) = 0.5*CTURB2*SQRT(TMP)
        enddo
      enddo
      !$OMP END PARALLEL DO

    endif

    ! *** *******************************************************************C
    !
    ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
    !
    !----------------------------------------------------------------------C
    if( ISWAVE >= 1 .and. .not. LSEDZLJ )then

      !$OMP PARALLEL DEFAULT(SHARED)
      !$OMP DO PRIVATE(ND,LF,LL,LP,L)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          TVAR3W(L) = TSX(LEC(L))
          TVAR3S(L) = TSY(LNC(L))
          TVAR3E(L) = TBX(LEC(L))
          TVAR3N(L) = TBY(LNC(L))
        enddo
      enddo
      !$OMP END DO

      !$OMP DO PRIVATE(ND,LF,LL,LP,L,TAUBC,TAUBC2,UTMP,VTMP,CURANG,TAUB2,TMP)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          if( LWVMASK(L) )then
            TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
            TAUBC = 0.5*SQRT(TAUBC2)   ! *** CURRENT ONLY
            CTAUC(L) = TAUBC
            UTMP = 0.5*STCUV(L)*( U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)) )+1.E-12
            VTMP = 0.5*STCUV(L)*( V(LN ,KSZV(LN) ) + V(L,KSZV(L)) )
            CURANG = ATAN2(VTMP,UTMP)
            TAUB2 = TAUBC*TAUBC +      (QQWV3(L)*QQWV3(L)) + 2.     *TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
            TAUB2 = MAX(TAUB2,0.)              ! *** CURRENT & WAVE
            QQ(L,0 )   = CTURB2*SQRT(TAUB2)  ! *** CELL CENTERED TURBULENT INTENSITY DUE TO CURRENTS & WAVES

            QQ(L,KC)   = 0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2)
          else
            TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2
            QQ(L,0 ) = 0.5*CTURB2*SQRT(TMP)

            TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2
            QQ(L,KC) = 0.5*CTURB2*SQRT(TMP)
          endif
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif

    ! *** *******************************************************************C
    !
    ! *** CALCULATE TURBULENT INTENSITY SQUARED
    !
    if( KC > 1 )then
      if( ISGOTM > 0 )then
        !$OMP PARALLEL DO DEFAULT(NONE) SHARED(LA,KC,QQ,QQSQR,TKE3D,DML,GL3D,HP) PRIVATE(L)
        do L = 2,LA
          QQ(L,0:KC)  = 2.*tke3d(L,0:KC)
          QQSQR(L,0)  = SQRT(QQ(L,0))
          DML(L,0:KC) = GL3D(L,0:KC)/HP(L)
        enddo
        !$OMP END PARALLEL DO
      else
        call CALQQ1
      endif
    endif
    TQQQ = TQQQ + (DSTIME(0)-TTDS)

    ! ****************************************************************************
    ! *** MPI communication
    call MPI_barrier(MPI_Comm_World, ierr)
    TTDS = DSTIME(0)
    if( KC > 1 )then
      call Communicate_QQ
    else
      call communicate_ghost_3d0(QQ)
      call Communicate_1D2(TBX, TBY)
      call Communicate_1D2(UV, VU)
    endif
    DSITIMING(7) = DSITIMING(7) + (DSTIME(0)-TTDS)
    ! ****************************************************************************
    ! *** COMPUTE THE SQRT OF THE TURBULENCE (M/S)
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,KS,LDMWET,LAWET,LWET,LLWET,LKWET,QQ,QQSQR) PRIVATE(ND,K,LF,LL,LP,L)
    do ND = 1,NDM  
      do K = 0,KS
        if( K == 0 )then
          LF = (ND-1)*LDMWET+1
          LL = MIN(LF+LDMWET-1,LAWET)
          do LP = LF,LL
            L = LWET(LP)
            QQSQR(L,0) = SQRT(QQ(L,0))  
          enddo
        else
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            QQSQR(L,K) = SQRT(QQ(L,K))  
          enddo
        endif
      enddo
    enddo 
    !$OMP END PARALLEL DO
    ! *** *******************************************************************C
    ! *** CALCULATE MEAN MASS TRANSPORT FIELD
    if( (ISSSMMT > 0 .or. ISWASP > 0) .and. RESSTEP > 0. )then
      NTMP = MOD(NCTBC,2)
      if( ISTL == 3 .and. NTMP == 0 ) CALL CALMMT
    endif

    ! *** *******************************************************************C
    !
    ! *** HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
    ! *** IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
    !
    !----------------------------------------------------------------------C
    if( NCTBC == NTSTBC )then
      NCTBC = 0
      ISTL = 2
      DELT = DT          ! *** NOT USED IN HDMT BUT INDICATIVE
      DELTD2 = 0.5*DT
      ROLD = 0.5
      RNEW = 0.5
      GOTO 500
    else
      NCTBC = NCTBC+1
      ISTL = 3
      DELT = DT2         ! *** NOT USED IN HDMT BUT INDICATIVE
      DELTD2 = DT
      ROLD = 0.
      RNEW = 1.
    endif

    ! *** *******************************************************************!
    ! *** *******************************************************************!
    ! *** HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
    ! *** *******************************************************************!
    ! *** *******************************************************************!

    ! *** WRITE TO TIME SERIES FILES
    CTIM = TIMESEC/TCON

    if( ISTMSR >= 1 )then
      if( NCTMSR >= NWTMSR )then
        call TMSR
        NCTMSR = 1
      else
        NCTMSR = NCTMSR+1
      endif
    endif

    ! *** *******************************************************************C
    ! *** CALCULATE MEAN MASS TRANSPORT FIELD
    if( (ISSSMMT > 0 .or. ISWASP > 0) .and. RESSTEP > 0 )then
      if( ISICM == 0 ) CALL CALMMT
    endif

    ! *** *******************************************************************C
    !
    ! *** ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
    !
    if( ISPD > 0 )then
      if( TIMEDAY >= LA_BEGTI0 .and. TIMEDAY <= LA_ENDTI0 )then
        TTDS = DSTIME(0)
        call DRIFTER_CALC
        TLRPD = TLRPD + (DSTIME(0)-TTDS)
      endif
    endif

    ! *** *******************************************************************!
    ! *** WRITE TO TIME VARYING GRAPHICS FILES
    if( HFREOUT == 1 )then
      call Gather_High_Frequency
      if( process_id == master_id )then
        do NS = 1,NSUBSET
          DEL = ABS(84600*(TIMEDAY-HFREDAY(NS)))
          if( DEL <= DELT .and. TIMEDAY <= HFREDAYEN(NS) )then
            if( NCDFOUT > 0 ) call nc_output_hf(0,ns)
            call HFREHYOUT(0,NS)
            call HFREWCOUT(0,NS)
            if( ISTRAN(8) >= 1 )then
              call HFREWQOUT(0,NS)
              if( ISRPEM > 0 ) CALL HFRERPEMOUT(0,NS)
            endif
            HFREDAY(NS) = HFREDAY(NS) + HFREMIN(NS)/1440
          endif
        enddo
      endif
    endif

    ! *** NETCDF OUTPUT
    if( NCDFOUT > 0 )then
      if( (TIMEDAY >= TBEGNCDF + NCDFSHOT*NCFREQ/1440) .and. (TIMEDAY <= TENDNCDF)  )then
        call Map_Write_NetCDF

        if( process_id == master_id )then
          call nc_output()
        endif
        NCDFSHOT = NCDFSHOT + 1.
      endif
    endif
    ! *** *******************************************************************C
    ! *** WRITE EFDC EXPLORER FORMAT OUTPUT
    if( ISPPH == 1 )then
      DTDYN = DT
      ! *** CHECK IF PROPWASH OUTPUT FORCES EE LINKAGE TO BE WRITTEN
      if( ISPROPWASH > 0 )then
        ! *** Determine if PW_Mesh was writting on any process
        if( iwrite_pwmesh > 0 )then
          if( TIMEDAY >= last_snapshot + freq_out_min/1440. )then
            iwrite_pwmesh = 2
          else
            iwrite_pwmesh = 0
          endif
        endif
      endif

      !CALL BCOUT_ACCUMULATE(DTDYN) @todo invoke and make work in MPI case
      if( ISTL == 3 .and. TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )then
        TTDS = DSTIME(0)
        call Map_Write_EE_Binary
        TMPIEE = TMPIEE + (DSTIME(0)-TTDS)

        if( process_id == master_id )then
          call EE_LINKAGE(0)
        endif
        if(iwrite_pwmesh > 0) last_snapshot = TIMEDAY
      endif
      if( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )then
        NSNAPSHOTS = NSNAPSHOTS + 1
      endif
    endif

    ! *** *******************************************************************!
    ! *** WRITE RESTART FILE EVERY ISRESTO REFERENCE PERIODS
    if( ABS(ISRESTO) >= 1 )then
      if( TIMEDAY >= TIME_RESTART )then
        call Restart_Out(TIME_RESTART, 0)
        if( ISTRAN(8) >= 1 )then
          call WQ_WCRST_OUT(TIME_RESTART)
          if( IWQBEN == 1 ) CALL WQSDRST_OUT(TIME_RESTART)
          if( ISRPEM > 0  ) CALL WQRPEMRST_Out(TIME_RESTART)
        endif
        TIME_RESTART = TIMEDAY + TIMERST
      endif
    endif

    ! *** *******************************************************************!
    ! *** RECORD TIME
    if( process_id == master_id )then
      if( TIMEDAY >= DAYNEXT )then
        TTDS = DSTIME(0) - T1TMP
        call TIMELOG(N,TIMEDAY,OUTDIR, TTDS)

        DAYNEXT = DAYNEXT + 1.
      endif
    endif


    ! *** *******************************************************************!
    if( process_id == master_id )then
      if( ISHOW > 0 ) CALL SHOWVAL
    endif
    ! *** *******************************************************************!

#ifdef _WIN  
    if( num_Processors == 1 )then     ! *** OMP Run
      if( KEY_PRESSED() )then
        if( ISEXIT() )then
          TTDS = DSTIME(0)
          call Map_Write_EE_Binary
          TMPIEE = TMPIEE + (DSTIME(0)-TTDS)

          call EE_LINKAGE(-1)

          EXIT                      ! *** Break out of time loop
        endif
      endif
    endif
#endif

    ! *** Update time counters
    if( TIMEDAY > HOURNEXT )then
      HOURNEXT = HOURNEXT + 1./24.
    endif
    if( TIMEDAY > HOUR06NEXT )then
      HOUR06NEXT = HOUR06NEXT + 6./24.
    endif
    if( TIMEDAY > HOUR12NEXT )then
      HOUR12NEXT = HOUR12NEXT + 12./24.
    endif

  enddo
  ! *** *******************************************************************!
  ! *** *******************************************************************!
  ! *** TIME LOOP COMPLETED
  ! *** *******************************************************************!
  ! *** *******************************************************************!

  THDMT = THDMT + (DSTIME(0)-T1TMP)    ! *** CPU time

  ! *** *******************************************************************C
  ! *** WRITE RESTART FILE
  if( ABS(ISRESTO) > 0 .and. TIMEDAY > TIME_RESTART-TIMERST*0.5 )then
    if( process_id == master_id ) WRITE(6,'(A,F12.4)') 'FINAL RESTART FILE: ',TIMEDAY
    call Restart_Out(TIME_RESTART,0)
    if( ISTRAN(8) >= 1 )then
      call WQ_WCRST_OUT(TIME_RESTART)
      if( IWQBEN == 1 ) CALL WQSDRST_OUT(TIME_RESTART)
      if( ISRPEM > 0  ) CALL WQRPEMRST_Out(TIME_RESTART)
    endif
  endif

  if( process_id == master_id )then
    ! *** *******************************************************************C
    ! *** OUTPUT COSMETIC VOLUME LOSSES FORM DRY CELLS
    if( NDRYSTP > 0 )then
      open(1,FILE = OUTDIR//'DRYLOSS.OUT')
      close(1,STATUS = 'DELETE')
      open(1,FILE = OUTDIR//'DRYLOSS.OUT')

      write(1,'(3A5,A14,A10)')'L','I','J','VOL WASTED','EQUIV HP'
      do L = 2,LA
        write(1,1993)L,IL(L),JL(L),VDWASTE(L),VDWASTE(L)/DXYP(L)
      enddo
1993  FORMAT(3I5,E14.6,F10.3)

      close(1)
    endif

    ! *** *******************************************************************C
    ! *** OUTPUT COURANT NUMBER DIAGNOSTICS
    if( ISINWV == 1 )then
      open(1,FILE = OUTDIR//'CFLMAX.OUT')
      close(1,STATUS = 'DELETE')
      open(1,FILE = OUTDIR//'CFLMAX.OUT')

      do L = 2,LA
        write(1,1991)IL(L),JL(L),(CFLUUU(L,K),K = 1,KC)
        write(1,1992)(CFLVVV(L,K),K = 1,KC)
        write(1,1992)(CFLWWW(L,K),K = 1,KC)
        write(1,1992)(CFLCAC(L,K),K = 1,KC)
      enddo

      close(1)
1991  FORMAT(2I5,12F7.2)
1992  FORMAT(10X,12F7.2)
    endif

    ! *** *******************************************************************C
    !
    ! *** OUTPUT FINAL FOOD CHAIN AVERAGING PERIOD
    if( ISTRAN(5) >= 1 .and. ISFDCH >= 1 ) CALL FOODCHAIN(1)

  endif ! End calculation on single processor

  END
