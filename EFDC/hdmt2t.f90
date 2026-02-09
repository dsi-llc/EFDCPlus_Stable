! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details EXECUTES THE FULL HYDRODYNAMIC AND MASS
!! TRANSPORT TIME INTEGRATION USING A TWO TIME LEVEL SCHEME
!
! @date    2015-06       Paul M. Craig     Implemented Sigma-Z (SGZ) in EE7.3
!!         2014-08-12    D H Chung         Set explicit precisions of integer & real
!!         2011-03       Paul M. Craig     Rewritten to F90 and added OMP
!!         2010-01       Chung Dang        Added the DSI version of Lagrangian Particle Tracking
!!         09-22-2004    Paul M. Craig     Merged DS and TT versions with the 06-04-2004 TT code
!!         05/02/2002    John Hamrick      Modified calculation of cell center bed stress (stored as QQ(l,0))
!!                                         for cells have source/sinks



SUBROUTINE HDMT2T

  ! *** *******************************************************************!
  use GLOBAL
  use Allocate_Initialize      
  use DRIFTER  ,only:DRIFTER_CALC
  use WINDWAVE ,only:WINDWAVEINIT,WINDWAVECAL,WINDWAVETUR
  use HIFREQOUT
  use RESTART_MODULE
  use EFDCOUT
  use CALCSERMOD,only: CALCSER
  use FIELDS
#ifndef GNU  
  USE IFPORT
#endif
  use WATERQUALITY,only:WQ3D
  use Variables_WQ
  use Variables_Propwash

  use Variables_MPI
  use OMP_LIB
  use MOD_NETCDF

  use MPI
  use Communicate_Ghost_Routines
  use Mod_Map_Write_EE_Binary
  use Mod_Map_Gather_Sort
  use Mod_Map_Write_NetCDF
  use MPI_All_Reduce

  use Mod_GOTM

  implicit none

  integer :: ILOGC, K, L, LP, LE, LW, LN, LS, LNW, LSE, LSW, ND, LF, LL, LG
  integer :: NLOOP, NTMP1, NTMP2, NMD, NDIFF, NS, NT

  real      :: TMP, HDFUFM, TAUBC2, TAUBC, UTMP, VTMP, CURANG, TAUB2, DELTD2, DZDDELT
  real      :: CTIM, WTM, WTMP, DELVOL, USGZ, VSGZ, AVTMP, ABTMP
  real(8)   :: DEL
  logical   :: LDEBUG
  character :: LFILENAME*30

  real,save,allocatable,dimension(:) :: WCOREW
  real,save,allocatable,dimension(:) :: WCORNS
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, T1TMP, TWAIT                 ! MODEL TIMING TEMPORARY VARIABLES
  real(RKD)           :: TIME_RESTART, DAYNEXT
  
  integer,save,allocatable,dimension(:) :: LCORNER
  integer,save,allocatable,dimension(:) :: LCORNWE
  integer,save,allocatable,dimension(:) :: LCORNSN

#ifdef _WIN
  logical, external :: KEY_PRESSED
  logical, external :: ISEXIT
#endif
  logical   :: RES, BSYNC
  real,save :: ADJC
  
  ! *** New variables for MPI
  integer :: ierr, i, j
  integer :: Start_Local_LA
  integer :: End_Local_LA

  real, dimension(2) :: DYN_IN, DYN_OUT
  
  logical :: lpmc = .false.
  
  ! *** ALLOCATE LOCAL ARRAYS
  if( .not. allocated(WCOREW) )then
    call AllocateDSI( WCOREW,   LCM,    0.0)
    call AllocateDSI( WCORNS,   LCM,    0.0)
    call AllocateDSI( LCORNER,  LCM,      0)
    call AllocateDSI( LCORNWE,  LCM,      0)
    call AllocateDSI( LCORNSN,  LCM,      0)

    if( KC > 1 )then
      ADJC = 1.0
    else
      ADJC = 1.1
    endif
    TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8
    DAYNEXT = TIMEDAY + 1.0
    TIME_RESTART = TIMEDAY + TIMERST
  endif

  if( process_id == master_id )then
    write(*,'(A)')'STARTING HDMT 2TL'
  endif
  T1TMP  = DSTIME(0)

  ISTL = 2
  IS2TL = 1
  FOURDPI = 4./PI
  RES = .TRUE.
  LDEBUG = .FALSE.

  ! *** INITIALIZE TIME FOR INITIAL CONDITION COMPUTATIONS
  TIMESEC = DBLE(TCON)*DBLE(TBEGIN)

  ! *** *******************************************************************!
  ! *** SET FLAGS FOR CORNER CELL BED STRESS CORRECTIONS
  if( ISCORTBC >= 1 )then
    ! *** SET FLAG FOR CELLS HAVING VOLUME SOURCE OR SINKS
    do L = 1,LC
      ISSBCP(L) = 0
    enddo

    do L = 2,LA
      if( RSSBCE(L) > 1.5)ISSBCP(L) = 1
      if( RSSBCW(L) > 1.5)ISSBCP(L) = 1
      if( RSSBCN(L) > 1.5)ISSBCP(L) = 1
      if( RSSBCS(L) > 1.5)ISSBCP(L) = 1
    enddo
  endif

  do L = 2,LA
    WCOREST(L) = 1.
    WCORWST(L) = 1.
    WCORNTH(L) = 1.
    WCORSTH(L) = 1.
  enddo

  ! *** *******************************************************************!
  ! *** REINITIALIZE VARIABLES
  if( ISRESTI == 0 )then
    do L = 2,LA
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
        H1PK(L,K) = HPK(L,K)
        U1(L,K) = U(L,K)
        V1(L,K) = V(L,K)
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


  ! *** *******************************************************************!
  ! *** CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
  ! *** CALCULATE VELOCITY GRADIENTS
  !----------------------------------------------------------------------!
  do L = 2,LA
    LSE = LSEC(L)
    LS  = LSC(L)
    LE  = LEC(L)
    LNW = LNWC(L)
    LW  = LWC(L)
    LN  = LNC(L)

    !UV(L)  = 0.25*( HP(LS  )*(U(LSE, KSZV(L))+U(LS,  KSZV(L)) ) + HP(L )*(U(LE, KSZV(L))+U(L, KSZV(L)) ) )*HVI(L)
    !U1V(L) = 0.25*( H1P(LS )*(U1(LSE,KSZV(L))+U1(LS, KSZV(L)) ) + H1P(L)*(U1(LE,KSZV(L))+U1(L,KSZV(L)) ) )*H1VI(L)
    !VU(L)  = 0.25*( HP(LWC(L) )*(V(LNW, KSZU(L))+V(LWC(L), KSZU(L)) ) + HP(L )*(V(LN,  KSZU(L))+V(L, KSZU(L)) ) )*HUI(L)
    !V1U(L) = 0.25*( H1P(LWC(L))*(V1(LNW,KSZU(L))+V1(LWC(L),KSZU(L)) ) + H1P(L)*(V1(LN, KSZU(L))+V1(L,KSZU(L)) ) )*H1UI(L)

    UV(L) = SVB(L)*( 0.5*HP(LS)*(U(LSE,KSZU(LSE)) + U(LS,KSZU(LS))) + 0.5*HP(L)*(U(LE,KSZU(LE)) + U(L,KSZU(L))) )*HVI(L)
    U1V(L) = SVB(L)*( 0.5*H1P(LS)*(U1(LSE,KSZU(LSE)) + U1(LS,KSZU(LS))) + 0.5*H1P(L)*(U1(LE,KSZU(LE)) + U1(L,KSZU(L))) )*H1VI(L)
    VU(L) = SUB(L)*( 0.5*HP(LW)*(V(LNW,KSZV(LNW)) + V(LW,KSZV(LW))) + 0.5*HP(L)*(V(LN,KSZV(LN)) + V(L,KSZV(L))) )*HUI(L)
    V1U(L) = SUB(L)*( 0.5*H1P(LW)*(V1(LNW,KSZV(LNW)) + V1(LW,KSZV(LW))) + 0.5*H1P(L)*(V1(LN,KSZV(LN)) + V1(L,KSZV(L))) )*H1UI(L)
  enddo

  ! *** MPI Communication of ghost cells
  ! ****************************************************************************
  call communicate_ghost_cells(UV,  'UV')
  call communicate_ghost_cells(U1V, 'U1V')
  call communicate_ghost_cells(VU,  'VU')
  call communicate_ghost_cells(V1U, 'V1U')
  ! ****************************************************************************

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
  DTDYN = DT  ! *** PMC - FOR INITIALIZATION
  if( ISRESTI == 0 .or. Restart_In_Ver < 1200 ) CALL CALTBXY

  ! *** *******************************************************************!
  ! *** CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
  if( ISHDMF >= 1 )then
    call CALHDMF
  endif

  ! *** *******************************************************************!
  ! *** CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
  !----------------------------------------------------------------------!
  N = -1
  call CALTSXY(1)

  ! *** *******************************************************************!
  ! *** SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  call CALTBXY

  if( ISRESTI == 0 .or. Restart_In_Ver < 1200 )then
    ! *** SET BOTTOM AND SURFACE STRESSES
    do L = 2,LA
      USGZ = U(L,KSZU(L))
      VSGZ = V(L,KSZV(L))
      TBX(L) = ( STBX(L)*SQRT(VU(L)*VU(L) + USGZ*USGZ) )*USGZ
      TBY(L) = ( STBY(L)*SQRT(UV(L)*UV(L) + VSGZ*VSGZ) )*VSGZ
    enddo
  endif
  
  ! *** MPI Communication routines
  ! ****************************************************************************
  call communicate_ghost_cells(TBX, 'TBX')
  call communicate_ghost_cells(TBY, 'TBY')
  ! ****************************************************************************

  N = 0

  !CALL CALTSXY

  !----------------------------------------------------------------------!
  ! *** SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
  do L = 2,LA
    HDFUFX(L) = 1.
    HDFUFY(L) = 1.
    HDFUF(L) = 1.
  enddo

  ! ***
  if( ISBSDFUF >= 1 )then
    HDFUFM = 1.E-12

    do L = 2,LA
      LS = LSC(L)
      HDFUFX(L) = HDFUFM+G*SUB(L)*HU(L)*(BELV(LWC(L))-BELV(L))*DXIU(L)
      HDFUFY(L) = HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)
    enddo

    do L = 2,LA
      if( HDFUFX(L)>0.0 )then
        HDFUFX(L) = TBX(L)/HDFUFX(L)
      else
        HDFUFX(L) = 1.0
      endif
      if( HDFUFY(L)>0.0 )then
        HDFUFY(L) = TBY(L)/HDFUFY(L)
      else
        HDFUFY(L) = 1.0
      endif
    enddo

  endif

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  if( ISWAVE == 0 .or. LSEDZLJ )then

    if( ISCORTBC == 0 )then
      ! *** NO CORNER CORRECTIONS - STANDARD APPROACH
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          TVAR3S(L) = TSY(LNC(L))
          TVAR3W(L) = TSX(LEC(L))
          TVAR3E(L) = TBX(LEC(L))
          TVAR3N(L) = TBY(LNC(L))
        enddo
      enddo

      ! *** COMPUTE CELL CENTERED BOTTOM AND SURFACE
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          QQ(L,0)  = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2 )
          QQ(L,KC) = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2 + (RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2 )
        enddo
      enddo

    else                             !IF(ISCORTBC >= 1 )then

      ! *** WITH CORNER CORRECTIONS
      do L = 2,LA
        if( ISSBCP(L) == 0 )then
          if( SUB(LEC(L)) < 0.5) WCOREST(L) = FSCORTBCV(L)
          if( SUB(L) < 0.5) WCORWST(L) = FSCORTBCV(L)
          if( SVB(LNC(L)) < 0.5) WCORNTH(L) = FSCORTBCV(L)
          if( SVB(L) < 0.5) WCORSTH(L) = FSCORTBCV(L)
        endif
      enddo

      do L = 2,LA
        WCOREW(L) = 1./(WCOREST(L)+WCORWST(L))
        WCORNS(L) = 1./(WCORNTH(L)+WCORSTH(L))
      enddo

      do L = 2,LA
        WCOREST(L) = WCOREST(L)*WCOREW(L)
        WCORWST(L) = WCORWST(L)*WCOREW(L)
        WCORNTH(L) = WCORNTH(L)*WCORNS(L)
        WCORSTH(L) = WCORSTH(L)*WCORNS(L)
      enddo

      do L = 2,LA
        TVAR3S(L) = TSY(LNC(L))
        TVAR3W(L) = TSX(LEC(L))
        TVAR3E(L) = TBX(LEC(L)   )
        TVAR3N(L) = TBY(LNC(L))
      enddo

      do L = 2,LA
        QQ(L,0 ) = CTURB2*SQRT(  (RSSBCE(L)*WCOREST(L)*TVAR3E(L) + RSSBCW(L)*WCORWST(L)*TBX(L))**2  + (RSSBCN(L)*WCORNTH(L)*TVAR3N(L) + RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
        QQ(L,KC) = 0.5*CTURB2*SQRT(  (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2  +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
      enddo

    endif
  endif    ! *** END OF ISWAVE = 0

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  !----------------------------------------------------------------------!
  if( ISWAVE >= 1 .and. .not. LSEDZLJ )then

    do L = 2,LA
      TVAR3S(L) = TSY(LNC(L))
      TVAR3W(L) = TSX(LEC(L))
      TVAR3E(L) = TBX(LEC(L)   )
      TVAR3N(L) = TBY(LNC(L))
    enddo

    do L = 2,LA
      TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
      TAUBC = 0.5*SQRT(TAUBC2)   ! *** CURRENT ONLY
      CTAUC(L) = TAUBC
      UTMP = 0.5*STCUV(L)*(U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)))+1.E-12
      VTMP = 0.5*STCUV(L)*(V(LN ,KSZV(LN )) + V(L,KSZV(L)))
      CURANG = ATAN2(VTMP,UTMP)
      TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
      TAUB2 = max(TAUB2,0.)          ! *** CURRENT & WAVE
      QQ(L,0 ) = CTURB2*SQRT(TAUB2)  ! *** CELL CENTERED TURBULENT INTENSITY DUE TO CURRENTS & WAVES
      QQ(L,KC) = 0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2)
    enddo
  endif
  QQSQR = SQRT(QQ)

  ! *** SET GRAIN STRESS
  if( (ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1) .and. .not. LSEDZLJ )then
    do L = 2,LA
      TAUBSED(L) = QQ(L,0)/CTURB2
      TAUBSND(L) = QQ(L,0)/CTURB2
    enddo
  endif

  ! *** *******************************************************************!
  ! ***  SET SWITCHES FOR TWO TIME LEVEL INTEGRATION
  DELT = DT
  DELTD2 = DT/2.

  ! *** *******************************************************************!
  ! *** BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION
  ! *** SET CYCLE COUNTER AND CALL TIMER
  NTIMER = 0
  N = 0

  ! *** EFDC_EXPLORER BEGIN BLOCK  RECORD TIME
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8

  if( process_id == master_id )then
    call TIMELOG(N, TIMEDAY, OUTDIR ,0._8)
  endif

  NTIMER = 1
  NINCRMT = 1
  NITER = 0
  NLOOP = 0

  ! *** ************************************************************************
  ! *** ************************************************************************
  ! *** BEGINNING OF THE MAIN TIME ITERATION LOOP FOR TWO TIME LEVEL SOLUTION
  ! *** ************************************************************************
  ! *** ************************************************************************
1001 continue
  if( TIMEDAY > TIMEEND ) GO TO 1000

  if( ISDYNSTP == 0 )then
    N = N + 1
  else
    NLOOP = NLOOP + 1
    if( NLOOP > NRAMPUP )then
      call CALSTEPD
    else
      DTDYN = DT
      NINCRMT = 1
    endif
    DELT    = DTDYN
    DELTD2  = DTDYN/2.
    N = N + NINCRMT
  endif
  NITER = NITER + 1

  TIMESEC = DBLE(DT)*DBLE(N)+DBLE(TCON)*DBLE(TBEGIN)
  TIMEDAY = TIMESEC/86400._8
  
  !print '(A,5I5,2F15.5)', ' HDMT 00 ', niter, n, process_id, NINCRMT, NSNAPSHOTS, SNAPSHOTS(NSNAPSHOTS), TIMEDAY   ! delme
  !print *, 'HDMT 12', niter, process_id   ! delme
  !if( NINCRMT > 15 )then
  !  pause 'pause: '  ! delme
  !endif
  
  
  if( ISDYNSTP == 0 )then
    ILOGC = ILOGC + 1
  else
    ILOGC = ILOGC + NINCRMT
  endif

  ! *** DSI BEGIN BLOCK
  if( N <= NLTS )then
    SNLT = 0.
  elseif( N > NLTS .and. N <= NTTS )then
    NTMP1 = N - NLTS
    NTMP2 = NTTS - NLTS + 1
    SNLT = FLOAT(NTMP1)/FLOAT(NTMP2)
  else
    SNLT = 1.
  endif

  ! *** SET THE GRAVITY ACCELERATION
  if( N <= NTSVB )then
    GP = GPO*(FLOAT(N)/FLOAT(NTSVB))
  else
    GP = GPO
  endif

  ! 2018-11-22, NTL: UPDATE TIME VARIABLE TOPOGRAPHY
  if( BATHY.IFLAG > 0 ) CALL UPDATETOPO(TIMEDAY,BELV,HP,HDRY)
  if( ROUGH.IFLAG > 0 ) CALL UPDATEFIELD(ROUGH,TIMEDAY,1,ZBR)

  ! *** *******************************************************************!
  ! *** UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS
  if( lpmc) print *, NIter, 'CALQVS'  ! delme
  call CALQVS

  ! *** *******************************************************************!
  ! *** CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
  if( KC > 1 )then
    TTDS = DSTIME(0)
    if( ISGOTM > 0 )then
      ! *** Update AQ with AV computed from GOTM - DTK
      call Advance_GOTM(ISTL)
                
      ! *** Apply maximum value if required
      if( ISAVBMX >= 1 )then  
        do ND = 1,NDM  
        
          do K = 1,KS  
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              AVTMP = AVMX*HPI(L)  
              ABTMP = ABMX*HPI(L)  
              AV(L,K) = min(AV(L,K),AVTMP)  
              AB(L,K) = min(AB(L,K),ABTMP)  
            enddo  
          enddo  
        enddo
      endif
      
      ! *** Compute the inverse of the average AV at U and V interfaces
      if( IGRIDV == 0 )then
        do ND = 1,NDM  
          do K = 1,KS  
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              AVUI(L,K) = 2.0/( AV(L,K) + AV(LWC(L),K) )
              AVVI(L,K) = 2.0/( AV(L,K) + AV(LSC(L),K) )
            enddo  
          enddo 
        enddo
      else
        do ND = 1,NDM  
          do K = 1,KS  
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              AVUI(L,K) = SUB(L)*(1. + SUB3D(L,K))/( AV(L,K)*HP(L) + SUB3D(L,K)*AV(LWC(L),K)*HP(LWC(L) ) )*HU(L)
              AVVI(L,K) = SVB(L)*(1. + SVB3D(L,K))/( AV(L,K)*HP(L) + SVB3D(L,K)*AV(LSC(L),K)*HP(LSC(L) ) )*HV(L)
            enddo  
          enddo 
        enddo
      endif
    else
      ! *** ORIGINAL EFDC MELLOR-YAMADA TURBULENCE CLOSURE
      call CALAVB
    endif
    TAVB = TAVB + (DSTIME(0)-TTDS)
  endif

  ! *** *******************************************************************!
  ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
  TTDS = DSTIME(0)
  if( ISWAVE > 2 .and. LSEDZLJ )then
    call WINDWAVECAL
  else
    if( ISWAVE == 1 ) CALL WAVEBL
    if( ISWAVE == 2 ) CALL WAVESXY
    if( ISWAVE >= 3 .and. NWSER > 0 ) CALL WINDWAVETUR
  endif
  TTBXY = TTBXY + (DSTIME(0)-TTDS)

  ! *** *******************************************************************!
  ! *** UPDATE TIME VARIABLE SURFACE WIND STRESSES
  TTDS = DSTIME(0)
  TMPITMP = 0.
  if( lpmc) print *, NIter, 'CALTSXY'  ! delme
  call CALTSXY(0)
  TTBXY = TTBXY + (DSTIME(0)-TTDS) - TMPITMP
  DSITIMING(9) = DSITIMING(9) + TMPITMP

  ! *** *******************************************************************!
  ! *** UPDATE TIME VARIABLE VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
  if( lpmc) print *, NIter, 'CALCSER'  ! delme
  call CALCSER
  call CALVEGSER
  
  PSERT(0) = 0.
  if( lpmc) print *, NIter, 'CALPSER'  ! delme
  if( NPSER >= 1 ) CALL CALPSER

  ! *** *******************************************************************!
  ! *** CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
  TTDS = DSTIME(0)
  if( lpmc) print *, NIter, 'CALEXP2T'  ! delme
  call CALEXP2T
  TCEXP = TCEXP + (DSTIME(0)-TTDS)

  ! *** *******************************************************************!
  ! *** SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
  TTDS = DSTIME(0)
  TMPITMP = 0.
  if( lpmc) print *, NIter, 'CALPUV2C'  ! delme
  call CALPUV2C 
  TPUV = TPUV + (DSTIME(0)-TTDS) - TMPITMP - TTWAIT          ! *** Includes CONGRAD time (TTWAIT is the wait time)
  DSITIMING(2) = DSITIMING(2)    + TMPITMP

  ! *** *******************************************************************!
  ! *** Advance internal variables
  UHDYF1 = UHDYF
  VHDXF1 = VHDXF
  UHDY1  = UHDY
  VHDX1  = VHDX
  U1     = U
  V1     = V
  W1     = W
  UCTR1  = UCTR
  VCTR1  = VCTR
  
  ! *** *******************************************************************!
  ! *** SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
  ! *** USING THE INTERNAL SHEARS STORED IN DU & DV
  TTDS = DSTIME(0)
  TMPITMP = 0.
  if( KC > 1 )then
    if( lpmc) print *, NIter, 'CALUVW'  ! delme
    call CALUVW 
  else
    UHDYF(:,1) = UHDYE(:)
    VHDXF(:,1) = VHDXE(:)
    UHDY(:,1)  = UHDYE(:)
    VHDX(:,1)  = VHDXE(:)
    U(:,1)     = UHDYE(:)*HUI(:)*DYIU(:)
    V(:,1)     = VHDXE(:)*HVI(:)*DXIV(:)
    W(:,1)     = 0.

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

  ! *** *******************************************************************!
  ! *** CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
  ! *** AT TIME LEVEL (iN+1)
  if( ISTRANACTIVE > 0 .or. ISGOTM > 0 )then
    if( lpmc)  print *, NIter, 'CALCONC'  ! delme
    call CALCONC
  else
    ! *** Call Propwash module if that option is selected
    if( propwash_on )then
      call Propwash_Calc_Sequence(0)
    endif
  endif

  ! *** Communicate variables with MPI calls before writing to files and looping computations
  ! ****************************************************************************
  call MPI_barrier(DSIcomm, ierr)
  TTDS = DSTIME(0)
  TMPITMP = 0.
  call Communicate_CON1
  TMPITMP = TMPITMP + DSTIME(0) - TTDS
  DSITIMING(6) = DSITIMING(6) + TMPITMP
  ! ****************************************************************************

  ! *** *******************************************************************!
  
  ! *** CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
  ! *** CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOUBLE TIME
  ! *** STEP TRANSPORT FIELD
  if( ISWQFLUX == 1 )then
    DTWQ = DELT
    if( ISICM >= 1 .or. ISRPEM > 0 .or. ISTRAN(4) >= 1 )then
      do L = 2,LA
        HWQ(L) = HP(L)
      enddo

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
    endif
    
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
      call MPI_barrier(DSIcomm, ierr)
      TTDS = DSTIME(0)
      call communicate_ghost_cells(WQV, NWQV)
      TMPITMP = DSTIME(0) - TTDS
      DSITIMING(6) = DSITIMING(6) + TMPITMP
    endif
  endif         ! *** END OF WQ SECTION

  ! *** *******************************************************************!
  ! ***  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING AN EQUATION OF STATE
  if( BSC > 1.E-6  )then
    if( lpmc) print *, NIter, 'CALBUOY'  ! delme
    call CALBUOY(.TRUE.)
  endif

  ! *** *******************************************************************!
  ! *** CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)
  if( ICALTB == 0 )then
    ! *** EFDC+ 10.0 and later approach to compute the shear velocities at the U and V faces

    ! *** Advance UV and VU
    U1V = UV
    V1U = VU
    
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)
      do LP = LF,LL
        L = LWET(LP)
        LSE = LSEC(L)
        LS  = LSC(L)
        LE  = LEC(L)
        LNW = LNWC(L)
        LW  = LWC(L)
        LN  = LNC(L)

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
  else
    ! *** LEGACY EFDC APPLICATIONS
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)

      do LP = LF,LL
        L = LWET(LP)
        LSE = LSEC(L)
        LS = LSC(L)
        LE = LEC(L)
        LNW = LNWC(L)
        LW = LWC(L)
        LN = LNC(L)

        UV(L) = ( ( HP(LSE)*U(LSE,KSZU(LSE)) + HP(LS)*U(LS,KSZU(LS)) ) + ( HP(LE)*U(LE,KSZU(LE)) + HP(L)*U(L,KSZU(L)) ) )/( 1.0 + SUB(LE) + SUB(LS) + SUB(LSE) )*HVI(L)
        VU(L) = ( ( HP(LNW)*V(LNW,KSZV(LNW)) + HP(LW)*V(LW,KSZV(LW)) ) + ( HP(LN)*V(LN,KSZV(LN)) + HP(L)*V(L,KSZV(L)) ) )/( 1.0 + SVB(LN) + SVB(LW) + SVB(LNW) )*HUI(L)
      enddo
    enddo

  endif

  ! *** *******************************************************************!
  ! *** CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
  ! *** AT TIME LEVEL (N)
  if( ISHDMF >= 1 )then
    TTDS = DSTIME(0)

    if( lpmc) print *, NIter, 'CALHDMF'  ! delme
    call CALHDMF
    
    THMDF = THMDF + (DSTIME(0)-TTDS)
  endif

  ! *** *******************************************************************!
  ! *** CALCULATE BOTTOM STRESS AT LEVEL (N+1)
  TTDS = DSTIME(0)
  call CALTBXY

  ! *** Set U and V at the upstream interface layer
  do L = 2,LA
    TVAR3W(L) = U(L,KSZU(L))
    TVAR3S(L) = V(L,KSZV(L))
  enddo
  
  if( ICALTB > 0 .and. ISAVCOMP == 0 )then
    TBX(:) = AVCON1*HUI(:)*TVAR3W(:)
    TBY(:) = AVCON1*HVI(:)*TVAR3S(:)
  else
    TBX(:) = ( STBX(:)*SQRT(VU(:)*VU(:) + TVAR3W(:)*TVAR3W(:)) )*TVAR3W(:)
    TBY(:) = ( STBY(:)*SQRT(UV(:)*UV(:) + TVAR3S(:)*TVAR3S(:)) )*TVAR3S(:)
  endif

  TTBXY = TTBXY + (DSTIME(0)-TTDS)
  
  ! *** *******************************************************************!
  ! *** SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
  if( ISBSDFUF >= 1 )then
    HDFUFM = 1.E-12

    do L = 2,LA
      LS = LSC(L)
      HDFUFX(L) = HDFUFM+G*SUB(L)*HU(L)*(BELV(LWC(L)) - BELV(L))*DXIU(L)
      HDFUFY(L) = HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )    - BELV(L))*DYIV(L)
    enddo

    do L = 2,LA
      if( HDFUFX(L)>0.0 )then
        HDFUFX(L) = TBX(L)/HDFUFX(L)
      else
        HDFUFX(L) = 1.0
      endif
      if( HDFUFY(L)>0.0 )then
        HDFUFY(L) = TBY(L)/HDFUFY(L)
      else
        HDFUFY(L) = 1.0
      endif
    enddo
  endif

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
  TTDS = DSTIME(0)
  if( ISWAVE == 0 .or. LSEDZLJ )then
    if( ISCORTBC == 0 )then
      TVAR3W(:) = TSX(LEC(:))
      TVAR3S(:) = TSY(LNC(:))
      TVAR3E(:) = TBX(LEC(:))
      TVAR3N(:) = TBY(LNC(:))

      do L = 2,LA
        TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2
        QQ(L,0)  = 0.5*CTURB2*SQRT(TMP)

        TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2
        QQ(L,KC) = 0.5*CTURB2*SQRT(TMP)
      enddo

    else    ! if( ISCORTBC >= 1 )then
      ! *** CORNER CORRECTIONS
      do L = 2,LA
        TVAR3S(L) = TSY(LNC(L))
        TVAR3W(L) = TSX(LEC(L))
        TVAR3E(L) = TBX(LEC(L))
        TVAR3N(L) = TBY(LNC(L))
      enddo

      do L = 2,LA
        WCOREST(L) = 1.
        WCORWST(L) = 1.
        WCORNTH(L) = 1.
        WCORSTH(L) = 1.
      enddo

      do L = 2,LA
        if( ISSBCP(L) == 0 )then
          if( SUB(LEC(L)) < 0.5)WCOREST(L) = FSCORTBCV(L)
          if( SUB(L) < 0.5)WCORWST(L) = FSCORTBCV(L)
          if( SVB(LNC(L)) < 0.5)WCORNTH(L) = FSCORTBCV(L)
          if( SVB(L) < 0.5)WCORSTH(L) = FSCORTBCV(L)
        endif
      enddo

      do L = 2,LA
        WCOREW(L) = 1./(WCOREST(L)+WCORWST(L))
        WCORNS(L) = 1./(WCORNTH(L)+WCORSTH(L))
      enddo

      do L = 2,LA
        WCOREST(L) = WCOREST(L)*WCOREW(L)
        WCORWST(L) = WCORWST(L)*WCOREW(L)
        WCORNTH(L) = WCORNTH(L)*WCORNS(L)
        WCORSTH(L) = WCORSTH(L)*WCORNS(L)
      enddo

      do L = 2,LA
        QQ(L,0 )   =     CTURB2*SQRT((RSSBCE(L)*WCOREST(L)*TVAR3E(L) + RSSBCW(L)*WCORWST(L)*TBX(L))**2 + (RSSBCN(L)*WCORNTH(L)*TVAR3N(L) + RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
        QQ(L,KC)   = 0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)            + RSSBCW(L)*TSX(L))**2            + (RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
      enddo
    endif  ! *** END OF ISCORTBC BLOCK
  endif    ! *** END OF ISWAVE = 0


  !----------------------------------------------------------------------!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
  if( ISWAVE >= 1 .and. .not. LSEDZLJ )then

    ! ***  STANDARD WAVE CALCULATIONS
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)

      do L = LF,LL
        TVAR3S(L) = TSY(LNC(L))
        TVAR3W(L) = TSX(LEC(L))
        TVAR3E(L) = TBX(LEC(L)   )
        TVAR3N(L) = TBY(LNC(L))
      enddo

      do L = LF,LL
        if( LWVMASK(L) )then
          TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
          TAUBC  = 0.5*SQRT(TAUBC2)   ! *** CURRENT ONLY
          CTAUC(L) = TAUBC
          UTMP   = 0.5*STCUV(L)*( U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)) )+1.E-12
          VTMP   = 0.5*STCUV(L)*( V(LNC(L),KSZV(LNC(L))) + V(L,KSZV(L)) )
          CURANG = ATAN2(VTMP,UTMP)
          TAUB2  = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
          TAUB2  = max(TAUB2,0.)              ! *** CURRENT & WAVE
          QQ(L,0 )   = CTURB2*SQRT(TAUB2)  ! *** CELL CENTERED TURBULENT INTENSITY DUE TO CURRENTS & WAVES
          QQ(L,KC)   = 0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2)
        else
          TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2
          QQ(L,0 ) = 0.5*CTURB2*SQRT(TMP)

          TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2
          QQ(L,KC) = 0.5*CTURB2*SQRT(TMP)
        endif
      enddo

    enddo  ! *** END OF DOMAIN

  endif    ! *** END OF ISWAVE>0

  ! *** *******************************************************************!
  ! *** CALCULATE TURBULENT INTENSITY SQUARED
  if( KC > 1 )then
    if( ISGOTM > 0 )then 
      do L = 2,LA
        QQ(L,0:KC)  = 2.*tke3d(L,0:KC)
        QQSQR(L,0)  = SQRT(QQ(L,0))
        DML(L,0:KC) = GL3D(L,0:KC)/HP(L)
      enddo
    else
      ! *** Original EFDC turbulent intensity
      if( lpmc) print *, NIter, 'CALQQ2T'  ! delme
      call CALQQ2T
    endif
  endif
  TQQQ = TQQQ + (DSTIME(0)-TTDS)
      
  ! ****************************************************************************
  ! *** MPI communication
  call MPI_barrier(DSIcomm, ierr)
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

  ! *** Compute the sqrt of the turbulence (m/s)
  ! $OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,KS,LDMWET,LAWET,LWET,LLWET,LKWET,QQ,QQSQR) PRIVATE(ND,K,LF,LL,LP,L)
  do ND = 1,NDM  
    do K = 0,KS
      if( K == 0 )then
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)
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
  ! $OMP END PARALLEL DO

  ! *** *******************************************************************!
  ! *** *******************************************************************!
  ! *** HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
  ! *** *******************************************************************!
  ! *** *******************************************************************!

  ! *** Write to time series files
  CTIM = TIMESEC/TCON
  if( ISTMSR >= 1 )then
    if( N >= NBTMSR .and. N <= NSTMSR )then
      if( NCTMSR >= NWTMSR )then

        if( process_id == master_id )then
          call TMSR !@todo add Map_TMSR routines so this will work on all processes?
        endif
        NDIFF = NWTMSR - NCTMSR
        NCTMSR = NINCRMT + NDIFF
      else
        NCTMSR = NCTMSR + NINCRMT
      endif
    endif
  endif

  ! *** *******************************************************************!
  ! *** CALCULATE MEAN MASS TRANSPORT FIELD
  if( (ISSSMMT > 0 .or. ISWASP > 0).and. RESSTEP > 0 ) CALL CALMMT

  ! *** *******************************************************************!
  ! *** ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
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

  ! *** *******************************************************************!
  ! *** WRITE EFDC EXPLORER FORMAT OUTPUT
  if( ISPPH == 1 )then
    ! *** CHECK IF PROPWASH OUTPUT FORCES EE LINKAGE TO BE WRITTEN
    if( ISPROPWASH > 0 )then
      ! *** Determine if PW_Mesh was writting on any process
      call DSI_All_Reduce(iwrite_pwmesh, i, MPI_SUM, TTDS, 0, TWAIT)

      if( i > 0 )then
        if( TIMEDAY >= last_snapshot + freq_out_min/1440. )then
          iwrite_pwmesh = 2
        else
          iwrite_pwmesh = 0
        endif
      endif
    endif
    
    !CALL BCOUT_ACCUMULATE(DTDYN) @todo invoke and make work in MPI case
    if( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )then
      TTDS = DSTIME(0)
      call Map_Write_EE_Binary
      TMPIEE = TMPIEE + (DSTIME(0)-TTDS)
      
      if( process_id == master_id )then
        call EE_LINKAGE(0)
      endif
     if( iwrite_pwmesh > 0 ) last_snapshot = TIMEDAY
    endif
  endif
  if( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )then
    NSNAPSHOTS = NSNAPSHOTS + 1
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
  if( ISHOW > 0 )then
    CALL SHOWVAL(BSYNC)

    if( BSYNC ) call MPI_barrier(DSIcomm, ierr)
  endif
  ! *** *******************************************************************!
  
# ifdef _WIN  
  if( num_Processors == 1 )then     ! *** OMP Run
    if( KEY_PRESSED() )then
      if( ISEXIT() )then
        TTDS = DSTIME(0)
        call Map_Write_EE_Binary
        TMPIEE = TMPIEE + (DSTIME(0)-TTDS)
      
        call EE_LINKAGE(-1)
        GOTO 1000
      endif
    endif
  endif
# endif
  
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
  
  GOTO 1001

  ! *** *******************************************************************!
  ! *** *******************************************************************!
  ! *** TIME LOOP COMPLETED
  ! *** *******************************************************************!
  ! *** *******************************************************************!
1000 continue
  THDMT = THDMT + (DSTIME(0)-T1TMP)

  ! *** *******************************************************************!
  ! *** WRITE RESTART FILE
  if( ABS(ISRESTO) > 0 .and. TIMEDAY > TIME_RESTART-TIMERST*0.5 )then
    if( process_id == master_id ) WRITE(6,'(A,F12.4)') 'FINAL RESTART FILE: ',TIMEDAY
    call Restart_Out(TIME_RESTART, 0)
    if( ISTRAN(8) >= 1 )then
      call WQ_WCRST_OUT(TIME_RESTART)
      if( IWQBEN == 1 ) CALL WQSDRST_OUT(TIME_RESTART)
      if( ISRPEM > 0  ) CALL WQRPEMRST_Out(TIME_RESTART)
    endif
  endif

  ! *** *******************************************************************!
  ! *** OUTPUT COURANT NUMBER DIAGNOSTICS
  if( process_id == master_id )then
    if( ISINWV == 1 )then
      open(1,FILE = OUTDIR//'CFLMAX.OUT')
      close(1,STATUS = 'DELETE')
      open(1,FILE = OUTDIR//'CFLMAX.OUT')

      do L = 2,LA
        write(1,1991) IL(L),JL(L),(CFLUUU(L,K),K = 1,KC)
        write(1,1992) (CFLVVV(L,K),K = 1,KC)
        write(1,1992) (CFLWWW(L,K),K = 1,KC)
        write(1,1992) (CFLCAC(L,K),K = 1,KC)
      enddo

      close(1)
1991  FORMAT(2I5,12F8.3)
1992  FORMAT(10X,12F8.3)
    endif
    ! *** *******************************************************************!
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

    ! *** *******************************************************************!
    ! *** OUTPUT FINAL FOOD CHAIN AVERAGING PERIOD
    if( ISTRAN(5) >= 1 .and. ISFDCH >= 1 )CALL FOODCHAIN(1)

  endif

  if( process_id == master_id )then
    ! *** *******************************************************************!
    close(90)
    close(98)

  endif

END
