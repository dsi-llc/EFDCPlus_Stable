! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details EXECUTES THE FULL HYDRODYNAMIC AND MASS
!! TRANSPORT TIME INTEGRATION USING A TWO TIME LEVEL SCHEME
!
! @date    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
!!         2014-08-12    D H CHUNG         SET EXPLICIT PRECISIONS OF INTEGER & REAL
!!         2011-03       Paul M. Craig     Rewritten to F90 and added OMP
!!         2010-01       Chung Dang        Added the DSI version of Lagrangian Particle Tracking
!!         09-22-2004    Paul M. Craig     Merged DS and TT versions with the 06-04-2004 TT code
!!         05/01/2002    John Hamrick      Modified calls to calbal and budget subroutines
!!                                            added calls to bal2t1, bal2t4, bal2t5
!!         05/02/2002    John Hamrick      Modified calculation of cell center bed stress (stored as QQ(l,0))
!!                                         for cells have source/sinks



SUBROUTINE HDMT2T

  ! *** *******************************************************************!
  USE GLOBAL
  USE DRIFTER  ,ONLY:DRIFTER_CALC
  USE WINDWAVE ,ONLY:WINDWAVEINIT,WINDWAVETUR,READWAVECELLS
  USE HIFREQOUT
  USE RESTART_MODULE
  USE EFDCOUT
  USE CALCSERMOD,ONLY: CALCSER
  USE FIELDS
  Use IFPORT
  USE WATERQUALITY,ONLY:WQ3D
  Use Variables_WQ
  Use Variables_Propwash

  USE Variables_MPI
  USE OMP_LIB
#ifdef NCOUT
  USE MOD_NETCDF
#endif
  USE MPI
  Use Communicate_Ghost_Routines
  Use Mod_Map_Write_EE_Binary
  Use Mod_Map_Gather_Sort
  Use Mod_Map_Write_NetCDF
  Use MPI_All_Reduce

  IMPLICIT NONE

  INTEGER :: ICALLTP, ILOGC, K, L, LP, LE, LW, LN, LS, LNW, LSE, LSW, ND, LF, LL, LG
  INTEGER :: NLOOP, NTMP1, NTMP2, NMD, NDIFF, NS, NT

  REAL      :: TMP, HDFUFM, TAUBC2, TAUBC, UTMP, VTMP, CURANG, TAUB2, DELTD2, DZDDELT
  REAL      :: CTIM, WTM, WTMP, DELVOL, USGZ, VSGZ
  REAL(8)   :: DEL
  LOGICAL   :: LDEBUG
  CHARACTER :: LFILENAME*30

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WCOREW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WCORNS
  REAL(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS, T1TMP, TWAIT                 ! MODEL TIMING TEMPORARY VARIABLES
  REAL(RKD)           :: TIMEHARM, TIMELAST, DAYNEXT

  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LCORNER
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LCORNWE
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LCORNSN

#ifdef _WIN
  LOGICAL, EXTERNAL :: KEY_PRESSED
  LOGICAL, EXTERNAL :: ISEXIT
#endif
  LOGICAL   :: RES
  REAL,SAVE :: ADJC
  
  ! *** New variables for MPI
  Integer :: ierr, i, j
  Integer :: Start_Local_LA
  Integer :: End_Local_LA

  REAL, DIMENSION(2) :: DYN_IN, DYN_OUT
  
  ! *** ALLOCATE LOCAL ARRAYS
  IF(  .NOT. ALLOCATED(WCOREW) )THEN
    ALLOCATE(WCOREW(LCM))
    ALLOCATE(WCORNS(LCM))
    ALLOCATE(LCORNER(LCM))
    ALLOCATE(LCORNWE(LCM))
    ALLOCATE(LCORNSN(LCM))
    WCOREW=0.0
    WCORNS=0.0
    LCORNER=0
    LCORNWE=0
    LCORNSN=0

    IF( KC > 1 )THEN
      ADJC = 1.0
    ELSE
      ADJC = 1.1
    ENDIF
    TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8
    DAYNEXT = TIMEDAY + 1.0
    TIMELAST = TIMEDAY + TIMERST
    TIMEHARM = DBLE(TCON)*DBLE(TBEGIN) + 300.   ! *** 05 MINUTE INCREMENTS
  ENDIF

  if( process_id == master_id )then
    WRITE(*,'(A)')'STARTING HDMT 2TL'
  endif
  T1TMP  = DSTIME(0)
  ICALLTP = 0

  ISTL=2
  IS2TL=1
  FOURDPI=4./PI
  RES = .TRUE.
  LDEBUG=.FALSE.

  ! *** INITIALIZE TIME FOR INITIAL CONDITION COMPUTATIONS
  TIMESEC = TCON*TBEGIN

  ! *** *******************************************************************!
  ! *** SET FLAGS FOR CORNER CELL BED STRESS CORRECTIONS
  IF( ISCORTBC >= 1 )THEN
    ! *** SET FLAG FOR CELLS HAVING VOLUME SOURCE OR SINKS
    DO L=1,LC
      ISSBCP(L)=0
    ENDDO

    DO L=2,LA
      IF( RSSBCE(L) > 1.5)ISSBCP(L)=1
      IF( RSSBCW(L) > 1.5)ISSBCP(L)=1
      IF( RSSBCN(L) > 1.5)ISSBCP(L)=1
      IF( RSSBCS(L) > 1.5)ISSBCP(L)=1
    ENDDO
  ENDIF

  DO L=2,LA
    WCOREST(L)=1.
    WCORWST(L)=1.
    WCORNTH(L)=1.
    WCORSTH(L)=1.
  ENDDO

  ! *** *******************************************************************!
  ! *** REINITIALIZE VARIABLES
  IF( ISRESTI == 0 )THEN
    DO L=2,LA
      H1P(L)=HP(L)
      H1U(L)=HU(L)
      H1UI(L)=HUI(L)
      H1V(L)=HV(L)
      H1VI(L)=HVI(L)
      UHDY1E(L)=UHDYE(L)
      VHDX1E(L)=VHDXE(L)
    ENDDO

    DO K=1,KC
      DO L=2,LA
        HPK(L,K) = HP(L)*DZC(L,K)
        H1PK(L,K) = HPK(L,K)
        U1(L,K)=U(L,K)
        V1(L,K)=V(L,K)
        UHDYF1(L,K)=UHDYF(L,K)
        VHDXF1(L,K)=VHDXF(L,K)
        UHDY1(L,K)=UHDY(L,K)
        VHDX1(L,K)=VHDX(L,K)
      ENDDO
    ENDDO
  ENDIF

  ! *** *******************************************************************!
  ! *** INITIALIZE COURANT NUMBER DIAGNOSTICS
  IF( ISINWV == 1 )THEN
    DO K=1,KC
      DO L=2,LA
        CFLUUU(L,K)=0.
        CFLVVV(L,K)=0.
        CFLWWW(L,K)=0.
        CFLCAC(L,K)=0.
      ENDDO
    ENDDO
  ENDIF
  ILOGC=0


  ! *** *******************************************************************!
  ! *** CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
  ! *** CALCULATE VELOCITY GRADIENTS
  !----------------------------------------------------------------------!
  DO L=2,LA
    LE=LEC(L)
    LN=LNC(L)
    LS=LSC(L)
    LNW=LNWC(L)
    LSE=LSEC(L)
    LSW=LSWC(L)

    UV(L) =0.25*( HP(LS  )*(U(LSE, KSZV(L))+U(LS,  KSZV(L)) ) + HP(L )*(U(LE, KSZV(L))+U(L, KSZV(L)) ) )*HVI(L)
    U1V(L)=0.25*( H1P(LS )*(U1(LSE,KSZV(L))+U1(LS, KSZV(L)) ) + H1P(L)*(U1(LE,KSZV(L))+U1(L,KSZV(L)) ) )*H1VI(L)
    VU(L) =0.25*( HP(LWC(L) )*(V(LNW, KSZU(L))+V(LWC(L), KSZU(L)) ) + HP(L )*(V(LN,  KSZU(L))+V(L, KSZU(L)) ) )*HUI(L)
    V1U(L)=0.25*( H1P(LWC(L))*(V1(LNW,KSZU(L))+V1(LWC(L),KSZU(L)) ) + H1P(L)*(V1(LN, KSZU(L))+V1(L,KSZU(L)) ) )*H1UI(L)
  ENDDO

  ! *** MPI Communication of ghost cells
#ifdef _MPI
  Call communicate_ghost_cells(UV,  'UV')
  Call communicate_ghost_cells(U1V, 'U1V')
  Call communicate_ghost_cells(VU,  'VU')
  Call communicate_ghost_cells(V1U, 'V1U')
#endif

  ! *** *******************************************************************!
  ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8

  IF( ISWAVE == 1 ) CALL WAVEBL
  IF( ISWAVE == 2 ) CALL WAVESXY
  IF( ISWAVE >= 3 .AND. NWSER > 0 )THEN
    CALL WINDWAVEINIT
    CALL WINDWAVETUR   !DHC FIRST CALL
    ! *** READ IN WAVE COMPUTATIONAL CELL LIST
    IF( IUSEWVCELLS /= 0 )THEN
      CALL READWAVECELLS
    ENDIF
  ENDIF

#ifdef NCOUT
  ! ** NETCDF INIT
  NC_DATESTAMP = ''
  IF( NCDFOUT > 0 )THEN
    ! *** Call subroutine to gather and map variables in prep for writing to the NetCDF file(s)
    IF( TIMEDAY >= TBEGNCDF ) Call Map_Write_NetCDF

    ! *** Only do NETCDF writing on the master
    IF( process_id == master_id )THEN
      CALL READCORN
      IF( TIMEDAY >= TBEGNCDF ) CALL nc_output()
    end if ! *** end on master

  ENDIF
#endif  
  
!> @todo modify writing out for hi frequency output with MPI
  if( process_id == master_id )THEN
    IF( HFREOUT == 1 )THEN
      DO NS=1,NSUBSET
#ifdef NCOUT          
        IF( NCDFOUT > 0 ) call nc_output_hf(1,ns)
#endif 
        CALL HFREHYOUT(1,NS)
        CALL HFREWCOUT(1,NS)
        IF( ISTRAN(8) >= 1 ) THEN
          CALL HFREWQOUT(1,NS)
          IF( ISRPEM > 0 ) CALL HFRERPEMOUT(1,NS)
        ENDIF
      ENDDO
    ENDIF
  endif

  ! *** *******************************************************************!
  ! *** FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  DTDYN = DT  ! *** PMC - FOR INITIALIZATION
  CALL CALTBXY(ISTL,IS2TL)

  ! *** *******************************************************************!
  ! *** CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
  IF( ISHDMF >= 1 ) CALL CALHDMF

  ! *** *******************************************************************!
  ! *** CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
  !----------------------------------------------------------------------!
  N=-1
  CALL CALTSXY

  ! *** *******************************************************************!
  ! *** SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  CALL CALTBXY(ISTL,IS2TL)

  ! *** SET BOTTOM AND SURFACE STRESSES
  DO L=2,LA
    USGZ = U(L,KSZU(L))
    VSGZ = V(L,KSZV(L))
    TBX(L) = ( STBX(L)*SQRT(VU(L)*VU(L) + USGZ*USGZ) )*USGZ
    TBY(L) = ( STBY(L)*SQRT(UV(L)*UV(L) + VSGZ*VSGZ) )*VSGZ
  ENDDO

  ! *** MPI Communication routines
#ifdef _MPI
  Call communicate_ghost_cells(TBX, 'TBX')
  Call communicate_ghost_cells(TBY, 'TBY')
#endif

  N=0

  CALL CALTSXY

  !----------------------------------------------------------------------!
  ! *** SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
  DO L=2,LA
    HDFUFX(L)=1.
    HDFUFY(L)=1.
    HDFUF(L)=1.
  ENDDO

  ! ***
  IF( ISBSDFUF >= 1 )THEN
    HDFUFM=1.E-12

    DO L=2,LA
      LS=LSC(L)
      HDFUFX(L)=HDFUFM+G*SUB(L)*HU(L)*(BELV(LWC(L))-BELV(L))*DXIU(L)
      HDFUFY(L)=HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )-BELV(L))*DYIV(L)
    ENDDO

    DO L=2,LA
      IF( HDFUFX(L)>0.0 )THEN
        HDFUFX(L)=TBX(L)/HDFUFX(L)
      ELSE
        HDFUFX(L)=1.0
      ENDIF
      IF( HDFUFY(L)>0.0 )THEN
        HDFUFY(L)=TBY(L)/HDFUFY(L)
      ELSE
        HDFUFY(L)=1.0
      ENDIF
    ENDDO

  ENDIF

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  IF( ISWAVE == 0 )THEN

    IF( ISCORTBC == 0 )THEN
      ! *** NO CORNER CORRECTIONS - STANDARD APPROACH
      !$OMP PARALLEL DEFAULT(SHARED)

      !$OMP DO PRIVATE(ND,LF,LL,LP,L)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)
          TVAR3S(L)=TSY(LNC(L))
          TVAR3W(L)=TSX(LEC(L))
          TVAR3E(L)=TBX(LEC(L))
          TVAR3N(L)=TBY(LNC(L))
        ENDDO
      ENDDO
      !$OMP END DO

      ! *** COMPUTE CELL CENTERED BOTTOM AND SURFACE
      !$OMP DO PRIVATE(ND,LF,LL,LP,L)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)
          QQ(L,0)  = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2 )
          QQ(L,KC) = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2 + (RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2 )
          QQSQR(L,0) = SQRT(QQ(L,0))
        ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

    ELSE                             !IF(ISCORTBC >= 1 )THEN

      ! *** WITH CORNER CORRECTIONS
      DO L=2,LA
        IF( ISSBCP(L) == 0 )THEN
          IF( SUB(LEC(L)) < 0.5) WCOREST(L)=FSCORTBCV(L)
          IF( SUB(L) < 0.5) WCORWST(L)=FSCORTBCV(L)
          IF( SVB(LNC(L)) < 0.5) WCORNTH(L)=FSCORTBCV(L)
          IF( SVB(L) < 0.5) WCORSTH(L)=FSCORTBCV(L)
        ENDIF
      ENDDO

      DO L=2,LA
        WCOREW(L)=1./(WCOREST(L)+WCORWST(L))
        WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))
      ENDDO

      DO L=2,LA
        WCOREST(L)=WCOREST(L)*WCOREW(L)
        WCORWST(L)=WCORWST(L)*WCOREW(L)
        WCORNTH(L)=WCORNTH(L)*WCORNS(L)
        WCORSTH(L)=WCORSTH(L)*WCORNS(L)
      ENDDO

      DO L=2,LA
        TVAR3S(L)=TSY(LNC(L))
        TVAR3W(L)=TSX(LEC(L))
        TVAR3E(L)=TBX(LEC(L)   )
        TVAR3N(L)=TBY(LNC(L))
      ENDDO

      DO L=2,LA
        QQ(L,0 )=CTURB2*SQRT(  (RSSBCE(L)*WCOREST(L)*TVAR3E(L) + RSSBCW(L)*WCORWST(L)*TBX(L))**2  + (RSSBCN(L)*WCORNTH(L)*TVAR3N(L) + RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
        QQ(L,KC)=0.5*CTURB2*SQRT(  (RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX(L))**2  +(RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
        QQSQR(L,0)=SQRT(QQ(L,0))
      ENDDO

    ENDIF
  ENDIF    ! *** END OF ISWAVE=0

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  !----------------------------------------------------------------------!
  IF( ISWAVE >= 1 )THEN

    DO L=2,LA
      TVAR3S(L)=TSY(LNC(L))
      TVAR3W(L)=TSX(LEC(L))
      TVAR3E(L)=TBX(LEC(L)   )
      TVAR3N(L)=TBY(LNC(L))
    ENDDO

    DO L=2,LA
      TAUBC2=(RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
      TAUBC=0.5*SQRT(TAUBC2)   ! *** CURRENT ONLY
      CTAUC(L)=TAUBC
      UTMP=0.5*STCUV(L)*(U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)))+1.E-12
      VTMP=0.5*STCUV(L)*(V(LN ,KSZV(LN )) + V(L,KSZV(L)))
      CURANG=ATAN2(VTMP,UTMP)
      TAUB2=TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
      TAUB2=MAX(TAUB2,0.)          ! *** CURRENT & WAVE
      QQ(L,0 )=CTURB2*SQRT(TAUB2)  ! *** CELL CENTERED TURBULENT INTENSITY DUE TO CURRENTS & WAVES
      QQ(L,KC)=0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2)
      QQSQR(L,0)=SQRT(QQ(L,0))
    ENDDO
  ENDIF

  ! ***  SET GRAIN STRESS
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    DO L=2,LA
      TAUBSED(L) = QQ(L,0)/CTURB2
      TAUBSND(L) = QQ(L,0)/CTURB2
    ENDDO
  ENDIF

  ! *** *******************************************************************!
  ! ***  SET SWITCHES FOR TWO TIME LEVEL INTEGRATION
  DELT=DT
  DELTD2=DT/2.
  DZDDELT=DZ/DELT

  ! *** *******************************************************************!
  ! *** BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION
  ! *** SET CYCLE COUNTER AND CALL TIMER
  NTIMER=0
  N=0

  ! *** EFDC_EXPLORER BEGIN BLOCK  RECORD TIME
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8

  if( process_id == master_id )THEN
    CALL TIMELOG(N,TIMEDAY,OUTDIR,0._8)
  End if

  NTIMER=1
  NINCRMT=1
  NITER=0
  NLOOP=0

  ! *** ************************************************************************
  ! *** ************************************************************************
  ! *** BEGINNING OF THE MAIN TIME ITERATION LOOP FOR TWO TIME LEVEL SOLUTION
  ! *** ************************************************************************
  ! *** ************************************************************************
1001 CONTINUE
  IF( TIMEDAY > TIMEEND ) GO TO 1000

  IF( ISDYNSTP == 0 )THEN
    N = N + 1
  ELSE
    NLOOP = NLOOP+1
    IF( NLOOP > NRAMPUP )THEN
      CALL CALSTEPD
    ELSE
      DTDYN = DT
      NINCRMT = 1
    ENDIF
    DELT    = DTDYN
    DELTD2  = DTDYN/2.
    DZDDELT = DZ/DTDYN
    N = N + NINCRMT
  ENDIF
  NITER = NITER + 1

  ETIMESEC = DT*FLOAT(N)
  ETIMEDAY = DT*FLOAT(N)/86400.
  TIMESEC = DBLE(DT)*DBLE(N)+DBLE(TCON)*DBLE(TBEGIN)
  TIMEDAY = TIMESEC/86400._8

  IF( ISDYNSTP == 0 )THEN
    ILOGC=ILOGC+1
  ELSE
    ILOGC=ILOGC+NINCRMT
  ENDIF

  ! *** DSI BEGIN BLOCK
  IF( N <= NLTS )THEN
    SNLT=0.
  ELSEIF( N > NLTS .AND. N <= NTTS )THEN
    NTMP1=N-NLTS
    NTMP2=NTTS-NLTS+1
    SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)
  ELSE
    SNLT=1.
  ENDIF

  ! *** SET THE GRAVITY ACCELERATION
  IF( N <= NTSVB )THEN
    GP=GPO*(FLOAT(N)/FLOAT(NTSVB))
  ELSE
    GP=GPO
  ENDIF

  ! 2018-11-22, NTL: UPDATE TIME VARIABLE TOPOGRAPHY
  IF( BATHY.IFLAG > 0 ) CALL UPDATETOPO(TIMEDAY,BELV,HP,HDRY)
  IF( ROUGH.IFLAG > 0 ) CALL UPDATEFIELD(ROUGH,TIMEDAY,1,ZBR)

  !----------------------------------------------------------------------!
  ! *** INITIALIZE TWO-TIME LEVEL BALANCES
  IF( ISBAL >= 1 ) CALL BAL2T1

  ! *** *******************************************************************!
  ! *** CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
  IF( KC > 1 )THEN
    TTDS = DSTIME(0)

    CALL CALAVB
    
    TAVB = TAVB + (DSTIME(0)-TTDS)
  ENDIF

  ! *** *******************************************************************!
  ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
  TTDS = DSTIME(0)
  IF( ISWAVE == 1 ) CALL WAVEBL
  IF( ISWAVE == 2 ) CALL WAVESXY
  IF( ISWAVE >= 3 .AND. NWSER > 0 ) CALL WINDWAVETUR   !DHC NEXT CALL
  TTBXY = TTBXY + (DSTIME(0)-TTDS)

  ! *** *******************************************************************!
  ! *** UPDATE TIME VARIABLE SURFACE WIND STRESSES
  TTDS = DSTIME(0)
  TMPITMP = 0.
  CALL CALTSXY
  TTBXY = TTBXY + (DSTIME(0)-TTDS) - TMPITMP
  DSITIMING(9) = DSITIMING(9) + TMPITMP

  ! *** *******************************************************************!
  ! *** CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
  TTDS = DSTIME(0)
  CALL CALEXP2T
  TCEXP = TCEXP + (DSTIME(0)-TTDS)

  ! *** *******************************************************************!
  ! *** UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS,
  ! *** VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
  CALL CALCSER(ISTL)
  CALL CALVEGSER
  
  CALL CALQVS(ISTL)
  PSERT(0) = 0.
  IF( NPSER >= 1 ) CALL CALPSER

  ! *** *******************************************************************!
  ! *** SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
  TTDS = DSTIME(0)
  TMPITMP = 0.
  CALL CALPUV2C
  TPUV = TPUV + (DSTIME(0)-TTDS) - TMPITMP - TTWAIT          ! *** Includes CONGRAD time (TTWAIT is the wait time)
  DSITIMING(2) = DSITIMING(2)    + TMPITMP

  ! *** *******************************************************************!
  ! *** ADVANCE INTERNAL VARIABLES
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,LP,L)
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        UHDYF1(L,K) = UHDYF(L,K)
        VHDXF1(L,K) = VHDXF(L,K)
        UHDY1(L,K)  = UHDY(L,K)
        VHDX1(L,K)  = VHDX(L,K)
        U1(L,K)     = U(L,K)
        V1(L,K)     = V(L,K)
        W1(L,K)     = W(L,K)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  ! *** *******************************************************************!
  ! *** SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
  ! *** USING THE INTERNAL SHEARS STORED IN DU & DV
  TTDS = DSTIME(0)
  TMPITMP = 0.
  IF( KC > 1 )THEN
    CALL CALUVW (ISTL)
  ELSE
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LA,LDM,UHDYE,UHDYF,UHDY,HUI,DYIU,U,VHDXE,VHDXF,VHDX,HVI,DXIV,V,W) PRIVATE(ND,LF,LL,L)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)

      DO L=LF,LL
        UHDYF(L,1) = UHDYE(L)
        UHDY(L,1)  = UHDYE(L)
        U(L,1)     = UHDYE(L)*HUI(L)*DYIU(L)
        VHDXF(L,1) = VHDXE(L)
        VHDX(L,1)  = VHDXE(L)
        V(L,1)     = VHDXE(L)*HVI(L)*DXIV(L)
        W(L,1)     = 0.
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    CALL CALUVW(ISTL)
  ENDIF
  TUVW = TUVW + (DSTIME(0)-TTDS) - TMPITMP
  DSITIMING(3) = DSITIMING(3) + TMPITMP

  ! *** *******************************************************************!
  ! *** CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
  ! *** AT TIME LEVEL (iN+1)
  TMPITMP = 0.
  IF( ISTRANACTIVE > 0 )then
    CALL CALCONC(ISTL,IS2TL)
  ELSE
    ! *** Call Propwash module if that option is selected
    IF( propwash_on )THEN
      Call Propwash_Calc_Sequence(0)
    END IF
  ENDIF
  
  ! *** Communicate variables with MPI calls before writing to files and looping computations
#ifdef _MPI
  Call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)
  Call Communicate_CON1(0)
  
  TMPITMP = DSTIME(0) - TTDS
  DSITIMING(6) = DSITIMING(6) + TMPITMP
#endif

  ! *** *******************************************************************!
  ! *** CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
  ! *** CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOUBLE TIME
  ! *** STEP TRANSPORT FIELD
  IF( ISWQFLUX == 1 )THEN
    IF( ISICM >= 1 .OR. ISRPEM > 0 .OR. ISTRAN(4) >= 1 )THEN
      DO L=2,LA
        HWQ(L) = HP(L)
      ENDDO

      ! *** ADD CHANNEL INTERACTIONS
      IF( MDCHH >= 1 )THEN
        DO NMD=1,MDCHH
          IF( MDCHTYP(NMD) == 1 )THEN
            HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))  +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))
            HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))  -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
          ENDIF
          IF( MDCHTYP(NMD) == 2 )THEN
            HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
            HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD)) - DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
          ENDIF
          IF( MDCHTYP(NMD) == 3 )THEN
            HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))
            HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD)) - DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))
            HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD)) - DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))
          ENDIF
        ENDDO
      ENDIF      ! *** END ADD CHANNEL INTERACTIONS
    ENDIF
    
    ! *** CALL WATER QAULITY KINETICS, SEDIMENT AND RPEM PROCESSES
    IF( ISTRAN(8) >= 1 ) CALL WQ3D
    IF( ISTRAN(4) >= 1 ) CALL CALSFT(ISTL,IS2TL)

    IF( ISICM >= 1 )THEN
      DO L=2,LA
        H2WQ(L) = HWQ(L)
      ENDDO
    ENDIF

    ! *** Communicate variables with MPI calls before writing to files and looping computations
#ifdef _MPI
    IF( ISTRAN(8) >= 1 )THEN
      Call MPI_barrier(MPI_Comm_World, ierr)
      TTDS = DSTIME(0)
      CALL communicate_ghost_cells(WQV, NWQV)
      TMPITMP = DSTIME(0) - TTDS
      DSITIMING(6) = DSITIMING(6) + TMPITMP
    ENDIF
#endif
  ENDIF         ! *** END OF WQ SECTION

  ! *** *******************************************************************!
  ! ***  UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING AN EQUATION OF STATE
  IF( BSC > 1.E-6  )THEN
    CALL CALBUOY(.TRUE.)
  End if

  ! *** CALL TWO-TIME LEVEL BALANCES
  IF( ISBAL >= 1 ) CALL BAL2T4

  ! *** *******************************************************************!
  ! *** CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)
  IF( ICALTB == 0 )THEN
    ! *** EFDC+ 10.0 AND LATER APPROACH TO COMPUTE THE SHEAR VELOCITIES AT THE U AND V FACES

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,LE,LS,LSE,LN,LW,LNW)
    DO ND=1,NDM
      LF=(ND-1)*LDMWET+1
      LL=MIN(LF+LDMWET-1,LAWET)

      ! *** ADVANCE THE VARIABLES
      DO LP=LF,LL
        L=LWET(LP)
        U1V(L) = UV(L)
        V1U(L) = VU(L)
      ENDDO
      
      !DIR$ NOFUSION
      DO LP=LF,LL
        L = LWET(LP)
        LSE = LSEC(L)
        LS  = LSC(L)
        LE  = LEC(L)
        LNW = LNWC(L)
        LW  = LWC(L)
        LN  = LNC(L)

        UV(L) = SVB(L)*( 0.5*HP(LS)*(U(LSE,KSZU(LSE)) + U(LS,KSZU(LS))) + 0.5*HP(L)*(U(LE,KSZU(LE)) + U(L,KSZU(L))) )*HVI(L)
        VU(L) = SUB(L)*( 0.5*HP(LW)*(V(LNW,KSZV(LNW)) + V(LW,KSZV(LW))) + 0.5*HP(L)*(V(LN,KSZV(LN)) + V(L,KSZV(L))) )*HUI(L)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ELSE
    ! *** LEGACY EFDC APPLICATIONS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,LE,LS,LSE,LN,LW,LNW)
    DO ND=1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)
        LSE=LSEC(L)
        LS=LSC(L)
        LE=LEC(L)
        LNW=LNWC(L)
        LW=LWC(L)
        LN=LNC(L)

        U1V(L)=UV(L)
        V1U(L)=VU(L)

        UV(L) = ( ( HP(LSE)*U(LSE,KSZU(LSE)) + HP(LS)*U(LS,KSZU(LS)) ) + ( HP(LE)*U(LE,KSZU(LE)) + HP(L)*U(L,KSZU(L)) ) )/( 1.0 + SUB(LE) + SUB(LS) + SUB(LSE) )*HVI(L)
        VU(L) = ( ( HP(LNW)*V(LNW,KSZV(LNW)) + HP(LW)*V(LW,KSZV(LW)) ) + ( HP(LN)*V(LN,KSZV(LN)) + HP(L)*V(L,KSZV(L)) ) )/( 1.0 + SVB(LN) + SVB(LW) + SVB(LNW) )*HUI(L)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  ENDIF

  ! *** *******************************************************************!
  ! *** CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
  ! *** AT TIME LEVEL (N)
  IF( ISHDMF >= 1 )THEN
    TTDS = DSTIME(0)

    CALL CALHDMF
    
    THMDF = THMDF + (DSTIME(0)-TTDS) - TMPITMP
  ENDIF

  ! *** *******************************************************************!
  ! *** CALCULATE BOTTOM STRESS AT LEVEL (N+1)
  TTDS = DSTIME(0)
  CALL CALTBXY(ISTL,IS2TL)

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,AVCON1,USGZ,VSGZ)       &
  !$OMP             SHARED(ICALTB,ISAVCOMP,NDM,LDMWET,LAWET,LWET,KSZ,KSZU,KSZV) &
  !$OMP             SHARED(AVCON,DZI,AVO,HUI,HVI,TBX,STBX,VU,U,TBY,STBY,UV,V)
  DO ND=1,NDM
    LF=(ND-1)*LDMWET+1
    LL=MIN(LF+LDMWET-1,LAWET)

    IF( ICALTB > 0 .AND. ISAVCOMP == 0 )THEN
      DO LP=LF,LL
        L=LWET(LP)
        USGZ = U(L,KSZU(L))
        VSGZ = V(L,KSZV(L))
        TBX(L) = AVCON1*HUI(L)*USGZ
        TBY(L) = AVCON1*HVI(L)*VSGZ
      ENDDO
    ELSE
      DO LP=LF,LL
        L = LWET(LP)
        USGZ = U(L,KSZU(L))
        VSGZ = V(L,KSZV(L))
        TBX(L) = ( STBX(L)*SQRT(VU(L)*VU(L) + USGZ*USGZ) )*USGZ
        TBY(L) = ( STBY(L)*SQRT(UV(L)*UV(L) + VSGZ*VSGZ) )*VSGZ
      ENDDO
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO
  TTBXY = TTBXY + (DSTIME(0)-TTDS)
  
  ! *** *******************************************************************!
  ! *** SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
  IF( ISBSDFUF >= 1 )THEN
    HDFUFM=1.E-12

    DO L=2,LA
      LS=LSC(L)
      HDFUFX(L) = HDFUFM+G*SUB(L)*HU(L)*(BELV(LWC(L)) - BELV(L))*DXIU(L)
      HDFUFY(L) = HDFUFM+G*SVB(L)*HV(L)*(BELV(LS )    - BELV(L))*DYIV(L)
    ENDDO

    DO L=2,LA
      IF( HDFUFX(L)>0.0 )THEN
        HDFUFX(L)=TBX(L)/HDFUFX(L)
      ELSE
        HDFUFX(L)=1.0
      ENDIF
      IF( HDFUFY(L)>0.0 )THEN
        HDFUFY(L)=TBY(L)/HDFUFY(L)
      ELSE
        HDFUFY(L)=1.0
      ENDIF
    ENDDO
  ENDIF

  ! *** *******************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
  TTDS = DSTIME(0)
  IF( ISWAVE == 0 )THEN
    IF( ISCORTBC == 0 )THEN
      ! ***  STANDARD CALCULATIONS - NO CORNER CORRECTS
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,L,K,LF,LL,LP,TMP)           &
      !$OMP             SHARED(NDM,LDMWET,LAWET,LWET,KC,LEC,LNC)             &
      !$OMP             SHARED(TVAR3S,TVAR3W,TVAR3E,TVAR3N,TBX,TBY,TSX,TSY)  &
      !$OMP             SHARED(CTURB2,RSSBCE,RSSBCW,RSSBCN,RSSBCS,QQ,QQSQR)
      DO ND=1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)

        DO LP=LF,LL
          L = LWET(LP)
          TVAR3W(L) = TSX(LEC(L))
          TVAR3S(L) = TSY(LNC(L))
          TVAR3E(L) = TBX(LEC(L))
          TVAR3N(L) = TBY(LNC(L))
        ENDDO

        DO LP=LF,LL
          L = LWET(LP)
          TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2
          QQ(L,0)  = 0.5*CTURB2*SQRT(TMP)

          TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2
          QQ(L,KC) = 0.5*CTURB2*SQRT(TMP)

          QQSQR(L,0)=SQRT(QQ(L,0))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ELSE    ! IF( ISCORTBC >= 1 )THEN
      ! *** CORNER CORRECTIONS
      DO L=2,LA
        TVAR3S(L)=TSY(LNC(L))
        TVAR3W(L)=TSX(LEC(L))
        TVAR3E(L)=TBX(LEC(L)   )
        TVAR3N(L)=TBY(LNC(L))
      ENDDO

      DO L=2,LA
        WCOREST(L)=1.
        WCORWST(L)=1.
        WCORNTH(L)=1.
        WCORSTH(L)=1.
      ENDDO

      DO L=2,LA
        IF( ISSBCP(L) == 0 )THEN
          IF( SUB(LEC(L)) < 0.5)WCOREST(L)=FSCORTBCV(L)
          IF( SUB(L) < 0.5)WCORWST(L)=FSCORTBCV(L)
          IF( SVB(LNC(L)) < 0.5)WCORNTH(L)=FSCORTBCV(L)
          IF( SVB(L) < 0.5)WCORSTH(L)=FSCORTBCV(L)
        ENDIF
      ENDDO

      DO L=2,LA
        WCOREW(L)=1./(WCOREST(L)+WCORWST(L))
        WCORNS(L)=1./(WCORNTH(L)+WCORSTH(L))
      ENDDO

      DO L=2,LA
        WCOREST(L)=WCOREST(L)*WCOREW(L)
        WCORWST(L)=WCORWST(L)*WCOREW(L)
        WCORNTH(L)=WCORNTH(L)*WCORNS(L)
        WCORSTH(L)=WCORSTH(L)*WCORNS(L)
      ENDDO

      DO L=2,LA
        QQ(L,0 )   =     CTURB2*SQRT((RSSBCE(L)*WCOREST(L)*TVAR3E(L) + RSSBCW(L)*WCORWST(L)*TBX(L))**2 + (RSSBCN(L)*WCORNTH(L)*TVAR3N(L) + RSSBCS(L)*WCORSTH(L)*TBY(L))**2)
        QQ(L,KC)   = 0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)            + RSSBCW(L)*TSX(L))**2            + (RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY(L))**2)
        QQSQR(L,0) = SQRT(QQ(L,0))
      ENDDO
    ENDIF  ! *** END OF ISCORTBC BLOCK
  ENDIF    ! *** END OF ISWAVE=0


  !----------------------------------------------------------------------!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
  IF( ISWAVE >= 1 )THEN

    ! ***  STANDARD WAVE CALCULATIONS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,TMP,TAUBC2,TAUBC,UTMP,VTMP,CURANG,TAUB2)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)

      DO L=LF,LL
        TVAR3S(L)=TSY(LNC(L))
        TVAR3W(L)=TSX(LEC(L))
        TVAR3E(L)=TBX(LEC(L)   )
        TVAR3N(L)=TBY(LNC(L))
      ENDDO

      DO L=LF,LL
        IF( LWVMASK(L) )THEN
          TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
          TAUBC  = 0.5*SQRT(TAUBC2)   ! *** CURRENT ONLY
          CTAUC(L) = TAUBC
          UTMP   = 0.5*STCUV(L)*( U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)) )+1.E-12
          VTMP   = 0.5*STCUV(L)*( V(LNC(L),KSZV(LNC(L))) + V(L,KSZV(L)) )
          CURANG = ATAN2(VTMP,UTMP)
          TAUB2  = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
          TAUB2  = MAX(TAUB2,0.)              ! *** CURRENT & WAVE
          QQ(L,0 )   = CTURB2*SQRT(TAUB2)  ! *** CELL CENTERED TURBULENT INTENSITY DUE TO CURRENTS & WAVES
          QQ(L,KC)   = 0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2)
          QQSQR(L,0) = SQRT(QQ(L,0))
        ELSE
          TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2
          QQ(L,0 ) = 0.5*CTURB2*SQRT(TMP)

          TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2
          QQ(L,KC) = 0.5*CTURB2*SQRT(TMP)

          QQSQR(L,0) = SQRT(QQ(L,0))
        ENDIF
      ENDDO

    ENDDO  ! *** END OF DOMAIN
    !$OMP END PARALLEL DO

  ENDIF    ! *** END OF ISWAVE>0

  ! *** *******************************************************************!
  ! *** CALCULATE TURBULENT INTENSITY SQUARED
  IF( KC > 1 )THEN
    CALL CALQQ2T
  ENDIF
  TQQQ = TQQQ + (DSTIME(0)-TTDS)
      
#ifdef _MPI
  Call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)
  Call Communicate_1D2(UV, VU)

  IF( BSC > 1E-6 )THEN
    ! *** Communicate updated B values
    Call communicate_ghost_cells(B)
    IF( ISTRAN(2) > 0 .AND. ISICE == 4 )THEN
      Call communicate_ghost_cells(RHOW)
    ENDIF
  ENDIF
  
  Call Communicate_1D2(TBX, TBY)

  TMPITMP = DSTIME(0) - TTDS
  DSITIMING(9) = DSITIMING(9) + TMPITMP

  TTDS = DSTIME(0)
  IF( KC > 1 )THEN
    Call Communicate_QQ
  ELSE
    Call Communicate_3D_0(QQ)
  ENDIF

  DSITIMING(7) = DSITIMING(7) + (DSTIME(0)-TTDS)
#endif

  ! *** *******************************************************************!
  ! *** *******************************************************************!
  ! *** HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
  ! *** *******************************************************************!
  ! *** *******************************************************************!
  !if( process_id == master_id )THEN
  ! *** WRITE TO TIME SERIES FILES
  IF( ISDYNSTP == 0 )THEN
    CTIM=DT*FLOAT(N)+TCON*TBEGIN
    CTIM=CTIM/TCON
  ELSE
    CTIM=TIMESEC/TCON
  ENDIF
  IF( ISTMSR >= 1 )THEN
    IF( N >= NBTMSR .AND. N <= NSTMSR )THEN
      IF( NCTMSR >= NWTMSR )THEN

        if( process_id == master_id )THEN
          CALL TMSR !@todo add Map_TMSR routines so this will work on all processes?
        endif
        NDIFF = NWTMSR - NCTMSR
        ICALLTP = 1
        NCTMSR = NINCRMT + NDIFF
      ELSE
        NCTMSR = NCTMSR + NINCRMT
      ENDIF
    ENDIF
  ENDIF

  ! *** *******************************************************************!
  ! *** OUTPUT ZERO DIMENSION VOLUME BALANCE
  if( process_id == master_id )THEN
    IF( ISDRY >= 1 .AND. ISDRY < 98 )THEN
      IF( ICALLTP == 1 .AND. DEBUG )THEN
        OPEN(1,FILE=OUTDIR//'ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
        DO LS=1,LORMAX
          IF( VOLZERD >= VOLSEL(LS) .AND. VOLZERD < VOLSEL(LS+1) )THEN
            WTM=VOLSEL(LS+1)-VOLZERD
            WTMP=VOLZERD-VOLSEL(LS)
            DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
            WTM=WTM/DELVOL
            WTMP=WTMP/DELVOL
            SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
            ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
          ENDIF
        ENDDO
        IF( ISDYNSTP == 0 )THEN
          CTIM=DT*FLOAT(N)+TCON*TBEGIN
          CTIM=CTIM/TCTMSR
        ELSE
          CTIM=TIMESEC/TCTMSR
        ENDIF
        WRITE(1,5304) CTIM,SELZERD,ASFZERD,VOLZERD,VETZERD
        CLOSE(1)
      ENDIF
5304  FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
    ENDIF
  end if

  ICALLTP=0

  ! *** *******************************************************************!
  ! *** CALCULATE MEAN MASS TRANSPORT FIELD
  IF( (ISSSMMT > 0 .OR. ISWASP > 0).AND. RESSTEP > 0 ) CALL CALMMT

  ! *** *******************************************************************!
  ! *** ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
  IF( ISPD > 0 )THEN
    IF( TIMEDAY >= LA_BEGTI0 .AND. TIMEDAY <= LA_ENDTI0 )THEN
      TTDS = DSTIME(0)
      CALL DRIFTER_CALC
      TLRPD = TLRPD + (DSTIME(0)-TTDS)
    ENDIF
  ENDIF

  ! *** CALL TWO-TIME LEVEL BALANCES
  IF( ISBAL >= 1 ) CALL BAL2T5

  ! *** *******************************************************************!
  ! *** PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS
  IF( ISHTA == 1 ) CALL CALHTA

  ! *** *******************************************************************!
  ! *** CALCULATE DISPERSION COEFFICIENTS
  IF( N >= NDISP )THEN
    IF( N >= NDISP .AND. NCTBC == 1 )THEN
      IF( ISDISP == 2 ) CALL CALDISP2
      IF( ISDISP == 3 ) CALL CALDISP3
    ENDIF
  ENDIF

  ! *** *******************************************************************!
  ! *** PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS
  IF( ISLSHA == 1 .AND. TIMESEC >= TIMEHARM )THEN
    CALL LSQHARM  !@TODO remove write statements
    TIMEHARM = TIMEHARM + 300.
  ENDIF

  ! *** *******************************************************************!
  ! *** WRITE TO TIME VARYING GRAPHICS FILES
  if( process_id == master_id )THEN
    !@TODO Ned to have an allgather statement here to write out files at EE interval
    IF( HFREOUT == 1 )THEN
      DO NS=1,NSUBSET
        DEL = ABS(84600*(TIMEDAY-HFREDAY(NS)))
        IF( DEL <= DELT .AND. TIMEDAY <= HFREDAYEN(NS) )THEN
#ifdef NCOUT          
          IF( NCDFOUT > 0 ) call nc_output_hf(0,ns)
#endif 
          CALL HFREHYOUT(0,NS)
          CALL HFREWCOUT(0,NS)
            IF( ISTRAN(8) >= 1 ) THEN
              CALL HFREWQOUT(0,NS)
              IF( ISRPEM > 0 ) CALL HFRERPEMOUT(0,NS)
            ENDIF
          HFREDAY(NS) = HFREDAY(NS) + HFREMIN(NS)/1440
        ENDIF
      ENDDO
    ENDIF
  end if

#ifdef NCOUT
    ! ** NETCDF OUTPUT
    IF( NCDFOUT > 0 )THEN
      IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) .AND. (TIMEDAY >= TBEGNCDF .AND. TIMEDAY <= TENDNCDF ) )THEN
        Call Map_Write_NetCDF
        
        if( process_id == master_id )THEN
            CALL nc_output()
        end if
      ENDIF
    ENDIF

#endif  

  ! *** *******************************************************************!
  ! *** WRITE EFDC EXPLORER FORMAT OUTPUT
  IF( ISPPH == 1 )THEN
    ! *** CHECK IF PROPWASH OUTPUT FORCES EE LINKAGE TO BE WRITTEN
    if( ISPROPWASH > 0 )then
      ! *** Determine if PW_Mesh was writting on any process
      CALL DSI_All_Reduce(iwrite_pwmesh, i, MPI_SUM, TTDS, 0, TWAIT)

      IF( i > 0 )THEN
        IF( TIMEDAY >= last_snapshot + freq_out_min/1440. )THEN
          iwrite_pwmesh = 2
        ELSE
          iwrite_pwmesh = 0
        ENDIF
      ENDIF
    ENDIF
    
    !CALL BCOUT_ACCUMULATE(DTDYN) @todo invoke and make work in MPI case
    IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )THEN
      TTDS = DSTIME(0)
      Call Map_Write_EE_Binary
      TMPIEE = TMPIEE + (DSTIME(0)-TTDS)
      
      if( process_id == master_id )THEN
        CALL EE_LINKAGE(0)
      endif
     IF( iwrite_pwmesh > 0 ) last_snapshot = TIMEDAY
    ENDIF
  ENDIF
  IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )THEN
    NSNAPSHOTS = NSNAPSHOTS + 1
  ENDIF
  
  ! *** *******************************************************************!
  ! *** WRITE TO TIME VARYING 3D HDF GRAPHICS FILES
  if( process_id == master_id )THEN
    IF( N == NC3DO .AND. IS3DO == 1 )THEN
      CALL OUT3D
      NC3DO=NC3DO+(NTSPTC/NP3DO)
    ENDIF
  end if

  ! *** *******************************************************************!
  ! *** WRITE RESTART FILE EVERY ISRESTO REFERENCE PERIODS
  IF( ISRESTO >= 1 )THEN
    IF( TIMEDAY >= TIMELAST )THEN
      CALL Restart_Out(0)
      IF( ISTRAN(8) >= 1 )THEN
        CALL WQ_WCRST_OUT
        IF( IWQBEN == 1 ) CALL WQSDRST_OUT
        IF( ISRPEM > 0  ) CALL WQRPEMRST_Out
      ENDIF
      TIMELAST = TIMEDAY + TIMERST
    ENDIF
  ENDIF

  ! *** *******************************************************************!
  ! *** RECORD TIME
  if( process_id == master_id )THEN
    IF( TIMEDAY >= DAYNEXT )THEN
      TTDS = DSTIME(0) - T1TMP
      CALL TIMELOG(N,TIMEDAY,OUTDIR, TTDS)
      
      DAYNEXT = DAYNEXT + 1.
    ENDIF
  end if

  
  ! *** *******************************************************************!
  if( process_id == master_id )THEN
    IF( ISHOW > 0 ) CALL SHOWVAL
  end if
  ! *** *******************************************************************!
  
#ifdef _WIN  
  if( num_Processors == 1 )THEN     ! *** OMP Run
    IF( KEY_PRESSED() )THEN
      IF( ISEXIT() )THEN
        TTDS = DSTIME(0)
        Call Map_Write_EE_Binary
        TMPIEE = TMPIEE + (DSTIME(0)-TTDS)
      
        CALL EE_LINKAGE(-1)
        GOTO 1000
      ENDIF
    ENDIF
  endif
#endif
  
  GOTO 1001

  ! *** *******************************************************************!
  ! *** *******************************************************************!
  ! *** TIME LOOP COMPLETED
  ! *** *******************************************************************!
  ! *** *******************************************************************!
1000 CONTINUE
  THDMT = THDMT + (DSTIME(0)-T1TMP)

  ! *** *******************************************************************!
  ! *** WRITE RESTART FILE
  IF( ABS(ISRESTO) > 0 .AND. TIMEDAY > TIMELAST-TIMERST*0.5 )THEN
    if( process_id == master_id ) WRITE(6,'(A,F12.4)') 'FINAL RESTART FILE: ',TIMEDAY
    CALL Restart_Out(0)
    IF( ISTRAN(8) >= 1 )THEN
      CALL WQ_WCRST_OUT
      IF( IWQBEN == 1 ) CALL WQSDRST_OUT
      IF( ISRPEM > 0  ) CALL WQRPEMRST_Out
    ENDIF
  ENDIF

  ! *** *******************************************************************!
  ! *** COMPLETE LEAST SQUARES HARMONIC ANALYSIS
  LSLSHA=1
  IF( ISLSHA == 1 ) CALL LSQHARM

  ! *** *******************************************************************!
  ! *** OUTPUT COURANT NUMBER DIAGNOSTICS
  if( process_id == master_id )THEN
    IF( ISINWV == 1 .AND. DEBUG )THEN
      OPEN(1,FILE=OUTDIR//'CFLMAX.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'CFLMAX.OUT')

      DO L=2,LA
        WRITE(1,1991) IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
        WRITE(1,1992) (CFLVVV(L,K),K=1,KC)
        WRITE(1,1992) (CFLWWW(L,K),K=1,KC)
        WRITE(1,1992) (CFLCAC(L,K),K=1,KC)
      ENDDO

      CLOSE(1)
1991  FORMAT(2I5,12F8.3)
1992  FORMAT(10X,12F8.3)
    ENDIF
    ! *** *******************************************************************!
    ! *** OUTPUT COSMETIC VOLUME LOSSES FORM DRY CELLS
    IF( NDRYSTP > 0 )THEN
      OPEN(1,FILE=OUTDIR//'DRYLOSS.OUT')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'DRYLOSS.OUT')

      WRITE(1,'(3A5,A14,A10)')'L','I','J','VOL WASTED','EQUIV HP'
      DO L=2,LA
        WRITE(1,1993)L,IL(L),JL(L),VDWASTE(L),VDWASTE(L)/DXYP(L)
      ENDDO
1993  FORMAT(3I5,E14.6,F10.3)

      CLOSE(1)
    ENDIF

    ! *** *******************************************************************!
    ! *** OUTPUT FINAL FOOD CHAIN AVERAGING PERIOD
    IF( ISTRAN(5) >= 1 .AND. ISFDCH >= 1 )CALL FOODCHAIN(1)

  end if
  ! *** *******************************************************************!
  ! *** OUTPUT FINAL MASS AND VOLUME BALANCES
  IF( ISBAL >= 1 )CALL BAL2T5

  if( process_id == master_id )THEN
    ! *** *******************************************************************!
    CLOSE(90)
    CLOSE(98)

  end if

END
