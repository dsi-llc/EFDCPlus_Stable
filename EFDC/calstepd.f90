! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CALSTEPD

  ! CHANGE RECORD
  ! *** SUBROUTINE CALSTEP ESTIMATE THE CURRENT MAXIMUM TIME STEP SIZE
  ! *** FORM LINEAR STABILITY CRITERIA AND A FACTOR OF SAFETY
  !
  !----------------------------------------------------------------------C
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP

  Use GLOBAL
  Use Allocate_Initialize
  
#ifdef _MPI  
  Use Variables_MPI
  USE MPI
  Use Variables_MPI_Mapping
  Use MPI_All_Reduce
#endif

  IMPLICIT NONE

  INTEGER :: ND, L, K, LF, LL, LE, LN, LP, KM, LLOC, ITRNTMP, NX, NDYN
  INTEGER :: NMD, LMDCHHT, LMDCHUT, LMDCHVT, ITMPR, LLOCOLD,  MINTYPE
  INTEGER,SAVE :: NUP
  !INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LDYN
  
  REAL      :: QUKTMP, QVKTMP, RTMPR, TESTTEMP, DTMAXX
  REAL      :: TMPUUU, DTTMP, TMPVVV, TMPVAL, THP1, THP2
  REAL      :: TOP, QXPLUS, QYPLUS, QZPLUS, QXMINS, QYMINS, QZMINS, QTOTAL, QSRC, BOT
  REAL      :: DTCOMP, DTWARN, DTDYNP
  REAL,SAVE :: HPLIM, HPLIMOLD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL2
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL3
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL4
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QSUBINN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QSUBOUT
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTHISTORY

  REAL(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS, TWAIT                 ! MODEL TIMING TEMPORARY VARIABLE
  real(rkd) :: dtold
  
  IF( .NOT. ALLOCATED(DTL1) )THEN
    Call AllocateDSI(DTL1,    LCM,  2.*TIDALP)
    Call AllocateDSI(DTL2,    LCM,  2.*TIDALP)
    Call AllocateDSI(DTL3,    LCM,  2.*TIDALP)
    Call AllocateDSI(DTL4,    LCM,  2.*TIDALP)
    Call AllocateDSI(QSUBINN, LCM,  KCM,  0.0)
    Call AllocateDSI(QSUBOUT, LCM,  KCM,  0.0)
    !Call AllocateDSI(LDYN,    LCM,  0)

    ALLOCATE(DTHISTORY(10))
   
    DTHISTORY = DT
    NUP = 0
  ENDIF

  ! ********************************************************************************
  
  ! *** SET WATER COLUMN CONSTITUENT TRANSPORT FLAG
  ITRNTMP=0
  DO NX=1,7
    ITRNTMP = ITRNTMP + ISTRAN(NX)
  ENDDO
  
  IF( NITER <= 1 ) DTDYN = DT

  DTMIN = DT
  DTMAXX = TIDALP

  ! *** DETERMINE SOURCE/SINKS FOR SUBGRID SCALE CHANNEL EXCHANGES
  IF( MDCHH >= 1 )THEN
    DO K=1,KC
      DO L=2,LA
        QSUBOUT(L,K)=0.0
        QSUBINN(L,K)=0.0
      ENDDO
    ENDDO

    DO K = 1,KC
      DO NMD = 1,MDCHH
        IF( LKSZ(L,K) )CYCLE
        LMDCHHT = LMDCHH(NMD)
        LMDCHUT = LMDCHU(NMD)
        LMDCHVT = LMDCHV(NMD)
        IF( MDCHTYP(NMD) == 1 )THEN
          QUKTMP = QCHANU(NMD)*DZC(LMDCHUT,K)
          QVKTMP = 0.
        ENDIF
        IF( MDCHTYP(NMD) == 2 )THEN
          QVKTMP = QCHANV(NMD)*DZC(LMDCHVT,K)
          QUKTMP = 0.
        ENDIF
        IF( MDCHTYP(NMD) == 3 )THEN
          QUKTMP = QCHANU(NMD)*DZC(LMDCHUT,K)
          QVKTMP = QCHANV(NMD)*DZC(LMDCHVT,K)
        ENDIF
        QSUBOUT(LMDCHHT,K) = QSUBOUT(LMDCHHT,K) + MIN(QUKTMP,0.) + MIN(QVKTMP,0.)
        QSUBINN(LMDCHHT,K) = QSUBINN(LMDCHHT,K) + MAX(QUKTMP,0.) + MAX(QVKTMP,0.)
        QSUBOUT(LMDCHUT,K) = QSUBOUT(LMDCHUT,K) - MAX(QUKTMP,0.)
        QSUBINN(LMDCHUT,K) = QSUBINN(LMDCHUT,K) - MIN(QUKTMP,0.)
        QSUBOUT(LMDCHVT,K) = QSUBOUT(LMDCHVT,K) - MAX(QVKTMP,0.)
        QSUBINN(LMDCHVT,K) = QSUBINN(LMDCHVT,K) - MIN(QVKTMP,0.)
      ENDDO
    ENDDO
  ENDIF

  ! *** Initialize
  DTL1 = DTMAXX
  IF( ITRNTMP > 0 )   DTL2 = DTMAXX
  IF( DTSSDHDT > 0. ) DTL4 = DTMAXX
  
  !$OMP PARALLEL DO DEFAULT(NONE)                                                               &
  !$OMP  SHARED(NDM, LDM, LA, LAWET, LWET, KC, ITRNTMP, LLWET, LKWET, LEC, LNC, IsGhost)        &
  !$OMP  SHARED(DTSSDHDT, DTMAXX, DXYP, DZC, UHDY2, VHDX2, W2, QSUM, DTL1, DTL2, DTL3, DTL4)    &
  !$OMP  SHARED(HDRY, HP, H1P, UHE, VHE, HUI, HVI, DXIU, DYIV, QSUBINN, QSUBOUT, DTDYN)         &
  !$OMP  PRIVATE(ND, L, K, LF, LL, LP, LE, LN, KM)                                              &
  !$OMP  PRIVATE(TMPUUU, DTTMP, TMPVVV, TMPVAL, TESTTEMP, TOP, QXPLUS, QYPLUS, QZPLUS)          &
  !$OMP  PRIVATE(QXMINS, QYMINS, QZMINS, QTOTAL, QSRC, BOT, THP1, THP2)
  DO ND=1,NDM
    LF = 2+(ND-1)*LDM
    LL = MIN(LF+LDM-1,LA)

    ! *** METHOD 1: COURANT–FRIEDRICHS–LEWY
    DO LP = 1,LAWET
      L = LWET(LP)   
      TMPVAL = 1.E-16
      TMPUUU = ABS(UHE(L))*HUI(L)*DXIU(L)
      TMPVVV = ABS(VHE(L))*HVI(L)*DYIV(L)
      TMPVAL = TMPUUU + TMPVVV
      IF( TMPVAL > 0.  ) DTL1(L)  = 1./TMPVAL
    ENDDO

    ! *** METHOD 2: POSITIVITY OF ADVECTED MATERIAL, DTL2
    IF( ITRNTMP >= 1 )THEN
      DO K=1,KC
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IF( IsGhost(L) ) CYCLE     ! *** Do not use ghost cells to set min DT
          
          LE = LEC(L)
          LN = LNC(L)
          KM = K - 1
          
          ! *** Volume of layer
          TOP    = H1P(L)*DXYP(L)*DZC(L,K)
          
          ! *** Positive Fluxes
          QXPLUS = UHDY2(LE,K)
          QXPLUS = MAX(QXPLUS,0.0)
          QYPLUS = VHDX2(LN,K)
          QYPLUS = MAX(QYPLUS,0.0)
          QZPLUS = W2(L,K)*DXYP(L)
          QZPLUS = MAX(QZPLUS,0.0)
          
          ! *** Negative fluxes
          QXMINS = UHDY2(L,K)
          QXMINS = -MIN(QXMINS,0.0)
          QYMINS = VHDX2(L,K)
          QYMINS = -MIN(QYMINS,0.0)
          QZMINS = W2(L,KM)*DXYP(L)
          QZMINS = -MIN(QZMINS,0.0)
          
          QTOTAL = QSUM(L,K) + QSUBOUT(L,K) + QSUBINN(L,K)
          QSRC   = -MIN(QTOTAL,0.0)
          BOT = QXPLUS + QYPLUS + QZPLUS + QXMINS + QYMINS + QZMINS + QSRC
          
          IF( BOT > 1.E-12 )THEN
            DTTMP = TOP/BOT
            DTL2(L) = MIN(DTL2(L),DTTMP)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! *** METHOD 3: IMPLICIT BOTTOM FRICTION AND ROTATIONAL ACCELERATION DAMPING
    !DO L=LF,LL
    !  IF( LACTIVE(L) )THEN
    !    TMPVAL=SUB(L)+SUB(LEC(L))+SVB(L)+SVB(LNC(L))
    !    IF( TMPVAL > 0.5 )THEN
    !      LN=LNC(L)
    !      TAUBC=QQ(L,0)/CTURB2
    !      UCTR=0.5*(U(L,1)+U(LEC(L),1))
    !      VCTR=0.5*(V(L,1)+V(LN,1))
    !      UHMAG=HP(L)*SQRT(UCTR*UCTR+VCTR*VCTR)
    !      IF( UHMAG > 0.0 )THEN
    !        FRIFRE=TAUBC/UHMAG
    !        FRIFRE2=FRIFRE*FRIFRE
    !        ACACTMP=(CAC(L,KC)*HPI(L)*DXYIP(L))**2
    !        IF( ACACTMP > FRIFRE2 )THEN
    !          DTTMP=2.*FRIFRE/(ACACTMP-FRIFRE2)
    !          DTL3(L)=MIN(DTL3(L),DTTMP)
    !        ENDIF
    !      ENDIF
    !    ENDIF
    !  ENDIF
    !ENDDO

    ! ***  METHOD 4: LIMIT RATE OF DEPTH CHANGE
    IF( DTSSDHDT > 0. )THEN
      DO LP = 1,LLWET(KC,ND)
        L = LKWET(LP,KC,ND)  
        IF( IsGhost(L) ) CYCLE     ! *** Do not use ghost cells to set min DT
        
        THP1 = HP(L) 
        IF( THP1 < HDRY/10.) THP1 = HDRY/10.
        THP2 = H1P(L) 
        IF( THP2 < HDRY/10.) THP2 = HDRY/10.
        TESTTEMP = MAX(ABS(THP1-THP2),1.E-06)
        TMPVAL   = DTDYN*HP(L)/TESTTEMP
        DTL4(L)  = DTSSDHDT*TMPVAL
      ENDDO
    ENDIF
  ENDDO   ! *** END OF DOMAIN
  !$OMP END PARALLEL DO

  ! *** CHOOSE THE MINIMUM OF THE THREE METHODS
  MINTYPE = 0

  L1LOC = MINLOC(DTL1,DIM=1)
  DTL1MN = DTL1(L1LOC)

  IF( ITRNTMP >= 1 )THEN
    L2LOC = MINLOC(DTL2,DIM=1)
    DTL2MN = DTL2(L2LOC)
  ENDIF

  !DTL3MN = 2.*DTMAXX
  !DO L=2,LA
  !  IF( DTL3MN > DTL3(L) )THEN
  !    DTL3MN = DTL3(L)
  !    L3LOC = L
  !  ENDIF
  !ENDDO

  DTL4MN = 2.*DTMAXX
  IF( DTSSDHDT > 0. )THEN
    L4LOC = MINLOC(DTL4,DIM=1)
    DTL4MN = DTL4(L4LOC)
  ENDIF

  ! *** FIND MINIMUM & APPLY A SAFETY FACTOR
  DTL1MN = DTL1MN*DTSSFAC
  DTTMP  = 2.*DTMAXX               ! *** DTTMP - MINIMUM TIME STEP WITHOUT SAFETY FACTOR
  IF( DTTMP > DTL1MN )THEN
    DTTMP  = DTL1MN
    DTCOMP = 1./DTSSFAC
    LLOC   = L1LOC
    MINTYPE = 1
  ENDIF

  IF( ITRNTMP >= 1 )THEN
    DTL2MN = DTL2MN*0.5             ! *** FIXED SAFETY FACTOR FOR ADVECTION OF 0.5
    IF( DTTMP > DTL2MN )THEN
      DTTMP  = DTL2MN
      DTCOMP = 2.
      LLOC   = L2LOC
      MINTYPE = 2
    ENDIF
  ENDIF
  
  IF( DTSSDHDT > 0. )THEN
    IF( DTTMP > DTL4MN )THEN
      DTTMP  = DTL4MN
      DTCOMP = 1.0
      LLOC   = L4LOC
      MINTYPE = 4
    ENDIF
  ENDIF
  
  LLOCOLD = LMINSTEP              ! *** Global Old LMINSTEP

  ! *** Synchronize all processes with minimum timestep and corresponding global LMINSTEP
  TMPVAL = MIN(DTTMP,DTMAX)
  LL = MAP2GLOBAL(LLOC).LG        ! *** Global LMINSTEP for current process

  CALL DSI_All_Reduce(TMPVAL, LL, DTTMP, LMINSTEP, MPI_MIN, TTDS, 1, TWAIT)
  DSITIMING(11) = DSITIMING(11) + TTDS

  DTCOMP = DTTMP*DTCOMP           ! *** Minimum delta T without any safety factors
  IF( DTCOMP < DTMIN )THEN
    WRITE(6,800) TIMEDAY, DTTMP, DTMIN, Map2Global(LLOC).IG, Map2Global(LLOC).JG, HP(LLOC)
    WRITE(6,801) Map2Global(L1LOC).IG, Map2Global(L1LOC).JG, DTL1MN, HP(L1LOC)
    WRITE(6,802) Map2Global(L2LOC).IG, Map2Global(L2LOC).JG, DTL2MN, HP(L2LOC)
    IF( DTSSDHDT > 0. )THEN
      WRITE(6,804) Map2Global(L4LOC).IG, Map2Global(L4LOC).JG, DTL4MN
    ENDIF

    if( process_id == master_id )THEN
      OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
      WRITE(8,800) TIMEDAY, DTTMP, DTMIN, Map2Global(LLOC).IG, Map2Global(LLOC).JG, HP(LLOC)
      WRITE(8,801) Map2Global(L1LOC).IG, Map2Global(L1LOC).JG, DTL1MN, HP(L1LOC)
      WRITE(8,802) Map2Global(L2LOC).IG, Map2Global(L2LOC).JG, DTL2MN, HP(L2LOC)
      !WRITE(8,803)IL(L3LOC),JL(L3LOC),DTL3MN,HP(L3LOC)
      !WRITE(6,803)IL(L3LOC),JL(L3LOC),DTL3MN,HP(L3LOC)
      IF( DTSSDHDT > 0. )THEN
        WRITE(8,804) Map2Global(L4LOC).IG, Map2Global(L4LOC).JG, DTL4MN
      ENDIF
      CLOSE(8)
    End if !***End calculation on master process
    DTTMP = DTMIN

  ELSEIF( DTTMP < DTMIN )THEN
    DTWARN = DTTMP
    DTTMP  = DTMIN
  ELSE
    TMPVAL = DTTMP/DTMIN
    ITMPR = NINT(TMPVAL)
    RTMPR = FLOAT(ITMPR)
    IF( RTMPR < TMPVAL )THEN
      DTTMP = RTMPR*DTMIN
    ELSE
      DTTMP = (RTMPR-1.)*DTMIN
    ENDIF
  ENDIF

  ! *** FORCE TO MINIMUM TIME STEP ON STARTUP
  IF( NITER < 2 ) DTTMP = DTMIN

  ! *** RESTRICT INCREASE IN TIME STEP TO DTMIN
  !DTDYNP = DTHISTORY(10) + DTMIN
  DTDYNP = DTDYN + DTMIN               ! *** AT THIS POINT DTDYN IS THE PREVIOUS ITERATIONS TIME STEP
  IF( DTTMP > DTDYNP )THEN
    NUP = NUP + 1
    IF( NUP >= NUPSTEP )THEN
      DTTMP = DTDYNP
      NUP = 0
    ELSE
      DTTMP = DTDYN
    ENDIF
  ELSEIF( DTTMP < DTDYN )THEN
    IF( (NUPSTEP >= 25 .AND. TIMEDAY > TBEGIN+1.) .OR. 2.*DTTMP < DTDYN )THEN
      if( process_id == master_id )THEN
        ! *** REPORT WHEN THERE IS A DECREASE
        IF( 2.*DTTMP < DTDYN )THEN
          WRITE(6,'(A,I10,F14.5,2(A,F10.4,I7))') 'DTDYN REDUCED',NITER,TIMEDAY,'     [DT,L]   OLD: ',DTDYN,LLOCOLD,'     NEW: ',DTTMP,LMINSTEP
        ENDIF
        
        OPEN(9,FILE=OUTDIR//'TIME.LOG',POSITION='APPEND')
        WRITE(9,'(A,I10,F14.5,2(A,F10.4,I7,F8.3))') 'DYNAMIC STEP DOWN',NITER,TIMEDAY,'     [DT,L,HP]   OLD: ',DTDYN,LLOCOLD,-99.99,'     NEW: ',DTTMP,LMINSTEP,-99.99
        CLOSE(9)
      end if
    ENDIF

    NUP = 0
  ELSE
    DTTMP = DTDYN
  ENDIF
  DTDYN = MIN(DTTMP,DTMAX)

  ! *** SET INCREMENTAL INCREASE IN OUTPUT COUNTER
  NINCRMT = NINT(DTDYN/DTMIN)
  DTDYN   = FLOAT(NINCRMT)*DTMIN

  !IF( ITRNTMP == 0 ) DTL2MN = DTDYN

100 FORMAT(5I5,5F12.5,E13.5)
101 FORMAT(3I5,E13.5)
800 FORMAT('  TIME,DTDYN,DTMIN,I,J,HP = ',F12.5,2E12.4,2I7,F10.3)
801 FORMAT('  MOM  ADV,I,J,DT,HP = ',2I7,E13.4,F10.3)
802 FORMAT('  MASS ADV,I,J,DT,HP = ',2I7,E13.4,F10.3)
803 FORMAT('  CURV ACC,I,J,DT,HP = ',2I7,E13.4,F10.3)
804 FORMAT('  LIM DHDT,I,J,DT,HP = ',2I7,E13.4,F10.3)
880 FORMAT(3I5,8E13.4)
8899 FORMAT(' DT3 ERROR ',2I5,6E13.5)

  !DO K=10,2,-1
  !  DTHISTORY(K) = DTHISTORY(K-1)
  !ENDDO
  !DTHISTORY(1) = DTDYN
     
  ! *** *******************************************************************C
  RETURN
END
