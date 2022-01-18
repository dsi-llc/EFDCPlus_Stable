! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CALQVS(ISTL_)

  ! *** SUBROUTINE CALQVS UPDATES TIME VARIABLE VOLUME SOURCES

  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !--------------------------------------------------------------------------------------!
  !    2018-03       PAUL M. CRAIG     IMPLEMENTED WHOLE CHANNEL ELEVATION RATING CURVE (NQCTYP=2)
  !    2018-02       N T LAM           IMPLEMENTED TIME VARIABLE GATE STRUCTURES
  !    2015-11       D H CHUNG  &      IMPLEMENTED NEW HYDRAULIC STRUCTURES BASED ON DILL
  !                  PAUL M. CRAIG
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2015-01       PAUL M. CRAIG     ADDED FULLY COUPLED ICE SUB-MODEL
  !                  DANG H CHUNG
  !    2014-09       PAUL M. CRAIG     IMPLEMENTED NEW LOW CHORD BC TYPE
  !    2014-08       D H CHUNG         SET EXPLICIT PRECISIONS of INTEGER & REAL

  USE GLOBAL
  Use Variables_WQ
  USE HEAT_MODULE, ONLY:ICECOMP
  USE HYDSTRUCMOD
  USE DIFFUSER_MODULE, ONLY:JPEFDC
  USE FIELDS

  Use Variables_MPI
  Use Variables_MPI_Mapping

  IMPLICIT NONE

  INTEGER :: ISTL_
  INTEGER :: L, K, LL, NS, NCTMP, M1, M2, ITYP, NCTL, NCTLT, I, NC, IG
  INTEGER :: IU, JU, ID, JD, LU, LD, ND, LF, IPMC, MU1, MU2, MD1, MD2
  INTEGER :: NTMP, NWR, KU, KD, NJP, LJP, KTMP, ITMPD, NTT, LEVELFLAG, MWET
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:)   :: FLOWTYPE,  NUNIQQSER,  NWET
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:, :) :: UNIQQSER

  REAL :: QWRABS, QSERTCELL, CTIM, TDIFF, WTM1, WTM2, HUP, HDW, DELH, HDRY9, TF, HPICE
  REAL :: TMPVAL, QCTLMAX, HTMPD, TDIFFU, WTM1U, WTM2U, TDIFFD, WTM1D, WTM2D, CLEVAPTMP
  REAL :: TVW, TVA, DTV, DTVL, LAMBDA, TM, VPTG                                           !< RYAN-HARLEMAN
  REAL :: RPORTS, QVJPTMP, RAVAIL, QSUMIET, QEAVAIL, DIFQVOL, DTAGW, RIFTRL, QSUMTMP
  REAL :: QSUM1, QSUM2                                                                    !< FLOW ACCUMULATORS TO HANDLE SPECIFIED FLOWS BELOW ACTIVE BOTTOM LAYER (SGZ)
  REAL :: ZHU, ZHD, HVAL, CVAL, RVAL, DENAIR, FWUP, QSW, QAA

  REAL(RKD) :: SPLIT, QLAYER1, QLAYER2, DDELT
  
  REAL,SAVE :: ICESTEP,DELTICE1,DELTICE2,ICEDAY

  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QLAYER
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: QSFACTTOT
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: QSFACT
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: QSUM3                                !< FLOW ACCUMULATOR FOR IN/OUT FLOWS FOR A DRY CELL
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:)   :: QSUM4                                !< FLOW ACCUMULATOR FOR IN/OUT FLOWS FOR A DRY CELL

  LOGICAL,SAVE :: LGROUPS

  IF( .NOT. ALLOCATED(QLAYER) )THEN
    ! ** FIRST CALL  FOR TYP=3/4
    ALLOCATE(QLAYER(NQCTL,KCM))
    ALLOCATE(USCELL(NQCTL))
    ALLOCATE(DSCELL(NQCTL))
    ALLOCATE(QSFACTTOT(NQCTL))
    ALLOCATE(QSFACT(NQCTL))
    IF( NGRPID > 0 .AND. NQSIJ > 0 )THEN
      ALLOCATE(QSUM3(2,NGRPID))
      ALLOCATE(QSUM4(NGRPID))
      ALLOCATE(FLOWTYPE(NGRPID))
      ALLOCATE(NUNIQQSER(NGRPID))
      ALLOCATE(NWET(NGRPID))
      ALLOCATE(UNIQQSER(NQSIJ,NGRPID))

      FLOWTYPE = 0
      DO LL=1,NQSIJ
        FLOWTYPE(GRPID(LL)) = NQSMUL(LL)
      ENDDO

      ! *** LIST OF UNIQUE QSER SERIES FOR EACH GROUP
      UNIQQSER = 0
      NUNIQQSER = 0
      DO LL=1,NQSIJ
        DO I=1,NUNIQQSER(GRPID(LL))
          IF( NQSERQ(LL) == UNIQQSER(I,GRPID(LL)) ) EXIT
        ENDDO
        IF( I > NUNIQQSER(GRPID(LL)) )THEN
          NUNIQQSER(GRPID(LL)) = NUNIQQSER(GRPID(LL)) + 1
          UNIQQSER(NUNIQQSER(GRPID(LL)),GRPID(LL)) = NQSERQ(LL)
        ENDIF
      ENDDO

      NWET = 0
      QSUM3 = 0.
      QSUM4 = 0.
    ENDIF

    QLAYER=0.
    QSFACTTOT=0.
    QSFACT=0.

    LEVELFLAG = 0
    DO NCTL=1,NQCTL
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      USCELL(NCTL) = LU
      IF( ISRESTI == 0 )THEN
        SAVESUB(1,NCTL) = SUBO(LU)
        SAVESVB(1,NCTL) = SVBO(LU)
      ENDIF

      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      DSCELL(NCTL) = LC
      IF( ID > 0 .AND. JD > 0 )THEN
        LD = LIJ(ID,JD)
        DSCELL(NCTL) = LD
        IF( ISRESTI == 0 )THEN
          SAVESUB(2,NCTL) = SUBO(LD)
          SAVESVB(2,NCTL) = SVBO(LD)
        ENDIF
      ELSE
        LD = LC
      ENDIF

      ! *** RESET IF USING RESTART
      IF( ISRESTI /= 0 )THEN
        HRUO(LU) = SAVESUB(1,NCTL)*DYU(LU)*DXIU(LU)
        HRVO(LU) = SAVESVB(1,NCTL)*DXV(LU)*DYIV(LU)
        HRUO(LD) = SAVESUB(2,NCTL)*DYU(LD)*DXIU(LD)
        HRVO(LD) = SAVESVB(2,NCTL)*DXV(LD)*DYIV(LD)
      ENDIF

      IF( HYD_STR(NCTL).NQCTYP == 3 .OR. HYD_STR(NCTL).NQCTYP == 4 )THEN
        ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE DIFFERENCE (4) DEPENDANT FLOWS
        NCTLT = HYD_STR(NCTL).NQCTLQ  ! *** SERIES INDEX

        ! *** ADJUST THE FLOWS TO ZERO AT LOWER CHORD DEPTHS
        IF( HYD_STR(NCTL).NQCTYP == 3 )THEN
          DELH = HCTLUM(NCTLT)*(HYD_STR(NCTL).BQCLCE-BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU)

          ! *** LOOKUP THE FLOWS
          M1=0
          M2=1
        70 M1=M1+1
          M2=M2+1
          IF( M2 > MQCTL(NCTLT) )THEN

            if( process_id == master_id )THEN
              WRITE(6,*)' BAD QCTL BOUNDARY SPECIFICATION.  LOW CHORD OUT OF RANGE OF RATING CURVE.'
              WRITE(6,6667)NCTL,NCTLT,IU,JU,ID,JD
              WRITE(6,*)' LOW CHORD:  ',HYD_STR(NCTL).BQCLCE,DELH
              OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
              WRITE(8,*)' BAD QCTL BOUNDARY SPECIFICATION.  LOW CHORD OUT OF RANGE OF RATING CURVE.'
              WRITE(8,6667)NCTL,NCTLT,IU,JU,ID,JD
              WRITE(8,*)' LOW CHORD:  ',HYD_STR(NCTL).BQCLCE,DELH
              CLOSE(8)
              CALL STOPP('.')
            endif

          ENDIF
          IF( DELH >= HDIFCTL(M1,NCTLT) .AND. DELH <= HDIFCTL(M2,NCTLT) )THEN
            ! *** FOUND LOWER CHORD ELEVATION, DETERMINE THE FLOW OFFSET
            TDIFF = HDIFCTL(M2,NCTLT)-HDIFCTL(M1,NCTLT)
            WTM1 = (HDIFCTL(M2,NCTLT)-DELH)/TDIFF
            WTM2 = (DELH-HDIFCTL(M1,NCTLT))/TDIFF
            DO K=1,KC
              QLAYER(NCTL,K) = WTM1*QCTL(M1,1,K,NCTLT)+WTM2*QCTL(M2,1,K,NCTLT)
            ENDDO

            CYCLE
          ELSE
            GOTO 70
          ENDIF

          ! *** INITIALIZE IF WSEL>LOWCHORD
          IF( ISRESTI == 0 )THEN
            HUP=BELV(LU)+HP(LU)
            IF( HUP <= HYD_STR(NCTL).BQCLCE )NLOWCHORD(NCTL) = -(HYD_STR(NCTL).NQCMINS+1)
          ENDIF

        ENDIF
      ENDIF

    ENDDO

    LGROUPS = .FALSE.
    DO NCTL=1,NQCTL
      IF( HYD_STR(NCTL).NQCTYP == -2 )THEN
        LGROUPS = .TRUE.
        IF( HYD_STR(NCTL).QCTLGRP < 1  )THEN
          WRITE(6,'(" BAD GROUPID AT NCTL: ",I5,I10)') NCTL,HYD_STR(NCTL).QCTLGRP
          CALL STOPP('.')
        ENDIF
      ENDIF
    ENDDO

    ! *** ICE SETTINGS
    ICESTEP = 60.*15.  ! *** DURATION TO ACHEIVE "QUASI-STEADY STATE" CONDITIONS FOR ICE SURFACE TEMPERATURE
    DELTICE1  = 0.0
    DELTICE2  = 0.0
    ICEDAY = INT(TIMEDAY)
  ENDIF

  IF( ISTL_ == 2 )THEN
    IF( ISDYNSTP == 0 )THEN
      DELT=DT
    ELSE
      DELT=DTDYN
    END IF
  ELSE
    DELT=DT2
  ENDIF
  DDELT = DELT
  HDRY9=0.9*HDRY

  ! **  INITIALIZE NULL (0) FLOW SERIES
  GWSERT(0)=0.
  QWRSERT(0)=0.
  QSERTCELL=0.0
  DO K=1,KC
    QSERT(K,0)=0.
    QCTLT(K,0,1)=0.
    QCTLT(K,0,2)=0.
    QCTLTO(K,0)=0.
  ENDDO

  ! ** 2018-10-08, NTL: ADD TIME VARIABLE RATING CURVE
  ! IMPLEMENTATION: RATING CURVES ARE SUDDENLY CHANGED BETWEEN TWO TIME STEPS
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP < 5 )THEN
      IF( HSCTL(NCTL).ITYPE == 1 )THEN
        ! *** STRUCTURE IS CONTROLLED BY TIME-SERIES
        ID = HSCTL(NCTL).ID
        CALL CHANGE_RATING_CURVE_BY_TIME(QCTLSER(ID), TIMESEC, HYD_STR(NCTL).NQCTLQ)
      ELSEIF (HSCTL(NCTL).ITYPE == 2 .OR. HSCTL(NCTL).ITYPE == 3 )THEN
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
        CALL CHANGE_RATING_CURVE_BY_RULES(NCTL, HYD_STR(NCTL).NQCTLQ)
      ELSE
        ! *** STRUCTURE IS UNCONTROLLED, DO NOTHING
      ENDIF
    ENDIF
  ENDDO

  !$OMP PARALLEL DEFAULT(SHARED)

  IF( NGWSER >= 1 )THEN
    !$OMP SINGLE
    NCTMP = 3 + NDYM + NSED + NSND + NTOX
    DO NC=1,NCTMP
      GWCSERT(0,NC) = 0.
    ENDDO
    !$OMP END SINGLE

    IF( ISTRAN(5) > 0 )THEN
      !$OMP DO PRIVATE(ND,LF,LL,L,NC)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO NC=1,NCTMP
          DO L=LF,LL
            CONGW(L,NC)=0.0
          ENDDO
        END DO
      END DO
      !$OMP END DO
    ENDIF
  ENDIF

  ! **  INITIALIZE TOTAL FLOW SERIES
  !$OMP DO PRIVATE(ND,LF,LL,L)
  DO ND=1,NDM
    LF=2+(ND-1)*LDM
    LL=MIN(LF+LDM-1,LA)
    DO L=LF,LL
      QSUM1E(L)=QSUME(L)
      QSUME(L)=0.
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** SELECTIVE ZEROING
  IF( KC > 1 )THEN
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      !$OMP DO PRIVATE(ND,LF,LL,L,K)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          QSUM(L,KSZ(L)) = 0.
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF

    ! *** ZERO EVAP/RAINFALL/ICE
    !$OMP DO PRIVATE(ND,LF,LL,L)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)
      DO L=LF,LL
        QSUM(L,KC) = 0.
      ENDDO
    ENDDO
    !$OMP END DO

    ! *** ZERO ALL DEFINED BC'S
    !$OMP SINGLE
    DO NS=1,NBCS
      L=LBCS(NS)
      DO K=1,KC
        QSUM(L,K) = 0.
      ENDDO
    ENDDO
    !$OMP END SINGLE

  ELSE
    ! *** SINGLE LAYER
    !$OMP DO PRIVATE(ND,LF,LL,L)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)
      DO L=LF,LL
        QSUM(L,KC) = 0.
      ENDDO
    ENDDO
    !$OMP END DO
  ENDIF
  !$OMP END PARALLEL

  ! *************************************************************************************
  ! *** VOLUME SOURCE/SINK INTERPOLATION
  DO NS=1,NQSER
    IF( ISTL_ == 2 )THEN
      IF( ISDYNSTP == 0 )THEN
        CTIM = DT*(FLOAT(N)-0.5)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
      ELSE
        CTIM = TIMESEC/TCQSER(NS)
      ENDIF
    ELSE
      IF( ISDYNSTP == 0 )THEN
        CTIM=DT*FLOAT(N-1)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
      ELSE
        CTIM=TIMESEC/TCQSER(NS)
      ENDIF
    ENDIF
    M2=MTSQLAST(NS)
    DO WHILE (CTIM > QSER(NS).TIM(M2))
      M2=M2+1
      IF( M2 > MQSER(NS) )THEN
        M2=MQSER(NS)
        EXIT
      ENDIF
    END DO
    MTSQLAST(NS)=M2
    M1 = M2-1
    TDIFF = QSER(NS).TIM(M2)-QSER(NS).TIM(M1)
    WTM1 = (QSER(NS).TIM(M2)-CTIM)/TDIFF
    WTM2 = (CTIM-QSER(NS).TIM(M1))/TDIFF
    DO K=1,KC
      QSERT(K,NS) = WTM1*QSER(NS).VAL(M1,K) + WTM2*QSER(NS).VAL(M2,K)
    ENDDO
  ENDDO

  IF( N == 1 )THEN
    if( process_id == master_id )THEN
      OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
      DO LL=1,NQSIJ
        L=LQS(LL)
        ITYP=LCT(L)
        IF( ITYP <= 0 .OR. ITYP >= 8 )THEN
          WRITE(6,6111)LL,IQS(LL),JQS(LL)
          WRITE(8,6111)LL,IQS(LL),JQS(LL)
        ENDIF
      ENDDO
      CLOSE(8)
    endif
  ENDIF

  DO LL=1,NQSIJ
    NS=NQSERQ(LL)
    L=LQS(LL)
    QSUM1=0.
    QSUM2=0.
    DO K=1,KC
      ! *** APPLY MULTIPLIERS HERE TO CORRECT MASS BALANCE PROBLEMS
      IF( NQSMUL(LL) < 5 .OR. ISDRY == 0 )THEN
        ! *** FIXED BC GROUP DISTRIBUTION BY QFACTOR
        QSS(K,LL)      = QSSE(LL)   *RQSMUL(LL)*DZC(L,K)
        QSERCELL(K,LL) = QSERT(K,NS)*RQSMUL(LL)*QFACTOR(LL)

      ELSEIF( NQSMUL(LL) == 5 )THEN
        ! *** REDISTRIBUTION OF DRY CELL FLOWS
        QSS(K,LL)      = QSSE(LL)*DZC(L,K)
        QSERCELL(K,LL) = QSERT(K,NS)*QFACTOR(LL)

        IF( HP(L) < HDRY )THEN
          QSUM3(1,GRPID(LL)) = QSUM3(1,GRPID(LL)) + QSS(K,LL)
          QSUM3(2,GRPID(LL)) = QSUM3(2,GRPID(LL)) + QSERCELL(K,LL)
          QSS(K,LL) = 0.
          QSERCELL(K,LL) = 0.
        ENDIF

      ELSEIF( NQSMUL(LL) == 6 )THEN
        QSS(K,LL) = 0.
        QSERCELL(K,LL) = 0.
        CYCLE
      ENDIF
      QSUM(L,K) = QSUM(L,K) + QSS(K,LL) + QSERCELL(K,LL)

      IF( LKSZ(L,K) )THEN
        ! *** LAYER IS INACTIVE
        QSUM1 = QSUM1 + QSS(K,LL)
        QSUM2 = QSUM2 + QSERCELL(K,LL)
        QSS(K,LL)      = 0.
        QSERCELL(K,LL) = 0.
        QSUM(L,K)      = 0.
      ENDIF
    ENDDO

    ! *** ADD FLOWS BELOW BOTTOM ACTIVE LAYER BACK TO ACTIVE LAYERS
    IF( QSUM1 /= 0. )THEN
      ! *** CONSTANT FLOW
      DO K = KSZ(L),KC
        QSS(K,LL) = QSS(K,LL) + QSUM1*DZC(L,K)
        QSUM(L,K) = QSUM(L,K) + QSUM1*DZC(L,K)
      ENDDO
    ENDIF

    IF( QSUM2 /= 0. )THEN
      ! *** TIME VARIABLE
      DO K = KSZ(L),KC
        QSERCELL(K,LL) = QSERCELL(K,LL) + QSUM2*DZC(L,K)
        QSUM(L,K)      = QSUM(L,K)      + QSUM2*DZC(L,K)
      ENDDO
    ENDIF
  ENDDO

  ! *** REDISTRIBUTE FLOWS INTO OTHER WET CELLS IN THE SAME GROUP
  IF( NQSIJ > 0 )THEN
    DO IG=1,NGRPID
      IF( FLOWTYPE(IG) == 5 )THEN
        IF( QSUM3(1,IG) /= 0. .OR. QSUM3(2,IG) /= 0. )THEN
          ! *** FOUND SOME FLOW, SO SPLIT FLOWS BY CELL VOLUME
          QSUM1 = 0.
          QSUM2 = 0.
          DO LL=1,NQSIJ
            IF( GRPID(LL) == IG )THEN
              ! *** FOUND A CELL IN THE GROUP.  GET VOLUME
              L = LQS(LL)
              IF( HP(L) > HDRY .OR. ISDRY == 0 )THEN
                QSUM1 = QSUM1 + HP(L)*DXYP(L)            ! *** TOTAL WET VOLUME
              ELSE
                QSUM2 = QSUM2 + HP(L)*DXYP(L)            ! *** TOTAL DRY VOLUME
              ENDIF
            ENDIF
          ENDDO

          ! *** NOW REDISTRIBUTE IF WET CELLS AVAILABLE
          DO LL=1,NQSIJ
            IF( GRPID(LL) == IG )THEN
              ! *** TEST IF THERE ARE ANY WET CELLS IN GROUP
              L = LQS(LL)
              IF( QSUM1 > 0. )THEN
                ! *** WET CELLS ARE AVAILABLE.  DISTRIBUTE FLOWS INTO GROUP
                IF( HP(L) > HDRY .OR. ISDRY == 0 )THEN
                  SPLIT = HP(L)*DXYP(L)/QSUM1                 ! *** FLOW SPLIT
                  DO K = KSZ(L),KC
                    QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                    QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                    QSS(K,LL)      = QSS(K,LL)      + QLAYER1
                    QSERCELL(K,LL) = QSERCELL(K,LL) + QLAYER2
                    QSUM(L,K)      = QSUM(L,K)      + QLAYER1 + QLAYER2
                  ENDDO
                ENDIF
              ELSEIF( QSUM2 > 0. )THEN
                ! *** NO WET CELLS AVAILABLE.  JUST USE THE ORIGINAL FLOWS
                SPLIT = HP(L)*DXYP(L)/QSUM2                 ! *** FLOW SPLIT
                DO K = KSZ(L),KC
                  QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                  QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                  QSS(K,LL)      = QSS(K,LL)      + QLAYER1
                  QSERCELL(K,LL) = QSERCELL(K,LL) + QLAYER2
                  QSUM(L,K)      = QSUM(L,K)      + QLAYER1 + QLAYER2
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      ELSEIF( FLOWTYPE(IG) == 6 )THEN
        ! *** REDISTRIBUTE TOTAL FLOW INTO WET CELLS BY VOLUME

        ! *** GET VOLUMES
        QSUM1 = 0.
        QSUM2 = 0.
        QSUM3(:,IG) = 0.
        MWET = 0

        ! *** GET TOTAL CONSTANT FLOW (QSUM3(1,IG)) INTO THIS GROUP 
        DO LL=1,NQSIJ
          IF( GRPID(LL) == IG )THEN
            ! *** FOUND A CELL IN THE GROUP.  GET VOLUME
            NS=NQSERQ(LL)
            L = LQS(LL)
            IF( HP(L) >= HDRY .OR. ISDRY == 0 )THEN
              QSUM1 = QSUM1 + HP(L)*DXYP(L)            ! *** TOTAL WET VOLUME
              MWET = MWET + 1
            ELSE
              QSUM2 = QSUM2 + HP(L)*DXYP(L)            ! *** TOTAL DRY VOLUME
            ENDIF
            QSUM3(1,IG) = QSUM3(1,IG) + QSSE(LL)       ! *** SUM CONSTANT FLOW
            DO K = 1,KC
              QSS(K,LL)      = 0.
              QSERCELL(K,LL) = 0.
            ENDDO
          ENDIF
        ENDDO

        ! *** GET TOTAL TIME SERIES FLOW (QSUM3(2,IG)) INTO THIS GROUP 
        DO I=1,NUNIQQSER(IG)
          NS = UNIQQSER(I,IG)
          DO K = 1,KC
            QSUM3(2,IG) = QSUM3(2,IG) + QSERT(K,NS)  ! *** SUM TIME SERIES FLOWS
          ENDDO
        ENDDO

        ! *** NOW REDISTRIBUTE IF WET CELLS AVAILABLE
        IU = 0
        DO LL=1,NQSIJ
          IF( GRPID(LL) == IG )THEN
            ! *** TEST IF THERE ARE ANY WET CELLS IN GROUP
            L = LQS(LL)

            ! *** TEST IF UPDATE IS NEEDED
            IF( QSUM1 > 0. )THEN
              IF( MWET == NWET(IG) .AND. ABS(QSUM1-QSUM4(IG))/QSUM1 < 0.10 )THEN
                ! *** USE PREVIOUS TIME STEPS SPLITS IF LITTLE HAS CHANGED SINCE THE LAST UPDATED SPLITS
                SPLIT = QFACTOR(LL)
                DO K = KSZ(L),KC
                  QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                  QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                  QSS(K,LL)      = QLAYER1
                  QSERCELL(K,LL) = QLAYER2
                  QSUM(L,K)      = QSUM(L,K) + QLAYER1 + QLAYER2
                ENDDO
              ELSE
                ! *** WET CELLS ARE AVAILABLE.  DISTRIBUTE FLOWS INTO GROUP
                IF( HP(L) >= HDRY .OR. ISDRY == 0 )THEN
                  SPLIT = HP(L)*DXYP(L)/QSUM1                 ! *** FLOW SPLIT
                  DO K = KSZ(L),KC
                    QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                    QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                    QSS(K,LL)      = QLAYER1
                    QSERCELL(K,LL) = QLAYER2
                    QSUM(L,K)      = QSUM(L,K) + QLAYER1 + QLAYER2
                  ENDDO
                  QFACTOR(LL) = SPLIT
                ELSE
                  QFACTOR(LL) = 0.0
                ENDIF
                IU = 1
              ENDIF

            ELSEIF( QSUM2 > 0. )THEN
              ! *** NO WET CELLS AVAILABLE.  JUST USE THE DRY VOLUMES
              SPLIT = HP(L)*DXYP(L)/QSUM2                 ! *** FLOW SPLIT
              DO K = KSZ(L),KC
                QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                QSS(K,LL)      = QLAYER1
                QSERCELL(K,LL) = QLAYER2
                QSUM(L,K)      = QSUM(L,K) + QLAYER1 + QLAYER2
              ENDDO
              QFACTOR(LL) = SPLIT
            ENDIF
          ENDIF
        ENDDO

        ! *** REPORT WHEN CHANGING FLOW SPLITS
        IF( NITER == 1 .OR. IU == 1 )THEN
          PRINT '(A,F10.4,I8,3I5,2E14.6)', 'Change in flow splits: ',TIMEDAY, NITER, IG, NWET(IG), MWET, QSUM4(IG), QSUM1
        ENDIF

        ! *** SAVE CURRENT WET CELL CONDITIONS
        NWET(IG) = MWET
        IF( IU == 1 ) QSUM4(IG) = QSUM1

      ENDIF
    ENDDO
  ENDIF

  ! *************************************************************************************
  ! **  GROUNDWATER SOURCE/SINK INTERPOLATION
  IF( NGWSER >= 1 .AND. GWSP.IFLAG == 0 )THEN
    NCTMP = 3 + NDYM + NSED + NSND + NTOX
    DO NS=1,NGWSER
      IF( ISTL_ == 2 )THEN
        IF( ISDYNSTP == 0 )THEN
          CTIM=DT*(FLOAT(N)-0.5)/TCGWSER(NS) + TBEGIN*(TCON/TCGWSER(NS))
        ELSE
          CTIM=TIMESEC/TCGWSER(NS)
        ENDIF
      ELSE
        IF( ISDYNSTP == 0 )THEN
          CTIM=DT*FLOAT(N-1)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
        ELSE
          CTIM=TIMESEC/TCGWSER(NS)
        ENDIF
      ENDIF
      M2 = MTSGWLAST(NS)
      DO WHILE ( CTIM > TGWSER(M2,NS) )
        M2 = M2+1
        IF( M2 > MGWSER(NS) )THEN
          M2 = MGWSER(NS)
          EXIT
        ENDIF
      END DO
      MTSGWLAST(NS) = M2
      M1 = M2-1

      TDIFF = TGWSER(M2,NS) - TGWSER(M1,NS)
      WTM1 = (TGWSER(M2,NS)-CTIM)/TDIFF
      WTM2 = (CTIM-TGWSER(M1,NS))/TDIFF
      GWSERT(NS) = WTM1*GWSER(M1,NS) + WTM2*GWSER(M2,NS)
      DO NC=1,NCTMP
        GWCSERT(NS,NC) = WTM1*GWCSER(M1,NS,NC) + WTM2*GWCSER(M2,NS,NC)
        IF( GWCSERT(NS,NC) < 0.0 )THEN
          PRINT '(A,I10,F10.4,2I6)','GWCSERT(NS,NC) is < 0.0!', niter, timeday, ns, nc  
          PAUSE
        ENDIF
      END DO
    ENDDO

  ENDIF

  ! *************************************************************************************
  ! *** CONTROL STRUCTURES AND TIDAL INLETS

  ! -------------------------------------------------------------------------------------
  ! *** SET CURRENT TIME QSFACTOR BASED ON ELEVATIONS
  IF( LGROUPS )THEN
    QSFACTTOT=0.
    DO NCTL=1,NQCTL
      IF( HYD_STR(NCTL).NQCTYP == -2 )THEN
        IU = HYD_STR(NCTL).IQCTLU
        JU = HYD_STR(NCTL).JQCTLU
        LU = LIJ(IU,JU)
        IF( HP(LU) >= HDRY )THEN
          QSFACTTOT(HYD_STR(NCTL).QCTLGRP) = QSFACTTOT(HYD_STR(NCTL).QCTLGRP) + HYD_STR(NCTL).QCTLMU
        ENDIF
      ENDIF
    ENDDO

    ! *** NOW COMPUTE CURRENT TIME STEP FLOW SPLITS
    DO NCTL=1,NQCTL
      IF( HYD_STR(NCTL).NQCTYP == -2 )THEN
        IF( QSFACTTOT(HYD_STR(NCTL).QCTLGRP) > 0. )THEN
          QSFACT(NCTL) = HYD_STR(NCTL).QCTLMU/QSFACTTOT(HYD_STR(NCTL).QCTLGRP)
        ELSE
          QSFACT(NCTL) = 0.
        ENDIF
      ENDIF
    ENDDO

    ! *** UPSTREAM ELEVATION DEPENDANT FLOWS
    DO NCTL=1,NQCTL
      IF( HYD_STR(NCTL).NQCTYP == -2 )THEN
        NCTLT = HYD_STR(NCTL).NQCTLQ  ! *** SERIES INDEX
        IU = HYD_STR(NCTL).IQCTLU
        JU = HYD_STR(NCTL).JQCTLU
        LU = LIJ(IU,JU)
        HUP = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
        HUP = HUP*HCTLUM(NCTLT)

        IF( HP(LU) < HWET )THEN
          DO K=1,KC
            QCTLT(K,NCTL,1)=0.
          ENDDO
        ELSE
          ! *** SUFFICIENT DEPTH, GET FLOWS
          M1=0
          M2=1
      500 M1=M1+1
          M2=M2+1
          IF( M2 > MQCTL(NCTLT) )THEN
            WRITE(6,*) ' US ELEVATION CONTROL STRUCTURE TABLE OUT OF BOUNDS'
            WRITE(6,'(" NCTL,NCTLT,IU,JU = ",4I5)') NCTL,NCTLT,IU,JU
            WRITE(6,'(" WSEL,HU = ",2F10.3)') HUP,HP(LU)
            OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
            WRITE(8,*) ' US ELEVATION CONTROL STRUCTURE TABLE OUT OF BOUNDS'
            WRITE(8,'(" NCTL,NCTLT,IU,JU = ",4I5)') NCTL,NCTLT,IU,JU
            WRITE(8,'(" WSEL,HU = ",2F10.3)') HUP,HP(LU)
            CLOSE(8)
            CALL STOPP('.')
          ENDIF
          IF( HUP >= HDIFCTL(M1,NCTLT) .AND. HUP <= HDIFCTL(M2,NCTLT) )THEN
            TDIFF = HDIFCTL(M2,NCTLT) - HDIFCTL(M1,NCTLT)
            WTM1 = (HDIFCTL(M2,NCTLT) - HUP)/TDIFF
            WTM2 = (HUP - HDIFCTL(M1,NCTLT))/TDIFF
            DO K=1,KC
              QCTLT(K,NCTL,1) = WTM1*QCTL(M1,1,K,NCTLT) + WTM2*QCTL(M2,1,K,NCTLT)
              QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*QSFACT(NCTL)
            ENDDO
          ELSEIF( HUP < HDIFCTL(1,NCTLT) )THEN
            ! *** Water levels below rating curve.  Assume zero flow
            DO K = 1,KC
              QCTLT(K,NCTL,1) = 0.
            ENDDO
          ELSE
            GOTO 500
          ENDIF
        ENDIF

        ! *** APPLY RQCMUL FLOW SCALING
        DO K=1,KC
          QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*HYD_STR(NCTL).RQCMUL
        ENDDO
      ENDIF   ! **** HYD_STR(NCTL).NQCTYP == -2
    ENDDO     ! *** END OF NQCTL LOOP
  ENDIF

  ! -------------------------------------------------------------------------------------
  ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP >= -1 .AND. HYD_STR(NCTL).NQCTYP <= 1 )THEN
      ! *** UPSTREAM ELEVATION/STAGE (0 & 1) OR DEPTH (-1) DEPENDANT FLOWS
      NCTLT=HYD_STR(NCTL).NQCTLQ  ! *** SERIES INDEX

      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      HUP = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
      IF( HYD_STR(NCTL).NQCTYP == -1 ) HUP = HP(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU

      ID=HYD_STR(NCTL).IQCTLD
      JD=HYD_STR(NCTL).JQCTLD
      IF( ID <= 1 .OR. JD <= 1 )THEN
        LD = 1
        HDW = 0.
      ELSE
        LD = LIJ(ID,JD)
        HDW = HP(LD) + BELV(LD) + HCTLDA(NCTLT) + HYD_STR(NCTL).HQCTLD
        IF( HYD_STR(NCTL).NQCTYP == -1 )THEN
          ! *** UPSTREAM DEPTH ONLY
          HDW = 0.
        ENDIF
      ENDIF

      DELH = HCTLUM(NCTLT)*HUP - HCTLDM(NCTLT)*HDW
      IF( HYD_STR(NCTL).NQCTYP == 0 .AND. AQCTL(NCTLT) > 0.0 )THEN
        IF( HUP < AQCTL(NCTLT) ) DELH = -100.
      ENDIF
      IF( DELH <= 0. .OR. HP(LU) < HWET )THEN
        ! *** PREVENT REVERSE FLOW OR INSUFFICENT DEPTH
        DO K=1,KC
          QCTLT(K,NCTL,1) = 0.
        ENDDO
      ELSE
        ! *** SUFFICIENT DEPTH
        IF( HYD_STR(NCTL).NQCTYP == 1 ) DELH=SQRT(DELH)
        M1=0
        M2=1
    600 M1=M1+1
        M2=M2+1
        IF( M2 > MQCTL(NCTLT) )THEN
          WRITE(6,6666)
          WRITE(6,6667)NCTL,NCTLT,IU,JU,ID,JD
          WRITE(6,6668)HUP,HP(LU),HDW,HP(LD)
          OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
          WRITE(8,6666)
          WRITE(8,6667)NCTL,NCTLT,IU,JU,ID,JD
          WRITE(8,6668)HUP,HP(LU),HDW,HP(LD)
          CLOSE(8)
          CALL STOPP('.')
        ENDIF
        IF( DELH >= HDIFCTL(M1,NCTLT) .AND. DELH <= HDIFCTL(M2,NCTLT) )THEN
          TDIFF = HDIFCTL(M2,NCTLT) - HDIFCTL(M1,NCTLT)
          WTM1 = (HDIFCTL(M2,NCTLT) - DELH)/TDIFF
          WTM2 = (DELH - HDIFCTL(M1,NCTLT))/TDIFF
          DO K=1,KC
            QCTLT(K,NCTL,1) = WTM1*QCTL(M1,1,K,NCTLT) + WTM2*QCTL(M2,1,K,NCTLT)
          ENDDO
        ELSE
          GOTO 600
        ENDIF
      ENDIF

      IF( HYD_STR(NCTL).NQCTYP == 1 )THEN
        ! *** ACCELERATION TYPE
        IF( ISTL_ == 3 )THEN
          DO K=1,KC
            QCTLST(K,NCTL)=QCTLT(K,NCTL,1)
            TMPVAL=QCTLTO(K,NCTL) + DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
            QCTLT(K,NCTL,1)=TMPVAL/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
            QCTLTO(K,NCTL)=QCTLT(K,NCTL,1)
            QCTLSTO(K,NCTL)=QCTLST(K,NCTL)
          ENDDO
        ELSE
          DO K=1,KC
            QCTLST(K,NCTL)=QCTLT(K,NCTL,1)
            TMPVAL=QCTLTO(K,NCTL) + DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
            QCTLT(K,NCTL,1)=TMPVAL/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
            QCTLT(K,NCTL,1)=0.5*(QCTLT(K,NCTL,1)+QCTLTO(K,NCTL))
          ENDDO
        ENDIF
      ENDIF

      ! *** APPLY RQCMUL FLOW SCALING
      DO K=1,KC
        QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*HYD_STR(NCTL).RQCMUL
      ENDDO

    ENDIF
  ENDDO

  ! -------------------------------------------------------------------------------------
  ! *** UPSTREAM/DOWNSTREAM ELEVATION RATING CURVE  HYD_STR(NCTL).NQCTYP == 2
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP == 2 )THEN
      NCTLT=HYD_STR(NCTL).NQCTLQ

      ! *** UPSTREAM ELEVATION
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      HUP = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU

      ! *** CHECK FOR ACTIVE FLOW CONDITIONS
      IF( HUP < HDIFCTL(1,NCTLT) .OR. HP(LU) < HWET )THEN
        ! *** INSUFFICIENT DEPTHS
        DO K=1,KC
          QCTLT(K,NCTL,1) = 0.
        ENDDO
        LD = 1
      ELSE
        ! *** SUFFICIENT DEPTHS
        ! *** DOWNSTREAM ELEVATION
        ID=HYD_STR(NCTL).IQCTLD
        JD=HYD_STR(NCTL).JQCTLD
        LD=LIJ(ID,JD)
        HDW = HP(LD) + BELV(LD) + HCTLDA(NCTLT) + HYD_STR(NCTL).HQCTLD
        HTMPD=HDIFCTD(1,NCTLT)+0.001
        HDW=MAX(HDW,HTMPD)

        MU1=0
        MU2=1
        MD1=0
        MD2=1

        ! *** FIND UPSTREAM ELEVATION BRACKET
    700 MU1=MU1+1
        MU2=MU1+1
        IF( MU2 > MQCTL(NCTLT) )THEN
          WRITE(6,6676)
          WRITE(6,6677)NCTL,NCTLT,IU,JU,ID,JD
          WRITE(6,6678)HUP,HP(LU),HDW,HP(LD)
          WRITE(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
          WRITE(8,6676)
          WRITE(8,6677)NCTL,NCTLT,IU,JU,ID,JD
          WRITE(8,6678)HUP,HP(LU),HDW,HP(LD)
          WRITE(8,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          CLOSE(8)
          CALL STOPP('.')
        ENDIF
        IF( HUP >= HDIFCTL(MU1,NCTLT) .AND. HUP <= HDIFCTL(MU2,NCTLT) )THEN
          ! *** FOUND VALID UPSTREAM ELEVATION RANGE
          TDIFFU=HDIFCTL(MU2,NCTLT)-HDIFCTL(MU1,NCTLT)
          WTM1U=(HDIFCTL(MU2,NCTLT)-HUP)/TDIFFU
          WTM2U=(HUP-HDIFCTL(MU1,NCTLT))/TDIFFU
        ELSE
          GOTO 700
        ENDIF

        ! *** FIND DOWNSTREAM ELEVATION BRACKET
    750 MD1=MD1+1
        MD2=MD1+1
        IF( MD2 > MQCTL(NCTLT) )THEN
          WRITE(6,6686)
          WRITE(6,6687)NCTL,NCTLT,IU,JU,ID,JD
          WRITE(6,6688)HUP,HP(LU),HDW,HP(LD)
          WRITE(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
          WRITE(8,6686)
          WRITE(8,6687)NCTL,NCTLT,IU,JU,ID,JD
          WRITE(8,6688)HUP,HP(LU),HDW,HP(LD)
          WRITE(8,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          CLOSE(8)
          CALL STOPP('.')
        ENDIF
        IF( HDW >= HDIFCTD(MD1,NCTLT) .AND. HDW <= HDIFCTD(MD2,NCTLT) )THEN
          ! *** FOUND VALID DOWNSTREAM ELEVATION RANGE
          TDIFFD=HDIFCTD(MD2,NCTLT)-HDIFCTD(MD1,NCTLT)
          WTM1D=(HDIFCTD(MD2,NCTLT)-HDW)/TDIFFD
          WTM2D=(HDW-HDIFCTD(MD1,NCTLT))/TDIFFD
        ELSE
          GOTO 750
        ENDIF

        ! *** DETERMINE FLOWS BASED ON UPSTREAM AND DOWNSTREAM ELEVATIONS
        DO K=1,KC
          QCTLT(K,NCTL,1) = WTM1U*( WTM1D*QCTL(MU1,MD1,K,NCTLT) + WTM2D*QCTL(MU1,MD2,K,NCTLT) ) &
            + WTM2U*( WTM1D*QCTL(MU2,MD1,K,NCTLT) + WTM2D*QCTL(MU2,MD2,K,NCTLT) )
        ENDDO
      ENDIF

      ! *** APPLY RQCMUL FLOW SCALING
      DO K=1,KC
        QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*HYD_STR(NCTL).RQCMUL
      ENDDO
    ENDIF
  ENDDO

  ! -------------------------------------------------------------------------------------
  ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE WITH A LOWER CHORD ACTIVATION TOGGLE
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP == 3 .OR. HYD_STR(NCTL).NQCTYP == 4 )THEN
      ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE DIFFERENCE (4) DEPENDANT FLOWS
      NCTLT=HYD_STR(NCTL).NQCTLQ

      IU=HYD_STR(NCTL).IQCTLU
      JU=HYD_STR(NCTL).JQCTLU
      LU=LIJ(IU,JU)

      ID=HYD_STR(NCTL).IQCTLD
      JD=HYD_STR(NCTL).JQCTLD
      IF( ID == 0 .AND. JD == 0 )THEN
        ! *** INVALID SPECIFICATION
        CYCLE
      ENDIF
      LD=LIJ(ID,JD)

      ! *** SKIP WHEN INSUFFICENT WATER
      IF( HP(LU) < HWET )CYCLE

      HUP = HP(LU) + BELV(LU) ! *** WSEL
      IF( HUP <= HYD_STR(NCTL).BQCLCE .AND. (N < HYD_STR(NCTL).NQCMINS .OR. NLOWCHORD(NCTL) < 0 .OR. NLOWCHORD(NCTL) > HYD_STR(NCTL).NQCMINS/2) )THEN
        ! *** WATER BELOW LOW CHORD

        ! *** SET THE CELL FACE FLAGS - ALLOW FULL HYDRODYNAMICS
        IF( LOWCHORDU(NCTL) /= -9999. )THEN
          if( process_id == master_id )THEN
            WRITE(*,*)' *** LOWER CHORD BC: OFF   NCTL=',NCTL

            OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
            WRITE(8,*)' *** LOWER CHORD BC: OFF   NCTL=',NCTL
            WRITE(8,*)'                           TIMEDAY=',TIMEDAY
            WRITE(8,*)'                           LU,LD=',LU,LD
            WRITE(8,*)'                           HUP,HDP=',HUP,HP(LD)+BELV(LD)
            WRITE(8,*)' ***'
            CLOSE(8)
          end if

          DO K=1,KC
            QCTLT(K,NCTL,1)=0.
          ENDDO
          IF( ID > IU )THEN
            IF( SAVESUB(2,NCTL)>0.5 )THEN
              SUB(LD)   = 1.0
              SUBO(LD)  = 1.0
              UHDYE(LD) = UHDYE(LU)
              TBX(LD)   = TBX(LU)
              TSX(LD)   = TSX(LU)
            ENDIF
          ELSE
            IF( SAVESUB(1,NCTL)>0.5 )THEN
              SUB(LU)   = 1.0
              SUBO(LU)  = 1.0
              UHDYE(LU)  = UHDYE(LD)
              TBX(LU)   = TBX(LD)
              TSX(LU)   = TSX(LD)
            ENDIF
          ENDIF
          IF( JD > JU )THEN
            IF( SAVESVB(2,NCTL)>0.5 )THEN
              SVB(LD)   = 1.0
              SVBO(LD)  = 1.0
              VHDXE(LD) = VHDXE(LU)
              TBY(LD)   = TBY(LU)
              TSY(LD)   = TSY(LU)
            ENDIF
          ELSE
            IF( SAVESVB(1,NCTL)>0.5 )THEN
              SVB(LU)=1.0
              SVBO(LU)=1.0
              VHDXE(LU) = VHDXE(LD)
              TBY(LU)   = TBY(LD)
              TSY(LU)   = TSY(LD)
            ENDIF
          ENDIF
          LOWCHORDU(NCTL) = -9999.
          LOWCHORDV(NCTL) = -9999.
          NLOWCHORD(NCTL)=0
          IPMC=0
        ENDIF
        NLOWCHORD(NCTL)=NLOWCHORD(NCTL)-1
        IPMC=0
        CYCLE

      ELSE

        ! *** ELEVATION ABOVE LOW CHORD
        IF( LOWCHORDU(NCTL) == -9999. )THEN
          IF( NLOWCHORD(NCTL) < -HYD_STR(NCTL).NQCMINS )THEN
            ! *** SET THE CELL FACE FLAGS - BLOCK FULL HYDRODYNAMICS

            WRITE(*,*)' *** LOWER CHORD BC: ON    NCTL=',NCTL

            IF( ID > IU )THEN
              SUB(LD) = 0.0
              SUBO(LD) = 0.0
              ! *** SAVE THE FLOWS
              LOWCHORDU(NCTL) = ABS(UHDYE(LD))
              UHDYE(LD)=0.0
            ELSEIF( ID < IU )THEN
              SUB(LU)=0.0
              SUBO(LU)=0.0
              ! *** SAVE THE FLOWS
              LOWCHORDU(NCTL) = ABS(UHDYE(LU))
              UHDYE(LU)=0.0
            ELSE
              LOWCHORDU(NCTL) = 0.
            ENDIF
            IF( JD > JU )THEN
              SVB(LD)=0.0
              SVBO(LD)=0.0
              ! *** SAVE THE FLOWS
              LOWCHORDV(NCTL) = ABS(VHDXE(LD))
              VHDXE(LD)=0.0
            ELSEIF( JD < JU )THEN
              SVB(LU)=0.0
              SVBO(LU)=0.0
              ! *** SAVE THE FLOWS
              LOWCHORDV(NCTL) = ABS(VHDXE(LU))
              VHDXE(LU)=0.0
            ELSE
              LOWCHORDV(NCTL) = 0.
            ENDIF
            if( process_id == master_id )THEN
              OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
              WRITE(8,*)' *** LOWER CHORD BC: ON    NCTL=',NCTL
              WRITE(8,*)'                           TIMEDAY=',TIMEDAY
              WRITE(8,*)'                           LU,LD=',LU,LD
              WRITE(8,*)'                           HUP,HDP=',HUP,HP(LD)+BELV(LD)
              WRITE(8,*)'                           QU,QV=',LOWCHORDU(NCTL),LOWCHORDV(NCTL)
              WRITE(8,*)' ***'
              CLOSE(8)
            endif

            NLOWCHORD(NCTL)=0
            IPMC=0
          ELSE
            NLOWCHORD(NCTL) = NLOWCHORD(NCTL)-1
            CYCLE
          ENDIF
        ENDIF

        IF( HYD_STR(NCTL).NQCTYP == 3 )THEN
          ! *** UPSTREAM DEPTH ONLY
          HUP  = HP(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
          DELH = HCTLUM(NCTLT)*HUP
        ELSE
          ! *** ELEVATION DIFFERENCE
          HUP  = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
          HDW  = HP(LD) + BELV(LD) + HCTLDA(NCTLT) + HYD_STR(NCTL).HQCTLD
          DELH = HCTLUM(NCTLT)*HUP - HCTLDM(NCTLT)*HDW
        ENDIF
        DELH=MAX( DELH,HDIFCTL(1,NCTLT) )

        ! *** LOOKUP THE FLOWS
        M1=0
        M2=1
    800 M1=M1+1
        M2=M2+1
        IF( M2 > MQCTL(NCTLT) )THEN
          if( process_id == master_id )THEN
            WRITE(6,6666)
            WRITE(6,6667)NCTL,NCTLT,IU,JU,ID,JD
            WRITE(6,6668)HUP,HP(LU),HDW,HP(LD)
            OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
            WRITE(8,6666)
            WRITE(8,6667)NCTL,NCTLT,IU,JU,ID,JD
            WRITE(8,6668)HUP,HP(LU),HDW,HP(LD)
            CLOSE(8)
            CALL STOPP('.')
          end if
        ENDIF
        IF( DELH >= HDIFCTL(M1,NCTLT) .AND. DELH <= HDIFCTL(M2,NCTLT) )THEN
          TDIFF = HDIFCTL(M2,NCTLT) - HDIFCTL(M1,NCTLT)
          WTM1 = (HDIFCTL(M2,NCTLT) - DELH)/TDIFF
          WTM2 = (DELH - HDIFCTL(M1,NCTLT))/TDIFF
          DO K=1,KC
            QCTLT(K,NCTL,1) = WTM1*QCTL(M1,1,K,NCTLT) + WTM2*QCTL(M2,1,K,NCTLT)
            ! *** SUBTRACT LOW CHORD FLOWS
            QCTLT(K,NCTL,1) = MAX(QCTLT(K,NCTL,1) - QLAYER(NCTL,K),0.)
          ENDDO
          DO K=1,KC
            ! ***                                 ** TOTAL FLOW OPEN CHANNEL FLOW  **
            QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1) + DZC(LU,K)*(LOWCHORDU(NCTL)+LOWCHORDV(NCTL))
          ENDDO
        ELSE
          GOTO 800
        ENDIF

        NLOWCHORD(NCTL)=NLOWCHORD(NCTL)+1
      ENDIF
    ENDIF
  ENDDO   ! *** END OF LOW/HI CHORD CHECK

  ! -------------------------------------------------------------------------------------
  ! *** COMPUTED FLOWS USING EQUATIONS OR DO FINAL UPDATES FOR OTHER TYPES
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP > 4 )THEN
      ! *** COMPUTED FLOWS USING HYDRAULIC STRUCTURE EQUATIONS:  HYD_STR(NCTL).NQCTYP > 4
      CALL COMPUTE_HSFLOW(NCTL)

    ELSE
      ! *** ADJUST FOR INACTIVE LAYERS AND ADD FINAL LAYERS FLOWS TO QSUM
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)

      ! *** CHECK VALID LAYERS
      QSUM1 = 0.
      DO K=1,KSZ(LU)-1
        QSUM1 = QSUM1 + QCTLT(K,NCTL,1)
        QCTLT(K,NCTL,1) = 0.
      ENDDO

      ! *** ADD FLOWS BELOW BOTTOM ACTIVE LAYER BACK TO ACTIVE LAYERS
      IF( QSUM1 > 0. )THEN
        DO K = KSZ(LU),KC
          QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1) + QSUM1*DZC(LU,K)
        ENDDO
      ENDIF

      ! *** LIMIT OUTFLOWS TO AVAILABLE WATER
      QCTLMAX = (HP(LU) - HDRY)*DXYP(LU)/(DELT*DZI)
      IF( QCTLMAX < 0. ) QCTLMAX = 0.
      DO K=1,KC
        QCTLT(K,NCTL,1) = MIN(QCTLT(K,NCTL,1),QCTLMAX)
      ENDDO

      ! *** FINAL ADJUSTED LAYER FLOWS ADDED TO TOTAL FLOWS
      DO K = KSZ(LU),KC
        QSUM(LU,K) = QSUM(LU,K) - QCTLT(K,NCTL,1)
      ENDDO

      ! *** DOWNSTREAM
      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      IF( ID /= 0 .AND. JD /= 0 )THEN
        LD = LIJ(ID,JD)
        QSUM2 = 0.
        DO K = 1,KSZ(LD)-1
          QSUM2 = QSUM2 + QCTLT(K,NCTL,1)
        ENDDO

        ! *** SET THE RETURN FLOW VARIABLE FOR LOOKUP TABLE CONTROL TYPE STRUCTURES
        DO K=KSZ(LD),KC
          QCTLT(K,NCTL,2) = QCTLT(K,NCTL,1) + QSUM2*DZC(LD,K)
          QSUM(LD,K) = QSUM(LD,K) + QCTLT(K,NCTL,2)
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  ! *************************************************************************************

  NTMP = 3 + NDYM + NTOX + NSED + NSND
  IF( ISTRAN(8) > 0 ) NTMP = NTMP + NWQV

  ! *************************************************************************************
  ! **  FLOW WITHDRAWAL AND RETURN
  DO NC=1,NTMP
    CQWRSERT(0,NC)=0.
  ENDDO
  DO NS=1,NQWRSR
    IF( ISTL_ == 2 )THEN
      IF( ISDYNSTP == 0 )THEN
        CTIM=DT*(FLOAT(N)-0.5)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
      ELSE
        CTIM=TIMESEC/TCQWRSR(NS)
      ENDIF
    ELSE
      IF( ISDYNSTP == 0 )THEN
        CTIM=DT*FLOAT(N-1)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
      ELSE
        CTIM=TIMESEC/TCQWRSR(NS)
      ENDIF
    ENDIF
    M2=MTSWRLAST(NS)
    DO WHILE (CTIM > TQWRSER(M2,NS))
      M2=M2+1
      IF( M2 > MQWRSR(NS) )THEN
        M2=MQWRSR(NS)
        EXIT
      ENDIF
    END DO
    MTSWRLAST(NS)=M2
    M1 = M2-1
    TDIFF=TQWRSER(M2,NS)-TQWRSER(M1,NS)
    WTM1=(TQWRSER(M2,NS)-CTIM)/TDIFF
    WTM2=(CTIM-TQWRSER(M1,NS))/TDIFF
    QWRSERT(NS)=WTM1*QWRSER(M1,NS)+WTM2*QWRSER(M2,NS)
    DO NC=1,NTMP
      CQWRSERT(NS,NC)=WTM1*CQWRSER(M1,NS,NC)+WTM2*CQWRSER(M2,NS,NC)
    ENDDO
  ENDDO

  IF( NQWR > 0 )THEN
    DO NWR=1,NQWR
      NS=WITH_RET(NWR).NQWRSERQ
      IF(WRCTL(NWR).ITYPE > 0 )THEN
        ! *** Withdrawal/Return is controlled by operation rules
        LU = LIJ(WRCTL(NWR).IREFUP, WRCTL(NWR).JREFUP)
        ZHU = BELV(LU) + HP(LU)

        IF( WRCTL(NWR).ITYPE == 2 )THEN
          ! *** W/R CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
          LD = LIJ(WRCTL(NWR).IREFDN, WRCTL(NWR).JREFDN)
          ZHD = BELV(LD) + HP(LD)
          HVAL = ZHU - ZHD
        ELSE
          ! *** W/R IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
          HVAL = ZHU
        ENDIF
        CVAL = WRCTL(NWR).CUR.FLOW
        CALL PUMP_OPERATION_RULES(WRCTL(NWR), HVAL, CVAL, RVAL)
        QWRSERT(NS) = RVAL
      ENDIF

      ! *** Handle +/- Flows for Withdrawal/Return Structures
      IF( QWRSERT(NS) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=WITH_RET(NWR).IQWRU
        JU=WITH_RET(NWR).JQWRU
        KU=WITH_RET(NWR).KQWRU
        ID=WITH_RET(NWR).IQWRD
        JD=WITH_RET(NWR).JQWRD
        KD=WITH_RET(NWR).KQWRD
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        ID=WITH_RET(NWR).IQWRU
        JD=WITH_RET(NWR).JQWRU
        KD=WITH_RET(NWR).KQWRU
        IU=WITH_RET(NWR).IQWRD
        JU=WITH_RET(NWR).JQWRD
        KU=WITH_RET(NWR).KQWRD
        WITH_RET(NWR).QWR=0.  ! *** Only allow time variable flows when using! -W/R
      ENDIF
      LU=LIJ(IU,JU)
      LD=LIJ(ID,JD)
      QWRABS = ABS(QWRSERT(NS))

      QSUM(LU,KU) = QSUM(LU,KU) - WITH_RET(NWR).QWR-QWRABS
      QSUM(LD,KD) = QSUM(LD,KD) + WITH_RET(NWR).QWR+QWRABS
      WRCTL(NWR).CUR.FLOW = QWRSERT(NS)
    ENDDO
  ENDIF

  ! **  CALL JPEFDC AND PLACE JET-PLUME VOLUMES SOURCES
  IF( NQJPIJ > 0 .AND. N == 1 ) CALL JPEFDC
  IF( NQJPIJ > 0 .AND. ISTL_ == 3 )THEN
    IF( NUDJPC(1) == NUDJP(1) )THEN
      CALL JPEFDC
      NUDJPC(1)=1
    ELSE
      NUDJPC(1)=NUDJPC(1)+1
    ENDIF
  ENDIF
  IF( NQJPIJ > 0 .AND. IS2TIM >= 1 )THEN
    IF( NUDJPC(1) == NUDJP(1) )THEN
      CALL JPEFDC
      NUDJPC(1)=1
    ELSE
      NUDJPC(1)=NUDJPC(1)+1
    ENDIF
  ENDIF

  ! *** PLACE JET-PLUME VOLUME SOURCES
  IF( NQJPIJ > 0 )THEN
    DO NJP=1,NQJPIJ
      IF( ICALJP(NJP) == 1 )THEN
        ! *** DISCHARGE ONLY USING QSER AND STANDARD CONCENTRATION SERIES
        RPORTS = FLOAT(NPORTJP(NJP))
        LJP = LIJ(IQJP(NJP),JQJP(NJP))
        KTMP = KEFFJP(NJP)

        ! QVJPTMP = JETPLUME DISCHARGE PER PORT
        ! QQCJP   = CONSTANT FLOW PER CELL
        ! QJPENT  = AMBIENT WATER ENTRAINMENT VOLUMETRIC RATE
        ! QJPENTT = AMBIENT WATER ENTRAINMENT VOLUMETRIC RATE PLUS PLUME VOLUME
        QVJPTMP  = QQCJP(NJP)
        DO K=1,KC
          QVJPTMP = QVJPTMP + QSERT(K,NQSERJP(NJP))
        ENDDO

        ! *** SUBTRACT THE ENTRAINMENT FROM EACH LAYER
        DO K=1,KC
          QSUM(LJP,K) = QSUM(LJP,K) - RPORTS*QJPENT(K,NJP)
        ENDDO

        ! *** PLACE DISCHARGE AND TOTAL ENTRAINMENT AT EFFECTIVE LAYER
        QSUM(LJP,KTMP) = QSUM(LJP,KTMP) + RPORTS*(QVJPTMP+QJPENTT(NJP))
      ENDIF
      IF( ICALJP(NJP) == 2 )THEN
        ! *** WITHDRAWAL/RETURN TYPE USING W/R SERIES
        RPORTS = FLOAT(NPORTJP(NJP))
        LJP = LIJ(IQJP(NJP),JQJP(NJP))
        KTMP = KEFFJP(NJP)

        ! QVJPTMP = JETPLUME DISCHARGE PER PORT
        QVJPTMP = QWRCJP(NJP) + QWRSERT(NQWRSERJP(NJP))

        ! SUBTRACT ENTRAIMENT FROM EACH LAYER
        DO K=1,KC
          QSUM(LJP,K) = QSUM(LJP,K) - RPORTS*QJPENT(K,NJP)
        ENDDO

        ! PLACE DISCHARGE AND TOTAL ENTRAINMENT AT EFFECTIVE LOCATION
        QSUM(LJP,KTMP) = QSUM(LJP,KTMP) + RPORTS*(QVJPTMP+QJPENTT(NJP))

        ! REMOVE DISCHARGE FROM UPSTREAM INTAKE CELL
        LU = LIJ(IUPCJP(NJP),JUPCJP(NJP))
        KU = KUPCJP(NJP)
        QSUM(LU,KU) = QSUM(LU,KU) - RPORTS*QVJPTMP
      ENDIF
    ENDDO
  ENDIF

  ! *** COMPUTE ICE MET/GROWTH VOLUMETRIC RATES (M3/S)
  IF( ISICE > 2 .AND. N > 0  .AND. ( LCHECKICE .OR. LFRAZIL ) )THEN
    DELTICE1 = DELTICE1 + DELT
    DELTICE2 = DELTICE2 + DELT

    CALL ICECOMP(ISTL_, DELTICE2, ICESTEP)
    IF( DELTICE2 >= ICESTEP ) DELTICE2 = 0.0  ! *** RESET ICE COVER TIMER AND VOLUME ACCUMULATORy
  ENDIF

  ! *** UPDATE TIME VARIABLE GROUNDWATER FIELD (MUST INCLUDE ALL CELLS)
  IF( GWSP.IFLAG > 0 )THEN
    CALL UPDATEFIELD(GWSP,TIMEDAY,1,QGW)
  ENDIF

  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** DETERMIME GROUNDWATER FLUXES
  IF( ISGWIT == 2 .AND. GWSP.IFLAG == 0 )THEN
    IF( NGWSER == 1 )THEN
      IF( IGWSER(1) == 0 )THEN
        ! *** FLOW SERIES
        !$OMP DO PRIVATE(ND,LF,LL,L)
        DO ND=1,NDM
          LF=2+(ND-1)*LDM
          LL=MIN(LF+LDM-1,LA)
          DO L=LF,LL
            QGW(L) = GWFAC(L)*GWSERT(NGWSL(L))
          END DO
        ENDDO
        !$OMP END DO
      ELSE
        ! *** SEEPAGE VELOCITY SERIES
        !$OMP DO PRIVATE(ND,LF,LL,L)
        DO ND=1,NDM
          LF=2+(ND-1)*LDM
          LL=MIN(LF+LDM-1,LA)
          DO L=LF,LL
            QGW(L) = GWFAC(L)*DXYP(L)*GWSERT(NGWSL(L))
          ENDDO
        ENDDO
        !$OMP END DO
      ENDIF
    ELSE
      !$OMP DO PRIVATE(ND,LF,LL,L)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          IF( IGWSER(NGWSL(L)) == 0 )THEN
            ! *** FLOW SERIES
            QGW(L) = GWFAC(L)*GWSERT(NGWSL(L))
          ELSE
            ! *** SEEPAGE VELOCITY SERIES
            QGW(L) = GWFAC(L)*DXYP(L)*GWSERT(NGWSL(L))
          ENDIF
        END DO
      ENDDO
      !$OMP END DO
    ENDIF
  ENDIF

  ! *** ACCUMULATE GROUNDWATER FLUXES
  IF( ISGWIT > 0 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,L)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)
      DO L=LF,LL
        QSUM(L,KSZ(L)) = QSUM(L,KSZ(L)) + QGW(L)
      ENDDO
    ENDDO
    !$OMP END DO
  ENDIF

  ! *** ADD ICE MET/GROWTH VOLUMETRIC RATES (M3/S)
  IF( ISICE > 2 .AND. ( LCHECKICE .OR. LFRAZIL ) .AND. DELTICE1 >= 60. )THEN
    !$OMP DO PRIVATE(ND,LF,LL,L,HPICE)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)
      DO L=LF,LL

        ! *** DUE TO MACHINE PRECISION ISSUES, ONLY APPLY ICE MET/FREEZE FLOW RATES WHEN LARGE ENOUGH
        HPICE = ICEVOL(L)*DXYIP(L)
        IF( ABS(HPICE) > HP(L)*1.E-6 )THEN
          ! *** LIMIT ICE GROWTH TO AVAILABLE WATER
          IF( HPICE < 0. .AND. ISDRY > 0 )THEN
            IF( HP(L)+HPICE < HDRYICE )THEN
              IF( ICEDAY /= INT(TIMEDAY) )THEN  ! *** ONLY DISPLAY WARNING ONCE PER DAY
                WRITE(6,'(A,I10,I5,3F10.5)') 'LIMITING ICE PRODUCTION (NITER,L,HP,HPICE): ',NITER,L,HP(L),HPICE,HP(L)+HPICE
                ICEDAY = INT(TIMEDAY)
              ENDIF
              ICEVOL(L) = 0.0
            ENDIF
          ENDIF

          ICERATE(L) = ICEVOL(L)/DELT
          QSUM(L,KC) = QSUM(L,KC) + ICERATE(L)
          ICEVOL(L)  = 0.0
        ELSEIF( ICETHICK(L) < HPICE )THEN
          ICERATE(L) = ICEVOL(L)/DELT
          QSUM(L,KC) = QSUM(L,KC) + ICERATE(L)
          ICEVOL(L)  = 0.0
        ENDIF
      ENDDO
    ENDDO
    !$OMP END DO
    !$OMP SINGLE
    DELTICE1 = 0.
    !$OMP END SINGLE
  ENDIF

  ! *** EVAPORATION AND RAINFALL.  BOTH RAINT AND EVAPT ARE IN M/S
  IF( NASER > 0  .OR. EVAP.IFLAG > 0 .OR. RAIN.IFLAG > 0 )THEN
    IF( IEVAP > 1 .OR. ISTOPT(2) == 1 .OR. ISTOPT(2) == 5 )THEN
      ! *** COMPUTE PARTIAL PRESSURE OF WATER VAPOR AT WATER TEMPERATURE (mb)
      !$OMP DO PRIVATE(ND,LF,LL,L)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          SVPW(L) = (10.**((0.7859+0.03477*TEM(L,KC))/(1.+0.00412*TEM(L,KC))))
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF

    ! *** PRECIPITATION, TREAT AS ICE, IF BELOW FREEZING AND HAS EXISTING ICE COVER
    IF( ISTRAN(2) > 0 )THEN
      IF( ISICE > 2 )THEN
        TMPVAL = 999.8426/RHOI*DELT
        !$OMP DO PRIVATE(ND,LF,LL,L)
        DO ND=1,NDM
          LF=2+(ND-1)*LDM
          LL=MIN(LF+LDM-1,LA)

          IF( LFRAZIL .AND. IEVAP == 1 )THEN
            DO L=LF,LL
              IF( ICETHICK(L) > 0. )THEN
                ICETHICK(L) = ICETHICK(L) - EVAPT(L)*TMPVAL          ! *** SUBLIMATION
                IF( ICETHICK(L) < 0.  )THEN
                  ICETHICK(L) = 0.0
                ENDIF
                IF( ICETHICK(L) < MINICETHICK ) ICECELL(L) = .FALSE.
                EVAPT(L) = 0.0
              ENDIF
            ENDDO
          ENDIF

          DO L=LF,LL
            IF( TATMT(L) < 0. )THEN
              IF( ICECELL(L) .OR. .NOT. LMASKDRY(L) )THEN
                ICETHICK(L) = ICETHICK(L) + RAINT(L)*TMPVAL
                IF( ICETHICK(L) < MINICETHICK )THEN
                  ICECELL(L) = .FALSE.
                ELSE
                  ICECELL(L) = .TRUE.
                ENDIF
                RAINT(L) = 0.0
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        !$OMP END DO
      ENDIF
    ENDIF

    IF( IEVAP == 0 )THEN
      ! *** IGNORE EVAPORATION
      !$OMP DO PRIVATE(ND,LF,LL,L,K)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*RAINT(L)
        ENDDO
      ENDDO
      !$OMP END DO
    
    ELSEIF( IEVAP == 1 .OR. EVAP.IFLAG > 0 .OR. RAIN.IFLAG > 0 )THEN
      ! *** USE INPUT EVAPORATION
      !$OMP DO PRIVATE(ND,LF,LL,L,K)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        ENDDO
      ENDDO
      !$OMP END DO

    ELSEIF( IEVAP == 2 )THEN
      ! *** COMPUTE EVAPORATION (ORIGINAL EFDC APPROACH)
      !$OMP DO PRIVATE(ND,LF,LL,L)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          EVAPT(L) = CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW(L) - VPAT(L))/PATMT(L)
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        ENDDO
      ENDDO
      !$OMP END DO

    ELSEIF( IEVAP == 11 )THEN
      ! *** COMPUTE EVAPORATION: RYAN-HARLEMAN (FROM CE-QUAL-W2, COLE & WELLS, 2011)
      !$OMP DO PRIVATE(ND,LF,LL,L,TVW,TVA,DTV,DTVL,LAMBDA,TM,VPTG)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          TVA  = (TATMT(L) +273.0)/(1.0-0.378*VPAT(L)/PATMT(L))
          TVW  = (TEM(L,KC)+273.0)/(1.0-0.378*SVPW(L)/PATMT(L))
          DTV  = TVW-TVA
          DTVL =  0.0084*WINDST(L)**3
          IF( DTV < DTVL) DTV = DTVL
          LAMBDA = 3.59*DTV**0.3333333
          LAMBDA = LAMBDA+4.26*WINDST(L)

          TM       = (TEM(L,KC)+TDEWT(L))*0.5
          VPTG     =  0.35 + 0.015*TM + 0.0012*TM*TM
          EVAPT(L) = VPTG*(TEM(L,KC) - TDEWT(L))*LAMBDA/2.45E9
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        ENDDO
      ENDDO
      !$OMP END DO

    ELSEIF( IEVAP > 2 .AND. IEVAP < 11 )THEN
      ! *** COMPUTE EVAPORATION USING THE WIND FUNCTION APPROACH  f(W) = A+B*W+C*W**2
      ! *** THE f(W) UNITS HAVE BEEN CONVERTED TO M/S/millibar
      !$OMP DO PRIVATE(ND,LF,LL,L,CLEVAPTMP)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          CLEVAPTMP  = WINDFA + WINDFB*WINDST(L) + WINDFC*WINDST(L)**2
          EVAPT(L)   = CLEVAPTMP*(SVPW(L)-VPAT(L))
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        ENDDO
      ENDDO
      !$OMP END DO

    ELSEIF( IEVAP == 12 ) THEN
      ! *** Compute evaporation using the approach from Arifin et al., 2016
      ! *** This approach combines the adjusment for wind speed from Croley 1989 with the algorithsm in Quinn 1979
      !$OMP DO PRIVATE(ND,LF,LL,L,DENAIR,CLEVAPTMP,FWUP,QSW,QAA)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          DENAIR = PATMT(L)*100./287.05 /(TATMT(L) + 273.15)                          ! *** Density of air (kg/m^3)
          CLEVAPTMP = 8.*10E-4*DENAIR/RHOW(L,KC)
          FWUP = 1.607 + 0.92*1.2475*WINDST(L) - 0.28*(TATMT(L) - TEM(L,KC))          ! *** Adjusted wind speed factor
          QSW = 0.622*SVPW(L) /(PATMT(L) - 0.378*SVPW(L))                             ! *** The specific humidity of water
          QAA = 0.622*VPAT(L) /(PATMT(L) - 0.378*VPAT(L))                             ! *** The specific humidity of air
          
          EVAPT(L)   = CLEVAPTMP*FWUP*(QSW - QAA)
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        ENDDO
      ENDDO
      !$OMP END DO
      
    ENDIF

    ! *** REMOVE EVAP & RAIN FROM CELLS WITH ICE COVER
    IF( ISICE > 0 .AND. IEVAP > 1 )THEN
      !$OMP DO PRIVATE(ND,LF,LL,L,TVW,TVA,DTV,DTVL,LAMBDA,TM,VPTG)
      DO ND=1,NDM
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          IF( ICECELL(L) )THEN
            ! *** BACK OUT ANY ADDED FLOWS
            QSUM(L,KC) = QSUM(L,KC) - DXYP(L)*( RAINT(L)-EVAPT(L) )
            EVAPT(L)   = 0.
            RAINT(L)   = 0.
          ENDIF
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF

  ENDIF

  ! **  DETERMINE NET EXTERNAL VOLUME SOURCE/SINK
  !$OMP DO PRIVATE(ND,LF,LL,L,K)
  DO ND=1,NDM
    LF=2+(ND-1)*LDM
    LL=MIN(LF+LDM-1,LA)
    DO K=1,KC
      DO L=LF,LL
        QSUME(L) = QSUME(L) + QSUM(L,K)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** ADJUST VOLUME SOURCE AND SINKS
  IF( ISGWIE == 0 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,L,K,QEAVAIL,RAVAIL,TMPVAL,TF,I)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)
      DO L=LF,LL
        IF( QSUME(L) < 0. )THEN
          ! *** LIMIT WITHDRAWAL VOLUMES TO AVAILABLE WATER IN CELL
          IF( HP(L) < HWET )THEN
            IF( HP(L) <= HDRY9 )THEN
              ! *** TURN OFF ANY OUTFLOWS
              QGW(L)   = 0.
              EVAPT(L) = 0.
              RAINT(L) = 0.

              ! *** ADJUST LAYER BY LAYER FLOWS
              DO K=1,KC
                QSUM(L,K) = 0.0
              ENDDO
              QSUME(L)=0.

              ! *** REDUCE ICE THICKNESS
              IF( ISICE > 2 .AND. ( ISTL_ == 3 .OR. IS2TL > 0 ) )THEN
                IF( ICEVOL(L) < 0. )THEN
                  TMPVAL = ICEVOL(L)*999.8426/RHOI/DXYP(L)
                  ICETHICK(L) = ICETHICK(L) + TMPVAL
                  ICETHICK(L) = MAX(ICETHICK(L),0.)
                  IF( ICETHICK(L) < MINICETHICK )ICECELL(L) = .FALSE.
                  ICEVOL(L) = 0.
                ENDIF
              ENDIF
            ELSE
              ! *** CHECK AVAILBILITY OF WATER IN CELL
              QEAVAIL = DXYP(L)*(HP(L)-HDRY9)*DELTI
              IF( QEAVAIL < ABS(QSUME(L)) )THEN
                ! *** LIMIT WITHDRAWAL VOLUMES TO AVAILABLE WATER IN CELL
                QEAVAIL  = -QEAVAIL
                RAVAIL   = QEAVAIL/QSUME(L)

                QSUME(L) = QEAVAIL
                ! *** ADJUST LAYER BY LAYER FLOWS
                DO K=1,KC
                  QSUM(L,K) = QSUM(L,K)*RAVAIL
                ENDDO

                QGW(L)   = QGW(L)*RAVAIL
                EVAPT(L) = EVAPT(L)*RAVAIL
                RAINT(L) = RAINT(L)*RAVAIL

                ! *** REDUCE ICE THICKNESS
                IF( ISICE > 2 .AND. ( ISTL_ == 3 .OR. IS2TL > 0 ) )THEN
                  IF( ICEVOL(L) < 0. )THEN
                    ! *** FREEZING TEMPERATURE OF WATER
                    IF( ISTRAN(1) > 0 )THEN
                      IF( SAL(L,KC) < 35. )THEN
                        TF = -0.0545*SAL(L,KC)
                      ELSE
                        TF=-0.31462-0.04177*SAL(L,KC)-0.000166*SAL(L,KC)*SAL(L,KC)
                      ENDIF
                    ELSE
                      TF=0.0
                    ENDIF

                    ICEVOL(L) = ICEVOL(L)*RAVAIL
                    TMPVAL = ICEVOL(L)*999.8426/RHOI/DXYP(L)
                    IF( ISICE == 3 )THEN
                      ! *** REDUCE ICE COVER THICKNESS
                      ICETHICK(L) = ICETHICK(L) + TMPVAL
                      ICETHICK(L) = MAX(ICETHICK(L),0.)
                    ELSE
                      ! *** RESTORE THE FRAZIL ICE
                      IF( TMPVAL+ICETHICK(L) > 0. )THEN
                        ICETHICK(L) = ICETHICK(L) + TMPVAL
                      ELSE
                        FRAZILICE(L,KC) = FRAZILICE(L,KC) + TMPVAL
                        FRAZILICE(L,KC) = MAX(FRAZILICE(L,KC),0.)
                      ENDIF
                    ENDIF
                    IF( ICETHICK(L) < MINICETHICK )THEN
                      ICECELL(L) = .FALSE.
                    ELSE
                      ICECELL(L) = .TRUE.
                    ENDIF

                    ! *** RESTORE THE TEMPERATURE (HEAT)
                    TEM(L,KC) = TEM(L,KC)*RAVAIL + (TEM1(L,KC)-TF)*(1.-RAVAIL)
                  ENDIF
                ENDIF  ! *** End of ISICE > 2
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !$OMP END DO

  ELSE    ! *** IF( ISGWIE >= 1 )THEN

    ! *** ADJUST SOURCES AND SINKS ESTIMATING SURFACE AND GROUNDWATER
    ! *** AVAILABLE FOR EVAPOTRANSPIRATON AND INFILTRATION
    !$OMP DO PRIVATE(ND,LF,LL,L,K,DTAGW,RIFTRL,RAVAIL,QSUMIET,QEAVAIL,QSUMTMP,DIFQVOL)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)
      DO L=LF,LL
        EVAPSW(L) = 0.
        EVAPGW(L) = 0.

        IF( HP(L) > HDRY9 )THEN
          ! *** APPLY MAXIMUM ET
          EVAPSW(L) = EVAPT(L)*DXYP(L)      ! *** m3/s
          QGW(L) = 0.

          ! *** CALCULATE DEPTH OF ACTIVE GROUNDWATER ELEV BELOW SURFACE
          DTAGW=BELV(L)-AGWELV(L)
          IF( DTAGW > 0.0 )THEN

            ! *** INFLITRATION CAN OCCUR, CALCULATE LIMITING RATE TO BRING
            ! *** GW ELEV TO SOIL SURFACE
            RIFTRL=RNPOR*DTAGW*DELTI

            ! *** SET RIFTRL TO MIN OF LIMITING RATE OR ACTUAL RATE
            RIFTRL=MIN(RIFTRM,RIFTRL)

            ! *** ESTIMATE RATE BASED ON AVAILABLE SURFACE WATER
            RAVAIL=(H1P(L)-HDRY)*DELTI-EVAPT(L)

            ! *** SET RIFTRL TO MIN OF AVAILABLE RATE OR LIMITING RATE
            RIFTRL=MIN(RAVAIL,RIFTRL)

            ! *** CONVERT TO VOLUME FLOW UNITS
            QGW(L) = -RIFTRL*DXYP(L)           ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT)
          ENDIF

          ! *** ADJUST VOLUME OUTFLOWS OF WET CELLS
          IF( QSUME(L) < 0.0 .AND. HP(L) < HWET )THEN
            QSUMIET = EVAPSW(L)-QGW(L)
            QEAVAIL = DXYP(L)*(H1P(L)-HDRY)*DELTI-QSUMIET
            QEAVAIL = MAX(QEAVAIL,0.0)
            QEAVAIL = -QEAVAIL
            QSUMTMP = MAX(QSUME(L),QEAVAIL)

            ! *** UPDATE LAYER FLOWS
            DIFQVOL=QSUME(L)-QSUMTMP
            DO K=KSZ(L),KC
              QSUM(L,K) = QSUM(L,K) - DIFQVOL*DZC(L,K)
            ENDDO
            QSUME(L)=QSUMTMP
          ENDIF
        ELSE
          ! *** CELL IS 'DRY'
          QGW(L)  = 0.
          EVAPSW(L) = 0.
          QSUME(L) = MAX(QSUME(L),0.0)      ! *** ONLY ALLOW INFLOWS
          DO K=1,KC
            QSUM(L,K) = MAX(QSUM(L,K),0.0)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !$OMP END DO
  ENDIF
  !$OMP END PARALLEL

  ! **  WRITE DIAGNOSTIC FILE FOR VOLUME SOURCES,SINKS, ETC
  ITMPD=0
  IF( ISDIQ == 2 .AND. ISTL_ == 2 ) ITMPD=1
  IF( ISDIQ == 1 ) ITMPD=1
  NTT = 3 + NDYM + NTOX + NSED + NSND

6665 FORMAT(' US ELEVATION CONTROL STRUCTURE TABLE OUT OF BOUNDS ')
6666 FORMAT(' SINGLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS ')
6667 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
6668 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
6676 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, UP ')
6677 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
6678 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
6679 FORMAT(' HUF,HUL,HDF,HDL = ',4(2X,E12.4))
6686 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, DW ')
6687 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
6688 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
6111 FORMAT(' INVALID NQSIJ LOCATION, NQSIJ,I,J = ',3I5)
  RETURN

  END SUBROUTINE

