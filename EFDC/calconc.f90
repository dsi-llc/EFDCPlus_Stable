! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALCONC

  ! *** SUBROUTINE CALCULATES THE CONCENTRATION OF DISSOLVED AND
  ! *** SUSPENDED CONSTITUENTS, INCLUDING SALINITY, TEMPERATURE, DYE,
  ! *** AND SUSPENDED SEDIMENT AT TIME LEVEL (N+1). THE VALUE OF ISTL
  ! *** INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2015-01       PAUL M. CRAIG     Added fully coupled Ice Sub-model with Frazil Ice Transport
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2014-08       D H CHUNG         SET EXPLICIT PRECISIONS OF INTEGER & REAL
  !    2011-03       PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP
  !    2002-05       John Hamrick      Modified calls to calbal and budget subroutines
  !                                     added calls to bal2t2, bal2t3
  !------------------------------------------------------------------------------------------------!

  USE GLOBAL
  Use Budget
  Use Allocate_Initialize      
  Use Variables_Propwash
  USE OMP_LIB
  USE HEAT_MODULE, ONLY:CALHEAT
  
# ifdef _MPI  
  Use MPI
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  Use Communicate_Ghost_Routines
# endif

  IMPLICIT NONE

  INTEGER :: K, L, LP, NS, ND, IT, IW, I, J
  INTEGER :: NTMP, LF, LL, LE, LN, LG, IERR
  INTEGER, SAVE :: NICE, IFIRST, NANTIDIFF
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ISKIP

  REAL      :: CDTMP, RCDZKMK, RCDZKK
  REAL,SAVE :: SEDTIME

  REAL(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS, TTDS1, TTDS2, TMPCOMM, TWAIT   ! *** Model timing temporary variables
  REAL(RKD)           :: CCUBTMP, CCMBTMP                     ! *** Vertical diffusion temporary variables

  ! *** VERTICAL DIFFUSION VARIABLES
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: CCLBTMP
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: EEB
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: VCU

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TOXASM
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SEDASM
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SNDASM

  IF( .NOT. ALLOCATED(EEB) )THEN
    Call AllocateDSI(CCLBTMP,   LCM,  KCM,  0.0)
    Call AllocateDSI(EEB,       LCM,  KCM,  0.0)
    Call AllocateDSI(VCU,       LCM,  KCM,  0.0)
    
    Call AllocateDSI(TOXASM,  NTXM,     0.0)
    Call AllocateDSI(SEDASM,  NSCM2,    0.0)
    Call AllocateDSI(SNDASM,  NSNM2,    0.0)
    Call AllocateDSI(ISKIP, NACTIVEWC,    0)

    NICE = 0
    SEDTIME = 0.0
    IFIRST = 0
    
    NANTIDIFF = 0
    DO IW=1,NACTIVEWC
      IF( .NOT. (ISADAC(IACTIVEWC1(IW)) == 0 .OR. ISCDCA(IACTIVEWC1(IW)) == 1) )THEN 
        NANTIDIFF = NANTIDIFF + 1
      ENDIF
    ENDDO
  ENDIF

  DELT=DT2
  IF( ISTL == 2 )THEN
    IF( ISDYNSTP == 0 )THEN
      DELT=DT
    ELSE
      DELT=DTDYN
    END IF
  ENDIF

  ! *** MASS BALANCE
  IF( IS2TIM >= 1 )THEN
    IF( ISBAL >= 1 )THEN
      CALL BAL2T3A
    ENDIF
  ENDIF
  IT=1

  ! ***************************************************************************************
  ! *** 3D ADVECTI0N TRANSPORT CALCULATION, STANDARD TRANSPORT
  TTDS1 = DSTIME(0)
  TWAIT = 0.
  TMPITMP = 0.

  ! *** PRESPECIFY THE UPWIND CELLS FOR 3D ADVECTION
  !$OMP PARALLEL DO PRIVATE(ND,K,LP,L,LE,LN)
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)
        IF( UHDY2(L,K) >= 0.0 )THEN
          LUPU(L,K) = LWC(L)
        ELSE
          LUPU(L,K) = L
        END IF
        IF( VHDX2(L,K) >= 0.0 )THEN
          LUPV(L,K) = LSC(L)
        ELSE
          LUPV(L,K) = L
        END IF
      ENDDO
    ENDDO

    IF( KC > 1 )THEN
      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)
          IF( W2(L,K) >= 0. )THEN
            KUPW(L,K) = K
          ELSE
            KUPW(L,K) = K + 1
          END IF
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO
  
  ! *** EE7.2 - MUST ZERO ALL THREADS FOR EVERY INSTANCE DO TO POTENTIAL OF DIFFERENT THREADS MAY NOT BE ZEROED FOR CURRENT TIME
  ! *** ZERO DRY CELL FLUXES
  IF( LADRY > 0 )THEN
    DO LP=1,LADRY
      L=LDRY(LP)
      FUHUD(L,:,:)=0.
      FVHUD(L,:,:)=0.
      FWUU(L,:,:) =0.

      LN=LNC(L)
      LE=LEC(L)
      FUHUD(LEC(L),:,:)=0.
      FVHUD(LNC(L),:,:)=0.
    ENDDO
    
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LADRY
          L=LDRY(LP)
          FUHVD(L,K,ND)=0.
          FVHVD(L,K,ND)=0.
          UUUU(L,K,ND) =0.
          VVVV(L,K,ND) =0.
          DUU(L,K,ND)  =0.
          DVV(L,K,ND)  =0.
          POS(L,K,ND)  =0.
          WWWW(L,K,ND) =0.

          LN=LNC(L)
          LE=LEC(L)
          FUHVD(LE,K,ND)=0.
          FVHVD(LN,K,ND)=0.
          UUUU(LE,K,ND) =0.
          VVVV(LN,K,ND) =0.
          DUU(LE,K,ND)  =0.
          DVV(LN,K,ND)  =0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! *** 3D ADVECTI0N-DIFFUSION TRANSPORT FOR ALL WATER COLUMN CONSITUENTS (SEE VARINIT FOR WCV INITIALIZATION)
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(IW,IT,L) SCHEDULE(STATIC,1)
  DO IW=1,NACTIVEWC
    !$  IT = OMP_GET_THREAD_NUM() + 1

    CALL CALTRAN(IACTIVEWC1(IW), IACTIVEWC2(IW), WCV(IW).VAL0, WCV(IW).VAL1, IW, IT, WCV(IW).WCLIMIT, ISKIP(IW))

    IF( ISICE > 2 .AND. IACTIVEWC1(IW) == 8 .AND. IACTIVEWC2(IW) == MSVDOX .AND. ISKIP(IW) == 0 )THEN
      ! *** ZERO SURFACE MELT FLUX
      DO L=1,LC
        FQC(L,KC,IT) = 0.
      ENDDO
      
    ENDIF
  ENDDO
  !$OMP END DO
  
  ! *** APPLY ANTI-DIFFUSION AND FLUX CORRECTOR
  !$OMP DO PRIVATE(IW,IT) SCHEDULE(STATIC,1)
  DO IW=1,NACTIVEWC
    IF( ISKIP(IW) == 0 )THEN
      !$  IT = OMP_GET_THREAD_NUM() + 1
      CALL CALTRAN_AD(IACTIVEWC1(IW), IACTIVEWC2(IW), WCV(IW).VAL0, WCV(IW).VAL1, IW, IT)
    ENDIF
  ENDDO
  !$OMP END DO

# ifdef _MPI
  ! ****************************************************************************
  ! *** MPI communication for FUHUD, FVHUD & FWUU
  !$OMP SINGLE
  Call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)

  Call Communicate_CON2
  
  TMPITMP = DSTIME(0) - TTDS
  DSITIMING(6) = DSITIMING(6) + TMPITMP
  !$OMP END SINGLE
  ! ****************************************************************************
# endif

  !$OMP DO PRIVATE(IW,K,LP,L,CDTMP) SCHEDULE(STATIC,1)
  DO IW=1,NACTIVEWC
    IF( ISADAC(IACTIVEWC1(IW)) == 1 .AND. ISCDCA(IACTIVEWC1(IW)) == 0 )THEN 
      ! *** APPLY THE ANTI-DIFFUSIVE ADVECTION CALCULATION TO STANDARD DONOR CELL SCHEME
      IF( ISKIP(IW) == 0 )THEN
        DO K=1,KC  
          DO LP=1,LLWET(K,0)
            L=LKWET(LP,K,0)
            CDTMP             = WCV(IW).VAL0(L,K)*HPK(L,K) + DELT*( ( FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                           +        ( FWUU(L,K-1,IW)-FWUU(L,K,IW) )  )  
            WCV(IW).VAL0(L,K) = CDTMP*HPKI(L,K)
          ENDDO  
        ENDDO    
        
        ! *** RESET OPEN BC CONCENTRATIONS
        DO K = 1,KC
          DO LP=1,NBCSOP
            L = LOBCS(LP)
            WCV(IW).VAL0(L,K) = WQBCCON(LP,K,IW)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  
  IF( ISICE == 4 .AND. (IS2TL > 0 .OR. (IS2TL == 0 .AND. NCTBC /= NTSTBC)) )THEN
    IF( LFRAZIL )THEN
      ! *** FRAZIL ICE TRANSPORT
      CALL CALTRANICE(FRAZILICE,FRAZILICE1,1)
      NICE = NICE+1

    ELSEIF( NICE > 0 )THEN
      ! *** ADVANCE THE FRAZIL ICE VARIABLE
      FRAZILICE1 = FRAZILICE
      NICE = 0
    ENDIF
  ENDIF

  TTDS2 = DSTIME(0)
  Call MPI_barrier(MPI_Comm_World, ierr)
  TWAIT = TWAIT + (DSTIME(0)- TTDS2)

  TSADV = TSADV + (DSTIME(0)-TTDS1) - TMPITMP - TWAIT

  ! ******************************************************************************************
  ! *** VERTICAL DIFFUSION IMPLICIT HALF STEP CALCULATION
  IF( KC == 1 ) GOTO 1500

  TTDS1 = DSTIME(0)

  IF( LADRY > 0 )THEN
    DO LP=1,LADRY
      L=LDRY(LP)
      CCLBTMP(L,:) = 0.
      EEB(L,:) = 1.
      VCU(L,:) = 0.
    ENDDO
  ENDIF

  !$OMP PARALLEL DO PRIVATE(ND,LF,LL,LP,L,K,RCDZKK,CCUBTMP,CCMBTMP,RCDZKMK,IW) SCHEDULE(STATIC,1)
  DO ND=1,NDM
    LF = (ND-1)*LDMWET+1
    LL = MIN(LF+LDMWET-1,LAWET)

    ! -------------------------------------------------------------------------------
    ! *** COMPUTE THE DIFFUSIVE FLUXES

    ! *** BOTTOM LAYER
    DO LP=1,LLWET(KS,ND)
      L = LKWET(LP,KS,ND)
      RCDZKK  = -DELT*CDZKK(L,KSZ(L))
      CCUBTMP = RCDZKK*HPI(L)*AB(L,KSZ(L))
      CCMBTMP = 1._8-CCUBTMP
      EEB(L,KSZ(L)) = 1._8/CCMBTMP
      VCU(L,KSZ(L)) = CCUBTMP*EEB(L,KSZ(L))
    ENDDO

    ! *** MIDDLE LAYERS
    DO K=2,KS
      DO LP=1,LLWET(K-1,ND)
        L = LKWET(LP,K-1,ND)
        RCDZKMK      = -DELT*CDZKMK(L,K)
        RCDZKK       = -DELT*CDZKK(L,K)
        CCLBTMP(L,K) = RCDZKMK*HPI(L)*AB(L,K-1)
        CCUBTMP      = RCDZKK*HPI(L)*AB(L,K)
        CCMBTMP      = 1._8-CCLBTMP(L,K)-CCUBTMP
        EEB(L,K)     = 1._8/( CCMBTMP - CCLBTMP(L,K)*VCU(L,K-1) )
        VCU(L,K)     = CCUBTMP*EEB(L,K)
      ENDDO
    ENDDO

    ! *** TOP LAYER
    K=KC
    DO LP=1,LLWET(KS,ND)
      L = LKWET(LP,KS,ND)
      RCDZKMK      = -DELT*CDZKMK(L,K)
      CCLBTMP(L,K) = RCDZKMK*HPI(L)*AB(L,K-1)
      CCMBTMP      = 1._8-CCLBTMP(L,K)
      EEB(L,K)     = 1._8/( CCMBTMP - CCLBTMP(L,K)*VCU(L,K-1) )
    ENDDO

    ! -------------------------------------------------------------------------------
    ! *** APPLY THE DIFFUSION

    ! *** BOTTOM LAYER
    DO IW=1,NACTIVEWC
      IF( ISKIP(IW) == 1 ) CYCLE
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND)
        K = KSZ(L)
        WCV(IW).VAL0(L,K) = WCV(IW).VAL0(L,K)*EEB(L,K)
      ENDDO

      ! *** MIDDLE AND TOP LAYERS
      DO K=2,KC
        DO LP=1,LLWET(K-1,ND)
          L=LKWET(LP,K-1,ND)
          WCV(IW).VAL0(L,K) = ( WCV(IW).VAL0(L,K) - CCLBTMP(L,K)*WCV(IW).VAL0(L,K-1) )*EEB(L,K)
        ENDDO
      ENDDO

      ! *** FINAL PASS
      DO K=KS,1,-1
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          WCV(IW).VAL0(L,K) = WCV(IW).VAL0(L,K) - VCU(L,K)*WCV(IW).VAL0(L,K+1)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  TVDIF = TVDIF + (DSTIME(0)-TTDS1)
  ! *** END OF VERTICAL DIFFUSION STEP
1500 CONTINUE

  ! ***  COMPUTE TOTALS FOR SEDIMENT TRANSPORT (AFTER ADVECTION/DIFFUSION)
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    !$OMP PARALLEL DEFAULT(SHARED)
    ! *** ZERO SED SEDIMENT ACCUMULATION ARRAYS
    IF( ISTRAN(6) > 0 )THEN
      !$OMP DO PRIVATE(ND,L,K,LF,LL,LP,NS) SCHEDULE(STATIC,1)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1
        LL=MIN(LF+LDMWET-1,LAWET)

        ! *** WATER COLUMN
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)
            SEDT(L,K) = 0.
          ENDDO
        ENDDO

        DO NS=1,NSED2
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)
              SEDT(L,K) = SEDT(L,K) + SED(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF

    ! *** ZERO SND SEDIMENT ACCUMULATION ARRAYS
    IF( ISTRAN(7) > 0 )THEN
      !$OMP DO PRIVATE(ND,L,K,LF,LL,LP,NS) SCHEDULE(STATIC,1)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1
        LL=MIN(LF+LDMWET-1,LAWET)

        ! *** WATER COLUMN
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)
            SNDT(L,K) = 0.
          ENDDO
        ENDDO

        DO NS=1,NSND
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)
              SNDT(L,K) = SNDT(L,K) + SND(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF
    !$OMP END PARALLEL
  ENDIF

  ! ***************************************************************************************
  ! *** SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION
  IF( ISTRAN(2) >= 1 )THEN
    TTDS1 = DSTIME(0)
    CALL CALHEAT

    THEAT = THEAT + (DSTIME(0)-TTDS1)
  ENDIF

  ! *** APPLY DYE PROCESSES
  IF( ISTRAN(3) >= 1 ) CALL CALDYE

  ! **************************************************************************************************
  ! *** BOTTOM AND INTERNAL SEDIMENT AND TOXIC CONTAMINATION SOURCE-SINK CALCULATION

  IF( ISPROPWASH == 2 )THEN
    ! *** If using propeller momentum, compute propeller velocities at every timestep, but skip if computing below
    IF( .NOT. ((SEDTIME+DELT) >= SEDSTEP .AND. (ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) .AND. TIMEDAY >= SEDSTART) )THEN
      Call Propwash_Calc_Sequence(1)
    ENDIF
  ENDIF

  ! *** SEDIMENT AND TOXICS SETTLING,DEPOSITION,RESUSPENSION,ETC
  IF( ( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) .AND. TIMEDAY >= SEDSTART )THEN
    IF( IS2TIM >= 1 )THEN
      ! *** FOR TWO TIME LEVEL SIMULATION
      SEDTIME = SEDTIME + DELT
      IF( SEDTIME >= SEDSTEP )THEN
        DTSED = SEDTIME
        
        ! *** Call Propwash module if that option is selected
        IF( propwash_on .and. IFIRST > 0 )THEN
          Call Propwash_Calc_Sequence(0)
        END IF
        
        CALL SSEDTOX
        SEDTIME = 0.0
      ENDIF
    ELSE
      ! *** FOR THREE TIME LEVEL SIMULATION
      IF( NCTBC == 1 )THEN
        DTSED = FLOAT(NTSTBC)*DT
        
        ! *** Call Propwash module if that option is selected
        IF( propwash_on .and. IFIRST > 0 )THEN
          Call Propwash_Calc_Sequence(0)
        END IF
        
        CALL SSEDTOX
        SEDTIME = 0.0
      ENDIF
    ENDIF
  ENDIF
  IFIRST = 1

  ! *** OPTIONAL MASS BALANCE CALCULATION
  IF( IS2TIM == 0 )THEN
    IF( ISTL/=2 .AND. ISBAL >= 1 )THEN
      CALL CALBAL2
      CALL CALBAL3
      NTMP=MOD(NCTBC,2)
      IF( NTMP == 0 )THEN
        CALL CBALEV2
        CALL CBALEV3
      ELSE
        CALL CBALOD2
        CALL CBALOD3
      ENDIF
    ENDIF
  ENDIF

  ! *** CALLS TO TWO-TIME LEVEL BALANCES
  !
  IF( IS2TIM >= 1 )THEN
    IF( ISBAL >= 1 )THEN
      CALL BAL2T2
      CALL BAL2T3B(1)
    ENDIF
  ENDIF

  ! *** SEDIMENT BUDGET CALCULATION    (DLK 10/15)
  IF( IS2TIM == 0 )THEN
    IF( ISTL/=2 .AND. ISSBAL >= 1 )THEN
      CALL BUDGET2
      CALL BUDGET3
    ENDIF
  ENDIF

  ! *** 2014-09 - REMOVED DATA ASSIMILATION CODE (PMC)

  RETURN

END

