! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALFQC(ISTL_, IS2TL_, MVAR, MO, CON, CON1, IT)

  ! **  SUBROUTINE CALFQC CALCULATES MASS SOURCES AND SINKS ASSOCIATED  
  ! **  WITH CONSTANT AND TIME SERIES INFLOWS AND OUTFLOWS; CONTROL  
  ! **  STRUCTURE INFLOWS AND OUTFLOWS; WITHDRAWAL AND RETURN STRUCTURE  
  ! **  OUTFLOWS; AND  EMBEDED CHANNEL INFLOWS AND OUTFLOWS  
  !  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------!
  !    2015-11       PAUL M. CRAIG     ELIMINATED DUPLICATE 3TL/2TL FQC COMPUTATIONS
  !                                    ADDED THE NEW HYDRAULIC CONTROL STRUCTURES WITH
  !                                    RETURN FLOWS OPTIONALLY PUT IN DIFFRENT LAYERS
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2011-03       Paul M. Craig     Rewritten to F90
  !----------------------------------------------------------------------------------!

  USE GLOBAL 
  Use Variables_WQ
  
  use mpi
  use Variables_MPI
  
  IMPLICIT NONE
  
  ! *** Passed in variables
  INTEGER, INTENT(IN) :: ISTL_, IS2TL_, MVAR, MO, IT
  REAL, INTENT(INOUT)    :: CON(LCM,KCM), CON1(LCM,KCM)
  
  ! *** Local variables
  REAL    :: FQOUT(NQCTLM), QTOUT(NQCTLM)

  INTEGER :: M, L, LP, K, ID, JD, KD, NWR, IU, JU, KU, LU, NS, NC, NT
  INTEGER :: LD, NMD, NJP, LJP, KTMP, NTMP
  INTEGER :: NQSTMP,  NCSTMP,  NCTL
  INTEGER :: LMDCHHT,  LMDCHUT,  LMDCHVT
  
  REAL    :: QWRABS,RPORTS,QVJPTMP,QCJPTMP,QVJPENT,CONUP
  REAL    :: RQWD, QUKTMP, QVKTMP, TMPVAL, DOSAT

  LOGICAL :: USENQSIJ
  
  ! *** VARIABLES USED FOR ANTI-DIFFUSION CALCULATIONS
  ! ***
  ! *** FQCPAD  - MASS LOADINGS INTO MODEL DOMAIN FROM ALL BC TYPES FOR EACH CELL AND LAYER  (OUTFLOWS EXCLUDED)
  ! *** QSUMPAD - TOTAL FLOWS INTO MODEL DOMAIN FROM ALL BC TYPES FOR EACH CELL AND LAYER  (OUTFLOWS EXCLUDED)
  
  M = MO

  ! *** SELECTIVE ZEROING
  IF( KC > 1 )THEN
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=1,LC  
        FQC(L,KSZ(L),IT) = 0.  
      ENDDO  
    ENDIF

    ! *** ZERO EVAP/RAINFALL
    IF( MVAR == 2 )THEN        
      DO L=1,LC  
        FQC(L,KC,IT) = 0.  
      ENDDO  
      IF( ISADAC(MVAR) >= 2 )THEN
        DO L=1,LC  
          FQCPAD(L,KC,IT) = 0.  
        ENDDO  
      ENDIF
      IF( ISADAC(MVAR) > 0 )THEN
        DO L=1,LC  
          QSUMPAD(L,KC,IT) = 0.  
        ENDDO  
      ENDIF
    ENDIF
    
    ! *** ATMOSPHERIC TOXIC DEPOSITION
    IF( ISTRAN(5) > 0 .AND. MVAR == 5 )THEN
      NT = MO - (3 + NDYE)
      IF( (TOXDEP(NT).TXDRY+TOXDEP(NT).TXWET) /= 0. .OR. (TOXDEP(NT).ITXDRY+TOXDEP(NT).ITXWET) /= 0 )THEN   ! *** Zero if atmospheric deposition being used
        DO L=1,LC  
          FQC(L,KC,IT)  = 0.  
          CONQ(L,KC,IT) = 0.5*(3.*CON(L,KC)-CON1(L,KC))  
        ENDDO  
      ENDIF
    ENDIF
    
    ! *** ZERO ALL DEFINED BCs
    DO NS=1,NBCS
      L=LBCS(NS)
      DO K=1,KC
        FQC(L,K,IT) = 0.  
        CONQ(L,K,IT) = 0.
        FQCPAD(L,K,IT) = 0  
        QSUMPAD(L,K,IT) = 0.  
      ENDDO
    ENDDO

  ELSE
    ! *** KC=1 RESET ALL
    FQC(1:LCM,KC,IT) = 0.
    CONQ(1:LCM,KC,IT)=0.
    IF( ISADAC(MVAR) >= 2 )THEN
      FQCPAD(1:LCM,1:KCM,IT) = 0.
    ENDIF
    QSUMPAD(1:LCM,1:KCM,IT) = 0.
  ENDIF

  ! **  CORRECT CONCENTRATION FOR WATER SURFACE ELEVATION DATA ASSIMILATION
  !IF(ISWSEDA > 0 )THEN
  !  DO K=1,KC
  !    DO L=1,LC
  !      FQC(L,K)=FQC(L,K)+QWSEDA(L,K)*CONQ(L,K)
  !      QSUMPAD(L,K)=QSUMPAD(L,K)+MAX(QWSEDA(L,K),0.)
  !     QSUMNAD(L,K)=QSUMNAD(L,K)+MIN(QWSEDA(L,K),0.)
  !    ENDDO
  !  ENDDO
  !ENDIF

  ! *** HANDLE LEGACY WATER QUALITY MASS LOADING BC TYPE, I.E. IWQPSL = 0 OR IWQPSL = 1
  USENQSIJ = .TRUE.
  IF( MVAR == 8 .AND. IWQPSL /= 2 ) USENQSIJ = .FALSE.

  ! **************************************************************************************
  ! ***     INITIALIZE VOLUMETRIC SOURCE-SINK FLUXES AND AUXILLARY VARIABLES  
  ! **************************************************************************************

  ! *** CONQ INITIALIZATION: 3TL STANDARD TIME STEP
  IF( ISTL_ == 3 )THEN  
    ! *** INITIALIZE BOTTOM LAYER FOR GW INTERACTIONS
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=1,LC  
        CONQ(L,KSZ(L),IT) = CON(L,KSZ(L))  
      ENDDO  
    ENDIF

    ! *** INITIALIZE TOP LAYER FOR EVAP/RAINFALL
    IF( MVAR == 2 )THEN        
      DO L=1,LC  
        CONQ(L,KC,IT) = CON(L,KC)
      ENDDO  
    ENDIF
    
    ! *** INITIALIZE ALL DEFINED BCs
    DO NS=1,NBCS
      L=LBCS(NS)
      DO K=KSZ(L),KC
        CONQ(L,K,IT) = CON(L,K)  
      ENDDO  
    ENDDO  
  ENDIF  ! *** END OF CONQ INITIALIZATION: 3TL STANDARD TIME STEP

  ! *** CONQ INITIALIZATION: 3TL CORRECTION STEP
  IF( ISTL_ == 2 .AND. IS2TL_ == 0 )THEN  
    ! *** INITIALIZE BOTTOM LAYER FOR GW INTERACTIONS
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=1,LC  
        CONQ(L,KSZ(L),IT) = 0.5*(CON(L,KSZ(L))+CON1(L,KSZ(L)))  
      ENDDO  
    ENDIF

    ! *** INITIALIZE TOP LAYER FOR EVAP/RAINFALL
    IF( MVAR == 2 )THEN        
      DO L=1,LC  
        CONQ(L,KC,IT) = 0.5*(CON(L,KC)+CON1(L,KC))  
      ENDDO  
    ENDIF
        
    ! *** INITIALIZE ALL DEFINED BCs
    DO NS=1,NBCS
      L=LBCS(NS)
      DO K=KSZ(L),KC
        CONQ(L,K,IT) = 0.5*(CON(L,K)+CON1(L,K))  
      ENDDO  
    ENDDO  
  ENDIF   ! *** END OF CONQ INITIALIZATION: 3TL CORRECTION STEP
       
  ! *** CONQ INITIALIZATION: 2TL STANDARD TIME STEP
  IF( IS2TL_ == 1 )THEN  
    ! *** INITIALIZE BOTTOM LAYER FOR GW INTERACTIONS
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=1,LC  
        CONQ(L,KSZ(L),IT) = 0.5*(3.*CON(L,KSZ(L))-CON1(L,KSZ(L)))  
      ENDDO  
    ENDIF

    ! *** INITIALIZE TOP LAYER FOR EVAP/RAINFALL
    IF( MVAR == 2 )THEN        
      DO L=1,LC  
        CONQ(L,KC,IT) = 0.5*(3.*CON(L,KC)-CON1(L,KC))  
      ENDDO  
    ENDIF
    
    ! *** INITIALIZE ALL DEFINED BCs
    DO NS=1,NBCS
      L=LBCS(NS)
      DO K=KSZ(L),KC
        CONQ(L,K,IT) = 0.5*(3.*CON(L,K)-CON1(L,K))  
      ENDDO  
    ENDDO  
  ENDIF  ! *** END OF CONQ INITIALIZATION: 2TL STANDARD TIME STEP

  ! *** GET MASS OUT FOR CONTROL STRUCTURES
  IF( NQCTL > 0 )THEN
    DO NCTL=1,NQCTL  
      IU = HYD_STR(NCTL).IQCTLU  
      JU = HYD_STR(NCTL).JQCTLU  
      LU = LIJ(IU,JU)  
      FQOUT(NCTL) = 0.
      QTOUT(NCTL) = 0.
      DO K=KSZ(LU),KC  
        QTOUT(NCTL) = QTOUT(NCTL) + QCTLT(K,NCTL,1)
        FQOUT(NCTL) = FQOUT(NCTL) + QCTLT(K,NCTL,1)*CONQ(LU,K,IT)
      ENDDO  
    ENDDO  
  ENDIF
  
  ! *********************************************************************!
  ! *** STANDARD VOLUMETRICS SOURCE SINK LOCATIONS 
  ! *** EFDCPlus - 2015-11 - ELIMINATED DUPLICATE MASS LOADINGS SECTIONS FOR THE DIFFERENT TIME STEP OPTIONS 

  ! *** FLOW BOUNDARY CELLS
  IF( USENQSIJ )THEN
    DO NS=1,NQSIJ  
      L=LQS(NS)  
      NQSTMP=NQSERQ(NS)  
      NCSTMP=NCSERQ(NS,MVAR)  
      DO K=KSZ(L),KC  
        ! ***                             |----------------  IN  -----------------|   |-----------  OUT  ----------------|
        FQC(L,K,IT)     = FQC(L,K,IT)     + MAX(QSS(K,NS),     0.) *CQS(K,NS,M)       + MIN(QSS(K,NS),0.)     *CONQ(L,K,IT)  &
                                          + MAX(QSERCELL(K,NS),0.) *CSERT(K,NCSTMP,M) + MIN(QSERCELL(K,NS),0.)*CONQ(L,K,IT)  
        
        FQCPAD(L,K,IT)  = FQCPAD(L,K,IT)  + MAX(QSS(K,NS),0.)*CQS(K,NS,M) + MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,M)  
        QSUMPAD(L,K,IT) = QSUMPAD(L,K,IT) + MAX(QSS(K,NS),0.)             + MAX(QSERCELL(K,NS),0.)  
      ENDDO  
    ENDDO  
  ENDIF
    
  ! ***  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS
  IF( NQJPIJ > 0 )THEN  
    DO NJP=1,NQJPIJ  
      IF( ICALJP(NJP) == 1 )THEN  
        ! *** DISCHARGE ONLY USING QSER AND STANDARD CONCENTRATION SERIES
        RPORTS=FLOAT(NPORTJP(NJP))  
        LJP=LIJ(IQJP(NJP),JQJP(NJP))  
        KTMP=KEFFJP(NJP)  
        
        ! ***  QVJPTMP = Time series discharge from jet-plume  
        QVJPTMP = 0.  
        DO K=KSZ(LJP),KC  
          QVJPTMP = QVJPTMP + QSERT(K,NQSERJP(NJP))
        ENDDO  

        ! *** Remove mass due to plume
        QCJPTMP=0.  
        QVJPENT=0.  
        ! *** Remove entrainment flux and calculate total entraiment flux  
        DO K=KSZ(LJP),KC  
          FQC(LJP,K,IT)  = FQC(LJP,K,IT) - QJPENT(K,NJP)*CONQ(LJP,K,IT)*RPORTS    ! *** Remove mass due to entrainment into the plume (includes actual discharge flow)
          QCJPTMP        = QCJPTMP       + QJPENT(K,NJP)*CONQ(LJP,K,IT)           ! *** Sum the total mass removed for each port
          QVJPENT        = QVJPENT       + QJPENT(K,NJP)                          ! *** Sum the total flows
          !QSUMNAD(LJP,K,IT)=QSUMNAD(LJP,K,IT)-RPORTS*QJPENT(K,NJP)  
        ENDDO  

        ! *** Place jet flux and entrainment flux is effective layer  
        FQC(LJP,KTMP,IT)     =     FQC(LJP,KTMP,IT) + RPORTS*QCJPTMP + RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M) + RPORTS*QVJPTMP*CSERT(1,NCSERJP(NJP,MVAR),M)   ! *** Add entrained mass plus discharge mass into the effective layer 
        FQCPAD(LJP,KTMP,IT)  =  FQCPAD(LJP,KTMP,IT) + RPORTS*QCJPTMP + RPORTS*QQCJP(NJP)*CQCJP(1,NJP,M) + RPORTS*QVJPTMP*CSERT(1,NCSERJP(NJP,MVAR),M)  
        QSUMPAD(LJP,KTMP,IT) = QSUMPAD(LJP,KTMP,IT) + RPORTS*QVJPENT + RPORTS*QQCJP(NJP)+RPORTS*QVJPTMP  
      ENDIF  
      
      IF( ICALJP(NJP) == 2 )THEN  
        ! *** WITHDRAWAL/RETURN TYPE USING W/R SERIES
        RPORTS=FLOAT(NPORTJP(NJP))  
        LJP=LIJ(IQJP(NJP),JQJP(NJP))  
        KTMP=KEFFJP(NJP)  
        NS=NQWRSERJP(NJP)  
        LU=LIJ(IUPCJP(NJP),JUPCJP(NJP))  
        KU=KUPCJP(NJP)  
        CONUP=CONQ(LU,KU,IT)  
        QCJPTMP=0.  
        QVJPENT=0.  

        ! *** REMOVE ENTRAIMENT FLUX AND CALCULATE TOTAL ENTRAINMENT  
        DO K=KSZ(LJP),KC  
          FQC(LJP,K,IT)=FQC(LJP,K,IT)-RPORTS*QJPENT(K,NJP)*CONQ(LJP,K,IT)  
          QCJPTMP=QCJPTMP+QJPENT(K,NJP)*CONQ(LJP,K,IT)  
          QVJPENT=QVJPENT+QJPENT(K,NJP)  
          !QSUMNAD(LJP,K,IT)=QSUMNAD(LJP,K,IT)-RPORTS*QJPENT(K,NJP)  
        ENDDO  

        ! *** PLACE ENTRAINMENT, CONSTANT AND TIME SERIES FLUXES IN EFFECTIVE CELL  
        FQC(LJP,KTMP,IT)     =     FQC(LJP,KTMP,IT) + RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)+RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)  
        FQCPAD(LJP,KTMP,IT)  =  FQCPAD(LJP,KTMP,IT) + RPORTS*QCJPTMP+RPORTS*QWRCJP(NJP)*(CWRCJP(NJP,M)+CONUP)+RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M)+CONUP)  
        QSUMPAD(LJP,KTMP,IT) = QSUMPAD(LJP,KTMP,IT) + RPORTS*QVJPENT+RPORTS*QWRCJP(NJP)+RPORTS*QWRSERT(NS)  

        ! *** REMOVAL WITHDRAWAL FROM UPSTREAM CELL  
        FQC(LU,KU,IT) = FQC(LU,KU,IT) - RPORTS*QWRCJP(NJP)*CONUP-RPORTS*QWRSERT(NS)*CONUP  
        !QSUMNAD(LU,KU,IT)=QSUMNAD(LU,KU,IT)-RPORTS*QWRCJP(NJP)-RPORTS*QWRSERT(NS)  
      ENDIF  
    ENDDO  
  ENDIF  

  ! *** CONTROL STRUCTURES
  DO NCTL=1,NQCTL  
    RQWD=1.  
    IU = HYD_STR(NCTL).IQCTLU  
    JU = HYD_STR(NCTL).JQCTLU  
    LU = LIJ(IU,JU)  
    ID = HYD_STR(NCTL).IQCTLD  
    JD = HYD_STR(NCTL).JQCTLD  
    IF( ID == 0 .AND. JD == 0 )THEN  
      LD = LC  
      RQWD = 0.  
    ELSE  
      LD = LIJ(ID,JD)  
    ENDIF  
    DO K=KSZ(LU),KC  
      FQC(LU,K,IT) = FQC(LU,K,IT) - QCTLT(K,NCTL,1)*CONQ(LU,K,IT)
      IF( QTOUT(NCTL) > 1.E-6 )THEN
        ! *** INTO THE DOMAIN ONLY
        TMPVAL = RQWD*QCTLT(K,NCTL,2)/QTOUT(NCTL)*FQOUT(NCTL)
        FQC(LD,K,IT)     = FQC(LD,K,IT)    + TMPVAL
        FQCPAD(LD,K,IT)  = FQCPAD(LD,K,IT) + TMPVAL
        QSUMPAD(LD,K,IT) = QSUMPAD(LD,K,IT)+ RQWD*QCTLT(K,NCTL,2)  
        !QSUMNAD(L,K,IT)=QSUMNAD(L,K,IT)-QCTLT(K,NCTL,1)  
      ENDIF
    ENDDO  
  ENDDO  

  ! *** WITHDRAWAL CONCENTRATION RISE RETURN
  DO NWR=1,NQWR  
    ! *** Handle +/- Flows for Withdrawal/Return Structures
    NQSTMP=WITH_RET(NWR).NQWRSERQ  
    IF( QWRSERT(NQSTMP) >= 0. )THEN
      !! *** Original Withdrawal/Return
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
      WITH_RET(NWR).QWR=0.  ! *** Only allow time variable flows when using -W/R 
    ENDIF
    QWRABS = ABS(QWRSERT(NQSTMP))
    LU=LIJ(IU,JU)  
    LD=LIJ(ID,JD)  
    NCSTMP=WITH_RET(NWR).NQWRSERQ  

    FQC(LU,KU,IT)=   FQC(LU,KU,IT)   -(WITH_RET(NWR).QWR+QWRABS)*CONQ(LU,KU,IT)  
    FQC(LD,KD,IT)=   FQC(LD,KD,IT)   +WITH_RET(NWR).QWR*(CONQ(LU,KU,IT)+CQWR(NWR,M))+QWRABS*(CONQ(LU,KU,IT)+CQWRSERT(NCSTMP,M))  
    FQCPAD(LD,KD,IT)=FQCPAD(LD,KD,IT)+WITH_RET(NWR).QWR*(CONQ(LU,KU,IT)+CQWR(NWR,M))+QWRABS*(CONQ(LU,KU,IT)+CQWRSERT(NCSTMP,M))  
    QSUMPAD(LD,KD,IT)=QSUMPAD(LD,KD,IT)+WITH_RET(NWR).QWR+QWRABS  
    !QSUMNAD(LU,KU,IT)=QSUMNAD(LU,KU,IT)-(WITH_RET(NWR).QWR+QWRSERT(NQSTMP))  
  ENDDO  

  ! *** SUBGRID SCALE CHANNEL EXCHANGE
  IF( MDCHH >= 1 )THEN  
    DO K=1,KC  
      DO NMD=1,MDCHH  
        LMDCHHT=LMDCHH(NMD)  
        LMDCHUT=LMDCHU(NMD)  
        LMDCHVT=LMDCHV(NMD)  
        IF( MDCHTYP(NMD) == 1 )THEN  
          QUKTMP=QCHANU(NMD)*DZC(L,K)  
          QVKTMP=0.  
        ENDIF  
        IF( MDCHTYP(NMD) == 2 )THEN  
          QVKTMP=QCHANV(NMD)*DZC(L,K)  
          QUKTMP=0.  
        ENDIF  
        IF( MDCHTYP(NMD) == 3 )THEN  
          QUKTMP=QCHANU(NMD)*DZC(L,K)  
          QVKTMP=QCHANV(NMD)*DZC(L,K)  
        ENDIF  
        FQC(LMDCHHT,K,IT)=FQC(LMDCHHT,K,IT)+MAX(QUKTMP,0.)*CONQ(LMDCHUT,K,IT)+MIN(QUKTMP,0.)*CONQ(LMDCHHT,K,IT)  &
                                            +MAX(QVKTMP,0.)*CONQ(LMDCHVT,K,IT)+MIN(QVKTMP,0.)*CONQ(LMDCHHT,K,IT)  
        FQC(LMDCHUT,K,IT)=FQC(LMDCHUT,K,IT)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K,IT)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K,IT)  
        FQC(LMDCHVT,K,IT)=FQC(LMDCHVT,K,IT)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K,IT)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K,IT)  
      ENDDO  
    ENDDO  
  ENDIF  

  ! *** GROUNDWATER FLUXES
  IF( ISGWIT > 0. )THEN  
    IF( ISTRAN(5) > 0 .AND. MVAR == 5 )THEN
      NTMP = 3 + NDYM + NSED + NSND + NTOX
      DO NC=1,NTMP
        DO L=2,LA
          CONGW(L,NC) = GWCSERT(NGWSL(L),NC) 
        END DO
      END DO
    ENDIF

    ! *** Compute in/out fluxes but bypass sediment AND chemical loads.  Chemical loads handled in CALTOXB.
    IF( MVAR /= 5 .AND. MVAR /= 6 .AND. MVAR /= 7 )THEN
      DO L=2,LA
        IF( QGW(L) < 0. )THEN        ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
          FQC(L,KSZ(L),IT) = FQC(L,KSZ(L),IT) + QGW(L)*CONQ(L,KSZ(L),IT)  
        ELSE
          FQC(L,KSZ(L),IT) = FQC(L,KSZ(L),IT) + QGW(L)*GWCSERT(NGWSL(L),M)
        ENDIF
      ENDDO  
    ENDIF
  ENDIF  
  
  ! *** ADD DRY AND/OR WET ATMOSPHERIC DEPOSITION
  IF( ISTRAN(5) > 0 .AND. MVAR == 5 )THEN
    NT = MO - (3 + NDYE)
    
    ! *** DRY DEPOSITION
    IF( TOXDEP(NT).TXDRY /= 0. .AND. TOXDEP(NT).ITXDRY == 0 )THEN
      ! *** CONSTANT FLUX
      DO LP=1,LASED
        L = LSED(LP)
        !    MG/S    =      MG/S            MG/S/M^2     M^2
        FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXDRY*DXYP(L)
      ENDDO
    ELSEIF( TOXDEP(NT).ITXDRY == 1 )THEN
      ! *** TIME VARIABLE FLUX
      DO LP=1,LASED
        L = LSED(LP)
        !    MG/S    =      MG/S            MG/S/M^2       M^2
        FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXDRYCUR*DXYP(L)
      ENDDO
    ELSEIF( TOXDEP(NT).ITXDRY == 2 )THEN
      ! *** CONSTANT + TIME VARIABLE FLUX
      DO LP=1,LASED
        L = LSED(LP)
        !    MG/S    =      MG/S              MG/S/M^2           MG/S/M^2       M^2
        FQC(L,KC,IT) = FQC(L,KC,IT) + (TOXDEP(NT).TXDRYCUR + TOXDEP(NT).TXDRY)*DXYP(L)
      ENDDO
    ENDIF
    
    ! *** WET DEPOSITION
    IF( RAINTSUM0 > 0. )THEN
      IF( TOXDEP(NT).TXWET /= 0. .AND. TOXDEP(NT).ITXWET == 0 )THEN
        ! *** CONSTANT FLUX
        DO LP=1,LASED
          L = LSED(LP)
          ! ***                               MG/M3        M2      M/S
          FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXWET*DXYP(L)*RAINT(L)     ! *** MG/S
        ENDDO
      ELSEIF( TOXDEP(NT).ITXWET == 1 )THEN
        ! *** TIME VARIABLE FLUX
        DO LP=1,LASED
          L = LSED(LP)
          FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXWETCUR*DXYP(L)*RAINT(L) 
        ENDDO
      ELSEIF( TOXDEP(NT).ITXWET == 2 )THEN
        ! *** CONSTANT + TIME VARIABLE FLUX
        DO LP=1,LASED
          L = LSED(LP)
          FQC(L,KC,IT) = FQC(L,KC,IT) + (TOXDEP(NT).TXWETCUR + TOXDEP(NT).TXWET)*DXYP(L)*RAINT(L) 
        ENDDO
      ENDIF
    ENDIF
  ENDIF

  ! *** ICE MELT FOR WATER QUALITY CONSTITUENTS
  IF( MVAR == 8 .AND. ISICE > 2 )THEN  
    IF( MO == MSVDOX )THEN
      DOSAT = 14.   ! *** FRESH WATER D.O. SATURATION FOR NEAR FREEZING WATER (G/M3)
      DO L=2,LA
        IF( ICERATE(L) > 0. )THEN
          FQC(L,KC,IT)     = FQC(L,KC,IT)     + ICERATE(L)*DOSAT
          !FQCPAD(L,KC,IT)  = FQCPAD(L,KC,IT)  + ICERATE(L)*DOSAT
          !QSUMPAD(L,KC,IT) = QSUMPAD(L,KC,IT) + ICERATE(L)
        ENDIF
        ICERATE(L) = 0.
      ENDDO
    ENDIF
  ENDIF
  
  ! *** TEMPERATURE ADJUSTMENTS FOR RAINFALL & EVAPORATION
  IF( M == 2 .AND. RAINTSUM0 > 0. )THEN 
    IF( ISTOPT(2) == 0 .OR. ISTOPT(2) == 3 )THEN  
      ! *** CONSTANT TEMPERATURE
      DO L=2,LA  
        FQC(L,KC,IT) = FQC(L,KC,IT) + RAINT(L)*TEMO*DXYP(L)  
      ENDDO  
    ELSE
      ! *** SPATIALLY VARIABLE TEMPERATURES
      DO L=2,LA  
        FQC(L,KC,IT)     = FQC(L,KC,IT)     + RAINT(L)*TATMT(L)*DXYP(L)  
        FQCPAD(L,KC,IT)  = FQCPAD(L,KC,IT)  + RAINT(L)*TATMT(L)*DXYP(L)
        QSUMPAD(L,KC,IT) = QSUMPAD(L,KC,IT) + RAINT(L)*DXYP(L)
      ENDDO  
    ENDIF  
  ENDIF
  
  IF( M == 2 )THEN  
    IF( ISTOPT(2) == 0 )THEN  
      DO L=2,LA  
        FQC(L,KC,IT) = FQC(L,KC,IT) - EVAPSW(L)*CONQ(L,KC,IT)  
      ENDDO  
    ENDIF  
  ENDIF  

  RETURN  

END  

