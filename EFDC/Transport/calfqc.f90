! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALFQC(MVAR, MO, CON, CON1, IT)

  ! ***  SUBROUTINE CALFQC CALCULATES MASS SOURCES AND SINKS ASSOCIATED  
  ! ***  WITH CONSTANT AND TIME SERIES INFLOWS AND OUTFLOWS; CONTROL  
  ! ***  STRUCTURE INFLOWS AND OUTFLOWS; WITHDRAWAL AND RETURN STRUCTURE  
  ! ***  OUTFLOWS; AND  EMBEDED CHANNEL INFLOWS AND OUTFLOWS  
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

  use GLOBAL 
  use Variables_WQ
  
  use mpi
  use Variables_MPI
  
  implicit none
  
  ! *** Passed in variables
  integer, intent(IN) :: MVAR, MO, IT
  real, intent(INOUT) :: CON(LCM,KCM), CON1(LCM,KCM)
  
  ! *** Local variables
  real    :: FQOUT(NQCTLM), QTOUT(NQCTLM)

  integer :: M, L, LP, K, ID, JD, KD, NWR, IU, JU, KU, LU, NS, NC, NT
  integer :: LD, NMD, NJP, LJP, KTMP, NTMP
  integer :: NQSTMP,  NCSTMP,  NCTL
  integer :: LMDCHHT,  LMDCHUT,  LMDCHVT
  
  real    :: QWRABS,RPORTS,QVJPTMP,QCJPTMP,QVJPENT,CONUP
  real    :: RQWD, QUKTMP, QVKTMP, TMPVAL, DOSAT

  logical :: USENQSIJ
  
  ! *** VARIABLES USED FOR ANTI-DIFFUSION CALCULATIONS
  ! ***
  ! *** FQCPAD  - MASS LOADINGS INTO MODEL DOMAIN FROM ALL BC TYPES FOR EACH CELL AND LAYER  (OUTFLOWS EXCLUDED)
  ! *** QSUMPAD - TOTAL FLOWS INTO MODEL DOMAIN FROM ALL BC TYPES FOR EACH CELL AND LAYER  (OUTFLOWS EXCLUDED)
  
  M = MO

  ! *** SELECTIVE ZEROING
  if( KC > 1 )then
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 1,LC  
        FQC(L,KSZ(L),IT) = 0.  
      enddo  
    endif

    ! *** ZERO EVAP/RAINFALL
    if( MVAR == 2 )then        
      do L = 1,LC  
        FQC(L,KC,IT) = 0.  
      enddo  
      ! *** FQCPAD and QSUMPAD are not used in EFDC+.  PMC - Retaining the code for future use.
      !IF( ISADAC(MVAR) >= 2 )then
      !  do L = 1,LC  
      !    FQCPAD(L,KC,IT) = 0.  
      !  enddo  
      !ENDIF
      !IF( ISADAC(MVAR) > 0 )then
      !  do L = 1,LC  
      !    QSUMPAD(L,KC,IT) = 0.  
      !  enddo  
      !ENDIF
    endif
    
    ! *** ATMOSPHERIC TOXIC DEPOSITION
    if( ISTRAN(5) > 0 .and. MVAR == 5 )then
      NT = MO - (3 + NDYE)
      if( (TOXDEP(NT).TXDRY+TOXDEP(NT).TXWET) /= 0. .or. (TOXDEP(NT).ITXDRY+TOXDEP(NT).ITXWET) /= 0 )then   ! *** Zero if atmospheric deposition being used
        do L = 1,LC  
          FQC(L,KC,IT)  = 0.  
          CONQ(L,KC,IT) = 0.5*(3.*CON(L,KC)-CON1(L,KC))  
        enddo  
      endif
    endif
    
    ! *** ZERO ALL DEFINED BCs
    do NS = 1,NBCS
      L = LBCS(NS)
      do K = 1,KC
        FQC(L,K,IT) = 0.  
        CONQ(L,K,IT) = 0.
        ! FQCPAD(L,K,IT) = 0  
        ! QSUMPAD(L,K,IT) = 0.  
      enddo
    enddo
    if( HDRYMOVE > 0.0 )then
      do NS = 1,NQSIJ  
        L = LQSMOVED(NS)
        do K = 1,KC
          FQC(L,K,IT) = 0.  
          CONQ(L,K,IT) = 0.
        enddo
      enddo
    endif

  else
    ! *** KC = 1 RESET ALL
    FQC(1:LCM,KC,IT) = 0.
    CONQ(1:LCM,KC,IT) = 0.
    ! if( ISADAC(MVAR) >= 2 )then
    !   FQCPAD(1:LCM,1:KCM,IT) = 0.
    ! ENDIF
    ! QSUMPAD(1:LCM,1:KCM,IT) = 0.
  endif

  ! ***  CORRECT CONCENTRATION FOR WATER SURFACE ELEVATION DATA ASSIMILATION
  !IF(ISWSEDA > 0 )then
  !  do K = 1,KC
  !    do L = 1,LC
  !      FQC(L,K) = FQC(L,K)+QWSEDA(L,K)*CONQ(L,K)
  !      QSUMPAD(L,K) = QSUMPAD(L,K)+MAX(QWSEDA(L,K),0.)
  !     QSUMNAD(L,K) = QSUMNAD(L,K)+MIN(QWSEDA(L,K),0.)
  !    enddo
  !  enddo
  !ENDIF

  ! *** HANDLE LEGACY WATER QUALITY MASS LOADING BC TYPE, I.E. IWQPSL = 0 OR IWQPSL = 1
  USENQSIJ = .TRUE.
  if( MVAR == 8 .and. IWQPSL /= 2 ) USENQSIJ = .FALSE.

  ! **************************************************************************************
  ! ***     INITIALIZE VOLUMETRIC SOURCE-SINK FLUXES AND AUXILLARY VARIABLES  
  ! **************************************************************************************

  ! *** CONQ INITIALIZATION: 3TL STANDARD TIME STEP
  if( ISTL == 3 )then  
    ! *** INITIALIZE BOTTOM LAYER FOR GW INTERACTIONS
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 1,LC  
        CONQ(L,KSZ(L),IT) = CON(L,KSZ(L))  
      enddo  
    endif

    ! *** INITIALIZE TOP LAYER FOR EVAP/RAINFALL
    if( MVAR == 2 )then        
      do L = 1,LC  
        CONQ(L,KC,IT) = CON(L,KC)
      enddo  
    endif
    
    ! *** INITIALIZE ALL DEFINED BCs
    do NS = 1,NBCS
      L = LBCS(NS)
      do K = KSZ(L),KC
        CONQ(L,K,IT) = CON(L,K)  
      enddo  
    enddo  
  endif  ! *** END OF CONQ INITIALIZATION: 3TL STANDARD TIME STEP

  ! *** CONQ INITIALIZATION: 3TL CORRECTION STEP
  if( ISTL == 2 .and. IS2TL == 0 )then  
    ! *** INITIALIZE BOTTOM LAYER FOR GW INTERACTIONS
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 1,LC  
        CONQ(L,KSZ(L),IT) = 0.5*(CON(L,KSZ(L))+CON1(L,KSZ(L)))  
      enddo  
    endif

    ! *** INITIALIZE TOP LAYER FOR EVAP/RAINFALL
    if( MVAR == 2 )then        
      do L = 1,LC  
        CONQ(L,KC,IT) = 0.5*(CON(L,KC)+CON1(L,KC))  
      enddo  
    endif
        
    ! *** INITIALIZE ALL DEFINED BCs
    do NS = 1,NBCS
      L = LBCS(NS)
      do K = KSZ(L),KC
        CONQ(L,K,IT) = 0.5*(CON(L,K)+CON1(L,K))  
      enddo  
    enddo  
  endif   ! *** END OF CONQ INITIALIZATION: 3TL CORRECTION STEP
       
  ! *** CONQ INITIALIZATION: 2TL STANDARD TIME STEP
  if( IS2TL == 1 )then  
    ! *** INITIALIZE BOTTOM LAYER FOR GW INTERACTIONS
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 1,LC  
        CONQ(L,KSZ(L),IT) = 0.5*(3.*CON(L,KSZ(L))-CON1(L,KSZ(L)))  
      enddo  
    endif

    ! *** INITIALIZE TOP LAYER FOR EVAP/RAINFALL
    if( MVAR == 2 )then        
      do L = 1,LC  
        CONQ(L,KC,IT) = 0.5*(3.*CON(L,KC)-CON1(L,KC))  
      enddo  
    endif
    
    ! *** INITIALIZE ALL DEFINED BCs
    do NS = 1,NBCS
      L = LBCS(NS)
      do K = KSZ(L),KC
        CONQ(L,K,IT) = 0.5*(3.*CON(L,K)-CON1(L,K))  
      enddo  
    enddo  
  endif  ! *** END OF CONQ INITIALIZATION: 2TL STANDARD TIME STEP

  ! *** GET MASS OUT FOR CONTROL STRUCTURES
  if( NQCTL > 0 )then
    do NCTL = 1,NQCTL  
      IU = HYD_STR(NCTL).IQCTLU  
      JU = HYD_STR(NCTL).JQCTLU  
      LU = LIJ(IU,JU)  
      FQOUT(NCTL) = 0.
      QTOUT(NCTL) = 0.
      do K = KSZ(LU),KC  
        QTOUT(NCTL) = QTOUT(NCTL) + QCTLT(K,NCTL,1)
        FQOUT(NCTL) = FQOUT(NCTL) + QCTLT(K,NCTL,1)*CONQ(LU,K,IT)
      enddo  
    enddo  
  endif
  
  ! *********************************************************************!
  ! *** STANDARD VOLUMETRICS SOURCE SINK LOCATIONS 
  ! *** EFDCPlus - 2015-11 - ELIMINATED DUPLICATE MASS LOADINGS SECTIONS FOR THE DIFFERENT TIME STEP OPTIONS 

  ! *** FLOW BOUNDARY CELLS
  if( USENQSIJ )then
    do NS = 1,NQSIJ  
      L = BCPS(NS).L  
      NQSTMP = BCPS(NS).NQSERQ  
      NCSTMP = BCPS(NS).NCSERQ(MVAR)  
      do K = KSZ(L),KC  
        ! ***                             |----------------  IN  -----------------|   |-----------  OUT  ----------------|
        FQC(L,K,IT)     = FQC(L,K,IT)     + max(QSS(K,NS),     0.) *CQS(K,NS,M)       + min(QSS(K,NS),0.)     *CONQ(L,K,IT)  &
                                          + max(QSERCELL(K,NS),0.) *CSERT(K,NCSTMP,M) + min(QSERCELL(K,NS),0.)*CONQ(L,K,IT)  
        
        ! FQCPAD(L,K,IT)  = FQCPAD(L,K,IT)  + max(QSS(K,NS),0.)*CQS(K,NS,M) + max(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,M)  
        ! QSUMPAD(L,K,IT) = QSUMPAD(L,K,IT) + max(QSS(K,NS),0.)             + max(QSERCELL(K,NS),0.)  
      enddo  
    enddo  
  else
    ! *** Handle WQ mass loading withdrawals
    do NS = 1,NQSIJ  
      L = BCPS(NS).L  
      NQSTMP = BCPS(NS).NQSERQ  
      NCSTMP = BCPS(NS).NCSERQ(MVAR)  
      do K = KSZ(L),KC  
        ! ***                             |-----------  OUT  ----------------|
        FQC(L,K,IT)     = FQC(L,K,IT)     + min(QSS(K,NS),0.)     *CONQ(L,K,IT)  &
                                          + min(QSERCELL(K,NS),0.)*CONQ(L,K,IT)  
      enddo  
    enddo
  endif
    
  ! ***  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS
  if( NQJPIJ > 0 )then  
    do NJP = 1,NQJPIJ  
      if( JET_PLM(NJP).ICALJP == 1 )then  
        ! *** DISCHARGE ONLY USING QSER AND STANDARD CONCENTRATION SERIES
        RPORTS = FLOAT(JET_PLM(NJP).NPORTJP)  
        LJP = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)  
        KTMP = KEFFJP(NJP)  
        
        ! ***  QVJPTMP = Time series discharge from jet-plume  
        QVJPTMP = 0.  
        do K = KSZ(LJP),KC  
          QVJPTMP = QVJPTMP + QSERT(K,JET_PLM(NJP).NQSERJP)
        enddo  

        ! *** Remove entrainment mass and calculate total entraiment flux  
        QCJPTMP = 0.  
        do K = KSZ(LJP),KC  
          FQC(LJP,K,IT)  = FQC(LJP,K,IT) - QJPENT(K,NJP)*CONQ(LJP,K,IT)*RPORTS    ! *** Remove mass due to entrainment into the plume
          QCJPTMP        = QCJPTMP       + QJPENT(K,NJP)*CONQ(LJP,K,IT)           ! *** Sum the total mass removed for each port
        enddo  

        ! *** Place jet flux and entrainment flux is effective layer  
        ! ***                                   Entrainment                        Constant Flow                                    Time Varying Flow
        FQC(LJP,KTMP,IT) = FQC(LJP,KTMP,IT) + RPORTS*QCJPTMP + RPORTS*JET_PLM(NJP).QQCJP*JET_PLM(NJP).CQCJP(1,M) + RPORTS*QVJPTMP*CSERT(1,JET_PLM(NJP).NCSERJP(MVAR),M)
      endif  
      
      if( JET_PLM(NJP).ICALJP == 2 )then  
        ! *** JET-PLUME type using withdrawal/return series
        RPORTS  = FLOAT(JET_PLM(NJP).NPORTJP)  
        LJP     = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)  
        KTMP    = KEFFJP(NJP)  
        NS      = JET_PLM(NJP).NQWRSERJP  
        LU      = LIJ(JET_PLM(NJP).IUPCJP,JET_PLM(NJP).JUPCJP)  
        KU      = JET_PLM(NJP).KUPCJP  
        CONUP   = CONQ(LU,KU,IT)  

        ! *** Remove entraiment mass and calculate total entrainment  
        QCJPTMP = 0.  
        do K = KSZ(LJP),KC  
          FQC(LJP,K,IT) = FQC(LJP,K,IT) - QJPENT(K,NJP)*CONQ(LJP,K,IT)*RPORTS    ! *** Remove mass due to entrainment into the plume  
          QCJPTMP = QCJPTMP             + QJPENT(K,NJP)*CONQ(LJP,K,IT)           ! *** Sum the total mass removed for each port
        enddo  

        ! *** Place entrainment, constant and time series fluxes in effective cell
        ! ***                                   Entrainment                         Constant Flow                                        Time Varying Flow
        FQC(LJP,KTMP,IT) = FQC(LJP,KTMP,IT) + RPORTS*QCJPTMP + RPORTS*JET_PLM(NJP).QWRCJP*(JET_PLM(NJP).CWRCJP(M) + CONUP) + RPORTS*QWRSERT(NS)*(CQWRSERT(NS,M) + CONUP)  

        ! *** Removal withdrawal from upstream cell  
        ! ***                                     Constant Flow               Time Varying Flow
        FQC(LU,KU,IT) = FQC(LU,KU,IT) - RPORTS*JET_PLM(NJP).QWRCJP*CONUP - RPORTS*QWRSERT(NS)*CONUP  
      endif  
    enddo  
  endif  

  ! *** CONTROL STRUCTURES
  do NCTL = 1,NQCTL  
    RQWD = 1.  
    IU = HYD_STR(NCTL).IQCTLU  
    JU = HYD_STR(NCTL).JQCTLU  
    LU = LIJ(IU,JU)  
    ID = HYD_STR(NCTL).IQCTLD  
    JD = HYD_STR(NCTL).JQCTLD  
    if( ID == 0 .and. JD == 0 )then  
      LD = LC  
      RQWD = 0.  
    else  
      LD = LIJ(ID,JD)  
    endif  
    do K = KSZ(LU),KC  
      FQC(LU,K,IT) = FQC(LU,K,IT) - QCTLT(K,NCTL,1)*CONQ(LU,K,IT)
      if( QTOUT(NCTL) > 1.E-6 )then
        ! *** INTO THE DOMAIN ONLY
        TMPVAL = RQWD*QCTLT(K,NCTL,2)/QTOUT(NCTL)*FQOUT(NCTL)
        FQC(LD,K,IT)     = FQC(LD,K,IT)    + TMPVAL
      endif
    enddo  
  enddo  

  ! *** WITHDRAWAL CONCENTRATION RISE RETURN
  do NWR = 1,NQWR  
    ! *** Handle +/- Flows for Withdrawal/Return Structures
    NQSTMP = WITH_RET(NWR).NQWRSERQ  
    if( QWRSERT(NQSTMP) >= 0. )then
      !! *** Original Withdrawal/Return
      IU = WITH_RET(NWR).IQWRU  
      JU = WITH_RET(NWR).JQWRU  
      KU = WITH_RET(NWR).KQWRU  
      ID = WITH_RET(NWR).IQWRD  
      JD = WITH_RET(NWR).JQWRD  
      KD = WITH_RET(NWR).KQWRD  
    else
      ! *** Reverse Flow Withdrawal/Return
      ID = WITH_RET(NWR).IQWRU  
      JD = WITH_RET(NWR).JQWRU  
      KD = WITH_RET(NWR).KQWRU  
      IU = WITH_RET(NWR).IQWRD  
      JU = WITH_RET(NWR).JQWRD  
      KU = WITH_RET(NWR).KQWRD 
      WITH_RET(NWR).QWR = 0.  ! *** Only allow time variable flows when using -W/R 
    endif
    QWRABS = ABS(QWRSERT(NQSTMP))
    LU = LIJ(IU,JU)  
    LD = LIJ(ID,JD)  
    NCSTMP = WITH_RET(NWR).NQWRSERQ  

    FQC(LU,KU,IT)=   FQC(LU,KU,IT) - (WITH_RET(NWR).QWR+QWRABS)*CONQ(LU,KU,IT)  
    FQC(LD,KD,IT)=   FQC(LD,KD,IT) + WITH_RET(NWR).QWR*(CONQ(LU,KU,IT)+CQWR(NWR,M))+QWRABS*(CONQ(LU,KU,IT)+CQWRSERT(NCSTMP,M))  
  enddo  

  ! *** SUBGRID SCALE CHANNEL EXCHANGE
  if( MDCHH >= 1 )then  
    do K = 1,KC  
      do NMD = 1,MDCHH  
        LMDCHHT = LMDCHH(NMD)  
        LMDCHUT = LMDCHU(NMD)  
        LMDCHVT = LMDCHV(NMD)  
        if( MDCHTYP(NMD) == 1 )then  
          QUKTMP = QCHANU(NMD)*DZC(L,K)  
          QVKTMP = 0.  
        endif  
        if( MDCHTYP(NMD) == 2 )then  
          QVKTMP = QCHANV(NMD)*DZC(L,K)  
          QUKTMP = 0.  
        endif  
        if( MDCHTYP(NMD) == 3 )then  
          QUKTMP = QCHANU(NMD)*DZC(L,K)  
          QVKTMP = QCHANV(NMD)*DZC(L,K)  
        endif  
        FQC(LMDCHHT,K,IT) = FQC(LMDCHHT,K,IT)+MAX(QUKTMP,0.)*CONQ(LMDCHUT,K,IT)+MIN(QUKTMP,0.)*CONQ(LMDCHHT,K,IT)  &
                                            +MAX(QVKTMP,0.)*CONQ(LMDCHVT,K,IT)+MIN(QVKTMP,0.)*CONQ(LMDCHHT,K,IT)  
        FQC(LMDCHUT,K,IT) = FQC(LMDCHUT,K,IT)-MAX(QUKTMP,0.)*CONQ(LMDCHUT,K,IT)-MIN(QUKTMP,0.)*CONQ(LMDCHHT,K,IT)  
        FQC(LMDCHVT,K,IT) = FQC(LMDCHVT,K,IT)-MAX(QVKTMP,0.)*CONQ(LMDCHVT,K,IT)-MIN(QVKTMP,0.)*CONQ(LMDCHHT,K,IT)  
      enddo  
    enddo  
  endif  

  ! *** GROUNDWATER FLUXES
  if( ISGWIT > 0. )then  
    if( ISTRAN(5) > 0 .and. MVAR == 5 )then
      NTMP = 3 + NDYM + NTOX + NSED + NSND
      do NC = 1,NTMP
        do L = 2,LA
          CONGW(L,NC) = GWCSERT(NGWSL(L),NC) 
        enddo
      enddo
    endif

    ! *** Compute in/out fluxes but bypass sediment AND chemical loads.  Chemical loads handled in CALTOXB.
    if( MVAR /= 5 .and. MVAR /= 6 .and. MVAR /= 7 )then
      do L = 2,LA
        if( QGW(L) < 0. )then        ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
          FQC(L,KSZ(L),IT) = FQC(L,KSZ(L),IT) + QGW(L)*CONQ(L,KSZ(L),IT)  
        else
          FQC(L,KSZ(L),IT) = FQC(L,KSZ(L),IT) + QGW(L)*GWCSERT(NGWSL(L),M)
        endif
      enddo  
    endif
  endif  
  
  ! *** ADD DRY AND/OR WET ATMOSPHERIC DEPOSITION
  if( ISTRAN(5) > 0 .and. MVAR == 5 )then
    NT = MO - (3 + NDYE)
    
    ! *** DRY DEPOSITION
    if( TOXDEP(NT).TXDRY /= 0. .and. TOXDEP(NT).ITXDRY == 0 )then
      ! *** CONSTANT FLUX
      do LP = 1,LASED
        L = LSED(LP)
        !    MG/S    =      MG/S            MG/S/M^2     M^2
        FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXDRY*DXYP(L)
      enddo
    elseif( TOXDEP(NT).ITXDRY == 1 )then
      ! *** TIME VARIABLE FLUX
      do LP = 1,LASED
        L = LSED(LP)
        !    MG/S    =      MG/S            MG/S/M^2       M^2
        FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXDRYCUR*DXYP(L)
      enddo
    elseif( TOXDEP(NT).ITXDRY == 2 )then
      ! *** CONSTANT + TIME VARIABLE FLUX
      do LP = 1,LASED
        L = LSED(LP)
        !    MG/S    =      MG/S              MG/S/M^2           MG/S/M^2       M^2
        FQC(L,KC,IT) = FQC(L,KC,IT) + (TOXDEP(NT).TXDRYCUR + TOXDEP(NT).TXDRY)*DXYP(L)
      enddo
    endif
    
    ! *** WET DEPOSITION
    if( RAINTSUM0 > 0. )then
      if( TOXDEP(NT).TXWET /= 0. .and. TOXDEP(NT).ITXWET == 0 )then
        ! *** CONSTANT FLUX
        do LP = 1,LASED
          L = LSED(LP)
          ! ***                               MG/M3        M2      M/S
          FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXWET*DXYP(L)*RAINT(L)     ! *** MG/S
        enddo
      elseif( TOXDEP(NT).ITXWET == 1 )then
        ! *** TIME VARIABLE FLUX
        do LP = 1,LASED
          L = LSED(LP)
          FQC(L,KC,IT) = FQC(L,KC,IT) + TOXDEP(NT).TXWETCUR*DXYP(L)*RAINT(L) 
        enddo
      elseif( TOXDEP(NT).ITXWET == 2 )then
        ! *** CONSTANT + TIME VARIABLE FLUX
        do LP = 1,LASED
          L = LSED(LP)
          FQC(L,KC,IT) = FQC(L,KC,IT) + (TOXDEP(NT).TXWETCUR + TOXDEP(NT).TXWET)*DXYP(L)*RAINT(L) 
        enddo
      endif
    endif
  endif

  ! *** ICE MELT FOR WATER QUALITY CONSTITUENTS
  if( MVAR == 8 .and. ISICE > 2 )then  
    if( MO == MSVDOX )then
      DOSAT = 14.   ! *** FRESH WATER D.O. SATURATION FOR NEAR FREEZING WATER (G/M3)
      do L = 2,LA
        if( ICERATE(L) > 0. )then
          FQC(L,KC,IT)     = FQC(L,KC,IT)     + ICERATE(L)*DOSAT
          ! FQCPAD(L,KC,IT)  = FQCPAD(L,KC,IT)  + ICERATE(L)*DOSAT
          ! QSUMPAD(L,KC,IT) = QSUMPAD(L,KC,IT) + ICERATE(L)
        endif
        ICERATE(L) = 0.
      enddo
    endif
  endif
  
  ! *** TEMPERATURE ADJUSTMENTS FOR RAINFALL & EVAPORATION
  if( M == 2 .and. RAINTSUM0 > 0. )then 
    if( ISTOPT(2) == 0 )then  
      ! *** CONSTANT TEMPERATURE
      do L = 2,LA  
        FQC(L,KC,IT) = FQC(L,KC,IT) + RAINT(L)*TEMO*DXYP(L)  
      enddo  
    else
      ! *** SPATIALLY VARIABLE TEMPERATURES
      do L = 2,LA  
        FQC(L,KC,IT)     = FQC(L,KC,IT)     + RAINT(L)*TATMT(L)*DXYP(L)  
        ! FQCPAD(L,KC,IT)  = FQCPAD(L,KC,IT)  + RAINT(L)*TATMT(L)*DXYP(L)
        ! QSUMPAD(L,KC,IT) = QSUMPAD(L,KC,IT) + RAINT(L)*DXYP(L)
      enddo  
    endif  
  endif
  
  if( M == 2 )then  
    if( ISTOPT(2) == 0 )then  
      do L = 2,LA  
        FQC(L,KC,IT) = FQC(L,KC,IT) - EVAPSW(L)*CONQ(L,KC,IT)  
      enddo  
    endif  
  endif  

  return

END  

