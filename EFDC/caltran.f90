! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @DETAILS DESCRIBE THE SUBROUTINE OR WHATEVEVER
! @DATE 10/22020
! @AUTHOR PAUL 
  
  SUBROUTINE CALTRAN (ISTL_, IS2TL_, MVAR, MO, CON, CON1, IW, IT, WCCUTOFF, ISKIP)  

  ! **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE  
  ! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSTITUENT M LEADING TO  
  ! **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES  
  ! **  THE NUMBER OF TIME LEVELS IN THE STEP  

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2015-01       PAUL M. CRAIG     Added fully coupled Ice Sub-model with Frazil Ice Transport
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH

  USE GLOBAL 
  Use Variables_MPI
  !  USE OMP_LIB
#ifdef _MPI
  Use MPI
#endif

  IMPLICIT NONE  
    
  ! *** Passed in variables
  INTEGER, INTENT(IN)    :: ISTL_, IS2TL_, MVAR, MO, IW, IT
  INTEGER, INTENT(INOUT) :: ISKIP
  REAL, INTENT(IN)       :: WCCUTOFF
  REAL, INTENT(INOUT)    :: CON(LCM,KCM), CON1(LCM,KCM)  
  
  ! *** Local variables
  REAL    :: DDELT, DDELTA, DDELTD2, WCSUM, BCSUM, WCMAX
  REAL    :: CTMP, CBT, AUHU, AVHV, UTERM, VTERM, WTERM  
  REAL    :: CBSTMP, CBWTMP, CBETMP, CBNTMP, UHU, VHV, AWW, WW  
  REAL    :: CWMAX, CEMAX, CSMAX, CNMAX, CMAXT  
  REAL    :: CWMIN, CEMIN, CSMIN, CNMIN, CMINT  
  
  INTEGER :: M, iunit, II, JJ, I, J, NS, ITRANFLOC  
  INTEGER :: ISUD, K, KMAX, NSID, IOBC, NMNLOD, NCELL, NPAR
  INTEGER :: LP, L, LN, LS, LE, LW, LSE, LNW, LL, LMAX

  ! *** ZERO ANY INITIALIZED CONCENTRATIONS BELOW BOTTOM ACTIVE LAYER
  IF( IGRIDV > 0 .AND. N < 5 )THEN
    DO L=2,LA
      DO K=1,KSZ(L)-1
        CON(L,K)  = 0.0
        CON1(L,K) = 0.0
      ENDDO
    ENDDO
  ENDIF
  
  ! *** SET UP FLOC TRANSPORT
  ITRANFLOC=0
  
  ISUD=1  
  IF( ISDYNSTP == 0 )THEN 
    ! *** FIXED DELTA T
    DDELT=DT2  
    DDELTA=DT2  
    IF( ISCDCA(MVAR) == 2 ) DDELTA=DT   ! *** Central Differencing (3TL) [EXPERIMENTAL]
    DDELTD2=DT
    IF( ISCDCA(MVAR) == 1 ) ISUD=0 
    IF( ISTL_ /= 3 )THEN  
      DDELT=DT  
      DDELTA=DT  
      DDELTD2=0.5*DT  
      IF( IS2TIM == 0 ) ISUD=0          ! *** 3TL CORRECTOR TIME STEP (ISTL=2)
    ENDIF  
  ELSE  
    ! *** DYNAMIC DELTA T
    DDELT=DTDYN  
    DDELTA=DTDYN  
    DDELTD2=0.5*DTDYN  
  END IF  
 
  M = MO  
  IF( IS2TL_ == 1 )THEN  
    ! *** ADVANCE CONCENTRATIONS BEFORE THE 2TL TRANSPORT CALCULATIONS
    ! *** SKIP UPDATING VARIABLES IF ALREADY COMPLETED BEFORE THIS STEP
    IF( MVAR /= 8 )THEN
      CON1(:,:) = CON(:,:)
    ENDIF  
  ENDIF  

  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)  
    DO K=1,KC  
      WQBCCON(IOBC,K,IW)  = CON(L,K)  
      WQBCCON1(IOBC,K,IW) = CON1(L,K)  
    ENDDO  
  ENDDO  
  
  ! **  CALCULATED EXTERNAL SOURCES AND SINKS  
  CALL CALFQC(ISTL_, IS2TL_, MVAR, M, CON, CON1, IT)  
  
  ! *** SKIP TRANSPORT IF NOTHING TO TRANSPORT
  ISKIP = 0
  IF( WCCUTOFF > 0.0 )THEN
    
    ! *** SUM ALL DEFINED BCs
    BCSUM = 0.0
    DO NS=1,NBCS
      L=LBCS(NS)
      DO K=1,KC
        BCSUM = BCSUM + FQC(L,K,IT)
      ENDDO
    ENDDO
    
    IF( BCSUM <= 0.0 )THEN
      WCSUM = 0.0
      !WCMAX = -1E32
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          WCSUM = WCSUM + CON(L,K)  
          !IF( CON(L,K) > WCMAX )THEN
          !  WCMAX = CON(L,K)
          !  LMAX = L
          !  KMAX = K
          !ENDIF
        ENDDO  
      ENDDO
      IF( WCSUM < WCCUTOFF )THEN
        ISKIP = 1
        RETURN
      ENDIF
    ENDIF
  ENDIF
  
  ! **  SELECT TRANSPORT OPTION, ISPLIT=1 FOR HORIZONTAL-VERTICAL  
  ! **  OPERATOR SPLITTING  
  ! **  BEGIN COMBINED ADVECTION SCHEME  
  ! **  ADVECTIVE FLUX CALCULATION  
 
  IF( ISTL_ == 2 ) GOTO 300  
  IF( ISCDCA(MVAR) == 0 ) GOTO 300   ! *** Upwind Differencing  (3TL)
  IF( ISCDCA(MVAR) == 1 ) GOTO 400   ! *** Central Differencing (3TL)  
  IF( ISCDCA(MVAR) == 2 ) GOTO 350   ! *** Upwind Differencing  (3TL) [EXPERIMENTAL]
  ! *** UHDY2 AND VHDX2 ARE LAYER FLOWS (SIGMA-Z VERSION)

  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED  
  ! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY  

  300 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      FUHUD(L,K,IW) = UHDY2(L,K)*CON1(LUPU(L,K),K)  
      FVHUD(L,K,IW) = VHDX2(L,K)*CON1(LUPV(L,K),K)  
    ENDDO  
  ENDDO
  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        FWUU(L,K,IW) = W2(L,K)*CON1(L,KUPW(L,K))  
      ENDDO
    ENDDO
  ENDIF  
  GOTO 500  

  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! **  AVERAGED BETWEEN  (N-1) AND (N+1) AND ADVECTED FIELD AVERAGED  
  ! **  BETWEEN AT (N-1) AND (N) IF ISTL 3 ONLY  

  350 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      CONTD(L,K,IT) = 0.5*(CON(L,K) + CON1(L,K))+DDELT*0.5*FQC(L,K,IT)*DXYIP(L)/H2P(L)  
    ENDDO  
  ENDDO  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      FUHUD(L,K,IW) = UHDY2(L,K)*CONTD(LUPU(L,K),K,IT)  
      FVHUD(L,K,IW) = VHDX2(L,K)*CONTD(LUPV(L,K),K,IT)  
    ENDDO  
  ENDDO  
  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        FWUU(L,K,IW) = W2(L,K)*CONTD(L,KUPW(L,K),IT)  
      ENDDO  
    ENDDO  
  ENDIF  
  GOTO 500  

  ! **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT  
  ! **  AVERAGED BETWEEN (N+1) AND (N-1) AND TRANSPORTED FIELD AT (N)  

  400 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      FUHUD(L,K,IW) = 0.5*UHDY2(L,K)*(CON(L,K) + CON(LWC(L),K))  
      FVHUD(L,K,IW) = 0.5*VHDX2(L,K)*(CON(L,K) + CON(LSC(L),K))  
    ENDDO  
  ENDDO  
  DO K=1,KC  
    DO LL=1,NCBS  
      L=LCBS(LL)  
      LN=LNC(L)  
      IF( VHDX2(LN,K) < 0.) FVHUD(LN,K,IW) = VHDX2(LN,K)*CON1(LN,K)  
    ENDDO  
    DO LL=1,NCBW  
      L=LCBW(LL)  
      IF( UHDY2(LEC(L),K) < 0.) FUHUD(LEC(L),K,IW) = UHDY2(LEC(L),K)*CON1(LEC(L),K)  
    ENDDO  
    DO LL=1,NCBE  
      L=LCBE(LL)  
      IF( UHDY2(L,K) > 0.) FUHUD(L,K,IW) = UHDY2(L,K)*CON1(LWC(L),K)  
    ENDDO  
    DO LL=1,NCBN  
      L=LCBN(LL)  
      LS =LSC(L)  
      IF( VHDX2(L,K) > 0.) FVHUD(L,K,IW) = VHDX2(L,K)*CON1(LS,K)  
    ENDDO  
  ENDDO  
  DO K=1,KS  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0) 
      FWUU(L,K,IW) = 0.5*W2(L,K)*(CON(L,K+1) + CON(L,K))  
    ENDDO  
  ENDDO  

  ! **  STANDARD ADVECTION CALCULATION  
  500 CONTINUE  
 
  ! *** CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX (PMC MOVED)  
  IF( ISHDMF == 2 .OR. ISHDMF == 4 ) CALL CALDIFF (CON1,IW)
    
  ! *** BEGIN IF ON TRANSPORT OPTION CHOICE  
  ! *** IF ISACAC EQ 0 INCLUDE FQC MASS SOURCES IN UPDATE  
  IF( ISCDCA(MVAR) == 0 )THEN
    ! *** Upwind Differencing (3TL & 2TL)  
    IF( ISTL_ == 2 )THEN  
      DO K=1,KC
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          CD(L,K,IT) = CON1(L,K)*H1PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                       FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                    + (FWUU(L,K-1,IW)-FWUU(L,K,IW)) ) 
        ENDDO  
      ENDDO

      IF( ISFCT(MVAR) >= 1 .AND. ISADAC(MVAR) > 0 )THEN 
        CON2(:,:,IW) = MAX(CON1(:,:),0.0)
      ENDIF  
    
    ELSE  ! *** IF ISTL = 3
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          CD(L,K,IT) = CON1(L,K)*H2PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                       FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                    + (FWUU(L,K-1,IW)-FWUU(L,K,IW))  )  
        ENDDO  
      ENDDO
      
      IF( ISFCT(MVAR) >= 1 .AND. ISADAC(MVAR) > 0 )THEN
        CON2(:,:,IW) = MAX(CON(:,:),0.0)
      ENDIF  
    ENDIF     ! *** ENDIF ON TIME LEVEL CHOICE FOR ISCDCA=0  
  
    IF( IS2TL_ == 0  .AND.  ISUD == 1 )THEN  
      ! *** ADVANCE CON1 TO CON (3TL)
      CON1(:,:) = CON(:,:)
    ENDIF  

    ! *** UPDATE NEW CONCENTRATIONS  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CON(L,K) = CD(L,K,IT)*HPKI(L,K)
      ENDDO
    ENDDO

  ELSE  ! *** ELSE ON TRANSPORT OPTION CHOICE: ISCDCA(MVAR)/=0
   
    ! *** Central Differencing (3TL)
    ! *** Experimental Upwind Differencing (3TL)  
    IF( ISTL_ == 2 )THEN  
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          CD(L,K,IT) = CON1(L,K)*H1PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                    &
                                                     FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                  + (FWUU(L,K-1,IW)-FWUU(L,K,IW)) )  
        ENDDO  
      ENDDO
      
    ELSE   ! *** ISTL == 3  
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0) 
          CD(L,K,IT) = CON1(L,K)*H2PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                       FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                    + (FWUU(L,K-1,IW)-FWUU(L,K,IW)) )  
        ENDDO  
      ENDDO  
    ENDIF   ! *** ENDIF ON TIME LEVEL CHOICE FOR ISCDCA/=0  

    ! *** SAVE CONCENTRATION FOR FLUX CORRECTOR
    IF( ISFCT(MVAR) >= 1 )THEN  
      DO K=1,KC
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          CON2(L,K,IW) = MAX(CON1(L,K),0.0)
        ENDDO
      ENDDO
    ENDIF  

    IF( ISUD == 1 )THEN  
      ! *** ADVANCE CON1 TO CON (3TL SOLUTION)
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          CON1(L,K) = CON(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  
      
    ! *** COMPUTE THE CURRENT CONCENTRATIONS
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CON(L,K) = CD(L,K,IT)*HPKI(L,K)
      ENDDO  
    ENDDO  
    
  ENDIF ! *** END OF TRANSPORT OPTION CHOICE  

  ! *** RESTORE ORIGINAL CONCENTRATIONS PRIOR TO APPLYING OPEN BC'S - 2TL & 3TL
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)  
    DO K=1,KC  
      CON1(L,K) = WQBCCON1(IOBC,K,IW)
    ENDDO  
  ENDDO  

  ! ******************************************************************************************
  ! *** APPLY OPEN BOUNDARY CONDITIONS, BASED ON DIRECTION OF FLOW  
  ! *** SOUTH OPEN BC, WITHOUT FLOCS
  IF( NCBS > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBS  
        NSID=NCSERS(LL,MVAR)  
        L=LCBS(LL) 
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        LN=LNC(L)  
        IF( VHDX2(LN,K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
                                    CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K) - FVHUD(LNC(L),K,IW))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2 ) CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K) - FVHUD(LNC(L),K,IW))*DXYIP(L)*HPKI(L,K)  
            IF( ISCDCA(MVAR) == 2 ) CTMP = 0.5*(CON1(L,K) + CON(L,K)) + 0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)  &  
                                         + DDELT*(0.5*VHDX2(LN,K)*(CON1(L,K) + CON(L,K))-FVHUD(LNC(L),K,IW))*DXYIP(L)*HPKI(L,K)
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          IF( M == 1 .AND. BSC > 0.0 )THEN
            ! *** LIMIT CONCENTRATIONS TO MAXIMUM BC CONCENTRATIONS AT BOTTOM LAYER (SALINITY ONLY)
            CBSTMP = CBS(LL,1,M) + CSERT(1,NSID,M)  
            IF( CON(L,K) > CBSTMP )THEN
              CON(L,K) = CBSTMP  
            ENDIF
          ENDIF  
          CLOS(LL,K,M) = CON(L,K)
          NLOS(LL,K,M) = NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT = WTCI(K,1)*CBS(LL,1,M) + WTCI(K,2)*CBS(LL,2,M) + CSERT(K,NSID,M)  
          NMNLOD = NITER - NLOS(LL,K,M)  
          IF( NMNLOD >= NTSCRS(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBSTMP = CLOS(LL,K,M) + (CBT-CLOS(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRS(LL))  
            CON(L,K) = MAX(CBSTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** WEST OPEN BC, WITHOUT FLOCS  
  IF( NCBW > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBW  
        NSID=NCSERW(LL,MVAR)  
        L=LCBW(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        IF( UHDY2(LEC(L),K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
                                    CTMP = CON1(L,K) + DDELT*(UHDY2(LEC(L),K)*CON1(L,K) - FUHUD(LEC(L),K,IW))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2 ) CTMP = CON1(L,K) + DDELT*(UHDY2(LEC(L),K)*CON1(L,K) - FUHUD(LEC(L),K,IW))*DXYIP(L)*HPKI(L,K)
            IF( ISCDCA(MVAR) == 2 ) CTMP = 0.5*(CON1(L,K) + CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)         &  
                                          + DDELT*(0.5*UHDY2(LEC(L),K)*(CON1(L,K) + CON(L,K))-FUHUD(LEC(L),K,IW))*DXYIP(L)*HPKI(L,K) 
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBWTMP = CBW(LL,1,M) + CSERT(1,NSID,M)  
          IF( M == 1 .AND. BSC > 0.0 .AND. CON(L,K) > CBWTMP ) CON(L,K) = CBWTMP    ! *** PREVENT INVERTED DENSITY PROFILES AT OPEN BOUNDARY
          CLOW(LL,K,M) = CON(L,K)  
          NLOW(LL,K,M) = NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT = WTCI(K,1)*CBW(LL,1,M) + WTCI(K,2)*CBW(LL,2,M) + CSERT(K,NSID,M)  
          NMNLOD = NITER - NLOW(LL,K,M)  
          IF( NMNLOD >= NTSCRW(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBWTMP=CLOW(LL,K,M)+(CBT-CLOW(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRW(LL))  
            CON(L,K) = MAX(CBWTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** EAST OPEN BC, WITHOUT FLOCS  
  IF( NCBE > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBE  
        NSID = NCSERE(LL,MVAR)  
        L = LCBE(LL)
        NMNLOD = -1
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) ) CYCLE 
        IF( UHDY2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
                                    CTMP = CON1(L,K) + DDELT*(FUHUD(L,K,IW) - UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2 ) CTMP = CON1(L,K) + DDELT*(FUHUD(L,K,IW) - UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)  
            IF( ISCDCA(MVAR) == 2 ) CTMP = 0.5*(CON1(L,K) + CON(L,K)) + 0.5*(CON1(L,K) - CON(L,K))*H2P(L)*HPKI(L,K)+DDELT*(FUHUD(L,K,IW) &  
                                          -0.5*UHDY2(L,K)*(CON1(L,K) + CON(L,K)))*DXYIP(L)*HPKI(L,K)
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBETMP = CBE(LL,1,M) + CSERT(1,NSID,M)  
          IF( M == 1 .AND. BSC > 0.0 .AND. CON(L,K) > CBETMP ) CON(L,K) = CBETMP    ! *** PREVENT INVERTED DENSITY PROFILES AT OPEN BOUNDARY
          CLOE(LL,K,M) = CON(L,K)  
          NLOE(LL,K,M) = NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT = WTCI(K,1)*CBE(LL,1,M) + WTCI(K,2)*CBE(LL,2,M) + CSERT(K,NSID,M)  
          NMNLOD = NITER - NLOE(LL,K,M)  
          IF( NMNLOD >= NTSCRE(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBETMP = CLOE(LL,K,M) + (CBT-CLOE(LL,K,M)) * FLOAT(NMNLOD)/FLOAT(NTSCRE(LL))  
            CON(L,K) = MAX(CBETMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF
    
  ! *** NORTH OPEN BC, WITHOUT FLOCS  
  IF( NCBN > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBN  
        NSID=NCSERN(LL,MVAR)  
        L=LCBN(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        LS=LSC(L)  
        IF( VHDX2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
                                    CTMP = CON1(L,K) + DDELT*(FVHUD(L,K,IW) - VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2 ) CTMP = CON1(L,K) + DDELT*(FVHUD(L,K,IW) - VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K) 
            IF( ISCDCA(MVAR) == 2 ) CTMP = 0.5*(CON1(L,K) + CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2PK(L,K)*HPKI(L,K) + DDELT*(FVHUD(L,K,IW) &  
                                          -0.5*VHDX2(L,K)*(CON1(L,K) + CON(L,K)))*DXYIP(L)*HPKI(L,K) 
          ENDIF  
          CON(L,K) = MAX(CTMP,0.)
          CBNTMP = CBN(LL,1,M) + CSERT(1,NSID,M)  
          IF( M == 1 .AND. BSC > 0.0 .AND. CON(L,K) > CBNTMP ) CON(L,K) = CBNTMP    ! *** PREVENT INVERTED DENSITY PROFILES AT OPEN BOUNDARY
          CLON(LL,K,M) = CON(L,K)  
          NLON(LL,K,M) = NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT = WTCI(K,1)*CBN(LL,1,M) + WTCI(K,2)*CBN(LL,2,M) + CSERT(K,NSID,M)  
          NMNLOD = NITER - NLON(LL,K,M)  
          IF( NMNLOD >= NTSCRN(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBNTMP=CLON(LL,K,M)+(CBT-CLON(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRN(LL))  
            CON(L,K) = MAX(CBNTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    
  
  ! ****************************************************************************************
 
  ! *** ZERO HEAT FLUXES
  IF( MVAR == 2 )THEN        
    ! *** ZERO EVAP/RAINFALL
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
  
  !  $ print *, ' CALTRAN-Leaving: Thread, MVAR, MO',IT,MVAR,MO  
  
  RETURN  
END  
