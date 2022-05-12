! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! *** SUBROUTINE CALTRAN_AD CALCULATES THE ANTIDIFFION AND THE FLUX CORRECT, IF REQUESTED
  
  SUBROUTINE CALTRAN_AD (ISTL_, IS2TL_, MVAR, MO, CON, CON1, IW, IT)  

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------------!
  !    2020-01       PAUL M. CRAIG     Split CALTRAN to two components for MPI

  USE GLOBAL 
  Use Variables_MPI
#ifdef _MPI
  Use MPI
#endif

  IMPLICIT NONE  
    
  ! *** Passed in variables
  INTEGER, INTENT(IN) :: ISTL_, IS2TL_, MVAR, MO, IW, IT  
  REAL, INTENT(INOUT) :: CON(LCM,KCM), CON1(LCM,KCM)  
  
  ! *** Local variables
  REAL    :: BSMALL, DDELT, DDELTA 
  REAL    :: CTMP, CBT, AUHU, AVHV, UTERM, VTERM, WTERM  
  REAL    :: CBSTMP, CBWTMP, CBETMP, CBNTMP, UHU, VHV, AWW, WW  
  REAL    :: CWMAX, CEMAX, CSMAX, CNMAX, CMAXT  
  REAL    :: CWMIN, CEMIN, CSMIN, CNMIN, CMINT  
  
  INTEGER :: M, iunit, II, JJ, I, J
  INTEGER :: K, NSID, IOBC, NMNLOD
  INTEGER :: LP,L, LN, LS, LE, LW, LSE, LNW, LL, ITRANFLOC  

  ! ****************************************************************************************
  
  ! **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATIONS WITH FLUX CORRECTOR  
  IF( ISADAC(MVAR) == 0 ) RETURN    ! *** SKIP IF NOT USING ANTI-DIFFUSION
  IF( ISCDCA(MVAR) == 1 ) RETURN    ! *** SKIP IF CENTRAL DIFFERENCE (3TL)

  BSMALL=1.0E-6  
  IF( ISDYNSTP == 0 )THEN 
    ! *** FIXED DELTA T
    DDELT = DT2  
    DDELTA = DT2  
    IF( ISCDCA(MVAR) == 2 ) DDELTA = DT   ! *** Central Differencing (3TL) [EXPERIMENTAL]
    IF( ISTL_/=3 )THEN  
      DDELT=DT  
      DDELTA=DT  
    ENDIF  
  ELSE  
    ! *** DYNAMIC DELTA T
    DDELT=DTDYN  
    DDELTA=DTDYN  
  END IF  
 
  ! ----------------------------------------------------------------------------------------
  ! **  STANDARD ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION  
  
  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)
    DO K=1,KC  
      WQBCCON(IOBC,K,IW) = CON(L,K)  
    ENDDO  
  ENDDO  

  ! *** GET ONLY POSITIVE CONCENTRATIONS FOR ANTI-DIFUSION CALCULATIONS
  DO K=1,KC
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      POS(L,K,IT) = MAX(CON(L,K),0.) 
    ENDDO  
  ENDDO  
   
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0) 
      UUUU(L,K,IT) = U2(L,K)*( POS(L,K,IT) - POS(LWC(L),K,IT) )*DXIU(L)  
      VVVV(L,K,IT) = V2(L,K)*( POS(L,K,IT) - POS(LSC(L),K,IT) )*DYIV(L)  
    ENDDO  
  ENDDO

  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        WWWW(L,K,IT) = W2(L,K)*( POS(L,K+1,IT) - POS(L,K,IT) )*HPI(L)*DZIG(L,K)
      ENDDO  
    ENDDO  
  ENDIF

  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      LE  = LEC(L)
      LW  = LWC(L)
      LN  = LNC(L)  
      LS  = LSC(L)  
      LNW = LNWC(L)  
      LSE = LSEC(L)
      
      ! *** U COMPONENTS
      AUHU = ABS(UHDY2(L,K))  
      UTERM = AUHU*( POS(L,K,IT) - POS(LW,K,IT) )  
      IF( UHDY2(L,K) >= 0.0 )THEN  
        ! *** QUANTITIES FROM WEST CELL
        UTERM = UTERM - 0.5*DDELTA*UHDY2(L,K)*( VVVV(LNW,K,IT)+VVVV(LW,K,IT) + WWWW(LW,K,IT)+WWWW(LW,K-1,IT) + UUUU(L,K,IT)+UUUU(LW,K,IT) )
      ELSE  
        ! *** QUANTITIES FROM CURRENT CELL
        UTERM = UTERM - 0.5*DDELTA*UHDY2(L,K)*( VVVV(LN,K,IT) +VVVV(L,K,IT)  + WWWW(L,K,IT) +WWWW(L,K-1,IT)  + UUUU(L,K,IT)+UUUU(LE,K,IT) )  
      ENDIF  
      UHU = UTERM/( POS(L,K,IT) + POS(LW,K,IT) + BSMALL )  
      FUHUD(L,K,IW) = MAX(UHU,0.)*POS(LW,K,IT) + MIN(UHU,0.)*POS(L,K,IT)

      ! *** V COMPONENTS
      AVHV = ABS(VHDX2(L,K))  
      VTERM = AVHV*( POS(L,K,IT) - POS(LS,K, IT) )  
      IF( VHDX2(L,K) >= 0.0 )THEN
        ! *** QUANTITIES FROM SOUTH CELL
        VTERM = VTERM - 0.5*DDELTA*VHDX2(L,K)*( UUUU(LS,K,IT)+UUUU(LSE,K,IT) + WWWW(LS,K,IT)+WWWW(LS,K-1,IT) + VVVV(LS,K,IT)+VVVV(L,K,IT) )  
      ELSE  
        ! *** QUANTITIES FROM CURRENT CELL
        VTERM = VTERM - 0.5*DDELTA*VHDX2(L,K)*( UUUU(L,K,IT) +UUUU(LE,K,IT)  + WWWW(L,K,IT) +WWWW(L ,K-1,IT) + VVVV(LN,K,IT)+VVVV(L,K,IT) )  
      ENDIF  
      VHV = VTERM/( POS(L,K,IT) + POS(LS,K ,IT) + BSMALL )  
      FVHUD(L,K,IW) = MAX(VHV,0.)*POS(LS,K ,IT) + MIN(VHV,0.)*POS(L,K,IT)  
    ENDDO  
  ENDDO  

  IF( KC > 1 )THEN
    DO K=1,KS
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)
        LE = LEC(L)
        LN = LNC(L)  
        AWW   = ABS(W2(L,K))  
        WTERM = AWW*(POS(L,K+1,IT)-POS(L,K,IT))  
        IF( W2(L,K) >= 0.0 )THEN  
          WTERM = WTERM-0.5*DDELTA*W2(L,K)*(UUUU(L,K,IT)   + UUUU(LE,K,IT)   + VVVV(L,K,IT)   + VVVV(LN,K,IT)   + WWWW(L,K,IT)+WWWW(L,K-1,IT))
        ELSE  
          WTERM = WTERM-0.5*DDELTA*W2(L,K)*(UUUU(L,K+1,IT) + UUUU(LE,K+1,IT) + VVVV(L,K+1,IT) + VVVV(LN,K+1,IT) + WWWW(L,K,IT)+WWWW(L,K+1,IT))
        ENDIF  
        WW = WTERM/( POS(L,K+1,IT) + POS(L,K,IT) + BSMALL )  
        FWUU(L,K,IW) = MAX(WW,0.)*POS(L,K,IT) + MIN(WW,0.)*POS(L,K+1,IT)  
      ENDDO  
    ENDDO  
  ENDIF

  ! ----------------------------------------------------------------------------------------
  ! ** ZERO ANTIDIFFUSIVE FLUXES FOR SELECTED CELLS
  
  ! ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR SOURCE CELLS  
  IF( ISADAC(MVAR) == 1 )THEN  
    ! *** ANTI-DIFFUSION TURNED OFF FOR SOURCE CELLS
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        IF( ABS(QSUM(L,K)) > 1.E-8 )THEN
          FUHUD(L     , K,IW) = 0.  
          FUHUD(LEC(L), K,IW) = 0.  
          FVHUD(L     , K,IW) = 0.  
          FVHUD(LNC(L), K,IW) = 0.  
          FWUU(L,K  ,IW)  = 0.  
          FWVV(L,K-1,IT)  = 0.  
        ENDIF  
      ENDDO  
    ENDDO  
  ENDIF  
  
  ! ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR OPEN BOUNDARY CELLS  
  DO K=1,KC  
    DO LL=1,NCBS  
      L = LCBS(LL)  
      FVHUD(LNC(L),K,IW) = 0.0  
    ENDDO  
    DO LL=1,NCBW  
      L = LCBW(LL)  
      FUHUD(LEC(L),K,IW) = 0.0  
    ENDDO  
    DO LL=1,NCBE  
      L = LCBE(LL)  
      FUHUD(L,     K,IW) = 0.0  
    ENDDO  
    DO LL=1,NCBN  
      L = LCBN(LL)  
      FVHUD(L,     K,IW) = 0.0  
    ENDDO  
  ENDDO  
  
  ! ----------------------------------------------------------------------------------------
  ! **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS  
  IF( ISFCT(MVAR) /= 0 )THEN  
  
    ! **  DETERMINE MAX AND MIN CONCENTRATIONS  
    DO K=1,KC
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CONTMX(L,K,IT) = MAX(POS(L,K,IT),ABS(CON2(L,K,IW)))  
        CONTMN(L,K,IT) = MIN(POS(L,K,IT),ABS(CON2(L,K,IW)))  
      ENDDO  
    ENDDO

    DO LP=1,LAWET
      L = LWET(LP) 
      CMAX(L,KSZ(L),IT) = MAX(CONTMX(L,KSZ(L),IT), CONTMX(L,KSZ(L)+1,IT))  
      CMAX(L,KC    ,IT) = MAX(CONTMX(L,KS    ,IT), CONTMX(L,KC      ,IT))  
      CMIN(L,KSZ(L),IT) = MIN(CONTMN(L,KSZ(L),IT), CONTMN(L,KSZ(L)+1,IT))  
      CMIN(L,KC    ,IT) = MIN(CONTMN(L,KS    ,IT), CONTMN(L,KC      ,IT))  
    ENDDO  
    
    DO K=2,KS  
      DO LP=1,LLWET(K-1,0)
        L = LKWET(LP,K-1,0)  
        CMAXT        = MAX(CONTMX(L,K-1,IT),CONTMX(L,K+1,IT))
        CMAX(L,K,IT) = MAX(CONTMX(L,K,IT),CMAXT)  
        CMINT        = MIN(CONTMN(L,K-1,IT),CONTMN(L,K+1,IT))  
        CMIN(L,K,IT) = MIN(CONTMN(L,K,IT),CMINT)  
      ENDDO  
    ENDDO  
    
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        LE = LEC(L)
        LS = LSC(L)  
        LN = LNC(L)  
        LW = LWC(L)
        
        CWMAX = SUB3D(L,K) *CONTMX(LW,K,IT)  
        CEMAX = SUB3D(LE,K)*CONTMX(LE,K,IT)  
        CSMAX = SVB3D(L,K) *CONTMX(LS,K,IT)  
        CNMAX = SVB3D(LN,K)*CONTMX(LN,K,IT)  
        CMAXT = MAX(CNMAX,CEMAX)  
        CMAXT = MAX(CMAXT,CSMAX)  
        CMAXT = MAX(CMAXT,CWMAX)  
        CMAX(L,K,IT) = MAX(CMAX(L,K,IT),CMAXT)  
        
        CWMIN = SUB3D(L,K) *CONTMN(LW,K,IT) + 1.E+6*(1.-SUB3D(L,K))
        CEMIN = SUB3D(LE,K)*CONTMN(LE,K,IT) + 1.E+6*(1.-SUB3D(LE,K))  
        CSMIN = SVB3D(L,K) *CONTMN(LS,K,IT) + 1.E+6*(1.-SVB3D(L,K))
        CNMIN = SVB3D(LN,K)*CONTMN(LN,K,IT) + 1.E+6*(1.-SVB3D(LN,K))  
        CMINT = MIN(CNMIN,CEMIN)  
        CMINT = MIN(CMINT,CSMIN)  
        CMINT = MIN(CMINT,CWMIN)  
        CMIN(L,K,IT) = MIN(CMIN(L,K,IT),CMINT)  
      ENDDO  
    ENDDO  
  
    ! **  SEPARATE POSITIVE AND NEGATIVE FLUXES  
    ! **  PUT NEGATIVE FLUXES INTO FUHVD, FVHVD, AND FWVV  
    DO K=1,KC
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        FUHVD(L,K,IT) = MIN(FUHUD(L,K,IW),0.)
        FUHUD(L,K,IW) = MAX(FUHUD(L,K,IW),0.)
        FVHVD(L,K,IT) = MIN(FVHUD(L,K,IW),0.)
        FVHUD(L,K,IW) = MAX(FVHUD(L,K,IW),0.)  
      ENDDO  
    ENDDO

    IF( KC > 1 )THEN  
      DO K=1,KS
        DO LP=1,LLWET(K,0)
          L = LKWET(LP,K,0)  
          FWVV(L,K,IT) = MIN(FWUU(L,K,IW),0.)  
          FWUU(L,K,IW) = MAX(FWUU(L,K,IW),0.)  
        ENDDO  
      ENDDO  
    ENDIF

    ! **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS  
    ! **  LOAD INTO DUU AND DVV, THEN ZERO VALUES AT BOUNDARIES  
    DO K=1,KC
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        LE = LEC(L)
        LN = LNC(L)
        ! ***                              MAX              MIN              MAX              MIN                  MAX              MIN
        DUU(L,K,IT) = DDELT*( DXYIP(L)*( FUHUD(L,K,IW)  - FUHVD(LE,K,IT) + FVHUD(L,K,IW)  - FVHVD(LN,K,IT) ) + ( FWUU(L,K-1,IW) - FWVV(L,K,IT))   )*HPKI(L,K)
        DVV(L,K,IT) = DDELT*( DXYIP(L)*( FUHUD(LE,K,IW) - FUHVD(L,K,IT)  + FVHUD(LN,K,IW) - FVHVD(L,K,IT)  ) + ( FWUU(L,K,IW)   - FWVV(L,K-1,IT)) )*HPKI(L,K)
      ENDDO  
    ENDDO  

    ! *** ZERO SELECTED FLUX CORRECTORS
    DO K=1,KC  
      DO IOBC=1,NBCSOP  
        L=LOBCS(IOBC)  
        DUU(L,K,IT) = 0.  
        DVV(L,K,IT) = 0.  
      ENDDO  
    END DO  
  
    ! **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DUU AND DVV  
    DO K=1,KC
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        IF( DUU(L,K,IT) > 0. ) DUU(L,K,IT) = (CMAX(L,K,IT)-POS(L,K,IT))/(DUU(L,K,IT)+BSMALL)  
        DUU(L,K,IT) = MIN(DUU(L,K,IT),1.)  
        IF( DVV(L,K,IT) > 0. ) DVV(L,K,IT) = (POS(L,K,IT)-CMIN(L,K,IT))/(DVV(L,K,IT)+BSMALL)  
        DVV(L,K,IT) = MIN(DVV(L,K,IT),1.)  
      ENDDO  
    ENDDO  
  
    ! **  LIMIT FLUXES  
    DO K=1,KC
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        LS = LSC(L)
        LW = LWC(L)
        FUHUD(L,K,IW) = MIN(DVV(LW,K,IT),DUU(L,K,IT))*FUHUD(L,K,IW) + MIN(DUU(LW,K,IT),DVV(L,K,IT))*FUHVD(L,K,IT)  
        FVHUD(L,K,IW) = MIN(DVV(LS,K,IT),DUU(L,K,IT))*FVHUD(L,K,IW) + MIN(DUU(LS,K,IT),DVV(L,K,IT))*FVHVD(L,K,IT)  
      ENDDO  
    ENDDO  

    DO K=1,KS
      DO LP=1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        FWUU(L,K,IW) = MIN(DVV(L,K,IT),DUU(L,K+1,IT))*FWUU(L,K,IW) + MIN(DUU(L,K,IT),DVV(L,K+1,IT))*FWVV(L,K,IT)  
      ENDDO  
    ENDDO  
    
    ! *** APPLY OPEN BOUNDARYS 
    DO LL=1,NBCSOP
      L=LOBCS(LL)
      DO K=1,KS  
        FWUU(L,K,IW) = 0.0
      ENDDO  
    ENDDO 
    DO LL=1,NBCSOP2
      L=LOBCS2(LL)
      DO K=1,KS  
        FWUU(L,K,IW) = 0.0
      ENDDO  
    ENDDO 

  ENDIF  ! *** END OF FLUX CORRECTOR SECTION
  
  ! *** APPLY THE ANTI-DIFFUSIVE ADVECTION CALCULATION
  ! *** Moved to CALCONC
  
  ! ----------------------------------------------------------------------------------------

  RETURN  
END  
