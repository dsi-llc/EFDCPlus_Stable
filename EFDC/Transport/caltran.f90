! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @DETAILS DESCRIBE THE SUBROUTINE OR WHATEVEVER
! @DATE 10/22020
! @AUTHOR PAUL 
  
  SUBROUTINE CALTRAN (MVAR, MO, CON, CON1, IW, IT, WCCUTOFF, ISKIP)  

  ! ***  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE  
  ! ***  TRANSPORT OF DISSOLVED OR SUSPENDED CONSTITUENT M LEADING TO  
  ! ***  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES  
  ! ***  THE NUMBER OF TIME LEVELS IN THE STEP  

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2015-01       PAUL M. CRAIG     Added fully coupled Ice Sub-model with Frazil Ice Transport
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH

  use GLOBAL 
  
  implicit none  
    
  ! *** Passed in variables
  integer, intent(IN)    :: MVAR, MO, IW, IT
  integer, intent(INOUT) :: ISKIP
  real, intent(IN)       :: WCCUTOFF
  real, intent(INOUT)    :: CON(LCM,KCM), CON1(LCM,KCM)  
  
  ! *** Local variables
  real    :: DDELT, DDELTA, DDELTD2, WCSUM, BCSUM, WCVMAX
  real    :: CTMP, CBT, AUHU, AVHV, UTERM, VTERM, WTERM  
  real    :: CBSTMP, CBWTMP, CBETMP, CBNTMP, UHU, VHV, AWW, WW  
  real    :: CWMAX, CEMAX, CSMAX, CNMAX, CMAXT  
  real    :: CWMIN, CEMIN, CSMIN, CNMIN, CMINT  
  
  integer :: M, iunit, II, JJ, I, J, NS, ITRANFLOC  
  integer :: ISUD, K, KMAX, NSID, IOBC, NMNLOD, NCELL, NPAR
  integer :: LP, L, LN, LS, LE, LW, LSE, LNW, LL, LMAX
  
  ! *** Zero any initialized concentrations below bottom active layer
  if( IGRIDV > 0 .and. N < 5 )then
    do L = 2,LA
      do K = 1,KSZ(L)-1
        CON(L,K)  = 0.0
        CON1(L,K) = 0.0
      enddo
    enddo
  endif
  
  ! *** Set up floc transport
  ITRANFLOC = 0
  
  ISUD = 1  
  if( ISDYNSTP == 0 )then 
    ! *** FIXED DELTA T
    DDELT = DT2  
    DDELTA = DT2  
    DDELTD2 = DT
    if( ISTL /= 3 )then  
      DDELT = DT  
      DDELTA = DT  
      DDELTD2 = 0.5*DT  
      if( IS2TIM == 0 ) ISUD = 0          ! *** 3TL CORRECTOR TIME STEP (ISTL = 2)
    endif  
  else  
    ! *** Dynamic delta T
    DDELT = DTDYN  
    DDELTA = DTDYN  
    DDELTD2 = 0.5*DTDYN  
  endif  
 
  M = MO  
  
  if( IS2TL == 1 )then  
    ! *** ADVANCE CONCENTRATIONS BEFORE THE 2TL TRANSPORT CALCULATIONS
    ! *** SKIP UPDATING VARIABLES IF ALREADY COMPLETED BEFORE THIS STEP
    if( MVAR /= 8 )then
      CON1(:,:) = CON(:,:)
    endif  
  endif  
  
  if( IW <= NACTIVEWC )then
    ! *** save OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
    do IOBC = 1,NBCSOP  
      L = LOBCS(IOBC)  
      do K = 1,KC  
        WQBCCON(IOBC,K,IW)  = CON(L,K)  
        WQBCCON1(IOBC,K,IW) = CON1(L,K)  
      enddo  
    enddo  
  
    ! ***  CALCULATED EXTERNAL SOURCES AND SINKS  
    call CALFQC(MVAR, M, CON, CON1, IT)  
  endif
  
  ! *** SKIP TRANSPORT IF NOTHING TO TRANSPORT
  ISKIP = 0
  WCVMAX = 0.0
  if( WCCUTOFF > 0.0 )then
    
    ! *** SUM ALL DEFINED BCs
    BCSUM = 0.0
    do NS = 1,NBCS
      L = LBCS(NS)
      do K = 1,KC
        BCSUM = BCSUM + FQC(L,K,IT)
      enddo
    enddo
    
    if( BCSUM <= 0.0 )then
      ! *** No loadings.  See if WC has any constituent to transport
      ISKIP = 1
      do K = 1,KC  
        do LP = 1,LLWET(K,0)
          L = LKWET(LP,K,0)
          WCVMAX = max(WCVMAX, CON(L,K))
          if( WCVMAX >= WCCUTOFF )then
            ISKIP = 0
            exit
          endif
        enddo  
        if( ISKIP == 0 ) exit
      enddo
      
      ! *** If minimum concentration not found.  Skip the transport for this consituent 
      if( ISKIP == 1 )then
        !PRINT '(3I6,2E12.4)', NITER, MVAR, MO, WCVMAX, WCCUTOFF   ! DELME
        return
      endif
    endif
  endif
      
  ! *****************************************************************************************
  ! *** Begin constituent transport
  
  ! *** Upwind Differencing
  do K = 1,KC  
    do LP = 1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      FUHUD(L,K,IW) = UHDY2(L,K)*CON1(LUPU(L,K),K)  
      FVHUD(L,K,IW) = VHDX2(L,K)*CON1(LUPV(L,K),K)  
    enddo  
  enddo
  if( KC > 1 )then  
    do K = 1,KS  
      do LP = 1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        FWUU(L,K,IW) = W2(L,K)*CON1(L,KUPW(L,K))  
      enddo
    enddo
  endif  
 
  ! *** Calculate and add horizontal diffusion flux
  if( (ISHDMF == 2 .or. ISHDMF == 4) .and. IW <= NACTIVEWC ) CALL CALDIFF (CON1,IW)
    
  ! *** Upwind Differencing (3TL & 2TL)  
  if( ISTL == 2 )then  
    do K = 1,KC
      do LP = 1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CD(L,K,IT) = CON1(L,K)*H1PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                     FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                  + (FWUU(L,K-1,IW)-FWUU(L,K,IW)) ) 
      enddo  
    enddo

    if( ISFCT(MVAR) >= 1 .and. ISADAC(MVAR) > 0 .and. IW <= NACTIVEWC )then 
      CON2(:,:,IW) = max(CON1(:,:),0.0)                      ! *** Save previous concentration for flux corrector
    endif  
    
  else  ! *** IF ISTL = 3
    do K = 1,KC  
      do LP = 1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CD(L,K,IT) = CON1(L,K)*H2PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                     FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                  + (FWUU(L,K-1,IW)-FWUU(L,K,IW))  )  
      enddo  
    enddo
      
    if( ISFCT(MVAR) >= 1 .and. ISADAC(MVAR) > 0 .and. IW <= NACTIVEWC )then
      CON2(:,:,IW) = max(CON(:,:),0.0)
    endif  
      
  endif   ! *** For ISTL = 2 and ISTL = 3
  
  if( IS2TL == 0 .and. ISUD == 1 )then  
    ! *** ADVANCE CON1 TO CON (3TL)
    CON1(:,:) = CON(:,:)
  endif  

  ! *** UPDATE NEW CONCENTRATIONS  
  do K = 1,KC  
    do LP = 1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      CON(L,K) = CD(L,K,IT)*HPKI(L,K)
    enddo
  enddo

  ! *** OPEN BOUNDARIES.  BYPASS FOR TURBULENCE ADVECTION (IW = 0)
  if( IW <= NACTIVEWC )then
    ! *** RESTORE ORIGINAL CONCENTRATIONS PRIOR TO APPLYING OPEN BC'S - 2TL & 3TL
    do IOBC = 1,NBCSOP
      L = LOBCS(IOBC)
      do K = 1,KC
        CON1(L,K) = WQBCCON1(IOBC,K,IW)
      enddo
    enddo

    ! ******************************************************************************************
    ! *** Apply open boundary conditions, based on direction of flow
    
    ! *** SOUTH OPEN BC
    if( NCBS > 0 )then
      do K = 1,KC
        do LL = 1,NCBS
          NSID = NCSERS(LL,MVAR)
          L = LCBS(LL)
          LN = LNC(L)
          if( LKSZ(L,K) .or. .not. LMASKDRY(L) ) CYCLE
          
          if( VHDX2(LN,K) <= 0. )then
            ! *** Flowing out of domain
            CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K) - FVHUD(LN,K,IW))*DXYIP(L)*HPKI(L,K)
            CON(L,K) = max(CTMP  ,0.)

            CLOS(LL,K,M) = CON(L,K)
            NLOS(LL,K,M) = NITER
          else
            ! *** Flowing into domain
            CBT = WTCI(K,1)*CBS(LL,1,M) + WTCI(K,2)*CBS(LL,2,M) + CSERT(K,NSID,M)
            
            NMNLOD = NITER - NLOS(LL,K,M)
            if( NMNLOD >= NTSCRS(LL) )then
              CON(L,K) = CBT
            else
              CBSTMP = CLOS(LL,K,M) + (CBT-CLOS(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRS(LL))
              CON(L,K) = max(CBSTMP,0.)
            endif
          endif
          if( ISUD == 1 ) CON1(L,K) = CON(L,K)
        enddo
      enddo
    endif

    ! *** WEST OPEN BC
    if( NCBW > 0 )then
      do K = 1,KC
        do LL = 1,NCBW
          NSID = NCSERW(LL,MVAR)
          L = LCBW(LL)
          LE = LEC(L)
          if( LKSZ(L,K) .or. .not. LMASKDRY(L) ) CYCLE
          
          if( UHDY2(LE,K) <= 0. )then
            ! *** Flowing out of domain
            CTMP = CON1(L,K) + DDELT*(UHDY2(LE,K)*CON1(L,K) - FUHUD(LE,K,IW))*DXYIP(L)*HPKI(L,K)
            CON(L,K) = max(CTMP  ,0.)
            
            CLOW(LL,K,M) = CON(L,K)
            NLOW(LL,K,M) = NITER
          else
            ! *** Flowing into domain
            CBT = WTCI(K,1)*CBW(LL,1,M) + WTCI(K,2)*CBW(LL,2,M) + CSERT(K,NSID,M)

            NMNLOD = NITER - NLOW(LL,K,M)
            if( NMNLOD >= NTSCRW(LL) )then
              CON(L,K) = CBT
            else
              CBWTMP = CLOW(LL,K,M)+(CBT-CLOW(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRW(LL))
              CON(L,K) = max(CBWTMP,0.)
            endif
          endif
          if( ISUD == 1 ) CON1(L,K) = CON(L,K)
        enddo
      enddo
    endif

    ! *** EAST OPEN BC
    if( NCBE > 0 )then
      do K = 1,KC
        do LL = 1,NCBE
          NSID = NCSERE(LL,MVAR)
          L = LCBE(LL)
          if( LKSZ(L,K) .or. .not. LMASKDRY(L) ) CYCLE
          
          if( UHDY2(L,K) >= 0. )then
            ! *** Flowing out of domain
            CTMP = CON1(L,K) + DDELT*(FUHUD(L,K,IW) - UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
            CON(L,K) = max(CTMP  ,0.)

            CLOE(LL,K,M) = CON(L,K)
            NLOE(LL,K,M) = NITER
          else
            ! *** Flowing into domain
            CBT = WTCI(K,1)*CBE(LL,1,M) + WTCI(K,2)*CBE(LL,2,M) + CSERT(K,NSID,M)
            
            NMNLOD = NITER - NLOE(LL,K,M)
            if( NMNLOD >= NTSCRE(LL) )then
              CON(L,K) = CBT
            else
              CBETMP = CLOE(LL,K,M) + (CBT-CLOE(LL,K,M)) * FLOAT(NMNLOD)/FLOAT(NTSCRE(LL))
              CON(L,K) = max(CBETMP,0.)
            endif
          endif
          if( ISUD == 1 ) CON1(L,K) = CON(L,K)
        enddo
      enddo
    endif

    ! *** NORTH OPEN BC
    if( NCBN > 0 )then
      do K = 1,KC
        do LL = 1,NCBN
          NSID = NCSERN(LL,MVAR)
          L = LCBN(LL)
          if( LKSZ(L,K) .or. .not. LMASKDRY(L) ) CYCLE

          if( VHDX2(L,K) >= 0. )then
            ! *** Flowing out of domain
            CTMP = CON1(L,K) + DDELT*(FVHUD(L,K,IW) - VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
            CON(L,K) = max(CTMP,0.)
            
            CLON(LL,K,M) = CON(L,K)
            NLON(LL,K,M) = NITER
          else
            ! *** Flowing into domain
            CBT = WTCI(K,1)*CBN(LL,1,M) + WTCI(K,2)*CBN(LL,2,M) + CSERT(K,NSID,M)

            NMNLOD = NITER - NLON(LL,K,M)
            if( NMNLOD >= NTSCRN(LL) )then
              CON(L,K) = CBT
            else
              CBNTMP = CLON(LL,K,M)+(CBT-CLON(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRN(LL))
              CON(L,K) = max(CBNTMP,0.)
            endif
          endif
          if( ISUD == 1 ) CON1(L,K) = CON(L,K)
        enddo
      enddo
    endif
  endif    ! *** END OF IW > 0 BLOCK
  
  ! ****************************************************************************************
 
  ! *** ZERO HEAT FLUXES
  if( MVAR == 2 )then        
    ! *** ZERO EVAP/RAINFALL
    do L = 1,LC  
      FQC(L,KC,IT) = 0.  
    enddo  
  endif
  
  !  $ print *, ' CALTRAN-Leaving: Thread, MVAR, MO',IT,MVAR,MO  

  ! *** ZERO HEAT FLUXES
  if( MVAR == 2 )then        
    ! *** ZERO EVAP/RAINFALL
    !DO L = 1,LC  
    !  FQC(L,KC,IT) = 0.  
    !ENDDO  
  endif
  
  return 
END  
