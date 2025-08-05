! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALTRANICE(CON, CON1, IT)  

  ! ***  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE  
  ! ***  TRANSPORT OF FRAZIL ICE OR ANY OTHER ICE IMPACTED TRANSPORT
  ! ***  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES  
  ! ***  THE NUMBER OF TIME LEVELS IN THE STEP  

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------------!
  ! *** 2015-01     PAUL M. CRAIG      Added fully coupled Ice Sub-model with Frazil Ice Transport
  ! *** 2014-09     PAUL M. CRAIG      ADDED THE LWET BYPASS APPROACH

  use GLOBAL  
  !  use OMP_LIB
  
  implicit none  
    
  ! *** Passed in variables
  integer, intent(IN) :: IT  
  real, intent(INOUT)    :: CON(LCM,KCM),CON1(LCM,KCM)  
  ! *** Local variables
  real,save,allocatable,dimension(:,:) :: FLUX
  real    :: DDELT  
  real    :: RDZIC, CTMP, WVEL, CLEFT, CRIGHT 
  
  integer :: LP,L, LN, LS, LL, K
  
  if( .not. allocated(FLUX) )then
    allocate(FLUX(LCM,KCM))
    FLUX = 0.0
  endif

  !ISUD = 1  
  if( ISDYNSTP == 0 )then
    ! *** FIXED DELTA T
    if( ISTL == 3 )then  
      DDELT = DT2
    else
      DDELT = DT
    endif
  else  
    DDELT = DTDYN   ! *** DYNAMIC DELTA T
  endif  
    
 ! *** ADVANCE THE CONSTITUENT
  CON1 = CON 
    
  ! *** SKIP BOUNDARY LOADINGS FOR FRAZIL TRANSPORT (MAY BE ADDED FOR FUTURE VERSIONS)
  !IF( .FALSE. )then
  !  ! *** save OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  !  do IOBC = 1,NBCSOP  
  !    L = LOBCS(IOBC)  
  !    do K = 1,KC  
  !      WQBCCON1(IOBC,K,IT) = CON1(L,K)  
  !      WQBCCON(IOBC,K,IT)  = CON(L,K)  
  !    enddo  
  !  enddo  
  !
  !  ! ***  CALCULATED EXTERNAL SOURCES AND SINKS  
  !  call CALFQC (ISTL,IS2TL_,MVAR,10,CON,CON1,IT)  
  !ENDIF

  ! ***  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! ***  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED  
  ! ***  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY  

  ! *** COMPUTE FLUXES
  do K = 1,KC
    do LP = 1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      FUHUD(L,K,IT) = UHDY2(L,K)*CON1(LUPU(L,K),K)
      FVHUD(L,K,IT) = VHDX2(L,K)*CON1(LUPV(L,K),K)  
    enddo  
  enddo
  
  if( KC > 1 )then  
    do K = 1,KS  
      do LP = 1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        FWUU(L,K,IT) = W2(L,K)*CON1(L,KUPW(L,K))  
      enddo  
    enddo  
  endif  
  
  ! *** CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX (PMC MOVED)  
  !IF(ISHDMF == 2 ) CALL CALDIFF (CON1,IT)  
    
  ! *** Upwind Differencing (3TL & 2TL)  
  if( ISTL == 2 )then  
    do K = 1,KC  
      RDZIC = DZIC(L,K)  
      do LP = 1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CD(L,K,IT) = CON1(L,K)*H1P(L)+DDELT*( ( FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT))*DXYIP(L) &
                                          + ( FWUU(L,K-1,IT)-FWUU(L,K,IT))*RDZIC )  
      enddo  
    enddo    
  else  ! *** IF ISTL = 3
    do K = 1,KC  
      RDZIC = DZIC(L,K)  
      do LP = 1,LLWET(K,0)
        L = LKWET(LP,K,0)  
        CD(L,K,IT) = CON1(L,K)*H2P(L)+DDELT*( ( FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT))*DXYIP(L) &
                                          + ( FWUU(L,K-1,IT)-FWUU(L,K,IT))*RDZIC )  
      enddo  
    enddo
  endif
  
  ! *** UPDATE NEW CONCENTRATIONS  
  do K = 1,KC  
    do LP = 1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      CON(L,K) = CD(L,K,IT)*HPI(L) 
    enddo

    ! *** RESET OPEN BC CONCENTRATIONS  
    !DO IOBC = 1,NBCSOP  
    !  L = LOBCS(IOBC)  
    !  CON(L,K) = WQBCCON(IOBC,K,IT)  
    !ENDDO  
  enddo  

  !IF( ISUD == 1 )then
  !  ! *** RESTORE ORIGINAL CONCENTRATIONS PRIOR TO APPLYING OPEN BC'S  
  !  do K = 1,KC  
  !    do IOBC = 1,NBCSOP  
  !      L = LOBCS(IOBC)  
  !      CON1(L,K) = WQBCCON(IOBC,K,IT)  
  !    enddo  
  !  enddo  
  !ENDIF  

  ! ******************************************************************************************
  ! *** COMPUTE THE FRAZIL ICE RISE
  
  ! *** BOTTOM LAYER
  do LP = 1,LAWET
    K = KSZ(L)
    L = LWET(LP) 
    WVEL   = DDELT*HPI(L)*DZIC(L,K)
    CLEFT  = 1.+RISEVEL*WVEL
    CRIGHT = CON(L,K)
    CON(L,K) = CRIGHT/CLEFT
    FLUX(L,K) = RISEVEL*CON(L,K) 
  enddo

  ! *** ALL OTHER LAYERS
  do K = 2,KC
    do LP = 1,LLWET(K,0)
      L = LKWET(LP,K,0)  
      WVEL   = DDELT*HPI(L)*DZIC(L,K)
      CLEFT  = 1.+RISEVEL*WVEL
      CRIGHT = CON(L,K)+FLUX(L,K-1)*WVEL
      CON(L,K) = CRIGHT/CLEFT
      FLUX(L,K) = RISEVEL*CON(L,K)
    enddo
  enddo

  ! *** ACCUMULATE ICE COVER
  do LP = 1,LAWET
    L = LWET(LP)
    ICETHICK(L) = ICETHICK(L) + FLUX(L,KC)*DDELT/(HP(L)*DZC(L,KC))
  enddo
  
  ! ******************************************************************************************
  ! *** APPLY OPEN BOUNDARY CONDITIONS, BASED ON DIRECTION OF FLOW  

  ! *** SOUTH OPEN BC
  if( NCBS > 0 )then
    do K = 1,KC  
      do LL = 1,NCBS  
        L = LCBS(LL)  
        LN = LNC(L)  
        if( VHDX2(LN,K) <= 0. )then  
          ! *** FLOWING OUT OF DOMAIN  
          if( ISTL == 2 )then  
            CTMP = CON1(L,K)+DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPI(L)  
          else  
            CTMP = CON1(L,K)+DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPI(L)  
            CON1(L,K) = CON(L,K)  
          endif  
          CON(L,K) = MAX(CTMP  ,0.)
        endif  
      enddo  
    enddo  
  endif    

  ! *** WEST OPEN BC 
  if( NCBW > 0 )then
    do K = 1,KC  
      do LL = 1,NCBW  
        L = LCBW(LL)  
        if( UHDY2(LEC(L),K) <= 0. )then  
          ! *** FLOWING OUT OF DOMAIN  
          if( ISTL == 2 )then  
            CTMP = CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPI(L)  
          else  
            CTMP = CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPI(L)  
            CON1(L,K) = CON(L,K)  
          endif  
          CON(L,K) = MAX(CTMP  ,0.)
        endif  
      enddo  
    enddo  
  endif    

  ! *** EAST OPEN BC
  if( NCBE > 0 )then
    do K = 1,KC  
      do LL = 1,NCBE  
        L = LCBE(LL)  
        if( UHDY2(L,K) >= 0. )then  
          ! *** FLOWING OUT OF DOMAIN  
          if( ISTL == 2 )then  
            CTMP = CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
          else  
            CTMP = CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
            CON1(L,K) = CON(L,K)  
          endif  
          CON(L,K) = MAX(CTMP  ,0.)
        endif  
      enddo  
    enddo  
  endif
    
  ! *** NORTH OPEN BC 
  if( NCBN > 0 )then
    do K = 1,KC  
      do LL = 1,NCBN  
        L = LCBN(LL)  
        LS = LSC(L)  
        if( VHDX2(L,K) >= 0. )then  
          ! *** FLOWING OUT OF DOMAIN  
          if( ISTL == 2 )then  
            CTMP = CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
          else  
            CTMP = CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
            CON1(L,K) = CON(L,K)  
          endif  
          CON(L,K) = MAX(CTMP  ,0.)
        endif  
      enddo  
    enddo  
  endif

  ! ----------------------------------------------------------------------------------------
  ! *** CALTRANICE EXIT 

  !  $ print *, ' CALTRAN-Leaving: Thread, MVAR, MO',IT,MVAR,MO  

  return 
END  
