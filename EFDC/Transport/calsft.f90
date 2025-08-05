! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALSFT
  !
  ! CHANGE RECORD
  ! ***  SUBROUTINE CALSFT CALCULATES THE TRANSPORT OF SHELL FISH LARVAE
  ! ***  AT TIME LEVEL (N+1).
  ! ***  CALLED ONLY ON ODD THREE TIME LEVEL STEPS  (PMC - NO, CALLED IN BOTH HDMT & HDMT2T)
  !

  use GLOBAL
  use Variables_WQ,ONLY : HWQ, H2WQ
  
  ! *** DSI BEGIN BLOCK
  implicit none

  integer :: ISDARK,ITIME,K,L,LN,LF,LL,ND,ISKIP
  real :: TIME,RTIME,TIMTMP,FANGTMP,UTMP,VTMP,VELEKB,VELNKB,CURANG,ANGDIF,RCDZKMK,CCLBTMP
  real :: RABOVE,HABOVE,TMPVAL,WMAXX,DZCIT,TMPVAL1,CDYETMP,RCDZKK,CCUBTMP,CCMBTMP,EEB

  real,save,allocatable,dimension(:)   :: HWQI
  real,save,allocatable,dimension(:)   :: WTFKB
  real,save,allocatable,dimension(:)   :: WTFKC
  real,save,allocatable,dimension(:,:) :: WWQ         

  if( .not. allocated(WTFKB) )then
    allocate(HWQI(LCM))
    allocate(WTFKB(KCM))
    allocate(WTFKC(KCM))
    allocate(WWQ(LCM,0:KCM))
    ! *** ZERO LOCAL ARRAYS
    HWQI  = 0.0
    WTFKB = 0.0
    WTFKC = 0.0
    WWQ   = 0.0
  endif
  ! *** DSI END BLOCK
  !
  !PMC      DELT = DT2
  ! *** PMC
  if( ISTL == 2 )then
    if( ISDYNSTP == 0 )then
      DELT = DT
    else
      DELT = DTDYN
    endif
  else
    DELT = DT2
  endif
  ! *** PMC
  !
  ! ***  UPDATED TIME SERIES CONCENTRATION BOUNDARY CONDITIONS
  ! ***  DETERMINE IF CURRENT TIME STEP IS DURING DAYLIGHT OR DARKNESS
  !
  if( ISSFLDN >= 1 )then
    ISDARK = 1
    TIME = TIMESEC/86400.
    ITIME = INT(TIME)
    RTIME = FLOAT(ITIME)
    TIMTMP = TIME-RTIME
    if( TIMTMP >= TSRSF .and. TIMTMP <= TSSSF) ISDARK = 0
  endif
  !
  ! ***  DETERMINE IF LOCAL CONDITIONS ARE EBB OR FLOOD
  !
  if( ISSFLFE >= 1 )then
    if( KC == 1 )then
      WTFKB(1) = 1.
      WTFKC(1) = 0.
    endif
    if( KC == 2 )then
      WTFKB(1) = 1.0
      WTFKC(1) = 0.0
      WTFKB(2) = 0.0
      WTFKC(2) = 1.0
    endif
    if( KC == 3 )then
      do K = 1,KC
        WTFKB(K) = FLOAT(KC-K)/FLOAT(KS)
        WTFKC(K) = 1.0-WTFKB(K)
      enddo
    endif
  !
  ! ***  SET SWITCHES TO EBB
  !
    do K = 1,KC
      do L = 2,LA
        UUU(L,K) = 0.
        VVV(L,K) = 1.
      enddo
    enddo
  !
  ! ***  RESET SWITCHES FOR FLOOD
  !
    do K = 1,KC
      do L = 2,LA
        LN = LNC(L)
        FANGTMP = ACCWFLD(L,1)*WTFKB(K)+ACCWFLD(L,2)*WTFKC(K)
        UTMP = 0.5*STCUV(L)*(U(LEC(L),K)+U(L,K))
        VTMP = 0.5*STCUV(L)*(V(LN ,K)+V(L,K))
        VELEKB = CUE(L)*UTMP+CVE(L)*VTMP+1.E-12
        VELNKB = CUN(L)*UTMP+CVN(L)*VTMP
        CURANG = ATAN2(VELNKB,VELEKB)
        ANGDIF = ABS(FANGTMP-CURANG)
        if( ANGDIF < 1.5708 )then
          UUU(L,K) = 1.
          VVV(L,K) = 0.
        endif
      enddo
    enddo
  endif
  !
  ! ***  SET UP ADVECTION FIELD
  ! ***  SET ATTACHED TO BOTTOM AND NO ADVECTIVE TRANSPORT IN BOTTOM
  ! ***  LAYER DURING EBB IF APPROPRIATE
  !
  if( ISSFLFE >= 1 )then
    if( SFNTBET < 1. )then
      do L = 2,LA
        K = KSZ(L)
        UHDY2(L,K) = UUU(L,K)*UHDY2(L,K)+SFNTBET*VVV(L,K)*UHDY2(L,K)
        VHDX2(L,K) = UUU(L,K)*VHDX2(L,K)+SFNTBET*VVV(L,K)*VHDX2(L,K)
        U2(L,K) = UUU(L,K)*U2(L,K)+SFNTBET*VVV(L,K)*U2(L,K)
        V2(L,K) = UUU(L,K)*V2(L,K)+SFNTBET*VVV(L,K)*V2(L,K)
      enddo
    endif
  endif

  ! *** COMPUTE SHELLFISH LARVAE ADVECTION
  call CALTRAN (4,4,SFL,SFL2,1,1,0.0,ISKIP)
  !CALL CALTRWQ (4,0,SFL,SFL2)   ! PMC
  !
  ! ***  SET UP VERTICAL MIGRATION AND SETTLING BEHAVIOR
  ! ***  INITIALIZE VERTICAL VELOCTIY TO TIME DEPENDENT SETTLING VELOCITY
  !
  do K = 1,KS
    do L = 2,LA
      WWQ(L,K) = -WSFLSTT
    enddo
  enddo
  do L = 2,LA
    WWQ(L,KC) = 0.
    WWQ(L,0) = 0.
  enddo
  if( ISSFLFE >= 1 .and. ISSFLDN >= 1 )then
  !
  ! ***  DAYLIGHT CONDITIONS
  !
    if( ISDARK == 0 )then
      do K = 1,KS
        RABOVE = FLOAT(KC-K)/DZI
        do L = 2,LA
  !
  ! ***   DETERMINE DISTANCE TO SURFACE
  !
          HABOVE = RABOVE*HWQ(L)
          if( UUU(L,K) > 0. )then
  !
  ! ***    FLOOD CONDITION : SWIM UP TO MIN DIST BELOW SURFACE
  !
            if( HABOVE > DSFLMNT) WWQ(L,K) = WSFLSMT
          else
  !
  ! ***    EBB CONDITION : CONTINUE TO SINK OR SWIM UP TO MAX DIST BL SURF
  !
            if( HABOVE > DSFLMXT) WWQ(L,K) = WSFLSMT
          endif
        enddo
      enddo
    endif
  !
  ! ***  DARK CONDITIONS
  !
    if( ISDARK == 1 )then
      do K = 1,KS
        do L = 2,LA
  !
  ! ***   FLOOD CONDITION : SWIM UP TO  SURFACE
  !
          WWQ(L,K) = VVV(L,K)*WWQ(L,K)+UUU(L,K)*WSFLSMT
        enddo
      enddo
    endif
  endif
  if( SFATBTT > 0. )then
    do L = 2,LA
      WWQ(L,0) = -WSFLSTT
    enddo
  endif
  !
  ! ***  CALCULATE NET VERTICAL SWIMING OR SETTLING
  !
  if( WSFLSMT == 0. ) GOTO 100
  !
  ! ***  LIMIT VERTICAL SETTLING AND/OR SWIMMING FOR STABILITY
  !
  do K = 0,KS
    do L = 2,LA
      WWW(L,K) = MIN(WWQ(L,K),0.)
      WWW(L,K) = ABS(WWW(L,K))
      WWQ(L,K) = MAX(WWQ(L,K),0.)
    enddo
  enddo
  TMPVAL = 0.25/(DELT*DZI)
  do K = 1,KS
    do L = 2,LA
      WMAXX = TMPVAL*HWQ(L)
      WWW(L,K) = MIN(WWW(L,K),WMAXX)
      WWQ(L,K) = MIN(WWQ(L,K),WMAXX)
      WWQ(L,K) = WWQ(L,K)-WWW(L,K)
    enddo
  enddo
  do K = 1,KS
    do L = 2,LA
      FWU(L,K) = MAX(WWQ(L,K),0.)*SFL(L,K) &
          +MIN(WWQ(L,K),0.)*SFL(L,K+1)
    enddo
  enddo
  if( SFATBTT > 0. )then
    do L = 2,LA
      SFLSBOT(L) = SFLSBOT(L)-DELT*FWU(L,0)
    enddo
  endif
  do K = 1,KC
    do L = 2,LA
      SFL(L,K) = SFL(L,K) + DELT*(FWU(L,K-1)-FWU(L,K))*DZIC(L,K)/HWQ(L)
    enddo
  enddo
  GOTO 200
    100 continue
  !
  ! ***  FULLY IMPLICIT SETTLING IF SWIMMING IS ZERO EVERYWHERE
  ! ***  FULLY IMPLICIT SETTLING IN SURFACE LAYER
  !
  TMPVAL = DELT*WSFLSTT
  DZCIT = TMPVAL*DZIC(L,KC)
  do L = 2,LA
    TMPVAL1 = DZCIT/HWQ(L)
    SFL(L,KC) = SFL(L,KC)/(1.+TMPVAL1)
  enddo
  !
  ! ***  FULLY IMPLICIT SETTLING IN REMAINING LAYERS
  !
  if( KC > 1 )then
    do K = KS,1,-1
      if( K < KSZ(L) ) CYCLE
      DZCIT = TMPVAL*DZIC(L,K)
      do L = 2,LA
        TMPVAL1 = DZCIT/HWQ(L)
        SFL(L,K) = (SFL(L,K)+TMPVAL1*SFL(L,K+1))/(1.+TMPVAL1)
      enddo
    enddo
  endif
  if( SFATBTT > 0. )then
    do L = 2,LA
      SFLSBOT(L) = SFLSBOT(L)+TMPVAL*SFL(L,KSZ(L))
    enddo
  endif
    200 continue
  do L = 2,LA
    FWU(L,0) = 0.
    WWQ(L,0) = 0.
    WWW(L,0) = 0.
  enddo
  !
  ! ***  CALCULATE LINEAR DECAY
  !
  if( RKDSFLT >= 0. )then
    CDYETMP = 1./(1.+DELT*RKDSFLT)
    do K = 1,KC
      do L = 2,LA
        SFL(L,K) = CDYETMP*SFL(L,K)
      enddo
    enddo
  endif
  if( KC == 1 ) GOTO 2000
  !
  ! ***  VERTICAL DIFFUSION CALCULATION
  !
  do L = 2,LA
    HWQI(L) = 1./HWQ(L)
  enddo
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = MIN(LF+LDM-1,LA)
    do L = LF,LL
      RCDZKK = -DELT*CDZKK(L,KSZ(L))
      CCUBTMP = RCDZKK*HWQI(L)*AB(L,KSZ(L))
      CCMBTMP = 1.-CCUBTMP
      EEB = 1./CCMBTMP
      CU1(L,KSZ(L)) = CCUBTMP*EEB
      SFL(L,KSZ(L)) = SFL(L,KSZ(L))*EEB
    enddo
  enddo
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = MIN(LF+LDM-1,LA)
    do K = 2,KS
      do L = LF,LL
        if( K < KSZ(L)+1 ) CYCLE
        RCDZKMK = -DELT*CDZKMK(L,K)
        RCDZKK = -DELT*CDZKK(L,K)
        CCLBTMP = RCDZKMK*HWQI(L)*AB(L,K-1)
        CCUBTMP = RCDZKK*HWQI(L)*AB(L,K)
        CCMBTMP = 1.-CCLBTMP-CCUBTMP
        EEB = 1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
        CU1(L,K) = CCUBTMP*EEB
        SFL(L,K) = (SFL(L,K)-CCLBTMP*SFL(L,K-1))*EEB
      enddo
    enddo
  enddo
  K = KC
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = MIN(LF+LDM-1,LA)
    do L = LF,LL
      RCDZKMK = -DELT*CDZKMK(L,K)
      CCLBTMP = RCDZKMK*HWQI(L)*AB(L,K-1)
      CCMBTMP = 1.-CCLBTMP
      EEB = 1./(CCMBTMP-CCLBTMP*CU1(L,K-1))
      SFL(L,K) = (SFL(L,K)-CCLBTMP*SFL(L,K-1))*EEB
    enddo
  enddo
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = MIN(LF+LDM-1,LA)
    do K = KC-1,1,-1
      do L = LF,LL
        SFL(L,K) = SFL(L,K)-CU1(L,K)*SFL(L,K+1)
      enddo
    enddo
  enddo
  !
  ! ***  UPDATE SHELL FISH LARVAE CONCENTRATIONS
  !
    2000 continue
  do K = 1,KC
    do L = 2,LA
      SFL2(L,K) = SFL(L,K)
    enddo
  enddo
  return

END

