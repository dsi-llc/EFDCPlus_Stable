! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALEXP2T  

  ! *** *******************************************************************!
  ! ***  SUBROUTINE CALEXP2T CALCULATES EXPLICIT MOMENTUM EQUATION TERMS  
  ! ***  USING A TWO TIME LEVEL SCHEME  
  
  ! *** VARIABLES   DESCRIPTION                                 UNITS
  ! *** FUHU,FVHU   GROSS MOMENTUM, U COMPONENTS                M4/S2
  ! *** FVHV,FUHV   GROSS MOMENTUM, V COMPONENTS                M4/S2
  ! *** FX, FY      INTERNAL MODE FORCING BY LAYER              M4/S2
  ! *** FBBX, FBBY  INTERNAL MODE BOUYANCY FORCING BY LAYER     M4/S2
  ! *** FCAX,FCAY   CORIOLIS FORCING BY LAYER                   M4/S2
  ! *** DU, DV      INTERNAL SHEARS BY LAYER                    M2/S2
  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! *** 2018-01     PAUL M. CRAIG      CORRECTED WITHDRAWAL/RETURN MOMENTUM OPTION
  ! ***                                  ADDED MOMENTUM OPTION FOR STANDARD FLOW BC
  ! *** 2016-02     PAUL M. CRAIG      UPDATED SIGMA-Z (SGZ) FOR EE8.0 
  ! *** 2015-06     PAUL M. CRAIG      IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! *** 2014-09     PAUL M. CRAIG      ADDED THE LWET BYPASS APPROACH
  ! *** 2014-01     Paul M. Craig      Fixed the TOT FXVEGE/FYVEGE when partial vegetation penetration using VEGK
  ! *** 2011-03     John Hamrick/Scott James  
  ! ***                                Fixed the Vegetative Resistance from AVG to TOT using FKC
  ! *** 2011-03     Paul M. Craig      Converted to F90, added OMP
  ! *** 2010-10     Scott James        Added MHK
  ! *** 2008-12     SANG YUK/PMC       Corrected The Explicit Internal Buoyancy Forcings
  ! *** 11/07/2001  John Hamrick       Added body forces fbodyfx and fbodyfy to external momentum equations
  ! *** 11/14/2001  John Hamrick       Corrected orientation of momentum fluxes from sinks and source
  ! *** 01/02/2002  John Hamrick       Corrected 2 layer (kc = -2) curvature acceleration correction
  ! *** 01/11/2002  John Hamrick       Added ick2cor,ck2uum,ck2vvm,ck2uvm,ck2uuc,ck2vvc,ck2uvc,ck2fcx,
  ! ***                                  ck2fcy to generalize two layer momentum flux and curvature 
  ! ***                                  acceleration correction
  ! *** 01/15/2002  John Hamrick       Modified calculation of coriolis-curvature accelerations at tidal
  ! ***                                  open boundaries
  ! *** 01/23/2002  John Hamrick       Added virtual momentum sources and sinks for subgrid scale channel
  ! ***                                  interactions, including local variables tmpvec1,tmpvec2,qmcsinkx,
  ! ***                                  qmcsinky,qmcsourx,qmsoury
  ! *** 03/19/2002  John Hamrick       Added dry cell bypass and consistent initialization of dry values

  use GLOBAL  
  use FIELDS
  use CYCLONE
  use Variables_Propwash
  use Variables_Ship
  use Mod_Active_Ship
  
  implicit none
  integer :: L, LP, K, LN, LS, ID, JD, KD, NWR, IU, JU, KU, LU, NS, LNW, LSE, LL, IT, I
  integer :: LEE, LWW, LNN, LSS
  integer :: LD, NMD, LHOST, LCHNU, LW, LE, LCHNV, ND, LF
  integer,save :: NSTB

  real :: TMPANG, WUU, WVV, CACSUMT, CFEFF, VEAST2, VWEST2, FCORE, FCORW
  real :: UNORT1, USOUT1, UNORT2, USOUT2, FCORN, FCORS
  real :: UTMPATV, VTMPATU, UMAGTMP, VMAGTMP, DZPU, DZPV
  real :: TMPVAL, WVFACT, DETH, CI11H, CI12H, CI22H, DETU
  real :: CI11V, CI12V, CI21V, CI22V, CI21H, CI12U, CI21U, CI22U, DETV, CI11U
  real :: UHC, UHB, VHC, VHB, UHC1, UHB1, VHC1, VHB1, UHC2, UHB2, VHC2, VHB2
  real :: UHB1MX, UHB1MN, VHC1MX, VHC1MN, UHC1MX, UHC1MN, VHB1MX
  real :: VHB1MN, UHB2MX, UHB2MN, VHC2MX, VHC2MN, UHC2MX, UHC2MN, VHB2MX
  real :: VHB2MN, BOTT, QMF, QUMF, VEAST1, VWEST1, QWRABS, VDIR

  real,save,allocatable,dimension(:)   :: CACSUM  
  real,save,allocatable,dimension(:,:) :: DZPC
  real,save,allocatable,dimension(:)   :: TMPVEC1  
  real,save,allocatable,dimension(:)   :: TMPVEC2  
  real,save,allocatable,dimension(:,:) :: FUHJ  
  real,save,allocatable,dimension(:,:) :: FVHJ  
  real,save,allocatable,dimension(:,:) :: QMCSINKX  
  real,save,allocatable,dimension(:,:) :: QMCSINKY  
  real,save,allocatable,dimension(:,:) :: QMCSOURX  
  real,save,allocatable,dimension(:,:) :: QMCSOURY  
  real,save,allocatable,dimension(:,:) :: BK
  
  if( .not. allocated(TMPVEC1) )then
    allocate(CACSUM(NTHREADS))  
    allocate(FUHJ(LCM,KCM))  
    allocate(FVHJ(LCM,KCM))  
    allocate(TMPVEC1(KCM))  
    allocate(TMPVEC2(KCM))
    allocate(BK(KCM,NDM))
    FUHJ = 0.
    FVHJ = 0.
    TMPVEC1 = 0.
    TMPVEC2 = 0.
    NSTB = 0
    CACSUM = 0.
    BK = 0.
    if( MDCHH >= 1 )then
      allocate(QMCSINKX(LCM,KCM))  
      allocate(QMCSINKY(LCM,KCM))  
      allocate(QMCSOURX(LCM,KCM))  
      allocate(QMCSOURY(LCM,KCM))  
      QMCSINKX = 0.
      QMCSINKY = 0.
      QMCSOURX = 0.
      QMCSOURY = 0.
    endif
    if( ISPNHYDS > 0 )then
      allocate(DZPC(LCM,KCM))
      DZPC = 0.
    endif
  endif

  if( ISDYNSTP == 0 )then  
    DELT = DT  
  else  
    DELT = DTDYN  
  endif  

  if( IS2TIM == 2 )then  
    DELT = 0.5*DT  
  endif  
 
  DELTI = 1./DELT  
 
  ! *** WAVE RAMPUP FACTOR
  if( ISWAVE == 2 .or. ISWAVE == 4 )then
    if( N < NTSWV )then  
      TMPVAL = FLOAT(N)/FLOAT(NTSWV)  
      WVFACT = 0.5-0.5*COS(PI*TMPVAL)  
    else  
      WVFACT = 1.0  
    endif  
  endif
  CACSUMT = 0.

  ! *** *******************************************************************!  
  ! *** ZERO NEWLY DRY CELL SHEARS/MOMENTUM
  if( LADRY > 0 )then
    do LP = 1,LADRY
      L = LDRY(LP)  
      FCAXE(L) = 0.  
      FCAYE(L) = 0.  
      FXE(L) = 0.  
      FYE(L) = 0.  
    enddo

    do K = 1,KC
      do LP = 1,LADRY
        L = LDRY(LP)  
        FUHU(L,K) = 0.  
        FVHU(L,K) = 0.  
        FUHV(L,K) = 0.  
        FVHV(L,K) = 0.  
        FWU(L,K) = 0.
        FWV(L,K) = 0.
        CAC(L,K) = 0.0
        FCAX(L,K) = 0.
        FCAY(L,K) = 0.
        FX(L,K) = 0.
        FY(L,K) = 0.
        FBBX(L,K) = 0.
        FBBY(L,K) = 0.
        DU(L,K) = 0.0  
        DV(L,K) = 0.0 
        
        ! *** TWO LAYER ROTATIONAL EFFECTS OR WITHDRAWAL/RETURN
        FUHJ(L,K) = 0.  
        FVHJ(L,K) = 0.  
      enddo
    enddo
    
    if( ISVEG > 0 )then
      do LP = 1,LADRY
        L = LDRY(LP)  
        FXVEGE(L) = 0.
        FYVEGE(L) = 0.
      enddo
      do K = 1,KC
        do LP = 1,LADRY
          L = LDRY(LP)  
          FXVEG(L,K) = 0.0
          FYVEG(L,K) = 0.0
        enddo
      enddo
    endif
  endif

  ! *** *******************************************************************!  
  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** *******************************************************************!  
  !  
  ! ***  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS  
  !  
  !----------------------------------------------------------------------!  
  !$OMP DO PRIVATE(ND,LF,LL,LP,L)
  do ND = 1,NDM  
    LF = (ND-1)*LDMWET+1  
    LL = MIN(LF+LDMWET-1,LAWET)

    do LP = LF,LL
      L = LWET(LP)  
      FCAXE(L) = 0.  
      FCAYE(L) = 0.  
      FXE(L) = 0.  
      FYE(L) = 0.  
    enddo  
  enddo
  !$OMP END DO

  !----------------------------------------------------------------------!  
  if( IS2LMC /= 1 )then 
    ! *** STANDARD CALCULATION. WITHOUT MOMENTUM-CURVATURE CORRECTION  
    !$OMP DO PRIVATE(ND,K,L,LF,LL,LP,LE,LN,LS,LW,UHC,UHB,VHC,VHB)                
    do ND = 1,NDM 
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)  
          LS = LSC(L)
          LW = LWC(L)

          ! *** U COMPONENTS  
          UHB = 0.5*(UHDY(L,K) + UHDY(LE,K))                                  ! *** m3/s
          VHC = 0.5*(VHDX(L,K) + VHDX(LW,K))                                  ! *** m3/s 

          ! ***      |-- EAST FLOWING --|  |-- WEST FLOWING --|
          FUHU(L,K) = MAX(UHB,0.)*U(L,K)  + MIN(UHB,0.)*U(LE,K)               ! *** m4/s2 
          ! ***      |-- NORTH FLOWING -|  |-- SOUTH FLOWING -|
          FVHU(L,K) = MAX(VHC,0.)*U(LS,K) + MIN(VHC,0.)*U(L,K)                ! *** m4/s2 

          ! *** V COMPONENTS
          VHB = 0.5*(VHDX(L,K) + VHDX(LN,K))                                  ! *** m3/s 
          UHC = 0.5*(UHDY(L,K) + UHDY(LS,K))                                  ! *** m3/s 

          ! ***      |-- NORTH FLOWING -|  |-- SOUTH FLOWING -|
          FVHV(L,K) = MAX(VHB,0.)*V(L,K)  + MIN(VHB,0.)*V(LN,K)               ! *** m4/s2
          ! ***      |-- EAST FLOWING --|  |-- WEST FLOWING --|
          FUHV(L,K) = MAX(UHC,0.)*V(LW,K) + MIN(UHC,0.)*V(L,K)                ! *** m4/s2 

        enddo
      enddo  
    enddo
    !$OMP END DO

  elseif( IS2LMC == 1 .and. KC == 2 )then  
    ! *** CALCULATION FOR MOMENTUM-CURVATURE CORRECTION (TWO LAYER ONLY)
    !$OMP DO PRIVATE(ND,K,L,LP,LF,LL,LE,LN,LS,LW,UHC1,UHB1,VHC1,VHB1) &
    !$OMP    PRIVATE(UHC2,UHB2,VHC2,VHB2,UHB1MX,UHB1MN,VHC1MX)  &
    !$OMP    PRIVATE(VHC1MN,UHC1MX,UHC1MN,VHB1MX,VHB1MN,UHB2MX,UHB2MN)  &
    !$OMP    PRIVATE(VHC2MX,VHC2MN,UHC2MX,UHC2MN,VHB2MX,VHB2MN,BOTT) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)

      do LP = LF,LL
        L = LWET(LP)
        LE = LEC(L)
        LN = LNC(L)  
        LS = LSC(L)  
        LW = LWC(L)
        UHC1 = 0.5*(UHDYF(L,1)+UHDYF(LS,1))  
        UHB1 = 0.5*(UHDYF(L,1)+UHDYF(LE,1))  
        VHC1 = 0.5*(VHDXF(L,1)+VHDXF(LW,1))  
        VHB1 = 0.5*(VHDXF(L,1)+VHDXF(LN,1))  
        UHC2 = 0.5*(UHDYF(L,2)+UHDYF(LS,2))  
        UHB2 = 0.5*(UHDYF(L,2)+UHDYF(LE,2))  
        VHC2 = 0.5*(VHDXF(L,2)+VHDXF(LW,2))  
        VHB2 = 0.5*(VHDXF(L,2)+VHDXF(LN,2))  
   
        UHB1MX = 0.  
        UHB1MN = 0.  
        VHC1MX = 0.  
        VHC1MN = 0.  
        UHC1MX = 0.  
        UHC1MN = 0.  
        VHB1MX = 0.  
        VHB1MN = 0.  
        UHB2MX = 0.  
        UHB2MN = 0.  
        VHC2MX = 0.  
        VHC2MN = 0.  
        UHC2MX = 0.  
        UHC2MN = 0.  
        VHB2MX = 0.  
        VHB2MN = 0.  
   
        BOTT = ABS(UHB1*U(L,1))  
        if( BOTT > 0.0)UHB1MX = 1.+CK2UUM*(UHB2-UHB1)*(U(L,2)  -U(L,1))  /UHB1*U(L,1)  
        BOTT = ABS(UHB1*U(LE,1))  
        if( BOTT > 0.0)UHB1MN = 1.+CK2UUM*(UHB2-UHB1)*(U(LE,2)-U(LE,1))/UHB1*U(LE,1)  
        BOTT = ABS(VHC1*U(LS,1))  
        if( BOTT > 0.0)VHC1MX = 1.+CK2UVM*(VHC2-VHC1)*(U(LS,2) -U(LS,1)) /VHC1* U(LS,1)  
        BOTT = ABS(VHC1*U(L,1))  
        if( BOTT > 0.0)VHC1MN = 1.+CK2UVM*(VHC2-VHC1)*(U(L,2)  -U(L,1))  /VHC1*U(L,1)  
        BOTT = ABS(UHC1*V(LW,1))  
        if( BOTT > 0.0)UHC1MX = 1.+CK2UVM*(UHC2-UHC1)*(V(LW,2)-V(LW,1))/UHC1*V(LW,1)  
        BOTT = ABS(UHC1*V(L,1))  
        if( BOTT > 0.0)UHC1MN = 1.+CK2UVM*(UHC2-UHC1)*(V(L,2)  -V(L,1))  /UHC1*V(L,1)  
        BOTT = ABS(VHB1*V(L,1))  
        if( BOTT > 0.0)VHB1MX = 1.+CK2VVM*(VHB2-VHB1)*(V(L,2)  -V(L,1))  /VHB1*V(L,1)  
        BOTT = ABS(VHB1*V(LN,1))  
        if( BOTT > 0.0)VHB1MN = 1.+CK2VVM*(VHB2-VHB1)*(V(LN,2) -V(LN,1)) /VHB1* V(LN,1)  
        BOTT = ABS(UHB2*U(L,2))  
        if( BOTT > 0.0)UHB2MX = 1.+CK2UUM*(UHB2-UHB1)*(U(L,2)  -U(L,1))  /UHB2*U(L,2)  
        BOTT = ABS(UHB2*U(LE,2))  
        if( BOTT > 0.0)UHB2MN = 1.+CK2UUM*(UHB2-UHB1)*(U(LE,2)-U(LE,1))/UHB2*U(LE,2)  
        BOTT = ABS(VHC2*U(LS,2))  
        if( BOTT > 0.0)VHC2MX = 1.+CK2UVM*(VHC2-VHC1)*(U(LS,2) -U(LS,1)) /VHC2* U(LS,2)  
        BOTT = ABS(VHC2*U(L,2))  
        if( BOTT > 0.0)VHC2MN = 1.+CK2UVM*(VHC2-VHC1)*(U(L,2)  -U(L,1))  /VHC2*U(L,2)  
        BOTT = ABS(UHC2*V(LW,2))  
        if( BOTT > 0.0)UHC2MX = 1.+CK2UVM*(UHC2-UHC1)*(V(LW,2)-V(LW,1))/UHC2*V(LW,2)  
        BOTT = ABS(UHC2*V(L,2))  
        if( BOTT > 0.0)UHC2MN = 1.+CK2UVM*(UHC2-UHC1)*(V(L,2)  -V(L,1))  /UHC2*V(L,2)  
        BOTT = ABS(VHB2*V(L,2))  
        if( BOTT > 0.0)VHB2MX = 1.+CK2VVM*(VHB2-VHB1)*(V(L,2)  -V(L,1))  /VHB2*V(L,2)  
        BOTT = ABS(VHB2*V(LN,2))  
        if( BOTT > 0.0)VHB2MN = 1.+CK2VVM*(VHB2-VHB1)*(V(LN,2) -V(LN,1)) /VHB2* V(LN,2)  
   
        FUHU(L,1) = UHB1MX*MAX(UHB1,0.)*U(L,1)  +UHB1MN*MIN(UHB1,0.)*U(LE,1)  
        FVHU(L,1) = VHC1MX*MAX(VHC1,0.)*U(LS,1) +VHC1MN*MIN(VHC1,0.)*U(L,1)  
        FUHV(L,1) = UHC1MX*MAX(UHC1,0.)*V(LW,1)+UHC1MN*MIN(UHC1,0.)*V(L,1)  
        FVHV(L,1) = VHB1MX*MAX(VHB1,0.)*V(L,1)  +VHB1MN*MIN(VHB1,0.)*V(LN,1)  
        FUHJ(L,1) = 0.  
        FVHJ(L,1) = 0.  
        FUHU(L,2) = UHB2MX*MAX(UHB2,0.)*U(L,2)  +UHB2MN*MIN(UHB2,0.)*U(LE,2)  
        FVHU(L,2) = VHC2MX*MAX(VHC2,0.)*U(LS,2) +VHC2MN*MIN(VHC2,0.)*U(L,2)  
        FUHV(L,2) = UHC2MX*MAX(UHC2,0.)*V(LW,2)+UHC2MN*MIN(UHC2,0.)*V(L,2)  
        FVHV(L,2) = VHB2MX*MAX(VHB2,0.)*V(L,2)  +VHB2MN*MIN(VHB2,0.)*V(LN,2)  
        FUHJ(L,2) = 0.  
        FVHJ(L,2) = 0.  
      enddo
    enddo  ! *** END OF DOMAIN
    !$OMP END DO

  endif
  
  !$OMP SINGLE
  ! *** ADD WITHDRAWAL/RETURN FLOW MOMENTUM FLUXES
  do NWR = 1,NQWR  
    if( ABS(WITH_RET(NWR).NQWRMFU) > 0 )then  
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS = WITH_RET(NWR).NQWRSERQ  
      if( QWRSERT(NS) >= 0. )then
        ! *** Original Withdrawal/Return
        IU = WITH_RET(NWR).IQWRU  
        JU = WITH_RET(NWR).JQWRU  
        KU = WITH_RET(NWR).KQWRU
        VDIR = 1.
      else
        ! *** Reverse Flow Withdrawal/Return
        IU = WITH_RET(NWR).IQWRD  
        JU = WITH_RET(NWR).JQWRD  
        KU = WITH_RET(NWR).KQWRD 
        VDIR = -1.
      endif
      LU = LIJ(IU,JU)  

      QWRABS = ABS(QWRSERT(NS))
      QMF = WITH_RET(NWR).QWR + QWRABS 
      QUMF = VDIR*QMF*( QMF/( HPK(LU,KU)*WITH_RET(NWR).BQWRMFU ) )   ! *** M4/S2
      if( WITH_RET(NWR).NQWRMFU ==  1 ) FUHJ(LU     ,KU) = -QUMF  
      if( WITH_RET(NWR).NQWRMFU ==  2 ) FVHJ(LU     ,KU) = -QUMF  
      if( WITH_RET(NWR).NQWRMFU ==  3 ) FUHJ(LEC(LU),KU) = -QUMF  
      if( WITH_RET(NWR).NQWRMFU ==  4 ) FVHJ(LNC(LU),KU) = -QUMF  
      if( WITH_RET(NWR).NQWRMFU == -1 ) FUHJ(LU     ,KU) = QUMF  
      if( WITH_RET(NWR).NQWRMFU == -2 ) FVHJ(LU     ,KU) = QUMF  
      if( WITH_RET(NWR).NQWRMFU == -3 ) FUHJ(LEC(LU),KU) = QUMF  
      if( WITH_RET(NWR).NQWRMFU == -4 ) FVHJ(LNC(LU),KU) = QUMF  
    endif  
    if( ABS(WITH_RET(NWR).NQWRMFD) > 0 )then  
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS = WITH_RET(NWR).NQWRSERQ  
      if( QWRSERT(NS) >= 0. )then
        ! *** Original Withdrawal/Return
        ID = WITH_RET(NWR).IQWRD  
        JD = WITH_RET(NWR).JQWRD  
        KD = WITH_RET(NWR).KQWRD  
        VDIR = 1.
      else
        ! *** Reverse Flow Withdrawal/Return
        ID = WITH_RET(NWR).IQWRU  
        JD = WITH_RET(NWR).JQWRU  
        KD = WITH_RET(NWR).KQWRU 
        VDIR = -1.
      endif
      LD = LIJ(ID,JD)

      QWRABS = ABS(QWRSERT(NS))
      QMF = WITH_RET(NWR).QWR + QWRABS 
      QUMF = VDIR*QMF*( QMF/( HPK(LD,KD)*WITH_RET(NWR).BQWRMFD ) )   ! *** M4/S2
      if( WITH_RET(NWR).NQWRMFD ==  1 ) FUHJ(LD     ,KD) = -QUMF  
      if( WITH_RET(NWR).NQWRMFD ==  2 ) FVHJ(LD     ,KD) = -QUMF  
      if( WITH_RET(NWR).NQWRMFD ==  3 ) FUHJ(LEC(LD),KD) = -QUMF  
      if( WITH_RET(NWR).NQWRMFD ==  4 ) FVHJ(LNC(LD),KD) = -QUMF  
      if( WITH_RET(NWR).NQWRMFD == -1 ) FUHJ(LD     ,KD) = QUMF  
      if( WITH_RET(NWR).NQWRMFD == -2 ) FVHJ(LD     ,KD) = QUMF  
      if( WITH_RET(NWR).NQWRMFD == -3 ) FUHJ(LEC(LD),KD) = QUMF  
      if( WITH_RET(NWR).NQWRMFD == -4 ) FVHJ(LNC(LD),KD) = QUMF  
    endif  
  enddo  
  
  ! *** ADD QSER MOMENTUM FLUXES

  do LL = 1,NQSIJ
    if( ABS(BCFL(LL).NQSMF) > 0 .and. BCFL(LL).NQSMF /= 5 )then  
      L = BCFL(LL).L
      do K = KSZ(L),KC
        ! *** Handle reversing flows in/out of domain
        if( QSERCELL(K,LL) >= 0. )then
          VDIR = 1.
        else
          VDIR = -1.
        endif
        !LD = LIJ(BCFL(LL).I,BCFL(LL).J)
        
        QMF = ABS(QSERCELL(K,LL))
        QUMF = VDIR*QMF*( QMF/( HPK(L,K)*BCFL(LL).QWIDTH ) )   ! *** M4/S2
        if( BCFL(LL).NQSMF ==  1 ) FUHJ(L     ,K) = -QUMF  
        if( BCFL(LL).NQSMF ==  2 ) FVHJ(L     ,K) = -QUMF  
        if( BCFL(LL).NQSMF ==  3 ) FUHJ(LEC(L),K) = -QUMF  
        if( BCFL(LL).NQSMF ==  4 ) FVHJ(LNC(L),K) = -QUMF  
        if( BCFL(LL).NQSMF == -1 ) FUHJ(L     ,K) = QUMF  
        if( BCFL(LL).NQSMF == -2 ) FVHJ(L     ,K) = QUMF  
        if( BCFL(LL).NQSMF == -3 ) FUHJ(LEC(L),K) = QUMF  
        if( BCFL(LL).NQSMF == -4 ) FVHJ(LNC(L),K) = QUMF  
      enddo
    endif  
  enddo

  ! *** Add propeller efflux velocity momentum
  if( ISPROPWASH == 2 .and. NACTIVESHIPS > 0 )then
    call add_ship_momentum(FUHJ, FVHJ)
  endif
  !$OMP END SINGLE

  !----------------------------------------------------------------------!  
  !  
  ! *** COMPUTE VERTICAL ACCELERATIONS
  if( KC > 1 )then
    !$OMP DO PRIVATE(ND,K,LP,L,LW,LS,WUU,WVV)
    do ND = 1,NDM  
      do K = 1,KS
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LW = LWC(L)
          LS = LSC(L)
          WUU = 0.5*DXYU(L)*(W(L,K) + W(LW,K))  
          WVV = 0.5*DXYV(L)*(W(L,K) + W(LS,K))  
          FWU(L,K) = MAX(WUU,0.)*U(L,K) + MIN(WUU,0.)*U(L,K+1)
          FWV(L,K) = MAX(WVV,0.)*V(L,K) + MIN(WVV,0.)*V(L,K+1)
        enddo  
      enddo
    enddo  
    !$OMP END DO
  
    ! *** APPLY OPEN BOUNDARYS
    !$OMP SINGLE
    do LL = 1,NBCSOP2
      L = LOBCS2(LL)
      do K = 1,KS  
        FWU(L,K) = 0.0
        FWV(L,K) = 0.0
      enddo  
    enddo 
    !$OMP END SINGLE 
  endif
  
  ! **********************************************************************!  
  ! *** BLOCK MOMENTUM FLUX ON LAND SIDE OF TRIANGULAR CELLS  
  if( ITRICELL > 0 )then
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          FUHU(L,K) = STCUV(L)*FUHU(L,K)
          FVHV(L,K) = STCUV(L)*FVHV(L,K)
        enddo
      enddo
    enddo   ! ***  END OF DOMAIN
    !$OMP END DO
  endif
  
  ! *** *******************************************************************!  
  ! *** CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS  
  !$OMP SINGLE
  CACSUM = 0. 
  CFMAX = CF  
  !$OMP END SINGLE
  
  if( ISCURVATURE )then

    if( ISDCCA == 0 )then  
      ! *** STANDARD CALCULATIONS, NO DIAGNOSTICS
      !$OMP DO PRIVATE(ND,K,LP,L,LE,LN)
      do ND = 1,NDM  

        do K = 1,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            LE = LEC(L)
            LN = LNC(L)
            CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)
            CACSUM(ND) = CACSUM(ND)+CAC(L,K)
          enddo
        enddo
      enddo   ! ***  END OF DOMAIN
      !$OMP END DO
      
    else  
      ! *** STANDARD CALCULATIONS, WITH DIAGNOSTICS
      !$OMP SINGLE
      IT = 1
      do K = 1,KC  
        do L = 2,LA
          LE = LEC(L)
          LN = LNC(L)  
          CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)  
          CFEFF = ABS(CAC(L,K))*DXYIP(L)*HPI(L)  
          CFMAX = MAX(CFMAX,CFEFF)  
          CACSUM(IT) = CACSUM(IT)+CAC(L,K) 
        enddo  
      enddo  
      !$OMP END SINGLE
    endif  

    ! *** ENSURE FCAY & FCAX ARE RESET
    !$OMP SINGLE
    CACSUMT = ABS(SUM(CACSUM(:)))
    !$OMP END SINGLE
    
    if( CACSUMT < 1.E-7 )then
      !$OMP DO PRIVATE(ND,K,LP,L)
      do ND = 1,NDM  
        do K = 1,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            FCAX(L,K) = 0.
            FCAY(L,K) = 0.
          enddo
        enddo
      enddo  
      !$OMP END DO
    endif
        
  endif

  ! *** *******************************************************************!  
  !  
  ! ***  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS  
  !  
  !----------------------------------------------------------------------!  

  ! ***  STANDARD CALCULATION
  if( CACSUMT > 1.E-7 )then
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW,LNW,LSE)
    do ND = 1,NDM
      do K = 1,KC  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)  
          LS = LSC(L)  
          LW = LWC(L)
          LNW = LNWC(L)  
          LSE = LSEC(L)  
          FCAX(L,K) = 0.25*SCAX(L)*( CAC(L,K)*(V(LN,K)+V(L,K)) + CAC(LW,K)*(V(LNW,K)+V(LW,K)) )  
          FCAY(L,K) = 0.25*SCAY(L)*( CAC(L,K)*(U(LE,K)+U(L,K)) + CAC(LS,K)*(U(LSE,K)+U(LS,K)) )  
        enddo
      enddo  
    enddo  
    !$OMP END DO

  endif

  !----------------------------------------------------------------------!  
  !  
  ! *** CALCULATION FOR MOMENTUM-CURVATURE CORRECTION  
  if( IS2LMC == 1 .and. CACSUMT > 1.E-7 )then  
    !$OMP SINGLE
    do LP = 1,LAWET
      L = LWET(LP)
      LE = LEC(L)
      LN = LNC(L)  
      LS = LSC(L)  
      LW = LWC(L)
      LNW = LNWC(L)  
      LSE = LSEC(L)  

      VEAST1 = V(LN,1)+V(L,1)  
      VWEST1 = V(LNW,1)+V(LW,1)  
      VEAST2 = V(LN,2)+V(L,2)  
      VWEST2 = V(LNW,2)+V(LW,2)  
        
      FCORE = CK2FCX*(CAC(L,2) -CAC(L,1)) *(VEAST2-VEAST1)  
      FCORW = CK2FCX*(CAC(LW,2)-CAC(LW,1))*(VWEST2-VWEST1)
        
      FCAX(L,1) = 0.25*SCAX(L)*(CAC(L,1)*VEAST1+FCORE +CAC(LW,1)     *VWEST1+FCORW)  
      FCAX(L,2) = 0.25*SCAX(L)*(CAC(L,2)*VEAST2+FCORE +CAC(LWC(LW),2)*VWEST2+FCORW)

      UNORT1 = U(LE,1)+U(L,1)  
      USOUT1 = U(LSE,1)+U(LS,1)  
      UNORT2 = U(LE,2)+U(L,2)  
      USOUT2 = U(LSE,2)+U(LS,2)  
      FCORN = CK2FCY*(CAC(L ,2)-CAC(L ,1))*(UNORT2-UNORT1)  
      FCORS = CK2FCY*(CAC(LS,2)-CAC(LS,1))*(USOUT2-USOUT1)  

      FCAY(L,1) = 0.25*SCAY(L)*(CAC(L,1)*UNORT1+FCORN +CAC(LS,1)*USOUT1+FCORS)  
      FCAY(L,2) = 0.25*SCAY(L)*(CAC(L,2)*UNORT2+FCORN +CAC(LS,2)*USOUT2+FCORS)  
    enddo  
    !$OMP END SINGLE
  endif  
  
  !----------------------------------------------------------------------!  
  ! ***  MODIFICATION FOR TYPE 2 OPEN BOUNDARIES  
  !$OMP SINGLE
  do LL = 1,NPBW  
    if( ISPBW(LL) == 2 .or. ISPBW(LL) == 5 )then  
      L = LPBW(LL)+1  
      LN = LNC(L)  
      do K = KSZ(L),KC  
        FCAX(L,K) = 0.5*SCAX(L)*CAC(L,K)*(V(LN,K)+V(L,K))  
      enddo  
    endif  
  enddo  

  do LL = 1,NPBE  
    if( ISPBE(LL) == 2 .or. ISPBE(LL) == 5 )then  
      L = LPBE(LL)  
      LNW = LNWC(L)  
      do K = KSZ(L),KC  
        FCAX(L,K) = 0.5*SCAX(L)*CAC(LWC(L),K)*(V(LNW,K)+V(LWC(L),K))  
      enddo  
    endif  
  enddo  

  do LL = 1,NPBS  
    if( ISPBS(LL) == 2 .or. ISPBS(LL) == 5 )then  
      L = LNC(LPBS(LL))  
      do K = KSZ(L),KC  
        FCAY(L,K) = 0.5*SCAY(L)*CAC(L,K)*(U(LEC(L),K)+U(L,K))  
      enddo  
    endif  
  enddo  

  do LL = 1,NPBN  
    if( ISPBN(LL) == 2 .or. ISPBN(LL) == 5 )then  
      L = LPBN(LL)  
      LS = LSC(L)  
      LSE = LSEC(L)  
      do K = KSZ(L),KC  
        FCAY(L,K) = 0.5*SCAY(L)*CAC(LS,K)*(U(LSE,K)+U(LS,K))  
      enddo  
    endif  
  enddo  
  !$OMP END SINGLE

  !----------------------------------------------------------------------!  
  !  INITIALIZE EXTERNAL FORCINGS WITH GROSS MOMENTUM
  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW)
  do ND = 1,NDM
    do K = 1,KC
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        LE = LEC(L)
        LN = LNC(L) 
        LS = LSC(L) 
        LW = LWC(L)
        FX(L,K) = FSGZU(L,K)*( FUHU(L,K)-FUHU(LW,K) + FVHU(LN,K)-FVHU(L,K) ) + FUHJ(L,K)   ! ***  M4/S2
        FY(L,K) = FSGZV(L,K)*( FVHV(L,K)-FVHV(LS,K) + FUHV(LE,K)-FUHV(L,K) ) + FVHJ(L,K)   ! ***  M4/S2
      enddo  
    enddo
  enddo
  !$OMP END DO

  ! *** TREAT BC'S NEAR EDGES
  !$OMP SINGLE
  do LL = 1,NBCS
    ! *** BC CELL
    L = LBCS(LL)
    do K = KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K) + FUHJ(L,K)   ! ***  M4/S2
      FY(L,K) = SAAY(L)*FY(L,K) + FVHJ(L,K)   ! ***  M4/S2
    enddo

    ! *** EAST/WEST ADJACENT CELL
    L = MAX(1,LBERC(LL))
    do K = KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K) + FUHJ(L,K)   ! ***  M4/S2
    enddo

    ! *** NORTH/SOUTH ADJACENT CELL
    L = MAX(1,LBNRC(LL))
    do K = KSZ(L),KC
      FY(L,K) = SAAY(L)*FY(L,K) + FVHJ(L,K)   ! ***  M4/S2
    enddo
    
  enddo
  
  ! *** Selective zero
  if( ISPROPWASH == 2 .and. NACTIVESHIPS > 0 )then
    do I = 1, total_ships
      if( all_ships(i).efflux_vel > 0.0 )then
        L = all_ships(i).pos.cell
        FUHJ(L,:)      = 0.0
        FUHJ(LEC(L),:) = 0.0
        FVHJ(L,:)      = 0.0
        FVHJ(LNC(L),:) = 0.0
      endif
      if( all_ships(i).ship.num_fixed_cells > 0 )then
        do NS = 1, all_ships(i).ship.num_fixed_cells
          L = all_ships(i).ship.fixed_cells(NS)
          FUHJ(L,:)      = 0.0
          FUHJ(LEC(L),:) = 0.0
          FVHJ(L,:)      = 0.0
          FVHJ(LNC(L),:) = 0.0
        enddo
      endif
    enddo
  endif
  
  !$OMP END SINGLE
      
  ! *** *******************************************************************!  
  !  
  ! ***  ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS  
  !  
  !----------------------------------------------------------------------!  
  if( ISVEG > 0 )then  

    ! *** ADD IN MHK DEVICES, IF NEEDED
    !$OMP SINGLE
    if( LMHK ) CALL MHKPWRDIS
    !$OMP END SINGLE

    if( ISDRY > 0 .and. LADRY > 0 )then
      ! *** UPDATE ACTIVE CELL BY LAYER LIST
      !$OMP DO PRIVATE(ND,K,LN,LP,L)
      do ND = 1,NDM  
        do K = 1,KC  
          LN = 0
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            if( LVEG(L) )then
              LN = LN+1
              LKVEG(LN,K,ND) = L
            endif
          enddo
          LLVEG(K,ND) = LN
        enddo
      enddo
      !$OMP END DO
    endif
    
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LW,LE,LS,LN,LNW,LSE) &
    !$OMP    PRIVATE(UTMPATV,VTMPATU,UMAGTMP,VMAGTMP)
    do ND = 1,NDM  
      do LP = 1,LLVEG(KC,ND)
        L = LKVEG(LP,KC,ND) 
        FXVEGE(L) = 0.  
        FYVEGE(L) = 0.  
      enddo  
      
      do K = 1,KC
        do LP = 1,LLVEG(K,ND)
          L = LKVEG(LP,K,ND)
          LW = LWC(L)    !west cell
          LE = LEC(L)    !east cell
          LS = LSC(L)    !south cell
          LN = LNC(L)    !north cell
          LNW = LNWC(L)  !northwest cell
          LSE = LSEC(L)  !southeast cell
          UTMPATV = 0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))       !u-velocity at v face
          VTMPATU = 0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))       !v-velocity at u face
          UMAGTMP = SQRT( U(L,K)*U(L,K) +   VTMPATU*VTMPATU )    !u-face velocity vector magnitude
          VMAGTMP = SQRT( UTMPATV*UTMPATV + V(L,K)*V(L,K) )      !v-face velocity vector magnitude

          !FXVEG/FYVEG come from CALTBXY unitless, but they are really just a form of the drag coefficient with terms accounting for the area density
          !FXVEG/FYVEG only change inasmuch as the water depth changes and are zero in layers not penetrated by vegetation
          !FXVEG/FYVEG are C_d(N/L^2) 
          !FXVEG/FYVEG are now multiplied by the cell area and cell-averaged velocity
          !FXVEG/FYVEG are C_d(N/L^2)A|q|
          FXVEG(L,K) = UMAGTMP*SUB3D(L,K)*FXVEG(L,K)   ![m/s] q_xC_d
          FYVEG(L,K) = VMAGTMP*SVB3D(L,K)*FYVEG(L,K)   ![m/s] q_yC_d
            
          !FXVEG/FXVEGE are multiplied by the local velocity to yield units of [m^4/s^2]
          !FXVEG/FXVEGE are added to the body forces as C_d(N/L^2)A|q|q
          FXVEGE(L) = FXVEGE(L)+FXVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FUHDXE 
          FYVEGE(L) = FYVEGE(L)+FYVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FVHDYE
        enddo  
      enddo  

      ! *** ADD VEGETATIVE DRAG TO INTERNAL MODE SHEAR
      if(ISVEG == 1 )then
        do K = 1,KC
          do LP = 1,LLVEG(K,ND)
            L = LKVEG(LP,K,ND)
            if( (K-KSZ(L)+1) > INT(VEGK(L)+1.) ) CYCLE
            ! ***               (        m/s         )  m/s    m2 
            FX(L,K) = FX(L,K) + (FXVEG(L,K)-FXVEGE(L))*U(L,K)*DXYU(L) ! *** (m^4/s^2)
            FY(L,K) = FY(L,K) + (FYVEG(L,K)-FYVEGE(L))*V(L,K)*DXYV(L) ! *** (m^4/s^2)
          enddo  
        enddo  
      else
        do K = 1,KC
          do LP = 1,LLVEG(K,ND)
            L = LKVEG(LP,K,ND)
            if( (K-KSZ(L)+1) > INT(VEGK(L)+1.) ) CYCLE
            ! ***                m3/s       m/s   
            FX(L,K) = FX(L,K) + FXVEG(L,K)*U(L,K) ! *** (m^4/s^2)
            FY(L,K) = FY(L,K) + FYVEG(L,K)*V(L,K) ! *** (m^4/s^2)
          enddo  
        enddo
      endif


      ! *** CONVERT THE AVG FXVEGE/FYVEGE TO TOTAL FXVEGE/FYVEGE
      ! *** Calculate vegetative dissipation for FUHDYE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
      do LP = 1,LLVEG(KC,ND)
        L = LKVEG(LP,KC,ND) 
        FXVEGE(L) = FXVEGE(L)*HUI(L)*VEGK(L)   ! *** Calculate vegetative dissipation for FUHDYE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
        FYVEGE(L) = FYVEGE(L)*HVI(L)*VEGK(L)   ! *** Calculate vegetative dissipation for FVHDXE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
      enddo 
      
      if( LMHK )then 
        do LP = 1,LLVEG(KC,ND)
          L = LKVEG(LP,KC,ND) 
          FXVEGE(L) = FXVEGE(L) + FXMHKE(L)   ! Add MHK to vegetative dissipation in FUHDYE for momentum conservation in CALPUV multiply by HUI took place in MHKPWRDIS
          FYVEGE(L) = FYVEGE(L) + FYMHKE(L)   ! Add MHK to vegetative dissipation in FVHDXE for momentum conservation in CALPUV multiply by HVI took place in MHKPWRDIS
          FXVEGE(L) = FXVEGE(L) + FXSUPE(L)   ! Add MHK support to vegetative dissipation in FUHDYE for momentum conservation in CALPUV
          FYVEGE(L) = FYVEGE(L) + FYSUPE(L)   ! Add MHK support to vegetative dissipation in FVHDXE for momentum conservation in CALPUV
        enddo
      endif
    enddo  ! *** END OF DOMAIN
    !$OMP END DO

  endif  

  ! *** *******************************************************************!  
  !  
  ! ***  ADD HORIZONTAL MOMENTUM DIFFUSION TO ADVECTIVE ACCELERATIONS  
  !  
  !----------------------------------------------------------------------!  
  if( ISHDMF >= 1 )then

    ! *** Modification for open boundaries for large cell aspect ratios
    !$OMP SINGLE
    do LL = 1,NPBW
      if( ISPBW(LL) == 0 .or. ISPBW(LL) == 1 .or. ISPBW(LL) == 4 )then
        L   = LEC(LPBW(LL))
        LE  = LEC(L)
        LEE = LEC(LE)
        ! *** Both adjacent and connected cell must both meet aspect criteria
        if( MAX(DXP(LE),DYP(LE)) / MIN(DXP(LE),DYP(LE)) > 4.0 )then
          if( MAX(DXP(LEE),DYP(LEE)) / MIN(DXP(LEE),DYP(LEE)) > 4.0 )then
            do K = KSZ(LE),KC
              FMDUY(LE,K) = 0.0
            enddo
            do K = KSZ(LEE),KC
              FMDUY(LEE,K) = 0.0
            enddo
          endif
        endif
      endif
    enddo
  
    do LL = 1,NPBE
      if( ISPBE(LL) == 0 .or. ISPBE(LL) == 1 .or. ISPBE(LL) == 4 )then
        L   = LPBE(LL)
        LW  = LWC(L)
        LWW = LEC(LW)
        ! *** Both adjacent and connected cell must both meet aspect criteria
        if( MAX(DXP(LW),DYP(LW)) / MIN(DXP(LW),DYP(LW)) > 4.0 )then
          if( MAX(DXP(LWW),DYP(LWW)) / MIN(DXP(LWW),DYP(LWW)) > 4.0 )then
            do K = KSZ(LW),KC
              FMDUY(LW,K) = 0.0
            enddo
            do K = KSZ(LWW),KC
              FMDUY(LWW,K) = 0.0
            enddo
          endif
        endif
      endif
    enddo
  
    do LL = 1,NPBS
      if( ISPBS(LL) == 0 .or. ISPBS(LL) == 1 .or. ISPBS(LL) == 4 )then
        L = LPBS(LL)
        LN = LNC(L)
        LNN = LEC(LN)
        ! *** Both adjacent and connected cell must both meet aspect criteria
        if( MAX(DXP(LN),DYP(LN)) / MIN(DXP(LN),DYP(LN)) > 4.0 )then
          if( MAX(DXP(LNN),DYP(LNN)) / MIN(DXP(LNN),DYP(LNN)) > 4.0 )then
            do K = KSZ(LN),KC
              FMDUY(LN,K) = 0.0
            enddo
            do K = KSZ(LNN),KC
              FMDUY(LNN,K) = 0.0
            enddo
          endif
        endif
      endif
    enddo
  
    do LL = 1,NPBN
      if( ISPBN(LL) == 0 .or. ISPBN(LL) == 1 .or. ISPBN(LL) == 4 )then
        L = LPBN(LL)
        LS = LSC(L)
        LSS = LEC(LS)
        ! *** Both adjacent and connected cell must both meet aspect criteria
        if( MAX(DXP(LS),DYP(LS)) / MIN(DXP(LS),DYP(LS)) > 4.0 )then
          if( MAX(DXP(LSS),DYP(LSS)) / MIN(DXP(LSS),DYP(LSS)) > 4.0 )then
            do K = KSZ(LS),KC
              FMDUY(LS,K) = 0.0
            enddo
            do K = KSZ(LSS),KC
              FMDUY(LSS,K) = 0.0
            enddo
          endif
        endif
      endif
    enddo
    !$OMP END SINGLE
    
    !$OMP DO PRIVATE(ND,K,LP,L) 
    do ND = 1,NDM  
      do K = 1,KC  
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          FX(L,K) = FX(L,K) - SDX(L)*( FMDUX(L,K) + FMDUY(L,K) )
          FY(L,K) = FY(L,K) - SDY(L)*( FMDVX(L,K) + FMDVY(L,K) )  
        enddo  
      enddo
    enddo
    !$OMP END DO

  endif  

  ! *** *******************************************************************!  
  !  
  ! *** ADD BODY FORCE TO ADVECTIVE ACCELERATIONS
  ! *** DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF = 1
  ! *** DISTRIBUTE OVER SURFACE LAYER IF ISBODYF = 2
  !  
  !----------------------------------------------------------------------!  
  if( ISBODYF == 1 )then  
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM  
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          FX(L,K) = FX(L,K) - DYU(L)*HU(L)*FBODYFX(L,K)
          FY(L,K) = FY(L,K) - DXV(L)*HV(L)*FBODYFY(L,K)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif  

  if( ISBODYF == 2 )then  
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      do LP = LF,LL
        L = LWET(LP)
        FX(L,KC) = FX(L,KC) - DZIC(L,KC)*DYU(L)*HU(L)*FBODYFX(L,KC)
        FY(L,KC) = FY(L,KC) - DZIC(L,KC)*DXV(L)*HV(L)*FBODYFY(L,KC)
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif

  ! *** *******************************************************************!  
  ! *** ADD EXPLICIT NONHYDROSTATIC PRESSURE   
  if( KC > 1 .and. ISPNHYDS >= 1 )then  
    !$OMP DO PRIVATE(ND,L,LP,K,TMPVAL) 
    do ND = 1,NDM
      do LP = 1,LLWET(KC,ND)
        L = LKWET(LP,KC,ND)  
        TMPVAL = 2./(DZC(L,KSZ(L))+DZC(L,KSZ(L)+1))  
        DZPC(L,1) = TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))  
      enddo
    enddo  
    !$OMP END DO

    !$OMP DO PRIVATE(ND,L,LP,TMPVAL) 
    do ND = 1,NDM  
      do LP = 1,LLWET(KC,ND)
        L = LKWET(LP,KC,ND)  
        TMPVAL = 2./(DZC(L,KC)+DZC(L,KC-1))  
        DZPC(L,KC) = TMPVAL*(PNHYDS(L,KC)-PNHYDS(L,KC-1))  
      enddo
    enddo  
    !$OMP END DO

    if( KC >= 3 )then  
      !$OMP DO PRIVATE(ND,K,LP,L,TMPVAL) 
      do ND = 1,NDM  
        do K = 2,KS  
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            TMPVAL = 2./(DZC(L,K+1)+2.*DZC(L,K)+DZC(L,K-1) )  
            DZPC(L,K) = TMPVAL*(PNHYDS(L,K+1)-PNHYDS(L,K-1))  
          enddo  
        enddo  
      enddo
      !$OMP END DO
    endif  

   !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,DZPU,DZPV)
    do ND = 1,NDM  
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LW = LWC(L)
          LS = LSC(L)  
          DZPU = 0.5*(DZPC(L,K)+DZPC(LW,K))  
          DZPV = 0.5*(DZPC(L,K)+DZPC(LS,K))  
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*( HU(L)*(PNHYDS(L,K)-PNHYDS(LW,K) ) - ( BELV(L)-BELV(LW) + ZZ(L ,K)*(HP(L)-HP(LW)) )*DZPU )
          FY(L,K) = FY(L,K) + SVB3D(L,K)*DXV(L)*( HV(L)*(PNHYDS(L,K)-PNHYDS(LS,K) ) - ( BELV(L)-BELV(LS) + ZZ(LS,K)*(HP(L)-HP(LS)) )*DZPV )
        enddo
      enddo  
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif

  if( IATMP > 0 .or. PRESS.IFLAG > 0 .or. ICYCLONE > 0 )then
    ! *** ADD AIR PRESSURE GRADIENT (PRESS.IFLAG ASSUMES A SPATIALLY VARYING PRESSURE FIELD)
    !$OMP DO PRIVATE(ND,K,LP,L,LS,LW)
    do ND = 1,NDM  
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LW = LWC(L)
          LS = LSC(L)
          ! ***                            M       M         M2/S2        
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*HU(L)*( ATMP(L)-ATMP(LW) )
          FY(L,K) = FY(L,K) + SVB3D(L,K)*DXV(L)*HV(L)*( ATMP(L)-ATMP(LS) )
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif

  !----------------------------------------------------------------------!  
  ! ***  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.  
  if( ISWAVE == 2 .or. ISWAVE == 4 )then

    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L) 
    do ND = 1,NDM  
      LF = (ND-1)*NWVCELLS+1  
      LL = MIN(LF+NWVCELLS-1,NWVCELLS)
      
      do K = 1,KC
        do LP = LF,LL
          L = LWVCELL(LP)  
          if( LKSZ(L,K) ) CYCLE  
          FX(L,K) = FX(L,K)+SUB3D(L,K)*WVFACT*SAAX(L)*FXWAVE(L,K)
          FY(L,K) = FY(L,K)+SVB3D(L,K)*WVFACT*SAAY(L)*FYWAVE(L,K)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif  
  
  
  ! *** Filter cells which have been impacted by newly wet cells
  !if( ISDRY > 0 )then
  !  !$OMP DO PRIVATE(ND,K,LP,L)
  !  do ND = 1,NDM   
  !    do K = 1,KC                delme - uncomment this block to match 3TL filtering
  !      do LP = 1,LLWET(K,ND)
  !        L = LKWET(LP,K,ND)
  !        if( NWET(L) > 10 ) CYCLE
  !        
  !        FX(L,K) = FX1(L,K) + SIGN(MAX(1E-8,MIN( ABS(FX(L,K)-FX1(L,K)), ABS(2.0*FX1(L,K)) )), FX(L,K)-FX1(L,K) )
  !        FY(L,K) = FY1(L,K) + SIGN(MAX(1E-8,MIN( ABS(FY(L,K)-FY1(L,K)), ABS(2.0*FY1(L,K)) )), FY(L,K)-FY1(L,K) )
  !        FX1(L,K) = FX(L,K)
  !        FY1(L,K) = FY(L,K)
  !      enddo
  !    enddo
  !  enddo   ! *** END OF DOMAIN
  !  !$OMP END DO
  !endif
  
  ! *** *******************************************************************!  
  ! ***  CALCULATE TOTAL EXTERNAL ACCELERATIONS  
  !$OMP DO PRIVATE(ND,L,K,LP)
  do ND = 1,NDM  
    do K = 1,KC
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        FCAXE(L) = FCAXE(L) + FCAX(L,K)*SGZU(L,K)  
        FCAYE(L) = FCAYE(L) + FCAY(L,K)*SGZV(L,K)  
        FXE(L) = FXE(L) + FX(L,K)*SGZU(L,K)  
        FYE(L) = FYE(L) + FY(L,K)*SGZV(L,K)  
      enddo  
    enddo
  enddo   ! *** END OF DOMAIN
  !$OMP END DO

  ! *** *******************************************************************!  
  ! ***  COMPLETE CALCULATION OF INTERNAL MODE ADVECTIVE ACCELERATIONS  
  if( KC > 1 )then

    ! *** *******************************************************************!  
    ! *** ADD VERTICAL MOMENTUM FLUX INTO HORIZONTAL MOMENTUM FLUX
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
    do ND = 1,NDM  
      do K = 1,KC
        do LP = 1,LLWETZ(K,ND)
          L = LKWETZ(LP,K,ND)
          FX(L,K) = FX(L,K) + SAAX(L)*( FWU(L,K)-FWU(L,K-1) )*DZIC(L,K)
          FY(L,K) = FY(L,K) + SAAY(L)*( FWV(L,K)-FWV(L,K-1) )*DZIC(L,K)
        enddo
      enddo  
    enddo   ! *** END OF DOMAIN
    !$OMP END DO

    ! *** *******************************************************************!  
    ! ***  ADD SUBGRID SCALE CHANNEL VIRTURAL MOMENTUM SOURCES AND SINKS  
    if( MDCHH >= 1 .and. ISCHAN == 3 )then  
      !$OMP DO PRIVATE(ND,LF,LL,L,LP,K)
      do ND = 1,NDM  
        LF = (ND-1)*LDMWET+1  
        LL = MIN(LF+LDMWET-1,LAWET)

        do K = 1,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            QMCSOURX(L,K) = 0.  
            QMCSOURY(L,K) = 0.  
            QMCSINKX(L,K) = 0.  
            QMCSINKY(L,K) = 0.  
          enddo  
        enddo  
      enddo  
      !$OMP END DO

      !$OMP SINGLE
      do NMD = 1,MDCHH  

        LHOST = LMDCHH(NMD)  
        LCHNU = LMDCHU(NMD)  
        LCHNV = LMDCHV(NMD)  

        DETH = CUE(LHOST)*CVN(LHOST)-CUN(LHOST)*CVE(LHOST)  
        CI11H = CVN(LHOST)/DETH  
        CI12H = -CUN(LHOST)/DETH  
        CI21H = -CVE(LHOST)/DETH  
        CI22H = CUE(LHOST)/DETH  

        DETU = CUE(LCHNU)*CVN(LCHNU)-CUN(LCHNU)*CVE(LCHNU)  
        CI11U = CVN(LCHNU)/DETU  
        CI12U = -CUN(LCHNU)/DETU  
        CI21U = -CVE(LCHNU)/DETU  
        CI22U = CUE(LCHNU)/DETU  

        DETV = CUE(LCHNV)*CVN(LCHNV)-CUN(LCHNV)*CVE(LCHNV)  
        CI11V = CVN(LCHNV)/DETV  
        CI12V = -CUN(LCHNV)/DETV  
        CI21V = -CVE(LCHNV)/DETV  
        CI22V = CUE(LCHNV)/DETV  

        ! *** X-DIRECTION CHANNEL  
        if( MDCHTYP(NMD) == 1 )then  
          if( QCHANU(NMD) > 0.0 )then  
            do K = 1,KC  
              QMCSINKX(LCHNU,K) = QMCSINKX(LCHNU,K)-0.5*DZC(L,K)*QCHANU(NMD)*(U(LCHNU,K)+U(LCHNU+1,K))  
              QMCSINKY(LCHNU,K) = QMCSINKY(LCHNU,K)-0.5*DZC(L,K)*QCHANU(NMD)*(V(LCHNU,K)+V(LNC(LCHNU),K))  
            enddo  
            do K = 1,KC  
              TMPVEC1(K) = CUE(LCHNU)*QMCSINKX(LCHNU,K)+CVE(LCHNU)*QMCSINKY(LCHNU,K)  
              TMPVEC2(K) = CUN(LCHNU)*QMCSINKX(LCHNU,K)+CVN(LCHNU)*QMCSINKY(LCHNU,K)  
            enddo  
            do K = 1,KC  
              QMCSOURX(LHOST,K) = QMCSOURX(LHOST,K)+CI11H*TMPVEC1(K)+CI12H*TMPVEC2(K)  
              QMCSOURY(LHOST,K) = QMCSOURY(LHOST,K)+CI21H*TMPVEC1(K)+CI22H*TMPVEC2(K)  
            enddo  
          else  
            do K = 1,KC  
              QMCSINKX(LHOST,K) = QMCSINKX(LHOST,K)+0.5*DZC(L,K)*QCHANU(NMD)*(U(LHOST,K)+U(LHOST+1,K))  
              QMCSINKY(LHOST,K) = QMCSINKY(LCHNU,K)+0.5*DZC(L,K)*QCHANU(NMD)*(V(LHOST,K)+V(LNC(LHOST),K))  
            enddo  
            do K = 1,KC  
              TMPVEC1(K) = CUE(LHOST)*QMCSINKX(LHOST,K)+CVE(LHOST)*QMCSINKY(LHOST,K)  
              TMPVEC2(K) = CUN(LHOST)*QMCSINKX(LCHNU,K)+CVN(LHOST)*QMCSINKY(LHOST,K)  
            enddo  
            do K = 1,KC  
              QMCSOURX(LCHNU,K) = QMCSOURX(LCHNU,K)-CI11U*TMPVEC1(K)-CI12U*TMPVEC2(K)  
              QMCSOURY(LCHNU,K) = QMCSOURY(LCHNU,K)-CI21U*TMPVEC1(K)-CI22U*TMPVEC2(K)  
            enddo  
          endif  
        endif  

        ! *** Y-DIRECTION CHANNEL  
        if( MDCHTYP(NMD) == 2 )then  
          if( QCHANV(NMD) > 0.0 )then  
            do K = 1,KC  
              QMCSINKX(LCHNV,K) = QMCSINKX(LCHNV,K)-0.5*DZC(L,K)*QCHANV(NMD)*(U(LCHNV,K)+U(LCHNV+1,K))  
              QMCSINKY(LCHNV,K) = QMCSINKY(LCHNV,K)-0.5*DZC(L,K)*QCHANV(NMD)*(V(LCHNV,K)+V(LNC(LCHNV),K))  
            enddo  
            do K = 1,KC  
              TMPVEC1(K) = CUE(LCHNV)*QMCSINKX(LCHNV,K)+CVE(LCHNV)*QMCSINKY(LCHNV,K)  
              TMPVEC2(K) = CUN(LCHNV)*QMCSINKX(LCHNV,K)+CVN(LCHNV)*QMCSINKY(LCHNV,K)  
            enddo  
            do K = 1,KC  
              QMCSOURX(LHOST,K) = QMCSOURX(LHOST,K)+CI11H*TMPVEC1(K)+CI12H*TMPVEC2(K)  
              QMCSOURY(LHOST,K) = QMCSOURY(LHOST,K)+CI21H*TMPVEC1(K)+CI22H*TMPVEC2(K)  
            enddo  
          else  
            do K = 1,KC  
              QMCSINKX(LHOST,K) = QMCSINKX(LHOST,K)+0.5*DZC(L,K)*QCHANV(NMD)*(U(LHOST,K)+U(LHOST+1,K))  
              QMCSINKY(LHOST,K) = QMCSINKY(LCHNV,K)+0.5*DZC(L,K)*QCHANV(NMD)*(V(LHOST,K)+V(LNC(LHOST),K))  
            enddo  
            do K = 1,KC  
              TMPVEC1(K) = CUE(LHOST)*QMCSINKX(LHOST,K)+CVE(LHOST)*QMCSINKY(LHOST,K)  
              TMPVEC2(K) = CUN(LHOST)*QMCSINKX(LCHNU,K)+CVN(LHOST)*QMCSINKY(LHOST,K)  
            enddo  
            do K = 1,KC  
              QMCSOURX(LCHNV,K) = QMCSOURX(LCHNV,K)-CI11V*TMPVEC1(K)-CI12V*TMPVEC2(K)  
              QMCSOURY(LCHNV,K) = QMCSOURY(LCHNV,K)-CI21V*TMPVEC1(K)-CI22V*TMPVEC2(K)  
            enddo  
          endif  
        endif  

      enddo  

      do K = 1,KC
        do L = 2,LA
          LE = LEC(L)
          LN = LNC(L)  
          if( QMCSOURX(L,K) /= 0.0 )then  
            TMPVAL = SUB3D(L,K)+SUB(LE)  
            TMPVAL = MAX(TMPVAL,1.0)  
            FX(L,K)  = FX(L,K)  - SUB3D(L,K)*QMCSOURX(L,K)/TMPVAL  
            FX(LE,K) = FX(LE,K) - SUB3D(LE,K)*QMCSOURX(L,K)/TMPVAL  
          endif  
          if( QMCSOURY(L,K) /= 0.0 )then  
            TMPVAL = SVB(L)+SVB(LN)  
            TMPVAL = MAX(TMPVAL,1.0)  
            FY(L,K)  = FY(L,K)  - SVB3D(L ,K)*QMCSOURX(L,K)/TMPVAL  
            FY(LN,K) = FY(LN,K) - SVB3D(LN,K)*QMCSOURX(L,K)/TMPVAL  
          endif  
          if( QMCSINKX(L,K) /= 0.0 )then  
            TMPVAL = SUB3D(L,K)+SUB(LE)  
            TMPVAL = MAX(TMPVAL,1.0)  
            FX(L,K)  = FX(L,K)  - SUB3D(L,K)*QMCSINKX(L,K)/TMPVAL  
            FX(LE,K) = FX(LE,K) - SUB3D(LE,K)*QMCSINKX(L,K)/TMPVAL  
          endif  
          if( QMCSINKY(L,K) /= 0.0 )then  
            TMPVAL = SVB(L)+SVB(LNC(L))  
            TMPVAL = MAX(TMPVAL,1.0)  
            FY(L,K)  = FY(L,K)  - SVB3D(L,K)*QMCSINKX(L,K)/TMPVAL  
            FY(LN,K) = FY(LN,K) - SVB3D(LN,K)*QMCSINKX(L,K)/TMPVAL  
          endif  
        enddo  
      enddo  
      !$OMP END SINGLE

    endif  ! *** END OF CHANNEL MODIFIER SECTION
    
  endif    ! *** END OF KC > 1 SECTION
  
  ! *** *******************************************************************!  
  !  
  ! ***  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR  
  ! ***  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP  
  ! ***  SBX = SBX*0.5*DYU & SBY = SBY*0.5*DXV  
  !  
  !----------------------------------------------------------------------!  

  ! *** ORIGINAL  
  if( BSC > 1.E-6 .and. KC > 1 )then

    if( IGRIDV == 1 )then
      ! *** SIGMA-ZED BOUYANCY SHEARS
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      do ND = 1,NDM  
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            LS = LSC(L) 
            LW = LWC(L)
            FBBX(L,K) = SUB3D(L,K)*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*SGZU(L,K+1) + (B(L,K)-B(LW,K))*SGZU(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = SVB3D(L,K)*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*SGZV(L,K+1) + (B(L,K)-B(LS,K))*SGZV(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          enddo  
        enddo  
      enddo
      !$OMP END DO
      
    elseif( IGRIDV > 1 )then
      ! *** SIGMA-ZED BOUYANCY SHEARS
      
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      do ND = 1,NDM
        ! *** ALL LAYERS
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND) 
            LS = LSC(L) 
            LW = LWC(L)
            FBBX(L,K) = SUB3D(L,K)*SBX(L)*GP*HU(L)*( BW(L,K+1)*HPW(L)*SGZW(L,K+1) - BE(LW,K+1)*HPE(LW)*SGZE(LW,K+1) + BW(L,K)*HPW(L)*SGZW(L,K) - BE(LW,K)*HPE(LW)*SGZE(LW,K)   &
                                                  - (BW(L,K+1)-BW(L,K)+BE(LW,K+1) - BE(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = SVB3D(L,K)*SBY(L)*GP*HV(L)*( BS(L,K+1)*HPS(L)*SGZS(L,K+1) - BN(LS,K+1)*HPN(LS)*SGZN(LS,K+1) + BS(L,K)*HPS(L)*SGZS(L,K) - BN(LS,K)*HPN(LS)*SGZN(LS,K)   &
                                                  - (BS(L,K+1)-BS(L,K)+BN(LS,K+1) - BN(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          enddo  
        enddo
      enddo
      !$OMP END DO
      
    elseif( IINTPG == 0 )then  
      ! *** IINTPG = 0  
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LS,LW) 
      do ND = 1,NDM  
        LF = (ND-1)*LDMWET+1  
        LL = MIN(LF+LDMWET-1,LAWET)
      
        do K = 1,KS  
          do LP = LF,LL
            L = LWET(LP)  
            LS = LSC(L)  
            LW = LWC(L)
            FBBX(L,K) = SBX(L)*GP*HU(L)*(HU(L)*( (B(L,K+1)-B(LW,K+1))*DZCK(K+1) + (B(L,K)-B(LW,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*(BELV(L)-BELV(LW)+Z(L,K)*(HP(L)-HP(LW))) )  
            FBBY(L,K) = SBY(L)*GP*HV(L)*(HV(L)*( (B(L,K+1)-B(LS,K+1))*DZCK(K+1) + (B(L,K)-B(LS,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*(BELV(L)-BELV(LS)+Z(L,K)*(HP(L)-HP(LS))) )  
          enddo  
        enddo  
      enddo
      !$OMP END DO

    elseif( IINTPG == 1. )then
      ! *** JACOBIAN
      
      if( KC <= 2 )then
        ! *** KC  <= 2 LAYERS
        
        !$OMP DO PRIVATE(ND,LF,LL,L,LP,K,LS,LW) 
        do ND = 1,NDM  
          LF = (ND-1)*LDMWET+1  
          LL = MIN(LF+LDMWET-1,LAWET)
      
          do K = 1,KS
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              LS = LSC(L)
              LW = LWC(L)
              FBBX(L,K) = SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K)-B(LW,K))*DZC(L,K) )-( B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K) )*(BELV(L)-BELV(LW)+Z(L,K)*(HP(L)-HP(LW))) )
              FBBY(L,K) = SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(L,K+1)+(B(L,K)-B(LS,K))*DZC(L,K) )-( B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K) )*(BELV(L)-BELV(LS)+Z(L,K)*(HP(L)-HP(LS))) )
            enddo
          enddo

        enddo   ! *** END OF DOMAIN
        !$OMP END DO
      
      else
        ! *** KC  > 2 LAYERS
        
        ! *** BOTTOM LAYER
        !$OMP DO PRIVATE(ND,LF,LL,L,LP,K,LS,LW) 
        do ND = 1,NDM  
          LF = (ND-1)*LDMWET+1  
          LL = MIN(LF+LDMWET-1,LAWET)
      
          do LP = LF,LL
            L = LWET(LP)  
            LW = LWC(L)
            LS = LSC(L)
            K = KSZ(L)
            FBBX(L,K) = SBX(L)*GP*HU(L)*( 0.5*HU(L)*( 0.5*(B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+2.0*(B(L,K  )-B(LW,K  ))*DZC(L,K  ) )-0.25*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.50*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW))) )
                      
            FBBY(L,K) = SBY(L)*GP*HV(L)*( 0.5*HV(L)*( 0.5*(B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+2.0*(B(L,K  )-B(LS ,K  ))*DZC(L,K  ) )-0.25*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))  &
                        *(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))  -0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))  &
                        *(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS))) )
          enddo
        enddo   ! *** END OF DOMAIN
        !$OMP END DO
      
        ! *** LAYER AT KC-1
        !$OMP DO PRIVATE(ND,LF,LL,L,LP,K,LS,LW) 
        do ND = 1,NDM  
          LF = (ND-1)*LDMWET+1  
          LL = MIN(LF+LDMWET-1,LAWET)
      
          K = KS
          do LP = LF,LL
            L = LWET(LP)  
            LW = LWC(L)
            LS = LSC(L)
            FBBX(L,K) = SBX(L)*GP*HU(L)*( 0.5*HU(L)*( 2.0*(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LW,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) )-0.50*(B(L,K+1)-B(L,K+1)+B(LW,K+1)-B(LW,K+1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW))) - 0.25*(B(L,K)-B(L,K-1)+B(LW,K)-B(LW,K-1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
                      
            FBBY(L,K) = SBY(L)*GP*HV(L)*( 0.5*HV(L)*( 2.0*(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LS ,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) )-0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))  &
                        *(BELV(L)-BELV(LS)+Z(L,K  )*(HP(L)-HP(LS)))   - 0.25*(B(L,K)-B(L,K-1)+B(LS,K)-B(LS ,K-1))   &
                        *(BELV(L)-BELV(LS)+Z(L,K-1)*(HP(L)-HP(LS))) )
          enddo
        enddo   ! *** END OF DOMAIN
        !$OMP END DO

        if( KC > 3 )then
        
          ! *** MIDDLE LAYERS

          !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
          do ND = 1,NDM  
            do K = 2,KS-1
              do LP = 1,LLWET(K-1,ND)
                L = LKWET(LP,K-1,ND) 
                LW = LWC(L)
                LS = LSC(L)
                FBBX(L,K) = SBX(L)*GP*HU(L)*( 0.5*HU(L)*( 0.5*(B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LW,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) ) - 0.25*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))  &
                           *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.50*(B(L,K+1)-B(L,K  )+B(LW,K+1)-B(LW,K  )) & 
                           *(BELV(L)-BELV(LW)+Z(L,K  )*(HP(L)-HP(LW)))-0.25*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1)) &
                           *(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
                         
                FBBY(L,K) = SBY(L)*GP*HV(L)*( 0.5*HV(L)*( 0.5*(B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LS ,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) )-0.25*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))    &
                           *(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))  &
                           *(BELV(L)-BELV(LS)+Z(L,K  )*(HP(L)-HP(LS)))-0.25*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))  &
                           *(BELV(L)-BELV(LS)+Z(L,K-1)*(HP(L)-HP(LS)))  )
              enddo
            enddo

          enddo   ! *** END OF DOMAIN
          !$OMP END DO

        endif     ! *** END OF KC > 3
      endif       ! *** END OF KC > 2
      ! *** END OF JACOBIAN SECTION

    elseif( IINTPG == 2 )then  
      ! *** FINITE VOLUME  
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      do ND = 1,NDM  
        do K = 1,KS  
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND)  
            LW = LWC(L)
            LS = LSC(L)
              
            FBBX(L,K) = SBX(L)*GP*HU(L)*( ( HP(L)*B(L,K+1)-HP(LW)*B(LW,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LW)*B(LW,K  ) )*DZC(L,K  ) )  & 
                       -SBX(L)*GP*(BELV(L)-BELV(LW))*( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LW)*B(LW,K+1)-HP(LW)*B(LW,K) )                  &
                       -SBX(L)*GP*(HP(L)-HP(LW))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K)+HP(LW)*ZZ(L,K+1)*B(LW,K+1)-HP(LW)*ZZ(L,K)*B(LW,K) )
                     
            FBBY(L,K) = SBY(L)*GP*HV(L)*( ( HP(L)*B(L,K+1)-HP(LS )*B(LS ,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LS )*B(LS ,K  ) )*DZC(L,K  ) )  &
                       -SBY(L)*GP*(BELV(L)-BELV(LS ))*( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LS)*B(LS ,K+1)-HP(LS)*B(LS ,K) )                    &
                       -SBY(L)*GP*(HP(L)-HP(LS ))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K) +HP(LS)*ZZ(L,K+1)*B(LS ,K+1)-HP(LS)*ZZ(L,K)*B(LS ,K) )
          enddo  
        enddo  
    
      enddo   ! *** END OF DOMAIN
      !$OMP END DO

    endif  
    
    ! *** BLOCKED LAYER FACE OPTION
    if( NBLOCKED > 0 .and. N > 1 )then 
      !$OMP SINGLE
      do LP = 1,NBLOCKED
        L = LBLOCKED(LP)
        
        if( KSZ(L) == KC ) CYCLE
        
        ! *** LAYER U BLOCKING
        if( BLDRAFTU(LP)+BLSILLU(LP) > 0.0 )then
          ! *** LAYER SPECIFIC DEPTH AVERAGED FLOWS (M3/S)
          if( SUB(L) == 0.0 ) CYCLE

          FBBX(L,:) = 0.0
          do K = KBBU(LP),KTBU(LP)-1
            LW = LWC(L)
            
            ! *** ADJUST ELEVATIONS AND DEPTHS
            if( BLSILLU(LP) > 0.0 )then
              ! *** SILL HEIGHT DICTATES BOTTOM ELEVATIONS
              BELVW(L)  = BELV(L) + BLSILLU(LP)
              BELVE(LW) = BELVW(L)
              HPW(L)  = HP(L)  - BLSILLU(LP)
              HPE(LW) = HP(LW) - BLSILLU(LP)
            else
              if( IGRIDV == 0 )then
                BELVW(L)  = BELV(L)
                BELVE(LW) = BELV(LW)
              endif
              HPW(L)  = HP(L)  - BLDRAFTU(LP)
              HPE(LW) = HP(LW) - BLDRAFTU(LP)
            endif

            FBBX(L,K) = SUB3D(L,K)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*SGZU(L,K+1) + (B(L,K)-B(LW,K))*SGZU(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
          enddo  
        endif
      
        ! *** LAYER V BLOCKING
        if( BLDRAFTV(LP)+BLSILLV(LP) > 0.0 )then
          ! *** LAYER SPECIFIC DEPTH AVERAGED FLOWS (M3/S)
          if( SVB(L) == 0.0 ) CYCLE

          FBBY(L,:) = 0.0
          do K = KBBV(LP),KTBV(LP)-1
            LS = LSC(L)  

            ! *** ADJUST ELEVATIONS AND DEPTHS
            if( BLSILLV(LP) > 0.0 )then
              ! *** SILL HEIGHT DICTATES BOTTOM ELEVATIONS
              BELVS(L)  = BELV(L) + BLSILLV(LP)
              BELVN(LS) = BELVS(L)
              HPS(L)  = HP(L)  - BLSILLV(LP)
              HPN(LS) = HP(LS) - BLSILLV(LP)
            else
              if( IGRIDV == 0 )then
                BELVS(L)  = BELV(L)
                BELVN(LS) = BELV(LS)
              endif
              HPS(L)  = HP(L)  - BLDRAFTV(LP)
              HPN(LS) = HP(LS) - BLDRAFTV(LP)
            endif

            FBBY(L,K) = SVB3D(L,K)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*SGZV(L,K+1) + (B(L,K)-B(LS,K))*SGZV(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          enddo  
        endif
      enddo
      !$OMP END SINGLE
    endif
  
  endif  ! *** END OF BOUYANCY 
  
  ! *** *******************************************************************!
  !
  ! ***  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
  !
  !----------------------------------------------------------------------!
  if( KC > 1 )then
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      
      ! *** COMPUTE THE INTERNAL SHEARS FOR THE LOWER LAYERS
      do K = 1,KS
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          DU(L,K) = CDZFU(L,K)*( HU(L)*( U(L,K+1)-U(L,K) )*DELTI + DXYIU(L)*( FCAX(L,K+1)-FCAX(L,K)   + FBBX(L,K) + SNLT*(FX(L,K)-FX(L,K+1)) ) )   ! *** M2/S2
          DV(L,K) = CDZFV(L,K)*( HV(L)*( V(L,K+1)-V(L,K) )*DELTI + DXYIV(L)*( FCAY(L,K)  -FCAY(L,K+1) + FBBY(L,K) + SNLT*(FY(L,K)-FY(L,K+1)) ) )   ! *** M2/S2
        enddo
      enddo

      ! *** ADD WIND SHEAR TO THE KC/KS INTERFACE
      if( NWSER > 0 )then
        do LP = 1,LLWET(KS,ND)
          L = LKWET(LP,KS,ND) 
          DU(L,KS) = DU(L,KS) - CDZUU(L,KS)*TSX(L)
          DV(L,KS) = DV(L,KS) - CDZUV(L,KS)*TSY(L)
        enddo
      endif
      
      ! *** Apply masks
      do K = 1,KS
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          DU(L,K) = SUB(L)*DU(L,K)
          DV(L,K) = SVB(L)*DV(L,K)
        enddo
      enddo      
    enddo    ! *** END OF DOMAIN
    !$OMP END DO

  endif
  !$OMP END PARALLEL

  ! *** *******************************************************************!
  return
  
END SUBROUTINE CALEXP2T 
