! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALEXP

  ! *** *******************************************************************!
  ! ***  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS  
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
  ! *** 2015-12     PAUL M. CRAIG      ADOPTED AQEA ISHDMF>0 FOR 3TL
  ! *** 2015-06     PAUL M. CRAIG      IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! *** 2015-02     Paul M. Craig      UPDATED OMP AND ADDED LDRY/LWET
  ! *** 2014-01     Paul M. Craig      Fixed the TOT FXVEGE/FYVEGE when partial vegetation penetration using VEGK
  ! *** 2012-09     Dang H Chung       Added OMP
  ! *** 2011-03     John Hamrick/Scott James
  ! ***                                Fixed the Vegetative Resistance from AVG to TOT using FKC
  ! *** 2010-10     Scott James        Added MHK
  !
  !----------------------------------------------------------------------C
  !
  ! *** *******************************************************************C
  !
  use GLOBAL
  use Allocate_Initialize
  use FIELDS
  use CYCLONE
  
  use Variables_MPI
  use Variables_Propwash
  use Variables_Ship
  use Mod_Active_Ship

  implicit none
  
  integer :: I, IOBC, LU, NS, ID,  JD, KD, LD, K,  L,  LL, LW, LNW, LSE, LE, LP  
  integer :: LEE, LWW, LNN, LSS
  integer :: ND, LF, NWR, IU, JU, KU, LN, LS   
  
  real :: QMF, QUMF, WUU, VTMPATU, UTMPATV, UMAGTMP, VMAGTMP, CACSUMT                                                                      
  real :: WVFACT, WVV, CFEFF, TMPVAL, DZPU, DZPV, FMDUY_TMP, FMDVX_TMP                                                       
  real :: VHC, VHB, DELTD2, UHC, UHB, QWRABS, VDIR                                                            

  real,save,allocatable,dimension(:)   :: CACSUM  
  real,save,allocatable,dimension(:,:) :: DZPC
  real,save,allocatable,dimension(:,:) :: FUHJ
  real,save,allocatable,dimension(:,:) :: FVHJ
  
  if( .not. allocated(DZPC) )then
    allocate(CACSUM(NTHREADS))  
    allocate(FUHJ(LCM,KCM))
    allocate(FVHJ(LCM,KCM))
    allocate(DZPC(LCM,KCM))
    CACSUM = 0.0
    FUHJ = 0.0
    FVHJ = 0.0
    DZPC = 0.0
  endif
  
  ! *** *******************************************************************C
  DELT = DT2
  DELTD2 = DT
  if( ISTL == 2 )then
    DELT = DT
    DELTD2 = 0.5*DT
  endif

  DELTI = 1./DELT

  if( N == 1 .and. DEBUG )then  
    open(1,FILE = OUTDIR//'MFLUX.DIA')  
    close(1,STATUS = 'DELETE')  
  endif  

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
        FUHU(L,K) = 0.0  
        FVHU(L,K) = 0.0  
        FUHV(L,K) = 0.0  
        FVHV(L,K) = 0.0  
        FWU(L,K)  = 0.0
        FWV(L,K)  = 0.0
        CAC(L,K)  = 0.0
        FCAX(L,K) = 0.0
        FCAY(L,K) = 0.0
        FX(L,K)   = 0.0
        FY(L,K)   = 0.0
        FBBX(L,K) = 0.0
        FBBY(L,K) = 0.0
        DU(L,K)   = 0.0  
        DV(L,K)   = 0.0  
        
        ! *** TWO LAYER ROTATIONAL EFFECTS OR WITHDRAWAL/RETURN
        FUHJ(L,K) = 0.0  
        FVHJ(L,K) = 0.0
        
        TVAR2E(L,K) = 0.0
        TVAR2N(L,K) = 0.0
      enddo
    enddo
    
    if( ISVEG > 0 )then
      do LP = 1,LADRY
        L = LDRY(LP)  
        FXVEGE(L) = 0.0
        FYVEGE(L) = 0.0
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
  ! *** *******************************************************************!  
  ! *** INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
  !----------------------------------------------------------------------!
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
  do ND = 1,NDM  
    LF = (ND-1)*LDMWET+1  
    LL = min(LF+LDMWET-1,LAWET)

    do LP = LF,LL
      L = LWET(LP)  
      FCAXE(L) = 0.0  
      FCAYE(L) = 0.0  
      FXE(L) = 0.0  
      FYE(L) = 0.0  
    enddo  
  enddo
  !$OMP END PARALLEL DO
  
  ! *** *******************************************************************C
  ! *** SELECT ADVECTIVE FLUX FORM
  ! ***

  if( ISTL == 2 )then
  
    ! *** THREE TIME LEVEL CORRECTOR STEP
    ! *** CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
    ! *** AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N

    !$OMP PARALLEL DEFAULT(SHARED)
    
    !$OMP DO PRIVATE(ND,K,LF,LL,L)
    do ND = 1,NDM
      LF = 2 + (ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do K = 1,KC
        do L = LF,LL
          TVAR2E(L,K) = 0.5*(UHDY(L,K) + UHDY1(L,K))
          TVAR2N(L,K) = 0.5*(VHDX(L,K) + VHDX1(L,K))
        enddo
      enddo
    enddo
    !$OMP END DO

    !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW,UHC,UHB,VHC,VHB)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)
          LS = LSC(L)
          LW = LWC(L)
          
          ! *** U COMPONENTS  
          UHB = 0.5*(TVAR2E(L,K) + TVAR2E(LE,K))
          VHC = 0.5*(TVAR2N(L,K) + TVAR2N(LW,K))

          ! ***       |-- EAST FLOWING --|    |-- WEST FLOWING --|
          FUHU(L,K) = (MAX(UHB,0.)*U1(L,K)  + min(UHB,0.)*U1(LE,K))
          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHU(L,K) = (MAX(VHC,0.)*U1(LS,K) + min(VHC,0.)*U1(L,K))

          ! *** V COMPONENTS
          VHB = 0.5*(TVAR2N(L,K) + TVAR2N(LN,K))
          UHC = 0.5*(TVAR2E(L,K) + TVAR2E(LS,K))
          
          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHV(L,K) = (MAX(VHB,0.)*V1(L,K)  + min(VHB,0.)*V1(LN, K))
          ! ***       |-- EAST  FLOWING --|   |-- WEST  FLOWING --|
          FUHV(L,K) = (MAX(UHC,0.)*V1(LW,K) + min(UHC,0.)*V1(L,K))
        enddo
      enddo
    enddo
    !$OMP END DO

    if( KC > 1 )then
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,WUU,WVV)
      do ND = 1,NDM
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L  = LKWET(LP,K,ND)  
            LW = LWC(L)
            LS = LSC(L)
            WUU = 0.25*DXYU(L)*(W2(L,K) + W2(LW,K))
            WVV = 0.25*DXYV(L)*(W2(L,K) + W2(LS,K))
            FWU(L,K) = max(WUU,0.)*U1(L,K) + min(WUU,0.)*U1(L,K+1)
            FWV(L,K) = max(WVV,0.)*V1(L,K) + min(WVV,0.)*V1(L,K+1)
          enddo
        enddo
      enddo
      !$OMP END DO
    endif
    !$OMP END PARALLEL
  
  elseif( ISTL /= 2 .and. ISCDMA == 0 )then

    ! *** THREE TIME LEVEL (LEAP-FROG) STEP
    ! *** WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
    ! *** UPWIND DIFFERENCE MOMENTUM ADVECTION 

    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW,UHC,UHB,VHC,VHB)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)
          LS = LSC(L)
          LW = LWC(L)

          ! *** U COMPONENTS  
          UHB = 0.5*(UHDY(L,K) + UHDY(LE,K))
          VHC = 0.5*(VHDX(L,K) + VHDX(LW,K))
          
          ! ***       |-- EAST FLOWING --|    |-- WEST FLOWING --|
          FUHU(L,K) = (MAX(UHB,0.)*U1(L,K)  + min(UHB,0.)*U1(LE,K))  
          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHU(L,K) = (MAX(VHC,0.)*U1(LS,K) + min(VHC,0.)*U1(L,K))
          
          ! *** V COMPONENTS
          VHB = 0.5*(VHDX(L,K)+VHDX(LN,K))
          UHC = 0.5*(UHDY(L,K)+UHDY(LS,K))

          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHV(L,K) = (MAX(VHB,0.)*V1(L,K)  + min(VHB,0.)*V1(LN, K))
          ! ***       |-- EAST  FLOWING --|   |-- WEST  FLOWING --|
          FUHV(L,K) = (MAX(UHC,0.)*V1(LW,K) + min(UHC,0.)*V1(L,K))
        enddo
      enddo
    enddo
    !$OMP END DO
    
    !----------------------------------------------------------------------!  
    ! *** COMPUTE VERTICAL ACCELERATIONS
    if( KC > 1 )then
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,WUU,WVV)
      do ND = 1,NDM

        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L  = LKWET(LP,K,ND)  
            LW = LWC(L)
            LS = LSC(L)
            WUU = 0.5*DXYU(L)*(W(L,K)+W(LW,K))
            WVV = 0.5*DXYV(L)*(W(L,K)+W(LS,K))
            FWU(L,K) = max(WUU,0.)*U1(L,K) + min(WUU,0.)*U1(L,K+1)
            FWV(L,K) = max(WVV,0.)*V1(L,K) + min(WVV,0.)*V1(L,K+1)
          enddo
        enddo
      enddo
      !$OMP END DO
    endif
    !$OMP END PARALLEL
  
  elseif( ISTL /= 2 .and. ISCDMA == 1 )then
    
    ! *** THREE TIME LEVEL (LEAP-FROG) STEP
    ! *** AT (N) AND TRANSPORTED FIELD AT (N)
    ! *** CENTRAL DIFFERENCE MOMENTUM ADVECTION
  
    !$OMP PARALLEL DEFAULT(SHARED) 
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)
          LS = LSC(L)
          LW = LWC(L)
          
          ! *** U COMPONENTS
          FUHU(L,K) = 0.25*(UHDY(L,K) + UHDY(LE,K))*(U(L,K)+U(LE,K))
          FVHU(L,K) = 0.25*(VHDX(L,K) + VHDX(LW,K))*(U(L,K)+U(LS,K))

          ! *** V COMPONENTS
          FVHV(L,K) = 0.25*(VHDX(L,K) + VHDX(LN,K))*(V(L,K)+V(LN,K))
          FUHV(L,K) = 0.25*(UHDY(L,K) + UHDY(LS,K))*(V(L,K)+V(LW,K))
        enddo
      enddo

      if( KC > 1 )then
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            LW = LWC(L)
            LS = LSC(L)
            FWU(L,K) = 0.25*DXYU(L)*(W(L,K)+W(LW,K))*(U(L,K+1)+U(L,K))
            FWV(L,K) = 0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
          enddo
        enddo
      endif
    enddo   ! ***  END OF DOMAIN
    !$OMP END DO
    !$OMP END PARALLEL   
  
  elseif( ISTL /= 2 .and. ISCDMA == 2 )then
    
    ! *** THREE TIME LEVEL (LEAP-FROG) STEP
    ! *** FIRST HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
    ! *** WITH TRANSPORT AT (N-1/2) AND TRANSPORTED FIELD AT (N-1)
    ! *** SECOND HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
    ! *** WITH TRANSPORT AT (N+1/2) AND TRANSPORTED FIELD AT (N)
    
    !$OMP PARALLEL DEFAULT(SHARED) 
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW,UHC,UHB,VHC,VHB,WUU,WVV)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)  
          U2(L,K) = U1(L,K)+U(L,K)
          V2(L,K) = V1(L,K)+V(L,K)
        enddo
      enddo

      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)
          LS = LSC(L)
          LW = LWC(L)

          ! *** U COMPONENTS
          UHB = 0.25*(UHDY(L,K)+UHDY(LE,K))
          VHC = 0.25*(VHDX(L,K)+VHDX(LW,K))
          
          FUHU(L,K) = max(UHB,0.)*U2(L,K)  + min(UHB,0.)*U2(LE,K)
          FVHU(L,K) = max(VHC,0.)*U2(LS,K) + min(VHC,0.)*U2(L,K)

          ! *** V COMPONENTS
          VHB = 0.25*(VHDX(L,K)+VHDX(LN,K))
          UHC = 0.25*(UHDY(L,K)+UHDY(LS,K))

          FVHV(L,K) = max(VHB,0.)*V2(L,K)  + min(VHB,0.)*V2(LN,K)
          FUHV(L,K) = max(UHC,0.)*V2(LW,K) + min(UHC,0.)*V2(L,K)
        enddo
      enddo
  
      if( KC > 1 )then
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L  = LKWET(LP,K,ND)  
            LW = LWC(L)
            LS = LSC(L)
            WUU = 0.25*DXYU(L)*(W(L,K)+W(LW,K))
            WVV = 0.25*DXYV(L)*(W(L,K)+W(LS,K))
            FWU(L,K) = max(WUU,0.)*U2(L,K) + min(WUU,0.)*U2(L,K+1)
            FWV(L,K) = max(WVV,0.)*V2(L,K) + min(WVV,0.)*V2(L,K+1)
          enddo
        enddo
      endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
 
  endif ! *** NUMERICAL SCHEMES
  
  ! *** Zero vertical momentum for open boundaries
  do IOBC = 1,NBCSOP2  
    L = LOBCS2(IOBC)  
    FWU(L,:) = 0.0
    FWV(L,:) = 0.0
  enddo  
    
  ! *** Add withdrawal/return flow momentum fluxes
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
    if( ABS(BCPS(LL).NQSMF) > 0 .and. BCPS(LL).NQSMF /= 5 )then  
      L = BCPS(LL).L
      do K = KSZ(L),KC
        ! *** Handle reversing flows in/out of domain
        if( QSERCELL(K,LL) >= 0. )then
          VDIR = 1.
        else
          VDIR = -1.
        endif
        !LD = LIJ(BCPS(LL).I,BCPS(LL).J)
        
        QMF = ABS(QSERCELL(K,LL))
        QUMF = VDIR*QMF*( QMF/( HPK(L,K)*BCPS(LL).QWIDTH ) )   ! *** M4/S2
        if( BCPS(LL).NQSMF ==  1 ) FUHJ(L     ,K) = -QUMF  
        if( BCPS(LL).NQSMF ==  2 ) FVHJ(L     ,K) = -QUMF  
        if( BCPS(LL).NQSMF ==  3 ) FUHJ(LEC(L),K) = -QUMF  
        if( BCPS(LL).NQSMF ==  4 ) FVHJ(LNC(L),K) = -QUMF  
        if( BCPS(LL).NQSMF == -1 ) FUHJ(L     ,K) = QUMF  
        if( BCPS(LL).NQSMF == -2 ) FVHJ(L     ,K) = QUMF  
        if( BCPS(LL).NQSMF == -3 ) FUHJ(LEC(L),K) = QUMF  
        if( BCPS(LL).NQSMF == -4 ) FVHJ(LNC(L),K) = QUMF  
      enddo
    endif  
  enddo

  ! *** Add propeller efflux velocity momentum
  if( ISPROPWASH == 2 .and. NACTIVESHIPS > 0 )then
    call add_ship_momentum(FUHJ, FVHJ)
  endif
  
  !$OMP PARALLEL DEFAULT(SHARED)     
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
  
  ! *** *******************************************************************C
  !
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
            CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K) + V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K) + U(L,K))*DXDJ(L) )*HP(L)
            CACSUM(ND) = CACSUM(ND)+CAC(L,K)
          enddo
        enddo
      enddo   ! ***  END OF DOMAIN
      !$OMP END DO
      
    else
      ! *** STANDARD CALCULATIONS, WITH DIAGNOSTICS
      !$OMP SINGLE
      do K = 1,KC
        do L = 2,LA
          LE = LEC(L)
          LN = LNC(L)
          CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K) + V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K) + U(L,K))*DXDJ(L) )*HP(L)
          CFEFF = ABS(CAC(L,K))*DXYIP(L)*HPI(L)
          CFMAX = max(CFMAX,CFEFF)
          CACSUM(1) = CACSUM(1)+CAC(L,K)
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

  ! *** *******************************************************************C
  !
  ! *** CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS
  !
  !----------------------------------------------------------------------C

  ! ***  STANDARD CALCULATION
  if( CACSUMT > 1.E-7 )then
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LE,LN,LS,LW,LNW,LSE)
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
          FCAX(L,K) = ROLD*FCAX(L,K) + 0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K) + V(L,K)) + CAC(LW,K)*(V(LNW,K) + V(LW,K)))
          FCAY(L,K) = ROLD*FCAY(L,K) + 0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(LE,K) + U(L,K)) + CAC(LS,K)*(U(LSE,K) + U(LS,K)))
        enddo
      enddo
    enddo
    !$OMP END DO

    !----------------------------------------------------------------------C
    !
    ! *** MODIFICATION FOR TYPE 2 OPEN BOUNDARIES
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

  endif    ! *** END OF CACSUMT > 1.E-7

  !----------------------------------------------------------------------C
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
    L = max(1,LBERC(LL))
    do K = KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K) + FUHJ(L,K)   ! ***  M4/S2
    enddo

    ! *** NORTH/SOUTH ADJACENT CELL
    L = max(1, LBNRC(LL))
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
  ! *** ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS
  !
  !----------------------------------------------------------------------!
  if( ISVEG > 0 )then
    ! *** ADD IN MHK DEVICES, IF NEEDED
    !$OMP SINGLE
    if( LMHK )CALL MHKPWRDIS
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
          FXVEG(L,K) = UMAGTMP*SUB3D(L,K)*FXVEG(L,K)  ![m/s] q_xC_d
          FYVEG(L,K) = VMAGTMP*SVB3D(L,K)*FYVEG(L,K)  ![m/s] q_yC_d

          !FXVEG/FXVEGE are multiplied by the local velocity to yield units of [m^4/s^2]
          !FXVEG/FXVEGE are added to the body forces as C_d(N/L^2)A|q|q
          FXVEGE(L) = FXVEGE(L) + FXVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FUHDXE
          FYVEGE(L) = FYVEGE(L) + FYVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FVHDYE
        enddo
      enddo

      ! *** ADD VEGETATIVE DRAG TO INTERNAL MODE SHEAR AND EXTERNAL MODE MOMENTUM
      if(ISVEG == 1 )then
        do K = 1,KC  
          do LP = 1,LLVEG(K,ND)
            L = LKVEG(LP,K,ND)
            if( (K-KSZ(L)+1) > INT(VEGK(L)+1.) ) CYCLE
            FX(L,K) = FX(L,K) + (FXVEG(L,K)-FXVEGE(L))*U(L,K)*DXYU(L) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) FXVEGE goes into FUHDXE for momentum conservation
            FY(L,K) = FY(L,K) + (FYVEG(L,K)-FYVEGE(L))*V(L,K)*DXYV(L) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) FYVEGE goes into FVHDYE for momentum conservation
          enddo
        enddo
      else
        do K = 1,KC
          do LP = 1,LLVEG(K,ND)
            L = LKVEG(LP,K,ND)
            if( (K-KSZ(L)+1) > INT(VEGK(L)+1.) ) CYCLE
            ! ***                m3/s      m/s   
            FX(L,K) = FX(L,K) + FXVEG(L,K)*U(L,K) ! *** (m^4/s^2)
            FY(L,K) = FY(L,K) + FYVEG(L,K)*V(L,K) ! *** (m^4/s^2)
          enddo  
        enddo
      endif
      
      ! *** CONVERT THE AVG FXVEGE/FYVEGE TO TOTAL FXVEGE/FYVEGE
      do LP = 1,LLVEG(KC,ND)
        L = LKVEG(LP,KC,ND) 
        FXVEGE(L) = FXVEGE(L)*HUI(L)*VEGK(L) !Calculate vegetative dissipation for FUHDYE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
        FYVEGE(L) = FYVEGE(L)*HVI(L)*VEGK(L) !Calculate vegetative dissipation for FVHDXE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
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
      
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif

  ! *** *******************************************************************!
  !
  ! *** ADD HORIZONTAL MOMENTUM DIFFUSION TO ADVECTIVE ACCELERATIONS
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
        if( max(DXP(LE),DYP(LE)) / min(DXP(LE),DYP(LE)) > 4.0 )then
          if( max(DXP(LEE),DYP(LEE)) / min(DXP(LEE),DYP(LEE)) > 4.0 )then
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
        if( max(DXP(LW),DYP(LW)) / min(DXP(LW),DYP(LW)) > 4.0 )then
          if( max(DXP(LWW),DYP(LWW)) / min(DXP(LWW),DYP(LWW)) > 4.0 )then
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
        if( max(DXP(LN),DYP(LN)) / min(DXP(LN),DYP(LN)) > 4.0 )then
          if( max(DXP(LNN),DYP(LNN)) / min(DXP(LNN),DYP(LNN)) > 4.0 )then
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
        if( max(DXP(LS),DYP(LS)) / min(DXP(LS),DYP(LS)) > 4.0 )then
          if( max(DXP(LSS),DYP(LSS)) / min(DXP(LSS),DYP(LSS)) > 4.0 )then
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
      
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW,FMDUY_TMP,FMDVX_TMP)
    do ND = 1,NDM  
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          ! *** THE FOLLOWING IS FROM AQEA 2015 CODE
          LS = LSC(L)
          LN = LNC(L)  
          LW = LWC(L)
          LE = LEC(L)
          ! *** CALCULATE DIFFUSIVE FLUX ON THE NORTHERN EXTERNAL BOUNDARY  
          if( SVB(LN) < 0.5  )then
            ! *** DXV1 IS 0 ON NORTHERN CLOSED BOUNDARY
            ! *** DYU1 IS CALCULATED USING THE U1 CLOSEST TO THE WALL AND 0
            FMDUY_TMP = DXU(L)*H1C(L)*AHC(L,K)*(-1.*U1(L,K)/DYU(L))
            FX(L,K) = FX(L,K) - SUB3D(L,K)*SDX(L)*( FMDUX(L,K) - FMDUX(LW,K) + FMDUY_TMP   - FMDUY(L,K) )
          else
            FX(L,K) = FX(L,K) - SUB3D(L,K)*SDX(L)*( FMDUX(L,K) - FMDUX(LW,K) + FMDUY(LN,K) - FMDUY(L,K) )
          endif
            
          ! *** CALCULATE DIFFUSIVE FLUX ON THE EASTERN EXTERNAL BOUNDARY  
          if( SUB(LE) < 0.5 )then
            ! *** DXV1 IS 0 IS CALCUALATED USING V1 CLOSEST TO THE WALL AND 0
            ! *** DYU1 IS 0 ON EASTERN CLOSED BOUNDARY
            FMDVX_TMP = (DYU(L))*H1C(L)*AHC(L,K)*((-1.*V1(L,K)/DXU(L)))
            FY(L,K) = FY(L,K) - SVB3D(L,K)*SDY(L)*( FMDVY(L,K) - FMDVY(LS,K) + FMDVX_TMP   - FMDVX(L,K) )
          else 
            FY(L,K) = FY(L,K) - SVB3D(L,K)*SDY(L)*( FMDVY(L,K) - FMDVY(LS,K) + FMDVX(LE,K) - FMDVX(L,K) )
          endif
        enddo  
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif
  
  ! *** *******************************************************************C
  !
  ! *** ADD BODY FORCE TO ADVECTIVE ACCELERATIONS
  ! *** DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF = 1
  ! *** DISTRIBUTE OVER SURFACE LAYER IF ISBODYF = 2
  !
  !----------------------------------------------------------------------C
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
      LL = min(LF+LDMWET-1,LAWET)
      do LP = LF,LL
        L = LWET(LP)
        FX(L,KC) = FX(L,KC) - DZIC(L,KC)*DYU(L)*HU(L)*FBODYFX(L,KC)
        FY(L,KC) = FY(L,KC) - DZIC(L,KC)*DXV(L)*HV(L)*FBODYFY(L,KC)
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif
  
  ! *** *******************************************************************C
  !
  ! *** ADD EXPLICIT NONHYDROSTATIC PRESSURE
  !
  if( KC > 1 .and. ISPNHYDS >= 1 )then
    !$OMP DO PRIVATE(ND,L,LP,K,TMPVAL) 
    do ND = 1,NDM  
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        TMPVAL    = 2./( DZC(L,KSZ(L)) + DZC(L,KSZ(L)+1) )
        DZPC(L,1) = TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))
      enddo
    enddo  
    !$OMP END DO
  
    !$OMP DO PRIVATE(ND,L,LP,TMPVAL) 
    do ND = 1,NDM  
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        TMPVAL     = 2./( DZC(L,KC) + DZC(L,KC-1) )
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
            TMPVAL    = 2./( DZC(L,K+1)+2.*DZC(L,K)+DZC(L,K-1) )
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
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*( HU(L)*(PNHYDS(L,K)-PNHYDS(LW,K) ) - ( BELV(L)-BELV(LW) + ZZ(L,K) *(HP(L)-HP(LW)) )*DZPU )
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
  !
  ! *** ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
  !
  if( ISWAVE == 2 .or. ISWAVE == 4 )then
  
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L) 
    do ND = 1,NDM  
      LF = (ND-1)*NWVCELLS+1  
      LL = min(LF+NWVCELLS-1,NWVCELLS)
      
      do K = 1,KC  
        do LP = LF,LL
          L = LWVCELL(LP)
          if( LKSZ(L,K) ) CYCLE  
          FX(L,K) = FX(L,K) + SUB3D(L,K)*WVFACT*SAAX(L)*FXWAVE(L,K)
          FY(L,K) = FY(L,K) + SVB3D(L,K)*WVFACT*SAAY(L)*FYWAVE(L,K)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  
  endif
  
  ! *** Filter cells which have been impacted by newly wet cells
  if( ISDRY > 0 )then
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM   
      do K = 1,KC  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          if( NWET(L) > 2*NTSTBC ) CYCLE
          
          FX(L,K) = FX1(L,K) + SIGN(MAX(1E-8,MIN( ABS(FX(L,K)-FX1(L,K)), ABS(2.0*FX1(L,K)) )), FX(L,K)-FX1(L,K) )
          FY(L,K) = FY1(L,K) + SIGN(MAX(1E-8,MIN( ABS(FY(L,K)-FY1(L,K)), ABS(2.0*FY1(L,K)) )), FY(L,K)-FY1(L,K) )
          FX1(L,K) = FX(L,K)
          FY1(L,K) = FY(L,K)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif
  
  ! *** *******************************************************************!
  !
  ! *** CALCULATE TOTAL EXTERNAL ACCELERATIONS
  !
  !----------------------------------------------------------------------!
  !$OMP DO PRIVATE(ND,K,LP,L)
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
    ! *** LIMIT THE VERTICAL MOMENTUM AT THE KSZ DIFFERENCES
  
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM  
      do K = 1,KC  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          FX(L,K) = FX(L,K) + SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(L,K)
          FY(L,K) = FY(L,K) + SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(L,K)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif
  
  ! *** *******************************************************************C
  !
  ! *** CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
  ! *** THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
  ! *** SBX = SBX*0.5*DYU & SBY = SBY*0.5*DXV
  !
  !----------------------------------------------------------------------C
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
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SUB3D(L,K)*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*SGZU(L,K+1) + (B(L,K)-B(LW,K))*SGZU(L,K) )                                &
                                                                                 - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SVB3D(L,K)*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*SGZV(L,K+1) + (B(L,K)-B(LS,K))*SGZV(L,K) )                                &
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
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SUB3D(L,K)*SBX(L)*GP*HU(L)*( BW(L,K+1)*HPW(L)*SGZW(L,K+1) - BE(LW,K+1)*HPE(LW)*SGZE(LW,K+1) + BW(L,K)*HPW(L)*SGZW(L,K) - BE(LW,K)*HPE(LW)*SGZE(LW,K)   &
                                                                        - (BW(L,K+1)-BW(L,K)+BE(LW,K+1) - BE(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SVB3D(L,K)*SBY(L)*GP*HV(L)*( BS(L,K+1)*HPS(L)*SGZS(L,K+1) - BN(LS,K+1)*HPN(LS)*SGZN(LS,K+1) + BS(L,K)*HPS(L)*SGZS(L,K) - BN(LS,K)*HPN(LS)*SGZN(LS,K)   &
                                                                        - (BS(L,K+1)-BS(L,K)+BN(LS,K+1) - BN(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          enddo  
        enddo
      enddo
      !$OMP END DO
      

    elseif( IINTPG == 0 )then  
      ! *** IINTPG = 0  
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LS,LW) 
      do ND = 1,NDM  
        do K = 1,KS  
        do LP = 1,LLWET(K,ND)
          L  = LKWET(LP,K,ND)  
            LS = LSC(L)
            LW = LWC(L)
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*DZCK(K+1) + (B(L,K)-B(LW,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*(BELV(L)-BELV(LW)+Z(L,K)*(HP(L)-HP(LW))) )
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZCK(K+1) + (B(L,K)-B(LS,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*(BELV(L)-BELV(LS)+Z(L,K)*(HP(L)-HP(LS))) )
            !FBBX(L,K) = SUBD(L)*FBBX(L,K)                                    
            !FBBY(L,K) = SVBD(L)*FBBY(L,K)                                    
          enddo
        enddo
      enddo
      !$OMP END DO  
  
    elseif( IINTPG == 1 )then
      ! *** JACOBIAN
      !$OMP SINGLE
      K = 1
      do L = 2,LA
        LW = LWC(L)
        LS = LSC(L)
        FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( 0.5*HU(L)*( (B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K  )-B(LW,K  ))*DZC(L,K  )+(B(L,K  )-B(LW,K  ))*DZC(L,K  ) ) &
                   -0.5*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))*(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.5*(B(L,K  )-B(L,K  )+B(LW,K  )-B(LW,K  ))*(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
  
        FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( 0.5*HV(L)*( (B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K  )-B(LS ,K  ))*DZC(L,K  )+(B(L,K  )-B(LS ,K  ))*DZC(L,K  ) ) &
                   -0.5*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.5*(B(L,K  )-B(L,K  )+B(LS ,K  )-B(LS ,K  ))*(BELV(L)-BELV(LS )+Z(L,K-1)*(HP(L)-HP(LS ))) )
      enddo
  
      if( KC > 2 )then
        K = KS
        do L = 2,LA
          LW = LWC(L)
          LS = LSC(L)
          FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( 0.5*HU(L)*( (B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K  )-B(LW,K  ))*DZC(L,K  )+(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) ) &
                     -0.5*(B(L,K+1)-B(L,K+1)+B(LW,K+1)-B(LW,K+1))*(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.5*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1))*(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
          FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( 0.5*HV(L)*( (B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K  )-B(LS ,K  ))*DZC(L,K  )+(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) ) &
                     -0.5*(B(L,K+1)-B(L,K+1)+B(LS ,K+1)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.5*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*(BELV(L)-BELV(LS )+Z(L,K-1)*(HP(L)-HP(LS ))) )
        enddo
      endif
  
      if( KC > 3 )then
        do K = 1,KS
          do L = 2,LA
            LW = LWC(L)
            LS = LSC(L)
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( 0.5*HU(L)*( (B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K  )-B(LW,K  ))*DZC(L,K  )+(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) ) &
                       -0.5*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))*(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.5*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1))*(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( 0.5*HV(L)*( (B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K  )-B(LS ,K  ))*DZC(L,K  )+(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) ) &
                       -0.5*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.5*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*(BELV(L)-BELV(LS )+Z(L,K-1)*(HP(L)-HP(LS ))) )
          enddo
        enddo
      endif
      !$OMP END SINGLE
        
    elseif( IINTPG == 2 )then
      ! *** FINITE VOLUME
      !$OMP SINGLE
      do K = 1,KS
        do L = 2,LA
          LW = LWC(L)
          LS = LSC(L)
          FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( ( HP(L)*B(L,K+1)-HP(LW)*B(LW,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LW)*B(LW,K  ) )*DZC(L,K  ) )-RNEW*SBX(L)*GP*(BELV(L)-BELV(LW)) &
                     *( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LW)*B(LW,K+1)-HP(LW)*B(LW,K) )   - RNEW*SBX(L)*GP*(HP(L)-HP(LW))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K)+HP(LW)*ZZ(L,K+1)*B(LW,K+1)-HP(LW)*ZZ(L,K)*B(LW,K) )
          FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( ( HP(L)*B(L,K+1)-HP(LS )*B(LS ,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LS )*B(LS ,K  ) )*DZC(L,K  ) )-RNEW*SBY(L)*GP*(BELV(L)-BELV(LS )) &
                     *( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LS)*B(LS ,K+1)-HP(LS)*B(LS ,K) ) - RNEW*SBY(L)*GP*(HP(L)-HP(LS ))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K)+HP(LS)*ZZ(L,K+1)*B(LS ,K+1)-HP(LS)*ZZ(L,K)*B(LS ,K) )
        enddo
      enddo
      !$OMP END SINGLE
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
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K)
    do ND = 1,NDM  
      ! *** COMPUTE THE INTERNAL SHEARS FOR THE LOWER LAYERS
      do K = 1,KS
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          DU(L,K) = CDZFU(L,K)*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI + DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)   + FBBX(L,K) + SNLT*(FX(L,K)-FX(L,K+1))) )
          DV(L,K) = CDZFV(L,K)*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI + DXYIV(L)*(FCAY(L,K)  -FCAY(L,K+1) + FBBY(L,K) + SNLT*(FY(L,K)-FY(L,K+1))) )
        enddo
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO
  endif

  ! *** ADD WIND SHEAR TO THE KC/KS INTERFACE
  if( (ISTL == 2 .and. NWSER > 0) .or. (ISTL == 2 .and. iGOTM_Test > 0 ) )then
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)

      do LP = LF,LL
        L = LWET(LP)  
        DU(L,KS) = DU(L,KS) - CDZUU(L,KS)*TSX(L)
        DV(L,KS) = DV(L,KS) - CDZUV(L,KS)*TSY(L)
      enddo
      
    enddo    ! *** END OF DOMAIN
    !$OMP END DO
  endif

  ! *** Apply masks
  !$OMP DO PRIVATE(ND,LP,K,L)
  do ND = 1,NDM  
    do K = 1,KS
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        DU(L,K) = SUB(L)*DU(L,K)
        DV(L,K) = SVB(L)*DV(L,K)
      enddo
    enddo
  enddo    ! *** END OF DOMAIN
  !$OMP END DO
  
  !$OMP END PARALLEL
       
  ! *** *******************************************************************!

  return
  
END SUBROUTINE CALEXP
