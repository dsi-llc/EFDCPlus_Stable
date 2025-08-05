! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALHDMF

  ! *** CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND
  ! *** DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING
  ! *** SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO

  ! *** ONLY VALID FOR ISHDMF >= 1
  !
  !----------------------------------------------------------------------------------------C
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------!
  !    2019-06       PAUL M. CRAIG     CHANGED APPROACH FOR WALL ROUGHNESS AND FMDUY/FMDVX
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2011-05       Paul M. Craig     Corrected DSQR equation from /4 to /2
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP
  !    2008-10       SANG YUK          CORRECTED THE DIFFUSIVE MOMENTUM FLUXES COMPUTATION
  !    2004-11       PAUL M. CRAIG     REWRITTEN AND RESTRUCTURED

  use GLOBAL
  use Allocate_Initialize

  implicit none

  integer :: L, LW, K, LL, NQSTMP, IU, JU, KU, NWR, ND, LN, LS, LP, LE, LG

  real      :: SLIPCO, TMPVAL, DSQR, WVFACT
  real      :: DTMPH, DTMPX, AHWVX, DX2DZBR, DY2DZBR, CSDRAG, SLIPFAC, FACES

  real,save,allocatable,dimension(:,:) :: AHEE
  real,save,allocatable,dimension(:,:) :: AHNN
  real,save,allocatable,dimension(:,:) :: SXY
  real,save,allocatable,dimension(:,:) :: DYU1
  real,save,allocatable,dimension(:,:) :: DYV1
  real,save,allocatable,dimension(:,:) :: DXU1
  real,save,allocatable,dimension(:,:) :: DXV1
  real,save,allocatable,dimension(:,:) :: FMDUX0, FMDUY0, FMDVY0, FMDVX0

  if( .not. allocated(AHEE) )then
    call AllocateDSI( AHEE,   LCM, KCM, 0.0)
    call AllocateDSI( AHNN,   LCM, KCM, 0.0)
    call AllocateDSI( SXY,    LCM, KCM, 0.0)
    call AllocateDSI( DYU1,   LCM, KCM, 0.0)
    call AllocateDSI( DYV1,   LCM, KCM, 0.0)
    call AllocateDSI( DXU1,   LCM, KCM, 0.0)
    call AllocateDSI( DXV1,   LCM, KCM, 0.0)
    call AllocateDSI( FMDUX0, LCM, KCM, 0.0)
    call AllocateDSI( FMDUY0, LCM, KCM, 0.0)
    call AllocateDSI( FMDVY0, LCM, KCM, 0.0)
    call AllocateDSI( FMDVX0, LCM, KCM, 0.0)
  endif
  
  ! *** SXX+SYY DEFINED AT CELL CENTERS AND STORED IN DXU1(L,K)
  SLIPCO = 1.
  if( AHD > 0.0 )then
    SLIPCO = 0.5/SQRT(AHD)
  endif

  ! ****************************************************************************
  ! *** RESET CELL WHEN INITIALLY DRY
  if( LADRY > 0 )then
    do K = 1,KC
      do LP = 1,LADRY
        L = LDRY(LP)
        AH(L,K)  = AHOXY(L)
        AHC(L,K) = AHOXY(L)
        DXU1(L,K) = 0.0
        DXV1(L,K) = 0.0
        DYU1(L,K) = 0.0
        DYV1(L,K) = 0.0
        FMDUX(L,K) = 0.0
        FMDUY(L,K) = 0.0
        FMDVY(L,K) = 0.0
        FMDVX(L,K) = 0.0
      enddo
    enddo
  endif

  !$OMP PARALLEL DEFAULT(SHARED)

  if( ISDRY > 0 .and. LADRY > 0 )then
    !$OMP DO PRIVATE(ND,K,LN,LP,L)
    do ND = 1,NDM
      do K = 1,KC
        LN = 0
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          if( LHDMF(L,K) )then
            LN = LN+1
            LKHDMF(LN,K,ND) = L
          endif
        enddo
        LLHDMF(K,ND) = LN    ! *** NUMBER OF WET HDMF CELLS FOR THE CURRENT LAYER
      enddo
    enddo
    !$OMP END DO
  endif

  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW,TMPVAL,DX2DZBR,DY2DZBR,CSDRAG,SLIPFAC,FACES)
  do ND = 1,NDM
    ! ***  CALCULATE HORIZONTAL VELOCITY SHEARS
    do K = 1,KC
      do LP = 1,LLHDMF(K,ND)
        L = LKHDMF(LP,K,ND)
        LE = LEC(L)
        LN = LNC(L)
        ! *** DXU1 = dU/dX, UNITS: 1/S
        DXU1(L,K) = ( U(LE,K) - U(L,K) )/DXP(L)
        ! *** DYV1 = dV/dY, UNITS: 1/S
        DYV1(L,K) = ( V(LN,K) - V(L,K) )/DYP(L)
      enddo
    enddo
    if( ISHDMF == 1 .or. ISHDMF == 2 )then
      ! *** HMD WITHOUT WALL EFFECTS

      ! *** DYU1 = dU/dY
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          LS = LSC(L)
          DYU1(L,K) = ( U(L,K) - U(LS,K) )/DYV(L)
        enddo
      enddo

      ! *** DXV1 = dV/dX
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          LW = LWC(L)
          DXV1(L,K) = ( V(L,K) - V(LW,K) )/DXU(L)
        enddo
      enddo

    else
      ! *** HMD WITH WALL EFFECTS

      ! *** DYU1 = dU/dY, DXV1 = dV/dX
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          FACES = SUB3D(L,K) + SUB3D(LEC(L),K) + SVB3D(L,K) + SVB3D(LNC(L),K)
          if( FACES > 3.0 )then
            ! *** OPEN WATER - AVERAGE SHEARS FROM BOTH SIDES OF THE CELL
            DYU1(L,K) = 0.5*( ( U(L,K) - U(LSC(L),K) )/DYV(L) + ( U(L,K) - U(LNC(L),K) )/DYV(LNC(L)) )      ! *** DYU1 = dU/dY
            DXV1(L,K) = 0.5*( ( V(L,K) - V(LWC(L),K) )/DXU(L) + ( V(L,K) - V(LEC(L),K) )/DXU(LEC(L)) )      ! *** DXV1 = dV/dX
          else
            ! *** NORTH/SOUTH WALLS
            FACES = SUB3D(L,K) + SUB3D(LEC(L),K)
            if( FACES > 0.5 )then
              if( SVB3D(L,K) > 0.5 .XOR. SVB3D(LNC(L),K) > 0.5 )then
                ! *** NORTH OR SOUTH FACE IS A WALL
                DY2DZBR = 1. + 0.5*DYP(L)/ZBRWALL
                CSDRAG  = 0.16/((LOG(DY2DZBR))**2)
                SLIPFAC  = SLIPCO*CSDRAG
                DYU1(L,K) = SLIPFAC*U(L,K)/DYP(L)
              endif
            endif

            ! *** EAST/WEST WALLS
            FACES = SVB3D(L,K) + SVB3D(LNC(L),K)
            if( FACES > 0.5 )then
              if( SUB3D(L,K) > 0.5 .XOR. SUB3D(LEC(L),K) > 0.5 )then
                ! *** EAST OR WEST FACE IS A WALL
                DY2DZBR = 1. + 0.5*DYP(L)/ZBRWALL
                CSDRAG  = 0.16/((LOG(DY2DZBR))**2)
                SLIPFAC  = SLIPCO*CSDRAG
                DXV1(L,K) = SLIPFAC*V(L,K)/DXP(L)
              endif
            endif
          endif
        enddo
      enddo
    endif
  enddo  ! *** END OF DOMAIN
  !$OMP END DO

  ! *** WITHDRAWAL/RETURN
  if( NQWR > 0 )then
    !$OMP SINGLE
    do NWR = 1,NQWR
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NQSTMP = WITH_RET(NWR).NQWRSERQ
      if( QWRSERT(NQSTMP) >= 0. )then
        ! *** Original Withdrawal/Return
        IU = WITH_RET(NWR).IQWRU
        JU = WITH_RET(NWR).JQWRU
        KU = WITH_RET(NWR).KQWRU
      else
        ! *** Reverse Flow Withdrawal/Return
        IU = WITH_RET(NWR).IQWRD
        JU = WITH_RET(NWR).JQWRD
        KU = WITH_RET(NWR).KQWRD
      endif
      DXU1(LIJ(IU,JU),KU) = 0.0
      DXV1(LIJ(IU,JU),KU) = 0.0
      DYU1(LIJ(IU,JU),KU) = 0.0
      DYV1(LIJ(IU,JU),KU) = 0.0
    enddo
    !$OMP END SINGLE
  endif

  ! *** SXY = dU/dY + dV/dX
  !$OMP DO PRIVATE(ND,K,LP,L)
  do ND = 1,NDM
    do K = 1,KC
      do LP = 1,LLHDMF(K,ND)
        L = LKHDMF(LP,K,ND)
        SXY(L,K) = DYU1(L,K) + DXV1(L,K)
      enddo
    enddo
  enddo
  !$OMP END DO

  !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,LN)  &
  !$OMP    PRIVATE(TMPVAL,DSQR,WVFACT,AHWVX,DTMPH,DTMPX)
  do ND = 1,NDM
    if( AHD > 0.0 )then
      ! *** CALCULATE SMAGORINSKY HORIZONTAL VISCOSITY
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          TMPVAL = AHDXY(L)*DXP(L)*DYP(L)
          DSQR = DXU1(L,K)*DXU1(L,K) + DYV1(L,K)*DYV1(L,K) + 0.5*SXY(L,K)*SXY(L,K)
          AH(L,K) = AHOXY(L) + TMPVAL*SQRT(DSQR)
        enddo
      enddo
    elseif( N < 10 .or. ISWAVE == 2 .or. ISWAVE == 4 )then
      ! *** ONLY NEED TO ASSIGN INITIALLY
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          AH(L,K) = AHOXY(L)
        enddo
      enddo
    endif

    ! ***  CALCULATE HORIZONTAL SMAG DIFFUSION DUE TO WAVE BREAKING
    if( ISWAVE == 2 .or. ISWAVE == 4 )then
      if( WVLSH > 0.0 .or. WVLSX > 0.0 )then
        if( ISWAVE == 2 .and. N < NTSWV )then
          TMPVAL = FLOAT(N)/FLOAT(NTSWV)
          WVFACT = 0.5-0.5*COS(PI*TMPVAL)
        else
          WVFACT = 1.0
        endif

        if( ISDRY > 0 )then
          do K = 1,KC
            do LP = 1,LLHDMF(K,ND)
              L = LKHDMF(LP,K,ND)
              if( LWVMASK(L) )then
                if( LMASKDRY(L) )then
                  if( WV(L).DISSIPA(K) > 0 )then
                    DTMPH = WV(L).DISSIPA(K)**0.3333
                  else
                    DTMPH = 0
                  endif
                  TMPVAL = 2.*PI/WV(L).FREQ     ! *** WAVE PERIOD
                  AHWVX = WVLSX*TMPVAL*TMPVAL
                  DTMPX = WV(L).DISSIPA(K)/HP(L)
                  AH(L,K) = AH(L,K)+WVFACT*(WVLSH*DTMPH*HP(L)+AHWVX*DTMPX)
                endif
              endif
            enddo
          enddo
        else
          do K = 1,KC
            do LP = 1,LLHDMF(K,ND)
              L = LKHDMF(LP,K,ND)
              if( LWVMASK(L) )then
                if( WV(L).DISSIPA(K) > 0 )then
                  DTMPH = WV(L).DISSIPA(K)**0.3333
                else
                  DTMPH = 0
                endif
                TMPVAL = 2.*PI/WV(L).FREQ
                AHWVX = WVLSX*TMPVAL*TMPVAL
                DTMPX = WV(L).DISSIPA(K)/HP(L)
                AH(L,K) = AH(L,K)+WVFACT*(WVLSH*DTMPH*HP(L)+AHWVX*DTMPX)
              endif
            enddo
          enddo
        endif
      endif
    endif
  enddo  ! *** END OF DOMAIN
  !$OMP END DO

  !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LW,LN)
  do ND = 1,NDM

    ! ***  CALCULATE DIFFUSIVE MOMENTUM FLUXES
    if( ISHDMF == 1 .or. ISHDMF == 2 )then
      ! *** NO WALL EFFECTS
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          LE = LEC(L)
          LN = LNC(L)
          LS = LSC(L)
          LW = LWC(L)
          FMDUX0(L,K) = ( DYP(L) *HP(L) *AH(L,K) *DXU1(L,K) - DYP(LW)*HP(LW)*AH(LW,K)*DXU1(LW,K) )*SUB(L)*SUB(LW)
          FMDUY0(L,K) = ( DXU(LN)*HU(LN)*AH(LN,K)*SXY(LN,K) - DXU(L) *HU(L) *AH(L,K) *SXY(L,K)   )*SVB(LW)*SVB(L)*SUB(LS)*SUB(L)
          
          FMDVY0(L,K) = ( DXP(L) *HP(L) *AH(L,K) *DYV1(L,K) - DXP(LS)*HP(LS)*AH(LS,K)*DYV1(LS,K) )*SVB(L)*SVB(LN)
          FMDVX0(L,K) = ( DYV(LE)*HV(LE)*AH(LE,K)*SXY(LE,K) - DYV(L) *HV(L) *AH(L,K) *SXY(L,K)   )*SVB(LW)*SVB(L)*SUB(LS)*SUB(L)
        enddo
      enddo
    else
      ! *** WITH WALL EFFECTS
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          LE = LEC(L)
          LN = LNC(L)
          LS = LSC(L)
          LW = LWC(L)
          FMDUX0(L,K) = ( DYP(L) *HP(L) *AH(L,K) *DXU1(L,K) - DYP(LW)*HP(LW)*AH(LW,K)*DXU1(LW,K) )
          FMDUY0(L,K) = (                                   - DXU(L) *HU(L) *AH(L,K) *SXY(L,K)   )
               
          FMDVY0(L,K) = ( DXP(L) *HP(L) *AH(L,K) *DYV1(L,K) - DXP(LS)*HP(LS)*AH(LS,K)*DYV1(LS,K) )
          FMDVX0(L,K) = (                                   - DYV(L) *HV(L) *AH(L,K) *SXY(L,K)   )
        enddo
      enddo
    endif
  enddo  ! *** END OF DOMAIN
  !$OMP END DO

  ! *** Advance the horizontal diffusion terms
  if( ISHDMFILTER > 0 .and. ((ISRESTI > 0 .and. NITER > 1) .or. NITER > 100) )then
    !$OMP DO PRIVATE(ND,K,LP,L) 
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L  = LKHDMF(LP,K,ND)  
          FMDUX(L,K) = FMDUX(L,K) + SIGN(MAX(1E-8,MIN( ABS(FMDUX0(L,K)-FMDUX(L,K)), ABS(0.25*FMDUX(L,K)) )), FMDUX0(L,K)-FMDUX(L,K) )                                             ! SQRT(FMDUX(L,K)*MAX(FMDUX0,1E-8))
          FMDUY(L,K) = FMDUY(L,K) + SIGN(MAX(1E-8,MIN( ABS(FMDUY0(L,K)-FMDUY(L,K)), ABS(0.25*FMDUY(L,K)) )), FMDUY0(L,K)-FMDUY(L,K) )                                             ! SQRT(FMDUX(L,K)*MAX(FMDUX0,1E-8))
          FMDVY(L,K) = FMDVY(L,K) + SIGN(MAX(1E-8,MIN( ABS(FMDVY0(L,K)-FMDVY(L,K)), ABS(0.25*FMDVY(L,K)) )), FMDVY0(L,K)-FMDVY(L,K) )                                             ! SQRT(FMDUX(L,K)*MAX(FMDUX0,1E-8))
          FMDVX(L,K) = FMDVX(L,K) + SIGN(MAX(1E-8,MIN( ABS(FMDVX0(L,K)-FMDVX(L,K)), ABS(0.25*FMDVX(L,K)) )), FMDVX0(L,K)-FMDVX(L,K) )                                             ! SQRT(FMDUX(L,K)*MAX(FMDUX0,1E-8))
        enddo
      enddo
    enddo  ! *** END OF DOMAIN  
    !$OMP END DO
  else
    !$OMP DO PRIVATE(ND,K,LP,L) 
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L  = LKHDMF(LP,K,ND)  
          FMDUX(L,K) = FMDUX0(L,K)
          FMDUY(L,K) = FMDUY0(L,K)
          FMDVY(L,K) = FMDVY0(L,K)
          FMDVX(L,K) = FMDVX0(L,K)
        enddo
      enddo
    enddo  ! *** END OF DOMAIN  
    !$OMP END DO
  endif

  if( ISDRY > 0 .and. NASPECT > 0 )then
    ! *** ZERO XY COMPONENT FOR CELLS WITH HIGH ASPECT RATIOS
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          if( LASPECT(L) )then
            FMDUY(L,K) = 0.
            FMDVX(L,K) = 0.
          endif
        enddo
      enddo
    enddo
    !$OMP END DO
  endif
  !$OMP END PARALLEL

  ! *** ZERO BOUNDARY CELL MOMENTUM DIFFUSION
  do LL = 1,NBCSOP
    L = LOBCS(LL)
    do K = 1,KC
      FMDUX(L,K) = 0.0
      FMDUY(L,K) = 0.0
      FMDVY(L,K) = 0.0
      FMDVX(L,K) = 0.0
    enddo
  enddo

  return

END
