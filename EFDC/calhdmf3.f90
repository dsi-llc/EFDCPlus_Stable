! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALHDMF3

  ! ***  SUBROUTINE CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND
  ! ***  DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING
  ! ***  SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO
  !
  ! *** *******************************************************************!
  ! *** SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS

  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! *** 2015-12     PAUL M. CRAIG      ADDED OMP AND UPDATED TO F90
  ! *** 2015-12     PAUL M. CRAIG      ADOPTED AQEA ISHDMF>0 FOR 3TL

  ! *** *******************************************************************!
  use GLOBAL  
  use Allocate_Initialize
  
  !  use OMP_LIB

  implicit none

  ! *** *******************************************************************!
  !
  ! **VARIABLE DEFINITION 
  integer :: L,  LF, LL, LP, LW, K, ND, LN, LS, LE, LNW, LNE, LSW, LSE
  real    :: TMPVAL, WVFACT, DTMPH, DTMPX, AHWVX, DTMP
  real    :: HPLW, HPLS, HPLSW
  
  real,save,allocatable,dimension(:,:) :: DYU1
  real,save,allocatable,dimension(:,:) :: DYV1
  real,save,allocatable,dimension(:,:) :: DXU1
  real,save,allocatable,dimension(:,:) :: DXV1
  real,save,allocatable,dimension(:,:) :: FMDUX0, FMDUY0, FMDVY0, FMDVX0
  real,save,allocatable,dimension(:)   :: HMC

  if( .not. allocated(DYU1) )then
    call AllocateDSI( DYU1,   LCM, KCM, 0.0)
    call AllocateDSI( DYV1,   LCM, KCM, 0.0)
    call AllocateDSI( DXU1,   LCM, KCM, 0.0)
    call AllocateDSI( DXV1,   LCM, KCM, 0.0)
    call AllocateDSI( FMDUX0, LCM, KCM, 0.0)
    call AllocateDSI( FMDUY0, LCM, KCM, 0.0)
    call AllocateDSI( FMDVY0, LCM, KCM, 0.0)
    call AllocateDSI( FMDVX0, LCM, KCM, 0.0)
    call AllocateDSI( HMC,    LCM, 0.0)

    FMDUX = 1E-8
    FMDUY = 1E-8
    FMDVY = 1E-8
    FMDVX = 1E-8
  endif
  
  if( ISDRY > 0 )then
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
          !FMDUX(L,K) = 1E-16
          !FMDUY(L,K) = 1E-16
          !FMDVY(L,K) = 1E-16
          !FMDVX(L,K) = 1E-16
        enddo  
      enddo 
    endif
  endif
  !IF( ISTL /= 3 ) return
  
  if( ISWAVE == 2 .or. ISWAVE == 4 )then
    if( WVLSH > 0.0 .or. WVLSX > 0.0 )then
      if( N < NTSWV )then
        TMPVAL = FLOAT(N)/FLOAT(NTSWV)
        WVFACT = 0.5-0.5*COS(PI*TMPVAL)
      else
        WVFACT = 1.0
      endif
      AHWVX = WVFACT*WVLSX*WVPRD*WVPRD
    endif
  endif
  
  ! ****************************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  
  if( ISDRY > 0 .and. LADRY > 0 )then
    !$OMP DO PRIVATE(ND,K,L,LN,LF,LL)
    do ND = 1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = MIN(LF+LDM-1,LA)
      do K = 1,KC  
        LN = 0
        do L = LF,LL
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
  
  !----------------------------------------------------------------------!
  ! ***  CALCULATE HORIZONTAL VELOCITY SHEARS
  !$OMP DO PRIVATE(ND,K,LP,L,LS,LN,LE,LW) 
  do ND = 1,NDM
    do K = 1,KC
      do LP = 1,LLHDMF(K,ND)
        L = LKHDMF(LP,K,ND)  
        LS = LSC(L)
        LN = LNC(L)

        LE = LEC(L)
        LW = LWC(L)

        ! *** SHEAR ACROSS CELL CENTERS
        DXU1(L,K) = (U1(LE,K)-U1(L,K))/DXP(L)
        DYV1(L,K) = (V1(LN,K)-V1(L,K))/DYP(L)
        
        ! *** HORIZONTAL X COMPONENT SHEAR AT SW CORNER
        DYU1(L,K) = 2.*SUBO(LS)*(U1(L,K)-U1(LS,K))/(DYU(L)+DYU(LS))
        if( SVBO(L) < 0.5 ) DYU1(L,K) = 2.0*U1(L,K)/DYU(L)

        ! *** HORIZONTAL Y COMPONENT SHEAR AT SW CORNER
        DXV1(L,K) = 2.*SVBO(LW)*(V1(L,K)-V1(LW,K))/(DXV(L)+DXV(LW))
        if( SUBO(L) < 0.5 ) DXV1(L,K) = 2.0*V1(L,K)/DXV(L)
        
      enddo
    enddo
  enddo   ! *** END OF DOMAIN 
  !$OMP END DO
  
  !----------------------------------------------------------------------!
  ! ***  CALCULATE HORIZONTAL VISCOSITY
  if( AHD > 0.0 )then
    !$OMP DO PRIVATE(ND,K,LP,L,LN,LS,LE,LW,LNW,LNE,LSW,LSE) 
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L   = LKHDMF(LP,K,ND)  
          LN  = LNC(L)
          LS  = LSC(L)
          LE  = LEC(L)
          LW  = LWC(L)
          LNW = LNWC(L)
          LNE = LNEC(L)
          LSW = LSWC(L)
          LSE = LSEC(L)

          ! *** CELL CENTROID
          AH(L,K)  = AHOXY(L) + AHDXY(L)*DXP(L)*DYP(L)                                                     &
                                  *SQRT( 2.*DXU1(L,K)*DXU1(L,K)  + 2.*DYV1(L,K)*DYV1(L,K)                  &
                                        + 0.0625*(DYU1(L,K) + DYU1(LN,K) + DYU1(LE,K) + DYU1(LNE,K)        &
                                        +         DXV1(L,K) + DXV1(LN,K) + DXV1(LE,K) + DXV1(LNE,K))**2. )
          ! *** SW CORNER
          AHC(L,K) = AHOXY(L) + AHDXY(L)*0.0625*( (DXP(L) + DXP(LW) + DXP(LS) + DXP(LSW))**2. )            &
                                  *SQRT(  0.125*(DXU1(L,K) + DXU1(LW,K) + DXU1(LS,K) + DXU1(LSW,K))**2.    &
                                        + 0.125*(DYV1(L,K) + DYV1(LW,K) + DYV1(LS,K) + DYV1(LSW,K))**2.    &
                                        + (DYU1(L,K)+ DXV1(L,K))**2. )
        enddo
      enddo
    enddo  ! *** END OF DOMAIN  
    !$OMP END DO
    

  elseif( N < 10 .or. ISWAVE == 2 .or. ISWAVE == 4 )then
    ! *** ONLY NEED TO ASSIGN INITIALLY
    !$OMP SINGLE
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          AH(L,K)  = AHOXY(L)
          AHC(L,K) = AHOXY(L)
        enddo
      enddo
    enddo
    !$OMP END SINGLE
  endif

  !----------------------------------------------------------------------!
  ! ***  CALCULATE HORIZONTAL DIFFUSION DUE TO WAVE BREAKING
  if( ISWAVE == 2 .or. ISWAVE == 4 )then
    if( WVLSH > 0.0 .or. WVLSX > 0.0 )then
      !$OMP DO PRIVATE(ND,K,LP,L,LN,LS,LW,LSW,LSE,DTMPH,DTMPX,DTMP,TMPVAL,HPLW,HPLS,HPLSW) 
      do ND = 1,NDM
        do K = 1,KC
          do LP = 1,LLHDMF(K,ND)
            L = LKHDMF(LP,K,ND)  
            DTMPH = (WVFACT*WV(L).DISSIPA(K))**0.3333
            DTMPX = WV(L).DISSIPA(K)/HP(L)
            AH(L,K) = AH(L,K) + WVLSH*DTMPH*HP(L) + AHWVX*DTMPX
          enddo
        enddo

        do K = 1,KC
          do LP = 1,LLHDMF(K,ND)
            L = LKHDMF(LP,K,ND)  
            LS = LSC(L)
            LW = LWC(L) 
            LSW = LSWC(L)
            
            TMPVAL = 1.+SUB(L)+SVB(L)+SUB(L)*SVB(L)
            HPLW   = SUB(L)*HP(LW)
            HPLS   = SVB(L)*HP(LS)
            HPLSW  = SUB(L)*HP(LSW)
            HMC(L)    = (HP(L)+HPLW+HPLS+HPLSW)/TMPVAL

            DTMP = 0.25*( WV(L).DISSIPA(K) + WV(LW).DISSIPA(K) + WV(LS).DISSIPA(K) + WV(LSW).DISSIPA(K) )

            DTMPH = (WVFACT*DTMP)**0.3333
            DTMPX = DTMP/HMC(L)

            AHC(L,K) = AHC(L,K) + WVLSH*DTMPH*HMC(L) + AHWVX*DTMPX
          enddo
        enddo
      enddo  ! *** END OF DOMAIN  
      !$OMP END DO
    endif
  endif

  ! *** SW CORNER AVERAGE DEPTHS (H1C)
  !$OMP DO PRIVATE(ND,LP,L,LS,LW,LSW,TMPVAL) 
  do ND = 1,NDM
    do LP = 1,LLHDMF(KC,ND)
      L   = LKHDMF(LP,KC,ND)  
      LS  = LSC(L)
      LW  = LWC(L)
      LSW = LSWC(L)

      H1C(L) = H1P(L)
      H1C(L) = H1C(L) + SUB(L)*H1P(LW)   ! *** WEST
      H1C(L) = H1C(L) + SVB(L)*H1P(LS)   ! *** SOUTH
      H1C(L) = H1C(L) + SUB(L)*H1P(LSW)  ! *** SOUTHWEST

      TMPVAL = 1.0/( 1.0 + SUB(L) + SUB(L) + 0.5*(SUB(L)+SVB(L)) )
      H1C(L) = TMPVAL*H1C(L)
    enddo
  enddo  ! *** END OF DOMAIN  
  !$OMP END DO

  !----------------------------------------------------------------------!
  ! ***  CALCULATE DIFFUSIVE MOMENTUM FLUXES
  !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
  do ND = 1,NDM
    do K = 1,KC
      do LP = 1,LLHDMF(K,ND)
        L  = LKHDMF(LP,K,ND)  
        LS = LSC(L)
        LW = LWC(L)

        FMDUX0(L,K) = 2.0* DYP(L)*H1P(L)*AH(L,K)*DXU1(L,K)
        FMDUY0(L,K) = 0.5*( DXU(L) + DXU(LS) )*H1C(L)*AHC(L,K)*( DYU1(L,K) + DXV1(L,K) )
        FMDVY0(L,K) = 2.0* DXP(L)*H1P(L)*AH(L,K)*DYV1(L,K)
        FMDVX0(L,K) = 0.5*( DYV(L) + DYV(LW) )*H1C(L)*AHC(L,K)*( DYU1(L,K) + DXV1(L,K) )
      enddo
    enddo
  enddo  ! *** END OF DOMAIN  
  !$OMP END DO
  
  ! *** Advance the horizontal diffusion termsy
  if( ISHDMFILTER > 0 .and. ((ISRESTI > 0 .and. NITER > 1) .or. NITER > 100) )then
    ! *** Apply a 25% growth rate limit on rising/falling HMD terms
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
    !$OMP DOPRIVATE(ND,K,LP,L)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          if( LASPECT(L) )then
            FMDUY(L,K) = 0.0
            FMDVX(L,K) = 0.0
          endif
        enddo
      enddo
    enddo  
    !$OMP END DO
  endif
  !$OMP END PARALLEL

  ! *** ZERO BOUNDARY CELL MOMENTUM DIFFUSION
  do LP = 1,NBCSOP
    L = LOBCS(LP)
    do K = 1,KC
      AHC(L,K) = 0.0
      FMDUX(L,K) = 0.0
      FMDUY(L,K) = 0.0
      FMDVY(L,K) = 0.0
      FMDVX(L,K) = 0.0
    enddo
  enddo

  return
  
END SUBROUTINE
