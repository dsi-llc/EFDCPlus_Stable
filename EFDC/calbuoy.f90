  ! ----------------------------------------------------------------------
  !   This file is a part of EFDC+
  !   Website:  https://eemodelingsystem.com/
  !   Repository: https://github.com/dsi-llc/EFDC_Plus.git
  ! ----------------------------------------------------------------------
  ! Copyright 2021-2024 DSI, LLC
  ! Distributed under the GNU GPLv2 License.
  ! ----------------------------------------------------------------------
  SUBROUTINE CALBUOY(UPDATE)

  ! ***********************************************************************!
  ! *** CALBUOY CALCULATES THE BUOYANCY USING MELLOR'S APPROXIMATION
  ! *** TO THE UNESCO EQUATION OF STATE
  ! *** MELLOR, G.L., J. ATM AND OCEAN TECH, VOL 8, P 609
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! ***  2015-06     PAUL M. CRAIG      IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  ! ***  2014-09     PAUL M. CRAIG      ADDED THE LWET BYPASS APPROACH
  ! ***  2014-08     D H CHUNG          SET EXPLICIT PRECISIONS OF INTEGER & REAL
  ! ***  2011-03     Paul M. Craig      Converted to F90, added OMP
  !

  use GLOBAL
  use Allocate_Initialize
  use GSW

  implicit none

  logical,intent(IN) :: UPDATE
  integer            :: NS, K, L, ND, LP, NN
  real(RKD), save    :: RHOO, ONED, RHO1
  real(RKD)          :: SSTMP, TTMP, RHTMP, TEM0
  real(RKD)          :: PSW        !< Pressure (dbar or m)
  real(RKD)          :: ps, sa, tm
  real(RKD), parameter :: p_ref = 0. ! *** dbar

  real(RKD),save,allocatable,dimension(:,:) :: PTEM   !< Potential Temperature   (degC)

  if( (ISGOTM > 0 .or. ISTRAN(2) > 0 ) .and. .not. allocated(PTEM) )then
    call AllocateDSI( PTEM,LCM, KCM, 0.0)
  endif

  ! *************************************************************************************
  if( IBSC == 1 )then
    ! *** DENSITY AS A LINEAR FUNCTION OF SALINITY ONLY.  FOR DIAGNOSTIC PURPOSES ONLY
    do K = 1,KC
      do L = 2,LA
        B(L,K) = 0.00075_8*SAL(L,K)
      enddo
    enddo
    return
  endif

  ! *************************************************************************************
  ! *** Density RHOO AT P = 0, S = 0, and T = TEMO.  Only compute once, since TEMO is a constant
  if( N <= 5 )then
    ONED = 1._8
    TEM0 = ABS(TEMO)
    RHOO = 999.842594 + 6.793952D-2*TEM0 - 9.095290D-3*TEM0*TEM0 + 1.001685D-4*TEM0*TEM0*TEM0 - 1.120083D-6*TEM0*TEM0*TEM0*TEM0 + 6.536332D-9*TEM0*TEM0*TEM0*TEM0*TEM0
  endif

  !$OMP PARALLEL DEFAULT(SHARED)

  ! *************************************************************************************
  ! *** SAVE THE CURRENT BUOYANCY BEFORE UPDATING
  if( UPDATE )then
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          B1(L,K) = B(L,K)
        enddo
      enddo
    enddo
    !$OMP END DO
  endif

  ! *************************************************************************************
  ! *** Density calculations

  !$OMP DO PRIVATE(ND, K, L, LP, NN, NS, SSTMP, TTMP, RHTMP, TEM0, RHO1, PSW, tm, ps, sa)
  do ND = 1,NDM

    ! *** Get temperature needed for density and buoyancy calculations
    if( ISTRAN(2) > 0 )then
      if( ISGOTM > 0 .or. IBSC == 2 )then
        ! *** Compute potential temperature from temperature and pressure
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            PSW = (1. - ZZ(L,K))*HP(L)        ! *** Pressue at the mid point of the layer [m]
            tm  = TEM(L,K)                    ! *** In-situ temperature [deg C]
            sa  = max(SAL(L,K),0.)         ! *** Absolute Salinity [g/kg] (Prevent negative value)
            PTEM(L,K) = gsw_pt_from_t (sa, tm, PSW, p_ref)
          enddo
        enddo
      else      
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L  = LKWET(LP,K,ND)
            PTEM(L,K) = TEM(L,K)            ! *** In-situ temperature [deg C]
          enddo
        enddo
      endif
    endif

    ! *** CASE: No salinity and no temperature
    if( ISTRAN(1) == 0 .and. ISTRAN(2) == 0 )then
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)

          B(L,K) = 0.0       ! *** Buoyancy [dimensionless]
        enddo
      enddo

      ! *** CASE: Salinity and no temperature
    elseif( ISTRAN(1) >= 1 .and. ISTRAN(2) == 0 )then
      TEM0 = ABS(TEMO)
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SSTMP = max(SAL(L,K),0.)

          RHO1 = RHOO + SSTMP*(0.824493 - 4.0899D-3*TEM0       &
            + 7.6438D-5*TEM0*TEM0                              &
            - 8.2467D-7*TEM0*TEM0*TEM0                         &
            + 5.3875D-9*TEM0*TEM0*TEM0*TEM0)                   &
            + SQRT(SSTMP)*SSTMP*(-5.72466D-3 + 1.0227D-4*TEM0  &
            - 1.6546D-6*TEM0*TEM0)                             &
            + 4.8314D-4*SSTMP*SSTMP

          RHOW(L,K) = RHO1                ! *** Density [kg/m^3]
          B(L,K) = (RHO1/RHOO)-1._8       ! *** Buoyancy [dimensionless]
        enddo
      enddo

      ! *** CASE: No salinity but with temperature
    elseif( ISTRAN(1) == 0 .and. ISTRAN(2) >= 1 )then
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          TTMP = PTEM(L,K)

          RHO1 = 999.842594 + 6.793952D-2*TTMP - 9.095290D-3*TTMP*TTMP  &
            + 1.001685D-4*TTMP*TTMP*TTMP                &
            - 1.120083D-6*TTMP*TTMP*TTMP*TTMP           &
            + 6.536332D-9*TTMP*TTMP*TTMP*TTMP*TTMP

          RHOW(L,K) = RHO1                ! *** Density [kg/m^3]
          B(L,K) = (RHO1/RHOO)-1._8       ! *** Buoyancy [dimensionless]
        enddo
      enddo

      ! *** CASE: Both salinity and temperature
    elseif( ISTRAN(1) >= 1 .and. ISTRAN(2) >= 1 )then
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SSTMP = max(SAL(L,K),0.)
          TTMP = PTEM(L,K)

          RHTMP = 999.842594 + 6.793952D-2*TTMP - 9.095290D-3*TTMP*TTMP  &
            + 1.001685D-4*TTMP*TTMP*TTMP                &
            - 1.120083D-6*TTMP*TTMP*TTMP*TTMP           &
            + 6.536332D-9*TTMP*TTMP*TTMP*TTMP*TTMP

          RHO1 = RHTMP + SSTMP*(0.824493 - 4.0899D-3*TTMP + 7.6438D-5*TTMP*TTMP   &
            - 8.2467D-7*TTMP*TTMP*TTMP                          &
            + 5.3875D-9*TTMP*TTMP*TTMP*TTMP)                    &
            + SQRT(SSTMP)*SSTMP*(-5.72466D-3 + 1.0227D-4*TTMP   &
            - 1.6546D-6*TTMP*TTMP)          &
            + 4.8314D-4*SSTMP*SSTMP

          RHOW(L,K) = RHO1                ! *** Density [kg/m^3]
          B(L,K) = (RHO1/RHOO)-1._8       ! *** Buoyancy [dimensionless]
        enddo
      enddo
    endif

    !------------------------------------------------------------------------------------
    ! *** Apply low sediment concentration correction to buoyancy
    if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then

      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          TVAR1S(L,K) = 0.
          TVAR1W(L,K) = 0.
        enddo
      enddo

      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED2
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              TVAR1S(L,K) = TVAR1S(L,K) + SDEN(NS)*SED(L,K,NS)
              TVAR1W(L,K) = TVAR1W(L,K) + (SSG(NS)-1.)*SDEN(NS)*SED(L,K,NS)
            enddo
          enddo
        enddo
      endif

      if( ISTRAN(7) >= 1 )then
        do NN = 1,NSND
          NS = NN+NSED
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              TVAR1S(L,K) = TVAR1S(L,K) + SDEN(NS)*SND(L,K,NN)
              TVAR1W(L,K) = TVAR1W(L,K) + (SSG(NS)-1.)*SDEN(NS)*SND(L,K,NN)
            enddo
          enddo
        enddo
      endif

      if( ISTRAN(1) == 0 .and. ISTRAN(2) == 0 )then
        ! *** Reset density
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            RHOW(L,K) = RHOO
          enddo
        enddo
      endif

      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          B(L,K) = B(L,K)*(1. - TVAR1S(L,K)) + TVAR1W(L,K)

          ! ***  Correction for sediment
          RHOW(L,K) = RHOW(L,K)*( 1. - TVAR1S(L,K) + TVAR1W(L,K) )
        enddo
      enddo

    endif
  enddo    ! *** End of domain loop
  !$OMP END DO
  !$OMP END PARALLEL

  ! *************************************************************************************
  return
  END
