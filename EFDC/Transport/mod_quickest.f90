! ==========================================================================================
! Module: MOD_QUICKEST
!
! Purpose: QUICKEST (Quadratic Upstream Interpolation for Convective Kinematics
!          with Estimated Streaming Terms) transport scheme for EFDC+
!
! References:
!  - Leonard, B.P. (1979). "A stable and accurate convective modelling procedure based
!    on quadratic upstream interpolation." Comp. Methods Appl. Mech. Eng., 19(1), 59-98.
!  - Leonard, B.P. (1991). "The ULTIMATE conservative difference scheme applied to
!    unsteady one-dimensional advection." Comp. Methods Appl. Mech. Eng., 88(1), 17-74.
!  - GOTM adv_center.f90 (P2_PDM limiter implementation)
!
! Author: Based on W2 documentation and GOTM implementation
! Date: 2025-12-19
! ==========================================================================================

MODULE MOD_QUICKEST
  USE GLOBAL
  IMPLICIT NONE

  CONTAINS

  ! ************************************************************************************
  ! FUNCTION: QUICKEST_FACE
  !
  ! Purpose: Computes face concentration using third-order QUICKEST scheme
  !          with ULTIMATE limiter for monotonicity
  !
  ! Inputs:
  !   Cu   - Upstream concentration (2 cells upstream of face)
  !   Cc   - Center concentration (1 cell upstream of face, interpolation center)
  !   Cd   - Downstream concentration (at or downstream of face)
  !   cfl  - Courant number at face (|U|*dt/dx, always positive)
  !
  ! Returns: Face concentration estimate
  !
  ! Algorithm: Based on Leonard (1979, 1991) and GOTM P2_PDM limiter
  ! ************************************************************************************
  REAL(RKD) FUNCTION QUICKEST_FACE(Cu, Cc, Cd, cfl)
    IMPLICIT NONE

    ! Arguments
    REAL(RKD), INTENT(IN) :: Cu    ! Upstream concentration
    REAL(RKD), INTENT(IN) :: Cc    ! Center concentration
    REAL(RKD), INTENT(IN) :: Cd    ! Downstream concentration
    REAL(RKD), INTENT(IN) :: cfl   ! Courant number (positive)

    ! Local variables
    REAL(RKD) :: deltaf    ! Concentration difference (downstream - center)
    REAL(RKD) :: deltafu   ! Concentration difference (center - upstream)
    REAL(RKD) :: ratio     ! Slope ratio r = deltafu / deltaf
    REAL(RKD) :: limiter   ! Flux limiter function (phi)
    REAL(RKD) :: x         ! QUICKEST parameter
    REAL(RKD) :: cfl_safe  ! Bounded CFL to prevent division by zero

    ! Small number to prevent division by zero
    REAL(RKD), PARAMETER :: EPSILON = 1.0E-30_RKD

    ! *** Compute concentration differences
    deltaf  = Cd - Cc
    deltafu = Cc - Cu

    ! *** Check for monotonicity (sign of gradients)
    !     If deltaf and deltafu have opposite signs or either is zero,
    !     fall back to upwind (limiter = 0)
    IF( deltaf * deltafu <= 0.0_RKD )THEN
      ! Local extremum or constant region - use first-order upwind
      QUICKEST_FACE = Cc
      RETURN
    ENDIF

    ! *** Compute slope ratio
    ratio = deltafu / (deltaf + EPSILON)

    ! *** Ensure CFL is in valid range (should be < 1 for stability)
    !     Bound CFL to prevent division by zero in limiter
    cfl_safe = MIN(MAX(cfl, EPSILON), 0.999_RKD)

    ! *** Compute ULTIMATE QUICKEST limiter (P2_PDM from GOTM)
    !     This is Leonard's (1991) ULTIMATE limiter applied to QUICKEST

    ! QUICKEST polynomial correction parameter
    x = (1.0_RKD - 2.0_RKD*cfl_safe) / 6.0_RKD

    ! Third-order QUICKEST limiter
    limiter = (0.5_RKD + x) + (0.5_RKD - x) * ratio

    ! Apply ULTIMATE monotonicity constraints
    ! Limiter must satisfy all three conditions to ensure TVD property:
    !   1. limiter <= 2*ratio/cfl        (prevents undershoot on upstream side)
    !   2. limiter <= 2/(1-cfl)          (prevents overshoot on downstream side)
    !   3. limiter as computed by QUICKEST formula
    limiter = MIN(2.0_RKD * ratio / (cfl_safe + EPSILON), limiter)
    limiter = MIN(limiter, 2.0_RKD / (1.0_RKD - cfl_safe))

    ! Ensure limiter is non-negative (positivity preserving)
    limiter = MAX(limiter, 0.0_RKD)

    ! *** Compute face concentration using limited QUICKEST scheme
    QUICKEST_FACE = Cc + 0.5_RKD * limiter * (1.0_RKD - cfl_safe) * deltaf

  END FUNCTION QUICKEST_FACE

  ! ************************************************************************************
  ! SUBROUTINE: QUICKEST_INIT
  !
  ! Purpose: Initialize QUICKEST module (called during model initialization)
  ! ************************************************************************************
  SUBROUTINE QUICKEST_INIT()
    IMPLICIT NONE

  END SUBROUTINE QUICKEST_INIT

END MODULE MOD_QUICKEST
