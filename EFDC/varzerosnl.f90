! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE VARZEROSNL

  !C *** THIS SUBROUTINE ZERO'S MANY ARRAYS AFTER ALLOCATION
  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 

  USE GLOBAL
  
  IMPLICIT NONE

  ! Begin SEDZLJ integer variables
  NCORENO = 0

  ! *** SCALARS
  HPMIN=0.25      ! *** Minimum depth to compute shears
  MAXDEPLIMIT=0.0 ! *** The maximum fraction of mass from bottom water column layer that can be deposited on active bed layer.
  RHO=1000.0      ! *** Density of water in kg/m^3.
  TACTM=0.0       ! *** Active layer multiplier
  WATERDENS=1.0   ! *** Density of water in g/cm^3
  
END SUBROUTINE VARZEROSNL
