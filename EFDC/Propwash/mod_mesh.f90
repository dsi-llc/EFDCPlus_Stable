! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
!> @details Contains subroutines to read propwash files
!> @author  Zander Mausolff
!---------------------------------------------------------------------------!
module mod_mesh
    
  use GLOBAL, only : RKD, RK4
  use mod_position
    
  implicit none

  private
    
  public :: mesh
    
  type, extends(position) :: mesh
    real(kind = rkd) :: area = 0.0    !< area of the cell in [m]
    real(kind = rkd) :: width = 0.0   !< width (radial direction) of the cell of the cell in [m]
    real(kind = rkd) :: length = 0.0  !< length (axial direction) of the cell of the cell in [m]
    integer          :: cell = 0      !< L index of the current position
  end type
    
    
end module mod_mesh