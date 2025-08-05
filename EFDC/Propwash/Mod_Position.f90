! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Mod_Position
    
    use GLOBAL, only : RKD, RK4
    use Mod_MPI_Helper_Functions
    
    implicit none
    
    private
    public :: position
    
    type :: position
        real (kind = RKD)   :: time   = 0.0       !< Current time we are at
        real (kind = RKD)   :: x_pos  = 0.0       !< x position [meters]
        real (kind = RKD)   :: y_pos  = 0.0       !< y position [meters]
        real (kind = RKD)   :: z_pos  = 0.0       !< z position [meters]
        real (kind = RKD)   :: var(3) = 0.0       !< Spatially dependent variable
        real (kind = RKD), allocatable :: ero(:)  !< Spatially dependent erosion rate by sediment class
    contains
    
        procedure, pass(self) :: write_out
        procedure, pass(self) :: interp_pos
        
    end type position
    
    
    contains
    
    !---------------------------------------------------------------------------!
    !< @details Interpolate between two coordinates
    !< @param[inout] self - return the interpolated value 
    !---------------------------------------------------------------------------!
    subroutine interp_pos(self)
    
        implicit none
        ! *** dummy variables
        class(position), intent(inout) :: self
        ! *** local variables
        
    end subroutine interp_pos
    
    !---------------------------------------------------------------------------!
    !< @details write out the x,y,z positions and dependent variable
    !< @param[in] self
    !< @param[in] unit_num - unit number to write to
    !< @param[in](optional) - counter if we are writing out multiple x,y,z positions
    !---------------------------------------------------------------------------!
    subroutine write_out(self, unit_num, cell, counter)
        
        implicit none
        
        class(position), intent(in) :: self
        integer, intent(in) :: unit_num
        integer, intent(in), optional :: cell, counter
        integer :: L
        
        if(present(cell) )then
          L = cell
        else
          L = -1
        endif
        
        write(unit_num,'(a,4f15.5,I8)') 'i | (x,y,z,v,L) ', self.x_pos, self.y_pos, self.z_pos, self.var(1), L
        
    end subroutine write_out
    
End module Mod_Position
