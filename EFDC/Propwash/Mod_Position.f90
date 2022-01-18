! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Mod_Position
    
    Use GLOBAL, only : RKD, RK4
    Use Mod_MPI_Helper_Functions
    
    implicit none
    
    private
    public :: position
    
    type :: position
        Real (kind = RKD)   :: time   = 0.0       !< Current time we are at
        Real (kind = RKD)   :: x_pos  = 0.0       !< x position [meters]
        Real (kind = RKD)   :: y_pos  = 0.0       !< y position [meters]
        Real (kind = RKD)   :: z_pos  = 0.0       !< z position [meters]
        Real (kind = RKD)   :: var(3) = 0.0       !< Spatially dependent variable
        Real (kind = RKD), allocatable :: ero(:)  !< Spatially dependent erosion rate by sediment class
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
        
        if(present(counter) )then
            write(unit_num,'(a,i8, 4f15.5,I8)') 'i | (x,y,z,v,L) ',counter, self.x_pos, self.y_pos, self.z_pos, self.var(1), L
        else
            write(unit_num,'(4f15.5,I8)') self.x_pos, self.y_pos, self.z_pos, self.var(1), L
        end if
        
    end subroutine write_out
    
End module Mod_Position