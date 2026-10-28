! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!< @details contains array of positions that makes up the 
!< @author Zander Mausolff
!---------------------------------------------------------------------------!
module mod_all_tracks
    
    use Mod_Position
    use Mod_Position_Cell
    use GLOBAL, only : RKD, RK4

    implicit none
    
    type all_tracks
    
        type(position_cell), allocatable, dimension(:) :: track_pos !< Array of all track positions for the ship
        integer :: num_positions = 0 !> keeps track of the number of positions the current track
        
    contains 

        procedure, pass(self) :: get_start_time
        procedure, pass(self) :: get_end_time
        procedure, pass(self) :: del_all_tracks

    end type

    contains
    !---------------------------------------------------------------------------!
    !< @details Determines the final time for a given track
    !< @author Zander Mausolff
    !---------------------------------------------------------------------------!
    subroutine get_start_time(self, start_time)
    
        implicit none
        
        ! *** Dummy variables
        class(all_tracks), intent(in) :: self
        real(kind = RKD), intent(inout) :: start_time
    
        start_time = self.track_pos(1).time
        
    end subroutine get_start_time
    
    !---------------------------------------------------------------------------!
    !< @details Determines the final time for a given track
    !< @author Zander Mausolff
    !---------------------------------------------------------------------------!
    subroutine get_end_time(self, final_time)
        
        implicit none
        
        ! *** Dummy variables
        class(all_tracks), intent(in) :: self
        real(kind = RKD), intent(inout) :: final_time !> final time for a given track
        
        
        final_time = self.track_pos( self.num_positions ).time
        
    
    end subroutine get_end_time
    
    !---------------------------------------------------------------------------!
    !< @details deallocates the all tracks array
    !< @author Zander Mausolff
    !---------------------------------------------------------------------------!
    subroutine del_all_tracks(self)
    
        implicit none
        ! ***
        class(all_tracks), intent(inout) :: self
        
        deallocate( self.track_pos)
        
    end subroutine del_all_tracks
    
end module 