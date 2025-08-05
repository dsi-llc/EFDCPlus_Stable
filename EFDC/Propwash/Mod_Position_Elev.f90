! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
module Mod_Position_Elev
    
    use GLOBAL, only : RKD, RK4
    use Mod_Position
    
    implicit none
    
    private
    public :: position_elev
    
    type, extends(position) :: position_elev
        
        real (kind = RKD)   :: b_elev  = 0.0       !< bottom elevation for the positions x,y coordinate[meters]
        real (kind = RKD)   :: w_depth = 0.0       !< water depths for the x,y coorinda [meters]
        
    end type position_elev
    
end module Mod_Position_Elev