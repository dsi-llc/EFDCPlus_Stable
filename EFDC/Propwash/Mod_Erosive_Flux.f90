! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Mod_Erosive_Flux
    
    Use GLOBAL, only : RKD, RK4
    
    private
    public :: erosive_flux
    
    type :: erosive_flux
    
        Real(kind = RKD), Allocatable, Dimension(:) :: sub_grid_flux 
        Real(kind = RKD) :: total_flux = 0.0
        
    end type erosive_flux
    
    ! *** interface for the constructor
    interface erosive_flux
        module procedure :: constructor_erosive_flux
    end interface erosive_flux
    
    contains
    
    ! *** construction definition
    type(erosive_flux) function constructor_erosive_flux(nx) result(self)
        ! *** Dummy variables
        integer :: nx
        
        ! *** Local variables
        
        ! *** setup the class
        allocate(self.sub_grid_flux (nx))
        self.sub_grid_flux = 0.0
        
        !self.total_flux = 0.0
        
    end function constructor_erosive_flux
    
End module Mod_Erosive_Flux