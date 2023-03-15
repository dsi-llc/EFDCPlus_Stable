! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Variables_Ship
    
  Use Mod_all_tracks
  Use Mod_Active_Ship
  use Mod_all_ship_tracks
    
  implicit none
    
  type :: track_ids
    integer :: prev = 1
    integer :: next = 1
  end type
  
  type(ship_type),       allocatable, dimension(:) :: all_new_ships
  type(active_ship),     allocatable, dimension(:) :: all_ships
  type(all_ship_tracks), allocatable, dimension(:) :: all_read_tracks

End module Variables_Ship
