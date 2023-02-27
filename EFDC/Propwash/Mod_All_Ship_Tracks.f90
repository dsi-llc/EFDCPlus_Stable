! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
module Mod_All_Ship_Tracks

  Use Mod_All_Tracks

  implicit none

  type All_Ship_Tracks

    type(all_tracks), allocatable, dimension(:) :: All_Ship_Tracks
    integer :: num_tracks = 0      !< total tracks for a given ship

  end type


end module mod_all_ship_tracks