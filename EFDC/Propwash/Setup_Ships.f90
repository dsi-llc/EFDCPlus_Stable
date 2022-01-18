! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
module Mod_Setup_Ships

contains

!---------------------------------------------------------------------------!
!> @details Sets up the ship objects using the read in prop data
!> @author Zander Mausolff
!---------------------------------------------------------------------------!
subroutine Setup_Ships(test_on)

  Use GLOBAL, only : timeday
  Use Variables_Propwash
  Use Variables_Ship

  Use Mod_Active_Ship
  Use Mod_All_Tracks
  Use Mod_Read_Propwash

  implicit none

  ! *** Dummy variables
  logical, optional, intent(in) :: test_on

  ! *** Local variables
  integer :: i, j, k, m
  integer :: total_tracks
  type(position_cell) :: start_position
  type(position_cell) :: new_position
  type(ship_type)     :: inp_ship
  type(active_ship)   :: new_ship
  real(rkd)           :: diam, dtime, dxx, dyy, distance, sternx, sterny

  allocate(all_ships(total_ships))

  do i = 1, total_ships
    if( all_read_tracks(i).num_tracks < 1 ) cycle
    
    ! *** Get the number of tracks for a particular ship
    if(present(test_on) )then
      ! *** The starting position of the ship should be the first track location
      start_position = position_cell (all_read_tracks(i).all_ship_tracks(1).track_pos(1).x_pos, &
                                      all_read_tracks(i).all_ship_tracks(1).track_pos(1).y_pos, &
                                      all_read_tracks(i).all_ship_tracks(1).track_pos(1).time,  &
                                      test_on)
    else
      ! *** The starting position of the ship should be the first track location
      start_position = position_cell (all_read_tracks(i).all_ship_tracks(1).track_pos(1).x_pos, &
                                      all_read_tracks(i).all_ship_tracks(1).track_pos(1).y_pos, &
                                      all_read_tracks(i).all_ship_tracks(1).track_pos(1).time )
    end if

    ! *** get the propeller object from the list of all ships
    inp_ship = all_new_ships(i)

    ! *** Instantiate new ship object
    new_ship = active_ship(inp_ship.mmsi,                      &
                           inp_ship.name,                      &
                           start_position,                     &
                           inp_ship,                           &
                           all_read_tracks(i).all_ship_tracks, &
                           all_read_tracks(i).num_tracks)

    ! *** Loop over the number of tracks for the ship
    do m = 1,new_ship.num_tracks
      ! *** Loop over the number of positions in a track
      do k = 1, new_ship.tracks(m).num_positions
        new_position = new_ship.tracks(m).track_pos(k)
        
        if( new_position.speed == -999. .and. k > 1 )then
          ! *** Compute speed from track points and time
          dtime = (new_position.time - new_ship.tracks(m).track_pos(k-1).time)*86400.   ! *** Seconds
          dxx   = (new_position.x_pos - new_ship.tracks(m).track_pos(k-1).x_pos)        ! *** meters
          dyy   = (new_position.y_pos - new_ship.tracks(m).track_pos(k-1).y_pos)        ! *** meters
          distance = sqrt(dxx**2 + dyy**2)
          dtime = max(dtime,1.)
          new_position.speed = min(distance/dtime,20.)                                  ! *** meters/second
        endif
        if( k == 2 )then
          if( new_ship.tracks(m).track_pos(k-1).speed == -999. ) new_ship.tracks(m).track_pos(k-1).speed = new_position.speed
        endif

        if( new_position.course == -999. .and. k > 1 )then
          ! *** Compute course from track points 
          dxx   = (new_position.x_pos - new_ship.tracks(m).track_pos(k-1).x_pos)        ! *** meters
          dyy   = (new_position.y_pos - new_ship.tracks(m).track_pos(k-1).y_pos)        ! *** meters
          new_position.course = atan2(dxx,dyy)                                          ! *** Compass angle towards in radians
        elseif( new_position.course /= -999 )then
          new_position.course = new_position.course / 180. * PI                         ! *** Convert to radians
        endif
        if( k == 2 )then
          if( new_ship.tracks(m).track_pos(k-1).course == -999. ) new_ship.tracks(m).track_pos(k-1).course = new_position.course
        endif
        
        if( new_position.heading == -999 .and. k > 1 )then
          ! *** Compute heading from track points 
          dxx   = (new_position.x_pos - new_ship.tracks(m).track_pos(k-1).x_pos)        ! *** meters
          dyy   = (new_position.y_pos - new_ship.tracks(m).track_pos(k-1).y_pos)        ! *** meters
          new_position.heading = atan2(dxx,dyy)                                         ! *** Compass angle towards in radians
        elseif( new_position.heading /= -999 )then
          new_position.heading = new_position.heading / 180. * PI                       ! *** Convert to radians
        endif
        if( k == 2 )then
          if( new_ship.tracks(m).track_pos(k-1).heading == -999. ) new_ship.tracks(m).track_pos(k-1).heading = new_position.heading
        endif
        
        if( new_position.draft == -999 .and. k > 1 )then
          ! *** Compute draft from ship data 
          new_position.draft = new_ship.ship.prop_diam + 1.0                            ! *** Keep propeller submerged
        endif
        if( k == 2 )then
          if( new_ship.tracks(m).track_pos(k-1).draft == -999. ) new_ship.tracks(m).track_pos(k-1).draft = new_position.draft
        endif
        
        if( new_position.power == -999 .and. k > 1 )then
          ! *** Test if RPS provided instead
          if( new_position.rps > 0.0 )then
              if( new_ship.ship.max_rps > 1. )then
                new_ship.ship.max_rps = new_position.rps
              endif
              new_position.power = -new_position.rps/new_ship.ship.max_rps
          else
            ! *** Missing power and rps
            if( new_position.speed < 1. )then
              new_position.power = 0.0
              new_position.rps   = 0.0
            else
              ! *** Compute applied power from ship data  delme !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delme
              new_position.power = -0.35       ! *** For now, just use 35% of max RPS if underway
              new_position.rps = new_position.power*new_ship.ship.max_rps
            endif
            ! delme
            ! delme
            ! delme
            ! delme
          endif
        endif
        if( k == 2 )then
          if( new_ship.tracks(m).track_pos(k-1).power == -999. ) new_ship.tracks(m).track_pos(k-1).power = new_position.power
          if( new_ship.tracks(m).track_pos(k-1).rps == -999. )   new_ship.tracks(m).track_pos(k-1).rps   = new_position.rps
        endif
                
        ! *** Apply offsets and conversions
        if( new_position.heading /= -999. )then
          sternx = new_position.x_pos - sin(new_position.heading)*new_ship.ship.ais_to_stern         ! *** Adjust for AIS antenna position relative to stern
          sterny = new_position.y_pos - cos(new_position.heading)*new_ship.ship.ais_to_stern         ! *** Adjust for AIS antenna position relative to stern

          new_position.x_pos = sternx + sin(new_position.heading)*new_ship.ship.dist_from_stern      ! *** Adjust from stern to center of propeller location
          new_position.y_pos = sterny + cos(new_position.heading)*new_ship.ship.dist_from_stern      ! *** Adjust from stern to center of propeller location

          new_position.heading = new_position.heading - pi/2.                                        ! *** Reverse direction from "to" to "from" (use compass orientation)
        endif

        new_position.draft = new_position.draft + new_ship.ship.prop_offset                          ! *** Depth of the propeller shaft (meters)
        new_position.draft = max(new_position.draft, 0.5*(0.0 + new_ship.ship.prop_diam))
  
        new_ship.tracks(m).track_pos(k) = new_position                                               ! *** Position now reflects propeller position
      enddo  ! *** End of number of points in track
    enddo    ! *** End of number of tracks
    
    ! *** Derived contant parameters
    new_ship.prop_radius = new_ship.ship.prop_diam / 2.
    if( new_ship.ship.ducted > 0 )then
      diam = new_ship.ship.prop_diam
    else
      diam = new_ship.ship.prop_diam*0.71                                                            ! *** Account for the jet contraction of an unducted propeller
    endif
    new_ship.prop_area   = PI*(0.5*diam)**2                                                          ! *** Single propeller area of effective diameter
    new_ship.prop_area   = new_ship.prop_area * new_ship.ship.num_props                              ! *** Total propeller area
    
    Allocate(new_ship.mesh_count(0:LCM))
    new_ship.mesh_count = 0
    
    ! *** Add to array of all ships
    all_ships(i) = new_ship

  end do


  ! *** No longer need data containing all ship tracks or all prop info
  deallocate(all_read_tracks)
  deallocate(all_new_ships)

  end subroutine Setup_Ships

end module Mod_Setup_Ships
