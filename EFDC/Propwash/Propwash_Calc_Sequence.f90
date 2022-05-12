! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
!> @details Calls subroutines to determine the  effect on the sediment
!! bed from the propwash
!> @author  Zander Mausolff & Paul Craig
!---------------------------------------------------------------------------!
Subroutine Propwash_Calc_Sequence(ieffluxonly)

  Use GLOBAL, only : timeday
  Use Variables_Propwash

  Use Mod_Active_Ship
  Use Variables_Ship

  implicit none

  ! *** Dummy variables
  integer, intent(in) :: ieffluxonly                           !< Used to bypass field and erosion calculations if only need efflux velocities for momentum field
  
  ! *** Local variables
  integer :: i, l
  integer :: prev_id, next_id !< keeps track of the array indices for the relevant ship tracks
  real(kind = rkd) :: start_time, end_time, effvel, ttds
  integer :: final_track_pos_id
  type(track_ids), allocatable, dimension(:), save :: pos_ids !> keeps track of which ship track we are in and which position within that track
  logical :: in_track, local_debug

  integer :: track_id
    
  REAL(RKD), EXTERNAL :: DSTIME

  ! *** One time processes
  if( .not. allocated(pos_ids) )then
    allocate(pos_ids(total_ships))       ! *** Note, all values in ids initialized to 1

    ! *** Determine the adjacent cells for each given valid cell
    Call Det_Adjacent_Cells
  endif
  
  if( (istran(6) + istran(7)) > 0 )then
    DTSEDJ = REAL(DTSED,8)
  else
    DTSEDJ = REAL(DELT,8)
  endif

#ifdef DEBUGGING
  local_debug = .true.
#else
  local_debug = .false.
#endif

  TTDS = DSTIME(0)             ! *** Track computational time used by propwash calculations
  freq_out_min  = 100000.
  iwrite_pwmesh = 0
  
  ! *** Perform the propwash related calculations over all the boats
  nactiveships = 0
  do i = 1, total_ships
    in_track = .FALSE.
    track_id = 1            ! delme - add this to the ship structure

    ! *** Update the time of the ship to the current time
    all_ships(i).pos.time = timeday
    all_ships(i).efflux_vel = 0.0

    ! *** see if the current time falls within any track
    Call all_ships(i).det_if_in_track(track_id, in_track)

    ! *** if the current time falls within a prescribed track, run through the propwash sequence
    if( in_track )then
      nactiveships = nactiveships + 1
      
      ! *** Determine which position within the track we are in
      Call all_ships(i).det_pos_in_track(track_id, pos_ids(i).prev, pos_ids(i).next)

      ! *** Interpolate position between track points
      !> @todo Set some criteria to avoid interpolation if the ship has not moved much or at all
      Call all_ships(i).interp_track(track_id, pos_ids(i).prev, pos_ids(i).next)

      if( all_ships(i).power == 0.0 .or. all_ships(i).pos.cell < 2 ) cycle   ! *** Skip if ship is outside the domain, docked or anchored
      
      if( ieffluxonly == 0 )then
        ! *** Setup 2d mesh for the velocity field
        call all_ships(i).setup_mesh(local_debug)

        ! *** Determine the erosive flux from the propwash velocity field
        call all_ships(i).calc_erosive_flux(local_debug)
      else
        ! *** Get the efflux velocity
        call all_ships(i).calc_velocity(10._RKD, 10._RKD, effvel, ieffluxonly)
      endif
    else
      ! *** Reset track counters
      pos_ids(i).prev = 1
      pos_ids(i).next = 2
    end if 
    
    ! *** Determine minimum output frequency
    if( all_ships(i).freq_out > 0. )then
      freq_out_min = min(all_ships(i).freq_out, freq_out_min)
    endif
    
  end do

  ! *** Sum to total erosion (g)
  if( nactiveships > 0 .and. ieffluxonly == 0 .and. (istran(6) > 0 .or. istran(7) > 0) )then
    do l=2,la
      prop_ero(l,0) = sum(prop_ero(l,1:nscm))    
    enddo
    if( icalc_bl > 0 )then
      do l=2,la
        prop_ero(l,0) = prop_ero(l,0) + sum(prop_bld(l,1:nscm))     ! ***  PROP_ERO(L,0) is used as a flag whether there is any 
        prop_bld(l,0) = sum(prop_bld(l,1:nscm))                     ! ***                erosion due to propwash in a cell.   
      enddo
    endif
  endif

  TPROPW = TPROPW + (DSTIME(0)-TTDS) 

End subroutine Propwash_Calc_Sequence

