! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !< @details  Module for active ship (SELF) for track, velocity, shear and
  !            Erosion rates
  !< @author  Paul Craig & Zander Mausolff
  !---------------------------------------------------------------------------!
Module Mod_Active_Ship

Use GLOBAL, only : RKD, RK4, PI

Use Mod_Position
Use Mod_Position_Cell
Use Mod_Ship
Use Mod_All_Tracks
Use Mod_Mesh
Use Variables_Propwash
Use Mod_MPI_Helper_Functions

implicit none

type ::active_ship

  integer (kind = RK4) :: mmsi       = 0        !< MMSI identifier of the ship
  character(len=100)   :: name       = ' '      !< Name of the ship
  type(position_cell)  :: pos                   !< Current position of the ship's propeller (x,y,z) coordinates
  type(position)       :: stern                 !< Current position of the ship's stern (x,y,z) coordinates
  type(position)       :: antenna               !< Current position of the ship AIS anntena (x,y,z) coordinates
  type(ship_type)      :: ship                  !< Ship properties

  real(kind = rkd)     :: rps        = 0.0      !< Current rps of the ship
  real(kind = rkd)     :: power      = 0.0      !< Current power of the ship [horsepower]
  real(kind = rkd)     :: speed      = 0.0      !< Current speed of the ship [m/s]
  real(kind = rkd)     :: course     = 0.0      !< Current course of the ship [radians]
  real(kind = rkd)     :: heading    = 0.0      !< Current heading of the ship [radians]
  real(kind = rkd)     :: draft      = 0.0      !< Current draught of the boat [meters]
  real(kind = rkd)     :: efflux_vel = 0.0      !< Current propeller efflux velocity [m/s].  Used for momentum addition

  real(kind = rkd)     :: prop_radius = 0.0      !< Derived parameter - Propeller radius [m]
  real(kind = rkd)     :: prop_area   = 0.0      !< Derived parameter - Total propeller area [m^2]
  real(kind = rkd)     :: max_bot_vel = 0.5      !< Maximum bottom velocity of last calculation

  integer (kind = RK4) :: bowcell    = 0        !< Current L index of the bow
  integer (kind = RK4) :: nsnap      = 0        !< Number of mesh snapshots written
  real(kind = rkd)     :: freq_out   = 0.0      !< Current subgrid mesh output frequncy [minutes]
  real(kind = rkd)     :: timesnap   = -99999.  !< Subgrid mesh output time [days]

  integer              :: num_tracks = 0        !< Number of tracks for this ship
  integer              :: itrack     = 0        !< Current track
  integer              :: ipos       = 0        !< Current track point
  integer (kind = RK4), allocatable :: mesh_count(:)           !< Count of number of subgrid cells are in a model cell

  type(all_tracks), allocatable, dimension(:) :: tracks        !< Array of all track positions for the ship
  type(mesh), allocatable, dimension(:)       :: axial_mesh    !< mesh pts along the axial line from the prop
  type(mesh), allocatable, dimension(:,:)     :: subgrid_mesh

contains
  ! *** define procedures that operate on the ship class
  procedure, pass(self) :: det_if_in_track                   !< determines the next ship position if in the track
  procedure, pass(self) :: det_pos_in_track
  procedure, pass(self) :: interp_track
  procedure, pass(self) :: setup_mesh
  procedure, pass(self) :: interp_depth_elev
  procedure, pass(self) :: interp_z
  procedure, pass(self) :: interp_bot_vel
  procedure, pass(self) :: delete_active_ship
  procedure, pass(self) :: write_xyz
  procedure, pass(self) :: calc_velocity
  procedure, pass(self) :: calc_erosive_flux
  procedure, pass(self) :: create_linkage
  end type active_ship

  ! *** Define interface for ship constructor
  interface active_ship
  module procedure :: constructor_active_ship
  end interface active_ship

contains

!---------------------------------------------------------------------------!
!< @details Constructor for the ship class
!---------------------------------------------------------------------------!
type(active_ship) function constructor_active_ship(inp_id, inp_name, inp_pos, inp_prop, inp_tracks, inp_num_tracks) result(self)

  ! *** Dummy variables
  integer  , intent(in) :: inp_id
  character(len=*)      , intent(in) :: inp_name
  type(position_cell)   , intent(in) :: inp_pos
  type(ship_type)       , intent(in) :: inp_prop
  type(all_tracks), Allocatable, Dimension(:), intent(in) :: inp_tracks !< Array of all track positions for the ship
  integer, intent(in) :: inp_num_tracks                                 !< number of tracks for this ship
  integer :: i,j

  ! *** Initialize based on inputed values
  self.mmsi = inp_id
  self.name = inp_name
  self.pos  = inp_pos
  self.ship = inp_prop
  self.tracks = inp_tracks
  self.num_tracks = inp_num_tracks

  ! *** Initialize to zero and set size of the 2d mesh array
  self.course  = 0.0
  self.heading = 0.0
  self.draft = 0.0

  ! *** Set the size of the meshs for the velocity field calculation

  allocate(self.axial_mesh(num_axial_elems))

  allocate(self.subgrid_mesh(num_radial_elems, num_axial_elems))

  do i=1,num_axial_elems
    do j=1,num_radial_elems
      allocate(self.subgrid_mesh(j,i).ero(nscm))
      self.subgrid_mesh(j,i).ero(:) = 0.0
    enddo
  enddo

end function constructor_active_ship

!---------------------------------------------------------------------------!
!< @details Interpolate between track points and gets ship heading
!< @param[inout] self - return the current location of the ship
!< @param[in] prev_id
!< @param[in] next_id
!---------------------------------------------------------------------------!
subroutine interp_track(self, track_id, prev_pos_id, next_pos_id, debug)

  Use Variables_MPI, only : mpi_log_unit

  implicit none

  ! *** Dummy variables
  class(active_ship), intent(inout) :: self
  integer, intent(in)               :: track_id
  integer, intent(in)               :: prev_pos_id
  integer, intent(in)               :: next_pos_id
  logical, optional, intent(in)     :: debug

  ! *** Local variables
  real(kind=rkd) :: dx, dy              !< Delta x and delta y (meters)
  real(kind=rkd) :: delta_t             !< time difference between two track positions
  real(kind=rkd) :: curr_delta_t        !< time difference between two track positions
  real(kind=rkd) :: x_0, y_0, x_1, y_1  !< Current and previous locations (meters)
  real(kind=rkd) :: x_b, y_b            !< Current bow location (meters)
  real(kind=rkd) :: radius, ratio       !< Radius from axial line to final point and ratio of distance between track points
  real(kind=rkd) :: p1, p2              !< Temporary power variables
  logical        :: no_heading_available

  ! *** Determine delta t between tracks (days)
  delta_t = (self.tracks(track_id).track_pos(next_pos_id).time) - (self.tracks(track_id).track_pos(prev_pos_id).time)
  delta_t = max(delta_t,1.15741E-05)    !*** 1.15741E-05 = 1/86400, one sec

  ! ***
  curr_delta_t = self.pos.time - self.tracks(track_id).track_pos(prev_pos_id).time

  ! *** x_0, y_0 are the previous ship position
  x_0 = self.tracks(track_id).track_pos(prev_pos_id).x_pos
  y_0 = self.tracks(track_id).track_pos(prev_pos_id).y_pos

  x_1 = self.tracks(track_id).track_pos(next_pos_id).x_pos
  y_1 = self.tracks(track_id).track_pos(next_pos_id).y_pos

  ! *** x_1,y_1 are the next ship positions
  ! *** Get changex,y directions
  dx = x_1 - x_0
  dy = y_1 - y_0
  ratio = curr_delta_t/delta_t

  ! *** Interpolate locations
  self.pos.x_pos = dx*ratio + x_0
  self.pos.y_pos = dy*ratio + y_0

  ! *** Cell Index.  For now, all propellers are assumed to be in the same cell
  self.pos.cell = get_cell(self.pos.cell, self.pos.x_pos, self.pos.y_pos)

  self.draft   = interp_value(self.tracks(track_id).track_pos(prev_pos_id).draft, &
                              self.tracks(track_id).track_pos(next_pos_id).draft, ratio)

  self.power   = interp_value(self.tracks(track_id).track_pos(prev_pos_id).power, &
                              self.tracks(track_id).track_pos(next_pos_id).power, ratio)

  self.speed   = interp_value(self.tracks(track_id).track_pos(prev_pos_id).speed, &
                              self.tracks(track_id).track_pos(next_pos_id).speed, ratio)

  self.heading = interp_angle(self.tracks(track_id).track_pos(prev_pos_id).heading, &
                              self.tracks(track_id).track_pos(next_pos_id).heading, ratio)

  self.course  =  interp_angle(self.tracks(track_id).track_pos(prev_pos_id).course, &
                               self.tracks(track_id).track_pos(next_pos_id).course, ratio)

  ! *** Deceleration
  if( self.tracks(track_id).track_pos(prev_pos_id).power < -1  .or. self.tracks(track_id).track_pos(next_pos_id).power < -1 )then
    self.heading = self.heading + pi                         !< Deceleration - Reverse heading to reflect propeller reversing
    
    if( self.tracks(track_id).track_pos(prev_pos_id).power > -1  .and. self.tracks(track_id).track_pos(next_pos_id).power < -1 )then
      ! *** Next point is declerating but previous point was not
      p1 = self.tracks(track_id).track_pos(prev_pos_id).power
      p2 = abs(self.tracks(track_id).track_pos(next_pos_id).power)/100.
    elseif( self.tracks(track_id).track_pos(prev_pos_id).power < -1  .and. self.tracks(track_id).track_pos(next_pos_id).power > -1 )then
      ! *** Previous point was declerating but next point is not
      p1 = abs(self.tracks(track_id).track_pos(prev_pos_id).power)/100.
      p2 = self.tracks(track_id).track_pos(next_pos_id).power
    else
      p1 = abs(self.tracks(track_id).track_pos(prev_pos_id).power)/100.
      p2 = abs(self.tracks(track_id).track_pos(next_pos_id).power)/100.
    endif
    
    self.power   = interp_value(p1, p2, ratio)
  endif
  
  if( self.power < -0.001 )then
    self.rps   = -self.power*self.ship.max_rps               !< Fraction of max RPS
    self.power = -1.                                         !< Set to inactive
  elseif( self.power > 0.001 )then
    self.power = self.power*self.ship.max_power              !< Fraction of max HP
    self.rps   = -1.                                         !< Set to inactive
  else
    self.power =  0.                                         !< Power off
    self.rps   =  0.                                         !< Power off
  endif

  if( self.tracks(track_id).track_pos(prev_pos_id).freq_out > 0. )then
    self.freq_out = self.tracks(track_id).track_pos(prev_pos_id).freq_out
  else
    self.freq_out = self.ship.freq_out
  endif

  if( ispropwash == 2 )then
    ! *** Get bow cell
    x_b = self.pos.x_pos + self.ship.length*cos(self.heading)
    y_b = self.pos.y_pos + self.ship.length*sin(self.heading)
    self.bowcell = get_cell(self.bowcell, x_b, y_b)
  endif

  ! *** optional debug write statements
  if(present( debug) )then
    write(mpi_log_unit, '(a,f15.5)') 'dx       ', dx
    write(mpi_log_unit, '(a,f15.5)') 'dy       ', dy
    write(mpi_log_unit, '(a,f15.5)') 'delta t  ', delta_t
    write(mpi_log_unit, '(a,f15.5)') 'curr dt  ', curr_delta_t
    write(mpi_log_unit, '(a,f15.5)') 'heading  ', self.heading
    write(mpi_log_unit, '(a,f15.5)') 'speed    ', self.speed
  end if

end subroutine interp_track

real(kind=rkd) function interp_value(value1, value2, factor)

  implicit none
  real(kind=rkd), intent(in)  :: value1,value2, factor

  interp_value =  value1 + factor*(value2 - value1)
  end function interp_value

  real(kind=rkd) function interp_angle(angle1, angle2, factor)

  implicit none
  real(kind=rkd), intent(in)  :: angle1, angle2, factor
  real(kind=rkd) :: sina, cosa

  sina = interp_value(sin(angle1), sin(angle2), factor)
  cosa = interp_value(cos(angle1), cos(angle2), factor)
  interp_angle = atan2(sina, cosa)
end function interp_angle

!---------------------------------------------------------------------------!
!< @details determines if we area track based on the current time
!< @param[in] self
!---------------------------------------------------------------------------!
subroutine det_if_in_track(self, track_id, in_track)

  implicit none

  ! *** Dummy variables
  class(active_ship), intent(inout) :: self
  integer, intent(inout) :: track_id
  logical, intent(inout) :: in_track

  ! *** Local variables
  integer :: i, j
  real(kind = RKD) :: start_time, end_time
  real(kind = RKD) :: prev_time, next_time

  ! *** Determine if the current time liesthe time region for a given track
  do i = 1, self.num_tracks

    ! *** get the start and end time for a given track
    call self.tracks(i).get_start_time(start_time)

    call self.tracks(i).get_end_time(end_time)

    ! *** Check if we are within this track, if so leave the loop
    if( start_time < self.pos.time .and. end_time > self.pos.time )then
      track_id = i
      in_track = .TRUE.
      return
    else
      in_track = .FALSE.
    end if
  end do

end subroutine det_if_in_track

!---------------------------------------------------------------------------!
!< @details Determines which positions within a track we are in
!< @param[in] self
!---------------------------------------------------------------------------!
subroutine det_pos_in_track(self, track_id, prev_pos_id, next_pos_id)

  implicit none

  ! *** Dummy variables
  class(active_ship), intent(inout) :: self
  integer, intent(in)    :: track_id
  integer, intent(inout) :: prev_pos_id
  integer, intent(inout) :: next_pos_id
  ! *** Local variables
  integer :: j
  real(kind = RKD) :: prev_time, next_time

  self.ipos = 0
  self.itrack = 0

  ! *** loop over the number of positions in the track
  do j = prev_pos_id, self.tracks(track_id).num_positions - 1

    ! *** Check if we are between ship tracks
    prev_time = self.tracks(track_id).track_pos(j).time
    next_time = self.tracks(track_id).track_pos(j + 1).time

    if( prev_time < self.pos.time .and. next_time >= self.pos.time  )then
      ! *** set the track and position indices
      prev_pos_id = j
      next_pos_id = j + 1
      self.ipos = j
      self.itrack = track_id
      exit
    end if

  end do

end subroutine det_pos_in_track

!---------------------------------------------------------------------------!
!< @details create 2d radial and axial mesh, projects onto the efdc (x,y)
!! grid coordinates
!< @param[in] self
!< @param[in] debug
!---------------------------------------------------------------------------!
subroutine setup_mesh(self, debug)

  use Variables_Propwash
  use XYIJCONV, only : BLOCKED

  implicit none
  ! *** Dummy variables
  class(active_ship), intent(inout) :: self
  logical, intent(in), optional :: debug !< prints out mesh for debugging
  character(30) :: filename

  ! *** Local variables
  integer :: i, j, k, L, Llast, iFlag, ipmc
  
  real(kind = rkd) :: x_o, y_o
  ! *** local variables
  real (kind = RKD) :: radial_sizing(num_radial_elems)
  real (kind = RKD) :: tmp_axial(num_axial_elems)

  real (kind = RKD) :: lower_rad, upper_rad, range_rad
  real (kind = RKD) :: lower_ax, upper_ax, range_ax
  real (kind = RKD) :: area, width_propellers, tan13, cone_max
  real (kind = RKD) :: width, length
  integer :: last_nodes(num_radial_elems), first_axial(num_radial_elems)

  ! *** Ignore track points outside domain
  if( self.pos.cell < 2 ) return

  width_propellers = real(self.ship.num_props-1)*self.ship.dist_between_props     ! *** Total width of all props between shafts

  lower_rad = -0.5*mesh_width*(self.ship.prop_diam + 0.5*width_propellers)        ! *** Width of mesh below ship centerline
  upper_rad =  0.5*mesh_width*(self.ship.prop_diam + 0.5*width_propellers)        ! *** Width of mesh above ship centerline
  range_rad = upper_rad - lower_rad

  lower_ax =  0.5*efflux_zone_mult*self.ship.prop_diam
  upper_ax =  mesh_length*self.ship.prop_diam
  range_ax = upper_ax - lower_ax

  width  = range_rad / num_radial_elems
  length = range_ax / num_axial_elems
  if( LSEDZLJ )then
    area = width*length*10000.                                                    ! *** Area in cm^2
  else                                                                            
    area = width*length                                                           ! *** Area in m^2
  endif
  width_propellers = 0.5*(self.ship.prop_diam + width_propellers)                 ! *** Half of the width of all props
  tan13 = 0.23087                                                                 ! *** Tan(13) velocity zone

  last_nodes = num_axial_elems
  first_axial = 0
  self.mesh_count = 0                                                             ! *** Zero subgrid cell count
  
  ! *** These loops are faster without OMP
  !print '(a,i12,a,f12.4,a,i6)', 'Mesh setup for MMSI = ', self.mmsi, ' @ ', timeday, ' L = ', self.pos.cell   ! delme
  
  ! *** Create x,y mesh
  do i = 1, num_axial_elems
    self.axial_mesh(i).x_pos = 0.0
    self.axial_mesh(i).y_pos = lower_ax +  range_ax * (i - 1) / (num_axial_elems - 1)
    do j = 1, num_radial_elems
      self.subgrid_mesh(j,i).x_pos = lower_rad + range_rad * (j - 1) / (num_radial_elems - 1)
      self.subgrid_mesh(j,i).y_pos = lower_ax +  range_ax * (i - 1) / (num_axial_elems - 1)
      ! *** constant area
      self.subgrid_mesh(j,i).area   = area
      self.subgrid_mesh(j,i).width  = width
      self.subgrid_mesh(j,i).length = length
      
      ! *** Check 13 degree cone
      cone_max = tan13 * self.subgrid_mesh(j,i).y_pos + width_propellers
      if( abs(self.subgrid_mesh(j,i).x_pos) > cone_max )then
        ! *** Deactivate points outside cone
        self.subgrid_mesh(j,i).x_pos = -9999.
        self.subgrid_mesh(j,i).x_pos = -9999.
        self.subgrid_mesh(j,i).area  = 0.
      endif
    end do
  end do

  ! *** Get index of first valid axial point for each radial
  do j = 1,num_radial_elems
    do i = 1,num_axial_elems
      if( abs(self.subgrid_mesh(j,i).x_pos + 9999.) > 1.0 )then
        first_axial(j) = i
        exit
      endif
    enddo
  enddo
  
  ! *** Translate and rotate points so the (x,y) coordinates align with the ship direction and location
  do i = 1, num_axial_elems

    x_o = self.axial_mesh(i).x_pos
    y_o = self.axial_mesh(i).y_pos

    self.axial_mesh(i).x_pos = x_o*sin(self.heading) - y_o*cos(self.heading) + self.pos.x_pos
    self.axial_mesh(i).y_pos = x_o*cos(self.heading) + y_o*sin(self.heading) + self.pos.y_pos

    ! *** Rotate the radial components
    do j = 1,num_radial_elems
      ! *** make some temporary variables to increase readability

      if( abs(self.subgrid_mesh(j,i).x_pos + 9999.) > 1.0 )then
        x_o = self.subgrid_mesh(j,i).x_pos
        y_o = self.subgrid_mesh(j,i).y_pos

        ! *** x rotation + translation about starting point
        self.subgrid_mesh(j,i).x_pos = x_o*sin(self.heading) - y_o*cos(self.heading) + self.pos.x_pos
        ! *** y rotation + translation about starting point
        self.subgrid_mesh(j,i).y_pos = x_o*cos(self.heading) + y_o*sin(self.heading) + self.pos.y_pos
      else
        self.subgrid_mesh(j,i).x_pos = -9999.
        self.subgrid_mesh(j,i).y_pos = -9999.
        self.subgrid_mesh(j,i).cell  = 0
      endif
    end do
  end do

  ! *** Assign L Index (OMP of LA loop in get_cell)
  do i = 1, num_axial_elems
    do j = 1,num_radial_elems
      if( abs(self.subgrid_mesh(j,i).x_pos + 9999.) < 1.0 ) cycle
      Llast = self.subgrid_mesh(j,i).cell
      self.subgrid_mesh(j,i).cell = get_cell(Llast, self.subgrid_mesh(j,i).x_pos, self.subgrid_mesh(j,i).y_pos)
      
      self.mesh_count(self.subgrid_mesh(j,i).cell) = self.mesh_count(self.subgrid_mesh(j,i).cell) + 1
    end do
  enddo

  ! *** 2021-04-27, NTL: Check blocking by dry land or cell masks
  do j = 1,num_radial_elems
    do i = 1, num_axial_elems - 1
      if( abs(self.subgrid_mesh(j,i).x_pos + 9999.) > 1.0 .and. abs(self.subgrid_mesh(j,i+1).x_pos + 9999.) > 1.0 ) then
        ! *** Check if the ray is blocked by dry land
        if( self.subgrid_mesh(j,i+1).cell == 0 )then
          last_nodes(j) = i
          exit
        endif
        if( self.subgrid_mesh(j,i).cell == 0 )then
          last_nodes(j) = i-1
          exit
        endif
        
        ! *** Check if the ray is blocked by any cell mask
        if( BLOCKED(self.subgrid_mesh(j,i).x_pos,   self.subgrid_mesh(j,i).y_pos,    &
                    self.subgrid_mesh(j,i+1).x_pos, self.subgrid_mesh(j,i+1).y_pos,  &
                    IL(self.subgrid_mesh(j,i).cell), JL(self.subgrid_mesh(j,i).cell)) )then
          last_nodes(j) = i
          exit
        endif
      endif
    enddo
  enddo
  
  ! *** deactivate the blocked grid nodes
  do j = 1,num_radial_elems
    if( last_nodes(j) < num_axial_elems )then
      do i = last_nodes(j)+1, num_axial_elems
        self.subgrid_mesh(j,i).cell  = 0
      enddo
    endif

    do i = 1, num_axial_elems
      if( self.subgrid_mesh(j,i).cell  == 0 )then
        self.subgrid_mesh(j,i).x_pos = -9999.
        self.subgrid_mesh(j,i).y_pos = -9999.
      endif
    enddo
  enddo

  ! *** Check for missing radials
  if( debug )then
    do j = 1,num_radial_elems
      i = first_axial(j)              ! *** First point in axial
      if( i < 1 .or. i > num_axial_elems )then
        write(901,'(a,i4,a,i10,f10.4, a,i12, a,2i5,2f12.2,f6.1)') 'No Valid Points for Radial: ', j, '  @ ', niter, timeday,   &
                                                                  '  for MMSI: ', self.mmsi,                                   &
                                                                  '  Track, Pt, X, Y Heading: ', self.itrack, self.ipos, self.pos.x_pos, self.pos.y_pos, self.heading
      endif
    enddo
  endif
  
  ! *** Final pass to remove axial points that start after the blockage
  ! *** Lower side
  do j = num_radial_elems/2-2,1,-1
    i = first_axial(j)              ! *** First point in axial
    if( i < 1 .or. i > num_axial_elems )then
      cycle
    endif
    if( abs(self.subgrid_mesh(j+1,i).x_pos + 9999.) < 1.0 )then
      do i = first_axial(j), num_axial_elems
        self.subgrid_mesh(j,i).x_pos = -9999.
        self.subgrid_mesh(j,i).y_pos = -9999.
        self.subgrid_mesh(j,i).cell  = 0
      enddo
    endif
  enddo
      
  ! *** Upper side
  do j = num_radial_elems/2+2,num_radial_elems
    i = first_axial(j)              ! *** First point in axial
    if( i < 1 .or. i > num_axial_elems )then
      cycle
    endif
    if( abs(self.subgrid_mesh(j-1,i).x_pos + 9999.) < 1.0  )then
      do i = first_axial(j), num_axial_elems
        self.subgrid_mesh(j,i).x_pos = -9999.
        self.subgrid_mesh(j,i).y_pos = -9999.
        self.subgrid_mesh(j,i).cell  = 0
      enddo
    endif
  enddo
  
  if( debug .and. (niter == 1 .or. mod(niter,100) == 0) )then
    do i = 1, num_axial_elems
      do j = 1, num_radial_elems
        call self.subgrid_mesh(j,i).write_out(800,self.subgrid_mesh(j,i).cell)
      end do
    end do
  end if

end subroutine setup_mesh

!---------------------------------------------------------------------------!
!< @details Interpolate the bed elevation and water surface elevation
!
!< @param[inout]  self
!---------------------------------------------------------------------------!
subroutine interp_depth_elev(self, cell, x_pos, y_pos, bed_elev, w_depth)

  Use Variables_Propwash

  Implicit None
  ! *** Dummy variables
  class(active_ship), intent(in)  :: self
  integer         , intent(in)    :: cell
  real(kind = rkd), intent(in)    :: x_pos
  real(kind = rkd), intent(in)    :: y_pos
  real(kind = rkd), intent(inout) :: bed_elev
  real(kind = rkd), intent(inout) :: w_depth

  ! *** Local variables
  integer :: L, LL
  integer :: i, j, k
  real(kind = rkd) :: bed_elev_1, weight, center_dist, zeta

  ! *** Initialize
  bed_elev_1 = 0.0
  weight     = 0.0
  zeta       = 0.0
  w_depth    = 0.0
  bed_elev   = 0.0
  do LL = 1, 9
    ! *** get the adjacent cell value
    !< @todo right now this assumes all propwash field liesa single cell based on the current ship position
    !! should modify to return the cell based on the x,y locations that are passed in
    L = adjacent_l(LL, cell)

    ! *** Only do if the cell is valid
    if( L >= 2 .AND. L <= LA )then
      ! *** INVERSE DISTANCE SQUARED interpolation, interpolate with center of surrounding cells
      center_dist = MAX( (x_pos - XCOR(L,5))**2 + (y_pos - YCOR(L,5))**2, 1D-8)

      bed_elev_1 = bed_elev_1 + BELV(L)/center_dist
      weight = weight + 1.0_rkd/center_dist

      ! *** for water surface interpolation
      zeta = zeta + (HP(L) + BELV(L))/center_dist
    end if
  end do
  if( weight > 0.0 )then
    bed_elev = bed_elev_1/weight
    zeta = zeta/weight

    ! *** Set the water surface elev
    w_depth = max(zeta - bed_elev_1/weight, hdry)
  endif

end subroutine interp_depth_elev

!---------------------------------------------------------------------------!
!< @details Interpolate the bed elevation and water surface elevation
!
!< @param[inout]  self
!---------------------------------------------------------------------------!
subroutine interp_z(self,cell, bed_elev, w_depth, z_loc)

  implicit none
  ! *** dummy variables
  class(active_ship), intent(in) :: self
  integer, intent(in) :: cell
  real(kind = rkd), intent(in) :: bed_elev
  real(kind = rkd), intent(in) :: w_depth
  real(kind = rkd), intent(inout) :: z_loc
  ! *** local variables
  integer :: k
  real(kind = rkd) :: z_sig

end subroutine interp_z

!---------------------------------------------------------------------------!
!< @details Interpolate bottom cell velocity
!
!< @param[inout]  self
!---------------------------------------------------------------------------!
subroutine interp_bot_vel(self, cell, x_pos, y_pos, u_vel, v_vel)

  implicit none
  ! *** dummy variables
  class(active_ship), intent(in) :: self
  integer, intent(in)     :: cell
  real(kind=rkd), intent(in) :: x_pos
  real(kind=rkd), intent(in) :: y_pos
  real(kind=rkd), intent(inout) :: u_vel
  real(kind=rkd), intent(inout) :: v_vel

  ! *** local variables
  integer :: ll, l_east, l_north, l
  integer :: k_bot, k_bot_cur
  real(kind = rkd) :: u_vel_b, v_vel_b
  real(kind = rkd) :: dist_center
  real(kind = rkd) :: xout, yout
  real(kind = rkd) :: vel_east, vel_north
  real(kind = rkd) :: vfactor, u_tmp, v_tmp
  real(kind = rkd) :: su_2, su_3, sv_2, sv_3

  VFACTOR = 0.5
  su_2 = 0.0
  su_3 = 0.0
  sv_2 = 0.0
  sv_3 = 0.0

  ! *** get bottom for current cell
  k_bot_cur = ksz(cell)

  ! *** go around 8 surrounding cells
  do ll = 1, 9

    l = adjacent_l(ll, cell)
    ! *** bottom layer
    k_bot = KSZ(l)
    ! *** make sure cell is valid
    if( l >= 2 .and. l <= LA )then

      ! *** get north and east cells
      l_east  = LEC(l)
      l_north = LNC(l)

      ! *** determine horizontal velocity at centroid bottom of the cell
      u_vel_b = 0.5*STCUV(l)*(RSSBCE(l)*U2(l_east, k_bot) + RSSBCW(l)*U2(l,k_bot) )
      v_vel_b = 0.5*STCUV(l)*(RSSBCN(l)*V2(l_north,k_bot) + RSSBCS(l)*V2(l,k_bot) )
      ! *** Determine distance from cell centroid
      dist_center = max( (x_pos - xcor(l,5))**2 + (y_pos - ycor(l,5))**2, 1d-8)

    else ! *** Edge cells, apply a zero face velocity

      u_tmp = slipfactor*U2(cell, k_bot_cur)
      v_tmp = slipfactor*V2(cell, k_bot_cur)

      IF( LL == 1 )THEN
        ! *** NORTHWEST
        XOUT = XCOR(cell,5) + ( DYP(cell)*CVE(cell) - DXP(cell)*CUE(cell) )
        YOUT = YCOR(cell,5) + ( DYP(cell)*CVN(cell) - DXP(cell)*CUN(cell) )
        dist_center = MAX((x_pos - XOUT)**2+(y_pos-YOUT)**2,1D-8)
      ELSEIF( LL == 2 )THEN
        ! *** NORTH FACE (MOVE WITH DRIFTER)
        XOUT = x_pos + 0.5*( DYP(cell)*CVE(cell) )
        YOUT = y_pos + 0.5*( DYP(cell)*CVN(cell) )
        dist_center = MAX( (x_pos-XOUT)**2 + (y_pos-YOUT)**2, 1D-8)
        v_tmp = -VFACTOR*MAX(V2(cell,k_bot_cur),0.0)
      ELSEIF( LL == 3 )THEN
        ! *** NORTHEAST
        XOUT = XCOR(cell,5) + ( DYP(cell)*CVE(cell) + DXP(cell)*CUE(cell) )
        YOUT = YCOR(cell,5) + ( DYP(cell)*CVN(cell) + DXP(cell)*CUN(cell) )
        dist_center = MAX((x_pos-XOUT)**2+(y_pos-YOUT)**2,1D-8)
      ELSEIF( LL == 4 )THEN
        ! *** WEST FACE (MOVE WITH DRIFTER)
        XOUT = x_pos - 0.5*( DXP(cell)*CUE(cell) )
        YOUT = y_pos - 0.5*( DXP(cell)*CUN(cell) )
        dist_center = MAX( (x_pos-XOUT)**2 + (y_pos-YOUT)**2, 1D-8)
        u_tmp = -VFACTOR*MIN(U2(LEC(cell),k_bot_cur),0.0)
      ELSEIF( LL == 6 )THEN
        ! *** EAST FACE (MOVE WITH DRIFTER)
        XOUT = x_pos + 0.5*( DXP(cell)*CUE(cell) )
        YOUT = y_pos + 0.5*( DXP(cell)*CUN(cell) )
        dist_center = MAX( (x_pos-XOUT)**2 + (y_pos-YOUT)**2, 1D-8)
        u_tmp = -VFACTOR*MAX(U2(cell,k_bot_cur),0.0)
      ELSEIF( LL == 7 )THEN
        ! *** SOUTHWEST
        XOUT = XCOR(cell,5) - ( DYP(cell)*CVE(cell) + DXP(cell)*CUE(cell) )
        YOUT = YCOR(cell,5) - ( DYP(cell)*CVN(cell) + DXP(cell)*CUN(cell) )
        dist_center = MAX((x_pos-XOUT)**2+(y_pos-YOUT)**2,1D-8)
      ELSEIF( LL == 8 )THEN
        ! *** SOUTH
        XOUT = x_pos - 0.5*( DYP(cell)*CVE(cell) )
        YOUT = y_pos - 0.5*( DYP(cell)*CVN(cell) )
        dist_center = MAX((x_pos-XOUT)**2+(y_pos-YOUT)**2,1D-8)
        v_tmp  = -VFACTOR*MIN(V2(LNC(cell),k_bot_cur),0.0)
      ELSEIF( LL == 9 )THEN
        ! *** SOUTHEAST
        XOUT = XCOR(cell,5) - ( DYP(cell)*CVE(cell) - DXP(cell)*CUE(cell) )
        YOUT = YCOR(cell,5) - ( DYP(cell)*CVN(cell) - DXP(cell)*CUN(cell) )
        dist_center = MAX((x_pos-XOUT)**2+(y_pos-YOUT)**2,1D-8)
      ENDIF

      l = cell
    end if

    ! *** Rotation?
    vel_east  = CUE(l)*u_tmp + CVE(l)*v_tmp
    vel_north = CUN(l)*u_tmp + CVN(l)*v_tmp

    ! *** Weighting
    su_2 = su_2 + vel_east / dist_center
    su_3 = su_3 + 1.0_rkd / dist_center

    sv_2 = sv_2 + vel_north / dist_center
    sv_3 = sv_3 + 1.0_rkd / dist_center

  end do

  ! *** Final interpolated values
  u_vel = su_2 / su_3
  v_vel = sv_2 / sv_3

end subroutine interp_bot_vel

!---------------------------------------------------------------------------!
!< @details Calculates propwash velocity for a given (ax,r) where ax is an
!           axial position and r is the radius away from that axial position
!< @author Zander Mausolff
!< @param[in]  ax - axial location
!< @param[in]  radius - radius away from the axial line
!< @param[out] velocity - returns the velocity determined from the emipirical
!! propwash velocity relationship
!---------------------------------------------------------------------------!
subroutine calc_velocity(self, ax, radius, bot_velocity, ieffluxonly)

  Use Mod_Ship

  implicit none

  ! *** dummy variables
  class(active_ship), intent(inout) :: self
  real(kind = rkd), intent(in)      :: ax               !< axial distance away from the starting point
  real(kind = rkd), intent(in)      :: radius           !< radius from the center point
  real(kind = rkd), intent(inout)   :: bot_velocity     !< radial velocity we want to calculate
  integer, intent(in)               :: ieffluxonly      !< Momentum only flag

  ! *** local
  integer          :: cell
  !< Ship position cell
  real(kind = rkd) :: efflux_vel, power
  real(kind = rkd) :: radius_vel_max                    !< location of the maximum radial velocity
  real(kind = rkd) :: radius_prop                       !< Radius of propeller
  real(kind = RKD) :: Aprime, Bprime, sigma             !< coeffs for vel_max_axial computation

  real(kind = RKD) :: vel_max_ax                        !< max axial velocity
  real(kind = RKD) :: EP, EPD, C1, rho

  real(kind = RKD) :: zone0, zone1, zone2               !< Propeller velocity regions

  ! *** Efflux zone - Maximum velocity at propeller
  if( self.rps >= 0.0 )then
    ! *** Based on angular speed
    !efflux_vel = 1.22 * self.rps**1.01 * self.ship.prop_diam**0.84 * self.ship.thrust_coeff**0.62   ! *** Hamill&Kee 2016
    efflux_vel  = 1.59 * self.rps * self.ship.prop_diam * sqrt(self.ship.thrust_coeff)               ! *** Jiang, etal 2019
    !efflux_vel = 1.33 * self.rps * self.ship.prop_diam * sqrt(self.ship.thrust_coeff)               ! *** Hamill 1987
  else
    ! *** Based on applied horsepower
    power = self.power / self.ship.num_props                                           ! *** Assume power is total ship power.  Calcs based on single prop
    if( self.ship.ducted == 0 )then
      ! *** Open propeller
      ! ***            HP              m/s   MPH/(m/s)
      EP = 23.57*power**0.974 - 2.3*(self.speed*2.237)**2 * sqrt(power)                ! *** Thrust Toutant Eq [lbs-force]  -  Ignores ambient current
      C1 = 0.71                                                                        ! *** Contraction coefficient from Maynord
    else
      ! *** Ducted propeller
      ! ***            HP              m/s   MPH/(m/s)
      EP = 31.82*power**0.974 - 5.4*(self.speed*2.237)**2 * sqrt(power)                ! *** Thrust [lbs-force]
      C1 = 1.0                                                                         ! *** Contraction coefficient from Maynord
    endif
    EP = max(EP, power)                                                                ! *** Minimum thrust.  Toutant breaks down for high speed, low power cases
    EP = 4.448*EP                                                                      ! *** Thrust [Newtons]

    ! *** Get water density
    if( bsc > 0.0 )then
      cell = self.pos.cell
      rho = RHOW(cell,KSZ(cell))
    else
      rho = 999.82                                                                     ! *** Water density @ 20 degC
    endif

    efflux_vel = 1.13/C1/self.ship.prop_diam * sqrt(EP/rho)
  endif
  self.efflux_vel = efflux_vel * efflux_mag_mult                                       ! *** efflux_mag_mult is a factor to account for subgrid turbulent losses
  if( ieffluxonly > 0 )then
    bot_velocity = efflux_vel
    return                                                                             ! *** Used if current iteration only needs momentum impacts
  endif

  vel_max_ax   = 0.0
  bot_velocity = 0.0

  ! *** cutoff distances for each flow regime
  zone0 = efflux_zone_mult * self.ship.prop_diam                ! *** Efflux zone
  zone1 = flow_est_zone1_mult * self.ship.prop_diam             ! *** Zone of flow establishment

  ! *** Velocities by Zone/Flow regime
  if( ax <= zone0 )then
    ! *** Efflux zone axial_mesh < 0.35Dp

    ! *** Max axial velocity
    vel_max_ax = efflux_vel

    ! *** radial velocity
    if( radius <= 0.75*self.ship.prop_diam )then
      bot_velocity = efflux_vel
    end if

  elseif( ax > zone0 .and. ax < zone1 )then
    ! *** Zone of flow estabishment  0.35Dp < axial distance < 3.25 Dp (Hamill & Kee 2016)
    vel_max_ax  = efflux_vel *(1.51 - 0.175 * (ax / self.ship.prop_diam) - 0.46 * self.ship.pitch_ratio)

    radius_prop = 0.5*self.ship.prop_diam
    radius_vel_max = 0.67 * (radius_prop - 0.5*self.ship.prop_hub_diam)                                    ! *** Radius of max vel at efflux plane (Rm0)

    ! *** determine radial velocity
    if( ax < radius_prop )then
      sigma = 0.5*radius_vel_max
    else
      sigma = 0.5*radius_vel_max + 0.075*(ax - radius_prop)
    end if
    bot_velocity = vel_max_ax * exp( -0.5*(radius - radius_vel_max)**2 / sigma**2 )

  else
    ! *** Zone of Established Flow
    Aprime = -11.4 * self.ship.thrust_coeff + 6.65 * self.ship.blade_area_ratio + 2.16 * self.ship.pitch_ratio;    ! eqn (2.32)
    if( Aprime > 1.0 )then
      Bprime = -(self.ship.thrust_coeff**(-0.216)) * self.ship.blade_area_ratio**1.024 / self.ship.pitch_ratio;    ! eqn (2.33)

      vel_max_ax = efflux_vel * Aprime * (ax / self.ship.prop_diam)**Bprime
    else
      ! *** If either BAR is not known or Aprime is too low, use Hashmi 1993
      vel_max_ax = efflux_vel * 0.638*exp(-0.097*(ax / self.ship.prop_diam))
    endif

    ! *** radial velocity calc
    bot_velocity = vel_max_ax * exp(-22.2*(radius/ax)**2)

  end if

  ! *** remove small velocities
  if( bot_velocity < epsilon )then
    bot_velocity = 0.0
  end if

end subroutine calc_velocity

!---------------------------------------------------------------------------!
!< @details Calculates the erosive flux from the propwash velocity field
!<          calculates single propeller velocities based on 2016 Hamill
!< @author Paul Craig, Luis Bastidas & Zander Mausolff
!< @param[inout] self
!---------------------------------------------------------------------------!
subroutine calc_erosive_flux(self, debug)

  implicit none

  ! *** Dummy variables
  class(active_ship), intent(inout) :: self
  logical, optional, intent(in)     :: debug

  ! *** Local variables
  integer :: i, j, k, m, q, isurf
  integer :: cell
  integer :: shear_counter
  real(kind = rkd) :: ang, ax, radius, bot_velocity, d1, d2, xp2, yp2, cos2, sin2, prop_off, vel2, velx, vely, xxx, yyy
  real(kind = rkd) :: bottom_roughness          !< Bottom roughness in meters
  real(kind = rkd) :: x_ax, y_ax, z_ax          !< (x,y,z) axial locations
  real(kind = rkd) :: x_bed, y_bed, z_bed       !< (x,y,z) locations at the bed
  real(kind = rkd) :: w_depth                   !< interpolated water depth
  real(kind = rkd) :: c_f                       !< shear friction factor
  real(kind = rkd) :: shear                     !< shear stress [dynes/cm^2]
  real(kind = rkd) :: shear_pa                  !< shear stress [Pa]
  real(kind = rkd) :: erosion                   !< erosion from propwash velocity field [g/cm^2]
  real(kind = rkd) :: prop_radius               !< Propeller radius [m]
  real(kind = rkd) :: elay(nscm)                !< erosion from propwash velocity field into suspension [g/cm^2]
  real(kind = rkd) :: ebld(nscm)                !< erosion from propwash velocity field into bedload [g/cm^2]
  real(kind = rkd) :: vb_max, vb_max1           !< Current iteration maximum bottom velocity

  character(31)    :: filename
  logical          :: bflag

  ! *** Ignore track points outside domain
  if( self.pos.cell < 2 ) return

  ! *** initialize
  elay = 0.0
  bot_velocity  = 0.0
  shear_counter = 0
  vb_max = 0.0
  vb_max1 = self.max_bot_vel

  ! ***
  prop_radius = self.prop_radius

  ! *** Get the Shear coefficient (Maynord 2000) based on depth at the propeller
  cell  = self.pos.cell
  x_bed = self.pos.x_pos
  y_bed = self.pos.y_pos
  call self.interp_depth_elev(cell, x_bed, y_bed, z_bed, w_depth)   ! *** Determines z_bed and water depth
  d2 = w_depth - self.draft
  d2 = max(d2,prop_radius)
  c_f = 0.01 * self.ship.prop_diam/d2                               ! *** Bottom friction coefficient - Use constant c_f (2021-06-28)

  ! *** Determine bottom shear stress at all intersection points
  !     of the bottom with propwash velocity field
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, m, isurf)   REDUCTION(MAX: vb_max)                                  &
  !$OMP                           PRIVATE(x_ax, y_ax, ax, prop_off, sin2, cos2, xp2, yp2, radius, d1, d2)           &
  !$OMP                           PRIVATE(vel2, erosion, shear, ang, velx, vely, ebld, xxx, yyy)                    &
  !$OMP                           FIRSTPRIVATE(elay, bot_velocity, shear_counter, vb_max1, prop_radius)             &
  !$OMP                           FIRSTPRIVATE(cell, x_bed, y_bed, z_bed, w_depth, c_f)                             &
  !$OMP                           SHARED(KSZ, RHOW, ISTRAN, LBED, LSEDZLJ, DTSEDJ, NSCM, PROP_ERO, PROP_BLD)        &
  !$OMP                           SHARED(self, num_axial_elems, num_radial_elems)
  do i = 1, num_axial_elems
    ! *** get axial location (x,y,z)
    x_ax = self.axial_mesh(i).x_pos         ! *** Model domain position X
    y_ax = self.axial_mesh(i).y_pos         ! *** Model domain position Y

    ! *** get axial distance away from ship prop
    ax = ((x_ax - self.pos.x_pos) **2 + (y_ax - self.pos.y_pos)**2)**(0.5)

    do j = 1, num_radial_elems
      if( self.subgrid_mesh(j,i).cell < 1 )then
        self.subgrid_mesh(j,i).z_pos = 0.0
        self.subgrid_mesh(j,i).var   = 0.0
        cycle
      endif

      ! *** get the cell index for the current (x,y)
      cell = self.subgrid_mesh(j,i).cell

      ! *** Determine the distance from the center axial line of the prop to the (x,y,z) that intersects the bed

      ! *** Get (x,y) of the horizontal plane
      x_bed = self.subgrid_mesh(j,i).x_pos
      y_bed = self.subgrid_mesh(j,i).y_pos

      ! *** interpolate the bed location based on the (x,y)
      call self.interp_depth_elev(cell, x_bed, y_bed, z_bed, w_depth)         ! *** Determines z_bed and water depth

      ! *** Vertical component of radius, i.e. depth at the mesh location
      d2 = w_depth - self.draft
      d2 = max(d2,prop_radius)

      ! *** Use superposition to calculate the total velocity for each point
      bot_velocity = 0.
      prop_off = -0.5*self.ship.dist_between_props*real(self.ship.num_props-1)
      do m=1,self.ship.num_props
        sin2 = prop_off*sin(self.heading)
        cos2 = prop_off*cos(self.heading)
        xp2 = x_ax + sin2
        yp2 = y_ax + cos2

        d1 = ( (x_bed - xp2)**2 + (y_bed - yp2)**2 )**0.5                     ! *** Horizontal component
        radius = sqrt( d1*d1 + d2*d2 )

        ! *** get the velocity at the sediment/water interface due to propwash
        call self.calc_velocity(ax, radius, vel2, 0)

        bot_velocity = bot_velocity + vel2

        prop_off = prop_off + self.ship.dist_between_props                    ! *** Offset for next propeller
      enddo

      self.subgrid_mesh(j,i).z_pos  = z_bed
      self.subgrid_mesh(j,i).var(1) = bot_velocity
      vb_max = max(vb_max, bot_velocity)

      ! *** Now bot_velocity is the magnitude of the velocity

      !< @todo should break up the velocity into x, y components and combine with ambient flow field
      !! Get ambient flow convert to units of [cm/s]

      !< @todo consider wave action

      erosion = 0.0

      ! *** Calculate SHEAR stress ONLY when there is a velocity value at the bed greater than cutoff
      if( bot_velocity > 1e-2 )then
        ! *** keep track of the number of times the shear is calculated
        shear_counter = shear_counter + 1

        ! *** Ambient velocities as of 2021-11-10 are ignored for analytical solution but used in for orginary sediment transport
        vel2 = bot_velocity                                                    ! *** Ignore ambient velocities

        ! *** Shear due to propwash velocities (Maynord 2000)
        !c_f = 0.01 * self.ship.prop_diam/d2                                  ! *** Bottom friction coefficient - Use constant c_f (2021-06-28)
        shear = 0.5 * Rhow(cell,ksz(cell)) * c_f * vel2**2                    ! *** Bottom shear in N/m**2        Using resulant magnitude vel2
        self.subgrid_mesh(j,i).var(2) = shear                                 ! *** Shear in Pascals
        !
        ! *** Determine erosion rate based on the shear stress
        if( ISTRAN(6)+ISTRAN(7) > 0 )then
          IF( LBED(cell) ) CYCLE                                              ! *** Bypass hard bottom cells

          ebld = 0.0
          elay = 0.0
          if( LSEDZLJ )then
            Call Calc_Prop_Erosion_SEDZLJ(cell, shear, elay, isurf)           ! *** Calculate erosion rates in g/cm^2

            elay(:) = elay(:)*self.subgrid_mesh(j,i).area                     ! *** Mass eroded for mesh cell per class (g)
            self.subgrid_mesh(j,i).ero(:) = elay(:)                           ! *** Total mass eroded for mesh cell per class (g)
            !ebld(:) = ebld(:)*self.subgrid_mesh(j,i).area                    ! *** Mass eroded into bedload layer for mesh cell per class (g)   ToDo
            !self.subgrid_mesh(j,i).ero(:) = (elay(:) + ebld(:))              ! *** Total mass eroded for mesh cell per class (g)                ToDo
          else
            shear = shear*0.001                                               ! *** Density normalized shear in m2/s2
            Call Calc_Prop_Erosion_Original(cell, shear, elay, ebld)          ! *** Calculate erosion rates
                                                                              ! *** elay in g/m^2 and ebld in g/m/s

            elay(:) = elay(:)*self.subgrid_mesh(j,i).area                     ! *** Mass eroded for mesh cell per class (g)
            
            ! *** Convert bedload to get total mass eroded for mesh cell per class (g)
            self.subgrid_mesh(j,i).ero(:) = (elay(:) + ebld(:)*self.subgrid_mesh(j,i).width*DTSEDJ)
          endif
        else
          self.subgrid_mesh(j,i).ero(:) = 0.0
          self.subgrid_mesh(j,i).var(3) = 0.0
        endif

      else
        self.subgrid_mesh(j,i).ero(:) = 0.0
        self.subgrid_mesh(j,i).var(2) = 0.0
        self.subgrid_mesh(j,i).var(3) = 0.0
      end if

    end do
  end do
  !$OMP END PARALLEL DO

  ! *** Accumulate the global propwash erosion (moved outside OMP due to "racing")
  do i = 1, num_axial_elems
    do j = 1, num_radial_elems
      if( self.subgrid_mesh(j,i).cell < 1 )then
        cycle
      endif

      cell = self.subgrid_mesh(j,i).cell          ! *** Gather erosion and put it into the cell the propwash component exists
      prop_ero(cell,1:nscm) = prop_ero(cell,1:nscm) + self.subgrid_mesh(j,i).ero(:)               ! *** (g)
      !prop_bld(cell,1:nscm) = prop_bld(cell,1:nscm) + ebld(1:nscm)/DBLE(self.mesh_count(cell))   ! *** Average unit discharge (g/m/s)
    enddo
  enddo
  
  self.max_bot_vel = 0.25 * vb_max          ! *** Fraction of maximum bottom velocity to be used to prevent double counting if ISPROPWASH = 2

  if( self.freq_out > 0. )then
    if( timeday >= self.timesnap )then
      if( self.power /= 0.0 .or. self.nsnap == 0 )then
        ! *** Output subgrid mesh data to binary format
        write(filename, '(A8,I9.9,A4)')    'pw_mesh_', self.mmsi, '.out'
        inquire(file = OUTDIR//filename, EXIST = bflag)
        if( .not. bflag )then
          Call create_linkage(self)
        endif
        close(801)
        open(unit = 801,file = OUTDIR//filename,status='OLD',position='APPEND',form='BINARY',shared)

        write(801) timeday
        write(801) real(self.pos.x_pos - center_x, 4)
        write(801) real(self.pos.y_pos - center_y, 4)
        write(801) real(self.stern.x_pos - center_x, 4)
        write(801) real(self.stern.y_pos - center_y, 4)
        write(801) real(self.antenna.x_pos - center_x, 4)
        write(801) real(self.antenna.y_pos - center_y,4)
        write(801) real(self.heading,4)

        write(801) int(self.pos.cell,4)

        write(801) ((real(self.subgrid_mesh(j,i).x_pos - center_x, 4), j = 1, num_radial_elems), i = 1, num_axial_elems)
        write(801) ((real(self.subgrid_mesh(j,i).y_pos - center_y, 4), j = 1, num_radial_elems), i = 1, num_axial_elems)
        write(801) ((real(self.subgrid_mesh(j,i).z_pos,4), j = 1, num_radial_elems), i = 1, num_axial_elems)

        do i = 1, num_axial_elems
          do j = 1, num_radial_elems
            if( self.subgrid_mesh(j,i).cell < 1 ) cycle
            erosion = sum(self.subgrid_mesh(j,i).ero(:))                                 ! *** Total mass eroded for subgrid area (g)
            self.subgrid_mesh(j,i).var(3) = erosion/DTSEDJ/self.subgrid_mesh(j,i).area   ! *** Mass erosion rate:  LSEDZLJ: g/cm^2/s,  Orig: g/m^2/s
            if( .not. LSEDZLJ )then
              self.subgrid_mesh(j,i).var(3) = self.subgrid_mesh(j,i).var(3)/10000.       ! *** Convert g/m^2/s to g/cm^2/s
            endif
          enddo
        enddo

        do k = 1, 3
          write(801) ((real(self.subgrid_mesh(j,i).var(k),4), j = 1, num_radial_elems), i = 1, num_axial_elems)
        enddo
        write(801) ((int(self.subgrid_mesh(j,i).cell,4), j = 1, num_radial_elems), i = 1, num_axial_elems)

        close(801)
        self.nsnap = self.nsnap + 1
        iwrite_pwmesh = 1                                                                ! *** Global flag to indicate a pw_mesh file has been updated
      endif

      ! *** Check on snapshot frequency option
      if( self.timesnap == -99999. )then
        self.timesnap = int(timeday*10000.)*0.0001
      endif
      self.timesnap = self.timesnap + self.freq_out/1440.

      ! *** Address periods of power = 0 within an active track
      if( self.timesnap < timeday )then
        self.timesnap = int(timeday*10000.)*0.0001 + self.freq_out/1440.
      endif
    endif
  end if

end subroutine calc_erosive_flux

!---------------------------------------------------------------------------!
!< @details writes out the (x,y,z coordinates of the vel field
!< @param[in] self
!---------------------------------------------------------------------------!
subroutine write_xyz(self, unit_num)

  implicit none
  ! *** Dummy variables
  class(active_ship), intent(inout) :: self
  integer, intent(in) :: unit_num
  ! *** Local variables
  integer :: i,j

  do i = 1, num_axial_elems
    write(unit_num,'(a)') ' '
    !call self.axial_mesh(i).write_out(unit_num)
    do j = 1, num_radial_elems
      !call self.radial_mesh(j,i).write_out(unit_num)

    end do
  end do

end subroutine write_xyz

!---------------------------------------------------------------------------!
!< @details deallocates ship arrays
!< @param[in] self
!---------------------------------------------------------------------------!
subroutine delete_active_ship(self)

  implicit none

  class(active_ship), intent(inout) :: self

  deallocate(self.tracks)

end subroutine delete_active_ship

subroutine create_linkage(self)

  implicit none

  class(active_ship), intent(inout) :: self
  character(32) :: filename
  integer :: verNum,headerSize,blockSize,varCount

  verNum = 10300
  varCount = 3           ! < bot_vel, shear, erosion
  headerSize = 17*4        ! (4*4 + 4 + 32 + 4*4)
  blockSize = 4*(10+(4+varCount)*num_axial_elems*num_radial_elems)

  write(filename, '(A8,I9.9,A4)')    'pw_mesh_',self.mmsi, '.out'
  open(unit = 801,file = OUTDIR//filename,status='UNKNOWN')
  close(801,status='DELETE')
  open(unit = 801,file = OUTDIR//filename,status='UNKNOWN',access='SEQUENTIAL',form='BINARY')
  filename = self.name
  write(801) int(verNum,4),int(headerSize,4),int(blockSize,4),int(varCount,4)     ! 4*4
  write(801) int(self.mmsi,4),filename                                            ! 4 + 32
  write(801) int(num_axial_elems, 4),int(num_radial_elems, 4)                     ! 2*4
  write(801) int(center_x, 4),int(center_y, 4)                                    ! 2*4
  close(801)

end subroutine create_linkage

End module Mod_Active_Ship
