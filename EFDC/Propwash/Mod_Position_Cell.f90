! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Mod_Position_Cell

  use GLOBAL, only : RKD, RK4
  use Mod_MPI_Helper_Functions
  use Mod_Position

  implicit none

  private
  public :: position_cell

  ! *** Extend the position type, to get (x,y,t)
  type, extends(position) :: position_cell

    Integer(kind = RK4) :: cell          = 0     !< L cell index the ship is in
    Integer(kind = RK4) :: sub_domain_id = 0     !< sub-domain the ship currently belongs in
    real(kind = rkd)    :: freq_out      = 0.0   !< If > 0, overwrites the ship defined output freqency
    real(kind = rkd)    :: speed         = 0.0   !< [meters/second]
    real(kind = rkd)    :: course        = 0.0   !< Compass direction the ship is moving toward
    real(kind = rkd)    :: heading       = 0.0   !< Compass direction the ship is pointed
    real(kind = rkd)    :: draft         = 0.0   !< Draft below water (used with prop_offset to set shaft depth) [meters]
    real(kind = rkd)    :: rps           = 0.0   !< Applied revolutions per second
    real(kind = rkd)    :: power         = 0.0   !< Applied power:  0.0 to  1.0:  Fraction of max_power
                                                 !<                 <0. to -1.0:  Fraction of max_rps
                                                 !<                 -999 Compute from ship parameters (delme TBD)
  contains

  procedure, pass(self) :: write_out_cell

  end type position_cell

  ! *** define interface for object constructor
  interface position_cell
    module procedure :: position_cell_constructor
  end interface position_cell

  contains
  !---------------------------------------------------------------------------!
  !< @details Cell position custom constructor
  !< @param[in] x_position
  !< @param[in] y_position
  !< @param[in] curr_time
  !< @param[out] self(position_cell) - implied return of type bound function
  !---------------------------------------------------------------------------!
  type(position_cell) function position_cell_constructor(x_position, y_position, curr_time, test) result(self)
  
    ! *** Dummy variables
    real(kind = RKD), intent(in) :: x_position
    real(kind = RKD), intent(in) :: y_position
    real(kind = RKD), intent(in) :: curr_time
    logical, optional, intent(in) :: test !> bypasses calls to cell location functions that require an actual EFDC+ geometry to be instantiated
    ! *** Local variables
    integer :: cell

    self.x_pos = x_position
    self.y_pos = y_position
    self.time  = curr_time

    ! *** Set cell values to zero for unit testing
    if(present(test) )then
      self.cell = 0
      self.sub_domain_id = 0
    else

      self.cell  = get_nearest_cell(x_position, y_position)

      ! *** Check that the cell returns is valid
      if(self.cell < 2  )then
        call STOPP('***Warning*** Cell index is invalid for this ships (x,y) location')
      else
        ! *** Deterine which process id the cell belongs in
        self.sub_domain_id = get_subdomain_id(self.cell)
      endif
    endif

  end function position_cell_constructor

  ! *** write out the position info
  subroutine write_out_cell(self, unit_num)

    implicit none

    class(position_cell), intent(in) :: self
    integer, optional, intent(in) :: unit_num

    if(present(unit_num) )then

      write(unit_num, '(a,f15.5)') 'x position:   ', self.x_pos
      write(unit_num, '(a,f15.5)') 'y position:   ', self.y_pos
      write(unit_num, '(a,f15.5)') 'time:         ', self.time
      write(unit_num, '(a,i4)')    'Subdomain ID: ', self.sub_domain_id
      write(unit_num, '(a,i8)')    'Cell index:   ', self.cell
    else
      ! *** write to some scratch unit
      write(884, '(a,f15.5)') 'x position:   ', self.x_pos
      write(884, '(a,f15.5)') 'y position:   ', self.y_pos
      write(884, '(a,f15.5)') 'time:         ', self.time
      write(884, '(a,i4)')    'Subdomain ID: ', self.sub_domain_id
      write(884, '(a,i8)')    'Cell index:   ', self.cell
    endif

  end subroutine write_out_cell

End module Mod_Position_Cell