! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  Module Mod_MPI_Helper_Functions

  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_Propwash
  use GLOBAL, only : KC, ISDRY, HDRY, HP
  use XYIJCONV

  implicit none

  private
  public :: get_nearest_cell, get_subdomain_id, r8_atan, Get_Cell


  !interface get_nearest_cell
  !    !procedure : get_nearest_cell_init
  !    procedure : get_nearest
  !end interface get_nearest_cell

  contains

  !---------------------------------------------------------------------------!
  !< @details determines the nearest cell (local subdomain index) to the current x,y position
  !< @author Zander Mausolff
  !< @param[in]  x_position
  !< @param[in]  y_position
  !< @param[out] cell_index
  !---------------------------------------------------------------------------!
  function get_nearest_cell(x_position, y_position) result(cell_index)
    ! *** Dummy variables
    real(kind = RKD), intent(in) :: x_position
    real(kind = RKD), intent(in) :: y_position

    ! *** Local variables
    real(kind = RKD), Allocatable, Dimension(:) :: Distance_to_Cells
    integer :: Min_L_loc_index(1)
    integer :: cell_index

    allocate(Distance_to_Cells(LA_Global))
    Distance_to_Cells = 0.0

    !   can constrain the movemet of a ship to only going to the adjacet cells
    ! *** Linear search for the nearest cell
    Distance_to_Cells(2:LA_Global) = SQRT((x_position - XCOR_Global(2:LA_Global,5))**2 + &
      (y_position - YCOR_Global(2:LA_Global,5))**2)
    ! *** Determine the minimum cell location based on the index of the distance to each cell
    Min_L_loc_index = minloc(Distance_to_Cells(2:LA_Global)) + 1

    ! *** Get the local cell index
    cell_index = Map2Local(Min_L_loc_index(1)).LL

    ! *** Check if there is water in this cell
    if( ISDRY > 0 .and. HP(cell_index) < HDRY  )then
      call STOPP("***WARNING*** the ship has entered a cell that has gone dry")
    endif

    deallocate(Distance_to_Cells)

  end function

  !-----------------------------------------------!
  !< @details Returns the subdomain process id for the current cell, excluding ghost cells
  !< @author Zander Mausolff
  !< @param[in]  cell_index
  !< @param[out] subdomain_id
  !---------------------------------------------------------------------------!
  pure function get_subdomain_id(cell_index) result (subdomain_id)
    ! *** Dummy variables
    integer, intent(in) :: cell_index

    ! *** Local variables
    integer :: subdomain_id

    subdomain_id = sorted_loc_to_glob(cell_index).process

  end function get_subdomain_id


  function r8_atan ( y, x )

    !*****************************************************************************80
    !
    !! R8_ATAN computes the inverse tangent of the ratio Y / X.
    !
    !  Discussion:
    !
    !    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
    !    the built in functions ATAN and ATAN2 already do.
    !
    !    However:
    !
    !    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
    !      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
    !      and [-PI,+PI] respectively;
    !
    !    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
    !     function by contrast always returns an angle in the first or fourth
    !     quadrants.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) Y, X, two quantities which represent the
    !    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
    !
    !    Output, real ( kind = 8 ) R8_ATAN, an angle between 0 and 2 * PI, whose
    !    tangent is (Y/X), and which lies in the appropriate quadrant so that
    !    the signs of its cosine and sine match those of X and Y.
    !
    implicit none

    real ( kind = 8 ) abs_x
    real ( kind = 8 ) abs_y
    real ( kind = 8 ) r8_atan
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) theta_0
    real ( kind = 8 ) value
    real ( kind = 8 ) x
    real ( kind = 8 ) y
    !
    !  Special cases:
    !
    if( x == 0.0D+00  )then

      if( 0.0D+00 < y  )then
        value = r8_pi / 2.0D+00
      elseif(  y < 0.0D+00  )then
        value = 3.0D+00 * r8_pi / 2.0D+00
      elseif(  y == 0.0D+00  )then
        value = 0.0D+00
      endif

    elseif(  y == 0.0D+00  )then

      if( 0.0D+00 < x  )then
        value = 0.0D+00
      elseif(  x < 0.0D+00  )then
        value = r8_pi
      endif
      !
      !  We assume that ATAN2 is correct when both arguments are positive.
      !
    else

      abs_y = abs ( y )
      abs_x = abs ( x )

      theta_0 = atan2 ( abs_y, abs_x )

      if( 0.0D+00 < x .and. 0.0D+00 < y  )then
        value = theta_0
      elseif(  x < 0.0D+00 .and. 0.0D+00 < y  )then
        value = r8_pi - theta_0
      elseif(  x < 0.0D+00 .and. y < 0.0D+00  )then
        value = r8_pi + theta_0
      elseif(  0.0D+00 < x .and. y < 0.0D+00  )then
        value = 2.0D+00 * r8_pi - theta_0
      endif

    endif

    r8_atan = value

    return
  end function r8_atan

  !---------------------------------------------------------------------------!
  !< @details Determines the cell (local subdomain index) to the current x,y position
  !<          Otherwize, returns 0 (not in domain)
  !< @author Paul M. Craig
  !< @param[in]  x_position, y_position in meters
  !< @param[in]  L                      Cell index
  !< @param[out] cell_index             Cell Index > 1 if valid
  !---------------------------------------------------------------------------!
  function get_cell(L, x_position, y_position) result(cell_index)
    ! *** Dummy variables
    integer,          intent(in) :: L
    real(kind = RKD), intent(in) :: x_position
    real(kind = RKD), intent(in) :: y_position
    integer                      :: cell_index, i, L1
    
    ! *** Test current cell
    if( L > 1 )then
      if( insidecell(L, x_position, y_position) )then
        cell_index = L
        return
      endif
      
      ! *** Search adjacent cells
      ! *** Order of cells
      ! ***   1  2  3
      ! ***   4  5  6
      ! ***   7  8  9
      do i = 1,9
        L1 = adjacent_l(i,L)
        if( L1 > 0 )then
          if( insidecell(L1, x_position, y_position) )then
            cell_index = L1
            return
          endif
        endif
      enddo
    
    endif
    
    ! *** Point not nearby, widen search
    cell_index = 0
    
    ! *** Still have not assigned L.  Search entire domain
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(L1) REDUCTION(MAX:cell_index)
    do L1 = 2,LA
      if( cell_index > 0 ) cycle       ! *** Once cell_index is defined fast forward as quickly as possible for all domains
      if( insidecell(L1, x_position, y_position) )then
        cell_index = L1
      endif
    enddo
    !$OMP END PARALLEL DO
    
    return

  end function get_cell

  End module Mod_MPI_Helper_Functions
