! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  !> @details Figures out the starting and ending valid L index ignoring the
  !> ghost cells
  !
  !---------------------------------------------------------------------------!
  !> @author Zander Mausolff
  !> @date 9/15/2019
  !---------------------------------------------------------------------------!

  Subroutine Get_LA_Local_No_Ghost()

  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Mapping
  use MPI

  implicit none

  ! *** Local variables
  integer :: I, J, L
  integer :: ierr

  L = 0
  do J = 3,JC-2
    do I = 3,IC-2
      if( LIJ(I,J) > 0 )then
        L = L + 1
      endif
    enddo
  enddo

  ! *** Get the number of active cells excluding anything in the ghost regions
  LA_Local_no_ghost = L

  End Subroutine Get_LA_Local_No_Ghost
