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
  !> @details Figures out the starting and ending valid L index ignoring the
  !> ghost cells
  !
  !---------------------------------------------------------------------------!
  !> @author Zander Mausolff
  !> @date 9/15/2019
  !---------------------------------------------------------------------------!

#ifdef _MPI
  Subroutine Get_LA_Local_No_Ghost()

  Use GLOBAL
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use MPI

  Implicit None

  ! *** Local variables
  Integer :: I, J, L
  Integer :: ierr

  L = 0
  DO J = 3,JC-2
    DO I = 3,IC-2
      IF( LIJ(I,J) > 0 )THEN
        L = L + 1
      ENDIF
    END DO
  END DO

  ! *** Get the number of active cells excluding anything in the ghost regions
  LA_Local_no_ghost = L

  End Subroutine Get_LA_Local_No_Ghost
#endif