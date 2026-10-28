! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Performs MPI Send/Receive with automatic handling of precision
! @author Paul M. Craig
! @date 2021-06-21
  
Module Mod_DSI_SendRecv

  use GLOBAL
  use MPI
  use Variables_MPI
  
  implicit none
  
  save

  Public :: DSI_SEND

  !***Set up the geneic interface for commuicating ghost cells
  !***This will select the the proper array size depending on what is passed in
  interface DSI_SEND

  Module Procedure DSI_SEND_Integer,   &
                   DSI_SEND_Real4,     &
                   DSI_SEND_Real8
  End interface

  Public :: DSI_RECV

  interface DSI_RECV

  Module Procedure DSI_RECV_Integer,   &
                   DSI_RECV_Real4,     &
                   DSI_RECV_Real8
  End interface

  ! *** Technical Note:  The STATUS variable returned by MPI_SEND and MPI_RECV is the number of BYTES sent or received.
  
  ! *** Next are each of the subroutines for the procedure
  Contains
 
SUBROUTINE DSI_SEND_Integer(ValIn, iLen, iDir)

  integer, Intent(in) :: ValIn(*)
  integer, Intent(in) :: iLen, iDir
  integer             :: IERR, status_message(MPI_Status_Size)

  if( num_Processors > 1 )then
    IERR = 0
    call MPI_SEND(ValIn, iLen, mpi_integer, iDir, process_id, DSIcomm, IERR)
  endif
  
END SUBROUTINE DSI_SEND_Integer

SUBROUTINE DSI_SEND_Real4(ValIn, iLen, iDir)

  real(4), Intent(in) :: ValIn(*)
  integer, Intent(in) :: iLen, iDir
  integer             :: IERR

  if( num_Processors > 1 )then
    IERR = 0
    call MPI_SEND(ValIn, iLen, mpi_real4, iDir, process_id, DSIcomm, IERR)
  endif
  
END SUBROUTINE DSI_SEND_Real4

SUBROUTINE DSI_SEND_Real8(ValIn, iLen, iDir)

  real(8), Intent(in) :: ValIn(*)
  integer, Intent(in) :: iLen, iDir
  integer             :: IERR

  if( num_Processors > 1 )then
    IERR = 0
    call MPI_SEND(ValIn, iLen, mpi_real8, iDir, process_id, DSIcomm, IERR)
  endif
  
END SUBROUTINE DSI_SEND_Real8

! ******************************************************************************************************************************
! @details Performs MPI Send with automatic handling of precision

SUBROUTINE DSI_RECV_Integer(Valout, iLen, iDir)

  integer, Intent(out) :: Valout(*)
  integer, Intent(in)  :: iLen, iDir
  integer             :: IERR, status_message(MPI_Status_Size)

  if( num_Processors > 1 )then
    IERR = 0
    call MPI_RECV(Valout, iLen, mpi_integer, iDir, iDir, DSIcomm, status_message, IERR)
  endif
  
END SUBROUTINE DSI_RECV_Integer

SUBROUTINE DSI_RECV_Real4(Valout, iLen, iDir)

  real(4), Intent(out) :: Valout(*)
  integer, Intent(in)  :: iLen, iDir
  integer             :: IERR, status_message(MPI_Status_Size)

  if( num_Processors > 1 )then
    IERR = 0
    call MPI_RECV(Valout, iLen, mpi_real4, iDir, iDir, DSIcomm, status_message, IERR)
  endif
  
END SUBROUTINE DSI_RECV_Real4

SUBROUTINE DSI_RECV_Real8(Valout, iLen, iDir)

  real(8), Intent(out) :: Valout(*)
  integer, Intent(in) :: iLen, iDir
  integer             :: IERR, status_message(MPI_Status_Size)

  if( num_Processors > 1 )then
    IERR = 0
    call MPI_RECV(Valout, iLen, mpi_real8, iDir, iDir, DSIcomm, status_message, IERR)
  endif
  
END SUBROUTINE DSI_RECV_Real8

END MODULE Mod_DSI_SendRecv
  