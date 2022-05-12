! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Performs MPI Send/Receive with automatic handling of precision
! @author Paul M. Craig
! @date 2021-06-21
  
Module Mod_DSI_SendRecv

  USE GLOBAL
  USE MPI
  Use Variables_MPI
  
  IMPLICIT NONE
  
  SAVE

  Public :: DSI_SEND

  !***Set up the geneic interface for commuicating ghost cells
  !***This will select the the proper array size depending on what is passed in
  Interface DSI_SEND

  Module Procedure DSI_SEND_Integer,   &
                   DSI_SEND_Real4,     &
                   DSI_SEND_Real8
  End interface

  Public :: DSI_RECV

  Interface DSI_RECV

  Module Procedure DSI_RECV_Integer,   &
                   DSI_RECV_Real4,     &
                   DSI_RECV_Real8
  End interface

  ! *** Technical Note:  The STATUS variable returned by MPI_SEND and MPI_RECV is the number of BYTES sent or received.
  
  ! *** Next are each of the subroutines for the procedure
  Contains
 
SUBROUTINE DSI_SEND_Integer(ValIn, iLen, iDir)

  Integer, Intent(in) :: ValIn(*)
  Integer, Intent(in) :: iLen, iDir
  INTEGER             :: IERR, status_message(MPI_Status_Size)

  IF( num_Processors > 1 )THEN
    IERR = 0
    CALL MPI_SEND(ValIn, iLen, mpi_integer, iDir, process_id, comm_2d, IERR)
  ENDIF
  
END SUBROUTINE DSI_SEND_Integer

SUBROUTINE DSI_SEND_Real4(ValIn, iLen, iDir)

  Real(4), Intent(in) :: ValIn(*)
  Integer, Intent(in) :: iLen, iDir
  INTEGER             :: IERR

  IF( num_Processors > 1 )THEN
    IERR = 0
    CALL MPI_SEND(ValIn, iLen, mpi_real4, iDir, process_id, comm_2d, IERR)
  ENDIF
  
END SUBROUTINE DSI_SEND_Real4

SUBROUTINE DSI_SEND_Real8(ValIn, iLen, iDir)

  Real(8), Intent(in) :: ValIn(*)
  Integer, Intent(in) :: iLen, iDir
  INTEGER             :: IERR

  IF( num_Processors > 1 )THEN
    IERR = 0
    CALL MPI_SEND(ValIn, iLen, mpi_real8, iDir, process_id, comm_2d, IERR)
  ENDIF
  
END SUBROUTINE DSI_SEND_Real8

! ******************************************************************************************************************************
! @details Performs MPI Send with automatic handling of precision

SUBROUTINE DSI_RECV_Integer(Valout, iLen, iDir)

  Integer, Intent(out) :: Valout(*)
  Integer, Intent(in)  :: iLen, iDir
  INTEGER             :: IERR, status_message(MPI_Status_Size)

  IF( num_Processors > 1 )THEN
    IERR = 0
    CALL MPI_RECV(Valout, iLen, mpi_integer, iDir, iDir, comm_2d, status_message, IERR)
  ENDIF
  
END SUBROUTINE DSI_RECV_Integer

SUBROUTINE DSI_RECV_Real4(Valout, iLen, iDir)

  Real(4), Intent(out) :: Valout(*)
  Integer, Intent(in)  :: iLen, iDir
  INTEGER             :: IERR, status_message(MPI_Status_Size)

  IF( num_Processors > 1 )THEN
    IERR = 0
    CALL MPI_RECV(Valout, iLen, mpi_real4, iDir, iDir, comm_2d, status_message, IERR)
  ENDIF
  
END SUBROUTINE DSI_RECV_Real4

SUBROUTINE DSI_RECV_Real8(Valout, iLen, iDir)

  Real(8), Intent(out) :: Valout(*)
  Integer, Intent(in) :: iLen, iDir
  INTEGER             :: IERR, status_message(MPI_Status_Size)

  IF( num_Processors > 1 )THEN
    IERR = 0
    CALL MPI_RECV(Valout, iLen, mpi_real8, iDir, iDir, comm_2d, status_message, IERR)
  ENDIF
  
END SUBROUTINE DSI_RECV_Real8

END MODULE Mod_DSI_SendRecv
  