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
! @details Initializes MPI environment, creates communicator shared by all processes
! @author  Zander Mausolff
! @date  5/22/2019
!---------------------------------------------------------------------------!
Subroutine Initialize_MPI

  use GLOBAL
  use MPI
  use Variables_MPI

  implicit none

  !***Dummy variables

  !***Local variables
  integer :: ierr     ! Error value needed
  integer :: required_thread
  integer :: provided_thread_sup
  !---------------------------------------------------------------------------!

  !***Required first call to start the environment
  !Call MPI_INIT(ierr)

  ! @todo look back at this when we need multithreading + MPI
  !***Ensure that only one thread will make MPI calls
  required_thread = MPI_Thread_Funneled
  call MPI_Init_thread( required_thread, provided_thread_sup, ierr)

  !***Determine the rank of the calling process in the communicator
  call MPI_COMM_RANK(MPI_Comm_World, process_id, ierr)

  !***Get total number of processors available
  call MPI_COMM_SIZE(MPI_Comm_World, num_Processors, ierr)

  !***Check thread level support
  !    if( provided_thread_sup .lt. required_thread )then
  !
  !         ! Insufficient support, degrade to 1 thread and warn the user
  !
  !         if( process_id .eq. 0 )then
  !            write(*,*) "Warning:  This MPI implementation provides ",   &
  !                        &         "insufficient threading support."
  !         endif
  !
  !         call omp_set_num_threads(1)
  !
  !    endif

  !  *** Define MPI real

!***Required first call to start the environment
    !Call MPI_INIT(ierr)
    
!***Determine the rank of the calling process in the communicator
    call MPI_COMM_RANK(MPI_Comm_World, process_id, ierr)

  !***Setup File for writing info out from each process
  call Setup_MPI_Debug_File

  call WriteBreak(mpi_log_unit)
  write (mpi_log_unit, '(a,I3)' ) 'Hello from process: ', process_id
  write (mpi_log_unit, '(a,I3)' ) 'The master process ID is:   ', master_id
  write (mpi_log_unit, '(a,I3)' ) 'Total # of processes: ', num_Processors
  call WriteBreak(mpi_log_unit)

  DoublePrecision = KIND(DXU) == 8
  
End Subroutine Initialize_MPI

