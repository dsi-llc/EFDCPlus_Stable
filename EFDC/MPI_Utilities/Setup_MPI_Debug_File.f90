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
  ! @details Sets up file for each process that is writes out various values
  !! for debugging purposes
  ! @date 8/16/2019
  ! @author Zander Mausolff
  !---------------------------------------------------------------------------!

  Subroutine Setup_MPI_Debug_File

  Use GLOBAL
  Use MPI
  Use Variables_MPI
  Use Variables_MPI_Write_Out
  USE IFPORT

  Implicit none

  ! *** Local
  integer :: ierr, RES
  character(24) :: mpi_filename

  ! *** Get unique unit number for writting out details of each MPI process
  mpi_log_unit = mpi_log_unit + process_id
  mpi_comm_unit = mpi_comm_unit + process_id
  mpi_mapping_unit = mpi_mapping_unit + process_id
  mpi_err_unit = mpi_err_unit + process_id

#ifdef _WIN
  OUTDIR='#output\'
#else
  OUTDIR='#output/'
#endif
  If( process_id == master_id )THEN
    RES = MAKEDIRQQ('#output')
    IF( .NOT. RES )THEN
      RES = GETLASTERRORQQ( )
      IF( RES == ERR$NOENT )THEN
        CALL STOPP('THE PATH FOR $OUTPUT IS NOT FOUND!')
      ELSEIF (RES == ERR$ACCES )THEN
        CALL STOPP('CANNOT CREATE THE FOLDER: PERMISSION DENIED!')
      ENDIF
    ENDIF
  endif
  Call MPI_barrier(MPI_Comm_World, ierr)

  ! *** Create output file for each processor to record inputted parameters, etc.
  write(mpi_filename, '(A13,I3.3,A4)')    'log_mpi_proc_',process_id, '.log'
  Open(unit = mpi_log_unit, status='replace', file = OUTDIR//mpi_filename)
  mpi_log_file = mpi_filename
  write(mpi_log_unit, '(33("*"),/,"*** EFDC+ Version: ",A10," ***",/,33("*"))') EFDC_VER
  
  ! *** Write the file out
  IF( MPI_DEBUG_FLAG == .TRUE. )THEN
    write(mpi_filename, '(A13,I3.3,A4)')  'comm_mpi_proc_',process_id, '.log'
    Open(unit = mpi_comm_unit, status='replace', file = OUTDIR//mpi_filename)
    write(mpi_comm_unit, '(33("*"),/,"*** EFDC+ Version: ",A10," ***",/,33("*"))') EFDC_VER
  endif

  ! *** Write the file out
  write(mpi_filename, '(A13,I3.3,A4)')    'map_mpi_proc_',process_id, '.log'
  Open(unit = mpi_mapping_unit, status='replace', file = OUTDIR//mpi_filename)
  write(mpi_mapping_unit, '(33("*"),/,"*** EFDC+ Version: ",A10," ***",/,33("*"))') EFDC_VER

  End subroutine Setup_MPI_Debug_File

  