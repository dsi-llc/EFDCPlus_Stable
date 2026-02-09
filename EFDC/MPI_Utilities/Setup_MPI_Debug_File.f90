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
  ! @details Sets up file for each process that is writes out various values
  !! for debugging purposes
  ! @date 8/16/2019
  ! @author Zander Mausolff
  !---------------------------------------------------------------------------!

  Subroutine Setup_MPI_Debug_File

  use GLOBAL
  use MPI
  use Variables_MPI
  use Variables_MPI_Write_Out
#ifndef GNU  
  USE IFPORT
#endif

  implicit none

  ! *** Local
  integer :: ierr, RES
  character(24) :: mpi_filename
  character*200 :: STR

  ! *** Get unique unit number for writting out details of each MPI process
  mpi_efdc_out_unit = mpi_efdc_out_unit + process_id
  mpi_log_unit = mpi_log_unit + process_id
  mpi_comm_unit = mpi_comm_unit + process_id
  mpi_mapping_unit = mpi_mapping_unit + process_id
  mpi_qdwaste_unit = mpi_qdwaste_unit + process_id
  mpi_error_unit = mpi_error_unit + process_id

#ifdef _WIN
  OUTDIR = '#output\'
#else
  OUTDIR = '#output/'
#endif
  if( process_id == master_id )then
#ifdef GNU
    RES = SYSTEM( 'mkdir -p ./' // trim(OUTDIR))
#else       
    RES = MAKEDIRQQ('#output')
    if( .not. RES )then
      RES = GETLASTERRORQQ( )
      if( RES == ERR$NOENT )then
        call STOPP('THE PATH FOR $OUTPUT IS NOT FOUND!')
      elseif( RES == ERR$ACCES )then
        call STOPP('CANNOT CREATE THE FOLDER: PERMISSION DENIED!')
      endif
    endif
#endif    
  endif
  call MPI_barrier(MPI_Comm_World, ierr)

  write(STR, '("*** EFDC+ Filename: ",A,",  Version: ",A10," ***")') TRIM(EFDC_EXE), EFDC_VER

  ! *** Create output file for each processor to record EFDC+ messages
  write(mpi_filename, '(A13,I3.3,A4)')    'EFDC_out_proc_', process_id, '.log'
  open(unit = mpi_efdc_out_unit, status = 'replace', file = OUTDIR//mpi_filename)
  mpi_efdc_out_file = mpi_filename
  write(mpi_efdc_out_unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))

  ! *** Create output file for each processor to record EFDC+ mpi messages
  write(mpi_filename, '(A13,I3.3,A4)')    'log_mpi_proc_', process_id, '.log'
  open(unit = mpi_log_unit, status = 'replace', file = OUTDIR//mpi_filename)
  mpi_log_file = mpi_filename
  write(mpi_log_unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
  
  ! *** Create output file for each processor to record EFDC+ errors
  write(mpi_filename, '(A15,I3.3,A4)')    'log_error_proc_', process_id, '.log'
  open(unit = mpi_error_unit, status = 'replace', file = OUTDIR//mpi_filename)
  mpi_error_file = mpi_filename
  write(mpi_error_unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
  
  ! *** Write the file out
  if( MPI_DEBUG_FLAG )then
    write(mpi_filename, '(A14,I3.3,A4)')  'comm_mpi_proc_', process_id, '.log'
    open(unit = mpi_comm_unit, status = 'replace', file = OUTDIR//mpi_filename)
    write(mpi_comm_unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
  endif

  ! *** Write the file out
  write(mpi_filename, '(A13,I3.3,A4)')    'map_mpi_proc_', process_id, '.log'
  open(unit = mpi_mapping_unit, status = 'replace', file = OUTDIR//mpi_filename)
  write(mpi_mapping_unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))

  End subroutine Setup_MPI_Debug_File

  