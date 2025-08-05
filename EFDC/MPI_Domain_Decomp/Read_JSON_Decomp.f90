! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  ! @details Reads the DECOMP.JNP file.  This file contains info
  !! on the width in the x-y directions should be set for the decomposition
  ! @author Zander Mausolff - Adopted from read_lorp.f90 from O'Donncha's repo
  ! @date 8/13/2019
  !
  Subroutine Read_JSON_Decomp

  use GLOBAL
  use Variables_MPI
  use Broadcast_Routines
  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  use MPI

  implicit none

  ! *** Local variables
  integer :: n_skip
  integer :: iso, i, j, it
  integer :: max_partition_size_write_out
  integer :: IERR                                      !< local MPI error flag
  
  ! *** Local variables for JSON reader
  type(fson_value), pointer :: json_data   !< Declare a pointer variables.  Always use a pointer with fson_value.
  
  character(len = 1024) :: strval, strval2
  
  ! *** Do reading only on the master process
  if( process_id == master_id )then
    !WRITE(6,'(A)') 'Reading JSON formatted DECOMP.JNP'
    
    ! *** Open up the decomp file formatted in json
    json_data => fson_parse("decomp.jnp")
        
    call fson_get(json_data, "number_i_subdomains", n_x_partitions)
    call fson_get(json_data, "number_j_subdomains", n_y_partitions)
    call fson_get(json_data, "number_active_subdomains", active_domains)

    ! *** Read in arrays
    call fson_get(json_data, "i_subdomain_widths", ic_decomp)
    call fson_get(json_data, "j_subdomain_widths", jc_decomp)
    
    ! *** QC
    if( sum(ic_decomp) /= IC )then
      print *, 'Sum of sub-domain IC''s does not equal global IC: ', sum(ic_decomp), IC
      pause
      stop
    endif
    if( sum(jc_decomp) /= JC )then
      print *, 'Sum of sub-domain JC''s does not equal global JC: ', sum(jc_decomp), JC
      pause
      stop
    endif
    
  endif

  ! *** Broadcast to all processors
  call Broadcast_Scalar(n_x_partitions, master_id)
  call Broadcast_Scalar(n_y_partitions, master_id)
  call Broadcast_Scalar(active_domains, master_id)

  if( process_id /= master_id )then
    ! *** Allocate ic_decomp and jc_decomp for other partitions
    allocate(ic_decomp(n_x_partitions))
    allocate(jc_decomp(n_y_partitions))
  endif
  
  ! *** Broadcast to all
  call Broadcast_Array(ic_decomp, master_id)
  call Broadcast_Array(jc_decomp, master_id)

  ! *** Write out to mpi file
  call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)') 'Reading DECOMP File'
  write(mpi_log_unit, '(a)') 'IC_DECOMP - partition widths in the x-direction: '
  write(mpi_log_unit, 950) (ic_decomp(i), i = 1,n_x_partitions)
  write(mpi_log_unit, '(a)') 'JC_DECOMP - partition heights in the y-direction:'
  write(mpi_log_unit, 950) (jc_decomp(j), j = 1,n_y_partitions)
  write(mpi_log_unit, '(a,I4)') '# of Partitions in X direction = ', n_x_partitions
  write(mpi_log_unit, '(a,I4)') '# of partitions in Y direction = ', n_y_partitions
  call WriteBreak(mpi_log_unit)

 950 Format (100(I5))     ! Formatting out the ic_decomp and jc_decomp arrays, max width of 100...
1000 Format('(a,I4)')

  call MPI_Barrier(DSIcomm, ierr)
 
End subroutine Read_JSON_Decomp
