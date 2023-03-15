! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!> @details Reads the JSON formatted DECOMP.inp file to get variables to
!! allocate domain decomposition related arrays
! 
!> @author Zander Mausolff
!> @date 2/6/2020
Subroutine Scan_JSON_Decomp

  Use Global
  Use Variables_MPI
  Use Broadcast_Routines
  Use fson
  Use mod_fson_value, only: fson_value_count, fson_value_get
#ifdef _MPI
  Use MPI
#endif
  Implicit None

  ! *** Local variables
  Integer :: idim !< Temporary value of the width in the x-direction as its being read in from DECOMP.JNP
  Integer :: jdim !< Temporary value of the width in the y-direction as its being read in from DECOMP.JNP
  Integer :: i, j, iso, it, nD, nelements
  Integer :: scan_decomp_unit = 200
  INTEGER(4) :: IERR                                    !< local MPI error flagreading
  Integer,Allocatable :: itmp(:)

  ! *** Local variables for JSON reader
  Type(fson_value), pointer :: json_data, array, item   !< Declare a pointer variables.  Always use a pointer with fson_value.
  character(len=1024) :: strval, strval2
  
  if( process_id == master_id )THEN
    WRITE(6,'(A)')'SCANNING INPUT FILE: DECOMP.JNP'
    
    ! *** Open up the decomp file formatted in json
    json_data => fson_parse("decomp.jnp")

    Call fson_get(json_data, "number_i_subdomains", n_x_partitions)
    Call fson_get(json_data, "number_j_subdomains", n_y_partitions)
    Call fson_get(json_data, "number_active_subdomains", active_domains)

    ! *** Read in arrays
    Call fson_get(json_data, "i_subdomain_widths", ic_decomp)
    Call fson_get(json_data, "j_subdomain_widths", jc_decomp)
    
    ! *** Get the active flags
    Call fson_get(json_data, "active_flag", itmp)               
    nelements = size(itmp)

    ! *** QC
    if( nelements /= n_x_partitions*n_y_partitions )then
      print *, 'Wrong Active Flag count:  nX*nY: ', n_x_partitions*n_y_partitions,', Entered: ', nelements
      stop
    endif
    
    max_width_x = maxval(ic_decomp, 1)
    max_width_y = maxval(jc_decomp, 1)
    max_width_x = max_width_x + 4
    max_width_y = max_width_y + 4

  End if
  
  ! *** Need to do broadcast with communicator MPI_Comm_World because 
  !     the comm_2d communicator has not been initiated yet
  Call MPI_BCAST(n_x_partitions, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)
  Call MPI_BARRIER(MPI_Comm_World, ierr)

  Call MPI_BCAST(n_y_partitions, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)
  Call MPI_BARRIER(MPI_Comm_World, ierr)

  Call MPI_BCAST(active_domains, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)
  Call MPI_BARRIER(MPI_Comm_World, ierr)

  Call MPI_BCAST(max_width_x, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)
  Call MPI_BARRIER(MPI_Comm_World, ierr)

  Call MPI_BCAST(max_width_y, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)
  Call MPI_BARRIER(MPI_Comm_World, ierr)

  Call MPI_BCAST(nelements, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)
  if( process_id /= master_id )then
    Allocate(itmp(nelements))
    itmp = 0
  endif
  Call MPI_BCAST(itmp, nelements, MPI_Integer, master_id, MPI_Comm_World, ierr)
  Call MPI_BARRIER(MPI_Comm_World, ierr)

  ! *** Set the active flags
  Allocate(decomp_active(0:n_x_partitions+1,0:n_y_partitions+1))
  Allocate(process_map(0:n_x_partitions+1,0:n_y_partitions+1))
  decomp_active = 0
  process_map = -1
  
  ! *** Retreive row by row
  it = 0
  nD = 0
  do j = 1,n_y_partitions
    do i = 1,n_x_partitions
      it = it + 1
      decomp_active(i,j) = itmp(it)
      
      if( itmp(it) == 1 )then
        process_map(i,j) = nD                !< Zero base
        nD = nD + 1
      endif
    enddo
  enddo
  
  ! *** Ensure consistency
  if( nD /= active_domains )then
    ! do something
  endif
  
  Return

End subroutine Scan_JSON_Decomp
