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
  ! @details Creates data structure that contains local i/j corresponding to
  !! global i/j and the process that cell resides in.  This gets communicated to
  !! all processes in this routine
  ! @author Zander Mausolff
  ! @date 12/4/2019
  !---------------------------------------------------------------------------!

  Subroutine Setup_Local_to_Global

  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Mapping
  use MPI

  implicit none

  integer :: ll, ii, i_tmp, j_tmp, i, j
  Integer(4), Allocatable, Dimension(:) :: recv_size_all !< Holds sizes each process will send for gathering IG2IL
  Integer(4), Allocatable, Dimension(:) :: disps_all !< Holds displacements array for gathering IG2IL onto all processes
  Integer(4) :: start
  Integer(4) :: send_size
  Integer(4) :: displacement
  integer(4) :: IERR     !< local MPI error flag
  integer :: l_sort
  !---------------------------------------------------------------------------!
  allocate(recv_size_all(num_Processors))
  allocate(disps_all(num_Processors))

  allocate(loc_to_glob((ic-4)*(jc-4)))
  loc_to_glob(:).process  =  0
  loc_to_glob(:).local_i  =  0
  loc_to_glob(:).local_j  =  0
  loc_to_glob(:).global_i =  0
  loc_to_glob(:).global_j =  0
  loc_to_glob(:).local_l  =  0

  allocate(all_loc_to_glob(ic_global*jc_global))
  all_loc_to_glob(:).process  = 0
  all_loc_to_glob(:).local_i  = 0
  all_loc_to_glob(:).local_j  = 0
  all_loc_to_glob(:).global_i = 0
  all_loc_to_glob(:).global_j = 0
  all_loc_to_glob(:).local_l  = 0
  all_loc_to_glob(:).global_l = 0

  allocate(sorted_loc_to_glob(LA_Global))
  sorted_loc_to_glob(:).process  = 0
  sorted_loc_to_glob(:).local_i  = 0
  sorted_loc_to_glob(:).local_j  = 0
  sorted_loc_to_glob(:).global_i = 0
  sorted_loc_to_glob(:).global_j = 0
  sorted_loc_to_glob(:).local_l  = 0
  sorted_loc_to_glob(:).global_l = 0

  if( active_domains == 1 )then
    ! *** Assign all_loc_to_glob data so the L indices are contiguous with respect to their global values
    Do i = 1, LA_Global-1
      l_sort = i + 1
      loc_to_glob(i).process = 0
      
      loc_to_glob(i).local_l = l_sort
      loc_to_glob(i).global_l = l_sort
      loc_to_glob(i).local_i = IL(l_sort)
      loc_to_glob(i).local_j = JL(l_sort)
      loc_to_glob(i).global_i = IL(l_sort)
      loc_to_glob(i).global_j = JL(l_sort)
      
      sorted_loc_to_glob(l_sort).process  = loc_to_glob(i).process
      sorted_loc_to_glob(l_sort).local_i  = loc_to_glob(i).local_i
      sorted_loc_to_glob(l_sort).local_j  = loc_to_glob(i).local_j
      sorted_loc_to_glob(l_sort).global_i = loc_to_glob(i).global_i
      sorted_loc_to_glob(l_sort).global_j = loc_to_glob(i).global_j
      sorted_loc_to_glob(l_sort).local_l  = loc_to_glob(i).local_l
      sorted_loc_to_glob(l_sort).global_l = loc_to_glob(i).global_l
    enddo
    i = 0 ! *** debugging only
    return
    
  endif
  
  ! *** Create local mapping from i/j to global i/j and tag the process that cell belongs to
  start = 0
  Do j = 1,jc_global
    Do i = 1,ic_global
      ! *** If we are in the subdomain select the corresponding global i
      if( IG2IL(i) > 2 .and. IG2IL(i) <= ic-2 )then                ! *** exclude ghost cells
        ! *** If we are in the subdomain select the corresponding global j
        if( JG2JL(j) > 2 .and. JG2JL(j) <= jc-2 )then              ! *** exclude ghost cells
          if( ijct_global(i,j) > 0 .and. ijct_global(i,j) < 9 )then
            start = start + 1
            ! *** ensure each i/j is mapped to a process
            loc_to_glob(start).process = process_id
            ! *** Local i/j
            loc_to_glob(start).local_i = IG2IL(i)
            loc_to_glob(start).local_j = JG2JL(j)
            ! *** Global i/j
            loc_to_glob(start).global_i = i
            loc_to_glob(start).global_j = j
            loc_to_glob(start).local_l  = lij(IG2IL(i), JG2JL(j))
            loc_to_glob(start).global_l = lij_global(i,j)
          endif
        endif
      endif
    Enddo
  Enddo

  ! *** First, need to make sure each process knows the size of all messages being sent
  send_size = start

  call MPI_Allgather(send_size, 1, MPI_Integer4, recv_size_all, 1, MPI_Integer4, comm_2d, ierr)

  ! *** Need to determine displacements for the global array about to be gathered
  disps_all(1) = 0
  do ii = 2, active_domains
    disps_all(ii) = disps_all(ii-1) + recv_size_all(ii-1)
  enddo
    
#ifdef GNU  
  ! *** Gather all onto single global data structure.  Doing it piecemeal to avoid setting up MPI data structure
  call MPI_AllGatherV(loc_to_glob.process, send_size, MPI_Integer4, &
    all_loc_to_glob.process, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
  ! *** Gather local I index
  call MPI_AllGatherV(loc_to_glob.local_i, send_size, MPI_Integer4, &
    all_loc_to_glob.local_i, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
  ! *** Gather local J index
  call MPI_AllGatherV(loc_to_glob.local_j, send_size, MPI_Integer4, &
    all_loc_to_glob.local_j, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
  ! *** Gather global I index
  call MPI_AllGatherV(loc_to_glob.global_i, send_size, MPI_Integer4, &
    all_loc_to_glob.global_i, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
  ! *** Gather global J index
  call MPI_AllGatherV(loc_to_glob.global_j, send_size, MPI_Integer4, &
    all_loc_to_glob.global_j, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
  ! *** Gather local L index
  call MPI_AllGatherV(loc_to_glob.local_l, send_size, MPI_Integer4, &
    all_loc_to_glob.local_l, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
  ! *** Gather global L index
  call MPI_AllGatherV(loc_to_glob.global_l, send_size, MPI_Integer4, &
    all_loc_to_glob.global_l, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, ierr)
#else
  ! *** Gather all onto single global data structure.  Doing it piecemeal to avoid setting up MPI data structure
  call MPI_AllGatherV(loc_to_glob.process, send_size, MPI_Integer4, &
    all_loc_to_glob.process, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
  ! *** Gather local I index
  call MPI_AllGatherV(loc_to_glob.local_i, send_size, MPI_Integer4, &
    all_loc_to_glob.local_i, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
  ! *** Gather local J index
  call MPI_AllGatherV(loc_to_glob.local_j, send_size, MPI_Integer4, &
    all_loc_to_glob.local_j, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
  ! *** Gather global I index
  call MPI_AllGatherV(loc_to_glob.global_i, send_size, MPI_Integer4, &
    all_loc_to_glob.global_i, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
  ! *** Gather global J index
  call MPI_AllGatherV(loc_to_glob.global_j, send_size, MPI_Integer4, &
    all_loc_to_glob.global_j, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
  ! *** Gather local L index
  call MPI_AllGatherV(loc_to_glob.local_l, send_size, MPI_Integer4, &
    all_loc_to_glob.local_l, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
  ! *** Gather global L index
  call MPI_AllGatherV(loc_to_glob.global_l, send_size, MPI_Integer4, &
    all_loc_to_glob.global_l, recv_size_all, disps_all, &
    MPI_Integer4, comm_2d, master_id, ierr)
#endif   
  !---------------------------------------------------------------------------!

  ! *** Sort the all_loc_to_glob data so the L indices are contiguous with respect to their global values
  Do i = 1, LA_Global-1
    ! *** Get the 'correct' global l
    l_sort = all_loc_to_glob(i).global_l
    ! *** Assign the global l value to a new array
    sorted_loc_to_glob(l_sort).process  = all_loc_to_glob(i).process
    sorted_loc_to_glob(l_sort).local_i  = all_loc_to_glob(i).local_i
    sorted_loc_to_glob(l_sort).local_j  = all_loc_to_glob(i).local_j
    sorted_loc_to_glob(l_sort).global_i = all_loc_to_glob(i).global_i
    sorted_loc_to_glob(l_sort).global_j = all_loc_to_glob(i).global_j
    sorted_loc_to_glob(l_sort).local_l  = all_loc_to_glob(i).local_l
    sorted_loc_to_glob(l_sort).global_l = all_loc_to_glob(i).global_l
  enddo


#ifdef DEBUGGING
  ! *** Global
  call WriteBreak(mpi_log_unit)
  if( MPI_Write_Flag )then
    write(mpi_log_unit,'(a)') 'Global mapping of ALL active cells and which process each cell belongs in'
    write(mpi_log_unit,'(a)') ' process,global l, local l, loc i, loc j, gl i, gl j'
    do i = 2, LA_Global
      write(mpi_log_unit,'(i6,a,i6,a,5i6)') sorted_loc_to_glob(i).process,' | ', i,' | ',sorted_loc_to_glob(i).local_l, sorted_loc_to_glob(i).local_i, &
        sorted_loc_to_glob(i).local_j, sorted_loc_to_glob(i).global_i, sorted_loc_to_glob(i).global_j
    enddo
    call WriteBreak(mpi_log_unit)
  endif
  
  write(mpi_log_unit,'(a,10i7)') 'Recv size all ', recv_size_all
  write(mpi_log_unit,'(a,2i7)')  'Sum recv ' , sum(recv_size_all), (ic_global*jc_global)
  write(mpi_log_unit,'(a,10i7)') 'Displacements in send ', disps_all
#endif

  End subroutine Setup_Local_to_Global
