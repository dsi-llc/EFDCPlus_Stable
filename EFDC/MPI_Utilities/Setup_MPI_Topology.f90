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
! @details Sets up communicator and MPI Topology routines.
!! The MPI topology routines setup a 'grid' for the processors
!! This grid allows each processor to know which other processors are nearby:
!! Example:
!!  |-----|-----|
!!  | P2  | P3  |
!!  |-----|-----|
!!  | P0  | P1  |
!!  |-----|-----|
!! The routines are intrinsic to MPI see
!! http://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-topo.html
! @author Zander Mausolff
! @date 8/14/2019
!---------------------------------------------------------------------------!

Subroutine Setup_MPI_Topology

  Use MPI
  Use GLOBAL
  Use Variables_MPI

  Implicit none

  ! *** Local
  Integer :: ierr, maxsize, nnodes, i, j, i1, edge, inode, iedge
  Integer,Allocatable,dimension(:) :: index
  Integer,Allocatable,dimension(:) :: NearestNeighbor
  Integer,Allocatable,dimension(:) :: edges
  
  Call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)') 'Setting up MPI Topology'
  write(mpi_log_unit, '(a)') 'This configures where process is in reference to one another'
  write(mpi_log_unit, '(a)') 'For example, with 4 processes the layout will look something like this'
  write(mpi_log_unit, '(a)') '|-----|-----|        |-----|-----|'
  write(mpi_log_unit, '(a)') '| P2  | P3  |        |(0,1)|(1,1)|'
  write(mpi_log_unit, '(a)') '|-----|-----| -----> |-----|-----|'
  write(mpi_log_unit, '(a)') '| P0  | P1  |        |(0,0)|(1,0)|'
  write(mpi_log_unit, '(a)') '|-----|-----|        |-----|-----|'
  write(mpi_log_unit, '(a)') ' Process IDs      Virtual Coordinates'

  ndim = 2 ! 2D topology is set
  is_periodic = .FALSE. !< no periodic boundary condition
  reorder     = .FALSE. !< Do not allow MPI to reoreder the location of processes for additional optimizations

  ! *** Set the number of partitions in each of the (2) directions
  dimensions(1) = n_x_partitions
  dimensions(2) = n_y_partitions
  
  ! *** Topology Option
  if( .false. )then
    ! *** Cartesian Topology
    ! *** MPI Tolopolgy routine, creates new communictor for 2d decomp
    Call MPI_Cart_Create(MPI_Comm_World, ndim, dimensions, is_periodic, reorder, comm_2d, ierr )

    ! *** Get rank of the process within the 'new' communicator comm_2d
    Call MPI_Comm_Rank(comm_2d,   process_id, ierr)  ! process_id becomes 'new id'

    ! *** Gets the location of each processes reference to one another
    Call MPI_Cart_Coords(comm_2d, process_id, 2, domain_coords)

    ! *** Get the East and West neighbors in addition to the north and south neighbors
    Call MPI_Cart_Shift(comm_2d, 0, 1, nbr_west,  nbr_east,  ierr)
    Call MPI_Cart_Shift(comm_2d, 1, 1, nbr_south, nbr_north, ierr)
  else
    ! *** Graph Topology
    ! *** Initial Graph Topology Arrays
    maxsize = n_x_partitions * n_y_partitions
    Allocate(index(maxsize))
    Allocate(NearestNeighbor(maxsize))
    Allocate(edges(maxsize*4))
    index = 0
    NearestNeighbor = 0
    edges = -1
    nnodes = active_domains
    
    ! *** Get Index List
    inode = 0
    iedge = 0
    do j = 1,n_y_partitions
      do i = 1,n_x_partitions
        if( decomp_active(i,j) == 1 )then
          ! *** Active domain.  Determine neghbors
          inode = inode + 1

          ! *** Count Nearest Neighbors
          if( decomp_active(i-1,j) == 1 )then                       ! *** West
            NearestNeighbor(inode) = NearestNeighbor(inode) + 1
            iedge = iedge + 1
            edges(iedge) = process_map(i-1,j)
          endif
          if( decomp_active(i+1,j) == 1 )then                       ! *** East
            NearestNeighbor(inode) = NearestNeighbor(inode) + 1
            iedge = iedge + 1
            edges(iedge) = process_map(i+1,j)
          endif
          if( decomp_active(i,j-1) == 1 )then                       ! *** South
            NearestNeighbor(inode) = NearestNeighbor(inode) + 1
            iedge = iedge + 1
            edges(iedge) = process_map(i,j-1)
          endif
          if( decomp_active(i,j+1) == 1 )then                       ! *** North
            NearestNeighbor(inode) = NearestNeighbor(inode) + 1
            iedge = iedge + 1
            edges(iedge) = process_map(i,j+1)
          endif
        endif
      enddo
    enddo
    
    ! *** Set the index array based on NearestNeighbor
    index(1) = NearestNeighbor(1)
    do i = 2,nnodes
      index(i) = index(i-1) + NearestNeighbor(i)
    enddo
    
    ! *** Create the Topology
    Call MPI_Graph_Create(MPI_Comm_World, nnodes, index, edges, .false., comm_2d, ierr)
    
    ! *** Get rank of the process within the 'new' communicator comm_2d
    Call MPI_Comm_Rank(comm_2d, process_id, ierr) 

    ! *** Initialize to Unconnected
    nbr_west  = -1
    nbr_east  = -1
    nbr_north = -1
    nbr_south = -1
    
    ! *** Determine Connections based on current rank, i.e. process_id
    do j = 1,n_y_partitions
      do i = 1,n_x_partitions
        if( process_map(i,j) == process_id )then
          ! *** Active domain.  Determine neghbors
          nbr_west  = process_map(i-1,j)
          nbr_east  = process_map(i+1,j)
          nbr_south  = process_map(i,j-1)
          nbr_north  = process_map(i,j+1)
          exit
        endif
      enddo
    enddo
  endif
  
  ! *** Set sub-domain coordinates
  domain_coords(1) = -1
  do j = 1,n_y_partitions
    do i = 1,n_x_partitions
      if( process_id == process_map(i,j) )then
        domain_coords(1) = i-1
        domain_coords(2) = j-1
        exit
      endif
    enddo
    if( domain_coords(1) /= -1 ) exit
  enddo
  
  write(mpi_log_unit, '(a,I3,a,i3,a)') 'Virtual Coordinates: (',domain_coords(1),',',domain_coords(2), ')'
  Call WriteInteger(dimensions(1), mpi_log_unit, '# I Partitions:     ')
  Call WriteInteger(dimensions(2), mpi_log_unit, '# J Partitions:     ')
  Call WriteInteger(nbr_west,  mpi_log_unit,  'West  neighbor:     ')
  Call WriteInteger(nbr_east,  mpi_log_unit,  'East  neighbor:     ')
  Call WriteInteger(nbr_north, mpi_log_unit,  'North neighbor:     ')
  Call WriteInteger(nbr_south, mpi_log_unit,  'South neighbor:     ')
  write(mpi_log_unit, '(a)') 'Note, (-1) values indicate there is no neighbor in the stated direction'
  write(mpi_log_unit, '(a)') 'Ending setup of MPI Topology'
  Call WriteBreak(mpi_log_unit)

  Call MPI_Barrier(comm_2d, ierr)

End subroutine Setup_MPI_Topology
