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
  ! Module: Variables_MPI
  !
  !> @details This module contains primary variables for MPI
  !
  !> @author Zander Mausolf
  !---------------------------------------------------------------------------!
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !---------------------------------------------------------------------------!
  !    2019-09       Zander     Moved variables into other modules for clarity
  !    2019-06       Zander     Added variables for MPI implementation
  !---------------------------------------------------------------------------!

  Module Variables_MPI

  use GLOBAL

  implicit none
  
  integer :: ic_pos, jc_pos !< starting positions within the global array in each I/J direction

  integer :: nincrmt_reduce
  real    :: dtdyn_reduce

  integer :: NPFOR_Readin

  integer, Allocatable, Dimension(:) :: recv_counts_3d !< Keeps track of the size of each processor is sending to the master for (L,K)
  integer, Allocatable, Dimension(:) :: recv_counts_1d !< Keeps track of the message size each process is sending (L) arrays

  integer :: stride
  integer :: send_size_1d !< size of the message for communicating the entire domain for (L) arrays


  integer :: LA_global  !< Global value of LA - the total number of cells
  integer :: LCM_Global !< Global LCM value as set in the SCANEFDC subroutine.
  integer :: ICM_Global !< gets original value of ICM as read in from file
  integer :: JCM_Global !< gets original value of JCM as read in from file

  integer :: End_global_LA   !< Ending global L index for a given subdomain
  integer :: Start_Global_LA !< Starting global L index for a given subdomain

  integer, Allocatable, Dimension(:) :: i_end_ghost_in_each !< keeps track of subdomains with ghost cells at end(east) side
  integer, Allocatable, Dimension(:) :: j_end_ghost_in_each
  integer, Allocatable, Dimension(:) :: i_start_ghost_in_each
  integer, Allocatable, Dimension(:) :: j_start_ghost_in_each

  character(12) :: filename_out

  integer :: process_id                 !< ID of processor
  integer :: num_Processors             !< Total number of processors
  integer :: master_id = 0              !< Defines the master node/process ID
                                        
  Character(24) :: mpi_log_file         !< General log file name for each processor
  Character(24) :: mpi_qdwaste          !< QDWASTE log file name for each processor
  Character(24) :: mpi_efdc_out_file    !< EFDC+ log file name for each processor
  Character(24) :: mpi_error_file       !< EFDC+ ERROR log file name for each processor
  
  ! *** Unit numbers for MPI log files, will be = (unit num + process_id)
  integer :: mpi_efdc_out_unit = 777    !< Unit so each process writes out its own EFDCLOG.OUT
  integer :: mpi_log_unit = 100         !< General log unit
  integer :: mpi_qdwaste_unit = 200     !< QDWASTE log unit
  integer :: mpi_comm_unit = 300        !< MPI setup debug file for the communication of ghost cell values
  integer :: mpi_mapping_unit = 500     !< Log of global to local cell mapping per process
  integer :: mpi_error_unit = 900       !< Log of error messages per process

  logical :: MPI_DEBUG_FLAG = .FALSE.   !< Boolean that turns on some print statements for debugging MPI code. This flag should be an input option at some point.
  logical :: MPI_Write_Flag = .FALSE.   

  integer, Allocatable, Dimension(:) :: displacements_L_index     !< Keeps track of displacements_L_index (in terms of MPI derived type) for ScatterV/GatherV processing
  integer, Allocatable, Dimension(:) :: recv_counts_array !< Counts the number of derived types each process recives

  ! *** For MPI Topology routines
  integer, dimension(2) :: domain_coords  !< Assumes 2D decomposition
  integer, dimension(2)  :: dimensions    !< Keeps track of dimensions of the problem
  logical, dimension(2)  :: is_periodic   !< Tells you if one of the directions is periodic
  logical :: reorder                      !< Boolean indicating if MPI should reorder and find a good way to assign
                                          !<       the process to the elements of the decomposition
  integer :: comm_2d                      !< New communicator for the MPI topology routines

  ! *** For new drifter routines
  integer :: NPD_TOT

  ! *** DECOMP.inp related variables
  integer :: n_ghost_rows = 2             !< Number of ghost rows/cols on each side of domain
  integer :: decomp_input_unit  = 123     !< unit for writing to DECOMP.inp
  integer :: decomp_output_unit = 124     !< unit for writing to DECOMP.out
                                          
  integer :: n_x_partitions               !< number of partitions in the x-direction
  integer :: n_y_partitions               !< number of partitions in the y-direction
  integer :: active_domains               !< total number of partitions in the domain

  integer, Allocatable, Dimension(:)   :: ic_decomp              !< Keeps track of the partition width  in the x-direction
  integer, Allocatable, Dimension(:)   :: jc_decomp              !< Keeps track of the partition height in the y-direction
  integer, Allocatable, Dimension(:,:) :: decomp_active          !< ic_decomp by jc_decomp array of flags, 0-Inactive, 1-Active
  integer, Allocatable, Dimension(:,:) :: process_map            !< ic_decomp by jc_decomp array of process_id's

  integer, Allocatable, Dimension(:) :: IB_Decomp, IE_Decomp     !< Starting and Ending global I for each partition in the x-direction
  integer, Allocatable, Dimension(:) :: JB_Decomp, JE_Decomp     !< Starting and Ending global J for each partition in the y-direction

  integer :: max_width_x                                         ! <== PNX
  integer :: max_width_y                                         ! <== PNY

  integer :: ic_global                                           !< IC value for the entire domain, this is what is originally read in from efdc.inp
  integer :: jc_global                                           !< JC value for the entire domain, this is what is originally read in from efdc.inp
  integer :: lc_global                                           !< LC value for the entire domain, this is what is originally read in from efdc.inp

  integer :: x_id
  integer :: y_id

  integer, Allocatable, Dimension(:) :: IL2IG
  integer, Allocatable, Dimension(:) :: JL2JG

  integer, Allocatable, Dimension(:) :: IG2IL
  integer, Allocatable, Dimension(:) :: JG2JL

  integer :: global_max_width_x ! <== GNX
  integer :: global_max_width_y ! <== GNY

  real(RKD), Allocatable, Dimension(:,:) :: XCOR_Global          !< 
  real(RKD), Allocatable, Dimension(:,:) :: YCOR_Global          !< 
  real(RKD), Allocatable, Dimension(:)   :: Area_Global          !< 
  
  integer, Allocatable, Dimension(:,:,:) :: LWDIR_Global         !< L index list for cells along each fetch for all cells in the global domain
  integer, Allocatable, Dimension(:)     :: UMASK_Global         !< Flag for U face mask on
  integer, Allocatable, Dimension(:)     :: VMASK_Global         !< Flag for V face mask on 
  real(RKD), Allocatable, Dimension(:,:) :: FWDIR_Global         !< Fetch length for wind sector midpoint for each wind direction
  
  ! *** Create type that helps in the global remapping
  Type mapping_lij
    integer :: PR      ! *** MPI process ID
    integer :: IL      ! *** Local I
    integer :: JL      ! *** Local J
    integer :: LL      ! *** Local L
    integer :: IG      ! *** Global I
    integer :: JG      ! *** Global J
    integer :: LG      ! *** Global L
  End Type

  type(mapping_lij), Allocatable, Dimension(:) :: Map2Global    !< Directly map L_local  to L_global
  type(mapping_lij), Allocatable, Dimension(:) :: Map2Local     !< Directly map L_global to L_local

  integer, Allocatable, Dimension(:,:) :: LIJ_Global            !< Keeps track of global active cell index based on I,J mapping
  integer, Allocatable, Dimension(:,:) :: IJCT_GLOBAL
  integer, Allocatable, Dimension(:,:) :: IJCTLT_GLOBAL
  
  !             nbr_north
  !          |-------------|
  !          |             |
  !          |             |
  ! nbr_west | process_id  |  nbr_east
  !          |             |
  !          |             |
  !          |-------------|
  !             nbr_south

  integer :: nbr_west  !< Index that keeps track of which process is the neighbor West of the local domain
  integer :: nbr_east  !< Index that keeps track of which process is the neighbor East of the local domain
  integer :: nbr_south !< Index that keeps track of which process is the neighbor South of the local domain
  integer :: nbr_north !< Index that keeps track of which process is the neighbor North of the local domain
  integer, Allocatable, Dimension(:,:,:) :: Comm_Cells     !< Sub-Domain active cell interface list
  integer, Allocatable, Dimension(:,:)   :: nComm_Cells    !< Sub-Domain active cell interface count

  integer :: lmap
  integer :: size_mpi  !< Number of MPI sub-domains a single domain is to write debugging information

  !------------------------------------------------------------------
  !***These variables are likely no longer needed
  ! *** For Task Division Routine
  integer, Allocatable, Dimension(:) :: num_task_per_process     !< The total nuber of tasks each process will take
  integer, Allocatable, Dimension(:) :: lower_task_bound         !< Array that keeps track of each processors lower bound constituent number it will calculate
  integer, Allocatable, Dimension(:) :: upper_task_bound         !< Array that keeps track of each processors upper bound constituent number it will calculate
  integer, Allocatable, Dimension(:) :: num_elem_to_each_array   !< Number of elements in each array
  integer :: new_three_dim_data_type                             !< New data type representing a contigous block of memory from a 3d array

  ! *** Timing variables
  Double Precision, Save :: start_time, end_time, total_time_calconc
  Double Precision, Save :: start_time_bcast,   end_time_bcast,   total_time_bcast
  Double Precision, Save :: start_time_scatter, end_time_scatter, total_time_scatter
  Double Precision, Save :: start_time_gather,  end_time_gather,  total_time_gather
  Double Precision, Save :: start_time_caltran, end_time_caltran, total_time_caltran
  !***End no longer needed
  !------------------------------------------------------------------
  
  logical :: DoublePrecision

  End module Variables_MPI
