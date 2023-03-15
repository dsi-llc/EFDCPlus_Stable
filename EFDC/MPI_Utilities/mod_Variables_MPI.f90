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

  Use Global

  Implicit None
  
  Integer :: ic_pos, jc_pos !< starting positions within the global array in each I/J direction

  Integer :: nincrmt_reduce
  Real    :: dtdyn_reduce

  Integer :: NPFOR_Readin

  Integer, Allocatable, Dimension(:) :: recv_counts_3d !< Keeps track of the size of each processor is sending to the master for (L,K)
  Integer, Allocatable, Dimension(:) :: recv_counts_1d !< Keeps track of the message size each process is sending (L) arrays

  Integer :: stride
  Integer :: send_size_1d !< size of the message for communicating the entire domain for (L) arrays


  Integer :: LA_global  !< Global value of LA - the total number of cells
  Integer :: LCM_Global !< Global LCM value as set in the SCANEFDC subroutine.
  Integer :: ICM_Global !< gets original value of ICM as read in from file
  Integer :: JCM_Global !< gets original value of JCM as read in from file

  Integer :: End_global_LA   !< Ending global L index for a given subdomain
  Integer :: Start_Global_LA !< Starting global L index for a given subdomain

  Integer, Allocatable, Dimension(:) :: i_end_ghost_in_each !< keeps track of subdomains with ghost cells at end(east) side
  Integer, Allocatable, Dimension(:) :: j_end_ghost_in_each
  Integer, Allocatable, Dimension(:) :: i_start_ghost_in_each
  Integer, Allocatable, Dimension(:) :: j_start_ghost_in_each

  Integer :: unit_efdc_out = 777 !< unit so each process writes out its own EFDC.OUT
  character(12) :: filename_out

  Integer :: process_id                 !< ID of processor
  Integer :: num_Processors             !< Total number of processors
  Integer :: master_id = 0              !< Defines the master node/process ID
                                        
  Character(22) :: mpi_log_file         !< Output log for MPI calculation
  Integer :: mpi_log_unit = 100         !< File unit that each processor can write to, will be = (unit num + process_id)
  Integer :: mpi_comm_unit = 300        !< Writes out the communication of ghost cell values
  Integer :: mpi_mapping_unit = 500     !< Log of global to local cell mapping per process
  Integer :: mpi_err_unit = 900         !< Log of error messages per process
  Logical :: MPI_DEBUG_FLAG = .FALSE.   !< Boolean that turns on some print statements for debugging MPI code. This flag should be an input option at some point.
  LOGICAL :: MPI_Write_Flag = .FALSE.   
  Integer :: MPI_Par_Flag = 0           !< Flag that indicates whether we are doing a MPI calculation or not

  Integer, Allocatable, Dimension(:) :: displacements_L_index     !< Keeps track of displacements_L_index (in terms of MPI derived type) for ScatterV/GatherV processing
  Integer, Allocatable, Dimension(:) :: recv_counts_array !< Counts the number of derived types each process recives

  ! *** For MPI Topology routines
  Integer, dimension(2) :: domain_coords  !< Assumes 2D decomposition
  Integer, dimension(2)  :: dimensions    !< keeps track of dimensions of the problem
  Logical, dimension(2)  :: is_periodic   !< tells you if one of the directions is periodic
  Logical :: reorder                      !< Boolean indicating if MPI should reorder and find a good way to assign
                                          !<       the process to the elements of the decomposition
  Integer :: comm_2d                      !< New communicator for the MPI topology routines

  ! *** For new drifter routines
  Integer :: NPD_TOT

  ! *** DECOMP.inp related variables
  Integer :: n_ghost_rows = 2             !< Number of ghost rows/cols on each side of domain
  Integer :: decomp_input_unit  = 123     !< unit for writing to DECOMP.inp
  Integer :: decomp_output_unit = 124     !< unit for writing to DECOMP.out
                                          
  Integer :: n_x_partitions               !< number of partitions in the x-direction
  Integer :: n_y_partitions               !< number of partitions in the y-direction
  Integer :: active_domains               !< total number of partitions in the domain

  Integer, Allocatable, Dimension(:)   :: ic_decomp              !< Keeps track of the partition width  in the x-direction
  Integer, Allocatable, Dimension(:)   :: jc_decomp              !< Keeps track of the partition height in the y-direction
  Integer, Allocatable, Dimension(:,:) :: decomp_active          !< ic_decomp by jc_decomp array of flags, 0-Inactive, 1-Active
  Integer, Allocatable, Dimension(:,:) :: process_map            !< ic_decomp by jc_decomp array of process_id's

  Integer, Allocatable, Dimension(:) :: IB_Decomp, IE_Decomp     !< Starting and Ending global I for each partition in the x-direction
  Integer, Allocatable, Dimension(:) :: JB_Decomp, JE_Decomp     !< Starting and Ending global J for each partition in the y-direction

  Integer :: max_width_x                                         ! <== PNX
  Integer :: max_width_y                                         ! <== PNY

  Integer :: ic_global                                           !< IC value for the entire domain, this is what is originally read in from efdc.inp
  Integer :: jc_global                                           !< JC value for the entire domain, this is what is originally read in from efdc.inp
  Integer :: lc_global                                           !< LC value for the entire domain, this is what is originally read in from efdc.inp

  Integer :: x_id
  Integer :: y_id

  Integer, Allocatable, Dimension(:) :: IL2IG
  Integer, Allocatable, Dimension(:) :: JL2JG

  Integer, Allocatable, Dimension(:) :: IG2IL
  Integer, Allocatable, Dimension(:) :: JG2JL

  Integer :: global_max_width_x ! <== GNX
  Integer :: global_max_width_y ! <== GNY

  Real(RKD), Allocatable, Dimension(:,:) :: XCOR_Global
  Real(RKD), Allocatable, Dimension(:,:) :: YCOR_Global
  Real(RKD), Allocatable, Dimension(:)   :: Area_Global
  Real(RKD), Allocatable, Dimension(:,:) :: FWDIR_Global
  Real(RKD), Allocatable, Dimension(:)   :: UMASK_Global
  Real(RKD), Allocatable, Dimension(:)   :: VMASK_Global
  
  ! *** Create type that helps in the global remapping
  Type mapping_lij
    Integer :: PR      ! *** MPI process ID
    Integer :: IL      ! *** Local I
    Integer :: JL      ! *** Local J
    Integer :: LL      ! *** Local L
    Integer :: IG      ! *** Global I
    Integer :: JG      ! *** Global J
    Integer :: LG      ! *** Global L
  End Type

  Type(mapping_lij), Allocatable, Dimension(:) :: Map2Global    !< Directly map L_local  to L_global
  Type(mapping_lij), Allocatable, Dimension(:) :: Map2Local     !< Directly map L_global to L_local

  Integer, Allocatable, Dimension(:,:) :: LIJ_Global            !< Keeps track of global active cell index based on I,J mapping
  Integer, Allocatable, Dimension(:,:) :: IJCT_GLOBAL
  Integer, Allocatable, Dimension(:,:) :: IJCTLT_GLOBAL
  
  !             nbr_north
  !          |-------------|
  !          |             |
  !          |             |
  ! nbr_west | process_id  |  nbr_east
  !          |             |
  !          |             |
  !          |-------------|
  !             nbr_south

  Integer :: nbr_west  !< Index that keeps track of which process is the neighbor West of the local domain
  Integer :: nbr_east  !< Index that keeps track of which process is the neighbor East of the local domain
  Integer :: nbr_south !< Index that keeps track of which process is the neighbor South of the local domain
  Integer :: nbr_north !< Index that keeps track of which process is the neighbor North of the local domain
  Integer, Allocatable, Dimension(:,:,:) :: Comm_Cells     !< Sub-Domain active cell interface list
  Integer, Allocatable, Dimension(:,:)   :: nComm_Cells    !< Sub-Domain active cell interface count

  Integer :: lmap
  Integer :: size_mpi  !< Number of MPI sub-domains a single domain is to write debugging information

  !------------------------------------------------------------------
  !***These variables are likely no longer needed
  ! *** For Task Division Routine
  Integer, Allocatable, Dimension(:) :: num_task_per_process     !< The total nuber of tasks each process will take
  Integer, Allocatable, Dimension(:) :: lower_task_bound         !< Array that keeps track of each processors lower bound constituent number it will calculate
  Integer, Allocatable, Dimension(:) :: upper_task_bound         !< Array that keeps track of each processors upper bound constituent number it will calculate
  Integer, Allocatable, Dimension(:) :: num_elem_to_each_array   !< Number of elements in each array
  Integer :: new_three_dim_data_type                             !< New data type representing a contigous block of memory from a 3d array

  ! *** Timing variables
  Double Precision, Save :: start_time, end_time, total_time_calconc
  Double Precision, Save :: start_time_bcast,   end_time_bcast,   total_time_bcast
  Double Precision, Save :: start_time_scatter, end_time_scatter, total_time_scatter
  Double Precision, Save :: start_time_gather,  end_time_gather,  total_time_gather
  Double Precision, Save :: start_time_caltran, end_time_caltran, total_time_caltran
  !***End no longer needed
  !------------------------------------------------------------------
  
  LOGICAL :: DoublePrecision

  End module Variables_MPI
