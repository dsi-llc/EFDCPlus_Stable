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
! @details Global arrays for mapping out boundary conditions and other 
!  mapping values
! @author Zander Mausolff
! @date 11/01/2019
!---------------------------------------------------------------------------!
Module Variables_MPI_Mapping

  implicit none
  
  integer, Allocatable, Dimension(:) :: IL_Global              ! *** Global I index based on Local L index
  integer, Allocatable, Dimension(:) :: JL_Global              ! *** Global J index based on Local L index
  integer, Allocatable, Dimension(:) :: IL_GL                  ! *** Global I index based on Global L index
  integer, Allocatable, Dimension(:) :: JL_GL                  ! *** Global J index based on Global L index
  
  integer :: La_local_no_ghost

  integer :: num_active_l_local
  
  integer,save :: NNoGhost !< Stores the number of ghost cells for a given subdomain
  integer,save,allocatable,dimension(:) :: NoGhost !< Array containing active cells, excluding ghost cells
  logical,save,allocatable,dimension(:) :: IsGhost
  Real,save,allocatable,dimension(:)    :: GhostMask

  integer, Allocatable, Dimension(:) :: lij_west_conn_outside !< Contains mapping of i/j outside of domain to LA values for connection arrays
  integer, Allocatable, Dimension(:) :: lij_east_conn_outside !< Contains mapping of i/j outside of domain to LA values for connection arrays

  integer :: offset_connectors_for_ns

  ! *** Create type that helps in the global remapping
  Type mapping_ij
    integer :: process   !< Process that the cell belongs to 
    integer :: local_i   !< Local I index for that subdomain
    integer :: local_j   !< Local J index for that subdomain
    integer :: global_i  !< Corresponding global I index for the given local I subdomain value
    integer :: global_j  !< Corresponding global J index for the given local I subdomain value
    integer :: local_l   !< Local L index for given subdomain
    integer :: global_l  !< Corresponding global L value for the given local L subdomain value. This is the global L that EE will understand
  End Type

  type(mapping_ij), Allocatable, Dimension(:) :: loc_to_glob !< Local type that contains local i/j and corresponding global i/j
  type(mapping_ij), Allocatable, Dimension(:) :: all_loc_to_glob !< Global type that contains ALL (even inactive cells) local i/j and corresponding global i/j
  type(mapping_ij), Allocatable, Dimension(:) :: active_loc_to_glob !< Maps Global --> ALL local values for only the active cells (globally sized)
  type(mapping_ij), Allocatable, Dimension(:) :: sorted_loc_to_glob !< Sorted array containing local i/j corresponding to global i/j.  These values are indexed off of the global L value

  Real, Allocatable, Dimension(:,:,:) :: GWCSER_Global !< GROUND WATER CONCENTRATION SERIES DATA
  integer, Allocatable, Dimension(:)  :: NGWSL_Global

  ! *** Global periodic boundary conditions
  integer :: NPBW_GL
  integer :: NPBE_GL
  integer :: NPBN_GL
  integer :: NPBS_GL

  ! *** Global periodic boundary conditions
  integer,allocatable,dimension (:) :: IPBS_GL
  integer,allocatable,dimension (:) :: IPBW_GL
  integer,allocatable,dimension (:) :: IPBE_GL
  integer,allocatable,dimension (:) :: IPBN_GL

  integer,allocatable,dimension (:) :: JPBS_GL
  integer,allocatable,dimension (:) :: JPBW_GL
  integer,allocatable,dimension (:) :: JPBE_GL
  integer,allocatable,dimension (:) :: JPBN_GL

  integer,allocatable,dimension (:) :: LPBS_GL
  integer,allocatable,dimension (:) :: LPBW_GL
  integer,allocatable,dimension (:) :: LPBE_GL
  integer,allocatable,dimension (:) :: LPBN_GL

  integer,allocatable,dimension (:) :: ISPBS_GL
  integer,allocatable,dimension (:) :: ISPBW_GL
  integer,allocatable,dimension (:) :: ISPBE_GL
  integer,allocatable,dimension (:) :: ISPBN_GL

  integer,allocatable,dimension (:) :: ISPRS_GL
  integer,allocatable,dimension (:) :: ISPRW_GL
  integer,allocatable,dimension (:) :: ISPRE_GL
  integer,allocatable,dimension (:) :: ISPRN_GL

  integer,allocatable,dimension (:) :: NPSERS_GL
  integer,allocatable,dimension (:) :: NPSERW_GL
  integer,allocatable,dimension (:) :: NPSERE_GL
  integer,allocatable,dimension (:) :: NPSERN_GL

  integer,allocatable,dimension (:) :: NPSERS1_GL
  integer,allocatable,dimension (:) :: NPSERW1_GL
  integer,allocatable,dimension (:) :: NPSERE1_GL
  integer,allocatable,dimension (:) :: NPSERN1_GL
  
  ! *** Local values are mapped in Map_OpenBC_Pressure
  real,allocatable,dimension(:,:)   :: PCBE_GL
  real,allocatable,dimension(:,:)   :: PCBN_GL
  real,allocatable,dimension(:,:)   :: PCBS_GL
  real,allocatable,dimension(:,:)   :: PCBW_GL
                                    
  real,allocatable,dimension(:,:)   :: PSBE_GL
  real,allocatable,dimension(:,:)   :: PSBN_GL
  real,allocatable,dimension(:,:)   :: PSBS_GL
  real,allocatable,dimension(:,:)   :: PSBW_GL

  real,allocatable,dimension (:)    :: TPCOORDS_GL
  real,allocatable,dimension (:)    :: TPCOORDW_GL
  real,allocatable,dimension (:)    :: TPCOORDE_GL
  real,allocatable,dimension (:)    :: TPCOORDN_GL
  ! *** End Map_OpenBC_Pressure

  ! *** Concentration global data arrays
  integer,allocatable,dimension (:)   :: ICBS_GL
  integer,allocatable,dimension (:)   :: JCBS_GL
  integer,allocatable,dimension (:)   :: NTSCRS_GL
  integer,allocatable,dimension (:,:) :: NCSERS_GL

  integer,allocatable,dimension (:)   :: ICBE_GL
  integer,allocatable,dimension (:)   :: JCBE_GL
  integer,allocatable,dimension (:)   :: NTSCRE_GL
  integer,allocatable,dimension (:,:) :: NCSERE_GL
  
  integer,allocatable,dimension (:)   :: ICBW_GL
  integer,allocatable,dimension (:)   :: JCBW_GL
  integer,allocatable,dimension (:)   :: NTSCRW_GL
  integer,allocatable,dimension (:,:) :: NCSERW_GL

  integer,allocatable,dimension (:)   :: ICBN_GL
  integer,allocatable,dimension (:)   :: JCBN_GL
  integer,allocatable,dimension (:)   :: NTSCRN_GL
  integer,allocatable,dimension (:,:) :: NCSERN_GL

  real,allocatable,dimension (:,:,:)  :: CBS_GL
  real,allocatable,dimension (:,:,:)  :: CBE_GL
  real,allocatable,dimension (:,:,:)  :: CBW_GL
  real,allocatable,dimension (:,:,:)  :: CBN_GL

  ! *** River global data arrays
  !integer,allocatable,dimension (:)   :: IQS_GL
  !integer,allocatable,dimension (:)   :: JQS_GL  delme
  !
  !integer,allocatable,dimension (:)   :: NQSMUL_GL
  !integer,allocatable,dimension (:)   :: NQSMF_GL
  !integer,allocatable,dimension (:)   :: NQSERQ_GL
  !integer,allocatable,dimension (:,:) :: NCSERQ_GL
  !integer,allocatable,dimension (:)   :: GRPID_GL
  !integer,allocatable,dimension (:)   :: LQS_GL

  !real,allocatable,dimension (:)      :: QSSE_GL
  !real,allocatable,dimension (:)      :: QFACTOR_GL
  !real,allocatable,dimension (:)      :: QWIDTH_GL
  !real,allocatable,dimension (:,:)    :: CQSE_GL
  ! *** End river global arrays

  integer :: MLTMSR_GL
  integer,allocatable,dimension (:) :: ILTMSR_GL
  integer,allocatable,dimension (:) :: JLTMSR_GL
  integer,allocatable,dimension (:) :: NTSSSS_GL
  integer,allocatable,dimension (:) :: MTMSRP_GL
  integer,allocatable,dimension (:) :: MTMSRC_GL
  integer,allocatable,dimension (:) :: MTMSRA_GL
  integer,allocatable,dimension (:) :: MTMSRUE_GL
  integer,allocatable,dimension (:) :: MTMSRUT_GL
  integer,allocatable,dimension (:) :: MTMSRU_GL
  integer,allocatable,dimension (:) :: MTMSRQE_GL
  integer,allocatable,dimension (:) :: MTMSRQ_GL
  integer,allocatable,dimension (:) :: MLTM_GL

  character(len = 20),allocatable,dimension (:):: CLTMSR_GL

End module Variables_MPI_Mapping
