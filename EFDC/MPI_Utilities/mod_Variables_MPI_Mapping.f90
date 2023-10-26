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
! @details Global arrays for mapping out boundary conditions and other 
!  mapping values
! @author Zander Mausolff
! @date 11/01/2019
!---------------------------------------------------------------------------!
Module Variables_MPI_Mapping

  Implicit None
  
  Integer, Allocatable, Dimension(:) :: IL_Global              ! *** Global I index based on Local L index
  Integer, Allocatable, Dimension(:) :: JL_Global              ! *** Global J index based on Local L index
  Integer, Allocatable, Dimension(:) :: IL_GL                  ! *** Global I index based on Global L index
  Integer, Allocatable, Dimension(:) :: JL_GL                  ! *** Global J index based on Global L index
  
  Integer :: La_local_no_ghost

  Integer :: num_active_l_local
  
  INTEGER,SAVE :: NNoGhost !< Stores the number of ghost cells for a given subdomain
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: NoGhost !< Array containing active cells, excluding ghost cells
  LOGICAL,SAVE,ALLOCATABLE,DIMENSION(:) :: IsGhost
  Real,SAVE,ALLOCATABLE,DIMENSION(:)    :: GhostMask

  Integer, Allocatable, Dimension(:) :: lij_west_conn_outside !< Contains mapping of i/j outside of domain to LA values for connection arrays
  Integer, Allocatable, Dimension(:) :: lij_east_conn_outside !< Contains mapping of i/j outside of domain to LA values for connection arrays

  Integer :: offset_connectors_for_ns

  ! *** Create type that helps in the global remapping
  Type mapping_ij
    Integer :: process   !< Process that the cell belongs to 
    Integer :: local_i   !< Local I index for that subdomain
    Integer :: local_j   !< Local J index for that subdomain
    Integer :: global_i  !< Corresponding global I index for the given local I subdomain value
    Integer :: global_j  !< Corresponding global J index for the given local I subdomain value
    Integer :: local_l   !< Local L index for given subdomain
    Integer :: global_l  !< Corresponding global L value for the given local L subdomain value. This is the global L that EE will understand
  End Type

  Type(mapping_ij), Allocatable, Dimension(:) :: loc_to_glob !< Local type that contains local i/j and corresponding global i/j
  Type(mapping_ij), Allocatable, Dimension(:) :: all_loc_to_glob !< Global type that contains ALL (even inactive cells) local i/j and corresponding global i/j
  Type(mapping_ij), Allocatable, Dimension(:) :: active_loc_to_glob !< Maps Global --> ALL local values for only the active cells (globally sized)
  Type(mapping_ij), Allocatable, Dimension(:) :: sorted_loc_to_glob !< Sorted array containing local i/j corresponding to global i/j.  These values are indexed off of the global L value

  Real, Allocatable, Dimension(:,:,:) :: GWCSER_Global !< GROUND WATER CONCENTRATION SERIES DATA
  Integer, Allocatable, Dimension(:)  :: NGWSL_Global

  ! *** Global periodic boundary conditions
  Integer :: NPBW_GL
  Integer :: NPBE_GL
  Integer :: NPBN_GL
  Integer :: NPBS_GL

  ! *** Global periodic boundary conditions
  INTEGER,ALLOCATABLE,DIMENSION (:) :: IPBS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: IPBW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: IPBE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: IPBN_GL

  INTEGER,ALLOCATABLE,DIMENSION (:) :: JPBS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: JPBW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: JPBE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: JPBN_GL

  INTEGER,ALLOCATABLE,DIMENSION (:) :: LPBS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: LPBW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: LPBE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: LPBN_GL

  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPBS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPBW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPBE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPBN_GL

  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPRS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPRW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPRE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ISPRN_GL

  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERN_GL

  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERS1_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERW1_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERE1_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NPSERN1_GL
  
  INTEGER,ALLOCATABLE,DIMENSION (:) :: TPCOORDS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: TPCOORDW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: TPCOORDE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: TPCOORDN_GL
 
  ! *** Local values are mapped in Map_OpenBC_Pressure
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PCBE_GL
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PCBN_GL
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PCBS_GL
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PCBW_GL
                                    
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PSBE_GL
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PSBN_GL
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PSBS_GL
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PSBW_GL
  ! *** End Map_OpenBC_Pressure

  ! *** Concentration global data arrays
  Integer,ALLOCATABLE,DIMENSION (:)   :: ICBS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: JCBS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NTSCRS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:,:) :: NCSERS_GL

  INTEGER,ALLOCATABLE,DIMENSION (:)   :: ICBE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: JCBE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NTSCRE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:,:) :: NCSERE_GL
  
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: ICBW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: JCBW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NTSCRW_GL
  INTEGER,ALLOCATABLE,DIMENSION (:,:) :: NCSERW_GL

  INTEGER,ALLOCATABLE,DIMENSION (:)   :: ICBN_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: JCBN_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NTSCRN_GL
  INTEGER,ALLOCATABLE,DIMENSION (:,:) :: NCSERN_GL

  REAL,ALLOCATABLE,DIMENSION (:,:,:)  :: CBS_GL
  REAL,ALLOCATABLE,DIMENSION (:,:,:)  :: CBE_GL
  REAL,ALLOCATABLE,DIMENSION (:,:,:)  :: CBW_GL
  REAL,ALLOCATABLE,DIMENSION (:,:,:)  :: CBN_GL

  ! *** River global data arrays
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: IQS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: JQS_GL

  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NQSMUL_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NQSMF_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: NQSERQ_GL
  INTEGER,ALLOCATABLE,DIMENSION (:,:) :: NCSERQ_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: GRPID_GL
  INTEGER,ALLOCATABLE,DIMENSION (:)   :: LQS_GL

  REAL,ALLOCATABLE,DIMENSION (:)      :: QSSE_GL
  REAL,ALLOCATABLE,DIMENSION (:)      :: QFACTOR_GL
  REAL,ALLOCATABLE,DIMENSION (:)      :: QWIDTH_GL
  REAL,ALLOCATABLE,DIMENSION (:,:)    :: CQSE_GL
  ! *** End river global arrays

  INTEGER :: MLTMSR_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: ILTMSR_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: JLTMSR_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: NTSSSS_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRP_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRC_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRA_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRUE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRUT_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRU_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRQE_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MTMSRQ_GL
  INTEGER,ALLOCATABLE,DIMENSION (:) :: MLTM_GL

  CHARACTER(len=20),ALLOCATABLE,DIMENSION (:):: CLTMSR_GL

End module Variables_MPI_Mapping
