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
! Module: Variables_MPI_Write_Out
!
!> @details Contains global variables that are written to the EE_*.OUT files
!
!> @author Zander Mausolff
!---------------------------------------------------------------------------!
! CHANGE RECORD
! DATE MODIFIED     BY               DESCRIPTION
!---------------------------------------------------------------------------!
!    2019-09       Zander    Added global variables that will be written out
!---------------------------------------------------------------------------!

Module Variables_MPI_Write_Out

  use GLOBAL

  implicit none

  Real, Allocatable, Dimension(:,:,:) :: RICEWHT_global

  real, target, allocatable, dimension(:) :: CUE_Global
  real, target, allocatable, dimension(:) :: CUN_Global
  real, target, allocatable, dimension(:) :: CVN_Global
  real, target, allocatable, dimension(:) :: CVE_Global
  
  real, target, allocatable, dimension(:,:)  :: DZC_Global
  real, target, allocatable, dimension(:)    :: TAUBSED_Global
  real, target, allocatable, dimension(:)    :: TAUBSND_Global
  real, target, allocatable, dimension(:)    :: TAUB_Global
  real, target, allocatable, dimension(:)    :: WNDVELE_Global
  real, target, allocatable, dimension(:)    :: WNDVELN_GLobal
  real, target, allocatable, dimension(:)    :: PATMT_Global
  
  Integer(IK4) :: JMN_Global, JMX_Global !> J min/max for the global grid
  Integer(IK4) :: IMN_Global, IMX_Global !> I min/max for the global grid
  
  integer, target, allocatable,dimension(:) :: LBCS_Global
  
  real(RK4),allocatable :: Values_Temp_Binary(:,:,:,:)
  
  Character*9 :: MPI_Outdir !< Directory within output that contains optional output logs for MPI run

  integer, Allocatable, Dimension(:) :: LWVCELL_Global
  logical, allocatable,dimension(:) :: LWVMASK_Global          !< MASK TO DETERMINE IF WAVE CALCUATIONS ARE ON/OFF FOR EACH CELL

  Real,Target, Allocatable, Dimension(:) :: TBY_Global
  Real,Target, Allocatable, Dimension(:) :: TBX_Global
  Real,Target, Allocatable, Dimension(:) :: TSY_Global
  Real,Target, Allocatable, Dimension(:) :: TSX_Global
  
  Real,Target, Allocatable, Dimension(:) :: TBY1_Global
  Real,Target, Allocatable, Dimension(:) :: TBX1_Global
  Real,Target, Allocatable, Dimension(:) :: TSY1_Global
  Real,Target, Allocatable, Dimension(:) :: TSX1_Global

  integer,Target, Allocatable, Dimension(:) :: LWC_Global
  integer,Target, Allocatable, Dimension(:) :: LEC_Global
  integer,Target, Allocatable, Dimension(:) :: LSC_Global
  integer,Target, Allocatable, Dimension(:) :: LNC_Global

  ! *** Hydrodynamics and other modules
  integer, Allocatable, Dimension(:)  :: KSZ_Global

  integer,Target, Allocatable, Dimension(:)  :: ISCDRY_Global
  integer,Target, Allocatable, Dimension(:)  :: NATDRY_Global
  integer,Target, Allocatable, Dimension(:)  :: IDRY_Global
  integer,Target, Allocatable, Dimension(:)  :: MVEG_Global

  Real,Target, Allocatable, Dimension(:)     :: BELV_Global
  Real,Target, Allocatable, Dimension(:)     :: DXP_Global
  Real,Target, Allocatable, Dimension(:)     :: DYP_Global
  Real,Target, Allocatable, Dimension(:)     :: HP_Global
  Real,Target, Allocatable, Dimension(:)     :: H1P_Global
  Real,Target, Allocatable, Dimension(:)     :: H2P_Global
  Real,Target, Allocatable, Dimension(:)     :: HWQ_Global
  Real,Target, Allocatable, Dimension(:)     :: H2WQ_Global
  Real,Target, Allocatable, Dimension(:)     :: ZBR_Global

  Real,Target, Allocatable, Dimension(:)     :: SHEAR_Global
  Real,Target, Allocatable, Dimension(:)     :: SHEAR_Global2
  Real,Target, Allocatable, Dimension(:)     :: SHEAR_Local
  Real,Target, Allocatable, Dimension(:)     :: SHEAR_Local2

  Real,Target, Allocatable, Dimension(:)     :: EVAPSW_Global
  Real,Target, Allocatable, Dimension(:)     :: EVAPGW_Global
  Real,Target, Allocatable, Dimension(:)     :: QGW_Global
  Real,Target, Allocatable, Dimension(:)     :: AGWELV_Global
                                            
  Real,Target, Allocatable, Dimension(:)     :: EVAPT_Global
  Real,Target, Allocatable, Dimension(:)     :: RAINT_Global
  
  Real,Target, Allocatable, Dimension(:)     :: RSSBCE_Global
  Real,Target, Allocatable, Dimension(:)     :: RSSBCW_Global
  Real,Target, Allocatable, Dimension(:)     :: RSSBCN_Global
  Real,Target, Allocatable, Dimension(:)     :: RSSBCS_Global

  Real,Target, Allocatable, Dimension(:)     :: UHDYE_Global
  Real,Target, Allocatable, Dimension(:)     :: UHDY1E_Global
  Real,Target, Allocatable, Dimension(:)     :: VHDXE_Global
  Real,Target, Allocatable, Dimension(:)     :: VHDX1E_Global
  Real,Target, Allocatable, Dimension(:)     :: SUB_Global
  Real,Target, Allocatable, Dimension(:)     :: SVB_Global

  Real,Target, Allocatable, Dimension(:,:)   :: U_Global
  Real,Target, Allocatable, Dimension(:,:)   :: V_Global
  Real,Target, Allocatable, Dimension(:,:)   :: U1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: V1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: W_Global

  Real,Target, Allocatable, Dimension(:,:)   :: QQ_Global
  Real,Target, Allocatable, Dimension(:,:)   :: QQ1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: QQL_Global
  Real,Target, Allocatable, Dimension(:,:)   :: QQL1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: DML_Global

  Real,Target, Allocatable, Dimension(:,:)   :: TKE3D_Global
  Real,Target, Allocatable, Dimension(:,:)   :: EPS3D_Global
  Real,Target, Allocatable, Dimension(:,:)   :: GL3D_Global

  Real,Target, Allocatable, Dimension(:)     :: QSUME_Global
  Real,Target, Allocatable, Dimension(:,:)   :: QSUM_Global

  Real,Target, Allocatable, Dimension(:,:)   :: VHDX2_Global
  Real,Target, Allocatable, Dimension(:,:)   :: UHDY2_Global

  Real,Target, Allocatable, Dimension(:)     :: TEMB_Global
  Real,Target, Allocatable, Dimension(:)     :: SHAD_Global
  Real,Target, Allocatable, Dimension(:)     :: ICETHICK_Global
  Real,Target, Allocatable, Dimension(:)     :: ICETEMP_Global
  
  Real,Target, Allocatable, Dimension(:,:)   :: SAL_Global
  Real,Target, Allocatable, Dimension(:,:)   :: TEM_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: DYE_Global
  Real,Target, Allocatable, Dimension(:,:)   :: SFL_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: TOX_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SED_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SND_Global

  Real,Target, Allocatable, Dimension(:,:)   :: SAL1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: TEM1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: DYE1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: SFL1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: TOX1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SED1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SND1_Global

  integer,Target, Allocatable, Dimension(:)  :: KBT_Global
  integer, Allocatable, Dimension(:)         :: BEDMAP_Global
  Real,Target, Allocatable, Dimension(:,:)   :: BDENBED_Global
  Real,Target, Allocatable, Dimension(:,:)   :: PORBED_Global

  Real,Target, Allocatable, Dimension(:,:)   :: HBED_Global
  Real,Target, Allocatable, Dimension(:,:)   :: VDRBED_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SEDB_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SNDB_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: TOXB_Global

  Real,Target, Allocatable, Dimension(:,:)   :: HBED1_Global
  Real,Target, Allocatable, Dimension(:,:)   :: VDRBED1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SEDB1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: SNDB1_Global
  Real,Target, Allocatable, Dimension(:,:,:) :: TOXB1_Global

  Real,Target, Allocatable, Dimension(:,:)   :: QSBDLDX_Global
  Real,Target, Allocatable, Dimension(:,:)   :: QSBDLDY_Global

  ! *** SEDZLJ
  integer,  Target, Allocatable, Dimension(:,:)   :: LAYERACTIVE_Global
  real(rkd),Target, Allocatable, Dimension(:)     :: TAU_Global
  real(rkd),Target, Allocatable, Dimension(:)     :: D50AVG_Global
  real(rkd),Target, Allocatable, Dimension(:,:)   :: BULKDENS_Global
  real(rkd),Target, Allocatable, Dimension(:,:)   :: TSED_Global
  real(rkd),Target, Allocatable, Dimension(:,:,:) :: PERSED_Global
  real(rkd),Target, Allocatable, Dimension(:,:)   :: ERO_SED_FLX_Global
  real(rkd),Target, Allocatable, Dimension(:,:)   :: DEP_SED_FLX_Global
  real(rkd),Target, Allocatable, Dimension(:,:)   :: CBL_Global
  real(rkd),Target, Allocatable, Dimension(:,:)   :: CBLTOX_Global
  
  ! *** Waves
  type(Wave), Target, Allocatable, Dimension(:) :: WV_Global
  Real,Target, Allocatable, Dimension(:) :: WV_HEIGHT_Global
  Real,Target, Allocatable, Dimension(:) :: WV_PERIOD_Global
  Real,Target, Allocatable, Dimension(:) :: WV_DIR_Global
  Real,Target, Allocatable, Dimension(:) :: WV_DISSIPA_Global
  
  Real,Target, Allocatable, Dimension(:,:) :: FXWAVE_Global
  Real,Target, Allocatable, Dimension(:,:) :: FYWAVE_Global
  Real,Target, Allocatable, Dimension(:)   :: QQWV3_Global
  Real,Target, Allocatable, Dimension(:,:) :: WVHUU_Global
  Real,Target, Allocatable, Dimension(:,:) :: WVHVV_Global
  Real,Target, Allocatable, Dimension(:,:) :: WVHUV_Global
  
  ! *** Water Quality Variables - Rooted Plant and Epiphyte (RPEM)
  Real,Target, Allocatable, Dimension(:)   :: WQRPS_Global
  Real,Target, Allocatable, Dimension(:)   :: WQRPR_Global
  Real,Target, Allocatable, Dimension(:)   :: WQRPE_Global
  Real,Target, Allocatable, Dimension(:)   :: WQRPD_Global
  logical,Allocatable,Dimension(:)         :: LMASKRPEM_Global
  
  ! *** Water Quality Variables - Sediment Diagenesis
  Real,Target, Allocatable, Dimension(:,:) :: SMPOP_Global
  Real,Target, Allocatable, Dimension(:,:) :: SMPON_Global
  Real,Target, Allocatable, Dimension(:,:) :: SMPOC_Global
  Real,Target, Allocatable, Dimension(:,:) :: SMDFN_Global
  Real,Target, Allocatable, Dimension(:,:) :: SMDFP_Global
  Real,Target, Allocatable, Dimension(:,:) :: SMDFC_Global
  
  Real,Target, Allocatable, Dimension(:)   :: SM1NH4_Global
  Real,Target, Allocatable, Dimension(:)   :: SM2NH4_Global
  Real,Target, Allocatable, Dimension(:)   :: SM1NO3_Global
  Real,Target, Allocatable, Dimension(:)   :: SM2NO3_Global
  Real,Target, Allocatable, Dimension(:)   :: SM1PO4_Global
  Real,Target, Allocatable, Dimension(:)   :: SM2PO4_Global
  Real,Target, Allocatable, Dimension(:)   :: SM1H2S_Global
  Real,Target, Allocatable, Dimension(:)   :: SM2H2S_Global
  Real,Target, Allocatable, Dimension(:)   :: SM1SI_Global
  Real,Target, Allocatable, Dimension(:)   :: SM2SI_Global
  Real,Target, Allocatable, Dimension(:)   :: SMPSI_Global
  Real,Target, Allocatable, Dimension(:)   :: SMBST_Global
  Real,Target, Allocatable, Dimension(:)   :: SMT_Global
  Real,Target, Allocatable, Dimension(:)   :: SMCSOD_Global
  Real,Target, Allocatable, Dimension(:)   :: SMNSOD_Global
  Real,Target, Allocatable, Dimension(:)   :: WQBFNH4_Global
  Real,Target, Allocatable, Dimension(:)   :: WQBFNO3_Global
  Real,Target, Allocatable, Dimension(:)   :: WQBFO2_Global
  Real,Target, Allocatable, Dimension(:)   :: WQBFCOD_Global
  Real,Target, Allocatable, Dimension(:)   :: WQBFPO4D_Global
  Real,Target, Allocatable, Dimension(:)   :: WQBFSAD_Global
  
  ! *** Water Quality Variables - Water Column
  Real,Target, Allocatable, Dimension(:,:,:) :: WQV_Global

  ! *** Boundard Conditions
  integer,Target, Allocatable, Dimension(:,:,:) :: NLOS_Global
  integer,Target, Allocatable, Dimension(:,:,:) :: NLOE_Global
  integer,Target, Allocatable, Dimension(:,:,:) :: NLOW_Global
  integer,Target, Allocatable, Dimension(:,:,:) :: NLON_Global

  Real,Target, Allocatable, Dimension(:,:,:)    :: CLOS_Global
  Real,Target, Allocatable, Dimension(:,:,:)    :: CLOE_Global
  Real,Target, Allocatable, Dimension(:,:,:)    :: CLOW_Global
  Real,Target, Allocatable, Dimension(:,:,:)    :: CLON_Global

  ! *** Temporary Arrays
  integer, Allocatable, Dimension(:)       :: I1D_Global           !<  Temporary variable for reading global IC fields
  integer, Allocatable, Dimension(:,:)     :: I2D_Global           !<  Temporary variable for reading global IC fields
  integer, Allocatable, Dimension(:,:,:)   :: I3D_Global           !<  Temporary variable for reading global IC fields
  real(RKD), Allocatable, Dimension(:)     :: R1D_Global           !<  Temporary variable for reading global IC fields
  real(RKD), Allocatable, Dimension(:,:)   :: R2D_Global           !<  Temporary variable for reading global IC fields
  real(RKD), Allocatable, Dimension(:,:,:) :: R3D_Global           !<  Temporary variable for reading global IC fields

  ! *** Propwash
  !REAL,target, allocatable, dimension(:,:,:) :: SDF_Global         !< Cohesive sediment concentration without fast settling classes     (g/m3), Also used for all sediments classes when NSEDFLUME > 0
  
End module Variables_MPI_Write_Out
