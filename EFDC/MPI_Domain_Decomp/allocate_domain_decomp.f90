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
  !> @details Allocated arrays needed for the domain decomposition routines
  !---------------------------------------------------------------------------!
  !> @author Zander Mausolff
  !> @date 9/9/2019
  !---------------------------------------------------------------------------!
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !---------------------------------------------------------------------------!
  !    2022-01      PMC      Updated array allocation and initialization
  !    2019-09      Zander   Allocates variables required for domain decomp
  !
  !---------------------------------------------------------------------------!
  Subroutine Allocate_Domain_Decomp

  use GLOBAL
  use Allocate_Initialize
  use Variables_WQ
  use WQ_RPEM_MODULE
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out

  implicit none

  ! *** Local variables
  integer :: L, NMAX

  ! *** Special Range Cases
  allocate(IL2IG(-2:global_max_width_x+2))
  allocate(JL2JG(-2:global_max_width_y+2))
  allocate(IG2IL(-2:global_max_width_x+2))
  allocate(JG2JL(-2:global_max_width_y+2))
  IL2IG = 0
  JL2JG = 0
  IG2IL = 0
  JG2JL = 0

  ! *** Character arrays
  allocate(CLTMSR_GL(MLTMSRM))
  !allocate(GRPID_GL(NQSIJM))

  ! *** Allocate general MPI arrays and initialize
  call AllocateDSI( lij_west_conn_outside, NPEWBP,  0)
  call AllocateDSI( lij_east_conn_outside, NPEWBP,  0)

  call AllocateDSI( GWCSER_Global,   NDGWSER, NGWSERM, NSTVM2, 0.0)
  call AllocateDSI( NGWSL_Global,    LCM_Global,    0)
                                                  
  call AllocateDSI( ISCDRY_Global,   LCM_Global,    0)
  call AllocateDSI( NATDRY_Global,   LCM_Global,    0)
  call AllocateDSI( IDRY_Global,     LCM_Global,    0)
                                                  
  call AllocateDSI( BELV_Global,     LCM_Global,  0.0)
  call AllocateDSI( DXP_Global,      LCM_Global,  0.0)
  call AllocateDSI( DYP_Global,      LCM_Global,  0.0)
  call AllocateDSI( HP_Global,       LCM_Global,  0.0)
  call AllocateDSI( H1P_Global,      LCM_Global,  0.0)
  call AllocateDSI( H2P_Global,      LCM_Global,  0.0)
  call AllocateDSI( HWQ_Global,      LCM_Global,  0.0)
  call AllocateDSI( H2WQ_Global,     LCM_Global,  0.0)
  call AllocateDSI( ZBR_Global,      LCM_Global,  0.0)
  call AllocateDSI( MVEG_Global,     LCM_Global,    0)
                                                  
  call AllocateDSI( SHEAR_Global,    LCM_Global,  0.0)
  call AllocateDSI( SHEAR_Global2,   LCM_Global,  0.0)
  call AllocateDSI( SHEAR_Local,     LCM,         0.0)    
  call AllocateDSI( SHEAR_Local2,    LCM,         0.0)    
                                                  
  call AllocateDSI( EVAPSW_Global,   LCM_Global,  0.0)
  call AllocateDSI( EVAPGW_Global,   LCM_Global,  0.0)
  call AllocateDSI( QGW_Global,      LCM_Global,  0.0)
  call AllocateDSI( AGWELV_Global,   LCM_Global,  0.0)
                                                  
  call AllocateDSI( EVAPT_Global,    LCM_Global,  0.0)
  call AllocateDSI( RAINT_Global,    LCM_Global,  0.0)
                                                  
  call AllocateDSI( RSSBCE_Global,   LCM_Global,  0.0)
  call AllocateDSI( RSSBCW_Global,   LCM_Global,  0.0)
  call AllocateDSI( RSSBCN_Global,   LCM_Global,  0.0)
  call AllocateDSI( RSSBCS_Global,   LCM_Global,  0.0)
                                                  
  call AllocateDSI( UHDYE_Global,    LCM_Global,  0.0)
  call AllocateDSI( UHDY1E_Global,   LCM_Global,  0.0)
  call AllocateDSI( VHDXE_Global,    LCM_Global,  0.0)
  call AllocateDSI( VHDX1E_Global,   LCM_Global,  0.0)
  call AllocateDSI( SUB_Global,      LCM_Global,  0.0)
  call AllocateDSI( SVB_Global,      LCM_Global,  0.0)

  call AllocateDSI( U_Global,        LCM_Global,  KCM, 0.0)
  call AllocateDSI( V_Global,        LCM_Global,  KCM, 0.0)
  call AllocateDSI( U1_Global,       LCM_Global,  KCM, 0.0)
  call AllocateDSI( V1_Global,       LCM_Global,  KCM, 0.0)
  call AllocateDSI( W_Global,        LCM_Global,  KCM, 0.0)

  call AllocateDSI( QQ_Global,       LCM_Global, -KCM, 0.0)
  call AllocateDSI( QQ1_Global,      LCM_Global, -KCM, 0.0)
  call AllocateDSI( QQL_Global,      LCM_Global, -KCM, 0.0)
  call AllocateDSI( QQL1_Global,     LCM_Global, -KCM, 0.0)
  call AllocateDSI( DML_Global,      LCM_Global, -KCM, 0.0)

  call AllocateDSI( QSUM_Global,     LCM_Global,  KCM, 0.0)
  call AllocateDSI( QSUME_Global,    LCM_Global,  0.0)

  call AllocateDSI( VHDX2_Global,    LCM_Global,  KCM, 0.0)
  call AllocateDSI( UHDY2_Global,    LCM_Global,  KCM, 0.0)

  call AllocateDSI( TBX_Global,      LCM_Global,  0.0)
  call AllocateDSI( TBY_Global,      LCM_Global,  0.0)
  call AllocateDSI( TSX_Global,      LCM_Global,  0.0)
  call AllocateDSI( TSY_Global,      LCM_Global,  0.0)

  call AllocateDSI( TBX1_Global,     LCM_Global,  0.0)
  call AllocateDSI( TBY1_Global,     LCM_Global,  0.0)
  call AllocateDSI( TSX1_Global,     LCM_Global,  0.0)
  call AllocateDSI( TSY1_Global,     LCM_Global,  0.0)

  if( ISGOTM > 0 )then
    call AllocateDSI( TKE3D_Global,  LCM_Global, -KCM, 0.0)
    call AllocateDSI( EPS3D_Global,  LCM_Global, -KCM, 0.0)
    call AllocateDSI( GL3D_Global,   LCM_Global, -KCM, 0.0)
  endif
  
  if( ISTRAN(1) > 0 )then
    call AllocateDSI( SAL_Global,    LCM_Global, KCM, 0.0)
    call AllocateDSI( SAL1_Global,   LCM_Global, KCM, 0.0)
  endif

  if( ISTRAN(2) > 0 )then
    call AllocateDSI( TEM_Global,    LCM_Global, KCM, 0.0)
    call AllocateDSI( TEM1_Global,   LCM_Global, KCM, 0.0)
    call AllocateDSI( TEMB_Global,   LCM_Global, 0.0)
    call AllocateDSI( SHAD_Global,   LCM_Global, 0.0)
    if( ISICE > 0 )then
      call AllocateDSI( ICETHICK_Global, LCM_Global, 0.0)
      call AllocateDSI( ICETEMP_Global,  LCM_Global, 0.0)
    endif
  endif
  if( ISTRAN(3) > 0 )then
    call AllocateDSI( DYE_Global,    LCM_Global, KCM, NDYM, 0.0)
    call AllocateDSI( DYE1_Global,   LCM_Global, KCM, NDYM, 0.0)
  endif
  if( ISTRAN(4) > 0 )then
    call AllocateDSI( SFL_Global,    LCM_Global, KCM, 0.0)
  endif

  if( ISTRAN(5) > 0 )then
    call AllocateDSI( TOX_Global,    LCM_Global, KCM, NTXM, 0.0)
    call AllocateDSI( TOX1_Global,   LCM_Global, KCM, NTXM, 0.0)
    call AllocateDSI( TOXB_Global,   LCM_Global, KBM, NTXM, 0.0)
    call AllocateDSI( TOXB1_Global,  LCM_Global, KBM, NTXM, 0.0)
  endif

  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    call AllocateDSI( HBED_Global,    LCM_Global, KBM, 0.0)
    call AllocateDSI( HBED1_Global,   LCM_Global, KBM, 0.0)

    call AllocateDSI( BEDMAP_Global,  LCM_Global, 0)
    call AllocateDSI( KBT_Global,     LCM_Global, 0)

    call AllocateDSI( BDENBED_Global, LCM_Global, KBM, 0.0)
    call AllocateDSI( PORBED_Global,  LCM_Global, KBM, 0.0)
    call AllocateDSI( VDRBED_Global,  LCM_Global, KBM, 0.0)
    call AllocateDSI( VDRBED1_Global, LCM_Global, KBM, 0.0)

    if( ISTRAN(7) > 0 .and. ICALC_BL > 0 .and. NSND > 0 )then
      call AllocateDSI( QSBDLDX_Global,  LCM_Global, NSND, 0.0)
      call AllocateDSI( QSBDLDY_Global,  LCM_Global, NSND, 0.0)
    elseif( ICALC_BL > 0 )then
      call AllocateDSI( QSBDLDX_Global,  LCM_Global, NSEDS, 0.0)
      call AllocateDSI( QSBDLDY_Global,  LCM_Global, NSEDS, 0.0)
    endif

    if( NSEDFLUME > 0 )then
      call AllocateDSI( LAYERACTIVE_Global, KB,         LCM_Global,   0)
      call AllocateDSI( ERO_SED_FLX_Global, LCM_Global, NSEDS2,     0.0)
      call AllocateDSI( DEP_SED_FLX_Global, LCM_Global, NSEDS2,     0.0)
      call AllocateDSI( TAU_Global,         LCM_Global, 0.0)
      call AllocateDSI( D50AVG_Global,      LCM_Global, 0.0)
      call AllocateDSI( BULKDENS_Global,    KB,         LCM_Global, 0.0)
      call AllocateDSI( TSED_Global,        KB,         LCM_Global, 0.0)
      call AllocateDSI( TSED0_Global,       KB,         LCM_Global, 0.0)
      call AllocateDSI( PERSED_Global,      NSEDS,      KB,         LCM_Global, 0.0)
      if( ICALC_BL > 0 )then                
        call AllocateDSI( CBL_Global,       LCM_Global, NSEDS, 0.0)
        if( ISTRAN(5) > 0 )then
          call AllocateDSI( CBLTOX_Global, LCM_Global, NTXM, 0.0)
        endif
      endif
    endif
  endif

  if( ISTRAN(6) > 0 )then
    call AllocateDSI( SED_Global,   LCM_Global, KCM, NSEDS2, 0.0)
    call AllocateDSI( SED1_Global,  LCM_Global, KCM, NSEDS2, 0.0)
    call AllocateDSI( SEDB_Global,  LCM_Global, KBM, NSEDS,  0.0)
    call AllocateDSI( SEDB1_Global, LCM_Global, KBM, NSEDS,  0.0)
    ! *** Propwash fast settling
    !Call AllocateDSI( SDF_Global,   LCM_Global, KCM, NSCM, 0.0)
  endif
  if( ISTRAN(7) > 0 )then
    call AllocateDSI( SND_Global,   LCM_Global, KCM, NSNM, 0.0)
    call AllocateDSI( SND1_Global,  LCM_Global, KCM, NSNM, 0.0)
    call AllocateDSI( SNDB_Global,  LCM_Global, KBM, NSNM, 0.0)
    call AllocateDSI( SNDB1_Global, LCM_Global, KBM, NSNM, 0.0)
  endif

  ! *** Water Quality
  if( ISTRAN(8) > 0 )then
    call AllocateDSI( WQV_Global,         LCM_Global, KCM, -NWQVM, 0.0)

    if( ISRPEM > 0 )then
      call AllocateDSI( WQRPS_Global,     LCM_Global, 0.0)
      call AllocateDSI( WQRPR_Global,     LCM_Global, 0.0)
      call AllocateDSI( WQRPE_Global,     LCM_Global, 0.0)
      call AllocateDSI( WQRPD_Global,     LCM_Global, 0.0)
      call AllocateDSI( LMASKRPEM_Global, LCM_Global, .false.)
    endif

    ! *** Sediment Diagenesis
    call AllocateDSI( SMPON_Global,    LCM_Global, NSMGM, 0.0 )
    call AllocateDSI( SMPOP_Global,    LCM_Global, NSMGM, 0.0 )
    call AllocateDSI( SMPOC_Global,    LCM_Global, NSMGM, 0.0 )
    call AllocateDSI( SMDFN_Global,    LCM_Global, NSMGM, 0.0 )
    call AllocateDSI( SMDFP_Global,    LCM_Global, NSMGM, 0.0 )
    call AllocateDSI( SMDFC_Global,    LCM_Global, NSMGM, 0.0 )

    call AllocateDSI( SM1NH4_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM2NH4_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM1NO3_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM2NO3_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM1PO4_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM2PO4_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM1H2S_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM2H2S_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SM1SI_Global,    LCM_Global, 0.0 )
    call AllocateDSI( SM2SI_Global,    LCM_Global, 0.0 )
    call AllocateDSI( SMPSI_Global,    LCM_Global, 0.0 )
    call AllocateDSI( SMBST_Global,    LCM_Global, 0.0 )
    call AllocateDSI( SMT_Global,      LCM_Global, 0.0 )
    call AllocateDSI( SMCSOD_Global,   LCM_Global, 0.0 )
    call AllocateDSI( SMNSOD_Global,   LCM_Global, 0.0 )
    call AllocateDSI( WQBFNH4_Global,  LCM_Global, 0.0 )
    call AllocateDSI( WQBFNO3_Global,  LCM_Global, 0.0 )
    call AllocateDSI( WQBFO2_Global,   LCM_Global, 0.0 )
    call AllocateDSI( WQBFCOD_Global,  LCM_Global, 0.0 )
    call AllocateDSI( WQBFPO4D_Global, LCM_Global, 0.0 )
    call AllocateDSI( WQBFSAD_Global,  LCM_Global, 0.0 )
  endif

  call AllocateDSI( WVHUU_Global,  LCM_Global,  KCM, 0.0)
  call AllocateDSI( WVHVV_Global,  LCM_Global,  KCM, 0.0)
  call AllocateDSI( WVHUV_Global,  LCM_Global,  KCM, 0.0)
  call AllocateDSI( QQWV3_Global,  LCM_Global,  0.0)
                                   
  call AllocateDSI( KSZ_Global,    LCM_Global,    1)
                                   
  call AllocateDSI( LWC_Global,    LCM_Global,    0)
  call AllocateDSI( LEC_Global,    LCM_Global,    0)
  call AllocateDSI( LSC_Global,    LCM_Global,    0)
  call AllocateDSI( LNC_Global,    LCM_Global,    0)

  call AllocateDSI( UMASK_Global, -LCM_Global,    0)
  call AllocateDSI( VMASK_Global, -LCM_Global,    0)

  if( ISWAVE > 0 )then

    call AllocateDSI( FXWAVE_Global,    LCM_Global, KCM, 0.0)
    call AllocateDSI( FYWAVE_Global,    LCM_Global, KCM, 0.0)
    call AllocateDSI( WV_HEIGHT_Global, LCM_Global, 0.0)
    call AllocateDSI( WV_PERIOD_Global, LCM_Global, 0.0)
    call AllocateDSI( WV_DIR_Global,    LCM_Global, 0.0)
    call AllocateDSI( WV_DISSIPA_Global,LCM_Global, 0.0)

    allocate(WV_Global(LCM_Global))
    do L = 1, LCM_Global
      call AllocateDSI( WV_Global(L).DISSIPA, KCM, 0.0)
    enddo
  endif

  call AllocateDSI( IJCT_GLOBAL,   global_max_width_x, global_max_width_y, 0)
  call AllocateDSI( IJCTLT_GLOBAL, global_max_width_x, global_max_width_y, 0)
  call AllocateDSI( LIJ_Global,    global_max_width_x, global_max_width_y, 0)

  ! *** Open Boundary Conditions
  call AllocateDSI( IPBS_GL,         NPBSM,   0)
  call AllocateDSI( IPBW_GL,         NPBWM,   0)
  call AllocateDSI( IPBE_GL,         NPBEM,   0)
  call AllocateDSI( IPBN_GL,         NPBNM,   0)

  call AllocateDSI( JPBS_GL,         NPBSM,   0)
  call AllocateDSI( JPBW_GL,         NPBWM,   0)
  call AllocateDSI( JPBE_GL,         NPBEM,   0)
  call AllocateDSI( JPBN_GL,         NPBNM,   0)

  call AllocateDSI( LPBS_GL,         NPBSM,   0)
  call AllocateDSI( LPBW_GL,         NPBWM,   0)
  call AllocateDSI( LPBE_GL,         NPBEM,   0)
  call AllocateDSI( LPBN_GL,         NPBNM,   0)

  call AllocateDSI( ISPBS_GL,        NPBSM,   0)
  call AllocateDSI( ISPBW_GL,        NPBWM,   0)
  call AllocateDSI( ISPBE_GL,        NPBEM,   0)
  call AllocateDSI( ISPBN_GL,        NPBNM,   0)

  call AllocateDSI( ISPRS_GL,        NPBSM,   0)
  call AllocateDSI( ISPRW_GL,        NPBWM,   0)
  call AllocateDSI( ISPRE_GL,        NPBEM,   0)
  call AllocateDSI( ISPRN_GL,        NPBNM,   0)

  call AllocateDSI( NPSERS_GL,       NPBSM,   0)
  call AllocateDSI( NPSERW_GL,       NPBWM,   0)
  call AllocateDSI( NPSERE_GL,       NPBEM,   0)
  call AllocateDSI( NPSERN_GL,       NPBNM,   0)

  call AllocateDSI( NPSERS1_GL,      NPBSM,   0)
  call AllocateDSI( NPSERW1_GL,      NPBWM,   0)
  call AllocateDSI( NPSERE1_GL,      NPBEM,   0)
  call AllocateDSI( NPSERN1_GL,      NPBNM,   0)

  call AllocateDSI( PCBS_GL,         NPBSM, MTM, 0.0)
  call AllocateDSI( PCBW_GL,         NPBWM, MTM, 0.0)
  call AllocateDSI( PCBE_GL,         NPBEM, MTM, 0.0)
  call AllocateDSI( PCBN_GL,         NPBNM, MTM, 0.0)

  call AllocateDSI( PSBS_GL,         NPBSM, MTM, 0.0)
  call AllocateDSI( PSBW_GL,         NPBWM, MTM, 0.0)
  call AllocateDSI( PSBE_GL,         NPBEM, MTM, 0.0)
  call AllocateDSI( PSBN_GL,         NPBNM, MTM, 0.0)

  call AllocateDSI( TPCOORDS_GL,     NPBSM,   0.0)
  call AllocateDSI( TPCOORDW_GL,     NPBWM,   0.0)
  call AllocateDSI( TPCOORDE_GL,     NPBEM,   0.0)
  call AllocateDSI( TPCOORDN_GL,     NPBNM,   0.0)

  ! *** Concentration global data arrays
  call AllocateDSI( ICBS_GL,         NBBSM,   0)
  call AllocateDSI( JCBS_GL,         NBBSM,   0)
  call AllocateDSI( NTSCRS_GL,       NBBSM,   0)
  call AllocateDSI( NCSERS_GL,       NBBSM,    NSTVM2,   0)
  call AllocateDSI( CBS_GL,          NBBSM, 2, NSTVM2, 0.0)
  call AllocateDSI( ICBW_GL,         NBBWM,   0)
  call AllocateDSI( JCBW_GL,         NBBWM,   0)
  call AllocateDSI( NTSCRW_GL,       NBBWM,   0)
  call AllocateDSI( NCSERW_GL,       NBBWM,    NSTVM2,   0)
  call AllocateDSI( CBW_GL,          NBBWM, 2, NSTVM2, 0.0)
  call AllocateDSI( ICBE_GL,         NBBEM,   0)
  call AllocateDSI( JCBE_GL,         NBBEM,   0)
  call AllocateDSI( NTSCRE_GL,       NBBEM,   0)
  call AllocateDSI( NCSERE_GL,       NBBEM,    NSTVM2,   0)
  call AllocateDSI( CBE_GL,          NBBEM, 2, NSTVM2, 0.0)
  call AllocateDSI( ICBN_GL,         NBBNM,   0)
  call AllocateDSI( JCBN_GL,         NBBNM,   0)
  call AllocateDSI( NTSCRN_GL,       NBBNM,   0)
  call AllocateDSI( NCSERN_GL,       NBBNM,    NSTVM2,   0)
  call AllocateDSI( CBN_GL,          NBBNM, 2, NSTVM2, 0.0)
  call AllocateDSI( ILTMSR_GL,       MLTMSRM, 0)
  call AllocateDSI( JLTMSR_GL,       MLTMSRM, 0)
  call AllocateDSI( NTSSSS_GL,       MLTMSRM, 0)

  call AllocateDSI( MTMSRP_GL,       MLTMSRM, 0)
  call AllocateDSI( MTMSRC_GL,       MLTMSRM, 0)
  call AllocateDSI( MTMSRA_GL,       MLTMSRM, 0)
  call AllocateDSI( MTMSRUE_GL,      MLTMSRM, 0)
  call AllocateDSI( MTMSRUT_GL,      MLTMSRM, 0)
  call AllocateDSI( MTMSRU_GL,       MLTMSRM, 0)
  call AllocateDSI( MTMSRQE_GL,      MLTMSRM, 0)
  call AllocateDSI( MTMSRQ_GL,       MLTMSRM, 0)
  call AllocateDSI( MLTM_GL,         MLTMSRM, 0)

  call AllocateDSI( NLOS_Global,     NBBSM, KCM, NSTVM2,   0)
  call AllocateDSI( NLOE_Global,     NBBEM, KCM, NSTVM2,   0)
  call AllocateDSI( NLOW_Global,     NBBWM, KCM, NSTVM2,   0)
  call AllocateDSI( NLON_Global,     NBBNM, KCM, NSTVM2,   0)

  call AllocateDSI( CLOS_Global,     NBBSM, KCM, NSTVM2, 0.0)
  call AllocateDSI( CLOE_Global,     NBBEM, KCM, NSTVM2, 0.0)
  call AllocateDSI( CLOW_Global,     NBBWM, KCM, NSTVM2, 0.0)
  call AllocateDSI( CLON_Global,     NBBNM, KCM, NSTVM2, 0.0)

  call AllocateDSI( CUE_Global,      LCM_Global, 0.0)
  call AllocateDSI( CUN_Global,      LCM_Global, 0.0)
  call AllocateDSI( CVN_Global,      LCM_Global, 0.0)
  call AllocateDSI( CVE_Global,      LCM_Global, 0.0)

  call AllocateDSI( DZC_Global,      LCM_Global, KCM, 0.0)

  ! *** only allocate if netcdf options turned on
  if(NCDFOUT > 0 )then
    call AllocateDSI( TAUBSED_Global,  LCM_GLOBAL, 0.0)
    call AllocateDSI( TAUBSND_Global,  LCM_GLOBAL, 0.0)
    call AllocateDSI( TAUB_Global,     LCM_GLOBAL, 0.0)
    call AllocateDSI( WNDVELE_Global,  LCM_GLOBAL, 0.0)
    call AllocateDSI( WNDVELN_GLobal,  LCM_GLOBAL, 0.0)
    call AllocateDSI( PATMT_Global,    LCM_GLOBAL, 1010.)
  endif

  If(MPI_Write_Flag )then
    write(mpi_log_unit, '(a)') 'Array sizes for domain decomposition:'
    write(mpi_log_unit, '(a, I6)') 'NPBEM    = ' , NPBEM
    write(mpi_log_unit, '(a, I6)') 'NPBNM    = ' , NPBNM
    write(mpi_log_unit, '(a, I6)') 'NPBSM    = ' , NPBSM
    write(mpi_log_unit, '(a, I6)') 'NQSIJM   = ' , NQSIJM
    write(mpi_log_unit, '(a, I6)') 'NSTVM    = ' , NSTVM
    write(mpi_log_unit, '(a, I6)') 'NBBNM    = ' , NBBNM
    write(mpi_log_unit, '(a, I6)') 'MLTMSRM  = ' , MLTMSRM
  Endif

  write(mpi_log_unit, '(a)') 'Finished allocating variables for domain decomp'
  call WriteBreak(mpi_log_unit)

  End subroutine Allocate_Domain_Decomp
