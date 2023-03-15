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

  Use GLOBAL
  Use Allocate_Initialize
  Use Variables_WQ
  Use WQ_RPEM_MODULE
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  
  Implicit none

  ! *** Local variables
  Integer :: L

  ! *** Special Range Cases
  Allocate(IL2IG(-2:global_max_width_x+2))
  Allocate(JL2JG(-2:global_max_width_y+2))
  Allocate(IG2IL(-2:global_max_width_x+2))
  Allocate(JG2JL(-2:global_max_width_y+2))
  IL2IG = 0
  JL2JG = 0
  IG2IL = 0
  JG2JL = 0

  ! *** Character arrays
  Allocate(CLTMSR_GL(MLTMSRM))
  Allocate(GRPID_GL(NQSIJM))
  
  ! *** Allocate general MPI arrays and initialize
  Call AllocateDSI( lij_west_conn_outside, NPEWBP,  0)
  Call AllocateDSI( lij_east_conn_outside, NPEWBP,  0)

  Call AllocateDSI( GWCSER_Global,   NDGWSER, NGWSERM, NSTVM2, 0.0)
  Call AllocateDSI( NGWSL_Global,    LCM_Global,   0)

  Call AllocateDSI( ISCDRY_Global,   LCM_Global,   0)
  Call AllocateDSI( NATDRY_Global,   LCM_Global,   0)
  Call AllocateDSI( IDRY_Global,     LCM_Global,   0)

  Call AllocateDSI( BELV_Global,     LCM_Global, 0.0)
  Call AllocateDSI( DXP_Global,      LCM_Global, 0.0)
  Call AllocateDSI( DYP_Global,      LCM_Global, 0.0)
  Call AllocateDSI( HP_Global,       LCM_Global, 0.0)
  Call AllocateDSI( H1P_Global,      LCM_Global, 0.0)
  Call AllocateDSI( HWQ_Global,      LCM_Global, 0.0)
  Call AllocateDSI( H2WQ_Global,     LCM_Global, 0.0)
  Call AllocateDSI( ZBR_Global,      LCM_Global, 0.0)
  Call AllocateDSI( MVEG_Global,     LCM_Global,   0)

  Call AllocateDSI( SHEAR_Global,    LCM_Global, 0.0)
  Call AllocateDSI( SHEAR_Global2,   LCM_Global, 0.0)
  Call AllocateDSI( SHEAR_Local,     LCM, 0.0)
  Call AllocateDSI( SHEAR_Local2,    LCM, 0.0)

  Call AllocateDSI( EVAPSW_Global,   LCM_Global, 0.0)
  Call AllocateDSI( EVAPGW_Global,   LCM_Global, 0.0)
  Call AllocateDSI( QGW_Global,      LCM_Global, 0.0)
  Call AllocateDSI( AGWELV_Global,   LCM_Global, 0.0)

  Call AllocateDSI( EVAPT_Global,    LCM_Global, 0.0)
  Call AllocateDSI( RAINT_Global,    LCM_Global, 0.0)

  Call AllocateDSI( RSSBCE_Global,   LCM_Global, 0.0)
  Call AllocateDSI( RSSBCW_Global,   LCM_Global, 0.0)
  Call AllocateDSI( RSSBCN_Global,   LCM_Global, 0.0)
  Call AllocateDSI( RSSBCS_Global,   LCM_Global, 0.0)

  Call AllocateDSI( UHDYE_Global,    LCM_Global, 0.0)
  Call AllocateDSI( UHDY1E_Global,   LCM_Global, 0.0)
  Call AllocateDSI( VHDXE_Global,    LCM_Global, 0.0)
  Call AllocateDSI( VHDX1E_Global,   LCM_Global, 0.0)
  Call AllocateDSI( SUB_Global,      LCM_Global, 0.0)
  Call AllocateDSI( SVB_Global,      LCM_Global, 0.0)

  Call AllocateDSI( U_Global,        LCM_Global, KCM, 0.0)
  Call AllocateDSI( V_Global,        LCM_Global, KCM, 0.0)
  Call AllocateDSI( U1_Global,       LCM_Global, KCM, 0.0)
  Call AllocateDSI( V1_Global,       LCM_Global, KCM, 0.0)
  Call AllocateDSI( W_Global,        LCM_Global, KCM, 0.0)

  Call AllocateDSI( QQ_Global,       LCM_Global, -KCM, 0.0)
  Call AllocateDSI( QQ1_Global,      LCM_Global, -KCM, 0.0)
  Call AllocateDSI( QQL_Global,      LCM_Global, -KCM, 0.0)
  Call AllocateDSI( QQL1_Global,     LCM_Global, -KCM, 0.0)
  Call AllocateDSI( DML_Global,      LCM_Global, -KCM, 0.0)

  Call AllocateDSI( QSUM_Global,     LCM_Global, KCM, 0.0)
  Call AllocateDSI( QSUME_Global,    LCM_Global, 0.0)

  Call AllocateDSI( VHDX2_Global,    LCM_Global, KCM, 0.0)
  Call AllocateDSI( UHDY2_Global,    LCM_Global, KCM, 0.0)

  IF( ISTRAN(1) > 0 )THEN
    Call AllocateDSI( SAL_Global,      LCM_Global, KCM, 0.0)
    Call AllocateDSI( SAL1_Global,     LCM_Global, KCM, 0.0)
  ENDIF

  IF( ISTRAN(2) > 0 )THEN
    Call AllocateDSI( TEM_Global,      LCM_Global, KCM, 0.0)
    Call AllocateDSI( TEM1_Global,     LCM_Global, KCM, 0.0)
    Call AllocateDSI( TEMB_Global,     LCM_Global, 0.0)
    Call AllocateDSI( SHAD_Global,     LCM_Global, 0.0)
    IF( ISICE > 0 )THEN
      Call AllocateDSI( ICETHICK_Global, LCM_Global, 0.0)
      Call AllocateDSI( ICETEMP_Global,  LCM_Global, 0.0)
    ENDIF
  ENDIF
  IF( ISTRAN(3) > 0 )THEN
    Call AllocateDSI( DYE_Global,      LCM_Global, KCM, NDYM, 0.0)
    Call AllocateDSI( DYE1_Global,     LCM_Global, KCM, NDYM, 0.0)
  ENDIF
  IF( ISTRAN(4) > 0 )THEN
    Call AllocateDSI( SFL_Global,      LCM_Global, KCM, 0.0)
  ENDIF

  IF( ISTRAN(5) > 0 )THEN
    Call AllocateDSI( TOX_Global,      LCM_Global, KCM, NTXM, 0.0)
    Call AllocateDSI( TOX1_Global,     LCM_Global, KCM, NTXM, 0.0)
    Call AllocateDSI( TOXB_Global,     LCM_Global, KBM, NTXM, 0.0)
    Call AllocateDSI( TOXB1_Global,    LCM_Global, KBM, NTXM, 0.0)
  ENDIF

  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    Call AllocateDSI( HBED_Global,     LCM_Global, KBM, 0.0)
    Call AllocateDSI( HBED1_Global,    LCM_Global, KBM, 0.0)

    Call AllocateDSI( BEDMAP_Global,   LCM_Global, 0)
    Call AllocateDSI( KBT_Global,      LCM_Global, 0)
    
    Call AllocateDSI( BDENBED_Global,  LCM_Global, KBM, 0.0)
    Call AllocateDSI( PORBED_Global,   LCM_Global, KBM, 0.0)
    Call AllocateDSI( VDRBED_Global,   LCM_Global, KBM, 0.0)
    Call AllocateDSI( VDRBED1_Global,  LCM_Global, KBM, 0.0)

    IF( ISTRAN(7) > 0 .AND. ICALC_BL > 0 .AND. NSND > 0 )THEN
      Call AllocateDSI( QSBDLDX_Global,  LCM_Global, NSND, 0.0)
      Call AllocateDSI( QSBDLDY_Global,  LCM_Global, NSND, 0.0)
    ELSEIF( ICALC_BL > 0 )THEN
      Call AllocateDSI( QSBDLDX_Global,  LCM_Global, NSCM, 0.0)
      Call AllocateDSI( QSBDLDY_Global,  LCM_Global, NSCM, 0.0)
    ENDIF

    IF( NSEDFLUME > 0 )THEN
      Call AllocateDSI( LAYERACTIVE_Global, KB, LCM_Global, 0)
      Call AllocateDSI( ETOTO_Global,    LCM_Global, 0.0)
      Call AllocateDSI( DEPO_Global,     LCM_Global, 0.0)
      Call AllocateDSI( TAU_Global,      LCM_Global, 0.0)
      Call AllocateDSI( D50AVG_Global,   LCM_Global, 0.0)
      Call AllocateDSI( BULKDENS_Global, KB, LCM_Global, 0.0)
      Call AllocateDSI( TSED_Global,     KB, LCM_Global, 0.0)
      Call AllocateDSI( PERSED_Global,   NSCM, KB, LCM_Global, 0.0)
      IF( ICALC_BL > 0 )THEN
        Call AllocateDSI( CBL_Global,    LCM_Global, NSCM, 0.0)
        IF( ISTRAN(5) > 0 )THEN
          Call AllocateDSI( CBLTOX_Global, LCM_Global, NTXM, 0.0)
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  
  IF( ISTRAN(6) > 0 )THEN
    Call AllocateDSI( SED_Global,      LCM_Global, KCM, NSCM2, 0.0)
    Call AllocateDSI( SED1_Global,     LCM_Global, KCM, NSCM2, 0.0)
    Call AllocateDSI( SEDB_Global,     LCM_Global, KBM, NSCM,  0.0)
    Call AllocateDSI( SEDB1_Global,    LCM_Global, KBM, NSCM,  0.0)
    ! *** Propwash fast settling
    !Call AllocateDSI( SDF_Global,      LCM_Global, KCM, NSCM, 0.0)
  ENDIF
  IF( ISTRAN(7) > 0 )THEN
    Call AllocateDSI( SND_Global,      LCM_Global, KCM, NSNM, 0.0)
    Call AllocateDSI( SND1_Global,     LCM_Global, KCM, NSNM, 0.0)
    Call AllocateDSI( SNDB_Global,     LCM_Global, KBM, NSNM, 0.0)
    Call AllocateDSI( SNDB1_Global,    LCM_Global, KBM, NSNM, 0.0)
  ENDIF

  ! *** Water Quality
  IF( ISTRAN(8) > 0 )then  
    Call AllocateDSI( WQV_Global,         LCM_Global, KCM, -NWQVM, 0.0)
    
    IF( ISRPEM > 0 )THEN
      Call AllocateDSI( WQRPS_Global,     LCM_Global, 0.0)
      Call AllocateDSI( WQRPR_Global,     LCM_Global, 0.0)
      Call AllocateDSI( WQRPE_Global,     LCM_Global, 0.0)
      Call AllocateDSI( WQRPD_Global,     LCM_Global, 0.0)
      Call AllocateDSI( LMASKRPEM_Global, LCM_Global, .false.)
    End if

    ! *** Sediment Diagenesis
    Call AllocateDSI( SMPON_Global,    LCM_Global, NSMGM, 0.0 )
    Call AllocateDSI( SMPOP_Global,    LCM_Global, NSMGM, 0.0 )
    Call AllocateDSI( SMPOC_Global,    LCM_Global, NSMGM, 0.0 )
    Call AllocateDSI( SMDFN_Global,    LCM_Global, NSMGM, 0.0 )
    Call AllocateDSI( SMDFP_Global,    LCM_Global, NSMGM, 0.0 )
    Call AllocateDSI( SMDFC_Global,    LCM_Global, NSMGM, 0.0 )

    Call AllocateDSI( SM1NH4_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM2NH4_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM1NO3_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM2NO3_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM1PO4_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM2PO4_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM1H2S_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM2H2S_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SM1SI_Global,    LCM_Global, 0.0 )
    Call AllocateDSI( SM2SI_Global,    LCM_Global, 0.0 )
    Call AllocateDSI( SMPSI_Global,    LCM_Global, 0.0 )
    Call AllocateDSI( SMBST_Global,    LCM_Global, 0.0 )
    Call AllocateDSI( SMT_Global,      LCM_Global, 0.0 )
    Call AllocateDSI( SMCSOD_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( SMNSOD_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( WQBFNH4_Global,  LCM_Global, 0.0 )
    Call AllocateDSI( WQBFNO3_Global,  LCM_Global, 0.0 )
    Call AllocateDSI( WQBFO2_Global,   LCM_Global, 0.0 )
    Call AllocateDSI( WQBFCOD_Global,  LCM_Global, 0.0 )
    Call AllocateDSI( WQBFPO4D_Global, LCM_Global, 0.0 )
    Call AllocateDSI( WQBFSAD_Global,  LCM_Global, 0.0 )
  ENDIF
  
  Call AllocateDSI( QQWV3_Global,    LCM_Global, 0.0)
  Call AllocateDSI( WVHUU_Global,    LCM_Global, KCM, 0.0)
  Call AllocateDSI( WVHVV_Global,    LCM_Global, KCM, 0.0)
  Call AllocateDSI( WVHUV_Global,    LCM_Global, KCM, 0.0)

  Call AllocateDSI( KSZ_Global,      LCM_Global,   1)

  Call AllocateDSI( LWC_Global,      LCM_Global,    0)
  Call AllocateDSI( LEC_Global,      LCM_Global,    0)
  Call AllocateDSI( LSC_Global,      LCM_Global,    0)
  Call AllocateDSI( LNC_Global,      LCM_Global,    0)
  
  Call AllocateDSI( UMASK_Global,    -LCM_Global, 0.0)
  Call AllocateDSI( VMASK_Global,    -LCM_Global, 0.0)
  
  IF( ISWAVE > 0 )THEN
  
    Call AllocateDSI( FXWAVE_Global,    LCM_Global, KCM, 0.0)
    Call AllocateDSI( FYWAVE_Global,    LCM_Global, KCM, 0.0)
    Call AllocateDSI( WV_HEIGHT_Global, LCM_Global, 0.0)
    Call AllocateDSI( WV_FREQ_Global,   LCM_Global, 0.0)
    Call AllocateDSI( WV_DIR_Global,    LCM_Global, 0.0)
    Call AllocateDSI( WV_DISSIPA_Global,LCM_Global, 0.0)
                                       
    Call AllocateDSI( TBX_Global,       LCM_Global, 0.0) 
    Call AllocateDSI( TBY_Global,       LCM_Global, 0.0)
    
    Allocate(WV_Global(LCM_Global))
    DO L=1, LCM_Global
      Call AllocateDSI( WV_Global(L).DISSIPA, KCM, 0.0)
    ENDDO
  ENDIF

  Call AllocateDSI( IJCT_GLOBAL,   global_max_width_x, global_max_width_y, 0)
  Call AllocateDSI( IJCTLT_GLOBAL, global_max_width_x, global_max_width_y, 0)
  Call AllocateDSI( LIJ_Global,    global_max_width_x, global_max_width_y, 0)

  ! *** Open Boundary Conditions
  Call AllocateDSI( IPBS_GL,         NPBSM,   0)
  Call AllocateDSI( IPBW_GL,         NPBWM,   0)
  Call AllocateDSI( IPBE_GL,         NPBEM,   0)
  Call AllocateDSI( IPBN_GL,         NPBNM,   0)
                                             
  Call AllocateDSI( JPBS_GL,         NPBSM,   0)
  Call AllocateDSI( JPBW_GL,         NPBWM,   0)
  Call AllocateDSI( JPBE_GL,         NPBEM,   0)
  Call AllocateDSI( JPBN_GL,         NPBNM,   0)

  Call AllocateDSI( LPBS_GL,         NPBSM,   0)
  Call AllocateDSI( LPBW_GL,         NPBWM,   0)
  Call AllocateDSI( LPBE_GL,         NPBEM,   0)
  Call AllocateDSI( LPBN_GL,         NPBNM,   0)

  Call AllocateDSI( ISPBS_GL,        NPBSM,   0)
  Call AllocateDSI( ISPBW_GL,        NPBWM,   0)
  Call AllocateDSI( ISPBE_GL,        NPBEM,   0)
  Call AllocateDSI( ISPBN_GL,        NPBNM,   0)
                                             
  Call AllocateDSI( ISPRS_GL,        NPBSM,   0)
  Call AllocateDSI( ISPRW_GL,        NPBWM,   0)
  Call AllocateDSI( ISPRE_GL,        NPBEM,   0)
  Call AllocateDSI( ISPRN_GL,        NPBNM,   0)
                                             
  Call AllocateDSI( NPSERS_GL,       NPBSM,   0)
  Call AllocateDSI( NPSERW_GL,       NPBWM,   0)
  Call AllocateDSI( NPSERE_GL,       NPBEM,   0)
  Call AllocateDSI( NPSERN_GL,       NPBNM,   0)
                                             
  Call AllocateDSI( NPSERS1_GL,      NPBSM,   0)
  Call AllocateDSI( NPSERW1_GL,      NPBWM,   0)
  Call AllocateDSI( NPSERE1_GL,      NPBEM,   0)
  Call AllocateDSI( NPSERN1_GL,      NPBNM,   0)

  Call AllocateDSI( TPCOORDS_GL,     NPBSM,   0)
  Call AllocateDSI( TPCOORDW_GL,     NPBWM,   0)
  Call AllocateDSI( TPCOORDE_GL,     NPBEM,   0)
  Call AllocateDSI( TPCOORDN_GL,     NPBNM,   0)

  Call AllocateDSI( PCBS_GL,         NPBSM, MTM, 0.0)
  Call AllocateDSI( PCBW_GL,         NPBWM, MTM, 0.0)
  Call AllocateDSI( PCBE_GL,         NPBEM, MTM, 0.0)
  Call AllocateDSI( PCBN_GL,         NPBNM, MTM, 0.0)

  Call AllocateDSI( PSBS_GL,         NPBSM, MTM, 0.0)
  Call AllocateDSI( PSBW_GL,         NPBWM, MTM, 0.0)
  Call AllocateDSI( PSBE_GL,         NPBEM, MTM, 0.0)
  Call AllocateDSI( PSBN_GL,         NPBNM, MTM, 0.0)

  ! *** River global data arrays
  Call AllocateDSI( IQS_GL,          NQSIJM,   0)
  Call AllocateDSI( JQS_GL,          NQSIJM,   0)

  Call AllocateDSI( QSSE_GL,         NQSIJM, 0.0)
  Call AllocateDSI( NQSMUL_GL,       NQSIJM,   0)
  Call AllocateDSI( NQSMF_GL,        NQSIJM,   0)
  Call AllocateDSI( NQSERQ_GL,       NQSIJM,   0)
  Call AllocateDSI( NCSERQ_GL,       NQSIJM, NSTVM2, 0)
  Call AllocateDSI( QWIDTH_GL,       NQSIJM, 0.0)
  Call AllocateDSI( QFACTOR_GL,      NQSIJM, 0.0)
  Call AllocateDSI( LQS_GL,          NQSIJM,   0)
  Call AllocateDSI( CQSE_GL,         NQSIJM, NSTVM2, 0.0)

  ! *** Concentration global data arrays
  Call AllocateDSI( ICBS_GL,         NBBSM,   0)
  Call AllocateDSI( JCBS_GL,         NBBSM,   0)
  Call AllocateDSI( NTSCRS_GL,       NBBSM,   0)
  Call AllocateDSI( NCSERS_GL,       NBBSM,    NSTVM2,   0)
  Call AllocateDSI( CBS_GL,          NBBSM, 2, NSTVM2, 0.0)
  Call AllocateDSI( ICBW_GL,         NBBWM,   0)
  Call AllocateDSI( JCBW_GL,         NBBWM,   0)
  Call AllocateDSI( NTSCRW_GL,       NBBWM,   0)
  Call AllocateDSI( NCSERW_GL,       NBBWM,    NSTVM2,   0)
  Call AllocateDSI( CBW_GL,          NBBWM, 2, NSTVM2, 0.0)
  Call AllocateDSI( ICBE_GL,         NBBEM,   0)
  Call AllocateDSI( JCBE_GL,         NBBEM,   0)
  Call AllocateDSI( NTSCRE_GL,       NBBEM,   0)
  Call AllocateDSI( NCSERE_GL,       NBBEM,    NSTVM2,   0)
  Call AllocateDSI( CBE_GL,          NBBEM, 2, NSTVM2, 0.0)
  Call AllocateDSI( ICBN_GL,         NBBNM,   0)
  Call AllocateDSI( JCBN_GL,         NBBNM,   0)
  Call AllocateDSI( NTSCRN_GL,       NBBNM,   0)
  Call AllocateDSI( NCSERN_GL,       NBBNM,    NSTVM2,   0)
  Call AllocateDSI( CBN_GL,          NBBNM, 2, NSTVM2, 0.0)
  Call AllocateDSI( ILTMSR_GL,       MLTMSRM, 0)
  Call AllocateDSI( JLTMSR_GL,       MLTMSRM, 0)
  Call AllocateDSI( NTSSSS_GL,       MLTMSRM, 0)

  Call AllocateDSI( MTMSRP_GL,       MLTMSRM, 0)
  Call AllocateDSI( MTMSRC_GL,       MLTMSRM, 0)
  Call AllocateDSI( MTMSRA_GL,       MLTMSRM, 0)
  Call AllocateDSI( MTMSRUE_GL,      MLTMSRM, 0)
  Call AllocateDSI( MTMSRUT_GL,      MLTMSRM, 0)
  Call AllocateDSI( MTMSRU_GL,       MLTMSRM, 0)
  Call AllocateDSI( MTMSRQE_GL,      MLTMSRM, 0)
  Call AllocateDSI( MTMSRQ_GL,       MLTMSRM, 0)
  Call AllocateDSI( MLTM_GL,         MLTMSRM, 0)

  Call AllocateDSI( NLOS_Global,     NBBSM, KCM, NSTVM2,   0)
  Call AllocateDSI( NLOE_Global,     NBBEM, KCM, NSTVM2,   0)
  Call AllocateDSI( NLOW_Global,     NBBWM, KCM, NSTVM2,   0)
  Call AllocateDSI( NLON_Global,     NBBNM, KCM, NSTVM2,   0)

  Call AllocateDSI( CLOS_Global,     NBBSM, KCM, NSTVM2, 0.0)
  Call AllocateDSI( CLOE_Global,     NBBEM, KCM, NSTVM2, 0.0)
  Call AllocateDSI( CLOW_Global,     NBBWM, KCM, NSTVM2, 0.0)
  Call AllocateDSI( CLON_Global,     NBBNM, KCM, NSTVM2, 0.0)

  ! *** Jet/Plume global mapping values, from C27 in efdc.inp
  IF( NQJPIJ > 0 )THEN
    Call AllocateDSI( ICALJP_GL,     NQJPM,   0)
    Call AllocateDSI( IQJP_GL,       NQJPM,   0)
    Call AllocateDSI( JQJP_GL,       NQJPM,   0)
    Call AllocateDSI( KQJP_GL,       NQJPM,   0)
    Call AllocateDSI( NPORTJP_GL,    NQJPM,   0)
      
    Call AllocateDSI( XJETL_GL,      NQJPM, 0.0)
    Call AllocateDSI( YJETL_GL,      NQJPM, 0.0)
    Call AllocateDSI( ZJET_GL,       NQJPM, 0.0)
    Call AllocateDSI( PHJET_GL,      NQJPM, 0.0)
    Call AllocateDSI( THJET_GL,      NQJPM, 0.0)
    Call AllocateDSI( DJET_GL,       NQJPM, 0.0)
    Call AllocateDSI( CFRD_GL,       NQJPM, 0.0)
    Call AllocateDSI( DJPER_GL,      NQJPM, 0.0)
  ENDIF

  Call AllocateDSI( CUE_Global,      LCM_Global, 0.0) 
  Call AllocateDSI( CUN_Global,      LCM_Global, 0.0) 
  Call AllocateDSI( CVN_Global,      LCM_Global, 0.0) 
  Call AllocateDSI( CVE_Global,      LCM_Global, 0.0) 

  ! *** only allocate if netcdf options turned on
  if(NCDFOUT > 0 )then
      Call AllocateDSI( DZC_Global,      LCM_Global, KCM, 0.0) 
      Call AllocateDSI( TAUBSED_Global,  LCM_GLOBAL, 0.0) 
      Call AllocateDSI( TAUBSND_Global,  LCM_GLOBAL, 0.0) 
      Call AllocateDSI( TAUB_Global,     LCM_GLOBAL, 0.0) 
      Call AllocateDSI( WNDVELE_Global,  LCM_GLOBAL, 0.0) 
      Call AllocateDSI( WNDVELN_GLobal,  LCM_GLOBAL, 0.0) 
      Call AllocateDSI( PATMT_Global,    LCM_GLOBAL, 1010.)
  endif

  If(MPI_Write_Flag == .TRUE. )THEN
    write(mpi_log_unit, '(a)') 'Array sizes for domain decomposition:'
    write(mpi_log_unit, '(a, I6)') 'NPBEM   =' , NPBEM
    write(mpi_log_unit, '(a, I6)') 'NPBNM   =' , NPBNM
    write(mpi_log_unit, '(a, I6)') 'NPBSM   =' , NPBSM
    write(mpi_log_unit, '(a, I6)') 'NQSIJM  =' , NQSIJM
    write(mpi_log_unit, '(a, I6)') 'NSTVM   =' , NSTVM
    write(mpi_log_unit, '(a, I6)') 'NBBNM   =' , NBBNM
    write(mpi_log_unit, '(a, I6)') 'MLTMSRM =' , MLTMSRM
  Endif

  write(mpi_log_unit, '(a)') 'Finished allocating variables for domain decomp'
  Call WriteBreak(mpi_log_unit)

  End subroutine Allocate_Domain_Decomp
