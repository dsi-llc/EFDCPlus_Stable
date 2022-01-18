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
!    2019-09      Zander   Allocates variables required for domain decomp
!
!---------------------------------------------------------------------------!
Subroutine Allocate_Domain_Decomp

  Use GLOBAL
  Use Variables_WQ
  Use WQ_RPEM_MODULE
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  
  Implicit none

  ! *** Read in variables

  ! *** Local variables
  Integer :: L

  
  Allocate(lij_west_conn_outside(NPEWBP))
  Allocate(lij_east_conn_outside(NPEWBP))

  Allocate(GWCSER_Global(NDGWSER,NGWSERM,NSTVM))
  Allocate(NGWSL_Global(LCM_Global))

  Allocate(ISCDRY_Global(LCM_Global))
  Allocate(NATDRY_Global(LCM_Global))
  Allocate(IDRY_Global(LCM_Global))

  Allocate(BELV_Global(LCM_Global))
  Allocate(DXP_Global(LCM_Global))
  Allocate(DYP_Global(LCM_Global))
  Allocate(HP_Global(LCM_Global))
  Allocate(H1P_Global(LCM_Global))
  Allocate(HWQ_Global(LCM_Global))
  Allocate(H2WQ_Global(LCM_Global))

  Allocate(SHEAR_Global(LCM_Global))
  Allocate(SHEAR_Global2(LCM_Global))
  Allocate(SHEAR_Local(LCM))
  Allocate(SHEAR_Local2(LCM))

  Allocate(EVAPSW_Global(LCM_Global))
  Allocate(EVAPGW_Global(LCM_Global))
  Allocate(QGW_Global(LCM_Global))
  Allocate(AGWELV_Global(LCM_Global))

  Allocate(EVAPT_Global(LCM_Global))
  Allocate(RAINT_Global(LCM_Global))

  Allocate(RSSBCE_Global(LCM_Global))
  Allocate(RSSBCW_Global(LCM_Global))
  Allocate(RSSBCN_Global(LCM_Global))
  Allocate(RSSBCS_Global(LCM_Global))

  Allocate(UHDYE_Global(LCM_Global))
  Allocate(UHDY1E_Global(LCM_Global))
  Allocate(VHDXE_Global(LCM_Global))
  Allocate(VHDX1E_Global(LCM_Global))
  Allocate(SUB_Global(LCM_Global))
  Allocate(SVB_Global(LCM_Global))

  Allocate(U_Global(LCM_Global, KCM))
  Allocate(V_Global(LCM_Global, KCM))
  Allocate(U1_Global(LCM_Global, KCM))
  Allocate(V1_Global(LCM_Global, KCM))
  Allocate(W_Global(LCM_Global, KCM))

  Allocate(QQ_Global(LCM_Global, 0:KCM))
  Allocate(QQ1_Global(LCM_Global, 0:KCM))
  Allocate(QQL_Global(LCM_Global, 0:KCM))
  Allocate(QQL1_Global(LCM_Global, 0:KCM))
  Allocate(DML_Global(LCM_Global, 0:KCM))

  Allocate(QSUM_Global(LCM_Global, KCM))
  Allocate(QSUME_Global(LCM_Global))

  Allocate(VHDX2_Global(LCM_Global, KCM))
  Allocate(UHDY2_Global(LCM_Global, KCM))

  IF( ISTRAN(1) > 0 )THEN
    Allocate(SAL_Global(LCM_Global, KCM))
    Allocate(SAL1_Global(LCM_Global, KCM))
  ENDIF

  IF( ISTRAN(2) > 0 )THEN
    Allocate(TEM_Global(LCM_Global, KCM))
    Allocate(TEM1_Global(LCM_Global, KCM))
    Allocate(TEMB_Global(LCM_Global))
    Allocate(SHAD_Global(LCM_Global))
    IF( ISICE > 0 )THEN
      Allocate(ICETHICK_Global(LCM_Global))
      Allocate(ICETEMP_Global(LCM_Global))
    ENDIF
  ENDIF
  IF( ISTRAN(3) > 0 )THEN
    Allocate(DYE_Global(LCM_Global, KCM, NDYM))
    Allocate(DYE1_Global(LCM_Global, KCM, NDYM))
  ENDIF
  IF( ISTRAN(4) > 0 )THEN
    Allocate(SFL_Global(LCM_Global, KCM))
  ENDIF

  IF( ISTRAN(5) > 0 )THEN
    Allocate(TOX_Global(LCM_Global, KCM, NTXM))
    Allocate(TOX1_Global(LCM_Global, KCM, NTXM))
    Allocate(TOXB_Global(LCM_Global, KBM, NTXM))
    Allocate(TOXB1_Global(LCM_Global, KBM, NTXM))
  ENDIF

  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    Allocate(HBED_Global(LCM_Global, KBM))
    Allocate(HBED1_Global(LCM_Global, KBM))

    Allocate(BEDMAP_Global(LCM_Global))
    Allocate(KBT_Global(LCM_Global))
    
    Allocate(BDENBED_Global(LCM_Global, KBM))
    Allocate(PORBED_Global(LCM_Global, KBM))
    Allocate(VDRBED_Global(LCM_Global, KBM))
    Allocate(VDRBED1_Global(LCM_Global, KBM))

    IF( ISTRAN(7) > 0 .AND. ICALC_BL > 0 .AND. NSND > 0 )THEN
      Allocate(QSBDLDX_Global(LCM_Global, NSND))
      Allocate(QSBDLDY_Global(LCM_Global, NSND))
    ELSEIF( ICALC_BL > 0 )THEN
      Allocate(QSBDLDX_Global(LCM_Global, NSCM))
      Allocate(QSBDLDY_Global(LCM_Global, NSCM))
    ENDIF

    IF( NSEDFLUME > 0 )THEN
      Allocate(LAYERACTIVE_Global(KB,LCM_Global))
      Allocate(TAU_Global(LCM_Global))
      Allocate(D50AVG_Global(LCM_Global))
      Allocate(BULKDENS_Global(KB,LCM_Global))
      Allocate(TSED_Global(KB,LCM_Global))
      Allocate(PERSED_Global(NSCM,KB,LCM_Global))
      Allocate(ETOTO_Global(LCM_Global))
      Allocate(DEPO_Global(LCM_Global))
      IF( ICALC_BL > 0 )THEN
        Allocate(CBL_Global(LCM_Global, NSCM))
        IF( ISTRAN(5) > 0 )THEN
          ALLOCATE(CBLTOX_Global(LCM_Global,NTXM))
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  
  IF( ISTRAN(6) > 0 )THEN
    Allocate(SED_Global(LCM_Global, KCM, NSCM))
    Allocate(SED1_Global(LCM_Global, KCM, NSCM))
    Allocate(SEDB_Global(LCM_Global, KBM, NSCM))
    Allocate(SEDB1_Global(LCM_Global, KBM, NSCM))
  ENDIF
  IF( ISTRAN(7) > 0 )THEN
    Allocate(SND_Global(LCM_Global, KCM, NSNM))
    Allocate(SND1_Global(LCM_Global, KCM, NSNM))
    Allocate(SNDB_Global(LCM_Global, KBM, NSNM))
    Allocate(SNDB1_Global(LCM_Global, KBM, NSNM))
  ENDIF

  ! *** Water Quality
  IF( ISTRAN(8) > 0 )then
    ALLOCATE(WQV_Global(LCM_Global, KCM, 0:NWQVM))
    
    IF( ISRPEM > 0 )THEN
      Allocate(WQRPS_Global(LCM_Global))
      Allocate(WQRPR_Global(LCM_Global))
      Allocate(WQRPE_Global(LCM_Global))
      Allocate(WQRPD_Global(LCM_Global))
      Allocate(LMASKRPEM_Global(LCM_Global))
    End if


    ! *** Sediment Diagenesis
    Allocate(SMPON_Global(LCM_Global,NSMGM) )
    Allocate(SMPOP_Global(LCM_Global,NSMGM) )
    Allocate(SMPOC_Global(LCM_Global,NSMGM) )
    Allocate(SMDFN_Global(LCM_Global,NSMGM) )
    Allocate(SMDFP_Global(LCM_Global,NSMGM) )
    Allocate(SMDFC_Global(LCM_Global,NSMGM) )

    Allocate(SM1NH4_Global  (LCM_Global) )
    Allocate(SM2NH4_Global  (LCM_Global) )
    Allocate(SM1NO3_Global  (LCM_Global) )
    Allocate(SM2NO3_Global  (LCM_Global) )
    Allocate(SM1PO4_Global  (LCM_Global) )
    Allocate(SM2PO4_Global  (LCM_Global) )
    Allocate(SM1H2S_Global  (LCM_Global) )
    Allocate(SM2H2S_Global  (LCM_Global) )
    Allocate(SM1SI_Global   (LCM_Global) )
    Allocate(SM2SI_Global   (LCM_Global) )
    Allocate(SMPSI_Global   (LCM_Global) )
    Allocate(SMBST_Global   (LCM_Global) )
    Allocate(SMT_Global     (LCM_Global) )
    Allocate(SMCSOD_Global  (LCM_Global) )
    Allocate(SMNSOD_Global  (LCM_Global) )
    Allocate(WQBFNH4_Global (LCM_Global) )
    Allocate(WQBFNO3_Global (LCM_Global) )
    Allocate(WQBFO2_Global  (LCM_Global) )
    Allocate(WQBFCOD_Global (LCM_Global) )
    Allocate(WQBFPO4D_Global(LCM_Global) )
    Allocate(WQBFSAD_Global( LCM_Global) )

  ENDIF
  
  Allocate(KSZ_Global(LCM_Global))

  Allocate(QQWV3_Global(LCM_Global))

  Allocate(WVHUU_Global(LCM_Global,KCM))
  Allocate(WVHVV_Global(LCM_Global,KCM))
  Allocate(WVHUV_Global(LCM_Global,KCM))

  Allocate(LWC_Global(LCM_Global))
  Allocate(LEC_Global(LCM_Global))
  Allocate(LSC_Global(LCM_Global))
  Allocate(LNC_Global(LCM_Global))
  
  Allocate(UMASK_Global(0:LCM_Global))
  Allocate(VMASK_Global(0:LCM_Global))
  
  IF( ISWAVE > 0 )THEN
  
    ALLOCATE(FXWAVE_Global(LCM_Global,KCM))
    ALLOCATE(FYWAVE_Global(LCM_Global,KCM))
    ALLOCATE(WV_Global       (LCM_Global))
    ALLOCATE(WV_HEIGHT_Global(LCM_Global))
    ALLOCATE(WV_FREQ_Global  (LCM_Global))
    ALLOCATE(WV_DIR_Global   (LCM_Global))
    
    Allocate(TBX_Global(LCM_Global)) 
    Allocate(TBY_Global(LCM_Global))
    
    DO L=1,LCM_Global
      ALLOCATE(WV_Global(L).DISSIPA(KCM))
    ENDDO
  ENDIF

  ALLOCATE(IJCT_GLOBAL  (global_max_width_x,  global_max_width_y))
  ALLOCATE(IJCTLT_GLOBAL(global_max_width_x,  global_max_width_y))

  Allocate(LIJ_Global(global_max_width_x, global_max_width_y) )

  ALLOCATE(IL2IG(-2:global_max_width_x+2))
  ALLOCATE(JL2JG(-2:global_max_width_y+2))
  ALLOCATE(IG2IL(-2:global_max_width_x+2))
  ALLOCATE(JG2JL(-2:global_max_width_y+2))

  ! *** Open Boundary Conditions
  ALLOCATE(IPBS_GL(NPBSM))
  ALLOCATE(IPBW_GL(NPBWM))
  ALLOCATE(IPBE_GL(NPBEM))
  ALLOCATE(IPBN_GL(NPBNM))

  ALLOCATE(JPBS_GL(NPBSM))
  ALLOCATE(JPBW_GL(NPBWM))
  ALLOCATE(JPBE_GL(NPBEM))
  ALLOCATE(JPBN_GL(NPBNM))

  ALLOCATE(LPBS_GL(NPBSM))
  ALLOCATE(LPBW_GL(NPBWM))
  ALLOCATE(LPBE_GL(NPBEM))
  ALLOCATE(LPBN_GL(NPBNM))

  ALLOCATE(ISPBS_GL(NPBSM))
  ALLOCATE(ISPBW_GL(NPBWM))
  ALLOCATE(ISPBE_GL(NPBEM))
  ALLOCATE(ISPBN_GL(NPBNM))

  ALLOCATE(NPSERS_GL(NPBSM))
  ALLOCATE(NPSERW_GL(NPBWM))
  ALLOCATE(NPSERE_GL(NPBEM))
  ALLOCATE(NPSERN_GL(NPBNM))

  ALLOCATE(NPSERS1_GL(NPBSM))
  ALLOCATE(NPSERW1_GL(NPBWM))
  ALLOCATE(NPSERE1_GL(NPBEM))
  ALLOCATE(NPSERN1_GL(NPBNM))

  ALLOCATE(TPCOORDS_GL(NPBSM))
  ALLOCATE(TPCOORDW_GL(NPBWM))
  ALLOCATE(TPCOORDE_GL(NPBEM))
  ALLOCATE(TPCOORDN_GL(NPBNM))

  ALLOCATE(PCBS_GL(NPBSM,MTM))
  ALLOCATE(PCBW_GL(NPBWM,MTM))
  ALLOCATE(PCBE_GL(NPBEM,MTM))
  ALLOCATE(PCBN_GL(NPBNM,MTM))

  ALLOCATE(PSBS_GL(NPBSM,MTM))
  ALLOCATE(PSBW_GL(NPBWM,MTM))
  ALLOCATE(PSBE_GL(NPBEM,MTM))
  ALLOCATE(PSBN_GL(NPBNM,MTM))

  ! *** River global data arrays
  ALLOCATE(IQS_GL(NQSIJM))
  ALLOCATE(JQS_GL(NQSIJM))

  ALLOCATE(QSSE_GL(NQSIJM))
  ALLOCATE(NQSMUL_GL(NQSIJM))
  ALLOCATE(NQSMF_GL(NQSIJM))
  ALLOCATE(NQSERQ_GL(NQSIJM))
  ALLOCATE(NCSERQ_GL(NQSIJM,NSTVM))
  ALLOCATE(QWIDTH_GL(NQSIJM))
  ALLOCATE(QFACTOR_GL(NQSIJM))
  ALLOCATE(GRPID_GL(NQSIJM))
  ALLOCATE(LQS_GL(NQSIJM))
  ALLOCATE(CQSE_GL(NQSIJM,NSTVM))

  ! *** Concentration global data arrays
  ALLOCATE(ICBS_GL(NBBSM))
  ALLOCATE(JCBS_GL(NBBSM))
  ALLOCATE(NTSCRS_GL(NBBSM))
  ALLOCATE(NCSERS_GL(NBBSM,NSTVM))
  ALLOCATE(CBS_GL(NBBSM,2,NSTVM))
  ALLOCATE(ICBW_GL(NBBWM))
  ALLOCATE(JCBW_GL(NBBWM))
  ALLOCATE(NTSCRW_GL(NBBWM))
  ALLOCATE(NCSERW_GL(NBBWM,NSTVM))
  ALLOCATE(CBW_GL(NBBWM,2,NSTVM))
  ALLOCATE(ICBE_GL(NBBEM))
  ALLOCATE(JCBE_GL(NBBEM))
  ALLOCATE(NTSCRE_GL(NBBEM))
  ALLOCATE(NCSERE_GL(NBBEM,NSTVM))
  ALLOCATE(CBE_GL(NBBEM,2,NSTVM))
  ALLOCATE(ICBN_GL(NBBNM))
  ALLOCATE(JCBN_GL(NBBNM))
  ALLOCATE(NTSCRN_GL(NBBNM))
  ALLOCATE(NCSERN_GL(NBBNM,NSTVM))
  ALLOCATE(CBN_GL(NBBNM,2,NSTVM))
  ALLOCATE(ILTMSR_GL(MLTMSRM))
  ALLOCATE(JLTMSR_GL(MLTMSRM))
  ALLOCATE(NTSSSS_GL(MLTMSRM))

  ALLOCATE(MTMSRP_GL(MLTMSRM))
  ALLOCATE(MTMSRC_GL(MLTMSRM))
  ALLOCATE(MTMSRA_GL(MLTMSRM))
  ALLOCATE(MTMSRUE_GL(MLTMSRM))
  ALLOCATE(MTMSRUT_GL(MLTMSRM))
  ALLOCATE(MTMSRU_GL(MLTMSRM))
  ALLOCATE(MTMSRQE_GL(MLTMSRM))
  ALLOCATE(MTMSRQ_GL(MLTMSRM))
  ALLOCATE(MLTM_GL(MLTMSRM))
  ALLOCATE(CLTMSR_GL(MLTMSRM))

  ALLOCATE(NLOS_Global(NBBSM,KCM,NSTVM))
  ALLOCATE(NLOE_Global(NBBEM,KCM,NSTVM))
  ALLOCATE(NLOW_Global(NBBWM,KCM,NSTVM))
  ALLOCATE(NLON_Global(NBBNM,KCM,NSTVM))

  ALLOCATE(CLOS_Global(NBBSM,KCM,NSTVM))
  ALLOCATE(CLOE_Global(NBBEM,KCM,NSTVM))
  ALLOCATE(CLOW_Global(NBBWM,KCM,NSTVM))
  ALLOCATE(CLON_Global(NBBNM,KCM,NSTVM))

  ! *** Jet/Plume global mapping values, from C27 in efdc.inp
  IF( NQJPIJ > 0 )THEN
    Allocate(ICALJP_GL(NQJPM))
    Allocate(IQJP_GL(NQJPM))
    Allocate(JQJP_GL(NQJPM))
    Allocate(KQJP_GL(NQJPM))
    Allocate(NPORTJP_GL(NQJPM))
      
    Allocate(XJETL_GL(NQJPM))
    Allocate(YJETL_GL(NQJPM))
    Allocate(ZJET_GL(NQJPM))
    Allocate(PHJET_GL(NQJPM))
    Allocate(THJET_GL(NQJPM))
    Allocate(DJET_GL(NQJPM))
    Allocate(CFRD_GL(NQJPM))
    Allocate(DJPER_GL(NQJPM))
      
    ICALJP_GL = 0
    IQJP_GL   = 0
    JQJP_GL   = 0
    KQJP_GL   = 0
    NPORTJP_GL = 0
  
    XJETL_GL = 0.0
    YJETL_GL = 0.0
    ZJET_GL  = 0.0
    PHJET_GL = 0.0
    THJET_GL = 0.0
    DJET_GL  = 0.0
    CFRD_GL  = 0.0
    DJPER_GL = 0.0
  ENDIF

  ! *** only allocate if netcdf options turned on
      Allocate(CUE_Global(LCM_Global)) 
      Allocate(CUN_Global(LCM_Global)) 
      Allocate(CVN_Global(LCM_Global)) 
      Allocate(CVE_Global(LCM_Global)) 
  if(NCDFOUT > 0 )then
      Allocate(DZC_Global(LCM_Global,KCM)) 
      Allocate(TAUBSED_Global(LCM_GLOBAL)) 
      Allocate(TAUBSND_Global(LCM_GLOBAL)) 
      Allocate(TAUB_Global   (LCM_GLOBAL)) 
      Allocate(WNDVELE_Global(LCM_GLOBAL)) 
      Allocate(WNDVELN_GLobal(LCM_GLOBAL)) 
      Allocate(PATMT_Global(LCM_GLOBAL))
  end if
  
  ! *********************************************************************************************************************
  ! *** Initialize arrays
      CUE_Global      = 0
      CUN_Global      = 0
      CVN_Global      = 0
      CVE_Global      = 0
  if(NCDFOUT > 0 )then    
      DZC_Global      = 0.0
      TAUBSED_Global  = 0.0
      TAUBSND_Global  = 0.0
      TAUB_Global     = 0.0
      WNDVELE_Global  = 0.0
      WNDVELN_GLobal  = 0.0
      PATMT_Global    = 1010.
  end if
  
  IJCT_GLOBAL = 0
  IJCTLT_GLOBAL = 0
  IL2IG = 0
  JL2JG = 0
  IG2IL = 0
  JG2JL = 0

  lij_west_conn_outside = 0
  lij_east_conn_outside = 0

  IPBS_GL   = 0
  IPBW_GL   = 0
  IPBE_GL   = 0
  IPBN_GL   = 0
  JPBS_GL   = 0
  JPBW_GL   = 0
  JPBE_GL   = 0
  JPBN_GL   = 0
  LPBS_GL   = 0
  LPBW_GL   = 0
  LPBE_GL   = 0
  LPBN_GL   = 0
  ISPBS_GL  = 0
  ISPBW_GL  = 0
  ISPBE_GL  = 0
  ISPBN_GL  = 0
  NPSERS_GL = 0
  NPSERW_GL = 0
  NPSERE_GL = 0
  NPSERN_GL = 0
  PCBE_GL   = 0
  PCBN_GL   = 0
  PCBS_GL   = 0
  PCBW_GL   = 0
  PSBE_GL   = 0
  PSBN_GL   = 0
  PSBS_GL   = 0
  PSBW_GL   = 0

  IQS_GL    = 0
  JQS_GL    = 0
  NQSMUL_GL = 0
  NQSMF_GL  = 0
  NQSERQ_GL = 0
  NCSERQ_GL = 0
  GRPID_GL  = 0
  LQS_GL    = 0
  
  QSSE_GL    = 0.0
  QWIDTH_GL  = 0.0
  QFACTOR_GL = 0.0
  CQSE_GL    = 0.0
  
  ICBS_GL   = 0
  JCBS_GL   = 0
  NTSCRS_GL = 0
  NCSERS_GL = 0
  CBS_GL    = 0.0
  
  ICBW_GL   = 0
  JCBW_GL   = 0
  NTSCRW_GL = 0
  NCSERW_GL = 0
  CBW_GL    = 0.0
  
  ICBE_GL   = 0
  JCBE_GL   = 0
  NTSCRE_GL = 0
  NCSERE_GL = 0
  CBE_GL    = 0.0
  
  ICBN_GL   = 0
  JCBN_GL   = 0
  NTSCRN_GL = 0
  NCSERN_GL = 0
  CBN_GL    = 0.0
  
  ILTMSR_GL = 0
  JLTMSR_GL = 0
  NTSSSS_GL = 0

  MTMSRP_GL  = 0
  MTMSRC_GL  = 0
  MTMSRA_GL  = 0
  MTMSRUE_GL = 0
  MTMSRUT_GL = 0
  MTMSRU_GL  = 0
  MTMSRQE_GL = 0
  MTMSRQ_GL  = 0
  MLTM_GL    = 0

  KSZ_Global   = 1
  BELV_Global  = 0.0
  DXP_Global   = 0.0
  DYP_Global   = 0.0
  HP_Global    = 0.0
  H1P_Global   = 0.0
  HWQ_Global   = 0.0
  H2WQ_Global  = 0.0

  U_Global     = 0.0
  V_Global     = 0.0
  W_Global     = 0.0
  QQ_Global    = 0.0
  QQL_Global   = 0.0
  DML_Global   = 0.0

  U1_Global    = 0.0
  V1_Global    = 0.0
  QQ1_Global   = 0.0
  QQL1_Global  = 0.0

  QSUM_Global  = 0.0
  QSUME_Global = 0.0
  VHDX2_Global = 0.0
  UHDY2_Global = 0.0
  SHEAR_Global = 0.0

  IF( ISTRAN(1) > 0 )THEN
    SAL_Global  = 0.0
    SAL1_Global = 0.0
  ENDIF
  
  IF( ISTRAN(2) > 0 )THEN
    TEM_Global   = 0.0
    TEM1_Global  = 0.0
    TEMB_Global  = 0.0
    SHAD_Global  = 0.0
    IF( ISICE > 0 )THEN
      ICETHICK_Global  = 0.0
      ICETEMP_Global   = 0.0
    ENDIF
  ENDIF
  
  IF( ISTRAN(3) > 0 )THEN
    DYE_Global = 0.0
    DYE1_Global = 0.0
  ENDIF
  
  IF( ISTRAN(4) > 0 )THEN
    SFL_Global = 0.0
  ENDIF
  
  IF( ISTRAN(5) > 0 )THEN
    TOX_Global   = 0.0
    TOX1_Global  = 0.0
    TOXB_Global  = 0.0
    TOXB1_Global = 0.0
  ENDIF

  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    HBED_Global     = 0.0
    HBED1_Global    = 0.0
                   
    BEDMAP_Global   = 0
    KBT_Global      = 0
                   
    BDENBED_Global  = 0.0
    PORBED_Global   = 0.0
    VDRBED_Global   = 0.0
    VDRBED1_Global  = 0.0

    IF( (ISTRAN(7) > 0 .AND. ICALC_BL > 0 .AND. NSND > 0 ) .OR. ICALC_BL > 0 )THEN
      QSBDLDX_Global = 0.0
      QSBDLDY_Global = 0.0
    ENDIF

    IF( NSEDFLUME > 0 )THEN
      LAYERACTIVE_Global = 0
      TAU_Global         = 0.0
      D50AVG_Global      = 0.0
      BULKDENS_Global    = 0.0
      TSED_Global        = 0.0
      PERSED_Global      = 0.0
      ETOTO_Global       = 0.0
      DEPO_Global        = 0.0
      IF( ICALC_BL > 0 )THEN
        CBL_Global       = 0.0
        IF( ISTRAN(5) > 0 )THEN
          CBLTOX_Global  = 0.0
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  IF( ISTRAN(6) > 0 )THEN
    SED_Global   = 0.0
    SED1_Global  = 0.0
    SEDB_Global  = 0.0
    SEDB1_Global = 0.0
  ENDIF
  
  IF( ISTRAN(7) > 0 )THEN
    SND_Global   = 0.0
    SND1_Global  = 0.0
    SNDB_Global  = 0.0
    SNDB1_Global = 0.0
  ENDIF

  NLOS_Global = 0
  NLOE_Global = 0
  NLOW_Global = 0
  NLON_Global = 0
  
  CLOS_Global = 0.0
  CLOE_Global = 0.0
  CLOW_Global = 0.0
  CLON_Global = 0.0

  If(MPI_Write_Flag == .TRUE. )THEN
    write(mpi_log_unit, '(a)') 'Array sizes for domain decomposition:'
    write(mpi_log_unit, '(a,I6)') 'NPBEM   =' , NPBEM
    write(mpi_log_unit, '(a,I6)') 'NPBNM   =' , NPBNM
    write(mpi_log_unit, '(a,I6)') 'NPBSM   =' , NPBSM
    write(mpi_log_unit, '(a,I6)') 'NQSIJM  =' , NQSIJM
    write(mpi_log_unit, '(a,I6)') 'NSTVM   =' , NSTVM
    write(mpi_log_unit, '(a,I6)') 'NBBNM   =' , NBBNM
    write(mpi_log_unit, '(a,I6)') 'MLTMSRM =' , MLTMSRM
  End if

  write(mpi_log_unit, '(a)') 'Finished allocating variables for domain decomp'
  Call WriteBreak(mpi_log_unit)

  End subroutine Allocate_Domain_Decomp
