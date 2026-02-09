! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  MODULE WQ_DIAGENESIS
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Module: Module name WQ_DIAGENESIS
  !
  !> @details 
  !
  !> @author 
  !> @date 02/2020
  !---------------------------------------------------------------------------!
  use GLOBAL    
  use INFOMOD
  use JULIANMOD
  
  use Variables_MPI
  use Broadcast_Routines
  use Variables_MPI_Mapping
  use Mod_Map_Global_to_Local
  use Variables_MPI_Write_Out
  
  use Variables_WQ

  implicit none 

  interface
  
  SUBROUTINE SOLVSMBE(SMV1,SMV2,SMA11,SMA22,SMA1,SMA2,SMB11,SMB22)     
    real, intent(IN)  :: SMA11, SMA22, SMA1, SMA2, SMB11, SMB22
    real, intent(OUT) :: SMV1, SMV2      
  END SUBROUTINE
  
  SUBROUTINE SEDFLUXNEW(L, SMSOD1, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                     &
                    SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                    SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                    SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                    RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                    CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, SEDFLUX)
  
    integer, intent(IN)  :: L
    real, intent(IN)     :: SMSOD1, SMCH4S, SMK1CH4, SMO2JC
    real, intent(IN)     :: SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM  
    real, intent(IN)     :: SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM
    real, intent(IN)     :: SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM
    real, intent(INOUT)  :: RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3
    real, intent(INOUT)  :: CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, SEDFLUX, RSMSS    
  END SUBROUTINE SEDFLUXNEW
    
  real FUNCTION ZBRENT(L, ISMERR, SMCH4S, SMK1CH4, SMO2JC,                                 &
                   SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                   SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                   SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                   RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                   CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS) 

    integer, intent(IN)  :: L
    integer, intent(OUT) :: ISMERR
    real, intent(IN)     :: SMCH4S, SMK1CH4, SMO2JC
    real, intent(IN)     :: SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM  
    real, intent(IN)     :: SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM
    real, intent(IN)     :: SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM
    real, intent(INOUT)  :: RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3
    real, intent(INOUT)  :: CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS   
  END FUNCTION 
                                          
  end interface
   
  contains
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SMRIN1
  !
  !> @details  READ THE MAIN CONTROL FILE OF SEDIMENT DIAGENESIS
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------! 
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SMRIN1
  !
  !> @details  READ THE MAIN CONTROL JNP FILE OF SEDIMENT DIAGENESIS
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------! 
  SUBROUTINE SMRIN1_JNP    
  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  use INFOMOD,only:SKIPCOM,READSTR
  
  Character(len = 79), allocatable :: TITLE(:)
  type(fson_value), Pointer :: json_data, item, zones
  
  integer :: IZN,ISMTSB,ISMTSE,ISMTSDT,M,IT,L,IM,I,J,IJC,ISMZX,NAL,LG
      
  Real,parameter :: SMCW2 = 2.739726E-5  !< cm/y to m/day
  real :: SMTHDD,SMTHDP,SMP1NH4,SMP2NH4,SMTHNH4,SMTHNO3,SMP1H2S,SMP2H2S,SMKD1HS,SMKP1HS,SMTHH2S
  real :: SMPOCR,SMKMH2S,SMKCH4,SMTHCH4,SMKSI,SMTHSI,TSMTSB,TSMTSE,SMTSDT,XSMK1H2S
  real :: STEMP,TT20,SUMNN,SUMNP,SUMNC
  
  Real,save,allocatable :: SMKPOC(:) !< Decay rate of POC at 20 degC in Layer 2 for G classes (1/day)
  Real,save,allocatable :: SMKPON(:) !< Decay rate of PON at 20 degC in Layer 2 for G classes (1/day)
  Real,save,allocatable :: SMKPOP(:) !< Decay rate of POP at 20 degC in Layer 2 for G classes (1/day)
  Real,save,allocatable :: SMTHKC(:) !< Constant for temperature adjustment for KPOC (unitless)
  Real,save,allocatable :: SMTHKN(:) !< Constant for temperature adjustment for KPON (unitless)
  Real,save,allocatable :: SMTHKP(:) !< Constant for temperature adjustment for KPOP (unitless)
  Real,allocatable      :: SMPON_TEMP(:)
  Real,allocatable      :: SMPOP_TEMP(:)
  Real,allocatable      :: SMPOC_TEMP(:)
  Real,allocatable      :: SMFNR_TEMP(:)
  Real,allocatable      :: SMFPR_TEMP(:)
  Real,allocatable      :: SMFCR_TEMP(:)
  Real,allocatable      :: SUMCBA(:)       !< Total of the reactive class splits for each algal class - Carbon
  Real,allocatable      :: SUMPBA(:)       !< Total of the reactive class splits for each algal class - Phosphorus
  Real,allocatable      :: SUMNBA(:)       !< Total of the reactive class splits for each algal class - Nitrogen
  
  if( .not. allocated(SMKPOC)  )then
    allocate(SMKPOC(NSMGM))
    allocate(SMKPON(NSMGM))
    allocate(SMKPOP(NSMGM))
    allocate(SMTHKC(NSMGM))
    allocate(SMTHKN(NSMGM))
    allocate(SMTHKP(NSMGM))
    allocate(SMPON_TEMP(NSMGM))
    allocate(SMPOP_TEMP(NSMGM))
    allocate(SMPOC_TEMP(NSMGM))
    allocate(SMFNR_TEMP(NSMGM))
    allocate(SMFPR_TEMP(NSMGM))
    allocate(SMFCR_TEMP(NSMGM))
    allocate(SUMCBA(NALGAEM))
    allocate(SUMPBA(NALGAEM))
    allocate(SUMNBA(NALGAEM))
    
    SMKPOC = 0.0
    SMKPON = 0.0
    SMKPOP = 0.0
    SMTHKC = 0.0
    SMTHKN = 0.0
    SMTHKP = 0.0
    SMPON_TEMP = 0.0
    SMPOP_TEMP = 0.0
    SMPOC_TEMP = 0.0
    SMFNR_TEMP = 0.0
    SMFPR_TEMP = 0.0
    SMFCR_TEMP = 0.0
    SUMCBA = 0.0
    SUMPBA = 0.0
    SUMNBA = 0.0     
  Endif
  
  if( process_id == master_id )then
    json_data => fson_parse("wq_3dsd.jnp")
    Write(*,'(A)') 'WQ: READING WQ_3DSD.JNP - MAIN DIAGENESIS CONTROL FILE'
         
    call fson_get(json_data, "title", TITLE)
    call fson_get(json_data, "initial_condition_option", ISMICI)
    call fson_get(json_data, "number_of_spatial_zones", ISMZ)
    call fson_get(json_data, "write_restart_option", ISMRST)
    call fson_get(json_data, "sediment_temperature_diffusion_coef", SMDIFT)
    call fson_get(json_data, "stoichiometric_coef_for_carbon_diagenesis.nitrification", SMO2NH4)
    call fson_get(json_data, "stoichiometric_coef_for_carbon_diagenesis.denitritrification", SMO2NO3)
    call fson_get(json_data, "stoichiometric_coef_for_carbon_diagenesis.H2S_oxidation", SMO2C)
    
    NSMZ = ISMZ
    SMDIFT = SMDIFT*8.64E4   !< Convert to m^2/day
    
    call fson_get(json_data, "fraction_from_algae.PON", SMFNBA)
    call fson_get(json_data, "fraction_from_algae.POP", SMFPBA)
    call fson_get(json_data, "fraction_from_algae.POC", SMFCBA)
    
    Do NAL = 1,ALG_COUNT
      SUMCBA(NAL) = SMFCBA(NAL,1) + SMFCBA(NAL,2) + SMFCBA(NAL,3)
      SUMPBA(NAL) = SMFPBA(NAL,1) + SMFPBA(NAL,2) + SMFPBA(NAL,3)
      SUMNBA(NAL) = SMFNBA(NAL,1) + SMFNBA(NAL,2) + SMFNBA(NAL,3)
      if( SUMCBA(NAL) < 0.9999 .or. SUMCBA(NAL) > 1.0001) CALL STOPP('ERROR!! SMFCBA(1)+SMFCBA(2)+SMFCBA(3) SHOULD BE 1')
      if( SUMPBA(NAL) < 0.9999 .or. SUMPBA(NAL) > 1.0001) CALL STOPP('ERROR!! SMFPBA(1)+SMFPBA(2)+SMFPBA(3) SHOULD BE 1')
      if( SUMNBA(NAL) < 0.9999 .or. SUMNBA(NAL) > 1.0001) CALL STOPP('ERROR!! SMFNBA(1)+SMFNBA(2)+SMFNBA(3) SHOULD BE 1')
    Enddo
    
    call fson_get(json_data, "particle_mixing.temperature_adjustment_for_particle_mixing_diffusion", SMTHDP)
    call fson_get(json_data, "particle_mixing.reference_concentration_for_GPOC1", SMPOCR)
    call fson_get(json_data, "particle_mixing.half_sat_const_for_oxygen", SMKMDP)
    call fson_get(json_data, "particle_mixing.min_diffusion_coef_for_particle_mixing", XSMDPMIN)
    call fson_get(json_data, "diffusion.temperature_adjustment_for_diffusion", SMTHDD)
    call fson_get(json_data, "diffusion.ratio_of_bioirrigation_to_bioturbation", SMRBIBT)
    SM1OKMDP = 1.0/SMKMDP
    
    call fson_get(json_data, "benthic_stress.activate_hysteresis_in_benthic_mixing", ISMHYST)
    call fson_get(json_data, "benthic_stress.first_order_decay_rate_for_accumulated_benthic_stress", SMKBST)
    call fson_get(json_data, "benthic_stress.critical_overlying_oxygen_conc_for_benthic_hysteresis", SMO2BS)
    call fson_get(json_data, "benthic_stress.duration_the_max_or_min_stress_retained", SMHYLAG)
    call fson_get(json_data, "benthic_stress.critical_hypoxia_duration", SMHYDUR)
    
    call fson_get(json_data, "layer_1.solid_concentration", SMM1)
    call fson_get(json_data, "layer_1.partition_coef_for_particulate_to_dissolved_NH4", SMP1NH4)
    call fson_get(json_data, "layer_1.partition_coef_for_H2S", SMP1H2S)
    call fson_get(json_data, "layer_1.nitrification.half_sat_const_for_ammonium", SMKMNH4)
    call fson_get(json_data, "layer_1.nitrification.half_sat_const_for_DO", SMKMO2N)
    call fson_get(json_data, "layer_1.nitrification.temperature_adjustment_for_KNH4", SMTHNH4)
    call fson_get(json_data, "layer_1.nitrification.temperature_adjustment_for_KNO3", SMTHNO3)
    
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.reaction_velocity_for_dissolved_sulfide_oxidation_at_20_degC", SMKD1HS)
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.reaction_velocity_for_particulate_sulfide_oxidation_at_20_degC", SMKP1HS)
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.temperature_adjustment_for_sulfide_oxidation", SMTHH2S)
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.const_to_normalize_sulfide_oxidation_rate_for_oxygen", SMKMH2S)
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.reaction_velocity_for_methane_at_20_degC", SMKCH4)
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.temperature_adjustment_for_methane_oxidation", SMTHCH4)
    call fson_get(json_data, "layer_1.sulfide_methane_oxidation.critical_salinity_for_producing_CH4_or_H2S", SMCSHSCH)
    
    call fson_get(json_data, "layer_1.PO4_sorption.critical_DO", SMCO2PO4)
    
    call fson_get(json_data, "layer_1.silica_sorption.enhancement_factor", SMDP1SI)
    call fson_get(json_data, "layer_1.silica_sorption.critical_DO", SMCO2SI)
    
    call fson_get(json_data, "layer_2.solid_concentration", SMM2)
    call fson_get(json_data, "layer_2.partition_coef_for_particulate_to_dissolved_NH4", SMP2NH4)
    call fson_get(json_data, "layer_2.partition_coef_for_particulate_to_dissolved_PO4", SMP2PO4)
    call fson_get(json_data, "layer_2.partition_coef_for_H2S", SMP2H2S)
    call fson_get(json_data, "layer_2.decay_rate.PON.decay_rate_at_20_degC", SMKPON)
    call fson_get(json_data, "layer_2.decay_rate.PON.temperature_adjustment", SMTHKN)
    call fson_get(json_data, "layer_2.decay_rate.POP.decay_rate_at_20_degC", SMKPOP)
    call fson_get(json_data, "layer_2.decay_rate.POP.temperature_adjustment", SMTHKP)
    call fson_get(json_data, "layer_2.decay_rate.POC.decay_rate_at_20_degC", SMKPOC)
    call fson_get(json_data, "layer_2.decay_rate.POC.temperature_adjustment", SMTHKC)
          
    call fson_get(json_data, "layer_2.particulate_biogenic_silica.detrital_flux_from_sources_other_than_algae", SMJDSI)
    call fson_get(json_data, "layer_2.particulate_biogenic_silica.dissolution_of_particulate_biogenic_silica.first_order_dissolution_rate_at_20_degC", SMKSI)
    call fson_get(json_data, "layer_2.particulate_biogenic_silica.dissolution_of_particulate_biogenic_silica.temperature_adjustment", SMTHSI)
    call fson_get(json_data, "layer_2.particulate_biogenic_silica.dissolution_of_particulate_biogenic_silica.half_sat_const", SMKMPSI)
    call fson_get(json_data, "layer_2.particulate_biogenic_silica.dissolution_of_particulate_biogenic_silica.silica_saturation_concentration_in_pore_water", SMSISAT)
    call fson_get(json_data, "layer_2.particulate_biogenic_silica.dissolution_of_particulate_biogenic_silica.partition_coefficient_for_silica", SMP2SI)

    call fson_get(json_data, "spatially_constant_initial_conditions.initial_accumulated_benthic_stress", SMBST(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.initial_sediment_temperature", SMT(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.PON", SMPON_TEMP)
    call fson_get(json_data, "spatially_constant_initial_conditions.POP", SMPOP_TEMP)
    call fson_get(json_data, "spatially_constant_initial_conditions.POC", SMPOC_TEMP)
    Do M = 1,3
      SMPON(1,M) = SMPON_TEMP(M)
      SMPOP(1,M) = SMPOP_TEMP(M)
      SMPOC(1,M) = SMPOC_TEMP(M)
    Enddo
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_1.NH4", SM1NH4(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_2.NH4", SM2NH4(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_2.NO3", SM2NO3(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_2.PO4", SM2PO4(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_2.H2S", SM2H2S(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_2.SUU", SMPSI(1))
    call fson_get(json_data, "spatially_constant_initial_conditions.layer_2.SAA", SM2SI(1))
          
    SMFD1NH4 = 1.0 / (1.0 + SMM1*SMP1NH4)
    SMFP1NH4 = 1.0 - SMFD1NH4
    SMFD2NH4 = 1.0 / (1.0 + SMM2*SMP2NH4)
    SMFP2NH4 = 1.0 - SMFD2NH4
    SMKMO2N  = SMKMO2N * 2.0
    SMFD2PO4 = 1.0 / (1.0 + SMM2*SMP2PO4)
    SMFP2PO4 = 1.0 - SMFD2PO4
    
    SMFD1H2S = 1.0 / (1.0 + SMM1*SMP1H2S)
    SMFP1H2S = 1.0 - SMFD1H2S
    SMFD2H2S = 1.0 / (1.0 + SMM2*SMP2H2S)
    SMFP2H2S = 1.0 - SMFD2H2S
    XSMK1H2S = (SMKD1HS*SMKD1HS*SMFD1H2S + SMKP1HS*SMKP1HS*SMFP1H2S) / (2.0*SMKMH2S)
          
    SMFD2SI = 1.0 / (1.0 + SMM2*SMP2SI)
    SMFP2SI = 1.0 - SMFD2SI

    call fson_get(json_data, "spatially_zones", zones)
    Do IZN = 1, fson_value_count(zones)
      ! *** Get the array item 
      item => fson_value_get(zones, IZN)
    
      call fson_get(item, "total_diagenesis_sediment_thickness", SMHSED(IZN))
      call fson_get(item, "sediment_burial_rate", SMW2(IZN))
      call fson_get(item, "diffusion_coef_in_pore_water", SMDD(IZN))
      call fson_get(item, "apparent_diffusion_coef_for_particle_mixing", SMDP(IZN))
      call fson_get(item, "optimal_reaction_velocity_for_nitrification_at_20_degC", SMKNH4(IZN))
      call fson_get(item, "factor_to_enhance_magnitude_of_SOD", SODMULT(IZN))
      call fson_get(item, "fraction_of_water_column_refractory.PON", SMFNR_TEMP)
      call fson_get(item, "fraction_of_water_column_refractory.POP", SMFPR_TEMP)
      call fson_get(item, "fraction_of_water_column_refractory.POC", SMFCR_TEMP)
      
      call fson_get(item, "layer_1.reaction_velocity_for_denitrification_at_20_degC", SMK1NO3(IZN))
      call fson_get(item, "layer_1.factor_to_enhance_sorption_of_PO4_when_DO_less_than_critical", SMDP1PO4(IZN))
      
      call fson_get(item, "layer_2.reaction_velocity_for_denitrification_at_20_degC", SMK2NO3(IZN))
      
      SMW2(IZN)     = SMW2(IZN)*SMCW2                            ! *** Convert to m/day
      SMDTOH(IZN)   = DTWQ/SMHSED(IZN)                           ! *** Assumes fixed DT
      SMHODT(IZN)   = SMHSED(IZN)/DTWQ     !< Fixed DELT only
      SMDP(IZN)     = SMDP(IZN) / (SMHSED(IZN)*SMPOCR+ 1.E-18)
      SMDD(IZN)     = SMDD(IZN) / (SMHSED(IZN)+ 1.E-18)
      SMKNH4(IZN)   = SMKNH4(IZN)*SMKNH4(IZN) * SMKMNH4
      SMK1NO3(IZN)  = SMK1NO3(IZN)*SMK1NO3(IZN)
      SM1DIFT(IZN)  = SMDIFT * SMDTOH(IZN)/(SMHSED(IZN)+ 1.E-18)
      SM2DIFT(IZN)  = 1.0 / (1.0 + SM1DIFT(IZN))
      SMW2DTOH(IZN) = 1.0 + SMW2(IZN)*SMDTOH(IZN)
      SMW2PHODT(IZN)= SMW2(IZN) + SMHODT(IZN)
      SMDPMIN(IZN)  = XSMDPMIN / (SMHSED(IZN)+ 1.E-18)
      Do M = 1,3
        SMFNR(IZN,M) = SMFNR_TEMP(M)
        SMFPR(IZN,M) = SMFPR_TEMP(M)
        SMFCR(IZN,M) = SMFCR_TEMP(M)
      Enddo
      SUMNN = SMFNR(IZN,1) + SMFNR(IZN,2) + SMFNR(IZN,3)
      SUMNP = SMFPR(IZN,1) + SMFPR(IZN,2) + SMFPR(IZN,3)
      SUMNC = SMFCR(IZN,1) + SMFCR(IZN,2) + SMFCR(IZN,3)
      if( SUMNN < 0.9999 .or. SUMNN > 1.0001) &
          call STOPP('ERROR!! SMFNR(I,1) + SMFNR(I,2) + SMFNR(I,3) SHOULD BE 1')
      if( SUMNP < 0.9999 .or. SUMNP > 1.0001) &
          call STOPP('ERROR!! SMFPR(I,1) + SMFPR(I,2) + SMFPR(I,3) SHOULD BE 1') 
      if( SUMNC < 0.9999 .or. SUMNC > 1.0001) &
          call STOPP('ERROR!! SMFCR(I,1) + SMFCR(I,2) + SMFCR(I,3) SHOULD BE 1')
      
    Enddo
        
    call fson_get(json_data, "writing_binary_benthic_flux_rates_output", ISSDBIN)
    call fson_get(json_data, "activate_diagnostic_output", ISMZB)
    call fson_get(json_data, "time_series_output.writing_begin", TSMTSB)
    call fson_get(json_data, "time_series_output.writing_end", TSMTSE)
    call fson_get(json_data, "time_series_output.writing_interval", SMTSDT)
    call fson_get(json_data, "time_series_output.number_of_locations", ISMTS)
  
  endif  ! *** End of master_id block
  
   ! *** Scalar Variables
  call Broadcast_Scalar(ISMZ,     master_id)
  call Broadcast_Scalar(NSMZ,     master_id)
  call Broadcast_Scalar(ISMICI,   master_id)
  call Broadcast_Scalar(ISMRST,   master_id)
  call Broadcast_Scalar(ISMHYST,  master_id)
  call Broadcast_Scalar(ISMZB,    master_id)

  call Broadcast_Scalar(ISMTS,    master_id)
  call Broadcast_Scalar(ISMTSB,   master_id)
  call Broadcast_Scalar(TSMTSB,   master_id)
  call Broadcast_Scalar(ISMTSE,   master_id)
  call Broadcast_Scalar(TSMTSE,   master_id)
  call Broadcast_Scalar(ISMTSDT,  master_id)
  call Broadcast_Scalar(SMTSDT,   master_id)
  call Broadcast_Scalar(ISSDBIN,  master_id)

  call Broadcast_Scalar(SMDIFT,   master_id)

  call Broadcast_Scalar(SMM1,     master_id)
  call Broadcast_Scalar(SMM2,     master_id)
  call Broadcast_Scalar(SMTHDD,   master_id)
  call Broadcast_Scalar(SMTHDP,   master_id)
  call Broadcast_Scalar(SMPOCR,   master_id)
  call Broadcast_Scalar(SMKMDP,   master_id)
  call Broadcast_Scalar(SMKBST,   master_id)
  call Broadcast_Scalar(XSMDPMIN, master_id)
  call Broadcast_Scalar(SMRBIBT,  master_id)

  call Broadcast_Scalar(SMO2BS,   master_id)
  call Broadcast_Scalar(SMHYLAG,  master_id)
  call Broadcast_Scalar(SMHYDUR,  master_id)
  call Broadcast_Scalar(SM1OKMDP, master_id)

  call Broadcast_Scalar(SMP1NH4,  master_id)
  call Broadcast_Scalar(SMP2NH4,  master_id)
  call Broadcast_Scalar(SMKMNH4,  master_id)
  call Broadcast_Scalar(SMKMO2N,  master_id)
  call Broadcast_Scalar(SMTHNH4,  master_id)
  call Broadcast_Scalar(SMTHNO3,  master_id)
  call Broadcast_Scalar(SMP2PO4,  master_id)
  call Broadcast_Scalar(SMCO2PO4, master_id)

  call Broadcast_Scalar(SMFD1NH4, master_id)
  call Broadcast_Scalar(SMFP1NH4, master_id)
  call Broadcast_Scalar(SMFD2NH4, master_id)
  call Broadcast_Scalar(SMFP2NH4, master_id)
  call Broadcast_Scalar(SMKMO2N,  master_id)
  call Broadcast_Scalar(SMFD2PO4, master_id)
  call Broadcast_Scalar(SMFP2PO4, master_id)

  call Broadcast_Scalar(SMP1H2S,  master_id)
  call Broadcast_Scalar(SMP2H2S,  master_id)
  call Broadcast_Scalar(SMKD1HS,  master_id)
  call Broadcast_Scalar(SMKP1HS,  master_id)
  call Broadcast_Scalar(SMTHH2S,  master_id)
  call Broadcast_Scalar(SMKMH2S,  master_id)
  call Broadcast_Scalar(SMKCH4,   master_id)
  call Broadcast_Scalar(SMTHCH4,  master_id)
  call Broadcast_Scalar(SMCSHSCH, master_id)

  call Broadcast_Scalar(SMFD1H2S, master_id)
  call Broadcast_Scalar(SMFP1H2S, master_id)
  call Broadcast_Scalar(SMFD2H2S, master_id)
  call Broadcast_Scalar(SMFP2H2S, master_id)
  call Broadcast_Scalar(XSMK1H2S, master_id)

  call Broadcast_Scalar(SMO2C,    master_id)
  call Broadcast_Scalar(SMO2NO3,  master_id)
  call Broadcast_Scalar(SMO2NH4,  master_id)

  call Broadcast_Scalar(SMKSI,    master_id)
  call Broadcast_Scalar(SMTHSI,   master_id)
  call Broadcast_Scalar(SMKMPSI,  master_id)
  call Broadcast_Scalar(SMSISAT,  master_id)
  call Broadcast_Scalar(SMP2SI,   master_id)
  call Broadcast_Scalar(SMDP1SI,  master_id)
  call Broadcast_Scalar(SMCO2SI,  master_id)
  call Broadcast_Scalar(SMJDSI,   master_id)

  call Broadcast_Scalar(SMFD2SI,  master_id)
  call Broadcast_Scalar(SMFP2SI,  master_id)

  ! *** Array Variables
  call Broadcast_Array(SMFNBA,    master_id)
  call Broadcast_Array(SMFPBA,    master_id)
  call Broadcast_Array(SMFCBA,    master_id)

  call Broadcast_Array(SMKPON,    master_id)
  call Broadcast_Array(SMKPOP,    master_id)
  call Broadcast_Array(SMKPOC,    master_id)

  call Broadcast_Array(SMTHKN,    master_id)
  call Broadcast_Array(SMTHKP,    master_id)
  call Broadcast_Array(SMTHKC,    master_id)

  call Broadcast_Array(SMPON,     master_id)
  call Broadcast_Array(SMPOP,     master_id)
  call Broadcast_Array(SMPOC,     master_id)

  call Broadcast_Array(SM1NH4,    master_id)
  call Broadcast_Array(SM2NH4,    master_id)
  call Broadcast_Array(SM2NO3,    master_id)
  call Broadcast_Array(SM2PO4,    master_id)
  call Broadcast_Array(SM2H2S,    master_id)
  call Broadcast_Array(SMPSI,     master_id)
  call Broadcast_Array(SM2SI,     master_id)
  call Broadcast_Array(SMBST,     master_id)
  call Broadcast_Array(SMT,       master_id)

  call Broadcast_Array(SMHSED,    master_id)
  call Broadcast_Array(SMW2,      master_id)
  call Broadcast_Array(SMDD,      master_id)
  call Broadcast_Array(SMDP,      master_id)
  call Broadcast_Array(SMKNH4,    master_id)
  call Broadcast_Array(SMK1NO3,   master_id)
  call Broadcast_Array(SMK2NO3,   master_id)
  call Broadcast_Array(SMDP1PO4,  master_id)
  call Broadcast_Array(SODMULT,   master_id)

  call Broadcast_Array(SMDTOH,    master_id)
  call Broadcast_Array(SMHODT,    master_id)
  call Broadcast_Array(SM1DIFT,   master_id)
  call Broadcast_Array(SM2DIFT,   master_id)
  call Broadcast_Array(SMW2DTOH,  master_id)
  call Broadcast_Array(SMW2PHODT, master_id)
  call Broadcast_Array(SMDPMIN,   master_id)

  call Broadcast_Array(SMFNR,     master_id)
  call Broadcast_Array(SMFPR,     master_id)
  call Broadcast_Array(SMFCR,     master_id)
  call Broadcast_Array(ISMZMAP,   master_id)

  ! *** SET UP LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY OVER -1OC TO 50OC
  WQTDsMIN = -10
  WQTDsMAX = +50
  STEMP = WQTDsMIN
  WQTDsINC = (WQTDsMAX-WQTDsMIN)/NWQTD

  if( process_id == master_id )then
    write(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for Sediment Diagenesis Temperature Dependency",                                                       &
                                             "WQTDsMIN = ", WQTDsMIN, "WQTDsMAX = ", WQTDsMAX, "WQTDsINC = ", WQTDsINC, "NWQTD = ", NWQTD
    write(2,'(A5,A10,20A10)') "IT", "TEMP", "SMTDND(1)", "SMTDND(2)", "SMTDND(3)", "SMTDPD(1)", "SMTDPD(2)", "SMTDPD(3)", "SMTDCD(1)", "SMTDCD(2)", "SMTDCD(3)",  &
                                            "SMTDDP",    "SMTDDD",    "SMTDNH4",   "SMTDNO3",   "SMK1H2S",   "SMTD1CH4",  "SMTD2CH4",  "SMTDSI"
  endif
  
  do IT = 1,NWQTD
    !STEMP = REAL(IT-1)*0.1 - 4.95
    TT20 = STEMP-20.0
    do M = 1,3
      SMTDND(IT,M) = SMKPON(M) * SMTHKN(M)**TT20
      SMTDPD(IT,M) = SMKPOP(M) * SMTHKP(M)**TT20
      SMTDCD(IT,M) = SMKPOC(M) * SMTHKC(M)**TT20
    enddo
    SMTDDP(IT)   = SMTHDP**TT20
    SMTDDD(IT)   = SMTHDD**TT20
    SMTDNH4(IT)  = SMTHNH4**TT20
    SMTDNO3(IT)  = SMTHNO3**TT20
    SMK1H2S(IT)  = XSMK1H2S * SMTHH2S**TT20
    SMTD1CH4(IT) = 0.97656**TT20 * 20.0
    SMTD2CH4(IT) = SMKCH4 * SMTHCH4**TT20
    SMTDSI(IT)   = SMKSI * SMTHSI**TT20

    if( process_id == master_id )then
      write(2,'(I5,F10.3,20F10.5)') IT, STEMP, SMTDND(IT,1), SMTDND(IT,2), SMTDND(IT,3), SMTDPD(IT,1), SMTDPD(IT,2), SMTDPD(IT,3), SMTDCD(IT,1), SMTDCD(IT,2), SMTDCD(IT,3),  &
                                               SMTDDP(IT)  , SMTDDD(IT)  , SMTDNH4(IT),  SMTDNO3(IT),  SMK1H2S(IT),  SMTD1CH4(IT), SMTD2CH4(IT), SMTDSI(IT)
    endif
    
    STEMP = STEMP + WQTDsINC
  enddo

  ! *** Spatially Variable Handling
  if( ISMICI /= 1 .and. ISMICI /= 2 )then
      
    do L = 2,LA
      do M = 1,NSMG
        SMPON(L,M) = SMPON(1,M)
        SMPOP(L,M) = SMPOP(1,M)
        SMPOC(L,M) = SMPOC(1,M)
      enddo
      SM1NH4(L) = SM1NH4(1)
      SM2NH4(L) = SM2NH4(1)
      SM2NO3(L) = SM2NO3(1)
      SM2PO4(L) = SM2PO4(1)
      SM2H2S(L) = SM2H2S(1)
      SMPSI(L)  = SMPSI(1)
      SM2SI(L)  = SM2SI(1)
      SMBST(L)  = SMBST(1)
      SMT(L)    = SMT(1)
    enddo
  endif

  ! *** Bed diagensis map
  do L = 2,LA
    ISMZMAP(L) = 1
  enddo

  ! *** READ IN MAPPING INFOR. FOR SPATIALLY-VARYING SED parameterS (UNIT #7).
  if( ISMZ  >  1 )then
    allocate(I1D_Global(LCM_Global))
    I1D_Global = 0

    if( process_id == master_id )then
      write(*,'(A)')' WQ: SD READING WQSDMAP.INP'
      open(1,FILE = 'wqsdmap.inp',STATUS = 'UNKNOWN')

      call SKIPCOM(1,'*',2)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      IM = 0
      
      IJC = IC_Global*JC_Global

      do M = 1,IJC
        read(1,*,end = 1111) I, J, ISMZX
        IM = IM + 1

        if( IJCT_Global(I,J) < 1 .or. IJCT_Global(I,J) > 8 .or. ISMZX > ISMZ )then
          PRINT*, 'I, J, IJCT(I,J) = ', I,J,IJCT_Global(I,J)
          call STOPP('ERROR!! INVALID (I,J) IN FILE WQSDMAP.INP')
        endif

        LG = LIJ_Global(I,J)
        I1D_Global(LG) = ISMZX                 ! *** ISMZMAP
      enddo
1111  continue

      if( IM /= (LA_Global-1) )then  
        PRINT *, 'ALL ACTIVE SED. CELLS SHOULD BE MAPPED FOR SED PAR.'
        call STOPP('ERROR!! NUMBER OF LINES IN FILE WQSDMAP.INP  = \ (LA-1)')
      endif
      close(1)
    endif
  
    call Broadcast_Array(I1D_Global, master_id)
    
    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        ISMZMAP(L) = I1D_Global(LG)
      endif
    enddo
    deallocate(I1D_Global)   
    
  endif
  
  END SUBROUTINE SMRIN1_JNP

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSDICI
  !
  !> @details READ IN SPATIALLY VARYING INITIAL CONDITIONS FOR SEDIMENT DIAGENESIS STATE VARIABLES 
  !---------------------------------------------------------------------------!
  ! 2013-03   Paul M. Craig    Restructed to F90, updated variables for OMP
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSDICI 
  
  implicit none

  integer :: I, J, M, NW, L, LG, MM
  
  real    :: XSM1NH4, XSM2NH4, XSM2NO3, XSM2PO4, XSM2H2S, XSMPSI, XSM2SI, XSMBST, XSMT
  
  real,save,allocatable,dimension(:) :: XSMPOC
  real,save,allocatable,dimension(:) :: XSMPON
  real,save,allocatable,dimension(:) :: XSMPOP
  
  character TITLE(3)*79,ICICONT*3

  if( .not. allocated(XSMPOC) )then
    allocate(XSMPOC(NSMGM))
    allocate(XSMPON(NSMGM))
    allocate(XSMPOP(NSMGM))
    XSMPOC = 0.0
    XSMPON = 0.0
    XSMPOP = 0.0
  endif

  write(*,'(A)')' WQ: SD READING WQSDICI.INP'
  open(1,FILE = 'wqsdici.inp',STATUS = 'OLD')

  write(2,999)
  
  call SKIPCOM(1,'*',2)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
    
  do M = 2, LA_Global
    read(1,*) I, J, (XSMPON(NW),NW = 1,NSMG),(XSMPOP(NW),NW = 1,NSMG),(XSMPOC(NW),NW = 1,NSMG),  &
                     XSM1NH4, XSM2NH4, XSM2NO3, XSM2PO4, XSM2H2S, XSMPSI, XSM2SI, XSMBST, XSMT
    if( IJCT_Global(I,J) < 1 .or. IJCT_Global(I,J) > 8 )then
      PRINT*, 'I, J, LINE# = ', I,J,M-1
      call STOPP('ERROR!! INVALID (I,J) IN FILE WQSDICI.INP')
    endif
    ! *** get the global L 
    LG = LIJ_Global(I,J)

    ! *** Set to global spatially dependent variables
    do NW = 1,NSMG
      SMPON_Global(LG,NW) = XSMPON(NW)
      SMPOP_Global(LG,NW) = XSMPOP(NW)
      SMPOC_Global(LG,NW) = XSMPOC(NW)
    enddo
    SM1NH4_Global(LG) = XSM1NH4 
    SM2NH4_Global(LG) = XSM2NH4
    SM2NO3_Global(LG) = XSM2NO3 
    SM2PO4_Global(LG) = XSM2PO4 
    SM2H2S_Global(LG) = XSM2H2S       
    SMPSI_Global(LG)  = XSMPSI 
    SM2SI_Global(LG)  = XSM2SI 
    SMBST_Global(LG)  = XSMBST 
    SMT_Global(LG)    = XSMT
               

    write(2,90) I,J, (SMPON_Global(LG,NW),NW = 1,NSMG), (SMPOP_Global(LG,NW),NW = 1,NSMG), (SMPOC_Global(LG,NW),NW = 1,NSMG),  &
                      SM1NH4_Global(LG), SM2NH4_Global(LG), SM2NO3_Global(LG), SM2PO4_Global(LG), SM2H2S_Global(LG),       &
                      SMPSI_Global(LG),  SM2SI_Global(LG),  SMBST_Global(LG),  SMT_Global(LG)
  enddo
  close(1)
  
    999 FORMAT(1X)
     50 FORMAT(A79)
     52 FORMAT(I7, 1X, A3)
     60 FORMAT(/, A24, I5, A24)
     84 FORMAT(3I5, 20F8.4, F8.2)
     90 FORMAT(2I5, 18E15.4)
  return

  END SUBROUTINE WQSDICI
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SMMBE
  !
  !> @details CONTROL SUBROUTINE (SMMBE) FOR SEDIMENT COMPONENT OF WATER QUALITY MODEL
  !---------------------------------------------------------------------------!
  ! ORGINALLY CODED BY K.-Y. PARK THEN OPTIMIZED AND MODIFIED BY J. M. HAMRICK
  ! 2013-03   Paul M. Craig    Restructed to F90, improved structure and added OMP
  !---------------------------------------------------------------------------!
  SUBROUTINE SMMBE
 
  implicit none

  integer :: ND, LF, LL, L, I, IZ, M, L1, ISMERR, NAL
  integer :: NERR1, NERR2, LIST1(100), LIST2(100)
  integer,save :: NDIAGEN = 0
  
  
  real :: SMBST1, RSMSS, CSODMSM 
  real :: SMCH4S, SMK1CH4, SMO2JC, SMSOD, CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM
  real :: SK1NH4SM, A1NH4SM, A2NH4SM, A22NH4SM, B1NH4SM, B2NH4SM 
  real :: SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM
  real :: SK1H2SSM, A1H2SSM, A2H2SSM, A22H2SSM, B1H2SSM, B2H2SSM
  real :: RSM1PO4,RSM2PO4,A11PO4SM,A22PO4SM,A1PO4SM,A2PO4SM,B11PO4SM,B22PO4SM
  real :: TIMTMP, SMDFNA, SMDFPA, SMDFCA
  real :: SMP1PO4, SMFD1PO4, SMFP1PO4, WQTT, SMP1SI, SMFD1SI, SMFP1SI, A1SISM, A2SISM, A11SISM, SMJ2SI, A22SISM
  real :: B11SISM, B22SISM, RSM1SI, RSM2SI
  real :: SMSOD1(100), SMSOD2(100)
  
  real, external :: ZBRENT

  ! *** FIRST-ORDER DECAY RATE FOR STRESS (/DAY)
  SMBST1 = 1.0 / (1.0 + SMKBST*DTWQ)        ! WQ VARIABLE DT
  NDIAGEN = NDIAGEN + 1
  
  NERR1 = 0
  NERR2 = 0
  LIST1 = 0
  LIST2 = 0 
  
  ! *** First call initializations
  if( NDIAGEN == 1 )then
    RSMSS = 0.
  endif
  
  ! ***   SMHSED = Total active sediment thickness (meters)
  ! ***     SMW2 = Sediment burial rate (m/day)
  ! ***     SMDD = Diffusion coefficient in pore water (m2/day)
  ! ***     SMDP = Apparent diffusion coefficient for particle mixing (m2/day)
  ! ***   SMKNH4 = Optimal reaction velocity for nitrification at 20 degC (m/day)
  ! ***  SMK1NO3 = Reaction velocity for denitrification in layer 1 at 20 degC (m/day)
  ! ***  SMK2NO3 = Reaction velocity for denitrification in layer 2 at 20 degC (m/day)
  ! *** SMDP1PO4 = Factor to enhance sorption of PO4 in layer 1 when DO is greater than DOcPO4 (unitless)
  ! ***  SODMULT = Factor to enhance magnitude of sediment oxygen demand (unitless)
  ! ***   SMDIFT = Sediment temperature diffusion rate (m2/day)
  ! *** XSMDPMIN = Minimum diffusion coefficient for particle mixing (m^2/d)

  ! *** Spatially variable parameters, updated for current DTWQ
  if( ISDYNSTP /= 0 )then
    do I = 1,ISMZ
      SMDTOH(I) = DTWQ/SMHSED(I)                              ! *** day/m
      SMHODT(I) = SMHSED(I)/DTWQ                              ! *** m/day 
      SM1DIFT(I) = SMDIFT * SMDTOH(I)/(SMHSED(I)+ 1.E-18)     ! *** dimensionless
      SM2DIFT(I) = 1.0 / (1.0 + SM1DIFT(I))                   ! *** dimensionless
      SMW2DTOH(I) = 1.0 + SMW2(I)*SMDTOH(I)                   ! *** dimensionless
      SMW2PHODT(I) = SMW2(I) + SMHODT(I)                      ! *** m/day
      SMDPMIN(I) = XSMDPMIN / (SMHSED(I)+ 1.E-18)             ! *** m/day
    enddo
  endif
  
  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND, LF, LL, L, IZ, M, NAL, TIMTMP, SMDFNA, SMDFPA, SMDFCA)  
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = min(LF+LDM-1,LA)

    ! *** SED TEMP., & FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY
    do L = LF,LL
      IZ = ISMZMAP(L)
      SMT(L) = (SMT(L) + SM1DIFT(IZ)*TEM(L,KSZ(L))) * SM2DIFT(IZ)
      ISMT(L) = NINT((SMT(L)-WQTDsMIN)/WQTDsINC)+1  ! *** DSI SINGLE! LINE
      if( ISMT(L) < 1 .or. ISMT(L) > NWQTD )then
        TIMTMP = TIMESEC/86400.
        
        if(process_id == master_id )then
            open(1,FILE = OUTDIR//'ERROR.LOG',POSITION = 'APPEND' ,STATUS = 'UNKNOWN')
            write(1,911) TIMTMP, L, IL(L), JL(L), TEM(L,KSZ(L)), SMT(L)
            close(1)
        endif
        
        PRINT *, 'L, TEM(L,KSZ(L)), SMT(L) = ', Map2Global(L).LG, TEM(L,KSZ(L)), SMT(L)
  
        ! ISMT(L) WAS SET EQUAL TO THE BOUNDS IF IT EXCEEDED THE BOUNDS, THUS
        ! THE MODEL IS NOW ALLOWED TO CONTINUE TO RUN.  THE USER SHOULD CHECK
        ! THE ERROR.LOG FILE FOR SEDIMENT TEMPERATURES OUT OF RANGE.
        !          STOP 'ERROR!! INVALID SEDIMENT TEMPERATURE'
  
        if( ISMT(L)  <  1) ISMT(L) = 1
        if( ISMT(L)  >  NWQTD) ISMT(L) = NWQTD
      endif
    enddo
    
    ! *** Computation of the depositional fluxes of POM
    do M = 1,NSMG
      do L = LF,LL
        SMDFNA = 0.0
        SMDFPA = 0.0
        SMDFCA = 0.0
        ! *** Algal Source Terms
        do NAL = 1,ALG_COUNT
          SMDFNA = SMDFNA + SMFNBA(NAL,M)*ALGAES(NAL).WQANCA*WQDFB(L,NAL)  ! *** RPON
          SMDFPA = SMDFPA + SMFPBA(NAL,M)*WQAPC(L)*WQDFB(L,NAL)            ! *** RPOP
          SMDFCA = SMDFCA + SMFCBA(NAL,M)*WQDFB(L,NAL)                     ! *** RPOC
        enddo                                                             
        ! *** Refractory POM source terms                                 
        SMDFN(L,M) = SMDFNA + SMFNR(ISMZMAP(L),M)*WQDFRN(L)                ! *** RPON
        SMDFP(L,M) = SMDFPA + SMFPR(ISMZMAP(L),M)*WQDFRP(L)                ! *** RPOP
        SMDFC(L,M) = SMDFCA + SMFCR(ISMZMAP(L),M)*WQDFRC(L)                ! *** RPOC
      enddo
    enddo
    
    ! *** Labile POM source terms, only G1 class
    do L = LF,LL
      SMDFN(L,1) = SMDFN(L,1) + WQDFLN(L)  ! *** LPON
      SMDFP(L,1) = SMDFP(L,1) + WQDFLP(L)  ! *** LPOP
      SMDFC(L,1) = SMDFC(L,1) + WQDFLC(L)  ! *** LPOC
    enddo
  enddo  ! *** END OF DOMAIN
  !$OMP END DO
  
  ! *** APPLY OPEN BOUNDARYS
  !OMP SINGLE
  do L1 = 1,NBCSOP
    L = LOBCS(L1)
    SMDFN(L,1) = 0.0
    SMDFP(L,1) = 0.0
    SMDFC(L,1) = 0.0
  enddo
  !OMP END SINGLE

  !$OMP DO PRIVATE(ND, LF, LL, L, M)        
   do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = min(LF+LDM-1,LA)

    !: SMW2 IN M/D,SMW2DTOH(IZ) = 1.0+SMW2*SMDTOH
    ! *** Solving mass balance equations for the concentration of POM 
    ! *** ADD SOURCE TERM FROM ALGAE AND POM DEPOSITION THEN ADJUST FOR BURIAL & DECAY
    do M = 1,NSMG
      do L = LF,LL
        ! ***         Organic M    Algal & POM Sources
        SMPOC(L,M) = (SMPOC(L,M) + SMDFC(L,M)*SMDTOH(ISMZMAP(L))) / (SMW2DTOH(ISMZMAP(L)) + SMTDCD(ISMT(L),M)*DTWQ + 1.E-18)
        SMPON(L,M) = (SMPON(L,M) + SMDFN(L,M)*SMDTOH(ISMZMAP(L))) / (SMW2DTOH(ISMZMAP(L)) + SMTDND(ISMT(L),M)*DTWQ + 1.E-18)
        SMPOP(L,M) = (SMPOP(L,M) + SMDFP(L,M)*SMDTOH(ISMZMAP(L))) / (SMW2DTOH(ISMZMAP(L)) + SMTDPD(ISMT(L),M)*DTWQ + 1.E-18)
      enddo
    enddo
    
    ! *** Computation of diagenesis flux
    do L = LF,LL
      SMDGFN(L) = SMHSED(ISMZMAP(L)) * (SMTDND(ISMT(L),1)*SMPON(L,1)  + SMTDND(ISMT(L),2)*SMPON(L,2))
      SMDGFP(L) = SMHSED(ISMZMAP(L)) * (SMTDPD(ISMT(L),1)*SMPOP(L,1)  + SMTDPD(ISMT(L),2)*SMPOP(L,2))
      SMDGFC(L) = SMHSED(ISMZMAP(L)) * (SMTDCD(ISMT(L),1)*SMPOC(L,1)  + SMTDCD(ISMT(L),2)*SMPOC(L,2))
  
      ! COMMON parameterS: SMBST1 = 1/(1+SMKBST*DTWQ),SM1OKMDP = 1/SMKMDP
      !: use SMTMP(L) TO STORE OLD SMBST(L)
      XSMO20(L) = max( WQV(L,KSZ(L),IDOX), 3.0 )
      SMTMP(L) = SMBST(L)
      ! *** Computation of the accumulated bentic stress
      if( XSMO20(L) < SMKMDP )then
        SMBST(L)  = (SMTMP(L) +DTWQ*(1.0-XSMO20(L)*SM1OKMDP)) * SMBST1
      else
        SMBST(L) = SMBST(L)*SMBST1
      endif
    enddo
  enddo  ! *** END OF DOMAIN
  !$OMP END DO
  
  ! *** APPLY OPEN BOUNDARYS
  !$OMP SINGLE
  do L1 = 1,NBCSOP
    L = LOBCS(L1)
    SMDGFN(L) = 0.0
    SMDGFP(L) = 0.0
    SMDGFC(L) = 0.0
    SMTMP(L) = 0.0
    SMBST(L) = 0.0
  enddo
  !$OMP END SINGLE

  !$OMP DO PRIVATE(ND, LF, LL, L, IZ)             
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = min(LF+LDM-1,LA)

    ! *** BENTHIC MIXING USING HYSTERESIS
    if( ISMHYST == 1 )then
      do L = LF,LL
        if( SCB(L) > 0.5 )then
          if( SMHYST(L) )then
            !IF(XSMO20(L) >= SMO2BS) ISMHYPD(L) = ISMHYPD(L) - 1    ! WQ VAR DT
            if( XSMO20(L) >= SMO2BS) SMHYPD(L) = SMHYPD(L) - DTWQ
            if( SMHYPD(L) <= 0. )then
              SMHYST(L) = .FALSE.
              SMHYPD(L) = 0.
            endif
            SMBST(L) = SMTMP(L)
          else
            !IF(XSMO20(L) < SMO2BS) ISMHYPD(L) = ISMHYPD(L) + 1    ! WQ VAR DT
            if( XSMO20(L) < SMO2BS) SMHYPD(L) = SMHYPD(L) + DTWQ
            if( SMHYPD(L) >= SMHYDUR )then
              SMHYST(L) = .TRUE.
              SMHYPD(L) = SMHYLAG
            endif
          endif
        endif
      enddo
      !ENDDO
    endif

    !: SMDPMIN(IZ) = SMDPMIN/SMHSED
    do L = LF,LL
      IZ = ISMZMAP(L)
      SMW12(L)  = SMDP(IZ)*SMTDDP(ISMT(L)) * SMPOC(L,1) * XSMO20(L) * (1.0-SMKBST*SMBST(L)) / (SMKMDP+XSMO20(L)+ 1.E-18) + SMDPMIN(IZ)
      SMKL12(L) = SMDD(IZ)*SMTDDD(ISMT(L)) + SMRBIBT*SMW12(L)
    enddo
  enddo  ! *** END OF DOMAIN
  !$OMP END DO
  
  ! *** APPLY OPEN BOUNDARYS
  !$OMP SINGLE
  do L1 = 1,NBCSOP
    L = LOBCS(L1)
    SMKL12(L) = 0.0
  enddo
  !$OMP END SINGLE

  !$OMP DO PRIVATE(ND, LF, LL, L, IZ, M, ISMERR, NAL)              &
  !$OMP    PRIVATE(RSMSS, SMCH4S, SMK1CH4, SMO2JC, SMSOD, CSODSM)  &
  !$OMP    PRIVATE(RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM)  &
  !$OMP    PRIVATE(SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM)  &
  !$OMP    PRIVATE(SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM)  &
  !$OMP    PRIVATE(SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM)  &
  !$OMP    PRIVATE(RSM1PO4,RSM2PO4,A11PO4SM,A22PO4SM,A1PO4SM,A2PO4SM,B11PO4SM,B22PO4SM)  &
  !$OMP    PRIVATE(TIMTMP, SMDFNA, SMDFPA, SMDFCA)  &
  !$OMP    PRIVATE(SMP1PO4, SMFD1PO4, SMFP1PO4, WQTT, SMP1SI, SMFD1SI, SMFP1SI, A1SISM, A2SISM, A11SISM, SMJ2SI, A22SISM)  &
  !$OMP    PRIVATE(B11SISM, B22SISM, RSM1SI, RSM2SI)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = min(LF+LDM-1,LA)

    ! *** Initialize conditional variables
    
    ! *** H2S for salinity > SMCSHSCH
    SK1H2SSM = 0.0
    A1H2SSM  = 0.0
    A2H2SSM  = 0.0
    A22H2SSM = 0.0
    B1H2SSM  = 0.0
    B2H2SSM  = 0.0
    
    ! *** Methane for salinity <= SMCSHSCH
    SMCH4S  = 0.0
    SMK1CH4 = 0.0

    ! *** Determine SOD, NH4, NO3 and DOC Fluxes and sediment concentrations
    do L = LF,LL
      if( SCB(L) > 0.5 )then
        IZ = ISMZMAP(L)
        
        ! *** Solving the mass-balance equations
        ! *** Ammonia
        SK1NH4SM = ( SMKNH4(IZ)*SMTDNH4(ISMT(L)) * XSMO20(L) ) / ( (SMKMO2N+XSMO20(L)+ 1.E-12) * (SMKMNH4+SM1NH4(L)) )
        A1NH4SM  = SMKL12(L)*SMFD1NH4 + SMW12(L)*SMFP1NH4 + SMW2(IZ)
        A2NH4SM  = SMKL12(L)*SMFD2NH4 + SMW12(L)*SMFP2NH4
        A22NH4SM = A2NH4SM + SMW2PHODT(IZ)
        B1NH4SM  = WQV(L,KSZ(L),INHX)
        B2NH4SM  = SMDGFN(L) + SMHODT(IZ)*SM2NH4(L)
      
        ! *** Nitrate
        SK1NO3SM = SMK1NO3(IZ)*SMTDNO3(ISMT(L))
        A1NO3SM  = SMKL12(L) + SMW2(IZ)
        A2NO3SM  = SMKL12(L)
        RK2NO3SM = SMK2NO3(IZ)*SMTDNO3(ISMT(L))
        A22NO3SM = A2NO3SM + SMW2PHODT(IZ) + RK2NO3SM
        B1NO3SM  = WQV(L,KSZ(L),INOX)
        B2NO3SM  = SMHODT(IZ)*SM2NO3(L)

        ! *** H2S/CH4
        SMO2JC = SMO2C*SMDGFC(L)
        if( SAL(L,KSZ(L)) > SMCSHSCH )then
          ! *** H2S for salinity > SMCSHSCH
          SK1H2SSM = SMK1H2S(ISMT(L)) * XSMO20(L)
          A1H2SSM  = SMKL12(L)*SMFD1H2S + SMW12(L)*SMFP1H2S + SMW2(IZ)
          A2H2SSM  = SMKL12(L)*SMFD2H2S + SMW12(L)*SMFP2H2S
          A22H2SSM = A2H2SSM + SMW2PHODT(IZ)
          B1H2SSM  = 0.0
          B2H2SSM  = SMHODT(IZ)*SM2H2S(L)
        else
          ! *** Methane for salinity <= SMCSHSCH
          SMCH4S  = (10.0 + HP(L) + SMHSED(IZ)) * SMTD1CH4(ISMT(L)) * SMKL12(L)
          SMK1CH4 = SMTD2CH4(ISMT(L))
        endif

        ! *** BACK SUBSTITUTION TO GET SURFACE MASS TRANSFER COEFFICIENT (SMSS)
        SMSOD = ZBRENT(L, ISMERR, SMCH4S, SMK1CH4, SMO2JC,                                &
                       SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                       SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                       SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                       SM1H2S(L), SM1NH4(L), SM1NO3(L), SM2H2S(L), SM2NH4(L) , SM2NO3(L), &
                       CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS)

        if( ISMERR > 0 )then
          if( ISMERR == 1 .and. NERR1 < 101 )then
            NERR1 = NERR1 + 1
            LIST1(NERR1) = L
            SMSOD1(NERR1) = SMSOD
          elseif( ISMERR == 2 .and. NERR2 < 101 )then
            NERR2 = NERR2 + 1
            LIST2(NERR2) = L
            SMSOD2(NERR2) = SMSOD
          endif
        endif
        SMSS(L)   = RSMSS                    ! *** SURFACE MASS TRANSFER COEFFICIENT (M/DAY)
        
        WQBFO2(L) = -SMSOD * SODMULT(IZ)     ! *** OXYGEN DEMAND (SOD) (G O2/M2/DAY)
        SMCSOD(L) = -CSODSM
        SMNSOD(L) = -RNSODSM
        SMJNIT(L) = RJNITSM
        SMJDEN(L) = RJDENSM

        ! *** Determine COD based on H2S or CH4
        SMJAQH2S(L) = AQJH2SSM + AQJCH4SM
        SMJGCH4(L)  = GJCH4SM
        WQBFNH4(L)  = SMSS(L) * (SMFD1NH4*SM1NH4(L) - WQV(L,KSZ(L),INHX))   ! *** Ammonium flux  (g NH4/m2/day)
        WQBFNO3(L)  = SMSS(L) * (SM1NO3(L) - WQV(L,KSZ(L),INOX))            ! *** Nitrate flux   (g NO3/m2/day)
        WQBFCOD(L)  = SMJAQH2S(L) - SMSS(L)*WQV(L,KSZ(L),ICOD)              ! *** COD flux       (g COD/m2/day)
        
      endif
    enddo

    ! *** Determine PO4 concentrations and water column fluxes using mass diffusion rate from ZBRENT (SMSS)
    do L = LF,LL
      if( SCB(L) > 0.5 )then
        if( XSMO20(L) < SMCO2PO4 )then
          SMP1PO4 = SMP2PO4 * SMDP1PO4(ISMZMAP(L))**(XSMO20(L)/(SMCO2PO4+ 1.E-18))
        else
          SMP1PO4 = SMP2PO4 * SMDP1PO4(ISMZMAP(L))
        endif
        SMFD1PO4 = 1.0 / (1.0 + SMM1*SMP1PO4)
        SMFP1PO4 = 1.0 - SMFD1PO4
        A1PO4SM = SMKL12(L)*SMFD1PO4 + SMW12(L)*SMFP1PO4 + SMW2(ISMZMAP(L))
        A2PO4SM = SMKL12(L)*SMFD2PO4 + SMW12(L)*SMFP2PO4
        A11PO4SM = SMSS(L)*SMFD1PO4 + A1PO4SM
        A22PO4SM = A2PO4SM + SMW2PHODT(ISMZMAP(L))
        B11PO4SM = SMSS(L) * WQPO4D(L,KSZ(L))
        B22PO4SM = SMDGFP(L) + SMHODT(ISMZMAP(L))*SM2PO4(L)

        call SOLVSMBE(RSM1PO4,RSM2PO4,A11PO4SM,A22PO4SM,A1PO4SM,A2PO4SM,B11PO4SM,B22PO4SM)

        SMD1PO4(L) = SMFD1PO4*RSM1PO4
        WQBFPO4D(L) = SMSS(L) * (SMD1PO4(L) - WQPO4D(L,KSZ(L)))
        SM1PO4(L) = RSM1PO4
        SM2PO4(L) = RSM2PO4
      endif
    enddo

    ! *** Determine SI concentrations and water column fluxes using mass diffusion rate from ZBRENT (SMSS)
    if( IWQSI == 1 )then
      do L = LF,LL
        if( SCB(L) > 0.5 )then
          SMDFSI(L) = 0.0
          do NAL = 1, NALGAE
            if( ALGAES(NAL).ISILICA > 0 ) SMDFSI(L) = SMDFSI(L) + (ALGAES(NAL).WQASC*WQDFB(L,NAL) + WQDFSI(L) + SMJDSI) * SMDTOH(ISMZMAP(L))
          enddo
          WQTT = DTWQ * SMTDSI(ISMT(L)) * (SMSISAT-SMFD2SI*SM2SI(L)) / (SMPSI(L)+SMKMPSI+ 1.E-18)
          SMPSI(L) = (SMPSI(L)+SMDFSI(L)) / (SMW2DTOH(ISMZMAP(L))+WQTT+ 1.E-18)
          if( XSMO20(L) < SMCO2SI )then
            SMP1SI = SMP2SI * SMDP1SI**(XSMO20(L)/(SMCO2SI+ 1.E-18))
          else
            SMP1SI = SMP2SI * SMDP1SI
          endif
          SMFD1SI = 1.0 / (1.0 + SMM1*SMP1SI)
          SMFP1SI = 1.0 - SMFD1SI
          A1SISM = SMKL12(L)*SMFD1SI + SMW12(L)*SMFP1SI + SMW2(ISMZMAP(L))
          A2SISM = SMKL12(L)*SMFD2SI + SMW12(L)*SMFP2SI
          A11SISM = SMSS(L)*SMFD1SI + A1SISM
          WQTT = SMTDSI(ISMT(L)) * SMPSI(L) * SMHSED(ISMZMAP(L)) / (SMPSI(L)+SMKMPSI+ 1.E-18)
          SMJ2SI = WQTT * SMSISAT
          A22SISM = A2SISM + SMW2PHODT(ISMZMAP(L)) + WQTT*SMFD2SI
          B11SISM = SMSS(L) * WQSAD(L,KSZ(L))
          B22SISM = SMHODT(ISMZMAP(L))*SM2SI(L) + SMJ2SI

          call SOLVSMBE(RSM1SI,RSM2SI,A11SISM,A22SISM,A1SISM,A2SISM,B11SISM,B22SISM)

          SMD1SI(L) = SMFD1SI*RSM1SI
          WQBFSAD(L) = SMSS(L) * (SMD1SI(L) - WQSAD(L,KSZ(L)))
          SM1SI(L)  = RSM1SI
          SM2SI(L)  = RSM2SI
        endif
      enddo
    endif

  enddo  ! *** END OF DOMAIN
  !$OMP END DO
  !$OMP END PARALLEL

  if(process_id == master_id )then
    if( NERR1 > 0 )then
      open(1,FILE = OUTDIR//'ZBRENT.LOG',STATUS = 'UNKNOWN', POSITION = 'APPEND')

      do ND = 1,NERR1
        L = LIST1(ND)
        if( L > 2 .and. L < LCM )then 
          write(1,401) ITNWQ, MAP2GLOBAL(L).LG, MAP2GLOBAL(L).IG, MAP2GLOBAL(L).JG, SMSOD1(ND), ' ROOT MUST BE BRACKETED FOR ZBRENT (SMMBE)'
        endif
      enddo

      close(1)
    endif
    if( NERR2 > 0 )then
      open(1,FILE = OUTDIR//'ZBRENT.LOG',STATUS = 'UNKNOWN', POSITION = 'APPEND')

      do ND = 1,NERR2
        L = LIST2(ND)
        if( L > 2 .and. L < LCM )then 
          write(1,401) ITNWQ, MAP2GLOBAL(L).LG, MAP2GLOBAL(L).IG, MAP2GLOBAL(L).JG, SMSOD2(ND), ' ZBRENT EXCEEDING MAXIMUM ITERATIONS (SMMBE)'
        endif
      enddo

      close(1)
    endif
  endif
              
  TIMTMP = TIMESEC/TCTMSR

  401 FORMAT(I8,3I5,E12.3,A36)
  911 FORMAT(/,'ERROR: TIME, L, I, J, TEM(L,KSZ(L)), SMT(L) = ', F10.5, 3I4, 2F10.4)

  return

  END SUBROUTINE SMMBE

  END MODULE WQ_DIAGENESIS
  
  
  ! *** EXTERNAL FUNCTIONS AND SUBROUTINES 
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SOLVSMBE
  !
  !> @details SOLVE 2X2 MATRIX 
  !---------------------------------------------------------------------------!
  ! 2013-03   Paul M. Craig    Restructed to F90, updated variables for OMP
  !---------------------------------------------------------------------------! 
  SUBROUTINE SOLVSMBE(SMV1,SMV2,SMA11,SMA22,SMA1,SMA2,SMB11,SMB22)  

  implicit none

  real, intent(IN)  :: SMA11, SMA22, SMA1, SMA2, SMB11, SMB22
  real, intent(OUT) :: SMV1, SMV2
  real :: SMA12, SMA21, SMDET

  SMA12 = -SMA2  
  SMA21 = -SMA1  
  SMDET = SMA11*SMA22 - SMA12*SMA21  
  if( SMDET == 0.0 )then  
    PRINT*, 'SINGULAR MATRIX: A11, A12, A21, A22, B11, B22'  
    PRINT*, SMA11,SMA12,SMA21,SMA22,SMB11,SMB22  
    call STOPP('')  
  endif  
  SMDET = 1.0 / SMDET  
  SMV1 = (SMB11*SMA22 - SMB22*SMA12) * SMDET  
  SMV2 = (SMB22*SMA11 - SMB11*SMA21) * SMDET  

  END SUBROUTINE SOLVSMBE
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SEDFLUXNEW
  !
  !> @details  SOLVE MASS-BALANCE EQ'S FOR NH4, NO3 & H2S/CH4 AND THEIR FLUXES. 
  !---------------------------------------------------------------------------!
  ! 
  !---------------------------------------------------------------------------! 
  SUBROUTINE SEDFLUXNEW(L, SMSOD1, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                       &
                      SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                      SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                      SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                      RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                      CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, SEDFLUX)
  use GLOBAL
  use Variables_WQ
  
  implicit none

  integer, intent(IN)  :: L
  real, intent(IN)     :: SMSOD1, SMCH4S, SMK1CH4, SMO2JC
  real, intent(IN)     :: SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM  
  real, intent(IN)     :: SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM
  real, intent(IN)     :: SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM
  real, intent(INOUT)  :: RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3
  real, intent(INOUT)  :: CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, SEDFLUX, RSMSS

  real :: RRNH4,RRNO3,RRH2S,SMTT1,SMTT2,SMTT3,SMSECH
  real :: A11NH4,B11NH4,B22NH4
  real :: A11NO3,B11NO3,B22NO3
  real :: A11H2S,B11H2S,B22H2S
  real :: CSODMSM, SMJ2H2S, SMSOD
  
  RSMSS = SMSOD1 / (XSMO20(L)+ 1.E-18)

  ! *** NH4
  RRNH4 = SK1NH4SM/(RSMSS+ 1.E-18)
  A11NH4 = RSMSS*SMFD1NH4 + A1NH4SM + RRNH4
  B11NH4 = RSMSS*B1NH4SM
  B22NH4 = B2NH4SM

  call SOLVSMBE(RSM1NH4,RSM2NH4,A11NH4,A22NH4SM,A1NH4SM,A2NH4SM,B11NH4,B22NH4)

  RJNITSM = RRNH4 * RSM1NH4
  
  ! *** NO3
  RRNO3 = SK1NO3SM/(RSMSS+ 1.E-18)
  A11NO3 = RSMSS + A1NO3SM + RRNO3
  B11NO3 = RJNITSM + RSMSS*B1NO3SM
  B22NO3 = B2NO3SM

  call SOLVSMBE(RSM1NO3,RSM2NO3,A11NO3,A22NO3SM,A1NO3SM,A2NO3SM,B11NO3,B22NO3)

  RJDENSM = RRNO3*RSM1NO3 + RK2NO3SM*RSM2NO3

  ! *** H2S/CH4: SMCH4S = 2*SMKL12*SMCH4S

  ! *** SMO2JC          - Max CH4 production
  ! *** SMO2NO3*RJDENSM - Nitrification
  ! *** SMJ2H2S         - Flux of methane in oxygen equivalent units (M/L2/T)
  SMJ2H2S = max(SMO2JC - SMO2NO3*RJDENSM, 0.0)
  if( SAL(L,KSZ(L)) > SMCSHSCH )then
    RRH2S = SK1H2SSM/(RSMSS+ 1.E-18)
    SMTT1 = RSMSS*SMFD1H2S
    A11H2S = SMTT1 + A1H2SSM + RRH2S
    B11H2S = B1H2SSM
    B22H2S = B2H2SSM + SMJ2H2S

    call SOLVSMBE(RSM1H2S,RSM2H2S,A11H2S,A22H2SSM,A1H2SSM,A2H2SSM,B11H2S,B22H2S)

    AQJH2SSM = SMTT1*RSM1H2S
    CSODSM = RRH2S*RSM1H2S
    AQJCH4SM = 0.0
    GJCH4SM = 0.0
  else
    CSODMSM = min( SQRT(SMCH4S*SMJ2H2S), SMJ2H2S )
    SMTT2 = SMK1CH4 / (RSMSS+ 1.E-18)
    if( SMTT2 < 80.0 )then
      SMTT3 = EXP(SMTT2)
      SMSECH = 2.0 / (SMTT3 + 1.0/SMTT3)
    else
      SMSECH = 0.0
    endif
    AQJCH4SM = CSODMSM*SMSECH
    CSODSM = CSODMSM - AQJCH4SM
    GJCH4SM = SMJ2H2S - CSODMSM
    AQJH2SSM = 0.0
  endif
  RNSODSM = SMO2NH4*RJNITSM
  SMSOD = CSODSM + RNSODSM
  SEDFLUX = SMSOD - SMSOD1

  return
  END SUBROUTINE SEDFLUXNEW
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Diagen_Log
  !
  !> @details  Write a report of diagenesis inputs to the WQ Log
  !---------------------------------------------------------------------------!
  SUBROUTINE DIAGEN_LOG
  
  
  
  END SUBROUTINE DIAGEN_LOG
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: FUNCTION ZBRENT
  !
  !> @details  USING BRENT'S METHOD (ZBRENT), FIND THE ROOT OF A FUNC 
  !            SEDFLUX KNOWN TO LIE BETWEEN RMIN & RMAX WITHIN AN ACCURACY 
  !            OF TOL (P. 253 IN NUMERICAL RECIPE). 
  !---------------------------------------------------------------------------!
  ! 2013-03           Paul M. Craig    Restructed to F90, updated variables for OMP
  !---------------------------------------------------------------------------! 
  FUNCTION ZBRENT(L, ISMERR, SMCH4S, SMK1CH4, SMO2JC,                                   &
                     SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                     SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                     SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                     RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                     CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS)
  
  implicit none
  
  integer, intent(IN)  :: L
  integer, intent(OUT) :: ISMERR
  real, intent(IN)     :: SMCH4S, SMK1CH4, SMO2JC
  real, intent(IN)     :: SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM  
  real, intent(IN)     :: SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM
  real, intent(IN)     :: SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM
  real, intent(INOUT)  :: RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3
  real, intent(INOUT)  :: CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS
  
  integer :: IZMAX,II
  real    :: EPS,TOL,RMIN,RMAX,ZBRENT
  real    :: A,B,C,D,E,FA,FB,FC,TOL1,XM
  real    :: P,Q,R,S

  external :: SEDFLUXNEW

  parameter (IZMAX = 100,EPS = 3.0E-8,TOL = 1.0E-5,RMIN = 1.0E-4,RMAX = 100.0)

  ISMERR = 0
  
  ! *** Bracket solution using RMIN and RMAX
  A = RMIN
  B = RMAX

  !FA = SEDFLUX(A,      RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3, RSMSS)
  call SEDFLUXNEW(L, A, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                              &
                  SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                  SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                  SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                  RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                  CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, FA)

  !FB = SEDFLUX(B,      RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3, RSMSS)
  call SEDFLUXNEW(L, B, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                              &
                  SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                  SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                  SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                  RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                  CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, FB)


  ZBRENT = 0
  if( FA*FB > 0.0 )then
    ISMERR = 1
    return
  endif
  
  ! *** Iterate over solution to converge on concentrations and fluxes
  FC = FB
  do II = 1,IZMAX
    if( FB*FC > 0.0 )then
      C = A
      FC = FA
      D = B-A
      E = D
    endif
    if( ABS(FC) < ABS(FB) )then
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
    endif
    TOL1 = 2.0*EPS*ABS(B) + 0.5*TOL
    XM = 0.5 * (C-B)
    
    ! *** Check for convergence
    if( ABS(XM) <= TOL1 .or. FB == 0.0 )then
      ZBRENT = B          ! *** Successful solution
      return
    endif
    
    if( ABS(E) >= TOL1 .and. ABS(FA) > ABS(FB) )then
      S = FB / FA
      if( A == C )then
        P = 2.0 * XM * S
        Q = 1.0 - S
      else
        Q = FA / FC
        R = FB / FC
        P = S * (2.0*XM*Q*(Q-R) - (B-A)*(R-1.0))
        Q = (Q-1.0) * (R-1.0) * (S-1.0)
      endif
      if( P > 0.0) Q = -Q
      P = ABS(P)
      if( 2.0*P  <  min(3.0*XM*Q-ABS(TOL1*Q), ABS(E*Q)) )then
        E = D
        D = P / Q
      else
        D = XM
        E = D
      endif
    else
      D = XM
      E = D
    endif
    A = B
    FA = FB
    if( ABS(D) > TOL1 )then
      B = B + D
    else
      B = B + SIGN(TOL1,XM)
    endif
    
    ! *** Calculate fluxes
    !FB = SEDFLUX(B, RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3, RSMSS)
    call SEDFLUXNEW(L, B, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                              &
                    SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                    SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                    SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                    RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                    CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, FB)

  enddo
  
  ! *** Solution did not converge.  Exit with an error flag
  ISMERR = 2
  ZBRENT = B

  return
  END

