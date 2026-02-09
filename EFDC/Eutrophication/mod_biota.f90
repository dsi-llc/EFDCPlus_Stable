! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE WQ_BIOTA
!--------------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!--------------------------------------------------------------------------------!
! Module: Module name WQ_BIOTA
!
!> @details Part of the eutrophication module and handles the following:
!              - Main biota input
!              - Vertical migration velocities of cyanobacteria cells/colonies
!              - Optional feedback of macrophyte growth and density on the 
!                hydrodynamics.  
!
!> @author Paul M. Craig, Duc Kien Tran
!> @date 06/2021, updated 11/2022
!--------------------------------------------------------------------------------!
  
  use GLOBAL    
  use Variables_WQ
  use Variables_MPI
  use WQ_ZOOPLANKTON
  
  implicit none

  real    :: BETAMAC_P = 1.0, BETAMAC_D = 0.39, CE4MAC = 1.05 ! *** James & O'Donncha (2019)
  
contains

!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine ALGAECONTROL
!
!> @details  READ THE ALGAE CONTROL JNP FILE FOR WATER QUALITY
!---------------------------------------------------------------------------!
!    I/O CONTROL VARIABLES
!    SPATIALLY AND TEMPORALLY CONSTANT parameterS
!---------------------------------------------------------------------------! 
SUBROUTINE ALGAECONTROL
  
  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  
  implicit none
  
  integer :: IZ, M, ITMP, NS, NAL, L, K, ISMOB
  real :: WTEMP, VISCW, BDLTMP
  type(fson_value), pointer :: json_data, phytogroups, item
  
  if( process_id == master_id  .and. NALGAE > 0 )then
    json_data => fson_parse("wq_biota.jnp")
    
    ! *** Get the phytoplankton group's parameters as an array
    call fson_get(json_data, "groups", phytogroups)
    
    ! *** Loop through each array item
    do NAL = 1, fson_value_count(phytogroups)
      ! *** Get the array item 
      item => fson_value_get(phytogroups, NAL)
      
      ! *** Lookup the values from the array
      call fson_get(item, "index",                                          ALGAES(NAL).IDN)            ! *** 
      call fson_get(item, "mobility_flag",                                  ISMOB)                      ! *** 
      call fson_get(item, "salinity_toxicity_flag",                         ALGAES(NAL).ISTOX)          ! *** 
      call fson_get(item, "winter_bloom_flag",                              ALGAES(NAL).ISBLOOM)        ! *** 
      call fson_get(item, "silica_active_flag",                             ALGAES(NAL).ISILICA)        ! *** 
      call fson_get(item, "settling_velocity",                              ALGAES(NAL).WQWS(1))        ! *** 
        
      if( ISMOB == 0 )then
        ALGAES(NAL).ISMOBILE = .FALSE.     
      else
        ALGAES(NAL).ISMOBILE = .TRUE.
      endif
      
      ! *** General parameters
      call fson_get(item, "carbon_to_chla_ratio",                           ALGAES(NAL).WQCHLA)         ! *** 
      call fson_get(item, "light_extinction_due_to_shading",                ALGAES(NAL).WQKEMAC)        ! *** Light_extinction_due_to_shading
      call fson_get(item, "O2_to_C_ratio",                                  ALGAES(NAL).WQALGOCR)       ! *** 
      call fson_get(item, "O2_to_N_ratio",                                  ALGAES(NAL).WQALGONT)       ! *** 
      call fson_get(item, "stoichiometry.nitrogen_to_carbon_ratio",         ALGAES(NAL).WQANCA)         ! *** 
      call fson_get(item, "stoichiometry.silica_to_carbon_ratio",           ALGAES(NAL).WQASC)          ! *** 
      call fson_get(item, "stoichiometry.factor_to_modify_C_to_P_ratio",    ALGAES(NAL).WQAPCM)         ! *** 
      
      ! *** Growth parameters
      call fson_get(item, "growth.max_growth_rate",                         ALGAES(NAL).WQPMA(1))       ! *** 
      call fson_get(item, "growth.photosynthesis_O2_to_C_ratio",            ALGAES(NAL).WQAOCRP)        ! *** 
      call fson_get(item, "growth.phosphorus_half_saturation",              ALGAES(NAL).WQKHPA)         ! *** 
      call fson_get(item, "growth.nitrogen_half_saturation",                ALGAES(NAL).WQKHNA)         ! *** 
      call fson_get(item, "growth.silica_half_saturation",                  ALGAES(NAL).WQKHS)          ! *** 
      call fson_get(item, "growth.CO2_half_saturation",                     ALGAES(NAL).WQKHCO2)        ! *** 
      call fson_get(item, "growth.halved_microsystis_growth_salinity",      ALGAES(NAL).WQSTOX)         ! *** 
      call fson_get(item, "growth.optimal_depth",                           ALGAES(NAL).WQDOP(1))       ! *** 
      call fson_get(item, "growth.optimal_growth.lower_temperature",        ALGAES(NAL).WQTM1)          ! *** 
      call fson_get(item, "growth.optimal_growth.upper_temperature",        ALGAES(NAL).WQTM2)          ! *** 
      call fson_get(item, "growth.optimal_growth.lower_coefficient",        ALGAES(NAL).WQKG1)          ! *** 
      call fson_get(item, "growth.optimal_growth.upper_coefficient",        ALGAES(NAL).WQKG2)          ! *** 
      call fson_get(item, "growth.plant_density_half_saturation_factor",    ALGAES(NAL).WQKBP(1))       ! *** 
      call fson_get(item, "growth.minimum_biomass",                         ALGAES(NAL).WQBMIN(1))      ! *** 
      call fson_get(item, "growth.maximum_biomass",                         ALGAES(NAL).WQBMAX)         ! *** 

      ! *** Metabolism parameters
      call fson_get(item, "metabolism.reference_temperature",               ALGAES(NAL).WQTR)           ! *** 
      call fson_get(item, "metabolism.reference_rate",                      ALGAES(NAL).WQBMRA(1))      ! *** 
      call fson_get(item, "metabolism.effect_of_temperature",               ALGAES(NAL).WQKTB)          ! *** 
      call fson_get(item, "metabolism.respiration_O2_to_C_ratio",           ALGAES(NAL).WQAOCRR)        ! *** 
      call fson_get(item, "metabolism.half_sat_DOC_excretion",              ALGAES(NAL).WQKHRA(1))      ! *** 
      call fson_get(item, "metabolism.fraction_of_carbon.DOC",              ALGAES(NAL).WQFCDB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_phosphorus.RPOP",         ALGAES(NAL).WQFPRB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_phosphorus.LPOP",         ALGAES(NAL).WQFPLB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_phosphorus.DOP",          ALGAES(NAL).WQFPDB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_phosphorus.PO4",          ALGAES(NAL).WQFPIB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_nitrogen.RPON",           ALGAES(NAL).WQFNRB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_nitrogen.LPON",           ALGAES(NAL).WQFNLB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_nitrogen.DON",            ALGAES(NAL).WQFNDB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_nitrogen.NH4",            ALGAES(NAL).WQFNIB)         ! *** 
      call fson_get(item, "metabolism.fraction_of_silica.SU",               ALGAES(NAL).WQFSPD)         ! *** 
      call fson_get(item, "metabolism.fraction_of_silica.SA",               ALGAES(NAL).WQFSID)         ! *** 
      
      ! *** Predation parameters
      call fson_get(item, "predation.predation_rate", ALGAES(NAL).WQPRRA(1))
      ! *** Set algal death rate if zooplankton is simulated
      if( IWQZPL == 1  ) ALGAES(NAL).WQDRA(1) = ALGAES(NAL).WQPRRA(1)
        
      call fson_get(item, "predation.optimal_predation.lower_temperature",  ALGAES(NAL).WQTP1)          ! *** 
      call fson_get(item, "predation.optimal_predation.upper_temperature",  ALGAES(NAL).WQTP2)          ! *** 
      call fson_get(item, "predation.optimal_predation.lower_coefficient",  ALGAES(NAL).WQKP1)          ! *** 
      call fson_get(item, "predation.optimal_predation.upper_coefficient",  ALGAES(NAL).WQKP2)          ! *** 
      call fson_get(item, "predation.fraction_of_carbon.RPOC",              ALGAES(NAL).WQFCRP)         ! *** 
      call fson_get(item, "predation.fraction_of_carbon.LPOC",              ALGAES(NAL).WQFCLP)         ! *** 
      call fson_get(item, "predation.fraction_of_carbon.DOC",               ALGAES(NAL).WQFCDP)         ! *** 
      call fson_get(item, "predation.fraction_of_phosphorus.RPOP",          ALGAES(NAL).WQFPRP)         ! *** 
      call fson_get(item, "predation.fraction_of_phosphorus.LPOP",          ALGAES(NAL).WQFPLP)         ! *** 
      call fson_get(item, "predation.fraction_of_phosphorus.DOP",           ALGAES(NAL).WQFPDP)         ! *** 
      call fson_get(item, "predation.fraction_of_phosphorus.PO4",           ALGAES(NAL).WQFPIP)         ! *** 
      call fson_get(item, "predation.fraction_of_nitrogen.RPON",            ALGAES(NAL).WQFNRP)         ! *** 
      call fson_get(item, "predation.fraction_of_nitrogen.LPON",            ALGAES(NAL).WQFNLP)         ! *** 
      call fson_get(item, "predation.fraction_of_nitrogen.DON",             ALGAES(NAL).WQFNDP)         ! *** 
      call fson_get(item, "predation.fraction_of_nitrogen.NH4",             ALGAES(NAL).WQFNIP)         ! *** 
      call fson_get(item, "predation.fraction_of_silica.SU",                ALGAES(NAL).WQFSPP)         ! *** 
      call fson_get(item, "predation.fraction_of_silica.SA",                ALGAES(NAL).WQFSID)         ! *** 
           
      call fson_get(item, "growth_velocity_limitation_option",              ALGAES(NAL).IWQVLIM)      
			CALL fson_get(item, "vel_lim.half_saturation_velocity",               ALGAES(NAL).WQKMV(1))       ! *** 
      call fson_get(item, "vel_lim.mininum_velocity",                       ALGAES(NAL).WQKMVMIN(1))    ! *** 
			CALL fson_get(item, "vel_lim.lf_param_a",                             ALGAES(NAL).WQKMVA(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_b",                             ALGAES(NAL).WQKMVB(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_c",                             ALGAES(NAL).WQKMVC(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_d",                             ALGAES(NAL).WQKMVD(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_e",                             ALGAES(NAL).WQKMVE(1))      ! *** 

      if( .not. ALGAES(NAL).ISMOBILE )then
        call fson_get(item, "layerthreshold",                               ALGAES(NAL).THRESHOLD)      ! *** Biomass concentration limit to force macrophyte mass into the next layer (mg C/L/m)
        call fson_get(item, "base_depth",                                   ALGAES(NAL).BASEDEPTH)      ! *** Distance above the bottom (i.e. depth) of the "base" of the macrophyte/periphyton growth (m)
        call fson_get(item, "maximum_length",                               ALGAES(NAL).MAXLENGTH)      ! *** Maximum length from the "base" to allow macrophyte/periphyton growth (m)
        call fson_get(item, "use_hydro_feedback",                           ALGAES(NAL).ISDRAG)      
        call fson_get(item, "laminar_flow_option",                          ALGAES(NAL).ISMACL)      

        call fson_get(item, "hydro_feedback.dragcoefficient",               ALGAES(NAL).DRAGCOEFF)      ! *** Vegetagive drag coefficient.  Includes stems and leaves
        call fson_get(item, "hydro_feedback.stemdiameter",                  ALGAES(NAL).MINDIAMETER)    ! *** Minimum stem diameter. As growth occurs, stemdiameter is updated based on % dry matter, stemdensity and concentration (m)
        call fson_get(item, "hydro_feedback.stemdensity",                   ALGAES(NAL).STEMDENSITY)    ! *** As growth occurs, stemdensity is constant (# stems/m2)
        call fson_get(item, "hydro_feedback.blocking_factor",               ALGAES(NAL).ALPMAC)         ! *** Macrophyte cell blocking factor
      endif      

      call fson_get(item, "variable_settling",                              ALGAES(NAL).ISVARSETTLE)    ! *** Option for variable algae settling.  
                                                                                                        ! ***    0 - use user specified settling/floating velocity, 
                                                                                                        ! ***    1 - Velocity varyies with time of day, 
                                                                                                        ! ***    2 - Velocity varies with time of day and available light, and 
                                                                                                        ! ***    3 - Dynamic Velocity (Visser, et.al. 1997)
      if( .not. ALGAES(NAL).ISMOBILE ) ALGAES(NAL).ISVARSETTLE = 0
      if( ALGAES(NAL).ISVARSETTLE > 0 )then 
        call fson_get(item, "amplitude",                                    ALGAES(NAL).AMPLITUDE)      ! *** Amplitude (m) of time varying settling                      (options 1 and 2)
        call fson_get(item, "phase_shift",                                  ALGAES(NAL).PHASESHIFT)     ! *** Phase shift (rad) of COS function of time varying settling  (options 1 and 2)
        call fson_get(item, "cell_density_min",                             ALGAES(NAL).DENSITYMIN)     ! *** Minimum cell density (kg/m3)                           (Visser)
        call fson_get(item, "cell_density_max",                             ALGAES(NAL).DENSITYMAX)     ! *** Maximum cell density (kg/m3)                           (Visser)
        call fson_get(item, "cell_density_init",                            ALGAES(NAL).DENSITYINI)     ! *** Initial cell density (kg/m3)                           (Visser)
        call fson_get(item, "cell_to_colony_ratio",                         ALGAES(NAL).CELL2COL)       ! *** Ratio of cell volume to colony volume                  (Visser)
        call fson_get(item, "cell_drag_coefficient",                        ALGAES(NAL).CELLDRAG)       ! *** Drag coefficient of cell/colony for Stokes             (Visser)
        call fson_get(item, "cell_radius" ,                                 ALGAES(NAL).CELLRAD)        ! *** Cell radius for Stokes  (m)                            (Visser)
        call fson_get(item, "cell_light_min",                               ALGAES(NAL).LIGHTMIN)       ! *** Minimum light for algal cell density increase (W/m2)   (Visser)
        call fson_get(item, "cell_light_sat",                               ALGAES(NAL).LIGHTSAT)       ! *** Light saturation constant for phytoplankton   (W/m2)   (Visser)
        call fson_get(item, "cell_increase_c1",                             ALGAES(NAL).C1)             ! *** Density increase multiplier   (s2/m3)                  (Visser)
        call fson_get(item, "cell_increase_c2",                             ALGAES(NAL).C2)             ! *** Density increase offset       (kg/m3/s2)               (Visser)
        call fson_get(item, "cell_decrease_f1",                             ALGAES(NAL).F1)             ! *** Density decrease slope        (/s)                     (Visser)
        call fson_get(item, "cell_decrease_f2",                             ALGAES(NAL).F2)             ! *** Density decrease intercept    (kg/m3/s)                (Visser)
        call fson_get(item, "minimum_depth",                                ALGAES(NAL).HVM)            ! *** Minimum water depth for active vertical migration (m) 
      endif

      ! *** Apply QC and special conditions
    enddo   
    
  Endif    ! *** End of master_id block
    
  ! *** Broadcast water quality global settings, not dependent on grid.
  
  ! *** Algae's parameters
  Do NAL = 1, NALGAE
     call Broadcast_Scalar(ALGAES(NAL).IDN,         master_id)
     call Broadcast_Scalar(ALGAES(NAL).ISMOBILE,    master_id)
     call Broadcast_Scalar(ALGAES(NAL).ISTOX,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).ISBLOOM,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).ISILICA,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQWS(1),     master_id)
          
     call Broadcast_Scalar(ALGAES(NAL).WQCHLA,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKEMAC,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQALGOCR,    master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQALGONT,    master_id)

     call Broadcast_Scalar(ALGAES(NAL).ISVARSETTLE, master_id)
     call Broadcast_Scalar(ALGAES(NAL).AMPLITUDE,   master_id)
     call Broadcast_Scalar(ALGAES(NAL).PHASESHIFT,  master_id)
     call Broadcast_Scalar(ALGAES(NAL).DENSITYMIN,  master_id)
     call Broadcast_Scalar(ALGAES(NAL).DENSITYMAX,  master_id)
     call Broadcast_Scalar(ALGAES(NAL).DENSITYINI,  master_id)
     call Broadcast_Scalar(ALGAES(NAL).CELL2COL,    master_id)
     call Broadcast_Scalar(ALGAES(NAL).CELLDRAG,    master_id)
     call Broadcast_Scalar(ALGAES(NAL).CELLRAD,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).LIGHTMIN,    master_id)
     call Broadcast_Scalar(ALGAES(NAL).LIGHTSAT,    master_id)
     call Broadcast_Scalar(ALGAES(NAL).C1,          master_id)
     call Broadcast_Scalar(ALGAES(NAL).C2,          master_id)
     call Broadcast_Scalar(ALGAES(NAL).F1,          master_id)
     call Broadcast_Scalar(ALGAES(NAL).F2,          master_id)
     call Broadcast_Scalar(ALGAES(NAL).HVM,         master_id)

     call Broadcast_Scalar(ALGAES(NAL).WQANCA,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQASC,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQAPCM,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQPMA(1),    master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQAOCRP,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKHPA,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKHNA,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKHS,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKHCO2,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQSTOX,      master_id)  
     call Broadcast_Scalar(ALGAES(NAL).WQDOP(1),    master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQTM1,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQTM2,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKG1,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKG2,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKBP(1),    master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQBMIN(1),   master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQBMAX,      master_id) 
     call Broadcast_Scalar(ALGAES(NAL).WQTR,        master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQBMRA(1),   master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKTB,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQAOCRR,     master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKHRA(1),   master_id) 
     call Broadcast_Scalar(ALGAES(NAL).WQFCDB,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFPRB,      master_id)     
     call Broadcast_Scalar(ALGAES(NAL).WQFPLB,      master_id)         
     call Broadcast_Scalar(ALGAES(NAL).WQFPDB,      master_id)          
     call Broadcast_Scalar(ALGAES(NAL).WQFPIB,      master_id)        
     call Broadcast_Scalar(ALGAES(NAL).WQFNRB,      master_id)   
     call Broadcast_Scalar(ALGAES(NAL).WQFNLB,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFNDB,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFNIB,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFSPD,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFSID,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQPRRA(1),   master_id)
     if( IWQZPL == 1 ) Call Broadcast_Scalar(ALGAES(NAL).WQDRA(1),   master_id)
     
     call Broadcast_Scalar(ALGAES(NAL).WQTP1,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQTP2,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKP1,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQKP2,       master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFCRP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFCLP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFCDP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFPRP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFPLP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFPDP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFPIP,      master_id)

     call Broadcast_Scalar(ALGAES(NAL).WQFNRP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFNLP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFNDP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFNIP,      master_id)                                            
     call Broadcast_Scalar(ALGAES(NAL).WQFSPP,      master_id)
     call Broadcast_Scalar(ALGAES(NAL).WQFSIP,      master_id)  
     
     call Broadcast_Scalar(ALGAES(NAL).IWQVLIM,     master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).WQKMV(1),    master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQKMVMIN(1), master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQKMVA(1),   master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQKMVB(1),   master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQKMVC(1),   master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQKMVD(1),   master_id)                                              
     call Broadcast_Scalar(ALGAES(NAL).WQKMVE(1),   master_id)                                              

     ! *** Macrophyte/Periphyton settings
     call Broadcast_Scalar(ALGAES(NAL).THRESHOLD,   master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).BASEDEPTH,   master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).MAXLENGTH,   master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).ISDRAG,      master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).ISMACL,      master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).DRAGCOEFF,   master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).MINDIAMETER, master_id)                                               
     call Broadcast_Scalar(ALGAES(NAL).STEMDENSITY, master_id)                   
     call Broadcast_Scalar(ALGAES(NAL).ALPMAC,      master_id)
  enddo
          
  ! *** macrophyte hydrodynamics feedback setting
  Do NAL = 1, NALGAE
    If( .not. ALGAES(NAL).ISMOBILE )then
      Do L = 2,LA
        DIAMETER_MAC(L,NAL) = ALGAES(NAL).MINDIAMETER
        BDLTMP = DIAMETER_MAC(L,NAL)*DIAMETER_MAC(L,NAL)*ALGAES(NAL).STEMDENSITY  !Normalized area of MACROALGAE in this cell (-)   
        MVEGZ(L,NAL) = max(1.E-18,(1. - ALGAES(NAL).ALPMAC*BDLTMP))               !Blockage ratio in Z (-)
        MVEGZ(1,NAL) = 1.
        MVEGZ(LC,NAL)= 1.
        MACAD(L,NAL) = ALGAES(NAL).STEMDENSITY*DIAMETER_MAC(L,NAL)  !a*d (m) used to calculate drag, see papers written on this James and O'Donncha EFM 2019 BDLPSQ
      enddo      
    endif
  enddo
  
  Do NAL = 1, NALGAE
    ALGAES(NAL).WQCHLA = 1.0/(ALGAES(NAL).WQCHLA + 1.E-12)
    ALGAES(NAL).WQSTOX = ALGAES(NAL).WQSTOX*ALGAES(NAL).WQSTOX
    ! *** Zeroing transport bypass flag of macro algae
    if( .not. ALGAES(NAL).ISMOBILE )then 
      ISTRWQ(19+NAL) = 0
    endif   
    
    call Broadcast_Scalar(ALGAES(NAL).WQCHLA,   master_id)
    call Broadcast_Scalar(ALGAES(NAL).WQSTOX,   master_id)
  enddo
  
  call Broadcast_Array(ISTRWQ,   master_id)
  
  WQAOCR = ALGAES(1).WQALGOCR
  WQAONT = ALGAES(1).WQALGONT
  call Broadcast_Scalar(WQAOCR,   master_id)
  call Broadcast_Scalar(WQAONT,   master_id)
  
  ! *** Deactivate silica flag if silica deactivated
  if( IWQSI == 0 )then
    do NAL = 1, NALGAE
      ALGAES(NAL).ISILICA = 0
    enddo
  endif
  
  !IF(  IWQBEN > 0  )then         ! *** Deprecated - Multiple alagal classes can interact with silica and no need to duplicate variables
  !  Do NAL = 1, NALGAE
  !    SMDWQANC(NAL) = ALGAES(NAL).WQANCA     
  !    if( ALGAES(NAL).ISILICA == 1 )then    
  !      SMWQASC = ALGAES(NAL).WQASC
  !      DIATOM  = ALGAES(NAL).IDN
  !    endif
  !  Enddo
  !  call Broadcast_Array (SMDWQANC,   master_id)
  !  call Broadcast_Scalar(SMWQASC,    master_id)
  !  call Broadcast_Scalar(DIATOM,     master_id)
  !End if

  ! **************************************************************************
  ! *** Set up look-up table for temperature dependency over -10 °C to 60 °C
  WTEMP = WQTDMIN
  
  do M = 1, NWQTD
    ! *** Reference temperature for algae growth
    do NAL = 1, NALGAE
      WQTDG(M,NAL) = 1.
      if( WTEMP < ALGAES(NAL).WQTM1 )then
        WQTDG(M,NAL) = EXP(-ALGAES(NAL).WQKG1*(WTEMP - ALGAES(NAL).WQTM1)*(WTEMP - ALGAES(NAL).WQTM1) )
      endif
      if( WTEMP > ALGAES(NAL).WQTM2 )then
        WQTDG(M,NAL) = EXP(-ALGAES(NAL).WQKG2*(WTEMP - ALGAES(NAL).WQTM2)*(WTEMP - ALGAES(NAL).WQTM2) )
      endif
    enddo
    
    ! *** Temperature adjustment for algae metabolism & Predation
    do NAL = 1, NALGAE 
      WQTDR(M,NAL) = EXP( ALGAES(NAL).WQKTB*(WTEMP - ALGAES(NAL).WQTR) ) 
    enddo
    
    ! *** Temperature adjustment for winter bloom.  Configure WQTP1/2 & WQKP1.1 so WQTDP is low/zero at low temps
    do NAL = 1,NALGAE
      WQTDP(M,NAL) = 1.
      if( ALGAES(NAL).ISBLOOM > 0 )then
        if( WTEMP < ALGAES(NAL).WQTP1 )then
          WQTDP(M,NAL) = EXP(-ALGAES(NAL).WQKP1*(WTEMP - ALGAES(NAL).WQTP1)*(WTEMP - ALGAES(NAL).WQTP1))
        endif
        if( WTEMP > ALGAES(NAL).WQTP2 )then
          WQTDP(M,NAL) = EXP(-ALGAES(NAL).WQKP2*(WTEMP - ALGAES(NAL).WQTP2)*(WTEMP - ALGAES(NAL).WQTP2))
        endif
      endif
    enddo
    
    WTEMP = WTEMP + WQTDINC
  enddo
  
  if( process_id == master_id )then
    write(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for Biota Temperature Dependency",                               &
                                             "WQTDMIN = ", WQTDMIN, "WQTDMAX = ", WQTDMAX, "WQTDINC = ", WQTDINC, "NWQTD = ", NWQTD
    write(2,'(/,A5,A10,30(" |",I2,3A10))') "IT", "TEMP", (NAL, "WQTDG", "WQTDR", "WQTDP", NAL = 1,NALGAE)
                                            
    do M = 1,NWQTD
      write(2,'(I5,F10.3,30(4X,3F10.5))') M, WQTDTEMP(M), (WQTDG(M,NAL), WQTDR(M,NAL), WQTDP(M,NAL), NAL = 1,NALGAE)
    enddo
  endif
  
  ! *** C44
  WQCHL(1,1) = 0.0
  if( ISWQLVL == 0 )then
    WQCHL(1,1) = WQV(1,1,1)*ALGAES(1).WQCHLA + WQV(1,1,2)*ALGAES(2).WQCHLA + WQV(1,1,3)*ALGAES(3).WQCHLA
  else
    do NAL = 1,NALGAE
      if( ALGAES(NAL).ISMOBILE )then
        ! *** Phytoplankton
        WQCHL(1,1) = WQCHL(1,1) + WQV(1,1,19+NAL)*ALGAES(NAL).WQCHLA
      endif
    enddo
  endif

  do K = 1,KC
    WQCHL(LC,K) = WQCHL(1,1)
  enddo
  do K = 1,KC
    do L = 2,LA
      WQCHL(L,K) = WQCHL(1,1)
    enddo
  enddo
  
  ! *** Read Zonal biota class values
  if( IWQZONES > 0 )then
    ! *** Defaults
    ! *** WQKMV    = 0.25     ! *** KMV velocity half-saturation
    ! *** WQKMVMIN = 0.15     ! *** Minimum velocity
    ! *** WQKBP    = 6.5      ! *** M density half saturation
    ! *** WQKMVA   = 1.0      ! *** Five-parameter logistic - A
    ! *** WQKMVB   = 12.0     ! *** Five-parameter logistic - B
    ! *** WQKMVC   = 0.3      ! *** Five-parameter logistic - C
    ! *** WQKMVD   = 0.25     ! *** Five-parameter logistic - D
    ! *** WQKMVE   = 2.0      ! *** Five-parameter logistic - E
    
    if( process_id == master_id )then
      write(*,'(A)')' WQ: WQ_BIO_ZONES.JNP'
      write(2,'(1/,A)')'Reading WQ_BIO_ZONES.JNP - Zonal kinetics for Biota parameters'
      json_data => fson_parse("wq_bio_zones.jnp")  
      call fson_get(json_data, "groups", phytogroups)  
      
      ! *** Macrophytes, periphyton and phytoplankton
      do NAL = 1, fson_value_count(phytogroups)
          item => fson_value_get(phytogroups, NAL)
          
          ! *** BIO = "index"
          call fson_get(item, "algae_dynamic.max_growth_rate",                 ALGAES(NAL).WQPMA)      ! *** 
          call fson_get(item, "algae_dynamic.predation(death)_rate",           ALGAES(NAL).WQPRRA)     ! *** 
          call fson_get(item, "algae_dynamic.metabolism_rate",                 ALGAES(NAL).WQBMRA)     ! *** 
          call fson_get(item, "growth_optimal_depth",                          ALGAES(NAL).WQDOP)      ! *** 
          call fson_get(item, "plant_density_half_saturation_factor",          ALGAES(NAL).WQKBP)      ! *** 
          call fson_get(item, "settling_velocity",                             ALGAES(NAL).WQWS)       ! *** 
          call fson_get(item, "minimum_biomass_zone",                          ALGAES(NAL).WQBMIN)

          call fson_get(item, "hydrolysis_constant_to_macrophytes",            ALGAES(NAL).WQKDCALM)   ! *** 
          call fson_get(item, "half_saturation_DOC_excretion",                 ALGAES(NAL).WQKHRA)     ! *** 
          call fson_get(item, "velocity_limitation.half_saturation_velocity",  ALGAES(NAL).WQKMV)      ! *** 
          call fson_get(item, "velocity_limitation.mininum_velocity",          ALGAES(NAL).WQKMVMIN)   ! *** 

          call fson_get(item, "velocity_limitation.logistic_function.param_a", ALGAES(NAL).WQKMVA)     ! *** 
          call fson_get(item, "velocity_limitation.logistic_function.param_b", ALGAES(NAL).WQKMVB)     ! *** 
          call fson_get(item, "velocity_limitation.logistic_function.param_c", ALGAES(NAL).WQKMVC)     ! *** 
          call fson_get(item, "velocity_limitation.logistic_function.param_d", ALGAES(NAL).WQKMVD)     ! *** 
          call fson_get(item, "velocity_limitation.logistic_function.param_e", ALGAES(NAL).WQKMVE)     ! *** 
      enddo
    endif   ! *** End of master_id block
    
    ! *** Broadcast
    do NAL = 1, NALGAE
      call Broadcast_Array(ALGAES(NAL).WQPMA,    master_id)
      call Broadcast_Array(ALGAES(NAL).WQPRRA,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQBMRA,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQDOP,    master_id)
      call Broadcast_Array(ALGAES(NAL).WQKBP,    master_id)
      call Broadcast_Array(ALGAES(NAL).WQWS,     master_id)
      call Broadcast_Array(ALGAES(NAL).WQBMIN,   master_id)
      
      call Broadcast_Array(ALGAES(NAL).WQKDCALM, master_id)
      call Broadcast_Array(ALGAES(NAL).WQKHRA,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQKMV,    master_id)
      call Broadcast_Array(ALGAES(NAL).WQKMVMIN, master_id)
      
      call Broadcast_Array(ALGAES(NAL).WQKMVA,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQKMVB,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQKMVC,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQKMVD,   master_id)
      call Broadcast_Array(ALGAES(NAL).WQKMVE,   master_id)
    enddo
    
  endif  

  ! *** Initialize 3D settling arrray with user defined settling rates
  do NAL = 1, NALGAE
    do L = 2, LA
      do K = 1, KC
        ALGAES(NAL).SETTLING(L,K) = ALGAES(NAL).WQWS(IWQZMAP(L,K))     ! *** 3D settling velocities (m/s)
        if( IVARSETL == 3 )then
          CELLDENSLIGHT(L,K,NAL) = ALGAES(NAL).DENSITYINI
          CELLDENS(L,K,NAL) = ALGAES(NAL).DENSITYINI
        endif
      enddo
    enddo
  enddo
  
END SUBROUTINE ALGAECONTROL
  
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine ALGAE_SETTLING
!
!> @details  Compute the algal settling velocities for class NAL
!
!---------------------------------------------------------------------------! 
SUBROUTINE ALGAE_SETTLING(NAL, L)
  
  integer, intent(IN) :: NAL, L
  integer :: K, ipmc
  real :: DEPTH, DENDELTA, EXTC, RADLAYER, VISCW
  real, external :: ViscosityW
  
  ! *** Set the settling rate (m/day) for the current cell and layer
  if( ALGAES(NAL).ISVARSETTLE == 1 )then
    ! *** User defined settling amplitude varying by time of day
    K = KC
    ALGAES(NAL).SETTLING(L,K) = ALGAES(NAL).AMPLITUDE*2.*PI*COS((TIMEDAY-DAYINT)*2.*PI + ALGAES(NAL).PHASESHIFT)
    do K = KSZ(L), KC-1
      ALGAES(NAL).SETTLING(L,K) = ALGAES(NAL).SETTLING(L,KC)
    enddo
    
  elseif( ALGAES(NAL).ISVARSETTLE == 2 )then
    ! *** User defined settling amplitude varying by time of day and available light (Belov & Giles, 1997)
    if( RADTOP(L,KC) >= ALGAES(NAL).LIGHTMIN )then
      do K = KSZ(L),KC
        DEPTH = 1. - ZZ(L,K)          ! *** Sigma at midpoint of layer
        DEPTH = DEPTH*HP(L)           ! *** Depth at midpoint of layer
        EXTC = RADKE(L,K)
        ALGAES(NAL).SETTLING(L,K) = ALGAES(NAL).AMPLITUDE*2.*PI*COS((TIMEDAY-DAYINT)*2.*PI + ALGAES(NAL).PHASESHIFT)*EXP(-EXTC*DEPTH)
      enddo
    else
      K = KC
      ALGAES(NAL).SETTLING(L,KC) = ALGAES(NAL).AMPLITUDE*2.*PI*COS((TIMEDAY-DAYINT)*2.*PI + ALGAES(NAL).PHASESHIFT)
      do K = KSZ(L),KC-1
        ALGAES(NAL).SETTLING(L,K) = ALGAES(NAL).SETTLING(L,KC)
      enddo
    endif
      
  elseif( ALGAES(NAL).ISVARSETTLE == 3 )then
    ! ***Variable settling using cell density (Visser, etal 1997)
    do K = KSZ(L), KC  
      !RADLAYER = min(0.5*(RADTOP(L,K) + RADBOT(L,K)),350.4)            ! *** Solar radiation (W/m2) at layer midpoint with a max of 1600 umol/m2/s 
      RADLAYER = 0.5*(RADTOP(L,K) + RADBOT(L,K))                        ! *** Solar radiation (W/m2) at layer midpoint
      if( RADLAYER >= ALGAES(NAL).LIGHTMIN )then
        DENDELTA = (ALGAES(NAL).C1*RADLAYER*EXP(-RADLAYER/ALGAES(NAL).LIGHTSAT) + ALGAES(NAL).C2)*DTWQ*86400.
        CELLDENSLIGHT(L,K,NAL) = CELLDENS(L,K,NAL) + DENDELTA                                         ! *** Save for density decrease
      else
        DENDELTA = (ALGAES(NAL).F1*(CELLDENSLIGHT(L,K,NAL) + 65.) + ALGAES(NAL).F2)*DTWQ*86400.       ! *** Rate of density decrease
      endif
      
      CELLDENS(L,K,NAL) = CELLDENS(L,K,NAL) + DENDELTA
      CELLDENS(L,K,NAL) = min(CELLDENS(L,K,NAL), ALGAES(NAL).DENSITYMAX)
      CELLDENS(L,K,NAL) = max(CELLDENS(L,K,NAL), ALGAES(NAL).DENSITYMIN)
      
      ! *** Compute settling velocity using Stokes
      VISCW = ViscosityW(TEM(L,KC))                                    ! *** Dynamic viscosity of water (kg/m/s)
      ALGAES(NAL).SETTLING(L,K) = 2.*G * ALGAES(NAL).CELLRAD**2 * ALGAES(NAL).CELL2COL * (CELLDENS(L,K,NAL) - RHOW(L,K)) / 9. / ALGAES(NAL).CELLDRAG / VISCW
      ALGAES(NAL).SETTLING(L,K) = ALGAES(NAL).SETTLING(L,K)*86400.     ! *** Convert to m/day
    enddo
    
  endif
    
END SUBROUTINE ALGAE_SETTLING
  
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine Feedback
!
!> @details Set parameters for vegetative drag to be used in CALEXP and CALTBXY
!---------------------------------------------------------------------------!
!> @author P.M. Craig
!---------------------------------------------------------------------------!
SUBROUTINE feedback
  
  integer :: NAL

  do NAL = 1,NALGAE
    !IF( .not. ALGAES(NAL).ISMOBILE .and. K == KSZ(L) )then
    !  ! *** Macrophytes and periphyton
    !  WQA6A(NAL) = ALGAES(NAL).WQFCDB + (1.0 - ALGAES(NAL).WQFCDB)*(ALGAES(NAL).WQKHRA(IZ)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)) 
    !  WQA6 = WQA6 + ( WQA6A(NAL)*WQBM(L,NAL) + ALGAES(NAL).WQFCDP*WQPR(L,NAL) )*WQVO(L,K,19+NAL) 
    !ENDIF
  enddo
    
END SUBROUTINE feedback
  
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine Macro_Veg
!
!> @details Grow and shrink macrophyte vegetation length and thickness
!
!---------------------------------------------------------------------------!
!> @author P.M. Craig
!---------------------------------------------------------------------------!
SUBROUTINE Macro_Veg(NAL)
  
  integer, intent(IN) :: NAL
    
  integer :: L, K, K1, K2, KTOP
  real :: LAYER_THRESHOLD, THICK, EXCESS_MASS, TVAL1, TVALHEI
    
  do L = 2,LA
    ! *** Set active layers               delme - TODO - handle growth downward from suspended base
    do K = KSZ(L),KC
      if( HP(L)*Z(L,K) >= ALGAES(NAL).BASEDEPTH )then
        LAYERTOP(NAL,L) = K                                                    ! *** Top active layer
        exit                                                                   ! *** Jump out of the layer loop
      endif
    enddo
    
    ! *** Set vegetation height/diameter
    MAC_CELL(L) = .FALSE.                                                      ! *** Reset total vegetation height flag
    HEIGHT_MAC(L,NAL) = 0.0
    LayerRatio_MAC(L,:,NAL) = 0.0
    DIAMETER_MAC(L,NAL) = ALGAES(NAL).MINDIAMETER
    
    KTOP = 0
    K2 = min(LAYERTOP(NAL,L)+1,KC)
    ! *** Loop over layers to redistribute vegetation mass
    do K = LAYERBOT(NAL,L), K2
      if( WQV(L,K,19+NAL) < 1E-8 ) CYCLE
      
      ! *** Allow for growth between layers, if allowed
      if( WQV(L,K,19+NAL) > ALGAES(NAL).THRESHOLD .and. K < KC )then
        ! *** Redistribute mass to other layers, i.e. "Grow" any excess mass into the layer above
        !LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K), 0.1)               ! *** Handle thin layers (g C/m3)
        TVAL1 = WQV(L,K,19+NAL) - ALGAES(NAL).THRESHOLD
        TVAL1 = TVAL1*HPK(L,K)                                                 ! *** Convert to mass of macrophyte carbon to move/grow from layer
        TVAL1 = TVAL1*HPKI(L,K+1)                                              ! *** Convert back to concentration based on layer thickness above.
        WQV(L,K,19+NAL)   = ALGAES(NAL).THRESHOLD
        WQV(L,K+1,19+NAL) = WQV(L,K+1,19+NAL) + TVAL1

        WQVO(L,K,19+NAL)   = ALGAES(NAL).THRESHOLD
        WQVO(L,K+1,19+NAL) = WQVO(L,K+1,19+NAL) + TVAL1
        !ELSE
          !! *** Macrophyte is at its full depth but mass is larger than threshold.
          !! *** Therefore, increase stem diameter to account for mass
          !THICK = Z(L,K) - Z(L,LAYERBOT(NAL,L)-1)
          !IF( THICK <= 0. ) THICK = 1.0                                             ! *** Handle single layer models
          !EXCESS_MASS = WQV(L,K,19+NAL) - ALGAES(NAL).THRESHOLD
          !
          !! *** Distribute excess mass
          !DO K1 = LAYERBOT(NAL,L),K
          !  LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K1), 0.1)             ! *** Handle thin layers
          !  if( K1 < LAYERTOP(NAL,L) )then
          !    WQV(L,K1,19+NAL)  = WQV(L,K1,19+NAL) + EXCESS_MASS*DZC(L,K1)/THICK    ! *** Update layer concentrations
          !    WQVO(L,K1,19+NAL) = WQV(L,K1,19+NAL)                                  ! *** Update layer concentrations
          !  endif
          !  DIAMETER_MAC(L,NAL) = DIAMETER_MAC(L,NAL)*EXCESS_MASS/LAYER_THRESHOLD/ALGAES(NAL).STEMDENSITY   ! *** Update stem diameter
          !  !IF( DIAMETER_MAC(L,NAL) < ALGAES(NAL).MINDIAMETER )then
          !  !  PRINT '(I6,2I5,F8.3,4F10.3)', NITER, NAL, L, HP(L), ALGAES(NAL).THRESHOLD, LAYER_THRESHOLD, DIAMETER_MAC(L,NAL), ALGAES(NAL).MINDIAMETER  ! DELME
          !  !ENDIF
          !  DIAMETER_MAC(L,NAL) = max(DIAMETER_MAC(L,NAL),ALGAES(NAL).MINDIAMETER)
          !ENDDO
          !exit
        !ENDIF
      endif
      KTOP = K
    enddo
    K1 = KTOP
    
    ! *** Update vegetation heights
    do K = K1, LAYERBOT(NAL,L), -1
      if( WQV(L,K,19+NAL) > 1E-8 )then
        !LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K), 0.1)             ! *** Handle thin layers (g C/m3)
        if( K >= LAYERBOT(NAL,L)+1 )then
          if( WQV(L,K-1,19+NAL) < ALGAES(NAL).THRESHOLD )then
            ! *** Move mass to layer below
            TVAL1 = ALGAES(NAL).THRESHOLD - WQV(L,K-1,19+NAL)
            TVAL1 = max(TVAL1,0.0)
            TVAL1 = TVAL1*HPK(L,K-1)                                           ! *** Convert to mass of macrophyte carbon to move/grow from layer
            TVAL1 = TVAL1*HPKI(L,K)                                            ! *** Convert back to concentration based on layer thickness above.
            WQV(L,K-1,19+NAL) = ALGAES(NAL).THRESHOLD
            WQV(L,K,19+NAL)   = WQV(L,K,19+NAL) - TVAL1
            if( WQV(L,K,19+NAL) < 1E-8 )then
              WQV(L,K,19+NAL) = 0.0
              KTOP = max(KTOP-1, LAYERBOT(NAL,L))
            endif
            
            WQVO(L,K-1,19+NAL) = ALGAES(NAL).THRESHOLD
            WQVO(L,K,19+NAL)   = WQVO(L,K,19+NAL) - TVAL1
            WQVO(L,K,19+NAL)   = max(WQVO(L,K,19+NAL), 0.0)
          endif
        endif
        
        ! *** Complete or partial penetration of current layer
        TVALHEI = min(WQV(L,K,19+NAL)/ALGAES(NAL).THRESHOLD, 1.0)              ! *** Height of macrophyte in current layer   DELME - /ALGAES(NAL).STEMDENSITY HERE AND EVERYWHERE ?
        HEIGHT_MAC(L,NAL) = HEIGHT_MAC(L,NAL) + TVALHEI*HPK(L,K)               ! *** Total height of macrophyte
        LayerRatio_MAC(L,K,NAL) = TVALHEI                                      ! *** Fraction of layer containing vegetation  
        if( LayerRatio_MAC(L,K,NAL) < 0.01 ) LayerRatio_MAC(L,K,NAL) = 0.0     ! *** Ignore very small heights for hydrodynamics
      endif
    enddo
  
    ! *** 
    if( KTOP == 0 )then
      KTOP = LAYERBOT(NAL,L)
      LAYERTOP(NAL,L) = KTOP
      K = KTOP
      !LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K), 0.1)               ! *** Handle thin layers (g C/m3)
      
      ! *** Partial penetration of current layer
      TVALHEI = min(WQV(L,K,19+NAL)/ALGAES(NAL).THRESHOLD, 1.0)                ! *** Height of macrophyte in current layer   DELME - /ALGAES(NAL).STEMDENSITY HERE AND EVERYWHERE ?
      HEIGHT_MAC(L,NAL) = HEIGHT_MAC(L,NAL) + TVALHEI*HPK(L,K)                 ! *** Total height of macrophyte
      LayerRatio_MAC(L,K,NAL) = TVALHEI                                        ! *** Fraction of layer containing vegetation  
      if( LayerRatio_MAC(L,K,NAL) < 0.01 ) LayerRatio_MAC(L,K,NAL) = 0.0       ! *** Ignore very small heights for hydrodynamics
    else
      LAYERTOP(NAL,L) = KTOP      
    endif
          
    if( HEIGHT_MAC(L,NAL) > 0.0 ) MAC_CELL(L) = .TRUE.                         ! *** Total vegetation height is used as a flag to compute turbulence
    
          ! *** ACCUMULATE ACTIVE VEGETATION LAYERS
          !VEGK(L) = VEGK(L) + HVGTC*HPKI(L,K)

  enddo        ! *** End of L index loop  
  
END SUBROUTINE Macro_Veg
  
END MODULE
  
