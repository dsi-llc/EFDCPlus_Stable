! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Variables_WQ
  use GLOBAL
  implicit none

  Character(len = 5) :: WQCONSTIT(50)
  Character(len = 5) :: ZOONAME(20)
  
  ! *** Variables for new Water Quality features
  integer, parameter :: MAXWQ = 50 !< Maximum number of WQ components
  integer, parameter :: NTSWQVM = 50    !< Maximum number of WQ time series
  integer :: NALGAE          !< Number of algal groups
  integer :: NALGAEM         !< Maximum of NALGAE or 1 for allocation of arrays
  integer :: ALG_COUNT       !< Number of algae available
  integer :: IWQZPL          !< Switch to Zooplankton module
  integer :: NZOOPL          !< Number of zooplankton groups modeling
  integer :: NWQVZ           !< Number of WQ components excluded zooplankton
  integer :: ISMOB(20)       !< Temporary mobile flag
  integer :: IVARSETL        !< Temporary algae variable settling flag
  
  ! *** list of index for WQ components
  integer :: IROC, ILOC, IDOC, IROP, ILOP, IDOP, IP4D, IRON, ILON, IDON, INHX,  &
             INOX, ISUU, ISAA, ICOD, IDOX, ITAM, IFCB, ICO2, ICHC, ICHD, ICHG
  
  ! *** Water Column Eutrophication variables
  
  ! *** Integer Variables
  integer :: NWQV            !< Number of water quality variables, including active and inactive consitutuents
  integer :: NWQVM           !< Maximum number of water quality variables, including active and inactive consitutuents
  integer :: IWQZONES = 0    !< Flag to activact the use of WQ zones
  integer :: IWQBEN = 0      !< Sediment diagensis flux option

  integer :: DIATOM          !< Constituent index of diatom using silica (Deprecated - Now multiple classes can interact with silica)
  integer :: ISSDBIN         !< 
  integer :: ISWQLVL         !< Global flag to active macrophyte and/or priphyton velocity limitation
  integer :: ITNWQ           !< Water quality kinetics counter
  integer :: IWQICI          !< Eutrophication initial condition option: 0 - Constant, 1 - spatially varying depth averaged (WQICI), 2 - Horizontal and vertical varying (WQ_WCRST)
  integer :: IWQM
  integer :: IWQNC
  integer :: IWQNPL
  integer :: IWQONC
  integer :: IWQORST
  integer :: IWQPSL
  integer :: IWQSI
  integer :: IWQSRP
  integer :: IWQSTOX
  integer :: IWQSUN
  integer :: IWQTS
  integer :: IWQTSB
  integer :: IWQTSDT
  integer :: IWQTSE

  integer :: NDWQPSR
  integer :: NDLTAVG
  integer :: NDLTCNT
  integer :: NDDOCNT
  integer :: NFIXED          !< Number of fixed biota classes.  Replaces old IDNOTRVA
  integer :: NPSTMSR         !< Number of mass loading time series
  integer :: NTSWQV          !< Number of WQ time series
  integer :: NWQCSRM
  integer :: NWQOBE
  integer :: NWQOBN
  integer :: NWQOBS
  integer :: NWQOBW
  integer :: NWQPS
  integer :: NWQPSM          !< Number of cells for point source mass loading
  integer :: NWQPSRM         !< Maximum number of cells for point source mass loading
  integer :: NWQTD           !< Number of temperature lookup table points
  integer :: NWQTDM          !< Maximum number of temperature lookup table points
  integer :: NWQTS           !< Number of WQ time series output locations (not used by EEMS)
  integer :: NWQTSM          !< Maximum number of WQ time series output locations (not used by EEMS)
  integer :: NWQZ            !< Number of water quality zones
  integer :: NWQZM           !< Maximum number of water quality zones
  integer :: IDOSFRM         !< Formulation for D.O. saturation: 0 = Garcia and Gordon (1992), 1 = Chapra et al. (1997), 2 = Genet et al. (1974)  
  integer :: IDOSELE         !< Elevation adjustment for D.O. saturation flag: 0 = Not used, 1 = Chapra et al. (1997), 2 = Zison et al. (1978)

  integer :: ISKINETICS(MAXWQ)    !< Eutrophication parameter activation flag
  integer :: ISTRWQ(MAXWQ)        !< Eutrophication parameter transport flag
  
  integer,allocatable,dimension(:)   :: ICPSL         !< 
  integer,allocatable,dimension(:)   :: JCPSL         !< 
  integer,allocatable,dimension(:)   :: KCPSL         !< WQ PSL

  integer,allocatable,dimension(:)   :: IWQCBE        !< East  open BC WQ list: I Index 
  integer,allocatable,dimension(:)   :: IWQCBN        !< North open BC WQ list: I Index
  integer,allocatable,dimension(:)   :: IWQCBS        !< South open BC WQ list: I Index
  integer,allocatable,dimension(:)   :: IWQCBW        !< West  open BC WQ list: I Index
  integer,allocatable,dimension(:,:) :: IWQOBE        !< East  open BC concentration series index
  integer,allocatable,dimension(:,:) :: IWQOBN        !< North open BC concentration series index
  integer,allocatable,dimension(:,:) :: IWQOBS        !< South open BC concentration series index
  integer,allocatable,dimension(:,:) :: IWQOBW        !< West  open BC concentration series index
  integer,allocatable,dimension(:,:) :: IWQPSC        !< 
  integer,allocatable,dimension(:,:) :: IWQPSV        !< 
  integer,allocatable,dimension(:,:) :: IWQZMAP       !< WQ zone loopup by L and K
  integer,allocatable,dimension(:)   :: ISUDPC        !< Water quality constituent activation flag

  integer,allocatable,dimension(:)   :: JWQCBE        !< East  open BC WQ list: J Index
  integer,allocatable,dimension(:)   :: JWQCBN        !< North open BC WQ list: J Index
  integer,allocatable,dimension(:)   :: JWQCBS        !< South open BC WQ list: J Index
  integer,allocatable,dimension(:)   :: JWQCBW        !< West  open BC WQ list: J Index

  integer  :: NWQCSR(NTSWQVM)

  ! *** Real Variables
  real :: BENDAY
  real :: DAYINT
  real :: DTWQ
  real :: DTWQO2
  real :: DTWQZO2
  real :: DOELEV    = 0.0
  real :: PARADJ    = 0.43   ! Photosynthetically active solar radiation fraction
  real :: WQAOCR    = 0.0
  real :: WQCIA     = 0.0
  real :: WQCIB     = 0.0
  real :: WQCIC     = 0.0
  real :: WQCIM     = 0.0
  real :: WQI0      = 0.0
  real :: WQI1      = 0.0
  real :: WQI2      = 0.0
  real :: WQI3      = 0.0
  real :: WQKECHL   = 0.0
  real :: WQKECHLE  = 0.0
  real :: WQKEPOM   = 0.0
  real :: WQKEDOM   = 0.0
  real :: WQKEPOC   = 0.0
  real :: WQKETSS   = 0.0
  real :: WQKINUPT  = 0.0
  real :: WQTDMIN   = 0.0
  real :: WQTDMAX   = 0.0
  real :: WQTDINC   = 0.0
  real :: WQTSDT    = 0.0
  real :: WQTRHDR   = 0.0
  real :: WQKTHDR   = 0.0
  real :: WQTRMNL   = 0.0
  real :: WQKTMNL   = 0.0
  real :: WQTNIT    = 0.0
  real :: WQKN1    = 0.0
  real :: WQKN2    = 0.0
  real :: WQKSU    = 0.0
  real :: WQTRSUA    = 0.0
  real :: WQKTSUA    = 0.0
  real :: WQBFTAM    = 0.0
  real :: WQTTAM    = 0.0
  real :: WQKTAM    = 0.0
  real :: WQTRCOD    = 0.0
  real :: WQKTCOD    = 0.0
  
  real,allocatable,dimension(:)     :: H2WQ          !< DEPTH IN WATER QUALITY MODEL
  real,allocatable,dimension(:)     :: HWQ           !< CELL CENTER DEPTH FOR WATER QUALITY
  real,allocatable,dimension(:)     :: HWQI          !< 1/HWQ

  real,allocatable,dimension(:)     :: SWQ           !< Temporary variable for salinity for current layer
  real,allocatable,dimension(:)     :: TWQ           !< Temporary variable for water temperature for current layer
  real,allocatable,dimension(:)     :: XSMO20
  real,allocatable,dimension(:)     :: VOLWQ         !< Temporary variable for inverse of layer volume (1/m3)
  real,allocatable,dimension(:)     :: WQAPC         !< Algal Phosphorus-to-Carbon ratio
  real,allocatable,dimension(:,:)   :: WQATM         !< Wet Depsition by WQ zone and WQ components
  real,target,allocatable,dimension(:) :: WQBFCOD    !< COD Sediment flux 
  real,target,allocatable,dimension(:) :: WQBFNH4    !< NH4 Sediment flux
  real,target,allocatable,dimension(:) :: WQBFNO3    !< NO3 Sediment flux
  real,target,allocatable,dimension(:) :: WQBFO2     !< SOD Sediment flux
  real,target,allocatable,dimension(:) :: WQBFPO4D   !< PO4 Sediment flux
  real,target,allocatable,dimension(:) :: WQBFSAD    !< Silica Sediment flux
  real,allocatable,dimension(:,:)   :: WQCHL         !< Total clhlorphyll in the water column
  real,allocatable,dimension(:)     :: WQDFLC        !< Depositional flux of LPOC 
  real,allocatable,dimension(:)     :: WQDFLN        !< Depositional flux of LPON
  real,allocatable,dimension(:)     :: WQDFLP        !< Depositional flux of LPOP
  real,allocatable,dimension(:)     :: WQDFRC        !< Depositional flux of RPOC 
  real,allocatable,dimension(:)     :: WQDFRN        !< Depositional flux of RPON
  real,allocatable,dimension(:)     :: WQDFRP        !< Depositional flux of RPOP
  real,allocatable,dimension(:)     :: WQDFSI        !< Depositional flux of SI
  real,allocatable,dimension(:)     :: WQKEB         !< Background extinction coefficient for WQ calculations
  real,allocatable,dimension(:,:,:) :: WQOBCE        !< Open BC concentration series
  real,allocatable,dimension(:,:,:) :: WQOBCN        !< Open BC concentration series
  real,allocatable,dimension(:,:,:) :: WQOBCS        !< Open BC concentration series
  real,allocatable,dimension(:,:,:) :: WQOBCW        !< Open BC concentration series
  real,allocatable,dimension(:,:)   :: WQPO4D        !< Bioavailable (i.e. dissolved) phase to total PO4
  real,allocatable,dimension(:,:)   :: WQSAD         !< Bioavailable (i.e. dissolved) phase to total SI
  real,allocatable,dimension(:)     :: WQTDTEMP      !< Uniform temperature lookup table for all temperature dependent WQ processes
  real,allocatable,dimension(:,:,:) :: WQWPSL        !< Mass loading of WQ parameters for point sources (only IWQPSL /= 2)
  real,allocatable,dimension(:,:)   :: WQWPSLC       !< Concentrations of WQ parameters for point sources


  ! *** Algae/Macrophyte variables
  type ALGAECLASS
    character(LEN = 20) :: PHYTONAME
    integer :: IDN  !< Algal ID group number, 1 - cyanobacteria, 2 - algae diatoms, 3 - greens algae, 4 - macroalgae
    integer :: ISTOX                                 !< Flag for salinity toxicity (use for Cyanobacteria)
    integer :: ISBLOOM                               !< Flag for bloom in winter (use for Diatoms)
    integer :: ISILICA                               !< Flag for silica active (use for Diatoms)
    integer :: IWQVLIM                               !< Velocity limitation options (use for Macrophytes)
    logical :: ISMOBILE                              !< Flag to indicate mobility
    
    real :: WQBMAX         !< Maximum biomass for macrophyte    (g C/m3)
    real :: WQTM1          !< Lower optimal temperature for algal group growth (degC)
    real :: WQTM2          !< Upper optimal temperature for algal group growth (degC)
    real :: WQKG1          !< Suboptimal temperature effect coef. for algal group growth
    real :: WQKG2          !< Superoptimal temperature effect coef. for algal group growth
    real :: WQTR           !< Reference temperature for algal group metabolism (degC)
    real :: WQKTB          !< Temperature effect coef. for algal group metabolism
    real :: WQTP1          !< Lower optimal temperature for algal group predation (degC)
    real :: WQTP2          !< Upper optimal temperature for algal group predation (degC)
    real :: WQKP1          !< Suboptimal temperature effect coef. for algal group predation
    real :: WQKP2          !< Superoptimal temperature effect coef. for algal group predation
    real :: WQKHNA         !< Nitrogen half-saturation for algal group        (mg/L)
    real :: WQKHPA         !< Phosphorus half-saturation for algal group      (mg/L)   
    real :: WQKHCO2        !< CO2 half-saturation consts for algal group
                           
    real :: WQCHLA         !< Carbon-to-Chlorophyll ratio for algal group (mg C/ug Chl)
    real :: WQFCRP         !< Fraction of predated Carbon produced as RPOC
    real :: WQFCLP         !< Fraction of predated Carbon produced as LPOC
    real :: WQFCDP         !< Fraction of predated Carbon produced as DOC
                           
    real :: WQFCDB         !< Fraction of basal metabolism exuded as DOC   
                           
    real :: WQFPRP         !< Fraction of predated Phosphorus produced as RPOP
    real :: WQFPLP         !< Fraction of predated Phosphorus produced as LPOP
    real :: WQFPDP         !< Fraction of predated Phosphorus produced as DOP
    real :: WQFPIP         !< Fraction of predated Phosphorus produced as Inorganic P
                           
    real :: WQFPRB         !< Fraction of metabolized Phosphorus produced as RPOP
    real :: WQFPLB         !< Fraction of metabolized Phosphorus produced as LPOP
    real :: WQFPDB         !< Fraction of metabolized Phosphorus produced as DOP
    real :: WQFPIB         !< Fraction of metabolized Phosphorus produced as P4T
                           
    real :: WQFNRP         !< Fraction of predated Nitrogen produced as RPON
    real :: WQFNLP         !< Fraction of predated Nitrogen produced as LPON
    real :: WQFNDP         !< Fraction of predated Nitrogen produced as DON
    real :: WQFNIP         !< Fraction of predated Nitrogen produced as Inorganic N
                           
    real :: WQFNRB         !< Fraction of metabolized Nitrogen produced as RPON
    real :: WQFNLB         !< Fraction of metabolized Nitrogen produced as LPON
    real :: WQFNDB         !< Fraction of metabolized Nitrogen produced as DON
    real :: WQFNIB         !< Fraction of metabolized Nitrogen produced as Inorganic N
                           
    real :: WQANCA         !< Nitrogen-to-carbon ratio for algal group (gN/gC)
                           
    real :: WQFSID         !< Fraction of metabolized Silica by algae (diatoms) produced as available Silica
    real :: WQFSPD         !< Fraction of metabolized Silica by algae (diatoms) produced as particulate biogenic Silica
    real :: WQFSIP         !< Fraction of predated Silica by algae (diatoms) produced as available Silica
    real :: WQFSPP         !< Fraction of predated Silica by algae (diatoms) produced as particulate biogenic Silica
                           
    real :: WQASC          !< Silica-to-carbon ratio for algae (diatoms)
    real :: WQKHS          !< Silica half-saturation for algae (diatoms) (mg/L)
    
    real :: WQC2DW         !< Carbon to Dry Weight Ratio
    real :: WQKEMAC        !< Light extinction coefficient due to macrophyte biomass   (m3/gC/m)
                           
    real :: WQSTOX         !< Salinity at which microcystis growth is halved (cyanobacteria)
    real :: WQAOCRP        !< Algae photosynthesis oxygen-to-carbon ratio (Macroalgae)
    real :: WQAOCRR        !< Algae respiration oxygen-to-carbon ratio (Macroalgae)
                           
    real :: WQALGOCR       !< Stoichiometric algae oxygen-to-carbon ratio (gO2/gC)
    real :: WQALGONT       !< Stoichiometric algae oxygen = to-nitrate ratio (gO2/gN)
    real :: WQAPCM         !< Factor to modify APC for macroalgae from 
     
    real :: THRESHOLD      !< Concentration factor to convert macrophyte concentration to plant height (g C/m2/(1 meter hieght (m)))
    real :: BASEDEPTH      !< Distance below the surface (i.e. depth) of the "base" of the macrophyte/periphyton growth (m)
    real :: MAXLENGTH      !< Maximum length from the "base" to allow macrophyte/periphyton growth (m)
    
    ! *** Hydrodynamic Feedback
    integer :: ISDRAG      !< Hydrodynamic feedback flag: 0 - None,  1 - Feedback activated
    integer :: ISMACL      !< Lamilar flow option: 0 - None,  1 - Included
    real :: DRAGCOEFF      !< Hydrodynamic drag coefficient for macrophyte stem and leaves
    real :: MINDIAMETER    !< Minimum diameter of each plant (m)
    real :: STEMHEIGHT     !< Total stem height of each plant (m).  Changes with growth
    real :: STEMDENSITY    !< Plant/stem density (#/m2)
    real :: ALPMAC         !< Plant/stem cell blocking factor
    
    ! *** Variable settling velocities
    integer :: ISVARSETTLE      !< Option for variable algae settling.  
                                !<    0 - use user specified settling/floating velocity, 
                                !<    1 - Velocity varyies with time of day, 
                                !<    2 - Velocity varies with time of day and available light, and 
                                !<    3 - Dynamic Velocity (Visser, et.al. 1997)
    real :: AMPLITUDE           !< Amplitude (m) of time varying settling (options 1 and 2)
    real :: PHASESHIFT          !< Phase shift for the daily calculations
    real :: DENSITYMIN          !< Minimum cell density (kg/m3)
    real :: DENSITYMAX          !< Maximum cell density (kg/m3)
    real :: DENSITYINI          !< Initial cell density (kg/m3)
    real :: CELL2COL            !< Ratio of algal cell volume to colony volume
    real :: CELLDRAG            !< Drag coefficient of algal cell/colony for Stokes
    real :: CELLRAD             !< Algal cell radius for Stokes  (m)
    real :: LIGHTMIN            !< Minimum light for algal cell density increase (W/m2)
    real :: LIGHTSAT            !< Light saturation constant for phytoplankton   (W/m2)
    real :: C1                  !< Algal cell density multiplier for I >= Ic     (s2/m3)
    real :: C2                  !< Algal cell density offset for I >= Ic         (kg/m3/s2)
    real :: F1                  !< Algal cell density slope for I < Ic           (/s)
    real :: F2                  !< Algal cell density intercept for I < Ic       (kg/m3/s)
    real :: HVM                 !< Minimum water depth for active vertical migration (m) 
    
    ! *** Zonal variables
    real,allocatable,dimension(:)     :: WQBMRA      !< Basal metabolism rate for algal group (1/day)
    real,allocatable,dimension(:)     :: WQDOP       !< Optimal depth (m) for algal group growth
    real,allocatable,dimension(:)     :: WQDRA       !< Maximum death rate for algal croup (1/day)
    real,allocatable,dimension(:)     :: WQKDCALM    !< Constant relating DOC hydrolysis rate to macroalgae
    real,allocatable,dimension(:)     :: WQKHRA      !< Half-sat. constant (gO2/m^3) for algal group DOC excretion
    real,allocatable,dimension(:)     :: WQPMA       !< Maximum growth rate for algal group (1/day)
    real,allocatable,dimension(:)     :: WQPRRA      !< Maximum predation rate on algal group (1/day)
    real,allocatable,dimension(:)     :: WQKBP       !< Growth limitation - Plant density half saturation factor   (G C/M2)
    real,allocatable,dimension(:)     :: WQKMV       !< Velocity limitation - half saturation velocity               (m/s)
    real,allocatable,dimension(:)     :: WQKMVMIN    !< Velocity limitation - Minimum velocity to apply limiation    (m/s)
    real,allocatable,dimension(:)     :: WQKMVA      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E
    real,allocatable,dimension(:)     :: WQKMVB      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    real,allocatable,dimension(:)     :: WQKMVC      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    real,allocatable,dimension(:)     :: WQKMVD      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    real,allocatable,dimension(:)     :: WQKMVE      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    real,allocatable,dimension(:)     :: WQWS        !< Settling velocity for algal group (m/day)
    real,allocatable,dimension(:)     :: WQBMIN      !< Minimum biomass for macrophyte    Input as g C/m2 but converted to g C/m3 using IC water depths

    ! *** 3D variables
    real,allocatable,dimension(:,:)   :: SETTLING    !< Settling velocity for algal group (m/day)
    
  end type 
  integer :: MACDRAG = 0
  type(ALGAECLASS),allocatable,dimension(:) :: ALGAES     !< ALGAL GROUPS
  real,allocatable,dimension(:,:,:) :: CELLDENSLIGHT      !< Last algal density when cell had light > critical light
  real,allocatable,dimension(:,:,:) :: CELLDENS           !< Current algal cell density
  
  real :: WQKRC       !< Minimum hydrolysis rate (1/day) of RPOC
  real :: WQKLC       !< Minimum hydrolysis rate (1/day) of LPOC
  real :: WQKRCALG    !< Constant relating RPOC hydrolysis rate to total Chl-
  real :: WQKLCALG    !< Constant relating LPOC hydrolysis rate to total Chl-
  real :: WQKDCALG    !< Constant relating DOC hydrolysis rate to total Chl-
  real :: WQAANOX     !< Ratio of denitrification rate to oxic DOC respiration rate
  
  real :: WQKRP       !< Minimum hydrolysis rate (1/day) of RPOP
  real :: WQKLP       !< Minimum hydrolysis rate (1/day) of LPOP
  real :: WQKDP       !< Minimum hydrolysis rate (1/day) of DOP
  real :: WQKRPALG    !< Constant relating hydrolysis rate of RPOP to algae
  real :: WQKLPALG    !< Constant relating hydrolysis rate of LPOP to algae
  real :: WQKDPALG    !< Constant relating hydrolysis rate of DOP to algae
  real :: WQCP1PRM    !< Constant used in determining algae Phos-to-Carbon ratio
  real :: WQCP2PRM    !< Constant used in determining algae Phos-to-Carbon ratio
  real :: WQCP3PRM    !< Constant used in determining algae Phos-to-Carbon ratio
  real :: WQKPO4P     !< Partition coefficient for sorbed/dissolved PO4
  
  real :: WQKRN       !< Minimum hydrolysis rate (1/day) of RPON
  real :: WQKLN       !< Minimum hydrolysis rate (1/day) of LPON
  real :: WQKDN       !< Minimum hydrolysis rate (1/day) of DON
  real :: WQKRNALG    !< Constant relating hydrolysis rate of RPON to algae
  real :: WQKLNALG    !< Constant relating hydrolysis rate of LPON to algae
  real :: WQKDNALG    !< Constant relating hydrolysis rate of DON to algae
  real :: WQKHNDO     !< Nitrification half-sat. constant for D.O.
  real :: WQKHNN      !< Nitrification half-sat. constant for NH4
  real :: WQANDC      !< Mass NO3 reduces per DOC oxidized (gN/gC)
  real :: WQNITM      !< Maximum nitrification rate (/day)
  real :: RNH4WQ_     !< Ammonia (for Current Layer)
  real :: RNO3WQ_     !< Nitrate (for Current Layer)
  
  real :: WQKSAP      !< Partition coef. for sorbed/dissolved SA
  real :: WQKHORDO    !< Oxic respiration half-sat. constant for D.O. (gO2/m^3)
  real :: WQKHDNN     !< Half-sat. constant for denitrification (gN/m^3)
  real :: WQAONT      !< Stoichiometric algae oxygen = to-nitrate ratio (gO2/gN)
  real :: WQFD        !< Solar radiation fraction daylength
  real :: WQISMIN     !< Minimum optimum solar radiation (Langley/day)
  real :: WQI0OPT     !< Optimal solar radiation
  real :: SOLSRDT     !< Daily average solar radiation interpolation for water quality
  real :: SOLFRDT     
  real :: WQKHBMF     !< D.O. concentration where TAM release is half the anoxic rate
  real :: WQTAMDMX    !< TAM solubility at anoxic conditions (mol/m^3)
  real :: WQKDOTAM    !< Constant relating TAM solubility to D.O.
  real :: WQKFCB      !< First-order fecal coliform bacteria decay rate (1/day)
  real :: WQTFCB      !< Temperature effect constant for KFCB decay rate
  real :: STEMFAC     !< SOD temperature factor, use 1.0 to ignore temperature effects  
  real :: WQTSB
  real :: WQTSE
  
  real :: WQHRAVG

  ! *** Integer Array Variables
  integer,allocatable,dimension(:)   :: IWQKA
  integer,allocatable,dimension(:,:) :: ICWQTS
  integer,allocatable,dimension(:,:) :: IBENMAP
  integer,allocatable,dimension(:)   :: NMLSER
  integer,allocatable,dimension(:)   :: MWQPTLT
  integer,allocatable,dimension(:)   :: IMWQZT
  integer,allocatable,dimension(:)   :: IMWQZT1
  integer,allocatable,dimension(:)   :: IMWQZT2
  integer,allocatable,dimension(:)   :: IWQT      !< Index for look-up table for temperature dependency
  
  ! *** Real Array Variables  
  real,allocatable,dimension(:,:)   :: WQPA       !< Net Growth rate for algal group (1/day)
  real,allocatable,dimension(:,:)   :: WQBM       !< Basal metabolism rate for algal group (1/day)
  real,allocatable,dimension(:,:)   :: WQPR       !< Predation rate on algal group (1/day)
  real,allocatable,dimension(:,:,:) :: WQBSETL    !< Vertical migration rates of phytoplankton between layers (settling or rising) (1/day)
  real,allocatable,dimension(:,:)   :: WQTDG      !< Temperature dependency of algal growth
  real,allocatable,dimension(:,:)   :: WQTDR      !< Algal metabolism/predation exponential function of a general algal group
  real,allocatable,dimension(:,:)   :: WQTDP      !< Algal predation exponential function of a general algal group - WQSKE2
  real,allocatable,dimension(:)     :: WQTDGP     !< Algal predation exponential function of a general algal group - WQSKE1
  real,allocatable,dimension(:)     :: WQTD1FCB   
  real,allocatable,dimension(:)     :: WQTD2FCB   
  real,allocatable,dimension(:)     :: UHEQ       
  real,allocatable,dimension(:)     :: WQHT       !< Fractional of Depth at the Top of the Layer
  real,allocatable,dimension(:,:,:) :: WQWDSL     
  real,allocatable,dimension(:,:)   :: DDOAVG     !< Diurnal domain analysis
  real,allocatable,dimension(:,:)   :: DDOMAX     !< Diurnal domain analysis
  real,allocatable,dimension(:,:)   :: DDOMIN     !< Diurnal domain analysis
  real,allocatable,dimension(:,:)   :: RLIGHTC    !< Light extinction analysis
  real,allocatable,dimension(:,:)   :: RLIGHTT    !< Light extinction analysis
  real,allocatable,dimension(:)     :: REAC       !< Global reaeration adjustment factor
  real,allocatable,dimension(:)     :: WQKDC      !< Minimum hydrolysis rate (1/day) of DOC
  real,allocatable,dimension(:)     :: WQKRO      !< Reaeration constant (3.933 for OConnor-Dobbins; 5.32 for Owen-Gibbs)
  real,allocatable,dimension(:)     :: WQKTR      !< Temperature rate constant for reaeration
  real,allocatable,dimension(:)     :: WQKHCOD    !< Oxygen half-saturation constant for COD decay (mg/L O2)
  real,allocatable,dimension(:)     :: WQKCD      !< COD decay rate (per day)
  real,allocatable,dimension(:)     :: WQTDHDR    !< Temperature effect on hydrolysis
  real,allocatable,dimension(:)     :: WQTDMNL    !< Temperature effect on mineralization
  real,allocatable,dimension(:)     :: WQKSUA     !< Temperature effect on PSI dissolution
  real,allocatable,dimension(:,:)   :: WQKCOD     !< Temperature effect on COD Oxidation    
  real,allocatable,dimension(:)     :: WQTDNIT    
  real,allocatable,dimension(:,:)   :: WQTDKR     !< Temperature effect Reaeration rate coefficient
  real,allocatable,dimension(:)     :: WQTDTAM    
  real,allocatable,dimension(:,:)   :: WQTAMP     
  real,allocatable,dimension(:)     :: WQWSLP     !< Settling velocity for refractory POM (m/day)
  real,allocatable,dimension(:)     :: WQWSRP     !< Settling velocity for refractory POM (m/day)
  real,allocatable,dimension(:)     :: WQWSS      !< Settling velocity for particles sorbed to TAM(m/day)
  real,allocatable,dimension(:,:,:) :: WQATML     
  real,allocatable,dimension(:)     :: XBENMUD    
  real,allocatable,dimension(:)     :: TSSRD      
  real,allocatable,dimension(:)     :: SOLFRD     
  real,allocatable,dimension(:)     :: SOLSRD     
  real,allocatable,dimension(:)     :: TCWQPSR    
  real,allocatable,dimension(:,:)   :: T_MLSER    
  real,allocatable,dimension(:,:,:) :: V_MLSER    
  real,allocatable,dimension(:,:)   :: WQPSSRT    
  real,allocatable,dimension(:)     :: WQI0BOT    
  real,allocatable,dimension(:)     :: DZWQ       !< Inverse of cell layer height
  real,allocatable,dimension(:)     :: WQKRPC     !< Hydrolysis rate of RPOC (1/day)
  real,allocatable,dimension(:)     :: WQKRPN     !< Hydrolysisrate of RPON (1/day)
  real,allocatable,dimension(:)     :: WQKRPP     !< Hydrolysisrate of RPOP (1/day)
  real,allocatable,dimension(:)     :: WQKLPC     !< Hydrolysis rate of LPOC (1/day)
  real,allocatable,dimension(:)     :: WQKLPN     !< Hydrolysisrate of LPON (1/day)
  real,allocatable,dimension(:)     :: WQKLPP     !< Hydrolysisrate of LPOP (1/day)  
  real,allocatable,dimension(:)     :: WQKHR      !< Heterotrophic respiration rate of DOC (1/day)
  real,allocatable,dimension(:)     :: WQDENIT    !< Denitrification rate of DOC (1/day)
  real,allocatable,dimension(:,:)   :: WQPN       !< Preference for ammonium uptake by algaes
  real,allocatable,dimension(:)     :: WQNIT      !< Nitrification rate (1/day)
  real,allocatable,dimension(:)     :: O2WQ       !< Dissolved oxygen concentration (mg/l)
  real,allocatable,dimension(:)     :: RNH4WQ     
  real,allocatable,dimension(:)     :: RNO3WQ     
  real,allocatable,dimension(:)     :: WQDOS      !< Saturation concentration of dissolved oxygen (gO2/m3)
  real,allocatable,dimension(:)     :: WQH10      
  real,allocatable,dimension(:,:)   :: XDOSAT     
  real,allocatable,dimension(:)     :: WQP19      
  real,allocatable,dimension(:)     :: WQRREA     !< Store reaeration rate for array out writing      
  real,allocatable,dimension(:)     :: WQKRDOS    !< Kr * DOs
  real,allocatable,dimension(:)     :: WQKK       !< Related to matrix K1
  real,allocatable,dimension(:)     :: WQRR       !< Accumulation of sources and sinks in the kinetic equation  
  real,allocatable,dimension(:,:)   :: WQRPSET    !< Settling rate of ROP
  real,allocatable,dimension(:,:)   :: WQLPSET    !< Settling rate of LOP
  real,allocatable,dimension(:,:)   :: WQWSSET    
  real,allocatable,dimension(:)     :: PO4DWQ     !< Phosphate (for Current Layer)
  real,allocatable,dimension(:)     :: RNH4NO3    !< Total Inorganic Nitrogen (for Current Layer)
  real,allocatable,dimension(:)     :: WQOBTOT    !< Total Algal Biomass (mg/l)
  real,allocatable,dimension(:)     :: WQKDON     !< Mineralization rate of DON (1/day)
  real,allocatable,dimension(:)     :: WQKDOP     !< Mineralization rate of DOP (1/day)  
  real,allocatable,dimension(:)     :: WQT10      
  real,allocatable,dimension(:)     :: WQT17      
  real,allocatable,dimension(:)     :: WQN17      
  real,allocatable,dimension(:)     :: WQO18      !< Matrix K1 in COD kinetic equation
  real,allocatable,dimension(:)     :: WQR20      
  real,allocatable,dimension(:)     :: SMAC       !< Flag to turn on (1.0) or off (0.0) fixed biota 
    
  ! *** Begin Dissolved Carbondioxide variables    VB
  real,allocatable,dimension(:)     :: CO2WQ
  real,allocatable,dimension(:)     :: CDOSATIDX
  real,allocatable,dimension(:)     :: WQCDOS
  real,allocatable,dimension(:,:)   :: WQITOP
  real,allocatable,dimension(:)     :: WQKRCDOS
  real,allocatable,dimension(:)     :: WQP22
  ! End Dissolved Carbondioxide variables  
  
  ! *** Macrophyte/Preiphyton variables
  logical,allocatable,dimension(:)   :: MAC_CELL  !< Flag indicating cell has macrophyte growth with hydrodynamic feedback
  integer,allocatable,dimension(:,:) :: LAYERBOT  !< Lowest/Bottom layer for macrophyte/periphyton growth
  integer,allocatable,dimension(:,:) :: LAYERTOP  !< Highest/Top   layer for macrophyte/periphyton growth

  real,target,allocatable,dimension(:,:)   :: HEIGHT_MAC      !< Current macrophyte height above the starting growth height (m)
  real,target,allocatable,dimension(:,:)   :: DIAMETER_MAC    !< Current macrophyte stem diameter (m)  
  real,target,allocatable,dimension(:,:,:) :: LayerRatio_MAC  !< Fraction of vegetation vertical extents in the layer
  
  real,allocatable,dimension(:,:) :: MACAD  !< Blockage ratio in Z (-)
  real,allocatable,dimension(:,:) :: MVEGZ  !< "a" - projected plant area per unit volume  (1/m)
  
  ! *** Sediment Diagenesis variables 
  integer :: ISMHYST
  integer :: ISMICI
  integer :: ISMORST
  integer :: ISMRST
  integer :: ISMTS
  integer :: ISMTSB
  integer :: ISMTSDT
  integer :: ISMTSE
  integer :: ISMZ
  integer :: ISMZB

  integer :: NSMG    !< Number of reactive glasses.  Hardwired to 3: G1, G2 and G3
  integer :: NSMGM   !< Maximum NSMG
  integer :: NSMZ    !< Number of sediment flux/diagenesis zones
  integer :: NSMZM   !< Maximum NSMZ

  real :: SMCSHSCH
  real :: SMFD1H2S
  real :: SMFD1NH4
  real :: SMO2NH4
  real :: SMO2NO3
  real :: SMPOCR

  real :: SMWQASC   !< Silica-to-carbon ratio for algae diatoms
  real :: SMDIFT    !< Diffusion coefficient for sediment temperature (m2/sec)
  real :: SMM1      !< Solid concentrations in Layer 1 (Kg/L)
  real :: SMM2      !< Solid concentrations in Layer 2 (Kg/L)   
  real :: SMKMDP    !< Particle mixing half-saturation constant for oxygen (mg/L)
  real :: SMKBST    !< First-order decay rate for accumulated benthic stress (1/day)
  real :: XSMDPMIN  !< Minimum diffusion coefficient for particle mixing (m^2/d)
  real :: SMRBIBT   !< Ratio of bio-irrigation to bioturbation (unitless)
  real :: SMO2BS    !< Critical overlying oxygen concentration below which benthic hysteresis occurs (mg/L)
  real :: SMHYLAG   !< Time duration for which the maximum or minimum stress is retained (days)
  real :: SMHYDUR   !< Critical hypoxia duration; if less than this value, no hysteresis occurs (days)
  real :: SMKMNH4   !< Nitrification half-sat. constant for ammonium (gN/m^3)
  real :: SMKMO2N   !< Nitrification half-sat. constant for dissolved oxygen (gO2/m^3)
  real :: SMP2PO4   !< Partition coefficient, ratio of particulate to dissolved PO4 in layer 2 (L/Kg)
  real :: SMCO2PO4  !< Critical dissolved oxygen for PO4 sorption (mg/L)
  real :: SMO2C     !< Stoichiometric coefficient for carbon diagenesis consumed by H2S oxidation (gO2/gC)
  real :: SMKMPSI   !< Silica dissolution half-saturation constant for PSi (g Si/m^3
  real :: SMSISAT   !< Saturation concentration of silica in pore water (g Si/m^3)
  real :: SMP2SI    !< Partition coefficient for Si in Layer 2, controls sorption of dissolved silica to solids (L/Kg)
  real :: SMDP1SI   !< Factor that enhances sorption of silica in layer 1 when D.O. exceeds DOcSi (unitless)
  real :: SMCO2SI   !< Critical dissolved oxygen for silica sorption in layer 1 (mg/L)
  real :: SMJDSI    !< Detrital flux of particulate biogenic silica from sources other than diatom algae (gSi/m^2/d)
  real :: SMFD2H2S  !< Dissolved fraction of H2S in Layer 2
  real :: SMFD2NH4  !< Dissolved fraction of NH4 in Layer 2
  real :: SMFD2PO4  !< Dissolved fraction of PO4 in Layer 2
  real :: SMFD2SI   !< Dissolved fraction of SI  in Layer 2
  real :: SMFP1H2S  !< Particulate fraction of H2S in Layer 1 
  real :: SMFP1NH4  !< Particulate fraction of NH4 in Layer 1 
  real :: SMFP2H2S  !< Particulate fraction of H2S in Layer 2 
  real :: SMFP2NH4  !< Particulate fraction of NH4 in Layer 2 
  real :: SMFP2PO4  !< Particulate fraction of PO4 in Layer 2 
  real :: SMFP2SI   !< Particulate fraction of SI  in Layer 2 
  real :: SM1OKMDP
  real :: WQTDsMIN
  real :: WQTDsMAX
  real :: WQTDsINC
  
  ! *** Integer Array Variables
  integer,allocatable,dimension(:)  :: ISMT       !<  Current Lookup Index Based On The Current Sediment Temperature
  integer,allocatable,dimension(:)  :: ISMZMAP
    
  ! *** Real Array Variables
  real,allocatable,dimension(:,:)   :: WQDFB      !< Coupling between water column and sediment diagenesis model
  real,allocatable,dimension(:)     :: SMD1PO4
  real,allocatable,dimension(:)     :: SMD1SI
  real,allocatable,dimension(:)     :: SMDD       !< Diffusion coefficient in pore water (m2/day)
  real,allocatable,dimension(:)     :: SMDFSI
  real,allocatable,dimension(:)     :: SMDGFC
  real,allocatable,dimension(:)     :: SMDGFN
  real,allocatable,dimension(:)     :: SMDGFP
  real,allocatable,dimension(:)     :: SMDP       !< Apparent diffusion coefficient for particle mixing (m2/day)
  real,allocatable,dimension(:)     :: SMDP1PO4   !< Factor to enhance sorption of PO4 in layer 1 when DO is less than DOcPO4 (unitless)
  real,allocatable,dimension(:)     :: SMDPMIN    !< Minimum diffusion coefficient for particle mixing divided by H2
  real,allocatable,dimension(:)     :: SMDTOH
  real,allocatable,dimension(:,:)   :: SMFCBA     !< fraction of POC from algae routed to G classes
  real,allocatable,dimension(:,:)   :: SMFCR      !< fraction of water column refractory POC routed to G-classes
  real,allocatable,dimension(:,:)   :: SMFNBA     !< fraction of PON from algae routed to G classes
  real,allocatable,dimension(:,:)   :: SMFNR      !< fraction of water column refractory PON routed to G-classes
  real,allocatable,dimension(:,:)   :: SMFPBA     !< fraction of POP from algae routed to G classes
  real,allocatable,dimension(:,:)   :: SMFPR      !< fraction of water column refractory POP routed to G-classes
  real,allocatable,dimension(:)     :: SMHODT     !< Benthic sediment depth divided by time step
  real,allocatable,dimension(:)     :: SMHYPD
  real,allocatable,dimension(:)     :: SMJAQH2S
  real,allocatable,dimension(:)     :: SMJDEN
  real,allocatable,dimension(:)     :: SMJGCH4
  real,allocatable,dimension(:)     :: SMJNIT 
  real,allocatable,dimension(:)     :: SMK1H2S
  real,allocatable,dimension(:)     :: SMK1NO3    !< Reaction velocity for denitrification in layer 1 at 20 degC (m/day)
  real,allocatable,dimension(:)     :: SMK2NO3    !< Reaction velocity for denitrification in layer 2 at 20 degC (m/day)
  real,allocatable,dimension(:)     :: SMKL12     !< Dissolved phase mixing coefficient
  real,allocatable,dimension(:)     :: SMKNH4     !< Optimal reaction velocity for nitrification at 20 degC (m/day)
  real,allocatable,dimension(:)     :: SMSS       !< Surface mass transfer coefficient
  real,allocatable,dimension(:)     :: SMTD1CH4
  real,allocatable,dimension(:)     :: SMTD2CH4
  real,allocatable,dimension(:,:)   :: SMTDCD
  real,allocatable,dimension(:)     :: SMTDDD     !< Constant for temperature adjustment for Dd raised to T-20
  real,allocatable,dimension(:)     :: SMTDDP     !< Constant for temperature adjustment for Dp raised to T-20
  real,allocatable,dimension(:,:)   :: SMTDND
  real,allocatable,dimension(:)     :: SMTDNH4    !< Constant for temperature adjustment for NH4 raised to T-20
  real,allocatable,dimension(:)     :: SMTDNO3    !< Constant for temperature adjustment for NO3 raised to T-20
  real,allocatable,dimension(:,:)   :: SMTDPD
  real,allocatable,dimension(:)     :: SMTDSI
  
  real,allocatable,dimension(:)     :: SMTMP
  real,allocatable,dimension(:)     :: SMW12      !< Particulate phase mixing coefficient
  real,allocatable,dimension(:)     :: SMW2       !< Sediment burial rate (cm/year)
  real,allocatable,dimension(:)     :: SMW2DTOH
  real,allocatable,dimension(:)     :: SMW2PHODT  
  real,allocatable,dimension(:)     :: SODMULT    !< Factor to enhance magnitude of sediment oxygen demand (unitless)
  real,allocatable,dimension(:)     :: SM1DIFT
  real,allocatable,dimension(:)     :: SM2DIFT
  real,allocatable,dimension(:)     :: SM2NO3     !< WATER QUALITY MODEL VARIABLE
  real,allocatable,dimension(:)     :: SM2PO4     !< WATER QUALITY MODEL VARIABLE
  real,allocatable,dimension(:)     :: SMHSED     !< WATER QUALITY MODEL VARIABLE

  real,allocatable,dimension(:,:)   :: SMPOC      !< Particulate organic Carbon     in the sediments
  real,allocatable,dimension(:,:)   :: SMPOP      !< Particulate organic Phosphorus in the sediments 
  real,allocatable,dimension(:,:)   :: SMPON      !< Particulate organic Nitrogen   in the sediments
  real,allocatable,dimension(:)     :: SMPSI      !< Sediment Silica

  real,target,allocatable,dimension(:)   :: SM1H2S !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM1NH4 !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM1NO3 !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM1PO4 !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM1SI  !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM2H2S !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM2NH4 !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SM2SI  !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SMBST  !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SMCSOD !< WATER QUALITY MODEL VARIABLE
    
  real,target,allocatable,dimension(:,:) :: SMDFC  !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:,:) :: SMDFN  !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:,:) :: SMDFP  !< WATER QUALITY MODEL VARIABLE
  real,target,allocatable,dimension(:)   :: SMNSOD !< WATER QUALITY MODEL VARIABLE
  logical,allocatable,dimension(:)       :: SMHYST !< WATER QUALITY MODEL VARIABLE
  
  real,target,allocatable,dimension(:)   :: SMT    !< Sediment temperature (degC)
    
  ! *** Zooplankton variables
  integer :: NZO                                 !< Index of zooplankton groups
  integer,allocatable,dimension(:)    :: IWQZT   !< Index for look-up table for temperature dependency
  real :: WQTDZMIN    !< Lower end of temperature for look-up table
  real :: WQTDZMAX    !< Upper end of temperature for look-up table
  real :: WTEMPZ
  real :: WQTDZINC    !< Increment of temperature for look-up table
  
  type ZOOPLGROUP
    character(LEN = 20) :: ZOONAME
    integer             :: IDZ
    integer             :: ISPREY         !< Flag to indicate prey or predator
    integer             :: ISPREDATOR     !< Flag to indicate prey or predator
    
    real :: CTZ         !< Carbon threshold for zooplankton grazing (gC/m3)
    real :: ANCZ        !< Nitrogen to Carbon ratio (gN/gC)
    real :: APCZ        !< Phosphorus to Carbon ratio (gP/gC)
    real :: ASCZ        !< Phosphorus to Carbon ratio (gSi/gC)
    real :: KHCZ        !< Prey density at which zooplankton grazing is halved (gC/m3)
    real :: ULZ         !< Utilization of labile particulate organic carbon
    real :: URZ         !< Utilization of refractory particulate organic carbon
    real :: UDZ         !< Utilization of Dissolved organic carbon
    real :: UZPL        !< Utilization of zooplankton as prey
    real :: KTGZ1       !< Effect of temperature blow optimal on grazing
    real :: KTGZ2       !< Effect of temperature above optimal on grazing
    real :: TMZG1       !< Lower end of optimal temperature for grazing
    real :: TMZG2       !< Upper end of optimal temperature for grazing
    real :: KTBZ        !< Effect of temperature on metabolism of zooplankton
    real :: TRZB        !< Reference temperature for zooplankton metabolism
    real :: KTPZ        !< Effect of temperature on predation of zooplankton
    real :: TRZP        !< Reference temperature for zooplankton predation
    real :: DOCRIT      !< Critical DO concentration below which zooplankton death occurs (g DO/m3)
    real :: DZEROZ      !< Zooplankton death at zero DO concentration
    real :: FCDDZ       !< Fraction of death carbon produced by zooplankton as DOC
    real :: FCLDZ       !< Fraction of death carbon produced by zooplankton as LPOC 
    real :: FCRDZ       !< Fraction of death carbon produced by zooplankton as LROC
    real :: FCDPZ       !< Fraction of predated carbon produced by zooplankton as DOC
    real :: FCLPZ       !< Fraction of predated carbon produced by zooplankton as LPOC
    real :: FCRPZ       !< Fraction of predated carbon produced by zooplankton as RPOC
    real :: FNLDZ       !< Fraction of death nitrogen produced by zooplankton as LPON
    real :: FNRDZ       !< Fraction of death nitrogen produced by zooplankton as RPON
    real :: FNDDZ       !< Fraction of death nitrogen produced by zooplankton as DON
    real :: FNIDZ       !< Fraction of death nitrogen produced by zooplankton as NH4
    real :: FNLPZ       !< Fraction of predated nitrogen produced by zooplankton as LPON
    real :: FNRPZ       !< Fraction of predated nitrogen produced by zooplankton as RPON
    real :: FNDPZ       !< Fraction of predated nitrogen produced by zooplankton as DON
    real :: FNIPZ       !< Fraction of predated nitrogen produced by zooplankton as NH4
    real :: FNDBZ       !< Fraction of basal metabolism nitrogen produced by zooplankton as DON
    real :: FNIBZ       !< Fraction of basal metabolism nitrogen produced by zooplankton as NH4
    real :: FPLDZ       !< Fraction of death phosphorus produced by zooplankton as LPOP
    real :: FPRDZ       !< Fraction of death phosphorus produced by zooplankton as RPOP
    real :: FPDDZ       !< Fraction of death phosphorus produced by zooplankton as DOP
    real :: FPIDZ       !< Fraction of death phosphorus produced by zooplankton as PO4
    real :: FPLPZ       !< Fraction of predated phosphorus produced by zooplankton as LPOP
    real :: FPRPZ       !< Fraction of predated phosphorus produced by zooplankton as RPOP
    real :: FPDPZ       !< Fraction of predated phosphorus produced by zooplankton as DOP
    real :: FPIPZ       !< Fraction of predated phosphorus produced by zooplankton as PO4
    real :: FPDBZ       !< Fraction of basal metabolism phosphorus produced by zooplankton as DOP
    real :: FPIBZ       !< Fraction of basal metabolism phosphorus produced by zooplankton as PO4
    real :: FSPDZ       !< Fraction of death Silica produced by zooplankton as SU
    real :: FSPPZ       !< Fraction of predated Silica produced by zooplankton as SU
    real :: FSADZ       !< Fraction of death Silica produced by zooplankton as SA
    real :: FSAPZ       !< Fraction of predated Silica produced by zooplankton as SA
    
    real,allocatable,dimension(:)   :: RMAXZ  !< Maximum ration of zooplankton (g prey C/g zooplankton C/d)
    real,allocatable,dimension(:)   :: BMRZ   !< Metabolic rate of zooplankton at reference temperature (1/day)
    real,allocatable,dimension(:)   :: PRRZ   !< Predation rate of zooplankton at reference temperature (1/day)   
    real,allocatable,dimension(:)   :: UBZ    !< Utilization of algaes by zooplankton
    real,allocatable,dimension(:)   :: WQGZ   !< Zooplankton growth rate (1/day)
    real,allocatable,dimension(:)   :: WQBZ   !< Zooplankton metabolic rate (1/day)
    real,allocatable,dimension(:)   :: WQPZ   !< Zooplankton predation rate (1/day)
    real,allocatable,dimension(:)   :: WQDZ   !< Zooplankton death rate (1/day)
  end type  
  
  type(ZOOPLGROUP),allocatable,dimension(:) :: ZOOPL
  
  ! *** Real Array Variables
  real,allocatable,dimension(:,:)   :: PRAZ     !< Prey available to zooplankton (gC/m3) 
  real,allocatable,dimension(:,:,:) :: BAZ      !< Available portion of algal group to zooplankton (gC/m3)
  real,allocatable,dimension(:,:)   :: DOCAZ    !< Available portion of dissolved organic carbon to zooplankton (gC/m3)
  real,allocatable,dimension(:,:)   :: RPOCAZ   !< Available portion of refractory particulate organic carbon to zooplankton (gC/m3)
  real,allocatable,dimension(:,:)   :: LPOCAZ   !< Available portion of labile particulate organic carbon to zooplankton (gC/m3)
  real,allocatable,dimension(:,:)   :: WQTDGZ   !< Temperature dependency of zooplankton grazing
  real,allocatable,dimension(:,:)   :: WQTDBZ   !< Exponential function of temperature for zooplankton metabolism
  real,allocatable,dimension(:,:)   :: WQTDPZ   !< Exponential function of temperature for zooplankton predation
  real,allocatable,dimension(:)     :: WQZKK    !< Related to matrix K1 in the kinetic equations
  real,allocatable,dimension(:)     :: WQZRR    !< Accumulation of sources and sinks in the kinetic equations
  real,allocatable,dimension(:,:,:) :: WQWPSZ   !< Point source loading
  real,allocatable,dimension(:,:)   :: FRLP     
  real,allocatable,dimension(:,:)   :: FRRP     
  real,allocatable,dimension(:,:)   :: FRSI     
  real,allocatable,dimension(:,:,:) :: SBZPAL   !< Total effect of zooplankton on phytoplankton
  real,allocatable,dimension(:,:)   :: SLPOCZ   !< Total effect of zooplanktons on LPOC
  real,allocatable,dimension(:,:)   :: SRPOCZ   !< Total effect of zooplanktons on RPOC
  real,allocatable,dimension(:,:)   :: SDOCZ    !< Total effect of zooplanktons on DOC
    
  real,allocatable,dimension(:,:)   :: SLPONZ   !< Total effect of zooplanktons on LPON
  real,allocatable,dimension(:,:)   :: SRPONZ   !< Total effect of zooplanktons on RPON
  real,allocatable,dimension(:,:)   :: SDONZ    !< Total effect of zooplanktons on DON
  real,allocatable,dimension(:,:)   :: SNH4Z    !< Total effect of zooplanktons on NH4
  real,allocatable,dimension(:,:)   :: SLPOPZ   !< Total effect of zooplanktons on LPOP
  real,allocatable,dimension(:,:)   :: SRPOPZ   !< Total effect of zooplanktons on RPOP
  real,allocatable,dimension(:,:)   :: SDOPZ    !< Total effect of zooplanktons on DOP
  real,allocatable,dimension(:,:)   :: SPO4Z    !< Total effect of zooplanktons on PO4
  real,allocatable,dimension(:,:)   :: SSUZ     !< Total effect of zooplanktons on SU
  real,allocatable,dimension(:,:)   :: SSAZ     !< Total effect of zooplanktons on SA
  real,allocatable,dimension(:,:)   :: SDOZ     !< Total effect of zooplanktons on DO
  
  Contains
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ_Allocate
  !
  !> @details  Allocating and Initializing of WQ array variables
  !---------------------------------------------------------------------------!
  Subroutine WQ_Allocate
    use Allocate_Initialize      
    integer :: NAL
    
    ! *** Call AllocateDSI integer arrays
    call AllocateDSI(IWQCBE,  NBBEM, 0)
    call AllocateDSI(IWQCBN,  NBBNM, 0)
    call AllocateDSI(IWQCBS,  NBBSM, 0)
    call AllocateDSI(IWQCBW,  NBBWM, 0)
    call AllocateDSI(IWQOBE,  NBBEM, NWQVM, 0)
    call AllocateDSI(IWQOBN,  NBBNM, NWQVM, 0)
    call AllocateDSI(IWQOBS,  NBBSM, NWQVM, 0)
    call AllocateDSI(IWQOBW,  NBBWM, NWQVM, 0)
    call AllocateDSI(IWQPSC,  LCM, KCM, 0)
    call AllocateDSI(IWQPSV,  LCM, KCM, 0)
    call AllocateDSI(IWQZMAP, LCM, KCM, 0)

    call AllocateDSI(JWQCBE,  NBBEM, 0)
    call AllocateDSI(JWQCBN,  NBBNM, 0)
    call AllocateDSI(JWQCBS,  NBBSM, 0)
    call AllocateDSI(JWQCBW,  NBBWM, 0)

    call AllocateDSI(ICPSL,   NWQPSM, 0)
    call AllocateDSI(JCPSL,   NWQPSM, 0)
    call AllocateDSI(KCPSL,   NWQPSM, 0)
    call AllocateDSI(MVPSL,   NWQPSM, 0)

    !Call AllocateDSI(NWQCSR, NTSWQVM, 0)  ! *** NWQCSR is used in SCANWQ before this subroutine is called

    ! *** Call AllocateDSI real arrays
    call AllocateDSI(TWQ,       LCM, 0.0)
    call AllocateDSI(SWQ,       LCM, 0.0)
    call AllocateDSI(VOLWQ,     LCM, 0.0)

    call AllocateDSI(WQAPC,     LCM, 0.0)
    call AllocateDSI(WQATM,     MAXWQ, NWQZM, 0.0)  
    call AllocateDSI(WQBFCOD,   LCM, 0.0)
    call AllocateDSI(WQBFNH4,   LCM, 0.0)
    call AllocateDSI(WQBFNO3,   LCM, 0.0)
    call AllocateDSI(WQBFO2,    LCM, 0.0)
    call AllocateDSI(WQBFPO4D,  LCM, 0.0)
    call AllocateDSI(WQBFSAD,   LCM, 0.0)
    call AllocateDSI(WQCHL,     LCM, KCM, 0.0)
    call AllocateDSI(WQDFLC,    LCM, 0.0)
    call AllocateDSI(WQDFLN,    LCM, 0.0)
    call AllocateDSI(WQDFLP,    LCM, 0.0)
    call AllocateDSI(WQDFRC,    LCM, 0.0)
    call AllocateDSI(WQDFRN,    LCM, 0.0)
    call AllocateDSI(WQDFRP,    LCM, 0.0)
    call AllocateDSI(WQDFSI,    LCM, 0.0)
    call AllocateDSI(WQOBCE,    NBBEM, 2, NWQVM, 0.0)
    call AllocateDSI(WQOBCN,    NBBNM, 2, NWQVM, 0.0)
    call AllocateDSI(WQOBCS,    NBBSM, 2, NWQVM, 0.0)
    call AllocateDSI(WQOBCW,    NBBWM, 2, NWQVM, 0.0)
    call AllocateDSI(WQPO4D,    LCM, KCM, 0.0)
    call AllocateDSI(WQSAD,     LCM, KCM, 0.0)
    call AllocateDSI(WQWPSL,    LCM,  KCM, NWQVM, 0.0)
    call AllocateDSI(XSMO20,    LCM, 0.0)

    ! *** Call AllocateDSI algae class variables
    allocate(ALGAES(NALGAEM))
    
    Do NAL = 1, NALGAEM
      ALGAES(NAL).IDN    = NAL
      call AllocateDSI(ALGAES(NAL).WQBMRA,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQDOP,    NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQDRA,    NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKDCALM, NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKHRA,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQPMA,    NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQPRRA,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQWS,     NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQBMIN,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKBP,    NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMV,    NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMVMIN, NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMVA,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMVB,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMVC,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMVD,   NWQZM, 0.0)
      call AllocateDSI(ALGAES(NAL).WQKMVE,   NWQZM, 0.0)

      call AllocateDSI(ALGAES(NAL).SETTLING, LCM, -KCM, 0.0)
    Enddo 
    
    call AllocateDSI(LAYERBOT,       NALGAEM,  LCM,      1)
    call AllocateDSI(LAYERTOP,       NALGAEM,  LCM,      KC)
    call AllocateDSI(HEIGHT_MAC,     LCM,     -NALGAEM,  0.)
    call AllocateDSI(DIAMETER_MAC,   LCM,      NALGAEM,  0.)
    call AllocateDSI(LayerRatio_MAC, LCM,      KCM,     NALGAEM,  0.)
    call AllocateDSI(MACAD,          LCM,      NALGAEM, 0.)
    call AllocateDSI(MVEGZ,          LCM,      NALGAEM, 0.)
    
    call AllocateDSI(IWQKA,    NWQZM,    0)
    call AllocateDSI(IBENMAP,  LCM,   2, 0)
    call AllocateDSI(NMLSER,   NWQPSRM,  0)
    call AllocateDSI(MWQPTLT,  NWQPSRM,  0)
    call AllocateDSI(IMWQZT,   LCM,   0)
    call AllocateDSI(IMWQZT1,  LCM,   0)
    call AllocateDSI(IMWQZT2,  LCM,   0)
    call AllocateDSI(IWQT,     LCM,   0)
    
    call AllocateDSI(WQPA,     LCM,   NALGAEM, 0.0)
    call AllocateDSI(WQBM,     LCM,   NALGAEM, 0.0)
    call AllocateDSI(WQPR,     LCM,   NALGAEM, 0.0)
    call AllocateDSI(WQBSETL,  LCM,   4,  NALGAEM, 0.0)
    call AllocateDSI(WQTDG,    NWQTDM, NALGAEM, 0.0)
    call AllocateDSI(WQTDR,    NWQTDM, NALGAEM, 0.0)
    call AllocateDSI(WQTDP,    NWQTDM, NALGAEM, 0.0)
    call AllocateDSI(WQTDGP,   NWQTDM, 0.0)
    call AllocateDSI(WQTD1FCB, NWQTDM, 0.0)
    call AllocateDSI(WQTD2FCB, NWQTDM, 0.0)
    call AllocateDSI(UHEQ,     LCM, 0.0)
    call AllocateDSI(WQHT,     KCM, 0.0)
    call AllocateDSI(WQWDSL,   LCM, KCM,  NWQVM,  0.0)  
    call AllocateDSI(DDOAVG,   LCM, KCM,  0.0)
    call AllocateDSI(DDOMAX,   LCM, KCM,  0.0)
    call AllocateDSI(DDOMIN,   LCM, KCM,  0.0)
    call AllocateDSI(RLIGHTC,  LCM, KCM,  0.0)
    call AllocateDSI(RLIGHTT,  LCM, KCM,  0.0)
    call AllocateDSI(REAC,     NWQZM,   0.0)
    call AllocateDSI(WQKDC,    NWQZM,   0.0)
    call AllocateDSI(WQKRO,    NWQZM,   0.0)
    call AllocateDSI(WQKTR,    NWQZM,   0.0)
    call AllocateDSI(WQKHCOD,  NWQZM,   0.0)
    call AllocateDSI(WQKCD,    NWQZM,   0.0)
    call AllocateDSI(WQTDHDR,  NWQTDM,  0.0)
    call AllocateDSI(WQTDMNL,  NWQTDM,  0.0)
    call AllocateDSI(WQKSUA,   NWQTDM,  0.0)
    call AllocateDSI(WQKCOD,   NWQTDM,  NWQZM, 0.0)
    call AllocateDSI(WQTDNIT,  NWQTDM,  0.0)
    call AllocateDSI(WQTDKR,   NWQTDM,  NWQZM, 0.0)
    call AllocateDSI(WQTDTAM,  NWQTDM,  0.0)
    call AllocateDSI(WQTAMP,   LCM,   KCM,   0.0)
    call AllocateDSI(WQWSLP,   NWQZM,   0.0)
    call AllocateDSI(WQWSRP,   NWQZM,   0.0)
    call AllocateDSI(WQWSS,    NWQZM,   0.0)
    call AllocateDSI(WQATML,   LCM,   KCM,  NWQVM, 0.0)
    call AllocateDSI(XBENMUD,  LCM,   0.0)
    call AllocateDSI(TSSRD,    NDASER,  0.0)
    call AllocateDSI(SOLFRD,   NDASER,  0.0)
    call AllocateDSI(SOLSRD,   NDASER,  0.0)
    call AllocateDSI(TCWQPSR,  NWQPSRM, 0.0)
    call AllocateDSI(T_MLSER,  NDWQPSR, NWQPSRM, 0.0)
    call AllocateDSI(V_MLSER,  NDWQPSR, NWQVM,  NWQPSRM, 0.0)
    call AllocateDSI(WQI0BOT,  LCM,     0.0)
    call AllocateDSI(DZWQ,     LCM,   0.0)
    call AllocateDSI(WQKRPC,   LCM,   0.0)
    call AllocateDSI(WQKRPN,   LCM,   0.0)
    call AllocateDSI(WQKRPP,   LCM,   0.0)
    call AllocateDSI(WQKLPC,   LCM,   0.0)
    call AllocateDSI(WQKLPN,   LCM,   0.0)
    call AllocateDSI(WQKLPP,   LCM,   0.0)
    call AllocateDSI(WQKHR,    LCM,   0.0)
    call AllocateDSI(WQDENIT,  LCM,   0.0)
    call AllocateDSI(WQPN,     LCM,   NALGAEM, 0.0)
    call AllocateDSI(WQNIT,    LCM,   0.0)    
    call AllocateDSI(O2WQ,     LCM,   0.0)    
    call AllocateDSI(RNH4WQ,   LCM,   0.0)    
    call AllocateDSI(RNO3WQ,   LCM,   0.0) 
    call AllocateDSI(WQDOS,    LCM,   0.0)    
    call AllocateDSI(WQH10,    LCM,   0.0)    
    call AllocateDSI(XDOSAT,   LCM,   KCM, 0.0)    
    call AllocateDSI(WQP19,    LCM,   0.0)    
    call AllocateDSI(WQRREA,   LCM,   0.0)    
    call AllocateDSI(WQKRDOS,  LCM,   0.0)    
    call AllocateDSI(WQKK,     LCM,   0.0)    
    call AllocateDSI(WQRR,     LCM,   0.0)      
    call AllocateDSI(WQRPSET,  LCM,   2,  0.0)    
    call AllocateDSI(WQLPSET,  LCM,   2,  0.0)    
    call AllocateDSI(WQWSSET,  LCM,   2,  0.0)    
    call AllocateDSI(PO4DWQ,   LCM,   0.0)    
    call AllocateDSI(RNH4NO3,  LCM,   0.0)    
    call AllocateDSI(WQOBTOT,  LCM,   0.0)    
    call AllocateDSI(WQKDON,   LCM,   0.0)    
    call AllocateDSI(WQKDOP,   LCM,   0.0)   
    call AllocateDSI(WQT10,    LCM,   0.0)    
    call AllocateDSI(WQT17,    LCM,   0.0)    
    call AllocateDSI(WQN17,    LCM,   0.0)    
    call AllocateDSI(WQO18,    LCM,   0.0)    
    call AllocateDSI(WQR20,    LCM,   0.0)    
    call AllocateDSI(SMAC,     LCM,   1.0)        !< Fixed biota flag set to "on"
    
    ! *** Special Cases
    allocate(ICWQTS(0:NWQVM,NWQTSM))
    allocate(WQV(LCM,KCM,0:NWQVM))
    allocate(WQVO(LCM,0:KCM,0:NWQVM))
    allocate(WQWPSLC(0:NWQPSM,NWQVM))
    allocate(WQPSSRT(NWQVM,0:NWQPSRM))
    ICWQTS  = 0
    WQV     = 0.0
    WQVO    = 0.0
    WQWPSLC = 0.0
    WQPSSRT = 0.0
    
    deallocate(WQKEB)
    call AllocateDSI(WQKEB, NWQZM, 0.0)

    ! *** Dissolved carbon dioxide variables
    call AllocateDSI(CO2WQ,     LCM, 0.0)
    call AllocateDSI(CDOSATIDX, LCM, 0.0)
    call AllocateDSI(WQCDOS,    LCM, 0.0)
    call AllocateDSI(WQITOP,    LCM, KCM, 0.0)
    call AllocateDSI(WQP22,     LCM, 0.0)
    call AllocateDSI(WQKRCDOS,  LCM, 0.0)

    do NAL = 1, NALGAE
      if( IVARSETL == 3 )then
        call AllocateDSI(CELLDENS,      LCM, KCM, NALGAE, 0.0)
        call AllocateDSI(CELLDENSLIGHT, LCM, KCM, NALGAE, 0.0)
        EXIT
      endif
    enddo
    
  End Subroutine WQ_Allocate

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SD_Allocate
  !
  !> @details  Allocating and Initializing of SD array variables
  !---------------------------------------------------------------------------!

  Subroutine SD_Allocate
    use Allocate_Initialize 
    
    ! *** Call AllocateDSI integer arrays
    call AllocateDSI(ISMT,     LCM, 0)
    call AllocateDSI(ISMZMAP,  LCM, 0)
    
    ! *** Call AllocateDSI real arrays
    call AllocateDSI(WQDFB,     LCM,  NALGAEM, 0.0)
    call AllocateDSI(SMD1PO4,   LCM,  0.0)
    call AllocateDSI(SMD1SI,    LCM,  0.0)
    call AllocateDSI(SMDD,      NSMZM,  0.0)
    call AllocateDSI(SMDFSI,    LCM,  0.0)
    call AllocateDSI(SMDGFC,    LCM,  0.0)
    call AllocateDSI(SMDGFN,    LCM,  0.0)
    call AllocateDSI(SMDGFP,    LCM,  0.0)
    call AllocateDSI(SMDP,      NSMZM,  0.0)
    call AllocateDSI(SMDP1PO4,  NSMZM,  0.0)
    call AllocateDSI(SMDPMIN,   NSMZM,  0.0)
    call AllocateDSI(SMDTOH,    NSMZM,  0.0)
    call AllocateDSI(SMFCBA,    NALGAEM, NSMGM, 0.0)
    call AllocateDSI(SMFNBA,    NALGAEM, NSMGM, 0.0)
    call AllocateDSI(SMFPBA,    NALGAEM, NSMGM, 0.0)
    call AllocateDSI(SMFCR,     NSMZM,  NSMGM, 0.0)
    call AllocateDSI(SMFNR,     NSMZM,  NSMGM, 0.0)
    call AllocateDSI(SMFPR,     NSMZM,  NSMGM, 0.0)
    call AllocateDSI(SMHODT,    NSMZM,  0.0)
    call AllocateDSI(SMHYPD,    LCM,  0.0)
    call AllocateDSI(SMJAQH2S,  LCM,  0.0)
    call AllocateDSI(SMJDEN,    LCM,  0.0)
    call AllocateDSI(SMJGCH4,   LCM,  0.0)
    call AllocateDSI(SMJNIT,    LCM,  0.0)
    call AllocateDSI(SMK1H2S,   NWQTDM, 0.0)
    call AllocateDSI(SMK1NO3,   NSMZM,  0.0)
    call AllocateDSI(SMK2NO3,   NSMZM,  0.0)
    call AllocateDSI(SMKL12,    LCM,    0.0)
    call AllocateDSI(SMKNH4,    NSMZM,  0.0)
    call AllocateDSI(SMSS,      LCM,  0.0)
    call AllocateDSI(SMTD1CH4,  NWQTDM, 0.0)
    call AllocateDSI(SMTD2CH4,  NWQTDM, 0.0)
    call AllocateDSI(SMTDCD,    NWQTDM, NSMGM, 0.0)
    call AllocateDSI(SMTDDD,    NWQTDM, 0.0)
    call AllocateDSI(SMTDDP,    NWQTDM, 0.0)
    call AllocateDSI(SMTDND,    NWQTDM, NSMGM, 0.0)
    call AllocateDSI(SMTDNH4,   NWQTDM, 0.0)
    call AllocateDSI(SMTDNO3,   NWQTDM, 0.0)
    call AllocateDSI(SMTDPD,    NWQTDM, NSMGM, 0.0)  
    call AllocateDSI(SMTDSI,    NWQTDM, 0.0) 
    call AllocateDSI(SMTMP,     LCM,    0.0)   
    call AllocateDSI(SMW12,     LCM,    0.0)
    call AllocateDSI(SMW2,      NSMZM,  0.0)
    call AllocateDSI(SMW2DTOH,  NSMZM,  0.0)
    call AllocateDSI(SMW2PHODT, NSMZM,  0.0)
    call AllocateDSI(SODMULT,   NSMZM,  0.0)
    call AllocateDSI(SM1DIFT,   NSMZM,  0.0)
    call AllocateDSI(SM2DIFT,   NSMZM,  0.0)
    call AllocateDSI(SM2NO3,    LCM,  0.0)
    call AllocateDSI(SM2PO4,    LCM,  0.0)
    call AllocateDSI(SMHSED,    NSMZM,  0.0)
    call AllocateDSI(SMPOP,     LCM,  NSMGM, 0.0) 
    call AllocateDSI(SM1H2S,    LCM,  0.0)
    call AllocateDSI(SM1NH4,    LCM,  0.0)
    call AllocateDSI(SM1NO3,    LCM,  0.0)
    call AllocateDSI(SM1PO4,    LCM,  0.0)
    call AllocateDSI(SM1SI,     LCM,  0.0)
    call AllocateDSI(SM2H2S,    LCM,  0.0)
    call AllocateDSI(SM2NH4,    LCM,  0.0)
    call AllocateDSI(SM2SI,     LCM,  0.0)
    call AllocateDSI(SMBST,     LCM,  0.0)
    call AllocateDSI(SMCSOD,    LCM,  0.0)
    call AllocateDSI(SMDFC,     LCM,  NSMGM, 0.0)
    call AllocateDSI(SMDFN,     LCM,  NSMGM, 0.0)
    call AllocateDSI(SMDFP,     LCM,  NSMGM, 0.0)
    call AllocateDSI(SMPOC,     LCM,  NSMGM, 0.0)
    call AllocateDSI(SMPON,     LCM,  NSMGM, 0.0)
    call AllocateDSI(SMNSOD,    LCM,  0.0)
    call AllocateDSI(SMPSI,     LCM,  0.0)
    call AllocateDSI(SMT,       LCM,  0.0)
       
    ! *** Special cases    
    allocate(SMHYST(LCM))
     
  End Subroutine SD_Allocate
    
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ3DINP
  !
  !> @details  Allocating and Initializing of SD array variables
  !---------------------------------------------------------------------------!

  Subroutine Zoo_Allocate
    use Allocate_Initialize      
    integer :: NZO
        
    ! *** Call AllocateDSI integer arrays
    call AllocateDSI(IWQZT,     LCM, 0)
    
    ! *** Call AllocateDSI algae class variables
    allocate(ZOOPL(NZOOPL))
    
    Do NZO = 1, NZOOPL
      call AllocateDSI(ZOOPL(NZO).RMAXZ, NWQZM,  0.0)
      call AllocateDSI(ZOOPL(NZO).BMRZ,  NWQZM,  0.0)
      call AllocateDSI(ZOOPL(NZO).PRRZ,  NWQZM,  0.0)
      call AllocateDSI(ZOOPL(NZO).UBZ,   NALGAEM, 0.0)
      call AllocateDSI(ZOOPL(NZO).WQGZ,  LCM,  0.0)
      call AllocateDSI(ZOOPL(NZO).WQBZ,  LCM,  0.0)
      call AllocateDSI(ZOOPL(NZO).WQPZ,  LCM,  0.0)
      call AllocateDSI(ZOOPL(NZO).WQDZ,  LCM,  0.0)
    Enddo 
     
    call AllocateDSI(PRAZ,    LCM,  NZOOPL, 0.0)
    call AllocateDSI(BAZ,     LCM,  NZOOPL, NALGAEM, 0.0)
    call AllocateDSI(DOCAZ,   LCM,  NZOOPL, 0.0)
    call AllocateDSI(RPOCAZ,  LCM,  NZOOPL, 0.0)
    call AllocateDSI(LPOCAZ,  LCM,  NZOOPL, 0.0)
    call AllocateDSI(WQTDGZ,  NWQTDM, NZOOPL, 0.0)
    call AllocateDSI(WQTDBZ,  NWQTDM, NZOOPL, 0.0)
    call AllocateDSI(WQTDPZ,  NWQTDM, NZOOPL, 0.0)
    call AllocateDSI(WQZKK,   LCM,  0.0)
    call AllocateDSI(WQZRR,   LCM,  0.0)
    call AllocateDSI(WQWPSZ,  LCM,  KCM,    NZOOPL, 0.0)
    call AllocateDSI(FRRP,    LCM,  NZOOPL, 0.0)
    call AllocateDSI(FRLP,    LCM,  NZOOPL, 0.0)
    call AllocateDSI(FRSI,    LCM,  NZOOPL, 0.0)
    call AllocateDSI(SBZPAL,  LCM,  KCM, NALGAEM, 0.0)
    call AllocateDSI(SLPOCZ,  LCM,  KCM, 0.0)
    call AllocateDSI(SRPOCZ,  LCM,  KCM, 0.0)
    call AllocateDSI(SDOCZ,   LCM,  KCM, 0.0)
    call AllocateDSI(SLPONZ,  LCM,  KCM, 0.0)
    call AllocateDSI(SRPONZ,  LCM,  KCM, 0.0)
    call AllocateDSI(SDONZ,   LCM,  KCM, 0.0)
    call AllocateDSI(SNH4Z,   LCM,  KCM, 0.0)
    call AllocateDSI(SLPOPZ,  LCM,  KCM, 0.0)
    call AllocateDSI(SRPOPZ,  LCM,  KCM, 0.0)
    call AllocateDSI(SDOPZ,   LCM,  KCM, 0.0)
    call AllocateDSI(SPO4Z,   LCM,  KCM, 0.0)
    call AllocateDSI(SSUZ,    LCM,  KCM, 0.0)
    call AllocateDSI(SSAZ,    LCM,  KCM, 0.0)
    call AllocateDSI(SDOZ,    LCM,  KCM, 0.0)
    
  End Subroutine Zoo_Allocate
  
  
End Module Variables_WQ
