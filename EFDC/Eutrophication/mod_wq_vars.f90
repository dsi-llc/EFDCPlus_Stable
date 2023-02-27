! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
Module Variables_WQ
  Use GLOBAL
  implicit none

  Character(len=5) :: WQCONSTIT(50)
  Character(len=5) :: ZOONAME(20)
  
  ! *** Variables for new Water Quality features
  INTEGER, PARAMETER :: MAXWQ = 50 !< Maximum number of WQ components
  INTEGER, PARAMETER :: NTSWQVM = 50    !< Maximum number of WQ time series
  INTEGER :: NALGAE          !< Number of algal groups
  INTEGER :: NALGAEM         !< Maximum of NALGAE or 1 for allocation of arrays
  INTEGER :: ALG_COUNT       !< Number of algae available
  INTEGER :: IWQZPL          !< Switch to Zooplankton module
  INTEGER :: NZOOPL          !< Number of zooplankton groups modeling
  INTEGER :: NWQVZ           !< Number of WQ components excluded zooplankton
  INTEGER :: ISMOB(20)       !< Temporaly mobile flag
  
  ! *** list of index for WQ components
  INTEGER :: IROC, ILOC, IDOC, IROP, ILOP, IDOP, IP4D, IRON, ILON, IDON, INHX,  &
             INOX, ISUU, ISAA, ICOD, IDOX, ITAM, IFCB, ICO2, ICHC, ICHD, ICHG
  
  ! *** Water Column Eutrophication variables
  
  ! *** Integer Variables
  INTEGER :: NWQV            !< Number of water quality variables, including active and inactive consitutuents
  INTEGER :: NWQVM           !< Maximum number of water quality variables, including active and inactive consitutuents
  INTEGER :: IWQZONES = 0    !< Flag to activact the use of WQ zones
  INTEGER :: IWQBEN = 0      !< Sediment diagensis flux option

  INTEGER :: DIATOM          !< Index of diatom when switching to silica
  INTEGER :: ISWQLVL         !< Global flag to active macrophyte and/or priphyton velocity limitation
  INTEGER :: ITNWQ           !< Water quality kinetics counter
  INTEGER :: IWQICI
  INTEGER :: IWQM
  INTEGER :: IWQNC
  INTEGER :: IWQNPL
  INTEGER :: IWQONC
  INTEGER :: IWQORST
  INTEGER :: IWQPSL
  INTEGER :: IWQRST
  INTEGER :: IWQSI
  INTEGER :: IWQSRP
  INTEGER :: IWQSTOX
  INTEGER :: IWQSUN
  INTEGER :: IWQTS
  INTEGER :: IWQTSB
  INTEGER :: IWQTSDT
  INTEGER :: IWQTSE


  INTEGER :: NDWQPSR
  INTEGER :: NFIXED          !< Number of fixed biota classes.  Replaces old IDNOTRVA
  INTEGER :: NTSWQV          !< Number of WQ time series
  INTEGER :: NWQCSRM
  INTEGER :: NWQOBE
  INTEGER :: NWQOBN
  INTEGER :: NWQOBS
  INTEGER :: NWQOBW
  INTEGER :: NWQPS
  INTEGER :: NWQPSM          !< Number of cells for point source mass loading
  INTEGER :: NWQPSRM         !< Maximum number of cells for point source mass loading
  INTEGER :: NWQTD           !< Number of temperature lookup table points
  INTEGER :: NWQTDM          !< Maximum number of temperature lookup table points
  INTEGER :: NWQTS           !< Number of WQ time series output locations (not used by EEMS)
  INTEGER :: NWQTSM          !< Maximum number of WQ time series output locations (not used by EEMS)
  INTEGER :: NWQZ            !< Number of water quality zones
  INTEGER :: NWQZM           !< Maximum number of water quality zones
  INTEGER :: IDOSFRM         !< Formulation for D.O. saturation: 0=Garcia and Gordon (1992), 1=Chapra et al. (1997), 2=Genet et al. (1974)  
  INTEGER :: IDOSELE         !< Elevation adjustment for D.O. saturation flag: 0=Not used, 1=Chapra et al. (1997), 2=Zison et al. (1978)

  INTEGER :: ISKINETICS(MAXWQ)    !< Eutrophication parameter activation flag
  INTEGER :: ISTRWQ(MAXWQ)        !< Eutrophication parameter transport flag
  
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IWQCBE        !< East  open BC WQ list: I Index 
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IWQCBN        !< North open BC WQ list: I Index
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IWQCBS        !< South open BC WQ list: I Index
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IWQCBW        !< West  open BC WQ list: I Index
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQOBE        !< East  open BC concentration series index
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQOBN        !< North open BC concentration series index
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQOBS        !< South open BC concentration series index
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQOBW        !< West  open BC concentration series index
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQPSC        !< 
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQPSV        !< 
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IWQZMAP       !< WQ zone loopup by L and K

  INTEGER,ALLOCATABLE,DIMENSION(:)     :: JWQCBE      !< East  open BC WQ list: J Index
  INTEGER,ALLOCATABLE,DIMENSION(:)     :: JWQCBN      !< North open BC WQ list: J Index
  INTEGER,ALLOCATABLE,DIMENSION(:)     :: JWQCBS      !< South open BC WQ list: J Index
  INTEGER,ALLOCATABLE,DIMENSION(:)     :: JWQCBW      !< West  open BC WQ list: J Index

  INTEGER  :: NWQCSR(NTSWQVM)

  ! *** Real Variables
  REAL :: DTWQ
  REAL :: DTWQO2
  REAL :: DTWQZO2
  REAL :: DOELEV    = 0.0
  REAL :: WQAOCR    = 0.0
  REAL :: WQCIA     = 0.0
  REAL :: WQCIB     = 0.0
  REAL :: WQCIC     = 0.0
  REAL :: WQCIM     = 0.0
  REAL :: WQI0      = 0.0
  REAL :: WQI1      = 0.0
  REAL :: WQI2      = 0.0
  REAL :: WQI3      = 0.0
  REAL :: WQKECHL   = 0.0
  REAL :: WQKECHLE  = 0.0
  REAL :: WQKEPOM   = 0.0
  REAL :: WQKEDOM   = 0.0
  REAL :: WQKEPOC   = 0.0
  REAL :: WQKETSS   = 0.0
  REAL :: WQKINUPT  = 0.0
  REAL :: WQTDMIN   = 0.0
  REAL :: WQTDMAX   = 0.0
  REAL :: WQTDINC   = 0.0
  REAL :: WQTSDT    = 0.0
  REAL :: WQTRHDR   = 0.0
  REAL :: WQKTHDR   = 0.0
  REAL :: WQTRMNL   = 0.0
  REAL :: WQKTMNL   = 0.0
  REAL :: WQTNIT    = 0.0
  REAL :: WQKN1    = 0.0
  REAL :: WQKN2    = 0.0
  REAL :: WQKSU    = 0.0
  REAL :: WQTRSUA    = 0.0
  REAL :: WQKTSUA    = 0.0
  REAL :: WQBFTAM    = 0.0
  REAL :: WQTTAM    = 0.0
  REAL :: WQKTAM    = 0.0
  REAL :: WQTRCOD    = 0.0
  REAL :: WQKTCOD    = 0.0
  
  REAL,ALLOCATABLE,DIMENSION(:)     :: SWQ           !< Temporary variable for salinity for current layer
  REAL,ALLOCATABLE,DIMENSION(:)     :: TWQ           !< Temporary variable for water temperature for current layer
  REAL,ALLOCATABLE,DIMENSION(:)     :: XSMO20
  REAL,ALLOCATABLE,DIMENSION(:)     :: VOLWQ         !< Temporary variable for inverse of layer volume (1/m3)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQAPC         !< Algal Phosphorus-to-Carbon ratio
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQATM         !< Wet Depsition by WQ zone and WQ components
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WQBFCOD    !< COD Sediment flux 
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WQBFNH4    !< NH4 Sediment flux
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WQBFNO3    !< NO3 Sediment flux
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WQBFO2     !< SOD Sediment flux
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WQBFPO4D   !< PO4 Sediment flux
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WQBFSAD    !< Silica Sediment flux
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQCHL         !< Total clhlorphyll in the water column
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFLC        !< Depositional flux of LPOC 
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFLN        !< Depositional flux of LPON
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFLP        !< Depositional flux of LPOP
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFRC        !< Depositional flux of RPOC 
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFRN        !< Depositional flux of RPON
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFRP        !< Depositional flux of RPOP
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDFSI        !< Depositional flux of SI
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKEB         !< Background extinction coefficient for WQ calculations
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQOBCE        !< Open BC concentration series
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQOBCN        !< Open BC concentration series
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQOBCS        !< Open BC concentration series
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQOBCW        !< Open BC concentration series
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQPO4D        !< Bioavailable (i.e. dissolved) phase to total PO4
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQSAD         !< Bioavailable (i.e. dissolved) phase to total SI
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTDTEMP      !< Uniform temperature lookup table for all temperature dependent WQ processes
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQWPSL        !< Mass loading of WQ parameters for point sources (only IWQPSL /= 2)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQWPSLC       !< Concentrations of WQ parameters for point sources


  ! *** Algae/Macrophyte variables
  TYPE ALGAECLASS
    CHARACTER(LEN = 20) :: PHYTONAME
    INTEGER :: IDN  !< Algal ID group number, 1 - cyanobacteria, 2 - algae diatoms, 3 - greens algae, 4 - macroalgae
    INTEGER :: ISTOX                                 !< Flag for salinity toxicity (use for Cyanobacteria)
    INTEGER :: ISBLOOM                               !< Flag for bloom in winter (use for Diatoms)
    INTEGER :: ISILICA                               !< Flag for silica active (use for Diatoms)
    INTEGER :: IWQVLIM                               !< Velocity limitation options (use for Macrophytes)
    
    LOGICAL :: ISMOBILE                              !< Flag to indicate mobility

    REAL,ALLOCATABLE,DIMENSION(:)     :: WQBMRA      !< Basal metabolism rate for algal group (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQDOP       !< Optimal depth (m) for algal group growth
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQDRA       !< Maximum death rate for algal croup (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKDCALM    !< Constant relating DOC hydrolysis rate to macroalgae
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKHRA      !< Half-sat. constant (gO2/m^3) for algal group DOC excretion
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQPMA       !< Maximum growth rate for algal group (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQPRRA      !< Maximum predation rate on algal group (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKBP       !< Growth limitation - Plant density half saturation factor   (G C/M2)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMV       !< Velocity limitation - half saturation velocity               (m/s)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMVMIN    !< Velocity limitation - Minimum velocity to apply limiation    (m/s)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMVA      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMVB      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMVC      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMVD      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQKMVE      !< Velocity limitation - A     Logistic funciton: Vel_limit = D + (A - D) / (1 + (Vel/C)**B)**E     
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQWS        !< Settling velocity for algal group (m/day)
    REAL,ALLOCATABLE,DIMENSION(:)     :: WQBMIN      !< Minimum biomass for macrophyte    Input as g C/m2 but converted to g C/m3 using IC water depths

    REAL :: WQBMAX         !< Maximum biomass for macrophyte    (g C/m3)
    REAL :: WQTM1          !< Lower optimal temperature for algal group growth (degC)
    REAL :: WQTM2          !< Upper optimal temperature for algal group growth (degC)
    REAL :: WQKG1          !< Suboptimal temperature effect coef. for algal group growth
    REAL :: WQKG2          !< Superoptimal temperature effect coef. for algal group growth
    REAL :: WQTR           !< Reference temperature for algal group metabolism (degC)
    REAL :: WQKTB          !< Temperature effect coef. for algal group metabolism
    REAL :: WQTP1          !< Lower optimal temperature for algal group predation (degC)
    REAL :: WQTP2          !< Upper optimal temperature for algal group predation (degC)
    REAL :: WQKP1          !< Suboptimal temperature effect coef. for algal group predation
    REAL :: WQKP2          !< Superoptimal temperature effect coef. for algal group predation
    REAL :: WQKHNA         !< Nitrogen half-saturation for algal group        (mg/L)
    REAL :: WQKHPA         !< Phosphorus half-saturation for algal group      (mg/L)   
    REAL :: WQKHCO2        !< CO2 half-saturation consts for algal group
                           
    REAL :: WQCHLA         !< Carbon-to-Chlorophyll ratio for algal group (mg C/ug Chl)
    REAL :: WQFCRP         !< Fraction of predated Carbon produced as RPOC
    REAL :: WQFCLP         !< Fraction of predated Carbon produced as LPOC
    REAL :: WQFCDP         !< Fraction of predated Carbon produced as DOC
                           
    REAL :: WQFCDB         !< Fraction of basal metabolism exuded as DOC   
                           
    REAL :: WQFPRP         !< Fraction of predated Phosphorus produced as RPOP
    REAL :: WQFPLP         !< Fraction of predated Phosphorus produced as LPOP
    REAL :: WQFPDP         !< Fraction of predated Phosphorus produced as DOP
    REAL :: WQFPIP         !< Fraction of predated Phosphorus produced as Inorganic P
                           
    REAL :: WQFPRB         !< Fraction of metabolized Phosphorus produced as RPOP
    REAL :: WQFPLB         !< Fraction of metabolized Phosphorus produced as LPOP
    REAL :: WQFPDB         !< Fraction of metabolized Phosphorus produced as DOP
    REAL :: WQFPIB         !< Fraction of metabolized Phosphorus produced as P4T
                           
    REAL :: WQFNRP         !< Fraction of predated Nitrogen produced as RPON
    REAL :: WQFNLP         !< Fraction of predated Nitrogen produced as LPON
    REAL :: WQFNDP         !< Fraction of predated Nitrogen produced as DON
    REAL :: WQFNIP         !< Fraction of predated Nitrogen produced as Inorganic N
                           
    REAL :: WQFNRB         !< Fraction of metabolized Nitrogen produced as RPON
    REAL :: WQFNLB         !< Fraction of metabolized Nitrogen produced as LPON
    REAL :: WQFNDB         !< Fraction of metabolized Nitrogen produced as DON
    REAL :: WQFNIB         !< Fraction of metabolized Nitrogen produced as Inorganic N
                           
    REAL :: WQANCA         !< Nitrogen-to-carbon ratio for algal group (gN/gC)
                           
    REAL :: WQFSID         !< Fraction of metabolized Silica by algae (diatoms) produced as available Silica
    REAL :: WQFSPD         !< Fraction of metabolized Silica by algae (diatoms) produced as particulate biogenic Silica
    REAL :: WQFSIP         !< Fraction of predated Silica by algae (diatoms) produced as available Silica
    REAL :: WQFSPP         !< Fraction of predated Silica by algae (diatoms) produced as particulate biogenic Silica
                           
    REAL :: WQASC          !< Silica-to-carbon ratio for algae (diatoms)
    REAL :: WQKHS          !< Silica half-saturation for algae (diatoms) (mg/L)
    
    REAL :: WQC2DW         !< Carbon to Dry Weight Ratio
    REAL :: WQKEMAC        !< Light extinction coefficient due to macrophyte biomass   (m3/gC/m)
                           
    REAL :: WQSTOX         !< Salinity at which microcystis growth is halved (cyanobacteria)
    REAL :: WQAOCRP        !< Algae photosynthesis oxygen-to-carbon ratio (Macroalgae)
    REAL :: WQAOCRR        !< Algae respiration oxygen-to-carbon ratio (Macroalgae)
                           
    REAL :: WQALGOCR       !< Stoichiometric algae oxygen-to-carbon ratio (gO2/gC)
    REAL :: WQALGONT       !< Stoichiometric algae oxygen=to-nitrate ratio (gO2/gN)
    REAL :: WQAPCM         !< Factor to modify APC for macroalgae from 
     
    REAL :: THRESHOLD      !< Concentration factor to convert macrophyte concentration to plant height (g C/m2/(1 meter hieght (m)))
    REAL :: BASEDEPTH      !< Distance below the surface (i.e. depth) of the "base" of the macrophyte/periphyton growth (m)
    REAL :: MAXLENGTH      !< Maximum length from the "base" to allow macrophyte/periphyton growth (m)
    
    ! *** Hydrodynamic Feedback
    INTEGER :: ISDRAG      !< Hydrodynamic feedback flag: 0 - None,  1 - Feedback activated
    REAL :: DRAGCOEFF      !< Hydrodynamic drag coefficient for macrophyte stem and leaves
    REAL :: MINDIAMETER    !< Minimum diameter of each plant (m)
    REAL :: STEMHEIGHT     !< Total stem height of each plant (m).  Changes with growth
    REAL :: STEMDENSITY    !< Plant/stem density (#/m2)
    
  ENDTYPE
  INTEGER :: MACDRAG = 0
  TYPE(ALGAECLASS),ALLOCATABLE,DIMENSION(:)    :: ALGAES     !< ALGAL GROUPS
  
  REAL :: WQKRC       !< Minimum hydrolysis rate (1/day) of RPOC
  REAL :: WQKLC       !< Minimum hydrolysis rate (1/day) of LPOC
  REAL :: WQKRCALG    !< Constant relating RPOC hydrolysis rate to total Chl-
  REAL :: WQKLCALG    !< Constant relating LPOC hydrolysis rate to total Chl-
  REAL :: WQKDCALG    !< Constant relating DOC hydrolysis rate to total Chl-
  REAL :: WQAANOX     !< Ratio of denitrification rate to oxic DOC respiration rate
  
  REAL :: WQKRP       !< Minimum hydrolysis rate (1/day) of RPOP
  REAL :: WQKLP       !< Minimum hydrolysis rate (1/day) of LPOP
  REAL :: WQKDP       !< Minimum hydrolysis rate (1/day) of DOP
  REAL :: WQKRPALG    !< Constant relating hydrolysis rate of RPOP to algae
  REAL :: WQKLPALG    !< Constant relating hydrolysis rate of LPOP to algae
  REAL :: WQKDPALG    !< Constant relating hydrolysis rate of DOP to algae
  REAL :: WQCP1PRM    !< Constant used in determining algae Phos-to-Carbon ratio
  REAL :: WQCP2PRM    !< Constant used in determining algae Phos-to-Carbon ratio
  REAL :: WQCP3PRM    !< Constant used in determining algae Phos-to-Carbon ratio
  REAL :: WQKPO4P     !< Partition coefficient for sorbed/dissolved PO4
  
  REAL :: WQKRN       !< Minimum hydrolysis rate (1/day) of RPON
  REAL :: WQKLN       !< Minimum hydrolysis rate (1/day) of LPON
  REAL :: WQKDN       !< Minimum hydrolysis rate (1/day) of DON
  REAL :: WQKRNALG    !< Constant relating hydrolysis rate of RPON to algae
  REAL :: WQKLNALG    !< Constant relating hydrolysis rate of LPON to algae
  REAL :: WQKDNALG    !< Constant relating hydrolysis rate of DON to algae
  REAL :: WQKHNDO     !< Nitrification half-sat. constant for D.O.
  REAL :: WQKHNN      !< Nitrification half-sat. constant for NH4
  REAL :: WQANDC      !< Mass NO3 reduces per DOC oxidized (gN/gC)
  REAL :: WQNITM      !< Maximum nitrification rate (/day)
  REAL :: RNH4WQ_     !< Ammonia (for Current Layer)
  REAL :: RNO3WQ_     !< Nitrate (for Current Layer)
  
  REAL :: WQKSAP      !< Partition coef. for sorbed/dissolved SA
  REAL :: WQKHORDO    !< Oxic respiration half-sat. constant for D.O. (gO2/m^3)
  REAL :: WQKHDNN     !< Half-sat. constant for denitrification (gN/m^3)
  REAL :: WQAONT      !< Stoichiometric algae oxygen=to-nitrate ratio (gO2/gN)
  REAL :: WQFD        !< Solar radiation fraction daylength
  REAL :: WQISMIN     !< Minimum optimum solar radiation (Langley/day)
  REAL :: WQI0OPT     !< Optimal solar radiation
  REAL :: SOLSRDT     !< Daily average solar radiation interpolation for water quality
  REAL :: SOLFRDT     
  REAL :: WQKHBMF     !< D.O. concentration where TAM release is half the anoxic rate
  REAL :: WQTAMDMX    !< TAM solubility at anoxic conditions (mol/m^3)
  REAL :: WQKDOTAM    !< Constant relating TAM solubility to D.O.
  REAL :: WQKFCB      !< First-order fecal coliform bacteria decay rate (1/day)
  REAL :: WQTFCB      !< Temperature effect constant for KFCB decay rate
  REAL :: STEMFAC     !< SOD temperature factor, use 1.0 to ignore temperature effects  
  REAL :: WQTSB
  REAL :: WQTSE
  
  REAL :: WQHRAVG

  ! *** Integer Array Variables
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IWQKA
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ICWQTS
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IBENMAP
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: MWQPSR
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: MWQPTLT
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IMWQZT
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IMWQZT1
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IMWQZT2
  INTEGER,ALLOCATABLE,DIMENSION(:)   :: IWQT      !< Index for look-up table for temperature dependency
  
  ! *** Real Array Variables  
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQPA       !< Net Growth rate for algal group (1/day)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQBM       !< Basal metabolism rate for algal group (1/day)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQPR       !< Predation rate on algal group (1/day)
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQBSET     !< Settling rate of a general algal group 
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDG      !< Temperature dependency of algal growth
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDR      !< Algal metabolism/predation exponential function of a general algal group
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDP      !< Algal predation exponential function of a general algal group - WQSKE2
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTDGP     !< Algal predation exponential function of a general algal group - WQSKE1
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTD1FCB   
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTD2FCB   
  REAL,ALLOCATABLE,DIMENSION(:)     :: UHEQ       
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQHT       !< Fractional of Depth at the Top of the Layer
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQWDSL     
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: DDOAVG     !< Diurnal domain analysis
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: DDOMAX     !< Diurnal domain analysis
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: DDOMIN     !< Diurnal domain analysis
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: RLIGHTC    !< Light extinction analysis
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: RLIGHTT    !< Light extinction analysis
  REAL,ALLOCATABLE,DIMENSION(:)     :: REAC       !< Global reaeration adjustment factor
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKDC      !< Minimum hydrolysis rate (1/day) of DOC
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKRO      !< Reaeration constant (3.933 for OConnor-Dobbins; 5.32 for Owen-Gibbs)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKTR      !< Temperature rate constant for reaeration
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKHCOD    !< Oxygen half-saturation constant for COD decay (mg/L O2)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKCD      !< COD decay rate (per day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTDHDR    !< Temperature effect on hydrolysis
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTDMNL    !< Temperature effect on mineralization
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKSUA     !< Temperature effect on PSI dissolution
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQKCOD     !< Temperature effect on COD Oxidation    
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTDNIT    
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDKR     !< Temperature effect Reaeration rate coefficient
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQTDTAM    
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTAMP     
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQWSLP     !< Settling velocity for refractory POM (m/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQWSRP     !< Settling velocity for refractory POM (m/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQWSS      !< Settling velocity for particles sorbed to TAM(m/day)
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQATML     
  REAL,ALLOCATABLE,DIMENSION(:)     :: XBENMUD    
  REAL,ALLOCATABLE,DIMENSION(:)     :: TSSRD      
  REAL,ALLOCATABLE,DIMENSION(:)     :: SOLFRD     
  REAL,ALLOCATABLE,DIMENSION(:)     :: SOLSRD     
  REAL,ALLOCATABLE,DIMENSION(:)     :: TCWQPSR    
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: TWQPSER    
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQPSSER    
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQPSSRT    
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQI0BOT    
  REAL,ALLOCATABLE,DIMENSION(:)     :: DZWQ       !< Inverse of cell layer height
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKRPC     !< Hydrolysis rate of RPOC (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKRPN     !< Hydrolysisrate of RPON (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKRPP     !< Hydrolysisrate of RPOP (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKLPC     !< Hydrolysis rate of LPOC (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKLPN     !< Hydrolysisrate of LPON (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKLPP     !< Hydrolysisrate of LPOP (1/day)  
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKHR      !< Heterotrophic respiration rate of DOC (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDENIT    !< Denitrification rate of DOC (1/day)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQPN       !< Preference for ammonium uptake by algaes
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQNIT      !< Nitrification rate (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: O2WQ       !< Dissolved oxygen concentration (mg/l)
  REAL,ALLOCATABLE,DIMENSION(:)     :: RNH4WQ     
  REAL,ALLOCATABLE,DIMENSION(:)     :: RNO3WQ     
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQDOS      !< Saturation concentration of dissolved oxygen (gO2/m3)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQH10      
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: XDOSAT     
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQP19      
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKRDOS    !< Kr * DOs
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKK       !< Related to matrix K1
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQRR       !< Accumulation of sources and sinks in the kinetic equation  
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQRPSET    !< Settling rate of ROP
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQLPSET    !< Settling rate of LOP
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQWSSET    
  REAL,ALLOCATABLE,DIMENSION(:)     :: PO4DWQ     !< Phosphate (for Current Layer)
  REAL,ALLOCATABLE,DIMENSION(:)     :: RNH4NO3    !< Total Inorganic Nitrogen (for Current Layer)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQOBTOT    !< Total Algal Biomass (mg/l)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKDON     !< Mineralization rate of DON (1/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKDOP     !< Mineralization rate of DOP (1/day)  
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQT10      
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQT17      
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQN17      
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQO18      !< Matrix K1 in COD kinetic equation
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQR20      
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMAC       !< Flag to turn on (1.0) or off (0.0) fixed biota 
    
  ! *** Begin Dissolved Carbondioxide variables    VB
  REAL,ALLOCATABLE,DIMENSION(:)     :: CO2WQ
  REAL,ALLOCATABLE,DIMENSION(:)     :: CDOSATIDX
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQCDOS
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQITOP
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQKRCDOS
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQP22
  ! End Dissolved Carbondioxide variables  
  
  ! *** Macrophyte/Preiphyton variables
  LOGICAL,ALLOCATABLE,DIMENSION(:)   :: MAC_CELL  !< Flag indicating cell has macrophyte growth with hydrodynamic feedback
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: LAYERBOT  !< Lowest/Bottom layer for macrophyte/periphyton growth
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: LAYERTOP  !< Highest/Top   layer for macrophyte/periphyton growth

  REAL,TARGET,ALLOCATABLE,DIMENSION(:,:)   :: HEIGHT_MAC      !< Current macrophyte height above the starting growth height (m)
  REAL,TARGET,ALLOCATABLE,DIMENSION(:,:)   :: DIAMETER_MAC    !< Current macrophyte stem diameter (m)  
  REAL,TARGET,ALLOCATABLE,DIMENSION(:,:,:) :: LayerRatio_MAC  !< Fraction of vegetation vertical extents in the layer
  
  ! *** Sediment Diagenesis variables 
  INTEGER :: ISMHYST
  INTEGER :: ISMICI
  INTEGER :: ISMORST
  INTEGER :: ISMRST
  INTEGER :: ISMTS
  INTEGER :: ISMTSB
  INTEGER :: ISMTSDT
  INTEGER :: ISMTSE
  INTEGER :: ISMZ
  INTEGER :: ISMZB

  INTEGER :: NSMG    !< Number of reactive glasses.  Hardwired to 3: G1, G2 and G3
  INTEGER :: NSMGM   !< Maximum NSMG
  INTEGER :: NSMZ    !< Number of sediment flux/diagenesis zones
  INTEGER :: NSMZM   !< Maximum NSMZ

  REAL :: SMCSHSCH
  REAL :: SMFD1H2S
  REAL :: SMFD1NH4
  REAL :: SMO2NH4
  REAL :: SMO2NO3
  REAL :: SMPOCR

  REAL :: SMWQASC   !< Silica-to-carbon ratio for algae diatoms
  REAL :: SMDIFT    !< Diffusion coefficient for sediment temperature (m2/sec)
  REAL :: SMM1      !< Solid concentrations in Layer 1 (Kg/L)
  REAL :: SMM2      !< Solid concentrations in Layer 2 (Kg/L)   
  REAL :: SMKMDP    !< Particle mixing half-saturation constant for oxygen (mg/L)
  REAL :: SMKBST    !< First-order decay rate for accumulated benthic stress (1/day)
  REAL :: XSMDPMIN  !< Minimum diffusion coefficient for particle mixing (m^2/d)
  REAL :: SMRBIBT   !< Ratio of bio-irrigation to bioturbation (unitless)
  REAL :: SMO2BS    !< Critical overlying oxygen concentration below which benthic hysteresis occurs (mg/L)
  REAL :: SMHYLAG   !< Time duration for which the maximum or minimum stress is retained (days)
  REAL :: SMHYDUR   !< Critical hypoxia duration; if less than this value, no hysteresis occurs (days)
  REAL :: SMKMNH4   !< Nitrification half-sat. constant for ammonium (gN/m^3)
  REAL :: SMKMO2N   !< Nitrification half-sat. constant for dissolved oxygen (gO2/m^3)
  REAL :: SMP2PO4   !< Partition coefficient, ratio of particulate to dissolved PO4 in layer 2 (L/Kg)
  REAL :: SMCO2PO4  !< Critical dissolved oxygen for PO4 sorption (mg/L)
  REAL :: SMO2C     !< Stoichiometric coefficient for carbon diagenesis consumed by H2S oxidation (gO2/gC)
  REAL :: SMKMPSI   !< Silica dissolution half-saturation constant for PSi (g Si/m^3
  REAL :: SMSISAT   !< Saturation concentration of silica in pore water (g Si/m^3)
  REAL :: SMP2SI    !< Partition coefficient for Si in Layer 2, controls sorption of dissolved silica to solids (L/Kg)
  REAL :: SMDP1SI   !< Factor that enhances sorption of silica in layer 1 when D.O. exceeds DOcSi (unitless)
  REAL :: SMCO2SI   !< Critical dissolved oxygen for silica sorption in layer 1 (mg/L)
  REAL :: SMJDSI    !< Detrital flux of particulate biogenic silica from sources other than diatom algae (gSi/m^2/d)
  REAL :: SMFD2H2S  !< Dissolved fraction of H2S in Layer 2
  REAL :: SMFD2NH4  !< Dissolved fraction of NH4 in Layer 2
  REAL :: SMFD2PO4  !< Dissolved fraction of PO4 in Layer 2
  REAL :: SMFD2SI   !< Dissolved fraction of SI  in Layer 2
  REAL :: SMFP1H2S  !< Particulate fraction of H2S in Layer 1 
  REAL :: SMFP1NH4  !< Particulate fraction of NH4 in Layer 1 
  REAL :: SMFP2H2S  !< Particulate fraction of H2S in Layer 2 
  REAL :: SMFP2NH4  !< Particulate fraction of NH4 in Layer 2 
  REAL :: SMFP2PO4  !< Particulate fraction of PO4 in Layer 2 
  REAL :: SMFP2SI   !< Particulate fraction of SI  in Layer 2 
  REAL :: SM1OKMDP
  REAL :: WQTDsMIN
  REAL :: WQTDsMAX
  REAL :: WQTDsINC
  
  ! *** Integer Array Variables
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: ISMT       !<  Current Lookup Index Based On The Current Sediment Temperature
  INTEGER,ALLOCATABLE,DIMENSION(:)  :: ISMZMAP
    
  ! *** Real Array Variables
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDWQANC   !< Nitrogen-to-carbon ratio for algae 
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQDFB      !< Coupling between water column and sediment diagenesis model
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMD1PO4
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMD1SI
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDD       !< Diffusion coefficient in pore water (m2/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDFSI
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDGFC
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDGFN
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDGFP
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDP       !< Apparent diffusion coefficient for particle mixing (m2/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDP1PO4   !< Factor to enhance sorption of PO4 in layer 1 when DO is less than DOcPO4 (unitless)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDPMIN    !< Minimum diffusion coefficient for particle mixing divided by H2
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMDTOH
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMFCBA     !< fraction of POC from algae routed to G classes
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMFCR      !< fraction of water column refractory POC routed to G-classes
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMFNBA     !< fraction of PON from algae routed to G classes
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMFNR      !< fraction of water column refractory PON routed to G-classes
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMFPBA     !< fraction of POP from algae routed to G classes
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMFPR      !< fraction of water column refractory POP routed to G-classes
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMHODT     !< Benthic sediment depth divided by time step
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMHYPD
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMJAQH2S
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMJDEN
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMJGCH4
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMJNIT 
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMK1H2S
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMK1NO3    !< Reaction velocity for denitrification in layer 1 at 20 degC (m/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMK2NO3    !< Reaction velocity for denitrification in layer 2 at 20 degC (m/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMKL12     !< Dissolved phase mixing coefficient
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMKNH4     !< Optimal reaction velocity for nitrification at 20 degC (m/day)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMSS       !< Surface mass transfer coefficient
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTD1CH4
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTD2CH4
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMTDCD
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTDDD     !< Constant for temperature adjustment for Dd raised to T-20
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTDDP     !< Constant for temperature adjustment for Dp raised to T-20
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMTDND
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTDNH4    !< Constant for temperature adjustment for NH4 raised to T-20
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTDNO3    !< Constant for temperature adjustment for NO3 raised to T-20
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMTDPD
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTDSI
  
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMTMP
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMW12      !< Particulate phase mixing coefficient
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMW2       !< Sediment burial rate (cm/year)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMW2DTOH
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMW2PHODT  
  REAL,ALLOCATABLE,DIMENSION(:)     :: SODMULT    !< Factor to enhance magnitude of sediment oxygen demand (unitless)
  REAL,ALLOCATABLE,DIMENSION(:)     :: SM1DIFT
  REAL,ALLOCATABLE,DIMENSION(:)     :: SM2DIFT
  REAL,ALLOCATABLE,DIMENSION(:)     :: SM2NO3     !< WATER QUALITY MODEL VARIABLE
  REAL,ALLOCATABLE,DIMENSION(:)     :: SM2PO4     !< WATER QUALITY MODEL VARIABLE
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMHSED     !< WATER QUALITY MODEL VARIABLE

  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMPOC      !< Particulate organic Carbon     in the sediments
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMPOP      !< Particulate organic Phosphorus in the sediments 
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SMPON      !< Particulate organic Nitrogen   in the sediments
  REAL,ALLOCATABLE,DIMENSION(:)     :: SMPSI      !< Sediment Silica

  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM1H2S !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM1NH4 !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM1NO3 !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM1PO4 !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM1SI  !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM2H2S !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM2NH4 !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SM2SI  !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SMBST  !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SMCSOD !< WATER QUALITY MODEL VARIABLE
    
  REAL,TARGET,ALLOCATABLE,DIMENSION(:,:) :: SMDFC  !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:,:) :: SMDFN  !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:,:) :: SMDFP  !< WATER QUALITY MODEL VARIABLE
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SMNSOD !< WATER QUALITY MODEL VARIABLE
  LOGICAL,ALLOCATABLE,DIMENSION(:)       :: SMHYST !< WATER QUALITY MODEL VARIABLE
  
  REAL,TARGET,ALLOCATABLE,DIMENSION(:)   :: SMT    !< Sediment temperature (degC)
    
  ! *** Zooplankton variables
  INTEGER :: NZO                                 !< Index of zooplankton groups
  INTEGER,ALLOCATABLE,DIMENSION(:)    :: IWQZT   !< Index for look-up table for temperature dependency
  REAL :: WQTDZMIN    !< Lower end of temperature for look-up table
  REAL :: WQTDZMAX    !< Upper end of temperature for look-up table
  REAL :: WTEMPZ
  REAL :: WQTDZINC    !< Increment of temperature for look-up table
  
  TYPE ZOOPLGROUP
    CHARACTER(LEN = 20) :: ZOONAME
    INTEGER             :: IDZ
    INTEGER             :: ISPREY         !< Flag to indicate prey or predator
    INTEGER             :: ISPREDATOR     !< Flag to indicate prey or predator
    
    REAL :: CTZ         !< Carbon threshold for zooplankton grazing (gC/m3)
    REAL :: ANCZ        !< Nitrogen to Carbon ratio (gN/gC)
    REAL :: APCZ        !< Phosphorus to Carbon ratio (gP/gC)
    REAL :: ASCZ        !< Phosphorus to Carbon ratio (gSi/gC)
    REAL :: KHCZ        !< Prey density at which zooplankton grazing is halved (gC/m3)
    REAL :: ULZ         !< Utilization of labile particulate organic carbon
    REAL :: URZ         !< Utilization of refractory particulate organic carbon
    REAL :: UDZ         !< Utilization of Dissolved organic carbon
    REAL :: UZPL        !< Utilization of zooplankton as prey
    REAL :: KTGZ1       !< Effect of temperature blow optimal on grazing
    REAL :: KTGZ2       !< Effect of temperature above optimal on grazing
    REAL :: TMZG1       !< Lower end of optimal temperature for grazing
    REAL :: TMZG2       !< Upper end of optimal temperature for grazing
    REAL :: KTBZ        !< Effect of temperature on metabolism of zooplankton
    REAL :: TRZB        !< Reference temperature for zooplankton metabolism
    REAL :: KTPZ        !< Effect of temperature on predation of zooplankton
    REAL :: TRZP        !< Reference temperature for zooplankton predation
    REAL :: DOCRIT      !< Critical DO concentration below which zooplankton death occurs (g DO/m3)
    REAL :: DZEROZ      !< Zooplankton death at zero DO concentration
    REAL :: FCDDZ       !< Fraction of death carbon produced by zooplankton as DOC
    REAL :: FCLDZ       !< Fraction of death carbon produced by zooplankton as LPOC 
    REAL :: FCRDZ       !< Fraction of death carbon produced by zooplankton as LROC
    REAL :: FCDPZ       !< Fraction of predated carbon produced by zooplankton as DOC
    REAL :: FCLPZ       !< Fraction of predated carbon produced by zooplankton as LPOC
    REAL :: FCRPZ       !< Fraction of predated carbon produced by zooplankton as RPOC
    REAL :: FNLDZ       !< Fraction of death nitrogen produced by zooplankton as LPON
    REAL :: FNRDZ       !< Fraction of death nitrogen produced by zooplankton as RPON
    REAL :: FNDDZ       !< Fraction of death nitrogen produced by zooplankton as DON
    REAL :: FNIDZ       !< Fraction of death nitrogen produced by zooplankton as NH4
    REAL :: FNLPZ       !< Fraction of predated nitrogen produced by zooplankton as LPON
    REAL :: FNRPZ       !< Fraction of predated nitrogen produced by zooplankton as RPON
    REAL :: FNDPZ       !< Fraction of predated nitrogen produced by zooplankton as DON
    REAL :: FNIPZ       !< Fraction of predated nitrogen produced by zooplankton as NH4
    REAL :: FNDBZ       !< Fraction of basal metabolism nitrogen produced by zooplankton as DON
    REAL :: FNIBZ       !< Fraction of basal metabolism nitrogen produced by zooplankton as NH4
    REAL :: FPLDZ       !< Fraction of death phosphorus produced by zooplankton as LPOP
    REAL :: FPRDZ       !< Fraction of death phosphorus produced by zooplankton as RPOP
    REAL :: FPDDZ       !< Fraction of death phosphorus produced by zooplankton as DOP
    REAL :: FPIDZ       !< Fraction of death phosphorus produced by zooplankton as PO4
    REAL :: FPLPZ       !< Fraction of predated phosphorus produced by zooplankton as LPOP
    REAL :: FPRPZ       !< Fraction of predated phosphorus produced by zooplankton as RPOP
    REAL :: FPDPZ       !< Fraction of predated phosphorus produced by zooplankton as DOP
    REAL :: FPIPZ       !< Fraction of predated phosphorus produced by zooplankton as PO4
    REAL :: FPDBZ       !< Fraction of basal metabolism phosphorus produced by zooplankton as DOP
    REAL :: FPIBZ       !< Fraction of basal metabolism phosphorus produced by zooplankton as PO4
    REAL :: FSPDZ       !< Fraction of death Silica produced by zooplankton as SU
    REAL :: FSPPZ       !< Fraction of predated Silica produced by zooplankton as SU
    REAL :: FSADZ       !< Fraction of death Silica produced by zooplankton as SA
    REAL :: FSAPZ       !< Fraction of predated Silica produced by zooplankton as SA
    
    REAL,ALLOCATABLE,DIMENSION(:)   :: RMAXZ  !< Maximum ration of zooplankton (g prey C/g zooplankton C/d)
    REAL,ALLOCATABLE,DIMENSION(:)   :: BMRZ   !< Metabolic rate of zooplankton at reference temperature (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)   :: PRRZ   !< Predation rate of zooplankton at reference temperature (1/day)   
    REAL,ALLOCATABLE,DIMENSION(:)   :: UBZ    !< Utilization of algaes by zooplankton
    REAL,ALLOCATABLE,DIMENSION(:)   :: WQGZ   !< Zooplankton growth rate (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)   :: WQBZ   !< Zooplankton metabolic rate (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)   :: WQPZ   !< Zooplankton predation rate (1/day)
    REAL,ALLOCATABLE,DIMENSION(:)   :: WQDZ   !< Zooplankton death rate (1/day)
  END TYPE 
  
  TYPE(ZOOPLGROUP),ALLOCATABLE,DIMENSION(:) :: ZOOPL
  
  ! *** Real Array Variables
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: PRAZ     !< Prey available to zooplankton (gC/m3) 
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: BAZ      !< Available portion of algal group to zooplankton (gC/m3)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: DOCAZ    !< Available portion of dissolved organic carbon to zooplankton (gC/m3)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: RPOCAZ   !< Available portion of refractory particulate organic carbon to zooplankton (gC/m3)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: LPOCAZ   !< Available portion of labile particulate organic carbon to zooplankton (gC/m3)
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDGZ   !< Temperature dependency of zooplankton grazing
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDBZ   !< Exponential function of temperature for zooplankton metabolism
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: WQTDPZ   !< Exponential function of temperature for zooplankton predation
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQZKK    !< Related to matrix K1 in the kinetic equations
  REAL,ALLOCATABLE,DIMENSION(:)     :: WQZRR    !< Accumulation of sources and sinks in the kinetic equations
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: WQWPSZ   !< Point source loading
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: FRLP     
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: FRRP     
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: FRSI     
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: SBZPAL   !< Total effect of zooplankton on phytoplankton
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SLPOCZ   !< Total effect of zooplanktons on LPOC
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SRPOCZ   !< Total effect of zooplanktons on RPOC
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SDOCZ    !< Total effect of zooplanktons on DOC
    
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SLPONZ   !< Total effect of zooplanktons on LPON
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SRPONZ   !< Total effect of zooplanktons on RPON
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SDONZ    !< Total effect of zooplanktons on DON
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SNH4Z    !< Total effect of zooplanktons on NH4
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SLPOPZ   !< Total effect of zooplanktons on LPOP
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SRPOPZ   !< Total effect of zooplanktons on RPOP
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SDOPZ    !< Total effect of zooplanktons on DOP
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SPO4Z    !< Total effect of zooplanktons on PO4
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SSUZ     !< Total effect of zooplanktons on SU
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SSAZ     !< Total effect of zooplanktons on SA
  REAL,ALLOCATABLE,DIMENSION(:,:)   :: SDOZ     !< Total effect of zooplanktons on DO
  
  Contains
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ_Allocate
  !
  !> @details  Allocating and Initializing of WQ array variables
  !---------------------------------------------------------------------------!
  Subroutine WQ_Allocate
    Use Allocate_Initialize      
    Integer :: NAL
    
    ! *** Call AllocateDSI integer arrays
    Call AllocateDSI(IWQCBE,  NBBEM, 0)
    Call AllocateDSI(IWQCBN,  NBBNM, 0)
    Call AllocateDSI(IWQCBS,  NBBSM, 0)
    Call AllocateDSI(IWQCBW,  NBBWM, 0)
    Call AllocateDSI(IWQOBE,  NBBEM, NWQVM, 0)
    Call AllocateDSI(IWQOBN,  NBBNM, NWQVM, 0)
    Call AllocateDSI(IWQOBS,  NBBSM, NWQVM, 0)
    Call AllocateDSI(IWQOBW,  NBBWM, NWQVM, 0)
    Call AllocateDSI(IWQPSC,  LCM, KCM, 0)
    Call AllocateDSI(IWQPSV,  LCM, KCM, 0)
    Call AllocateDSI(IWQZMAP, LCM, KCM, 0)

    Call AllocateDSI(JWQCBE,  NBBEM, 0)
    Call AllocateDSI(JWQCBN,  NBBNM, 0)
    Call AllocateDSI(JWQCBS,  NBBSM, 0)
    Call AllocateDSI(JWQCBW,  NBBWM, 0)

    Call AllocateDSI(ICPSL,   NWQPSM, 0)
    Call AllocateDSI(JCPSL,   NWQPSM, 0)
    Call AllocateDSI(KCPSL,   NWQPSM, 0)
    Call AllocateDSI(MVPSL,   NWQPSM, 0)
    Call AllocateDSI(MWQCTLT, NWQCSRM, NWQVM, 0)

    !Call AllocateDSI(NWQCSR, NTSWQVM, 0)  ! *** NWQCSR is used in SCANWQ before this subroutine is called

    ! *** Call AllocateDSI real arrays
    Call AllocateDSI(TWQ,       LCM, 0.0)
    Call AllocateDSI(SWQ,       LCM, 0.0)
    Call AllocateDSI(VOLWQ,     LCM, 0.0)

    Call AllocateDSI(WQAPC,     LCM, 0.0)
    Call AllocateDSI(WQATM,     MAXWQ, NWQZM, 0.0)  
    Call AllocateDSI(WQBFCOD,   LCM, 0.0)
    Call AllocateDSI(WQBFNH4,   LCM, 0.0)
    Call AllocateDSI(WQBFNO3,   LCM, 0.0)
    Call AllocateDSI(WQBFO2,    LCM, 0.0)
    Call AllocateDSI(WQBFPO4D,  LCM, 0.0)
    Call AllocateDSI(WQBFSAD,   LCM, 0.0)
    Call AllocateDSI(WQCHL,     LCM, KCM, 0.0)
    Call AllocateDSI(WQDFLC,    LCM, 0.0)
    Call AllocateDSI(WQDFLN,    LCM, 0.0)
    Call AllocateDSI(WQDFLP,    LCM, 0.0)
    Call AllocateDSI(WQDFRC,    LCM, 0.0)
    Call AllocateDSI(WQDFRN,    LCM, 0.0)
    Call AllocateDSI(WQDFRP,    LCM, 0.0)
    Call AllocateDSI(WQDFSI,    LCM, 0.0)
    Call AllocateDSI(WQOBCE,    NBBEM, 2, NWQVM, 0.0)
    Call AllocateDSI(WQOBCN,    NBBNM, 2, NWQVM, 0.0)
    Call AllocateDSI(WQOBCS,    NBBSM, 2, NWQVM, 0.0)
    Call AllocateDSI(WQOBCW,    NBBWM, 2, NWQVM, 0.0)
    Call AllocateDSI(WQPO4D,    LCM, KCM, 0.0)
    Call AllocateDSI(WQSAD,     LCM, KCM, 0.0)
    Call AllocateDSI(WQWPSL,    LCM,  KCM, NWQVM, 0.0)
    Call AllocateDSI(XSMO20,    LCM, 0.0)

    ! *** Call AllocateDSI algae class varibales
    Allocate(ALGAES(NALGAEM))
    
    Do NAL = 1, NALGAEM
      ALGAES(NAL).IDN    = NAL
      Call AllocateDSI(ALGAES(NAL).WQBMRA,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQDOP,    NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQDRA,    NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKDCALM, NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKHRA,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQPMA,    NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQPRRA,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQWS,     NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQBMIN,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKBP,    NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMV,    NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMVMIN, NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMVA,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMVB,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMVC,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMVD,   NWQZM, 0.0)
      Call AllocateDSI(ALGAES(NAL).WQKMVE,   NWQZM, 0.0)
    Enddo 
    
    Call AllocateDSI(LAYERBOT,       NALGAEM,  LCM,      1)
    Call AllocateDSI(LAYERTOP,       NALGAEM,  LCM,      KC)
    Call AllocateDSI(HEIGHT_MAC,     LCM,     -NALGAEM,  0.)
    Call AllocateDSI(DIAMETER_MAC,   LCM,      NALGAEM,  0.)
    Call AllocateDSI(LayerRatio_MAC, LCM,      KCM,     NALGAEM,  0.)
    
    Call AllocateDSI(IWQKA,    NWQZM,    0)
    Call AllocateDSI(IBENMAP,  LCM,   2, 0)
    Call AllocateDSI(MWQPSR,   NWQPSRM,  0)
    Call AllocateDSI(MWQPTLT,  NWQPSRM,  0)
    Call AllocateDSI(IMWQZT,   LCM,   0)
    Call AllocateDSI(IMWQZT1,  LCM,   0)
    Call AllocateDSI(IMWQZT2,  LCM,   0)
    Call AllocateDSI(IWQT,     LCM,   0)
    
    Call AllocateDSI(WQPA,     LCM,   NALGAEM, 0.0)
    Call AllocateDSI(WQBM,     LCM,   NALGAEM, 0.0)
    Call AllocateDSI(WQPR,     LCM,   NALGAEM, 0.0)
    Call AllocateDSI(WQBSET,   LCM,   2,  NALGAEM, 0.0)
    Call AllocateDSI(WQTDG,    NWQTDM, NALGAEM, 0.0)
    Call AllocateDSI(WQTDR,    NWQTDM, NALGAEM, 0.0)
    Call AllocateDSI(WQTDP,    NWQTDM, NALGAEM, 0.0)
    Call AllocateDSI(WQTDGP,   NWQTDM, 0.0)
    Call AllocateDSI(WQTD1FCB, NWQTDM, 0.0)
    Call AllocateDSI(WQTD2FCB, NWQTDM, 0.0)
    Call AllocateDSI(UHEQ,     LCM, 0.0)
    Call AllocateDSI(WQHT,     KCM, 0.0)
    Call AllocateDSI(WQWDSL,   LCM, KCM,  NWQVM,  0.0)  
    Call AllocateDSI(DDOAVG,   LCM, KCM,  0.0)
    Call AllocateDSI(DDOMAX,   LCM, KCM,  0.0)
    Call AllocateDSI(DDOMIN,   LCM, KCM,  0.0)
    Call AllocateDSI(RLIGHTC,  LCM, KCM,  0.0)
    Call AllocateDSI(RLIGHTT,  LCM, KCM,  0.0)
    Call AllocateDSI(REAC,     NWQZM,   0.0)
    Call AllocateDSI(WQKDC,    NWQZM,   0.0)
    Call AllocateDSI(WQKRO,    NWQZM,   0.0)
    Call AllocateDSI(WQKTR,    NWQZM,   0.0)
    Call AllocateDSI(WQKHCOD,  NWQZM,   0.0)
    Call AllocateDSI(WQKCD,    NWQZM,   0.0)
    Call AllocateDSI(WQTDHDR,  NWQTDM,  0.0)
    Call AllocateDSI(WQTDMNL,  NWQTDM,  0.0)
    Call AllocateDSI(WQKSUA,   NWQTDM,  0.0)
    Call AllocateDSI(WQKCOD,   NWQTDM,  NWQZM, 0.0)
    Call AllocateDSI(WQTDNIT,  NWQTDM,  0.0)
    Call AllocateDSI(WQTDKR,   NWQTDM,  NWQZM, 0.0)
    Call AllocateDSI(WQTDTAM,  NWQTDM,  0.0)
    Call AllocateDSI(WQTAMP,   LCM,   KCM,   0.0)
    Call AllocateDSI(WQWSLP,   NWQZM,   0.0)
    Call AllocateDSI(WQWSRP,   NWQZM,   0.0)
    Call AllocateDSI(WQWSS,    NWQZM,   0.0)
    Call AllocateDSI(WQATML,   LCM,   KCM,  NWQVM, 0.0)
    Call AllocateDSI(XBENMUD,  LCM,   0.0)
    Call AllocateDSI(TSSRD,    NDASER,  0.0)
    Call AllocateDSI(SOLFRD,   NDASER,  0.0)
    Call AllocateDSI(SOLSRD,   NDASER,  0.0)
    Call AllocateDSI(TCWQPSR,  NWQPSRM, 0.0)
    Call AllocateDSI(TWQPSER,  NDWQPSR, NWQPSRM, 0.0)
    Call AllocateDSI(WQPSSER,  NDWQPSR, NWQVM,  NWQPSRM, 0.0)
    Call AllocateDSI(WQI0BOT,  LCM,     0.0)
    Call AllocateDSI(DZWQ,     LCM,   0.0)
    Call AllocateDSI(WQKRPC,   LCM,   0.0)
    Call AllocateDSI(WQKRPN,   LCM,   0.0)
    Call AllocateDSI(WQKRPP,   LCM,   0.0)
    Call AllocateDSI(WQKLPC,   LCM,   0.0)
    Call AllocateDSI(WQKLPN,   LCM,   0.0)
    Call AllocateDSI(WQKLPP,   LCM,   0.0)
    Call AllocateDSI(WQKHR,    LCM,   0.0)
    Call AllocateDSI(WQDENIT,  LCM,   0.0)
    Call AllocateDSI(WQPN,     LCM,   NALGAEM, 0.0)
    Call AllocateDSI(WQNIT,    LCM,   0.0)    
    Call AllocateDSI(O2WQ,     LCM,   0.0)    
    Call AllocateDSI(RNH4WQ,   LCM,   0.0)    
    Call AllocateDSI(RNO3WQ,   LCM,   0.0) 
    Call AllocateDSI(WQDOS,    LCM,   0.0)    
    Call AllocateDSI(WQH10,    LCM,   0.0)    
    Call AllocateDSI(XDOSAT,   LCM,   KCM, 0.0)    
    Call AllocateDSI(WQP19,    LCM,   0.0)    
    Call AllocateDSI(WQKRDOS,  LCM,   0.0)    
    Call AllocateDSI(WQKK,     LCM,   0.0)    
    Call AllocateDSI(WQRR,     LCM,   0.0)      
    Call AllocateDSI(WQRPSET,  LCM,   2,  0.0)    
    Call AllocateDSI(WQLPSET,  LCM,   2,  0.0)    
    Call AllocateDSI(WQWSSET,  LCM,   2,  0.0)    
    Call AllocateDSI(PO4DWQ,   LCM,   0.0)    
    Call AllocateDSI(RNH4NO3,  LCM,   0.0)    
    Call AllocateDSI(WQOBTOT,  LCM,   0.0)    
    Call AllocateDSI(WQKDON,   LCM,   0.0)    
    Call AllocateDSI(WQKDOP,   LCM,   0.0)   
    Call AllocateDSI(WQT10,    LCM,   0.0)    
    Call AllocateDSI(WQT17,    LCM,   0.0)    
    Call AllocateDSI(WQN17,    LCM,   0.0)    
    Call AllocateDSI(WQO18,    LCM,   0.0)    
    Call AllocateDSI(WQR20,    LCM,   0.0)    
    Call AllocateDSI(SMAC,     LCM,   1.0)        !< Fixed biota flag set to "on"
    
    ! *** Special Cases
    Allocate(ICWQTS(0:NWQVM,NWQTSM))
    Allocate(WQV(LCM,KCM,0:NWQVM))
    Allocate(WQVO(LCM,0:KCM,0:NWQVM))
    Allocate(WQWPSLC(0:NWQPSM,NWQVM))
    Allocate(WQPSSRT(NWQVM,0:NWQPSRM))
    ICWQTS  = 0
    WQV     = 0.0
    WQVO    = 0.0
    WQWPSLC = 0.0
    WQPSSRT = 0.0
    
    Deallocate(WQKEB)
    Call AllocateDSI(WQKEB, NWQZM, 0.0)

    ! *** Dissolved carbon dioxide variables
    Call AllocateDSI(CO2WQ,     LCM, 0.0)
    Call AllocateDSI(CDOSATIDX, LCM, 0.0)
    Call AllocateDSI(WQCDOS,    LCM, 0.0)
    Call AllocateDSI(WQITOP,    LCM, KCM, 0.0)
    Call AllocateDSI(WQP22,     LCM, 0.0)
    Call AllocateDSI(WQKRCDOS,  LCM, 0.0)

  End Subroutine WQ_Allocate

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine SD_Allocate
  !
  !> @details  Allocating and Initializing of SD array variables
  !---------------------------------------------------------------------------!

  Subroutine SD_Allocate
    Use Allocate_Initialize 
    
    ! *** Call AllocateDSI integer arrays
    Call AllocateDSI(ISMT,     LCM, 0)
    Call AllocateDSI(ISMZMAP,  LCM, 0)
    
    ! *** Call AllocateDSI real arrays
    Call AllocateDSI(SMDWQANC,  NALGAEM, 0.0)
    Call AllocateDSI(WQDFB,     LCM,  NALGAEM, 0.0)
    Call AllocateDSI(SMD1PO4,   LCM,  0.0)
    Call AllocateDSI(SMD1SI,    LCM,  0.0)
    Call AllocateDSI(SMDD,      NSMZM,  0.0)
    Call AllocateDSI(SMDFSI,    LCM,  0.0)
    Call AllocateDSI(SMDGFC,    LCM,  0.0)
    Call AllocateDSI(SMDGFN,    LCM,  0.0)
    Call AllocateDSI(SMDGFP,    LCM,  0.0)
    Call AllocateDSI(SMDP,      NSMZM,  0.0)
    Call AllocateDSI(SMDP1PO4,  NSMZM,  0.0)
    Call AllocateDSI(SMDPMIN,   NSMZM,  0.0)
    Call AllocateDSI(SMDTOH,    NSMZM,  0.0)
    Call AllocateDSI(SMFCBA,    NALGAEM, NSMGM, 0.0)
    Call AllocateDSI(SMFNBA,    NALGAEM, NSMGM, 0.0)
    Call AllocateDSI(SMFPBA,    NALGAEM, NSMGM, 0.0)
    Call AllocateDSI(SMFCR,     NSMZM,  NSMGM, 0.0)
    Call AllocateDSI(SMFNR,     NSMZM,  NSMGM, 0.0)
    Call AllocateDSI(SMFPR,     NSMZM,  NSMGM, 0.0)
    Call AllocateDSI(SMHODT,    NSMZM,  0.0)
    Call AllocateDSI(SMHYPD,    LCM,  0.0)
    Call AllocateDSI(SMJAQH2S,  LCM,  0.0)
    Call AllocateDSI(SMJDEN,    LCM,  0.0)
    Call AllocateDSI(SMJGCH4,   LCM,  0.0)
    Call AllocateDSI(SMJNIT,    LCM,  0.0)
    Call AllocateDSI(SMK1H2S,   NWQTDM, 0.0)
    Call AllocateDSI(SMK1NO3,   NSMZM,  0.0)
    Call AllocateDSI(SMK2NO3,   NSMZM,  0.0)
    Call AllocateDSI(SMKL12,    LCM,    0.0)
    Call AllocateDSI(SMKNH4,    NSMZM,  0.0)
    Call AllocateDSI(SMSS,      LCM,  0.0)
    Call AllocateDSI(SMTD1CH4,  NWQTDM, 0.0)
    Call AllocateDSI(SMTD2CH4,  NWQTDM, 0.0)
    Call AllocateDSI(SMTDCD,    NWQTDM, NSMGM, 0.0)
    Call AllocateDSI(SMTDDD,    NWQTDM, 0.0)
    Call AllocateDSI(SMTDDP,    NWQTDM, 0.0)
    Call AllocateDSI(SMTDND,    NWQTDM, NSMGM, 0.0)
    Call AllocateDSI(SMTDNH4,   NWQTDM, 0.0)
    Call AllocateDSI(SMTDNO3,   NWQTDM, 0.0)
    Call AllocateDSI(SMTDPD,    NWQTDM, NSMGM, 0.0)  
    Call AllocateDSI(SMTDSI,    NWQTDM, 0.0) 
    Call AllocateDSI(SMTMP,     LCM,    0.0)   
    Call AllocateDSI(SMW12,     LCM,    0.0)
    Call AllocateDSI(SMW2,      NSMZM,  0.0)
    Call AllocateDSI(SMW2DTOH,  NSMZM,  0.0)
    Call AllocateDSI(SMW2PHODT, NSMZM,  0.0)
    Call AllocateDSI(SODMULT,   NSMZM,  0.0)
    Call AllocateDSI(SM1DIFT,   NSMZM,  0.0)
    Call AllocateDSI(SM2DIFT,   NSMZM,  0.0)
    Call AllocateDSI(SM2NO3,    LCM,  0.0)
    Call AllocateDSI(SM2PO4,    LCM,  0.0)
    Call AllocateDSI(SMHSED,    NSMZM,  0.0)
    Call AllocateDSI(SMPOP,     LCM,  NSMGM, 0.0) 
    Call AllocateDSI(SM1H2S,    LCM,  0.0)
    Call AllocateDSI(SM1NH4,    LCM,  0.0)
    Call AllocateDSI(SM1NO3,    LCM,  0.0)
    Call AllocateDSI(SM1PO4,    LCM,  0.0)
    Call AllocateDSI(SM1SI,     LCM,  0.0)
    Call AllocateDSI(SM2H2S,    LCM,  0.0)
    Call AllocateDSI(SM2NH4,    LCM,  0.0)
    Call AllocateDSI(SM2SI,     LCM,  0.0)
    Call AllocateDSI(SMBST,     LCM,  0.0)
    Call AllocateDSI(SMCSOD,    LCM,  0.0)
    Call AllocateDSI(SMDFC,     LCM,  NSMGM, 0.0)
    Call AllocateDSI(SMDFN,     LCM,  NSMGM, 0.0)
    Call AllocateDSI(SMDFP,     LCM,  NSMGM, 0.0)
    Call AllocateDSI(SMPOC,     LCM,  NSMGM, 0.0)
    Call AllocateDSI(SMPON,     LCM,  NSMGM, 0.0)
    Call AllocateDSI(SMNSOD,    LCM,  0.0)
    Call AllocateDSI(SMPSI,     LCM,  0.0)
    Call AllocateDSI(SMT,       LCM,  0.0)
       
    ! *** Special cases    
    ALLOCATE(SMHYST(LCM))
     
  End Subroutine SD_Allocate
    
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ3DINP
  !
  !> @details  Allocating and Initializing of SD array variables
  !---------------------------------------------------------------------------!

  Subroutine Zoo_Allocate
    Use Allocate_Initialize      
    Integer :: NZO
        
    ! *** Call AllocateDSI integer arrays
    Call AllocateDSI(IWQZT,     LCM, 0)
    
    ! *** Call AllocateDSI algae class varibales
    Allocate(ZOOPL(NZOOPL))
    
    Do NZO = 1, NZOOPL
      Call AllocateDSI(ZOOPL(NZO).RMAXZ, NWQZM,  0.0)
      Call AllocateDSI(ZOOPL(NZO).BMRZ,  NWQZM,  0.0)
      Call AllocateDSI(ZOOPL(NZO).PRRZ,  NWQZM,  0.0)
      Call AllocateDSI(ZOOPL(NZO).UBZ,   NALGAEM, 0.0)
      Call AllocateDSI(ZOOPL(NZO).WQGZ,  LCM,  0.0)
      Call AllocateDSI(ZOOPL(NZO).WQBZ,  LCM,  0.0)
      Call AllocateDSI(ZOOPL(NZO).WQPZ,  LCM,  0.0)
      Call AllocateDSI(ZOOPL(NZO).WQDZ,  LCM,  0.0)
    Enddo 
     
    Call AllocateDSI(PRAZ,    LCM,  NZOOPL, 0.0)
    Call AllocateDSI(BAZ,     LCM,  NZOOPL, NALGAEM, 0.0)
    Call AllocateDSI(DOCAZ,   LCM,  NZOOPL, 0.0)
    Call AllocateDSI(RPOCAZ,  LCM,  NZOOPL, 0.0)
    Call AllocateDSI(LPOCAZ,  LCM,  NZOOPL, 0.0)
    Call AllocateDSI(WQTDGZ,  NWQTDM, NZOOPL, 0.0)
    Call AllocateDSI(WQTDBZ,  NWQTDM, NZOOPL, 0.0)
    Call AllocateDSI(WQTDPZ,  NWQTDM, NZOOPL, 0.0)
    Call AllocateDSI(WQZKK,   LCM,  0.0)
    Call AllocateDSI(WQZRR,   LCM,  0.0)
    Call AllocateDSI(WQWPSZ,  LCM,  KCM,    NZOOPL, 0.0)
    Call AllocateDSI(FRRP,    LCM,  NZOOPL, 0.0)
    Call AllocateDSI(FRLP,    LCM,  NZOOPL, 0.0)
    Call AllocateDSI(FRSI,    LCM,  NZOOPL, 0.0)
    Call AllocateDSI(SBZPAL,  LCM,  KCM, NALGAEM, 0.0)
    Call AllocateDSI(SLPOCZ,  LCM,  KCM, 0.0)
    Call AllocateDSI(SRPOCZ,  LCM,  KCM, 0.0)
    Call AllocateDSI(SDOCZ,   LCM,  KCM, 0.0)
    Call AllocateDSI(SLPONZ,  LCM,  KCM, 0.0)
    Call AllocateDSI(SRPONZ,  LCM,  KCM, 0.0)
    Call AllocateDSI(SDONZ,   LCM,  KCM, 0.0)
    Call AllocateDSI(SNH4Z,   LCM,  KCM, 0.0)
    Call AllocateDSI(SLPOPZ,  LCM,  KCM, 0.0)
    Call AllocateDSI(SRPOPZ,  LCM,  KCM, 0.0)
    Call AllocateDSI(SDOPZ,   LCM,  KCM, 0.0)
    Call AllocateDSI(SPO4Z,   LCM,  KCM, 0.0)
    Call AllocateDSI(SSUZ,    LCM,  KCM, 0.0)
    Call AllocateDSI(SSAZ,    LCM,  KCM, 0.0)
    Call AllocateDSI(SDOZ,    LCM,  KCM, 0.0)
    
  End Subroutine Zoo_Allocate
  
  
End Module Variables_WQ
