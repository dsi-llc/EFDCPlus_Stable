! ----------------------------------------------------------------------
!   This file is a part of EFDC +
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------

REAL FUNCTION CALVOLTERM(HP, HPKCI, STCUV, UHE, VHE, HU, HV, TEM_KC, SAL_KC, TATMT, PATMT, WINDST,    &
                         VEL_MAX, DEP_MIN, MOLW_I, C_HE, AIRCON, TCOEFF, MULT, WCCON, KL_OPT, PFTW)

  use DIFFUSER_MODULE, only:FUNDEN
  implicit none  

  real, external :: ViscosityW

  ! *** Input Variables
  integer, intent(IN) :: KL_OPT
  real, intent(IN)    :: HP, HPKCI, HU, HV, UHE, VHE, STCUV
  real, intent(IN)    :: TEM_KC, SAL_KC
  real, intent(IN)    :: VEL_MAX, DEP_MIN, MOLW_I, AIRCON, TCOEFF, C_HE, WCCON
  real, intent(IN)    :: TATMT, PATMT, WINDST
  real, intent(IN)    :: MULT
  real, intent(IN)    :: PFTW

  ! *** Local Variables
  real  :: CDRAG, CONC, DA, DENA, DENW, DISFRAC, DW, HE, HPU, HPI
  real  :: K_L, K_G, KV, RA, R_CONST, SCA, SCW, TKA, TKW
  real  :: U_AVG, USTR, V_AVG, VEL_AVG, VISCA, VISCW, VOLTERM

  ! *** HP:          Water depth (m)
  ! *** HPKCI:       Top layer water thickness (m)
  ! *** STCUV:       Flag to zero shear for diagonal cells
  ! *** UHE/VHE:     Unit discharge in the u and v directions (m2/s)
  ! *** HU/HV:       Flow depths at cell interface (m)
  ! *** TEM_KC:      Surface temperature (deg c)
  ! *** VEL_MAX:     Maximum velocity for lake condition (m/s)
  ! *** DEP_MIN:     Minimum depth for lake condition (m)
  ! *** SAL_KC:      Surface salinity concentration (g/kg)
  ! *** MOLW_I:      1./molecular weight (g/mole)**0.66667
  ! *** TATMT:       Atmospheric temperature (deg c)
  ! *** PATMT:       Atmospheric pressure (atm)
  ! *** WINDST:      Wind speed (m/s)
  ! *** MULT:        Volatilization rate adjustment factor (-)
  ! *** AIRCON:      Atmospheric concentration (ug/l)
  ! *** TCOEFF:      Mass transfer temperature coefficient (-)
  ! *** C_HE:        Henry's constant
  ! *** WCCON:       Constituent concentration    (g/m3;mg/l)
  ! *** PFTW:        Solids partitioning fraction (-)
  ! *** KL_OPT:      Liquid film Coefficient option

  HPI = 1./HP
  R_CONST = 8.206E-5

  ! *** Calculate depth averaged velocity
  U_AVG = STCUV*UHE/HU
  V_AVG = STCUV*VHE/HV
  VEL_AVG = SQRT( U_AVG*U_AVG + V_AVG*V_AVG )

  TKW = TEM_KC + 273.15                                            ! *** Temperature (Kelvin)
  TKA = TATMT + 273.15                                             ! *** Temperature (Kelvin)
  if( VEL_AVG > VEL_MAX .or. HP < DEP_MIN )then
    ! *** RIVER/STREAM CONDITIONS: WATER TURBULENCE CONTROLLED GAS TRANSFER
    HPU = 3.45*VEL_AVG**2.5

    if( HP < 0.61 )then
      ! *** OWENS FORMULA FOR K_L (K_A)
      K_L = 5.35*VEL_AVG**0.667 * HPI**1.85                         ! *** Compute reaeration rate (1/day)
      K_L = K_L * HP / 86400.                                      ! *** Multiply reaeration rate with depth to get transfer rate, then convert from m/day to m/s
    elseif( HP > HPU )then
      ! *** O'CONNOR-DOBBINS
      DW = 22.E-9*MOLW_I                                           ! *** Diffusivity of toxic in water (m^2/s)
      K_L = (DW*VEL_AVG)**0.5 * HPI**1.5                           ! *** Compute reaeration rate (1/s), because all units are already in m and s, no need to /86400
      K_L = K_L * HP                                               ! *** Multiply reaeration rate with depth to get transfer rate (m/s)
    else
      ! *** CHURCHILL
      K_L = 5.026*VEL_AVG**0.969 * HPI**1.673                       ! *** Compute reaeration rate (1/day)
      K_L = K_L * HP / 86400.                                      ! *** Multiply reaeration rate with depth to get transfer rate, then convert from m/day to m/s
    endif
    K_G = 0.0011574                                                ! *** Set gas transfer coefficient (m/s)  [100 m/day]
  else
    ! *** LAKE/QUIESCENT CONDITIONS: AIR TURBULENCE CONTROLLED GAS TRANSFER

    ! *** SCHMIDT NUMBER FOR WATER (dimensionless)
    VISCW = ViscosityW(TEM_KC)                                     ! *** Viscosity of water (kg/m/s)
    DENW = FUNDEN(SAL_KC,0.0,TEM_KC)                               ! *** kg/m^3

    DW = 22.E-9*MOLW_I                                             ! *** Diffusivity of toxic in water (m^2/s)
    SCW = VISCW/DENW/DW                                            ! *** Schmidt Number

    ! *** SCHMIDT NUMBER FOR AIR (dimensionless)
    
    RA = 287.05                                                    ! *** Specific Gas Constant (m^2/s^2/K) (J/kg/K)  (J = kg*m^2/s^2)
    DENA = PATMT*100./RA/TKA                                       ! *** Density of air (kg/m^3)
    VISCA = 1.458E-6*TKA**1.5/(TKA+110.4)                          ! *** Viscosity of air (kg/m/s)
    DA = 1.9E-4*MOLW_I                                             ! *** Diffusivity of toxic in air (m^2/s)
    SCA = VISCA/DENA/DA                                            ! *** Schmidt Number

    ! *** COMPUTE GAS AND LIQUID TRANSFER COEFFICIENTS (M/S)
    CDRAG = ( 6.1 + 0.63*WINDST )*1.0E-04
    CDRAG = AMAX1(CDRAG,0.0011)
    USTR = WINDST*CDRAG**0.5                                       ! *** Shear velocity (m/s)
    
    ! *** CALCULATE LIQUID FILM COEFFICIENT
    if( KL_OPT == 2 )then
      ! *** CALCULATE K_L AS PER MACKAY AND YEUN (1983)
      if( USTR < 0.3 )then
        K_L = 0.0144*USTR**2.2/SQRT(SCW) + 1E-6
      else
        K_L = 0.00341*USTR/SQRT(SCW) + 1E-6
      endif
      K_G = 0.0462*USTR/SCA**0.6667 + 1E-3                         ! *** .0462 = k^0.33/LAMDA
    elseif( KL_OPT == 1 )then                                      
      ! *** CALCULATE K_L AS PER O'CONNOR                          
      K_L = 0.0462*USTR*(DENA/DENW)**.5/SCW**0.6667                ! *** .0462 = k^0.33/LAMDA
      K_G = 0.0462*USTR/SCA**0.6667                                ! *** .0462 = k^0.33/LAMDA
    endif                                                          
    K_G = AMAX1(K_G,1.0E-06)                                       ! *** Gas mass transfer (m/s)
  endif                                                            
  K_L = AMAX1(K_L,1.0E-9)                                          ! *** Liquid mass transfer (m/s)                                                            
  HE = C_HE/R_CONST/TKA                                            ! *** Henry's term (dimensionless)
  KV = 1.0 / ( 1.0/K_L + 1.0/(K_G*HE) )                            ! *** Mass transfer rate at 20 degC (m/s)
  KV = KV*TCOEFF**(TEM_KC-20.)                                     ! *** Mass transfer rate at water temperature (m/s)

  CONC = WCCON
  
  ! *** PFTW - Particulate fraction 
  DISFRAC = 1./(1. + PFTW)                                         ! *** Dissolved fraction
  CONC = WCCON * DISFRAC
  
  ! *** Compute volatization flux term ( mg/(l*s) )
  VOLTERM = MULT*KV*HPKCI*( CONC - AIRCON*0.001 )

  CALVOLTERM = VOLTERM
  
END FUNCTION CALVOLTERM


