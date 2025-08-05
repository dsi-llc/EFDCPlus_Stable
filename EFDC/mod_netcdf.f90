! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! Write NetCDF UGrid format
! 2021-01-18, Nghiem Tien Lam
module mod_netcdf
! DELME - 2022-06-17 - Not updated for the water column fast settling classes
    
    use GLOBAL
    use julianmod
    use hydstrucmod
    use Variables_WQ
    use shellfishmod
    use fields
    use cyclone
    use convertwgs84
    use netcdf
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    use Variables_MPI_Drifter
    use Mod_Map_Write_NetCDF

    implicit none

    integer(4), parameter        :: ZERO = 0, ONE = 1, FOUR = 4, NC_VAR_CNT = 155, LEN_NAME = 20
    real(4), parameter           :: MISSING_VALUE = -999.
    
    type nc_var
        integer(2)     :: id = 0
        character(20)  :: name,units
        integer(2)     :: data_type,dim_type,layer_dim,comp_dim
        character(100) :: long_name
        character(100) :: standard_name
        integer(4),allocatable,dimension(:) :: idx
    end type

    type nc_dataset
        integer(4) :: id = 0
        integer(2) :: isopen = 0
        integer(4) :: node_dim, col_dim, row_dim, lcm_dim, cnr_dim
        integer(4) :: nodecnt, cellcnt, colcnt, rowcnt
        integer(4) :: nc_msh(9), nc_sig, nc_kl, nc_kc, nc_lyr
        integer(4) :: nc_cuv, nc_rssbc, nc_seddia
        integer(4) :: cell_dim, time_dim, name_dim
        integer(4) :: kc_dim, dye_dim, tox_dim, nalg_dim, nzoo_dim, nshf_dim
        integer(4) :: kb_dim, sed_dim, snd_dim, sxd_dim, nseds2_dim
        integer(4) :: nqctl_dim, nqwr_dim, time_idx
        integer(4) :: nc_time
        integer(4) :: idx(NC_VAR_CNT)
    end type
        
    
    type(nc_dataset),allocatable :: nc(:)
    character(10)                :: nc_datestamp

    contains

    function define_nc_vars(nc) result(ncvar)
    type(nc_dataset),intent(in) :: nc
    type(nc_var)                :: ncvar(NC_VAR_CNT)
    ! *** VARIABLE DECLARATIONS
    ncvar(1:NC_VAR_CNT) = (/&
        nc_var(  1, 'BELV',                 'm',            0, 3, 0,         0,           'Bottom elevation', 'sea_floor_depth_below_reference_ellipsoid'), &
        nc_var(  2, 'HP',                   'm',            0, 3, 0,         0,           'Water depth', 'depth'), &
        nc_var(  3, 'WSEL',                 'm',            0, 3, 0,         0,           'Water surface elevation', 'water_surface_height_above_reference_datum'), &
        nc_var(  4, 'U',                    'm/s',          0, 4, nc.KC_DIM, 0,           'Eastward water velocity', 'sea_water_x_velocity'), &
        nc_var(  5, 'V',                    'm/s',          0, 4, nc.KC_DIM, 0,           'Northward water velocity', 'sea_water_y_velocity'), &
        nc_var(  6, 'W',                    'm/s',          0, 4, nc.KC_DIM, 0,           'Upward water velocity', 'upward_sea_water_velocity'), &
        nc_var(  7, 'DELT',                 's',            0, 1, 0,         0,           'Time step', ''), &
        nc_var(  8, 'total_shear',          'N/m2',         0, 3, 0,         0,           'Total bed shear stress', ''), &
        nc_var(  9, 'cohesive_shear',       'N/m2',         0, 3, 0,         0,           'Cohesive grain size shear stress', ''), &
        nc_var( 10, 'noncohesive_shear',    'N/m2',         0, 3, 0,         0,           'Non-cohesive grain size shear stress', ''), &
        nc_var( 11, 'current_shear',        'N/m2',         0, 3, 0,         0,           'Current-induced bed shear stress', ''), &
        nc_var( 12, 'wave_shear',           'N/m2',         0, 3, 0,         0,           'Wave-induced bed shear stress', ''), &
        nc_var( 13, 'wind_x',               'm/s',          0, 3, 0,         0,           'Eastward wind velocity', 'eastward_wind'), &
        nc_var( 14, 'wind_y',               'm/s',          0, 3, 0,         0,           'Northward wind velocity', 'northward_wind'), &
        nc_var( 15, 'wave_height',          'm',            0, 3, 0,         0,           'Significant wave height', 'sea_surface_wave_significant_height'), &
        nc_var( 16, 'wave_direction',       'degree',       0, 3, 0,         0,           'Wave direction', 'sea_surface_wave_from_direction'), &
        nc_var( 17, 'wave_period',          's',            0, 3, 0,         0,           'Wave period', 'sea_surface_wave_mean_period'), &
        nc_var( 18, 'wave_dissipation',     'W/m2',         0, 3, 0,         0,           'Wave energy dissipation', ''), &
        nc_var( 19, 'wave_SXX',             'W/m2',         0, 3, 0,         0,           'Wave radian stress Sxx', 'sea_surface_wave_xx_radiation_stress'), &
        nc_var( 20, 'wave_SXY',             'W/m2',         0, 3, 0,         0,           'Wave radian stress Sxy', 'sea_surface_wave_xy_radiation_stress'), &
        nc_var( 21, 'wave_SYY',             'W/m2',         0, 3, 0,         0,           'Wave radian stress Syy', 'sea_surface_wave_yy_radiation_stress'), &
        nc_var( 22, 'salinity',             'ppt',          0, 4, nc.KC_DIM, 0,           'Salinity', 'sea_water_salinity'), &
        nc_var( 23, 'temperature',          'degC',         0, 4, nc.KC_DIM, 0,           'Water temperature', 'sea_water_temperature'), &
        nc_var( 24, 'bed_temperature',      'degC',         0, 3, 0,         0,           'Bed temperature', 'temperature_in_ground'), &
        nc_var( 25, 'evaporation',          'mm/s',         0, 3, 0,         0,           'Evaporation', 'surface_water_evaporation_flux'), &
        nc_var( 26, 'rainfall',             'mm/day',       0, 3, 0,         0,           'Rainfall', 'rainfall_rate'), &
        nc_var( 27, 'dye',                  'mg/L',         0, 6, nc.KC_DIM, nc.DYE_DIM,  'Dye', ''), &
        nc_var( 28, 'shellfish_larvae',     'mg/L',         0, 4, nc.KC_DIM, 0,           'Shellfish larvae', ''), &
        nc_var( 29, 'toxics',               'ug/L',         0, 6, nc.KC_DIM, nc.TOX_DIM,  'Toxics', ''), &
        nc_var( 30, 'bed_toxics',           'mg/m2',        0, 6, nc.KB_DIM, nc.TOX_DIM,  'Bed toxics', ''), &
        nc_var( 31, 'bed_top',              'layer',        1, 3, 0,         0,           'Sediment bed top layer index', ''), &
        nc_var( 32, 'bed_thickness',        'm',            0, 3, nc.KB_DIM, 0,           'Bed sediment thickness', ''), &
        nc_var( 33, 'bed_wet_density',      'kg/m3',        0, 3, nc.KB_DIM, 0,           'Sediment bed wet density', ''), &
        nc_var( 34, 'bed_porosity',         '1',            0, 3, nc.KB_DIM, 0,           'Bed sediment porosity', ''), &
        nc_var( 35, 'bed_cohesive',         'g/m2',         0, 6, nc.KB_DIM, nc.SED_DIM,  'Bed cohesive sediment', ''), &
        nc_var( 36, 'bed_noncohesive',      'g/m2',         0, 6, nc.KB_DIM, nc.SND_DIM,  'Bed non-cohesive sediment', ''), &
        nc_var( 37, 'cohesive_sediment',    'mg/L',         0, 6, nc.KC_DIM, nc.SED_DIM,  'Cohesive sediment', ''), &
        nc_var( 38, 'noncohesive_sediment', 'mg/L',         0, 6, nc.KC_DIM, nc.SND_DIM,  'Non-cohesive sediment', ''), &
        nc_var( 39, 'active_layer',         '1',            1, 4, nc.KB_DIM, 0,           'Active layer', ''), &
        nc_var( 40, 'bed_tau',              'dynes/cm2',    0, 3, 0,         0,           'Shear stress', ''), &
        nc_var( 41, 'bed_mass',             'g/cm2',        0, 4, nc.KB_DIM, 0,           'Sediment mass', ''), &
        nc_var( 42, 'bed_dry_density',      'g/cm3',        0, 4, nc.KB_DIM, 0,           'Dry bulk density', ''), &
        nc_var( 43, 'bed_mass_fraction',    '1',            0, 6, nc.KB_DIM, nc.SXD_DIM,  'Sediment mass fraction', ''), &
        nc_var( 44, 'bed_d50',              'microns',      0, 3, 0,         0,           'Average particle size of bed surface', ''), &
        nc_var( 45, 'bed_erosion_rate',     'g/cm2/s',      0, 5, 0,         nc.nseds2_dim,  'Total erosion rate in the cell', ''), &
        nc_var( 46, 'bed_deposition_rate',  'g/cm2/s',      0, 5, 0,         nc.nseds2_dim,  'Total deposition rate in the cell', ''), &
        nc_var( 47, 'bedload_mass',         'g/m2',         0, 5, 0,         nc.SXD_DIM,  'Bedload mass for each size class', ''), &
        nc_var( 48, 'bedload_x',            'g/s',          0, 5, 0,         nc.SXD_DIM,  'Bedload transport in x-direction', ''), &
        nc_var( 49, 'bedload_y',            'g/s',          0, 5, 0,         nc.SXD_DIM,  'Bedload transport in y-direction', ''), &
        nc_var( 50, 'bedload_toxics',       'mg/m2',        0, 5, 0,         nc.TOX_DIM,  'Bedload toxic concentration', ''), &

        nc_var( 51, 'EVAPSW',               'm3/s',         0, 3, 0,         0,           'Groundwater transpiration', ''), &
        nc_var( 52, 'EVAPGW',               'm3/s',         0, 3, 0,         0,           'Groundwater evaporation', ''), &
        nc_var( 53, 'gw_flow',              'm3/s',         0, 3, 0,         0,           'Groundwater inflow', ''), &
        nc_var( 54, 'gw_level',             'm',            0, 3, 0,         0,           'Groundwater table', ''), &
        nc_var( 55, 'ice_thick',            'm',            0, 3, 0,         0,           'Ice thickness', 'floating_ice_thickness'), &
        nc_var( 56, 'ice_temp',             'degC',         0, 3, 0,         0,           'Ice temperature', ''), &

        nc_var( 57, 'RPOC',                 'mg/L',         0, 4, nc.KC_DIM, 0,           'Refractory particulate organic carbon', ''), &
        nc_var( 58, 'LPOC',                 'mg/L',         0, 4, nc.KC_DIM, 0,           'Labile particulate organic carbon', ''), &
        nc_var( 59, 'DOC',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Dissolved organic carbon', ''), &
        nc_var( 60, 'RPOP',                 'mg/L',         0, 4, nc.KC_DIM, 0,           'Refractory particulate organic phosphorus', ''), &
        nc_var( 61, 'LPOP',                 'mg/L',         0, 4, nc.KC_DIM, 0,           'Labile particulate organic phosphorus', ''), &
        nc_var( 62, 'DOP',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Dissolved organic phosphorus', ''), &
        nc_var( 63, 'P4D',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Total phosphate', ''), &
        nc_var( 64, 'RPON',                 'mg/L',         0, 4, nc.KC_DIM, 0,           'Refractory particulate organic nitrogen', ''), &
        nc_var( 65, 'LPON',                 'mg/L',         0, 4, nc.KC_DIM, 0,           'Labile particulate organic nitrogen', ''), &
        nc_var( 66, 'DON',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Dissolved organic nitrogen', ''), &
        nc_var( 67, 'NHX',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Ammonia nitrogen', ''), &
        nc_var( 68, 'NOX',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Nitrate nitrogen', ''), &
        nc_var( 69, 'SUU',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Particulate biogenic silica', ''), &
        nc_var( 70, 'SAA',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Available dissolved silica', ''), &
        nc_var( 71, 'COD',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Chemical oxygen demand', ''), &
        nc_var( 72, 'DOX',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Dissolved oxygen', ''), &
        nc_var( 73, 'TAM',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Total active metal', ''), &
        nc_var( 74, 'FCB',                  'MPN/100ml',    0, 4, nc.KC_DIM, 0,           'Fecal coliform bacteria', ''), &
        nc_var( 75, 'CO2',                  'mg/L',         0, 4, nc.KC_DIM, 0,           'Carbon dioxide', ''), &
        nc_var( 76, 'ALG',                  'mg C L-1',     0, 6, nc.KC_DIM, nc.NALG_DIM, 'Phytoplankton', ''), &
        nc_var( 77, 'ZOO',                  'mg C L-1',     0, 6, nc.KC_DIM, nc.NZOO_DIM, 'Zooplankton', ''), &

        nc_var( 78, 'SMPON1',               'g/m3',         0, 3, 0,         0,           'Particulate organic nitrogen in G-class 1', ''), &
        nc_var( 79, 'SMPON2',               'g/m3',         0, 3, 0,         0,           'Particulate organic nitrogen in G-class 2', ''), &
        nc_var( 80, 'SMPON3',               'g/m3',         0, 3, 0,         0,           'Particulate organic nitrogen in G-class 3', ''), &
        nc_var( 81, 'SMPOP1',               'g/m3',         0, 3, 0,         0,           'Particulate organic phosphorus in G-class 1', ''), &
        nc_var( 82, 'SMPOP2',               'g/m3',         0, 3, 0,         0,           'Particulate organic phosphorus in G-class 2', ''), &
        nc_var( 83, 'SMPOP3',               'g/m3',         0, 3, 0,         0,           'Particulate organic phosphorus in G-class 3', ''), &
        nc_var( 84, 'SMPOC1',               'g/m3',         0, 3, 0,         0,           'Particulate organic carbon in G-class 1', ''), &
        nc_var( 85, 'SMPOC2',               'g/m3',         0, 3, 0,         0,           'Particulate organic carbon in G-class 2', ''), &
        nc_var( 86, 'SMPOC3',               'g/m3',         0, 3, 0,         0,           'Particulate organic carbon in G-class 3', ''), &
        nc_var( 87, 'SMDFN1',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from PON into G1', ''), &
        nc_var( 88, 'SMDFN2',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from PON into G2', ''), &
        nc_var( 89, 'SMDFN3',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from PON into G3', ''), &
        nc_var( 90, 'SMDFP1',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from POP into G1', ''), &
        nc_var( 91, 'SMDFP2',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from POP into G2', ''), &
        nc_var( 92, 'SMDFP3',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from POP into G3', ''), &
        nc_var( 93, 'SMDFC1',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from POC into G1', ''), &
        nc_var( 94, 'SMDFC2',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from POC into G2', ''), &
        nc_var( 95, 'SMDFC3',               'g/m2/day',     0, 3, 0,         0,           'Flux to sediment bed from POC into G3', ''), &
        nc_var( 96, 'SM1NH4',               'g/m3',         0, 3, 0,         0,           'Ammonium-nitrogen concentration in layer 1', ''), &
        nc_var( 97, 'SM2NH4',               'g/m3',         0, 3, 0,         0,           'Ammonium-nitrogen concentration in layer 2', ''), &
        nc_var( 98, 'SM1NO3',               'g/m3',         0, 3, 0,         0,           'Nitrate-nitrogen concentration in layer 1', ''), &
        nc_var( 99, 'SM2NO3',               'g/m3',         0, 3, 0,         0,           'Nitrate-nitrogen concentration in layer 2', ''), &
        nc_var(100, 'SM1PO4',               'g/m3',         0, 3, 0,         0,           'Orthophosphate as phosphorus concentration in layer 1', ''), &
        nc_var(101, 'SM2PO4',               'g/m3',         0, 3, 0,         0,           'Orthophosphate as phosphorus concentration in layer 2', ''), &
        nc_var(102, 'SM1H2S',               'g/m3',         0, 3, 0,         0,           'Sulfide concentration in layer 1', ''), &
        nc_var(103, 'SM2H2S',               'g/m3',         0, 3, 0,         0,           'Sulfide concentration in layer 2', ''), &
        nc_var(104, 'SM1SI',                'g/m3',         0, 3, 0,         0,           'Available dissolved silica in layer 1', ''), &
        nc_var(105, 'SM2SI',                'g/m3',         0, 3, 0,         0,           'Available dissolved silica in layer 2', ''), &
        nc_var(106, 'SMPSI',                'g/m3',         0, 3, 0,         0,           'Biogenic particulate silica', ''), &
        nc_var(107, 'SMBST',                'day',          0, 3, 0,         0,           'Accumulated benthic stress', ''), &
        nc_var(108, 'SMT',                  'degC',         0, 3, 0,         0,           'Sediment temperature', ''), &
        nc_var(109, 'SMCSOD',               'g/m2/day',     0, 3, 0,         0,           'Carbonaceous sediment oxygen demands', ''), &
        nc_var(110, 'SMNSOD',               'g/m2/day',     0, 3, 0,         0,           'Nitrogenous sediment oxygen demands', ''), &
        nc_var(111, 'WQBFNH4',              'g/m2/day',     0, 3, 0,         0,           'Ammonium flux', ''), &
        nc_var(112, 'WQBFNO3',              'g/m2/day',     0, 3, 0,         0,           'Nitrate flux', ''), &
        nc_var(113, 'WQBFO2',               'g/m2/day',     0, 3, 0,         0,           'Sediment oxygen flux', ''), &
        nc_var(114, 'WQBFCOD',              'g/m2/day',     0, 3, 0,         0,           'COD flux', ''), &
        nc_var(115, 'WQBFPO4D',             'g/m2/day',     0, 3, 0,         0,           'Phosphate flux', ''), &
        nc_var(116, 'WQBFSAD',              'g/m2/day',     0, 3, 0,         0,           'Silica flux', ''), &

        nc_var(117, 'WQRPS',                'g C m-2',      0, 3, 0,         0,           'Shoot carbon', ''), &
        nc_var(118, 'WQRPR',                'g C m-2',      0, 3, 0,         0,           'Root carbon', ''), &
        nc_var(119, 'WQRPE',                'g C m-2',      0, 3, 0,         0,           'Epiphyte carbon', ''), &
        nc_var(120, 'WQRPD',                'g C m-2',      0, 3, 0,         0,           'Detritus carbon', ''), &

        nc_var(121, 'HS_IN',                'm3/s',         0, 8, nc.KC_DIM, nc.NQCTL_DIM,'Control structure inflow', ''), &
        nc_var(122, 'HS_OUT',               'm3/s',         0, 8, nc.KC_DIM, nc.NQCTL_DIM,'Control structure outflow', ''), &
        nc_var(123, 'HS_STATE',             '1',            1, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure operation state', ''), &
        nc_var(124, 'HS_UNITS',             '1',            1, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure operation units', ''), &
        nc_var(125, 'HS_CURVE',             '1',            1, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure operation rating curve', ''), &
        nc_var(126, 'HS_FLOW',              'm3/s',         0, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure flow', ''), &
        nc_var(127, 'HS_HEIGHT',            'm',            0, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure opening height', ''), &
        nc_var(128, 'HS_WIDTH',             'm',            0, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure opening width', ''), &
        nc_var(129, 'HS_SILL',              'm',            0, 7, 0,         nc.NQCTL_DIM,'Hydraulic structure sill', ''), &

        nc_var(130, 'WR_STATE',             '1',            1, 7, 0,         nc.NQWR_DIM, 'Withdrawal/return operation state', ''), &
        nc_var(131, 'WR_FLOW',              'm3/s',         0, 7, 0,         nc.NQWR_DIM, 'Withdrawal/return flow', ''), &

        nc_var(132, 'SF_WCK',               'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish carbon weight per grid cell layer', ''), &
        nc_var(133, 'SF_WIN',               'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish carbon weight per individual', ''), &
        nc_var(134, 'SF_LIK',               'mm',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish length', ''), &
        nc_var(135, 'SF_NIK',               '#',            0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish individuals per layer per grid cell', ''), &
        nc_var(136, 'SF_DIK',               '#',            0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Number of death shellfish individuals', ''), &
        nc_var(137, 'SF_HIK',               '#',            0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish harvested individuals', ''), &
        nc_var(138, 'SF_HWK',               'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish harvested weight', ''), &
        nc_var(139, 'SF_NPK',               'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish net production', ''), &
        nc_var(140, 'SF_PRK',               'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish reproductive tissue production', ''), &
        nc_var(141, 'SF_CRBK',              'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish cumulative reproductive biomass', ''), &
        nc_var(142, 'SF_SPNK',              'gC',           0, 6, nc.KC_DIM, nc.NSHF_DIM, 'Shellfish spawn', ''), &
        nc_var(143, 'SF_WDC',               'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish dry meat weight per grid cell', ''), &
        nc_var(144, 'SF_NIC',               '#',            0, 5, 0,         nc.NSHF_DIM, 'Shellfish individuals per grid cell', ''), &
        nc_var(145, 'SF_FR',                'L/#/hour',     0, 5, 0,         nc.NSHF_DIM, 'Shellfish filtration rate', ''), &
        nc_var(146, 'SF_BMG',               '1/day',        0, 5, 0,         nc.NSHF_DIM, 'Shellfish basal metabolic rate', ''), &
        nc_var(147, 'SF_RESP',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish respiration', ''), &
        nc_var(148, 'SF_GRAZ',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish ingestion', ''), &
        nc_var(149, 'SF_DEAT',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish Death', ''), &
        nc_var(150, 'SF_FECA',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish feces', ''), &
        nc_var(151, 'SF_URIN',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish urine', ''), &
        nc_var(152, 'SF_RPOC',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish RPOC', ''), &
        nc_var(153, 'SF_LPOC',              'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish LPOC', ''), &
        nc_var(154, 'SF_DOC',               'gC',           0, 5, 0,         nc.NSHF_DIM, 'Shellfish DOC', ''), &
        nc_var(155, 'air_pressure',        'hPa',           0, 3, 0,         0,           'Barometric pressure', 'air_pressure') /)
    end function

    subroutine set_nc_flags
        integer :: i
        !do i = 1, 50
        !    IS_NC_OUT(i) = 1
        !enddo
        if( ISWAVE == 0 ) IS_NC_OUT(3) = 0                                                 ! *** Waves
        if( NWSER == 0 .and. WIND.IFLAG == 0 .and. ICYCLONE == 0 ) IS_NC_OUT(4) = 0        ! *** Winds
        if( NASER == 0 .and. PRESS.IFLAG == 0 .and. ICYCLONE == 0 ) IS_NC_OUT(5) = 0       ! *** Barometric pressure    
        if( IEVAP <= 1 .and. RAIN.IFLAG == 0 ) IS_NC_OUT(6) = 0                            ! *** Rainfall
        if( IEVAP <= 1 .and. EVAP.IFLAG == 0 ) IS_NC_OUT(7) = 0                            ! *** Evaporation
        if( ISGWIE == 0 .and. GWSP.IFLAG == 0 ) IS_NC_OUT(8) = 0                           ! *** Groundwater
        if( NQCTL == 0 ) IS_NC_OUT(9) = 0                                                  ! *** Hydraulic Structures
        if( NQWR == 0 ) IS_NC_OUT(10) = 0                                                  ! *** Withdrawal/return
        if( NQJPIJ == 0 ) IS_NC_OUT(11) = 0                                                ! *** Jet plume
        if( NPD == 0 ) IS_NC_OUT(12) = 0                                                   ! *** Lagrangian Particle Tracks
        if( ISTRAN(2) == 0) IS_NC_OUT(13) = 0                                              ! *** Water Temperature
        if( ISTRAN(2) == 0 .or. TEMBO == 0.) IS_NC_OUT(14) = 0                             ! *** Bed temperature
        if( ISTRAN(2) == 0 .or. ISICE < 3 ) IS_NC_OUT(15) = 0                              ! *** Ice thickness and temperature
        if( ISTRAN(1) == 0 ) IS_NC_OUT(16) = 0                                             ! *** Salinity
        if( ISTRAN(3) == 0 .or. NDYE == 0 ) IS_NC_OUT(17) = 0                              ! *** Dye
        if( ISTRAN(6) == 0 ) IS_NC_OUT(18) = 0                                             ! *** Cohesive Sediments
        if( ISTRAN(7) == 0 ) IS_NC_OUT(19) = 0                                             ! *** Non-cohesive Sediments
        if( ISTRAN(6) == 0 .and. ISTRAN(7) == 0 .and. BATHY.IFLAG == 0 ) IS_NC_OUT(20) = 0 ! *** Bed Elevation
        if( ISTRAN(6) == 0 ) IS_NC_OUT(21) = 0                                             ! *** Bed cohesive sediment
        if( ISTRAN(7) == 0 ) IS_NC_OUT(22) = 0                                             ! *** Bed noncohesive sediment
        if( ICALC_BL == 0 .or. NSND == 0 ) IS_NC_OUT(23) = 0                               ! *** Bed loads
        if( ISTRAN(5) == 0 ) IS_NC_OUT(24) = 0                                             ! *** Water Columns Toxics
        if( ISTRAN(5) == 0 ) IS_NC_OUT(25) = 0                                             ! *** Bed Toxics
        do i = 1, 19
            if( ISTRAN(8) == 0 .or. ISTRWQ(i) == 0 ) IS_NC_OUT(i+25) = 0                   ! *** Water quality
        enddo
        IS_NC_OUT(44) = 0                                                                   ! *** CO2 is not activated yet!
        if( ISTRAN(8) == 0 .or. NALGAE == 0 ) IS_NC_OUT(45) = 0                            ! *** Phytoplankton
        if( ISTRAN(8) == 0 .or. NZOOPL == 0 ) IS_NC_OUT(46) = 0                            ! *** Zooplankton
        if( ISTRAN(8) == 0 .or. ISRPEM == 0 ) IS_NC_OUT(47) = 0                            ! *** RPEM
        if( ISTRAN(8) == 0 .or. ISFFARM == 0 ) IS_NC_OUT(48) = 0                           ! *** Shellfish Farm
        if( ISTRAN(4) == 0) IS_NC_OUT(49) = 0                                              ! *** Shellfish larvae
        if( ISTRAN(8) == 0 .or. IWQBEN == 0 .or. ISSDBIN == 0) IS_NC_OUT(50) = 0           ! *** Sediment Diagenesis
    end subroutine
    
    function nc_create_file(filename) result(nc)
    character(*),intent(in) :: filename
    type(nc_var) :: ncvar(NC_VAR_CNT)
    type(nc_dataset) :: nc
    character(20) :: bdate,zonestr*3
    real(8)   :: xm(1),ym(1),xll(1),yll(1), sumdzc
    real(RKD), allocatable, dimension(:,:) :: xc,yc
    real(4), allocatable :: lon(:),lat(:),lonc(:),latc(:),sigma(:),kcsgz(:),dzsgz(:,:),array2d(:,:)
    integer, allocatable :: I_indices(:),J_indices(:),L_indices(:)
    integer, allocatable, dimension(:,:) :: num, nv
    integer :: status,colcnt,rowcnt,lyrcnt,nsxd,i,j,k,k1,l,c,r,m,dim_mesh,mesh(1),empty(0),id_crs
        
    ! *** ENTER DEFINE MODE
    status = check_err(NF90_CREATE(filename, NF90_CLOBBER+NF90_NETCDF4, nc.id),'nc_create');!NF90_HDF5
    if(status /= NF90_NOERR )then
        return
    endif

    BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)
    ! *** modified to handle global
    if( LSEDZLJ )then
        NSXD = NSEDS
    else
        NSXD = NSED + NSND
    endif

    nc.cellcnt = LA_Global - 1
    nc.colcnt = IMX_Global - IMN_Global + 1
    nc.rowcnt = JMX_Global - JMN_Global + 1
    
    allocate(xc(nc.rowcnt+1,nc.colcnt+1), yc(nc.rowcnt+1,nc.colcnt+1), num(nc.rowcnt+1,nc.colcnt+1), nv(4,nc.cellcnt))
    allocate(I_indices(nc.cellcnt),J_indices(nc.cellcnt),L_indices(nc.cellcnt))
    xc(:,:) = 0.
    yc(:,:) = 0.
    num(:,:) = 0

    call readcorn
    call utmpars

    !open(100,file = 'mesh.txt',status = 'unknown')
    !write(*,*) 'I = ',IMN_Global,IMX_Global
    !write(*,*) 'J = ',JMN_Global,JMX_Global
    !write(100,*) nc.colcnt,nc.rowcnt,LA_Global
    do L = 2,LA_Global
        c = sorted_loc_to_glob(L).global_i - IMN_Global + 1
        r = sorted_loc_to_glob(L).global_j - JMN_Global + 1
        I = c + IMN_Global - 1
        J = r + JMN_Global - 1
        xc(r,c) = xc(r,c) + XCR(1,j,i)
        yc(r,c) = yc(r,c) + YCR(1,j,i)
        num(r,c) = num(r,c) + 1
        xc(r+1,c) = xc(r+1,c) + XCR(2,j,i)
        yc(r+1,c) = yc(r+1,c) + YCR(2,j,i)
        num(r+1,c) = num(r+1,c) + 1
        xc(r+1,c+1) = xc(r+1,c+1) + XCR(3,j,i)
        yc(r+1,c+1) = yc(r+1,c+1) + YCR(3,j,i)
        num(r+1,c+1) = num(r+1,c+1) + 1
        xc(r,c+1) = xc(r,c+1) + XCR(4,j,i)
        yc(r,c+1) = yc(r,c+1) + YCR(4,j,i)
        num(r,c+1) = num(r,c+1) + 1    
        !write(100,'(2I4,8F12.3)') i,j,(XCR(k,j,i),YCR(k,j,i),k = 1,4)
        L_indices(L-1) = L
        I_indices(L-1) = I
        J_indices(L-1) = J
    enddo

    nc.nodecnt = 0
    do c = 1,nc.colcnt+1
        do r = 1,nc.rowcnt+1
            if( num(r,c) > 0 )then
                xm(1) = xc(r,c) / num(r,c)
                ym(1) = yc(r,c) / num(r,c)
                if( ISWGS84 == 0 .and. UTMZ /= 0 )then
                    call utmr_wgs84(xm,ym,xll,yll)
                    xc(r,c) = xll(1)
                    yc(r,c) = yll(1)
                else
                    xc(r,c) = xm(1)
                    yc(r,c) = ym(1)
                endif
                nc.nodecnt = nc.nodecnt + 1
                num(r,c) = nc.nodecnt
                !write(100,'(2I4,I6,2F12.3,2F15.10)') c,r,num(r,c),xm(1),ym(1),xll(1),yll(1)
            endif
        enddo
    enddo

    allocate(lon(nc.nodecnt),lat(nc.nodecnt),lonc(nc.cellcnt), latc(nc.cellcnt))
    do c = 1,nc.colcnt+1
        do r = 1,nc.rowcnt+1
            k = num(r,c)
            if( k > 0 )then
                lon(k) = xc(r,c)
                lat(k) = yc(r,c)
                !write(100,'(2I4,I6,2F15.10)') c,r,k,lon(k),lat(k)
            endif
        enddo
    enddo

    do L = 2,LA_Global
        I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
        J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
        k = L - 1
        nv(1, k) = num(j,i)
        nv(2, k) = num(j,i+1)
        nv(3, k) = num(j+1,i+1)
        nv(4, k) = num(j+1,i)
        lonc(k) = 0.
        latc(k) = 0.
        m = 0
        do c = 1,4
            r = nv(c, k)
            !if( r > 0 )then
                lonc(k) = lonc(k) + lon(r)
                latc(k) = latc(k) + lat(r)
                m = m + 1
            !endif
        enddo
        if( m > 0 )then
            lonc(k) = lonc(k) / m
            latc(k) = latc(k) / m
        endif
        !write(100,'(2I6,3I4,2F15.10,4I6)') k,L,I,J,m,lonc(k),latc(k), (nv(c, k), c = 1,4)
    enddo
    deallocate(xc, yc, num)
    !close(100)

    ! *** DEFINE DIMENSIONS
    !mesh(1) = 1
    status = check_err(nf90_def_dim(nc.id, 'time', NF90_UNLIMITED, nc.time_dim),'time_dim');
    status = check_err(nf90_def_dim(nc.id, 'Mesh2D_nNodes', nc.nodecnt, nc.node_dim),'node_dim');
    status = check_err(nf90_def_dim(nc.id, 'Mesh2D_nFaces', nc.cellcnt, nc.cell_dim),'cell_dim');
    status = check_err(nf90_def_dim(nc.id, 'FOUR', FOUR, nc.cnr_dim),'cnr_dim');
    status = check_err(nf90_def_dim(nc.id, 'ONE', ONE, dim_mesh),'mesh_dim');
    status = check_err(nf90_def_dim(nc.id, 'LCM', LCM_Global-1, nc.lcm_dim),'lcm_dim');

    status = check_err(nf90_def_dim(nc.id, 'IC', nc.colcnt, nc.col_dim),'col_dim');
    status = check_err(nf90_def_dim(nc.id, 'JC', nc.rowcnt, nc.row_dim),'row_dim');
    status = check_err(nf90_def_dim(nc.id, 'KC', KC, nc.kc_dim),'kc_dim');
    status = check_err(nf90_def_dim(nc.id, 'KB', KB, nc.kb_dim),'kb_dim');
    status = check_err(nf90_def_dim(nc.id, 'NDYE', NDYE, nc.dye_dim),'dye_dim');
    status = check_err(nf90_def_dim(nc.id, 'NTOX', NTOX, nc.tox_dim),'tox_dim');
    status = check_err(nf90_def_dim(nc.id, 'NSED', NSED, nc.sed_dim),'sed_dim');
    status = check_err(nf90_def_dim(nc.id, 'NSND', NSND, nc.snd_dim),'snd_dim');
    status = check_err(nf90_def_dim(nc.id, 'NSXD', NSXD, nc.sxd_dim),'sxd_dim');
    status = check_err(nf90_def_dim(nc.id, 'NSEDS2', NSEDS2, nc.nseds2_dim),'nseds2_dim');
    status = check_err(nf90_def_dim(nc.id, 'NALG', NALGAE, nc.nalg_dim),'nalg_dim');
    status = check_err(nf90_def_dim(nc.id, 'NZOO', NZOOPL, nc.nzoo_dim),'nzoo_dim');
    status = check_err(nf90_def_dim(nc.id, 'NSHF', NSF, nc.nshf_dim),'nshf_dim');    
    status = check_err(nf90_def_dim(nc.id, 'NQCTL', NQCTL, nc.nqctl_dim),'nqctl_dim');    
    status = check_err(nf90_def_dim(nc.id, 'NQWR', NQWR, nc.nqwr_dim),'nqwr_dim');        
    !status = check_err(nf90_def_dim(nc.id, 'LPT_NPD', NPD, nc.lpt_npd_dim),'lpt_npd_dim');
    !status = check_err(nf90_def_dim(nc.id, 'LPT_TIME', TIMECNT, nc.lpt_time_dim),'lpt_tim_dim');

    ! *** DEFINE VARIABLES
    status = check_err(nf90_def_var(ncid = nc.id,name = 'crs',xtype = NF90_INT,varid = id_crs),'def_crs')
    status = nf90_put_att(nc.id,id_crs,'grid_mapping_name','latitude_longitude')
    status = nf90_put_att(nc.id,id_crs,'longitude_of_prime_meridian',0.0)
    status = nf90_put_att(nc.id,id_crs,'semi_major_axis', 6378137.0 )
    status = nf90_put_att(nc.id,id_crs,'inverse_flattening',298.257223563)
    status = nf90_put_att(nc.id,id_crs,'epsg_code','EPSG:4326')
    
    status = check_err(nf90_def_var(nc.id,'Mesh2D',NF90_INT,dim_mesh, nc.nc_msh(1)),'def_mesh')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'long_name','mesh topology')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'cf_role','mesh_topology')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'topology_dimension', 2 )
    status = nf90_put_att(nc.id,nc.nc_msh(1),'node_coordinates','Mesh2D_node_x Mesh2D_node_y')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'face_coordinates','Mesh2D_face_x Mesh2D_face_y')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'face_node_connectivity','Mesh2D_face_nodes')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'layer_dimension','KC')
    !status = nf90_put_att(nc.id,nc.nc_msh(1),'node_dimension','NODE')
    !status = nf90_put_att(nc.id,nc.nc_msh(1),'face_dimension','CELL')
    !status = nf90_put_att(nc.id,nc.nc_msh(1),'max_face_nodes_dimension','FOUR')
    !status = nf90_put_att(nc.id,nc.nc_msh(1),'edge_node_connectivity','ne')

    status = check_err(nf90_def_var(nc.id,'Mesh2D_face_nodes',NF90_INT,(/ nc.cnr_dim, nc.cell_dim /),nc.nc_msh(2)),'def_nv')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(2), 1, 1, deflev)
    !status = nf90_put_att(nc.id,nc.nc_msh(2),'standard_name','face_node_connnectivity')
    status = nf90_put_att(nc.id,nc.nc_msh(2),'long_name','nodes surrounding element')
    status = nf90_put_att(nc.id,nc.nc_msh(2),'cf_role','face_node_connnectivity')
    status = nf90_put_att(nc.id,nc.nc_msh(2),'start_index', 1 )
    !status = nf90_put_att(nc.id,nc.nc_msh(2),'mesh','Mesh2D')
    status = nf90_put_att(nc.id,nc.nc_msh(2),'units','nondimensional')
    !status = nf90_put_att(nc.id,nc.nc_msh(2),'location','face')
    status = nf90_put_att(nc.id,nc.nc_msh(2), '_FillValue', ZERO)

    status = check_err(nf90_def_var(nc.id,'Mesh2D_node_x',NF90_FLOAT,(/ nc.node_dim /),nc.nc_msh(3)),'def_lon')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(3), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(3),'standard_name','longitude')
    status = nf90_put_att(nc.id,nc.nc_msh(3),'long_name','x-coordinate of mesh nodes')
    status = nf90_put_att(nc.id,nc.nc_msh(3),'units','degrees_east')
    status = nf90_put_att(nc.id,nc.nc_msh(3),'axis','X')
    !status = nf90_put_att(nc.id,nc.nc_msh(3),'mesh','Mesh2D')
    !status = nf90_put_att(nc.id,nc.nc_msh(3),'location','node')
    
    status = check_err(nf90_def_var(nc.id,'Mesh2D_node_y',NF90_FLOAT,(/ nc.node_dim /),nc.nc_msh(4)),'def_lat')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(4), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(4),'standard_name','latitude')
    status = nf90_put_att(nc.id,nc.nc_msh(4),'long_name','y-coordinate of mesh nodes')
    status = nf90_put_att(nc.id,nc.nc_msh(4),'units','degrees_north')
    status = nf90_put_att(nc.id,nc.nc_msh(4),'axis','Y')
    !status = nf90_put_att(nc.id,nc.nc_msh(4),'mesh','Mesh2D')
    !status = nf90_put_att(nc.id,nc.nc_msh(4),'location','node')

    status = check_err(nf90_def_var(nc.id,'Mesh2D_face_x',NF90_FLOAT,(/ nc.cell_dim /),nc.nc_msh(5)),'def_lonc')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(5), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(5),'standard_name','longitude')
    status = nf90_put_att(nc.id,nc.nc_msh(5),'long_name','x-coordinate of cell centroid')
    status = nf90_put_att(nc.id,nc.nc_msh(5),'units','degrees_east')
    status = nf90_put_att(nc.id,nc.nc_msh(5),'axis','X')
    !status = nf90_put_att(nc.id,nc.nc_msh(5),'mesh','Mesh2D')
    !status = nf90_put_att(nc.id,nc.nc_msh(5),'location','face')

    status = check_err(nf90_def_var(nc.id,'Mesh2D_face_y',NF90_FLOAT,(/ nc.cell_dim /),nc.nc_msh(6)),'def_latc')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(6), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(6),'standard_name','latitude')
    status = nf90_put_att(nc.id,nc.nc_msh(6),'long_name','y-coordinate of cell centroid')
    status = nf90_put_att(nc.id,nc.nc_msh(6),'units','degrees_north')
    status = nf90_put_att(nc.id,nc.nc_msh(6),'axis','Y')
    !status = nf90_put_att(nc.id,nc.nc_msh(6),'mesh','Mesh2D')
    !status = nf90_put_att(nc.id,nc.nc_msh(6),'location','face')

    status = check_err(nf90_def_var(nc.id,'cell_id',NF90_UINT,(/ nc.cell_dim /),nc.nc_msh(7)),'def_cell_L')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(7), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(7),'long_name','cell index')
    status = nf90_put_att(nc.id,nc.nc_msh(7),'cf_role','timeseries_id')
    status = nf90_put_att(nc.id,nc.nc_msh(7),'start_index', 2 )
    status = nf90_put_att(nc.id,nc.nc_msh(7),'mesh','Mesh2D')
    status = nf90_put_att(nc.id,nc.nc_msh(7),'location','face')
    
    if( IS_NC_OUT(1) > 0 )then
        status = check_err(nf90_def_var(nc.id,'col_id',NF90_UINT,(/ nc.cell_dim /),nc.nc_msh(8)),'def_cell_I')
        status = nf90_def_var_deflate(nc.id, nc.nc_msh(8), 1, 1, deflev)
        status = nf90_put_att(nc.id,nc.nc_msh(8),'long_name','cell column index')
        status = nf90_put_att(nc.id,nc.nc_msh(8),'start_index', 3 )

        status = check_err(nf90_def_var(nc.id,'row_id',NF90_UINT,(/ nc.cell_dim /),nc.nc_msh(9)),'def_cell_J')
        status = nf90_def_var_deflate(nc.id, nc.nc_msh(9), 1, 1, deflev)
        status = nf90_put_att(nc.id,nc.nc_msh(9),'long_name','cell row index')    
        status = nf90_put_att(nc.id,nc.nc_msh(9),'start_index', 3 )
    endif
    if( IGRIDV > 0 )then
        status = check_err(nf90_def_var(nc.id,'layers',NF90_BYTE,(/ nc.cell_dim /),nc.nc_kc),'def_kc')
        status = nf90_put_att(nc.id,nc.nc_kc,'long_name','number of vertical layers')
        status = nf90_put_att(nc.id, nc.nc_kc, 'coordinates','face_y face_x')
        status = nf90_put_att(nc.id, nc.nc_kc, 'mesh','Mesh2D')
        status = nf90_put_att(nc.id, nc.nc_kc, 'location','face')

        status = check_err(nf90_def_var(nc.id,'bottom_layer',NF90_BYTE,(/ nc.cell_dim /),nc.nc_kl),'def_kl')
        status = nf90_put_att(nc.id,nc.nc_kl,'long_name','index of bottom layers')
        status = nf90_put_att(nc.id, nc.nc_kl, 'coordinates','face_y face_x')
        status = nf90_put_att(nc.id, nc.nc_kl, 'mesh','Mesh2D')
        status = nf90_put_att(nc.id, nc.nc_kl, 'location','face')

        status = check_err(nf90_def_var(nc.id,'Mesh2D_layers',NF90_FLOAT,(/ nc.kc_dim, nc.cell_dim/), nc.nc_sig),'def_sigma')
        status = nf90_put_att(nc.id,nc.nc_sig,'standard_name','ocean_sigma_coordinate')
        status = nf90_put_att(nc.id,nc.nc_sig,'long_name','sigma at layer midpoints')
        status = nf90_put_att(nc.id,nc.nc_sig,'units','sigma_level')
        status = nf90_put_att(nc.id,nc.nc_sig,'positive','up')
        status = nf90_put_att(nc.id,nc.nc_sig,'formula_terms','sigma: Mesh2D_layers eta: WSEL depth: Bottom')
        status = nf90_put_att(nc.id, nc.nc_sig, 'coordinates','face_y face_x')

        status = check_err(nf90_def_var(nc.id,'layer_thickness',NF90_FLOAT,(/ nc.kc_dim, nc.cell_dim/), nc.nc_lyr),'def_layer')
        status = nf90_put_att(nc.id,nc.nc_lyr,'long_name','water layer thickness')
        status = nf90_put_att(nc.id,nc.nc_lyr,'units','1')
        status = nf90_put_att(nc.id, nc.nc_lyr, 'coordinates','face_y face_x')
        status = nf90_put_att(nc.id, nc.nc_lyr, 'mesh','Mesh2D')
        status = nf90_put_att(nc.id, nc.nc_lyr, 'location','face')
    else
        status = check_err(nf90_def_var(nc.id,'Mesh2D_layers',NF90_FLOAT,(/ nc.kc_dim /), nc.nc_sig),'def_sigma')
        status = nf90_put_att(nc.id,nc.nc_sig,'standard_name','ocean_sigma_coordinate')
        status = nf90_put_att(nc.id,nc.nc_sig,'long_name','sigma at layer midpoints')
        status = nf90_put_att(nc.id,nc.nc_sig,'units','sigma_level')
        status = nf90_put_att(nc.id,nc.nc_sig,'positive','up')
        status = nf90_put_att(nc.id,nc.nc_sig,'formula_terms','sigma: Mesh2D_layers eta: WSEL depth: BELV')
        !status = nf90_put_att(nc.id,nc.nc_sig,'_CoordinateZisPositive','up')
        !status = nf90_put_att(nc.id,nc.nc_sig,'_CoordinateTransformType','Vertical')
        !status = nf90_put_att(nc.id,nc.nc_sig,'_CoordinateAxisType','GeoZ')
        !status = nf90_put_att(nc.id,nc.nc_sig,'_CoordinateAxes','sigma')

        status = check_err(nf90_def_var(nc.id,'layer_thickness',NF90_FLOAT,(/ nc.kc_dim /), nc.nc_lyr),'def_layer')
        status = nf90_put_att(nc.id,nc.nc_lyr,'long_name','water layer thickness')
        status = nf90_put_att(nc.id,nc.nc_lyr,'units','1')
    endif

    status = check_err(nf90_def_var(nc.id,'CUV',NF90_FLOAT,(/ nc.cnr_dim, nc.cell_dim/), nc.nc_cuv),'def_CUV')
    status = nf90_put_att(nc.id, nc.nc_cuv,'long_name','cell rotation')

    status = check_err(nf90_def_var(nc.id,'RSSBC',NF90_FLOAT,(/ nc.cnr_dim, nc.lcm_dim/), nc.nc_rssbc),'def_RSSBC')
    status = nf90_put_att(nc.id,  nc.nc_rssbc,'long_name','velocity edge weighting')

    status = check_err(nf90_def_var(nc.id,'SEDDIA',NF90_FLOAT,(/ nc.sxd_dim/), nc.nc_seddia),'def_SEDDIA')
    status = nf90_put_att(nc.id, nc.nc_seddia,'long_name','sediment grain size')
    status = nf90_put_att(nc.id, nc.nc_seddia,'units','micrometers')

    ncvar = define_nc_vars(nc)
    !call set_nc_flags
    call def_time_vars(nc, 0, ncvar)
    
    ! *** ASSIGN GLOBAL ATTRIBUTES
    write(zonestr,'(I3)') abs(utmz)
    status = check_err(nf90_put_att(nc.id,nf90_global,'Conventions','CF-1.11, UGRID-1.0'),'att_cf')
    status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'base_date',BDATE),'att_base')
    status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'Project',trim(proj)),'att_prj')
    if( utmz > 0 )then
        status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'utm_zone','UTM Zone '//trim(zonestr)//' Northern Hemisphere'),'')
    else
        status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'utm_zone','UTM Zone '//trim(zonestr)//' Southern Hemisphere'),'')
    endif
    status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'institution','DSI, LLC'),'att_dsi')

    ! *** LEAVE DEFINE MODE
    status = check_err(nf90_enddef(nc.id),'end_def');

    status = check_err(nf90_put_var(nc.id, nc.nc_msh(2), nv),'nv')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(3), lon),'node_x')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(4), lat),'node_y')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(5), lonc),'face_x')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(6), latc),'face_y')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(7), L_indices),'L')

    if( IS_NC_OUT(1) > 0 )then
        status = check_err(nf90_put_var(nc.id, nc.nc_msh(8), I_indices),'I')
        status = check_err(nf90_put_var(nc.id, nc.nc_msh(9), J_indices),'J')
    endif
    
    deallocate(lon, lat, lonc, latc, nv, L_indices, I_indices, J_indices)

    if( IGRIDV > 0 )then ! Sigma-Zed level
        allocate(kcsgz(nc.cellcnt),dzsgz(KC,nc.cellcnt),array2d(KC,nc.cellcnt))
        kcsgz = KMINV
        dzsgz = MISSING_VALUE
        array2d = MISSING_VALUE
        do L = 2,LA_Global ! *** Modified to handle global
            I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
            J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
            k1 = KSZ_Global(L)
            kcsgz(L-1) = KC - k1 + 1    ! Number of vertical layers
            sumdzc = 0.0
            do k = k1,KC
                sumdzc = sumdzc + DZC(L,K)
                dzsgz(k,L-1) = (sumdzc - 0.5*DZC(L,K))
                array2d(k,L-1) = DZC(L,K)
            enddo
        enddo
        status = check_err(nf90_put_var(nc.id, nc.nc_kl, real(KSZ_Global(2:LA_Global),4)),'kl')
        status = check_err(nf90_put_var(nc.id, nc.nc_kc, kcsgz),'kc')
        status = check_err(nf90_put_var(nc.id, nc.nc_sig, dzsgz),'dz')
        status = check_err(nf90_put_var(nc.id, nc.nc_lyr, array2d),'lyr')
        deallocate(kcsgz,dzsgz,array2d)
    else  ! *** Standard sigma stretched level
        allocate(sigma(KC))
        sumdzc = 0.0
        do k = 1,KC
            sumdzc = sumdzc + DZCK(K)
            sigma(k) = sumdzc - 0.5*DZCK(K)
        enddo
        status = check_err(nf90_put_var(nc.id, nc.nc_sig, sigma),'sig')
        status = check_err(nf90_put_var(nc.id, nc.nc_lyr, real(DZCK(1:KC),4)),'lyr')
        deallocate(sigma)
    endif

    allocate(array2d(4,LCM_Global))
    do L = 2,LA_GLobal
        array2d(1,L-1) = CUE_Global(L)
        array2d(2,L-1) = CVE_Global(L)
        array2d(3,L-1) = CUN_Global(L)
        array2d(4,L-1) = CVN_Global(L)
    enddo
    status = nf90_put_var(nc.id, nc.nc_cuv, array2d(:,1:nc.cellcnt), (/ 1, 1 /), (/ 4, nc.cellcnt /))

    do L = 2,LCM_Global
        array2d(1,L-1) = RSSBCE_GLobal(L)
        array2d(2,L-1) = RSSBCW_GLobal(L)
        array2d(3,L-1) = RSSBCN_GLobal(L)
        array2d(4,L-1) = RSSBCS_GLobal(L)
    enddo
    status = nf90_put_var(nc.id, nc.nc_rssbc, array2d, (/ 1, 1 /), (/ 4, LCM_Global-1 /))

    if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
        status = nc_write_array1d(nc, nc.nc_seddia, SEDDIA, 1, NSXD, 1., 'SEDDIA')
    endif

    !status = nc_write_cell_1d(nc.id, ix, nc.nc_idx(1), BELV_Global, 1., 'Zb')

    deallocate(array2d)
    end function
    
    function nc_create_hf_file(is,filename) result(nc)
    type(nc_dataset) :: nc
    integer,intent(in) :: is
    character(*),intent(in) :: filename
    type(nc_var) :: ncvar(NC_VAR_CNT)
    character(20) :: bdate,zonestr*3
    character(LEN_NAME) :: location
    real(8)   :: xm(1),ym(1),xll(1),yll(1), sumdzc
    real(4), allocatable :: lon(:),lat(:),array2d(:,:)
    integer, allocatable :: I_indices(:),J_indices(:),L_indices(:)
    integer :: status,pntcnt,nsxd,i,j,k,k1,l,c,r,m,empty(0)
    
    ! *** ENTER DEFINE MODE
    status = check_err(NF90_CREATE(filename, NF90_CLOBBER+NF90_NETCDF4, nc.id),'nc_create');!NF90_HDF5
    if(status /= NF90_NOERR )then
        return
    endif

    ! *** modified to handle global
    BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)
    if( LSEDZLJ )then
        NSXD = NSEDS
    else
        NSXD = NSED + NSND
    endif

    pntcnt = NPNT(IS)

    ! *** DEFINE DIMENSIONS
    status = check_err(nf90_def_dim(nc.id,'TIME', NF90_UNLIMITED, nc.time_dim),'time_dim');
    status = check_err(nf90_def_dim(nc.id,'LOCATION', pntcnt, nc.cell_dim),'pnt_dim');
    status = check_err(nf90_def_dim(nc.id,'NAME', LEN_NAME, nc.name_dim),'name_dim');        

    status = check_err(nf90_def_dim(nc.id,'KC', KC, nc.kc_dim),'kc_dim');
    status = check_err(nf90_def_dim(nc.id,'KB', KB, nc.kb_dim),'kb_dim');
    status = check_err(nf90_def_dim(nc.id,'NDYE', NDYE, nc.dye_dim),'dye_dim');
    status = check_err(nf90_def_dim(nc.id,'NTOX', NTOX, nc.tox_dim),'tox_dim');
    status = check_err(nf90_def_dim(nc.id,'NSED', NSED, nc.sed_dim),'sed_dim');
    status = check_err(nf90_def_dim(nc.id,'NSND', NSND, nc.snd_dim),'snd_dim');
    status = check_err(nf90_def_dim(nc.id,'NSXD', NSXD, nc.sxd_dim),'sxd_dim');
    status = check_err(nf90_def_dim(nc.id,'NSEDS2', NSEDS2, nc.nseds2_dim),'nseds2_dim');
    status = check_err(nf90_def_dim(nc.id,'NALG', NALGAE, nc.nalg_dim),'nalg_dim');
    status = check_err(nf90_def_dim(nc.id,'NZOO', NZOOPL, nc.nzoo_dim),'nzoo_dim');
    status = check_err(nf90_def_dim(nc.id, 'NSHF', NSF, nc.nshf_dim),'nshf_dim');    
    status = check_err(nf90_def_dim(nc.id,'NQCTL', NQCTL, nc.nqctl_dim),'nqctl_dim');    
    status = check_err(nf90_def_dim(nc.id,'NQWR', NQWR, nc.nqwr_dim),'nqwr_dim');        

    ! *** DEFINE VARIABLES
    status = check_err(nf90_def_var(nc.id,'cell',NF90_UINT,(/ nc.cell_dim /),nc.nc_msh(1)),'def_cell_L')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(1), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(1),'long_name','cell index')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'cf_role','timeseries_id')
    status = nf90_put_att(nc.id,nc.nc_msh(1),'start_index',(/ 2 /))
    
    status = check_err(nf90_def_var(nc.id,'col',NF90_UINT,(/ nc.cell_dim /),nc.nc_msh(2)),'def_cell_I')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(2), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(2),'long_name','cell column index')
    status = nf90_put_att(nc.id,nc.nc_msh(2),'start_index',(/ 3 /))

    status = check_err(nf90_def_var(nc.id,'row',NF90_UINT,(/ nc.cell_dim /),nc.nc_msh(3)),'def_cell_J')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(3), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(3),'long_name','cell row index')
    status = nf90_put_att(nc.id,nc.nc_msh(3),'start_index',(/ 3 /))
        
    status = check_err(nf90_def_var(nc.id,'lon',NF90_FLOAT,(/ nc.cell_dim /),nc.nc_msh(4)),'def_cell_x')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(4), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(4),'standard_name','longitude')
    status = nf90_put_att(nc.id,nc.nc_msh(4),'long_name','longitude')
    status = nf90_put_att(nc.id,nc.nc_msh(4),'units','degrees_east')

    status = check_err(nf90_def_var(nc.id,'lat',NF90_FLOAT,(/ nc.cell_dim /),nc.nc_msh(5)),'def_cell_y')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(5), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(5),'standard_name','latitude')
    status = nf90_put_att(nc.id,nc.nc_msh(5),'long_name','latitude')
    status = nf90_put_att(nc.id,nc.nc_msh(5),'units','degrees_north')

    status = check_err(nf90_def_var(nc.id,'name',NF90_CHAR,(/ nc.name_dim, nc.cell_dim /),nc.nc_msh(6)),'def_lat')
    status = nf90_def_var_deflate(nc.id, nc.nc_msh(6), 1, 1, deflev)
    status = nf90_put_att(nc.id,nc.nc_msh(6),'standard_name','name')
    status = nf90_put_att(nc.id,nc.nc_msh(6),'long_name','name')
        
    ncvar = define_nc_vars(nc)
    !call set_nc_flags
    call def_time_vars(nc, is, ncvar)

    ! *** ASSIGN GLOBAL ATTRIBUTES
    write(zonestr,'(I3)') abs(utmz)
    status = check_err(nf90_put_att(nc.id,nf90_global,'Conventions','CF-1.11'),'att_cf')
    status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'Base_date',BDATE),'att_base')
    status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'Project',trim(proj)),'att_prj')
    if( utmz > 0 )then
        status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'utm_zone','UTM Zone '//trim(zonestr)//' Northern Hemisphere'),'')
    else
        status = check_err(nf90_put_att(nc.id,NF90_GLOBAL,'utm_zone','UTM Zone '//trim(zonestr)//' Southern Hemisphere'),'')
    endif

    ! *** LEAVE DEFINE MODE
    status = check_err(nf90_enddef(nc.id),'end_def');
 
    allocate(lon(pntcnt),lat(pntcnt),I_indices(pntcnt),J_indices(pntcnt),L_indices(pntcnt))
    lon = MISSING_VALUE
    lat = MISSING_VALUE
    if( IJHFRE(IS) == 0 )then
        do i = 1, NPNT(IS)
            xm(1) = HFREGRP(IS).XCEL(i)
            ym(1) = HFREGRP(IS).YCEL(i)
            if( ISWGS84 == 0 .and. UTMZ /= 0 )then
                call utmr_wgs84(xm, ym, xll, yll)
                lon(i) = xll(1)
                lat(i) = yll(1)
            else
                lon(i) = xm(1)
                lat(i) = ym(1)
            endif
        enddo
    endif
    
    do k = 1, NPNT(IS)
        i = HFREGRP(is).ICEL(k)
        j = HFREGRP(is).JCEL(k)
        I_indices(k) = i
        J_indices(k) = j
        L_indices(k) = LIJ_Global(i,j)
    enddo
    
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(1), L_indices),'L')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(2), I_indices),'I')    !HFREGRP(is).ICEL
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(3), J_indices),'J')    !HFREGRP(is).JCEL
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(4), lon),'node_x')
    status = check_err(nf90_put_var(nc.id, nc.nc_msh(5), lat),'lat')
    do k = 1, NPNT(IS)
        location = trim(HFREGRP(is).NAME(k))    !//CHAR(0)
        status = check_err(nf90_put_var(nc.id, nc.nc_msh(6), location, (/ 1, k /), (/ LEN_NAME, 1 /)),'name')
    enddo

    deallocate(lon, lat, I_indices, J_indices, L_indices)
    end function 
    
    subroutine def_time_vars(nc, ix, ncvar)
    type(nc_dataset),intent(inout) :: nc
    type(nc_var),intent(in) :: ncvar(:)
    integer,intent(in) :: ix
    integer :: i, status 
    character(20) :: bdate
    
    BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)

    status = check_err(nf90_def_var(nc.id,'time',NF90_DOUBLE, nc.time_dim, nc.nc_time),'def_time')
    status = nf90_put_att(nc.id,nc.nc_time,'standard_name','time')
    status = nf90_put_att(nc.id,nc.nc_time,'long_name','time')
    status = nf90_put_att(nc.id,nc.nc_time,'units','days since '//trim(bdate))
    status = nf90_put_att(nc.id,nc.nc_time,'base_date',trim(bdate))
    status = nf90_put_att(nc.id,nc.nc_time,'calendar','julian')
    status = nf90_put_att(nc.id,nc.nc_time,'axis','T')
    status = nf90_put_att(nc.id,nc.nc_time,'format','modified julian day (MJD)')
    !status = nf90_put_att(nc.id,nc.nc_time,'time_zone','UTC')

    !call define_vars()
    do i = 1,7
        call def_var(nc, ncvar, i)
    enddo

    if( ISSPH(8) >= 1 )then
        if( IS_NC_OUT(2) > 0 )then
            call def_var(nc, ncvar, 8)    ! *** TOTAL BED SHEAR STRESS

            if( ISBEDSTR == 1 .and. .not. LSEDZLJ )then
                call def_var(nc, ncvar, 9)    ! *** Stress from non-cohesive components
            endif

            call def_var(nc, ncvar, 12)    ! *** Bed Shear due to Waves Only
            call def_var(nc, ncvar, 11)    ! *** Shear due to Current Only
            if( ISWAVE >= 3 )then
                do i = 15,17
                    call def_var(nc, ncvar, i)
                enddo
                if( ISWAVE == 4 )then
                    do i = 18,21
                        call def_var(nc, ncvar, i)
                    enddo
                endif
            endif
        endif

        if( IS_NC_OUT(16) > 0 ) call def_var(nc, ncvar, 22) ! *** Salinity

        if( IS_NC_OUT(13) > 0 ) call def_var(nc, ncvar, 23) ! *** Temperature
        if( IS_NC_OUT(14) > 0 ) call def_var(nc, ncvar, 24) ! *** Bed temperature
        if( IS_NC_OUT(7) > 0 )  call def_var(nc, ncvar, 25) ! *** Evaporation
        if( IS_NC_OUT(6) > 0 )  call def_var(nc, ncvar, 26) ! *** Rainfall
        if( IS_NC_OUT(17) > 0 ) call def_var(nc, ncvar, 27) ! *** Dye
        if( IS_NC_OUT(49) > 0 ) call def_var(nc, ncvar, 28) ! *** Shellfish larvae

        if( IS_NC_OUT(24) > 0 ) call def_var(nc, ncvar, 29) ! *** Toxics
        if( IS_NC_OUT(25) > 0 ) call def_var(nc, ncvar, 30) ! *** Bed toxics

        if( IS_NC_OUT(21) > 0 .or. IS_NC_OUT(22) > 0 )then
            call def_var(nc, ncvar, 31)     ! *** WRITE THE TOP LAYER INDEX
            do i = 32,34
                call def_var(nc, ncvar, i)  ! *** Bed thickness, wet_density, porosity
            enddo
        endif
        if( IS_NC_OUT(18) > 0 ) call def_var(nc, ncvar, 37)     ! *** Cohesive sediment concentration
        if( IS_NC_OUT(19) > 0 ) call def_var(nc, ncvar, 38)     ! *** Noncohesive sediment concentration
        if( IS_NC_OUT(21) > 0 ) call def_var(nc, ncvar, 35)     ! *** Bed cohesive sediment
        if( IS_NC_OUT(22) > 0 ) call def_var(nc, ncvar, 36)     ! *** Bed noncohesive sediment
        if( IS_NC_OUT(23) > 0 )then
            call def_var(nc, ncvar, 48)     ! *** Bedload_x
            call def_var(nc, ncvar, 49)     ! *** Bedload_y
        endif
        
        if( IS_NC_OUT(8) > 0 )then
            do i = 51,54
                call def_var(nc, ncvar, i)  ! *** Groundwater
            enddo
        endif

        if( IS_NC_OUT(15) > 0 )then
            call def_var(nc, ncvar, 55)     ! *** Ice thickness
            call def_var(nc, ncvar, 56)     ! *** Ice temperature
        endif
    endif

    if( LSEDZLJ )then
        call def_var(nc, ncvar, 39)     ! *** This is used to write out the active layers instead of using the LAYERACTIVE array.

        ! *** REAL*4  - Global arrays use reversed indicies for MPI mapping.  Writes the same order as OMP original
        do i = 40,46
            call def_var(nc, ncvar, i)
        enddo

        if( ICALC_BL > 0 )then
            do i = 47,49
                call def_var(nc, ncvar, i)
            enddo
        endif

        !if( IS_NC_OUT(25) > 0 ) call def_var(nc, ncvar, 30) ! *** Bed toxics
        if( IS_NC_OUT(25) > 0 .and. ICALC_BL > 0 ) call def_var(nc, ncvar, 50) ! *** Bedload toxic concentration (mg/m^2)

    else
        if( ISBEXP >= 1 .and. KB > 1  .and. ISSPH(8) < 1 )then
            if( IS_NC_OUT(21) > 0 .or. IS_NC_OUT(22) > 0 )then
                do i = 31,34
                    call def_var(nc, ncvar, i)  ! *** Bed top, thickness, wet_density, porosity
                enddo
            endif
            if( IS_NC_OUT(21) > 0 ) call def_var(nc, ncvar, 35) ! *** Bed cohesive
            if( IS_NC_OUT(22) > 0 ) call def_var(nc, ncvar, 36) ! *** Bed noncohesive
            if( IS_NC_OUT(25) > 0 ) call def_var(nc, ncvar, 30) ! *** Bed toxics
        endif
    endif

    do i = 57,74
        if( IS_NC_OUT(i-31) > 0 ) call def_var(nc, ncvar, i)  ! *** Water quality
    enddo
    if( IS_NC_OUT(45) > 0 ) call def_var(nc, ncvar, 76)     ! *** Algal groups
    if( IS_NC_OUT(46) > 0 ) call def_var(nc, ncvar, 77)     ! *** Zooplankton groups
    if( IS_NC_OUT(50) > 0 )then
        do i = 78,116
            call def_var(nc, ncvar, i)                      ! *** Sediment Diagenesis
        enddo
    endif
    if( IS_NC_OUT(47) > 0 )then
        do i = 117,120
            call def_var(nc, ncvar, i)                      ! *** RPEM
        enddo
    endif
    if( IS_NC_OUT(48) > 0 )then
        do i = 132,154
            call def_var(nc, ncvar, i)                      ! *** Shellfish Farm
        enddo
    endif

    if( IS_NC_OUT(9) > 0 )then
        call def_var(nc, ncvar, 121)                        ! *** Hydraulic Structures
        call def_var(nc, ncvar, 122)
        if(NQCTLSER > 0 .or. NQCTRULES > 0 )then
            do i = 123,129
                if(i /= 124) call def_var(nc, ncvar, i)
            enddo
        endif
    endif
    if( IS_NC_OUT(10) > 0 )then
        call def_var(nc, ncvar, 130)                        ! *** Withdrawal/return
        call def_var(nc, ncvar, 131)
    endif

    ! Global variables for WNDVELE, WNDVELN, PATMT are not yet available
    if( IS_NC_OUT(4) > 0 )then
        do i = 13, 14
            call def_var(nc, ncvar, i)                      ! *** Winds
        enddo
    endif
    if( IS_NC_OUT(5) > 0 )then
        call def_var(nc, ncvar, 155)                        ! *** Barometric pressure
    endif
    
    end subroutine 
    

    subroutine def_var(nc, ncvar, i)
    type(nc_dataset),intent(inout) :: nc
    type(nc_var),intent(in) :: ncvar(:)
    integer,intent(in) ::   i
    integer :: id, xtype, lyr_dim, cmp_dim, status

    lyr_dim = ncvar(i).layer_dim
    cmp_dim = ncvar(i).comp_dim
    
    if( ncvar(i).data_type == 1 )then
        xtype = NF90_SHORT
    else
        xtype = NF90_FLOAT
    endif

    !if( ncvar(i).id > 0 )then
        select case (ncvar(i).dim_type)
        case (1)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ nc.TIME_DIM /), id)
        case (2)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ nc.CELL_DIM /), id)
        case (3)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ nc.CELL_DIM, nc.TIME_DIM /), id)
        case (4)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ nc.CELL_DIM, lyr_dim, nc.TIME_DIM /), id)
        case (5)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ nc.CELL_DIM, cmp_dim, nc.TIME_DIM /), id)
        case (6)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ nc.CELL_DIM, lyr_dim, cmp_dim, nc.TIME_DIM /), id)
        case (7)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ cmp_dim, nc.TIME_DIM /), id)
        case (8)
            status = nf90_def_var(nc.id, ncvar(i).name, xtype, (/ lyr_dim, cmp_dim, nc.TIME_DIM /), id)
        end select
        status = nf90_def_var_deflate(nc.id, id, 1, 1, deflev)
        if(len_trim(ncvar(i).standard_name) > 0 )then
            status = nf90_put_att(nc.id, id, 'standard_name', ncvar(i).standard_name)
        endif
        status = nf90_put_att(nc.id, id, 'long_name', ncvar(i).long_name)
        status = nf90_put_att(nc.id, id, 'units', ncvar(i).units)
        select case (ncvar(i).dim_type)
        case (1)
            status = nf90_put_att(nc.id, id, 'coordinates','time')
        case (2)
            status = nf90_put_att(nc.id, id, 'coordinates','Mesh2D_face_y Mesh2D_face_x')
        case (3)
            status = nf90_put_att(nc.id, id, 'coordinates','time Mesh2D_face_y Mesh2D_face_x')
        case (4)
            status = nf90_put_att(nc.id, id, 'coordinates','time Mesh2D_layers Mesh2D_face_y Mesh2D_face_x')
        case (5)
            !status = nf90_put_att(nc.id, id, 'coordinates','time class Mesh2D_face_y Mesh2D_face_x')
            status = nf90_put_att(nc.id, id, 'coordinates','time Mesh2D_face_y Mesh2D_face_x')
        case (6)
            !status = nf90_put_att(nc.id, id, 'coordinates','time class Mesh2D_layers Mesh2D_face_y Mesh2D_face_x')
            status = nf90_put_att(nc.id, id, 'coordinates','time Mesh2D_layers Mesh2D_face_y Mesh2D_face_x')
        case (7)
            !status = nf90_put_att(nc.id, id, 'coordinates','time class')
            status = nf90_put_att(nc.id, id, 'coordinates','time')
        case (8)
            !status = nf90_put_att(nc.id, id, 'coordinates','time class sigma')
            status = nf90_put_att(nc.id, id, 'coordinates','time Mesh2D_layers')
        end select        
        if( ncvar(i).dim_type >= 2 .and. ncvar(i).dim_type <= 6 )then
            status = nf90_put_att(nc.id, id, 'mesh', 'Mesh2D')
            status = nf90_put_att(nc.id, id, 'location', 'face')
            status = nf90_put_att(nc.id, id, 'grid_mapping','crs')
        endif
        status = nf90_put_att(nc.id, id, 'type', ncvar(i).dim_type)
        if( ncvar(i).data_type == 1 )then
            status = nf90_put_att(nc.id, id, '_FillValue', ZERO)
        else
            status = nf90_put_att(nc.id, id, '_FillValue', MISSING_VALUE)
        endif
        nc.idx(i) = id
    !endif
    end subroutine

    integer function nc_open_file(nc_id, filename)
        character(*),intent(in) :: filename
        integer :: nc_id, status
    
        status = check_err(NF90_OPEN(filename, NF90_WRITE + NF90_SHARE, nc_id),'nc_open');
        if(status /= NF90_NOERR )then
            print *, nf90_strerror(status)
        endif
        nc_open_file = status
    end function
    
    subroutine nc_close_file(nc_id)
        integer :: nc_id,status
        if( ncdfout > 1 )then
            status = check_err(nf90_close(nc_id),'close')
            if(status /= NF90_NOERR )then
                print * ,'cannot close nc file!'
                return
            endif
        endif
    end subroutine

    integer function check_err(status,msg)
    integer,intent(in) :: status
    character(*),intent(in) :: msg
    check_err = status
    if( status .ne. NF90_NOERR )then
        print *, msg//': '//trim(nf90_strerror(status))
        call stopp('.')
    endif
    end function

    subroutine readcorn
    ! *** use GLOBAL ,only:XCR,YCR,UCOR,IC,JC
    use infomod,only:readstr
    implicit none
    integer(4)::i,j,m
    character*80 :: str*200

    if(allocated(xcr)) return
    open(ucor,file = 'corners.inp',action = 'READ')
    str = readstr(ucor)
    allocate(xcr(4,jc_global,ic_global),ycr(4,jc_global,ic_global))
    xcr = missing_value
    ycr = missing_value
    do while(.TRUE.)
        read(ucor,*,end = 100,err = 998) i,j,(xcr(m,j,i),ycr(m,j,i),m = 1,4)
    enddo
100 close(ucor)

    return
998 close(ucor)
    call stopp('CORNERS.INP READING ERROR!')
    end subroutine
    
    ! *** WRITES 1D SCALAR FIELD TO A NETCDF FILE
    integer function nc_write_array1d(nc, id, array1d, offset, count, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer(4),intent(in) :: id
    integer,intent(in) :: count, offset
    real,intent(in) :: array1d(:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:) :: values
    integer :: i,j,status,str1d(2),cnt1d(2)

    allocate(values(count))

    do i = 1,count
        j = i + offset
        values(i) = factor * array1d(j)
    enddo

    str1d = (/ 1, nc.time_idx/)
    cnt1d = (/count, 1/)
    status = nf90_put_var(nc.id, id, values, str1d, cnt1d)
  

    deallocate(values)
    nc_write_array1d = status
    status = check_err(status, 'put_'//msg)
    end function
    
    ! *** WRITES 1D SCALAR FIELD TO A NETCDF FILE
    integer function nc_write_cell_1d(nc, is, ix, array1d, factor, msg)!1, LA_Global-1
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix
    real,intent(in) :: array1d(:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:) :: values
    integer :: i,j,ip,L, count, status,str1d(2),cnt1d(2)

    if( is > 0 )then
        count = NPNT(is)
        allocate(values(count))
        do ip = 1, count
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            values(ip) = factor * array1d(L)
        enddo
    else
        count = LA_Global-1
        allocate(values(count))
        do i = 1,count
            L = i + 1
            values(i) = factor * array1d(L)
        enddo
    endif
    
    str1d = (/ 1, nc.time_idx /)
    cnt1d = (/ count, 1 /)
    status = nf90_put_var(nc.id, nc.idx(ix), values, str1d, cnt1d)

    deallocate(values)
    nc_write_cell_1d = status
    status = check_err(status, 'put_'//msg)
    end function

    integer function nc_write_cell_2d(nc, is, ix, array2d, nk, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix, nk
    real,intent(in) :: array2d(:,:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:,:) :: values
    integer :: i,j,k,L,ip,count,status,str2d(3),cnt2d(3)

    if( is > 0 )then
        count = NPNT(is)
        allocate(values(count, nk))
        do ip = 1, count
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            do k = 1,nk
                values(ip, k) = factor * array2d(L, k)
            enddo
        enddo
    else
        count = LA_Global-1
        allocate(values(count, nk))
        do i = 1, count
            L = i + 1
            do k = 1,nk
                values(i, k) = factor * array2d(L, k)
            enddo
        enddo
    endif

    str2d = (/ 1, 1, nc.time_idx/)
    cnt2d = (/ count, nk, 1/)
    status = nf90_put_var(nc.id, nc.idx(ix), values, str2d, cnt2d)

    deallocate(values)
    nc_write_cell_2d = status
    status = check_err(status, 'put_'//msg)
    end function

    integer function nc_write_wc_2d(nc, is, ix, array2d, nk, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix, nk
    real,intent(in) :: array2d(:,:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:,:) :: values
    integer :: i,j,k,k1,L,ip, count,status,str2d(3),cnt2d(3)

    if( is > 0 )then
        count = NPNT(is)
        allocate(values(count, nk))
        do ip = 1, count
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            k1 = KSZ_Global(L)
            do k = k1,KC
                values(ip, k) = factor * array2d(L, k)
            enddo
        enddo
    else
        count = LA_Global-1
        allocate(values(count, nk))
        values = MISSING_VALUE
        do i = 1, count
            L = i + 1
            k1 = KSZ_Global(L)
            do k = k1,KC
                values(i, k) = factor * array2d(L, k)
            enddo
        enddo
    endif
    
    str2d = (/ 1, 1, nc.time_idx/)
    cnt2d = (/ count, nk, 1/)
    status = nf90_put_var(nc.id, nc.idx(ix), values, str2d, cnt2d)

    deallocate(values)
    nc_write_wc_2d = status
    status = check_err(status, 'put_'//msg)
    end function

    integer function nc_write_wc_3d(nc, is, ix, array3d, nk, ncmp, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix, nk, ncmp
    real,intent(in) :: array3d(:,:,:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:,:,:) :: values
    integer :: i,j,k,k1,L,m,ip, count,status,str3d(4),cnt3d(4)

    if( is > 0 )then
        count = NPNT(is)
        allocate(values(count, nk, ncmp))
        values = MISSING_VALUE
        do ip = 1, count
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            k1 = KSZ_Global(L)
            do k = k1,KC
                do m = 1,ncmp
                    values(ip, k, m) = factor * array3d(L, k, m)
                enddo
            enddo
        enddo
    else
        count = LA_Global-1
        allocate(values(count, nk, ncmp))
        values = MISSING_VALUE
        do i = 1, count
            L = i + 1
            k1 = KSZ_Global(L)
            do k = k1,KC
                do m = 1,ncmp
                    values(i, k, m) = factor * array3d(L, k, m)
                enddo
            enddo
        enddo
    endif
    
    str3d = (/ 1, 1, 1, nc.time_idx/)
    cnt3d = (/ count, nk, ncmp, 1/)
    status = nf90_put_var(nc.id, nc.idx(ix), values, str3d, cnt3d)

    deallocate(values)
    nc_write_wc_3d = status
    status = check_err(status, 'put_'//msg)
    end function

    integer function nc_write_bed_2d(nc, is, ix, array2d, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix
    real,intent(in) :: array2d(:,:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:,:) :: values
    integer :: i,j,L,k,k1,m,ip, count,status,str2d(2),cnt2d(3)

    str2d = (/ 1, nc.time_idx/)
    if( is > 0 )then
        cnt2d = (/ NPNT(is), KB, 1/)
        allocate(values(NPNT(is), KB))
        values = MISSING_VALUE
        do ip = 1, NPNT(is)
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            do k = 1,KBT_Global(L)
              values(ip,k) = factor * array2d(L,k)
            enddo
        enddo
    else
        count = LA_Global-1
        cnt2d = (/count, KB, 1/)
        allocate(values(count, KB))
        values = MISSING_VALUE
        do i = 1, count
            L = i + 1
            do k = 1,KBT_Global(L)
              values(i,k) = factor * array2d(L,k)
            enddo
        enddo
    endif

    status = nf90_put_var(nc.id, nc.idx(ix), values, str2d, cnt2d)

    deallocate(values)
    nc_write_bed_2d = status
    status = check_err(status, 'put_'//msg)
    end function

    integer function nc_write_bed_3d(nc, is, ix, array3d, ncmp, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix, ncmp
    real,intent(in) :: array3d(:,:,:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:,:,:) :: values
    integer :: i,j,L,k,k1,m,ip, count,status,str3d(4),cnt3d(4)

    str3d = (/ 1, 1, 1, nc.time_idx/)
    if( is > 0 )then
        cnt3d = (/ NPNT(is), KB, ncmp, 1/)
        allocate(values(NPNT(is), KB, ncmp))
        values = MISSING_VALUE
        do ip = 1, NPNT(is)
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            do k = 1,KBT_Global(L)
                do m = 1,ncmp
                    values(ip, k, m) = factor * array3d(L, k, m)
                enddo
            enddo
        enddo
    else
        count = LA_Global-1
        cnt3d = (/ count, KB, ncmp, 1/)
        allocate(values(count, KB, ncmp))
        values = MISSING_VALUE
        do i = 1, count
            L = i + 1
            do k = 1,KBT_Global(L)
                do m = 1,ncmp
                    values(i, k, m) = factor * array3d(L, k, m)
                enddo
            enddo
        enddo
    endif    

    status = nf90_put_var(nc.id, nc.idx(ix), values, str3d, cnt3d)

    deallocate(values)
    nc_write_bed_3d = status
    status = check_err(status, 'put_'//msg)
    end function

    integer function nc_write_cell_3d(nc, is, ix, array3d, nk, ncmp, factor, msg)
    type(nc_dataset),intent(in) :: nc
    integer,intent(in) :: is, ix, nk, ncmp
    real,intent(in) :: array3d(:,:,:), factor
    character(*),intent(in) :: msg
    real(4), allocatable, dimension(:,:,:) :: values
    integer :: i,j,L,k,k1,m,ip, count,status,str3d(4),cnt3d(4)
    real(4) :: TMP

    str3d = (/ 1, 1, 1, nc.time_idx/)
    if( is > 0 )then
        cnt3d = (/ NPNT(is), nk, ncmp, 1/)
        allocate(values(NPNT(is), nk, ncmp))
        do ip = 1, NPNT(is)
            i = HFREGRP(is).ICEL(ip)
            j = HFREGRP(is).JCEL(ip)
            L = LIJ_Global(i,j)
            do k = 1, nk
                do m = 1,ncmp
                    values(ip, k, m) = factor * array3d(L, k, m)
                enddo
            enddo
        enddo
    else
        count = LA_Global-1
        cnt3d = (/ count, nk, ncmp, 1/)
        allocate(values(count, nk, ncmp))
        do i = 1, count
            L = i + 1
            do k = 1, nk
                do m = 1,ncmp
                    values(i, k, m) = factor * array3d(L, k, m)
                enddo
            enddo
        enddo
    endif   

    status = nf90_put_var(nc.id, nc.idx(ix), values, str3d, cnt3d)

    deallocate(values)
    nc_write_cell_3d = status
    status = check_err(status, 'put_'//msg)
    
  end function

  subroutine nc_write(nc, is, timesec)
    
    type(nc_dataset),intent(inout) :: nc
    integer,intent(in) :: is
    double precision,intent(in) :: timesec
    integer :: status, i, L, k, k1, k2, m, LN
    integer(2),allocatable :: int1d(:),int2d(:,:)
    real,allocatable :: array1d(:), array2d(:,:), array3d(:,:,:)
    real(8) :: ti(1), xm(1),ym(1),xll(1),yll(1)
    real(4) :: utmp, vtmp, tmp, SHEAR, real1d(1)

    if(process_id /= master_id ) return

    nc.time_idx = nc.time_idx + 1
    ti = TIMEDAY

    if(is == 0) WRITE(*,'(A14,I2,A10,F15.3)') 'NETCDF OUTPUT ',is,' @ HOURS: ',TI(1)*24

    allocate(array1d(LA_Global))
    array1d = MISSING_VALUE

    status = nf90_put_var(nc.id, nc.nc_time, TI, (/nc.time_idx/), (/1/))
    real1d(1) = DELT
    status = nf90_put_var(nc.id, nc.idx(7), real1d, (/nc.time_idx/), (/1/))

    if( IS_NC_OUT(20) > 0 .or. nc.time_idx <= 1 )then
        status = nc_write_cell_1d(nc, is, 1, BELV_Global, 1., 'Zb')
    endif
    
    status = nc_write_cell_1d(nc, is, 2, HP_Global, 1., 'HP')

    do L = 2,LA_Global
        if( HP_Global(L) > HDRY )then
            array1d(L) = BELV_Global(L) + HP_Global(L)
        endif
    enddo
    status = nc_write_cell_1d(nc, is, 3, array1d, 1., 'WSEL')

    allocate(array3d(LA_Global, KC, 3))
    array3d = MISSING_VALUE
    do L = 2,LA_GLobal
        if( HP_GLobal(L) > HDRY )then
            k1 = KSZ_Global(L)
            !WRITE(EE_UNIT) ((REAL(SAL_Global(L,K),4), K = KSZ_Global(L),KC), L = 2,LA_Global)
            if( ROTA == 1 )then
                LN = LNC_GLobal(L)
                do k = k1,KC
                    utmp = 0.5*(RSSBCE_GLobal(L)*U_GLobal(L+1,K) + RSSBCW_GLobal(L)*U_GLobal(L,K))
                    vtmp = 0.5*(RSSBCN_GLobal(L)*V_GLobal(LN ,K) + RSSBCS_GLobal(L)*V_GLobal(L,K))
                    array3d(L, k, 1) = CUE_Global(L)*utmp + CVE_Global(L)*vtmp
                    array3d(L, k, 2) = CUN_Global(L)*utmp + CVN_Global(L)*vtmp
                    array3d(L, k, 3) = W_Global(L,K)
                enddo
            else
                do k = k1,KC
                    array3d(L, k, 1) = U_Global(L,K)
                    array3d(L, k, 2) = V_Global(L,K)
                    array3d(L, k, 3) = W_Global(L,K)
                enddo
            endif
        endif
    enddo
    status = nc_write_wc_2d(nc, is, 4, array3d(:,:,1), KC, 1., 'U')
    status = nc_write_wc_2d(nc, is, 5, array3d(:,:,2), KC, 1., 'V')
    status = nc_write_wc_2d(nc, is, 6, array3d(:,:,3), KC, 1., 'W')
    deallocate(array3d)

    if( ISSPH(8) >= 1 )then
        if( IS_NC_OUT(21) > 0 .or. IS_NC_OUT(22) > 0 )then
            ! *** WRITE THE TOP LAYER INDEX
            status = nc_write_cell_1d(nc, is, 31,float(KBT_Global) , 1., 'KBT')
            !allocate(int1d(LA_Global-1))
            !int1d = 0
            !do L = 2,LA_Global
            !    int1d(L-1) = KBT_Global(L)
            !enddo
            !status = nf90_put_var(nc.id, nc.idx(31), int1d, (/ 1, nc.time_idx/), (/ nc.cellcnt, 1/))
            !deallocate(int1d)
        endif

        if( IS_NC_OUT(2) > 0 )then
            ! *** TOTAL BED SHEAR STRESS
            do L = 2,LA_Global
                array1d(L) = SHEAR_Global(L)*RHOW(L,KSZ_Global(L)) ! Convert from m2/s2 to Pa
            enddo
            status = nc_write_cell_1d(nc, is, 8, array1d, 1., 'SHEAR')

            if( ISBEDSTR == 1 .and. .not. LSEDZLJ )then
                ! *** Stress from non-cohesive components
                status = nc_write_cell_1d(nc, is, 9, SHEAR_Global2, 1., 'SHEAR2')
            endif

            if( ISWAVE >= 1 )then
                ! *** Bed Shear due to Waves Only
                do L = 2,LA_Global
                  array1d(L) = QQWV3_Global(L)*RHOW(L,KSZ_Global(L)) ! Convert from m2/s2 to Pa
                enddo
                status = nc_write_cell_1d(nc, is, 12, array1d, 1., 'SHEAR2')
                ! *** Shear due to Current Only
                ! @todo - calculate using global values.  Need to get a global LEC_Global,LNC   delme
                do L = 2,LA_Global
                    SHEAR = ( RSSBCE_Global(L)*TBX_Global(LEC_Global(L)) + RSSBCW_Global(L)*TBX_Global(L) )**2  + &
                        ( RSSBCN_Global(L)*TBY_Global(LNC_Global(L)) + RSSBCS_Global(L)*TBY_Global(L) )**2
                    SHEAR = 0.5*SQRT(SHEAR)
                    array1d(L) = SHEAR*RHOW(L,KSZ_Global(L))            ! *** Bed Shear due to Current Only in Pa
                enddo
                status = nc_write_cell_1d(nc, is, 11, array1d, 1., 'SHEAR1')
                if( ISWAVE >= 3 )then
                    do L = 2,LA_Global
                        array1d(L) = WV_Global(L).HEIGHT
                    enddo
                    status = nc_write_cell_1d(nc, is, 15, array1d, 1., 'Hs')
                    do L = 2,LA_Global
                        array1d(L) = WV_Global(L).DIR
                    enddo
                    status = nc_write_cell_1d(nc, is, 16, array1d, 1., 'Dp')
                    do L = 2,LA_Global
                        array1d(L) = WV_Global(L).FREQ
                    enddo
                    status = nc_write_cell_1d(nc, is, 17, array1d, 1., 'Tp')
                    if( ISWAVE == 4 )then
                        do L = 2,LA_Global
                            array1d(L) = WV_Global(L).FREQ
                        enddo
                        status = nc_write_cell_1d(nc, is, 18, array1d, 1., 'Diss')
                        status = nc_write_cell_1d(nc, is, 19, WVHUU_Global(:,KC), 1., 'Sxx')
                        status = nc_write_cell_1d(nc, is, 20, WVHUV_Global(:,KC), 1., 'Sxy')
                        status = nc_write_cell_1d(nc, is, 21, WVHVV_Global(:,KC), 1., 'Syy')
                    endif
                endif
            endif
        endif

        if( IS_NC_OUT(16) > 0 ) status = nc_write_wc_2d(nc, is, 22, SAL_Global, KC, 1., 'Salt')     ! *** Salinity

        if( IS_NC_OUT(13) > 0 ) status = nc_write_wc_2d(nc, is, 23, TEM_Global, KC, 1., 'TEM')      ! *** Temperature
        if( IS_NC_OUT(14) > 0 ) status = nc_write_cell_1d(nc, is, 24, TEMB_Global, 1., 'TEMB')      ! *** Bed temperature
        if( IS_NC_OUT(7) > 0 )  status = nc_write_cell_1d(nc, is, 25, EVAPT_Global*1000, 1., 'EVAPT')    ! *** Evaporation
        if( IS_NC_OUT(6) > 0 )  status = nc_write_cell_1d(nc, is, 26, RAINT_Global*86400000, 1., 'RAINT')    ! *** Rainfall
        if( IS_NC_OUT(17) > 0 ) status = nc_write_wc_3d(nc, is, 27, DYE_Global, KC, NDYE, 1.0, 'DYE') ! *** Dye
        if( IS_NC_OUT(49) > 0 ) status = nc_write_wc_2d(nc, is, 28, SFL_Global, KC, 1., 'SFL')      ! *** Shellfish larvae

        if( IS_NC_OUT(24) > 0 ) status = nc_write_wc_3d(nc, is, 29, TOX_Global, KC, NTOX, 1.0, 'TOX') ! *** Toxics
        if( IS_NC_OUT(25) > 0 ) status = nc_write_bed_3d(nc, is, 30, TOXB_Global, NTOX, 1., 'TOXB') ! *** Bed toxics

        if( IS_NC_OUT(21) > 0 .or. IS_NC_OUT(22) > 0 )then
            status = nc_write_bed_2d(nc, is, 32, HBED_Global, 1.0, 'HBED')
            status = nc_write_bed_2d(nc, is, 33, BDENBED_Global, 1.0, 'BDENBED')
            status = nc_write_bed_2d(nc, is, 34, PORBED_Global, 1.0, 'PORBED')
        endif
        if( IS_NC_OUT(18) > 0 ) status = nc_write_wc_3d(nc, is, 37, SED_Global, KC, NSED, 1.0, 'SED') ! *** Cohesive sediment concentration
        if( IS_NC_OUT(19) > 0 ) status = nc_write_wc_3d(nc, is, 38, SND_GLobal, KC, NSND, 1.0, 'SND') ! *** Noncohesive sediment concentration
        if( IS_NC_OUT(21) > 0 ) status = nc_write_bed_3d(nc, is, 35, SEDB_Global, NSED, 1., 'SEDB')   ! *** Bed cohesive sediment
        if( IS_NC_OUT(22) > 0 ) status = nc_write_bed_3d(nc, is, 36, SNDB_Global, NSND, 1., 'SNDB')   ! *** Bed noncohesive sediment
        
        if( ICALC_BL > 0 .and. NSND > 0 )then
          if( IS_NC_OUT(23) > 0 )then
            status = nc_write_cell_2d(nc, is, 48, real(QSBDLDX_Global), NSND, 1., 'QSBDLDX') ! *** Bedload_x
            status = nc_write_cell_2d(nc, is, 49, real(QSBDLDY_Global), NSND, 1., 'QSBDLDY') ! *** Bedload_y
          endif
        endif
    
        if( IS_NC_OUT(8) > 0 )then
            status = nc_write_cell_1d(nc, is, 51, EVAPSW_Global, 1., 'EVAPSW')
            status = nc_write_cell_1d(nc, is, 52, EVAPGW_Global, 1., 'EVAPGW')
            status = nc_write_cell_1d(nc, is, 53, QGW_Global, 1., 'QGW')
            status = nc_write_cell_1d(nc, is, 54, AGWELV_Global, 1., 'AGWELV')
        endif

        if( IS_NC_OUT(15) > 0 )then
            ! *** Ice thickness
            status = nc_write_cell_1d(nc, is, 55, ICETHICK_Global, 1., 'ICETHICK')
            ! *** Ice temperature
            status = nc_write_cell_1d(nc, is, 56, ICETEMP_Global, 1., 'ICETEMP')
        endif
    endif
    deallocate(array1d)

    if( LSEDZLJ )then
        ! *** This is used to write out the active layers instead of using the LAYERACTIVE array.
        allocate(int2d(LA_Global-1,KB))
        int2d = 0
        ! *** LAYERACTIVE(KB,LCM) - This is = 1 when a bed layer (KB index) exists with mass
        do L = 2,LA_Global
            do K = 1,KB
                if( TSED_GLOBAL(K,L) > 1E-8 )then
                    int2d(L-1, k) = 1
                endif
            enddo
        enddo
        status = nf90_put_var(nc.id, nc.idx(39), int2d, (/1, 1, nc.time_idx/), (/KB, nc.cellcnt, 1/))
        deallocate(int2d)
        
        allocate(array2d(LA_Global, KB))
        array2d = 0.
        ! *** TAU(LCM)      - Shear Stress in dynes/cm^2
        status = nc_write_cell_1d(nc, is, 40, real(TAU_Global), 1., 'TAU')
        ! *** TSED(KB,LCM)  - This is the mass in g/cm^2 in each layer
        !status = nf90_put_var(nc.id, nc.idx(41), TSED_Global(:,2:LA_Global), (/1, 1, nc.time_idx/), (/ KB, LA_Global-1, 1/))
        do L = 2,LA_Global
            do K = 1,KB
                array2d(L, k) = TSED_Global(k, L)
            enddo
        enddo
        status = nc_write_cell_2d(nc, is, 41, array2d, KB, 1., 'TSED')
        ! *** BULKDENS(KB,LCM) - Dry Bulk density of each layer (g/cm^3)
        !status = nf90_put_var(nc.id, nc.idx(42), BULKDENS_Global(:,2:LA_Global), (/1, 1, nc.time_idx/), (/ KB, LA_Global-1, 1/))
        do L = 2,LA_Global
            do K = 1,KB
                array2d(L, k) = BULKDENS_Global(k, L)
            enddo
        enddo
        status = nc_write_cell_2d(nc, is, 42, array2d, KB, 1., 'BULKDENS')
        deallocate(array2d)
        ! *** PERSED(NSEDS,KB,LCM) - This is the mass percentage of each size class in a layer
        !status = nf90_put_var(nc.id, nc.idx(43), PERSED_Global(:,:,2:LA_Global), (/1, 1, 1, nc.time_idx/), (/NSEDS, KB, LA_Global-1, 1/))
        allocate(array3d(LA_Global, KB, NSEDS))
        do L = 2,LA_Global
            do K = 1,KB
                do i = 1,NSEDS
                    array3d(L, k, i) = PERSED_Global(i, k, L)
                enddo
            enddo
        enddo
        status = nc_write_cell_3d(nc, is, 43, array3d, KB, NSEDS, 1., 'PERSED')
        deallocate(array3d)
        ! *** D50AVG(LCM)   - Average particle size of bed surface (microns)
        status = nc_write_cell_1d(nc, is, 44, real(D50AVG_Global), 1., 'D50AVG')
        ! *** ADD IN DEP_SED_FLX_GLOBAL AND ERO_SED_FLX_GLOBAL HERE
        ! *** ADD IN DEP_SED_FLX_GLOBAL AND ERO_SED_FLX_GLOBAL HERE

        if( ICALC_BL > 0 )then
            ! *** Bedload mass inb g/m^2 for each size class
            status = nc_write_cell_2d(nc, is, 47, real(CBL_Global), NSEDS, 10000., 'CBL')
            ! *** Bedload flux in X direction (g/s)
            status = nc_write_cell_2d(nc, is, 48, real(QSBDLDX_Global), NSEDS, 1., 'QSBDLDX')
            ! *** Bedload flux in Y direction (g/s)
            status = nc_write_cell_2d(nc, is, 49, real(QSBDLDY_Global), NSEDS, 1., 'QSBDLDY')
        endif

        if( IS_NC_OUT(25) > 0 )then
            ! *** Bed toxics
            allocate(array3d(LA_Global, KB, NTOX))
            array3d = MISSING_VALUE
            do i = 1,NTOX
                do L = 2,LA_Global
                    do k = 1,KB
                        if( k < 3 .and. TSED_Global(K,L) > 0. .and. HBED_Global(L,KBT_Global(L)) > 0. )then  ! *** avoid division by HBED = 0 --- DKT
                            tmp = 0.01*TSED_Global(K,L)/BULKDENS_Global(K,L)             ! *** HBED(L,K)
                            tmp = tmp / HBED_Global(L,KBT_Global(L))
                            tmp = tmp * TOXB_Global(L,KBT_Global(L),i)
                            array3d(L, k, i) = tmp
                        elseif( TSED_Global(K,L) > 0. )then
                            array3d(L, k, i) = TOXB_Global(L,K,i)
                        else
                            array3d(L, k, i) = 0.0_4
                        endif
                    enddo
                enddo
            enddo
            status = nc_write_cell_3d(nc, is, 30, array3d, KB, NTOX, 1.0, 'TOXB')
            deallocate(array3d)
        endif
        if( IS_NC_OUT(25) > 0 .and. ICALC_BL > 0 )then
            ! *** Bedload toxic concentration (mg/m^2)
            status = nc_write_cell_2d(nc, is, 50, real(CBLTOX_Global), NTOX, 1., 'CBLTOX')
        endif
    else
        if( ISBEXP >= 1 .and. KB > 1 .and. ISSPH(8) < 1 )then
            if( IS_NC_OUT(21) > 0 .or. IS_NC_OUT(22) > 0 )then
                !allocate(int1d(LA_Global-1))
                !int1d = 0
                !do L = 2,LA_Global
                !    int1d(L-1) = KBT_Global(L)
                !enddo
                !status = nf90_put_var(nc.id, nc.idx(31), int1d, (/ 1, nc.time_idx/), (/ nc.cellcnt, 1/))
                !deallocate(int1d)                
                !status = nf90_put_var(nc.id, nc.idx(31), integer(KBT_Global(2:LA_Global),2), (/1, nc.time_idx/), (/nc.cellcnt, 1/))
                status = nc_write_cell_1d(nc, is, 31, float(KBT_Global), 1., 'KBT')
                status = nc_write_cell_2d(nc, is, 32, HBED_Global, KB, 1., 'HBED')
                status = nc_write_cell_2d(nc, is, 33, BDENBED_Global, KB, 1., 'BDENBED')
                status = nc_write_cell_2d(nc, is, 34, PORBED_Global, KB, 1., 'PORBED')
            endif
            if( IS_NC_OUT(21) > 0 ) status = nc_write_cell_3d(nc, is, 35, SEDB_Global, KB, NSED, 1., 'SEDB') ! *** Bed cohesive
            if( IS_NC_OUT(22) > 0 ) status = nc_write_cell_3d(nc, is, 36, SNDB_Global, KB, NSND, 1., 'SNDB') ! *** Bed noncohesive
            if( IS_NC_OUT(25) > 0 )then
                ! *** Bed toxics
                allocate(array2d(LA_Global, NTOX))
                array2d = MISSING_VALUE
                do i = 1,NTOX
                    do L = 2,LA_Global
                        array2d(L, i) = TOXB_Global(L,KBT_Global(L),i)
                    enddo
                enddo
                status = nc_write_wc_2d(nc, is, 30, array2d, NTOX, 1.0, 'TOXB')
                deallocate(array2d)
            endif
        endif
    endif

    do i = 1,18
        if( IS_NC_OUT(i+25) > 0 )then
            status = nc_write_wc_2d(nc, is, 56+i, WQV_Global(:,:,i), KC, 1., 'WQ')  ! *** Water quality
        endif
    enddo        
    if( IS_NC_OUT(45) > 0 )then
        k1 = 20
        k2 = 19 + NALGAE
        status = nc_write_wc_3d(nc, is, 76, WQV_Global(:,:,k1:k2), KC, NALGAE, 1., 'ALGAE') ! *** Algal groups
    endif
    if( IS_NC_OUT(46) > 0 )then     
        k1 = 19 + NALGAE + 1
        k2 = k1 + NZOOPL - 1
        status = nc_write_wc_3d(nc, is, 77, WQV_Global(:,:,k1:k2), KC, NZOOPL, 1., 'ZOOP')  ! *** Zooplankton groups
    endif
    if( IS_NC_OUT(50) > 0 )then                ! *** Sediment Diagenesis
        status = nc_write_cell_1d(nc, is,  78, SMPON_Global(:,1), 1., 'SMPON1')
        status = nc_write_cell_1d(nc, is,  79, SMPON_Global(:,2), 1., 'SMPON2')
        status = nc_write_cell_1d(nc, is,  80, SMPON_Global(:,3), 1., 'SMPON3')
        status = nc_write_cell_1d(nc, is,  81, SMPOP_Global(:,1), 1., 'SMPOP1')
        status = nc_write_cell_1d(nc, is,  82, SMPOP_Global(:,2), 1., 'SMPOP2')
        status = nc_write_cell_1d(nc, is,  83, SMPOP_Global(:,3), 1., 'SMPOP3')
        status = nc_write_cell_1d(nc, is,  84, SMPOC_Global(:,1), 1., 'SMPOC1')
        status = nc_write_cell_1d(nc, is,  85, SMPOC_Global(:,2), 1., 'SMPOC2')
        status = nc_write_cell_1d(nc, is,  86, SMPOC_Global(:,3), 1., 'SMPOC3')
        status = nc_write_cell_1d(nc, is,  87, SMDFN_Global(:,1), 1., 'SMDFN1')
        status = nc_write_cell_1d(nc, is,  88, SMDFN_Global(:,2), 1., 'SMDFN2')
        status = nc_write_cell_1d(nc, is,  89, SMDFN_Global(:,3), 1., 'SMDFN3')
        status = nc_write_cell_1d(nc, is,  90, SMDFP_Global(:,1), 1., 'SMDFP1')
        status = nc_write_cell_1d(nc, is,  91, SMDFP_Global(:,2), 1., 'SMDFP2')
        status = nc_write_cell_1d(nc, is,  92, SMDFP_Global(:,3), 1., 'SMDFP3')
        status = nc_write_cell_1d(nc, is,  93, SMDFC_Global(:,1), 1., 'SMDFC1')
        status = nc_write_cell_1d(nc, is,  94, SMDFC_Global(:,2), 1., 'SMDFC2')
        status = nc_write_cell_1d(nc, is,  95, SMDFC_Global(:,3), 1., 'SMDFC3')
        status = nc_write_cell_1d(nc, is,  96, SM1NH4_Global(:), 1., 'SM1NH4')
        status = nc_write_cell_1d(nc, is,  97, SM2NH4_Global(:), 1., 'SM2NH4')
        status = nc_write_cell_1d(nc, is,  98, SM1NO3_Global(:), 1., 'SM1NO3')
        status = nc_write_cell_1d(nc, is,  99, SM2NO3_Global(:), 1., 'SM2NO3')
        status = nc_write_cell_1d(nc, is, 100, SM1PO4_Global(:), 1., 'SM1PO4')
        status = nc_write_cell_1d(nc, is, 101, SM2PO4_Global(:), 1., 'SM2PO4')
        status = nc_write_cell_1d(nc, is, 102, SM1H2S_Global(:), 1., 'SM1H2S')
        status = nc_write_cell_1d(nc, is, 103, SM2H2S_Global(:), 1., 'SM2H2S')
        status = nc_write_cell_1d(nc, is, 104, SM1SI_Global(:), 1., 'SM1SI')
        status = nc_write_cell_1d(nc, is, 105, SM2SI_Global(:), 1., 'SM2SI')
        status = nc_write_cell_1d(nc, is, 106, SMPSI_Global(:), 1., 'SMPSI')
        status = nc_write_cell_1d(nc, is, 107, SMBST_Global(:), 1., 'SMBST')
        status = nc_write_cell_1d(nc, is, 108, SMT_Global(:), 1., 'SMT')
        status = nc_write_cell_1d(nc, is, 109, SMCSOD_Global(:), 1., 'SMCSOD')
        status = nc_write_cell_1d(nc, is, 110, SMNSOD_Global(:), 1., 'SMNSOD')
        status = nc_write_cell_1d(nc, is, 111, WQBFNH4_Global(:), 1., 'WQBFNH4')
        status = nc_write_cell_1d(nc, is, 112, WQBFNO3_Global(:), 1., 'WQBFNO3')
        status = nc_write_cell_1d(nc, is, 113, WQBFO2_Global(:), 1., 'WQBFO2')
        status = nc_write_cell_1d(nc, is, 114, WQBFCOD_Global(:), 1., 'WQBFCOD')
        status = nc_write_cell_1d(nc, is, 115, WQBFPO4D_Global(:), 1., 'WQBFPO4D')
        status = nc_write_cell_1d(nc, is, 116, WQBFSAD_Global(:), 1., 'WQBFSAD')        
    endif
    if( IS_NC_OUT(47) > 0 )then    ! *** RPEM
        status = nc_write_cell_1d(nc, is, 117, WQRPS_Global(:), 1., 'WQRPS')
        status = nc_write_cell_1d(nc, is, 118, WQRPR_Global(:), 1., 'WQRPR')
        status = nc_write_cell_1d(nc, is, 119, WQRPE_Global(:), 1., 'WQRPE')
        status = nc_write_cell_1d(nc, is, 120, WQRPD_Global(:), 1., 'WQRPD')
    endif
    if( IS_NC_OUT(48) > 0 )then    ! *** Shellfish Farm
        allocate(array3d(LA_Global, KC, NSF))
        array3d = MISSING_VALUE
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).CELL_C(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 132, array3d, KC, NSF, 1.0, 'CELL_C')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).INDI_C(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 133, array3d, KC, NSF, 1.0, 'INDI_C')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).SHLEN(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 134, array3d, KC, NSF, 1.0, 'SHLEN')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).INDI(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 135, array3d, KC, NSF, 1.0, 'INDI')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).INDI_D(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 136, array3d, KC, NSF, 1.0, 'INDI_D')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).INDI_H(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 137, array3d, KC, NSF, 1.0, 'INDI_H')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).HARV_C(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 138, array3d, KC, NSF, 1.0, 'HARV_C')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).NP(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 139, array3d, KC, NSF, 1.0, 'NP')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).PR(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 140, array3d, KC, NSF, 1.0, 'PR')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).CRB(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 141, array3d, KC, NSF, 1.0, 'CRB')
        do i = 1,NSF
            do L = 2,LA_Global
                do k = 1,KC
                    array3d(L, k, i) = SF(i).SPAWN(L,K)
                enddo
            enddo
        enddo
        status = nc_write_wc_3d(nc, is, 142, array3d, KC, NSF, 1.0, 'SPAWN')
        deallocate(array3d)

        allocate(array2d(LA_Global, NSF))
        array2d = MISSING_VALUE

        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).WQAQ(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 143, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).WQEA(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 144, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).FR(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 145, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).BMG(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 146, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_RESPI(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 147, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_GRAZI(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 148, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_DEATH(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 149, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_FECAL(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 150, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_URINE(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 151, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_RPOC(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 152, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_LPOC(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 153, array2d, NSF, 1.0, 'WQAQ')
        do i = 1,NSF
            do L = 2,LA_Global
                array2d(L, i) = SF(i).B_DOC(L)
            enddo
        enddo
        status = nc_write_wc_2d(nc, is, 154, array2d, NSF, 1.0, 'WQAQ')
        deallocate(array2d)
    endif
    
    if( IS_NC_OUT(9) > 0 )then     ! *** Hydraulic Structures
        status = nf90_put_var(nc.id, nc.idx(121), real(QCTLT(:, :, 1),4), (/ 1, 1, nc.time_idx/), (/ KC, NQCTL, 1/))
        status = nf90_put_var(nc.id, nc.idx(122), real(QCTLT(:, :, 2),4), (/ 1, 1, nc.time_idx/), (/ KC, NQCTL, 1/))

        if(NQCTLSER > 0 .or. NQCTRULES > 0 )then
            allocate(int1d(NQCTL))
            int1d = 0
            do i = 1,NQCTL
                int1d(i) = HSCTL_GL(i).CUR.STATE
            enddo
            status = nf90_put_var(nc.id, nc.idx(123), int1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            !do i = 1,NQCTL
            !    int1d(i) = HSCTL_GL(i).CUR.UNITS
            !enddo
            !status = nf90_put_var(nc.id, nc.idx(124), int1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            do i = 1,NQCTL
                int1d(i) = HSCTL_GL(i).CUR.ID
            enddo
            status = nf90_put_var(nc.id, nc.idx(125), int1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            deallocate(int1d)
            allocate(array1d(NQCTL))
            do i = 1,NQCTL
                array1d(i) = HSCTL_GL(i).CUR.FLOW
            enddo
            status = nf90_put_var(nc.id, nc.idx(126), array1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            do i = 1,NQCTL
                array1d(i) = HSCTL_GL(i).CUR.HEIGHT
            enddo
            status = nf90_put_var(nc.id, nc.idx(127), array1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            do i = 1,NQCTL
                array1d(i) = HSCTL_GL(i).CUR.WIDTH
            enddo
            status = nf90_put_var(nc.id, nc.idx(128), array1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            do i = 1,NQCTL
                array1d(i) = HSCTL_GL(i).CUR.SILL
            enddo
            status = nf90_put_var(nc.id, nc.idx(129), array1d, (/ 1, nc.time_idx/), (/ NQCTL, 1/))
            deallocate(array1d)
        endif
    endif
    if( IS_NC_OUT(10) > 0 )then    ! *** Withdrawal/return
        allocate(int1d(NQWR))
        int1d = 0
        do i = 1,NQWR
            int1d(i) = WITH_RET_CTL(i).CUR.STATE
        enddo
        status = nf90_put_var(nc.id, nc.idx(130), int1d, (/ 1, nc.time_idx/), (/ NQWR, 1/))
        deallocate(int1d)
        allocate(array1d(NQWR))
        do i = 1,NQWR
            array1d(i) = WITH_RET_CTL(i).CUR.FLOW
        enddo
        status = nf90_put_var(nc.id, nc.idx(131), array1d, (/ 1, nc.time_idx/), (/ NQWR, 1/))
        deallocate(array1d)
    endif

    if( IS_NC_OUT(4) > 0 )then     ! *** Winds
        status = nc_write_cell_1d(nc, is, 13, WNDVELE_Global, 1., 'Wx')
        status = nc_write_cell_1d(nc, is, 14, WNDVELN_Global, 1., 'Wy')
    endif
    if( IS_NC_OUT(5) > 0 )then     ! *** Barometric pressure
        status = nc_write_cell_1d(nc, is, 155, PATMT_Global, 1., 'Pa')
    endif
    
    end subroutine

    subroutine nc_output()
    real(8) :: timehour
    character(8) :: datestamp
    character(80) :: filename
    integer :: status
    logical :: file_exists

    if(process_id /= master_id ) return
    
    if( .not.allocated(nc)) allocate(nc(0:NSUBSET)) 

    if( issglfil == 1 )then
        ! *** only create the file at the first call TOGREGOR this subroutine
        filename = './/#output//efdc_out.nc'
        INQUIRE(FILE = filename, EXIST = file_exists)
        if( (nc(0).isopen == 1) .and. (file_exists) )then
            status = nc_open_file(nc(0).id, filename)
        else
            nc(0) = nc_create_file(filename)
            nc(0).time_idx = 0
            nc(0).isopen = 1
        endif
    elseif( issglfil == 0 )then
        timehour = timeday*24
        call togregor(datestamp,timeday)
        filename = './/#output//efdc_'//datestamp//'.nc'
        INQUIRE(FILE = filename, EXIST = file_exists)
        if((nc_datestamp == datestamp) .and. (file_exists) )then
            status = nc_open_file(nc(0).id, filename)
        else
            if( nc(0).isopen == 1 )then
                call nc_write(nc(0), 0, timesec)
                call nc_close_file(nc(0).id)  
                if( timeday >= tendncdf) return
            endif
            nc(0) = nc_create_file(filename)
            nc(0).time_idx = 0
            nc(0).isopen = 1
        endif
        nc_datestamp = datestamp
    endif 
    
    call nc_write(nc(0), 0, timesec)
    call nc_close_file(nc(0).id)  
    
    end subroutine
    
    subroutine nc_output_hf(ncal,is)
    integer,intent(in) :: ncal,is
    real(8) :: timehour
    character(3) :: nstr
    character(80) :: filename
    integer :: status
    logical :: file_exists

    if(process_id /= master_id ) return
    
    if( .not.allocated(nc)) allocate(nc(0:NSUBSET)) 
    
    write(nstr,'(I3.3)') is
    filename  = OUTDIR//'efdc_hf_'//trim(nstr)//'.nc'
    
    INQUIRE(FILE = filename, EXIST = file_exists)
    if( (nc(is).isopen == 1) .and. (file_exists) )then
        status = nc_open_file(nc(is).id, filename)  !ncal = 1
    else
        nc(is) = nc_create_hf_file(is, filename)    !ncal = 0
        nc(is).time_idx = 0
        nc(is).isopen = 1
    endif
    
    call nc_write(nc(is), is, timesec)
    call nc_close_file(nc(is).id)  
    end subroutine

  end module
