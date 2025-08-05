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
!< @details Define type that characterizes a ship and its propeller
!! When an instance of this type is created it initializes all starting values
!! to zero
!< @author  Zander Mausolff & Paul Craig
!---------------------------------------------------------------------------!
Module Mod_Ship

  use GLOBAL, only : rkd, rk4, ik4

  implicit none
    
  private
    
  public :: ship_type
    
  ! *** 
  type :: ship_type    
    integer                :: mmsi               = 0        !< Identifier of the ship
    character(len = 100)     :: name               = ' '      !< Name of the ship
    integer (kind = ik4)   :: power_source       = 0        !< Type of power source to use: 0 - Horsepower, 1 - Propeller RPM
    real    (kind = rkd)   :: length             = 0.0      !< Length of the ship  [meters]
    real    (kind = rkd)   :: beam               = 0.0      !< Max width of the ship  [meters]
    real    (kind = rkd)   :: max_draft          = 0.0      !< Max draught of the ship [meters]
    real    (kind = rkd)   :: min_draft          = 0.0      !< Min draught of the ship [meters]
    real    (kind = rkd)   :: max_power          = 0.0      !< Max power of the ship [horsepower]
    real    (kind = rkd)   :: max_rps            = 0.0      !< Max revolutions per second
    real    (kind = rkd)   :: ais_to_stern       = 0.0      !< Distance from AIS antenna to stern of ship [meters]
    real    (kind = rkd)   :: dist_between_props = 0.0      !< Distance between propellers [meters]
    real    (kind = rkd)   :: dist_from_stern    = 0.0      !< Horizontal distance from the stern to the propeller(s) [meters]
    real    (kind = rkd)   :: prop_offset        = 0.0      !< Distance of propeller shaft from draft [meters]
                                                            
    integer (kind = ik4)   :: num_props          = 0        !< Number of propellers   
    integer (kind = ik4)   :: num_blades         = 0        !< z   : Number of blades
    integer (kind = ik4)   :: ducted             = 0        !< Flag for propeller nozzles, 0-None, 1-Ducted
    real    (kind = rkd)   :: thrust_coeff       = 0.0      !< Ct  : Propeller Thrust Coefficient
    real    (kind = rkd)   :: prop_diam          = 0.0      !< D   : Propeller diameter [meters]
    real    (kind = rkd)   :: prop_hub_diam      = 0.0      !< Dh  : Propeller hub diameter [meters]
    real    (kind = rkd)   :: blade_area_ratio   = 0.0      !< BAR : Blade area ratio
    real    (kind = rkd)   :: pitch_ratio        = 0.0      !< PD  : Blade pitch ratio

    real    (kind = rkd)   :: freq_out           = 0.0      !< Subgrid mesh output frequncy [minutes]
    
    integer (kind = ik4), allocatable :: mesh_count(:)       !< Count of number of subgrid cells are in a model cell

    integer (kind = ik4)   :: num_fixed_cells     = 0        !< Fix number and cells for test cases  
    integer (kind = ik4), allocatable :: fixed_cells(:)      !< Cell ID for fixed cells
    real    (kind = RK4), allocatable :: fixed_frac(:)       !< Fraction of energy per cell

  contains 
    
    procedure, pass(self) :: write_out !< writes out the basic information to a file
        
  end type ship_type
     
  contains
  !! *** Create constructor interface
  !interface propeller
  !    module procedure :: constructor_propeller
  !end interface propeller
  !
  !contains
  !
  !! *** Create constructor to initialize propeller objects
  !type(ship) function constructor_propeller() result(self)
  !
  !    self.mmsi = 0
  !    self.name = ' '
  !    self.thrust_coeff     = 0.0
  !    self.max_rps          = 0.0
  !    self.prop_diam        = 0.0
  !    self.prop_hub_diam    = 0.0
  !    self.blade_area_ratio = 0.0
  !    self.pitch_ratio      = 0.0
  !    self.num_blades       = 0
  !    
  !end function constructor_propeller
    
  !---------------------------------------------------------------------------!
  !< @details writes the prop info read in to a file
  !< @param[in] unit_num
  !< @author  Zander Mausolff
  !---------------------------------------------------------------------------!
  subroutine write_out(self, unit_num)
    
    implicit none
        
    class(ship_type), intent(inout) :: self !< unit number of the file to write out to
    integer, intent(in) :: unit_num 
        
    write(unit_num, '(a)')           '-----------------------'
    !write(unit_num, '(a,i10)')       'Ship MMSI:             ', self.mmsi
    !write(unit_num, '(a,a) ')        'Ship Name:             ', trim(self.name)     
    write(unit_num, '(a,f10.3) ')    'Thrust Coeff:          ', self.thrust_coeff
    if( self.power_source == 0 )then
      write(unit_num, '(a,f10.3) ')  'Max Power:             ', self.max_power         
    elseif( self.power_source == 1 )then
      write(unit_num, '(a,f10.3) ')  'Max RPS Prop:          ', self.max_rps         
    endif
    write(unit_num, '(a,f10.3) ')    'Prop Diameter [m]:     ', self.prop_diam        
    write(unit_num, '(a,f10.3) ')    'Prop Hub Diameter [m]: ', self.prop_hub_diam    
    write(unit_num, '(a,f10.3) ')    'Blade Area Ratio:      ', self.blade_area_ratio 
    write(unit_num, '(a,f10.3) ')    'Pitch Ratio:           ', self.pitch_ratio      
    write(unit_num, '(a,i4) ')       'Number Blades:         ', self.num_blades     

  end subroutine write_out
    
End module Mod_Ship
