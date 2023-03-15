! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!-------------------------------------------------------------------------------------------!

! Module: Propeller Wash: Variables_Propwash
!
!> @authors Paul M. Craig with Zander Mausolff and Luis Bastidas
!-------------------------------------------------------------------------------------------!
! CHANGE RECORD
! DATE MODIFIED   BY                DESCRIPTION
!-------------------------------------------------------------------------------------------!
!    2020-01      Zander Mausolff   Initially added basic propwash structure
!    2020-12      Paul M. Craig     Finished propeller wash with linkage to the SEDZLJ
!                                   and toxics module
!    2021-12      Paul M. Craig     Added propeller jet efflux momentum to EFDC+ flow field
!-------------------------------------------------------------------------------------------!
  Module Variables_Propwash

  Use GLOBAL    
    
  implicit none
    
  integer, parameter :: prop_log_unit = 668
  
  integer :: ISPROPWASH                                !< Primary flag for propwash options
  integer :: total_ships                               !< Total number of ships to be simulated
                                                       
  Logical :: propwash_on = .FALSE.                     !< Decides if we are executing the propwash module

  integer, allocatable, dimension(:,:) :: adjacent_l
  integer(kind = rk4) :: num_radial_elems              !< Number of propwash mesh points perpindicular to axis of propeller
  integer(kind = rk4) :: num_axial_elems               !< Number of propwash mesh points along axis of propeller
                                                       
  real(kind = rkd) :: TPROPW                           !< Run time logger
  real(kind = rkd) :: mesh_width    = 30.              !< Number of propeller diameters wide
  real(kind = rkd) :: mesh_Length   = 60.              !< Number of propeller diameters long
  real(kind = rkd) :: efflux_zone_mult    = 0.35       !< multiplier for the size of the efflux zone
  real(kind = rkd) :: flow_est_zone1_mult = 3.5        !< multiplier for the size of the flow establishment zone (single prop or no-influnce-zone)
  real(kind = rkd) :: flow_est_zone2_mult = 14.        !< multiplier for the size of the flow establishment zone (dual prop influence zone)
  real(kind = rkd) :: efflux_mag_mult     = 0.75       !< multiplier to adjust computed efflux velocities when ISPROPWASH = 2 (usual lower to address small scale turbulent eddies)
  real(kind = rkd) :: fraction_fast       = 0.         !< Fraction of the resuspended sediments from propwash that will have a faster settling speed than the standard class settling
  real(kind = rkd) :: fast_multiplier     = 30.        !< Fast class settling velocity multiplier 
  
  real(kind = rkd), parameter :: epsilon = 1.e-3       !< some threshold value
  real(kind = rkd), parameter :: pi_prop = 3.1415926535897932

  ! *** Timing variables
  integer(kind = rk4) :: nactiveships = 0              !< Number of active ships at current time
  integer(kind = rk4) :: iwrite_pwmesh = 0             !< Flag to indicate whether a PW_MESH file was written to for this time step
  real(kind = rkd)    :: freq_out_min = 0.             !< Minimum frequency output of PW_MESH files
  real(kind = rkd)    :: last_snapshot                 !< Timeday of the last snapshot written
  
  ! *** Total erosion
  real(kind = rkd), allocatable :: prop_ero(:,:)       !< Integrated erosion into suspended load over entire cell due to propwash from all ships, by sediment class (g/cm^2)
  real(kind = rkd), allocatable :: prop_bld(:,:)       !< Integrated erosion into bedload load over entire cell due to propwash from all ships, by sediment class (g/cm^2)

  integer :: IWC2BED(100)        !< Bed index for all WC cohesives
  integer :: IBED2WC(100)        !< WC index for all bed cohesives
  integer :: IFASTCLASS(100)     !< Flag to indicate whether a specific sediment class will be split into slow/fast settling
  
Contains
  

End Module Variables_Propwash   
