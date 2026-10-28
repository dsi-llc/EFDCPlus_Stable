! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Holds new variables for domain decomposition for the Drifter routines
! @author Zander Mausolff
! @date 9/4/0219

Module Variables_MPI_Drifter

    use GLOBAL
    
    implicit none
    
    logical, Allocatable :: mask_lla(:)
  
    integer :: local_tot_drifters_send
        
    ! *** identifies if a subdomain shoul be recieving data for the drifter routines
    integer :: send_east  !< Number of drifters that should be sent in the  east  direction
    integer :: send_west  !< Number of drifters that should be sent in the  west  direction
    integer :: send_north !< Number of drifters that should be sent in the  north direction
    integer :: send_south !< Number of drifters that should be sent in the  south direction
    
    integer :: recv_east  !< Number of drifters that should be recvd in from the east  direction
    integer :: recv_west  !< Number of drifters that should be recvd in from the west  direction
    integer :: recv_north !< Number of drifters that should be recvd in from the north direction
    integer :: recv_south !< Number of drifters that should be recvd in from the south direction
    
    integer, Allocatable :: LLA_Process(:)
    
    integer :: num_drifters_send_west !< keeps track of the number of drifters to send West
    integer :: num_drifters_send_east
    integer :: num_drifters_send_north
    integer :: num_drifters_send_south
    
    integer :: global_max_drifters_to_comm !< keeps track of the max number of drifters to keep track of amongs all processes
    
    integer, Allocatable :: drifter_ids_send_west(:) !< Contains the drifters that will be communicated from the 
                                                     !! current domain to the the respective direction indicated 
                                                     !! by this variable name.  E.g. this array will send from the current --> West subdomain
    integer, Allocatable :: drifter_ids_send_east(:)
    integer, Allocatable :: drifter_ids_send_north(:)
    integer, Allocatable :: drifter_ids_send_south(:)

    integer, Allocatable :: drifter_ids_recv_west (:)
    integer, Allocatable :: drifter_ids_recv_east (:)
    integer, Allocatable :: drifter_ids_recv_north(:)
    integer, Allocatable :: drifter_ids_recv_south(:)
    
    integer, Allocatable :: gl_drift_mapping(:)
    
    ! *** Global data read in from the input file
    real(RKD), save, allocatable :: LA_BEGTI_Global(:)  
    real(RKD), save, allocatable :: LA_ENDTI_Global(:)  
    integer, Allocatable         :: LLA_Global(:)
    integer, save, allocatable   :: LA_GRP_Global(:)
                                 

End Module Variables_MPI_Drifter
