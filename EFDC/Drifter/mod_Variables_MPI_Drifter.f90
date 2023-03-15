! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Holds new variables for domain decomposition for the Drifter routines
! @author Zander Mausolff
! @date 9/4/0219

Module Variables_MPI_Drifter

    Use Global
    
    Implicit None
    
    Logical, Allocatable :: mask_lla(:)
  
    Integer :: local_tot_drifters_send
        
    ! *** identifies if a subdomain shoul be recieving data for the drifter routines
    Integer :: send_east  !< Number of drifters that should be sent in the  east  direction
    Integer :: send_west  !< Number of drifters that should be sent in the  west  direction
    Integer :: send_north !< Number of drifters that should be sent in the  north direction
    Integer :: send_south !< Number of drifters that should be sent in the  south direction
    
    Integer :: recv_east  !< Number of drifters that should be recvd in from the east  direction
    Integer :: recv_west  !< Number of drifters that should be recvd in from the west  direction
    Integer :: recv_north !< Number of drifters that should be recvd in from the north direction
    Integer :: recv_south !< Number of drifters that should be recvd in from the south direction
    
    Integer, Allocatable :: LLA_Process(:)
    
    Integer :: num_drifters_send_west !< keeps track of the number of drifters to send West
    Integer :: num_drifters_send_east
    Integer :: num_drifters_send_north
    Integer :: num_drifters_send_south
    
    Integer :: global_max_drifters_to_comm !< keeps track of the max number of drifters to keep track of amongs all processes
    
    Integer, Allocatable :: drifter_ids_send_west(:) !< Contains the drifters that will be communicated from the 
                                                     !! current domain to the the respective direction indicated 
                                                     !! by this variable name.  E.g. this array will send from the current --> West subdomain
    Integer, Allocatable :: drifter_ids_send_east(:)
    Integer, Allocatable :: drifter_ids_send_north(:)
    Integer, Allocatable :: drifter_ids_send_south(:)

    Integer, Allocatable :: drifter_ids_recv_west (:)
    Integer, Allocatable :: drifter_ids_recv_east (:)
    Integer, Allocatable :: drifter_ids_recv_north(:)
    Integer, Allocatable :: drifter_ids_recv_south(:)
    
    Integer, Allocatable :: gl_drift_mapping(:)
    
    ! *** Global data read in from the input file
    REAL(RKD), SAVE, ALLOCATABLE :: LA_BEGTI_Global(:)  
    REAL(RKD), SAVE, ALLOCATABLE :: LA_ENDTI_Global(:)  
    Integer, Allocatable         :: LLA_Global(:)
    INTEGER, SAVE, ALLOCATABLE   :: LA_GRP_Global(:)
                                 

End Module Variables_MPI_Drifter
