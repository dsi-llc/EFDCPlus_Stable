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
! @details Creates an array of L values for only the active cells, excluding 
!! ghost cells
! @author Paul Craig/Zander Mausolff
! @date 1/6/2019
!---------------------------------------------------------------------------!  
Subroutine Create_List_No_Ghost_Cells
    
  use GLOBAL
  use Variables_MPI_Mapping
  use MPI

  implicit none

  ! *** Local variables
  integer :: I, J, L

  if(.not.allocated(NoGhost) )then
    allocate(NoGhost(LCM))
    allocate(IsGhost(LCM))
    allocate(GhostMask(LCM))
  endif

  ! *** Save the list of active cells, excluding ghost cells
  IsGhost = .TRUE.
  GhostMask = 0.0
  NoGhost = 0
  NNoGhost = 0
  do J = 3,JC-2
    do I = 3,IC-2
      L = LIJ(I,J)
      if( L > 0 )then ! *** Only records the active cells
        NNoGhost = NNoGhost + 1
        NoGhost(NNoGhost) = L
        IsGhost(L) = .FALSE.
        GhostMask(L) = 1.0
      endif
    enddo
  enddo
  
  ! *** Get the number of active cells excluding anything in the ghost regions
  LA_Local_no_ghost = NNoGhost
  
End Subroutine Create_List_No_Ghost_Cells
