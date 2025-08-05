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
! Module: Mod_Gather_Soln
!
!> @details 
!
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 1/8/2020
!---------------------------------------------------------------------------!
    
Module Mod_Gather_Soln

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    
    implicit none

    save
    
    Public :: Gather_Soln
    
    interface Gather_Soln
    
        Module Procedure Gather_1D_Real, &
                         Gather_1D_Real_RK8, &
                         Gather_1D_Int,  &
                         Gather_2D_Real, &
                         Gather_2D_Real_RK8, &
                         Gather_2D_Int,  &
                         Gather_3D_Real, &
                         Gather_3D_Real_RK8, &
                         Gather_3D_Int
    
        
    end interface

    Contains
    
!---------------------------------------------------------------------------!
! @details Gather entire spatial solution for 1D array of Reals
!
!> @param[in] size_local_1d: size of local soln and corresponding mapping arrays
!> @param[in] Soln_Local_1D: collapsed 1d solution excluding ghost cells
!> @param[in] size_global_1d: size of the global solution excluding ghost cells
!> @param[inout] Soln_Global_1D: global solution to be 
!> @param[in] Mapping_Local: mapping of local L to global L values
!> @param[inout] Mapping_Global: contains local to global mapping from each subdomain
!    
! @author Zander Mausolff
! @date 1/8/2020
!---------------------------------------------------------------------------!

Subroutine Gather_1D_Real(size_local_1d,  Soln_Local_1D, &
                          size_global_1d, Soln_Global_1D, Mapping_Local, Mapping_Global)

    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer, Intent(in)    :: size_local_1d 
    real(4),    Intent(in)    :: Soln_Local_1D(size_local_1d )
    integer, Intent(in)    :: size_global_1d
    real(4),    Intent(inout) :: Soln_Global_1D(size_global_1d) 
    integer, Intent(in)    :: Mapping_Local(size_local_1d)
    integer, Intent(inout) :: Mapping_Global(size_global_1d)
    
    !***Local variables
    integer :: ierr
    integer :: i
    
    If(.not.allocated(recv_counts_1d) )then
      allocate(recv_counts_1d(active_domains))
    endif
    
    !***we only want to send the number of active cells in the 1D array
    send_size_1d = num_active_l_local ! we do not want to send the ghost cells.. so ommit them.
    
    recv_counts_1d(:) = 0
    
    !***make sure each process has the same recv counts array for the 1d case
    call MPI_Allgather(send_size_1d,   1, MPI_Int, recv_counts_1d, 1, MPI_Int, DSIcomm, ierr)
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_1d,   mpi_real, &                                        ! Send buff
                     Soln_Global_1D, recv_counts_1d, displacements_L_index, mpi_real, master_id,  &     ! Recv buff
                     DSIcomm, ierr)
    
    !***Gather the mapping onto a global array
    call MPI_Gatherv(Mapping_Local,  send_size_1d,   MPI_Integer, &                                     ! Send buff
                     Mapping_Global, recv_counts_1d, displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
    
End Subroutine Gather_1D_Real

! *****************************************************************************************************
! *****************************************************************************************************
!---------------------------------------------------------------------------!
! @details Gather entire spatial solution for 1D array of Reals
!
!> @param[in] size_local_1d: size of local soln and corresponding mapping arrays
!> @param[in] Soln_Local_1D: collapsed 1d solution excluding ghost cells
!> @param[in] size_global_1d: size of the global solution excluding ghost cells
!> @param[inout] Soln_Global_1D: global solution to be 
!> @param[in] Mapping_Local: mapping of local L to global L values
!> @param[inout] Mapping_Global: contains local to global mapping from each subdomain
!    
! @author Zander Mausolff
! @date 1/8/2020
!---------------------------------------------------------------------------!

Subroutine Gather_1D_Real_RK8(size_local_1d,  Soln_Local_1D, &
                              size_global_1d, Soln_Global_1D, Mapping_Local, Mapping_Global)

    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer,  Intent(in)    :: size_local_1d 
    real(rkd),Intent(in)    :: Soln_Local_1D(size_local_1d )
    integer,  Intent(in)    :: size_global_1d
    real(rkd),Intent(inout) :: Soln_Global_1D(size_global_1d) 
    integer,  Intent(in)    :: Mapping_Local(size_local_1d)
    integer,  Intent(inout) :: Mapping_Global(size_global_1d)
    
    !***Local variables
    integer :: ierr
    integer :: i
    
    If(.not.allocated(recv_counts_1d) )then
      allocate(recv_counts_1d(active_domains))
    endif
    
    !***we only want to send the number of active cells in the 1D array
    send_size_1d = num_active_l_local ! we do not want to send the ghost cells.. so ommit them.
    
    recv_counts_1d(:) = 0
    
    !***make sure each process has the same recv counts array for the 1d case
    call MPI_Allgather(send_size_1d,   1, MPI_Int, recv_counts_1d, 1, MPI_Int, DSIcomm, ierr)
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_1d,   mpi_real8, &                                       ! Send buff
                     Soln_Global_1D, recv_counts_1d, displacements_L_index, mpi_real8, master_id,  &    ! Recv buff
                     DSIcomm, ierr)
    
    !***Gather the mapping onto a global array
    call MPI_Gatherv(Mapping_Local,  send_size_1d,   MPI_Integer, &                                     ! Send buff
                     Mapping_Global, recv_counts_1d, displacements_l_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
        
 End Subroutine Gather_1D_Real_RK8
!---------------------------------------------------------------------------!
! @details Gather entire spatial solution for 1D array of Integers
!
!> @param[in] size_local_1d: size of local soln and corresponding mapping arrays
!> @param[in] Soln_Local_1D: collapsed 1d solution excluding ghost cells
!> @param[in] size_global_1d: size of the global solution excluding ghost cells
!> @param[inout] Soln_Global_1D: global solution to be 
!> @param[in] Mapping_Local: mapping of local L to global L values
!> @param[inout] Mapping_Global: contains local to global mapping from each subdomain
!    
! @author Zander Mausolff
! @date 1/8/2020
!---------------------------------------------------------------------------!
Subroutine Gather_1D_int(size_local_1d,  Soln_Local_1D, &
                         size_global_1d, Soln_Global_1D, Mapping_Local, Mapping_Global)

    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Passed in
    integer, Intent(in)    :: size_local_1d 
    integer, Intent(in)    :: Soln_Local_1D(size_local_1d )
    integer, Intent(in)    :: size_global_1d
    integer, Intent(inout) :: Soln_Global_1D(size_global_1d)
    integer, Intent(in)    :: Mapping_Local(size_local_1d)
    integer, Intent(inout) :: Mapping_Global(size_global_1d)
    
    !***Local variables
    integer :: ierr
    integer :: i
    
    If(.not.allocated(recv_counts_1d) )then
      allocate(recv_counts_1d(active_domains))
    endif
    
    !***we only want to send the number of active cells in the 1D array
    send_size_1d = num_active_l_local ! we do not want to send the ghost cells.. so ommit them.
    
    recv_counts_1d(:) = 0
    
    !***make sure each process has the same recv counts array for the 1d case
    call MPI_Allgather(send_size_1d,   1, MPI_Int, recv_counts_1d, 1, MPI_Int, DSIcomm, ierr)
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_1d,   MPI_Integer, &                                     ! Send buff
                     Soln_Global_1D, recv_counts_1d, displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
    
    !***Gather the mapping onto a global array
    call MPI_Gatherv(Mapping_Local,  send_size_1d,   MPI_Integer, &                                     ! Send buff
                     Mapping_Global, recv_counts_1d, displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
    
End Subroutine Gather_1D_int
                         
!---------------------------------------------------------------------------!
! @details Gather entire spatial solution for (L,H) arrays
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size: 
!> @param[in] Soln_Local_1D: 
!> @param[in] first_dim_size_global
!> @param[in] second_dim_size_1D
!> @param[inout] Soln_Global_1D: 
!> @param[in] Mapping_Local: mapping of local L to global L values
!> @param[inout] Mapping_Global: contains local to global mapping from each subdomain
!    
! @author Zander Mausolff
! @date 1/8/2020
!
!---------------------------------------------------------------------------!

    Subroutine Gather_2D_Real(first_dim_size,        second_dim_size,    Soln_Local_1D, &
                              first_dim_size_global, second_dim_size_1D, Soln_Global_1D, Mapping_Local, Mapping_Global)
    
    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer, Intent(in)    :: first_dim_size
    integer, Intent(in)    :: second_dim_size
    real(4), intent(in)    :: Soln_Local_1D(first_dim_size*second_dim_size)
    integer, intent(in)    :: first_dim_size_global
    integer, Intent(in)    :: second_dim_size_1D
    real(4), intent(inout) :: Soln_Global_1D(first_dim_size_global*second_dim_size_1D) !< assumes LCM global sizing
    integer, intent(in)    :: Mapping_Local(first_dim_size*second_dim_size_1D)
    integer, intent(inout) :: Mapping_Global(first_dim_size_global*second_dim_size_1D)
    
    !***Local variables
    integer :: ierr
    integer :: send_size_3d  !< Size of message for communicating the entire domain for (L,K) arrays
    integer :: i
    
    If(.not.allocated(recv_counts_3d) )then
      allocate(recv_counts_3d(active_domains))
    endif
    
    !***we only want to send the number of active cells in the 1D array
    send_size_3d = num_active_l_local! we do not want to send the ghost cells.. so ommit them.
    
    recv_counts_3d(:) = 0

    !***make sure each process has the same recv counts array for the 3D case
    call MPI_Allgather(send_size_3d,   1, MPI_Int, recv_counts_3d, 1, MPI_Int, DSIcomm, ierr)
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D, send_size_3d, mpi_real, &                                            ! Send buff
                     Soln_Global_1D, recv_counts_3d , displacements_L_index, mpi_real, master_id,  &     ! Recv buff
                     DSIcomm, ierr)
    
    !***Gather the mapping onto a global array
    call MPI_Gatherv(Mapping_Local, send_size_3d, MPI_Integer, &                                         ! Send buff
                     Mapping_Global, recv_counts_3d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
    
End Subroutine Gather_2D_Real

!---------------------------------------------------------------------------!
! @details Gather entire spatial solution for (L,H) arrays
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size: 
!> @param[in] Soln_Local_1D: 
!> @param[in] first_dim_size_global
!> @param[in] second_dim_size_1D
!> @param[inout] Soln_Global_1D: 
!> @param[in] Mapping_Local: mapping of local L to global L values
!> @param[inout] Mapping_Global: contains local to global mapping from each subdomain
!    
! @author Zander Mausolff
! @date 1/8/2020
!
!---------------------------------------------------------------------------!

    Subroutine Gather_2D_Real_RK8(first_dim_size,        second_dim_size,    Soln_Local_1D, &
                                  first_dim_size_global, second_dim_size_1D, Soln_Global_1D, Mapping_Local, Mapping_Global)
    
    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer,  Intent(in)    :: first_dim_size
    integer,  Intent(in)    :: second_dim_size
    real(rkd),Intent(in)    :: Soln_Local_1D(first_dim_size*second_dim_size)
    integer,  Intent(in)    :: first_dim_size_global
    integer,  Intent(in)    :: second_dim_size_1D
    real(rkd),Intent(inout) :: Soln_Global_1D(first_dim_size_global*second_dim_size_1D) !< assumes LCM global sizing
    integer,  Intent(in)    :: Mapping_Local(first_dim_size*second_dim_size_1D)
    integer,  Intent(inout) :: Mapping_Global(first_dim_size_global*second_dim_size_1D)
    
    !***Local variables
    integer :: ierr
    integer :: send_size_3d  !< Size of message for communicating the entire domain for (L,K) arrays
    integer :: i
    
    If(.not.allocated(recv_counts_3d) )then
      allocate(recv_counts_3d(active_domains))
    endif
    
    !***we only want to send the number of active cells in the 1D array
    send_size_3d = num_active_l_local! we do not want to send the ghost cells.. so ommit them.
    
    recv_counts_3d(:) = 0

    !***make sure each process has the same recv counts array for the 3D case
    call MPI_Allgather(send_size_3d,   1, MPI_Int, &
                       recv_counts_3d, 1, MPI_Int, DSIcomm, ierr)
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D, send_size_3d, mpi_real8, &                                           ! Send buff
                     Soln_Global_1D, recv_counts_3d , displacements_L_index, mpi_real8, master_id,  &    ! Recv buff
                     DSIcomm, ierr)
    
    !***Gather the mapping onto a global array
    call MPI_Gatherv(Mapping_Local, send_size_3d, MPI_Integer, &                                         ! Send buff
                     Mapping_Global, recv_counts_3d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
    
End Subroutine Gather_2D_Real_RK8
    
!---------------------------------------------------------------------------!
! @details Gather entire spatial solution for (L,H) arrays
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size: 
!> @param[in] Soln_Local_1D: 
!> @param[in] first_dim_size_global
!> @param[in] second_dim_size_1D
!> @param[inout] Soln_Global_1D: 
!> @param[in] Mapping_Local: mapping of local L to global L values
!> @param[inout] Mapping_Global: contains local to global mapping from each subdomain
!    
! @author Zander Mausolff
! @date 1/8/2020
!
!---------------------------------------------------------------------------!

    Subroutine Gather_2D_Int(first_dim_size,        second_dim_size,    Soln_Local_1D, &
                             first_dim_size_global, second_dim_size_1D, Soln_Global_1D, Mapping_Local, Mapping_Global)
    
    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer, Intent(in)    :: first_dim_size
    integer, Intent(in)    :: second_dim_size
    integer, intent(in)    :: Soln_Local_1D(first_dim_size*second_dim_size)
    integer, intent(in)    :: first_dim_size_global
    integer, Intent(in)    :: second_dim_size_1D
    integer, intent(inout) :: Soln_Global_1D(first_dim_size_global*second_dim_size_1D) !< assumes LCM global sizing
    integer, intent(in)    :: Mapping_Local(first_dim_size*second_dim_size_1D)
    integer, intent(inout) :: Mapping_Global(first_dim_size_global*second_dim_size_1D)
    
    !***Local variables
    integer :: ierr
    integer :: send_size_3d  !< Size of message for communicating the entire domain for (L,K) arrays
    integer :: i
    
    If(.not.allocated(recv_counts_3d) )then
      allocate(recv_counts_3d(active_domains))
    endif
    
    !***we only want to send the number of active cells in the 1D array
    send_size_3d = num_active_l_local! we do not want to send the ghost cells.. so ommit them.
    
    recv_counts_3d(:) = 0

    !***make sure each process has the same recv counts array for the 3D case
    call MPI_Allgather(send_size_3d, 1, MPI_Int, recv_counts_3d, 1, MPI_Int, DSIcomm, ierr)
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_3d, MPI_Integer, &                                        ! Send buff
                     Soln_Global_1D, recv_counts_3d , displacements_L_index, MPI_Integer, master_id, &   ! Recv buff
                     DSIcomm, ierr)
                     
    !***Gather the mapping onto a global array
    call MPI_Gatherv(Mapping_Local,  send_size_3d, MPI_Integer, &                                        ! Send buff
                     Mapping_Global, recv_counts_3d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
    
    End Subroutine Gather_2D_Int
!---------------------------------------------------------------------------!  
! Subroutine: Gather_3D_Real
!
!> @details Gather entire spatial solution for (L,H,N) arrays
!
!> @param[in] first_dim_size
!> @param[in] second_dim_size
!> @param[in] third_dim_size
!> @param[in] Soln_Local_1D
!> @param[in] Mapping_Local
!> @param[in] first_dim_size_gl
!> @param[in] second_dim_size_gl
!> @param[in] third_dim_size_gl
!> @param[inout] Soln_Global_1D
!> @param[inout] Mapping_Global    
!
!> @author Zander Mausolff
! @date 1/8/2020
!---------------------------------------------------------------------------!
    
Subroutine Gather_3D_Real(first_dim_size,    second_dim_size,    third_dim_size,    Soln_Local_1D,  Mapping_Local, &
                          first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_1D, Mapping_Global)
   
    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer, intent(in)    :: first_dim_size
    integer, Intent(in)    :: second_dim_size
    integer, intent(in)    :: third_dim_size
    real(4), intent(in)    :: Soln_Local_1D(first_dim_size*second_dim_size*third_dim_size)
    integer, intent(in)    :: Mapping_Local(first_dim_size*second_dim_size*third_dim_size)
    integer, intent(in)    :: first_dim_size_gl
    integer, Intent(in)    :: second_dim_size_gl
    integer, intent(in)    :: third_dim_size_gl
    real(4), intent(inout) :: Soln_Global_1D(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl) !< assumes LCM global sizing
    integer, intent(inout) :: Mapping_Global(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl)
    
    !***Local variables
    integer :: ierr     !< local MPI error flag
    integer :: i, send_size_4d
    integer, Allocatable, Dimension(:) :: recv_counts_4d
    
    if(.not.Allocated(recv_counts_4d) )then
        allocate(recv_counts_4d(active_domains))
    endif
    
    !***Length of message size
    send_size_4d    = num_active_l_local ! we do not want to send the ghost cells.. so ommit them.

    !***make sure each process has the same recv counts array for the 3D case
    call MPI_Allgather(send_size_4d,   1, MPI_Int, &
                       recv_counts_4d, 1, MPI_Int, DSIcomm, ierr)
                       
                        
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_4d, mpi_real, &                                           ! Send buff
                     Soln_Global_1D, recv_counts_4d , displacements_L_index, mpi_real, master_id,  &     ! Recv buff
                     DSIcomm, ierr)
                     
    !***Gather the mapping onto a global array                   
    call MPI_Gatherv(Mapping_Local,  send_size_4d, MPI_Integer, &                                        ! Send buff
                     Mapping_Global, recv_counts_4d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)                  
                     
End Subroutine Gather_3D_Real
                          
!---------------------------------------------------------------------------!  
! Subroutine: Gather_3D_Real
!
!> @details Gather entire spatial solution for (L,H,N) arrays
!
!> @param[in] first_dim_size
!> @param[in] second_dim_size
!> @param[in] third_dim_size
!> @param[in] Soln_Local_1D
!> @param[in] Mapping_Local
!> @param[in] first_dim_size_gl
!> @param[in] second_dim_size_gl
!> @param[in] third_dim_size_gl
!> @param[inout] Soln_Global_1D
!> @param[inout] Mapping_Global    
!
!> @author Zander Mausolff
! @date 1/8/2020
!---------------------------------------------------------------------------!
    
Subroutine Gather_3D_Real_RK8(first_dim_size,    second_dim_size,    third_dim_size,    Soln_Local_1D,  Mapping_Local, &
                              first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_1D, Mapping_Global)
   
    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer,  intent(in)    :: first_dim_size
    integer,  Intent(in)    :: second_dim_size
    integer,  intent(in)    :: third_dim_size
    real(rkd),intent(in)    :: Soln_Local_1D(first_dim_size*second_dim_size*third_dim_size)
    integer,  intent(in)    :: Mapping_Local(first_dim_size*second_dim_size*third_dim_size)
    integer,  intent(in)    :: first_dim_size_gl
    integer,  Intent(in)    :: second_dim_size_gl
    integer,  intent(in)    :: third_dim_size_gl
    real(rkd),intent(inout) :: Soln_Global_1D(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl) !< assumes LCM global sizing
    integer,  intent(inout) :: Mapping_Global(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl)
    
    !***Local variables
    integer :: ierr     !< local MPI error flag
    integer :: i, send_size_4d
    integer, Allocatable, Dimension(:) :: recv_counts_4d
    
    if(.not.Allocated(recv_counts_4d) )then
        allocate(recv_counts_4d(active_domains))
    endif
    
    !***Length of message size
    send_size_4d    = num_active_l_local ! we do not want to send the ghost cells.. so ommit them.

    !***make sure each process has the same recv counts array for the 3D case
    call MPI_Allgather(send_size_4d,   1, MPI_Int, &
                       recv_counts_4d, 1, MPI_Int, DSIcomm, ierr)
                       
                        
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_4d, mpi_real8, &                                          ! Send buff
                     Soln_Global_1D, recv_counts_4d , displacements_L_index, mpi_real8, master_id,  &    ! Recv buff
                     DSIcomm, ierr)
                     
    !***Gather the mapping onto a global array                   
    call MPI_Gatherv(Mapping_Local,  send_size_4d, MPI_Integer, &                                        ! Send buff
                     Mapping_Global, recv_counts_4d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)                  
                     
End Subroutine Gather_3D_Real_RK8
                              
!---------------------------------------------------------------------------!  
! Subroutine: Gather_3D_Int
!
!> @details Gather entire spatial solution for (L,H,N) arrays
!
!> @param[in] first_dim_size
!> @param[in] second_dim_size
!> @param[in] third_dim_size
!> @param[in] Soln_Local_1D
!> @param[in] Mapping_Local
!> @param[in] first_dim_size_gl
!> @param[in] second_dim_size_gl
!> @param[in] third_dim_size_gl
!> @param[inout] Soln_Global_1D
!> @param[inout] Mapping_Global    
!
!> @author Zander Mausolff
! @date 1/8/2020
!---------------------------------------------------------------------------!
    
Subroutine Gather_3D_Int(first_dim_size,    second_dim_size,    third_dim_size,    Soln_Local_1D,  Mapping_Local, &
                         first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_1D, Mapping_Global)
   
    use GLOBAL
    use Variables_MPI
    use MPI
    
    implicit none
    
    !***Read in
    integer, intent(in)    :: first_dim_size
    integer, Intent(in)    :: second_dim_size
    integer, intent(in)    :: third_dim_size
    integer, intent(in)    :: Soln_Local_1D(first_dim_size*second_dim_size*third_dim_size)
    integer, intent(in)    :: Mapping_Local(first_dim_size*second_dim_size*third_dim_size)
    integer, intent(in)    :: first_dim_size_gl
    integer, Intent(in)    :: second_dim_size_gl
    integer, intent(in)    :: third_dim_size_gl
    integer, intent(inout) :: Soln_Global_1D(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl) !< assumes LCM global sizing
    integer, intent(inout) :: Mapping_Global(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl)
    
    !***Local variables
    integer :: ierr     !< local MPI error flag
    integer :: i, send_size_4d
    integer, Allocatable, Dimension(:) :: recv_counts_4d
    
    if(.not.Allocated(recv_counts_4d) )then
        allocate(recv_counts_4d(active_domains))
    endif
    
    !***Length of message size
    send_size_4d    = num_active_l_local ! we do not want to send the ghost cells.. so ommit them.

    !***make sure each process has the same recv counts array for the 3D case
    call MPI_Allgather(send_size_4d,   1, MPI_Int, recv_counts_4d, 1, MPI_Int, DSIcomm, ierr)
                        
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_4d, MPI_Integer, &                                        ! Send buff
                     Soln_Global_1D, recv_counts_4d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)
                     
    !***Gather the mapping onto a global array                   
    call MPI_Gatherv(Mapping_Local,  send_size_4d, MPI_Integer, &                                        ! Send buff
                     Mapping_Global, recv_counts_4d , displacements_L_index, MPI_Integer, master_id,  &  ! Recv buff
                     DSIcomm, ierr)                  
                     
End Subroutine Gather_3D_Int

End Module Mod_Gather_Soln
