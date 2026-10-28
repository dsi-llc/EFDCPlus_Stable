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
! Module: Mod_Map_Soln
!
!> @details Module that contains interface for mapping active cells in a soln
!!  and collapsing to a 1D array.  It also determines how the  
!
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 1/8/2020
!---------------------------------------------------------------------------!
    
Module Mod_Map_Soln

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    
    implicit none

    save
    
    Public :: Map_Soln
    
    interface Map_Soln
    
        Module Procedure Map_1D_Real, &
                         Map_1D_Real_RK8, &
                         Map_1D_Int,  &
                         Map_2D_Real, &
                         Map_2D_Real_RK8, &
                         Map_2D_Int,  &
                         Map_3D_Real, &
                         Map_3D_Real_RK8, &
                         Map_3D_Int
    
    end interface

    Contains
    
!---------------------------------------------------------------------------!
! Subroutine: Map_1D_Real
!
!> @details Maps local values to 1 contiguous 1d array of REAL elements
!
!> @param[in] size_soln_local: size of the local Soln (containing ghost cells)
!> @param[in] Soln_Local: Local 1D solution array
!> @param[in] size_soln_local_1D: size of  1D soln (without ghost cells)
!> @param[inout] Soln_Local_1D: Local 1D solution excluding ghost cells
!> @param[inout] Map_Local_to_Global_1D: global mapping of local L values
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 9/10/2019
!---------------------------------------------------------------------------!
    Subroutine Map_1D_Real(size_soln_local, Soln_Local,&
                           size_soln_local_1D, Soln_Local_1D, &
                           Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI

    implicit none

    ! *** Passed in variables
    integer, Intent(in)    :: size_soln_local
    real(4), Intent(in)       :: Soln_Local(size_soln_local) 
    integer, Intent(in)    :: size_soln_local_1D
    real(4), Intent(inout)    :: Soln_Local_1D(size_soln_local_1D)            
    integer, Intent(inout) :: Map_Local_to_Global_1D(size_soln_local_1D)   

    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process

    if(.not.Allocated(all_local_LA_array) )then
      allocate(all_local_LA_array(active_domains))
    endif

    If(.not.allocated(displacements_L_index) )then
      allocate(displacements_L_index(active_domains))
    Endif

    all_local_LA_array(:) = 0
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
      ! *** Get first 2 rows
      do i = 3, ic - 2
        ! *** Map to the correct l index
        l  = lij(i,j)
        
        ! *** Get global L information
        i_gl = IL2IG(i)
        j_gl = JL2JG(j)
        l_gl = lij_global(i_gl,j_gl)
        
        ! *** Build up array containing active cells excluding ghost cells
        if( l > 0 )then
          iii = iii + 1
          Soln_Local_1D(iii) = Soln_Local(l)   ! *** Setup 1D aray of the solution
          Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
        endif
      enddo
    enddo

    num_active_l_local = iii        ! ***  iii is the count the number of active cells, excluding ghost
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, all_local_LA_array, 1, MPI_Int, DSIcomm, ierr)

    displacements_L_index(1) = 0 ! Want to start at L = 2
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
      displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo

    End Subroutine Map_1D_Real
! *****************************************************************************************************
! *****************************************************************************************************
    !---------------------------------------------------------------------------!
! Subroutine: Map_1D_Real
!
!> @details Maps local values to 1 contiguous 1d array of REAL elements
!
!> @param[in] size_soln_local: size of the local Soln (containing ghost cells)
!> @param[in] Soln_Local: Local 1D solution array
!> @param[in] size_soln_local_1D: size of  1D soln (without ghost cells)
!> @param[inout] Soln_Local_1D: Local 1D solution excluding ghost cells
!> @param[inout] Map_Local_to_Global_1D: global mapping of local L values
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 9/10/2019
!---------------------------------------------------------------------------!
    Subroutine Map_1D_Real_RK8(size_soln_local, Soln_Local,&
                           size_soln_local_1D, Soln_Local_1D, &
                           Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI

    implicit none

    ! *** Passed in variables
    integer, Intent(in)    :: size_soln_local
    real(rkd), Intent(in)       :: Soln_Local(size_soln_local) 
    integer, Intent(in)    :: size_soln_local_1D
    real(rkd), Intent(inout)    :: Soln_Local_1D(size_soln_local_1D)            
    integer, Intent(inout) :: Map_Local_to_Global_1D(size_soln_local_1D)   

    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process

    if(.not.Allocated(all_local_LA_array) )then
      allocate(all_local_LA_array(active_domains))
    endif

    If(.not.allocated(displacements_L_index) )then
      allocate(displacements_L_index(active_domains))
    Endif

    all_local_LA_array(:) = 0
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
      ! *** Get first 2 rows
      do i = 3, ic - 2
        ! *** Map to the correct l index
        l  = lij(i,j)
        ! *** Get global L information
        i_gl = IL2IG(i)
        j_gl = JL2JG(j)
        l_gl = lij_global(i_gl,j_gl)
        ! *** Build up array containing active cells excluding ghost cells
        if( l > 0 )then
          iii = iii + 1
          Soln_Local_1D(iii) = Soln_Local(l) !***Setup 1D aray of the solution
          Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
        endif
      enddo
    enddo

    num_active_l_local = iii
    ! ***  iii will count the number of

    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, all_local_LA_array, 1, MPI_Int, DSIcomm, ierr)

    displacements_L_index(1) = 0 ! Want to start at L = 2
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
      displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo

    End Subroutine Map_1D_Real_RK8
!---------------------------------------------------------------------------!
! Subroutine: Map_1D_Int
!
!> @details Maps local values to 1 contiguous 1d array of REAL elements
!
!> @param[in] size_soln_local: size of the local Soln (containing ghost cells)
!> @param[in] Soln_Local: Local 1D solution array
!> @param[in] size_soln_local_1D: size of  1D soln (without ghost cells)
!> @param[inout] Soln_Local_1D: Local 1D solution excluding ghost cells
!> @param[inout] Map_Local_to_Global_1D: global mapping of local L values
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!  
    Subroutine Map_1D_Int(size_soln_local, Soln_Local,&
                           size_soln_local_1D, Soln_Local_1D, &
                           Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI
    
    implicit none
    
    ! *** Passed in variables
    integer, Intent(in)    :: size_soln_local
    integer, Intent(in)    :: Soln_Local(size_soln_local) 
    integer, Intent(in)    :: size_soln_local_1D
    integer, Intent(inout) :: Soln_Local_1D(size_soln_local_1D)           
    integer, Intent(inout) :: Map_Local_to_Global_1D(size_soln_local_1D)   
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
      allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
      allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
      ! *** Get first 2 rows
      do i = 3, ic - 2
        ! *** Map to the correct l index
        l  = lij(i,j)
        ! *** Get global L information
        i_gl = IL2IG(i)
        j_gl = JL2JG(j)
        l_gl = lij_global(i_gl,j_gl)
        ! *** Build up array containing active cells excluding ghost cells
        if( l > 0 )then
          iii = iii + 1
          Soln_Local_1D(iii) = Soln_Local(l)   ! *** Setup 1D aray of the solution
          Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
        endif
      enddo
    enddo
    
    num_active_l_local = iii
    ! ***  iii will count the number of
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, all_local_LA_array, 1, MPI_Int, DSIcomm, ierr)
    
    displacements_L_index(1) = 0 ! Want to start at L = 2
    
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
      displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo

    End Subroutine Map_1D_Int

! *****************************************************************************************************
! *****************************************************************************************************
!---------------------------------------------------------------------------!  
! Subroutine: Map_2D_real
!
!> @details Maps local values to 1 contiguous 1d array 
!
!> @param[in] first_dim_size_local
!> @param[in] second_dim_size_local
!> @param[in] Soln_Local
!> @param[in] second_sim_size_1D
!> @param[inout] Soln_Local_1D
!> @param[inout] Soln_Global: Global array containing entire spatial solution
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 9/10/2019
!---------------------------------------------------------------------------!
Subroutine Map_2D_Real(first_dim_size_local, second_dim_size_local, Soln_Local, &
					   second_sim_size_1D, Soln_Local_1D, Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping
    use MPI
    
    implicit none
    
    ! *** Passed in variables
    integer, Intent(in)    :: first_dim_size_local
    integer, Intent(in)    :: second_dim_size_local
    real(4), Intent(in)    :: Soln_Local(first_dim_size_local, second_dim_size_local)      !< spatial solution local to a process
    integer, Intent(in)    :: second_sim_size_1D
    real(4), Intent(inout) :: Soln_Local_1D(LA_Local_no_ghost*second_sim_size_1D)          !< global solution for entire spatial domain
    integer, Intent(inout) :: Map_Local_to_Global_1D(LA_Local_no_ghost*second_sim_size_1D) !< mapping of local to global values
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl 
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
        allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
        allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
        ! *** Get first 2 rows
        do i = 3, ic - 2
        ! *** Grab the layers
            do k = 1, second_sim_size_1D
                ! *** Map to the correct l index
                l  = lij(i,j)
                ! *** Get global L information
                i_gl = IL2IG(i)
                j_gl = JL2JG(j)
                l_gl = lij_global(i_gl,j_gl)
                ! *** Build up array containing active cells excluding ghost cells
                if( l > 0 )then
                    iii = iii + 1
                    Soln_Local_1D(iii) = Soln_Local(l,k) ! *** Setup 1D aray of the solution
                    Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
                endif
            enddo
        enddo
    enddo 
    
    num_active_l_local = iii
    ! ***  iii will count the number of 
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, &
                       all_local_LA_array, 1, MPI_Int, &
                       DSIcomm, ierr)

    displacements_L_index(1) = 0 ! Want to start at L = 2
    
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
        displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo
    
  End Subroutine Map_2D_Real
! *****************************************************************************************************
! *****************************************************************************************************
!---------------------------------------------------------------------------!  
! Subroutine: Map_2D_real
!
!> @details Maps local values to 1 contiguous 1d array 
!
!> @param[in] first_dim_size_local
!> @param[in] second_dim_size_local
!> @param[in] Soln_Local
!> @param[in] second_sim_size_1D
!> @param[inout] Soln_Local_1D
!> @param[inout] Soln_Global: Global array containing entire spatial solution
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 9/10/2019
!---------------------------------------------------------------------------!
Subroutine Map_2D_Real_RK8(first_dim_size_local, second_dim_size_local, Soln_Local, &
					   second_sim_size_1D, Soln_Local_1D, Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping
    use MPI
    
    implicit none
    
    ! *** Passed in variables
    integer, Intent(in)    :: first_dim_size_local
    integer, Intent(in)    :: second_dim_size_local
    real(rkd), Intent(in)       :: Soln_Local(first_dim_size_local, second_dim_size_local)    !< spatial solution local to a process
    integer, Intent(in)    :: second_sim_size_1D
    real(rkd), Intent(inout)    :: Soln_Local_1D(LA_Local_no_ghost*second_sim_size_1D)    !< global solution for entire spatial domain
    integer, Intent(inout) :: Map_Local_to_Global_1D(LA_Local_no_ghost*second_sim_size_1D) !< mapping of local to global values
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl 
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
        allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
        allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
        ! *** Get first 2 rows
        do i = 3, ic - 2
        ! *** Grab the layers
            do k = 1, second_sim_size_1D
                ! *** Map to the correct l index
                l  = lij(i,j)
                ! *** Get global L information
                i_gl = IL2IG(i)
                j_gl = JL2JG(j)
                l_gl = lij_global(i_gl,j_gl)
                ! *** Build up array containing active cells excluding ghost cells
                if( l > 0 )then
                    iii = iii + 1
                    Soln_Local_1D(iii) = Soln_Local(l,k) ! *** Setup 1D aray of the solution
                    Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
                endif
            enddo
        enddo
    enddo 
    
    num_active_l_local = iii
    ! ***  iii will count the number of 
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, &
                       all_local_LA_array, 1, MPI_Int, &
                       DSIcomm, ierr)

    displacements_L_index(1) = 0 ! Want to start at L = 2
    
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
        displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo
    

    End Subroutine Map_2D_Real_RK8
!---------------------------------------------------------------------------!  
! Subroutine: Map_2D_Int
!
!> @details Maps local values to 1 contiguous 1d array 
!
!> @param[in] first_dim_size_local
!> @param[in] second_dim_size_local
!> @param[in] Soln_Local
!> @param[in] second_sim_size_1D
!> @param[inout] Soln_Local_1D
!> @param[inout] Soln_Global: Global array containing entire spatial solution
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 9/10/2019
!---------------------------------------------------------------------------!
Subroutine Map_2D_Int(first_dim_size_local, second_dim_size_local, Soln_Local, &
					   second_sim_size_1D, Soln_Local_1D, Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI
    
    implicit none
    
    ! *** Passed in variables
    integer, Intent(in)    :: first_dim_size_local
    integer, Intent(in)    :: second_dim_size_local
    integer, Intent(in)    :: Soln_Local(first_dim_size_local, second_dim_size_local)      !< spatial solution local to a process
    integer, Intent(in)    :: second_sim_size_1D
    integer, Intent(inout) :: Soln_Local_1D(LA_Local_no_ghost*second_sim_size_1D)          !< global solution for entire spatial domain
    integer, Intent(inout) :: Map_Local_to_Global_1D(LA_Local_no_ghost*second_sim_size_1D) !< mapping of local to global values
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl 
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
        allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
        allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
        ! *** Get first 2 rows
        do i = 3, ic - 2
        ! *** Grab the layers
            do k = 1, second_sim_size_1D
                ! *** Map to the correct l index
                l  = lij(i,j)
                ! *** Get global L information
                i_gl = IL2IG(i)
                j_gl = JL2JG(j)
                l_gl = lij_global(i_gl,j_gl)
                ! *** Build up array containing active cells excluding ghost cells
                if( l > 0 )then
                    iii = iii + 1
                    Soln_Local_1D(iii) = Soln_Local(l,k) ! *** Setup 1D aray of the solution
                    Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
                endif
            enddo
        enddo
    enddo 
    
    num_active_l_local = iii
    ! ***  iii will count the number of 
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, &
                       all_local_LA_array, 1, MPI_Int, &
                       DSIcomm, ierr)

    displacements_L_index(1) = 0 ! Want to start at L = 2
    
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
        displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo
    

    End Subroutine Map_2D_Int
    
!---------------------------------------------------------------------------!  
! Subroutine: Map_3D_Real
!
!> @details Maps local values to 1 contiguous 1d array 
!
!> @param[in]  first_dim_size: 
!> @param[in]  second_dim_size: 
!> @param[in]  third_dim_size: 
!> @param[in]  Soln_Local: 
!> @param[in[  second_dim_size_1D: 
!> @param[in]  third_dim_size_1D:
!> @param[inout] Soln_Local_1D: 
!> @param[inout] Soln_Global: Global array containing entire spatial solution
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 1/8/2020
!---------------------------------------------------------------------------!
Subroutine Map_3D_Real(first_dim_size, second_dim_size, third_dim_size, Soln_Local, &
	       second_dim_size_1D, third_dim_size_1D, Soln_Local_1D, Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI
    
    implicit none
    
    ! *** Passed in variables

    integer, Intent(in)    :: first_dim_size  
    integer, Intent(in)    :: second_dim_size
    integer, Intent(in)    :: third_dim_size
    real(4), Intent(in)       :: Soln_Local(first_dim_size, second_dim_size, third_dim_size) 
    integer, Intent(in)    :: second_dim_size_1D
    integer, Intent(in)    :: third_dim_size_1D
    real(4), Intent(inout)    :: Soln_Local_1D(LA_Local_no_ghost*second_dim_size_1D*third_dim_size_1D) 
    integer, Intent(inout) :: Map_Local_to_Global_1D(LA_Local_no_ghost*second_dim_size_1D*third_dim_size_1D) 
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl, nn
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
        allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
        allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
        ! *** Get first 2 rows
        do i = 3, ic - 2
        ! *** Grab the layers
            ! *** Map to the correct l index
            l  = lij(i,j)
            ! *** Get global L information
            i_gl = IL2IG(i)
            j_gl = JL2JG(j)
            l_gl = lij_global(i_gl,j_gl)
            do k = 1, second_dim_size_1D
                do nn = 1, third_dim_size_1D
                    ! *** Build up array containing active cells excluding ghost cells
                    if( l > 0 )then
                        iii = iii + 1
                        Soln_Local_1D(iii) = Soln_Local(l,k,nn) ! *** Setup 1D aray of the solution
                        Map_Local_to_Global_1D(iii) = l_gl      ! *** Get corresponding global l value for a given local
                    endif
                enddo ! n loop
            enddo ! k loop
        enddo ! i loop
    enddo ! j loop
    
    num_active_l_local = iii
    ! ***  iii will count the number of 
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, &
                       all_local_LA_array, 1, MPI_Int, &
                       DSIcomm, ierr)
    
    displacements_L_index(1) = 0 ! Want to start at L = 2
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
        displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo


End Subroutine Map_3D_Real
!---------------------------------------------------------------------------!  
! Subroutine: Map_3D_Real
!
!> @details Maps local values to 1 contiguous 1d array 
!
!> @param[in]  first_dim_size: 
!> @param[in]  second_dim_size: 
!> @param[in]  third_dim_size: 
!> @param[in]  Soln_Local: 
!> @param[in[  second_dim_size_1D: 
!> @param[in]  third_dim_size_1D:
!> @param[inout] Soln_Local_1D: 
!> @param[inout] Soln_Global: Global array containing entire spatial solution
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 1/8/2020
!---------------------------------------------------------------------------!
Subroutine Map_3D_Real_RK8(first_dim_size, second_dim_size, third_dim_size, Soln_Local, &
	       second_dim_size_1D, third_dim_size_1D, Soln_Local_1D, Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI
    
    implicit none
    
    ! *** Passed in variables

    integer, Intent(in)    :: first_dim_size  
    integer, Intent(in)    :: second_dim_size
    integer, Intent(in)    :: third_dim_size
    real(rkd), Intent(in)       :: Soln_Local(first_dim_size, second_dim_size, third_dim_size) 
    integer, Intent(in)    :: second_dim_size_1D
    integer, Intent(in)    :: third_dim_size_1D
    real(rkd), Intent(inout)    :: Soln_Local_1D(LA_Local_no_ghost*second_dim_size_1D*third_dim_size_1D) 
    integer, Intent(inout) :: Map_Local_to_Global_1D(LA_Local_no_ghost*second_dim_size_1D*third_dim_size_1D) 
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl, nn
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
        allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
        allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
        ! *** Get first 2 rows
        do i = 3, ic - 2
        ! *** Grab the layers
            ! *** Map to the correct l index
            l  = lij(i,j)
            ! *** Get global L information
            i_gl = IL2IG(i)
            j_gl = JL2JG(j)
            l_gl = lij_global(i_gl,j_gl)
            do k = 1, second_dim_size_1D
                do nn = 1, third_dim_size_1D
                    ! *** Build up array containing active cells excluding ghost cells
                    if( l > 0 )then
						iii = iii + 1
                        Soln_Local_1D(iii) = Soln_Local(l,k,nn) !***Setup 1D aray of the solution
                        Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
                    endif
                enddo ! n loop
            enddo ! k loop
        enddo ! i loop
    enddo ! j loop
    
    num_active_l_local = iii
    ! ***  iii will count the number of 
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, &
                       all_local_LA_array, 1, MPI_Int, &
                       DSIcomm, ierr)
    
    displacements_L_index(1) = 0 ! Want to start at L = 2
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
        displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo

End Subroutine Map_3D_Real_RK8
!---------------------------------------------------------------------------!  
! Subroutine: Map_3D_Int
!
!> @details Maps local values to 1 contiguous 1d array 
!
!> @param[in]  first_dim_size: 
!> @param[in]  second_dim_size: 
!> @param[in]  third_dim_size: 
!> @param[in]  Soln_Local: 
!> @param[in[  second_dim_size_1D: 
!> @param[in]  third_dim_size_1D:
!> @param[inout] Soln_Local_1D: 
!> @param[inout] Soln_Global: Global array containing entire spatial solution
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 1/8/2020
!---------------------------------------------------------------------------!
Subroutine Map_3D_Int(first_dim_size, second_dim_size, third_dim_size, Soln_Local, &
	       second_dim_size_1D, third_dim_size_1D, Soln_Local_1D, Map_Local_to_Global_1D)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use MPI
    
    implicit none
    
    ! *** Passed in variables

    integer, Intent(in)    :: first_dim_size  
    integer, Intent(in)    :: second_dim_size
    integer, Intent(in)    :: third_dim_size
    integer, Intent(in)    :: Soln_Local(first_dim_size, second_dim_size, third_dim_size) 
    integer, Intent(in)    :: second_dim_size_1D
    integer, Intent(in)    :: third_dim_size_1D
    integer, Intent(inout) :: Soln_Local_1D(LA_Local_no_ghost*second_dim_size_1D*third_dim_size_1D) 
    integer, Intent(inout) :: Map_Local_to_Global_1D(LA_Local_no_ghost*second_dim_size_1D*third_dim_size_1D) 
    
    ! *** Local variables
    integer :: i, j, k, kk, l, iii, i_gl, j_gl, l_gl, nn
    integer :: ierr
    integer, Allocatable, Dimension(:) :: all_local_LA_array !< Will contain the local LA value from each process
    
    if(.not.Allocated(all_local_LA_array) )then
        allocate(all_local_LA_array(active_domains))
    endif
    
    If(.not.allocated(displacements_L_index) )then
        allocate(displacements_L_index(active_domains))
    Endif
    
    all_local_LA_array(:) = 0
    
    ! *** Create 1D array of the solution excluding the ghost values
    iii = 0
    do j = 3, jc - 2
        ! *** Get first 2 rows
        do i = 3, ic - 2
        ! *** Grab the layers
            ! *** Map to the correct l index
            l  = lij(i,j)
            ! *** Get global L information
            i_gl = IL2IG(i)
            j_gl = JL2JG(j)
            l_gl = lij_global(i_gl,j_gl)
            do k = 1, second_dim_size_1D
                do nn = 1, third_dim_size_1D
                    ! *** Build up array containing active cells excluding ghost cells
                    if( l > 0 )then
						iii = iii + 1
                        Soln_Local_1D(iii) = Soln_Local(l,k,nn) !***Setup 1D aray of the solution
                        Map_Local_to_Global_1D(iii) = l_gl   ! *** Get corresponding global l value for a given local
                    endif
                enddo ! n loop
            enddo ! k loop
        enddo ! i loop
    enddo ! j loop
    
    num_active_l_local = iii
    ! ***  iii will count the number of 
    
    ! *** Collect the LA values from each process into an array so we can figure out the displacements
    !   for future collective operations
    call MPI_Allgather(num_active_l_local, 1, MPI_Int, &
                       all_local_LA_array, 1, MPI_Int, &
                       DSIcomm, ierr)
    
    displacements_L_index(1) = 0 ! Want to start at L = 2
    ! *** Build up the displacements array on each process
    do i = 2, active_domains
        displacements_L_index(i) = displacements_L_index(i-1) + all_local_LA_array(i-1)
    enddo


End Subroutine Map_3D_Int

End Module Mod_Map_Soln
