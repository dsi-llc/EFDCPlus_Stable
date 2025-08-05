! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------

Module Mod_Map_Gather_Sort

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_MapGatherSort
    
    implicit none

    save
    
    Public :: Map_Gather_Sort
    
    interface Map_Gather_Sort
    
        Module Procedure Map_Gather_Sort_Real_1D, &
                         Map_Gather_Sort_Real_2D, &
                         Map_Gather_Sort_Real_3D, &   
                         Map_Gather_Sort_Real_1D_RK8, &
                         Map_Gather_Sort_Real_2D_RK8, &
                         Map_Gather_Sort_Real_3D_RK8, &
                         Map_Gather_Sort_Int_1D, &
                         Map_Gather_Sort_Int_2D, &
                         Map_Gather_Sort_Int_3D
    
    end interface

    Contains
    
    
!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Real_1D
!
!> @details 
!> @param[in]    size_local_1d
!> @param[in]    Initial_Soln_Local
!> @param[in]    size_global_1d
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!    
Subroutine Map_Gather_Sort_Real_1D(size_local_1d, Initial_Soln_Local, size_global_1d, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in) :: size_local_1d
    real(4),    Intent(in) :: Initial_Soln_Local(size_local_1d)
    integer, Intent(in) :: size_global_1d
    real(4), Intent(inout) :: Final_Global_Soln(size_global_1d)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Soln_Local_1D(LA_Local_no_ghost))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost))
    allocate(Soln_Global_1D(size_global_1d))
    allocate(Global_Local_L_to_Global(size_global_1d))
    
    Soln_Local_1D  = 0.0 
    Soln_Global_1D = 0.0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(size_local_1d, Initial_Soln_Local, la_local_no_ghost, Soln_Local_1D, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(la_local_no_ghost, Soln_Local_1D, size_global_1d, &
                    Soln_Global_1D, Map_Local_L_to_Global, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(size_global_1d, size_global_1d, Soln_Global_1D, Global_Local_L_to_Global, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Soln_Local_1D)   
    deallocate(Soln_Global_1D) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
    
End Subroutine Map_Gather_Sort_Real_1D

!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Real_2D
!
!> @details 
!> @param[in]  first_dim_size
!> @param[in]  second_dim_size
!> @param[in]  Initial_Soln_Local
!> @param[in]  first_dim_size_gl
!> @param[in]  second_dim_size_gl
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!
Subroutine Map_Gather_Sort_Real_2D(first_dim_size,    second_dim_size,    Initial_Soln_Local, &
                                   first_dim_size_gl, second_dim_size_gl, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in)    :: first_dim_size
    integer, Intent(in)    :: second_dim_size
    real(4), Intent(in)    :: Initial_Soln_Local(first_dim_size, second_dim_size)
    integer, Intent(in)    :: first_dim_size_gl
    integer, Intent(in)    :: second_dim_size_gl
    real(4), Intent(inout) :: Final_Global_Soln(first_dim_size_gl, second_dim_size_gl)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Soln_Local_1D(LA_Local_no_ghost*second_dim_size))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost*second_dim_size))
    
    allocate(Soln_Global_1D(first_dim_size_gl*second_dim_size_gl))
    allocate(Global_Local_L_to_Global(first_dim_size_gl*second_dim_size_gl))
    
    Soln_Local_1D  = 0.0 
    Soln_Global_1D = 0.0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(first_dim_size, second_dim_size, Initial_Soln_Local, &
                  second_dim_size_gl, Soln_Local_1D, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(la_local_no_ghost, second_dim_size,    Soln_Local_1D, &
                     first_dim_size_gl, second_dim_size_gl, Soln_Global_1D, &
                     Map_Local_L_to_Global, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(first_dim_size_gl, second_dim_size_gl, Soln_Global_1D, Global_Local_L_to_Global, &
                            first_dim_size_gl, second_dim_size_gl, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Soln_Local_1D)   
    deallocate(Soln_Global_1D) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
End Subroutine Map_Gather_Sort_Real_2D

!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Real_3D
!
!> @details 
!> @param[in]    size_local_1d
!> @param[in]    Initial_Soln_Local
!> @param[in]    size_global_1d
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!    
Subroutine Map_Gather_Sort_Real_3D(first_dim_size, second_dim_size, third_dim_size, Initial_Soln_Local, &
                                   first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    integer, Intent(in) :: third_dim_size
    real(4),    Intent(in) :: Initial_Soln_Local(first_dim_size, second_dim_size, third_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl 
    integer, Intent(in) :: third_dim_size_gl
    real(4), Intent(inout) :: Final_Global_Soln(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Soln_Local_1D(LA_Local_no_ghost*second_dim_size*third_dim_size))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost*second_dim_size*third_dim_size))
    
    allocate(Soln_Global_1D(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl))
    allocate(Global_Local_L_to_Global(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl))
    
    Soln_Local_1D  = 0.0 
    Soln_Global_1D = 0.0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(first_dim_size, second_dim_size, third_dim_size, Initial_Soln_Local, &
                  second_dim_size_gl,third_dim_size_gl, Soln_Local_1D, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(first_dim_size, second_dim_size, third_dim_size, Soln_Local_1D, Map_Local_L_to_Global, &
                     first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, &
                     Soln_Global_1D, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_1D, Global_Local_L_to_Global, &
                            first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Soln_Local_1D)   
    deallocate(Soln_Global_1D) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
    
End Subroutine Map_Gather_Sort_Real_3D
!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Real_1D_RK8
!
!> @details 
!> @param[in]    size_local_1d
!> @param[in]    Initial_Soln_Local
!> @param[in]    size_global_1d
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!    
Subroutine Map_Gather_Sort_Real_1D_RK8(size_local_1d, Initial_Soln_Local, size_global_1d, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in) :: size_local_1d
    real(rkd),    Intent(in) :: Initial_Soln_Local(size_local_1d)
    integer, Intent(in) :: size_global_1d
    real(rkd), Intent(inout) :: Final_Global_Soln(size_global_1d)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Soln_Local_1D_RK8(LA_Local_no_ghost))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost))
    allocate(Soln_Global_1D_RK8(size_global_1d))
    allocate(Global_Local_L_to_Global(size_global_1d))
    
    Soln_Local_1D_RK8  = 0.0 
    Soln_Global_1D_RK8 = 0.0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(size_local_1d, Initial_Soln_Local, la_local_no_ghost, Soln_Local_1D_RK8, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(la_local_no_ghost, Soln_Local_1D_RK8, size_global_1d, &
                    Soln_Global_1D_RK8, Map_Local_L_to_Global, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(size_global_1d, size_global_1d, Soln_Global_1D_RK8, Global_Local_L_to_Global, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Soln_Local_1D_RK8)   
    deallocate(Soln_Global_1D_RK8) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
    
End Subroutine Map_Gather_Sort_Real_1D_RK8

!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Real_2D_RK8
!
!> @details 
!> @param[in]  first_dim_size
!> @param[in]  second_dim_size
!> @param[in]  Initial_Soln_Local
!> @param[in]  first_dim_size_gl
!> @param[in]  second_dim_size_gl
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!
Subroutine Map_Gather_Sort_Real_2D_RK8(first_dim_size, second_dim_size, Initial_Soln_Local, &
                        first_dim_size_gl, second_dim_size_gl, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    real(rkd),    Intent(in) :: Initial_Soln_Local(first_dim_size, second_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl
    real(rkd), Intent(inout) :: Final_Global_Soln(first_dim_size_gl, second_dim_size_gl)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Soln_Local_1D_RK8(LA_Local_no_ghost*second_dim_size))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost*second_dim_size))
    allocate(Soln_Global_1D_RK8(first_dim_size_gl*second_dim_size_gl))
    allocate(Global_Local_L_to_Global(first_dim_size_gl*second_dim_size_gl))
    
    Soln_Local_1D_RK8  = 0.0 
    Soln_Global_1D_RK8 = 0.0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(first_dim_size, second_dim_size, Initial_Soln_Local, &
                  second_dim_size_gl, Soln_Local_1D_RK8, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(la_local_no_ghost, second_dim_size, Soln_Local_1D_RK8, &
                     first_dim_size_gl, second_dim_size_gl, &
                     Soln_Global_1D_RK8, Map_Local_L_to_Global, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(first_dim_size_gl, second_dim_size_gl, Soln_Global_1D_RK8, Global_Local_L_to_Global, &
                            first_dim_size_gl, second_dim_size_gl, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Soln_Local_1D_RK8)   
    deallocate(Soln_Global_1D_RK8) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
End Subroutine Map_Gather_Sort_Real_2D_RK8

!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Real_3D_RK8
!
!> @details 
!> @param[in]    size_local_1d
!> @param[in]    Initial_Soln_Local
!> @param[in]    size_global_1d
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!    
Subroutine Map_Gather_Sort_Real_3D_RK8(first_dim_size, second_dim_size, third_dim_size, Initial_Soln_Local, &
                                   first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    integer, Intent(in) :: third_dim_size
    real(rkd),    Intent(in) :: Initial_Soln_Local(first_dim_size, second_dim_size, third_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl 
    integer, Intent(in) :: third_dim_size_gl
    real(rkd), Intent(inout) :: Final_Global_Soln(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Soln_Local_1D_RK8(LA_Local_no_ghost*second_dim_size*third_dim_size))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost*second_dim_size*third_dim_size))
    
    allocate(Soln_Global_1D_RK8(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl))
    allocate(Global_Local_L_to_Global(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl))
    
    Soln_Local_1D_RK8  = 0.0 
    Soln_Global_1D_RK8 = 0.0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(first_dim_size, second_dim_size, third_dim_size, Initial_Soln_Local, &
                  second_dim_size_gl,third_dim_size_gl, Soln_Local_1D_RK8, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(first_dim_size, second_dim_size, third_dim_size, Soln_Local_1D_RK8, Map_Local_L_to_Global, &
                     first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, &
                     Soln_Global_1D_RK8, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_1D_RK8, Global_Local_L_to_Global, &
                            first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Soln_Local_1D_RK8)   
    deallocate(Soln_Global_1D_RK8) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
    
End Subroutine Map_Gather_Sort_Real_3D_RK8

!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Int_1D
!
!> @details 
!> @param[in]    size_local_1d
!> @param[in]    Initial_Soln_Local
!> @param[in]    size_global_1d
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!    
Subroutine Map_Gather_Sort_Int_1D(size_local_1d, Initial_Soln_Local, size_global_1d, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in)    :: size_local_1d
    integer, Intent(in)    :: Initial_Soln_Local(size_local_1d)
    integer, Intent(in)    :: size_global_1d
    integer, Intent(inout) :: Final_Global_Soln(size_global_1d)
    
    ! *** Local variables

    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Int_Soln_Local_1D(LA_Local_no_ghost))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost))
    
    allocate(Int_Soln_Global_1D(size_global_1d))
    allocate(Global_Local_L_to_Global(size_global_1d))
    
    Int_Soln_Local_1D  = 0
    Int_Soln_Global_1D = 0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(size_local_1d, Initial_Soln_Local, la_local_no_ghost, Int_Soln_Local_1D, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(la_local_no_ghost, Int_Soln_Local_1D, size_global_1d, &
                    Int_Soln_Global_1D, Map_Local_L_to_Global, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(size_global_1d, size_global_1d, Int_Soln_Global_1D, Global_Local_L_to_Global, Final_Global_Soln)
    Endif
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Int_Soln_Local_1D)   
    deallocate(Int_Soln_Global_1D) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
    
End Subroutine Map_Gather_Sort_Int_1D
!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Int_2D
!
!> @details 
!> @param[in]  first_dim_size
!> @param[in]  second_dim_size
!> @param[in]  Initial_Soln_Local
!> @param[in]  first_dim_size_gl
!> @param[in]  second_dim_size_gl
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!
Subroutine Map_Gather_Sort_Int_2D(first_dim_size,    second_dim_size,    Initial_Soln_Local, &
                                  first_dim_size_gl, second_dim_size_gl, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in)    :: first_dim_size
    integer, Intent(in)    :: second_dim_size
    integer, Intent(in)    :: Initial_Soln_Local(first_dim_size, second_dim_size)
    integer, Intent(in)    :: first_dim_size_gl
    integer, Intent(in)    :: second_dim_size_gl
    integer, Intent(inout) :: Final_Global_Soln(first_dim_size_gl, second_dim_size_gl)
    
    ! *** Local variables
    allocate(Int_Soln_Local_1D(LA_Local_no_ghost*second_dim_size))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost*second_dim_size))
    
    allocate(Int_Soln_Global_1D(first_dim_size_gl*second_dim_size_gl))
    allocate(Global_Local_L_to_Global(first_dim_size_gl*second_dim_size_gl))
    
    
    Int_Soln_Local_1D  = 0
    Int_Soln_Global_1D = 0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0
    
    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(first_dim_size,     second_dim_size,   Initial_Soln_Local, &
                  second_dim_size_gl, Int_Soln_Local_1D, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(la_local_no_ghost,  second_dim_size,    Int_Soln_Local_1D,  &
                     first_dim_size_gl,  second_dim_size_gl, Int_Soln_Global_1D, &
                     Map_Local_L_to_Global, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(first_dim_size_gl, second_dim_size_gl, Int_Soln_Global_1D, Global_Local_L_to_Global, &
                            first_dim_size_gl, second_dim_size_gl, Final_Global_Soln)
    Endif
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Int_Soln_Local_1D)   
    deallocate(Int_Soln_Global_1D) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
End Subroutine Map_Gather_Sort_Int_2D
!---------------------------------------------------------------------------!  
! Subroutine: Map_Gather_Sort_Int_3D
!
!> @details 
!> @param[in]    size_local_1d
!> @param[in]    Initial_Soln_Local
!> @param[in]    size_global_1d
!> @param[inout] Final_Global_Soln
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!    
Subroutine Map_Gather_Sort_Int_3D(first_dim_size,    second_dim_size,    third_dim_size,    Initial_Soln_Local, &
                                  first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Final_Global_Soln)

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    
    ! ***
    use Mod_Map_Soln
    use Mod_Gather_Soln
    use Mod_Sort_Global_Soln
    ! *** 
    
    implicit none
    
    ! *** Passed in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    integer, Intent(in) :: third_dim_size
    integer, Intent(in) :: Initial_Soln_Local(first_dim_size, second_dim_size, third_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl 
    integer, Intent(in) :: third_dim_size_gl
    integer, Intent(inout) :: Final_Global_Soln(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl)
    
    ! *** Local variables
    
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    ! *** Allocate the temporary arrays used for the collapsing and remapping
    allocate(Int_Soln_Local_1D(LA_Local_no_ghost*second_dim_size*third_dim_size))
    allocate(Map_Local_L_to_Global(LA_Local_no_ghost*second_dim_size*third_dim_size))
    
    allocate(Int_Soln_Global_1D(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl))
    allocate(Global_Local_L_to_Global(first_dim_size_gl*second_dim_size_gl*third_dim_size_gl))
    
    Int_Soln_Local_1D  = 0 
    Int_Soln_Global_1D = 0
    Map_Local_L_to_Global = 0
    Global_Local_L_to_Global = 0

    ! *** Collapse the solution to 1D, excluding the ghost cells and get the local to global mapping
    call Map_Soln(first_dim_size, second_dim_size, third_dim_size, Initial_Soln_Local, &
                  second_dim_size_gl,third_dim_size_gl, Int_Soln_Local_1D, Map_Local_L_to_Global)
    
    ! *** Gather local solutions and mappings onto a single array 
    call Gather_Soln(first_dim_size, second_dim_size, third_dim_size, Int_Soln_Local_1D, Map_Local_L_to_Global, &
                     first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, &
                     Int_Soln_Global_1D, Global_Local_L_to_Global)
    
    If( process_id == master_id  )then 
      ! *** Sort the global array so L indexing is correct
      call Sort_Global_Soln(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Int_Soln_Global_1D, Global_Local_L_to_Global, &
                            first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Final_Global_Soln)
    Endif
        
    
    ! *** Deallocate so the other routines can set their own sizes
    deallocate(Int_Soln_Local_1D)   
    deallocate(Int_Soln_Global_1D) 
    deallocate(Map_Local_L_to_Global)
    deallocate(Global_Local_L_to_Global)
    
    
End Subroutine Map_Gather_Sort_Int_3D
End Module Mod_Map_Gather_Sort
