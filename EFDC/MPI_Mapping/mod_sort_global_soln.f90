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
    
Module Mod_Sort_Global_Soln

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    
    implicit none

    save
    
    Public :: Sort_Global_Soln
    
    interface Sort_Global_Soln
    
        Module Procedure Sort_1D_Real, &
                         Sort_1D_Real_RK8, &
                         Sort_1D_Int,  &
                         Sort_2D_Real, &
                         Sort_2D_Real_RK8, &
                         Sort_2D_Int,  &
                         Sort_3D_Real, &
                         Sort_3D_Real_RK8, &
                         Sort_3D_Int
    
    end interface

    Contains
    
!---------------------------------------------------------------------------!  
! Subroutine: Sort_1D_Real
!
!> @details Sorts the 1d solution array of Reals back to the global 1D (L) 
!
!> @param[in]  loop_bound_1d 
!> @param[in]  size_global_1d
!> @param[in]  Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in]  Mapping: Contains mapping of L values 
!> @param[out] Soln_Global_Remapped: Remapped global solution (ordered)
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!
Subroutine Sort_1D_Real(loop_bound_1d, size_global_1d, Soln_Global_1D, Mapping, Soln_Global_Remapped) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none
    
    ! *** Read in
    integer, Intent(in) :: loop_bound_1d
    integer, Intent(in) :: size_global_1d
    real(4), Intent(in)    :: Soln_Global_1D(size_global_1d)
    integer, Intent(in) :: Mapping(size_global_1d)
    real(4), Intent(inout) :: Soln_Global_Remapped(size_global_1d)
    
    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl

    ii = 0
    do l = 1, loop_bound_1d
        ii = ii + 1
        ! *** Get the global L
        l_gl = Mapping(l)
        ! *** Make sure we get valid mapping value
        if( l_gl > 0 )then
            ! *** Map back to the global (L,K) space
            Soln_Global_Remapped(l_gl) = Soln_Global_1D(ii)
        endif
    enddo
    
End subroutine Sort_1D_Real
!---------------------------------------------------------------------------!  
! Subroutine: Sort_1D_Real_RK8
!
!> @details Sorts the 1d solution array of Reals back to the global 1D (L) 
!
!> @param[in]  loop_bound_1d 
!> @param[in]  size_global_1d
!> @param[in]  Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in]  Mapping: Contains mapping of L values 
!> @param[out] Soln_Global_Remapped: Remapped global solution (ordered)
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!
Subroutine Sort_1D_Real_RK8(loop_bound_1d, size_global_1d, Soln_Global_1D, Mapping, Soln_Global_Remapped) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none
    
    ! *** Read in
    integer, Intent(in) :: loop_bound_1d
    integer, Intent(in) :: size_global_1d
    real(rkd), Intent(in)    :: Soln_Global_1D(size_global_1d)
    integer, Intent(in) :: Mapping(size_global_1d)
    real(rkd), Intent(inout) :: Soln_Global_Remapped(size_global_1d)
    
    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl

    ii = 0
    do l = 1, loop_bound_1d
        ii = ii + 1
        ! *** Get the global L
        l_gl = Mapping(l)
        ! *** Make sure we get valid mapping value
        if( l_gl > 0 )then
            ! *** Map back to the global (L,K) space
            Soln_Global_Remapped(l_gl) = Soln_Global_1D(ii)
        endif
    enddo
    
End subroutine Sort_1D_Real_RK8
!---------------------------------------------------------------------------!  
! Subroutine: Sort_1D_Int
!
!> @details Sorts the 1d solution array of Integers back to the global 1D (L) 
!
!> @param[in]  loop_bound_1d 
!> @param[in]  size_global_1d
!> @param[in]  Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in]  Mapping: Contains mapping of L values 
!> @param[out] Soln_Global_Remapped: Remapped global solution (ordered)
!
!> @author Zander Mausolff
!> @date 1/8/2019
!---------------------------------------------------------------------------!
Subroutine Sort_1D_Int(loop_bound_1d, size_global_1d, Soln_Global_1D, Mapping, Soln_Global_Remapped) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none
    
    ! *** Read in
    integer, Intent(in)    :: loop_bound_1d
    integer, Intent(in)    :: size_global_1d
    integer, Intent(in)    :: Soln_Global_1D(size_global_1d)
    integer, Intent(in)    :: Mapping(size_global_1d)
    integer, Intent(inout) :: Soln_Global_Remapped(size_global_1d)
    
    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl

    ii = 0
    do l = 1, loop_bound_1d
        ii = ii + 1
        ! *** Get the global L
        l_gl = Mapping(l)
        ! *** Make sure we get valid mapping value
        if( l_gl > 0  )then
            ! *** Map back to the global (L,K) space
            Soln_Global_Remapped(l_gl) = Soln_Global_1D(ii)
        endif
    enddo
    
End subroutine Sort_1D_Int
!---------------------------------------------------------------------------!  
! Subroutine: Sort_2D_Real(k_max, k_size, Soln_Local, Soln_Global)
!
!> @details Sorts the 1d solution array back to the global (L,K) 
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size
!> @param[in] Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in] Mapping: Contains mapping of L values 
!> @param[in] first_dim_size_gl
!> @param[in] second_dim_size_gl
!> @param[inout] Soln_Global_3D: Correctly ordered global array containing 
!
!> @author Zander Mausolff
!> @date 1/8/2020                   
!---------------------------------------------------------------------------!
Subroutine Sort_2D_Real(first_dim_size, second_dim_size, Soln_Global_1D, Mapping, &
                        first_dim_size_gl, second_dim_size_gl, Soln_Global_3D) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none
    
    ! *** Read in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    real(4), Intent(in)    :: Soln_Global_1D(first_dim_size*second_dim_size)
    integer, Intent(in) :: Mapping(first_dim_size*second_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl
    real(4), Intent(inout) :: Soln_Global_3D(first_dim_size_gl, second_dim_size_gl)

    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl
    
    ii = 0
    do l = 2, first_dim_size_gl
        do k = 1, second_dim_size_gl
            ii = ii + 1
            ! *** Get the global L
            l_gl = Mapping(ii)
            ! *** Make sure we get valid mapping value
            if(l_gl > 0 )then
                ! *** Map back to the global (L,K) space
                Soln_Global_3D(l_gl, k) = Soln_Global_1D(ii)
            endif
        enddo
    enddo
    
    
End subroutine Sort_2D_Real
!---------------------------------------------------------------------------!  
! Subroutine: Sort_2D_Real_RK8(k_max, k_size, Soln_Local, Soln_Global)
!
!> @details Sorts the 1d solution array back to the global (L,K) 
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size
!> @param[in] Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in] Mapping: Contains mapping of L values 
!> @param[in] first_dim_size_gl
!> @param[in] second_dim_size_gl
!> @param[inout] Soln_Global_3D: Correctly ordered global array containing 
!
!> @author Zander Mausolff
!> @date 1/8/2020                   
!---------------------------------------------------------------------------!
Subroutine Sort_2D_Real_RK8(first_dim_size, second_dim_size, Soln_Global_1D, Mapping, &
                        first_dim_size_gl, second_dim_size_gl, Soln_Global_3D) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none
    
    ! *** Read in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    real(rkd), Intent(in)    :: Soln_Global_1D(first_dim_size*second_dim_size)
    integer, Intent(in) :: Mapping(first_dim_size*second_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl
    real(rkd), Intent(inout) :: Soln_Global_3D(first_dim_size_gl, second_dim_size_gl)

    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl
    
    ii = 0
    do l = 2, first_dim_size_gl
        do k = 1, second_dim_size_gl
            ii = ii + 1
            ! *** Get the global L
            l_gl = Mapping(ii)
            ! *** Make sure we get valid mapping value
            if(l_gl > 0 )then
                ! *** Map back to the global (L,K) space
                Soln_Global_3D(l_gl, k) = Soln_Global_1D(ii)
            endif
        enddo
    enddo
    
    
End subroutine Sort_2D_Real_RK8
!---------------------------------------------------------------------------!  
! Subroutine: Sort_2D_Int(k_max, k_size, Soln_Local, Soln_Global)
!
!> @details Sorts the 1d solution array back to the global (L,K) 
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size
!> @param[in] Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in] Mapping: Contains mapping of L values 
!> @param[in] first_dim_size_gl
!> @param[in] second_dim_size_gl
!> @param[inout] Soln_Global_3D: Correctly ordered global array containing 
!
!> @author Zander Mausolff
!> @date 1/8/2020                  
!---------------------------------------------------------------------------!
Subroutine Sort_2D_Int(first_dim_size, second_dim_size, Soln_Global_1D, Mapping, &
                        first_dim_size_gl, second_dim_size_gl, Soln_Global_3D) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none
    
    ! *** Read in
    integer, Intent(in) :: first_dim_size
    integer, Intent(in) :: second_dim_size
    integer, Intent(in)    :: Soln_Global_1D(first_dim_size*second_dim_size)
    integer, Intent(in) :: Mapping(first_dim_size*second_dim_size)
    integer, Intent(in) :: first_dim_size_gl
    integer, Intent(in) :: second_dim_size_gl
    integer, Intent(inout) :: Soln_Global_3D(first_dim_size_gl, second_dim_size_gl)

    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl
    
    ii = 0
    do l = 2, first_dim_size_gl
        do k = 1, second_dim_size_gl
            ii = ii + 1
            ! *** Get the global L
            l_gl = Mapping(ii)
            ! *** Make sure we get valid mapping value
            if(l_gl > 0 )then
                ! *** Map back to the global (L,K) space
                Soln_Global_3D(l_gl, k) = Soln_Global_1D(ii)
            endif
        enddo
    enddo
    
    
End subroutine Sort_2D_Int
!---------------------------------------------------------------------------!  
! Subroutine: Sort_3D_Real
!
!> @details Sorts the 1d solution array back to the global (L,K) 
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size: 
!> @param[in] third_dim_size:
!> @param[in] Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in] Mapping: Contains mapping of L values 
!> @param[in] first_dim_size_gl: 
!> @param[in] second_dim_size_gl: 
!> @param[in] third_dim_size_gl:
!> @param[out] Soln_Global_3D: Correctly ordered global array containing 
!
!> @author Zander Mausolff
!> @date 1/8/2020                       
!---------------------------------------------------------------------------!
Subroutine Sort_3D_Real(first_dim_size, second_dim_size, third_dim_size, Soln_Global_1D, Mapping, &
            first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_3D) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none

    
    ! *** Read in
    integer, Intent(in) :: first_dim_size
	Integer, Intent(in) :: second_dim_size
    integer, Intent(in) :: third_dim_size
    real(4), Intent(in)    :: Soln_Global_1D(first_dim_size*second_dim_size*third_dim_size)
    integer, Intent(in) :: Mapping(first_dim_size*second_dim_size*third_dim_size)
    integer, Intent(in) :: first_dim_size_gl
	Integer, Intent(in) :: second_dim_size_gl
    integer, Intent(in) :: third_dim_size_gl
    real(4), Intent(inout) :: Soln_Global_3D(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl)
    
    
    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl, nn
    
    ii = 0
    do l = 2, first_dim_size_gl
        do k = 1, second_dim_size_gl
            do nn = 1, third_dim_size_gl
				ii = ii + 1
				!***Get the global L
				l_gl = Mapping(ii)
                ! *** Make sure we get valid mapping value
                if(l_gl > 0 )then
                    ! *** Map back to the global (L,K) space
                    Soln_Global_3D(l_gl, k, nn) = Soln_Global_1D(ii)
                endif
            enddo ! nn loop
        enddo ! k loop
    enddo ! l loop
    
    
End subroutine Sort_3D_Real
!---------------------------------------------------------------------------!  
! Subroutine: Sort_3D_Real_RK8
!
!> @details Sorts the 1d solution array back to the global (L,K) 
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size: 
!> @param[in] third_dim_size:
!> @param[in] Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in] Mapping: Contains mapping of L values 
!> @param[in] first_dim_size_gl: 
!> @param[in] second_dim_size_gl: 
!> @param[in] third_dim_size_gl:
!> @param[out] Soln_Global_3D: Correctly ordered global array containing 
!
!> @author Zander Mausolff
!> @date 1/8/2020                       
!---------------------------------------------------------------------------!
Subroutine Sort_3D_Real_RK8(first_dim_size, second_dim_size, third_dim_size, Soln_Global_1D, Mapping, &
            first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_3D) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none

    
    ! *** Read in
    integer, Intent(in) :: first_dim_size
	Integer, Intent(in) :: second_dim_size
    integer, Intent(in) :: third_dim_size
    real(rkd), Intent(in)    :: Soln_Global_1D(first_dim_size*second_dim_size*third_dim_size)
    integer, Intent(in) :: Mapping(first_dim_size*second_dim_size*third_dim_size)
    integer, Intent(in) :: first_dim_size_gl
	Integer, Intent(in) :: second_dim_size_gl
    integer, Intent(in) :: third_dim_size_gl
    real(rkd), Intent(inout) :: Soln_Global_3D(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl)
    
    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl, nn
    
    ii = 0
    do l = 2, first_dim_size_gl
        do k = 1, second_dim_size_gl
            do nn = 1, third_dim_size_gl
				ii = ii + 1
				!***Get the global L
				l_gl = Mapping(ii)
                ! *** Make sure we get valid mapping value
                if(l_gl > 0 )then
                    ! *** Map back to the global (L,K) space
                    Soln_Global_3D(l_gl, k, nn) = Soln_Global_1D(ii)
                endif
            enddo ! nn loop
        enddo ! k loop
    enddo ! l loop
    
    
End subroutine Sort_3D_Real_RK8
!---------------------------------------------------------------------------!  
! Subroutine: Sort_3D_Int
!
!> @details Sorts the 1d solution array back to the global (L,K) 
!
!> @param[in] first_dim_size: 
!> @param[in] second_dim_size: 
!> @param[in] third_dim_size:
!> @param[in] Soln_Global_1D: 1D array of the global values (unordered)
!> @param[in] Mapping: Contains mapping of L values 
!> @param[in] first_dim_size_gl: 
!> @param[in] second_dim_size_gl: 
!> @param[in] third_dim_size_gl:
!> @param[out] Soln_Global_3D: Correctly ordered global array containing 
!
!> @author Zander Mausolff
!> @date 1/8/2020                       
!---------------------------------------------------------------------------!
Subroutine Sort_3D_Int(first_dim_size, second_dim_size, third_dim_size, Soln_Global_1D, Mapping, &
            first_dim_size_gl, second_dim_size_gl, third_dim_size_gl, Soln_Global_3D) 

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Write_Out
    use Variables_MPI_Mapping

    implicit none

    
    ! *** Read in
    integer, Intent(in) :: first_dim_size
	Integer, Intent(in) :: second_dim_size
    integer, Intent(in) :: third_dim_size
    integer, Intent(in) :: Soln_Global_1D(first_dim_size*second_dim_size*third_dim_size)
    integer, Intent(in) :: Mapping(first_dim_size*second_dim_size*third_dim_size)
    integer, Intent(in) :: first_dim_size_gl
	Integer, Intent(in) :: second_dim_size_gl
    integer, Intent(in) :: third_dim_size_gl
    integer, Intent(inout) :: Soln_Global_3D(first_dim_size_gl, second_dim_size_gl, third_dim_size_gl)
    
    
    ! *** Local variables
    integer :: i, j ,k, l, ii, l_gl, nn
    
    ii = 0
    do l = 2, first_dim_size_gl
        do k = 1, second_dim_size_gl
            do nn = 1, third_dim_size_gl
				ii = ii + 1
				!***Get the global L
				l_gl = Mapping(ii)
                ! *** Make sure we get valid mapping value
                if(l_gl > 0 )then
                    ! *** Map back to the global (L,K) space
                    Soln_Global_3D(l_gl, k, nn) = Soln_Global_1D(ii)
                endif
            enddo ! nn loop
        enddo ! k loop
    enddo ! l loop
    
    
End subroutine Sort_3D_Int

End Module Mod_Sort_Global_Soln