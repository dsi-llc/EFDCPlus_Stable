! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!  
! Module: Variables_MPI_MapGatherSort
!
!> @details 
!
!> @author Zander Mausolff
!> @date 1/9/2019
!---------------------------------------------------------------------------! 
Module Variables_MPI_MapGatherSort

    use GLOBAL

    implicit none
    
    real(4), Allocatable, Dimension(:)    :: Soln_Local_1D !< Stores collapsed to 1D arrays of the local soln 
    real(4), Allocatable, Dimension(:)    :: Soln_Global_1D !< Stores collapsed 1D arrays containing the global soln (unsorted)
    integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global !< Stores collapsed mapping of local L to global L 
    integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global !< Stores (from all processes) collapsed mapping of local L to global L

    integer, Allocatable, Dimension(:)    :: Int_Soln_Local_1D !< Stores collapsed to 1D arrays of the local soln 
    integer, Allocatable, Dimension(:)    :: Int_Soln_Global_1D !< Stores collapsed 1D arrays containing the global soln (unsorted)
 
    real(rkd), Allocatable, Dimension(:)    :: Soln_Local_1D_RK8 !< Stores collapsed to 1D arrays of the local soln 
    real(rkd), Allocatable, Dimension(:)    :: Soln_Global_1D_RK8 !< Stores collapsed 1D arrays containing the global soln (unsorted)

    ! *** Real types
    Type pointer_soln_1D_Real !< Points to 1D arrays that are going to be written out
        real(4), Pointer :: val(:) => null()
    End type
    
    Type pointer_soln_2D_Real !< Points to 2D arrays that are going to be written out
        real(4), Pointer :: val(:,:) => null()
    End type
    
    Type pointer_soln_3D_Real !< Points to 2D arrays that are going to be written out
        real(4), Pointer :: val(:,:,:) => null()
    End type
    
    ! *** Real(8) types
    Type pointer_soln_1D_Real_RK8 !< Points to 1D arrays that are going to be written out
        real(rkd), Pointer :: val(:) => null()
    End type
    
    Type pointer_soln_2D_Real_RK8 !< Points to 2D arrays that are going to be written out
        real(rkd), Pointer :: val(:,:) => null()
    End type
    
    Type pointer_soln_3D_Real_RK8 !< Points to 2D arrays that are going to be written out
        real(rkd), Pointer :: val(:,:,:) => null()
    End type
    
    ! *** Integer types
    Type pointer_soln_1D_Int !< Points to 1D arrays that are going to be written out
        integer, Pointer :: val(:) => null()
    End type
    
    Type pointer_soln_2D_Int !< Points to 2D arrays that are going to be written out
        integer, Pointer :: val(:,:) => null()
    End type
    
    Type pointer_soln_3D_Int !< Points to 2D arrays that are going to be written out
        integer, Pointer :: val(:,:,:) => null()
    End type
    
    ! *** Real(4)
    type(pointer_soln_1D_Real), Dimension(100) :: Local_Arrays_to_Write_1D_Real
    type(pointer_soln_1D_Real), Dimension(100) :: Global_Arrays_to_Write_1D_Real
    
    type(pointer_soln_2D_Real), Dimension(100) :: Local_Arrays_to_Write_2D_Real
    type(pointer_soln_2D_Real), Dimension(100) :: Global_Arrays_to_Write_2D_Real
    
    type(pointer_soln_3D_Real), Dimension(100) :: Local_Arrays_to_Write_3D_Real
    type(pointer_soln_3D_Real), Dimension(100) :: Global_Arrays_to_Write_3D_Real
    
    ! *** Real(8)
    type(pointer_soln_1D_Real_RK8), Dimension(100) :: Local_Arrays_to_Write_1D_Real_RK8
    type(pointer_soln_1D_Real_RK8), Dimension(100) :: Global_Arrays_to_Write_1D_Real_RK8
    
    type(pointer_soln_2D_Real_RK8), Dimension(100) :: Local_Arrays_to_Write_2D_Real_RK8
    type(pointer_soln_2D_Real_RK8), Dimension(100) :: Global_Arrays_to_Write_2D_Real_RK8
    
    type(pointer_soln_3D_Real_RK8), Dimension(100) :: Local_Arrays_to_Write_3D_Real_RK8
    type(pointer_soln_3D_Real_RK8), Dimension(100) :: Global_Arrays_to_Write_3D_Real_RK8
    
    ! *** Integer
    type(pointer_soln_1D_Int), Dimension(100) :: Local_Arrays_to_Write_1D_Int
    type(pointer_soln_1D_Int), Dimension(100) :: Global_Arrays_to_Write_1D_Int
    
    type(pointer_soln_2D_Int), Dimension(100) :: Local_Arrays_to_Write_2D_Int
    type(pointer_soln_2D_Int), Dimension(100) :: Global_Arrays_to_Write_2D_Int
    
    type(pointer_soln_3D_Int), Dimension(100) :: Local_Arrays_to_Write_3D_Int
    type(pointer_soln_3D_Int), Dimension(100) :: Global_Arrays_to_Write_3D_Int
    
    integer, Dimension(100) :: Dim_Array_Written_Out

End Module Variables_MPI_MapGatherSort