! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!  
! Module: Mod_Assign_Loc_Glob_For_Write
!
!> @details Handles assigning pointer arrays for local and global solutions
!! so they can be processed by the map/gather/sort routines in prep for writing
!! out to EE_*.OUT binary files.  It handles both real and integer arrays 
!! automaticall through the module interface.
!
!> @param[in]  index
!> @param[in]  first_dim_local
!> @param[in]  Soln_Local(first_dim_local)
!> @param[in]  first_dim_gl
!> @param[in]  Soln_Global(first_dim_gl)
!
!> @author Zander Mausolff
!> @date 1/9/2019
!---------------------------------------------------------------------------!
Module Mod_Assign_Loc_Glob_For_Write

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_MapGatherSort
    
    implicit none

    save
    
    Public :: Assign_Loc_Glob_For_Write
    
    interface Assign_Loc_Glob_For_Write
    
        Module Procedure Assign_Loc_Glob_For_Write_1D_Real, &
                         Assign_Loc_Glob_For_Write_1D_Real_RK8, &
                         Assign_Loc_Glob_For_Write_1D_Int,  &
                         Assign_Loc_Glob_For_Write_2D_Real, &
                         Assign_Loc_Glob_For_Write_2D_Real_RK8, &
                         Assign_Loc_Glob_For_Write_2D_Int,  &                         
                         Assign_Loc_Glob_For_Write_3D_Real, &
                         Assign_Loc_Glob_For_Write_3D_Real_RK8, &                    
                         Assign_Loc_Glob_For_Write_3D_Int
    
        
    end interface

    Contains
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_1D_Real
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  Soln_Local(first_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_1D_Real(index, first_dim_local, Soln_Local, &
                                                 first_dim_gl, Soln_Global)
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in) :: index
        integer, Intent(in) :: first_dim_local
        real(4), Target, Intent(inout) :: Soln_Local(first_dim_local)
        integer, Intent(in) :: first_dim_gl
        real(4),  Target,   Intent(inout) :: Soln_Global(first_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_1D_Real(index).val  => Soln_Local
        Global_Arrays_to_Write_1D_Real(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 1
        
    End Subroutine Assign_Loc_Glob_For_Write_1D_Real
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_1D_Real
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  Soln_Local(first_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_1D_Real_RK8(index, first_dim_local, Soln_Local, &
                                                 first_dim_gl, Soln_Global)
        use GLOBAL
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in) :: index
        integer, Intent(in) :: first_dim_local
        real(rkd), Target, Intent(inout) :: Soln_Local(first_dim_local)
        integer, Intent(in) :: first_dim_gl
        real(rkd),  Target,   Intent(inout) :: Soln_Global(first_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_1D_Real_RK8(index).val  => Soln_Local
        Global_Arrays_to_Write_1D_Real_RK8(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 100 + rkd
        
    End Subroutine Assign_Loc_Glob_For_Write_1D_Real_RK8
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_1D_Int
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  Soln_Local(first_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_1D_Int(index, first_dim_local, Soln_Local, &
                                                 first_dim_gl, Soln_Global)
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in) :: index
        integer, Intent(in) :: first_dim_local
        integer,Target, Intent(inout) :: Soln_Local(first_dim_local)
        integer, Intent(in) :: first_dim_gl
        integer,Target,  Intent(inout) :: Soln_Global(first_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_1D_Int(index).val  => Soln_Local
        Global_Arrays_to_Write_1D_Int(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 10
        
    End Subroutine Assign_Loc_Glob_For_Write_1D_Int
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_2D_Real
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  second_dim_local
    !> @param[in]  Soln_Local(first_dim_local, second_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  second_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl, second_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_2D_Real(index, first_dim_local, second_dim_local, Soln_Local, &
                                                 first_dim_gl, second_dim_gl, Soln_Global)
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in) :: index
        integer, Intent(in) :: first_dim_local
        integer, Intent(in) :: second_dim_local
        real(4), Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local)
        integer, Intent(in) :: first_dim_gl
        integer, Intent(in) :: second_dim_gl
        real(4), Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_2D_Real(index).val  => Soln_Local
        Global_Arrays_to_Write_2D_Real(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 2
    
    End Subroutine Assign_Loc_Glob_For_Write_2D_Real
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_2D_Real
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  second_dim_local
    !> @param[in]  Soln_Local(first_dim_local, second_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  second_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl, second_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_2D_Real_RK8(index, first_dim_local, second_dim_local, Soln_Local, &
                                                 first_dim_gl, second_dim_gl, Soln_Global)
        use GLOBAL
                                        
        implicit none
        
        ! *** Read in 
        integer, Intent(in) :: index
        integer, Intent(in) :: first_dim_local
        integer, Intent(in) :: second_dim_local
        real(rkd), Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local)
        integer, Intent(in) :: first_dim_gl
        integer, Intent(in) :: second_dim_gl
        real(rkd), Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_2D_Real_RK8(index).val  => Soln_Local
        Global_Arrays_to_Write_2D_Real_RK8(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 200 + rkd
    
    End Subroutine Assign_Loc_Glob_For_Write_2D_Real_RK8
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_2D_Int
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  second_dim_local
    !> @param[in]  Soln_Local(first_dim_local, second_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  second_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl, second_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_2D_Int(index, first_dim_local, second_dim_local, Soln_Local, &
                                                 first_dim_gl, second_dim_gl, Soln_Global)
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in) :: index
        integer, Intent(in) :: first_dim_local
        integer, Intent(in) :: second_dim_local
        integer, Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local)
        integer, Intent(in) :: first_dim_gl
        integer, Intent(in) :: second_dim_gl
        integer, Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_2D_Int(index).val  => Soln_Local
        Global_Arrays_to_Write_2D_Int(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 20
    
    End Subroutine Assign_Loc_Glob_For_Write_2D_Int
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_3D_Real
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  second_dim_local
    !> @param[in]  third_dim_local
    !> @param[in]  Soln_Local(first_dim_local, second_dim_local, third_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  second_dim_gl
    !> @param[in]  third_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_3D_Real(index, first_dim_local, second_dim_local, third_dim_local, Soln_Local, &
                                                 first_dim_gl, second_dim_gl, third_dim_gl, Soln_Global)
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in)         :: index
        integer, Intent(in)         :: first_dim_local
        integer, Intent(in)         :: second_dim_local
        integer, Intent(in)         :: third_dim_local
        real(4), Target, Intent(inout) :: Soln_Local(first_dim_local, second_dim_local, third_dim_local)
        integer, Intent(in)         :: first_dim_gl
        integer, Intent(in)         :: second_dim_gl
        integer, Intent(in)         :: third_dim_gl
        real(4), Target, Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_3D_Real(index).val  => Soln_Local
        Global_Arrays_to_Write_3D_Real(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 3
    
    End Subroutine Assign_Loc_Glob_For_Write_3D_Real
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_3D_Real
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  second_dim_local
    !> @param[in]  third_dim_local
    !> @param[in]  Soln_Local(first_dim_local, second_dim_local, third_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  second_dim_gl
    !> @param[in]  third_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_3D_Real_RK8(index, first_dim_local, second_dim_local, third_dim_local, Soln_Local, &
                                                 first_dim_gl, second_dim_gl, third_dim_gl, Soln_Global)
        use GLOBAL
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in)         :: index
        integer, Intent(in)         :: first_dim_local
        integer, Intent(in)         :: second_dim_local
        integer, Intent(in)         :: third_dim_local
        real(rkd), Target, Intent(inout) :: Soln_Local(first_dim_local, second_dim_local, third_dim_local)
        integer, Intent(in)         :: first_dim_gl
        integer, Intent(in)         :: second_dim_gl
        integer, Intent(in)         :: third_dim_gl
        real(rkd), Target, Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_3D_Real_RK8(index).val  => Soln_Local
        Global_Arrays_to_Write_3D_Real_RK8(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 300 + rkd
    
    End Subroutine Assign_Loc_Glob_For_Write_3D_Real_RK8
    !---------------------------------------------------------------------------!  
    ! Subroutine: Assign_Loc_Glob_For_Write_3D_Int
    !
    !> @details 
    !
    !> @param[in]  index
    !> @param[in]  first_dim_local
    !> @param[in]  second_dim_local
    !> @param[in]  third_dim_local
    !> @param[in]  Soln_Local(first_dim_local, second_dim_local, third_dim_local)
    !> @param[in]  first_dim_gl
    !> @param[in]  second_dim_gl
    !> @param[in]  third_dim_gl
    !> @param[in]  Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
    !
    !> @author Zander Mausolff
    !> @date 1/9/2019
    !---------------------------------------------------------------------------!
    Subroutine Assign_Loc_Glob_For_Write_3D_Int(index, first_dim_local, second_dim_local, third_dim_local, Soln_Local, &
                                                 first_dim_gl, second_dim_gl, third_dim_gl, Soln_Global)
                                                 
        implicit none
        
        ! *** Read in 
        integer, Intent(in)            :: index
        integer, Intent(in)            :: first_dim_local
        integer, Intent(in)            :: second_dim_local
        integer, Intent(in)            :: third_dim_local
        integer,Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local, third_dim_local)
        integer, Intent(in)            :: first_dim_gl
        integer, Intent(in)            :: second_dim_gl
        integer, Intent(in)            :: third_dim_gl
        integer,Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_3D_Int(index).val  => Soln_Local
        Global_Arrays_to_Write_3D_Int(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 30
    
    End Subroutine Assign_Loc_Glob_For_Write_3D_Int
    
End Module Mod_Assign_Loc_Glob_For_Write