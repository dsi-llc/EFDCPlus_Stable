! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

    Use Global
    Use Variables_MPI
    Use Variables_MPI_Mapping
    Use Variables_MPI_MapGatherSort
    
    Implicit None

    Save
    
    Public :: Assign_Loc_Glob_For_Write
    
    Interface Assign_Loc_Glob_For_Write
    
        Module Procedure Assign_Loc_Glob_For_Write_1D_Real, &
                         Assign_Loc_Glob_For_Write_1D_Real_RK8, &
                         Assign_Loc_Glob_For_Write_1D_Int,  &
                         Assign_Loc_Glob_For_Write_2D_Real, &
                         Assign_Loc_Glob_For_Write_2D_Real_RK8, &
                         Assign_Loc_Glob_For_Write_2D_Int,  &                         
                         Assign_Loc_Glob_For_Write_3D_Real, &
                         Assign_Loc_Glob_For_Write_3D_Real_RK8, &                    
                         Assign_Loc_Glob_For_Write_3D_Int
    
        
    End Interface

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
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in) :: index
        Integer, Intent(in) :: first_dim_local
        Real(4), Target, Intent(inout) :: Soln_Local(first_dim_local)
        Integer, Intent(in) :: first_dim_gl
        Real(4),  Target,   Intent(inout) :: Soln_Global(first_dim_gl)
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
        Use Global
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in) :: index
        Integer, Intent(in) :: first_dim_local
        Real(rkd), Target, Intent(inout) :: Soln_Local(first_dim_local)
        Integer, Intent(in) :: first_dim_gl
        Real(rkd),  Target,   Intent(inout) :: Soln_Global(first_dim_gl)
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
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in) :: index
        Integer, Intent(in) :: first_dim_local
        Integer,Target, Intent(inout) :: Soln_Local(first_dim_local)
        Integer, Intent(in) :: first_dim_gl
        Integer,Target,  Intent(inout) :: Soln_Global(first_dim_gl)
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
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in) :: index
        Integer, Intent(in) :: first_dim_local
        Integer, Intent(in) :: second_dim_local
        Real(4), Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local)
        Integer, Intent(in) :: first_dim_gl
        Integer, Intent(in) :: second_dim_gl
        Real(4), Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl)
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
        Use Global
                                        
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in) :: index
        Integer, Intent(in) :: first_dim_local
        Integer, Intent(in) :: second_dim_local
        Real(rkd), Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local)
        Integer, Intent(in) :: first_dim_gl
        Integer, Intent(in) :: second_dim_gl
        Real(rkd), Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl)
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
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in) :: index
        Integer, Intent(in) :: first_dim_local
        Integer, Intent(in) :: second_dim_local
        Integer, Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local)
        Integer, Intent(in) :: first_dim_gl
        Integer, Intent(in) :: second_dim_gl
        Integer, Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl)
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
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in)         :: index
        Integer, Intent(in)         :: first_dim_local
        Integer, Intent(in)         :: second_dim_local
        Integer, Intent(in)         :: third_dim_local
        Real(4), Target, Intent(inout) :: Soln_Local(first_dim_local, second_dim_local, third_dim_local)
        Integer, Intent(in)         :: first_dim_gl
        Integer, Intent(in)         :: second_dim_gl
        Integer, Intent(in)         :: third_dim_gl
        Real(4), Target, Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
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
        Use Global
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in)         :: index
        Integer, Intent(in)         :: first_dim_local
        Integer, Intent(in)         :: second_dim_local
        Integer, Intent(in)         :: third_dim_local
        Real(rkd), Target, Intent(inout) :: Soln_Local(first_dim_local, second_dim_local, third_dim_local)
        Integer, Intent(in)         :: first_dim_gl
        Integer, Intent(in)         :: second_dim_gl
        Integer, Intent(in)         :: third_dim_gl
        Real(rkd), Target, Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
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
                                                 
        Implicit None
        
        ! *** Read in 
        Integer, Intent(in)            :: index
        Integer, Intent(in)            :: first_dim_local
        Integer, Intent(in)            :: second_dim_local
        Integer, Intent(in)            :: third_dim_local
        Integer,Target,  Intent(inout) :: Soln_Local(first_dim_local, second_dim_local, third_dim_local)
        Integer, Intent(in)            :: first_dim_gl
        Integer, Intent(in)            :: second_dim_gl
        Integer, Intent(in)            :: third_dim_gl
        Integer,Target,  Intent(inout) :: Soln_Global(first_dim_gl, second_dim_gl, third_dim_gl)
        ! *** Local variables
    
        Local_Arrays_to_Write_3D_Int(index).val  => Soln_Local
        Global_Arrays_to_Write_3D_Int(index).val => Soln_Global
        Dim_Array_Written_Out(index) = 30
    
    End Subroutine Assign_Loc_Glob_For_Write_3D_Int
    
End Module Mod_Assign_Loc_Glob_For_Write