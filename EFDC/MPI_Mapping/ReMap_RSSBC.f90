! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!  
!                     EFDC+ Developed by DSI, LLC. 
!---------------------------------------------------------------------------!  
! Subroutine: Map_RSSBC
!
!> @details Maps the local RSSBC* values to global ones in the SETBCS routine
!
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!---------------------------------------------------------------------------!  
! CHANGE RECORD  
! DATE MODIFIED     BY               DESCRIPTION
!---------------------------------------------------------------------------!
!    2019-10      Zander     Maps the local RSSBC* values to global ones 
!                           
!---------------------------------------------------------------------------!
Subroutine ReMap_RSSBC

    Use GLOBAL
    Use Variables_MPI
    Use Variables_MPI_Write_Out
    Use Variables_MPI_Mapping
    Use Mod_Map_Soln
    Use Mod_Gather_Soln
    Use Mod_Sort_Global_Soln
    
    Implicit none

    ! *** Passed in
    
    ! *** Local Variables
    Real,    Allocatable, Dimension(:) :: Local_1D
    Real,    Allocatable, Dimension(:) :: Global_1D
    Integer, Allocatable, Dimension(:) :: Mapping_Local
    Integer, Allocatable, Dimension(:) :: Mapping_Global
    Integer :: num_active_l_local_1d
    Integer :: l
    
    ! *** Allocate local variables
    Allocate(Local_1D(LA_Local_no_ghost))
    Allocate(Global_1D(LCM_Global))
    Allocate(Mapping_Local(LA_Local_no_ghost))
    Allocate(Mapping_Global(LCM_Global))
    
    ! *** Initialize to zero
    Local_1D  = 0.0
    Global_1D = 0.0
    Mapping_Local = 0
    Mapping_Global = 0

    ! *** East (RSSBCE)

    ! *** Map to 1d array containing l values excluding ghost cells and build up array mapping to global l values
    Call Map_Soln(LCM, RSSBCE, LA_Local_no_ghost, Local_1D, Mapping_Local)
    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
    Call Gather_Soln(LA_Local_no_ghost, Local_1D, LCM_Global, Global_1D, Mapping_Local, Mapping_Global ) 
    ! *** Sort the global array so L indexing is correct
    Call Sort_Global_Soln(LCM_Global, LCM_Global, Global_1D, Mapping_Global, RSSBCE_Global)

    Local_1D  = 0.0
    Global_1D = 0.0
    Mapping_Local = 0
    Mapping_Global = 0
    
    ! *** West (RSSBCW)
    ! *** Map to 1d array containing l values excluding ghost cells and build up array mapping to global l values
    Call Map_Soln(LCM, RSSBCW, LA_Local_no_ghost, Local_1D, Mapping_Local)
    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
    Call Gather_Soln(LA_Local_no_ghost, Local_1D, LCM_Global, Global_1D, Mapping_Local, Mapping_Global ) 
    ! *** Sort the global array so L indexing is correct
    Call Sort_Global_Soln(LCM_Global, LCM_Global, Global_1D, Mapping_Global, RSSBCW_Global)
    
    Local_1D  = 0.0
    Global_1D = 0.0
    Mapping_Local = 0
    Mapping_Global = 0
    
    ! *** North (RSSBCN)
    ! *** Map to 1d array containing l values excluding ghost cells and build up array mapping to global l values
    Call Map_Soln(LCM, RSSBCN, LA_Local_no_ghost, Local_1D, Mapping_Local)
    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
    Call Gather_Soln(LA_Local_no_ghost, Local_1D, LCM_Global, Global_1D, Mapping_Local, Mapping_Global ) 
    ! *** Sort the global array so L indexing is correct
    Call Sort_Global_Soln(LCM_Global, LCM_Global, Global_1D, Mapping_Global, RSSBCN_Global)

    Local_1D  = 0.0
    Global_1D = 0.0
    Mapping_Local  = 0
    Mapping_Global = 0
    
    ! *** South (RSSBCS)
    ! *** Map to 1d array containing l values excluding ghost cells and build up array mapping to global l values
    Call Map_Soln(LCM, RSSBCS, LA_Local_no_ghost, Local_1D, Mapping_Local)
    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
    Call Gather_Soln(LA_Local_no_ghost, Local_1D, LCM_Global, Global_1D, Mapping_Local, Mapping_Global ) 
    ! *** Sort the global array so L indexing is correct
    Call Sort_Global_Soln(LCM_Global, LCM_Global, Global_1D, Mapping_Global, RSSBCS_Global)
    
    ! *** End mapping to global
    Deallocate(Local_1D)
    Deallocate(Global_1D)
    Deallocate(Mapping_Local)
    Deallocate(Mapping_Global)
    
End subroutine ReMap_RSSBC
