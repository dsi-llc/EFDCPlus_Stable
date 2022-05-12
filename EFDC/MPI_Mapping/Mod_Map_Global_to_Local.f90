! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! *****************************************************************************************************
! *****************************************************************************************************
! Module: Map_Global_to_Local
!
!> @details Maps global L array (LCM_Global) to local L array (LCM)
!
!> @param[inout] Array
!
!---------------------------------------------------------------------------!
!> @author Paul M. Craig
!> @date   2020-02-26
!---------------------------------------------------------------------------!
Module Mod_Map_Global_to_Local

  Use GLOBAL
  Use Variables_MPI
  Use Broadcast_Routines
    
  Implicit None

  Save
    
  Public :: Map_Global_to_Local
    
  Interface Map_Global_to_Local
    
    Module Procedure Map_Global_to_Local_RK4, &
                     Map_Global_to_Local_RK8, &
                     Map_Global_to_Local_I4,  &
                     Map_Global_to_Local_L4,  &
                     Map_Global_to_Local_R2D, &
                     Map_Global_to_Local_R2D_RK8, &
                     Map_Global_to_Local_I2D, &
                     Map_Global_to_Local_GRP, &
                     Map_Global_to_Local_GRP_RK8
    
  End Interface

  Contains
    
  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_RK4(Array)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Real(4), Intent(inout), Allocatable, Dimension(:) :: Array
        
    ! *** Local variables
    Integer :: L, LG
    Real(4), Allocatable, Dimension(:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global))
    ArrayG(:) = 0.
    
    do LG = 1,LCM_Global
      ArrayG(LG) = Array(LG)
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM))
    Array = 0.
    
    ! *** Now map global to local
    do LG = 1,LCM_Global
      L = Map2Local(LG).LL
      if( L > 0 )then
        Array(L) = ArrayG(LG)
      endif
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_RK4
    
  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_RK8(Array)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Real(8), Intent(inout), Allocatable, Dimension(:) :: Array
        
    ! *** Local variables
    Integer :: L, LG
    Real(8), Allocatable, Dimension(:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global))
    ArrayG(:) = 0.
    
    do LG = 1,LCM_Global
      ArrayG(LG) = Array(LG)
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM))
    Array = 0.
    
    ! *** Now map global to local
    do LG = 1,LCM_Global
      L = Map2Local(LG).LL
      if( L > 0 )then
        Array(L) = ArrayG(LG)
      endif
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_RK8

    
  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_I4(Array)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Integer(4), Intent(inout), Allocatable, Dimension(:) :: Array
        
    ! *** Local variables
    Integer :: L, LG
    Integer(4), Allocatable, Dimension(:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global))
    ArrayG(:) = 0.
    
    do LG = 1,LCM_Global
      ArrayG(LG) = Array(LG)
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM))
    Array = 0
    
    ! *** Now map global to local
    do LG = 1,LCM_Global
      L = Map2Local(LG).LL
      if( L > 0 )then
        Array(L) = ArrayG(LG)
      endif
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_I4

    
  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_L4(Array)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Logical, Intent(inout), Allocatable, Dimension(:) :: Array
        
    ! *** Local variables
    Integer :: L, LG
    Logical, Allocatable, Dimension(:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global))
    ArrayG(:) = 0.
    
    do LG = 1,LCM_Global
      ArrayG(LG) = Array(LG)
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM))
    Array = 0.
    
    ! *** Now map global to local
    do LG = 1,LCM_Global
      L = Map2Local(LG).LL
      if( L > 0 )then
        Array(L) = ArrayG(LG)
      endif
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_L4

  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_R2D(Array,KMAX)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Integer, Intent(in) :: KMAX
    Real(4), Intent(inout), Allocatable, Dimension(:,:) :: Array
        
    ! *** Local variables
    Integer :: L, LG, K
    Real(4), Allocatable, Dimension(:,:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global,KMAX))
    ArrayG = 0.
    
    do K = 1,KMAX
      do LG = 1,LCM_Global
        ArrayG(LG,K) = Array(LG,K)
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM,KMAX))
    Array = 0.
    
    ! *** Now map global to local
    do K = 1,KMAX
      do LG = 1,LCM_Global
        L = Map2Local(LG).LL
        if( L > 0 )then
          Array(L,K) = ArrayG(LG,K)
        endif
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_R2D

  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_R2D_RK8(Array,KMAX)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Integer, Intent(in) :: KMAX
    Real(8), Intent(inout), Allocatable, Dimension(:,:) :: Array
        
    ! *** Local variables
    Integer :: L, LG, K
    Real(8), Allocatable, Dimension(:,:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global,KMAX))
    ArrayG = 0.
    
    do K = 1,KMAX
      do LG = 1,LCM_Global
        ArrayG(LG,K) = Array(LG,K)
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM,KMAX))
    Array = 0.
    
    ! *** Now map global to local
    do K = 1,KMAX
      do LG = 1,LCM_Global
        L = Map2Local(LG).LL
        if( L > 0 )then
          Array(L,K) = ArrayG(LG,K)
        endif
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_R2D_RK8
  
  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_I2D(Array,KMAX)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Integer, Intent(in) :: KMAX
    Integer, Intent(inout), Allocatable, Dimension(:,:) :: Array
        
    ! *** Local variables
    Integer :: L, LG, K
    Integer, Allocatable, Dimension(:,:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global,KMAX))
    ArrayG = 0
    
    do K = 1,KMAX
      do LG = 1,LCM_Global
        ArrayG(LG,K) = Array(LG,K)
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM,KMAX))
    Array = 0.
    
    ! *** Now map global to local
    do K = 1,KMAX
      do LG = 1,LCM_Global
        L = Map2Local(LG).LL
        if( L > 0 )then
          Array(L,K) = ArrayG(LG,K)
        endif
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_I2D

  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_GRP(Array,KMAX,NGRP)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Integer, Intent(in) :: KMAX, NGRP
    Real(4), Intent(inout), Allocatable, Dimension(:,:,:) :: Array
        
    ! *** Local variables
    Integer :: L, LG, NS, K
    Real(4), Allocatable, Dimension(:,:,:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global,KMAX,NGRP))
    ArrayG = 0.
    
    do NS = 1,NGRP
      do K = 1,KMAX
        do LG = 1,LCM_Global
          ArrayG(LG,K,NS) = Array(LG,K,NS)
        enddo
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM,KMAX,NGRP))
    Array = 0.
    
    ! *** Now map global to local
    do NS = 1,NGRP
      do K = 1,KMAX
        do LG = 1,LCM_Global
          L = Map2Local(LG).LL
          if( L > 0 )then
            Array(L,K,NS) = ArrayG(LG,K,NS)
          endif
        enddo
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_GRP
  !---------------------------------------------------------------------------!  
  Subroutine Map_Global_to_Local_GRP_RK8(Array,KMAX,NGRP)

    Use GLOBAL
    Use Variables_MPI
    
    Implicit none
    
    ! *** Passed in variables
    Integer, Intent(in) :: KMAX, NGRP
    Real(8), Intent(inout), Allocatable, Dimension(:,:,:) :: Array
        
    ! *** Local variables
    Integer :: L, LG, NS, K
    Real(8), Allocatable, Dimension(:,:,:) :: ArrayG
    
    ! *** First Get the global array from the master_id
    Call Broadcast_Array(Array, master_id)

    ! *** Setup temporary arrays for transfer
    Allocate(ArrayG(LCM_Global,KMAX,NGRP))
    ArrayG = 0.
    
    do NS = 1,NGRP
      do K = 1,KMAX
        do LG = 1,LCM_Global
          ArrayG(LG,K,NS) = Array(LG,K,NS)
        enddo
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(Array)
    
    ! *** Reallocate to local LCM
    Allocate(Array(LCM,KMAX,NGRP))
    Array = 0.
    
    ! *** Now map global to local
    do NS = 1,NGRP
      do K = 1,KMAX
        do LG = 1,LCM_Global
          L = Map2Local(LG).LL
          if( L > 0 )then
            Array(L,K,NS) = ArrayG(LG,K,NS)
          endif
        enddo
      enddo
    enddo
    
    ! *** Release current array memory
    deallocate(ArrayG)

  End Subroutine Map_Global_to_Local_GRP_RK8
  
  
End Module Mod_Map_Global_to_Local
  