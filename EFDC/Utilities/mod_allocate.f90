! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Performs allocation with initialization
! @author Paul M. Craig
! @date 2020-10-31
Module Allocate_Initialize

  Use GLOBAL, only : RKD, RK4
  IMPLICIT NONE
  
  Public :: AllocateDSI

  Interface AllocateDSI

  Module Procedure AllocateDSI_Integer1,   &
                   AllocateDSI_Integer2,   &
                   AllocateDSI_Integer3,   &
                   AllocateDSI_Logical1,   &
                   AllocateDSI_Logical2,   &
                   AllocateDSI_Real1,      &
                   AllocateDSI_Real2,      &
                   AllocateDSI_Real3,      &
                   AllocateDSI_Real4,      &
                   AllocateDSI_RKD1,       &
                   AllocateDSI_RKD2,       &
                   AllocateDSI_RKD3,       &
                   AllocateDSI_RKD4
  End interface

  ! *** Next are each of the subroutines for the procedure
  Contains
 
! *****************************************************************************************
SUBROUTINE AllocateDSI_Integer1(ArrIn, size1, iVal)

  Integer, Allocatable, Intent(inout) :: ArrIn(:)
  Integer, Intent(in)    :: size1, ival
  
  Allocate(ArrIn(size1))
  
  ArrIn = iVal
  
END SUBROUTINE AllocateDSI_Integer1

SUBROUTINE AllocateDSI_Integer2(ArrIn, size1, size2, iVal)

  Integer, Allocatable, Intent(inout) :: ArrIn(:,:)
  Integer, Intent(in)    :: size1, size2, ival
  
  Allocate(ArrIn(size1,size2))
  
  ArrIn = iVal
  
END SUBROUTINE AllocateDSI_Integer2

SUBROUTINE AllocateDSI_Integer3(ArrIn, size1, size2, size3, iVal)

  Integer, Allocatable, Intent(inout) :: ArrIn(:,:,:)
  Integer, Intent(in)    :: size1, size2, size3, ival
  
  Allocate(ArrIn(size1,size2,size3))
  
  ArrIn = iVal
  
END SUBROUTINE AllocateDSI_Integer3

! *****************************************************************************************
SUBROUTINE AllocateDSI_Logical1(ArrIn, size1, lVal)

  Logical, Allocatable, Intent(inout) :: ArrIn(:)
  Integer, Intent(in)    :: size1
  Logical, Intent(in)    :: lval
  
  Allocate(ArrIn(size1))
  
  ArrIn = lVal
  
END SUBROUTINE AllocateDSI_Logical1

SUBROUTINE AllocateDSI_Logical2(ArrIn, size1, size2, lval)

  Logical, Allocatable, Intent(inout) :: ArrIn(:,:)
  Integer, Intent(in)    :: size1, size2
  Logical,    Intent(in)    :: lval
  
  Allocate(ArrIn(size1,size2))
  
  ArrIn = lval
  
END SUBROUTINE AllocateDSI_Logical2

! *****************************************************************************************
SUBROUTINE AllocateDSI_Real1(ArrIn, size1, val)

  Real(RK4), Allocatable, Intent(inout) :: ArrIn(:)
  Integer, Intent(in) :: size1
  Real,    Intent(in) :: val
  
  Allocate(ArrIn(size1))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_Real1

SUBROUTINE AllocateDSI_Real2(ArrIn, size1, size2, val)

  Real(RK4), Allocatable, Intent(inout) :: ArrIn(:,:)
  Integer, Intent(in)    :: size1, size2
  Real,    Intent(in)    :: val
  
  Allocate(ArrIn(size1,size2))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_Real2

SUBROUTINE AllocateDSI_Real3(ArrIn, size1, size2, size3, val)

  Real(RK4), Allocatable, Intent(inout) :: ArrIn(:,:,:)
  Integer, Intent(in)    :: size1, size2, size3
  Real,    Intent(in)    :: val
  
  Allocate(ArrIn(size1,size2,size3))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_Real3

SUBROUTINE AllocateDSI_Real4(ArrIn, size1, size2, size3, size4, val)

  Real(RK4), Allocatable, Intent(inout) :: ArrIn(:,:,:,:)
  Integer, Intent(in)    :: size1, size2, size3, size4
  Real,    Intent(in)    :: val
  
  Allocate(ArrIn(size1,size2,size3, size4))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_Real4

! *****************************************************************************************
SUBROUTINE AllocateDSI_RKD1(ArrIn, size1, val)

  Real(RKD), Allocatable, Intent(inout) :: ArrIn(:)
  Integer, Intent(in) :: size1
  Real,    Intent(in) :: val
  
  Allocate(ArrIn(size1))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_RKD1

SUBROUTINE AllocateDSI_RKD2(ArrIn, size1, size2, val)

  Real(RKD), Allocatable, Intent(inout) :: ArrIn(:,:)
  Integer, Intent(in)    :: size1, size2
  Real,    Intent(in)    :: val
  
  Allocate(ArrIn(size1,size2))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_RKD2

SUBROUTINE AllocateDSI_RKD3(ArrIn, size1, size2, size3, val)

  Real(RKD), Allocatable, Intent(inout) :: ArrIn(:,:,:)
  Integer, Intent(in)    :: size1, size2, size3
  Real,    Intent(in)    :: val
  
  Allocate(ArrIn(size1,size2,size3))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_RKD3

SUBROUTINE AllocateDSI_RKD4(ArrIn, size1, size2, size3, size4, val)

  Real(RKD), Allocatable, Intent(inout) :: ArrIn(:,:,:,:)
  Integer, Intent(in)    :: size1, size2, size3, size4
  Real,    Intent(in)    :: val
  
  Allocate(ArrIn(size1,size2,size3, size4))
  
  ArrIn = Val
  
END SUBROUTINE AllocateDSI_RKD4

END MODULE Allocate_Initialize
