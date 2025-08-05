! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)

  ! ***  FROM NUMERICAL RECIPES
  ! CHANGE RECORD
  !
  ! 2014-08           D H CHUNG        SET EXPLICIT PRECISIONS OF INTEGER & REAL
  implicit none

  integer :: M,N,I,J,JJ
  integer :: MP,NP
  
  real :: U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),S
  real,save,allocatable :: TMP(:)
  
  if( .not. allocated(TMP) )then
    allocate(TMP(N))
    TMP = 0.
  endif
  
  do 12 J = 1,N
    S = 0.
    if( W(J) /= 0. )then
      do 11 I = 1,M
        S = S+U(I,J)*B(I)
      11 continue
      S = S/W(J)
    endif
    TMP(J) = S
  12 continue

  do 14 J = 1,N
    S = 0.
    do 13 JJ = 1,N
      S = S+V(J,JJ)*TMP(JJ)
    13   continue
    X(J) = S
  14 continue
  deallocate(TMP)

END SUBROUTINE

