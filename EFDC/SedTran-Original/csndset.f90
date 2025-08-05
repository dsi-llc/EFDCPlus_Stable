! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION CSNDSET(SND,SDEN,IOPT)
  !
  ! CHANGE RECORD
  ! ***  CALCULATES HINDERED SETTLING CORRECTION FOR CLASS NS NONCOHESIVE
  ! ***  SEDIMENT

  implicit none
  
  integer :: IOPT
  real :: CSNDSET,SND,SDEN,ROPT

  ROPT = FLOAT(IOPT)
  CSNDSET = (1.-SDEN*SND)**ROPT

END FUNCTION

