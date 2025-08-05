! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION SETSTVEL(D,SSG)

  ! ***  NONCOHEASIVE SEDIMENT SETTLING AND SHIELDS CRITERIA
  ! ***  USING VAN RIJN'S EQUATIONS

  ! CHANGE RECORD

  implicit none

  real :: SETSTVEL,D,SSG,VISC,GP,GPD,SQGPD,RD,WSET,TMP

  VISC = 1.D-6
  GP = (SSG-1.)*9.82
  GPD = GP*D
  SQGPD = SQRT(GPD)
  RD = SQGPD*D/VISC

  ! ***  SETTLING VELOCITY

  if( D < 1.0E-4 )then
    WSET = SQGPD*RD/18.
  endif
  if( D >= 1.0E-4 .and. D < 1.E-3 )then
    TMP = SQRT(1.+0.01*RD*RD)-1.
    WSET = 10.0*SQGPD*TMP/RD
  endif
  if( D >= 1.E-3 )then
    WSET = 1.1*SQGPD
  endif
  SETSTVEL = WSET

END FUNCTION

