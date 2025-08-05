! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SETSHLD(TSC,THETA,D,SSG,DSR,USC)

  ! CHANGE RECORD

  ! ***  NONCOHESIVE SEDIMENT SETTLING AND SHIELDS CRITERIA
  ! ***  USING VAN RIJN'S EQUATIONS

  implicit none

  real :: TSC,THETA,D,SSG,DSR,USC,VISC,GP,TMP,GPD

  if( D <= 1E-8 )then
    TSC = 0.0000055
    USC = 0.0023
    return
  endif

  VISC = 1.D-6
  GP = (SSG-1.)*9.82
  TMP = GP/(VISC*VISC)
  DSR = D*(TMP**0.333333)
  GPD = GP*D

  ! ***  SHIELDS
  if(DSR <= 4.0 )then
    THETA = 0.24/DSR
  endif
  if(DSR > 4.0 .and. DSR <= 10.0 )then
    THETA = 0.14/(DSR**0.64)
  endif
  if(DSR > 10.0 .and. DSR <= 20.0 )then
    THETA = 0.04/(DSR**0.1)
  endif
  if(DSR > 20.0 .and. DSR <= 150.0 )then
    THETA = 0.013*(DSR**0.29)
  endif
  if(DSR > 150.0 )then
    THETA = 0.055
  endif
  TSC = GPD*THETA
  USC = SQRT(TSC)

END SUBROUTINE

