! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION CSEDRESS(DENBULK,WRSPO,VDRO,VDR,VDRC,IOPT)

  ! CHANGE RECORD

  ! ***  CALCULATES SURFACE EROSION RATE OF COHESIVE
  ! ***  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
  ! ***  IOPT = 1  BASED ON
  ! **
  ! ***  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT ERODIBILITY
  ! ***  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
  ! ***  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
  ! ***  IOPT = 2  BASED ON J. M. HAMRICK'S MODIFICATION OF
  ! **
  ! ***  SANFORD, L.P., AND J. P. Y. MAA, 2001: A UNIFIED EROSION FORMULATION
  ! ***  FOR FINE SEDIMENT, MARINE GEOLOGY, 179, 9-23.
  !
  ! ***  IOPT = 4 & 5 BASED ON J. M. HAMRICK'S parameterIZATION OF SEDFLUME
  !     TEST DATA
  ! **

  implicit none

  integer :: IOPT
  real :: CSEDRESS,DENBULK,WRSPO,VDRO,VDR,VDRC
  real :: BULKDEN,TMP,TMPVAL,FACTOR

  if( IOPT == 1 )then
    BULKDEN = 0.001*DENBULK    ! *** TO PREVENT CORRUPTING THE DENBULK VARIABLE
    if( BULKDEN <= 1.065 )then
      CSEDRESS = 0.62
    else
      TMP = 0.198/(BULKDEN-1.0023)
      TMP = EXP(TMP)
      CSEDRESS = 6.4E-4*(10.**TMP)
    endif
  elseif( IOPT == 2 )then
    CSEDRESS = WRSPO*(1.+VDRO)/(1.+VDR)
  elseif( IOPT == 3 )then
    CSEDRESS = WRSPO*(1.+VDRO)/(1.+VDRC)
  elseif( IOPT == 4 )then
    TMPVAL = (1.+VDRO)/(1.+VDR)
    FACTOR = EXP(-TMPVAL)
    CSEDRESS = FACTOR*WRSPO*(1.+VDRO)/(1.+VDR)
  elseif( IOPT == 5 )then
    TMPVAL = (1.+VDRO)/(1.+VDRC)
    FACTOR = EXP(-TMPVAL)
    CSEDRESS = FACTOR*WRSPO*(1.+VDRO)/(1.+VDRC)
  elseif( IOPT >= 99 )then
    CSEDRESS = WRSPO
  endif

END FUNCTION

