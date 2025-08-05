! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION CSEDTAUB(DENBULK,IOPT)

  ! ***  CALCULATES CRITICAL STRESS FOR BULK OR MASS EROSION OF COHESIVE
  ! ***  SEDIMENT AS A FUNCTION OF BED BULK DENSITY
  ! **
  ! ***  IOPT = 1  BASED ON
  ! ***  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
  ! ***  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
  ! ***  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
  ! **
  ! ***  IOPT = 2  BASED ON
  ! ***  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT DRODIBILITY
  ! ***  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
  ! ***  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661

  ! CHANGE RECORD

  implicit none
  
  integer :: IOPT
  real    :: CSEDTAUB, DENBULK, BULKDEN

  if( IOPT == 1 )then
    BULKDEN = 0.001*DENBULK  ! *** PMC Changed to prevent
    if( BULKDEN <= 1.013 )then
      CSEDTAUB = 0.0
    else
      CSEDTAUB = 0.001*(9.808*BULKDEN-9.934)
    endif
  elseif( IOPT == 2 )then
    BULKDEN = 0.001*DENBULK  ! *** PMC Changed to prevent
    if( BULKDEN <= 1.013 )then
      CSEDTAUB = 0.0
    else
      CSEDTAUB = 0.001*(9.808*BULKDEN-9.934)
    endif
  else
    call STOPP('CSEDTAUB: BAD SEDIMENT CRITICAL STRESS OPTION! STOPPING!')
  endif

END FUNCTION

