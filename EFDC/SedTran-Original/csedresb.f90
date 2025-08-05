! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION CSEDRESB(DENBULK,WRSPO,VDRO,VDR,VDRC,IOPT)

  ! ***  CALCULATES BULK EROSION RATE OF COHESIVE
  ! ***  SEDIMENT AS A FUNCTION OF BED BULK DENSITY AND OTHER VARIABLES
  ! ***  CURRENT OPTIONS SHOULD NOT BE USED
  ! ***  IOPT = 1  BASED ON
  ! **
  ! ***  HWANG, K. N., AND A. J. MEHTA, 1989: FINE SEDIMENT ERODIBILITY
  ! ***  IN LAKE OKEECHOBEE FLORIDA. COASTAL AND OCEANOGRAPHIC ENGINEERING
  ! ***  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661
  ! ***  IOPT = 2  BASED ON J. M. HAMRICK'S MODIFICATION OF
  ! **
  ! ***  SANFORD, L.P., AND J. P. Y. MAA, 2001: A UNIFIED EROSION FORMULATI
  ! ***  FOR FINE SEDIMENT, MARINE GEOLOGY, 179, 9-23.
  ! CHANGE RECORD

  implicit none
  
  integer :: IOPT
  real    :: CSEDRESB,DENBULK,WRSPO,VDRO,VDR,VDRC

  if( IOPT >= 1 ) CSEDRESB = 0.0
  !
  !        DENBULK = 0.001*DENBULK
  !          CSEDRESB = 0.62
  !        else
  !          TMP = 0.198/(DENBULK-1.0023)
  !          CSEDRESB = 6.4E-4*(10.**TMP)
  !

END FUNCTION

