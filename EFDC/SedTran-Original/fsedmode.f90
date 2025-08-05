! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION FSEDMODE(WS, USTOT, USGRN, RSNDM, ISNDM1, ISNDM2, IMODE)

  ! *** FSEDMODE SETS BEDLOAD (IMODE = 1) AND SUSPENDED LOAD (IMODE = 2)
  ! *** TRANSPORT FRACTIONS (DIMENSIONLESS)
  
  ! ***  WS     - SETTING VELOCITY (M/S)
  ! ***  USTOT  - HYDRODYNAMIC SHEAR VELOCITY (M/S)
  ! ***  USGRN  - GRAIN SHEAR VELOCITY (M/S)
  ! ***  RSNDM  - VALUE OF USTAR/WSET FOR BINARY SWITCH BETWEEN BEDLOAD AND SUSPENDED LOAD (DIMENSIONLESS)

  implicit none

  integer :: ISNDM1,ISNDM2,IMODE
  real   :: FSEDMODE, WS, USTOT, USGRN, RSNDM, US, USDWS, TMPVAL

  if( WS == 0. )then
    FSEDMODE = 0.
    return
  endif

  ! *** CHOOSE BETWEEEN TOTAL STRESS SHEAR VELOCITY AND GRAIN STRESS
  ! *** SHEAR VELOCITY
  if( ISNDM2 == 0 )then
    US = USTOT
  else
    US = USGRN
  endif
  USDWS = US/WS

  ! *** CHOOSE BETWEEN MODE OPTIONS

  if( ISNDM1 == 0 )then
    ! *** ISNDM1 = 0 SET BOTH BEDLOAD AND SUSPENDED LOAD FRACTIONS TO 1.0
    FSEDMODE = 1.0

  elseif( ISNDM1 == 1 )then
    ! *** ISNDM1 = 1 SET BEDLOAD FRACTION TO 1. use BINARY RELATIONSHIP FOR SUSPENDED
    FSEDMODE = 0.
    if( IMODE == 1 )then
      FSEDMODE = 1.0
    else
      if( USDWS >= RSNDM ) FSEDMODE = 1.
    endif

  elseif( ISNDM1 == 2 )then
    ! *** ISNDM1 = 2 SET BEDLOAD FRACTION TO 1, use LINEAR RELATIONSHIP FOR SUSPENDED
    if( IMODE == 1 )then
      FSEDMODE = 1.0
    else
      TMPVAL = ((USDWS)-0.4)/9.6
      TMPVAL = MIN(TMPVAL,1.0)
      TMPVAL = MAX(TMPVAL,0.0)
      FSEDMODE = TMPVAL
    endif

  elseif( ISNDM1 == 3 )then
    ! *** ISNDM1 = 3 use BINARY RELATIONSHIP FOR BEDLOAD AND SUSPENDED LOAD
    FSEDMODE = 0.
    if( IMODE == 1 )then
      if( USDWS < RSNDM ) FSEDMODE = 1.
    else
      if( USDWS >= RSNDM ) FSEDMODE = 1.
    endif

  elseif( ISNDM1 == 4 )then
    ! *** ISNDM1 = 4 use LINEAR RELATIONSHIP FOR BEDLOAD AND SUSPENDED LOAD
    TMPVAL = ((USDWS)-0.4)/9.6
    TMPVAL = MIN(TMPVAL,1.0)
    TMPVAL = MAX(TMPVAL,0.0)
    if( IMODE == 1 )then
      FSEDMODE = 1.-TMPVAL
    else
      FSEDMODE = TMPVAL
    endif
  endif

END FUNCTION
