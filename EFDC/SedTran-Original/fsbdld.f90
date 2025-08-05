! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION FSBDLD(DIASED, GPDIASED, D50, DEP, PEXP, PHID, CSHIELDS, SBDLDP, ISOPT)

  ! ***  CALCULATES DIMENSIONLESS BED LOAD TRANSPORT COEFFICIENT (PHI) 
  ! ***  ISOPT = 0  Constant value
  ! ***          1  Van Rijn (1984)
  ! ***          2  Modified Enguland-Hansen
  ! ***          3  Wu, et.al.
  ! ***  FSBDLD = PHI (dimensionless)

  implicit none

  integer,intent(IN) :: ISOPT
  real,intent(IN)    :: DIASED, GPDIASED, D50, DEP, PEXP, PHID, CSHIELDS, SBDLDP
  real :: FSBDLD
  real :: RD, TMP, TMP1, TMP2

  if( ISOPT == 0 )then
    ! ***  ISOPT = 0  USER SPECIFIED
    FSBDLD = SBDLDP
  
  elseif( ISOPT == 1 )then
    ! ***  ISOPT = 1  BASED ON
    ! ***
    ! ***  VAN RIJN, L. C., 1984: SEDIMENT TRANSPORT, PART I: BED
    ! ***  LOAD TRANSPORT, J. HYDRAULIC ENGINEERING, 110, 1431-1455.
    RD  = DIASED*SQRT(GPDIASED)*1.E6
    RD  = (1./RD)**0.2
    TMP = CSHIELDS**2.1
    FSBDLD = 0.053*RD/TMP

  elseif( ISOPT == 2 )then
    ! ***  ISOPT = 2  BASED ON MODIFIED ENGULAND-HANSEN FORMULA
    ! ***  REFERENCE TO BE ADDED
    TMP1   = (DEP/D50)**0.33333
    TMP2   = (PEXP/PHID)**1.125
    FSBDLD = 2.0367*TMP1*TMP2

  elseif( ISOPT == 3 )then
    ! ***  ISOPT = 3  BASED ON WU, WANG AND JIA, J. HYDR. RES. V38, 2000
    TMP1   = 0.03*((PHID/PEXP)**0.6)
    TMP2   = TMP1**2.2
    FSBDLD = 0.0053/TMP2
  endif

END FUNCTION

