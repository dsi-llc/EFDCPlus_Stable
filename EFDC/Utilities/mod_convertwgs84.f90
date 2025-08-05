! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE CONVERTWGS84
  ! *** CONVERT GEOGRAPHIC COORDINATES TO UTM AND REVERSE
  ! *** AUTHOR: DH CHUNG
  ! *** START : 2008
  ! *** UPDATE: 2016-06-16

  use GLOBAL,only:RKD,PI,HEMI,UTMZ

  implicit none

  ! *** HEMI = HEMISPHERE  (1:NORTH/2SOUTH)
  ! *** UTMZ = ZONE NUMBER (1:60)
  ! *** GEOGRAPHY SYSTEM:
  ! *** 1.WGS84/NAD83 /2.GRS80 /3.WGS72

  character(20):: GEOSYS = '1.WGS84/NAD83'
  real(RKD):: R_MJR      !SEMI MAJOR AXIS
  real(RKD):: R_MNR      !SEMI MINOR AXIS
  real(RKD):: SCLF       !SCALE ALONG CENTRAL MERIDIAN
  real(RKD):: FE         !X OFFSET IN METER
  real(RKD):: FN         !Y OFFSET IN METER
  real(RKD):: LAM0       !CENTRAL MERIDIAN OF ZONE    [DEGREE]
  real(RKD):: FLA,E2,EP2

  contains

  SUBROUTINE UTMPARS
  ! *** UTM parameterS

  if( HEMI==1 )then
    FN = 0
  elseif(HEMI==2 )then
    FN = 1.D7
  endif

  FE = 500000._RKD
  SCLF = 0.9996_RKD

  if( GEOSYS(1:1)=='1' )then
    R_MJR = 6378137._RKD 
    FLA  = 1._RKD/298.257223563_RKD
  elseif(GEOSYS(1:1)=='2' )then
    R_MJR = 6378137._RKD
    FLA  = 1._RKD/298.257222101_RKD
  elseif(GEOSYS(1:1)=='3' )then
    R_MJR = 6378135._RKD
    FLA  = 1._RKD/298.26_RKD
  endif

  E2 = 2*FLA-FLA**2
  EP2 = E2/(1-E2)
  R_MNR = R_MJR*SQRT(1-E2)
  LAM0 = 6*(UTMZ-30)-3    !CENTRAL MERIDIAN IN DEGREE

  END SUBROUTINE

  SUBROUTINE UTM_WGS84(LON,LAT,XUTM,YUTM)
  ! *********************************************************
  ! *** INPUT:
  ! LAT   = LATITUDE OF POINTS          [DECIMAL DEGREE]
  ! LON   = LONGITUDE OF POINTS         [DECIMAL DEGREE]
  ! *** OUTPUT:
  ! XUTM = ABSCISSA OF POINTS: XUTM    [M] 
  ! YUTM = ORDINATE OF POINTS: YUTM    [M]

  real(RKD)  ,intent(IN )::LON(:),LAT(:)
  real(RKD)  ,intent(OUT)::XUTM(:),YUTM(:)
  real(RKD):: F2,AP,BP,CP,DP,EP,L0
  real(RKD),dimension(SIZE(LON))::A2,A4,A6,B1,B3,B5,DLAM,RN,LEN1,TAPH,ETA2,LAM,PHI

  L0  = PI*LAM0/180
  LAM = PI*LON/180
  PHI = PI*LAT/180
  DLAM= LAM-L0

  RN = R_MJR/SQRT(1-E2*SIN(PHI)**2)
  F2 = (R_MJR-R_MNR)/(R_MJR+R_MNR)

  AP = R_MJR*(1-F2+(5._RKD/4)*(F2**2-F2**3)+(81._RKD/64)*(F2**4-F2**5))
  BP = (3._RKD/2)*R_MJR*F2*(1-F2+(7._RKD/8)*(F2**2-F2**3)+(55._RKD/64)*(F2**4-F2**5))
  CP = (15._RKD/16)*R_MJR*F2**2*(1-F2+(3._RKD/4)*(F2**2-F2**3))
  DP = (35._RKD/48)*R_MJR*F2**3*(1-F2+(11._RKD/16)*(F2**2-F2**3))
  EP = (315._RKD/51)*(R_MJR*F2**4)*(1-F2)
  LEN1 = AP*PHI-BP*SIN(2*PHI)+CP*SIN(4*PHI)-DP*SIN(6*PHI)+EP*SIN(8*PHI)
  ETA2 = EP2*COS(PHI)**2
  TAPH = TAN(PHI)**2
  A2 = RN/2*SIN(PHI)*COS(PHI)
  A4 = RN/24*SIN(PHI)*COS(PHI)**3*(5-TAPH+9*ETA2+4*ETA2**2)
  A6 = RN/720*SIN(PHI)*COS(PHI)**5*(61-58*TAPH+TAPH**2+270*ETA2-330*ETA2*TAPH)
  B1 = RN*COS(PHI)
  B3 = RN/6*COS(PHI)**3*(1-TAPH+ETA2)
  B5 = RN/120*COS(PHI)**5*(5-18*TAPH+TAPH**2+14*ETA2-58*ETA2*TAPH+13*ETA2**2-64*ETA2**2*TAPH)

  YUTM = SCLF*(LEN1+A2*DLAM**2+A4*DLAM**4+A6*DLAM**6)+FN
  XUTM = SCLF*(B1*DLAM+B3*DLAM**3+B5*DLAM**5)+FE

  END SUBROUTINE

  SUBROUTINE UTMR_WGS84(XUTM,YUTM,LON,LAT)
  ! *** INPUT:
  ! XUTM           [M] 
  ! YUTM           [M]
  ! *** OUTPUT:
  ! LON  = LONGITUDE (DECIMAL DEGREE)
  ! LAT  = LATITUDE  (DECIMAL DEGREE)

  real(RKD), intent(IN )::XUTM(:),YUTM(:)
  real(RKD), intent(OUT)::LON(:),LAT(:)
  real(RKD)::E1,J1,J2,J3,J4,L0
  real(RKD),dimension(SIZE(XUTM))::RM,MU,BX,ETAX2,VX2,NX,TX,A2,A4,A6,B1,B3,B5,Y

  L0  = PI*LAM0/180
  RM  = YUTM/SCLF
  MU = RM/(R_MJR*(1-E2/4-3*E2**2/64 - 5*E2**3/256))

  E1 = (1-SQRT(1-E2))/(1+SQRT(1-E2))
  J1 = 3*E1/2 - 27*E1**3/32 
  J2 = 21*E1**2/16 - 55*E1**4/32
  J3 = 151*E1**3/96
  J4 = 1097*E1**4/512
  BX = MU+J1*SIN(2*MU)+J2*SIN(4*MU)+J3*SIN(6*MU)+J4*SIN(8*MU)
  ETAX2 = EP2*COS(BX)**2
  VX2 = 1+ETAX2
  NX  = R_MJR/SQRT(1-E2*SIN(BX)**2)
  TX  = TAN(BX)
  A2  = -VX2*TX/(2*NX**2)
  A4  = -A2/(12*NX**2)*(5+3*TX**2+ETAX2-9*ETAX2*TX**2-4*ETAX2**2)
  A6  = -A2/(360*NX**4)*(61-90*TX**2+45*TX**4-46*ETAX2)

  B1 = 1/(NX*COS(BX))
  B3  = -B1/(6*NX**2)*(1+2*TX**2+ETAX2)
  B5  = -B1/(120*NX**4)*(5+28*TX**2+24*TX**4+6*ETAX2+8*ETAX2*TX**2)
  Y   = (XUTM-FE)/SCLF
  LAT = BX+A2*Y**2+A4*Y**4+A6*Y**6
  LON = B1*Y+B3*Y**3+B5*Y**5 +L0
  LAT = LAT*180/PI ! TO DEGREE
  LON = LON*180/PI ! TO DEGREE

  END SUBROUTINE

END MODULE
