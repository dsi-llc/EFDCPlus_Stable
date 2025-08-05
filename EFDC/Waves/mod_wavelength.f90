! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE WAVELENGTH
  !Author: Dang Huu Chung

  use GLOBAL,only:RKD
  implicit none

  real(RKD),PRIVATE,parameter :: G = 9.81
  real(RKD),PRIVATE,parameter :: PI = 3.14159265358979

  contains

  FUNCTION DISRELATION(RLS,TP,HD,U,PHI) RESULT(FWL)
    !Dispersion Relation:          FWL = 0
    !RLS: Wave length              [m]
    !TP : Wave period              [s]
    !HD : Water depth              [m]
    !U  : Depth-average velocity   [m/s]
    !PHI: Angle of (wave,current)  [Rad]
    real(RKD) :: FWL,RLS,TP,HD,U,PHI
    FWL = (RLS/TP-U*COS(PHI))-SQRT(G*RLS/2._8/PI*TANH(2._8*PI*HD/RLS))
  END FUNCTION

  RECURSIVE SUBROUTINE BISEC(FUN,A0,B0,TOL,TP,HD,U,PHI,X)
    real(RKD),intent(IN)  :: A0,B0,TOL,TP,HD,U,PHI
    real(RKD),external    :: FUN
    real(RKD),intent(OUT) :: X
    integer               :: IST
    real(RKD)             :: A,B,FA,FB,FX
    A = A0
    B = B0
    IST = 1
    FA = FUN(A,TP,HD,U,PHI)
    FB = FUN(B,TP,HD,U,PHI)
    if( FA*FB < 0 )then
      X = 0.5*(A+B)
      FX = FUN(X,TP,HD,U,PHI)
      if( FA*FX < 0 )then
        B = X
      else
        A = X
      endif
      if( ABS(A-B) <= TOL )then
        return
      else
        call BISEC(FUN,A,B,TOL,TP,HD,U,PHI,X)
      endif
    elseif( FA == 0 )then
      X = A
    elseif( FB == 0 )then
      X = B
    else
      IST = -1
      call STOPP('DISPERSION RELATION: FA.FB>0')
    endif
  END SUBROUTINE BISEC
  
  FUNCTION RTBIS(FUNC,X1,X2,XACC,TP,HD,U,PHI)
    real(RKD) :: RTBIS,X1,X2,XACC,TP,HD,U,PHI
    real(RKD),external :: FUNC
    integer   :: J,JMAX
    real(RKD) :: DX,F,FMID,XMID
    parameter (JMAX = 40)
  
    FMID = FUNC(X2,TP,HD,U,PHI)
    F = FUNC(X1,TP,HD,U,PHI)
    if( F*FMID >= 0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
    if( F<0. )then
      RTBIS = X1
      DX = X2-X1
    else
      RTBIS = X2
      DX = X1-X2
    endif
    do J = 1,JMAX
      DX = DX*.5
      XMID = RTBIS+DX
      FMID = FUNC(XMID,TP,HD,U,PHI)
      if( FMID <= 0.)RTBIS = XMID
      if( ABS(DX)<XACC .or. FMID == 0. ) return
    enddo
    PAUSE 'TOO MANY BISECTIONS IN RTBIS'
  END FUNCTION

END MODULE


   
