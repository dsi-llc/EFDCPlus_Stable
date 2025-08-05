! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION FPROBDEP(TAUDDD,TAUBBB)

  ! *** *******************************************************************C
  !
  ! ***  FPROBDEP CALCULATES PROBABILITY OF DEPOSITION USING PROBABILITY
  ! ***  INTEGRAL
  !
  ! ***  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
  !
  ! ***  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
  !
  !----------------------------------------------------------------------C
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY                 DATE APPROVED    BY
  !
  !----------------------------------------------------------------------C
  !
  ! *** *******************************************************************C
  !

  implicit none

  integer :: INEG
  real   :: FPROBDEP,YVAL,XVAL,POLYX,EXPY,FUNY,TMPVAL
  real, intent(IN) :: TAUDDD,TAUBBB


  ! *** EVALUATION ASSUMES TAUBBB > TAUDDD
  YVAL = 2.04*LOG(0.25*((TAUBBB/TAUDDD)-1.)*EXP(1.27*TAUDDD))
  INEG = 0
  if( YVAL < 0.0 )then
    INEG = 1
    YVAL = ABS(YVAL)
  endif
  XVAL = 1.0/(1.0+0.3327*YVAL)
  POLYX = XVAL*(0.4632-0.1202*XVAL+0.9373*XVAL*XVAL)
  EXPY = -0.5*YVAL*YVAL
  FUNY = 0.3989*EXP(EXPY)
  TMPVAL = 1.0-FUNY*POLYX
  if( INEG == 1 ) TMPVAL = 1.0-TMPVAL
  FPROBDEP = 1.0-TMPVAL

  ! ***  NOTES
  !     0.3989 = 1/SQRT(2*PI)

END FUNCTION
