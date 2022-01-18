! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SETOPENBC2

  ! CHANGE RECORD
  ! ** SUBROUTINE SETOBC SETS OPEN BOUNDARY CONDITIONS FOR
  !    CALPUV2T & CALPUV2C   AND  CALPUV9 & CALPUV9C
  !
  ! *** MODIFIED BY PAUL M. CRAIG

  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: LL,L
  
  ! *** SET OPEN BOUNDARY SURFACE ELEVATIONS FOR ADJACENT CELLS (OPTIONS 0 AND 3)
  ! *** IF NOT RADIATION SEPARATION SPECIFIED

  DO LL=1,NPBW
    IF( ISPBW(LL) == 0 .OR. ISPBW(LL) == 3 )THEN
      L=LPBW(LL)
      CW(LEC(L))=0.
    ENDIF
  ENDDO

  DO LL=1,NPBE
    IF( ISPBE(LL) == 0 .OR. ISPBE(LL) == 3 )THEN
      L=LPBE(LL)
      CE(LWC(L))=0.
    ENDIF
  ENDDO

  DO LL=1,NPBS
    IF( ISPBS(LL) == 0 .OR. ISPBS(LL) == 3 )THEN
      L=LPBS(LL)
      CS(LNC(L))=0.
    ENDIF
  ENDDO

  DO LL=1,NPBN
    IF( ISPBN(LL) == 0 .OR. ISPBN(LL) == 3 )THEN
      L=LPBN(LL)
      CN(LSC(L))=0.
    ENDIF
  ENDDO

  RETURN
END
