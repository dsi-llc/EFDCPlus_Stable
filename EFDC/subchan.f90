! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SUBCHAN(QCHANUT,QCHANVT,IACTIVE,DE_T)

  ! CHANGE RECORD
  ! ** SUBROUTINE SUBCHAN CALCULATES SUBGRID CHANNEL INTERACTIONS AND IS
  ! ** CALLED FROM CALPUV2C AND CALPUV9C

  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: IACTIVE(NCHANM)
  INTEGER :: NMD,LHOST,LCHNU,LCHNV,IHCHMX,IHCHMN
  
  REAL :: QCHANUT(NCHANM), QCHANVT(NCHANM)
  REAL :: HCHNMX, HCHNMN, RLAMN, RLAMO
  REAL :: DE_T, SRFHOST, SRFCHAN, WCHAN, RLCHN, HCHAN, TMPVAL

  RLAMN=QCHERR  
  RLAMO=1.-RLAMN  

  DO NMD=1,MDCHH
    CCCCHU(NMD)=0.0
    CCCCHV(NMD)=0.0
    CCCCHH(NMD)=0.0
    LHOST=LMDCHH(NMD)
    LCHNU=LMDCHU(NMD)
    LCHNV=LMDCHV(NMD)

    ! *** X-DIRECTION CHANNEL
    IF( MDCHTYP(NMD) == 1 )THEN
      IACTIVE(NMD)=0
      SRFHOST=H1P(LHOST)+BELV(LHOST)
      SRFCHAN=H1P(LCHNU)+BELV(LCHNU)
      IF( SRFCHAN > SRFHOST )THEN
        IF( H1P(LCHNU) > HDRY )THEN
          IACTIVE(NMD)=1
        ENDIF
      ENDIF
      IF( SRFHOST > SRFCHAN )THEN
        IF( H1P(LHOST) > HDRY )THEN
          IACTIVE(NMD)=1
        ENDIF
      ENDIF
      IF( HP(LHOST) <= 0.0 .OR. HP(LCHNU) <= 0.0 )THEN
        IF( IACTIVE(NMD) == 1 )THEN
          IACTIVE(NMD)=0
        ENDIF
      ENDIF
      IF( IACTIVE(NMD) == 1 )THEN
        WCHAN=DXP(LCHNU)
        RLCHN=0.5*DYP(LCHNU)+CHANLEN(NMD)
        HCHAN=0.5*DYP(LCHNU)*H1P(LCHNU)+CHANLEN(NMD)*H1P(LHOST)
        HCHAN=HCHAN/RLCHN
        IF( HCHAN > 0. )THEN
          TMPVAL=CHANFRIC(NMD)*DE_T/(HCHAN*HCHAN*WCHAN)
          CCCCHU(NMD)=1./(1.+TMPVAL*ABS(QCHANUT(NMD)))
          CCCCHV(NMD)=DE_T*HCHAN*WCHAN/RLCHN
        ENDIF
      ENDIF
    ENDIF

    ! *** Y-DIRECTION CHANNEL
    IF( MDCHTYP(NMD) == 2 )THEN
      IHCHMX=0
      IHCHMN=0
      IACTIVE(NMD)=0
      SRFHOST=H1P(LHOST)+BELV(LHOST)
      SRFCHAN=H1P(LCHNV)+BELV(LCHNV)
      IF( SRFCHAN > SRFHOST )THEN
        IF( H1P(LCHNV) > HDRY )THEN
          HCHNMX=-H1P(LCHNV)
          IHCHMX=1
          IACTIVE(NMD)=1
        ENDIF
      ENDIF
      IF( SRFHOST > SRFCHAN )THEN
        IF( H1P(LHOST) > HDRY )THEN
          HCHNMN=H1P(LHOST)
          IHCHMN=1
          IACTIVE(NMD)=1
        ENDIF
      ENDIF
      IF( HP(LHOST) <= 0.0 .OR. HP(LCHNU) <= 0.0 )THEN
        IF( IACTIVE(NMD) == 1 )THEN
          IACTIVE(NMD)=0
        ENDIF
      ENDIF
      IF( IACTIVE(NMD) == 1 )THEN
        WCHAN=DYP(LCHNV)
        RLCHN=0.5*DXP(LCHNV)+CHANLEN(NMD)
        HCHAN=0.5*DXP(LCHNV)*H1P(LCHNV)+CHANLEN(NMD)*H1P(LHOST)
        HCHAN=HCHAN/RLCHN
        IF( IHCHMX == 1 ) HCHAN=MAX(HCHAN,HCHNMX)
        IF( IHCHMN == 1 ) HCHAN=MIN(HCHAN,HCHNMN)
        IF( HCHAN > 0. )THEN
          TMPVAL=CHANFRIC(NMD)*DE_T/(HCHAN*HCHAN*WCHAN)
          CCCCHU(NMD)=1./(1.+TMPVAL*ABS(QCHANVT(NMD)))
          CCCCHV(NMD)=DE_T*HCHAN*WCHAN/RLCHN
        ENDIF
      ENDIF
    ENDIF
    
    CC(LHOST)=CC(LHOST)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
    CCCCHH(NMD)        =G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
    
    IF( MDCHTYP(NMD) == 1 )THEN
      CC(LCHNU)=CC(LCHNU)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
      TMPVAL=G*(RLAMO+RLAMN*CCCCHU(NMD))*QCHANUT(NMD) - G*RLAMN*RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST) - P1(LCHNU))
      FP(LHOST)=FP(LHOST)+TMPVAL
      FP(LCHNU)=FP(LCHNU)-TMPVAL
    ENDIF
    
    IF( MDCHTYP(NMD) == 2 )THEN
      CC(LCHNV)=CC(LCHNV)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
      TMPVAL=G*(RLAMO+RLAMN*CCCCHU(NMD))*QCHANVT(NMD) - G*RLAMN*RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST) - P1(LCHNV))
      FP(LHOST)=FP(LHOST)+TMPVAL
      FP(LCHNV)=FP(LCHNV)-TMPVAL
    ENDIF
    
  ENDDO

  RETURN
END