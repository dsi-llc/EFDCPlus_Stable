! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SUBCHAN(QCHANUT,QCHANVT,IACTIVE,DE_T)

  ! CHANGE RECORD
  ! *** SUBROUTINE SUBCHAN CALCULATES SUBGRID CHANNEL INTERACTIONS AND IS
  ! *** CALLED FROM CALPUV2C AND CALPUV9C

  use GLOBAL
  implicit none
  
  integer :: IACTIVE(NCHANM)
  integer :: NMD,LHOST,LCHNU,LCHNV,IHCHMX,IHCHMN
  
  real :: QCHANUT(NCHANM), QCHANVT(NCHANM)
  real :: HCHNMX, HCHNMN, RLAMN, RLAMO
  real :: DE_T, SRFHOST, SRFCHAN, WCHAN, RLCHN, HCHAN, TMPVAL

  RLAMN = QCHERR  
  RLAMO = 1.-RLAMN  

  do NMD = 1,MDCHH
    CCCCHU(NMD) = 0.0
    CCCCHV(NMD) = 0.0
    CCCCHH(NMD) = 0.0
    LHOST = LMDCHH(NMD)
    LCHNU = LMDCHU(NMD)
    LCHNV = LMDCHV(NMD)

    ! *** X-DIRECTION CHANNEL
    if( MDCHTYP(NMD) == 1 )then
      IACTIVE(NMD) = 0
      SRFHOST = H1P(LHOST)+BELV(LHOST)
      SRFCHAN = H1P(LCHNU)+BELV(LCHNU)
      if( SRFCHAN > SRFHOST )then
        if( H1P(LCHNU) > HDRY )then
          IACTIVE(NMD) = 1
        endif
      endif
      if( SRFHOST > SRFCHAN )then
        if( H1P(LHOST) > HDRY )then
          IACTIVE(NMD) = 1
        endif
      endif
      if( HP(LHOST) <= 0.0 .or. HP(LCHNU) <= 0.0 )then
        if( IACTIVE(NMD) == 1 )then
          IACTIVE(NMD) = 0
        endif
      endif
      if( IACTIVE(NMD) == 1 )then
        WCHAN = DXP(LCHNU)
        RLCHN = 0.5*DYP(LCHNU)+CHANLEN(NMD)
        HCHAN = 0.5*DYP(LCHNU)*H1P(LCHNU)+CHANLEN(NMD)*H1P(LHOST)
        HCHAN = HCHAN/RLCHN
        if( HCHAN > 0. )then
          TMPVAL = CHANFRIC(NMD)*DE_T/(HCHAN*HCHAN*WCHAN)
          CCCCHU(NMD) = 1./(1.+TMPVAL*ABS(QCHANUT(NMD)))
          CCCCHV(NMD) = DE_T*HCHAN*WCHAN/RLCHN
        endif
      endif
    endif

    ! *** Y-DIRECTION CHANNEL
    if( MDCHTYP(NMD) == 2 )then
      IHCHMX = 0
      IHCHMN = 0
      IACTIVE(NMD) = 0
      SRFHOST = H1P(LHOST)+BELV(LHOST)
      SRFCHAN = H1P(LCHNV)+BELV(LCHNV)
      if( SRFCHAN > SRFHOST )then
        if( H1P(LCHNV) > HDRY )then
          HCHNMX = -H1P(LCHNV)
          IHCHMX = 1
          IACTIVE(NMD) = 1
        endif
      endif
      if( SRFHOST > SRFCHAN )then
        if( H1P(LHOST) > HDRY )then
          HCHNMN = H1P(LHOST)
          IHCHMN = 1
          IACTIVE(NMD) = 1
        endif
      endif
      if( HP(LHOST) <= 0.0 .or. HP(LCHNU) <= 0.0 )then
        if( IACTIVE(NMD) == 1 )then
          IACTIVE(NMD) = 0
        endif
      endif
      if( IACTIVE(NMD) == 1 )then
        WCHAN = DYP(LCHNV)
        RLCHN = 0.5*DXP(LCHNV)+CHANLEN(NMD)
        HCHAN = 0.5*DXP(LCHNV)*H1P(LCHNV)+CHANLEN(NMD)*H1P(LHOST)
        HCHAN = HCHAN/RLCHN
        if( IHCHMX == 1 ) HCHAN = MAX(HCHAN,HCHNMX)
        if( IHCHMN == 1 ) HCHAN = MIN(HCHAN,HCHNMN)
        if( HCHAN > 0. )then
          TMPVAL = CHANFRIC(NMD)*DE_T/(HCHAN*HCHAN*WCHAN)
          CCCCHU(NMD) = 1./(1.+TMPVAL*ABS(QCHANVT(NMD)))
          CCCCHV(NMD) = DE_T*HCHAN*WCHAN/RLCHN
        endif
      endif
    endif
    
    CC(LHOST) = CC(LHOST)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
    CCCCHH(NMD)         = G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
    
    if( MDCHTYP(NMD) == 1 )then
      CC(LCHNU) = CC(LCHNU)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
      TMPVAL = G*(RLAMO+RLAMN*CCCCHU(NMD))*QCHANUT(NMD) - G*RLAMN*RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST) - P1(LCHNU))
      FP(LHOST) = FP(LHOST)+TMPVAL
      FP(LCHNU) = FP(LCHNU)-TMPVAL
    endif
    
    if( MDCHTYP(NMD) == 2 )then
      CC(LCHNV) = CC(LCHNV)+G*RLAMN*RLAMN*CCCCHU(NMD)*CCCCHV(NMD)
      TMPVAL = G*(RLAMO+RLAMN*CCCCHU(NMD))*QCHANVT(NMD) - G*RLAMN*RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST) - P1(LCHNV))
      FP(LHOST) = FP(LHOST)+TMPVAL
      FP(LCHNV) = FP(LCHNV)-TMPVAL
    endif
    
  enddo

  return
END
