! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE BANKEROSED

  ! ***  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
  !
  ! ***  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
  !
  !----------------------------------------------------------------------C
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY                 DATE APPROVED    BY
  !----------------------------------------------------------------------C
  !
  ! *** *******************************************************************C
  !
  ! ***  SUBROUTINE BANKEROSED CALCULATES SEDIMENT TRANSPORT DUE TO BANK
  ! ***  EROSION.  TRANSPORT IS FROM BANK BED CELL TO CHANNEL BED AND
  ! ***  WATER COLUMN CELLS
  !
  ! *** *******************************************************************C

  use GLOBAL

  implicit none

  integer :: L,NS,NT,LBANK,LCHAN,K,NX,NP
  real    :: TIME,BKEROBKB,BKEROCHW,BKEROCHB,WVEL,TMPVAL

  ! *** *******************************************************************C
  !
  if( ISDYNSTP == 0 )then
    TIME = (DT*FLOAT(N)+TCON*TBEGIN)/TCON
  else
    TIME = TIMESEC/TCON
  endif

  ! *** *******************************************************************C
  !
  !  INITIALIZE BANK EROSION VARIABLES
  !
  do L = 1,LC
    QSBDTOPBEBKB(L) = 0.
    QSBDTOPBECHB(L) = 0.
    QSBDTOPBECHW(L) = 0.
    QWBDTOPBEBKB(L) = 0.
    QWBDTOPBECHB(L) = 0.
    QWBDTOPBECHW(L) = 0.
  enddo
  !
  do NS = 1,NSED
  do L = 1,LC
    SEDFBEBKB(L,NS) = 0.
    SEDFBECHB(L,NS) = 0.
    SEDFBECHW(L,NS) = 0.
  enddo
  enddo
  !
  do NS = 1,NSND
  do L = 1,LC
    SNDFBEBKB(L,NS) = 0.
    SNDFBECHB(L,NS) = 0.
    SNDFBECHW(L,NS) = 0.
  enddo
  enddo
  !
  do NT = 1,NTOX
  do L = 1,LC
    TOXFBEBKB(L,NT) = 0.
    TOXFBECHB(L,NT) = 0.
    TOXFBECHW(L,NT) = 0.
  enddo
  enddo
  !
  ! *** *******************************************************************C
  !
  !  LOAD SEDIMENT FLUXES
  !
  if( ISTRAN(6) > 0 )then
    do NS = 1,NSED
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        BKEROBKB = VFRBED(LBANK,K,NS)*FBESER(NP) &
                 *BESERT(NBESERN(NP))
        BKEROCHW = FWCBESERT(NBESERN(NP))
        BKEROCHB = 1.-FWCBESERT(NBESERN(NP))
  ! *** *******************************************************************C
  ! HQI Change to restrict bank erosion to cells with available sediment
  ! RM, 01/08/07
  !            SEDFBEBKB(NP,NS) = BKEROBKB*DXYIP(LBANK)
  !            SEDFBECHB(NP,NS) = -BKEROCHB*BKEROBKB*DXYIP(LCHAN)
  !            SEDFBECHW(NP,NS) = BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        if( (DELT*BKEROBKB*DXYIP(LBANK)) > SEDB(LBANK,KBT(LBANK),NS) )then

          SEDFBEBKB(NP,NS) = 0.
          SEDFBECHB(NP,NS) = 0.
          SEDFBECHW(NP,NS) = 0.
        else
          SEDFBEBKB(NP,NS) = BKEROBKB*DXYIP(LBANK)
          SEDFBECHB(NP,NS) = -BKEROCHB*BKEROBKB*DXYIP(LCHAN)
          SEDFBECHW(NP,NS) = BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        endif
  ! End HQI Change
  ! *** *******************************************************************C
      enddo
    enddo
  endif
  !
  if( ISTRAN(7) > 0 )then
    do NX = 1,NSND
    NS = NSED+NX
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        BKEROBKB = VFRBED(LBANK,K,NS)*FBESER(NP) &
                 *BESERT(NBESERN(NP))
        BKEROCHW = FWCBESERT(NBESERN(NP))
        BKEROCHB = 1.-FWCBESERT(NBESERN(NP))
  ! *** *******************************************************************C
  ! HQI Change to restrict bank erosion to cells with available sediment
  ! RM, 01/08/07
  !     SNDFBEBKB(NP,NX) = BKEROBKB*DXYIP(LBANK)
  !     SNDFBECHB(NP,NX) = -BKEROCHB*BKEROBKB*DXYIP(LCHAN)
  !     SNDFBECHW(NP,NX) = BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        if( (DELT*BKEROBKB*DXYIP(LBANK)) > SNDB(LBANK,KBT(LBANK),NX) )then

          SNDFBEBKB(NP,NX) = 0.
          SNDFBECHB(NP,NX) = 0.
          SNDFBECHW(NP,NX) = 0.
        else
          SNDFBEBKB(NP,NX) = BKEROBKB*DXYIP(LBANK)
          SNDFBECHB(NP,NX) = -BKEROCHB*BKEROBKB*DXYIP(LCHAN)
          SNDFBECHW(NP,NX) = BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        endif
  ! End HQI Change
  ! *** *******************************************************************C
      enddo
    enddo
  endif
  !
  ! *** *******************************************************************C
  !
  !  UPDATE BED AND WATER COLUMN SEDIMENT CONCENTRATION
  !
  if( ISTRAN(6) > 0 )then
    do NS = 1,NSED
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        WVEL = DELT*HPI(LCHAN)*DZIC(LCHAN,KSZ(LCHAN))
        SEDB(LBANK,KBT(LBANK),NS) = SEDB(LBANK,KBT(LBANK),NS) &
                                 -DELT*SEDFBEBKB(NP,NS)
        SEDB(LCHAN,KBT(LCHAN),NS) = SEDB(LCHAN,KBT(LCHAN),NS) &
                                 -DELT*SEDFBECHB(NP,NS)
        SED(LCHAN,KSZ(LCHAN),NS) = SED(LCHAN,1,NS)+WVEL*SEDFBECHW(NP,NS)
      enddo
    enddo
  endif
  !
  if( ISTRAN(7) > 0 )then
    do NS = 1,NSND
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        WVEL = DELT*HPI(LCHAN)*DZIC(LCHAN,KSZ(LCHAN))
        SNDB(LBANK,KBT(LBANK),NS) = SNDB(LBANK,KBT(LBANK),NS) &
                                 -DELT*SNDFBEBKB(NP,NS)
        SNDB(LCHAN,KBT(LCHAN),NS) = SNDB(LCHAN,KBT(LCHAN),NS) &
                                 -DELT*SNDFBECHB(NP,NS)
        SND(LCHAN,KSZ(LCHAN),NS) = SND(LCHAN,1,NS)+WVEL*SNDFBECHW(NP,NS)
      enddo
    enddo
  endif
  !
  ! *** *******************************************************************C
  !
  !  CALCULATE SEDIMENT AND WATER VOLUME FLUXES FROM BANK
  !
  !  COHESIVE
  !
  if( ISTRAN(6) > 0 )then
  if( IBMECH == 1 .and. SEDVRDT < 0.0 )then
  !
    do NS = 1,NSED
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        QSBDTOPBEBKB(NP) = QSBDTOPBEBKB(NP) + 0.001*SEDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
        QWBDTOPBEBKB(NP) = QWBDTOPBEBKB(NP) + 0.001*VDRBED(LBANK,K)*SEDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
    enddo
    enddo
  !
  else
  !
    do NS = 1,NSED
      DSEDGMM = 1./(1.E6*SSG(NS))
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        QSBDTOPBEBKB(NP) = QSBDTOPBEBKB(NP) + DSEDGMM*SEDFBEBKB(NP,NS)
        QWBDTOPBEBKB(NP) = QWBDTOPBEBKB(NP) + DSEDGMM*VDRBED(LBANK,K)*SEDFBEBKB(NP,NS)
      enddo
    enddo
  !
  endif
  endif
  !
  !  NONCOHESIVE
  !
  if( ISTRAN(7) > 0 )then
  if( IBMECH == 1 .and. SEDVRDT < 0.0 )then
  !
    do NS = 1,NSND
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        QSBDTOPBEBKB(NP) = QSBDTOPBEBKB(NP) + 0.001*SNDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
        QWBDTOPBEBKB(NP) = QWBDTOPBEBKB(NP) + 0.001*VDRBED(LBANK,K)*SNDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
      enddo
    enddo
  !
  else
  !
    do NS = 1,NSND
      DSEDGMM = 1./(1.E6*SSG(NS+NSED))
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K = KBT(LBANK)
        QSBDTOPBEBKB(NP) = QSBDTOPBEBKB(NP) + DSEDGMM*SNDFBEBKB(NP,NS)
        QWBDTOPBEBKB(NP) = QWBDTOPBEBKB(NP) + DSEDGMM*VDRBED(LBANK,K)*SNDFBEBKB(NP,NS)
      enddo
    enddo
  !
  endif
  endif
  !
  ! *** *******************************************************************C
  !
  !  CALCULATE SEDIMENT AND WATER VOLUME FLUXES FOR CHANNEL BED AND
  !  WATER COLUMN
  !
  do NP = 1,NBEPAIR
    LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
    LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
    BKEROCHW = FWCBESERT(NBESERN(NP))
    BKEROCHB = 1.-FWCBESERT(NBESERN(NP))
    TMPVAL = DXYP(LBANK)*DXYIP(LCHAN)
    QSBDTOPBECHB(NP) = -TMPVAL*BKEROCHB*QSBDTOPBEBKB(NP)
    QWBDTOPBECHB(NP) = -TMPVAL*BKEROCHB*QWBDTOPBEBKB(NP)
    QSBDTOPBECHW(NP) = TMPVAL*BKEROCHW*QSBDTOPBEBKB(NP)
    QWBDTOPBECHW(NP) = TMPVAL*BKEROCHW*QWBDTOPBEBKB(NP)
  enddo
  !
  ! *** *******************************************************************C
  !
  return
END
