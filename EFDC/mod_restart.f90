! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE RESTART_MODULE
! ** THIS MODULE INCLUDES THE SUBROUTINES FOR
! ** WRITTING AND READING RESTART FILES:
! ** Restart_Out, WSMRST, WQ_WCRST_OUT, RESTMOD
! ** Restart_In, RESTIN2, RESTIN10, WQSDRST_OUT, WQ_WCRST_IN
! ** ALL THE BINARY OUTPUT FILES WILL BE APPENDED
! ** UPDATED: 2013-06-23 BY DH CHUNG
USE GLOBAL
Use WQ_DIAGENESIS
Use Allocate_Initialize      
USE Variables_WQ, only:WQCONSTIT
USE JULIANMOD
USE IFPORT
USE INFOMOD,ONLY:NUMCOL,SKIPCOM,READSTR
Use MPI
Use Variables_MPI
Use Mod_DSI_SendRecv
Use Broadcast_Routines
Use Variables_MPI_Write_Out
Use Variables_MPI_Mapping

IMPLICIT NONE

CONTAINS

SUBROUTINE Restart_Out(IRSTYP)

  ! ** GENERATE A RESTART FILE FOR USE BY EFDC+

  INTEGER,INTENT(IN) :: IRSTYP
  INTEGER :: I, IDUM, JDUM, K, L, LG, LL, M, MD, NMD, NP, NS, NT, NX, NHR, VER, NCTL, IMASK
  INTEGER :: NSCM0

  REAL :: RVAL, TMP1, TMP2, TMP3, TMP4, SURF, HB1, EL1
  CHARACTER(2) :: SNHR

  ! *** Build the global arrays needed for the RESTART.OUT file
  CALL Gather_Restart_Arrays

  ! *** Only the master process is needed after this point
  if( process_id /= master_id ) return

  IF( ISGREGOR == 1 ) CALL TOGREGOR(SDATE,TIMEDAY)

  IF( IRSTYP == 0 )THEN
    PRINT *,'Restart Snapshot @ Timeday: ',SNGL(TIMEDAY)
    IF( ISGREGOR == 1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = 'RESTART'//'_'//TRIM(SDATE)//'_'//SNHR//'.OUT'
    ELSE
      RESFILE = 'RESTART.OUT'
    ENDIF
    
    OPEN(99,FILE=OUTDIR//TRIM(RESFILE),STATUS='REPLACE')
  ENDIF
  IF( IRSTYP == 1 )THEN
    OPEN(99,FILE=OUTDIR//'CRASHST.OUT',STATUS='REPLACE')
  ENDIF

  ! *****************************************************************************************************
  ! *** BEGIN RESTART.OUT
  VER = 1020                                      ! *** Changed for Version 10.2
  WRITE(99,'(I20,4X,F12.4,2I10)') N, TIMEDAY, VER
  WRITE(99,'(9I6)') ISCO(0:8)
  DO LG=2,LA_Global
    WRITE(99,906) HP_Global(LG), H1P_Global(LG), HWQ_Global(LG), H2WQ_Global(LG), BELV_Global(LG)
    WRITE(99,906) UHDYE_Global(LG), UHDY1E_Global(LG), VHDXE_Global(LG), VHDX1E_Global(LG)
    WRITE(99,913) ISCDRY_Global(LG), NATDRY_Global(LG), IDRY_Global(LG), SUB_Global(LG), SVB_Global(LG)
    WRITE(99,906) (U_Global(LG,K),K=1,KC)
    WRITE(99,906) (U1_Global(LG,K),K=1,KC)
    WRITE(99,906) (V_Global(LG,K),K=1,KC)
    WRITE(99,906) (V1_Global(LG,K),K=1,KC)
    WRITE(99,906) (QQ_Global(LG,K),K=0,KC)
    WRITE(99,906) (QQ1_Global(LG,K),K=0,KC)
    WRITE(99,906) (QQL_Global(LG,K),K=0,KC)
    WRITE(99,906) (QQL1_Global(LG,K),K=0,KC)
    WRITE(99,906) (DML_Global(LG,K),K=0,KC)
    IF( ISTRAN(2) >= 1 .AND. ISICE > 2 )THEN
      WRITE(99,906) ICETHICK_Global(LG), ICETEMP_Global(LG)
    ENDIF
    IF( ISTRAN(1) >= 1 .AND. ISCO(1) == 1 )THEN
      WRITE(99,906) (SAL_Global(LG,K),K=1,KC)
      WRITE(99,906) (SAL1_Global(LG,K),K=1,KC)
    ENDIF
    IF( ISTRAN(2) >= 1 .AND. ISCO(2) == 1 )THEN
      WRITE(99,906) TEMB_Global(LG),(MAX(TEM_Global(LG,K),0.0),K=1,KC)
      WRITE(99,906) (MAX(TEM1_Global(LG,K),0.0),K=1,KC)
    ENDIF
    IF( ISTRAN(3) >= 1 .AND. ISCO(3) == 1 )THEN
      DO MD=1,NDYE
        WRITE(99,906) (DYE_Global(LG,K,MD),K=1,KC)
        WRITE(99,906) (DYE1_Global(LG,K,MD),K=1,KC)
      ENDDO
    ENDIF
    IF( ISTRAN(4) >= 1 .AND. ISCO(4) == 1 )THEN
      WRITE(99,906) SFLSBOT(LG),(SFL_Global(LG,K),K=1,KC)
      WRITE(99,906) SFLSBOT(LG),(SFL2(LG,K),K=1,KC)
    ENDIF
    IF( ISTRAN(5) >= 1 .AND. ISCO(5) == 1 )THEN
      DO NT=1,NTXM
        WRITE(99,906) (TOXB_Global(LG,K,NT),K=1,KB)
        WRITE(99,906) (TOX_Global(LG,K,NT),K=1,KC)
        WRITE(99,906) (TOXB1_Global(LG,K,NT),K=1,KB)
        WRITE(99,906) (TOX1_Global(LG,K,NT),K=1,KC)
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 .AND. ISCO(6) == 1 )THEN
      DO NS=1,NSCM
        WRITE(99,906) (SEDB_Global(LG,K,NS),K=1,KB)
        WRITE(99,906) (SED_Global(LG,K,NS),K=1,KC)
        WRITE(99,906) (SEDB1_Global(LG,K,NS),K=1,KB)
        WRITE(99,906) (SED1_Global(LG,K,NS),K=1,KC)
      ENDDO
      IF( NSCM2 > NSCM )THEN
        IF( LG == 2 )THEN
          WRITE(99,908) NSCM, NSCM2                 ! *** Make NSCM2 available before propwash read in
        ENDIF
        DO NS=NSCM+1,NSCM2
          WRITE(99,906) (SED_Global(LG,K,NS),K=1,KC)
          WRITE(99,906) (SED1_Global(LG,K,NS),K=1,KC)
        ENDDO
      ENDIF
    ENDIF
    IF( ISTRAN(7) >= 1 .AND. ISCO(7) == 1 )THEN
      DO NS=1,NSNM
        WRITE(99,906) (SNDB_Global(LG,K,NS),K=1,KB)
        WRITE(99,906) (SND_Global(LG,K,NS),K=1,KC)
        WRITE(99,906) (SNDB1_Global(LG,K,NS),K=1,KB)
        WRITE(99,906) (SND1_Global(LG,K,NS),K=1,KC)
      ENDDO
    ENDIF
    IF( (ISTRAN(6) >= 1 .AND. ISCO(6) == 1 ) .OR. (ISTRAN(7) >= 1 .AND. ISCO(7) == 1 ) )THEN
      WRITE(99,906) (HBED_Global(LG,K),K=1,KB)
      WRITE(99,906) (HBED1_Global(LG,K),K=1,KB)
      WRITE(99,906) (VDRBED_Global(LG,K),K=1,KB)
      WRITE(99,906) (VDRBED1_Global(LG,K),K=1,KB)
    ENDIF
  ENDDO

  ! *** BOUNDARY CONDITIONS
  MD = 0
  NT = 0
  NS = 0
  NX = 0
  DO I = 1,7
    IF( ISTRAN(I) > 0 .AND. ISCO(I) > 0 )THEN
100   IF( I == 3 )THEN
        MD = MD + 1
        M = MSVDYE(MD)
      ELSEIF( I == 4 )THEN
        M = 3 + NDYE
      ELSEIF( I == 5 )THEN
        NT = NT + 1
        M = MSVTOX(NT)
      ELSEIF( I == 6 )THEN
        NS = NS + 1
        M = MSVSED(NS)
      ELSEIF( I == 7 )THEN
        NX = NX + 1
        M = MSVSND(NX)
      ELSE
        M = I
      ENDIF

      DO LL=1,NPBS_GL
        DO K=1,KC
          NLOS_Global(LL,K,M) = MAX( NLOS_Global(LL,K,M)-NITER, -NTSCRS_GL(LL) )
        ENDDO
        WRITE(99,908) (NLOS_Global(LL,K,M),K=1,KC)
        WRITE(99,906) (CLOS_Global(LL,K,M),K=1,KC)
      ENDDO
      DO LL=1,NPBW_GL
        DO K=1,KC
          NLOW_Global(LL,K,M) = MAX( NLOW_Global(LL,K,M)-NITER, -NTSCRW_GL(LL) )
        ENDDO
        WRITE(99,908) (NLOW_Global(LL,K,M),K=1,KC)
        WRITE(99,906) (CLOW_Global(LL,K,M),K=1,KC)
      ENDDO
      DO LL=1,NPBE_GL
        DO K=1,KC
          NLOE_Global(LL,K,M) = MAX( NLOE_Global(LL,K,M)-NITER, -NTSCRE_GL(LL) )
        ENDDO
        WRITE(99,908) (NLOE_Global(LL,K,M),K=1,KC)
        WRITE(99,906) (CLOE_Global(LL,K,M),K=1,KC)
      ENDDO
      DO LL=1,NPBN_GL
        DO K=1,KC
          NLON_Global(LL,K,M) = MAX( NLON_Global(LL,K,M)-NITER, -NTSCRN_GL(LL) )
        ENDDO
        WRITE(99,908) (NLON_Global(LL,K,M),K=1,KC)
        WRITE(99,906) (CLON_Global(LL,K,M),K=1,KC)
      ENDDO
      NSCM0 = MAX(1,NSED)                        ! *** Temporary variable since propwash also uses NSCM
      IF( I == 3 .AND. MD < NDYE ) GOTO 100      ! *** LOOP BACK AND WRITE NEXT DYE CLASS
      IF( I == 5 .AND. NT < NTXM ) GOTO 100      ! *** LOOP BACK AND WRITE NEXT TOX CLASS
      IF( I == 6 .AND. NS < NSCM0) GOTO 100      ! *** LOOP BACK AND WRITE NEXT SED CLASS
      IF( I == 7 .AND. NX < NSNM ) GOTO 100      ! *** LOOP BACK AND WRITE NEXT SND CLASS
    ENDIF
  ENDDO

  DO LG=2,LA_Global
    WRITE(99,906) QSUME_Global(LG), (QSUM_Global(LG,K),K=1,KC)
  ENDDO

  IF( MDCHH >= 1 )THEN
    DO NMD=1,MDCHH
      WRITE(99,'(6I5,2X,E17.8,2X,E17.8)') IMDCHH(NMD), JMDCHH(NMD), IMDCHU(NMD), JMDCHU(NMD), IMDCHV(NMD), JMDCHV(NMD), QCHANU(NMD), QCHANV(NMD)
    ENDDO
  ENDIF

  IF( ISGWIE >= 1 )THEN
    if( num_Processors > 1 ) Call STOPP('MPI does not support ISGWIE > 0 option')
    DO L=2,LA_Global
      WRITE(99,906) AGWELV(L), AGWELV1(L)
    ENDDO
  ENDIF

  ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE WITH A LOWER CHORD ACTIVATION TOGGLE
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP == 3 .OR. HYD_STR(NCTL).NQCTYP == 4 )THEN
      if( num_Processors > 1 ) Call STOPP('MPI does not support NQCTYP = 3 or 4 option')
      ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE (4) DEPENDANT FLOWS
      WRITE(99,906) SAVESUB(1, NCTL), SAVESVB(1, NCTL), SAVESUB(2, NCTL), SAVESVB(2, NCTL)
      WRITE(99,906) LOWCHORDU(NCTL), LOWCHORDV(NCTL)
      WRITE(99,908) NCTL, NLOWCHORD(NCTL)
    ENDIF
  ENDDO
906 FORMAT(100E17.8)
908 FORMAT(100I10)
912 FORMAT(3I5,12F7.3,5(15X,12F7.3))
913 FORMAT(3I10,2F7.3)

  CLOSE(99)

  ! *** SEDZLJ
  IF( LSEDZLJ  )THEN
    IF( ISGREGOR == 1 )THEN
      RSTFILE = 'SEDBED_HOT'//'_'//TRIM(SDATE)//'_'//SNHR//'.OUT'
    ELSE
      RSTFILE = 'SEDBED_HOT.OUT'
    ENDIF
    OPEN(99,FILE=OUTDIR//TRIM(RSTFILE),FORM='FORMATTED',STATUS='REPLACE')

    WRITE(99,34569) ((LAYERACTIVE_Global(K,LG),K=1,KB),LG=2,LA_Global)
    WRITE(99,34569) (KBT_Global(LG),LG=2,LA_Global)
    WRITE(99,34567) (D50AVG_Global(LG),LG=2,LA_Global)
    WRITE(99,34568) ((BULKDENS_Global(K,LG),K=1,KB),LG=2,LA_Global)
    WRITE(99,34568) ((TSED_Global(K,LG),K=1,KB),LG=2,LA_Global)
    WRITE(99,34568) (((PERSED_Global(NS,K,LG),NS=1,NSCM),K=1,KB),LG=2,LA_Global)

    IF( ICALC_BL > 0 )THEN
      WRITE(99,34568) ((CBL_Global(LG,NS),NS=1,NSCM),LG=2,LA_Global)
      IF( ISTRAN(5) > 0 ) WRITE(99,34568) ((CBLTOX_Global(LG,NT),NT=1,NTOX),LG=2,LA_Global)
    ENDIF

    CLOSE(99)
34567 FORMAT(E17.9)
34568 FORMAT(6E17.9)
34569 FORMAT(8I8)
  END IF

  ! *** LPT RESTART
  IF( ISPD > 0 )THEN
    IF( ISGREGOR == 1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = 'DRIFTER'//'_'//TRIM(SDATE)//'_'//SNHR//'.RST'
    ELSE
      RESFILE = 'DRIFTER.RST'
    ENDIF

    OPEN(99,FILE=OUTDIR//TRIM(RESFILE),STATUS='REPLACE')

    WRITE(99,'(50I2)') JSPD(1:NPD)
    WRITE(99,'(30I3)') KLA(1:NPD)
    WRITE(99,'(2(I8,2F15.5,F10.5))') (LLA(NP), XLA(NP), YLA(NP), ZLA(NP),NP=1,NPD)

    CLOSE(99,STATUS='KEEP')
  ENDIF

  ! *** SPECIAL FILES
  !IF( ISWAVE >= 1 )THEN               THIS FILE WAS NOT READ IN Restart_In
  !  OPEN(1,FILE=OUTDIR//'WVQWCP.OUT',STATUS='UNKNOWN')
  !  CLOSE(1, STATUS='DELETE')
  !  OPEN(1,FILE=OUTDIR//'WVQWCP.OUT',STATUS='UNKNOWN')
  !  DO L=2,LA_Global
  !    WRITE(1,'(2I5,2X,6E13.4)') IL(L), JL(L), QQWV1(L), QQWV2(L), QQWV3(L), QQWC(L), QQWCR(L), QQ_Global(L,0)
  !  ENDDO
  !  CLOSE(1)
  !ENDIF

  ! **  OUTPUT SALINITY AND TEMPERATURE DATA ASSIMILATION
  IF( NLCDA > 0 )THEN
    ! *** DEPRECATED
  ENDIF

  IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. ISDTXBUG == 1 .AND. TIMEDAY >= TIMEEND )THEN
    OPEN(1,FILE=OUTDIR//'BEDRST.SED',STATUS='REPLACE')
    WRITE(1,111)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (SEDB_Global(LG,K,1),K=1,KB)
      IF( NSCM > 1 )THEN
        DO NX=2,NSCM
          WRITE(1,102)(SEDB_Global(LG,K,NX),K=1,KB)
        END DO
      ENDIF
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.SND',STATUS='REPLACE')
    WRITE(1,112)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (SNDB_Global(LG,K,1),K=1,KB)
      IF( NSNM > 1 )THEN
        DO NX=2,NSNM
          WRITE(1,102)(SNDB_Global(LG,K,NX),K=1,KB)
        END DO
      ENDIF
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.VDR',STATUS='REPLACE')
    WRITE(1,113)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (VDRBED_Global(LG,K),K=1,KB)
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.POR',STATUS='REPLACE')
    WRITE(1,114)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (PORBED_Global(LG,K),K=1,KB)
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.ZHB',STATUS='REPLACE')
    WRITE(1,115)
    DO LG=2,LA_Global
      HB1 = 0.0
      DO K=1,KBT(LG)
        HB1 = HB1 + HBED_Global(LG,K)       ! *** HBEDA(LG)
      END DO
      EL1 = BELV_Global(LG) - HB1           ! *** ZELBEDA(LG)

      WRITE(1,101) IL_GL(LG), JL_GL(LG), EL1, HB1, (HBED_Global(LG,K),K=1,KB)
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.BDN',STATUS='REPLACE')
    WRITE(1,116)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (BDENBED_Global(LG,K),K=1,KB)
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.ELV',STATUS='REPLACE')
    WRITE(1,117)
    RVAL = 0.
    TMP1 = 0.
    TMP2 = 0.
    TMP3 = 0.
    TMP4 = 0.
    DO LG=2,LA_Global
      HB1 = 0.0
      DO K=1,KBT(LG)
        HB1 = HB1 + HBED_Global(LG,K)       ! *** HBEDA(LG)
      END DO
      EL1 = BELV_Global(LG) - HB1           ! *** ZELBEDA(LG)

      RVAL = RVAL + 1.
      TMP1 = TMP1 + EL1
      TMP2 = TMP2 + HB1
      TMP3 = TMP3 + BELV_Global(LG)
      TMP4 = TMP4 + HP_Global(LG)
      SURF = HP_Global(LG) + BELV_Global(LG)
      WRITE(1,101) IL_GL(LG), JL_GL(LG), EL1, HB1, BELV_Global(LG), HP_Global(LG), SURF
    ENDDO
    TMP1 = TMP1/RVAL
    TMP2 = TMP2/RVAL
    TMP3 = TMP3/RVAL
    TMP4 = TMP4/RVAL
    IDUM = 0
    JDUM = 0
    WRITE(1,101) IDUM, JDUM, TMP1, TMP2, TMP3, TMP4
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'WATRST.SED',STATUS='REPLACE')
    WRITE(1,118)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (SED_Global(LG,K,1),K=1,KC)
      IF( NSCM > 1 )THEN
        DO NX=2,NSCM
          WRITE(1,102)(SED_Global(LG,K,NX),K=1,KC)
        END DO
      ENDIF
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'WATRST.SND',STATUS='REPLACE')
    WRITE(1,119)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), (SND_Global(LG,K,1),K=1,KC)
      IF( NSNM > 1 )THEN
        DO NX=2,NSNM
          WRITE(1,102)(SND_Global(LG,K,NX),K=1,KC)
        END DO
      ENDIF
    ENDDO
    CLOSE(1)

    OPEN(1,FILE=OUTDIR//'BEDRST.BDL',STATUS='REPLACE')
    WRITE(1,120)
    DO LG=2,LA_Global
      WRITE(1,101) IL_GL(LG), JL_GL(LG), QSBDLDX_Global(LG,1),QSBDLDY_Global(LG,1)
      IF( NSNM > 1 )THEN
        DO NX=2,NSNM
          WRITE(1,102)QSBDLDX_Global(LG,NX),QSBDLDY_Global(LG,NX)
        END DO
      ENDIF
    ENDDO
    CLOSE(1)

    IF( ISTRAN(5) > 0 )THEN
      OPEN(1,FILE=OUTDIR//'BEDRST.TOX',STATUS='REPLACE')
      DO NT=1,NTXM
        WRITE(1,121)NT
        DO LG=2,LA_Global
          WRITE(1,101) IL_GL(LG), JL_GL(LG), (TOXB_Global(LG,K,NT),K=1,KB)
        ENDDO
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF
101 FORMAT(2I5,50E13.5)
102 FORMAT(10X,50E13.5)
111 FORMAT('   IL   JL    SEDBT(K=1,KB)')
112 FORMAT('   IL   JL    SNDBT(K=1,KB)')
113 FORMAT('   IL   JL    VRDBED(K=1,KB)')
114 FORMAT('   IL   JL    PORBED(K=1,KB)')
115 FORMAT('   IL   JL    ZBEDB        HBEDT        HBED(K=1,KB)')
116 FORMAT('   IL   JL    BDENBED(K=1,KB)')
117 FORMAT('   IL   JL    ZBEDB        HBEDT        BELV', '        HWCOL        SELV')
118 FORMAT('   IL   JL    SEDT(K=1,KC)')
119 FORMAT('   IL   JL    SNDT(K=1,KC)')
120 FORMAT('   IL   JL    QSBDLDX      QSBDLDY')
121 FORMAT('   IL   JL    TOXB(K=1,KB,NT)  NT = ',I5)

END SUBROUTINE Restart_Out


SUBROUTINE Restart_In(OPT)
  ! ** ADDED CODE TO PROPERLY INITIAL RESTART INPUT FOR DRYING AND WETTING
  ! ** SUBROUTINE Restart_In READS A RESTART FILE EXPORTED BY EFDC_DSI
  ! ** USING TRUE RESTART FILE IN EFDC_DSI (ISRESTI == 1/-1/11)
  ! ** ISRESTI=1/-1,Restart_In_Ver=720: RESTART FILES (.INP) BY EFDC_7.2
  ! **              Restart_In_Ver=710: RESTART FILES (.INP) BY EFDC_7.0-7.1

  INTEGER,INTENT(IN),OPTIONAL :: OPT
  INTEGER :: ISBELVC, I, K, M, MD, NS, NT, NX, NMD, ITMP1, JTMP1, NCOL, NCTL, ISCOCHK(0:10), IMASK
  INTEGER :: ITMP2, JTMP2, ITMP3, JTMP3, LS, LN, LDUM, IDUM, JDUM, IDFLAG, ICORDRY, IU, JU
  INTEGER :: III, JJJ, L, LL, LG, LU, ID, JD, LD

  REAL :: BELTMP, TMPVAL, RDZC, HDRY2, DHPDT, SUBW, SUBE, SVBS, SVBN, RDRY
  REAL :: SUBL, SVBL

  REAL,ALLOCATABLE,DIMENSION(:) :: TDUMMY
  CHARACTER(100) :: STR

  ALLOCATE(TDUMMY(KCM))
  TDUMMY=0.

  RSTFIRST_WS  = 0
  RSTFIRST_VEL = 0
  RSTFIRST_WC  = 0
  RSTFIRST_WQ  = 0
  RSTFIRST_SEDZLJ = 0
  RSTFIRST_RPEM = 0
  RSTFIRST_BED  = 0
  RSTFIRST_SD   = 0
  RSTFIRST_BC   = 0

  if( process_id == master_id )then
    WRITE(*,'(A)')'READING RESTART FILE: RESTART.INP'
    OPEN(UINP,FILE='restart.inp',STATUS='UNKNOWN')
    
    ISBELVC=0

    READ(UINP,'(A)',ERR=1000) STR
    NCOL = NUMCOL(STR)
    IF( NCOL == 3 )THEN
      READ(STR,*,ERR=1000) NRESTART, TBEGINC, Restart_In_Ver
    ELSE
      READ(STR,*,ERR=1000) NRESTART, TBEGINC
      Restart_In_Ver = 710
    ENDIF

    IF( PRESENT(OPT) )THEN
      CLOSE(UINP)
      RETURN
    ENDIF

    IF( Restart_In_Ver <= 710 .AND. ISTRAN(2) > 0 .AND. ISCI(2) >= 1 )THEN
      WRITE(*,'(A)')'READING RESTART FILE: TEMP.RST'
      OPEN(111,FILE='temp.rst',STATUS='UNKNOWN')
      
      DO L=2,LA_Global
        READ(111,*) LDUM, IDUM, JDUM, (TDUMMY(K),K=1,KC),TEMB_Global(L)
      ENDDO
      CLOSE(111)

    ELSEIF( Restart_In_Ver >= 720 )THEN
      READ(UINP,*) ISCOCHK(0:8)   ! *** Read flags indicating what is in the file
    ENDIF

    ! *****************************************************************************************************
    DO LG=2,LA_Global
      IF( ISRESTIOPT == 0 )THEN
        ! *** Standard format
        READ(UINP,*,ERR=1001) HP_Global(LG), H1P_Global(LG), HWQ_Global(LG), H2WQ_Global(LG), BELV_Global(LG)
      ELSE
        ! *** Allow the bathymetry to use the values in BELV_Global.INP instead of the restart file
        READ(UINP,*,ERR=1002)HP_Global(LG), H1P_Global(LG), HWQ_Global(LG), H2WQ_Global(LG), BELTMP
        IF( BELTMP /= BELV_Global(LG) )THEN
          ISBELVC = 1
          WRITE(6,600) IL_GL(LG), IL_GL(LG), BELTMP, BELV_Global(LG)
          WRITE(8,600) IL_GL(LG), JL_GL(LG), BELTMP, BELV_Global(LG)
          HP_Global(LG)   = HP_Global(LG)   + BELTMP - BELV_Global(LG)
          H1P_Global(LG)  = H1P_Global(LG)  + BELTMP - BELV_Global(LG)
          HWQ_Global(LG)  = HWQ_Global(LG)  + BELTMP - BELV_Global(LG)
          H2WQ_Global(LG) = H2WQ_Global(LG) + BELTMP - BELV_Global(LG)
        ENDIF
      ENDIF
      IF( HP_Global(LG) < 0.0 .OR. H1P_Global(LG) < 0.0 )THEN
        WRITE(6,9696) LG, IL_GL(LG), JL_GL(LG), HP_Global(LG), H1P_Global(LG)
        CALL STOPP('.')
      ELSEIF( HP_Global(LG) < 0.0001 .OR. H1P_Global(LG) < 0.0001 )THEN
        WRITE(6,9698) LG, IL(LG), JL(LG), HP_Global(LG), H1P_Global(LG)
        HP_Global(LG)    = HDRY*0.9
        H1P_Global(LG) = HDRY*0.9
      ENDIF
      
      READ(UINP,*,ERR=1003) UHDYE_Global(LG), UHDY1E_Global(LG), VHDXE_Global(LG), VHDX1E_Global(LG)

      IF( Restart_In_Ver > 1000 )THEN
        READ(UINP,*) ISCDRY_Global(LG), NATDRY_Global(LG), IDRY_Global(LG), SUB_Global(LG), SVB_Global(LG)
      ENDIF

      READ(UINP,*,ERR=1004) (U_Global(LG,K),K=1,KC)
      READ(UINP,*,ERR=1005) (U1_Global(LG,K),K=1,KC)
      READ(UINP,*,ERR=1006) (V_Global(LG,K),K=1,KC)
      READ(UINP,*,ERR=1007) (V1_Global(LG,K),K=1,KC)
      READ(UINP,*,ERR=1008) (QQ_Global(LG,K),K=0,KC)
      READ(UINP,*,ERR=1009) (QQ1_Global(LG,K),K=0,KC)
      READ(UINP,*,ERR=1010) (QQL_Global(LG,K),K=0,KC)
      READ(UINP,*,ERR=1011) (QQL1_Global(LG,K),K=0,KC)
      READ(UINP,*,ERR=1012) (DML_Global(LG,K),K=0,KC)

      IF( ISCOCHK(2) == 1 .AND. Restart_In_Ver >= 720 .AND. ISICE > 2 )THEN
        IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 )THEN
          READ(UINP,*,ERR=1035) ICETHICK_Global(LG), ICETEMP_Global(LG)                  ! ** EFDC_7.3
        ELSE
          READ(UINP,*,ERR=1035) TMPVAL, TMPVAL
        ENDIF
      ENDIF

      ! *** REMOVED THE ISTRAN CHECK SINCE IT IS NOT REQUIRED JUST TO LOAD THE RESTART FILE.  (EE7.2)
      IF( ISCOCHK(1) == 1 )THEN
        IF( ISTRAN(1) > 0 .AND. ISCI(1) > 0 )THEN
          READ(UINP,*,ERR=1013) (SAL_Global(LG,K),K=1,KC)
          READ(UINP,*,ERR=1014) (SAL1_Global(LG,K),K=1,KC)
        ELSE
          READ(UINP,*,ERR=1013) (TMPVAL,K=1,KC)
          READ(UINP,*,ERR=1014) (TMPVAL,K=1,KC)
        ENDIF
      ENDIF

      IF( ISCOCHK(2) == 1 )THEN
        IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 )THEN
          IF( Restart_In_Ver >= 720 )THEN
            READ(UINP,*,ERR=1015) TEMB_Global(LG),(TEM_Global(LG,K),K=1,KC)              ! ** EFDC_7.2
          ELSEIF (Restart_In_Ver <= 710 )THEN
            READ(UINP,*,ERR=1015) (TEM_Global(LG,K),K=1,KC)                              ! ** EFDC_7.0-7.1
          ENDIF
          READ(UINP,*,ERR=1016) (TEM1_Global(LG,K),K=1,KC)
        ELSE
          IF( Restart_In_Ver >= 720 )THEN
            READ(UINP,*,ERR=1015) TMPVAL,(TMPVAL,K=1,KC)
          ELSEIF (Restart_In_Ver <= 710 )THEN
            READ(UINP,*,ERR=1015) (TMPVAL,K=1,KC)
          ENDIF
          READ(UINP,*,ERR=1016) (TMPVAL,K=1,KC)
        ENDIF
      ENDIF

      IF( ISCOCHK(3) == 1 )THEN
        IF( ISTRAN(3) > 0 .AND. ISCI(3) > 0 )THEN
          DO MD=1,NDYE
            READ(UINP,*,ERR=1017) (DYE_Global(LG,K,MD),K=1,KC)
            READ(UINP,*,ERR=1018) (DYE1_Global(LG,K,MD),K=1,KC)
          ENDDO
        ELSE
          DO MD=1,NDYE
            READ(UINP,*,ERR=1017) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1018) (TMPVAL,K=1,KC)
          ENDDO
        ENDIF
      ENDIF

      IF( ISCOCHK(4) == 1 )THEN
        ! *** delme - not MPI
        IF( ISTRAN(4) > 0 .AND. ISCI(4) > 0 )THEN
          READ(UINP,*,ERR=1021) SFLSBOT(LG),(SFL(LG,K),K=1,KC)
          READ(UINP,*,ERR=1022) SFLSBOT(LG),(SFL2(LG,K),K=1,KC)
        ELSE
          READ(UINP,*,ERR=1021) TMPVAL,(TMPVAL,K=1,KC)
          READ(UINP,*,ERR=1022) TMPVAL,(TMPVAL,K=1,KC)
        ENDIF
      ENDIF

      IF( ISCOCHK(5) == 1 )THEN
        IF( ISTRAN(5) > 0 .AND. ISCI(5) > 0 )THEN
          DO NT=1,NTXM
            READ(UINP,*,ERR=1019) (TOXB_Global(LG,K,NT),K=1,KB)
            READ(UINP,*,ERR=1019) (TOX_Global(LG,K,NT),K=1,KC)
            READ(UINP,*,ERR=1019) (TOXB1_Global(LG,K,NT),K=1,KB)
            READ(UINP,*,ERR=1019) (TOX1_Global(LG,K,NT),K=1,KC)
          ENDDO
        ELSE
          DO NT=1,NTXM
            READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1019) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1020) (TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1020) (TMPVAL,K=1,KC)
          ENDDO
        ENDIF
      ENDIF

      IDFLAG = 0
      IF( ISCOCHK(6) == 1 )THEN
        IF( ISTRAN(6) > 0 .AND. ISCI(6) > 0 )THEN
          DO NS=1,NSCM
            READ(UINP,*,ERR=1019) (SEDB_Global(LG,K,NS),K=1,KB)
            READ(UINP,*,ERR=1019) (SED_Global(LG,K,NS),K=1,KC)
            READ(UINP,*,ERR=1019) (SEDB1_Global(LG,K,NS),K=1,KB)
            READ(UINP,*,ERR=1019) (SED1_Global(LG,K,NS),K=1,KC)
          ENDDO
          IF( NSCM2 > NSCM )THEN
            IF( LG == 2 )THEN
              READ(UINP,*,ERR=1019) NSCM, NSCM2                 ! *** Make NSCM2 available before propwash read in
            ENDIF
            DO NS=NSCM+1,NSCM2
              READ(UINP,*,ERR=1019) (SED_Global(LG,K,NS),K=1,KC)
              READ(UINP,*,ERR=1019) (SED1_Global(LG,K,NS),K=1,KC)
            ENDDO
          ENDIF
        ELSE
          DO NS=1,NSCM
            READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1019) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1020) (TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1020) (TMPVAL,K=1,KC)
          ENDDO
          IF( NSCM2 > NSCM )THEN
            IF( LG == 2 )THEN
              READ(UINP,*,ERR=1019) III, III                 ! *** Make NSCM2 available before propwash read in
            ENDIF
            DO NS=NSCM+1,NSCM2
              READ(UINP,*,ERR=1019) (TMPVAL,K=1,KC)
              READ(UINP,*,ERR=1019) (TMPVAL,K=1,KC)
            ENDDO
          ENDIF
        ENDIF
        IDFLAG = 1
      ENDIF

      IF( ISCOCHK(7) == 1 )THEN
        IF( ISTRAN(7) > 0 .AND. ISCI(7) > 0 )THEN
          DO NS=1,NSNM
            READ(UINP,*,ERR=1019) (SNDB_Global(LG,K,NS),K=1,KB)
            READ(UINP,*,ERR=1019) (SND_Global(LG,K,NS),K=1,KC)
            READ(UINP,*,ERR=1019) (SNDB1_Global(LG,K,NS),K=1,KB)
            READ(UINP,*,ERR=1019) (SND1_Global(LG,K,NS),K=1,KC)
          ENDDO
        ELSE
          DO NS=1,NSNM
            READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1019) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1020) (TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1020) (TMPVAL,K=1,KC)
          ENDDO
        ENDIF
        IDFLAG = 1
      ENDIF

      IF( IDFLAG == 1 )THEN
        IF( ( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 ) .AND. ( ISCI(6) > 0 .OR. ISCI(7) > 0 ) )THEN
          READ(UINP,*,ERR=1019) (HBED_Global(LG,K),K=1,KB)
          READ(UINP,*,ERR=1019) (HBED1_Global(LG,K),K=1,KB)
          READ(UINP,*,ERR=1019) (VDRBED_Global(LG,K),K=1,KB)
          READ(UINP,*,ERR=1019) (VDRBED1_Global(LG,K),K=1,KB)
        ELSE
          READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
          READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
          READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
          READ(UINP,*,ERR=1019) (TMPVAL,K=1,KB)
        ENDIF
      ENDIF
    ENDDO   ! *** END OF MAIN DOMAIN LOOP
    
    ! *** Check minimums
    DO LG=2,LA_Global
      DO K=1,KC
        IF( QQ_GLOBAL(LG,K) < QQMIN )THEN
          QQ_GLOBAL(LG,K) = QQMIN
        ENDIF
        IF( QQ1_GLOBAL(LG,K) < QQMIN )THEN
          QQ1_GLOBAL(LG,K) = QQMIN
        ENDIF
      ENDDO
      DO K=1,KS
        IF( QQL_GLOBAL(LG,K) < QQLMIN )THEN
          QQL_GLOBAL(LG,K) = QQLMIN
        ENDIF
        IF( DML_GLOBAL(LG,K) < DMLMIN )THEN
          DML_GLOBAL(LG,K) = DMLMIN
        ENDIF
      ENDDO
    ENDDO

  endif     ! *** End of master_id block

  ! *****************************************************************************************************
  ! *** Communicate primary global scalars
  Call Broadcast_Scalar(NRESTART,         master_id)
  Call Broadcast_Scalar(TBEGINC,          master_id)
  Call Broadcast_Scalar(Restart_In_Ver,   master_id)

  ! *** Communicate primary global arrays
  Call Broadcast_Array(HP_Global,         master_id)
  Call Broadcast_Array(H1P_Global,        master_id)
  Call Broadcast_Array(HWQ_Global,        master_id)
  Call Broadcast_Array(H2WQ_Global,       master_id)
  Call Broadcast_Array(BELV_Global,       master_id)

  Call Broadcast_Array(UHDYE_Global,      master_id)
  Call Broadcast_Array(UHDY1E_Global,     master_id)
  Call Broadcast_Array(VHDXE_Global,      master_id)
  Call Broadcast_Array(VHDX1E_Global,     master_id)
  Call Broadcast_Array(SUB_Global,        master_id)
  Call Broadcast_Array(SVB_Global,        master_id)

  Call Broadcast_Array(U_Global,          master_id)
  Call Broadcast_Array(U1_Global,         master_id)
  Call Broadcast_Array(V_Global,          master_id)
  Call Broadcast_Array(V1_Global,         master_id)
  Call Broadcast_Array(QQ_Global,         master_id)
  Call Broadcast_Array(QQ1_Global,        master_id)
  Call Broadcast_Array(QQL_Global,        master_id)
  Call Broadcast_Array(QQL1_Global,       master_id)
  Call Broadcast_Array(DML_Global,        master_id)

  ! *** Conditional Broadcasts
  Call Broadcast_Array(ISCOCHK,           master_id)
  IF( ISTRAN(1) > 0 .AND. ISCI(1) > 0 .AND. ISCOCHK(1) == 1 )THEN
    Call Broadcast_Array(SAL_Global,      master_id)
    Call Broadcast_Array(SAL1_Global,     master_id)
  ENDIF

  IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 .AND. ISCOCHK(2) == 1 )THEN
    Call Broadcast_Array(TEM_Global,      master_id)
    Call Broadcast_Array(TEM1_Global,     master_id)
    Call Broadcast_Array(TEMB_Global,     master_id)
    IF( ISICE > 0 )THEN
      Call Broadcast_Array(ICETHICK_Global, master_id)
      Call Broadcast_Array(ICETEMP_Global,  master_id)
    ENDIF
  ENDIF

  IF( ISTRAN(3) > 0 .AND. ISCI(3) > 0 .AND. ISCOCHK(3) == 1 )THEN
    Call Broadcast_Array(DYE_Global,      master_id)
    Call Broadcast_Array(DYE1_Global,     master_id)
  ENDIF

  IF( ISTRAN(4) > 0 .AND. ISCI(4) > 0 .AND. ISCOCHK(4) == 1 )THEN
    ! PLACE HOLDER
  ENDIF

  IF(  ISTRAN(5) > 0 .AND. ISCI(5) > 0 .AND. ISCOCHK(5) == 1 )THEN
    Call Broadcast_Array(TOXB_Global,     master_id)
    Call Broadcast_Array(TOXB1_Global,    master_id)
    Call Broadcast_Array(TOX_Global,      master_id)
    Call Broadcast_Array(TOX1_Global,     master_id)
  ENDIF

  IDFLAG = 0
  IF( ISTRAN(6) > 0 .AND. ISCI(6) > 0 .AND. ISCOCHK(6) == 1 )THEN
    Call Broadcast_Array(SEDB_Global,     master_id)
    Call Broadcast_Array(SEDB1_Global,    master_id)
    Call Broadcast_Array(SED_Global,      master_id)
    Call Broadcast_Array(SED1_Global,     master_id)
    IDFLAG = 1
  ENDIF

  IF( ISTRAN(7) > 0 .AND. ISCI(7) > 0 .AND. ISCOCHK(7) == 1 )THEN
    Call Broadcast_Array(SNDB_Global,     master_id)
    Call Broadcast_Array(SNDB1_Global,    master_id)
    Call Broadcast_Array(SND_Global,      master_id)
    Call Broadcast_Array(SND1_Global,     master_id)
    IDFLAG = 1
  ENDIF

  IF( IDFLAG == 1 )THEN
    Call Broadcast_Array(HBED_Global,     master_id)
    Call Broadcast_Array(HBED1_Global,    master_id)
    Call Broadcast_Array(VDRBED_Global,   master_id)
    Call Broadcast_Array(VDRBED1_Global,  master_id)
  ENDIF

  ! *** Populate local arrays
  DO LG=2,LA_Global
    LL = Map2Local(LG).LL
    IF( LL > 0 )THEN
      BELV(LL)   = BELV_Global(LG)
      HP(LL)     = HP_Global(LG)
      H1P(LL)    = H1P_Global(LG)
      HWQ(LL)    = HWQ_Global(LG)
      H2WQ(LL)   = H2WQ_Global(LG)
      UHDYE(LL)  = UHDYE_Global(LG)
      UHDY1E(LL) = UHDY1E_Global(LG)
      VHDXE(LL)  = VHDXE_Global(LG)
      VHDX1E(LL) = VHDX1E_Global(LG)
      IF( Restart_In_Ver > 1000 )THEN
        SUB(LL)  = SUB_Global(LG)
        SVB(LL)  = SVB_Global(LG)
        ISCDRY(LL) = ISCDRY_Global(LG)
        NATDRY(LL) = NATDRY_Global(LG)
        IF( IDRY_Global(LG) == 0 )THEN
          LMASKDRY(LL) = .TRUE.
        ELSE
          LMASKDRY(LL) = .FALSE.
        ENDIF
      ENDIF
      DO K = 1,KC
        U(LL,K)  = U_Global(LG,K)
        V(LL,K)  = V_Global(LG,K)
        U1(LL,K) = U1_Global(LG,K)
        V1(LL,K) = V1_Global(LG,K)
      ENDDO
      DO K = 0,KC
        QQ(LL,K)   = MAX(QQ_Global(LG,K),QQMIN)
        QQ1(LL,K)  = MAX(QQ1_Global(LG,K),QQMIN)
        QQL(LL,K)  = MAX(QQL_Global(LG,K),QQLMIN)
        QQL1(LL,K) = MAX(QQL1_Global(LG,K),QQLMIN)
        DML(LL,K)  = MAX(DML_Global(LG,K),DMLMIN)
      ENDDO

      IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 .AND. ISCOCHK(2) == 1 .AND. ISICE > 2 )THEN
        ICETHICK(LL) = ICETHICK_Global(LG)
        ICETEMP(LL)  = ICETEMP_Global(LG)
        IF( ICETHICK(LL) > 0.0 ) ICECELL(LL) = .TRUE.
      ENDIF

      IF( ISTRAN(1) > 0 .AND. ISCI(1) > 0 .AND. ISCOCHK(1) == 1 )THEN
        SAL(LL,1:KC)  = SAL_Global(LG,1:KC)
        SAL1(LL,1:KC) = SAL1_Global(LG,1:KC)
      ENDIF

      IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 .AND. ISCOCHK(2) == 1 )THEN
        TEM(LL,1:KC)  = TEM_Global(LG,1:KC)
        TEM1(LL,1:KC) = TEM1_Global(LG,1:KC)
        TEMB(LL)      = TEMB_Global(LG)

        ! *** ENSURE VALID VALUES
        !TEM(LL,1:KC)  = MAX(MIN(TEM(LL,1:KC),70.),-5.)                  ! *** Allow negative temperatures for ice formation process
        !TEM1(LL,1:KC) = MAX(MIN(TEM1(LL,1:KC),70.),-5.)
        !TEMB(LL)      = MAX(MIN(TEMB(LL),50.),0.)
      ENDIF

      IF( ISTRAN(3) > 0 .AND. ISCI(3) > 0 .AND. ISCOCHK(3) == 1 )THEN
        DO MD=1,NDYE
          DYE(LL,1:KC,MD)   = DYE_Global(LG,1:KC,MD)
          DYE1(LL,1:KC,MD)  = DYE1_Global(LG,1:KC,MD)
        ENDDO
      ENDIF

      IF( ISTRAN(4) > 0 .AND. ISCI(4) > 0 .AND. ISCOCHK(4) == 1 )THEN
        ! *** delme - not MPI
      ENDIF

      IF(  ISTRAN(5) > 0 .AND. ISCI(5) > 0 .AND. ISCOCHK(5) == 1 )THEN
        DO NT=1,NTXM
          TOX(LL,1:KC,NT)   = TOX_Global(LG,1:KC,NT)
          TOX1(LL,1:KC,NT)  = TOX1_Global(LG,1:KC,NT)
          TOXB(LL,1:KB,NT)  = TOXB_Global(LG,1:KB,NT)
          TOXB1(LL,1:KB,NT) = TOXB1_Global(LG,1:KB,NT)
        ENDDO
      ENDIF

      IDFLAG = 0
      IF( ISTRAN(6) > 0 .AND. ISCI(6) > 0 .AND. ISCOCHK(6) == 1 )THEN
        DO NS=1,NSCM
          SED(LL,1:KC,NS)   = SED_Global(LG,1:KC,NS)
          SED1(LL,1:KC,NS)  = SED1_Global(LG,1:KC,NS)
          IF( NSEDFLUME == 0 .OR. (NSEDFLUME == 3 .AND. IHTSTRT == 0) )THEN
            SEDB(LL,1:KB,NS)  = SEDB_Global(LG,1:KB,NS)
            SEDB1(LL,1:KB,NS) = SEDB1_Global(LG,1:KB,NS)
          ENDIF
        ENDDO
        IF( NSCM2 > NSCM )THEN
          DO NS=NSCM+1,NSCM2
            SED(LL,1:KC,NS)   = SED_Global(LG,1:KC,NS)
            SED1(LL,1:KC,NS)  = SED1_Global(LG,1:KC,NS)
          ENDDO
        ENDIF
        IDFLAG = 1
      ENDIF

      IF( ISTRAN(7) > 0 .AND. ISCI(7) > 0 .AND. ISCOCHK(7) == 1 )THEN
        DO NS=1,NSNM
          SND(LL,1:KC,NS)   = SND_Global(LG,1:KC,NS)
          SND1(LL,1:KC,NS)  = SND1_Global(LG,1:KC,NS)
          IF( NSEDFLUME == 0 .OR. (NSEDFLUME == 3 .AND. IHTSTRT == 0) )THEN
            SNDB(LL,1:KB,NS)  = SNDB_Global(LG,1:KB,NS)
            SNDB1(LL,1:KB,NS) = SNDB1_Global(LG,1:KB,NS)
          ENDIF
        ENDDO
        IDFLAG = 1
      ENDIF

      IF( NSEDFLUME == 0 .OR. (NSEDFLUME == 3 .AND. IHTSTRT == 0) )THEN
        IF( IDFLAG == 1 )THEN
          HBED(LL,1:KB)    = HBED_Global(LG,1:KB)
          HBED1(LL,1:KB)   = HBED1_Global(LG,1:KB)
          VDRBED(LL,1:KB)  = VDRBED_Global(LG,1:KB)
          VDRBED1(LL,1:KB) = VDRBED1_Global(LG,1:KB)
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  ! *****************************************************************************************************
  ! *** BOUNDARY CONDITIONS
  if( process_id == master_id )then
    MD = 0
    NT = 0
    NS = 0
    NX = 0
    DO I = 1,7
      IF( ISCOCHK(I) == 1 )THEN
100     IF( I == 3 )THEN
          MD = MD + 1
          M = MSVDYE(MD)
        ELSEIF( I == 4 )THEN
          M = 3 + NDYE
        ELSEIF( I == 5 )THEN
          NT = NT + 1
          M = MSVTOX(NT)
        ELSEIF( I == 6 )THEN
          NS = NS + 1
          M = MSVSED(NS)
        ELSEIF( I == 7 )THEN
          NX = NX + 1
          M = MSVSND(NX)
        ELSE
          M = I
        ENDIF

        IF( ISTRAN(I) > 0 .AND. ISCI(I) > 0 )THEN
          DO LL=1,NPBS_GL
            READ(UINP,*,ERR=1023) (NLOS_Global(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1024) (CLOS_Global(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NPBW_GL
            READ(UINP,*,ERR=1025) (NLOW_Global(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1026) (CLOW_Global(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NPBE_GL
            READ(UINP,*,ERR=1027) (NLOE_Global(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1028) (CLOE_Global(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NPBN_GL
            READ(UINP,*,ERR=1029) (NLON_Global(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1030) (CLON_Global(LL,K,M),K=1,KC)
          ENDDO
        ELSE
          DO LL=1,NPBS_GL
            READ(UINP,*,ERR=1023) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1024) (TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NPBW_GL
            READ(UINP,*,ERR=1025) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1026) (TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NPBE_GL
            READ(UINP,*,ERR=1027) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1028) (TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NPBN_GL
            READ(UINP,*,ERR=1029) (TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1030) (TMPVAL,K=1,KC)
          ENDDO
        ENDIF

        IF( I == 3 .AND. MD < NDYE ) GOTO 100      ! *** LOOP BACK AND READ NEXT DYE CLASS
        IF( I == 5 .AND. NT < NTXM ) GOTO 100      ! *** LOOP BACK AND READ NEXT TOX CLASS
        IF( I == 6 .AND. NS < NSCM ) GOTO 100      ! *** LOOP BACK AND READ NEXT SED CLASS
        IF( I == 7 .AND. NX < NSNM ) GOTO 100      ! *** LOOP BACK AND READ NEXT SND CLASS
      ENDIF
    ENDDO
  endif

  ! *****************************************************************************************************
  ! *** Communcate boundary condition arrays

  ! *** Build local arrays
  IF( NPBS_GL > 0 )THEN
    Call Broadcast_Array(NLOS_Global, master_id)
    Call Broadcast_Array(CLOS_Global, master_id)
    NLOS = 0
    CLOS = 0.

    LL = 0
    DO I = 1,NPBS_GL
      III = IG2IL(IPBS_GL(I))                    ! Get local I value
      JJJ = JG2JL(JPBS_GL(I))                    ! Get local J value
      IF(  III > 0 .AND. III < IC )THEN
        IF(  JJJ > 0 .AND. JJJ < JC )THEN
          LL = LL + 1
          NLOS(LL,1:KC,1:NSTVM) = NLOS_Global(I,1:KC,1:NSTVM)
          CLOS(LL,1:KC,1:NSTVM) = CLOS_Global(I,1:KC,1:NSTVM)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  IF( NPBE_GL > 0 )THEN
    Call Broadcast_Array(NLOE_Global, master_id)
    Call Broadcast_Array(CLOE_Global, master_id)
    NLOE = 0
    CLOE = 0.

    LL = 0
    DO I = 1,NPBE_GL
      III = IG2IL(IPBE_GL(I))                    ! Get local I value
      JJJ = JG2JL(JPBE_GL(I))                    ! Get local J value
      IF(  III > 0 .AND. III < IC )THEN
        IF(  JJJ > 0 .AND. JJJ < JC )THEN
          LL = LL + 1
          NLOE(LL,1:KC,1:NSTVM) = NLOE_Global(I,1:KC,1:NSTVM)
          CLOE(LL,1:KC,1:NSTVM) = CLOE_Global(I,1:KC,1:NSTVM)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  IF( NPBW_GL > 0 )THEN
    Call Broadcast_Array(NLOW_Global, master_id)
    Call Broadcast_Array(CLOW_Global, master_id)
    NLOW = 0
    CLOW = 0.

    LL = 0
    DO I = 1,NPBW_GL
      III = IG2IL(IPBW_GL(I))                    ! Get local I value
      JJJ = JG2JL(JPBW_GL(I))                    ! Get local J value
      IF(  III > 0 .AND. III < IC )THEN
        IF(  JJJ > 0 .AND. JJJ < JC )THEN
          LL = LL + 1
          NLOW(LL,1:KC,1:NSTVM) = NLOW_Global(I,1:KC,1:NSTVM)
          CLOW(LL,1:KC,1:NSTVM) = CLOW_Global(I,1:KC,1:NSTVM)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  IF( NPBN_GL > 0 )THEN
    Call Broadcast_Array(NLON_Global, master_id)
    Call Broadcast_Array(CLON_Global, master_id)
    NLON = 0
    CLON = 0.

    LL = 0
    DO I = 1,NPBN_GL
      III = IG2IL(IPBN_GL(I))                    ! Get local I value
      JJJ = JG2JL(JPBN_GL(I))                    ! Get local J value
      IF(  III > 0 .AND. III < IC )THEN
        IF(  JJJ > 0 .AND. JJJ < JC )THEN
          LL = LL + 1
          NLON(LL,1:KC,1:NSTVM) = NLON_Global(I,1:KC,1:NSTVM)
          CLON(LL,1:KC,1:NSTVM) = CLON_Global(I,1:KC,1:NSTVM)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  ! *** BOUNDARY FLOWS
  if( process_id == master_id )then
    DO L=2,LA_Global
      READ(UINP,*,ERR=1031) QSUME_Global(L), (QSUM_Global(L,K),K=1,KC)
    ENDDO
  endif
  Call Broadcast_Array(QSUME_Global, master_id)
  Call Broadcast_Array(QSUM_Global,  master_id)

  ! *** Populate local arrays
  DO LG=2,LA_Global
    LL = Map2Local(LG).LL
    IF( LL > 0 )THEN
      QSUME(LL)     = QSUME_Global(LG)
      QSUM(LL,1:KC) = QSUM_Global(LG,1:KC)
    ENDIF
  ENDDO

  IF( MDCHH >= 1 )THEN   ! *** delme - Not supported by MPI
    if( num_Processors > 1 ) Call STOPP('MPI does not support Channel Modifiers')
    DO NMD=1,MDCHH
      READ(UINP,*,ERR=1032) ITMP1,JTMP1,ITMP2,JTMP2,ITMP3,JTMP3,QCHANU(NMD),QCHANV(NMD)
    ENDDO
  ELSE
    DO NMD=1,MDCHH
      QCHANU(NMD)=0.
      QCHANV(NMD)=0.
      QCHANUN(NMD)=0.
      QCHANVN(NMD)=0.
    ENDDO
  ENDIF

  IF( ISGWIE >= 1 )THEN
    if( num_Processors > 1 ) Call STOPP('MPI does not support ISGWIE > 0 option')
    DO L=2,LA_Global
      READ(UINP,*,ERR=1033) AGWELV(L), AGWELV1(L)
    ENDDO
  ENDIF

  ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE WITH A LOWER CHORD ACTIVATION TOGGLE
  DO NCTL=1,NQCTL
    IF( HYD_STR(NCTL).NQCTYP == 3 .OR. HYD_STR(NCTL).NQCTYP == 4 )THEN
      if( num_Processors > 1 ) Call STOPP('MPI does not support NQCTYP = 3 or 4 option')
      ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE (4) DEPENDANT FLOWS
      READ(1,*,ERR=1034) LL,NLOWCHORD(NCTL)
      READ(1,*,ERR=1034) SAVESUB(1,NCTL),SAVESVB(1,NCTL),SAVESUB(2,NCTL),SAVESVB(2,NCTL)
      READ(1,*,ERR=1034) LOWCHORDU(NCTL),LOWCHORDV(NCTL)

      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)

      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      IF( ID == 0 .AND. JD == 0 )THEN
        ! *** INVALID SPECIFICATION
        CYCLE
      ENDIF
      LD=LIJ(ID,JD)

      IF( NLOWCHORD(NCTL) < 1 )THEN
        IF( ID > IU )THEN
          IF( SAVESUB(2,NCTL) > 0.5 )THEN
            SUB(LD)   = 1.0
            SUBO(LD)  = 1.0
          ENDIF
        ELSE
          IF( SAVESUB(1,NCTL) > 0.5 )THEN
            SUB(LU)   = 1.0
            SUBO(LU)  = 1.0
          ENDIF
        ENDIF
        IF( JD > JU )THEN
          IF( SAVESVB(2,NCTL) > 0.5 )THEN
            SVB(LD)   = 1.0
            SVBO(LD)  = 1.0
          ENDIF
        ELSE
          IF( SAVESVB(1,NCTL) > 0.5 )THEN
            SVB(LU)=1.0
            SVBO(LU)=1.0
          ENDIF
        ENDIF
      ELSE
        ! *** SET THE CELL FACE FLAGS - BLOCK FULL HYDRODYNAMICS
        IF( ID > IU )THEN
          SUB(LD)=0.0
          SUBO(LD)=0.0
          UHDYE(LD)=0.0
        ELSE
          SUB(LU)=0.0
          SUBO(LU)=0.0
          UHDYE(LU)=0.0
        ENDIF
        IF( JD > JU )THEN
          SVB(LD)=0.0
          SVBO(LD)=0.0
          VHDXE(LD)=0.0
        ELSE
          SVB(LU)=0.0
          SVBO(LU)=0.0
          VHDXE(LU)=0.0
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  ! *** CLOSE MAIN RESTART.OUT FILE
  if( process_id == master_id ) CLOSE(UINP)

  ! *****************************************************************************************************
  ! *** Compute local intermediate arrays

6666 FORMAT(3I10,F12.6)
6667 FORMAT(7I5,2X,E12.4,2X,E12.4)
  DO K=1,KC
    SAL(1,K)=0.
    SAL1(1,K)=0.
    TEM(1,K)=0.
    TEM1(1,K)=0.
    DYE(1,K,1:NDYE)=0.
    DYE1(1,K,1:NDYE)=0.
    SFL(1,K)=0.
    SFL2(1,K)=0.
    IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
      SED(1,K,1:NSED)=0.
      SED1(1,K,1:NSED)=0.
      SND(1,K,1:NSND)=0.
      SND1(1,K,1:NSND)=0.
      IF( ISTRAN(5) > 0 )THEN
        TOX(1,K,1:NTOX)=0.
        TOX1(1,K,1:NTOX)=0.
      ENDIF
    ENDIF

    VHDX(1,K)=0.
    UHDY(1,K)=0.
    VHDXF(1,K)=0.
    UHDYF(1,K)=0.
    VHDX1(1,K)=0.
    UHDY1(1,K)=0.
    VHDXF1(1,K)=0.
    UHDYF1(1,K)=0.

    SAL(LC,K)=0.
    SAL1(LC,K)=0.
    TEM(LC,K)=0.
    TEM1(LC,K)=0.
    DYE(LC,K,1:NDYE)=0.
    DYE1(LC,K,1:NDYE)=0.
    SFL(LC,K)=0.
    SFL2(LC,K)=0.
    IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
      SED(LC,K,1:NSED)=0.
      SED1(LC,K,1:NSED)=0.
      SND(LC,K,1:NSND)=0.
      SND1(LC,K,1:NSND)=0.
      IF( ISTRAN(5) > 0 )THEN
        TOX(LC,K,1:NTOX)=0.
        TOX1(LC,K,1:NTOX)=0.
      ENDIF
    ENDIF

    VHDX(LC,K)=0.
    UHDY(LC,K)=0.
    VHDXF(LC,K)=0.
    UHDYF(LC,K)=0.
    VHDX1(LC,K)=0.
    UHDY1(LC,K)=0.
    VHDXF1(LC,K)=0.
    UHDYF1(LC,K)=0.
  ENDDO

  DO L=2,LA
    UHDYE(L)  = SUB(L)*UHDYE(L)
    UHDY1E(L) = SUB(L)*UHDY1E(L)
    VHDXE(L)  = SVB(L)*VHDXE(L)
    VHDX1E(L) = SVB(L)*VHDX1E(L)
  ENDDO
  DO K=1,KC
    DO L=2,LA
      U(L,K)  = SUB(L)*U(L,K)
      U1(L,K) = SUB(L)*U1(L,K)
      V(L,K)  = SVB(L)*V(L,K)
      V1(L,K) = SVB(L)*V1(L,K)
    ENDDO
  ENDDO
  DO L=2,LA
    LS=LSC(L)
    H1U(L) = 0.5*(DXP(L)*DYP(L)*H1P(L)+DXP(LWC(L))*DYP(LWC(L))*H1P(LWC(L)))/(DXU(L)*DYU(L))
    H1V(L) = 0.5*(DXP(L)*DYP(L)*H1P(L)+DXP(LSC(L))*DYP(LSC(L))*H1P(LSC(L)))/(DXV(L)*DYV(L))
    P1(L)  = G*(H1P(L)+BELV(L))
    HU(L)  = 0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LWC(L))*DYP(LWC(L))*HP(LWC(L)))/(DXU(L)*DYU(L))
    HV(L)  = 0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LSC(L))*DYP(LSC(L))*HP(LSC(L)))/(DXV(L)*DYV(L))
    P(L)   = G*(HP(L)+BELV(L))
    HPI(L) = 1./HP(L)
    HUI(L) = 1./HU(L)
    HVI(L) = 1./HV(L)
    H1UI(L) = 1./H1U(L)
    H1VI(L) = 1./H1V(L)
  ENDDO
  H1U(1) = H1U(2)
  H1V(1) = H1V(2)
  P1(1) = P1(2)
  HU(1) = HU(2)
  HV(1) = HV(2)
  P(1) = P(2)
  HPI(1) = 1./HP(2)
  HUI(1) = 1./HU(2)
  HVI(1) = 1./HV(2)
  H1UI(1) = 1./H1U(2)
  H1VI(1) = 1./H1V(2)
  H1U(LC) = H1U(LA)
  H1V(LC) = H1V(LA)
  P1(LC) = P1(LA)
  HU(LC) = HU(LA)
  HV(LC) = HV(LA)
  P(LC) = P(LA)
  HPI(LC) = 1./HP(LA)
  HUI(LC) = 1./HU(LA)
  HVI(LC) = 1./HV(LA)
  H1UI(LC) = 1./H1U(LA)
  H1VI(LC) = 1./H1V(LA)

  DO K=1,KC
    DO L=2,LA
      UHDYF1(L,K) = DYU(L)*H1U(L)*U1(L,K)
      VHDXF1(L,K) = DXV(L)*H1V(L)*V1(L,K)
      UHDY1(L,K)  = UHDYF1(L,K)*DZC(L,K)
      VHDX1(L,K)  = VHDXF1(L,K)*DZC(L,K)

      UHDYF(L,K)  = DYU(L)*HU(L)*U(L,K)
      VHDXF(L,K)  = DXV(L)*HV(L)*V(L,K)
      UHDY(L,K)   = UHDYF(L,K)*DZC(L,K)
      VHDX(L,K)   = VHDXF(L,K)*DZC(L,K)

      SAL(L,K)  = MAX(SAL(L,K),0.)
      SAL1(L,K) = MAX(SAL1(L,K),0.)

      TEM(L,K)  = MAX(TEM(L,K),0.)
      TEM1(L,K) = MAX(TEM1(L,K),0.)

      DO MD = 1,NDYE
        DYE(L,K,MD)  = MAX(DYE(L,K,MD),0.)
        DYE1(L,K,MD) = MAX(DYE1(L,K,MD),0.)
      ENDDO

      DO NT = 1,NTOX
        TOX(L,K,NT)  = MAX(TOX(L,K,NT),0.)
        TOX1(L,K,NT) = MAX(TOX1(L,K,NT),0.)
      ENDDO

      DO NS = 1,NSED
        SED(L,K,NS)  = MAX(SED(L,K,NS),0.)
        SED1(L,K,NS) = MAX(SED1(L,K,NS),0.)
      ENDDO

      DO NX = 1,NSND
        SND(L,K,NX)  = MAX(SND(L,K,NX),0.)
        SND1(L,K,NX) = MAX(SND1(L,K,NX),0.)
      ENDDO

    ENDDO
  ENDDO

  ! *** CORRECT FOR CHANGED BOTTOM ELEV
  Call Broadcast_Scalar(ISBELVC,  master_id)
  IF( ISRESTIOPT == 0 .AND. ISBELVC == 1 )THEN
    DO L=2,LA
      UHE(L)=0.
      VHE(L)=0.
    ENDDO
    DO K=1,KC
      DO L=2,LA
        UHE(L)=UHE(L)+UHDYF1(L,K)
        VHE(L)=VHE(L)+VHDXF1(L,K)
      ENDDO
    ENDDO
    DO L=2,LA
      IF( UHE(L) /= 0. )THEN
        TMPVAL=UHDY1E(L)/UHE(L)
        DO K=1,KC
          U1(L,K)=TMPVAL*U1(L,K)
        ENDDO
      ENDIF
      IF( VHE(L) /= 0. )THEN
        TMPVAL=VHDX1E(L)/VHE(L)
        DO K=1,KC
          V1(L,K)=TMPVAL*V1(L,K)
        ENDDO
      ENDIF
    ENDDO
    DO L=2,LA
      UHE(L)=0.
      VHE(L)=0.
    ENDDO
    DO K=1,KC
      DO L=2,LA
        UHE(L)=UHE(L)+UHDYF(L,K)
        VHE(L)=VHE(L)+VHDXF(L,K)
      ENDDO
    ENDDO

    DO L=2,LA
      IF( UHE(L) /= 0. )THEN
        TMPVAL=UHDYE(L)/UHE(L)
        DO K=1,KC
          U(L,K)=TMPVAL*U(L,K)
        ENDDO
      ENDIF
      IF( VHE(L) /= 0. )THEN
        TMPVAL=VHDXE(L)/VHE(L)
        DO K=1,KC
          V(L,K)=TMPVAL*V(L,K)
        ENDDO
      ENDIF
    ENDDO

    DO K=1,KC
      DO L=2,LA
        UHDYF1(L,K) = DYU(L)*H1U(L)*U1(L,K)
        VHDXF1(L,K) = DXV(L)*H1V(L)*V1(L,K)
        UHDY1(L,K)  = UHDYF1(L,K)*DZC(L,K)
        VHDX1(L,K)  = VHDXF1(L,K)*DZC(L,K)

        UHDYF(L,K) = DYU(L)*HU(L)*U(L,K)
        VHDXF(L,K) = DXV(L)*HV(L)*V(L,K)
        UHDY(L,K)  = UHDYF(L,K)*DZC(L,K)
        VHDX(L,K)  = VHDXF(L,K)*DZC(L,K)
      ENDDO
    ENDDO
  ENDIF
  N=0

  IF( ISDRY == 0 )THEN
    CALL CALTSXY
    CALL CALQVS (2)
  ENDIF

  DO K=1,KS
    DO L=2,LA
      IF( LKSZ(L,K) )CYCLE
      RDZC = DZC(L,K)
      W(L,K)  = SWB(L)*( W(L,K-1)  - RDZC*(UHDYF(LEC(L),K) - UHDYF(L,K) - UHDYE(LEC(L)) + UHDYE(L) &
        + VHDXF(LNC(L),K) - VHDXF(L,K) - VHDXE(LNC(L)) + VHDXE(L))*DXYIP(L)) &
        + SWB(L)*( QSUM(L,K) - RDZC*QSUME(L) )*DXYIP(L)

      W1(L,K) = SWB(L)*( W1(L,K-1) - RDZC*(UHDYF1(LEC(L),K) - UHDYF1(L,K) - UHDY1E(LEC(L)) + UHDY1E(L) &
        + VHDXF1(LNC(L),K) - VHDXF1(L,K) - VHDX1E(LNC(L)) + VHDX1E(L))*DXYIP(L)) &
        + SWB(L)*( QSUM(L,K) - RDZC*QSUME(L) )*DXYIP(L)
    ENDDO
  ENDDO

  ! **  SET DRYING AND WETTING FLAGS
  IF( ISDRY > 0 .AND. ISDRY < 99 )THEN
    DO L=2,LA
      ISCDRY(L)=0
      LS=LSC(L)
      LN=LNC(L)
      IF( HP(L) <= HDRY )THEN
        ISCDRY(L)=1
        SUB(L)=0.
        SVB(L)=0.
        SUB(LEC(L))=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(LEC(L))=0.
        SBY(LN)=0.
      ENDIF
    ENDDO
  ENDIF

  ! *****************************************************************************************************
  IF( ISDRY == 99 .AND. Restart_In_Ver < 1010 )THEN
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING RESTART FILE: RSTWD.INP'
      OPEN(1,FILE='rstwd.inp',STATUS='UNKNOWN')

      DO LG=2,LA_Global
        READ(1,*) LDUM, IDUM, JDUM, ISCDRY_Global(LG), NATDRY_Global(LG), IDRY_Global(LG), SUB_Global(LG), SVB_Global(LG)
      ENDDO

      CLOSE(1)
    endif

    ! *** Populate local arrays
    DO LG=2,LA_Global
      LL = Map2Local(LG).LL
      IF( LL > 0 )THEN
        SUB(LL)  = SUB_Global(LG)
        SVB(LL)  = SVB_Global(LG)
        ISCDRY(LL) = ISCDRY_Global(LG)
        NATDRY(LL) = NATDRY_Global(LG)
        IF( IDRY_Global(LG) == 0 )THEN
          LMASKDRY(LL) = .TRUE.
        ELSE
          LMASKDRY(LL) = .FALSE.
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  RETURN

101 FORMAT(I5)
102 FORMAT(3I5,12F8.2)

  ! **  WRITE READ ERRORS ON RESTART
1000 WRITE(6,2000)
  CALL STOPP('.')

1001 WRITE(6,2001)L
  CALL STOPP('.')

1002 WRITE(6,2002)L
  CALL STOPP('.')

1003 WRITE(6,2003)L
  CALL STOPP('.')

1004 WRITE(6,2004)L
  CALL STOPP('.')

1005 WRITE(6,2005)L
  CALL STOPP('.')

1006 WRITE(6,2006)L
  CALL STOPP('.')

1007 WRITE(6,2007)L
  CALL STOPP('.')

1008 WRITE(6,2008)L
  CALL STOPP('.')

1009 WRITE(6,2009)L
  CALL STOPP('.')

1010 WRITE(6,2010)L
  CALL STOPP('.')

1011 WRITE(6,2011)L
  CALL STOPP('.')

1012 WRITE(6,2012)L
  CALL STOPP('.')

1013 WRITE(6,2013)L
  CALL STOPP('.')

1014 WRITE(6,2014)L
  CALL STOPP('.')

1015 WRITE(6,2015)L
  CALL STOPP('.')

1016 WRITE(6,2016)L
  CALL STOPP('.')

1017 WRITE(6,2017)L
  CALL STOPP('.')

1018 WRITE(6,2018)L
  CALL STOPP('.')

1019 WRITE(6,2019)L
  CALL STOPP('.')

1020 WRITE(6,2020)L
  CALL STOPP('.')

1021 WRITE(6,2021)L
  CALL STOPP('.')

1022 WRITE(6,2022)L
  CALL STOPP('.')

1023 WRITE(6,2023)L
  CALL STOPP('.')

1024 WRITE(6,2024)L
  CALL STOPP('.')

1025 WRITE(6,2025)L
  CALL STOPP('.')

1026 WRITE(6,2026)L
  CALL STOPP('.')

1027 WRITE(6,2027)L
  CALL STOPP('.')

1028 WRITE(6,2028)L
  CALL STOPP('.')

1029 WRITE(6,2029)L
  CALL STOPP('.')

1030 WRITE(6,2030)L
  CALL STOPP('.')

1031 WRITE(6,2031)L
  CALL STOPP('.')

1032 WRITE(6,2032)NMD
  CALL STOPP('.')

1033 WRITE(6,2033)L
  CALL STOPP('.')

1034 WRITE(6,2034)NCTL
  CALL STOPP('.')

1035 WRITE(6,2035)L
  CALL STOPP('.')

600 FORMAT(2X,'I,J,BELVOLD,BELVNEW',2I5,2F12.2)
2000 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1000')
2001 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1001 L =',I6)
2002 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1002 L =',I6)
2003 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1003 L =',I6)
2004 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1004 L =',I6)
2005 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1005 L =',I6)
2006 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1006 L =',I6)
2007 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1007 L =',I6)
2008 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1008 L =',I6)
2009 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1009 L =',I6)
2010 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1010 L =',I6)
2011 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1011 L =',I6)
2012 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1012 L =',I6)
2013 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1013 L =',I6)
2014 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1014 L =',I6)
2015 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1015 L =',I6)
2016 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1016 L =',I6)
2017 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1017 L =',I6)
2018 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1018 L =',I6)
2019 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1019 L =',I6)
2020 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1020 L =',I6)
2021 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1021 L =',I6)
2022 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1022 L =',I6)
2023 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1023 L =',I6)
2024 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1024 L =',I6)
2025 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1025 L =',I6)
2026 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1026 L =',I6)
2027 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1027 L =',I6)
2028 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1028 L =',I6)
2029 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1029 L =',I6)
2030 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1030 L =',I6)
2031 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1031 L =',I6)
2032 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1032 NMD =',I6)
2033 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1033 L =',I6)
2034 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1034 NCTL =',I6)
2035 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1035 L =',I6)
9696 FORMAT('  NEGATIVE DEPTH RESTART, L,I,J,HP,H1P = ',3I7,2F10.4)
9698 FORMAT('  ZERO DEPTH RESTART (WARN), L,I,J,HP,H1P = ',3I7,2F10.4)

  RETURN
  
END SUBROUTINE Restart_In

SUBROUTINE WQSDRST_IN
  ! ** READ ICS FROM RESTART FILE FROM INSMRST.
  ! ** THIS FILE IS ONLY USED WHEN iICI=2 (WQ_3DSD.INP)
  ! ** CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE
  ! ** THE ASCII FILE INSTEAD.
  !
  LOGICAL FEXIST
  INTEGER :: NREST, M, L, LG, NW
  REAL    :: RSTTIME

  INQUIRE(FILE='wqsdrst.bin', EXIST=FEXIST)
    
  IF( .NOT. FEXIST )THEN
    ! *** FORMATTED
    WRITE(*,'(A)')' WQ: READING WQSDRST.INP'
    OPEN(1,FILE='wqsdrst.inp',STATUS='UNKNOWN')

    CALL SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES
    DO M=2,LA_Global
      READ(1,*) L,(SMPON_Global(L,NW),NW=1,NSMG), (SMPOP_Global(L,NW),NW=1,NSMG), (SMPOC_Global(L,NW),NW=1,NSMG),  &
                   SM1NH4_Global(L), SM2NH4_Global(L), SM2NO3_Global(L), SM2PO4_Global(L), SM2H2S_Global(L),       &
                   SMPSI_Global(L),  SM2SI_Global(L),  SMBST_Global(L),  SMT_Global(L)
    ENDDO
  ELSE
    ! *** BINARY
    WRITE(*,'(A)')' WQ: READING WQSDRST.BIN'
    OPEN(1,FILE='wqsdrst.bin',STATUS='UNKNOWN',FORM='BINARY')
    
    READ(1) NREST, RSTTIME
    WRITE(0,911) NREST, RSTTIME
911 FORMAT(' READING BINARY WQSDRST.BIN FILE ...    NN, TIME = ', I10, F11.5)
    DO M=2,LA_Global
      READ(1) L
      READ(1) (SMPON_Global(L,NW),NW=1,NSMG), (SMPOP_Global(L,NW),NW=1,NSMG), (SMPOC_Global(L,NW),NW=1,NSMG),  &
               SM1NH4_Global(L), SM2NH4_Global(L), SM2NO3_Global(L), SM2PO4_Global(L), SM2H2S_Global(L),       &
               SMPSI_Global(L),  SM2SI_Global(L),  SMBST_Global(L),  SMT_Global(L)
    ENDDO
  ENDIF

  CLOSE(1)

 90 FORMAT(I5, 18E12.4)
999 FORMAT(1X)

END SUBROUTINE WQSDRST_IN

SUBROUTINE WQSDRST_OUT
  ! ** WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT ISMORST.
  INTEGER      :: L, NW, NHR
  CHARACTER(2) :: SNHR

  if( process_id == master_id )then
    IF( ISGREGOR ==  1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = OUTDIR//'WQSDRST_'//TRIM(SDATE)//'_'//SNHR//'.OUT'
    ELSE
      RESFILE = OUTDIR//'WQSDRST.OUT'
    ENDIF
    OPEN(1,FILE=RESFILE,STATUS='REPLACE')

    ! *** BEGIN WQSDRST
    WRITE(1,101) N, TIMEDAY
    WRITE(1,888)
    DO L=2,LA_Global
      WRITE(1,90) L,(SMPON_Global(L,NW),NW=1,NSMG), (SMPOP_Global(L,NW),NW=1,NSMG), (SMPOC_Global(L,NW),NW=1,NSMG),  &
                     SM1NH4_Global(L), SM2NH4_Global(L), SM2NO3_Global(L), SM2PO4_Global(L), SM2H2S_Global(L),       &
                     SMPSI_Global(L),  SM2SI_Global(L),  SMBST_Global(L),  SMT_Global(L)
    ENDDO
    CLOSE(1)
  endif
  
 90 FORMAT(I5, 1P, 18E12.4)
101 FORMAT('C  SM RESTART FILE TIME STEP, TIME = ',I10,F13.5)
888 FORMAT('C    L',                                                      &
           '       GPON1       GPON2       GPON3       GPOP1    GPOP2',  &
           '       GPOP3       GPOC1       GPOC2       GPOC3    G1NH4',  &
           '       G2NH4       G2NO3       G2PO4       G2H2S    GPSI',   &
           '        G2SI        GBST         GT')

END SUBROUTINE WQSDRST_OUT

SUBROUTINE WQ_WCRST_IN
  ! ** INITIAL CONDITION FILE
  ! ** READ ICS FROM RESTART FILE FROM INWQRST(WQ_3DWC.INP: IWQICI == 2 ).

  INTEGER        :: L, LG, K, NREST, NWQV0, M, NW, LL, KK, NCOL
  INTEGER(IK8)   :: NNN
  REAL           :: RSTTIME
  CHARACTER(100) :: STR
  LOGICAL FEXIST

  ! CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE THE ASCII FILE INSTEAD.
  INQUIRE(FILE='wq_wcrst.bin', EXIST=FEXIST)

  NCOL  = 0
  IF( .NOT. FEXIST )THEN
    ! *** FORMATTED
    WRITE(*,'(A)')' WQ: RESTART: WQ_WCRST.INP'
    OPEN(1,FILE='wq_wcrst.inp',STATUS='UNKNOWN')

    ! *** SKIP THE FIRST 2 LINES
    READ(1,999)
    READ(1,999)

    READ(1,'(A)') STR
    NCOL = INDEX(STR,'RADIATION')
    IF( NCOL > 0 )THEN
      NCOL = INDEX(STR,'=') + 1
      !READ(STR,*) VER,NNN,RSTTIME
      READ(STR(NCOL:100),* ) WQI1, WQI2 ,WQI3

      ! *** SKIP THE HEADER LINE
      READ(1,999)
    ENDIF

    DO LG=2,LA_Global
      DO K=1,KC
        READ(1,*) L, M, (WQV_Global(LG,K,NW),NW=1,NWQV)
      ENDDO
    ENDDO
      
  ELSE
    ! *** BINARY
    WRITE(*,'(A)')' WQ: RESTART: WQ_WCRST.BIN'
    OPEN(1,FILE='wq_wcrst.bin',STATUS='UNKNOWN',FORM='BINARY')
      
    READ(1) NREST, RSTTIME
    WRITE(0,911) NREST, RSTTIME
911 FORMAT(' READING BINARY WQ_WCRST.BIN FILE ...    NN, TIME = ', I7, F11.5)

    DO L=2,LA_Global
      DO KK=1,KC
        READ(1) LG, K
        READ(1) (WQV_Global(LG,K,NW),NW=1,NWQV)
      ENDDO
    ENDDO
  ENDIF
    
  CLOSE(1)
  
 90 FORMAT(2I5, 21E12.4)
999 FORMAT(1X)

END SUBROUTINE WQ_WCRST_IN

SUBROUTINE WQ_WCRST_OUT
  ! ** WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT WQ_WCRST.

  INTEGER :: NWQV0, L, K, NW, NHR, VER
  CHARACTER(2) :: SNHR

  if( process_id == master_id )then
    IF( ISGREGOR ==  1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = OUTDIR//'WQ_WCRST_'//TRIM(SDATE)//'_'//SNHR//'.OUT'
    ELSE
      RESFILE = OUTDIR//'WQ_WCRST.OUT'
    ENDIF
    OPEN(1,FILE=RESFILE,STATUS='REPLACE')

    ! *** BEGIN WQ_WCRST (ASCII FILE)
    VER = 8400
    WRITE(1,'(A)') 'C WQ_WCRST.INP - WQ INITIAL CONDITIONS'
    WRITE(1,101) VER, N, TIMEDAY
    WRITE(1,102) WQI1, WQI2, WQI3
    WRITE(1,103) WQCONSTIT(1:NWQV)
    
    DO L=2,LA_Global
      DO K=1,KC
        WRITE(1,104) L,K,(WQV_Global(L,K,NW),NW=1,NWQV)
      ENDDO
    ENDDO
    CLOSE(1)
  endif
  
101 FORMAT('CC  WQ RESTART FILE TIME STEP, TIME = ',2I10,F12.4)
102 FORMAT('CC  PREVIOUS DAYS SOLAR RADIATION   = ',3F12.5)
103 FORMAT('C   L    K',50A14)
104 FORMAT(2I5, 50E14.6)

END SUBROUTINE WQ_WCRST_OUT

SUBROUTINE WQRPEMRST_Out

  INTEGER :: L,NHR
  CHARACTER(2) :: SNHR

  ! **********************************************************************C
  ! ** RESTART OUTPUT
  if( process_id == master_id )then
    IF( ISGREGOR ==  1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = OUTDIR//'WQRPEMRST_'//TRIM(SDATE)//'_'//SNHR//'.OUT'
    ELSE
      RESFILE = OUTDIR//'WQRPEMRST_Out.OUT'
    ENDIF
    OPEN(1,FILE=RESFILE,STATUS='REPLACE')

    WRITE(1,110)TIMEDAY
    DO L=2,LA_Global
      IF( LMASKRPEM_Global(L) )THEN
        WRITE(1,111) L, WQRPS_Global(L), WQRPR_Global(L), WQRPE_Global(L), WQRPD_Global(L)
      ENDIF
    ENDDO
    CLOSE(1)
  endif
  
110 FORMAT(F12.4)
111 FORMAT(I6,6E13.5)

END SUBROUTINE

SUBROUTINE Setup_Continuation_Files(RESTARTF)
  ! ** EDIT EFDC.INP FOR RESTART
  ! ** COPY RESTART FILES INTO PROJECT FOLDER
  ! ** CHANGE THE FOLLOWING PARAMETERS:
  ! ** C2: ISRESTI =1 EFDC RUNS RESTART
  ! **     ICONTINUE=0 EFDC WRITES OUT FILES FROM BEGINING AS USUAL
  ! **     ICONTINUE=1 EFDC READS OUT FILES & WRITES THE NEXT OUTPUTS AT THE TRUE
  ! **                LOCATION RIGHT AFTER THE TBEGINC
  ! ** C7  NTC CHANGES TO THE TOTAL NUMBER OF DAYS IN EFDC.INP
  ! ** C8: TBEGIN STILL KEEPS THE PREVIOUS VALUE   IN EFDC.INP
  ! ** RESET ISCI(1:8),ISCO(1:8) & SOME PARAMETERS
  CHARACTER(*) :: RESTARTF
  LOGICAL(4) :: RESLOG

  INTEGER :: SLEN,L
  Character(6) :: copy_command
  CHARACTER*3  FLABEL*12
  CHARACTER(40) :: STR*200,RSTFILE1,RSTFILE2,RSTFILE3,RSTFILE4,RSTFILE5
#ifdef LINUX
  copy_command = 'cp'
#else
  copy_command = 'copy'
#endif
  ! ** COPY FILES
  IF( ISGREGOR ==  1 )THEN
    RESTARTF = ADJUSTL(RESTARTF)
    SLEN = LEN_TRIM(RESTARTF)
    DO L=1,SLEN
      IF( RESTARTF(L:L) == '_') EXIT
    ENDDO
    FLABEL = RESTARTF(L:SLEN)
  ELSE
    FLABEL = ''
  ENDIF

#ifdef LINUX
  RSTFILE1 = './'//OUTDIR//'RESTART'//TRIM(FLABEL)//'.OUT'
  RSTFILE2 = './'//OUTDIR//'RSTWD'//TRIM(FLABEL)//'.OUT'
  RSTFILE3 = './'//OUTDIR//'WQSDRST'//TRIM(FLABEL)//'.OUT'
  RSTFILE4 = './'//OUTDIR//'WQ_WCRST'//TRIM(FLABEL)//'.OUT'
  RSTFILE5 = './'//OUTDIR//'SEDBED_HOT'//TRIM(FLABEL)//'.OUT'
#else
  RSTFILE1 = OUTDIR//'RESTART'//TRIM(FLABEL)//'.OUT'
  RSTFILE2 = OUTDIR//'RSTWD'//TRIM(FLABEL)//'.OUT'
  RSTFILE3 = OUTDIR//'WQSDRST'//TRIM(FLABEL)//'.OUT'
  RSTFILE4 = OUTDIR//'WQ_WCRST'//TRIM(FLABEL)//'.OUT'
  RSTFILE5 = OUTDIR//'SEDBED_HOT'//TRIM(FLABEL)//'.OUT'
#endif

  STR =TRIM(copy_command)//' '//trim(RSTFILE1)//' restart.inp'
  WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: RESTART.INP'
  RESLOG = SYSTEMQQ(TRIM(STR))

  IF( ISDRY > 0 .AND. Restart_In_Ver < 1000 )THEN
    WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: RSTWD.INP'
    STR =TRIM(copy_command)//' '//trim(RSTFILE2)//' rstwd.inp'
    RESLOG = SYSTEMQQ(TRIM(STR))
  ENDIF

  IF( LSEDZLJ )THEN
    WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: SEDBED_HOT.SDF'
    STR =TRIM(copy_command)//' '//trim(RSTFILE5)//' sedbed_hot.sdf'
    RESLOG = SYSTEMQQ(TRIM(STR))
  ENDIF

  ! ** WQ: THIS MUST COME AFTER READING EFDC.INP **
  IF( ISTRAN(8) >= 1 )THEN
    ! ** wq_3dwc.jnp
    WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: WQ_WCRST.INP'
    STR =TRIM(copy_command)//' '//trim(RSTFILE4)//' wq_wcrst.inp'
    RESLOG = SYSTEMQQ(trim(str))

    IF( IWQBEN == 1 .AND. ISMICI == 2 )THEN
      ! ** wq_3dsd.jnp
      WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: WQSDRST.INP'
      STR =TRIM(copy_command)//' '//trim(RSTFILE3)//' wqsdrst.inp'
      RESLOG = SYSTEMQQ(TRIM(STR))
    ENDIF

  ENDIF

  CALL Restart_In(1)  !READ NRESTART, TBEGINC, Restart_In_Ver ONLY

END SUBROUTINE Setup_Continuation_Files

SUBROUTINE Gather_Restart_Arrays

  !---------------------------------------------------------------------------!
  ! *** Gathers values for all of RESTART variables
  !---------------------------------------------------------------------------!
  Use Mod_Assign_Loc_Glob_For_Write
  Use Mod_Map_Write_EE_Binary

  Implicit None

  ! *** Local variables
  Integer :: i, III, IERR, j, JJJ, K, L, LL, LG, NS
  Integer,Target, Allocatable, Dimension(:) :: IDRY
  Real,Target, Allocatable, Dimension(:)    :: ICE

  ! *** Clear out the pointer array
  j = 0

  ! *** 1D Section
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(HP,1), HP, size(HP_Global,1), HP_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(H1P,1), H1P, size(H1P_Global,1), H1P_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(HWQ,1), HWQ, size(HWQ_Global,1), HWQ_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(H2WQ,1), H2WQ, size(H2WQ_Global,1), H2WQ_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(BELV,1), BELV, size(BELV_Global,1), BELV_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(UHDYE,1), UHDYE, size(UHDYE_Global,1), UHDYE_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(UHDY1E,1), UHDY1E, size(UHDY1E_Global,1), UHDY1E_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(VHDXE,1), VHDXE, size(VHDXE_Global,1), VHDXE_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(VHDX1E,1), VHDX1E, size(VHDX1E_Global,1), VHDX1E_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(ISCDRY,1), ISCDRY, size(ISCDRY_Global,1), ISCDRY_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(NATDRY,1), NATDRY, size(NATDRY_Global,1), NATDRY_Global)

  ALLOCATE(IDRY(LCM))
  IDRY = 0
  DO L = 2,LA
    IF( .NOT. LMASKDRY(L) ) IDRY(L) = 1
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(IDRY,1), IDRY, size(IDRY_Global,1), IDRY_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SUB,1), SUB, size(SUB_Global,1), SUB_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SVB,1), SVB, size(SVB_Global,1), SVB_Global)

  ! *** 2D Section
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(U,1), size(U,2), U, size(U_Global,1), size(U_Global,2), U_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(U1,1), size(U1,2), U1, size(U1_Global,1), size(U1_Global,2), U1_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(V,1), size(V,2), V, size(V_Global,1), size(V_Global,2), V_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(V1,1), size(V1,2), V1, size(V1_Global,1), size(V1_Global,2), V1_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(QQ,1), size(QQ,2), QQ, size(QQ_Global,1), size(QQ_Global,2), QQ_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(QQ1,1), size(QQ1,2), QQ1, size(QQ1_Global,1), size(QQ1_Global,2), QQ1_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(QQL,1), size(QQL,2), QQL, size(QQL_Global,1), size(QQL_Global,2), QQL_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(QQL1,1), size(QQL1,2), QQL1, size(QQL1_Global,1), size(QQL1_Global,2), QQL1_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(DML,1), size(DML,2), DML, size(DML_Global,1), size(DML_Global,2), DML_Global)

  ! *** Conditional Arrays Section
  IF( ISTRAN(1) >= 1 .AND. ISCO(1) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SAL,1), size(SAL,2), SAL, size(SAL_Global,1), size(SAL_Global,2), SAL_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SAL1,1), size(SAL1,2), SAL1, size(SAL1_Global,1), size(SAL1_Global,2), SAL1_Global)
  ENDIF

  IF( ISTRAN(2) >= 1 .AND. ISCO(2) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TEM,1), size(TEM,2), TEM, size(TEM_Global,1), size(TEM_Global,2), TEM_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TEM1,1), size(TEM1,2), TEM1, size(TEM1_Global,1), size(TEM1_Global,2), TEM1_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TEMB,1), TEMB, size(TEMB_Global,1), TEMB_Global)

    IF( ISICE > 2 )THEN
      ALLOCATE(ICE(LCM))
      ICE = 0.0

      ! *** Get ice thickness total
      DO L = 2,LA
        ICE(L) = ICETHICK(L) + ICEVOL(L)*999.8426/RHOI*DXYIP(L)
        IF( ISICE == 4 ) ICE(L) = ICE(L) + SUM(FRAZILICE(L,KSZ(L):KC))
      ENDDO

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(ICE,1), ICE, size(ICETHICK_Global,1), ICETHICK_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(ICETEMP,1), ICETEMP, size(ICETEMP_Global,1), ICETEMP_Global)
    ENDIF

  ENDIF

  IF( ISTRAN(3) >= 1 .AND. ISCO(3) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(DYE,1),        size(DYE,2),        size(DYE,3),        DYE, &
                                      size(DYE_Global,1), size(DYE_Global,2), size(DYE_Global,3), DYE_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(DYE1,1), size(DYE1,2), size(DYE1,3), DYE1, &
                                      size(DYE1_Global,1), size(DYE1_Global,2), size(DYE1_Global,3), DYE1_Global)
  ENDIF

  IF( ISTRAN(4) >= 1 .AND. ISCO(4) == 1 )THEN
    ! *** TODO
  ENDIF

  IF( ISTRAN(5) >= 1 .AND. ISCO(5) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOX,1),        size(TOX,2),        size(TOX,3),        TOX, &
                                      size(TOX_Global,1), size(TOX_Global,2), size(TOX_Global,3), TOX_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOX1,1),        size(TOX1,2),        size(TOX1,3),        TOX1, &
                                      size(TOX1_Global,1), size(TOX1_Global,2), size(TOX1_Global,3), TOX1_Global)
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOXB,1),        size(TOXB,2),        size(TOXB,3),        TOXB, &
                                      size(TOXB_Global,1), size(TOXB_Global,2), size(TOXB_Global,3), TOXB_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOXB1,1),        size(TOXB1,2),        size(TOXB1,3),        TOXB1, &
                                      size(TOXB1_Global,1), size(TOXB1_Global,2), size(TOXB1_Global,3), TOXB1_Global)
  ENDIF

  IF( ISTRAN(6) >= 1 .AND. ISCO(6) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SED,1),        size(SED,2),        size(SED,3),        SED, &
                                      size(SED_Global,1), size(SED_Global,2), size(SED_Global,3), SED_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SED1,1),        size(SED1,2),        size(SED1,3),        SED1, &
                                      size(SED1_Global,1), size(SED1_Global,2), size(SED1_Global,3), SED1_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SEDB,1),        size(SEDB,2),        size(SEDB,3),        SEDB, &
                                      size(SEDB_Global,1), size(SEDB_Global,2), size(SEDB_Global,3), SEDB_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SEDB1,1),        size(SEDB1,2),        size(SEDB1,3),        SEDB1, &
                                      size(SEDB1_Global,1), size(SEDB1_Global,2), size(SEDB1_Global,3), SEDB1_Global)
  ENDIF

  IF( ISTRAN(7) >= 1 .AND. ISCO(7) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SND,1),        size(SND,2),        size(SND,3),        SND, &
                                      size(SND_Global,1), size(SND_Global,2), size(SND_Global,3), SND_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SND1,1),        size(SND1,2),        size(SND1,3),        SND1, &
                                      size(SND1_Global,1), size(SND1_Global,2), size(SND1_Global,3), SND1_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SNDB,1),        size(SNDB,2),        size(SNDB,3),        SNDB, &
                                      size(SNDB_Global,1), size(SNDB_Global,2), size(SNDB_Global,3), SNDB_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SNDB1,1),        size(SNDB1,2),        size(SNDB1,3),        SNDB1, &
                                      size(SNDB1_Global,1), size(SNDB1_Global,2), size(SNDB1_Global,3), SNDB1_Global)
  ENDIF

  IF( (ISTRAN(6) >= 1 .AND. ISCO(6) == 1 ) .OR. (ISTRAN(7) >= 1 .AND. ISCO(7) == 1 ) )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(HBED,1), size(HBED,2), HBED, size(HBED_Global,1), size(HBED_Global,2), HBED_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(HBED1,1), size(HBED1,2), HBED1, size(HBED1_Global,1), size(HBED1_Global,2), HBED1_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(VDRBED,1), size(VDRBED,2), VDRBED, size(VDRBED_Global,1), size(VDRBED_Global,2), VDRBED_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(VDRBED1,1), size(VDRBED1,2), VDRBED1, size(VDRBED1_Global,1), size(VDRBED1_Global,2), VDRBED1_Global)
  ENDIF

  IF( ISTRAN(8) >= 1 .AND. ISCO(8) == 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(WQV,1),        size(WQV,2),        size(WQV,3),        WQV, &
                                      size(WQV_Global,1), size(WQV_Global,2), size(WQV_Global,3), WQV_Global)
    
    IF( IWQBEN == 1 )THEN
      !CALL WQSDRST_OUT
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SMPON,1), size(SMPON,2), SMPON, size(SMPON_Global,1), size(SMPON_Global,2), SMPON_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SMPOP,1), size(SMPOP,2), SMPOP, size(SMPOP_Global,1), size(SMPOP_Global,2), SMPOP_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SMPOC,1), size(SMPOC,2), SMPOC, size(SMPOC_Global,1), size(SMPOC_Global,2), SMPOC_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SM1NH4,1), SM1NH4, size(SM1NH4_Global,1), SM1NH4_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SM2NH4,1), SM2NH4, size(SM2NH4_Global,1), SM2NH4_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SM2NO3,1), SM2NO3, size(SM2NO3_Global,1), SM2NO3_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SM2PO4,1), SM2PO4, size(SM2PO4_Global,1), SM2PO4_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SM2H2S,1), SM2H2S, size(SM2H2S_Global,1), SM2H2S_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SMPSI,1), SMPSI, size(SMPSI_Global,1), SMPSI_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SM2SI,1), SM2SI, size(SM2SI_Global,1), SM2SI_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SMBST,1), SMBST, size(SMBST_Global,1), SMBST_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SMT,1), SMT, size(SMT_Global,1), SMT_Global)
    ENDIF
    
    IF( ISRPEM > 0  )THEN
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(WQRPS,1), WQRPS, size(WQRPS_Global,1), WQRPS_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(WQRPR,1), WQRPR, size(WQRPR_Global,1), WQRPR_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(WQRPE,1), WQRPE, size(WQRPE_Global,1), WQRPE_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(WQRPD,1), WQRPD, size(WQRPD_Global,1), WQRPD_Global)
    ENDIF

  ENDIF

  ! *** Boundary conditions
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(QSUME,1), QSUME, size(QSUME_Global,1), QSUME_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(QSUM,1), size(QSUM,2), QSUM, size(QSUM_Global,1), size(QSUM_Global,2), QSUM_Global)

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  ! *** Build local open boundary arrays
  IF( NPBS_GL > 0 )THEN
    Call Communicate_Open_BC(NBBSM, NPBS_GL, NPBS, IPBS_GL, IPBS, JPBS_GL, JPBS, NLOS_Global, CLOS_Global, NLOS, CLOS)
  ENDIF
  IF( NPBE_GL > 0 )THEN
    Call Communicate_Open_BC(NBBEM, NPBE_GL, NPBE, IPBE_GL, IPBE, JPBE_GL, JPBE, NLOE_Global, CLOE_Global, NLOE, CLOE)
  ENDIF
  IF( NPBW_GL > 0 )THEN
    Call Communicate_Open_BC(NBBWM, NPBW_GL, NPBW, IPBW_GL, IPBW, JPBW_GL, JPBW, NLOW_Global, CLOW_Global, NLOW, CLOW)
  ENDIF
  IF( NPBN_GL > 0 )THEN
    Call Communicate_Open_BC(NBBNM, NPBN_GL, NPBN, IPBN_GL, IPBN, JPBN_GL, JPBN, NLON_Global, CLON_Global, NLON, CLON)
  ENDIF

  ! *** SEDZLJ
  IF( LSEDZLJ  )THEN
    Call Communicate_SEDZLJ
  ENDIF

END SUBROUTINE Gather_Restart_Arrays

SUBROUTINE Communicate_Open_BC(NBBM, NPB_GL, NPB, IPB_GL, IPB, JPB_GL, JPB, NLO_Global, CLO_Global, NLO, CLO)
  !------------------------------------------------------------------------------------!
  ! *** Communicates the open boundary concentration rate of change controls
  !------------------------------------------------------------------------------------!

  ! *** Declare passed variables
  Integer, Intent(inout) :: NLO(NBBM,KCM,NSTVM), NLO_Global(NBBM,KCM,NSTVM)
  Real,    Intent(inout) :: CLO(NBBM,KCM,NSTVM), CLO_Global(NBBM,KCM,NSTVM)
  Integer, Intent(in)    :: NBBM, NPB_GL, NPB
  Integer, Intent(in)    :: IPB_GL(NBBM), IPB(NBBM), JPB_GL(NBBM), JPB(NBBM)

  ! *** Local variables
  Integer :: i, III, j, JJJ, K, L, LL, LG, nDim, NS
  INTEGER :: IERR, status_message(MPI_Status_Size)

  nDim = NBBM*KCM*NSTVM
  ALLOCATE(I1D_Global(nDim))
  Call AllocateDSI(R1D_Global, nDim, 0.0)

  NLO_Global = -999      ! *** Initialize to missing
  I1D_Global  = -999     ! *** Initialize to missing
  R1D_Global  = -999.0   ! *** Initialize to missing

  ! *** Update global lists
  DO j = 1,NPB_GL
    DO I = 1,NPB
      ! *** Match local BC to Global BC
      if( IL2IG(IPB(I)) == IPB_GL(j) .and.  JL2JG(JPB(I)) == JPB_GL(j) )then
        DO NS = 1,NSTVM
          DO K = 1,KC
            NLO_Global(J,K,NS) = NLO(I,K,NS)
            CLO_Global(J,K,NS) = CLO(I,K,NS)
          ENDDO
        ENDDO
      endif
    ENDDO
  ENDDO

  if( process_id /= master_id )then
    ! *** Collapse to 1D
    LL = 0
    DO NS = 1,NSTVM
      DO K = 1,KC
        DO j = 1,NPB_GL
          LL = LL + 1
          I1D_Global(LL) = NLO_Global(J,K,NS)
          R1D_Global(LL) = CLO_Global(J,K,NS)
        ENDDO
      ENDDO
    ENDDO

    ! *** Send updated global lists
    CALL DSI_SEND(I1D_Global, nDim, master_id)
    CALL DSI_SEND(R1D_Global, nDim, master_id)
  else
    ! *** Master Process.  Gather the boundary lists and accumulate into the final global list
    DO I = 1,num_Processors-1
      CALL DSI_RECV(I1D_Global, nDim, I)
      CALL DSI_RECV(R1D_Global, nDim, I)

      ! *** Integrate arrays into global list
      LL = 0
      DO NS = 1,NSTVM
        DO K = 1,KC
          DO j = 1,NPB_GL
            LL = LL + 1
            IF( I1D_Global(LL) /= -999 )THEN
              NLO_Global(J,K,NS) = I1D_Global(LL)
              CLO_Global(J,K,NS) = R1D_Global(LL)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! *** QC - Ensure all of the values have been set
    DO NS = 1,NSTVM
      DO K = 1,KC
        DO j = 1,NPB_GL
          IF( NLO_Global(J,K,NS) == -999 )THEN
            Call STOPP('Bad open BC communication: Missing values')
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  endif

  DEALLOCATE(I1D_Global)
  DEALLOCATE(R1D_Global)

END SUBROUTINE Communicate_Open_BC

Subroutine Communicate_SEDZLJ
  !------------------------------------------------------------------------------------!
  ! *** Communicates the open boundary concentration rate of change controls
  !------------------------------------------------------------------------------------!

  Use Mod_Assign_Loc_Glob_For_Write
  Use Mod_Map_Write_EE_Binary

  Implicit None

  ! *** Local
  Integer :: i, j, L, LG, k, nn
  Integer,   Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_ActLay
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_TSED
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_BULKDENS
  Real(rkd), Target, Allocatable, Dimension(:,:,:) :: Reverse_Temp_3D

  Integer,   Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_ActLay
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_TSED
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_BULKDENS
  Real(rkd), Target, Allocatable, Dimension(:,:,:) :: Gl_Reverse_Temp_3D

  Allocate(Reverse_Temp_2D_ActLay(LCM,KB))
  Allocate(Reverse_Temp_2D_TSED(LCM,KB))
  Allocate(Reverse_Temp_2D_BULKDENS(LCM,KB))
  Allocate(Reverse_Temp_3D(LCM,KB,NSCM))

  Allocate(GL_Reverse_Temp_2D_ActLay(LCM_Global,KB))
  Allocate(Gl_Reverse_Temp_2D_TSED(LCM_Global,KB))
  Allocate(Gl_Reverse_Temp_2D_BULKDENS(LCM_Global,KB))
  Allocate(GL_Reverse_Temp_3D(LCM_Global,KB,NSCM))

  Reverse_Temp_2D_ActLay      = 0
  Reverse_Temp_2D_TSED        = 0.0
  Reverse_Temp_2D_BULKDENS    = 0.0
  Reverse_Temp_3D             = 0.0
  Gl_Reverse_Temp_2D_ActLay   = 0
  Gl_Reverse_Temp_2D_TSED     = 0.0
  Gl_Reverse_Temp_2D_BULKDENS = 0.0
  Gl_Reverse_Temp_3D          = 0.0

  j = 0

  ! ********************************************************************************************************************
  ! *** Integer

  ! *** Reverse Order of Array (TSED)
  DO L=2,LA
    DO K=1,KB
      Reverse_Temp_2D_ActLay(L,K) = LAYERACTIVE(K,L)
    ENDDO
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_ActLay,1), size(Reverse_Temp_2D_ActLay,2), Reverse_Temp_2D_ActLay, &
    size(Gl_Reverse_Temp_2D_ActLay,1),size(Gl_Reverse_Temp_2D_ActLay,2), Gl_Reverse_Temp_2D_ActLay)

  ! ********************************************************************************************************************
  ! *** REAL*8

  ! *** Reverse Order of Array (TSED)
  DO L=2,LA
    DO K=1,KB
      Reverse_Temp_2D_TSED(L,K) = TSED(K,L)
    ENDDO
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_TSED,1), size(Reverse_Temp_2D_TSED,2), Reverse_Temp_2D_TSED, &
    size(Gl_Reverse_Temp_2D_TSED,1),size(Gl_Reverse_Temp_2D_TSED,2), Gl_Reverse_Temp_2D_TSED)

  ! *** Reverse Order of Array (BULKDENS)
  DO L=2,LA
    DO K=1,KB
      Reverse_Temp_2D_BULKDENS(L,K) = BULKDENS(K,L)
    ENDDO
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_BULKDENS,1), size(Reverse_Temp_2D_BULKDENS,2), Reverse_Temp_2D_BULKDENS, &
    size(Gl_Reverse_Temp_2D_BULKDENS,1),size(Gl_Reverse_Temp_2D_BULKDENS,2), Gl_Reverse_Temp_2D_BULKDENS)

  ! *** Reverse Order of Array (PERSED)
  DO L=2,LA
    DO K=1,KB
      DO nn=1,NSCM
        Reverse_Temp_3D(L,K,nn) = PERSED(nn,K,L)
      ENDDO
    ENDDO
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_3D,1), size(Reverse_Temp_3D,2), size(Reverse_Temp_3D,3), Reverse_Temp_3D, &
    size(Gl_Reverse_Temp_3D,1),size(Gl_Reverse_Temp_3D,2),size(Gl_Reverse_Temp_3D,3), Gl_Reverse_Temp_3D)


  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(D50AVG,1), D50AVG, size(D50AVG_Global,1), D50AVG_Global)

  IF( ICALC_BL > 0 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(CBL,1), size(CBL,2), CBL, size(CBL_Global,1),size(CBL_Global,2), CBL_Global)

    IF( ISTRAN(5) > 0 )THEN
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(CBLTOX,1), size(CBLTOX,2), CBLTOX, size(CBLTOX_Global,1),size(CBLTOX_Global,2), CBLTOX_Global)
    ENDIF
  ENDIF

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  ! *** Reverse global array order
  DO LG=2,LA_Global
    Do k = 1, KB
      LAYERACTIVE_Global(K,LG) = Gl_Reverse_Temp_2D_ActLay(LG,K)
    End do
  End do

  ! *** Reverse global array order
  DO LG=2,LA_Global
    Do k = 1, KB
      TSED_Global(K,LG) = Gl_Reverse_Temp_2D_TSED(LG,K)
    End do
  End do

  ! *** Reverse global array order
  DO LG = 2,LA_Global
    Do k = 1, KB
      BULKDENS_Global(K,LG) = Gl_Reverse_Temp_2D_BULKDENS(LG,K)
    End do
  End do

  ! *** Reverse global array order
  DO LG = 2,LA_Global
    DO nn = 1, NSCM
      Do k = 1, KB
        PERSED_Global(nn,k,LG) = Gl_Reverse_Temp_3D(LG,K,nn)
      End do
    End do
  End do

  Deallocate(Reverse_Temp_2D_ActLay)
  Deallocate(Reverse_Temp_2D_TSED)
  Deallocate(Reverse_Temp_2D_BULKDENS)
  Deallocate(Reverse_Temp_3D)
  Deallocate(Gl_Reverse_Temp_2D_ActLay)
  Deallocate(Gl_Reverse_Temp_2D_TSED)
  Deallocate(Gl_Reverse_Temp_2D_BULKDENS)
  Deallocate(Gl_Reverse_Temp_3D)

End Subroutine Communicate_SEDZLJ

END MODULE
