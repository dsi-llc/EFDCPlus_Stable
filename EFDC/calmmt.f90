! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CALMMT

  ! CHANGE RECORD
  ! **  SUBROUTINE CALMMTF CALCULATES THE MEAN MASS TRANSPORT FIELD

  USE GLOBAL
  Use Variables_MPI

  INTEGER      :: L, K, NSN, NT, NS, NMD, NWR, LN, LS, LSW, LL, LT, I, J, ITMP, NSC, LE, LW, MD, IFLAG
  REAL         :: UTMP1, VTMP1, UTMP, VTMP, FLTWT, TMPVAL, HPLW, HPLS, HPLSW, HMC
  REAL, SAVE   :: MMTCNT

  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: VPX
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: VPY
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: VPZ

  ! *** ALLOCATE ARRAYS
  IF( .NOT. ALLOCATED(VPX) )THEN
    ALLOCATE(VPX(LCM,0:KCM))
    ALLOCATE(VPY(LCM,0:KCM))
    ALLOCATE(VPZ(LCM,KCM))
    VPX = 0.0
    VPY = 0.0
    VPZ = 0.0

    MMTCNT = 0.
    NMMT = 1

    IFLAG = 0
    IF( RESSTEP > TIDALP ) IFLAG = 1
    CALL INITIALIZE_RESIDUALS(IFLAG)
    
    ! *** Initiate and write the first ouput to WASP linkage file 
#ifdef WASPOUT
    IF( ISWASP == 17 ) CALL WASP8HYDRO
#endif
  ENDIF

  ! *** Select time step depending on whether a static or dynamic time step is used
  IF( ISDYNSTP == 0 )THEN
    DELT = DT
  ELSE
    DELT = DTDYN
  ENDIF

  ! **  INITIALIZE CE-QUAL-ICM INTERFACE
  IF( ISICM >= 1 .AND. JSWASP == 1 ) CALL CEQICM

  ! *** ACCUMULATE MEAN MASS TRANPORT
  DO L = 2,LA
    LN = LNC(L)
    LE = LEC(L)
    HLPF(L)     = HLPF(L) + HP(L)*DELT
    QSUMELPF(L) = QSUMELPF(L) + QSUME(L)*DELT
    UTMP1    = 0.5*(UHDYE(LE) + UHDYE(L))/(DYP(L)*HP(L))
    VTMP1    = 0.5*(VHDXE(LN) + VHDXE(L))/(DXP(L)*HP(L))
    UTMP     = CUE(L)*UTMP1 + CVE(L)*VTMP1
    VTMP     = CUN(L)*UTMP1 + CVN(L)*VTMP1
    UELPF(L) = UELPF(L) + UTMP*DELT
    VELPF(L) = VELPF(L) + VTMP*DELT
    RAINLPF(L) = RAINLPF(L) + DXYP(L)*RAINT(L)*DELT
  ENDDO

  IF( ISGWIE == 0 )THEN
    DO L = 2,LA
      EVPSLPF(L) = EVPSLPF(L) + DXYP(L)*EVAPT(L)*DELT
    ENDDO
  ELSE
    DO L = 2,LA
      EVPSLPF(L) = EVPSLPF(L) + EVAPSW(L)*DELT         ! *** m3
      EVPGLPF(L) = EVPGLPF(L) + EVAPGW(L)*DELT         ! *** m3
      RINFLPF(L) = RINFLPF(L) - QGW(L)*DELT            ! *** m3
      GWLPF(L)   = GWLPF(L) + AGWELV(L)*DELT           ! *** m3
    ENDDO
  ENDIF

  IF( ISTRAN(5) > 0 )THEN
    DO NT = 1,NTOX
      DO K = 1,KB
        DO L = 2,LA
          TOXBLPF(L,K,NT) = TOXBLPF(L,K,NT) + TOXB(L,K,NT)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(6) > 0 )THEN
    DO NSC = 1,NSED
      DO K = 1,KB
        DO L = 2,LA
          SEDBLPF(L,K,NSC) = SEDBLPF(L,K,NSC) + SEDB(L,K,NSC)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    DO NSN = 1,NSND
      DO K = 1,KB
        DO L = 2,LA
          SNDBLPF(L,K,NSN) = SNDBLPF(L,K,NSN) + SNDB(L,K,NSN)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISWASP == 99 .OR. ISICM >= 1 )THEN
    ! *** ICM OR RCA LINKAGE
    DO K = 1,KS
      DO L = 2,LA
        ABEFF(L,K) = ABEFF(L,K) + AB(L,K)*(SAL(L,K+1)-SAL(L,K))*DELT
        ABLPF(L,K) = ABLPF(L,K) + (AB(L,K)*HP(L))*DELT
        WLPF(L,K)  = WLPF(L,K)  + W(L,K)*DELT
      ENDDO
    ENDDO

    IF( RESSTEP > TIDALP )THEN
      ! DELME - NOT TESTED
      DO K = 1,KS
        DO L = 2,LA
          WIRT(L,K)  = WIRT(L,K)  + DT*W(L,K)
          WTLPF(L,K) = WTLPF(L,K) + DT*(FLOAT(NMMT)-0.5)*W(L,K)
        ENDDO
      ENDDO
    ENDIF
  ELSE
    ! *** WASP LINKAGE
    DO K = 1,KS
      DO L = 2,LA
        ABEFF(L,K) = ABEFF(L,K) + AB(L,K)*(SAL(L,K+1)-SAL(L,K))*DELT
        ABLPF(L,K) = ABLPF(L,K) + AB(L,K)*DELT
        WLPF(L,K)  = WLPF(L,K)  + W(L,K)*DELT
      ENDDO
    ENDDO

    IF( RESSTEP > TIDALP )THEN
      ! DELME - NOT TESTED
      DO K = 1,KS
        DO L = 2,LA
          WIRT(L,K)  = WIRT(L,K)  + DT*W(L,K)
          WTLPF(L,K) = WTLPF(L,K) + DT*(FLOAT(NMMT)-0.5)*W(L,K)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  DO K = 1,KC
    DO L = 2,LA
      LS = LSC(L)
      LW = LWC(L)
      AHULPF(L,K) = AHULPF(L,K) + 0.5*(AH(L,K) + AH(LW,K))*DELT
      AHVLPF(L,K) = AHVLPF(L,K) + 0.5*(AH(L,K) + AH(LS,K))*DELT
      SALLPF(L,K) = SALLPF(L,K) + SAL(L,K)*DELT
      TEMLPF(L,K) = TEMLPF(L,K) + TEM(L,K)*DELT
      SFLLPF(L,K) = SFLLPF(L,K) + SFL(L,K)*DELT
      UHLPF(L,K)  = UHLPF(L,K) + UHDY2(L,K)/DYU(L)*DELT
      VHLPF(L,K)  = VHLPF(L,K) + VHDX2(L,K)/DXV(L)*DELT
      QSUMLPF(L,K) = QSUMLPF(L,K) + QSUM(L,K)*DELT
    ENDDO
  ENDDO

  IF( RESSTEP > TIDALP )THEN
    ! *** DELME - NOT TESTED
    DO K = 1,KC
      DO L = 2,LA
        LS = LSC(L)
        LW = LWC(L)
        ULPF(L,K)   = ULPF(L,K) + U(L,K)*DELT
        VLPF(L,K)   = VLPF(L,K) + V(L,K)*DELT

        UIRT(L,K)   = UIRT(L,K) + U(L,K)*DELT
        UTLPF(L,K)  = UTLPF(L,K) + (FLOAT(NMMT)-0.5)*U(L,K)

        VIRT(L,K)   = VIRT(L,K) + V(L,K)*DELT
        VTLPF(L,K)  = VTLPF(L,K) + (FLOAT(NMMT)-0.5)*V(L,K)
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(3) > 0 )THEN
    DO MD = 1,NDYE
      DO K = 1,KC
        DO L = 1,LC
          DYELPF(L,K,MD) = DYELPF(L,K,MD) + DYE(L,K,MD)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(5) > 0 )THEN
    DO NT = 1,NTOX
      DO K = 1,KC
        DO L = 2,LA
          TOXLPF(L,K,NT) = TOXLPF(L,K,NT) + TOX(L,K,NT)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(6) > 0 )THEN
    DO NSC = 1,NSED
      DO K = 1,KC
        DO L = 2,LA
          SEDLPF(L,K,NSC) = SEDLPF(L,K,NSC) + SED(L,K,NSC)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    DO NSN = 1,NSND
      DO K = 1,KC
        DO L = 2,LA
          SNDLPF(L,K,NSN) = SNDLPF(L,K,NSN) + SND(L,K,NSN)*DELT
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(5) > 0 )THEN
    DO NT = 1,NTOX
      DO NS = 1,NSED+NSND
        DO K = 1,KC
          DO L = 1,LC
            TXPFLPF(L,K,NS,NT) = TXPFLPF(L,K,NS,NT) + TOXPFW(L,K,NS,NT)*DELT
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! *** SOURCE SINK
  DO NS = 1,NQSER
    DO K = 1,KC
      QSRTLPP(K,NS) = QSRTLPP(K,NS) + MAX(QSERT(K,NS),0.)*DELT
      QSRTLPN(K,NS) = QSRTLPN(K,NS) + MIN(QSERT(K,NS),0.)*DELT
    ENDDO
  ENDDO
  DO NS = 1,NQCTL
    DO K = 1,KC
      QCTLTLP(K,NS) = QCTLTLP(K,NS) + QCTLT(K,NS,1)*DELT
    ENDDO
  ENDDO
  DO NMD = 1,MDCHH
    QCHNULP(NMD) = QCHNULP(NMD) + QCHANU(NMD)*DELT
    QCHNVLP(NMD) = QCHNVLP(NMD) + QCHANV(NMD)*DELT
  ENDDO
  DO NWR = 1,NQWR
    QWRSERTLP(NWR) = QWRSERTLP(NWR) + QWRSERT(NWR)*DELT
  ENDDO

  ! ***
  IF( RESSTEP > TIDALP )THEN
    ! DELME - NOT TESTED
    DO K = 1,KS
      DO L = 2,LA
        LS = LSC(L)
        LW = LWC(L)
        VPX(L,K) = VPX(L,K) + 0.25*(V(L,K+1) + V(L,K))*(WIRT(L,K) + WIRT(LS,K))
        VPY(L,K) = VPY(L,K) + 0.25*(W(L,K)   + W(LW,K))*(UIRT(L,K+1) + UIRT(L,K))
      ENDDO
    ENDDO
    DO K = 1,KC
      DO L = 2,LA
        LS = LSC(L)
        LW = LWC(L)
        VPZ(L,K) = VPZ(L,K) + 0.25*(U(L,K) + U(LS,K))*(VIRT(L,K) + VIRT(LW,K))
      ENDDO
    ENDDO
  ENDIF
  ! *****************************************************************************************************************

  ! *** ACCUMULATE THE DELTA T
  MMTCNT = MMTCNT + DELT
  IF( ISICM >= 1 ) MMTCNT = MMTCNT + DELT

  ! *****************************************************************************************************************
  ! **  CHECK FOR END OF FILTER
  IF( TIMEDAY >= WASPTIME(NTIME) )THEN

    ! **  COMPLETE THE FILTERING
    FLTWT = 1./MMTCNT
    IF( ISICM >= 1 ) FLTWT = 2.*FLTWT

    ! *** FINALIZE FOR CASE WHEN MEAN MASS TRANPORT INTERVAL IS < REFERENCE PERIOD
    DO L = 2,LA
      HLPF(L)     = FLTWT*HLPF(L)
      QSUMELPF(L) = FLTWT*QSUMELPF(L)
      UELPF(L)    = FLTWT*UELPF(L)
      VELPF(L)    = FLTWT*VELPF(L)
      RAINLPF(L)  = FLTWT*RAINLPF(L)
    ENDDO

    IF( ISGWIE == 0 )THEN
      DO L = 2,LA
        EVPSLPF(L) = FLTWT*EVPSLPF(L)
      ENDDO
    ELSE
      DO L = 2,LA
        EVPSLPF(L) = FLTWT*EVPSLPF(L)
        EVPGLPF(L) = FLTWT*EVPGLPF(L)
        RINFLPF(L) = FLTWT*RINFLPF(L)
        GWLPF(L)   = FLTWT*GWLPF(L)
      ENDDO
    ENDIF

    IF( ISTRAN(5) > 0 )THEN
      DO NT = 1,NTOX
        DO K = 1,KB
          DO L = 2,LA
            TOXBLPF(L,K,NT) = TOXBLPF(L,K,NT)*FLTWT
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( ISTRAN(6) > 0 )THEN
      DO NSC = 1,NSED
        DO K = 1,KB
          DO L = 2,LA
            SEDBLPF(L,K,NSC) = SEDBLPF(L,K,NSC)*FLTWT
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( ISTRAN(7) > 0 )THEN
      DO NSN = 1,NSND
        DO K = 1,KB
          DO L = 2,LA
            SNDBLPF(L,K,NSN) = SNDBLPF(L,K,NSN)*FLTWT
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DO K = 1,KS
      DO L = 2,LA
        ABEFF(L,K) = FLTWT*ABEFF(L,K)
        ABLPF(L,K) = FLTWT*ABLPF(L,K)
        WLPF(L,K)  = FLTWT*WLPF(L,K)
      ENDDO
    ENDDO

    DO K = 1,KC
      DO L = 2,LA
        AHULPF(L,K)  = FLTWT*AHULPF(L,K)
        AHVLPF(L,K)  = FLTWT*AHVLPF(L,K)
        SALLPF(L,K)  = FLTWT*SALLPF(L,K)
        TEMLPF(L,K)  = FLTWT*TEMLPF(L,K)
        SFLLPF(L,K)  = FLTWT*SFLLPF(L,K)
        UHLPF(L,K)   = FLTWT*UHLPF(L,K)
        VHLPF(L,K)   = FLTWT*VHLPF(L,K)
        QSUMLPF(L,K) = FLTWT*QSUMLPF(L,K)
      ENDDO
    ENDDO

    IF( ISTRAN(3) > 0 )THEN
      DO MD = 1,NDYE
        DO K = 1,KC
          DO L = 1,LC
            DYELPF(L,K,MD) = DYELPF(L,K,MD)*FLTWT
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( ISTRAN(6) > 0 )THEN
      DO NSC = 1,NSED
        DO K = 1,KC
          DO L = 2,LA
            SEDLPF(L,K,NSC) = SEDLPF(L,K,NSC)*FLTWT
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( ISTRAN(7) > 0 )THEN
      DO NSN = 1,NSND
        DO K = 1,KC
          DO L = 2,LA
            SNDLPF(L,K,NSN) = SNDLPF(L,K,NSN)*FLTWT
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( ISTRAN(5) > 0 )THEN
      DO NT = 1,NTOX
        DO K = 1,KC
          DO L = 2,LA
            TOXLPF(L,K,NT) = TOXLPF(L,K,NT)*FLTWT
          ENDDO
        ENDDO
      ENDDO
      DO NT = 1,NTOX
        DO NS = 1,NSED+NSND
          DO K = 1,KC
            DO L = 1,LC
              TXPFLPF(L,K,NS,NT) = TXPFLPF(L,K,NS,NT)*FLTWT
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! *** SOURCES & SINKS
    DO NS = 1,NQSER
      DO K = 1,KC
        QSRTLPP(K,NS) = FLTWT*QSRTLPP(K,NS)
        QSRTLPN(K,NS) = FLTWT*QSRTLPN(K,NS)
      ENDDO
    ENDDO
    DO NS = 1,NQCTL
      DO K = 1,KC
        QCTLTLP(K,NS) = FLTWT*QCTLTLP(K,NS)
      ENDDO
    ENDDO
    DO NMD = 1,MDCHH
      QCHNULP(NMD) = FLTWT*QCHNULP(NMD)
      QCHNVLP(NMD) = FLTWT*QCHNVLP(NMD)
    ENDDO
    DO NWR = 1,NQWR
      QWRSERTLP(NWR) = FLTWT*QWRSERTLP(NWR)
    ENDDO

    IF( RESSTEP > TIDALP )THEN
      ! DELME - NOT TESTED

      ! *** POTENTIAL TRANSPORT VELOCITY
      DO K = 1,KS
        DO L = 2,LA
          LS = LSC(L)
          LW = LWC(L)
          VPX(L,K) = VPX(L,K) - 0.25*(VTLPF(L,K+1) + VTLPF(L,K)) *(WLPF(L,K)   + WLPF(LS,K))
          VPY(L,K) = VPY(L,K) - 0.25*(WTLPF(L,K)   + WTLPF(LW,K))*(ULPF(L,K+1) + ULPF(L,K))
        ENDDO
      ENDDO
      DO K = 1,KC
        DO L = 2,LA
          LS = LSC(L)
          LSW = LSWC(L)
          LW = LWC(L)
          VPZ(L,K) = VPZ(L,K) - 0.25*(UTLPF(L,K) + UTLPF(LS,K))*(VLPF(L,K) + VLPF(LW,K))
          TMPVAL = 1.+SUB(L)+SVB(L)+SUB(L)*SVB(L)
          HPLW   = SUB(L)*HP(LW)
          HPLS   = SVB(L)*HP(LS)
          HPLSW  = SUB(L)*HP(LSW)
          HMC    = (HP(L)+HPLW+HPLS+HPLSW)/TMPVAL
          VPZ(L,K) = VPZ(L,K)*HMC*SUB(L)*SUB(LS)*SVB(L)*SVB(LW)
        ENDDO
      ENDDO
      DO K = 1,KC
        DO L = 2,LA
          LS = LSC(L)
          LN = LNC(L)
          LE = LEC(L)
          UVPT(L,K) = (VPZ(LN,K)-VPZ(L,K))/DYU(L) - DZIC(L,K)*(VPY(L,K)-VPY(L,K-1))      ! DELME DZIC???
          VVPT(L,K) = DZIC(L,K)*(VPX(L,K)-VPX(L,K-1)) - (VPZ(LE,K)-VPZ(L,K))/DXV(L)
        ENDDO
      ENDDO
      DO K = 1,KS
        DO L = 2,LA
          LS = LSC(L)
          LN = LNC(L)
          LE = LEC(L)
          WVPT(L,K) = (VPY(LE,K)-VPY(L,K))/DXP(L)-(VPX(LN,K)-VPX(L,K) )/DYP(L)
        ENDDO
      ENDDO
    ENDIF

    ! *** RESIDUAL FLOWS ACROSS OPEN BCs
    QXW = 0.
    QXWVP = 0.
    DO K = 1,KC
      DO LL = 1,NPBW
        L = LPBW(LL)
        LE = LEC(L)
        QXW   = QXW   + UHLPF(LE,K)*DYU(LE)
        QXWVP = QXWVP + UVPT(LE,K)*DZC(L,K)*DYU(LE)
      ENDDO
    ENDDO

    QXE = 0.
    QXEVP = 0.
    DO K = 1,KC
      DO LL = 1,NPBE
        L = LPBE(LL)
        QXE = QXE+UHLPF(L,K)*DYU(L)
        QXEVP = QXEVP+UVPT(L,K)*DZC(L,K)*DYU(L)
      ENDDO
    ENDDO

    QYS = 0.
    QYSVP = 0.
    DO K = 1,KC
      DO LL = 1,NPBS
        L = LPBS(LL)
        LN = LNC(L)
        QYS = QYS+VHLPF(LN,K)*DXV(LN)
        QYSVP = QYSVP+VVPT(LN,K)*DZC(L,K)*DXV(LN)
      ENDDO
    ENDDO

    QYN = 0.
    QYNVP = 0.
    DO K = 1,KC
      DO LL = 1,NPBN
        L = LPBN(LL)
        QYN   = QYN   + VHLPF(L,K)*DXV(L)
        QYNVP = QYNVP + VVPT(L,K)*DZC(L,K)*DXV(L)
      ENDDO
    ENDDO

    ! *****************************************************************************************************************
    ! *** Process output options
    ! **  OUTPUT RESIDUAL TRANSPORT TO FILE RESTRAN.OUT
    IF( ISSSMMT == 1 .OR. ISWASP > 0 .OR. (ISSSMMT == 2 .AND. TIMEDAY >= FLOAT(NTC)*TIDALP) )THEN

      if( process_id == master_id )THEN
        IF( ISRESTR == 1 )THEN
          IF( JSRESTR == 1 )THEN
            OPEN(98,FILE = OUTDIR//'RESTRAN.OUT',STATUS = 'UNKNOWN')
            CLOSE(98,STATUS = 'DELETE')
            OPEN(98,FILE = OUTDIR//'RESTRAN.OUT',STATUS = 'UNKNOWN')
            JSRESTR = 0
          ELSE
            OPEN(98,FILE = OUTDIR//'RESTRAN.OUT',POSITION = 'APPEND',STATUS = 'UNKNOWN')
          ENDIF
          IF( RESSTEP < TIDALP )THEN
            DO LT = 2,LALT
              I = ILLT(LT)
              J = JLLT(LT)
              L = LIJ(I,J)
              WRITE(98,907)HP(L),HLPF(L),QSUMELPF(L)
              WRITE(98,907)(UHLPF(L,K),K = 1,KC)
              WRITE(98,907)(VHLPF(L,K),K = 1,KC)
              WRITE(98,907)(AHULPF(L,K),K = 1,KC)
              WRITE(98,907)(AHVLPF(L,K),K = 1,KC)
              WRITE(98,907)(SALLPF(L,K),K = 1,KC)
              WRITE(98,907)(ABLPF(L,K),K = 1,KS)
              WRITE(98,907)(ABEFF(L,K),K = 1,KS)
            ENDDO
          ELSE
            DO LT = 2,LALT
              I = ILLT(LT)
              J = JLLT(LT)
              L = LIJ(I,J)
              WRITE(98,907)HP(L),HLPF(L),QSUMELPF(L)
              WRITE(98,907)(UHLPF(L,K),K = 1,KC)
              WRITE(98,907)(VHLPF(L,K),K = 1,KC)
              WRITE(98,907)(VPZ(L,K),K = 1,KC)
              WRITE(98,907)(AHULPF(L,K),K = 1,KC)
              WRITE(98,907)(AHVLPF(L,K),K = 1,KC)
              WRITE(98,907)(SALLPF(L,K),K = 1,KC)
              WRITE(98,907)(VPX(L,K),K = 1,KS)
              WRITE(98,907)(VPY(L,K),K = 1,KS)
              WRITE(98,907)(ABLPF(L,K),K = 1,KS)
            ENDDO
          ENDIF
          CLOSE(98)
        ENDIF
907     FORMAT(12E12.4)
#ifdef WASPOUT
        ! **  OUTPUT TO WASP COMPATIABLE FILES
        !IF( ISWASP == 4  ) CALL WASP4             delme - untested
        !IF( ISWASP == 5  ) CALL WASP5
        !IF( ISWASP == 6  ) CALL WASP6
        IF( ISWASP == 7  ) CALL WASP7
        IF( ISWASP == 17 ) CALL WASP8HYDRO
        !IF( ISRCA  >= 1  ) CALL RCAHQ
        !IF( ISICM  >= 1  ) CALL CEQICM
#endif
        ! **  WRITE GRAPHICS FILES FOR RESIDUAL VARIABLES
        ! **  RESIDUAL SALINITY CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
        IF( ISRSPH(1) == 1 .AND. ISTRAN(1) >= 1 )THEN
          CALL RSALPLTH(1,SALLPF)
        ENDIF
        IF( ISRSPH(2) == 1 .AND. ISTRAN(2) >= 1 )THEN
          CALL RSALPLTH(2,TEMLPF)
        ENDIF
        IF( ISRSPH(3) == 1 .AND. ISTRAN(3) >= 1 )THEN
          DO MD = 1,NDYE
            CALL RSALPLTH(3,DYELPF(1,1,MD))
          ENDDO
        ENDIF
        IF( ISRSPH(4) == 1 .AND. ISTRAN(4) >= 1 )THEN
          CALL RSALPLTH(4,SFLLPF)
        ENDIF

        IF( ISTRAN(6) > 0 )THEN
          DO K = 2,KB
            DO L = 2,LA
              SEDBTLPF(L,K) = 0.
            ENDDO
          ENDDO
          DO K = 1,KC
            DO L = 2,LA
              SEDTLPF(L,K) = 0.
            ENDDO
          ENDDO
        ENDIF

        IF( ISTRAN(7) > 0 )THEN
          DO K = 2,KB
            DO L = 2,LA
              SNDBTLPF(L,K) = 0.
            ENDDO
          ENDDO
          DO K = 1,KC
            DO L = 2,LA
              SNDTLPF(L,K) = 0.
            ENDDO
          ENDDO
        ENDIF

        IF( ISRSPH(5) == 1 .AND. ISTRAN(5) >= 1 )THEN
          DO K = 1,KC
            DO L = 2,LA
              TVAR1S(L,K) = TOXLPF(L,K,1)
            ENDDO
          ENDDO
          DO NT = 1,NTOX
            CALL RSALPLTH(5,TVAR1S)   ! *** DELME - WHY IS TOXLPF NOT PASSED WITHOUT COPYING THE ARRAY?
          ENDDO

        ENDIF

        IF( ISTRAN(6) > 0 )THEN
          DO NS = 1,NSED
            DO K = 1,KB
              DO L = 2,LA
                SEDBTLPF(L,K) = SEDBTLPF(L,K) + SEDBLPF(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
          DO NS = 1,NSED
            DO K = 1,KC
              DO L = 2,LA
                SEDTLPF(L,K) = SEDTLPF(L,K) + SEDLPF(L,K,NS)
              ENDDO
            ENDDO
          ENDDO

          IF( ISRSPH(6) == 1 .AND. ISTRAN(6) >= 1 )THEN
            DO NSC = 1,NSED
              CALL RSALPLTH(6,SEDTLPF)
            ENDDO
          ENDIF
        ENDIF

        IF( ISTRAN(7) > 0 )THEN
          DO NS = 1,NSND
            DO K = 1,KB
              DO L = 2,LA
                SNDBTLPF(L,K) = SNDBTLPF(L,K)+SNDBLPF(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
          DO NS = 1,NSND
            DO K = 1,KC
              DO L = 2,LA
                SNDTLPF(L,K) = SNDTLPF(L,K) + SNDLPF(L,K,NS)
              ENDDO
            ENDDO
          ENDDO

          IF( ISRSPH(7) == 1 .AND. ISTRAN(7) >= 1 )THEN
            DO NSN = 1,NSND
              CALL RSALPLTH(7,SNDTLPF)    ! ***  DELME THIS WRITE THE SAME ARRAY NSND TIMES?
            ENDDO
          ENDIF
        ENDIF

        ! **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:    RVELPLTH
        IF( ISRVPH >= 1 ) CALL RVELPLTH

        ! **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:  RSURFPLT
        IF( ISRPPH == 1 ) CALL RSURFPLT

        ! **  RESIDUAL SALINITY AND VERTICAL MASS DIFFUSIVITY CONTOURING IN 3 VERTICAL PLANES:  RSALPLTV
        DO ITMP = 1,7
          IF( ISRSPV(ITMP) >= 1 ) CALL RSALPLTV(ITMP)
        ENDDO

        ! **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND
        ! **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:  RVELPLTV
        IF( ISRVPV >= 1 ) CALL RVELPLTV

        ! **  RESIDUAL 3D SCALAR AND VECTOR OUTPUT FILES
        IF( ISR3DO >= 1 ) CALL ROUT3D
      end if   ! *** End of master process

      IFLAG = 0
      IF( RESSTEP > TIDALP ) IFLAG = 1
      CALL INITIALIZE_RESIDUALS(IFLAG)

      MMTCNT = 0.
      NMMT = 0
    ENDIF   ! *** END OF OUTPUT PROCESSING
    NTIME = NTIME + 1
  ENDIF   ! *** END OF MMTCNT >= RESSTEP

  IF( ISICM >= 1 )THEN
    NMMT = NMMT + 2
  ELSE
    NMMT = NMMT + 1
  ENDIF

  RETURN

END

SUBROUTINE INITIALIZE_RESIDUALS(IFLAG)

  USE GLOBAL
  Use Variables_MPI

  INTEGER :: IFLAG
  INTEGER :: L, K, NSN, NT, NS, NMD, NWR, LN, LS, LSW, LL, LT, I, J, ITMP, NSC, LE, LW, MD
  REAL    :: UTMP1, VTMP1, UTMP, VTMP, FLTWT, TMPVAL, HPLW, HPLS, HPLSW, HMC

  ! *** INITIALIZE FOR CASE WHEN MEAN MASS TRANPORT INTERVAL IS >= REFERENCE PERIOD
  DO L = 1,LC
    HLPF(L) = 0.
    QSUMELPF(L) = 0.
    UELPF(L) = 0.
    VELPF(L) = 0.
    RAINLPF(L) = 0.
    EVPSLPF(L) = 0.
    EVPGLPF(L) = 0.
    RINFLPF(L) = 0.
    GWLPF(L) = 0.
  ENDDO

  IF( ISTRAN(6) > 0 )THEN
    DO NSC = 1,NSED
      DO K = 1,KB
        DO L = 1,LC
          SEDBLPF(L,K,NSC) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    DO NSN = 1,NSND
      DO K = 1,KB
        DO L = 1,LC
          SNDBLPF(L,K,NSN) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    DO NT = 1,NTOX
      DO K = 1,KB
        DO L = 1,LC
          TOXBLPF(L,K,NT) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  HLPF(1) = HMIN
  HLPF(LC) = HMIN
  DO K = 1,KS
    DO L = 1,LC
      ABEFF(L,K) = 0.
      ABLPF(L,K) = 0.
      WLPF(L,K) = 0.
    ENDDO
  ENDDO

  DO K = 1,KC
    DO L = 1,LC
      AHULPF(L,K) = 0.
      AHVLPF(L,K) = 0.
      SALLPF(L,K) = 0.
      TEMLPF(L,K) = 0.
      SFLLPF(L,K) = 0.
      UHLPF(L,K) = 0.
      VHLPF(L,K) = 0.
      QSUMLPF(L,K) = 0.
    ENDDO
  ENDDO

  IF( ISTRAN(3) > 0 )THEN
    DO MD = 1,NDYE
      DO K = 1,KC
        DO L = 1,LC
          DYELPF(L,K,MD) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(6) > 0 )THEN
    DO NSC = 1,NSED
      DO K = 1,KC
        DO L = 1,LC
          SEDLPF(L,K,NSC) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    DO NSN = 1,NSND
      DO K = 1,KC
        DO L = 1,LC
          SNDLPF(L,K,NSN) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    DO NT = 1,NTOX
      DO K = 1,KC
        DO L = 1,LC
          TOXLPF(L,K,NT) = 0.
        ENDDO
      ENDDO
    ENDDO

    DO NT = 1,NTOX
      DO NS = 1,NSED+NSND
        DO K = 1,KC
          DO L = 1,LC
            TXPFLPF(L,K,NS,NT) = 0.
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! *** SOURCE/SINK TERMS
  DO NS = 1,NQSER
    DO K = 1,KC
      QSRTLPP(K,NS) = 0.
      QSRTLPN(K,NS) = 0.
    ENDDO
  ENDDO
  DO NS = 1,NQCTL
    DO K = 1,KC
      QCTLTLP(K,NS) = 0.
    ENDDO
  ENDDO
  DO NMD = 1,MDCHH
    QCHNULP(NMD) = 0.
    QCHNVLP(NMD) = 0.
  ENDDO
  DO NWR = 1,NQWR
    QWRSERTLP(NWR) = 0.
  ENDDO

  IF( IFLAG > 0 )THEN    ! *** RESSTEP > TIDALP
    DO K = 1,KS
      DO L = 1,LC
        WTLPF(L,K) = 0.
        WIRT(L,K) = 0.
      ENDDO
    ENDDO

    DO K = 1,KC
      DO L = 1,LC
        UIRT(L,K) = 0.
        ULPF(L,K) = 0.
        UTLPF(L,K) = 0.
        VIRT(L,K) = 0.
        VLPF(L,K) = 0.
        VTLPF(L,K) = 0.
      ENDDO
    ENDDO
  ENDIF

  RETURN

  END SUBROUTINE INITIALIZE_RESIDUALS
