! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CALMMT

  ! CHANGE RECORD
  ! ***  SUBROUTINE CALMMTF CALCULATES THE MEAN MASS TRANSPORT FIELD

  use GLOBAL
  use Variables_MPI

  integer      :: L, K, NSN, NT, NS, NMD, NWR, LN, LS, LSW, LL, LT, I, J, ITMP, NSC, LE, LW, MD, IFLAG
  real         :: UTMP1, VTMP1, UTMP, VTMP, FLTWT, TMPVAL, HPLW, HPLS, HPLSW, HMC
  real         :: QXE
  real :: QXEVP
  real :: QXW
  real :: QXWVP
  real :: QYN
  real :: QYNVP
  real :: QYS
  real :: QYSVP
  real, save   :: MMTCNT

  real,save,allocatable,dimension(:,:)   :: VPX
  real,save,allocatable,dimension(:,:)   :: VPY
  real,save,allocatable,dimension(:,:)   :: VPZ

  ! *** ALLOCATE ARRAYS
  if( .not. allocated(VPX) )then
    allocate(VPX(LCM,0:KCM))
    allocate(VPY(LCM,0:KCM))
    allocate(VPZ(LCM,KCM))
    VPX = 0.0
    VPY = 0.0
    VPZ = 0.0

    MMTCNT = 0.
    NMMT = 1

    IFLAG = 0
    if( RESSTEP > TIDALP ) IFLAG = 1
    call INITIALIZE_RESIDUALS(IFLAG)
    
    ! *** Initiate and write the first ouput to WASP linkage file 
#ifdef WASPOUT
    if( ISWASP == 17 ) CALL WASP8HYDRO
#endif
  endif

  ! *** Select time step depending on whether a static or dynamic time step is used
  if( ISDYNSTP == 0 )then
    DELT = DT
  else
    DELT = DTDYN
  endif

  ! ***  INITIALIZE CE-QUAL-ICM INTERFACE
  if( ISICM >= 1 .and. JSWASP == 1 ) CALL CEQICM

  ! *** ACCUMULATE MEAN MASS TRANPORT
  do L = 2,LA
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
  enddo

  if( ISGWIE == 0 )then
    do L = 2,LA
      EVPSLPF(L) = EVPSLPF(L) + DXYP(L)*EVAPT(L)*DELT
    enddo
  else
    do L = 2,LA
      EVPSLPF(L) = EVPSLPF(L) + EVAPSW(L)*DELT         ! *** m3
      EVPGLPF(L) = EVPGLPF(L) + EVAPGW(L)*DELT         ! *** m3
      RINFLPF(L) = RINFLPF(L) - QGW(L)*DELT            ! *** m3
      GWLPF(L)   = GWLPF(L) + AGWELV(L)*DELT           ! *** m3
    enddo
  endif

  if( ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      do K = 1,KB
        do L = 2,LA
          TOXBLPF(L,K,NT) = TOXBLPF(L,K,NT) + TOXB(L,K,NT)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(6) > 0 )then
    do NSC = 1,NSED2
      do K = 1,KB
        do L = 2,LA
          SEDBLPF(L,K,NSC) = SEDBLPF(L,K,NSC) + SEDB(L,K,NSC)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(7) > 0 )then
    do NSN = 1,NSND
      do K = 1,KB
        do L = 2,LA
          SNDBLPF(L,K,NSN) = SNDBLPF(L,K,NSN) + SNDB(L,K,NSN)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISWASP == 99 .or. ISICM >= 1 )then
    ! *** ICM OR RCA LINKAGE
    do K = 1,KS
      do L = 2,LA
        ABEFF(L,K) = ABEFF(L,K) + AB(L,K)*(SAL(L,K+1)-SAL(L,K))*DELT
        ABLPF(L,K) = ABLPF(L,K) + (AB(L,K)*HP(L))*DELT
        WLPF(L,K)  = WLPF(L,K)  + W(L,K)*DELT
      enddo
    enddo

    if( RESSTEP > TIDALP )then
      ! DELME - NOT TESTED
      do K = 1,KS
        do L = 2,LA
          WIRT(L,K)  = WIRT(L,K)  + DT*W(L,K)
          WTLPF(L,K) = WTLPF(L,K) + DT*(FLOAT(NMMT)-0.5)*W(L,K)
        enddo
      enddo
    endif
  else
    ! *** WASP LINKAGE
    do K = 1,KS
      do L = 2,LA
        ABEFF(L,K) = ABEFF(L,K) + AB(L,K)*(SAL(L,K+1)-SAL(L,K))*DELT
        ABLPF(L,K) = ABLPF(L,K) + AB(L,K)*DELT
        WLPF(L,K)  = WLPF(L,K)  + W(L,K)*DELT
      enddo
    enddo

    if( RESSTEP > TIDALP )then
      ! DELME - NOT TESTED
      do K = 1,KS
        do L = 2,LA
          WIRT(L,K)  = WIRT(L,K)  + DT*W(L,K)
          WTLPF(L,K) = WTLPF(L,K) + DT*(FLOAT(NMMT)-0.5)*W(L,K)
        enddo
      enddo
    endif
  endif

  do K = 1,KC
    do L = 2,LA
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
    enddo
  enddo

  if( RESSTEP > TIDALP )then
    ! *** DELME - NOT TESTED
    do K = 1,KC
      do L = 2,LA
        LS = LSC(L)
        LW = LWC(L)
        ULPF(L,K)   = ULPF(L,K) + U(L,K)*DELT
        VLPF(L,K)   = VLPF(L,K) + V(L,K)*DELT

        UIRT(L,K)   = UIRT(L,K) + U(L,K)*DELT
        UTLPF(L,K)  = UTLPF(L,K) + (FLOAT(NMMT)-0.5)*U(L,K)

        VIRT(L,K)   = VIRT(L,K) + V(L,K)*DELT
        VTLPF(L,K)  = VTLPF(L,K) + (FLOAT(NMMT)-0.5)*V(L,K)
      enddo
    enddo
  endif

  if( ISTRAN(3) > 0 )then
    do MD = 1,NDYE
      do K = 1,KC
        do L = 1,LC
          DYELPF(L,K,MD) = DYELPF(L,K,MD) + DYE(L,K,MD)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      do K = 1,KC
        do L = 2,LA
          TOXLPF(L,K,NT) = TOXLPF(L,K,NT) + TOX(L,K,NT)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(6) > 0 )then
    do NSC = 1,NSED2
      do K = 1,KC
        do L = 2,LA
          SEDLPF(L,K,NSC) = SEDLPF(L,K,NSC) + SED(L,K,NSC)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(7) > 0 )then
    do NSN = 1,NSND
      do K = 1,KC
        do L = 2,LA
          SNDLPF(L,K,NSN) = SNDLPF(L,K,NSN) + SND(L,K,NSN)*DELT
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      do NS = 1,NSED2+NSND
        do K = 1,KC
          do L = 1,LC
            TXPFLPF(L,K,NS,NT) = TXPFLPF(L,K,NS,NT) + TOXPFW(L,K,NS,NT)*DELT
          enddo
        enddo
      enddo
    enddo
  endif

  ! *** SOURCE SINK
  do NS = 1,NQSER
    do K = 1,KC
      QSRTLPP(K,NS) = QSRTLPP(K,NS) + MAX(QSERT(K,NS),0.)*DELT
      QSRTLPN(K,NS) = QSRTLPN(K,NS) + MIN(QSERT(K,NS),0.)*DELT
    enddo
  enddo
  do NS = 1,NQCTL
    do K = 1,KC
      QCTLTLP(K,NS) = QCTLTLP(K,NS) + QCTLT(K,NS,1)*DELT
    enddo
  enddo
  do NMD = 1,MDCHH
    QCHNULP(NMD) = QCHNULP(NMD) + QCHANU(NMD)*DELT
    QCHNVLP(NMD) = QCHNVLP(NMD) + QCHANV(NMD)*DELT
  enddo
  do NWR = 1,NQWR
    QWRSERTLP(NWR) = QWRSERTLP(NWR) + QWRSERT(NWR)*DELT
  enddo

  ! ***
  if( RESSTEP > TIDALP )then
    ! DELME - NOT TESTED
    do K = 1,KS
      do L = 2,LA
        LS = LSC(L)
        LW = LWC(L)
        VPX(L,K) = VPX(L,K) + 0.25*(V(L,K+1) + V(L,K))*(WIRT(L,K) + WIRT(LS,K))
        VPY(L,K) = VPY(L,K) + 0.25*(W(L,K)   + W(LW,K))*(UIRT(L,K+1) + UIRT(L,K))
      enddo
    enddo
    do K = 1,KC
      do L = 2,LA
        LS = LSC(L)
        LW = LWC(L)
        VPZ(L,K) = VPZ(L,K) + 0.25*(U(L,K) + U(LS,K))*(VIRT(L,K) + VIRT(LW,K))
      enddo
    enddo
  endif
  ! *****************************************************************************************************************

  ! *** ACCUMULATE THE DELTA T
  MMTCNT = MMTCNT + DELT
  if( ISICM >= 1 ) MMTCNT = MMTCNT + DELT

  ! *****************************************************************************************************************
  ! ***  CHECK FOR END OF FILTER
  if( TIMEDAY >= WASPTIME(NTIME) )then

    ! ***  COMPLETE THE FILTERING
    FLTWT = 1./MMTCNT
    if( ISICM >= 1 ) FLTWT = 2.*FLTWT

    ! *** FINALIZE FOR CASE WHEN MEAN MASS TRANPORT INTERVAL IS < REFERENCE PERIOD
    do L = 2,LA
      HLPF(L)     = FLTWT*HLPF(L)
      QSUMELPF(L) = FLTWT*QSUMELPF(L)
      UELPF(L)    = FLTWT*UELPF(L)
      VELPF(L)    = FLTWT*VELPF(L)
      RAINLPF(L)  = FLTWT*RAINLPF(L)
    enddo

    if( ISGWIE == 0 )then
      do L = 2,LA
        EVPSLPF(L) = FLTWT*EVPSLPF(L)
      enddo
    else
      do L = 2,LA
        EVPSLPF(L) = FLTWT*EVPSLPF(L)
        EVPGLPF(L) = FLTWT*EVPGLPF(L)
        RINFLPF(L) = FLTWT*RINFLPF(L)
        GWLPF(L)   = FLTWT*GWLPF(L)
      enddo
    endif

    if( ISTRAN(5) > 0 )then
      do NT = 1,NTOX
        do K = 1,KB
          do L = 2,LA
            TOXBLPF(L,K,NT) = TOXBLPF(L,K,NT)*FLTWT
          enddo
        enddo
      enddo
    endif

    if( ISTRAN(6) > 0 )then
      do NSC = 1,NSED2
        do K = 1,KB
          do L = 2,LA
            SEDBLPF(L,K,NSC) = SEDBLPF(L,K,NSC)*FLTWT
          enddo
        enddo
      enddo
    endif

    if( ISTRAN(7) > 0 )then
      do NSN = 1,NSND
        do K = 1,KB
          do L = 2,LA
            SNDBLPF(L,K,NSN) = SNDBLPF(L,K,NSN)*FLTWT
          enddo
        enddo
      enddo
    endif

    do K = 1,KS
      do L = 2,LA
        ABEFF(L,K) = FLTWT*ABEFF(L,K)
        ABLPF(L,K) = FLTWT*ABLPF(L,K)
        WLPF(L,K)  = FLTWT*WLPF(L,K)
      enddo
    enddo

    do K = 1,KC
      do L = 2,LA
        AHULPF(L,K)  = FLTWT*AHULPF(L,K)
        AHVLPF(L,K)  = FLTWT*AHVLPF(L,K)
        SALLPF(L,K)  = FLTWT*SALLPF(L,K)
        TEMLPF(L,K)  = FLTWT*TEMLPF(L,K)
        SFLLPF(L,K)  = FLTWT*SFLLPF(L,K)
        UHLPF(L,K)   = FLTWT*UHLPF(L,K)
        VHLPF(L,K)   = FLTWT*VHLPF(L,K)
        QSUMLPF(L,K) = FLTWT*QSUMLPF(L,K)
      enddo
    enddo

    if( ISTRAN(3) > 0 )then
      do MD = 1,NDYE
        do K = 1,KC
          do L = 1,LC
            DYELPF(L,K,MD) = DYELPF(L,K,MD)*FLTWT
          enddo
        enddo
      enddo
    endif

    if( ISTRAN(6) > 0 )then
      do NSC = 1,NSED2
        do K = 1,KC
          do L = 2,LA
            SEDLPF(L,K,NSC) = SEDLPF(L,K,NSC)*FLTWT
          enddo
        enddo
      enddo
    endif

    if( ISTRAN(7) > 0 )then
      do NSN = 1,NSND
        do K = 1,KC
          do L = 2,LA
            SNDLPF(L,K,NSN) = SNDLPF(L,K,NSN)*FLTWT
          enddo
        enddo
      enddo
    endif

    if( ISTRAN(5) > 0 )then
      do NT = 1,NTOX
        do K = 1,KC
          do L = 2,LA
            TOXLPF(L,K,NT) = TOXLPF(L,K,NT)*FLTWT
          enddo
        enddo
      enddo
      do NT = 1,NTOX
        do NS = 1,NSED2+NSND
          do K = 1,KC
            do L = 1,LC
              TXPFLPF(L,K,NS,NT) = TXPFLPF(L,K,NS,NT)*FLTWT
            enddo
          enddo
        enddo
      enddo
    endif

    ! *** SOURCES & SINKS
    do NS = 1,NQSER
      do K = 1,KC
        QSRTLPP(K,NS) = FLTWT*QSRTLPP(K,NS)
        QSRTLPN(K,NS) = FLTWT*QSRTLPN(K,NS)
      enddo
    enddo
    do NS = 1,NQCTL
      do K = 1,KC
        QCTLTLP(K,NS) = FLTWT*QCTLTLP(K,NS)
      enddo
    enddo
    do NMD = 1,MDCHH
      QCHNULP(NMD) = FLTWT*QCHNULP(NMD)
      QCHNVLP(NMD) = FLTWT*QCHNVLP(NMD)
    enddo
    do NWR = 1,NQWR
      QWRSERTLP(NWR) = FLTWT*QWRSERTLP(NWR)
    enddo

    if( RESSTEP > TIDALP )then
      ! DELME - NOT TESTED

      ! *** POTENTIAL TRANSPORT VELOCITY
      do K = 1,KS
        do L = 2,LA
          LS = LSC(L)
          LW = LWC(L)
          VPX(L,K) = VPX(L,K) - 0.25*(VTLPF(L,K+1) + VTLPF(L,K)) *(WLPF(L,K)   + WLPF(LS,K))
          VPY(L,K) = VPY(L,K) - 0.25*(WTLPF(L,K)   + WTLPF(LW,K))*(ULPF(L,K+1) + ULPF(L,K))
        enddo
      enddo
      do K = 1,KC
        do L = 2,LA
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
        enddo
      enddo
      do K = 1,KC
        do L = 2,LA
          LS = LSC(L)
          LN = LNC(L)
          LE = LEC(L)
          UVPT(L,K) = (VPZ(LN,K)-VPZ(L,K))/DYU(L) - DZIC(L,K)*(VPY(L,K)-VPY(L,K-1))      ! DELME DZIC???
          VVPT(L,K) = DZIC(L,K)*(VPX(L,K)-VPX(L,K-1)) - (VPZ(LE,K)-VPZ(L,K))/DXV(L)
        enddo
      enddo
      do K = 1,KS
        do L = 2,LA
          LS = LSC(L)
          LN = LNC(L)
          LE = LEC(L)
          WVPT(L,K) = (VPY(LE,K)-VPY(L,K))/DXP(L)-(VPX(LN,K)-VPX(L,K) )/DYP(L)
        enddo
      enddo
    endif

    ! *** RESIDUAL FLOWS ACROSS OPEN BCs
    QXW = 0.
    QXWVP = 0.
    do K = 1,KC
      do LL = 1,NPBW
        L = LPBW(LL)
        LE = LEC(L)
        QXW   = QXW   + UHLPF(LE,K)*DYU(LE)
        QXWVP = QXWVP + UVPT(LE,K)*DZC(L,K)*DYU(LE)
      enddo
    enddo

    QXE = 0.
    QXEVP = 0.
    do K = 1,KC
      do LL = 1,NPBE
        L = LPBE(LL)
        QXE = QXE+UHLPF(L,K)*DYU(L)
        QXEVP = QXEVP+UVPT(L,K)*DZC(L,K)*DYU(L)
      enddo
    enddo

    QYS = 0.
    QYSVP = 0.
    do K = 1,KC
      do LL = 1,NPBS
        L = LPBS(LL)
        LN = LNC(L)
        QYS = QYS+VHLPF(LN,K)*DXV(LN)
        QYSVP = QYSVP+VVPT(LN,K)*DZC(L,K)*DXV(LN)
      enddo
    enddo

    QYN = 0.
    QYNVP = 0.
    do K = 1,KC
      do LL = 1,NPBN
        L = LPBN(LL)
        QYN   = QYN   + VHLPF(L,K)*DXV(L)
        QYNVP = QYNVP + VVPT(L,K)*DZC(L,K)*DXV(L)
      enddo
    enddo

    ! *****************************************************************************************************************
    ! *** Process output options
    ! ***  OUTPUT RESIDUAL TRANSPORT TO FILE RESTRAN.OUT
    if( ISSSMMT == 1 .or. ISWASP > 0 .or. (ISSSMMT == 2 .and. TIMEDAY >= FLOAT(NTC)*TIDALP) )then

      if( process_id == master_id )then
        if( ISRESTR == 1 )then
          if( JSRESTR == 1 )then
            open(98,FILE = OUTDIR//'RESTRAN.OUT',STATUS = 'UNKNOWN')
            close(98,STATUS = 'DELETE')
            open(98,FILE = OUTDIR//'RESTRAN.OUT',STATUS = 'UNKNOWN')
            JSRESTR = 0
          else
            open(98,FILE = OUTDIR//'RESTRAN.OUT',POSITION = 'APPEND',STATUS = 'UNKNOWN')
          endif
          if( RESSTEP < TIDALP )then
            do LT = 2,LALT
              I = ILLT(LT)
              J = JLLT(LT)
              L = LIJ(I,J)
              write(98,907)HP(L),HLPF(L),QSUMELPF(L)
              write(98,907)(UHLPF(L,K),K = 1,KC)
              write(98,907)(VHLPF(L,K),K = 1,KC)
              write(98,907)(AHULPF(L,K),K = 1,KC)
              write(98,907)(AHVLPF(L,K),K = 1,KC)
              write(98,907)(SALLPF(L,K),K = 1,KC)
              write(98,907)(ABLPF(L,K),K = 1,KS)
              write(98,907)(ABEFF(L,K),K = 1,KS)
            enddo
          else
            do LT = 2,LALT
              I = ILLT(LT)
              J = JLLT(LT)
              L = LIJ(I,J)
              write(98,907)HP(L),HLPF(L),QSUMELPF(L)
              write(98,907)(UHLPF(L,K),K = 1,KC)
              write(98,907)(VHLPF(L,K),K = 1,KC)
              write(98,907)(VPZ(L,K),K = 1,KC)
              write(98,907)(AHULPF(L,K),K = 1,KC)
              write(98,907)(AHVLPF(L,K),K = 1,KC)
              write(98,907)(SALLPF(L,K),K = 1,KC)
              write(98,907)(VPX(L,K),K = 1,KS)
              write(98,907)(VPY(L,K),K = 1,KS)
              write(98,907)(ABLPF(L,K),K = 1,KS)
            enddo
          endif
          close(98)
        endif
907     FORMAT(12E12.4)
#ifdef WASPOUT
        ! ***  OUTPUT TO WASP COMPATIABLE FILES
        !IF( ISWASP == 4  ) CALL WASP4             delme - untested
        !IF( ISWASP == 5  ) CALL WASP5
        !IF( ISWASP == 6  ) CALL WASP6
        if( ISWASP == 7  ) CALL WASP7
        if( ISWASP == 17 ) CALL WASP8HYDRO
        !IF( ISRCA  >= 1  ) CALL RCAHQ
        !IF( ISICM  >= 1  ) CALL CEQICM
#endif
        ! ***  WRITE GRAPHICS FILES FOR RESIDUAL VARIABLES
        ! ***  RESIDUAL SALINITY CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
        if( ISRSPH(1) == 1 .and. ISTRAN(1) >= 1 )then
          call RSALPLTH(1,SALLPF)
        endif
        if( ISRSPH(2) == 1 .and. ISTRAN(2) >= 1 )then
          call RSALPLTH(2,TEMLPF)
        endif
        if( ISRSPH(3) == 1 .and. ISTRAN(3) >= 1 )then
          do MD = 1,NDYE
            call RSALPLTH(3,DYELPF(1,1,MD))
          enddo
        endif
        if( ISRSPH(4) == 1 .and. ISTRAN(4) >= 1 )then
          call RSALPLTH(4,SFLLPF)
        endif

        if( ISTRAN(6) > 0 )then
          do K = 2,KB
            do L = 2,LA
              SEDBTLPF(L,K) = 0.
            enddo
          enddo
          do K = 1,KC
            do L = 2,LA
              SEDTLPF(L,K) = 0.
            enddo
          enddo
        endif

        if( ISTRAN(7) > 0 )then
          do K = 2,KB
            do L = 2,LA
              SNDBTLPF(L,K) = 0.
            enddo
          enddo
          do K = 1,KC
            do L = 2,LA
              SNDTLPF(L,K) = 0.
            enddo
          enddo
        endif

        if( ISRSPH(5) == 1 .and. ISTRAN(5) >= 1 )then
          do K = 1,KC
            do L = 2,LA
              TVAR1S(L,K) = TOXLPF(L,K,1)
            enddo
          enddo
          do NT = 1,NTOX
            call RSALPLTH(5,TVAR1S)   ! *** DELME - WHY IS TOXLPF NOT PASSED WITHOUT COPYING THE ARRAY?
          enddo

        endif

        if( ISTRAN(6) > 0 )then
          do NS = 1,NSED2
            do K = 1,KB
              do L = 2,LA
                SEDBTLPF(L,K) = SEDBTLPF(L,K) + SEDBLPF(L,K,NS)
              enddo
            enddo
          enddo
          do NS = 1,NSED2
            do K = 1,KC
              do L = 2,LA
                SEDTLPF(L,K) = SEDTLPF(L,K) + SEDLPF(L,K,NS)
              enddo
            enddo
          enddo

          if( ISRSPH(6) == 1 .and. ISTRAN(6) >= 1 )then
            do NSC = 1,NSED2
              call RSALPLTH(6,SEDTLPF)
            enddo
          endif
        endif

        if( ISTRAN(7) > 0 )then
          do NS = 1,NSND
            do K = 1,KB
              do L = 2,LA
                SNDBTLPF(L,K) = SNDBTLPF(L,K)+SNDBLPF(L,K,NS)
              enddo
            enddo
          enddo
          do NS = 1,NSND
            do K = 1,KC
              do L = 2,LA
                SNDTLPF(L,K) = SNDTLPF(L,K) + SNDLPF(L,K,NS)
              enddo
            enddo
          enddo

          if( ISRSPH(7) == 1 .and. ISTRAN(7) >= 1 )then
            do NSN = 1,NSND
              call RSALPLTH(7,SNDTLPF)    ! ***  DELME THIS WRITE THE SAME ARRAY NSND TIMES?
            enddo
          endif
        endif

        ! ***  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:    RVELPLTH
        !if( ISRVPH >= 1 ) CALL RVELPLTH

        ! ***  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:  RSURFPLT
        !if( ISRPPH == 1 ) CALL RSURFPLT

        ! ***  RESIDUAL SALINITY AND VERTICAL MASS DIFFUSIVITY CONTOURING IN 3 VERTICAL PLANES:  RSALPLTV
        !do ITMP = 1,7
        !  if( ISRSPV(ITMP) >= 1 ) CALL RSALPLTV(ITMP)
        !enddo

        ! ***  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND
        ! ***  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:  RVELPLTV
        !if( ISRVPV >= 1 ) CALL RVELPLTV

        ! ***  RESIDUAL 3D SCALAR AND VECTOR OUTPUT FILES
        !if( ISR3DO >= 1 ) CALL ROUT3D
      endif   ! *** End of master process

      IFLAG = 0
      if( RESSTEP > TIDALP ) IFLAG = 1
      call INITIALIZE_RESIDUALS(IFLAG)

      MMTCNT = 0.
      NMMT = 0
    endif   ! *** END OF OUTPUT PROCESSING
    NTIME = NTIME + 1
  endif   ! *** END OF MMTCNT >= RESSTEP

  if( ISICM >= 1 )then
    NMMT = NMMT + 2
  else
    NMMT = NMMT + 1
  endif

  return

END

SUBROUTINE INITIALIZE_RESIDUALS(IFLAG)

  use GLOBAL
  use Variables_MPI

  integer :: IFLAG
  integer :: L, K, NSN, NT, NS, NMD, NWR, LN, LS, LSW, LL, LT, I, J, ITMP, NSC, LE, LW, MD
  real    :: UTMP1, VTMP1, UTMP, VTMP, FLTWT, TMPVAL, HPLW, HPLS, HPLSW, HMC

  ! *** INITIALIZE FOR CASE WHEN MEAN MASS TRANPORT INTERVAL IS >= REFERENCE PERIOD
  do L = 1,LC
    HLPF(L) = 0.
    QSUMELPF(L) = 0.
    UELPF(L) = 0.
    VELPF(L) = 0.
    RAINLPF(L) = 0.
    EVPSLPF(L) = 0.
    EVPGLPF(L) = 0.
    RINFLPF(L) = 0.
    GWLPF(L) = 0.
  enddo

  if( ISTRAN(6) > 0 )then
    do NSC = 1,NSED2
      do K = 1,KB
        do L = 1,LC
          SEDBLPF(L,K,NSC) = 0.
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(7) > 0 )then
    do NSN = 1,NSND
      do K = 1,KB
        do L = 1,LC
          SNDBLPF(L,K,NSN) = 0.
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(7) > 0 )then
    do NT = 1,NTOX
      do K = 1,KB
        do L = 1,LC
          TOXBLPF(L,K,NT) = 0.
        enddo
      enddo
    enddo
  endif

  HLPF(1) = HMIN
  HLPF(LC) = HMIN
  do K = 1,KS
    do L = 1,LC
      ABEFF(L,K) = 0.
      ABLPF(L,K) = 0.
      WLPF(L,K) = 0.
    enddo
  enddo

  do K = 1,KC
    do L = 1,LC
      AHULPF(L,K) = 0.
      AHVLPF(L,K) = 0.
      SALLPF(L,K) = 0.
      TEMLPF(L,K) = 0.
      SFLLPF(L,K) = 0.
      UHLPF(L,K) = 0.
      VHLPF(L,K) = 0.
      QSUMLPF(L,K) = 0.
    enddo
  enddo

  if( ISTRAN(3) > 0 )then
    do MD = 1,NDYE
      do K = 1,KC
        do L = 1,LC
          DYELPF(L,K,MD) = 0.
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(6) > 0 )then
    do NSC = 1,NSED2
      do K = 1,KC
        do L = 1,LC
          SEDLPF(L,K,NSC) = 0.
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(7) > 0 )then
    do NSN = 1,NSND
      do K = 1,KC
        do L = 1,LC
          SNDLPF(L,K,NSN) = 0.
        enddo
      enddo
    enddo
  endif

  if( ISTRAN(7) > 0 )then
    do NT = 1,NTOX
      do K = 1,KC
        do L = 1,LC
          TOXLPF(L,K,NT) = 0.
        enddo
      enddo
    enddo

    do NT = 1,NTOX
      do NS = 1,NSED2+NSND
        do K = 1,KC
          do L = 1,LC
            TXPFLPF(L,K,NS,NT) = 0.
          enddo
        enddo
      enddo
    enddo
  endif

  ! *** SOURCE/SINK TERMS
  do NS = 1,NQSER
    do K = 1,KC
      QSRTLPP(K,NS) = 0.
      QSRTLPN(K,NS) = 0.
    enddo
  enddo
  do NS = 1,NQCTL
    do K = 1,KC
      QCTLTLP(K,NS) = 0.
    enddo
  enddo
  do NMD = 1,MDCHH
    QCHNULP(NMD) = 0.
    QCHNVLP(NMD) = 0.
  enddo
  do NWR = 1,NQWR
    QWRSERTLP(NWR) = 0.
  enddo

  if( IFLAG > 0 )then    ! *** RESSTEP > TIDALP
    do K = 1,KS
      do L = 1,LC
        WTLPF(L,K) = 0.
        WIRT(L,K) = 0.
      enddo
    enddo

    do K = 1,KC
      do L = 1,LC
        UIRT(L,K) = 0.
        ULPF(L,K) = 0.
        UTLPF(L,K) = 0.
        VIRT(L,K) = 0.
        VLPF(L,K) = 0.
        VTLPF(L,K) = 0.
      enddo
    enddo
  endif

  return

  END SUBROUTINE INITIALIZE_RESIDUALS
