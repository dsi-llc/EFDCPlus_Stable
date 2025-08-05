! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE WASP6

  ! CHANGE RECORD
  ! == == == ==  == 
  ! REVISIONS:
  ! == == == ==  == 
  !  M. MORTON 06/06/94: THIS VERSION WRITES DISPERSION TO THE WASPDH.OUT
  !  M. MORTON 06/07/94: WRITES HYDRODYNAMIC INFORMATION AND DISPERSION TO
  !        DATA GROUP B use WASPB.MRM  (DO NOT use WASPB.OUT)
  !        DATA GROUP C use WASPC.OUT
  !        DATA GROUP D use WASPD.MRM  (DO NOT use WASPD.OUT)
  ! == == == == ==  = 
  ! ***  SUBROUTINE WASP5 WRITES OUTPUT FILES PROVIDING ADVECTIVE AND
  ! ***  DIFFUSIVE TRANSPORT FIELDS FOR THE WASP6 WATER QUALITY MODEL
  !
  use GLOBAL

  integer :: L,K,I,J,LT,LCLTM2,LTYPE,KWASP,LBELOW,LWSPTMP,NTEX,NORSH,NORSV,KMUL,LWASPW,LW,NCHNH,NCHNV,ISTMP
  integer :: KMUL1,KMUL2,KMUL3,NBRK,IBPTMP,NINQ,NOQSH,NOQSV,NOQS,LL,LN,NORS,LSLT,LWASPS,LS,NBRKQ,LDTM,LUTM
  integer :: NTEXX,NJUN,NCHN,NODYN,LCHN,LCELL,IZERO,LCHNUM,IMTMP,IPTMP,JMTMP,JPTMP,KPTMP,LCELTMP,LE1
  integer(IK4),save,allocatable,dimension(:) :: LDTMP
  integer(IK4),save,allocatable,dimension(:) :: LUTMP

  real    :: SVPT,SCALR,WSS1,WSS2,WSS3,VOLUME,DXYSUM,UNITY,ADDL,TSTOP,TSTART,TSMALL,TZERO,TENDHYD
  real    :: T1,T2,T3,T4,T5,T6,D1,D2,D3,D4,D5,D6,ADDLW,ADDLS,DTWASP,RMNDUM,RLENTH,WIDTH,VELTMP,DUMVOL,DEPTMP
  real    :: VOLTMP,RZERO,FLOWX,FLOWY,FLOWZ,QQSUM,VOLUM,DEPTH,VELX,VELY,VELZ,VELMAG
  real(RK4),save,allocatable,dimension(:) :: QTMP

  character*50 TITLEB,TITLEC

  if( .not. allocated(LDTMP) )then
    allocate(LDTMP((KCM+1)*LCM))
    allocate(LUTMP((KCM+1)*LCM))
    allocate(QTMP((KCM+1)*LCM))
    LDTMP = 0
    LUTMP = 0
    QTMP = 0.0
  endif
  TITLEB = 'DATA GROUP B: EXCHANGE COEFFICIENTS'
  TITLEC = 'DATA GROUP C: VOLUMES'
  !
  ! ***  WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ***  THE VALUE OF X IN THE F10.X FORMATS MAY NEED TO BE CHANGED
  ! ***  FROM PROBLEM TO PROBLEM.  A PRELIMINARY RUN USING E10.3
  ! ***  CAN BE USED TO SPEED THE ADJUSTMENT
  ! ***  READ CONTROL DATA FOR WRITING TO WASP COMPATIBLE FILES
  !
  SVPT = 1.
  if( NTSMMT < NTSPTC)SVPT = 0.
  if( JSWASP == 1 )then
    write(*,'(A)')'READING EFDC.WSP'
    open(1,FILE = 'EFDC.WSP',STATUS = 'UNKNOWN')
    read(1,1)
    read(1,1)
    read(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
    read(1,1)
    read(1,1)
    read(1,*) NRFLD,SCALR,CONVR,ISNKH
    read(1,1)
    read(1,1)
    read(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD
    read(1,1)
    read(1,1)
    read(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3
    close(1)
  endif
      1 FORMAT (80X)
  !
  ! ***  WRITE HORIZONTAL POSITION AND LAYER FILE WASPP.OUT
  ! ***  WRITE INITIAL VOLUME FILE WASPC.OUT
  ! ***  FILE WASPC.OUT IS CONSISTENT WITH DATA GROUP C SPECIFICATIONS
  ! ***  ON PAGE 11 OF THE WASP5.1 MANUAL PART B, SEPT 1993
  ! ***  FILE WASPP.OUT DEFINES THE LAYER (1 IS SURFACE WATER LAYER, WITH
  ! ***  LAYER NUMBERING INCREASING WITH DEPTH IN WATER COLUMN) AND
  ! ***  HORIZONTAL POSITIONS IN LON,LAT OR UTME, UTMN OF THE WATER
  ! ***  QUALITY (LONG TERM TRANSPORT) CELLS OR SEGEMENTS
  !
  if( JSWASP == 1 )then
    open(90,FILE = OUTDIR//'wasp\WASPP.OUT',STATUS = 'UNKNOWN')
    open(93,FILE = OUTDIR//'wasp\WASPC.OUT',STATUS = 'UNKNOWN')
    close(90,STATUS = 'DELETE')
    close(93,STATUS = 'DELETE')
    open(90,FILE = OUTDIR//'wasp\WASPP.OUT',STATUS = 'UNKNOWN')
    open(93,FILE = OUTDIR//'wasp\WASPC.OUT',STATUS = 'UNKNOWN')
  !
  !       IVOPT = 2
  !       IBEDV = 0
  !
    write(93,1031)IVOPT,IBEDV,TDINTS,TITLEC
  !
  !       SCALV = 1.
  !       CONVV = 1.
  !
    write(93,1032)SCALV,CONVV
  !
  !       VMULT = 0.
  !       VEXP = 0.
  !       DMULT = 0.
  !       DEXP = 0.
  !
    LCLTM2 = LCLT-2
    LWASP = 0
    if( KC > 1 )then
      LTYPE = 1
      KWASP = 1
      do LT = 2,LALT
        LWASP = LWASP+1
        LBELOW = LWASP+LCLTM2
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        DMULT = HLPF(L)*DZC(L,KC)
        VOLUME = DXYP(L)*HLPF(L)*DZC(L,KC)
        write(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
        write(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, &
            DMULT,DEXP
      enddo
      LTYPE = 2
      do K = KS,2,-1
        KWASP = KC-K+1
        do LT = 2,LALT
          LWASP = LWASP+1
          LBELOW = LWASP+LCLTM2
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          DMULT = HLPF(L)*DZC(L,K)
          VOLUME = DXYP(L)*HLPF(L)*DZC(L,K)
          write(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
          write(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, &
              DMULT,DEXP
        enddo
      enddo
    endif
    LTYPE = 2
    if( KC == 1 ) LTYPE = 1
    KWASP = KC
    do LT = 2,LALT
      LWASP = LWASP+1
  !
  !        LBELOW = 0
  !
      LBELOW = LWASP+LCLTM2
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      DMULT = HLPF(L)*DZC(L,KSZ(L))
      VOLUME = DXYP(L)*HLPF(L)*DZC(L,KSZ(L))
      write(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
      write(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, &
          DMULT,DEXP
    enddo
    LTYPE = 3
    KWASP = KC+1
    DXYSUM = 0.
    LWSPTMP = LWASP+1
    do LT = 2,LALT
      LWSPTMP = LWSPTMP+1
    enddo
  !
  ! THE FOLLOWING THE LOWER BENTHIC LAYER.  ALL UPPER BENTHIC LAYER SEGMEN
  ! HAVE THIS LAYER IMMEDIATELY BELOW THEM:
  !
    do LT = 2,LALT
      LWASP = LWASP+1
      LBELOW = LWSPTMP
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      DXYSUM = DXYSUM+DXYP(L)
      VOLUME = DXYP(L)*DEPSED
      write(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
      write(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, &
          DEPSED,DEXP
    enddo
  !
  ! NEXT DO THE LOWER BENTHIC LAYER:
  !
    LTYPE = 4
    KWASP = KC+2
    LWASP = LWASP+1
    LBELOW = 0
    DMULT = DEPSED
    VOLUME = DXYSUM*DEPSED
    write(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
    write(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, &
        DMULT,DEXP
    close(90)
    close(93)
  endif
   1001 FORMAT(4I5,2F10.4)
   1031 FORMAT(2I5,F10.4,10X,A50)
   1032 FORMAT(2F10.4)
  !
  ! FORMAT 1033 AS COMMENTED OUT IS TROUBLESOME ... BETTER CHANGE SHOWN
  !
   1033 FORMAT(3I10,F10.1,4F10.3)
  !
  ! ***  WRITE DIFFUSIVE AND DISPERSIVE TRANSPORT FILE WASPB.OUT
  ! ***  FILE WASPB.OUT IS CONSISTENT WITH DATA GROUP B SPECIFICATIONS
  ! ***  ON PAGE 8 OF THE WASP5.1 MANUAL PART B, SEPT 1993
  !
  if( JSWASP == 1 )then
    open(91,FILE = OUTDIR//'wasp\WASPB.OUT',STATUS = 'UNKNOWN')
    close(91,STATUS = 'DELETE')
    open(91,FILE = OUTDIR//'wasp\WASPB.OUT',STATUS = 'UNKNOWN')
  !
  !       NRFLD = 1
  !
    write(91,1011)NRFLD,TITLEB
    NTEX = NTS/NTSMMT
  !
  !       SCALR = 1.
  !       CONVR = 1.
  !
    write(91,1012)NTEX,SCALR,CONVR
    close(91)
    open(91,FILE = OUTDIR//'wasp\WASPB.OUT',POSITION = 'APPEND' ,STATUS = 'UNKNOWN')
    LCLTM2 = LCLT-2
    NORSH = 0
    NORSV = 0
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      NORSH = NORSH+INT(SUBO(L))+INT(SVBO(L))
      NORSV = NORSV+INT(SPB(L))
    enddo
    NORS = ISNKH*KC*NORSH+KS*NORSV
    write(91,1013)NORS
    if( ISNKH == 1 )then
      UNITY = 1.
      do K = KC,1,-1
        KMUL = KC-K
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SUB(L) == 1. )then
            LWASP = LT-1+KMUL*LCLTM2
            LWASPW = LWASP-1
            LW = LWC(L)
            ADDLW = DYU(L)*AHULPF(L,K)*DZC(L,K)*0.5*(HLPF(L) &
                +HLPF(LW))*DXIU(L)
            write(91,1014) ADDLW,UNITY,LWASPW,LWASP
          endif
        enddo
      enddo
    endif
    if( ISNKH == 1 )then
      UNITY = 1.
      do K = KC,1,-1
        KMUL = KC-K
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SVB(L) == 1. )then
            LWASP = LT-1+KMUL*LCLTM2
            LSLT = LSCLT(LT)
            LWASPS = LSLT-1+KMUL*LCLTM2
            LS = LSC(L)
            ADDLS = DXV(L)*AHVLPF(L,K)*DZC(L,K)*0.5*(HLPF(L) &
                +HLPF(LS))*DYIV(L)
            write(91,1014) ADDLS,UNITY,LWASPS,LWASP
          endif
        enddo
      enddo
    endif
    if( KC > 1 )then
      UNITY = 1.
      do K = KS,1,-1
        KMUL1 = KS-K
        KMUL2 = KMUL1+1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SPB(L) == 1. )then
            LWASP = LT-1+KMUL1*LCLTM2
            LBELOW = LT-1+KMUL2*LCLTM2
            ADDL = DXYP(L)*ABLPF(L,K)*DZIG(L,K)
            write(91,1014) ADDL,UNITY,LWASP,LBELOW
          endif
        enddo
      enddo
    endif
    NBRK = 6
    write(91,1015)NBRK
    
    TSTOP = TIMESEC
    TSTART = TSTOP-DT*FLOAT(NTSMMT)
  
    TSTOP = TSTOP/86400.
    TSTART = TSTART/86400.
    TSMALL = 1.E-5
    D1 = 0.
    T1 = 0.-2*TSMALL
    D2 = 0.
    T2 = TSTART-TSMALL
    D3 = 1.
    T3 = TSTART+TSMALL
    D4 = 1.
    T4 = TSTOP-TSMALL
    D5 = 0.
    T5 = TSTOP+TSMALL
    D6 = 0.
    T6 = 2*TSMALL+(DT*FLOAT(NTS)+TBEGIN*TCON)/86400.
    write(91,1016)D1,T1,D2,T2,D3,T3,D4,T4
    write(91,1016)D5,T5,D6,T6
    close(91)
  !
  ! ***  ADD PORE WATER EXCHANGE FIELD ON LAST CALL
  !
    open(91,FILE = OUTDIR//'wasp\WASPB.OUT',POSITION = 'APPEND' ,STATUS = 'UNKNOWN')
    NTEX = 1
    SCALR = 1.
    CONVR = 1.
    write(91,1012)NTEX,SCALR,CONVR
    NORSV = 0
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      NORSV = NORSV+INT(SPB(L))
    enddo
    write(91,1013)NORSV
    if( KC >= 1 )then
      KMUL2 = KC+1
      UNITY = 1.
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        if( SPB(L) == 1. )then
          LWASP = LT-1+KC*LCLTM2
          LBELOW = LT-1+KMUL2*LCLTM2
          ADDL = 2.*DXYP(L)*SEDIFF/DEPSED
          write(91,1014) ADDL,UNITY,LWASP,LBELOW
        endif
      enddo
    endif
    NBRK = 6
    write(91,1015)NBRK
    if( ISDYNSTP == 0 )then
      TSTOP = DT*FLOAT(N)+TCON*TBEGIN
      TSTART = TSTOP-DT*FLOAT(NTSMMT)
    else
      TSTOP = TIMESEC
      TSTART = TSTOP-DT*FLOAT(NTSMMT)
    endif
    TSTOP = TSTOP/86400.
    TSTART = TSTART/86400.
    TSMALL = 1.E-5
    D1 = 0.
    T1 = 0.-2*TSMALL
    D2 = 0.
    T2 = TSTART-TSMALL
    D3 = 1.
    T3 = TSTART+TSMALL
    D4 = 1.
    T4 = TSTOP-TSMALL
    D5 = 0.
    T5 = TSTOP+TSMALL
    D6 = 0.
    T6 = 2*TSMALL+(DT*FLOAT(NTS)+TBEGIN*TCON)/86400.
    write(91,1016)D1,T1,D2,T2,D3,T3,D4,T4
    write(91,1016)D5,T5,D6,T6
    IBPTMP = 0
    write(91,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
        IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
        IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
        IBPTMP,IBPTMP,IBPTMP,IBPTMP
    close(91)
  endif
   1011 FORMAT(I5,10X,A50)
   1012 FORMAT(I5,2F10.4)
   1013 FORMAT(I5)
   1014 FORMAT(2E10.3,2I5)
   1015 FORMAT(I5)
   1016 FORMAT(4(E10.3,F10.5))
   1017 FORMAT(16I5)
  !
  ! ***  WRITE ADVECTIVE TRANSPORT FILE WASPD.OUT
  ! ***  FILE WASPD.OUT IS CONSISTENT WITH DATA GROUP D.1 SPECIFICATIONS
  ! ***  ON PAGE 13 OF THE WASP5.1 MANUAL PART B, SEPT 1993
  ! ***  THIS FILE IS WRITTEN ONLY IF ISWASPD = 1
  !!!!!!!!!!!CHANGES ON NEXT 2 LINES
  !
  if( ISWASPD == 1 )then
    if( JSWASP == 1 )then
      open(92,FILE = OUTDIR//'wasp\WASPD.OUT',STATUS = 'UNKNOWN')
      close(92,STATUS = 'DELETE')
      open(92,FILE = OUTDIR//'wasp\WASPD.OUT',STATUS = 'UNKNOWN')
  !
  !       IQOPT = 1
  !       NFIELD = 1
  !
      write(92,1021)IQOPT,NFIELD,HYDFIL
      NINQ = NTS/NTSMMT
  !
  !       SCALQ = 1
  !       CONVQ = 1
  !
      write(92,1022)NINQ,SCALQ,CONVQ
      close(92)
    endif
    open(92,FILE = OUTDIR//'wasp\WASPD.OUT',POSITION = 'APPEND' ,STATUS = 'UNKNOWN')
    LCLTM2 = LCLT-2
    NOQSH = 0
    NOQSV = 0
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
  !
  !!!!!!!!!!!!!!!CHANGES ON NEXT 3 LINES
  !
      NOQSH = NOQSH+INT(SUBO(L))+INT(SVBO(L))
      if( IJCTLT(I+1,J) == 8 ) NOQSH = NOQSH+1
      if( IJCTLT(I,J+1) == 8 ) NOQSH = NOQSH+1
      NOQSV = NOQSV+INT(SWB(L))
    enddo
    NOQS = KC*NOQSH+KS*NOQSV
    write(92,1023)NOQS
    LL = 0
    do K = KC,1,-1
      KMUL = KC-K
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
  !
  !!!!!!!!!!!!!!CHANGES ON NEXT 15 LINES
  !
        if( SUBO(L) == 1. )then
          LL = LL+1
          LDTMP(LL) = LT-1+KMUL*LCLTM2
          LUTMP(LL) = LDTMP(LL)-1
          if( IJCTLT(I-1,J) == 8 ) LUTMP(LL) = 0
          QTMP(LL) = DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(L,K)
        endif
        if( IJCTLT(I+1,J) == 8 )then
          if( SUBO(LEC(L)) == 1. )then
            LL = LL+1
            LDTMP(LL) = 0
            LUTMP(LL) = LT-1+KMUL*LCLTM2
            QTMP(LL) = DYU(LEC(L))*(UHLPF(LEC(L),K)+SVPT*UVPT(LEC(L),K))*DZC(L,K)
          endif
        endif
      enddo
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
  !
  !!!!!!!!!!!!!!CHANGES ON NEXT 16 LINES
  !
        if( SVBO(L) == 1. )then
          LL = LL+1
          LSLT = LSCLT(LT)
          LDTMP(LL) = LT-1+KMUL*LCLTM2
          LUTMP(LL) = LSLT-1+KMUL*LCLTM2
          if( IJCTLT(I,J-1) == 8 ) LUTMP(LL) = 0
          QTMP(LL) = DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(L,K)
        endif
        if( IJCTLT(I,J+1) == 8 )then
          LN = LNC(L)
          if( SVBO(LN) == 1 )then
            LL = LL+1
            LDTMP(LL) = 0
            LUTMP(LL) = LT-1+KMUL*LCLTM2
            QTMP(LL) = DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(L,K)
          endif
        endif
      enddo
    enddo
    if( KC > 1 )then
      do K = KS,1,-1
        KMUL1 = KS-K
        KMUL2 = KMUL1+1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SWB(L) == 1. )then
            LL = LL+1
            LUTMP(LL) = LT-1+KMUL1*LCLTM2
            LDTMP(LL) = LT-1+KMUL2*LCLTM2
            QTMP(LL) = -DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
          endif
        enddo
      enddo
    endif
    do L = 1,LL,4
      LE1 = LEC(LEC(L))
      write(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L), &
          QTMP(LEC(L)),LUTMP(LEC(L)),LDTMP(LEC(L)), &
          QTMP(LE1),LUTMP(LE1),LDTMP(LE1), &
          QTMP(LEC(LE1)),LUTMP(LEC(LE1)),LDTMP(LEC(LE1))
    enddo
    NBRKQ = 6
    write(92,1025)NBRKQ
    write(92,1026)D1,T1,D2,T2,D3,T3,D4,T4
    write(92,1026)D5,T5,D6,T6
    close(92)
  !
  !!!!!!!!!!!CHANGES ON NEXT 2 LINES
  !
  endif
   1021 FORMAT(2I5,A12)
   1022 FORMAT(I5,2F10.4)
   1023 FORMAT(I5)
   1024 FORMAT(4(E10.3,2I5))
   1025 FORMAT(I5)
   1026 FORMAT(4(2F10.5))
  if( JSWASP  ==  1 )then
    open(92,FILE = OUTDIR//'wasp\WASPD.MRM',STATUS = 'UNKNOWN')
    write(92,2020) IQOPT,NFIELD,HYDFIL
    NINQ = 0
    SCALQ = 1.0
    CONVQ = 1.0/86400.0
  !
  ! DATA BLOCK D.1 (ADVECTIVE FLOWS) IS NOT NEEDED SINCE HYD FILE IS USED:
  ! DATA BLOCK D.2 (PORE WATER FLOWS) NOT NEEDED:
  !
    write(92,2022) NINQ,SCALQ,CONVQ
  !
  ! DATA BLOCK D.3 (SEDIMENT #1 TRANSPORT FIELD):
  !
    NINQ = 1
    write(92,2023) NINQ,SCALQ,CONVQ
    if( KC > 1 )then
      do K = KS,0,-1
        KMUL1 = KS-K
        KMUL2 = KMUL1+1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SWB(L) == 1. )then
            LL = LL+1
            LUTMP(LL) = LT-1+KMUL1*LCLTM2
            LDTMP(LL) = LT-1+KMUL2*LCLTM2
  !
  ! QTMP ARRAY WILL HOLD THE PLAN VIEW AREA OF EACH CELL:
  !
            QTMP(LL)= DXYP(L)
          endif
        enddo
      enddo
    endif
    write(92,2030) LL
    do L = 1,LL,4
      LE1 = LEC(LEC(L))
      write(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L), &
          QTMP(LEC(L)),LUTMP(LEC(L)),LDTMP(LEC(L)), &
          QTMP(LE1),LUTMP(LE1),LDTMP(LE1), &
          QTMP(LEC(LE1)),LUTMP(LEC(LE1)),LDTMP(LEC(LE1))
    enddo
    NBRKQ = 2
    T1 = 1.0
    T2 = 366.0
    write(92,2030) NBRKQ
    write(92,2031) WSS1,T1,WSS1,T2
  !
  ! DATA BLOCK D.4 (SEDIMENT #2 TRANSPORT FIELD):
  !
    NINQ = 1
    write(92,2024) NINQ,SCALQ,CONVQ
    write(92,2030) LL
    do L = 1,LL,4
      LE1 = LEC(LEC(L))
      write(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L), &
          QTMP(LEC(L)),LUTMP(LEC(L)),LDTMP(LEC(L)), &
          QTMP(LE1),LUTMP(LE1),LDTMP(LE1), &
          QTMP(LEC(LE1)),LUTMP(LEC(LE1)),LDTMP(LEC(LE1))
    enddo
    NBRKQ = 2
    T1 = 1.0
    T2 = 366.0
    write(92,2030) NBRKQ
    write(92,2031) WSS2,T1,WSS2,T2
  !
  ! DATA BLOCK D.5 (SEDIMENT #3 TRANSPORT FIELD):
  !
    NINQ = 1
    write(92,2025) NINQ,SCALQ,CONVQ
    write(92,2030) LL
    do L = 1,LL,4
      LE1 = LEC(LEC(L))
      write(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L), &
          QTMP(LEC(L)),LUTMP(LEC(L)),LDTMP(LEC(L)), &
          QTMP(LE1),LUTMP(LE1),LDTMP(LE1), &
          QTMP(LEC(LE1)),LUTMP(LEC(LE1)),LDTMP(LEC(LE1))
    enddo
    NBRKQ = 2
    T1 = 1.0
    T2 = 366.0
    write(92,2030) NBRKQ
    write(92,2031) WSS3,T1,WSS3,T2
  !
  ! ADD SYSTEM BYPASS ARRAY TO BOTTOM OF DATA GROUP D:
  !
    write(92,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
        IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
        IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
        IBPTMP,IBPTMP,IBPTMP,IBPTMP
    close(92)
  endif
   2020 FORMAT(2I5,A12,'  DATA GROUP D: FLOWS')
   2021 FORMAT(I5,2E10.3,'  DATA BLOCK D.1 ADVECTIVE FLOWS')
   2022 FORMAT(I5,2E10.3,'  DATA BLOCK D.2 PORE WATER FLOWS')
   2023 FORMAT(I5,2E10.3,'  DATA BLOCK D.3 SEDIMENT #1 TRANSPORT FIELD')
   2024 FORMAT(I5,2E10.3,'  DATA BLOCK D.4 SEDIMENT #2 TRANSPORT FIELD')
   2025 FORMAT(I5,2E10.3,'  DATA BLOCK D.5 SEDIMENT #3 TRANSPORT FIELD')
   2030 FORMAT(I5)
   2031 FORMAT(2(E10.3,F10.5))
  !
  ! ***  WRITE TO EXTERNAL HYDRO FILE WASPDH.OUT AND DIAGNOSTIC VERSION
  ! ***  OF SAME FILE WASPDHD.OUT
  !
  if( JSWASP == 1 )then
    open(90,FILE = OUTDIR//'wasp\WASPDHD.OUT',STATUS = 'UNKNOWN')
    if( IQOPT == 3 ) open(94,FILE = OUTDIR//'wasp\WASPDH.OUT'  ,STATUS = 'UNKNOWN')
    if( IQOPT == 4 ) open(95,FILE = OUTDIR//'wasp\WASPDHU.OUT' ,STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')
    open(96,FILE = OUTDIR//'wasp\WASPB.MRM',STATUS = 'UNKNOWN')
    close(90,STATUS = 'DELETE')
    if( IQOPT == 3 ) close(94,STATUS = 'DELETE')
    if( IQOPT == 4 ) close(95,STATUS = 'DELETE')
    close(96,STATUS = 'DELETE')
    open(90,FILE = OUTDIR//'wasp\WASPDHD.OUT',STATUS = 'UNKNOWN')
    if( IQOPT == 3 ) open(94,FILE = OUTDIR//'wasp\WASPDH.OUT'  ,STATUS = 'UNKNOWN')
    if( IQOPT == 4 ) open(95,FILE = OUTDIR//'wasp\WASPDHU.OUT' ,STATUS = 'UNKNOWN',  FORM = 'UNFORMATTED')
    open(96,FILE = OUTDIR//'wasp\WASPB.MRM',STATUS = 'UNKNOWN')
    write(96,1011) NRFLD,TITLEB
    NTEXX = 1
    write(96,1012) NTEXX,SCALR,CONVR
  !
  ! WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 1, DATA OPTIONS:
  !  NJUN = NUMBER OF SEGMENTS CONNECTED BY FLOWS FROM THE HYD. FILE
  !  NCHN = NUMBER OF INTERFACIAL FLOW PAIRS FROM THE HYD. FILE
  !  DTWASP = WASP5 TIME STEP (SECONDS)
  !  TZERO = BEGIN TIME STEP FOR HYD. FILE (SECONDS)
  !  TENDHYD = END TIME STEP FOR HYD. FILE (SECONDS)
  !  ISTMP = CONTROL SWITCH, 0 = TIME VARIABLE SEGMENT DEPTHS AND VELOCITIES
  !          ARE READ; 1 = TIME VARIABLE SEGMENT DEPTHS AND VELOCITIES ARE N
  !          READ.
  !        NCHNC(KL) = 0
  !         LCHNC(KL,M) = 0
  !
    NJUN = KC*(LCLT-2)
    NCHNH = 0
    NCHNV = 0
  !
  !!!!!!!!!CHANGES NEXT 13 LINES
  !
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      NCHNH = NCHNH+INT(SUBO(L))
      if( IJCTLT(I+1,J) == 8 )then
        if( SUBO(LEC(L)) == 1.) NCHNH = NCHNH+1
      endif
      NCHNH = NCHNH+INT(SVBO(L))
      if( IJCTLT(I,J+1) == 8 )then
        if( SVBO(LNC(L)) == 1.) NCHNH = NCHNH+1
      endif
      NCHNV = NCHNV+INT(SWB(L))
    enddo
    NCHN = KC*NCHNH+KS*NCHNV
    ISTMP = 0
    NODYN = 1
    NODYN = NODYN
    DTWASP = DT * FLOAT(NTSMMT)
    TZERO = TBEGIN*TCON
    TENDHYD = TZERO+NTS*DT
    write(90,901)NJUN,NCHN
    if( IQOPT == 3 )then
      write(94,941) NJUN,NCHN, DTWASP, TZERO,TENDHYD,ISTMP
    endif
    if( IQOPT == 4 )then
      write(95) NJUN,NCHN, DTWASP, TZERO,TENDHYD,ISTMP
    endif
    write(96,1013) NCHN
  !
  ! ***  CHANNEL DATA
  ! WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 2, SEGMENT INTERFACE PAIRS:
  !   WASP EXPECTS TO SEE BOUNDARY SEGMENTS DESIGNATED AS "0".
  !
    RMNDUM = 0.
    LCHN = 0
    do K = KC,1,-1
      KMUL = KC-K
  !
  !!!!!!!!!!!!!!!CHANGES ON NEXT 38 LINES
  !
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        if( SUBO(L) == 1. )then
          LDTM = LT-1+KMUL*LCLTM2
          LUTM = LDTM-1
          if( IJCTLT(I-1,J) == 8 ) LUTM = 0
          RLENTH = DXU(L)
          WIDTH = DYU(L)
          LCHN = LCHN+1
  !
  !             LCHNC(LDTM,NCHNC(LDTM)) = LCHN
  !             LCHNC(LUTM,NCHNC(LUTM)) = LCHN
  !
          if( ISDHD == 1 ) WRITE(90,902)LCHN,RLENTH,WIDTH, &
              RMNDUM,LUTM,LDTM
          if( IQOPT == 3 ) WRITE(94,941) LUTM,LDTM
          if( IQOPT == 4 ) WRITE(95) LUTM,LDTM
          write(96,1014) UNITY,UNITY,LUTM,LDTM
        endif
        if( IJCTLT(I+1,J) == 8 )then
          if( SUBO(LEC(L)) == 1. )then
            LDTM = 0
            LUTM = LT-1+KMUL*LCLTM2
            RLENTH = DXU(LEC(L))
            WIDTH = DYU(LEC(L))
            LCHN = LCHN+1
  !
  !               LCHNC(LDTM,NCHNC(LDTM)) = LCHN
  !               LCHNC(LUTM,NCHNC(LUTM)) = LCHN
  !
            if( ISDHD  ==  1) WRITE(90,902) LCHN,RLENTH,WIDTH, &
                RMNDUM,LUTM,LDTM
            if( IQOPT == 3 ) WRITE(94,941) LUTM,LDTM
            if( IQOPT == 4 ) WRITE(95) LUTM,LDTM
            UNITY = 1.0
            write(96,1014) UNITY,UNITY,LUTM,LDTM
          endif
        endif
      enddo
  !
  !!!!!!!!!CHANGES NEXT 41 LINES
  !
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        if( SVBO(L) == 1. )then
          LSLT = LSCLT(LT)
          LDTM = LT-1+KMUL*LCLTM2
          LUTM = LSLT-1+KMUL*LCLTM2
          if( IJCTLT(I,J-1) == 8 ) LUTM = 0
          RLENTH = DYV(L)
          WIDTH = DXV(L)
          LCHN = LCHN+1
  !
  !             LCHNC(LDTM,NCHNC(LDTM)) = LCHN
  !             LCHNC(LUTM,NCHNC(LUTM)) = LCHN
  !
          if( ISDHD  ==  1) WRITE(90,902) LCHN,RLENTH,WIDTH, &
              RMNDUM,LUTM,LDTM
          if( IQOPT == 3 ) WRITE(94,941) LUTM,LDTM
          if( IQOPT == 4 ) WRITE(95) LUTM,LDTM
          write(96,1014) UNITY,UNITY,LUTM,LDTM
        endif
        if( IJCTLT(I,J+1) == 8 )then
          LN = LNC(L)
          if( SVBO(LN) == 1. )then
            LSLT = LSCLT(LT)
            LDTM = 0
            LUTM = LT-1+KMUL*LCLTM2
            RLENTH = DYV(LN)
            WIDTH = DXV(LN)
            LCHN = LCHN+1
  !
  !               LCHNC(LDTM,NCHNC(LDTM)) = LCHN
  !               LCHNC(LUTM,NCHNC(LUTM)) = LCHN
  !
            if( ISDHD  ==  1) WRITE(90,902) LCHN,RLENTH,WIDTH, &
                RMNDUM,LUTM,LDTM
            if( IQOPT == 3 ) WRITE(94,941) LUTM,LDTM
            if( IQOPT == 4 ) WRITE(95) LUTM,LDTM
            write(96,1014) UNITY,UNITY,LUTM,LDTM
          endif
        endif
      enddo
    enddo
    if( KC > 1 )then
      do K = KS,1,-1
        KMUL1 = KS-K
        KMUL2 = KMUL1+1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SWB(L) == 1. )then
            LUTM = LT-1+KMUL1*LCLTM2
            LDTM = LT-1+KMUL2*LCLTM2
            RLENTH = HLPF(L)*DZG(L,K)
            WIDTH = SQRT(DXYP(L))
            LCHN = LCHN+1
  !
  !               LCHNC(LDTM,NCHNC(LDTM)) = LCHN
  !               LCHNC(LUTM,NCHNC(LUTM)) = LCHN
  !
            write(90,902)LCHN,RLENTH,WIDTH,RMNDUM,LUTM,LDTM
            if( IQOPT == 3 ) WRITE(94,941)LUTM,LDTM
            if( IQOPT == 4 ) WRITE(95) LUTM,LDTM
            write(96,1014) UNITY,UNITY,LUTM,LDTM
          endif
        enddo
      enddo
  !
  ! WRITE OUT TIME SERIES OF ZERO DISPERSION COEFFICIENTS:
  !
      D1 = 0.0
      T1 = TZERO/TCON
      D2 = 0.0
      T2 = TENDHYD/TCON
      NBRKQ = 2
      write(96,905) NBRKQ
      write(96,1016) D1,T1, D2,T2
  !
  ! FOR EXCHANGE BETWEEN THE LOWER WATER SURFACE LAYER AND THE UPPER
  ! BENTHIC LAYER, DO THE FOLLOWING:
  !
      write(96,1012) NTEXX,SCALR,CONVR
      NTEXX = 0
      do K = 1,1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SWB(L) == 1. )then
            NTEXX = NTEXX+1
          endif
        enddo
      enddo
      write(96,1013) NTEXX
      do K = 1,1
        KMUL1 = KS-K
        KMUL2 = KMUL1+1
        KMUL3 = KMUL2+1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SWB(L) == 1. )then
            LUTM = LT-1+KMUL2*LCLTM2
            LDTM = LT-1+KMUL3*LCLTM2
            write(96,1014) DXYP(L),DEPSED,LUTM,LDTM
          endif
        enddo
      enddo
  !
  ! WRITE OUT TIME SERIES OF WATER-BENTHIC EXCHANGE DISPERSION COEFFICIENT
  !
      D1 = SEDIFF
      T1 = TZERO/TCON
      D2 = SEDIFF
      T2 = TENDHYD/TCON
      NBRKQ = 2
      write(96,905) NBRKQ
      write(96,1016) D1,T1, D2,T2
      IBPTMP = 0
      write(96,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
          IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
          IBPTMP,IBPTMP,IBPTMP,IBPTMP, &
          IBPTMP,IBPTMP,IBPTMP,IBPTMP
    endif
  !
  ! ***  JUNCTION DATA WITH INITIAL CONDITIONS
  ! WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 3, INITIAL SEGMENT PROPERTIE
  !
    VELTMP = 0.
    DUMVOL = 0.
    do K = KC,1,-1
      KMUL = KC-K
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        LCELL = LT-1+KMUL*LCLTM2
        DEPTMP = HLPF(L)*DZC(L,K)
        VOLTMP = DEPTMP*DXYP(L)
        if( ISDHD  ==  1) WRITE(90,904) LCELL,VOLTMP,I,J
        if( IQOPT == 3 ) WRITE(94,9440) VOLTMP,DUMVOL,DEPTMP,VELTMP
        if( IQOPT == 4 ) WRITE(95) VOLTMP,DEPTMP,VELTMP
      enddo
    enddo
    close(90)
    if( IQOPT == 3 ) close(94)
    if( IQOPT == 4 ) close(95)
    close(96)
  endif
  !
  ! ***  WRITE TIME STEP, VOLUME AND FLOW DATA
  !
  open(90,FILE = OUTDIR//'wasp\WASPDHD.OUT',POSITION = 'APPEND'   ,STATUS = 'UNKNOWN')
  if( IQOPT == 3 )then
    open(94,FILE = OUTDIR//'wasp\WASPDH.OUT',POSITION = 'APPEND'  ,STATUS = 'UNKNOWN')
  endif
  if( IQOPT == 4 )then
    open(95,FILE = OUTDIR//'wasp\WASPDHU.OUT',POSITION = 'APPEND' ,STATUS = 'UNKNOWN',  FORM = 'UNFORMATTED')
  endif
  LCLTM2 = LCLT-2
  IZERO = 0
  RZERO = 0
  IZERO = IZERO
  RZERO = RZERO
  !
  ! WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 4, BQ(J) FLOW IN INTERFACE
  ! PAIR "J":
  !
  LCHNUM = 0
  do K = KC,1,-1
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
  !
  ! +++++ FOLLOWING LINES BY M. MORTON TO INPUT DISPERSION TO HYD FILE:
  !
      ADDLW = 0.0
      if( SUB(L) == 1. )then
        LW = LWC(L)
        ADDLW = DYU(L)*AHULPF(L,K)*DZC(L,K)*0.5*(HLPF(L) &
            +HLPF(LW))*DXIU(L)
      endif
  !
  ! +++++ ABOVE ADDED BY M. MORTON
  !!!!!!!!!CHANGES NEXT 12 LINES
  !
      if( SUBO(L) == 1. )then
        FLOWX = DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(L,K)
        IMTMP = I-1
        LCHNUM = LCHNUM+1
        if( ISDHD  ==  1) WRITE(90,944) FLOWX,IMTMP,I,J,K
        if( IQOPT == 3 ) WRITE(94,946) FLOWX, ADDLW
        if( IQOPT == 4 ) WRITE(95) FLOWX, ADDLW
      endif
      if( IJCTLT(I+1,J) == 8 )then
        if( SUBO(LEC(L)) == 1. )then
          FLOWX = DYU(LEC(L))*(UHLPF(LEC(L),K)+SVPT*UVPT(LEC(L),K))*DZC(L,K)
          IPTMP = I+1
          LCHNUM = LCHNUM+1
          if( ISDHD  ==  1) WRITE(90,944) LCHNUM,FLOWX,I,IPTMP,J,K
          if( IQOPT == 3 ) WRITE(94,946) FLOWX, ADDLW
          if( IQOPT == 4 ) WRITE(95) FLOWX, ADDLW
        endif
      endif
    enddo
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
  !
  ! +++++ FOLLOWING LINES BY M. MORTON TO INPUT DISPERSION TO HYD FILE:
  !
      ADDLS = 0.0
      if( SVB(L) == 1. )then
        LS = LSC(L)
        ADDLS = DXV(L)*AHVLPF(L,K)*DZC(L,K)*0.5*(HLPF(L) &
            +HLPF(LS))*DYIV(L)
      endif
  !
  ! +++++ ABOVE ADDED BY M. MORTON
  !!!!!!!!CHANGES NEXT 13 LINES
  !
      if( SVBO(L) == 1. )then
        FLOWY = DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(L,K)
        JMTMP = J-1
        LCHNUM = LCHNUM+1
        if( ISDHD  ==  1) WRITE(90,944) LCHNUM,FLOWY,I,JMTMP,J,K
        if( IQOPT == 3 ) WRITE(94,946) FLOWY, ADDLS
        if( IQOPT == 4 ) WRITE(95) FLOWY, ADDLS
      endif
      if( IJCTLT(I,J+1) == 8 )then
        LN = LNC(L)
        if( SVBO(LN) == 1. )then
          FLOWY = DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(L,K)
          JPTMP = J+1
          LCHNUM = LCHNUM+1
          if( ISDHD  ==  1) WRITE(90,944) LCHNUM,FLOWY,I,J,JPTMP,K
          if( IQOPT == 3 ) WRITE(94,946) FLOWY, ADDLS
          if( IQOPT == 4 ) WRITE(95) FLOWY, ADDLS
        endif
      endif
    enddo
  enddo
  if( KC > 1 )then
    do K = KS,1,-1
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
  !
  ! +++++ FOLLOWING LINES BY M. MORTON TO INPUT DISPERSION TO HYD FILE:
  !
        ADDL = 0.0
        if( SPB(L) == 1. )then
          ADDL = DXYP(L)*ABLPF(L,K)*DZIG(L,K)
        endif
  !
  ! +++++ ABOVE ADDED BY M. MORTON
  !
        if( SWB(L) == 1 )then
          FLOWZ = -DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
          KPTMP = K+1
          LCHNUM = LCHNUM+1
          if( ISDHD  ==  1) WRITE(90,944) LCHNUM,FLOWZ,I,J,K,KPTMP
          if( IQOPT == 3 ) WRITE(94,946) FLOWZ, ADDL
          if( IQOPT == 4 ) WRITE(95) FLOWZ, ADDL
        endif
      enddo
    enddo
  endif
  !
  ! WRITE WASP5 HYDRODYNAMIC FILE DATA RECORD 5, SEGMENT PROPERTIES:
  !
  QQSUM = 0.
  LCELTMP = 0
  do K = KC,1,-1
    do LT = 2,LALT
      LCELTMP = LCELTMP+1
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      LN = LNC(L)
      VOLUM = DXYP(L)*HLPF(L)*DZC(L,K)
      DEPTH = HLPF(L)*DZC(L,K)
      VELX = 0.5*(UHLPF(L,K)+SVPT*UVPT(L,K) &
          +UHLPF(LEC(L),K)+SVPT*UVPT(LEC(L),K))/HLPF(L)
      VELY = 0.5*(VHLPF(L,K)+SVPT*VVPT(L,K) &
          +VHLPF(LN,K)+SVPT*VVPT(LN,K))/HLPF(L)
      VELZ = 0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1) &
          +WLPF(L,K)+SVPT*WVPT(L,K))
      VELMAG = SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
      if( ISDHD  ==  1) WRITE(90,902) LCELTMP,VOLUM,I,J,K
      if( IQOPT == 3 ) WRITE(94,946) VOLUM,QQSUM,DEPTH,VELMAG
      if( IQOPT == 4 ) WRITE(95) VOLUM, DEPTH, VELMAG
    enddo
  enddo
  close(90)
  if( IQOPT == 3 ) close(94)
  if( IQOPT == 4 ) close(95)
    901 FORMAT(2I5,E12.4,4I5,E12.4)
    902 FORMAT(I5,2X,3F20.8,3I5)
    903 FORMAT(3E12.4,2I5)
    904 FORMAT(I5,2X,F20.8,10I5)
    905 FORMAT(I5)
    906 FORMAT(5E12.4)
    941 FORMAT(2I5,3F20.8,I5)
    942 FORMAT(3E12.4,2I5)
    943 FORMAT(3E12.4,2I5)
    944 FORMAT(I5,2X,F20.8,10I5)
   9440 FORMAT(4F20.8)
    945 FORMAT(I5)
    946 FORMAT(4E17.9)
  JSWASP = 0
  return

END

