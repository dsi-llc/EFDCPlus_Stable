! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CEQICM

  ! ***  SUBROUTINE FOR INTERFACING CE-QUAL-ICM MODEL

  ! CHANGE RECORD

  use GLOBAL
  use INFOMOD,ONLY : SKIPCOM, READSTR
  use Variables_WQ,ONLY : HWQ, H2WQ
  
  implicit none
  
  integer :: NSKIP,IDMPCL,JDMPCL,M,IDUM,JDUM,MDUM,NS,NCTL,NWR,NMD,K,KICM,L,LM
  integer :: NFICMS,LTMP,NMID,NQSTMP,LE,LN,MM
  real    :: TIME,SVPT,TIMMID,TAVGTMP,TMPVAL
  character*80 STR*200

  real,save,allocatable,dimension(:)    :: QINRCA
  real,save,allocatable,dimension(:)    :: TMPICMF
  real,save,allocatable,dimension(:,:)  :: QINTFL
  
  integer,save,allocatable,dimension(:) :: IDRICM

  if( .not. allocated(QINRCA) )then
    allocate(QINRCA(NQSIJM))
    allocate(QINTFL(NQINFLM,KCM))
    allocate(TMPICMF(2*LCM*KCM))
    allocate(IDRICM(2*LCM*KCM))
    QINRCA = 0.
    QINTFL = 0.
    TMPICMF = 0.
    IDRICM = 0
    LFEFDC = 0
  endif

  ! ***  READ CONTROL FILES AND WRITE TIME INVARIANT FILES ON FIRST ENTRY
  if( JSWASP == 0 ) GOTO 1000
  TIME = TIMESEC/86400.              ! *** Time is in days

  !
  ! ***  READ MAIN CONTROL FILE
  !
  write(*,'(A)')'READING EFDC.ICM'
  open(1,FILE = 'efdc.icm',STATUS = 'UNKNOWN')
  STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
  read(1,*)ISDICM,IAUXICM,ISTICM,IDMPCL,JDMPCL,NCICM,NFICM
  close(1)

  ! ***  READ CELL MAPPING FILE
  write(*,'(A)')'READING EFDC_C_ICM.INP'
  open(1,FILE = 'efdc_c_icm.inp',STATUS = 'UNKNOWN')
  STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
  do M = 1,NCICM
    read(1,*)IDUM,JDUM,LCEFDC(M),KCEFDC(M),MDUM
  enddo
  close(1)

  ! ***  READ FLOW MAPPING FILE
  write(*,'(A)')'READING EFDC_F_ICM.INP'
  open(1,FILE = 'efdc_f_icm.inp',STATUS = 'UNKNOWN')
  STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
  do M = 1,NFICM
    read(1,*) IDUM, JDUM, LFEFDC(M), KFEFDC(M), MDUM, IDRICM(M)
  enddo
  close(1)

  ! ***  WRITE I,J INDICES DEFINING FLOWS BETWEEN ARBITARY CELLS
  ! ***  (POSTIVE FLOW DIRECTION DEFINED FROM FIRST TO SECOND I,J PAIR)
  if( IAUXICM >= 1 )then
    open(1,FILE = OUTDIR//'flwmap.inp',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'flwmap.inp',STATUS = 'UNKNOWN')
    NINTFL = 0
    if( NQSIJ > 0 )then
      do NS = 1,NQSIJ
        NINTFL = NINTFL+1
        write(1,101)IDMPCL,JDMPCL,BCPS(NS).I,BCPS(NS).J
      enddo
    endif
    if( NQCTL > 0 )then
      do NCTL = 1,NQCTL
        NINTFL = NINTFL+1
        if( HYD_STR(NCTL).IQCTLD > 0 )then
          write(1,101) HYD_STR(NCTL).IQCTLU, HYD_STR(NCTL).JQCTLU, &
                       HYD_STR(NCTL).IQCTLD, HYD_STR(NCTL).JQCTLD
        else
          write(1,101) HYD_STR(NCTL).IQCTLU,HYD_STR(NCTL).JQCTLU, &
              IDMPCL,JDMPCL
        endif
      enddo
    endif
    if( NQWR > 0 )then
      do NWR = 1,NQWR
        NINTFL = NINTFL+1
        write(1,101)WITH_RET(NWR).IQWRU,WITH_RET(NWR).JQWRU,WITH_RET(NWR).IQWRD,WITH_RET(NWR).JQWRD
      enddo
    endif
    if( MDCHH > 0 )then
      do NMD = 1,MDCHH
        if( IMDCHU(NMD) > 1 )then
          NINTFL = NINTFL+1
         write(1,101)IMDCHU(NMD),JMDCHU(NMD),IMDCHH(NMD),JMDCHH(NMD)
        endif
      enddo
      do NMD = 1,MDCHH
        if( IMDCHV(NMD) > 1 )then
          NINTFL = NINTFL+1
         write(1,101)IMDCHV(NMD),JMDCHV(NMD),IMDCHH(NMD),JMDCHH(NMD)
        endif
      enddo
    endif
    close(1)
  endif

  ! ***  WRITE EXTERNAL INFLOW LOCATIONS
  open(1,FILE = OUTDIR//'INFLOWIJ.DAT',STATUS = 'UNKNOWN')
  close(1,STATUS = 'DELETE')
  open(1,FILE = OUTDIR//'INFLOWIJ.DAT',STATUS = 'UNKNOWN')
  write(1,110)
  write(1,111)
  if( NQSIJ >= 1 )then
    do NS = 1,NQSIJ
      write(1,112)NS,BCPS(NS).I,BCPS(NS).J
    enddo
  endif
  close(1)

  ! ***  INITIALIZE INTERFACE FILES
  open(1,FILE = OUTDIR//'efdchyd.inp',FORM = 'UNFORMATTED',STATUS = 'UNKNOWN')
  close(1,STATUS = 'DELETE')
  open(1,FILE = OUTDIR//'efdchyd.inp',FORM = 'UNFORMATTED',STATUS = 'UNKNOWN')
  if( ISDICM == 1 )then
    open(2,FILE = OUTDIR//'EFDCHYD.ASC',STATUS = 'UNKNOWN')
    close(2,STATUS = 'DELETE')
    open(2,FILE = OUTDIR//'EFDCHYD.ASC',STATUS = 'UNKNOWN')
  endif

  ! ***  WRITE SIGMA GRID FRACTIONAL LAYER THICKNESS (SURFACE DOWN FOR ICM)
  write(1) (DZC(L,K),K = KC,1,-1)
  if( ISDICM == 1 )then
    write(2,2001)TIME
    write(2,2002)
    do K = KC,1,-1
      KICM = KC-K+1
      write(2,293)KICM,DZC(L,K),K
    enddo
  endif
  !
  ! ***  WRITE CELL HORIZONTAL SURFACE AREAS (LOOP OVER SURFACE LAYER CELLS
  !
  write(1) (DXYP(L),L = 2,LA)
  if( ISDICM == 1 )then
    write(2,2003)
    do L = 2,LA
      LM = LWC(L)
      write(2,293)LM,DXYP(L),L
    enddo
  endif
  !
  ! ***  WRITE CELL HORIZONTAL DEPTH AVERAGE FLOW FACE AREAS
  !
  NFICMS = NFICM/KC
  do M = 1,NFICMS
    L = LFEFDC(M)
    TMPICMF(M) = 0.0  ! *** DSI SINGLE LINE
    if( IDRICM(M) == 1 ) TMPICMF(M) = DYU(L)*HMU(L)
    if( IDRICM(M) == 2 ) TMPICMF(M) = DXV(L)*HMV(L)
  enddo
  write(1) (TMPICMF(M),M = 1,NFICMS)
  if( ISDICM == 1 )then
    write(2,2004)
    do M = 1,NFICMS
      write(2,293)M,TMPICMF(M),IDRICM(M),LFEFDC(M)
    enddo
  endif
  !
  ! ***  WRITE CELL HORIZONTAL X DIRECTION LENGTHS
  !
  write(1) (DXP(L),L = 2,LA)
  if( ISDICM == 1 )then
    write(2,2005)
    do L = 2,LA
      LM = LWC(L)
      write(2,293)LM,DXP(L),L
    enddo
  endif
  !
  ! ***  WRITE CELL HORIZONTAL Y DIRECTION LENGTHS
  !
  write(1) (DYP(L),L = 2,LA)
  if( ISDICM == 1 )then
    write(2,2006)
    do L = 2,LA
      LM = LWC(L)
      write(2,293)LM,DYP(L),L
    enddo
  endif
  !
  ! ***  WRITE CELL HORIZONTAL DEPTH AVERAGE INITIAL VOLUMES
  !
  do L = 2,LA
    TVAR3C(LWC(L)) = DXYP(L)*H2WQ(L)
    HTMP(L) = H2WQ(L)
  enddo
  LTMP = LA-1
  write(1) (TVAR3C(L),L = 1,LTMP)
  if( ISDICM == 1 )then
    write(2,2007)
    do L = 2,LA
      LM = LWC(L)
      write(2,293)LM,TVAR3C(LM),L
    enddo
  endif
  close(1)
  if( ISDICM == 1 ) close(2)
  !
  ! ***  INITIALIZE OTHER FILES TO RECEIVE TIME VARYING DATA
  !
  if( IAUXICM == 1 )then
    open(1,FILE = OUTDIR//'INFLOW.DAT',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'hydrlgy.inp',FORM = 'UNFORMATTED' &
           ,STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    if( ISDICM == 1 )then
      open(1,FILE = OUTDIR//'HYDRLGY.ASC',STATUS = 'UNKNOWN')
      close(1,STATUS = 'DELETE')
    endif
    open(1,FILE = OUTDIR//'intflw.inp',FORM = 'UNFORMATTED' &
           ,STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    if( ISDICM == 1 )then
      open(1,FILE = OUTDIR//'INTFLW.ASC',STATUS = 'UNKNOWN')
      close(1,STATUS = 'DELETE')
    endif
  endif
  if( ISTICM == 1 )then
    open(1,FILE = OUTDIR//'salicm.inp',FORM = 'UNFORMATTED' &
           ,STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'temicm.inp',FORM = 'UNFORMATTED' &
           ,STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
  endif
  if( ISTICM == 1 .and. ISDICM == 1 )then
    open(1,FILE = OUTDIR//'SALICM.ASC',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'TEMICM.ASC',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
  endif
  open(1,FILE = OUTDIR//'efdcflw.inp',FORM = 'UNFORMATTED' &
           ,STATUS = 'UNKNOWN')
  close(1,STATUS = 'DELETE')
  open(1,FILE = OUTDIR//'efdcrme.inp',FORM = 'UNFORMATTED' &
           ,STATUS = 'UNKNOWN')
  close(1,STATUS = 'DELETE')
  if( ISDICM == 1 )then
    open(1,FILE = OUTDIR//'EFDCFLW.ASC',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'EFDCRME.ASC',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
  endif
  !
  ! ***  WRITE SUMMARY DATA TO LOG FILE
  !
  open(1,FILE = OUTDIR//'EFDCICM.LOG',STATUS = 'UNKNOWN')
  close(1,STATUS = 'DELETE')
  open(1,FILE = OUTDIR//'EFDCICM.LOG',STATUS = 'UNKNOWN')
  write(1,2001) TIME
  write(1,102)IC,JC,KC
  write(1,103)NINTFL
  write(1,104)IDMPCL,JDMPCL
  write(1,2002)
  write(1,2003)
  write(1,2004)
  write(1,2005)
  write(1,2006)
  write(1,2007)
  close(1)
  JSWASP = 0
  return
   1000 continue
  !
  ! ***  SET VECTOR POTENTIAL TRANSPORT SWITCH TO INCLUDE VECTOR POTENTIAL
  !
  SVPT = 0.
  if( NTSMMT >= NTSPTC) SVPT = 1.
  TIME = TIMESEC/86400.  

  NMID = N-(NTSMMT/2)
  TIMMID = (DT*FLOAT(NMID)+TCON*TBEGIN)/86400.
  !
  ! ***  WRITE TIME AT END OF AVERAGING PERIOD TO EFDCICM.LOG
  !
  open(1,FILE = OUTDIR//'EFDCICM.LOG',STATUS = 'UNKNOWN' &
       ,POSITION = 'APPEND')
  write(1,106)TIME
  write(1,2008)
  write(1,2009)
  write(1,2010)
  write(1,2011)
  close(1)
  !
  ! ***  WRITE INFLOWS AT END OF AVERAGING PERIOD TO INFLOW.DAT
  !
  if( IAUXICM >= 1 )then
    open(1,FILE = OUTDIR//'INFLOW.DAT',STATUS = 'UNKNOWN' &
          ,POSITION = 'APPEND')
    do NS = 1,NQSIJ
      QINRCA(NS) = 0.
    enddo
    do NS = 1,NQSIJ
      NQSTMP = BCPS(NS).NQSERQ
      if( NQSTMP > 0 )then
        do K = 1,KC
          QINRCA(NS) = QINRCA(NS)+QSRTLPP(K,NQSTMP)+MAX(QSS(K,NS),0.)
        enddo
      else
        do K = 1,KC
          QINRCA(NS) = QINRCA(NS)+MAX(QSS(K,NS),0.)
        enddo
      endif
    enddo
    write(1,120)TIME,(QINRCA(NS),NS = 1,NQSIJ)
    close(1)
  endif
  !
  ! ***  WRITE INTERNAL FLOWS TO INTFLW.INP
  !
  if( IAUXICM >= 1 )then
    open(1,FILE = OUTDIR//'intflw.inp',FORM = 'UNFORMATTED' &
          ,STATUS = 'UNKNOWN', POSITION = 'APPEND')
    if( ISDICM == 1 )then
      open(2,FILE = OUTDIR//'INTFLW.ASC',STATUS = 'UNKNOWN' &
            ,POSITION = 'APPEND')
      write(2,106)TIME
    endif
    do K = 1,KC
      do NS = 1,NCRCA1
        QINTFL(NS,K) = 0.
      enddo
    enddo
    NINTFL = 0
    if( NQSIJ > 0 )then
      do NS = 1,NQSIJ
        NINTFL = NINTFL+1
        NQSTMP = BCPS(NS).NQSERQ
        if( NQSTMP > 0 )then
          do K = 1,KC
            QINTFL(NINTFL,K) = QINTFL(NINTFL,K)+QSRTLPN(K,NQSTMP) &
                +MIN(QSS(K,NS),0.)
          enddo
        else
          do K = 1,KC
            QINTFL(NINTFL,K) = QINTFL(NINTFL,K)+MIN(QSS(K,NS),0.)
          enddo
        endif
      enddo
    endif
    if( NQCTL > 0 )then
      do NCTL = 1,NQCTL
        NINTFL = NINTFL+1
        do K = 1,KC
          QINTFL(NINTFL,K) = QINTFL(NINTFL,K)+QCTLTLP(K,NCTL)
        enddo
      enddo
    endif
    if( NQWR > 0 )then
      do NWR = 1,NQWR
        NQSTMP = WITH_RET(NWR).NQWRSERQ
        NINTFL = NINTFL+1
        do K = 1,KC
          if( K == WITH_RET(NWR).KQWRU )then
            QINTFL(NINTFL,K) = QINTFL(NINTFL,K)+QWRSERTLP(NQSTMP) &
                +WITH_RET(NWR).QWR
          endif
        enddo
      enddo
    endif
    if( MDCHH > 0 )then
      do NMD = 1,MDCHH
        if( IMDCHU(NMD) > 1 )then
          NINTFL = NINTFL+1
          do K = 1,KC
            QINTFL(NINTFL,K) = QINTFL(NINTFL,K)+QCHNULP(NMD)*DZC(L,K)
          enddo
        endif
      enddo
      do NMD = 1,MDCHH
        if( IMDCHV(NMD) > 1 )then
          NINTFL = NINTFL+1
          do K = 1,KC
            QINTFL(NINTFL,K) = QINTFL(NINTFL,K)+QCHNULP(NMD)*DZC(L,K)
          enddo
        endif
      enddo
    endif
    write(1) QINTFL
    if( ISDICM == 1 )then
      write(2,213)
      do NS = 1,NINTFL
        write(2,211)NS,(QINTFL(NS,K),K = 1,KC)
      enddo
    endif
    close(1)
    if( ISDICM == 1 ) close(2)
  endif
  !
  ! ***  WRITE HYDRODYNAMIC DATA TO EFDCHDY.INP
  !
  open(1,FILE = OUTDIR//'efdchyd.inp',FORM = 'UNFORMATTED' &
      ,STATUS = 'UNKNOWN',  POSITION = 'APPEND')
  write(1) TIME
  open(3,FILE = OUTDIR//'efdcflw.inp',FORM = 'UNFORMATTED' &
      ,STATUS = 'UNKNOWN',  POSITION = 'APPEND')
  write(3) TIME
  open(13,FILE = OUTDIR//'efdcrme.inp',FORM = 'UNFORMATTED' &
      ,STATUS = 'UNKNOWN',  POSITION = 'APPEND')
  write(13) TIME
  if( ISDICM == 1 )then
    open(2,FILE = OUTDIR//'EFDCHYD.ASC',STATUS = 'UNKNOWN' &
          ,POSITION = 'APPEND')
    write(2,106) TIME
    open(4,FILE = OUTDIR//'EFDCFLW.ASC',STATUS = 'UNKNOWN' &
          ,POSITION = 'APPEND')
    write(4,106) TIME
    open(14,FILE = OUTDIR//'EFDCRME.ASC',STATUS = 'UNKNOWN' &
           ,POSITION = 'APPEND')
    write(14,106) TIME
  endif
  TAVGTMP = FLOAT(NTSMMT)*DT
  !
  ! ***  WRITE CELL HORIZONTAL DEPTH AVERAGE VOLUMES
  !
  do L = 1,LC
    TVAR3E(L) = 0.
    TVAR3W(L) = 0.
    TVAR3S(L) = 0.
    TVAR3N(L) = 0.
    TVAR3C(L) = 0.
  enddo
  do K = 1,KC
    do L = 2,LA
      TVAR3W(L) = TVAR3W(L)+DYU(L)*UHLPF(L,K)
      TVAR2W(L,K) = DYU(L)*UHLPF(L,K)
    enddo
  enddo
  do K = 1,KC
    do L = 2,LA
      LE = LEC(L)
      TVAR3E(L) = TVAR3E(L)+DYU(LE)*UHLPF(LE,K)
      TVAR2E(L,K) = DYU(LE)*UHLPF(LE,K)
    enddo
  enddo
  do K = 1,KC
    do L = 2,LA
      TVAR3S(L) = TVAR3S(L)+DXV(L)*VHLPF(L,K)
      TVAR2S(L,K) = DXV(L)*VHLPF(L,K)
    enddo
  enddo
  do K = 1,KC
    do L = 2,LA
      LN = LNC(L)
      TVAR3N(L) = TVAR3N(L)+DXV(LN)*VHLPF(LN,K)
      TVAR2N(L,K) = DXV(LN)*VHLPF(LN,K)
    enddo
  enddo
  do K = 1,KC
    do L = 2,LA
      TVAR3N(L) = TVAR3N(L)*DZC(L,K)
      TVAR3S(L) = TVAR3S(L)*DZC(L,K)
      TVAR3E(L) = TVAR3E(L)*DZC(L,K)
      TVAR3W(L) = TVAR3W(L)*DZC(L,K)
    enddo
  enddo
  do K = 1,KC
    do L = 2,LA
      TVAR2N(L,K) = TVAR2N(L,K)*DZC(L,K)
      TVAR2S(L,K) = TVAR2S(L,K)*DZC(L,K)
      TVAR2E(L,K) = TVAR2E(L,K)*DZC(L,K)
      TVAR2W(L,K) = TVAR2W(L,K)*DZC(L,K)
    enddo
  enddo
  TMPVAL = DT*FLOAT(NTSMMT)
  do L = 2,LA
    WLPF(L,0) = 0.
    WLPF(L,KC) = 0.
    TVAR2C(L,0) = 0.
  enddo
  do K = 1,KS
    do L = 2,LA
      WLPF(L,K) = WLPF(L,K-1)+TVAR2W(L,K)-TVAR2E(L,K) &
          +TVAR2S(L,K)-TVAR2N(L,K) &
          +QSUMLPF(L,K)
    enddo
  enddo
  do L = 2,LA
    TVAR3C(L) = WLPF(L,KS)+TVAR2W(L,KC)-TVAR2E(L,KC) &
        +TVAR2S(L,KC)-TVAR2N(L,KC) &
        +QSUMLPF(L,KC)
  enddo
  do K = 1,KS
    do L = 2,LA
      TVAR2C(L,K) = TVAR2C(L,K-1)+DZC(L,K)*TVAR3C(L)
    enddo
  enddo
  do K = 1,KS
    do L = 2,LA
      WLPF(L,K) = WLPF(L,K)-TVAR2C(L,K)
    enddo
  enddo
  do L = 2,LA
    TVAR3C(L) = TMPVAL*TVAR3C(L)
  enddo
  do L = 2,LA
    TVAR3C(L) = DXYP(L)*HTMP(L)+TVAR3C(L)
  enddo
  do L = 2,LA
    HTMP(L) = TVAR3C(L)*DXYIP(L)
  enddo
  write(1) (TVAR3C(L),L = 2,LA)
  do L = 2,LA
    QSUMLPF(L,KC) = QSUMLPF(L,KC)-RAINLPF(L)+EVPSLPF(L)
  enddo
  do K = KC,1,-1
    write(3) (QSUMLPF(L,K),L = 2,LA)
  enddo
  if( ISDICM == 1 )then
    write(2,2008)
    write(4,401)
    do L = 2,LA
      LM = LWC(L)
      TMPVAL = DXYP(L)*HWQ(L)
      write(2,2294)LM,L,HP(L),HWQ(L),HTMP(L)
      write(4,402)LM,L,(QSUMLPF(L,K),K = 1,KC)
  !
  !             RVAL2 = -DYU(LP)*UHLPF(LP,K)*DZC(L,K)
  !             RVAL4 = -DXV(LN)*VHLPF(LN,K)*DZC(L,K)
  !
    enddo
  endif
  do L = 2,LA
    TVAR3C(L) = RAINLPF(L)-EVPSLPF(L)
  enddo
  write(13) (TVAR3C(L),L = 2,LA)
  do L = 2,LA
    LM = LWC(L)
    write(14,402)LM,L,TVAR3C(L)
  enddo
    222 FORMAT(' ERROR ',2I5,6F12.2)
    401 FORMAT(/,'  LICM     L    QSUMLPF(L,K) K = 1,KC',/)
    402 FORMAT(2I6,12E13.5)
   2294 FORMAT(2I6,4F12.5)
  !
  ! ***  WRITE HORIZONTAL FLOWS
  !
  do M = 1,NFICM
    L = LFEFDC(M)
    K = KFEFDC(M)
    if( IDRICM(M) == 1 )then
      TMPICMF(M) = DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(L,K)
    endif
    if( IDRICM(M) == 2 )then
      TMPICMF(M) = DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(L,K)
    endif
  enddo
  write(1) (TMPICMF(M),M = 1,NFICM)
  if( ISDICM == 1 )then
    write(2,2009)
    do M = 1,NFICM
      L = LFEFDC(M)
      K = KFEFDC(M)
      write(2,293)M,TMPICMF(M),IDRICM(M),L,K
    enddo
  endif
  !
  ! ***  WRITE VERTICAL DIFFUSIVITY
  !
  MM = 0
  do K = KS,1,-1
    do L = 2,LA
      MM = MM+1
      TMPICMF(MM) = ABLPF(L,K)
    enddo
  enddo
  write(1) (TMPICMF(M),M = 1,MM)
  if( ISDICM == 1 )then
    write(2,2010)
    MM = 0
    do K = KS,1,-1
      do L = 2,LA
        MM = MM+1
        write(2,293)MM,ABLPF(L,K),L,K
      enddo
    enddo
  endif
  !
  ! ***  DETERMINATION OF VERTICAL VELOCITY AND VOLUME FLUX
  ! ***  LOAD NET DEPTH INTEGRATED INFLOWS INTO TVAR3E AND
  ! ***  OUT FLOWS INTO TVAR3N
  !       WLPF(L,0) = 0.
  !        TVAR1N(L,K) = 0.
  !        TVAR1S(L,K) = 0.
  !        TVAR1E(L,K) = 0.
  !        TVAR1W(L,K) = 0.
  ! ***  CALCULATE QZ (STORED IN WLPF(L,K)
  ! ***  WRITE VERTICAL VOLUME FLUX
  !
  MM = 0
  do K = KS,1,-1
    do L = 2,LA
      MM = MM+1
      TMPICMF(MM) = WLPF(L,K)
    enddo
  enddo
  write(1) (TMPICMF(M),M = 1,MM)
  if( ISDICM == 1 )then
    write(2,2011)
    MM = 0
    do K = KS,1,-1
      do L = 2,LA
        MM = MM+1
        write(2,293)MM,WLPF(L,K),L,K
      enddo
    enddo
    MM = 0
    do L = 2,LA
      write(2,293)MM,WLPF(L,KC),L,KC
    enddo
  endif
  close(1)
  close(3)
  close(13)
  if( ISDICM == 1 ) close(2)
  if( ISDICM == 1 ) close(4)
  if( ISDICM == 1 ) close(14)
  !
  ! ***  WRITE SALINITY
  !
  if( ISTICM == 1 )then
    open(1,FILE = OUTDIR//'salicm.inp',FORM = 'UNFORMATTED' &
          ,POSITION = 'APPEND')
    MM = 0
    do K = KS,1,-1
      do L = 2,LA
        MM = MM+1
        TMPICMF(MM) = SALLPF(L,K)
      enddo
    enddo
    write(1) (TMPICMF(M),M = 1,MM)
    close(1)
  endif
  if( ISTICM == 1 .and. ISDICM == 1 )then
    open(1,FILE = OUTDIR//'SALICM.ASC',POSITION = 'APPEND')
    write(1,106) TIME
    write(1,2012)
    MM = 0
    do K = KS,1,-1
      do L = 2,LA
        MM = MM+1
        write(2,293)MM,SALLPF(L,K),L,K
      enddo
    enddo
    close(1)
  endif
  !
  ! ***  WRITE TEMPERATURE
  !
  if( ISTICM == 1 )then
    open(1,FILE = OUTDIR//'temicm.inp',FORM = 'UNFORMATTED' &
          ,POSITION = 'APPEND')
    MM = 0
    do K = KS,1,-1
      do L = 2,LA
        MM = MM+1
        TMPICMF(MM) = TEMLPF(L,K)
      enddo
    enddo
    write(1) (TMPICMF(M),M = 1,MM)
    close(1)
  endif
  if( ISTICM == 1 .and. ISDICM == 1 )then
    open(1,FILE = OUTDIR//'TEMICM.ASC',POSITION = 'APPEND')
    write(1,106) TIME
    write(1,2013)
    MM = 0
    do K = KS,1,-1
      do L = 2,LA
        MM = MM+1
        write(2,293)MM,TEMLPF(L,K),L,K
      enddo
    enddo
    close(1)
  endif
  !
  ! ***  WRITE HYDROLOGY
  !
  if( IAUXICM == 1 )then
    open(1,FILE = OUTDIR//'hydrlgy.inp',FORM = 'UNFORMATTED' &
          ,POSITION = 'APPEND')
    write(1)TIME
    MM = 0
    do L = 2,LA
      MM = MM+1
      TVAR3C(MM) = RAINLPF(L)
    enddo
    write(1) (TVAR3C(M),M = 1,MM)
    MM = 0
    do L = 2,LA
      MM = MM+1
      TVAR3C(MM) = EVPSLPF(L)
    enddo
    write(1) (TVAR3C(M),M = 1,MM)
    MM = 0
    do L = 2,LA
      MM = MM+1
      TVAR3C(MM) = RINFLPF(L)
    enddo
    write(1) (TVAR3C(M),M = 1,MM)
    close(1)
    if( ISDICM == 1 )then
      open(2,FILE = OUTDIR//'HYDRLGY.ASC',POSITION = 'APPEND')
      write(2,106)TIME
      write(2,212)
      do L = 2,LA
        write(2,200)L,IL(L),JL(L),RAINLPF(L),EVPSLPF(L),EVPGLPF(L), &
            RINFLPF(L),GWLPF(L)
      enddo
      close(2)
    endif
  endif
    100 FORMAT(120X)
    101 FORMAT(4I10)
    102 FORMAT(/,' NROW,NCOL,NLAYR = ',3I10/)
    103 FORMAT(/,' NO INTERNAL FLOWS, NINTFL (LINES) IN FLWMAP.INP = ', &
      I10/)
    104 FORMAT(/,' ROW, COLUMN INDICES OF DUMP CELL = ',2I10/)
    105 FORMAT(/,' SIMULATION STARTING TIME IN DAYS = ',F12.6/)
    106 FORMAT(/,' TIME IN DAYS AT END OF AVERAGING PERIOD = ',F12.6/)
    110 FORMAT(' LOCATION OF INFLOWS ',/)
    111 FORMAT(' INFLOW #   ROW INDEX   COLUMN INDEX ',/)
    112 FORMAT(2X,I5,7X,I5,7X,I5)
    120 FORMAT(F12.6,13F12.4)
    200 FORMAT(3I5,6E14.6)
    201 FORMAT(' L,I(ROW),J(COL),QX(I,J,K),K = 1,KC ',/)
    202 FORMAT(' L,I(ROW),J(COL),QY(I,J,K),K = 1,KC ',/)
    203 FORMAT(' L,I(ROW),J(COL),QZ(I,J,K),K = 1,KS ',/)
    204 FORMAT(' L,I(ROW),J(COL),AX(I,J,K),K = 1,KC ',/)
    205 FORMAT(' L,I(ROW),J(COL),AY(I,J,K),K = 1,KC ',/)
    206 FORMAT(' L,I(ROW),J(COL),AZ(I,J,K),K = 1,KS ',/)
    207 FORMAT(' L,I(ROW),J(COL),SELS(I,J),SELE(I,J),DSEL(I,J) ',/)
    208 FORMAT(' L,I(ROW),J(COL),SAL(I,J,K),K = 1,KC ',/)
    209 FORMAT(' L,I(ROW),J(COL),TEM(I,J,K),K = 1,KC ',/)
    210 FORMAT(//)
    211 FORMAT(I5,2X,6E15.6)
    212 FORMAT(' L,I(ROW),J(ROW),RAINLPF(I,J),EVPSLPF(I,J),EVPGLPF(I,J), &
      RINFLPF(I,J),GWLPF(I,J) ',/)
    213 FORMAT(' NQINTFL,QINTFL ',/)
    215 FORMAT(' L,I(ROW),J(COL),SURFELV START AVG INTERVAL',/)
    216 FORMAT(' L,I(ROW),J(COL),DEL SURFELV OVER INTERVAL',/)
    291 FORMAT(I8,F8.4)
    292 FORMAT(I8,E13.5)
    293 FORMAT(I8,E13.5,3I8)
    294 FORMAT(I8,E13.5,E13.5,3I8)
   2001 FORMAT(/,' TIME AT ICM INTERFACE INITIALIZATION = ',F12.4,/)
   2002 FORMAT(/,' SIGMA LAYER FRACTIONAL THICKNESS: KICM, DZ, KEFDC',/)
   2003 FORMAT(/,' HORIZONTAL CELL SURFACE AREAS, TOP LAYER : LICM, AREA' &
      ,' LEFDC',/)
   2004 FORMAT(/,' DEPTH INTEGRATED HORIZONTAL FLOW FACE AREAS: NF, AREA' &
      ,' IDIR, LEFDC',/)
   2005 FORMAT(/,' X DIRECTION CELL LENTHS, TOP LAYER: LICM, DX, LEFDC',/)
   2006 FORMAT(/,' Y DIRECTION CELL LENTHS, TOP LAYER: LICM, DY, LEFDC',/)
   2007 FORMAT(/,' INITIAL DEPTH AVERGE CELL VOLUMES: LICM, VOL, LEFDC',/)
   2008 FORMAT(/,' DEPTH AVERGE CELL VOLUMES: LICM, VOL, LEFDC',/)
   2009 FORMAT(/,' HORIZONTAL FLOWS: NF, FLOW, IDIR, LEFDC, KEFDC'/)
   2010 FORMAT(/,' VERTICAL DIFFUSIVITY: ND, DIFF, LEFDC, KEFDC',/)
   2011 FORMAT(/,' VERTICAL SIGMA FLOW: ND, FLOWZ, LEFDC, KEFDC',/)
   2012 FORMAT(/,' HYDRO SALINITY: NICM, SAL, LEFDC, KEFDC',/)
   2013 FORMAT(/,' HYDRO SALINITY: NICM, TEM, LEFDC, KEFDC',/)
  return
END

