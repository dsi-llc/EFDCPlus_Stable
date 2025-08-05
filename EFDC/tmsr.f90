! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE TMSR

  ! ***  SUBROUTINE TMSR WRITES TIME SERIES FILES FOR SURFACE ELEVATION,
  ! ***  VELOCITY, CONCENTRATION, AND VOLUME SOURCES AT SPECIFIED
  ! ***  (I,J) POINTS

  ! CHANGE RECORD

  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out
  
  implicit none
  
  integer :: MLTM,MSDIG,MTMP,MFDIG,NX,NT,NTSSS,MTSSS,L,NMD,LMDCHHT
  integer :: LMDCHUT,LMDCHVT,MTSCC,I,J,LN,K,KTMP,KTMP1,NSXD,NDOC,MD
  
  real    :: TAUW1,TAUW2,TIME,SEDBTMP,SNDBTMP,PPTMP,TMPVAL,UTMP1,VTMP1,UTMP,VTMP
  real    :: TBEAST,TBNORT,TAUBDYN,TAUB2,CURANG,TAUTOT,VELMAG,TMPDRAG,TAUBSEDDYN
  real    :: TAUBSNDDYN,QBEDLOADX,QBEDLOADY,CQBEDLOADXT,CQBEDLOADYT,RUVTMP,QRNT

  character*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7, &
      TITLE11,TITLE12,TITLE13,TITLE14,TITLE15,TITLE16,TITLE17, &
      TITLE18,TITLE19,TFNTXWT,TFNTXWF,TFNTXWC,TFNTXWP, &
      TFNTXBT,TFNTXBF,TFNTXBC,TFNTXBP
  character*10 CTUNIT
  character*1 CZTT(0:9)
  character*1 CCHTMF,CCHTMS
  character*2 CNTOX(25),CNSND(10),CNSBL(10)
  
  real,allocatable,dimension(:) :: ATMPP
  real,allocatable,dimension(:) :: BTMP          

  real,allocatable,dimension(:) :: PORH
  real,allocatable,dimension(:) :: QCHANUIJ
  real,allocatable,dimension(:) :: QCHANVIJ
  real,allocatable,dimension(:) :: SEDSND
  real,allocatable,dimension(:) :: SEDSNDB
  real,allocatable,dimension(:) :: TXBC
  real,allocatable,dimension(:) :: TXBF
  real,allocatable,dimension(:) :: TXBP
  real,allocatable,dimension(:) :: TXBT
  real,allocatable,dimension(:) :: TXWC
  real,allocatable,dimension(:) :: TXWF
  real,allocatable,dimension(:) :: TXWP

  character* 2,save,allocatable,dimension(:)     :: CNTMSR
  character*20,save,allocatable,dimension(:)     :: FNAVB
  character*20,save,allocatable,dimension(:)     :: FNAVV
  character*20,save,allocatable,dimension(:)     :: FNBED
  character*20,save,allocatable,dimension(:)     :: FNDOX
  character*20,save,allocatable,dimension(:)     :: FNDYE
  character*20,save,allocatable,dimension(:)     :: FNNHX
  character*20,save,allocatable,dimension(:)     :: FNQ3D
  character*20,save,allocatable,dimension(:)     :: FNQQE
  character*20,save,allocatable,dimension(:)     :: FNSAL
  character*22,save,allocatable,dimension(:,:)   :: FNSBL
  character*20,save,allocatable,dimension(:)     :: FNSED
  character*20,save,allocatable,dimension(:)     :: FNSEL
  character*20,save,allocatable,dimension(:)     :: FNSFL
  character*22,save,allocatable,dimension(:,:)   :: FNSND
  character*20,save,allocatable,dimension(:)     :: FNTEM
  character*20,save,allocatable,dimension(:)     :: FNTOC
  character*22,save,allocatable,dimension(:,:)   :: FNTOX
  character*23,save,allocatable,dimension(:,:)   :: FNTXBC
  character*23,save,allocatable,dimension(:,:)   :: FNTXBF
  character*23,save,allocatable,dimension(:,:)   :: FNTXBP
  character*23,save,allocatable,dimension(:,:)   :: FNTXBT
  character*23,save,allocatable,dimension(:,:)   :: FNTXWC
  character*23,save,allocatable,dimension(:,:)   :: FNTXWF
  character*23,save,allocatable,dimension(:,:)   :: FNTXWP
  character*23,save,allocatable,dimension(:,:)   :: FNTXWT
  character*20,save,allocatable,dimension(:)     :: FNU3D
  character*20,save,allocatable,dimension(:)     :: FNUVE
  character*20,save,allocatable,dimension(:)     :: FNUVT
  character*20,save,allocatable,dimension(:)     :: FNV3D

  allocate(ATMPP(KCM))
  allocate(BTMP(KCM))     
  allocate(PORH(KBM))
  allocate(QCHANUIJ(LCM))
  allocate(QCHANVIJ(LCM))
  allocate(SEDSND(KCM))
  allocate(SEDSNDB(KBM))
  allocate(TXBC(KBM))
  allocate(TXBF(KBM))
  allocate(TXBP(KBM))
  allocate(TXBT(KBM))
  allocate(TXWC(KCM))
  allocate(TXWF(KCM))
  allocate(TXWP(KCM))

  ! *** ALLOCATE LOCAL ARRAYS
  if( .not. allocated(CNTMSR) )then
    allocate(CNTMSR(MLTMSRM))
    allocate(FNAVB(MLTMSRM))
    allocate(FNAVV(MLTMSRM))
    allocate(FNBED(MLTMSRM))
    allocate(FNDOX(MLTMSRM))
    allocate(FNDYE(MLTMSRM))
    allocate(FNNHX(MLTMSRM))
    allocate(FNQ3D(MLTMSRM))
    allocate(FNQQE(MLTMSRM))
    allocate(FNSAL(MLTMSRM))
    allocate(FNSBL(MLTMSRM,NSNM))
    allocate(FNSED(MLTMSRM))
    allocate(FNSEL(MLTMSRM))
    allocate(FNSFL(MLTMSRM))
    allocate(FNSND(MLTMSRM,NSNM))
    allocate(FNTEM(MLTMSRM))
    allocate(FNTOC(MLTMSRM))
    allocate(FNTOX(MLTMSRM,NTXM))
    allocate(FNTXBC(MLTMSRM,NTXM))
    allocate(FNTXBF(MLTMSRM,NTXM))
    allocate(FNTXBP(MLTMSRM,NTXM))
    allocate(FNTXBT(MLTMSRM,NTXM))
    allocate(FNTXWC(MLTMSRM,NTXM))
    allocate(FNTXWF(MLTMSRM,NTXM))
    allocate(FNTXWP(MLTMSRM,NTXM))
    allocate(FNTXWT(MLTMSRM,NTXM))
    allocate(FNU3D(MLTMSRM))
    allocate(FNUVE(MLTMSRM))
    allocate(FNUVT(MLTMSRM))
    allocate(FNV3D(MLTMSRM))
    
    CNTMSR = ' '
    FNAVB = ' '
    FNAVV = ' '
    FNBED = ' '
    FNDOX = ' '
    FNDYE = ' '
    FNNHX = ' '
    FNQ3D = ' '
    FNQQE = ' '
    FNSAL = ' '
    FNSBL = ' '
    FNSED = ' '
    FNSEL = ' '
    FNSFL = ' '
    FNSND = ' '
    FNTEM = ' '
    FNTOC = ' '
    FNU3D = ' '
    FNUVE = ' '
    FNUVT = ' '
    FNV3D = ' '
  endif
  
  ! *** INITIALIZE LOCAL ARRAYS
  PORH = 0.
  QCHANUIJ = 0.
  QCHANVIJ = 0.
  SEDSND = 0.
  SEDSNDB = 0.
  TXBC = 0.
  TXBF = 0.
  TXBP = 0.
  TXBT = 0.
  TXWC = 0.
  TXWF = 0.
  TXWP = 0.

  if( JSTMSR /= 1 ) GOTO 300
  !
  if( MLTMSR > MLTMSRM )then
    WRITE (6,600)
    STOP
  endif
  if( MLTMSR > 99 )then
    WRITE (6,601)
    STOP
  endif
  !
    600 FORMAT(' NUMBER OF TIME SER LOC, MLTMSR, EXCEEDS DIM, MLTMSRM')
    601 FORMAT(' NUMBER OF TIME SERIES LOCATIONS EXCEED 99')
  !
        TAUW1 = 0.
        TAUW2 = 0.
  !
  CZTT(0) = '0'
  CZTT(1) = '1'
  CZTT(2) = '2'
  CZTT(3) = '3'
  CZTT(4) = '4'
  CZTT(5) = '5'
  CZTT(6) = '6'
  CZTT(7) = '7'
  CZTT(8) = '8'
  CZTT(9) = '9'
  !
  do MLTM = 1,MLTMSR_GL
      MSDIG = MOD(MLTM,10)
      MTMP = MLTM-MSDIG
      MFDIG = MTMP/10
      CCHTMF = CZTT(MFDIG)
      CCHTMS = CZTT(MSDIG)
      CNTMSR(MLTM)= CCHTMF // CCHTMS
  enddo
  !
  if( TCTMSR == 1.) CTUNIT = 'SECONDS'
  if( TCTMSR == 60.) CTUNIT = 'MINUTES'
  if( TCTMSR == 3600.) CTUNIT = 'HOURS'
  if( TCTMSR == 86400.) CTUNIT = 'DAYS'
  !
   CNTOX( 1)= '01'
   CNTOX( 2)= '02'
   CNTOX( 3)= '03'
   CNTOX( 4)= '04'
   CNTOX( 5)= '05'
   CNTOX( 6)= '06'
   CNTOX( 7)= '07'
   CNTOX( 8)= '08'
   CNTOX( 9)= '09'
   CNTOX(10)= '10'
   CNTOX(11)= '11'
   CNTOX(12)= '12'
   CNTOX(13)= '13'
   CNTOX(14)= '14'
   CNTOX(15)= '15'
   CNTOX(16)= '16'
   CNTOX(17)= '17'
   CNTOX(18)= '18'
   CNTOX(19)= '19'
   CNTOX(20)= '20'
   CNTOX(21)= '21'
   CNTOX(22)= '22'
   CNTOX(23)= '23'
   CNTOX(24)= '24'
   CNTOX(25)= '25'
  !
   CNSND( 1)= '01'
   CNSND( 2)= '02'
   CNSND( 3)= '03'
   CNSND( 4)= '04'
   CNSND( 5)= '05'
   CNSND( 6)= '06'
   CNSND( 7)= '07'
   CNSND( 8)= '08'
   CNSND( 9)= '09'
   CNSND(10)= '10'
  !
   CNSBL( 1)= '01'
   CNSBL( 2)= '02'
   CNSBL( 3)= '03'
   CNSBL( 4)= '04'
   CNSBL( 5)= '05'
   CNSBL( 6)= '06'
  !

  ! ***  WRITE HEADINGS
  !
  TITLE1 = ' SALINITY (PSU) TIME SERIES, K = 1,KC'
  TITLE2 = ' TEMPERATURE (DEG C) TIME SERIES, K = 1,KC'
  TITLE3 = ' DYE CONC (KG/M**3) TIME SERIES, K = 1,KC'
  TITLE4 = ' SED CONC (MG/LITER) TIME SERIES, K = 1,KC'
  TITLE5 = ' TOXIC CONC (M/TOT VOL - UG/LITER) 1-4 BED,5-8 WC'
  TITLE6 = ' VISCOSITY (CM**2/S) TIME SERIES, K = 1,KS'
  TITLE7 = ' DIFFUSIVITY (CM**2/S) TIME SERIES, K = 1,KS'
  TITLE11 = ' SURFACE ELEVATION & DEPTH (METERS) TIME SERIES'
  TITLE12 = ' EXT MODE E,N VEL (CM/S) TBX TBY TB TBCG TBNG (CM/S)**2'
  TITLE13 = ' EXT MODE U,V TRANSPORT (M**3/S) TIME SERIES'
  TITLE14 = ' INT MODE EAST VEL (CM/S) TIME SERIES, K = 1,KC'
  TITLE15 = ' INT MODE NORTH VEL (CM/S) TIME SERIES, K = 1,KC'
  TITLE16 = ' EXT MODE VOLUME S/S (M**3/S) TIME SERIES'
  TITLE17 = ' INT MODE VOL S/S (M**3/S) TIME SERIES, K = 1,KC'
  TITLE18 = ' SED BED LOAD QSX QSY (GM/S) CQSX CQSY (MG/L) '
  TITLE19 = ' BED TOP  KBT HBED(KBT) HBED(KBT-1) VOIDR FRAC SED/SND'
  TFNTXWT = ' TOTAL TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
  TFNTXWF = ' FREE DIS TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
  TFNTXWC = ' DOC COMP TOXIC CONC WATER COL (M/TOT VOL), UG/LITER'
  TFNTXWP = ' TOT PART TOXIC CONC WATER COL (M/M), UG/GM'
  TFNTXBT = ' TOTAL TOXIC CONC SED BED (M/TOT VOL), UG/LITER'
  TFNTXBF = ' FREE DIS TOXIC CONC SED BED (M/PORE VOL), UG/LITER'
  TFNTXBC = ' DOC COMP TOXIC CONC SED BED (M/PORE VOL), UG/LITER'
  TFNTXBP = ' TOT PART TOXIC CONC SED BED (M/M), UG/GM'
  !
  if( ISTMSR == 2 )then
  do MLTM = 1,MLTMSR
    if( MTMSRC(MLTM) == 1 )then
      if( ISTRAN(1) >= 1 )then
        FNSAL(MLTM) = OUTDIR//'SALTS' // CNTMSR(MLTM) // '.OUT'
      endif
      if( ISTRAN(2) >= 1 )then
        FNTEM(MLTM) = OUTDIR//'TEMTS' // CNTMSR(MLTM) // '.OUT'
      endif
      if( ISTRAN(3) >= 1 )then
        FNDYE(MLTM) = OUTDIR//'DYETS' // CNTMSR(MLTM) // '.OUT'
      endif
      if( ISTRAN(4) >= 1 )then
        FNSFL(MLTM) = OUTDIR//'SFLTS' // CNTMSR(MLTM) // '.OUT'
      endif
      if( ISTRAN(6) >= 1 )then
        FNSED(MLTM) = OUTDIR//'SEDTS' // CNTMSR(MLTM) // '.OUT'
      endif
  !          if( ISTRAN(7) >= 1 )then
  !            FNSND(MLTM) = OUTDIR//'SNDTS' // CNTMSR(MLTM) // '.OUT'
  !          endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
        FNSND(MLTM,NX) = OUTDIR//'SND' // CNSND(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNSBL(MLTM,NX) = OUTDIR//'SBL' // CNSBL(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
        enddo
      endif
      if( ISTRAN(8) >= 1 )then
        FNDOX(MLTM) = OUTDIR//'DOXTS' // CNTMSR(MLTM) // '.OUT'
        FNTOC(MLTM) = OUTDIR//'TOCTS' // CNTMSR(MLTM) // '.OUT'
        FNNHX(MLTM) = OUTDIR//'NHXTS' // CNTMSR(MLTM) // '.OUT'
      endif
      if( ISTRAN(5) >= 1 )then
        do NT = 1,NTOX
        FNTXWT(MLTM,NT) = OUTDIR//'TXWT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXWF(MLTM,NT) = OUTDIR//'TXWF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXWC(MLTM,NT) = OUTDIR//'TXWC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXWP(MLTM,NT) = OUTDIR//'TXWP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXBT(MLTM,NT) = OUTDIR//'TXBT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXBF(MLTM,NT) = OUTDIR//'TXBF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXBC(MLTM,NT) = OUTDIR//'TXBC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        FNTXBP(MLTM,NT) = OUTDIR//'TXBP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
        enddo
      endif
    endif
    if( MTMSRA(MLTM) == 1 )then
      FNAVV(MLTM) = OUTDIR//'AVVTS' // CNTMSR(MLTM) // '.OUT'
      FNAVB(MLTM) = OUTDIR//'AVBTS' // CNTMSR(MLTM) // '.OUT'
    endif
    if( MTMSRP(MLTM) == 1 )then
      FNSEL(MLTM) = OUTDIR//'SELTS' // CNTMSR(MLTM) // '.OUT'
    endif
    if( MTMSRUE(MLTM) == 1 )then
      FNUVE(MLTM) = OUTDIR//'UVETS' // CNTMSR(MLTM) // '.OUT'
    endif
    if( MTMSRUT(MLTM) == 1 )then
      FNUVT(MLTM) = OUTDIR//'UVTTS' // CNTMSR(MLTM) // '.OUT'
    endif
    if( MTMSRU(MLTM) == 1 )then
      FNU3D(MLTM) = OUTDIR//'U3DTS' // CNTMSR(MLTM) // '.OUT'
      FNV3D(MLTM) = OUTDIR//'V3DTS' // CNTMSR(MLTM) // '.OUT'
    endif
    if( MTMSRQE(MLTM) == 1 )then
      FNQQE(MLTM) = OUTDIR//'QQETS' // CNTMSR(MLTM) // '.OUT'
    endif
    if( MTMSRQ(MLTM) == 1 )then
      FNQ3D(MLTM) = OUTDIR//'Q3DTS' // CNTMSR(MLTM) // '.OUT'
    endif
  enddo
  JSTMSR = 0
  endif
  !
  if( JSTMSR == 0 ) GOTO 300
  !
  do MLTM = 1,MLTMSR
    if( MTMSRC(MLTM) == 1 )then
      if( ISTRAN(1) >= 1 )then
        FNSAL(MLTM) = OUTDIR//'SALTS' // CNTMSR(MLTM) // '.OUT'
        open(11,FILE = FNSAL(MLTM),STATUS = 'UNKNOWN')
        close(11,STATUS = 'DELETE')
        open(11,FILE = FNSAL(MLTM),STATUS = 'UNKNOWN')
        WRITE (11,100) TITLE1
        WRITE (11,101) CLTMSR(MLTM)
        WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (11,102) CTUNIT
        close(11)
      endif
      if( ISTRAN(2) >= 1 )then
        FNTEM(MLTM) = OUTDIR//'TEMTS' // CNTMSR(MLTM) // '.OUT'
        open(21,FILE = FNTEM(MLTM),STATUS = 'UNKNOWN')
        close(21,STATUS = 'DELETE')
        open(21,FILE = FNTEM(MLTM),STATUS = 'UNKNOWN')
        WRITE (21,100) TITLE2
        WRITE (21,101) CLTMSR(MLTM)
        WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (21,102) CTUNIT
        close(21)
      endif
      if( ISTRAN(3) >= 1 )then
        FNDYE(MLTM) = OUTDIR//'DYETS' // CNTMSR(MLTM) // '.OUT'
        open(31,FILE = FNDYE(MLTM),STATUS = 'UNKNOWN')
        close(31,STATUS = 'DELETE')
        open(31,FILE = FNDYE(MLTM),STATUS = 'UNKNOWN')
        WRITE (31,100) TITLE3
        WRITE (31,101) CLTMSR(MLTM)
        WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (31,102) CTUNIT
        close(31)
      endif
      if( ISTRAN(4) >= 1 )then
        FNDYE(MLTM) = OUTDIR//'SFLTS' // CNTMSR(MLTM) // '.OUT'
        open(31,FILE = FNSFL(MLTM),STATUS = 'UNKNOWN')
        close(31,STATUS = 'DELETE')
        open(31,FILE = FNSFL(MLTM),STATUS = 'UNKNOWN')
        WRITE (31,100) TITLE3
        WRITE (31,101) CLTMSR(MLTM)
        WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (31,102) CTUNIT
        close(31)
      endif
      if( ISTRAN(6) >= 1 )then
        FNSED(MLTM) = OUTDIR//'SEDTS' // CNTMSR(MLTM) // '.OUT'
        open(41,FILE = FNSED(MLTM),STATUS = 'UNKNOWN')
        close(41,STATUS = 'DELETE')
        open(41,FILE = FNSED(MLTM),STATUS = 'UNKNOWN')
        WRITE (41,100) TITLE4
        WRITE (41,101) CLTMSR(MLTM)
        WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (41,102) CTUNIT
        close(41)
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
        FNSND(MLTM,NX) = OUTDIR//'SND'// CNSND(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
        open(41,FILE = FNSND(MLTM,NX),STATUS = 'UNKNOWN')
        close(41,STATUS = 'DELETE')
        open(41,FILE = FNSND(MLTM,NX),STATUS = 'UNKNOWN')
        WRITE (41,100) TITLE4
        WRITE (41,101) CLTMSR(MLTM)
        WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (41,102) CTUNIT
        close(41)
        enddo
        do NX = 1,NSND
        FNSBL(MLTM,NX) = OUTDIR//'SBL'// CNSBL(NX) // 'TS' // CNTMSR(MLTM) // '.OUT'
        open(41,FILE = FNSBL(MLTM,NX),STATUS = 'UNKNOWN')
        close(41,STATUS = 'DELETE')
        open(41,FILE = FNSBL(MLTM,NX),STATUS = 'UNKNOWN')
        WRITE (41,100) TITLE18
        WRITE (41,101) CLTMSR(MLTM)
        WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (41,102) CTUNIT
        close(41)
        enddo
  !            FNSND(MLTM) = OUTDIR//'SNDTS' // CNTMSR(MLTM) // '.OUT'
  !            open(41,FILE = FNSND(MLTM),STATUS = 'UNKNOWN')
  !            close(41,STATUS = 'DELETE')
  !            open(41,FILE = FNSND(MLTM),STATUS = 'UNKNOWN')
  !            WRITE (41,100) TITLE4
  !            WRITE (41,101) CLTMSR(MLTM)
  !            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
  !            WRITE (41,102) CTUNIT
  !            close(41)
      endif
      if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
        FNBED(MLTM) = OUTDIR//'BEDTS' // CNTMSR(MLTM) // '.OUT'
        open(41,FILE = FNBED(MLTM),STATUS = 'UNKNOWN')
        close(41,STATUS = 'DELETE')
        open(41,FILE = FNBED(MLTM),STATUS = 'UNKNOWN')
        WRITE (41,100) TITLE19
        WRITE (41,101) CLTMSR(MLTM)
        WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (41,102) CTUNIT
        close(41)
      endif
      if( ISTRAN(8) >= 1 )then
        FNDOX(MLTM) = OUTDIR//'DOXTS' // CNTMSR(MLTM) // '.OUT'
        open(41,FILE = FNDOX(MLTM),STATUS = 'UNKNOWN')
        close(41,STATUS = 'DELETE')
        open(41,FILE = FNDOX(MLTM),STATUS = 'UNKNOWN')
        WRITE (41,100) TITLE4
        WRITE (41,101) CLTMSR(MLTM)
        WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (41,102) CTUNIT
        close(41)
        FNTOC(MLTM) = OUTDIR//'TOCTS' // CNTMSR(MLTM) // '.OUT'
        open(42,FILE = FNTOC(MLTM),STATUS = 'UNKNOWN')
        close(42,STATUS = 'DELETE')
        open(42,FILE = FNTOC(MLTM),STATUS = 'UNKNOWN')
        WRITE (42,100) TITLE4
        WRITE (42,101) CLTMSR(MLTM)
        WRITE (42,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (42,102) CTUNIT
        close(42)
        FNNHX(MLTM) = OUTDIR//'NHXTS' // CNTMSR(MLTM) // '.OUT'
        open(43,FILE = FNNHX(MLTM),STATUS = 'UNKNOWN')
        close(43,STATUS = 'DELETE')
        open(43,FILE = FNNHX(MLTM),STATUS = 'UNKNOWN')
        WRITE (43,100) TITLE4
        WRITE (43,101) CLTMSR(MLTM)
        WRITE (43,103)ILTMSR(MLTM),JLTMSR(MLTM)
        WRITE (43,102) CTUNIT
        close(43)
      endif
      if( ISTRAN(5) >= 1 )then
        do NT = 1,NTOX
          FNTOX(MLTM,NT) = OUTDIR//'TOX' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTOX(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTOX(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TITLE5
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXWT(MLTM,NT) = OUTDIR//'TXWT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXWT(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXWT(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXWT
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXWF(MLTM,NT) = OUTDIR//'TXWF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXWF(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXWF(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXWF
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXWC(MLTM,NT) = OUTDIR//'TXWC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXWC(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXWC(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXWC
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXWP(MLTM,NT) = OUTDIR//'TXWP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXWP(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXWP(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXWP
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXBT(MLTM,NT) = OUTDIR//'TXBT' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXBT(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXBT(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXBT
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXBF(MLTM,NT) = OUTDIR//'TXBF' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXBF(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXBF(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXBF
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXBC(MLTM,NT) = OUTDIR//'TXBC' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXBC(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXBC(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXBC
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
          FNTXBP(MLTM,NT) = OUTDIR//'TXBP' // CNTOX(NT) // 'TS' // CNTMSR(MLTM) // '.OUT'
          open(51,FILE = FNTXBP(MLTM,NT),STATUS = 'UNKNOWN')
          close(51,STATUS = 'DELETE')
          open(51,FILE = FNTXBP(MLTM,NT),STATUS = 'UNKNOWN')
          WRITE (51,100) TFNTXBP
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          close(51)
        enddo
      endif
    endif    ! ***  MTMSRC(MLTM) == 1

    if( MTMSRA(MLTM) == 1 )then
      FNAVV(MLTM) = OUTDIR//'AVVTS' // CNTMSR(MLTM) // '.OUT'
      open(61,FILE = FNAVV(MLTM),STATUS = 'UNKNOWN')
      close(61,STATUS = 'DELETE')
      open(61,FILE = FNAVV(MLTM),STATUS = 'UNKNOWN')
      WRITE (61,100) TITLE6
      WRITE (61,101) CLTMSR(MLTM)
      WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (61,102) CTUNIT
      close(61)
      
      FNAVB(MLTM) = OUTDIR//'AVBTS' // CNTMSR(MLTM) // '.OUT'
      open(71,FILE = FNAVB(MLTM),STATUS = 'UNKNOWN')
      close(71,STATUS = 'DELETE')
      open(71,FILE = FNAVB(MLTM),STATUS = 'UNKNOWN')
      WRITE (71,100) TITLE7
      WRITE (71,101) CLTMSR(MLTM)
      WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (71,102) CTUNIT
      close(71)
    endif
    
    if( MTMSRP(MLTM) == 1 )then
      FNSEL(MLTM) = OUTDIR//'SELTS' // CNTMSR(MLTM) // '.OUT'
      open(11,FILE = FNSEL(MLTM),STATUS = 'UNKNOWN')
      close(11,STATUS = 'DELETE')
      open(11,FILE = FNSEL(MLTM),STATUS = 'UNKNOWN')
      WRITE (11,100) TITLE11
      WRITE (11,101) CLTMSR(MLTM)
      WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (11,102) CTUNIT
      close(11)
    endif
    if( MTMSRUE(MLTM) == 1 )then
      FNUVE(MLTM) = OUTDIR//'UVETS' // CNTMSR(MLTM) // '.OUT'
      open(21,FILE = FNUVE(MLTM),STATUS = 'UNKNOWN')
      close(21,STATUS = 'DELETE')
      open(21,FILE = FNUVE(MLTM),STATUS = 'UNKNOWN')
      WRITE (21,100) TITLE12
      WRITE (21,101) CLTMSR(MLTM)
      WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (21,102) CTUNIT
      close(21)
    endif
    if( MTMSRUT(MLTM) == 1 )then
      FNUVT(MLTM) = OUTDIR//'UVTTS' // CNTMSR(MLTM) // '.OUT'
      open(31,FILE = FNUVT(MLTM),STATUS = 'UNKNOWN')
      close(31,STATUS = 'DELETE')
      open(31,FILE = FNUVT(MLTM),STATUS = 'UNKNOWN')
      WRITE (31,100) TITLE13
      WRITE (31,101) CLTMSR(MLTM)
      WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (31,102) CTUNIT
      close(31)
    endif
    if( MTMSRU(MLTM) == 1 )then
      FNU3D(MLTM) = OUTDIR//'U3DTS' // CNTMSR(MLTM) // '.OUT'
      open(41,FILE = FNU3D(MLTM),STATUS = 'UNKNOWN')
      close(41,STATUS = 'DELETE')
      open(41,FILE = FNU3D(MLTM),STATUS = 'UNKNOWN')
      WRITE (41,100) TITLE14
      WRITE (41,101) CLTMSR(MLTM)
      WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (41,102) CTUNIT
      close(41)
      FNV3D(MLTM) = OUTDIR//'V3DTS' // CNTMSR(MLTM) // '.OUT'
      open(51,FILE = FNV3D(MLTM),STATUS = 'UNKNOWN')
      close(51,STATUS = 'DELETE')
      open(51,FILE = FNV3D(MLTM),STATUS = 'UNKNOWN')
      WRITE (51,100) TITLE15
      WRITE (51,101) CLTMSR(MLTM)
      WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (51,102) CTUNIT
      close(51)
    endif
    if( MTMSRQE(MLTM) == 1 )then
      FNQQE(MLTM) = OUTDIR//'QQETS' // CNTMSR(MLTM) // '.OUT'
      open(61,FILE = FNQQE(MLTM),STATUS = 'UNKNOWN')
      close(61,STATUS = 'DELETE')
      open(61,FILE = FNQQE(MLTM),STATUS = 'UNKNOWN')
      WRITE (61,100) TITLE16
      WRITE (61,101) CLTMSR(MLTM)
      WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (61,102) CTUNIT
      close(61)
    endif
    if( MTMSRQ(MLTM) == 1 )then
      FNQ3D(MLTM) = OUTDIR//'Q3DTS' // CNTMSR(MLTM) // '.OUT'
      open(71,FILE = FNQ3D(MLTM),STATUS = 'UNKNOWN')
      close(71,STATUS = 'DELETE')
      open(71,FILE = FNQ3D(MLTM),STATUS = 'UNKNOWN')
      WRITE (71,100) TITLE17
      WRITE (71,101) CLTMSR(MLTM)
      WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
      WRITE (71,102) CTUNIT
      close(71)
    endif
  enddo
  !
  JSTMSR = 0
  !
  !----------------------------------------------------------------------C
  !
    300 continue
  !
  !----------------------------------------------------------------------C
  !
  TIME = TIMESEC/TCTMSR

  !
  FOURDPI = 4./PI
  !
  ! ***  STEP CURRENT TIME INTERVALS FOR WRITE SCENARIOS
  !
  do NTSSS = 1,NTSSTSP
   do MTSSS = 1,MTSSTSP(NTSSS)
    if( TIME >= TSSTRT(MTSSS,NTSSS) )then
    if( TIME <= TSSTOP(MTSSS,NTSSS) )then
      MTSCUR(NTSSS) = MTSSS
    endif
    endif
   enddo
  enddo
  !
  if( MDCHH >= 1 )then
    do L = 2,LA
      QCHANUIJ(L) = 0.0
      QCHANVIJ(L) = 0.0
    enddo
    do NMD = 1,MDCHH
      LMDCHHT = LMDCHH(NMD)
      LMDCHUT = LMDCHU(NMD)
      LMDCHVT = LMDCHV(NMD)
      if( MDCHTYP(NMD) == 1 )then
        QCHANUIJ(LMDCHHT) = QCHANUIJ(LMDCHHT)+QCHANU(NMD)
        QCHANUIJ(LMDCHUT) = QCHANUIJ(LMDCHUT)-QCHANU(NMD)
      endif
      if( MDCHTYP(NMD) == 2 )then
        QCHANVIJ(LMDCHHT) = QCHANVIJ(LMDCHHT)+QCHANV(NMD)
        QCHANVIJ(LMDCHVT) = QCHANVIJ(LMDCHVT)-QCHANV(NMD)
      endif
      if( MDCHTYP(NMD) == 3 )then
        QCHANUIJ(LMDCHHT) = QCHANUIJ(LMDCHHT)+QCHANU(NMD)
        QCHANUIJ(LMDCHUT) = QCHANUIJ(LMDCHUT)-QCHANU(NMD)
        QCHANVIJ(LMDCHHT) = QCHANVIJ(LMDCHHT)+QCHANV(NMD)
        QCHANVIJ(LMDCHVT) = QCHANVIJ(LMDCHVT)-QCHANV(NMD)
      endif
    enddo
  endif
  !
  do MLTM = 1,MLTMSR
   NTSSS = NTSSSS(MLTM)
   MTSCC = MTSCUR(NTSSS)
   if( TIME >= TSSTRT(MTSCC,NTSSS) )then
   if( TIME <= TSSTOP(MTSCC,NTSSS) )then
    I = ILTMSR(MLTM)
    J = JLTMSR(MLTM)
    L = LIJ(I,J)
    LN = LNC(L)
    if( MTMSRC(MLTM) == 1 )then
      if( ISTRAN(1) >= 1 )then
        open(11,FILE = FNSAL(MLTM),POSITION = 'APPEND')
        write(11,201)TIME,(SAL(L,K),K = 1,KC)
        close(11)
      endif
      if( ISTRAN(2) >= 1 )then
        open(21,FILE = FNTEM(MLTM),POSITION = 'APPEND')
        write(21,201)TIME,(TEM(L,K),K = 1,KC)
        close(21)
      endif
      if( ISTRAN(3) >= 1 )then
        MD = 1  ! DELME - HARDWIRED
        open(31,FILE = FNDYE(MLTM),POSITION = 'APPEND')
        write(31,201)TIME,(DYE(L,K,MD),K = 1,KC)
        close(31)
      endif
      if( ISTRAN(4) >= 1 )then
        open(31,FILE = FNSFL(MLTM),POSITION = 'APPEND')
        write(31,201)TIME,(SFL(L,K),K = 1,KC)
        close(31)
      endif
      if( ISTRAN(6) >= 1 )then
        open(41,FILE = FNSED(MLTM),POSITION = 'APPEND')
        if( ISNDAL == 2 )then
          SEDBTMP = SEDBT(L,KBT(L))+SEDBT(L,KBT(L)-1)
        else
          SEDBTMP = SEDBT(L,KBT(L))
        endif
        write(41,201)TIME,SEDBTMP,(SEDT(L,K),K = 1,KC)
        close(41)
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          open(41,FILE = FNSND(MLTM,NX),POSITION = 'APPEND')
          if( ISNDAL == 2 )then
            SNDBTMP = SNDB(L,KBT(L),NX)+SNDB(L,KBT(L)-1,NX)
          else
            SNDBTMP = SNDB(L,KBT(L),NX)
          endif
          write(41,201)TIME,SNDBTMP,(SND(L,K,NX),K = 1,KC), &
                                   SNDEQSAV(L,NX)
          close(41)
        enddo
        do NX = 1,NSND
          open(41,FILE = FNSBL(MLTM,NX),POSITION = 'APPEND')
  !            CQBEDLOADX = 0.
  !            CQBEDLOADY = 0.
          if( UHDYE(L) /= 0.0 ) CQBEDLOADX(L,NX) = QSBDLDX(L,NX)/UHDYE(L)
          if( VHDXE(L) /= 0.0 ) CQBEDLOADY(L,NX) = QSBDLDY(L,NX)/VHDXE(L)
          write(41,201)TIME,QSBDLDX(L,NX),QSBDLDY(L,NX), &
                 CQBEDLOADX(L,NX),CQBEDLOADY(L,NX),SNDFBL(L,NX)
          close(41)
        enddo
  !            open(41,FILE = FNSND(MLTM,NX),POSITION = 'APPEND')
  !            write(41,201)TIME,SNDBT(L,KBT(L)),(SNDT(L,K),K = 1,KC)
  !            close(41)
      endif
      if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
        open(41,FILE = FNBED(MLTM),POSITION = 'APPEND')
        KTMP = KBT(L)
        KTMP1 = KBT(L)-1
        KTMP1 = max(KTMP1,1)
        NSXD = NSED+NSND
  !            write(41,221)TIME,KTMP,HBED(L,KTMP),HBED(L,KTMP1),
        write(41,221)TIME,KTMP,HBED(L,KTMP),HBED(L,KTMP1),HBEDA(L), &
                VDRBED(L,KTMP),(VFRBED(L,KTMP,NX),NX = 1,NSXD)
        close(41)
      endif
      if( ISTRAN(8) >= 1 )then
        open(41,FILE = FNDOX(MLTM),POSITION = 'APPEND')
        write(41,201)TIME,(WQV(L,K,19),K = 1,KC)
        close(41)
        open(42,FILE = FNTOC(MLTM),POSITION = 'APPEND')
        write(42,201)TIME,(WQV(L,K,6),K = 1,KC)
        close(42)
        open(43,FILE = FNNHX(MLTM),POSITION = 'APPEND')
        write(43,201)TIME,(WQV(L,K,14),K = 1,KC)
        close(43)
      endif
      if( ISTRAN(5) >= 1 )then
        do NT = 1,NTOX
  !
          open(51,FILE = FNTOX(MLTM,NT),POSITION = 'APPEND')
          open(52,FILE = FNTXWT(MLTM,NT),POSITION = 'APPEND')
          open(53,FILE = FNTXWF(MLTM,NT),POSITION = 'APPEND')
          open(54,FILE = FNTXWC(MLTM,NT),POSITION = 'APPEND')
          open(55,FILE = FNTXWP(MLTM,NT),POSITION = 'APPEND')
          open(56,FILE = FNTXBT(MLTM,NT),POSITION = 'APPEND')
          open(57,FILE = FNTXBF(MLTM,NT),POSITION = 'APPEND')
          open(58,FILE = FNTXBC(MLTM,NT),POSITION = 'APPEND')
          open(59,FILE = FNTXBP(MLTM,NT),POSITION = 'APPEND')
  !
          NDOC = NSED+NSND+1
          do K = 1,KC
            SEDSND(K) = SEDT(L,K)+SNDT(L,K)
            TXWF(K) = TOXFDFW(L,K,NT)*TOX(L,K,NT)
            TXWC(K) = TOXCDFW(L,K,NT)*TOX(L,K,NT)
            TXWP(K) = TOXPFTW(L,K,NT)*TOX(L,K,NT)
          enddo
          do K = 1,KB
            SEDSNDB(K) = SEDBT(L,K)+SNDBT(L,K)
            TXBF(K) = 0.
            TXBC(K) = 0.
            TXBP(K) = 0.
            TXBT(K) = 0.
          enddo
          do K = 1,KBT(L)
            PORH(K) = 1.0/PORBED(L,K)
            TXBF(K) = TOXFDFB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
            TXBC(K) = TOXCDFB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
            TXBP(K) = TOXPFTB(L,K,NT)*TOXB(L,K,NT)/HBED(L,K)
            TXBT(K) = TOXB(L,K,NT)/HBED(L,K)
          enddo

          K = KBT(L)
          write(51,201)TIME,TXBT(K),TXBF(K),TXBC(K),TXBP(K),TOX(L,KSZ(L),NT),TXWF(1),TXWC(1),TXWP(1)
          do K = 1,KC
            if( SEDSND(K) > 0.0)TXWP(K) = 1000.*TXWP(K)/SEDSND(K)
          enddo

          do K = 1,KBT(L)
            TXBP(K) = TXBP(K)*HBED(L,K)
            TXBF(K) = TXBF(K)*PORH(K)
            TXBC(K) = TXBC(K)*PORH(K)
            if( SEDSNDB(K) > 0.0)TXBP(K) = 1000.*TXBP(K)/SEDSNDB(K)
          enddo
          write(52,201)TIME,(TOX(L,K,NT),K = 1,KC)
          write(53,201)TIME,(TXWF(K),K = 1,KC)
          write(54,201)TIME,(TXWC(K),K = 1,KC)
          write(55,201)TIME,(TXWP(K),K = 1,KC)
          write(56,201)TIME,(TXBT(K),K = 1,KB)
          write(57,201)TIME,(TXBF(K),K = 1,KB)
          write(58,201)TIME,(TXBC(K),K = 1,KB)
          write(59,201)TIME,(TXBP(K),K = 1,KB)
  !
          close(51)
          close(52)
          close(53)
          close(54)
          close(55)
          close(56)
          close(57)
          close(58)
          close(59)
  !
        enddo
  !            do NT = 1,NTOX
  !            open(51,FILE = FNTOXPF(MLTM,NT),POSITION = 'APPEND')
  !            write(51,201)TIME,TOXPFTB(L,KBT(L),NT),
  !     &                  (TOXPFTW(L,K,NT),K = 1,KC)
  !            close(51)
  !            enddo
  !            do NT = 1,NTOX
  !            open(51,FILE = FNTOXFD(MLTM,NT),POSITION = 'APPEND')
  !            write(51,201)TIME,TOXPFTB(L,KBT(L),NT),
  !     &                  (TOXPFTW(L,K,NT),K = 1,KC)
  !            close(51)
  !            enddo
      endif
    endif    ! ***  MTMSRC(MLTM) == 1
    
    if( MTMSRA(MLTM) == 1 )then
      open(61,FILE = FNAVV(MLTM),POSITION = 'APPEND')
      open(71,FILE = FNAVB(MLTM),POSITION = 'APPEND')
      do K = 1,KS
        ATMPP(K) = 10000.*AV(L,K)*HP(L)
      enddo
      write(61,201)TIME,(ATMPP(K),K = 1,KS)
      do K = 1,KS
        ATMPP(K) = 10000.*AB(L,K)*HP(L)
      enddo
      write(71,201)TIME,(ATMPP(K),K = 1,KS)
      close(61)
      close(71)
    endif
    if( MTMSRP(MLTM) == 1 )then
      open(11,FILE = FNSEL(MLTM),POSITION = 'APPEND')
      PPTMP = HP(L)+BELV(L)
      TMPVAL = VDWASTE(L)/DXYP(L)
      if( (ISTRAN(6)+ISTRAN(7)) > 0 )then
        write(11,201)TIME,PPTMP,HP(L),BELV(L),HBEDA(L),ZELBEDA(L),TMPVAL,VDWASTE(L)
      else
        write(11,201)TIME,PPTMP,HP(L),BELV(L),0.,0.,TMPVAL,VDWASTE(L)
      endif
      close(11)
    endif
    if( MTMSRUE(MLTM) == 1 )then
      open(21,FILE = FNUVE(MLTM),POSITION = 'APPEND')
      UTMP1 = 50.*(UHDYE(LEC(L))+UHDYE(L))/(DYP(L)*HP(L))
      VTMP1 = 50.*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
      if( SPB(L) == 0 )then
        UTMP1 = 2.*UTMP1
        VTMP1 = 2.*VTMP1
      endif
      UTMP = CUE(L)*UTMP1+CVE(L)*VTMP1
      VTMP = CUN(L)*UTMP1+CVN(L)*VTMP1
      UTMP1 = 5000.*(TBX(LEC(L))+TBX(L))
      VTMP1 = 5000.*(TBY(LN)+TBY(L))
      TBEAST = CUE(L)*UTMP1+CVE(L)*VTMP1
      TBNORT = CUN(L)*UTMP1+CVN(L)*VTMP1
      TAUBDYN = 10000.*TAUB(L)
      TAUB2 = TAUBDYN*TAUBDYN
      if( ISWAVE > 0 )then
        CURANG = ATAN2(VTMP,UTMP)
        TAUW1 = 10000.*QQWV1(L)
        TAUW2 = 10000.*QQWV2(L)
        TAUB2 = TAUB2 + (TAUW2*TAUW2) + 2.*TAUBDYN*TAUW2*COS(CURANG-WV(L).DIR)
      endif
      TAUB2 = max(TAUB2,0.)
      TAUTOT = SQRT(TAUB2)
      VELMAG = UTMP*UTMP+VTMP*VTMP
      TMPDRAG = 0.0
      
      if( (ISTRAN(6)+ISTRAN(7)) > 0 )then
        TAUBSEDDYN = 10000.*TAUBSED(L)
        TAUBSNDDYN = 10000.*TAUBSND(L)
        if( VELMAG > 0.0)TMPDRAG = TAUTOT/VELMAG
        write(21,201)TIME,UTMP,VTMP,TBEAST,TBNORT,TAUBDYN,TAUBSEDDYN,TAUBSNDDYN
      else
        write(21,201)TIME,UTMP,VTMP,TBEAST,TBNORT,TAUBDYN,0.,0.
      endif
      close(21)
    endif
    if( MTMSRUT(MLTM) == 1 )then
      open(31,FILE = FNUVT(MLTM),POSITION = 'APPEND')
      QBEDLOADX = 0.
      QBEDLOADY = 0.
      CQBEDLOADXT = 0.
      CQBEDLOADYT = 0.
      do NX = 1,NSND
        QBEDLOADX = QBEDLOADX+QSBDLDX(L,NX)
        QBEDLOADY = QBEDLOADY+QSBDLDY(L,NX)
        CQBEDLOADXT = CQBEDLOADXT+CQBEDLOADX(L,NX)
        CQBEDLOADYT = CQBEDLOADYT+CQBEDLOADY(L,NX)
      enddo
  !          if( UHDYE(L) /= 0.0)CQBEDLOADX = QBEDLOADX/UHDYE(L)
  !          if( VHDXE(L) /= 0.0)CQBEDLOADY = QBEDLOADY/VHDXE(L)
      write(31,201)TIME,UHDYE(L),VHDXE(L),QBEDLOADX,QBEDLOADY, &
                   CQBEDLOADXT,CQBEDLOADYT
      close(31)
    endif
    if( MTMSRU(MLTM) == 1 )then
      open(41,FILE = FNU3D(MLTM),POSITION = 'APPEND')
      open(51,FILE = FNV3D(MLTM),POSITION = 'APPEND')
      RUVTMP = 50.
      if( SPB(L) == 0 ) RUVTMP = 100.
      do K = 1,KC
       UTMP1 = RUVTMP*(U(LEC(L),K)+U(L,K))
       VTMP1 = RUVTMP*(V(LN,K)+V(L,K))
       ATMPP(K) = CUE(L)*UTMP1+CVE(L)*VTMP1
       BTMP(K) = CUN(L)*UTMP1+CVN(L)*VTMP1
      enddo
      write(41,201)TIME,(ATMPP(K),K = 1,KC)
      write(51,201)TIME,(BTMP(K),K = 1,KC)
      close(41)
      close(51)
    endif
    if( MTMSRQE(MLTM) == 1 )then
      open(61,FILE = FNQQE(MLTM),POSITION = 'APPEND')
      QRNT = DXYP(L)*RAINT(L)
      write(61,201)TIME,QSUME(L),QRNT,EVAPSW(L),EVAPGW(L),QGW(L),QCHANUIJ(L),QCHANVIJ(L),QDWASTE(L),VDWASTE(L)
      close(61)
    endif
    if( MTMSRQ(MLTM) == 1 )then
      open(71,FILE = FNQ3D(MLTM),POSITION = 'APPEND')
      WRITE (71,201)TIME,(QSUM(L,K),K = 1,KC)
      close(71)
    endif
   endif
   endif
  enddo
  !
  ! *** *******************************************************************C
  !
    100 FORMAT(A80)
    101 FORMAT('  AT LOCATION  ',A20)
    102 FORMAT('  TIME IN FIRST COLUMN HAS UNITS OF ',A10)
    103 FORMAT('  CELL I,J = ',2I5)
    201 FORMAT(F12.5,12E12.4)
    221 FORMAT(F12.5,I5,14E12.4)
  !
  ! *** *******************************************************************C
  !
  return
END
