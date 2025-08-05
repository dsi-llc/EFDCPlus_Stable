  ! ----------------------------------------------------------------------
  !   This file is a part of EFDC+
  !   Website:  https://eemodelingsystem.com/
  !   Repository: https://github.com/dsi-llc/EFDC_Plus.git
  ! ----------------------------------------------------------------------
  ! Copyright 2021-2024 DSI, LLC
  ! Distributed under the GNU GPLv2 License.
  ! ----------------------------------------------------------------------
  SUBROUTINE INPUT()

  ! *** SUBROUTINE INPUT READS ALL INPUT DATA EXCEPT DATA IN LXLY.INP,
  ! *** MASK.INP AND RESTART.INP

  !----------------------------------------------------------------------!
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3

  use GLOBAL
  use Variables_WQ
  use RESTART_MODULE, only:Setup_Continuation_Files
  use INFOMOD ,only:SKIPCOM, READSTR
  use HYDSTRUCMOD
  use SHELLFISHMOD, only: NSF, ISFFARM, NSFCELLS, READ_SHELLFISH_JSON
  use Variables_Propwash
  use FIELDS
  use CYCLONE
  use Allocate_Initialize

  ! *** New for MPI
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out
  use Broadcast_Routines
  use Mod_GOTM
  ! *** End new for MPI

#ifndef GNU  
  use IFPORT
#endif

  implicit none

  character*120 TEXT
  character*200 STR
  character*10  CDUM
  character*3   NCARD
  character     ADUMMY*5,RESTARTF*50,STRC*650

  real(RK4) :: SEEPRATE(1000)
  real,allocatable,dimension(:) :: RMULADS
  real,allocatable,dimension(:) :: ADDADS
  real,allocatable,dimension(:) :: BOFFMHK, BOFFSUP, TOFFMHK, TOFFSUP
  real,allocatable,dimension(:,:) :: PFX2
  real,allocatable,dimension(:,:) :: QSERSM,RULES

  integer,allocatable,dimension(:) :: IPARTSP   ! *** NEW BEDFORD
  integer,allocatable,dimension(:) :: IDX
  integer,allocatable,dimension(:,:) :: IBLTMP

  integer :: ITURB = 0, IPMC, IS, NS, NX, NT, NW, NA, NI, NWR
  integer :: L, LP, LW, LE, LS, LN, LG, LU, LD
  integer :: IDUM, NDUM, LCM2T, I, J, K, M, KDUM, IMDXDY, NPFORN
  integer :: ISO, ITIDASM, NPFOR, NPFORS, NPFORW, NPFORE, ITSSS
  integer :: NDUM1, NDUM2, NTMP, ITYPE, IFLAG
  integer :: MTMP
  integer :: MMAX, MS, MMIN, NJP
  integer :: JSFDCH, IDUMMY, NPP, MD, MU
  integer :: MTSSS, IACROSS, JCTMP, JACROSS, JT, JF, JLAST, NMD
  integer :: IFIRST, ILAST, IT, NP, LT, MVEGIJT, NMDXDY, INITTEMP
  integer :: ITMP, JTMP, LTMP, LL, ITMPU, JTMPU, ITMPD, JTMPD, ID, JD, NSEEPCLASSES
  integer :: IZONE, LDUM, JDUM, IVAL, ISALTYP, IREAD, KBINPUT, ISCOLD
  integer :: ISTYP, ICHGQS, NCTMP
  integer :: NFLAGPWR, ITMPPMX, NPMXZ, NPMXPTS, PMIXSF, NZ
  integer :: NPBPH, IASERVER, NC, ISWDINT
  integer :: ISQCTRL(6)
  integer :: IBEDBDNU
  integer :: IBEDDDNU
  integer :: IBEDLAYU

  real*8 :: TMPDATE, TIMEQCTLSERIES
  real   :: ADMAX, ADMIN, AHMAX, AHMIN
  real   :: DXIJ, DYIJ, HIJ, BELVIJ, ZBRIJ, RVALUE, PSERTMP, DIAMMHK
  real   :: DXYCVT, HADADJ, RAD, AMP, T1, T2, TMPAMP, TMPPHS, BOTTOM
  real   :: TOP1, TOP2
  real   :: SEDVDRT, DSTR, USTR, DUM, POR, QSERTMP
  real   :: RMDX, RMDY, CVTFACX, CVTFACY, FBODY1, FBODY2, BDLTMP
  real   :: RMULADJ, ADDADJ, RMULADJS, ADDADJS, PSERTMPS, CSERTMP
  real   :: QCTLTMP, SOLRCVT, CLDCVT, WINDSCT, RMULADJCOV, RMULADJTHK
  real   :: TDELTA, TMPCVT, WSEL, TMP, TOFFSET

  real,external :: SETSTVEL, PARSE_REAL
  logical :: BFLAG
  logical(4) :: RES

  integer :: RI_GLOBAL, RJ_GLOBAL
  integer :: IGL, JGL
  integer :: i_loc, j_loc, l_local, nn, ierr
  integer, allocatable, dimension(:) :: LPMXZ_Global

  real,allocatable,dimension(:)   :: CONINIT
  real,allocatable,dimension(:)   :: CONBINIT
  real,allocatable,dimension(:)   :: CQSTMP
  real,allocatable,dimension(:,:) :: CPFAM0        !< Tides
  real,allocatable,dimension(:,:) :: CPFAM1        !< Tides
  real,allocatable,dimension(:,:) :: CPFAM2        !< Tides
  real,allocatable,dimension(:,:) :: SPFAM0
  real,allocatable,dimension(:,:) :: SPFAM1
  real,allocatable,dimension(:,:) :: SPFAM2

  call AllocateDSI( RMULADS,  NSTM2,  0.0)
  call AllocateDSI( ADDADS,   NSTM2,  0.0)
  call AllocateDSI( IPARTSP,  NTXM,     0)
  call AllocateDSI( CONINIT,  KCM,    0.0)
  call AllocateDSI( CONBINIT, KBM,    0.0)
  call AllocateDSI( CQSTMP,   NSTVM2, 0.0)
  call AllocateDSI( QSERSM,   NDQSER, KCM, 0.0)

  G = 9.81
  PI = 3.1415926535898
  PI2 = 2.*PI
2 FORMAT(A80)

  NCARD = '1'
  if( process_id == master_id )then
    ! *** READ MAIN INPUT FILE EFDC.INP
    write(*,'(A)')'READING THE MAIN EFDC CONTROL FILE: EFDC.INP'
    open(1,FILE = 'efdc.inp',STATUS = 'UNKNOWN')

    !1**  READ TITLE CARD
    call SEEK('C1',0)
    read(1,2) TITLE
  endif
  call Broadcast_Scalar(TITLE,     master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,2) TITLE

  !C1A**  READ MODE OPTIONS
  NCARD = '1A'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C1A',0)
    read(1,*,IOSTAT = ISO) IS2TIM, IGRIDH, IGRIDV, KMINV, SGZHPDELTA  !, ISWGS84        ! NTL: Waiting for the implementation of Geographic Coordinate
  endif

  call Broadcast_Scalar(ISO,        master_id)
  call Broadcast_Scalar(IS2TIM,     master_id)
  call Broadcast_Scalar(IGRIDH,     master_id)
  call Broadcast_Scalar(IGRIDV,     master_id)
  call Broadcast_Scalar(KMINV,      master_id)
  call Broadcast_Scalar(SGZHPDELTA, master_id)
  !Call Broadcast_Scalar(ISWGS84, master_id)        ! NTL: Waiting for the implementation of Geographic Coordinate

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) IS2TIM,IGRIDH,IGRIDV,KMINV,SGZHPDELTA  !, ISWGS84        ! NTL: Waiting for the implementation of Geographic Coordinate
  if( ISO > 0 ) GOTO 100

  !C2**  READ RESTART AND DIAGNOSTIC SWITCHES
  ! *** ********************************************************
  NCARD = '2'
  if( process_id == master_id )then
    call SEEK('C2',0)
    read(1,*,IOSTAT = ISO) ISRESTI, ISRESTO, ISRESTR, ISGREGOR, ISLOG, ISDIVEX, ISNEGH, ISMMC, ldum, ICONTINUE, ISHOW
  endif

  call Broadcast_Scalar(ISO,      master_id)
  call Broadcast_Scalar(ISRESTI,  master_id)
  call Broadcast_Scalar(ISRESTO,  master_id)
  call Broadcast_Scalar(ISRESTR,  master_id)
  call Broadcast_Scalar(ISGREGOR, master_id)
  call Broadcast_Scalar(ISLOG,    master_id)
  call Broadcast_Scalar(ISDIVEX,  master_id)
  call Broadcast_Scalar(ISNEGH,   master_id)
  call Broadcast_Scalar(ISMMC,    master_id)
  call Broadcast_Scalar(ICONTINUE,master_id)
  call Broadcast_Scalar(ISHOW,    master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) ISRESTI, ISRESTO, ISRESTR, ISGREGOR, ISLOG, ISDIVEX, ISNEGH, ISMMC, ldum, ICONTINUE, ISHOW
  if( ISO > 0 ) GOTO 100

  ! *** HANDLE BATHYMETRY ADJUSTMENTS
  ISRESTIOPT = 0
  if( ISRESTI == -1 )then
    ISRESTIOPT = 1
    ISRESTI = 1
  endif

  ! *** ********************************************************
  RESTARTF = " "          ! *** Initialize
  Restart_In_Ver = -1     ! *** Initialize

  if( ISRESTI == 1 .and. ICONTINUE == 1 )then
    ! *** FOR CONTINUATION MODE:
    if( process_id == master_id )then
      call SEEK('C2A',0)
      read(1,*,IOSTAT = ISO) Restart_In_Ver
      read(1,*,IOSTAT = ISO) RESTARTF
    endif
  endif
  call Broadcast_Scalar(ISO,            master_id)
  call Broadcast_Scalar(Restart_In_Ver, master_id)
  call Broadcast_Scalar(RESTARTF,       master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) Restart_In_Ver, RESTARTF
  if( ISO > 0 ) GOTO 100

  if( ISMMC < 0 )then
    DEBUG = .TRUE.
    ISMMC = 0
    if( process_id == master_id ) WRITE(*,'(A)') 'DEBUG ON'
  else
    DEBUG = .FALSE.
    if( process_id == master_id ) WRITE(*,'(A)') 'DEBUG OFF'
  endif

  !C3**  READ RELAXATION parameterS
  NCARD = '3'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C3',0)
    read(1,*,IOSTAT = ISO) RP, RSQM, ITERM, IRVEC, IATMP, IWDRAG, ITERHPM, ldum, ISDSOLV, tmp
  endif

  call Broadcast_Scalar(RP     , master_id)
  call Broadcast_Scalar(RSQM   , master_id)
  call Broadcast_Scalar(ITERM  , master_id)
  call Broadcast_Scalar(IRVEC  , master_id)
  call Broadcast_Scalar(IATMP  , master_id)
  call Broadcast_Scalar(IWDRAG , master_id)
  call Broadcast_Scalar(ITERHPM, master_id)
  call Broadcast_Scalar(ISDSOLV, master_id)

  if( IRVEC /= 0 .and. IRVEC /= 9 ) CALL STOPP('INVALID IRVEC')
  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) RP, RSQM, ITERM, IRVEC, IATMP, IWDRAG, ITERHPM, ldum, ISDSOLV, tmp
  if( ISO > 0 ) GOTO 100

  !C4**  READ LONGTERM MASS TRANSPORT INTEGRATION ONLY SWITCHES
  NCARD = '4'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C4',0)
    read(1,*,IOSTAT = ISO) ISLTMT, ISSSMMT, RESSTEP
  endif

  call Broadcast_Scalar(ISO ,    master_id)
  call Broadcast_Scalar(ISLTMT , master_id)
  call Broadcast_Scalar(ISSSMMT, master_id)
  call Broadcast_Scalar(RESSTEP, master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) ISLTMT, ISSSMMT, RESSTEP
  if( ISO > 0 ) GOTO 100

  !C5**  READ MOMENTUM ADVECTION AND DIFFUSION SWITCHES AND MISC
  NCARD = '5'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C5',0)
    read(1,*,IOSTAT = ISO) ISCDMA, ISHDMF, ldum, ISWASP, ISDRY, ICALTB, ISRLID, ISVEG, ISVEGL, ISITB, IHMDSUB, IINTPG, ISHDMFILTER

    if( MACDRAG > 0 .and. ISVEG == 0 )  ISVEG = 2
  endif

  call Broadcast_Scalar(ISCDMA ,     master_id)
  call Broadcast_Scalar(ISHDMF ,     master_id)
  call Broadcast_Scalar(ISWASP ,     master_id)
  call Broadcast_Scalar(ISDRY  ,     master_id)
  call Broadcast_Scalar(ICALTB ,     master_id)
  call Broadcast_Scalar(ISRLID ,     master_id)
  call Broadcast_Scalar(ISVEG  ,     master_id)
  call Broadcast_Scalar(ISVEGL ,     master_id)
  call Broadcast_Scalar(ISITB  ,     master_id)
  call Broadcast_Scalar(IHMDSUB,     master_id)
  call Broadcast_Scalar(IINTPG ,     master_id)
  call Broadcast_Scalar(ISHDMFILTER, master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) ISCDMA, ISHDMF, ldum, ISWASP, ISDRY, ICALTB, ISRLID, ISVEG, ISVEGL, ISITB, IHMDSUB, IINTPG, ISHDMFILTER
  if( ISO > 0 ) GOTO 100

  IDRYTBP = 0
  if( ISDRY < 0 )then
    ISDRY = ABS(ISDRY)
    IDRYTBP = 1
  endif
  if( ISVEG < 1 )ISITB = 0
  if( ISWASP == 99)ISICM = 1
  if( ISRLID == 1) ISDRY = -1
  if( ISWASP == 10)ISRCA = 1
  JSWAVE = 0
  ! PMC      IS1DCHAN = 0
  ! PMC  if( ISCDMA == 10) IS1DCHAN = 1                              Y

  !C6**  DISSOLVED AND SUSPENDED CONSTITUENT TRANSPORT SWITCHES
  NCARD = '6'

  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C6',0)
    do NS = 0,8
      read(1,*,IOSTAT = ISO) ISTRAN(NS), ISTOPT(NS), ldum, ISADAC(NS), ISFCT(NS), ldum, ldum, ldum, ISCI(NS), ISCO(NS)
      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) ISTRAN(NS), ISTOPT(NS), ldum, ISADAC(NS), ISFCT(NS), ldum, ldum, ldum, ISCI(NS), ISCO(NS)
      if( ISO > 0 ) GOTO 100
    enddo
  endif

  call Broadcast_Array(ISTRAN, master_id)
  call Broadcast_Array(ISTOPT, master_id)
  call Broadcast_Array(ISADAC, master_id)
  call Broadcast_Array(ISFCT , master_id)
  call Broadcast_Array(ISCI  , master_id)
  call Broadcast_Array(ISCO  , master_id)
  NFLTMT = 1

  if( ISTRAN(8) >= 1 .and. ISTRAN(2) == 0 )then
    PRINT *,'*** WARNING: TEMPERATURE SHOULD BE ACTIVATED FOR WQ CALCULATIONS!'
#ifdef GNU  
    call SLEEP(5)       !< 5 seconds
#else
    call SLEEPQQ(5000)  !< 5000 milliseconds
#endif 
  endif

  ! *** ********************************************************
  if( ISRESTI == 1 .and. ICONTINUE == 1 )then
    ! *** MODEL RUN CONTINUATION
    if( process_id == master_id )then
      close(1)
      call Setup_Continuation_Files(RESTARTF)           ! GENERATE RESTART FILES, OTHERWISE EFDC LOADS EXISTING ONES
      open(1,FILE = 'efdc.inp',STATUS = 'UNKNOWN')
    endif
    call Broadcast_Scalar(NRESTART,         master_id)
    call Broadcast_Scalar(TBEGINC,          master_id)
    call Broadcast_Scalar(Restart_In_Ver,   master_id)
  endif

  !  *** DEACTIVATE ANY UNUSED OPTIONS
  if( ISTRAN(2) < 1 ) ISTOPT(2) = 0

  ! *** SET TRANSPORT FLAG
  ISTRANACTIVE = 0
  do I = 1,8
    if( ISTRAN(I) > 0 ) ISTRANACTIVE = 1
  enddo

  !C7**  READ TIME-RELATED INTEGER parameterS
  NCARD = '7'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C7',0)
    read(1,*,IOSTAT = ISO) NTC, NTSPTC, NLTC, NTTC, ldum, NTSTBC, ldum, NTCVB, NTSMMT, ldum, NDRYSTP, NRAMPUP, NUPSTEP

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) NTC, NTSPTC, NLTC, NTTC, ldum, NTSTBC, ldum, NTCVB, NTSMMT, ldum, NDRYSTP, NRAMPUP, NUPSTEP
    ! *** Not Used:  NTCPP
    ! *** Not Used:  NTCNB
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(NTC    , master_id)
  call Broadcast_Scalar(NTSPTC , master_id)
  call Broadcast_Scalar(NLTC   , master_id)
  call Broadcast_Scalar(NTTC   , master_id)
  call Broadcast_Scalar(NTSTBC , master_id)
  call Broadcast_Scalar(NTCVB  , master_id)
  call Broadcast_Scalar(NTSMMT , master_id)
  call Broadcast_Scalar(NDRYSTP, master_id)
  call Broadcast_Scalar(NRAMPUP, master_id)
  call Broadcast_Scalar(NUPSTEP, master_id)

  if( NRAMPUP < 1 ) NRAMPUP = 1
  if( NUPSTEP < 2 ) NUPSTEP = 2

  NDRYSTP = ABS(NDRYSTP)   ! *** IF > 0, EFDC+ WILL WASTE ISOLATED CELL WATER AFTER NDRYSTP STEPS

  !C8**  READ TIME-RELATED REAL parameterS
  NCARD = '8'
  ! *** *******************************************************
  if( process_id == master_id )then
    call SEEK('C8',0)
    read(1,*,IOSTAT = ISO) TCON, TBEGIN, TIDALP, CF, ISCORV, ISDCCA, ISCFL, ISCFLM, DTSSFAC, DTSSDHDT, DTMAX
  endif

  call Broadcast_Scalar(ISO     , master_id)
  call Broadcast_Scalar(TCON    , master_id)
  call Broadcast_Scalar(TBEGIN  , master_id)
  call Broadcast_Scalar(TIDALP  , master_id)
  call Broadcast_Scalar(CF      , master_id)
  call Broadcast_Scalar(ISCORV  , master_id)
  call Broadcast_Scalar(ISDCCA  , master_id)
  call Broadcast_Scalar(ISCFL   , master_id)
  call Broadcast_Scalar(ISCFLM  , master_id)
  call Broadcast_Scalar(DTSSFAC , master_id)
  call Broadcast_Scalar(DTSSDHDT, master_id)
  call Broadcast_Scalar(DTMAX   , master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) TCON, TBEGIN, TIDALP, CF, ISCORV, ISDCCA, ISCFL, ISCFLM, DTSSFAC, DTSSDHDT, DTMAX
  if( ISO > 0 ) GOTO 100

  if( ISRESTI == 1 .and. ICONTINUE == 1 )then
    ! *** FOR CONTINUATION MODE (TBEGINC IS SET FROM THE RESTART FILE)
    NTC = NTC + INT(TBEGIN - TBEGINC)
    TBEGIN = TBEGINC
  endif

  if( DTSSFAC > 0.0 )then
    ISDYNSTP = 1
    DT = TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
    if( DTMAX <= DT*2. ) DTMAX = 3600.
  else
    ISDYNSTP = 0
  endif
  if( IS2TIM == 0 ) ISDYNSTP = 0

  ! *** Calculating number of WASP outputs and snapshot time
  if( RESSTEP > 0 )then
    NWASPOUT = NTC*TIDALP/RESSTEP + 1
    call AllocateDSI(WASPTIME, NWASPOUT, 0.)
    do IT = 1,NWASPOUT
      WASPTIME(IT) = TBEGIN + IT*RESSTEP/DBLE(86400.)
    enddo
  endif

  !C9**  READ SPACE RELATED AND SMOOTHING parameterS
  NCARD = '9'
  ! *** *******************************************************
  if( process_id == master_id )then
    call SEEK('C9',0)
    read(1,*,IOSTAT = ISO) IC_GLOBAL, JC_GLOBAL, LC_GLOBAL, LVC, ldum, NDM, LDM, ISMASK, NBLOCKED, ISCONNECT, ldum, ldum, tmp, tmp
  endif

  call Broadcast_Scalar(IC_GLOBAL, master_id)
  call Broadcast_Scalar(JC_GLOBAL, master_id)
  call Broadcast_Scalar(LC_GLOBAL, master_id)

  call Broadcast_Scalar(ISO      , master_id)
  call Broadcast_Scalar(LVC      , master_id)
  call Broadcast_Scalar(NDM      , master_id)
  call Broadcast_Scalar(LDM      , master_id)
  call Broadcast_Scalar(ISMASK   , master_id)
  call Broadcast_Scalar(NBLOCKED , master_id)
  call Broadcast_Scalar(ISCONNECT, master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) IC_GLOBAL, JC_GLOBAL, LC_GLOBAL, LVC, ldum, NDM, LDM, ISMASK, NBLOCKED, ISCONNECT, ldum, ldum, tmp, tmp
  if( ISO > 0 ) GOTO 100

  LA_Global = LC_Global - 1

  IS2LMC = 0
  if( KC < 0 )then
    KC = -KC
    IS2LMC = 1
  endif

  ! *** DOMAIN DECOMPOSITION CHECKS FOR HORIZONTAL LOOPS
  ! *** MAKE CONSISTENT WITH NEW OMP APPROACH
  ! @todo remove this NDM approach
  LCM2T = LC_GLOBAL - 2

  LCM2T = LC - 2
  NDM = NTHREADS
  LDM = INT(FLOAT(LCM2T)/FLOAT(NTHREADS)) + 1

  if( KC >= 2) ISITB = 0

  ! *** Remapping of IC and JC for domain decomposition with MPI

  ! *** Get the local x,y-ids for identifying where we are on 'node map'
  !   the + 1 is for fortran indexing as the id's will start with 0 otherwise
  x_id = domain_coords(1) + 1
  y_id = domain_coords(2) + 1

  ! *** Set IC and JC for the the subdomain
  IC = IC_DECOMP(x_id)
  if( n_x_partitions > 1 )then
    if( process_map(x_id-1,y_id) /= -1 )  IC = IC + n_ghost_rows   ! *** West Edge
    if( process_map(x_id+1,y_id) /= -1 )  IC = IC + n_ghost_rows   ! *** East Edge
  endif

  JC = JC_DECOMP(y_id)
  if( n_y_partitions > 1 )then
    if( process_map(x_id,y_id-1) /= -1 )  JC = JC + n_ghost_rows   ! *** South Edge
    if( process_map(x_id,y_id+1) /= -1 )  JC = JC + n_ghost_rows   ! *** North Edge
  endif
  call MPI_Barrier(comm_2d, ierr)

  call Parent_Grid

  global_max_width_x = ic_global + 1 + 4
  global_max_width_y = jc_global + 1 + 4

  ! *** After this point IC and JC are Local!

  ! *** write out local IC and JCs
  call WriteBreak(mpi_log_unit)
  write(mpi_log_unit,'(A)') 'After Parent Grid Mapping in INPUT routine '
  write(mpi_log_unit,'(A,I5)') 'IC: ',IC
  write(mpi_log_unit,'(A,I5)') 'JC: ',JC
  write(mpi_log_unit,'(A,I5)') 'Global Max Width x-direction + 4 for ghost: ', global_max_width_x
  write(mpi_log_unit,'(A,I5)') 'Global Max Width x-direction + 4 for ghost: ', global_max_width_y
  call WriteBreak(mpi_log_unit)

  !C9A**  READ VERTICAL SPACE RELATED  parameterS
  NCARD = '9A'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C9A',0)
    read(1,*,IOSTAT = ISO) KC, ldum, ldum, tmp, tmp, ldum
  endif

  call Broadcast_Scalar(ISO    , master_id)
  call Broadcast_Scalar(KC     , master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) KC, ldum, ldum, tmp, tmp, ldum
  if( ISO > 0 ) GOTO 100

  !C10*  READ LAYER THICKNESS IN VERTICAL
  NCARD = '10'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C10',0)
    do K = 1,KC
      read(1,*,IOSTAT = ISO) KDUM, DZCK(K)
    enddo
  endif

  call Broadcast_Scalar(ISO , master_id)
  call Broadcast_Array(DZCK, master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) (DZCK(K),K = 1,KC)
  if( ISO > 0 ) GOTO 100

  !C11*  READ GRID, ROUGHNESS, MASKING AND DEPTH parameterS
  NCARD = '11'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C11',0)
    read(1,*,IOSTAT = ISO) DX, DY, DXYCVT, IMDXDY, ZBRADJ, ZBRCVRT, HMIN, HADADJ, HCVRT, HDRY, HWET, BELADJ, BELCVRT, HDRYMOVE
  endif

  call Broadcast_Scalar(DX     ,  master_id)
  call Broadcast_Scalar(DY     ,  master_id)
  call Broadcast_Scalar(DXYCVT ,  master_id)
  call Broadcast_Scalar(IMDXDY ,  master_id)
  call Broadcast_Scalar(ZBRADJ ,  master_id)
  call Broadcast_Scalar(ZBRCVRT,  master_id)
  call Broadcast_Scalar(HMIN   ,  master_id)
  call Broadcast_Scalar(HADADJ ,  master_id)
  call Broadcast_Scalar(HCVRT  ,  master_id)
  call Broadcast_Scalar(HDRY   ,  master_id)
  call Broadcast_Scalar(HWET   ,  master_id)
  call Broadcast_Scalar(BELADJ ,  master_id)
  call Broadcast_Scalar(BELCVRT,  master_id)
  call Broadcast_Scalar(HDRYMOVE, master_id)

  write(mpi_efdc_out_unit,1002) NCARD
  write(mpi_efdc_out_unit,*) DX, DY, DXYCVT, IMDXDY, ZBRADJ, ZBRCVRT, HMIN, HADADJ, HCVRT, HDRY, HWET, BELADJ, BELCVRT, HDRYMOVE
  if( ISO > 0 ) GOTO 100

  HDRYICE = 0.91*HDRY
  HDRYWAV = 1.2*HWET

  !C11A* READ TWO-LAYER MOMENTUM FLUX AND CURVATURE ACCELERATION
  !     CORRECTION FACTORS
  NCARD = '11A'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C11A',0)
    read(1,*,IOSTAT = ISO) ICK2COR, CK2UUM, CK2VVM, CK2UVM, CK2UUC, CK2VVC, CK2UVC, CK2FCX, CK2FCY
    write(mpi_efdc_out_unit,1002) NCARD

    write(mpi_efdc_out_unit,*) ICK2COR, CK2UUM, CK2VVM, CK2UVM, CK2UUC, CK2VVC, CK2UVC, CK2FCX, CK2FCY
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(ICK2COR, master_id)
  call Broadcast_Scalar(CK2UUM , master_id)
  call Broadcast_Scalar(CK2VVM , master_id)
  call Broadcast_Scalar(CK2UVM , master_id)
  call Broadcast_Scalar(CK2UUC , master_id)
  call Broadcast_Scalar(CK2VVC , master_id)
  call Broadcast_Scalar(CK2UVC , master_id)
  call Broadcast_Scalar(CK2FCX , master_id)
  call Broadcast_Scalar(CK2FCY , master_id)

  if( ICK2COR >= 1 )then
    IS2LMC = ICK2COR
  endif

  !C11B* READ CORNER CELL BOTTOM STRESS CORRECTION OPTIONS
  NCARD = '11B'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C11B',0)
    read(1,*,IOSTAT = ISO) ISCORTBC, ISCORTBCD, FSCORTBC

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ISCORTBC,ISCORTBCD,FSCORTBC
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(ISCORTBC, master_id)
  call Broadcast_Scalar(ISCORTBCD,master_id)
  call Broadcast_Scalar(FSCORTBC, master_id)

  !C12*  READ TURBULENT DIFFUSION parameterS
  NCARD = '12'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C12',0)
    read(1,*,IOSTAT = ISO) AHO, AHD, AVO, ABO, AVMX, ABMX, VISMUD, AVCON, ZBRWALL

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) AHO, AHD, AVO, ABO, AVMX, ABMX, VISMUD, AVCON, ZBRWALL
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(AHO    , master_id)
  call Broadcast_Scalar(AHD    , master_id)
  call Broadcast_Scalar(AVO    , master_id)
  call Broadcast_Scalar(ABO    , master_id)
  call Broadcast_Scalar(AVMX   , master_id)
  call Broadcast_Scalar(ABMX   , master_id)
  call Broadcast_Scalar(VISMUD , master_id)
  call Broadcast_Scalar(AVCON  , master_id)
  call Broadcast_Scalar(ZBRWALL, master_id)

  ISAVCOMP = 1
  if( AVCON == 0. ) ISAVCOMP = 0

  !C12A*  READ TURBULENCE CLOSURE OPTIONS
  NCARD = '12A'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C12A',0)
    read(1,*,IOSTAT = ISO) ISTOPT(0), ISSQL, ISAVBMX, ISFAVB, ISINWV, ISLLIM, IFPROX, XYRATIO, BC_EDGEFACTOR

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ISTOPT(0), ISSQL, ISAVBMX, ISFAVB, ISINWV, ISLLIM, IFPROX, XYRATIO, BC_EDGEFACTOR
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(ISTOPT(0)    , master_id)
  call Broadcast_Scalar(ISSQL        , master_id)
  call Broadcast_Scalar(ISAVBMX      , master_id)
  call Broadcast_Scalar(ISFAVB       , master_id)
  call Broadcast_Scalar(ISINWV       , master_id)
  call Broadcast_Scalar(ISLLIM       , master_id)
  call Broadcast_Scalar(IFPROX       , master_id)
  call Broadcast_Scalar(XYRATIO      , master_id)
  call Broadcast_Scalar(BC_EDGEFACTOR, master_id)

  if( BC_EDGEFACTOR < 0 ) BC_EDGEFACTOR = 0.0
  if( BC_EDGEFACTOR > 1 ) BC_EDGEFACTOR = 1.0

  !C12B*  READ GOTM TURBULENCE OPTIONS
  NCARD = '12B'

  if( process_id == master_id )then
    call SEEK('C12B',0)
    read(1,*,IOSTAT = ISO) ISGOTM, IFRICTION, ICALNN, ICALSS, CHARNOCK

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*)  ISGOTM, IFRICTION, ICALNN, ICALSS, CHARNOCK
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(ISGOTM,        master_id)
  call Broadcast_Scalar(IFRICTION,     master_id)
  call Broadcast_Scalar(ICALNN,        master_id)
  call Broadcast_Scalar(ICALSS,        master_id)
  call Broadcast_Scalar(CHARNOCK,      master_id)

  ! *** Initialize GOTM variables
  if( ISGOTM > 0 ) CALL Init_GOTM

  ! *** DELME IGOTMTEST
  INQUIRE(FILE = 'GOTMTEST.INP',EXIST = RES)
  if( RES )then
    if( process_id == master_id )then
      open(99,FILE = 'GOTMTEST.INP',STATUS = 'UNKNOWN')
      STR = READSTR(99)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(99,*)  iGOTM_Test, SurfaceShearX, SurfaceShearY, HeatFlux
      close(99)
    endif
    call Broadcast_Scalar(iGOTM_Test,    master_id)
    call Broadcast_Scalar(SurfaceShearX, master_id)
    call Broadcast_Scalar(SurfaceShearY, master_id)
    call Broadcast_Scalar(HeatFlux,      master_id)

    write(mpi_efdc_out_unit,*)  ' FILE = GOTMTEST.INP'
    write(mpi_efdc_out_unit,*)  iGOTM_Test, SurfaceShearX, SurfaceShearY, HeatFlux
  endif
  ! *** DELME

  !C13*  READ TURBULENCE CLOSURE parameterS
  NCARD = '13'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C13',0)
    ! *** PMC - CTE2 NOT USED
    read(1,*,IOSTAT = ISO) VKC, CTURB, CTURB2B, CTE1, CTE2, CTE3, CTE4, CTE5, RIQMAX, QQMIN, QQLMIN, DMLMIN

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) VKC, CTURB, CTURB2B, CTE1, CTE2, CTE3, CTE4, CTE5, RIQMAX, QQMIN, QQLMIN, DMLMIN
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(VKC      , master_id)
  call Broadcast_Scalar(CTURB    , master_id)
  call Broadcast_Scalar(CTURB2B  , master_id)
  call Broadcast_Scalar(CTE1     , master_id)
  call Broadcast_Scalar(CTE2     , master_id)
  call Broadcast_Scalar(CTE3     , master_id)
  call Broadcast_Scalar(CTE4     , master_id)
  call Broadcast_Scalar(CTE5     , master_id)
  call Broadcast_Scalar(RIQMAX   , master_id)
  call Broadcast_Scalar(QQMIN    , master_id)
  call Broadcast_Scalar(QQLMIN   , master_id)
  call Broadcast_Scalar(DMLMIN   , master_id)

  !C14*  READ TIDAL & ATMOSPHERIC FORCING, GROUND WATER AND SUBGRID CHANNEL parameterS
  NCARD = '14'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C14',0)
    read(1,*,IOSTAT = ISO) MTIDE, NWSER, NASER, ISGWIT, ISCHAN, ISWAVE, ITIDASM, ISPERC, ISBODYF, ISPNHYDS, ISPROPWASH

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) MTIDE, NWSER, NASER, ISGWIT, ISCHAN, ISWAVE, ITIDASM, ISPERC, ISBODYF, ISPNHYDS, ISPROPWASH
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(MTIDE     , master_id)
  call Broadcast_Scalar(NWSER     , master_id)
  call Broadcast_Scalar(NASER     , master_id)
  call Broadcast_Scalar(ISGWIT    , master_id)
  call Broadcast_Scalar(ISCHAN    , master_id)
  call Broadcast_Scalar(ISWAVE    , master_id)
  call Broadcast_Scalar(ITIDASM   , master_id)
  call Broadcast_Scalar(ISPERC    , master_id)
  call Broadcast_Scalar(ISBODYF   , master_id)
  call Broadcast_Scalar(ISPNHYDS  , master_id)
  call Broadcast_Scalar(ISPROPWASH, master_id)

  ISWCBL = 0
  ISWVSD = 0
  if( ISPERC > 0 ) ISGWIT = 3

  ! *** INITIALIZE PROPWASH CALCULATIONS FLAG
  if( ISPROPWASH > 0 ) propwash_on = .TRUE.

  !C14A* READ SAND GRAIN NIKURADSE ROUGHNESS
  KSW = 0.00001 * 2.5  ! *** DEFAULT IS 10 MICRON

  if( ISWAVE >= 1 )then
    NCARD = '14A'
    if( process_id == master_id )then
      call SEEK('C14A',0)
      read(1,*,IOSTAT = ISO) KSW, IUSEWVCELLS, IFWAVE, SWANGRP, ISSTEAD
      if( ISO > 0 ) GOTO 100
    endif

    call Broadcast_Scalar(KSW        , master_id)
    call Broadcast_Scalar(IUSEWVCELLS, master_id)
    call Broadcast_Scalar(IFWAVE     , master_id)
    call Broadcast_Scalar(SWANGRP    , master_id)
    call Broadcast_Scalar(ISSTEAD    , master_id)

    if( ISWAVE == 1 .or. ISWAVE == 2 .or. ISWAVE == 4 )then
      NCARD = '14B'
      if( process_id == master_id )then
        call SEEK('C14B',0)
        read(1,*,IOSTAT = ISO) ISWRSR, ISWRSI, WVDISV, WVLSH, WVLSX, ISWVSD, WVLCAL, NTSWV, ISWCBL, ISDZBR
        if( ISO > 0 ) GOTO 100
      endif

      call Broadcast_Scalar(ISWRSR, master_id)
      call Broadcast_Scalar(ISWRSI, master_id)
      call Broadcast_Scalar(WVDISV, master_id)
      call Broadcast_Scalar(WVLSH , master_id)
      call Broadcast_Scalar(WVLSX , master_id)
      call Broadcast_Scalar(ISWVSD, master_id)
      call Broadcast_Scalar(WVLCAL, master_id)
      call Broadcast_Scalar(NTSWV , master_id)
      call Broadcast_Scalar(ISWCBL, master_id)
      call Broadcast_Scalar(ISDZBR, master_id)

      NTSWV = MAX(NTSWV,1)
    endif
  endif

  !C14C TIME & SPACE VARYING FORCINGS     2018-10-12, NTL:
  NCARD = '14C'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C14C',0)
    read(1,*,IOSTAT = ISO) BATHY.IFLAG, ROUGH.IFLAG, VEGE.IFLAG,  GWSP.IFLAG,  WIND.IFLAG, PRESS.IFLAG,  &
      RAIN.IFLAG,  EVAP.IFLAG,  SHELT.IFLAG, SHADE.IFLAG, SNOW.IFLAG, ICETHK.IFLAG, &
      SEDZLJER.IFLAG

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) BATHY.IFLAG, ROUGH.IFLAG, VEGE.IFLAG,  GWSP.IFLAG,  WIND.IFLAG, PRESS.IFLAG,  &
      RAIN.IFLAG,  EVAP.IFLAG,  SHELT.IFLAG, SHADE.IFLAG, SNOW.IFLAG, ICETHK.IFLAG, &
      SEDZLJER.IFLAG
  endif

  call Broadcast_Scalar(BATHY.IFLAG  , master_id)
  call Broadcast_Scalar(ROUGH.IFLAG , master_id)
  call Broadcast_Scalar(VEGE.IFLAG  , master_id)
  call Broadcast_Scalar(GWSP.IFLAG  , master_id)
  call Broadcast_Scalar(WIND.IFLAG  , master_id)
  call Broadcast_Scalar(PRESS.IFLAG , master_id)
  call Broadcast_Scalar(RAIN.IFLAG  , master_id)
  call Broadcast_Scalar(EVAP.IFLAG  , master_id)
  call Broadcast_Scalar(SHELT.IFLAG , master_id)
  call Broadcast_Scalar(SHADE.IFLAG , master_id)
  call Broadcast_Scalar(SNOW.IFLAG , master_id)
  call Broadcast_Scalar(ICETHK.IFLAG, master_id)
  call Broadcast_Scalar(SEDZLJER.IFLAG, master_id)

  !C14C  PARAMETRIC CYCLONE TIME & SPACE VARYING FORCING     2021-04-23, NTL:
  NCARD = '14D'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C14D',0)
    read(1,*,IOSTAT = ISO) ICYCLONE, RHO_A, PINF, THETAMAX, WDRAG1, WDRAG2, CDRAG1, CDRAG2

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ICYCLONE, RHO_A, PINF, THETAMAX, WDRAG1, WDRAG2, CDRAG1, CDRAG2
  endif

  call Broadcast_Scalar(ICYCLONE, master_id)
  call Broadcast_Scalar(RHO_A,    master_id)
  call Broadcast_Scalar(PINF,     master_id)
  call Broadcast_Scalar(THETAMAX, master_id)
  call Broadcast_Scalar(WDRAG1,   master_id)
  call Broadcast_Scalar(WDRAG2,   master_id)
  call Broadcast_Scalar(CDRAG1,   master_id)
  call Broadcast_Scalar(CDRAG2,   master_id)

  if( MTIDE > 0 )then
    !C15*  READ PERIODIC FORCING (TIDAL) CONSTITUENT SYMBOLS AND PERIODS
    NCARD = '15'
    if( process_id == master_id )then
      call SEEK('C15',0)
      do M = 1,MTIDE
        read(1,*,IOSTAT = ISO) SYMBOL(M),TCP(M)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) SYMBOL(M),TCP(M)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(SYMBOL, master_id)
    call Broadcast_Array(TCP,    master_id)

  endif

  !C16*  READ SURFACE ELEVATION OR PRESSURE BOUNDARY CONDITION parameterS
  NCARD = '16'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C16',0)
    read(1,*,IOSTAT = ISO) NPBS, NPBW, NPBE, NPBN, NPFOR_Readin, NPFORT, NPSER, PDGINIT

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) NPBS, NPBW, NPBE, NPBN, NPFOR_Readin, NPFORT, NPSER, PDGINIT
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(NPBS        , master_id)
  call Broadcast_Scalar(NPBW        , master_id)
  call Broadcast_Scalar(NPBE        , master_id)
  call Broadcast_Scalar(NPBN        , master_id)
  call Broadcast_Scalar(NPFOR_Readin, master_id)
  call Broadcast_Scalar(NPFORT      , master_id)
  call Broadcast_Scalar(NPSER       , master_id)
  call Broadcast_Scalar(PDGINIT     , master_id)

  if( NPFOR_Readin > 0 )then
    ! ***  READ PERIODIC FORCING (TIDAL) SURFACE ELEVATION OR
    NCARD = '17'
    call AllocateDSI( CPFAM0,   NPFORM, MTM, 0.0)
    call AllocateDSI( CPFAM1,   NPFORM, MTM, 0.0)
    call AllocateDSI( CPFAM2,   NPFORM, MTM, 0.0)
    call AllocateDSI( SPFAM0,   NPFORM, MTM, 0.0)
    call AllocateDSI( SPFAM1,   NPFORM, MTM, 0.0)
    call AllocateDSI( SPFAM2,   NPFORM, MTM, 0.0)
    call AllocateDSI( PFX2,     NPFORM, MTM, 0.0)

    if( process_id == master_id )then
      call SEEK('C17',0)
      do NP = 1,NPFOR_Readin
        do M = 1,MTIDE
          if( NPFORT == 0 )then
            read(1,*,IOSTAT = ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
            write(mpi_efdc_out_unit,1002) NCARD
            write(mpi_efdc_out_unit,*) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M)
            if( ISO > 0 ) GOTO 100
          elseif( NPFORT >= 1 )then
            read(1,*,IOSTAT = ISO) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M)
            RAD = PI2*PFPH(NP,M)/TCP(M)
            CPFAM0(NP,M) = PFAM(NP,M)*COS(RAD)
            SPFAM0(NP,M) = PFAM(NP,M)*SIN(RAD)
            write(mpi_efdc_out_unit,1002) NCARD
            write(mpi_efdc_out_unit,*) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M), CPFAM0(NP,M), SPFAM0(NP,M)
            if( ISO > 0 ) GOTO 100
            read(1,*,IOSTAT = ISO) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M)
            RAD = PI2*PFPH(NP,M)/TCP(M)
            CPFAM1(NP,M) = PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
            SPFAM1(NP,M) = PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
            write(mpi_efdc_out_unit,1002) NCARD
            write(mpi_efdc_out_unit,*) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M), CPFAM1(NP,M), SPFAM1(NP,M)
            CPFAM2(NP,M) = 0.0
            SPFAM2(NP,M) = 0.0
          elseif( NPFORT == 2 )then
            read(1,*,IOSTAT = ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),PFX2(NP,M)
            RAD = PI2*PFPH(NP,M)/TCP(M)
            if( PFX2(NP,M)>0.0 )then
              CPFAM2(NP,M) = PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
              SPFAM2(NP,M) = PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
            else
              CPFAM2(NP,M) = 0.
              SPFAM2(NP,M) = 0.
            endif
          endif

        enddo
      enddo

    endif

    if( NPFORT == 0 )then
      call Broadcast_Array(PFAM,   master_id)
      call Broadcast_Array(PFPH,   master_id)
    elseif( NPFORT >= 1 )then
      call Broadcast_Array(PFAM,   master_id)
      call Broadcast_Array(PFPH,   master_id)
      call Broadcast_Array(CPFAM0, master_id)
      call Broadcast_Array(SPFAM0, master_id)
      call Broadcast_Array(CPFAM1, master_id)
      call Broadcast_Array(SPFAM1, master_id)
      call Broadcast_Array(CPFAM2, master_id)
      call Broadcast_Array(SPFAM2, master_id)
    elseif( NPFORT == 2 )then
      call Broadcast_Array(PFAM,   master_id)
      call Broadcast_Array(PFPH,   master_id)
      call Broadcast_Array(PFX2,   master_id)
      call Broadcast_Array(CPFAM2, master_id)
      call Broadcast_Array(SPFAM2, master_id)
    endif

  endif

  if( NPBS > 0 )then
    !C18*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
    !     ON SOUTH OPEN BOUNDARIES
    NCARD = '18'
    if( process_id == master_id )then
      call SEEK('C18',0)
      if( NPFORT == 0 )then
        do L = 1,NPBS
          ! *** Added global arrays for MPI
          read(1,*,IOSTAT = ISO) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L)
          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L)
          if( ISO > 0 ) GOTO 100

          do M = 1,MTIDE
            if( NPFORS == 0) EXIT
            RAD = PI2*PFPH(NPFORS,M)/TCP(M)
            AMP = G*PFAM(NPFORS,M)
            PCBS_GL(L,M) = AMP*COS(RAD)
            PSBS_GL(L,M) = AMP*SIN(RAD)
          enddo
        enddo

      elseif( NPFORT == 1 )then
        do L = 1,NPBS
          read(1,*,IOSTAT = ISO) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L), NPSERS1_GL(L), TPCOORDS_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L), NPSERS1_GL(L), TPCOORDS_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORS == 0) EXIT
            PCBS_GL(L,M) = CPFAM0(NPFORS,M)+TPCOORDS_GL(L)*CPFAM1(NPFORS,M) + TPCOORDS_GL(L)*TPCOORDS_GL(L)*CPFAM2(NPFORS,M)
            PSBS_GL(L,M) = SPFAM0(NPFORS,M)+TPCOORDS_GL(L)*SPFAM1(NPFORS,M) + TPCOORDS_GL(L)*TPCOORDS_GL(L)*SPFAM2(NPFORS,M)
            TMPAMP = SQRT(PCBS(L,M)*PCBS(L,M)+PSBS_GL(L,M)*PSBS_GL(L,M))
            TMPPHS = ATAN2(PSBS_GL(L,M),PCBS_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBS_GL(L,M) = G*PCBS_GL(L,M)
            PSBS_GL(L,M) = G*PSBS_GL(L,M)
          enddo
        enddo

      elseif( NPFORT == 2 )then
        do L = 1,NPBS
          read(1,*,IOSTAT = ISO) IPBS_GL(L),JPBS_GL(L),ISPBS_GL(L), ISPRS_GL(L), NPFORS,NPSERS_GL(L),NPSERS1_GL(L), TPCOORDS_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L), NPSERS1_GL(L), TPCOORDS_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORS == 0) EXIT
            BOTTOM = PFX2(NPFORS,M)*(1.0-PFX2(NPFORS,M))
            TOP1 = TPCOORDS_GL(L)*PFX2(NPFORS,M)*(TPCOORDS_GL(L)-PFX2(NPFORS,M))
            TOP2 = TPCOORDS_GL(L)*(1.0-TPCOORDS_GL(L))
            if( BOTTOM == 0.0 )then
              TOP1 = TPCOORDS_GL(L)
              TOP2 = TPCOORDS_GL(L)*TPCOORDS_GL(L)
            else
              TOP1 = TOP1/BOTTOM
              TOP2 = TOP2/BOTTOM
            endif
            PCBS_GL(L,M) = CPFAM0(NPFORS,M)+TOP1*CPFAM1(NPFORS,M)+TOP2*CPFAM2(NPFORS,M)
            PSBS_GL(L,M) = SPFAM0(NPFORS,M)+TOP1*SPFAM1(NPFORS,M)+TOP2*SPFAM2(NPFORS,M)
            TMPAMP = SQRT(PCBS_GL(L,M)*PCBS_GL(L,M)+PSBS_GL(L,M)*PSBS_GL(L,M))
            TMPPHS = ATAN2(PSBS_GL(L,M),PCBS_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBS_GL(L,M) = G*PCBS_GL(L,M)
            PSBS_GL(L,M) = G*PSBS_GL(L,M)
          enddo
        enddo
      endif
2068  FORMAT(I4,3X,A2,5X,E14.4,3E14.5,5X,2I5)
2069  FORMAT(I4,3X,A2,5X,2E14.4,5X,2I5)

    endif

    call Broadcast_Array(IPBS_GL,     master_id)
    call Broadcast_Array(JPBS_GL,     master_id)
    call Broadcast_Array(ISPBS_GL,    master_id)
    call Broadcast_Array(ISPRS_GL,    master_id)
    call Broadcast_Array(NPSERS_GL,   master_id)
    call Broadcast_Array(PCBS_GL,     master_id)
    call Broadcast_Array(PSBS_GL,     master_id)
    if( NPFORT > 0 )then
      call Broadcast_Array(NPSERS1_GL,  master_id)
      call Broadcast_Array(TPCOORDS_GL, master_id)
    endif
  endif

  if( NPBW > 0 )then
    !C19*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
    !     ON WEST OPEN BOUNDARIES
    NCARD = '19'
    if( process_id == master_id )then

      call SEEK('C19',0)
      if( NPFORT == 0 )then
        do L = 1,NPBW
          read(1,*,IOSTAT = ISO) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L)
          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORW == 0) EXIT
            RAD = PI2*PFPH(NPFORW,M)/TCP(M)
            AMP = G*PFAM(NPFORW,M)
            PCBW_GL(L,M) = AMP*COS(RAD)
            PSBW_GL(L,M) = AMP*SIN(RAD)
          enddo
        enddo

      elseif( NPFORT == 1 )then
        do L = 1,NPBW
          read(1,*,IOSTAT = ISO) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORW == 0) EXIT
            PCBW_GL(L,M) = CPFAM0(NPFORW,M)+TPCOORDW_GL(L)*CPFAM1(NPFORW,M)+TPCOORDW_GL(L)*TPCOORDW_GL(L)*CPFAM2(NPFORW,M)
            PSBW_GL(L,M) = SPFAM0(NPFORW,M)+TPCOORDW_GL(L)*SPFAM1(NPFORW,M)+TPCOORDW_GL(L)*TPCOORDW_GL(L)*SPFAM2(NPFORW,M)
            TMPAMP = SQRT(PCBW_GL(L,M)*PCBW_GL(L,M)+PSBW_GL(L,M)*PSBW_GL(L,M))
            TMPPHS = 0.0
            if( TMPAMP>0.0) TMPPHS = ATAN2(PSBW_GL(L,M),PCBW_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBW_GL(L,M) = G*PCBW_GL(L,M)
            PSBW_GL(L,M) = G*PSBW_GL(L,M)
          enddo
        enddo

      elseif( NPFORT == 2 )then
        do L = 1,NPBW
          read(1,*,IOSTAT = ISO) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORW == 0) EXIT
            BOTTOM = PFX2(NPFORW,M)*(1.0-PFX2(NPFORW,M))
            TOP1 = TPCOORDW_GL(L)*PFX2(NPFORW,M)*(TPCOORDW_GL(L)-PFX2(NPFORW,M))
            TOP2 = TPCOORDW_GL(L)*(1.0-TPCOORDW_GL(L))
            if( BOTTOM == 0.0 )then
              TOP1 = TPCOORDW_GL(L)
              TOP2 = TPCOORDW_GL(L)*TPCOORDW_GL(L)
            else
              TOP1 = TOP1/BOTTOM
              TOP2 = TOP2/BOTTOM
            endif
            PCBW_GL(L,M) = CPFAM0(NPFORW,M)+TOP1*CPFAM1(NPFORW,M)+TOP2*CPFAM2(NPFORW,M)
            PSBW_GL(L,M) = SPFAM0(NPFORW,M)+TOP1*SPFAM1(NPFORW,M)+TOP2*SPFAM2(NPFORW,M)
            TMPAMP = SQRT(PCBW_GL(L,M)*PCBW_GL(L,M)+PSBW_GL(L,M)*PSBW_GL(L,M))
            TMPPHS = 0.0
            if( TMPAMP>0.0) TMPPHS = ATAN2(PSBW_GL(L,M),PCBW_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBW_GL(L,M) = G*PCBW_GL(L,M)
            PSBW_GL(L,M) = G*PSBW_GL(L,M)
          enddo
        enddo
      endif
    endif

    call Broadcast_Array(IPBW_GL,     master_id)
    call Broadcast_Array(JPBW_GL,     master_id)
    call Broadcast_Array(ISPBW_GL,    master_id)
    call Broadcast_Array(ISPRW_GL,    master_id)
    call Broadcast_Array(NPSERW_GL,   master_id)
    call Broadcast_Array(PCBW_GL,     master_id)
    call Broadcast_Array(PSBW_GL,     master_id)
    if( NPFORT > 0 )then
      call Broadcast_Array(NPSERW1_GL,  master_id)
      call Broadcast_Array(TPCOORDW_GL, master_id)
    endif
  endif

  if( NPBE > 0 )then
    !C20*  READ PERIODIC FORCING (TIDAL)ELEVATION BOUNDARY CONDTIONS
    !     ON EAST OPEN BOUNDARIES
    NCARD = '20'
    if( process_id == master_id )then
      call SEEK('C20',0)
      if( NPFORT == 0 )then
        do L = 1,NPBE
          read(1,*,IOSTAT = ISO) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORE == 0) EXIT
            RAD = PI2*PFPH(NPFORE,M)/TCP(M)
            AMP = G*PFAM(NPFORE,M)
            PCBE_GL(L,M) = AMP*COS(RAD)
            PSBE_GL(L,M) = AMP*SIN(RAD)
          enddo
        enddo

      elseif( NPFORT == 1 )then
        do L = 1,NPBE
          read(1,*,IOSTAT = ISO) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORE == 0) EXIT
            PCBE_GL(L,M) = CPFAM0(NPFORE,M)+TPCOORDE_GL(L)*CPFAM1(NPFORE,M)+TPCOORDE_GL(L)*TPCOORDE_GL(L)*CPFAM2(NPFORE,M)
            PSBE_GL(L,M) = SPFAM0(NPFORE,M)+TPCOORDE_GL(L)*SPFAM1(NPFORE,M)+TPCOORDE_GL(L)*TPCOORDE_GL(L)*SPFAM2(NPFORE,M)
            TMPAMP = SQRT(PCBE_GL(L,M)*PCBE_GL(L,M)+PSBE_GL(L,M)*PSBE_GL(L,M))
            TMPPHS = 0.0
            if( TMPAMP>0.0) TMPPHS = ATAN2(PSBE_GL(L,M),PCBE_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBE_GL(L,M) = G*PCBE_GL(L,M)
            PSBE_GL(L,M) = G*PSBE_GL(L,M)
          enddo
        enddo

      elseif( NPFORT == 2 )then
        do L = 1,NPBE
          read(1,*,IOSTAT = ISO) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORE == 0) EXIT
            BOTTOM = PFX2(NPFORE,M)*(1.0-PFX2(NPFORE,M))
            TOP1 = TPCOORDE_GL(L)*PFX2(NPFORE,M)*(TPCOORDE_GL(L)-PFX2(NPFORE,M))
            TOP2 = TPCOORDE_GL(L)*(1.0-TPCOORDE_GL(L))
            if( BOTTOM == 0.0 )then
              TOP1 = TPCOORDE_GL(L)
              TOP2 = TPCOORDE_GL(L)*TPCOORDE_GL(L)
            else
              TOP1 = TOP1/BOTTOM
              TOP2 = TOP2/BOTTOM
            endif
            PCBE_GL(L,M) = CPFAM0(NPFORE,M)+TOP1*CPFAM1(NPFORE,M)+TOP2*CPFAM2(NPFORE,M)
            PSBE_GL(L,M) = SPFAM0(NPFORE,M)+TOP1*SPFAM1(NPFORE,M)+TOP2*SPFAM2(NPFORE,M)
            TMPAMP = SQRT(PCBE_GL(L,M)*PCBE_GL(L,M)+PSBE_GL(L,M)*PSBE_GL(L,M))
            TMPPHS = 0.0
            if( TMPAMP>0.0) TMPPHS = ATAN2(PSBE_GL(L,M),PCBE_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBE_GL(L,M) = G*PCBE_GL(L,M)
            PSBE_GL(L,M) = G*PSBE_GL(L,M)
          enddo
        enddo
      endif
    endif

    call Broadcast_Array(IPBE_GL,     master_id)
    call Broadcast_Array(JPBE_GL,     master_id)
    call Broadcast_Array(ISPBE_GL,    master_id)
    call Broadcast_Array(ISPRE_GL,    master_id)
    call Broadcast_Array(NPSERE_GL,   master_id)
    call Broadcast_Array(PCBE_GL,     master_id)
    call Broadcast_Array(PSBE_GL,     master_id)
    if( NPFORT > 0 )then
      call Broadcast_Array(NPSERE1_GL,  master_id)
      call Broadcast_Array(TPCOORDE_GL, master_id)
    endif
  endif

  if( NPBN > 0 )then
    !C21*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
    !     ON NORTH OPEN BOUNDARIES
    NCARD = '21'
    if( process_id == master_id )then

      call SEEK('C21',0)
      if( NPFORT == 0 )then
        do L = 1,NPBN
          read(1,*,IOSTAT = ISO) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORN == 0) EXIT
            RAD = PI2*PFPH(NPFORN,M)/TCP(M)
            AMP = G*PFAM(NPFORN,M)
            PCBN_GL(L,M) = AMP*COS(RAD)
            PSBN_GL(L,M) = AMP*SIN(RAD)
          enddo
        enddo

      elseif( NPFORT >= 1 )then
        do L = 1,NPBN
          read(1,*,IOSTAT = ISO) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORN == 0) EXIT
            PCBN_GL(L,M) = CPFAM0(NPFORN,M)+TPCOORDN_GL(L)*CPFAM1(NPFORN,M)+TPCOORDN_GL(L)*TPCOORDN_GL(L)*CPFAM2(NPFORN,M)
            PSBN_GL(L,M) = SPFAM0(NPFORN,M)+TPCOORDN_GL(L)*SPFAM1(NPFORN,M)+TPCOORDN_GL(L)*TPCOORDN_GL(L)*SPFAM2(NPFORN,M)
            TMPAMP = SQRT(PCBN_GL(L,M)*PCBN_GL(L,M)+PSBN_GL(L,M)*PSBN_GL(L,M))
            TMPPHS = 0.0
            if( TMPAMP>0.0) TMPPHS = ATAN2(PSBN_GL(L,M),PCBN_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBN_GL(L,M) = G*PCBN_GL(L,M)
            PSBN_GL(L,M) = G*PSBN_GL(L,M)
          enddo
        enddo

      elseif( NPFORT == 2 )then
        do L = 1,NPBN
          read(1,*,IOSTAT = ISO) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)
          if( ISO > 0 ) GOTO 100
          do M = 1,MTIDE
            if( NPFORN == 0) EXIT
            BOTTOM = PFX2(NPFORN,M)*(1.0-PFX2(NPFORN,M))
            TOP1 = TPCOORDN_GL(L)*PFX2(NPFORN,M)*(TPCOORDN_GL(L)-PFX2(NPFORN,M))
            TOP2 = TPCOORDN_GL(L)*(1.0-TPCOORDN_GL(L))
            if( BOTTOM == 0.0 )then
              TOP1 = TPCOORDN_GL(L)
              TOP2 = TPCOORDN_GL(L)*TPCOORDN_GL(L)
            else
              TOP1 = TOP1/BOTTOM
              TOP2 = TOP2/BOTTOM
            endif
            PCBN_GL(L,M) = CPFAM0(NPFORN,M)+TOP1*CPFAM1(NPFORN,M)+TOP2*CPFAM2(NPFORN,M)
            PSBN_GL(L,M) = SPFAM0(NPFORN,M)+TOP1*SPFAM1(NPFORN,M)+TOP2*SPFAM2(NPFORN,M)
            TMPAMP = SQRT(PCBN_GL(L,M)*PCBN_GL(L,M)+PSBN_GL(L,M)*PSBN_GL(L,M))
            TMPPHS = 0.0
            if( TMPAMP>0.0) TMPPHS = ATAN2(PSBN_GL(L,M),PCBN_GL(L,M))
            TMPPHS = TMPPHS*TCP(M)/PI2
            if( TMPPHS<0.0)TMPPHS = TMPPHS+TCP(M)
            PCBN_GL(L,M) = G*PCBN_GL(L,M)
            PSBN_GL(L,M) = G*PSBN_GL(L,M)
          enddo
        enddo
      endif
    endif

    call Broadcast_Array(IPBN_GL,     master_id)
    call Broadcast_Array(JPBN_GL,     master_id)
    call Broadcast_Array(ISPBN_GL,    master_id)
    call Broadcast_Array(ISPRN_GL,    master_id)
    call Broadcast_Array(NPSERN_GL,   master_id)
    call Broadcast_Array(PCBN_GL,     master_id)
    call Broadcast_Array(PSBN_GL,     master_id)
    if( NPFORT > 0 )then
      call Broadcast_Array(NPSERN1_GL,  master_id)
      call Broadcast_Array(TPCOORDN_GL, master_id)
    endif
  endif

  !22*  READ NUM OF SEDIMENT AMD TOXICS AND NUM OF CONCENTRATION TIME SERIES
  NCARD = '22'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C22',0)
    read(1,*,IOSTAT = ISO) NDYE,NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),NCSER(4),NCSER(5),NCSER(6),NCSER(7),ldum

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) NDYE,NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),NCSER(4),NCSER(5),NCSER(6),NCSER(7),ldum
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(NDYE  , master_id)
  call Broadcast_Scalar(NTOX  , master_id)
  call Broadcast_Scalar(NSED  , master_id)
  call Broadcast_Scalar(NSND  , master_id)
  call Broadcast_Array(NCSER(1:7), master_id)

  ! *** REMOVE UNUSED SETTINGS TO ALLOW FOR SELECTIVE ALLOCATIONS
  if( ISTRAN(6) < 1)  NSED = 0
  if( ISTRAN(7) < 1 ) NSND = 0
  if( NSED == 0 .and. NSND == 0 ) ISTRAN(5) = 0
  if( ISTRAN(5) < 1 ) NTOX = 0
  if( NTOX == 0 ) NCSER(5) = 0
  if( NSED == 0 ) NCSER(6) = 0
  if( NSND == 0 ) NCSER(7) = 0

  !C23*  READ VELOCITY, VOL SOUR/SINK, FLOW CONTROL, & WITHDRAW/RETURN DATA
  NCARD = '23'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C23',0)
    read(1,*,IOSTAT = ISO) NQSIJ, NQJPIJ, NQSER, NQCTL, NQCTLT, NHYDST, NQWR, NQWRSR, ISDIQ, NQCTLSER, NQCTRULES

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) NQSIJ, NQJPIJ, NQSER, NQCTL, NQCTLT, NHYDST, NQWR, NQWRSR, ISDIQ, NQCTLSER, NQCTRULES
    if( ISO > 0 ) GOTO 100
  endif
  ! *** Broadcast in Scan_EFDC

  NGRPID = 0
  if( NQSIJ > 0 )then
    !C24*  READ VOLUME SOURCE/SINK LOCATIONS, MAGNITUDES, & VOL & CONC SERIES
    NCARD = '24'
    if( process_id == master_id )then
      call SEEK('C24',0)
      do NS = 1,NQSIJ
        read(1,*,IOSTAT = ISO) BCFL_GL(NS).I, BCFL_GL(NS).J, BCFL_GL(NS).QSSE, BCFL_GL(NS).NQSMUL, BCFL_GL(NS).NQSMF, BCFL_GL(NS).NQSERQ,   BCFL_GL(NS).NCSERQ(1), BCFL_GL(NS).NCSERQ(2), BCFL_GL(NS).NCSERQ(3), BCFL_GL(NS).NCSERQ(4), BCFL_GL(NS).NCSERQ(5),  &
          BCFL_GL(NS).NCSERQ(6), BCFL_GL(NS).NCSERQ(7), BCFL_GL(NS).QWIDTH,   BCFL_GL(NS).QFACTOR,  BCFL_GL(NS).GRPID

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) BCFL_GL(NS).I, BCFL_GL(NS).J, BCFL_GL(NS).QSSE, BCFL_GL(NS).NQSMUL, BCFL_GL(NS).NQSMF, BCFL_GL(NS).NQSERQ,   BCFL_GL(NS).NCSERQ(1), BCFL_GL(NS).NCSERQ(2), BCFL_GL(NS).NCSERQ(3), BCFL_GL(NS).NCSERQ(4), BCFL_GL(NS).NCSERQ(5),  &
          BCFL_GL(NS).NCSERQ(6), BCFL_GL(NS).NCSERQ(7), BCFL_GL(NS).QWIDTH,   BCFL_GL(NS).QFACTOR,  BCFL_GL(NS).GRPID
        if( ISO > 0 ) GOTO 100

      enddo
    endif

    !C25*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS
    !     SAL,TEM,DYE,SFL,TOX(1 TO NOTX)
    NCARD = '25'
    if( process_id == master_id )then
      call SEEK('C25',0)
      MMAX = 3 + NDYM + NTOX
      do NS = 1,NQSIJ
        read(1,*,IOSTAT = ISO) (BCFL_GL(NS).CQSE(M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (BCFL_GL(NS).CQSE(M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    !C26*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS
    !     SED(1 TO NSED),SND(1 TO NSND)
    NCARD = '26'

    if( process_id == master_id )then
      call SEEK('C26',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do NS = 1,NQSIJ
        read(1,*,IOSTAT = ISO) (BCFL_GL(NS).CQSE(M),M = MMIN,MMAX)
        write(mpi_efdc_out_unit,1002) NCARD

        write(mpi_efdc_out_unit,*) (BCFL_GL(NS).CQSE(M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    
    do NS = 1, NQSIJ
      call Broadcast_Scalar( BCFL_GL(NS).I       , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).J       , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).QSSE    , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).NQSMUL  , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).NQSMF   , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).NQSERQ  , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).QWIDTH  , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).QFACTOR , master_id)
      call Broadcast_Scalar( BCFL_GL(NS).GRPID   , master_id)
                                                 
      call Broadcast_Array( BCFL_GL(NS).NCSERQ   , master_id)
      call Broadcast_Array( BCFL_GL(NS).CQSE     , master_id)
    enddo

  endif

  if( NQJPIJ > 0 )then
    !C27*  READ JET-PLUME SOURCE LOCATIONS AND parameterS
    NCARD = '27'
    if( process_id == master_id )then
      call SEEK('C27',0)
      do NJP = 1,NQJPIJ
        read(1,*,IOSTAT = ISO) IDUM, JET_PLM_GL(NJP).ICALJP, JET_PLM_GL(NJP).IQJP,  JET_PLM_GL(NJP).JQJP, JET_PLM_GL(NJP).KQJP,  JET_PLM_GL(NJP).NPORTJP,   &
          JET_PLM_GL(NJP).XJETL,  JET_PLM_GL(NJP).YJETL, JET_PLM_GL(NJP).ZJET, JET_PLM_GL(NJP).PHJET, JET_PLM_GL(NJP).THJET,     &
          JET_PLM_GL(NJP).DJET,   JET_PLM_GL(NJP).CFRD,  JET_PLM_GL(NJP).DJPER

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) IDUM, JET_PLM_GL(NJP).ICALJP, JET_PLM_GL(NJP).IQJP,  JET_PLM_GL(NJP).JQJP, JET_PLM_GL(NJP).KQJP,  JET_PLM_GL(NJP).NPORTJP,   &
          JET_PLM_GL(NJP).XJETL,  JET_PLM_GL(NJP).YJETL, JET_PLM_GL(NJP).ZJET, JET_PLM_GL(NJP).PHJET, JET_PLM_GL(NJP).THJET,     &
          JET_PLM_GL(NJP).DJET,   JET_PLM_GL(NJP).CFRD,  JET_PLM_GL(NJP).DJPER
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    !C28*  READ JET-PLUME SOURCE LOCATIONS AND parameterS
    NCARD = '28'
    if( process_id == master_id )then

      call SEEK('C28',0)
      do NJP = 1,NQJPIJ
        read(1,*,IOSTAT = ISO) IDUM, JET_PLM_GL(NJP).NJEL,   JET_PLM_GL(NJP).NJPMX, JET_PLM_GL(NJP).ISENT,  JET_PLM_GL(NJP).ISTJP,  NDUM, JET_PLM_GL(NJP).IOUTJP,  &
          JET_PLM_GL(NJP).NZPRJP, JET_PLM_GL(NJP).ISDJP, JET_PLM_GL(NJP).IUPCJP, JET_PLM_GL(NJP).JUPCJP, JET_PLM_GL(NJP).KUPCJP

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) IDUM, JET_PLM_GL(NJP).NJEL,   JET_PLM_GL(NJP).NJPMX, JET_PLM_GL(NJP).ISENT,  JET_PLM_GL(NJP).ISTJP,  NDUM, JET_PLM_GL(NJP).IOUTJP,  &
          JET_PLM_GL(NJP).NZPRJP, JET_PLM_GL(NJP).ISDJP, JET_PLM_GL(NJP).IUPCJP, JET_PLM_GL(NJP).JUPCJP, JET_PLM_GL(NJP).KUPCJP
        if( ISO > 0 ) GOTO 100
        if( NJP == 1 ) NUDJP = NDUM
      enddo

    endif

    !C29*  READ ADDITIONAL JET-PLUME parameterS
    NCARD = '29'
    if( process_id == master_id )then
      call SEEK('C29',0)
      do NJP = 1,NQJPIJ
        read(1,*,IOSTAT = ISO) IDUM, JET_PLM_GL(NJP).QQCJP,      JET_PLM_GL(NJP).NQSERJP,    JET_PLM_GL(NJP).NQWRSERJP,  JET_PLM_GL(NJP).NCSERJP(1), JET_PLM_GL(NJP).NCSERJP(2),   &
          JET_PLM_GL(NJP).NCSERJP(3), JET_PLM_GL(NJP).NCSERJP(4), JET_PLM_GL(NJP).NCSERJP(5), JET_PLM_GL(NJP).NCSERJP(6), JET_PLM_GL(NJP).NCSERJP(7)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) IDUM, JET_PLM_GL(NJP).QQCJP, JET_PLM_GL(NJP).NQSERJP, JET_PLM_GL(NJP).NQWRSERJP, JET_PLM_GL(NJP).NCSERJP(1:7)

        if( ISO > 0 ) GOTO 100

        if( JET_PLM_GL(NJP).ICALJP == 2 )then
          JET_PLM_GL(NJP).QWRCJP = JET_PLM_GL(NJP).QQCJP
          JET_PLM_GL(NJP).QQCJP = 0.
        else
          JET_PLM_GL(NJP).QWRCJP = 0.
        endif
      enddo
    endif

    !C30*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT
    !     JET-PLUME SOURCES SAL,TEM,DYE,SFL,TOX(1 TO NOTX)
    NCARD = '30'
    if( process_id == master_id )then
      call SEEK('C30',0)
      MMAX = 3 + NDYM + NTOX
      do NJP = 1,NQJPIJ
        read(1,*,IOSTAT = ISO) (CQSTMP(M),M = 1,MMAX)
        write(mpi_efdc_out_unit,1002) NCARD

        write(mpi_efdc_out_unit,*) (CQSTMP(M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100

        if( JET_PLM_GL(NJP).ICALJP == 1 )then
          do MS = 1,MMAX
            JET_PLM_GL(NJP).CWRCJP(MS) = 0.0
            do K = 1,KC
              JET_PLM_GL(NJP).CQCJP(K,MS) = CQSTMP(MS)
            enddo
          enddo
        else
          do MS = 1,MMAX
            JET_PLM_GL(NJP).CWRCJP(MS) = CQSTMP(MS)
            do K = 1,KC
              JET_PLM_GL(NJP).CQCJP(K,MS) = 0.0
            enddo
          enddo
        endif
      enddo

    endif

    !C31*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT
    !     JET-PLUME SOURCES SED(1 TO NSED),SND(1 TO NSND)
    NCARD = '31'
    if( process_id == master_id )then
      call SEEK('C31',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      if( ISTRAN(8) > 0 )then
        MMAX = MMAX + NWQV
      endif
      do NJP = 1,NQJPIJ
        read(1,*,IOSTAT = ISO) (CQSTMP(M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CQSTMP(M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100

        MS = MAX(MMIN-1, 0)
        if( JET_PLM_GL(NJP).ICALJP == 1 )then
          do MS = MMIN,MMAX
            JET_PLM_GL(NJP).CWRCJP(MS) = 0.0
            do K = 1,KC
              JET_PLM_GL(NJP).CQCJP(K,MS) = CQSTMP(MS)
            enddo
          enddo

          do M = 1,NWQV
            if( ISKINETICS(M) > 0 )then
              MS = MS + 1
              do K = 1,KC
                JET_PLM_GL(NJP).CQCJP(K,MS) = CQSTMP(M + MMIN - 1)        ! *** Only used constituents that will be simulated
              enddo
            endif
          enddo
        else
          do MS = MMIN,MMAX
            JET_PLM_GL(NJP).CWRCJP(MS) = CQSTMP(MS)
            do K = 1,KC
              JET_PLM_GL(NJP).CQCJP(K,MS) = 0.
            enddo
          enddo
        endif
      enddo

    endif

    call Broadcast_Scalar( NUDJP,     master_id)
    do NJP = 1,NQJPIJ
      call Broadcast_Scalar( JET_PLM_GL(NJP).ICALJP,    master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).IQJP,      master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).JQJP,      master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).KQJP,      master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).IUPCJP,    master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).JUPCJP,    master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).KUPCJP,    master_id)

      call Broadcast_Scalar( JET_PLM_GL(NJP).CFRD,      master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).DJET,      master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).DJPER,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).PHJET,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).QQCJP,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).QWRCJP,    master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).THJET,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).XJETL,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).YJETL,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).ZJET,      master_id)

      call Broadcast_Scalar( JET_PLM_GL(NJP).IOUTJP,    master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).ISENT,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).ISTJP,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).ISDJP,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).NJEL,      master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).NJPMX,     master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).NPORTJP,   master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).NQSERJP,   master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).NQWRSERJP, master_id)
      call Broadcast_Scalar( JET_PLM_GL(NJP).NZPRJP,    master_id)

      call Broadcast_Array( JET_PLM_GL(NJP).NCSERJP,    master_id)
      call Broadcast_Array( JET_PLM_GL(NJP).CQCJP,      master_id)
      call Broadcast_Array( JET_PLM_GL(NJP).CWRCJP,     master_id)
    enddo

  endif

  if( NQCTL > 0 )then
    !C32*  READ SURF ELEV OR PRESS DEPENDENT FLOW CONTROL STRUCTURE INFO
    NCARD = '32'

    if( process_id == master_id )then
      call SEEK('C32',0)
      do NC = 1,NQCTL
        read(1,*,IOSTAT = ISO) HYD_STR_GL(NC).IQCTLU, HYD_STR_GL(NC).JQCTLU, HYD_STR_GL(NC).IQCTLD, HYD_STR_GL(NC).JQCTLD,  HYD_STR_GL(NC).NQCTYP, HYD_STR_GL(NC).NQCTLQ,  HYD_STR_GL(NC).NQCMUL,    &
          HYD_STR_GL(NC).HQCTLU, HYD_STR_GL(NC).HQCTLD, HYD_STR_GL(NC).QCTLMU, HYD_STR_GL(NC).BQCLCE,  HYD_STR_GL(NC).NQCMINS, HYD_STR_GL(NC).HS_FACTOR, &
          HYD_STR_GL(NC).MMASKS, HYD_STR_GL(NC).TRANSIT, HYD_STR_GL(NC).QCTLGRP

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) HYD_STR_GL(NC).IQCTLU, HYD_STR_GL(NC).JQCTLU, HYD_STR_GL(NC).IQCTLD, HYD_STR_GL(NC).JQCTLD,  HYD_STR_GL(NC).NQCTYP, HYD_STR_GL(NC).NQCTLQ,  HYD_STR_GL(NC).NQCMUL,    &
          HYD_STR_GL(NC).HQCTLU, HYD_STR_GL(NC).HQCTLD, HYD_STR_GL(NC).QCTLMU, HYD_STR_GL(NC).BQCLCE,  HYD_STR_GL(NC).NQCMINS, HYD_STR_GL(NC).HS_FACTOR, &
          HYD_STR_GL(NC).MMASKS, HYD_STR_GL(NC).TRANSIT, HYD_STR_GL(NC).QCTLGRP
        if( ISO > 0 ) GOTO 100

        if( HYD_STR_GL(NC).NQCTYP > 4 )then
          if( HYD_STR_GL(NC).HS_FACTOR <= 0  )then
            PRINT *,' *** BAD HS_FACTOR: NC, HS_FACTOR = ',NC, HYD_STR_GL(NC).HS_FACTOR
            call STOPP('')
          endif
        endif

      enddo
    endif

    do NC = 1,NQCTL
      call Broadcast_Scalar(HYD_STR_GL(NC).IQCTLU     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).JQCTLU     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).IQCTLD     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).JQCTLD     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).NQCTYP     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).NQCTLQ     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).NQCMUL     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).HQCTLU     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).HQCTLD     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).QCTLMU     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).BQCLCE     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).NQCMINS    , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).HS_FACTOR  , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).MMASKS     , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).TRANSIT    , master_id)
      call Broadcast_Scalar(HYD_STR_GL(NC).QCTLGRP    , master_id)
    enddo

    do NC = 1,NQCTL
      do K = 1,KC
        QCTLTO(K,NC) = 0.
        QCTLT(K,NC,1) = 0.
        QCTLT(K,NC,2) = 0.
      enddo
    enddo

    !C32A*  READ THE EQUATION parameterS FOR EACH OF THE HYDRAULIC STRUCTURE EQUATIONS
    if( NHYDST > 0 )then
      NCARD = '32A'
      if( process_id == master_id )then
        call SEEK('C32A',0)
        do L = 1,NHYDST
          ! *** COMPUTE FLOWS USING HYDRAULIC STRUCTURE EQUATIONS
          read(1,*,IOSTAT = ISO) NS, NX, HS_REVERSE(L), HS_XSTYPE(L), HS_WIDTH(L), HS_HEIGHT(L), HS_LENGTH(L), HS_MANN(L), HS_ANGLE(L), HS_USELEV(L), HS_DSELEV(L),  &
            HS_COEFF(L,1), HS_COEFF(L,2), HS_COEFF(L,3), HS_COEFF(L,4)

          call HYDSTRUCT_CHECK(NX,L)  ! *** CHECK STRUCTURE DEFINITIONS

          ! *** Use THE SAME INVERT ELEVATION FOR US/DS FOR SLUICE GATES AND WEIRS
          if( NX > 5 )then
            HS_DSELEV(L) = HS_USELEV(L)
          endif

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) NS, NX, HS_REVERSE(L), HS_XSTYPE(L), HS_WIDTH(L), HS_HEIGHT(L), HS_LENGTH(L), HS_MANN(L), HS_ANGLE(L), HS_USELEV(L), HS_DSELEV(L), &
            HS_COEFF(L,1), HS_COEFF(L,2), HS_COEFF(L,3), HS_COEFF(L,4)
          HS_ANGLE(L) = HS_ANGLE(L)*PI/180.0D0
        enddo
      endif

      call Broadcast_Scalar(NS, master_id)
      call Broadcast_Scalar(NX, master_id)
      call Broadcast_Array(HS_REVERSE, master_id)
      call Broadcast_Array(HS_XSTYPE , master_id)
      call Broadcast_Array(HS_WIDTH  , master_id)
      call Broadcast_Array(HS_HEIGHT , master_id)
      call Broadcast_Array(HS_LENGTH , master_id)
      call Broadcast_Array(HS_MANN   , master_id)
      call Broadcast_Array(HS_ANGLE  , master_id)
      call Broadcast_Array(HS_USELEV , master_id)
      call Broadcast_Array(HS_DSELEV , master_id)
      call Broadcast_Array(HS_COEFF  , master_id)
      call Broadcast_Array(HS_COEFF  , master_id)
      call Broadcast_Array(HS_COEFF  , master_id)
      call Broadcast_Array(HS_COEFF  , master_id)
    endif

    !C32B*  READ THE CONTROL INFO FOR ALL HYDRAULIC STRUCTURES
    NCARD = '32B'
    if( process_id == master_id )then
      if( NQCTLSER > 0 .or. NQCTRULES > 0 )then
        call SEEK('C32B',0)
        do NC = 1,NQCTL
          HSCTL_GL(NC).ITYPE = 0
          HSCTL_GL(NC).ID = 0
          HSCTL_GL(NC).SUBID = 0

          read(1,*,IOSTAT = ISO) HSCTL_GL(NC).ITYPE, HSCTL_GL(NC).ID, ITMPU, JTMPU, ITMPD, JTMPD, &
            HSCTL_GL(NC).CUR.STATE, HSCTL_GL(NC).CUR.HEIGHT, HSCTL_GL(NC).CUR.WIDTH, &
            HSCTL_GL(NC).CUR.SILL,  HYD_STR_GL(NC).NWRGRP

          HSCTL_GL(NC).IREFUP = ITMPU
          HSCTL_GL(NC).JREFUP = JTMPU
          HSCTL_GL(NC).IREFDN = ITMPD
          HSCTL_GL(NC).JREFDN = JTMPD

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*) HSCTL_GL(NC).ITYPE, HSCTL_GL(NC).ID, &
            HSCTL_GL(NC).IREFUP, HSCTL_GL(NC).JREFUP, HSCTL_GL(NC).IREFDN, HSCTL_GL(NC).JREFDN, &
            HSCTL_GL(NC).CUR.STATE, HSCTL_GL(NC).CUR.HEIGHT, HSCTL_GL(NC).CUR.WIDTH, &
            HSCTL_GL(NC).CUR.SILL,  HYD_STR_GL(NC).NWRGRP
        enddo
      else
        do NC = 1,NQCTL
          HSCTL_GL(NC).ITYPE = 0
          HSCTL_GL(NC).ID = 0

          HSCTL_GL(NC).IREFUP = 0
          HSCTL_GL(NC).JREFUP = 0
          HSCTL_GL(NC).IREFDN = 0
          HSCTL_GL(NC).JREFDN = 0

          HSCTL_GL(NC).CUR.STATE = 0
          HSCTL_GL(NC).CUR.HEIGHT = 0.
          HSCTL_GL(NC).CUR.WIDTH = 0.
          HSCTL_GL(NC).CUR.SILL = 0.
          HSCTL_GL(NC).CUR.FLOW = 0.
        enddo
      endif
    endif

    Do L = 1, NQCTL
      call Broadcast_Scalar(HSCTL_GL(L).ITYPE     , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).ID        , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).IREFUP    , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).JREFUP    , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).IREFDN    , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).JREFDN    , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).CUR.STATE , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).CUR.HEIGHT, master_id)
      call Broadcast_Scalar(HSCTL_GL(L).CUR.WIDTH , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).CUR.SILL  , master_id)
      call Broadcast_Scalar(HSCTL_GL(L).CUR.FLOW  , master_id)
      call Broadcast_Scalar(HYD_STR_GL(L).NWRGRP   , master_id)
    enddo

  endif

  if( NQWR > 0 )then
    !C33*  READ FLOW WITHDRAWAL, HEAT OR MATERIAL ADDITION, FLOW RETURN DATA
    NCARD = '33'
    if( process_id == master_id )then
      call SEEK('C33',0)
      !J = 0
      do NWR = 1,NQWR
        read(1,*,IOSTAT = ISO) WITH_RET_GL(NWR).IQWRU,    WITH_RET_GL(NWR).JQWRU,   WITH_RET_GL(NWR).KQWRU,   WITH_RET_GL(NWR).IQWRD,      &
          WITH_RET_GL(NWR).JQWRD,    WITH_RET_GL(NWR).KQWRD,   WITH_RET_GL(NWR).QWR,     WITH_RET_GL(NWR).NQWRSERQ,   &
          WITH_RET_GL(NWR).NQWRMFU,  WITH_RET_GL(NWR).NQWRMFD, WITH_RET_GL(NWR).BQWRMFU, WITH_RET_GL(NWR).BQWRMFD,    &
          WITH_RET_GL(NWR).ANGWRMFD, WITH_RET_GL(NWR).GROUPID

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) WITH_RET_GL(NWR).IQWRU,   WITH_RET_GL(NWR).JQWRU,    WITH_RET_GL(NWR).KQWRU,    WITH_RET_GL(NWR).IQWRD,   WITH_RET_GL(NWR).JQWRD,    &
          WITH_RET_GL(NWR).KQWRD,   WITH_RET_GL(NWR).QWR,      WITH_RET_GL(NWR).NQWRSERQ, WITH_RET_GL(NWR).NQWRMFU, WITH_RET_GL(NWR).NQWRMFD,  &
          WITH_RET_GL(NWR).BQWRMFU, WITH_RET_GL(NWR).BQWRMFD,  WITH_RET_GL(NWR).ANGWRMFD
        if( ISO > 0 ) GOTO 100

        !J = MAX(J,WITH_RET_GL(NWR).GROUPID)
        !IF( WITH_RET_GL(NWR).GROUPID < 1 )then
        !  J = J + 1
        !  WITH_RET_GL(NWR).GROUPID = J
        !ENDIF
      enddo

      !C33A*  READ FLOW WITHDRAWAL/RETURN CONTROL
      NCARD = '33A'
      call SEEK('C33A',0)
      do NWR = 1,NQWR
        WITH_RET_CTL_GL(NWR).ITYPE = 0
        WITH_RET_CTL_GL(NWR).ID = WITH_RET_GL(NWR).NQWRSERQ
        WITH_RET_CTL_GL(NWR).SUBID = 0

        read(1,*,IOSTAT = ISO) WITH_RET_CTL_GL(NWR).ITYPE, WITH_RET_CTL_GL(NWR).ID, ITMPU, JTMPU, ITMPD, JTMPD, &
          WITH_RET_CTL_GL(NWR).CUR.STATE, WITH_RET_CTL_GL(NWR).CUR.FLOW

        WITH_RET_CTL_GL(NWR).IREFUP = ITMPU
        WITH_RET_CTL_GL(NWR).JREFUP = JTMPU
        WITH_RET_CTL_GL(NWR).IREFDN = ITMPD
        WITH_RET_CTL_GL(NWR).JREFDN = JTMPD

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) WITH_RET_CTL_GL(NWR).ITYPE,  WITH_RET_CTL_GL(NWR).ID,     WITH_RET_CTL_GL(NWR).IREFUP,    WITH_RET_CTL_GL(NWR).JREFUP,   &
          WITH_RET_CTL_GL(NWR).IREFDN, WITH_RET_CTL_GL(NWR).JREFDN, WITH_RET_CTL_GL(NWR).CUR.STATE, WITH_RET_CTL_GL(NWR).CUR.FLOW
      enddo

    endif

    do NWR = 1,NQWR
      call Broadcast_Scalar(WITH_RET_GL(NWR).IQWRU    , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).JQWRU    , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).KQWRU    , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).IQWRD    , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).JQWRD    , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).KQWRD    , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).QWR      , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).NQWRSERQ , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).NQWRMFU  , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).NQWRMFD  , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).BQWRMFU  , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).BQWRMFD  , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).ANGWRMFD , master_id)
      call Broadcast_Scalar(WITH_RET_GL(NWR).GROUPID  , master_id)

      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).ITYPE     , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).ID        , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).IREFUP    , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).JREFUP    , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).IREFDN    , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).JREFDN    , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).CUR.STATE , master_id)
      call Broadcast_Scalar(WITH_RET_CTL_GL(NWR).CUR.FLOW  , master_id)
    enddo

    !C34*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '34'
    if( process_id == master_id )then
      call SEEK('C34',0)
      MMAX = 3 + NDYM + NTOX
      do NWR = 1,NQWR
        read(1,*,IOSTAT = ISO) (CQWR(NWR,MS),MS = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CQWR(NWR,MS),MS = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    !C35*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES
    !     SED(1 TO NSED),SND(1 TO NSND)
    NCARD = '35'
    if( process_id == master_id )then
      call SEEK('C35',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do NWR = 1,NQWR
        read(1,*,IOSTAT = ISO) (CQWR(NWR,MS),MS = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CQWR(NWR,MS),MS = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

  endif

  call Broadcast_Array(CQWR, master_id)

  if( NSED > 0 .or. NSND > 0 )then
    !C36*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
    NCARD = '36'
    if( process_id == master_id )then
      call SEEK('C36',0)
      read(1,*,IOSTAT = ISO)ISEDINT,ISEDBINT,NSEDFLUME,ISMUD,ISBEDMAP,ISEDVW,ISNDVW,KB,ISDTXBUG

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*)ISEDINT,ISEDBINT,NSEDFLUME,ISMUD,ISBEDMAP,ISEDVW,ISNDVW,KB,ISDTXBUG
      if( ISO > 0 ) GOTO 100
    endif

    call Broadcast_Scalar(ISEDINT,   master_id)
    call Broadcast_Scalar(ISEDBINT,  master_id)
    call Broadcast_Scalar(NSEDFLUME, master_id)
    call Broadcast_Scalar(ISMUD    , master_id)
    call Broadcast_Scalar(ISBEDMAP , master_id)
    call Broadcast_Scalar(ISEDVW   , master_id)
    call Broadcast_Scalar(ISNDVW   , master_id)
    call Broadcast_Scalar(KB       , master_id)
    call Broadcast_Scalar(ISDTXBUG , master_id)

    ! *** FORCE READ IF SPATIALLY VARYING SEDIMENTS
    if( NSEDFLUME == 3 .and. ISEDINT < 2 ) ISEDINT = 3

    !C36A*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
    NCARD = '36A'
    if( process_id == master_id )then
      call SEEK('C36A',0)
      read(1,*,IOSTAT = ISO)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO
      if( ISO > 0 ) GOTO 100
    endif

    call Broadcast_Scalar(ISBEDSTR , master_id)
    call Broadcast_Scalar(ISBSDIAM , master_id)
    call Broadcast_Scalar(ISBSDFUF , master_id)
    call Broadcast_Scalar(COEFTSBL , master_id)
    call Broadcast_Scalar(VISMUDST , master_id)
    call Broadcast_Scalar(ISBKERO  , master_id)

    !C36B*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
    NCARD = '36B'
    if( process_id == master_id )then
      call SEEK('C36B',0)
      if( NSED>0 .or. NSND > 0 )then
        read(1,*,IOSTAT = ISO) ldum,ISNDAL,IALTYP,IALSTUP,ISEDEFF,HBEDAL,COEHEFF,COEHEFF2

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) ldum,ISNDAL,IALTYP,IALSTUP,ISEDEFF,HBEDAL,COEHEFF,COEHEFF2
        if( ISO > 0 ) GOTO 100
      endif
      ! *** FORCE BED ARMORING OPTION "INITIALIZATION AT STARTUP" TO BE OFF
      if( ISRESTI > 0 ) IALSTUP = 0
    endif

    call Broadcast_Scalar(ISNDAL   , master_id)
    call Broadcast_Scalar(IALTYP   , master_id)
    call Broadcast_Scalar(IALSTUP  , master_id)
    call Broadcast_Scalar(ISEDEFF  , master_id)
    call Broadcast_Scalar(HBEDAL   , master_id)
    call Broadcast_Scalar(COEHEFF  , master_id)
    call Broadcast_Scalar(COEHEFF2 , master_id)

    !C37*  BED MECHANICAL PROPERTIES parameter SET 1
    NCARD = '37'
    if( process_id == master_id )then
      call SEEK('C37',0)
      if( NSED > 0 .or. NSND > 0 )then
        read(1,*,IOSTAT = ISO) SEDSTEP, SEDSTART, IBMECH, IMORPH, HBEDMAX, BEDPORC, SEDMDMX, SEDMDMN, SEDVDRD, SEDVDRM, SEDVRDT

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) SEDSTEP, SEDSTART, IBMECH, IMORPH, HBEDMAX, BEDPORC, SEDMDMX, SEDMDMN, SEDVDRD, SEDVDRM, SEDVRDT
        if( ISO > 0 ) GOTO 100

        if( SEDSTART <= TBEGIN ) SEDSTART = TBEGIN
      else
        BEDPORC = 0.4
        SEDVDRD = BEDPORC/(1.-BEDPORC)
        SEDVDRM = SEDVDRD
        SEDVDRT = 0.0
      endif
    endif

    call Broadcast_Scalar(SEDSTEP  , master_id)
    call Broadcast_Scalar(SEDSTART , master_id)
    call Broadcast_Scalar(IBMECH   , master_id)
    call Broadcast_Scalar(IMORPH   , master_id)
    call Broadcast_Scalar(HBEDMAX  , master_id)
    call Broadcast_Scalar(BEDPORC  , master_id)
    call Broadcast_Scalar(SEDMDMX  , master_id)
    call Broadcast_Scalar(SEDMDMN  , master_id)
    call Broadcast_Scalar(SEDVDRD  , master_id)
    call Broadcast_Scalar(SEDVDRM  , master_id)
    call Broadcast_Scalar(SEDVRDT  , master_id)

    ! *** Uniform porosities and void ratios
    if( IBMECH == 0 )then
      SEDVDRD  = BEDPORC/(1.-BEDPORC)
      SEDVDRM  = SEDVDRD
      SEDVRDT  = 0.0
      VDRBED   = SEDVDRD
      VDRBED1  = SEDVDRD
      PORBED   = BEDPORC
      PORBED1  = BEDPORC
      BEDBINIT = SEDVDRD
    endif

    SNDVDRD = BEDPORC/(1.-BEDPORC)
    do NS = 1,NSED
      VDRDEPO(NS) = SEDVDRD
    enddo
    do NS = 1,NSND
      NX = NS + NSED
      VDRDEPO(NX) = SNDVDRD
    enddo
  endif

  if( NSED > 0 )then
    !C39*  READ COHESIVE SEDIMENT parameter SET 1 REPEAT DATA LINE NSED TIMES
    NCARD = '39'

    if( process_id == master_id )then
      call SEEK('C39',0)
      if( NSED > 0 )then
        do NS = 1,NSED
          read(1,*,IOSTAT = ISO)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS),ISPROBDEP(NS)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISPROBDEP(NS)
          if( ISO > 0 ) GOTO 100
          SEDDIA(NS) = 0.
        enddo
      endif
    endif

    call Broadcast_Array(SEDO      , master_id)
    call Broadcast_Array(SEDBO     , master_id)
    call Broadcast_Array(SDEN      , master_id)
    call Broadcast_Array(SSG       , master_id)
    call Broadcast_Array(WSEDO     , master_id)
    call Broadcast_Array(SEDN      , master_id)
    call Broadcast_Array(SEXP      , master_id)
    call Broadcast_Array(TAUD      , master_id)
    call Broadcast_Array(ISPROBDEP , master_id)

    !C40*  READ COHESIVE SEDIMENT parameter SET 2 REPEAT DATA LINE NSED TIMES
    NCARD = '40'
    if( process_id == master_id )then
      call SEEK('C40',0)
      do NS = 1,NSED
        read(1,*,IOSTAT = ISO) IWRSP(NS), IWRSPB(NS), WRSPO(NS), TAUR(NS), TAUN(NS), TEXP(NS), VDRRSPO(NS), COSEDHID(NS), SEDS(NS).WCLIMIT

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) IWRSP(NS), IWRSPB(NS), WRSPO(NS), TAUR(NS), TAUN(NS), TEXP(NS), VDRRSPO(NS), COSEDHID(NS), SEDS(NS).WCLIMIT
        if( ISO > 0 ) GOTO 100

        if( NS == 1 .and. IWRSP(NS) == 999 )then
          write(*,'(A)')'READING TAU_CRIT_COH.INP'
          open(1001,FILE = 'tau_crit_coh.inp',STATUS = 'OLD')
          do L = 2, 4393
            read(1001,*,IOSTAT = ISO) (TAUCRCOH(L,K),K = 1,10)
          enddo
          close(1001)
        endif
        if( ISO > 0 ) GOTO 100
        ISNDEQ(NS) = 0

      enddo

    endif

    call Broadcast_Array(IWRSP    , master_id)
    call Broadcast_Array(IWRSPB   , master_id)
    call Broadcast_Array(WRSPO    , master_id)
    call Broadcast_Array(TAUR     , master_id)
    call Broadcast_Array(TAUN     , master_id)
    call Broadcast_Array(TEXP     , master_id)
    call Broadcast_Array(VDRRSPO  , master_id)
    call Broadcast_Array(COSEDHID , master_id)

  endif

  TAUCMIN = 1000.
  ICALC_BL = 0
  SSG = 2.65         ! *** DEFAULT GRAIN DENSITY USING QUARTZ
  if( NSND > 0 )then
    !C41*  READ NONCOHESIVE SEDIMENT parameter SET 1 REPEAT DATA LINE NSND TIMES
    NCARD = '41'

    if( process_id == master_id )then
      call SEEK('C41',0)
      do NX = 1,NSND
        NS = NSED + NX
        read(1,*,IOSTAT = ISO) SEDO(NS), SEDBO(NS), SDEN(NS), SSG(NS), SEDDIA(NS), WSEDO(NS), SEDN(NS), SEXP(NS), TAUD(NS), ISEDSCOR(NS), SEDS(NS).WCLIMIT

        ! *** IF SETTLING VELOCITY IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
        if( WSEDO(NS) < 0.0 )then
          WSEDO(NS) = SETSTVEL(SEDDIA(NS),SSG(NS))
        endif

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) SEDO(NS), SEDBO(NS), SDEN(NS), SSG(NS), SEDDIA(NS), WSEDO(NS), SEDN(NS), SEXP(NS), TAUD(NS), ISEDSCOR(NS), SEDS(NS).WCLIMIT
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    call Broadcast_Array(SEDO     , master_id)
    call Broadcast_Array(SEDBO    , master_id)
    call Broadcast_Array(SDEN     , master_id)
    call Broadcast_Array(SSG      , master_id)
    call Broadcast_Array(SEDDIA   , master_id)
    call Broadcast_Array(WSEDO    , master_id)
    call Broadcast_Array(SEDN     , master_id)
    call Broadcast_Array(SEXP     , master_id)
    call Broadcast_Array(TAUD     , master_id)
    call Broadcast_Array(ISEDSCOR , master_id)
    do NS = 1,NSEDS
      call Broadcast_Scalar(SEDS(NS).WCLIMIT , master_id)
    enddo

    !C42*  READ NONCOHESIVE SEDIMENT parameter SET 2 REPEAT DATA LINE NSND TIMES
    NCARD = '42'
    if( process_id == master_id )then
      call SEEK('C42',0)
      do NX = 1,NSND
        NS = NSED + NX
        read(1,*,IOSTAT = ISO)ISNDEQ(NS),ISBDLD(NS),TAUR(NS),TAUN(NS),TCSHIELDS(NS),ISLTAUC(NS),IBLTAUC(NS),IROUSE(NX),ISNDM1(NX),ISNDM2(NX),RSNDM(NX)

        ! IF TAUR(NS) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
        !       TAUR:     CRITICAL SHIELDS STRESS IN (M/S)**2   (ISNDEQ = 2)
        !       TAUN:     EQUAL TO TAUR FOR NONCHOESIVE SED TRANS  (ISNDEQ = 2)
        !       TEXP:     CRITICAL SHIELDS parameter  (ISNDEQ = 2)
        DSTR = 0.0
        USTR = 0.0

        ! *** IF TAUR(NS) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
        if( TAUR(NS) < 0.0 )then
          call SETSHLD(TAUR(NS),TCSHIELDS(NS),SEDDIA(NS),SSG(NS),DSTR,USTR)
          TAUN(NS) = TAUR(NS)
        endif
        TAUCMIN = MIN(TAUCMIN,TAUR(NS))

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*)ISNDEQ(NS),TAUR(NS),TAUN(NS),TCSHIELDS(NS),SEDDIA(NS),SSG(NS),DSTR,USTR
        if( ISO > 0 ) GOTO 100

        IWRSP(NS) = 0
        WRSPO(NS) = 0
      enddo

    endif

    call Broadcast_Array(ISNDEQ         , master_id)
    call Broadcast_Array(ISBDLD         , master_id)
    call Broadcast_Array(TAUR           , master_id)
    call Broadcast_Array(TAUN           , master_id)
    call Broadcast_Array(TCSHIELDS      , master_id)
    call Broadcast_Array(SEDDIA         , master_id)
    call Broadcast_Array(ISLTAUC        , master_id)
    call Broadcast_Array(IBLTAUC        , master_id)
    call Broadcast_Array(IROUSE         , master_id)
    call Broadcast_Array(ISNDM1         , master_id)
    call Broadcast_Array(ISNDM2         , master_id)
    call Broadcast_Array(RSNDM          , master_id)
    call Broadcast_Array(SSG            , master_id)
    call Broadcast_Scalar(DSTR          , master_id)
    call Broadcast_Scalar(USTR          , master_id)
    call Broadcast_Scalar(TAUCMIN       , master_id)

    !C42A*  READ NONCOHESIVE SEDIMENT BED LOAD parameterS
    NCARD = '42A'
    if( process_id == master_id )then
      call SEEK('C42A',0)
      do NX = 1,NSND
        read(1,*,IOSTAT = ISO)ICALC_BL,SBDLDA(NX),SBDLDB(NX),SBDLDG1(NX),SBDLDG2(NX),SBDLDG3(NX),SBDLDG4(NX),SBDLDP(NX),ISBLFUC,BLBSNT

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*)ICALC_BL,SBDLDA(NX),SBDLDB(NX),SBDLDG1(NX),SBDLDG2(NX),SBDLDG3(NX),SBDLDG4(NX),SBDLDP(NX),ISBLFUC,BLBSNT
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    call Broadcast_Scalar(ICALC_BL , master_id)
    call Broadcast_Scalar(ISBLFUC  , master_id)
    call Broadcast_Scalar(BLBSNT   , master_id)
    call Broadcast_Array(SBDLDA    , master_id)
    call Broadcast_Array(SBDLDB    , master_id)
    call Broadcast_Array(SBDLDG1   , master_id)
    call Broadcast_Array(SBDLDG2   , master_id)
    call Broadcast_Array(SBDLDG3   , master_id)
    call Broadcast_Array(SBDLDG4   , master_id)
    call Broadcast_Array(SBDLDP    , master_id)

  endif

  if( NTOX > 0 )then
    !C43A*  READ TOXIC CONTAMINANT INITIAL CONDITIONS
    NCARD = '43A'
    if( process_id == master_id )then
      call SEEK('C43A',0)
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO) NDUM, ITXINT(NT), ITXBDUT(NT), TOXINTW(NT), TOXINTB(NT)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) NDUM, ITXINT(NT), ITXBDUT(NT), TOXINTW(NT), TOXINTB(NT)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    call Broadcast_Array(ITXINT , master_id)
    call Broadcast_Array(ITXBDUT, master_id)
    call Broadcast_Array(TOXINTW, master_id)
    call Broadcast_Array(TOXINTB, master_id)

    !C43B*  READ TOXIC KINETIC OPTION FLAGS
    NCARD = '43B'
    if( process_id == master_id )then
      call SEEK('C43B',0)
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO) NDUM, ITOXKIN(1,NT) ,ITOXKIN(2,NT), ITOXKIN(3,NT), ITOXKIN(4,NT), ITOXKIN(5,NT), ITOXKIN(6,NT), TOXS(NT).WCLIMIT

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) NDUM, ITOXKIN(1,NT) ,ITOXKIN(2,NT), ITOXKIN(3,NT), ITOXKIN(4,NT), ITOXKIN(5,NT), ITOXKIN(6,NT), TOXS(NT).WCLIMIT
        if( ISO > 0 ) GOTO 100
      enddo

    endif

    call Broadcast_Array(ITOXKIN, master_id)
    do NT = 1,NTOX
      call Broadcast_Scalar(TOXS(NT).WCLIMIT, master_id)
    enddo

    !C43C*  READ TOXIC TIMING AND VOLATILIZATION SWITCHES
    NCARD = '43C'
    if( process_id == master_id )then
      call SEEK('C43C',0)
      read(1,*,IOSTAT = ISO) TOXSTEPW, TOXSTEPB
      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) TOXSTEPW, TOXSTEPB
      if( ISO > 0 ) GOTO 100
      if( TOXSTEPW < SEDSTEP ) TOXSTEPW = SEDSTEP    ! *** TOXIC KINETICS CANNOT OPERATE ON A FINER TIMESCALE THAN SEDIMENTS
      if( TOXSTEPB < SEDSTEP ) TOXSTEPB = SEDSTEP    ! *** TOXIC BED MIXING AND DIFFUSION CANNOT OPERATE ON A FINER TIMESCALE THAN SEDIMENTS
      if( ISTRAN(5) > 0 .and. ISTRAN(2) == 0 .and. ( IVOLTEMP < 0 .or. IVOLTEMP-1 > NCSER(2) ) ) CALL STOPP('IVOLTEMP IS OUT OF RANGE')
    endif

    call Broadcast_Scalar(TOXSTEPW      , master_id)
    call Broadcast_Scalar(TOXSTEPB      , master_id)

    !C43D*  READ TOXIC BULK DECAY AND BIODEGRADATION parameterS
    NCARD = '43D'
    if( process_id == master_id )then
      call SEEK('C43D',0)
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO) NDUM,TOXS(NT).BLK_KW,   TOXS(NT).BLK_KB,   TOXS(NT).BLK_MXD, TOXS(NT).BIO_KW, TOXS(NT).BIO_KB, TOXS(NT).BIO_MXD, &
          TOXS(NT).BIO_Q10W, TOXS(NT).BIO_Q10B, TOXS(NT).BIO_TW,   TOXS(NT).BIO_TB
        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*)  NDUM,TOXS(NT).BLK_KW,   TOXS(NT).BLK_KB,   TOXS(NT).BLK_MXD, TOXS(NT).BIO_KW, TOXS(NT).BIO_KB, TOXS(NT).BIO_MXD, &
          TOXS(NT).BIO_Q10W, TOXS(NT).BIO_Q10B, TOXS(NT).BIO_TW,  TOXS(NT).BIO_TB
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    do NT = 1,NTOX
      call Broadcast_Scalar( TOXS(NT).BLK_KW    , master_id)
      call Broadcast_Scalar( TOXS(NT).BLK_KB    , master_id)
      call Broadcast_Scalar( TOXS(NT).BLK_MXD   , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_KW    , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_KB    , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_MXD   , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_Q10W  , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_Q10B  , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_TW    , master_id)
      call Broadcast_Scalar( TOXS(NT).BIO_TB    , master_id)
    enddo

    !C43E*  READ TOXIC VOLATILIZATION parameterS
    NCARD = '43E'
    if( process_id == master_id )then
      call SEEK('C43E',0)
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO) NDUM, TOXS(NT).VOL.KL_OPT, TOXS(NT).VOL.MW, TOXS(NT).VOL.HE, TOXS(NT).VOL.TCOEFF, TOXS(NT).VOL.AIRCON, TOXS(NT).VOL.MULT
        write(mpi_efdc_out_unit,1002) NCARD

        write(mpi_efdc_out_unit,*) NDUM, TOXS(NT).VOL.KL_OPT, TOXS(NT).VOL.MW, TOXS(NT).VOL.HE, TOXS(NT).VOL.TCOEFF, TOXS(NT).VOL.AIRCON, TOXS(NT).VOL.MULT
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    do NT = 1,NTOX
      call Broadcast_Scalar(TOXS(NT).VOL.KL_OPT , master_id)
      call Broadcast_Scalar(TOXS(NT).VOL.MW     , master_id)
      call Broadcast_Scalar(TOXS(NT).VOL.HE     , master_id)
      call Broadcast_Scalar(TOXS(NT).VOL.TCOEFF , master_id)
      call Broadcast_Scalar(TOXS(NT).VOL.AIRCON , master_id)
      call Broadcast_Scalar(TOXS(NT).VOL.MULT   , master_id)
    enddo

    ! *** Delme - Photolysis   RKTOXP(NT),SKTOXP(NT)

    !C44*  READ TOXIC CONTAMINANT parameterS: SORBTION
    NCARD = '44'
    if( process_id == master_id )then
      call SEEK('C44',0)
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO) NDUM, ISTOC(NT), DIFTOX(NT), DIFTOXS(NT), PDIFTOX(NT), DPDIFTOX(NT)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,* )NDUM, ISTOC(NT), DIFTOX(NT), DIFTOXS(NT), PDIFTOX(NT), DPDIFTOX(NT)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    call Broadcast_Array(ISTOC    , master_id)
    call Broadcast_Array(DIFTOX   , master_id)
    call Broadcast_Array(DIFTOXS  , master_id)
    call Broadcast_Array(PDIFTOX  , master_id)
    call Broadcast_Array(DPDIFTOX , master_id)

    ! *** Set up spatially varying particle mixing flags
    do NT = 1,NTOX
      ISPMXZ(NT) = 0
      if( PDIFTOX(NT) < 0.0 ) ISPMXZ(NT) = 1

      ISDIFBW(NT) = 0
      if( DIFTOXS(NT) < 0.0 )then
        DIFTOXS(NT) = ABS(DIFTOXS(NT))
        ISDIFBW(NT) = 1
      endif
    enddo

    !C45*  READ TOXIC CONTAMINANT-SEDIMENT INTERACTION parameterS
    NCARD = '45'
    if( process_id == master_id )then
      call SEEK('C45',0)
      do NT = 1,NTOX
        if( NSED > 0 )then
          do NS = 1,NSED
            read(1,*,IOSTAT = ISO) NDUM1, NDUM2, ITXPARW(NS,NT), CONPARW(NS,NT), TOXPARW(1,NS,NT), TOXPARB(1,NS,NT)

            write(mpi_efdc_out_unit,1002) NCARD
            write(mpi_efdc_out_unit,*) NDUM1, NDUM2, ITXPARW(NS,NT), CONPARW(NS,NT), TOXPARW(1,NS,NT), TOXPARB(1,NS,NT)
            if( ISO > 0 ) GOTO 100
          enddo
        endif
        if( NSND > 0 )then
          do NX = 1,NSND
            NS = NSED + NX
            read(1,*,IOSTAT = ISO) NDUM1, NDUM2, ITXPARW(NS,NT), CONPARW(NS,NT), TOXPARW(1,NS,NT), TOXPARB(1,NS,NT)

            write(mpi_efdc_out_unit,1002) NCARD
            write(mpi_efdc_out_unit,*) NDUM1, NDUM2, ITXPARW(NS,NT), CONPARW(NS,NT), TOXPARW(1,NS,NT), TOXPARB(1,NS,NT)
            if( ISO > 0 ) GOTO 100
          enddo
        endif
      enddo
    endif

    call Broadcast_Array(ITXPARW, master_id)
    call Broadcast_Array(CONPARW, master_id)
    call Broadcast_Array(TOXPARW, master_id)
    call Broadcast_Array(TOXPARB, master_id)

    !C45A*  READ TOXIC CONTAMINANT ORGANIC CARBON parameterS
    NCARD = '45A'
    if( process_id == master_id )then
      call SEEK('C45A',0)
      if( NTOX > 0 )then
        read(1,*,IOSTAT = ISO)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
        if( ISTRAN(5) == 0 )then
          ISTDOCW = 0
          ISTPOCW = 0
          ISTDOCB = 0
          ISTPOCB = 0
        endif
        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
        if( ISO > 0 ) GOTO 100
      endif
    endif

    call Broadcast_Scalar(ISTDOCW,  master_id)
    call Broadcast_Scalar(ISTPOCW,  master_id)
    call Broadcast_Scalar(ISTDOCB,  master_id)
    call Broadcast_Scalar(ISTPOCB,  master_id)
    call Broadcast_Scalar(STDOCWC,  master_id)
    call Broadcast_Scalar(STPOCWC,  master_id)
    call Broadcast_Scalar(STDOCBC,  master_id)
    call Broadcast_Scalar(STPOCBC,  master_id)

    !C45B* READ TOXIC CONTAMINANT-ORGANIC CARBON INTERACTION parameterS
    NCARD = '45B'
    if( process_id == master_id )then
      call SEEK('C45B',0)
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO) NDUM1, NDUM2, ITXPARWC(1,NT), CONPARWC(1,NT), TOXPARWC(1,NT), TOXPARBC(1,NT)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) NDUM1, NDUM2, ITXPARWC(1,NT), CONPARWC(1,NT), TOXPARWC(1,NT), TOXPARBC(1,NT)
        if( ISO > 0 ) GOTO 100

        read(1,*,IOSTAT = ISO) NDUM1, NDUM2, ITXPARWC(2,NT), CONPARWC(2,NT), TOXPARWC(2,NT), TOXPARBC(2,NT)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) NDUM1, NDUM2, ITXPARWC(2,NT), CONPARWC(2,NT), TOXPARWC(2,NT), TOXPARBC(2,NT)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    TOXPARWC = ABS(TOXPARWC)    ! *** Remove any legacy flags
    TOXPARBC = ABS(TOXPARBC)    ! *** Remove any legacy flags

    call Broadcast_Array(ITXPARWC,  master_id)
    call Broadcast_Array(CONPARWC,  master_id)
    call Broadcast_Array(TOXPARWC,  master_id)
    call Broadcast_Array(TOXPARBC,  master_id)

    !C45C* READ TOXIC CONTAMINANT-ORGANIC CARBON WATER COLUMN POC FRACTIONS
    NCARD = '45C'
    if( process_id == master_id )then
      call SEEK('C45C',0)
      write(mpi_efdc_out_unit,1002) NCARD
      NTMP = NSED+NSND
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO)NDUM1,(FPOCWST(NS,NT),NS = 1,NTMP)

        write(mpi_efdc_out_unit,*)NDUM1,(FPOCWST(NS,NT),NS = 1,NTMP)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(FPOCWST, master_id)

    ! RESET INORGANIC SEDIMENT PARTITION COEFFICIENTS BASED ON
    ! FRACTION OF POC ASSIGNED TO EACH SEDIMENT CLASS IN WATER COLUMN
    do NT = 1,NTOX
      if( ISTOC(NT) >= 2 )then
        if( NSED > 0 )then
          do NS = 1,NSED
            ITXPARW(NS,NT) = 0
            CONPARW(NS,NT) = 0.
            TOXPARW(2:LA,NS,NT) = TOXPARWC(2,NT)
          enddo
        endif
        if( NSND > 0 )then
          do NX = 1,NSND
            NS = NSED+NX
            ITXPARW(NS,NT) = 0
            CONPARW(NS,NT) = 0.
            TOXPARW(2:LA,NS,NT) = TOXPARWC(2,NT)
          enddo
        endif
      endif
    enddo

    !C45D* READ TOXIC CONTAMINANT-ORGANIC CARBON SED BED POC FRACTIONS
    NCARD = '45D'
    if( process_id == master_id )then
      call SEEK('C45D',0)
      write(mpi_efdc_out_unit,1002) NCARD
      NTMP = NSED+NSND
      do NT = 1,NTOX
        read(1,*,IOSTAT = ISO)NDUM1,(FPOCBST(NS,NT),NS = 1,NTMP)

        if( ISO > 0 ) GOTO 100
        write(mpi_efdc_out_unit,*)NDUM1,(FPOCBST(NS,NT),NS = 1,NTMP)
      enddo
    endif
    call Broadcast_Array(FPOCBST, master_id)

    !C45E* READ TOXIC CONTAMINANT ATMOSPHERIC DEPOSITIONS SETTINGS
    NCARD = '45E'
    if( process_id == master_id )then
      call SEEK('C45E',0)
      write(mpi_efdc_out_unit,1002) NCARD

      do NT = 1,NTOX
        read(1,*,ERR = 100) NDUM, TOXDEP(NT).ITXDRY, TOXDEP(NT).TXDRY, TOXDEP(NT).ITXWET, TOXDEP(NT).TXWET

        TOXDEP(NT).TXDRY = TOXDEP(NT).TXDRY/86400.   ! *** CONVERT FROM MG/M2/DAY TO MG/M2/SEC
      enddo
    endif
    do NT = 1,NTOX
      call Broadcast_Scalar(TOXDEP(NT).ITXDRY, master_id)
      call Broadcast_Scalar(TOXDEP(NT).TXDRY,  master_id)
      call Broadcast_Scalar(TOXDEP(NT).ITXWET, master_id)
      call Broadcast_Scalar(TOXDEP(NT).TXWET,  master_id)
    enddo

    ! RESET INORGANIC SEDIMENT PARTITION COEFFICIENTS BASED ON
    ! FRACTION OF POC ASSIGNED TO EACH SEDIMENT CLASS IN SEDIMENT BED
    do NT = 1,NTOX
      if( ISTOC(NT) >= 2 )then
        if( NSED > 0 )then
          do NS = 1,NSED
            TOXPARB(2:LA,NS,NT) = TOXPARBC(2,NT)
          enddo
        endif
        if( NSND > 0 )then
          do NX = 1,NSND
            NS = NSED+NX
            TOXPARB(2:LA,NS,NT) = TOXPARBC(2,NT)
          enddo
        endif
      endif
    enddo
  endif  ! *** END OF NTOX>0 BLOCK

  !C46*  READ BUOYANCY, TEMPERATURE, DYE DATA AND CONCENTRATION BC DATA
  NCARD = '46'
  IBSC = 0
  if( process_id == master_id )then
    call SEEK('C46',0)
    read(1,*,IOSTAT = ISO) BSC, IBSC, TEMO, tmp, tmp, NCBS, NCBW, NCBE, NCBN, IVOLTEMP, VOL_VEL_MAX, VOL_DEP_MIN

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) BSC, IBSC, TEMO, tmp, tmp, NCBS, NCBW, NCBE, NCBN, IVOLTEMP, VOL_VEL_MAX, VOL_DEP_MIN
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(BSC         , master_id)
  call Broadcast_Scalar(IBSC        , master_id)
  call Broadcast_Scalar(TEMO        , master_id)
  call Broadcast_Scalar(NCBS        , master_id)
  call Broadcast_Scalar(NCBW        , master_id)
  call Broadcast_Scalar(NCBE        , master_id)
  call Broadcast_Scalar(NCBN        , master_id)
  call Broadcast_Scalar(IVOLTEMP    , master_id)
  call Broadcast_Scalar(VOL_VEL_MAX , master_id)
  call Broadcast_Scalar(VOL_DEP_MIN , master_id)

  ! *** ENSURE CONSISTENCY OF TEMPERATURE OPTIONS
  if( ISTRAN(2) > 0 )then
    IVOLTEMP = 0    ! *** FORCE THE use OF COMPUTED TEMPERATURES
  elseif( IVOLTEMP > 1 )then
    if( ISTRAN(3) == 0 .and. ISTRAN(5) == 0 ) IVOLTEMP = 0       ! *** NEITHER VOLALILIZATION OR TOXICS ARE ACTIVE, SO DISABLE READING/PROCESSING IVOLTEMP
    if( NCSER(2) < IVOLTEMP-1 ) CALL STOPP('Temperature option IVOLTEMP points to an invalid temperature time series!')
  endif

  ! *** DISABLE BOUYANCY IF ALL CONSTITUENTS IMPACTING DENSITY ARE OFF (7.2)
  if( BSC > 0. .and. (ISTRAN(1) < 1 .and. ISTRAN(2) < 1 .and.  .not. (ISTRAN(6) > 0 .or. ISTRAN(7) > 0)) )then
    BSC = 0.
    IBSC = 0
  endif

  if( TEMO < 0.0 )then
    TEMO = ABS(TEMO)
    INITTEMP = 1
  else
    INITTEMP = 0
  endif

  if( ISICE > 0 )then
    !C46A*   READ ICE EFFECTS
    NCARD = '46A'
    if( process_id == master_id )then
      call SEEK('C46A',0)
      read(1,*,IOSTAT = ISO) ISICE, NISER, TEMPICE, CDICE, ICETHMX, RICETHK0

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) ISICE, NISER, TEMPICE, CDICE, ICETHMX, RICETHK0
      if( ISO > 0 ) GOTO 100
    endif
  endif
  call Broadcast_Scalar(ISICE   , master_id)
  call Broadcast_Scalar(NISER   , master_id)
  call Broadcast_Scalar(TEMPICE , master_id)
  call Broadcast_Scalar(CDICE   , master_id)
  call Broadcast_Scalar(ICETHMX , master_id)
  call Broadcast_Scalar(RICETHK0, master_id)

  RHOI = 917.
  RISEVEL = 0.01
  if( ISICE > 2 )then
    !C46B*   READ ICE MODULE parameterS
    NCARD = '46B'
    if( process_id == master_id )then
      call SEEK('C46B',0)
      read(1,*,IOSTAT = ISO) HWI, ICEK, ALBEDOI, BETAI, GAMMAI, MINICETHICK, RHOI, ISRHEVAP, RISEVEL, MELTFACTOR, AFWI, BFWI, CFWI
      write(mpi_efdc_out_unit,1002) NCARD

      write(mpi_efdc_out_unit,*) HWI, ICEK, ALBEDOI, BETAI, GAMMAI, MINICETHICK, ICETHMX, RHOI, ISRHEVAP, RISEVEL, MELTFACTOR, AFWI, BFWI, CFWI
      if( ISO > 0 ) GOTO 100
    endif
    call Broadcast_Scalar(HWI           , master_id)
    call Broadcast_Scalar(ICEK          , master_id)
    call Broadcast_Scalar(ALBEDOI       , master_id)
    call Broadcast_Scalar(BETAI         , master_id)
    call Broadcast_Scalar(GAMMAI        , master_id)
    call Broadcast_Scalar(MINICETHICK   , master_id)
    call Broadcast_Scalar(ICETHMX       , master_id)
    call Broadcast_Scalar(RHOI          , master_id)
    call Broadcast_Scalar(ISRHEVAP      , master_id)
    call Broadcast_Scalar(RISEVEL       , master_id)
    call Broadcast_Scalar(MELTFACTOR    , master_id)
    call Broadcast_Scalar(AFWI          , master_id)
    call Broadcast_Scalar(BFWI          , master_id)
    call Broadcast_Scalar(CFWI          , master_id)

    if( CFWI <= 0. )then
      AFWI = 9.2    ! *** WIND FUNCTION A, Edinger, et. al. (1974)
      BFWI = 0.46   ! *** WIND FUNCTION B, Edinger, et. al. (1974)
      CFWI = 2.0    ! *** WIND FUNCTION C, Edinger, et. al. (1974)
    else
      ! *** Units Conversion = 1 mmhg = 1.33322 millibars
      AFWI = AFWI*1.33322
      BFWI = BFWI*1.33322
    endif

    ! *** FRAZIL ICE TRANSPORT DEFAULTS
    ISFCT(10)  = 0

  endif
  if( ISICE == 2) NISER = 1

  if( NASER > 0 )then
    ! *** READING TWO NEW CARDS FOR EFDC_073 INSTAED OF ASER.INP FOR 072
    DS_LAT = 0.0
    DS_LONG = 0.0
    COMPUTESOLRAD = .FALSE.

    NCARD = '46C'
    if( process_id == master_id )then
      call SEEK('C46C',0)
      read(1,*,IOSTAT = ISO) DS_LONG, DS_LAT, COMPUTESOLRAD, USESHADE, IEVAP, WINDFA, WINDFB, WINDFC, PBLZ, TEM_HRZ, COARE_NITS

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) DS_LONG, DS_LAT, COMPUTESOLRAD, USESHADE, IEVAP, WINDFA, WINDFB, WINDFC, PBLZ, TEM_HRZ, COARE_NITS
      if( ISO > 0 ) GOTO 100
    endif

    call Broadcast_Scalar(DS_LONG       , master_id)
    call Broadcast_Scalar(DS_LAT        , master_id)
    call Broadcast_Scalar(COMPUTESOLRAD , master_id)
    call Broadcast_Scalar(USESHADE      , master_id)
    call Broadcast_Scalar(IEVAP         , master_id)
    call Broadcast_Scalar(WINDFA        , master_id)
    call Broadcast_Scalar(WINDFB        , master_id)
    call Broadcast_Scalar(WINDFC        , master_id)
    call Broadcast_Scalar(PBLZ          , master_id)
    call Broadcast_Scalar(TEM_HRZ       , master_id)
    call Broadcast_Scalar(COARE_NITS    , master_id)

    AFW = 9.2        ! *** Edinger, et. al. (1974)  W/m2/mmHg
    BFW = 0.46       ! *** Edinger, et. al. (1974)  W/m2/mmHg
    CFW = 2.0        ! *** Edinger, et. al. (1974)  W/m2/mmHg

    ! *** CONVERT WIND FACTOR COEFFICIENTS FROM W/M2/MILLIBAR TO M/S/MILLIBAR
    ! *** Latent Heat of Evaporation = 2259 KJ/KG
    ! *** Density of Water           = 1000 kg/m3
    WINDFA = WINDFA*4.42674E-10
    WINDFB = WINDFB*4.42674E-10
    WINDFC = WINDFC*4.42674E-10

    NCARD = '46D'
    if( process_id == master_id )then
      call SEEK('C46D',0)
      ! *** TBEDIT deprecated.  Changed to TEMBO
      ! *** DABEDT deprecated.  Changed to TEMTHKO
      read(1,*,IOSTAT = ISO) IASWRAD, REVC, RCHC, ISVHEAT, SWRATNF, SWRATNS, FSWRATF, TEMTHKO, TEMBO, HTBED1, HTBED2, WQKETSS, FSOLRADMIN
      write(mpi_efdc_out_unit,1002) NCARD

      write(mpi_efdc_out_unit,*) IASWRAD, REVC, RCHC, ISVHEAT, SWRATNF, SWRATNS, FSWRATF, TEMTHKO, TEMBO, HTBED1, HTBED2, WQKETSS, FSOLRADMIN
      if( ISO > 0 ) GOTO 100
    endif
    WQKEB(1) = SWRATNF
    If( ISTRAN(6) == 0 .and. ISTRAN(7) == 0  )then
      WQKETSS = 0.0
    endif

    call Broadcast_Scalar(IASWRAD   , master_id)
    call Broadcast_Scalar(REVC      , master_id)
    call Broadcast_Scalar(RCHC      , master_id)
    call Broadcast_Scalar(ISVHEAT   , master_id)
    call Broadcast_Scalar(SWRATNF   , master_id)
    call Broadcast_Scalar(SWRATNS   , master_id)
    call Broadcast_Scalar(FSWRATF   , master_id)
    call Broadcast_Scalar(TEMTHKO   , master_id)
    call Broadcast_Scalar(TEMBO     , master_id)
    call Broadcast_Scalar(HTBED1    , master_id)
    call Broadcast_Scalar(HTBED2    , master_id)
    call Broadcast_Scalar(WQKEB(1)  , master_id)
    call Broadcast_Scalar(WQKETSS   , master_id)

  endif

  ! *** READ DYE TYPES
  if( ISTRAN(3) > 0 .and. NDYE > 0 )then
    NCARD = '46E'
    if( process_id == master_id )then
      call SEEK('C46E',0)
      do MD = 1,NDYE
        read(1,*,IOSTAT = ISO) NDUM1, DYES(MD).ITYPE, DYES(MD).KRATE0, DYES(MD).KRATE1, DYES(MD).TADJ, DYES(MD).TREF, &
          DYES(MD).SETTLE, DYES(MD).ICFLAG, DYES(MD).IC, DYES(MD).WCLIMIT
        if( ISO > 0 ) GOTO 100

        write(mpi_efdc_out_unit,*) NDUM1, DYES(MD).ITYPE, DYES(MD).KRATE0, DYES(MD).KRATE1, DYES(MD).TADJ, DYES(MD).TREF, &
          DYES(MD).SETTLE, DYES(MD).ICFLAG, DYES(MD).IC, DYES(MD).WCLIMIT
      enddo
    endif

    !C46F*  READ DYE VOLATILIZATION parameterS
    NCARD = '46F'
    if( process_id == master_id )then
      call SEEK('C46F',0)
      do MD = 1,NDYE
        read(1,*,IOSTAT = ISO) NDUM, DYES(MD).VOL.KL_OPT, DYES(MD).VOL.MW, DYES(MD).VOL.HE, DYES(MD).VOL.TCOEFF, DYES(MD).VOL.AIRCON, DYES(MD).VOL.MULT
        write(mpi_efdc_out_unit,1002) NCARD

        write(mpi_efdc_out_unit,*) NDUM, DYES(MD).VOL.KL_OPT, DYES(MD).VOL.MW, DYES(MD).VOL.HE, DYES(MD).VOL.TCOEFF, DYES(MD).VOL.AIRCON, DYES(MD).VOL.MULT
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    Do MD = 1,NDYE
      call Broadcast_Scalar(DYES(MD).ITYPE      , master_id)
      call Broadcast_Scalar(DYES(MD).KRATE0     , master_id)
      call Broadcast_Scalar(DYES(MD).KRATE1     , master_id)
      call Broadcast_Scalar(DYES(MD).TADJ       , master_id)
      call Broadcast_Scalar(DYES(MD).TREF       , master_id)
      call Broadcast_Scalar(DYES(MD).SETTLE     , master_id)
      call Broadcast_Scalar(DYES(MD).ICFLAG     , master_id)
      call Broadcast_Scalar(DYES(MD).IC         , master_id)
      call Broadcast_Scalar(DYES(MD).WCLIMIT    , master_id)
      call Broadcast_Scalar(DYES(MD).VOL.KL_OPT , master_id)
      call Broadcast_Scalar(DYES(MD).VOL.MW     , master_id)
      call Broadcast_Scalar(DYES(MD).VOL.HE     , master_id)
      call Broadcast_Scalar(DYES(MD).VOL.TCOEFF , master_id)
      call Broadcast_Scalar(DYES(MD).VOL.AIRCON , master_id)
      call Broadcast_Scalar(DYES(MD).VOL.MULT   , master_id)
    enddo


    !C43C*  READ TOXIC TIMING AND VOLATILIZATION SWITCHES
    NCARD = '46G'
    if( process_id == master_id )then
      call SEEK('C46G',0)
      read(1,*,IOSTAT = ISO) DYESTEPW

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) DYESTEPW
      if( ISO > 0 ) GOTO 100
    endif
  endif

  ! ****************************************************************************C
  if( NCBS > 0 )then

    !C47*  READ LOCATIONS OF CONC BC'S ON SOUTH BOUNDARIES
    NCARD = '47'
    if( process_id == master_id )then
      call SEEK('C47',0)
      do L = 1,NCBS
        read(1,*,IOSTAT = ISO) ICBS_GL(L),JCBS_GL(L),NTSCRS_GL(L),NCSERS_GL(L,1),NCSERS_GL(L,2),NCSERS_GL(L,3),NCSERS_GL(L,4),NCSERS_GL(L,5),NCSERS_GL(L,6),NCSERS_GL(L,7)
        write(mpi_efdc_out_unit,1002) NCARD

        write(mpi_efdc_out_unit,*) ICBS_GL(L),JCBS_GL(L),NTSCRS_GL(L),NCSERS_GL(L,1),NCSERS_GL(L,2),NCSERS_GL(L,3),NCSERS_GL(L,4),NCSERS_GL(L,5),NCSERS_GL(L,6),NCSERS_GL(L,7)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(ICBS_GL  , master_id)
    call Broadcast_Array(JCBS_GL  , master_id)
    call Broadcast_Array(NTSCRS_GL, master_id)
    call Broadcast_Array(NCSERS_GL, master_id)

    !C48*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '48'
    if( process_id == master_id )then
      call SEEK('C48',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBS
        read(1,*,IOSTAT = ISO) (CBS_GL(L,1,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBS_GL(L,1,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBS_GL, master_id)

    !C49*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '49'
    if( process_id == master_id )then
      call SEEK('C49',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBS
        read(1,*,IOSTAT = ISO) (CBS_GL(L,1,M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBS_GL(L,1,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBS_GL, master_id)

    !C50*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '50'
    if( process_id == master_id )then
      call SEEK('C50',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBS
        read(1,*,IOSTAT = ISO) (CBS_GL(L,2,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBS_GL(L,2,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBS_GL, master_id)

    !C51*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '51'
    if( process_id == master_id )then
      call SEEK('C51',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBS
        read(1,*,IOSTAT = ISO) (CBS_GL(L,2,M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBS_GL(L,2,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBS_GL, master_id)

  endif


  if( NCBW > 0 )then
    !C52*  READ LOCATIONS OF CONC BC'S ON WEST BOUNDARIES
    NCARD = '52'
    if( process_id == master_id )then
      call SEEK('C52',0)
      do L = 1,NCBW
        read(1,*,IOSTAT = ISO) ICBW_GL(L),JCBW_GL(L),NTSCRW_GL(L),NCSERW_GL(L,1),NCSERW_GL(L,2),NCSERW_GL(L,3),NCSERW_GL(L,4),NCSERW_GL(L,5),NCSERW_GL(L,6),NCSERW_GL(L,7)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) ICBW_GL(L),JCBW_GL(L),NTSCRW_GL(L),NCSERW_GL(L,1),NCSERW_GL(L,2),NCSERW_GL(L,3),NCSERW_GL(L,4),NCSERW_GL(L,5),NCSERW_GL(L,6),NCSERW_GL(L,7)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(ICBW_GL, master_id)
    call Broadcast_Array(JCBW_GL, master_id)
    call Broadcast_Array(NTSCRW_GL, master_id)
    call Broadcast_Array(NCSERW_GL, master_id)

    !C53*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '53'
    if( process_id == master_id )then
      call SEEK('C53',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBW
        read(1,*,IOSTAT = ISO) (CBW_GL(L,1,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBW_GL(L,1,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBW_GL, master_id)

    !C54*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '54'
    if( process_id == master_id )then
      call SEEK('C54',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBW
        read(1,*,IOSTAT = ISO) (CBW_GL(L,1,M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBW_GL(L,1,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBW_GL, master_id)

    !C55*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '55'
    if( process_id == master_id )then
      call SEEK('C55',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBW
        read(1,*,IOSTAT = ISO) (CBW_GL(L,2,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBW_GL(L,2,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBW_GL, master_id)

    !C56*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '56'
    if( process_id == master_id )then
      call SEEK('C56',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBW
        read(1,*,IOSTAT = ISO) (CBW_GL(L,2,M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBW_GL(L,2,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBW_GL, master_id)

  endif

  if( NCBE > 0 )then
    !C57*  READ LOCATIONS OF CONC BC'S ON EAST BOUNDARIES
    NCARD = '57'
    if( process_id == master_id )then
      call SEEK('C57',0)
      do L = 1,NCBE
        read(1,*,IOSTAT = ISO) ICBE_GL(L),JCBE_GL(L),NTSCRE_GL(L),NCSERE_GL(L,1),NCSERE_GL(L,2),NCSERE_GL(L,3),NCSERE_GL(L,4),NCSERE_GL(L,5),NCSERE_GL(L,6),NCSERE_GL(L,7)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) ICBE_GL(L),JCBE_GL(L),NTSCRE_GL(L),NCSERE_GL(L,1),NCSERE_GL(L,2),NCSERE_GL(L,3),NCSERE_GL(L,4),NCSERE_GL(L,5),NCSERE_GL(L,6),NCSERE_GL(L,7)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(ICBE_GL,   master_id)
    call Broadcast_Array(JCBE_GL,   master_id)
    call Broadcast_Array(NCSERE_GL, master_id)
    call Broadcast_Array(NTSCRE_GL,  master_id)

    !C58*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '58'
    if( process_id == master_id )then
      call SEEK('C58',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBE
        read(1,*,IOSTAT = ISO) (CBE_GL(L,1,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBE_GL(L,1,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBE_GL, master_id)

    !C59*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '59'
    if( process_id == master_id )then
      call SEEK('C59',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBE
        read(1,*,IOSTAT = ISO) (CBE_GL(L,1,M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBE_GL(L,1,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBE_GL, master_id)

    !C60*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '60'
    if( process_id == master_id )then
      call SEEK('C60',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBE
        read(1,*,IOSTAT = ISO) (CBE_GL(L,2,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBE_GL(L,2,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBE_GL, master_id)

    !C61*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '61'
    if( process_id == master_id )then
      call SEEK('C61',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBE
        read(1,*,IOSTAT = ISO) (CBE_GL(L,2,M),M = MMIN,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBE_GL(L,2,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBE_GL, master_id)

  endif    ! *** NCBE>0

  if( NCBN > 0 )then
    !C62*  READ LOCATIONS OF CONC BC'S ON NORTH BOUNDARIES
    NCARD = '62'
    if( process_id == master_id )then
      call SEEK('C62',0)
      do L = 1,NCBN
        read(1,*,IOSTAT = ISO) ICBN_GL(L),JCBN_GL(L),NTSCRN_GL(L),NCSERN_GL(L,1),NCSERN_GL(L,2),NCSERN_GL(L,3),NCSERN_GL(L,4),NCSERN_GL(L,5),NCSERN_GL(L,6),NCSERN_GL(L,7)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) ICBN_GL(L),JCBN_GL(L),NTSCRN_GL(L),NCSERN_GL(L,1),NCSERN_GL(L,2),NCSERN_GL(L,3),NCSERN_GL(L,4),NCSERN_GL(L,5),NCSERN_GL(L,6),NCSERN_GL(L,7)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(ICBN_GL, master_id)
    call Broadcast_Array(JCBN_GL, master_id)
    call Broadcast_Array(NTSCRN_GL, master_id)
    call Broadcast_Array(NCSERN_GL, master_id)

    !C63*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '63'
    if( process_id == master_id )then
      call SEEK('C63',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBN
        read(1,*,IOSTAT = ISO) (CBN_GL(L,1,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBN_GL(L,1,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo

    endif
    call Broadcast_Array(CBN_GL, master_id)

    !C64*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '64'
    if( process_id == master_id )then
      call SEEK('C64',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBN
        read(1,*,IOSTAT = ISO) (CBN_GL(L,1,M),M = MMIN,MMAX)
        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBN_GL(L,1,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo

    endif
    call Broadcast_Array(CBN_GL, master_id)

    !C65*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD = '65'
    if( process_id == master_id )then
      call SEEK('C65',0)
      MMAX = 3 + NDYM + NTOX
      do L = 1,NCBN
        read(1,*,IOSTAT = ISO) (CBN_GL(L,2,M),M = 1,MMAX)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBN_GL(L,2,M),M = 1,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array(CBN_GL, master_id)

    !C66*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD = '66'
    if( process_id == master_id )then
      call SEEK('C66',0)
      MMIN = MMAX+1
      MMAX = MMAX+NSED+NSND
      do L = 1,NCBN
        read(1,*,IOSTAT = ISO) (CBN_GL(L,2,M),M = MMIN,MMAX)
        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*) (CBN_GL(L,2,M),M = MMIN,MMAX)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    call Broadcast_Array(CBN_GL, master_id)

  endif    ! *** NCBN>0

  !C66A*  READ CONCENTRATION DATA ASSIMILATION parameterS
  NCARD = '66A'

  !C66B*  READ CONCENTRATION DATA ASSIMILATION LOCATIONS AND SERIES IDENTIFIERS
  NCARD = '66B'

  !C67*  READ NEUTRALLY BUOYANT PARTICLE DRIFTER DATA
  NCARD = '67'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C67',0)
    read(1,*,IOSTAT = ISO) ISPD

    NPD = 0
    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ISPD
    ! ***
    if( ISO > 0 ) GOTO 100
  endif

  call Broadcast_Scalar(ISPD      , master_id)

  !C68*  READ NEUTRALLY BUOYANT PARTICLE INITIAL POSITIONS     DELME - These are not used in EFDC+  See mod_drifter.f90
  NCARD = '68'

  !C69*  CONSTANTS FOR LONGITUDE AND LATITUDE OF CELL CENTERS (NOT USED)
  NCARD = '69'

  !C70*  CONTROLS FOR WRITING ASCII OR BINARY DUMP FILES  (NOT USED)
  NCARD = '70'

  !C71*  CONTROLS FOR HORIZONTAL PLANE SCALAR FIELD CONTOURING
  NCARD = '71'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C71',0)
    do NS = 1,7
      read(1,*,IOSTAT = ISO) ISSPH(NS),ldum,ISRSPH(NS),ISPHXY(NS)

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) ISSPH(NS),ldum,ISRSPH(NS),ISPHXY(NS)
    enddo
    if( ISO > 0 ) GOTO 100

  endif
  call Broadcast_Array(ISSPH , master_id)
  call Broadcast_Array(ISRSPH, master_id)
  call Broadcast_Array(ISPHXY, master_id)

  ! *** SET WATER COLUMN LINKAGE FLAG IF ANY CONSTITUENT OUTPUT IS ENABLED
  do NS = 1,7
    if( ISTRAN(NS) >= 1 ) ISSPH(8) = 1
  enddo
  if( ISSPH(4) >= 1 ) ISSPH(8) = 1

  !C71A*  CONTROLS FOR HORIZONTAL PLANE SEDIMENT BED PROPERTIES
  NCARD = '71A'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C71A',0)
    read(1,*,IOSTAT = ISO) ITMP,ISBEXP,NPBPH

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ISBEXP,NPBPH
    if( ISO > 0 ) GOTO 100
  endif
  call Broadcast_Scalar(ISBEXP , master_id)
  call Broadcast_Scalar(NPBPH  , master_id)

  !C71B*  CONTROLS FOR FOOD CHAIN MODEL OUTPUT
  NCARD = '71B'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C71B',0)
    read(1,*,IOSTAT = ISO) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
    if( ISO > 0 ) GOTO 100
  endif
  call Broadcast_Scalar(ISFDCH , master_id)
  call Broadcast_Scalar(NFDCHZ , master_id)
  call Broadcast_Scalar(HBFDCH , master_id)
  call Broadcast_Scalar(TFCAVG , master_id)

  !C72*  CONTROLS FOR EFDC_EXPLORER LINKAGE AND SURFACE ELEVATION RESIDUAL OUTPUT
  JSFDCH = 1
  NCARD = '72'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C72',0)
    read(1,*,IOSTAT = ISO) ISPPH,NPPPH,ldum,ldum

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) ISPPH,NPPPH,ldum,ldum
    if( ISO > 0 ) GOTO 100
    if( ISPPH < 0 ) ISPPH = 1
  endif
  call Broadcast_Scalar(ISPPH  , master_id)
  call Broadcast_Scalar(NPPPH  , master_id)

  !C73*  CONTROLS FOR HORIZONTAL PLANE VELOCITY PLOTTING
  NCARD = '73'

  !C74*  CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING
  NCARD = '74'

  !C75*  MORE CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING
  NCARD = '75'

  !C76*  I,J LOCATIONS DEFINING VERTICAL PLANE FOR CONTOURING
  NCARD = '76'

  !C77*  CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
  NCARD = '77'

  !C78*   MORE CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
  NCARD = '78'

  !C79*   MORE CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
  NCARD = '79'

  !C80*  CONTROLS FOR 3D FIELD OUTPUT
  NCARD = '80'

  !C81* OUTPUT ACTIVATION AND SCALES FOR 3D FIELD OUTPUT
  NCARD = '81'

  !C82* INPLACE HARMONIC ANALYSIS parameterS
  NCARD = '82'

  !C83* HARMONIC ANALYSIS LOCATIONS AND SWITCHES
  NCARD = '83'

  !C84* CONTROLS FOR WRITING TO TIME SERIES FILES
  NCARD = '84'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C84',0)
    read(1,*,IOSTAT = ISO)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR
    if( ISO > 0 ) GOTO 100
  endif
  call Broadcast_Scalar(ISTMSR   , master_id)
  call Broadcast_Scalar(MLTMSR   , master_id)
  call Broadcast_Scalar(NWTMSR   , master_id)
  call Broadcast_Scalar(NBTMSR   , master_id)
  call Broadcast_Scalar(NSTMSR   , master_id)
  call Broadcast_Scalar(NTSSTSP  , master_id)
  call Broadcast_Scalar(TCTMSR   , master_id)

  JSTMSR = 1
  NCTMSR = 1

  if( NTSSTSP > 0 )then
    NCARD = '85'
    if( process_id == master_id )then
      call SEEK('C85',0)
      do ITSSS = 1,NTSSTSP
        read(1,*,IOSTAT = ISO)IDUM,MTSSTSP(ITSSS)
        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*)IDUM,MTSSTSP(ITSSS)
        if( ISO > 0 ) GOTO 100
      enddo
    endif

    call Broadcast_Array(MTSSTSP, master_id)

    NCARD = '86'
    if( process_id == master_id )then
      call SEEK('C86',0)
      do ITSSS = 1,NTSSTSP
        do MTSSS = 1,MTSSTSP(ITSSS)
          read(1,*,IOSTAT = ISO)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),TSSTOP(MTSSS,ITSSS)

          write(mpi_efdc_out_unit,1002) NCARD
          write(mpi_efdc_out_unit,*)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),TSSTOP(MTSSS,ITSSS)
          if( ISO > 0 ) GOTO 100
        enddo
      enddo
    endif
    call Broadcast_Array( TSSTRT, master_id)
    call Broadcast_Array( TSSTOP, master_id)

  endif

  if( MLTMSR > 0 )then
    NCARD = '87'
    if( process_id == master_id )then
      call SEEK('C87',0)
      do M = 1,MLTMSR
        read(1,*,IOSTAT = ISO)ILTMSR_GL(M),JLTMSR_GL(M),NTSSSS_GL(M),MTMSRP_GL(M),MTMSRC_GL(M),MTMSRA_GL(M),&
          MTMSRUE_GL(M),MTMSRUT_GL(M),MTMSRU_GL(M),MTMSRQE_GL(M),MTMSRQ_GL(M),CLTMSR_GL(M)

        write(mpi_efdc_out_unit,1002) NCARD
        write(mpi_efdc_out_unit,*)ILTMSR_GL(M),JLTMSR_GL(M),NTSSSS_GL(M),MTMSRP_GL(M),MTMSRC_GL(M),MTMSRA_GL(M),&
          MTMSRUE_GL(M),MTMSRUT_GL(M),MTMSRU_GL(M),MTMSRQE_GL(M),MTMSRQ_GL(M),CLTMSR_GL(M)
        if( ISO > 0 ) GOTO 100
      enddo
    endif
    call Broadcast_Array( ILTMSR_GL,  master_id)
    call Broadcast_Array( JLTMSR_GL,  master_id)
    call Broadcast_Array( NTSSSS_GL,  master_id)
    call Broadcast_Array( MTMSRP_GL,  master_id)
    call Broadcast_Array( MTMSRC_GL,  master_id)
    call Broadcast_Array( MTMSRA_GL,  master_id)
    call Broadcast_Array( MTMSRUE_GL, master_id)
    call Broadcast_Array( MTMSRUT_GL, master_id)
    call Broadcast_Array( MTMSRU_GL , master_id)
    call Broadcast_Array( MTMSRQE_GL, master_id)
    call Broadcast_Array( MTMSRQ_GL , master_id)
    !Call Broadcast_Array(CLTMSR , master_id)  ! *** Character    DELME - TODO
  endif

  !C88 ** HIGH FREQUENCY OUTPUT
  NCARD = '88'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C88',0)
    read(1,*,IOSTAT = ISO) HFREOUT

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) HFREOUT
    if( ISO > 0 ) GOTO 100
  endif
  call Broadcast_Scalar(HFREOUT, master_id)

  ! *** NETCDF GENERATION CONTROLS
  NCARD = '91'
  ! *** ********************************************************
  if( process_id == master_id )then
    call SEEK('C91',0)
    read(1,'(A)',IOSTAT = ISO) STR
    !                                       NOT      NOT
    read(STR,*,ERR = 100) NCDFOUT,DEFLEV,ROTA,BLK,UTMZ,HREST,BASEDATE,BASETIME,PROJ
    if( UTMZ < 0 )then
      UTMZ = ABS(UTMZ)
      HEMI = 2
    else
      HEMI = 1
    endif

    write(mpi_efdc_out_unit,1002) NCARD
    write(mpi_efdc_out_unit,*) NCDFOUT,DEFLEV,ROTA,BLK,HEMI,UTMZ,HREST,BASEDATE,BASETIME,PROJ
  endif

  call Broadcast_Scalar( NCDFOUT   , master_id)
  call Broadcast_Scalar( DEFLEV    , master_id)
  call Broadcast_Scalar( ROTA      , master_id)
  call Broadcast_Scalar( BLK       , master_id)
  call Broadcast_Scalar( HEMI      , master_id)
  call Broadcast_Scalar( UTMZ      , master_id)
  call Broadcast_Scalar( HREST     , master_id)
  call Broadcast_Scalar( BASEDATE  , master_id)
  call Broadcast_Scalar( BASETIME  , master_id)
  call Broadcast_Scalar( PROJ      , master_id)

  IS_NC_OUT = 0
  if( NCDFOUT > 0 )then
    NCARD = '91A'
    if( process_id == master_id )then
      call SEEK('C91A',0)
      read(1,*,ERR = 100) ISSGLFIL,TBEGNCDF,TENDNCDF,NCFREQ

      NCARD = '91B'
      call SEEK('C91B',0)
      read(1,*,ERR = 100) (IS_NC_OUT(i), i = 1,25)

      write(mpi_efdc_out_unit,1002) NCARD
      write(mpi_efdc_out_unit,*) IS_NC_OUT(1:25)

    endif

    call Broadcast_Scalar( ISSGLFIL, master_id)
    call Broadcast_Scalar( TBEGNCDF, master_id)
    call Broadcast_Scalar( TENDNCDF, master_id)
    call Broadcast_Scalar( NCFREQ,   master_id)

    call Broadcast_Array( IS_NC_OUT, master_id)

  endif

  NTS = INT8(NTC)*NTSPTC
  NSVSFP = 0

  DT = TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)

  GOTO 2000


  ! *** ********************************************************
  ! *** WRITE INPUT ERROR MESSAGES AND TERMINATE RUN
100 WRITE(6,1001)NCARD
  write(mpi_error_unit,1001)NCARD
  write(mpi_error_unit,1001)NCARD

  call STOPP('')

2000 continue

  ! *** READ CELL TYPES FROM FILES CELL.INP
  if( process_id == master_id )then
    write(*,'(A)')'READING CELL.INP'
    open(1,FILE = 'cell.inp',STATUS = 'UNKNOWN')

    ! *** SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT
    STRC = READSTR(1)
    read(STRC,*)JCTMP

    if( JCTMP /= JC_GLOBAL )then

      ! ***   READ OLD FILE FORMAT
      JACROSS = JC_Global
      if(JC_Global > 640) JACROSS = 640
      do JT = 1,JC_GLOBAL,JACROSS
        JF = JT
        JLAST = JT + JACROSS-1
        if( JLAST > JC_GLOBAL) JLAST = JC_GLOBAL
        do I = 1,IC_GLOBAL
          read(1,6,IOSTAT = ISO) (IJCT_GLOBAL(I,J),J = JF,JLAST)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
        enddo
      enddo

    else

      if( IC_Global > 640 )then
        IACROSS = 640
        do IT = 1,IC_Global,IACROSS
          IFIRST = IT
          ILAST = IT + IACROSS - 1
          if( ILAST>IC_Global) ILAST = IC_GLOBAL

          do J = JC_GLOBAL,1,-1
            read(1,66,IOSTAT = ISO) ADUMMY, (IJCT_GLOBAL(I,J),I = IFIRST,ILAST)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
          enddo
        enddo
      else
        IFIRST = 1
        ILAST = IC_GLOBAL
        do J = JC_GLOBAL,1,-1
          read(1,66,IOSTAT = ISO) ADUMMY, (IJCT_GLOBAL(I,J),I = IFIRST,ILAST)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
        enddo
      endif
    endif
    close(1)
  endif
  call Broadcast_Array(IJCT_GLOBAL, master_id)

8 FORMAT ('   CELL TYPE ARRAY,J = ',I5,2X,'TO J = ',I5,//)

  !----------------------------------------------------------------------C
  ! *** ********************************************************
  if( process_id == master_id )then

    write(*,'(A)')'READING CELLLT.INP'
    open(1,FILE = 'celllt.inp',STATUS = 'UNKNOWN')

    ! *** SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT
    STRC = READSTR(1)
    read(STRC,*)JCTMP

    if( JCTMP /= JC_GLOBAL )then
      ! ***   READ OLD FILE FORMAT
      JACROSS = JC_Global
      if( JC_Global > 640 ) JACROSS = 640
      do JT = 1,JC_Global,JACROSS
        JF = JT
        JLAST = JT + JACROSS - 1
        if( JLAST>JC_Global) JLAST = JC_Global
        WRITE (7,8) JF, JLAST
        do I = 1,IC_Global
          read(1,6,IOSTAT = ISO) (IJCTLT_Global(I,J),J = JF,JLAST)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELLLT.INP')
          WRITE (7,16) (IJCTLT_Global(I,J),J = JF,JLAST)
        enddo
        write(mpi_efdc_out_unit,15)
      enddo

    else
      ! ***   READ NEW FILE FORMAT
      if( IC_GLOBAL > 640 )then
        IACROSS = 640
        do IT = 1,IC_GLOBAL,IACROSS
          IFIRST = IT
          ILAST = IT+IACROSS-1
          if( ILAST>IC_GLOBAL) ILAST = IC_GLOBAL
          WRITE (mpi_efdc_out_unit,'(2I5)') IFIRST, ILAST
          do J = JC_GLOBAL,1,-1
            read(1,66,IOSTAT = ISO) ADUMMY, (IJCTLT_GLOBAL(I,J),I = IFIRST,ILAST)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELLLT.INP')
            WRITE (mpi_efdc_out_unit,166) ADUMMY, (IJCTLT_GLOBAL(I,J),I = IFIRST,ILAST)
          enddo
          write(mpi_efdc_out_unit,15)
        enddo
      else
        IFIRST = 1
        ILAST = IC_GLOBAL
        WRITE (mpi_efdc_out_unit,88) IFIRST, ILAST
        do J = JC_GLOBAL,1,-1
          read(1,66,IOSTAT = ISO) ADUMMY, (IJCTLT_GLOBAL(I,J),I = IFIRST,ILAST)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELLLT.INP')
          WRITE (mpi_efdc_out_unit,166) ADUMMY, (IJCTLT_GLOBAL(I,J),I = IFIRST,ILAST)
        enddo
        write(mpi_efdc_out_unit,15)
      endif
    endif
    close(1)

88  FORMAT (//,' *** CELLMAP',//,'   CELLLT TYPE ARRAY,  I = ',I5,' TO I = ',I5,//)
  endif

  call Broadcast_Array(IJCTLT_GLOBAL, master_id)

  ! *** ********************************************************
  ! ***  GENERATE CELL MAPPINGS
  call Child_Grid

  ! *** This is a domain specific cell map
  call CELLMAP

  ! *** FORMAT STATEMENTS FOR EFDC.INP
15 FORMAT (/)
6 FORMAT (640I1)
16 FORMAT (1X,640I1)
66 FORMAT (A5,10000I1)
166 FORMAT (1X,A5,10000I1)

  ! *** READ CURVILINEAR-ORTHOGONAL OR VARIABLE CELL DATA FROM FILE DXDY.INP

  ! *** INITIALIZE CELL DIMENSIONS TO CONSTANT CARTESIAN OR DUMMY VALUES
  do L = 1,LC
    DXP(L) = DX*DXYCVT
    DYP(L) = DY*DXYCVT
    ZBR(L) = ZBRADJ
  enddo

  ! *** READ IN DX, DY, DEPTH AND BOTTOM ELEVATION AT CELL CENTERS OF VARIABLE CELLS
  LMHK = .FALSE.
  if( process_id == master_id )then
    write(*,'(A)')'READING DXDY.INP'
    open(1,FILE = 'dxdy.inp',STATUS = 'UNKNOWN')
    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR = READSTR(1)
  endif

  HP_Global(:) = HMIN

  IMN = 1000000
  JMN = 1000000
  IMX = 0
  JMX = 0
  ! *** Global versions
  IMN_Global = 1000000
  JMN_Global = 1000000
  IMX_Global = 0
  JMX_Global = 0

  ! *** READ DXDY WITHOUT VEGETATION
  if( process_id == master_id )then
    MVEGIJT = 0
    do LP = 2,LA_Global
      if( ISVEG > 0 .and. MACDRAG == 0 )then
        read(1,*,IOSTAT = ISO) IGL, JGL, DXIJ, DYIJ, HIJ, BELVIJ, ZBRIJ, MVEGIJT
      else
        read(1,*,IOSTAT = ISO) IGL, JGL, DXIJ, DYIJ, HIJ, BELVIJ, ZBRIJ
      endif
      if( IGL < 2 .or. IGL > IC_GLOBAL-1 ) ISO = 3
      if( JGL < 2 .or. JGL > JC_GLOBAL-1 ) ISO = 4
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DXDY.INP')

      LG = LIJ_GLOBAL(IGL,JGL)

      DXP_Global(LG)  = DXIJ
      DYP_Global(LG)  = DYIJ
      HP_Global(LG)   = HADADJ + HCVRT*HIJ
      HP_Global(LG)   = MAX(HP_Global(LG), HMIN)
      BELV_GLOBAL(LG) = BELADJ + BELCVRT*BELVIJ
      ZBR_GLOBAL(LG)  = ZBRADJ + ZBRCVRT*ZBRIJ
      MVEG_GLOBAL(LG) = MVEGIJT

      ! *** get global versions as they are needed for the NETCDF writing out
      IMN_Global = MIN(IMN_Global, IGL)
      JMN_Global = MIN(JMN_Global, JGL)
      IMX_Global = MAX(IMX_Global, IGL)
      JMX_Global = MAX(JMX_Global, JGL)
    enddo

    close(1)
  endif

  call Broadcast_Array(DXP_Global,  master_id)
  call Broadcast_Array(DYP_Global,  master_id)
  call Broadcast_Array(HP_Global,   master_id)
  call Broadcast_Array(BELV_Global, master_id)
  call Broadcast_Array(ZBR_Global,  master_id)
  call Broadcast_Array(MVEG_Global, master_id)

  call Broadcast_Scalar(IMN_Global, master_id)
  call Broadcast_Scalar(JMN_Global, master_id)
  call Broadcast_Scalar(IMX_Global, master_id)
  call Broadcast_Scalar(JMX_Global, master_id)

  ! *** Map to Local Domain
  do LG = 2,LA_Global
    L = Map2Local(LG).LL
    if( L > 1 )then
      DXP(L) = DXP_Global(LG)
      DYP(L) = DYP_Global(LG)
      HP(L)  = HP_Global(LG)
      HMP(L) = HP_Global(LG)
      BELV(L)  = BELV_Global(LG)
      BELV1(L) = BELV_Global(LG)
      ZBR(L)   = ZBR_Global(LG)
      MVEGL(L) = MVEG_Global(LG)

      ! ***  MHK
      if( MVEGL(L) > 90 )then
        LMHK = .TRUE.
        ITURB = ITURB + 1
        IJLTURB(ITURB,1) = Map2Local(LG).IL
        IJLTURB(ITURB,2) = Map2Local(LG).JL
        IJLTURB(ITURB,3) = L
      endif
    endif
  enddo

  ! *** OPEN FILE MODDXDY.INP TO MODIFY INPUT VALUES OF DX AND DY
  if( IMDXDY > 0 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING MODDXDY.INP'
      open(1,FILE = 'moddxdy.inp',STATUS = 'UNKNOWN')
      ! *** SKIP OVER TITLE AND HEADER LINES
      STR = READSTR(1)
      read(1,*) NMDXDY
    Endif

    call Broadcast_Scalar(NMDXDY, master_id)

    if( NMDXDY >= 1 )then
      do NMD = 1, NMDXDY
        if( process_id == master_id )then
          read(1,*) ITMP, JTMP, RMDX, RMDY
        Endif

        call Broadcast_Scalar(ITMP, master_id)
        call Broadcast_Scalar(JTMP, master_id)
        call Broadcast_Scalar(RMDX, master_id)
        call Broadcast_Scalar(RMDY, master_id)

        ! *** Get local i,j
        I = IG2IL(ITMP)
        J = JG2JL(JTMP)
        ! *** Only select cells in the current domain
        if( I > 0 .and. I <= IC )then
          if( J > 0 .and. J <= JC )then
            LTMP = LIJ(ITMP, JTMP) ! Map to local L value

            DXP(LTMP) = RMDX*DXP(LTMP)
            DYP(LTMP) = RMDY*DYP(LTMP)
          endif
        endif

      enddo
    endif

    ! *** Only have the master close the file
    if( process_id == master_id )then
      close(1)
    endif

  endif

  ! *** OPEN FILE MODCHAN.INP TO INSERT SUBGRID CHANNELS INTO HOST CELLS
  MDCHH = 0
  if( ISCHAN > 0 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING MODCHAN.INP'
      open(1,FILE = 'modchan.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      if( ISCHAN == 1 )then
        read(1,*) MDCHH, ldum, ldum
        read(1,*) ldum, ldum, QCHERR
        if( MDCHH >= 1 )then
          do NMD = 1,MDCHH
            read(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD)
            QCHANU(NMD) = 0.
            QCHANUN(NMD) = 0.
            QCHANV(NMD) = 0.
            QCHANVN(NMD) = 0.
          enddo
        endif
      endif
      if( ISCHAN == 2 )then
        read(1,*) MDCHH, ldum, ldum
        read(1,*) ldum, ldum, QCHERR
        if( MDCHH >= 1 )then
          do NMD = 1,MDCHH
            read(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD),CHANLEN(NMD),PMDCH(NMD)
            QCHANU(NMD) = 0.
            QCHANUN(NMD) = 0.
            QCHANV(NMD) = 0.
            QCHANVN(NMD) = 0.
          enddo
        endif
      endif
      close(1)
      if( MDCHH >= 1 )then
        do NMD = 1,MDCHH
          LMDCHH(NMD) = LIJ(IMDCHH(NMD),JMDCHH(NMD))
          if( IMDCHU(NMD) == 1 .and. JMDCHU(NMD) == 1 )then
            LMDCHU(NMD) = 1
          else
            LMDCHU(NMD) = LIJ(IMDCHU(NMD),JMDCHU(NMD))
          endif
          if( IMDCHV(NMD) == 1 .and. JMDCHV(NMD) == 1 )then
            LMDCHV(NMD) = 1
          else
            LMDCHV(NMD) = LIJ(IMDCHV(NMD),JMDCHV(NMD))
          endif
        enddo
      endif

    endif

    call Broadcast_Array(MDCHTYP , master_id)
    call Broadcast_Array(IMDCHH  , master_id)
    call Broadcast_Array(JMDCHH  , master_id)
    call Broadcast_Array(IMDCHU  , master_id)
    call Broadcast_Array(JMDCHU  , master_id)
    call Broadcast_Array(IMDCHV  , master_id)
    call Broadcast_Array(JMDCHV  , master_id)
    call Broadcast_Array(CHANLEN , master_id)
    call Broadcast_Array(PMDCH   , master_id)

  endif

  ! *** ENABLE SPATIALLY VARIABLE BACKGROUND AHO
  if( AHO < 0. )then
    AHMAX = ABS(AHO)
    AHMIN = 1.0E32
    do L = 2,LA
      AHOXY(L) = ABS(AHO)*DXP(L)*DYP(L)
      AHMAX    = MAX(AHOXY(L),AHMAX)
      AHMIN    = MIN(AHOXY(L),AHMIN)
    enddo

    if( process_id == master_id )then
      PRINT '(A,2E16.5)','VARIABLE AHO USED (MIN,MAX): ', AHMIN, AHMAX
    endif

  else
    AHOXY = AHO
  endif

  ! *** SMAGORINSKY AND BACKGROUND DIFFUSIVITY
  ! *** Constant and/or default values
  AHO = ABS(AHO)

  ! *** ENABLE SPATIALLY VARIABLE SMAGORINSKY AND BACKGROUND DIFFUSIVITY
  if( AHD < 0. )then
    call AllocateDSI(R2D_Global, LCM_Global, 2, 0.0)

    AHMAX = -1.0E32
    AHMIN = 1.0E32
    ADMAX = -1.0E32
    ADMIN = 1.0E32

    AHDXY = ABS(AHD)       ! *** Default value

    if( process_id == master_id )then
      write(*,'(A)')'READING AHMAP.INP'
      open(1,FILE = 'ahmap.inp')

      STR = READSTR(1)

      do LL = 2,LA_Global
        read(1,*,END = 200) LG, ITMP, JTMP, T1, T2

        R2D_Global(LG,1) = T1
        R2D_Global(LG,2) = T2
      enddo
200   close(1)
    endif

    call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        if( R2D_Global(LG,1) < 0. )then
          AHOXY(L) = ABS(R2D_Global(LG,1))*DXP(L)*DYP(L)
        else
          AHOXY(L) = R2D_Global(LG,1)
        endif
        AHDXY(L) = R2D_Global(LG,2)

        AHMAX = MAX(AHOXY(L),AHMAX)
        AHMIN = MIN(AHOXY(L),AHMIN)
        ADMAX = MAX(AHDXY(L),ADMAX)
        ADMIN = MIN(AHDXY(L),ADMIN)
      endif
    enddo
    deallocate(R2D_Global)

    AHD = ABS(AHD)
    PRINT '(A,4E12.4)','VARIABLE AHO & AHD USED (MIN,MAX): ',AHMIN,AHMAX,ADMIN,ADMAX
  else
    AHDXY = AHD       ! *** Constant value
  endif

  ! *** VERTICAL VISCOSITY AND DIFFUSIVITY
  ! *** Constant and/or default values
  AVO = ABS(AVO)
  ABO = ABS(ABO)
  AVOXY = ABS(AVO)
  AVBXY = ABS(ABO)

  ! *** ENABLE SPATIALLY VARIABLE VERTICAL VISCOSITY AND DIFFUSIVITY
  if( AVO < 0. )then
    call AllocateDSI(R2D_Global, LCM_Global, 2, 0.0)

    AVMAX = -1.0E32
    AVMIN = 1.0E32
    ABMAX = -1.0E32
    ABMIN = 1.0E32

    if( process_id == master_id )then
      write(*,'(A)')'READING AVMAP.INP'
      open(1,FILE = 'avmap.inp')

      STR = READSTR(1)

      do LL = 2,LA_Global
        Read(1,*,END = 300)  LG, ID, JD, T1, T2

        R2D_Global(LG,1) = T1   ! *** AVOXY
        R2D_Global(LG,2) = T2   ! *** AVBXY
      enddo ! *** end loop over the file
300   close(1)
    endif ! *** end on master

    call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        AVOXY(L) = R2D_Global(LG,1)
        AVBXY(L) = R2D_Global(LG,2)

        AHMAX = MAX(AVOXY(L), AVMAX)
        AHMIN = MIN(AVOXY(L), AVMIN)
        ADMAX = MAX(AVBXY(L), ABMAX)
        ADMIN = MIN(AVBXY(L), ABMIN)
      endif

      PRINT '(A,4E12.4)','VARIABLE AVO & ABO USED (MIN,MAX): ',AHMIN,AHMAX,ADMIN,ADMAX
    enddo
    deallocate(R2D_Global)
  endif

  ! *** OPEN FILE CHANSEC.INP FOR 1-D CHANNEL CROSS SECTION DATA
  ! *** REMOVED 2004-09-19  PMC

  ! *** OPEN FILE GWATER.INP TO SPECIFY GROUNDWATER INTERACTION
  ! *** BY INFILTRATION AND EVAPOTRANSPIRATION
  ISGWIE = 0
  RNPOR = 1.E-12
  if( ISGWIT == 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING GWATER.INP'
      open(1,FILE = 'gwater.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      read(1,*) ISGWIE
      if( ISGWIE >= 1 )then
        read(1,*) DAGWZ,RNPOR,RIFTRM
      else
        DAGWZ = 0.0
        RNPOR = 1.E-12
        RIFTRM = 0.0
      endif
      close(1)
    endif

    call Broadcast_Scalar(DAGWZ, master_id)
    call Broadcast_Scalar(RNPOR, master_id)
    call Broadcast_Scalar(RIFTRM,master_id)

  endif

339 FORMAT(2I5,6F14.5)

  ! *** OPEN FILE FBODY.INP TO READ IN SPATIALLY VARYING BODY FORCES
  if( ISBODYF >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING FBODY.INP'
      open(1,FILE = 'fbody.inp',STATUS = 'UNKNOWN')
      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      read(1,*)CVTFACX,CVTFACY
    endif

    call Broadcast_Scalar(CVTFACX, master_id)
    call Broadcast_Scalar(CVTFACY, master_id)

    do LL = 2,LA_Global
      if( process_id == master_id )then
        read(1,*)ITMP,JTMP,FBODY1,FBODY2
      endif

      call Broadcast_Scalar(ITMP,   master_id)
      call Broadcast_Scalar(JTMP,   master_id)
      call Broadcast_Scalar(FBODY1, master_id)
      call Broadcast_Scalar(FBODY2, master_id)

      ! *** Get local i,j
      ID = IG2IL(ITMP)
      JD = JG2JL(JTMP)

      if(id > 0 .and. id <= ic )then
        if(jd > 0 .and. jd <= jc )then

          L = LIJ(id, jd)
          if( ISBODYF == 1 )then
            do K = 1,KC
              FBODYFX(L,K) = CVTFACX*FBODY1
              FBODYFY(L,K) = CVTFACY*FBODY2
            enddo
          endif
          if( ISBODYF == 2 )then
            do K = 1,KC-1
              FBODYFX(L,K) = 0.0
              FBODYFY(L,K) = 0.0
            enddo
            FBODYFX(L,KC) = CVTFACX*FBODY1
            FBODYFY(L,KC) = CVTFACY*FBODY2
          endif
        endif
      endif
    enddo
    do K = 1,KC
      FBODYFX(1,K) = 0.0
      FBODYFY(1,K) = 0.0
      FBODYFX(LC,K) = 0.0
      FBODYFY(LC,K) = 0.0
    enddo

    if( process_id == master_id )then
      close(1)
    endif

  endif

  ! *** Read shellfish farm configuration
  if( ISFFARM > 0 .and. NSF > 0 )then
    call READ_SHELLFISH_JSON(LCM,KC)
  endif

  ! *** OPEN FILE GWMAP.INP TO SPECIFY GROUNDWATER INTERACTION BY AMBIENT GROUNDWATER FLOW
  if( ISGWIT > 1 )then
    ISGWIE = 0

    if( ISGWIT == 3 )then
      if( process_id == master_id )then
        write(*,'(A)')'READING GWSEEP.INP'
        open(1,FILE = 'gwseep.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)
        read(1,*)NSEEPCLASSES
        do M = 1,NSEEPCLASSES
          read(1,*)IDUM,SEEPRATE(M)
        enddo
        close(1)
      endif

      call Broadcast_Scalar(NSEEPCLASSES, master_id)
      call Broadcast_Array(SEEPRATE,      master_id)
    endif

    ! *** Read GWMAP
    allocate(I2D_Global(LCM_Global,3), R1D_Global(LCM_Global))
    I2D_Global = 0
    R1D_Global = 0.

    if( process_id == master_id )then
      write(*,'(A)')'READING GWMAP.INP'
      open(1,FILE = 'gwmap.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      do LL = 2,LA_Global
        !READ(1,*) ITMP,JTMP,IZONE,RVALUE
        read(1,*) (I2D_Global(LL,K),K = 1,3), R1D_Global(LL)
      enddo
      close(1)
    endif

    call Broadcast_Array(I2D_Global, master_id)
    call Broadcast_Array(R1D_Global, master_id)

    ! *** Map to Local Domain
    do LL = 2,LA_GLOBAL
      ID = I2D_Global(LL,1)
      JD = I2D_Global(LL,2)
      LG = LIJ_Global(ID,JD)
      L = Map2Local(LG).LL
      if( L > 1 )then
        IZONE  = I2D_Global(LL,3)
        RVALUE = R1D_Global(LL)
        if( ISGWIT == 3 )then
          if( IZONE > NSEEPCLASSES )then

            write(mpi_log_unit,*)'BAD SEEPAGE CLASS AT I,J = ',id,jd
            call STOPP('')
          endif
          QGW(L) = SEEPRATE(IZONE)*RVALUE
        else
          NGWSL(L) = IZONE
          GWFAC(L) = RVALUE
        endif
      endif
    enddo

    deallocate(I2D_Global)
    deallocate(R1D_Global)
  endif

  ! *** READ IN SPATIALLY VARYING SEDIMENT ROUGHNESS HEIGHT FOR
  ! *** DETERMINING GRAIN STRESS
  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
    if( ISBEDSTR == 3 )then
      call AllocateDSI(R1D_Global, LCM_Global, 0.0)

      if( process_id == master_id )then
        write(*,'(A)')'READING SEDROUGH.INP'
        open(1,FILE = 'sedrough.inp')
        STR = READSTR(1)

        do LG = 2, LA_Global ! @todo make sure this is global
          read(1,*) LDUM, IDUM, JDUM, R1D_Global(LG)                    ! *** ZBRSED
          if( R1D_Global(LG) <= 0.0 )then
            STOP ' BAD SEDIMENT ROUGHNESS IN SEDROUGH.INP'
          endif
        enddo
        close(1)
      endif

      call Broadcast_Array(R1D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          ZBRSED(L) = R1D_Global(LG)
        endif
      enddo

      deallocate(R1D_Global)

    endif

  endif

  ! *** Spatially varying partition coefficients (Water column)
  if( ISTRAN(5) > 0 .and. NSEDS > 0 .and. NTOX > 0 )then
    ! *** Read in spatially varying Kd's/Koc's if any of the sediments or toxics are flagged
    if( ANY(TOXPARW < 0.0) )then
      TOXPARW = ABS(TOXPARW)

      ! *** Spatially variable
      call AllocateDSI(R2D_Global, LCM_Global, NSEDS, 0.0)

      if( process_id == master_id )then
        write(*,'(A)')'READING PARTITIONW.INP'
        open(1,FILE = 'partitionw.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)                                                           ! *** Read the header
      endif

      do NT = 1,NTOX
        if( process_id == master_id )then
          read(1,*) TEXT                                                           ! *** Toxic ID
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(L,NS),NS = 1,NSEDS)    ! *** TOXPARW
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PARTITIONW.INP')
          enddo
        endif
        call Broadcast_Array(R2D_Global, master_id)

        ! *** Map to Local Domain
        do NS = 1,NSEDS
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              TOXPARW(L,NS,NT) = ABS(R2D_Global(LG,NS))
            endif
          enddo
        enddo
      enddo

      if( process_id == master_id )then
        close(1)
      endif
      deallocate(R2D_Global)
    else
      ! *** Spatially constant
      do NT = 1,NTOX
        do NS = 1,NSEDS
          TOXPARW(2:LA,NS,NT) = ABS(TOXPARW(1,NS,NT))
        enddo
      enddo
    endif
  endif

  ! *** Spatially varying partition coefficients (Bed)
  if( ISTRAN(5) > 0 .and. NSEDS > 0 .and. NTOX > 0 )then
    ! *** Read in spatially varying Kd's/Koc's if any of the sediments or toxics are flagged
    if( ANY(TOXPARB < 0.0) )then
      TOXPARB = ABS(TOXPARB)

      ! *** Spatially variable
      call AllocateDSI(R2D_Global, LCM_Global, NSEDS, 0.0)

      if( process_id == master_id )then
        write(*,'(A)')'READING PARTITIONB.INP'
        open(1,FILE = 'partitionb.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)
      endif

      do NT = 1,NTOX
        if( process_id == master_id )then
          read(1,*) TEXT                                                           ! *** Toxic ID
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(L,NS),NS = 1,NSEDS)    ! *** TOXPARB
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PARTITIONB.INP')
          enddo
        endif

        call Broadcast_Array(R2D_Global, master_id)

        ! *** Map to Local Domain
        do NS = 1,NSEDS
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              TOXPARB(L,NS,NT) = R2D_Global(LG,NS)
            endif
          enddo
        enddo
      enddo

      if( process_id == master_id )then
        close(1)
      endif

      deallocate(R2D_Global)

    else
      ! *** Spatially constant
      do NT = 1,NTOX
        do NS = 1,NSEDS
          TOXPARW(2:LA,NS,NT) = ABS(TOXPARW(1,NS,NT))
        enddo
      enddo
    endif
  endif

  ! *** OPEN FILE DOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** DISSOLVED ORGANIC CARBON IN WATER COLUMN
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) == 1 .or. ISTOC(NT) == 2 ) IVAL = 1
  enddo

  if( IVAL == 1 )then
    if( ISTDOCW == 1 )then
      call AllocateDSI(R2D_Global, LCM_Global, KCM, 0.0)

      if( process_id == master_id )then
        write(*,'(A)')'READING DOCW.INP'
        open(1,FILE = 'docw.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)
        read(1,*) ISALTYP

        if( ISALTYP == 0 )then
          do LG = 2,LA_Global
            read(1,*,IOSTAT = ISO) (R2D_Global(LG,K),K = 1,KC)                     ! *** STDOCW
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCW.INP')
          enddo
        else
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(L,K),K = 1,KC)    ! *** STDOCW
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCW.INP')
          enddo
        endif
        close(1)
      endif

      call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          STDOCW(L,:) = R2D_Global(LG,:)
        endif
      enddo

      deallocate(R2D_Global)
    endif

  endif

  ! *** OPEN FILE POCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON IN WATER COLUMN
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) == 1 ) IVAL = 1
  enddo
  if( IVAL == 1 )then
    if( ISTPOCW == 1 )then
      call AllocateDSI(R2D_Global, LCM_Global, KCM, 0.0)

      if( process_id == master_id )then
        write(*,'(A)')'READING POCW.INP'
        open(1,FILE = 'pocw.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)
        read(1,*) ISALTYP

        if( ISALTYP == 0 )then
          do LG = 2,LA_Global
            read(1,*,IOSTAT = ISO) (R2D_Global(LG,K),K = 1,KC)                    ! *** STPOCW
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCW.INP')
          enddo
        else
          do LG = 2,LA_Global
            read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(LG,K),K = 1,KC)  ! *** STPOCW
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCW.INP')
          enddo
        endif
        close(1)
      endif

      call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          STPOCW(L,:) = R2D_Global(LG,:)
        endif
      enddo
      deallocate(R2D_Global)
    endif

    ! *** PREPROCESS STPOCW TO PREVENT ANY DIVISION BY ZERO
    do L = 2,LC-1
      do K = 1,KC
        if( STPOCW(L,K) <= 0.0 ) STPOCW(L,K) = 1E-12
      enddo
    enddo

  endif

  ! *** OPEN FILE FPOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS
  ! *** IN WATER COLUMN
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) == 2 .or. ISTOC(NT) == 3 ) IVAL = 1
  enddo

  if( IVAL == 1 )then
    if( ISTPOCW == 3 )then
      NX = NSED + NSND
      allocate(R3D_Global(LCM_Global,KCM,NX))
      R3D_Global = 0.0

      if( process_id == master_id )then
        write(*,'(A)')'READING FPOCW.INP'
        open(1,FILE = 'fpocw.inp',STATUS = 'UNKNOWN')

        do NS = 1,NX
          STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
          read(1,*) ISALTYP

          if( ISALTYP == 0 )then
            do LG = 2,LA_Global
              read(1,*,IOSTAT = ISO) (R3D_Global(LG,K,NS),K = 1,KC)                                    ! *** STFPOCW
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCW.INP')
            enddo
          else
            do LG = 2, LA_Global
              read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(R3D_Global(LG,K,NS),K = 1,KC)                      ! *** STFPOCW
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCW.INP')
            enddo
          endif
        enddo
        close(1)
      endif

      call Broadcast_Array(R3D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          STFPOCW(L,:,1:NX) = R3D_Global(LG,:,1:NX)
        endif
      enddo
      deallocate(R3D_Global)
    endif
  endif

  ! *** OPEN FILE DOCB.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** DISSOLVED ORGANIC CARBON IN SEDIMENT BED
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) == 1 .or. ISTOC(NT) == 2 ) IVAL = 1
  enddo
  if( IVAL == 1 )then
    if( ISTDOCB == 1 )then
      call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)
      STDOCB(:,:) = 0.0

      if( process_id == master_id )then
        write(*,'(A)')'READING DOCB.INP'
        open(1,FILE = 'docb.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)
        read(1,*) ISALTYP, IREAD, KBINPUT

        if( IREAD == 0 )then
          if( ISALTYP == 0 )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) R2D_Global(L,1)                                ! *** STDOCB
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              do K = 2,KB
                R2D_Global(L,K) = R2D_Global(L,1)
              enddo
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, R2D_Global(L,1)                    ! *** STDOCB
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              do K = 2,KB
                R2D_Global(L,K) = R2D_Global(L,1)
              enddo
            enddo
          endif
        endif

        if( IREAD == 1 )then
          if( ISALTYP == 0 )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) (R2D_Global(L,K),K = 1,KB)                        ! *** STDOCB
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(LDUM,K),K = 1,KB)   ! *** STDOCB
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
            enddo
          endif
        endif

        if( IREAD == 2 )then
          if( ISALTYP == 0 )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) (R2D_Global(L,K),K = 1,KBINPUT)                   ! *** STDOCB
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              do K = KBINPUT,KB
                R2D_Global(L,K) = R2D_Global(L,KBINPUT)
              enddo
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(L,K),K = 1,KBINPUT) ! *** STDOCB
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              do K = KBINPUT,KB
                R2D_Global(L,K) = R2D_Global(L,KBINPUT)
              enddo
            enddo
          endif
        endif

        close(1)
      endif

      call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          STDOCB(L,:) = R2D_Global(LG,:)
        endif
      enddo
      deallocate(R2D_Global)
    endif
  endif

  ! *** OPEN FILE POCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON IN BED
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) == 1 ) IVAL = 1
  enddo
  if( IVAL == 1 )then
    if( ISTPOCB == 1 )then
      if( process_id == master_id )then
        write(*,'(A)')'READING POCB.INP'
        open(1,FILE = 'pocb.inp',STATUS = 'UNKNOWN')
        STR = READSTR(1)
        read(1,*) ISALTYP, IREAD, KBINPUT
      endif

      call Broadcast_Scalar(ISALTYP, master_id)
      call Broadcast_Scalar(IREAD, master_id)
      call Broadcast_Scalar(KBINPUT, master_id)


      STPOCB(L,K) = 0.0

      if( IREAD == 0 )then
        ! ***
        if( ISALTYP == 0 )then
          if( num_Processors > 1 )then
            write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
            pause
          endif

          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) STPOCB(L,1)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            do K = 2,KB
              STPOCB(L,K) = STPOCB(L,1)
            enddo
          enddo
        else
          if( num_Processors > 1 )then
            write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
            pause
          endif
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,STPOCB(L,1)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            do K = 2,KB
              STPOCB(L,K) = STPOCB(L,1)
            enddo
          enddo
        endif
      endif

      if( IREAD == 1 )then
        ! ***
        if( ISALTYP == 0 )then
          if( num_Processors > 1 )then
            write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
            pause
          endif
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) (STPOCB(L,K),K = 1,KB)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
          enddo
        else
          ! *** delme - this needs to be more general
          call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)
          STDOCB(:,:) = 0.0

          if( process_id == master_id )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM,IDUM,JDUM,(STPOCB(1,K),K = 1,KB)
              R2D_Global(LDUM,:) = STPOCB(1,:)
              if( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
            enddo
            STPOCB(1,:) = 0.0
          endif

          call Broadcast_Array(R2D_Global, master_id)

          ! *** Map to Local Domain
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              STPOCB(L,:) = R2D_Global(LG,:)
            endif
          enddo
          deallocate(R2D_Global)
        endif
      endif

      if( IREAD == 2 )then
        ! ***
        if( num_Processors > 1 )then
          write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
          pause
        endif

        if( ISALTYP == 0 )then
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) (STPOCB(L,K),K = 1,KBINPUT)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            do K = KBINPUT,KB
              STPOCB(L,K) = STPOCB(L,KBINPUT)
            enddo
          enddo
        else
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(STPOCB(L,K),K = 1,KBINPUT)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            do K = KBINPUT,KB
              STPOCB(L,K) = STPOCB(L,KBINPUT)
            enddo
          enddo
        endif
      endif

      if( process_id == master_id )then
        close(1)
      endif

    endif
  endif

  ! *** OPEN FILE FPOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS
  ! *** IN BED
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) == 2 .or. ISTOC(NT) == 3) IVAL = 1
  enddo
  if( IVAL == 1 )then
    if( ISTPOCB == 3 )then
      if( process_id == master_id )then
        write(*,'(A)')'READING FPOCB.INP'
        open(1,FILE = 'fpocb.inp',STATUS = 'UNKNOWN')
        do NS = 1,NSED+NSND
          STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
          read(1,*) ISALTYP,IREAD,KBINPUT

          if( IREAD == 0 )then
            if( ISALTYP == 0 )then
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO) STFPOCB(L,1,NS)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                do K = 2,KB
                  STFPOCB(L,K,NS) = STFPOCB(L,1,NS)
                enddo
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,STFPOCB(L,1,NS)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                do K = 2,KB
                  STFPOCB(L,K,NS) = STFPOCB(L,1,NS)
                enddo
              enddo
            endif
          endif
          if( IREAD == 1 )then
            if( ISALTYP == 0 )then
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO) (STFPOCB(L,K,NS),K = 1,KB)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(STFPOCB(L,K,NS),K = 1,KB)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
              enddo
            endif
          endif
          if( IREAD == 2 )then
            if( ISALTYP == 0 )then
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO) (STFPOCB(L,K,NS),K = 1,KBINPUT)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                do K = KBINPUT,KB
                  STFPOCB(L,K,NS) = STFPOCB(L,KBINPUT,NS)
                enddo
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(STFPOCB(L,K,NS),K = 1,KBINPUT,NS)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                do K = KBINPUT,KB
                  STFPOCB(L,K,NS) = STFPOCB(L,KBINPUT,NS)
                enddo
              enddo
            endif
          endif
        enddo
        close(1)
      endif

      call Broadcast_Array(STFPOCB, master_id)

    endif
  endif

  ! *** *******************************************************************C
  !###########################################################################
  ! HQI Change to include sptially varying, but time constant bulk foc
  ! FPOCB  - Bulk foc from data
  ! PFPOCB - Pseudo foc from data, to be used for all partitioning calculations
  ! RM, 02/29/04
  ! *** *******************************************************************C

  ! *** OPEN FILE FOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON IN BED AND PSEUDO-POC IN BED
  if( ISTPOCB == 4 )then
    call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)

    do K = 1,KB
      do L = 2,LA
        FPOCB(L,K) = 0.0
      enddo
    enddo

    if( process_id == master_id )then
      write(*,'(A)')'READING FOCB.INP'
      open(1,FILE = 'focb.inp',STATUS = 'UNKNOWN')

      STR = READSTR(1)
      read(1,*) ISALTYP, IREAD, KBINPUT

      do LG = 2,La_Global
        read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(LG,K),K = 1,KBINPUT)    ! *** FPOCB
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FOCB.INP')
      enddo
      close(1)
    endif

    call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        do K = 1,KBINPUT
          FPOCB(L,K) = R2D_Global(LG,K)/1000000.
        enddo
        do K = KBINPUT+1,KB
          FPOCB(L,K) = FPOCB(L,KBINPUT)
        enddo
      endif
    enddo

    R2D_Global = 0.0
    if( process_id == master_id )then
      write(*,'(A)')'READING PSEUDO_FOCB.INP'
      open(1,FILE = 'pseudo_focb.inp',STATUS = 'UNKNOWN')

      STR = READSTR(1)
      read(1,*) ISALTYP, IREAD, KBINPUT

      do LG = 2,LA_Global
        read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (R2D_Global(LG,K),K = 1,KBINPUT)     ! *** PFPOCB
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSEUDO_FOCB.INP')
      enddo
      close(1)
    endif

    call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        do K = 1,KBINPUT
          PFPOCB(L,K) = R2D_Global(LG,K)/1000000.
        enddo
        do K = KBINPUT+1,KB
          PFPOCB(L,K) = PFPOCB(L,KBINPUT)
        enddo
      endif
    enddo
    close(1)

    deallocate(R2D_Global)

  endif

  if( ISVHEAT > 0 .and. .not. (ISTOPT(2) == 1 .and. IASWRAD == 3) ) ISVHEAT = 0

  ! *** Assigning background light extinction to SWRATNF when using fast/slow Extinction coefficients option
  do L = 2,LA
    if( ISTRAN(8) > 0 )then
      SVKEBACK(L) = WQKEB(1)
    else
      SVKEBACK(L) = SWRATNF
    endif
  enddo

  if( ISTRAN(2) > 0 .and. IASWRAD == 3 .and. ISVHEAT > 0 )then
    call AllocateDSI(R2D_Global, LCM_Global, 3, 0.0)

    ! *** DEFAULT VALUES
    do L = 2,LA
      SVREVC(L) = 0.001*ABS(REVC)
      SVRCHC(L) = 0.001*ABS(RCHC)
      if( REVC < 0. ) LSVHTWINDE(L) = .TRUE.
      if( RCHC < 0. ) LSVHTWINDC(L) = .TRUE.
    enddo

    if( process_id == master_id )then
      write(*,'(A)')'READING SVHTFACT.INP'
      open(1,FILE = 'svhtfact.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      ! *** LOOP OVER THE DOMAIN AND ONLY SET THE CELLS THAT ARE IN THE FILE
      do LG = 2,LA_Global
        read(1,*,END = 101) LL, ITMP, JTMP, R2D_Global(LG,1), R2D_Global(LG,2), R2D_Global(LG,3)
      enddo
101   close(1)
    endif

    call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        LSVHTWINDE(L) = R2D_Global(LG,1) < 0.
        LSVHTWINDC(L) = R2D_Global(LG,2) < 0.
        SVREVC(L)     = 0.001*ABS(R2D_Global(LG,1))
        SVRCHC(L)     = 0.001*ABS(R2D_Global(LG,2))
        SVKEBACK(L)   = R2D_Global(LG,3)
      endif
    enddo

    deallocate(R2D_Global)

    call Broadcast_Array(LSVHTWINDE , master_id)
    call Broadcast_Array(LSVHTWINDC , master_id)
    call Broadcast_Array(SVREVC     , master_id)
    call Broadcast_Array(SVRCHC     , master_id)
    call Broadcast_Array(SVKEBACK   , master_id)

  endif

  ! *** READ IN INITIAL SALINITY, TEMPERATURE, DYE, SED, SND, TOX
  ! *** FOR COLD STARTS FORM FILE XXXX.INP
  ! *** SALINITY
  if( ISTRAN(1) >= 1 )then
      do K = 1,KC
        do L = 2,LA
          SALINIT(L,K) = 0.
        enddo
      enddo

      if( ISRESTI == 0 .or. (ISRESTI >= 1 .and. ISCI(1) == 0) .or. (ISTOPT(1) > 1))  then
        if( ISTOPT(1) >= 1 )then
          if( process_id == master_id )then
            write(*,'(A)')'READING SALT.INP'
            open(1,FILE = 'salt.inp',STATUS = 'UNKNOWN')

            ! ***   SKIP OVER TITLE AND AND HEADER LINES
            STR = READSTR(1)
            read(1,*)ISALTYP

            if( ISALTYP == 0 )then
              do LG = 2, LA_GLOBAL
                read(1,*,IOSTAT = ISO) (SAL_Global(LG,K),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SALT.INP')
              enddo
            else
              do L = 2, LA_GLOBAL
                read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(CONINIT(K),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SALT.INP')

                LG = LIJ_Global(IDUM,JDUM)
                SAL_Global(LG,:) = CONINIT(:)
              enddo
            endif
            close(1)
          endif

          call Broadcast_Array(SAL_Global, master_id)

          ! *** Map to Local Domain
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              SALINIT(L,:) = SAL_Global(LG,:)
            endif
          enddo
        endif
      endif

  endif

  ! *** TEMPERATURE
  if( ISTRAN(2) >= 1 )then
      TEM_Global = TEMO
      do K = 1,KC
        do L = 2,LA
          TEMINIT(L,K) = TEMO
        enddo
      enddo

      if( ISRESTI == 0 .or. (ISRESTI >= 1 .and. ISCI(2) == 0) .or. (ISTOPT(2) > 9) ) then
        if( ISTOPT(2) >= 1 .or. INITTEMP > 0 )then
          if( process_id == master_id )then
            write(*,'(A)')'READING TEMP.INP'
            open(1,FILE = 'temp.inp',STATUS = 'UNKNOWN')
            ! ***   SKIP OVER TITLE AND AND HEADER LINES
            STR = READSTR(1)
            read(1,*) ISALTYP

            if( ISALTYP == 0 )then
              do LG = 2,LA_Global
                read(1,*,IOSTAT = ISO) (TEM_Global(LG,K),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TEMP.INP')
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONINIT(K),K = 1,KC)

                LG = LIJ_Global(IDUM,JDUM)
                TEM_Global(LG,:) = CONINIT(:)
              enddo
            endif
            close(1)
          endif
          call Broadcast_Array(TEM_Global, master_id)

          ! *** Map to Local Domain
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              TEMINIT(L,:) = TEM_Global(LG,:)
            endif
          enddo
        endif
      endif
  endif

  ! *** DYE
  if( ISTRAN(3) >= 1 )then
    IFLAG = SUM(DYES(:).ICFLAG)                                                      ! *** Spatially varying flag
    BFLAG = ISRESTI == 0 .or. (ISRESTI == 1 .and. ISCI(3) == 0) .or. ISTOPT(3) > 1   ! *** Cold start flag
    if( BFLAG )then
      if( IFLAG > 0 )then
        if( process_id == master_id )then
          write(*,'(A)')'READING DYE.INP'
          open(1,FILE = 'dye.inp',STATUS = 'UNKNOWN')

          do MD = 1,NDYE
            ! ***   SKIP OVER TITLE AND AND HEADER LINES
            STR = READSTR(1)
            read(1,*)ISALTYP

            if( ISALTYP == 0 )then
              do L = LG,LA_Global
                read(1,*,IOSTAT = ISO) (DYE_Global(LG,K,MD),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DYE.INP')
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(CONINIT(K),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DYE.INP')

                LG = LIJ_Global(IDUM,JDUM)
                DYE_Global(LG,:,MD) = CONINIT(:)
              enddo
            endif
          enddo
          close(1)
        endif
        call Broadcast_Array(DYE_Global, master_id)

        ! *** Map to Local Domain
        do MD = 1,NDYE
          if( BFLAG .and. DYES(MD).ICFLAG > 0 )then
            ! *** Update IC for this class
            do LG = 2,LA_GLOBAL
              L = Map2Local(LG).LL
              if( L > 1 )then
                DYEINIT(L,:,MD) = DYE_Global(LG,:,MD)
              endif
            enddo
          else
            ! *** Spatially Constant
            DYEINIT(:,:,MD) = DYES(MD).IC
          endif
        enddo
      else
        ! *** Spatially Constant
        do MD = 1,NDYE
          DYEINIT(:,:,MD) = DYES(MD).IC
        enddo
      endif

    endif

  endif    ! *** END OF DYE.INP LOADING

  ! *** SFL
  if( ISRESTI == 0 .and. ISTRAN(4) >= 1 )then
    if( ISTOPT(4) >= 1 )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SFL.INP'
        open(1,FILE = 'sfl.inp',STATUS = 'UNKNOWN')

        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        read(1,*) ISALTYP

        if( ISALTYP == 0 )then
          do LG = 2,LA_Global
            read(1,*,IOSTAT = ISO) (SFL_Global(LG,K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFL.INP')
          enddo
        else
          do L = 2,LA_Global
            read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONINIT(K),K = 1,KC)

            LG = LIJ_Global(IDUM,JDUM)
            SFL_Global(LG,:) = CONINIT(:)
          enddo
        endif
        close(1)
      endif

      call Broadcast_Array(SFL_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          SFLINIT(L,:) = SFL_Global(LG,:)
        endif
      enddo

    endif
  endif

  ! ****************************************************************************************
  ! *** SEDIMENTS AND SEDIMENT BED FILES

  ! *** INITIALIZE HARD-BOTTOM, IF USED
  LASED = 0
  LDMSED = 0
  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    allocate(BEDMAP(LCM))
    allocate(LSED(LCM))
    allocate(LBED(LCM))

    BEDMAP = 1                ! *** BED PROCESSES ACTIVE
    LBED = .FALSE.
    do L = 2,LA
      LASED = LASED + 1
      LSED(LASED) = L
    enddo

    if( ISBEDMAP > 0 )then
      ! *** READ USER SPECIFIED SEDIMENT ACTIVE CELL LIST
      if( process_id == master_id )then
        write(*,'(A)')'READING BEDMAP.INP'
        open(1,FILE = 'bedmap.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
        do L = 2,LA_GLOBAL
          read(1,*,IOSTAT = ISO,END = 400) IDUM,JDUM,NDUM1
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDMAP.INP')

          LG = LIJ_Global(IDUM,JDUM)
          BEDMAP_Global(LG) = NDUM1
        enddo
400     close(1)
      endif
      call Broadcast_Array(BEDMAP_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          BEDMAP(L) = BEDMAP_Global(LG)
        endif
      enddo

      LASED = 0
      LSED = 0
      LBED = .TRUE.
      do L = 2,LA
        if( BEDMAP(L) > 0 )then
          LASED = LASED + 1
          LSED(LASED) = L
          LBED(L) = .FALSE.
        endif
      enddo
    endif

    LDMSED = INT(LASED/NTHREADS)+1

    ! *** OPEN FILE SEDBLBC.INP TO READ IN SEDIMENT BEDLOAD OUTFLOW
    ! *** OR RECIRCULATION BOUNDARY CONDITIONS
    NSBDLDBC = 0
    if( (ICALC_BL == 1 .and. .not. LSEDZLJ) .or. (ICALC_BL >= 2 .and. LSEDZLJ) )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SEDBLBC.INP'
        open(1,FILE = 'sedblbc.inp',STATUS = 'UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        read(1,*) NDUM

        allocate(IBLTMP(5,NDUM))
        IBLTMP = 0

        do NS = 1,NDUM
          read(1,*) (IBLTMP(I,NS),I = 1,5)
        enddo
        close(1)
      endif
      call Broadcast_Scalar(NDUM, master_id)

      if( process_id /= master_id )then
        allocate(IBLTMP(5,NDUM))
        IBLTMP = 0
      endif
      call Broadcast_Array(IBLTMP, master_id)

      do NS = 1,NDUM
        ! *** Check if in current node
        LG = LIJ_Global(IBLTMP(1,NS),IBLTMP(2,NS))
        LU = Map2Local(LG).LL
        if( LU > 0 )then
          ! *** ENSURE US CELLS ARE PART OF THE ACTIVE MAP, OTHERWISE DISREGARD
          do LL= 1,LASED
            if( LSED(LL) == LU ) EXIT
          enddo
          if( LL > LASED ) CYCLE

          NSBDLDBC = NSBDLDBC + 1
          ISDBLDIR(NSBDLDBC) = ABS(IBLTMP(5,NS))
          LSBLBCU(NSBDLDBC) = LU

          if( IBLTMP(3,NS) > 0 .and. IBLTMP(4,NS) > 0 )then
            LG = LIJ_Global(IBLTMP(3,NS),IBLTMP(4,NS))
            LD = Map2Local(LG).LL
            LSBLBCD(NSBDLDBC) = LD
          else
            LSBLBCD(NSBDLDBC) = 0
          endif
        endif
      enddo

      do NS = 1,NSBDLDBC
        ISDBLDIR(NS) = ABS(ISDBLDIR(NS))
        if( ISDBLDIR(NS) == 1 )then
          ! *** EAST-WEST
          if( SUB(LSBLBCU(NS)) < 0.5 .and. SUB(LEC(LSBLBCU(NS))) > 0.5 ) ISDBLDIR(NS) = -1
        endif
        if( ISDBLDIR(NS) == 2 )then
          ! *** NORTH-SOUTH
          if( SVB(LSBLBCU(NS)) < 0.5 .and. SVB(LEC(LSBLBCU(NS))) > 0.5 ) ISDBLDIR(NS) = -2
        endif
      enddo

    endif
  endif

  ! ****************************************************************************************
  ! *** TOXICS
  if( ISTRAN(5) >= 1 )then
    do NT = 1,NTOX
      do K = 1,KC
        do L = 2,LA
          TOXINIT(L,K,NT) = TOXINTW(NT)
        enddo
      enddo
    enddo
    do NT = 1,NTOX
      do K = 1,KB
        do L = 2,LA
          TOXBINIT(L,K,NT) = TOXINTB(NT)
        enddo
      enddo
    enddo

    ! *** SPATIALLY VARYING WATER COLUMN INITIAL CONDITIONS
    ISCOLD = 1
    if( (ISRESTI >= 1 .and. ISCI(5) > 0) .or. ISLTMT > 0 ) ISCOLD = 0
    if( ISCOLD == 1 )then
      if( ITXINT(1) == 2 .or. ITXINT(1) == 3 )then
        if( process_id == master_id )then
          write(*,'(A)')'READING TOXW.INP'
          open(1,FILE = 'toxw.inp',STATUS = 'UNKNOWN')

          do NT = 1,NTOX
            STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
            read(1,*) ISALTYP, ITOXWU(NT)

            if( ISALTYP == 0 )then
              do LG = 2,LA_Global
                read(1,*,IOSTAT = ISO) (TOX_Global(LG,K,NT),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXW.INP')
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONINIT(K),K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXW.INP')

                LG = LIJ_Global(IDUM,JDUM)
                TOX_Global(LG,:,NT) = CONINIT(:)
              enddo
            endif
          enddo
          close(1)
        endif

        call Broadcast_Array(ITOXWU,     master_id)
        call Broadcast_Array(TOX_Global, master_id)

        ! *** Map to Local Domain
        do NT = 1,NTOX
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              TOXINIT(L,:,NT) = TOX_Global(LG,:,NT)
            endif
          enddo
        enddo

        ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
        if( process_id == master_id )then
          write(*,'(A)')'READING TOXB.INP'
          open(1,FILE = 'toxb.inp',STATUS = 'UNKNOWN')

          do NT = 1,NTOX
            STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
            ! *** ITOXBU is not used
            read(1,*) ISALTYP,ITOXBU(NT),KBINPUT

            do K = 1,KB
              do LG = 2,LA_Global
                TOXB_Global(LG,K,NT) = 0.0
              enddo
            enddo

            if( ISALTYP == 0 )then
              do LG = 2,LA_Global
                read(1,*,IOSTAT = ISO) (TOXB_Global(LG,K,NT),K = 1,KBINPUT)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXB.INP')
              enddo
            else
              do L = 2,LA_Global
                read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONBINIT(K),K = 1,KBINPUT)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXB.INP')

                LG = LIJ_Global(IDUM,JDUM)
                TOXB_Global(LG,1:KBINPUT,NT) = CONBINIT(1:KBINPUT)
              enddo
            endif

            ! *** FILL ANY MISSING LAYERS
            do K = KBINPUT+1,KB
              TOXB_Global(LG,K,NT) = TOXB_Global(LG,KBINPUT,NT)
            enddo
          enddo
          close(1)
        endif
        call Broadcast_Array(ITOXBU,      master_id)
        call Broadcast_Array(TOXB_Global, master_id)

        ! *** Map to Local Domain
        do NT = 1,NTOX
          do LG = 2,LA_GLOBAL
            L = Map2Local(LG).LL
            if( L > 1 )then
              TOXBINIT(L,:,NT) = TOXB_Global(LG,:,NT)
            endif
          enddo
        enddo
      endif
    endif
  endif

  ! ****************************************************************************************
  ! *** COHESIVE/SEDZLJ SEDIMENTS
  if( ISTRAN(6)  >= 1 )then
    do NS = 1,NSED
      do K = 1,KC
        do L = 2,LA
          SEDINIT(L,K,NS) = SEDO(NS)
        enddo
      enddo
    enddo
    do NS = 1,NSED
      do K = 1,KB
        do L = 2,LA
          SEDBINIT(L,K,NS) = SEDBO(NS)
        enddo
      enddo
    enddo

    ! *** SPATIALLY VARYING WATER COLUMN INITIAL CONDITIONS
    ISCOLD = 1
    if( (ISRESTI >= 1 .and. ISCI(6) > 0) .or. ISLTMT > 0 ) ISCOLD = 0
    if( ISCOLD == 1 .and. (ISEDINT == 1 .or. ISEDINT == 3) )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SEDW.INP'
        open(1,FILE = 'sedw.inp',STATUS = 'UNKNOWN')

        do NS = 1,NSED
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR = READSTR(1)
          read(1,*)ISALTYP

          if( ISALTYP == 0 )then
            do LG = 2,LA_Global
              read(1,*,IOSTAT = ISO) (SED_Global(LG,K,NS),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDW.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONINIT(K),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDW.INP')

              LG = LIJ_Global(IDUM,JDUM)
              SED_Global(LG,:,NS) = CONINIT(:)
            enddo
          endif
        enddo
        close(1)
      endif
      call Broadcast_Array(SED_Global, master_id)

      ! *** Map to Local Domain
      do NS = 1,NSED
        do LG = 2,LA_GLOBAL
          L = Map2Local(LG).LL
          if( L > 1 )then
            SEDINIT(L,:,NS) = SED_Global(LG,:,NS)
          endif
        enddo
      enddo
    endif

    ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
    if( NSEDFLUME == 0 )then
      ISCOLD = 1
      if( (ISRESTI >= 1 .and. ISCI(6) > 0) .or. ISLTMT > 0 ) ISCOLD = 0
    elseif( NSEDFLUME == 3 .and. IHTSTRT == 0 )then
      ISCOLD = 1
    else
      ISCOLD = 0
    endif

    if( ISCOLD == 1 .and. (ISEDINT == 2 .or. ISEDINT == 3) )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SEDB.INP'
        open(1,FILE = 'sedb.inp',STATUS = 'UNKNOWN')

        do NS = 1,NSED
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR = READSTR(1)
          read(1,*) ISALTYP, ISEDBU(NS), KBINPUT

          do K = 1,KB
            do LG = 2,LA_Global
              SEDB_Global(LG,K,NS) = 0.0
            enddo
          enddo

          if( ISALTYP == 0 )then
            do LG = 2,LA_Global
              read(1,*,IOSTAT = ISO) (SEDB_Global(LG,K,NS),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDB.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONBINIT(K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDB.INP')

              LG = LIJ_Global(IDUM,JDUM)
              SEDB_Global(LG,1:KBINPUT,NS) = CONBINIT(1:KBINPUT)
            enddo
          endif
        enddo
        close(1)
      endif
      call Broadcast_Array(ISEDBU,      master_id)
      call Broadcast_Array(SEDB_Global, master_id)

      ! *** Map to Local Domain
      do NS = 1,NSED
        do LG = 2,LA_GLOBAL
          L = Map2Local(LG).LL
          if( L > 1 )then
            SEDBINIT(L,:,NS) = SEDB_Global(LG,:,NS)
          endif
        enddo
      enddo
    endif
  endif   ! *** END OF NON-COHESIVE CONFIGURATION

  ! ****************************************************************************************
  ! *** NON-COHESIVE SEDIMENTS
  if( ISTRAN(7) >= 1 )then
    do NX = 1,NSND
      NS = NSED + NX
      do K = 1,KC
        do L = 2,LA
          SNDINIT(L,K,NX) = SEDO(NS)
        enddo
      enddo
    enddo
    do NX = 1,NSND
      NS = NSED + NX
      do K = 1,KB
        do L = 2,LA
          SNDBINIT(L,K,NX) = SEDBO(NS)
        enddo
      enddo
    enddo

    ! *** SPATIALLY VARYING WATER COLUMN INITIAL CONDITIONS
    ISCOLD = 1
    if( (ISRESTI >= 1 .and. ISCI(7) > 0) .or. ISLTMT > 0 ) ISCOLD = 0
    if( ISCOLD == 1 .and. (ISEDINT == 1 .or. ISEDINT == 3) )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SNDW.INP'
        open(1,FILE = 'sndw.inp',STATUS = 'UNKNOWN')

        do NX = 1,NSND
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR = READSTR(1)
          read(1,*)ISALTYP
          if( ISALTYP == 0 )then
            do LG = 2,LA_Global
              read(1,*,IOSTAT = ISO) (SND_Global(LG,K,NX),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDW.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONINIT(K),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDW.INP')

              LG = LIJ_Global(IDUM,JDUM)
              SND_Global(LG,:,NX) = CONINIT(:)
            enddo
          endif
        enddo
        close(1)
      endif
      call Broadcast_Array(SND_Global, master_id)

      ! *** Map to Local Domain
      do NX = 1,NSND
        do LG = 2,LA_GLOBAL
          L = Map2Local(LG).LL
          if( L > 1 )then
            SNDINIT(L,:,NX) = SND_Global(LG,:,NX)
          endif
        enddo
      enddo

    endif

    ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
    if( NSEDFLUME == 0 )then
      ISCOLD = 1
      if( (ISRESTI >= 1 .and. ISCI(7) > 0) .or. ISLTMT > 0 ) ISCOLD = 0
    else
      ISCOLD = 0
    endif

    if( ISCOLD == 1 .and. (ISEDINT == 2 .or. ISEDINT == 3) )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SNDB.INP'
        open(1,FILE = 'sndb.inp',STATUS = 'UNKNOWN')
        do NX = 1,NSND
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR = READSTR(1)
          read(1,*) ISALTYP,ISNDBU(NX),KBINPUT

          do K = 1,KB
            do LG = 2,LC_Global
              SNDB_Global(LG,K,NX) = 0.0
            enddo
          enddo

          if( ISALTYP == 0 )then
            do LG = 2,LA_Global
              read(1,*,IOSTAT = ISO) (SNDB_Global(LG,K,NX),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDB.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) LDUM, IDUM, JDUM, (CONBINIT(K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDB.INP')

              LG = LIJ_Global(IDUM,JDUM)
              SNDB_Global(LG,1:KBINPUT,NX) = CONBINIT(1:KBINPUT)
            enddo
          endif
        enddo
        close(1)
      endif

      call Broadcast_Array(ISNDBU,      master_id)
      call Broadcast_Array(SNDB_Global, master_id)

      ! *** Map to Local Domain
      do NX = 1,NSND
        do LG = 2,LA_GLOBAL
          L = Map2Local(LG).LL
          if( L > 1 )then
            SNDBINIT(L,:,NX) = SNDB_Global(LG,:,NX)
          endif
        enddo
      enddo
    endif
  endif   ! *** END OF NON-COHESIVE CONFIGURATION

  ! ****************************************************************************************
  !  ** SEDIMENT BED LAYER CONFIGURATION
  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
    if( .not. LSEDZLJ )then
      ISCOLD = 1
      if( ISLTMT > 0 )then
        ISCOLD = 0
      elseif( ISRESTI >= 1 )then
        if( ISTRAN(6) > 0 .and. ISCI(6) > 0 ) ISCOLD = 0
        if( ISTRAN(7) > 0 .and. ISCI(7) > 0 ) ISCOLD = 0
      endif
    elseif( NSEDFLUME == 3 .and. IHTSTRT == 0 )then
      ISCOLD = 1
    else
      ISCOLD = 0
    endif

    if( ISCOLD == 1 .and. (ISEDINT == 2 .or. ISEDINT == 3) )then
      !  ** BED LAYER THICKNESS
      call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)

      if( process_id == master_id )then
        write(*,'(A)')'READING BEDLAY.INP'
        open(1,FILE = 'bedlay.inp',STATUS = 'UNKNOWN')

        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        read(1,*) IBEDLAYU, ISALTYP, KBINPUT
        if( IBEDLAYU > 0 )then
          do K = 1,KB
            do L = 2,LA_Global
              R2D_Global(L,K) = 0.0
            enddo
          enddo
          if( ISALTYP == 0 )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) (R2D_Global(L,K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDLAY.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(R2D_Global(L,K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDLAY.INP')
            enddo
          endif
        endif
        close(1)
      endif   ! *** End master process

      call Broadcast_Scalar(IBEDLAYU,  master_id)
      call Broadcast_Array(R2D_Global, master_id)

      if( IBEDLAYU == 1 ) TMPCVT = 0.001    ! *** THICKNESS IN MM
      if( IBEDLAYU == 2 ) TMPCVT = 0.01     ! *** THICKNESS IN CM
      if( IBEDLAYU == 3 ) TMPCVT = 1.0      ! *** THICKNESS IN M

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          BEDLINIT(L,:) = TMPCVT*R2D_Global(LG,:)     ! *** THICKNESS IN M
        endif
      enddo

      ! *** BED LAYER BULK DENSITY
      R2D_Global = 0.0
      if( process_id == master_id )then
        write(*,'(A)')'READING BEDBDN.INP'
        open(1,FILE = 'bedbdn.inp',STATUS = 'UNKNOWN')

        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        read(1,*) IBEDBDNU, ISALTYP, KBINPUT

        if( IBEDBDNU > 0 )then
          do K = 1,KB
            do L = 2,LA_Global
              R2D_Global(L,K) = 0.0
            enddo
          enddo
          if( ISALTYP == 0 )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) (R2D_Global(L,K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDBDN.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(R2D_Global(L,K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDBDN.INP')
            enddo
          endif
        endif
        close(1)
      endif   ! *** End master process

      call Broadcast_Scalar(IBEDBDNU,  master_id)
      call Broadcast_Array(R2D_Global, master_id)

      TMPCVT = 0.0                          ! *** INVALID UNITS
      if( IBEDBDNU == 1 ) TMPCVT = 1.0      ! *** BULK DENSITY IN KG/M**3
      if( IBEDBDNU == 2 ) TMPCVT = 1000.0   ! *** BULK DENSITY IN GM/CM**3

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          BEDBINIT(L,:) = TMPCVT*R2D_Global(LG,:)     ! *** BULK DENSITY IN KG/M**3
        endif
      enddo

      ! *** BED VOID RATIO
      R2D_Global = 0.0
      if( process_id == master_id )then
        write(*,'(A)')'READING BEDDDN.INP'
        open(1,FILE = 'bedddn.inp',STATUS = 'UNKNOWN')

        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        read(1,*) IBEDDDNU, ISALTYP, KBINPUT
        if( IBEDDDNU > 0 )then
          do K = 1,KB
            do L = 2,LA_Global
              R2D_Global(L,K) = 0.0
            enddo
          enddo
          if( ISALTYP == 0 )then
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO) (R2D_Global(L,K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDDDN.INP')
            enddo
          else
            do L = 2,LA_Global
              read(1,*,IOSTAT = ISO)LDUM,IDUM,JDUM,(R2D_Global(L,K),K = 1,KBINPUT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDDDN.INP')
            enddo
          endif
        endif
        close(1)
      endif   ! *** End master process

      call Broadcast_Scalar(IBEDDDNU,  master_id)
      call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          BEDDINIT(L,:) = R2D_Global(LG,:)
        endif
      enddo
      deallocate(R2D_Global)

      ! *** INPUT IS DRY DENSITY (G/CM3)
      if( IBEDDDNU <= 1 )then
        do K = 1,KB
          do L = 2,LA
            POR = 0.001*(BEDBINIT(L,K) - BEDDINIT(L,K))
            BEDBINIT(L,K) = POR/(1.0 - POR)
          enddo
        enddo
      endif
      ! *** INPUT IS POROSITY
      if( IBEDDDNU == 2 )then
        do K = 1,KB
          do L = 2,LA
            BEDDINIT(L,K) = BEDDINIT(L,K)/(1. - BEDDINIT(L,K))
          enddo
        enddo
      endif
      ! *** ELSE INPUT IS VOID RATIO (IBEDDDNU == 3)

      ! *** IF MASS FRACTION IS INPUT THEN ENSURE SUM OF 1.0 THEN COMPUTE MASS
      IFLAG = 0
      do NS = 1,NSED
        if( ISEDBINT == 1 .and. ISEDBU(NS) /= 1 ) IFLAG = 1
      enddo
      do NX = 1,NSND
        if( ISEDBINT == 1 .and. ISNDBU(NX) /= 1 ) IFLAG = 1
      enddo

      if( IFLAG == 0 )then
        ! *** CHECK FRACTIONS
        do L = 2,LA
          do K = 1,KB
            T1 = 0.0
            do NS = 1,NSED
              T1 = T1 + SEDBINIT(L,K,NS)
            enddo
            do NX = 1,NSND
              T1 = T1 + SNDBINIT(L,K,NX)
            enddo

            if( ABS(1.0 - T1) > 1E-6 .and. T1 > 0.0 )then
              ! *** FRACTIONS DO NOT ADD UP.  SET FRACTIONS
              do NS = 1,NSED
                SEDBINIT(L,K,NS) = SEDBINIT(L,K,NS)/T1
              enddo
              do NX = 1,NSND
                SNDBINIT(L,K,NX) = SNDBINIT(L,K,NX)/T1
              enddo
            endif
          enddo
        enddo
      endif

    elseif( ISCOLD == 1 .and. ISEDINT < 2 )then

      ! *** LAYER THICKNESSES  (M)
      BEDLINIT = 0.0
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do L = 2,LA
              BEDLINIT(L,K) = BEDLINIT(L,K) + SDEN(NS)*SEDB(L,K,NS)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED+NX
          do K = 1,KB
            do L = 2,LA
              BEDLINIT(L,K) = BEDLINIT(L,K)+SDEN(NS)*SNDB(L,K,NX)
            enddo
          enddo
        enddo
      endif
    endif    ! *** END OF COLD START

  endif      ! *** END OF SEDIMENT BED CONFIGURATION

  ! ****************************************************************************************************************************
  ! *** External read forcings

  ! *** Read in open boundary surface elevation time series from the file PSER.INP
  if( NPSER >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING PSER.INP'
      open(1,FILE = 'pser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      do NS = 1,NPSER
        read(1,*,IOSTAT = ISO) ITYPE, TSPS(NS).NREC, TSPS(NS).TMULT, TSPS(NS).TOFFSET, RMULADJ, ADDADJ, PSERZDF(NS), INTPSER(NS)
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
        PSERZDF(NS) = G*PSERZDF(NS)
        if( ITYPE == 1 )then
          read(1,*,IOSTAT = ISO) RMULADJS,ADDADJS,PSERZDS(NS)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
          PSERZDS(NS) = G*PSERZDS(NS)
        else
          RMULADJS = 0
          ADDADJS = 0.
          PSERZDS(NS) = 0.
        endif
        if( ITYPE == 0 )then
          do M = 1,TSPS(NS).NREC
            read(1,*,IOSTAT = ISO) TSPS(NS).TIM(M), PSERTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
            TSPS(NS).TIM(M) = TSPS(NS).TIM(M) + TSPS(NS).TOFFSET
            TSPS(NS).VAL(M,1) = G*(PSERTMP + ADDADJ)*RMULADJ                   ! *** m2/s2
            TSPS(NS).VAL(M,2) = 0.0
          enddo
        elseif( ITYPE == 1 )then
          do M = 1,TSPS(NS).NREC
            read(1,*,IOSTAT = ISO) TSPS(NS).TIM(M), PSERTMP, PSERTMPS
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
            TSPS(NS).TIM(M) = TSPS(NS).TIM(M)+TSPS(NS).TOFFSET
            TSPS(NS).VAL(M,1) = G*(PSERTMP  + ADDADJ)*RMULADJ                  ! *** m2/s2
            TSPS(NS).VAL(M,2) = G*(PSERTMPS + ADDADJS)*RMULADJS                ! *** m2/s2
          enddo
        endif

        ! *** Compute the average pressure value for this series.  Used for radiation BC's
        do M = 2,TSPS(NS).NREC
          TDELTA = TSPS(NS).TIM(M) - TSPS(NS).TIM(M-1)
          PSERAVG(NS,1) = PSERAVG(NS,1) + TSPS(NS).VAL(M,1)*TDELTA
          PSERAVG(NS,2) = PSERAVG(NS,2) + TSPS(NS).VAL(M,2)*TDELTA
        enddo
        PSERAVG(NS,1) = PSERAVG(NS,1)/(TSPS(NS).TIM(TSPS(NS).NREC) - TSPS(NS).TIM(1))
        PSERAVG(NS,2) = PSERAVG(NS,2)/(TSPS(NS).TIM(TSPS(NS).NREC) - TSPS(NS).TIM(1))
      enddo
      close(1)
    endif

    call Broadcast_Array(PSERZDS, master_id)
    call Broadcast_Array(PSERZDF, master_id)
    call Broadcast_Array(INTPSER, master_id)
    call Broadcast_Array(PSERAVG, master_id)

    do NS = 1,NPSER
      call Broadcast_Scalar(TSPS(NS).NREC ,    master_id)
      call Broadcast_Scalar(TSPS(NS).TMULT ,   master_id)
      call Broadcast_Scalar(TSPS(NS).TOFFSET , master_id)
      call Broadcast_Array(TSPS(NS).TIM,       master_id)
      call Broadcast_Array(TSPS(NS).VAL,       master_id)
    enddo

  endif
6776 FORMAT(A20)

  ! *** Read in volumetric source or river inflow time series from the file QSER.INP
  if( NQSER >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING QSER.INP'
      open(1,FILE = 'qser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      do NS = 1,NQSER
        read(1,*,IOSTAT = ISO) ISTYP, TSFL(NS).NREC, TSFL(NS).TMULT, TSFL(NS).TOFFSET, RMULADJ, ADDADJ, ICHGQS
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')

        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')
          do M = 1,TSFL(NS).NREC
            read(1,*,IOSTAT = ISO)TSFL(NS).TIM(M),QSERTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')
            TSFL(NS).TIM(M) = TSFL(NS).TIM(M)+TSFL(NS).TOFFSET
            QSERTMP = RMULADJ*(QSERTMP+ADDADJ)
            if( ICHGQS == 1)  QSERTMP = MAX(QSERTMP,0.0)
            if( ICHGQS == -1) QSERTMP = MIN(QSERTMP,0.0)
            do K = 1,KC
              TSFL(NS).VAL(M,K) = QSERTMP*WKQ(K)
            enddo
          enddo
        else
          do M = 1,TSFL(NS).NREC
            read(1,*,IOSTAT = ISO)TSFL(NS).TIM(M),(TSFL(NS).VAL(M,K), K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')
            TSFL(NS).TIM(M) = TSFL(NS).TIM(M)+TSFL(NS).TOFFSET
            do K = 1,KC
              TSFL(NS).VAL(M,K) = RMULADJ*(TSFL(NS).VAL(M,K) + ADDADJ)
              if( ICHGQS == 1)  TSFL(NS).VAL(M,K) = MAX(TSFL(NS).VAL(M,K),0.0)
              if( ICHGQS == -1) TSFL(NS).VAL(M,K) = MIN(TSFL(NS).VAL(M,K),0.0)
            enddo
          enddo
        endif
      enddo
      close(1)
    endif

    do NS = 1,NQSER
      call Broadcast_Scalar(TSFL(NS).NREC ,    master_id)
      call Broadcast_Scalar(TSFL(NS).TMULT ,   master_id)
      call Broadcast_Scalar(TSFL(NS).TOFFSET , master_id)
      call Broadcast_Array(TSFL(NS).TIM,       master_id)
      call Broadcast_Array(TSFL(NS).VAL,       master_id)
    enddo

  endif
2222 FORMAT(2I5,F12.7,F12.4)

  ! *** Read in flow withdrawl-return flow and concentration rise time series from the file QWRS.INP
  if( NQWRSR >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING QWRS.INP'
      open(1,FILE = 'qwrs.inp',STATUS = 'UNKNOWN')

      NCTMP = 3 + NDYM + NSED + NSND + NTOX
      if( ISTRAN(8) > 0 ) NCTMP = NCTMP + NWQV    ! *** Includes inactive WQ parameters also

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      do NS = 1,NQWRSR
        read(1,*,IOSTAT = ISO) ISTYP, TSWR(NS).NREC, TSWR(NS).TMULT, TSWR(NS).TOFFSET, RMULADJ, ADDADJ
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QWRS.INP')
        if( ISTYP == 0 )then
          ! *** FLOW ONLY.  NO RISE/FALL
          do NC = 1,NCTMP
            do M = 1,TSWR(NS).NREC
              TSWR(NS).VAL(M,NC) = 0.
            enddo
          enddo
          do M = 1,TSWR(NS).NREC
            read(1,*,IOSTAT = ISO) TSWR(NS).TIM(M), TSWR(NS).VAL(M,0)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QWRS.INP')
            TSWR(NS).TIM(M)   = TSWR(NS).TIM(M) + TSWR(NS).TOFFSET
            TSWR(NS).VAL(M,0) = (RMULADJ*(TSWR(NS).VAL(M,0)+ADDADJ))
          enddo
        else
          ! *** FLOW WITH RISE/FALL
          do M = 1,TSWR(NS).NREC
            read(1,*,IOSTAT = ISO) TSWR(NS).TIM(M), TSWR(NS).VAL(M,0),(TSWR(NS).VAL(M,NC),NC = 1,NCTMP)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QWRS.INP')
            TSWR(NS).TIM(M)   = TSWR(NS).TIM(M) + TSWR(NS).TOFFSET
            TSWR(NS).VAL(M,0) = (RMULADJ*(TSWR(NS).VAL(M,0) + ADDADJ))
          enddo
        endif
      enddo
      close(1)

    endif

    do NS = 1,NQWRSR
      call Broadcast_Scalar(TSWR(NS).NREC ,    master_id)
      call Broadcast_Scalar(TSWR(NS).TMULT ,   master_id)
      call Broadcast_Scalar(TSWR(NS).TOFFSET , master_id)
      call Broadcast_Array(TSWR(NS).TIM,       master_id)
      call Broadcast_Array(TSWR(NS).VAL,       master_id)
    enddo

  endif

  ! *** Read in groundwater inflow/outflow and concentration time from the file GWSER.INP
  if( ISGWIT == 2 )then
    ! *** SKIP OVER TITLE AND AND HEADER LINES
    NCTMP = 3 + NDYM + NTOX + NSED + NSND

    ! *** READ IN GW CONCENTRATIONS
    if( process_id == master_id )then
      write(*,'(A)')'READING GWSER.INP'
      open(1,FILE = 'gwser.inp',STATUS = 'UNKNOWN')    ! TODO - TSR

      STR = READSTR(1)
      read(1,*)NGWSER
      if( NGWSER > 0 )then
        do NS = 1,NGWSER
          read(1,*,IOSTAT = ISO) MGWSER(NS), TCGWSER(NS), TAGWSER(NS), RMULADJ, ADDADJ, IGWSER(NS)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE GWSER.INP')
          do M = 1,MGWSER(NS)
            read(1,*,IOSTAT = ISO)TGWSER(M,NS),GWSER(M,NS),(GWCSER(M,NS,NC),NC = 1,NCTMP)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE GWSER.INP')
            TGWSER(M,NS) = TGWSER(M,NS)+TAGWSER(NS)
            GWSER(M,NS)  = RMULADJ*(GWSER(M,NS) + ADDADJ)
          enddo
        enddo
      endif
      close(1)
    endif !***End on master

    call Broadcast_Array(MGWSER, master_id)
    call Broadcast_Array(TCGWSER, master_id)
    call Broadcast_Array(TAGWSER, master_id)
    call Broadcast_Array(IGWSER, master_id)
    call Broadcast_Array(GWSER,  master_id)
    call Broadcast_Array(GWCSER,  master_id)
    call Broadcast_Array(TGWSER, master_id)

  endif

  ! *** Read in open boundary or volumetric source salinity time serie from the file SSER.INP
  if( ISTRAN(1) >= 1 .and. NCSER(1) >= 1 )then
    NC = 1
    if( process_id == master_id )then
      write(*,'(A)')'READING SSER.INP'
      open(1,FILE = 'sser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      NC = 1
      do NS = 1,NCSER(NC)
        read(1,*,IOSTAT = ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TOFFSET, RMULADJ, ADDADJ
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
          do M = 1,TSSAL(NS).NREC
            read(1,*,IOSTAT = ISO) TSSAL(NS).TIM(M),CSERTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
            TSSAL(NS).TIM(M) = TSSAL(NS).TIM(M) + TOFFSET
            do K = 1,KC
              TSSAL(NS).VAL(M,K)  = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
            enddo
          enddo
        else
          do M = 1,TSSAL(NS).NREC
            read(1,*,IOSTAT = ISO) TSSAL(NS).TIM(M),(TSSAL(NS).VAL(M,K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
            TSSAL(NS).TIM(M) = TSSAL(NS).TIM(M) + TOFFSET
            do K = 1,KC
              TSSAL(NS).VAL(M,K) = RMULADJ*(TSSAL(NS).VAL(M,K) + ADDADJ)
            enddo
          enddo
        endif
      enddo
      close(1)
    endif


    call Broadcast_Array(MCSER,  master_id)
    call Broadcast_Array(TCCSER, master_id)

    do NS = 1,NCSER(NC)
      call Broadcast_Array(TSSAL(NS).VAL, master_id)
      call Broadcast_Array(TSSAL(NS).TIM, master_id)
    enddo

  endif

  ! *** Read in open boundary or volumetric source temperature time series from the file TSER.INP
  if( ( ISTRAN(2) >= 1 .and. NCSER(2) >= 1 ) .or. IVOLTEMP > 1 )then

    NC = 2
    if( process_id == master_id )then
      write(*,'(A)')'READING TSER.INP'
      open(1,FILE = 'tser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      do NS = 1,NCSER(NC)
        read(1,*,IOSTAT = ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TOFFSET, RMULADJ, ADDADJ
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSTEM(NS).TIM(M),CSERTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
            TSTEM(NS).TIM(M) = TSTEM(NS).TIM(M) + TOFFSET
            do K = 1,KC
              TSTEM(NS).VAL(M,K) = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
            enddo
          enddo
        else
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSTEM(NS).TIM(M),(TSTEM(NS).VAL(M,K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
            TSTEM(NS).TIM(M) = TSTEM(NS).TIM(M) + TOFFSET
            do K = 1,KC
              TSTEM(NS).VAL(M,K) = RMULADJ*(TSTEM(NS).VAL(M,K) + ADDADJ)
            enddo
          enddo
        endif
      enddo
      close(1)
    endif

    call Broadcast_Array(MCSER,  master_id)
    call Broadcast_Array(TCCSER, master_id)
    do NS = 1,NCSER(NC)
      call Broadcast_Array(TSTEM(NS).TIM, master_id)
      call Broadcast_Array(TSTEM(NS).VAL, master_id)
    enddo

  endif

  ! *** Read in open boundary or volumetric source dye concentration time series from the file DSER.INP
  if( ISTRAN(3) >= 1 .and. NCSER(3) >= 1 )then
    NC = 3
    if( process_id == master_id )then
      write(*,'(A)')'READING DSER.INP'
      open(1,FILE = 'dser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      do NS = 1,NCSER(NC)
        read(1,*,IOSTAT = ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TOFFSET, RMULADJ, ADDADJ
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSDYE(NS,1).TIM(M),CSERTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
            TSDYE(NS,1).TIM(M) = TSDYE(NS,1).TIM(M) + TOFFSET
            do K = 1,KC
              TSDYE(NS,1).VAL(M,K) = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
            enddo
            do MD = 2,NDYE
              TSDYE(NS,MD).TIM(M) = TSDYE(NS,1).TIM(M)
              read(1,*,IOSTAT = ISO) CSERTMP
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
              do K = 1,KC
                TSDYE(NS,MD).VAL(M,K) = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
              enddo
            enddo
          enddo
        else
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSDYE(NS,1).TIM(M),(TSDYE(NS,1).VAL(M,K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
            TSDYE(NS,1).TIM(M) = TSDYE(NS,1).TIM(M) + TOFFSET
            do K = 1,KC
              TSDYE(NS,1).VAL(M,K) = RMULADJ*(TSDYE(NS,1).VAL(M,K) + ADDADJ)
            enddo
            do MD = 2,NDYE
              TSDYE(NS,MD).TIM(M) = TSDYE(NS,1).TIM(M)
              read(1,*,IOSTAT = ISO) (TSDYE(NS,MD).VAL(M,K), K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
              do K = 1,KC
                TSDYE(NS,MD).VAL(M,K) = RMULADJ*(TSDYE(NS,MD).VAL(M,K) + ADDADJ)
              enddo
            enddo
          enddo
        endif
      enddo
      close(1)
    endif !*** end master

    do NS = 1,NCSER(NC)
      call Broadcast_Scalar(MCSER(NS,NC),  master_id)
      call Broadcast_Scalar(TCCSER(NS,NC), master_id)
      do MD = 1,NDYE
        call Broadcast_Array(TSDYE(NS,MD).TIM, master_id)
        call Broadcast_Array(TSDYE(NS,MD).VAL, master_id)
      enddo
    enddo

  endif

  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SHELL FISH LARVAE
  ! *** TIME SERIES FROM THE FILE SFSER.INP
  if(ISTRAN(4) >= 1 .and. NCSER(4) >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING SFSER.INP'
      open(1,FILE = 'sfser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      NC = 4
      do NS = 1,NCSER(NC)
        read(1,*,IOSTAT = ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TOFFSET,RMULADJ,ADDADJ
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSSFL(NS).TIM(M),CSERTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
            TSSFL(NS).TIM(M) = TSSFL(NS).TIM(M) + TOFFSET
            do K = 1,KC
              TSSFL(NS).VAL(M,K) = (RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            enddo
          enddo
        else
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSSFL(NS).TIM(M),(TSSFL(NS).VAL(M,K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
            TSSFL(NS).TIM(M) = TSSFL(NS).TIM(M) + TOFFSET
            do K = 1,KC
              TSSFL(NS).VAL(M,K) = RMULADJ*(TSSFL(NS).VAL(M,K)+ADDADJ)
            enddo
          enddo
        endif
      enddo
      close(1)
    endif !*** end master

    do NS = 1,NCSER(NC)
      call Broadcast_Scalar(MCSER(NS,NC),  master_id)
      call Broadcast_Scalar(TCCSER(NS,NC), master_id)
      call Broadcast_Array(TSSFL(NS).TIM,  master_id)
      call Broadcast_Array(TSSFL(NS).VAL,  master_id)
    enddo

  endif

  ! *** Read in open boundary or volumetric source cohesive sediment concentration time series from the file SDSER.INP
  if( ISTRAN(6) >= 1 .and. NSED > 0 )then
    NC = 6
    if( NCSER(NC) >= 1 )then

      if( process_id == master_id )then
        write(*,'(A)')'READING SDSER.INP'
        open(1,FILE = 'sdser.inp',STATUS = 'UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        do NS = 1,NCSER(NC)
          read(1,*,IOSTAT = ISO)ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TOFFSET, RMULADS(1), ADDADS(1)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
          if( NSED>1 )then
            do NT = 2,NSED
              read(1,*,IOSTAT = ISO) RMULADS(NT),ADDADS(NT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
            enddo
          endif
          if( ISTYP == 1 )then
            read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
            do M = 1,MCSER(NS,NC)
              read(1,*,IOSTAT = ISO) TSSED(NS,1).TIM(M),CSERTMP
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
              TSSED(NS,1).TIM(M) = TSSED(NS,1).TIM(M) + TOFFSET
              do K = 1,KC
                TSSED(NS,1).VAL(M,K) = ( RMULADS(1)*( CSERTMP + ADDADS(1) ) )*WKQ(K)
              enddo
              do NT = 2,NSED
                !NTT = NT-1
                TSSED(NS,NT).TIM(M) = TSSED(NS,1).TIM(M)
                read(1,*,IOSTAT = ISO) CSERTMP
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
                do K = 1,KC
                  TSSED(NS,NT).VAL(M,K)  = ( RMULADS(NT)*( CSERTMP + ADDADS(NT) ) )*WKQ(K)
                enddo
              enddo
            enddo
          else
            do M = 1,MCSER(NS,NC)
              read(1,*,IOSTAT = ISO) TSSED(NS,1).TIM(M),(TSSED(NS,1).VAL(M,K),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
              TSSED(NS,1).TIM(M) = TSSED(NS,1).TIM(M) + TOFFSET
              do K = 1,KC
                TSSED(NS,1).VAL(M,K) = RMULADS(1)*( TSSED(NS,1).VAL(M,K) + ADDADS(1) )
              enddo
              do NT = 2,NSED
                !NTT = NT-1
                TSSED(NS,NT).TIM(M) = TSSED(NS,1).TIM(M)
                read(1,*,IOSTAT = ISO)(TSSED(NS,NT).VAL(M,K), K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
                do K = 1,KC
                  TSSED(NS,NT).VAL(M,K) = RMULADS(NT)*( TSSED(NS,NT).VAL(M,K) + ADDADS(NT) )
                enddo
              enddo
            enddo
          endif
        enddo
        close(1)
      endif

      do NS = 1,NCSER(NC)
        call Broadcast_Scalar(MCSER(NS,NC),       master_id)
        call Broadcast_Scalar(TCCSER(NS,NC),      master_id)
        do NT = 1,NSED
          call Broadcast_Array(TSSED(NS,NT).TIM, master_id)
          call Broadcast_Array(TSSED(NS,NT).VAL, master_id)
        enddo
      enddo
    endif
  endif

  ! *** Read in open boundary or volumetric source noncohesive sediment concentration time series from the file SNSER.INP
  if( ISTRAN(7) >= 1 .and. NSND > 0 )then
    NC = 7
    if( NCSER(NC) >= 1 )then
      if( process_id == master_id )then
        write(*,'(A)')'READING SNSER.INP'
        open(1,FILE = 'snser.inp',STATUS = 'UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)

        do NS = 1,NCSER(NC)
          read(1,*,IOSTAT = ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TOFFSET,RMULADS(1),ADDADS(1)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
          if( NSND > 1 )then
            do NT = 2,NSND
              read(1,*,IOSTAT = ISO)RMULADS(NT),ADDADS(NT)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
            enddo
          endif
          if( ISTYP == 1 )then
            read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
            do M = 1,MCSER(NS,NC)
              read(1,*,IOSTAT = ISO) TSSND(NS,1).TIM(M),CSERTMP
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
              TSSND(NS,1).TIM(M) = TSSND(NS,1).TIM(M) + TOFFSET
              do K = 1,KC
                TSSND(NS,1).VAL(M,K) = (RMULADS(1)*(CSERTMP+ADDADS(1)))*WKQ(K)
              enddo
              do NT = 2,NSND
                !NTT = NT-1
                TSSND(NS,NT).TIM(M) = TSSND(NS,1).TIM(M)
                read(1,*,IOSTAT = ISO)CSERTMP
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
                do K = 1,KC
                  TSSND(NS,NT).VAL(M,K)  = (RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
                enddo
              enddo
            enddo
          else
            do M = 1,MCSER(NS,NC)
              read(1,*,IOSTAT = ISO) TSSND(NS,1).TIM(M),(TSSND(NS,1).VAL(M,K),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
              TSSND(NS,1).TIM(M) = TSSND(NS,1).TIM(M) + TOFFSET
              do K = 1,KC
                TSSND(NS,1).VAL(M,K) = RMULADS(1)*(TSSND(NS,1).VAL(M,K) + ADDADS(1))
              enddo
              do NT = 2,NSND
                !NTT = NT-1
                TSSND(NS,NT).TIM(M) = TSSND(NS,1).TIM(M)
                read(1,*,IOSTAT = ISO)(TSSND(NS,NT).VAL(M,K), K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
                do K = 1,KC
                  TSSND(NS,NT).VAL(M,K)  = RMULADS(NT)*(TSSND(NS,NT).VAL(M,K)+ADDADS(NT))
                enddo
              enddo
            enddo
          endif
        enddo
        close(1)
      endif

      do NS = 1,NCSER(NC)
        call Broadcast_Scalar(MCSER(NS,NC),       master_id)
        call Broadcast_Scalar(TCCSER(NS,NC),      master_id)
        do NT = 1,NSND
          call Broadcast_Array(TSSND(NS,NT).TIM, master_id)
          call Broadcast_Array(TSSND(NS,NT).VAL, master_id)
        enddo
      enddo
    endif
  endif

  ! *** Read in open boundary or volumetric source toxic contaminant concentration time series from the file TXSER.INP
  if( ISTRAN(5) >= 1 .and. NTOX > 0 )then
    NC = 5
    if( NCSER(NC) >= 1 )then
      if( process_id == master_id )then
        write(*,'(A)')'READING TXSER.INP'
        open(1,FILE = 'txser.inp',STATUS = 'UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        do NS = 1,NCSER(NC)
          read(1,*,IOSTAT = ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TOFFSET,RMULADJ,ADDADJ
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
          if( ISTYP == 1 )then
            read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
            do M = 1,MCSER(NS,NC)
              read(1,*,IOSTAT = ISO) TSTOX(NS,1).TIM(M),CSERTMP
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
              TSTOX(NS,1).TIM(M) = TSTOX(NS,1).TIM(M) + TOFFSET
              do K = 1,KC
                TSTOX(NS,1).VAL(M,K) = (RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
              enddo
              do NT = 2,NTOX
                !NTT = NT-1
                TSTOX(NS,NT).TIM(M) = TSTOX(NS,1).TIM(M)
                read(1,*,IOSTAT = ISO) CSERTMP
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
                do K = 1,KC
                  TSTOX(NS,NT).VAL(M,K) = ( RMULADJ*( CSERTMP + ADDADJ) )*WKQ(K)
                enddo
              enddo
            enddo
          else
            do M = 1,MCSER(NS,NC)
              read(1,*,IOSTAT = ISO) TSTOX(NS,1).TIM(M),(TSTOX(NS,1).VAL(M,K), K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
              TSTOX(NS,1).TIM(M) = TSTOX(NS,1).TIM(M) + TOFFSET
              do K = 1,KC
                TSTOX(NS,1).VAL(M,K) = RMULADJ*(TSTOX(NS,1).VAL(M,K) + ADDADJ)
              enddo
              do NT = 2,NTOX
                !NTT = NT-1
                TSTOX(NS,NT).TIM(M) = TSTOX(NS,1).TIM(M)
                read(1,*,IOSTAT = ISO)(TSTOX(NS,NT).VAL(M,K), K = 1,KC)
                if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
                do K = 1,KC
                  TSTOX(NS,NT).VAL(M,K) = RMULADJ*(TSTOX(NS,NT).VAL(M,K) + ADDADJ)
                enddo
              enddo
            enddo
          endif
        enddo
        close(1)
      endif

      do NS = 1,NCSER(NC)
        call Broadcast_Scalar(MCSER(NS,NC),       master_id)
        call Broadcast_Scalar(TCCSER(NS,NC),      master_id)
        do NT = 1,NTOX
          call Broadcast_Array(TSTOX(NS,NT).TIM, master_id)
          call Broadcast_Array(TSTOX(NS,NT).VAL, master_id)
        enddo
      enddo

    endif

    ! *** TIME SERIES DRY DEPOSITION BEING USED
    if( SUM(TOXDEP(:).ITXDRY) > 0 )then
      write(*,'(A)')'READING TXDRY.INP'
      open(1,FILE = 'txdry.inp',STATUS = 'UNKNOWN')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*,IOSTAT = ISO) ISTYP, TXDRYSER(1).NREC, T1, T2, RMULADJ, ADDADJ
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')

      ! *** IGNORE ISTYP, ALWAYS use ONE VALUE
      read(1,*)                                                     ! *** SKIP SPLITS
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXDRY.INP')

      do M = 1,TXDRYSER(1).NREC
        read(1,*,IOSTAT = ISO) TXDRYSER(1).TIM(M), CSERTMP
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXDRY.INP')

        TXDRYSER(1).TIM(M)   = (TXDRYSER(1).TIM(M) + T2)*T1         ! *** TIMES ARE IN SECONDS
        TXDRYSER(1).VAL(M,1) = ((CSERTMP + ADDADJ)*RMULADJ)/86400.  ! *** MG/M2/SEC
        TOXDEP(1).ITDRY = 1                                         ! *** SET INITIAL TIME INDEX
        do NT = 2,NTOX
          read(1,*,IOSTAT = ISO) CSERTMP
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXDRY.INP')
          TXDRYSER(1).VAL(M,NT) = (CSERTMP + ADDADJ)*RMULADJ
          TOXDEP(NT).ITDRY = 1                                     ! *** SET INITIAL TIME INDEX
        enddo
      enddo
      close(1)

    endif

    ! *** TIME SERIES WET DEPOSITION BEING USED
    if( SUM(TOXDEP(:).ITXWET) > 0 )then
      write(*,'(A)')'READING TXWET.INP'
      open(1,FILE = 'txwet.inp',STATUS = 'UNKNOWN')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*,IOSTAT = ISO) ISTYP, TXWETSER(1).NREC, T1, T2, RMULADJ, ADDADJ
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')

      ! *** IGNORE ISTYP, ALWAYS use ONE VALUE
      read(1,*)                                                     ! *** SKIP SPLITS
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXWET.INP')

      do M = 1,TXWETSER(1).NREC
        read(1,*,IOSTAT = ISO) TXWETSER(1).TIM(M), CSERTMP
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXWET.INP')

        TXWETSER(1).TIM(M)   = (TXWETSER(1).TIM(M) + T2)*T1         ! *** TIMES ARE IN SECONDS
        TXWETSER(1).VAL(M,1) = (CSERTMP + ADDADJ)*RMULADJ           ! *** MG/M3
        TOXDEP(1).ITWET = 1                                         ! *** SET INITIAL TIME INDEX

        do NT = 2,NTOX
          read(1,*,IOSTAT = ISO) CSERTMP
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXWET.INP')
          TXWETSER(1).VAL(M,NT) = (CSERTMP + ADDADJ)*RMULADJ
          TOXDEP(NT).ITWET = 1                                      ! *** SET INITIAL TIME INDEX
        enddo
      enddo
      close(1)
    endif
  endif

  ! *** Read in free surface elevation or pressure controlled flow specification from the file QCTL.INP
  if( NQCTLT >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING QCTL.INP'
      open(1,FILE = 'qctl.inp',STATUS = 'UNKNOWN')
      if( DEBUG )then
        open(99,FILE = OUTDIR//'qctlck.inp',STATUS = 'UNKNOWN')
        close(99,STATUS = 'DELETE')
        open(99,FILE = OUTDIR//'qctlck.inp',STATUS = 'UNKNOWN')
      endif

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      do NS = 1,NQCTLT
        read(1,*, IOSTAT = ISO)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
        if( DEBUG )then
          write(99,991)NS
          write(99,992)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
        endif
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
        if( ISTYP == 0 )then
          do M = 1,MQCTL(NS)
            read(1,*,IOSTAT = ISO) HDIFCTL(M,NS),(QCTL(M,1,K,NS),K = 1,KC)
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
            do K = 1,KC
              QCTL(M,1,K,NS) = RMULADJ*(QCTL(M,1,K,NS)+ADDADJ)
            enddo
          enddo
        endif
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
          do M = 1,MQCTL(NS)
            read(1,*,IOSTAT = ISO) HDIFCTL(M,NS),QCTLTMP
            if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
            do K = 1,KC
              QCTL(M,1,K,NS) = RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
            enddo
          enddo
        endif
        if( ISTYP == 2 )then
          do MD = 1,MQCTL(NS)
            do MU = 1,MQCTL(NS)
              read(1,*,IOSTAT = ISO) HDIFCTL(MU,NS),HDIFCTD(MD,NS),(QCTL(MU,MD,K,NS),K = 1,KC)
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
              do K = 1,KC
                QCTL(MU,MD,K,NS) = RMULADJ*(QCTL(MU,MD,K,NS)+ADDADJ)
              enddo
            enddo
          enddo
        endif
        if( ISTYP == 3 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
          do MD = 1,MQCTL(NS)
            do MU = 1,MQCTL(NS)
              read(1,*,IOSTAT = ISO)HDIFCTL(MU,NS),HDIFCTD(MD,NS),QCTLTMP
              if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
              do K = 1,KC
                QCTL(MU,MD,K,NS) = RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
              enddo
            enddo
          enddo
        endif
        if( DEBUG )then
          if( ISTYP <= 1 )then
            do M = 1,MQCTL(NS)
              write(99,993)M,HDIFCTL(M,NS),(QCTL(M,1,K,NS),K = 1,KC)
            enddo
          endif
          if( ISTYP >= 2 )then
            do MD = 1,MQCTL(NS)
              do MU = 1,MQCTL(NS)
                write(99,994)MU,MD,HDIFCTL(MU,NS),HDIFCTD(MD,NS),(QCTL(MU,MD,K,NS),K = 1,KC)
              enddo
            enddo
          endif
        endif
      enddo
      close(1)
      if( DEBUG)CLOSE(99)
    endif

    call Broadcast_Array(MQCTL   , master_id)
    call Broadcast_Array(HCTLUA  , master_id)
    call Broadcast_Array(HCTLUM  , master_id)
    call Broadcast_Array(HCTLDA  , master_id)
    call Broadcast_Array(HCTLDM  , master_id)
    call Broadcast_Array(AQCTL   , master_id)
    call Broadcast_Array(HDIFCTL , master_id)
    call Broadcast_Array(HDIFCTD , master_id)
    call Broadcast_Array(QCTL    , master_id)

  endif

  ! *** Read control time-series from the file QCTLSER.INP
  if( NQCTLSER > 0 )then
    allocate(QCTLSER(NQCTLSER))

    if( process_id == master_id )then
      write(*,'(A)')'READING QCTLSER.INP'
      open(1,FILE = 'qctlser.inp',STATUS = 'UNKNOWN')
      if( DEBUG )then
        open(99,FILE = OUTDIR//'qctlser.log',STATUS = 'REPLACE')
      endif

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
    endif

    do NS = 1,NQCTLSER
      if( process_id == master_id )then
        read(1,*, IOSTAT = ISO) QCTLSER(NS).ITYPE, NX, QCTLSER(NS).TMUL, QCTLSER(NS).TADD, (ISQCTRL(M),M = 1,6)
        QCTLSER(NS).PARAM = 0
        QCTLSER(NS).COUNT = NX
        IVAL = 1
        do M = 1,6
          if(ISQCTRL(M) > 0) QCTLSER(NS).PARAM = QCTLSER(NS).PARAM + IVAL
          IVAL = 2*IVAL
        enddo
        if( DEBUG )then
          write(99,991) NS
          write(99,992) QCTLSER(NS).ITYPE, NX, QCTLSER(NS).TMUL, QCTLSER(NS).TADD, QCTLSER(NS).PARAM
        endif
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTLSER.INP')
      endif
      call Broadcast_Scalar(QCTLSER(NS).ITYPE   , master_id)
      call Broadcast_Scalar(QCTLSER(NS).TMUL    , master_id)
      call Broadcast_Scalar(QCTLSER(NS).TADD    , master_id)
      call Broadcast_Scalar(QCTLSER(NS).PARAM   , master_id)
      call Broadcast_Scalar(QCTLSER(NS).COUNT   , master_id)

      NX = QCTLSER(NS).COUNT
      allocate(QCTLSER(NS).TIME(NX))
      allocate(QCTLSER(NS).ID(NX))
      allocate(QCTLSER(NS).HEIGHT(NX))
      allocate(QCTLSER(NS).WIDTH(NX))
      allocate(QCTLSER(NS).SILL(NX))
      !ALLOCATE(QCTLSER(NS).NUM(NX))
      allocate(QCTLSER(NS).FLOW(NX))
      !ALLOCATE(QCTLSER(NS).RATE(NX))

      if( process_id == master_id )then
        do M = 1,NX
          read(1,*,IOSTAT = ISO) TIMEQCTLSERIES, QCTLSER(NS).HEIGHT(M), QCTLSER(NS).WIDTH(M), QCTLSER(NS).SILL(M), &
            NDUM, QCTLSER(NS).FLOW(M), QCTLSER(NS).ID(M)    !,QCTLSER(NS).RATE(M)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTLSER.INP')
          QCTLSER(NS).TIME(M) = TIMEQCTLSERIES*QCTLSER(NS).TMUL + QCTLSER(NS).TADD
        enddo
      endif
      call Broadcast_Array(QCTLSER(NS).TIME     , master_id)
      call Broadcast_Array(QCTLSER(NS).ID       , master_id)
      call Broadcast_Array(QCTLSER(NS).HEIGHT   , master_id)
      call Broadcast_Array(QCTLSER(NS).WIDTH    , master_id)
      call Broadcast_Array(QCTLSER(NS).SILL     , master_id)
      !Call Broadcast_Array(QCTLSER(NS).NUM      , master_id)  ! NOT USED
      call Broadcast_Array(QCTLSER(NS).FLOW     , master_id)
    enddo
    if( process_id == master_id )then
      close(1)
      if( DEBUG ) close(99)
    endif


  endif

  ! *** Read control rules from the file QCTRULES.INP
  if( NQCTRULES > 0 )then
    allocate(QCTRULES(NQCTRULES))

    if( process_id == master_id )then
      write(*,'(A)')'READING QCTRULES.INP'
      open(1,FILE = 'qctrules.inp',STATUS = 'UNKNOWN')

      open(99,FILE = OUTDIR//'QCTRULES.log',STATUS = 'UNKNOWN')
      close(99,STATUS = 'DELETE')
      open(99,FILE = OUTDIR//'QCTRULES.log',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      do NS = 1,NQCTRULES
        STR = READSTR(1)       ! *** Skip all comments to the type of control flags

        read(1,*, IOSTAT = ISO) NX,NX,(ISQCTRL(M),M = 1,6)  !,QCTRULES(NS).PARAM
        QCTRULES(NS).PARAM = 0
        IVAL = 1
        do M = 1,6
          if(ISQCTRL(M) > 0) QCTRULES(NS).PARAM = QCTRULES(NS).PARAM + IVAL
          IVAL = 2*IVAL
        enddo
        write(99,*) NS,NX,QCTRULES(NS).PARAM
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTRULES.INP')

        QCTRULES(NS).NTRGON = 0
        QCTRULES(NS).NTROFF = 0

        STR = READSTR(1)       ! *** Skip intermediate comments

        allocate(RULES(NX,9))
        allocate(IDX(NX))
        do M = 1,NX
          read(1,*,IOSTAT = ISO) RULES(M,1), ISTYP, RULES(M,4), RULES(M,5), RULES(M,6), NDUM, RULES(M,8), RULES(M,9), IVAL
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTRULES.INP')
          IDX(M) = M
          RULES(M,2) = ISTYP
          RULES(M,3) = IVAL
          RULES(M,7) = NDUM
          if( ISTYP == 0 )then
            QCTRULES(NS).NTROFF = QCTRULES(NS).NTROFF + 1
          else
            QCTRULES(NS).NTRGON = QCTRULES(NS).NTRGON + 1
          endif
        enddo
        ! SORT TRIGGER LEVELS DESCENDING
        do M = 1,NX-1
          do L  = M+1,NX
            if(RULES(IDX(L),1) > RULES(IDX(M),1) )then
              IVAL = IDX(L)
              IDX(L) = IDX(M)
              IDX(M) = IVAL
            endif
          enddo
        enddo

        ! *** Setup list of triggers to open the gates/pumps
        allocate(QCTRULES(NS).TRGON(QCTRULES(NS).NTRGON))
        L = 0
        do M = 1,NX
          if( RULES(M,2) == 1. )then
            L = L + 1
            QCTRULES(NS).TRGON(L).LEVEL = RULES(M,1)             ! TRIGGER LEVEL
            QCTRULES(NS).TRGON(L).STATE = INT(RULES(M,2))        ! 0 = OFF, 1 = ON
            QCTRULES(NS).TRGON(L).ID = INT(RULES(M,3))           ! INDEX OF RATING CURVE/MASK FOR CONTROL parameterS
            QCTRULES(NS).TRGON(L).HEIGHT = RULES(M,4)            ! OPENING HEIGHT (M), FOR UPWARD OPENING
            QCTRULES(NS).TRGON(L).WIDTH = RULES(M,5)             ! OPENING WIDTH (M), FOR SIDEWARD OPENING
            QCTRULES(NS).TRGON(L).SILL = RULES(M,6)              ! SILL LEVEL CHANGE (M), FOR DOWNWARD OPENING
            !QCTRULES(NS).TRGON(L).UNITS = INT(RULES(M,7))       ! NUMBER OF GATES, PUMP UNITS   not used
            QCTRULES(NS).TRGON(L).FLOW = RULES(M,8)              ! FLOW RATE (M3/S), FOR PUMPS
            QCTRULES(NS).TRGON(L).RATE = RULES(M,9)              ! RATE OF OPENING/CLOSING

            write(99,*) L,QCTRULES(NS).TRGON(L).LEVEL, QCTRULES(NS).TRGON(L).STATE, &
              QCTRULES(NS).TRGON(L).ID,    QCTRULES(NS).TRGON(L).HEIGHT, QCTRULES(NS).TRGON(L).WIDTH, &
              QCTRULES(NS).TRGON(L).SILL,  QCTRULES(NS).TRGON(L).FLOW,   QCTRULES(NS).TRGON(L).RATE

          endif
        enddo

        ! *** Setup list of triggers to close the gates/pumps
        allocate(QCTRULES(NS).TROFF(QCTRULES(NS).NTROFF))
        L = 0
        do M = NX,1,-1
          if( RULES(M,2) == 0. )then
            L = L + 1
            QCTRULES(NS).TROFF(L).LEVEL = RULES(M,1)             ! TRIGGER LEVEL
            QCTRULES(NS).TROFF(L).STATE = INT(RULES(M,2))        ! 0 = OFF, 1 = ON
            QCTRULES(NS).TROFF(L).ID = INT(RULES(M,3))           ! INDEX OF RATING CURVE/MASK FOR CONTROL parameterS
            QCTRULES(NS).TROFF(L).HEIGHT = RULES(M,4)            ! OPENING HEIGHT (M), FOR UPWARD OPENING
            QCTRULES(NS).TROFF(L).WIDTH = RULES(M,5)             ! OPENING WIDTH (M), FOR SIDEWARD OPENING
            QCTRULES(NS).TROFF(L).SILL = RULES(M,6)              ! SILL LEVEL CHANGE (M), FOR DOWNWARD OPENING
            !QCTRULES(NS).TROFF(L).UNITS = INT(RULES(M,7))       ! NUMBER OF GATES, PUMP UNITS   not used
            QCTRULES(NS).TROFF(L).FLOW = RULES(M,8)              ! FLOW RATE (M3/S), FOR PUMPS
            QCTRULES(NS).TROFF(L).RATE = RULES(M,9)              ! RATE OF OPENING/CLOSING
            write(99,*) L,QCTRULES(NS).TROFF(L).LEVEL, QCTRULES(NS).TROFF(L).STATE, &
              QCTRULES(NS).TROFF(L).ID,    QCTRULES(NS).TROFF(L).HEIGHT, QCTRULES(NS).TROFF(L).WIDTH, &
              QCTRULES(NS).TROFF(L).SILL,  QCTRULES(NS).TROFF(L).FLOW,   QCTRULES(NS).TROFF(L).RATE
          endif
        enddo
        deallocate(RULES,IDX)
      enddo
      close(1)
      close(99)
    endif ! *** end on master


    do NS = 1,NQCTRULES
      call Broadcast_Scalar(QCTRULES(NS).PARAM  , master_id)

      call Broadcast_Scalar(QCTRULES(NS).NTRGON  , master_id)
      ! *** only allocate on other processes besides master
      if( process_id /= master_id )then
        allocate(QCTRULES(NS).TRGON(QCTRULES(NS).NTRGON))
      endif

      do M = 1,QCTRULES(NS).NTRGON
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).LEVEL  , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).STATE  , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).ID     , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).HEIGHT , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).WIDTH  , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).SILL   , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).FLOW   , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TRGON(M).RATE   , master_id)
      enddo

      call Broadcast_Scalar(QCTRULES(NS).NTROFF   , master_id)
      ! *** only allocate on other processes besides master
      if(process_id /= master_id )then
        allocate(QCTRULES(NS).TROFF(QCTRULES(NS).NTROFF))
      endif

      do M = 1,QCTRULES(NS).NTROFF
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).LEVEL  , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).STATE  , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).ID     , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).HEIGHT , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).WIDTH  , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).SILL   , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).FLOW   , master_id)
        call Broadcast_Scalar(QCTRULES(NS).TROFF(M).RATE   , master_id)
      enddo
    enddo
  endif

  ! *** Q/C CHECK CONTROL DATA
  if( NQCTL > 0 .and. (NQCTLSER > 0 .or. NQCTRULES > 0) )then
    do L = 1,NQCTL
      IVAL = HSCTL_GL(L).ID
      if( HSCTL_GL(L).ITYPE == 1 )then
        ! *** STRUCTURE IS CONTROLLED BY TIME-SERIES
        if( IVAL < 1 .or. IVAL > NQCTLSER )then
          call STOPP('DATA ERROR: TIME-SERIES INDEX FOR HYDRAULIC STRUCTURE CONTROL IS INVALID!')
        endif

      elseif( HSCTL_GL(L).ITYPE == 2 .or. HSCTL_GL(L).ITYPE == 3 )then
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
        if( IVAL < 1 .or. IVAL > NQCTRULES )then
          call STOPP('DATA ERROR: RULES INDEX FOR HYDRAULIC STRUCTURE CONTROL IS INVALID!')
        endif
        if( (QCTRULES(IVAL).PARAM .and. 32) == 32 )then
          ! *** SET OF LOOKUP TABLE
          do I = 1,QCTRULES(IVAL).NTRGON
            ID = QCTRULES(IVAL).TRGON(I).ID
            if( ID < 1 .or. ID > NQCTLT) CALL STOPP('DATA ERROR: INVALID INDEX FOR LOOKUP TABLE DEFINED IN CONTROL RULES')
          enddo
          do I = 1,QCTRULES(IVAL).NTROFF
            ID = QCTRULES(IVAL).TROFF(I).ID
            if( ID < 1 .or. ID > NQCTLT) CALL STOPP('DATA ERROR: INVALID INDEX FOR LOOKUP TABLE DEFINED IN CONTROL RULES')
          enddo
        elseif( QCTRULES(IVAL).PARAM > 0 )then
          ! *** GATE OPENING
        else
          ! *** PUMP FLOWS
        endif
        if( HSCTL_GL(L).ITYPE == 2 )then
          ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
          !LU = LIJ(HSCTL_GL(NCTL).IREFUP, HSCTL_GL(NCTL).JREFUP)
        elseif( HSCTL_GL(L).ITYPE == 3 )then
          ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
          !LU = LIJ(HSCTL_GL(NCTL).IREFUP, HSCTL_GL(NCTL).JREFUP)
          !LD = LIJ(HSCTL_GL(NCTL).IREFDN, HSCTL_GL(NCTL).JREFDN)
        endif
      else
        ! *** STRUCTURE IS UNCONTROLLED
      endif
    enddo
  endif

  if( NQWR > 0 .and. NQCTRULES > 0 )then
    do L = 1,NQWR
      IVAL = WITH_RET_CTL_GL(L).ID
      if( WITH_RET_CTL_GL(L).ITYPE == 1 .or. WITH_RET_CTL_GL(L).ITYPE == 2 )then
        ! *** W/R IS CONTROLLED BY OPERATION RULES
        if( IVAL < 1 .or. IVAL > NQCTRULES )then
          call STOPP('DATA ERROR: RULES INDEX FOR WITHDRAWAL/RETURN CONTROL IS INVALID!')
        endif
      endif
    enddo
  endif

991 FORMAT(/,'CONTROL TABLE NS  = ',I5,/)
992 FORMAT(2I5,7F10.4)
993 FORMAT(I5,11F10.4)
994 FORMAT(2I5,11F10.4)
1001 FORMAT(/,'READ ERROR FROM FILE EFDC.INP ON CARD ',A3/)
1002 FORMAT(/,'INPUT ECHO NCARD = ',A/)

  ! *** Initialization of atmospheric data
  do L = 1,LC
    PATMT(L) = 1000.
    TATMT(L) = 0.
    RAINT(L) = 0.
    EVAPT(L) = 0.
    SOLSWRT(L) = 1.  ! *** Address SUNDAY.INP Option
    CLOUDT(L) = 0.
    RHAT(L) = 0.
    VPAT(L) = 0.
    CLEVAP(L) = 0.
    CCNHTT(L) = 0.
  enddo

  ! *** Read in atmospheric data time series from the file ASER.INP
  LDAYLIGHT = .FALSE.
  LRAIN = .FALSE.
  if( NASER > 0 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING ASER.INP'
      open(1,FILE = 'aser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      call SKIPCOM(1,'*')
      read(1,'(A)')TEXT
      IASERVER = PARSE_REAL(TEXT)*1000.
      STR = READSTR(1)

      if( IASERVER < 7300 )then
        IEVAP = -1
      endif

      write(*,'(A,I5)')'  NUMBER OF ATMOSPHERIC SERIES = ',NASER
      write(*,'(A,L3)')'  COMPUTESOLRAD = ',COMPUTESOLRAD
      write(*,'(A,F10.2)')'  DS_LONG = ',DS_LONG
      write(*,'(A,F10.2)')'  DS_LAT = ',DS_LAT

      do NA = 1,NASER
        read(1,*,IOSTAT = ISO) M, TSATM(NA).TMULT, TSATM(NA).TOFFSET, IRELH(NA), RAINCVT, EVAPCVT, SOLRCVT, CLDCVT
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE ASER.INP')
        if( NA == 1 .and. IEVAP == -1 )then
          IEVAP = 1                     ! *** Use ASER DATA
          if( EVAPCVT < 0 ) IEVAP = 2   ! *** LEGACY INPUT TO COMPUTE EVAPORATION USING EFDC ORIGINAL APPROACH
        elseif( IEVAP == 0 )then
          ! *** IGNORE EVAPORATION BUT ALWAYS INCLUDE RAINFALL FROM ASER
          EVAPCVT = 0.
        endif
        write(*,'(A,I5,I10,F6.1)')'  SERIES, NPTS  = ',NA,TSATM(NA).NREC

        ! *** These parameters are read in for every series but only the last is actually used
        if( IASERVER < 7300 )then
          read(1,*,IOSTAT = ISO) IASWRAD, REVC, RCHC, SWRATNF, SWRATNS, FSWRATF, TEMTHKO, TEMBO, HTBED1, HTBED2
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE ASER.INP')
        endif

        do M = 1,TSATM(NA).NREC
          read(1,*,IOSTAT = ISO) TSATM(NA).TIM(M), (TSATM(NA).VAL(M,I),I = 1,7)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE ASER.INP')
        enddo

        ! *** APPLY CONVERSIONS
        do M = 1,TSATM(NA).NREC
          TSATM(NA).TIM(M)  = TSATM(NA).TIM(M)+TSATM(NA).TOFFSET
          TSATM(NA).VAL(M,4)= RAINCVT*TSATM(NA).VAL(M,4)
          TSATM(NA).VAL(M,5)= EVAPCVT*TSATM(NA).VAL(M,5)
          TSATM(NA).VAL(M,6)= SOLRCVT*TSATM(NA).VAL(M,6)
          TSATM(NA).VAL(M,7)=  CLDCVT*TSATM(NA).VAL(M,7)
        enddo
      enddo
      close(1)
    endif

    call Broadcast_Scalar(IASERVER, master_id)

    if( IASERVER < 7300 )then
      call Broadcast_Scalar(IASWRAD, master_id)
      call Broadcast_Scalar(REVC,    master_id)
      call Broadcast_Scalar(RCHC,    master_id)
      call Broadcast_Scalar(SWRATNF, master_id)
      call Broadcast_Scalar(SWRATNS, master_id)
      call Broadcast_Scalar(FSWRATF, master_id)
      call Broadcast_Scalar(TEMTHKO, master_id)
      call Broadcast_Scalar(TEMBO,   master_id)
      call Broadcast_Scalar(HTBED1,  master_id)
      call Broadcast_Scalar(HTBED2,  master_id)
    endif

    do NA = 1,NASER
      call Broadcast_Scalar(IRELH(NA)        , master_id)
      call Broadcast_Scalar(TSATM(NA).NREC   , master_id)
      call Broadcast_Scalar(TSATM(NA).TMULT   , master_id)
      call Broadcast_Scalar(TSATM(NA).TOFFSET , master_id)
      call Broadcast_Array(TSATM(NA).TIM     , master_id)
      call Broadcast_Array(TSATM(NA).VAL     , master_id)
    enddo

    ! *** RHO    = 1000.0  Density (kg / m^3)
    ! *** CP     = 4179.0  Specific Heat (J / kg / degC)
    ! *** 0.2393E-6 = 1/RHO/CP
  endif

  if( NASER > 1 )then
    call AllocateDSI(R2D_Global, LCM_Global, NASER, 0.0)

    if( process_id == master_id )then
      write(*,'(A)')'  READING ATMMAP.INP'
      open(1,FILE = 'atmmap.inp',STATUS = 'UNKNOWN')
      STR = READSTR(1)

      read(1,*) i
      if( I /= NATMMAP ) CALL STOPP('BAD NATMMAP')      ! *** BROADCAST IN SCAN_EFDC
    endif

    do NA = 1,NATMMAP
      if( process_id == master_id )then
        read(1,*) TATMMAPBEG(NA), TATMMAPEND(NA)
        STR = READSTR(1)

        do L = 2,LA_Global
          read(1,*) ID,JD,(atmwht_temp(nn),NN = 1,NASER)

          LG = LIJ_GLOBAL(ID,JD)
          R2D_Global(LG,:) = atmwht_temp(:)
        enddo
      endif

      call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          ATMWHT(:,L,NA) = R2D_Global(LG,:)
        endif
      enddo

    enddo

    ! *** send to all processes
    call Broadcast_Array(TATMMAPBEG, master_id)
    call Broadcast_Array(TATMMAPEND, master_id)

    if( process_id == master_id )then
      close(1)
    endif

    deallocate(R2d_Global)
  endif

  ! *** Read in bed temperatures from the file
  TEMB(1) = ABS(TEMBO)
  TEMB(LC) = TEMB(1)

  TEMB1(1)  = TEMB(1)
  TEMB1(LC) = TEMB(1)
  if( ISRESTI == 0 )then
    do L = 2,LA
      TEMB(L)  = TEMB(1)
      TEMB1(L) = TEMB(1)
    enddo
  endif

  ! *** Read in wind series from the file WSER.INP
  do L = 2,LA
    WINDST(L) = 0.
    TSX(L) = 0.
    TSY(L) = 0.
  enddo

  if( NWSER > 0 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING WSER.INP'
      open(1,FILE = 'wser.inp',STATUS = 'UNKNOWN')

      ! *** PERIOD TO TURN OFF WIND SHEAR (DEPRECATED IN EE7.3 SINCE ADDITION OF ISICE)
      write(*,'(A,I5)')  '  NUMBER OF WIND SERIES = ',NWSER

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      call SKIPCOM(1,'*')

      ! *** LOOP OVER EACH TIME SERIES
      STR = READSTR(1)

      do NW = 1,NWSER
        read(1,*,IOSTAT = ISO) M, TSWND(NW).TMULT, TSWND(NW).TOFFSET, WINDSCT, ISWDINT, WINDH(NW)
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE WSER.INP')

        if( WINDH(NW) <= 0.1) WINDH(NW) = 2.0     ! *** FIXING THE WIND SPEED MEASUREMENT HEIGHT (I.E. use THE WIND SPEED AS ENTERED)
        write(*,'(A,I5,I10,F6.1)')'  SERIES, NPTS, ANEMOMETER HEIGHT (m) = ',NW,TSWND(NW).NREC,WINDH(NW)

        do M = 1,TSWND(NW).NREC
          read(1,*,IOSTAT = ISO) TSWND(NW).TIM(M), (TSWND(NW).VAL(M,I),I = 1,2) !TWSER(M,NW),WINDS(M,NW),WINDD(M,NW)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE WSER.INP')
        enddo

        ! *** APPLY FACTORS AND OFFSETS AND COMPLETE TSWND DATA STRUCTURE
        do M = 1,TSWND(NW).NREC
          TSWND(NW).TIM(M) = TSWND(NW).TIM(M)+TSWND(NW).TOFFSET
        enddo
        if( ISWDINT <= 1 )then
          do M = 1,TSWND(NW).NREC
            TSWND(NW).VAL(M,1) = WINDSCT*TSWND(NW).VAL(M,1)
          enddo
        endif
        if( ISWDINT == 1 )then
          do M = 1,TSWND(NW).NREC
            if( TSWND(NW).VAL(M,2) <= 180.0 )then
              TSWND(NW).VAL(M,2) = TSWND(NW).VAL(M,2) + 180.
              if( TSWND(NW).VAL(M,2) == 360.) TSWND(NW).VAL(M,2) = 0.
            else
              TSWND(NW).VAL(M,2) = TSWND(NW).VAL(M,2) - 180.
              if( TSWND(NW).VAL(M,2) == 360.0) TSWND(NW).VAL(M,2) = 0.
            endif
          enddo
        endif
        if( ISWDINT == 2 )then
          do M = 1,TSWND(NW).NREC
            TSWND(NW).VAL(M,1) = WINDSCT*TSWND(NW).VAL(M,1)
            TSWND(NW).VAL(M,2) = WINDSCT*TSWND(NW).VAL(M,2)
          enddo
        endif
      enddo
      close(1)
    endif

    do NW = 1,NWSER
      call Broadcast_Scalar(WINDH(NW)         , master_id)
      call Broadcast_Scalar(TSWND(NW).NREC    , master_id)
      call Broadcast_Scalar(TSWND(NW).TMULT   , master_id)
      call Broadcast_Scalar(TSWND(NW).TOFFSET , master_id)
      call Broadcast_Array(TSWND(NW).TIM      , master_id)
      call Broadcast_Array(TSWND(NW).VAL      , master_id)
    enddo

  endif

  if( NWSER > 1 )then
    call AllocateDSI(R2D_Global, LCM_Global, NWSER, 0.0)

    if( process_id == master_id )then
      write(*,'(A)')'  READING WNDMAP.INP'
      open(1,FILE = 'wndmap.inp',STATUS = 'UNKNOWN', SHARED)

      STR = READSTR(1)
      read(1,*) i
      if( I /= NWNDMAP ) CALL STOPP('BAD NWNDMAP')      ! *** BROADCAST IN SCAN_EFDC
    endif

    do NW = 1,NWNDMAP

      if( process_id == master_id )then
        read(1,*) TWNDMAPBEG(NW), TWNDMAPEND(NW)
        STR = READSTR(1)

        do L = 2,LA_Global
          read(1,*) ID,JD,(WNDWHT_TEMP(nn),NN = 1,NWSER)

          LG = LIJ_GLOBAL(ID,JD)
          R2D_Global(LG,:) = WNDWHT_TEMP(:)
        enddo
      endif

      call Broadcast_Scalar(TWNDMAPBEG(NW), master_id)
      call Broadcast_Scalar(TWNDMAPEND(NW), master_id)
      call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
          WNDWHT(:,L,NW) = R2D_Global(LG,:)
        endif
      enddo
    enddo

    if( process_id == master_id )then
      close(1)
    endif
    deallocate(R2d_Global)
  endif

  ! *** Read in externally specified ice cover from the file ISER.INP
  if( ISICE == 1 .and. NISER >= 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING ISER.INP'
      open(1,FILE = 'iser.inp')
      STR = READSTR(1)

      do NS = 1,NISER
        RICECOVT(NS) = 0.
        RICETHKT(NS) = 0.
        MITLAST(NS)  = 2

        read(1,*,IOSTAT = ISO) M,TSICE(NS).TMULT,TOFFSET,RMULADJCOV
        if( ISO > 0 ) CALL STOPP('ISER.INP: READING ERROR')

        do M = 1,TSICE(NS).NREC
          ! *** TSICE(NS).VAL(M,1) ice on/off flag.
          ! *** TSICE(NS).VAL(M,2) is ice thickness, a legacy approach that was never used in EFDC
          read(1,*,IOSTAT = ISO) TSICE(NS).TIM(M),TSICE(NS).VAL(M,1)

          if( TSICE(NS).VAL(M,1) > 0.0 )then    ! *** VAL(M,2) is not used, for display only
            TSICE(NS).VAL(M,2) = RICETHK0       ! *** VAL(M,2) is not used, for display only
          else                                  ! *** VAL(M,2) is not used, for display only
            TSICE(NS).VAL(M,2) = 0.0            ! *** VAL(M,2) is not used, for display only
          endif                                 ! *** VAL(M,2) is not used, for display only
          if( ISO > 0 ) CALL STOPP('ISER.INP: READING ERROR')
        enddo

        do M = 1,TSICE(NS).NREC
          TSICE(NS).TIM(M)   = TSICE(NS).TIM(M) + TOFFSET
          TSICE(NS).VAL(M,1) = RMULADJCOV*TSICE(NS).VAL(M,1)
        enddo

      enddo

      close(1)
    endif

  elseif( ISICE == 2 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING ISTAT.INP'
      open(1,FILE = 'istat.inp')
      STR = READSTR(1)
      NS = 1    ! *** ALWAYS NISER = 1
      MITLAST(NS) = 2

      read(1,*,IOSTAT = ISO) M,TSICE(NS).TMULT,TOFFSET
      if( ISO > 0 ) CALL STOPP('ISTAT.INP: READING ERROR')

      do M = 1,TSICE(NS).NREC
        read(1,*,IOSTAT = ISO) TSICE(NS).TIM(M), TSICE(NS).VAL(M,1)  ! *** Fraction of ice cover: 0 to 1
        if( ISO > 0 ) CALL STOPP('ISTAT.INP: READING ERROR')
        TSICE(NS).TIM(M) = TSICE(NS).TIM(M) + TOFFSET
        if( TSICE(NS).VAL(M,1) > 1.0 ) TSICE(NS).VAL(M,1) = 1.0
      enddo
    endif

    ! *** send to all other processes
    call Broadcast_Array(MITLAST, master_id)
    do NS = 1, NISER
      call Broadcast_Scalar(TSICE(NS).NREC    , master_id)
      call Broadcast_Scalar(TSICE(NS).TMULT   , master_id)

      call Broadcast_Array(TSICE(NS).TIM, master_id)
      call Broadcast_Array(TSICE(NS).VAL, master_id)
    enddo
  endif

  if( ISICE == 1 .and. NISER > 1 )then
    ! *** Read ice time series weighting
    if(process_id == master_id )then
      write(*,'(A)')'READING ICEMAP.INP'
      open(1,FILE = 'icemap.inp')
      STR = READSTR(1)
      read(1,*) NICEMAP

      do NI = 1,NICEMAP
        read(1,*)TICEMAPBEG(NI),TICEMAPEND(NI)
        STR = READSTR(1)
        do L = 2,LA_Global
          read(1,*) ID,JD,(RICEWHT_global(NI,L,N),N = 1,NISER)
        enddo
      enddo
    endif

    call Broadcast_Array(TICEMAPBEG, master_id)
    call Broadcast_Array(TICEMAPEND, master_id)
    call Broadcast_Array(RICEWHT_global, master_id)

    ! *** Map to local
    do LG = 2, LA_GLOBAL
      L = MAP2LOCAL(LG).LL
      if( L > 1 )then
        RICEWHT(:,L,:) = RICEWHT_GLOBAL(:,LG,:)
      endif
    enddo

    if( process_id == master_id )then
      close(1)
    endif
  endif

  ! *** End of ICE *******************************

  ! *** Read in shell fish larave behavior data from the file SFBSER.INP
  if( ISTRAN(4) >= 1 )then

    if( process_id == master_id )then
      write(*,'(A)')'READING SFBSER.INP'
      open(1,FILE = 'sfbser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      read(1,*,IOSTAT = ISO) MSFSER,TCSFSER,TASFSER,TSRSF,TSSSF,ISSFLDN,ISSFLFE,SFLKILL
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFBSER.INP')
      do M = 1,MSFSER
        read(1,*,IOSTAT = ISO) TSFSER(M),RKDSFL(M),WSFLST(M),WSFLSM(M),DSFLMN(M),DSFLMX(M),SFNTBE(M),SFATBT(M)
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFBSER.INP')
      enddo
      close(1)
    endif

    call Broadcast_Scalar(MSFSER,  master_id)
    call Broadcast_Scalar(TCSFSER, master_id)
    call Broadcast_Scalar(TASFSER, master_id)
    call Broadcast_Scalar(TSRSF,   master_id)
    call Broadcast_Scalar(TSSSF,   master_id)
    call Broadcast_Scalar(ISSFLDN, master_id)
    call Broadcast_Scalar(ISSFLFE, master_id)
    call Broadcast_Scalar(SFLKILL, master_id)

    call Broadcast_Array(TSFSER, master_id)
    call Broadcast_Array(RKDSFL, master_id)
    call Broadcast_Array(WSFLST, master_id)
    call Broadcast_Array(WSFLSM, master_id)
    call Broadcast_Array(DSFLMN, master_id)
    call Broadcast_Array(DSFLMX, master_id)
    call Broadcast_Array(SFNTBE, master_id)
    call Broadcast_Array(SFATBT, master_id)

  endif

  ! ********************************************************************************************************************
  ! *** Read vegetation data from VEGE.INP and VEGSER.INP

  if( ISVEG > 0 .and. MACDRAG == 0 )then
    if( process_id == master_id )then
      write(*,'(A)')'READING VEGE.INP'
      open(1,FILE = 'vege.inp',STATUS = 'UNKNOWN',SHARED)
      STR = READSTR(1)
      read(1,*) MVEGTYP, MVEGOW, NVEGSER, UVEGSCL
      if( UVEGSCL <= 0. ) UVEGSCL = 1.E-12

      ! *** RDLPSQ - Stem density (#/m2)
      ! *** BPVEG  - Stem diameter (m)
      ! *** HPVEG  - Height of vegetation (m)
      ! *** ALPVEG - Form factor (dimensionless)
      ! *** SCVEG  - Drag coefficient (dimensionless)
      do M = 1,MVEGTYP
        read(1,*,ERR = 3120)IDUM, NVEGSERV(M), RDLPSQ(M), BPVEG(M), HPVEG(M), ALPVEG(M), SCVEG(M)
        BDLTMP   = BPVEG(M)*BPVEG(M)*RDLPSQ(M)          ! *** "ad" - dimensionless population density
        PVEGZ(M) = MAX(1.E-18,(1. - ALPVEG(M)*BDLTMP))  ! *** Turbulence/canopy adjustment factor to increase drag and turbulence. 0-No adjustment, 1-Infinite drag (invalid) (dimensionless)
        BDLPSQ(M) = BPVEG(M)*RDLPSQ(M)                  ! *** "a" - projected plant area per unit volume  (1/m)
      enddo
      close(1)
    endif !***end on master

    call Broadcast_Scalar(MVEGTYP, master_id)
    call Broadcast_Scalar(MVEGOW , master_id)
    call Broadcast_Scalar(NVEGSER, master_id)
    call Broadcast_Scalar(UVEGSCL, master_id)

    call Broadcast_Array(PVEGZ,    master_id)
    call Broadcast_Array(BDLPSQ ,  master_id)
    call Broadcast_Array(NVEGSERV, master_id)
    call Broadcast_Array(RDLPSQ  , master_id)
    call Broadcast_Array(BPVEG   , master_id)
    call Broadcast_Array(HPVEG   , master_id)
    call Broadcast_Array(ALPVEG  , master_id)
    call Broadcast_Array(SCVEG   , master_id)

    if( NVEGSER > 0 )then
      do M = 1,NVEGSER
        MVEGTLAST(M) = 1
      enddo
    endif


  GOTO 3124
  
3120 call STOPP('READ ERROR FOR FILE VEGE.INP')
3122 call STOPP('READ ERROR FOR FILE MHK.INP')

3124 continue

  if( NVEGSER >= 1 )then
    if(process_id == master_id )then
      write(*,'(A)')'READING VEGSER.INP'
      open(1,FILE = 'vegser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      do NS = 1,NVEGSER
        read(1,*,IOSTAT = ISO) MVEGSER(NS),TCVEGSER(NS),TAVEGSER(NS)
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE VEGSER.INP')
        do M = 1,MVEGSER(NS)
          read(1,*,IOSTAT = ISO)TVEGSER(M,NS),VEGSERR(M,NS),VEGSERB(M,NS),VEGSERH(M,NS)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE VEGSER.INP')
          TVEGSER(M,NS) = TVEGSER(M,NS)+TAVEGSER(NS)
        enddo
      enddo
      close(1)
    endif

    call Broadcast_Array(MVEGSER , master_id)
    call Broadcast_Array(TCVEGSER, master_id)
    call Broadcast_Array(TAVEGSER, master_id)
    call Broadcast_Array(TVEGSER,  master_id)
    call Broadcast_Array(VEGSERR,  master_id)
    call Broadcast_Array(VEGSERB,  master_id)
    call Broadcast_Array(VEGSERH,  master_id)

    ! *** REINITIALIZE CLASSES HAVING TIME SERIES INFORMATION
    do M = 1,MVEGTYP
      if( NVEGSERV(M) > 0 )then
        NS = NVEGSERV(M)
        RDLPSQ(M) = VEGSERR(1,NS)
        BPVEG(M) = VEGSERB(1,NS)
        HPVEG(M) = VEGSERH(1,NS)
        BDLTMP = BPVEG(M)*BPVEG(M)*RDLPSQ(M)
        PVEGZ(M) = 1.-ALPVEG(M)*BDLTMP
        BDLPSQ(M) = BPVEG(M)*RDLPSQ(M)
      endif
    enddo
  endif
  GOTO 7122
7120 call STOPP('READ ERROR FOR FILE VEGSER.INP')
7122 continue

  ! ********************************************************************************************************************
    ! *** Read Marine Hydrokinetics (MHK) file
    if( LMHK )then
      if( process_id == master_id )then
        write(*,'(A)')'READING MHK.INP'
        open(1,FILE = 'mhk.inp',STATUS = 'UNKNOWN')

        STR = READSTR(1)

        read(1,*,ERR = 3122)MHKTYP,NFLAGPWR,UPSTREAM,OUTPUTFLAG
      endif ! *** end on master

      ! *** send to all
      call Broadcast_Scalar(MHKTYP,    master_id)
      call Broadcast_Scalar(NFLAGPWR,  master_id)
      call Broadcast_Scalar(UPSTREAM,  master_id)
      call Broadcast_Scalar(OUTPUTFLAG,master_id)


      allocate(BOFFMHK(MHKTYP),BOFFSUP(MHKTYP))
      allocate(TOFFMHK(MHKTYP),TOFFSUP(MHKTYP))

      if( process_id == master_id )then

        if( NFLAGPWR == 1 )then
          do M = 1,MHKTYP
            read(1,*,ERR = 3122)WIDTHMHK(M),WIDTHSUP(M), BOFFMHK(M),BOFFSUP(M),TOFFMHK(M),TOFFSUP(M),CTMHK(M),CDSUP(M),VMINCUT(M),VMAXCUT(M),DENMHK(M)
            CTMHK(M) = CTMHK(M)*DENMHK(M)
            CDSUP(M) = CDSUP(M)*DENMHK(M)

            do L = 2,LA
              if( M+90 == MVEGL(L) )then
                ZMINMHK(M,L) = BELV(L)+BOFFMHK(M)
                ZMAXMHK(M,L) = BELV(L)+TOFFMHK(M)
                ZMINSUP(M,L) = BELV(L)+BOFFSUP(M)
                ZMAXSUP(M,L) = BELV(L)+TOFFSUP(M)
                DIAMMHK = ZMAXMHK(M,L)-ZMINMHK(M,L)
                if( DIAMMHK<0.0 )then   !error check
                  PRINT *,'MHK ZMIN > ZMAX'
                  call STOPP('.')
                endif
              endif
            enddo
          enddo
          read(1,*) !skip the header line
          read(1,*)BETAMHK_D,BETAMHK_P,CE4MHK,PB_COEF
        endif
        ! *** send to all processes
        call Broadcast_Array(WIDTHMHK, master_id)
        call Broadcast_Array(WIDTHSUP, master_id)
        call Broadcast_Array(BOFFMHK,  master_id)
        call Broadcast_Array(BOFFSUP,  master_id)
        call Broadcast_Array(TOFFMHK,  master_id)
        call Broadcast_Array(TOFFSUP,  master_id)
        call Broadcast_Array(CTMHK,    master_id)
        call Broadcast_Array(CDSUP,    master_id)
        call Broadcast_Array(VMINCUT,  master_id)
        call Broadcast_Array(VMAXCUT,  master_id)
        call Broadcast_Array(DENMHK,   master_id)

        call Broadcast_Array(ZMINMHK, master_id)
        call Broadcast_Array(ZMAXMHK, master_id)
        call Broadcast_Array(ZMINSUP, master_id)
        call Broadcast_Array(ZMAXSUP, master_id)

        call Broadcast_Scalar(BETAMHK_D,master_id)
        call Broadcast_Scalar(BETAMHK_P,master_id)
        call Broadcast_Scalar(CE4MHK,   master_id)
        call Broadcast_Scalar(PB_COEF,  master_id)

      elseif( NFLAGPWR == 2 )then
        PRINT*,'Not available yet'

      elseif( NFLAGPWR == 3 )then !FFP input style
        if(process_id == master_id )then
          read(1,*,ERR = 3122)WIDTHMHK(M),WIDTHSUP(M),BOFFMHK(M),HEIGHTMHK(M),HEIGHTSUP(M),REFELEV(M),CTMHK(M),CDSUP(M),VMINCUT(M),VMAXCUT(M),DENMHK(M)
          CTMHK(M) = CTMHK(M)*DENMHK(M)
          CDSUP(M) = CDSUP(M)*DENMHK(M)
          do L = 2,LA
            if( M+90 == MVEGL(L) )then
              ZMINMHK(M,L) = BELV(L)+REFELEV(M)
              ZMAXMHK(M,L) = BELV(L)+REFELEV(M)+HEIGHTMHK(M)
              ZMINSUP(M,L) = BELV(L)
              ZMAXSUP(M,L) = BELV(L)+REFELEV(M)+HEIGHTSUP(M)
              DIAMMHK = ZMAXMHK(M,L)-ZMINMHK(M,L)
              if( DIAMMHK<0.0 )then
                PRINT*,'MHK ZMIN > ZMAX'
                call STOPP('')
              endif
            endif
          enddo

          read(1,*) !skip the header line
          read(1,*)BETAMHK_D,BETAMHK_P,CE4MHK,PB_COEF
        endif

        call Broadcast_Array(WIDTHMHK, master_id)
        call Broadcast_Array(WIDTHSUP, master_id)
        call Broadcast_Array(BOFFMHK,  master_id)
        call Broadcast_Array(HEIGHTMHK,master_id)
        call Broadcast_Array(HEIGHTSUP,master_id)
        call Broadcast_Array(REFELEV,  master_id)
        call Broadcast_Array(CTMHK,    master_id)
        call Broadcast_Array(CDSUP,    master_id)
        call Broadcast_Array(VMINCUT,  master_id)
        call Broadcast_Array(VMAXCUT,  master_id)
        call Broadcast_Array(DENMHK,   master_id)

        call Broadcast_Array(ZMINMHK, master_id)
        call Broadcast_Array(ZMAXMHK, master_id)
        call Broadcast_Array(ZMINSUP, master_id)
        call Broadcast_Array(ZMAXSUP, master_id)

        call Broadcast_Scalar(BETAMHK_D,master_id)
        call Broadcast_Scalar(BETAMHK_P,master_id)
        call Broadcast_Scalar(CE4MHK,   master_id)
        call Broadcast_Scalar(PB_COEF,  master_id)

      endif

      if(process_id == master_id )then
        close(1)
      endif

      do L = 2,LA
        if( MVEGL(L)>90 )then
          if( (DIAMMHK>DXP(L) .or. DIAMMHK>DYP(L)) .and. DENMHK(MVEGL(L)-90)>1.0 )then
            PRINT*,'MHK DIAMETER EXCEEDS CELL SIZE'
            PRINT*,'AND DENSITY >= 1'
            call STOPP('')
          endif
        endif
      enddo

    endif
    !!!END SCJ BLOCK

  endif

  ! ********************************************************************************************************************
  ! *** Read bank erosion map and time series files
  if( (ISTRAN(6)>0 .or. ISTRAN(7) > 0 ) .and. ISBKERO >= 1 )then
    if(process_id == master_id )then
      write(*,'(A)')'READING BEMAP.INP'
      open(1,FILE = 'bemap.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      read(1,*)NBEPAIR,NBESER

      do NS = 1,NBEPAIR
        read(1,*)IBANKBE(NS),JBANKBE(NS),ICHANBE(NS),JCHANBE(NS),NBESERN(NS),FBESER(NS)
      enddo

      close(1)
    endif
    ! *** send to all
    call Broadcast_Scalar(NBEPAIR,master_id)
    call Broadcast_Scalar(NBESER, master_id)

    call Broadcast_Array(IBANKBE, master_id)
    call Broadcast_Array(JBANKBE, master_id)
    call Broadcast_Array(ICHANBE, master_id)
    call Broadcast_Array(JCHANBE, master_id)
    call Broadcast_Array(NBESERN, master_id)
    call Broadcast_Array(FBESER , master_id)

  endif


  if( (ISTRAN(6)>0 .or. ISTRAN(7) > 0 ) .and. ISBKERO >= 1 .and. NBESER > 0 )then
    if(process_id == master_id )then
      write(*,'(A)')'READING BESER.INP'
      open(1,FILE = 'beser.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      do NS = 1,NBESER

        read(1,*)MBESER(NS),TCBESER(NS),TABESER(NS),RMULADJ,ADDADJ
        MBETLAST(NS) = 1

        do M = 1,MBESER(NS)
          read(1,*)TBESER(M,NS),BESER(M,NS),FWCBESER(M,NS)
          TBESER(M,NS) = TBESER(M,NS)+TABESER(NS)
          BESER(M,NS) = RMULADJ*(BESER(M,NS)+ADDADJ)
        enddo

      enddo

      close(1)
    endif
    ! *** Send to all
    call Broadcast_Array(MBESER, master_id)
    call Broadcast_Array(TCBESER, master_id)
    call Broadcast_Array(TABESER, master_id)
    call Broadcast_Array(TBESER, master_id)
    call Broadcast_Array(BESER, master_id)
    call Broadcast_Array(FWCBESER, master_id)

  endif

  ! ********************************************************************************************************************
  ! *** READ ZONALLY VARYING SEDIMENT BED PARTICLE MIXING
  !----------------------------------------------------------------------C
  ITMPPMX = 0
  do NT = 1,NTOX
    if( ISPMXZ(NT) == 1 )ITMPPMX = 1
  enddo

  if( ISTRAN(5) > 0 .and. ITMPPMX == 1 )then
    if(process_id == master_id )then
      write(*,'(A)')'READING PARTMIX.INP'
      open(1,FILE = 'partmix.inp')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      !#######################################################################
      !     HQI change to input multiplication scale factor for particle mixing rate
      !     RM 10/06/05
      read(1,*)NPMXZ,NPMXPTS,PMIXSF

      !#######################################################################
      do NZ = 1,NPMXZ
        do NP = 1,NPMXPTS
          read(1,*) PMXDEPTH(NP,NZ), PMXCOEF(NP,NZ)
          !#######################################################################
          !     HQI change to input multiplication scale factor for particle mixing rate
          !     RM 10/06/05
          PMXCOEF(NP,NZ) = PMIXSF*PMXCOEF(NP,NZ)
          !#######################################################################
        enddo
      enddo

      close(1)
    endif

    ! *** send to all processes
    call Broadcast_Scalar(NPMXZ,  master_id)
    call Broadcast_Scalar(NPMXPTS,master_id)
    call Broadcast_Scalar(PMIXSF, master_id)

    call Broadcast_Array(PMXDEPTH, master_id)
    call Broadcast_Array(PMXCOEF , master_id)

    if(process_id == master_id )then
      write(*,'(A)')'READING PMXMAP.INP'
      open(1,FILE = 'pmxmap.inp')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)

      do L = 2,LC_Global-1
        read(1,*)LDUM,IDUM,JDUM,LPMXZ_Global(L)
      enddo

      close(1)
    endif
    ! *** send to all processes
    call Broadcast_Array(LPMXZ_Global, master_id)

    do LG = 2 , LC_GLOBAL - 1
      L = MAP2LOCAL(LG).LL
      if( L > 1 )then
        LPMXZ(L) = LPMXZ_GLOBAL(LG)
      endif
    enddo

  endif

  ! ********************************************************************************************************************
  ! *** Spatially/temporaly varying wind fields
  call ReadFields()     ! 2018-10-12, NTL: Read time & space varying fields
  call ReadCyclones()   ! 2021-03-31, NTL: Read cyclone tracks

  END

  ! ********************************************************************************************************************
  ! *** DSI UTILITIES
  FUNCTION PARSE_REAL(INLINE)

  use GLOBAL,only: IK4
  implicit none
  integer(IK4) :: I1,I2,ILEN,IPOS,ILEN2
  real :: PARSE_REAL
  character*(*) INLINE
  character*15  CVAL

  ILEN = LEN_TRIM(INLINE)
  PARSE_REAL = 0.
  do I1 = 1,ILEN
    if( INLINE(I1:I1) == ':' )then
      do IPOS = I1+1,ILEN
        if( INLINE(IPOS:IPOS) /= ' ')EXIT
      enddo
      if( IPOS>ILEN) return
      CVAL = INLINE(IPOS:ILEN)
      ILEN2 = LEN_TRIM(CVAL)
      do I2 = 1,ILEN2
        if( CVAL(I2:I2) == ' ' .or. CVAL(I2:I2) == ',' .or. I2 == ILEN2 )then
          read(CVAL(1:I2),'(F12.1)',ERR = 999)PARSE_REAL
          return
        endif
      enddo
    endif
  enddo
999 call STOPP(' ERROR PARSING REAL')

  END FUNCTION

  FUNCTION PARSE_LOGICAL(INLINE)

  use GLOBAL,only: IKV
  implicit none

  integer :: ILEN,IC,JC,IPOS,ILEN2
  character*(*) INLINE
  character*12  CVAL
  logical       PARSE_LOGICAL

  ILEN = LEN_TRIM(INLINE)
  do IC = 1,ILEN
    if( INLINE(IC:IC) == ':' )then
      do IPOS = IC+1,ILEN
        if( INLINE(IPOS:IPOS) /= ' ')EXIT
      enddo
      if( IPOS>ILEN ) return
      CVAL = INLINE(IPOS:ILEN)
      ILEN2 = LEN_TRIM(CVAL)
      do JC = 1,ILEN2
        if( CVAL(JC:JC) == ' ' .or. CVAL(JC:JC) == ',' .or. JC == ILEN2 )then
          if( CVAL(1:1) == 'T' .or. CVAL(1:1) == 'Y' )then
            PARSE_LOGICAL = .TRUE.
            return
          endif
        endif
      enddo
    endif
  enddo
  PARSE_LOGICAL = .FALSE.

  END FUNCTION

  SUBROUTINE PARSESTRING(string,substring)
  character(*),intent(inout) :: string
  character(*) :: substring
  integer :: posBT(2), pos

  posBT(1) = index(string,' ')      ! Blank
  posBT(2) = index(string,char(9))  ! Tab
  pos = minval( posBT, posBT > 0 )
  substring = string(1:pos)
  string = adjustl(string(pos+1:))
  end SUBROUTINE

  SUBROUTINE READ_SUBSET
  use GLOBAL
  use INFOMOD ,ONLY:SKIPCOM, READSTR
  use XYIJCONV,ONLY:XY2IJ
  use Variables_MPI
  use Broadcast_Routines

  implicit none
  character(256) :: string
  character(20) :: substring
  integer :: ISO,NS,IS,NP

  If( process_id == master_id )then
    open(1,FILE = 'subset.inp',action = 'read')
    call SKIPCOM(1,'*')
    read(1,*,IOSTAT = ISO) NSUBSET
    if( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
  endif

  call Broadcast_Scalar(NSUBSET, master_id)

  allocate(HFREGRP(NSUBSET))
  allocate(IJHFRE(NSUBSET),HFREDAYBG(NSUBSET),HFREDAYEN(NSUBSET))
  allocate(HFREDUR(NSUBSET),HFREDAY(NSUBSET),HFREMIN(NSUBSET),NPNT(NSUBSET))

  do NS = 1,NSUBSET
    If( process_id == master_id )then
      call SKIPCOM(1,'*')
      read(1,*,IOSTAT = ISO) IS,IJHFRE(IS),HFREDAYBG(IS),HFREDUR(IS),HFREMIN(IS),NPNT(IS)
      if( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
    endif

    call Broadcast_Scalar(IS,            master_id)
    call Broadcast_Scalar(NPNT(IS),      master_id)

    ALLOCATE (HFREGRP(IS).ICEL(NPNT(IS)),HFREGRP(IS).JCEL(NPNT(IS)))
    ALLOCATE (HFREGRP(IS).XCEL(NPNT(IS)),HFREGRP(IS).YCEL(NPNT(IS)))
    ALLOCATE (HFREGRP(IS).NAME(NPNT(IS)))

    If( process_id == master_id )then
      do NP = 1,NPNT(IS)
        call SKIPCOM(1,'*')
        read(1,'(a)') string
        string = adjustl(trim(string))
        if( IJHFRE(IS) == 1 )then
          !READ(1,*,advance = 'no',IOSTAT = ISO) HFREGRP(IS)%ICEL(NP),HFREGRP(IS)%JCEL(NP)
          call PARSESTRING(string, substring)
          read(substring,*) HFREGRP(IS).ICEL(NP)
          call PARSESTRING(string, substring)
          read(substring,*) HFREGRP(IS).JCEL(NP)
        else
          !READ(1,*,advance = 'no',IOSTAT = ISO) HFREGRP(IS)%XCEL(NP),HFREGRP(IS)%YCEL(NP)
          call PARSESTRING(string, substring)
          read(substring,*) HFREGRP(IS).XCEL(NP)
          call PARSESTRING(string, substring)
          read(substring,*) HFREGRP(IS).YCEL(NP)
        endif
        if( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
        string = TRIM(string)
        if( string(1:1) == '!') string = string(2:20)
        HFREGRP(IS).NAME(NP) = TRIM(string)
      enddo
    endif

  enddo

  call Broadcast_Array(IJHFRE,    master_id)
  call Broadcast_Array(HFREDAYBG, master_id)
  call Broadcast_Array(HFREDUR,   master_id)
  call Broadcast_Array(HFREMIN,   master_id)


  do NS = 1,NSUBSET
    call Broadcast_Array(HFREGRP(NS)%ICEL, master_id)
    call Broadcast_Array(HFREGRP(NS)%JCEL, master_id)
    call Broadcast_Array(HFREGRP(NS)%XCEL, master_id)
    call Broadcast_Array(HFREGRP(NS)%YCEL, master_id)
  enddo

  If( process_id == master_id ) close(1)


  END SUBROUTINE
