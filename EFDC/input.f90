! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE INPUT(TITLE)

  ! *** SUBROUTINE INPUT READS ALL INPUT DATA EXCEPT DATA IN LXLY.INP,
  ! *** MASK.INP AND RESTART.INP

  !----------------------------------------------------------------------!
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3

  Use GLOBAL
  Use Variables_WQ
  Use RESTART_MODULE, ONLY:Setup_Continuation_Files
  Use INFOMOD ,ONLY:SKIPCOM, READSTR
  USE HYDSTRUCMOD
  Use SHELLFISHMOD, ONLY: NSF, ISFFARM, NSFCELLS, READ_SHELLFISH_JSON
  Use Variables_Propwash
  Use FIELDS
  Use CYCLONE
  Use Allocate_Initialize

  ! *** New for MPI
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  Use Broadcast_Routines
  ! *** End new for MPI

  Use IFPORT

  IMPLICIT NONE

  CHARACTER*80 TEXT, TITLE, STR*200
  CHARACTER*10 CDUM
  CHARACTER*3  NCARD
  CHARACTER    ADUMMY*5,RESTARTF*50,STRC*650

  REAL(RK4) :: SEEPRATE(1000)
  REAL,ALLOCATABLE,DIMENSION(:) :: RMULADS
  REAL,ALLOCATABLE,DIMENSION(:) :: ADDADS
  REAL,ALLOCATABLE,DIMENSION(:) :: BOFFMHK, BOFFSUP, TOFFMHK, TOFFSUP
  REAL,ALLOCATABLE,DIMENSION(:,:) :: PFX2
  REAL,ALLOCATABLE,DIMENSION(:,:) :: QSERSM,RULES

  INTEGER,ALLOCATABLE,DIMENSION(:) :: IPARTSP   ! *** NEW BEDFORD
  INTEGER,ALLOCATABLE,DIMENSION(:) :: IDX
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: IBLTMP

  INTEGER :: ITURB=0, IPMC, IS, NS, NX, NT, NW, NA, NI, NWR
  INTEGER :: L, LP, LW, LE, LS, LN, LG, LU, LD
  INTEGER :: IDUM, NDUM, LCM2T, I, J, K, M, KDUM, IMDXDY, NPFORN
  INTEGER :: ISO, ITIDASM, NPFOR, NPFORS, NPFORW, NPFORE, ITSSS
  INTEGER :: NDUM1, NDUM2, NTMP, ITYPE, IFLAG
  INTEGER :: MTMP
  INTEGER :: MMAX, MS, MMIN
  INTEGER :: JSFDCH, IDUMMY, NPP, MD, MU
  INTEGER :: MTSSS, IACROSS, JCTMP, JACROSS, JT, JF, JLAST, NMD
  INTEGER :: IFIRST, ILAST, IT, NP, LT, MVEGIJT, NMDXDY, INITTEMP
  INTEGER :: ITMP, JTMP, LTMP, LL, ITMPU, JTMPU, ITMPD, JTMPD, ID, JD, NSEEPCLASSES
  INTEGER :: IZONE, LDUM, JDUM, IVAL, ISALTYP, IREAD, KBINPUT, ISCOLD, ISEDINIT
  INTEGER :: ISTYP, ISMOOTH, ICHGQS, NCTMP
  INTEGER :: NFLAGPWR, ITMPPMX, NPMXZ, NPMXPTS, PMIXSF, NZ
  INTEGER :: NPBPH, IASERVER, NC
  INTEGER :: ISQCTRL(6)

  REAL*8 :: TMPDATE
  REAL   :: ADMAX, ADMIN, AHMAX, AHMIN
  REAL   :: DXIJ, DYIJ, HIJ, BELVIJ, ZBRIJ, RVALUE, PSERTMP, DIAMMHK
  REAL   :: DXYCVT, HADADJ, RAD, AMP, T1, T2, TMPAMP, TMPPHS, BOTTOM
  REAL   :: TOP1, TOP2 ! , QSSE,  Removed Zander for MPI Domain decomp
  REAL   :: SEDVDRT, DSTR, USTR, DUM, QSERTMP
  REAL   :: RMDX, RMDY, CVTFACX, CVTFACY, FBODY1, FBODY2, BDLTMP
  REAL   :: RMULADJ, ADDADJ, RMULADJS, ADDADJS, PSERTMPS, CSERTMP
  REAL   :: QCTLTMP, SOLRCVT, CLDCVT, WINDSCT, RMULADJCOV, RMULADJTHK
  REAL   :: TDELTA
  REAL,EXTERNAL :: SETSTVEL, PARSE_REAL
  LOGICAL :: BFLAG

  INTEGER :: RI_GLOBAL, RJ_GLOBAL
  INTEGER :: IGL, JGL
  Integer :: i_loc, j_loc, l_local, nn, ierr
  Integer, allocatable, dimension(:) :: LPMXZ_Global
  
  REAL,ALLOCATABLE,DIMENSION(:) :: CONINIT
  REAL,ALLOCATABLE,DIMENSION(:) :: CONBINIT

  Call AllocateDSI(RMULADS,  NSTM,   0.0)      
  Call AllocateDSI(ADDADS,   NSTM,   0.0)      
  Call AllocateDSI(QSERSM,   NDQSER, KCM, 0.0)
  Call AllocateDSI(PFX2,     NPFORM, MTM, 0.0)
  Call AllocateDSI(IPARTSP,  NTXM,     0)
  Call AllocateDSI(CONINIT,  KCM,    0.0)
  Call AllocateDSI(CONBINIT, KBM,    0.0)

  G=9.81
  PI=3.1415926535898
  PI2=2.*PI
2 FORMAT(A80)

  if( process_id == master_id )THEN
    ! *** READ MAIN INPUT FILE EFDC.INP
    WRITE(*,'(A)')'READING THE MAIN EFDC CONTROL FILE: EFDC.INP'
    OPEN(1,FILE='efdc.inp',STATUS='UNKNOWN')

    !1**  READ TITLE CARD
    NCARD='1'
    CALL SEEK('C1')
    READ(1,2) TITLE
    WRITE(7,1002)NCARD
    WRITE(7,2) TITLE
  end if
  Call Broadcast_Scalar(TITLE,     master_id)

  !C1A**  READ MODE OPTIONS
  NCARD='1A'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C1A')
    READ(1,*,IOSTAT=ISO) IS2TIM, IGRIDH, IGRIDV, KMINV, SGZHPDELTA
    
    WRITE(7,1002)NCARD
    WRITE(7,*) IS2TIM,IGRIDH,IGRIDV,KMINV,SGZHPDELTA
    IF( ISO > 0 ) GOTO 100
  end if

  Call Broadcast_Scalar(IS2TIM,     master_id)
  Call Broadcast_Scalar(IGRIDH,     master_id)
  Call Broadcast_Scalar(IGRIDV,     master_id)
  Call Broadcast_Scalar(KMINV,      master_id)
  Call Broadcast_Scalar(SGZHPDELTA, master_id)

  !C2**  READ RESTART AND DIAGNOSTIC SWITCHES
  ! *** ********************************************************
  if( process_id == master_id )THEN
    NCARD='2'
    CALL SEEK('C2')
    READ(1,*,IOSTAT=ISO) ISRESTI, ISRESTO, ISRESTR, ISGREGOR, ISLOG, ISDIVEX, ISNEGH, ISMMC, ISBAL, ICONTINUE, ISHOW
  end if

  Call Broadcast_Scalar(ISRESTI,  master_id)
  Call Broadcast_Scalar(ISRESTO,  master_id)
  Call Broadcast_Scalar(ISRESTR,  master_id)
  Call Broadcast_Scalar(ISGREGOR, master_id)
  Call Broadcast_Scalar(ISLOG,    master_id)
  Call Broadcast_Scalar(ISDIVEX,  master_id)
  Call Broadcast_Scalar(ISNEGH,   master_id)
  Call Broadcast_Scalar(ISMMC,    master_id)
  Call Broadcast_Scalar(ISBAL,    master_id)
  Call Broadcast_Scalar(ICONTINUE,master_id)
  Call Broadcast_Scalar(ISHOW,    master_id)

  ! *** HANDLE BATHYMETRY ADJUSTMENTS
  ISRESTIOPT = 0
  IF( ISRESTI == -1 )THEN
    ISRESTIOPT = 1
    ISRESTI = 1
  ENDIF

  IF( ISRESTI == 1 .AND. ICONTINUE == 1 )THEN
    ! ** FOR CONTINUATION MODE:
    if( process_id == master_id )THEN
      CALL SEEK('C2A')
      READ(1,*,IOSTAT=ISO) Restart_In_Ver
      READ(1,*,IOSTAT=ISO) RESTARTF
    end if
  ENDIF
  
  ! *** ********************************************************
  if( process_id == master_id )THEN
    WRITE(7,1002)NCARD
    WRITE(7,*) ISRESTI, ISRESTO, ISRESTR, ISGREGOR, ISLOG, ISDIVEX, ISNEGH, ISMMC, ISBAL, ICONTINUE, ISHOW
    IF( ISO > 0 ) GOTO 100
  End if

  Call Broadcast_Scalar(ISRESTI,  master_id)
  Call Broadcast_Scalar(ISRESTO,  master_id)
  Call Broadcast_Scalar(ISRESTR,  master_id)
  Call Broadcast_Scalar(ISGREGOR, master_id)
  Call Broadcast_Scalar(ISLOG,    master_id)
  Call Broadcast_Scalar(ISDIVEX,  master_id)
  Call Broadcast_Scalar(ISNEGH,   master_id)
  Call Broadcast_Scalar(ISMMC,    master_id)
  Call Broadcast_Scalar(ISBAL,    master_id)
  Call Broadcast_Scalar(ICONTINUE,master_id)
  Call Broadcast_Scalar(ISHOW,    master_id)

  IF( ISMMC < 0 )THEN
    DEBUG = .TRUE.
    ISMMC = 0
    if( process_id == master_id ) WRITE(*,'(A)') 'DEBUG ON'
  ELSE
    DEBUG = .FALSE.
    if( process_id == master_id ) WRITE(*,'(A)') 'DEBUG OFF'
  ENDIF

  !C3**  READ RELAXATION PARAMETERS
  NCARD='3'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C3')
    READ(1,*,IOSTAT=ISO) RP, RSQM, ITERM, IRVEC, IATMP, IWDRAG, ITERHPM, IDRYCK, ISDSOLV, FILT3TL
    WRITE(7,1002)NCARD
    WRITE(7,*) RP, RSQM, ITERM, IRVEC, IATMP, IWDRAG, ITERHPM, IDRYCK, ISDSOLV, FILT3TL
    IF( ISO > 0 ) GOTO 100
    IF( IRVEC /= 0 .AND. IRVEC /= 9 ) CALL STOPP('INVALID IRVEC')
  end if

  Call Broadcast_Scalar(RP     , master_id)
  Call Broadcast_Scalar(RSQM   , master_id)
  Call Broadcast_Scalar(ITERM  , master_id)
  Call Broadcast_Scalar(IRVEC  , master_id)
  Call Broadcast_Scalar(IATMP  , master_id)
  Call Broadcast_Scalar(IWDRAG , master_id)
  Call Broadcast_Scalar(ITERHPM, master_id)
  Call Broadcast_Scalar(IDRYCK , master_id)
  Call Broadcast_Scalar(ISDSOLV, master_id)
  Call Broadcast_Scalar(FILT3TL, master_id)

  !C4**  READ LONGTERM MASS TRANSPORT INTEGRATION ONLY SWITCHES
  NCARD='4'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C4')
    READ(1,*,IOSTAT=ISO) ISLTMT, ISSSMMT, RESSTEP
    WRITE(7,1002)NCARD
    WRITE(7,*) ISLTMT, ISSSMMT, RESSTEP
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ISLTMT , master_id)
  Call Broadcast_Scalar(ISSSMMT, master_id)
  Call Broadcast_Scalar(RESSTEP, master_id)

  !C5**  READ MOMENTUM ADVECTION AND DIFFUSION SWITCHES AND MISC
  NCARD='5'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C5')
    READ(1,*,IOSTAT=ISO) ISCDMA, ISHDMF, ISDISP, ISWASP, ISDRY, ICALTB, ISRLID, ISVEG, ISVEGL, ISITB, IHMDSUB, IINTPG
    WRITE(7,1002)NCARD
    WRITE(7,*) ISCDMA, ISHDMF, ISDISP, ISWASP, ISDRY, ICALTB, ISRLID, ISVEG, ISVEGL, ISITB, IHMDSUB, IINTPG
    ! *** Not Used: IDRYCK

    IF( ISO > 0 ) GOTO 100
    IF( MACDRAG > 0 .AND. ISVEG == 0 )  ISVEG = 2
  endif

  Call Broadcast_Scalar(ISCDMA , master_id)
  Call Broadcast_Scalar(ISHDMF , master_id)
  Call Broadcast_Scalar(ISDISP , master_id)
  Call Broadcast_Scalar(ISWASP , master_id)
  Call Broadcast_Scalar(ISDRY  , master_id)
  Call Broadcast_Scalar(ICALTB , master_id)
  Call Broadcast_Scalar(ISRLID , master_id)
  Call Broadcast_Scalar(ISVEG  , master_id)
  Call Broadcast_Scalar(ISVEGL , master_id)
  Call Broadcast_Scalar(ISITB  , master_id)
  Call Broadcast_Scalar(IHMDSUB, master_id)
  Call Broadcast_Scalar(IINTPG , master_id)
  IDRYTBP=0
  IF( ISDRY < 0 )THEN
    ISDRY=ABS(ISDRY)
    IDRYTBP=1
  ENDIF
  IF( ISVEG < 1 )ISITB=0
  IF( ISWASP == 99)ISICM=1
  IF( ISRLID == 1) ISDRY=-1
  IF( ISWASP == 10)ISRCA=1
  JSWAVE=0
  ! PMC      IS1DCHAN=0
  ! PMC  IF( ISCDMA == 10) IS1DCHAN=1                              Y
  ISCOSMIC=0

  !C6**  DISSOLVED AND SUSPENDED CONSTITUENT TRANSPORT SWITCHES
  NCARD='6'

  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C6')
    DO NS=0,8
      READ(1,*,IOSTAT=ISO) ISTRAN(NS), ISTOPT(NS), ISCDCA(NS), ISADAC(NS), ISFCT(NS), ISPLIT(NS), ISADAH(NS), ISADAV(NS), ISCI(NS), ISCO(NS)
      IF( ISCDCA(NS) >= 4) ISCOSMIC=1
      WRITE(7,1002)NCARD
      WRITE(7,*) ISTRAN(NS), ISTOPT(NS), ISCDCA(NS), ISADAC(NS), ISFCT(NS), ISPLIT(NS), ISADAH(NS), ISADAV(NS), ISCI(NS), ISCO(NS)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  end if

  Call Broadcast_Array(ISTRAN, master_id)
  Call Broadcast_Array(ISTOPT, master_id)
  Call Broadcast_Array(ISCDCA, master_id)
  Call Broadcast_Array(ISADAC, master_id)
  Call Broadcast_Array(ISFCT , master_id)
  Call Broadcast_Array(ISPLIT, master_id)
  Call Broadcast_Array(ISADAH, master_id)
  Call Broadcast_Array(ISADAV, master_id)
  Call Broadcast_Array(ISCI  , master_id)
  Call Broadcast_Array(ISCO , master_id)

  IF( ISTRAN(8) >= 1 .AND. ISTRAN(2) == 0 )THEN
    PRINT *,'*** WARNING: TEMPERATURE SHOULD BE ACTIVATED FOR WQ CALCULATIONS!'
    CALL SLEEPQQ(5000)
  ENDIF

  ! *** ********************************************************
  IF( ISRESTI == 1 .AND. ICONTINUE == 1 )THEN
    ! *** MODEL RUN CONTINUATION
    if( process_id == master_id )THEN
      CLOSE(1)
      CALL Setup_Continuation_Files(RESTARTF)           ! GENERATE RESTART FILES, OTHERWISE EFDC LOADS EXISTING ONES
      OPEN(1,FILE='efdc.inp',STATUS='UNKNOWN')
    endif
    Call Broadcast_Scalar(NRESTART,         master_id)
    Call Broadcast_Scalar(TBEGINC,          master_id)
    Call Broadcast_Scalar(Restart_In_Ver,   master_id)
  ENDIF

  !  *** DEACTIVATE ANY UNUSED OPTIONS
  IF( ISTRAN(2) < 1 ) ISTOPT(2) = 0

  ! *** SET TRANSPORT FLAG
  ISTRANACTIVE=0
  DO I=1,8
    IF( ISTRAN(I) > 0 ) ISTRANACTIVE = 1
  ENDDO

  !C7**  READ TIME-RELATED INTEGER PARAMETERS
  NCARD='7'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C7')
    READ(1,*,IOSTAT=ISO) NTC, NTSPTC, NLTC, NTTC, NTCPP, NTSTBC, NTCNB, NTCVB, NTSMMT, NFLTMT, NDRYSTP, NRAMPUP, NUPSTEP

    WRITE(7,1002)NCARD
    WRITE(7,*) NTC, NTSPTC, NLTC, NTTC, NTCPP, NTSTBC, NTCNB, NTCVB, NTSMMT, NFLTMT, NDRYSTP, NRAMPUP, NUPSTEP
    ! *** Not Used:  NTCPP
    ! *** Not Used:  NTCNB
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(NTC    , master_id)
  Call Broadcast_Scalar(NTSPTC , master_id)
  Call Broadcast_Scalar(NLTC   , master_id)
  Call Broadcast_Scalar(NTTC   , master_id)
  Call Broadcast_Scalar(NTCPP  , master_id)
  Call Broadcast_Scalar(NTSTBC , master_id)
  Call Broadcast_Scalar(NTCNB  , master_id)
  Call Broadcast_Scalar(NTCVB  , master_id)
  Call Broadcast_Scalar(NTSMMT , master_id)
  Call Broadcast_Scalar(NFLTMT , master_id)
  Call Broadcast_Scalar(NDRYSTP, master_id)
  Call Broadcast_Scalar(NRAMPUP, master_id)
  Call Broadcast_Scalar(NUPSTEP, master_id)

  IF( NRAMPUP < 1 ) NRAMPUP = 1
  IF( NUPSTEP < 2 ) NUPSTEP = 2

  NDRYSTP = ABS(NDRYSTP)   ! *** IF > 0, EFDC+ WILL WASTE ISOLATED CELL WATER AFTER NDRYSTP STEPS

  !C8**  READ TIME-RELATED REAL PARAMETERS
  NCARD='8'
  ! *** *******************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C8')
    READ(1,*,IOSTAT=ISO) TCON, TBEGIN, TIDALP, CF, ISCORV, ISDCCA, ISCFL, ISCFLM, DTSSFAC, DTSSDHDT, DTMAX

    WRITE(7,1002)NCARD
    WRITE(7,*) TCON, TBEGIN, TIDALP, CF, ISCORV, ISDCCA, ISCFL, ISCFLM, DTSSFAC, DTSSDHDT, DTMAX
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(TCON    , master_id)
  Call Broadcast_Scalar(TBEGIN  , master_id)
  Call Broadcast_Scalar(TIDALP  , master_id)
  Call Broadcast_Scalar(CF      , master_id)
  Call Broadcast_Scalar(ISCORV  , master_id)
  Call Broadcast_Scalar(ISDCCA  , master_id)
  Call Broadcast_Scalar(ISCFL   , master_id)
  Call Broadcast_Scalar(ISCFLM  , master_id)
  Call Broadcast_Scalar(DTSSFAC , master_id)
  Call Broadcast_Scalar(DTSSDHDT, master_id)
  Call Broadcast_Scalar(DTMAX   , master_id)

  IF( ISRESTI == 1 .AND. ICONTINUE == 1 )THEN
    ! ** FOR CONTINUATION MODE (TBEGINC IS SET FROM THE RESTART FILE)
    NTC = NTC + INT(TBEGIN - TBEGINC)
    TBEGIN = TBEGINC
  ENDIF

  IF( DTSSFAC > 0.0 )THEN
    ISDYNSTP = 1
    DT = TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
    IF( DTMAX <= DT*2. ) DTMAX = 3600.
  ELSE
    ISDYNSTP = 0
  ENDIF
  IF( IS2TIM == 0 ) ISDYNSTP = 0
  
  ! *** Calculating number of WASP outputs and snapshot time
  IF( RESSTEP > 0 ) THEN
    NWASPOUT = NTC*TIDALP/RESSTEP + 1
    CALL AllocateDSI(WASPTIME, NWASPOUT, 0.)
    DO IT = 1,NWASPOUT
      WASPTIME(IT) = TBEGIN + IT*RESSTEP/DBLE(86400.)
    ENDDO
  ENDIF
  !C9**  READ SPACE RELATED AND SMOOTHING PARAMETERS
  NCARD='9'
  ! *** *******************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C9')
    READ(1,*,IOSTAT=ISO) IC_GLOBAL, JC_GLOBAL, LC_GLOBAL, LVC, ISCLO, NDM, LDM, ISMASK, NBLOCKED, ISCONNECT, NSHMAX, NSBMAX, WSMH, WSMB
    WRITE(7,1002)NCARD
    WRITE(7,*) IC_GLOBAL, JC_GLOBAL, LC_GLOBAL, LVC, ISCLO, NDM, LDM, ISMASK, NBLOCKED, ISCONNECT, NSHMAX, NSBMAX, WSMH, WSMB
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(IC_GLOBAL, master_id)
  Call Broadcast_Scalar(JC_GLOBAL, master_id)
  Call Broadcast_Scalar(LC_GLOBAL, master_id)

  Call Broadcast_Scalar(LVC      , master_id)
  Call Broadcast_Scalar(ISCLO    , master_id)
  Call Broadcast_Scalar(NDM      , master_id)
  Call Broadcast_Scalar(LDM      , master_id)
  Call Broadcast_Scalar(ISMASK   , master_id)
  Call Broadcast_Scalar(NBLOCKED,  master_id)
  Call Broadcast_Scalar(ISCONNECT, master_id)
  Call Broadcast_Scalar(NSHMAX   , master_id)
  Call Broadcast_Scalar(NSBMAX   , master_id)
  Call Broadcast_Scalar(WSMH     , master_id)
  Call Broadcast_Scalar(WSMB     , master_id)
  
  LA_Global = LC_Global - 1
  
  IS2LMC=0
  IF( KC<0 )THEN
    KC=-KC
    IS2LMC=1
  ENDIF

  ! *** DOMAIN DECOMPOSITION CHECKS FOR HORIZONTAL LOOPS
  ! *** MAKE CONSISTENT WITH NEW OMP APPROACH
  ! @todo remove this NDM approach
  LCM2T = LC_GLOBAL - 2

  LCM2T = LC - 2
  NDM = NTHREADS
  LDM = INT(FLOAT(LCM2T)/FLOAT(NTHREADS)) + 1

  IF( KC >= 2) ISITB = 0

  ! *** Remapping of IC and JC for domain decomposition with MPI
#ifdef _MPI
  ! *** Get the local x,y-ids for identifying where we are on 'node map'
  !   the + 1 is for fortran indexing as the id's will start with 0 otherwise
  x_id = domain_coords(1) + 1
  y_id = domain_coords(2) + 1
  
  ! *** Set IC and JC for the the subdomain
  IC = IC_DECOMP(x_id)
  IF( n_x_partitions > 1 )THEN
    IF( process_map(x_id-1,y_id) /= -1 )  IC = IC + n_ghost_rows   ! *** West Edge
    IF( process_map(x_id+1,y_id) /= -1 )  IC = IC + n_ghost_rows   ! *** East Edge
  ENDIF
  
  JC = JC_DECOMP(y_id)
  IF( n_y_partitions > 1 )THEN
    IF( process_map(x_id,y_id-1) /= -1 )  JC = JC + n_ghost_rows   ! *** South Edge
    IF( process_map(x_id,y_id+1) /= -1 )  JC = JC + n_ghost_rows   ! *** North Edge
  ENDIF
  Call MPI_Barrier(comm_2d, ierr)

  Call Parent_Grid

  global_max_width_x = ic_global + 1 + 4 
  global_max_width_y = jc_global + 1 + 4 

  ! *** After this point IC and JC are Local!

  ! *** write out local IC and JCs
  Call WriteBreak(mpi_log_unit)
  write(mpi_log_unit,'(A)') 'After Parent Grid Mapping in INPUT routine '
  write(mpi_log_unit,'(A,I5)') 'IC: ',IC
  write(mpi_log_unit,'(A,I5)') 'JC: ',JC
  write(mpi_log_unit,'(A,I5)') 'Global Max Width x-direction + 4 for ghost: ', global_max_width_x
  write(mpi_log_unit,'(A,I5)') 'Global Max Width x-direction + 4 for ghost: ', global_max_width_y
  Call WriteBreak(mpi_log_unit)

#else ! No MPI
  IC = IC_GLOBAL
  JC = JC_GLOBAL
  DO I = 1, IC_GLOBAL
    DO J = 1, JC_GLOBAL
      IG2IL(I) = I
      JG2JL(J) = J
    END DO
  END DO
#endif

  !C9A**  READ VERTICAL SPACE RELATED  PARAMETERS
  NCARD='9A'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C9A')
    READ(1,*,IOSTAT=ISO) KC, KSIG, ISETGVC, SELVREF, BELVREF, ISGVCCK
    
    WRITE(7,1002)NCARD
    WRITE(7,*) KC, KSIG, ISETGVC, SELVREF, BELVREF, ISGVCCK
    ! *** Not Used: KSIG
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(KC     , master_id)
  Call Broadcast_Scalar(KSIG   , master_id)
  Call Broadcast_Scalar(ISETGVC, master_id)
  Call Broadcast_Scalar(SELVREF, master_id)
  Call Broadcast_Scalar(BELVREF, master_id)
  Call Broadcast_Scalar(ISGVCCK, master_id)

  !C10*  READ LAYER THICKNESS IN VERTICAL
  NCARD='10'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C10')
    DO K=1,KC
      READ(1,*,IOSTAT=ISO) KDUM, DZCK(K)
      WRITE(7,1002)NCARD
      WRITE(7,*)KDUM,DZCK(K)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  endif

  Call Broadcast_Array(DZCK, master_id)

  !C11*  READ GRID, ROUGHNESS, MASKING AND DEPTH PARAMETERS
  NCARD='11'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C11')
    READ(1,*,IOSTAT=ISO) DX, DY, DXYCVT, IMDXDY, ZBRADJ, ZBRCVRT, HMIN, HADADJ, HCVRT, HDRY, HWET, BELADJ, BELCVRT

    WRITE(7,1002)NCARD
    WRITE(7,*) DX, DY, DXYCVT, IMDXDY, ZBRADJ, ZBRCVRT, HMIN, HADADJ, HCVRT, HDRY, HWET, BELADJ, BELCVRT
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(DX     , master_id)
  Call Broadcast_Scalar(DY     , master_id)
  Call Broadcast_Scalar(DXYCVT , master_id)
  Call Broadcast_Scalar(IMDXDY , master_id)
  Call Broadcast_Scalar(ZBRADJ , master_id)
  Call Broadcast_Scalar(ZBRCVRT, master_id)
  Call Broadcast_Scalar(HMIN   , master_id)
  Call Broadcast_Scalar(HADADJ , master_id)
  Call Broadcast_Scalar(HCVRT  , master_id)
  Call Broadcast_Scalar(HDRY   , master_id)
  Call Broadcast_Scalar(HWET   , master_id)
  Call Broadcast_Scalar(BELADJ , master_id)
  Call Broadcast_Scalar(BELCVRT, master_id)

  HDRYICE = 0.91*HDRY
  HDRYWAV = 1.2*HWET

  !C11A* READ TWO-LAYER MOMENTUM FLUX AND CURVATURE ACCELERATION
  !     CORRECTION FACTORS
  NCARD='11A'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C11A')
    READ(1,*,IOSTAT=ISO) ICK2COR, CK2UUM, CK2VVM, CK2UVM, CK2UUC, CK2VVC, CK2UVC, CK2FCX, CK2FCY
    WRITE(7,1002)NCARD

    WRITE(7,*) ICK2COR, CK2UUM, CK2VVM, CK2UVM, CK2UUC, CK2VVC, CK2UVC, CK2FCX, CK2FCY
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ICK2COR, master_id)
  Call Broadcast_Scalar(CK2UUM , master_id)
  Call Broadcast_Scalar(CK2VVM , master_id)
  Call Broadcast_Scalar(CK2UVM , master_id)
  Call Broadcast_Scalar(CK2UUC , master_id)
  Call Broadcast_Scalar(CK2VVC , master_id)
  Call Broadcast_Scalar(CK2UVC , master_id)
  Call Broadcast_Scalar(CK2FCX , master_id)
  Call Broadcast_Scalar(CK2FCY , master_id)

  IF( ICK2COR >= 1 )THEN
    IS2LMC=ICK2COR
  END IF

  !C11B* READ CORNER CELL BOTTOM STRESS CORRECTION OPTIONS
  NCARD='11B'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C11B')
    READ(1,*,IOSTAT=ISO) ISCORTBC, ISCORTBCD, FSCORTBC
    
    WRITE(7,1002)NCARD
    WRITE(7,*) ISCORTBC,ISCORTBCD,FSCORTBC
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ISCORTBC, master_id)
  Call Broadcast_Scalar(ISCORTBCD,master_id)
  Call Broadcast_Scalar(FSCORTBC, master_id)
  
  !C12*  READ TURBULENT DIFFUSION PARAMETERS
  NCARD='12'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C12')
    READ(1,*,IOSTAT=ISO) AHO, AHD, AVO, ABO, AVMX, ABMX, VISMUD, AVCON, ZBRWALL
  
    WRITE(7,1002)NCARD
    WRITE(7,*) AHO, AHD, AVO, ABO, AVMX, ABMX, VISMUD, AVCON, ZBRWALL
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(AHO    , master_id)
  Call Broadcast_Scalar(AHD    , master_id)
  Call Broadcast_Scalar(AVO    , master_id)
  Call Broadcast_Scalar(ABO    , master_id)
  Call Broadcast_Scalar(AVMX   , master_id)
  Call Broadcast_Scalar(ABMX   , master_id)
  Call Broadcast_Scalar(VISMUD , master_id)
  Call Broadcast_Scalar(AVCON  , master_id)
  Call Broadcast_Scalar(ZBRWALL, master_id)

  ISAVCOMP = 1
  IF( AVCON == 0. ) ISAVCOMP = 0

  !C12A*  READ TURBULENCE CLOSURE OPTIONS
  NCARD='12A'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C12A')
    READ(1,*,IOSTAT=ISO) ISTOPT(0), ISSQL, ISAVBMX, ISFAVB, ISINWV, ISLLIM, IFPROX, XYRATIO, BC_EDGEFACTOR

    WRITE(7,1002)NCARD
    WRITE(7,*) ISTOPT(0), ISSQL, ISAVBMX, ISFAVB, ISINWV, IFPROX, XYRATIO, BC_EDGEFACTOR
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ISTOPT(0)    , master_id)
  Call Broadcast_Scalar(ISSQL        , master_id)
  Call Broadcast_Scalar(ISAVBMX      , master_id)
  Call Broadcast_Scalar(ISFAVB       , master_id)
  Call Broadcast_Scalar(ISINWV       , master_id)
  Call Broadcast_Scalar(IFPROX       , master_id)
  Call Broadcast_Scalar(XYRATIO      , master_id)
  Call Broadcast_Scalar(BC_EDGEFACTOR, master_id)

  IF( BC_EDGEFACTOR < 0 ) BC_EDGEFACTOR = 0.0
  IF( BC_EDGEFACTOR > 1 ) BC_EDGEFACTOR = 1.0

  !C13*  READ TURBULENCE CLOSURE PARAMETERS
  NCARD='13'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C13')
    ! *** PMC - CTE2 NOT USED
    READ(1,*,IOSTAT=ISO) VKC, CTURB, CTURB2B, CTE1, CTE2, CTE3, CTE4, CTE5, RIQMAX, QQMIN, QQLMIN, DMLMIN

    WRITE(7,1002)NCARD
    WRITE(7,*) VKC, CTURB, CTURB2B, CTE1, CTE2, CTE3, CTE4, CTE5, RIQMAX, QQMIN, QQLMIN, DMLMIN
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(VKC      , master_id)
  Call Broadcast_Scalar(CTURB    , master_id)
  Call Broadcast_Scalar(CTURB2B  , master_id)
  Call Broadcast_Scalar(CTE1     , master_id)
  Call Broadcast_Scalar(CTE2     , master_id)
  Call Broadcast_Scalar(CTE3     , master_id)
  Call Broadcast_Scalar(CTE4     , master_id)
  Call Broadcast_Scalar(CTE5     , master_id)
  Call Broadcast_Scalar(RIQMAX   , master_id)
  Call Broadcast_Scalar(QQMIN    , master_id)
  Call Broadcast_Scalar(QQLMIN   , master_id)
  Call Broadcast_Scalar(DMLMIN   , master_id)

  !C14*  READ TIDAL & ATMOSPHERIC FORCING, GROUND WATER AND SUBGRID CHANNEL PARAMETERS
  NCARD='14'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C14')
    READ(1,*,IOSTAT=ISO) MTIDE, NWSER, NASER, ISGWIT, ISCHAN, ISWAVE, ITIDASM, ISPERC, ISBODYF, ISPNHYDS, ISPROPWASH

    WRITE(7,1002)NCARD
    WRITE(7,*) MTIDE, NWSER, NASER, ISGWIT, ISCHAN, ISWAVE, ITIDASM, ISPERC, ISBODYF, ISPNHYDS, ISPROPWASH
    IF( ISO > 0 ) GOTO 100
  end if

  Call Broadcast_Scalar(MTIDE     , master_id)
  Call Broadcast_Scalar(NWSER     , master_id)
  Call Broadcast_Scalar(NASER     , master_id)
  Call Broadcast_Scalar(ISGWIT    , master_id)
  Call Broadcast_Scalar(ISCHAN    , master_id)
  Call Broadcast_Scalar(ISWAVE    , master_id)
  Call Broadcast_Scalar(ITIDASM   , master_id)
  Call Broadcast_Scalar(ISPERC    , master_id)
  Call Broadcast_Scalar(ISBODYF   , master_id)
  Call Broadcast_Scalar(ISPNHYDS  , master_id)
  Call Broadcast_Scalar(ISPROPWASH, master_id)

  ISWCBL = 0
  ISWVSD = 0
  IF( ISPERC > 0 ) ISGWIT = 3
  
  ! *** INITIALIZE PROPWASH CALCULATIONS FLAG
  IF( ISPROPWASH > 0 ) propwash_on = .TRUE.

  !C14A* READ SAND GRAIN NIKURADSE ROUGHNESS
  KSW = 0.00001 * 2.5  ! *** DEFAULT IS 10 MICRON

  IF( ISWAVE >= 1 )THEN
    NCARD='14A'
    if( process_id == master_id )THEN
      CALL SEEK('C14A')
      READ(1,*,IOSTAT=ISO) KSW, IUSEWVCELLS, IFWAVE, SWANGRP, ISSTEAD
      IF( ISO > 0 ) GOTO 100
    endif

    Call Broadcast_Scalar(KSW        , master_id)
    Call Broadcast_Scalar(IUSEWVCELLS, master_id)
    Call Broadcast_Scalar(IFWAVE     , master_id)
    Call Broadcast_Scalar(SWANGRP    , master_id)
    Call Broadcast_Scalar(ISSTEAD    , master_id)

    IF( ISWAVE == 1 .OR. ISWAVE == 2 .OR. ISWAVE == 4 )THEN
      NCARD='14B'
      if( process_id == master_id )THEN
        CALL SEEK('C14B')
        READ(1,*,IOSTAT=ISO) ISWRSR, ISWRSI, WVDISV, WVLSH, WVLSX, ISWVSD, WVLCAL, NTSWV, ISWCBL, ISDZBR
        IF( ISO > 0 ) GOTO 100
      endif

      Call Broadcast_Scalar(ISWRSR, master_id)
      Call Broadcast_Scalar(ISWRSI, master_id)
      Call Broadcast_Scalar(WVDISV, master_id)
      Call Broadcast_Scalar(WVLSH , master_id)
      Call Broadcast_Scalar(WVLSX , master_id)
      Call Broadcast_Scalar(ISWVSD, master_id)
      Call Broadcast_Scalar(WVLCAL, master_id)
      Call Broadcast_Scalar(NTSWV , master_id)
      Call Broadcast_Scalar(ISWCBL, master_id)
      Call Broadcast_Scalar(ISDZBR, master_id)

      NTSWV = MAX(NTSWV,1)
    ENDIF
  ENDIF

  !C14C TIME & SPACE VARYING FORCINGS     2018-10-12, NTL:
  NCARD='14C'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C14C')
    READ(1,*,IOSTAT=ISO) BATHY.IFLAG, ROUGH.IFLAG, VEGE.IFLAG,  GWSP.IFLAG,  WIND.IFLAG, PRESS.IFLAG,  &
                         RAIN.IFLAG,  EVAP.IFLAG,  SHELT.IFLAG, SHADE.IFLAG, SNOW.IFLAG, ICETHK.IFLAG, &
                         SEDZLJER.IFLAG

    WRITE(7,1002)NCARD
    WRITE(7,*) BATHY.IFLAG, ROUGH.IFLAG, VEGE.IFLAG,  GWSP.IFLAG,  WIND.IFLAG, PRESS.IFLAG,  &
               RAIN.IFLAG,  EVAP.IFLAG,  SHELT.IFLAG, SHADE.IFLAG, SNOW.IFLAG, ICETHK.IFLAG, &
               SEDZLJER.IFLAG
  endif

  Call Broadcast_Scalar(BATHY.IFLAG  , master_id)
  Call Broadcast_Scalar(ROUGH.IFLAG , master_id)
  Call Broadcast_Scalar(VEGE.IFLAG  , master_id)
  Call Broadcast_Scalar(GWSP.IFLAG  , master_id)
  Call Broadcast_Scalar(WIND.IFLAG  , master_id)
  Call Broadcast_Scalar(PRESS.IFLAG , master_id)
  Call Broadcast_Scalar(RAIN.IFLAG  , master_id)
  Call Broadcast_Scalar(EVAP.IFLAG  , master_id)
  Call Broadcast_Scalar(SHELT.IFLAG , master_id)
  Call Broadcast_Scalar(SHADE.IFLAG , master_id)
  Call Broadcast_Scalar(SNOW.IFLAG , master_id)
  Call Broadcast_Scalar(ICETHK.IFLAG, master_id)
  Call Broadcast_Scalar(SEDZLJER.IFLAG, master_id)

  !C14C  PARAMETRIC CYCLONE TIME & SPACE VARYING FORCING     2021-04-23, NTL:
  NCARD='14D'
  ! *** ********************************************************
  IF( process_id == master_id )THEN
    CALL SEEK('C14D')
    READ(1,*,IOSTAT=ISO) ICYCLONE,RHO_A,PINF,THETAMAX,WDRAG1,WDRAG2,CDRAG1,CDRAG2

    WRITE(7,1002)NCARD
    WRITE(7,*) ICYCLONE,RHO_A,PINF,THETAMAX,WDRAG1,WDRAG2,CDRAG1,CDRAG2
  ENDIF
  
  Call Broadcast_Scalar(ICYCLONE, master_id)
  Call Broadcast_Scalar(RHO_A,    master_id)
  Call Broadcast_Scalar(PINF,     master_id)
  Call Broadcast_Scalar(THETAMAX, master_id)
  Call Broadcast_Scalar(WDRAG1,   master_id)
  Call Broadcast_Scalar(WDRAG2,   master_id)
  Call Broadcast_Scalar(CDRAG1,   master_id)
  Call Broadcast_Scalar(CDRAG2,   master_id)
  
  IF( MTIDE > 0 )THEN
    !C15*  READ PERIODIC FORCING (TIDAL) CONSTITUENT SYMBOLS AND PERIODS
    NCARD='15'
    if( process_id == master_id )THEN
      CALL SEEK('C15')
      DO M=1,MTIDE
        READ(1,*,IOSTAT=ISO) SYMBOL(M),TCP(M)

        WRITE(7,1002)NCARD
        WRITE(7,*) SYMBOL(M),TCP(M)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    !all Broadcast_Array(SYMBOL, master_id)    ! DELME - TODO
    Call Broadcast_Array(TCP, master_id)

  ENDIF

  !C16*  READ SURFACE ELEVATION OR PRESSURE BOUNDARY CONDITION PARAMETERS
  NCARD='16'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C16')
    READ(1,*,IOSTAT=ISO) NPBS, NPBW, NPBE, NPBN, NPFOR_Readin, NPFORT, NPSER, PDGINIT

    WRITE(7,1002)NCARD
    WRITE(7,*) NPBS, NPBW, NPBE, NPBN, NPFOR_Readin, NPFORT, NPSER, PDGINIT
    IF( ISO > 0 ) GOTO 100

  endif

  Call Broadcast_Scalar(NPBS   , master_id)
  Call Broadcast_Scalar(NPBW   , master_id)
  Call Broadcast_Scalar(NPBE   , master_id)
  Call Broadcast_Scalar(NPBN   , master_id)
  Call Broadcast_Scalar(NPFOR_Readin , master_id)
  Call Broadcast_Scalar(NPFORT , master_id)
  Call Broadcast_Scalar(NPSER  , master_id)
  Call Broadcast_Scalar(PDGINIT, master_id)

  IF( NPFOR_Readin > 0 )THEN
    !C17*  READ PERIODIC FORCING (TIDAL) SURFACE ELEVATION OR
    NCARD='17'
    if( process_id == master_id )THEN
      CALL SEEK('C17')
      DO NP=1,NPFOR_Readin
        DO M=1,MTIDE
          IF( NPFORT == 0 )THEN
            READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
            WRITE(7,1002)NCARD
            WRITE(7,*) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M)
            IF( ISO > 0 ) GOTO 100
          ELSEIF( NPFORT >= 1 )THEN
            READ(1,*,IOSTAT=ISO) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M)
            RAD=PI2*PFPH(NP,M)/TCP(M)
            CPFAM0(NP,M)=PFAM(NP,M)*COS(RAD)
            SPFAM0(NP,M)=PFAM(NP,M)*SIN(RAD)
            WRITE(7,1002)NCARD
            WRITE(7,*) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M), CPFAM0(NP,M), SPFAM0(NP,M)
            IF( ISO > 0 ) GOTO 100
            READ(1,*,IOSTAT=ISO) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M)
            RAD=PI2*PFPH(NP,M)/TCP(M)
            CPFAM1(NP,M)=PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
            SPFAM1(NP,M)=PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
            WRITE(7,1002)NCARD
            WRITE(7,*) NDUM, CDUM, PFAM(NP,M), PFPH(NP,M), CPFAM1(NP,M), SPFAM1(NP,M)
            CPFAM2(NP,M)=0.0
            SPFAM2(NP,M)=0.0
          ELSEIF( NPFORT == 2 )THEN
            READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),PFX2(NP,M)
            RAD=PI2*PFPH(NP,M)/TCP(M)
            IF( PFX2(NP,M)>0.0 )THEN
              CPFAM2(NP,M)=PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
              SPFAM2(NP,M)=PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
            ELSE
              CPFAM2(NP,M)=0.
              SPFAM2(NP,M)=0.
            ENDIF
          ENDIF

        ENDDO
      ENDDO

    endif

    IF( NPFORT == 0 )THEN
      Call Broadcast_Array(PFAM,   master_id)
      Call Broadcast_Array(PFPH,   master_id)
    ELSEIF( NPFORT >= 1 )THEN
      Call Broadcast_Array(PFAM,   master_id)
      Call Broadcast_Array(PFPH,   master_id)
      Call Broadcast_Array(CPFAM0, master_id)
      Call Broadcast_Array(SPFAM0, master_id)
      Call Broadcast_Array(CPFAM1, master_id)
      Call Broadcast_Array(SPFAM1, master_id)
      Call Broadcast_Array(CPFAM2, master_id)
      Call Broadcast_Array(SPFAM2, master_id)
    ELSEIF( NPFORT == 2 )THEN
      Call Broadcast_Array(PFAM,   master_id)
      Call Broadcast_Array(PFPH,   master_id)
      Call Broadcast_Array(PFX2,   master_id)
      Call Broadcast_Array(CPFAM2, master_id)
      Call Broadcast_Array(SPFAM2, master_id)
    ENDIF

  ENDIF

  IF( NPBS > 0 )THEN
    !C18*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
    !     ON SOUTH OPEN BOUNDARIES
    NCARD='18'
    if( process_id == master_id )THEN
      CALL SEEK('C18')
      IF( NPFORT == 0 )THEN
        DO L=1,NPBS
          ! *** Added global arrays for MPI
          READ(1,*,IOSTAT=ISO) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L)
          IF( ISO > 0 ) GOTO 100

          DO M=1,MTIDE
            IF( NPFORS == 0) EXIT
            RAD=PI2*PFPH(NPFORS,M)/TCP(M)
            AMP=G*PFAM(NPFORS,M)
            PCBS_GL(L,M)=AMP*COS(RAD)
            PSBS_GL(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 1 )THEN
        DO L=1,NPBS
          READ(1,*,IOSTAT=ISO) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L), NPSERS1_GL(L), TPCOORDS_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L), NPSERS1_GL(L), TPCOORDS_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORS == 0) EXIT
            PCBS_GL(L,M)=CPFAM0(NPFORS,M)+TPCOORDS_GL(L)*CPFAM1(NPFORS,M) + TPCOORDS_GL(L)*TPCOORDS_GL(L)*CPFAM2(NPFORS,M)
            PSBS_GL(L,M)=SPFAM0(NPFORS,M)+TPCOORDS_GL(L)*SPFAM1(NPFORS,M) + TPCOORDS_GL(L)*TPCOORDS_GL(L)*SPFAM2(NPFORS,M)
            TMPAMP=SQRT(PCBS(L,M)*PCBS(L,M)+PSBS_GL(L,M)*PSBS_GL(L,M))
            TMPPHS=ATAN2(PSBS_GL(L,M),PCBS_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBS_GL(L,M)=G*PCBS_GL(L,M)
            PSBS_GL(L,M)=G*PSBS_GL(L,M)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 2 )THEN
        DO L=1,NPBS
          READ(1,*,IOSTAT=ISO) IPBS_GL(L),JPBS_GL(L),ISPBS_GL(L), ISPRS_GL(L), NPFORS,NPSERS_GL(L),NPSERS1_GL(L), TPCOORDS_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBS_GL(L), JPBS_GL(L), ISPBS_GL(L), ISPRS_GL(L), NPFORS, NPSERS_GL(L), NPSERS1_GL(L), TPCOORDS_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORS == 0) EXIT
            BOTTOM=PFX2(NPFORS,M)*(1.0-PFX2(NPFORS,M))
            TOP1=TPCOORDS_GL(L)*PFX2(NPFORS,M)*(TPCOORDS_GL(L)-PFX2(NPFORS,M))
            TOP2=TPCOORDS_GL(L)*(1.0-TPCOORDS_GL(L))
            IF( BOTTOM == 0.0 )THEN
              TOP1=TPCOORDS_GL(L)
              TOP2=TPCOORDS_GL(L)*TPCOORDS_GL(L)
            ELSE
              TOP1=TOP1/BOTTOM
              TOP2=TOP2/BOTTOM
            ENDIF
            PCBS_GL(L,M)=CPFAM0(NPFORS,M)+TOP1*CPFAM1(NPFORS,M)+TOP2*CPFAM2(NPFORS,M)
            PSBS_GL(L,M)=SPFAM0(NPFORS,M)+TOP1*SPFAM1(NPFORS,M)+TOP2*SPFAM2(NPFORS,M)
            TMPAMP=SQRT(PCBS_GL(L,M)*PCBS_GL(L,M)+PSBS_GL(L,M)*PSBS_GL(L,M))
            TMPPHS=ATAN2(PSBS_GL(L,M),PCBS_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBS_GL(L,M)=G*PCBS_GL(L,M)
            PSBS_GL(L,M)=G*PSBS_GL(L,M)
          ENDDO
        ENDDO
      ENDIF
2068  FORMAT(I4,3X,A2,5X,E14.4,3E14.5,5X,2I5)
2069  FORMAT(I4,3X,A2,5X,2E14.4,5X,2I5)

    endif
    
    Call Broadcast_Array(IPBS_GL,     master_id)
    Call Broadcast_Array(JPBS_GL,     master_id)
    Call Broadcast_Array(ISPBS_GL,    master_id)
    Call Broadcast_Array(ISPRS_GL,    master_id)
    Call Broadcast_Array(NPSERS_GL,   master_id)
    Call Broadcast_Array(PCBS_GL,     master_id)
    Call Broadcast_Array(PSBS_GL,     master_id)
    IF( NPFORT > 0 )THEN
      Call Broadcast_Array(NPSERS1_GL,  master_id)
      Call Broadcast_Array(TPCOORDS_GL, master_id)
    ENDIF
  ENDIF

  IF( NPBW > 0 )THEN
    !C19*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
    !     ON WEST OPEN BOUNDARIES
    NCARD='19'
    if( process_id == master_id )THEN

      CALL SEEK('C19')
      IF( NPFORT == 0 )THEN
        DO L=1,NPBW
          READ(1,*,IOSTAT=ISO) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L)
          WRITE(7,1002)NCARD
          WRITE(7,*) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORW == 0) EXIT
            RAD=PI2*PFPH(NPFORW,M)/TCP(M)
            AMP=G*PFAM(NPFORW,M)
            PCBW_GL(L,M)=AMP*COS(RAD)
            PSBW_GL(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 1 )THEN
        DO L=1,NPBW
          READ(1,*,IOSTAT=ISO) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORW == 0) EXIT
            PCBW_GL(L,M)=CPFAM0(NPFORW,M)+TPCOORDW_GL(L)*CPFAM1(NPFORW,M)+TPCOORDW_GL(L)*TPCOORDW_GL(L)*CPFAM2(NPFORW,M)
            PSBW_GL(L,M)=SPFAM0(NPFORW,M)+TPCOORDW_GL(L)*SPFAM1(NPFORW,M)+TPCOORDW_GL(L)*TPCOORDW_GL(L)*SPFAM2(NPFORW,M)
            TMPAMP=SQRT(PCBW_GL(L,M)*PCBW_GL(L,M)+PSBW_GL(L,M)*PSBW_GL(L,M))
            TMPPHS=0.0
            IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBW_GL(L,M),PCBW_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBW_GL(L,M)=G*PCBW_GL(L,M)
            PSBW_GL(L,M)=G*PSBW_GL(L,M)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 2 )THEN
        DO L=1,NPBW
          READ(1,*,IOSTAT=ISO) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBW_GL(L), JPBW_GL(L), ISPBW_GL(L), ISPRW_GL(L), NPFORW, NPSERW_GL(L), NPSERW1_GL(L), TPCOORDW_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORW == 0) EXIT
            BOTTOM=PFX2(NPFORW,M)*(1.0-PFX2(NPFORW,M))
            TOP1=TPCOORDW_GL(L)*PFX2(NPFORW,M)*(TPCOORDW_GL(L)-PFX2(NPFORW,M))
            TOP2=TPCOORDW_GL(L)*(1.0-TPCOORDW_GL(L))
            IF( BOTTOM == 0.0 )THEN
              TOP1=TPCOORDW_GL(L)
              TOP2=TPCOORDW_GL(L)*TPCOORDW_GL(L)
            ELSE
              TOP1=TOP1/BOTTOM
              TOP2=TOP2/BOTTOM
            ENDIF
            PCBW_GL(L,M)=CPFAM0(NPFORW,M)+TOP1*CPFAM1(NPFORW,M)+TOP2*CPFAM2(NPFORW,M)
            PSBW_GL(L,M)=SPFAM0(NPFORW,M)+TOP1*SPFAM1(NPFORW,M)+TOP2*SPFAM2(NPFORW,M)
            TMPAMP=SQRT(PCBW_GL(L,M)*PCBW_GL(L,M)+PSBW_GL(L,M)*PSBW_GL(L,M))
            TMPPHS=0.0
            IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBW_GL(L,M),PCBW_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBW_GL(L,M)=G*PCBW_GL(L,M)
            PSBW_GL(L,M)=G*PSBW_GL(L,M)
          ENDDO
        ENDDO
      ENDIF
    endif  
    
    Call Broadcast_Array(IPBW_GL,     master_id)
    Call Broadcast_Array(JPBW_GL,     master_id)
    Call Broadcast_Array(ISPBW_GL,    master_id)
    Call Broadcast_Array(ISPRW_GL,    master_id)
    Call Broadcast_Array(NPSERW_GL,   master_id)
    Call Broadcast_Array(PCBW_GL,     master_id)
    Call Broadcast_Array(PSBW_GL,     master_id)
    IF( NPFORT > 0 )THEN
      Call Broadcast_Array(NPSERW1_GL,  master_id)
      Call Broadcast_Array(TPCOORDW_GL, master_id)
    ENDIF
  ENDIF

  IF( NPBE > 0 )THEN
    !C20*  READ PERIODIC FORCING (TIDAL)ELEVATION BOUNDARY CONDTIONS
    !     ON EAST OPEN BOUNDARIES
    NCARD='20'
    if( process_id == master_id )THEN
      CALL SEEK('C20')
      IF( NPFORT == 0 )THEN
        DO L=1,NPBE
          READ(1,*,IOSTAT=ISO) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORE == 0) EXIT
            RAD=PI2*PFPH(NPFORE,M)/TCP(M)
            AMP=G*PFAM(NPFORE,M)
            PCBE_GL(L,M)=AMP*COS(RAD)
            PSBE_GL(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 1 )THEN
        DO L=1,NPBE
          READ(1,*,IOSTAT=ISO) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORE == 0) EXIT
            PCBE_GL(L,M)=CPFAM0(NPFORE,M)+TPCOORDE_GL(L)*CPFAM1(NPFORE,M)+TPCOORDE_GL(L)*TPCOORDE_GL(L)*CPFAM2(NPFORE,M)
            PSBE_GL(L,M)=SPFAM0(NPFORE,M)+TPCOORDE_GL(L)*SPFAM1(NPFORE,M)+TPCOORDE_GL(L)*TPCOORDE_GL(L)*SPFAM2(NPFORE,M)
            TMPAMP=SQRT(PCBE_GL(L,M)*PCBE_GL(L,M)+PSBE_GL(L,M)*PSBE_GL(L,M))
            TMPPHS=0.0
            IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBE_GL(L,M),PCBE_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBE_GL(L,M)=G*PCBE_GL(L,M)
            PSBE_GL(L,M)=G*PSBE_GL(L,M)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 2 )THEN
        DO L=1,NPBE
          READ(1,*,IOSTAT=ISO) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBE_GL(L), JPBE_GL(L), ISPBE_GL(L), ISPRE_GL(L), NPFORE, NPSERE_GL(L), NPSERE1_GL(L), TPCOORDE_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORE == 0) EXIT
            BOTTOM=PFX2(NPFORE,M)*(1.0-PFX2(NPFORE,M))
            TOP1=TPCOORDE_GL(L)*PFX2(NPFORE,M)*(TPCOORDE_GL(L)-PFX2(NPFORE,M))
            TOP2=TPCOORDE_GL(L)*(1.0-TPCOORDE_GL(L))
            IF( BOTTOM == 0.0 )THEN
              TOP1=TPCOORDE_GL(L)
              TOP2=TPCOORDE_GL(L)*TPCOORDE_GL(L)
            ELSE
              TOP1=TOP1/BOTTOM
              TOP2=TOP2/BOTTOM
            ENDIF
            PCBE_GL(L,M)=CPFAM0(NPFORE,M)+TOP1*CPFAM1(NPFORE,M)+TOP2*CPFAM2(NPFORE,M)
            PSBE_GL(L,M)=SPFAM0(NPFORE,M)+TOP1*SPFAM1(NPFORE,M)+TOP2*SPFAM2(NPFORE,M)
            TMPAMP=SQRT(PCBE_GL(L,M)*PCBE_GL(L,M)+PSBE_GL(L,M)*PSBE_GL(L,M))
            TMPPHS=0.0
            IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBE_GL(L,M),PCBE_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBE_GL(L,M)=G*PCBE_GL(L,M)
            PSBE_GL(L,M)=G*PSBE_GL(L,M)
          ENDDO
        ENDDO
      ENDIF
    endif

    Call Broadcast_Array(IPBE_GL,     master_id)
    Call Broadcast_Array(JPBE_GL,     master_id)
    Call Broadcast_Array(ISPBE_GL,    master_id)
    Call Broadcast_Array(ISPRE_GL,    master_id)
    Call Broadcast_Array(NPSERE_GL,   master_id)
    Call Broadcast_Array(PCBE_GL,     master_id)
    Call Broadcast_Array(PSBE_GL,     master_id)
    IF( NPFORT > 0 )THEN
      Call Broadcast_Array(NPSERE1_GL,  master_id)
      Call Broadcast_Array(TPCOORDE_GL, master_id)
    ENDIF
  ENDIF

  IF( NPBN > 0 )THEN
    !C21*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS
    !     ON NORTH OPEN BOUNDARIES
    NCARD='21'
    if( process_id == master_id )THEN

      CALL SEEK('C21')
      IF( NPFORT == 0 )THEN
        DO L=1,NPBN
          READ(1,*,IOSTAT=ISO) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORN == 0) EXIT
            RAD=PI2*PFPH(NPFORN,M)/TCP(M)
            AMP=G*PFAM(NPFORN,M)
            PCBN_GL(L,M)=AMP*COS(RAD)
            PSBN_GL(L,M)=AMP*SIN(RAD)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT >= 1 )THEN
        DO L=1,NPBN
          READ(1,*,IOSTAT=ISO) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORN == 0) EXIT
            PCBN_GL(L,M)=CPFAM0(NPFORN,M)+TPCOORDN_GL(L)*CPFAM1(NPFORN,M)+TPCOORDN_GL(L)*TPCOORDN_GL(L)*CPFAM2(NPFORN,M)
            PSBN_GL(L,M)=SPFAM0(NPFORN,M)+TPCOORDN_GL(L)*SPFAM1(NPFORN,M)+TPCOORDN_GL(L)*TPCOORDN_GL(L)*SPFAM2(NPFORN,M)
            TMPAMP=SQRT(PCBN_GL(L,M)*PCBN_GL(L,M)+PSBN_GL(L,M)*PSBN_GL(L,M))
            TMPPHS=0.0
            IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBN_GL(L,M),PCBN_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBN_GL(L,M)=G*PCBN_GL(L,M)
            PSBN_GL(L,M)=G*PSBN_GL(L,M)
          ENDDO
        ENDDO
        
      ELSEIF( NPFORT == 2 )THEN
        DO L=1,NPBN
          READ(1,*,IOSTAT=ISO) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)

          WRITE(7,1002)NCARD
          WRITE(7,*) IPBN_GL(L), JPBN_GL(L), ISPBN_GL(L), ISPRN_GL(L), NPFORN, NPSERN_GL(L), NPSERN1_GL(L), TPCOORDN_GL(L)
          IF( ISO > 0 ) GOTO 100
          DO M=1,MTIDE
            IF( NPFORN == 0) EXIT
            BOTTOM=PFX2(NPFORN,M)*(1.0-PFX2(NPFORN,M))
            TOP1=TPCOORDN_GL(L)*PFX2(NPFORN,M)*(TPCOORDN_GL(L)-PFX2(NPFORN,M))
            TOP2=TPCOORDN_GL(L)*(1.0-TPCOORDN_GL(L))
            IF( BOTTOM == 0.0 )THEN
              TOP1=TPCOORDN_GL(L)
              TOP2=TPCOORDN_GL(L)*TPCOORDN_GL(L)
            ELSE
              TOP1=TOP1/BOTTOM
              TOP2=TOP2/BOTTOM
            ENDIF
            PCBN_GL(L,M)=CPFAM0(NPFORN,M)+TOP1*CPFAM1(NPFORN,M)+TOP2*CPFAM2(NPFORN,M)
            PSBN_GL(L,M)=SPFAM0(NPFORN,M)+TOP1*SPFAM1(NPFORN,M)+TOP2*SPFAM2(NPFORN,M)
            TMPAMP=SQRT(PCBN_GL(L,M)*PCBN_GL(L,M)+PSBN_GL(L,M)*PSBN_GL(L,M))
            TMPPHS=0.0
            IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBN_GL(L,M),PCBN_GL(L,M))
            TMPPHS=TMPPHS*TCP(M)/PI2
            IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
            PCBN_GL(L,M)=G*PCBN_GL(L,M)
            PSBN_GL(L,M)=G*PSBN_GL(L,M)
          ENDDO
        ENDDO
      ENDIF
    endif

    Call Broadcast_Array(IPBN_GL,     master_id)
    Call Broadcast_Array(JPBN_GL,     master_id)
    Call Broadcast_Array(ISPBN_GL,    master_id)
    Call Broadcast_Array(ISPRN_GL,    master_id)
    Call Broadcast_Array(NPSERN_GL,   master_id)
    Call Broadcast_Array(PCBN_GL,     master_id)
    Call Broadcast_Array(PSBN_GL,     master_id)
    IF( NPFORT > 0 )THEN
      Call Broadcast_Array(NPSERN1_GL,  master_id)
      Call Broadcast_Array(TPCOORDN_GL, master_id)
    ENDIF
  ENDIF

  !22*  READ NUM OF SEDIMENT AMD TOXICS AND NUM OF CONCENTRATION TIME SERIES
  NCARD='22'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C22')
    READ(1,*,IOSTAT=ISO) NDYE,NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),NCSER(4),NCSER(5),NCSER(6),NCSER(7),ISSBAL

    WRITE(7,1002)NCARD
    WRITE(7,*) NDYE,NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),NCSER(4),NCSER(5),NCSER(6),NCSER(7),ISSBAL
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(NDYE  , master_id)
  Call Broadcast_Scalar(NTOX  , master_id)
  Call Broadcast_Scalar(NSED  , master_id)
  Call Broadcast_Scalar(NSND  , master_id)
  Call Broadcast_Scalar(ISSBAL, master_id)
  Call Broadcast_Array(NCSER(1:7), master_id)

  ! *** REMOVE UNUSED SETTINGS TO ALLOW FOR SELECTIVE ALLOCATIONS
  IF(  ISTRAN(6) < 1)  NSED = 0
  IF(  ISTRAN(7) < 1 ) NSND = 0
  IF(  NSED == 0 .AND. NSND == 0 ) ISTRAN(5) = 0
  IF(  ISTRAN(5) < 1 ) NTOX = 0
  IF( NTOX == 0 ) NCSER(5) = 0
  IF( NSED == 0 ) NCSER(6) = 0
  IF( NSND == 0 ) NCSER(7) = 0

  IF( ISTRAN(6) == 0 .AND. ISTRAN(7) == 0 )THEN
    ISSBAL = 0  ! *** PMC SINGLE LINE
  ENDIF

  !!22B*  BIOLOGICAL MODEL SETTINGS
  !NCARD='22B'
  !! *** ********************************************************
  !if( process_id == master_id )THEN
  !  CALL SEEK('C22B')
  !  READ(1,*,IOSTAT=ISO) NSF,ISFFARM,NSFCELLS
  !
  !  WRITE(7,1002)NCARD
  !  WRITE(7,*) NSF,ISFFARM,NSFCELLS
  !  IF( ISO > 0 ) GOTO 100
  !endif
  !
  !Call Broadcast_Scalar(NSF     , master_id)
  !Call Broadcast_Scalar(ISFFARM , master_id)
  !Call Broadcast_Scalar(NSFCELLS, master_id)

  !C23*  READ VELOCITY, VOL SOUR/SINK, FLOW CONTROL, & WITHDRAW/RETURN DATA
  NCARD='23'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C23')
    READ(1,*,IOSTAT=ISO) NQSIJ, NQJPIJ, NQSER, NQCTL, NQCTLT, NHYDST, NQWR, NQWRSR, ISDIQ, NQCTLSER, NQCTRULES

    WRITE(7,1002)NCARD
    WRITE(7,*) NQSIJ, NQJPIJ, NQSER, NQCTL, NQCTLT, NHYDST, NQWR, NQWRSR, ISDIQ, NQCTLSER, NQCTRULES
    IF( ISO > 0 ) GOTO 100
  endif
  ! *** Broadcast in Scan_EFDC

  NGRPID = 0
  IF( NQSIJ > 0 )THEN
    !C24*  READ VOLUME SOURCE/SINK LOCATIONS, MAGNITUDES, & VOL & CONC SERIES
    NCARD='24'
    if( process_id == master_id )THEN
      CALL SEEK('C24')
      DO NS=1,NQSIJ
        READ(1,*,IOSTAT=ISO) IQS_GL(NS), JQS_GL(NS), QSSE_GL(NS), NQSMUL_GL(NS), NQSMF_GL(NS), NQSERQ_GL(NS),   NCSERQ_GL(NS,1), NCSERQ_GL(NS,2), NCSERQ_GL(NS,3), NCSERQ_GL(NS,4), NCSERQ_GL(NS,5),  &
                                                                                               NCSERQ_GL(NS,6), NCSERQ_GL(NS,7), QWIDTH_GL(NS),   QFACTOR_GL(NS),  GRPID_GL(NS)

        WRITE(7,1002)NCARD
        WRITE(7,*)           IQS_GL(NS), JQS_GL(NS), QSSE_GL(NS), NQSMUL_GL(NS), NQSMF_GL(NS), NQSERQ_GL(NS),   NCSERQ_GL(NS,1), NCSERQ_GL(NS,2), NCSERQ_GL(NS,3), NCSERQ_GL(NS,4), NCSERQ_GL(NS,5),  &
                                                                                               NCSERQ_GL(NS,6), NCSERQ_GL(NS,7), QWIDTH_GL(NS),   QFACTOR_GL(NS),  GRPID_GL(NS)
        IF( ISO > 0 ) GOTO 100
        
        ! *** MOVED TO MAP_RIVER
        !DO K=1,KC
        !  QSS(K,NS) = QSSE(NS)*DZCK(K)
        !ENDDO
        !
        !! *** GET NUMBER IF GROUPS
        !IF( GRPID(NS) > NGRPID ) NGRPID = GRPID(NS)
      ENDDO
    endif

    Call Broadcast_Array(IQS_GL, master_id)
    Call Broadcast_Array(JQS_GL, master_id)
    Call Broadcast_Array(QSSE_GL, master_id)
    Call Broadcast_Array(NQSMUL_GL, master_id)
    Call Broadcast_Array(NQSMF_GL  , master_id)
    Call Broadcast_Array(NQSERQ_GL , master_id)
    Call Broadcast_Array(NCSERQ_GL, master_id)
    Call Broadcast_Array(QWIDTH_GL , master_id)
    Call Broadcast_Array(QFACTOR_GL, master_id)
    Call Broadcast_Array(GRPID_GL, master_id)

    !C25*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS
    !     SAL,TEM,DYE,SFL,TOX(1 TO NOTX)
    NCARD='25'
    if( process_id == master_id )THEN
      CALL SEEK('C25')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NQSIJ
        READ(1,*,IOSTAT=ISO) (CQSE_GL(L,M),M=1,MMAX)

        WRITE(7,1002)NCARD
        WRITE(7,*) (CQSE_GL(L,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(CQSE_GL, master_id)

    !C26*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS
    !     SED(1 TO NSED),SND(1 TO NSND)
    NCARD='26'

    if( process_id == master_id )THEN
      CALL SEEK('C26')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NQSIJ
        READ(1,*,IOSTAT=ISO) (CQSE_GL(L,M),M=MMIN,MMAX)
        WRITE(7,1002)NCARD

        WRITE(7,*) (CQSE_GL(L,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
  ENDIF

  Call Broadcast_Array(CQSE_GL, master_id)

  IF( NQJPIJ > 0 )THEN
    !C27*  READ JET/PLUME SOURCE LOCATIONS AND PARAMETERS
    NCARD='27'
    if( process_id == master_id )THEN
      CALL SEEK('C27')
      DO L=1,NQJPIJ
        READ(1,*,IOSTAT=ISO) IDUM,ICALJP(L),IQJP_GL(L),JQJP_GL(L),KQJP_GL(L),NPORTJP_GL(L),XJETL_GL(L),&
                             YJETL_GL(L),ZJET_GL(L),PHJET_GL(L),THJET_GL(L),DJET_GL(L),CFRD_GL(L),DJPER_GL(L)

        WRITE(7,1002)NCARD
        WRITE(7,*)IDUM,ICALJP(L),IQJP_GL(L),JQJP_GL(L),KQJP_GL(L),NPORTJP_GL(L),XJETL_GL(L),&
                       YJETL_GL(L),ZJET_GL(L),PHJET_GL(L),THJET_GL(L),DJET_GL(L),CFRD_GL(L),DJPER_GL(L)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(ICALJP, master_id)
    Call Broadcast_Array(IQJP_GL  , master_id)
    Call Broadcast_Array(JQJP_GL  , master_id)
    Call Broadcast_Array(KQJP_GL  , master_id)
    Call Broadcast_Array(NPORTJP_GL, master_id)
    Call Broadcast_Array(XJETL_GL , master_id)
    Call Broadcast_Array(YJETL_GL , master_id)
    Call Broadcast_Array(ZJET_GL  , master_id)
    Call Broadcast_Array(PHJET_GL  , master_id)
    Call Broadcast_Array(THJET_GL  , master_id)
    Call Broadcast_Array(DJET_GL  , master_id)
    Call Broadcast_Array(CFRD_GL   , master_id)
    Call Broadcast_Array(DJPER_GL  , master_id)


    !C28*  READ JET/PLUME SOURCE LOCATIONS AND PARAMETERS
    NCARD='28'
    if( process_id == master_id )THEN

      CALL SEEK('C28')
      DO L=1,NQJPIJ
        READ(1,*,IOSTAT=ISO) IDUM,NJEL(L),NJPMX(L),ISENT(L),ISTJP(L),NUDJP(L),IOUTJP(L),NZPRJP(L),ISDJP(L),IUPCJP(L),JUPCJP(L),KUPCJP(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IDUM,NJEL(L),NJPMX(L),ISENT(L),ISTJP(L),NUDJP(L),IOUTJP(L),NZPRJP(L),ISDJP(L),IUPCJP(L),JUPCJP(L),KUPCJP(L)
        IF( ISO > 0 ) GOTO 100
      ENDDO

    endif

    Call Broadcast_Array(NJEL  , master_id)
    Call Broadcast_Array(NJPMX , master_id)
    Call Broadcast_Array(ISENT , master_id)
    Call Broadcast_Array(ISTJP , master_id)
    Call Broadcast_Array(NUDJP , master_id)
    Call Broadcast_Array(IOUTJP, master_id)
    Call Broadcast_Array(NZPRJP, master_id)
    Call Broadcast_Array(ISDJP , master_id)
    Call Broadcast_Array(IUPCJP, master_id)
    Call Broadcast_Array(JUPCJP, master_id)
    Call Broadcast_Array(KUPCJP, master_id)

    !C29*  READ ADDITIONAL JET/PLUME PARAMETERS
    NCARD='29'
    if( process_id == master_id )THEN
      CALL SEEK('C29')
      DO L=1,NQJPIJ
        READ(1,*,IOSTAT=ISO) IDUM,QQCJP(L),NQSERJP(L),NQWRSERJP(L),NCSERJP(L,1),NCSERJP(L,2),NCSERJP(L,3),NCSERJP(L,4),NCSERJP(L,5),NCSERJP(L,6),NCSERJP(L,7)

        WRITE(7,1002)NCARD
        WRITE(7,*) IDUM,QQCJP(L),NQSERJP(L),NQWRSERJP(L),NCSERJP(L,1),NCSERJP(L,2),NCSERJP(L,3),NCSERJP(L,4),NCSERJP(L,5),NCSERJP(L,6),NCSERJP(L,7)
        
        NUDJPC(L) = 1
        IF( ISO > 0 ) GOTO 100

        IF( ICALJP(L) == 2 )THEN
          QWRCJP(L)=QQCJP(L)
          QQCJP(L)=0.
        ELSE
          QWRCJP(L)=0.
        ENDIF
      ENDDO
    endif

    Call Broadcast_Array(QQCJP,     master_id)
    Call Broadcast_Array(NQSERJP,   master_id)
    Call Broadcast_Array(NQWRSERJP, master_id)
    Call Broadcast_Array(NCSERJP  , master_id)

    IF( NQJPIJ > 1 )THEN
      DO L=2,NQJPIJ
        NUDJP(L)=NUDJP(1)
      ENDDO
    ENDIF


    !C30*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT
    !     JET/PLUME SOURCES SAL,TEM,DYE,SFL,TOX(1 TO NOTX)
    NCARD='30'
    if( process_id == master_id )THEN
      CALL SEEK('C30')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NQJPIJ
        READ(1,*,IOSTAT=ISO) (CQSE(M),M=1,MMAX)
        WRITE(7,1002)NCARD
        
        WRITE(7,*) (CQSE(M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
        
        IF( ICALJP(L) == 1 )THEN
          DO MS=1,MMAX
            CWRCJP(L,MS)=0.0
            DO K=1,KC
              CQCJP(K,L,MS)=CQSE(MS)
            ENDDO
          ENDDO
        ELSE
          DO MS=1,MMAX
            CWRCJP(L,MS)=CQSE(MS)
            DO K=1,KC
              CQCJP(K,L,MS)=0.0
            ENDDO
          ENDDO
        ENDIF
      ENDDO

    endif

    Call Broadcast_Array(CQCJP, master_id)
    Call Broadcast_Array(CWRCJP, master_id)
    Call Broadcast_Array(CQSE, master_id)

    !C31*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT
    !     JET/PLUME SOURCES SED(1 TO NSED),SND(1 TO NSND)
    NCARD='31'
    if( process_id == master_id )THEN
      CALL SEEK('C31')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      IF( ISTRAN(8) > 0 )THEN
        MMAX = MMAX + NWQV
      ENDIF
      DO L=1,NQJPIJ
        READ(1,*,IOSTAT=ISO) (CQSE(M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CQSE(M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
        
        MS = MAX(MMIN-1, 0)
        IF( ICALJP(L) == 1 )THEN
          DO M = 1,NWQV
            IF( ISKINETICS(M) > 0 )THEN
              MS = MS + 1
              DO K=1,KC
                CQCJP(K,L,MS) = CQSE(M + MMIN - 1)        ! *** Only used constituents that will be simulated
              ENDDO
            ENDIF
          ENDDO
        ELSE
          DO MS=MMIN,MMAX
            CWRCJP(L,MS)=CQSE(MS)
            DO K=1,KC
              CQCJP(K,L,MS)=0.
            ENDDO
          ENDDO
        ENDIF
      ENDDO

    endif
    Call Broadcast_Array(CQCJP, master_id)
    Call Broadcast_Array(CWRCJP, master_id)
    Call Broadcast_Array(CQSE, master_id)

  ENDIF

  IF( NQCTL > 0 )THEN
    !C32*  READ SURF ELEV OR PRESS DEPENDENT FLOW CONTROL STRUCTURE INFO
    NCARD='32'

    if( process_id == master_id )THEN
      CALL SEEK('C32')
      DO NC=1,NQCTL
        READ(1,*,IOSTAT=ISO) HYD_STR_GL(NC).IQCTLU, HYD_STR_GL(NC).JQCTLU, HYD_STR_GL(NC).IQCTLD, HYD_STR_GL(NC).JQCTLD,  HYD_STR_GL(NC).NQCTYP, HYD_STR_GL(NC).NQCTLQ,  HYD_STR_GL(NC).NQCMUL,    &
                             HYD_STR_GL(NC).HQCTLU, HYD_STR_GL(NC).HQCTLD, HYD_STR_GL(NC).QCTLMU, HYD_STR_GL(NC).QCTLGRP, HYD_STR_GL(NC).BQCLCE, HYD_STR_GL(NC).NQCMINS, HYD_STR_GL(NC).HS_FACTOR, &
                             HYD_STR_GL(NC).MMASKS, HYD_STR_GL(NC).HS_TRANSITION

        WRITE(7,1002)NCARD
        WRITE(7,*) HYD_STR_GL(NC).IQCTLU, HYD_STR_GL(NC).JQCTLU, HYD_STR_GL(NC).IQCTLD, HYD_STR_GL(NC).JQCTLD,  HYD_STR_GL(NC).NQCTYP, HYD_STR_GL(NC).NQCTLQ,  HYD_STR_GL(NC).NQCMUL,    &
                   HYD_STR_GL(NC).HQCTLU, HYD_STR_GL(NC).HQCTLD, HYD_STR_GL(NC).QCTLMU, HYD_STR_GL(NC).QCTLGRP, HYD_STR_GL(NC).BQCLCE, HYD_STR_GL(NC).NQCMINS, HYD_STR_GL(NC).HS_FACTOR, &
                   HYD_STR_GL(NC).MMASKS, HYD_STR_GL(NC).HS_TRANSITION
        IF( ISO > 0 ) GOTO 100
        
        IF( HYD_STR_GL(NC).NQCTYP > 4 )THEN
          IF( HYD_STR_GL(NC).HS_FACTOR <= 0  )THEN
            PRINT *,' *** BAD HS_FACTOR: NC, HS_FACTOR = ',NC, HYD_STR_GL(NC).HS_FACTOR
            CALL STOPP('')
          ENDIF
        ENDIF

      ENDDO
    endif

    DO NC=1,NQCTL
      Call Broadcast_Scalar(HYD_STR_GL(NC).IQCTLU     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).JQCTLU     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).IQCTLD     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).JQCTLD     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).NQCTYP     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).NQCTLQ     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).NQCMUL     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).HQCTLU     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).HQCTLD     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).QCTLMU     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).QCTLGRP    , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).BQCLCE     , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).NQCMINS    , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).HS_FACTOR  , master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).HS_TRANSITION, master_id)
      Call Broadcast_Scalar(HYD_STR_GL(NC).MMASKS     , master_id)
    ENDDO
    
    DO NC=1,NQCTL
      DO K=1,KC
        QCTLTO(K,NC)=0.
        QCTLT(K,NC,1)=0.
        QCTLT(K,NC,2)=0.
      ENDDO
    ENDDO
    
    !C32A*  READ THE EQUATION PARAMETERS FOR EACH OF THE HYDRAULIC STRUCTURE EQUATIONS
    IF( NHYDST > 0 )THEN
      NCARD='32A'
      if( process_id == master_id )THEN
        CALL SEEK('C32A')
        DO L=1,NHYDST
          ! *** COMPUTE FLOWS USING HYDRAULIC STRUCTURE EQUATIONS
          READ(1,*,IOSTAT=ISO) NS, NX, HS_REVERSE(L), HS_XSTYPE(L), HS_WIDTH(L), HS_HEIGHT(L), HS_LENGTH(L), HS_MANN(L), HS_ANGLE(L), HS_USELEV(L), HS_DSELEV(L),  &
                                       HS_COEFF(L,1), HS_COEFF(L,2), HS_COEFF(L,3), HS_COEFF(L,4)

          CALL HYDSTRUCT_CHECK(NX,L)  ! *** CHECK STRUCTURE DEFINITIONS

          ! *** USE THE SAME INVERT ELEVATION FOR US/DS FOR SLUICE GATES AND WEIRS
          IF( NX > 5 )THEN
            HS_DSELEV(L) = HS_USELEV(L)
          ENDIF

          WRITE(7,1002)NCARD
          WRITE(7,*) NS, NX, HS_REVERSE(L), HS_XSTYPE(L), HS_WIDTH(L), HS_HEIGHT(L), HS_LENGTH(L), HS_MANN(L), HS_ANGLE(L), HS_USELEV(L), HS_DSELEV(L), &
                             HS_COEFF(L,1), HS_COEFF(L,2), HS_COEFF(L,3), HS_COEFF(L,4)
          HS_ANGLE(L) = HS_ANGLE(L)*PI/180.0D0
        ENDDO
      ENDIF

      Call Broadcast_Scalar(NS, master_id)
      Call Broadcast_Scalar(NX, master_id)
      Call Broadcast_Array(HS_REVERSE, master_id)
      Call Broadcast_Array(HS_XSTYPE , master_id)
      Call Broadcast_Array(HS_WIDTH  , master_id)
      Call Broadcast_Array(HS_HEIGHT , master_id)
      Call Broadcast_Array(HS_LENGTH , master_id)
      Call Broadcast_Array(HS_MANN   , master_id)
      Call Broadcast_Array(HS_ANGLE  , master_id)
      Call Broadcast_Array(HS_USELEV , master_id)
      Call Broadcast_Array(HS_DSELEV , master_id)
      Call Broadcast_Array(HS_COEFF  , master_id)
      Call Broadcast_Array(HS_COEFF  , master_id)
      Call Broadcast_Array(HS_COEFF  , master_id)
      Call Broadcast_Array(HS_COEFF  , master_id)
    endif

    !C32B*  READ THE CONTROL INFO FOR ALL HYDRAULIC STRUCTURES
    NCARD='32B'
    if( process_id == master_id )THEN
      IF( NQCTLSER > 0 .OR. NQCTRULES > 0 )THEN
        CALL SEEK('C32B')
        DO NC=1,NQCTL
          HSCTL_GL(NC).ITYPE = 0
          HSCTL_GL(NC).ID = 0
          HSCTL_GL(NC).SUBID = 0

          READ(1,*,IOSTAT=ISO) HSCTL_GL(NC).ITYPE, HSCTL_GL(NC).ID, ITMPU, JTMPU, ITMPD, JTMPD, &
                               HSCTL_GL(NC).CUR.STATE, HSCTL_GL(NC).CUR.HEIGHT, HSCTL_GL(NC).CUR.WIDTH, &
                               HSCTL_GL(NC).CUR.SILL,  NDUM,                    HSCTL_GL(NC).CUR.FLOW
          
          HSCTL_GL(NC).IREFUP = ITMPU
          HSCTL_GL(NC).JREFUP = JTMPU
          HSCTL_GL(NC).IREFDN = ITMPD
          HSCTL_GL(NC).JREFDN = JTMPD

          WRITE(7,1002)NCARD
          WRITE(7,*) HSCTL_GL(NC).ITYPE, HSCTL_GL(NC).ID, &
                     HSCTL_GL(NC).IREFUP, HSCTL_GL(NC).JREFUP, HSCTL_GL(NC).IREFDN, HSCTL_GL(NC).JREFDN, &
                     HSCTL_GL(NC).CUR.STATE, HSCTL_GL(NC).CUR.HEIGHT, HSCTL_GL(NC).CUR.WIDTH, &
                     HSCTL_GL(NC).CUR.SILL,  HSCTL_GL(NC).CUR.FLOW
        ENDDO
      ELSE
        DO NC=1,NQCTL
          HSCTL_GL(NC).ITYPE = 0
          HSCTL_GL(NC).ID = 0
                    
          HSCTL_GL(NC).IREFUP = 0
          HSCTL_GL(NC).JREFUP = 0
          HSCTL_GL(NC).IREFDN = 0
          HSCTL_GL(NC).JREFDN = 0
          
          HSCTL_GL(NC).CUR.STATE = 0.
          HSCTL_GL(NC).CUR.HEIGHT = 0.
          HSCTL_GL(NC).CUR.WIDTH = 0.
          HSCTL_GL(NC).CUR.SILL = 0.
          HSCTL_GL(NC).CUR.FLOW = 0.
        ENDDO
      ENDIF
    endif

    Do L=1, NQCTL
      Call Broadcast_Scalar(HSCTL_GL(L).ITYPE     , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).ID        , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).IREFUP    , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).JREFUP    , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).IREFDN    , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).JREFDN    , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).CUR.STATE , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).CUR.HEIGHT, master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).CUR.WIDTH , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).CUR.SILL  , master_id)
      !Call Broadcast_Scalar(HSCTL_GL(L).CUR.UNITS , master_id)
      Call Broadcast_Scalar(HSCTL_GL(L).CUR.FLOW  , master_id)
    End do

  ENDIF

  IF( NQWR > 0 )THEN
    !C33*  READ FLOW WITHDRAWAL, HEAT OR MATERIAL ADDITION, FLOW RETURN DATA
    NCARD='33'
    if( process_id == master_id )THEN
      CALL SEEK('C33')
      DO NWR=1,NQWR
        READ(1,*,IOSTAT=ISO) WITH_RET_GL(NWR).IQWRU,   WITH_RET_GL(NWR).JQWRU,    WITH_RET_GL(NWR).KQWRU,   WITH_RET_GL(NWR).IQWRD, WITH_RET_GL(NWR).JQWRD, WITH_RET_GL(NWR).KQWRD,  &
                             WITH_RET_GL(NWR).QWR,     WITH_RET_GL(NWR).NQWRSERQ, WITH_RET_GL(NWR).NQWRMFU, WITH_RET_GL(NWR).NQWRMFD,                                       &
                             WITH_RET_GL(NWR).BQWRMFU, WITH_RET_GL(NWR).BQWRMFD,  WITH_RET_GL(NWR).ANGWRMFD

        WRITE(7,1002)NCARD
        WRITE(7,*) WITH_RET_GL(NWR).IQWRU,   WITH_RET_GL(NWR).JQWRU,    WITH_RET_GL(NWR).KQWRU,   WITH_RET_GL(NWR).IQWRD, WITH_RET_GL(NWR).JQWRD, WITH_RET_GL(NWR).KQWRD,  &
                   WITH_RET_GL(NWR).QWR,     WITH_RET_GL(NWR).NQWRSERQ, WITH_RET_GL(NWR).NQWRMFU, WITH_RET_GL(NWR).NQWRMFD,                                       &
                   WITH_RET_GL(NWR).BQWRMFU, WITH_RET_GL(NWR).BQWRMFD,  WITH_RET_GL(NWR).ANGWRMFD
        IF( ISO > 0 ) GOTO 100
      ENDDO

      !C33A*  READ FLOW WITHDRAWAL/RETURN CONTROL
      NCARD='33A'
      CALL SEEK('C33A')
      DO NWR=1,NQWR
        WRCTL_GL(NWR).ITYPE = 0
        WRCTL_GL(NWR).ID = WITH_RET_GL(NWR).NQWRSERQ
        WRCTL_GL(NWR).SUBID = 0
        
        READ(1,*,IOSTAT=ISO) WRCTL_GL(NWR).ITYPE, WRCTL_GL(NWR).ID, ITMPU, JTMPU, ITMPD, JTMPD, &
                             WRCTL_GL(NWR).CUR.STATE, WRCTL_GL(NWR).CUR.FLOW

        WRCTL_GL(NWR).IREFUP = ITMPU
        WRCTL_GL(NWR).JREFUP = JTMPU
        WRCTL_GL(NWR).IREFDN = ITMPD
        WRCTL_GL(NWR).JREFDN = JTMPD

        WRITE(7,1002)NCARD
        WRITE(7,*) WRCTL_GL(NWR).ITYPE, WRCTL_GL(NWR).ID, WRCTL_GL(NWR).IREFUP, WRCTL_GL(NWR).JREFUP, WRCTL_GL(NWR).IREFDN, WRCTL_GL(NWR).JREFDN, &
                   WRCTL_GL(NWR).CUR.STATE, WRCTL_GL(NWR).CUR.FLOW
      ENDDO

    endif

    DO NWR=1,NQWR
      Call Broadcast_Scalar(WITH_RET_GL(NWR).IQWRU    , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).JQWRU    , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).KQWRU    , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).IQWRD    , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).JQWRD    , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).KQWRD    , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).QWR      , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).NQWRSERQ , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).NQWRMFU  , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).NQWRMFD  , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).BQWRMFU  , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).BQWRMFD  , master_id)
      Call Broadcast_Scalar(WITH_RET_GL(NWR).ANGWRMFD , master_id)
      
      Call Broadcast_Scalar(WRCTL_GL(NWR).ITYPE     , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).ID        , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).IREFUP    , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).JREFUP    , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).IREFDN    , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).JREFDN    , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).CUR.STATE , master_id)
      Call Broadcast_Scalar(WRCTL_GL(NWR).CUR.FLOW  , master_id)
    ENDDO

    !C34*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='34'
    if( process_id == master_id )THEN
      CALL SEEK('C34')
      MMAX = 3 + NDYM + NTOX
      DO NWR=1,NQWR
        READ(1,*,IOSTAT=ISO) (CQWR(NWR,MS),MS=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CQWR(NWR,MS),MS=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    !C35*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES
    !     SED(1 TO NSED),SND(1 TO NSND)
    NCARD='35'
    if( process_id == master_id )THEN
      CALL SEEK('C35')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO NWR=1,NQWR
        READ(1,*,IOSTAT=ISO) (CQWR(NWR,MS),MS=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CQWR(NWR,MS),MS=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    ENDIF

  endif

  Call Broadcast_Array(CQWR, master_id)

  IF( NSED > 0 .OR. NSND > 0 )THEN
    !C36*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
    NCARD='36'
    if( process_id == master_id )THEN
      CALL SEEK('C36')
      READ(1,*,IOSTAT=ISO)ISEDINT,ISEDBINT,NSEDFLUME,ISMUD,ISBEDMAP,ISEDVW,ISNDVW,KB,ISDTXBUG

      WRITE(7,1002)NCARD
      WRITE(7,*)ISEDINT,ISEDBINT,NSEDFLUME,ISMUD,ISBEDMAP,ISEDVW,ISNDVW,KB,ISDTXBUG
      IF( ISO > 0 ) GOTO 100
    endif

    Call Broadcast_Scalar(ISEDINT, master_id)
    Call Broadcast_Scalar(ISEDBINT, master_id)
    Call Broadcast_Scalar(NSEDFLUME, master_id)
    Call Broadcast_Scalar(ISMUD    , master_id)
    Call Broadcast_Scalar(ISBEDMAP , master_id)
    Call Broadcast_Scalar(ISEDVW   , master_id)
    Call Broadcast_Scalar(ISNDVW   , master_id)
    Call Broadcast_Scalar(KB       , master_id)
    Call Broadcast_Scalar(ISDTXBUG , master_id)

    ! *** FORCE READ IF SPATIALLY VARYING SEDIMENTS
    IF( NSEDFLUME == 3 .AND. ISEDINT < 2 ) ISEDINT = 3

    !C36A*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
    NCARD='36A'
    if( process_id == master_id )THEN
      CALL SEEK('C36A')
      READ(1,*,IOSTAT=ISO)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO

      WRITE(7,1002)NCARD
      WRITE(7,*)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO
      IF( ISO > 0 ) GOTO 100
    endif

    Call Broadcast_Scalar(ISBEDSTR, master_id)
    Call Broadcast_Scalar(ISBSDIAM, master_id)
    Call Broadcast_Scalar(ISBSDFUF, master_id)
    Call Broadcast_Scalar(COEFTSBL, master_id)
    Call Broadcast_Scalar(VISMUDST, master_id)
    Call Broadcast_Scalar(ISBKERO , master_id)

    !C36B*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS
    NCARD='36B'
    if( process_id == master_id )THEN
      CALL SEEK('C36B')
      IF( NSED>0 .OR. NSND > 0 )THEN
        READ(1,*,IOSTAT=ISO)ISEDAL,ISNDAL,IALTYP,IALSTUP,ISEDEFF,HBEDAL,COEHEFF,COEHEFF2

        WRITE(7,1002)NCARD
        WRITE(7,*)ISEDAL,ISNDAL,IALTYP,IALSTUP,HBEDAL,COEHEFF,COEHEFF2
        IF( ISO > 0 ) GOTO 100
      ENDIF
      ! *** FORCE BED ARMORING OPTION "INITIALIZATION AT STARTUP" TO BE OFF
      IF( ISRESTI > 0 ) IALSTUP = 0
    endif

    Call Broadcast_Scalar(ISEDAL  , master_id)
    Call Broadcast_Scalar(ISNDAL  , master_id)
    Call Broadcast_Scalar(IALTYP  , master_id)
    Call Broadcast_Scalar(IALSTUP , master_id)
    Call Broadcast_Scalar(HBEDAL  , master_id)
    Call Broadcast_Scalar(COEHEFF , master_id)
    Call Broadcast_Scalar(COEHEFF2, master_id)

    !C37*  BED MECHANICAL PROPERTIES PARAMETER SET 1
    NCARD='37'
    if( process_id == master_id )THEN
      CALL SEEK('C37')
      IF( NSED>0 .OR. NSND > 0 )THEN
        READ(1,*,IOSTAT=ISO) SEDSTEP, SEDSTART, IBMECH, IMORPH, HBEDMAX, BEDPORC, SEDMDMX, SEDMDMN, SEDVDRD, SEDVDRM, SEDVRDT

        WRITE(7,1002)NCARD
        WRITE(7,*) SEDSTEP,SEDSTART,IBMECH,IMORPH,HBEDMAX,BEDPORC,SEDMDMX,SEDMDMN,SEDVDRD,SEDVDRM,SEDVRDT
        IF( ISO > 0 ) GOTO 100
        
        IF( SEDSTART <= TBEGIN ) SEDSTART = TBEGIN
      ELSE
        BEDPORC = 0.4
        SEDVDRD = 1.*(1-BEDPORC)
        SEDVDRM = SEDVDRD
        SEDVDRT = 0.0
      ENDIF
    endif

    Call Broadcast_Scalar(SEDSTEP , master_id)
    Call Broadcast_Scalar(SEDSTART, master_id)
    Call Broadcast_Scalar(IBMECH  , master_id)
    Call Broadcast_Scalar(IMORPH  , master_id)
    Call Broadcast_Scalar(HBEDMAX , master_id)
    Call Broadcast_Scalar(BEDPORC , master_id)
    Call Broadcast_Scalar(SEDMDMX , master_id)
    Call Broadcast_Scalar(SEDMDMN , master_id)
    Call Broadcast_Scalar(SEDVDRD , master_id)
    Call Broadcast_Scalar(SEDVDRM , master_id)
    Call Broadcast_Scalar(SEDVRDT , master_id)

    IF( IBMECH == 0 )THEN
      SNDVDRD = BEDPORC/(1.-BEDPORC)
      SEDVDRM = SEDVDRD
    END IF

    SNDVDRD = BEDPORC/(1.-BEDPORC)
    DO NS=1,NSED
      VDRDEPO(NS)=SEDVDRD
    ENDDO
    DO NS=1,NSND
      NX = NS + NSED
      VDRDEPO(NX) = SNDVDRD
    ENDDO

    !C38*  BED MECHANICAL PROPERTIES PARAMETER SET 2
    NCARD='38'
    if( process_id == master_id )THEN
      CALL SEEK('C38')
      READ(1,*,IOSTAT=ISO)IBMECHK,BMECH1,BMECH2,BMECH3,BMECH4,BMECH5,BMECH6
      WRITE(7,1002)NCARD
      
      WRITE(7,*)IBMECHK,BMECH1,BMECH2,BMECH3,BMECH4,BMECH5,BMECH6
      IF( ISO > 0 ) GOTO 100
    endif

    Call Broadcast_Scalar(IBMECHK, master_id)
    Call Broadcast_Scalar(BMECH1 , master_id)
    Call Broadcast_Scalar(BMECH2 , master_id)
    Call Broadcast_Scalar(BMECH3 , master_id)
    Call Broadcast_Scalar(BMECH4 , master_id)
    Call Broadcast_Scalar(BMECH5 , master_id)
    Call Broadcast_Scalar(BMECH6 , master_id)

  ENDIF

  IF( NSED > 0 )THEN
    !C39*  READ COHESIVE SEDIMENT PARAMETER SET 1 REPEAT DATA LINE NSED TIMES
    NCARD='39'

    if( process_id == master_id )THEN
      CALL SEEK('C39')
      HADJ=0.0
      IF( NSED > 0 )THEN
        DO NS=1,NSED
          READ(1,*,IOSTAT=ISO)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS),ISPROBDEP(NS)

          WRITE(7,1002)NCARD
          WRITE(7,*)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISPROBDEP(NS)
          IF( ISO > 0 ) GOTO 100
          SEDDIA(NS) = 0.
          HADJ = SEDN(1)
        ENDDO
        IF( HADJ<HWET)HADJ=HWET  ! *** PMC-PROVIDE MORE CONTROL FOR !MORPH CHANGE LIMITS
      ENDIF
    endif

    Call Broadcast_Array(SEDO        , master_id)
    Call Broadcast_Array(SEDBO       , master_id)
    Call Broadcast_Array(SDEN        , master_id)
    Call Broadcast_Array(SSG         , master_id)
    Call Broadcast_Array(WSEDO       , master_id)
    Call Broadcast_Array(SEDN        , master_id)
    Call Broadcast_Array(SEXP        , master_id)
    Call Broadcast_Array(TAUD        , master_id)
    Call Broadcast_Array(ISPROBDEP   , master_id)

    !C40*  READ COHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSED TIMES
    NCARD='40'
    if( process_id == master_id )THEN
      CALL SEEK('C40')
      DO NS=1,NSED
        READ(1,*,IOSTAT=ISO)IWRSP(NS),IWRSPB(NS),WRSPO(NS),TAUR(NS),TAUN(NS),TEXP(NS),VDRRSPO(NS),COSEDHID(NS)

        WRITE(7,1002)NCARD
        WRITE(7,*)IWRSP(NS),IWRSPB(NS),WRSPO(NS),TAUR(NS),TAUN(NS),TEXP(NS),VDRRSPO(NS),COSEDHID(NS)
        IF( ISO > 0 ) GOTO 100

        IF( NS == 1 .AND. IWRSP(NS) == 999 )THEN
          WRITE(*,'(A)')'READING TAU_CRIT_COH.INP'
          OPEN(1001,FILE='tau_crit_coh.inp',STATUS='OLD')
          DO L = 2, 4393
            READ(1001,*,IOSTAT=ISO) (TAUCRCOH(L,K),K=1,10)
          ENDDO
          CLOSE(1001)
        ENDIF
        IF( ISO > 0 ) GOTO 100
        ISNDEQ(NS)=0

      ENDDO

    endif

    Call Broadcast_Array(IWRSP   , master_id)
    Call Broadcast_Array(IWRSPB   , master_id)
    Call Broadcast_Array(WRSPO    , master_id)
    Call Broadcast_Array(TAUR     , master_id)
    Call Broadcast_Array(TAUN     , master_id)
    Call Broadcast_Array(TEXP     , master_id)
    Call Broadcast_Array(VDRRSPO  , master_id)
    Call Broadcast_Array(COSEDHID , master_id)

  ENDIF

  TAUCMIN = 1000.
  ICALC_BL = 0
  SSG = 2.65         ! *** DEFAULT GRAIN DENSITY USING QUARTZ
  IF( NSND > 0 )THEN
    !C41*  READ NONCOHESIVE SEDIMENT PARAMETER SET 1 REPEAT DATA LINE NSND TIMES
    NCARD='41'

    if( process_id == master_id )THEN
      CALL SEEK('C41')
      DO NX=1,NSND
        NS=NX+NSED
        READ(1,*,IOSTAT=ISO)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),SEDDIA(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS)

        ! *** IF SETTLING VELOCITY IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
        IF( WSEDO(NS) < 0.0 )THEN
          WSEDO(NS) = SETSTVEL(SEDDIA(NS),SSG(NS))
        ENDIF

        WRITE(7,1002)NCARD
        WRITE(7,*)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),SEDDIA(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(SEDO     , master_id)
    Call Broadcast_Array(SEDBO    , master_id)
    Call Broadcast_Array(SDEN     , master_id)
    Call Broadcast_Array(SSG      , master_id)
    Call Broadcast_Array(SEDDIA   , master_id)
    Call Broadcast_Array(WSEDO    , master_id)
    Call Broadcast_Array(SEDN     , master_id)
    Call Broadcast_Array(SEXP     , master_id)
    Call Broadcast_Array(TAUD     , master_id)
    Call Broadcast_Array(ISEDSCOR , master_id)

    !C42*  READ NONCOHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSND TIMES
    NCARD='42'
    if( process_id == master_id )THEN
      CALL SEEK('C42')
      DO NX=1,NSND
        NS=NX+NSED
        READ(1,*,IOSTAT=ISO)ISNDEQ(NS),ISBDLD(NS),TAUR(NS),TAUN(NS),TCSHIELDS(NS),ISLTAUC(NS),IBLTAUC(NS),IROUSE(NX),ISNDM1(NX),ISNDM2(NX),RSNDM(NX)

        ! IF TAUR(NS) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
        !       TAUR:     CRITICAL SHIELDS STRESS IN (M/S)**2   (ISNDEQ=2)
        !       TAUN:     EQUAL TO TAUR FOR NONCHOESIVE SED TRANS  (ISNDEQ=2)
        !       TEXP:     CRITICAL SHIELDS PARAMETER  (ISNDEQ=2)
        DSTR=0.0
        USTR=0.0

        ! *** IF TAUR(NS) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
        IF( TAUR(NS)<0.0 )THEN
          CALL SETSHLD(TAUR(NS),TCSHIELDS(NS),SEDDIA(NS),SSG(NS),DSTR,USTR)
          TAUN(NS)=TAUR(NS)
        ENDIF
        TAUCMIN = MIN(TAUCMIN,TAUR(NS))

        WRITE(7,1002)NCARD
        WRITE(7,*)ISNDEQ(NS),TAUR(NS),TAUN(NS),TCSHIELDS(NS),SEDDIA(NS),SSG(NS),DSTR,USTR
        IF( ISO > 0 ) GOTO 100
        
        IWRSP(NS) = 0
        WRSPO(NS) = 0
      ENDDO

    endif

    Call Broadcast_Array(ISNDEQ         , master_id)
    Call Broadcast_Array(ISBDLD         , master_id)
    Call Broadcast_Array(TAUR           , master_id)
    Call Broadcast_Array(TAUN           , master_id)
    Call Broadcast_Array(TCSHIELDS      , master_id)
    Call Broadcast_Array(ISLTAUC        , master_id)
    Call Broadcast_Array(IBLTAUC        , master_id)
    Call Broadcast_Array(IROUSE         , master_id)
    Call Broadcast_Array(TCSHIELDS      , master_id)
    Call Broadcast_Array(ISNDM1         , master_id)
    Call Broadcast_Array(ISNDM2         , master_id)
    Call Broadcast_Array(RSNDM          , master_id)
    Call Broadcast_Array(SSG            , master_id)
    Call Broadcast_Scalar(DSTR          , master_id)
    Call Broadcast_Scalar(USTR          , master_id)
    Call Broadcast_Scalar(TAUCMIN       , master_id)

    !C42A*  READ NONCOHESIVE SEDIMENT BED LOAD PARAMETERS
    NCARD='42A'
    if( process_id == master_id )THEN
      CALL SEEK('C42A')
      DO NS=1,NSND
        READ(1,*,IOSTAT=ISO)ICALC_BL,SBDLDA(NS),SBDLDB(NS),SBDLDG1(NS),SBDLDG2(NS),SBDLDG3(NS),SBDLDG4(NS),SBDLDP(NS),ISBLFUC,BLBSNT

        WRITE(7,1002)NCARD
        WRITE(7,*)ICALC_BL,SBDLDA(NS),SBDLDB(NS),SBDLDG1(NS),SBDLDG2(NS),SBDLDG3(NS),SBDLDG4(NS),SBDLDP(NS),ISBLFUC,BLBSNT
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Scalar(ICALC_BL , master_id)
    Call Broadcast_Scalar(ISBLFUC  , master_id)
    Call Broadcast_Scalar(BLBSNT   , master_id)
    Call Broadcast_Array(SBDLDA    , master_id)
    Call Broadcast_Array(SBDLDB    , master_id)
    Call Broadcast_Array(SBDLDG1   , master_id)
    Call Broadcast_Array(SBDLDG2   , master_id)
    Call Broadcast_Array(SBDLDG3   , master_id)
    Call Broadcast_Array(SBDLDG4   , master_id)
    Call Broadcast_Array(SBDLDP    , master_id)

  ENDIF

  IF( NTOX > 0 )THEN
    !C43A*  READ TOXIC CONTAMINANT INITIAL CONDITIONS
    NCARD='43A'
    if( process_id == master_id )THEN
      CALL SEEK('C43A')
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM,ITXINT(NT),ITXBDUT(NT),TOXINTW(NT),TOXINTB(NT)
        
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM,ITXINT(NT),ITXBDUT(NT),TOXINTW(NT),TOXINTB(NT)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(ITXINT , master_id)
    Call Broadcast_Array(ITXBDUT, master_id)
    Call Broadcast_Array(TOXINTW, master_id)
    Call Broadcast_Array(TOXINTB, master_id)

    !C43B*  READ TOXIC KINETIC OPTION FLAGS
    NCARD='43B'
    if( process_id == master_id )THEN
      CALL SEEK('C43B')
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO) NDUM,ITOXKIN(1,NT),ITOXKIN(2,NT),ITOXKIN(3,NT),ITOXKIN(4,NT),ITOXKIN(5,NT)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) NDUM,ITOXKIN(1,NT),ITOXKIN(2,NT),ITOXKIN(3,NT),ITOXKIN(4,NT),ITOXKIN(5,NT)
        IF( ISO > 0 ) GOTO 100
      ENDDO

    endif

    Call Broadcast_Array(ITOXKIN, master_id)
    Call Broadcast_Array(ITOXKIN, master_id)
    Call Broadcast_Array(ITOXKIN, master_id)
    Call Broadcast_Array(ITOXKIN, master_id)
    Call Broadcast_Array(ITOXKIN, master_id)

    !C43C*  READ TOXIC TIMING AND VOLATILIZATION SWITCHES
    NCARD='43C'
    if( process_id == master_id )THEN
      CALL SEEK('C43C')
      READ(1,*,IOSTAT=ISO) TOXSTEPW, TOXSTEPB, TOX_VEL_MAX, TOX_DEP_MAX,ITOXTEMP,TOXTEMP
      WRITE(7,1002)NCARD
      WRITE(7,*) TOXSTEPW, TOXSTEPB, TOX_VEL_MAX, TOX_DEP_MAX,ITOXTEMP,TOXTEMP
      IF( ISO > 0 ) GOTO 100
      IF( TOXSTEPW < SEDSTEP ) TOXSTEPW = SEDSTEP    ! *** TOXIC KINETICS CANNOT OPERATE ON A FINER TIMESCALE THAN SEDIMENTS
      IF( TOXSTEPB < SEDSTEP ) TOXSTEPB = SEDSTEP    ! *** TOXIC BED MIXING AND DIFFUSION CANNOT OPERATE ON A FINER TIMESCALE THAN SEDIMENTS
      IF( ISTRAN(5) > 0 .AND. ISTRAN(2) == 0 .AND. ( ITOXTEMP < 0 .OR. ITOXTEMP-1 > NCSER(2) ) ) CALL STOPP('ITOXTEMP IS OUT OF RANGE')

    endif

    Call Broadcast_Scalar(TOXSTEPW      , master_id)
    Call Broadcast_Scalar(TOXSTEPB      , master_id)
    Call Broadcast_Scalar(TOX_VEL_MAX   , master_id)
    Call Broadcast_Scalar(TOX_DEP_MAX   , master_id)
    Call Broadcast_Scalar(ITOXTEMP      , master_id)
    Call Broadcast_Scalar(TOXTEMP       , master_id)

    !C43D*  READ TOXIC BULK DECAY AND BIODEGRADATION PARAMETERS
    NCARD='43D'
    if( process_id == master_id )THEN
      CALL SEEK('C43D')
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO) NDUM,TOX_BLK_KW(NT), TOX_BLK_KB(NT), TOX_BLK_MXD(NT), TOX_BIO_KW(NT), TOX_BIO_KB(NT), TOX_BIO_MXD(NT), &
          TOX_BIO_Q10W(NT), TOX_BIO_Q10B(NT), TOX_BIO_TW(NT), TOX_BIO_TB(NT)
        WRITE(7,1002)NCARD
        WRITE(7,*)  NDUM,TOX_BLK_KW(NT), TOX_BLK_KB(NT), TOX_BIO_MXD(NT), TOX_BIO_KW(NT), TOX_BIO_KB(NT), TOX_BIO_MXD(NT), &
          TOX_BIO_Q10W(NT), TOX_BIO_Q10B(NT), TOX_BIO_TW(NT), TOX_BIO_TB(NT)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif


    Call Broadcast_Array(TOX_BLK_KW     , master_id)
    Call Broadcast_Array(TOX_BLK_KB     , master_id)
    Call Broadcast_Array(TOX_BIO_MXD    , master_id)
    Call Broadcast_Array(TOX_BIO_KW      , master_id)
    Call Broadcast_Array(TOX_BIO_KB      , master_id)
    Call Broadcast_Array(TOX_BIO_MXD     , master_id)
    Call Broadcast_Array(TOX_BIO_Q10W   , master_id)
    Call Broadcast_Array(TOX_BIO_Q10B   , master_id)
    Call Broadcast_Array(TOX_BIO_TW     , master_id)
    Call Broadcast_Array(TOX_BIO_TB   , master_id)

    !C43D*  READ TOXIC VOLATILIZATION PARAMETERS
    NCARD='43E'
    if( process_id == master_id )THEN
      CALL SEEK('C43E')
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO) NDUM,TOX_MW(NT), TOX_HE(NT), TOX_KV_TCOEFF(NT), TOX_ATM(NT), TOX_ADJ(NT)
        WRITE(7,1002)NCARD
        
        WRITE(7,*) NDUM,TOX_MW(NT), TOX_HE(NT), TOX_KV_TCOEFF(NT), TOX_ATM(NT), TOX_ADJ(NT)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif


    Call Broadcast_Array(TOX_MW         , master_id)
    Call Broadcast_Array(TOX_HE         , master_id)
    Call Broadcast_Array(TOX_KV_TCOEFF  , master_id)
    Call Broadcast_Array(TOX_ATM        , master_id)
    Call Broadcast_Array(TOX_ADJ        , master_id)

    ! *** Delme - Photolysis   RKTOXP(NT),SKTOXP(NT)

    !C44*  READ TOXIC CONTAMINANT PARAMETERS: SORBTION
    NCARD='44'
    if( process_id == master_id )THEN
      CALL SEEK('C44')
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO) NDUM, ISTOC(NT), DIFTOX(NT), DIFTOXS(NT), PDIFTOX(NT), DPDIFTOX(NT)

        WRITE(7,1002)NCARD
        WRITE(7,* )NDUM, ISTOC(NT), DIFTOX(NT), DIFTOXS(NT), PDIFTOX(NT), DPDIFTOX(NT)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(ISTOC   , master_id)
    Call Broadcast_Array(DIFTOX  , master_id)
    Call Broadcast_Array(DIFTOXS  , master_id)
    Call Broadcast_Array(PDIFTOX  , master_id)
    Call Broadcast_Array(DPDIFTOX , master_id)

    DO NT=1,NTOX
      ISPMXZ(NT)=0
      IF( PDIFTOX(NT) < 0.0 ) ISPMXZ(NT)=1
      ISDIFBW(NT)=0
      IF( DIFTOXS(NT)<0.0 )THEN
        DIFTOXS(NT)=ABS(DIFTOXS(NT))
        ISDIFBW(NT)=1
      ENDIF
    ENDDO

    !C45*  READ TOXIC CONTAMINANT-SEDIMENT INTERACTION PARAMETERS
    NCARD='45'
    if( process_id == master_id )THEN
      CALL SEEK('C45')
      DO NT=1,NTOX
        IF( NSED > 0 )THEN
          DO NS=1,NSED
            READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
            
            WRITE(7,1002)NCARD
            WRITE(7,*)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
            IF( ISO > 0 ) GOTO 100
          ENDDO
        ENDIF
        IF( NSND > 0 )THEN
          DO NX=1,NSND
            NS=NX+NSED
            READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
            
            WRITE(7,1002)NCARD
            WRITE(7,*)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
            IF( ISO > 0 ) GOTO 100
          ENDDO
        ENDIF
      ENDDO
    endif

    Call Broadcast_Array(ITXPARW, master_id)
    Call Broadcast_Array(TOXPARW, master_id)
    Call Broadcast_Array(CONPARW, master_id)
    Call Broadcast_Array(ITXPARB, master_id)
    Call Broadcast_Array(TOXPARB, master_id)
    Call Broadcast_Array(CONPARB, master_id)

    !C45A*  READ TOXIC CONTAMINANT ORGANIC CARBON PARAMETERS
    NCARD='45A'
    if( process_id == master_id )THEN
      CALL SEEK('C45A')
      IF( NTOX > 0 )THEN
        READ(1,*,IOSTAT=ISO)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
        IF( ISTRAN(5) == 0 )THEN
          ISTDOCW=0
          ISTPOCW=0
          ISTDOCB=0
          ISTPOCB=0
        ENDIF
        WRITE(7,1002)NCARD
        WRITE(7,*)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
        IF( ISO > 0 ) GOTO 100
      ENDIF
    endif

    Call Broadcast_Scalar(ISTDOCW,  master_id)
    Call Broadcast_Scalar(ISTPOCW,  master_id)
    Call Broadcast_Scalar(ISTDOCB,  master_id)
    Call Broadcast_Scalar(ISTPOCB,  master_id)
    Call Broadcast_Scalar(STDOCWC,  master_id)
    Call Broadcast_Scalar(STPOCWC,  master_id)
    Call Broadcast_Scalar(STDOCBC,  master_id)
    Call Broadcast_Scalar(STPOCBC,  master_id)

    !C45B* READ TOXIC CONTAMINANT-ORGANIC CARBON INTERACTION PARAMETERS
    NCARD='45B'
    if( process_id == master_id )THEN
      CALL SEEK('C45B')
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARWC(1,NT),TOXPARWC(1,NT),CONPARWC(1,NT),ITXPARBC(1,NT),TOXPARBC(1,NT),CONPARBC(1,NT)
        
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM1,NDUM2,ITXPARWC(1,NT),TOXPARWC(1,NT),CONPARWC(1,NT),ITXPARBC(1,NT),TOXPARBC(1,NT),CONPARBC(1,NT)
        IF( ISO > 0 ) GOTO 100
        
        READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARWC(2,NT),TOXPARWC(2,NT),CONPARWC(2,NT),ITXPARBC(2,NT),TOXPARBC(2,NT),CONPARBC(2,NT)
        
        WRITE(7,1002)NCARD
        WRITE(7,*)NDUM1,NDUM2,ITXPARWC(2,NT),TOXPARWC(2,NT),CONPARWC(2,NT),ITXPARBC(2,NT),TOXPARBC(2,NT),CONPARBC(2,NT)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(ITXPARWC,  master_id)
    Call Broadcast_Array(TOXPARWC,  master_id)
    Call Broadcast_Array(CONPARWC,  master_id)
    Call Broadcast_Array(ITXPARBC,  master_id)
    Call Broadcast_Array(TOXPARBC,  master_id)
    Call Broadcast_Array(CONPARBC,  master_id)

    !C45C* READ TOXIC CONTAMINANT-ORGANIC CARBON WATER COLUMN POC FRACTIONS
    NCARD='45C'
    if( process_id == master_id )THEN
      CALL SEEK('C45C')
      WRITE(7,1002)NCARD
      NTMP=NSED+NSND
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM1,(FPOCWST(NS,NT),NS=1,NTMP)
        
        WRITE(7,*)NDUM1,(FPOCWST(NS,NT),NS=1,NTMP)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(FPOCWST, master_id)
    
    ! RESET INORGANIC SEDIMENT PARTITION COEFFICIENTS BASED ON
    ! FRACTION OF POC ASSIGNED TO EACH SEDIMENT CLASS IN WATER COLUMN
    DO NT=1,NTOX
      IF( ISTOC(NT) >= 2 )THEN
        IF( NSED > 0 )THEN
          DO NS=1,NSED
            ITXPARW(NS,NT)=0
            TOXPARW(NS,NT)=TOXPARWC(2,NT)
            CONPARW(NS,NT)=0.
          ENDDO
        ENDIF
        IF( NSND > 0 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            ITXPARW(NS,NT)=0
            TOXPARW(NS,NT)=TOXPARWC(2,NT)
            CONPARW(NS,NT)=0.
          ENDDO
        ENDIF
      ENDIF
    ENDDO

    !C45D* READ TOXIC CONTAMINANT-ORGANIC CARBON SED BED POC FRACTIONS
    NCARD='45D'
    if( process_id == master_id )THEN
      CALL SEEK('C45D')
      WRITE(7,1002)NCARD
      NTMP = NSED+NSND
      DO NT=1,NTOX
        READ(1,*,IOSTAT=ISO)NDUM1,(FPOCBST(NS,NT),NS=1,NTMP)
        
        IF( ISO > 0 ) GOTO 100
        WRITE(7,*)NDUM1,(FPOCBST(NS,NT),NS=1,NTMP)
      ENDDO
    endif
    Call Broadcast_Array(FPOCBST, master_id)

    !C45E* READ TOXIC CONTAMINANT ATMOSPHERIC DEPOSITIONS SETTINGS
    NCARD='45E'
    if( process_id == master_id )THEN
      CALL SEEK('C45E')
      WRITE(7,1002)NCARD

      DO NT=1,NTOX
        READ(1,*,ERR=100) NDUM, TOXDEP(NT).ITXDRY, TOXDEP(NT).TXDRY, TOXDEP(NT).ITXWET, TOXDEP(NT).TXWET

        TOXDEP(NT).TXDRY = TOXDEP(NT).TXDRY/86400.   ! *** CONVERT FROM MG/M2/DAY TO MG/M2/SEC
      ENDDO
    endif
    DO NT=1,NTOX
      Call Broadcast_Scalar(TOXDEP(NT).ITXDRY, master_id)
      Call Broadcast_Scalar(TOXDEP(NT).TXDRY,  master_id)
      Call Broadcast_Scalar(TOXDEP(NT).ITXWET, master_id)
      Call Broadcast_Scalar(TOXDEP(NT).TXWET,  master_id)
    ENDDO
    
    ! RESET INORGANIC SEDIMENT PARTITION COEFFICIENTS BASED ON
    ! FRACTION OF POC ASSIGNED TO EACH SEDIMENT CLASS IN SEDIMENT BED
    DO NT=1,NTOX
      IF( ISTOC(NT) >= 2 )THEN
        IF( NSED > 0 )THEN
          DO NS=1,NSED
            ITXPARB(NS,NT)=0
            TOXPARB(NS,NT)=TOXPARBC(2,NT)
            CONPARB(NS,NT)=0
          ENDDO
        ENDIF
        IF( NSND > 0 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            ITXPARB(NS,NT)=0
            TOXPARB(NS,NT)=TOXPARBC(2,NT)
            CONPARB(NS,NT)=0
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF  ! *** END OF NTOX>0 BLOCK

  !C46*  READ BUOYANCY, TEMPERATURE, DYE DATA AND CONCENTRATION BC DATA
  NCARD='46'
  if( process_id == master_id )THEN
    CALL SEEK('C46')
    READ(1,*,IOSTAT=ISO) BSC, TEMO, HEQT, RKDYE, NCBS, NCBW, NCBE, NCBN

    WRITE(7,1002)NCARD
    WRITE(7,*) BSC, TEMO, HEQT, RKDYE, NCBS, NCBW, NCBE, NCBN
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(BSC   , master_id)
  Call Broadcast_Scalar(TEMO  , master_id)
  Call Broadcast_Scalar(HEQT  , master_id)
  Call Broadcast_Scalar(RKDYE , master_id)
  Call Broadcast_Scalar(NCBS  , master_id)
  Call Broadcast_Scalar(NCBW  , master_id)
  Call Broadcast_Scalar(NCBE  , master_id)
  Call Broadcast_Scalar(NCBN  , master_id)

  IF( BSC == 2. )THEN
    BSC = 1.
    IBSC = 1
  ELSE
    IBSC = 0
  ENDIF
  
  ! *** DISABLE BOUYANCY IF ALL CONSTITUENTS IMPACTING DENSITY ARE OFF (7.2)
  IF( BSC > 0. .AND. (ISTRAN(1) < 1 .AND. ISTRAN(2) < 1 .AND.  .NOT. (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)) )THEN
    BSC = 0.
    IBSC = 0
  ENDIF
  
  IF( TEMO < 0.0 )THEN
    TEMO = ABS(TEMO)
    INITTEMP = 1
  ELSE
    INITTEMP = 0
  ENDIF

  IF( ISICE > 0 )THEN
    !C46A*   READ ICE EFFECTS
    NCARD='46A'
    if( process_id == master_id )THEN
      CALL SEEK('C46A')
      READ(1,*,IOSTAT=ISO)ISICE,NISER,TEMPICE,CDICE,ICETHMX,RICETHK0
      
      WRITE(7,1002)NCARD
      WRITE(7,*)ISICE,NISER,TEMPICE,CDICE,ICETHMX,RICETHK0
      IF( ISO > 0 ) GOTO 100
    ENDIF
  endif
  Call Broadcast_Scalar(ISICE   , master_id)
  Call Broadcast_Scalar(NISER   , master_id)
  Call Broadcast_Scalar(TEMPICE , master_id)
  Call Broadcast_Scalar(CDICE   , master_id)
  Call Broadcast_Scalar(ICETHMX , master_id)
  Call Broadcast_Scalar(RICETHK0, master_id)

  RHOI = 917.
  RISEVEL = 0.01
  IF( ISICE > 2 )THEN
    !C46B*   READ ICE MODULE PARAMETERS
    NCARD='46B'
    if( process_id == master_id )THEN
      CALL SEEK('C46B')
      READ(1,*,IOSTAT=ISO) HWI,ICEK,ALBEDOI,BETAI,GAMMAI,MINICETHICK,RHOI,ISRHEVAP,RISEVEL,MELTFACTOR,AFWI,BFWI,CFWI
      WRITE(7,1002) NCARD
      
      WRITE(7,*) HWI,ICEK,ALBEDOI,BETAI,GAMMAI,MINICETHICK,ICETHMX,RHOI,ISRHEVAP,RISEVEL,MELTFACTOR,AFWI,BFWI,CFWI
      IF( ISO > 0 ) GOTO 100
    endif
    Call Broadcast_Scalar(HWI           , master_id)
    Call Broadcast_Scalar(ICEK          , master_id)
    Call Broadcast_Scalar(ALBEDOI       , master_id)
    Call Broadcast_Scalar(BETAI         , master_id)
    Call Broadcast_Scalar(GAMMAI        , master_id)
    Call Broadcast_Scalar(MINICETHICK   , master_id)
    Call Broadcast_Scalar(ICETHMX       , master_id)
    Call Broadcast_Scalar(RHOI          , master_id)
    Call Broadcast_Scalar(ISRHEVAP      , master_id)
    Call Broadcast_Scalar(RISEVEL       , master_id)
    Call Broadcast_Scalar(MELTFACTOR    , master_id)
    Call Broadcast_Scalar(AFWI          , master_id)
    Call Broadcast_Scalar(BFWI          , master_id)
    Call Broadcast_Scalar(CFWI          , master_id)

    IF( CFWI <= 0. )THEN
      AFWI = 9.2    ! *** WIND FUNCTION A, Edinger, et. al. (1974)
      BFWI = 0.46   ! *** WIND FUNCTION B, Edinger, et. al. (1974)
      CFWI = 2.0    ! *** WIND FUNCTION C, Edinger, et. al. (1974)
    ELSE
      ! *** Units Conversion = 1 mmhg = 1.33322 millibars
      AFWI = AFWI*1.33322
      BFWI = BFWI*1.33322
    ENDIF

    ! *** FRAZIL ICE TRANSPORT DEFAULTS
    ISCDCA(10) = 0  ! *** UPWIND DIFFERENCE (3TL ONLY)
    ISFCT(10)  = 0

  ENDIF
  IF( ISICE == 2) NISER = 1

  IF( NASER > 0 )THEN
    ! ** READING TWO NEW CARDS FOR EFDC_073 INSTAED OF ASER.INP FOR 072
    DS_LAT=0.0
    DS_LONG=0.0
    COMPUTESOLRAD=.FALSE.

    NCARD='46C'
    if( process_id == master_id )THEN
      CALL SEEK('C46C')
      READ(1,*,IOSTAT=ISO) DS_LONG, DS_LAT, COMPUTESOLRAD, USESHADE, IEVAP, WINDFA, WINDFB, WINDFC
      
      WRITE(7,1002) NCARD
      WRITE(7,*) DS_LONG, DS_LAT, COMPUTESOLRAD, USESHADE, IEVAP, WINDFA, WINDFB, WINDFC
      IF( ISO > 0 ) GOTO 100
    endif

    Call Broadcast_Scalar(DS_LONG       , master_id)
    Call Broadcast_Scalar(DS_LAT        , master_id)
    Call Broadcast_Scalar(COMPUTESOLRAD , master_id)
    Call Broadcast_Scalar(USESHADE      , master_id)
    Call Broadcast_Scalar(IEVAP         , master_id)
    Call Broadcast_Scalar(WINDFA        , master_id)
    Call Broadcast_Scalar(WINDFB        , master_id)
    Call Broadcast_Scalar(WINDFC        , master_id)

    AFW = 9.2        ! *** Edinger, et. al. (1974)  W/m2/mmHg
    BFW = 0.46       ! *** Edinger, et. al. (1974)  W/m2/mmHg
    CFW = 2.0        ! *** Edinger, et. al. (1974)  W/m2/mmHg

    ! *** CONVERT WIND FACTOR COEFFICIENTS FROM W/M2/MILLIBAR TO M/S/MILLIBAR
    ! *** Latent Heat of Evaporation = 2259 KJ/KG
    ! *** Density of Water           = 1000 kg/m3
    WINDFA = WINDFA*4.42674E-10
    WINDFB = WINDFB*4.42674E-10
    WINDFC = WINDFC*4.42674E-10

    NCARD='46D'
    if( process_id == master_id )THEN
      CALL SEEK('C46D')
      ! *** TBEDIT deprecated.  Changed to TEMBO
      ! *** DABEDT deprecated.  Changed to TEMTHKO
      READ(1,*,IOSTAT=ISO) IASWRAD, REVC, RCHC, ISVHEAT, SWRATNF, SWRATNS, FSWRATF, TEMTHKO, TEMBO, HTBED1, HTBED2, WQKEB(1), WQKETSS
      WRITE(7,1002) NCARD
      
      WRITE(7,*) IASWRAD, REVC, RCHC, ISVHEAT, SWRATNF, SWRATNS, FSWRATF, TEMTHKO, TEMBO, HTBED1, HTBED2, WQKEB(1), WQKETSS
      IF( ISO > 0 ) GOTO 100
    endif

    Call Broadcast_Scalar(IASWRAD   , master_id)
    Call Broadcast_Scalar(REVC      , master_id)
    Call Broadcast_Scalar(RCHC      , master_id)
    Call Broadcast_Scalar(ISVHEAT   , master_id)
    Call Broadcast_Scalar(SWRATNF   , master_id)
    Call Broadcast_Scalar(SWRATNS   , master_id)
    Call Broadcast_Scalar(FSWRATF   , master_id)
    Call Broadcast_Scalar(TEMTHKO    , master_id)
    Call Broadcast_Scalar(TEMBO    , master_id)
    Call Broadcast_Scalar(HTBED1    , master_id)
    Call Broadcast_Scalar(HTBED2    , master_id)
    Call Broadcast_Scalar(WQKEB(1)  , master_id)
    Call Broadcast_Scalar(WQKETSS   , master_id)

  ENDIF

  ! *** READ DYE TYPES
  IF( ISTRAN(3) > 0 .AND. NDYE > 0 )THEN
    NCARD='46E'
    if( process_id == master_id )THEN
      CALL SEEK('C46E')
      DO MD=1,NDYE
        READ(1,*,IOSTAT=ISO) NDUM1, DYES(MD).ITYPE, DYES(MD).KRATE0, DYES(MD).KRATE1, DYES(MD).TADJ, DYES(MD).TREF, DYES(MD).SETTLE
        IF( ISO > 0 ) GOTO 100
        
        WRITE(7,*) NDUM1, DYES(MD).ITYPE, DYES(MD).KRATE0, DYES(MD).KRATE1, DYES(MD).TADJ, DYES(MD).TREF, DYES(MD).SETTLE
      ENDDO
    endif

    Do MD=1,NDYE
      Call Broadcast_Scalar(DYES(MD).ITYPE    , master_id)
      Call Broadcast_Scalar(DYES(MD).KRATE0   , master_id)
      Call Broadcast_Scalar(DYES(MD).KRATE1   , master_id)
      Call Broadcast_Scalar(DYES(MD).TADJ     , master_id)
      Call Broadcast_Scalar(DYES(MD).TREF     , master_id)
      Call Broadcast_Scalar(DYES(MD).SETTLE   , master_id)
    ENDDO

  ENDIF

  ! ****************************************************************************C
  IF( NCBS > 0 )THEN

    !C47*  READ LOCATIONS OF CONC BC'S ON SOUTH BOUNDARIES
    NCARD='47'
    if( process_id == master_id )THEN
      CALL SEEK('C47')
      DO L=1,NCBS
        READ(1,*,IOSTAT=ISO) ICBS_GL(L),JCBS_GL(L),NTSCRS_GL(L),NCSERS_GL(L,1),NCSERS_GL(L,2),NCSERS_GL(L,3),NCSERS_GL(L,4),NCSERS_GL(L,5),NCSERS_GL(L,6),NCSERS_GL(L,7)
        WRITE(7,1002)NCARD
        
        WRITE(7,*) ICBS_GL(L),JCBS_GL(L),NTSCRS_GL(L),NCSERS_GL(L,1),NCSERS_GL(L,2),NCSERS_GL(L,3),NCSERS_GL(L,4),NCSERS_GL(L,5),NCSERS_GL(L,6),NCSERS_GL(L,7)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(ICBS_GL  , master_id)
    Call Broadcast_Array(JCBS_GL  , master_id)
    Call Broadcast_Array(NTSCRS_GL, master_id)
    Call Broadcast_Array(NCSERS_GL, master_id)

    !C48*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='48'
    if( process_id == master_id )THEN
      CALL SEEK('C48')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBS
        READ(1,*,IOSTAT=ISO) (CBS_GL(L,1,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBS_GL(L,1,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBS_GL, master_id)

    !C49*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='49'
    if( process_id == master_id )THEN
      CALL SEEK('C49')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBS
        READ(1,*,IOSTAT=ISO) (CBS_GL(L,1,M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBS_GL(L,1,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBS_GL, master_id)

    !C50*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='50'
    if( process_id == master_id )THEN
      CALL SEEK('C50')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBS
        READ(1,*,IOSTAT=ISO) (CBS_GL(L,2,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBS_GL(L,2,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBS_GL, master_id)

    !C51*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='51'
    if( process_id == master_id )THEN
      CALL SEEK('C51')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBS
        READ(1,*,IOSTAT=ISO) (CBS_GL(L,2,M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBS_GL(L,2,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBS_GL, master_id)

  ENDIF


  IF( NCBW > 0 )THEN
    !C52*  READ LOCATIONS OF CONC BC'S ON WEST BOUNDARIES
    NCARD='52'
    if( process_id == master_id )THEN
      CALL SEEK('C52')
      DO L=1,NCBW
        READ(1,*,IOSTAT=ISO) ICBW_GL(L),JCBW_GL(L),NTSCRW_GL(L),NCSERW_GL(L,1),NCSERW_GL(L,2),NCSERW_GL(L,3),NCSERW_GL(L,4),NCSERW_GL(L,5),NCSERW_GL(L,6),NCSERW_GL(L,7)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) ICBW_GL(L),JCBW_GL(L),NTSCRW_GL(L),NCSERW_GL(L,1),NCSERW_GL(L,2),NCSERW_GL(L,3),NCSERW_GL(L,4),NCSERW_GL(L,5),NCSERW_GL(L,6),NCSERW_GL(L,7)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(ICBW_GL, master_id)
    Call Broadcast_Array(JCBW_GL, master_id)
    Call Broadcast_Array(NTSCRW_GL, master_id)
    Call Broadcast_Array(NCSERW_GL, master_id)

    !C53*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='53'
    if( process_id == master_id )THEN
      CALL SEEK('C53')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBW
        READ(1,*,IOSTAT=ISO) (CBW_GL(L,1,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBW_GL(L,1,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBW_GL, master_id)

    !C54*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='54'
    if( process_id == master_id )THEN
      CALL SEEK('C54')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBW
        READ(1,*,IOSTAT=ISO) (CBW_GL(L,1,M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBW_GL(L,1,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBW_GL, master_id)

    !C55*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='55'
    if( process_id == master_id )THEN
      CALL SEEK('C55')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBW
        READ(1,*,IOSTAT=ISO) (CBW_GL(L,2,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBW_GL(L,2,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBW_GL, master_id)

    !C56*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='56'
    if( process_id == master_id )THEN
      CALL SEEK('C56')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBW
        READ(1,*,IOSTAT=ISO) (CBW_GL(L,2,M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBW_GL(L,2,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBW_GL, master_id)

  ENDIF

  IF( NCBE > 0 )THEN
    !C57*  READ LOCATIONS OF CONC BC'S ON EAST BOUNDARIES
    NCARD='57'
    if( process_id == master_id )THEN
      CALL SEEK('C57')
      DO L=1,NCBE
        READ(1,*,IOSTAT=ISO) ICBE_GL(L),JCBE_GL(L),NTSCRE_GL(L),NCSERE_GL(L,1),NCSERE_GL(L,2),NCSERE_GL(L,3),NCSERE_GL(L,4),NCSERE_GL(L,5),NCSERE_GL(L,6),NCSERE_GL(L,7)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) ICBE_GL(L),JCBE_GL(L),NTSCRE_GL(L),NCSERE_GL(L,1),NCSERE_GL(L,2),NCSERE_GL(L,3),NCSERE_GL(L,4),NCSERE_GL(L,5),NCSERE_GL(L,6),NCSERE_GL(L,7)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(ICBE_GL,   master_id)
    Call Broadcast_Array(JCBE_GL,   master_id)
    Call Broadcast_Array(NCSERE_GL, master_id)
    Call Broadcast_Array(NTSCRE_GL,  master_id)

    !C58*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='58'
    if( process_id == master_id )THEN
      CALL SEEK('C58')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBE
        READ(1,*,IOSTAT=ISO) (CBE_GL(L,1,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBE_GL(L,1,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBE_GL, master_id)

    !C59*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='59'
    if( process_id == master_id )THEN
      CALL SEEK('C59')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBE
        READ(1,*,IOSTAT=ISO) (CBE_GL(L,1,M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBE_GL(L,1,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBE_GL, master_id)
    
    !C60*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='60'
    if( process_id == master_id )THEN
      CALL SEEK('C60')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBE
        READ(1,*,IOSTAT=ISO) (CBE_GL(L,2,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBE_GL(L,2,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBE_GL, master_id)

    !C61*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='61'
    if( process_id == master_id )THEN
      CALL SEEK('C61')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBE
        READ(1,*,IOSTAT=ISO) (CBE_GL(L,2,M),M=MMIN,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBE_GL(L,2,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBE_GL, master_id)

  ENDIF    ! *** NCBE>0

  IF( NCBN > 0 )THEN
    !C62*  READ LOCATIONS OF CONC BC'S ON NORTH BOUNDARIES
    NCARD='62'
    if( process_id == master_id )THEN
      CALL SEEK('C62')
      DO L=1,NCBN
        READ(1,*,IOSTAT=ISO) ICBN_GL(L),JCBN_GL(L),NTSCRN_GL(L),NCSERN_GL(L,1),NCSERN_GL(L,2),NCSERN_GL(L,3),NCSERN_GL(L,4),NCSERN_GL(L,5),NCSERN_GL(L,6),NCSERN_GL(L,7)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) ICBN_GL(L),JCBN_GL(L),NTSCRN_GL(L),NCSERN_GL(L,1),NCSERN_GL(L,2),NCSERN_GL(L,3),NCSERN_GL(L,4),NCSERN_GL(L,5),NCSERN_GL(L,6),NCSERN_GL(L,7)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(ICBN_GL, master_id)
    Call Broadcast_Array(JCBN_GL, master_id)
    Call Broadcast_Array(NTSCRN_GL, master_id)
    Call Broadcast_Array(NCSERN_GL, master_id)

    !C63*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='63'
    if( process_id == master_id )THEN
      CALL SEEK('C63')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBN
        READ(1,*,IOSTAT=ISO) (CBN_GL(L,1,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBN_GL(L,1,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO

    endif
    Call Broadcast_Array(CBN_GL, master_id)

    !C64*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='64'
    if( process_id == master_id )THEN
      CALL SEEK('C64')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBN
        READ(1,*,IOSTAT=ISO) (CBN_GL(L,1,M),M=MMIN,MMAX)
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBN_GL(L,1,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO

    endif
    Call Broadcast_Array(CBN_GL, master_id)

    !C65*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)
    NCARD='65'
    if( process_id == master_id )THEN
      CALL SEEK('C65')
      MMAX = 3 + NDYM + NTOX
      DO L=1,NCBN
        READ(1,*,IOSTAT=ISO) (CBN_GL(L,2,M),M=1,MMAX)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBN_GL(L,2,M),M=1,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(CBN_GL, master_id)

    !C66*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES
    !     SED(1 TO NSED),SND(1,NSND)
    NCARD='66'
    if( process_id == master_id )THEN
      CALL SEEK('C66')
      MMIN=MMAX+1
      MMAX=MMAX+NSED+NSND
      DO L=1,NCBN
        READ(1,*,IOSTAT=ISO) (CBN_GL(L,2,M),M=MMIN,MMAX)
        WRITE(7,1002)NCARD
        WRITE(7,*) (CBN_GL(L,2,M),M=MMIN,MMAX)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(CBN_GL, master_id)

  ENDIF    ! *** NCBN>0

  !C66A*  READ CONCENTRATION DATA ASSIMILATION PARAMETERS
  NCARD='66A'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C66A')
    READ(1,*,IOSTAT=ISO) NLCDA,TSCDA,(ISCDA(K),K=1,7)
    WRITE(7,1002)NCARD
    WRITE(7,*) NLCDA,TSCDA,(ISCDA(K),K=1,7)
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(NLCDA , master_id)
  Call Broadcast_Scalar(TSCDA , master_id)
  Call Broadcast_Array(ISCDA  , master_id)

  IF( NLCDA > 0 )THEN
    !C66B*  READ CONCENTRATION DATA ASSIMILATION LOCATIONS AND
    !      SERIES IDENTIFIERS
    NCARD='66B'
    if( process_id == master_id )THEN
      CALL SEEK('C66B')
      WRITE(7,1002)NCARD
      DO L=1,NLCDA
        READ(1,*,IOSTAT=ISO) ITPCDA(L),ICDA(L),JCDA(L),ICCDA(L),JCCDA(L),(NCSERA(L,K),K=1,7)
        WRITE(7,*)           ITPCDA(L),ICDA(L),JCDA(L),ICCDA(L),JCCDA(L),(NCSERA(L,K),K=1,7)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(ITPCDA, master_id)
    Call Broadcast_Array(ICDA, master_id)
    Call Broadcast_Array(JCDA, master_id)
    Call Broadcast_Array(ICCDA, master_id)
    Call Broadcast_Array(JCCDA, master_id)
    Call Broadcast_Array(NCSERA, master_id)

  ENDIF

  !C67*  READ NEUTRALLY BUOYANT PARTICLE DRIFTER DATA
  NCARD='67'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C67')
    READ(1,*,IOSTAT=ISO) ISPD
    
    NPD = 0
    WRITE(7,1002)NCARD
    WRITE(7,*) ISPD
    ! ***
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ISPD      , master_id)

  IF( NPD > 0 )THEN
    !C68*  READ NEUTRALLY BUOYANT PARTICLE INITIAL POSITIONS     DELME - These are not used in EFDC+  See mod_drifter.f90
    NCARD='68'
    if( process_id == master_id )THEN
      CALL SEEK('C68')
      !DO NP=1,NPD
      !  READ(1,*,IOSTAT=ISO) RI_GLOBAL, RJ_GLOBAL, RK(NP)
      !  WRITE(7,1002)NCARD
      !  WRITE(7,*) RI_GLOBAL,RJ_GLOBAL,RK(NP)
      !ENDDO
    endif


    !C69*  CONSTANTS FOR LONGITUDE AND LATITUDE OF CELL CENTERS (NOT USED)
    NCARD='69'
    ! *** ********************************************************
    if( process_id == master_id )THEN
      CALL SEEK('C69')
      !READ(1,*,IOSTAT=ISO) CDLON1,CDLON2,CDLON3,CDLAT1,CDLAT2,CDLAT3
      !
      !WRITE(7,1002)NCARD
      !WRITE(7,*) CDLON1,CDLON2,CDLON3,CDLAT1,CDLAT2,CDLAT3
      IF( ISO > 0 ) GOTO 100
    endif
  ENDIF

  !C70*  CONTROLS FOR WRITING ASCII OR BINARY DUMP FILES  (NOT USED)
  NCARD='70'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C70')
    READ(1,*,IOSTAT=ISO)ISDUMP,ISADMP,NSDUMP,TSDUMP,TEDUMP,ISDMPP,ISDMPU,ISDMPW,ISDMPT,IADJDMP

    WRITE(7,1002)NCARD
    WRITE(7,*)ISDUMP,ISADMP,NSDUMP,TSDUMP,TEDUMP,ISDMPP,ISDMPU,ISDMPW,ISDMPT,IADJDMP
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ISDUMP , master_id)
  Call Broadcast_Scalar(ISADMP , master_id)
  Call Broadcast_Scalar(NSDUMP , master_id)
  Call Broadcast_Scalar(TSDUMP , master_id)
  Call Broadcast_Scalar(TEDUMP , master_id)
  Call Broadcast_Scalar(ISDMPP , master_id)
  Call Broadcast_Scalar(ISDMPU , master_id)
  Call Broadcast_Scalar(ISDMPW , master_id)
  Call Broadcast_Scalar(ISDMPT , master_id)
  Call Broadcast_Scalar(IADJDMP, master_id)

  JSDUMP=1
  NCDUMP=1

  !C71*  CONTROLS FOR HORIZONTAL PLANE SCALAR FIELD CONTOURING
  NCARD='71'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C71')
    DO NS=1,7
      READ(1,*,IOSTAT=ISO) ISSPH(NS),NPSPH(NS),ISRSPH(NS),ISPHXY(NS)
      
      WRITE(7,1002)NCARD
      WRITE(7,*) ISSPH(NS),NPSPH(NS),ISRSPH(NS),ISPHXY(NS)
    ENDDO
    IF( ISO > 0 ) GOTO 100

  endif
  Call Broadcast_Array(ISSPH , master_id)
  Call Broadcast_Array(NPSPH , master_id)
  Call Broadcast_Array(ISRSPH, master_id)
  Call Broadcast_Array(ISPHXY, master_id)

  ! *** SET WATER COLUMN LINKAGE FLAG IF ANY CONSTITUENT OUTPUT IS ENABLED
  DO NS=1,7
    IF( ISTRAN(NS) >= 1 )ISSPH(8)=1
  ENDDO
  IF( ISSPH(4) >= 1 )ISSPH(8)=1

  !C71A*  CONTROLS FOR HORIZONTAL PLANE SEDIMENT BED PROPERTIES
  NCARD='71A'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C71A')
    READ(1,*,IOSTAT=ISO) ITMP,ISBEXP,NPBPH

    WRITE(7,1002)NCARD
    WRITE(7,*) ISBEXP,NPBPH
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(ISBEXP , master_id)
  Call Broadcast_Scalar(NPBPH  , master_id)

  JSBPH=1
  JSBPHA=1

  !C71B*  CONTROLS FOR FOOD CHAIN MODEL OUTPUT
  NCARD='71B'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C71B')
    READ(1,*,IOSTAT=ISO) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
    WRITE(7,1002)NCARD
    WRITE(7,*) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(ISFDCH , master_id)
  Call Broadcast_Scalar(NFDCHZ , master_id)
  Call Broadcast_Scalar(HBFDCH , master_id)
  Call Broadcast_Scalar(TFCAVG , master_id)

  !C72*  CONTROLS FOR EFDC_EXPLORER LINKAGE AND SURFACE ELEVATION RESIDUAL OUTPUT
  JSFDCH=1
  NCARD='72'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C72')
    READ(1,*,IOSTAT=ISO) ISPPH,NPPPH,ISRPPH,IPPHXY

    WRITE(7,1002)NCARD
    WRITE(7,*) ISPPH,NPPPH,ISRPPH,IPPHXY
    IF( ISO > 0 ) GOTO 100
    IF( ISPPH < 0 ) ISPPH = 1
  endif
  Call Broadcast_Scalar(ISPPH  , master_id)
  Call Broadcast_Scalar(NPPPH  , master_id)
  Call Broadcast_Scalar(ISRPPH , master_id)
  Call Broadcast_Scalar(IPPHXY , master_id)

  !C73*  CONTROLS FOR HORIZONTAL PLANE VELOCITY PLOTTING
  NCARD='73'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C73')
    READ(1,*,IOSTAT=ISO) ISVPH,NPVPH,ISRVPH,IVPHXY

    WRITE(7,1002)NCARD
    WRITE(7,*) ISVPH,NPVPH,ISRVPH,IVPHXY
    IF( ISO > 0 ) GOTO 100
  endif

  Call Broadcast_Scalar(ISVPH  , master_id)
  Call Broadcast_Scalar(NPVPH  , master_id)
  Call Broadcast_Scalar(ISRVPH , master_id)
  Call Broadcast_Scalar(IVPHXY , master_id)

  !C74*  CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING
  NCARD='74'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C74')
    READ(1,*,IOSTAT=ISO) ISECSPV,NPSPV(1),ISSPV(1),ISRSPV(1),ISHPLTV(1)

    WRITE(7,1002)NCARD
    WRITE(7,*) ISECSPV,NPSPV(1),ISSPV(1),ISRSPV(1),ISHPLTV(1)
    SHPLTV(1)=FLOAT(ISHPLTV(1))
    SBPLTV(1)=1.0-SHPLTV(1)
    DO NS=2,7
      READ(1,*,IOSTAT=ISO) IDUMMY,NPSPV(NS),ISSPV(NS),ISRSPV(NS),ISHPLTV(NS)

      WRITE(7,1002)NCARD
      WRITE(7,*) IDUMMY,NPSPV(NS),ISSPV(NS),ISRSPV(NS),ISHPLTV(NS)
      SHPLTV(NS)=FLOAT(ISHPLTV(NS))
      SBPLTV(NS)=1.0-SHPLTV(NS)
    ENDDO
    IF( ISO > 0 ) GOTO 100
  endif


  Call Broadcast_Array(NPSPV  , master_id)
  Call Broadcast_Array(ISSPV  , master_id)
  Call Broadcast_Array(ISRSPV , master_id)
  Call Broadcast_Array(ISHPLTV, master_id)

  IF( ISECSPV > 0 )THEN
    !C75*  MORE CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING
    NCARD='75'
    if( process_id == master_id )THEN
      CALL SEEK('C75')
      DO IS=1,ISECSPV
        READ(1,*,IOSTAT=ISO) DUM,NIJSPV(IS),CCTITLE(10+IS)

        WRITE(7,1002)NCARD
        WRITE(7,*) DUM,NIJSPV(IS),CCTITLE(10+IS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(NIJSPV, master_id)

    !C76*  I,J LOCATIONS DEFINING VERTICAL PLANE FOR CONTOURING
    NCARD='76'
    if( process_id == master_id )THEN
      CALL SEEK('C76')
      DO IS=1,ISECSPV
        DO NPP=1,NIJSPV(IS)
          READ(1,*,IOSTAT=ISO) DUM,ISPV(NPP,IS),JSPV(NPP,IS)

          WRITE(7,1002)NCARD
          WRITE(7,*) DUM,ISPV(NPP,IS),JSPV(NPP,IS)
          IF( ISO > 0 ) GOTO 100
        ENDDO
      ENDDO
    endif
    Call Broadcast_Array(ISPV, master_id)
    Call Broadcast_Array(JSPV, master_id)

  ENDIF

  !C77*  CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
  NCARD='77'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C77')
    READ(1,*,IOSTAT=ISO) ISECVPV,NPVPV,ISVPV,ISRVPV
    
    WRITE(7,1002)NCARD
    WRITE(7,*) ISECVPV,NPVPV,ISVPV,ISRVPV
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(ISECVPV, master_id)
  Call Broadcast_Scalar(NPVPV  , master_id)
  Call Broadcast_Scalar(ISVPV  , master_id)
  Call Broadcast_Scalar(ISRVPV , master_id)

  IF( ISECVPV > 0 )THEN
    !C78*   MORE CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
    NCARD='78'
    if( process_id == master_id )THEN
      CALL SEEK('C78')
      DO IS=1,ISECVPV
        READ(1,*,IOSTAT=ISO) DUM,NIJVPV(IS),ANGVPV(IS),CVTITLE(10+IS)
        WRITE(7,1002)NCARD
        WRITE(7,*) DUM,NIJVPV(IS),ANGVPV(IS),CVTITLE(10+IS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif


    Call Broadcast_Array(NIJVPV , master_id)
    Call Broadcast_Array(ANGVPV , master_id)

    !C79*   MORE CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
    NCARD='79'
    if( process_id == master_id )THEN
      CALL SEEK('C79')
      DO IS=1,ISECVPV
        DO NPP=1,NIJVPV(IS)
          READ(1,*,IOSTAT=ISO) DUM,IVPV(NPP,IS),JVPV(NPP,IS)
          WRITE(7,1002)NCARD
          WRITE(7,*) DUM,IVPV(NPP,IS),JVPV(NPP,IS)
          IF( ISO > 0 ) GOTO 100
        ENDDO
      ENDDO
    endif
    Call Broadcast_Array(IVPV, master_id)
    Call Broadcast_Array(JVPV, master_id)

  ENDIF

  !C80*  CONTROLS FOR 3D FIELD OUTPUT
  NCARD='80'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C80')
    READ(1,*,IOSTAT=ISO)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
    
    WRITE(7,1002)NCARD
    WRITE(7,*)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(IS3DO  , master_id)
  Call Broadcast_Scalar(ISR3DO , master_id)
  Call Broadcast_Scalar(NP3DO  , master_id)
  Call Broadcast_Scalar(KPC    , master_id)
  Call Broadcast_Scalar(NWGG   , master_id)
  Call Broadcast_Scalar(I3DMIN , master_id)
  Call Broadcast_Scalar(I3DMAX , master_id)
  Call Broadcast_Scalar(J3DMIN , master_id)
  Call Broadcast_Scalar(J3DMAX , master_id)
  Call Broadcast_Scalar(I3DRW  , master_id)
  Call Broadcast_Scalar(SELVMAX, master_id)
  Call Broadcast_Scalar(BELVMIN, master_id)

  NCALL3D=0
  NRCAL3D=0

  !C81* OUTPUT ACTIVATION AND SCALES FOR 3D FIELD OUTPUT
  NCARD='81'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C81')
    READ(1,*,IOSTAT=ISO)CDUM,IS3DUUU,JS3DUUU,UUU3DMA,UUU3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DUUU,JS3DUUU,UUU3DMA,UUU3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DVVV,JS3DVVV,VVV3DMA,VVV3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DVVV,JS3DVVV,VVV3DMA,VVV3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DWWW,JS3DWWW,WWW3DMA,WWW3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DWWW,JS3DWWW,WWW3DMA,WWW3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DSAL,JS3DSAL,SAL3DMA,SAL3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DSAL,JS3DSAL,SAL3DMA,SAL3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DTEM,JS3DTEM,TEM3DMA,TEM3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DTEM,JS3DTEM,TEM3DMA,TEM3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DDYE,JS3DDYE,DYE3DMA,DYE3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DDYE,JS3DDYE,DYE3DMA,DYE3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DSED,JS3DSED,SED3DMA,SED3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DSED,JS3DSED,SED3DMA,SED3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DSND,JS3DSND,SND3DMA,SND3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DSND,JS3DSND,SND3DMA,SND3DMI
    IF( ISO > 0 ) GOTO 100
    READ(1,*,IOSTAT=ISO)CDUM,IS3DTOX,JS3DTOX,TOX3DMA,TOX3DMI
    WRITE(7,1002)NCARD
    WRITE(7,*)CDUM,IS3DTOX,JS3DTOX,TOX3DMA,TOX3DMI
    IF( ISO > 0 ) GOTO 100

  endif


  Call Broadcast_Scalar(IS3DUUU, master_id)
  Call Broadcast_Scalar(JS3DUUU, master_id)
  Call Broadcast_Scalar(UUU3DMA, master_id)
  Call Broadcast_Scalar(UUU3DMI, master_id)
  Call Broadcast_Scalar(IS3DVVV, master_id)
  Call Broadcast_Scalar(JS3DVVV, master_id)
  Call Broadcast_Scalar(VVV3DMA, master_id)
  Call Broadcast_Scalar(VVV3DMI, master_id)
  Call Broadcast_Scalar(IS3DWWW, master_id)
  Call Broadcast_Scalar(JS3DWWW, master_id)
  Call Broadcast_Scalar(WWW3DMA, master_id)
  Call Broadcast_Scalar(WWW3DMI, master_id)
  Call Broadcast_Scalar(IS3DSAL, master_id)
  Call Broadcast_Scalar(JS3DSAL, master_id)
  Call Broadcast_Scalar(SAL3DMA, master_id)
  Call Broadcast_Scalar(SAL3DMI, master_id)
  Call Broadcast_Scalar(IS3DTEM, master_id)
  Call Broadcast_Scalar(JS3DTEM, master_id)
  Call Broadcast_Scalar(TEM3DMA, master_id)
  Call Broadcast_Scalar(TEM3DMI, master_id)
  Call Broadcast_Scalar(IS3DDYE, master_id)
  Call Broadcast_Scalar(JS3DDYE, master_id)
  Call Broadcast_Scalar(DYE3DMA, master_id)
  Call Broadcast_Scalar(DYE3DMI, master_id)
  Call Broadcast_Scalar(IS3DSED, master_id)
  Call Broadcast_Scalar(JS3DSED, master_id)
  Call Broadcast_Scalar(SED3DMA, master_id)
  Call Broadcast_Scalar(SED3DMI, master_id)
  Call Broadcast_Scalar(IS3DSND, master_id)
  Call Broadcast_Scalar(JS3DSND, master_id)
  Call Broadcast_Scalar(SND3DMA, master_id)
  Call Broadcast_Scalar(SND3DMI, master_id)
  Call Broadcast_Scalar(IS3DTOX, master_id)
  Call Broadcast_Scalar(JS3DTOX, master_id)
  Call Broadcast_Scalar(TOX3DMA, master_id)
  Call Broadcast_Scalar(TOX3DMI, master_id)



  IF( ISTRAN(1) < 1 )THEN
    IS3DSAL = 0
    JS3DSAL = 0
  ENDIF
  IF( ISTRAN(2) < 1 )THEN
    IS3DTEM = 0
    JS3DTEM = 0
  ENDIF
  IF( ISTRAN(3) < 1 )THEN
    IS3DDYE = 0
    JS3DDYE = 0
  ENDIF
  IF( ISTRAN(5) < 1 )THEN
    IS3DTOX = 0
    JS3DTOX = 0
  ENDIF
  IF( ISTRAN(6) < 1 )THEN
    IS3DSED = 0
    JS3DSED = 0
  ENDIF
  IF( ISTRAN(7) < 1 )THEN
    IS3DSND = 0
    JS3DSND = 0
  ENDIF

  !C82* INPLACE HARMONIC ANALYSIS PARAMETERS
  NCARD='82'
  ! *** ********************************************************
  if( process_id == master_id )THEN

    CALL SEEK('C82')
    READ(1,*,IOSTAT=ISO) ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
    
    WRITE(7,1002)NCARD
    WRITE(7,*) ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(ISLSHA  , master_id)
  Call Broadcast_Scalar(MLLSHA  , master_id)
  Call Broadcast_Scalar(NTCLSHA , master_id)
  Call Broadcast_Scalar(ISLSTR  , master_id)
  Call Broadcast_Scalar(ISHTA   , master_id)

  IF( MLLSHA > 0 )THEN
    !C83* HARMONIC ANALYSIS LOCATIONS AND SWITCHES
    NCARD='83'
    if( process_id == master_id )THEN
      CALL SEEK('C83')
      DO M=1,MLLSHA
        READ(1,*,IOSTAT=ISO) ILLSHA(M),JLLSHA(M),LSHAP(M),LSHAB(M),LSHAUE(M),LSHAU(M),CLSL(M)
        
        WRITE(7,1002)NCARD
        WRITE(7,*) ILLSHA(M),JLLSHA(M),LSHAP(M),LSHAB(M),LSHAUE(M),LSHAU(M),CLSL(M)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(ILLSHA , master_id)
    Call Broadcast_Array(JLLSHA , master_id)
    Call Broadcast_Array(LSHAP  , master_id)
    Call Broadcast_Array(LSHAB  , master_id)
    Call Broadcast_Array(LSHAUE , master_id)
    Call Broadcast_Array(LSHAU  , master_id)
  ENDIF

  !C84* CONTROLS FOR WRITING TO TIME SERIES FILES
  NCARD='84'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C84')
    READ(1,*,IOSTAT=ISO)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR

    WRITE(7,1002)NCARD
    WRITE(7,*)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(ISTMSR   , master_id)
  Call Broadcast_Scalar(MLTMSR   , master_id)
  Call Broadcast_Scalar(NBTMSR   , master_id)
  Call Broadcast_Scalar(NSTMSR   , master_id)
  Call Broadcast_Scalar(NWTMSR   , master_id)
  Call Broadcast_Scalar(NTSSTSP  , master_id)
  Call Broadcast_Scalar(TCTMSR   , master_id)

  JSTMSR=1
  NCTMSR=1
  JSHYDOUT=1
  NCHYDOUT=1

  IF( NTSSTSP > 0 )THEN
    NCARD='85'
    if( process_id == master_id )THEN
      CALL SEEK('C85')
      DO ITSSS=1,NTSSTSP
        READ(1,*,IOSTAT=ISO)IDUM,MTSSTSP(ITSSS)
        WRITE(7,1002)NCARD
        WRITE(7,*)IDUM,MTSSTSP(ITSSS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif

    Call Broadcast_Array(MTSSTSP, master_id)

    NCARD='86'
    if( process_id == master_id )THEN
      CALL SEEK('C86')
      DO ITSSS=1,NTSSTSP
        DO MTSSS=1,MTSSTSP(ITSSS)
          READ(1,*,IOSTAT=ISO)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),TSSTOP(MTSSS,ITSSS)

          WRITE(7,1002)NCARD
          WRITE(7,*)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),TSSTOP(MTSSS,ITSSS)
          IF( ISO > 0 ) GOTO 100
        ENDDO
      ENDDO
    endif
    Call Broadcast_Array( TSSTRT, master_id)
    Call Broadcast_Array( TSSTOP, master_id)

  ENDIF

  IF( MLTMSR > 0 )THEN
    NCARD='87'
    if( process_id == master_id )THEN
      CALL SEEK('C87')
      DO M=1,MLTMSR
        READ(1,*,IOSTAT=ISO)ILTMSR_GL(M),JLTMSR_GL(M),NTSSSS_GL(M),MTMSRP_GL(M),MTMSRC_GL(M),MTMSRA_GL(M),&
          MTMSRUE_GL(M),MTMSRUT_GL(M),MTMSRU_GL(M),MTMSRQE_GL(M),MTMSRQ_GL(M),CLTMSR_GL(M)
        
        WRITE(7,1002)NCARD
        WRITE(7,*)ILTMSR_GL(M),JLTMSR_GL(M),NTSSSS_GL(M),MTMSRP_GL(M),MTMSRC_GL(M),MTMSRA_GL(M),&
          MTMSRUE_GL(M),MTMSRUT_GL(M),MTMSRU_GL(M),MTMSRQE_GL(M),MTMSRQ_GL(M),CLTMSR_GL(M)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    endif
    Call Broadcast_Array(ILTMSR_GL, master_id)
    Call Broadcast_Array(JLTMSR_GL, master_id)
    Call Broadcast_Array(NTSSSS_GL, master_id)
    Call Broadcast_Array(MTMSRP_GL, master_id)
    Call Broadcast_Array(MTMSRC_GL, master_id)
    Call Broadcast_Array(MTMSRA_GL, master_id)
    Call Broadcast_Array(MTMSRUE_GL, master_id)
    Call Broadcast_Array(MTMSRUT_GL, master_id)
    Call Broadcast_Array(MTMSRU_GL , master_id)
    Call Broadcast_Array(MTMSRQE_GL, master_id)
    Call Broadcast_Array(MTMSRQ_GL , master_id)
    !Call Broadcast_Array(CLTMSR , master_id)  ! *** Character    DELME - TODO
  ENDIF

  !C88 ** HIGH FREQUENCY OUTPUT
  NCARD='88'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C88')
    READ(1,*,IOSTAT=ISO) HFREOUT
    
    WRITE(7,1002)NCARD
    WRITE(7,*) HFREOUT
    IF( ISO > 0 ) GOTO 100
  endif
  Call Broadcast_Scalar(HFREOUT, master_id)

  !NCARD='88'
  !CALL SEEK('C88')
  !READ(1,*,IOSTAT=ISO)ISVSFP,MDVSFP,MLVSFP,TMVSFP,TAVSFP
  !WRITE(7,1002)NCARD
  !WRITE(7,*)ISVSFP,MDVSFP,MLVSFP,TMVSFP,TAVSFP
  !IF(ISO > 0 ) GOTO 100
  JSVSFP=1

#ifdef NCOUT   
  ! *** NETCDF GENERATION CONTROLS
  NCARD='91'
  ! *** ********************************************************
  if( process_id == master_id )THEN
    CALL SEEK('C91')
    READ(1,'(A)',IOSTAT=ISO) STR
    !                                       NOT      NOT
    READ(STR,*,ERR=100) NCDFOUT,DEFLEV,ROTA,BLK,UTMZ,HREST,BASEDATE,BASETIME,PROJ
    IF( UTMZ < 0 )THEN
      UTMZ = ABS(UTMZ)
      HEMI = 2
    ELSE
      HEMI = 1
    ENDIF

    WRITE(7,1002)NCARD
    WRITE(7,*) NCDFOUT,DEFLEV,ROTA,BLK,HEMI,UTMZ,HREST,BASEDATE,BASETIME,PROJ
  endif

  Call Broadcast_Scalar(NCDFOUT   , master_id)
  Call Broadcast_Scalar(DEFLEV    , master_id)
  Call Broadcast_Scalar(ROTA      , master_id)
  Call Broadcast_Scalar(BLK       , master_id)
  Call Broadcast_Scalar(HEMI      , master_id)
  Call Broadcast_Scalar(UTMZ      , master_id)
  Call Broadcast_Scalar(HREST     , master_id)
  Call Broadcast_Scalar(BASEDATE  , master_id)
  Call Broadcast_Scalar(BASETIME  , master_id)
  Call Broadcast_Scalar(PROJ      , master_id)

  IF( NCDFOUT > 0 )THEN
    NCARD='91A'
    if( process_id == master_id )THEN
      CALL SEEK('C91A')
      READ(1,*,ERR=100) ISSGLFIL,TBEGNCDF,TENDNCDF

      NCARD='91B'
      CALL SEEK('C91B')
      READ(1,*,ERR=100) ISNCDF
      DO NS=1,8
        IF( ISTRAN(NS) == 0 ) ISNCDF(NS) = 0
      ENDDO

      IF(ISPD < 1 )   ISNCDF(9)  = 0
      IF(NWSER == 0 .AND. WIND.IFLAG == 0 .AND. ICYCLONE == 0) ISNCDF(11) = 0
      IF(ISWAVE == 0 ) ISNCDF(12) = 0

      WRITE(7,1002)NCARD
      WRITE(7,*) ISNCDF(1:12)

    endif

    Call Broadcast_Scalar(ISSGLFIL, master_id)
    Call Broadcast_Scalar(TBEGNCDF, master_id)
    Call Broadcast_Scalar(TENDNCDF, master_id)
    
    Call Broadcast_Array(ISNCDF, master_id)
    
  ENDIF
#endif

  NTS=INT8(NTC)*NTSPTC
  NBVSFP=NTC*NTSPTC
  NSVSFP=0

  DT = TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)

  GOTO 2000


  ! *** ********************************************************
  if( process_id == master_id )THEN

    ! *** WRITE INPUT ERROR MESSAGES AND TERMINATE RUN
100 WRITE(6,1001)NCARD
    WRITE(8,1001)NCARD
    WRITE(7,1001)NCARD

    CALL STOPP('')
  endif

2000 CONTINUE

  ! *** ********************************************************
  if( process_id == master_id )THEN
    ! *** NOW REWIND UNIT 1 & READ IN AS CHARACTER TO WRITE TO UNIT 7
    REWIND (1)
21  READ(1,22,END=24) TEXT
    WRITE (7,23) TEXT
    GOTO 21
24  CONTINUE
    CLOSE(1)
22  FORMAT (A80)
23  FORMAT (1X,A80)
  endif

  ! *** READ CELL TYPES FROM FILES CELL.INP
  if( process_id == master_id )THEN
    WRITE(*,'(A)')'READING CELL.INP'
    OPEN(1,FILE='cell.inp',STATUS='UNKNOWN')

    ! *** SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT
    STRC=READSTR(1)
    READ(STRC,*)JCTMP

    IF( JCTMP /= JC_GLOBAL )THEN

      ! ***   READ OLD FILE FORMAT
      JACROSS = JC_Global
      IF(JC_Global > 640) JACROSS = 640
      DO JT=1,JC_GLOBAL,JACROSS
        JF=JT
        JLAST = JT + JACROSS-1
        IF( JLAST > JC_GLOBAL) JLAST = JC_GLOBAL
        DO I=1,IC_GLOBAL 
          READ(1,6,IOSTAT=ISO) (IJCT_GLOBAL(I,J),J=JF,JLAST)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
        ENDDO
      ENDDO

    ELSE

      IF( IC_Global > 640 )THEN
        IACROSS = 640
        DO IT=1,IC_Global,IACROSS
          IFIRST = IT
          ILAST = IT + IACROSS - 1
          IF( ILAST>IC_Global) ILAST = IC_GLOBAL
        
          DO J=JC_GLOBAL,1,-1
            READ(1,66,IOSTAT=ISO) ADUMMY, (IJCT_GLOBAL(I,J),I=IFIRST,ILAST)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
          ENDDO
        ENDDO
      ELSE
        IFIRST = 1
        ILAST = IC_GLOBAL
        DO J=JC_GLOBAL,1,-1
          READ(1,66,IOSTAT=ISO) ADUMMY, (IJCT_GLOBAL(I,J),I=IFIRST,ILAST)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
        ENDDO
      ENDIF
    ENDIF
    CLOSE(1)
  ENDIF
  Call Broadcast_Array(IJCT_GLOBAL, master_id)
  
8 FORMAT ('   CELL TYPE ARRAY,J=',I5,2X,'TO J=',I5,//)

  !----------------------------------------------------------------------C
  ! *** ********************************************************
  if( process_id == master_id )THEN

    WRITE(*,'(A)')'READING CELLLT.INP'
    OPEN(1,FILE='celllt.inp',STATUS='UNKNOWN')

    ! *** SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT
    STRC=READSTR(1)
    READ(STRC,*)JCTMP

    IF( JCTMP /= JC_GLOBAL )THEN
      ! ***   READ OLD FILE FORMAT
      JACROSS = JC_Global
      IF( JC_Global > 640 ) JACROSS = 640
      DO JT=1,JC_Global,JACROSS
        JF = JT
        JLAST = JT + JACROSS - 1
        IF( JLAST>JC_Global) JLAST = JC_Global
        WRITE (7,8) JF, JLAST
        DO I=1,IC_Global
          READ(1,6,IOSTAT=ISO) (IJCTLT_Global(I,J),J=JF,JLAST)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELLLT.INP')
          WRITE (7,16) (IJCTLT_Global(I,J),J=JF,JLAST)
        ENDDO
        WRITE(7,15)
      ENDDO

    ELSE
      ! ***   READ NEW FILE FORMAT
      IF( IC_GLOBAL > 640 )THEN
        IACROSS = 640
        DO IT=1,IC_GLOBAL,IACROSS
          IFIRST = IT
          ILAST = IT+IACROSS-1
          IF( ILAST>IC_GLOBAL) ILAST = IC_GLOBAL
          WRITE (7,'(2I5)') IFIRST, ILAST
          DO J=JC_GLOBAL,1,-1
            READ(1,66,IOSTAT=ISO) ADUMMY, (IJCTLT_GLOBAL(I,J),I=IFIRST,ILAST)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELLLT.INP')
            WRITE (7,166) ADUMMY, (IJCTLT_GLOBAL(I,J),I=IFIRST,ILAST)
          ENDDO
          WRITE(7,15)
        ENDDO
      ELSE
        IFIRST = 1
        ILAST = IC_GLOBAL
        WRITE (7,88) IFIRST, ILAST
        DO J=JC_GLOBAL,1,-1
          READ(1,66,IOSTAT=ISO) ADUMMY, (IJCTLT_GLOBAL(I,J),I=IFIRST,ILAST)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELLLT.INP')
          WRITE (7,166) ADUMMY, (IJCTLT_GLOBAL(I,J),I=IFIRST,ILAST)
        ENDDO
        WRITE(7,15)
      ENDIF
    ENDIF
    CLOSE(1)

88  FORMAT ('   CELLLT TYPE ARRAY,I=',I5,' TO I=',I5,//)
  endif
  
  Call Broadcast_Array(IJCTLT_GLOBAL, master_id)

  ! *** ********************************************************
  ! ***  GENERATE CELL MAPPINGS
  ! *** If running and MPI calculation call childgrid
#ifdef _MPI
  Call Child_Grid

#else
  DO I = 1, IC
    IL2IG(I) = I
  END DO
  DO J = 1, JC
    JL2JG(J) = J
  END DO
#endif

  ! *** This is a modified cell map
  Call CELLMAP

  Call Write_Mapping_3Dgraphics !***Moved from end of CELLMAP to here

  ! *** FORMAT STATEMENTS FOR EFDC.INP
 15 FORMAT (/)
  6 FORMAT (640I1)
 16 FORMAT (1X,640I1)
 66 FORMAT (A5,10000I1)
166 FORMAT (1X,A5,10000I1)

  ! *** READ CURVILINEAR-ORTHOGONAL OR VARIABLE CELL DATA FROM FILE DXDY.INP

  ! *** INITIALIZE CELL DIMENSIONS TO CONSTANT CARTESIAN OR DUMMY VALUES
  DO L=1,LC
    DXP(L) = DX*DXYCVT
    DYP(L) = DY*DXYCVT
    ZBR(L) = ZBRADJ
  ENDDO

  ! *** READ IN DX, DY, DEPTH AND BOTTOM ELEVATION AT CELL CENTERS OF VARIABLE CELLS
  LMHK=.FALSE.
  if( process_id == master_id )THEN
    WRITE(*,'(A)')'READING DXDY.INP'
    OPEN(1,FILE='dxdy.inp',STATUS='UNKNOWN')
    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR=READSTR(1)
  end if

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
  
  IF( ISVEG /= 1 )THEN
    ! *** READ DXDY WITHOUT VEGETATION
    DO LL=2,LA_GLOBAL
      if( process_id == master_id )THEN
        READ(1,*,IOSTAT=ISO) IGL, JGL, DXIJ, DYIJ, HIJ, BELVIJ, ZBRIJ
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DXDY.INP')
        
        ! *** get global versions as they are needed for the NETCDF writing out
        IMN_Global = MIN(IMN_Global, IGL) 
        JMN_Global = MIN(JMN_Global, JGL) 
        IMX_Global = MAX(IMX_Global, IGL) 
        JMX_Global = MAX(JMX_Global, JGL) 
        
      end if

      ! *** Added MPI Domain decomp
      Call Broadcast_Scalar(IGL,   master_id)
      Call Broadcast_Scalar(JGL,   master_id)
      Call Broadcast_Scalar(DXIJ,  master_id)
      Call Broadcast_Scalar(DYIJ,  master_id)
      Call Broadcast_Scalar(HIJ,   master_id)
      Call Broadcast_Scalar(BELVIJ,master_id)
      Call Broadcast_Scalar(ZBRIJ, master_id)

      LG = LIJ_Global(IGL, JGL)
      DXP_Global(LG) = DXIJ
      DYP_Global(LG) = DYIJ
      HP_Global(LG) = HADADJ + HCVRT*HIJ
      HP_Global(LG) = MAX(HP_Global(LG), HMIN)

      ! *** Get local i,j
      I = IG2IL(igl)
      J = JG2JL(jgl)

      IF( I > 0 .AND. I <= IC )THEN
        IF( J > 0 .AND. J <= JC )THEN
          L = LIJ(I,J)
          IF( L > 0 )THEN
          
            IMN = MIN(IMN,I)
            JMN = MIN(JMN,J)
            IMX = MAX(IMX,I)
            JMX = MAX(JMX,J)
          
            DXP(L) = DXYCVT*DXIJ
            DYP(L) = DXYCVT*DYIJ
            HMP(L) = HADADJ + HCVRT*HIJ
            HMP(L) = MAX(HMP(L), HMIN)

            BELV(L)  = BELADJ + BELCVRT*BELVIJ
            BELV1(L) = BELADJ + BELCVRT*BELVIJ
            ZBR(L)   = ZBRADJ + ZBRCVRT*ZBRIJ
          End if
        End if
      End if
    ENDDO
      
  ELSE
    ! *** READ DXDY WITH VEGETATION
    DO LL=2,LA_GLOBAL
      if( process_id == master_id )THEN
        READ(1,*,IOSTAT=ISO) IGL, JGL, DXIJ, DYIJ, HIJ, BELVIJ, ZBRIJ, MVEGIJT  ! !SCJ adding MHK devices when MVEGIJT>90
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DXDY.INP')
        
        ! *** get global versions as they are needed for the NETCDF writing out
        IMN_Global = MIN(IMN_Global, IGL) 
        JMN_Global = MIN(JMN_Global, JGL) 
        IMX_Global = MAX(IMX_Global, IGL) 
        JMX_Global = MAX(JMX_Global, JGL) 
        
      endif

      Call Broadcast_Scalar(IGL,     master_id)
      Call Broadcast_Scalar(JGL,     master_id)
      Call Broadcast_Scalar(DXIJ,    master_id)
      Call Broadcast_Scalar(DYIJ,    master_id)
      Call Broadcast_Scalar(HIJ,     master_id)
      Call Broadcast_Scalar(BELVIJ,  master_id)
      Call Broadcast_Scalar(ZBRIJ,   master_id)
      Call Broadcast_Scalar(MVEGIJT, master_id)

      ! *** Make sure the master gets HMP global exactly from input file
      LG = LIJ_Global(IGL,JGL)
      DXP_Global(LG) = DXIJ
      DYP_Global(LG) = DYIJ
      HP_Global(LG) = HADADJ + HCVRT*HIJ
      HP_Global(LG) = MAX(HP_Global(LG), HMIN)
        
      ! *** Get local i,j
      I = IG2IL(igl)
      J = JG2JL(jgl)
        
      IF( I > 0 .AND. I <= IC )THEN
        IF( J > 0 .AND. J <= JC )THEN
          L = LIJ(I,J)
          IF( L > 0 )THEN
            L=LIJ(I,J)
            IMN = MIN(IMN,I)
            JMN = MIN(JMN,J)
            IMX = MAX(IMX,I)
            JMX = MAX(JMX,J)
            DXP(L) = DXYCVT*DXIJ
            DYP(L) = DXYCVT*DYIJ
            HMP(L) = HADADJ + HCVRT*HIJ
            HMP(L) = MAX(HMP(L),HMIN)
        
            BELV(L)  = BELADJ + BELCVRT*BELVIJ
            BELV1(L) = BELADJ + BELCVRT*BELVIJ
            ZBR(L)   = ZBRADJ + ZBRCVRT*ZBRIJ
            MVEGL(L) = MVEGIJT
              
            ! ***  MHK
            IF( MVEGIJT > 90 )THEN
              LMHK = .TRUE.
              ITURB = ITURB+1
              IJLTURB(ITURB,1) = I
              IJLTURB(ITURB,2) = J
              IJLTURB(ITURB,3) = LIJ(I,J)
            ENDIF
          End if
        End if
      End if
    ENDDO
  ENDIF
    
  if( process_id == master_id )THEN
    CLOSE(1)
  end if

  ! *** OPEN FILE MODDXDY.INP TO MODIFY INPUT VALUES OF DX AND DY
  IF( IMDXDY > 0 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING MODDXDY.INP'
      OPEN(1,FILE='moddxdy.inp',STATUS='UNKNOWN')
      ! *** SKIP OVER TITLE AND HEADER LINES
      STR=READSTR(1)
      READ(1,*) NMDXDY
    Endif

    Call Broadcast_Scalar(NMDXDY, master_id)

    IF( NMDXDY >= 1 )THEN
      DO NMD=1, NMDXDY
        if( process_id == master_id )THEN
          READ(1,*) ITMP, JTMP, RMDX, RMDY
        Endif

        Call Broadcast_Scalar(ITMP, master_id)
        Call Broadcast_Scalar(JTMP, master_id)
        Call Broadcast_Scalar(RMDX, master_id)
        Call Broadcast_Scalar(RMDY, master_id)

        ! *** Get local i,j
        I = IG2IL(ITMP)
        J = JG2JL(JTMP)
        ! *** Only select cells in the current domain
        IF( I > 0 .AND. I <= IC )THEN
          IF( J > 0 .AND. J <= JC )THEN
            LTMP = LIJ(ITMP, JTMP) ! Map to local L value

            DXP(LTMP) = RMDX*DXP(LTMP)
            DYP(LTMP) = RMDY*DYP(LTMP)
          End if
        End if

      ENDDO
    ENDIF

    ! *** Only have the master close the file
    if( process_id == master_id )THEN
      CLOSE(1)
    endif

  ENDIF

  ! *** OPEN FILE MODCHAN.INP TO INSERT SUBGRID CHANNELS INTO HOST CELLS
  MDCHH=0
  IF( ISCHAN > 0 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING MODCHAN.INP'
      OPEN(1,FILE='modchan.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      IF( ISCHAN == 1 )THEN
        READ(1,*) MDCHH,MDCHHD,MDCHHD2
        READ(1,*) MDCHITM,MDCHHQ,QCHERR
        IF( MDCHH >= 1 )THEN
          DO NMD=1,MDCHH
            READ(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD)
            QCHANU(NMD)=0.
            QCHANUN(NMD)=0.
            QCHANV(NMD)=0.
            QCHANVN(NMD)=0.
          ENDDO
        ENDIF
      ENDIF
      IF( ISCHAN == 2 )THEN
        READ(1,*) MDCHH,MDCHHD,MDCHHD2
        READ(1,*) MDCHITM,MDCHHQ,QCHERR
        IF( MDCHH >= 1 )THEN
          DO NMD=1,MDCHH
            READ(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD),CHANLEN(NMD),PMDCH(NMD)
            QCHANU(NMD)=0.
            QCHANUN(NMD)=0.
            QCHANV(NMD)=0.
            QCHANVN(NMD)=0.
          ENDDO
        ENDIF
      ENDIF
      CLOSE(1)
      IF( MDCHH >= 1 )THEN
        DO NMD=1,MDCHH
          LMDCHH(NMD)=LIJ(IMDCHH(NMD),JMDCHH(NMD))
          IF( IMDCHU(NMD) == 1 .AND. JMDCHU(NMD) == 1 )THEN
            LMDCHU(NMD)=1
          ELSE
            LMDCHU(NMD)=LIJ(IMDCHU(NMD),JMDCHU(NMD))
          ENDIF
          IF( IMDCHV(NMD) == 1 .AND. JMDCHV(NMD) == 1 )THEN
            LMDCHV(NMD)=1
          ELSE
            LMDCHV(NMD)=LIJ(IMDCHV(NMD),JMDCHV(NMD))
          ENDIF
        ENDDO
      ENDIF

    endif

    Call Broadcast_Array(MDCHTYP , master_id)
    Call Broadcast_Array(IMDCHH  , master_id)
    Call Broadcast_Array(JMDCHH  , master_id)
    Call Broadcast_Array(IMDCHU  , master_id)
    Call Broadcast_Array(JMDCHU  , master_id)
    Call Broadcast_Array(IMDCHV  , master_id)
    Call Broadcast_Array(JMDCHV  , master_id)
    Call Broadcast_Array(CHANLEN , master_id)
    Call Broadcast_Array(PMDCH   , master_id)

  ENDIF

  ! *** ENABLE SPATIALLY VARIABLE BACKGROUND AHO
  IF( AHO < 0. )THEN
    AHMAX = ABS(AHO)
    AHMIN = 1.0E32
    DO L=2,LA
      AHOXY(L) = ABS(AHO)*DXP(L)*DYP(L)
      AHMAX    = MAX(AHOXY(L),AHMAX)
      AHMIN    = MIN(AHOXY(L),AHMIN)
    ENDDO

    if( process_id == master_id )then
      PRINT '(A,2E16.5)','VARIABLE AHO USED (MIN,MAX): ', AHMIN, AHMAX
    endif

  ELSE
    AHOXY = AHO
  ENDIF

  ! *** SMAGORINSKY AND BACKGROUND DIFFUSIVITY
  ! *** Constant and/or default values
  AHO = ABS(AHO)         

  ! *** ENABLE SPATIALLY VARIABLE SMAGORINSKY AND BACKGROUND DIFFUSIVITY
  IF( AHD < 0. )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, 2, 0.0)
    
    AHMAX = -1.0E32
    AHMIN = 1.0E32
    ADMAX = -1.0E32
    ADMIN = 1.0E32
    
    AHDXY = ABS(AHD)       ! *** Default value
    
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING AHMAP.INP'
      OPEN(1,FILE='ahmap.inp')

      STR = READSTR(1)

      DO LL=2,LA_Global
        READ(1,*,END=200) LG, ITMP, JTMP, T1, T2
        
        R2D_Global(LG,1) = T1
        R2D_Global(LG,2) = T2
      ENDDO
200   CLOSE(1) 
    endif
    
    Call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        IF( R2D_Global(LG,1) < 0. )THEN
          AHOXY(L) = ABS(R2D_Global(LG,1))*DXP(L)*DYP(L)
        ELSE
          AHOXY(L) = R2D_Global(LG,1)
        ENDIF
        AHDXY(L) = R2D_Global(LG,2)
        
        AHMAX = MAX(AHOXY(L),AHMAX)
        AHMIN = MIN(AHOXY(L),AHMIN)
        ADMAX = MAX(AHDXY(L),ADMAX)
        ADMIN = MIN(AHDXY(L),ADMIN)
      ENDIF
    ENDDO
    Deallocate(R2D_Global)

    AHD = ABS(AHD)
    PRINT '(A,4E12.4)','VARIABLE AHO & AHD USED (MIN,MAX): ',AHMIN,AHMAX,ADMIN,ADMAX
  ELSE
    AHDXY = AHD       ! *** Constant value
  ENDIF
  
  ! *** VERTICAL VISCOSITY AND DIFFUSIVITY
  ! *** Constant and/or default values
  AVO = ABS(AVO)
  ABO = ABS(ABO)
  AVOXY = ABS(AVO)
  AVBXY = ABS(ABO)

  ! *** ENABLE SPATIALLY VARIABLE VERTICAL VISCOSITY AND DIFFUSIVITY
  IF( AVO < 0. )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, 2, 0.0)

    AVMAX = -1.0E32
    AVMIN = 1.0E32
    ABMAX = -1.0E32
    ABMIN = 1.0E32

    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING AVMAP.INP'
      OPEN(1,FILE='avmap.inp')
      
      STR = READSTR(1)

      DO LL=2,LA_Global
        Read(1,*,END=300)  LG, ID, JD, T1, T2

        R2D_Global(LG,1) = T1   ! *** AVOXY
        R2D_Global(LG,2) = T2   ! *** AVBXY
      ENDDO ! *** end loop over the file
300   CLOSE(1) 
    endif ! *** end on master
    
    Call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        AVOXY(L) = R2D_Global(LG,1)
        AVBXY(L) = R2D_Global(LG,2)
        
        AHMAX = MAX(AVOXY(L), AVMAX)
        AHMIN = MIN(AVOXY(L), AVMIN)
        ADMAX = MAX(AVBXY(L), ABMAX)
        ADMIN = MIN(AVBXY(L), ABMIN)
      ENDIF
      
      PRINT '(A,4E12.4)','VARIABLE AVO & ABO USED (MIN,MAX): ',AHMIN,AHMAX,ADMIN,ADMAX
    ENDDO
    Deallocate(R2D_Global)
  ENDIF
  
  ! *** OPEN FILE CHANSEC.INP FOR 1-D CHANNEL CROSS SECTION DATA
  ! *** REMOVED 2004-09-19  PMC

  ! *** OPEN FILE GWATER.INP TO SPECIFY GROUNDWATER INTERACTION
  ! *** BY INFILTRATION AND EVAPOTRANSPIRATION
  ISGWIE=0
  ! NDL ADDED 2018-03-21 TO FIXED ISSUE DIVIDE BY REZO RNPOR IN LINE 903 OF CALPUV2C.F90
  RNPOR=1.E-12
  IF( ISGWIT == 1 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING GWATER.INP'
      OPEN(1,FILE='gwater.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      READ(1,*) ISGWIE
      IF( ISGWIE >= 1 )THEN
        READ(1,*) DAGWZ,RNPOR,RIFTRM
      ELSE
        DAGWZ=0.0
        RNPOR=1.E-12
        RIFTRM=0.0
      ENDIF
      CLOSE(1)
    endif

    Call Broadcast_Scalar(DAGWZ, master_id)
    Call Broadcast_Scalar(RNPOR, master_id)
    Call Broadcast_Scalar(RIFTRM,master_id)

  ENDIF

339 FORMAT(2I5,6F14.5)

  ! *** OPEN FILE FBODY.INP TO READ IN SPATIALLY VARYING BODY FORCES
  IF( ISBODYF >= 1 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING FBODY.INP'
      OPEN(1,FILE='fbody.inp',STATUS='UNKNOWN')
      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      READ(1,*)CVTFACX,CVTFACY
    endif

    Call Broadcast_Scalar(CVTFACX, master_id)
    Call Broadcast_Scalar(CVTFACY, master_id)

    DO LL=2,LA_Global
      if( process_id == master_id )THEN
        READ(1,*)ITMP,JTMP,FBODY1,FBODY2
      endif

      Call Broadcast_Scalar(ITMP,   master_id)
      Call Broadcast_Scalar(JTMP,   master_id)
      Call Broadcast_Scalar(FBODY1, master_id)
      Call Broadcast_Scalar(FBODY2, master_id)

      ! *** Get local i,j
      ID = IG2IL(ITMP)
      JD = JG2JL(JTMP)

      if(id > 0 .and. id <= ic )THEN
        if(jd > 0 .and. jd <= jc )THEN

          L=LIJ(id, jd)
          IF( ISBODYF == 1 )THEN
            DO K=1,KC
              FBODYFX(L,K)=CVTFACX*FBODY1
              FBODYFY(L,K)=CVTFACY*FBODY2
            ENDDO
          ENDIF
          IF( ISBODYF == 2 )THEN
            DO K=1,KC-1
              FBODYFX(L,K)=0.0
              FBODYFY(L,K)=0.0
            ENDDO
            FBODYFX(L,KC)=CVTFACX*FBODY1
            FBODYFY(L,KC)=CVTFACY*FBODY2
          ENDIF
        endif
      endif
    ENDDO
    DO K=1,KC
      FBODYFX(1,K)=0.0
      FBODYFY(1,K)=0.0
      FBODYFX(LC,K)=0.0
      FBODYFY(LC,K)=0.0
    END DO

    if( process_id == master_id )THEN
      CLOSE(1)
    endif

  ENDIF

  ! *** OPEN FILE GWMAP.INP TO SPECIFY GROUNDWATER INTERACTION BY AMBIENT GROUNDWATER FLOW
  IF( ISGWIT > 1 )THEN
    ISGWIE = 0

    IF( ISGWIT == 3 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING GWSEEP.INP'
        OPEN(1,FILE='gwseep.inp',STATUS='UNKNOWN')

        STR=READSTR(1)
        READ(1,*)NSEEPCLASSES
        DO M=1,NSEEPCLASSES
          READ(1,*)IDUM,SEEPRATE(M)
        ENDDO
        CLOSE(1)
      endif
      
      Call Broadcast_Scalar(NSEEPCLASSES, master_id)
      Call Broadcast_Array(SEEPRATE,      master_id)
    ENDIF

    ! *** Read GWMAP
    ALLOCATE(I2D_Global(LCM_Global,3), R1D_Global(LCM_Global))
    I2D_Global = 0
    R1D_Global = 0.
    
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING GWMAP.INP'
      OPEN(1,FILE='gwmap.inp',STATUS='UNKNOWN')
      
      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)

      DO LL=2,LA_Global
        !READ(1,*) ITMP,JTMP,IZONE,RVALUE
        READ(1,*) (I2D_Global(LL,K),K=1,3), R1D_Global(LL)
      ENDDO
      CLOSE(1)
    endif
    
    Call Broadcast_Array(I2D_Global, master_id)
    Call Broadcast_Array(R1D_Global, master_id)

    ! *** Map to Local Domain
    DO LL=2,LA_GLOBAL
      ID = I2D_Global(LL,1)
      JD = I2D_Global(LL,2)
      LG = LIJ_Global(ID,JD)
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        IZONE  = I2D_Global(LL,3)
        RVALUE = R1D_Global(LL)
        IF( ISGWIT == 3 )THEN
          IF( IZONE > NSEEPCLASSES )THEN

            WRITE(mpi_log_unit,*)'BAD SEEPAGE CLASS AT I,J=',id,jd
            CALL STOPP('')
          ENDIF
          QGW(L) = SEEPRATE(IZONE)*RVALUE
        ELSE
          NGWSL(L) = IZONE
          GWFAC(L) = RVALUE
        ENDIF
      ENDIF
    ENDDO

    Deallocate(I2D_Global)
    Deallocate(R1D_Global)
  ENDIF

  ! *** READ IN SPATIALLY VARYING SEDIMENT ROUGHNESS HEIGHT FOR
  ! *** DETERMINING GRAIN STRESS
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    IF( ISBEDSTR == 3 )THEN
      Call AllocateDSI(R1D_Global, LCM_Global, 0.0)
      
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SEDROUGH.INP'
        OPEN(1,FILE='sedrough.inp')
        STR=READSTR(1)

        DO LG=2, LA_Global ! @todo make sure this is global
          READ(1,*) LDUM, IDUM, JDUM, R1D_Global(LG)                    ! *** ZBRSED
          IF( R1D_Global(LG) <= 0.0 )THEN
            STOP ' BAD SEDIMENT ROUGHNESS IN SEDROUGH.INP'
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      
      Call Broadcast_Array(R1D_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          ZBRSED(L) = R1D_Global(LG) 
        ENDIF
      ENDDO
      
      Deallocate(R1D_Global)
      
    ENDIF
    
  ENDIF

  ! *** OPEN FILE DOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** DISSOLVED ORGANIC CARBON IN WATER COLUMN
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2 ) IVAL=1
  ENDDO
  
  IF( IVAL == 1 )THEN
    IF( ISTDOCW == 1 )THEN
      Call AllocateDSI(R2D_Global, LCM_Global, KCM, 0.0)
      
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING DOCW.INP'
        OPEN(1,FILE='docw.inp',STATUS='UNKNOWN')
        
        STR=READSTR(1)
        READ(1,*) ISALTYP
      
        IF( ISALTYP == 0 )THEN
          DO LG=2,LA_Global
            READ(1,*,IOSTAT=ISO) (R2D_Global(LG,K),K=1,KC)                     ! *** STDOCW
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCW.INP')
          ENDDO
        ELSE
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (R2D_Global(L,K),K=1,KC)    ! *** STDOCW
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCW.INP')
          ENDDO
        ENDIF
        CLOSE(1) 
      endif

      Call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          STDOCW(L,:) = R2D_Global(LG,:)
        ENDIF
      ENDDO

      Deallocate(R2D_Global)
    ENDIF
    
  ENDIF

  ! *** OPEN FILE POCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON IN WATER COLUMN
  IVAL = 0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1 ) IVAL = 1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCW == 1 )THEN
      Call AllocateDSI(R2D_Global, LCM_Global, KCM, 0.0)
      
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING POCW.INP'
        OPEN(1,FILE='pocw.inp',STATUS='UNKNOWN')
        
        STR=READSTR(1)
        READ(1,*) ISALTYP

        IF( ISALTYP == 0 )THEN
          DO LG=2,LA_Global
            READ(1,*,IOSTAT=ISO) (R2D_Global(LG,K),K=1,KC)                    ! *** STPOCW
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCW.INP')
          ENDDO
        ELSE
          DO LG=2,LA_Global
            READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (R2D_Global(LG,K),K=1,KC)  ! *** STPOCW
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCW.INP')
          ENDDO
        ENDIF
        CLOSE(1)
      endif
      
      Call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          STPOCW(L,:) = R2D_Global(LG,:)
        ENDIF
      ENDDO
      Deallocate(R2D_Global)
    ENDIF

    ! *** PREPROCESS STPOCW TO PREVENT ANY DIVISION BY ZERO
    DO L=2,LC-1
      DO K=1,KC
        IF( STPOCW(L,K) <= 0.0 ) STPOCW(L,K) = 1E-12
      ENDDO
    ENDDO

  ENDIF

  ! *** OPEN FILE FPOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS
  ! *** IN WATER COLUMN
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 2 .OR. ISTOC(NT) == 3 ) IVAL=1
  ENDDO
  
  IF( IVAL == 1 )THEN
    IF( ISTPOCW == 3 )THEN
      NX = NSED + NSND
      Allocate(R3D_Global(LCM_Global,KCM,NX))
      R3D_Global = 0.0
      
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING FPOCW.INP'
        OPEN(1,FILE='fpocw.inp',STATUS='UNKNOWN')
        
        DO NS=1,NX
          STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
          READ(1,*) ISALTYP

          IF( ISALTYP == 0 )THEN
            DO LG=2,LA_Global
              READ(1,*,IOSTAT=ISO) (R3D_Global(LG,K,NS),K=1,KC)                                    ! *** STFPOCW
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCW.INP')
            ENDDO
          ELSE
            DO LG=2, LA_Global
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(R3D_Global(LG,K,NS),K=1,KC)                      ! *** STFPOCW
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCW.INP')
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif

      Call Broadcast_Array(R3D_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          STFPOCW(L,:,1:NX) = R3D_Global(LG,:,1:NX)
        ENDIF
      ENDDO
      Deallocate(R3D_Global)
    ENDIF
  ENDIF

  ! *** OPEN FILE DOCB.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT
  ! *** DISSOLVED ORGANIC CARBON IN SEDIMENT BED
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2 ) IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTDOCB == 1 )THEN
      Call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)
      STDOCB(:,:) = 0.0
      
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING DOCB.INP'
        OPEN(1,FILE='docb.inp',STATUS='UNKNOWN')
        
        STR=READSTR(1)
        READ(1,*) ISALTYP, IREAD, KBINPUT
      
        IF( IREAD == 0 )THEN
          IF( ISALTYP == 0 )THEN
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) R2D_Global(L,1)                                ! *** STDOCB
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              DO K=2,KB
                R2D_Global(L,K) = R2D_Global(L,1)
              ENDDO
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, R2D_Global(L,1)                    ! *** STDOCB
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              DO K=2,KB
                R2D_Global(L,K) = R2D_Global(L,1)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF( IREAD == 1 )THEN
          IF( ISALTYP == 0 )THEN
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) (R2D_Global(L,K),K=1,KB)                        ! *** STDOCB
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (R2D_Global(LDUM,K),K=1,KB)   ! *** STDOCB
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
            ENDDO
          ENDIF
        ENDIF

        IF( IREAD == 2 )THEN
          IF( ISALTYP == 0 )THEN
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) (R2D_Global(L,K),K=1,KBINPUT)                   ! *** STDOCB
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              DO K=KBINPUT,KB
                R2D_Global(L,K) = R2D_Global(L,KBINPUT)
              ENDDO
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (R2D_Global(L,K),K=1,KBINPUT) ! *** STDOCB
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DOCB.INP')
              DO K=KBINPUT,KB
                R2D_Global(L,K) = R2D_Global(L,KBINPUT)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        CLOSE(1)
      endif

      Call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          STDOCB(L,:) = R2D_Global(LG,:)
        ENDIF
      ENDDO
      Deallocate(R2D_Global)
    ENDIF
  ENDIF

  ! *** OPEN FILE POCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON IN BED
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1 ) IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCB == 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING POCB.INP'
        OPEN(1,FILE='pocb.inp',STATUS='UNKNOWN')
        STR=READSTR(1)
        READ(1,*) ISALTYP, IREAD, KBINPUT
      endif

      Call Broadcast_Scalar(ISALTYP, master_id)
      Call Broadcast_Scalar(IREAD, master_id)
      Call Broadcast_Scalar(KBINPUT, master_id)


      STPOCB(L,K)=0.0

      IF( IREAD == 0 )THEN
        ! ***
        IF( ISALTYP == 0 )THEN
          if( num_Processors > 1 )then
            write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
            pause
          endif
          
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO) STPOCB(L,1)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            DO K=2,KB
              STPOCB(L,K)=STPOCB(L,1)
            ENDDO
          ENDDO
        ELSE
          if( num_Processors > 1 )then
            write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
            pause
          endif
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STPOCB(L,1)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            DO K=2,KB
              STPOCB(L,K)=STPOCB(L,1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      
      IF( IREAD == 1 )THEN
        ! ***
        IF( ISALTYP == 0 )THEN
          if( num_Processors > 1 )then
            write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
            pause
          endif
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO) (STPOCB(L,K),K=1,KB)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
          ENDDO
        ELSE
          ! *** delme - this needs to be more general
          Call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)
          STDOCB(:,:) = 0.0
          
          if( process_id == master_id )then
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM,IDUM,JDUM,(STPOCB(1,K),K=1,KB)
              R2D_Global(LDUM,:) = STPOCB(1,:)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
            ENDDO
            STPOCB(1,:) = 0.0
          endif
          
          Call Broadcast_Array(R2D_Global, master_id)

          ! *** Map to Local Domain
          DO LG=2,LA_GLOBAL
            L = Map2Local(LG).LL
            IF( L > 1 )THEN
              STPOCB(L,:) = R2D_Global(LG,:)
            ENDIF
          ENDDO
          Deallocate(R2D_Global)
        ENDIF
      ENDIF
      
      IF( IREAD == 2 )THEN
        ! ***
        if( num_Processors > 1 )then
          write(*,*) 'WARNING: this POCB.INP format is not compatible with MPI version of EFDC+'
          pause
        endif

        IF( ISALTYP == 0 )THEN
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO) (STPOCB(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            DO K=KBINPUT,KB
              STPOCB(L,K)=STPOCB(L,KBINPUT)
            ENDDO
          ENDDO
        ELSE
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STPOCB(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE POCB.INP')
            DO K=KBINPUT,KB
              STPOCB(L,K)=STPOCB(L,KBINPUT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      if( process_id == master_id )THEN
        CLOSE(1)
      endif

    ENDIF
  ENDIF

  ! *** OPEN FILE FPOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS
  ! *** IN BED
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 2 .OR. ISTOC(NT) == 3) IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCB == 3 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING FPOCB.INP'
        OPEN(1,FILE='fpocb.inp',STATUS='UNKNOWN')
        DO NS=1,NSED+NSND
          STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
          READ(1,*) ISALTYP,IREAD,KBINPUT

          IF( IREAD == 0 )THEN
            IF( ISALTYP == 0 )THEN
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) STFPOCB(L,1,NS)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                DO K=2,KB
                  STFPOCB(L,K,NS)=STFPOCB(L,1,NS)
                ENDDO
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STFPOCB(L,1,NS)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                DO K=2,KB
                  STFPOCB(L,K,NS)=STFPOCB(L,1,NS)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF( IREAD == 1 )THEN
            IF( ISALTYP == 0 )THEN
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) (STFPOCB(L,K,NS),K=1,KB)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STFPOCB(L,K,NS),K=1,KB)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
              ENDDO
            ENDIF
          ENDIF
          IF( IREAD == 2 )THEN
            IF( ISALTYP == 0 )THEN
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) (STFPOCB(L,K,NS),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                DO K=KBINPUT,KB
                  STFPOCB(L,K,NS)=STFPOCB(L,KBINPUT,NS)
                ENDDO
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STFPOCB(L,K,NS),K=1,KBINPUT,NS)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FPOCB.INP')
                DO K=KBINPUT,KB
                  STFPOCB(L,K,NS)=STFPOCB(L,KBINPUT,NS)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        CLOSE(1)
      end if

      Call Broadcast_Array(STFPOCB, master_id)

    ENDIF
  ENDIF

  ! *** *******************************************************************C
  !###########################################################################
  ! HQI Change to include sptially varying, but time constant bulk foc
  ! FPOCB  - Bulk foc from data
  ! PFPOCB - Pseudo foc from data, to be used for all partitioning calculations
  ! RM, 02/29/04
  ! *** *******************************************************************C

  ! *** OPEN FILE FOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT
  ! *** PARTICULATE ORGANIC CARBON IN BED AND PSEUDO-POC IN BED
  IF( ISTPOCB == 4 )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)
    
    DO K=1,KB
      DO L=2,LA
        FPOCB(L,K)=0.0
      ENDDO
    ENDDO
      
    if( process_id == master_id )THEN               
      WRITE(*,'(A)')'READING FOCB.INP'                  
      OPEN(1,FILE='focb.inp',STATUS='UNKNOWN')
      
      STR=READSTR(1)
      READ(1,*) ISALTYP, IREAD, KBINPUT
      
      DO LG=2,La_Global
        READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (R2D_Global(LG,K),K=1,KBINPUT)    ! *** FPOCB
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE FOCB.INP')
      ENDDO
      CLOSE(1)
    end if
    
    Call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        DO K=1,KBINPUT
          FPOCB(L,K) = R2D_Global(LG,K)/1000000.
        ENDDO
        DO K=KBINPUT+1,KB
          FPOCB(L,K)=FPOCB(L,KBINPUT)
        ENDDO
      ENDIF
    ENDDO
      
    R2D_Global = 0.0
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING PSEUDO_FOCB.INP' 
      OPEN(1,FILE='pseudo_focb.inp',STATUS='UNKNOWN')
      
      STR=READSTR(1)
      READ(1,*) ISALTYP, IREAD, KBINPUT
      
      DO LG=2,LA_Global
        READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (R2D_Global(LG,K),K=1,KBINPUT)     ! *** PFPOCB
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSEUDO_FOCB.INP')
      ENDDO
      CLOSE(1)
    end if
    
    Call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        DO K=1,KBINPUT
          PFPOCB(L,K) = R2D_Global(LG,K)/1000000.
        ENDDO
        DO K=KBINPUT+1,KB
          PFPOCB(L,K)=PFPOCB(L,KBINPUT)
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)

    Deallocate(R2D_Global)
      
  ENDIF

  IF( ISVHEAT > 0 .AND. ISTOPT(2) /= 5 ) ISVHEAT = 0
  DO L=2,LA
    IF( ISTOPT(2) == 1 )THEN
      SVKEBACK(L) = SWRATNF
    ELSE
      SVKEBACK(L) = WQKEB(1)
    ENDIF
  ENDDO

  IF( ISTRAN(2) > 0 .AND. ISTOPT(2) == 5 .AND. ISVHEAT > 0 )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, 3, 0.0)
    
    ! *** DEFAULT VALUES
    DO L=2,LA
      SVREVC(L) = 0.001*ABS(REVC)
      SVRCHC(L) = 0.001*ABS(RCHC)
      IF( REVC < 0. ) LSVHTWINDE(L) = .TRUE.
      IF( RCHC < 0. ) LSVHTWINDC(L) = .TRUE.
    ENDDO

    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING SVHTFACT.INP'
      OPEN(1,FILE='svhtfact.inp',STATUS='UNKNOWN')  

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)

      ! *** LOOP OVER THE DOMAIN AND ONLY SET THE CELLS THAT ARE IN THE FILE
      DO LG=2,LA_Global
        READ(1,*,END=101) LL, ITMP, JTMP, R2D_Global(LG,1), R2D_Global(LG,2), R2D_Global(LG,3)
      ENDDO
101   CLOSE(1)
    endif

    Call Broadcast_Array(R2D_Global, master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        LSVHTWINDE(L) = R2D_Global(LG,1) < 0.
        LSVHTWINDC(L) = R2D_Global(LG,2) < 0.
        SVREVC(L)     = 0.001*ABS(R2D_Global(LG,1))
        SVRCHC(L)     = 0.001*ABS(R2D_Global(LG,2))
        SVKEBACK(L)   = R2D_Global(LG,3)
      ENDIF
    ENDDO
      
    Deallocate(R2D_Global)
      
    Call Broadcast_Array(LSVHTWINDE , master_id)
    Call Broadcast_Array(LSVHTWINDC , master_id)
    Call Broadcast_Array(SVREVC     , master_id)
    Call Broadcast_Array(SVRCHC     , master_id)
    Call Broadcast_Array(SVKEBACK   , master_id)

  ENDIF

  ! *** READ IN INITIAL SALINITY, TEMPERATURE, DYE, SED, SND, TOX
  ! *** FOR COLD STARTS FORM FILE XXXX.INP
  ! *** SALINITY
  DO K=1,KC
    DO L=2,LA
      SALINIT(L,K)=0.
    ENDDO
  ENDDO

  IF( ISTRAN(1) >= 1 .AND. (ISRESTI == 0 .OR. (ISRESTI >= 1 .AND. ISCI(1) == 0) .OR. (ISTOPT(1) > 1))  )THEN
    IF( ISTOPT(1) >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SALT.INP'
        OPEN(1,FILE='salt.inp',STATUS='UNKNOWN')

        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR=READSTR(1)
        READ(1,*)ISALTYP

        IF( ISALTYP == 0 )THEN
          DO LG=2, LA_GLOBAL
            READ(1,*,IOSTAT=ISO) (SAL_Global(LG,K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SALT.INP')
          ENDDO
        ELSE
          DO L=2, LA_GLOBAL
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(CONINIT(K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SALT.INP')

            LG = LIJ_Global(IDUM,JDUM)
            SAL_Global(LG,:) = CONINIT(:)
          ENDDO
        ENDIF
        CLOSE(1)
      End if

      Call Broadcast_Array(SAL_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          SALINIT(L,:) = SAL_Global(LG,:)
        ENDIF
      ENDDO
    ENDIF

  ENDIF

  ! *** TEMPERATURE
  TEM_Global = TEMO
  DO K=1,KC
    DO L=2,LA
      TEMINIT(L,K) = TEMO
    ENDDO
  ENDDO

  IF( ISTRAN(2) >= 1 .AND. (ISRESTI == 0 .OR. (ISRESTI >= 1 .AND. ISCI(2) == 0) .OR. (ISTOPT(2) > 9) ) )THEN
    IF( ISTOPT(2) >= 1 .OR. INITTEMP > 0 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING TEMP.INP'
        OPEN(1,FILE='temp.inp',STATUS='UNKNOWN')
        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        READ(1,*) ISALTYP

        IF( ISALTYP == 0 )THEN
          DO LG=2,LA_Global
            READ(1,*,IOSTAT=ISO) (TEM_Global(LG,K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TEMP.INP')
          ENDDO
        ELSE
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONINIT(K),K=1,KC)

            LG = LIJ_Global(IDUM,JDUM)
            TEM_Global(LG,:) = CONINIT(:)
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
      Call Broadcast_Array(TEM_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          TEMINIT(L,:) = TEM_Global(LG,:)
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! *** DYE
  IF( ISTRAN(3) >= 1 .AND. (ISRESTI == 0 .OR. (ISRESTI >= 1 .AND. ISCI(3) == 0) .OR. (ISTOPT(3) > 1)) )THEN
    IF( ISTOPT(3) >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING DYE.INP'
        OPEN(1,FILE='dye.inp',STATUS='UNKNOWN')

        DO MD = 1,NDYE
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*)ISALTYP

          IF( ISALTYP == 0 )THEN
            DO L=LG,LA_Global
              READ(1,*,IOSTAT=ISO) (DYE_Global(LG,K,MD),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DYE.INP')
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(CONINIT(K),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DYE.INP')

              LG = LIJ_Global(IDUM,JDUM)
              DYE_Global(LG,:,MD) = CONINIT(:)
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      end if
      Call Broadcast_Array(DYE_Global, master_id)

      ! *** Map to Local Domain
      DO MD = 1,NDYE
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            DYEINIT(L,:,MD) = DYE_Global(LG,:,MD)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ENDIF

  ! *** SFL
  IF( ISRESTI == 0 .AND. ISTRAN(4) >= 1 )THEN
    IF( ISTOPT(4) >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SFL.INP'
        OPEN(1,FILE='sfl.inp',STATUS='UNKNOWN')

        ! ***   SKIP OVER TITLE AND AND HEADER LINES
        STR = READSTR(1)
        READ(1,*) ISALTYP

        IF( ISALTYP == 0 )THEN
          DO LG=2,LA_Global
              READ(1,*,IOSTAT=ISO) (SFL_Global(LG,K),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFL.INP')          
          ENDDO
        ELSE
          DO L=2,LA_Global
            READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONINIT(K),K=1,KC)
          
            LG = LIJ_Global(IDUM,JDUM)
            SFL_Global(LG,:) = CONINIT(:)
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF

      Call Broadcast_Array(SFL_Global, master_id)
      
      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          SFLINIT(L,:) = SFL_Global(LG,:)
        ENDIF
      ENDDO

    ENDIF
  ENDIF

  ! ****************************************************************************************
  ! *** SEDIMENTS AND SEDIMENT BED FILES
  
  ! *** INITIALIZE HARD-BOTTOM, IF USED
  LASED = 0
  LDMSED = 0
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    ALLOCATE(BEDMAP(LCM))
    ALLOCATE(LSED(LCM))
    ALLOCATE(LBED(LCM))

    BEDMAP = 1                ! *** BED PROCESSES ACTIVE
    LBED = .FALSE.
    DO L=2,LA
      LASED = LASED + 1
      LSED(LASED) = L
    ENDDO

    IF( ISBEDMAP > 0 )THEN
      ! *** READ USER SPECIFIED SEDIMENT ACTIVE CELL LIST
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING BEDMAP.INP'
        OPEN(1,FILE='bedmap.inp',STATUS='UNKNOWN')

        STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
        DO L=2,LA_GLOBAL
          READ(1,*,IOSTAT=ISO,END=400) IDUM,JDUM,NDUM1
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDMAP.INP')
          
          LG = LIJ_Global(IDUM,JDUM)
          BEDMAP_Global(LG) = NDUM1
        ENDDO
400     CLOSE(1)
      endif
      Call Broadcast_Array(BEDMAP_Global, master_id)
      
      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          BEDMAP(L) = BEDMAP_Global(LG)
        ENDIF
      ENDDO

      LASED = 0
      LSED = 0
      LBED = .TRUE.
      DO L=2,LA
        IF( BEDMAP(L) > 0 )THEN
          LASED = LASED + 1
          LSED(LASED) = L
          LBED(L) = .FALSE.
        ENDIF
      ENDDO
    ENDIF
    
    LDMSED = INT(LASED/NTHREADS)+1
    
    ! *** OPEN FILE SEDBLBC.INP TO READ IN SEDIMENT BEDLOAD OUTFLOW
    ! *** OR RECIRCULATION BOUNDARY CONDITIONS
    NSBDLDBC = 0
    IF( (ICALC_BL == 1 .AND. .NOT. LSEDZLJ) .OR. (ICALC_BL >= 2 .AND. LSEDZLJ) )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SEDBLBC.INP'
        OPEN(1,FILE='sedblbc.inp',STATUS='UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR=READSTR(1)
        READ(1,*) NDUM
        
        ALLOCATE(IBLTMP(5,NDUM))
        IBLTMP = 0
        
        DO NS=1,NDUM
          READ(1,*) (IBLTMP(I,NS),I=1,5)
        ENDDO
        CLOSE(1)
      endif
      call Broadcast_Scalar(NDUM, master_id)
      
      if( process_id /= master_id )THEN
        ALLOCATE(IBLTMP(5,NDUM))
        IBLTMP = 0
      endif
      call Broadcast_Array(IBLTMP, master_id)
      
      DO NS=1,NDUM
        ! *** Check if in current node
        LG = LIJ_Global(IBLTMP(1,NS),IBLTMP(2,NS))
        LU = Map2Local(LG).LL
        IF( LU > 0 )THEN
          ! *** ENSURE US CELLS ARE PART OF THE ACTIVE MAP, OTHERWISE DISREGARD
          DO LL= 1,LASED
            IF( LSED(LL) == LU ) EXIT
          ENDDO  
          IF( LL > LASED ) CYCLE

          NSBDLDBC = NSBDLDBC + 1
          ISDBLDIR(NSBDLDBC) = ABS(IBLTMP(5,NS))
          LSBLBCU(NSBDLDBC) = LU
            
          IF( IBLTMP(3,NS) > 0 .AND. IBLTMP(4,NS) > 0 )THEN
            LG = LIJ_Global(IBLTMP(3,NS),IBLTMP(4,NS))
            LD = Map2Local(LG).LL
            LSBLBCD(NSBDLDBC) = LD
          ELSE
            LSBLBCD(NSBDLDBC) = 0
          ENDIF
        ENDIF
      ENDDO
      
      DO NS=1,NSBDLDBC
        ISDBLDIR(NS) = ABS(ISDBLDIR(NS))
        IF( ISDBLDIR(NS) == 1 )THEN
          ! *** EAST-WEST
          IF( SUB(LSBLBCU(NS)) < 0.5 .AND. SUB(LEC(LSBLBCU(NS))) > 0.5 ) ISDBLDIR(NS) = -1
        ENDIF
        IF( ISDBLDIR(NS) == 2 )THEN
          ! *** NORTH-SOUTH
          IF( SVB(LSBLBCU(NS)) < 0.5 .AND. SVB(LEC(LSBLBCU(NS))) > 0.5 ) ISDBLDIR(NS) = -2
        ENDIF
      ENDDO

    ENDIF
  ENDIF

  ! ****************************************************************************************
  ! *** TOXICS
  IF( ISTRAN(5) >= 1 )THEN
    DO NT=1,NTOX
      DO K=1,KC
        DO L=2,LA
          TOXINIT(L,K,NT)=TOXINTW(NT)
        ENDDO
      ENDDO
    ENDDO
    DO NT=1,NTOX
      DO K=1,KB
        DO L=2,LA
          TOXBINIT(L,K,NT)=TOXINTB(NT)
        ENDDO
      ENDDO
    ENDDO

    ! *** SPATIALLY VARYING WATER COLUMN INITIAL CONDITIONS
    ISCOLD = 1 
    IF( (ISRESTI >= 1 .AND. ISCI(5) > 0) .OR. ISLTMT > 0 ) ISCOLD = 0
    IF( ISCOLD == 1 )THEN
      IF( ITXINT(1) == 2 .OR. ITXINT(1) == 3 )THEN
        if( process_id == master_id )THEN
          WRITE(*,'(A)')'READING TOXW.INP'
          OPEN(1,FILE='toxw.inp',STATUS='UNKNOWN')
  
          DO NT=1,NTOX
            STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
            READ(1,*) ISALTYP, ITOXWU(NT)
              
            IF( ISALTYP == 0 )THEN
              DO LG=2,LA_Global
                READ(1,*,IOSTAT=ISO) (TOX_Global(LG,K,NT),K=1,KC)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXW.INP')
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONINIT(K),K=1,KC)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXW.INP')

                LG = LIJ_Global(IDUM,JDUM)
                TOX_Global(LG,:,NT) = CONINIT(:)
              ENDDO
            ENDIF
          ENDDO
          CLOSE(1)
        endif
        
        Call Broadcast_Array(ITOXWU,     master_id)
        Call Broadcast_Array(TOX_Global, master_id)

        ! *** Map to Local Domain
        DO NT = 1,NTOX
          DO LG=2,LA_GLOBAL
            L = Map2Local(LG).LL
            IF( L > 1 )THEN
              TOXINIT(L,:,NT) = TOX_Global(LG,:,NT)
            ENDIF
          ENDDO
        ENDDO

        ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
        if( process_id == master_id )THEN
          WRITE(*,'(A)')'READING TOXB.INP'
          OPEN(1,FILE='toxb.inp',STATUS='UNKNOWN')
  
            DO NT=1,NTOX
              STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
              ! *** ITOXBU is not used
              READ(1,*) ISALTYP,ITOXBU(NT),KBINPUT
                
              DO K=1,KB
                DO LG=2,LA_Global
                  TOXB_Global(LG,K,NT)=0.0
                ENDDO
              ENDDO
              
              IF( ISALTYP == 0 )THEN
                DO LG=2,LA_Global
                  READ(1,*,IOSTAT=ISO) (TOXB_Global(LG,K,NT),K=1,KBINPUT)
                  IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXB.INP')
                ENDDO
              ELSE
                DO L=2,LA_Global
                  READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONBINIT(K),K=1,KBINPUT)
                  IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TOXB.INP')
            
                  LG = LIJ_Global(IDUM,JDUM)
                  TOXB_Global(LG,1:KBINPUT,NT) = CONBINIT(1:KBINPUT)
                ENDDO
              ENDIF
              
              ! *** FILL ANY MISSING LAYERS
              DO K = KBINPUT+1,KB
                TOXB_Global(LG,K,NT) = TOXB_Global(LG,KBINPUT,NT)
              ENDDO
            ENDDO
          CLOSE(1)
        end if
        Call Broadcast_Array(ITOXBU,      master_id)
        Call Broadcast_Array(TOXB_Global, master_id)

        ! *** Map to Local Domain
        DO NT = 1,NTOX
          DO LG=2,LA_GLOBAL
            L = Map2Local(LG).LL
            IF( L > 1 )THEN
              TOXBINIT(L,:,NT) = TOXB_Global(LG,:,NT)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDIF

  ! ****************************************************************************************
  ! *** COHESIVE/SEDZLJ SEDIMENTS
  IF( ISTRAN(6)  >= 1 )THEN
    DO NS=1,NSED
      DO K=1,KC
        DO L=2,LA
          SEDINIT(L,K,NS)=SEDO(NS)
        ENDDO
      ENDDO
    ENDDO
    DO NS=1,NSED
      DO K=1,KB
        DO L=2,LA
          SEDBINIT(L,K,NS)=SEDBO(NS)
        ENDDO
      ENDDO
    ENDDO

    ! *** SPATIALLY VARYING WATER COLUMN INITIAL CONDITIONS
    ISEDINIT = 0                        ! *** SPATIALLY VARYING WC OFF
    IF( ISEDINT == 1 ) ISEDINIT = 1     ! *** SPATIALLY VARYING WC ON
    IF( ISEDINT == 3 ) ISEDINIT = 1

    ISCOLD = 1 
    IF( (ISRESTI >= 1 .AND. ISCI(6) > 0) .OR. ISLTMT > 0 ) ISCOLD = 0
    IF( ISCOLD == 1 .AND. ISEDINIT >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SEDW.INP'
        OPEN(1,FILE='sedw.inp',STATUS='UNKNOWN')
          
        DO NS=1,NSED
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*)ISALTYP

          IF( ISALTYP == 0 )THEN
            DO LG=2,LA_Global
              READ(1,*,IOSTAT=ISO) (SED_Global(LG,K,NS),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDW.INP')
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONINIT(K),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDW.INP')

              LG = LIJ_Global(IDUM,JDUM)
              SED_Global(LG,:,NS) = CONINIT(:)
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      Call Broadcast_Array(SED_Global, master_id)

      ! *** Map to Local Domain
      DO NS = 1,NSED
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            SEDINIT(L,:,NS) = SED_Global(LG,:,NS)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
    ISEDINIT = 0                        ! *** SPATIALLY VARYING BED OFF
    IF( ISEDINT == 2 ) ISEDINIT = 1     ! *** SPATIALLY VARYING BED ON
    IF( ISEDINT == 3 ) ISEDINIT = 1

    IF( NSEDFLUME == 0 )THEN
      ISCOLD = 1   
      IF( (ISRESTI >= 1 .AND. ISCI(6) > 0) .OR. ISLTMT > 0 ) ISCOLD = 0
    ELSEIF( NSEDFLUME == 3 .AND. IHTSTRT == 0 )THEN
      ISCOLD = 1
      ISEDINIT = 1
    ELSE
      ISCOLD = 0
    ENDIF
      
    IF( ISCOLD == 1 .AND. ISEDINIT >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SEDB.INP'
        OPEN(1,FILE='sedb.inp',STATUS='UNKNOWN')
          
        DO NS=1,NSED
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*) ISALTYP, ISEDBU(NS), KBINPUT
            
          DO K=1,KB
            DO LG=2,LA_Global
              SEDB_Global(LG,K,NS)=0.0
            ENDDO
          ENDDO
            
          IF( ISALTYP == 0 )THEN
            DO LG=2,LA_Global
              READ(1,*,IOSTAT=ISO) (SEDB_Global(LG,K,NS),K=1,KBINPUT)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDB.INP')
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONBINIT(K),K=1,KBINPUT)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SEDB.INP')
          
              LG = LIJ_Global(IDUM,JDUM)
              SEDB_Global(LG,1:KBINPUT,NS) = CONBINIT(1:KBINPUT)
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      Call Broadcast_Array(ISEDBU,      master_id)
      Call Broadcast_Array(SEDB_Global, master_id)

      ! *** Map to Local Domain
      DO NS = 1,NSED
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            SEDBINIT(L,:,NS) = SEDB_Global(LG,:,NS)
          ENDIF
        ENDDO
      ENDDO
        
    ENDIF
  ENDIF   ! *** END OF NON-COHESIVE CONFIGURATION

  ! ****************************************************************************************
  ! *** NON-COHESIVE SEDIMENTS
  IF( ISTRAN(7) >= 1 )THEN
    DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KC
        DO L=2,LA
          SNDINIT(L,K,NX) = SEDO(NS)
        ENDDO
      ENDDO
    ENDDO
    DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KB
        DO L=2,LA
          SNDBINIT(L,K,NX) = SEDBO(NS)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** SPATIALLY VARYING WATER COLUMN INITIAL CONDITIONS
    ISEDINIT = 0                        ! *** SPATIALLY VARYING WC OFF
    IF( ISEDINT == 1 ) ISEDINIT = 1     ! *** SPATIALLY VARYING WC ON
    IF( ISEDINT == 3 ) ISEDINIT = 1
    
    ISCOLD = 1 
    IF( (ISRESTI >= 1 .AND. ISCI(7) > 0) .OR. ISLTMT > 0 ) ISCOLD = 0
    IF( ISCOLD == 1 .AND. ISEDINIT >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SNDW.INP'
        OPEN(1,FILE='sndw.inp',STATUS='UNKNOWN')
          
        DO NX=1,NSND
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*)ISALTYP
          IF( ISALTYP == 0 )THEN
            DO LG=2,LA_Global
              READ(1,*,IOSTAT=ISO) (SND_Global(LG,K,NX),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDW.INP')
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONINIT(K),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDW.INP')

              LG = LIJ_Global(IDUM,JDUM)
              SND_Global(LG,:,NX) = CONINIT(:)
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      Call Broadcast_Array(SND_Global, master_id)

      ! *** Map to Local Domain
      DO NX = 1,NSND
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            SNDINIT(L,:,NX) = SND_Global(LG,:,NX)
          ENDIF
        ENDDO
      ENDDO

    ENDIF

    ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
    ISEDINIT = 0                        ! *** SPATIALLY VARYING BED OFF
    IF( ISEDINT == 2 ) ISEDINIT = 1     ! *** SPATIALLY VARYING BED ON
    IF( ISEDINT == 3 ) ISEDINIT = 1

    IF( NSEDFLUME == 0 )THEN
      ISCOLD = 1   
      IF( (ISRESTI >= 1 .AND. ISCI(7) > 0) .OR. ISLTMT > 0 ) ISCOLD = 0
    ELSE
      ISCOLD = 0
    ENDIF
      
    IF( ISCOLD == 1 .AND. ISEDINIT >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SNDB.INP'
        OPEN(1,FILE='sndb.inp',STATUS='UNKNOWN')
        DO NX=1,NSND
          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*) ISALTYP,ISNDBU(NX),KBINPUT
            
          DO K=1,KB
            DO LG=2,LC_Global
              SNDB_Global(LG,K,NX) = 0.0
            ENDDO
          ENDDO

          IF( ISALTYP == 0 )THEN
            DO LG=2,LA_Global
              READ(1,*,IOSTAT=ISO) (SNDB_Global(LG,K,NX),K=1,KBINPUT)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDB.INP')
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*,IOSTAT=ISO) LDUM, IDUM, JDUM, (CONBINIT(K),K=1,KBINPUT)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNDB.INP')
          
              LG = LIJ_Global(IDUM,JDUM)
              SNDB_Global(LG,1:KBINPUT,NX) = CONBINIT(1:KBINPUT)
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif

      Call Broadcast_Array(ISNDBU,      master_id)
      Call Broadcast_Array(SNDB_Global, master_id)

      ! *** Map to Local Domain
      DO NX = 1,NSND
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            SNDBINIT(L,:,NX) = SNDB_Global(LG,:,NX)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF   ! *** END OF NON-COHESIVE CONFIGURATION

  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    IFLAG = 0
    DO NS=1,NSED
      IF( ISEDBINT == 1 .AND. ISEDBU(NS) /= 1 ) IFLAG = 1
    ENDDO
    DO NX=1,NSND
      IF( ISEDBINT == 1 .AND. ISNDBU(NX) /= 1 ) IFLAG = 1
    ENDDO

    IF( IFLAG == 0 )THEN
      ! *** CHECK FRACTIONS
      DO L=2,LA
        DO K=1,KB
          T1 = 0.0
          DO NS=1,NSED
            T1 = T1 + SEDBINIT(L,K,NS)
          ENDDO
          DO NX=1,NSND
            T1 = T1 + SNDBINIT(L,K,NX)
          ENDDO

          IF( ABS(1.0 - T1) > 1E-6 .AND. T1 > 0.0 )THEN
            ! *** FRACTIONS DO NOT ADD UP.  SET FRACTIONS
            DO NS=1,NSED
              SEDBINIT(L,K,NS) = SEDBINIT(L,K,NS)/T1
            ENDDO
            DO NX=1,NSND
              SNDBINIT(L,K,NX) = SNDBINIT(L,K,NX)/T1
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! ****************************************************************************************
  !  ** SEDIMENT BED LAYER CONFIGURATION
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    ! *** SPATIALLY VARYING SEDIMENT BED INITIAL CONDITIONS
    ISEDINIT = 0                        ! *** SPATIALLY VARYING BED OFF
    IF( ISEDINT == 2 ) ISEDINIT = 1     ! *** SPATIALLY VARYING BED ON
    IF( ISEDINT == 3 ) ISEDINIT = 1

    IF( NSEDFLUME == 0 )THEN
      ISCOLD = 1
      IF( ISLTMT > 0 )THEN
        ISCOLD = 0
      ELSEIF( ISRESTI >= 1 )THEN 
       IF( ISTRAN(6) > 0 .AND. ISCI(6) > 0 ) ISCOLD = 0
       IF( ISTRAN(7) > 0 .AND. ISCI(7) > 0 ) ISCOLD = 0
      ENDIF
    ELSEIF( NSEDFLUME == 3 .AND. IHTSTRT == 0 )THEN
      ISCOLD = 1
      ISEDINIT = 1
    ELSE
      ISCOLD = 0
    ENDIF
      
    IF( ISCOLD == 1 .AND. ISEDINIT >= 1 )THEN
      IF( IBMECH /= 9 )THEN
        !  ** BED LAYER THICKNESS
        Call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)

        if( process_id == master_id )THEN
          WRITE(*,'(A)')'READING BEDLAY.INP'
          OPEN(1,FILE='bedlay.inp',STATUS='UNKNOWN')

          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*) IBEDLAYU, ISALTYP, KBINPUT
          IF( IBEDLAYU > 0 )THEN
            DO K=1,KB
              DO L=2,LA_Global
                R2D_Global(L,K) = 0.0
              ENDDO
            ENDDO
            IF( ISALTYP == 0 )THEN
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) (R2D_Global(L,K),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDLAY.INP')
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(R2D_Global(L,K),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDLAY.INP')
              ENDDO
            ENDIF
          ENDIF
          CLOSE(1)
        endif   ! *** End master process

        Call Broadcast_Scalar(IBEDLAYU,  master_id)
        Call Broadcast_Array(R2D_Global, master_id)

        ! *** Map to Local Domain
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            BEDLINIT(L,:) = R2D_Global(LG,:)
          ENDIF
        ENDDO
        DEALLOCATE(R2D_Global)

        ! *** BED LAYER BULK DENSITY
        Call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)

        if( process_id == master_id )THEN
          WRITE(*,'(A)')'READING BEDBDN.INP'
          OPEN(1,FILE='bedbdn.inp',STATUS='UNKNOWN')

          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*) IBEDBDNU, ISALTYP, KBINPUT
        
          IF( IBEDBDNU > 0 )THEN
            DO K=1,KB
              DO L=2,LA_Global
                R2D_Global(L,K)=0.0
              ENDDO
            ENDDO
            IF( ISALTYP == 0 )THEN
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) (R2D_Global(L,K),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDBDN.INP')
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(R2D_Global(L,K),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDBDN.INP')
              ENDDO
            ENDIF
          ENDIF
          CLOSE(1)
        endif   ! *** End master process

        Call Broadcast_Scalar(IBEDBDNU,  master_id)
        Call Broadcast_Array(R2D_Global, master_id)

        ! *** Map to Local Domain
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            BEDBINIT(L,:) = R2D_Global(LG,:)
          ENDIF
        ENDDO
        DEALLOCATE(R2D_Global)

        ! *** BED LAYER DRY DENSITY, POROSITY OR VOID RATIO
        Call AllocateDSI(R2D_Global, LCM_Global, KBM, 0.0)

        if( process_id == master_id )THEN
          WRITE(*,'(A)')'READING BEDDDN.INP'
          OPEN(1,FILE='bedddn.inp',STATUS='UNKNOWN')

          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*) IBEDDDNU, ISALTYP, KBINPUT
          IF( IBEDDDNU > 0 )THEN
            DO K=1,KB
              DO L=2,LA_Global
                R2D_Global(L,K)=0.0
              ENDDO
            ENDDO
            IF( ISALTYP == 0 )THEN
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO) (R2D_Global(L,K),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDDDN.INP')
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(R2D_Global(L,K),K=1,KBINPUT)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE BEDDDN.INP')
              ENDDO
            ENDIF
          ENDIF
          CLOSE(1)
        end if

        Call Broadcast_Scalar(IBEDDDNU,  master_id)
        Call Broadcast_Array(R2D_Global, master_id)
      
        ! *** Map to Local Domain
        DO LG=2,LA_GLOBAL
          L = Map2Local(LG).LL
          IF( L > 1 )THEN
            BEDDINIT(L,:) = R2D_Global(LG,:)
          ENDIF
        ENDDO
        DEALLOCATE(R2D_Global)
      ENDIF   ! *** End IBMECH /= 9

      ! *** CONSOLIDATION MAP
      IF( IBMECH == 9 .AND. .NOT. LSEDZLJ )THEN
        if( process_id == master_id )THEN
          WRITE(*,'(A)')'READING CONSOLMAP.INP'
          OPEN(1,FILE='consolmap.inp',STATUS='UNKNOWN')

          ! ***   SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*)ISALTYP
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) LCONSOL(L)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CONSOLMAP.INP')
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,LCONSOL(L)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CONSOLMAP.INP')
            ENDDO
          ENDIF
          CLOSE(1)
        end if

        Call Broadcast_Array(LCONSOL, master_id)
      ENDIF
    ENDIF    ! *** END OF COLD START
  ENDIF      ! *** END OF SEDIMENT BED CONFIGURATION

  ! *** READ IN OPEN BOUNDARY SURFACE ELEVATION TIME SERIES FROM THE
  ! *** FILE PSER.INP
  IF( NPSER >= 1 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING PSER.INP'
      OPEN(1,FILE='pser.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      DO NS=1,NPSER
        READ(1,*,IOSTAT=ISO) ITYPE,MPSER(NS),TCPSER(NS),TAPSER(NS),RMULADJ,ADDADJ,PSERZDF(NS),INTPSER(NS)
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
        PSERZDF(NS)=G*PSERZDF(NS)
        IF( ITYPE == 1 )THEN
          READ(1,*,IOSTAT=ISO) RMULADJS,ADDADJS,PSERZDS(NS)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
          PSERZDS(NS)=G*PSERZDS(NS)
        ELSE
          RMULADJS=0
          ADDADJS=0.
          PSERZDS(NS)=0.
        ENDIF
        IF( ITYPE == 0 )THEN
          DO M=1,MPSER(NS)
            READ(1,*,IOSTAT=ISO)TPSER(M,NS),PSERTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
            TPSER(M,NS) = TPSER(M,NS)+TAPSER(NS)
            PSER(M,NS)  = G*(PSERTMP+ADDADJ)*RMULADJ                  ! *** m2/s2
            PSERS(M,NS) = 0.0
          ENDDO
        ELSEIF( ITYPE == 1 )THEN
          DO M = 1,MPSER(NS)
            READ(1,*,IOSTAT = ISO)TPSER(M,NS),PSERTMP,PSERTMPS
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE PSER.INP')
            TPSER(M,NS) = TPSER(M,NS)+TAPSER(NS)
            PSER(M,NS)  = G*(PSERTMP+ADDADJ)*RMULADJ                  ! *** m2/s2
            PSERS(M,NS) = G*(PSERTMPS+ADDADJS)*RMULADJS               ! *** m2/s2
          ENDDO
        ENDIF

        ! *** Compute the average pressure value for this series.  Used for radiation BC's
        DO M = 2,MPSER(NS)
          TDELTA = TPSER(M,NS) - TPSER(M-1,NS)
          PSERAVG(NS,1) = PSERAVG(NS,1) + PSER(M,NS)*TDELTA
          PSERAVG(NS,2) = PSERAVG(NS,2) + PSERS(M,NS)*TDELTA
        ENDDO
        PSERAVG(NS,1) = PSERAVG(NS,1)/(TPSER(MPSER(NS),NS) - TPSER(1,NS))
        PSERAVG(NS,2) = PSERAVG(NS,2)/(TPSER(MPSER(NS),NS) - TPSER(1,NS))
      ENDDO
      CLOSE(1)
    endif

    Call Broadcast_Array(MPSER,   master_id)
    Call Broadcast_Array(TCPSER,  master_id)
    Call Broadcast_Array(TAPSER,  master_id)
    Call Broadcast_Scalar(RMULADJ,master_id)
    Call Broadcast_Scalar(ADDADJ, master_id)

    Call Broadcast_Array(PSERZDS, master_id)
    Call Broadcast_Array(PSERZDF, master_id)
    Call Broadcast_Array(INTPSER, master_id)

    Call Broadcast_Array(TPSER,   master_id)
    Call Broadcast_Array(PSER,    master_id)
    Call Broadcast_Array(PSERAVG, master_id)
    Call Broadcast_Array(PSERS,   master_id)

  ENDIF
6776 FORMAT(A20)

  ! *** READ IN VOLUMETRIC SOURCE OR RIVER INFLOW TIME SERIES FROM THE
  ! *** FILE QSER.INP
  IF( NQSER >= 1 )THEN
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING QSER.INP'
      OPEN(1,FILE='qser.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      DO NS=1,NQSER
        ISMOOTH=0
        READ(1,*,IOSTAT=ISO) ISTYP, MQSER(NS), TCQSER(NS), TAQSER(NS), RMULADJ, ADDADJ, ICHGQS
        IF( MQSER(NS)<0 )THEN
          ISMOOTH=1
          MQSER(NS)=-MQSER(NS)
        ENDIF
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')

        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')
          DO M=1,MQSER(NS)
            READ(1,*,IOSTAT=ISO)QSER(NS).TIM(M),QSERTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')
            QSER(NS).TIM(M) = QSER(NS).TIM(M)+TAQSER(NS)
            QSERTMP = RMULADJ*(QSERTMP+ADDADJ)
            IF( ICHGQS == 1)  QSERTMP = MAX(QSERTMP,0.0)
            IF( ICHGQS == -1) QSERTMP = MIN(QSERTMP,0.0)
            DO K=1,KC
              QSER(NS).VAL(M,K) = QSERTMP*WKQ(K)
            ENDDO
          ENDDO
        ELSE
          DO M=1,MQSER(NS)
            READ(1,*,IOSTAT=ISO)QSER(NS).TIM(M),(QSER(NS).VAL(M,K), K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QSER.INP')
            QSER(NS).TIM(M) = QSER(NS).TIM(M)+TAQSER(NS)
            DO K=1,KC
              QSER(NS).VAL(M,K) = RMULADJ*(QSER(NS).VAL(M,K)+ADDADJ)
              IF( ICHGQS == 1)  QSER(NS).VAL(M,K) = MAX(QSER(NS).VAL(M,K),0.0)
              IF( ICHGQS == -1) QSER(NS).VAL(M,K) = MIN(QSER(NS).VAL(M,K),0.0)
            ENDDO
          ENDDO
        ENDIF

        IF( ISMOOTH == 1 )THEN
          DO K=1,KC
            DO M=1,MQSER(NS)
              QSERSM(M,K) = QSER(NS).VAL(M,K)
            ENDDO
          ENDDO
          DO K=1,KC
            DO M=2,MQSER(NS)-1
              QSER(NS).VAL(M,K) = 0.25*(QSERSM(M-1,K)+QSERSM(M+1,K))+0.5*QSERSM(M,K)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    endif

    call Broadcast_Array(MQSER    , master_id)
    call Broadcast_Array(TCQSER    , master_id)
    call Broadcast_Array(TAQSER    , master_id)
    
    DO NS=1,NQSER
      Call Broadcast_Scalar(QSER(NS).NREC , master_id)
      call Broadcast_Array(QSER(NS).TIM   , master_id)
      call Broadcast_Array(QSER(NS).VAL   , master_id)
    ENDDO      
    
  ENDIF
2222 FORMAT(2I5,F12.7,F12.4)

  ! *** READ IN FLOW WITHDRAWL-RETURN FLOW AND CONCENTRATION RISE
  ! *** TIME SERIES FROM THE FILE QWRS.INP
  IF( NQWRSR >= 1 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING QWRS.INP'
      OPEN(1,FILE='qwrs.inp',STATUS='UNKNOWN')

      NCTMP = 3 + NDYM + NSED + NSND + NTOX
      IF( ISTRAN(8) > 0 ) NCTMP = NCTMP + NWQV    ! *** Includes inactive WQ parameters also

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      DO NS=1,NQWRSR
        READ(1,*,IOSTAT=ISO) ISTYP, MQWRSR(NS), TCQWRSR(NS), TAQWRSR(NS), RMULADJ, ADDADJ
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QWRS.INP')
        IF( ISTYP == 0 )THEN
          ! *** FLOW ONLY.  NO RISE/FALL
          DO NC=1,NCTMP
            DO M=1,MQWRSR(NS)
              CQWRSER(M,NS,NC)=0.
            ENDDO
          ENDDO
          DO M=1,MQWRSR(NS)
            READ(1,*,IOSTAT=ISO)TQWRSER(M,NS),QWRSER(M,NS)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QWRS.INP')
            TQWRSER(M,NS)=TQWRSER(M,NS)+TAQWRSR(NS)
            QWRSER(M,NS)=(RMULADJ*(QWRSER(M,NS)+ADDADJ))
          ENDDO
        ELSE
          ! *** FLOW WITH RISE/FALL
          DO M=1,MQWRSR(NS)
            READ(1,*,IOSTAT=ISO) TQWRSER(M,NS), QWRSER(M,NS),(CQWRSER(M,NS,NC),NC=1,NCTMP)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QWRS.INP')
            TQWRSER(M,NS) = TQWRSER(M,NS) + TAQWRSR(NS)
            QWRSER(M,NS) = (RMULADJ*(QWRSER(M,NS) + ADDADJ))
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)

    endif
    
    Call Broadcast_Array(MQWRSR,  master_id)
    Call Broadcast_Array(TCQWRSR, master_id)
    Call Broadcast_Array(TAQWRSR, master_id)
    Call Broadcast_Array(TQWRSER, master_id)
    Call Broadcast_Array(QWRSER,  master_id)
    Call Broadcast_Array(CQWRSER, master_id)
    
  ENDIF

  ! *** READ IN GROUNDWATER INFLOW/OUTFLOW AND CONCENTRATION TIME
  ! *** SERIES FROM THE FILE GWSER.INP
  IF( ISGWIT == 2 )THEN
    ! *** SKIP OVER TITLE AND AND HEADER LINES
    NCTMP = 3 + NDYM + NSED + NSND + NTOX

    ! *** READ IN GW CONCENTRATIONS
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING GWSER.INP'
      OPEN(1,FILE='gwser.inp',STATUS='UNKNOWN')

      STR=READSTR(1)
      READ(1,*)NGWSER
      IF( NGWSER > 0 )THEN
        DO NS=1,NGWSER
          READ(1,*,IOSTAT=ISO) MGWSER(NS), TCGWSER(NS), TAGWSER(NS), RMULADJ, ADDADJ, IGWSER(NS)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE GWSER.INP')
          DO M=1,MGWSER(NS)
            READ(1,*,IOSTAT=ISO)TGWSER(M,NS),GWSER(M,NS),(GWCSER(M,NS,NC),NC=1,NCTMP)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE GWSER.INP')
            TGWSER(M,NS) = TGWSER(M,NS)+TAGWSER(NS)
            GWSER(M,NS)  = RMULADJ*(GWSER(M,NS) + ADDADJ)
          ENDDO
        ENDDO
      ENDIF
      CLOSE(1)
    endif !***End on master

    Call Broadcast_Array(MGWSER, master_id)
    Call Broadcast_Array(TCGWSER, master_id)
    Call Broadcast_Array(TAGWSER, master_id)
    Call Broadcast_Array(IGWSER, master_id)
    Call Broadcast_Array(GWSER,  master_id)
    Call Broadcast_Array(TGWSER, master_id)

  ENDIF

  ! *** READ IN SPATIAL MAPS AND TIME SERIES FOR EXTERNAL SPECIFICATION OF
  ! *** PARTICULATE ORGANIC CARBON FOR USE IN TOXIC CONTAMINANT SORPTION
  ! *** DISSOLVED ORGANIC CARBON
  ! *** SKIP OVER TITLE AND AND HEADER LINES
  !        READ(1,*)
  !        READ(1,*)NOCSER
  !            READ(1,*,IOSTAT=ISO)MOCSER(NS),TCOCSER(NS),TAOCSER(NS),
  !              READ(1,*,IOSTAT=ISO)TOCSER(M,NS),DOCWSER(M,NS),
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SALINITY TIME SERIES
  ! *** FROM THE FILE SSER.INP
  IF( ISTRAN(1) >= 1 .AND. NCSER(1) >= 1 )THEN
    NC = 1
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING SSER.INP'
      OPEN(1,FILE='sser.inp',STATUS='UNKNOWN')
                                                                  
      ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
      STR=READSTR(1)
      NC = 1
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TACSER(NS,NC), RMULADJ, ADDADJ
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSAL(NS).TIM(M),CSERTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
            TSSAL(NS).TIM(M) = TSSAL(NS).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSSAL(NS).VAL(M,K) =(RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO) TSSAL(NS).TIM(M),(TSSAL(NS).VAL(M,K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SSER.INP')
            TSSAL(NS).TIM(M) = TSSAL(NS).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSSAL(NS).VAL(M,K) = RMULADJ*(TSSAL(NS).VAL(M,K) + ADDADJ)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    endif

    
    Call Broadcast_Array(MCSER,  master_id)
    Call Broadcast_Array(TCCSER, master_id)
    Call Broadcast_Array(TACSER, master_id)
    
    DO NS=1,NCSER(NC)
        Call Broadcast_Array(TSSAL(NS).VAL, master_id)
        Call Broadcast_Array(TSSAL(NS).TIM, master_id)
    End do
    
  ENDIF

  !========================================================================
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE TEMPERATURE TIME
  ! *** SERIES FROM THE FILE TSER.INP
  IF( ( ISTRAN(2) >= 1 .AND. NCSER(2) >= 1 ) .OR. ( ISTRAN(5) > 0 .AND. ITOXTEMP > 1 ) )THEN

    NC = 2
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING TSER.INP'
      OPEN(1,FILE='tser.inp',STATUS='UNKNOWN')
                                                                  
      ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
      STR=READSTR(1)
      
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TACSER(NS,NC), RMULADJ, ADDADJ
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSTEM(NS).TIM(M),CSERTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
            TSTEM(NS).TIM(M) = TSTEM(NS).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSTEM(NS).VAL(M,K) = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
           READ(1,*,IOSTAT=ISO) TSTEM(NS).TIM(M),(TSTEM(NS).VAL(M,K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TSER.INP')
            TSTEM(NS).TIM(M) = TSTEM(NS).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSTEM(NS).VAL(M,K) = RMULADJ*(TSTEM(NS).VAL(M,K) + ADDADJ)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    end if

   Call Broadcast_Array(MCSER,  master_id)
   Call Broadcast_Array(TCCSER, master_id)
   Call Broadcast_Array(TACSER, master_id)
      
    DO NS=1,NCSER(NC)
      Call Broadcast_Array(TSTEM(NS).TIM, master_id)
      Call Broadcast_Array(TSTEM(NS).VAL, master_id)
    ENDDO
  
  ENDIF

  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE DYE CONCENTRATION
  ! *** TIME SERIES FROM THE FILE DSER.INP
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE DYE CONCENTRATION
  ! *** TIME SERIES FROM THE FILE DSER.INP
  IF( ISTRAN(3) >= 1 .AND. NCSER(3) >= 1 )THEN
    NC = 3
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING DSER.INP'
      OPEN(1,FILE='dser.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TACSER(NS,NC), RMULADJ, ADDADJ
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSDYE(NS,1).TIM(M),CSERTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
            TSDYE(NS,1).TIM(M) = TSDYE(NS,1).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSDYE(NS,1).VAL(M,K) = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
            ENDDO
            DO MD=2,NDYE
              TSDYE(NS,MD).TIM(M) = TSDYE(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO) CSERTMP
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
              DO K=1,KC
                TSDYE(NS,MD).VAL(M,K) = (RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSDYE(NS,1).TIM(M),(TSDYE(NS,1).VAL(M,K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE DSER.INP')
            TSDYE(NS,1).TIM(M) = TSDYE(NS,1).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSDYE(NS,1).VAL(M,K) = RMULADJ*(TSDYE(NS,1).VAL(M,K) + ADDADJ)
            ENDDO
            DO MD=2,NDYE
              TSDYE(NS,MD).TIM(M) = TSDYE(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO) (TSDYE(NS,MD).VAL(M,K), K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
              DO K=1,KC
                TSDYE(NS,MD).VAL(M,K) = RMULADJ*(TSDYE(NS,MD).VAL(M,K) + ADDADJ)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    End if !*** end master
    
    DO NS=1,NCSER(NC)
      Call Broadcast_Scalar(MCSER(NS,NC),  master_id)
      Call Broadcast_Scalar(TCCSER(NS,NC), master_id)
      Call Broadcast_Scalar(TACSER(NS,NC), master_id)
      DO MD=1,NDYE
        Call Broadcast_Array(TSDYE(NS,MD).TIM, master_id)
        Call Broadcast_Array(TSDYE(NS,MD).VAL, master_id)
      ENDDO
    ENDDO
  
  ENDIF

  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SHELL FISH LARVAE
  ! *** TIME SERIES FROM THE FILE SFSER.INP
  IF(ISTRAN(4) >= 1 .AND. NCSER(4) >= 1 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING SFSER.INP'
      OPEN(1,FILE='sfser.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      NC = 4
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSFL(NS).TIM(M),CSERTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
            TSSFL(NS).TIM(M) = TSSFL(NS).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSSFL(NS).VAL(M,K)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSFL(NS).TIM(M),(TSSFL(NS).VAL(M,K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
            TSSFL(NS).TIM(M) = TSSFL(NS).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSSFL(NS).VAL(M,K)=RMULADJ*(TSSFL(NS).VAL(M,K)+ADDADJ)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    End if !*** end master
  ENDIF

  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE COHESIVE SEDIMENT
  ! *** CONCENTRATION TIME SERIES FROM THE FILE SDSER.INP
  IF( ISTRAN(6) >= 1 .AND. NSED > 0 )THEN
    NC = 6  !MSVSED(1)
    IF( NCSER(NC) >= 1 )THEN

      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SDSER.INP'
        OPEN(1,FILE='sdser.inp',STATUS='UNKNOWN')
                                                                  
        ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
        STR=READSTR(1)
        DO NS=1,NCSER(NC)
          READ(1,*,IOSTAT=ISO)ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TACSER(NS,NC), RMULADS(1), ADDADS(1)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
          IF( NSED>1 )THEN
            DO NT=2,NSED
              READ(1,*,IOSTAT=ISO) RMULADS(NT),ADDADS(NT)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
            ENDDO
          ENDIF
          IF( ISTYP == 1 )THEN
            READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
            DO M=1,MCSER(NS,NC)
              READ(1,*,IOSTAT=ISO) TSSED(NS,1).TIM(M),CSERTMP
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
              TSSED(NS,1).TIM(M) = TSSED(NS,1).TIM(M) + TACSER(NS,NC)
              DO K=1,KC
               TSSED(NS,1).VAL(M,K) = ( RMULADS(1)*( CSERTMP + ADDADS(1) ) )*WKQ(K)
              ENDDO
              DO NT=2,NSED
                !NTT=NT-1
                TSSED(NS,NT).TIM(M)=TSSED(NS,1).TIM(M)
                READ(1,*,IOSTAT=ISO) CSERTMP
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
                DO K=1,KC
                  TSSED(NS,NT).VAL(M,K) =( RMULADS(NT)*( CSERTMP + ADDADS(NT) ) )*WKQ(K)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO M=1,MCSER(NS,NC)
              READ(1,*,IOSTAT=ISO) TSSED(NS,1).TIM(M),(TSSED(NS,1).VAL(M,K),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
              TSSED(NS,1).TIM(M) = TSSED(NS,1).TIM(M) + TACSER(NS,NC)
              DO K=1,KC
                TSSED(NS,1).VAL(M,K) = RMULADS(1)*( TSSED(NS,1).VAL(M,K) + ADDADS(1) )
              ENDDO
              DO NT=2,NSED
                !NTT=NT-1
                TSSED(NS,NT).TIM(M)=TSSED(NS,1).TIM(M)
                READ(1,*,IOSTAT=ISO)(TSSED(NS,NT).VAL(M,K), K=1,KC)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SDSER.INP')
                DO K=1,KC
                  TSSED(NS,NT).VAL(M,K) = RMULADS(NT)*( TSSED(NS,NT).VAL(M,K) + ADDADS(NT) )
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      
      DO NS=1,NCSER(NC)
        Call Broadcast_Scalar(MCSER(NS,NC),       master_id)
        Call Broadcast_Scalar(TCCSER(NS,NC),      master_id)
        Call Broadcast_Scalar(TACSER(NS,NC),      master_id)
        DO NT=1,NSED
          Call Broadcast_Array(TSSED(NS,NT).TIM, master_id)
          Call Broadcast_Array(TSSED(NS,NT).VAL, master_id)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! *** CHECK SEDIMENT SERIES
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE NONCOHESIVE SEDIMENT
  ! *** CONCENTRATION TIME SERIES FROM THE FILE SNSER.INP
  IF( ISTRAN(7) >= 1 .AND. NSND > 0 )THEN
    NC = 7  ! MSVSND(1)
    IF( NCSER(NC) >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SNSER.INP'
        OPEN(1,FILE='snser.inp',STATUS='UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR=READSTR(1)
      
        DO NS=1,NCSER(NC)
          READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADS(1),ADDADS(1)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
          IF( NSND > 1 )THEN
            DO NT=2,NSND
              READ(1,*,IOSTAT=ISO)RMULADS(NT),ADDADS(NT)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
            ENDDO
          ENDIF
          IF( ISTYP == 1 )THEN
            READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
            DO M=1,MCSER(NS,NC)
              READ(1,*,IOSTAT=ISO) TSSND(NS,1).TIM(M),CSERTMP
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
              TSSND(NS,1).TIM(M) = TSSND(NS,1).TIM(M) + TACSER(NS,NC)
              DO K=1,KC
                TSSND(NS,1).VAL(M,K)=(RMULADS(1)*(CSERTMP+ADDADS(1)))*WKQ(K)
              ENDDO
              DO NT=2,NSND
                !NTT=NT-1
                TSSND(NS,NT).TIM(M)=TSSND(NS,1).TIM(M)
                READ(1,*,IOSTAT=ISO)CSERTMP
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
                DO K=1,KC
                  TSSND(NS,NT).VAL(M,K) =(RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO M=1,MCSER(NS,NC)
              READ(1,*,IOSTAT=ISO) TSSND(NS,1).TIM(M),(TSSND(NS,1).VAL(M,K),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
              TSSND(NS,1).TIM(M) = TSSND(NS,1).TIM(M) + TACSER(NS,NC)
              DO K=1,KC
                TSSND(NS,1).VAL(M,K) = RMULADS(1)*(TSSND(NS,1).VAL(M,K) + ADDADS(1))
              ENDDO
              DO NT=2,NSND
                !NTT=NT-1
                TSSND(NS,NT).TIM(M)=TSSND(NS,1).TIM(M)
                READ(1,*,IOSTAT=ISO)(TSSND(NS,NT).VAL(M,K), K=1,KC)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SNSER.INP')
                DO K=1,KC
                  TSSND(NS,NT).VAL(M,K) =RMULADS(NT)*(TSSND(NS,NT).VAL(M,K)+ADDADS(NT))
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      
      DO NS=1,NCSER(NC)
        Call Broadcast_Scalar(MCSER(NS,NC),       master_id)
        Call Broadcast_Scalar(TCCSER(NS,NC),      master_id)
        Call Broadcast_Scalar(TACSER(NS,NC),      master_id)
        DO NT=1,NSND
          Call Broadcast_Array(TSSND(NS,NT).TIM, master_id)
          Call Broadcast_Array(TSSND(NS,NT).VAL, master_id)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE TOXIC CONTAMINANT
  ! *** CONCENTRATION TIME SERIES FROM THE FILE TXSER.INP
  IF( ISTRAN(5) >= 1 .AND. NTOX > 0 )THEN
    NC = 5
    IF( NCSER(NC) >= 1 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING TXSER.INP'
        OPEN(1,FILE='txser.inp',STATUS='UNKNOWN')

        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR=READSTR(1)
        DO NS=1,NCSER(NC)
          READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
          IF( ISTYP == 1 )THEN
            READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
            DO M=1,MCSER(NS,NC)
              READ(1,*,IOSTAT=ISO) TSTOX(NS,1).TIM(M),CSERTMP
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
              TSTOX(NS,1).TIM(M) = TSTOX(NS,1).TIM(M) + TACSER(NS,NC)
              DO K=1,KC
                TSTOX(NS,1).VAL(M,K) = (RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
              ENDDO
              DO NT=2,NTOX
                !NTT=NT-1
                TSTOX(NS,NT).TIM(M)=TSTOX(NS,1).TIM(M)
                READ(1,*,IOSTAT=ISO) CSERTMP
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
                DO K=1,KC
                  TSTOX(NS,NT).VAL(M,K) = ( RMULADJ*( CSERTMP + ADDADJ) )*WKQ(K)
                ENDDO
              ENDDO
            ENDDO
          ELSE
            DO M=1,MCSER(NS,NC)
              READ(1,*,IOSTAT=ISO) TSTOX(NS,1).TIM(M),(TSTOX(NS,1).VAL(M,K), K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
              TSTOX(NS,1).TIM(M) = TSTOX(NS,1).TIM(M) + TACSER(NS,NC)
              DO K=1,KC
                TSTOX(NS,1).VAL(M,K) = RMULADJ*(TSTOX(NS,1).VAL(M,K) + ADDADJ)
              ENDDO
              DO NT=2,NTOX
                !NTT=NT-1
                TSTOX(NS,NT).TIM(M)=TSTOX(NS,1).TIM(M)
                READ(1,*,IOSTAT=ISO)(TSTOX(NS,NT).VAL(M,K), K=1,KC)
                IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')
                DO K=1,KC
                  TSTOX(NS,NT).VAL(M,K)=RMULADJ*(TSTOX(NS,NT).VAL(M,K) + ADDADJ)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      endif
      
      DO NS=1,NCSER(NC)
        Call Broadcast_Scalar(MCSER(NS,NC),       master_id)
        Call Broadcast_Scalar(TCCSER(NS,NC),      master_id)
        Call Broadcast_Scalar(TACSER(NS,NC),      master_id)
        DO NT=1,NTOX
          Call Broadcast_Array(TSTOX(NS,NT).TIM, master_id)
          Call Broadcast_Array(TSTOX(NS,NT).VAL, master_id)
        ENDDO
      ENDDO
      
    ENDIF

    ! *** TIME SERIES DRY DEPOSITION BEING USED
    IF( SUM(TOXDEP(:).ITXDRY) > 0 )THEN
      WRITE(*,'(A)')'READING TXDRY.INP'
      OPEN(1,FILE='txdry.inp',STATUS='UNKNOWN')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*,IOSTAT=ISO) ISTYP, TXDRYSER(1).NREC, T1, T2, RMULADJ, ADDADJ
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')

      ! *** IGNORE ISTYP, ALWAYS USE ONE VALUE
      READ(1,*)                                                     ! *** SKIP SPLITS
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXDRY.INP')

      DO M=1,TXDRYSER(1).NREC
        READ(1,*,IOSTAT=ISO) TXDRYSER(1).TIM(M), CSERTMP
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXDRY.INP')

        TXDRYSER(1).TIM(M)   = (TXDRYSER(1).TIM(M) + T2)*T1         ! *** TIMES ARE IN SECONDS
        TXDRYSER(1).VAL(M,1) = ((CSERTMP + ADDADJ)*RMULADJ)/86400.  ! *** MG/M2/SEC
        TOXDEP(1).ITDRY = 1                                         ! *** SET INITIAL TIME INDEX
        DO NT=2,NTOX
          READ(1,*,IOSTAT=ISO) CSERTMP
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXDRY.INP')
          TXDRYSER(1).VAL(M,NT) = (CSERTMP + ADDADJ)*RMULADJ
          TOXDEP(NT).ITDRY = 1                                     ! *** SET INITIAL TIME INDEX
        ENDDO
      ENDDO
      CLOSE(1)

    ENDIF

    ! *** TIME SERIES WET DEPOSITION BEING USED
    IF( SUM(TOXDEP(:).ITXWET) > 0 )THEN
      WRITE(*,'(A)')'READING TXWET.INP'
      OPEN(1,FILE='txwet.inp',STATUS='UNKNOWN')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*,IOSTAT=ISO) ISTYP, TXWETSER(1).NREC, T1, T2, RMULADJ, ADDADJ
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXSER.INP')

      ! *** IGNORE ISTYP, ALWAYS USE ONE VALUE
      READ(1,*)                                                     ! *** SKIP SPLITS
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXWET.INP')

      DO M=1,TXWETSER(1).NREC
        READ(1,*,IOSTAT=ISO) TXWETSER(1).TIM(M), CSERTMP
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXWET.INP')

        TXWETSER(1).TIM(M)   = (TXWETSER(1).TIM(M) + T2)*T1         ! *** TIMES ARE IN SECONDS
        TXWETSER(1).VAL(M,1) = (CSERTMP + ADDADJ)*RMULADJ           ! *** MG/M3
        TOXDEP(1).ITWET = 1                                         ! *** SET INITIAL TIME INDEX

        DO NT=2,NTOX
          READ(1,*,IOSTAT=ISO) CSERTMP
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE TXWET.INP')
          TXWETSER(1).VAL(M,NT) = (CSERTMP + ADDADJ)*RMULADJ
          TOXDEP(NT).ITWET = 1                                      ! *** SET INITIAL TIME INDEX
        ENDDO
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF

  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SHELL FISH LARVAE
  ! *** TIME SERIES FROM THE FILE SFSER.INP
  IF( ISTRAN(4) >= 1 .AND. NCSER(4) >= 1 )THEN
    WRITE(*,'(A)')'READING SFSER.INP'
    OPEN(1,FILE='sfser.inp',STATUS='UNKNOWN')

    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR=READSTR(1)
    NC=7
    DO NS=1,NCSER(NC)
      READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
        DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO) TSSFL(NS).TIM(M),CSERTMP
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
          TSSFL(NS).TIM(M) = TSSFL(NS).TIM(M) + TACSER(NS,NC)
          DO K=1,KC
            TSSFL(NS).VAL(M,K)=(RMULADJ*(CSERTMP + ADDADJ))*WKQ(K)
          ENDDO
        ENDDO
      ELSE
        DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO) TSSFL(NS).TIM(M),(TSSFL(NS).VAL(M,K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFSER.INP')
          TSSFL(NS).TIM(M) = TSSFL(NS).TIM(M) + TACSER(NS,NC)
          DO K=1,KC
            TSSFL(NS).VAL(M,K) = RMULADJ*(TSSFL(NS).VAL(M,K) + ADDADJ)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
  ! *** READ IN SHELLFISH FARM
  IF((ISFFARM > 0) .AND. (NSF > 0) )THEN
    CALL READ_SHELLFISH_JSON(LCM,KC)    
  ENDIF

  ! *** READ IN FREE SURFACE ELEVATION OR PRESSURE CONTROLLED FLOW
  ! *** SPECIFICATION FROM THE FILE QCTL.INP
  ! *** THE FLOW IS GIVEN BY:
  !         FREE SURFACE
  !         FREE SURFACE
  !        FLOW=0
  !      ELSE
  !        ENTER QCTL(M,K,NS) VS HDIFCTL(M,NS) TABLE WITH DELH TO GIVE
  IF( NQCTLT >= 1 )THEN
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING QCTL.INP'
      OPEN(1,FILE='qctl.inp',STATUS='UNKNOWN')
      IF( DEBUG )THEN
        OPEN(99,FILE=OUTDIR//'qctlck.inp',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE=OUTDIR//'qctlck.inp',STATUS='UNKNOWN')
      ENDIF

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      DO NS=1,NQCTLT
        READ(1,*, IOSTAT=ISO)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
        IF( DEBUG )THEN
          WRITE(99,991)NS
          WRITE(99,992)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
        ENDIF
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
        IF( ISTYP == 0 )THEN
          DO M=1,MQCTL(NS)
            READ(1,*,IOSTAT=ISO) HDIFCTL(M,NS),(QCTL(M,1,K,NS),K=1,KC)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
            DO K=1,KC
              QCTL(M,1,K,NS)=RMULADJ*(QCTL(M,1,K,NS)+ADDADJ)
            ENDDO
          ENDDO
        ENDIF
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
          DO M=1,MQCTL(NS)
            READ(1,*,IOSTAT=ISO) HDIFCTL(M,NS),QCTLTMP
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
            DO K=1,KC
              QCTL(M,1,K,NS)=RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
            ENDDO
          ENDDO
        ENDIF
        IF( ISTYP == 2 )THEN
          DO MD=1,MQCTL(NS)
            DO MU=1,MQCTL(NS)
              READ(1,*,IOSTAT=ISO) HDIFCTL(MU,NS),HDIFCTD(MD,NS),(QCTL(MU,MD,K,NS),K=1,KC)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
              DO K=1,KC
                QCTL(MU,MD,K,NS)=RMULADJ*(QCTL(MU,MD,K,NS)+ADDADJ)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF( ISTYP == 3 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
          DO MD=1,MQCTL(NS)
            DO MU=1,MQCTL(NS)
              READ(1,*,IOSTAT=ISO)HDIFCTL(MU,NS),HDIFCTD(MD,NS),QCTLTMP
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTL.INP')
              DO K=1,KC
                QCTL(MU,MD,K,NS)=RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF( DEBUG )THEN
          IF( ISTYP <= 1 )THEN
            DO M=1,MQCTL(NS)
              WRITE(99,993)M,HDIFCTL(M,NS),(QCTL(M,1,K,NS),K=1,KC)
            ENDDO
          ENDIF
          IF( ISTYP >= 2 )THEN
            DO MD=1,MQCTL(NS)
              DO MU=1,MQCTL(NS)
                WRITE(99,994)MU,MD,HDIFCTL(MU,NS),HDIFCTD(MD,NS),(QCTL(MU,MD,K,NS),K=1,KC)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      CLOSE(1)
      IF( DEBUG)CLOSE(99)
    endif

    Call Broadcast_Array(MQCTL   , master_id)
    Call Broadcast_Array(HCTLUA  , master_id)
    Call Broadcast_Array(HCTLUM  , master_id)
    Call Broadcast_Array(HCTLDA  , master_id)
    Call Broadcast_Array(HCTLDM  , master_id)
    Call Broadcast_Array(AQCTL   , master_id)
    Call Broadcast_Array(HDIFCTL , master_id)
    Call Broadcast_Array(HDIFCTD , master_id)
    Call Broadcast_Array(QCTL    , master_id)

  ENDIF

  ! *** READ CONTROL TIME-SERIES
  IF( NQCTLSER > 0 )THEN
    ALLOCATE(QCTLSER(NQCTLSER))
    
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING QCTLSER.INP'
      OPEN(1,FILE='qctlser.inp',STATUS='UNKNOWN')
      IF( DEBUG )THEN
        OPEN(99,FILE=OUTDIR//'qctlser.log',STATUS='REPLACE')
      ENDIF

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
    endif
    
    DO NS=1,NQCTLSER
      if( process_id == master_id )then
        READ(1,*, IOSTAT=ISO) QCTLSER(NS).ITYPE, NX, QCTLSER(NS).TMUL, QCTLSER(NS).TADD, (ISQCTRL(M),M=1,6)
        QCTLSER(NS).PARAM = 0
        QCTLSER(NS).COUNT = NX
        IVAL = 1
        DO M=1,6
          IF(ISQCTRL(M) > 0) QCTLSER(NS).PARAM = QCTLSER(NS).PARAM + IVAL
          IVAL = 2*IVAL
        ENDDO
        IF( DEBUG )THEN
          WRITE(99,991) NS
          WRITE(99,992) QCTLSER(NS).ITYPE, NX, QCTLSER(NS).TMUL, QCTLSER(NS).TADD, QCTLSER(NS).PARAM
        ENDIF
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTLSER.INP')
      endif
      Call Broadcast_Scalar(QCTLSER(NS).ITYPE   , master_id)
      Call Broadcast_Scalar(QCTLSER(NS).TMUL    , master_id)
      Call Broadcast_Scalar(QCTLSER(NS).TADD    , master_id)
      Call Broadcast_Scalar(QCTLSER(NS).PARAM   , master_id)
      Call Broadcast_Scalar(QCTLSER(NS).COUNT   , master_id)
      
      NX = QCTLSER(NS).COUNT
      ALLOCATE(QCTLSER(NS).TIME(NX))
      ALLOCATE(QCTLSER(NS).ID(NX))
      ALLOCATE(QCTLSER(NS).HEIGHT(NX))
      ALLOCATE(QCTLSER(NS).WIDTH(NX))
      ALLOCATE(QCTLSER(NS).SILL(NX))
      !ALLOCATE(QCTLSER(NS).NUM(NX))
      ALLOCATE(QCTLSER(NS).FLOW(NX))
      !ALLOCATE(QCTLSER(NS).RATE(NX))
      
      if( process_id == master_id )then
        DO M=1,NX
          READ(1,*,IOSTAT=ISO) QCTLTMP, QCTLSER(NS).HEIGHT(M), QCTLSER(NS).WIDTH(M), QCTLSER(NS).SILL(M), &
                               NDUM,    QCTLSER(NS).FLOW(M),  QCTLSER(NS).ID(M)    !,QCTLSER(NS).RATE(M)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTLSER.INP')
          QCTLSER(NS).TIME(M) = QCTLTMP*QCTLSER(NS).TMUL + QCTLSER(NS).TADD
        ENDDO
      endif
      Call Broadcast_Array(QCTLSER(NS).TIME     , master_id)
      Call Broadcast_Array(QCTLSER(NS).ID       , master_id)
      Call Broadcast_Array(QCTLSER(NS).HEIGHT   , master_id)
      Call Broadcast_Array(QCTLSER(NS).WIDTH    , master_id)
      Call Broadcast_Array(QCTLSER(NS).SILL     , master_id)
      !Call Broadcast_Array(QCTLSER(NS).NUM      , master_id)  ! NOT USED
      Call Broadcast_Array(QCTLSER(NS).FLOW     , master_id)
    ENDDO
    if( process_id == master_id )then
      CLOSE(1)
      IF( DEBUG ) CLOSE(99)
    endif
    
    
  ENDIF

  ! *** READ CONTROL RULES
  IF( NQCTRULES > 0 )THEN
    ALLOCATE(QCTRULES(NQCTRULES))

    if( process_id == master_id )then
        WRITE(*,'(A)')'READING QCTRULES.INP'
        OPEN(1,FILE='qctrules.inp',STATUS='UNKNOWN')

        OPEN(99,FILE=OUTDIR//'QCTRULES.log',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE=OUTDIR//'QCTRULES.log',STATUS='UNKNOWN')
        
        ! *** SKIP OVER TITLE AND AND HEADER LINES
        DO NS=1,NQCTRULES
          STR=READSTR(1)       ! *** Skip all comments to the type of control flags

          READ(1,*, IOSTAT=ISO) NX,NX,(ISQCTRL(M),M=1,6)  !,QCTRULES(NS).PARAM
          QCTRULES(NS).PARAM = 0
          IVAL = 1
          DO M=1,6
            IF(ISQCTRL(M) > 0) QCTRULES(NS).PARAM = QCTRULES(NS).PARAM + IVAL
            IVAL = 2*IVAL
          ENDDO
          WRITE(99,*) NS,NX,QCTRULES(NS).PARAM
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTRULES.INP')
        
          QCTRULES(NS).NTRGON = 0
          QCTRULES(NS).NTROFF = 0
        
          STR=READSTR(1)       ! *** Skip intermediate comments

          ALLOCATE(RULES(NX,9))
          ALLOCATE(IDX(NX))
          DO M=1,NX
            READ(1,*,IOSTAT=ISO) RULES(M,1), ISTYP, RULES(M,4), RULES(M,5), RULES(M,6), NDUM, RULES(M,8), RULES(M,9), IVAL
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE QCTRULES.INP')
            IDX(M) = M
            RULES(M,2) = ISTYP
            RULES(M,3) = IVAL
            RULES(M,7) = NDUM
            IF( ISTYP == 0 )THEN
              QCTRULES(NS).NTROFF = QCTRULES(NS).NTROFF + 1
            ELSE
              QCTRULES(NS).NTRGON = QCTRULES(NS).NTRGON + 1
            ENDIF
          ENDDO
          ! SORT TRIGGER LEVELS DESCENDING
          DO M=1,NX-1
            DO L =M+1,NX
              IF(RULES(IDX(L),1) > RULES(IDX(M),1) )THEN
                IVAL = IDX(L)
                IDX(L) = IDX(M)
                IDX(M) = IVAL
              ENDIF
            ENDDO
          ENDDO

          ! *** Setup list of triggers to open the gates/pumps
          ALLOCATE(QCTRULES(NS).TRGON(QCTRULES(NS).NTRGON))
          L = 0
          DO M=1,NX
            IF( RULES(M,2) == 1. )THEN
              L = L + 1
              QCTRULES(NS).TRGON(L).LEVEL = RULES(M,1)             ! TRIGGER LEVEL
              QCTRULES(NS).TRGON(L).STATE = INT(RULES(M,2))        ! 0 = OFF, 1 = ON
              QCTRULES(NS).TRGON(L).ID = INT(RULES(M,3))           ! INDEX OF RATING CURVE/MASK FOR CONTROL PARAMETERS
              QCTRULES(NS).TRGON(L).HEIGHT = RULES(M,4)            ! OPENING HEIGHT (M), FOR UPWARD OPENING
              QCTRULES(NS).TRGON(L).WIDTH = RULES(M,5)             ! OPENING WIDTH (M), FOR SIDEWARD OPENING
              QCTRULES(NS).TRGON(L).SILL = RULES(M,6)              ! SILL LEVEL CHANGE (M), FOR DOWNWARD OPENING
              !QCTRULES(NS).TRGON(L).UNITS = INT(RULES(M,7))       ! NUMBER OF GATES, PUMP UNITS   not used
              QCTRULES(NS).TRGON(L).FLOW = RULES(M,8)              ! FLOW RATE (M3/S), FOR PUMPS    
              QCTRULES(NS).TRGON(L).RATE = RULES(M,9)              ! RATE OF OPENING/CLOSING
              
              WRITE(99,*) L,QCTRULES(NS).TRGON(L).LEVEL, QCTRULES(NS).TRGON(L).STATE, &
                            QCTRULES(NS).TRGON(L).ID,    QCTRULES(NS).TRGON(L).HEIGHT, QCTRULES(NS).TRGON(L).WIDTH, &
                            QCTRULES(NS).TRGON(L).SILL,  QCTRULES(NS).TRGON(L).FLOW,   QCTRULES(NS).TRGON(L).RATE
                          
            ENDIF
          ENDDO

          ! *** Setup list of triggers to close the gates/pumps
          ALLOCATE(QCTRULES(NS).TROFF(QCTRULES(NS).NTROFF))
          L = 0
          DO M=NX,1,-1
            IF( RULES(M,2) == 0. )THEN
              L = L + 1
              QCTRULES(NS).TROFF(L).LEVEL = RULES(M,1)             ! TRIGGER LEVEL
              QCTRULES(NS).TROFF(L).STATE = INT(RULES(M,2))        ! 0 = OFF, 1 = ON
              QCTRULES(NS).TROFF(L).ID = INT(RULES(M,3))           ! INDEX OF RATING CURVE/MASK FOR CONTROL PARAMETERS
              QCTRULES(NS).TROFF(L).HEIGHT = RULES(M,4)            ! OPENING HEIGHT (M), FOR UPWARD OPENING
              QCTRULES(NS).TROFF(L).WIDTH = RULES(M,5)             ! OPENING WIDTH (M), FOR SIDEWARD OPENING
              QCTRULES(NS).TROFF(L).SILL = RULES(M,6)              ! SILL LEVEL CHANGE (M), FOR DOWNWARD OPENING
              !QCTRULES(NS).TROFF(L).UNITS = INT(RULES(M,7))       ! NUMBER OF GATES, PUMP UNITS   not used
              QCTRULES(NS).TROFF(L).FLOW = RULES(M,8)              ! FLOW RATE (M3/S), FOR PUMPS
              QCTRULES(NS).TROFF(L).RATE = RULES(M,9)              ! RATE OF OPENING/CLOSING
              WRITE(99,*) L,QCTRULES(NS).TROFF(L).LEVEL, QCTRULES(NS).TROFF(L).STATE, &
                            QCTRULES(NS).TROFF(L).ID,    QCTRULES(NS).TROFF(L).HEIGHT, QCTRULES(NS).TROFF(L).WIDTH, &
                            QCTRULES(NS).TROFF(L).SILL,  QCTRULES(NS).TROFF(L).FLOW,   QCTRULES(NS).TROFF(L).RATE
            ENDIF
          ENDDO
          DEALLOCATE(RULES,IDX)
        ENDDO
        CLOSE(1)
        CLOSE(99)
    end if ! *** end on master

    
    DO NS=1,NQCTRULES
      Call Broadcast_Scalar(QCTRULES(NS).PARAM  , master_id)

      Call Broadcast_Scalar(QCTRULES(NS).NTRGON  , master_id)
      ! *** only allocate on other processes besides master
      if( process_id /= master_id )then
        ALLOCATE(QCTRULES(NS).TRGON(QCTRULES(NS).NTRGON))
      end if
      
      DO M=1,QCTRULES(NS).NTRGON
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).LEVEL  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).STATE  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).ID     , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).HEIGHT , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).WIDTH  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).SILL   , master_id)
        !Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).UNITS  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).FLOW   , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TRGON(M).RATE   , master_id)
      ENDDO
        
      Call Broadcast_Scalar(QCTRULES(NS).NTROFF   , master_id)
      ! *** only allocate on other processes besides master
      if(process_id /= master_id )then
          ALLOCATE(QCTRULES(NS).TROFF(QCTRULES(NS).NTROFF))
      end if
      
      DO M=1,QCTRULES(NS).NTROFF
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).LEVEL  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).STATE  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).ID     , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).HEIGHT , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).WIDTH  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).SILL   , master_id)
        !Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).UNITS  , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).FLOW   , master_id)
        Call Broadcast_Scalar(QCTRULES(NS).TROFF(M).RATE   , master_id)
      ENDDO
    ENDDO
  ENDIF

  ! *** CHECK CONTROL DATA
  IF( NQCTL > 0 .AND. (NQCTLSER > 0 .OR. NQCTRULES > 0) )THEN
    DO L=1,NQCTL
      IVAL = HSCTL_GL(L).ID
      IF( HSCTL_GL(L).ITYPE == 1 )THEN
        ! *** STRUCTURE IS CONTROLLED BY TIME-SERIES
        IF( IVAL < 1 .OR. IVAL > NQCTLSER )THEN
          CALL STOPP('DATA ERROR: TIME-SERIES INDEX FOR HYDRAULIC STRUCTURE CONTROL IS INVALID!')
        ENDIF
        !IF (QCTLSER(IVAL).TIME(1) > TBEGIN .OR. QCTLSER(IVAL).TIME(QCTLSER(IVAL).COUNT) < TEND) CALL STOPP('DATA ERROR'
      ELSEIF (HSCTL_GL(L).ITYPE == 2 .OR. HSCTL_GL(L).ITYPE == 3 )THEN
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
        IF( IVAL < 1 .OR. IVAL > NQCTRULES )THEN
          CALL STOPP('DATA ERROR: RULES INDEX FOR HYDRAULIC STRUCTURE CONTROL IS INVALID!')
        ENDIF
        IF( (QCTRULES(IVAL).PARAM .AND. 32) == 32 )THEN
          ! *** SET OF LOOKUP TABLE
          DO I=1,QCTRULES(IVAL).NTRGON
            ID = QCTRULES(IVAL).TRGON(I).ID
            IF( ID < 1 .OR. ID > NQCTLT) CALL STOPP('DATA ERROR: INVALID INDEX FOR LOOKUP TABLE DEFINED IN CONTROL RULES')
          ENDDO
          DO I=1,QCTRULES(IVAL).NTROFF
            ID = QCTRULES(IVAL).TROFF(I).ID
            IF( ID < 1 .OR. ID > NQCTLT) CALL STOPP('DATA ERROR: INVALID INDEX FOR LOOKUP TABLE DEFINED IN CONTROL RULES')
          ENDDO
        ELSEIF (QCTRULES(IVAL).PARAM > 0 )THEN
          ! *** GATE OPENING
        ELSE
          ! *** PUMP FLOWS
        ENDIF
        IF( HSCTL_GL(L).ITYPE == 2 )THEN
          ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
          !LU = LIJ(HSCTL_GL(NCTL).IREFUP, HSCTL_GL(NCTL).JREFUP)
        ELSEIF (HSCTL_GL(L).ITYPE == 3 )THEN
          ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
          !LU = LIJ(HSCTL_GL(NCTL).IREFUP, HSCTL_GL(NCTL).JREFUP)
          !LD = LIJ(HSCTL_GL(NCTL).IREFDN, HSCTL_GL(NCTL).JREFDN)
        ENDIF
      ELSE
        ! *** STRUCTURE IS UNCONTROLLED
      ENDIF
    ENDDO
  ENDIF
  
  IF( NQWR > 0 .AND. NQCTRULES > 0 )THEN
    DO L=1,NQWR
      IVAL = WRCTL_GL(L).ID
      IF( WRCTL_GL(L).ITYPE == 1 .OR. WRCTL_GL(L).ITYPE == 2 )THEN
        ! *** W/R IS CONTROLLED BY OPERATION RULES
        IF( IVAL < 1 .OR. IVAL > NQCTRULES )THEN
          CALL STOPP('DATA ERROR: RULES INDEX FOR WITHDRAWAL/RETURN CONTROL IS INVALID!')
        ENDIF
      ENDIF
    ENDDO
  ENDIF

991 FORMAT(/,'CONTROL TABLE NS =',I5,/)
992 FORMAT(2I5,7F10.4)
993 FORMAT(I5,11F10.4)
994 FORMAT(2I5,11F10.4)
1001 FORMAT(/,'READ ERROR FROM FILE EFDC.INP ON CARD ',A3/)
1002 FORMAT(/,'INPUT ECHO NCARD = ',A/)

  DO L=1,LC
    PATMT(L)=1000.
    TATMT(L)=0.
    RAINT(L)=0.
    EVAPT(L)=0.
    SOLSWRT(L)=1.  ! *** Address SUNDAY.INP Option
    CLOUDT(L)=0.
    RHAT(L)=0.
    VPAT(L)=0.
    CLEVAP(L)=0.
    CCNHTT(L)=0.
  ENDDO
  
  ! *** DEACTIVATE THE ATMOSPHERIC FILE IF TEMPERATURE IS NOT SIMUATED
  !IF( ISTRAN(2) == 0 )NASER = 0  2015-06-22 DEPRECATED TO ALLOW EVAP AND RAIN
  LDAYLIGHT = .FALSE.
  LRAIN = .FALSE.
  IF( NASER > 0 )THEN
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING ASER.INP'
      OPEN(1,FILE='aser.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      CALL SKIPCOM(1,'*')
      READ(1,'(A)')TEXT
      IASERVER=PARSE_REAL(TEXT)*1000.
      STR=READSTR(1)

      IF( IASERVER < 7300 )THEN
        IEVAP = -1
      ENDIF

      WRITE(*,'(A,I5)')'  NUMBER OF ATMOSPHERIC SERIES=',NASER
      WRITE(*,'(A,L3)')'  COMPUTESOLRAD=',COMPUTESOLRAD
      WRITE(*,'(A,F10.2)')'  DS_LONG=',DS_LONG
      WRITE(*,'(A,F10.2)')'  DS_LAT=',DS_LAT

      DO NS=1,NASER
        READ(1,*,IOSTAT=ISO) M,TCASER(NS),TAASER(NS),IRELH(NS),RAINCVT,EVAPCVT,SOLRCVT,CLDCVT
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE ASER.INP')
        IF( NS == 1 .AND. IEVAP == -1 )THEN
          IEVAP = 1                     ! *** USE ASER DATA
          IF( EVAPCVT < 0 ) IEVAP = 2   ! *** LEGACY INPUT TO COMPUTE EVAPORATION USING EFDC ORIGINAL APPROACH
        ELSEIF( IEVAP == 0 )THEN
          ! *** IGNORE EVAPORATION BUT ALWAYS INCLUDE RAINFALL FROM ASER
          EVAPCVT = 0.
        ENDIF
        WRITE(*,'(A,I5,I10,F6.1)')'  SERIES, NPTS =',NS,TSATM(NS).NREC

        ! *** These parameters are read in for every series but only the last is actually used
        IF( IASERVER < 7300 )THEN
          READ(1,*,IOSTAT=ISO) IASWRAD,REVC,RCHC,SWRATNF,SWRATNS,FSWRATF,TEMTHKO,TEMBO,HTBED1,HTBED2
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE ASER.INP')
        ENDIF

        DO M=1,TSATM(NS).NREC
          READ(1,*,IOSTAT=ISO) TSATM(NS).TIM(M),(TSATM(NS).VAL(M,I),I=1,7)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE ASER.INP')
        ENDDO
        
        ! *** APPLY CONVERSIONS
        DO M=1,TSATM(NS).NREC
          TSATM(NS).TIM(M)  = TSATM(NS).TIM(M)+TAASER(NS)
          TSATM(NS).VAL(M,4)= RAINCVT*TSATM(NS).VAL(M,4)
          TSATM(NS).VAL(M,5)= EVAPCVT*TSATM(NS).VAL(M,5)
          TSATM(NS).VAL(M,6)= SOLRCVT*TSATM(NS).VAL(M,6)
          TSATM(NS).VAL(M,7)=  CLDCVT*TSATM(NS).VAL(M,7)
        ENDDO
      ENDDO
      CLOSE(1)
    endif

    call Broadcast_Scalar(IASERVER, master_id)
    
    if( IASERVER < 7300 )THEN
      Call Broadcast_Scalar(IASWRAD, master_id)
      Call Broadcast_Scalar(REVC,    master_id)
      Call Broadcast_Scalar(RCHC,    master_id)
      Call Broadcast_Scalar(SWRATNF, master_id)
      Call Broadcast_Scalar(SWRATNS, master_id)
      Call Broadcast_Scalar(FSWRATF, master_id)
      Call Broadcast_Scalar(TEMTHKO, master_id)
      Call Broadcast_Scalar(TEMBO,   master_id)
      Call Broadcast_Scalar(HTBED1,  master_id)
      Call Broadcast_Scalar(HTBED2,  master_id)
    end if
    
    DO NS=1,NASER
      Call Broadcast_Scalar(TCASER(NS)     , master_id)
      Call Broadcast_Scalar(TAASER(NS)     , master_id)
      Call Broadcast_Scalar(IRELH(NS)      , master_id)
      Call Broadcast_Scalar(TSATM(NS).NREC , master_id)
      Call Broadcast_Array(TSATM(NS).TIM   , master_id)
      Call Broadcast_Array(TSATM(NS).VAL   , master_id)
    ENDDO
      
    ! *** RHO    = 1000.0  Density (kg / m^3)
    ! *** CP     = 4179.0  Specific Heat (J / kg / degC)
    ! *** 0.2393E-6 = 1/RHO/CP
  ENDIF

  IF( NASER > 1 )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, NASER, 0.0)
    
    IF( process_id == master_id )THEN
      WRITE(*,'(A)')'  READING ATMMAP.INP'
      OPEN(1,FILE='atmmap.inp',STATUS='UNKNOWN')
      STR = READSTR(1)
      
      READ(1,*) i
      IF( I /= NATMMAP ) CALL STOPP('BAD NATMMAP')      ! *** BROADCAST IN SCAN_EFDC
    ENDIF    

    DO NA=1,NATMMAP
      if( process_id == master_id )THEN
          READ(1,*) TATMMAPBEG(NA), TATMMAPEND(NA)
          STR=READSTR(1)
          
          DO L=2,LA_Global
              READ(1,*) ID,JD,(atmwht_temp(nn),NN=1,NASER)
              
              LG = LIJ_GLOBAL(ID,JD)
              R2D_Global(LG,:) = atmwht_temp(:)
          ENDDO
      end if

      Call Broadcast_Array(R2D_Global, master_id)
      
      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
           ATMWHT(:,L,NA) = R2D_Global(LG,:)
        ENDIF
      ENDDO

    ENDDO
    
    ! *** send to all processes
    Call Broadcast_Array(TATMMAPBEG, master_id)
    Call Broadcast_Array(TATMMAPEND, master_id)
    
    if( process_id == master_id )THEN
      CLOSE(1)
    endif
    
    DEALLOCATE(R2d_Global)
  ENDIF
  
  ! *** BED TEMPERATURE
  TEMB(1) = ABS(TEMBO)
  TEMB(LC) = TEMB(1)

  TEMB1(1)  = TEMB(1)
  TEMB1(LC) = TEMB(1)
  IF( ISRESTI == 0 )THEN
    DO L=2,LA
      TEMB(L)  = TEMB(1)
      TEMB1(L) = TEMB(1)
    ENDDO
  ENDIF

  ! *** READ IN ABOVE WATER SURFACE WIND TIME SERIES FROM THE
  ! *** FILE WSER.INP
  DO L=2,LA
    WINDST(L)=0.
    TSX(L)=0.
    TSY(L)=0.
  ENDDO

  IF( NWSER > 0 )THEN
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING WSER.INP'
      OPEN(1,FILE='wser.inp',STATUS='UNKNOWN')

      ! *** PERIOD TO TURN OFF WIND SHEAR (DEPRECATED IN EE7.3 SINCE ADDITION OF ISICE)
      WRITE(*,'(A,I5)')  '  NUMBER OF WIND SERIES=',NWSER

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      CALL SKIPCOM(1,'*')

      ! *** LOOP OVER EACH TIME SERIES
      STR=READSTR(1)
      
      DO NS=1,NWSER
        READ(1,*,IOSTAT=ISO) M, TCWSER(NS), TAWSER(NS), WINDSCT, ISWDINT(NS), WINDH(NS)
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE WSER.INP')

        IF( WINDH(NS) <= 0.1) WINDH(NS) = 2.0     ! *** FIXING THE WIND SPEED MEASUREMENT HEIGHT (I.E. USE THE WIND SPEED AS ENTERED)
        WRITE(*,'(A,I5,I10,F6.1)')'  SERIES, NPTS, ANEMOMETER HEIGHT (m)=',NS,TSWND(NS).NREC,WINDH(NS)

        DO M=1,TSWND(NS).NREC
          READ(1,*,IOSTAT=ISO) TSWND(NS).TIM(M), (TSWND(NS).VAL(M,I),I=1,2) !TWSER(M,NS),WINDS(M,NS),WINDD(M,NS)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE WSER.INP')
        ENDDO

        ! *** APPLY FACTORS AND OFFSETS AND COMPLETE TSWND DATA STRUCTURE
        DO M=1,TSWND(NS).NREC
          TSWND(NS).TIM(M)=TSWND(NS).TIM(M)+TAWSER(NS)
        ENDDO
        IF( ISWDINT(NS) <= 1 )THEN
          DO M=1,TSWND(NS).NREC
            TSWND(NS).VAL(M,1)=WINDSCT*TSWND(NS).VAL(M,1)
          ENDDO
        ENDIF
        IF( ISWDINT(NS) == 1 )THEN
          DO M=1,TSWND(NS).NREC
            IF( TSWND(NS).VAL(M,2) <= 180. )THEN
              TSWND(NS).VAL(M,2)=TSWND(NS).VAL(M,2)+180.
              IF( TSWND(NS).VAL(M,2) == 360.) TSWND(NS).VAL(M,2)=0.
            ELSE
              TSWND(NS).VAL(M,2)=TSWND(NS).VAL(M,2)-180.
              IF( TSWND(NS).VAL(M,2) == 360.) TSWND(NS).VAL(M,2)=0.
            ENDIF
          ENDDO
        ENDIF
        IF( ISWDINT(NS) == 2 )THEN
          DO M=1,TSWND(NS).NREC
            TSWND(NS).VAL(M,1)=WINDSCT*TSWND(NS).VAL(M,1)
            TSWND(NS).VAL(M,2)=WINDSCT*TSWND(NS).VAL(M,2)
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    endif
    
    DO NS=1,NWSER
      Call Broadcast_Scalar(TSWND(NS).NREC, master_id)
      Call Broadcast_Scalar(TCWSER(NS)   , master_id)
      Call Broadcast_Scalar(TAWSER(NS)   , master_id)
      Call Broadcast_Scalar(ISWDINT(NS)  , master_id)
      Call Broadcast_Scalar(WINDH(NS)    , master_id)
      Call Broadcast_Array(TSWND(NS).TIM , master_id)
      Call Broadcast_Array(TSWND(NS).VAL , master_id)
    ENDDO

  ENDIF

  IF( NWSER > 1 )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, NWSER, 0.0)

    if( process_id == master_id )THEN
      WRITE(*,'(A)')'  READING WNDMAP.INP'
      OPEN(1,FILE='wndmap.inp',STATUS='UNKNOWN', SHARED)
      
      STR=READSTR(1)
      READ(1,*) i
      IF( I /= NWNDMAP ) CALL STOPP('BAD NWNDMAP')      ! *** BROADCAST IN SCAN_EFDC
    endif

    DO NW=1,NWNDMAP

      IF( process_id == master_id )THEN
        READ(1,*) TWNDMAPBEG(NW), TWNDMAPEND(NW)
        STR=READSTR(1)

        DO L=2,LA_Global
          READ(1,*) ID,JD,(WNDWHT_TEMP(nn),NN=1,NWSER)
          
          LG = LIJ_GLOBAL(ID,JD)
          R2D_Global(LG,:) = WNDWHT_TEMP(:)
        ENDDO
      ENDIF
      
      Call Broadcast_Scalar(TWNDMAPBEG(NW), master_id)
      Call Broadcast_Scalar(TWNDMAPEND(NW), master_id)
      Call Broadcast_Array(R2D_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
           WNDWHT(:,L,NW) = R2D_Global(LG,:)
        ENDIF
      ENDDO
    ENDDO

    if( process_id == master_id )THEN
      CLOSE(1)
    endif
    DEALLOCATE(R2d_Global)
  ENDIF

    ! **  READ IN EXTERNALLY SPECIFIED ICE COVER INFORMATION
    ! **  FROM THE FILE ISER.INP
    !----------------------------------------------------------------------C
    IF( ISICE == 1 .AND. NISER >= 1 )THEN
        if( process_id == master_id )THEN
            WRITE(*,'(A)')'READING ISER.INP'
            OPEN(1,FILE='iser.inp')
            STR=READSTR(1)
            
            DO NS=1,NISER
              RICECOVT(NS)=0.
              RICETHKT(NS)=0.
              MITLAST(NS) =2
            
              READ(1,*,IOSTAT=ISO) M,TCISER(NS),TAISER(NS),RMULADJCOV
              IF( ISO > 0 ) CALL STOPP('ISER.INP: READING ERROR')
            
              DO M=1,TSICE(NS).NREC
                ! *** TSICE(NS).VAL(M,1) is the fraction of ice coverage.  
                ! *** Ice thickness (RICETHKS from ICECOVER.INP) was a legacy approach that was not used in EFDC.  Removed in EFDC+
                READ(1,*,IOSTAT=ISO) TSICE(NS).TIM(M),TSICE(NS).VAL(M,1)
            
                IF( TSICE(NS).VAL(M,1) > 1.0 )THEN
                  TSICE(NS).VAL(M,2) = 1.0
                ELSE
                  TSICE(NS).VAL(M,2) = 0.0
                ENDIF
                IF( ISO > 0 ) CALL STOPP('ISER.INP: READING ERROR')
              ENDDO
            
              DO M=1,TSICE(NS).NREC
                TSICE(NS).TIM(M)   = TSICE(NS).TIM(M) + TAISER(NS)
                TSICE(NS).VAL(M,1) = RMULADJCOV*TSICE(NS).VAL(M,1)
              ENDDO
            
            ENDDO
            
            CLOSE(1)
        end if
        
        ! *** send to all other processes
        Call Broadcast_Array(MITLAST, master_id)
        Call Broadcast_Array(TCISER,  master_id)
        Call Broadcast_Array(TAISER,  master_id)
        
        DO NS=1, NISER
          Call Broadcast_Array(TSICE(NS).TIM, master_id)
          Call Broadcast_Array(TSICE(NS).VAL, master_id)
        END DO
        
    ELSEIF( ISICE == 2 )THEN
      if( process_id == master_id )THEN
        WRITE(*,'(A)')'READING ISTAT.INP'
        OPEN(1,FILE='istat.inp')
        STR=READSTR(1)
        NS = 1    ! ** ALWAYS NISER = 1
        MITLAST(NS) = 2
        
        READ(1,*,IOSTAT=ISO) M,TCISER(NS),TAISER(NS)
        IF( ISO > 0 ) CALL STOPP('ISTAT.INP: READING ERROR')
        
        DO M=1,TSICE(NS).NREC
          READ(1,*,IOSTAT=ISO) TSICE(NS).TIM(M), TSICE(NS).VAL(M,1)  ! *** Fraction of ice cover: 0 to 1
          IF( ISO > 0 ) CALL STOPP('ISTAT.INP: READING ERROR')
          TSICE(NS).TIM(M) = TSICE(NS).TIM(M) + TAISER(NS)
          if( TSICE(NS).VAL(M,1) > 1.0 ) TSICE(NS).VAL(M,1) = 1.0
        ENDDO
      END IF
      
      ! *** send to all other processes
      Call Broadcast_Array(MITLAST, master_id)
      Call Broadcast_Array(TCISER,  master_id)
      Call Broadcast_Array(TAISER,  master_id)
      
      DO NS=1, NISER
        Call Broadcast_Array(TSICE(NS).TIM, master_id)
        Call Broadcast_Array(TSICE(NS).VAL, master_id)
      end do
      
    ENDIF

    !----------------------------------------------------------------------C
    IF( ISICE == 1 .AND. NISER > 1 )THEN
      ! *** Read ice time series weighting
      if(process_id == master_id )then
          WRITE(*,'(A)')'READING ICEMAP.INP'
          OPEN(1,FILE='icemap.inp')
          STR=READSTR(1)
          READ(1,*) NICEMAP
        
          DO NI=1,NICEMAP
            READ(1,*)TICEMAPBEG(NI),TICEMAPEND(NI)
            STR=READSTR(1)
            DO L=2,LA_Global
              READ(1,*) ID,JD,(RICEWHT_global(NI,L,N),N=1,NISER)
            ENDDO
          ENDDO
      end if
        
      Call Broadcast_Array(TICEMAPBEG, master_id)
      Call Broadcast_Array(TICEMAPEND, master_id)
      Call Broadcast_Array(RICEWHT_global, master_id)
        
      ! *** Map to local 
      DO LG = 2, LA_GLOBAL
        L = MAP2LOCAL(LG).LL
        IF( L > 1 )THEN
          RICEWHT(:,L,:) = RICEWHT_GLOBAL(:,LG,:)
        END IF
      END DO
        
      IF( process_id == master_id )THEN
        CLOSE(1)
      END IF
    ENDIF

    ! *** End of ICE *******************************

    ! *** READ IN SHELL FISH LARAVE BEHAVIOR DATA
    ! *** FROM THE FILE SFBSER.INP
    IF( ISTRAN(4) >= 1 )THEN
    
      IF( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SFBSER.INP'
        OPEN(1,FILE='sfbser.inp',STATUS='UNKNOWN')
          
        ! *** SKIP OVER TITLE AND AND HEADER LINES
        STR=READSTR(1)
          
        READ(1,*,IOSTAT=ISO) MSFSER,TCSFSER,TASFSER,TSRSF,TSSSF,ISSFLDN,ISSFLFE,SFLKILL
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFBSER.INP')
        DO M=1,MSFSER
          READ(1,*,IOSTAT=ISO) TSFSER(M),RKDSFL(M),WSFLST(M),WSFLSM(M),DSFLMN(M),DSFLMX(M),SFNTBE(M),SFATBT(M)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE SFBSER.INP')
        ENDDO
        CLOSE(1)
      END IF
      
      Call Broadcast_Scalar(MSFSER,  master_id)
      Call Broadcast_Scalar(TCSFSER, master_id)
      Call Broadcast_Scalar(TASFSER, master_id)
      Call Broadcast_Scalar(TSRSF,   master_id)
      Call Broadcast_Scalar(TSSSF,   master_id)
      Call Broadcast_Scalar(ISSFLDN, master_id)
      Call Broadcast_Scalar(ISSFLFE, master_id)
      Call Broadcast_Scalar(SFLKILL, master_id)
      
      Call Broadcast_Array(TSFSER, master_id)
      Call Broadcast_Array(RKDSFL, master_id)
      Call Broadcast_Array(WSFLST, master_id)
      Call Broadcast_Array(WSFLSM, master_id)
      Call Broadcast_Array(DSFLMN, master_id)
      Call Broadcast_Array(DSFLMX, master_id)
      Call Broadcast_Array(SFNTBE, master_id)
      Call Broadcast_Array(SFATBT, master_id)
      
    ENDIF

    ! *** READ VEGETATION DATA FROM VEGE.INP AND VEGSER.INP

  IF( ISVEG == 1 )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING VEGE.INP'
      OPEN(1,FILE='vege.inp',STATUS='UNKNOWN',SHARED)
      STR=READSTR(1)
      READ(1,*) MVEGTYP, MVEGOW, NVEGSER, UVEGSCL
      IF( UVEGSCL <= 0. ) UVEGSCL = 1.E-12
      
      DO M=1,MVEGTYP
        READ(1,*,ERR=3120)IDUM, NVEGSERV(M), RDLPSQ(M), BPVEG(M), HPVEG(M), ALPVEG(M), SCVEG(M)
        BDLTMP   = BPVEG(M)*BPVEG(M)*RDLPSQ(M)          ! *** "ad" - dimensionless population density
        !PVEGX(M) = 1. - BETVEG(M)*BDLTMP               ! *** Not used
        !PVEGY(M) = 1. - GAMVEG(M)*BDLTMP               ! *** Not used
        PVEGZ(M) = MAX(1.E-18,(1. - ALPVEG(M)*BDLTMP))  ! *** Turbulence/canopy adjustment factor to increase drag and turbulence. 0-No adjustment, 1-Infinite drag (invalid) (dimensionless)
        BDLPSQ(M) = BPVEG(M)*RDLPSQ(M)                  ! *** "a" - projected plant area per unit volume  (1/m)
      ENDDO
      CLOSE(1)
    end if !***end on master

    Call Broadcast_Scalar(MVEGTYP, master_id)
    Call Broadcast_Scalar(MVEGOW , master_id)
    Call Broadcast_Scalar(NVEGSER, master_id)
    Call Broadcast_Scalar(UVEGSCL, master_id)

    !Call Broadcast_Array(PVEGX,    master_id)
    !Call Broadcast_Array(PVEGY,    master_id)
    Call Broadcast_Array(PVEGZ,    master_id)
    Call Broadcast_Array(BDLPSQ ,  master_id)
    Call Broadcast_Array(NVEGSERV, master_id)
    Call Broadcast_Array(RDLPSQ  , master_id)
    Call Broadcast_Array(BPVEG   , master_id)
    Call Broadcast_Array(HPVEG   , master_id)
    Call Broadcast_Array(ALPVEG  , master_id)
    !Call Broadcast_Array(BETVEG  , master_id)
    !Call Broadcast_Array(GAMVEG  , master_id)
    Call Broadcast_Array(SCVEG   , master_id)

    IF( NVEGSER > 0 )THEN
      DO M=1,NVEGSER
        MVEGTLAST(M)=1
      ENDDO
    ENDIF

    
      ! *** Read MHK file
      IF( LMHK )THEN !MHK devices
        if( process_id == master_id )THEN
            WRITE(*,'(A)')'READING MHK.INP'
            OPEN(1,FILE='mhk.inp',STATUS='UNKNOWN')
            
            STR=READSTR(1)
            
            READ(1,*,ERR=3122)MHKTYP,NFLAGPWR,UPSTREAM,OUTPUTFLAG
        end if ! *** end on master
        
        ! *** send to all
        Call Broadcast_Scalar(MHKTYP,    master_id)
        Call Broadcast_Scalar(NFLAGPWR,  master_id)
        Call Broadcast_Scalar(UPSTREAM,  master_id)
        Call Broadcast_Scalar(OUTPUTFLAG,master_id)
        
        
        ALLOCATE(BOFFMHK(MHKTYP),BOFFSUP(MHKTYP))
        ALLOCATE(TOFFMHK(MHKTYP),TOFFSUP(MHKTYP))
        
        if( process_id == master_id )THEN

            IF( NFLAGPWR == 1 )THEN
              DO M=1,MHKTYP
                READ(1,*,ERR=3122)WIDTHMHK(M),WIDTHSUP(M), BOFFMHK(M),BOFFSUP(M),TOFFMHK(M),TOFFSUP(M),CTMHK(M),CDSUP(M),VMINCUT(M),VMAXCUT(M),DENMHK(M)
                CTMHK(M)=CTMHK(M)*DENMHK(M)
                CDSUP(M)=CDSUP(M)*DENMHK(M)

                DO L=2,LA
                  IF( M+90 == MVEGL(L) )THEN
                    ZMINMHK(M,L)=BELV(L)+BOFFMHK(M)
                    ZMAXMHK(M,L)=BELV(L)+TOFFMHK(M)
                    ZMINSUP(M,L)=BELV(L)+BOFFSUP(M)
                    ZMAXSUP(M,L)=BELV(L)+TOFFSUP(M)
                    DIAMMHK=ZMAXMHK(M,L)-ZMINMHK(M,L)
                    IF( DIAMMHK<0.0 )THEN   !error check
                      PRINT *,'MHK ZMIN > ZMAX'
                      CALL STOPP('.')
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
              READ(1,*) !skip the header line
              READ(1,*)BETAMHK_D,BETAMHK_P,CE4MHK,PB_COEF
          end if
          ! *** send to all processes
          Call Broadcast_Array(WIDTHMHK, master_id)
          Call Broadcast_Array(WIDTHSUP, master_id)
          Call Broadcast_Array(BOFFMHK,  master_id)
          Call Broadcast_Array(BOFFSUP,  master_id)
          Call Broadcast_Array(TOFFMHK,  master_id)
          Call Broadcast_Array(TOFFSUP,  master_id)
          Call Broadcast_Array(CTMHK,    master_id)
          Call Broadcast_Array(CDSUP,    master_id)
          Call Broadcast_Array(VMINCUT,  master_id)
          Call Broadcast_Array(VMAXCUT,  master_id)
          Call Broadcast_Array(DENMHK,   master_id)

          Call Broadcast_Array(ZMINMHK, master_id)
          Call Broadcast_Array(ZMAXMHK, master_id)
          Call Broadcast_Array(ZMINSUP, master_id)
          Call Broadcast_Array(ZMAXSUP, master_id)
          
          Call Broadcast_Scalar(BETAMHK_D,master_id)
          Call Broadcast_Scalar(BETAMHK_P,master_id)
          Call Broadcast_Scalar(CE4MHK,   master_id)
          Call Broadcast_Scalar(PB_COEF,  master_id)
          
        ELSEIF( NFLAGPWR == 2 )THEN
          PRINT*,'Not available yet'

        ELSEIF( NFLAGPWR == 3 )THEN !FFP input style
          if(process_id == master_id )then
              READ(1,*,ERR=3122)WIDTHMHK(M),WIDTHSUP(M),BOFFMHK(M),HEIGHTMHK(M),HEIGHTSUP(M),REFELEV(M),CTMHK(M),CDSUP(M),VMINCUT(M),VMAXCUT(M),DENMHK(M)
              CTMHK(M)=CTMHK(M)*DENMHK(M)
              CDSUP(M)=CDSUP(M)*DENMHK(M)
              DO L=2,LA
                IF( M+90 == MVEGL(L) )THEN
                  ZMINMHK(M,L)=BELV(L)+REFELEV(M)
                  ZMAXMHK(M,L)=BELV(L)+REFELEV(M)+HEIGHTMHK(M)
                  ZMINSUP(M,L)=BELV(L)
                  ZMAXSUP(M,L)=BELV(L)+REFELEV(M)+HEIGHTSUP(M)
                  DIAMMHK=ZMAXMHK(M,L)-ZMINMHK(M,L)
                  IF( DIAMMHK<0.0 )THEN
                    PRINT*,'MHK ZMIN > ZMAX'
                    CALL STOPP('')
                  ENDIF
                ENDIF
              ENDDO

              READ(1,*) !skip the header line
              READ(1,*)BETAMHK_D,BETAMHK_P,CE4MHK,PB_COEF
          end if
          
          Call Broadcast_Array(WIDTHMHK, master_id)
          Call Broadcast_Array(WIDTHSUP, master_id)
          Call Broadcast_Array(BOFFMHK,  master_id)
          Call Broadcast_Array(HEIGHTMHK,master_id)
          Call Broadcast_Array(HEIGHTSUP,master_id)
          Call Broadcast_Array(REFELEV,  master_id)
          Call Broadcast_Array(CTMHK,    master_id)
          Call Broadcast_Array(CDSUP,    master_id)
          Call Broadcast_Array(VMINCUT,  master_id)
          Call Broadcast_Array(VMAXCUT,  master_id)
          Call Broadcast_Array(DENMHK,   master_id)
          
          Call Broadcast_Array(ZMINMHK, master_id)
          Call Broadcast_Array(ZMAXMHK, master_id)
          Call Broadcast_Array(ZMINSUP, master_id)
          Call Broadcast_Array(ZMAXSUP, master_id)
          
          Call Broadcast_Scalar(BETAMHK_D,master_id)
          Call Broadcast_Scalar(BETAMHK_P,master_id)
          Call Broadcast_Scalar(CE4MHK,   master_id)
          Call Broadcast_Scalar(PB_COEF,  master_id)
          
        ENDIF
        
        if(process_id == master_id )then
            CLOSE(1)
        end if
        
        DO L=2,LA
          IF( MVEGL(L)>90 )THEN
            IF( (DIAMMHK>DXP(L) .OR. DIAMMHK>DYP(L)) .AND. DENMHK(MVEGL(L)-90)>1.0 )THEN
              PRINT*,'MHK DIAMETER EXCEEDS CELL SIZE'
              PRINT*,'AND DENSITY >= 1'
              CALL STOPP('')
            ENDIF
          ENDIF
        ENDDO
        
      ENDIF
      !!!END SCJ BLOCK

    ENDIF

    GOTO 3124
3120 CALL STOPP('READ ERROR FOR FILE VEGE.INP')
3122 CALL STOPP('READ ERROR FOR FILE MHK.INP')

3124 CONTINUE

    IF( NVEGSER >= 1 )THEN
      if(process_id == master_id )then
          WRITE(*,'(A)')'READING VEGSER.INP'
          OPEN(1,FILE='vegser.inp',STATUS='UNKNOWN')

          ! *** SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          DO NS=1,NVEGSER
            READ(1,*,IOSTAT=ISO) MVEGSER(NS),TCVEGSER(NS),TAVEGSER(NS)
            IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE VEGSER.INP')
            DO M=1,MVEGSER(NS)
              READ(1,*,IOSTAT=ISO)TVEGSER(M,NS),VEGSERR(M,NS),VEGSERB(M,NS),VEGSERH(M,NS)
              IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE VEGSER.INP')
              TVEGSER(M,NS)=TVEGSER(M,NS)+TAVEGSER(NS)
            ENDDO
          ENDDO
          CLOSE(1)
      end if
      
      Call Broadcast_Array(MVEGSER , master_id)
      Call Broadcast_Array(TCVEGSER, master_id)
      Call Broadcast_Array(TAVEGSER, master_id)
      Call Broadcast_Array(TVEGSER,  master_id)
      Call Broadcast_Array(VEGSERR,  master_id)
      Call Broadcast_Array(VEGSERB,  master_id)
      Call Broadcast_Array(VEGSERH,  master_id)
      
      ! *** REINITIALIZE CLASSES HAVING TIME SERIES INFORMATION
      DO M=1,MVEGTYP
        IF( NVEGSERV(M) > 0 )THEN
          NS=NVEGSERV(M)
          RDLPSQ(M)=VEGSERR(1,NS)
          BPVEG(M)=VEGSERB(1,NS)
          HPVEG(M)=VEGSERH(1,NS)
          BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
          !PVEGX(M)=1.-BETVEG(M)*BDLTMP   ! *** Not used
          !PVEGY(M)=1.-BETVEG(M)*BDLTMP   ! *** Not used
          PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
          BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
        ENDIF
      ENDDO
    ENDIF
    GOTO 7122
7120 CALL STOPP('READ ERROR FOR FILE VEGSER.INP')
7122 CONTINUE

    ! *** *******************************************************************C
    ! *** READ BANK EROSION MAP AND TIME SERIES FILE
    IF( (ISTRAN(6)>0 .OR. ISTRAN(7) > 0 ) .AND. ISBKERO >= 1 )THEN
      if(process_id == master_id )then
          WRITE(*,'(A)')'READING BEMAP.INP'
          OPEN(1,FILE='bemap.inp',STATUS='UNKNOWN')
          
          ! *** SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          READ(1,*)NBEPAIR,NBESER
          
          DO NS=1,NBEPAIR
            READ(1,*)IBANKBE(NS),JBANKBE(NS),ICHANBE(NS),JCHANBE(NS),NBESERN(NS),FBESER(NS)
          ENDDO
          
          CLOSE(1)
      end if
      ! *** send to all
      Call Broadcast_Scalar(NBEPAIR,master_id)
      Call Broadcast_Scalar(NBESER, master_id)
      
      Call Broadcast_Array(IBANKBE, master_id)
      Call Broadcast_Array(JBANKBE, master_id)
      Call Broadcast_Array(ICHANBE, master_id)
      Call Broadcast_Array(JCHANBE, master_id)
      Call Broadcast_Array(NBESERN, master_id)
      Call Broadcast_Array(FBESER , master_id)
      
    ENDIF


    IF( (ISTRAN(6)>0 .OR. ISTRAN(7) > 0 ) .AND. ISBKERO >= 1 .AND. NBESER > 0 )THEN
      if(process_id == master_id )then
          WRITE(*,'(A)')'READING BESER.INP'
          OPEN(1,FILE='beser.inp',STATUS='UNKNOWN')

          ! *** SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)

          DO NS=1,NBESER

            READ(1,*)MBESER(NS),TCBESER(NS),TABESER(NS),RMULADJ,ADDADJ
            MBETLAST(NS)=1

            DO M=1,MBESER(NS)
              READ(1,*)TBESER(M,NS),BESER(M,NS),FWCBESER(M,NS)
              TBESER(M,NS)=TBESER(M,NS)+TABESER(NS)
              BESER(M,NS)=RMULADJ*(BESER(M,NS)+ADDADJ)
            ENDDO

          ENDDO

          CLOSE(1)
      end if
      ! *** Send to all
      Call Broadcast_Array(MBESER, master_id)
      Call Broadcast_Array(TCBESER, master_id)
      Call Broadcast_Array(TABESER, master_id)
      Call Broadcast_Array(TBESER, master_id)
      Call Broadcast_Array(BESER, master_id)
      Call Broadcast_Array(FWCBESER, master_id)
      
    ENDIF

    ! *** *******************************************************************C
    ! *** READ ZONALLY VARYING SEDIMENT BED PARTICLE MIXING
    !----------------------------------------------------------------------C
    ITMPPMX=0
    DO NT=1,NTOX
      IF( ISPMXZ(NT) == 1 )ITMPPMX=1
    ENDDO

    IF( ISTRAN(5) > 0 .AND. ITMPPMX == 1 )THEN
      if(process_id == master_id )then
           WRITE(*,'(A)')'READING PARTMIX.INP'
           OPEN(1,FILE='partmix.inp')
           
           ! *** SKIP OVER TITLE AND AND HEADER LINES
           STR=READSTR(1)
           
           !#######################################################################
           !     HQI change to input multiplication scale factor for particle mixing rate
           !     RM 10/06/05
           READ(1,*)NPMXZ,NPMXPTS,PMIXSF
           
           !#######################################################################
           DO NZ=1,NPMXZ
             DO NP=1,NPMXPTS
               READ(1,*)PMXDEPTH(NP,NZ),PMXCOEF(NP,NZ)
               !#######################################################################
               !     HQI change to input multiplication scale factor for particle mixing rate
               !     RM 10/06/05
               PMXCOEF(NP,NZ) = PMIXSF*PMXCOEF(NP,NZ)
               !#######################################################################
             ENDDO
           ENDDO
           
           CLOSE(1)
        end if
        
        ! *** send to all processes
        Call Broadcast_Scalar(NPMXZ,  master_id)
        Call Broadcast_Scalar(NPMXPTS,master_id)
        Call Broadcast_Scalar(PMIXSF, master_id)
        
        Call Broadcast_Array(PMXDEPTH, master_id)
        Call Broadcast_Array(PMXCOEF , master_id)
        
      if(process_id == master_id )then
          WRITE(*,'(A)')'READING PMXMAP.INP'
          OPEN(1,FILE='pmxmap.inp')
          
          ! *** SKIP OVER TITLE AND AND HEADER LINES
          STR=READSTR(1)
          
          DO L=2,LC_Global-1
            READ(1,*)LDUM,IDUM,JDUM,LPMXZ_Global(L)
          ENDDO
          
          CLOSE(1)
       end if
       ! *** send to all processes
       Call Broadcast_Array(LPMXZ_Global, master_id)
       
       DO LG = 2 , LC_GLOBAL - 1
          L = MAP2LOCAL(LG).LL
          IF( L > 1 )THEN
            LPMXZ(L) = LPMXZ_GLOBAL(LG)
          END IF
       END DO
       
    ENDIF

  call ReadFields()     ! 2018-10-12, NTL: Read time & space varying fields
  call ReadCyclones()   ! 2021-03-31, NTL: Read cyclone tracks
  
END

! ***************************************************************************************
! *** DSI UTILITIES
FUNCTION PARSE_REAL(INLINE)

  USE GLOBAL,ONLY: IK4
  IMPLICIT NONE
  INTEGER(IK4) :: I1,I2,ILEN,IPOS,ILEN2
  REAL :: PARSE_REAL
  CHARACTER*(*) INLINE
  CHARACTER*15  CVAL

  ILEN=LEN_TRIM(INLINE)
  PARSE_REAL=0.
  DO I1=1,ILEN
    IF( INLINE(I1:I1) == ':' )THEN
      DO IPOS=I1+1,ILEN
        IF( INLINE(IPOS:IPOS)/=' ')EXIT
      ENDDO
      IF( IPOS>ILEN)RETURN
      CVAL=INLINE(IPOS:ILEN)
      ILEN2=LEN_TRIM(CVAL)
      DO I2=1,ILEN2
        IF( CVAL(I2:I2) == ' ' .OR. CVAL(I2:I2) == ',' .OR. I2 == ILEN2 )THEN
          READ(CVAL(1:I2),'(F12.1)',ERR=999)PARSE_REAL
          RETURN
        ENDIF
      ENDDO
    ENDIF
  ENDDO
999 CALL STOPP(' ERROR PARSING REAL')
    
END FUNCTION

FUNCTION PARSE_LOGICAL(INLINE)

  USE GLOBAL,ONLY: IKV
  IMPLICIT NONE

  INTEGER :: ILEN,IC,JC,IPOS,ILEN2
  CHARACTER*(*) INLINE
  CHARACTER*12  CVAL
  LOGICAL       PARSE_LOGICAL

  ILEN=LEN_TRIM(INLINE)
  DO IC=1,ILEN
    IF( INLINE(IC:IC) == ':' )THEN
      DO IPOS=IC+1,ILEN
        IF( INLINE(IPOS:IPOS)/=' ')EXIT
      ENDDO
      IF( IPOS>ILEN)RETURN
      CVAL=INLINE(IPOS:ILEN)
      ILEN2=LEN_TRIM(CVAL)
      DO JC=1,ILEN2
        IF( CVAL(JC:JC) == ' ' .OR. CVAL(JC:JC) == ',' .OR. JC == ILEN2 )THEN
          IF( CVAL(1:1) == 'T' .OR. CVAL(1:1) == 'Y' )THEN
            PARSE_LOGICAL=.TRUE.
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  PARSE_LOGICAL=.FALSE.

END FUNCTION

SUBROUTINE PARSESTRING(string,substring) 
    character(*),intent(inout) :: string
    character(len=25) :: substring
    integer :: posBT(2), pos

    posBT(1) = index(string,' ')      ! Blank
    posBT(2) = index(string,char(9))  ! Tab
    pos = minval( posBT, posBT > 0 )
    substring = string(1:pos)
    string = adjustl(string(pos+1:))
end SUBROUTINE
    
SUBROUTINE READ_SUBSET
    USE GLOBAL
    Use INFOMOD ,ONLY:SKIPCOM, READSTR
    USE XYIJCONV,ONLY:XY2IJ
    Use Variables_MPI
    Use Broadcast_Routines
    
    IMPLICIT NONE
    character(256) :: string
    character(25) :: substring
    integer :: ISO,NS,IS,NP
     
    If( process_id == master_id )THEN
      OPEN(1,FILE='subset.inp',action='read')
      CALL SKIPCOM(1,'*')
      READ(1,*,IOSTAT=ISO) NSUBSET
      IF( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
    End if
    
    Call Broadcast_Scalar(NSUBSET, master_id)
    
    ALLOCATE(HFREGRP(NSUBSET))
    ALLOCATE(IJHFRE(NSUBSET),HFREDAYBG(NSUBSET),HFREDAYEN(NSUBSET))
    ALLOCATE(HFREDUR(NSUBSET),HFREDAY(NSUBSET),HFREMIN(NSUBSET),NPNT(NSUBSET))
    
    DO NS=1,NSUBSET
      If( process_id == master_id )Then
        CALL SKIPCOM(1,'*')
        READ(1,*,IOSTAT=ISO) IS,IJHFRE(IS),HFREDAYBG(IS),HFREDUR(IS),HFREMIN(IS),NPNT(IS)
        IF( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
      End if
          
      Call Broadcast_Scalar(IS,            master_id)
      Call Broadcast_Scalar(NPNT(IS),      master_id)
      
      ALLOCATE (HFREGRP(IS)%ICEL(NPNT(IS)),HFREGRP(IS)%JCEL(NPNT(IS)))
      ALLOCATE (HFREGRP(IS)%XCEL(NPNT(IS)),HFREGRP(IS)%YCEL(NPNT(IS)))
      ALLOCATE (HFREGRP(IS)%NAME(NPNT(IS)))
      
      If( process_id == master_id )Then
        DO NP=1,NPNT(IS)
          CALL SKIPCOM(1,'*')
          read(1,'(a)') string
          string = adjustl(trim(string))
          IF( IJHFRE(IS) == 1 )THEN
            !READ(1,*,advance='no',IOSTAT=ISO) HFREGRP(IS)%ICEL(NP),HFREGRP(IS)%JCEL(NP)
            call PARSESTRING(string, substring)
            read(substring,*) HFREGRP(IS)%ICEL(NP)
            call PARSESTRING(string, substring)
            read(substring,*) HFREGRP(IS)%JCEL(NP)
          ELSE
            !READ(1,*,advance='no',IOSTAT=ISO) HFREGRP(IS)%XCEL(NP),HFREGRP(IS)%YCEL(NP)
            call PARSESTRING(string, substring)
            read(substring,*) HFREGRP(IS)%XCEL(NP)
            call PARSESTRING(string, substring)
            read(substring,*) HFREGRP(IS)%YCEL(NP)
          ENDIF
          IF( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
          string = TRIM(string)
          IF (string(1:1)=='!') string = string(2:20)
          HFREGRP(IS)%NAME(NP) = TRIM(string)
        ENDDO
      End if
      
    ENDDO
 
    Call Broadcast_Array(IJHFRE,    master_id)
    Call Broadcast_Array(HFREDAYBG, master_id)
    Call Broadcast_Array(HFREDUR,   master_id)
    Call Broadcast_Array(HFREMIN,   master_id)

      
    DO NS=1,NSUBSET
      Call Broadcast_Array(HFREGRP(NS)%ICEL, master_id)
      Call Broadcast_Array(HFREGRP(NS)%JCEL, master_id)
      Call Broadcast_Array(HFREGRP(NS)%XCEL, master_id)
      Call Broadcast_Array(HFREGRP(NS)%YCEL, master_id)
    ENDDO
    
    If( process_id == master_id ) CLOSE(1)
    
    
END SUBROUTINE    
