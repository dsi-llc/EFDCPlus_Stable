! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  MODULE EFDCOUT

#ifndef GNU  
  USE IFPORT
#endif
  use GLOBAL
  use DRIFTER, only:DRIFTER_OUT
  use HYDSTRUCMOD
  use SHELLFISHMOD
  use FIELDS
  use Variables_WQ
  use WQ_RPEM_MODULE
  
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out

  implicit none

  ! *** DO NOT CHANGE THESE PRECISIONS AS THEY ARE LINKED TO EFDC_EXPLORER
  real(RKD),PRIVATE,save :: EETIME, SUMTIME

  integer(IK4),PRIVATE,save :: NSXD,NBCCELLS,NCELLLIST,INITBCOUT = 0
  integer(IK4),PRIVATE,save :: IWQ(40)
  integer(IK4),PRIVATE,save :: NSEDSTEPS, NBEDSTEPS, NWQVAR, CELL3D
  integer(IK4),PRIVATE,save,allocatable,dimension(:) :: BCCELLS

  real,PRIVATE,save,allocatable,dimension(:,:) :: BCQSUM  ! @todo define!!!
  real,PRIVATE,save,allocatable,dimension(:,:) :: BCUHDY2 ! @todo define!!!
  real,PRIVATE,save,allocatable,dimension(:,:) :: BCVHDX2 ! @todo define!!!
  real,PRIVATE,save,allocatable,dimension(:)   :: BCQSUME ! @todo define!!!

  character(30) :: FILENAME

  Character(1) :: EE_PARALLEL_ID !< This is used to write out an integer based on the processor that is writing out

  integer :: EE_UNIT = 95
  integer :: IERRio, NERR, IORIGIN
  contains

  SUBROUTINE EE_LINKAGE(JSEXPLORER)
  !------------------------------------------------------------------------
  ! ***  SUBROUTINE EE_LINKAGE (OLD EEXPOUT.FOR) WRITES BINARY OUTPUT FILES:
  ! ***    EE_HYD    - WATER DEPTH AND VELOCITY
  ! ***    EE_WC     - WATER COLUMN AND TOP LAYER OF SEDIMENTS
  ! ***    EE_BC     - EFDC COMPUTED BOUNDARY FLOWS
  ! ***    EE_BED    - SEDIMENT BED LAYER INFORMATION
  ! ***    EE_WQ     - WATER QUALITY INFORMATION FOR THE WATER COLUMN
  ! ***    EE_SD     - SEDIMENT DIAGENSIS INFORMATION
  ! ***    EE_RPEM   - ROOTED PLANT AND EPIPHYTE MODEL
  ! ***    EE_ARRAYS - GENERAL/USER DEFINED ARRAY DUMP. LINKED TO
  ! ***                EFDC_EXPLORER FOR DISPLAY
  ! ***    EE_SEDZLJ - SEDIMENT BED DATA FOR SEDZLJ SUB-MODEL
  !------------------------------------------------------------------------
  integer(IK4),intent(IN) :: JSEXPLORER
  integer :: I, IW, J, NS, NX, NT
  
  if( ISDYNSTP == 0 )then
    DELT = DT
  else
    DELT = DTDYN
  endif

  if( JSEXPLORER == 1 )then
    call HEADEROUT
  elseif( JSEXPLORER == -1 )then
    ! *** FORCE ALL OUTPUT
    NSEDSTEPS = 32000
    NBEDSTEPS = 32000
  endif

  ! *** SET TIME EE LINKAGE FILES
  EETIME = TIMESEC
  EETIME = EETIME/86400._8

  call WSOUT
  call VELOUT
  
  ! *** disable writing ou to EE_BC.OUT with MPI
  if(num_Processors > 1 )then
      
  else
    call BCOUT
  endif

  if( ISSPH(8) >= 1 ) CALL WCOUT

  if( LSEDZLJ )then
    call SEDZLJOUT
    
    !*** Initilize deposition and erosion fluxes for each output frequency
    DEP_SED_FLX = 0.0
    ERO_SED_FLX = 0.0
  else
    if( ISBEXP >= 1 .and. KB > 1 )then
      if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 ) CALL BEDOUT
    endif
  endif

  if( ISTRAN(8) > 0 )then
    call WQOUT
    if( IWQBEN > 0 .and. ISSDBIN /= 0 ) CALL SDOUT(JSEXPLORER)
    if( ISRPEM > 0) CALL RPEMOUT(JSEXPLORER)
    if( ISFFARM > 0) CALL SHELLFISHOUT()
  endif

  if( ISINWV == 2 .and. JSEXPLORER /= 1 ) CALL ARRAYSOUT
  if( ISPD > 0 .and. JSEXPLORER == -1 ) CALL DRIFTER_OUT(.TRUE.)

  END SUBROUTINE

  SUBROUTINE HEADEROUT
  integer(IK4) :: NACTIVE,VER,HSIZE,BSIZE,I,NS,MW,NW,L,K,ITYPE,ITIMEVAR,LL,N2D,N3D
  character(8) :: ARRAYNAME

  if( iswqlvl == 0 .and. NFIXED > 0 )then
    NWQVAR = NWQVAR + 1
  else
    NWQVAR = NWQV              ! *** Number of water quality variables no longer changes 
  endif
  NACTIVE = LA_Global - 1

  ! *** NUMBER OF 3D CELLS
  CELL3D = 0
  do L = 2, LA_Global !***Modified for Global
    CELL3D = CELL3D + (KC - KSZ_Global(L) + 1)
  enddo

  ! *** WATER DEPTHS
  FILENAME = OUTDIR//'EE_WS.OUT'
  VER   = 8400
  HSIZE = 6*4
  BSIZE = (2 + 3 + NACTIVE)*4
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
  close(EE_UNIT,STATUS = 'DELETE')
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
  write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
  write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(NACTIVE,4)
  close(EE_UNIT,STATUS = 'KEEP')

  ! *** VELOCITY
  FILENAME = OUTDIR//'EE_VEL.OUT'
  VER   = 8400
  HSIZE = (9 + 4*LCM_Global)*4 !***Modified for MPI
  BSIZE = (2 + 2 + 3*CELL3D)*4
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
  close(EE_UNIT,STATUS = 'DELETE')
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
  write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
  write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(NACTIVE,4)
  write(EE_UNIT) REAL(RSSBCE_Global,4)
  write(EE_UNIT) REAL(RSSBCW_Global,4)
  write(EE_UNIT) REAL(RSSBCS_Global,4)
  write(EE_UNIT) REAL(RSSBCN_Global,4)
  close(EE_UNIT,STATUS = 'KEEP')

  if( ISSPH(8) >= 1 )then
    ! *** WATER COLUMN AND TOP LAYER OF SEDIMENT
    FILENAME = OUTDIR//'EE_WC.OUT'
    NSXD = NSED+NSND
    VER = 11700
    HSIZE = (29 + NSXD + 1)*4      ! *** +1 is for NSED2
    BSIZE = BLOCKWC(CELL3D)

    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
    close(EE_UNIT,STATUS = 'DELETE')
    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
    write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
    write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
    write(EE_UNIT) (INT(ISTRAN(I),4),I = 1,7)
    write(EE_UNIT) INT(NDYE,4),INT(NSED,4),INT(NSND,4),INT(NTOX,4),INT(NSED2,4)
    write(EE_UNIT) INT(ISWAVE,4),INT(ISBEDSTR,4),LOGICAL(LSEDZLJ,4),INT(ICALC_BL,4),REAL(TEMBO,4)
    write(EE_UNIT) INT(IEVAP,4),INT(ISGWIE,4),INT(ISICE,4)
    
    do NS = 1,NSXD
      write(EE_UNIT) REAL(SEDDIA(NS),4)
    enddo
    close(EE_UNIT,STATUS = 'KEEP')
  endif

  ! *** SEDFLUME MODEL RESULTS
  if( LSEDZLJ )then
    ! *** SEDIMENT BED
    FILENAME = OUTDIR//'EE_SEDZLJ.OUT'
    VER = 12000
    HSIZE = 23*4
    BSIZE = BLOCKSEDZLJ(CELL3D)

    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
    close(EE_UNIT,STATUS = 'DELETE')
    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
    write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
    write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
    write(EE_UNIT) (INT(ISTRAN(I),4),I = 1,7)
    write(EE_UNIT) INT(NSEDS,4),INT(ITBM,4),INT(NSICM,4),INT(NTOX,4),INT(ICALC_BL,4),INT(NSEDS2,4)
    close(EE_UNIT,STATUS = 'KEEP')
    NBEDSTEPS = ISBEXP - 1            ! *** The very first bed snapshot will be saved
  else
    ! *** SEDIMENT BED LAYERS FOR ORIGINAL SEDIMENT TRANSPORT APPROACH
    if( ISBEXP >= 1 .and. KB > 1 )then
      if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
        FILENAME = OUTDIR//'EE_BED.OUT'
        VER = 8400
        HSIZE = (18 + NSXD)*4
        BSIZE = BLOCKBED(CELL3D)

        open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
        close(EE_UNIT,STATUS = 'DELETE')
        open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
        write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
        write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
        write(EE_UNIT) (INT(ISTRAN(I),4),I = 1,7)
        write(EE_UNIT) INT(NSED,4),INT(NSND,4),INT(NTOX,4)

        do NS = 1,NSXD
          write(EE_UNIT) REAL(SEDDIA(NS),4)
        enddo
        close(EE_UNIT,STATUS = 'KEEP')
        NBEDSTEPS = ISBEXP - 1            ! *** The very bed first snapshot will be saved, regardless of the skip count
      endif
    endif
  endif

  ! *** WATER QUALITY MODEL (HEM3D) RESULTS
  if( ISTRAN(8) > 0 )then
    FILENAME = OUTDIR//'EE_WQ.OUT'

    IWQ = 0
    BSIZE = 0

    do NW = 1,NWQVAR
      IWQ(NW) = ISKINETICS(NW)
      if( IWQ(NW) > 0) BSIZE = BSIZE + 1
    enddo
    
    VER = 10300                 !< version 10.3
    HSIZE = (13 + NWQVAR)*4
    BSIZE = (2 + BSIZE*CELL3D)*4
    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
    close(EE_UNIT,STATUS = 'DELETE')
    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
    write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
    write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
    write(EE_UNIT) INT(NALGAE,4),INT(NZOOPL,4),INT(NWQVAR,4)
    write(EE_UNIT) (INT(IWQ(NW),4),NW = 1,NWQVAR)
    close(EE_UNIT,STATUS = 'KEEP')

    ! *** save SEDIMENT DIAGENESIS RESULTS
    if( ISSDBIN /= 0 )then
      FILENAME = OUTDIR//'EE_SD.OUT'
      VER = 8400
      HSIZE = 8*4
      BSIZE = (2 + (6*3 + 21)*NACTIVE)*4
      open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
      close(EE_UNIT,STATUS = 'DELETE')
      open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
      write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
      close(EE_UNIT,STATUS = 'KEEP')
      NSEDSTEPS = -1
    endif

    if( ISRPEM > 0 )then
      FILENAME = OUTDIR//'EE_RPEM.OUT'
      VER = 8400
      HSIZE = (10 + NACTIVE)*4
      BSIZE = 0
      do L = 2,LA_Global
        if( LMASKRPEM_Global(L) .or. INITRPEM == 0 )then
          BSIZE = BSIZE + 1
        endif
      enddo
      BSIZE = (2 + 4*BSIZE)*4

      open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
      close(EE_UNIT,STATUS = 'DELETE')
      open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
      write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
      write(EE_UNIT) INT(NRPEM,4),INT(INITRPEM,4)
      write(EE_UNIT) (LOGICAL(LMASKRPEM_Global(L),4),L = 2,LA_Global)
      close(EE_UNIT,STATUS = 'KEEP')
    endif

    if( ISFFARM > 0 )then
      FILENAME = OUTDIR//'EE_SHF.OUT'
      VER = 10000
      HSIZE = (10 + NACTIVE)*4
      N2D = 0
      N3D = 0
      do L = 2,LA_Global
        if( FARMCELL(L) > 0 )then
          N2D = N2D + 1
          N3D = N3D + (KC - KSZ_Global(L) + 1)
        endif
      enddo
      BSIZE = (2 + NSF*(12*N2D + 11*N3D))*4

      open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
      close(EE_UNIT,STATUS = 'DELETE')
      open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      write(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
      write(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(NACTIVE,4)
      write(EE_UNIT) INT(NSF,4),INT(NSFCELLS,4),INT(NSFHARV,4)
      write(EE_UNIT) (INT(FARMCELL(L),4),L = 2,LA_Global)
      close(EE_UNIT,STATUS = 'KEEP')
    endif

  endif

  ! *** BC OUTPUT: QSUM
  !CALL BCOUT_INITIALIZE !***Initializes some arrays that hold some info later on.. not clearly explained

  ! *** Set file name and write out to #output directory
  FILENAME = OUTDIR//'EE_BC.OUT'
  ! *** Setup binary file to write to
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
  close(EE_UNIT,STATUS = 'DELETE')
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)

  ! *** ONLY OUTPUT CELLS THAT HAVE BC'S
  if( KC > 1 )then  ! *** If more than one layer
    NBCCELLS = NBCS + NBCSOP  ! Number of boundary condition cells?
    NBCCELLS = NBCCELLS + LA_Global                                      ! *** EVAP/RAINFALL/ICE
    if( NGWSER > 0 .or. ISGWIT /= 0 ) NBCCELLS = NBCCELLS + LA_Global    ! *** GROUNDWATER

    allocate(BCCELLS(NBCCELLS)) !*** @todo what is this

    ! *** BUILD CELL LIST FOR SELECTIVE QSUM
    NCELLLIST = 0
    ! *** BC CELLS
    do LL = 1, NBCS
      NCELLLIST = NCELLLIST+1
      BCCELLS(NCELLLIST) = LBCS(LL)
    enddo

    ! *** OPEN BC CELLS
    do LL = 1, NBCSOP
      NCELLLIST = NCELLLIST+1
      BCCELLS(NCELLLIST) = LOBCS(LL)
    enddo
  else
    ! *** SINGLE LAYER
    NBCCELLS = LA_Global
    NCELLLIST = 0
    allocate(BCCELLS(1))
  endif
  ! *** Setup values for writing to header

  VER   = 8400 !*** version number, @todo should this be updated??
  HSIZE = (20 + NCELLLIST + NBCS + NPBS + NPBW + NPBE + NPBN)*4 !***Header size, why is there a +20?
  BSIZE = BLOCKBC(CELL3D) !***Block size? what is cell3d supposed to do?

  write(EE_UNIT) INT(VER,4), INT(HSIZE,4), INT(BSIZE,4)
  write(EE_UNIT) INT(IC_Global,4), INT(JC_Global,4), INT(KC,4), INT(NACTIVE,4) !***Modified for global
  write(EE_UNIT) INT(NBCS,4), INT(NBCCELLS,4), INT(NCELLLIST,4)
  write(EE_UNIT) (INT(BCCELLS(L),4),L = 1,NCELLLIST)
  write(EE_UNIT) INT(NPBS,4),  INT(NPBW,4), INT(NPBE,4),    INT(NPBN,4)
  write(EE_UNIT) INT(NQCTL,4), INT(NQWR,4), INT(NQCTLSER,4),INT(NQCTRULES,4)
  write(EE_UNIT) INT(NGWSER,4),INT(ISGWIT,4)
  write(EE_UNIT) (INT(LBCS(L),4),L = 1,NBCS) 
  write(EE_UNIT) (INT(LPBS(L),4),L = 1,NPBS)
  write(EE_UNIT) (INT(LPBW(L),4),L = 1,NPBW)
  write(EE_UNIT) (INT(LPBE(L),4),L = 1,NPBE)
  write(EE_UNIT) (INT(LPBN(L),4),L = 1,NPBN)
  close(EE_UNIT,STATUS = 'KEEP')

  ! *** End BC OUT

  if( ISINWV == 2 )then
    FILENAME = OUTDIR//'EE_ARRAYS.OUT'
    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
    close(EE_UNIT,STATUS = 'DELETE')
    open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)

    ! *** EE LINKAGE VERSION
    VER = 7300
    write(EE_UNIT) INT(VER,4)

    ! *** NUMBER OF TIME VARYING ARRAYS
    NS = 0
    if (IS_ARRAY_OUT(2) == 1 ) NS = NS + 9
    if (IS_ARRAY_OUT(3) == 1 ) NS = NS + 4
    if (IS_ARRAY_OUT(4) == 1 ) then
      NS = NS + 2
      if( ISTOPT(2) == 1 .or. ISTOPT(2) == 2 ) NS = NS + 3
    endif
    if (IS_ARRAY_OUT(5) == 1 ) NS = NS + 1 + 3*NALGAE
    write(EE_UNIT) INT(NS,4)

    ! FLAGS: ARRAY TYPE, TIME VARIABLE
    ! ARRAY TYPE:    0 = L            DIM'D
    !                1 = L,KC         DIM'D
    !                2 = L,0:KC       DIM'D
    !                3 = L,KB         DIM'D
    !                4 = L,KC,NCLASS  DIM'D
    ! TIME VARIABLE: 0 = NOT CHANGING
    !                1 = TIME VARYING

    ! *****************************************************************
    ! *** THE FOLLOWING ARRAYS ARE TIME INVARIANT SO ONLY WRITTEN ONCE
    ! *****************************************************************
    if( IS_ARRAY_OUT(1) == 1 )then
      ITIMEVAR = 0
      ITYPE = 0
      write(EE_UNIT) INT(ITYPE,4),INT(ITIMEVAR,4)
      ARRAYNAME = 'SBX'
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(SBXO(L),4)
      enddo
      
      ITYPE = 0
      write(EE_UNIT) INT(ITYPE,4),INT(ITIMEVAR,4)
      ARRAYNAME = 'SBY'
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(SBYO(L),4)
      enddo
      
      ITYPE = 1
      write(EE_UNIT) INT(ITYPE,4),INT(ITIMEVAR,4)
      ARRAYNAME = 'DZC'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(DZC(L,K),4)
        enddo
      enddo
    endif

    close(EE_UNIT,STATUS = 'KEEP')

  endif
  END SUBROUTINE

  SUBROUTINE WSOUT

  implicit none

  ! *** OUTPUT WATER DEPTH
  integer(IK4) :: VER, HSIZE, BSIZE, I, J, LG, LL, L2
  integer(IK4) :: L, ITMP, LTMP, NN, NS, ISTAT
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME

  FILENAME = OUTDIR//'EE_WS.OUT'
  EE_UNIT = 96                        ! *** EE_WS.OUT

  ! *** If we are restarting and continuing a run
  if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_WS == 0 )then
    RSTFIRST_WS = 1
    write(*,'(A)')'READING TO STARTING TIME FOR WS'
    FSIZE = FILESIZE(FILENAME)

    open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
    read(EE_UNIT) VER,HSIZE,BSIZE
    if( VER /= 8400 )then
      write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
      call STOPP('.')
    endif
    OFFSET = HSIZE + 4

    NS = 0
    do while(OFFSET < FSIZE)
#ifdef GNU
      CALL FSEEK(EE_UNIT,OFFSET,0, ISTAT)

      if(IS_IOSTAT_END(ISTAT)) then
        write(EE_UNIT) INT(N+NRESTART,4), EETIME,REAL(DELT,4), INT(LMINSTEP,4)
        exit
      endif
#else
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      if( EOF(EE_UNIT) )then
        write(EE_UNIT) INT(N+NRESTART,4), EETIME,REAL(DELT,4), INT(LMINSTEP,4)
        exit
      endif
#endif       
      if( ISTAT /= 0 ) exit

      read(EE_UNIT) PTIME

      if( DEBUG )then
        NS = NS+1
        if( NS == 8 )then
          write(*,'(F10.3)')PTIME
          NS = 0
        else
          !write(*,'(F10.3,\)')PTIME
	  WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
        endif
      endif
      if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

      OFFSET = OFFSET + BSIZE
    enddo

    if( DEBUG ) WRITE(*,'(" ")')
    write(*,'(A)')'FINSIHED READING WS'

#ifdef GNU
    call FSEEK(EE_UNIT,8,1, ISTAT)
#else
    ISTAT = FSEEK(EE_UNIT,8,1)
#endif    
  else !***Normal write out
    ! *** Writing out some header information for the current time step
    NERR = 0
    IORIGIN = 1
    100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
    write(EE_UNIT) INT(N+NRESTART,4), EETIME, REAL(DELT,4),INT(LMINSTEP,4)
  endif

  if( ISRESTI == 0 .or. ICONTINUE == 0) NRESTART = 0

  ! *** Check for anti-reflection type BC's and replace the interior cell depths with the boundary values
  do LL = 1,NPBW_GL
    if( ISPBW_GL(LL) == 4 .or. ISPBW_GL(LL) == 5 )then
      I = IG2IL(IPBW_GL(LL))
      J = JG2JL(JPBW_GL(LL))
      if( I > 0 .and. I <= IC_GLOBAL )then
        if( J > 0 .and. J <= JC_GLOBAL )then
          LG = LIJ_GLOBAL(I,J) 
          L2 = LEC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        endif
      endif
    endif
  enddo
  do LL = 1,NPBE_GL
    if( ISPBE_GL(LL) == 4 .or. ISPBE_GL(LL) == 5 )then
      I = IG2IL(IPBE_GL(LL))
      J = JG2JL(JPBE_GL(LL))
      if( I > 0 .and. I <= IC_GLOBAL )then
        if( J > 0 .and. J <= JC_GLOBAL )then
          LG = LIJ_GLOBAL(I,J) 
          L2 = LWC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        endif
      endif
    endif
  enddo
  do LL = 1,NPBS_GL
    if( ISPBS_GL(LL) == 4 .or. ISPBS_GL(LL) == 5 )then
      I = IG2IL(IPBS_GL(LL))
      J = JG2JL(JPBS_GL(LL))
      if( I > 0 .and. I <= IC_GLOBAL )then
        if( J > 0 .and. J <= JC_GLOBAL )then
          LG = LIJ_GLOBAL(I,J) 
          L2 = LNC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        endif
      endif
    endif
  enddo
  do LL = 1,NPBN_GL
    if( ISPBN_GL(LL) == 4 .or. ISPBN_GL(LL) == 5 )then
      I = IG2IL(IPBN_GL(LL))
      J = JG2JL(JPBN_GL(LL))
      if( I > 0 .and. I <= IC_GLOBAL )then
        if( J > 0 .and. J <= JC_GLOBAL )then
          LG = LIJ_GLOBAL(I,J) 
          L2 = LSC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        endif
      endif
    endif
  enddo
 
  write(EE_UNIT) (REAL(HP_Global(L),4), L = 2,LA_Global)

  FLUSH(EE_UNIT)
  close(EE_UNIT,STATUS = 'KEEP')
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_WS.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif

END SUBROUTINE

SUBROUTINE VELOUT

  ! *** WATER DEPTH AND VELOCITY OUTPUT
  integer(IK4) :: VER,HSIZE,BSIZE
  integer(IK4) :: L,K,ITMP,NN,NS,ISTAT
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME

  FILENAME = OUTDIR//'EE_VEL.OUT'
  EE_UNIT = 97                        ! *** EE_VEL.OUT
  if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_VEL == 0 )then
    RSTFIRST_VEL = 1
    write(*,'(A)')'READING TO STARTING TIME FOR VEL'
    FSIZE = FILESIZE(FILENAME)
    open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
    read(EE_UNIT) VER,HSIZE,BSIZE
    if( VER /= 8400 )then
      write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
      call STOPP('.')
    endif
    OFFSET = HSIZE + 4

    NS = 0
    do while(OFFSET < FSIZE)
#ifdef GNU
      call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

      if(IS_IOSTAT_END(ISTAT)) then
        write(EE_UNIT) INT(N+NRESTART,4),EETIME,REAL(DELT,4)
        exit
      endif
#else
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      if( EOF(EE_UNIT) )then
        write(EE_UNIT) INT(N+NRESTART,4),EETIME,REAL(DELT,4)
        exit
      endif
#endif 
      if( ISTAT /= 0 ) exit

      read(EE_UNIT) PTIME

      if( DEBUG )then
        NS = NS+1
        if( NS == 8 )then
          write(*,'(F10.3)')PTIME
          NS = 0
        else
          WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
        endif
      endif
      if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

      OFFSET = OFFSET + BSIZE
    enddo
    if( DEBUG ) WRITE(*,'(" ")')
    write(*,'(A)')'FINISHED READING VEL'

#ifdef GNU
    CALL FSEEK(EE_UNIT,4,1, ISTAT)
#else
    ISTAT = FSEEK(EE_UNIT,4,1)
#endif 
  else
    NERR = 0
    IORIGIN = 2
    100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
    write(EE_UNIT) INT(N+NRESTART,4), EETIME, REAL(DELT,4)
  endif

  if( ISRESTI == 0 .or. ICONTINUE == 0) NRESTART = 0
  
  write(EE_UNIT) ((REAL(U_Global(L,K),4), K = KSZ_Global(L), KC), L = 2, LA_Global)
  write(EE_UNIT) ((REAL(V_Global(L,K),4), K = KSZ_Global(L), KC), L = 2, LA_Global)
  write(EE_UNIT) ((REAL(W_Global(L,K),4), K = KSZ_Global(L), KC), L = 2, LA_Global)

  FLUSH(EE_UNIT)
  close(EE_UNIT,STATUS = 'KEEP')
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_VEL.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif
  
END SUBROUTINE

SUBROUTINE WCOUT

  ! *** WATER COLUMN OF CONSTITUENTS OUTPUT
  integer(IK4) :: VER, I, ITMP, NSED4, NSND4, NTOX4
  integer(IK4) :: NS, HSIZE, BSIZE, ISTAT
  integer(IK4) :: L, K, NT, NX, MD
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP, SHEAR
  real(RKD)    :: PTIME
  logical      :: LTMP

  FILENAME = OUTDIR//'EE_WC.OUT'
  EE_UNIT = 98                        ! *** EE_WC.OUT
  if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_WC == 0 )then
    ! *** Run continuation
    RSTFIRST_WC = 1
    write(*,'(A)')'READING TO STARTING TIME FOR WC'
    FSIZE = FILESIZE(FILENAME)
    open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
    read(EE_UNIT) VER,HSIZE,BSIZE
    if( VER /= 11700 )then
      write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
      call STOPP('.')
    endif
    OFFSET = HSIZE

    NS = 0
    do while(OFFSET < FSIZE)
#ifdef GNU
      call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

      if(IS_IOSTAT_END(ISTAT)) then
        write(EE_UNIT) EETIME
        exit
      endif
#else
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      if( EOF(EE_UNIT) )then
        write(EE_UNIT) EETIME
        exit
      endif
#endif   
      if( ISTAT /= 0 ) exit

      read(EE_UNIT) PTIME

      if( DEBUG )then
        NS = NS+1
        if( NS == 8 )then
          write(*,'(F10.3)')PTIME
          NS = 0
        else
          WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
        endif
      endif
      if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

      OFFSET = OFFSET + BSIZE
    enddo
    if( DEBUG ) WRITE(*,'(" ")')
    write(*,'(A)')'FINISHED READING WC'

  else
    NERR = 0
    IORIGIN = 3
    !PRINT *, 'OPEN: NITER = ', NITER
    100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
    write(EE_UNIT) EETIME
  endif
  
  ! *** WRITE THE TOP LAYER INDEX
  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    write(EE_UNIT) (INT(KBT_Global(L),4), L = 2,LA_Global)
  endif

  ! *** TOTAL BED SHEAR STRESS
  do L = 2,LA_Global
    write(EE_UNIT) REAL(SHEAR_Global(L),4)
  enddo

  if( ISBEDSTR == 1 .and. .not. LSEDZLJ )then
    ! *** Stress from non-cohesive components
    do L = 2,LA_Global
      write(EE_UNIT) REAL(SHEAR_Global2(L),4)
    enddo
  endif

  if( ISWAVE >= 1 )then
    if( LSEDZLJ )then
      QQWV3_Global = QQWV3_Global* 0.1 / (WATERDENS*1000.)   ! *** Convert from dyne/cm2 to Pa and density normalize
    endif
    write(EE_UNIT) (REAL(QQWV3_Global(L),4), L = 2,LA_Global)  ! *** Bed Shear due to Waves Only
    ! *** Shear due to Current Only
    if( LSEDZLJ )then
      do L = 2,LA_Global
        SHEAR = SHEAR_Global(L) - QQWV3_Global(L)
        !SHEAR = max(SHEAR,0.0)
        write(EE_UNIT) SHEAR             ! *** Bed Shear due to Current Only
      enddo
    else
      do L = 2,LA_Global
        SHEAR = ( RSSBCE_Global(L)*TBX_Global(LEC_Global(L)) + RSSBCW_Global(L)*TBX_Global(L) )**2  + &
                ( RSSBCN_Global(L)*TBY_Global(LNC_Global(L)) + RSSBCS_Global(L)*TBY_Global(L) )**2
        SHEAR = 0.5*SQRT(SHEAR)
        write(EE_UNIT) SHEAR             ! *** Bed Shear due to Current Only
      enddo
    endif
    if( ISWAVE >= 3 )then
      write(EE_UNIT) (REAL(WV_HEIGHT_Global(L),4), L = 2,LA_Global)
      write(EE_UNIT) (REAL(WV_PERIOD_Global(L),4), L = 2,LA_Global)
      write(EE_UNIT) (REAL(WV_DIR_Global(L),4), L = 2,LA_Global)
      if( ISWAVE == 4 )then
        write(EE_UNIT) (REAL(WV_DISSIPA_Global(L),4), L = 2,LA_Global)      ! *** DISSIPATION
        write(EE_UNIT) (REAL(WVHUU_Global(L,KC),4), L = 2,LA_Global)        ! *** SXX (M3/S2)
        write(EE_UNIT) (REAL(WVHVV_Global(L,KC),4), L = 2,LA_Global)        ! *** SYY (M3/S2)
        write(EE_UNIT) (REAL(WVHUV_Global(L,KC),4), L = 2,LA_Global)        ! *** SXY (M3/S2)
      endif
    endif
  endif

  if( ISTRAN(1) >= 1 ) WRITE(EE_UNIT) ((REAL(SAL_Global(L,K),4), K = KSZ_Global(L),KC), L = 2,LA_Global)

  if( ISTRAN(2) >= 1 )then
    write(EE_UNIT) ((REAL(TEM_Global(L,K),4), K = KSZ_Global(L),KC), L = 2,LA_Global)
    if( TEMBO > 0.) WRITE(EE_UNIT) (REAL(TEMB_Global(L),4), L = 2,LA_Global)
    if( IEVAP > 1 )then
      write(EE_UNIT) (REAL(EVAPT_Global(L),4), L = 2,LA_Global)
      write(EE_UNIT) (REAL(RAINT_Global(L),4), L = 2,LA_Global)
    endif
  endif

  if( ISTRAN(3) >= 1 ) WRITE(EE_UNIT) (((REAL(DYE_Global(L,K,MD),4), K = KSZ_Global(L),KC), L = 2,LA_Global), MD = 1,NDYE)

  if( ISTRAN(4) >= 1 ) WRITE(EE_UNIT) ((REAL(SFL_Global(L,K),4), K = KSZ_Global(L),KC), L = 2,LA_Global)

  if( ISTRAN(5) >= 1 )then
    write(EE_UNIT) ((REAL(TOXB_Global(L,KBT_Global(L),NT),4), L = 2,LA_Global), NT = 1,NTOX)
    
    write(EE_UNIT) (((REAL(TOX_Global(L,K,NT),4), K = KSZ_Global(L),KC), L = 2,LA_Global), NT = 1,NTOX)
  endif

  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
    write(EE_UNIT) (REAL(BELV_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(HBED_Global(L,KBT_Global(L)),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(BDENBED_Global(L,KBT_Global(L)),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(PORBED_Global(L,KBT_Global(L)),4), L = 2,LA_Global)
    if( ISTRAN(6) >= 1 ) WRITE(EE_UNIT) ((REAL(SEDB_Global(L,KBT_Global(L),NS),4), L = 2,LA_Global), NS = 1,NSED)
    if( ISTRAN(7) >= 1 ) WRITE(EE_UNIT) ((REAL(SNDB_Global(L,KBT_Global(L),NX),4), L = 2,LA_Global), NX = 1,NSND)
    write(EE_UNIT) ((REAL(0.0,4), L = 2,LA_Global), NS = 1,NSED+NSND)
    if( ISTRAN(6) >= 1 ) WRITE(EE_UNIT) (((REAL(SED_Global(L,K,NS),4), K = KSZ_Global(L),KC), L = 2,LA_Global), NS = 1,NSED2)
    if( ISTRAN(7) >= 1 )then
      write(EE_UNIT) (((REAL(SND_GLobal(L,K,NX),4), K = KSZ_Global(L),KC), L = 2,LA_Global), NX = 1,NSND)
      if( ICALC_BL > 0 .and. NSND > 0 )then
        ! *** CALCULATE EQUIVALENT CONCENTRATIONS
        write(EE_UNIT) ((REAL(QSBDLDX_Global(L,NX),4), L = 2,LA_Global), NX = 1,NSND)
        write(EE_UNIT) ((REAL(QSBDLDY_Global(L,NX),4), L = 2,LA_Global), NX = 1,NSND)
      endif
    endif
  elseif( BATHY.IFLAG > 0 )then
      write(EE_UNIT) (REAL(BELV_Global(L),4), L = 2,LA_Global)
  endif
  if( ISGWIE > 0 )then
    write(EE_UNIT) (REAL(EVAPSW_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(EVAPGW_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(QGW_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(AGWELV_Global(L),4), L = 2,LA_Global)
  endif

  if( ISTRAN(2) > 0 .and. ISICE  >= 3 )then
    write(EE_UNIT) (REAL(ICETHICK_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(ICETEMP_Global(L),4), L = 2,LA_Global)
  endif

  FLUSH(EE_UNIT)
  close(EE_UNIT,STATUS = 'KEEP')
  !PRINT *, 'CLOSE: NITER = ', NITER
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_WC.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif

END SUBROUTINE

INTEGER FUNCTION BLOCKWC(CELL3D)

  INTEGER:: DSIZE, CELL2D, CELL3D

  CELL2D = LA_Global - 1

  DSIZE = 2 !**EETIME
  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    DSIZE = DSIZE + CELL2D                      ! TOP LAYER
  endif

  ! *** READ THE WATER COLUMN AND TOP LAYER OF SEDIMENT DATA, IF NEEDED
  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    if( LSEDZLJ )then
      DSIZE = DSIZE + CELL2D                  ! SHEAR
    elseif( ISBEDSTR >= 1 )then
      DSIZE = DSIZE + CELL2D                  ! TAUBSED
      if( ISBEDSTR == 1 )then
        DSIZE = DSIZE + CELL2D                ! TAUBSND
      endif
    else
      ! *** TOTAL BED SHEAR STRESS
      DSIZE = DSIZE + CELL2D                  ! TAUB
    endif
  else
    ! *** TOTAL BED SHEAR STRESS
    DSIZE = DSIZE + CELL2D                    ! SHEAR
  endif
  if( ISWAVE >= 1 )then
    DSIZE = DSIZE + 2*CELL2D                  ! QQWV3, SHEAR
    if( ISWAVE >= 3 )then
      DSIZE = DSIZE + 3*CELL2D                ! WV(L).HEIGHT, WV.FREQ
      if( ISWAVE == 4 )then
        DSIZE = DSIZE + 4*CELL2D              ! WV.DISSIPA, WVHUU, WVHVV, WVHUV
      endif
    endif
  endif
  if( ISTRAN(1) >= 1 ) DSIZE = DSIZE + CELL3D   ! SAL_Global
  if( ISTRAN(2) >= 1 )then
    DSIZE = DSIZE + CELL3D                      ! TEM_Global
    if( TEMBO > 0.) DSIZE = DSIZE + CELL2D      ! TEMB
    if( IEVAP > 1 )then                         
      DSIZE = DSIZE + 2 * CELL2D                ! EVAPT, RAINT
    endif
  endif
  if( ISTRAN(3) >= 1 ) DSIZE = DSIZE + NDYE*CELL3D   ! DYE
  if( ISTRAN(4) >= 1 ) DSIZE = DSIZE + CELL3D        ! SFL
  if( ISTRAN(5) >= 1 )then
    DSIZE = DSIZE + NTOX*CELL2D               ! TOXB
    DSIZE = DSIZE + NTOX*CELL3D               ! TOX
  endif
  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
    DSIZE = DSIZE + 4*CELL2D                  ! BELV, HBED, BDENBED, PORBED
    if( ISTRAN(6) == 1 )then
      DSIZE = DSIZE + 2*NSED*CELL2D           ! SEDB, VFRBED
      DSIZE = DSIZE + NSED2*CELL3D            ! SED
    endif
    if( ISTRAN(7) >= 1 )then
      DSIZE = DSIZE + 2*NSND*CELL2D           ! SNDB, VFRBED
      DSIZE = DSIZE + NSND*CELL3D             ! SND
      if( ICALC_BL > 0 .and. NSND > 0 )then
        DSIZE = DSIZE + 2*NSND*CELL2D         ! QSBDLDX, QSBDLDY
      endif
    endif
  elseif( BATHY.IFLAG > 0 )then
      DSIZE = DSIZE + CELL2D                  ! BELV
  endif
  ! *** NEW BLOCK
  if( ISGWIE > 0 )then
    DSIZE = DSIZE + 4*CELL2D                  ! EVAPSW, EVAPGW, QGW, AGWELV
  endif
  if(ISTRAN(2) > 0 .and. ISICE  >= 3 )then
    DSIZE = DSIZE + 2*CELL2D                  ! ICETHICK, ICETEMP
  endif
  BLOCKWC = 4*DSIZE                             ! SINGLE PRECISION
  END FUNCTION

  SUBROUTINE SEDZLJOUT

  integer(IK4) :: VER, HSIZE, BSIZE, ISTAT
  integer(IK4) :: L, NS, NT, K, ITMP
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP, SURFACE
  real(RKD)    :: PTIME

  NBEDSTEPS = NBEDSTEPS + 1
  if( NBEDSTEPS >= ISBEXP )then
    FILENAME = OUTDIR//'EE_SEDZLJ.OUT'
    EE_UNIT = 99                        ! *** EE_SEDZLJ.OUT
    if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_SEDZLJ == 0 )then
      RSTFIRST_SEDZLJ = 1
      write(*,'(A)')'READING TO STARTING TIME FOR SEDZLJ'
      FSIZE = FILESIZE(FILENAME)
      open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
      read(EE_UNIT) VER,HSIZE,BSIZE
      if( VER /= 12000 )then
        write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
        call STOPP('.')
      endif
      OFFSET = HSIZE

      NS = 0
      do while(OFFSET < FSIZE)
#ifdef GNU
        call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

        if(IS_IOSTAT_END(ISTAT)) then
          write(EE_UNIT) EETIME
          exit
        endif
#else
        ISTAT = FSEEK(EE_UNIT,OFFSET,0)

        if( EOF(EE_UNIT) )then
          write(EE_UNIT) EETIME
          exit
        endif
#endif 
        if( ISTAT /= 0 ) exit

        read(EE_UNIT) PTIME

        if( DEBUG )then
          NS = NS+1
          if( NS == 8 )then
            write(*,'(F10.3)')PTIME
            NS = 0
          else
            WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
          endif
        endif
        if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

        OFFSET = OFFSET + BSIZE
      enddo
      if( DEBUG ) WRITE(*,'(" ")')
      write(*,'(A)')'FINISHED READING SEDZLJ'

    else
      NERR = 0
      IORIGIN = 5
      100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
      write(EE_UNIT) EETIME
    endif
    if( DTSED == 0. ) DTSED = DT

    ! *** This is used to write out the active layers instead of using the LAYERACTIVE array. 
    do L = 2,LA_Global
      do K = 1,KB
        if( TSED_GLOBAL(K,L) > 1E-8 )then
          write(EE_UNIT) INT(1,4)                                                    ! *** LAYERACTIVE(KB,LCM) - This is = 1 when a bed layer (KB index) exists with mass
        else
          write(EE_UNIT) INT(0,4)
        endif
      enddo
    enddo
  
    ! *** REAL*4  - Global arrays use reversed indicies for MPI mapping.  Writes the same order as OMP original
    write(EE_UNIT) (REAL(TAU_Global(L),4), L = 2,LA_Global)                                      ! *** TAU(LCM)      - Shear Stress in dynes/cm^2
    write(EE_UNIT) ((REAL(TSED_Global(K,L),4), K = 1,KB), L = 2,LA_Global)                       ! *** TSED(KB,LCM)  - This is the mass in g/cm^2 in each layer
    write(EE_UNIT) ((REAL(BULKDENS_Global(K,L),4), K = 1,KB), L = 2,LA_Global)                   ! *** BULKDENS(KB,LCM) - Dry Bulk density of each layer (g/cm^3)
    write(EE_UNIT) (((REAL(PERSED_Global(NS,K,L),4), K = 1,KB), L = 2,LA_Global), NS = 1,NSEDS)  ! *** PERSED(NSEDS,KB,LCM) - This is the mass percentage of each size class in a layer
    write(EE_UNIT) (REAL(D50AVG_Global(L),4), L = 2,LA_Global)                                   ! *** D50AVG(LCM)   - Average particle size of bed surface (microns)
                                                                                                 
    write(EE_UNIT) ((REAL(ERO_SED_FLX_Global(L,NS),4), L = 2,LA_Global), NS = 1,NSEDS2)              ! *** ERO_SED_FLX(LCM)    - Total erosion rate in the cell g/cm^2/s
    write(EE_UNIT) ((REAL(DEP_SED_FLX_Global(L,NS),4), L = 2,LA_Global), NS = 1,NSEDS2)              ! *** DEP_SED_FLX(LCM)     - Total deposition rate in the cell g/cm^2/s
                                                                                      
    if( ICALC_BL > 0 )then                                                              
      write(EE_UNIT) ((REAL(CBL_Global(L,NS)*10000.,4), L = 2,LA_Global), NS = 1,NSEDS)          ! *** Bedload mass inb g/m^2 for each size class
      write(EE_UNIT) ((REAL(QSBDLDX_Global(L,NS),4), L = 2,LA_Global), NS = 1,NSEDS)             ! *** Bedload flux in X direction (g/s)
      write(EE_UNIT) ((REAL(QSBDLDY_Global(L,NS),4), L = 2,LA_Global), NS = 1,NSEDS)             ! *** Bedload flux in Y direction (g/s)
    endif

    if( ISTRAN(5) > 0 )then
      do NT = 1,NTOX
        do L = 2,LA_Global
          do K = 1,KB
            if( K < 3 .and. TSED_Global(K,L) > 0. .and. HBED_Global(L,KBT_Global(L)) > 0. )then  ! *** avoid division by HBED = 0 --- DKT
              TMP = 0.01*TSED_Global(K,L)/BULKDENS_Global(K,L)             ! *** HBED(L,K)
              TMP = TMP/HBED_Global(L,KBT_Global(L))
              TMP = TMP*TOXB_Global(L,KBT_Global(L),NT)
              write(EE_UNIT) TMP
            elseif( TSED_Global(K,L) > 0. )then
              write(EE_UNIT) REAL(TOXB_Global(L,K,NT),4)
            else
              write(EE_UNIT) 0.0_4
            endif
          enddo
        enddo
      enddo

      if( ICALC_BL > 0 )then
        write(EE_UNIT) ((REAL(CBLTOX_Global(L,NT),4), L = 2,LA_Global), NT = 1,NTOX)    ! *** Bedload toxic concentration (mg/m^2)
      endif
    endif

    FLUSH(EE_UNIT)
    close(EE_UNIT,STATUS = 'KEEP')
    
    NBEDSTEPS = 0
  endif
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_SEDZLJ.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif
  
END SUBROUTINE

INTEGER FUNCTION BLOCKSEDZLJ(CELL3D)

  INTEGER:: DSIZE, CELL2D, CELL3D

  CELL2D = LA_Global - 1
  DSIZE = 2                             ! *** EETIME

  ! *** INTEGER*4
  DSIZE = DSIZE + KB*CELL2D             ! *** LAYERACTIVE

  ! *** REAL*4
  DSIZE = DSIZE + 2*CELL2D              ! *** TAU,D50AVG
  DSIZE = DSIZE + 2*KB*CELL2D           ! *** TSED,BULKDENS
  DSIZE = DSIZE + KB*NSCM*CELL2D        ! *** PERSED
  DSIZE = DSIZE + 2*NSEDS2*CELL2D       ! *** ERO_SED_FLX, DEP_SED_FLX
  if( ICALC_BL > 0 )then
    DSIZE = DSIZE + 3*NSEDS*CELL2D      ! *** CBL,QSBDLDX,QSBDLDY
  endif
  if( ISTRAN(5) >= 1 )then
    DSIZE = DSIZE + NTOX*KB*CELL2D      ! *** TOXB
    if( ICALC_BL > 0 )then
      DSIZE = DSIZE + NTOX*CELL2D       ! *** CBLTOX
    endif
  endif
  BLOCKSEDZLJ = 4*DSIZE                 ! *** CONVERT TO BYTES
  END FUNCTION

  SUBROUTINE BEDOUT

  integer(IK4) :: VER, HSIZE, BSIZE, ISTAT
  integer(IK4) :: I, NS, L, K, NX, NT, ITMP
  integer(IK8) :: FSIZE,  OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME

  ! *** SEDIMENT BED LAYERS
  NBEDSTEPS = NBEDSTEPS + 1
  if( NBEDSTEPS >= ISBEXP )then
    FILENAME = OUTDIR//'EE_BED.OUT'
    if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_BED == 0 )then
      RSTFIRST_BED = 1
      write(*,'(A)')'READING TO STARTING TIME FOR BED'
      FSIZE = FILESIZE(FILENAME)
      open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
      read(EE_UNIT) VER,HSIZE,BSIZE
      if( VER /= 8400 )then
        write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
        call STOPP('.')
      endif
      OFFSET = HSIZE

      NS = 0
      do while(OFFSET < FSIZE)
#ifdef GNU
        call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

        if(IS_IOSTAT_END(ISTAT)) then
          write(EE_UNIT) EETIME
          exit
        endif
#else
        ISTAT = FSEEK(EE_UNIT,OFFSET,0)

        if( EOF(EE_UNIT) )then
          write(EE_UNIT) EETIME
          exit
        endif
#endif 
        if( ISTAT /= 0 ) exit

        read(EE_UNIT) PTIME

        if( DEBUG )then
          NS = NS+1
          if( NS == 8 )then
            write(*,'(F10.3)')PTIME
            NS = 0
          else
            WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
          endif
        endif
        if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

        OFFSET = OFFSET + BSIZE
      enddo
      if( DEBUG ) WRITE(*,'(" ")')
      write(*,'(A)')'FINISHED READING BED'
    else
      NERR = 0
      IORIGIN = 4
      100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
      write(EE_UNIT) EETIME
    endif

    write(EE_UNIT) (INT(KBT_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) ((REAL(HBED_Global(L,K),4), K = 1,KB), L = 2,LA_Global)
    write(EE_UNIT) ((REAL(BDENBED_Global(L,K),4), K = 1,KB), L = 2,LA_Global)
    write(EE_UNIT) ((REAL(PORBED_Global(L,K),4), K = 1,KB), L = 2,LA_Global)

    if( ISTRAN(6) >= 1 )then
      do NS = 1,NSED
        write(EE_UNIT) ((REAL(SEDB_Global(L,K,NS),4), K = 1,KB), L = 2,LA_Global)
      enddo
    endif
    if( ISTRAN(7) >= 1 )then
      do NX = 1,NSND
        NS = NSED+NX
        write(EE_UNIT) ((REAL(SNDB_Global(L,K,NX),4), K = 1,KB), L = 2,LA_Global)
      enddo
    endif
    if( ISTRAN(5) >= 1 )then
      do NT = 1,NTOX
        write(EE_UNIT) ((REAL(TOXB_Global(L,K,NT),4), K = 1,KB), L = 2,LA_Global)
      enddo
    endif
    FLUSH(EE_UNIT)
    close(EE_UNIT,STATUS = 'KEEP')
    
    NBEDSTEPS = 0
  endif
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_BED.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif

END SUBROUTINE

INTEGER FUNCTION BLOCKBED(CELL3D)
  
  INTEGER:: DSIZE, CELL2D, CELL3D
  
  CELL2D = LA_Global - 1
  DSIZE = 2 !**EETIME
  DSIZE = DSIZE + CELL2D
  DSIZE = DSIZE + 3*KB*CELL2D
  if( ISTRAN(6) >= 1 )then
    DSIZE = DSIZE + NSED*KB*CELL2D
  endif
  if( ISTRAN(7) >= 1 )then
    DSIZE = DSIZE + NSND*KB*CELL2D
  endif
  if( ISTRAN(5) >= 1 )then
    DSIZE = DSIZE + NTOX*KB*CELL2D
  endif
  BLOCKBED = 4*DSIZE
  
  END FUNCTION

  SUBROUTINE SDOUT(JSEXPLORER)

  integer(IK4) :: VER,HSIZE,BSIZE,ISTAT
  integer(IK4) :: L,K,JSEXPLORER,ITMP,NS
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME

  NSEDSTEPS = NSEDSTEPS+1
  if( NSEDSTEPS >= ABS(ISSDBIN) .or. JSEXPLORER == 1 )then

    FILENAME = OUTDIR//'EE_SD.OUT'
    EE_UNIT = 101                      ! *** EE_SD.OUT
    if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_SD == 0 )then
      RSTFIRST_SD = 1
      write(*,'(A)')'READING TO STARTING TIME FOR SD'
      FSIZE = FILESIZE(FILENAME)
      open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
      read(EE_UNIT) VER,HSIZE,BSIZE
      if( VER /= 8400 )then
        write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
        call STOPP('.')
      endif
      OFFSET = HSIZE

      NS = 0
      do while(OFFSET < FSIZE)
#ifdef GNU
        call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

        if(IS_IOSTAT_END(ISTAT)) then
          write(EE_UNIT) EETIME
          exit
        endif
#else
        ISTAT = FSEEK(EE_UNIT,OFFSET,0)

        if( EOF(EE_UNIT) )then
          write(EE_UNIT) EETIME
          exit
        endif
#endif 
        if( ISTAT /= 0 ) exit

        read(EE_UNIT) PTIME

        if( DEBUG )then
          NS = NS+1
          if( NS == 8 )then
            write(*,'(F10.3)')PTIME
            NS = 0
          else
            WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
          endif
        endif
        if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

        OFFSET = OFFSET + BSIZE
      enddo
      if( DEBUG ) WRITE(*,'(" ")')
      write(*,'(A)')'FINISHED READING SD'

    else
      NERR = 0
      IORIGIN = 7
      100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
      write(EE_UNIT) EETIME
    endif

    write(EE_UNIT) ((REAL(SMPON_Global(L,K),4), L = 2,LA_Global), K = 1,3)
    write(EE_UNIT) ((REAL(SMPOP_Global(L,K),4), L = 2,LA_Global), K = 1,3)
    write(EE_UNIT) ((REAL(SMPOC_Global(L,K),4), L = 2,LA_Global), K = 1,3)
    write(EE_UNIT) ((REAL(SMDFN_Global(L,K),4), L = 2,LA_Global), K = 1,3)
    write(EE_UNIT) ((REAL(SMDFP_Global(L,K),4), L = 2,LA_Global), K = 1,3)
    write(EE_UNIT) ((REAL(SMDFC_Global(L,K),4), L = 2,LA_Global), K = 1,3)
    write(EE_UNIT) (REAL(SM1NH4_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM2NH4_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM1NO3_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM2NO3_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM1PO4_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM2PO4_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM1H2S_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM2H2S_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM1SI_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SM2SI_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SMPSI_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SMBST_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SMT_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SMCSOD_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(SMNSOD_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(WQBFNH4_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(WQBFNO3_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(WQBFO2_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(WQBFCOD_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(WQBFPO4D_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) (REAL(WQBFSAD_Global(L),4), L = 2,LA_Global)

    FLUSH(EE_UNIT)
    close(EE_UNIT,STATUS = 'KEEP')
    NSEDSTEPS = 0
  endif
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_SD.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif
  
  END SUBROUTINE

  SUBROUTINE RPEMOUT(JSEXPLORER)

  integer(IK4) :: VER,HSIZE,BSIZE,ISTAT
  integer(IK4) :: L,JSEXPLORER,NS,ITMP
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME
  logical(4)   :: LMASK

  ! *** RPEM
  if( ISRPEM > 0 )then
    ! *** IF JSEXPLORER = 1 THEN WRITE THE ARRAYS (I.E. IC'S)
    NRPEMSTEPS = NRPEMSTEPS+1

    if( NRPEMSTEPS >= NRPEMEE .or. JSEXPLORER == 1 )then
      FILENAME = OUTDIR//'EE_RPEM.OUT'
      EE_UNIT = 102                      ! *** EE_RPEM.OUT
      if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_RPEM == 0 )then
        RSTFIRST_RPEM = 1
        write(*,'(A)')'READING TO STARTING TIME FOR RPEM'
        FSIZE = FILESIZE(FILENAME)
        open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
        read(EE_UNIT) VER,HSIZE,BSIZE
        if( VER /= 8400 )then
          write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
          call STOPP('.')
        endif
        OFFSET = HSIZE

        NS = 0
        do while(OFFSET < FSIZE)
#ifdef GNU
          call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

          if(IS_IOSTAT_END(ISTAT)) then
            write(EE_UNIT) EETIME
            exit
          endif
#else
          ISTAT = FSEEK(EE_UNIT,OFFSET,0)

          if( EOF(EE_UNIT) )then
            write(EE_UNIT) EETIME
            exit
          endif
#endif 
          if( ISTAT /= 0 ) exit

          read(EE_UNIT) PTIME

          if( DEBUG )then
            NS = NS+1
            if( NS == 8 )then
              write(*,'(F10.3)')PTIME
              NS = 0
            else
              WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
            endif
          endif
          if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

          OFFSET = OFFSET + BSIZE
        enddo
        if( DEBUG ) WRITE(*,'(" ")')
        write(*,'(A)')'FINISHED READING RPEM'

      else
        NERR = 0
        IORIGIN = 9
        100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
        write(EE_UNIT) EETIME
      endif

      do L = 2,LA_Global
        if( LMASKRPEM_Global(L) .or. INITRPEM == 0 )then
          write(EE_UNIT) REAL(WQRPS_Global(L),4)
        endif
      enddo
      do L = 2,LA_Global
        if( LMASKRPEM_Global(L) .or. INITRPEM == 0 )then
          write(EE_UNIT) REAL(WQRPR_Global(L),4)
        endif
      enddo
      do L = 2,LA_Global
        if( LMASKRPEM_Global(L) .or. INITRPEM == 0 )then
          write(EE_UNIT) REAL(WQRPE_Global(L),4)
        endif
      enddo
      do L = 2,LA_Global
        if( LMASKRPEM_Global(L) .or. INITRPEM == 0 )then
          write(EE_UNIT) REAL(WQRPD_Global(L),4)
        endif
      enddo

      FLUSH(EE_UNIT)
      close(EE_UNIT,STATUS = 'KEEP')
      NRPEMSTEPS = 0
    endif
  endif
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_RPEM.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif
  
END SUBROUTINE

SUBROUTINE WQOUT

  integer(IK4) :: VER,HSIZE,BSIZE,ISTAT
  integer(IK4) :: L,K,MW,NW,ITMP,NS
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: WQ
  real(RKD)    :: PTIME

  ! *** ISWQLVL = 0
  ! ***  1) CHC - cyanobacteria
  ! ***  2) CHD - diatom algae
  ! ***  3) CHG - green algae
  ! ***  4) ROC - refractory particulate organic carbon
  ! ***  5) LOC - labile particulate organic carbon
  ! ***  6) DOC - dissolved organic carbon
  ! ***  7) ROP - refractory particulate organic phosphorus
  ! ***  8) LOP - labile particulate organic phosphorus
  ! ***  9) DOP - dissolved organic phosphorus
  ! *** 10) P4D - total phosphate
  ! *** 11) RON - refractory particulate organic nitrogen 22) macroalgae
  ! *** 12) LON - labile particulate organic nitrogen
  ! *** 13) DON - dissolved organic nitrogen
  ! *** 14) NHX - ammonia nitrogen
  ! *** 15) NOX - nitrate nitrogen
  ! *** 16) SUU - particulate biogenic silica
  ! *** 17) SAA - dissolved available silica
  ! *** 18) COD - chemical oxygen demand
  ! *** 19) DOX - dissolved oxygen
  ! *** 20) TAM - total active metal
  ! *** 21) FCB - fecal coliform bacteria

  ! *** ISWQLVL = 1
  ! ***  1) ROC - Refractory particulate organic carbon
  ! ***  2) LOC - Labile particulate organic carbon
  ! ***  3) DOC - Dissolved organic carbon
  ! ***  4) ROP - Refractory particulate organic phosphorus
  ! ***  5) LOP - Labile particulate organic phosphorus
  ! ***  6) DOP - Dissolved organic phosphorus
  ! ***  7) P4D - Total phosphate
  ! ***  8) RON - Refractory particulate organic nitrogen
  ! ***  9) LON - Labile particulate organic nitrogen
  ! *** 10) DON - Dissolved organic nitrogen
  ! *** 11) NHX - Ammonia nitrogen
  ! *** 12) NOX - Nitrate nitrogen
  ! *** 13) SUU - Particulate biogenic silica
  ! *** 14) SAA - Dissolved available silica
  ! *** 15) COD - Chemical oxygen demand
  ! *** 16) DOX - Dissolved oxygen
  ! *** 17) TAM - Total active metal
  ! *** 18) FCB - Fecal coliform bacteria
  ! *** 19) CO2 - Dissolved carbon dioxide
  ! *** 20) BIO - Any number of biota classes (NALGAE), phytoplankton, macrophytes and zooplankton
  ! ***           Can include green, blue-green, cyanobacteria, etc.  

  if( ISRESTI == 1 )then
    IWQ = 0
    do MW = 1,NWQV
      IWQ(MW) = ISKINETICS(MW)
    enddo
  endif

  FILENAME = OUTDIR//'EE_WQ.OUT'
  EE_UNIT = 100                      ! *** EE_WQ.OUT
  if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_WQ == 0 )then
    RSTFIRST_WQ = 1
    write(*,'(A)')'READING TO STARTING TIME FOR WQ'
    FSIZE = FILESIZE(FILENAME)
    open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
    read(EE_UNIT) VER,HSIZE,BSIZE
    if( VER /= 10300 )then
      write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
      call STOPP('.')
    endif
    OFFSET = HSIZE

    NS = 0
    do while(OFFSET < FSIZE)
#ifdef GNU
      call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

      if(IS_IOSTAT_END(ISTAT)) then
        write(EE_UNIT) EETIME
        exit
      endif
#else
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      if( EOF(EE_UNIT) )then
        write(EE_UNIT) EETIME
        exit
      endif
#endif 
      if( ISTAT /= 0 ) exit

      read(EE_UNIT) PTIME

      if( DEBUG )then
        NS = NS+1
        if( NS == 8 )then
          write(*,'(F10.3)')PTIME
          NS = 0
        else
          WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
        endif
      endif
      if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

      OFFSET = OFFSET + BSIZE
    enddo
    if( DEBUG ) WRITE(*,'(" ")')
    write(*,'(A)')'FINISHED READING WQ'

  else
    NERR = 0
    IORIGIN = 6
    100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
    write(EE_UNIT) EETIME
  endif

  if( ISWQLVL == 0 )then
    ! *** Starts from ROC (4 - 22) as algae are moved to the end
    do NW = 4,22 
      if( IWQ(NW) > 0 )then
        do L = 2,LA_Global
          do K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            write(EE_UNIT) WQ
          enddo
        enddo
      endif
    enddo
    
    ! *** Algae (1,2,3) are now moved to the end
    do NW = 1,3
      if( IWQ(NW) > 0 )then
        do L = 2,LA_Global
          do K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            write(EE_UNIT) WQ
          enddo
        enddo
      endif
    enddo
    
    ! *** 23) macroalgae
    if( nfixed > 0 )then
      NW = 23 
      if( IWQ(NW) > 0 )then
        do L = 2,LA_Global
          do K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            write(EE_UNIT) WQ
          enddo
        enddo
      endif
    endif
    
    ! *** Zooplankton
    if( IWQZPL > 0 )then
      do NW = 1,NZOOPL
        do L = 2,LA_Global
          do K = KSZ_Global(L),KC
            if( NFIXED > 0 )then
              WQ = WQV_Global(L,K,22+NW)
            else
              WQ = WQV_Global(L,K,21+NW)
            endif
            write(EE_UNIT) WQ
          enddo
        enddo
      enddo
    endif  
    
  else
    ! *** New DSI standard WQ
    do NW = 1,NWQV
      if( IWQ(NW) > 0 )then
        do L = 2,LA_Global
          do K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            write(EE_UNIT) WQ
          enddo
        enddo
      endif
    enddo
  endif
  
  FLUSH(EE_UNIT)
  close(EE_UNIT,STATUS = 'KEEP')

  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_WQ.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif
  
END SUBROUTINE

SUBROUTINE BCOUT
  ! *** OUTPUT QSUM AND OPEN BC FLOWS
  integer(IK4) :: VER, HSIZE, BSIZE, ISTAT
  integer(IK4) :: L, K, IBC, LL, LQ, LE, LN, ITMP, NS, NWR
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME

  FILENAME = OUTDIR//'EE_BC.OUT'
  EE_UNIT = 104                      ! *** EE_BC.OUT

  if( ISRESTI /= 0 .and. ICONTINUE == 1 .and. RSTFIRST_BC == 0 )then
    RSTFIRST_BC = 1
    write(*,'(A)')'READING TO STARTING TIME FOR BC'
    FSIZE = FILESIZE(FILENAME)
    open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)
    read(EE_UNIT) VER,HSIZE,BSIZE
    if( VER /= 8400 )then
      write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
      call STOPP('.')
    endif
    OFFSET = HSIZE

    NS = 0
    do while(OFFSET < FSIZE)
#ifdef GNU
      call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

      if(IS_IOSTAT_END(ISTAT)) then
        write(EE_UNIT) EETIME
        exit
      endif
#else
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      if( EOF(EE_UNIT) )then
        write(EE_UNIT) EETIME
        exit
      endif
#endif 
      if( ISTAT /= 0 ) exit

      read(EE_UNIT) PTIME

      if( DEBUG )then
        NS = NS+1
        if( NS == 8 )then
          write(*,'(F10.3)')PTIME
          NS = 0
        else
          WRITE(*,'(F10.3)',ADVANCE='NO') PTIME
        endif
      endif
      if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

      OFFSET = OFFSET + BSIZE
    enddo
    if( DEBUG ) WRITE(*,'(" ")')
    write(*,'(A)')'FINSIHED READING BC'

  else
    NERR = 0
    IORIGIN = 0
    100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
    write(EE_UNIT) EETIME
  endif

  if( ISRESTI == 0 .or. ICONTINUE == 0) NRESTART = 0

  ! *** GET AVERAGE FLOWS FOR SNAPSHOT
  ! @todo put bcout_average routine back in and make it worth with MPI
  !CALL BCOUT_AVERAGE

  ! *** OUTPUT SELECTIVE QSUM
  if( KC > 1 )then

    write(EE_UNIT) (REAL(QSUM_Global(L,KC),4), L = 2,LA_Global)

    do IBC = 1, NBCS
      L = LBCS(IBC) !
      write(EE_UNIT) (REAL(QSUM_Global(L,K),4),K = 1,KC)
    enddo
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      write(EE_UNIT) (REAL(QSUM_Global(L, KSZ_Global(L)),4), L = 2,LA_Global)
    endif
  else
    ! *** SINGLE LAYER
    do L = 2,LA_Global
      write(EE_UNIT) REAL(QSUME_Global(L),4)
    enddo
  endif

  ! ***  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  do LL = 1, NPBS
    LQ = LPBS(LL)
    LN = LNC_Global(LQ)
    write(EE_UNIT) (REAL(VHDX2_Global(LN,K),4),K = 1,KC)
  enddo
  do LL = 1, NPBW
    LQ = LPBW(LL)
    LE = LEC_Global(LQ)
    write(EE_UNIT) (REAL(UHDY2_Global(LE,K),4),K = 1,KC)
  enddo
  do LL = 1, NPBE
    LQ = LPBE(LL)
    write(EE_UNIT) (REAL(UHDY2_Global(LQ,K),4),K = 1,KC)
  enddo
  do LL = 1, NPBN
    LQ = LPBN(LL)
    write(EE_UNIT) (REAL(VHDX2_Global(LQ,K),4),K = 1,KC)
  enddo

  ! @todo verify this will work for MPI
  ! *** save FOR INFORMATION ONLY.   QCTLT AND HSCTL.CUR.FLOW ALREADY INCLUDED IN QSUM
  if( NQCTL > 0 )then

    write(EE_UNIT) (((REAL(QCTLT(K,L,NS),4), K = 1,KC), NS = 1,2), L = 1,NQCTL) !< do not need a global one since this was not decomposed

    if(NQCTLSER > 0 .or. NQCTRULES > 0 )then
      write(EE_UNIT) (INT(HSCTL(L).CUR.STATE), L = 1,NQCTL)
      write(EE_UNIT) (INT(1), L = 1,NQCTL)                     ! HSCTL(L).CUR.UNITS not used
      write(EE_UNIT) (INT(HSCTL(L).CUR.ID), L = 1,NQCTL)
      write(EE_UNIT) (REAL(HSCTL(L).CUR.FLOW,4), L = 1,NQCTL)
      write(EE_UNIT) (REAL(HSCTL(L).CUR.HEIGHT,4), L = 1,NQCTL)
      write(EE_UNIT) (REAL(HSCTL(L).CUR.WIDTH,4), L = 1,NQCTL)
      write(EE_UNIT) (REAL(HSCTL(L).CUR.SILL,4), L = 1,NQCTL)
    endif
  endif

  ! @todo verify this will work for MPI
  ! *** save FOR INFORMATION ONLY.   WITH_RET_CTL ALREADY INCLUDED IN QSUM
  if( NQWR > 0 )then
    write(EE_UNIT) (INT(WITH_RET_CTL(NWR).CUR.STATE),   NWR = 1,NQWR)
    write(EE_UNIT) (REAL(WITH_RET_CTL(NWR).CUR.FLOW,4), NWR = 1,NQWR)
  endif

  close(EE_UNIT,STATUS = 'KEEP')

  !SUMTIME = 0.
  !BCQSUM = 0.
  !BCQSUME = 0.
  !BCUHDY2 = 0.
  !BCVHDX2 = 0
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_BC.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif

  END SUBROUTINE

  SUBROUTINE BCOUT_ACCUMULATE(DELT)
  ! *** ACCUMULATE BOUNDARY INFLOWS/OUTFLOWS
  real,intent(IN) :: DELT
  integer(IK4) :: L,K,LL,LQ,LE,LN,IBC

  if( INITBCOUT == 0 ) CALL BCOUT_INITIALIZE

  SUMTIME = SUMTIME + DELT

  ! *** SUM SELECTIVE QSUM
  if( KC > 1 )then
    do L = 2,LA_Global
      BCQSUM(L,KC) = BCQSUM(L,KC) + QSUM(L,KC)*DELT
    enddo
    do IBC = 1,NBCS
      L = LBCS(IBC)
      do K = 1,KS
        BCQSUM(L,K) = BCQSUM(L,K) + QSUM(L,K)*DELT
      enddo
    enddo
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 2,LA_Global
        K = KSZ_Global(L)
        BCQSUM(L,K) = BCQSUM(L,K) + QSUM(L,K)*DELT
      enddo
    endif
  else
    ! *** SINGLE LAYER
    do L = 2,LA_Global
      BCQSUME(L) = BCQSUME(L) + QSUME(L)*DELT
    enddo
  endif

  ! ***  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  do LL = 1,NPBS
    LQ = LPBS(LL)
    LN = LNC_Global(LQ)
    do K = 1,KC
      BCVHDX2(LN,K) = BCVHDX2(LN,K) + VHDX2(LN,K)*DELT
    enddo
  enddo
  do LL = 1,NPBW
    LQ = LPBW(LL)
    LE = LEC_Global(LQ)
    do K = 1,KC
      BCUHDY2(LE,K) = BCUHDY2(LE,K) + UHDY2(LE,K)*DELT
    enddo
  enddo
  do LL = 1,NPBE
    LQ = LPBE(LL)
    do K = 1,KC
      BCUHDY2(LQ,K) = BCUHDY2(LQ,K) + UHDY2(LQ,K)*DELT
    enddo
  enddo
  do LL = 1,NPBN
    LQ = LPBN(LL)
    do K = 1,KC
      BCVHDX2(LQ,K) = BCVHDX2(LQ,K) + VHDX2(LQ,K)*DELT
    enddo
  enddo

  return

  END SUBROUTINE

  SUBROUTINE BCOUT_AVERAGE
  ! *** COMPUTE AVEARGE BOUNDARY INFLOWS/OUTFLOWS

  integer(IK4) :: L,K,LL,LQ,LE,LN,IBC

  if( SUMTIME <= 0. ) return

  ! *** SUM SELECTIVE QSUM
  if( KC > 1 )then
    do L = 2,LA_Global
      BCQSUM(L,KC) = BCQSUM(L,KC)/SUMTIME
    enddo
    do IBC = 1,NBCS
      L = LBCS(IBC)
      do K = 1,KS
        BCQSUM(L,K) = BCQSUM(L,K)/SUMTIME
      enddo
    enddo
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 2,LA_Global
        K = KSZ_Global(L)
        BCQSUM(L,K) = BCQSUM(L,K)/SUMTIME
      enddo
    endif
  else
    ! *** SINGLE LAYER
    do L = 2,LA_Global
      BCQSUME(L) = BCQSUME(L)/SUMTIME
    enddo
  endif

  ! ***  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  do LL = 1,NPBS
    LQ = LPBS(LL)
    LN = LNC_Global(LQ)
    do K = 1,KC
      BCVHDX2(LN,K) = BCVHDX2(LN,K)/SUMTIME
    enddo
  enddo
  do LL = 1,NPBW
    LQ = LPBW(LL)
    LE = LEC_Global(LQ)
    do K = 1,KC
      BCUHDY2(LE,K) = BCUHDY2(LE,K)/SUMTIME
    enddo
  enddo
  do LL = 1,NPBE
    LQ = LPBE(LL)
    do K = 1,KC
      BCUHDY2(LQ,K) = BCUHDY2(LQ,K)/SUMTIME
    enddo
  enddo
  do LL = 1,NPBN
    LQ = LPBN(LL)
    do K = 1,KC
      BCVHDX2(LQ,K) = BCVHDX2(LQ,K)/SUMTIME
    enddo
  enddo

  return

  END SUBROUTINE

  SUBROUTINE BCOUT_INITIALIZE

  ! *** SETUP ACCUMULATION ARRAYS - @todo specify what the hell each is accumulating so a human can read...

  allocate(BCQSUM(LCM_Global,  KCM))  ! *** Modified so LCM is global
  allocate(BCUHDY2(LCM_Global, KCM))  ! *** Modified so LCM is global
  allocate(BCVHDX2(LCM_Global, KCM))  ! *** Modified so LCM is global
  allocate(BCQSUME(LCM_Global))       ! *** Modified so LCM is global

  SUMTIME = 0.
  BCQSUM = 0.
  BCUHDY2 = 0.
  BCVHDX2 = 0.
  BCQSUME = 0.

  INITBCOUT = 1

  END SUBROUTINE

  SUBROUTINE BCOUT_ZERO
  ! *** COMPUTE AVEARGE BOUNDARY INFLOWS/OUTFLOWS

  integer(IK4) :: L,K,LL,LQ,LE,LN,IBC

  SUMTIME = 0.

  ! *** SUM SELECTIVE QSUM
  if( KC > 1 )then
    do L = 2,LA_Global
      BCQSUM(L,KC) = 0.
    enddo
    do IBC = 1,NBCS
      L = LBCS(IBC)
      do K = 1,KS
        BCQSUM(L,K) = 0.
      enddo
    enddo
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      do L = 2,LA_Global
        K = KSZ_Global(L)
        BCQSUM(L,K) = 0.
      enddo
    endif
  else
    ! *** SINGLE LAYER
    BCQSUME = 0.
  endif

  ! ***  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  do LL = 1,NPBS
    LQ = LPBS(LL)
    LN = LNC_Global(LQ)
    do K = 1,KC
      BCVHDX2(LN,K) = 0.
    enddo
  enddo
  do LL = 1,NPBW
    LQ = LPBW(LL)
    LE = LEC_Global(LQ)
    do K = 1,KC
      BCUHDY2(LE,K) = 0.
    enddo
  enddo
  do LL = 1,NPBE
    LQ = LPBE(LL)
    do K = 1,KC
      BCUHDY2(LQ,K) = 0.
    enddo
  enddo
  do LL = 1,NPBN
    LQ = LPBN(LL)
    do K = 1,KC
      BCVHDX2(LQ,K) = 0.
    enddo
  enddo

END SUBROUTINE

INTEGER FUNCTION BLOCKBC(CELL3D)
  integer :: DSIZE,CELL2D,CELL3D
  CELL2D = LA_Global - 1
  DSIZE = 2                 ! *** EETIME
  if( KC > 1 )then
    DSIZE = DSIZE + CELL2D
    DSIZE = DSIZE + KC*NBCS
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      DSIZE = DSIZE + CELL2D
    endif
  else
    DSIZE = DSIZE + CELL2D  ! *** SINGLE LAYER
  endif

  ! ***  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  DSIZE = DSIZE + KC*NPBS
  DSIZE = DSIZE + KC*NPBW
  DSIZE = DSIZE + KC*NPBE
  DSIZE = DSIZE + KC*NPBN
  if( NQCTL > 0 )then
    DSIZE = DSIZE + NQCTL*2*KC
    if(NQCTLSER > 0 .or. NQCTRULES > 0 )then
      DSIZE = DSIZE + 7*NQCTL
    endif
  endif
  if( NQWR > 0 )then
    DSIZE = DSIZE + 2*NQWR
  endif
  BLOCKBC = DSIZE*4
END FUNCTION

SUBROUTINE SHELLFISHOUT()

  ! *** SHELLFISH OUTPUT
  integer(IK4) :: VER,HSIZE,BSIZE
  integer(IK4) :: I,L,K,ITMP,NN,NS,ISTAT
  integer(IK8) :: FSIZE, OFFSET
  real(RK4)    :: TMP
  real(RKD)    :: PTIME

  FILENAME = OUTDIR//'EE_SHF.OUT'
  EE_UNIT = 103                      ! *** EE_SHF.OUT
  if( ISRESTI /= 0 .and. ICONTINUE == 1 )then
    write(*,'(A)')'READING TO STARTING TIME FOR SHELLFISH'
    FSIZE = FILESIZE(FILENAME)
    open(EE_UNIT,FILE = FILENAME,ACTION = 'READWRITE',STATUS = 'OLD',FORM = FMT_BINARY,SHARED)     
    read(EE_UNIT) VER,HSIZE,BSIZE
    if( VER /= 8400 )then
      write(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'      ! *** Run continuation can only work with the latest version
      call STOPP('.')
    endif
    OFFSET = HSIZE
    
    NS = 0
    do while(OFFSET < FSIZE)
#ifdef GNU
      call FSEEK(EE_UNIT,OFFSET,0, ISTAT)

      if(IS_IOSTAT_END(ISTAT)) then
        write(EE_UNIT) EETIME
        exit
      endif
#else
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      if( EOF(EE_UNIT) )then
        write(EE_UNIT) EETIME
        exit
      endif
#endif  
      if( ISTAT /= 0 ) exit
        
      read(EE_UNIT) PTIME
      
      if( DEBUG )then
        NS = NS+1
        if( NS == 8 )then
          write(*,'(F10.3)')PTIME
          NS = 0
        else
          write(*,'(F10.3)',ADVANCE='NO') PTIME
        endif
      endif
      if( ABS(PTIME-TIMEDAY) <= 1E-4 .or. PTIME > TIMEDAY ) exit

      OFFSET = OFFSET + BSIZE
    enddo
    if( DEBUG ) WRITE(*,'(" ")')
    write(*,'(A)')'FINISHED READING SHELLFISH'
  else
    NERR = 0
    IORIGIN = 8
    100 open(EE_UNIT,FILE = FILENAME,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)
    write(EE_UNIT) EETIME
  endif
  
  if( ISRESTI == 0 .or. ICONTINUE == 0) NRESTART = 0 
  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).CELL_C(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI_C(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).SHLEN(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI_D(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI_H(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
    do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).HARV_C(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).NP(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
    do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).PR(L,K),4), K = KSZ(L),KC)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).CRB(L,K),4), K = KSZ(L),KC)
    enddo
  enddo
    do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).SPAWN(L,K),4), K = KSZ(L),KC)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).WQAQ(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).WQEA(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).FR(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).BMG(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_RESPI(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_GRAZI(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_DEATH(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_FECAL(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_URINE(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_RPOC(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_LPOC(L),4)
    enddo
  enddo  
  do I = 1,NSF
    do L = 2,LA
      if(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_DOC(L),4)
    enddo
  enddo

  FLUSH(EE_UNIT)
  close(EE_UNIT,STATUS = 'KEEP') 
  
  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_SHF.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif
  
END SUBROUTINE

SUBROUTINE GRIDOUT
  use INFOMOD,only:READSTR
  implicit none
  integer(IK4) :: NACTIVE,VER,HSIZE,BSIZE,K,L,CELL3D

  NACTIVE = LA_Global-1

  FILENAME = OUTDIR//'EE_GRD.OUT'
  VER   = 8400
  HSIZE = 8*4
  BSIZE = 0
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN')
  close(EE_UNIT,STATUS = 'DELETE')
  open(EE_UNIT,FILE = FILENAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY,SHARED)
  write(EE_UNIT) VER,HSIZE,BSIZE
  write(EE_UNIT) INT(IC,4),INT(JC,4),INT(KC,4),NACTIVE,IGRIDV

  !XCR(1:4,JMN:JMX,IMN:IMX),YCR(1:4,JMN:JMX,IMN:IMX)
  if( IGRIDV > 0 )then ! Sigma-Zed level
    write(EE_UNIT) (INT(KSZ_Global(L),4), L = 2,LA_Global)
    write(EE_UNIT) ((REAL(DZC(L,K),4), K = 1,KC), L = 2,LA_Global)
  else  ! *** Standard sigma stretched level
    write(EE_UNIT) (REAL(DZCK(K),4), K = 1,KC)
  endif
  close(EE_UNIT,STATUS = 'KEEP')

END SUBROUTINE

SUBROUTINE ARRAYSOUT

  ! *** TIME VARIABLE USER SPECIFIED ARRAYS IN THIS SECTION

  ! FLAGS: ARRAY TYPE, TIME VARIABLE
  ! ARRAY TYPE:    0 = L            DIM'D
  !                1 = L,KC         DIM'D
  !                2 = L,0:KC       DIM'D
  !                3 = L,KB         DIM'D
  !                4 = L,KC,NCLASS  DIM'D
  ! TIME VARIABLE: 0 = NOT CHANGING
  !                1 = TIME VARYING

  integer(IK4) :: K,L,ITYPE,ITIMEVAR,NAL
  real(RK4)    :: ZERO
  real(RK4)    :: TMPVAL
  character*8  ARRAYNAME
  character*2  SNUM

  if( ISINWV == 2 )then
    ZERO = 0.0

    FILENAME = OUTDIR//'EE_ARRAYS.OUT'
    EE_UNIT  = 105                      ! *** EE_ARRAYS.OUT

    NERR = 0
    IORIGIN = 10
    100 open(EE_UNIT,FILE = FILENAME,POSITION = 'APPEND',STATUS = 'UNKNOWN',FORM = FMT_BINARY,SHARED,ERR = 999,IOSTAT = IERRio)

    ! *** Time varibale flag
    ITIMEVAR = 1
    
    ! *** Hydrodynamic varibales
    if( IS_ARRAY_OUT(2) == 1 )then
      ! *** Water pressure
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'P'
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(P(L),4)
      enddo
      
      ! *** Average over depth flow
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'UHDYE' 
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(UHDYE(L),4)
      enddo
          
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'VHDXE' 
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(VHDXE(L),4)
      enddo
      
      ! *** Internal shear
      ITYPE = 1
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'DU'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(DU(L,K),4)
        enddo
      enddo
    
      ITYPE = 1
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'DV'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(DV(L,K),4)
        enddo
      enddo
      
      ! *** Bed shear
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'TBX' 
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(TBX(L),4)
      enddo
      
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'TBY' 
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(TBY(L),4)
      enddo
      
      ! *** Surface shear
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'TSX' 
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(TSX(L),4)
      enddo
      
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'TSY' 
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(TSY(L),4)
      enddo
    endif
    
    ! *** Turbulence quantities
    if( IS_ARRAY_OUT(3) == 1 )then
      ! *** Turbulence intensity or kinetic energy (m2/s2)
      ITYPE = 2
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'QQ'
      write(EE_UNIT) ARRAYNAME
      do K = 0,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(QQ(L,K),4)
        enddo
      enddo
      ! *** Turbulence length-scale (m)
      ITYPE = 2
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'DML'
      write(EE_UNIT) ARRAYNAME
      do K = 0,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(DML(L,K),4)
        enddo
      enddo    
      ! *** Vertical eddy visocity (m/s)
      ITYPE = 1
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'AV'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(AV(L,K),4)
        enddo
      enddo
      ! *** Vertical eddy diffusivity (m/s)
      ITYPE = 1
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'AB'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          write(EE_UNIT) REAL(AB(L,K),4)
        enddo
      enddo
    endif
        
    ! *** Thermal quantities
    if( IS_ARRAY_OUT(4) == 1 )then
      ! *** Solar radiation at top layers    
      ITYPE = 1
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'RADTOP'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          TMPVAL = RADTOP(L,K)
          write(EE_UNIT) TMPVAL
        enddo
      enddo

      ! *** Solar radiation at bottom layers 
      ITYPE = 1
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'RADBOT'
      write(EE_UNIT) ARRAYNAME
      do K = 1,KC
        do L = 2,LA_Global
          TMPVAL = RADBOT(L,K)
          write(EE_UNIT) TMPVAL
        enddo
      enddo
      if( ISTOPT(2) == 1 .or. ISTOPT(2) == 2) then
        ! *** Long wave solar radiation heat flux        
        ITYPE = 0
        write(EE_UNIT) ITYPE, ITIMEVAR
        ARRAYNAME = 'HBLW' 
        write(EE_UNIT) ARRAYNAME
        do L = 2,LA_Global
          write(EE_UNIT) REAL(HW_OUT(L),4)
        enddo
        ! *** Latent heat flux        
        ITYPE = 0
        write(EE_UNIT) ITYPE, ITIMEVAR
        ARRAYNAME = 'HBCV' 
        write(EE_UNIT) ARRAYNAME
        do L = 2,LA_Global
          write(EE_UNIT) REAL(HL_OUT(L),4)
        enddo
        ! *** Sensible heat flux        
        ITYPE = 0
        write(EE_UNIT) ITYPE, ITIMEVAR
        ARRAYNAME = 'HBEV' 
        write(EE_UNIT) ARRAYNAME
        do L = 2,LA_Global
          write(EE_UNIT) REAL(HS_OUT(L),4)
        enddo
      endif
    endif

    ! *** Water Quality quantities
    if( IS_ARRAY_OUT(5) == 1 )then
      ! *** Reareation rate (1/day)
      ITYPE = 0
      write(EE_UNIT) ITYPE, ITIMEVAR
      ARRAYNAME = 'WQREA'
      write(EE_UNIT) ARRAYNAME
      do L = 2,LA_Global
        write(EE_UNIT) REAL(WQRREA(L),4)
      enddo
      
      !*** Algae Kinetic rate
      ITYPE = 0
      do NAL = 1, NALGAE
        write(SNUM,'(I2.2)') NAL
        ! *** Growth rate
        write(EE_UNIT) ITYPE, ITIMEVAR
        ARRAYNAME = 'WQPA_'//SNUM
        write(EE_UNIT) ARRAYNAME
        do L = 2,LA_Global
          write(EE_UNIT) REAL(WQPA(L,NAL),4)
        enddo
        ! *** Metabolism rate
        write(EE_UNIT) ITYPE, ITIMEVAR
        ARRAYNAME = 'WQBM_'//SNUM
        write(EE_UNIT) ARRAYNAME
        do L = 2,LA_Global
          write(EE_UNIT) REAL(WQBM(L,NAL),4)
        enddo
        ! *** Death/predation rate
        write(EE_UNIT) ITYPE, ITIMEVAR
        ARRAYNAME = 'WQPR_'//SNUM
        write(EE_UNIT) ARRAYNAME
        do L = 2,LA_Global
          write(EE_UNIT) REAL(WQPR(L,NAL),4)
        enddo
      enddo
    endif
    

    !IF( .FALSE. )then
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'P'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(P(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FUHDYE' 
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FUHDYE(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FVHDXE' 
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FVHDXE(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FXE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    TMPVAL = FXE(L)
    !    write(EE_UNIT) TMPVAL
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FYE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    TMPVAL = FYE(L)
    !    write(EE_UNIT) TMPVAL
    !  enddo
    !
    !ENDIF
    !
    !IF( (ISTRAN(1) + ISTRAN(2)) > 0 )then
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FPGXE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FPGXE(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FPGYE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FPGYE(L),4)
    !  enddo
    !ENDIF
    !
    !IF( ISCURVATURE )then
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FCAXE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FCAXE(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FCAYE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FCAYE(L),4)
    !  enddo
    !ENDIF
    !
    !IF( ISHDMF >= 1 )then
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FMDUX'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = FMDUX(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FMDUY'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = FMDUY(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FMDVX'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = FMDVX(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'FMDVY'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = FMDVY(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !ENDIF
    !
    !IF( ISTRAN(2) > 0 )then
    !  if( IASWRAD == 3 )then
    !    ITYPE = 1
    !    write(EE_UNIT) ITYPE, ITIMEVAR
    !    ARRAYNAME = 'RADKE'
    !    write(EE_UNIT) ARRAYNAME
    !    do K = 1,KC
    !      do L = 2,LA_Global
    !        TMPVAL = RADKE(L,K)
    !        write(EE_UNIT) TMPVAL
    !      enddo
    !    enddo
    !  endif
    !  
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'RADTOP'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = RADTOP(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !  
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'RADNET'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = RADNET(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !  
    !  ITYPE = 1
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'RADBOT'
    !  write(EE_UNIT) ARRAYNAME
    !  do K = 1,KC
    !    do L = 2,LA_Global
    !      TMPVAL = RADBOT(L,K)
    !      write(EE_UNIT) TMPVAL
    !    enddo
    !  enddo
    !  
    !ENDIF
    !
    !IF( ISTRAN(8) > 0 )then
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'WQREA(1/day)'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(WQRREA(L),4)
    !  enddo
    !ENDIF
    !
    !! *** Miscellaneous but useful arrays
    !if( .false. )then
    !  ! *** SURFACE HEAT FLUX ARRAYS
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HBLW'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TEST1(1,L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HBCD'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TEST1(2,L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HBCV'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TEST1(3,L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HBEV'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TEST1(4,L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HBNT'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TEST1(5,L),4)
    !  enddo
    !  ! *** END OF SURFACE HEAT FLUX ARRAYS
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'WINDST'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(WINDST(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'WINDCD10'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(WINDCD10(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'TSX'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TSX(L),4)
    !  enddo
    !
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'TSY'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(TSY(L),4)
    !  enddo
    !endif
    !
    !IF( ISTRAN(2) > 0 .and. ISTOPT(2) == 2 )then
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'TS_COARE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(CDCOARE(L),4)
    !  enddo
    !  
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HS_COARE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(HSCOARE(L),4)
    !  enddo
    !  
    !  ITYPE = 0
    !  write(EE_UNIT) ITYPE, ITIMEVAR
    !  ARRAYNAME = 'HL_COARE'
    !  write(EE_UNIT) ARRAYNAME
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(HLCOARE(L),4)
    !  enddo
    !ENDIF
    
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'U'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(U(L,K),4)
    !  enddo
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'V'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(V(L,K),4)
    !  enddo
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'TBX'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TBX(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'TBY'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TBY(L),4)
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'B'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(B(L,K),4)
    !  enddo
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'U1'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(U1(L,K),4)
    !  enddo
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'V1'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(V1(L,K),4)
    !  enddo
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'V1U'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(V1U(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'U1V'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(U1V(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'H1U'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(H1U(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'H1V'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(H1V(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'STBX'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(STBX(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'STBY'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(STBY(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'LMASKDRY'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  if( LMASKDRY(L) )then
    !    write(EE_UNIT) REAL(1.0,4)
    !  else
    !    write(EE_UNIT) REAL(0.0,4)
    !  endif
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'TBX'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TBX(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'TBY'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TBY(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'QSBDLDX'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(QSBDLDX(L,2),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'QSBDLDY'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(QSBDLDY(L,2),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'BLDELTA'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL( ( QSBDLDX(L,2)-QSBDLDX(LEC_Global(L),2) + QSBDLDY(L,2)-QSBDLDY(LNC(L),2) ),4)
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FXMHK'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    if( K < KSZ_Global(L) )then
    !      write(EE_UNIT) REAL(-999,4)
    !    else
    !      write(EE_UNIT) REAL(FXMHK(L,k),4)
    !    endif
    !  enddo
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FYMHK'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    if( K < KSZ_Global(L) )then
    !      write(EE_UNIT) REAL(-999.,4)
    !    else
    !      write(EE_UNIT) REAL(FYMHK(L,K),4)
    !    endif
    !  enddo
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FXMHKE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(FXMHKE(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FYMHKE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(FYMHKE(L),4)
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FXSUP'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    if( K < KSZ_Global(L) )then
    !      write(EE_UNIT) REAL(-999,4)
    !    else
    !      write(EE_UNIT) REAL(FXSUP(L,k),4)
    !    endif
    !  enddo
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FYSUP'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    if( K < KSZ_Global(L) )then
    !      write(EE_UNIT) REAL(-999.,4)
    !    else
    !      write(EE_UNIT) REAL(FYSUP(L,K),4)
    !    endif
    !  enddo
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FXSUPE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(FXSUPE(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FYSUPE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(FYSUPE(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'ET'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TEST1(1,L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'CSHE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TEST1(2,L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'ICETEMP'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(TEST1(3,L),4)
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FRAZILICE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FRAZILICE(L,K),4)
    !  enddo
    !ENDDO

    !ITYPE = 2
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'W2'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 0,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(W2(L,K),4)
    !  enddo
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'UV'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(UV(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'VU'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(VU(L),4)
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FXVEG'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FXVEG(L,K),4)
    !  enddo
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FYVEG'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO K = 1,KC
    !  do L = 2,LA_Global
    !    write(EE_UNIT) REAL(FYVEG(L,K),4)
    !  enddo
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FXVEGE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(FXVEGE(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT) ITYPE, ITIMEVAR
    !ARRAYNAME = 'FYVEGE'
    !WRITE(EE_UNIT) ARRAYNAME
    !DO L = 2,LA_Global
    !  write(EE_UNIT) REAL(FYVEGE(L),4)
    !ENDDO

    ! FLUSH(EE_UNIT)
    close(EE_UNIT,STATUS = 'KEEP')

  endif

  return
  
999 if( IERRio == 30 .and. NERR < 3 )then
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_ARRAYS.OUT) error.  Try # ',nerr
      call SLEEP(5)
      if( NERR == 3 )then
        close(EE_UNIT,STATUS = 'KEEP')
      endif
      GOTO 100
    else
      call FILE_ERROR
      call STOPP('.')
    endif

END SUBROUTINE

INTEGER(IK8) FUNCTION FILESIZE(FNAME)
#ifndef GNU  
  USE IFPORT
#endif
  character*(*), intent(IN) :: FNAME
  integer(IK4) :: ISTAT
  integer(IK4) :: FINFO(12)

  ISTAT = STAT(FNAME,FINFO)     !  ISTAT = FSTAT(IUNIT,FINFO)
  if( ISTAT == 0 )then
    FILESIZE = FINFO(8)
  else
    FILESIZE = 0
  endif

  END FUNCTION
  
  SUBROUTINE FILE_ERROR
    
    SELECT CASE ( IORIGIN )
    CASE (0)
      write(*,'(a)') ' Error tying to open EE_BC.OUT'

    CASE (1)
      write(*,'(a)') ' Error tying to open EE_WS.OUT'

    CASE (2)
      write(*,'(a)') ' Error tying to open EE_VEL.OUT'

    CASE (3)
      write(*,'(a)') ' Error tying to open EE_WC.OUT'

    CASE (4)
      write(*,'(a)') ' Error tying to open EE_BED.OUT'

    CASE (5)
      write(*,'(a)') ' Error tying to open EE_SEDZLJ.OUT'

    CASE (6)
      write(*,'(a)') ' Error tying to open EE_WQ.OUT'

    CASE (7)
      write(*,'(a)') ' Error tying to open EE_SD.OUT'

    CASE (8)
      write(*,'(a)') ' Error tying to open EE_SHF.OUT'

    CASE (9)
      write(*,'(a)') ' Error tying to open EE_RPEM.OUT'

    CASE (10)
      write(*,'(a)') ' Error tying to open EE_ARRAYS.OUT'

    CASE DEFAULT
      write(*,'(a)') ' Error tying to open EE linkage files'
    
    END SELECT  
    
  END SUBROUTINE FILE_ERROR

  END MODULE
