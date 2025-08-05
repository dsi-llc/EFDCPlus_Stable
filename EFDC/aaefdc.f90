! ----------------------------------------------------------------------!
!   EFDC+: Multifunctional surface water modeling system
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------!
! Copyright 2021-2025 DSI, LLC
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------!
!
!  RELEASE:         EFDCPlus_12.3
!                   Domain Decomposition with MPI
!                   Propeller Wash 
!                   New WQ kinetics with user defined algal groups and zooplankton
!                   Dynamic time stepping adjustment for 3TL solution
!                   General Ocean Turbulence Model (GOTM)
!                   SIGMA-Zed (SGZ) Vertical Layering
!
!  DATE:            2025-05-01
!  BY:              DSI, LLC
!                   EDMONDS, WASHINGTON  98020
!                   USA
!                   WEBSITE:  https://www.eemodelingsystem.com/
!  LEAD DEVELOPER:  PAUL M. CRAIG

PROGRAM EFDC
  !
  ! *** WELCOME TO THE ENVIRONMENTAL FLUID DYNAMICS COMPUTER CODE PLUS (EFDC+)
  ! **
  ! *** ORIGINALLY DEVELOPED BY JOHN M. HAMRICK WHILE AT
  ! ***                         VIRGINIA INSTITUTE OF MARINE SCIENCE
  ! ***                         SCHOOL OF MARINE SCIENCE, THE COLLEGE OF
  ! ***                         WILLIAM AND MARY, GLOUCESTER POINT, VA 23062
  ! **
  ! *** THIS VERSION OF EFDC (EFDC+) HAS BEEN SIGNIFICANTLY ENHANCED AND
  ! *** UPDATED AND IS NOW MAINTAINED BY DSI, LLC, EDMONDS, WA.
  ! **
  ! *** EFDC SOLVES THE 3D REYNOLDS AVERAGED NAVIER-STOKES
  ! *** EQUATIONS (WITH HYDROSTATIC AND BOUSINESSQ APPROXIMATIONS) AND
  ! *** TRANSPORT EQUATIONS FOR TURBULENT INTENSITY, TURBULENT
  ! *** INTENSITY X LENGTH SCALE, SALINITY (OR WATER VAPOR CONTENT),
  ! *** TEMPERATURE, AN INERT TRACER (CALLED DYE), A DYNAMICALLY ACTIVE
  ! *** SUSPENDED SETTLING PARTICLE FIELD (CALLED SEDIMENT).  A FREE
  ! *** SURFACE OR RIGID LID IS PRESENT ON THE VERTICAL BOUNDARY Z = 1
  ! *** IN THE SIGMA STRETCHED VERTICAL COORDINATE.  THE HORIZONTAL
  ! *** COORDINATE SYSTEM IS CURVILINEAR AND ORTHOGONAL.
  ! *** THE NUMERICAL SOLUTION SCHEME IS ON A SPATIALLY STAGGERED MAC
  ! *** OR C GRID.
  ! *** SPATIAL SOLUTION OF THE EXTERNAL MODE FOR THE FREE SURFACE
  ! *** ELEVATION OR KINEMATIC PRESSURE UNDER THE RIGID LID IS BY
  ! *** CONJUGATE GRADIENT SOLUTION OF A PSEUDO-HEMHOLTZ EQUATION.
  ! *** THE INTERNAL SOLUTION IS IMPLICIT FOR THE VERTICAL SHEAR OR
  ! *** VELOCITY STRUCTURE.
  ! *** A NUMBER OF OPTIONS ARE AVAILABLE FOR REPRESENTING THE ADVECTIVE
  ! *** TRANSPORT TERMS IN THE MOMENTUM AND SCALAR TRANSPORT EQUATIONS.
  ! **
  ! *** PRIMARY DOCUMENTATION INCLUDES:
  !     HAMRICK, J. M., 1992:  A THREE-DIMENSIONAL ENVIRONMENTAL
  !     FLUID DYNAMICS COMPUTER CODE: THEORETICAL AND COMPUTATIONAL
  !     ASPECTS. THE COLLEGE OF WILLIAM AND MARY, VIRGINIA INSTITUTE
  !     OF MARINE SCIENCE, SPECIAL REPORT 317, 63 PP.
  !
  !     A COMPREHENSIVE EFDC+ THEORY KNOWLEDGE BASE HAS BEEN DEVELOPED BY
  !     DSI AND CAN BE FOUND AT:
  !     https://eemodelingsystem.atlassian.net/wiki/spaces/ETG/overview
  !
  ! *** DSI ASSUMES NO LIABILITY FOR use OF THIS CODE FOR ENVIRONMENTAL
  ! *** AND ENGINEERING STUDIES.
  ! **

  !----------------------------------------------------------------------!
  ! EFDC+ CHANGE RECORD
  ! DATE MODIFIED    BY                DESCRIPTION
  !----------------------------------------------------------------------!
  !    2002-02       John H. Hamrick   Added original sediment and toxics transport
  !    2004-06       Dudley Benton     Added dynamic memory allocation
  !                  Paul M. Craig
  !    2006-11       Paul M. Craig     Added Rooted Plant and Epiphyte Module
  !    2008-1O       Craig Jones       Added SEDZLJ and MHK working with SNL
  !                  Scott James
  !                  Paul M. Craig
  !    2009-03       Dang H. Chung     Added wind waves with sediment interaction with wave generated currents
  !                  Paul M. Craig
  !    2009-06       Paul M. Craig     Added new Lagrangian Particle Tracking Module
  !                  Dang H. Chung
  !    2010-01       Paul M. Craig     Added concentration/flow based BC's for the ICM WQ module
  !    2014-03       Paul M. Craig     Added heat coupled ice transport
  !    2011-06       Paul M. Craig     Added OpenMP implementation for most modules
  !                  Dang H. Chung
  !    2015-1O       Paul M. Craig     Implemented Sigma-Z
  !    2018-02       Paul M. Craig     Added unlimited dye classes with 0th and 1st order decay with temperature dependency
  !    2018-09       Paul M. Craig     Added chemical (toxics) volatilization and degradation
  !    2018-1O       Zander Mausolff   Implemented domain decomposition using OpenMPI
  !                  Paul M. Craig
  !    2020-06       Nghiem T. Lam     Added shellfish farm module
  !    2020-12       Paul M. Craig     Added propeller wash with linkage to the SEDZLJ and toxics module
  !                  Zander Mausolff
  !    2021-08       Duc Kien Tran     Rewrote WQ kinetics to allow unlimited phytoplankton and macrophyte classes, with unlimited zooplankton classes
  !                  Paul M. Craig
  !    2021-12       Paul M. Craig     Added propeller jet efflux momentum to EFDC+ flow field
  !    2022-01       Paul M. Craig     Added cohesive mass erosion classes for propwash induced resuspension

  use GLOBAL
  use OMP_LIB
  use Allocate_Initialize      
#ifndef GNU  
  USE IFPORT
#endif
  use XYIJCONV,only:XY2IJ
  use RESTART_MODULE
  use EFDCOUT
  use INFOMOD, only:SKIPCOM,READSTR
  use DRIFTER, only:DRIFTER_INP,AREA_CENTRD
  use FIELDS
  use WATERQUALITY, only:WQ3DINP
  use SHELLFISHMOD, only:ISFFARM,FREE_SHELLFISH
  use Cyclone
  use CONVERTWGS84
  use DRIFTER, only:close_nc_lpt

  ! *** MPI  modules
  use MPI
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out
  use MPI_All_Reduce
  use Broadcast_Routines
  use Communicate_Ghost_Routines
  use Mod_Map_Write_EE_Binary
  use Mod_Map_Gather_Sort

  ! *** Propwash
  use Variables_Propwash
  use Mod_Read_Propwash 
  use Mod_Setup_Ships

  implicit none

  real(RKD), allocatable,dimension(:,:) :: SHOTS
  integer,   allocatable,dimension(:,:) :: IFDCH
  integer,   allocatable,dimension(:)   :: IDX

  real, allocatable,dimension(:) :: THICK

  real :: TSHIFT, DET, XLNUTME, YLTUTMN, TA1, TA2, FORCSUM, ZERO, MAXTHICK, DEPTHMIN, MAXTHICK_local
  real :: CCUE, CCVE, CCUN, CCVN, TMPVAL, TMPCOR, ANG1, ANG2, ANG, DETTMP, DZPC
  real :: ANGTMP1, ANGTMP2, DYUP1, DYUM1, DDXDDDY, DDYDDDX, DXVLN, DXVLS, C
  real :: TMP, C1, DZGTMP, FRACK
  real :: GRADW, GRADE, GRADS, GRADN, MHKE, MHKE_Global, MHKL, MHKL_Global
  real(8) :: XUTM(1), YUTM(1), XLL(1), YLL(1)
  real :: ACS2
  real :: AS2
  real :: AC2
  real :: TNT            !< 
  real :: TPN

  real(RK4) :: CPUTIME(2)
  real(RKD) :: T0, T1, T2, T3, T4, T9, DELSNAP, TCPU, TTDS, TWAIT
  real(RKD), external :: DSTIME

  integer :: COUNT, NSNAPMAX, NT, iStatus, NS, NTMP, NSHOTS, ISNAP, NFDCHIJ
  integer :: ITMPVAL, IFIRST, ILAST, ISO, JDUMY, IISTMP, LPBTMP, LBELMIN
  integer :: L, K, KK, I, J, IP, IS, M, LL, LP, LF, LT, LN, LS, NX, LW, KM, LE, LG, MW, ND, IYEAR, NP, NC, MD
  integer :: IU, JU, KU, LU, ID, JD, KD, LD, NWR, NJP
  
  integer(IK4) :: IERROR
  integer(IK4) :: IRET
  integer(IK8) :: NREST,NN
  integer(1)   :: VERSION

  logical(4) :: RES, LDEBUG, BFLAG

  character*24  :: mpi_filename
  character*19  :: MPI_LIST
  character*20  :: BUFFER
  character*120 :: LINE
  character*8   :: DAY
  character*200 :: STR
  
  ! *** MPI Variables
  integer(4) :: IERR     !< local MPI error flag
  integer    :: IIN, JIN
  Double Precision :: starting_time, ending_time
  
  ! *** When updating the version date also update the EXE date in the project settings
  EFDC_VER = '2025-05-01'
  
  ! *** GET START TIME OF ELAPSED TIME COUNTER
  TCPU = DTIME(CPUTIME)
  TIME_START = DSTIME(1)
  VERSION = 3

  IERR = 0
#ifdef DEBUGGING
  call sleep(10)
#endif

  ! ****************************************************************************
  ! *** Get EFDC+ version
  NTHREADS = 1

#ifdef GNU
  COUNT = IARGC()
  call GETARG(0, STR)
#else
  COUNT = NARGS()
  call GETARG(0, STR, iStatus)
#endif 

  I = INDEX(STR,"\",.TRUE.)
  if( I == 0 )then
    I = INDEX(STR,":")
  endif
  EFDC_EXE = STR(I+1:LEN(STR))

  process_id = 0
  master_id  = 0
  num_Processors = 1

  ! ****************************************************************************
  ! *** Initialize MPI communicator
  call Initialize_MPI

  ! ***  Write to screen how many processes/domain are being used for each process
  write(*,11) num_Processors, process_id
11 FORMAT( '***********************************************************************************',/, &
           '***  This EFDCPlus run is using:',I6,' MPI domain(s).  Preparing domain:',i6,'  ***',/, &
           '***********************************************************************************',/)

  if( process_id == master_id )then

    ! *** GET THE COMMAND LINE ARGUMENTS, IF ANY
    if( COUNT > 1 )then
      ! *** ARGUMENT 1
#ifdef GNU
      call GETARG(1, BUFFER)
#else
      call GETARG(1, BUFFER, iStatus)
#endif          
      !$  if( Buffer(1:3) == '-NT' .or. Buffer(1:3) == '-nt' )then
      !$    read(Buffer(4:10),*) NTHREADS
      !$  endif
      if( COUNT > 2 )then
        ! *** ARGUMENT 2
#ifdef GNU
        call GETARG(2, BUFFER)
#else
        call GETARG(2, BUFFER, iStatus)
#endif            
        !$  if( Buffer(1:3) == '-NT' .or. Buffer(1:3) == '-nt' )then
        !$    read(Buffer(4:10),*) NTHREADS
        !$  endif
      endif
    endif

    if( VERSION <= 43 .and. VERSION >= 41 .and. NTHREADS > 2 ) NTHREADS = 2

    call WELCOME

  endif !*** End calculation on master process

  !print *, 'here'
  !pause
  call INITFIELDS()
  
  ! *** Send off the number of threads  read in from the command line so 
  ! that the OMP_SET_NUM_THREADS is set properly on all processes
  call MPI_Bcast(NTHREADS, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)

  ! *** BEGIN OMP SECTION
  !$OMP PARALLEL
  !$OMP MASTER
  !$    NT = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
  !$    if( NT < NTHREADS ) NTHREADS = NT
  !$    call OMP_SET_NUM_THREADS(NTHREADS)
  
  If( process_id == master_id )then
    !$    write(*,1) num_Processors, NTHREADS

1   FORMAT( '***********************************************************************************',/, &
            '***  This EFDC+ run will use:',I6,' MPI domain(s) and ',I2,' thread(s) per domain   ***',/, &
            '***********************************************************************************',/)
  endif

#ifdef _WIN
  OUTDIR = '#output\'
  MPI_Outdir = 'mpi_logs\'
#else
  OUTDIR = '#output/'
  MPI_Outdir = 'mpi_logs/'
#endif
#ifdef GNU  
  FMT_BINARY = 'UNFORMATTED'
#else
  FMT_BINARY = 'BINARY'
#endif 

  ! *** Always create the new file
  call MPI_barrier(MPI_Comm_World, ierr)
  
  If( process_id == master_id )then !***Start calculation on master process
    write(STR, '("*** EFDC+ Filename: ",A,",  Version: ",A10," ***")') TRIM(EFDC_EXE), EFDC_VER

    ! *** OPEN OUTPUT FILES
    open(9,FILE = OUTDIR//'TIME.LOG', STATUS = 'REPLACE')
    write(9, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
    close(9)

    ! *** DELETE THE FOLLOWING
    open(1,FILE = OUTDIR//'DRYWET.LOG', STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'SEDIAG.OUT', STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'CFL.OUT',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'NEGSEDSND.OUT', STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'ERROR.LOG', STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')

    ! *** PAUSE THE SCREEN TO ALLOW THE USER TO REVIEW THE NUMBER OF THREADS
#ifdef GNU  
    call SLEEP(3)       !< 3 seconds
#else
    call SLEEPQQ(3000)  !< 3000 milliseconds
#endif    

  endif !***end calculation on master process

  ! ****************************************************************************
  call MPI_BCAST(NTHREADS, 1, MPI_Int, master_id, MPI_Comm_World,ierr)
  call MPI_BARRIER(MPI_Comm_World, ierr)
  ! ****************************************************************************

  ! *** CALL INPUT SUBROUTINE
  NRESTART = 0
  call VARINIT

  ! *** Allocate variables for domain decomposition
  call Allocate_Domain_Decomp

  call INPUT
  
  ! *** Writes out the variables read in by the input on each process
  !Call VerifyInput      ! *** Used to test the input/broadcast process.  Activate for testing.

  If( process_id == master_id )then
    ! *** Handle Runtime Flag
    open(1,FILE = '0run',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = '0run',STATUS = 'UNKNOWN')
    close(1)
  endif !***End calculation on master process

  ! *** SET TIME RELATED parameterS
  ! *** THE parameter NTC = NUMBER OF TIME CYCLES, CONTROLS THE LENGTH OF RUN (NUMBER OF TIME STEPS)
  ! *** STARTING AND ENDING DATES, IN DAYS
  TIMEDAY  = DBLE(TCON)*DBLE(TBEGIN)/86400._8
  TIMEEND  = TIMEDAY + (DBLE(TIDALP)*DBLE(NTC) + DBLE(DT))/86400._8
  HOURNEXT   = TIMEDAY + 1./24.
  HOUR06NEXT = TIMEDAY + 6./24.
  HOUR12NEXT = TIMEDAY + 12./24.

  ! *** MODEL RUN TIME ACCUMULATORS
  !TCYCLE = 0.0
  TLRPD = 0.0
  THDMT = 0.0
  TVDIF = 0.0
  TCGRS = 0.0
  TTSED = 0.0
  TSSTX = 0.0
  TCONG = 0.0
  TSADV = 0.0
  TPUV = 0.0
  TCEXP = 0.0
  TAVB = 0.0
  THEAT = 0.0
  THMDF = 0.0
  TMISC = 0.0
  TUVW = 0.0
  TQCTL = 0.0
  TQQQ = 0.0
  TMPIEE = 0.0
  TMPIGH = 0.0
  TMPITMP = 0.0
  TWQKIN = 0.0
  TWQRPEM = 0.0
  TWQSED = 0.0
  DSITIMING = 0.0
  TTWAIT = 0.0
  TPROPW = 0.0
  
  CFMAX = CF
  PI = ACOS(-1.0)

  ! ***  NTC:     Number of reference time periods in run
  ! ***  NTSPTC:  Number of time steps per reference time period
  ! ***  TCON:    Conversion multiplier to change tbegin to seconds
  ! ***  TBEGIN:  Time origin of run
  ! ***  TIDALP:  Reference time period in sec (ie 44714.16s or 86400s)
  TPN = REAL(NTSPTC)
  NTS = INT8(NTC)*NTSPTC/NFLTMT    ! *** Total # of time steps for entire simulation (typically NFLTMT = 1)
  NLTS = NTSPTC*NLTC               ! *** # Transition Step to Completely linear
  NTTS = NTSPTC*NTTC               ! *** # Transition Step to Fully Non-linear
  NTTS = NTTS+NLTS                 ! *** Total # of Steps to End of Transition
  SNLT = 0.
  NCTBC = 1
  NTSVB = NTCVB*NTSPTC             ! *** Variable Bouyancy
  NTSPTC2 = 2*NTSPTC/NFLTMT        ! *** Twice # of time steps per reference time period
  NSHOWR = 0
  NSHOWC = 0
  do NS = 1,NASER
    MTSALAST(NS) = 2
  enddo
  do NS = 1,NWSER
    MTSWLAST(NS) = 2
  enddo
  do NS = 1,NPSER
    MTSPLAST(NS) = 2
  enddo
  do NS = 1,NQSER
    MTSQLAST(NS) = 2
  enddo
  do NS = 1,NQWRSR
    MTSWRLAST(NS) = 2
  enddo
  do NS = 1,NGWSER
    MTSGWLAST(NS) = 2
  enddo
  do NC = 1,8
    do NN = 1,NCSER(NC)
      MTSCLAST(NN,NC) = 2
    enddo
  enddo
  MSFTLST = 1

  ! *** EFDC_EXPLORER PRIMARY OUTPUT FREQUENCY DEFINITION TABLE
  ! *** REMOVED THE OBSOLETE SUBROUTINE SURFPLT
  NSNAPSHOTS = 0
  NSHOTS = 0

  if( HFREOUT == 1 )then
    if( process_id == master_id )then
      write(*,'(A)')'HIGH FREQUENCY TIME SERIES ARE BEING USED'
      write(*,'(A)')'  READING SUBSET.INP'
      open(1,FILE = 'subset.inp',action = 'read')
      call SKIPCOM(1,'*',2)
      read(1,*,IOSTAT = ISO) NSUBSET
      if( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
    endif
    call Broadcast_Scalar(NSUBSET, master_id)
    
    allocate(HFREGRP(NSUBSET))
    call AllocateDSI( IJHFRE,    NSUBSET, 0)
    call AllocateDSI( NPNT,      NSUBSET, 0)

    call AllocateDSI( HFREDAYEN, NSUBSET, 0.)
    call AllocateDSI( HFREDAYBG, NSUBSET, 0.)
    call AllocateDSI( HFREDUR,   NSUBSET, 0.)
    call AllocateDSI( HFREDAY,   NSUBSET, 0.)
    call AllocateDSI( HFREMIN,   NSUBSET, 0.)

    if( process_id == master_id )then
      do NS = 1,NSUBSET
        call SKIPCOM(1,'*',2)
        read(1,*,IOSTAT = ISO) IS, IJHFRE(IS), HFREDAYBG(IS), HFREDUR(IS), HFREMIN(IS), NPNT(IS)
        if( ISO > 0 )then
          call STOPP('SUBSET.INP: READING ERROR!')
        endif

        call AllocateDSI( HFREGRP(IS).ICEL, NPNT(IS), 0)
        call AllocateDSI( HFREGRP(IS).JCEL, NPNT(IS), 0)
        call AllocateDSI( HFREGRP(IS).XCEL, NPNT(IS), 0.)
        call AllocateDSI( HFREGRP(IS).YCEL, NPNT(IS), 0.)
        allocate(HFREGRP(IS).NAME(NPNT(IS)))

        do NP = 1,NPNT(IS)
          call SKIPCOM(1,'*',2)
          read(1,'(a)') STR
          STR = adjustl(trim(STR))
          if( IJHFRE(IS) == 1 )then
            call PARSESTRING(STR, BUFFER)
            read(BUFFER,*) HFREGRP(IS).ICEL(NP)
            call PARSESTRING(STR, BUFFER)
            read(BUFFER,*) HFREGRP(IS).JCEL(NP)
          else
            call PARSESTRING(STR, BUFFER)
            read(BUFFER,*) HFREGRP(IS).XCEL(NP)
            call PARSESTRING(STR, BUFFER)
            read(BUFFER,*) HFREGRP(IS).YCEL(NP)
          endif
          if( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
          STR = TRIM(STR)
          if( STR(1:1) == '!') STR = STR(2:20)
          HFREGRP(IS).NAME(NP) = TRIM(STR)          
        enddo
      enddo
      close(1)
    endif    ! *** end calculation on master process
    
    call Broadcast_Array(IJHFRE,    master_id)
    call Broadcast_Array(HFREDAYBG, master_id)
    call Broadcast_Array(HFREDUR,   master_id)
    call Broadcast_Array(HFREMIN,   master_id)
    call Broadcast_Array(NPNT,      master_id)

    if( process_id /= master_id )then
      do NS = 1,NSUBSET
        call AllocateDSI( HFREGRP(NS).ICEL, NPNT(NS),   0)
        call AllocateDSI( HFREGRP(NS).JCEL, NPNT(NS),   0)
        call AllocateDSI( HFREGRP(NS).XCEL, NPNT(NS), 0.0)
        call AllocateDSI( HFREGRP(NS).YCEL, NPNT(NS), 0.0)
      enddo
    endif
    
    do NS = 1,NSUBSET
      call Broadcast_Array(HFREGRP(NS).ICEL, master_id)
      call Broadcast_Array(HFREGRP(NS).JCEL, master_id)
      call Broadcast_Array(HFREGRP(NS).XCEL, master_id)
      call Broadcast_Array(HFREGRP(NS).YCEL, master_id)
    enddo
    
    ! *** UPDATE LOCATIONS TO LOCAL DOMAIN, REMOVING ANY POINTS NOT IN THE CURRENT PROCESS
    if( .not. allocated(XCOR) ) CALL AREA_CENTRD
    if( process_id == master_id )then
      do NS = 1,NSUBSET
        if( IJHFRE(NS) == 0 )then
          ! *** CONVERT IJ TO LOCAL
          call XY2IJ(HFREGRP(NS), 1)          ! *** Get I & J from X and Y
        endif
        HFREDAY(NS)   = HFREDAYBG(NS)
        HFREDAYEN(NS) = HFREDAYBG(NS) + HFREDUR(NS)/24.
      enddo 
    endif
  endif

  ! ***  TCON:    CONVERSION MULTIPLIER TO CHANGE TBEGIN TO SECONDS
  ! ***  TBEGIN:  TIME ORIGIN OF RUN
  ! ***  TIDALP:  REFERENCE TIME PERIOD IN SEC (IE 44714.16S OR 86400S)
  if( ISPPH == 0 )then
  
  elseif( ISPPH == 2 )then
    NSNAPSHOTS = 2
    SNAPSHOTS(1) = (TBEGIN*TCON + TIDALP*NTC)/86400. - 0.001
    SNAPSHOTS(2) = SNAPSHOTS(1) + 0.001
  else
    ! *** ISPPH == 1 OR 100

    DELSNAP = DBLE(TIDALP)/DBLE(NPPPH)/DBLE(86400.)
    T1 = DBLE(TBEGIN)*DBLE(TCON)/DBLE(86400.)
    T2 = (DBLE(TBEGIN)*DBLE(TCON) + DBLE(TIDALP)*DBLE(NTC))/DBLE(86400.)  ! *** TIMEDAY AT END OF RUN
    T0 = T1
    T9 = T2 + 1.E-6                                                     ! *** TIMEDAY AT END OF RUN PLUS PRECISION TOLERANCE
    NN = (T2-T1)/DELSNAP+2

    ! *** HIGH FREQUENCY SNAPSHOTS
    if( ISPPH == 100 )then
      write(*,'(A)')'HIGH FREQ SNAPSHOTS USED'
      NS = 0
      write(*,'(A)')'  READING SNAPSHOTS.INP'
      open(1,FILE = 'snapshots.inp',STATUS = 'UNKNOWN',ERR = 999)

      call SKIPCOM(1,'*')
      read(1,*)NSHOTS
      allocate(SHOTS(NSHOTS+1,3))
      SHOTS = 0.0
      do NS = 1,NSHOTS
        read(1,*)(SHOTS(NS,K),K = 1,3)
        SHOTS(NS,3) = SHOTS(NS,3)/1440._8
        SHOTS(NS,2) = SHOTS(NS,2)/24._8
      enddo
999   NSHOTS = NS-1
      if( NSHOTS < 0 )NSHOTS = 0
      close(1)

      SHOTS(NSHOTS+1,1) = T2

      do NS = 1,NSHOTS
        NN = NN + NINT(SHOTS(NS,2)/SHOTS(NS,3),8) + 1
      enddo
    endif

    ! *** Reallocate the SNAPSHOTS array
    NSNAPMAX = NN+2
100 DEALLOCATE(SNAPSHOTS)

    call AllocateDSI( SNAPSHOTS, NSNAPMAX, 0.0)

    ! *** BUILD THE SNAPSHOT DATES
    ISNAP = 1
    if( NSHOTS > 0 )then
      T4 = SHOTS(ISNAP,1)+SHOTS(ISNAP,2)  ! *** ENDING TIME FOR HIGH FREQ
      do while (T4 < T0)
        ISNAP = ISNAP+1
        if( ISNAP > NSHOTS )then
          NSHOTS = 0
          exit
        endif
        T4 = SHOTS(ISNAP,1) + SHOTS(ISNAP,2)
      enddo
    endif

    NSNAPSHOTS = 0
    T1 = T0
    do while ( T1 <= T9 )
      if( ISNAP <= NSHOTS )then
        if( T1 > SHOTS(ISNAP,1) )then
          ! *** ENSURE NO OVERLAPPING PERIODS
          T3 = SHOTS(ISNAP,1)
          if( NSNAPSHOTS >= 1 )then
            do while (T3 < SNAPSHOTS(NSNAPSHOTS))
              T3 = T3+SHOTS(ISNAP,3)
            enddo
          endif

          T4 = SHOTS(ISNAP,1)+SHOTS(ISNAP,2)  ! *** ENDING TIME FOR HI!GH FREQ
          if( T4 > T0 .and. T3 < T9 )then
            ! *** VALID PERIOD, SO BUILD HF SNAPSHOTS

            ! *** CHECK 1ST AND LAST HF PERIOD
            if( T4 > T9 )T4 = T2
            do while ( T3 <= T0 )
              T3 = T3+SHOTS(ISNAP,3)
            enddo

            ! *** BUILD HF SNAPSHOTS
            do while ( T3 <= T4 )
              NSNAPSHOTS = NSNAPSHOTS+1
              SNAPSHOTS(NSNAPSHOTS) = T3
              T3 = T3+SHOTS(ISNAP,3)
            enddo
          endif
          ISNAP = ISNAP+1  ! *** INCREMENT TO THE NEXT PERIOD

          ! *** Synch up the regular intervals
          if( NSNAPSHOTS > 0 )then
            do while (T1 < SNAPSHOTS(NSNAPSHOTS))
              T1 = T1+DELSNAP
            enddo
          endif

          if( T1 < SHOTS(ISNAP,1) .and. T1 >= T0 .and. T1 <= T9 )then
            NSNAPSHOTS = NSNAPSHOTS+1
            if( NSNAPSHOTS > NSNAPMAX )then
              NSNAPMAX = NSNAPMAX + 10
              GOTO 100
            endif
            SNAPSHOTS(NSNAPSHOTS) = T1
          endif
        elseif( T1 >= T0 .and. T1 <= T9 )then
          NSNAPSHOTS = NSNAPSHOTS+1
          if( NSNAPSHOTS > NSNAPMAX )then
            NSNAPMAX = NSNAPMAX + 10
            GOTO 100
          endif
          SNAPSHOTS(NSNAPSHOTS) = T1
        endif
      elseif( T1 >= T0 .and. T1 <= T9 )then
        NSNAPSHOTS = NSNAPSHOTS+1
        if( NSNAPSHOTS > NSNAPMAX )then
          NSNAPMAX = NSNAPMAX + 10
          GOTO 100
        endif
        SNAPSHOTS(NSNAPSHOTS) = T1
      endif
      T1 = T1+DELSNAP
    enddo
  endif

  if( ISPPH == 100 ) ISPPH = 1  !DHC: 2013-08-15

  if( process_id == master_id )then
    if( MPI_DEBUG_FLAG )then
      write(*,'(A,I8)')'NSNAPSHOTS = ',  NSNAPSHOTS
      write(mpi_efdc_out_unit,*)'NSNAPSHOTS = ',  NSNAPSHOTS
      do I = 1,NSNAPSHOTS
        write(mpi_efdc_out_unit,*)'SNAPSHOT: ',I,SNAPSHOTS(I)
      enddo
    endif
  endif !*** End calc on master process

  ! *** ENSURE LAST SNAPSHOT SPANS THE END OF THE MODEL RUN
  SNAPSHOTS(NSNAPSHOTS+1) = T2 + 1.
  NSNAPMAX = NSNAPSHOTS
  NSNAPSHOTS = min(NSNAPSHOTS,2)        ! *** NSNAPSHOTS = 1 IS THE INITIAL CONDITION.  SET TO FIRST MODEL RESULTS SNAPSHOT

  ! *** SET CONTROLS FOR WRITING TO FILTERED, AVERAGED OR RESIDUAL
  ! *** 2D SCALAR CONTOURING AND 2D VELOCITY VECTOR PLOTTING FILES
  ! *** RESIDUAL SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATION
  ! *** CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
  do N = 1,7
    if( ISRSPH(N) >= 1 ) JSRSPH(N) = 1
  enddo

  ! *** Set controls for writing to drifter output
  ! *** 
  if( ISRESTR == 1 ) JSRESTR = 1  ! *** Mean mass transport
  JSWASP = 0
  if( ISWASP >= 1 ) JSWASP = 1

  ! *** SET SOME CONSTANTS
  CTURB2 = CTURB**0.667
  CTURB3 = CTURB**0.333
  KS     = KC-1
  if( KS == 0 ) KS = 1
  DZI = REAL(KC)
  DT  = TIDALP*REAL(NFLTMT)/REAL(NTSPTC)
  DTI = 1./DT
  DT2 = 2.*DT
  DTMIN  = DT
  AVCON1 = 2.*(1.-AVCON)*DZI*AVO
  G    = 9.81
  GPO  = G*BSC
  GI   = 1./G
  PI   = 3.1415926535898
  PI2  = 2.*PI
  TIMERST = (TIDALP*ABS(ISRESTO))/86400.

  ! *** SET WEIGHTS FOR SALINITY AND TEMPERATURE BOUNDARY INTERPOLATION
  if( KC > 1 )then
    do K = 1,KC
      WTCI(K,1) = REAL(K-KC)/REAL(1-KC)
      WTCI(K,2) = REAL(K-1)/REAL(KC-1)
    enddo
  else
    WTCI(1,1) = 0.5
    WTCI(1,2) = 0.5
  endif

  ! *** INITIALIZE ARRAYS
  call AINIT

  ! *** READ IN XLON AND YLAT OR UTME AND UTMN OF CELL CENTERS OF
  ! *** CURVILINEAR PORTION OF THE  GRID
  ISHELTERVARY = 0

  ierr = 0
  ! ****************************************************************************
  ! *** MPI communication
  call MPI_barrier(DSIcomm, ierr)
  call communicate_ghost_cells(DXU, 'DXU')
  call communicate_ghost_cells(DYU, 'DYU')
  call communicate_ghost_cells(DXV, 'DXV')
  call communicate_ghost_cells(DYV, 'DYV')
  call communicate_ghost_cells(HMU, 'HMU')
  call communicate_ghost_cells(HMV, 'HMV')
  ! ****************************************************************************

  ! *** Read LXLY.INP and set curvature flag
  ISCURVATURE = .FALSE.
  call AllocateDSI( R2D_Global, LCM_Global, 8, 0.0)
    
  if( process_id == master_id )then
    write(*,'(A)')'READING LXLY.INP'
    open(1,FILE = 'lxly.inp',STATUS = 'UNKNOWN')
      
    call SKIPCOM(1,'*')

    call utmpars
      
    ! *** Domain Centroid
    Center_X = 0.
    Center_Y = 0.
      
    ! *** READ LXLY
    do LL = 2,LA_Global
      if( ISCORV == 1 )then
        ! *** READ LXLY WITH SPATIALLY VARYING CORIOLIS
        ! ***                              1        2     3     4     5     6       7      8
        ! ***              IIN, JIN, XLNUTME, YLTUTMN, CCUE, CCVE, CCUN, CCVN, TMPVAL, TMPCOR
        read(1,*,ERR = 3000) IIN, JIN, (R2D_Global(1,J),J = 1,8)
      else
        ! *** READ LXLY WITHOUT SPATIALLY VARYING CORIOLIS
        read(1,*,ERR = 3000) IIN, JIN, (R2D_Global(1,J),J = 1,7)
        R2D_Global(1,8) = CF
      endif
        
      LG = LIJ_GLOBAL(IIN,JIN)
      do J = 1,8
        R2D_Global(LG,J) = R2D_Global(1,J)
      enddo
        
      if( ISWGS84 == 1 )then
        ! *** 2022-08-10, NTL: Convert Lon/Lat to UTM coordinates
        XLL(1) = R2D_Global(1,1)
        YLL(1) = R2D_Global(1,2)
        call UTM_WGS84(XLL, YLL, XUTM, YUTM)
        R2D_Global(1,1) = XUTM(1)
        R2D_Global(1,2) = YUTM(1)
      endif
        
      Center_X = Center_X + R2D_Global(1,1)
      Center_Y = Center_Y + R2D_Global(1,2)
    enddo
    close(1)
      
    Center_X = Center_X / REAL(LA_Global-1,8)
    Center_Y = Center_Y / REAL(LA_Global-1,8)
    Center_X = INT(Center_X)
    Center_Y = INT(Center_Y)
            
    do LG = 2,LA_GLOBAL
      ANG1 = ATAN2(R2D_Global(LG,5), R2D_Global(LG,3))
      ANG2 = ATAN2(-R2D_Global(LG,4),R2D_Global(LG,6))
      ANG = 0.5*(ANG1 + ANG2)
      if( SIGN(1.,ANG1) /= SIGN(1.,ANG2) )then
        if( ABS(ANG1) > (1.57) .or. ABS(ANG2) > (1.57) )then
          ! *** HANDLE THE DISCONTINUITY AT THE 180 DEGREES ANGLE
          ANG = ANG + ACOS(-1.0)
        endif
      endif
          
      CUE_Global(LG) = COS(ANG)
      CVE_Global(LG) = -SIN(ANG)
      CUN_Global(LG) = SIN(ANG)
      CVN_Global(LG) = COS(ANG)
    enddo
  endif   ! *** End of master process
    
  call MPI_Barrier(DSIcomm, ierr)
  call Broadcast_Array(R2D_Global, master_id)
  call Broadcast_Scalar(Center_X,  master_id)
  call Broadcast_Scalar(Center_Y,  master_id)
  call Broadcast_Array(CUE_Global, master_id)  
  call Broadcast_Array(CVE_Global, master_id)  
  call Broadcast_Array(CUN_Global, master_id)  
  call Broadcast_Array(CVN_Global, master_id)  
    
  FORCSUM = 0.

  ! *** Map to Local Domain
  do LG = 2,LA_GLOBAL
    FORCSUM = FORCSUM + R2D_Global(LG,8)
      
    L = Map2Local(LG).LL
    if( L > 1 )then
      DLON(L) = R2D_Global(LG,1)      ! *** This is easting coordinate in meters
      DLAT(L) = R2D_Global(LG,2)      ! *** This is northing coordinate in meters
          
      CUE(L) = CUE_Global(LG)
      CVE(L) = CVE_Global(LG)
      CUN(L) = CUN_Global(LG)
      CVN(L) = CVN_Global(LG)
      FCORC(L) = R2D_Global(LG,8)
          
      WINDSTKA(L) = R2D_Global(LG,7)
      if( L > 2 .and. WINDSTKA(L-1) /= WINDSTKA(L) ) ISHELTERVARY = 1
          
      DETTMP = 1./( CUE(L)*CVN(L) - CUN(L)*CVE(L) )
      if( DETTMP == 0.0 )then
        write(6,6262)
        write(6,6263) IL(L), JL(L)
        call STOPP('.')
      endif
    endif
  enddo
  FCORC(1) = FCORC(2)
  FCORC(LC) = FCORC(LA)
  if( FORCSUM > 1.0E-6 ) ISCURVATURE = .TRUE.

  deallocate(R2D_Global)
    
6262 FORMAT('  SINGULAR INVERSE TRANSFORM FROM E,N TO CURV X,Y')
6263 FORMAT('  I,J  = ',2I10/)


  ! *** COMPUTE CELL AREAS AND CENTROIDS
  if( ISWAVE >= 3 .or. ISPD > 0 .or. (ISWAVE >= 1 .and. ISWAVE <= 2 .and. IFWAVE == 1 .and. SWANGRP == 0) )then
    ! *** COMPUTE CELL AREAS AND CENTROIDS, REQUIRES CORNERS.INP FILE
    if( .not. allocated(XCOR) ) CALL AREA_CENTRD
  endif

  ! *** READ IN SEDFLUME DATA
  if( LSEDZLJ )then
    call SEDIC
  endif
  
  ! *** QC option for moving QSER series with depths are small 
  if( HDRYMOVE > 0.0 .and. ISDRY == 0 ) HDRYMOVE = 0.0                ! *** Prevent QSER cells from moving when ISDRY = 0
  if( HDRYMOVE > 0.0 .and. HDRYMOVE <= HDRY ) HDRYMOVE = 1.5*HDRY     ! *** Don't allow HDRYMOVE to be <= HDRY   
  
  GOTO 3002
3000 call STOPP('READ ERROR FOR FILE LXLY.INP')
3002 continue

  ZERO = 0.
  if( DEBUG )then
    write(mpi_log_unit,'(a)'), ' *** Writing LIJMAP.OUT'
    do L = 2,LA
      write(mpi_log_unit,1113) L,IL(L),JL(L),ZERO
    enddo
  endif
1112 FORMAT (2I5,2F12.4,6F12.7)
1113 FORMAT (3I5,F10.2)

  ! *** SET CORNER CELL STRESS CORRECTION
  do L = 2,LA
    FSCORTBCV(L) = 0.0
  enddo

  if( ISCORTBC >= 1 )then
    do L = 2,LA
      FSCORTBCV(L) = FSCORTBC
    enddo
  endif

  if( ISCORTBC == 2 .and. DEBUG )then
    if( process_id == master_id )then
      write(*,'(A)')'READING CORNERC.INP'
      open(1,FILE = 'cornerc.inp')
      
      call SKIPCOM(1,'*')
      read(1,*)NTMP

      do NT = 1,NTMP
        read(1,*)I,J,TMPVAL
        L = LIJ(I,J)
        FSCORTBCV(L) = TMPVAL
      enddo
      close(1)
    endif

    call Broadcast_Array(FSCORTBCV, master_id)
  endif

  ! *** READ SPATIAL AVERAGING MAP FOR FOOD CHAIN MODEL OUTPUT
  if( ISFDCH == 1 )then
    do L = 1,LC
      MFDCHZ(L) = 0
    enddo
    write(*,'(A)')'READING FOODCHAIN.INP'
    open(1,FILE = 'foodchain.inp')

    call SKIPCOM(1,'*')

    read(1,*)NFDCHIJ
    do LT = 1,NFDCHIJ
      read(1,*)I,J,ITMPVAL
      L = LIJ(I,J)
      MFDCHZ(L) = ITMPVAL
    enddo
    close(1)

  elseif( ISFDCH == 2 )then
    ! *** NEW BEDFORD APPLICATION
    call AllocateDSI( IFDCH, ICM, JCM, 0)

    write(*,'(A)')'READING FOODCHAIN.INP'
    open(1,FILE = 'foodchain.inp',STATUS = 'OLD')

    call SKIPCOM(1,'*')

    IFIRST = 1
    ILAST = IC
    do J = JC,1,-1
      read(1,66,IOSTAT = ISO)C,(IFDCH(I,J),I = IFIRST,ILAST)
      if( ISO > 0 )CALL STOPP('READ ERROR FOR FILE FOODCHAIN.INP')
      WRITE (7,166)JDUMY,(IFDCH(I,J),I = IFIRST,ILAST)
    enddo
    do L = 2,LA
      I = IL(L)
      J = JL(L)
      MFDCHZ(L) = IFDCH(I,J)
      if( NCORENO(I,J) == 7 )MFDCHZ(L) = 0
      if( MFDCHZ(L) > 1 )then
        do K = 1,KB
          STDOCB(L,K) = STDOCB(L,K)*0.10
        enddo
      endif
    enddo
    close(1)
66  FORMAT (I3,2X,113I1)
166 FORMAT (1X,I3,2X,113I1)
  endif

  ! *** READ IN COUNTER CLOCKWISE ANGLE FROM EAST SPECIFYING PRINCIPAL FLOOD FLOW DIRECTION
  if( ISTRAN(4) >= 1 .and. ISSFLFE >= 1 )then
    write(*,'(A)')'READING FLDANG.INP'
    open(1,FILE = 'fldang.inp',STATUS = 'UNKNOWN')
    do LL = 2,LA
      read(1,*,ERR = 3130)I,J,ANGTMP1,ANGTMP2
      L = LIJ(I,J)
      ACCWFLD(L,1) = 0.0174533*ANGTMP1
      ACCWFLD(L,2) = 0.0174533*ANGTMP2
    enddo
    close(1)
  endif

  GOTO 3132

3130 call STOPP('READ ERROR FOR FILE FLDANG.INP')
3132 continue

  !---------------------------------------------------------------------------!
  ! *** SET BOUNDARY CONDITION SWITCHES
  call SETBCS

  ! ****************************************************************************
  ! *** MPI communication
  ! *** Communicate SUB/SVB switches
  call communicate_ghost_cells(SUB, 'SUB')
  call communicate_ghost_cells(SVB, 'SVB')
  call communicate_ghost_cells(RSSBCE, 'RSSBCE')
  call communicate_ghost_cells(RSSBCW, 'RSSBCW')
  call communicate_ghost_cells(RSSBCN, 'RSSBCN')
  call communicate_ghost_cells(RSSBCS, 'RSSBCS')
  ! ****************************************************************************
  
  ! *** ALLOCATE MEMORY FOR VARIABLE TO STORE CONCENTRATIONS AT OPEN BOUNDARIES
  call AllocateDSI( WQBCCON,  NBCSOP, KCM, NACTIVEWC, 0.0)
  call AllocateDSI( WQBCCON1, NBCSOP, KCM, NACTIVEWC, 0.0)

  ! *** READ THE DRIFTER DATA
  if( ISPD > 0 ) CALL DRIFTER_INP

  ! *** CALCUATE CURVATURE METRICS (NEW ADDITION)
  do L = 1,LC
    DYDI(L) = 0.
    DXDJ(L) = 0.
  enddo

  ! ***DYDI
  TMPVAL = 0.
  do L = 2,LA
    I = IL(L)
    J = JL(L)

    ! *** KEEP ZERO'S FOR DYDI AT DOMAIN BORDERS
    if( I < 2 .or. I > IC-1 ) CYCLE
    if( J < 2 .or. J > JC-1 ) CYCLE

    ! *** DSI - CHANGED CELL TYPE 5 TO 8
    if( IJCT(I-1,J) >= 1 .and. IJCT(I-1,J) <= 8 )then
      if( IJCT(I+1,J) >= 1 .and. IJCT(I+1,J) <= 8 )then
        DYDI(L) = DYU(LEC(L))-DYU(L)
      else
        DDYDDDX = 2.*(DYP(L)-DYP(LWC(L)))/(DXP(L)+DXP(LWC(L)))
        DYUP1 = DYP(L)+0.5*DDYDDDX*DXP(L)
        DYDI(L) = DYUP1-DYU(L)
      endif
    else
      if( IJCT(I+1,J) >= 1 .and. IJCT(I+1,J) <= 8 )then
        DDYDDDX = 2.*(DYP(LEC(L))-DYP(L))/(DXP(LEC(L))+DXP(L))
        DYUM1 = DYP(L)-0.5*DDYDDDX*DXP(L)
        DYDI(L) = DYU(L)-DYUM1
      else
        DYDI(L) = 0.0
      endif
    endif
    if( DYDI(L) > 1.E-7 )then
      TMPVAL = 1.
    endif
  enddo

  ! ***DXDJ
  do L = 2,LA
    LN = LNC(L)
    LS = LSC(L)
    I = IL(L)
    J = JL(L)

    ! *** KEEP ZERO'S FOR DXDJ AT DOMAIN BORDERS
    if( I < 2 .or. I > IC-1 ) CYCLE
    if( J < 2 .or. J > JC-1 ) CYCLE

    ! *** DSI - CHANGED CELL TYPE 5 TO 8
    if( IJCT(I,J-1) >= 1 .and. IJCT(I,J-1) <= 8 )then
      if( IJCT(I,J+1) >= 1 .and. IJCT(I,J+1) <= 8 )then
        DXDJ(L) = DXV(LN)-DXV(L)
      else
        DDXDDDY = 2.*(DXP(L)-DXP(LS))/(DYP(L)+DYP(LS))
        DXVLN = DXP(L)+0.5*DDXDDDY*DYP(L)
        DXDJ(L) = DXVLN-DXV(L)
      endif
    else
      if( IJCT(I,J+1) >= 1 .and. IJCT(I,J+1) <= 8 )then
        DDXDDDY = 2.*(DXP(LN)-DXP(L))/(DYP(LN)+DYP(L))
        DXVLS = DXP(L)-0.5*DDXDDDY*DYP(L)
        DXDJ(L) = DXV(L)-DXVLS
      else
        DXDJ(L) = 0.0
      endif
    endif
    if( DXDJ(L) > 1.E-7 )then
      TMPVAL = 1.
    endif
  enddo

  ! ****************************************************************************
  ! *** MPI communication
  call communicate_ghost_cells(DYDI, 'DYDI')
  call communicate_ghost_cells(DXDJ, 'DXDJ')
  ! ****************************************************************************
  
  ! *** SETUP UP EDGE OF HARD BOTTOM REGIONS TO HANDLE BEDLOAD TRYING TO EXIT THE HARD BOTTOM REGION
  if( ISBEDMAP > 0 )then
    if( (ISTRAN(7) > 0 .and. ICALC_BL > 0) .or. (NSEDFLUME > 0 .and. ICALC_BL > 0) )then
      call AllocateDSI( IDX, LCM, 0)

      ! *** WEST
      BEDEDGEW.NEDGE = 0
      do LP = 1,LASED
        L = LSED(LP)
        LW = LWC(L)
        if( LW < LA .and. SUBO(L) > 0.5 .and. LBED(LW) )then
          BEDEDGEW.NEDGE = BEDEDGEW.NEDGE + 1
          IDX(BEDEDGEW.NEDGE) = L
        endif
      enddo
      if( BEDEDGEW.NEDGE > 0 )then
        call AllocateDSI( BEDEDGEW.LEDGE, BEDEDGEW.NEDGE, 0)
        do LP = 1,BEDEDGEW.NEDGE
          BEDEDGEW.LEDGE(LP) = IDX(LP)
        enddo
      endif

      ! *** EAST
      IDX = 0
      BEDEDGEE.NEDGE = 0
      do LP = 1,LASED
        L = LSED(LP)
        LE = LEC(L)
        if( LE < LA .and. SUBO(LE) > 0.5 .and. LBED(LE) )then
          BEDEDGEE.NEDGE = BEDEDGEE.NEDGE + 1
          IDX(BEDEDGEE.NEDGE) = L
        endif
      enddo
      if( BEDEDGEE.NEDGE > 0 )then
        call AllocateDSI( BEDEDGEE.LEDGE, BEDEDGEE.NEDGE, 0)
        do LP = 1,BEDEDGEE.NEDGE
          BEDEDGEE.LEDGE(LP) = IDX(LP)
        enddo
      endif

      ! *** NORTH
      IDX = 0
      BEDEDGEN.NEDGE = 0
      do LP = 1,LASED
        L = LSED(LP)
        LN = LNC(L)
        if( LN < LA .and. SVBO(LN) > 0.5 .and. LBED(LN) )then
          BEDEDGEN.NEDGE = BEDEDGEN.NEDGE + 1
          IDX(BEDEDGEN.NEDGE) = L
        endif
      enddo
      if( BEDEDGEN.NEDGE > 0 )then
        call AllocateDSI( BEDEDGEN.LEDGE, BEDEDGEN.NEDGE, 0)
        do LP = 1,BEDEDGEN.NEDGE
          BEDEDGEN.LEDGE(LP) = IDX(LP)
        enddo
      endif

      ! *** SOUTH
      IDX = 0
      BEDEDGES.NEDGE = 0
      do LP = 1,LASED
        L = LSED(LP)
        LS = LSC(L)
        if( LS < LA .and. SVBO(L) > 0.5 .and. LBED(LS) )then
          BEDEDGES.NEDGE = BEDEDGES.NEDGE + 1
          IDX(BEDEDGES.NEDGE) = L
        endif
      enddo
      if( BEDEDGES.NEDGE > 0 )then
        call AllocateDSI( BEDEDGES.LEDGE, BEDEDGES.NEDGE, 0)
        do LP = 1,BEDEDGES.NEDGE
          BEDEDGES.LEDGE(LP) = IDX(LP)
        enddo
      endif
      deallocate(IDX)
    endif
  endif

  ! **********************************************************************************
  ! *** EFDC CELL MAP LOGGING
  write(mpi_efdc_out_unit,'(A)') 'EFDC Cell Bottom Elevation Gradients'
  write(mpi_efdc_out_unit,'(3A6,4A10,4A14)') 'I','J','L','BELV','HP','DX','DY','GRADW','GRADE','GRADS','GRADN'
  do L = 2,LA
    GRADW = -999.
    GRADE = -999.
    GRADS = -999.
    GRADN = -999.
    if( SUBO(L) > 0.      ) GRADW = ( BELV(LWC(L))+HP(LWC(L) ) - ( BELV(L)+HP(L) )           )/DXU(L)
    if( SUBO(LEC(L)) > 0. ) GRADE = ( BELV(L)+HP(L)            - ( BELV(LEC(L))+HP(LEC(L)) ) )/DXU(LEC(L))
    if( SVBO(L) > 0.      ) GRADS = ( BELV(LSC(L))+HP(LSC(L) ) - ( BELV(L)+HP(L) )           )/DYV(L)
    if( SVBO(LNC(L)) > 0. ) GRADN = ( BELV(L)+HP(L)            - ( BELV(LNC(L))+HP(LNC(L)) ) )/DYV(LNC(L))
    write(mpi_efdc_out_unit,'(3I6,4F10.3,4F14.5)') IL(L),JL(L),L,BELV(L),HP(L),DXP(L),DYP(L),GRADW,GRADE,GRADS,GRADN
  enddo

  ! **********************************************************************************
  ! *** ACTIVE CELL LISTS
  call AllocateDSI( LLWET, KCM, -NDM,    0)
  call AllocateDSI( LKWET, LCM,  KCM, -NDM, 0)
  LLWET = 0
  LKWET = 0
  if( ISHDMF > 0 )then
    call AllocateDSI( LHDMF,  LCM,  KCM, .false.)
    call AllocateDSI( LLHDMF, KCM, -NDM,       0)
    call AllocateDSI( LKHDMF, LCM,  KCM,    -NDM,  0)

    NHDMF = LA-1
  endif

  ! *** ENSURE SUM OF DZC = 1.0000 TO MACHINE PRECISION
  DZPC = 0.0
  do K = 1,KC
    DZPC = DZPC + DZCK(K)
  enddo
  do K = 1,KC
    DZCK(K) = DZCK(K)/DZPC
  enddo

  ! ***************************************************************************
  ! *** BEGIN SIGMA-Z VARIABLE INITIALIZATION (SGZ)
  ! ***
  ! *** Use ORIGINAL DEPTH FROM DXDY (HMP) SO KSZ'S ARE CONSISTENT FOR COLD START,
  ! ***   RESTART AND CONTINUATION RUNS
  if( IGRIDV > 0 )then
    do L = 2,LA
      HMP(L) = HP(L) + SGZHPDELTA
    enddo
  endif

  MAXTHICK_local = -1
  do L = 2,LA
    if( HMP(L) > MAXTHICK_local ) MAXTHICK_local = HMP(L)
  enddo
    
  ! *** Should use MAXTHICK global to compute KSZ
  call DSI_All_Reduce(MAXTHICK_local, MAXTHICK, MPI_Max, TTDS, 1, TWAIT)

  ! *** INITIALIZE BOTTOM LAYERS
  if( IGRIDV == 0 )then
    ! *** Initialize global arrays for reporting
    do LG = 2,LA_Global
      KSZ_Global(LG) = 1
      DZC_Global(LG,:) = DZCK(K) 
    enddo
    
  else

    ! *** Read SGZ bottom active layer for Sigma-Z cells
    call AllocateDSI( R1D_Global, KCM, 0.0)
      
      if( process_id == master_id )then
        write(*,'(A)')'READING SGZLAYER.INP'
        open(1,FILE = 'sgzlayer.inp',STATUS = 'UNKNOWN')
        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES

      do LL = 2,LA_global
        if( IGRIDV == 1 )then
          read(1,*,END = 1000) IIN, JIN, K
        else
          read(1,*,END = 1000) IIN, JIN, K, (R1D_Global(KK),KK = 1,KC)
        endif
          
          LG = LIJ_Global(IIN,JIN)
          if( LG > 0 )then
            if( K < 1 )then
              write(6,'(A,3I5)')'*** ERROR: BAD KSZ IN SGZLAYER.INP, I,J,KSZ',I,J,K
              write(mpi_error_unit,'(A,3I5)')'*** ERROR: BAD KSZ IN SGZLAYER.INP, I,J,KSZ',I,J,K
              K = 1
            endif
            KSZ_Global(LG) = K
            DZPC = SUM(R1D_Global(:))                 ! *** Ensure the DZC's sum to 1.0
            if( DZPC <= 0.0 ) DZPC = 1.0
            DZC_Global(LG,:) = R1D_Global(:)/DZPC     ! *** Ensure the DZC's sum to 1.0
          else
            write(6,'(A,3I5)') '*** ERROR: BAD L INDEX IN SGZLAYER.INP, I,J,L',I,J,LG
            write(mpi_error_unit,'(A,3I5)') '*** ERROR: BAD L INDEX IN SGZLAYER.INP, I,J,L',I,J,LG
          endif
        enddo

1000    close(1)
      endif

      call MPI_Barrier(DSIcomm, ierr)
      call Broadcast_Array(DZC_Global, master_id)
      call Broadcast_Array(KSZ_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        I = Map2Local(LG).IG
        J = Map2Local(LG).JG
        KSZ(L) = KSZ_Global(LG)
        if( IGRIDV > 1 )then
          do K = 1,KC
            DZC(L,K) = DZC_Global(LG,K)
            if( K < KSZ(L) .and. DZC(L,K) /= 0.0 )then
              write(6,'(A,2I5,I8,I5)') '*** ERROR: DZC(L,K) FROM SGZLAYER.INP - NOT 0.0  @  I,J,L,K: ',I,J,LG,K
              call STOPP('.')
            endif
            if( K >= KSZ(L) .and. DZC(L,K) < 5E-5 )then
              write(6,'(A,2I5,I8,I5)') '*** ERROR: DZC(L,K) FROM SGZLAYER.INP - TOO SMALL  @   I,J,L,K: ',I,J,LG,K
              call STOPP('.')
            endif
          enddo
        endif
      endif
    enddo
      
    deallocate(R1D_Global)

    ! *** Communicate the local to global arrays
    call Assign_Loc_Glob_For_Write(1, size(KSZ,1), KSZ, size(KSZ_Global), KSZ_Global)
    call Assign_Loc_Glob_For_Write(2, size(DZC,1), size(DZC,2), DZC, size(DZC_Global,1), size(DZC_Global,2), DZC_Global)
    call Handle_Calls_MapGatherSort(2)

    if( process_id == master_id )then
      ! *** EXPORT THE LAYERING, AS COMPUTED OR ENTERED.
      open(1,FILE = OUTDIR//'SGZLAYER.OUT')
      write(1,'(A)') 'C ** SGZLAYER.OUT - LIST OF THE BOTTOMMOST ACTIVE LAYER'
      write(1,'(A)') 'C **    IGRIDV =  1 - BOTTOM LAYER PROVIDED,  2 - UNIFORM LAYERING DETERMINED FROM DZCK,  3 - UNIFORM LAYERING DETERMINED BY USER DEFINED THICKNESSES'
      write(1,'(A)') 'C **  '
      write(1,'(A,1000I10)') 'C **     L     I     J  KSZ', (K,K = 1,KC) 
      do LG = 2,LA_Global
        write(1,'(3X,I7,2I6,I5,1000F11.7)') LG, IL_GL(LG), JL_GL(LG), KSZ_Global(LG), (DZC_Global(LG,K),K = 1,KC)
      enddo
      close(1)
    endif

  endif

  ! *** SET VERTICAL GRID DEPENDENT ARRAYS AND LAYER MASKS
  LSGZU = .FALSE.
  LSGZV = .FALSE.
  do L = 2,LA
    LW = LWC(L)
    LE = LEC(L)
    LS = LSC(L)
    LN = LNC(L)

    DZPC = 0.
    do K = KSZ(L),KC
      DZPC = DZPC + DZCK(K)
    enddo

    do K = 1,KC
      SUB3D(L,K) = SUB(L)
      SVB3D(L,K) = SVB(L)

      ! *** INACTIVE LAYER MASK FOR CURRENT CELL
      if( K < KSZ(L) )then
        LKSZ(L,K) = .TRUE.
        DZC(L,K)  = 0.0
      elseif( IGRIDV > 0 )then
        ! *** SIGMA-ZED DZC
        LKSZ(L,K) = .FALSE.
        if( IGRIDV == 1 ) DZC(L,K) = DZCK(K)/DZPC
      else
        ! *** SIGMA STRETCH DZC
        LKSZ(L,K) = .FALSE.
        DZC(L,K)  = DZCK(K)/DZPC
      endif

      ! *** U FACE
      if( SUBO(L) > 0. )then
        KM = max(KSZ(LW), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR U FACE
        if( K >= KM )then
          LSGZU(L,K) = .TRUE.
        else
          SUB3D(L,K) = 0.0
        endif
      endif

      ! *** V FACE
      if( SVBO(L) > 0. )then
        KM = max(KSZ(LS), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR V FACE
        if( K >= KM )then
          LSGZV(L,K) = .TRUE.
        else
          SVB3D(L,K) = 0.0
        endif
      endif

      ! *** SET FACE/LAYER FLAGS
      SUB3D(L,K) = SUB3D(L,K)*SUB(L)
      SVB3D(L,K) = SVB3D(L,K)*SVB(L)
      SUB3DO(L,K) = SUB3D(L,K)
      SVB3DO(L,K) = SVB3D(L,K)
    enddo

    if( SUBO(L) < 0.5 )then
      KSZU(L) = KSZ(L)
    else
      KSZU(L) = max(KSZ(L),KSZ(LW))
    endif
    if( SVBO(L) < 0.5 )then
      KSZV(L) = KSZ(L)
    else
      KSZV(L) = max(KSZ(L),KSZ(LS))
    endif

    do K = KSZ(L),KC
      if( DZC(L,K) < 1.E-5 )then
        write(*,'(A,3I5,F10.4,I5)') 'ERROR: BAD LAYER THICKNESS FOR L,K,KSZ,HP = ',L,K,KSZ(L),HMP(L),process_id
        write(mpi_error_unit,'(A,3I5,F10.4)') 'ERROR: BAD LAYER THICKNESS FOR L,K,KSZ,HP = ',L,K,KSZ(L),HMP(L)
        call STOPP('.')
      endif
      DZIC(L,K) = 1./DZC(L,K)
    enddo

    DZIG(L,0) = 0.
    DZIG(L,KC) = 0.
    do K = KSZ(L),KS
      DZG(L,K)  = 0.5*(DZC(L,K)+DZC(L,K+1))
      DZIG(L,K) = 1./DZG(L,K)
    enddo

    CDZKMK(L,KSZ(L)) = 0.
    do K = KSZ(L)+1,KC
      CDZKMK(L,K) = DZIG(L,K-1)*DZIC(L,K)
    enddo
    do K = KSZ(L),KS
      CDZKK(L,K) = DZIC(L,K)*DZIG(L,K)
      CDZKKP(L,K) = DZIG(L,K)*DZIC(L,K+1)
    enddo
    CDZKK(L,KC) = 0.

    Z(L,0) = 0.
    do K = KSZ(L),KC
      Z(L,K)  = Z(L,K-1) +     DZC(L,K)   ! *** TOP OF LAYER Z
      ZZ(L,K) = Z(L,K)   - 0.5*DZC(L,K)   ! *** MID LAYER Z

      ! *** WALL PROXIMITY FUNCTION
      if( IFPROX == 0 ) FPROX(L,K) = 0.
      if( K /= KC )then
        ! *** Mellor Yamada, 1982 (parabolic shape)
        if( IFPROX == 1 ) FPROX(L,K) = 1./(VKC*Z(L,K))**2 + 1./(VKC*(1. - Z(L,K)))**2
        ! Blumberg et al., 1992 (correction for open channel flow)
        if( IFPROX == 2 ) FPROX(L,K) = (1./(VKC*Z(L,K))**2) + CTE5*(1./(VKC*(1.-Z(L,K)))**2)/(CTE4+0.00001)
        ! *** Burchard et al., 1998 (symmetric linear shape)
        if( IFPROX == 3 ) FPROX(L,K) = 1./(VKC*MIN(Z(L,K),(1. - Z(L,K))))**2
        ! *** Burchard et al., 2001 (simulation of an infinitely deep basin)
        if( IFPROX == 4 ) FPROX(L,K) = 1./(VKC*(1. - Z(L,K)))**2
      endif
    enddo

  enddo
  
  ! ****************************************************************************
  ! *** MPI communication
  call Communicate_Ghost_LCM0(SUB3DO)         ! *** Handles special case of having a zero as the starting index for LCM
  call Communicate_Ghost_LCM0(SVB3DO)         ! *** Handles special case of having a zero as the starting index for LCM
  call Communicate_Ghost_LCM0(SUB3D)          ! *** Handles special case of having a zero as the starting index for LCM
  call Communicate_Ghost_LCM0(SVB3D)          ! *** Handles special case of having a zero as the starting index for LCM
  call communicate_ghost_cells(KSZU, 'KSZU')
  call communicate_ghost_cells(KSZV, 'KSZV')
  call communicate_ghost_cells(LSGZU)
  call communicate_ghost_cells(LSGZV)
  ! ****************************************************************************
    
  ! *** CALCULATE CONSTANT HORIZONTAL SPATIAL ARRAYS
  do L = 2,LA
    DXYU(L)  = DXU(L)*DYU(L)
    DXYV(L)  = DXV(L)*DYV(L)
    DXYP(L)  = STCAP(L)*DXP(L)*DYP(L)
    DXIU(L)  = 1./DXU(L)
    DYIU(L)  = 1./DYU(L)
    DXIV(L)  = 1./DXV(L)
    DYIV(L)  = 1./DYV(L)
    DXYIP(L) = 1./(STCAP(L)*DXP(L)*DYP(L))
    DXYIU(L) = 1./(DXU(L)*DYU(L))
    DXYIV(L) = 1./(DXV(L)*DYV(L))
  enddo
  
  ! *** Read cell lists for comparison of single domain runs to MPI multiple domain runs
  if( active_domains == 1 )then
    INQUIRE(FILE = 'domain_list_000.txt', EXIST = RES)
    if( RES )then
      ! *** Sub-domain cell lists are available.  Read them
      open(5000, file = 'domain_list_000.txt')    
      read(5000,*) size_mpi
      close(5000)
      
      ! *** Loop over each sub-domain and get active list
      do i = 0, size_mpi-1
        write(mpi_list, '(A12,I3.3,A4)') 'domain_list_', i, '.txt'
        open(5000, file = mpi_list, status = 'unknown')
        
        read(5000,*) nd                ! *** Skip it
        do lg = 2, la_global
          read(5000,*) lp, nd, ls
          domain_list(lg,nd+1) = ls
        enddo
        close(5000)
      enddo
    endif
  else
    ! *** Save the current sub-domains cell list
    write(mpi_list, '(A12,I3.3,A4)') 'domain_list_', process_id, '.txt'
    open(5000+process_id, file = OUTDIR//mpi_list, status = 'unknown')
    
    write(5000+process_id,'(i10)') active_domains
    do lg = 2, la_global
      domain_list(lg,1) = lg
      write(5000+process_id,'(3i6)') lg, process_id, map2local(lg).ll
    enddo
    close(5000)
  endif
  
  ! *** SET ACTIVE CELL LIST (WITHOUT WET/DRY CELL CONSIDERATION)
  do ND = 1,NDM
    LF = 2 + (ND-1)*LDM
    LL = min(LF+LDM-1,LA)
    do K = 1,KC
      LN = 0
      do L = LF,LL
        if( K < KSZ(L) )CYCLE
        LN = LN+1
        LKWET(LN,K,ND) = L
      enddo
      LLWET(K,ND) = LN
    enddo
  enddo

  ! *** Limit BC layers to active layers
  do NWR = 1,NQWR
    IU = WITH_RET(NWR).IQWRU
    JU = WITH_RET(NWR).JQWRU
    LU = LIJ(IU,JU)
    WITH_RET(NWR).KQWRU = max(WITH_RET(NWR).KQWRU, KSZ(LU))
    
    ID = WITH_RET(NWR).IQWRD
    JD = WITH_RET(NWR).JQWRD
    LD = LIJ(ID,JD)
    if( LD > 1 )then
      WITH_RET(NWR).KQWRD = max(WITH_RET(NWR).KQWRD, KSZ(LD))  
    endif
  enddo
  do NJP = 1,NQJPIJ
    LU = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)
    JET_PLM(NJP).KUPCJP = max(JET_PLM(NJP).KUPCJP, KSZ(LU))
  enddo
  
  If( process_id == master_id )then
    if( IGRIDV > 0 )then
      if( DEBUG ) WRITE(6,'(A)') 'TOTAL NUMBER OF ACTIVE COMPUTATIONAL CELLS PER LAYER'
      write(mpi_efdc_out_unit,'(//,A)') 'TOTAL NUMBER OF ACTIVE COMPUTATIONAL CELLS PER LAYER'
      do K = KC,1,-1
        if( DEBUG ) WRITE(6,8000) '  Layer, LLWET:',K,SUM(LLWET(K,1:NDM))
        write(mpi_efdc_out_unit,8000) '  Layer, LLWET:',K,SUM(LLWET(K,1:NDM))
      enddo
      write(mpi_efdc_out_unit,8000) ' '
    endif
  endif !***End calc on master process

8000 FORMAT(A,I5,3I10)

  if( ISHDMF > 0 )then
    LLHDMF = LLWET
    LKHDMF = LKWET
  endif

  ! **********************************************************************************

  ! *** READ RESTART CONDITIONS OR INITIALIZE SCALAR FIELDS
  !     ISRESTI == 10 READS AND OLD RESTART FILE GENERATED BY
  !     PRE SEPTEMBER 8, 1992 VERSIONS OF EFDC.FOR
  if( ISRESTI >= 1 ) CALL Restart_In

  ! *** INITIALIZE SALINITY FIELD IF NOT READ IN FROM RESTART FILE
  if( ISTRAN(1) >= 1 .and. (ISRESTI == 0 .or.  &
     (ISRESTI >= 1 .and. ISCI(1) == 0 )  .or.  &
     (ISTOPT(1) > 1)) )then  ! *** PMC SINGLE LINE - FORCE IC
    if( ISTOPT(1) >= 1 )then
      NREST = 0
      do K = 1,KC
        do L = 2,LA
          SAL(L,K) = SALINIT(L,K)
          SAL1(L,K) = SALINIT(L,K)
        enddo
      enddo

      M = 1
      do K = 1,KC
        do LL = 1,NCBS
          L = LCBS(LL)
          CLOS(LL,K,M) = SALINIT(L,K)
          NLOS(LL,K,M) = 0
          if( NCSERS(LL,1) == 0 )then
            SAL(L,K) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
            SAL1(L,K) = SAL(L,K)
          endif
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBW
          L = LCBW(LL)
          CLOW(LL,K,M) = SALINIT(L,K)
          NLOW(LL,K,M) = 0
          if( NCSERW(LL,1) == 0 )then
            SAL(L,K) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
            SAL1(L,K) = SAL(L,K)
          endif
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBE
          L = LCBE(LL)
          CLOE(LL,K,M) = SALINIT(L,K)
          NLOE(LL,K,M) = 0
          if( NCSERE(LL,1) == 0 )then
            SAL(L,K) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
            SAL1(L,K) = SAL(L,K)
          endif
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBN
          L = LCBN(LL)
          CLON(LL,K,M) = SALINIT(L,K)
          NLON(LL,K,M) = 0
          if( NCSERN(LL,1) == 0 )then
            SAL(L,K) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
            SAL1(L,K) = SAL(L,K)
          endif
        enddo
      enddo

    endif
  endif
9101 FORMAT(I5)
9102 FORMAT(3I5,12F8.2)

  ! *** INITIALIZE TEMP FIELD IF NOT READ IN FROM RESTART FILE
  if( ISTRAN(2) >= 1 .and. (ISRESTI == 0                      .or.  &
    (ISRESTI >= 1 .and. ISCI(2) == 0 ) .or.  &
    (ISTOPT(2) > 9)) )then  ! *** PMC SINGLE LINE - FORCE IC
    ! *** SPATIALLY VARYING TEMPERATURE FIELD
    NREST = 0
    do K = 1,KC
      do L = 2,LA
        TEM(L,K)  = TEMINIT(L,K)
        TEM1(L,K) = TEM(L,K)
      enddo
    enddo

    M = 2
    do K = 1,KC
      do LL = 1,NCBS
        L = LCBS(LL)
        CLOS(LL,K,M) = TEMINIT(L,K)
        NLOS(LL,K,M) = 0
        if( NCSERS(LL,2) == 0 )then
          TEM(L,K) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
          TEM1(L,K) = TEM(L,K)
        endif
      enddo
    enddo
    do K = 1,KC
      do LL = 1,NCBW
        L = LCBW(LL)
        CLOW(LL,K,M) = TEMINIT(L,K)
        NLOW(LL,K,M) = 0
        if( NCSERW(LL,2) == 0 )then
          TEM(L,K) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
          TEM1(L,K) = TEM(L,K)
        endif
      enddo
    enddo
    do K = 1,KC
      do LL = 1,NCBE
        L = LCBE(LL)
        CLOE(LL,K,M) = TEMINIT(L,K)
        NLOE(LL,K,M) = 0
        if( NCSERE(LL,2) == 0 )then
          TEM(L,K) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
          TEM1(L,K) = TEM(L,K)
        endif
      enddo
    enddo
    do K = 1,KC
      do LL = 1,NCBN
        L = LCBN(LL)
        CLON(LL,K,M) = TEMINIT(L,K)
        NLON(LL,K,M) = 0
        if( NCSERN(LL,2) == 0 )then
          TEM(L,K) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
          TEM1(L,K) = TEM(L,K)
        endif
      enddo
    enddo

  endif
  
  ! *** INITIALIZE TEMPERATURE BC IF NOT READ IN FROM RESTART FILE
  !     AND CONSTANT INITIAL CONDITION IS USED
  if( ISRESTI == 0 .and. ISTRAN(2) >= 1 )then
    if( ISTOPT(2) == 0 )then
      ! *** CONSTANT TEMPERATURE FIELD
      M = 2
      do K = 1,KC
        do LL = 1,NCBS
          CLOS(LL,K,M) = TEMO
          NLOS(LL,K,M) = 0
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBW
          CLOW(LL,K,M) = TEMO
          NLOW(LL,K,M) = 0
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBE
          CLOE(LL,K,M) = TEMO
          NLOE(LL,K,M) = 0
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBN
          CLON(LL,K,M) = TEMO
          NLON(LL,K,M) = 0
        enddo
      enddo
    endif
  endif

  ! *** RESET IC OPTION, IF USED
  if( ISTOPT(2) > 9) ISTOPT(2) = ISTOPT(2)-10 ! PMC SINGLE LINE

  !
  ! *** INITIALIZE DYE FIELD 
  if( ISTRAN(3) >= 1 )then
    BFLAG = ISRESTI == 0 .or. (ISRESTI == 1 .and. ISCI(3) == 0) .or. ISTOPT(3) > 1   ! *** Cold start flag
    if( BFLAG )then
      do MD = 1,NDYE
        do K = 1,KC
          do L = 2,LA
            DYE(L,K,MD)  = DYEINIT(L,K,MD)
            DYE1(L,K,MD) = DYE(L,K,MD)
          enddo
        enddo

        M = 2 + MD
        do K = 1,KC
          do LL = 1,NCBS
            L = LCBS(LL)
            CLOS(LL,K,M) = DYEINIT(L,K,MD)
            NLOS(LL,K,M) = 0
            if( NCSERS(LL,3) == 0 )then
              DYE(L,K,MD) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              DYE1(L,K,MD) = DYE(L,K,MD)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            L = LCBW(LL)
            CLOW(LL,K,M) = DYEINIT(L,K,MD)
            NLOW(LL,K,M) = 0
            if( NCSERW(LL,3) == 0 )then
              DYE(L,K,MD) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              DYE1(L,K,MD) = DYE(L,K,MD)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            L = LCBE(LL)
            CLOE(LL,K,M) = DYEINIT(L,K,MD)
            NLOE(LL,K,M) = 0
            if( NCSERE(LL,3) == 0 )then
              DYE(L,K,MD) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              DYE1(L,K,MD) = DYE(L,K,MD)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            L = LCBN(LL)
            CLON(LL,K,M) = DYEINIT(L,K,MD)
            NLON(LL,K,M) = 0
            if( NCSERN(LL,3) == 0 )then
              DYE(L,K,MD) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              DYE1(L,K,MD) = DYE(L,K,MD)
            endif
          enddo
        enddo
      enddo
    endif
  endif

  ! *** INITIALIZE SFL if( ISRESTI == 0.AND ISTRAN(4) >= 1 )
  if( ISRESTI == 0 .and. ISTRAN(4) >= 1 )then
    if( ISTOPT(4) == 11 )then
      do K = 1,KC
        do L = 1,LC
          SFL(L,K) = SFLINIT(L,K)
          SFL2(L,K) = SFLINIT(L,K)
        enddo
      enddo

      M = 3 + NDYM
      do K = 1,KC
        do LL = 1,NCBS
          L = LCBS(LL)
          CLOS(LL,K,M) = SFLINIT(L,K)
          NLOS(LL,K,M) = 0
          SFL(L,K) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
          SFL2(L,K) = SFL(L,K)
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBW
          L = LCBW(LL)
          CLOW(LL,K,M) = SFLINIT(L,K)
          NLOW(LL,K,M) = 0
          SFL(L,K) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
          SFL2(L,K) = SFL(L,K)
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBE
          L = LCBE(LL)
          CLOE(LL,K,M) = SFLINIT(L,K)
          NLOE(LL,K,M) = 0
          SFL(L,K) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
          SFL2(L,K) = SFL(L,K)
        enddo
      enddo
      do K = 1,KC
        do LL = 1,NCBN
          L = LCBN(LL)
          CLON(LL,K,M) = SFLINIT(L,K)
          NLON(LL,K,M) = 0
          SFL(L,K) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
          SFL2(L,K) = SFL(L,K)
        enddo
      enddo
    endif
  endif

  ! *** INITIALIZE TOX AND BC IF NOT READ IN FROM RESTART FILE
  ! *** AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(5) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      if( ITXINT(NT) == 1 .or. ITXINT(NT) == 3 )then
        do K = 1,KC
          do L = 2,LA
            TOX(L,K,NT) = TOXINIT(L,K,NT)
            TOX1(L,K,NT) = TOX(L,K,NT)
          enddo
        enddo

        M = MSVTOX(NT)
        do K = 1,KC
          do LL = 1,NCBS
            L = LCBS(LL)
            CLOS(LL,K,M) = TOXINIT(L,K,NT)
            NLOS(LL,K,M) = 0
            if( NCSERS(LL,4) == 0 )then
              TOX(L,K,NT) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              TOX1(L,K,NT) = TOX(L,K,NT)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            L = LCBW(LL)
            CLOW(LL,K,M) = TOXINIT(L,K,NT)
            NLOW(LL,K,M) = 0
            if( NCSERW(LL,5) == 0 )then
              TOX(L,K,NT) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              TOX1(L,K,NT) = TOX(L,K,NT)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            L = LCBE(LL)
            CLOE(LL,K,M) = TOXINIT(L,K,NT)
            NLOE(LL,K,M) = 0
            if( NCSERE(LL,5) == 0 )then
              TOX(L,K,NT) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              TOX1(L,K,NT) = TOX(L,K,NT)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            L = LCBN(LL)
            CLON(LL,K,M) = TOXINIT(L,K,NT)
            NLON(LL,K,M) = 0
            if( NCSERN(LL,5) == 0 )then
              TOX(L,K,NT) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              TOX1(L,K,NT) = TOX(L,K,NT)
            endif
          enddo
        enddo

      endif
    enddo
  endif

  ! *** INITIALIZE TOX BC IF NOT READ IN FROM RESTART FILE
  ! *** AND CONSTANT INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(5) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      if( ITXINT(NT) == 0 .or. ITXINT(NT) == 2 )then
        M = MSVTOX(NT)
        do K = 1,KC
          do LL = 1,NCBS
            CLOS(LL,K,M) = TOXINTW(NT)
            NLOS(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            CLOW(LL,K,M) = TOXINTW(NT)
            NLOW(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            CLOE(LL,K,M) = TOXINTW(NT)
            NLOE(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            CLON(LL,K,M) = TOXINTW(NT)
            NLON(LL,K,M) = 0
          enddo
        enddo
      endif
    enddo
  endif

  ! *** INITIALIZE TOX BED IF NOT READ IN FROM RESTART FILE
  ! *** AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(5) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      if( ITXINT(NT) == 2 .or. ITXINT(NT) == 3 )then
        do K = 1,KB
          do L = 2,LA
            TOXB(L,K,NT) = TOXBINIT(L,K,NT)
            TOXB1(L,K,NT) = TOXB(L,K,NT)
          enddo
        enddo
      endif
    enddo
  endif

  ! *** INITIALIZE SEDIMENTS FOR WATER COLUMN AND BC'S, IF NOT READ IN FROM RESTART FILE
  ! *** AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(6) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(6) > 0 )then
    if( ISEDINT == 1 .or. ISEDINT == 3 )then
      do NS = 1,NSED
        do K = 1,KC
          do L = 2,LA
            SED(L,K,NS) = SEDINIT(L,K,NS)
            SED1(L,K,NS) = SED(L,K,NS)
          enddo
        enddo

        M = MSVSED(NS)
        do K = 1,KC
          do LL = 1,NCBS
            L = LCBS(LL)
            CLOS(LL,K,M) = SEDINIT(L,K,NS)
            NLOS(LL,K,M) = 0
            if( NCSERS(LL,6) == 0 )then
              SED(L,K,NS) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              SED1(L,K,NS) = SED(L,K,NS)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            L = LCBW(LL)
            CLOW(LL,K,M) = SEDINIT(L,K,NS)
            NLOW(LL,K,M) = 0
            if( NCSERW(LL,6) == 0 )then
              SED(L,K,NS) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              SED1(L,K,NS) = SED(L,K,NS)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            L = LCBE(LL)
            CLOE(LL,K,M) = SEDINIT(L,K,NS)
            NLOE(LL,K,M) = 0
            if( NCSERE(LL,6) == 0 )then
              SED(L,K,NS) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              SED1(L,K,NS) = SED(L,K,NS)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            L = LCBN(LL)
            CLON(LL,K,M) = SEDINIT(L,K,NS)
            NLON(LL,K,M) = 0
            if( NCSERN(LL,6) == 0 )then
              SED(L,K,NS) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              SED1(L,K,NS) = SED(L,K,NS)
            endif
          enddo
        enddo

      enddo

    endif
  endif
  if( ISTRAN(6) > 0 ) DEALLOCATE(SEDINIT)

  ! *** INITIALIZE SEDIMENT WATER COLUMN BC'S, IF NOT READ IN FROM RESTART FILE AND
  ! *** CONSTANT INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(6) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(6) > 0 )then
    if( ISEDINT == 0 .or. ISEDINT == 2 )then
      do NS = 1,NSED
        M = MSVSED(NS)
        do K = 1,KC
          do LL = 1,NCBS
            CLOS(LL,K,M) = SEDO(NS)
            NLOS(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            CLOW(LL,K,M) = SEDO(NS)
            NLOW(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            CLOE(LL,K,M) = SEDO(NS)
            NLOE(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            CLON(LL,K,M) = SEDO(NS)
            NLON(LL,K,M) = 0
          enddo
        enddo
      enddo
    endif
  endif

  ! *** INITIALIZE SED BED IF NOT READ IN FROM RESTART FILE
  ! *** AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(6) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(6) > 0 )then
    if( ISEDINT == 2 .or. ISEDINT == 3 )then
      do NS = 1,NSED
        do K = 1,KB
          do L = 2,LA
            SEDB(L,K,NS) = SEDBINIT(L,K,NS)
            SEDB1(L,K,NS) = SEDB(L,K,NS)
          enddo
        enddo
      enddo
    endif
  endif

  ! *** INITIALIZE SND AND BC IF NOT READ IN FROM RESTART FILE
  ! *** AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(7) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(7) > 0 )then
    if( ISEDINT == 1 .or. ISEDINT == 3 )then
      do NS = 1,NSND
        do K = 1,KC
          do L = 2,LA
            SND(L,K,NS) = SNDINIT(L,K,NS)
            SND1(L,K,NS) = SND(L,K,NS)
          enddo
        enddo

        M = MSVSND(NS)
        do K = 1,KC
          do LL = 1,NCBS
            L = LCBS(LL)
            CLOS(LL,K,M) = SNDINIT(L,K,NS)
            NLOS(LL,K,M) = 0
            if( NCSERS(LL,7) == 0 )then
              SND(L,K,NS) = WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              SND1(L,K,NS) = SND(L,K,NS)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            L = LCBW(LL)
            CLOW(LL,K,M) = SNDINIT(L,K,NS)
            NLOW(LL,K,M) = 0
            if( NCSERW(LL,7) == 0 )then
              SND(L,K,NS) = WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              SND1(L,K,NS) = SND(L,K,NS)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            L = LCBE(LL)
            CLOE(LL,K,M) = SNDINIT(L,K,NS)
            NLOE(LL,K,M) = 0
            if( NCSERE(LL,7) == 0 )then
              SND(L,K,NS) = WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              SND1(L,K,NS) = SND(L,K,NS)
            endif
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            L = LCBN(LL)
            CLON(LL,K,M) = SNDINIT(L,K,NS)
            NLON(LL,K,M) = 0
            if( NCSERN(LL,7) == 0 )then
              SND(L,K,NS) = WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              SND1(L,K,NS) = SND(L,K,NS)
            endif
          enddo
        enddo

      enddo
    endif
  endif
  if( ISTRAN(7) > 0 ) DEALLOCATE(SNDINIT)

  ! *** INITIALIZE SND BC IF NOT READ IN FROM RESTART FILE AND
  ! *** CONSTANT INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(7) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(7) > 0 )then
    if( ISEDINT == 0 .or. ISEDINT == 2 )then
      do NX = 1,NSND
        NS = NSED + NX
        M = MSVSND(NX)
        do K = 1,KC
          do LL = 1,NCBS
            CLOS(LL,K,M) = SEDO(NS)
            NLOS(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBW
            CLOW(LL,K,M) = SEDO(NS)
            NLOW(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBE
            CLOE(LL,K,M) = SEDO(NS)
            NLOE(LL,K,M) = 0
          enddo
        enddo
        do K = 1,KC
          do LL = 1,NCBN
            CLON(LL,K,M) = SEDO(NS)
            NLON(LL,K,M) = 0
          enddo
        enddo
      enddo
    endif
  endif

  ! *** INITIALIZE SND BED, IF NOT READ IN FROM RESTART FILE
  ! *** AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP = 1
  if( ISRESTI == 0 ) IISTMP = 0
  if( ISRESTI >= 1 .and. ISCI(7) == 0 ) IISTMP = 0
  if( IISTMP == 0 .and. ISTRAN(7) > 0 )then
    if( ISEDINT == 2 .or. ISEDINT == 3 )then
      do NX = 1,NSND
        do K = 1,KB
          do L = 2,LA
            SNDB(L,K,NX) = SNDBINIT(L,K,NX)
            SNDB1(L,K,NX) = SNDB(L,K,NX)
          enddo
        enddo
      enddo
    endif
  endif

  ! *** INITIALIZE SEDIMENT BED
  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 ) CALL BEDINIT
  if( ISTRAN(6) > 0 ) DEALLOCATE(SEDBINIT)
  if( ISTRAN(7) > 0 ) DEALLOCATE(SNDBINIT)

  ! *** SET THE WET CELL LIST
  LAWET = 0
  LADRY = 0
  do L = 2,LA
    LAWET = LAWET+1
    LWET(LAWET) = L
  enddo
  LDMWET = INT(LAWET/NTHREADS)+1
  LDMDRY = 0

  ! *** ACTIVATE DYE TRACER CONTINUITY CHECK
  if( ISMMC == 1 .and. ISTRAN(3) > 0 )then
    do MD = 1,NDYE
      do K = 1,KC
        do L = 1,LC
          DYE(L,K,MD) = 1.
          DYE1(L,K,MD) = 1.
        enddo
      enddo

      M = 2 + MD
      do K = 1,KC
        do LL = 1,NCBS
          CLOS(LL,K,M) = 1.
          NLOS(LL,K,M) = 0
        enddo
        do LL = 1,NCBW
          CLOW(LL,K,M) = 1.
          NLOW(LL,K,M) = 0
        enddo
        do LL = 1,NCBE
          CLOE(LL,K,M) = 1.
          NLOE(LL,K,M) = 0
        enddo
        do LL = 1,NCBN
          CLON(LL,K,M) = 1.
          NLON(LL,K,M) = 0
        enddo
      enddo
    enddo
  endif
  
  ! *** SET WATER COLUMN BYPASS LIMITS
  IP = 0
  if( ISTRAN(1) > 0 )then
    IP = IP + 1
  endif
  if( ISTRAN(2) > 0 )then
    IP = IP + 1
  endif
  if( ISTRAN(3) > 0 )then
    do MD = 1,NDYE
      IP = IP + 1
      WCV(IP).WCLIMIT = DYES(MD).WCLIMIT
    enddo
  endif
  if( ISTRAN(4) > 0 )then
    IP = IP + 1
  endif
  if( ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      IP = IP + 1
      WCV(IP).WCLIMIT = TOXS(NT).WCLIMIT
    enddo
  endif
  if( ISTRAN(6) > 0 )then
    do NS = 1,NSED
      IP = IP + 1
      WCV(IP).WCLIMIT = SEDS(NS).WCLIMIT
    enddo
  endif
  if( ISTRAN(7) > 0 )then
    do NS = 1,NSND
      IP = IP + 1
      WCV(IP).WCLIMIT = SEDS(NSED+NS).WCLIMIT
    enddo
  endif
  if( ISTRAN(8) > 0 )then
    do MW = 1,NWQV
      if( ISTRWQ(MW) == 1 )then
        IP = IP + 1
        WCV(IP).WCLIMIT = 0.0
      endif
    enddo
  endif
  
  ! *** SET 3D CELL FACE CONSTANTS
  FSGZU = 1.
  FSGZV = 1.
  do L = 2,LA
    LW = LWC(L)
    LS = LSC(L)

    do K = 1,KC
      ! *** U FACE
      if( SUBO(L) > 0. )then
        KM = max(KSZ(LW), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR U FACE
        if( K >= KM )then
          if( IGRIDV > 0 )then
            if( KSZ(LW) > KSZ(L) )then
              SGZU(L,K)  = DZC(LW,K)
            else
              SGZU(L,K)  = DZC(L,K)
            endif
          else
            SGZU(L,K)  = max(DZC(LW,K),DZC(L,K))
          endif
        endif
      endif

      ! *** V FACE
      if( SVBO(L) > 0. )then
        KM = max(KSZ(LS), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR V FACE
        if( K >= KM )then
          if( IGRIDV > 0 )then
            if( KSZ(LS) > KSZ(L) )then
              SGZV(L,K)  = DZC(LS,K)
            else
              SGZV(L,K)  = DZC(L,K)
            endif
          else
            SGZV(L,K)  = max(DZC(LS,K),DZC(L,K))
          endif
        endif
      endif
    enddo

    if( IGRIDV > 1 )then
      if( SUBO(L) > 0.5 )then
        TMP = SUM(SGZU(L,1:KC))
        do K = 1,KC
          SGZU(L,K) = SGZU(L,K) / TMP
        enddo
      endif

      if( SVBO(L) > 0.5 )then
        TMP = SUM(SGZV(L,1:KC))
        do K = 1,KC
          SGZV(L,K) = SGZV(L,K) / TMP
        enddo
      endif

      ! *** QC
      if( SUB(L) > 0. .and. ABS(1.0 - SUM(SGZU(L,1:KC))) > 1.E-5 )then
        PRINT *,'BAD SGZU:  L, KSZ, SUM(SGZU(L,1:KC))',L,KSZ(L),SUM(SGZU(L,1:KC))
        do K = 1,KC
          PRINT *,K,SGZU(L,K),DZC(L,K)
        enddo
        call STOPP('BAD SGZU')
      endif
      if( SVB(L) > 0. .and. ABS(1.0 - SUM(SGZV(L,1:KC))) > 1.E-5 )then
        PRINT *,'BAD SGZV:  L, KSZ, SUM(SGZV(L,1:KC))',L,KSZ(L),SUM(SGZV(L,1:KC))
        do K = 1,KC
          PRINT *,K,SGZU(L,K),DZC(L,K)
        enddo
        call STOPP('BAD SGZV')
      endif
      ! *** END QC
    endif

    ! *** CELL FACE/INTERFACE CONSTANTS
    do K = KSZU(L),KS
      DZGU(L,K) = 0.5*(SGZU(L,K)+SGZU(L,K+1))
      if( SUBO(L) > 0. )then
        CDZFU(L,K) = SGZU(L,K)*SGZU(L,K+1)/(SGZU(L,K)+SGZU(L,K+1))
        CDZUU(L,K) = -SGZU(L,  K)/(SGZU(L,K)+SGZU(L,K+1))
        CDZLU(L,K) = -SGZU(L,K+1)/(SGZU(L,K)+SGZU(L,K+1))
      endif
    enddo
    do K = KSZV(L),KS
      DZGV(L,K) = 0.5*(SGZV(L,K)+SGZV(L,K+1))
      if( SVBO(L) > 0. )then
        CDZFV(L,K) = SGZV(L,K)*SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
        CDZUV(L,K) = -SGZV(L,K)  /(SGZV(L,K)+SGZV(L,K+1))
        CDZLV(L,K) = -SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
      endif
    enddo

    ! *** U FACE
    if( SUBO(L) > 0. )then
      CDZRU(L,KSZU(L)) = SGZU(L,KSZU(L))-1.
      CDZDU(L,KSZU(L)) = SGZU(L,KSZU(L))
      do K = KSZU(L)+1,KS
        KM = max(KSZ(LW), KSZ(L))      ! *** MINIMUM ACTIVE LAYERS FOR U FACE
        if( K >= KM )then
          CDZRU(L,K) = CDZRU(L,K-1)+SGZU(L,K)
          CDZDU(L,K) = CDZDU(L,K-1)+SGZU(L,K)
        endif
      enddo
    endif

    ! *** V FACE
    if( SVBO(L) > 0. )then
      CDZRV(L,KSZV(L)) = SGZV(L,KSZV(L))-1.
      CDZDV(L,KSZV(L)) = SGZV(L,KSZV(L))
      do K = KSZV(L)+1,KS
        KM = max(KSZ(LS), KSZ(L))      ! *** MINIMUM ACTIVE LAYERS FOR V FACE
        if( K >= KM )then
          CDZRV(L,K) = CDZRV(L,K-1)+SGZV(L,K)
          CDZDV(L,K) = CDZDV(L,K-1)+SGZV(L,K)
        endif
      enddo
    endif

    do K = 1,KS
      ! *** U FACE
      if( SUBO(L) > 0. )then
        CDZRU(L,K) = CDZRU(L,K)*DZGU(L,K)*CDZLU(L,KSZU(L))
        CDZMU(L,K) = 0.5*SGZU(L,K)*SGZU(L,K+1)
      endif

      ! *** V FACE
      if( SVBO(L) > 0. )then
        CDZRV(L,K) = CDZRV(L,K)*DZGV(L,K)*CDZLV(L,KSZV(L))
        CDZMV(L,K) = 0.5*SGZV(L,K)*SGZV(L,K+1)
      endif
    enddo

  enddo

  ! *** SET LIST OF CELLS WHOSE ACTIVE LAYER COUNT IS 1
  LDMSGZ1 = 0
  if( IGRIDV > 0 )then
    call AllocateDSI( LSGZ1, LCM, 0)
    do L = 2,LA
      if( KSZ(L) == KC )then
        LASGZ1 = LASGZ1 + 1
        LSGZ1(LASGZ1) = L
      endif
    enddo
    LDMSGZ1 = INT(LASGZ1/NDM) + 1
  endif

  ! ***  TREAT SGZ CELLS WHOSE KSZ = KC (i.e. KMIN = 1)
  call AllocateDSI( LLWETZ, KCM, -NDM,       0)
  call AllocateDSI( LKWETZ, LCM,  KCM, -NDM, 0)

  ! ****************************************************************************
  ! *** MPI communication - Communicate face values
  call communicate_ghost_cells(SGZU)
  call communicate_ghost_cells(SGZV) 
  
  call communicate_ghost_cells(CDZDU)
  call communicate_ghost_cells(CDZDV)
  
  call communicate_ghost_cells(CDZFU)
  call communicate_ghost_cells(CDZFV)
  
  call communicate_ghost_cells(CDZLU)
  call communicate_ghost_cells(CDZLV)

  call communicate_ghost_cells(CDZRU)
  call communicate_ghost_cells(CDZRV)
  
  call communicate_ghost_cells(CDZUU)
  call communicate_ghost_cells(CDZUV)
  ! ****************************************************************************

  do ND = 1,NDM
    do K = 1,KS
      LLWETZ(K,ND) = LLWET(K,ND)
      do LP = 1,LLWET(K,ND)
        LKWETZ(LP,K,ND) = LKWET(LP,K,ND)
      enddo
    enddo

    LLWETZ(KC,ND) = LLWET(KS,ND)
    do LP = 1,LLWET(KS,ND)
      LKWETZ(LP,KC,ND) = LKWET(LP,KS,ND)
    enddo
  enddo

  ! *** CALCULATE WET/DRY HORIZONTAL SPATIAL ARRAYS
  do L = 2,LA
    HRU(L)   = SUB(L)*HMU(L)*DYU(L)*DXIU(L)
    HRV(L)   = SVB(L)*HMV(L)*DXV(L)*DYIV(L)
    HRUO(L)  = SUBO(L)*DYU(L)*DXIU(L)
    HRVO(L)  = SVBO(L)*DXV(L)*DYIV(L)
    SBX(L)   = 0.5*SUB(L)*DYU(L)
    SBY(L)   = 0.5*SVB(L)*DXV(L)
    SBXO(L)  = 0.5*SUBO(L)*DYU(L)
    SBYO(L)  = 0.5*SVBO(L)*DXV(L)
  enddo
  
  ! ****************************************************************************
  ! *** MPI communication - Communicate face values
  call communicate_ghost_cells(SBX, 'SBX')
  call communicate_ghost_cells(SBY, 'SBY')
  call communicate_ghost_cells(SBY, 'SBXO')
  call communicate_ghost_cells(SBY, 'SBYO')
  call communicate_ghost_cells(SBY, 'HRUO')
  call communicate_ghost_cells(SBY, 'HRVO')
  ! ****************************************************************************

  ! *** Determine FSGZU/FSGZV for gross momentum
  do L = 2,LA
    LW = LWC(L)
    LS = LSC(L)
    do K = 1,KC
      if( SGZU(L,K) > 0. )then
        FSGZU(L,K) = 1./SGZU(L,K)
      else
        FSGZU(L,K) = 0.0
      endif
      if( SGZV(L,K) > 0. )then
        FSGZV(L,K) = 1./SGZV(L,K)
      else
        FSGZV(L,K) = 0.0
      endif
    enddo
  enddo

  ! *** Third pass at cell constants
  do L = 2,LA
    LW = LWC(L)
    LE = LEC(L)
    LS = LSC(L)
    LN = LNC(L)

    if( IGRIDV > 0 )then
      ! *** Cell interface metrics
      FRACK = min(0.1*MAXTHICK/REAL(KC), 0.25)

      BELVW(L) = BELV(L)
      BELVE(L) = BELV(L)
      BELVS(L) = BELV(L)
      BELVN(L) = BELV(L)

      KSZW(L) = KSZ(L)
      KSZE(L) = KSZ(L)
      KSZS(L) = KSZ(L)
      KSZN(L) = KSZ(L)

      do K = KSZ(L),KC
        SGZW(L,K) = DZC(L,K)
        SGZE(L,K) = DZC(L,K)
        SGZS(L,K) = DZC(L,K)
        SGZN(L,K) = DZC(L,K)

        SGZKW(K,L) = DZC(L,K)
        SGZKE(K,L) = DZC(L,K)
        SGZKS(K,L) = DZC(L,K)
        SGZKN(K,L) = DZC(L,K)

        ! *** Top of layer
        ZW(L,K) = Z(L,K)
        ZE(L,K) = Z(L,K)
        ZS(L,K) = Z(L,K)
        ZN(L,K) = Z(L,K)

        ! *** Middle of layer
        ZZW(K,L) = ZZ(L,K)
        ZZE(K,L) = ZZ(L,K)
        ZZS(K,L) = ZZ(L,K)
        ZZN(K,L) = ZZ(L,K)
      enddo

      if( SUBO(L) > 0. )then
        KSZW(L) = KSZU(L)
        if( KSZ(LW) > KSZ(L) )then
          BELVW(L) = BELV(LW)
          !BELVW(L) = max(BELV(LW)-FRACK,BELV(L))
          do K = 1,KC
            SGZW(L,K)  = DZC(LW,K)
            SGZKW(K,L) = DZC(LW,K)
            ZW(L,K)    = Z(LW,K)
            ZZW(K,L)   = ZZ(LW,K)
          enddo
        endif
        if( IGRIDV > 1 .and. KSZ(LW) == KSZ(L) .and. KSZ(L) < KC-KMINV+1 )then
          if( BELV(L) >= BELV(LW) )then
            LL = L
          else
            LL = LW
          endif
          BELVW(L) = BELV(LL)
          do K = 1,KC
            SGZW(L,K)  = DZC(LL,K)
            SGZKW(K,L) = DZC(LL,K)
            ZW(L,K)    = Z(LL,K)
            ZZW(K,L)   = ZZ(LL,K)
          enddo
        endif
      endif

      if( SUBO(LE) > 0. )then
        KSZE(L) = KSZU(LE)
        if( KSZ(LE) > KSZ(L) )then
          BELVE(L) = BELV(LE)
          !BELVE(L) = max(BELV(LE)-FRACK,BELV(L))
          do K = 1,KC
            SGZE(L,K)  = DZC(LE,K)
            SGZKE(K,L) = DZC(LE,K)
            ZE(L,K)    = Z(LE,K)
            ZZE(K,L)   = ZZ(LE,K)
          enddo
        elseif( IGRIDV > 1 .and. KSZ(LE) == KSZ(L) .and. KSZ(L) < KC-KMINV+1 )then
          if( BELV(L) >= BELV(LE) )then
            LL = L
          else
            LL = LE
          endif
          BELVE(L) = BELV(LL)
          do K = 1,KC
            SGZE(L,K)  = DZC(LL,K)
            SGZKE(K,L) = DZC(LL,K)
            ZE(L,K)    = Z(LL,K)
            ZZE(K,L)   = ZZ(LL,K)
          enddo
        endif
      endif

      if( SVBO(L) > 0. )then
        KSZS(L) = KSZV(L)
        if( KSZ(LS) > KSZ(L) )then
          BELVS(L) = BELV(LS)
          !BELVS(L) = max(BELV(LS)-FRACK,BELV(L))
          do K = 1,KC
            SGZS(L,K)  = DZC(LS,K)
            SGZKS(K,L) = DZC(LS,K)
            ZS(L,K)    = Z(LS,K)
            ZZS(K,L)   = ZZ(LS,K)
          enddo
        elseif( IGRIDV > 1 .and. KSZ(LS) == KSZ(L) .and. KSZ(L) < KC-KMINV+1 )then
          if( BELV(L) >= BELV(LS) )then
            LL = L
          else
            LL = LS
          endif
          BELVS(L) = BELV(LL)
          do K = 1,KC
            SGZS(L,K)  = DZC(LL,K)
            SGZKS(K,L) = DZC(LL,K)
            ZS(L,K)    = Z(LL,K)
            ZZS(K,L)   = ZZ(LL,K)
          enddo
        endif
      endif

      if( SVBO(LN) > 0. )then
        KSZN(L) = KSZV(LN)
        if( KSZ(LN) > KSZ(L) )then
          BELVN(L) = BELV(LN)
          !BELVN(L) = max(BELV(LN)-FRACK,BELV(L))
          do K = 1,KC
            SGZN(L,K)  = DZC(LN,K)
            SGZKN(K,L) = DZC(LN,K)
            ZN(L,K)    = Z(LN,K)
            ZZN(K,L)   = ZZ(LN,K)
          enddo
        elseif( IGRIDV > 1 .and. KSZ(LN) == KSZ(L) .and. KSZ(L) < KC-KMINV+1 )then
          if( BELV(L) >= BELV(LN) )then
            LL = L
          else
            LL = LN
          endif
          BELVN(L) = BELV(LL)
          do K = 1,KC
            SGZN(L,K)  = DZC(LL,K)
            SGZKN(K,L) = DZC(LL,K)
            ZN(L,K)    = Z(LL,K)
            ZZN(K,L)   = ZZ(LL,K)
          enddo
        endif
      endif
    endif

    do K = KSZ(L),KS
      if( SUBO(L) > 0.0 )then
        if( LSGZU(L,K) )then
          DZGTMP = 0.5*(SGZU(L,K)+SGZU(L,K+1))
          DZGTMP = 1./DZGTMP
          DZIGSD4U(L,K) = 0.25*DZGTMP*DZGTMP
        else
          DZIGSD4U(L,K) = 0.0
        endif
      else
        DZIGSD4U(L,K) = 0.25*DZIG(L,K)*DZIG(L,K)
      endif

      if( SVBO(L) > 0.0 )then
        if( LSGZV(L,K) )then
          DZGTMP = 0.5*(SGZV(L,K)+SGZV(L,K+1))
          DZGTMP = 1./DZGTMP
          DZIGSD4V(L,K) = 0.25*DZGTMP*DZGTMP
        else
          DZIGSD4V(L,K) = 0.
        endif
      else
        DZIGSD4V(L,K) = 0.25*DZIG(L,K)*DZIG(L,K)
      endif
    enddo

  enddo

  ! *** Zero momentum switches for cells along active sub-domain boundaries. (MPI)
  do L = 1, LC  
    if( IL(L) == IC )then
      FSGZU(L,:) = 0.0
      FSGZV(L,:) = 0.0
    endif
    if( JL(L) == JC )then
      FSGZU(L,:) = 0.0
      FSGZV(L,:) = 0.0
    endif
  enddo
  
  ! ****************************************************************************
  ! *** MPI communication - Communicate face values
  call communicate_ghost_cells(SBX, 'SBX')
  call communicate_ghost_cells(SBY, 'SBY')
  if( IGRIDV > 0 )then
    call communicate_ghost_cells(BELVW, 'BELVW')
    call communicate_ghost_cells(BELVE, 'BELVE')
    call communicate_ghost_cells(BELVS, 'BELVS')
    call communicate_ghost_cells(BELVN, 'BELVN')
  
    call communicate_ghost_cells(KSZW, 'KSZW')
    call communicate_ghost_cells(KSZE, 'KSZE')
    call communicate_ghost_cells(KSZS, 'KSZS')
    call communicate_ghost_cells(KSZN, 'KSZN')
  
    call communicate_ghost_cells(SGZW)
    call communicate_ghost_cells(SGZE)
    call communicate_ghost_cells(SGZS)
    call communicate_ghost_cells(SGZN)
  
    call communicate_ghost_3d0(ZW)
    call communicate_ghost_3d0(ZE)
    call communicate_ghost_3d0(ZS)
    call communicate_ghost_3d0(ZN)
  endif
  ! ****************************************************************************

  if( IGRIDV > 0 )then
    do L = 2,LA
      HMP(L) = HMP(L) - SGZHPDELTA

      ! *** SET FACE DEPTHS
      HPW(L) = HP(L) + BELV(L) - BELVW(L)
      HPE(L) = HP(L) + BELV(L) - BELVE(L)
      HPS(L) = HP(L) + BELV(L) - BELVS(L)
      HPN(L) = HP(L) + BELV(L) - BELVN(L)
    enddo
  endif

  ! *** COMPUTE HU/HV FOR INITIAL CONDITIONS
  if( ISRESTI == 0 )then
    do L = 2,LA
      LW = LWC(L)
      LS = LSC(L)
    
      if( KSZ(LW) > KSZ(L) )then
        HU(L) = max( 0.5*HPK(L,KSZ(LW)), HP(LW)*(1.+DZC(L,KSZ(LW))*0.1) )
      elseif( KSZ(LW) < KSZ(L) )then
        HU(L) = max( 0.5*HPK(LW,KSZ(L)), HP(L)*(1.+DZC(LW,KSZ(L))*0.1) )
      else
        HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )/(DXYP(L) + DXYP(LW))
      endif
    
      if( KSZ(LS) > KSZ(L) )then
        HV(L) = max( 0.5*HPK(L,KSZ(LS)), HP(LS)*(1.+DZC(L,KSZ(LS))*0.1) )
      elseif( KSZ(LS) < KSZ(L) )then
        HV(L) = max( 0.5*HPK(LS,KSZ(L)), HP(L)*(1.+DZC(LS,KSZ(L))*0.1) )
      else
        HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )/(DXYP(L) + DXYP(LS))
      endif
    
      HPI(L) = 1./HP(L)
      HUI(L) = 1./HU(L)
      HVI(L) = 1./HV(L)
    
    enddo
  endif

  ! *** INITIALIZE BUOYANCY AND EQUATION OF STATE - ALL CELLS
  N = 0

  ! *** Make sure to check IF density effects are turned on and bouyancy is intialized correctly
  if( BSC > 1.E-6 )then
    call CALBUOY(.FALSE.)

    ! ****************************************************************************
    ! *** MPI communication
    call communicate_ghost_cells(B)
    if( ISTRAN(2) > 0 .and. ISICE == 4 )then
      call communicate_ghost_cells(RHOW)
    endif
    ! ****************************************************************************
  endif

  ! *** COMPUTATIONAL CELL LIST BY SUB-DOMAIN AND LAYER (WET/DRY CONSIDERED)
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = min(LF+LDM-1,LA)
    do K = 1,KC
      LN = 0
      do L = LF,LL
        if( LKSZ(L,K) )CYCLE
        if( LMASKDRY(L) )then
          LN = LN+1
          LKWET(LN,K,ND) = L
        endif
      enddo
      LLWET(K,ND) = LN
    enddo
  enddo

  ! *** COMPUTATIONAL CELL LIST FOR ENTIRE DOMAIN, i.e. ND = 0, BY LAYER (WET/DRY CONSIDERED)
  do K = 1,KC
    LN = 0
    do L = 2,LA
      if( LKSZ(L,K) )CYCLE
      if( LMASKDRY(L) )then
        LN = LN+1
        LKWET(LN,K,0) = L  ! *** Wet Cell for Layer K
      endif
    enddo
    LLWET(K,0) = LN        ! *** Total Wet Cells for Layer K
  enddo

  do K = 1,KC
    do L = 2,LA
      HPK(L,K) = HP(L)*DZC(L,K)
    enddo
  enddo
  if( ISRESTI > 0 )then
    do K = 1,KC
      do L = 2,LA
        H1PK(L,K) = H1P(L)*DZC(L,K)
        H2PK(L,K) = H1PK(L,K)
      enddo
    enddo
  endif

  TMPVAL = max(HMIN,.01)
  HPKI = 1./TMPVAL
  do K = 1,KC
    do L = 2,LA
      if( K < KSZ(L) )CYCLE
      HPKI(L,K) = 1./HPK(L,K)
    enddo
  enddo

  if( ISHDMF > 0 )then
    ! *** READ THE HMD SUBSET LIST, IF NEEDED
    if( IHMDSUB > 0 )then
      write(*,'(A)')'READING MAPHMD.INP'
      open(1,FILE = 'maphmd.inp',STATUS = 'UNKNOWN')
      LINE = READSTR(1) ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*,IOSTAT = ISO) NHDMF
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPHMD.INP')
      do NP = 1,NHDMF
        read(1,*,IOSTAT = ISO) L,I,J
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPHMD.INP')
        do K = KSZ(L),KC
          LHDMF(L,K) = .TRUE.
        enddo
      enddo
      close(1)
      PRINT '(A,I10)','  NUMBER OF HMD CELLS FOUND: ',NHDMF
    else
      do L = 2,LA
        do K = KSZ(L),KC
          LHDMF(L,K) = .TRUE.
        enddo
      enddo
    endif

    do ND = 1,NDM
      do K = 1,KC
        LN = 0
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          if( LHDMF(L,K) )then
            LN = LN+1
            LKHDMF(LN,K,ND) = L
          endif
        enddo
        LLHDMF(K,ND) = LN     ! *** NUMBER OF WET HDMF CELLS FOR THE CURRENT LAYER
      enddo
    enddo
  endif

  ! *** END OF SIGMA-Z VARIABLE INITIALIZATION (SGZ)
  ! ***************************************************************************

  ! *** DSI BEGIN SEEPAGE
  if( ISGWIT == 3 )then
    do L = 2,LA
      QGW(L) = QGW(L)*DXYP(L)                ! *** m3/s
    enddo
  endif

  ! *** LARGE ASPECT RATIO ASSIGMENTS
  NASPECT = 0
  if( XYRATIO > 1.1 )then
    call AllocateDSI( LASPECT, LCM, .false.)
    do L = 2,LA
      if( DXP(L) > XYRATIO*DYP(L) .or. DYP(L) > XYRATIO*DXP(L) )then
        NASPECT = NASPECT+1
        LASPECT(L) = .TRUE.
      endif
    enddo
  endif

  ! *** INITIALIZE ELEVATION OF ACTIVE GROUNDWATER ZONE FOR COLD START
  if( ISGWIE >= 1 .and. ISRESTI == 0 )then
    do L = 2,LA
      if( HP(L) > HDRY )then
        AGWELV(L) = BELV(L)
      else
        AGWELV(L) = BELAGW(L)
      endif
    enddo
    do L = 2,LA
      AGWELV1(L) = AGWELV(L)
      AGWELV2(L) = AGWELV(L)
    enddo
    open(1,FILE = OUTDIR//'GWELV.OUT',STATUS = 'UNKNOWN')
    write(1,5400)
    write(1,5402)
    do L = 2,LA
      write(1,5401)IL(L),JL(L),BELV(L),BELAGW(L),AGWELV(L)
    enddo
    close(1)
  endif
5400 FORMAT('   I   J    BELELV      BELAGW        AGWELV',/)
5401 FORMAT(1X,2I5,2X,F10.5,2X,F10.5,2X,F10.5)
5402 FORMAT(/)

  ! *** CALCULATE CONSTANT C ARRAYS FOR EXTERNAL P SOLUTION
  ! *** HRU = SUB*HMU*DYU/DXU & HRV = SVB*HMV*DXV/DYV
  ! *** DXYIP = 1/(DXP*DYP)
  if( IRVEC /= 9 )then
    do L = 2,LA
      CC(L) = 1.
      CCC(L) = 1.
    enddo
    if( ISRLID == 1 )then
      do L = 2,LA
        CC(L) = 0.
        CCC(L) = 0.
        if( SPB(L) == 0. ) CC(L) = 1.
        if( SPB(L) == 0. ) CCC(L) = 1.
      enddo
    endif
    do L = 2,LA
      LE = LEC(L)
      LN = LNC(L)
      C1 = -G*DT*DT*SPB(L)*DXYIP(L)
      CS(L) = C1*HRV(L)
      CW(L) = C1*HRU(L)
      CE(L) = C1*HRU(LE)
      CN(L) = C1*HRV(LN)
      CC(L) = CC(L)-CS(L)-CW(L)-CE(L)-CN(L)
      CCS(L) = 0.25*CS(L)
      CCW(L) = 0.25*CW(L)
      CCE(L) = 0.25*CE(L)
      CCN(L) = 0.25*CN(L)
      CCC(L) = CCC(L)-CCS(L)-CCW(L)-CCE(L)-CCN(L)
      CCCI(L) = 1./CCC(L)
    enddo
  endif

  ! *** *********************************************************************

  ! *** INITIALIZE WATER QUALITY MODEL AND READ INPUT
  if( ISTRAN(8) >= 1 ) CALL WQ3DINP

  ! *** Added call to some shell fish initialization
  if( ISFFARM > 0 .and. NSF > 0 )then
      do L = 2,LA
        call SHELLFISH_REDIST(L) 
        call SHELLFISH_LENGTH(L)
      enddo
  endif

  ! *** RELEASE ARRAYS NOT USED
  if( ISWASP /= 17 )then
    deallocate(DYEINIT)
    deallocate(SALINIT)
    deallocate(TEMINIT)
    deallocate(SFLINIT)
  endif

  ! *** *********************************************************************
  ! *** Call Propwash initialization if specified
  if( propwash_on )then
    
    ! *** Need to get the cell areas 
    if( .not. allocated(XCOR) ) CALL AREA_CENTRD
    
#ifdef DEBUGGING
    LDEBUG = .true.
#else
    LDEBUG = .false.
#endif

    ! *** Read the propwash values
    if( process_id == master_id )then
      call Read_propwash_config(LDEBUG)
    
      call Read_Ship_Data(LDEBUG)
    
      call Read_Ship_Tracks(LDEBUG)
    endif
    
    call Broadcast_Propwash
    
    ! *** Setup all of the ships
    call Setup_Ships(LDEBUG)
    
  endif
  
  ! *** Initialize propwash linkage, even if not used
  if( istran(6) > 0 .or. istran(7) > 0 )then
    call AllocateDSI( prop_ero, lcm, -NSEDS, 0.0)
    call AllocateDSI( prop_bld, lcm, -NSEDS, 0.0)
  endif  
  
  ! *** INITIALIZE EFDC EXPLORER OUTPUT (SKIP IF CONTINUATION)
  NITER = 0
  if( ISPPH == 1 .and. (ISRESTI == 0 .or. (ISRESTI /= 0 .and. ICONTINUE == 0)) )then
    ! ****************************************************************************
    ! *** Get local LA
    call Get_LA_Local_No_Ghost
    
    if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
        if( IHTSTRT > 0 )then
            ! *** do nothing
    
        else! KBT_Global is not read in from a file so we need to create it from the local KBTs
            ! *** Map the local KBTs to a global one
            call Map_Gather_Sort_Int_1D(size(KBT,1), KBT, size(KBT_Global,1), KBT_Global)
        endif
    endif
  
    call Map_Write_EE_Binary
    ! ****************************************************************************

    if( process_id == master_id )then
      call EE_LINKAGE(1)
    endif
  endif

  ! *** Close files no longer needed to avoid drive unavailable errors
  close(mpi_efdc_out_unit)  ! *** EFDC+ log file File
  close(mpi_log_unit)       ! *** MPI/EFDC Log File
  close(mpi_error_unit)     ! *** EFDC Error log file
  close(mpi_mapping_unit)   ! *** MPI Mapping File

  ! *** Setup processor specific optional log files
  if( ISDRY > 0 .and. NDRYSTP > 0 )then
    write(mpi_filename, '(A17,I3.3,A4)')    'log_qdwaste_proc_',process_id, '.log'
    Open(unit = mpi_qdwaste_unit, status = 'replace', file = OUTDIR//mpi_filename)
    mpi_qdwaste = mpi_filename
    write(STR, '("*** EFDC+ Filename: ",A,",  Version: ",A10," ***")') TRIM(EFDC_EXE), EFDC_VER
    write(mpi_qdwaste_unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
    close(mpi_qdwaste_unit)
  endif
  
  LB = 2                    ! *** LB is a globally declared integer that can be hardwired for dubugging a certain cell
  
  ! *** *********************************************************************
  ! *** SELECT FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION OR
  ! *** LONG-TERM MASS TRANSPORT CALCULATION  (DISABLED)
  NITERAT = 0
  
  if( IS2TIM == 0 ) CALL HDMT
  if( IS2TIM >= 1 ) CALL HDMT2T
  !IF(ISLTMT >= 1 ) CALL LTMT
  
  ! *** ENSURE FINAL SNAPSHOT HAS BEEN WRITTEN
  if( ISPPH == 1 .and. NSNAPSHOTS == NSNAPMAX )then
    !CALL BCOUT_ACCUMULATE(DTDYN) !@todo need to modify to owrk with mpi
    call Map_Write_EE_Binary
    
    If( process_id == master_id )then
      call EE_LINKAGE(0)
    endif

  endif

  ! *** *********************************************************************
  ! *** WRITE END OF RUN SUMMARIES TO SCREEN AND LOG FILES
  open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')

  ! *** REPORT THE MHK SUMMARY
  if( LMHK )then
    write(mpi_efdc_out_unit,'(//"***********************************************************************")')
    write(mpi_efdc_out_unit,'("  MHK SUMMARY")')
    MHKL = SUM(ESUP(1:TCOUNT,2:LA))
    write(mpi_efdc_out_unit,'("  SUPPORT ENERGY LOSS  ",E14.6," MW-HR")') MHKL
    MHKE = SUM(EMHK(1:TCOUNT,2:LA))
    write(mpi_efdc_out_unit,'("  MHK ENERGY OUTPUT    ",E14.6," MW-HR")') MHKE
    write(mpi_efdc_out_unit,'("***********************************************************************"//)')
  endif

  ! *** REPORT ANY SEDZLJ WARNINGS
  if( NWARNING > 0 )then
    write(mpi_efdc_out_unit,'(/A,I10)')   'SEDZLJ WARNING.  DURING THE SIMULATION THE ACTUAL BED SHEAR STRESS EXCEEDED'
    write(mpi_efdc_out_unit,'(A,I10,A/)') '                 THE RANGE OF DEFINED SHEAR CATAGORIES IN THE SEDFLUME DATA',NWARNING,' TIMES'
  endif
  write(mpi_efdc_out_unit,2) num_Processors, NTHREADS

  ! *** Sum information over all processes
  call DSI_All_Reduce(NWARNING, COUNT, MPI_SUM, TTDS, 1, TWAIT)
  if( LMHK )then
    call DSI_All_Reduce(MHKL, MHKL_Global, MPI_SUM, TTDS, 1, TWAIT)
    call DSI_All_Reduce(MHKE, MHKE_Global, MPI_SUM, TTDS, 1, TWAIT)
  endif

  If( process_id == master_id )then
    ! *** DISPLAY THE MHK SUMMARY
    if( LMHK )then
      write(6,'("SUPPORT ENERGY LOSS    ",F10.4," MW-HR")') MHKL_Global
      write(6,'("MHK ENERGY OUTPUT      ",F10.4," MW-HR")') MHKE_Global
      write(6,'(A)')'See EFDCLOG.OUT File'
    endif
    
    ! *** DISPLAY ANY SEDZLJ WARNINGS
    if( COUNT > 0 )then
      write(6,'(/A,I10)')   'SEDZLJ WARNING.  DURING THE SIMULATION THE ACTUAL BED SHEAR STRESS EXCEEDED'
      write(6,'(A,I10,A/)') '                 THE RANGE OF DEFINED SHEAR CATAGORIES IN THE SEDFLUME DATA',COUNT,' TIMES'
    endif
    
    ! *** PRINT THE NUMBER OF THREADS
2   FORMAT( '********************************************************************************',/, &
            '***     This EFDCPlus Run used:',I6,' MPI processor(s) and ',I2,' thread(s)      ***',/, &
            '********************************************************************************',/)
    write(6,2) num_Processors, NTHREADS
  endif          ! *** End of Master Thread block

  ! *** OUTPUT TIMING (OUTPUT TO TIME.LOG (UNIT 9) USED BY EFDC_EXPLORER)
  TIME_END = DSTIME(1)
  TIME_END = (TIME_END-TIME_START)/3600.

  ! *** DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS,BUT MAY NOT BE
  ! *** SUPPORTED ON OTHER SYSTEMS.
  TCPU = DTIME(CPUTIME)
  TCPU = TCPU/3600.    !/NTHREADS
  CPUTIME(1) = CPUTIME(1)/3600.   !/NTHREADS
  CPUTIME(2) = CPUTIME(2)/3600.   !/NTHREADS

  THDMT   = THDMT - TAVB - TTBXY - TCEXP - TPUV - TUVW - THMDF - TQQQ - TSADV - TVDIF - THEAT - TTSED - TSSTX - TLRPD - TWQKIN - TWQSED - TWQRPEM - TPROPW 
  TCONG   = TCONG/3600.             ! *** INCLUDED IN TPUV
  THDMT   = THDMT/3600.
  TPUV    = TPUV/3600.
  TTSED   = (TTSED-TSSTX)/3600.
  TSSTX   = TSSTX/3600.
  TCEXP   = TCEXP/3600.
  TAVB    = TAVB/3600.
  TUVW    = TUVW/3600.
  TQQQ    = TQQQ/3600.
  TTBXY   = TTBXY/3600.
  THEAT   = THEAT/3600.
  TVDIF   = TVDIF/3600.
  TSADV   = TSADV/3600.
  TLRPD   = TLRPD/3600.
  THMDF   = THMDF/3600.
  TMPIEE  = TMPIEE/3600.
  TWQKIN  = TWQKIN/3600.
  TWQRPEM = TWQRPEM/3600.
  TWQSED  = TWQSED/3600.
  TPROPW  = TPROPW/3600.

  DSITIMING(:)  = DSITIMING(:)/3600.
  DSITIMING(10) = SUM(DSITIMING(1:9))
  DSITIMING(11) = DSITIMING(11) + DSITIMING(12)
  TMPIGH        = DSITIMING(10) + DSITIMING(11)
  THDMT         = THDMT - TMPIGH - TMPIEE
  
  If( process_id == master_id )then ! *** Show only on master process
     call close_nc_lpt()
    !------------------------------------------------------------------------
    ! *** WRITE THE TIME SUMMARY TO THE SCREEN
    if( LSEDZLJ )then
      write(6,1995) THDMT    ,TTSED , TSSTX
    else
      write(6,1996) THDMT    ,TTSED , TSSTX
    endif
    write(6,1997) TPUV      ,TCONG
    write(6,1998) TCEXP     ,TAVB
    write(6,1999) TUVW      ,TQQQ
    write(6,2000) TTBXY     ,THEAT
    write(6,2001) TLRPD     ,TSADV
    write(6,2002) THMDF     ,TVDIF
    if( ISTRAN(8) > 0 )then
      write(6,2003) TWQKIN  ,TWQRPEM
      write(6,2004) TWQSED  ,0.
    endif
    if( ISPROPWASH > 0 )then
      write(6,2005) TPROPW  ,0.
    endif
    ! ****************************************************************************
    write(6,2006) TMPIEE    ,TMPIGH
    ! ****************************************************************************
    write(6,2007) CPUTIME(1),CPUTIME(2)
    write(6,2008) TIME_END  ,TCPU
    write(6,'(//,20(A8,F14.5,/))') (DSITIME(I),DSITIMING(I),I = 1,13)

    !------------------------------------------------------------------------
    ! *** WRITE THE TIME SUMMARY TO THE EFDC LOG
    if( LSEDZLJ )then
      write(mpi_efdc_out_unit,1995)THDMT    ,TTSED , TSSTX
    else
      write(mpi_efdc_out_unit,1996)THDMT    ,TTSED, TSSTX
    endif
    write(mpi_efdc_out_unit,1997) TPUV      ,TCONG
    write(mpi_efdc_out_unit,1998) TCEXP     ,TAVB
    write(mpi_efdc_out_unit,1999) TUVW      ,TQQQ
    write(mpi_efdc_out_unit,2000) TTBXY     ,THEAT
    write(mpi_efdc_out_unit,2001) TLRPD     ,TSADV
    write(mpi_efdc_out_unit,2002) THMDF     ,TVDIF
    write(mpi_efdc_out_unit,2003) TWQKIN    ,TWQRPEM
    write(mpi_efdc_out_unit,2004) TWQSED    ,0.
    write(mpi_efdc_out_unit,2005) TPROPW    ,0.
    write(mpi_efdc_out_unit,2006) TMPIEE    ,TMPIGH
    write(mpi_efdc_out_unit,2007) CPUTIME(1),CPUTIME(2)
    write(mpi_efdc_out_unit,2008) TIME_END,  TCPU

    !------------------------------------------------------------------------
1995 FORMAT('***TIMING (HOURS)',/, &
      'T HDMT ONLY  = ',F8.4,'  T SEDZLJ     = ',F8.4, F8.4)
1996 FORMAT('***TIMING (HOURS)',/, &
      'T HDMT ONLY  = ',F8.4,'  T SSEDTOX    = ',F8.4, F8.4)
1997 FORMAT('T CALPUV     = ',F8.4,'  T CONG GRAD  = ',F8.4)
1998 FORMAT('T EXPLICIT   = ',F8.4,'  T CALC AV    = ',F8.4)
1999 FORMAT('T CALC UVW   = ',F8.4,'  T TURB QQQ   = ',F8.4)
2000 FORMAT('T T&B SHEAR  = ',F8.4,'  T HEAT PRCS  = ',F8.4)
2001 FORMAT('T PART TRK   = ',F8.4,'  T ADV TRANSP = ',F8.4)
2002 FORMAT('T HORIZ DIF  = ',F8.4,'  T VERT DFUSN = ',F8.4)
2003 FORMAT('WQ KINETICS  = ',F8.4,'  WQ RPEM      = ',F8.4)
2004 FORMAT('WQ DIAGEN    = ',F8.4,'  NOT USED     = ',F8.4)
2005 FORMAT('T PROPWASH   = ',F8.4,'  NOT USED     = ',F8.4)
2006 FORMAT('MPI EE GATH  = ',F8.4,'  MPI COMMUNIC = ',F8.4)
2007 FORMAT('CPU USER     = ',F8.4,'  CPU SYSTEM   = ',F8.4)
2008 FORMAT('ELAPSED TIME = ',F8.4,'  CPU TIME     = ',F8.4)

    !------------------------------------------------------------------------
    ! *** WRITE THE TIME SUMMARY TO THE TIME LOG
    call TIMELOG(N,TIMEDAY,OUTDIR,TIME_END*3600._8)

    open(9,FILE = OUTDIR//'TIME.LOG',POSITION = 'APPEND')
    write(9,2) num_Processors, NTHREADS
    if( LSEDZLJ )then
      write(9,6995)THDMT    ,TTSED, TSSTX
    else
      write(9,6996)THDMT    ,TTSED, TSSTX
    endif
    write(9,6997) TPUV      ,TCONG
    write(9,6998) TCEXP     ,TAVB
    write(9,6999) TUVW      ,TQQQ
    write(9,7000) TTBXY     ,THEAT
    write(9,7001) TLRPD     ,TSADV
    write(9,7002) THMDF     ,TVDIF
    write(9,7003) TWQKIN    ,TWQRPEM
    write(9,7004) TWQSED    ,0.
    write(9,7005) TPROPW    ,0.
    write(9,7006) TMPIEE    ,TMPIGH
    write(9,7007) CPUTIME(1),CPUTIME(2)
    write(9,7008) TIME_END  ,TCPU
    write(9,'(//,20(A8,F14.5,/))') (DSITIME(I),DSITIMING(I),I = 1,16)

6995 FORMAT('T HDMT ONLY   = ',F14.4,'  T SEDZLJ      = ',F14.4,F14.4)
6996 FORMAT('T HDMT ONLY   = ',F14.4,'  T SSEDTOX     = ',F14.4,F14.4)
6997 FORMAT('T CALPUV      = ',F14.4,'  T CONG GRAD   = ',F14.4)
6998 FORMAT('T EXPLICIT    = ',F14.4,'  T CALC AV     = ',F14.4)
6999 FORMAT('T CALC UVW    = ',F14.4,'  T TURB QQQ    = ',F14.4)
7000 FORMAT('T T&B SHEAR   = ',F14.4,'  T HEAT PRCS   = ',F14.4)
7001 FORMAT('T PART TRK    = ',F14.4,'  T ADV TRANSP  = ',F14.4)
7002 FORMAT('T HORIZ DIFF  = ',F14.4,'  T VERT DFUSN  = ',F14.4)
7003 FORMAT('WQ KINETICS   = ',F14.4,'  WQ RPEM       = ',F14.4)
7004 FORMAT('WQ DIAGEN     = ',F14.4,'  NOT USED      = ',F14.4)
7005 FORMAT('T PROPWASH    = ',F14.4,'  NOT USED      = ',F14.4)
7006 FORMAT('MPI EE GATH   = ',F14.4,'  MPI COMMUNIC  = ',F14.4)
7007 FORMAT('CPU USER      = ',F14.4,'  CPU SYSTEM    = ',F14.4)
7008 FORMAT('ELAPSED TIME  = ',F14.4,'  CPU TIME      = ',F14.4)
    !------------------------------------------------------------------------

    ! *** CLOSE OUTPUT  FILES
    close(mpi_efdc_out_unit)  ! *** EFDC+ log file File
    close(9)

    ! *** CSV file for timing evaluations
    open(9,FILE = OUTDIR//'TIME.CSV',STATUS = 'REPLACE')
    write(9,'(A,/)') TITLE
    write(9,'("MPI processors = ,",I10,",, OMP Threads = ,",I10,///)') ,num_Processors, NTHREADS
    if( LSEDZLJ )then
      write(9,8995)THDMT    ,TTSED, TSSTX, DSITIME(8), DSITIMING(8)
    else
      write(9,8996)THDMT    ,TTSED, TSSTX, DSITIME(8), DSITIMING(8)
    endif
    write(9,8997) TPUV      ,TCONG,   DSITIME(1), DSITIMING(1)
    write(9,8998) TCEXP     ,TAVB,    DSITIME(2), DSITIMING(2)
    write(9,8999) TUVW      ,TQQQ,    DSITIME(3), DSITIMING(3)
    write(9,9000) TTBXY     ,THEAT,   DSITIME(4), DSITIMING(4)
    write(9,9001) TLRPD     ,TSADV,   DSITIME(5), DSITIMING(5)
    write(9,9002) THMDF     ,TVDIF,   DSITIME(6), DSITIMING(6)
    write(9,9003) TWQKIN    ,TWQRPEM, DSITIME(7), DSITIMING(7)
    write(9,9004) TWQSED    ,0.,      DSITIME(9), DSITIMING(9)
    write(9,9005) TMPIEE    ,TMPIGH,  DSITIME(10), DSITIMING(10)
    write(9,9006) CPUTIME(1),CPUTIME(2), DSITIME(11), DSITIMING(11)
    write(9,9007) TIME_END  ,TCPU,    DSITIME(12), DSITIMING(12)
    
    close(9)

8995 FORMAT('T HDMT ONLY   = ,',F14.5,',,  T SEDZLJ      = ,',F14.5,',',F14.5,',,',A,',',F14.7)
8996 FORMAT('T HDMT ONLY   = ,',F14.5,',,  T SSEDTOX     = ,',F14.5,',',F14.5,',,',A,',',F14.7)
8997 FORMAT('T CALPUV      = ,',F14.5,',,  T CONG GRAD   = ,',F14.5,',,,',A,',',F14.7)
8998 FORMAT('T EXPLICIT    = ,',F14.5,',,  T CALC AV     = ,',F14.5,',,,',A,',',F14.7)
8999 FORMAT('T CALC UVW    = ,',F14.5,',,  T TURB QQQ    = ,',F14.5,',,,',A,',',F14.7)
9000 FORMAT('T T&B SHEAR   = ,',F14.5,',,  T HEAT PRCS   = ,',F14.5,',,,',A,',',F14.7)
9001 FORMAT('T PART TRK    = ,',F14.5,',,  T ADV TRANSP  = ,',F14.5,',,,',A,',',F14.7)
9002 FORMAT('T HORIZ DIFF  = ,',F14.5,',,  T VERT DFUSN  = ,',F14.5,',,,',A,',',F14.7)
9003 FORMAT('WQ KINETICS   = ,',F14.5,',,  WQ RPEM       = ,',F14.5,',,,',A,',',F14.7)
9004 FORMAT('WQ DIAGEN     = ,',F14.5,',,  NOT USED      = ,',F14.5,',,,',A,',',F14.7)
9005 FORMAT('MPI EE GATH   = ,',F14.5,',,  MPI COMMUNIC  = ,',F14.5,',,,',A,',',F14.7)
9006 FORMAT('CPU USER      = ,',F14.5,',,  CPU SYSTEM    = ,',F14.5,',,,',A,',',F14.7)
9007 FORMAT('ELAPSED TIME  = ,',F14.5,',,  CPU TIME      = ,',F14.5,',,,',A,',',F14.7)
    !------------------------------------------------------------------------
    
    ! *** CLOSE FIELDS AND SHELLFISH
    call FREEFIELDS()
    !IF( ISFFARM > 0) CALL FREE_SHELLFISH()

#ifdef WASPOUT
    ! *** Move the WASP Linkage File to the output folder, if needed.
    if( ISWASP > 0 .or. ISRCA > 0 .or. ISICM > 0 )then
      INQUIRE(FILE = HYDFIL,EXIST = RES)
      if( RES )then
        ! *** LINKAGE FILE
        if( IHL_HANDLE /= 0 )then
          call hlclose(IHL_HANDLE,IERROR)
          IHL_HANDLE = 0
        endif
        TITLE = OUTDIR//'wasp\'//HYDFIL
        INQUIRE(FILE = TITLE,EXIST = RES)
        if( RES )then
          I = DELFILESQQ (TITLE)
          if( I == 1)RES = RENAMEFILEQQ(HYDFIL,TITLE)
        else
          RES = RENAMEFILEQQ(HYDFIL,TITLE)
        endif

        ! *** LOG/DEBUG FILE
        INQUIRE(FILE = OUTDIR//'wasp\log.txt',EXIST = RES)
        if( RES )then
          I = DELFILESQQ (OUTDIR//'wasp\log.txt')
          if( I == 1)RES = RENAMEFILEQQ('log.txt',OUTDIR//'wasp\log.txt')
        else
          RES = RENAMEFILEQQ('log.txt',OUTDIR//'wasp\log.txt')
        endif
      endif
    endif
#endif

    ! *** Handle Runtime Flag
    open(1,FILE = '0run',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')

  endif !***End calculation on master process

  ! ************************************************************************************************
  ! *** Report timing for every MPI processor
  open(mpi_log_unit,FILE = OUTDIR//mpi_log_file,POSITION = 'APPEND')
  
  write(mpi_log_unit,'(/,80("*"),//,A,/)') TITLE
  write(mpi_log_unit,'("MPI Processes = ",I10,"   OMP Threads = ",I10,/,"MPI Process   = ",I10,//)') ,num_Processors, NTHREADS, process_id
  if( LSEDZLJ )then
    write(mpi_log_unit,6995)THDMT    ,TTSED, TSSTX
  else
    write(mpi_log_unit,6996)THDMT    ,TTSED, TSSTX
  endif
  write(mpi_log_unit,6997) TPUV      ,TCONG
  write(mpi_log_unit,6998) TCEXP     ,TAVB
  write(mpi_log_unit,6999) TUVW      ,TQQQ
  write(mpi_log_unit,7000) TTBXY     ,THEAT
  write(mpi_log_unit,7001) TLRPD     ,TSADV
  write(mpi_log_unit,7002) THMDF     ,TVDIF
  write(mpi_log_unit,7003) TWQKIN    ,TWQRPEM
  write(mpi_log_unit,7004) TWQSED    ,0.
  write(mpi_log_unit,7005) TPROPW    ,0.
  write(mpi_log_unit,7006) TMPIEE    ,TMPIGH
  write(mpi_log_unit,7007) CPUTIME(1),CPUTIME(2)
  write(mpi_log_unit,7008) TIME_END  ,TCPU
  write(mpi_log_unit,'(//,20(A8,F14.5,/))') (DSITIME(I),DSITIMING(I),I = 1,16)

  close(mpi_log_unit)
  
  ! *** Do some additional processing of the timing routines
  ! *** Fine the Min/Max to aid in load balancing 
  call Report_Max_Min_Timing(TPUV,    'CALPUV  ')
  call Report_Max_Min_Timing(TCEXP,   'EXPLICIT')
  call Report_Max_Min_Timing(TAVB,    'CALC AV ')
  call Report_Max_Min_Timing(TQQQ,    'TURB QQQ')
  call Report_Max_Min_Timing(TTBXY,   'TB SHEAR')
  call Report_Max_Min_Timing(TSADV,   'ADV TRAN')
  call Report_Max_Min_Timing(THMDF,   'CALHDMF ')
  call Report_Max_Min_Timing(TVDIF,   'VERTDIFF')
  call Report_Max_Min_Timing(THEAT,   'CALHEAT ')
  call Report_Max_Min_Timing(TTSED,   'SED TRAN')
  call Report_Max_Min_Timing(TSSTX,   'TOX PROC')
  call Report_Max_Min_Timing(TLRPD,   'PARTTRAC')
  call Report_Max_Min_Timing(TWQKIN,  'WQ KIN  ')
  call Report_Max_Min_Timing(TWQSED,  'WQ DIAGN')
  call Report_Max_Min_Timing(TWQRPEM, 'RPEM    ')

  ! *** COMMUNICATION TIMES
  if( process_id == master_id )then
    open(9,FILE = OUTDIR//'TIME.LOG',POSITION = 'APPEND')
    write(9,'(//,"*** MPI COMMUNICATION TIMES ***",//)')
    close(9)
  endif
  
  call Report_Max_Min_Timing(DSITIMING(1),DSITIME(1))
  call Report_Max_Min_Timing(DSITIMING(2),DSITIME(2))
  call Report_Max_Min_Timing(DSITIMING(3),DSITIME(3))
  call Report_Max_Min_Timing(DSITIMING(4),DSITIME(4))
  call Report_Max_Min_Timing(DSITIMING(5),DSITIME(5))
  call Report_Max_Min_Timing(DSITIMING(6),DSITIME(6))
  call Report_Max_Min_Timing(DSITIMING(7),DSITIME(7))
  call Report_Max_Min_Timing(DSITIMING(8),DSITIME(8))
  call Report_Max_Min_Timing(DSITIMING(9),DSITIME(9))
  call Report_Max_Min_Timing(DSITIMING(10),DSITIME(10))
  call Report_Max_Min_Timing(DSITIMING(11),DSITIME(11))
  

  ! *** End the mpi calculation

#ifdef DEBUGGING
  if( process_id == 0 )then
    Pause
  endif
#endif 

  call MPI_Finalize(ierr)

  END

#ifdef GNU  
  SUBROUTINE STOPP(MSG)
#else
  SUBROUTINE STOPP(MSG, IOPEN)
#endif   
  use GLOBAL
  use Variables_MPI
  
  ! *** STOP WITH A PAUSE TO ALLOW USERS TO SEE MESSAGE IF EFDC LAUNCHED BY EE
  character(LEN = *),  intent(IN) :: MSG
#ifndef GNU  
  integer, OPTIONAL, intent(IN) :: IOPEN

  if( .not. PRESENT(IOPEN) )then
#endif   
    open(mpi_error_unit,FILE = OUTDIR//mpi_error_file,POSITION = 'APPEND')
#ifndef GNU  
  endif
#endif   
  

  if( LEN_TRIM(MSG) < 1 )then
    if( process_id == master_id) PRINT '("EFDCPlus Stopped: ",A)', 'SEE #OUTPUT\log_error_proc_xxx FOR MORE INFORMATION'
    write(mpi_error_unit,'("EFDCPlus Stopped: ",A)')
  else
    if( process_id == master_id) PRINT '("EFDCPlus Stopped: ",A)', MSG
    write(mpi_error_unit,'("EFDCPlus Stopped: ",A)') MSG
  endif

  close(mpi_error_unit)
  
  call MPI_Finalize(ierr)

#ifdef DEBUGGING
  Pause
#endif
  
  ERROR STOP
  END SUBROUTINE STOPP
  ! *** Useful Regular Expressions
  ! ***
  ! *** Adds single space around equals sign:     (?<![= ]) = (?! )    
