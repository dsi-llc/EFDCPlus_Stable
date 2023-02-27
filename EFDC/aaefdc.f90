! ----------------------------------------------------------------------!
!   EFDC+: Multifunctional surface water modeling system
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------!
! Copyright 2021-2022 DSI, LLC
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
!  RELEASE:         EFDCPlus_11.5
!                   Domain Decomposition with MPI
!                   Propeller Wash 
!                   New WQ kinetics with user defined algal groups and zooplankton
!                   SIGMA-Zed (SGZ) Vertical Layering
!  DATE:            2022-10-06
!  BY:              DSI, LLC
!                   EDMONDS, WASHINGTON  98020
!                   USA
!                   WEBSITE:  https://www.eemodelingsystem.com/
!  LEAD DEVELOPER:  PAUL M. CRAIG

PROGRAM EFDC
  !
  ! **  WELCOME TO THE ENVIRONMENTAL FLUID DYNAMICS COMPUTER CODE PLUS (EFDC+)
  ! **
  ! **  ORIGINALLY DEVELOPED BY JOHN M. HAMRICK WHILE AT
  ! **                          VIRGINIA INSTITUTE OF MARINE SCIENCE
  ! **                          SCHOOL OF MARINE SCIENCE, THE COLLEGE OF
  ! **                          WILLIAM AND MARY, GLOUCESTER POINT, VA 23062
  ! **
  ! **  THIS VERSION OF EFDC (EFDC+) HAS BEEN SIGNIFICANTLY ENHANCED AND
  ! **  UPDATED AND IS NOW MAINTAINED BY DSI, LLC, EDMONDS, WA.
  ! **
  ! **  EFDC SOLVES THE 3D REYNOLDS AVERAGED NAVIER-STOKES
  ! **  EQUATIONS (WITH HYDROSTATIC AND BOUSINESSQ APPROXIMATIONS) AND
  ! **  TRANSPORT EQUATIONS FOR TURBULENT INTENSITY, TURBULENT
  ! **  INTENSITY X LENGTH SCALE, SALINITY (OR WATER VAPOR CONTENT),
  ! **  TEMPERATURE, AN INERT TRACER (CALLED DYE), A DYNAMICALLY ACTIVE
  ! **  SUSPENDED SETTLING PARTICLE FIELD (CALLED SEDIMENT).  A FREE
  ! **  SURFACE OR RIGID LID IS PRESENT ON THE VERTICAL BOUNDARY Z=1
  ! **  IN THE SIGMA STRETCHED VERTICAL COORDINATE.  THE HORIZONTAL
  ! **  COORDINATE SYSTEM IS CURVILINEAR AND ORTHOGONAL.
  ! **  THE NUMERICAL SOLUTION SCHEME IS ON A SPATIALLY STAGGERED MAC
  ! **  OR C GRID.
  ! **  SPATIAL SOLUTION OF THE EXTERNAL MODE FOR THE FREE SURFACE
  ! **  ELEVATION OR KINEMATIC PRESSURE UNDER THE RIGID LID IS BY
  ! **  CONJUGATE GRADIENT SOLUTION OF A PSEUDO-HEMHOLTZ EQUATION.
  ! **  THE INTERNAL SOLUTION IS IMPLICIT FOR THE VERTICAL SHEAR OR
  ! **  VELOCITY STRUCTURE.
  ! **  A NUMBER OF OPTIONS ARE AVAILABLE FOR REPRESENTING THE ADVECTIVE
  ! **  TRANSPORT TERMS IN THE MOMENTUM AND SCALAR TRANSPORT EQUATIONS.
  ! **
  ! **  PRIMARY DOCUMENTATION INCLUDES:
  !     HAMRICK, J. M., 1992:  A THREE-DIMENSIONAL ENVIRONMENTAL
  !     FLUID DYNAMICS COMPUTER CODE: THEORETICAL AND COMPUTATIONAL
  !     ASPECTS. THE COLLEGE OF WILLIAM AND MARY, VIRGINIA INSTITUTE
  !     OF MARINE SCIENCE, SPECIAL REPORT 317, 63 PP.
  !
  !     A COMPREHENSIVE EFDC+ THEORY KNOWLEDGE BASE HAS BEEN DEVELOPED BY
  !     DSI AND CAN BE FOUND AT:
  !     https://eemodelingsystem.atlassian.net/wiki/spaces/ETG/overview
  !
  ! **  DSI ASSUMES NO LIABILITY FOR USE OF THIS CODE FOR ENVIRONMENTAL
  ! **  AND ENGINEERING STUDIES.
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
  !    2021-08       Tran D. Kien      Rewrote WQ kinetics to allow unlimited phytoplankton and macrophyte classes, with unlimited zooplankton classes
  !                  Paul M. Craig
  !    2021-12       Paul M. Craig     Added propeller jet efflux momentum to EFDC+ flow field
  !    2022-01       Paul M. Craig     Added cohesive mass erosion classes for propwash induced resuspension

  USE GLOBAL
  USE OMP_LIB
  Use Allocate_Initialize      
  USE IFPORT
  USE XYIJCONV,ONLY:XY2IJ
  USE RESTART_MODULE
  USE EFDCOUT
  USE INFOMOD, ONLY:SKIPCOM,READSTR
  USE DRIFTER, ONLY:DRIFTER_INP,AREA_CENTRD
  USE FIELDS
  USE WATERQUALITY, ONLY:WQ3DINP
  USE SHELLFISHMOD, ONLY:ISFFARM,FREE_SHELLFISH
  USE Cyclone
#ifdef NCOUT
  USE DRIFTER, ONLY:close_nc_lpt
#endif
  ! *** MPI  modules
  USE MPI
  USE Variables_MPI
  USE Variables_MPI_Mapping
  USE Variables_MPI_Write_Out
  USE MPI_All_Reduce
  USE Broadcast_Routines
  Use Communicate_Ghost_Routines
  Use Mod_Map_Write_EE_Binary
  Use Mod_Map_Gather_Sort

  ! *** Propwash
  Use Variables_Propwash
  Use Mod_Read_Propwash 
  Use Mod_Setup_Ships

  IMPLICIT NONE

  CHARACTER*80 TITLE

  REAL(RKD), ALLOCATABLE,DIMENSION(:,:) :: SHOTS
  INTEGER,   ALLOCATABLE,DIMENSION(:,:) :: IFDCH
  INTEGER,   ALLOCATABLE,DIMENSION(:)   :: IDX

  REAL, ALLOCATABLE,DIMENSION(:) :: THICK

  REAL :: TSHIFT, DET, XLNUTME, YLTUTMN, TA1, TA2, FORCSUM, ZERO, MAXTHICK, DEPTHMIN, MAXTHICK_local
  REAL :: CCUE, CCVE, CCUN, CCVN, TMPVAL, TMPCOR, ANG1, ANG2, ANG, DETTMP, DZPC
  REAL :: ANGTMP1, ANGTMP2, DYUP1, DYUM1, DDXDDDY, DDYDDDX, DXVLN, DXVLS, C
  REAL :: TMP, BELMIN, VOLLDRY, WTM, WTMP, DELVOL, C1, ETMP, DZGTMP, FRACK
  REAL :: GRADW, GRADE, GRADS, GRADN

  REAL(RK4) :: CPUTIME(2)
  REAL(RKD) :: T0, T1, T2, T3, T4, T9, DELSNAP, TCPU, TTDS, TWAIT
  REAL(RKD), EXTERNAL :: DSTIME

  INTEGER :: COUNT, NSNAPMAX, NT, iStatus, NS, NTMP, NSHOTS, ISNAP, NFDCHIJ
  INTEGER :: ITMPVAL, IFIRST, ILAST, ISO, JDUMY, IISTMP, LPBTMP, LBELMIN
  INTEGER :: L, K, I, J, IP, IS, M, LL, LP, LF, LT, LN, LS, NX, LW, KM, LE, LG, ND, IYEAR, NP, NC, MD

  INTEGER(IK4) :: IERROR, NWR, IU, JU, LU, ID, JD, LD, NJP
  INTEGER(IK4) :: IRET
  INTEGER(IK8) :: NREST,NN
  INTEGER(1)   :: VERSION

  LOGICAL LDEBUG
  LOGICAL(4) :: RES

  CHARACTER*20  BUFFER
  CHARACTER*120 LINE
  CHARACTER*8   DAY
  CHARACTER*200 STR

  ! *** MPI Variables
  INTEGER(4) :: IERR     !< local MPI error flag
  Integer    :: IIN, JIN
  Integer    :: size_mpi
  Double Precision :: starting_time, ending_time
  
  ! *** When updating the version data also update the EXE date in the project settings
  EFDC_VER = '2022-10-06'
  
  IERR = 0
#ifdef DEBUGGING
  !call sleep(20)
#endif

#ifdef _MPI
  ! *** Initialize MPI communicator
  Call Initialize_MPI

  ! ***  Write to screen how manh processes/domains are being used for each process
  WRITE(*,11) num_Processors, process_id
11 FORMAT( '***********************************************************************************',/, &
           '***  This EFDCPlus run is using:',I6,' MPI domains(s). Preparing domain:',i6,'  ***',/, &
           '***********************************************************************************',/)
#else
  process_id = 0
  master_id  = 0
  num_Processors   = 1
#endif

  ! *** GET START TIME OF ELAPSED TIME COUNTER
  TIME_START = DSTIME(1)
  VERSION = 3

  IF( process_id == master_id )THEN

    COUNT = NARGS()
    NTHREADS = 1

    CALL GETARG(0, STR, iStatus)
    I = INDEX(STR,"\",.TRUE.)
    IF( I == 0 )THEN
      I = INDEX(STR,":")
    ENDIF
    EFDC_EXE = STR(I+1:LEN(STR))
    
    ! *** GET THE COMMAND LINE ARGUMENTS, IF ANY
    IF( COUNT > 1 )THEN
      ! *** ARGUMENT 1
      CALL GETARG(1, BUFFER, iStatus)
      !$  IF( Buffer(1:3) == '-NT' .OR. Buffer(1:3) == '-nt' )THEN
      !$    READ(Buffer(4:10),*) NTHREADS
      !$  ENDIF
      IF( COUNT > 2 )THEN
        ! *** ARGUMENT 2
        CALL GETARG(2, BUFFER, iStatus)
        !$  IF( Buffer(1:3) == '-NT' .OR. Buffer(1:3) == '-nt' )THEN
        !$    READ(Buffer(4:10),*) NTHREADS
        !$  ENDIF
      ENDIF
    ENDIF

    IF( VERSION <= 43 .AND. VERSION >= 41 .AND. NTHREADS > 2 ) NTHREADS = 2

    CALL WELCOME

  ENDIF !*** End calculation on master process

  CALL INITFIELDS()
  
  ! *** Send off the number of threads  read in from the command line so 
  ! that the OMP_SET_NUM_THREADS is set properly on all processes
  Call MPI_Bcast(NTHREADS, 1, MPI_Integer, master_id, MPI_Comm_World, ierr)

  ! *** BEGIN OMP SECTION
  !$OMP PARALLEL
  !$OMP MASTER
  !$    NT = OMP_GET_NUM_THREADS()
  !$OMP END MASTER
  !$OMP END PARALLEL
  !$    IF( NT < NTHREADS ) NTHREADS = NT
  !$    CALL OMP_SET_NUM_THREADS(NTHREADS)
  
  If( process_id == master_id )THEN
    !$    WRITE(*,1) num_Processors, NTHREADS

1   FORMAT( '***********************************************************************************',/, &
            '***  This EFDC+ run will use:',I6,' MPI domains(s) and ',I2,' thread(s) per domain  ***',/, &
            '***********************************************************************************',/)
  ENDIF

#ifdef _WIN
  OUTDIR='#output\'
  MPI_Outdir = 'mpi_logs\'
#else
  OUTDIR='#output/'
  MPI_Outdir = 'mpi_logs/'
#endif

  ! *** Create output file for each processor to record inputted parameters, etc.
  Call MPI_barrier(MPI_Comm_World, ierr)

  unit_efdc_out = process_id + unit_efdc_out
  Write(filename_out, '(A5,I3.3,A4)') 'EFDC_',process_id,'.out'

  ! *** Always create the new file
  Call MPI_barrier(MPI_Comm_World, ierr)
  Open(unit_efdc_out, FILE=OUTDIR//filename_out,STATUS='replace')
  If( process_id == master_id )THEN !***Start calculation on master process
    write(STR, '("*** EFDC+ Filename: ",A,",  Version: ",A10," ***")') TRIM(EFDC_EXE), EFDC_VER

    ! **  OPEN OUTPUT FILES
    OPEN(7,FILE=OUTDIR//'EFDC.OUT', STATUS='REPLACE')
    write(7, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))

    OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT', STATUS='REPLACE', SHARED)
    write(8, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))

    OPEN(9,FILE=OUTDIR//'TIME.LOG', STATUS='REPLACE')
    write(9, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
    CLOSE(9)

    ! *** DELETE THE FOLLOWING
    OPEN(1,FILE=OUTDIR//'DRYWET.LOG', STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'SEDIAG.OUT', STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'CFL.OUT',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'NEGSEDSND.OUT', STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'ERROR.LOG', STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')

    ! *** PAUSE THE SCREEN TO ALLOW THE USER TO REVIEW THE NUMBER OF THREADS
    CALL SLEEPQQ(3000)

  ENDIF !***end calculation on master process

#ifdef _MPI
  call MPI_BCAST(NTHREADS, 1, MPI_Int, master_id, MPI_Comm_World,ierr)
  call MPI_BARRIER(MPI_Comm_World, ierr)
#endif 

  ! **  CALL INPUT SUBROUTINE
  NRESTART = 0
  CALL VARINIT

  ! *** Allocate variables for domain decomposition
  Call Allocate_Domain_Decomp

  CALL INPUT(TITLE)
  
  ! *** Writes out the variables read in by the input on each process
  !Call VerifyInput      ! *** Used to test the input/broadcast process.  Activate for testing.

  If( process_id == master_id )THEN
    ! *** Handle Runtime Flag
    OPEN(1,FILE='0run',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE='0run',STATUS='UNKNOWN')
    CLOSE(1)
  ENDIF !***End calculation on master process

  ! *** SET TIME RELATED PARAMETERS
  ! *** THE PARAMETER NTC=NUMBER OF TIME CYCLES, CONTROLS THE LENGTH OF RUN (NUMBER OF TIME STEPS)
  ! *** STARTING AND ENDING DATES, IN DAYS
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8
  TIMEEND = TIMEDAY + (DBLE(TIDALP)*DBLE(NTC) + DBLE(DT))/86400._8

  ! *** MODEL RUN TIME ACCUMULATORS
  !TCYCLE=0.0
  TLRPD = 0.0
  THDMT = 0.0
  TVDIF = 0.0
  TCGRS = 0.0
  TTSED = 0.0
  TSSTX = 0.0
  TCONG = 0.0
  TSADV = 0.0
  TRELAXV = 0.0
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
  NBAN = 49

  ! ***  NTC:     Number of reference time periods in run
  ! ***  NTSPTC:  Number of time steps per reference time period
  ! ***  TCON:    Conversion multiplier to change tbegin to seconds
  ! ***  TBEGIN:  Time origin of run
  ! ***  TIDALP:  Reference time period in sec (ie 44714.16s or 86400s)
  ! ***  NFLTMT:  Number of sub-periods within a reference period, typically = 1 (Research)
  TPN = REAL(NTSPTC)
  NTS = INT8(NTC)*NTSPTC/NFLTMT    ! *** Total # of time steps for entire simulation (typically NFLTMT = 1)
  NLTS = NTSPTC*NLTC               ! *** # Transition Step to Completely linear
  NTTS = NTSPTC*NTTC               ! *** # Transition Step to Fully Non-linear
  NTTS = NTTS+NLTS                 ! *** Total # of Steps to End of Transition
  SNLT = 0.
  NCTBC = 1
  NPRINT = 1
  NTSVB = NTCVB*NTSPTC             ! *** Variable Bouyancy
  ITRMAX = 0
  ITRMIN = 1000
  ERRMAX = 1E-9
  ERRMIN = 1000.
  NBAL = 1
  NBALE = 1
  NBALO = 1
  NBUD = 1
  NHAR = 1
  NTSPTC2 = 2*NTSPTC/NFLTMT        ! *** Twice # of time steps per reference time period
  NDISP = NTS-NTSPTC+2
  NSHOWR = 0
  NSHOWC = 0
  DO NS=1,NASER
    MTSALAST(NS) = 2
  ENDDO
  DO NS=1,NWSER
    MTSWLAST(NS) = 2
  ENDDO
  DO NS=1,NPSER
    MTSPLAST(NS) = 2
  ENDDO
  DO NS=1,NQSER
    MTSQLAST(NS) = 2
  ENDDO
  DO NS=1,NQWRSR
    MTSWRLAST(NS) = 2
  ENDDO
  DO NS=1,NGWSER
    MTSGWLAST(NS) = 2
  ENDDO
  DO NC=1,8
    DO NN=1,NCSER(NC)
      MTSCLAST(NN,NC) = 2
    ENDDO
  ENDDO
  MSFTLST=1

  ! **  EFDC_EXPLORER PRIMARY OUTPUT FREQUENCY DEFINITION TABLE
  ! **  REMOVED THE OBSOLETE SUBROUTINE SURFPLT
  NSNAPSHOTS=0
  NSHOTS=0

  IF( HFREOUT == 1 )THEN
    IF( process_id == master_id )THEN
      WRITE(*,'(A)')'HIGH FREQUENCY TIME SERIES ARE BEING USED'
      WRITE(*,'(A)')'  READING SUBSET.INP'
      OPEN(1,FILE='subset.inp',action='read')
      CALL SKIPCOM(1,'*',2)
      READ(1,*,IOSTAT=ISO) NSUBSET
      IF( ISO > 0 ) CALL STOPP('SUBSET.INP: READING ERROR!')
    endif
    Call Broadcast_Scalar(NSUBSET, master_id)
    
    Allocate(HFREGRP(NSUBSET))
    Call AllocateDSI( IJHFRE,    NSUBSET, 0)
    Call AllocateDSI( NPNT,      NSUBSET, 0)

    Call AllocateDSI( HFREDAYEN, NSUBSET, 0.)
    Call AllocateDSI( HFREDAYBG, NSUBSET, 0.)
    Call AllocateDSI( HFREDUR,   NSUBSET, 0.)
    Call AllocateDSI( HFREDAY,   NSUBSET, 0.)
    Call AllocateDSI( HFREMIN,   NSUBSET, 0.)

    IF( process_id == master_id )THEN
      DO NS=1,NSUBSET
        CALL SKIPCOM(1,'*',2)
        READ(1,*,IOSTAT=ISO) IS, IJHFRE(IS), HFREDAYBG(IS), HFREDUR(IS), HFREMIN(IS), NPNT(IS)
        IF( ISO > 0 )THEN
          CALL STOPP('SUBSET.INP: READING ERROR!')
        ENDIF

        Call AllocateDSI( HFREGRP(IS).ICEL, NPNT(IS), 0)
        Call AllocateDSI( HFREGRP(IS).JCEL, NPNT(IS), 0)
        Call AllocateDSI( HFREGRP(IS).XCEL, NPNT(IS), 0.)
        Call AllocateDSI( HFREGRP(IS).YCEL, NPNT(IS), 0.)

        DO NP=1,NPNT(IS)
          CALL SKIPCOM(1,'*',2)
          IF( IJHFRE(IS) == 1 )THEN
            READ(1,*,IOSTAT=ISO) HFREGRP(IS).ICEL(NP), HFREGRP(IS).JCEL(NP)
          ELSE
            READ(1,*,IOSTAT=ISO) HFREGRP(IS).XCEL(NP), HFREGRP(IS).YCEL(NP)
          ENDIF
          IF( ISO > 0 )THEN
            CALL STOPP('SUBSET.INP: READING ERROR!')
          ENDIF
        ENDDO
      ENDDO
      CLOSE(1)
    endif    ! *** end calculation on master process
    
    Call Broadcast_Array(IJHFRE,    master_id)
    Call Broadcast_Array(HFREDAYBG, master_id)
    Call Broadcast_Array(HFREDUR,   master_id)
    Call Broadcast_Array(HFREMIN,   master_id)
    Call Broadcast_Array(NPNT,      master_id)

    if( process_id /= master_id )THEN
      DO NS=1,NSUBSET
        Call AllocateDSI( HFREGRP(NS).ICEL, NPNT(NS),   0)
        Call AllocateDSI( HFREGRP(NS).JCEL, NPNT(NS),   0)
        Call AllocateDSI( HFREGRP(NS).XCEL, NPNT(NS), 0.0)
        Call AllocateDSI( HFREGRP(NS).YCEL, NPNT(NS), 0.0)
      ENDDO
    endif
    
    DO NS=1,NSUBSET
      Call Broadcast_Array(HFREGRP(NS).ICEL, master_id)
      Call Broadcast_Array(HFREGRP(NS).JCEL, master_id)
      Call Broadcast_Array(HFREGRP(NS).XCEL, master_id)
      Call Broadcast_Array(HFREGRP(NS).YCEL, master_id)
    ENDDO
    
    ! *** UPDATE LOCATIONS TO LOCAL DOMAIN, REMOVING ANY POINTS NOT IN THE CURRENT PROCESS
    DO NS=1,NSUBSET
      IF( IJHFRE(NS) == 0 )THEN
        ! *** CONVERT IJ TO LOCAL
        CALL XY2IJ(HFREGRP(NS), 1)          ! *** Get I & J from X and Y
      ENDIF
      HFREDAY(NS)   = HFREDAYBG(NS)
      HFREDAYEN(NS) = HFREDAYBG(NS) + HFREDUR(NS)/24.
      
      ! *** REMOVE POINTS THAT ARE NOT IN DOMAIN
      IP = 0
      DO NP = 1,NPNT(NS)
        IP = IP + 1
        IF( IJHFRE(NS) == 0 )THEN
          ! *** ICEL AND JCEL ARE LOCAL DOMAIN INDICIES
          IF( HFREGRP(NS).ICEL(NP) > 1 )THEN
            L = LIJ(HFREGRP(NS).ICEL(NP), HFREGRP(NS).JCEL(NP))       ! *** Point is in the domain
          ELSE
            L = 0                                                     ! *** Point is outside the domain
          ENDIF
        ELSE
          ! *** ICEL AND JCEL ARE GLOBAL DOMAIN INDICIES
          LG = LIJ_Global(HFREGRP(NS).ICEL(NP), HFREGRP(NS).JCEL(NP))
          L = Map2Local(LG).LL
        ENDIF
      
        ! *** Check if valid.  Otherwise skip point
        IF( L > 1 )THEN
          ! *** Point in domain
          HFREGRP(NS).ICEL(NP) = IL(L)
          HFREGRP(NS).JCEL(NP) = JL(L)
          
          IF( NP > IP )THEN
            HFREGRP(NS).ICEL(IP) = HFREGRP(NS).ICEL(NP)
            HFREGRP(NS).JCEL(IP) = HFREGRP(NS).JCEL(NP)
            HFREGRP(NS).XCEL(IP) = HFREGRP(NS).XCEL(NP)
            HFREGRP(NS).YCEL(IP) = HFREGRP(NS).YCEL(NP)
          ENDIF
        ELSE
          IP = IP - 1
        ENDIF
      ENDDO
      NPNT(NS) = IP 
    ENDDO
  ENDIF

  ! ***  TCON:    CONVERSION MULTIPLIER TO CHANGE TBEGIN TO SECONDS
  ! ***  TBEGIN:  TIME ORIGIN OF RUN
  ! ***  TIDALP:  REFERENCE TIME PERIOD IN SEC (IE 44714.16S OR 86400S)
  IF( ISPPH == 0 )THEN
    NCPPH=0
    JSPPH=0
  ELSEIF( ISPPH == 2 )THEN
    NCPPH=NTS-(NTSPTC-(NTSPTC/NPPPH))/NFLTMT
    JSPPH=1
    NSNAPSHOTS=2
    SNAPSHOTS(1) = (TBEGIN*TCON + TIDALP*NTC)/86400. - 0.001
    SNAPSHOTS(2) = SNAPSHOTS(1) + 0.001
  ELSE
    ! *** ISPPH == 1 OR 100
    NCPPH=NTSPTC/NPPPH/NFLTMT
    JSPPH=1

    DELSNAP = DBLE(TIDALP)/DBLE(NPPPH)/DBLE(86400.)
    T1 = DBLE(TBEGIN)*DBLE(TCON)/DBLE(86400.)
    T2 = (DBLE(TBEGIN)*DBLE(TCON) + DBLE(TIDALP)*DBLE(NTC))/DBLE(86400.)  ! *** TIMEDAY AT END OF RUN
    T0 = T1
    T9 = T2 + 1.E-6                                                     ! *** TIMEDAY AT END OF RUN PLUS PRECISION TOLERANCE
    NN = (T2-T1)/DELSNAP+2

    ! *** HIGH FREQUENCY SNAPSHOTS
    IF( ISPPH == 100 )THEN
      WRITE(*,'(A)')'HIGH FREQ SNAPSHOTS USED'
      NS=0
      WRITE(*,'(A)')'  READING SNAPSHOTS.INP'
      OPEN(1,FILE='snapshots.inp',STATUS='UNKNOWN',ERR=999)

      CALL SKIPCOM(1,'*')
      READ(1,*)NSHOTS
      ALLOCATE(SHOTS(NSHOTS+1,3))
      SHOTS = 0.0
      DO NS=1,NSHOTS
        READ(1,*)(SHOTS(NS,K),K=1,3)
        SHOTS(NS,3) = SHOTS(NS,3)/1440._8
        SHOTS(NS,2) = SHOTS(NS,2)/24._8
      ENDDO
999   NSHOTS=NS-1
      IF( NSHOTS < 0 )NSHOTS=0
      CLOSE(1)

      SHOTS(NSHOTS+1,1)=T2

      DO NS=1,NSHOTS
        NN = NN + NINT(SHOTS(NS,2)/SHOTS(NS,3),8) + 1
      ENDDO
    ENDIF

    ! *** Reallocate the SNAPSHOTS array
    NSNAPMAX = NN+2
100 DEALLOCATE(SNAPSHOTS)

    Call AllocateDSI( SNAPSHOTS, NSNAPMAX, 0.0)

    ! *** BUILD THE SNAPSHOT DATES
    ISNAP = 1
    IF( NSHOTS > 0 )THEN
      T4 = SHOTS(ISNAP,1)+SHOTS(ISNAP,2)  ! *** ENDING TIME FOR HIGH FREQ
      DO WHILE (T4 < T0)
        ISNAP=ISNAP+1
        IF( ISNAP > NSHOTS )THEN
          NSHOTS=0
          EXIT
        ENDIF
        T4 = SHOTS(ISNAP,1) + SHOTS(ISNAP,2)
      ENDDO
    ENDIF

    NSNAPSHOTS = 0
    T1 = T0
    DO WHILE ( T1 <= T9 )
      IF( ISNAP <= NSHOTS )THEN
        IF( T1 > SHOTS(ISNAP,1) )THEN
          ! *** ENSURE NO OVERLAPPING PERIODS
          T3=SHOTS(ISNAP,1)
          IF( NSNAPSHOTS >= 1 )THEN
            DO WHILE (T3 < SNAPSHOTS(NSNAPSHOTS))
              T3=T3+SHOTS(ISNAP,3)
            ENDDO
          ENDIF

          T4=SHOTS(ISNAP,1)+SHOTS(ISNAP,2)  ! *** ENDING TIME FOR HI!GH FREQ
          IF( T4 > T0 .AND. T3 < T9 )THEN
            ! *** VALID PERIOD, SO BUILD HF SNAPSHOTS

            ! *** CHECK 1ST AND LAST HF PERIOD
            IF( T4 > T9 )T4=T2
            DO WHILE ( T3 <= T0 )
              T3=T3+SHOTS(ISNAP,3)
            ENDDO

            ! *** BUILD HF SNAPSHOTS
            DO WHILE ( T3 <= T4 )
              NSNAPSHOTS=NSNAPSHOTS+1
              SNAPSHOTS(NSNAPSHOTS)=T3
              T3=T3+SHOTS(ISNAP,3)
            ENDDO
          ENDIF
          ISNAP=ISNAP+1  ! *** INCREMENT TO THE NEXT PERIOD

          ! *** Synch up the regular intervals
          IF( NSNAPSHOTS > 0 )THEN
            DO WHILE (T1 < SNAPSHOTS(NSNAPSHOTS))
              T1=T1+DELSNAP
            ENDDO
          ENDIF

          IF( T1 < SHOTS(ISNAP,1) .AND. T1 >= T0 .AND. T1 <= T9 )THEN
            NSNAPSHOTS=NSNAPSHOTS+1
            IF( NSNAPSHOTS > NSNAPMAX )THEN
              NSNAPMAX = NSNAPMAX + 10
              GOTO 100
            ENDIF
            SNAPSHOTS(NSNAPSHOTS)=T1
          ENDIF
        ELSEIF( T1 >= T0 .AND. T1 <= T9 )THEN
          NSNAPSHOTS=NSNAPSHOTS+1
          IF( NSNAPSHOTS > NSNAPMAX )THEN
            NSNAPMAX = NSNAPMAX + 10
            GOTO 100
          ENDIF
          SNAPSHOTS(NSNAPSHOTS)=T1
        ENDIF
      ELSEIF( T1 >= T0 .AND. T1 <= T9 )THEN
        NSNAPSHOTS=NSNAPSHOTS+1
        IF( NSNAPSHOTS > NSNAPMAX )THEN
          NSNAPMAX = NSNAPMAX + 10
          GOTO 100
        ENDIF
        SNAPSHOTS(NSNAPSHOTS)=T1
      ENDIF
      T1=T1+DELSNAP
    ENDDO
  ENDIF

  IF( ISPPH == 100 ) ISPPH = 1  !DHC: 2013-08-15

  IF( process_id == master_id )THEN
    if(MPI_DEBUG_FLAG == .TRUE. )THEN
      WRITE(*,'(A,I8)')'NSNAPSHOTS=',  NSNAPSHOTS
      WRITE(mpi_log_unit,*)'NSNAPSHOTS=',  NSNAPSHOTS
      DO I=1,NSNAPSHOTS
        WRITE(mpi_log_unit,*)'SNAPSHOT: ',I,SNAPSHOTS(I)
      ENDDO
    endif
  ENDIF !*** End calc on master process

  ! *** ENSURE LAST SNAPSHOT SPANS THE END OF THE MODEL RUN
  SNAPSHOTS(NSNAPSHOTS+1) = T2 + 1.
  NSNAPMAX = NSNAPSHOTS
  NSNAPSHOTS = MIN(NSNAPSHOTS,2)        ! *** NSNAPSHOTS=1 IS THE INITIAL CONDITION.  SET TO FIRST MODEL RESULTS SNAPSHOT

  ! **  THREE-DIMENSIONAL HDF FORMAT GRAPHICS FILES: SUBROUTINE OUT3D
  IF( IS3DO == 1 )THEN
    NC3DO=NTS-(NTSPTC-(NTSPTC/NP3DO))/NFLTMT
  ENDIF

  ! **  SET CONTROLS FOR WRITING TO FILTERED, AVERAGED OR RESIDUAL
  ! **  2D SCALAR CONTOURING AND 2D VELOCITY VECTOR PLOTTING FILES
  ! **  RESIDUAL SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATION
  ! **  CONTOURING IN HORIZONTAL: SUBROUTINE RSALPLTH
  DO N=1,7
    IF( ISRSPH(N) >= 1 ) JSRSPH(N)=1
  ENDDO

  ! **  RESIDUAL VELOCITY VECTOR PLOTTING IN HORIZONTAL PLANES:
  ! **  SUBROUTINE RVELPLTH
  IF( ISRVPH >= 1 ) JSRVPH=1

  ! **  RESIDUAL SURFACE ELEVATION PLOTTING IN HORIZONTAL PLANES:
  ! **  SUBROUTINE RVELPLTH
  IF( ISRPPH >= 1 ) JSRPPH=1

  ! **  RESIDUAL SCALAR FIELD CONTOURING IN VERTICAL
  ! **  PLANES: SUBROUTINE RSALPLTV
  DO N=1,7
    IF( ISRSPV(N) >= 1 ) JSRSPV(N)=1
  ENDDO

  ! **  RESIDUAL NORMAL AND TANGENTIAL VELOCITY CONTOURING AND AND
  ! **  TANGENTIAL VELOCITY VECTOR PLOTTING IN VERTICAL PLANES:
  ! **  SUBROUTINE RVELPLTV
  IF( ISRVPV >= 1 ) JSRVPV=1

  ! **  SET CONTROLS FOR WRITING TO DRIFTER, HARMONIC ANALYSIS,
  ! **  RESIDUAL TRANSPORT, AND BLANCE OUTPUT FILES
  NCPD=1
  JSLSHA=1
  IF( ISLSHA == 1 )THEN
    LSLSHA=0
    NCLSHA=NTS-NTCLSHA*NTSPTC
  ENDIF
  IF( ISRESTR == 1 ) JSRESTR=1
  JSWASP=0
  IF( ISWASP >= 1 ) JSWASP=1
  IF( ISBAL >= 1 )THEN
    JSBAL=1
    JSBALO=1
    JSBALE=1
  ENDIF
  JSSBAL=1

  ! **  SET SOME CONSTANTS
  JSTBXY = 0
  CTURB2 = CTURB**0.667
  CTURB3 = CTURB**0.333
  KS     = KC-1
  IF( KS == 0 ) KS = 1
  DZI = REAL(KC)
  DZ  = 1./DZI
  DZS = DZ*DZ
  DT  = TIDALP*REAL(NFLTMT)/REAL(NTSPTC)
  DTI = 1./DT
  DT2 = 2.*DT
  DTMIN  = DT
  AVCON1 = 2.*(1.-AVCON)*DZI*AVO
  G    = 9.81
  GPO  = G*BSC
  GI   = 1./G
  GID2 = .5*GI
  PI   = 3.1415926535898
  PI2  = 2.*PI
  TCVP = 0.0625*TIDALP/PI
  TIMERST = (TIDALP*ISRESTO)/86400.

  ! **  SET CONSTANTS FOR M2 TIDAL CYCLE HARMONIC ANALYSIS
  IF( ISHTA > 0 )THEN
    AC=0.
    AS=0.
    ACS=0.
    TSHIFT=(TBEGIN*TCON/DT)+REAL(NTC-2)*NTSPTC
    DO N=1,NTSPTC
      TNT     = REAL(N)+TSHIFT
      NP      = NTSPTC+N
      WC(N)   = COS(2.*PI*TNT/TPN)          ! DELME - WC and WS are not properly dimensioned.  WS is used for NTSPTC and NTSPTC2
      WS(N)   = SIN(2.*PI*TNT/TPN)
      WC(NP)  = WC(N)
      WS(NP)  = WS(N)
      AC      = AC + 2.*WC(N)*WC(N)
      AS      = AS + 2.*WS(N)*WS(N)
      ACS     = 0.
      WC2(N)  = COS(4.*PI*TNT/TPN)
      WS2(N)  = SIN(4.*PI*TNT/TPN)
      WC2(NP) = WC2(N)
      WS2(NP) = WC2(N)
      AC2     = AC2 + 2.*WC2(N)*WC2(N)
      AS2     = AS2 + 2.*WS2(N)*WS2(N)
      ACS2    = 0.
    ENDDO
    DET=AC*AS-ACS*ACS
    AS=AS/DET
    AC=AC/DET
    ACS=ACS/DET
    DET=AC2*AS2-ACS2*ACS2
    AS2=AS2/DET
    AC2=AC2/DET
    ACS2=ACS2/DET
  ENDIF

  ! **  SET WEIGHTS FOR SALINITY AND TEMPERATURE BOUNDARY INTERPOLATION
  IF( KC > 1 )THEN
    DO K=1,KC
      WTCI(K,1)=REAL(K-KC)/REAL(1-KC)
      WTCI(K,2)=REAL(K-1)/REAL(KC-1)
    ENDDO
  ELSE
    WTCI(1,1)=0.5
    WTCI(1,2)=0.5
  ENDIF

  ! **  INITIALIZE ARRAYS
  CALL AINIT

  ! **  READ IN XLON AND YLAT OR UTME AND UTMN OF CELL CENTERS OF
  ! **  CURVILINEAR PORTION OF THE  GRID
  ISHELTERVARY = 0

  ierr = 0
  Call MPI_Barrier(MPI_Comm_World, ierr)
  IF( ISCLO == 1 )THEN
    Call AllocateDSI( R2D_Global, LCM_Global, 8, 0.0)
    
    IF( process_id == master_id )THEN
      WRITE(*,'(A)')'READING LXLY.INP'
      OPEN(1,FILE='lxly.inp',STATUS='UNKNOWN')
      
      CALL SKIPCOM(1,'*')

      ! *** Domain Centroid
      Center_X = 0.
      Center_Y = 0.
      
      ! *** READ LXLY
      DO LL=2,LA_Global
        IF( ISCORV == 1 )THEN
          ! *** READ LXLY WITH SPATIALLY VARYING CORIOLIS
          !                                 1       2       3     4     5     6      7       8
          !READ(1,*,ERR=3000) IIN, JIN, XLNUTME, YLTUTMN, CCUE, CCVE, CCUN, CCVN, TMPVAL, TMPCOR   ! *** IIN, JIN for MPI domain decomp
          READ(1,*,ERR=3000) IIN, JIN, (R2D_Global(1,J),J=1,8)
        ELSE
          ! *** READ LXLY WITHOUT SPATIALLY VARYING CORIOLIS
          !READ(1,*,ERR=3000) IIN, JIN, XLNUTME, YLTUTMN, CCUE, CCVE, CCUN, CCVN, TMPVAL           ! *** IIN, JIN for MPI domain decomp
          READ(1,*,ERR=3000) IIN, JIN, (R2D_Global(1,J),J=1,7)
          R2D_Global(1,8) = CF
        ENDIF
        
        LG = LIJ_GLOBAL(IIN,JIN)
        DO J=1,8
          R2D_Global(LG,J) = R2D_Global(1,J)
        ENDDO
        
        Center_X = Center_X + R2D_Global(1,1)
        Center_Y = Center_Y + R2D_Global(1,2)
      ENDDO
      CLOSE(1)
      
      Center_X = Center_X / REAL(LA_Global-1,8)
      Center_Y = Center_Y / REAL(LA_Global-1,8)
      Center_X = INT(Center_X)
      Center_Y = INT(Center_Y)
            
      DO LG=2,LA_GLOBAL
        ANG1 = ATAN2(R2D_Global(LG,5), R2D_Global(LG,3))
        ANG2 = ATAN2(-R2D_Global(LG,4),R2D_Global(LG,6))
        ANG = 0.5*(ANG1 + ANG2)
        IF( SIGN(1.,ANG1) /= SIGN(1.,ANG2) )THEN
          IF( ABS(ANG1) > (1.57) .OR. ABS(ANG2) > (1.57) )THEN
            ! *** HANDLE THE DISCONTINUITY AT THE 180 DEGREES ANGLE
            ANG = ANG + ACOS(-1.0)
          ENDIF
        ENDIF
          
        CUE_Global(LG) = COS(ANG)
        CVE_Global(LG) = -SIN(ANG)
        CUN_Global(LG) = SIN(ANG)
        CVN_Global(LG) = COS(ANG)
      ENDDO
    ENDIF
    
    Call MPI_Barrier(comm_2d, ierr)
    Call Broadcast_Array(R2D_Global, master_id)
    Call Broadcast_Scalar(Center_X,  master_id)
    Call Broadcast_Scalar(Center_Y,  master_id)
    Call Broadcast_Array(CUE_Global, master_id)  
    Call Broadcast_Array(CVE_Global, master_id)  
    Call Broadcast_Array(CUN_Global, master_id)  
    Call Broadcast_Array(CVN_Global, master_id)  
    
    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        DLON(L) = R2D_Global(LG,1)
        DLAT(L) = R2D_Global(LG,2)
          
        CUE(L) = CUE_Global(LG)
        CVE(L) = CVE_Global(LG)
        CUN(L) = CUN_Global(LG)
        CVN(L) = CVN_Global(LG)
        FCORC(L) = R2D_Global(LG,8)
          
        WINDSTKA(L) = R2D_Global(LG,7)
        IF( L > 2 .AND. WINDSTKA(L-1) /= WINDSTKA(L) ) ISHELTERVARY = 1
          
        DETTMP = 1./( CUE(L)*CVN(L) - CUN(L)*CVE(L) )
        IF( DETTMP == 0.0 )THEN
          WRITE(6,6262)
          WRITE(6,6263) IL(L), JL(L)
          CALL STOPP('.')
        ENDIF
      ENDIF
    ENDDO
    DEALLOCATE(R2D_Global)
    
6262 FORMAT('  SINGULAR INVERSE TRANSFORM FROM E,N TO CURV X,Y')
6263 FORMAT('  I,J =',2I10/)
  ENDIF


  ! *** COMPUTE CELL AREAS AND CENTROIDS
  IF( ISWAVE >= 3 .OR. ISPD > 0 .OR. (ISWAVE >= 1 .AND. ISWAVE <= 2 .AND. IFWAVE == 1 .AND. SWANGRP == 0) )THEN
    ! *** COMPUTE CELL AREAS AND CENTROIDS, REQUIRES CORNERS.INP FILE
    IF( .NOT. ALLOCATED(XCOR) ) CALL AREA_CENTRD
  ENDIF

  ! *** READ IN SEDFLUME DATA
  IF( LSEDZLJ )THEN
    CALL SEDIC
  ENDIF

  ! *** SET CURVATURE FLAG
  ISCURVATURE = .FALSE.
  FORCSUM = 0.
  DO L=2,LA
    FORCSUM=FORCSUM + FCORC(L)
  ENDDO
  IF( FORCSUM > 1.0E-6 ) ISCURVATURE=.TRUE.

  FCORC(1)=FCORC(2)
  FCORC(LC)=FCORC(LA)

  GOTO 3002
3000 CALL STOPP('READ ERROR FOR FILE LXLY.INP')
3002 CONTINUE

  ZERO=0.
  If( process_id == master_id )THEN
    IF( DEBUG )THEN
      OPEN(1,FILE=OUTDIR//'LIJMAP.OUT',STATUS='UNKNOWN')
      DO L=2,LA
        WRITE(1,1113)L,IL(L),JL(L),ZERO
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF !***End on master process
1112 FORMAT (2I5,2F12.4,6F12.7)
1113 FORMAT (3I5,F10.2)

  ! **  SET CORNER CELL STRESS CORRECTION
  DO L=2,LA
    FSCORTBCV(L)=0.0
  ENDDO

  IF( ISCORTBC >= 1 )THEN
    DO L=2,LA
      FSCORTBCV(L)=FSCORTBC
    ENDDO
  ENDIF

  IF( ISCORTBC == 2 .AND. DEBUG )THEN
    IF( process_id == master_id )THEN
      WRITE(*,'(A)')'READING CORNERC.INP'
      OPEN(1,FILE='cornerc.inp')
      
      CALL SKIPCOM(1,'*')
      READ(1,*)NTMP

      DO NT=1,NTMP
        READ(1,*)I,J,TMPVAL
        L=LIJ(I,J)
        FSCORTBCV(L)=TMPVAL
      ENDDO
      CLOSE(1)
    endif

    Call Broadcast_Array(FSCORTBCV, master_id)
  ENDIF

  ! **  READ SPATIAL AVERAGING MAP FOR FOOD CHAIN MODEL OUTPUT
  IF( ISFDCH == 1 )THEN
    DO L=1,LC
      MFDCHZ(L)=0
    ENDDO
    WRITE(*,'(A)')'READING FOODCHAIN.INP'
    OPEN(1,FILE='foodchain.inp')

    CALL SKIPCOM(1,'*')

    READ(1,*)NFDCHIJ
    DO LT=1,NFDCHIJ
      READ(1,*)I,J,ITMPVAL
      L=LIJ(I,J)
      MFDCHZ(L)=ITMPVAL
    ENDDO
    CLOSE(1)

  ELSEIF( ISFDCH == 2 )THEN
    ! *** NEW BEDFORD APPLICATION
    Call AllocateDSI( IFDCH, ICM, JCM, 0)

    WRITE(*,'(A)')'READING FOODCHAIN.INP'
    OPEN(1,FILE='foodchain.inp',STATUS='OLD')

    CALL SKIPCOM(1,'*')

    IFIRST=1
    ILAST=IC
    DO J=JC,1,-1
      READ(1,66,IOSTAT=ISO)C,(IFDCH(I,J),I=IFIRST,ILAST)
      IF( ISO > 0 )CALL STOPP('READ ERROR FOR FILE FOODCHAIN.INP')
      WRITE (7,166)JDUMY,(IFDCH(I,J),I=IFIRST,ILAST)
    END DO
    DO L=2,LA
      I=IL(L)
      J=JL(L)
      MFDCHZ(L)=IFDCH(I,J)
      IF( NCORENO(I,J) == 7 )MFDCHZ(L)=0
      IF( MFDCHZ(L) > 1 )THEN
        DO K=1,KB
          STDOCB(L,K)=STDOCB(L,K)*0.10
        END DO
      END IF
    END DO
    CLOSE(1)
66  FORMAT (I3,2X,113I1)
166 FORMAT (1X,I3,2X,113I1)
  ENDIF

  ! **  READ IN COUNTER CLOCKWISE ANGLE FROM EAST SPECIFYING PRINCIPAL FLOOD FLOW DIRECTION
  IF( ISTRAN(4) >= 1 .AND. ISSFLFE >= 1 )THEN
    WRITE(*,'(A)')'READING FLDANG.INP'
    OPEN(1,FILE='fldang.inp',STATUS='UNKNOWN')
    DO LL=2,LA
      READ(1,*,ERR=3130)I,J,ANGTMP1,ANGTMP2
      L=LIJ(I,J)
      ACCWFLD(L,1)=0.0174533*ANGTMP1
      ACCWFLD(L,2)=0.0174533*ANGTMP2
    ENDDO
    CLOSE(1)
  ENDIF

  GOTO 3132

3130 CALL STOPP('READ ERROR FOR FILE FLDANG.INP')
3132 CONTINUE

  !---------------------------------------------------------------------------!
  ! **  SET BOUNDARY CONDITION SWITCHES
  CALL SETBCS

#ifdef _MPI
  Call communicate_ghost_cells(DXU, 'DXU')
  Call communicate_ghost_cells(DYU, 'DYU')
  Call communicate_ghost_cells(DXV, 'DXV')
  Call communicate_ghost_cells(DYV, 'DYV')
  Call communicate_ghost_cells(HMU, 'HMU')
  Call communicate_ghost_cells(HMV, 'HMV')

  Call communicate_ghost_cells(SUBO, 'SUBO')
  Call communicate_ghost_cells(SVBO, 'SVBO')
  Call communicate_ghost_cells(SUB,  'SUB')
  Call communicate_ghost_cells(SVB,  'SVB')
  Call communicate_ghost_cells(SAAX, 'SAAX')
  Call communicate_ghost_cells(SAAY, 'SAAY')
  Call communicate_ghost_cells(SCAX, 'SCAX')
  Call communicate_ghost_cells(SCAY, 'SCAY')
  Call communicate_ghost_cells(SDX,  'SDX')
  Call communicate_ghost_cells(SDY,  'SDY')
  Call communicate_ghost_cells(SWB,  'SWB')
#endif

  !---------------------------------------------------------------------------!

  ! *** ALLOCATE MEMORY FOR VARIABLE TO STORE CONCENTRATIONS AT OPEN BOUNDARIES
  Call AllocateDSI( WQBCCON,  NBCSOP, KCM, NACTIVEWC, 0.0)
  Call AllocateDSI( WQBCCON1, NBCSOP, KCM, NACTIVEWC, 0.0)
  WQBCCON  = 0.0
  WQBCCON1 = 0.0

  ! *** READ THE DRIFTER DATA
  IF( ISPD > 0 ) CALL DRIFTER_INP

  ! **  CALCUATE CURVATURE METRICS (NEW ADDITION)
  DO L=1,LC
    DYDI(L)=0.
    DXDJ(L)=0.
  ENDDO

  ! ** DYDI
  TMPVAL = 0.
  DO L=2,LA
    I=IL(L)
    J=JL(L)

    ! *** KEEP ZERO'S FOR DYDI AT DOMAIN BORDERS
    IF( I < 2 .OR. I > IC-1 ) CYCLE
    IF( J < 2 .OR. J > JC-1 ) CYCLE

    ! *** DSI - CHANGED CELL TYPE 5 TO 8
    IF( IJCT(I-1,J) >= 1 .AND. IJCT(I-1,J) <= 8 )THEN
      IF( IJCT(I+1,J) >= 1 .AND. IJCT(I+1,J) <= 8 )THEN
        DYDI(L) = DYU(LEC(L))-DYU(L)
      ELSE
        DDYDDDX=2.*(DYP(L)-DYP(LWC(L)))/(DXP(L)+DXP(LWC(L)))
        DYUP1=DYP(L)+0.5*DDYDDDX*DXP(L)
        DYDI(L)=DYUP1-DYU(L)
      ENDIF
    ELSE
      IF( IJCT(I+1,J) >= 1 .AND. IJCT(I+1,J) <= 8 )THEN
        DDYDDDX=2.*(DYP(LEC(L))-DYP(L))/(DXP(LEC(L))+DXP(L))
        DYUM1=DYP(L)-0.5*DDYDDDX*DXP(L)
        DYDI(L)=DYU(L)-DYUM1
      ELSE
        DYDI(L)=0.0
      ENDIF
    ENDIF
    IF( DYDI(L) > 1.E-7 )THEN
      TMPVAL=1.
    ENDIF
  ENDDO

  ! ** DXDJ
  DO L=2,LA
    LN=LNC(L)
    LS=LSC(L)
    I=IL(L)
    J=JL(L)

    ! *** KEEP ZERO'S FOR DXDJ AT DOMAIN BORDERS
    IF( I < 2 .OR. I > IC-1 ) CYCLE
    IF( J < 2 .OR. J > JC-1 ) CYCLE

    ! *** DSI - CHANGED CELL TYPE 5 TO 8
    IF( IJCT(I,J-1) >= 1 .AND. IJCT(I,J-1) <= 8 )THEN
      IF( IJCT(I,J+1) >= 1 .AND. IJCT(I,J+1) <= 8 )THEN
        DXDJ(L)=DXV(LN)-DXV(L)
      ELSE
        DDXDDDY=2.*(DXP(L)-DXP(LS))/(DYP(L)+DYP(LS))
        DXVLN=DXP(L)+0.5*DDXDDDY*DYP(L)
        DXDJ(L)=DXVLN-DXV(L)
      ENDIF
    ELSE
      IF( IJCT(I,J+1) >= 1 .AND. IJCT(I,J+1) <= 8 )THEN
        DDXDDDY=2.*(DXP(LN)-DXP(L))/(DYP(LN)+DYP(L))
        DXVLS=DXP(L)-0.5*DDXDDDY*DYP(L)
        DXDJ(L)=DXV(L)-DXVLS
      ELSE
        DXDJ(L)=0.0
      ENDIF
    ENDIF
    IF( DXDJ(L) > 1.E-7 )THEN
      TMPVAL = 1.
    ENDIF
  ENDDO

  ! *** SETUP UP EDGE OF HARD BOTTOM REGIONS TO HANDLE BEDLOAD TRYING TO EXIT THE HARD BOTTOM REGION
  IF( ISBEDMAP > 0 )THEN
    IF( (ISTRAN(7) > 0 .AND. ICALC_BL > 0) .OR. (NSEDFLUME > 0 .AND. ICALC_BL > 0) )THEN
      Call AllocateDSI( IDX, LCM, 0)

      ! *** WEST
      BEDEDGEW.NEDGE = 0
      DO LP=1,LASED
        L = LSED(LP)
        LW = LWC(L)
        IF( LW < LA .AND. SUBO(L) > 0.5 .AND. LBED(LW) )THEN
          BEDEDGEW.NEDGE = BEDEDGEW.NEDGE + 1
          IDX(BEDEDGEW.NEDGE) = L
        ENDIF
      ENDDO
      IF( BEDEDGEW.NEDGE > 0 )THEN
        Call AllocateDSI( BEDEDGEW.LEDGE, BEDEDGEW.NEDGE, 0)
        DO LP=1,BEDEDGEW.NEDGE
          BEDEDGEW.LEDGE(LP) = IDX(LP)
        ENDDO
      ENDIF

      ! *** EAST
      IDX = 0
      BEDEDGEE.NEDGE = 0
      DO LP=1,LASED
        L = LSED(LP)
        LE = LEC(L)
        IF( LE < LA .AND. SUBO(LE) > 0.5 .AND. LBED(LE) )THEN
          BEDEDGEE.NEDGE = BEDEDGEE.NEDGE + 1
          IDX(BEDEDGEE.NEDGE) = L
        ENDIF
      ENDDO
      IF( BEDEDGEE.NEDGE > 0 )THEN
        Call AllocateDSI( BEDEDGEE.LEDGE, BEDEDGEE.NEDGE, 0)
        DO LP=1,BEDEDGEE.NEDGE
          BEDEDGEE.LEDGE(LP) = IDX(LP)
        ENDDO
      ENDIF

      ! *** NORTH
      IDX = 0
      BEDEDGEN.NEDGE = 0
      DO LP=1,LASED
        L = LSED(LP)
        LN = LNC(L)
        IF( LN < LA .AND. SVBO(LN) > 0.5 .AND. LBED(LN) )THEN
          BEDEDGEN.NEDGE = BEDEDGEN.NEDGE + 1
          IDX(BEDEDGEN.NEDGE) = L
        ENDIF
      ENDDO
      IF( BEDEDGEN.NEDGE > 0 )THEN
        Call AllocateDSI( BEDEDGEN.LEDGE, BEDEDGEN.NEDGE, 0)
        DO LP=1,BEDEDGEN.NEDGE
          BEDEDGEN.LEDGE(LP) = IDX(LP)
        ENDDO
      ENDIF

      ! *** SOUTH
      IDX = 0
      BEDEDGES.NEDGE = 0
      DO LP=1,LASED
        L = LSED(LP)
        LS = LSC(L)
        IF( LS < LA .AND. SVBO(L) > 0.5 .AND. LBED(LS) )THEN
          BEDEDGES.NEDGE = BEDEDGES.NEDGE + 1
          IDX(BEDEDGES.NEDGE) = L
        ENDIF
      ENDDO
      IF( BEDEDGES.NEDGE > 0 )THEN
        Call AllocateDSI( BEDEDGES.LEDGE, BEDEDGES.NEDGE, 0)
        DO LP=1,BEDEDGES.NEDGE
          BEDEDGES.LEDGE(LP) = IDX(LP)
        ENDDO
      ENDIF
      DEALLOCATE(IDX)
    ENDIF
  ENDIF

  If( process_id == master_id )THEN
    ! **********************************************************************************
    ! *** EFDC CELL MAP LOGGING
    WRITE(8,'(3A6,4A10,4A14)') 'I','J','L','BELV','HP','DX','DY','GRADW','GRADE','GRADS','GRADN'
    DO L=2,LA
      GRADW = -999.
      GRADE = -999.
      GRADS = -999.
      GRADN = -999.
      IF( SUBO(L) > 0.      ) GRADW = ( BELV(LWC(L))+HP(LWC(L) ) - ( BELV(L)+HP(L) )           )/DXU(L)
      IF( SUBO(LEC(L)) > 0. ) GRADE = ( BELV(L)+HP(L)            - ( BELV(LEC(L))+HP(LEC(L)) ) )/DXU(LEC(L))
      IF( SVBO(L) > 0.      ) GRADS = ( BELV(LSC(L))+HP(LSC(L) ) - ( BELV(L)+HP(L) )           )/DYV(L)
      IF( SVBO(LNC(L)) > 0. ) GRADN = ( BELV(L)+HP(L)            - ( BELV(LNC(L))+HP(LNC(L)) ) )/DYV(LNC(L))
      WRITE(8,'(3I6,4F10.3,4F14.5)') IL(L),JL(L),L,BELV(L),HP(L),DXP(L),DYP(L),GRADW,GRADE,GRADS,GRADN
    ENDDO
  ENDIF !***End calc on master process

  ! **********************************************************************************
  ! *** ACTIVE CELL LISTS
  Call AllocateDSI( LLWET, KCM, -NDM,    0)
  Call AllocateDSI( LKWET, LCM,  KCM, -NDM, 0)
  LLWET=0
  LKWET=0
  IF( ISHDMF > 0 )THEN
    Call AllocateDSI( LHDMF,  LCM,  KCM, .false.)
    Call AllocateDSI( LLHDMF, KCM, -NDM,       0)
    Call AllocateDSI( LKHDMF, LCM,  KCM,    -NDM,  0)

    NHDMF = LA-1
  ENDIF

  ! *** ENSURE SUM OF DZC = 1.0000 TO MACHINE PRECISION
  DZPC=0.
  DO K=1,KC
    DZPC = DZPC+DZCK(K)
  ENDDO
  DZPC = DZPC - 1.0
  DO K=1,KC
    DZCK(K) = DZCK(K) + DZPC/DZI
  ENDDO

  ! ***************************************************************************
  ! *** BEGIN SIGMA-Z VARIABLE INITIALIZATION (SGZ)
  ! ***
  ! *** USE ORIGINAL DEPTH FROM DXDY (HMP) SO KSZ'S ARE CONSISTENT FOR COLD START,
  ! ***   RESTART AND CONTINUATION RUNS
  IF( IGRIDV > 0 )THEN
    DO L=2,LA
      HMP(L) = HMP(L) + SGZHPDELTA
    ENDDO
  ENDIF

  ! *** INITIALIZE BOTTOM LAYERS
  IF( IGRIDV > 0 )THEN
    ! *** READ SGZ BOTTOM ACTIVE LAYER FOR SIGMA-Z CELLS
    IF( IGRIDV == 1 )THEN

      IF( process_id == master_id )THEN
        WRITE(*,'(A)')'READING SGZLAYER.INP'
        OPEN(1,FILE='sgzlayer.inp',STATUS='UNKNOWN')
        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES

        DO LL=2,LA_global
          READ(1,*,END=1000) IIN, JIN, K

          LG = LIJ_Global(IIN,JIN)
          IF( LG > 0 )THEN
            IF( K < 1 )THEN
              WRITE(6,'(A,3I5)')'*** ERROR: BAD KSZ IN SGZLAYER.INP, I,J,KSZ',I,J,K
              K = 1
            ENDIF
            KSZ_Global(LG) = K
          ELSE
              WRITE(6,'(A,3I5)') '*** ERROR: BAD L INDEX IN SGZLAYER.INP, I,J,L',I,J,LG
          ENDIF
        ENDDO

1000    CLOSE(1)
      endif

      Call MPI_Barrier(comm_2d, ierr)
      Call Broadcast_Array(KSZ_Global, master_id)

      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          KSZ(L) = KSZ_Global(LG)
        ENDIF
      ENDDO
    ENDIF

    MAXTHICK_local = -1
    DO L=2,LA
      IF( HMP(L) > MAXTHICK_local ) MAXTHICK_local = HMP(L)
    ENDDO
    
    ! *** Should use MAXTHICK global to compute KSZ - DKT
    Call DSI_All_Reduce(MAXTHICK_local, MAXTHICK, MPI_Max, TTDS, 1, TWAIT)

    ! *** COMPUTE NOMINAL LAYER THICKNESSES
    IF( IGRIDV == 2 )THEN
      KSZ = 1
      Call AllocateDSI( THICK, KCM, 0.0)

      DEPTHMIN = 0.
      DO K=KC,1,-1
        IF( K > (KC-KMINV) )DEPTHMIN = DEPTHMIN + DZCK(K)*MAXTHICK
        THICK(K) = DZCK(K)*MAXTHICK
      ENDDO

      ! *** ASSIGN LAYER THICKNESSES
      DO L=2,LA
        IF( HMP(L) < DEPTHMIN )THEN
          ! *** DEFAULT LAYERING
          DZPC = 0.
          KSZ(L) = KC - KMINV + 1
          DO K=KSZ(L),KC
            DZPC = DZPC + DZCK(K)
          ENDDO
          DO K=KSZ(L),KC
            DZC(L,K)  = DZCK(K)/DZPC    ! *** STANDARD DZC
          ENDDO
        ELSE
          ! *** SGZ LAYERING
          TMP = 0.
          TMPCOR = 0.
          DO K=KC,KSZ(L),-1
            TMP = TMP + THICK(K)
            IF( TMP > HMP(L)+0.00001 )THEN
              ! *** DEFINE REMAINING BOTTOM FRACTION
              TMP = TMP - THICK(K)
              DZC(L,K) = (HMP(L)-TMP)/HMP(L)
              KSZ(L) = K
              IF( DZC(L,K) < DZCK(K)*0.2 .OR. DZC(L,K)*HMP(L) < HDRY )THEN
                DZC(L,K+1) = DZC(L,K+1)+DZC(L,K)
                KSZ(L) = KSZ(L)+1
                DZC(L,K) = 0.0
              ENDIF
              EXIT
            ELSE
              DZC(L,K) = THICK(K)/HMP(L)
            ENDIF
            TMPCOR = TMPCOR+DZC(L,K)
          ENDDO
        ENDIF

        ! *** QC
        TMP = ABS(1.0 - SUM(DZC(L,1:KC)))
        IF( ABS(1.0 - SUM(DZC(L,1:KC))) > 2.E-6 )THEN  ! Changed from 1.E-6 to 2.E-6 for the Lake Tahoe - DKT
          TMP = SUM(DZC(L,1:KC))
          PRINT *,' BAD DZC FOR CELL: ',L,TMP
          CALL STOPP('.')
        ENDIF
      ENDDO
    ENDIF

    IF( process_id == master_id )THEN
      ! *** EXPORT THE LAYERING, AS COMPUTED OR ENTERED.
      OPEN(1,FILE=OUTDIR//'SGZLAYER.OUT')
      WRITE(1,'(A)') 'C ** SGZLAYER.OUT - LIST OF THE BOTTOMMOST ACTIVE LAYER'
      WRITE(1,'(A)') 'C **    IGRIDV =  1 - PROVIDED BY SGZLAYER.INP,  2 - KSZ(L) CALCULATED'
      WRITE(1,'(A)') 'C **  '
      WRITE(1,'(A)') 'C **     I     J   KSZ'
      DO LG=2,LA_Global
        WRITE(1,'(4X,I7,2I6,I5)') LG, IL_GL(LG), JL_GL(LG), KSZ_Global(LG)
      ENDDO
      CLOSE(1)
    ENDIF

  ENDIF

  ! **  SET VERTICAL GRID DEPENDENT ARRAYS AND LAYER MASKS
  LSGZU=.FALSE.
  LSGZV=.FALSE.
  DO L=2,LA
    LW = LWC(L)
    LE = LEC(L)
    LS = LSC(L)
    LN = LNC(L)

    DZPC=0.
    DO K=KSZ(L),KC
      DZPC = DZPC + DZCK(K)
    ENDDO

    DO K=1,KC
      SUB3D(L,K) = SUB(L)
      SVB3D(L,K) = SVB(L)

      ! *** INACTIVE LAYER MASK FOR CURRENT CELL
      IF( K < KSZ(L) )THEN
        DZC(L,K)  = 0.0
        LKSZ(L,K) = .TRUE.
      ELSEIF( IGRIDV > 0 )THEN
        ! *** SIGMA-ZED DZC
        LKSZ(L,K) = .FALSE.
        IF( IGRIDV == 1 ) DZC(L,K) = DZCK(K)/DZPC
      ELSE
        ! *** SIGMA STRETCH DZC
        DZC(L,K)  = DZCK(K)/DZPC
        LKSZ(L,K) = .FALSE.
      ENDIF

      ! *** U FACE
      IF( SUBO(L) > 0. )THEN
        KM=MAX(KSZ(LW), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR U FACE
        IF( K >= KM )THEN
          LSGZU(L,K) = .TRUE.
        ELSE
          SUB3D(L,K)=0.0
        ENDIF
      ENDIF

      ! *** V FACE
      IF( SVBO(L) > 0. )THEN
        KM=MAX(KSZ(LS), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR V FACE
        IF( K >= KM )THEN
          LSGZV(L,K) = .TRUE.
        ELSE
          SVB3D(L,K)=0.0
        ENDIF
      ENDIF

      ! *** SET FACE/LAYER FLAGS
      SUB3D(L,K) = SUB3D(L,K)*SUB(L)
      SVB3D(L,K) = SVB3D(L,K)*SVB(L)
      SUB3DO(L,K) = SUB3D(L,K)
      SVB3DO(L,K) = SVB3D(L,K)
    ENDDO

    IF( SUBO(L) < 0.5 )THEN
      KSZU(L)=KSZ(L)
    ELSE
      KSZU(L)=MAX(KSZ(L),KSZ(LW))
    ENDIF
    IF( SVBO(L) < 0.5 )THEN
      KSZV(L)=KSZ(L)
    ELSE
      KSZV(L)=MAX(KSZ(L),KSZ(LS))
    ENDIF

    DO K=KSZ(L),KC
      IF( DZC(L,K) < 1.E-4 )THEN
        WRITE(*,'(A,3I5,F10.4,I5)') 'ERROR: BAD LAYER THICKNESS FOR L,K,KSZ,HP = ',L,K,KSZ(L),HMP(L),process_id
        WRITE(mpi_log_unit,'(A,3I5,F10.4)') 'ERROR: BAD LAYER THICKNESS FOR L,K,KSZ,HP = ',L,K,KSZ(L),HMP(L)
        CALL STOPP('.')
      ENDIF
      DZIC(L,K) = 1./DZC(L,K)
    ENDDO

    DZIG(L,0) = 0.
    DZIG(L,KC) = 0.
    DO K=KSZ(L),KS
      DZG(L,K)  = 0.5*(DZC(L,K)+DZC(L,K+1))
      DZIG(L,K) = 1./DZG(L,K)
    ENDDO

    CDZKMK(L,KSZ(L)) = 0.
    DO K=KSZ(L)+1,KC
      CDZKMK(L,K) = DZIG(L,K-1)*DZIC(L,K)
    ENDDO
    DO K=KSZ(L),KS
      CDZKK(L,K) = DZIC(L,K)*DZIG(L,K)
      CDZKKP(L,K) = DZIG(L,K)*DZIC(L,K+1)
    ENDDO
    CDZKK(L,KC)=0.

    Z(L,0)=0.
    DO K=KSZ(L),KC
      Z(L,K)  = Z(L,K-1) +     DZC(L,K)   ! *** TOP OF LAYER Z
      ZZ(L,K) = Z(L,K)   - 0.5*DZC(L,K)   ! *** MID LAYER Z

      ! *** WALL PROXIMITY FUNCTION
      IF( IFPROX == 0 ) FPROX(L,K)=0.
      IF( K /= KC )THEN
        IF( IFPROX == 1 ) FPROX(L,K) =  1./(VKC*Z(L,K)*(1.-Z(L,K)))**2
        IF( IFPROX == 2 ) FPROX(L,K) = (1./(VKC*Z(L,K))**2)            +CTE5*(1./(VKC*(1.-Z(L,K)))**2)/(CTE4+0.00001)
      ENDIF
    ENDDO

  ENDDO

#ifdef _MPI
  Call communicate_ghost_cells(DZC)
  Call communicate_ghost_cells(DZG)
  Call communicate_ghost_cells(KSZU, 'KSZU')
  Call communicate_ghost_cells(KSZV, 'KSZV')
  Call communicate_ghost_cells(CDZKMK)
  Call communicate_ghost_cells(CDZKK)
  Call communicate_ghost_cells(CDZKKP)
  Call communicate_3d_0(DZIC)
  Call Communicate_3D_0(DZIG)
  Call Communicate_3D_0(Z)
  Call Communicate_3D_0(ZZ)
#endif

  ! *** SET ACTIVE CELL LIST (WITHOUT WET/DRY CELL CONSIDERATION)
  DO ND=1,NDM
    LF = 2 + (ND-1)*LDM
    LL = MIN(LF+LDM-1,LA)
    DO K=1,KC
      LN=0
      DO L=LF,LL
        IF( K < KSZ(L) )CYCLE
        LN = LN+1
        LKWET(LN,K,ND)=L
      ENDDO
      LLWET(K,ND)=LN
    ENDDO
  ENDDO

  ! *** Limit BC layers to active layers
  DO NWR=1,NQWR
    IU = WITH_RET(NWR).IQWRU
    JU = WITH_RET(NWR).JQWRU
    LU = LIJ(IU,JU)
    WITH_RET(NWR).KQWRU = MAX(WITH_RET(NWR).KQWRU, KSZ(LU))
    
    ID = WITH_RET(NWR).IQWRD
    JD = WITH_RET(NWR).JQWRD
    LD = LIJ(ID,JD)
    IF( LD > 1 )THEN
      WITH_RET(NWR).KQWRD = MAX(WITH_RET(NWR).KQWRD, KSZ(LD))  
    ENDIF
  ENDDO
  DO NJP=1,NQJPIJ
    LU = LIJ(IQJP(NJP),JQJP(NJP))
    KUPCJP(NJP) = MAX(KUPCJP(NJP), KSZ(LU))
  ENDDO
  
  If( process_id == master_id )THEN
    IF( IGRIDV > 0 )THEN
      IF( DEBUG ) WRITE(6,'(A)') 'TOTAL NUMBER OF ACTIVE COMPUTATIONAL CELLS PER LAYER'
      WRITE(8,'(//,A)') 'TOTAL NUMBER OF ACTIVE COMPUTATIONAL CELLS PER LAYER'
      DO K=KC,1,-1
        IF( DEBUG ) WRITE(6,8000) '  Layer, LLWET:',K,SUM(LLWET(K,1:NDM))
        WRITE(8,8000) '  Layer, LLWET:',K,SUM(LLWET(K,1:NDM))
      ENDDO
      WRITE(8,8000) ' '
    ENDIF
  ENDIF !***End calc on master process

8000 FORMAT(A,I5,3I10)

  IF( ISHDMF > 0 )THEN
    LLHDMF = LLWET
    LKHDMF = LKWET
  ENDIF

  ! **********************************************************************************

  ! **  READ RESTART CONDITIONS OR INITIALIZE SCALAR FIELDS
  !     ISRESTI == 10 READS AND OLD RESTART FILE GENERATED BY
  !     PRE SEPTEMBER 8, 1992 VERSIONS OF EFDC.FOR
  IF( ISRESTI >= 1 ) CALL Restart_In

  ! *** INITIALIZE SALINITY FIELD IF NOT READ IN FROM RESTART FILE
  IF( ISTRAN(1) >= 1 .AND. (ISRESTI == 0 .OR.  &
     (ISRESTI >= 1 .AND. ISCI(1) == 0 )  .OR.  &
     (ISTOPT(1) > 1)) )THEN  ! *** PMC SINGLE LINE - FORCE IC
    IF( ISTOPT(1) >= 1 )THEN
      NREST=0
      DO K=1,KC
        DO L=2,LA
          SAL(L,K)=SALINIT(L,K)
          SAL1(L,K)=SALINIT(L,K)
        ENDDO
      ENDDO

      M = 1
      DO K=1,KC
        DO LL=1,NCBS
          L=LCBS(LL)
          CLOS(LL,K,M)=SALINIT(L,K)
          NLOS(LL,K,M)=0
          IF( NCSERS(LL,1) == 0 )THEN
            SAL(L,K)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
            SAL1(L,K)=SAL(L,K)
          ENDIF
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBW
          L=LCBW(LL)
          CLOW(LL,K,M)=SALINIT(L,K)
          NLOW(LL,K,M)=0
          IF( NCSERW(LL,1) == 0 )THEN
            SAL(L,K)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
            SAL1(L,K)=SAL(L,K)
          ENDIF
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBE
          L=LCBE(LL)
          CLOE(LL,K,M)=SALINIT(L,K)
          NLOE(LL,K,M)=0
          IF( NCSERE(LL,1) == 0 )THEN
            SAL(L,K)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
            SAL1(L,K)=SAL(L,K)
          ENDIF
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBN
          L=LCBN(LL)
          CLON(LL,K,M)=SALINIT(L,K)
          NLON(LL,K,M)=0
          IF( NCSERN(LL,1) == 0 )THEN
            SAL(L,K)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
            SAL1(L,K)=SAL(L,K)
          ENDIF
        ENDDO
      ENDDO

    ENDIF
  ENDIF
9101 FORMAT(I5)
9102 FORMAT(3I5,12F8.2)

  ! **  INITIALIZE TEMP FIELD IF NOT READ IN FROM RESTART FILE
  IF( ISTRAN(2) >= 1 .AND. (ISRESTI == 0                      .OR.  &
    (ISRESTI >= 1 .AND. ISCI(2) == 0 ) .OR.  &
    (ISTOPT(2) > 9)) )THEN  ! *** PMC SINGLE LINE - FORCE IC
    ! *** SPATIALLY VARYING TEMPERATURE FIELD
    NREST=0
    DO K=1,KC
      DO L=2,LA
        TEM(L,K)  = TEMINIT(L,K)
        TEM1(L,K) = TEM(L,K)
      ENDDO
    ENDDO

    M = 2
    DO K=1,KC
      DO LL=1,NCBS
        L=LCBS(LL)
        CLOS(LL,K,M)=TEMINIT(L,K)
        NLOS(LL,K,M)=0
        IF( NCSERS(LL,2) == 0 )THEN
          TEM(L,K)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
          TEM1(L,K)=TEM(L,K)
        ENDIF
      ENDDO
    ENDDO
    DO K=1,KC
      DO LL=1,NCBW
        L=LCBW(LL)
        CLOW(LL,K,M)=TEMINIT(L,K)
        NLOW(LL,K,M)=0
        IF( NCSERW(LL,2) == 0 )THEN
          TEM(L,K)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
          TEM1(L,K)=TEM(L,K)
        ENDIF
      ENDDO
    ENDDO
    DO K=1,KC
      DO LL=1,NCBE
        L=LCBE(LL)
        CLOE(LL,K,M)=TEMINIT(L,K)
        NLOE(LL,K,M)=0
        IF( NCSERE(LL,2) == 0 )THEN
          TEM(L,K)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
          TEM1(L,K)=TEM(L,K)
        ENDIF
      ENDDO
    ENDDO
    DO K=1,KC
      DO LL=1,NCBN
        L=LCBN(LL)
        CLON(LL,K,M)=TEMINIT(L,K)
        NLON(LL,K,M)=0
        IF( NCSERN(LL,2) == 0 )THEN
          TEM(L,K)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
          TEM1(L,K)=TEM(L,K)
        ENDIF
      ENDDO
    ENDDO

  ENDIF

  ! **  INITIALIZE TEMPERATURE BC IF NOT READ IN FROM RESTART FILE
  !     AND CONSTANT INITIAL CONDITION IS USED
  IF( ISRESTI == 0 .AND. ISTRAN(2) >= 1 )THEN
    IF( ISTOPT(2) == 0 )THEN
      ! *** CONSTANT TEMPERATURE FIELD
      M = 2
      DO K=1,KC
        DO LL=1,NCBS
          CLOS(LL,K,M)=TEMO
          NLOS(LL,K,M)=0
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBW
          CLOW(LL,K,M)=TEMO
          NLOW(LL,K,M)=0
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBE
          CLOE(LL,K,M)=TEMO
          NLOE(LL,K,M)=0
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBN
          CLON(LL,K,M)=TEMO
          NLON(LL,K,M)=0
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! *** RESET IC OPTION, IF USED
  IF( ISTOPT(2) > 9) ISTOPT(2) = ISTOPT(2)-10 ! PMC SINGLE LINE
  !
  ! **  INITIALIZE DYE FIELD IF NOT READ IN FROM RESTART FILE
  !
  IF( ISTRAN(3) >= 1 .AND. (ISRESTI == 0                      .OR.  &
    (ISRESTI >= 1 .AND. ISCI(3) == 0 ) .OR.  &
    (ISTOPT(3) > 1)) )THEN  ! *** PMC SINGLE LINE - FORCE IC
    IF( ISTOPT(3) >= 1 )THEN
      ! *** SPATIALLY VARIABLE DYE FIELD
      NREST=0
      DO MD=1,NDYE
        DO K=1,KC
          DO L=2,LA
            DYE(L,K,MD)=DYEINIT(L,K,MD)
            DYE1(L,K,MD)=DYE(L,K,MD)
          ENDDO
        ENDDO

        M = 2 + MD
        DO K=1,KC
          DO LL=1,NCBS
            L=LCBS(LL)
            CLOS(LL,K,M)=DYEINIT(L,K,MD)
            NLOS(LL,K,M)=0
            IF( NCSERS(LL,3) == 0 )THEN
              DYE(L,K,MD)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              DYE1(L,K,MD)=DYE(L,K,MD)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            L=LCBW(LL)
            CLOW(LL,K,M)=DYEINIT(L,K,MD)
            NLOW(LL,K,M)=0
            IF( NCSERW(LL,3) == 0 )THEN
              DYE(L,K,MD)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              DYE1(L,K,MD)=DYE(L,K,MD)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            L=LCBE(LL)
            CLOE(LL,K,M)=DYEINIT(L,K,MD)
            NLOE(LL,K,M)=0
            IF( NCSERE(LL,3) == 0 )THEN
              DYE(L,K,MD)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              DYE1(L,K,MD)=DYE(L,K,MD)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            L=LCBN(LL)
            CLON(LL,K,M)=DYEINIT(L,K,MD)
            NLON(LL,K,M)=0
            IF( NCSERN(LL,3) == 0 )THEN
              DYE(L,K,MD)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              DYE1(L,K,MD)=DYE(L,K,MD)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! **  INITIALIZE DYE BC IF NOT READ IN FROM RESTART FILE
  ! **  AND CONSTANT INITIAL CONDITIONS ARE USED
  IF( (ISRESTI == 0 .AND. ISTRAN(3) >= 1 ) .OR.  &
    (ISRESTI >= 1 .AND. ISCI(3) == 0 ) )THEN     ! *** PMC SINGLE LINE
    IF( ISTOPT(3) == 0 )THEN
      DO MD=1,NDYE
        M = 2 + MD
        DO K=1,KC
          DO LL=1,NCBS
            CLOS(LL,K,M)=0.
            NLOS(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            CLOW(LL,K,M)=0.
            NLOW(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            CLOE(LL,K,M)=0.
            NLOE(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            CLON(LL,K,M)=0.
            NLON(LL,K,M)=0
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! **  INITIALIZE SFL IF( ISRESTI == 0.AND ISTRAN(4) >= 1 )
  IF( ISRESTI == 0 .AND. ISTRAN(4) >= 1 )THEN
    IF( ISTOPT(4) == 11 )THEN
      DO K=1,KC
        DO L=1,LC
          SFL(L,K)=SFLINIT(L,K)
          SFL2(L,K)=SFLINIT(L,K)
        ENDDO
      ENDDO

      M = 3 + NDYM
      DO K=1,KC
        DO LL=1,NCBS
          L=LCBS(LL)
          CLOS(LL,K,M)=SFLINIT(L,K)
          NLOS(LL,K,M)=0
          SFL(L,K)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
          SFL2(L,K)=SFL(L,K)
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBW
          L=LCBW(LL)
          CLOW(LL,K,M)=SFLINIT(L,K)
          NLOW(LL,K,M)=0
          SFL(L,K)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
          SFL2(L,K)=SFL(L,K)
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBE
          L=LCBE(LL)
          CLOE(LL,K,M)=SFLINIT(L,K)
          NLOE(LL,K,M)=0
          SFL(L,K)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
          SFL2(L,K)=SFL(L,K)
        ENDDO
      ENDDO
      DO K=1,KC
        DO LL=1,NCBN
          L=LCBN(LL)
          CLON(LL,K,M)=SFLINIT(L,K)
          NLON(LL,K,M)=0
          SFL(L,K)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
          SFL2(L,K)=SFL(L,K)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! **  INITIALIZE TOX AND BC IF NOT READ IN FROM RESTART FILE
  ! **  AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(5) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(5) > 0 )THEN
    DO NT=1,NTOX
      IF( ITXINT(NT) == 1 .OR. ITXINT(NT) == 3 )THEN
        DO K=1,KC
          DO L=2,LA
            TOX(L,K,NT)=TOXINIT(L,K,NT)
            TOX1(L,K,NT)=TOX(L,K,NT)
          ENDDO
        ENDDO

        M = MSVTOX(NT)
        DO K=1,KC
          DO LL=1,NCBS
            L=LCBS(LL)
            CLOS(LL,K,M)=TOXINIT(L,K,NT)
            NLOS(LL,K,M)=0
            IF( NCSERS(LL,4) == 0 )THEN
              TOX(L,K,NT)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              TOX1(L,K,NT)=TOX(L,K,NT)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            L=LCBW(LL)
            CLOW(LL,K,M)=TOXINIT(L,K,NT)
            NLOW(LL,K,M)=0
            IF( NCSERW(LL,5) == 0 )THEN
              TOX(L,K,NT)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              TOX1(L,K,NT)=TOX(L,K,NT)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            L=LCBE(LL)
            CLOE(LL,K,M)=TOXINIT(L,K,NT)
            NLOE(LL,K,M)=0
            IF( NCSERE(LL,5) == 0 )THEN
              TOX(L,K,NT)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              TOX1(L,K,NT)=TOX(L,K,NT)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            L=LCBN(LL)
            CLON(LL,K,M)=TOXINIT(L,K,NT)
            NLON(LL,K,M)=0
            IF( NCSERN(LL,5) == 0 )THEN
              TOX(L,K,NT)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              TOX1(L,K,NT)=TOX(L,K,NT)
            ENDIF
          ENDDO
        ENDDO

      ENDIF
    ENDDO
  ENDIF

  ! **  INITIALIZE TOX BC IF NOT READ IN FROM RESTART FILE
  ! **  AND CONSTANT INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(5) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(5) > 0 )THEN
    DO NT=1,NTOX
      IF( ITXINT(NT) == 0 .OR. ITXINT(NT) == 2 )THEN
        M = MSVTOX(NT)
        DO K=1,KC
          DO LL=1,NCBS
            CLOS(LL,K,M)=TOXINTW(NT)
            NLOS(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            CLOW(LL,K,M)=TOXINTW(NT)
            NLOW(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            CLOE(LL,K,M)=TOXINTW(NT)
            NLOE(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            CLON(LL,K,M)=TOXINTW(NT)
            NLON(LL,K,M)=0
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! **  INITIALIZE TOX BED IF NOT READ IN FROM RESTART FILE
  ! **  AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(5) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(5) > 0 )THEN
    DO NT=1,NTOX
      IF( ITXINT(NT) == 2 .OR. ITXINT(NT) == 3 )THEN
        DO K=1,KB
          DO L=2,LA
            TOXB(L,K,NT)=TOXBINIT(L,K,NT)
            TOXB1(L,K,NT)=TOXB(L,K,NT)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  ! **  INITIALIZE SED AND BC IF NOT READ IN FROM RESTART FILE
  ! **  AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(6) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(6) > 0 )THEN
    IF( ISEDINT == 1 .OR. ISEDINT == 3 )THEN
      DO NS=1,NSED
        DO K=1,KC
          DO L=2,LA
            SED(L,K,NS)=SEDINIT(L,K,NS)
            SED1(L,K,NS)=SED(L,K,NS)
          ENDDO
        ENDDO

        M = MSVSED(NS)
        DO K=1,KC
          DO LL=1,NCBS
            L=LCBS(LL)
            CLOS(LL,K,M)=SEDINIT(L,K,NS)
            NLOS(LL,K,M)=0
            IF( NCSERS(LL,6) == 0 )THEN
              SED(L,K,NS)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              SED1(L,K,NS)=SED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            L=LCBW(LL)
            CLOW(LL,K,M)=SEDINIT(L,K,NS)
            NLOW(LL,K,M)=0
            IF( NCSERW(LL,6) == 0 )THEN
              SED(L,K,NS)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              SED1(L,K,NS)=SED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            L=LCBE(LL)
            CLOE(LL,K,M)=SEDINIT(L,K,NS)
            NLOE(LL,K,M)=0
            IF( NCSERE(LL,6) == 0 )THEN
              SED(L,K,NS)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              SED1(L,K,NS)=SED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            L=LCBN(LL)
            CLON(LL,K,M)=SEDINIT(L,K,NS)
            NLON(LL,K,M)=0
            IF( NCSERN(LL,6) == 0 )THEN
              SED(L,K,NS)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              SED1(L,K,NS)=SED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO

      ENDDO

    ENDIF
  ENDIF
  IF( ISTRAN(6) > 0 ) DEALLOCATE(SEDINIT)

  ! **  INITIALIZE SED BC IF NOT READ IN FROM RESTART FILE AND
  ! **  CONSTANT INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(6) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(6) > 0 )THEN
    IF( ISEDINT == 0 .OR. ISEDINT == 2 )THEN
      DO NS=1,NSED
        M = MSVSED(NS)
        DO K=1,KC
          DO LL=1,NCBS
            CLOS(LL,K,M)=SEDO(NS)
            NLOS(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            CLOW(LL,K,M)=SEDO(NS)
            NLOW(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            CLOE(LL,K,M)=SEDO(NS)
            NLOE(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            CLON(LL,K,M)=SEDO(NS)
            NLON(LL,K,M)=0
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! **  INITIALIZE SED BED IF NOT READ IN FROM RESTART FILE
  ! **  AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(6) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(6) > 0 )THEN
    IF( ISEDINT == 2 .OR. ISEDINT == 3 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            SEDB(L,K,NS) = SEDBINIT(L,K,NS)
            SEDB1(L,K,NS) = SEDB(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! **  INITIALIZE SND AND BC IF NOT READ IN FROM RESTART FILE
  ! **  AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(7) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(7) > 0 )THEN
    IF( ISEDINT == 1 .OR. ISEDINT == 3 )THEN
      DO NS=1,NSND
        DO K=1,KC
          DO L=2,LA
            SND(L,K,NS)=SNDINIT(L,K,NS)
            SND1(L,K,NS)=SND(L,K,NS)
          ENDDO
        ENDDO

        M = MSVSND(NS)
        DO K=1,KC
          DO LL=1,NCBS
            L=LCBS(LL)
            CLOS(LL,K,M)=SNDINIT(L,K,NS)
            NLOS(LL,K,M)=0
            IF( NCSERS(LL,7) == 0 )THEN
              SND(L,K,NS)=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
              SND1(L,K,NS)=SND(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            L=LCBW(LL)
            CLOW(LL,K,M)=SNDINIT(L,K,NS)
            NLOW(LL,K,M)=0
            IF( NCSERW(LL,7) == 0 )THEN
              SND(L,K,NS)=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
              SND1(L,K,NS)=SND(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            L=LCBE(LL)
            CLOE(LL,K,M)=SNDINIT(L,K,NS)
            NLOE(LL,K,M)=0
            IF( NCSERE(LL,7) == 0 )THEN
              SND(L,K,NS)=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
              SND1(L,K,NS)=SND(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            L=LCBN(LL)
            CLON(LL,K,M)=SNDINIT(L,K,NS)
            NLON(LL,K,M)=0
            IF( NCSERN(LL,7) == 0 )THEN
              SND(L,K,NS)=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
              SND1(L,K,NS)=SND(L,K,NS)
            ENDIF
          ENDDO
        ENDDO

      ENDDO
    ENDIF
  ENDIF
  IF( ISTRAN(7) > 0 ) DEALLOCATE(SNDINIT)

  ! **  INITIALIZE SND BC IF NOT READ IN FROM RESTART FILE AND
  ! **  CONSTANT INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(7) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(7) > 0 )THEN
    IF( ISEDINT == 0 .OR. ISEDINT == 2 )THEN
      DO NX=1,NSND
        NS = NSED + NX
        M = MSVSND(NX)
        DO K=1,KC
          DO LL=1,NCBS
            CLOS(LL,K,M)=SEDO(NS)
            NLOS(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBW
            CLOW(LL,K,M)=SEDO(NS)
            NLOW(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBE
            CLOE(LL,K,M)=SEDO(NS)
            NLOE(LL,K,M)=0
          ENDDO
        ENDDO
        DO K=1,KC
          DO LL=1,NCBN
            CLON(LL,K,M)=SEDO(NS)
            NLON(LL,K,M)=0
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! **  INITIALIZE SND BED IF NOT READ IN FROM RESTART FILE
  ! **  AND VARIABLE INITIAL CONDITIONS ARE USED
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(7) == 0 ) IISTMP=0
  IF( IISTMP == 0 .AND. ISTRAN(7) > 0 )THEN
    IF( ISEDINT == 2 .OR. ISEDINT == 3 )THEN
      DO NX=1,NSND
        DO K=1,KB
          DO L=2,LA
            SNDB(L,K,NX) = SNDBINIT(L,K,NX)
            SNDB1(L,K,NX) = SNDB(L,K,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! *** INITIALIZE SEDIMENT BED
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) CALL BEDINIT
  IF( ISTRAN(6) > 0 ) DEALLOCATE(SEDBINIT)
  IF( ISTRAN(7) > 0 ) DEALLOCATE(SNDBINIT)

  ! *** SET THE WET CELL LIST
  LAWET = 0
  LADRY=0
  DO L=2,LA
    LAWET = LAWET+1
    LWET(LAWET)=L
  ENDDO
  LDMWET = INT(LAWET/NTHREADS)+1
  LDMDRY = 0

  ! **  ACTIVATE DYE TRACER CONTINUITY CHECK
  IF( ISMMC == 1 .AND. ISTRAN(3) > 0 )THEN
    DO MD=1,NDYE
      DO K=1,KC
        DO L=1,LC
          DYE(L,K,MD)=1.
          DYE1(L,K,MD)=1.
        ENDDO
      ENDDO

      M = 2 + MD
      DO K=1,KC
        DO LL=1,NCBS
          CLOS(LL,K,M)=1.
          NLOS(LL,K,M)=0
        ENDDO
        DO LL=1,NCBW
          CLOW(LL,K,M)=1.
          NLOW(LL,K,M)=0
        ENDDO
        DO LL=1,NCBE
          CLOE(LL,K,M)=1.
          NLOE(LL,K,M)=0
        ENDDO
        DO LL=1,NCBN
          CLON(LL,K,M)=1.
          NLON(LL,K,M)=0
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! *** SET 3D CELL FACE CONSTANTS
  FSGZU=1.
  FSGZV=1.
  DO L=2,LA
    LW=LWC(L)
    LS=LSC(L)


    DO K=1,KC

      ! *** U FACE
      IF( SUBO(L) > 0. )THEN
        KM=MAX(KSZ(LW), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR U FACE
        IF( K >= KM )THEN
          IF( IGRIDV > 0 )THEN
            IF( KSZ(LW) > KSZ(L) )THEN
              SGZU(L,K)  = DZC(LW,K)
            ELSE
              SGZU(L,K)  = DZC(L,K)
            ENDIF
          ELSE
            SGZU(L,K)  = MAX(DZC(LW,K),DZC(L,K))
          ENDIF
        ENDIF
      ENDIF

      ! *** V FACE
      IF( SVBO(L) > 0. )THEN
        KM=MAX(KSZ(LS), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR V FACE
        IF( K >= KM )THEN
          IF( IGRIDV > 0 )THEN
            IF( KSZ(LS) > KSZ(L) )THEN
              SGZV(L,K)  = DZC(LS,K)
            ELSE
              SGZV(L,K)  = DZC(L,K)
            ENDIF
          ELSE
            SGZV(L,K)  = MAX(DZC(LS,K),DZC(L,K))
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    IF( IGRIDV == 2 )THEN
      IF( SUBO(L) > 0.5 )THEN
        TMP = SUM(SGZU(L,1:KC))
        DO K=1,KC
          SGZU(L,K) = SGZU(L,K) / TMP
        ENDDO
      ENDIF

      IF( SVBO(L) > 0.5 )THEN
        TMP = SUM(SGZV(L,1:KC))
        DO K=1,KC
          SGZV(L,K) = SGZV(L,K) / TMP
        ENDDO
      ENDIF

      ! *** QC
      IF( SUB(L) > 0. .AND. ABS(1.0 - SUM(SGZU(L,1:KC))) > 1.E-5 )THEN
        PRINT *,'BAD SGZU:  L, KSZ, SUM(SGZU(L,1:KC))',L,KSZ(L),SUM(SGZU(L,1:KC))
        DO K=1,KC
          PRINT *,K,SGZU(L,K),DZC(L,K)
        ENDDO
        CALL STOPP('BAD SGZU')
      ENDIF
      IF( SVB(L) > 0. .AND. ABS(1.0 - SUM(SGZV(L,1:KC))) > 1.E-5 )THEN
        PRINT *,'BAD SGZV:  L, KSZ, SUM(SGZV(L,1:KC))',L,KSZ(L),SUM(SGZV(L,1:KC))
        DO K=1,KC
          PRINT *,K,SGZU(L,K),DZC(L,K)
        ENDDO
        CALL STOPP('BAD SGZV')
      ENDIF
      ! *** END QC
    ENDIF

    ! *** CELL FACE/INTERFACE CONSTANTS
    DO K=KSZU(L),KS
      DZGU(L,K) = 0.5*(SGZU(L,K)+SGZU(L,K+1))
      IF( SUBO(L) > 0. )THEN
        CDZFU(L,K) = SGZU(L,K)*SGZU(L,K+1)/(SGZU(L,K)+SGZU(L,K+1))
        CDZUU(L,K) = -SGZU(L,  K)/(SGZU(L,K)+SGZU(L,K+1))
        CDZLU(L,K) = -SGZU(L,K+1)/(SGZU(L,K)+SGZU(L,K+1))
      ENDIF
    ENDDO
    DO K=KSZV(L),KS
      DZGV(L,K) = 0.5*(SGZV(L,K)+SGZV(L,K+1))
      IF( SVBO(L) > 0. )THEN
        CDZFV(L,K) = SGZV(L,K)*SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
        CDZUV(L,K) = -SGZV(L,K)  /(SGZV(L,K)+SGZV(L,K+1))
        CDZLV(L,K) = -SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
      ENDIF
    ENDDO

    ! *** U FACE
    IF( SUBO(L) > 0. )THEN
      CDZRU(L,KSZU(L)) = SGZU(L,KSZU(L))-1.
      CDZDU(L,KSZU(L)) = SGZU(L,KSZU(L))
      DO K=KSZU(L)+1,KS
        KM = MAX(KSZ(LW), KSZ(L))      ! *** MINIMUM ACTIVE LAYERS FOR U FACE
        IF( K >= KM )THEN
          CDZRU(L,K) = CDZRU(L,K-1)+SGZU(L,K)
          CDZDU(L,K) = CDZDU(L,K-1)+SGZU(L,K)
        ENDIF
      ENDDO
    ENDIF

    ! *** V FACE
    IF( SVBO(L) > 0. )THEN
      CDZRV(L,KSZV(L)) = SGZV(L,KSZV(L))-1.
      CDZDV(L,KSZV(L)) = SGZV(L,KSZV(L))
      DO K=KSZV(L)+1,KS
        KM = MAX(KSZ(LS), KSZ(L))      ! *** MINIMUM ACTIVE LAYERS FOR V FACE
        IF( K >= KM )THEN
          CDZRV(L,K) = CDZRV(L,K-1)+SGZV(L,K)
          CDZDV(L,K) = CDZDV(L,K-1)+SGZV(L,K)
        ENDIF
      ENDDO
    ENDIF

    DO K=1,KS
      ! *** U FACE
      IF( SUBO(L) > 0. )THEN
        CDZRU(L,K) = CDZRU(L,K)*DZGU(L,K)*CDZLU(L,KSZU(L))
        CDZMU(L,K) = 0.5*SGZU(L,K)*SGZU(L,K+1)
      ENDIF

      ! *** V FACE
      IF( SVBO(L) > 0. )THEN
        CDZRV(L,K) = CDZRV(L,K)*DZGV(L,K)*CDZLV(L,KSZV(L))
        CDZMV(L,K) = 0.5*SGZV(L,K)*SGZV(L,K+1)
      ENDIF
    ENDDO

  ENDDO

  ! *** SET LIST OF CELLS WHOSE ACTIVE LAYER COUNT IS 1
  LDMSGZ1 = 0
  IF( IGRIDV > 0 )THEN
    Call AllocateDSI( LSGZ1, LCM, 0)
    DO L=2,LA
      IF( KSZ(L) == KC )THEN
        LASGZ1 = LASGZ1+1
        LSGZ1(LASGZ1)=L
      ENDIF
    ENDDO
    LDMSGZ1 = INT(LASGZ1/NDM)+1
  ENDIF

  ! ***  TREAT SGZ CELLS WHOSE KSZ=KC (i.e. KMIN=1)
  Call AllocateDSI( LLWETZ, KCM, -NDM,       0)
  Call AllocateDSI( LKWETZ, LCM,  KCM, -NDM, 0)

  DO ND=1,NDM
    DO K=1,KS
      LLWETZ(K,ND) = LLWET(K,ND)
      DO LP=1,LLWET(K,ND)
        LKWETZ(LP,K,ND) = LKWET(LP,K,ND)
      ENDDO
    ENDDO

    LLWETZ(KC,ND) = LLWET(KS,ND)
    DO LP=1,LLWET(KS,ND)
      LKWETZ(LP,KC,ND) = LKWET(LP,KS,ND)
    ENDDO
  ENDDO

  ! *** CALCULATE CONSTANT HORIZONTAL SPATIAL ARRAYS
  DO L=2,LA
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
    HRU(L)   = SUB(L)*HMU(L)*DYU(L)*DXIU(L)
    HRV(L)   = SVB(L)*HMV(L)*DXV(L)*DYIV(L)
    HRUO(L)  = SUBO(L)*DYU(L)*DXIU(L)
    HRVO(L)  = SVBO(L)*DXV(L)*DYIV(L)
    SBX(L)   = 0.5*SUB(L)*DYU(L)
    SBY(L)   = 0.5*SVB(L)*DXV(L)
    SBXO(L)  = 0.5*SUBO(L)*DYU(L)
    SBYO(L)  = 0.5*SVBO(L)*DXV(L)
  ENDDO

  ! *** DETERMINE FSGZU/FSGZV FOR GROSS MOMENTUM
  DO L=2,LA
    LW=LWC(L)
    LS=LSC(L)
    DO K=1,KC
      IF( SGZU(L,K) > 0. )THEN
        FSGZU(L,K) = 1./SGZU(L,K)
      ELSE
        FSGZU(L,K) = 0.
      ENDIF
      IF( SGZV(L,K) > 0. )THEN
        FSGZV(L,K) = 1./SGZV(L,K)
      ELSE
        FSGZV(L,K) = 0.
      ENDIF
    ENDDO
  ENDDO

  ! *** THIRD PASS AT CELL CONSTANTS
  DO L=2,LA
    LW=LWC(L)
    LE=LEC(L)
    LS=LSC(L)
    LN=LNC(L)

    IF( IGRIDV > 0 )THEN
      ! *** CELL INTERFACE METRICS
      FRACK = MIN(0.1*MAXTHICK*DZ,.25)

      BELVW(L) = BELV(L)
      BELVE(L) = BELV(L)
      BELVS(L) = BELV(L)
      BELVN(L) = BELV(L)

      KSZW(L) = KSZ(L)
      KSZE(L) = KSZ(L)
      KSZS(L) = KSZ(L)
      KSZN(L) = KSZ(L)

      DO K=KSZ(L),KC
        SGZW(L,K) = DZC(L,K)
        SGZE(L,K) = DZC(L,K)
        SGZS(L,K) = DZC(L,K)
        SGZN(L,K) = DZC(L,K)

        SGZKW(K,L) = DZC(L,K)
        SGZKE(K,L) = DZC(L,K)
        SGZKS(K,L) = DZC(L,K)
        SGZKN(K,L) = DZC(L,K)

        ! *** TOP OF LAYER
        ZW(L,K) = Z(L,K)
        ZE(L,K) = Z(L,K)
        ZS(L,K) = Z(L,K)
        ZN(L,K) = Z(L,K)

        ! *** MIDDLE OF LAYER
        ZZW(K,L) = ZZ(L,K)
        ZZE(K,L) = ZZ(L,K)
        ZZS(K,L) = ZZ(L,K)
        ZZN(K,L) = ZZ(L,K)
      ENDDO

      IF( SUBO(L) > 0. )THEN
        KSZW(L) = KSZU(L)
        IF( KSZ(LW) > KSZ(L) )THEN
          BELVW(L) = BELV(LW)
          !BELVW(L) = MAX(BELV(LW)-FRACK,BELV(L))
          DO K=1,KC
            SGZW(L,K)  = DZC(LW,K)
            SGZKW(K,L) = DZC(LW,K)
            ZW(L,K)    = Z(LW,K)
            ZZW(K,L)   = ZZ(LW,K)
          ENDDO
        ENDIF
        IF( IGRIDV == 2 .AND. KSZ(LW) == KSZ(L) .AND. KSZ(L) < KC-KMINV+1 )THEN
          IF( BELV(L) >= BELV(LW) )THEN
            LL = L
          ELSE
            LL = LW
          ENDIF
          BELVW(L) = BELV(LL)
          DO K=1,KC
            SGZW(L,K)  = DZC(LL,K)
            SGZKW(K,L) = DZC(LL,K)
            ZW(L,K)    = Z(LL,K)
            ZZW(K,L)   = ZZ(LL,K)
          ENDDO
        ENDIF
      ENDIF

      IF( SUBO(LE) > 0. )THEN
        KSZE(L)=KSZU(LE)
        IF( KSZ(LE) > KSZ(L) )THEN
          BELVE(L) = BELV(LE)
          !BELVE(L) = MAX(BELV(LE)-FRACK,BELV(L))
          DO K=1,KC
            SGZE(L,K)  = DZC(LE,K)
            SGZKE(K,L) = DZC(LE,K)
            ZE(L,K)    = Z(LE,K)
            ZZE(K,L)   = ZZ(LE,K)
          ENDDO
        ELSEIF( IGRIDV == 2 .AND. KSZ(LE) == KSZ(L) .AND. KSZ(L) < KC-KMINV+1 )THEN
          IF( BELV(L) >= BELV(LE) )THEN
            LL = L
          ELSE
            LL = LE
          ENDIF
          BELVE(L) = BELV(LL)
          DO K=1,KC
            SGZE(L,K)  = DZC(LL,K)
            SGZKE(K,L) = DZC(LL,K)
            ZE(L,K)    = Z(LL,K)
            ZZE(K,L)   = ZZ(LL,K)
          ENDDO
        ENDIF
      ENDIF

      IF( SVBO(L) > 0. )THEN
        KSZS(L)=KSZV(L)
        IF( KSZ(LS) > KSZ(L) )THEN
          BELVS(L) = BELV(LS)
          !BELVS(L) = MAX(BELV(LS)-FRACK,BELV(L))
          DO K=1,KC
            SGZS(L,K)  = DZC(LS,K)
            SGZKS(K,L) = DZC(LS,K)
            ZS(L,K)    = Z(LS,K)
            ZZS(K,L)   = ZZ(LS,K)
          ENDDO
        ELSEIF( IGRIDV == 2 .AND. KSZ(LS) == KSZ(L) .AND. KSZ(L) < KC-KMINV+1 )THEN
          IF( BELV(L) >= BELV(LS) )THEN
            LL = L
          ELSE
            LL = LS
          ENDIF
          BELVS(L) = BELV(LL)
          DO K=1,KC
            SGZS(L,K)  = DZC(LL,K)
            SGZKS(K,L) = DZC(LL,K)
            ZS(L,K)    = Z(LL,K)
            ZZS(K,L)   = ZZ(LL,K)
          ENDDO
        ENDIF
      ENDIF

      IF( SVBO(LN) > 0. )THEN
        KSZN(L)=KSZV(LN)
        IF( KSZ(LN) > KSZ(L) )THEN
          BELVN(L) = BELV(LN)
          !BELVN(L) = MAX(BELV(LN)-FRACK,BELV(L))
          DO K=1,KC
            SGZN(L,K)  = DZC(LN,K)
            SGZKN(K,L) = DZC(LN,K)
            ZN(L,K)    = Z(LN,K)
            ZZN(K,L)   = ZZ(LN,K)
          ENDDO
        ELSEIF( IGRIDV == 2 .AND. KSZ(LN) == KSZ(L) .AND. KSZ(L) < KC-KMINV+1 )THEN
          IF( BELV(L) >= BELV(LN) )THEN
            LL = L
          ELSE
            LL = LN
          ENDIF
          BELVN(L) = BELV(LL)
          DO K=1,KC
            SGZN(L,K)  = DZC(LL,K)
            SGZKN(K,L) = DZC(LL,K)
            ZN(L,K)    = Z(LL,K)
            ZZN(K,L)   = ZZ(LL,K)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    DO K=KSZ(L),KS
      IF( SUBO(L) > 0.0 )THEN
        IF( LSGZU(L,K) )THEN
          DZGTMP = 0.5*(SGZU(L,K)+SGZU(L,K+1))
          DZGTMP = 1./DZGTMP
          DZIGSD4U(L,K)=0.25*DZGTMP*DZGTMP
        ELSE
          DZIGSD4U(L,K)=0.0
        ENDIF
      ELSE
        DZIGSD4U(L,K)=0.25*DZIG(L,K)*DZIG(L,K)
      ENDIF

      IF( SVBO(L) > 0.0 )THEN
        IF( LSGZV(L,K) )THEN
          DZGTMP = 0.5*(SGZV(L,K)+SGZV(L,K+1))
          DZGTMP = 1./DZGTMP
          DZIGSD4V(L,K)=0.25*DZGTMP*DZGTMP
        ELSE
          DZIGSD4V(L,K)=0.
        ENDIF
      ELSE
        DZIGSD4V(L,K)=0.25*DZIG(L,K)*DZIG(L,K)
      ENDIF
    ENDDO

  ENDDO

#ifdef _MPI
  Call communicate_ghost_cells(SGZV)
  Call communicate_ghost_cells(SGZU)
  IF( IGRIDV > 0 .OR. ISBLOCKED > 0 )THEN
    Call communicate_ghost_cells(BELVW, 'BELVW')
    Call communicate_ghost_cells(BELVE, 'BELVE')
    Call communicate_ghost_cells(BELVS, 'BELVS')
    Call communicate_ghost_cells(BELVN, 'BELVN')
    Call communicate_ghost_cells(KSZW, 'KSZW')
    Call communicate_ghost_cells(KSZE, 'KSZE')
    Call communicate_ghost_cells(KSZS, 'KSZS')
    Call communicate_ghost_cells(KSZN, 'KSZN')
    Call communicate_ghost_cells(SGZW)
    Call communicate_ghost_cells(SGZE)
    Call communicate_ghost_cells(SGZS)
    Call communicate_ghost_cells(SGZN)
    !Call communicate_ghost_cells(SGZKW)   Reversed indicies.  Uses 0:KCM,LCM
    !Call communicate_ghost_cells(SGZKE)
    !Call communicate_ghost_cells(SGZKS)
    !Call communicate_ghost_cells(SGZKN)
    Call Communicate_3D_0(ZW)
    Call Communicate_3D_0(ZE)
    Call Communicate_3D_0(ZS)
    Call Communicate_3D_0(ZN)
    !Call communicate_ghost_cells(ZZW)   Reversed indicies.  Uses 0:KCM,LCM
    !Call communicate_ghost_cells(ZZE)
    !Call communicate_ghost_cells(ZZS)
    !Call communicate_ghost_cells(ZZN)
  ENDIF
#endif

  IF( IGRIDV > 0 )THEN
    DO L=2,LA
      HMP(L) = HMP(L)-SGZHPDELTA

      ! *** SET FACE DEPTHS
      HPW(L) = HP(L)+BELV(L) - BELVW(L)
      HPE(L) = HP(L)+BELV(L) - BELVE(L)
      HPS(L) = HP(L)+BELV(L) - BELVS(L)
      HPN(L) = HP(L)+BELV(L) - BELVN(L)
    ENDDO
  ENDIF

  ! *** COMPUTE HU/HV FOR INITIAL CONDITIONS
  DO L=2,LA
    LW=LWC(L)
    LS=LSC(L)

    IF( KSZ(LW) > KSZ(L) )THEN
      HU(L) = MAX( 0.5*HPK(L,KSZ(LW)), HP(LW)*(1.+DZC(L,KSZ(LW))*0.1) )
    ELSEIF( KSZ(LW) < KSZ(L) )THEN
      HU(L) = MAX( 0.5*HPK(LW,KSZ(L)), HP(L)*(1.+DZC(LW,KSZ(L))*0.1) )
    ELSE
      HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )/(DXYP(L) + DXYP(LW))
    ENDIF

    IF( KSZ(LS) > KSZ(L) )THEN
      HV(L) = MAX( 0.5*HPK(L,KSZ(LS)), HP(LS)*(1.+DZC(L,KSZ(LS))*0.1) )
    ELSEIF( KSZ(LS) < KSZ(L) )THEN
      HV(L) = MAX( 0.5*HPK(LS,KSZ(L)), HP(L)*(1.+DZC(LS,KSZ(L))*0.1) )
    ELSE
      HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )/(DXYP(L) + DXYP(LS))
    ENDIF

    HPI(L) = 1./HP(L)
    HUI(L) = 1./HU(L)
    HVI(L) = 1./HV(L)

  ENDDO

  ! **  INITIALIZE BUOYANCY AND EQUATION OF STATE - ALL CELLS
  N=0

  ! *** Make sure to check IF density effects are turned on and bouyancy is intialized correctly
  IF( BSC > 1.E-6 )THEN
    CALL CALBUOY(.FALSE.)
#ifdef _MPI
    ! *** Communicate updated B values
    Call communicate_ghost_cells(B)
    IF( ISTRAN(2) > 0 .AND. ISICE == 4 )THEN
      Call communicate_ghost_cells(RHOW)
    ENDIF
#endif
  ENDIF

  ! *** COMPUTATIONAL CELL LIST BY SUB-DOMAIN AND LAYER (WET/DRY CONSIDERED)
  DO ND=1,NDM
    LF=2+(ND-1)*LDM
    LL=MIN(LF+LDM-1,LA)
    DO K=1,KC
      LN=0
      DO L=LF,LL
        IF( LKSZ(L,K) )CYCLE
        IF( LMASKDRY(L) )THEN
          LN = LN+1
          LKWET(LN,K,ND) = L
        ENDIF
      ENDDO
      LLWET(K,ND) = LN
    ENDDO
  ENDDO

  ! *** COMPUTATIONAL CELL LIST FOR ENTIRE DOMAIN, i.e. ND=0, BY LAYER (WET/DRY CONSIDERED)
  DO K=1,KC
    LN=0
    DO L=2,LA
      IF( LKSZ(L,K) )CYCLE
      IF( LMASKDRY(L) )THEN
        LN = LN+1
        LKWET(LN,K,0) = L  ! *** Wet Cell for Layer K
      ENDIF
    ENDDO
    LLWET(K,0) = LN        ! *** Total Wet Cells for Layer K
  ENDDO

  DO K=1,KC
    DO L=2,LA
      HPK(L,K) = HP(L)*DZC(L,K)
    ENDDO
  ENDDO
  IF( ISRESTI > 0 )THEN
    DO K=1,KC
      DO L=2,LA
        H1PK(L,K) = H1P(L)*DZC(L,K)
        H2PK(L,K) = H1PK(L,K)
      ENDDO
    ENDDO
  ENDIF

  TMPVAL = MAX(HMIN,.01)
  HPKI = 1./TMPVAL
  DO K=1,KC
    DO L=2,LA
      IF( K < KSZ(L) )CYCLE
      HPKI(L,K) = 1./HPK(L,K)
    ENDDO
  ENDDO

  IF( ISHDMF > 0 )THEN
    ! *** READ THE HMD SUBSET LIST, IF NEEDED
    IF( IHMDSUB > 0 )THEN
      WRITE(*,'(A)')'READING MAPHMD.INP'
      OPEN(1,FILE='maphmd.inp',STATUS='UNKNOWN')
      LINE=READSTR(1) ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*,IOSTAT=ISO) NHDMF
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPHMD.INP')
      DO NP=1,NHDMF
        READ(1,*,IOSTAT=ISO) L,I,J
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPHMD.INP')
        DO K=KSZ(L),KC
          LHDMF(L,K) = .TRUE.
        ENDDO
      ENDDO
      CLOSE(1)
      PRINT '(A,I10)','  NUMBER OF HMD CELLS FOUND: ',NHDMF
    ELSE
      DO L=2,LA
        DO K=KSZ(L),KC
          LHDMF(L,K) = .TRUE.
        ENDDO
      ENDDO
    ENDIF

    DO ND=1,NDM
      DO K=1,KC
        LN=0
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)
          IF( LHDMF(L,K) )THEN
            LN = LN+1
            LKHDMF(LN,K,ND) = L
          ENDIF
        ENDDO
        LLHDMF(K,ND) = LN     ! *** NUMBER OF WET HDMF CELLS FOR THE CURRENT LAYER
      ENDDO
    ENDDO
  ENDIF

  ! *** END OF SIGMA-Z VARIABLE INITIALIZATION (SGZ)
  ! ***************************************************************************

  ! *** DSI BEGIN SEEPAGE
  IF( ISGWIT == 3 )THEN
    DO L=2,LA
      QGW(L) = QGW(L)*DXYP(L)                ! *** m3/s
    ENDDO
  ENDIF

  ! **  DEACTIVATE DRY CELLS
6902 FORMAT('  DRYING AT N,I,J =',I10,2I6,'  H,H1,H2 =',3(2X,E12.4))

  ! **  INITIALIZE ZERO DIMENSION VOLUME BALANCE
  IF( ISDRY >= 1 .AND. ISDRY <= 98 )THEN

    If( process_id == master_id )THEN
      OPEN(1,FILE=OUTDIR//'ZVOLBAL.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'AVSEL.OUT',STATUS='UNKNOWN')
    ENDIF

    LPBTMP=0
    DO L=2,LA
      TVAR3C(L)=0
      IF( SPB(L) == 0 )THEN
        LPBTMP=LPBTMP+1
        TVAR3C(L)=1
      ENDIF
      LORDER(L)=0
    ENDDO
    TVAR3C(1)=1
    TVAR3C(LC)=1
    LORMAX=LC-2-LPBTMP
    DO LS=1,LORMAX
      BELMIN=100000.
      DO L=2,LA
        IF( SPB(L) /= 0 .AND. TVAR3C(L) /= 1 )THEN
          IF( BELV(L) < BELMIN )THEN
            LBELMIN=L
            BELMIN=BELV(L)
          ENDIF
        ENDIF
      ENDDO
      LORDER(LS)=LBELMIN
      TVAR3C(LBELMIN)=1
    ENDDO

    If( process_id == master_id )THEN
      WRITE(1,5300)
    ENDIF

    LS=1
    L=LORDER(LS)
    BELSURF(LS)=BELV(L)
    ASURFEL(LS)=DXYP(L)
    VOLSEL(LS)=0.

    If( process_id == master_id )THEN
      WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
    ENDIF

    DO LS=2,LORMAX
      L=LORDER(LS)
      BELSURF(LS)=BELV(L)
      ASURFEL(LS)=ASURFEL(LS-1)+DXYP(L)
      VOLSEL(LS)=VOLSEL(LS-1)+0.5*(BELSURF(LS)-BELSURF(LS-1))* &
        (ASURFEL(LS)+ASURFEL(LS-1))
      If( process_id == master_id )THEN
        WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
      ENDIF
    ENDDO
    LS=LORMAX+1
    BELSURF(LS)=BELV(L)+10.0
    ASURFEL(LS)=ASURFEL(LS-1)
    VOLSEL(LS)=VOLSEL(LS-1)+0.5*(BELSURF(LS)-BELSURF(LS-1))* &
      (ASURFEL(LS)+ASURFEL(LS-1))
    If( process_id == master_id )THEN
      WRITE(1,5301)LS,BELSURF(LS),ASURFEL(LS),VOLSEL(LS)
    endif
    VOLZERD=0.
    VOLLDRY=0.
    DO L=2,LA
      IF( SPB(L) /= 0 )THEN
        VOLZERD=VOLZERD+DXYP(L)*HP(L)
        IF( HP(L) > HDRY) VOLLDRY=VOLLDRY+DXYP(L)*HP(L)
      ENDIF
    ENDDO
    DO LS=1,LORMAX
      IF( VOLZERD >= VOLSEL(LS) .AND. VOLZERD < VOLSEL(LS+1) )THEN
        WTM=VOLSEL(LS+1)-VOLZERD
        WTMP=VOLZERD-VOLSEL(LS)
        DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
        WTM=WTM/DELVOL
        WTMP=WTMP/DELVOL
        SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
        ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
      ENDIF
    ENDDO
    VETZERD=VOLZERD

    If( process_id == master_id )THEN
      WRITE(1,5302)
      WRITE(1,5303) SELZERD,ASFZERD,VOLZERD,VOLLDRY
      CLOSE(1)
    ENDIF

  ENDIF
5300 FORMAT('   M    BELSURF     ASURFEL        VOLSEL',/)
5301 FORMAT(1X,I5,2X,F10.5,2X,E12.4,2X,E12.4)
5302 FORMAT(/)
5303 FORMAT(2X,F10.5,3(2X,E12.4))

  ! *** LARGE ASPECT RATIO ASSIGMENTS
  NASPECT = 0
  IF( XYRATIO > 1.1 )THEN
    Call AllocateDSI( LASPECT, LCM, .false.)
    DO L=2,LA
      IF( DXP(L) > XYRATIO*DYP(L) .OR. DYP(L) > XYRATIO*DXP(L) )THEN
        NASPECT = NASPECT+1
        LASPECT(L) = .TRUE.
      ENDIF
    ENDDO
  ENDIF

  ! **  INITIALIZE ELEVATION OF ACTIVE GROUNDWATER ZONE FOR COLD START
  IF( ISGWIE >= 1 .AND. ISRESTI == 0 )THEN
    DO L=2,LA
      IF( HP(L) > HDRY )THEN
        AGWELV(L)=BELV(L)
      ELSE
        IF( BELAGW(L) < SELZERD )THEN
          AGWELV(L)=SELZERD
          AGWELV(L)=MIN(AGWELV(L),BELV(L))
        ELSE
          AGWELV(L)=BELAGW(L)
        ENDIF
      ENDIF
    ENDDO
    DO L=2,LA
      AGWELV1(L)=AGWELV(L)
      AGWELV2(L)=AGWELV(L)
    ENDDO
    OPEN(1,FILE=OUTDIR//'GWELV.OUT',STATUS='UNKNOWN')
    WRITE(1,5400)
    WRITE(1,5402)
    DO L=2,LA
      WRITE(1,5401)IL(L),JL(L),BELV(L),BELAGW(L),AGWELV(L)
    ENDDO
    CLOSE(1)
  ENDIF
5400 FORMAT('   I   J    BELELV      BELAGW     ', &
    '   AGWELV',/)
5401 FORMAT(1X,2I5,2X,F10.5,2X,F10.5,2X,F10.5)
5402 FORMAT(/)

  ! **  CALCULATE CONSTANT C ARRAYS FOR EXTERNAL P SOLUTION
  ! **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV
  ! **  DXYIP=1/(DXP*DYP)
  IF( IRVEC /= 9 )THEN
    DO L=2,LA
      CC(L)=1.
      CCC(L)=1.
    ENDDO
    IF( ISRLID == 1 )THEN
      DO L=2,LA
        CC(L)=0.
        CCC(L)=0.
        IF( SPB(L) == 0. ) CC(L)=1.
        IF( SPB(L) == 0. ) CCC(L)=1.
      ENDDO
    ENDIF
    DO L=2,LA
      LE=LEC(L)
      LN=LNC(L)
      C1=-G*DT*DT*SPB(L)*DXYIP(L)
      CS(L)=C1*HRV(L)
      CW(L)=C1*HRU(L)
      CE(L)=C1*HRU(LE)
      CN(L)=C1*HRV(LN)
      CC(L)=CC(L)-CS(L)-CW(L)-CE(L)-CN(L)
      CCI(L)=1./CC(L)
      CCS(L)=0.25*CS(L)
      CCW(L)=0.25*CW(L)
      CCE(L)=0.25*CE(L)
      CCN(L)=0.25*CN(L)
      CCC(L)=CCC(L)-CCS(L)-CCW(L)-CCE(L)-CCN(L)
      CCCI(L)=1./CCC(L)
    ENDDO
  ENDIF

  ! *** *********************************************************************

  ! *** INITIALIZE WATER QUALITY MODEL AND READ INPUT
  IF( ISTRAN(8) >= 1 ) CALL WQ3DINP

  ! *** Added call to some shell fish initialization
  IF(  ISFFARM > 0 .AND. NSF > 0 )THEN
      DO L=2,LA
        CALL SHELLFISH_REDIST(L) 
        CALL SHELLFISH_LENGTH(L)
      ENDDO
  ENDIF

  ! *** RELEASE ARRAYS NOT USED
  IF( ISWASP /= 17 )THEN
    DEALLOCATE(DYEINIT)
    DEALLOCATE(SALINIT)
    DEALLOCATE(TEMINIT)
    DEALLOCATE(SFLINIT)
  ENDIF

  IF( ISR3DO > 0 .OR. IS3DO > 0 )THEN
    BELVMIN = 1.E32
    SELVMAX = -1.E32
    DO L=2,LA
      IF( ISDRY == 0 .OR. HP(L) >= HDRY ) SELVMAX = MAX(SELVMAX,BELV(L)+HP(L))
      BELVMIN = MIN(BELVMIN,BELV(L))
    ENDDO
    DZPC=(SELVMAX-BELVMIN)/FLOAT(KC)

    ZP(0)=BELVMIN
    DO K=1,KC
      ZP(K) = ZP(K-1)+DZPC
    ENDDO
    DO K=1,KC
      ZZP(K) = 0.5*(ZP(K)+ZP(K-1))
    ENDDO
  ENDIF

#ifdef _MPI
  Call Broadcast_Scalar(ISCURVATURE, master_id )
  Call communicate_ghost_cells(DYDI, 'DYDI')
  Call communicate_ghost_cells(DXDJ, 'DXDJ')
  Call Communicate_3d_Real_Zero_LCM(SUB3DO) ! Handles special case of having a zero as the starting index
  Call Communicate_3d_Real_Zero_LCM(SVB3DO) ! Handles special case of having a zero as the starting index
  Call communicate_ghost_cells(SUB3D)
  Call communicate_ghost_cells(SVB3D)
  !Call communicate_ghost_cells(LKSZ)
  !Call communicate_ghost_cells(LSGZU)
  !Call communicate_ghost_cells(LSGZV)
  Call communicate_ghost_cells(FSGZU)
  Call communicate_ghost_cells(FSGZV)
#endif


  ! *** *********************************************************************
  ! *** Call Propwash initialization if specified
  IF( propwash_on )THEN
    
    ! *** Need to get the cell areas 
    IF( .NOT. ALLOCATED(XCOR) ) CALL AREA_CENTRD
    
#ifdef DEBUGGING
    LDEBUG = .true.
#else
    LDEBUG = .false.
#endif

    ! *** Read the propwash values
    if( process_id == master_id )then
      Call Read_propwash_config(LDEBUG)
    
      Call Read_Ship_Data(LDEBUG)
    
      Call Read_Ship_Tracks(LDEBUG)
    endif
    
    Call Broadcast_Propwash
    
    ! *** Setup all of the ships
    Call Setup_Ships(LDEBUG)
    
  end if
  
  ! *** Initialize propwash linkage, even if not used
  if( istran(6) > 0 .or. istran(7) > 0 )then
    Call AllocateDSI( prop_ero, lcm, -nscm, 0.0)
    Call AllocateDSI( prop_bld, lcm, -nscm, 0.0)
  endif  
  
  ! *** INITIALIZE EFDC EXPLORER OUTPUT (SKIP IF CONTINUATION)
  NITER = 0
  IF( ISPPH == 1 .AND. (ISRESTI == 0 .OR. (ISRESTI /= 0 .AND. ICONTINUE == 0)))THEN
#ifdef _MPI
    ! *** Get local LA
    Call Get_LA_Local_No_Ghost
    
    IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
        if( IHTSTRT > 0 )then
            ! *** do nothing
    
        else! KBT_Global is not read in from a file so we need to create it from the local KBTs
            ! *** Map the local KBTs to a global one
            Call Map_Gather_Sort_Int_1D(size(KBT,1), KBT, size(KBT_Global,1), KBT_Global)
        endif
    End if
  
    Call Map_Write_EE_Binary
#endif

    IF( process_id == master_id )THEN
      CALL EE_LINKAGE(1)
    ENDIF
  ENDIF

  If( process_id == master_id )THEN
    ! *** CLOSE FILES NO LONGER NEEDED TO AVOID DRIVE UNAVAILABLE ERRORS
    CLOSE(7)             ! *** EFDC.OUT
    CLOSE(8)             ! *** EFDCLOG.OUT
  ENDIF
  CLOSE(mpi_log_unit)    ! *** MPI Log File

  LB = 2                 ! *** LB is a globally declared integer that can be hardwired for dubugging a certain cell
  
  ! *** *********************************************************************
  ! *** SELECT FULL HYDRODYNAMIC AND MASS TRANSPORT CALCULATION OR
  ! *** LONG-TERM MASS TRANSPORT CALCULATION  (DISABLED)
  NITERAT = 0

  IF( IS2TIM == 0 ) CALL HDMT
  IF( IS2TIM >= 1 ) CALL HDMT2T
  !IF(ISLTMT >= 1 ) CALL LTMT
  
  ! *** ENSURE FINAL SNAPSHOT HAS BEEN WRITTEN
  IF( ISPPH == 1 .AND. NSNAPSHOTS == NSNAPMAX )THEN
    !CALL BCOUT_ACCUMULATE(DTDYN) !@todo need to modify to owrk with mpi
    Call Map_Write_EE_Binary
    
    If( process_id == master_id )THEN
      CALL EE_LINKAGE(0)
    ENDIF

  ENDIF

  ! *** *********************************************************************
  ! *** WRITE END OF RUN SUMMARIES TO SCREEN AND LOG FILES
  If( process_id == master_id )THEN !***Calc on master process
    CLOSE(8)
    OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')

    ! *** DISPLAY THE MHK SUMMARY
    IF( LMHK )THEN
      ETMP=SUM(ESUP(1:TCOUNT,2:LA))
      WRITE(6,'("SUPPORT ENERGY LOSS    ",F10.4," MW-HR")')ETMP
      ETMP=SUM(EMHK(1:TCOUNT,2:LA))
      WRITE(6,'("MHK ENERGY OUTPUT      ",F10.4," MW-HR")')ETMP
      WRITE(6,'(A)')'See EFDCLOG.OUT File'

      ! *** WRITE TO LOG FILE
      WRITE(8,'(//"***********************************************************************")')
      WRITE(8,'("  MHK SUMMARY")')
      ETMP=SUM(ESUP(1:TCOUNT,2:LA))
      WRITE(8,'("  SUPPORT ENERGY LOSS  ",E14.6," MW-HR")')ETMP
      ETMP=SUM(EMHK(1:TCOUNT,2:LA))
      WRITE(8,'("  MHK ENERGY OUTPUT    ",E14.6," MW-HR")')ETMP
      WRITE(8,'("***********************************************************************"//)')
    ENDIF

    ! *** DISPLAY ANY SEDZLJ WARNINGS
    IF( NWARNING > 0 )THEN
      WRITE(6,'(/A,I10)')   'SEDZLJ WARNING.  DURING THE SIMULATION THE ACTUAL BED SHEAR STRESS EXCEEDED'
      WRITE(6,'(A,I10,A/)') '                 THE RANGE OF DEFINED SHEAR CATAGORIES IN THE SEDFLUME DATA',NWARNING,' TIMES'
      WRITE(8,'(/A,I10)')   'SEDZLJ WARNING.  DURING THE SIMULATION THE ACTUAL BED SHEAR STRESS EXCEEDED'
      WRITE(8,'(A,I10,A/)') '                 THE RANGE OF DEFINED SHEAR CATAGORIES IN THE SEDFLUME DATA',NWARNING,' TIMES'
    ENDIF

    ! *** PRINT THE NUMBER OF THREADS
2   FORMAT( '********************************************************************************',/, &
            '***     This EFDCPlus Run used:',I6,' MPI processor(s) and ',I2,' thread(s)      ***',/, &
            '********************************************************************************',/)
    WRITE(6,2) num_Processors, NTHREADS
    WRITE(8,2) num_Processors, NTHREADS
  endif          ! *** End of Master Thread block

  ! *** OUTPUT TIMING (OUTPUT TO TIME.LOG (UNIT 9) USED BY EFDC_EXPLORER)
  TIME_END = DSTIME(1)
  TIME_END = (TIME_END-TIME_START)/3600.

  ! *** DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS,BUT MAY NOT BE
  ! *** SUPPORTED ON OTHER SYSTEMS.
  TCPU = DTIME(CPUTIME)
  TCPU = TCPU/3600./NTHREADS
  CPUTIME(1) = CPUTIME(1)/3600./NTHREADS
  CPUTIME(2) = CPUTIME(2)/3600./NTHREADS

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
  

  If( process_id == master_id )THEN ! *** Show only on master process
#ifdef NCOUT
     call close_nc_lpt()
#endif
    !------------------------------------------------------------------------
    ! *** WRITE THE TIME SUMMARY TO THE SCREEN
    IF( LSEDZLJ )THEN
      WRITE(6,1995) THDMT    ,TTSED , TSSTX
    ELSE
      WRITE(6,1996) THDMT    ,TTSED , TSSTX
    ENDIF
    WRITE(6,1997) TPUV      ,TCONG
    WRITE(6,1998) TCEXP     ,TAVB
    WRITE(6,1999) TUVW      ,TQQQ
    WRITE(6,2000) TTBXY     ,THEAT
    WRITE(6,2001) TLRPD     ,TSADV
    WRITE(6,2002) THMDF     ,TVDIF
    IF( ISTRAN(8) > 0 )THEN
      WRITE(6,2003) TWQKIN  ,TWQRPEM
      WRITE(6,2004) TWQSED  ,0.
    ENDIF
    IF( ISPROPWASH > 0 )THEN
      WRITE(6,2005) TPROPW  ,0.
    ENDIF
#ifdef _MPI
    WRITE(6,2006) TMPIEE    ,TMPIGH
#endif
    WRITE(6,2007) CPUTIME(1),CPUTIME(2)
    WRITE(6,2008) TIME_END  ,TCPU
    WRITE(6,'(//,20(A8,F14.5,/))') ((DSITIME(I),DSITIMING(I)),I=1,13)

    !------------------------------------------------------------------------
    ! *** WRITE THE TIME SUMMARY TO THE EFDC LOG
    IF( LSEDZLJ )THEN
      WRITE(8,1995)THDMT    ,TTSED , TSSTX
    ELSE
      WRITE(8,1996)THDMT    ,TTSED, TSSTX
    ENDIF
    WRITE(8,1997) TPUV      ,TCONG
    WRITE(8,1998) TCEXP     ,TAVB
    WRITE(8,1999) TUVW      ,TQQQ
    WRITE(8,2000) TTBXY     ,THEAT
    WRITE(8,2001) TLRPD     ,TSADV
    WRITE(8,2002) THMDF     ,TVDIF
    WRITE(8,2003) TWQKIN    ,TWQRPEM
    WRITE(8,2004) TWQSED    ,0.
    WRITE(8,2005) TPROPW    ,0.
    WRITE(8,2006) TMPIEE    ,TMPIGH
    WRITE(8,2007) CPUTIME(1),CPUTIME(2)
    WRITE(8,2008) TIME_END,  TCPU

    !------------------------------------------------------------------------
1995 FORMAT('***TIMING (HOURS)',/, &
      'T HDMT ONLY =',F8.4,'  T SEDZLJ    =',F8.4, F8.4)
1996 FORMAT('***TIMING (HOURS)',/, &
      'T HDMT ONLY =',F8.4,'  T SSEDTOX   =',F8.4, F8.4)
1997 FORMAT('T CALPUV    =',F8.4,'  T CONG GRAD =',F8.4)
1998 FORMAT('T EXPLICIT  =',F8.4,'  T CALC AV   =',F8.4)
1999 FORMAT('T CALC UVW  =',F8.4,'  T TURB QQQ  =',F8.4)
2000 FORMAT('T T&B SHEAR =',F8.4,'  T HEAT PRCS =',F8.4)
2001 FORMAT('T PART TRK  =',F8.4,'  T ADV TRANSP=',F8.4)
2002 FORMAT('T HORIZ DIF =',F8.4,'  T VERT DFUSN=',F8.4)
2003 FORMAT('WQ KINETICS =',F8.4,'  WQ RPEM     =',F8.4)
2004 FORMAT('WQ DIAGEN   =',F8.4,'  NOT USED    =',F8.4)
2005 FORMAT('T PROPWASH  =',F8.4,'  NOT USED    =',F8.4)
2006 FORMAT('MPI EE GATH =',F8.4,'  MPI COMMUNIC=',F8.4)
2007 FORMAT('CPU USER    =',F8.4,'  CPU SYSTEM  =',F8.4)
2008 FORMAT('ELAPSED TIME=',F8.4,'  CPU TIME    =',F8.4)

    !------------------------------------------------------------------------
    ! *** WRITE THE TIME SUMMARY TO THE TIME LOG
    CALL TIMELOG(N,TIMEDAY,OUTDIR,TIME_END*3600._8)

    OPEN(9,FILE=OUTDIR//'TIME.LOG',POSITION='APPEND')
    WRITE(9,2) num_Processors, NTHREADS
    IF( LSEDZLJ )THEN
      WRITE(9,6995)THDMT    ,TTSED, TSSTX
    ELSE
      WRITE(9,6996)THDMT    ,TTSED, TSSTX
    ENDIF
    WRITE(9,6997) TPUV      ,TCONG
    WRITE(9,6998) TCEXP     ,TAVB
    WRITE(9,6999) TUVW      ,TQQQ
    WRITE(9,7000) TTBXY     ,THEAT
    WRITE(9,7001) TLRPD     ,TSADV
    WRITE(9,7002) THMDF     ,TVDIF
    WRITE(9,7003) TWQKIN    ,TWQRPEM
    WRITE(9,7004) TWQSED    ,0.
    WRITE(9,7005) TWQSED    ,0.
    WRITE(9,7006) TMPIEE    ,TMPIGH
    WRITE(9,7007) CPUTIME(1),CPUTIME(2)
    WRITE(9,7008) TIME_END  ,TCPU
    WRITE(9,'(//,20(A8,F14.5,/))') ((DSITIME(I),DSITIMING(I)),I=1,16)

6995 FORMAT('T HDMT ONLY  =',F14.4,'  T SEDZLJ     =',F14.4,F14.4)
6996 FORMAT('T HDMT ONLY  =',F14.4,'  T SSEDTOX    =',F14.4,F14.4)
6997 FORMAT('T CALPUV     =',F14.4,'  T CONG GRAD  =',F14.4)
6998 FORMAT('T EXPLICIT   =',F14.4,'  T CALC AV    =',F14.4)
6999 FORMAT('T CALC UVW   =',F14.4,'  T TURB QQQ   =',F14.4)
7000 FORMAT('T T&B SHEAR  =',F14.4,'  T HEAT PRCS  =',F14.4)
7001 FORMAT('T PART TRK   =',F14.4,'  T ADV TRANSP =',F14.4)
7002 FORMAT('T HORIZ DIFF =',F14.4,'  T VERT DFUSN =',F14.4)
7003 FORMAT('WQ KINETICS  =',F14.4,'  WQ RPEM      =',F14.4)
7004 FORMAT('WQ DIAGEN    =',F14.4,'  NOT USED     =',F14.4)
7005 FORMAT('T PROPWASH   =',F14.4,'  NOT USED     =',F14.4)
7006 FORMAT('MPI EE GATH  =',F14.4,'  MPI COMMUNIC =',F14.4)
7007 FORMAT('CPU USER     =',F14.4,'  CPU SYSTEM   =',F14.4)
7008 FORMAT('ELAPSED TIME =',F14.4,'  CPU TIME     =',F14.4)
    !------------------------------------------------------------------------

    ! *** CLOSE OUTPUT  FILES
    CLOSE(8)
    CLOSE(9)

    ! *** CSV file for timing evaluations
    OPEN(9,FILE=OUTDIR//'TIME.CSV',STATUS='REPLACE')
    WRITE(9,'(A,/)') TITLE
    WRITE(9,'("MPI processors = ,",I10,",, OMP Threads = ,",I10,///)') ,num_Processors, NTHREADS
    IF( LSEDZLJ )THEN
      WRITE(9,8995)THDMT    ,TTSED, TSSTX, DSITIME(8), DSITIMING(8)
    ELSE
      WRITE(9,8996)THDMT    ,TTSED, TSSTX, DSITIME(8), DSITIMING(8)
    ENDIF
    WRITE(9,8997) TPUV      ,TCONG,   DSITIME(1), DSITIMING(1)
    WRITE(9,8998) TCEXP     ,TAVB,    DSITIME(2), DSITIMING(2)
    WRITE(9,8999) TUVW      ,TQQQ,    DSITIME(3), DSITIMING(3)
    WRITE(9,9000) TTBXY     ,THEAT,   DSITIME(4), DSITIMING(4)
    WRITE(9,9001) TLRPD     ,TSADV,   DSITIME(5), DSITIMING(5)
    WRITE(9,9002) THMDF     ,TVDIF,   DSITIME(6), DSITIMING(6)
    WRITE(9,9003) TWQKIN    ,TWQRPEM, DSITIME(7), DSITIMING(7)
    WRITE(9,9004) TWQSED    ,0.,      DSITIME(9), DSITIMING(9)
    WRITE(9,9005) TMPIEE    ,TMPIGH,  DSITIME(10), DSITIMING(10)
    WRITE(9,9006) CPUTIME(1),CPUTIME(2), DSITIME(11), DSITIMING(11)
    WRITE(9,9007) TIME_END  ,TCPU,    DSITIME(12), DSITIMING(12)
    
    CLOSE(9)

8995 FORMAT('T HDMT ONLY  =,',F14.5,',,  T SEDZLJ     =,',F14.5,',',F14.5,',,',A,',',F14.7)
8996 FORMAT('T HDMT ONLY  =,',F14.5,',,  T SSEDTOX    =,',F14.5,',',F14.5,',,',A,',',F14.7)
8997 FORMAT('T CALPUV     =,',F14.5,',,  T CONG GRAD  =,',F14.5,',,,',A,',',F14.7)
8998 FORMAT('T EXPLICIT   =,',F14.5,',,  T CALC AV    =,',F14.5,',,,',A,',',F14.7)
8999 FORMAT('T CALC UVW   =,',F14.5,',,  T TURB QQQ   =,',F14.5,',,,',A,',',F14.7)
9000 FORMAT('T T&B SHEAR  =,',F14.5,',,  T HEAT PRCS  =,',F14.5,',,,',A,',',F14.7)
9001 FORMAT('T PART TRK   =,',F14.5,',,  T ADV TRANSP =,',F14.5,',,,',A,',',F14.7)
9002 FORMAT('T HORIZ DIFF =,',F14.5,',,  T VERT DFUSN =,',F14.5,',,,',A,',',F14.7)
9003 FORMAT('WQ KINETICS  =,',F14.5,',,  WQ RPEM      =,',F14.5,',,,',A,',',F14.7)
9004 FORMAT('WQ DIAGEN    =,',F14.5,',,  NOT USED     =,',F14.5,',,,',A,',',F14.7)
9005 FORMAT('MPI EE GATH  =,',F14.5,',,  MPI COMMUNIC =,',F14.5,',,,',A,',',F14.7)
9006 FORMAT('CPU USER     =,',F14.5,',,  CPU SYSTEM   =,',F14.5,',,,',A,',',F14.7)
9007 FORMAT('ELAPSED TIME =,',F14.5,',,  CPU TIME     =,',F14.5,',,,',A,',',F14.7)
    !------------------------------------------------------------------------
    
    ! *** CLOSE FIELDS AND SHELLFISH
    CALL FREEFIELDS()
    !IF( ISFFARM > 0) CALL FREE_SHELLFISH()

#ifdef WASPOUT
    ! *** Move the WASP Linkage File to the output folder, if needed.
    IF( ISWASP > 0 .OR. ISRCA > 0 .OR. ISICM > 0 )THEN
      INQUIRE(FILE=HYDFIL,EXIST=RES)
      IF( RES )THEN
        ! *** LINKAGE FILE
        IF( IHL_HANDLE /= 0 )THEN
          CALL hlclose(IHL_HANDLE,IERROR)
          IHL_HANDLE = 0
        ENDIF
        TITLE = OUTDIR//'wasp\'//HYDFIL
        INQUIRE(FILE=TITLE,EXIST=RES)
        IF( RES )THEN
          I=DELFILESQQ (TITLE)
          IF( I == 1)RES=RENAMEFILEQQ(HYDFIL,TITLE)
        ELSE
          RES=RENAMEFILEQQ(HYDFIL,TITLE)
        ENDIF

        ! *** LOG/DEBUG FILE
        INQUIRE(FILE=OUTDIR//'wasp\log.txt',EXIST=RES)
        IF( RES )THEN
          I=DELFILESQQ (OUTDIR//'wasp\log.txt')
          IF( I == 1)RES=RENAMEFILEQQ('log.txt',OUTDIR//'wasp\log.txt')
        ELSE
          RES=RENAMEFILEQQ('log.txt',OUTDIR//'wasp\log.txt')
        ENDIF
      ENDIF
    ENDIF
#endif

    ! *** Handle Runtime Flag
    OPEN(1,FILE='0run',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')

  ENDIF !***End calculation on master process

  ! ************************************************************************************************
  ! *** Report timing for every MPI processor
  OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
  
  Call WriteBreak(mpi_log_unit)
  WRITE(mpi_log_unit,'(A,/)') TITLE
  WRITE(mpi_log_unit,'("MPI Processes = ,",I10,",, OMP Threads = ,",I10,/,"MPI Process = ,",I10,//)') ,num_Processors, NTHREADS, process_id
  IF( LSEDZLJ )THEN
    WRITE(mpi_log_unit,6995)THDMT    ,TTSED, TSSTX
  ELSE
    WRITE(mpi_log_unit,6996)THDMT    ,TTSED, TSSTX
  ENDIF
  WRITE(mpi_log_unit,6997) TPUV      ,TCONG
  WRITE(mpi_log_unit,6998) TCEXP     ,TAVB
  WRITE(mpi_log_unit,6999) TUVW      ,TQQQ
  WRITE(mpi_log_unit,7000) TTBXY     ,THEAT
  WRITE(mpi_log_unit,7001) TLRPD     ,TSADV
  WRITE(mpi_log_unit,7002) THMDF     ,TVDIF
  WRITE(mpi_log_unit,7003) TWQKIN    ,TWQRPEM
  WRITE(mpi_log_unit,7004) TWQSED    ,0.
  WRITE(mpi_log_unit,7005) TWQSED    ,0.
  WRITE(mpi_log_unit,7006) TMPIEE    ,TMPIGH
  WRITE(mpi_log_unit,7007) CPUTIME(1),CPUTIME(2)
  WRITE(mpi_log_unit,7008) TIME_END  ,TCPU
  WRITE(mpi_log_unit,'(//,20(A8,F14.5,/))') ((DSITIME(I),DSITIMING(I)),I=1,16)

  CLOSE(mpi_log_unit)
  
  ! *** Do some additional processing of the timing routines
  ! *** Fine the Min/Max to aid in load balancing 
  Call Report_Max_Min_Timing(TPUV,    'CALPUV  ')
  Call Report_Max_Min_Timing(TCEXP,   'EXPLICIT')
  Call Report_Max_Min_Timing(TAVB,    'CALC AV ')
  Call Report_Max_Min_Timing(TQQQ,    'TURB QQQ')
  Call Report_Max_Min_Timing(TTBXY,   'TB SHEAR')
  Call Report_Max_Min_Timing(TSADV,   'ADV TRAN')
  Call Report_Max_Min_Timing(THMDF,   'CALHDMF ')
  Call Report_Max_Min_Timing(TVDIF,   'VERTDIFF')
  Call Report_Max_Min_Timing(THEAT,   'CALHEAT ')
  Call Report_Max_Min_Timing(TTSED,   'SED TRAN')
  Call Report_Max_Min_Timing(TSSTX,   'TOX PROC')
  Call Report_Max_Min_Timing(TLRPD,   'PARTTRAC')
  Call Report_Max_Min_Timing(TWQKIN,  'WQ KIN  ')
  Call Report_Max_Min_Timing(TWQSED,  'WQ DIAGN')
  Call Report_Max_Min_Timing(TWQRPEM, 'RPEM    ')

  ! *** COMMUNICATION TIMES
  if( process_id == master_id )then
    OPEN(9,FILE=OUTDIR//'TIME.LOG',POSITION='APPEND')
    WRITE(9,'(//,"*** MPI COMMUNICATION TIMES ***",//)')
    CLOSE(9)
  endif
  
  Call Report_Max_Min_Timing(DSITIMING(1),DSITIME(1))
  Call Report_Max_Min_Timing(DSITIMING(2),DSITIME(2))
  Call Report_Max_Min_Timing(DSITIMING(3),DSITIME(3))
  Call Report_Max_Min_Timing(DSITIMING(4),DSITIME(4))
  Call Report_Max_Min_Timing(DSITIMING(5),DSITIME(5))
  Call Report_Max_Min_Timing(DSITIMING(6),DSITIME(6))
  Call Report_Max_Min_Timing(DSITIMING(7),DSITIME(7))
  Call Report_Max_Min_Timing(DSITIMING(8),DSITIME(8))
  Call Report_Max_Min_Timing(DSITIMING(9),DSITIME(9))
  Call Report_Max_Min_Timing(DSITIMING(10),DSITIME(10))
  Call Report_Max_Min_Timing(DSITIMING(11),DSITIME(11))
  
  !CLOSE(mpi_log_unit)

  ! *** End the mpi calculation

#ifdef DEBUGGING
  if( process_id == 0 )then
    Pause
  endif
#endif 

  Call MPI_Finalize(ierr)

  END

  SUBROUTINE STOPP(MSG)
  ! *** STOP WITH A PAUSE TO ALLOW USERS TO SEE MESSAGE IF EFDC LAUNCHED BY EE
  CHARACTER(*),INTENT(IN) :: MSG

  IF( LEN_TRIM(MSG) < 1 )THEN
    PRINT '("EFDCPlus Stopped: ",A)', 'SEE #OUTPUT\EFDC_xxxx.log FOR MORE INFORMATION'
  ELSE
    PRINT '("EFDCPlus Stopped: ",A)', MSG
  ENDIF

  PAUSE

  Call MPI_Finalize(ierr)
  STOP

  END SUBROUTINE STOPP
  ! *** Useful Regular Expressions
  ! ***
  ! *** Adds single space around equals sign:     (?<![= ])=(?! )    
