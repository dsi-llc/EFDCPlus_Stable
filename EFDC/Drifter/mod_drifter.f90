! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
!> @details THIS MODULE COMPLETELY REPLACES THE PREVIOUS VERSIONS OF 
!! PARTICLE TRACKING IN EFDC. THE CARDS C67 AND C68 IN THE EFDC.INP FILE WERE 
!! LEFT INTACT TO PROVIDE COMPATIBILITY WITH OLDER VERSIONS OF EFDC.
!! 2018-11  MODIFIED BY PAUL M. CRAIG.  MAJOR REWRITE OF 3D LAGRANGIAN 
!! VELOCITY FIELD AND THE CONTAINER FUNCTION
!
!> @author Zander Mausolff (Updated for MPI)
!> @date 2/17/2020
!---------------------------------------------------------------------------!  
MODULE DRIFTER

  use GLOBAL
  use OMP_LIB
  use XYIJCONV
  use INFOMOD,only:SKIPCOM
#ifndef GNU  
  USE IFPORT
#endif
  use MPI
  use Variables_MPI_Drifter
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out
  use Mod_Gather_Soln
  use Mod_Gather_Drifter
  use Mod_Communicate_Drifters
  use Broadcast_Routines

  use netcdf
  use CONVERTWGS84

  implicit none

  real(RKD),save             :: DAYNEXT
  integer,save               :: NPDAY, NPDAY_MAX
  integer,save,allocatable   :: NPLIST(:)

  real(RKD),pointer,PRIVATE  :: ZCTR(:,:)
  real(RKD),save,allocatable :: LA_BEGTI(:)
  real(RKD),save,allocatable :: LA_ENDTI(:)
  real,save,allocatable      :: GRPWS(:)

  integer,save               :: NSTOPPED(10)
  integer,save,allocatable   :: LA_GRP(:)
  integer,save,allocatable   :: EW_Connectors(:)
  integer,save,allocatable   :: NS_Connectors(:)
  integer,save,allocatable   :: BEDFIX(:)           ! OPTION TO FIX A PARTICLE ON BED AFTER DEPOSITION ON BED
                             
  integer,save,allocatable   :: LCTLU(:)
  integer,save,allocatable   :: LCTLD(:)
  integer,save,allocatable   :: LWRU(:)
  integer,save,allocatable   :: LWRD(:)
                             
  logical,save,allocatable   :: ACTIVE(:)
                             
  real(RKD)                  :: DELTD, DIFFH, DIFFV
  real(8)                    :: EETIME              ! USING REAL*8 EXPLICITY FOR EE LINKAGE REQUIREMENTS
  integer(IK4)               :: NGRP
  integer(IK4)               :: ADJVERP             ! OPTION FOR ADJUSTING VERTICAL POSITION IN CASE OF FULL 3D
  integer(IK4), save         :: NTIMES

  integer :: EE_UNIT = 95
  
  real(RKD), external :: DSTIME
  integer(4) :: nc_lpt(0:6), lpt_time_dim, lpt_npd_dim, lpt_time_idx

  contains

  !---------------------------------------------------------------------------!
  !>@details SOLVE DIFFERENTIAL EQS. FOR (X,Y,Z):
  !! DX = U.DELTD+RAN.SQRT(2EH.DELTD)
  !! DY = V.DELTD+RAN.SQRT(2EH.DELTD)
  !! DZ = W.DELTD+RAN.SQRT(2EV.DELTD)
  !! U(L,K),V(L,K),W(L,K),K = 1,KC,L = 2:LA    CURRENT TIME
  !! U1(L,K),V1(L,K),W1(L,K),K = 1,KC,L = 2:LA PREVIOUS TIME
  !! N: TIME STEP
  !---------------------------------------------------------------------------!

SUBROUTINE DRIFTER_CALC

  implicit none

  ! *** Local variables
  integer(IK4) :: NPP, NP, NW, IT, L
  integer(4)   :: VER, HSIZE, BSIZE, ISOIL        ! *** USING INTEGER*4 EXPLICITY FOR EE LINKAGE REQUIREMENTS

  real(RKD) :: KDX1, KDX2, KDX3, KDX4
  real(RKD) :: KDY1, KDY2, KDY3, KDY4
  real(RKD) :: KDZ1, KDZ2, KDZ3, KDZ4
  real(RKD) :: XLA1, YLA1, ZLA1
  real(RKD) :: U2NP, V2NP, W2NP
  real(RKD) :: DAHX1, DAHY1, DAVZ1, DAHX2, DAHY2, DAVZ2
  real(RKD) :: DXYMIN, KWEIGHT, TODAY
  real(RKD) :: TTDS

  real(RKD),save :: TIMENEXT

  logical(IK4) :: BEDGEMOVE, LFORCE
  character*80 :: TITLE,METHOD

  ! *** New variables for MPI
  integer :: ierr, ii
  integer :: in_west, in_east, in_north, in_south
  integer :: send_tag
  integer :: recv_tag
  integer :: status_msg(MPI_Status_Size)
  real(RKD) :: dxymin_global
  
  ! *** End new variables for MPI

  ! *** Decide if a dynamic or constant time step is used
  if( ISDYNSTP == 0 )then ! *** Constant time step
    DELTD = DT
  else ! *** Dynamic time step
    DELTD = DTDYN
  endif
  
  ! *** FIRST CALL AT GLOBAL RELEASE TIME 
  if( LA_FIRSTCALL == 1 )then
    do L = 2,LA
      ZCTR(L,0:KC) = ZZ(L,0:KC)
      ZCTR(L,KC+1) = Z(L,KC)            ! *** Z is the top layer of the interface (dimensionless)
    enddo

    ! *** Find the min cell dimension for the x and y directions.  Used for EE linkage scaling.
    DXYMIN = MIN(MINVAL(DXP_Global(2:LA_Global)), MINVAL(DYP_Global(2:LA_Global))) ! *** DXP cell dimension in x direction @ center
    if( DXYMIN < 0.2 )then
      XYZSCL = 1000
    else
      XYZSCL = 100
    endif

    ! *** MAKE SURE THE FILE IS NEW - this is the binary file read by EE to visualize drifters
    if(process_id == master_id )then ! *** Only make the master process look at the file
      call create_nc_lpt()
      
      open(ULGR,FILE = OUTDIR//'EE_DRIFTER.OUT',STATUS = 'UNKNOWN',FORM = FMT_BINARY)
      close(ULGR,STATUS = 'DELETE')
      open(ULGR,FILE = OUTDIR//'EE_DRIFTER.OUT',ACTION = 'WRITE',FORM = FMT_BINARY)
      ISOIL = 0
      if( ANY(ISOILSPI == 1) ) ISOIL = 1
      VER = 8400
      HSIZE = 6*4
      
      ! *** Write out header information for the drifter binary file
      write(ULGR) VER, HSIZE
      write(ULGR) INT(NPD,4), INT(KC,4), XYZSCL                 ! *** NPD --> is global
      write(ULGR) ISOIL
        
      FLUSH(ULGR)
      close(ULGR,STATUS = 'KEEP')
    endif

    ! *** LA_FREQ - output frequency for the lagrangian calculation
    TIMENEXT = TIMEDAY + LA_FREQ + 0.000001_8
    LA_FIRSTCALL = 0
    NTIMES = 0

    ! *** If we are not reading from a restart file proceed
    if( ISRESTI == 0 )then
      ! *** Call routine to initalize some stuff for the oil spill module
      if( ANY(ISOILSPI == 1) ) CALL INIT_OIL                    ! *** GET: DVOL, DARE

      ! *** Loop over the number of drifters
      do NP = 1,NPD
        ! *** If the cell has not been initalized and the current time is within the lagrangian particle tracking time interval read in.
        if( JSPD(NP) == 1 .and. TIMEDAY >= LA_BEGTI(NP) .and. TIMEDAY <= LA_ENDTI(NP) )then
          ! *** Only process drifters that are active in a given subdomain
          if( LLA_Process(NP) == process_id )then
            call CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)     ! *** Only NP is used for initialization
          endif
          JSPD(NP) = 0                                          ! *** Particle is now initialized
        endif
      enddo
    else ! *** Read from the restart file
      !< @todo Update to read the restart file
      ! *** Warning not working for MPI
      write(*,*) '*** WARNING *** Restart functionality for Lagrangian Particle Tracking module is not available'
      call DRIFTER_RST

      ! *** INITIALIZE BOTTOM ELEVATION AND DEPTH
      do NP = 1,NPD
        if( LLA(NP) > 1 )then
          call DRF_DEPTH(LLA(NP),NP,BELVLA(NP),HPLA(NP))
        endif
      enddo

    endif

    ! *** INITIALIZE DAYNEXT VALUE
    DAYNEXT = DBLE(INT(TIMEDAY)) - 1E-5
    NSTOPPED = 0

    ! *** Initialize to zero the send/receive arrays only once
    drifter_ids_send_west  = 0
    drifter_ids_send_east  = 0
    drifter_ids_send_north = 0
    drifter_ids_send_south = 0
  
    drifter_ids_recv_west  = 0
    drifter_ids_recv_east  = 0
    drifter_ids_recv_north = 0
    drifter_ids_recv_south = 0
  endif ! *** End first call to this routine
  
  ! *** BUILD DAY LIST OF DRIFTERS
  if( TIMEDAY > DAYNEXT )then
      
    ! *** REPORT ANY DRIFTERS THAT MAY HAVE STOPPED FOR WHATEVER REASON
    if( SUM(NSTOPPED) > 0 )then
      if( NSTOPPED(1) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS EXPIRED DUE TO EVAPOR. & BIODEG.:   ',NSTOPPED(1), ' ON PROCESS: ', process_id
      if( NSTOPPED(2) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS STOPPED DUE TO DEPOSITION:          ',NSTOPPED(2), ' ON PROCESS: ', process_id
      if( NSTOPPED(3) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO HYDRAULIC STRUCTURE: ',NSTOPPED(3), ' ON PROCESS: ', process_id
      if( NSTOPPED(4) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO OPEN BOUNDARY:       ',NSTOPPED(4), ' ON PROCESS: ', process_id
      if( NSTOPPED(5) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO WITHDAWAL:           ',NSTOPPED(5), ' ON PROCESS: ', process_id
      if( NSTOPPED(6) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO WITHDRAWAL/RETURN:   ',NSTOPPED(6), ' ON PROCESS: ', process_id
      if( NSTOPPED(7) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS STOPPED DUE TO CELL DRYING:         ',NSTOPPED(7), ' ON PROCESS: ', process_id
      NSTOPPED = 0
    endif
    call MPI_barrier(MPI_Comm_World, ierr)

    TODAY = DBLE(INT(TIMEDAY))
    DAYNEXT = TODAY + 1.
    
    ! *** Update the active drifters for the next day
    NPDAY = 0
    do NP = 1,NPD
      ! *** Exclude drifters not in current domain
      if( LLA_Process(NP) /= process_id )then
        NPP = 0 ! FOR DEBUGGING
        CYCLE
      endif
      
      ! *** Exclude stopped drifters
      if( JSPD(NP) == 0 .and. LLA_Global(NP) < 2 )then
        NPP = 0 ! FOR DEBUGGING
        CYCLE
      endif

      ! *** Exclude drifters whose timing is not in the current time window
      ! *** LA_BEGTI(NP) < LA_BEGTI0:  Released before start of LPT computations
      ! *** LA_BEGTI(NP) > DAYNEXT:    Released after next day
      ! *** LA_ENDTI(NP) < TODAY:      Stopped tracking before current day
      if( LA_BEGTI(NP) < LA_BEGTI0 .or. LA_BEGTI(NP) > DAYNEXT .or. LA_ENDTI(NP) < TODAY )then
        NPP = 0 ! FOR DEBUGGING
        CYCLE
      endif
      
      NPDAY = NPDAY + 1
      NPLIST(NPDAY) = NP
    enddo
    PRINT '(A,I8,A,I8)','# DRIFTERS ACTIVE FOR THE NEXT PERIOD:         ',NPDAY, ' ON PROCESS: ', process_id

    call MPI_Allreduce(NPDAY, NPDAY_MAX, 1, MPI_Integer, MPI_Max, comm_2d, ierr)
    
  endif

  ! *** NEXT CALL --------------------------------------------------------------
  DAHX1 = 0
  DAHY1 = 0
  DAVZ1 = 0
  DAHX2 = 0
  DAHY2 = 0
  DAVZ2 = 0
  MOC   = 0
  IT = 1
  DIFFH = SQRT(2.*LA_HORDIF*DELTD)
  DIFFV = SQRT(2.*LA_VERDIF*DELTD)

  ! *** These keep track of the number of drifters to communicate between processes
  num_drifters_send_west  = 0
  num_drifters_send_east  = 0
  num_drifters_send_north = 0
  num_drifters_send_south = 0
  
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NPP,NP,BEDGEMOVE,XLA1,YLA1,ZLA1)          &
  !$OMP  PRIVATE(U2NP,V2NP,W2NP,IT,NW,KWEIGHT)                                        &
  !$OMP  PRIVATE(KDX1,KDY1,KDZ1,KDX2,KDY2,KDZ2,KDX3,KDY3,KDZ3,KDX4,KDY4,KDZ4)         &
  !$OMP  FIRSTPRIVATE(DAHX1,DAHY1,DAVZ1,DAHX2,DAHY2,DAVZ2)
  do NPP = 1,NPDAY ! *** Loop over the number of active drifters for this time step
    !$ IT = OMP_GET_THREAD_NUM() + 1
  
    ! *** GET CURRENT ACTIVE DRIFTER
    NP = NPLIST(NPP)

    ! *** If we are outside of the time window we should be tracking this drifter then move onto next drifter
    if( TIMEDAY < LA_BEGTI(NP) .or. TIMEDAY > LA_ENDTI(NP) )then
      if( JSPD(NP) == 0 )then
        ! *** Drifter has been initialized.  Update LLA arrays
        if( LLA(NP) > 1 .and. TIMEDAY > LA_ENDTI(NP) )then
          LLA_Global(NP) = 1
        else
          LLA_Global(NP) = 0
        endif
        LLA(NP) = 0
      endif
      CYCLE
    endif
    
    ! *** Drifter is not active yet and the time is above the beginning of the drifter activation time
    if( JSPD(NP) == 1 .and. TIMEDAY >= LA_BEGTI(NP) )then
      ! *** Initialize Drifter
      call CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)
      JSPD(NP) = 0
    endif

    ! *** Call some subroutines specific to the oil spill module
    if( ISOILSPI(LA_GRP(NP)) == 1 .and. LLA(NP) >= 2 )then
      call OIL_PROC(NP)            ! *** CALCULATE DVOL(NP)
      if( DVOL(NP) <= 1D-9 )then
        LLA(NP) = 0
        LLA_Global(NP) = 1
        if( NPD < 10000 .or. DEBUG )then
          PRINT '(A41,I8)','DRIFTER EXPIRED DUE TO EVAPOR. & BIODEG.:',NP
        else
          NSTOPPED(1) = NSTOPPED(1) + 1
        endif
        CYCLE
      endif
    endif

    ! *** If this particle is fixed in a bed and within the active domain
    if( BEDFIX(LA_GRP(NP)) == 1 .and. LLA(NP) >= 2  )then
      ! *** CHECK FOR DEPOSITION
      if( (ZLA(NP) - BELVLA(NP)) <= 1D-6  )then
        LLA(NP) = 0
        LLA_Global(NP) = 1

        if( NPD < 10000 .or. DEBUG )then
          PRINT '(A38,I8)','DRIFTER HAS STOPPED DUE TO DEPOSITION:',NP
        else
          NSTOPPED(2) = NSTOPPED(2) + 1
        endif

        CYCLE
      endif
      
      ! *** CHECK FOR DRY CELL
      if( ISDRY > 0 )then
        if( .not. LMASKDRY(LLA(NP))  )then
          LLA(NP) = 0
          LLA_Global(NP) = 1
          if( NPD < 10000 .or. DEBUG )then
            PRINT '(A38,I8)','DRIFTER HAS STOPPED DUE TO CELL DRYING:',NP
          else
            NSTOPPED(7) = NSTOPPED(7) + 1
          endif

          CYCLE
        endif
      endif
    endif
    
    ! *** If the cell index is less than 2 the drifter has stopped (0), left the domain (0) or exited the sub-domain (1)
    if( LLA(NP) < 2 ) CYCLE

    ! *** If the cell has dried,
    ! *** checks the water surface elevation against the depth at which we have decided a cell goes dry --> somewhat redundant?
    if( ISDRY > 0 .and. HP(LLA(NP)) < HDRY ) CYCLE
    
    ! *** Set the coordinates for the current cell
    XLA1 = XLA(NP)
    YLA1 = YLA(NP)
    ZLA1 = ZLA(NP)

    NW = LA_GRP(NP) ! *** Set the group for this drifter, this was read in by DRIFTER_INP
    BEDGEMOVE = .FALSE.
    
    ! *** APPLY DEPTH CHANGE ON VERTICAL POSITION
    if( LA_ZCAL == 1 .and. ADJVERP == 1 )then
      ZLA1 = ZLA1 + (HP(LLA(NP))-H1P(LLA(NP)))*(ZLA1-BELVLA(NP))/HPLA(NP)
    endif

    ! *** If horitzontal diffusivity is turned off in EE
    if( LA_DIFOP /= 0 .or. ISHDMF == 0 )then
      KDX1 = 0.
      KDY1 = 0.
      KDZ1 = 0.
      KDX2 = 0.
      KDY2 = 0.
      KDZ2 = 0.
      KDX3 = 0.
      KDY3 = 0.
      KDZ3 = 0.
      KWEIGHT = 1.
      GOTO 4
    else
      KWEIGHT = 6.0
    endif
    ! ***************************************************************************************
    ! *** RUNGE-KUTTA
    ! *** PASS 1
    call DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    if( LA_DIFOP == 0 .and. ISHDMF > 0 ) CALL DIFGRAD(LLA(NP),KLA(NP),NP,DAHX1,DAHY1,DAVZ1)
    KDX1 = DELTD*(U2NP + DAHX1)
    KDY1 = DELTD*(V2NP + DAHY1)
    XLA(NP)  = XLA1 + 0.5*KDX1
    YLA(NP)  = YLA1 + 0.5*KDY1
    if( LA_ZCAL == 1 )then
      KDZ1 = DELTD*(W2NP - GRPWS(NW) + DAVZ1)
      ZLA(NP) = ZLA1 + 0.5*KDZ1
    else
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    endif
    call CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)

    if( LLA(NP)<2 .or. BEDGEMOVE) CYCLE

    ! *** PASS 2
    call DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    if( LA_DIFOP == 0 .and. ISHDMF>0 ) CALL DIFGRAD(LLA(NP),KLA(NP),NP,DAHX2,DAHY2,DAVZ2)
    KDX2 = DELTD*(U2NP + DAHX2)
    KDY2 = DELTD*(V2NP + DAHY2)
    XLA(NP) = XLA1 + 0.5*KDX2
    YLA(NP) = YLA1 + 0.5*KDY2
    if( LA_ZCAL == 1 )then
      KDZ2 = DELTD*(W2NP - GRPWS(NW) + DAVZ2)
      ZLA(NP) = ZLA1 + 0.5*KDZ2
    else
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    endif
    call CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)

    if( LLA(NP)<2 .or. BEDGEMOVE) CYCLE

    ! *** PASS 3
    call DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    if( LA_DIFOP == 0 .and. ISHDMF>0 ) CALL DIFGRAD(LLA(NP), KLA(NP), NP, DAHX2, DAHY2, DAVZ2)
    KDX3 = DELTD*(U2NP + DAHX2)
    KDY3 = DELTD*(V2NP + DAHY2)
    XLA(NP) = XLA1 + KDX3
    YLA(NP) = YLA1 + KDY3
    if( LA_ZCAL == 1 )then
      KDZ3 = DELTD*(W2NP - GRPWS(NW) + DAVZ2)
      ZLA(NP) = ZLA1 + KDZ3
    else
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    endif
    call CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)

    if( LLA(NP)<2 .or. BEDGEMOVE ) CYCLE

    ! *** PASS 4
4   continue
    call DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    if( LA_DIFOP == 0 .and. ISHDMF>0 ) CALL DIFGRAD(LLA(NP),KLA(NP),NP,DAHX2,DAHY2,DAVZ2)
    KDX4 = DELTD*(U2NP + DAHX2)
    KDY4 = DELTD*(V2NP + DAHY2)
    XLA(NP) = XLA1 + (KDX1 + 2.0*KDX2 + 2.0*KDX3 + KDX4)/KWEIGHT
    YLA(NP) = YLA1 + (KDY1 + 2.0*KDY2 + 2.0*KDY3 + KDY4)/KWEIGHT
    if( LA_ZCAL == 1 )then
      KDZ4 = DELTD*(W2NP-GRPWS(NW)+DAVZ2)
      ZLA(NP) = ZLA1 + (KDZ1 + 2.0*KDZ2 + 2.0*KDZ3 + KDZ4)/KWEIGHT
    else
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    endif

    ! *** WIND-DRIFT COMPONENT
    if( IOSWD > 0 )then
      if( DLA(NP) < 0.5 ) CALL WIND_DRIFT(LLA(NP), KLA(NP), NP)
    endif

    ! *** RANDOM-WALK COMPONENT
    if( LA_PRAN > 0 ) CALL RANDCAL(LLA(NP), KLA(NP), NP)
    
    call CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)
    
    if( LLA(NP) >= 2 .and. LLA(NP) <= LA .and. ISOILSPI(LA_GRP(NP)) == 1  )then
      MOC(LLA(NP)) = MOC(LLA(NP)) + DVOL(NP)*DDEN(LA_GRP(NP))                ! TOTAL MASS OF OIL (KG) AT A CELL
    endif

    ! *** Check if the drifter has moved into the ghost cell boundary
    call Check_Drifter_In_Ghost(NP)

  enddo
  !$OMP END PARALLEL DO
  
  ! *** Communicate ghost info if a drifter entered the ghost cell region
  ! *** Get the max number of drifters to communicate amongst all directions
  local_tot_drifters_send = num_drifters_send_west + num_drifters_send_east + num_drifters_send_north + num_drifters_send_south
  
  ! *** Get global max number of drifters to communicate amongst all processes
  call MPI_barrier(MPI_Comm_World, ierr)
  
  TTDS = DSTIME(0)
  global_max_drifters_to_comm = 0
  call MPI_Allreduce(local_tot_drifters_send, global_max_drifters_to_comm, 1, MPI_Integer, MPI_Max, comm_2d, ierr)
  DSITIMING(11) = DSITIMING(11) + (DSTIME(0)-TTDS)

  ! *** Only need to communicate if there are drifters in any of the ghost cells of any subdomain
  if( global_max_drifters_to_comm > 0 )then
    ! *** If a drifter is leaving a subdomain let the subdomain know it will be receiving data
    call Notify_Receiving_Domains
    
    ! *** Communicate all drifter info
    call Communicate_Drifters(LLA_Global, 1)

    !Call Communicate_Drifters(JSPD)  ! by definition, if the drifter is communicated JSPD = 0
     
    call Communicate_Drifters(XLA)
     
    call Communicate_Drifters(YLA)

    call Communicate_Drifters(ZLA)

    call Communicate_Drifters(DLA)

    call Communicate_Drifters(HPLA)
     
    call Communicate_Drifters(KLA)
     
    call Communicate_Drifters(BELVLA)

    ! *** Zero send list for next iteration
    drifter_ids_send_west(1:global_max_drifters_to_comm)  = 0
    drifter_ids_send_east(1:global_max_drifters_to_comm)  = 0
    drifter_ids_send_north(1:global_max_drifters_to_comm) = 0
    drifter_ids_send_south(1:global_max_drifters_to_comm) = 0

    drifter_ids_recv_west(1:global_max_drifters_to_comm)  = 0
    drifter_ids_recv_east(1:global_max_drifters_to_comm)  = 0
    drifter_ids_recv_north(1:global_max_drifters_to_comm) = 0
    drifter_ids_recv_south(1:global_max_drifters_to_comm) = 0
  endif ! *** End sequence that communicates drifters 
  
  ! *** WRITE THE CURRENT TRACK POSITION
  if( TIMEDAY >= TIMENEXT )then   ! .or. TIMEDAY+1E-5 >= TIMEEND )then
    call DRIFTER_OUT(.FALSE.)
    TIMENEXT = TIMENEXT + LA_FREQ
  endif
  
END SUBROUTINE DRIFTER_CALC
!---------------------------------------------------------------------------!
!> @details Lets other process know that it will be receiving drifter data
!  and how many drifters 
!
!> @author Zander Mausolff 
!
!---------------------------------------------------------------------------! 
Subroutine Notify_Receiving_Domains

	implicit none

  ! *** Local variables
  integer :: status_msg(MPI_Status_Size)
  integer :: ierr, ii, iii, iFound, NP
  integer :: total_send_drifters
     
  ! *** Determine the number of drifters to send in each direction
  send_east   = num_drifters_send_east
  send_west   = num_drifters_send_west
  send_north  = num_drifters_send_north
  send_south  = num_drifters_send_south

  ! *** These hold the number of drifters that each subdomain will recieve
  recv_east  = 0 
  recv_west  = 0 
  recv_north = 0 
  recv_south = 0 
     
  total_send_drifters = send_east + send_west + send_north + send_south
     
  ! *** Let each process know if it is going to be recieving data

  ! *** West/East send/recv
  call MPI_Sendrecv(drifter_ids_send_west, global_max_drifters_to_comm, MPI_Integer, nbr_west, 0, &
                    drifter_ids_recv_east, global_max_drifters_to_comm, MPI_Integer, nbr_east, 0, &
                    comm_2d, status_msg, ierr)
     
  ! *** East/West send/recv
  call MPI_Sendrecv(drifter_ids_send_east, global_max_drifters_to_comm, MPI_Integer, nbr_east, 0, &
                    drifter_ids_recv_west, global_max_drifters_to_comm, MPI_Integer, nbr_west, 0, &
                    comm_2d, status_msg, ierr)

  ! ***
  call MPI_Sendrecv(drifter_ids_send_north, global_max_drifters_to_comm, MPI_Integer, nbr_north, 0, &
                    drifter_ids_recv_south, global_max_drifters_to_comm, MPI_Integer, nbr_south, 0, &
                    comm_2d, status_msg, ierr)     

  ! ***
  call MPI_Sendrecv(drifter_ids_send_south, global_max_drifters_to_comm, MPI_Integer, nbr_south, 0, &
                    drifter_ids_recv_north, global_max_drifters_to_comm, MPI_Integer, nbr_north, 0, &
                    comm_2d, status_msg, ierr)
     
  recv_east  = minloc(drifter_ids_recv_east,1)  - 1
  recv_west  = minloc(drifter_ids_recv_west,1)  - 1
  recv_north = minloc(drifter_ids_recv_north,1) - 1
  recv_south = minloc(drifter_ids_recv_south,1) - 1
     
  ! *** Update NPDAY for recevied drifters.  First check for drifters already in NPDAY list
  do ii = 1,recv_east
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_east(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_east(ii))
    endif
  enddo
    
  do ii = 1,recv_west
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_west(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_west(ii))
    endif
  enddo
    
  do ii = 1,recv_north
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_north(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_north(ii))
    endif
  enddo

  do ii = 1,recv_south
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_south(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_south(ii))
    endif
  enddo
  
End subroutine Notify_Receiving_Domains


!---------------------------------------------------------------------------!
!> @details Reads the drifter restart file
!> @todo Enable remapping for domain decomposition
!---------------------------------------------------------------------------!  
SUBROUTINE DRIFTER_RST
  ! *** READ LPT RESTART FILE
  integer :: NP
  logical FEXIST

  write(*,'(A)')'READING RESTART FILE: DRIFTER.RST'

  ! CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE THE ASCII FILE INSTEAD.
  RESFILE = 'drifter.rst'
  INQUIRE(FILE = RESFILE, EXIST = FEXIST)

  if( FEXIST )then
    open(99,FILE = RESFILE)

    read(99,'(50I2)') JSPD(1:NPD)
    read(99,'(30I3)') KLA(1:NPD)
    read(99,'(2(I8,2F15.5,F10.5))') (LLA_Global(NP), XLA(NP), YLA(NP), ZLA(NP),NP = 1,NPD)

    close(99)
  endif

END SUBROUTINE DRIFTER_RST
!---------------------------------------------------------------------------!
!> @details Writes to the EE_drifter.out binary file for EE linkage
!
!> @author Zander Mausolff (Updated for MPI)
!> @date 2/17/2020
!---------------------------------------------------------------------------!  
SUBROUTINE DRIFTER_OUT(FORCEEXIT)
  real(RKD)           :: TTDS, TTDS1                     ! MODEL TIMING TEMPORARY VARIABLES

  logical, intent(IN) :: FORCEEXIT

  ! *** Local variables
  logical :: LFORCE
  integer(IK4) :: NP, NACT
  integer(4) :: IX,IY,IZ,IV     ! *** USING INTEGER*4 EXPLICITY FOR EE LINKAGE REQUIREMENTS

  ! *** All processes must send their data
  TTDS = DSTIME(0)
  
  call Gather_Drifter_Arrays
  
  TTDS1 = DSTIME(0) - TTDS
  TMPIEE = TMPIEE + TTDS1
  TLRPD = TLRPD - TTDS1

  ! *** Only proceed with master process
  if( process_id == master_id )then

    if( FORCEEXIT )then
      ! *** USER INITIATED EXIT.  WRITE FINAL POSITION OF EVERY DRIFTER
      PRINT *,'COMPLETE DRIFTER SNAPSHOT (USER): ', TIMEDAY
      LFORCE = .TRUE.
      NACT = NPDAY
    else
      if( TIMEDAY+LA_FREQ >= TIMEEND )then
        ! *** WRITE FINAL POSITION OF EVERY DRIFTER AT END OF RUN
        LFORCE = .TRUE.
        PRINT *,'COMPLETE DRIFTER SNAPSHOT (EoR): ', TIMEDAY
        NACT = NPD
      else
        LFORCE = .FALSE.
        NACT = 0
        do NP = 1,NPD
          if( LLA_Global(NP) >= 1 )then
            NACT = NACT + 1
          endif
        enddo
      endif
    endif

    call write_nc_lpt()
    
    ! *** WRITE TO THE FILE
    open(ULGR,FILE = OUTDIR//'EE_DRIFTER.OUT',STATUS = 'UNKNOWN',FORM = FMT_BINARY,POSITION = 'APPEND')
    EETIME = TIMEDAY

    write(ULGR) EETIME

    NTIMES = NTIMES + 1

    write(ULGR) INT(NACT,4)

    if( ANY(ISOILSPI == 1) )then
      do NP = 1,NPD
        if( LLA_Global(NP) >= 1 .or. LFORCE )then
          IX = NINT(XYZSCL*XLA(NP))                                     ! *** SCALED
          IY = NINT(XYZSCL*YLA(NP))                                     ! *** SCALED
          IZ = NINT(XYZSCL*ZLA(NP))                                     ! *** SCALED
          IV = NINT(1D6*DVOL(NP))                                       ! *** IN CUBIC CM
          write(ULGR) INT(NP,4), INT(LLA_Global(NP),4), IX, IY, IZ, IV
        endif
      enddo
    else
      do NP = 1,NPD
        if( LLA_Global(NP) >= 1 .or. LFORCE )then
          IX = NINT(XYZSCL*XLA(NP))                                     ! *** SCALED
          IY = NINT(XYZSCL*YLA(NP))                                     ! *** SCALED
          IZ = NINT(XYZSCL*ZLA(NP))                                     ! *** SCALED
          write(ULGR) INT(NP,4), INT(LLA_Global(NP),4), IX, IY, IZ      
        endif
      enddo
    endif
    FLUSH(ULGR)

    close(ULGR,STATUS = 'KEEP')
  endif
  
  ! *** Write out the last position then set to inactive for future writes
  do NP = 1,NPD
    if( LLA_Global(NP) == 1 )then
      LLA_Global(NP) = 0                  
    endif
  enddo
      
END SUBROUTINE DRIFTER_OUT
!---------------------------------------------------------------------------!
!
!> @details READING INPUT DATA OF INITIAL LOCATIONS OF DRIFTERS and mapping them out
!! to subdomains if necessary
!!OUTPUT: NPD,XLA,YLA,ZLA,NP = 1:NPD, LA_BEGTI, LA_ENDTI, LA_FREQ,LANDT\
!
! Below are some options that are read in for the drifters
! ISPD:    0 Drifter calculations not performed
!          2 Lagrangian Particle Tracks Computed using Runge-Kutta 4
! ZOPT:    Lagrangian Particle Tracking Option in the Z direction:
!          0 Depths are Fixed at the Initial Seeding Depth
!          1 Fully 3D Lagrangian Neutrally Buoyant Particles
! PRAN:    Random Walk Option to add a random movement
!          0 No random component
!          1 Random Walk, Horizontal ONLY
!          2 Random Walk, Vertical   ONLY
!          3 Random Walk, 3D Random  FULL
! DIFOP    Option For Random Walk Diffusivity
!          0 use AH(L,K) and AV(L,K) from EFDC computations (Horizontal diffusion should be turned on in EFDC)
!          1 use HORDIF  VERDIF from this file
! HORDIF   Horizontal Diffusivity (m^2/s) in case of DIFOP = 1
!          Horizontal Factor for derivative in case of DIFOP = 0
! VERDIF   Vertical Diffusivity   (m^2/s) in case of DIFOP = 1
!          Vertical Factor for derivative in case of DIFOP = 0
! DEPOP:   Option for specifying initial vertical position in input file
!          0 Elevations are specified
!          1 Depths are specified
!---------------------------------------------------------------------------!
SUBROUTINE DRIFTER_INP
  ! ********************************************************************
  implicit none

  ! *** Local variables

  integer :: I, J, NP, N, K, NPN, IN, JN, IS, JS, L, NC
  real(RKD) :: RANVAL
  !REAL(8),external :: DRAND   !IT NEEDS THIS STATEMENT IN CASE OF implicit none @todo Zander 2/14/20 - seems like it is not?
  character(200)  :: STR
  real(RKD)  :: Distance_to_Cells(LA_Global)

  integer :: Min_L_loc_index(1)
  integer :: I_global, J_global, L_local, LLO, NI, i_local, j_local
  integer, Allocatable, Dimension(:) :: new_seed
  integer :: seed_size !< size of the array that represents the random number seed
  integer :: LG

  ! *** Do reading only on the master process
  if(process_id == master_id )then
    write(*,'(A)')'READING DRIFTER.INP'
    open(ULOC,FILE = 'drifter.inp',ACTION = 'READ')
    call SKIPCOM(ULOC,'*')
    read(ULOC,*) LA_ZCAL, LA_PRAN, LA_DIFOP, LA_HORDIF, LA_VERDIF, DEPOP, ADJVERP, SLIPFACTOR, IOSWD, OSWDA, OSWDB
    call SKIPCOM(ULOC,'*')
    read(ULOC,*) LA_BEGTI0, LA_ENDTI0, LA_FREQ
    call SKIPCOM(ULOC,'*')
    read(ULOC,*) NPD,NGRP
    write(*,'(A41,I8)')'DRIFTER: NUMBER OF DRIFTERS INITIALZED: ',NPD
  endif

  ! *** Broadcast read in variables to all processes
  call Broadcast_Scalar(LA_ZCAL,    master_id)
  call Broadcast_Scalar(LA_PRAN,    master_id)
  call Broadcast_Scalar(LA_DIFOP,   master_id)
  call Broadcast_Scalar(LA_HORDIF,  master_id)
  call Broadcast_Scalar(LA_VERDIF,  master_id)
  call Broadcast_Scalar(DEPOP,      master_id)
  call Broadcast_Scalar(ADJVERP,    master_id)
  call Broadcast_Scalar(SLIPFACTOR, master_id)
  call Broadcast_Scalar(IOSWD,      master_id)
  call Broadcast_Scalar(OSWDA,      master_id)
  call Broadcast_Scalar(OSWDB,      master_id)

  call Broadcast_Scalar(LA_BEGTI0,  master_id)
  call Broadcast_Scalar(LA_ENDTI0,  master_id)
  call Broadcast_Scalar(LA_FREQ,    master_id)
  call Broadcast_Scalar(NPD,        master_id)
  call Broadcast_Scalar(NGRP,       master_id)

  ! *** Drifter output frequency
  LA_FREQ = LA_FREQ/1440.
  if( LA_FREQ <= 0. ) LA_FREQ = DBLE(TIDALP)/DBLE(NPPPH)/DBLE(86400.)

  ! NG      : Order number of group
  ! ISOILSPI: = 1 for oil spill; = 0 for normal drifter as usual
  ! GRPWS   : Settling velocity of drifter [m/s]
  ! GVOL    : Total volume of oil spill for a group [m^3]
  ! DDEN    : Density of an oil dritfer    [kg/m^3]
  ! DRAT    : First order biodegradation rate of an oil drifter [1/day,> = 0]
  ! DTEM    : Biodegradation Rate Reference Temperature [deg C]
  ! DVAP    : Vapor pressure of an oil drifter     [Pa]
  ! DVMO    : Molar volume of an oil drifter at STP [m^3/mol]

  ! *** ALLOCATE THE DRIFTER ARRAYS
  allocate(XLA(NPD), YLA(NPD), ZLA(NPD), DLA(NPD))
  allocate(LLA(NPD), KLA(NPD), HPLA(NPD), BELVLA(NPD))
  allocate(LA_BEGTI(NPD),LA_ENDTI(NPD), LA_GRP(NPD), JSPD(NPD))
  allocate(ACTIVE(NPD), NPLIST(NPD))

  allocate(ZCTR(LCM, 0:KC+1))
  allocate(NS_Connectors(LCM), EW_Connectors(LCM))
  allocate(BEDFIX(NGRP))
  allocate(DARE(NGRP), DTEM(NGRP), GRPWS(NGRP), ISOILSPI(NGRP), DVOL0(NGRP))
  allocate(DVAP(NGRP), DVMO(NGRP), DDEN(NGRP), DRAT(NGRP), GVOL(NGRP))

  allocate(MOC(LA))

  ! *** New variables for MPI - global arrays for the initial read and for writing out
  allocate(LLA_Global(NPD))

  !< @todo resize these to reduce memory usage
  allocate(drifter_ids_send_west(NPD))
  allocate(drifter_ids_send_east(NPD))
  allocate(drifter_ids_send_north(NPD))
  allocate(drifter_ids_send_south(NPD))

  allocate(drifter_ids_recv_east (NPD))
  allocate(drifter_ids_recv_west (NPD))
  allocate(drifter_ids_recv_north(NPD))
  allocate(drifter_ids_recv_south(NPD))

  allocate(LLA_Process(NPD))

  ! *** INIITALIZE THE ARRAYS (NPD)
  XLA = 0.0
  YLA = 0.0
  ZLA = 0.0
  DLA = 0.0
  HPLA = 0.0
  BELVLA = 0.0
  LA_BEGTI = 0.0
  LA_ENDTI = 0.0
  LA_GRP = 0
  LLA = 0
  KLA = 0
  JSPD= 0
  MOC = 0
  ACTIVE = .TRUE.
  NPLIST = 0


  ! *** INIITALIZE THE ARRAYS (NGRP)
  GRPWS = 0.0
  DARE = 0.0
  DTEM = 0.0
  DVAP = 0.0
  DVMO = 0.0
  DDEN = 0.0
  DRAT = 0.0
  GVOL = 0.0
  DVOL0 = 0
  ISOILSPI = 0

  ZCTR = 0.0

  ! *** end initialize new MPI arrays

  ! *** Do reading only on master
  if(process_id == master_id )then
    call SKIPCOM(ULOC,'*')
    do N = 1,NGRP
      read(ULOC,'(A)') STR
      read(STR,*) K,ISOILSPI(K),GRPWS(K),BEDFIX(K)
      if( ISOILSPI(K) > 0 )then
        read(STR,*) K,ISOILSPI(K),GRPWS(K),BEDFIX(K),GVOL(K),DDEN(K),DRAT(K),DTEM(K),DVAP(K),DVMO(K)
        DRAT(K) = DRAT(K)/86400.
      endif
    enddo
  endif

  ! *** Broadcast arrays just read in by the master process
  call Broadcast_Array(ISOILSPI, master_id)
  call Broadcast_Array(GRPWS   , master_id)
  call Broadcast_Array(BEDFIX  , master_id)
  call Broadcast_Array(GVOL    , master_id)
  call Broadcast_Array(DDEN    , master_id)
  call Broadcast_Array(DRAT    , master_id)
  call Broadcast_Array(DTEM    , master_id)
  call Broadcast_Array(DVAP    , master_id)
  call Broadcast_Array(DVMO    , master_id)

  ! *** INITIALIZE FOR THE FIRST CALL
  LA_FIRSTCALL = 1
  JSPD = 1           ! *** JSPD = 1 IF NOT INITIALIZED,  JSPD = 0 AFTER CELL IS INTIALIZED
  
  NS_Connectors = 0  ! *** Default - No connection
  do NPN = 1,NPNSBP
    IN = INPNS(NPN)
    JN = JNPNS(NPN)
    IS = ISPNS(NPN)
    JS = JSPNS(NPN)
    NS_Connectors(LIJ(IN,JN)) = 1
    NS_Connectors(LIJ(IS,JS)) = 2
  enddo

  EW_Connectors = 0  ! *** Default - No connection
  do NPN = 1,NPEWBP
    IN = IWPEW(NPN)
    JN = JWPEW(NPN)
    IS = IEPEW(NPN)
    JS = JEPEW(NPN)
    EW_Connectors(LIJ(IN,JN)) = 1
    EW_Connectors(LIJ(IS,JS)) = 2
  enddo

  ! *** Do reading only on master
  if(process_id == master_id )then
    call SKIPCOM(ULOC,'*')
    if( DEPOP == 1 )then
      do NP = 1,NPD
        ! *** Read Depths
        read(ULOC,*,err = 999) XLA(NP), YLA(NP), DLA(NP), LA_BEGTI(NP), LA_ENDTI(NP), LA_GRP(NP)
        DLA(NP) = MAX(DLA(NP),0.0_RKD)
        if( LA_GRP(NP) < 1 ) LA_GRP(NP) = 1
      enddo
    else
      do NP = 1,NPD
        ! *** Read Elevations
        read(ULOC,*,err = 999) XLA(NP), YLA(NP), ZLA(NP), LA_BEGTI(NP), LA_ENDTI(NP), LA_GRP(NP)
        if( LA_GRP(NP) < 1 ) LA_GRP(NP) = 1
      enddo
    endif
    close(ULOC)
  endif

  ! *** Broadcast the arrays to all processes
  call Broadcast_Array(XLA, master_id)
  call Broadcast_Array(YLA, master_id)
  if( DEPOP == 1  )then
    call Broadcast_Array(DLA, master_id)
  else
    call Broadcast_Array(ZLA, master_id)
  endif
  call Broadcast_Array(LA_BEGTI, master_id)
  call Broadcast_Array(LA_ENDTI, master_id)
  call Broadcast_Array(LA_GRP,   master_id)

  ! ***************************************************************************************************************
  
  ! *** Initialize new MPI global arrays
  drifter_ids_send_west  = 0
  drifter_ids_send_east  = 0
  drifter_ids_send_north = 0
  drifter_ids_send_south = 0

  drifter_ids_recv_east  = 0
  drifter_ids_recv_west  = 0
  drifter_ids_recv_north = 0
  drifter_ids_recv_south = 0

  LLA_Process = -1               ! *** Not Current Processor

  ! *** Now we need to map global values read in to local values
  do NP = 1, NPD
    ! *** First need to determine which cell we are in, drifter values are read in based on x,y,z
    LLA(NP) = -1
    LLA_Global(NP) = -1
    
    if( LA_BEGTI(NP) >= LA_BEGTI0 .and. LA_BEGTI(NP) <= LA_ENDTI0 )then
      ! *** Get Local L
      do J = 3,JC-2
        do I = 3,IC-2
          L = LIJ(I,J)
          if( L > 0 )then
            ! *** Test if Drifter is inside the cell
            if( INSIDECELL(L, XLA(NP), YLA(NP)) )then
              LLA(NP) = L
              ! *** LLA_Global(NP) is initialized only when activated
              LLA_Process(NP) = Map2Global(L).PR
              Exit
            endif
          endif
        enddo
        if( LLA(NP) > 1 ) EXIT
      enddo
    endif
  enddo

  ! *** Write out information regaring the local mapping of drifters
  if(MPI_DEBUG_FLAG )then
    write(mpi_efdc_out_unit,'(a,i8)') 'Total number of local drifters: ', NPD

    call writebreak(mpi_efdc_out_unit)
    write(mpi_efdc_out_unit,'(a)') 'Writing out local drifter values'
    write(mpi_efdc_out_unit,'(a)') '               NP,      LLA, LLA_GL,  Proc ID, XLA(NP),      YLA(NP),          DLA(NP),    LA_BEGTI(NP), LA_ENDTI(NP),  LA_GRP(NP)'
    ! *** Write out - TESTING/DEBUGGING only
    do np = 1, NPD

      ! *** If the cell is inside the subdomain
      if( DEPOP == 1  )then
        write(mpi_efdc_out_unit,'(a,4i8,5f15.4,i8)') 'Drifter: ', NP, LLA(NP), LLA_Global(NP), LLA_Process(NP), XLA(NP), YLA(NP), ZLA(NP), LA_BEGTI(NP),LA_ENDTI(NP),LA_GRP(NP)
      else
        write(mpi_efdc_out_unit,'(a,4i8,5f15.4,i8)') 'Drifter: ', NP, LLA(NP), LLA_Global(NP), LLA_Process(NP), XLA(NP), YLA(NP), DLA(NP), LA_BEGTI(NP),LA_ENDTI(NP),LA_GRP(NP)
      endif

    enddo
  endif

  !< @todo Improve random number generation
  ! *** This is where the random number is first set.
  ! *** This should be adjusted to instead call the random number seed on all processes
  ! *** and then call the random number.  That way we can make repeatable events with pseudo random numbers

  ! *** Setup the random seed
  seed_size = 3
  allocate(new_seed(seed_size))
  new_seed(:) = process_id + 1
  call RANDOM_SEED(put = new_seed(1:seed_size))

  if( LA_PRAN > 0  )then
    call RANDOM_NUMBER(RANVAL)  ! *** Replaced 2/14/2020 --> RANVAL = DRAND(1)
  endif

  ! *** SETUP L LOOKUPS FOR HYDRAULIC STRUCTURES
  allocate(LCTLU(MAX(NQCTL,1)),LCTLD(MAX(NQCTL,1)))
  LCTLU = 0
  LCTLD = 0
  do NC = 1,NQCTL
    LCTLU(NC) = LIJ(HYD_STR(NC).IQCTLU, HYD_STR(NC).JQCTLU)
    if( HYD_STR(NC).IQCTLD > 0 .and. HYD_STR(NC).JQCTLD > 0 )then
      LCTLD(NC) = LIJ(HYD_STR(NC).IQCTLD, HYD_STR(NC).JQCTLD)
    endif
  enddo

  ! *** SETUP L LOOKUPS FOR HYDRAULIC STRUCTURES
  allocate(LWRU(MAX(NQWR,1)), LWRD(MAX(NQWR,1)))
  LWRU = 0
  LWRD = 0
  do NC = 1,NQWR
    LWRU(NC) = LIJ(WITH_RET(NC).IQWRU, WITH_RET(NC).JQWRU)
    if( WITH_RET(NC).IQWRD > 0 .and. WITH_RET(NC).JQWRD > 0 )then
      LWRD(NC) = LIJ(WITH_RET(NC).IQWRD, WITH_RET(NC).JQWRD)
    endif
  enddo

  return
  999 call STOPP('DRIFTER.INP READING ERROR!')

END SUBROUTINE DRIFTER_INP

!---------------------------------------------------------------------------!
!> @details DETERMINING LLA,KLA,BELVLA,HPLA FOR THE FIRST CALL
!! UPDATING XLA,YLA,LLA,KLA,BELVLA,HPLA FOR THE NEXT CALL
!! FOR EACH DRIFTER (XLA,YLA,ZLA)
!! BY FINDING THE NEAREST CELL CENTROID
!! THEN EXPANDING TO THE NEIGHBOUR CELLS
!! HP(LIJ(I,J))     : WATER DEPTH = WATER SUR. - BELV
!! BELV(LIJ(I,J))   : BOTTOM ELEVATION OF A CELL
!! BELVLA           : BED ELEVATION AT DRIFTER NI POSITION
!! HPLA             : WATER DEPTH AT DRIFTER NI POSITION
!! DLON(L),L = 2:LA ? : CELL CENTROID XCEL = XCOR(L,5)
!! DLAT(L),L = 2:LA ? : CELL CENTROID YCEL = YCOR(L,5)
!
! @param[inout] BEDGEMOVE
! @param[in] XLA2 - x coordinates
! @param[in] YLA2 - y coordinates
! @param[in] ZLA2 - z coordinates
! @param[in] NI   - Current drifter particle
!
! Some other notes:
! *** INPUT:
! *** IF DEPOP = 0: XLA,YLA,ZLA,XCOR(L,5),YCOR(L,5),BELV,HP
! *** IF DEPOP = 1: XLA,YLA,    XCOR(L,5),YCOR(L,5),BELV,HP,DLA
! *** OUTPUT:
! ***             XLA,YLA,LLA(NP),KLA(NP),BELVLA(NP),HPLA(NP)
!---------------------------------------------------------------------------!

SUBROUTINE CONTAINER(BEDGEMOVE, XLA2, YLA2, ZLA2, NI)

  implicit none

  logical(4),intent(INOUT) :: BEDGEMOVE 
  real(RKD) ,intent(IN)    :: XLA2 
  real(RKD) ,intent(IN)    :: YLA2 
  real(RKD) ,intent(IN)    :: ZLA2 
  integer,intent(IN)       :: NI   

  ! *** Local variables
  integer :: NPSTAT, NWR, NS, NPN, LLO, LI, NCTL, ipmc
  integer :: LMILOC(1), K, L, ILN, JLN, LE, LN, LM, LLL, LLA2
  integer :: LL, LU, LD
  real(RKD)  :: RADLA(LA), VELEK, VELNK, VPRO, XXO, YYO, DMIN, D, X0, Y0, X3, Y3, XP, YP, UTMPB, VTMPB, OFF, PMCTIME
  logical(4) :: LINSIDE

  ! *** DETERMINE THE NEAREST CELL CENTROID
  if( JSPD(NI) == 1 )then
    ! *** GET CELL FOR THE FIRST CALL
    ! *** Determines distance from the current drifter position to the centroid of all cells
    RADLA(2:LA) = SQRT((XLA(NI)-XCOR(2:LA,5))**2+(YLA(NI)-YCOR(2:LA,5))**2)
    
    ! *** Finds the closest cell by getting the minimum value
    LMILOC = MINLOC(RADLA(2:LA))
    ILN = IL(LMILOC(1)+1)    !I OF THE NEAREST CELL FOR DRIFTER
    JLN = JL(LMILOC(1)+1)    !J OF THE NEAREST CELL FOR DRIFTER
    LLO = LIJ(ILN,JLN)
        
    ! *** Initialization of ZLA/DLA, BELVLA, HPLA, KLA parameters on the first call
    
    ! *** CHECK IF DRIFTER IS INSIDE THE CELL
    LINSIDE = .FALSE.
    if( INSIDECELL(LLO,XLA(NI),YLA(NI)) )then
      LINSIDE = .TRUE.
    endif
    
    if( .not. LINSIDE .and. LWC(LLO) > 0 .and. LWC(LLO) <= LA )then
      if( INSIDECELL(LWC(LLO),XLA(NI),YLA(NI)) )then
        LINSIDE = .TRUE.
        LLO = LWC(LLO)
      endif
    endif
    
    if( .not. LINSIDE .and. LEC(LLO) > 0 .and. LEC(LLO) <= LA )then
      if( INSIDECELL(LEC(LLO),XLA(NI),YLA(NI)) )then
        LINSIDE = .TRUE.
        LLO = LEC(LLO)
      endif
    endif
    
    if( .not. LINSIDE .and. LSC(LLO) > 0 .and. LSC(LLO) <= LA )then
      if( INSIDECELL(LSC(LLO),XLA(NI),YLA(NI)) )then
        LINSIDE = .TRUE.
        LLO = LSC(LLO)
      endif
    endif
    
    if( .not. LINSIDE .and. LNC(LLO) > 0 .and. LNC(LLO) <= LA )then
      if( INSIDECELL(LNC(LLO),XLA(NI),YLA(NI)) )then
        LINSIDE = .TRUE.
        LLO = LNC(LLO)
      endif
    endif
    
    ! *** DRIFTER IS IN THE CELL.  DETERMINE ALL OF THE SETTINGS
    LLA(NI) = LLO

    call DRF_DEPTH(LLA(NI),NI,BELVLA(NI),HPLA(NI))

    call DRF_LAYER(LLA(NI),BELVLA(NI),HPLA(NI),KLA(NI),ZLA(NI))

    !THE FIRST CALL: CONVERT DLA TO DLA/ZLA
    if( DEPOP == 1 )then
      ZLA(NI) = HPLA(NI) + BELVLA(NI)-DLA(NI)
    else
      DLA(NI) = MAX(HPLA(NI)+BELVLA(NI)-ZLA(NI),0._8)
    endif
      
    ! *** Settings based on if oil spill calculations are being considers
    if( ISOILSPI(LA_GRP(NI)) == 1 .and. DDEN(LA_GRP(NI)) < 1000. .and. GRPWS(LA_GRP(NI)) == 0. )then
      ! *** FORCE OIL TO SURFACE IF NOT USING SETTLING/RISING VELOCITIES
      DLA(NI) = 0.005
      ZLA(NI) = HPLA(NI) + BELVLA(NI) - DLA(NI)
    endif

    LLA_Global(NI) = Map2Global(LLO).LG

    return
  endif

  ! *************************************  NORMAL PROCESSING  *********************************************
  ! *** CHECK IF OLD AND NEW POSITIONS ARE IN THE SAME CELL
  if( INSIDECELL(LLA(NI), XLA(NI), YLA(NI)) )then
    call DRF_DEPTH(LLA(NI), NI, BELVLA(NI), HPLA(NI))
    call DRF_LAYER(LLA(NI), BELVLA(NI), HPLA(NI), KLA(NI), ZLA(NI))
    
    ! *** CHECK FOR BOUNDARY CONDITIONS
    if( MOD(NITER,3) == 0 ) CALL CHECK_BCS
    
    return
  endif

  ! *************************************  OUTSIDE CELL  *********************************************
  ! *** Either the drifter is outside the domain or has changed cells.  Process the change.
  ! *** Save the drifter information before any changes by EDGEMOVE
  ILN = IL(LLA(NI))
  JLN = JL(LLA(NI))
  LLO = LLA(NI)
  XXO = XLA(NI)
  YYO = YLA(NI)

  ! *** DETERMINE THE CELL CONTAINING THE DRIFTER WITHIN 9 CELLS: LLA(NI)
  NPSTAT  = 0
  LINSIDE = .FALSE.
  LI = 0
  K  = KLA(NI)   ! *** NEEDED FOR PARTIAL LAYER MASKS

  ! *** Loop over the adjacent cells
  do LL = 1,9
    L = LADJ(LL,LLO)
    if( L < 2 .or. L > LA .or. LL == 5 ) CYCLE

    if( INSIDECELL(L,XLA(NI),YLA(NI)) )then
      ! *** PARTICLE IS INSIDE CURRENT CELL FOR NEXT CALL
      ! *** DEALING WITH THE WALLS
      LINSIDE = .TRUE.
      LI = LL
      if( LL == 1 )then
        ! *** NORTHWEST
        if( (SUB(LLO) > 0.5 .and. SVB(L) > 0.5) .or. (SVB(LADJ(2,LLO)) > 0.5 .and. SUB(LADJ(2,LLO)) > 0.5) )then
          if( ISOK(LL, LLO, K, L) )then
            LLA(NI) = L
            NPSTAT  = 1
          endif
        endif
      elseif( LL == 2 .and. SVB(LADJ(2,LLO)) > 0.5 )then
        ! *** NORTH
        if( ISOK(LL, LLO, K, L) )then
          LLA(NI) = L
          NPSTAT  = 1
        endif
      elseif( LL == 3 )then  ! .and. (SVB(LADJ(2,LLO)) > 0.5 .or. SUB(LADJ(6,LLO)) > 0.5) )then
        ! *** NORTHEAST
        if( (SUB(LADJ(6,LLO)) > 0.5 .and. SVB(LADJ(3,LLO)) > 0.5) .or. (SVB(LADJ(2,LLO)) > 0.5 .and. SUB(LADJ(3,LLO)) > 0.5) )then
          if( ISOK(LL, LLO, K, L) )then
            LLA(NI) = L
            NPSTAT  = 1
          endif
        endif
      elseif( LL == 4 .and. SUB(LLO) > 0.5 )then
        ! *** WEST
        if( ISOK(LL, LLO, K, L) )then
          LLA(NI) = L
          NPSTAT  = 1
        endif
      elseif( LL == 6 .and. SUB(LADJ(6,LLO)) > 0.5 )then
        ! *** EAST
        if( ISOK(LL, LLO, K, L) )then
          LLA(NI) = L
          NPSTAT  = 1
        endif
      elseif( LL == 7 )then  ! .and. (SVB(L) > 0.5 .or. SUB(L) > 0.5) )then
        ! *** SOUTHWEST
        if( (SUB(LLO) > 0.5 .and. SVB(LADJ(4,LLO)) > 0.5) .or. (SVB(LLO) > 0.5 .and. SUB(LADJ(8,LLO)) > 0.5) )then
          if( ISOK(LL, LLO, K, L) )then
            LLA(NI) = L
            NPSTAT  = 1
          endif
        endif
      elseif( LL == 8 .and. SVB(LLO) > 0.5 )then
        ! *** SOUTH
        if( ISOK(LL, LLO, K, L) )then
          LLA(NI) = L
          NPSTAT  = 1
        endif
      elseif( LL == 9 )then  ! .and. (SVB(L) > 0.5 .or. SUB(LE) > 0.5) )then
        ! *** SOUTHEAST
        if( (SUB(LADJ(6,LLO)) > 0.5 .and. SVB(LADJ(6,LLO)) > 0.5) .or. (SVB(LLO) > 0.5 .and. SUB(LADJ(9,LLO)) > 0.5) )then
          if( ISOK(LL, LLO, K, L) )then
            LLA(NI) = L
            NPSTAT  = 1
          endif
        endif
      endif

      EXIT

    else
      ! *** NOT IN CELL L, CHECK SPECIAL CONDITIONS
      
      ! Possible drifter got moved into a cell with a connector?
      if( NS_Connectors(LLO) >= 1 .and. NS_Connectors(L) >= 1 )then        ! *** N-S connection
        if( (L == LNC(LLO) .and. SVB3D(L,K) > 0.5) .or. (L == LSC(LLO) .and. SVB3D(LLO,K) > 0.5) )then
          VELEK = CUE(LLO)*U(LLO,KC) + CVE(LLO)*V(LLO,KC)
          VELNK = CUN(LLO)*U(LLO,KC) + CVN(LLO)*V(LLO,KC)
          VPRO  = VELEK*(XCOR(L,5)-XCOR(LLO,5)) + VELNK*(YCOR(L,5)-YCOR(LLO,5))
          if( VPRO > 0 )then
            LLA(NI) = L
            XLA(NI) = XCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            YLA(NI) = YCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            NPSTAT  = 1
            LLA_Global(NI) = Map2Global(LLA(NI)).LG
            EXIT
          endif
        endif
      elseif( EW_Connectors(LLO) >= 1 .and. EW_Connectors(L) >= 1 )then    ! *** E-W connection
        if( (L == LEC(LLO) .and. SUB3D(L,K) > 0.5) .or. (L == LWC(LLO) .and. SUB3D(LLO,K) > 0.5) )then
          VELEK = CUE(LLO)*U(LLO,KC) + CVE(LLO)*V(LLO,KC)
          VELNK = CUN(LLO)*U(LLO,KC) + CVN(LLO)*V(LLO,KC)
          VPRO  = VELEK*(XCOR(L,5)-XCOR(LLO,5)) + VELNK*(YCOR(L,5)-YCOR(LLO,5))
          if( VPRO > 0 )then
            LLA(NI) = L
            XLA(NI) = XCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            YLA(NI) = YCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            NPSTAT  = 1
            LLA_Global(NI) = Map2Global(LLA(NI)).LG
            EXIT
          endif
        endif
      endif
    endif

  enddo
  
  ! *** CHECK IF THE PARTICLE IS INSIDE THE MODEL DOMAIN
  if( NPSTAT == 0 )then
    ! *** CHECK IF THE PARTICLES SHOULD EXIT DOMAIN
    call CHECK_BCS

    ! *** IF LLO IS NOT A BC CELL THEN CONTINUE PROCESSING
    if( LLA(NI) > 1 )then
      ! *** UPDATE POSITION BASED ON EDGE VELOCITIES (HALF OF CELL CENTERED)
      X0 = XLA(NI)
      Y0 = YLA(NI)

      ! **********************************************************************************
      D = 1E32
      if( LINSIDE .and. LI > 0 )then
        ! *** PARTICLE IS INSIDE AN INVALID CELL.  MOVE TO VALID
        LLL = LLO
        if( LI == 2 )then
          OFF = 0.01*DYP(LLO)
          call DIST2LINE(LLO,2,X0,Y0,1,OFF,D,X3,Y3)     ! *** NORTH EDGE

        elseif( LI == 4 )then
          OFF = 0.01*DXP(LLO)
          call DIST2LINE(LLO,1,X0,Y0,1,OFF,D,X3,Y3)     ! *** WEST EDGE

        elseif( LI == 6 )then
          OFF = 0.01*DXP(LLO)
          call DIST2LINE(LLO,3,X0,Y0,1,OFF,D,X3,Y3)     ! *** EAST EDGE

        elseif( LI == 8 )then
          OFF = 0.01*DYP(LLO)
          call DIST2LINE(LLO,4,X0,Y0,1,OFF,D,X3,Y3)     ! *** SOUTH EDGE

        elseif( LL == 1 )then
          ! *** NORTHWEST
          if( SUB3D(LLO,K) > 0.5 )then
            L = LADJ(4,LLO)
            OFF = 0.01*DYP(L)
            call DIST2LINE(L,2,X0,Y0,1,OFF,D,X3,Y3)     ! *** NORTH EDGE
            LLL = L
          elseif( SVB3D(LADJ(2,LLO),K) > 0.5 )then
            L = LADJ(2,LLO)
            OFF = 0.01*DXP(L)
            call DIST2LINE(L,1,X0,Y0,1,OFF,D,X3,Y3)     ! *** WEST EDGE
            LLL = L
          endif
        elseif( LL == 3 )then
          ! *** NORTHEAST
          if( SUB3D(LADJ(6,LLO),K) > 0.5 )then
            L = LADJ(6,LLO)
            OFF = 0.01*DYP(L)
            call DIST2LINE(L,2,X0,Y0,1,OFF,D,X3,Y3)     ! *** NORTH EDGE
            LLL = L
          elseif( SVB3D(LADJ(2,LLO),K) > 0.5 )then
            L = LADJ(2,LLO)
            OFF = 0.01*DXP(L)
            call DIST2LINE(L,3,X0,Y0,1,OFF,D,X3,Y3)     ! *** EAST EDGE
            LLL = L
          endif
        elseif( LL == 7 )then
          ! *** SOUTHWEST
          if( SUB3D(LLO,K) > 0.5 )then
            L = LADJ(4,LLO)
            OFF = 0.01*DYP(L)
            call DIST2LINE(L,4,X0,Y0,1,OFF,D,X3,Y3)     ! *** SOUTH EDGE
            LLL = L
          elseif( SVB3D(LLO,K) > 0.5 )then
            L = LADJ(8,LLO)
            OFF = 0.01*DXP(L)
            call DIST2LINE(L,1,X0,Y0,1,OFF,D,X3,Y3)     ! *** WEST EDGE
            LLL = L
          endif
        elseif( LL == 9 )then
          ! *** SOUTHEAST
          if( SUB3D(LADJ(6,LLO),K) > 0.5 )then
            L = LADJ(6,LLO)
            OFF = 0.01*DYP(L)
            call DIST2LINE(L,4,X0,Y0,1,OFF,D,X3,Y3)     ! *** SOUTH EDGE
            LLL = L
          elseif( SVB3D(LLO,K) > 0.5 )then
            L = LADJ(8,LLO)
            OFF = 0.01*DXP(L)
            call DIST2LINE(L,3,X0,Y0,1,OFF,D,X3,Y3)     ! *** EAST EDGE
            LLL = L
          endif
        endif

        if( ABS(D) < 1.E32 )then
          if( INSIDECELL(LLL,X3,Y3) )then
            BEDGEMOVE = .TRUE.
            XLA(NI) = X3
            YLA(NI) = Y3
            LLA(NI) = LLL
          else
            XLA(NI) = XXO
            YLA(NI) = YYO
            LLA(NI) = LLO
          endif
        else
          XLA(NI) = XXO
          YLA(NI) = YYO
          LLA(NI) = LLO
        endif
        LLA_Global(NI) = Map2Global(LLA(NI)).LG

        call DRF_DEPTH(LLA(NI),NI,BELVLA(NI),HPLA(NI))
        call DRF_LAYER(LLA(NI),BELVLA(NI),HPLA(NI),KLA(NI),ZLA(NI))

        return
      endif

      ! **********************************************************************************
      ! *** PARTICLE IS OUTSIDE DOMAIN

      ! *** FIND CLOSEST FACE CONNECTED CELL (NEEDED TO ADDRESS SMALL GAPS BETWEEN CELLS)
      DMIN = 1E32
      LLL = 0
      do LL = 1,4
        if( LL == 1 .or. LL == 3 )then
          OFF = 0.01*DXP(LLO)
          call DIST2LINE(LLO,LL,X0,Y0,1,OFF,D,X3,Y3)
        else
          OFF =  0.01*DYP(LLO)
          call DIST2LINE(LLO,LL,X0,Y0,1,OFF,D,X3,Y3)
        endif
        if( ABS(D) < ABS(DMIN) .and. ABS(D) < (25.*OFF) )then
          DMIN = D
          LLL = LL
          XP = X3
          YP = Y3
        endif
      enddo

      if( LLL /= 0 )then
        ! *** FOUND AN EDGE, SEE IF IT IS ADJACENT TO AN ACTIVE CELL
        L = 0
        if( LLL == 2 )then
          L = LNC(LLO)             ! *** NORTH
          if( SVB3D(L,K) > 0.5 )then
            LLA(NI) = L
            NPSTAT = 1
          endif
        elseif( LLL == 3 )then
          L = LEC(LLO)             ! *** EAST
          if( SUB3D(L,K) > 0.5 )then
            LLA(NI) = L
            NPSTAT = 1
          endif
        elseif( LLL == 4 )then
          L = LLO                  ! *** SOUTH
          if( SVB3D(L,K) > 0.5 )then
            LLA(NI) = LSC(L)
            NPSTAT = 1
          endif
        elseif( LLL == 1 )then
          L = LLO                  ! *** WEST
          if( SUB3D(L,K) > 0.5 )then
            LLA(NI) = LWC(L)
            NPSTAT = 1
          endif
        endif

        if( NPSTAT == 0 )then
          ! *** MOVE THE PARTICLE BACK INTO THE ORIGINAL CELL
          LLA(NI) = LLO
          NPSTAT = 1
          XLA(NI) = XP
          YLA(NI) = YP
          BEDGEMOVE = .TRUE.

        endif
        LLA_Global(NI) = Map2Global(LLA(NI)).LG
      endif

      if( NPSTAT == 0 )then
        ! *** PARTICLE MUST BE IN ONE OF THE CORNER CELLS
        DMIN = 1E32
        LLL = 1
        do LL = 1,4
          D = SQRT((XLA(NI)-XCOR(LLO,LL))**2 + (YLA(NI)-YCOR(LLO,LL))**2)
          if( D < DMIN )then
            DMIN = D
            LLL = LL
          endif
        enddo

        ! *** PLACE THE PARTICLE IN THE ASSOCIATED CORNER
        BEDGEMOVE = .TRUE.
        NPSTAT = 1
        OFF = 0.01
        X0 = OFF*DXP(LLO)
        Y0 = OFF*DYP(LLO)
        if( LLL == 1 )then
          XLA(NI) = XCOR(LLO,LLL) - (-X0 * CUE(LLO) - Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (-X0 * CUN(LLO) - Y0 * CVN(LLO))
        elseif( LLL == 2 )then
          XLA(NI) = XCOR(LLO,LLL) - (-X0 * CUE(LLO) + Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (-X0 * CUN(LLO) + Y0 * CVN(LLO))
        elseif( LLL == 3 )then
          XLA(NI) = XCOR(LLO,LLL) - (X0 * CUE(LLO) + Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (X0 * CUN(LLO) + Y0 * CVN(LLO))
        elseif( LLL == 4 )then
          XLA(NI) = XCOR(LLO,LLL) - (X0 * CUE(LLO) - Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (X0 * CUN(LLO) - Y0 * CVN(LLO))
        else
          LLA(NI) = 1
          LLA_Global(NI) = 1
          NPSTAT = 0
        endif
      endif
    endif

  else    ! if( NPSTAT == 1 )then
    ! *** DRIFTER IS NOW IN A NEW CELL.  CHECK BOUNDARY CONDITIONS
    call CHECK_BCS
  endif

  ! *** DETERMINE BOTTOM ELEVATION AND TOTAL WATER DEPTH OF DRIFTERS FOR EVERY TIMESTEP
  if( LLA(NI) >= 2 )then
    LLA_Global(NI) = Map2Global(LLA(NI)).LG
    call DRF_DEPTH(LLA(NI),NI,BELVLA(NI),HPLA(NI))
    call DRF_LAYER(LLA(NI),BELVLA(NI),HPLA(NI),KLA(NI),ZLA(NI))
  endif

! ***************************************************************************
contains

!---------------------------------------------------------------------------!
!> @details Checks if drifter cell location is within a boundary condition
!! cell.  (For MPI) checks if the cell has entered the ghost region
!---------------------------------------------------------------------------!
SUBROUTINE CHECK_BCS

  implicit none
  real :: RAD
  
  ! *** LLA(NI) is the L index of the last valid cell the drifter was in.
  
  ! *** Checks if cell id is valid
  if( LLA(NI) < 2 .or. LLA(NI) > LA )then
    return
  endif

  ! *** Check if drifter is within the North, South, East, West open boundary (Only for drifters outside the domain)
  if( NBCSOP > 0 .and. NPSTAT == 0 )then
    if( ANY(LPBN == LLA(NI)) .or. ANY(LPBS == LLA(NI)) .or.  &
        ANY(LPBE == LLA(NI)) .or. ANY(LPBW == LLA(NI)) )   THEN

      call SET_DRIFTER_OUT(XLA2,YLA2,ZLA2)

      if( NPD < 10000 .or. DEBUG )then
        PRINT '(A36,I8)','OPEN BOUNDARY, DRIFTER IS OUTSIDE:',NI
      else
        NSTOPPED(4) = NSTOPPED(4) + 1
      endif
      return
      
    endif
  endif
  
  ! *** Check if drifter in outflow boundary condition
  if( NQSIJ > 0 )then
    if( ANY(BCFL(:).L == LLA(NI)) .and. QSUM(LLA(NI),KLA(NI)) < 0. )then
      ! *** If inside a cell then use cell dimensions to determine if should remove
      if( NPSTAT == 1 )then
        ! *** INSIDE CELL
        RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
      else
        RAD = 0.
      endif
            
      if( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )then
        call SET_DRIFTER_OUT(XLA2,YLA2,ZLA2)

        if( NPD < 10000 .or. DEBUG )then
          PRINT '(A36,I8)','WITHDAWAL CELL, DRIFTER IS OUTSIDE:',NI
        else
          NSTOPPED(5) = NSTOPPED(5) + 1
        endif
      endif    
      return
      
    endif
  endif
  
  ! *** Check if drifter is in a withdrawal cell
  if( NQWR > 0 )then
    if( ANY(LWRU == LLA(NI)) .or. ANY(LWRD == LLA(NI)) )then
      ! *** PMC - MODIFIED TO HANDLE +/- FLOWS FOR THE W/R BOUNDARY CONDITION
      LLA2 = LLA(NI)
      do NWR = 1,NQWR
        ! *** Handle +/- Flows for Withdrawal/Return Structures
        NS = WITH_RET(NWR).NQWRSERQ
        if( QWRSERT(NS) >= 0. )then
          ! *** Original Withdrawal/Return
          if( WITH_RET(NWR).IQWRU == IL(LLA2) .and. WITH_RET(NWR).JQWRU == JL(LLA2) .and. WITH_RET(NWR).KQWRU == KLA(NI) )then
            ! *** If inside a cell then use cell dimensions to determine if should remove
            if( NPSTAT == 1 )then
              ! *** INSIDE CELL
              RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
            else
              RAD = 0.
            endif
            
            if( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )then
              ! ***  RETURN DRIFTER
              LLA(NI)   = LIJ(WITH_RET(NWR).IQWRD,WITH_RET(NWR).JQWRD)
              LLA2      = LLA(NI)
              XLA(NI)   = XCOR(LLA2,5)
              YLA(NI)   = YCOR(LLA2,5)
              ZLA(NI)   = BELV(LLA2) + HP(LLA2)*ZZ(LLA2,WITH_RET(NWR).KQWRD)
              KLA(NI)   = WITH_RET(NWR).KQWRD
              BEDGEMOVE = .TRUE.
            endif
            EXIT
          endif
        else
          ! *** Reverse Flow Withdrawal/Return
          if( WITH_RET(NWR).IQWRD == IL(LLA2) .and. WITH_RET(NWR).JQWRD == JL(LLA2) .and. WITH_RET(NWR).KQWRD == KLA(NI) )then
            ! *** If inside a cell then use cell dimensions to determine if should remove
            if( NPSTAT == 1 )then
              ! *** INSIDE CELL
              RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
            else
              RAD = 0.
            endif
            
            if( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )then
              ! ***  RETURN DRIFTER
              LLA(NI)   = LIJ(WITH_RET(NWR).IQWRU,WITH_RET(NWR).JQWRU)
              LLA2      = LLA(NI)
              XLA(NI)   = XCOR(LLA2,5)
              YLA(NI)   = YCOR(LLA2,5)
              ZLA(NI)   = BELV(LLA2) + HP(LLA2)*ZZ(LLA2,WITH_RET(NWR).KQWRU)
              KLA(NI)   = WITH_RET(NWR).KQWRU
              BEDGEMOVE = .TRUE.
            endif
            EXIT
          endif
        endif
      enddo

      ! ***
      if( LLA(NI) == 1 )then
        ! *** Should never get here
        if( NPD < 10000 .or. DEBUG )then
          PRINT '(A36,I8)','WITHDRAWAL/RETURN, DRIFTER IS OUTSIDE:',NI
        else
          NSTOPPED(6) = NSTOPPED(6) + 1
        endif
      endif
      return
      
    endif
  endif     ! *** END OF NQWR

  if( NQCTL > 0 )then
    if( ANY(LCTLU == LLA(NI)) )then
      ! *** HYDRAULIC STRUCTURE.  RETURN DRIFTER TO DOWNSTREAM CELL, IF ANY
      LLA2 = LLA(NI)
      do NCTL = 1,NQCTL
        LU = LIJ(HYD_STR(NCTL).IQCTLU, HYD_STR(NCTL).JQCTLU)
        if( LU == LLA2 )then
          ! *** If inside a cell then use cell dimensions to determine if should remove
          if( NPSTAT == 1 )then
            ! *** INSIDE CELL
            RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
          else
            RAD = 0.
          endif
            
          if( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )then
            LD = LIJ(MAX(HYD_STR(NCTL).IQCTLD,1), MAX(HYD_STR(NCTL).JQCTLD,1))
            if( LD >= 2 )then
              LLA(NI) = LD
              XLA(NI) = XCOR(LLA(NI),5)
              YLA(NI) = YCOR(LLA(NI),5)
              ZLA(NI) = BELV(LD) + HP(LD)
              BEDGEMOVE = .TRUE.
              EXIT
            else
              if( NPD < 10000 .or. DEBUG )then
                PRINT '(A40,I8)','HYDRAULIC STRUCTURE, DRIFTER IS OUTSIDE:',NI
              else
                NSTOPPED(3) = NSTOPPED(3) + 1
              endif
              call SET_DRIFTER_OUT(XLA2, YLA2, ZLA2)
              EXIT
            endif
          endif
        endif
      enddo
    endif
    return
  endif
  
END SUBROUTINE CHECK_BCS

!---------------------------------------------------------------------------!
!> @details Sets drifter cell to '1' so that it is effectively invalid
!! and updates the x,y,z locations of the drifter
!> @param[in] XLA2
!> @param[in] YLA2
!> @param[in] ZLA2
!---------------------------------------------------------------------------!
SUBROUTINE SET_DRIFTER_OUT(XLA2,YLA2,ZLA2)

  implicit none

  real(RKD) ,intent(IN)   :: XLA2 !<
  real(RKD) ,intent(IN)   :: YLA2 !<
  real(RKD) ,intent(IN)   :: ZLA2 !<

  XLA(NI) = XLA2
  YLA(NI) = YLA2
  ZLA(NI) = ZLA2
  LLA(NI) = 1
  LLA_Global(NI) = 1

END SUBROUTINE SET_DRIFTER_OUT

!---------------------------------------------------------------------------!
!> @details Checks if layers are open better "open" faces and will allow drifters
!           to move between cells
!---------------------------------------------------------------------------!
LOGICAL FUNCTION ISOK(LL, LLO, K, L)

  implicit none
  
  integer, intent(IN) :: LL, LLO, K, L
  integer :: LLB
  logical :: FLAG
  
  if( NBLOCKED > 0 )then
    FLAG = .FALSE.
    do LLB = 1,NBLOCKED
      if( LBLOCKED(LLB) == L .or. LBLOCKED(LLB) == LLO )then
        if( LL == 1 )then
          ! *** NORTHWEST
          if( (SUB3D(LLO,K) > 0.5 .and. SVB3D(L,K) > 0.5) .or. (SVB3D(LADJ(2,LLO),K) > 0.5 .and. SUB3D(LADJ(2,LLO),K) > 0.5) )then
            FLAG = .TRUE.
          endif
        elseif( LL == 2 .and. SVB3D(LADJ(2,LLO),K) > 0.5 )then
          ! *** NORTH
            FLAG = .TRUE.
        elseif( LL == 3 )then  ! .and. (SVB3D(LADJ(2,LLO),K) > 0.5 .or. SUB3D(LADJ(6,LLO),K) > 0.5) )then
          ! *** NORTHEAST
          if( (SUB3D(LADJ(6,LLO),K) > 0.5 .and. SVB3D(LADJ(3,LLO),K) > 0.5) .or. (SVB3D(LADJ(2,LLO),K) > 0.5 .and. SUB3D(LADJ(3,LLO),K) > 0.5) )then
            FLAG = .TRUE.
          endif
        elseif( LL == 4 .and. SUB3D(LLO,K) > 0.5 )then
          ! *** WEST
            FLAG = .TRUE.
        elseif( LL == 6 .and. SUB3D(LADJ(6,LLO),K) > 0.5 )then
          ! *** EAST
            FLAG = .TRUE.
        elseif( LL == 7 )then  ! .and. (SVB3D(L,K) > 0.5 .or. SUB3D(L,K) > 0.5) )then
          ! *** SOUTHWEST
          if( (SUB3D(LLO,K) > 0.5 .and. SVB3D(LADJ(4,LLO),K) > 0.5) .or. (SVB3D(LLO,K) > 0.5 .and. SUB3D(LADJ(8,LLO),K) > 0.5) )then
            FLAG = .TRUE.
          endif
        elseif( LL == 8 .and. SVB3D(LLO,K) > 0.5 )then
          ! *** SOUTH
            FLAG = .TRUE.
        elseif( LL == 9 )then  ! .and. (SVB3D(L,K) > 0.5 .or. SUB3D(LE,K) > 0.5) )then
          ! *** SOUTHEAST
          if( (SUB3D(LADJ(6,LLO),K) > 0.5 .and. SVB3D(LADJ(6,LLO),K) > 0.5) .or. (SVB3D(LLO,K) > 0.5 .and. SUB3D(LADJ(9,LLO),K) > 0.5) )then
            FLAG = .TRUE.
          endif
        endif
        ISOK = FLAG
        return
      endif
    enddo
  endif
  
  ! *** ALLOW TRANSPORT FOR ALL LAYERS EXCEPT FOR DRIFTERS < KSZ(LLO)-2
  if( K >= KSZ(L)-1 )then
    ISOK = .TRUE.
  else
    ISOK = .FALSE.
  endif

END FUNCTION ISOK

END SUBROUTINE CONTAINER

!---------------------------------------------------------------------------!
!> @details CALCULATING VELOCITY COMPONENTS AT DRIFTER LOCATION BY USING
!! INVERSE DISTANCE POWER 2 INTERPOLATION FOR VELOCITY COMPONENTS
!! AT THE CENTROID OF THE CELLS
!
!> @param[in] LNI
!> @param[in] KNI
!> @param[in] NI
!> @param[out] U2NI
!> @param[out] V2NI
!> @param[out] W2NI
!---------------------------------------------------------------------------!

SUBROUTINE DRF_VELOCITY(LNI,KNI,NI,U2NI,V2NI,W2NI)

  integer,intent(IN )    :: LNI !< Cell the drifter belongs to
  integer,intent(IN )    :: KNI !< Layer the drifter belongs to
  integer,intent(IN )    :: NI  !< Drifter number
  real(RKD) ,intent(OUT) :: U2NI !< Interpoled U velocity
  real(RKD) ,intent(OUT) :: V2NI !< Interpoled V velocity
  real(RKD) ,intent(OUT) :: W2NI !< Interpoled W velocity
  ! *** Local variables
  integer :: L,LE,LN,K1,K2,K3,KZ1,KZ2,LL
  real(RKD) :: RAD2,SU1,SU2,SU3,SV1,SV2,SV3,SW1,SW2,SW3,XOUT,YOUT
  real(RKD) :: UTMPB,VTMPB,UTMPB1,VTMPB1,WTMPB,WTMPB1
  real(RKD) :: VELEK,VELNK,VELEK1,VELNK1,ZSIG,DZCTR,DZTOP
  real(RKD) :: UKB,UKT,VKB,VKT,UKB1,UKT1,VKB1,VKT1,VFACTOR,VE,UE

  ! *** ORDER OF CELLS
  ! ***   1  2  3
  ! ***   4  5  6
  ! ***   7  8  9
  VFACTOR = 0.5
  SU1 = 0
  SU2 = 0
  SU3 = 0
  SV1 = 0
  SV2 = 0
  SV3 = 0
  SW1 = 0
  SW2 = 0
  SW3 = 0

  ! *** UPDATE VERTICAL POSITION FOR TOP/BOTTOM LAYERS
  ZSIG = (ZLA(NI)-BELVLA(NI))/HPLA(NI)
  ZSIG = MAX(Z(LNI,0),MIN(Z(LNI,KC),ZSIG))
  ZLA(NI) = ZSIG*HPLA(NI)+BELVLA(NI)
  if( ZSIG >= ZZ(LNI,KNI) )then
    ! *** DRIFTER IS ABOVE THE LAYER MID DEPTH
    K1 = KNI
    K2 = MIN(KNI+1,KC)
    KZ1= KNI
    KZ2= KNI+1
  else
    ! *** DRIFTER IS BELOW THE LAYER MID DEPTH
    K1 = MAX(KSZ(LNI),KNI-1)
    K2 = KNI
    KZ1= KNI-1
    KZ2= KNI
  endif

  ! *** Go over the 8 cells around the active
  do LL = 1,9
    L = LADJ(LL,LNI)
    
    ! *** Drifter is active and within a valid cell
    if( L >= 2 .and. L <= LA )then
      !IF( KSZ(L) > KNI ) CYCLE
      LE = LEC(L)
      LN = LNC(L)

      ! *** CALCULATING HORIZONTAL VELOCITY COMPONENTS AT CENTROID-BOTTOM
      UKB  = 0.5*STCUV(L)*( RSSBCE(L)*U2(LE,K1) + RSSBCW(L)*U2(L,K1) )
      VKB  = 0.5*STCUV(L)*( RSSBCN(L)*V2(LN,K1) + RSSBCS(L)*V2(L,K1) )

      ! *** CALCULATING HORIZONTAL VELOCITY COMPONENTS AT CENTROID-TOP
      if( K2 < KSZ(L) )then
        K3 = MIN(K2+1,KSZ(L))      ! *** Use BOTTOM ACTIVE LAYER IF ONE LAYER HIGHER
      else
        K3 = K2
      endif
      UKT  = 0.5*STCUV(L)*( RSSBCE(L)*U2(LE,K3) + RSSBCW(L)*U2(L,K3) )
      VKT  = 0.5*STCUV(L)*( RSSBCN(L)*V2(LN,K3) + RSSBCS(L)*V2(L,K3) )

      ! *** DISTANCE FROM CENTROID SQUARED
      RAD2 = MAX((XLA(NI)-XCOR(L,5))**2 + (YLA(NI)-YCOR(L,5))**2,1D-8)

      DZCTR = ZCTR(L,KZ2)-ZCTR(L,KZ1)
      if( DZCTR > 1D-8 )then
        ! *** INTERPOLATION FOR HORIZONTAL VELOCITY COMPONENT
        UTMPB  = (UKT -UKB )*(ZSIG-ZCTR(L,KZ1))/DZCTR+UKB
        VTMPB  = (VKT -VKB )*(ZSIG-ZCTR(L,KZ1))/DZCTR+VKB
      else
        UTMPB  = UKT
        VTMPB  = VKT
      endif

      DZTOP = Z(L,KNI)-Z(L,KNI-1)
      if( DZTOP > 1D-8 )then
        ! *** INTERPOLATION FOR VERTICAL VELOCITY COMPONENT
        WTMPB  = (W2(L,KNI)-W2(L,KNI-1))*(ZSIG-Z(L,KNI-1))/DZTOP + W2(L,KNI-1)
      else
        WTMPB  = W2(L,KNI)
      endif

    else
      ! *** EDGE CELLS, APPLY A ZERO FACE VELOCITIES

      ! *** ORDER OF CELLS
      ! ***   1  2  3
      ! ***   4  5  6
      ! ***   7  8  9

      UTMPB = SLIPFACTOR*U2(LNI,KNI)
      VTMPB = SLIPFACTOR*V2(LNI,KNI)
      WTMPB = SLIPFACTOR*W2(LNI,KNI)

      ! *** PMC NOTE - THESE COORDINATES COULD BE CALCULATED ONCE AND SAVED FOR LATER USE
      if( LL == 1 )then
        ! *** NORTHWEST
        XOUT = XCOR(LNI,5) + ( DYP(LNI)*CVE(LNI) - DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) + ( DYP(LNI)*CVN(LNI) - DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      elseif( LL == 2 )then
        ! *** NORTH FACE (MOVE WITH DRIFTER)
        XOUT = XLA(NI) + 0.5*( DYP(LNI)*CVE(LNI) )
        YOUT = YLA(NI) + 0.5*( DYP(LNI)*CVN(LNI) )
        RAD2 = MAX( (XLA(NI)-XOUT)**2 + (YLA(NI)-YOUT)**2, 1D-8)
        VTMPB = -VFACTOR*MAX(V2(LNI,KNI),0.0)
      elseif( LL == 3 )then
        ! *** NORTHEAST
        XOUT = XCOR(LNI,5) + ( DYP(LNI)*CVE(LNI) + DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) + ( DYP(LNI)*CVN(LNI) + DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      elseif( LL == 4 )then
        ! *** WEST FACE (MOVE WITH DRIFTER)
        XOUT = XLA(NI) - 0.5*( DXP(LNI)*CUE(LNI) )
        YOUT = YLA(NI) - 0.5*( DXP(LNI)*CUN(LNI) )
        RAD2 = MAX( (XLA(NI)-XOUT)**2 + (YLA(NI)-YOUT)**2, 1D-8)
        UTMPB = -VFACTOR*MIN(U2(LEC(LNI),KNI),0.0)
      elseif( LL == 6 )then
        ! *** EAST FACE (MOVE WITH DRIFTER)
        XOUT = XLA(NI) + 0.5*( DXP(LNI)*CUE(LNI) )
        YOUT = YLA(NI) + 0.5*( DXP(LNI)*CUN(LNI) )
        RAD2 = MAX( (XLA(NI)-XOUT)**2 + (YLA(NI)-YOUT)**2, 1D-8)
        UTMPB = -VFACTOR*MAX(U2(LNI,KNI),0.0)
      elseif( LL == 7 )then
        ! *** SOUTHWEST
        XOUT = XCOR(LNI,5) - ( DYP(LNI)*CVE(LNI) + DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) - ( DYP(LNI)*CVN(LNI) + DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      elseif( LL == 8 )then
        ! *** SOUTH  (MOVE WITH DRIFTER)
        !XOUT = XCOR(LNI,5) - ( DYP(LNI)*CVE(LNI) )
        !YOUT = YCOR(LNI,5) - ( DYP(LNI)*CVN(LNI) )
        XOUT = XLA(NI) - 0.5*( DYP(LNI)*CVE(LNI) )
        YOUT = YLA(NI) - 0.5*( DYP(LNI)*CVN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
        VTMPB  = -VFACTOR*MIN(V2(LNC(LNI),KNI),0.0)
      elseif( LL == 9 )then
        ! *** SOUTHEAST
        XOUT = XCOR(LNI,5) - ( DYP(LNI)*CVE(LNI) - DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) - ( DYP(LNI)*CVN(LNI) - DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      endif
      L = LNI
    endif

    ! *** ROTATION
    VELEK  = CUE(L)*UTMPB  + CVE(L)*VTMPB
    VELNK  = CUN(L)*UTMPB  + CVN(L)*VTMPB

    ! *** APPLY WEIGHTING
    SU2 = SU2 + VELEK/RAD2
    SU3 = SU3 + 1._8/RAD2
    SV2 = SV2 + VELNK/RAD2
    SV3 = SV3 + 1._8/RAD2
    SW2 = SW2 + WTMPB/RAD2
    SW3 = SW3 + 1._8/RAD2

  enddo

  U2NI = SU2/SU3
  V2NI = SV2/SV3
  W2NI = SW2/SW3

END SUBROUTINE DRF_VELOCITY

SUBROUTINE RANDCAL(L,K,NP)
  integer,intent(IN) :: L,K,NP
  !REAL(8),external :: DRAND -Return double-precision random numbers in the range 0.0 to 1.0, not including the end points.
  real(RKD) :: COEF

  ! *** HORIZONTAL
  if( LA_PRAN == 1 .or. LA_PRAN == 3 )then
    if( LA_DIFOP == 0 )then
      COEF = SQRT(2.*AH(L,K)*DELTD)*LA_HORDIF
    else
      COEF = DIFFH
    endif
#ifdef GNU  
    XLA(NP) = XLA(NP) + (2.*RAND(0)-1.)*COEF
    YLA(NP) = YLA(NP) + (2.*RAND(0)-1.)*COEF
#else
    XLA(NP) = XLA(NP) + (2.*DRAND(0)-1.)*COEF
    YLA(NP) = YLA(NP) + (2.*DRAND(0)-1.)*COEF
#endif    
  endif

  ! *** VERTICAL (IF LA_ZCAL = 1 THEN FULLY 3D)
  if( LA_PRAN >= 2 .and. LA_ZCAL == 1 )then
    if( LA_DIFOP == 0 )then
      COEF = SQRT(2.*AV(L,K)*DELTD)*LA_VERDIF
    else
      COEF = DIFFV
    endif
#ifdef GNU  
    ZLA(NP) = ZLA(NP) + (2*RAND(0)-1)*COEF
#else
    ZLA(NP) = ZLA(NP) + (2*DRAND(0)-1)*COEF
#endif 
  endif
  
END SUBROUTINE

! **** APPLY WIND SHEAR FOR DRIFTERS AT THE SURFACE
SUBROUTINE WIND_DRIFT(L,K,NP)

  integer,intent(IN) :: L,K,NP
  real(RKD) :: WSPD, COEF, XCOEF, YCOEF
  
  if( IOSWD == 1 )then
    ! *** APPLY SHEAR TO THE SURFACE OIL BASED ON 
    WSPD = SQRT(WNDVELE(L)**2 + WNDVELN(L)**2)       ! *** 10m winds
    COEF = OSWDA * WSPD + OSWDB 
  
    XLA(NP) = XLA(NP) + COEF * WNDVELE(L) * DELTD 
    YLA(NP) = YLA(NP) + COEF * WNDVELN(L) * DELTD
    
  elseif( IOSWD == 2 )then
    XCOEF = OSWDA * WNDVELE(L) + OSWDB
    YCOEF = OSWDA * WNDVELN(L) + OSWDB
  
    XLA(NP) = XLA(NP) + XCOEF * WNDVELE(L) * DELTD
    YLA(NP) = YLA(NP) + YCOEF * WNDVELN(L) * DELTD
  endif
  
END SUBROUTINE WIND_DRIFT

!---------------------------------------------------------------------------!
!> @details INTERPOLATION OF THE TOTAL WATER DEPTH AND BOTTOM ELEVATION
!! FOR THE DRIFTER NI AT EACH TIME INSTANT AND EACH LOCATION
!
!> @param[in]  LNI
!> @param[in]  NI
!> @param[out] BELVNI
!> @param[out] HPNI
!---------------------------------------------------------------------------!
SUBROUTINE DRF_DEPTH(LNI, NI, BELVNI, HPNI)

  implicit none

  integer,intent(IN) :: LNI        !< Cell index 'L' for the drifter
  integer,intent(IN) :: NI         !< drifter number
  real(RKD),intent(OUT) :: BELVNI  !< Interpolated bed elevation
  real(RKD),intent(OUT) :: HPNI    !< Interpolated water surface elevation

  ! *** Local variables
  integer :: L,LL
  real(RKD) :: BELVNI1,BELVNI2,RAD2,ZETA

  BELVNI1 = 0.
  BELVNI2 = 0.
  ZETA = 0.

  do LL = 1,9
    L = LADJ(LL,LNI)
    if( L >= 2 .and. L <= LA )then
      RAD2 = MAX((XLA(NI)-XCOR(L,5))**2 + (YLA(NI)-YCOR(L,5))**2,1D-8)         ! *** INVERSE DISTANCE SQUARED
      !RAD2 = SQRT(MAX((XLA(NI)-XCOR(L,5))**2 + (YLA(NI)-YCOR(L,5))**2,1D-8))  ! *** LINEAR INTERPOLATION
      BELVNI1 = BELVNI1 + BELV(L)/RAD2
      BELVNI2 = BELVNI2 + 1._8/RAD2
      ZETA = ZETA + (HP(L)+BELV(L))/RAD2
    endif
  enddo
  BELVNI = BELVNI1/BELVNI2        ! *** INTERPOLATED BELV
  ZETA = ZETA/BELVNI2             ! *** INTERPOLATED WSEL
  HPNI = MAX(ZETA - BELVNI,HDRY)

END SUBROUTINE DRF_DEPTH

!---------------------------------------------------------------------------!
!> @details Determines the drifters vertical location and layer it belongs in
!! based on interpolated values
!
!> @param[in]  LNI
!> @param[in]  BELVNI
!> @param[in]  HPNI
!> @param[out] KLN
!> @param[inout] ZLN
!---------------------------------------------------------------------------!
SUBROUTINE DRF_LAYER(LNI,BELVNI,HPNI,KLN,ZLN)

  implicit none
  ! *** Dummy variables
  integer,intent(IN)    :: LNI    !< Current cell 'L' index the drifter belongs to
  real(RKD), intent(IN) :: BELVNI !< Interpolated bed elevation
  real(RKD), intent(IN) :: HPNI   !< Interpolated water surface elevation
  integer,intent(OUT)   :: KLN    !< Layer this drifter belongs in
  real(RKD), intent(INOUT) :: ZLN !< Vertical location of the drifter based on interpolated values

  ! *** Local variables
  integer :: K
  real(RKD) :: ZSIG

  ! *** Check if the bottom most active layer (for the current cell) is equal to the top layer
  if( KSZ(LNI) == KC )then  ! *** Alberta
    KLN = KC
    ZSIG = (ZLN-BELVNI)/HPNI
    ZSIG = MAX(Z(LNI,0),MIN(Z(LNI,KC),ZSIG))
    ZLN  = ZSIG*HPNI + BELVNI
  elseif( LNI >= 2 .and. KC > 1 )then
    ZSIG = (ZLN - BELVNI)/HPNI
    ZSIG = MAX( Z(LNI,0), MIN(Z(LNI,KC),ZSIG) )
    ZLN = ZSIG*HPNI + BELVNI
    ! *** From the current active layer for cell the drifter belongs to the top layer
    do K = KSZ(LNI),KC
      ! *** Based on the new vertical height determine the layer the drifter belongs in
      if( ZSIG  >= Z(LNI,K-1) .and. ZSIG  <= Z(LNI,K) )then
        KLN = K
        EXIT
      endif
    enddo
  endif

END SUBROUTINE DRF_LAYER

!---------------------------------------------------------------------------!
!> @details CALCULATING GRADIENT COMPONENTS OF DIFFUSION COEFFICIENT
!! AT DRIFTER LOCATION BY BACKWARD SPACE METHOD
!
!> @param[in]   LNI
!> @param[in]   KNI
!> @param[in]   NI
!> @param[out]  DAHX
!> @param[out]  DAHY
!> @param[out]  DAVZ

!---------------------------------------------------------------------------!

SUBROUTINE DIFGRAD(LNI, KNI, NI, DAHX, DAHY, DAVZ)   ! **************************

  integer,intent(IN )    :: LNI !< Cell index the drifter is in
  integer,intent(IN )    :: KNI !< Layer index the drifer is in
  integer,intent(IN )    :: NI !< Drifter number
  real(RKD) ,intent(OUT) :: DAHX !< Horizontal diffusion coeffecient in x-direction
  real(RKD) ,intent(OUT) :: DAHY !< horizontal diffusion coeffecient in y-direction
  real(RKD) ,intent(OUT) :: DAVZ !< horizontal diffusion coeffecient in z-direction

  ! *** Local variables
  real(RKD) :: DAHXC, DAHYC, ZSIG, WEIC, DAVB, DAVT

  DAHXC = SUB3D(LNI,KNI)*( AH(LNI,KNI) - AH(LWC(LNI),KNI) )/DXU(LNI)
  DAHYC = SVB3D(LNI,KNI)*( AH(LNI,KNI) - AH(LSC(LNI),KNI) )/DYV(LNI)

  ! ***  HORIZONTAL DIFFUSION COMPONENT (M/S)
  DAHX = CUE(LNI)*DAHXC + CVE(LNI)*DAHYC
  DAHY = CUN(LNI)*DAHXC + CVN(LNI)*DAHYC

  ! *** APPLY A FACTOR TO GET THE HORIZONTAL DIFFUSION COMPONENT (M/S)
  DAHX = DAHX*LA_HORDIF
  DAHY = DAHY*LA_HORDIF

  if( KC <= 2 .or. LA_ZCAL == 0 )then
    DAVZ = 0.
  else
    ! *** AV IS THE VERTICAL DIFFUSION COEFFICIENT AT THE TOP OF THE LAYER
    ! *** AV(L,KC) = 0 AND AV(L,0) = 0
    ZSIG = (ZLA(NI)-BELVLA(NI))/HPLA(NI)
    ZSIG = MAX(Z(LNI,0),MIN(Z(LNI,KC),ZSIG))
    ZLA(NI) = ZSIG*HPLA(NI)+BELVLA(NI)
    WEIC = Z(LNI,KC)-(Z(LNI,KNI)-ZSIG)*DZIC(LNI,KNI)  !WEIGHTING COEFFICIENT

    if( KNI == KC )then
      DAVT = 0.
      DAVB = (AV(LNI,KS)-AV(LNI,KS-1))/HPK(LNI,KNI)
    elseif( KNI == 1 )then
      DAVT = (AV(LNI,KNI+1)-AV(LNI,KNI))/HPK(LNI,KNI)
      DAVB = 0.
    else
      DAVT = 2._RKD*(AV(LNI,KNI+1)-AV(LNI,KNI  ))/(DZC(LNI,KNI+1)+DZC(LNI,KNI  ))/HP(LNI)
      DAVB = 2._RKD*(AV(LNI,KNI  )-AV(LNI,KNI-1))/(DZC(LNI,KNI  )+DZC(LNI,KNI-1))/HP(LNI)
    endif
    DAVZ = (DAVT-DAVB)*WEIC+DAVB

    ! *** APPLY A VERTICAL GRADIENT FACTOR
    DAVZ = DAVZ*LA_VERDIF

  endif

END SUBROUTINE DIFGRAD

SUBROUTINE OIL_PROC(NP)
  ! *** CALCULATING  OIL REDUCTION DUE TO:
  ! *** EVAPORATION AND BIODEGRADATION
  ! *** OUTPUT: DVOL(NP)
  integer,intent(IN ) :: NP
  integer :: NG,LNP
  real(RKD) :: DFV,RKA,TK,RGAS,DTHE
  real(RKD) :: HENRY,OILMAS,BIODEG,VWIND,TDEP,TNP
  DATA RGAS /8.314_RKD/

  ! *** OIL EVAPORATION
  ! *** BASED ON THE PAPER:
  ! *** EVAPORATION RATE OF SPILLS OF HYDROCARBONS AND PETROLEUM MIXTURES
  ! *** BY WARREN STLVER AND DONALD MACKAY
  ! *** VAPOR PRESSURE IS REFERED FROM: 2013 CRUDE CHARACTERISTICS
  ! *** MOLAR VOLUME IS REFERED FROM: PhD THESIS BY NA DU
  !    INVESTIGATION OF HYDROGENATED AND FLUORINATED SURFACTANT
  !    BASED-SYSTEMS FOR THE DESIGN OF POROUS SILICA MATERIALS

  NG  = LA_GRP(NP)
  LNP = LLA(NP)

  ! *** LIMIT EVAPORATION, ONLY WHEN DRIFTER IS NEAR SURFACE
  if( DLA(NP) < 0.1 .and. DVAP(NG) > 0. )then
    if( ICECELL(LNP) )then
      VWIND = 0.001
    else
      VWIND  = WINDST(LNP)
      VWIND = MAX(VWIND,0.001)
    endif
    RKA = 0.00125_8*VWIND                ! *** MODIFIED MASS TRANSFER COEFFICENT: RKA = 0.00125_8*VWIND [M/S]
    DTHE = RKA*DARE(NG)*DELTD/DVOL0(NG)  ! *** EVAPORATIVE EXPOSURE FOR A DRIFTER (DIMENSIONLESS)
    if( ISTRAN(2) > 0 )then
      TK = TEM(LNP,KC)+273.15
    else
      TK = DTEM(NG)+273.15
    endif
    HENRY = DVAP(NG)*DVMO(NG)/(RGAS*TK)  ! *** HENRY'S LAW: DIMENSIONLESS, DVMO = 0.000194-0.000293 M3/MOL
    DFV = HENRY*DTHE                     ! *** STIVER AND MACKAY,         DFV: RATIO OF VAPOR VOLUME AND INITIAL VOLUME
    ! *** 1-Fv = (1-Fm)(RHOo/RHOf)   DFV ASSUMES NO CHANGE IN DENSITY (TO BE UPDATED LATER)
    DVOL(NP) = DVOL(NP)*(1.-DFV)
  endif

  ! *** OIL BIODEGRADATION
  if( DARE(NG)>0. )then
    ! *** DETERMINE TEMPERATURE DEPENDENCY
    if( ISTRAN(2) > 0 )then
      TNP = TEM(LNP,KLA(NP))-DTEM(NG)
      TDEP = EXP(-0.01*TNP*TNP)
    else
      TDEP = 1.0
    endif
    BIODEG = EXP(-TDEP*DRAT(NG)*DELTD)                     ! *** FRACTION REMAINING AT END OF TIME STEP
    OILMAS = DVOL(NP)*DDEN(NG)*BIODEG                      ! *** MASS OF AN OIL DRIFTER IN KG
    DVOL(NP) = OILMAS/DDEN(NG)                             ! *** DVOL UPDATE
  endif


END SUBROUTINE

SUBROUTINE INIT_OIL
  integer :: NG,NP
  integer :: NUMD(NGRP)
  real(RKD)   :: XMAXG,XMING,YMAXG,YMING,RAD
  real(RKD)   :: GARE(NGRP)

  ! *** CALCULATING SURFACE AREA OF OIL SPILL FOR EACH GROUP
  GARE = 0  ! TOTAL SURFACE AREA OF OIL FOR EACH GROUP
  DARE = 0  ! SURFACE AREA OF OIL DRIFTER FOR EACH GROUP
  NUMD = 0  ! NUMBER OF OIL DRIFTERS FOR EACH GROUP

  do NG = 1,NGRP
    if( ISOILSPI(NG) == 0 ) CYCLE

    XMING = MINVAL(XLA,LA_GRP == NG)
    XMAXG = MAXVAL(XLA,LA_GRP == NG)
    YMING = MINVAL(YLA,LA_GRP == NG)
    YMAXG = MAXVAL(YLA,LA_GRP == NG)
    RAD   = MAX(0.25*(XMAXG-XMING+YMAXG-YMING),1.0)
    NUMD(NG) = SUM(LA_GRP,LA_GRP == NG)/NG
    GARE(NG) = PI*RAD**2
    DARE(NG) = GARE(NG)/NUMD(NG)
    DVOL0(NG) = GVOL(NG)/NUMD(NG)
  enddo

  ! *** OIL: NPD,NGRP,LA_GRP(NP)
  allocate(DVOL(NPD))
  DVOL = 0
  do NP = 1,NPD
    if( ISOILSPI(LA_GRP(NP)) == 0 ) CYCLE
    DVOL(NP) = DVOL0(LA_GRP(NP))           !INITIAL VOLUME OF A DRIFTER

  enddo

END SUBROUTINE

!---------------------------------------------------------------------------!
!> @details Checks if a drifter has moved into the ghost cell region of a
!! given subdomain
!
!> @author Zander
!
!---------------------------------------------------------------------------!
Subroutine Check_Drifter_In_Ghost(ID)

  implicit none

  ! *** Dummy variables
  integer, intent(in) :: ID

  ! *** Local variables
  integer :: drifter_local_l !< Cell L index the drifter is currently in
  integer :: local_i, local_j, index, LG
  integer :: process_id_to_send

  ! *** Checks if cell id is valid and active
  if( LLA(ID) > 1 .and. LLA(ID) < LC .and. JSPD(ID) == 0 )then
    ! *** get the local drifter
    drifter_local_l = LLA(ID)
    
    ! *** Get the local i/j from the L index ofthe current drifter
    local_i = IL(drifter_local_l)  ! Map2Global(drifter_local_l).IL
    local_j = JL(drifter_local_l)  ! Map2Global(drifter_local_l).JL

    ! *** Test if in ghost cells
    if( local_i < 3 .or. local_i > IC-2 .or. local_j < 3 .or. local_j > jc-2 )then
      ! *** Drifter is in a ghost cell
      
      ! *** Get the global L value
      LG = Map2Global(drifter_local_l).LG
      process_id_to_send = sorted_loc_to_glob(LG).process
      
      ! *** Determine where this process is in reference to the current one
      if( nbr_west == process_id_to_send )then
        ! *** Increment the number of drifters to commmunicate
        num_drifters_send_west = num_drifters_send_west + 1

        drifter_ids_send_west(num_drifters_send_west) = ID

      elseif( nbr_east == process_id_to_send )then
        ! *** Increment the number of drifters to commmunicate
        num_drifters_send_east = num_drifters_send_east + 1
        
        drifter_ids_send_east(num_drifters_send_east) = ID

      elseif( nbr_north == process_id_to_send )then
        ! *** Increment the number of drifters to commmunicate
        num_drifters_send_north = num_drifters_send_north + 1

        drifter_ids_send_north(num_drifters_send_north) = ID

      elseif( nbr_south == process_id_to_send )then
        ! *** Increment the number of drifters to commmunicate
        num_drifters_send_south = num_drifters_send_south + 1

        drifter_ids_send_south(num_drifters_send_south) = ID
      endif
    endif  
  endif

End Subroutine Check_Drifter_In_Ghost

SUBROUTINE Gather_Drifter_Arrays

  !---------------------------------------------------------------------------!
  ! *** Gathers values for all of Drifter variables needed for EE linkage
  !---------------------------------------------------------------------------!
  use Mod_Assign_Loc_Glob_For_Write
  use Mod_Map_Write_EE_Binary

  implicit none

  ! *** Local variables
  integer :: ierr, ii, LG, NP, global_id

  ! *** Global arrays before they are sorted
  integer :: NPD_Write
  integer,   Save, Allocatable :: map_id(:)
  integer,   Save, Allocatable :: LLA_Unsorted_Global(:)
  real(RKD), Save, Allocatable :: XLA_Unsorted_Global(:)
  real(RKD), Save, Allocatable :: YLA_Unsorted_Global(:)
  real(RKD), Save, Allocatable :: ZLA_Unsorted_Global(:)
  logical,   Save, Allocatable :: mask(:)                         !< mask used to exclude deactivated cells or ones that have been communicated to another process

  ! *** Arrays to be allocated and initialized during the "PACK" process
  integer,   Allocatable :: map_tmp(:)
  integer,	 Allocatable :: lla_tmp(:)
  real(RKD), Allocatable :: xla_tmp(:)
  real(RKD), Allocatable :: yla_tmp(:)
  real(RKD), Allocatable :: zla_tmp(:)

  ! *** FIRST CALL INITIALIZATION  
  if( .not. allocated(map_id) )then
    ! *** Allocate Arrays
    allocate(gl_drift_mapping(NPD))

    allocate(map_id(NPD))
    allocate(LLA_Unsorted_Global(NPD))
    allocate(XLA_Unsorted_Global(NPD))
    allocate(YLA_Unsorted_Global(NPD))
    allocate(ZLA_Unsorted_Global(NPD))
    allocate(mask(NPD))
    
    !Allocate(map_id(NPD))
    !Allocate(lla_tmp(NPD))
    !Allocate(xla_tmp(NPD))
    !Allocate(yla_tmp(NPD))
    !Allocate(zla_tmp(NPD))

    ! *** Zero Arrays
    gl_drift_mapping = 0

    map_id              = 0
    LLA_Unsorted_Global = 0
    XLA_Unsorted_Global = 0.0
    YLA_Unsorted_Global = 0.0
    ZLA_Unsorted_Global = 0.0
    mask                = .false.

    !map_tmp = 0
    !lla_tmp = 0
    !xla_tmp = 0.0
    !yla_tmp = 0.0
    !zla_tmp = 0.0
  endif
  
  ! *** Convert local LLA to have global L values
  Do NP = 1,NPD
    map_id(NP) = 0
    gl_drift_mapping(NP) = 0
    
    ! *** Check if this drifter belongs to the current subdomain
    if( LLA_Process(NP) /= process_id )then
      mask(NP) = .FALSE.
      cycle
    endif
  
    ! *** Check for inactive drifters
    if( jspd(NP) == 1 )then
      mask(NP) = .FALSE.
      cycle
    endif

    ! *** Check for stopped drifters
    if( jspd(NP) == 0 .and. LLA_Global(NP) < 1 )then
      mask(NP) = .FALSE.
      cycle
    endif

    ! *** Setup mask
    mask(NP) = .TRUE.
    map_id(NP) = NP
  enddo

  ! *** This removes parts of each local array that does not have any drifter info in it.  Makes a coniguous array to make the gathering easier
  lla_tmp = PACK(LLA_Global, mask)
  xla_tmp = PACK(xla, mask)
  yla_tmp = PACK(yla, mask)
  zla_tmp = PACK(zla, mask)
  map_tmp = PACK(map_id, mask)

  ! *** Gather up LLA, use generic Gather_Soln routine so we can get the mapping
  call Gather_Drifter(size(lla_tmp), lla_tmp, NPD, LLA_Unsorted_Global)
  
  ! *** Need to get the mapping so we can put drifters in the correct global position
  call Gather_Drifter(size(map_tmp), map_tmp, NPD, gl_drift_mapping)
  
  ! *** Gather just the XLA array
  call Gather_Drifter(size(xla_tmp), xla_tmp, NPD, XLA_Unsorted_Global)
  
  ! *** Gather just the YLA array
  call Gather_Drifter(size(yla_tmp), yla_tmp, NPD, YLA_Unsorted_Global)
  
  ! *** Gather just the ZLA array
  call Gather_Drifter(size(zla_tmp), zla_tmp, NPD, ZLA_Unsorted_Global)

  if( process_id == master_id )then
    ! *** Sort so the global indexing matches when writing out
    do NP = 1,NPD

      ! *** Get the correct global drifter id
      global_id = gl_drift_mapping(NP)

      if( global_id > 0 .and. global_id <= NPD )then
        LLA_Global(global_id) =  LLA_Unsorted_Global(NP)
        
        XLA(global_id) =  XLA_Unsorted_Global(NP)
        YLA(global_id) =  YLA_Unsorted_Global(NP)
        ZLA(global_id) =  ZLA_Unsorted_Global(NP)
      !else
      !  LLA(NP) = 1
      !  XLA(NP) = 0.0
      !  YLA(NP) = 0.0
      !  ZLA(NP) = 0.0
      endif
    enddo
  endif

  ! *** Deallocate to save memory as this is used infrequently
  deallocate(lla_tmp)
  deallocate(xla_tmp)
  deallocate(yla_tmp)
  deallocate(zla_tmp)
  deallocate(map_tmp)
  
END SUBROUTINE Gather_Drifter_Arrays

! *** Binary search of a sorted list
INTEGER FUNCTION Search_List(NN, NLIST, iTarget)

  integer, intent(IN) :: NN, iTarget
  integer, intent(IN) :: NLIST(NPD)
  
  ! *** Local Variables
  integer :: I1, I2, I3, nnn
  
  nnn = 0
  I1 = 1
  I3 = NN 
  do while (I1 <= I3)
    nnn = nnn + 1
    I2 = (I1 + I3) / 2
    
    if( NLIST(I2) > iTarget )then
      I3 = I2 - 1
    elseif( NLIST(I2) < iTarget )then
      I1 = I2 + 1
    else
      ! *** Found the target
      Search_List = I2
      return
    endif
  enddo
  Search_List = -1
  
End FUNCTION Search_List


! *** Insert a value into a sorted list
SUBROUTINE Insert_Value(NN, NLIST, iTarget)

  integer, intent(IN)    :: iTarget
  integer, intent(INOUT) :: NN, NLIST(NPD)
  
  ! *** Local Variables
  integer :: I1, I2, I3
  
  if( NN == NPD ) return
  
  ! *** Find the insertion point
  I2 = 1
  do while (NLIST(I2) < iTarget)
    I2 = I2 + 1
    if( I2 > NN )then
      ! *** Add at the end
      NN = NN + 1
      NLIST(NN) = iTarget
      return
    endif
  enddo
  
  ! *** Shift the remaining values and insert new value
  I1 = I2
  do I3 = NN,I1,-1
    NLIST(I3+1) = NLIST(I3)
  enddo
  NN = NN + 1
  NLIST(I2) = iTarget
  
End SUBROUTINE Insert_Value
    subroutine create_nc_lpt()
        real(4), parameter :: MISSING_VALUE = -999.
        integer(4), parameter :: TIMECNT = NF90_UNLIMITED
        character(80) :: filename
        character(20) :: bdate
        integer :: status

        if( IS_NC_OUT(12) == 0 ) return

        lpt_time_idx = 0
    
        BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)
        filename  = OUTDIR//'efdc_drifters.nc'
        
            ! *** ENTER DEFINE MODE
        status = NF90_CREATE(filename, NF90_CLOBBER+NF90_NETCDF4, nc_lpt(0))    !NF90_HDF5
        if(status /= NF90_NOERR )then
            return
        endif
        status = nf90_def_dim(nc_lpt(0), 'TRACKS', NPD, lpt_npd_dim)
        status = nf90_def_dim(nc_lpt(0), 'TIME', TIMECNT, lpt_time_dim)
        
        status = nf90_def_var(nc_lpt(0),'time',NF90_DOUBLE,(/ lpt_time_dim /), nc_lpt(1))
        status = nf90_put_att(nc_lpt(0), nc_lpt(1),'standard_name','time')
        status = nf90_put_att(nc_lpt(0),nc_lpt(1),'long_name','time')
        status = nf90_put_att(nc_lpt(0),nc_lpt(1),'units','days since '//trim(bdate))
        status = nf90_put_att(nc_lpt(0),nc_lpt(1),'calendar','julian')
        status = nf90_put_att(nc_lpt(0),nc_lpt(1),'axis','T')
        status = nf90_put_att(nc_lpt(0),nc_lpt(1),'format','modified julian day (MJD)')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(1),'time_zone','UTC')            

        status = nf90_def_var(nc_lpt(0),'lat',nf90_float,(/lpt_npd_dim,lpt_time_dim/),nc_lpt(3))
        status = nf90_def_var_deflate(nc_lpt(0), nc_lpt(3), 1, 1, deflev)
        status = nf90_put_att(nc_lpt(0),nc_lpt(3),'standard_name','latitude')
        status = nf90_put_att(nc_lpt(0),nc_lpt(3),'long_name','drifter latitude')
        status = nf90_put_att(nc_lpt(0),nc_lpt(3),'units','degrees_north')
        status = nf90_put_att(nc_lpt(0),nc_lpt(3),'valid_min','-90')
        status = nf90_put_att(nc_lpt(0),nc_lpt(3),'valid_max','90')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(3),'grid_mapping','crs')
        status = nf90_put_att(nc_lpt(0),nc_lpt(3),'_FillValue',missing_value)

        status = nf90_def_var(nc_lpt(0),'lon',nf90_float,(/lpt_npd_dim,lpt_time_dim/),nc_lpt(2))
        status = nf90_def_var_deflate(nc_lpt(0), nc_lpt(2), 1, 1, deflev)
        status = nf90_put_att(nc_lpt(0),nc_lpt(2),'standard_name','longitude')
        status = nf90_put_att(nc_lpt(0),nc_lpt(2),'long_name','drifter longitude')
        status = nf90_put_att(nc_lpt(0),nc_lpt(2),'units','degrees_east')
        status = nf90_put_att(nc_lpt(0),nc_lpt(2),'valid_min','-180')
        status = nf90_put_att(nc_lpt(0),nc_lpt(2),'valid_max','180')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(2),'grid_mapping','crs')
        status = nf90_put_att(nc_lpt(0),nc_lpt(2),'_FillValue',missing_value)

        status = nf90_def_var(nc_lpt(0),'elev',nf90_float,(/lpt_npd_dim,lpt_time_dim/), nc_lpt(4))
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'standard_name','altitude')
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'long_name','drifter elevation')
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'units','m')
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'positive','up')
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'axis','Z')
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'coordinates','time lon lat')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'grid_mapping','crs')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'_CoordinateZisPositive','up')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'_CoordinateTransformType','Vertical')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'_CoordinateAxisType','GeoZ')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'_CoordinateAxes','lev')
        status = nf90_put_att(nc_lpt(0),nc_lpt(4),'_FillValue',MISSING_VALUE)

        if( any(isoilspi == 1) )then
            status = nf90_def_var(nc_lpt(0),'oil_vol',nf90_float,(/lpt_npd_dim,lpt_time_dim/), nc_lpt(5))
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'standard_name','')
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'long_name','drifter oil volume')
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'units','m3')
            !status = nf90_put_att(nc_lpt(0),nc_lpt(5),'grid_mapping','crs')
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'coordinates','time lon lat elev')
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'_FillValue',MISSING_VALUE)

            status = nf90_def_var(nc_lpt(0),'oil_mass',nf90_float,(/lpt_npd_dim,lpt_time_dim/), nc_lpt(6))
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'standard_name','')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'long_name','drifter oil mass')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'units','mg/L')
            !status = nf90_put_att(nc_lpt(0),nc_lpt(6),'grid_mapping','crs')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'coordinates','time lon lat elev')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'_FillValue',MISSING_VALUE)
        endif
        
        ! *** ASSIGN GLOBAL ATTRIBUTES
        status = nf90_put_att(nc_lpt(0),nf90_global,'Conventions','CF-1.4')
        status = nf90_put_att(nc_lpt(0),NF90_GLOBAL,'Base_date',BDATE)

        ! *** LEAVE DEFINE MODE
        status = nf90_enddef(nc_lpt(0))
    end subroutine 
    
    subroutine close_nc_lpt()
        integer :: status
        if( IS_NC_OUT(12) == 0 ) return
        status = nf90_close(nc_lpt(0))
        if(status /= NF90_NOERR )then
            print * ,'cannot close efdc_drifters.nc!'
            return
        endif
    end subroutine
    
    subroutine write_nc_lpt()
        real(4), parameter :: MISSING_VALUE = -999.
        real(8) :: ti(1), xm(1),ym(1),xll(1),yll(1)
        real(4), allocatable :: array2d(:,:)
        integer :: i,status,str1d(2),cnt1d(2)

        if( IS_NC_OUT(12) == 0 ) return
        allocate(array2d(2,NPD))
        array2d = MISSING_VALUE
        do i = 1,NPD
            if( LLA_Global(i) < 2 .or. LLA_Global(i) > LA_Global) cycle ! *** Modified for global
            if( ISWGS84 == 0 .and. UTMZ /= 0 )then
                xm(1) = XLA(i)
                ym(1) = YLA(i)
                call utmr_wgs84(xm,ym,xll,yll)
                array2d(1,i) = xll(1)
                array2d(2,i) = yll(1)
            else
                array2d(1,i) = XLA(i)
                array2d(2,i) = YLA(i)
            endif
        enddo
        if(process_id /= master_id ) return

        lpt_time_idx = lpt_time_idx + 1
        ti = TIMEDAY
        str1d = (/1, lpt_time_idx/)
        cnt1d = (/NPD, 1/)

        status = nf90_put_var(nc_lpt(0), nc_lpt(1), TI, (/lpt_time_idx/), (/1/))
        status = nf90_put_var(nc_lpt(0), nc_lpt(2), array2d(1,:), str1d, cnt1d)      ! XLA
        status = nf90_put_var(nc_lpt(0), nc_lpt(3), array2d(2,:), str1d, cnt1d)      ! YLA
        status = nf90_put_var(nc_lpt(0), nc_lpt(4), real(ZLA), str1d, cnt1d)         ! ZLA 
        if( any(isoilspi == 1) )then
            status = nf90_put_var(nc_lpt(0), nc_lpt(5), real(DVOL), str1d, cnt1d)    ! DVOL
            status = nf90_put_var(nc_lpt(0), nc_lpt(6), real(MOC), str1d, cnt1d)     ! MOC
        endif
        deallocate(array2d)
    end subroutine
END MODULE
