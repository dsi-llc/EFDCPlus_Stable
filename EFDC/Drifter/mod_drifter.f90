! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

  USE GLOBAL
  USE OMP_LIB
  USE XYIJCONV
  USE INFOMOD,ONLY:SKIPCOM
  USE IFPORT

  Use MPI
#ifdef _MPI
  Use Variables_MPI_Drifter
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  Use Mod_Gather_Soln
  Use Mod_Gather_Drifter
  Use Mod_Communicate_Drifters
  Use Broadcast_Routines
#endif
#ifdef NCOUT
  use netcdf
  use CONVERTWGS84
#endif

  IMPLICIT NONE

  REAL(RKD),SAVE             :: DAYNEXT
  INTEGER,SAVE               :: NPDAY, NPDAY_MAX
  INTEGER,SAVE,ALLOCATABLE   :: NPLIST(:)
  INTEGER,SAVE,ALLOCATABLE   :: LADJ(:,:)

  REAL(RKD),POINTER,PRIVATE  :: ZCTR(:,:)
  REAL(RKD),SAVE,ALLOCATABLE :: LA_BEGTI(:)
  REAL(RKD),SAVE,ALLOCATABLE :: LA_ENDTI(:)
  REAL,SAVE,ALLOCATABLE      :: GRPWS(:)

  INTEGER,SAVE               :: NSTOPPED(10)
  INTEGER,SAVE,ALLOCATABLE   :: LA_GRP(:)
  INTEGER,SAVE,ALLOCATABLE   :: EW_Connectors(:)
  INTEGER,SAVE,ALLOCATABLE   :: NS_Connectors(:)
  INTEGER,SAVE,ALLOCATABLE   :: BEDFIX(:)           ! OPTION TO FIX A PARTICLE ON BED AFTER DEPOSITION ON BED
                             
  INTEGER,SAVE,ALLOCATABLE   :: LCTLU(:)
  INTEGER,SAVE,ALLOCATABLE   :: LCTLD(:)
  INTEGER,SAVE,ALLOCATABLE   :: LWRU(:)
  INTEGER,SAVE,ALLOCATABLE   :: LWRD(:)
                             
  LOGICAL,SAVE,ALLOCATABLE   :: ACTIVE(:)
                             
  REAL(RKD)                  :: DELTD, DIFFH, DIFFV
  REAL(8)                    :: EETIME              ! USING REAL*8 EXPLICITY FOR EE LINKAGE REQUIREMENTS
  INTEGER(IK4)               :: NGRP
  INTEGER(IK4)               :: ADJVERP             ! OPTION FOR ADJUSTING VERTICAL POSITION IN CASE OF FULL 3D
  INTEGER(IK4), SAVE         :: NTIMES

  INTEGER :: EE_UNIT = 95
  
  REAL(RKD), EXTERNAL :: DSTIME
#ifdef NCOUT
  integer(4) :: nc_lpt(0:6), lpt_time_dim, lpt_npd_dim, lpt_time_idx
#endif

  CONTAINS

  !---------------------------------------------------------------------------!
  !>@details SOLVE DIFFERENTIAL EQS. FOR (X,Y,Z):
  !! DX=U.DELTD+RAN.SQRT(2EH.DELTD)
  !! DY=V.DELTD+RAN.SQRT(2EH.DELTD)
  !! DZ=W.DELTD+RAN.SQRT(2EV.DELTD)
  !! U(L,K),V(L,K),W(L,K),K=1,KC,L=2:LA    CURRENT TIME
  !! U1(L,K),V1(L,K),W1(L,K),K=1,KC,L=2:LA PREVIOUS TIME
  !! N: TIME STEP
  !---------------------------------------------------------------------------!

SUBROUTINE DRIFTER_CALC

  Implicit none

  ! *** Local variables
  INTEGER(IK4) :: NPP, NP, NW, IT, L
  INTEGER(4)   :: VER, HSIZE, BSIZE, ISOIL        ! *** USING INTEGER*4 EXPLICITY FOR EE LINKAGE REQUIREMENTS

  REAL(RKD) :: KDX1, KDX2, KDX3, KDX4
  REAL(RKD) :: KDY1, KDY2, KDY3, KDY4
  REAL(RKD) :: KDZ1, KDZ2, KDZ3, KDZ4
  REAL(RKD) :: XLA1, YLA1, ZLA1
  REAL(RKD) :: U2NP, V2NP, W2NP
  REAL(RKD) :: DAHX1, DAHY1, DAVZ1, DAHX2, DAHY2, DAVZ2
  REAL(RKD) :: DXYMIN, KWEIGHT, TODAY
  REAL(RKD) :: TTDS

  REAL(RKD),SAVE :: TIMENEXT

  LOGICAL(IK4) :: BEDGEMOVE, LFORCE
  CHARACTER*80 :: TITLE,METHOD

  ! *** New variables for MPI
  Integer :: ierr, ii
  Integer :: in_west, in_east, in_north, in_south
  Integer :: send_tag
  Integer :: recv_tag
  Integer :: status_msg(MPI_Status_Size)
  Real(RKD) :: dxymin_global
  
  ! *** End new variables for MPI

  ! *** Decide if a dynamic or constant time step is used
  IF( ISDYNSTP == 0 )THEN ! *** Constant time step
    DELTD=DT
  ELSE ! *** Dynamic time step
    DELTD=DTDYN
  ENDIF
  
  ! *** FIRST CALL AT GLOBAL RELEASE TIME 
  IF( LA_FIRSTCALL == 1 )THEN
    DO L=2,LA
      ZCTR(L,0:KC) = ZZ(L,0:KC)
      ZCTR(L,KC+1) = Z(L,KC)            ! *** Z is the top layer of the interface (dimensionless)
    ENDDO

    ! *** Find the min cell dimension for the x and y directions.  Used for EE linkage scaling.
    DXYMIN = MIN(MINVAL(DXP_Global(2:LA_Global)), MINVAL(DYP_Global(2:LA_Global))) ! *** DXP cell dimension in x direction @ center
    IF( DXYMIN < 0.2 )THEN
      XYZSCL = 1000
    ELSE
      XYZSCL = 100
    ENDIF

    ! *** MAKE SURE THE FILE IS NEW - this is the binary file read by EE to visualize drifters
    if(process_id == master_id )then ! *** Only make the master process look at the file
#ifdef NCOUT
      call create_nc_lpt()
#endif      
      
      OPEN(ULGR,FILE=OUTDIR//'EE_DRIFTER.OUT',STATUS='UNKNOWN',FORM='BINARY')
      CLOSE(ULGR,STATUS='DELETE')
      OPEN(ULGR,FILE=OUTDIR//'EE_DRIFTER.OUT',ACTION='WRITE',FORM='BINARY')
      ISOIL = 0
      IF( ANY(ISOILSPI == 1) ) ISOIL = 1
      VER = 8400
      HSIZE = 6*4
      
      ! *** Write out header information for the drifter binary file
      WRITE(ULGR) VER, HSIZE
      WRITE(ULGR) INT(NPD,4), INT(KC,4), XYZSCL                 ! *** NPD --> is global
      WRITE(ULGR) ISOIL
        
      FLUSH(ULGR)
      CLOSE(ULGR,STATUS='KEEP')
    end if

    ! *** LA_FREQ - output frequency for the lagrangian calculation
    TIMENEXT = TIMEDAY + LA_FREQ + 0.000001_8
    LA_FIRSTCALL = 0
    NTIMES = 0

    ! *** If we are not reading from a restart file proceed
    IF( ISRESTI == 0 )THEN
      ! *** Call routine to initalize some stuff for the oil spill module
      IF( ANY(ISOILSPI == 1) ) CALL INIT_OIL                    ! *** GET: DVOL, DARE

      ! *** Loop over the number of drifters
      DO NP=1,NPD
        ! *** If the cell has not been initalized and the current time is within the lagrangian particle tracking time interval read in.
        IF( JSPD(NP) == 1 .AND. TIMEDAY >= LA_BEGTI(NP) .AND. TIMEDAY <= LA_ENDTI(NP) )THEN
          ! *** Only process drifters that are active in a given subdomain
          if( LLA_Process(NP) == process_id )then
            CALL CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)     ! *** Only NP is used for initialization
          end if
          JSPD(NP) = 0                                          ! *** Particle is now initialized
        ENDIF
      ENDDO
    ELSE ! *** Read from the restart file
      !< @todo Update to read the restart file
      ! *** Warning not working for MPI
      write(*,*) '*** WARNING *** Restart functionality for Lagrangian Particle Tracking module is not available'
      CALL DRIFTER_RST

      ! *** INITIALIZE BOTTOM ELEVATION AND DEPTH
      DO NP=1,NPD
        IF( LLA(NP) > 1 )THEN
          CALL DRF_DEPTH(LLA(NP),NP,BELVLA(NP),HPLA(NP))
        ENDIF
      ENDDO

    ENDIF

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
  ENDIF ! *** End first call to this routine
  
  ! *** BUILD DAY LIST OF DRIFTERS
  IF( TIMEDAY > DAYNEXT )THEN
      
    ! *** REPORT ANY DRIFTERS THAT MAY HAVE STOPPED FOR WHATEVER REASON
    IF( SUM(NSTOPPED) > 0 )THEN
      IF( NSTOPPED(1) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS EXPIRED DUE TO EVAPOR. & BIODEG.:   ',NSTOPPED(1), ' ON PROCESS: ', process_id
      IF( NSTOPPED(2) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS STOPPED DUE TO DEPOSITION:          ',NSTOPPED(2), ' ON PROCESS: ', process_id
      IF( NSTOPPED(3) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO HYDRAULIC STRUCTURE: ',NSTOPPED(3), ' ON PROCESS: ', process_id
      IF( NSTOPPED(4) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO OPEN BOUNDARY:       ',NSTOPPED(4), ' ON PROCESS: ', process_id
      IF( NSTOPPED(5) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO WITHDAWAL:           ',NSTOPPED(5), ' ON PROCESS: ', process_id
      IF( NSTOPPED(6) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS OUTSIDE DUE TO WITHDRAWAL/RETURN:   ',NSTOPPED(6), ' ON PROCESS: ', process_id
      IF( NSTOPPED(7) > 0 ) PRINT '(A,I8,A,I8)','# DRIFTERS STOPPED DUE TO CELL DRYING:         ',NSTOPPED(7), ' ON PROCESS: ', process_id
      NSTOPPED = 0
    ENDIF
    Call MPI_barrier(MPI_Comm_World, ierr)

    TODAY = DBLE(INT(TIMEDAY))
    DAYNEXT = TODAY + 1.
    
    ! *** Update the active drifters for the next day
    NPDAY = 0
    DO NP=1,NPD
      ! *** Exclude drifters not in current domain
      if( LLA_Process(NP) /= process_id )then
        NPP = 0 ! FOR DEBUGGING
        CYCLE
      endif
      
      ! *** Exclude stopped drifters
      IF( JSPD(NP) == 0 .AND. LLA_Global(NP) < 2 )THEN
        NPP = 0 ! FOR DEBUGGING
        CYCLE
      ENDIF

      ! *** Exclude drifters whose timing is not in the current time window
      ! *** LA_BEGTI(NP) < LA_BEGTI0:  Released before start of LPT computations
      ! *** LA_BEGTI(NP) > DAYNEXT:    Released after next day
      ! *** LA_ENDTI(NP) < TODAY:      Stopped tracking before current day
      IF( LA_BEGTI(NP) < LA_BEGTI0 .OR. LA_BEGTI(NP) > DAYNEXT .OR. LA_ENDTI(NP) < TODAY )THEN
        NPP = 0 ! FOR DEBUGGING
        CYCLE
      ENDIF
      
      NPDAY = NPDAY + 1
      NPLIST(NPDAY) = NP
    ENDDO
    PRINT '(A,I8,A,I8)','# DRIFTERS ACTIVE FOR THE NEXT PERIOD:         ',NPDAY, ' ON PROCESS: ', process_id

    Call MPI_Allreduce(NPDAY, NPDAY_MAX, 1, MPI_Integer, MPI_Max, comm_2d, ierr)
    
  ENDIF

  ! *** NEXT CALL --------------------------------------------------------------
  DAHX1=0
  DAHY1=0
  DAVZ1=0
  DAHX2=0
  DAHY2=0
  DAVZ2=0
  MOC  =0
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
  DO NPP=1,NPDAY ! *** Loop over the number of active drifters for this time step
    !$ IT = OMP_GET_THREAD_NUM() + 1
  
    ! *** GET CURRENT ACTIVE DRIFTER
    NP = NPLIST(NPP)

    ! *** If we are outside of the time window we should be tracking this drifter then move onto next drifter
    IF( TIMEDAY < LA_BEGTI(NP) .OR. TIMEDAY > LA_ENDTI(NP) )THEN
      IF( JSPD(NP) == 0 )THEN
        ! *** Drifter has been initialized.  Update LLA arrays
        IF( LLA(NP) > 1 .AND. TIMEDAY > LA_ENDTI(NP) )THEN
          LLA_Global(NP) = 1
        ELSE
          LLA_Global(NP) = 0
        ENDIF
        LLA(NP) = 0
      ENDIF
      CYCLE
    ENDIF
    
    ! *** Drifter is not active yet and the time is above the beginning of the drifter activation time
    IF( JSPD(NP) == 1 .AND. TIMEDAY >= LA_BEGTI(NP) )THEN
      ! *** Initialize Drifter
      CALL CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)
      JSPD(NP) = 0
    ENDIF

    ! *** Call some subroutines specific to the oil spill module
    IF( ISOILSPI(LA_GRP(NP)) == 1 .AND. LLA(NP) >= 2 )THEN
      CALL OIL_PROC(NP)            ! *** CALCULATE DVOL(NP)
      IF( DVOL(NP) <= 1D-9 )THEN
        LLA(NP) = 0
        LLA_Global(NP) = 1
        IF( NPD < 10000 .OR. DEBUG )THEN
          PRINT '(A41,I8)','DRIFTER EXPIRED DUE TO EVAPOR. & BIODEG.:',NP
        ELSE
          NSTOPPED(1) = NSTOPPED(1) + 1
        ENDIF
        CYCLE
      ENDIF
    ENDIF

    ! *** If this particle is fixed in a bed and within the active domain
    IF( BEDFIX(LA_GRP(NP)) == 1 .AND. LLA(NP) >= 2  )THEN
      ! *** CHECK FOR DEPOSITION
      IF( (ZLA(NP) - BELVLA(NP)) <= 1D-6  )THEN
        LLA(NP) = 0
        LLA_Global(NP) = 1

        IF( NPD < 10000 .OR. DEBUG )THEN
          PRINT '(A38,I8)','DRIFTER HAS STOPPED DUE TO DEPOSITION:',NP
        ELSE
          NSTOPPED(2) = NSTOPPED(2) + 1
        ENDIF

        CYCLE
      ENDIF
      
      ! *** CHECK FOR DRY CELL
      IF( ISDRY > 0 )THEN
        IF( .NOT. LMASKDRY(LLA(NP))  )THEN
          LLA(NP) = 0
          LLA_Global(NP) = 1
          IF( NPD < 10000 .OR. DEBUG )THEN
            PRINT '(A38,I8)','DRIFTER HAS STOPPED DUE TO CELL DRYING:',NP
          ELSE
            NSTOPPED(7) = NSTOPPED(7) + 1
          ENDIF

          CYCLE
        ENDIF
      ENDIF
    ENDIF
    
    ! *** If the cell index is less than 2 the drifter has stopped (0), left the domain (0) or exited the sub-domain (1)
    IF( LLA(NP) < 2 ) CYCLE

    ! *** If the cell has dried,
    ! *** checks the water surface elevation against the depth at which we have decided a cell goes dry --> somewhat redundant?
    IF( ISDRY > 0 .AND. HP(LLA(NP)) < HDRY ) CYCLE
    
    ! *** Set the coordinates for the current cell
    XLA1 = XLA(NP)
    YLA1 = YLA(NP)
    ZLA1 = ZLA(NP)

    NW = LA_GRP(NP) ! *** Set the group for this drifter, this was read in by DRIFTER_INP
    BEDGEMOVE = .FALSE.
    
    ! *** APPLY DEPTH CHANGE ON VERTICAL POSITION
    IF( LA_ZCAL == 1 .AND. ADJVERP == 1 )THEN
      ZLA1 = ZLA1 + (HP(LLA(NP))-H1P(LLA(NP)))*(ZLA1-BELVLA(NP))/HPLA(NP)
    ENDIF

    ! *** If horitzontal diffusivity is turned off in EE
    IF( LA_DIFOP /= 0 .OR. ISHDMF == 0 )THEN
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
    ELSE
      KWEIGHT = 6.0
    ENDIF
    ! ***************************************************************************************
    ! *** RUNGE-KUTTA
    ! *** PASS 1
    CALL DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    IF( LA_DIFOP == 0 .AND. ISHDMF > 0 ) CALL DIFGRAD(LLA(NP),KLA(NP),NP,DAHX1,DAHY1,DAVZ1)
    KDX1 = DELTD*(U2NP + DAHX1)
    KDY1 = DELTD*(V2NP + DAHY1)
    XLA(NP)  = XLA1 + 0.5*KDX1
    YLA(NP)  = YLA1 + 0.5*KDY1
    IF( LA_ZCAL == 1 )THEN
      KDZ1 = DELTD*(W2NP - GRPWS(NW) + DAVZ1)
      ZLA(NP) = ZLA1 + 0.5*KDZ1
    ELSE
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    ENDIF
    CALL CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)

    IF( LLA(NP)<2 .OR. BEDGEMOVE) CYCLE

    ! *** PASS 2
    CALL DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    IF( LA_DIFOP == 0 .AND. ISHDMF>0 ) CALL DIFGRAD(LLA(NP),KLA(NP),NP,DAHX2,DAHY2,DAVZ2)
    KDX2 = DELTD*(U2NP + DAHX2)
    KDY2 = DELTD*(V2NP + DAHY2)
    XLA(NP) = XLA1 + 0.5*KDX2
    YLA(NP) = YLA1 + 0.5*KDY2
    IF( LA_ZCAL == 1 )THEN
      KDZ2 = DELTD*(W2NP - GRPWS(NW) + DAVZ2)
      ZLA(NP) = ZLA1 + 0.5*KDZ2
    ELSE
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    ENDIF
    CALL CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)

    IF( LLA(NP)<2 .OR. BEDGEMOVE) CYCLE

    ! *** PASS 3
    CALL DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    IF( LA_DIFOP == 0 .AND. ISHDMF>0 ) CALL DIFGRAD(LLA(NP), KLA(NP), NP, DAHX2, DAHY2, DAVZ2)
    KDX3 = DELTD*(U2NP + DAHX2)
    KDY3 = DELTD*(V2NP + DAHY2)
    XLA(NP) = XLA1 + KDX3
    YLA(NP) = YLA1 + KDY3
    IF( LA_ZCAL == 1 )THEN
      KDZ3 = DELTD*(W2NP - GRPWS(NW) + DAVZ2)
      ZLA(NP) = ZLA1 + KDZ3
    ELSE
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    ENDIF
    CALL CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)

    IF( LLA(NP)<2 .OR. BEDGEMOVE ) CYCLE

    ! *** PASS 4
4   CONTINUE
    CALL DRF_VELOCITY(LLA(NP),KLA(NP),NP,U2NP,V2NP,W2NP)

    IF( LA_DIFOP == 0 .AND. ISHDMF>0 ) CALL DIFGRAD(LLA(NP),KLA(NP),NP,DAHX2,DAHY2,DAVZ2)
    KDX4 = DELTD*(U2NP + DAHX2)
    KDY4 = DELTD*(V2NP + DAHY2)
    XLA(NP) = XLA1 + (KDX1 + 2.0*KDX2 + 2.0*KDX3 + KDX4)/KWEIGHT
    YLA(NP) = YLA1 + (KDY1 + 2.0*KDY2 + 2.0*KDY3 + KDY4)/KWEIGHT
    IF( LA_ZCAL == 1 )THEN
      KDZ4 = DELTD*(W2NP-GRPWS(NW)+DAVZ2)
      ZLA(NP) = ZLA1 + (KDZ1 + 2.0*KDZ2 + 2.0*KDZ3 + KDZ4)/KWEIGHT
    ELSE
      ZLA(NP) = HPLA(NP) + BELVLA(NP) - DLA(NP)
    ENDIF

    ! *** WIND-DRIFT COMPONENT
    IF( IOSWD > 0 )THEN
      IF( DLA(NP) < 0.5 ) CALL WIND_DRIFT(LLA(NP), KLA(NP), NP)
    ENDIF

    ! *** RANDOM-WALK COMPONENT
    IF( LA_PRAN > 0 ) CALL RANDCAL(LLA(NP), KLA(NP), NP)
    
    CALL CONTAINER(BEDGEMOVE, XLA1, YLA1, ZLA1, NP)
    
    IF(LLA(NP) >=2 .AND. LLA(NP) <= LA .AND. ISOILSPI(LA_GRP(NP)) == 1  )THEN
      MOC(LLA(NP)) = MOC(LLA(NP)) + DVOL(NP)*DDEN(LA_GRP(NP))                ! TOTAL MASS OF OIL (KG) AT A CELL
    ENDIF

    ! *** Check if the drifter has moved into the ghost cell boundary
    Call Check_Drifter_In_Ghost(NP)

  ENDDO
  !$OMP END PARALLEL DO
  
  ! *** Communicate ghost info if a drifter entered the ghost cell region
  ! *** Get the max number of drifters to communicate amongst all directions
  local_tot_drifters_send = num_drifters_send_west + num_drifters_send_east + num_drifters_send_north + num_drifters_send_south
  
  ! *** Get global max number of drifters to communicate amongst all processes
  Call MPI_barrier(MPI_Comm_World, ierr)
  
  TTDS = DSTIME(0)
  global_max_drifters_to_comm = 0
  Call MPI_Allreduce(local_tot_drifters_send, global_max_drifters_to_comm, 1, MPI_Integer, MPI_Max, comm_2d, ierr)
  DSITIMING(11) = DSITIMING(11) + (DSTIME(0)-TTDS)

  ! *** Only need to communicate if there are drifters in any of the ghost cells of any subdomain
  if( global_max_drifters_to_comm > 0 )then
    ! *** If a drifter is leaving a subdomain let the subdomain know it will be receiving data
    Call Notify_Receiving_Domains
    
    ! *** Communicate all drifter info
    Call Communicate_Drifters(LLA_Global, 1)

    !Call Communicate_Drifters(JSPD)  ! by definition, if the drifter is communicated JSPD = 0
     
    Call Communicate_Drifters(XLA)
     
    Call Communicate_Drifters(YLA)

    Call Communicate_Drifters(ZLA)

    Call Communicate_Drifters(DLA)

    Call Communicate_Drifters(HPLA)
     
    Call Communicate_Drifters(KLA)
     
    Call Communicate_Drifters(BELVLA)

    ! *** Zero send list for next iteration
    drifter_ids_send_west(1:global_max_drifters_to_comm)  = 0
    drifter_ids_send_east(1:global_max_drifters_to_comm)  = 0
    drifter_ids_send_north(1:global_max_drifters_to_comm) = 0
    drifter_ids_send_south(1:global_max_drifters_to_comm) = 0

    drifter_ids_recv_west(1:global_max_drifters_to_comm)  = 0
    drifter_ids_recv_east(1:global_max_drifters_to_comm)  = 0
    drifter_ids_recv_north(1:global_max_drifters_to_comm) = 0
    drifter_ids_recv_south(1:global_max_drifters_to_comm) = 0
  end if ! *** End sequence that communicates drifters 
  
  ! *** WRITE THE CURRENT TRACK POSITION
  IF( TIMEDAY >= TIMENEXT )THEN   ! .OR. TIMEDAY+1E-5 >= TIMEEND )THEN
    CALL DRIFTER_OUT(.FALSE.)
    TIMENEXT = TIMENEXT + LA_FREQ
  ENDIF
  
END SUBROUTINE DRIFTER_CALC
!---------------------------------------------------------------------------!
!> @details Lets other process know that it will be receiving drifter data
!  and how many drifters 
!
!> @author Zander Mausolff 
!
!---------------------------------------------------------------------------! 
Subroutine Notify_Receiving_Domains

	Implicit None

  ! *** Local variables
  Integer :: status_msg(MPI_Status_Size)
  Integer :: ierr, ii, iii, iFound, NP
  Integer :: total_send_drifters
     
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
  Call MPI_Sendrecv(drifter_ids_send_west, global_max_drifters_to_comm, MPI_Integer, nbr_west, 0, &
                    drifter_ids_recv_east, global_max_drifters_to_comm, MPI_Integer, nbr_east, 0, &
                    comm_2d, status_msg, ierr)
     
  ! *** East/West send/recv
  Call MPI_Sendrecv(drifter_ids_send_east, global_max_drifters_to_comm, MPI_Integer, nbr_east, 0, &
                    drifter_ids_recv_west, global_max_drifters_to_comm, MPI_Integer, nbr_west, 0, &
                    comm_2d, status_msg, ierr)

  ! ***
  Call MPI_Sendrecv(drifter_ids_send_north, global_max_drifters_to_comm, MPI_Integer, nbr_north, 0, &
                    drifter_ids_recv_south, global_max_drifters_to_comm, MPI_Integer, nbr_south, 0, &
                    comm_2d, status_msg, ierr)     

  ! ***
  Call MPI_Sendrecv(drifter_ids_send_south, global_max_drifters_to_comm, MPI_Integer, nbr_south, 0, &
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
      Call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_east(ii))
    endif
  enddo
    
  do ii = 1,recv_west
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_west(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      Call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_west(ii))
    endif
  enddo
    
  do ii = 1,recv_north
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_north(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      Call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_north(ii))
    endif
  enddo

  do ii = 1,recv_south
    iii = Search_List(NPDAY, NPLIST, drifter_ids_recv_south(ii))
    if( iii < 1 )then
      ! *** Drifter not found.  Add to list
      Call Insert_Value(NPDAY, NPLIST, drifter_ids_recv_south(ii))
    endif
  enddo
  
End subroutine Notify_Receiving_Domains


!---------------------------------------------------------------------------!
!> @details Reads the drifter restart file
!> @todo Enable remapping for domain decomposition
!---------------------------------------------------------------------------!  
SUBROUTINE DRIFTER_RST
  ! *** READ LPT RESTART FILE
  INTEGER :: NP
  LOGICAL FEXIST

  WRITE(*,'(A)')'READING RESTART FILE: DRIFTER.RST'

  ! CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE THE ASCII FILE INSTEAD.
  RESFILE = 'drifter.rst'
  INQUIRE(FILE=RESFILE, EXIST=FEXIST)

  IF( FEXIST )THEN
    OPEN(99,FILE=RESFILE)

    READ(99,'(50I2)') JSPD(1:NPD)
    READ(99,'(30I3)') KLA(1:NPD)
    READ(99,'(2(I8,2F15.5,F10.5))') (LLA_Global(NP), XLA(NP), YLA(NP), ZLA(NP),NP=1,NPD)

    CLOSE(99)
  ENDIF

END SUBROUTINE DRIFTER_RST
!---------------------------------------------------------------------------!
!> @details Writes to the EE_drifter.out binary file for EE linkage
!
!> @author Zander Mausolff (Updated for MPI)
!> @date 2/17/2020
!---------------------------------------------------------------------------!  
SUBROUTINE DRIFTER_OUT(FORCEEXIT)
  REAL(RKD)           :: TTDS, TTDS1                     ! MODEL TIMING TEMPORARY VARIABLES

  LOGICAL, INTENT(IN) :: FORCEEXIT

  ! *** Local variables
  LOGICAL :: LFORCE
  INTEGER(IK4) :: NP, NACT
  INTEGER(4) :: IX,IY,IZ,IV     ! *** USING INTEGER*4 EXPLICITY FOR EE LINKAGE REQUIREMENTS

  ! *** All processes must send their data
  TTDS = DSTIME(0)
  
  Call Gather_Drifter_Arrays
  
  TTDS1 = DSTIME(0) - TTDS
  TMPIEE = TMPIEE + TTDS1
  TLRPD = TLRPD - TTDS1

  ! *** Only proceed with master process
  if( process_id == master_id )then

    IF( FORCEEXIT )THEN
      ! *** USER INITIATED EXIT.  WRITE FINAL POSITION OF EVERY DRIFTER
      PRINT *,'COMPLETE DRIFTER SNAPSHOT (USER): ', TIMEDAY
      LFORCE = .TRUE.
      NACT = NPDAY
    ELSE
      IF( TIMEDAY+LA_FREQ >= TIMEEND )THEN
        ! *** WRITE FINAL POSITION OF EVERY DRIFTER AT END OF RUN
        LFORCE = .TRUE.
        PRINT *,'COMPLETE DRIFTER SNAPSHOT (EoR): ', TIMEDAY
        NACT = NPD
      ELSE
        LFORCE = .FALSE.
        NACT = 0
        DO NP=1,NPD
          IF( LLA_Global(NP) >= 1 )THEN
            NACT = NACT + 1
          ENDIF
        ENDDO
      ENDIF
    ENDIF

#ifdef NCOUT
    call write_nc_lpt()
#endif
    
    ! *** WRITE TO THE FILE
    OPEN(ULGR,FILE=OUTDIR//'EE_DRIFTER.OUT',STATUS='UNKNOWN',FORM='BINARY',POSITION='APPEND')
    EETIME = TIMEDAY

    WRITE(ULGR) EETIME

    NTIMES = NTIMES + 1

    WRITE(ULGR) INT(NACT,4)

    IF( ANY(ISOILSPI == 1) )THEN
      DO NP=1,NPD
        IF( LLA_Global(NP) >= 1 .OR. LFORCE )THEN
          IX = NINT(XYZSCL*XLA(NP))                                     ! *** SCALED
          IY = NINT(XYZSCL*YLA(NP))                                     ! *** SCALED
          IZ = NINT(XYZSCL*ZLA(NP))                                     ! *** SCALED
          IV = NINT(1D6*DVOL(NP))                                       ! *** IN CUBIC CM
          WRITE(ULGR) INT(NP,4), INT(LLA_Global(NP),4), IX, IY, IZ, IV
        ENDIF
      ENDDO
    ELSE
      DO NP=1,NPD
        IF( LLA_Global(NP) >= 1 .OR. LFORCE )THEN
          IX = NINT(XYZSCL*XLA(NP))                                     ! *** SCALED
          IY = NINT(XYZSCL*YLA(NP))                                     ! *** SCALED
          IZ = NINT(XYZSCL*ZLA(NP))                                     ! *** SCALED
          WRITE(ULGR) INT(NP,4), INT(LLA_Global(NP),4), IX, IY, IZ      
        ENDIF
      ENDDO
    ENDIF
    FLUSH(ULGR)

    CLOSE(ULGR,STATUS='KEEP')
  endif
  
  ! *** Write out the last position then set to inactive for future writes
  DO NP=1,NPD
    IF( LLA_Global(NP) == 1 )THEN
      LLA_Global(NP) = 0                  
    ENDIF
  ENDDO
      
END SUBROUTINE DRIFTER_OUT
!---------------------------------------------------------------------------!
!
!> @details READING INPUT DATA OF INITIAL LOCATIONS OF DRIFTERS and mapping them out
!! to subdomains if necessary
!!OUTPUT: NPD,XLA,YLA,ZLA,NP=1:NPD, LA_BEGTI, LA_ENDTI, LA_FREQ,LANDT\
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
!          0 Use AH(L,K) and AV(L,K) from EFDC computations (Horizontal diffusion should be turned on in EFDC)
!          1 Use HORDIF  VERDIF from this file
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
  Implicit none

  ! *** Local variables

  INTEGER :: I, J, NP, N, K, NPN, IN, JN, IS, JS, L, NC
  REAL(RKD) :: RANVAL
  !REAL(8),EXTERNAL :: DRAND   !IT NEEDS THIS STATEMENT IN CASE OF IMPLICIT NONE @todo Zander 2/14/20 - seems like it is not?
  CHARACTER(200)  :: STR
  REAL(RKD)  :: Distance_to_Cells(LA_Global)

  INTEGER :: Min_L_loc_index(1)
  Integer :: I_global, J_global, L_local, LLO, NI, i_local, j_local
  Integer, Allocatable, Dimension(:) :: new_seed
  Integer :: seed_size !< size of the array that represents the random number seed
  Integer :: LG

  ! *** Do reading only on the master process
  if(process_id == master_id )then
    WRITE(*,'(A)')'READING DRIFTER.INP'
    OPEN(ULOC,FILE='drifter.inp',ACTION='READ')
    CALL SKIPCOM(ULOC,'*')
    READ(ULOC,*) LA_ZCAL, LA_PRAN, LA_DIFOP, LA_HORDIF, LA_VERDIF, DEPOP, ADJVERP, SLIPFACTOR, IOSWD, OSWDA, OSWDB
    CALL SKIPCOM(ULOC,'*')
    READ(ULOC,*) LA_BEGTI0, LA_ENDTI0, LA_FREQ
    CALL SKIPCOM(ULOC,'*')
    READ(ULOC,*) NPD,NGRP
    WRITE(*,'(A41,I8)')'DRIFTER: NUMBER OF DRIFTERS INITIALZED: ',NPD
  end if

  ! *** Broadcast read in variables to all processes
  Call Broadcast_Scalar(LA_ZCAL,    master_id)
  Call Broadcast_Scalar(LA_PRAN,    master_id)
  Call Broadcast_Scalar(LA_DIFOP,   master_id)
  Call Broadcast_Scalar(LA_HORDIF,  master_id)
  Call Broadcast_Scalar(LA_VERDIF,  master_id)
  Call Broadcast_Scalar(DEPOP,      master_id)
  Call Broadcast_Scalar(ADJVERP,    master_id)
  Call Broadcast_Scalar(SLIPFACTOR, master_id)
  Call Broadcast_Scalar(IOSWD,      master_id)
  Call Broadcast_Scalar(OSWDA,      master_id)
  Call Broadcast_Scalar(OSWDB,      master_id)

  Call Broadcast_Scalar(LA_BEGTI0,  master_id)
  Call Broadcast_Scalar(LA_ENDTI0,  master_id)
  Call Broadcast_Scalar(LA_FREQ,    master_id)
  Call Broadcast_Scalar(NPD,        master_id)
  Call Broadcast_Scalar(NGRP,       master_id)

  ! *** Drifter output frequency
  LA_FREQ = LA_FREQ/1440.
  IF( LA_FREQ <= 0. ) LA_FREQ = DBLE(TIDALP)/DBLE(NPPPH)/DBLE(86400.)

  ! NG      : Order number of group
  ! ISOILSPI: = 1 for oil spill; = 0 for normal drifter as usual
  ! GRPWS   : Settling velocity of drifter [m/s]
  ! GVOL    : Total volume of oil spill for a group [m^3]
  ! DDEN    : Density of an oil dritfer    [kg/m^3]
  ! DRAT    : First order biodegradation rate of an oil drifter [1/day,>=0]
  ! DTEM    : Biodegradation Rate Reference Temperature [deg C]
  ! DVAP    : Vapor pressure of an oil drifter     [Pa]
  ! DVMO    : Molar volume of an oil drifter at STP [m^3/mol]

  ! *** ALLOCATE THE DRIFTER ARRAYS
  ALLOCATE(XLA(NPD), YLA(NPD), ZLA(NPD), DLA(NPD))
  ALLOCATE(LLA(NPD), KLA(NPD), HPLA(NPD), BELVLA(NPD))
  ALLOCATE(LA_BEGTI(NPD),LA_ENDTI(NPD), LA_GRP(NPD), JSPD(NPD))
  ALLOCATE(ACTIVE(NPD), NPLIST(NPD))
  ALLOCATE(LADJ(9,LCM))

  ALLOCATE(ZCTR(LCM, 0:KC+1))
  ALLOCATE(NS_Connectors(LCM), EW_Connectors(LCM))
  ALLOCATE(BEDFIX(NGRP))
  ALLOCATE(DARE(NGRP), DTEM(NGRP), GRPWS(NGRP), ISOILSPI(NGRP), DVOL0(NGRP))
  ALLOCATE(DVAP(NGRP), DVMO(NGRP), DDEN(NGRP), DRAT(NGRP), GVOL(NGRP))

  ALLOCATE(MOC(LA))

  ! *** New variables for MPI - global arrays for the initial read and for writing out
  Allocate(LLA_Global(NPD))

  !< @todo resize these to reduce memory usage
  Allocate(drifter_ids_send_west(NPD))
  Allocate(drifter_ids_send_east(NPD))
  Allocate(drifter_ids_send_north(NPD))
  Allocate(drifter_ids_send_south(NPD))

  Allocate(drifter_ids_recv_east (NPD))
  Allocate(drifter_ids_recv_west (NPD))
  Allocate(drifter_ids_recv_north(NPD))
  Allocate(drifter_ids_recv_south(NPD))

  Allocate(LLA_Process(NPD))

  ! *** INIITALIZE THE ARRAYS (NPD)
  XLA=0.0
  YLA=0.0
  ZLA=0.0
  DLA=0.0
  HPLA=0.0
  BELVLA=0.0
  LA_BEGTI=0.0
  LA_ENDTI=0.0
  LA_GRP=0
  LLA = 0
  KLA = 0
  JSPD= 0
  MOC = 0
  ACTIVE = .TRUE.
  NPLIST = 0


  ! *** INIITALIZE THE ARRAYS (NGRP)
  GRPWS=0.0
  DARE=0.0
  DTEM=0.0
  DVAP=0.0
  DVMO=0.0
  DDEN=0.0
  DRAT=0.0
  GVOL=0.0
  DVOL0=0
  ISOILSPI=0

  ZCTR=0.0

  ! *** end initialize new MPI arrays

  ! *** Do reading only on master
  if(process_id == master_id )then
    CALL SKIPCOM(ULOC,'*')
    DO N=1,NGRP
      READ(ULOC,'(A)') STR
      READ(STR,*) K,ISOILSPI(K),GRPWS(K),BEDFIX(K)
      IF( ISOILSPI(K) > 0 )THEN
        READ(STR,*) K,ISOILSPI(K),GRPWS(K),BEDFIX(K),GVOL(K),DDEN(K),DRAT(K),DTEM(K),DVAP(K),DVMO(K)
        DRAT(K)=DRAT(K)/86400.
      ENDIF
    ENDDO
  end if

  ! *** Broadcast arrays just read in by the master process
  Call Broadcast_Array(ISOILSPI, master_id)
  Call Broadcast_Array(GRPWS   , master_id)
  Call Broadcast_Array(BEDFIX  , master_id)
  Call Broadcast_Array(GVOL    , master_id)
  Call Broadcast_Array(DDEN    , master_id)
  Call Broadcast_Array(DRAT    , master_id)
  Call Broadcast_Array(DTEM    , master_id)
  Call Broadcast_Array(DVAP    , master_id)
  Call Broadcast_Array(DVMO    , master_id)

  ! *** INITIALIZE FOR THE FIRST CALL
  LA_FIRSTCALL = 1
  JSPD = 1           ! *** JSPD = 1 IF NOT INITIALIZED,  JSPD = 0 AFTER CELL IS INTIALIZED
  
  NS_Connectors = 0  ! *** Default - No connection
  DO NPN=1,NPNSBP
    IN = INPNS(NPN)
    JN = JNPNS(NPN)
    IS = ISPNS(NPN)
    JS = JSPNS(NPN)
    NS_Connectors(LIJ(IN,JN)) = 1
    NS_Connectors(LIJ(IS,JS)) = 2
  ENDDO

  EW_Connectors = 0  ! *** Default - No connection
  DO NPN=1,NPEWBP
    IN = IWPEW(NPN)
    JN = JWPEW(NPN)
    IS = IEPEW(NPN)
    JS = JEPEW(NPN)
    EW_Connectors(LIJ(IN,JN)) = 1
    EW_Connectors(LIJ(IS,JS)) = 2
  ENDDO

  ! *** OBTAIN ACTIVE CELLS SURROUNDING EACH CELL
  LADJ = 1
  DO L=2,LA
    ! *** ORDER OF CELLS
    ! ***   1  2  3
    ! ***   4  5  6
    ! ***   7  8  9
    LADJ(5,L) = L
    IF( SUBO(L)      + SVBO(LNWC(L)) > 1.5 .OR. SVBO(LNC(L)) + SUBO(LNC(L))  > 1.5 ) LADJ(1,L) = LNWC(L)
    IF( SVBO(LNC(L)) + SUBO(LNEC(L)) > 1.5 .OR. SUBO(LEC(L)) + SVBO(LNEC(L)) > 1.5 ) LADJ(3,L) = LNEC(L)
    IF( SUBO(L)      + SVBO(LWC(L))  > 1.5 .OR. SVBO(L)      + SUBO(LSC(L))  > 1.5 ) LADJ(7,L) = LSWC(L)
    IF( SVBO(L)      + SUBO(LSEC(L)) > 1.5 .OR. SUBO(LEC(L)) + SVBO(LEC(L))  > 1.5 ) LADJ(9,L) = LSEC(L)

    IF( SVBO(LNC(L)) > 0.5 ) LADJ(2,L) = LNC(L)
    IF( SUBO(L) > 0.5 )      LADJ(4,L) = LWC(L)
    IF( SUBO(LEC(L)) > 0.5 ) LADJ(6,L) = LEC(L)
    IF( SVBO(L) > 0.5 )      LADJ(8,L) = LSC(L)
  ENDDO

  ! *** Do reading only on master
  if(process_id == master_id)then
    CALL SKIPCOM(ULOC,'*')
    IF( DEPOP == 1 )THEN
      DO NP=1,NPD
        ! *** Read Depths
        READ(ULOC,*,ERR=999) XLA(NP), YLA(NP), DLA(NP), LA_BEGTI(NP), LA_ENDTI(NP), LA_GRP(NP)
        DLA(NP) = MAX(DLA(NP),0.0_RKD)
        IF( LA_GRP(NP) < 1 ) LA_GRP(NP) = 1
      ENDDO
    ELSE
      DO NP=1,NPD
        ! *** Read Elevations
        READ(ULOC,*,ERR=999) XLA(NP), YLA(NP), ZLA(NP), LA_BEGTI(NP), LA_ENDTI(NP), LA_GRP(NP)
        IF( LA_GRP(NP) < 1 ) LA_GRP(NP) = 1
      ENDDO
    ENDIF
    CLOSE(ULOC)
  end if

  ! *** Broadcast the arrays to all processes
  Call Broadcast_Array(XLA, master_id)
  Call Broadcast_Array(YLA, master_id)
  if( DEPOP == 1  )then
    Call Broadcast_Array(DLA, master_id)
  else
    Call Broadcast_Array(ZLA, master_id)
  end if
  Call Broadcast_Array(LA_BEGTI, master_id)
  Call Broadcast_Array(LA_ENDTI, master_id)
  Call Broadcast_Array(LA_GRP,   master_id)

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
    
    IF( LA_BEGTI(NP) >= LA_BEGTI0 .AND. LA_BEGTI(NP) <= LA_ENDTI0 )THEN
      ! *** Get Local L
      DO J = 3,JC-2
        DO I = 3,IC-2
          L = LIJ(I,J)
          if( L > 0 )then
            ! *** Test if Drifter is inside the cell
            IF( INSIDECELL(L, XLA(NP), YLA(NP)) )THEN
              LLA(NP) = L
              ! *** LLA_Global(NP) is initialized only when activated
              LLA_Process(NP) = Map2Global(L).PR
              Exit
            ENDIF
          ENDIF
        ENDDO
        IF( LLA(NP) > 1 ) EXIT
      enddo
    endif
  end do

  ! *** Write out information regaring the local mapping of drifters
  if(MPI_DEBUG_FLAG == .TRUE. )then
    write(mpi_log_unit,'(a,i8)') 'Total number of local drifters: ', NPD

    Call writebreak(mpi_log_unit)
    write(mpi_log_unit,'(a)') 'Writing out local drifter values'
    write(mpi_log_unit,'(a)') '               NP,      LLA, LLA_GL,  Proc ID, XLA(NP),      YLA(NP),          DLA(NP),    LA_BEGTI(NP), LA_ENDTI(NP),  LA_GRP(NP)'
    ! *** Write out - TESTING/DEBUGGING only
    do np = 1, NPD

      ! *** If the cell is inside the subdomain
      if( DEPOP == 1  )then
        write(mpi_log_unit,'(a,4i8,5f15.4,i8)') 'Drifter: ', NP, LLA(NP), LLA_Global(NP), LLA_Process(NP), XLA(NP), YLA(NP), ZLA(NP), LA_BEGTI(NP),LA_ENDTI(NP),LA_GRP(NP)
      else
        write(mpi_log_unit,'(a,4i8,5f15.4,i8)') 'Drifter: ', NP, LLA(NP), LLA_Global(NP), LLA_Process(NP), XLA(NP), YLA(NP), DLA(NP), LA_BEGTI(NP),LA_ENDTI(NP),LA_GRP(NP)
      end if

    end do
  end if

  !< @todo Improve random number generation
  ! *** This is where the random number is first set.
  ! *** This should be adjusted to instead call the random number seed on all processes
  ! *** and then call the random number.  That way we can make repeatable events with pseudo random numbers

  ! *** Setup the random seed
  seed_size = 3
  Allocate(new_seed(seed_size))
  new_seed(:) = process_id + 1
  Call RANDOM_SEED(put = new_seed(1:seed_size))

  IF( LA_PRAN > 0  )then
    Call RANDOM_NUMBER(RANVAL)  ! *** Replaced 2/14/2020 --> RANVAL = DRAND(1)
  end if

  ! *** SETUP L LOOKUPS FOR HYDRAULIC STRUCTURES
  ALLOCATE(LCTLU(MAX(NQCTL,1)),LCTLD(MAX(NQCTL,1)))
  LCTLU = 0
  LCTLD = 0
  DO NC=1,NQCTL
    LCTLU(NC) = LIJ(HYD_STR(NC).IQCTLU, HYD_STR(NC).JQCTLU)
    IF( HYD_STR(NC).IQCTLD > 0 .AND. HYD_STR(NC).JQCTLD > 0 )THEN
      LCTLD(NC) = LIJ(HYD_STR(NC).IQCTLD, HYD_STR(NC).JQCTLD)
    ENDIF
  ENDDO

  ! *** SETUP L LOOKUPS FOR HYDRAULIC STRUCTURES
  ALLOCATE(LWRU(MAX(NQWR,1)), LWRD(MAX(NQWR,1)))
  LWRU = 0
  LWRD = 0
  DO NC=1,NQWR
    LWRU(NC) = LIJ(WITH_RET(NC).IQWRU, WITH_RET(NC).JQWRU)
    IF( WITH_RET(NC).IQWRD > 0 .AND. WITH_RET(NC).JQWRD > 0 )THEN
      LWRD(NC) = LIJ(WITH_RET(NC).IQWRD, WITH_RET(NC).JQWRD)
    ENDIF
  ENDDO

  RETURN
  999 CALL STOPP('DRIFTER.INP READING ERROR!')

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
!! DLON(L),L=2:LA ? : CELL CENTROID XCEL = XCOR(L,5)
!! DLAT(L),L=2:LA ? : CELL CENTROID YCEL = YCOR(L,5)
!
! @param[inout] BEDGEMOVE
! @param[in] XLA2 - x coordinates
! @param[in] YLA2 - y coordinates
! @param[in] ZLA2 - z coordinates
! @param[in] NI   - Current drifter particle
!
! Some other notes:
! *** INPUT:
! *** IF DEPOP=0: XLA,YLA,ZLA,XCOR(L,5),YCOR(L,5),BELV,HP
! *** IF DEPOP=1: XLA,YLA,    XCOR(L,5),YCOR(L,5),BELV,HP,DLA
! *** OUTPUT:
! ***             XLA,YLA,LLA(NP),KLA(NP),BELVLA(NP),HPLA(NP)
!---------------------------------------------------------------------------!

SUBROUTINE CONTAINER(BEDGEMOVE, XLA2, YLA2, ZLA2, NI)

  Implicit None

  LOGICAL(4),INTENT(INOUT) :: BEDGEMOVE 
  REAL(RKD) ,INTENT(IN)    :: XLA2 
  REAL(RKD) ,INTENT(IN)    :: YLA2 
  REAL(RKD) ,INTENT(IN)    :: ZLA2 
  INTEGER,INTENT(IN)       :: NI   

  ! *** Local variables
  INTEGER :: NPSTAT, NWR, NS, NPN, LLO, LI, NCTL, ipmc
  INTEGER :: LMILOC(1), K, L, ILN, JLN, LE, LN, LM, LLL, LLA2
  INTEGER :: LL, LU, LD
  REAL(RKD)  :: RADLA(LA), VELEK, VELNK, VPRO, XXO, YYO, DMIN, D, X0, Y0, X3, Y3, XP, YP, UTMPB, VTMPB, OFF, PMCTIME
  LOGICAL(4) :: LINSIDE

  ! *** DETERMINE THE NEAREST CELL CENTROID
  IF( JSPD(NI) == 1 )THEN
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
    IF( INSIDECELL(LLO,XLA(NI),YLA(NI)) )THEN
      LINSIDE = .TRUE.
    ENDIF
    
    IF( .NOT. LINSIDE .AND. LWC(LLO) > 0 .AND. LWC(LLO) <= LA )THEN
      IF( INSIDECELL(LWC(LLO),XLA(NI),YLA(NI)) )THEN
        LINSIDE = .TRUE.
        LLO = LWC(LLO)
      ENDIF
    ENDIF
    
    IF( .NOT. LINSIDE .AND. LEC(LLO) > 0 .AND. LEC(LLO) <= LA )THEN
      IF( INSIDECELL(LEC(LLO),XLA(NI),YLA(NI)) )THEN
        LINSIDE = .TRUE.
        LLO = LEC(LLO)
      ENDIF
    ENDIF
    
    IF( .NOT. LINSIDE .AND. LSC(LLO) > 0 .AND. LSC(LLO) <= LA )THEN
      IF( INSIDECELL(LSC(LLO),XLA(NI),YLA(NI)) )THEN
        LINSIDE = .TRUE.
        LLO = LSC(LLO)
      ENDIF
    ENDIF
    
    IF( .NOT. LINSIDE .AND. LNC(LLO) > 0 .AND. LNC(LLO) <= LA )THEN
      IF( INSIDECELL(LNC(LLO),XLA(NI),YLA(NI)) )THEN
        LINSIDE = .TRUE.
        LLO = LNC(LLO)
      ENDIF
    ENDIF
    
    ! *** DRIFTER IS IN THE CELL.  DETERMINE ALL OF THE SETTINGS
    LLA(NI) = LLO

    CALL DRF_DEPTH(LLA(NI),NI,BELVLA(NI),HPLA(NI))

    CALL DRF_LAYER(LLA(NI),BELVLA(NI),HPLA(NI),KLA(NI),ZLA(NI))

    !THE FIRST CALL: CONVERT DLA TO DLA/ZLA
    IF( DEPOP == 1 )THEN
      ZLA(NI) = HPLA(NI) + BELVLA(NI)-DLA(NI)
    ELSE
      DLA(NI) = MAX(HPLA(NI)+BELVLA(NI)-ZLA(NI),0._8)
    ENDIF
      
    ! *** Settings based on if oil spill calculations are being considers
    IF( ISOILSPI(LA_GRP(NI)) == 1 .AND. DDEN(LA_GRP(NI)) < 1000. .AND. GRPWS(LA_GRP(NI)) == 0. )THEN
      ! *** FORCE OIL TO SURFACE IF NOT USING SETTLING/RISING VELOCITIES
      DLA(NI) = 0.005
      ZLA(NI) = HPLA(NI) + BELVLA(NI) - DLA(NI)
    ENDIF

    LLA_Global(NI) = Map2Global(LLO).LG

    RETURN
  ENDIF

  ! *************************************  NORMAL PROCESSING  *********************************************
  ! *** CHECK IF OLD AND NEW POSITIONS ARE IN THE SAME CELL
  IF( INSIDECELL(LLA(NI), XLA(NI), YLA(NI)) )THEN
    CALL DRF_DEPTH(LLA(NI), NI, BELVLA(NI), HPLA(NI))
    CALL DRF_LAYER(LLA(NI), BELVLA(NI), HPLA(NI), KLA(NI), ZLA(NI))
    
    ! *** CHECK FOR BOUNDARY CONDITIONS
    IF( MOD(NITER,3) == 0 ) CALL CHECK_BCS
    
    RETURN
  ENDIF

  ! *************************************  OUTSIDE CELL  *********************************************
  ! *** EITHER THE DRIFTER IS OUTSIDE THE DOMAIN OR HAS CHANGED CELLS.  PROCESS THE CHANGE
  ! *** SAVE THE DRIFTER INFORMATION BEFORE ANY CHANGES BY EDGEMOVE
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
  DO LL = 1,9
    L = LADJ(LL,LLO)
    IF( L < 2 .OR. L > LA .OR. LL == 5 ) CYCLE

    IF( INSIDECELL(L,XLA(NI),YLA(NI)) )THEN
      ! *** PARTICLE IS INSIDE CURRENT CELL FOR NEXT CALL
      ! *** DEALING WITH THE WALLS
      LINSIDE = .TRUE.
      LI = LL
      IF( LL == 1 )THEN
        ! *** NORTHWEST
        IF( (SUB(LLO) > 0.5 .AND. SVB(L) > 0.5) .OR. (SVB(LADJ(2,LLO)) > 0.5 .AND. SUB(LADJ(2,LLO)) > 0.5) )THEN
          IF( ISOK(LL, LLO, K, L) )THEN
            LLA(NI) = L
            NPSTAT  = 1
          ENDIF
        ENDIF
      ELSEIF( LL == 2 .AND. SVB(LADJ(2,LLO)) > 0.5 )THEN
        ! *** NORTH
        IF( ISOK(LL, LLO, K, L) )THEN
          LLA(NI) = L
          NPSTAT  = 1
        ENDIF
      ELSEIF( LL == 3 )THEN  ! .AND. (SVB(LADJ(2,LLO)) > 0.5 .OR. SUB(LADJ(6,LLO)) > 0.5) )THEN
        ! *** NORTHEAST
        IF( (SUB(LADJ(6,LLO)) > 0.5 .AND. SVB(LADJ(3,LLO)) > 0.5) .OR. (SVB(LADJ(2,LLO)) > 0.5 .AND. SUB(LADJ(3,LLO)) > 0.5) )THEN
          IF( ISOK(LL, LLO, K, L) )THEN
            LLA(NI) = L
            NPSTAT  = 1
          ENDIF
        ENDIF
      ELSEIF( LL == 4 .AND. SUB(LLO) > 0.5 )THEN
        ! *** WEST
        IF( ISOK(LL, LLO, K, L) )THEN
          LLA(NI) = L
          NPSTAT  = 1
        ENDIF
      ELSEIF( LL == 6 .AND. SUB(LADJ(6,LLO)) > 0.5 )THEN
        ! *** EAST
        IF( ISOK(LL, LLO, K, L) )THEN
          LLA(NI) = L
          NPSTAT  = 1
        ENDIF
      ELSEIF( LL == 7 )THEN  ! .AND. (SVB(L) > 0.5 .OR. SUB(L) > 0.5) )THEN
        ! *** SOUTHWEST
        IF( (SUB(LLO) > 0.5 .AND. SVB(LADJ(4,LLO)) > 0.5) .OR. (SVB(LLO) > 0.5 .AND. SUB(LADJ(8,LLO)) > 0.5) )THEN
          IF( ISOK(LL, LLO, K, L) )THEN
            LLA(NI) = L
            NPSTAT  = 1
          ENDIF
        ENDIF
      ELSEIF( LL == 8 .AND. SVB(LLO) > 0.5 )THEN
        ! *** SOUTH
        IF( ISOK(LL, LLO, K, L) )THEN
          LLA(NI) = L
          NPSTAT  = 1
        ENDIF
      ELSEIF( LL == 9 )THEN  ! .AND. (SVB(L) > 0.5 .OR. SUB(LE) > 0.5) )THEN
        ! *** SOUTHEAST
        IF( (SUB(LADJ(6,LLO)) > 0.5 .AND. SVB(LADJ(6,LLO)) > 0.5) .OR. (SVB(LLO) > 0.5 .AND. SUB(LADJ(9,LLO)) > 0.5) )THEN
          IF( ISOK(LL, LLO, K, L) )THEN
            LLA(NI) = L
            NPSTAT  = 1
          ENDIF
        ENDIF
      ENDIF

      EXIT

    ELSE
      ! *** NOT IN CELL L, CHECK SPECIAL CONDITIONS
      
      ! Possible drifter got moved into a cell with a connector?
      IF( NS_Connectors(LLO) >= 1 .AND. NS_Connectors(L) >= 1 )THEN        ! *** N-S connection
        IF( (L == LNC(LLO) .AND. SVB3D(L,K) > 0.5) .OR. (L == LSC(LLO) .AND. SVB3D(LLO,K) > 0.5) )THEN
          VELEK = CUE(LLO)*U(LLO,KC) + CVE(LLO)*V(LLO,KC)
          VELNK = CUN(LLO)*U(LLO,KC) + CVN(LLO)*V(LLO,KC)
          VPRO  = VELEK*(XCOR(L,5)-XCOR(LLO,5)) + VELNK*(YCOR(L,5)-YCOR(LLO,5))
          IF( VPRO > 0 )THEN
            LLA(NI) = L
            XLA(NI) = XCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            YLA(NI) = YCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            NPSTAT  = 1
            LLA_Global(NI) = Map2Global(LLA(NI)).LG
            EXIT
          ENDIF
        ENDIF
      ELSEIF( EW_Connectors(LLO) >= 1 .AND. EW_Connectors(L) >= 1 )THEN    ! *** E-W connection
        IF( (L == LEC(LLO) .AND. SUB3D(L,K) > 0.5) .OR. (L == LWC(LLO) .AND. SUB3D(LLO,K) > 0.5) )THEN
          VELEK = CUE(LLO)*U(LLO,KC) + CVE(LLO)*V(LLO,KC)
          VELNK = CUN(LLO)*U(LLO,KC) + CVN(LLO)*V(LLO,KC)
          VPRO  = VELEK*(XCOR(L,5)-XCOR(LLO,5)) + VELNK*(YCOR(L,5)-YCOR(LLO,5))
          IF( VPRO > 0 )THEN
            LLA(NI) = L
            XLA(NI) = XCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            YLA(NI) = YCOR(L,5)   ! *** Since drifter not in a cell, move to destination cell centroid
            NPSTAT  = 1
            LLA_Global(NI) = Map2Global(LLA(NI)).LG
            EXIT
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  ENDDO
  
  ! *** CHECK IF THE PARTICLE IS INSIDE THE MODEL DOMAIN
  IF( NPSTAT == 0 )THEN
    ! *** CHECK IF THE PARTICLES SHOULD EXIT DOMAIN
    CALL CHECK_BCS

    ! *** IF LLO IS NOT A BC CELL THEN CONTINUE PROCESSING
    IF( LLA(NI) > 1 )THEN
      ! *** UPDATE POSITION BASED ON EDGE VELOCITIES (HALF OF CELL CENTERED)
      X0 = XLA(NI)
      Y0 = YLA(NI)

      ! **********************************************************************************
      D = 1E32
      IF( LINSIDE .AND. LI > 0 )THEN
        ! *** PARTICLE IS INSIDE AN INVALID CELL.  MOVE TO VALID
        LLL = LLO
        IF( LI == 2 )THEN
          OFF = 0.01*DYP(LLO)
          CALL DIST2LINE(LLO,2,X0,Y0,1,OFF,D,X3,Y3)     ! *** NORTH EDGE

        ELSEIF( LI == 4 )THEN
          OFF = 0.01*DXP(LLO)
          CALL DIST2LINE(LLO,1,X0,Y0,1,OFF,D,X3,Y3)     ! *** WEST EDGE

        ELSEIF( LI == 6 )THEN
          OFF = 0.01*DXP(LLO)
          CALL DIST2LINE(LLO,3,X0,Y0,1,OFF,D,X3,Y3)     ! *** EAST EDGE

        ELSEIF( LI == 8 )THEN
          OFF = 0.01*DYP(LLO)
          CALL DIST2LINE(LLO,4,X0,Y0,1,OFF,D,X3,Y3)     ! *** SOUTH EDGE

        ELSEIF( LL == 1 )THEN
          ! *** NORTHWEST
          IF( SUB3D(LLO,K) > 0.5 )THEN
            L = LADJ(4,LLO)
            OFF = 0.01*DYP(L)
            CALL DIST2LINE(L,2,X0,Y0,1,OFF,D,X3,Y3)     ! *** NORTH EDGE
            LLL = L
          ELSEIF( SVB3D(LADJ(2,LLO),K) > 0.5 )THEN
            L = LADJ(2,LLO)
            OFF = 0.01*DXP(L)
            CALL DIST2LINE(L,1,X0,Y0,1,OFF,D,X3,Y3)     ! *** WEST EDGE
            LLL = L
          ENDIF
        ELSEIF( LL == 3 )THEN
          ! *** NORTHEAST
          IF( SUB3D(LADJ(6,LLO),K) > 0.5 )THEN
            L = LADJ(6,LLO)
            OFF = 0.01*DYP(L)
            CALL DIST2LINE(L,2,X0,Y0,1,OFF,D,X3,Y3)     ! *** NORTH EDGE
            LLL = L
          ELSEIF( SVB3D(LADJ(2,LLO),K) > 0.5 )THEN
            L = LADJ(2,LLO)
            OFF = 0.01*DXP(L)
            CALL DIST2LINE(L,3,X0,Y0,1,OFF,D,X3,Y3)     ! *** EAST EDGE
            LLL = L
          ENDIF
        ELSEIF( LL == 7 )THEN
          ! *** SOUTHWEST
          IF( SUB3D(LLO,K) > 0.5 )THEN
            L = LADJ(4,LLO)
            OFF = 0.01*DYP(L)
            CALL DIST2LINE(L,4,X0,Y0,1,OFF,D,X3,Y3)     ! *** SOUTH EDGE
            LLL = L
          ELSEIF( SVB3D(LLO,K) > 0.5 )THEN
            L = LADJ(8,LLO)
            OFF = 0.01*DXP(L)
            CALL DIST2LINE(L,1,X0,Y0,1,OFF,D,X3,Y3)     ! *** WEST EDGE
            LLL = L
          ENDIF
        ELSEIF( LL == 9 )THEN
          ! *** SOUTHEAST
          IF( SUB3D(LADJ(6,LLO),K) > 0.5 )THEN
            L = LADJ(6,LLO)
            OFF = 0.01*DYP(L)
            CALL DIST2LINE(L,4,X0,Y0,1,OFF,D,X3,Y3)     ! *** SOUTH EDGE
            LLL = L
          ELSEIF( SVB3D(LLO,K) > 0.5 )THEN
            L = LADJ(8,LLO)
            OFF = 0.01*DXP(L)
            CALL DIST2LINE(L,3,X0,Y0,1,OFF,D,X3,Y3)     ! *** EAST EDGE
            LLL = L
          ENDIF
        ENDIF

        IF( ABS(D) < 1.E32 )THEN
          IF( INSIDECELL(LLL,X3,Y3) )THEN
            BEDGEMOVE = .TRUE.
            XLA(NI) = X3
            YLA(NI) = Y3
            LLA(NI) = LLL
          ELSE
            XLA(NI) = XXO
            YLA(NI) = YYO
            LLA(NI) = LLO
          ENDIF
        ELSE
          XLA(NI) = XXO
          YLA(NI) = YYO
          LLA(NI) = LLO
        ENDIF
        LLA_Global(NI) = Map2Global(LLA(NI)).LG

        CALL DRF_DEPTH(LLA(NI),NI,BELVLA(NI),HPLA(NI))
        CALL DRF_LAYER(LLA(NI),BELVLA(NI),HPLA(NI),KLA(NI),ZLA(NI))

        RETURN
      ENDIF

      ! **********************************************************************************
      ! *** PARTICLE IS OUTSIDE DOMAIN

      ! *** FIND CLOSEST FACE CONNECTED CELL (NEEDED TO ADDRESS SMALL GAPS BETWEEN CELLS)
      DMIN = 1E32
      LLL = 0
      DO LL=1,4
        IF( LL == 1 .OR. LL == 3 )THEN
          OFF = 0.01*DXP(LLO)
          CALL DIST2LINE(LLO,LL,X0,Y0,1,OFF,D,X3,Y3)
        ELSE
          OFF =  0.01*DYP(LLO)
          CALL DIST2LINE(LLO,LL,X0,Y0,1,OFF,D,X3,Y3)
        ENDIF
        IF( ABS(D) < ABS(DMIN) .AND. ABS(D) < (25.*OFF) )THEN
          DMIN = D
          LLL = LL
          XP = X3
          YP = Y3
        ENDIF
      ENDDO

      IF( LLL /= 0 )THEN
        ! *** FOUND AN EDGE, SEE IF IT IS ADJACENT TO AN ACTIVE CELL
        L = 0
        IF( LLL == 2 )THEN
          L = LNC(LLO)             ! *** NORTH
          IF( SVB3D(L,K) > 0.5 )THEN
            LLA(NI) = L
            NPSTAT = 1
          ENDIF
        ELSEIF( LLL == 3 )THEN
          L = LEC(LLO)             ! *** EAST
          IF( SUB3D(L,K) > 0.5 )THEN
            LLA(NI) = L
            NPSTAT = 1
          ENDIF
        ELSEIF( LLL == 4 )THEN
          L = LLO                  ! *** SOUTH
          IF( SVB3D(L,K) > 0.5 )THEN
            LLA(NI) = LSC(L)
            NPSTAT = 1
          ENDIF
        ELSEIF( LLL == 1 )THEN
          L = LLO                  ! *** WEST
          IF( SUB3D(L,K) > 0.5 )THEN
            LLA(NI) = LWC(L)
            NPSTAT = 1
          ENDIF
        ENDIF

        IF( NPSTAT == 0 )THEN
          ! *** MOVE THE PARTICLE BACK INTO THE ORIGINAL CELL
          LLA(NI) = LLO
          NPSTAT = 1
          XLA(NI) = XP
          YLA(NI) = YP
          BEDGEMOVE = .TRUE.

        ENDIF
        LLA_Global(NI) = Map2Global(LLA(NI)).LG
      ENDIF

      IF( NPSTAT == 0 )THEN
        ! *** PARTICLE MUST BE IN ONE OF THE CORNER CELLS
        DMIN = 1E32
        LLL = 1
        DO LL=1,4
          D = SQRT((XLA(NI)-XCOR(LLO,LL))**2 + (YLA(NI)-YCOR(LLO,LL))**2)
          IF( D < DMIN )THEN
            DMIN = D
            LLL = LL
          ENDIF
        ENDDO

        ! *** PLACE THE PARTICLE IN THE ASSOCIATED CORNER
        BEDGEMOVE = .TRUE.
        NPSTAT = 1
        OFF = 0.01
        X0 = OFF*DXP(LLO)
        Y0 = OFF*DYP(LLO)
        IF( LLL == 1 )THEN
          XLA(NI) = XCOR(LLO,LLL) - (-X0 * CUE(LLO) - Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (-X0 * CUN(LLO) - Y0 * CVN(LLO))
        ELSEIF( LLL == 2 )THEN
          XLA(NI) = XCOR(LLO,LLL) - (-X0 * CUE(LLO) + Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (-X0 * CUN(LLO) + Y0 * CVN(LLO))
        ELSEIF( LLL == 3 )THEN
          XLA(NI) = XCOR(LLO,LLL) - (X0 * CUE(LLO) + Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (X0 * CUN(LLO) + Y0 * CVN(LLO))
        ELSEIF( LLL == 4 )THEN
          XLA(NI) = XCOR(LLO,LLL) - (X0 * CUE(LLO) - Y0 * CVE(LLO))
          YLA(NI) = YCOR(LLO,LLL) - (X0 * CUN(LLO) - Y0 * CVN(LLO))
        ELSE
          LLA(NI) = 1
          LLA_Global(NI) = 1
          NPSTAT = 0
        ENDIF
      ENDIF
    ENDIF

  ELSE    ! IF( NPSTAT == 1 )THEN
    ! *** DRIFTER IS NOW IN A NEW CELL.  CHECK BOUNDARY CONDITIONS
    CALL CHECK_BCS
  ENDIF

  ! *** DETERMINE BOTTOM ELEVATION AND TOTAL WATER DEPTH OF DRIFTERS FOR EVERY TIMESTEP
  IF( LLA(NI) >= 2 )THEN
    LLA_Global(NI) = Map2Global(LLA(NI)).LG
    CALL DRF_DEPTH(LLA(NI),NI,BELVLA(NI),HPLA(NI))
    CALL DRF_LAYER(LLA(NI),BELVLA(NI),HPLA(NI),KLA(NI),ZLA(NI))
  ENDIF

! ***************************************************************************
CONTAINS

!---------------------------------------------------------------------------!
!> @details Checks if drifter cell location is within a boundary condition
!! cell.  (For MPI) checks if the cell has entered the ghost region
!---------------------------------------------------------------------------!
SUBROUTINE CHECK_BCS

  Implicit None
  REAL :: RAD
  
  ! *** LLA(NI) is the L index of the last valid cell the drifter was in.
  
  ! *** Checks if cell id is valid
  IF( LLA(NI) < 2 .OR. LLA(NI) > LA )THEN
    RETURN
  ENDIF

  ! *** Check if drifter is within the North, South, East, West open boundary (Only for drifters outside the domain)
  IF( NBCSOP > 0 .AND. NPSTAT == 0 )THEN
    IF( ANY(LPBN == LLA(NI)) .OR. ANY(LPBS == LLA(NI)) .OR.  &
        ANY(LPBE == LLA(NI)) .OR. ANY(LPBW == LLA(NI)) )   THEN

      CALL SET_DRIFTER_OUT(XLA2,YLA2,ZLA2)

      IF( NPD < 10000 .OR. DEBUG )THEN
        PRINT '(A36,I8)','OPEN BOUNDARY, DRIFTER IS OUTSIDE:',NI
      ELSE
        NSTOPPED(4) = NSTOPPED(4) + 1
      ENDIF
      RETURN
      
    ENDIF
  ENDIF
  
  ! *** Check if drifter in outflow boundary condition
  IF( NQSIJ > 0 )THEN
    IF( ANY(LQS == LLA(NI)) .AND. QSUM(LLA(NI),KLA(NI)) < 0. )THEN
      ! *** If inside a cell then use cell dimensions to determine if should remove
      IF( NPSTAT == 1 )THEN
        ! *** INSIDE CELL
        RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
      ELSE
        RAD = 0.
      ENDIF
            
      IF( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )THEN
        CALL SET_DRIFTER_OUT(XLA2,YLA2,ZLA2)

        IF( NPD < 10000 .OR. DEBUG )THEN
          PRINT '(A36,I8)','WITHDAWAL CELL, DRIFTER IS OUTSIDE:',NI
        ELSE
          NSTOPPED(5) = NSTOPPED(5) + 1
        ENDIF
      ENDIF    
      RETURN
      
    ENDIF
  ENDIF
  
  ! *** Check if drifter is in a withdrawal cell
  IF( NQWR > 0 )THEN
    IF( ANY(LWRU == LLA(NI)) .OR. ANY(LWRD == LLA(NI)) )THEN
      ! *** PMC - MODIFIED TO HANDLE +/- FLOWS FOR THE W/R BOUNDARY CONDITION
      LLA2 = LLA(NI)
      DO NWR = 1,NQWR
        ! *** Handle +/- Flows for Withdrawal/Return Structures
        NS = WITH_RET(NWR).NQWRSERQ
        IF( QWRSERT(NS) >= 0. )THEN
          ! *** Original Withdrawal/Return
          IF( WITH_RET(NWR).IQWRU == IL(LLA2) .AND. WITH_RET(NWR).JQWRU == JL(LLA2) .AND. WITH_RET(NWR).KQWRU == KLA(NI) )THEN
            ! *** If inside a cell then use cell dimensions to determine if should remove
            IF( NPSTAT == 1 )THEN
              ! *** INSIDE CELL
              RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
            ELSE
              RAD = 0.
            ENDIF
            
            IF( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )THEN
              ! ***  RETURN DRIFTER
              LLA(NI)   = LIJ(WITH_RET(NWR).IQWRD,WITH_RET(NWR).JQWRD)
              LLA2      = LLA(NI)
              XLA(NI)   = XCOR(LLA2,5)
              YLA(NI)   = YCOR(LLA2,5)
              ZLA(NI)   = BELV(LLA2) + HP(LLA2)*ZZ(LLA2,WITH_RET(NWR).KQWRD)
              KLA(NI)   = WITH_RET(NWR).KQWRD
              BEDGEMOVE = .TRUE.
            ENDIF
            EXIT
          ENDIF
        ELSE
          ! *** Reverse Flow Withdrawal/Return
          IF( WITH_RET(NWR).IQWRD == IL(LLA2) .AND. WITH_RET(NWR).JQWRD == JL(LLA2) .AND. WITH_RET(NWR).KQWRD == KLA(NI) )THEN
            ! *** If inside a cell then use cell dimensions to determine if should remove
            IF( NPSTAT == 1 )THEN
              ! *** INSIDE CELL
              RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
            ELSE
              RAD = 0.
            ENDIF
            
            IF( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )THEN
              ! ***  RETURN DRIFTER
              LLA(NI)   = LIJ(WITH_RET(NWR).IQWRU,WITH_RET(NWR).JQWRU)
              LLA2      = LLA(NI)
              XLA(NI)   = XCOR(LLA2,5)
              YLA(NI)   = YCOR(LLA2,5)
              ZLA(NI)   = BELV(LLA2) + HP(LLA2)*ZZ(LLA2,WITH_RET(NWR).KQWRU)
              KLA(NI)   = WITH_RET(NWR).KQWRU
              BEDGEMOVE = .TRUE.
            ENDIF
            EXIT
          ENDIF
        ENDIF
      ENDDO

      ! ***
      IF( LLA(NI) == 1 )THEN
        ! *** Should never get here
        IF( NPD < 10000 .OR. DEBUG )THEN
          PRINT '(A36,I8)','WITHDRAWAL/RETURN, DRIFTER IS OUTSIDE:',NI
        ELSE
          NSTOPPED(6) = NSTOPPED(6) + 1
        ENDIF
      ENDIF
      RETURN
      
    ENDIF
  ENDIF     ! *** END OF NQWR

  IF( NQCTL > 0 )THEN
    IF( ANY(LCTLU == LLA(NI)) )THEN
      ! *** HYDRAULIC STRUCTURE.  RETURN DRIFTER TO DOWNSTREAM CELL, IF ANY
      LLA2 = LLA(NI)
      DO NCTL=1,NQCTL
        LU = LIJ(HYD_STR(NCTL).IQCTLU, HYD_STR(NCTL).JQCTLU)
        IF( LU == LLA2 )THEN
          ! *** If inside a cell then use cell dimensions to determine if should remove
          IF( NPSTAT == 1 )THEN
            ! *** INSIDE CELL
            RAD = SQRT( (XLA(NI)-XCOR(LLA2,5))**2 + (YLA(NI)-YCOR(LLA2,5))**2 )
          ELSE
            RAD = 0.
          ENDIF
            
          IF( RAD < MIN(DXP(LLA2),DYP(LLA2))*0.25 )THEN
            LD = LIJ(MAX(HYD_STR(NCTL).IQCTLD,1), MAX(HYD_STR(NCTL).JQCTLD,1))
            IF( LD >= 2 )THEN
              LLA(NI) = LD
              XLA(NI) = XCOR(LLA(NI),5)
              YLA(NI) = YCOR(LLA(NI),5)
              ZLA(NI) = BELV(LD) + HP(LD)
              BEDGEMOVE = .TRUE.
              EXIT
            ELSE
              IF( NPD < 10000 .OR. DEBUG )THEN
                PRINT '(A40,I8)','HYDRAULIC STRUCTURE, DRIFTER IS OUTSIDE:',NI
              ELSE
                NSTOPPED(3) = NSTOPPED(3) + 1
              ENDIF
              CALL SET_DRIFTER_OUT(XLA2, YLA2, ZLA2)
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    RETURN
  ENDIF
  
END SUBROUTINE CHECK_BCS

!---------------------------------------------------------------------------!
!> @details Sets drifter cell to '1' so that it is effectively invalid
!! and updates the x,y,z locations of the drifter
!> @param[in] XLA2
!> @param[in] YLA2
!> @param[in] ZLA2
!---------------------------------------------------------------------------!
SUBROUTINE SET_DRIFTER_OUT(XLA2,YLA2,ZLA2)

  Implicit None

  REAL(RKD) ,INTENT(IN)   :: XLA2 !<
  REAL(RKD) ,INTENT(IN)   :: YLA2 !<
  REAL(RKD) ,INTENT(IN)   :: ZLA2 !<

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

  Implicit None
  
  INTEGER, INTENT(IN) :: LL, LLO, K, L
  INTEGER :: LLB
  LOGICAL :: FLAG
  
  IF( NBLOCKED > 0 )THEN
    FLAG = .FALSE.
    DO LLB = 1,NBLOCKED
      IF( LBLOCKED(LLB) == L .OR. LBLOCKED(LLB) == LLO )THEN
        IF( LL == 1 )THEN
          ! *** NORTHWEST
          IF( (SUB3D(LLO,K) > 0.5 .AND. SVB3D(L,K) > 0.5) .OR. (SVB3D(LADJ(2,LLO),K) > 0.5 .AND. SUB3D(LADJ(2,LLO),K) > 0.5) )THEN
            FLAG = .TRUE.
          ENDIF
        ELSEIF( LL == 2 .AND. SVB3D(LADJ(2,LLO),K) > 0.5 )THEN
          ! *** NORTH
            FLAG = .TRUE.
        ELSEIF( LL == 3 )THEN  ! .AND. (SVB3D(LADJ(2,LLO),K) > 0.5 .OR. SUB3D(LADJ(6,LLO),K) > 0.5) )THEN
          ! *** NORTHEAST
          IF( (SUB3D(LADJ(6,LLO),K) > 0.5 .AND. SVB3D(LADJ(3,LLO),K) > 0.5) .OR. (SVB3D(LADJ(2,LLO),K) > 0.5 .AND. SUB3D(LADJ(3,LLO),K) > 0.5) )THEN
            FLAG = .TRUE.
          ENDIF
        ELSEIF( LL == 4 .AND. SUB3D(LLO,K) > 0.5 )THEN
          ! *** WEST
            FLAG = .TRUE.
        ELSEIF( LL == 6 .AND. SUB3D(LADJ(6,LLO),K) > 0.5 )THEN
          ! *** EAST
            FLAG = .TRUE.
        ELSEIF( LL == 7 )THEN  ! .AND. (SVB3D(L,K) > 0.5 .OR. SUB3D(L,K) > 0.5) )THEN
          ! *** SOUTHWEST
          IF( (SUB3D(LLO,K) > 0.5 .AND. SVB3D(LADJ(4,LLO),K) > 0.5) .OR. (SVB3D(LLO,K) > 0.5 .AND. SUB3D(LADJ(8,LLO),K) > 0.5) )THEN
            FLAG = .TRUE.
          ENDIF
        ELSEIF( LL == 8 .AND. SVB3D(LLO,K) > 0.5 )THEN
          ! *** SOUTH
            FLAG = .TRUE.
        ELSEIF( LL == 9 )THEN  ! .AND. (SVB3D(L,K) > 0.5 .OR. SUB3D(LE,K) > 0.5) )THEN
          ! *** SOUTHEAST
          IF( (SUB3D(LADJ(6,LLO),K) > 0.5 .AND. SVB3D(LADJ(6,LLO),K) > 0.5) .OR. (SVB3D(LLO,K) > 0.5 .AND. SUB3D(LADJ(9,LLO),K) > 0.5) )THEN
            FLAG = .TRUE.
          ENDIF
        ENDIF
        ISOK = FLAG
        RETURN
      ENDIF
    ENDDO
  ENDIF
  
  ! *** ALLOW TRANSPORT FOR ALL LAYERS EXCEPT FOR DRIFTERS < KSZ(LLO)-2
  IF( K >= KSZ(L)-1 )THEN
    ISOK = .TRUE.
  ELSE
    ISOK = .FALSE.
  ENDIF

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

  INTEGER,INTENT(IN )    :: LNI !< Cell the drifter belongs to
  INTEGER,INTENT(IN )    :: KNI !< Layer the drifter belongs to
  INTEGER,INTENT(IN )    :: NI  !< Drifter number
  REAL(RKD) ,INTENT(OUT) :: U2NI !< Interpoled U velocity
  REAL(RKD) ,INTENT(OUT) :: V2NI !< Interpoled V velocity
  REAL(RKD) ,INTENT(OUT) :: W2NI !< Interpoled W velocity
  ! *** Local variables
  INTEGER :: L,LE,LN,K1,K2,K3,KZ1,KZ2,LL
  REAL(RKD) :: RAD2,SU1,SU2,SU3,SV1,SV2,SV3,SW1,SW2,SW3,XOUT,YOUT
  REAL(RKD) :: UTMPB,VTMPB,UTMPB1,VTMPB1,WTMPB,WTMPB1
  REAL(RKD) :: VELEK,VELNK,VELEK1,VELNK1,ZSIG,DZCTR,DZTOP
  REAL(RKD) :: UKB,UKT,VKB,VKT,UKB1,UKT1,VKB1,VKT1,VFACTOR,VE,UE

  ! *** ORDER OF CELLS
  ! ***   1  2  3
  ! ***   4  5  6
  ! ***   7  8  9
  VFACTOR = 0.5
  SU1=0
  SU2=0
  SU3=0
  SV1=0
  SV2=0
  SV3=0
  SW1=0
  SW2=0
  SW3=0

  ! *** UPDATE VERTICAL POSITION FOR TOP/BOTTOM LAYERS
  ZSIG = (ZLA(NI)-BELVLA(NI))/HPLA(NI)
  ZSIG = MAX(Z(LNI,0),MIN(Z(LNI,KC),ZSIG))
  ZLA(NI) = ZSIG*HPLA(NI)+BELVLA(NI)
  IF( ZSIG >= ZZ(LNI,KNI) )THEN
    ! *** DRIFTER IS ABOVE THE LAYER MID DEPTH
    K1 = KNI
    K2 = MIN(KNI+1,KC)
    KZ1= KNI
    KZ2= KNI+1
  ELSE
    ! *** DRIFTER IS BELOW THE LAYER MID DEPTH
    K1 = MAX(KSZ(LNI),KNI-1)
    K2 = KNI
    KZ1= KNI-1
    KZ2= KNI
  ENDIF

  ! *** Go over the 8 cells around the active
  DO LL=1,9
    L = LADJ(LL,LNI)
    
    ! *** Drifter is active and within a valid cell
    IF( L >= 2 .AND. L <= LA )THEN
      !IF( KSZ(L) > KNI ) CYCLE
      LE = LEC(L)
      LN = LNC(L)

      ! *** CALCULATING HORIZONTAL VELOCITY COMPONENTS AT CENTROID-BOTTOM
      UKB  = 0.5*STCUV(L)*( RSSBCE(L)*U2(LE,K1) + RSSBCW(L)*U2(L,K1) )
      VKB  = 0.5*STCUV(L)*( RSSBCN(L)*V2(LN,K1) + RSSBCS(L)*V2(L,K1) )

      ! *** CALCULATING HORIZONTAL VELOCITY COMPONENTS AT CENTROID-TOP
      IF( K2 < KSZ(L) )THEN
        K3 = MIN(K2+1,KSZ(L))      ! *** USE BOTTOM ACTIVE LAYER IF ONE LAYER HIGHER
      ELSE
        K3 = K2
      ENDIF
      UKT  = 0.5*STCUV(L)*( RSSBCE(L)*U2(LE,K3) + RSSBCW(L)*U2(L,K3) )
      VKT  = 0.5*STCUV(L)*( RSSBCN(L)*V2(LN,K3) + RSSBCS(L)*V2(L,K3) )

      ! *** DISTANCE FROM CENTROID SQUARED
      RAD2 = MAX((XLA(NI)-XCOR(L,5))**2 + (YLA(NI)-YCOR(L,5))**2,1D-8)

      DZCTR = ZCTR(L,KZ2)-ZCTR(L,KZ1)
      IF( DZCTR > 1D-8 )THEN
        ! *** INTERPOLATION FOR HORIZONTAL VELOCITY COMPONENT
        UTMPB  = (UKT -UKB )*(ZSIG-ZCTR(L,KZ1))/DZCTR+UKB
        VTMPB  = (VKT -VKB )*(ZSIG-ZCTR(L,KZ1))/DZCTR+VKB
      ELSE
        UTMPB  = UKT
        VTMPB  = VKT
      ENDIF

      DZTOP = Z(L,KNI)-Z(L,KNI-1)
      IF( DZTOP > 1D-8 )THEN
        ! *** INTERPOLATION FOR VERTICAL VELOCITY COMPONENT
        WTMPB  = (W2(L,KNI)-W2(L,KNI-1))*(ZSIG-Z(L,KNI-1))/DZTOP + W2(L,KNI-1)
      ELSE
        WTMPB  = W2(L,KNI)
      ENDIF

    ELSE
      ! *** EDGE CELLS, APPLY A ZERO FACE VELOCITIES

      ! *** ORDER OF CELLS
      ! ***   1  2  3
      ! ***   4  5  6
      ! ***   7  8  9

      UTMPB = SLIPFACTOR*U2(LNI,KNI)
      VTMPB = SLIPFACTOR*V2(LNI,KNI)
      WTMPB = SLIPFACTOR*W2(LNI,KNI)

      ! *** PMC NOTE - THESE COORDINATES COULD BE CALCULATED ONCE AND SAVED FOR LATER USE
      IF( LL == 1 )THEN
        ! *** NORTHWEST
        XOUT = XCOR(LNI,5) + ( DYP(LNI)*CVE(LNI) - DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) + ( DYP(LNI)*CVN(LNI) - DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      ELSEIF( LL == 2 )THEN
        ! *** NORTH FACE (MOVE WITH DRIFTER)
        XOUT = XLA(NI) + 0.5*( DYP(LNI)*CVE(LNI) )
        YOUT = YLA(NI) + 0.5*( DYP(LNI)*CVN(LNI) )
        RAD2 = MAX( (XLA(NI)-XOUT)**2 + (YLA(NI)-YOUT)**2, 1D-8)
        VTMPB = -VFACTOR*MAX(V2(LNI,KNI),0.0)
      ELSEIF( LL == 3 )THEN
        ! *** NORTHEAST
        XOUT = XCOR(LNI,5) + ( DYP(LNI)*CVE(LNI) + DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) + ( DYP(LNI)*CVN(LNI) + DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      ELSEIF( LL == 4 )THEN
        ! *** WEST FACE (MOVE WITH DRIFTER)
        XOUT = XLA(NI) - 0.5*( DXP(LNI)*CUE(LNI) )
        YOUT = YLA(NI) - 0.5*( DXP(LNI)*CUN(LNI) )
        RAD2 = MAX( (XLA(NI)-XOUT)**2 + (YLA(NI)-YOUT)**2, 1D-8)
        UTMPB = -VFACTOR*MIN(U2(LEC(LNI),KNI),0.0)
      ELSEIF( LL == 6 )THEN
        ! *** EAST FACE (MOVE WITH DRIFTER)
        XOUT = XLA(NI) + 0.5*( DXP(LNI)*CUE(LNI) )
        YOUT = YLA(NI) + 0.5*( DXP(LNI)*CUN(LNI) )
        RAD2 = MAX( (XLA(NI)-XOUT)**2 + (YLA(NI)-YOUT)**2, 1D-8)
        UTMPB = -VFACTOR*MAX(U2(LNI,KNI),0.0)
      ELSEIF( LL == 7 )THEN
        ! *** SOUTHWEST
        XOUT = XCOR(LNI,5) - ( DYP(LNI)*CVE(LNI) + DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) - ( DYP(LNI)*CVN(LNI) + DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      ELSEIF( LL == 8 )THEN
        ! *** SOUTH  (MOVE WITH DRIFTER)
        !XOUT = XCOR(LNI,5) - ( DYP(LNI)*CVE(LNI) )
        !YOUT = YCOR(LNI,5) - ( DYP(LNI)*CVN(LNI) )
        XOUT = XLA(NI) - 0.5*( DYP(LNI)*CVE(LNI) )
        YOUT = YLA(NI) - 0.5*( DYP(LNI)*CVN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
        VTMPB  = -VFACTOR*MIN(V2(LNC(LNI),KNI),0.0)
      ELSEIF( LL == 9 )THEN
        ! *** SOUTHEAST
        XOUT = XCOR(LNI,5) - ( DYP(LNI)*CVE(LNI) - DXP(LNI)*CUE(LNI) )
        YOUT = YCOR(LNI,5) - ( DYP(LNI)*CVN(LNI) - DXP(LNI)*CUN(LNI) )
        RAD2 = MAX((XLA(NI)-XOUT)**2+(YLA(NI)-YOUT)**2,1D-8)
      ENDIF
      L = LNI
    ENDIF

    ! ** ROTATION
    VELEK  = CUE(L)*UTMPB  + CVE(L)*VTMPB
    VELNK  = CUN(L)*UTMPB  + CVN(L)*VTMPB

    ! *** APPLY WEIGHTING
    SU2 = SU2 + VELEK/RAD2
    SU3 = SU3 + 1._8/RAD2
    SV2 = SV2 + VELNK/RAD2
    SV3 = SV3 + 1._8/RAD2
    SW2 = SW2 + WTMPB/RAD2
    SW3 = SW3 + 1._8/RAD2

  ENDDO

  U2NI = SU2/SU3
  V2NI = SV2/SV3
  W2NI = SW2/SW3

END SUBROUTINE DRF_VELOCITY

SUBROUTINE RANDCAL(L,K,NP)
  INTEGER,INTENT(IN) :: L,K,NP
  !REAL(8),EXTERNAL :: DRAND -Return double-precision random numbers in the range 0.0 to 1.0, not including the end points.
  REAL(RKD) :: COEF

  ! *** HORIZONTAL
  IF( LA_PRAN == 1 .OR. LA_PRAN == 3 )THEN
    IF( LA_DIFOP == 0 )THEN
      COEF = SQRT(2.*AH(L,K)*DELTD)*LA_HORDIF
    ELSE
      COEF = DIFFH
    ENDIF
    XLA(NP) = XLA(NP) + (2.*DRAND(0)-1.)*COEF
    YLA(NP) = YLA(NP) + (2.*DRAND(0)-1.)*COEF
  ENDIF

  ! *** VERTICAL (IF LA_ZCAL=1 THEN FULLY 3D)
  IF( LA_PRAN >= 2 .AND. LA_ZCAL == 1 )THEN
    IF( LA_DIFOP == 0 )THEN
      COEF = SQRT(2.*AV(L,K)*DELTD)*LA_VERDIF
    ELSE
      COEF = DIFFV
    ENDIF
    ZLA(NP) = ZLA(NP) + (2*DRAND(0)-1)*COEF
  ENDIF
  
END SUBROUTINE

! **** APPLY WIND SHEAR FOR DRIFTERS AT THE SURFACE
SUBROUTINE WIND_DRIFT(L,K,NP)

  INTEGER,INTENT(IN) :: L,K,NP
  REAL(RKD) :: WSPD, COEF, XCOEF, YCOEF
  
  IF( IOSWD == 1 )THEN
    ! *** APPLY SHEAR TO THE SURFACE OIL BASED ON 
    WSPD = SQRT(WNDVELE(L)**2 + WNDVELN(L)**2)
    COEF = OSWDA * WSPD + OSWDB 
  
    XLA(NP) = XLA(NP) + COEF * WNDVELE(L) * DELTD 
    YLA(NP) = YLA(NP) + COEF * WNDVELN(L) * DELTD
    
  ELSEIF( IOSWD == 2 )THEN
    XCOEF = OSWDA * WNDVELE(L) + OSWDB
    YCOEF = OSWDA * WNDVELN(L) + OSWDB
  
    XLA(NP) = XLA(NP) + XCOEF * WNDVELE(L) * DELTD
    YLA(NP) = YLA(NP) + YCOEF * WNDVELN(L) * DELTD
  ENDIF
  
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

  Implicit None

  INTEGER,INTENT(IN) :: LNI        !< Cell index 'L' for the drifter
  INTEGER,INTENT(IN) :: NI         !< drifter number
  REAL(RKD),INTENT(OUT) :: BELVNI  !< Interpolated bed elevation
  REAL(RKD),INTENT(OUT) :: HPNI    !< Interpolated water surface elevation

  ! *** Local variables
  INTEGER :: L,LL
  REAL(RKD) :: BELVNI1,BELVNI2,RAD2,ZETA

  BELVNI1 = 0.
  BELVNI2 = 0.
  ZETA = 0.

  DO LL=1,9
    L = LADJ(LL,LNI)
    IF( L >= 2 .AND. L <= LA )THEN
      RAD2 = MAX((XLA(NI)-XCOR(L,5))**2 + (YLA(NI)-YCOR(L,5))**2,1D-8)         ! *** INVERSE DISTANCE SQUARED
      !RAD2 = SQRT(MAX((XLA(NI)-XCOR(L,5))**2 + (YLA(NI)-YCOR(L,5))**2,1D-8))  ! *** LINEAR INTERPOLATION
      BELVNI1 = BELVNI1 + BELV(L)/RAD2
      BELVNI2 = BELVNI2 + 1._8/RAD2
      ZETA = ZETA + (HP(L)+BELV(L))/RAD2
    ENDIF
  ENDDO
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

  Implicit None
  ! *** Dummy variables
  INTEGER,INTENT(IN)    :: LNI    !< Current cell 'L' index the drifter belongs to
  REAL(RKD), INTENT(IN) :: BELVNI !< Interpolated bed elevation
  REAL(RKD), INTENT(IN) :: HPNI   !< Interpolated water surface elevation
  INTEGER,INTENT(OUT)   :: KLN    !< Layer this drifter belongs in
  REAL(RKD), INTENT(INOUT) :: ZLN !< Vertical location of the drifter based on interpolated values

  ! *** Local variables
  INTEGER :: K
  REAL(RKD) :: ZSIG

  ! *** Check if the bottom most active layer (for the current cell) is equal to the top layer
  IF( KSZ(LNI) == KC )THEN  ! *** Alberta
    KLN=KC
    ZSIG = (ZLN-BELVNI)/HPNI
    ZSIG = MAX(Z(LNI,0),MIN(Z(LNI,KC),ZSIG))
    ZLN  = ZSIG*HPNI + BELVNI
  ELSEIF( LNI >= 2 .AND. KC > 1 )THEN
    ZSIG = (ZLN - BELVNI)/HPNI
    ZSIG = MAX( Z(LNI,0), MIN(Z(LNI,KC),ZSIG) )
    ZLN = ZSIG*HPNI + BELVNI
    ! *** From the current active layer for cell the drifter belongs to the top layer
    DO K=KSZ(LNI),KC
      ! *** Based on the new vertical height determine the layer the drifter belongs in
      IF( ZSIG  >= Z(LNI,K-1) .AND. ZSIG  <= Z(LNI,K) )THEN
        KLN = K
        EXIT
      ENDIF
    ENDDO
  ENDIF

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

  INTEGER,INTENT(IN )    :: LNI !< Cell index the drifter is in
  INTEGER,INTENT(IN )    :: KNI !< Layer index the drifer is in
  INTEGER,INTENT(IN )    :: NI !< Drifter number
  REAL(RKD) ,INTENT(OUT) :: DAHX !< Horizontal diffusion coeffecient in x-direction
  REAL(RKD) ,INTENT(OUT) :: DAHY !< horizontal diffusion coeffecient in y-direction
  REAL(RKD) ,INTENT(OUT) :: DAVZ !< horizontal diffusion coeffecient in z-direction

  ! *** Local variables
  REAL(RKD) :: DAHXC, DAHYC, ZSIG, WEIC, DAVB, DAVT

  DAHXC = SUB3D(LNI,KNI)*( AH(LNI,KNI) - AH(LWC(LNI),KNI) )/DXU(LNI)
  DAHYC = SVB3D(LNI,KNI)*( AH(LNI,KNI) - AH(LSC(LNI),KNI) )/DYV(LNI)

  ! ***  HORIZONTAL DIFFUSION COMPONENT (M/S)
  DAHX = CUE(LNI)*DAHXC + CVE(LNI)*DAHYC
  DAHY = CUN(LNI)*DAHXC + CVN(LNI)*DAHYC

  ! *** APPLY A FACTOR TO GET THE HORIZONTAL DIFFUSION COMPONENT (M/S)
  DAHX = DAHX*LA_HORDIF
  DAHY = DAHY*LA_HORDIF

  IF( KC <= 2 .OR. LA_ZCAL == 0 )THEN
    DAVZ = 0.
  ELSE
    ! *** AV IS THE VERTICAL DIFFUSION COEFFICIENT AT THE TOP OF THE LAYER
    ! *** AV(L,KC) = 0 AND AV(L,0) = 0
    ZSIG = (ZLA(NI)-BELVLA(NI))/HPLA(NI)
    ZSIG = MAX(Z(LNI,0),MIN(Z(LNI,KC),ZSIG))
    ZLA(NI)=ZSIG*HPLA(NI)+BELVLA(NI)
    WEIC = Z(LNI,KC)-(Z(LNI,KNI)-ZSIG)*DZIC(LNI,KNI)  !WEIGHTING COEFFICIENT

    IF( KNI == KC )THEN
      DAVT = 0.
      DAVB = (AV(LNI,KS)-AV(LNI,KS-1))/HPK(LNI,KNI)
    ELSEIF (KNI == 1 )THEN
      DAVT = (AV(LNI,KNI+1)-AV(LNI,KNI))/HPK(LNI,KNI)
      DAVB = 0.
    ELSE
      DAVT = 2._RKD*(AV(LNI,KNI+1)-AV(LNI,KNI  ))/(DZC(LNI,KNI+1)+DZC(LNI,KNI  ))/HP(LNI)
      DAVB = 2._RKD*(AV(LNI,KNI  )-AV(LNI,KNI-1))/(DZC(LNI,KNI  )+DZC(LNI,KNI-1))/HP(LNI)
    ENDIF
    DAVZ = (DAVT-DAVB)*WEIC+DAVB

    ! *** APPLY A VERTICAL GRADIENT FACTOR
    DAVZ = DAVZ*LA_VERDIF

  ENDIF

END SUBROUTINE DIFGRAD

SUBROUTINE OIL_PROC(NP)
  ! ** CALCULATING  OIL REDUCTION DUE TO:
  ! ** EVAPORATION AND BIODEGRADATION
  ! ** OUTPUT: DVOL(NP)
  INTEGER,INTENT(IN ) :: NP
  INTEGER :: NG,LNP
  REAL(RKD) :: DFV,RKA,TK,RGAS,DTHE
  REAL(RKD) :: HENRY,OILMAS,BIODEG,VWIND,TDEP,TNP
  DATA RGAS /8.314_RKD/

  ! ** OIL EVAPORATION
  ! ** BASED ON THE PAPER:
  ! ** EVAPORATION RATE OF SPILLS OF HYDROCARBONS AND PETROLEUM MIXTURES
  ! ** BY WARREN STLVER AND DONALD MACKAY
  ! ** VAPOR PRESSURE IS REFERED FROM: 2013 CRUDE CHARACTERISTICS
  ! ** MOLAR VOLUME IS REFERED FROM: PhD THESIS BY NA DU
  !    INVESTIGATION OF HYDROGENATED AND FLUORINATED SURFACTANT
  !    BASED-SYSTEMS FOR THE DESIGN OF POROUS SILICA MATERIALS

  NG  = LA_GRP(NP)
  LNP = LLA(NP)

  ! *** LIMIT EVAPORATION, ONLY WHEN DRIFTER IS NEAR SURFACE
  IF( DLA(NP)<0.1 .AND. DVAP(NG)>0. )THEN
    IF( ICECELL(LNP) )THEN
      VWIND = 0.001
    ELSE
      VWIND  = WINDST(LNP)
      VWIND = MAX(VWIND,0.001)
    ENDIF
    RKA = 0.00125_8*VWIND                ! *** MODIFIED MASS TRANSFER COEFFICENT: RKA = 0.00125_8*VWIND [M/S]
    DTHE = RKA*DARE(NG)*DELTD/DVOL0(NG)  ! *** EVAPORATIVE EXPOSURE FOR A DRIFTER (DIMENSIONLESS)
    IF( ISTRAN(2) > 0 )THEN
      TK = TEM(LNP,KC)+273.15
    ELSE
      TK = DTEM(NG)+273.15
    ENDIF
    HENRY = DVAP(NG)*DVMO(NG)/(RGAS*TK)  ! *** HENRY'S LAW: DIMENSIONLESS, DVMO = 0.000194-0.000293 M3/MOL
    DFV = HENRY*DTHE                     ! *** STIVER AND MACKAY,         DFV: RATIO OF VAPOR VOLUME AND INITIAL VOLUME
    ! *** 1-Fv = (1-Fm)(RHOo/RHOf)   DFV ASSUMES NO CHANGE IN DENSITY (TO BE UPDATED LATER)
    DVOL(NP) = DVOL(NP)*(1.-DFV)
  ENDIF

  ! *** OIL BIODEGRADATION
  IF( DARE(NG)>0. )THEN
    ! *** DETERMINE TEMPERATURE DEPENDENCY
    IF( ISTRAN(2) > 0 )THEN
      TNP = TEM(LNP,KLA(NP))-DTEM(NG)
      TDEP = EXP(-0.01*TNP*TNP)
    ELSE
      TDEP=1.0
    ENDIF
    BIODEG = EXP(-TDEP*DRAT(NG)*DELTD)                     ! *** FRACTION REMAINING AT END OF TIME STEP
    OILMAS = DVOL(NP)*DDEN(NG)*BIODEG                      ! *** MASS OF AN OIL DRIFTER IN KG
    DVOL(NP) = OILMAS/DDEN(NG)                             ! *** DVOL UPDATE
  ENDIF


END SUBROUTINE

SUBROUTINE INIT_OIL
  INTEGER :: NG,NP
  INTEGER :: NUMD(NGRP)
  REAL(RKD)   :: XMAXG,XMING,YMAXG,YMING,RAD
  REAL(RKD)   :: GARE(NGRP)

  ! ** CALCULATING SURFACE AREA OF OIL SPILL FOR EACH GROUP
  GARE = 0  ! TOTAL SURFACE AREA OF OIL FOR EACH GROUP
  DARE = 0  ! SURFACE AREA OF OIL DRIFTER FOR EACH GROUP
  NUMD = 0  ! NUMBER OF OIL DRIFTERS FOR EACH GROUP

  DO NG=1,NGRP
    IF( ISOILSPI(NG) == 0 ) CYCLE

    XMING = MINVAL(XLA,LA_GRP == NG)
    XMAXG = MAXVAL(XLA,LA_GRP == NG)
    YMING = MINVAL(YLA,LA_GRP == NG)
    YMAXG = MAXVAL(YLA,LA_GRP == NG)
    RAD   = MAX(0.25*(XMAXG-XMING+YMAXG-YMING),1.0)
    NUMD(NG) = SUM(LA_GRP,LA_GRP == NG)/NG
    GARE(NG) = PI*RAD**2
    DARE(NG) = GARE(NG)/NUMD(NG)
    DVOL0(NG) = GVOL(NG)/NUMD(NG)
  ENDDO

  ! ** OIL: NPD,NGRP,LA_GRP(NP)
  ALLOCATE(DVOL(NPD))
  DVOL = 0
  DO NP=1,NPD
    IF( ISOILSPI(LA_GRP(NP)) == 0 ) CYCLE
    DVOL(NP) = DVOL0(LA_GRP(NP))           !INITIAL VOLUME OF A DRIFTER

  ENDDO

END SUBROUTINE

!---------------------------------------------------------------------------!
!> @details Checks if a drifter has moved into the ghost cell region of a
!! given subdomain
!
!> @author Zander
!
!---------------------------------------------------------------------------!
Subroutine Check_Drifter_In_Ghost(ID)

  Implicit none

  ! *** Dummy variables
  Integer, intent(in) :: ID

  ! *** Local variables
  Integer :: drifter_local_l !< Cell L index the drifter is currently in
  Integer :: local_i, local_j, index, LG
  Integer :: process_id_to_send

  ! *** Checks if cell id is valid and active
  IF( LLA(ID) > 1 .AND. LLA(ID) < LC .AND. JSPD(ID) == 0 )THEN
    ! *** get the local drifter
    drifter_local_l = LLA(ID)
    
    ! *** Get the local i/j from the L index ofthe current drifter
    local_i = IL(drifter_local_l)  ! Map2Global(drifter_local_l).IL
    local_j = JL(drifter_local_l)  ! Map2Global(drifter_local_l).JL

    ! *** Test if in ghost cells
    IF( local_i < 3 .OR. local_i > IC-2 .OR. local_j < 3 .OR. local_j > jc-2 )THEN
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

      else if( nbr_south == process_id_to_send )then
        ! *** Increment the number of drifters to commmunicate
        num_drifters_send_south = num_drifters_send_south + 1

        drifter_ids_send_south(num_drifters_send_south) = ID
      endif
    ENDIF  
  ENDIF

End Subroutine Check_Drifter_In_Ghost

SUBROUTINE Gather_Drifter_Arrays

  !---------------------------------------------------------------------------!
  ! *** Gathers values for all of Drifter variables needed for EE linkage
  !---------------------------------------------------------------------------!
  Use Mod_Assign_Loc_Glob_For_Write
  Use Mod_Map_Write_EE_Binary

  Implicit None

  ! *** Local variables
  Integer :: ierr, ii, LG, NP, global_id

  ! *** Global arrays before they are sorted
  Integer :: NPD_Write
  Integer,   Save, Allocatable :: map_id(:)
  Integer,   Save, Allocatable :: LLA_Unsorted_Global(:)
  Real(RKD), Save, Allocatable :: XLA_Unsorted_Global(:)
  Real(RKD), Save, Allocatable :: YLA_Unsorted_Global(:)
  Real(RKD), Save, Allocatable :: ZLA_Unsorted_Global(:)
  Logical,   Save, Allocatable :: mask(:)                         !< mask used to exclude deactivated cells or ones that have been communicated to another process

  ! *** Arrays to be allocated and initialized during the "PACK" process
  Integer,   Allocatable :: map_tmp(:)
  Integer,	 Allocatable :: lla_tmp(:)
  Real(RKD), Allocatable :: xla_tmp(:)
  Real(RKD), Allocatable :: yla_tmp(:)
  Real(RKD), Allocatable :: zla_tmp(:)

  ! *** FIRST CALL INITIALIZATION  
  IF( .NOT. ALLOCATED(map_id) )THEN
    ! *** Allocate Arrays
    Allocate(gl_drift_mapping(NPD))

    Allocate(map_id(NPD))
    Allocate(LLA_Unsorted_Global(NPD))
    Allocate(XLA_Unsorted_Global(NPD))
    Allocate(YLA_Unsorted_Global(NPD))
    Allocate(ZLA_Unsorted_Global(NPD))
    Allocate(mask(NPD))
    
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
  ENDIF
  
  ! *** Convert local LLA to have global L values
  Do NP = 1,NPD
    map_id(NP) = 0
    gl_drift_mapping(NP) = 0
    
    ! *** Check if this drifter belongs to the current subdomain
    if( LLA_Process(NP) /= process_id )then
      mask(NP) = .FALSE.
      cycle
    end if
  
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
  End do

  ! *** This removes parts of each local array that does not have any drifter info in it.  Makes a coniguous array to make the gathering easier
  lla_tmp = PACK(LLA_Global, mask)
  xla_tmp = PACK(xla, mask)
  yla_tmp = PACK(yla, mask)
  zla_tmp = PACK(zla, mask)
  map_tmp = PACK(map_id, mask)

  ! *** Gather up LLA, use generic Gather_Soln routine so we can get the mapping
  Call Gather_Drifter(size(lla_tmp), lla_tmp, NPD, LLA_Unsorted_Global)
  
  ! *** Need to get the mapping so we can put drifters in the correct global position
  Call Gather_Drifter(size(map_tmp), map_tmp, NPD, gl_drift_mapping)
  
  ! *** Gather just the XLA array
  Call Gather_Drifter(size(xla_tmp), xla_tmp, NPD, XLA_Unsorted_Global)
  
  ! *** Gather just the YLA array
  Call Gather_Drifter(size(yla_tmp), yla_tmp, NPD, YLA_Unsorted_Global)
  
  ! *** Gather just the ZLA array
  Call Gather_Drifter(size(zla_tmp), zla_tmp, NPD, ZLA_Unsorted_Global)

  if( process_id == master_id )then
    ! *** Sort so the global indexing matches when writing out
    DO NP = 1,NPD

      ! *** Get the correct global drifter id
      global_id = gl_drift_mapping(NP)

      IF( global_id > 0 .and. global_id <= NPD )THEN
        LLA_Global(global_id) =  LLA_Unsorted_Global(NP)
        
        XLA(global_id) =  XLA_Unsorted_Global(NP)
        YLA(global_id) =  YLA_Unsorted_Global(NP)
        ZLA(global_id) =  ZLA_Unsorted_Global(NP)
      !else
      !  LLA(NP) = 1
      !  XLA(NP) = 0.0
      !  YLA(NP) = 0.0
      !  ZLA(NP) = 0.0
      END IF
    END DO
  endif

  ! *** Deallocate to save memory as this is used infrequently
  Deallocate(lla_tmp)
  Deallocate(xla_tmp)
  Deallocate(yla_tmp)
  Deallocate(zla_tmp)
  Deallocate(map_tmp)
  
END SUBROUTINE Gather_Drifter_Arrays

! *** Binary search of a sorted list
INTEGER FUNCTION Search_List(NN, NLIST, iTarget)

  INTEGER, INTENT(IN) :: NN, iTarget
  INTEGER, INTENT(IN) :: NLIST(NPD)
  
  ! *** Local Variables
  Integer :: I1, I2, I3, nnn
  
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

  INTEGER, INTENT(IN)    :: iTarget
  INTEGER, INTENT(INOUT) :: NN, NLIST(NPD)
  
  ! *** Local Variables
  Integer :: I1, I2, I3
  
  IF( NN == NPD ) RETURN
  
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
  DO I3 = NN,I1,-1
    NLIST(I3+1) = NLIST(I3)
  ENDDO
  NN = NN + 1
  NLIST(I2) = iTarget
  
End SUBROUTINE Insert_Value
#ifdef NCOUT
    subroutine create_nc_lpt()
        real(4), parameter :: MISSING_VALUE=-999.
        integer(4), parameter :: TIMECNT=NF90_UNLIMITED
        character(80) :: filename
        character(20) :: bdate
        integer :: status

        if( ISNCDF(9) == 0 ) return

        lpt_time_idx = 0
    
        BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)
        filename  = OUTDIR//'efdc_drifters.nc'
        
            ! ** ENTER DEFINE MODE
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
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'grid_mapping','crs')
        !status = nf90_put_att(nc_lpt(0),nc_lpt(4),'coordinates','trk_lat trk_lon')
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
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'coordinates','trk_time trk_lon trk_lat trk_lev')
            status = nf90_put_att(nc_lpt(0),nc_lpt(5),'_FillValue',MISSING_VALUE)

            status = nf90_def_var(nc_lpt(0),'oil_mass',nf90_float,(/lpt_npd_dim,lpt_time_dim/), nc_lpt(6))
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'standard_name','')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'long_name','drifter oil mass')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'units','mg/L')
            !status = nf90_put_att(nc_lpt(0),nc_lpt(6),'grid_mapping','crs')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'coordinates','trk_time trk_lon trk_lat trk_lev')
            status = nf90_put_att(nc_lpt(0),nc_lpt(6),'_FillValue',MISSING_VALUE)
        endif
        
        ! ** ASSIGN GLOBAL ATTRIBUTES
        status = nf90_put_att(nc_lpt(0),nf90_global,'Conventions','CF-1.4')
        status = nf90_put_att(nc_lpt(0),NF90_GLOBAL,'Base_date',BDATE)

        ! ** LEAVE DEFINE MODE
        status = nf90_enddef(nc_lpt(0))
    end subroutine 
    
    subroutine close_nc_lpt()
        integer :: status
        if( ISNCDF(9) == 0 ) return
        status = nf90_close(nc_lpt(0))
        if(status /= NF90_NOERR )then
            print * ,'cannot close efdc_drifters.nc!'
            return
        endif
    end subroutine
    
    subroutine write_nc_lpt()
        real(4), parameter :: MISSING_VALUE=-999.
        real(8) :: ti(1), xm(1),ym(1),xll(1),yll(1)
        real(4), allocatable :: array2d(:,:)
        integer :: i,status,str1d(2),cnt1d(2)

        if( ISNCDF(9) == 0 ) return
        allocate(array2d(2,NPD))
        array2d = MISSING_VALUE
        do i=1,NPD
            if( LLA_Global(i) < 2 .OR. LLA_Global(i) > LA_Global) cycle ! *** Modified for global
            xm(1)=XLA(i)
            ym(1)=YLA(i)
            call utmr_wgs84(xm,ym,xll,yll)
            array2d(1,i) = xll(1)
            array2d(2,i) = yll(1)
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
#endif
END MODULE
