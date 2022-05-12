! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
REAL*8 FUNCTION DSTIME(IOPTION)  
  ! *** Generic Function to Provide Model Timing
  ! *** DSTIME returns back the number of seconds since some event
  USE GLOBAL,ONLY:NTHREADS,RKD,IK4
  USE IFPORT
  USE MPI
  
  INTEGER(IK4), INTENT(IN) :: IOPTION
  
  REAL                :: TARR(2), TPMC
  REAL,STATIC         :: LASTTIME
  REAL(RKD)           :: CPUSEC
  
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:) :: STARTEDD  

  IF( .NOT. ALLOCATED(STARTEDD) )THEN
    ALLOCATE(STARTEDD(0:3))
    STARTEDD = 0.
    
    ! *** ELAPSED TIME FROM January 1, 1970
    STARTEDD(1) = RTC()
  ENDIF
    
  ! *** CHANGE THE FOLLOWING LINES TO CORRESPOND TO THE PLATORM AND COMPILER
  IF( IOPTION == 0 )THEN
    ! *** MPI Timing in seconds
    DSTIME = MPI_WTIME()
    !! *** TOTAL RUN CPU TIME
    !CALL CPU_TIME(TPMC)
    !DSTIME = DBLE(TPMC)/DBLE(NTHREADS)
    
  ELSEIF( IOPTION == 1 )THEN
    ! *** ELAPSED TIME FROM January 1, 1970
    CPUSEC=RTC()
    DSTIME=CPUSEC-STARTEDD(1)
    
  ELSEIF( IOPTION == 2 )THEN
    ! *** TOTAL RUN CPU TIME
    TPMC=DTIME(TARR)
    DSTIME=DBLE(TPMC)/DBLE(NTHREADS)
    
    IF( DSTIME < 0.0 )THEN
      DSTIME=ABS(DSTIME)
    ENDIF
    IF( DSTIME < LASTTIME)DSTIME = LASTTIME
    LASTTIME = DSTIME
    
  ELSEIF( IOPTION == 3 )THEN
    ! *** RUN SPECIFIC CPU TIME ("USER TIME")
    TPMC=DTIME(TARR)
    
    DSTIME=DBLE(TARR(1))/DBLE(NTHREADS)
    IF( DSTIME < LASTTIME)DSTIME = LASTTIME
    LASTTIME = DSTIME
    
  ELSEIF( IOPTION == 4 )THEN
    ! *** MPI Timing in seconds
    DSTIME = MPI_WTIME()
  ENDIF

  RETURN

END FUNCTION
