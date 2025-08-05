! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
REAL*8 FUNCTION DSTIME(IOPTION)  
  ! *** Generic Function to Provide Model Timing
  ! *** DSTIME returns back the number of seconds since some event
  use GLOBAL,only:NTHREADS,RKD,IK4
#ifndef GNU  
  use IFPORT
#endif
  use OMP_LIB
  use MPI
  
  integer(IK4), intent(IN) :: IOPTION
  
  real                :: TARR(2), TPMC
  real,STATIC         :: LASTTIME
  real(RKD)           :: CPUSEC
  
  real(RKD),save,allocatable,dimension(:) :: STARTEDD  

  DSTIME = 0.0
  
  if( .not. allocated(STARTEDD) )then
    allocate(STARTEDD(0:3))
    STARTEDD = 0.
    
#ifdef GNU
    call CPU_TIME(TPMC)
    STARTEDD(1) = TPMC
#else
    ! *** ELAPSED TIME FROM January 1, 1970
    STARTEDD(1) = RTC()
#endif     
  endif
  
#ifdef GNU  
    ! *** TOTAL RUN CPU TIME
    call CPU_TIME(TPMC)
    DSTIME = DBLE(TPMC)-STARTEDD(1)
    !DSTIME = DBLE(TPMC)/DBLE(NTHREADS)
#else
  ! *** CHANGE THE FOLLOWING LINES TO CORRESPOND TO THE PLATORM AND COMPILER
  if( IOPTION == 0 )then
    ! *** MPI Timing in seconds
#ifdef DEBUGGING
    DSTIME = MPI_WTIME()
#else
    DSTIME = omp_get_wtime()
#endif
    !! *** TOTAL RUN CPU TIME
    !CALL CPU_TIME(TPMC)
    !DSTIME = DBLE(TPMC)/DBLE(NTHREADS)
    
  elseif( IOPTION == 1 )then
    ! *** ELAPSED TIME FROM January 1, 1970
    CPUSEC = RTC()
    DSTIME = CPUSEC - STARTEDD(1)
    
  elseif( IOPTION == 2 )then
    ! *** TOTAL RUN CPU TIME
    TPMC = DTIME(TARR)
    DSTIME = DBLE(TPMC)/DBLE(NTHREADS)
    
    if( DSTIME < 0.0 )then
      DSTIME = ABS(DSTIME)
    endif
    if( DSTIME < LASTTIME)DSTIME = LASTTIME
    LASTTIME = DSTIME
    
  elseif( IOPTION == 3 )then
    ! *** RUN SPECIFIC CPU TIME ("USER TIME")
    TPMC = DTIME(TARR)
    
    DSTIME = DBLE(TARR(1))/DBLE(NTHREADS)
    if( DSTIME < LASTTIME ) DSTIME = LASTTIME
    LASTTIME = DSTIME
    
  elseif( IOPTION == 4 )then
    ! *** MPI Timing in seconds
#ifdef DEBUGGING
    DSTIME = MPI_WTIME()
#else
    DSTIME = omp_get_wtime()
#endif
  endif
#endif 
  return

END FUNCTION
