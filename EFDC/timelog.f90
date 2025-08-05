! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE TIMELOG(N,DAYJUL,OUTDIR1,TTDS) 
   
  ! *** WRITE OUT MODEL TIME STEP AND SUN/PC SYSTEM CLOCK TIME TO TIME.LOG FILE

  use GLOBAL,only:IK8,NITER
  
  implicit none
  
  integer(IK8) :: N
  DOUBLE PRECISION, intent(IN) :: DAYJUL, TTDS
  character*8,      intent(IN) :: OUTDIR1
  character*8 MRMDATE,MRMTIME*10  
  
  ! *** WRITE OUT MODEL TIME STEP AND SYSTEM CLOCK TIME TO TIME.LOG  
  call DATE_AND_TIME(MRMDATE,MRMTIME)  
    
  open(9,FILE = OUTDIR1//'TIME.LOG',POSITION = 'APPEND')  
  write(9,100) NITER, DAYJUL, MRMDATE, MRMTIME, TTDS/3600. 
  close(9)

  100 FORMAT(' ','NITER  = ',I12,5X,'TIMEDAY = ',F12.4,5X,'DATE = ',A8,5X,'TIME = ',A10,5X,'Elapsed Time (hrs) = ',F12.5)  
  
END SUBROUTINE
 
