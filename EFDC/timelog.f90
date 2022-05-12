! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE TIMELOG(N,DAYJUL,OUTDIR1,TTDS) 
   
  ! *** WRITE OUT MODEL TIME STEP AND SUN/PC SYSTEM CLOCK TIME TO TIME.LOG FILE

  USE GLOBAL,ONLY:IK8,NITER
  
  IMPLICIT NONE
  
  INTEGER(IK8) :: N
  DOUBLE PRECISION, INTENT(IN) :: DAYJUL, TTDS
  CHARACTER*8,      INTENT(IN) :: OUTDIR1
  CHARACTER*8 MRMDATE,MRMTIME*10  
  
  ! *** WRITE OUT MODEL TIME STEP AND SYSTEM CLOCK TIME TO TIME.LOG  
  CALL DATE_AND_TIME(MRMDATE,MRMTIME)  
    
  OPEN(9,FILE=OUTDIR1//'TIME.LOG',POSITION='APPEND')  
  WRITE(9,100) NITER, DAYJUL, MRMDATE, MRMTIME, TTDS/3600. 
  CLOSE(9)

  100 FORMAT(' ','NITER =',I12,5X,'TIMEDAY = ',F12.4,5X,'DATE = ',A8,5X,'TIME = ',A10,5X,'Elapsed Time (hrs) = ',F12.5)  
  
END SUBROUTINE
 
