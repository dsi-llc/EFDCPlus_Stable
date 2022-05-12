! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! *** THE FUNCTIONS IN THIS FILE WERE ADDED TO REPLACE THE C LIB FUNCTIONS OF KBHIT AND ATEXIT
! *** THE NEW APPROACH FOR KEYBOARD RUN CONTROLS ALLOWS FOR 32BIT AND/OR 64BIT COMPILATIONS
!
! CHANGE RECORD  
! DATE MODIFIED     BY               DESCRIPTION
!----------------------------------------------------------------------!
! 2011-03           PAUL M. CRAIG    ADDED NEW APPROACH FOR KEYBOARD RUN CONTROLS
!                                    TO ELIMINATE INCOMPATIBILITIES WITH 64BIT  

#ifdef _WIN  
! *****************************************************************************
LOGICAL FUNCTION KEY_PRESSED()  
  ! *** DETERMINES IF THE USER HAS PRESSED A KEY

  USE IFCORE
  
  KEY_PRESSED = PEEKCHARQQ ( )
  
END FUNCTION KEY_PRESSED

! *****************************************************************************
LOGICAL FUNCTION ISEXIT()
  ! *** DETERMINES IF THE KEYBOARD HIT RESUMES (I1 /= I2) OR TERMINATES (I1 == I2) THE RUN
  
  USE IFCORE
  USE GLOBAL,ONLY:IK4
  
  INTEGER(IK4) :: I1, I2
  CHARACTER*1 :: KEY
  
  KEY = GETCHARQQ()
  I1=ICHAR(KEY)
  
  WRITE(*,'(A)')'PROGRAM PAUSED BY USER'  
  WRITE(*,'(A)')'  EFDC_DSI: TO EXIT PRESS THE SAME KEY: ' // KEY
  WRITE(*,'(A)')'  EFDC_DSI: TO CONTINUE RUN PRESS ANY OTHER KEY'  

  KEY = GETCHARQQ()
  I2=ICHAR(KEY)
  
  IF( I1 /= I2 )THEN
    ISEXIT = .FALSE.
  ELSE
    ISEXIT = .TRUE.
  ENDIF
  
END FUNCTION ISEXIT

! *****************************************************************************
SUBROUTINE QUIT

  USE IFCORE
  CHARACTER KEY*1  

  WRITE(6,'(/,A28,$)')'TAP ANY KEY TO EXIT EFDC_DSI'

  KEY = GETCHARQQ()

  RETURN
  
  END SUBROUTINE
#endif

