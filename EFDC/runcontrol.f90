! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
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

  use IFCORE
  
  KEY_PRESSED = PEEKCHARQQ ( )
  
END FUNCTION KEY_PRESSED

! *****************************************************************************
LOGICAL FUNCTION ISEXIT()
  ! *** DETERMINES IF THE KEYBOARD HIT RESUMES (I1 /= I2) OR TERMINATES (I1 == I2) THE RUN
  
  use IFCORE
  use GLOBAL,only:IK4
  
  integer(IK4) :: I1, I2
  character*1 :: KEY
  
  KEY = GETCHARQQ()
  I1 = ICHAR(KEY)
  
  write(*,'(A)')'PROGRAM PAUSED BY USER'  
  write(*,'(A)')'  EFDC+: TO EXIT PRESS THE SAME KEY: ' // KEY
  write(*,'(A)')'  EFDC+: TO CONTINUE RUN PRESS ANY OTHER KEY'  

  KEY = GETCHARQQ()
  I2 = ICHAR(KEY)
  
  if( I1 /= I2 )then
    ISEXIT = .FALSE.
  else
    ISEXIT = .TRUE.
  endif
  
END FUNCTION ISEXIT

! *****************************************************************************
SUBROUTINE QUIT

  use IFCORE
  character KEY*1  

  write(6,'(/,A28,$)')'TAP ANY KEY TO EXIT EFDC_DSI'

  KEY = GETCHARQQ()

  return
  
  END SUBROUTINE
#endif

