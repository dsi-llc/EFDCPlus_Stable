! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE INFOMOD
!Author: Dang Huu Chung

use GLOBAL,only:IK4
implicit none
contains

FUNCTION READSTR(UNIT) RESULT(STR)
  integer(IK4),intent(IN) :: UNIT
  character(200) :: STR
  integer(IK4) :: ISTR,I
  do while (.TRUE.)
    read(UNIT,'(A)',err = 1000,end = 1010) STR
    STR = ADJUSTL(STR)
    ISTR = ICHAR(STR(1:1))
    I = 1
    do while (ISTR == 9) 
      I = I+1
      ISTR = ICHAR(STR(I:I))
    enddo
    SELECT CASE (ISTR)
    CASE (45,46,48:57)  !CHARACTER = -, ., 0:9
      BACKSPACE UNIT
      return
    END SELECT
  enddo
  return
  
1000 call STOPP('READ ERROR!')
1010 call STOPP('END OF FILE BEFORE EXPECTED!')
 
END FUNCTION

SUBROUTINE SKIPCOM(IUNIT,CC,IUOUT)
  integer(IK4),  intent(IN) :: IUNIT    
  integer(IK4),  intent(IN),OPTIONAL :: IUOUT
  character(1),intent(IN) :: CC
  character(250) :: LINE,COMM*1(4)
  integer(IK4)   :: I,ISTR
  DATA COMM /'C','c','*','#'/
 
  do while(.TRUE.)
    read(IUNIT, '(A)', end = 999) LINE      
    if( PRESENT(IUOUT)) WRITE(IUOUT,'(A)') LINE
    LINE = ADJUSTL(LINE)
    ISTR = ICHAR(LINE(1:1))
    I = 1
    do while (ISTR == 9) 
      I = I+1
      ISTR = ICHAR(LINE(I:I))
    enddo
    if( LINE(I:I) == CC .or. ANY(COMM == LINE(I:I)) )then
      CYCLE
    else
      BACKSPACE(IUNIT)
      exit
    endif
  enddo
  999 RETURN
END SUBROUTINE
 
FUNCTION FINDSTR(STR,SS,NCOL) RESULT(COLM)
  character(*)  :: STR,SS
  character(10) :: SSN(NCOL)
  integer(IK4)  :: M,COLM,NCOL,NL
  COLM = 0
  read(STR,*,end = 100) (SSN(M),M = 1,NCOL)
  100 continue
  do M = 1,NCOL
    SSN(M) = ADJUSTL(SSN(M))
    NL = INDEX(SSN(M),SS)
    if( NL > 0 )then
      COLM = M
      return
    endif
  enddo
END FUNCTION

FUNCTION NUMCOL(STR) RESULT(NC)

  integer(IK4) :: M,NC,NL
  character(*) :: STR,STR1*200

  STR1 = ADJUSTL(STR)
  NL = LEN_TRIM(STR1)
  if( NL == 0 )then
    NC = 0
    return
  endif
  NC = 1
  do M = 2,NL
    if( STR1(M:M) == '' .and. STR1(M-1:M-1)/='' )then
      NC = NC+1
    endif
  enddo
END FUNCTION
 
END MODULE
 

