! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE JULIANMOD
  ! *** AUTHOR: DANG HUU CHUNG
  ! *** DATE  : 2010-06-26
  use GLOBAL  

  implicit none
   ! *** JULIAN TIME
    integer(IK4),PRIVATE  :: DA0             !BASE DAY
    integer(IK4),PRIVATE  :: MN0             !BASE MONTH
    integer(IK4),PRIVATE  :: YR0             !BASE YEAR
    integer(IK4),PRIVATE  :: HH0             !BASE HH
    integer(IK4),PRIVATE  :: MM0             !BASE MM
    integer(IK4),PRIVATE  :: SS0             !BASE SS 
    real(8),   PRIVATE  :: HR0
  
  contains

  FUNCTION JULIDAY(MM,ID,IYYY) RESULT(JULDAY)
    ! *** FROM A BOOK OF NUMERICAL METHODS
    integer(IK4) :: JULDAY,ID,IYYY,MM,IGREG
    parameter (IGREG = 15+31*(10+12*1582))
    integer(IK4) :: JA,JM,JY
    JY = IYYY
    if( JY == 0) STOP 'JULDAY: THERE IS NO YEAR ZERO'
    if( JY < 0) JY = JY+1
    if( MM > 2 )then
      JM = MM+1
    else
      JY = JY-1
      JM = MM+13
    endif
    JULDAY = 365*JY+INT(0.25D0*JY+2000D0)+INT(30.6001D0*JM)+ID+1718995
    if( ID+31*(MM+12*IYYY).GE.IGREG )then
      JA = INT(0.01D0*JY)
      JULDAY = JULDAY+2-JA+INT(0.25D0*JA)
    endif
    return
  END FUNCTION

  SUBROUTINE CALDAT(JULIAN,MM,ID,IYYY)
    ! *** FROM A BOOK OF NUMERICAL METHODS
    integer(IK4) :: ID,IYYY,JULIAN,MM,IGREG
    parameter (IGREG = 2299161)
    integer(IK4) :: JA,JALPHA,JB,JC,JD,JE
    if( JULIAN.GE.IGREG )then
      JALPHA = INT(((JULIAN-1867216)-0.25D0)/36524.25D0)
      JA = JULIAN+1+JALPHA-INT(0.25D0*JALPHA)
    elseif( JULIAN < 0 )then
      JA = JULIAN+36525*(1-JULIAN/36525)
    else
      JA = JULIAN
    endif
    JB = JA+1524
    JC = INT(6680.0D0+((JB-2439870)-122.1D0)/365.25D0)
    JD = 365*JC+INT(0.25D0*JC)
    JE = INT((JB-JD)/30.6001D0)
    ID = JB-JD-INT(30.6001D0*JE)
    MM = JE-1
    if( MM > 12 )MM = MM-12
    IYYY = JC-4715
    if( MM > 2 )IYYY = IYYY-1
    if( IYYY <= 0)IYYY = IYYY-1
    if( JULIAN < 0 ) IYYY = IYYY-100*(1-JULIAN/36525)
    return
  END SUBROUTINE

   SUBROUTINE DATEPRO(STR,YR,MN,DD)
   
     character(*),intent(IN) :: STR
     character(LEN = LEN_TRIM(STR)) :: SS
     integer(IK4),intent(OUT) :: DD,MN,YR
     integer(IK4) :: M,NL
     
     SS = STR
     NL = LEN_TRIM(SS)
     do M = 1,NL
       if( SS(M:M) == '-' .or. SS(M:M) == ':' .or. SS(M:M) == '/' )then
         SS(M:M)  = ''
       endif
     enddo
     read(SS,*) YR,MN,DD
     
   END SUBROUTINE
 
  SUBROUTINE TOGREGOR(SDATETIME,DAYTIME)
  
    real(8),intent(IN) :: DAYTIME
    character(*),intent(OUT) :: SDATETIME
    integer(IK4) :: JD0,JULIAN,MN,DD,YR 

    call DATEPRO(BASETIME,HH0,MM0,SS0) 
    call DATEPRO(BASEDATE,YR0,MN0,DA0) 

    HR0 = HH0+MM0/60._8+SS0/3600._8
    JD0 = JULIDAY(MN0,DA0,YR0)         
   
    JULIAN = IDINT(DAYTIME+JD0)
    call CALDAT(JULIAN,MN,DD,YR)

    write(SDATETIME,'(I4,I2.2,I2.2)') YR,MN,DD
  
  END SUBROUTINE
 
END MODULE
