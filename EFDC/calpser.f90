! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALPSER

  ! CHANGE RECORD
  ! *** SUBROUTINE CALPSER UPDATES TIME VARIABLE SURFACE ELEVATION
  ! *** BOUNDARY CONDITIONS

  use GLOBAL
  implicit none
  
  integer   :: NS,M2,M1
  real(RKD) :: TIME,TDIFF
  real      :: WTM1,WTM2
  real      :: Y1,Y2,Y3,Y4
  
  PSERT(0) = 0.
  
  do NS = 1,NPSER
    TIME = TIMESEC/DBLE(TSPS(NS).TMULT)

    M2 = MTSPLAST(NS)
    do while (TIME > TSPS(NS).TIM(M2) )
      M2 = M2+1
      if( M2 > TSPS(NS).NREC )then
        M2 = TSPS(NS).NREC
        EXIT
      endif
    enddo
    MTSPLAST(NS) = M2  
    M1 = M2-1
    
    TDIFF = TSPS(NS).TIM(M2) - TSPS(NS).TIM(M1)
    
    if( INTPSER(NS) == 0 .or. M1 < 2 .or. M2 > TSPS(NS).NREC-1 )then
      ! *** LINEAR INTERPOLATION
      WTM1 = (TSPS(NS).TIM(M2) - TIME)/TDIFF 
      WTM2 = (TIME - TSPS(NS).TIM(M1))/TDIFF
      PSERT(NS)  = WTM1*TSPS(NS).VAL(M1,1) + WTM2*TSPS(NS).VAL(M2,1) + PDGINIT  ! *** ADD OFFSET    (m2/s2)
      PSERST(NS) = WTM1*TSPS(NS).VAL(M1,2) + WTM2*TSPS(NS).VAL(M2,2) + PDGINIT  ! *** ADD OFFSET    (m2/s2)
    else
      ! *** CATMULL–ROM SPLINE
      WTM1 = (TIME - TSPS(NS).TIM(M1))/TDIFF
      Y1 = TSPS(NS).VAL(M1-1,1)
      Y2 = TSPS(NS).VAL(M1,1)
      Y3 = TSPS(NS).VAL(M2,1)
      Y4 = TSPS(NS).VAL(M2+1,1)
      PSERT(NS) = 0.5*( (2.*Y2) + (-Y1 + Y3)*WTM1 + (2.*Y1 - 5.*Y2 + 4.*Y3 - Y4)*WTM1**2 + (-Y1 + 3.*Y2 - 3.*Y3 + Y4)*WTM1**3 )
      
      ! *** Uncomment out the following for QC tests of the spline interpolation
      !if(ns==3) WRITE(100+NS,'(F12.6,F10.3,2I10,F8.5,2F10.3)') TIME,PSERT(NS)/g, M1,M2,WTM1, TSPS(NS).VAL(M1,1)/g, TSPS(NS).VAL(M2,1)/g                        ! DELME
      !IF( MOD(NITER,1000)==0 .and. NS==3 ) WRITE(*,'(F12.6,F10.3,2I10,F8.5,2F10.3)') TIME,PSERT(NS)/g, M1,M2,WTM1, TSPS(NS).VAL(M1,1)/g, TSPS(NS).VAL(M2,1)/g  ! DELME
      
      Y1 = TSPS(NS).VAL(M1-1,2)
      Y2 = TSPS(NS).VAL(M1,2)
      Y3 = TSPS(NS).VAL(M2,2)
      Y4 = TSPS(NS).VAL(M2+1,2)
      PSERST(NS) = 0.5*( (2.*Y2) + (-Y1 + Y3)*WTM1 + (2.*Y1 - 5.*Y2 + 4.*Y3 - Y4)*WTM1**2 + (-Y1 + 3.*Y2 - 3.*Y3 + Y4)*WTM1**3 )

    endif
  enddo
  
  return
  
END

