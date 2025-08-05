! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALVEGSER()

  ! *** SUBROUTINE CALVEGSR UPDATES TIME VARIABLE VEGETATION RESISTANCE
  ! *** parameterS

  !  NVEGSER = NUMBER OF VEGETATION TIME SERIES
  !  NVEGSERV(NVEGTPM) = TIME SERIES ID FOR SPECIFIC VEGETATION CLASS
  !  MVEGTLAST(NVEGSERM) = PLACE HOLDER IN INTERPOLATION TABLE
  !  TCVEGSER(NVEGSERM) = TIME CONVERSION FACTOR FOR TIME VARIABLE
  !  TVEGSER(NDVEGSER,NVEGSERM) = TIME OF DATA
  !  VEGSERRT(NVEGSERM) = CURRENT VALUE OF RDLPSQ
  !  VEGSERBT(NVEGSERM) = CURRENT VALUE OF BPVEG
  !  VEGSERHT(NVEGSERM) = CURRENT VALUE OF HPVEG
  !  VEGSERR(NDVEGSER,NVEGSERM) = TIME VARYING VALUES OF RDLPSQ
  !  VEGSERB(NDVEGSER,NVEGSERM) = TIME VARYING VALUES OF BPVEG
  !  VEGSERH(NDVEGSER,NVEGSERM) = TIME VARYING VALUES OF HPVEG

  ! CHANGE RECORD

  use GLOBAL
  
  implicit none
  
  integer :: NS,M1,M2,M,NSTMP
  real   :: TIME,TDIFF,WTM1,WTM2,BDLTMP

  if( NVEGSER > 0 )then
    do NS = 1,NVEGSER
      TIME = TIMESEC/TCVEGSER(NS)  

      M1 = MVEGTLAST(NS)
    100     continue
      M2 = M1+1
      if( TIME > TVEGSER(M2,NS) )then
        M1 = M2
        GOTO 100
      else
        MVEGTLAST(NS) = M1
      endif
      TDIFF = TVEGSER(M2,NS)-TVEGSER(M1,NS)
      WTM1 = (TVEGSER(M2,NS)-TIME)/TDIFF
      WTM2 = (TIME-TVEGSER(M1,NS))/TDIFF
      VEGSERRT(NS) = WTM1*VEGSERR(M1,NS)+WTM2*VEGSERR(M2,NS)
      VEGSERBT(NS) = WTM1*VEGSERB(M1,NS)+WTM2*VEGSERB(M2,NS)
      VEGSERHT(NS) = WTM1*VEGSERH(M1,NS)+WTM2*VEGSERH(M2,NS)
    enddo
    do M = 1,MVEGTYP
      NSTMP = NVEGSERV(M)
      if( NSTMP > 0 )then
        RDLPSQ(M) = VEGSERRT(NSTMP)
        BPVEG(M) = VEGSERBT(NSTMP)
        HPVEG(M) = VEGSERHT(NSTMP)
        BDLTMP = BPVEG(M)*BPVEG(M)*RDLPSQ(M)    ! *** Population density (dimensionless)
        PVEGZ(M) = 1.-ALPVEG(M)*BDLTMP
        BDLPSQ(M) = BPVEG(M)*RDLPSQ(M)
      endif
    enddo
  endif
  return
END

