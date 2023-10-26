! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEDZLJ_SHEAR
  
  ! *** CALCULATES Wave and Current Shear Stress Based on Log-Law of Cristofferson Jonsson
  ! *** 
  ! *** ORIGINAL:  May 24, 2007
  ! ***  Craig Jones and Scott James
  ! *** REVISED: SIGMA-ZED AND OMP - 2016-11-07
  ! ***  Paul M. Craig

  USE GLOBAL
  USE WINDWAVE
  Use Allocate_Initialize      
  
  IMPLICIT NONE
  
  INTEGER :: L,ND,LF,LL,LP
  INTEGER :: M1,M2
  INTEGER :: FZONE
  
  REAL(RKD)  :: MMW, SIGMAWV, JJW, KN, FWVTPP, ABM, WA, WR, UWBM, WAMP
  REAL(RKD)  :: VELMAG, VELANG, DELW, APROUGH, KNFACTOR
  REAL(RKD)  :: UTMP, VTMP
  REAL(RKD)  :: WVLENGTH, WVANGLE, WFTIM, EXCURSION
  REAL(RKD)  :: FC1, FC2, FWINDSQ, FC, FWW, FWVHT, SHEAR, SHEARC, SHEARW, GROWTH
  REAL(RKD)  :: FWINDS, FWINDD
  REAL(RKD)  :: TDIFF, WTM1, WTM2, AVGDEPTH
  
  REAL(RKD) ,SAVE,ALLOCATABLE,DIMENSION(:) :: WVFREQ
  REAL(RKD) ,SAVE,ALLOCATABLE,DIMENSION(:) :: WVORBIT
  REAL(RKD) ,SAVE,ALLOCATABLE,DIMENSION(:) :: ZBTEMP

  IF( .NOT. ALLOCATED(WVFREQ) )THEN
    Call AllocateDSI(WVFREQ,  LCM, 0.0)
    Call AllocateDSI(WVORBIT, LCM, 0.0)
    Call AllocateDSI(ZBTEMP,  LCM, 0.0)
  ENDIF
  GROWTH = 0.1_8
  
  ! ***************************************************************************

  IF( TAUCONST == 0 )THEN
    ! *** Constant Tau not used.  Compute spatially variable Tau

    ! ***************************************************************************
    IF( ISWAVE > 0 .AND. ISWNWAVE == 0 )THEN
      ! *** Using EFDC+ windwave parameters
      CALL WINDWAVECAL    ! *** Compute all wave parameters from the fetch and wind fields
      
      ! *** Assign to SEDZLJ variables
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND, LF, LL, LP, L)   &
      !$OMP             SHARED(NDM, LDMSED, LASED, LSED, WV, WVFREQ, WVORBIT, UWVSQ, WVANG)
      DO ND = 1,NDM  
        LF = (ND-1)*LDMSED+1  
        LL = MIN(LF+LDMSED-1,LASED)
        DO LP = LF,LL
          L = LSED(LP)
          WVFREQ(L)  = WV(L).FREQ
          WVORBIT(L) = WV(L).UDEL
          UWVSQ(L)   = WV(L).UDEL*WV(L).UDEL
          WVANG(L)   = WV(L).DIR          
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      
    ELSEIF( ISWAVE > 0 .AND. ISWNWAVE == 1 )THEN
      ! Wind Wave Fetch
      
      ! Convert wind input into current wind info for wind-driven wave calcs
      WFTIM = TIMESEC/TCWSER(1)    ! *** 3tl and 2tl both use TIMESEC
      M2 = MTSWLAST(1)
      DO WHILE (WFTIM > TSWND(1).TIM(M2))
        M2 = M2+1
        IF( M2 > TSWND(1).NREC )THEN
          M2 = TSWND(1).NREC
          EXIT
        ENDIF
      END DO
      MTSWLAST(1) = M2 
      M1 = M2-1
      TDIFF = TSWND(1).TIM(M2)-TSWND(1).TIM(M1)  
      WTM1 = (TSWND(1).TIM(M2)-WFTIM)/TDIFF  
      WTM2 = (WFTIM-TSWND(1).TIM(M1))/TDIFF 
      FWINDS = WTM1*TSWND(1).VAL(M1,1)+WTM2*TSWND(1).VAL(M2,1)

      IF( FWINDS > 1.0 )THEN
        ! *** Only compute wind wave if wind speed > 1 m/s 
        IF( ABS(TSWND(1).VAL(M1,2)-TSWND(1).VAL(M2,2))<180.0 )THEN
          FWINDD = WTM1*TSWND(1).VAL(M1,2)+WTM2*TSWND(1).VAL(M2,2)
        ELSE
          IF( TSWND(1).VAL(M1,2) > TSWND(1).VAL(M2,2) )THEN
              FWINDD = WTM1*TSWND(1).VAL(M1,2)+WTM2*(TSWND(1).VAL(M2,2)+360.0)
          ELSE
              FWINDD = WTM1*(TSWND(1).VAL(M1,2)+360)+WTM2*TSWND(1).VAL(M2,2)
          ENDIF
          IF( FWINDD >= 360.0)FWINDD = FWINDD-360.0 
        ENDIF
        ! Convert wind into direction it is blowing "from"
        IF( FWINDD <= 180.0 )THEN  
          FWINDD = FWINDD+180.0  
          IF( FWINDD == 360.0)FWINDD = 0.0         
        ELSE  
          FWINDD = FWINDD-180.0  
          IF( FWINDD == 360.0)FWINDD = 0.0 
        ENDIF
        
        ! *** Calculate which of the 8 wind zones (FWZONE) the wind is coming from
        ! *** Also the Waveangle CCW from East.  Waveangle is center of sector.
        IF( FWINDD >= 337.5 .OR. FWINDD<22.5 )THEN
          FZONE = 1
          WVANGLE = 4.712
        ELSEIF( FWINDD >= 22.5 .AND. FWINDD<67.5 )THEN
          FZONE = 2
          WVANGLE = 3.927
        ELSEIF( FWINDD >= 67.5 .AND. FWINDD<112.5 )THEN
          FZONE = 3
          WVANGLE = 3.142
        ELSEIF( FWINDD >= 112.5 .AND. FWINDD<157.5 )THEN
          FZONE = 4
          WVANGLE = 2.356
        ELSEIF( FWINDD >= 157.5 .AND. FWINDD<202.5 )THEN
          FZONE = 5
          WVANGLE = 1.571
        ELSEIF( FWINDD >= 202.5 .AND. FWINDD<247.5 )THEN
          FZONE = 6
          WVANGLE = 0.7854
        ELSEIF( FWINDD >= 247.5 .AND. FWINDD<292.5 )THEN
          FZONE = 7
          WVANGLE = 0.
        ELSEIF( FWINDD >= 292.5 .AND. FWINDD<337.5 )THEN
          FZONE = 8
          WVANGLE = 5.4978
        ENDIF
        FWINDSQ = FWINDS*FWINDS

        !Calculate Domain Average Depth
        !Needs to be calculated along each fetch
        !This is sufficient for small systems
        AVGDEPTH = SUM(HP(2:LA))/FLOAT(LA-1)
        
        !Calculate wave height, period, orbital velocity, and length
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,FWVHT,FC1,FC2,FWVTPP,WVLENGTH)
        DO ND = 1,NDM  
          LF = (ND-1)*LDMSED+1  
          LL = MIN(LF+LDMSED-1,LASED)
          DO LP = LF,LL
            L = LSED(LP)
            FC1   = (FWINDSQ/9.8)*0.283*TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75)
            FC2   = TANH(0.0125*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.42/TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75))   
            FWVHT = MIN(HP(L),FC1*FC2)
              
            FC1       = (FWINDS/9.8)*7.54*TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375)
            FC2       = TANH(0.077*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.25/TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375))   
            FWVTPP    = FC1*FC2
            WVFREQ(L) = 2.0*PI/FWVTPP
              
            WVLENGTH   = FWVTPP*SQRT(9.8*HP(L)/FC2)
            WVLENGTH   = MAX(1.0,WVLENGTH)
            WVORBIT(L) = MAX(0.01,PI*FWVHT/(FWVTPP*SINH(HP(L)*2.0*PI/WVLENGTH)))

            WVANG(L)   = WVANGLE
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ELSE
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND = 1,NDM  
          LF = (ND-1)*LDMSED+1  
          LL = MIN(LF+LDMSED-1,LASED)
          DO LP = LF,LL
            L = LSED(LP)
            WVFREQ(L)  = 0.0
            WVORBIT(L) = 0.0
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
      
    ELSEIF( ISWAVE > 0 .AND. ISWNWAVE == 2 )THEN
      ! *** Read in EFDC STWAVE Wave Field
      NWVCOUNT = NWVCOUNT + DTSEDJ/3600
      IF( NWVCOUNT == STWVTIM )THEN
        NWVCOUNT = 0
        STINC = STINC+1

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,FWVHT,FC1,FC2,WVLENGTH,EXCURSION)
        DO ND = 1,NDM  
          LF = (ND-1)*LDMSED+1  
          LL = MIN(LF+LDMSED-1,LASED)
          DO LP = LF,LL
            L = LSED(LP)
            IF( STINC > STWVNUM ) EXIT
            IF( STWVTP(L,STINC) > 0.0 )THEN
              WVFREQ(L)  = 2.0*PI/STWVTP(L,STINC)
              FWVHT      = MIN(HP(L),STWVHT(L,STINC))
              FC1        = 9.8*STWVTP(L,STINC)**2
              FC2        = 2.*PI*SQRT(TANH(4.*PI**2*HP(L)/(STWVTP(L,STINC)**2*9.8)))
              WVLENGTH   = FC1/FC2
              EXCURSION  = FWVHT/(2.*SINH((2.*PI/WVLENGTH)*HP(L)))
              WVORBIT(L) = EXCURSION*WVFREQ(L)
              WVORBIT(L) = MAX(0.01,WVORBIT(L))
              WVANG(L)   = STWVDR(L,STINC)
            ELSE
              WVFREQ(L)  = 0.0
              WVORBIT(L) = 0.0
              WVANG(L)   = 0.0
            ENDIF
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
      ENDIF
    ENDIF    
         
    ! **************************************************************************************************
    ! *** Begin shear stress calculations
    KNFACTOR = 1.

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, KN, UTMP, VTMP, VELMAG, FC, VELANG)   &
    !$OMP                             PRIVATE(FWW, SIGMAWV, MMW, JJW, DELW, APROUGH, SHEAR, SHEARC, SHEARW)
    DO ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      ! *** Set up roughness, velocities, and angles 
      IF( ZBSKIN <= 0. )THEN
        ! *** VARIABLE BOTTOM SKIN FRICTION.  ZBTEMP(L) = ZBR(L)  [meters]
        DO LP = LF,LL
          L = LSED(LP)
          ZBTEMP(L) = MAX(D50AVG(L)*1.E-6,1.E-6)         ! *** Convert microns to meters
        ENDDO
      ELSE
        DO LP = LF,LL
          L = LSED(LP)
          ZBTEMP(L) = ZBSKIN*1.E-6
        ENDDO
      ENDIF
      
      DO LP = LF,LL
        L = LSED(LP)
        IF( HP(L) < HPMIN )THEN
          CYCLE    ! *** Leave the shear unchanged from the last time HP >= HPMIN
        ENDIF
        
        KN = KNFACTOR*ZBTEMP(L)                                                     ! *** Nikuradse roughness (m)
        
        ! *** Calculate Cell Centroid Depth Averaged Velocity Magnitude in cm/s
        UTMP = 100.0*STCUV(L)*( UHE(LEC(L))+UHE(L) )/( SUB(LEC(L))*HU(LEC(L))+HU(L) ) + 1.0E-12
        VTMP = 100.0*STCUV(L)*( VHE(LNC(L))+VHE(L) )/( SVB(LNC(L))*HV(LNC(L))+HV(L) )
        VELMAG = SQRT(UTMP**2 + VTMP**2)                                            ! *** Current velocity (cm/s)
          
        ! *** Calculate Initial Friction Factors (Updated using Parker, 2004 approach)
        FC = ( 0.42/LOG(11.*HP(L)/(2.0*ZBTEMP(L))) )**2                             ! *** Friction factor for current only (dimensionless)
          
        ! *** Compute bed shear stresses
        IF( ISWNWAVE == 0 .AND. UWVSQ(L) == 0.0 .OR. ISWAVE ==  0 )THEN
          ! *** Compute bed shear due to currents only (dynes/cm2)
          SHEAR = FC*VELMAG**2
          QQWV3(L) = 0.0
        ELSE
          ! *** Calculate combined wave friction factor using Christoffersen & Jonsson (1985)
            
          ! *** Calculate Current Angle CCW From X axis
          IF( UTMP > 0.0 .AND. VTMP > 0.0 )THEN
            VELANG = ATAN(VTMP/UTMP)
          ELSEIF( UTMP < 0.0 )THEN
            VELANG = ATAN(VTMP/UTMP)+PI
          ELSEIF( UTMP > 0.0 .AND. VTMP < 0.0 )THEN
            VELANG = 2*PI + ATAN(VTMP/UTMP)
          ELSEIF( UTMP == 0.0 )THEN
            VELANG = SIGN(0.5*PI,VTMP)
          ENDIF
          SHEARC = FC*VELMAG**2                                                     ! *** Shear due to current only (cm2/s2)

          ! *** Calculate wave friction factor
          FWW = 2.0*(0.0747*(KN*WVFREQ(L)/WVORBIT(L)))**0.66                            ! *** Pure wave friction factor.  Equation shown bottom of pg 401, Beta = 0.0747
          SIGMAWV = FC/FWW*(VELMAG/(WVORBIT(L)*100.0))**2                               ! *** Eq 3.8
          MMW = SQRT( 1.0+SIGMAWV**2+2.0*SIGMAWV*ABS(COS(VELANG-WVANG(L))) )            ! *** Eq 3.10   Adjusting for wave/current angles
          JJW = WVORBIT(L)/(KN*WVFREQ(L))*SQRT(MMW*FWW/2.0)                             ! *** Eq 4.12
          FWW = MMW*0.15/JJW                                                            ! *** Update the adjusted wave friction coefficient. End of Step 3 Pg 412

          ! *** Calculate wave boundary layer info
          DELW = KN*0.273*SQRT(JJW)                                                     ! *** Eq 4.11 Wave boundary layer thickness (m)
          APROUGH = KNFACTOR*DELW*EXP(-5.62*DELW/KN*SQRT(SIGMAWV/MMW))                  ! *** Eq 4.23 Apparent roughness  (m)
          
          ! *** Calculate new current friction factor accounting for waves
          FC = 2.0*(1.0/(2.38*LOG(KNFACTOR*HP(L)/(2.718*KN))-2.38*LOG(APROUGH/KN)))**2  ! *** Eq 4.25 1st iteration of the combined wave and current friction factor (dimensionless)
              
          ! *** Iterate once more
          SIGMAWV = FC/FWW*(VELMAG/(WVORBIT(L)*100.0))**2                               ! *** Eq 3.8
          MMW = SQRT(1.0+SIGMAWV**2+2.0*SIGMAWV*ABS(COS(VELANG-WVANG(L))))              ! *** Eq 3.10   Adjusting for wave/current angles
          JJW = WVORBIT(L)/(KN*WVFREQ(L))*SQRT(MMW*FWW/2.0)                             ! *** Eq 4.12
          FWW = MMW*0.15/JJW
              
          ! *** Calculate total wave and current shear stress (dynes/cm^2)
          SHEARW = 0.5*FWW*MMW*WVORBIT(L)**2                                            ! *** Final wave only shear stress   (m2/s2)
          SHEARW = 10000.0*SHEARW                                                       ! *** Final wave only shear stress   (cm2/s2)
          
          QQWV3(L) = SHEARW                                                         ! *** Shear due to wave action only  (cm2/s2)
          QQWCR(L) = SHEARC                                                         ! *** Shear due to current only      (cm2/s2)
          SHEAR = SHEARC + SHEARW                                                   ! *** Shear due to current and waves (cm2/s2)
          
        ENDIF

        ! *** Limit the rate of growth to 10% of previous value for bed shear (current and wave)
        IF( SHEAR > TAU(L)*(1.+GROWTH) )THEN
          TAU(L) = TAU(L) + GROWTH*(SHEAR-TAU(L))
          IF( QQWV3(L) > 0.0 )THEN
            QQWV3(L) = TAU(L) - SHEARC 
          ENDIF
        ELSE
          TAU(L) = SHEAR
        ENDIF
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    
    ! ******************************
     
  ELSE
    ! *** Set constant tau is TAUCONST (dynes/cm^2) is greater than 0
    DO L = 2,LA
      TAU(L) = TAUCONST
      TAUB(L) = 0.1*TAU(L)
    ENDDO
  ENDIF
  
  ! *************************
  RETURN
  
END SUBROUTINE SEDZLJ_SHEAR
