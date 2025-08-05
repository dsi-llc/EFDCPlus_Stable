! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEDZLJ_SHEAR
  
  ! *** CALCULATES Wave and Current Shear Stress Based on Log-Law of Cristofferson Jonsson
  ! *** 
  ! *** ORIGINAL:  May 24, 2007
  ! ***  Craig Jones and Scott James
  ! *** REVISED: SIGMA-ZED AND OMP - 2016-11-07
  ! ***  Paul M. Craig

  use GLOBAL
  use WINDWAVE
  use Allocate_Initialize      
  
  implicit none
  
  integer :: L,ND,LF,LL,LP
  integer :: M1,M2
  integer :: FZONE
  
  real(RKD)  :: MMW, SIGMAWV, JJW, KN, FWVTPP, ABM, WA, WR, UWBM, WAMP
  real(RKD)  :: VELMAG, VELANG, DELW, APROUGH, KNFACTOR
  real(RKD)  :: UTMP, VTMP
  real(RKD)  :: WVLENGTH, WVANGLE, WFTIM, EXCURSION
  real(RKD)  :: FC1, FC2, FWINDSQ, FC, FWW, FWVHT, SHEAR, SHEARC, SHEARW, GROWTH
  real(RKD)  :: FWINDS, FWINDD
  real(RKD)  :: TDIFF, WTM1, WTM2, AVGDEPTH
  
  real(RKD) ,save,allocatable,dimension(:) :: WVFREQ
  real(RKD) ,save,allocatable,dimension(:) :: WVORBIT
  real(RKD) ,save,allocatable,dimension(:) :: ZBTEMP
  real(RKD) ,save,allocatable,dimension(:) :: WVANG

  if( .not. allocated(WVFREQ) )then
    call AllocateDSI(WVFREQ,  LCM, 0.0)
    call AllocateDSI(WVORBIT, LCM, 0.0)
    call AllocateDSI(ZBTEMP,  LCM, 0.0)
    call AllocateDSI(WVANG,  LCM, 0.0)
  endif
  GROWTH = 0.1_8
  
  ! ***************************************************************************

  if( TAUCONST == 0 )then
    ! *** Constant Tau not used.  Compute spatially variable Tau

    ! ***************************************************************************
    if( ISWAVE > 0 .and. ISWNWAVE == 0 )then
      ! *** Using EFDC+ windwave parameters
      call WINDWAVECAL    ! *** Compute all wave parameters from the fetch and wind fields
      ! *** Assign to SEDZLJ variables
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND, LF, LL, LP, L)   &
      !$OMP             SHARED(NDM, LDMSED, LASED, LSED, WV, WVFREQ, WVORBIT, UWVSQ, WVANG)
      do ND = 1,NDM  
        LF = (ND-1)*LDMSED+1  
        LL = min(LF+LDMSED-1,LASED)
        do LP = LF,LL
          L = LSED(LP)
          WVFREQ(L)  = WV(L).FREQ
          WVORBIT(L) = WV(L).UDEL
          UWVSQ(L)   = WV(L).UDEL*WV(L).UDEL
          WVANG(L)   = WV(L).DIR          
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    elseif( ISWAVE > 0 .and. ISWNWAVE == 1 )then
      ! Wind Wave Fetch
      
      ! Convert wind input into current wind info for wind-driven wave calcs   (delme - use CALSTXY winds)
      WFTIM = TIMESEC/TSWND(1).TMULT    ! *** 3tl and 2tl both use TIMESEC
      M2 = MTSWLAST(1)
      do while (WFTIM > TSWND(1).TIM(M2))
        M2 = M2+1
        if( M2 > TSWND(1).NREC )then
          M2 = TSWND(1).NREC
          EXIT
        endif
      enddo
      MTSWLAST(1) = M2 
      M1 = M2-1
      TDIFF = TSWND(1).TIM(M2)-TSWND(1).TIM(M1)  
      WTM1 = (TSWND(1).TIM(M2)-WFTIM)/TDIFF  
      WTM2 = (WFTIM-TSWND(1).TIM(M1))/TDIFF 
      FWINDS = WTM1*TSWND(1).VAL(M1,1)+WTM2*TSWND(1).VAL(M2,1)

      if( FWINDS > 1.0 )then
        ! *** Only compute wind wave if wind speed > 1 m/s 
        if( ABS(TSWND(1).VAL(M1,2)-TSWND(1).VAL(M2,2))<180.0 )then
          FWINDD = WTM1*TSWND(1).VAL(M1,2)+WTM2*TSWND(1).VAL(M2,2)
        else
          if( TSWND(1).VAL(M1,2) > TSWND(1).VAL(M2,2) )then
              FWINDD = WTM1*TSWND(1).VAL(M1,2)+WTM2*(TSWND(1).VAL(M2,2)+360.0)
          else
              FWINDD = WTM1*(TSWND(1).VAL(M1,2)+360)+WTM2*TSWND(1).VAL(M2,2)
          endif
          if( FWINDD >= 360.0)FWINDD = FWINDD-360.0 
        endif
        ! Convert wind into direction it is blowing "from"
        if( FWINDD <= 180.0 )then  
          FWINDD = FWINDD+180.0  
          if( FWINDD == 360.0)FWINDD = 0.0         
        else  
          FWINDD = FWINDD-180.0  
          if( FWINDD == 360.0)FWINDD = 0.0 
        endif
        
        ! *** Calculate which of the 8 wind zones (FWZONE) the wind is coming from
        ! *** Also the Waveangle CCW from East.  Waveangle is center of sector.
        if( FWINDD >= 337.5 .or. FWINDD<22.5 )then
          FZONE = 1
          WVANGLE = 4.712
        elseif( FWINDD >= 22.5 .and. FWINDD<67.5 )then
          FZONE = 2
          WVANGLE = 3.927
        elseif( FWINDD >= 67.5 .and. FWINDD<112.5 )then
          FZONE = 3
          WVANGLE = 3.142
        elseif( FWINDD >= 112.5 .and. FWINDD<157.5 )then
          FZONE = 4
          WVANGLE = 2.356
        elseif( FWINDD >= 157.5 .and. FWINDD<202.5 )then
          FZONE = 5
          WVANGLE = 1.571
        elseif( FWINDD >= 202.5 .and. FWINDD<247.5 )then
          FZONE = 6
          WVANGLE = 0.7854
        elseif( FWINDD >= 247.5 .and. FWINDD<292.5 )then
          FZONE = 7
          WVANGLE = 0.
        elseif( FWINDD >= 292.5 .and. FWINDD<337.5 )then
          FZONE = 8
          WVANGLE = 5.4978
        endif
        FWINDSQ = FWINDS*FWINDS

        !Calculate Domain Average Depth
        !Needs to be calculated along each fetch
        !This is sufficient for small systems
        AVGDEPTH = SUM(HP(2:LA))/FLOAT(LA-1)
        
        !Calculate wave height, period, orbital velocity, and length
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,FWVHT,FC1,FC2,FWVTPP,WVLENGTH)
        do ND = 1,NDM  
          LF = (ND-1)*LDMSED+1  
          LL = min(LF+LDMSED-1,LASED)
          do LP = LF,LL
            L = LSED(LP)
            FC1   = (FWINDSQ/9.8)*0.283*TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75)
            FC2   = TANH(0.0125*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.42/TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75))   
            FWVHT = min(HP(L),FC1*FC2)
              
            FC1       = (FWINDS/9.8)*7.54*TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375)
            FC2       = TANH(0.077*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.25/TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375))   
            FWVTPP    = FC1*FC2
            WVFREQ(L) = 2.0*PI/FWVTPP
              
            WVLENGTH   = FWVTPP*SQRT(9.8*HP(L)/FC2)
            WVLENGTH   = max(1.0,WVLENGTH)
            WVORBIT(L) = max(0.01,PI*FWVHT/(FWVTPP*SINH(HP(L)*2.0*PI/WVLENGTH)))

            WVANG(L)   = WVANGLE
          enddo
        enddo
        !$OMP END PARALLEL DO
      else
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM  
          LF = (ND-1)*LDMSED+1  
          LL = min(LF+LDMSED-1,LASED)
          do LP = LF,LL
            L = LSED(LP)
            WVFREQ(L)  = 0.0
            WVORBIT(L) = 0.0
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      
    elseif( ISWAVE > 0 .and. ISWNWAVE == 2 )then
      ! *** Read in EFDC STWAVE Wave Field
      NWVCOUNT = NWVCOUNT + DTSEDJ/3600
      if( NWVCOUNT == STWVTIM )then
        NWVCOUNT = 0
        STINC = STINC+1

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,FWVHT,FC1,FC2,WVLENGTH,EXCURSION)
        do ND = 1,NDM  
          LF = (ND-1)*LDMSED+1  
          LL = min(LF+LDMSED-1,LASED)
          do LP = LF,LL
            L = LSED(LP)
            if( STINC > STWVNUM ) EXIT
            if( STWVTP(L,STINC) > 0.0 )then
              WVFREQ(L)  = 2.0*PI/STWVTP(L,STINC)
              FWVHT      = min(HP(L),STWVHT(L,STINC))
              FC1        = 9.8*STWVTP(L,STINC)**2
              FC2        = 2.*PI*SQRT(TANH(4.*PI**2*HP(L)/(STWVTP(L,STINC)**2*9.8)))
              WVLENGTH   = FC1/FC2
              EXCURSION  = FWVHT/(2.*SINH((2.*PI/WVLENGTH)*HP(L)))
              WVORBIT(L) = EXCURSION*WVFREQ(L)
              WVORBIT(L) = max(0.01,WVORBIT(L))
              WVANG(L)   = STWVDR(L,STINC)
            else
              WVFREQ(L)  = 0.0
              WVORBIT(L) = 0.0
              WVANG(L)   = 0.0
            endif
          enddo
        enddo
        !$OMP END PARALLEL DO
        
      endif
    endif    
         
    ! **************************************************************************************************
    ! *** Begin shear stress calculations
    KNFACTOR = 1.

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, KN, UTMP, VTMP, VELMAG, FC, VELANG)   &
    !$OMP                             PRIVATE(FWW, SIGMAWV, MMW, JJW, DELW, APROUGH, SHEAR, SHEARC, SHEARW)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      ! *** Set up roughness, velocities, and angles 
      if( ZBSKIN <= 0. )then
        ! *** VARIABLE BOTTOM SKIN FRICTION.  ZBTEMP(L) = ZBR(L)  [meters]
        do LP = LF,LL
          L = LSED(LP)
          ZBTEMP(L) = max(D50AVG(L)*1.E-6,1.E-6)         ! *** Convert microns to meters
        enddo
      else
        do LP = LF,LL
          L = LSED(LP)
          ZBTEMP(L) = ZBSKIN*1.E-6
        enddo
      endif
      
      do LP = LF,LL
        L = LSED(LP)
        if( HP(L) < HPMIN )then
          CYCLE    ! *** Leave the shear unchanged from the last time HP >= HPMIN
        endif
        
        KN = KNFACTOR*ZBTEMP(L)                                                     ! *** Nikuradse roughness (m)
        
        ! *** Calculate Cell Centroid Depth Averaged Velocity Magnitude in cm/s
        UTMP = 100.0*STCUV(L)*( UHE(LEC(L))+UHE(L) )/( SUB(LEC(L))*HU(LEC(L))+HU(L) ) + 1.0E-12
        VTMP = 100.0*STCUV(L)*( VHE(LNC(L))+VHE(L) )/( SVB(LNC(L))*HV(LNC(L))+HV(L) )
        VELMAG = SQRT(UTMP**2 + VTMP**2)                                            ! *** Current velocity (cm/s)
          
        ! *** Calculate Initial Friction Factors (Updated using Parker, 2004 approach)
        FC = ( 0.42/LOG(11.*HP(L)/(2.0*ZBTEMP(L))) )**2                             ! *** Friction factor for current only (dimensionless)
          
        ! *** Compute bed shear stresses
        if( (ISWNWAVE == 0 .and. WVORBIT(L) < 1.D-6) .or. ISWAVE ==  0 )then
          ! *** Compute bed shear due to currents only (dynes/cm2)
          SHEAR = FC*VELMAG**2
          QQWV3(L) = 0.0
        else
          ! *** Calculate combined wave friction factor using Christoffersen & Jonsson (1985)
            
          ! *** Calculate Current Angle CCW From X axis
          VELANG = 0.0
          if( UTMP > 0.0 .and. VTMP > 0.0 )then
            VELANG = ATAN(VTMP/UTMP)
          elseif( UTMP < 0.0 .and. VTMP > 0.0 )then
            VELANG = ATAN(VTMP/UTMP)+PI
          elseif( UTMP > 0.0 .and. VTMP < 0.0 )then
            VELANG = 2*PI + ATAN(VTMP/UTMP)
          elseif( UTMP < 0.0 .and. VTMP < 0.0 )then
            VELANG = PI + ATAN(VTMP/UTMP)
          endif
          SHEARC = FC*VELMAG**2                                                         ! *** Shear due to current only (cm2/s2)

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
          
          QQWV3(L) = SHEARW                                                             ! *** Shear due to wave action only  (cm2/s2)
          QQWCR(L) = SHEARC                                                             ! *** Shear due to current only      (cm2/s2)
          SHEAR = SHEARC + SHEARW                                                       ! *** Shear due to current and waves (cm2/s2)
          
        endif

        ! *** Limit the rate of growth to 10% of previous value for bed shear (current and wave)
        if( SHEAR > TAU(L)*(1.+GROWTH) )then
          TAU(L) = TAU(L) + GROWTH*(SHEAR-TAU(L))
          if( QQWV3(L) > 0.0 )then
            QQWV3(L) = TAU(L) - SHEARC 
          endif
        else
          TAU(L) = SHEAR
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! ******************************
     
  else
    ! *** Set constant tau is TAUCONST (dynes/cm^2) is greater than 0
    do L = 2,LA
      TAU(L) = TAUCONST
      TAUB(L) = 0.1*TAU(L)
    enddo
  endif
  
  ! *************************
  return
  
END SUBROUTINE SEDZLJ_SHEAR
