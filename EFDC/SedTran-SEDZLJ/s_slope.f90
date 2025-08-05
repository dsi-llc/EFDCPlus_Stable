! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEDZLJ_SLOPE
  
  use GLOBAL
  implicit none
  
  ! CALCULATES slope change on shear
  ! 
  ! REVISION DATE :  Nov 19, 2009
  ! Craig Jones and Scott James
  ! *** ***********************************************************************
  ! Check to see if we've set a constant Tau
  integer :: L, NS
  real :: DZBETA, DZTHETA, COSB, SINB, COST, TANT, MINTAU
  real :: ALPHA, AA, BB, FCO, FB, FG, TEMP1, TUNEP, TUNER
  
  DOUBLE PRECISION ::  XDIST,YDIST,DELZX,DELZY,VELMAG,DELX,DELY

  AA = 1.0             !Geometry factor (Lick, 2009, Figure 3.25)
  BB = 1.0             !Geometry factor (Lick, 2009, Figure 3.25)
  ALPHA = ATAN(AA/BB)  !Geometry factor (Lick, 2009, p.96)
  FCO = 0.0            !Cohesive force = 7 N/m^2 * D50(K)^2, but zero for non-cohesives (Lick, 2009, p.89)
  FB = 0.0             !Binding force, but zero for non-cohesives (Lick, 2009, p.96)
  FG = 1.0             !Gravitational force (Lick, 2008, Eqn. 3.8), arbitrarily set to 1 because it is irrelevant for non-cohesives
  TUNEP = 1.0          !Bedload pitch tuning factor
  TUNER = 1.0          !Bedload roll  tuning factor
  MINTAU = 10000.
  do NS = 1,NSEDS
    if( D50(NS) >= BEDLOAD_CUTOFF ) MINTAU = MIN(TCRE(NS),MINTAU)
  enddo
  
  do L = 2,LA
    if( TAU(L) > MINTAU )then
      ! *** For Beta Angle "Pitch" and Theta Angle "Roll"
      if( SUB(LEC(L)) /= 0. .and. SUB(L) /= 0. .and. SUB(LWC(L)) /= 0. )then      ! *** Central differencing
        DELZX = BELV(LEC(L))-BELV(LWC(L))                      
        XDIST = 0.5*DXP(LEC(L))+DXP(L)+0.5*DXP(LWC(L))                            ! *** Across 3 cells
      elseif( SUB(LEC(L)) == 0. .and. SUB(L) /= 0. .and. SUB(LWC(L)) /= 0. )then  ! *** Forward differencing
        DELZX = BELV(L)-BELV(LWC(L))
        XDIST = 0.5*(DXP(L)+DXP(LWC(L)))                                          ! *** Across 2 cells
      elseif( SUB(LEC(L)) /= 0. .and. SUB(L) /= 0. .and. SUB(LWC(L)) == 0. )then
        DELZX = BELV(LEC(L))-BELV(L)
        XDIST = 0.5*(DXP(LEC(L))+DXP(L))                                            ! *** Across 2 cells
      else
        DELZX = 0.0                                                               ! *** No active cells on either side
        XDIST = 1.0
      endif
    
      if( SVB(LNC(L))/=0. .and. SVB(L) /= 0. .and. SVB(LSC(L)) /= 0. )then        ! *** Central differencing
        DELZY = BELV(LNC(L))-BELV(LSC(L))
        YDIST = 0.5*DYP(LNC(L))+DYP(L)+0.5*DYP(LSC(L))                              ! *** Across 3 cells
      elseif( SVB(LNC(L)) == 0. .and. SVB(L) /= 0. .and. SVB(LSC(L)) /= 0. )then  ! *** Forward differencing 
        DELZY = BELV(L)-BELV(LSC(L))
        YDIST = 0.5*(DYP(L)+DYP(LSC(L)))                                            ! *** Across 2 cells
      elseif( SVB(LNC(L)) /= 0. .and. SVB(L) /= 0. .and. SVB(LSC(L)) == 0. )then  ! *** Forward differencing
        DELZY = BELV(LNC(L))-BELV(L)
        YDIST = 0.5*(DYP(LNC(L))+DYP(L))                                            ! *** Across 2 cells
      else
        DELZY = 0.0                                                               ! *** No active cells on either side
        YDIST = 1.0
      endif
    
      VELMAG = SQRT(U(L,KSZU(L))**2 + V(L,KSZV(L))**2)  ! *** Local flow speed
      if( VELMAG > 0.0 )then                            ! *** Velocity directional vectors
        DELX = U(L,KSZU(L))/VELMAG
        DELY = V(L,KSZV(L))/VELMAG  
      else           
        DELX = 0.0
        DELY = 0.0
      endif  
    
      DZBETA  = DELX/XDIST*DELZX + DELY/YDIST*DELZY  ! *** Change in bed elevation for the velocity pitch angle
      DZTHETA = DELY/XDIST*DELZX + DELX/YDIST*DELZY  ! *** Change in bed elevation fot the velocity roll  angle
      COSB = 1.0/SQRT(1.0+DZBETA**2)                 ! *** Cosine of pitch angle
      SINB = DZBETA/SQRT(1.0+DZBETA**2)              ! *** Sine of pitch angle
      COST = 1.0/SQRT(1.0+DZTHETA**2)                ! *** Cosine of roll angle
      TANT = DZTHETA                                 ! *** Tangent of roll angle
      TEMP1 = (1.0+(FCO+FB)/(FG*COSB*COST))**2-(TANT/TAN(ALPHA))**2 
    
      ! *** Slope erosion from Willy Lick Work
      if( TEMP1 > 0.0 )then                          ! *** Calculate shear stress and erosion scaling factor
        SH_SCALE(L) = MAX(0.1,BB/AA*SINB+COSB*COST*SQRT(TEMP1)) ! *** Lick, 2009, Eqn.3.36
      else
        SH_SCALE(L) = 0.1                              ! *** Minimum value
      endif
    
      ! *** Slope bedload factors from Lesser et al. summary
      ALPHA_PX(L) = 1.0-TUNEP*(TAN(ALPHA)/(COS(ATAN(DELZX/XDIST))*(TAN(ALPHA)-DELZX/XDIST))-1.0) ! *** Bedload velocity pitch angle x-correction
      ALPHA_RX(L,1:NSEDS) = -TUNER*SQRT(TCRE(1:NSEDS)/TAU(L))*DELZY/YDIST                          ! *** Bedload velocity roll  angle x-correction
      ALPHA_PY(L) = 1.0-TUNEP*(TAN(ALPHA)/(COS(ATAN(DELZY/YDIST))*(TAN(ALPHA)-DELZY/YDIST))-1.0) ! *** Bedload velocity pitch angle y-correction
      ALPHA_RY(L,1:NSEDS) = -TUNER*SQRT(TCRE(1:NSEDS)/TAU(L))*DELZX/XDIST                          ! *** Bedload velocity roll  angle y-correction
    else
      SH_SCALE(L) = 1.0
      ALPHA_PX(L) = 1.0
      ALPHA_RX(L,1:NSEDS) = 0.
      ALPHA_PY(L) = 1.0
      ALPHA_RY(L,1:NSEDS) = 0.
    endif
  enddo
  
  return

END SUBROUTINE SEDZLJ_SLOPE
