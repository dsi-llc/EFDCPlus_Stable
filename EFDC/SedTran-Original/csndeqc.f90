! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
FUNCTION CSNDEQC(IOPT,SNDDIA,SSG,WS,TAUR,TAUB,D50,SIGPHI,SNDDMX,VDR,ISNDAL)  

  ! *** CALCULATES NEAR BED REFERENCE CONCENTRATION FOR NONCOHESIVE SEDIMENT  
  ! 
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-03-02        PAUL M. CRAIG    FIXED CASE WHEN WS <= 0
  !                   JOHN HAMRICK     IF USTAR < WS SET EQUILBRIUM NEAR BED CONCENTRATION TO ZERO
  !                                    TRANSPORT VIA BED LOAD IMPLEMENTED IN SSEDTOX FOR 
  !                                    USTAR**2 > CRITICAL SHIELDS AND USTAR < WS

  !     SNDDIA = SAND GRAIN DIAMETER 
  !     SSG    = SAND GRAIN SPECIFIC GRAVITY 
  !     WS     = SAND GRAIN SETTLING VELOCITY
  !     TAUR   = WATER DENSITY NORMALIZED CRITICAL SHIELDS STRESS
  !     TAUB   = WATER DENSITY NORMALIZED BED STRESS
  !     SSG    = SAND SPECIFIC GRAVITY 
  !     SIGPHI = PHI SIZE STANDARD DEVIATION
  !     SNDDMX = D90 SEDIMENT DIAMETER OR MAX SEDIMENT DIAMETER
  !
  
  implicit none
  integer :: IOPT,ISNDAL
  real :: CSNDEQC,SNDDIA,SSG,WS,TAUR,TAUB,D50,SIGPHI,SNDDMX,VDR,USTAR
  real :: REY,DFAC,RLAM,VAL,TMP,TAURS,REY3,RATIO,TMPVAL
  
  USTAR = SQRT(TAUB)
  
  if( USTAR < WS )then
     CSNDEQC = 0.
     
  elseif( IOPT == 1 )then  
    ! *** IOPT = 1  BASED ON  
    ! *** 
    ! *** GARCIA, M., AND G. PARKER, 1991: ENTRAINMENT OF BED SEDIMENT  
    ! *** INTO SUSPENSION, J. HYDRAULIC ENGINEERING, 117, 414-435.  

    if( WS > 0. )then
      REY = 1.E6*SNDDIA*SQRT( 9.8*(SSG-1.)*SNDDIA )   !EQ 42
      REY = REY**0.6                                  !SEE EQ 43
      DFAC = 1.
      if( ISNDAL >= 1) DFAC = (SNDDIA/D50)**0.2       !SEE EQ 43
      RLAM = 1.-0.29*SIGPHI                           !EQ 51
      VAL = DFAC*RLAM*REY*USTAR/WS                    !Z IN EQ 43
      VAL = 1.3E-7*(VAL**5)                           !TOP OF EQ 45
      TMP = VAL/(1+3.33*VAL)                          !EQ 45
      CSNDEQC = 1.E6*SSG*TMP                          !CONVERT TO MASS CONC
    else
      CSNDEQC = 0.
    endif

  elseif( IOPT == 2 )then  
    ! *** IOPT = 2  BASED ON  
    ! *** 
    ! *** SMITH, J. D., AND S. R. MCLEAN, 1977: SPATIALLY AVERAGED FLOW  
    ! *** OVER A WAVY SURFACE, J. GEOPHYSICAL RESEARCH, 82, 1735-1746.  

    VAL = 2.4E-3*( (TAUB/TAUR)-1. )  
    VAL = MAX(VAL,0.)  
    TMP = 0.65*VAL/(1.+VAL)  
    CSNDEQC = 1.E6*SSG*TMP  

  elseif( IOPT == 3 )then
    ! *** IOPT = 3  BASED ON  
    ! *** 
    ! *** VAN RIJN, L. C., 1984: SEDIMENT TRANSPORT, PART II: SUSPENDED  
    ! *** LOAD TRANSPORT, J. HYDRAULIC ENGINEERING, 110, 1623-1641.  

    if( WS > 0. )then
      REY = 1.E4*SNDDIA*( (9.8*(SSG-1.))**0.333 )  
      if( REY <= 10. ) TAURS = (4.*WS/REY)**2  
      if( REY  > 10. ) TAURS = 0.16*WS*WS                      ! *** Corrected 2021-06 from 0.016.  0.16 = 0.4^2 from VanRijn 1984
      REY3 = REY**0.3  
      VAL = (TAUB/TAURS)-1.  
      VAL = MAX(VAL,0.)  
      VAL = VAL**1.5  
      RATIO = SNDDIA/(3.*SNDDMX)  
      TMP = 0.015*RATIO*VAL/REY3  
      CSNDEQC = 1.E6*SSG*TMP  
    else
      CSNDEQC = 0.
    endif
    
  elseif( IOPT == 4 )then

    ! *** IOPT = 4  BASED ON 
    ! **
    ! *** J.M. HAMRICK'S parameterIZATION OF SEDFLUME DATA NO CRITICAL STRESS
    if( WS > 0. )then
      REY = 1.E4*SNDDIA*( (9.8*(SSG-1.))**0.333 )
      TMPVAL = SQRT(TAUB)/WS
      REY3 = REY**1.333
      TMPVAL = REY3*TMPVAL
      VAL = (TMPVAL-1.0)**5.
      VAL = 4.E-9*VAL
      CSNDEQC = 1.E6*SSG*VAL/(1.+VDR)
    else
      CSNDEQC = 0.
    endif

  elseif( IOPT == 5 )then
    ! *** IOPT = 5  BASED ON 
    ! **
    ! *** J.M. HAMRICK'S parameterIZATION OF SEDFLUME DATA WITH CRITICAL STRESS
    if( WS > 0. )then
      REY = 1.E4*SNDDIA*( (9.8*(SSG-1.))**0.333 )
      TMPVAL = SQRT(TAUB)/WS
      REY3 = REY**1.333
      TMPVAL = REY3*TMPVAL
      VAL = 0.0
      if( TMPVAL > 1.0 ) VAL = (TMPVAL-1.0)**5.
      VAL = 4.E-9*VAL
      CSNDEQC = 1.E6*SSG*VAL/(1.+VDR)
    else
      CSNDEQC = 0.
    endif
  else
    ! *** BAD OPTION
    call STOPP('BAD CSNDEQC OPTION') 
  endif

  return 
  END  

