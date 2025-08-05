! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!< @details Calculates erosion rate (g/m2/s) for a given shear stress (TAUP (Pascals))
!< @author  Paul Craig
!< @param[in]  L
!< @param[in]  TAUP
!< @param[out] ELAY, SURFACE
!---------------------------------------------------------------------------!
Subroutine Calc_Prop_Erosion_SEDZLJ(L, TAUP, ELAY, SURFACE)

  use GLOBAL

  implicit none

  ! *** Dummy variables
  integer, intent(in)           :: L           !> cell index
  real(kind = rkd), intent(in)  :: TAUP        !> shear stress                  N/m**2
  real(kind = rkd), intent(out) :: ELAY(NSEDS)  !> calculated erosion per class  gm/cm**2
  integer, intent(out) :: SURFACE
  
  ! *** Local variables
  integer :: K, NS, SLLN, NACTLAY, NT, K1, SURFOLD
  integer :: NSC0, NSC1, NTAU0, NTAU1, ICORE

  real(RKD) :: CSEDSS,SQR2PI
  real(RKD) :: D50TMPP, D50AVGL
  real(RKD) :: DEP
  real(RKD) :: EBD, ERO
  real(RKD) :: ERATEMOD
  real(RKD) :: NSCTOT
  real(RKD) :: ONE = 1.0
  real(RKD) :: PFY
  real(RKD) :: PX
  real(RKD) :: PY
  real(RKD) :: SEDFLUX
  real(RKD) :: SN00
  real(RKD) :: SN01
  real(RKD) :: SN10
  real(RKD) :: SN11
  real(RKD) :: TEMP,TEMP1,TEMP2
  real(RKD) :: TACT,TSUM,TAUDYNE
  real(RKD) :: TAUCRIT
  real(RKD) :: VZDIF
  real(RKD) :: WDTDZ
  real(RKD) :: ERATEMAX

  real(RKD) ,dimension(NSEDS) :: QBFLUX
  real(RKD) ,dimension(NSEDS) :: CSEDVR
  real(RKD) ,dimension(NSEDS) :: CTB
  real(RKD) ,dimension(NSEDS) :: DEPBL
  real(RKD) ,dimension(NSEDS) :: DEPTSS
  real(RKD) ,dimension(NSEDS) :: ETOT
  real(RKD) ,dimension(NSEDS) :: PROB
  real(RKD) ,dimension(NSEDS) :: PROBVR
  real(RKD) ,dimension(NSEDS) :: SMASS
  real(RKD) ,dimension(NSEDS) :: TTEMP
  real(RKD) ,dimension(KB)   :: INITMASS

  real(RKD) ,dimension(2)    :: NSCD
  real(RKD) ,dimension(2)    :: TAUDD

  ! *** HARD BOTTOM BYPASS HANDLED IN CALLING ROUTINE

  ! *** ********************************************************************************************************************************************
  ! *** Get things set up for Erosion Calculations
  SURFACE = KB+1
  do K = 2,KB                       ! *** Ignore active layer
    if( TSED(K,L) > 0.001 )then   ! *** Ignore layers than are less than 0.1 mm (typically)
      SURFACE = K
      EXIT
    endif
  enddo

  TAUDYNE = 10.*TAUP              ! *** Convert to dynes
  if( SURFACE > KB )then
    ! *** Set a tiny concentration as a flag to be used in SedTran integration if shear stress is high
    WHERE( TAUDYNE >= TAUCRITE(1:NSEDS) )
      ELAY(1:NSEDS)  = 1.E-12
    ENDWHERE
    return                       ! *** Insufficient sediments
  endif

  ! *** *******************************************************************
  ICORE = NCORENO(IL(L),JL(L))

  ! *** Calculate erosion/deposition for top layer for all cells
  ERO = 0.0                       ! *** Initialize total erosion for the cell
  INITMASS(1:KB) = TSED(1:KB,L)   ! *** Save the starting sediment mass by layers

  ! *** Calculate Average particle size of surface layer so we can calculate
  ! *** active layer unit mass
  D50AVGL = SUM(PERSED(1:NSEDS,SURFACE,L)*D50(1:NSEDS))             ! *** Calculate local d50 at sediment bed surface

  ! *** Calculate TAUCRIT Based on the Average Particle Size of Surface
  ! *** Then calculate the Active Layer unit mass (TACT) from it.
  ! *** Ta =  Tam * Davg * (Tau/Taucr)
  if( LAYERACTIVE(SURFACE,L) < 2 )then
    ! Identify Size Class interval to use for Taucrit erosion calculation
    if(D50AVGL < SCND(1) )then
      NSCD(1) = SCND(1)
      NSCD(2) = SCND(2)
      NSC0 = 1
      NSC1 = 2
      D50AVGL = SCND(1)                                             ! *** Prevent division (s_shear) by zero when there is no sediment in the layer
    elseif( D50AVGL >= SCND(NSICM) )then
      NSCD(1) = SCND(NSICM-1)
      NSCD(2) = SCND(NSICM)
      NSC0 = NSICM-1
      NSC1 = NSICM
    else
      do NS = 1,NSICM-1
        if( D50AVGL >= SCND(NS) .and. D50AVGL < SCND(NS+1) )then
          NSCD(1) = SCND(NS)
          NSCD(2) = SCND(NS+1)
          NSC0 = NS
          NSC1 = NS+1
          EXIT
        endif
      enddo
    endif

    TAUCRIT = TAUCRITE(NSC0)+(TAUCRITE(NSC1)-TAUCRITE(NSC0))/(NSCD(2)-NSCD(1))*(D50AVGL-NSCD(1))
  else
    ! *** IN-PLACE SEDIMENTS
    TAUCRIT = TAUCOR(SURFACE,L)
  endif

  ! *** Compute the requried active layer thickness (cm)
  if( TAUDYNE < TAUCRIT )then
    return
  else
    TACT = TACTM*D50AVGL*(TAUDYNE/TAUCRIT)*(BULKDENS(1,L)/10000.0)
  endif

  ! *** ********************************************************************************************************************************************
  ! *** Now calculate the Erosion Rate
  K = SURFACE

  ! *** Find upper and lower limits of size classes on mean bed diameter
  if( (D50AVGL+1E-6) < SCND(1) )then
    NS = 1
    NSCD(1) = SCND(NS)
    NSCD(2) = SCND(NS+1)
    NSC0 = NS
    NSC1 = NS+1
    D50AVGL = SCND(1)
  elseif( (D50AVGL-1E-6) > SCND(NSICM) )then
    NS = NSICM - 1
    NSCD(1) = SCND(NS)
    NSCD(2) = SCND(NS+1)
    NSC0 = NS
    NSC1 = NS+1
    D50AVGL = SCND(NSICM)
  else
    do NS = 1,NSICM-1
      if( D50AVGL >= SCND(NS) .and. D50AVGL < SCND(NS+1) )then
        NSCD(1) = SCND(NS)
        NSCD(2) = SCND(NS+1)
        NSC0 = NS
        NSC1 = NS+1
        EXIT
      endif
    enddo
  endif

  ! *** Calculate TAUCRIT Based on the D50 of the bed or from Sedflume Data
  if( LAYERACTIVE(SURFACE,L) < 2 )then                ! *** For active/deposited layers
    TAUCRIT = TAUCRITE(NSC0) + (TAUCRITE(NSC1)-TAUCRITE(NSC0))/(NSCD(2)-NSCD(1))*(D50AVGL-NSCD(1)) !interpolation
    TAUCOR(K,L) = TAUCRIT
  else
    ! *** SEDFlume data (depth interpolation)
    SN01 = TSED(K,L)/TSED0(K,L)                       ! *** Weighting factor 1 for interpolation
    SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)          ! *** Weighting factor 2
    TAUCRIT = SN01*TAUCOR(K,L) + SN11*TAUCOR(K+1,L)
  endif

  ! *** Check if the shear is greater than critical shears.  If not, exit erosion FUNCTION
  if( TAUDYNE < TAUCRIT ) return

  ! *** Now, calculate erosion rates

  ! *** Find the upper and lower limits of the Shear Stress for the interpolation
  if( NSEDFLUME == 1 )then

    if( TAUDYNE >= TAULOC(ITBM) )then
      TAUDD(1) = TAULOC(ITBM-1)
      TAUDD(2) = TAULOC(ITBM)
      NTAU0 = ITBM-1
      NTAU1 = ITBM

    elseif( TAUDYNE < TAULOC(1) )then
      TAUDD(1) = TAULOC(1)
      TAUDD(2) = TAULOC(2)
      NTAU0 = 1
      NTAU1 = 2
    else
      do NS = 1,ITBM-1
        if( TAUDYNE >= TAULOC(NS) .and. TAUDYNE < TAULOC(NS+1) )then
          TAUDD(1) = TAULOC(NS)
          TAUDD(2) = TAULOC(NS+1)
          NTAU0 = NS
          NTAU1 = NS+1
          EXIT
        endif
      enddo
    endif

  endif

  ! *** Interpolate the erosion rates for shear stress and depth.
  ! *** This utilizes normal sedflume data for deeper layers.
  if( LAYERACTIVE(SURFACE, L) == 2 )then
    ! *** Calculate erosion rates of deeper layers (SEDFlume data)
    if( NSEDFLUME == 1 )then
      SN00 = (TAUDD(2)-TAUDYNE)/(TAUDD(2)-TAUDD(1)) ! *** weighting factor 1 for interpolation
      SN10 = (TAUDD(1)-TAUDYNE)/(TAUDD(1)-TAUDD(2)) ! *** weighting factor 2
      SN01 = TSED(K,L)/TSED0(K,L)                   ! *** weighting factor 3
      SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)      ! *** weighting factor 4

      if( K+1 <= KB )then  ! *** Maximum erosion rate
        ERATEMAX = ( SN00*EXP(SN11*LOG(ERATE(K+1,L,NTAU0))+SN01*LOG(ERATE(K,L,NTAU0))) &
          + SN10*EXP(SN11*LOG(ERATE(K+1,L,NTAU1))+SN01*LOG(ERATE(K,L,NTAU1))) )*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))
      else                 ! *** Do not allow erosion through the bottom layer
        ERATEMAX = ( SN00*EXP(SN11*LOG(1.0E-9)+SN01*LOG(ERATE(K,L,NTAU0))) &
          + SN10*EXP(SN11*LOG(1.0E-9)+SN01*LOG(ERATE(K,L,NTAU1))) )*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))
      endif
    else
      if( TAUDYNE > TAUCOR(K,L) )then                                                  ! *** Check that the applied shear exceeds the critical shear stress for this layer
        ! *** Erosion rate values (cm/s) computed by equation assume shear in Pascals so convert dynes
        SN00 = EA(ICORE,K)*TAUP**EN(ICORE,K)                                           ! *** Erosion rate (cm/s) of the top layer

        if( K+1 <= KB )then
          SN10 = EA(ICORE,K+1)*TAUP**EN(ICORE,K+1)                                     ! *** Erosion rate (cm/s) of the layer below
        else
          SN10 = 0.0                                                                   ! *** Modeled erosion rate in limited by bottom
        endif

        SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)                                       ! *** Mass weighting factor
        ERATEMAX = ((SN10-SN00)*SN11 + SN00)*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))       ! *** linear interpolation for remaining mass in current layer    (g/cm2/s)
        ERATEMAX = MIN(ERATEMAX,MAXRATE(ICORE,K))                                      ! *** Limit erosion rate
      else
        ERATEMAX = 0.0
      endif
    endif
  else
    ! *** For Active and deposited sediment layers
    ! *** The erosion rate for these layers is determined from
    ! *** Sedflume experiments and is based on average particle Size (D50AVG)
    ! ***
    NSCTOT = NSCD(2)-NSCD(1)                                                         ! *** difference in interpolant size class
    D50TMPP = D50AVGL-NSCD(1)                                                      ! *** difference from local size class and lower interpolant
    if( NSEDFLUME == 1 )then
      SN00 = (TAUDD(2)-TAUDYNE)/(TAUDD(2)-TAUDD(1))                                  ! *** weighting factor 1 for interpolation
      SN10 = (TAUDD(1)-TAUDYNE)/(TAUDD(1)-TAUDD(2))                                  ! *** weigthing factor 2
      SN01 = D50TMPP/NSCTOT                                                          ! *** weighting factor 3
      SN11 = (NSCTOT-D50TMPP)/NSCTOT                                                 ! *** weighting factor 4
      ERATEMAX = (SN00*EXP(SN11*LOG(ERATEND(NSC0,NTAU0)) + SN01*LOG(ERATEND(NSC1,NTAU0))) + SN10*EXP(SN11*LOG(ERATEND(NSC0,NTAU1)) +  &   ! *** log-linear interpolation
        SN01*LOG(ERATEND(NSC1,NTAU1))))*BULKDENS(K,L)*SQRT(1./SH_SCALE(L))
    else
      ! *** Erosion rate values (cm/s) computed by equation assume shear in Pascals so convert dynes
      SN00 = ACTDEPA(NSC0)*TAUP**ACTDEPN(NSC0)                                       ! *** Erosion rate 1 (cm/s)
      SN10 = ACTDEPA(NSC1)*TAUP**ACTDEPN(NSC1)                                       ! *** Erosion rate 2 (cm/s)
      SN11 = D50TMPP/NSCTOT                                                          ! *** Weighting factor
      ERATEMAX = ((SN10-SN00)*SN11 + SN00)*BULKDENS(K,L)*SQRT(1./SH_SCALE(L))        ! *** linear interpolation around size class (g/cm2/s)
      ERATEMAX = MIN(ERATEMAX,ACTDEPMAX(NSC0))                                       ! *** Limit erosion rate
    endif
  endif

  ! *** Sort out Thicknesses and Erosion Rates
  EBD = ERATEMAX*DTSEDJ                                                              ! *** Maximum mass potentially eroded this time step for this layer (g/cm^2)

  ! *** If the shear stress is less than the critical shear stress for a
  ! *** particular size class, then it is not eroded from the bed.

  ! *** Conservation of sediment mass will be addressed on a full cell basis in SEDZLJ
  ! *** ELAY(NS) = Total erosion at this cell of size class NS
  ! *** ERO      = Total erosion at this cell
  do NS = 1,NSEDS
    if( TAUDYNE >= TCRE(ns) )then
      ELAY(NS)  = PERSED(NS,K,L)*EBD
    else
      ELAY(NS)  = 0.0
    endif
  enddo
  
  !ERO = SUM(ELAY(1:NSEDS))                                                            ! *** Total erosion from the layer   (g/cm^2)

  return

  END Subroutine Calc_Prop_Erosion_SEDZLJ

  
!---------------------------------------------------------------------------!
!< @details  CALCULATES EROSION RATE (GM/M2/S) FOR A GIVEN
!! SHEAR STRESS (TAUP (Pascals))
!< @author  Paul Craig
!< @param[in] L
!< @param[in] TAUP
!< @param[out] ELAY   Erosion into suspended load by class
!< @param[out] EBLD   Erosion into bedload by class
!---------------------------------------------------------------------------!
Subroutine Calc_Prop_Erosion_Original(L, TAUP, ELAY, EBLD)

  ! ***
  use GLOBAL

  implicit none

  ! *** Dummy variables
  integer, intent(in)           :: L           !> cell index
  real(kind = rkd), intent(in)  :: TAUP        !> density normalized shear stress               m2/s2
  real(kind = rkd), intent(out) :: ELAY(NSEDS)  !> calculated erosion into suspension per class  g/m2
  real(kind = rkd), intent(out) :: EBLD(NSEDS)  !> calculated erosion into bedload per class     g/m2
  
  integer :: NS, NX
  real :: CSHIELDS, FACBEDLP, FACSUSLP, ROUSEP, SHIELDS, SNDEQP, SIGP, SNDEQBP, TAUE, TAUBSNDP, TAURTMP
  real :: TMPSEDHID, USTARP, USTARSNDP, WESE, WESEMX, WSFAC, ZEQP, ZEQDP, ZEQDIP, ZEQMAX
  real :: BDLDE, BDLDTMPP, BDLDTMP, BDLDTMPA, BDLDTMPB, CSHIELDSC
  real :: BOT, TOP, TMPSTR, TMPVAL, DIASEDPROP, DSEDGMMPROP, GPDIASEDPROP
  real, external :: FSBDLD, FSEDMODE, CSNDZEQ, CSNDEQC
  
  ! *** Cohesive sediments
  if( ISTRAN(6) > 0 )then
    if( IWRSP(1) == 4 )then
      TMPSTR = 0.0
    else
      TMPSTR = 1.0
    endif    

    do NS = 1,NSED
      ! *** SET MAXIMUM EROSION RATE
      WESEMX = DELTI*SEDB(L,KBT(L),NS)

      if( TAUP > TAURB(L,KBT(L)) )then
        ! *** MASS/BULK EROSION
        WESE = WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
        WESE = MIN(WESE,WESEMX)
      else
        TAUE = 0.
        if( TAUP > TAURS(L,KBT(L)) )then
          ! *** SURFACE EROSION
          WESE = WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)

          ! *** SET NORMALIZING SHEAR BASED ON INPUT OPTION
          if( IWRSP(1) >= 99 )then
            ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU FROM SSCOHSEDPMAP
            TAURTMP = TAUNS(L,KBT(L))
          elseif( (IWRSP(1) >= 2) .and. (IWRSP(1) < 4) )then
            ! *** 1 HWANG AND METHA - LAKE OKEECHOBEE
            ! *** 2 HAMRICK'S MODIFICATION OF SANFORD AND MAA
            ! *** 3 SAME AS 2 EXCEPT VOID RATIO OF COHESIVE SEDIMENT FRACTION IS USED
            TAURTMP = TAUR(1)
          else
            ! *** Use DIRECTLY INPUT PROPERTIES
            TAURTMP = TAURS(L,KBT(L))
          endif

          TAUE = (TAUP-TMPSTR*TAURS(L,KBT(L)))/TAURTMP
          TAUE = MAX(TAUE,0.0)

          ! *** SET NON-COHESIVE HIDING FACTOR
          TMPSEDHID = 1.0
          if( ISTRAN(7) >= 1 .and. COSEDHID(1) /= 0.0 )then
            TMPSEDHID = (FRACCOH(L,KBT(L)))**COSEDHID(1)
          endif
              
          ! *** FINALIZE SURFACE EROSION RATE
          if( IWRSP(1) < 99 )then
            WESE = TMPSEDHID*WESE*( TAUE**TEXP(NS) )
          else
            ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU EXPONENT FROM SSCOHSEDPMAP
            WESE = TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
          endif
          WESE = MIN(WESE,WESEMX)
        else
          ! *** NO EROSION 
          WESE = 0.0
        endif
        ELAY(NS) = WESE                                            ! *** Erosion flux in g/m2/s
      endif
    enddo
  endif
  
  ! *** Non-Cohesive sediments
  if( ISTRAN(7) > 0 )then
    
    ! *** Erosion flux into the water column
    TAUBSNDP = TAUP                                                ! *** Propwash shear.  use total shear m2/s2
    USTARP = SQRT(TAUP)                                            ! *** USTAR m/s
    USTARSNDP = USTARP
    
    do NX = 1,NSND
      NS = NSED + NX
      
      DSEDGMMPROP  = 1./(1.E6*SSG(NS))                             ! *** Specific volume (m**3/g)
      DIASEDPROP   = SEDDIA(NS)                                    ! *** Sediment particle diameter (m)
      GPDIASEDPROP = G*(SSG(NS)-1.)*DIASEDPROP  
      
      ! *** SET THE ROUSE parameter USING TOTAL
      if( USTARP > 1E-12 )then  
        ROUSEP = WSETA(L,0,NS)/(VKC*USTARP)  
      else  
        ROUSEP = 250000.*WSETA(L,0,NS)  
      endif  
      
      ! *** FSEDMODE SETS BEDLOAD (IMODE = 1) AND SUSPENDED LOAD (IMODE = 2) TRANSPORT FRACTIONS
      FACSUSLP = FSEDMODE(WSETA(L,0,NS),USTARP,USTARSNDP,RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)  

      WESE = 0.
      if( ISNDEQ(NS) == 0 )then
        ! *** USER SPECIFIED EQUILIBRIUM CONCENTRATION
        ZEQP    = 0.05         ! *** Garcia & Parker
        ZEQDP   = ZEQP*DZIC(L,KSZ(L))  
        ZEQDIP  = 1./ZEQDP  
        SNDEQBP = TAUR(NS)
      else
        ! *** SET EQUILIBRUIM CONCENTRATION REFERENCE HEIGHT (DIMENSIONLESS)  
        ZEQP = CSNDZEQ(ISNDEQ(NS),DIASEDPROP,GPDIASEDPROP,TAUR(NS),TAUBSNDP,SEDDIA50(L,KBT(L)),HP(L),SSG(NS),WSETA(L,0,NS))
          
        ZEQMAX   = 0.5*DZC(L,KSZ(L))  
        ZEQP   = MIN(ZEQP,ZEQMAX)  
        ZEQDP  = ZEQP*DZIC(L,KSZ(L))  
        ZEQDIP = 1./ZEQDP  
        SIGP = SIGPHI(L,KBT(L))  

        ! *** SET EQUILIBRUIM CONCENTRATION  
        SNDEQBP = CSNDEQC(ISNDEQ(NS),DIASEDPROP,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUBSNDP,SEDDIA50(L,KBT(L)),SIGP,ZEQP,VDRBED(L,KBT(L)),ISNDAL)      !  

        ! *** GET CRITICAL SHIELDS STRESS FOR GARCIA AND PARKER  
        if( ISNDEQ(NS) == 1 )then  
          if( ISLTAUC(NS) <= 1 )then  
            CSHIELDS = TCSHIELDS(NS)  
          elseif( ISLTAUC(NS) == 2 )then  
            CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASEDPROP  
          elseif( ISLTAUC(NS) == 3 )then  
            CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)  
          else
            CSHIELDS = TCSHIELDS(NS)  
          endif  

          ! *** Modify crtical shear stress to account for cohesive fraction
          if( ISEDEFF == 2 )then  
            TMPVAL = 1. + (COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
            CSHIELDS = TMPVAL*CSHIELDS  
          endif  
          SHIELDS = TAUBSNDP/GPDIASEDPROP  
          if( SHIELDS < CSHIELDS ) SNDEQBP = 0.  
        endif  

        ! *** Modify crtical erosion rates to account for cohesive fraction
        if( ISEDEFF == 1 )then
          SNDEQBP = SNDEQBP*EXP(-COEHEFF*FRACCOH(L,KBT(L)))  
        endif
      endif
        
      if( ROUSEP < 0.999 .or. ROUSEP > 1.001 )then  
        TOP = (ZEQDP**(ROUSEP-1.))-1.  
        BOT = (1.-ROUSEP)*(ZEQDIP-1.)  
        SNDEQP = SNDEQBP*TOP/BOT  
        SNDEQP = FACSUSLP*VFRBED(L,KBT(L),NS)*MAX(SNDEQP,0.)  
      else  
        TOP = LOG(ZEQDIP)  
        BOT = (ZEQDIP-1.)  
        SNDEQP = SNDEQBP*TOP/BOT  
        SNDEQP = FACSUSLP*VFRBED(L,KBT(L),NS)*MAX(SNDEQP,0.)  
      endif  

      ! *** SET RESUSPENSION FLUX
      if( KSZ(L) == KC )then
        WSFAC = 2.*(1.+ROUSEP)/(2.+ROUSEP*(1.-ZEQP))               ! *** Single layer
      else
        WSFAC = 1.
      endif
      WESE = WSFAC*WSETA(L,0,NS)*SNDEQP                            ! *** Erosion flux in g/m2/s  
      ELAY(NS) = WESE*DTSEDJ                                       ! *** Erosion flux in g/m2
    enddo
    
    ! *** Erosion flux into bedload
    if( ICALC_BL > 0 )then
      do NX = 1,NSND
        NS = NSED + NX
      
        DSEDGMMPROP  = 1./(1.E6*SSG(NS))                           ! *** Specific volume (m**3/g)
        DIASEDPROP   = SEDDIA(NS)                                  ! *** Sediment particle diameter (m)
        GPDIASEDPROP = G*(SSG(NS)-1.)*DIASEDPROP                   ! *** "excess" density (m2/s2)
        BDLDTMP = SQRT(GPDIASEDPROP)*DIASEDPROP/DSEDGMMPROP

        CSHIELDS = TCSHIELDS(NS)
      
        ! *** Modify crtical shear stress to account for cohesive fraction
        if( ISEDEFF == 2 )then
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        endif
      
        !*** Compute the cell centered bedload transport rates using the specified option
        BDLDE = 0.0
        if( ISBDLD(NS) == 0 )then
          ! *** Calculate cell center transport rates using generic bed load equation
          BDLDTMPP = SBDLDP(NX)                                                          ! *** Constant Phi
          if( BDLDTMPP > 0.0 )then
            FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP ,USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
            SHIELDS = TAUBSNDP/GPDIASEDPROP
            if( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDS )then
              if( SBDLDA(NX) > 0.0 )then
                BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
              else
                BDLDTMPA = 1.0
              endif
              if( SBDLDB(NX) > 0.0 )then
                BDLDTMPB = (SBDLDG3(NX)*SQRT(SHIELDS) - SBDLDG4(NX)*SQRT(CSHIELDS))**SBDLDB(NX)
              else
                BDLDTMPB = 1.0
              endif
              BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA*BDLDTMPB    ! *** (g/m/s)
            endif
          endif

        elseif( ISBDLD(NS) == 1 )then
          ! *** Calculate cell center transport rates using Van Rijn bed load equation

          ! *** Activate Garcia & Parker armoring
          if( ISNDAL == 1 )then
            TMPVAL = LOG10(19.*DIASEDPROP/SEDDIA50(L,KBT(L)))
            TMPVAL = 1.66667/(TMPVAL**2)
            CSHIELDSC = CSHIELDS50(L)*TMPVAL
          else
            CSHIELDSC = TCSHIELDS(NS)
          endif

          BDLDTMPP = FSBDLD(DIASEDPROP, GPDIASEDPROP, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))   ! *** Calculate Phi
          if( ISNDAL == 1 )then
            BDLDTMPP = ((DIASEDPROP/SEDDIA50(L,KBT(L)))**0.3)*BDLDTMPP                   ! *** Adjust Phi using Garcia & Parker
          endif
        
          FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP, USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
          SHIELDS = TAUBSNDP/GPDIASEDPROP
          if( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDSC )then
            BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDSC)**SBDLDA(NX)
            ! ***      (-)     (-)                G/M/S     (-)     (-) 
            BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA               ! *** (g/m/s)
          endif

        elseif( ISBDLD(NS) == 2 )then
          ! *** Calculate cell center transport rates using Engelund-Hansen
          if( IBLTAUC(NS) == 2 )then
            CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASEDPROP
          elseif( IBLTAUC(NS) == 3 )then
            CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
          endif
        
          SHIELDS = TAUBSNDP/GPDIASEDPROP
          if( SHIELDS >= CSHIELDS )then
            BDLDTMPP = FSBDLD(DIASEDPROP, GPDIASEDPROP, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
            !IF( HGDH(L) > 0.0) BDLDTMPP = BDLDTMPP/(HGDH(L)**0.333)      Disable grain separation factor
            FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP, USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
            BDLDTMPA = (SBDLDG1(NX)*SHIELDS)**SBDLDA(NX)
            BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA               ! *** (g/m/s)
          endif

        elseif( ISBDLD(NS) == 3 )then
          ! *** Calculate cell center transport rates using Wu, Wang, and Jia
          BDLDTMPP = FSBDLD(DIASEDPROP, GPDIASEDPROP, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
          CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
          if( ISEDEFF == 2 )then
            TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
            CSHIELDS = TMPVAL*CSHIELDS
          endif
        
          FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP, USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
          SHIELDS = TAUBSNDP/GPDIASEDPROP
          if( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDS )then
            BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
            BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA               ! *** (g/m/s)
          endif

        endif  ! *** End of cell centered bedload transport calculations

        ! *** Modify crtical erosion rates to account for cohesive fraction
        if( ISEDEFF == 1 )then
          TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
          BDLDE = TMPVAL*BDLDE
        endif
        EBLD(NS) = BDLDE                                     ! *** Erosion into bedload due to propwash   (g/m/s)
      enddo
    endif
  endif
  
  return

  END Subroutine Calc_Prop_Erosion_Original

!---------------------------------------------------------------------------!
!< @details Add propeller energy to momentum terms in explicit solution
!---------------------------------------------------------------------------!
subroutine add_ship_momentum(FUHJ, FVHJ)
  
  use Variables_Propwash
  use Variables_Ship
  use Mod_Active_Ship    

  implicit none

  ! *** Dummy variables
  Real, intent(inout) :: FUHJ(LCM,KCM),  FVHJ(LCM,KCM)

  ! *** Local variables
  integer :: i, L, k, ip
  real    :: DZPU, PSGZU, PSGZV, UHC, UHC1, UHC2, VHB, VHB1, VHB2
  real    :: ANG, COS1, SIN1, PTOP, PBOT, PAREA, FRAC, RATIO, FRAC_POWER
    
  ip = 1
  do i = 1, total_ships
    if( all_ships(i).efflux_vel > 0.0 )then
      ! *** Add propeller energy
      do
        if( all_ships(i).ship.num_fixed_cells > 0 )then
          ! *** Special Case
          L = all_ships(i).ship.fixed_cells(ip)
          FRAC_POWER = all_ships(i).ship.fixed_frac(ip)
          ip = ip + 1
        else
          L = all_ships(i).pos.cell
          FRAC_POWER = 1.0
        endif
        ang = all_ships(i).heading - 0.5*PI         
        cos1 = cos(ang)
        sin1 = sin(ang)
        
        UHC = all_ships(i).efflux_vel*cos1
        VHB = all_ships(i).efflux_vel*sin1

        ! *** Velocity components rotated to grid
        VHB1 = CVN(L)*UHC + CVE(L)*VHB            
        UHC1 = CUN(L)*UHC + CUE(L)*VHB
      
        ! *** Propeller area and radius already accounts for jet contraction, if not ducted
        PAREA = all_ships(i).prop_area * FRAC_POWER
        
        ! *** Add propeller momentum to ambient conditions.  
        PTOP  = MAX((all_ships(i).draft - all_ships(i).prop_radius), 0.0)
        PBOT  = PTOP + all_ships(i).ship.prop_diam
        DZPU  = 0.0     ! *** Depth to bottom of layer
        RATIO = 0.0
        FRAC = 0.0
      
        ! *** Split the energy into the appropriate layers
        do K = KC,KSZ(L),-1
          DZPU = DZPU + HPK(L,K)
          if( DZPU < PTOP ) CYCLE

          if( FRAC == 0.0 )then
            ! *** Top fraction
            RATIO = (DZPU - PTOP)/all_ships(i).ship.prop_diam
            RATIO = MIN(RATIO,1.0)
            FRAC  = RATIO
          elseif( DZPU > PBOT )then
            ! *** Bottom fraction
            RATIO = 1.0 - FRAC
            FRAC = 1.0
          else
            ! *** Middle fractions
            RATIO = HPK(L,K)/all_ships(i).ship.prop_diam
            RATIO = MIN(RATIO,1.0)
            FRAC  = FRAC + RATIO
          endif
            
          ! *** Add momentum to layer
          PSGZU = 0.0
          if( SGZU(L,K) > 0.0 ) PSGZU = 1.0/SGZU(L,K)
          UHC2  = -RATIO*ABS(UHC1*PAREA)*UHC1
          if( SUB(L) < 0.5 .or. UHC1 > 0.0 )then
            ! *** Apply momentum to east face
            FUHJ(LEC(L),K) = PSGZU*UHC2
          else
            ! *** Apply momentum to west face
            FUHJ(L,K)      = PSGZU*UHC2
          endif
          
          PSGZV = 0.0
          if( SGZV(L,K) > 0.0 ) PSGZV = 1.0/SGZV(L,K)
          VHB2  = -RATIO*ABS(VHB1*PAREA)*VHB1
          if( SVB(L) < 0.5 .or. VHB1 > 0.0 )then
            ! *** Apply momentum to north face
            FVHJ(LNC(L),K) = PSGZV*VHB2
          else
            ! *** Apply momentum to south face
            FVHJ(L,K)      = PSGZV*VHB2
          endif
          if( FRAC >= 1.0 ) EXIT
        enddo    ! *** Layer Loop
        
        if( ip > all_ships(i).ship.num_fixed_cells ) EXIT
      enddo
    endif
  enddo        ! *** Vessel Loop

end subroutine add_ship_momentum
  
  
