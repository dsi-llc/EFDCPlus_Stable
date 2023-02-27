! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

  USE GLOBAL

  IMPLICIT NONE

  ! *** Dummy variables
  Integer, intent(in)           :: L           !> cell index
  Real(kind = rkd), intent(in)  :: TAUP        !> shear stress                  N/m**2
  Real(kind = rkd), intent(out) :: ELAY(NSCM)  !> calculated erosion per class  gm/cm**2
  integer, intent(out) :: SURFACE
  
  ! *** Local variables
  INTEGER :: K, NS, SLLN, NACTLAY, NT, K1, SURFOLD
  INTEGER :: NSC0, NSC1, NTAU0, NTAU1, ICORE

  REAL(RKD) :: CSEDSS,SQR2PI
  REAL(RKD) :: D50TMPP, D50AVGL
  REAL(RKD) :: DEP
  REAL(RKD) :: EBD, ERO
  REAL(RKD) :: ERATEMOD
  REAL(RKD) :: NSCTOT
  REAL(RKD) :: ONE=1.0
  REAL(RKD) :: PFY
  REAL(RKD) :: PX
  REAL(RKD) :: PY
  REAL(RKD) :: SEDFLUX
  REAL(RKD) :: SN00
  REAL(RKD) :: SN01
  REAL(RKD) :: SN10
  REAL(RKD) :: SN11
  REAL(RKD) :: TEMP,TEMP1,TEMP2
  REAL(RKD) :: TACT,TSUM,TAUDYNE
  REAL(RKD) :: TAUCRIT
  REAL(RKD) :: VZDIF
  REAL(RKD) :: WDTDZ
  REAL(RKD) :: ERATEMAX

  REAL(RKD) ,DIMENSION(NSCM) :: QBFLUX
  REAL(RKD) ,DIMENSION(NSCM) :: CSEDVR
  REAL(RKD) ,DIMENSION(NSCM) :: CTB
  REAL(RKD) ,DIMENSION(NSCM) :: DEPBL
  REAL(RKD) ,DIMENSION(NSCM) :: DEPTSS
  REAL(RKD) ,DIMENSION(NSCM) :: ETOT
  REAL(RKD) ,DIMENSION(NSCM) :: PROB
  REAL(RKD) ,DIMENSION(NSCM) :: PROBVR
  REAL(RKD) ,DIMENSION(NSCM) :: SMASS
  REAL(RKD) ,DIMENSION(NSCM) :: TTEMP
  REAL(RKD) ,DIMENSION(KB)   :: INITMASS

  REAL(RKD) ,DIMENSION(2)    :: NSCD
  REAL(RKD) ,DIMENSION(2)    :: TAUDD

  ! *** HARD BOTTOM BYPASS HANDLED IN CALLING ROUTINE

  ! *** ********************************************************************************************************************************************
  ! *** Get things set up for Erosion Calculations
  SURFACE = KB+1
  DO K=2,KB                       ! *** Ignore active layer
    IF( TSED(K,L) > 0.001 )THEN   ! *** Ignore layers than are less than 0.1 mm (typically)
      SURFACE = K
      EXIT
    ENDIF
  ENDDO

  TAUDYNE = 10.*TAUP              ! *** Convert to dynes
  IF( SURFACE > KB )THEN
    ! *** Set a tiny concentration as a flag to be used in SedTran integration if shear stress is high
    WHERE( TAUDYNE >= TAUCRITE(1:NSCM) )
      ELAY(1:NSCM)  = 1.E-12
    ENDWHERE
    RETURN                        ! *** Insufficient sediments
  ENDIF

  ! *** *******************************************************************
  ICORE = NCORENO(IL(L),JL(L))

  ! *** Calculate erosion/deposition for top layer for all cells
  ERO = 0.0                       ! *** Initialize total erosion for the cell
  INITMASS(1:KB) = TSED(1:KB,L)   ! *** Save the starting sediment mass by layers

  ! *** Calculate Average particle size of surface layer so we can calculate
  ! *** active layer unit mass
  D50AVGL = SUM(PERSED(1:NSCM,SURFACE,L)*D50(1:NSCM))             ! *** Calculate local d50 at sediment bed surface

  ! *** Calculate TAUCRIT Based on the Average Particle Size of Surface
  ! *** Then calculate the Active Layer unit mass (TACT) from it.
  ! *** Ta =  Tam * Davg * (Tau/Taucr)
  IF( LAYERACTIVE(SURFACE,L) < 2 )THEN
    ! Identify Size Class interval to use for Taucrit erosion calculation
    IF(D50AVGL < SCND(1) )THEN
      NSCD(1)=SCND(1)
      NSCD(2)=SCND(2)
      NSC0=1
      NSC1=2
      D50AVGL = SCND(1)                                             ! *** Prevent division (s_shear) by zero when there is no sediment in the layer
    ELSEIF( D50AVGL >= SCND(NSICM) )THEN
      NSCD(1)=SCND(NSICM-1)
      NSCD(2)=SCND(NSICM)
      NSC0=NSICM-1
      NSC1=NSICM
    ELSE
      DO NS=1,NSICM-1
        IF( D50AVGL >= SCND(NS) .AND. D50AVGL < SCND(NS+1) )THEN
          NSCD(1)=SCND(NS)
          NSCD(2)=SCND(NS+1)
          NSC0=NS
          NSC1=NS+1
          EXIT
        ENDIF
      ENDDO
    ENDIF

    TAUCRIT = TAUCRITE(NSC0)+(TAUCRITE(NSC1)-TAUCRITE(NSC0))/(NSCD(2)-NSCD(1))*(D50AVGL-NSCD(1))
  ELSE
    ! *** IN-PLACE SEDIMENTS
    TAUCRIT = TAUCOR(SURFACE,L)
  ENDIF

  ! *** Compute the requried active layer thickness (cm)
  IF( TAUDYNE < TAUCRIT )THEN
    RETURN
  ELSE
    TACT = TACTM*D50AVGL*(TAUDYNE/TAUCRIT)*(BULKDENS(1,L)/10000.0)
  ENDIF

  ! *** ********************************************************************************************************************************************
  ! *** Now calculate the Erosion Rate
  K = SURFACE

  ! *** Find upper and lower limits of size classes on mean bed diameter
  IF( (D50AVGL+1E-6) < SCND(1) )THEN
    NS = 1
    NSCD(1) = SCND(NS)
    NSCD(2) = SCND(NS+1)
    NSC0 = NS
    NSC1 = NS+1
    D50AVGL = SCND(1)
  ELSEIF( (D50AVGL-1E-6) > SCND(NSICM) )THEN
    NS = NSICM - 1
    NSCD(1) = SCND(NS)
    NSCD(2) = SCND(NS+1)
    NSC0 = NS
    NSC1 = NS+1
    D50AVGL = SCND(NSICM)
  ELSE
    DO NS=1,NSICM-1
      IF( D50AVGL >= SCND(NS) .AND. D50AVGL < SCND(NS+1) )THEN
        NSCD(1) = SCND(NS)
        NSCD(2) = SCND(NS+1)
        NSC0 = NS
        NSC1 = NS+1
        EXIT
      ENDIF
    ENDDO
  ENDIF

  ! *** Calculate TAUCRIT Based on the D50 of the bed or from Sedflume Data
  IF( LAYERACTIVE(SURFACE,L) < 2 )THEN                ! *** For active/deposited layers
    TAUCRIT = TAUCRITE(NSC0) + (TAUCRITE(NSC1)-TAUCRITE(NSC0))/(NSCD(2)-NSCD(1))*(D50AVGL-NSCD(1)) !interpolation
    TAUCOR(K,L) = TAUCRIT
  ELSE
    ! *** SEDFlume data (depth interpolation)
    SN01 = TSED(K,L)/TSED0(K,L)                       ! *** Weighting factor 1 for interpolation
    SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)          ! *** Weighting factor 2
    TAUCRIT = SN01*TAUCOR(K,L) + SN11*TAUCOR(K+1,L)
  ENDIF

  ! *** Check if the shear is greater than critical shears.  If not, exit erosion FUNCTION
  IF( TAUDYNE < TAUCRIT ) RETURN

  ! *** Now, calculate erosion rates

  ! *** Find the upper and lower limits of the Shear Stress for the interpolation
  IF( NSEDFLUME == 1 )THEN

    IF( TAUDYNE >= TAULOC(ITBM) )THEN
      TAUDD(1) = TAULOC(ITBM-1)
      TAUDD(2) = TAULOC(ITBM)
      NTAU0 = ITBM-1
      NTAU1 = ITBM

    ELSEIF( TAUDYNE < TAULOC(1) )THEN
      TAUDD(1) = TAULOC(1)
      TAUDD(2) = TAULOC(2)
      NTAU0 = 1
      NTAU1 = 2
    ELSE
      DO NS=1,ITBM-1
        IF( TAUDYNE >= TAULOC(NS) .AND. TAUDYNE < TAULOC(NS+1) )THEN
          TAUDD(1) = TAULOC(NS)
          TAUDD(2) = TAULOC(NS+1)
          NTAU0 = NS
          NTAU1 = NS+1
          EXIT
        ENDIF
      ENDDO
    ENDIF

  ENDIF

  ! *** Interpolate the erosion rates for shear stress and depth.
  ! *** This utilizes normal sedflume data for deeper layers.
  IF( LAYERACTIVE(SURFACE, L) == 2 )THEN
    ! *** Calculate erosion rates of deeper layers (SEDFlume data)
    IF( NSEDFLUME == 1 )THEN
      SN00 = (TAUDD(2)-TAUDYNE)/(TAUDD(2)-TAUDD(1)) ! *** weighting factor 1 for interpolation
      SN10 = (TAUDD(1)-TAUDYNE)/(TAUDD(1)-TAUDD(2)) ! *** weighting factor 2
      SN01 = TSED(K,L)/TSED0(K,L)                   ! *** weighting factor 3
      SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)      ! *** weighting factor 4

      IF( K+1 <= KB )THEN  ! *** Maximum erosion rate
        ERATEMAX=( SN00*EXP(SN11*LOG(ERATE(K+1,L,NTAU0))+SN01*LOG(ERATE(K,L,NTAU0))) &
          + SN10*EXP(SN11*LOG(ERATE(K+1,L,NTAU1))+SN01*LOG(ERATE(K,L,NTAU1))) )*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))
      ELSE                 ! *** Do not allow erosion through the bottom layer
        ERATEMAX = ( SN00*EXP(SN11*LOG(1.0E-9)+SN01*LOG(ERATE(K,L,NTAU0))) &
          + SN10*EXP(SN11*LOG(1.0E-9)+SN01*LOG(ERATE(K,L,NTAU1))) )*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))
      ENDIF
    ELSE
      IF( TAUDYNE > TAUCOR(K,L) )THEN                                                  ! *** Check that the applied shear exceeds the critical shear stress for this layer
        ! *** Erosion rate values (cm/s) computed by equation assume shear in Pascals so convert dynes
        SN00 = EA(ICORE,K)*TAUP**EN(ICORE,K)                                           ! *** Erosion rate (cm/s) of the top layer

        IF( K+1 <= KB )THEN
          SN10 = EA(ICORE,K+1)*TAUP**EN(ICORE,K+1)                                     ! *** Erosion rate (cm/s) of the layer below
        ELSE
          SN10 = 0.0                                                                   ! *** Modeled erosion rate in limited by bottom
        ENDIF

        SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)                                       ! *** Mass weighting factor
        ERATEMAX = ((SN10-SN00)*SN11 + SN00)*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))       ! *** linear interpolation for remaining mass in current layer    (g/cm2/s)
        ERATEMAX = MIN(ERATEMAX,MAXRATE(ICORE,K))                                      ! *** Limit erosion rate
      ELSE
        ERATEMAX = 0.0
      ENDIF
    ENDIF
  ELSE
    ! *** For Active and deposited sediment layers
    ! *** The erosion rate for these layers is determined from
    ! *** Sedflume experiments and is based on average particle Size (D50AVG)
    ! ***
    NSCTOT = NSCD(2)-NSCD(1)                                                         ! *** difference in interpolant size class
    D50TMPP = D50AVGL-NSCD(1)                                                      ! *** difference from local size class and lower interpolant
    IF( NSEDFLUME == 1 )THEN
      SN00 = (TAUDD(2)-TAUDYNE)/(TAUDD(2)-TAUDD(1))                                  ! *** weighting factor 1 for interpolation
      SN10 = (TAUDD(1)-TAUDYNE)/(TAUDD(1)-TAUDD(2))                                  ! *** weigthing factor 2
      SN01 = D50TMPP/NSCTOT                                                          ! *** weighting factor 3
      SN11 = (NSCTOT-D50TMPP)/NSCTOT                                                 ! *** weighting factor 4
      ERATEMAX = (SN00*EXP(SN11*LOG(ERATEND(NSC0,NTAU0)) + SN01*LOG(ERATEND(NSC1,NTAU0))) + SN10*EXP(SN11*LOG(ERATEND(NSC0,NTAU1)) +  &   ! *** log-linear interpolation
        SN01*LOG(ERATEND(NSC1,NTAU1))))*BULKDENS(K,L)*SQRT(1./SH_SCALE(L))
    ELSE
      ! *** Erosion rate values (cm/s) computed by equation assume shear in Pascals so convert dynes
      SN00 = ACTDEPA(NSC0)*TAUP**ACTDEPN(NSC0)                                       ! *** Erosion rate 1 (cm/s)
      SN10 = ACTDEPA(NSC1)*TAUP**ACTDEPN(NSC1)                                       ! *** Erosion rate 2 (cm/s)
      SN11 = D50TMPP/NSCTOT                                                          ! *** Weighting factor
      ERATEMAX = ((SN10-SN00)*SN11 + SN00)*BULKDENS(K,L)*SQRT(1./SH_SCALE(L))        ! *** linear interpolation around size class (g/cm2/s)
      ERATEMAX = MIN(ERATEMAX,ACTDEPMAX(NSC0))                                       ! *** Limit erosion rate
    ENDIF
  ENDIF

  ! *** Sort out Thicknesses and Erosion Rates
  EBD = ERATEMAX*DTSEDJ                                                              ! *** Maximum mass potentially eroded this time step for this layer (g/cm^2)

  ! *** If the shear stress is less than the critical shear stress for a
  ! *** particular size class, then it is not eroded from the bed.

  ! *** Conservation of sediment mass will be addressed on a full cell basis in SEDZLJ
  ! *** ELAY(NS) = Total erosion at this cell of size class NS
  ! *** ERO      = Total erosion at this cell
  DO NS = 1,NSCM
    IF( TAUDYNE >= TCRE(ns) )THEN
      ELAY(NS)  = PERSED(NS,K,L)*EBD
    ELSE
      ELAY(NS)  = 0.0
    ENDIF
  ENDDO
  
  !ERO = SUM(ELAY(1:NSCM))                                                            ! *** Total erosion from the layer   (g/cm^2)

  RETURN

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
  USE GLOBAL

  IMPLICIT NONE

  ! *** Dummy variables
  Integer, intent(in)           :: L           !> cell index
  Real(kind = rkd), intent(in)  :: TAUP        !> density normalized shear stress               m2/s2
  Real(kind = rkd), intent(out) :: ELAY(NSCM)  !> calculated erosion into suspension per class  g/m2
  Real(kind = rkd), intent(out) :: EBLD(NSCM)  !> calculated erosion into bedload per class     g/m2
  
  Integer :: NS, NX
  Real :: CSHIELDS, FACBEDLP, FACSUSLP, ROUSEP, SHIELDS, SNDEQP, SIGP, SNDEQBP, TAUE, TAUBSNDP, TAURTMP
  Real :: TMPSEDHID, USTARP, USTARSNDP, WESE, WESEMX, WSFAC, ZEQP, ZEQDP, ZEQDIP, ZEQMAX
  Real :: BDLDE, BDLDTMPP, BDLDTMP, BDLDTMPA, BDLDTMPB, CSHIELDSC
  Real :: BOT, TOP, TMPSTR, TMPVAL, DIASEDPROP, DSEDGMMPROP, GPDIASEDPROP
  REAL, EXTERNAL :: FSBDLD, FSEDMODE, CSNDZEQ, CSNDEQC
  
  ! *** Cohesive sediments
  if( ISTRAN(6) > 0 )then
    IF( IWRSP(1) == 4 )THEN
      TMPSTR=0.0
    ELSE
      TMPSTR=1.0
    ENDIF    

    DO NS = 1,NSED
      ! *** SET MAXIMUM EROSION RATE
      WESEMX = DELTI*SEDB(L,KBT(L),NS)

      IF( TAUP > TAURB(L,KBT(L)) )THEN
        ! *** MASS/BULK EROSION
        WESE = WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
        WESE = MIN(WESE,WESEMX)
      ELSE
        TAUE = 0.
        IF( TAUP > TAURS(L,KBT(L)) )THEN
          ! *** SURFACE EROSION
          WESE = WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)

          ! *** SET NORMALIZING SHEAR BASED ON INPUT OPTION
          IF( IWRSP(1) >= 99 )THEN
            ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU FROM SSCOHSEDPMAP
            TAURTMP = TAUNS(L,KBT(L))
          ELSEIF ((IWRSP(1) >= 2) .AND. (IWRSP(1) < 4) )THEN
            ! *** 1 HWANG AND METHA - LAKE OKEECHOBEE
            ! *** 2 HAMRICK'S MODIFICATION OF SANFORD AND MAA
            ! *** 3 SAME AS 2 EXCEPT VOID RATIO OF COHESIVE SEDIMENT FRACTION IS USED
            TAURTMP = TAUR(1)
          ELSE
            ! *** USE DIRECTLY INPUT PROPERTIES
            TAURTMP = TAURS(L,KBT(L))
          ENDIF

          TAUE = (TAUP-TMPSTR*TAURS(L,KBT(L)))/TAURTMP
          TAUE = MAX(TAUE,0.0)

          ! *** SET NON-COHESIVE HIDING FACTOR
          TMPSEDHID = 1.0
          IF( ISTRAN(7) >= 1 .AND. COSEDHID(1) /= 0.0 )THEN
            TMPSEDHID = (FRACCOH(L,KBT(L)))**COSEDHID(1)
          ENDIF
              
          ! *** FINALIZE SURFACE EROSION RATE
          IF( IWRSP(1) < 99 )THEN
            WESE = TMPSEDHID*WESE*( TAUE**TEXP(NS) )
          ELSE
            ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU EXPONENT FROM SSCOHSEDPMAP
            WESE = TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
          ENDIF
          WESE = MIN(WESE,WESEMX)
        ELSE
          ! *** NO EROSION 
          WESE = 0.0
        ENDIF
        ELAY(NS) = WESE                                            ! *** Erosion flux in g/m2/s
      ENDIF
    ENDDO
  endif
  
  ! *** Non-Cohesive sediments
  if( ISTRAN(7) > 0 )then
    
    ! *** Erosion flux into the water column
    TAUBSNDP = TAUP                                                ! *** Propwash shear.  Use total shear m2/s2
    USTARP = SQRT(TAUP)                                            ! *** USTAR m/s
    USTARSNDP = USTARP
    
    DO NX = 1,NSND
      NS = NSED + NX
      
      DSEDGMMPROP  = 1./(1.E6*SSG(NS))                             ! *** Specific volume (m**3/g)
      DIASEDPROP   = SEDDIA(NS)                                    ! *** Sediment particle diameter (m)
      GPDIASEDPROP = G*(SSG(NS)-1.)*DIASEDPROP  
      
      ! *** SET THE ROUSE PARAMETER USING TOTAL
      IF( USTARP > 1E-12 )THEN  
        ROUSEP = WSETA(L,0,NS)/(VKC*USTARP)  
      ELSE  
        ROUSEP = 250000.*WSETA(L,0,NS)  
      END IF  
      
      ! *** FSEDMODE SETS BEDLOAD (IMODE=1) AND SUSPENDED LOAD (IMODE=2) TRANSPORT FRACTIONS
      FACSUSLP = FSEDMODE(WSETA(L,0,NS),USTARP,USTARSNDP,RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)  

      WESE = 0.
      IF( ISNDEQ(NS) == 0 )THEN
        ! *** USER SPECIFIED EQUILIBRIUM CONCENTRATION
        ZEQP    = 0.05         ! *** Garcia & Parker
        ZEQDP   = ZEQP*DZIC(L,KSZ(L))  
        ZEQDIP  = 1./ZEQDP  
        SNDEQBP = TAUR(NS)
      ELSE
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
        IF( ISNDEQ(NS) == 1 )THEN  
          IF( ISLTAUC(NS) <= 1 )THEN  
            CSHIELDS = TCSHIELDS(NS)  
          ELSEIF( ISLTAUC(NS) == 2 )THEN  
            CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASEDPROP  
          ELSEIF( ISLTAUC(NS) == 3 )THEN  
            CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)  
          ELSE
            CSHIELDS = TCSHIELDS(NS)  
          ENDIF  

          ! *** Modify crtical shear stress to account for cohesive fraction
          IF( ISEDEFF == 2 )THEN  
            TMPVAL = 1. + (COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
            CSHIELDS = TMPVAL*CSHIELDS  
          ENDIF  
          SHIELDS = TAUBSNDP/GPDIASEDPROP  
          IF( SHIELDS < CSHIELDS ) SNDEQBP = 0.  
        ENDIF  

        ! *** Modify crtical erosion rates to account for cohesive fraction
        IF( ISEDEFF == 1 )THEN
          SNDEQBP = SNDEQBP*EXP(-COEHEFF*FRACCOH(L,KBT(L)))  
        ENDIF
      ENDIF
        
      IF( ROUSEP < 0.999 .OR. ROUSEP > 1.001 )THEN  
        TOP = (ZEQDP**(ROUSEP-1.))-1.  
        BOT = (1.-ROUSEP)*(ZEQDIP-1.)  
        SNDEQP = SNDEQBP*TOP/BOT  
        SNDEQP = FACSUSLP*VFRBED(L,KBT(L),NS)*MAX(SNDEQP,0.)  
      ELSE  
        TOP = LOG(ZEQDIP)  
        BOT = (ZEQDIP-1.)  
        SNDEQP = SNDEQBP*TOP/BOT  
        SNDEQP = FACSUSLP*VFRBED(L,KBT(L),NS)*MAX(SNDEQP,0.)  
      ENDIF  

      ! *** SET RESUSPENSION FLUX
      IF( KSZ(L) == KC )THEN
        WSFAC = 2.*(1.+ROUSEP)/(2.+ROUSEP*(1.-ZEQP))               ! *** Single layer
      ELSE
        WSFAC = 1.
      ENDIF
      WESE = WSFAC*WSETA(L,0,NS)*SNDEQP                            ! *** Erosion flux in g/m2/s  
      ELAY(NS) = WESE*DTSEDJ                                       ! *** Erosion flux in g/m2
    ENDDO
    
    ! *** Erosion flux into bedload
    IF( ICALC_BL > 0 )THEN
      DO NX = 1,NSND
        NS = NSED + NX
      
        DSEDGMMPROP  = 1./(1.E6*SSG(NS))                           ! *** Specific volume (m**3/g)
        DIASEDPROP   = SEDDIA(NS)                                  ! *** Sediment particle diameter (m)
        GPDIASEDPROP = G*(SSG(NS)-1.)*DIASEDPROP                   ! *** "excess" density (m2/s2)
        BDLDTMP = SQRT(GPDIASEDPROP)*DIASEDPROP/DSEDGMMPROP

        CSHIELDS = TCSHIELDS(NS)
      
        ! *** Modify crtical shear stress to account for cohesive fraction
        IF( ISEDEFF == 2 )THEN
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        ENDIF
      
        !*** Compute the cell centered bedload transport rates using the specified option
        BDLDE = 0.0
        IF( ISBDLD(NS) == 0 )THEN
          ! *** Calculate cell center transport rates using generic bed load equation
          BDLDTMPP = SBDLDP(NX)                                                          ! *** Constant Phi
          IF( BDLDTMPP > 0.0 )THEN
            FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP ,USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
            SHIELDS = TAUBSNDP/GPDIASEDPROP
            IF( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDS )THEN
              IF( SBDLDA(NX) > 0.0 )THEN
                BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
              ELSE
                BDLDTMPA = 1.0
              ENDIF
              IF( SBDLDB(NX) > 0.0 )THEN
                BDLDTMPB = (SBDLDG3(NX)*SQRT(SHIELDS) - SBDLDG4(NX)*SQRT(CSHIELDS))**SBDLDB(NX)
              ELSE
                BDLDTMPB = 1.0
              ENDIF
              BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA*BDLDTMPB    ! *** (g/m/s)
            ENDIF
          ENDIF

        ELSEIF( ISBDLD(NS) == 1 )THEN
          ! *** Calculate cell center transport rates using Van Rijn bed load equation

          ! *** Activate Garcia & Parker armoring
          IF( ISNDAL == 1 )THEN
            TMPVAL = LOG10(19.*DIASEDPROP/SEDDIA50(L,KBT(L)))
            TMPVAL = 1.66667/(TMPVAL**2)
            CSHIELDSC = CSHIELDS50(L)*TMPVAL
          ELSE
            CSHIELDSC = TCSHIELDS(NS)
          ENDIF

          BDLDTMPP = FSBDLD(DIASEDPROP, GPDIASEDPROP, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))   ! *** Calculate Phi
          IF( ISNDAL == 1 )THEN
            BDLDTMPP = ((DIASEDPROP/SEDDIA50(L,KBT(L)))**0.3)*BDLDTMPP                   ! *** Adjust Phi using Garcia & Parker
          ENDIF
        
          FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP, USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
          SHIELDS = TAUBSNDP/GPDIASEDPROP
          IF( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDSC )THEN
            BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDSC)**SBDLDA(NX)
            ! ***      (-)     (-)                G/M/S     (-)     (-) 
            BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA               ! *** (g/m/s)
          ENDIF

        ELSEIF( ISBDLD(NS) == 2 )THEN
          ! *** Calculate cell center transport rates using Engelund-Hansen
          IF( IBLTAUC(NS) == 2 )THEN
            CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASEDPROP
          ELSEIF( IBLTAUC(NS) == 3 )THEN
            CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
          ENDIF
        
          SHIELDS = TAUBSNDP/GPDIASEDPROP
          IF( SHIELDS >= CSHIELDS )THEN
            BDLDTMPP = FSBDLD(DIASEDPROP, GPDIASEDPROP, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
            !IF( HGDH(L) > 0.0) BDLDTMPP = BDLDTMPP/(HGDH(L)**0.333)      Disable grain separation factor
            FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP, USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
            BDLDTMPA = (SBDLDG1(NX)*SHIELDS)**SBDLDA(NX)
            BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA               ! *** (g/m/s)
          ENDIF

        ELSEIF( ISBDLD(NS) == 3 )THEN
          ! *** Calculate cell center transport rates using Wu, Wang, and Jia
          BDLDTMPP = FSBDLD(DIASEDPROP, GPDIASEDPROP, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
          CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
          IF( ISEDEFF == 2 )THEN
            TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
            CSHIELDS = TMPVAL*CSHIELDS
          ENDIF
        
          FACBEDLP = FSEDMODE(WSETA(L,0,NS), USTARP, USTARSNDP, RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
          SHIELDS = TAUBSNDP/GPDIASEDPROP
          IF( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDS )THEN
            BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
            BDLDE = FACBEDLP*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA               ! *** (g/m/s)
          ENDIF

        ENDIF  ! *** End of cell centered bedload transport calculations

        ! *** Modify crtical erosion rates to account for cohesive fraction
        IF( ISEDEFF == 1 )THEN
          TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
          BDLDE = TMPVAL*BDLDE
        ENDIF
        EBLD(NS) = BDLDE                                     ! *** Erosion into bedload due to propwash   (g/m/s)
      ENDDO
    ENDIF
  endif
  
  RETURN

  END Subroutine Calc_Prop_Erosion_Original

!---------------------------------------------------------------------------!
!< @details Add propeller energy to momentum terms in explicit solution
!---------------------------------------------------------------------------!
subroutine add_ship_momentum(FUHJ, FVHJ)
  
  Use Variables_Propwash
  Use Variables_Ship
  Use Mod_Active_Ship    

  implicit none

  ! *** Dummy variables
  Real, intent(inout) :: FUHJ(LCM,KCM),  FVHJ(LCM,KCM)

  ! *** Local variables
  integer :: i, L, k, ip
  real    :: DZPU, UHC, UHC1, UHC2, VHB, VHB1, VHB2
  real    :: ANG, COS1, SIN1, PTOP, PBOT, PAREA, FRAC, RATIO, FRAC_POWER
    
  ip = 1
  DO i = 1, total_ships
    IF( all_ships(i).efflux_vel > 0.0 )THEN
      ! *** Add propeller energy
      do
        IF( all_ships(i).ship.num_fixed_cells > 0 )THEN
          ! *** Special Case
          L = all_ships(i).ship.fixed_cells(ip)
          FRAC_POWER = all_ships(i).ship.fixed_frac(ip)
          ip = ip + 1
        ELSE
          L = all_ships(i).pos.cell
          FRAC_POWER = 1.0
        ENDIF
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
        DO K = KC,KSZ(L),-1
          DZPU = DZPU + HPK(L,K)
          IF( DZPU < PTOP ) CYCLE

          IF( FRAC == 0.0 )THEN
            ! *** Top fraction
            RATIO = (DZPU - PTOP)/all_ships(i).ship.prop_diam
            RATIO = MIN(RATIO,1.0)
            FRAC  = RATIO
          ELSEIF( DZPU > PBOT )THEN
            ! *** Bottom fraction
            RATIO = 1.0 - FRAC
            FRAC = 1.0
          ELSE
            ! *** Middle fractions
            RATIO = HPK(L,K)/all_ships(i).ship.prop_diam
            RATIO = MIN(RATIO,1.0)
            FRAC  = FRAC + RATIO
          ENDIF
            
          ! *** Add momentum to layer
          UHC2  = -RATIO*ABS(UHC1*PAREA)*UHC1
          IF( SUB(L) < 0.5 .OR. UHC1 > 0.0 )THEN
            ! *** Apply momentum to east face
            FUHJ(LEC(L),K) = FSGZU(L,K)*UHC2
          ELSE
            ! *** Apply momentum to west face
            FUHJ(L,K)      = FSGZU(L,K)*UHC2
          ENDIF
          
          VHB2  = -RATIO*ABS(VHB1*PAREA)*VHB1
          IF( SVB(L) < 0.5 .OR. VHB1 > 0.0 )THEN
            ! *** Apply momentum to north face
            FVHJ(LNC(L),K) = FSGZV(L,K)*VHB2
          ELSE
            ! *** Apply momentum to south face
            FVHJ(L,K)      = FSGZV(L,K)*VHB2
          ENDIF
          IF( FRAC >= 1.0 ) EXIT
        ENDDO    ! *** Layer Loop
        
        IF( ip > all_ships(i).ship.num_fixed_cells ) EXIT
      enddo
    ENDIF
  ENDDO        ! *** Vessel Loop

end subroutine add_ship_momentum
  
  
