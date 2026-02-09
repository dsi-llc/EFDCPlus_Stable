! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEDZLJ(L)
  
  ! CALCULATES DEPOSITION, ENTRAINMENT, NET FLUX, TOTAL THICKNESS,
  ! LAYER ORDER, AND COMPONENT LAYER FRACTION
  ! UCSB, Craig Jones and Wilbert Lick

  ! ORIGINAL:  May 24, 2006
  ! REVISION DATE :  August 29th, 2007
  !  Craig Jones and Scott James
  !  Updated to fix Active Layer issues

  ! 2016-11-07  Paul M. Craig  Updated for SIGMA-ZED AND OMP
  ! 2017-01-04  Paul M. Craig  Toxics and bedload mass balance updates
  ! 2021-06-07  Paul M. Craig  Updated for propwash
  ! 2022-01-07  Paul M. Craig  Added mass erosion classes for propwash

  use GLOBAL
  use FIELDS
  use Variables_Propwash
  
  implicit none
  
  integer :: K, NS, L, SURFACE, SLLN, NACTLAY, NT, K1, SURFOLD
  integer :: NSC0, NSC1, NTAU0, NTAU1, ICORE

  real(RKD) :: CSEDSS,SQR2PI
  real(RKD) :: D50TMPP
  real(RKD) :: DEP
  real(RKD) :: EBD
  real(RKD) :: ERATEMOD
  real(RKD) :: ERO
  real(RKD) :: NSCTOT
  real(RKD) :: ONE = 1.0_8
  real(RKD) :: PFY
  real(RKD) :: PX
  real(RKD) :: PY
  real(RKD) :: RATIOMASS
  real(RKD) :: SEDFLUX
  real(RKD) :: SN00
  real(RKD) :: SN01
  real(RKD) :: SN10
  real(RKD) :: SN11
  real(rkd) :: SUMPER
  real(RKD) :: TEMP, TEMP1, TEMP2
  real(RKD) :: TACT, TSUM
  real(RKD) :: TAUCRIT, TCRSUS_FAST
  real(RKD) :: VZDIF
  real(RKD) :: WDTDZ

  real(RKD) ,dimension(NSEDS)  :: CSEDVR
  real(RKD) ,dimension(NSEDS2) :: CTB
  real(RKD) ,dimension(NSEDS2) :: DEPBL
  real(RKD) ,dimension(NSEDS2) :: DEPTSS
  real(RKD) ,dimension(NSEDS2) :: DEPTSSB
  real(RKD) ,dimension(NSEDS2) :: ELAY
  real(RKD) ,dimension(NSEDS2) :: ETOT
  real(RKD) ,dimension(NSEDS2) :: PROB
  real(RKD) ,dimension(NSEDS)  :: PROBVR
  real(RKD) ,dimension(NSEDS2) :: QBFLUX
  real(RKD) ,dimension(NSEDS2) :: SMASS
  real(RKD) ,dimension(NSEDS)  :: TTEMP
  real(RKD) ,dimension(KBM,NSEDS)  :: PMC   ! DELME

  real(RKD) ,dimension(KB) :: INITMASS
  real(RKD) ,dimension(2)  :: NSCD
  real(RKD) ,dimension(2)  :: TAUDD
  
  logical :: NOPROPWASH
  
  ! *** *******************************************************************
  ICORE = NCORENO(IL(L),JL(L))

  ! *** CALCULATE EROSION/DEPOSITION FOR TOP LAYER FOR ALL CELLS
  CSEDVR(1:NSEDS)    = 0.0    ! *** Initialize van Rijn's equilibrium bedload concentration
  CTB(1:NSEDS2)      = 0.0    ! *** Initialize bottom layer concentration 
  DEPBL(1:NSEDS2)    = 0.0    ! *** Initialize deposition from bedload
  DEPTSS(1:NSEDS2)   = 0.0    ! *** Initialize deposition from suspended load
  DEPTSSB(1:NSEDS2)  = 0.0    ! *** Initialize deposition from suspended load collapsed into the bed class (propwash)
  ELAY(1:NSEDS)      = 0.0    ! *** Initialize top-layer erosion rates for each sediment size class
  ETOT(1:NSEDS2)     = 0.0    ! *** Initialize erosion rates for each size class and each cell
  PROB(1:NSEDS2)     = 0.0    ! *** Initialize probability of deposition of suspended load
  PROBVR(1:NSEDS)    = 0.0    ! *** 
  QBFLUX(1:NSEDS2)   = 0.0    ! *** Initialize suspended load flux
  SMASS(1:NSEDS2)    = 0.0    ! *** Initialize bedload sediment mass available for deposition  (g/cm^2)
  TTEMP(1:NSEDS)     = 0.0    ! *** 
  
  QBLFLUX(L,1:NSEDS) = 0.0    ! *** Initialize bedload flux
  
  DEP = 0.0                  ! *** Initialize total deposition for the cell
  EBD = 0.0                  ! *** Initialize non-propwash total erosion for the cell
  NOPROPWASH = .TRUE.
  
  INITMASS(1:KB) = TSED(1:KB,L)   ! *** Save the starting sediment mass by layers  (g/cm^2)
  DO K = 1,KB
    PMC(K,1:NSEDS) = PERSED(1:NSEDS,K,L)*INITMASS(K)   ! DELME
  ENDDO
  SQR2PI = 1. / SQRT(2.*PI)

  ! *** Convert Bottom Shear from (m/s)^2 to dynes/cm^2, if using shear from EFDC
  !     TAU(L) = TAUBSED(L)*1000.0*10.0
  
  ! *** Convert Bottom Concentrations and estimate bottom concentration for KC = 1
  if( KSZ(L) == KC )then    ! *** Alberta
    ! *** Estimate bottom concentration assuming exponential sediment concentration profile in single layer
    USTAR(L) = SQRT(TAU(L)/10000.0)                                     ! *** USTAR (m/s) = SQRT(Tau/RhoH2O) and RhoH2O is 1000 kg/cm^3. 
    do NS = 1,NSEDS2
      VZDIF = max(20.0,0.067*HPCM(L)*USTAR(L)*100._8)                   ! *** Convert units of USTAR to cm/s
      TEMP2 = HPCM(L)*DWS(NS)/VZDIF
      CTB(NS) = SED(L,KSZ(L),NS)*TEMP2*(ONE/(ONE-EXP(-TEMP2)))*1.0D-06  ! *** Convert bottom layer concentration from mg/L (g/m^3) to g/cm^3
    enddo
  else
    CTB(1:NSEDS2) = SED(L,KSZ(L),1:NSEDS2)*1.0D-06                        ! *** Convert bottom layer concentration from mg/L (g/m^3) to g/cm^3
  endif
  
  ! *** ********************************************************************************************************************************************
  ! CALCULATE DEPOSITION
  ! *** Temporarily calculate the sediment mass in the active layer (g/cm^3)
  ! *** so that percentages can be determined after deposition.
  
  ! *** Deposition from suspended load by sediment sizes
  do NS = 1,NSEDS
    ! *** If there is no sediment of that size available in the water column then it cannot be deposited
    if( CTB(NS) < 1.0E-20 ) CYCLE                                      

    ! *** Calculate probability for suspended load deposition
    ! *** Based on Gessler (1965) if D50 > 200um or Krone if < 200 um  (SNL Eqs 15 - 18)
    if( PROP_ERO(L,0) > 0. )then
      PROB(NS) = 0.0                                                ! *** Zero probability when propeller wash is active
      NOPROPWASH = .FALSE.
    elseif( D50(NS) >= BEDLOAD_CUTOFF )then
      PY  = 1.7544*(TCRSUS(NS)/(TAU(L) + 1.0D-18)-ONE)              ! *** 1.7544 = 1/sigma =  1./0.57 from Gessler
      if( PY >= 0.0 )then
        PFY = SQR2PI*EXP(-0.5d0*PY*PY)
        PX  = ONE/(ONE + 0.33267*PY)
        PROB(NS) = ONE - PFY*(0.43618*PX - 0.12016*PX*PX + 0.93729*PX**3)
      else
        PY  = DABS(PY)
        PFY = SQR2PI*EXP(-0.5d0*PY*PY)
        PX  = ONE/(ONE + 0.33267*PY)
        PROB(NS) = PFY*(0.43618*PX - 0.12016*PX*PX + 0.93729*PX**3)
      endif
    elseif( TAU(L) <= TCRSUS(NS) )then
      PROB(NS) = ONE - TAU(L)/(TCRSUS(NS))                          ! *** Krones deposition probability is calculated
    else
      PROB(NS) = 0.0                                                ! *** Zero probability when TAU > TCRSUS
    endif
     
    ! *** Calculate SMASS of sediment present in water and allow only that much to be deposited.
    SMASS(NS)  = SED(L,KSZ(L),NS)*1.0D-06*DZC(L,KSZ(L))*HPCM(L)*MAXDEPLIMIT         ! *** SMASS(NS) is the total sediment mass in the first layer.  It is calculated as a precaution so that no more than the total amount of mass in the first layer can deposit onto the sediment bed PERSED time step.
    DEPTSS(NS) = CTB(NS)*PROB(NS)*(DWS(NS)*DTSEDJ)                                  ! *** Deposition of a size class is equal to the probability of deposition times the settling rate times the time step times the sediment concentration
    DEPTSS(NS) = min(MAX(DEPTSS(NS),0.0),SMASS(NS))                                 ! *** Do not allow more sediment to deposit than is available in the water-column layer above the bed
    DEPTSSB(NS) = DEPTSS(NS)                                                        ! *** Non "fast" classes go straight to the bed
    DEP_SED_FLX(L,NS) = DEP_SED_FLX(L,NS) + DEPTSSB(NS)                             ! *** Accumulate total deposition for each origin class (g/cm^2)
  enddo
  
  ! *** Collapse propwash "fast" classes into the associated standard class
  if( NOPROPWASH )then
    if( NSEDS2 > NSEDS )then
      do NS = NSEDS+1, NSEDS2
        TCRSUS_FAST = TCRSUS(IWC2BED(NS))*fast_multiplier
        if( TAU(L) <= TCRSUS_FAST )then
          PROB(NS) = ONE - TAU(L)/TCRSUS_FAST                                         ! *** Krones deposition probability is calculated
          ! *** Calculate SMASS of sediment present in water and allow only that much to be deposited.
          SMASS(NS)  = SED(L,KSZ(L),NS)*1.0D-06*DZC(L,KSZ(L))*HPCM(L)*MAXDEPLIMIT     ! *** SMASS(NS) is the total sediment mass in the first layer.  It is calculated as a precaution so that no more than the total amount of mass in the first layer can deposit onto the sediment bed PERSED time step.
          DEPTSS(NS) = CTB(NS)*PROB(NS)*(DWS(NS)*DTSEDJ)                              ! *** Deposition of a size class is equal to the probability of deposition times the settling rate times the time step times the sediment concentration
          DEPTSS(NS) = min(MAX(DEPTSS(NS),0.0),SMASS(NS))                             ! *** Do not allow more sediment to deposit than is available in the water-column layer above the bed
          DEPTSSB(IWC2BED(NS))  = DEPTSSB(IWC2BED(NS)) + DEPTSS(NS)                   ! *** Collapsed "fast" classes to origin class
          DEP_SED_FLX(L,NS) = DEP_SED_FLX(L,NS) + DEPTSS(NS)                          ! *** Accumulate total deposition for each fast class (g/cm^2)   
        endif
      enddo
    endif
  endif
  
  
  ! *** Deposition from Bedload
  if( ICALC_BL > 0 .and. PROP_ERO(L,0) == 0.0 )then
    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE

      if( CBL(L,NS) < 1E-20 ) CYCLE
      SMASS(NS) = CBL(L,NS)                                   ! *** Local mass available for deposition is the bedload concentration times the bedload height
      CSEDSS = SMASS(NS)/(DWS(NS)*DTSEDJ)                     ! *** Concentration eroded into bedload from the last time step
      CSEDVR(NS) = (0.18*2.65*TRANS(L,NS)/DISTAR(NS)*0.65)    ! *** van Rijn's (1981, Eq. 21) equilibrium bedload concentration (here TRANS is the transport parameter for bedload calculations)
      if( CSEDVR(NS) <= 0.0 )then                             ! *** If there is no equilibrium bedload available
        PROBVR(NS) = 1.0                                      ! *** then deposition probability is unity
      else
        PROBVR(NS) = min(CSEDSS/CSEDVR(NS),1.0)               ! *** van Rijn probability of deposition from bedload
      endif
      if( CSEDSS <= 0.0 ) PROBVR(NS) = PROB(NS)               ! *** In case CBL = 0  The deposition from bedload reverts to Gessler's for that particle size.
      
      DEPBL(NS) = PROBVR(NS)*CBL(L,NS)*(DWS(NS)*DTSEDJ)       ! *** Calculate Bedload Deposition (DEPBL)  (g/cm^2)
                                                              ! *** Deposition from bedload is the van Rijn probability times bedload concentration time settling velocity times the time step
      DEPBL(NS) = min(MAX(DEPBL(NS),0.0),SMASS(NS))           ! *** Do not allow more bedload deposition than bedload mass available
    enddo
    DEP = SUM(DEPBL(1:NSEDS)) + SUM(DEPTSSB(1:NSEDS))         ! *** Total deposition, sum of bedload and suspended load deposition for all size classes
    DEP_SED_FLX(L,1:NSCM) = DEP_SED_FLX(L,1:NSCM) + DEPBL(1:NSCM) ! *** Total deposition flux with bedload and suspended load (g/cm^2)
  else
    DEP = SUM(DEPTSSB(1:NSEDS))                                ! *** Total deposition is the sum of suspended load deposition for all size classes
  endif
  
  ! *** TOTAL DEPOSITION
  !DEP_SED_FLX(L) = DEP/DTSEDJ                                 ! *** Calculate the deposition rate (g/cm^2/s)
  
  ! *** ********************************************************************************************************************************************
  ! *** Get things set up for Erosion Calculations
  ! *** Find the next layer (SLLN) of sediment below the top layer based on the LAYERACTIVE from the previous time step
  SLLN = KB
  if( LAYERACTIVE(2,L) == 1 )then   ! *** Check to ensure that there are at least 2 sediment layers present here
    SLLN = 2                        ! *** Layer number below active layer (always 3 to start after SEDIC, always 2 when there at least 2 layers)
  else
    do K = 3,KB
      if( LAYERACTIVE(K,L) > 0 .and. LAYERACTIVE(K-1,L) == 0 )then
        SLLN = K                      ! *** Next layer is renumbered as necessary (happens if active layer is eroded completely)
        EXIT
      endif
    enddo
  endif
  
  ! *** Calculate Average particle size of surface layer so we can calculate
  ! *** active layer unit mass
  if( LAYERACTIVE(1,L) == 1 )then
    SURFACE = 1                     ! *** Surface variable set to top layer
  else
    SURFACE = SLLN                  ! *** Otherwise the top layer is SLLN
  endif
  D50AVG(L) = SUM(PERSED(1:NSEDS,SURFACE,L)*D50(1:NSEDS))             ! *** Calculate local d50 at sediment bed surface
 
  ! *** Calculate TAUCRIT Based on the Average Particle Size of Surface
  ! *** Then calculate the Active Layer unit mass (TACT) from it.
  ! *** Ta =  Tam * Davg * (Tau/Taucr)
  TAUCRIT = 1.E6
  if( LAYERACTIVE(SURFACE,L) == 1 )then
    ! *** Active/deposited layer
    if(D50AVG(L) < SCND(1) )then
      NSCD(1) = SCND(1)
      NSCD(2) = SCND(2)
      NSC0 = 1
      NSC1 = 2
      D50AVG(L) = SCND(1)                                             ! *** Prevent division (s_shear) by zero when there is no sediment in the layer 
    elseif( D50AVG(L) >= SCND(NSICM) )then
      NSCD(1) = SCND(NSICM-1)
      NSCD(2) = SCND(NSICM)
      NSC0 = NSICM-1
      NSC1 = NSICM
    else  
      do NS = 1,NSICM-1
        if( D50AVG(L) >= SCND(NS) .and. D50AVG(L) < SCND(NS+1) )then
          NSCD(1) = SCND(NS)
          NSCD(2) = SCND(NS+1)
          NSC0 = NS
          NSC1 = NS+1
          EXIT
        endif
      enddo
    endif
  
    TAUCRIT = TAUCRITE(NSC0)+(TAUCRITE(NSC1)-TAUCRITE(NSC0))/(NSCD(2)-NSCD(1))*(D50AVG(L)-NSCD(1))
    TAUCOR(SURFACE,L) = TAUCRIT
  elseif( LAYERACTIVE(SURFACE,L) == 2 )then
    ! *** In-place sediment layer
    TAUCRIT = TAUCOR(SURFACE,L)
  endif
  
  ! *** Compute the requried active layer thickness (cm)
  if( TAU(L)/TAUCRIT < 1.0 )then
    TACT = TACTM*D50AVG(L)*(BULKDENS(1,L)/10000.0)
  else
    TACT = TACTM*D50AVG(L)*(TAU(L)/TAUCRIT)*(BULKDENS(1,L)/10000.0)
  endif

  ! *** This is where we determine if there is an un-erodeable size class.  
  ! *** If there is one, we need an active layer.
  NACTLAY = 0
  do K = 1,KB                            ! *** Search through all layers
    if( LAYERACTIVE(K,L) > 0 )then     ! *** If the layer is present
      do NS = 1,NSEDS
        ! *** TAU is > the critical shear stress for erosion based on the d50 but less than the critical shear stress for one or more size classes
        if( PERSED(NS,K,L) > 0.0 .and. TAU(L) < TCRE(NS) .and. TAU(L) > TAUCRIT )then !if the size class is present and the sediment bed is eroding and there is insufficient shear stress for suspension
          NACTLAY = 1                  ! *** There is an active layer (logical variable)
          LAYERACTIVE(1,L) = 1         ! *** There is an active layer
          EXIT
        endif
      enddo
    endif
  enddo

  ! *** If the layer exposed can erode, then 
  ! *** Use the active layer model, otherwise just put deposited material
  ! *** on top.  Also if there is a size class present in the top layer
  ! *** that will not erode, create an active layer.
  SURFOLD = 0
  if( TSED(1,L) > 0.0 .or. NACTLAY /= 0 )then      ! *** If there is mass in the active layer, we must go through the sorting routine
    ! *** No active layer for pure erosion (active layer needed for coarsening and deposition)
  
    ! *** Sort layers so that the active layer is always Ta thick.
    ! *** Recalculate the mass fractions after borrowing from lower layers
    if( LAYERACTIVE(1,L) == 1 )then
      if( TSED(1,L) > TACT )then                        ! *** At this point TSED does not include the deposited sediment for this time step            
        ! *** There is deposition over this time step.  Redistribute excess mass to the deposition layer (i.e. layer 2)
        FORALL(NS = 1:NSEDS)PERSED(NS,2,L) = (PERSED(NS,2,L)*TSED(2,L)+PERSED(NS,1,L)*(TSED(1,L)-TACT))/(TSED(2,L)+(TSED(1,L)-TACT))  ! *** Recalculate mass fractions
        LAYERACTIVE(2,L) = 1                            ! *** Ensure that the second layer logical is turned on
        SLLN = 2                                        ! *** Next lower layer is 2
        TSED(2,L) = TSED(2,L) + TSED(1,L) - TACT        ! *** Add layer unit mass in excess of active-layer unit mass to next lower layer
        TSED(1,L) = TACT                                ! *** Reset top layer unit mass to active layer unit mass
        
      elseif( TSED(1,L) < TACT .and. TSED(1,L)+TSED(SLLN,L) > TACT .and. TAU(L) > TAUCOR(SLLN,L) )then
        ! *** There is net erosion over this time step and there is sufficient sediment below to reconstitute the active layer
        FORALL(NS = 1:NSEDS) PERSED(NS,1,L) = (PERSED(NS,1,L)*TSED(1,L) + PERSED(NS,SLLN,L)*(TACT-TSED(1,L)))/TACT     ! *** Recalculate the mass fraction
        TSED(SLLN,L) = TSED(SLLN,L) - (TACT-TSED(1,L))  ! *** Borrow unit mass from lower layer
        TSED(1,L) = TACT                                ! *** Reset top layer unit mass to active layer unit mass
        
      elseif( TSED(1,L) < TACT .and. TSED(1,L)+TSED(SLLN,L) <= TACT .and. TAU(L) > TAUCOR(SLLN,L) )then
        ! *** There is net erosion over this time step and there is NOT sufficient sediment below to reconstitute the active layer
        FORALL(NS = 1:NSEDS)
          PERSED(NS,1,L) = (PERSED(NS,1,L)*TSED(1,L) + PERSED(NS,SLLN,L)*(TSED(SLLN,L)))/(TSED(1,L) + TSED(SLLN,L))  ! *** Recalculate the mass fraction
          PERSED(NS,SLLN,L) = 0.0                       ! *** No more sediment available in next lower layer
        ENDFORALL
        TSED(1,L) = TSED(1,L) + TSED(SLLN,L)            ! *** Add available sediment to layer
        SURFOLD = KBT(L)                                ! *** Save the residual layer
        TSED(SLLN,L) = 0.0                              ! *** No more sediment available below this
        LAYERACTIVE(SLLN,L) = 0                         ! *** Layer has been eliminated in the logical variable
        if( SLLN < KB )then
          SLLN = SLLN + 1                               ! *** Set next layer lower
          do while (TSED(SLLN,L) <= 0. )
            SLLN = SLLN + 1                                 
            if( SLLN > KB )then
              SLLN = KB                                 ! *** Do not allow specification of the next lower layer to be below the bottom sediment layer
              EXIT
            endif
          enddo
        endif
      endif
    endif
  endif

  ! *** Propwash
  if( NACTIVESHIPS > 0 )then                                   ! *** Propwash is active and there are active ships. 
    PROP_ERO(L,1:NSEDS) = PROP_ERO(L,1:NSEDS)/DXYP(L)/10000.   ! *** Convert mass from g to g/cm^2.  Ignore fast classes since not fractionated yet
  endif
  
  ! *** ********************************************************************************************************************************************
  ! *** Now calculate the Erosion Rates
  K1 = 0
  do K = 1,KB                                          ! *** Loop through all sediment layers so that they are properly eroded and sorted
    if( LAYERACTIVE(K,L) == 0 ) CYCLE                ! *** If the layer is gone don't consider it
    if( LAYERACTIVE(1,L) == 1 .and. K /= 1 )then
      if( PROP_ERO(L,0) == 0.0 )then
        EXIT                                         ! *** If it is depositional, there is no need to consider erosion.  Don't exit if propwash is active
      endif
    endif
    if( K > 1 )then
      if( LAYERACTIVE(K-1,L) > 0 )then
        EXIT                                         ! *** Exit loop if layer above still has mass (i.e. not eroded completely, so erosion flux satisfied)
      endif
    endif

    D50AVG(L) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))  ! *** Find mean diameter of Layer
    
    ! *** Find upper and lower limits of size classes on mean bed diameter
    if( (D50AVG(L)+1E-6) < SCND(1) )then
      if( TSED(k,L) > 1E-8 )then
        PRINT '("Limits!  L: ",I6,", K: ", I3,", BED MASS: ",E14.6,", COMPUTED D50: ",f10.4,", MIN GRAINSIZE: ",F10.4)', Map2Global(L).LG, K, TSED(K,L), D50AVG(L), SCND(1)
        if( TSED(K,L) > 1E-4 )then
          call STOPP('COMPUTED BED D50 < MINIMUM BED GRAINSIZE CLASS!')
        else
          PERSED(1:NSEDS,K,L) = 0.0
          PERSED(1,K,L) = 1.0
          NS = 1
          NSCD(1) = SCND(NS)
          NSCD(2) = SCND(NS+1)
          NSC0 = NS
          NSC1 = NS+1
          D50AVG(L) = SCND(1)
        endif
      else
        PERSED(1:NSEDS,K,L) = 1./FLOAT(NSEDS)
        D50AVG(L) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))
      endif
    elseif( (D50AVG(L)-1E-6) > SCND(NSICM) )then
      temp = sum(PERSED(1:NSEDS,K,L))   ! delme
      print '(a,f15.5,2i5,f12.8,15e14.6)', 'D50 too big', timeday, l, k, temp, (PERSED(ns,K,L), ns = 1, nseds), prop_ero(l,0)   ! delme
      pause
      if( TSED(k,L) > 1E-8 )then
        PRINT '("Limits!  L: ",I6,", K: ", I3,", BED MASS: ",E14.6,", COMPUTED D50: ",f10.4,", MAX GRAINSIZE: ",F10.4)', Map2Global(L).LG, K, TSED(K,L), D50AVG(L), SCND(NSICM)
        if( TSED(K,L) > 1E-4 )then
          call STOPP('COMPUTED BED D50 > MAXIMUM BED GRAINSIZE CLASS!')
        else
          PERSED(1:NSEDS,K,L) = 0.0
          PERSED(NSEDS,K,L) = 1.0
          NS = NSICM - 1
          NSCD(1) = SCND(NS)
          NSCD(2) = SCND(NS+1)
          NSC0 = NS
          NSC1 = NS+1
          D50AVG(L) = SCND(NSICM)
        endif
      else
        PERSED(1:NSEDS,K,L) = 1./FLOAT(NSEDS)
        D50AVG(L) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))
      endif
    else   
      do NS = 1,NSICM-1
        if( D50AVG(L) >= SCND(NS) .and. D50AVG(L) < SCND(NS+1) )then
          NSCD(1) = SCND(NS)
          NSCD(2) = SCND(NS+1)
          NSC0 = NS
          NSC1 = NS+1
          EXIT
        endif
      enddo
    endif
    
    !if( timeday > 59.9 )then
    !  print '(f15.5,i10,i6,i3,i5,f10.4,2e12.4)', timeday, niter, Map2Global(L).LG, K, LAYERACTIVE(K,L), HP(L), TSED(k,L), TSED0(K,L) ! delme
    !endif

    ! *** Calculate TAUCRIT Based on the D50 of the bed or from Sedflume Data
    if( LAYERACTIVE(K,L) == 1 )then
      ! *** For active layers
      TAUCRIT = TAUCRITE(NSC0) + (TAUCRITE(NSC1)-TAUCRITE(NSC0))/(NSCD(2)-NSCD(1))*(D50AVG(L)-NSCD(1)) !interpolation
      TAUCOR(K,L) = TAUCRIT
    elseif( LAYERACTIVE(K,L) == 2 )then         ! *** IC Bed sediments
      ! *** SEDFlume data (depth interpolation)
      SN01 = TSED(K,L)/TSED0(K,L)               ! *** Weighting factor 1 for interpolation
      SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)  ! *** Weighting factor 2
      TAUCRIT = SN01*TAUCOR(K,L) + SN11*TAUCOR(K+1,L)
    endif
    ERO = 0.                                    ! *** Total erosion from the layer
    
    ! *** Check if the shear is greater than critical shears.  If not, exit erosion loop
    if( TAU(L) < SH_SCALE(L)*TAUCRIT .and. PROP_ERO(L,0) == 0.0 )then
      EXIT
    endif
    
    ! *** Now, calculate erosion rates  ----------------------------------------------------------------------------------
    ! *** Find the upper and lower limits of the Shear Stress for the interpolation
    if( TAU(L) >= TAULOC(ITBM) )then
      if( NWARNING < 100 )then
        PRINT '(A,I6,I5,E12.4)','*** WARNING  TAU >= MAXIMUM TAUCORE',Map2Global(L).LG,ICORE,TAU(L)
      endif
      NWARNING = NWARNING+1
      if( (MOD(NWARNING,100) == 0 .or. NWARNING == 1) .and. NDM == 1 )then  ! *** Can't write to log when running multi-threaded 
        open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        write(mpi_efdc_out_unit,'(A,I6,I5,F10.3,E12.4)')'*** WARNING  TAU >= MAXIMUM TAUCORE:  L,ICORE,TIMEDAY,TAU',Map2Global(L).LG,ICORE,TIMEDAY,TAU(L)
        close(mpi_efdc_out_unit)
      endif
      TAUDD(1) = TAULOC(ITBM-1)
      TAUDD(2) = TAULOC(ITBM)
      NTAU0 = ITBM-1
      NTAU1 = ITBM
      
    elseif( TAU(L) < TAULOC(1) )then
      if( NWARNING < 100 )then
        PRINT '(A,I6,I5,E12.4)','*** WARNING  TAU < MINIMUM TAUCORE',Map2Global(L).LG,ICORE,TAU(L)
      endif
      NWARNING = NWARNING+1
      if( (MOD(NWARNING,100) == 0 .or. NWARNING == 1) .and. NDM == 1  )then  ! *** Can't write to log when running multi-threaded 
        close(mpi_efdc_out_unit)
        open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        write(mpi_efdc_out_unit,'(A,I6,I5,F10.3,E12.4)')'*** WARNING  TAU < MINIMUM TAUCORE:  L,ICORE,TIMEDAY,TAU',Map2Global(L).LG,ICORE,TIMEDAY,TAU(L)
        close(mpi_efdc_out_unit)
      endif
      TAUDD(1) = TAULOC(1)
      TAUDD(2) = TAULOC(2)
      NTAU0 = 1
      NTAU1 = 2
    else
      do NS = 1,ITBM-1
        if( TAU(L) >= TAULOC(NS) .and. TAU(L) < TAULOC(NS+1) )then
          TAUDD(1) = TAULOC(NS)
          TAUDD(2) = TAULOC(NS+1)
          NTAU0 = NS
          NTAU1 = NS+1
          EXIT
        endif
      enddo
    endif
    
    ! *** Interpolate the erosion rates for shear stress and depth.
    ! *** This utilizes normal sedflume data for deeper layers.
    if( LAYERACTIVE(K,L) == 2 )then 
      ! *** Calculate erosion rates of deeper layers (SEDFlume data)
      if( NSEDFLUME == 1 )then
        SN00 = (TAUDD(2)-TAU(L))/(TAUDD(2)-TAUDD(1)) ! *** weighting factor 1 for interpolation
        SN10 = (TAUDD(1)-TAU(L))/(TAUDD(1)-TAUDD(2)) ! *** weighting factor 2
        SN01 = TSED(K,L)/TSED0(K,L)                  ! *** weighting factor 3
        SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)     ! *** weighting factor 4
        
        if( K+1 <= KB )then  ! *** Modeled erosion rate
          ERATEMOD = ( SN00*EXP(SN11*LOG(ERATE(K+1,L,NTAU0))+SN01*LOG(ERATE(K,L,NTAU0))) &
                   + SN10*EXP(SN11*LOG(ERATE(K+1,L,NTAU1))+SN01*LOG(ERATE(K,L,NTAU1))) )*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))
        else                 ! *** Do not allow erosion through the bottom layer
          ERATEMOD = ( SN00*EXP(SN11*LOG(1.0E-9)+SN01*LOG(ERATE(K,L,NTAU0))) &
                     + SN10*EXP(SN11*LOG(1.0E-9)+SN01*LOG(ERATE(K,L,NTAU1))) )*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))
        endif
      else
       
        if( TAU(L) > TAUCOR(K,L) )then                                               ! *** Check that the applied shear exceeds the critical shear stress for this layer
          ! *** Erosion rate values (cm/s) computed by equation assume shear in Pascals so convert dynes
          SN00 = EA(ICORE,K)*(0.1*TAU(L))**EN(ICORE,K)                               ! *** Erosion rate (cm/s) of the top layer

          if( K+1 <= KB )then
            SN10 = EA(ICORE,K+1)*(0.1*TAU(L))**EN(ICORE,K+1)                         ! *** Erosion rate (cm/s) of the layer below
          else
            SN10 = 0.0                                                               ! *** Modeled erosion rate in limited by bottom
          endif

          SN11 = (TSED0(K,L)-TSED(K,L))/TSED0(K,L)                                   ! *** Mass weighting factor
          ERATEMOD = ((SN10-SN00)*SN11 + SN00)*BULKDENS(K,L)*SQRT(ONE/SH_SCALE(L))   ! *** linear interpolation for remaining mass in current layer    (g/cm2/s)
          ERATEMOD = min(ERATEMOD,MAXRATE(ICORE,K))                                  ! *** Limit erosion rate
        else
          ERATEMOD = 0.0
        endif
      endif
    elseif( LAYERACTIVE(K,L) == 1 )then
      ! *** For Layers One and Two (the newly deposited sediments)
      ! *** The erosion rate for these layers is determined from 
      ! *** Sedflume experiments and is based on average particle
      ! *** Size (D50AVG) 
      NSCTOT = NSCD(2)-NSCD(1)                                                       ! *** difference in interpolant size class
      D50TMPP = D50AVG(L)-NSCD(1)                                                    ! *** difference from local size class and lower interpolant
      if( NSEDFLUME == 1 )then
        SN00 = (TAUDD(2)-TAU(L))/(TAUDD(2)-TAUDD(1))                                 ! *** weighting factor 1 for interpolation
        SN10 = (TAUDD(1)-TAU(L))/(TAUDD(1)-TAUDD(2))                                 ! *** weigthing factor 2
        SN01 = D50TMPP/NSCTOT                                                        ! *** weighting factor 3
        SN11 = (NSCTOT-D50TMPP)/NSCTOT                                               ! *** weighting factor 4
        ERATEMOD = (SN00*EXP(SN11*LOG(ERATEND(NSC0,NTAU0)) + SN01*LOG(ERATEND(NSC1,NTAU0))) + SN10*EXP(SN11*LOG(ERATEND(NSC0,NTAU1)) +  &   ! *** log-linear interpolation
                             SN01*LOG(ERATEND(NSC1,NTAU1))))*BULKDENS(K,L)*SQRT(1./SH_SCALE(L))
      else
        ! *** Erosion rate values (cm/s) computed by equation assume shear in Pascals so convert dynes
        SN00 = ACTDEPA(NSC0)*(0.1*TAU(L))**ACTDEPN(NSC0)                             ! *** Erosion rate 1 (cm/s)
        SN10 = ACTDEPA(NSC1)*(0.1*TAU(L))**ACTDEPN(NSC1)                             ! *** Erosion rate 2 (cm/s)
        SN11 = D50TMPP/NSCTOT                                                        ! *** Weighting factor 
        ERATEMOD = ((SN10-SN00)*SN11 + SN00)*BULKDENS(K,L)*SQRT(1./SH_SCALE(L))      ! *** linear interpolation around size class (g/cm2/s)
        ERATEMOD = min(ERATEMOD,ACTDEPMAX(NSC0))                                     ! *** Limit erosion rate
      endif
    endif

    ! *** Sort out Thicknesses and Erosion Rates
    EBD = ERATEMOD*DTSEDJ                                                            ! *** Maximum mass potentially eroded this time step for this layer (g/cm^2)

    ! *** If the shear stress is less than the critical shear stress for a
    ! *** particular size class, then it is not eroded from the bed.

    ! *** Calculate New sediment mass (TTEMP) of each sediment in Bed (g/cm^2)
    ! *** Conservation of sediment mass, you can only erode as much as is there.
    ! *** ETOT(NS) = Total erosion at this cell of size class NS
    ! *** ERO_SED_FLX(L,NSEDS2) = Erosion by sediment class for this cell
    ! *** ELAY(NS) = Erosion from this layer of size class NS
    ! *** ERO      = Total erosion from layer 
    if( PROP_ERO(L,0) > 0.0 )then                                                    ! *** Propwash is active and there are active ships
      ! *** Active ship traffic with erosion
      PSUS(L,1:NSEDS) = 1.0                                                          ! *** Set probability of erosion mass into suspension

      ! *** Include erosion due to ambient currents for all propwash options
      WHERE( TAU(L) >= TCRE(1:NSEDS) )
        ELAY(1:NSEDS) = PERSED(1:NSEDS,K,L)*EBD                                      ! *** Compute erosion due to ambient currents (g/cm^2)
      ENDWHERE
      ELAY(1:NSEDS)  = ELAY(1:NSEDS) + PROP_ERO(L,1:NSEDS)                           ! *** Add propwash induced erosion
      TTEMP(1:NSEDS) = PERSED(1:NSEDS,K,L)*TSED(K,L) - ELAY(1:NSEDS)                 ! *** Remaining mass in layer for each size class
      EBD = SUM(ELAY)                                                                ! *** Updated total mass erosion
    else
      ! *** Standard erosion processing
      WHERE( TAU(L) >= TCRE(1:NSEDS) )
        ELAY(1:NSEDS)  = PERSED(1:NSEDS,K,L)*EBD
        ETOT(1:NSEDS)  = ETOT(1:NSEDS) + ELAY(1:NSEDS)
        TTEMP(1:NSEDS) = PERSED(1:NSEDS,K,L)*TSED(K,L) - ELAY(1:NSEDS)               ! *** Remaining mass in layer for each size class
      elseWHERE
        ETOT(1:NSEDS)  = 0.0
        ELAY(1:NSEDS)  = 0.0
        TTEMP(1:NSEDS) = PERSED(1:NSEDS,K,L)*TSED(K,L)
      ENDWHERE
    endif
    ! *** At this point, ELAY only includes the total potential mass to be eroded frm the bed. Propwash Fast classes not addressed yet
    
    ! *** Ensure sufficient mass in current layer, otherwise empty the current layer
    ! *** and reduce erosion for next layer
    do NS = 1,NSEDS
      if( TTEMP(NS) < 0.0 )then
        ! *** The mass by class is negative, so zero the mass for that class
        TTEMP(NS) = 0.0                                                              ! *** Set remaining mass to zero, i.e. class fully eroded from this layer
        ELAY(NS)  = PERSED(NS,K,L)*TSED(K,L)                                         ! *** Only allow available mass to erode
        if( PROP_ERO(L,0) > 0. )then
          ETOT(NS) = ELAY(NS)                                                        ! *** Reset total erosion to mass available
          if( ETOT(NS) < 0.0 )then
            ! *** Should never trigger this event
            print '(a,f10.4,i8,3i5,f9.5,6e12.4)', 'Bad mass accounting N, L, NS = ', timeday, NITER, L, K, NS, persed(ns,k,l), tsum, PROP_ERO(L,ns), EBD, tsed(k,l), elay(ns), etot(ns)
            ETOT(NS) = 0.0
            ELAY(NS) = 0.0
          endif
        else
          ETOT(NS)  = ETOT(NS) - PERSED(NS,K,L)*EBD + ELAY(NS)                       ! *** Recalculate total erosion
        endif
      endif
      
      ! *** Mass erosion for propwash
      if( PROP_ERO(L,NS) > 0.0 )then
        if( IBED2WC(NS) > NSEDS )then
          ! *** Split eroded mass into fast and normal WC classes
          ETOT(NS) = (1.-fraction_fast)*ELAY(NS)                                     ! *** Reset total erosion to mass available for the base class
          ETOT(IBED2WC(NS)) = fraction_fast*ELAY(NS)                                 ! *** Reset total erosion to mass available for the "fast" class
          PROP_ERO(L,NS) = PROP_ERO(L,NS) - ELAY(NS)                                 ! *** Reduce propwash erosion by the amount removed from current layer
        else
          ETOT(NS) = ELAY(NS)                                                        ! *** Reset total erosion to mass available  
          PROP_ERO(L,NS) = PROP_ERO(L,NS) - ELAY(NS)                                 ! *** Reduce propwash erosion by the amount removed from current layer
        endif
        PROP_ERO(L,NS) = max(PROP_ERO(L,NS),0.0)
      endif
    enddo
    ERO = SUM(ELAY(1:NSEDS))                                                         ! *** Actual final total erosion from the layer   (g/cm^2)
    
    ! *** Subtract total erosion from layer unit mass then Calculate new percentages
    TEMP = TSED(K,L) - ERO                                                           ! *** Eroded layer unit mass.  TSED already has deposition added.
    TSUM = sum(TTEMP(:))
    
    if( TEMP < 1e-12 .or. TSUM <= 0.0 )then                                          ! *** If the remaining mass in the layer is negative, set its mass to zero
      TSED(K,L) = 0.0                                                                ! *** This layer has no mass
      LAYERACTIVE(K,L) = 0                                                           ! *** This layer is absent
      PERSED(1:NSEDS,K,L) = 0.0                                                      ! *** Zero mass fractions
      
    elseif( TEMP < 1e-4 .and. TSUM < 1e-4 )then
      ! *** If the sediment layer contains very little mass, just use the TTEMP mass to recompute PERSED
      TSED(K,L) = sum(TTEMP(:))
      PERSED(1:NSEDS,K,L) = TTEMP(1:NSEDS)/TSED(K,L)                                 ! *** New mass fractions
    
    else 
      ! *** Standard erosion processing
      RATIOMASS = TEMP/TSUM
      TTEMP(1:NSEDS) = TTEMP(1:NSEDS)*RATIOMASS                                      ! *** Ensure mass balance to the precision of the compiler
      TSED(K,L) = TEMP                                                               ! *** New layer unit mass (g/cm^2)

      PERSED(1:NSEDS,K,L) = TTEMP(1:NSEDS)/TSED(K,L)                                 ! *** New mass fractions
      
      ! delme
      if( abs(RATIOMASS - ONE) > 1e-8 )then
        SUMPER = sum(PERSED(1:NSEDS,K,L))
        print '(a,f15.5,i10,2i5,6e18.10,f10.1,f10.7)', 'Bad PERSED: ',  TIMEDAY, NITER, L, K, SUMPER, INITMASS(K), ERO, TSED(K,L),   &
                                                       sum(TTEMP(:)), TEMP!, sum(prop_ero(l,1:nseds)) , PROP_ERO(L,0)/DXYP(L)/10000.  ! delme
        do ns = 1,nseds
          print '(a,f15.5,i10,2i5,4e18.10,i5)', '            ',  TIMEDAY, NITER, L, K, PERSED(NS,K,L), ELAY(NS), TTEMP(NS), PMC(K,NS), NS  ! delme
        enddo
        
        ! *** If the sediment layer contains very little mass, just use the TTEMP mass to recompute PERSED
        if( TSED(K,L) < 1e-3 )then
          TSED(K,L) = sum(TTEMP(:))
          PERSED(1:NSEDS,K,L) = TTEMP(1:NSEDS)/TSED(K,L)                             ! *** New mass fractions
        else
          print '(a,e16.8,a,e16.8)', 'Stopping.  Computed Total Sediment = ', TSED(K,L), ',  Sum of Class by Class = ', sum(TTEMP(:))
          call STOPP('.')
        endif
        pause   ! DELME
      endif
    endif                                                                            

  enddo   ! ALL_LAYERS                                                               
  ERO_SED_FLX(L,1:NSEDS2) = ERO_SED_FLX(L,1:NSEDS2) + ETOT(1:NSEDS2)                 ! *** Accumulate the erosion flux in cell by each sediment class (g/cm^2)
  
  ! *** Add deposition to top layer
  if( DEP > 0.0 )then                                                                  
    ! *** There is deposition, calculate the new layer unit mass and mass fractions
    LAYERACTIVE(1,L) = 1                                                             ! *** Top sediment bed layer exists (because there is deposition)
    TTEMP(1:NSEDS) = PERSED(1:NSEDS,1,L)*TSED(1,L)                                   ! *** Mass for each size class for bed before deposition
    TSED(1,L) = TSED(1,L) + DEP                                                      ! *** Add the deposited mass to the active layer 1
    
    if( ICALC_BL > 0 )then
      PERSED(1:NSEDS,1,L) = ( TTEMP(1:NSEDS) + DEPTSSB(1:NSEDS) + DEPBL(1:NSEDS) )/TSED(1,L)  ! *** Bedload possible
    else
      PERSED(1:NSEDS,1,L) = ( TTEMP(1:NSEDS) + DEPTSSB(1:NSEDS) )/TSED(1,L)          ! *** No bedload
    endif
  endif
  
  ! *** DETERMINE TOTAL SEDIMENT FLUX
  if( ICALC_BL > 0 )then
    do NS = 1,NSEDS
      TEMP = (ONE-PSUS(L,NS))*ETOT(NS)                      ! *** Calculate erosion into bedload             (g/cm^2)
      EBL(L,NS)   = TEMP*10000.                             ! *** Save erosion into bedload                  (g/m^2)
      DBL(L,NS)   = DEPBL(NS)*10000.                        ! *** Save deposition from bedload               (g/m^2)
      QBLFLUX(L,NS) = TEMP - DEPBL(NS)                      ! *** Flux in/out of bed for bedload             (g/cm^2)
    enddo
    do NS = 1,NSEDS2
      QBFLUX(NS)    = PSUS(L,NS)*ETOT(NS) - DEPTSS(NS)      ! *** Flux in/out of bed for suspended suspended load  (g/cm^2)
    enddo
  else
    QBFLUX(1:NSEDS2) = ETOT(1:NSEDS2) - DEPTSS(1:NSEDS2)       ! *** Net erosion (+) / Deposition (-) from TSS  (g/cm^2)
  endif

  ! *** Update sediment variables for use in the remaining EFDC routines
  FORALL(K = 1:KB) SEDDIA50(L,K) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))
  HBED(L,1:KB) = 0.01_8*TSED(1:KB,L)/BULKDENS(1:KB,L)                                                    ! *** HBED-Bed height (m)  TSED-sediment layer unit mass (g/cm^2)  BULKDENS-Dry Density of Sediment (g/cm^3)
                                                                                                         
  ! *** Water column fluxes and concentrations                                                           
  WDTDZ                 = DTSEDJ*HPKI(L,KSZ(L))                                                          ! *** Delta t over Delta z
  SEDF(L,0,1:NSEDS2)     = QBFLUX(1:NSEDS2)*10000._8/DTSEDJ                                              ! *** SEDF-Suspended Sediment flux (g/m^2/s), QBFLUX (g/cm^2)
  SED(L,KSZ(L),1:NSEDS2) = SED2(L,KSZ(L),1:NSEDS2) + (SEDF(L,0,1:NSEDS2)-SEDF(L,KSZ(L),1:NSEDS2))*WDTDZ  ! *** SED-Suspended sediment concentration (g/m^3)
  
  ! *** Check for negative concentrations
  K = KSZ(L)
  do NS = 1,NSEDS2
    if( SED(L,K,NS) < 0. )then
      if( SED(L,K,NS) < -0.001 )then
        open(1,FILE = OUTDIR//'NEGSEDSND.OUT',POSITION = 'APPEND')  
        write(1,"(' Warning: WC  SED < 0: TIME, NS, L, I, J, K, HP, NEGSED = ',F12.4,5I5,4E13.4)" ) TIMEDAY, NS, Map2Global(L).LG, Map2Global(L).IG, Map2Global(L).JG, K, HP(L), SED(L,K,NS)  
        PRINT "(' Warning: WC  SED < 0: TIME, NS, L, I, J, K, HP, NEGSED = ',F12.4,5I5,4E13.4)", TIMEDAY, NS, Map2Global(L).LG, Map2Global(L).IG, Map2Global(L).JG, K, HP(L), SED(L,K,NS)  
        close(1)
      endif  
      SED(L,K,NS) = 0.0    ! *** Continue with warning
    endif  
  enddo
  
  ! *** Set EFDC standard bed mass variables
  SEDBT(L,1:KB)        = TSED(1:KB,L)*10000.                                                      ! *** SEDBT-Total sediment mass (g/m^2) in a layer, TSED-sediment layer unit mass (g/cm^2)
  do K = 1,KB
    SEDB1(L,K,1:NSEDS) = SEDB(L,K,1:NSEDS)
    SEDB(L,K,1:NSEDS)  = SEDBT(L,K)*PERSED(1:NSEDS,K,L)                                           ! *** SEDB-Sediment mass (g/m^2) by class in each layer
  enddo

  ! *** Update the volumetric flux rates (m/s) for sediment and porewater.
  QSBDTOP(L) = 0.                                                   
  QWBDTOP(L) = 0.     
  if( ICALC_BL > 0 )then
    do NS = 1,NSEDS
      SEDFLUX = SEDF(L,0,NS) + ( EBL(L,NS)-DBL(L,NS) )/DTSEDJ                                     ! *** Total sediment flux           (g/m^2/s)
      QSBDTOP(L) = QSBDTOP(L) + SSGI(NS)*SEDFLUX                                                  ! *** Volume of sediment exchange   (m/s)
      QWBDTOP(L) = QWBDTOP(L) + SSGI(NS)*SEDFLUX*VDRBED(L,KBT(L))                                 ! *** Volume of water exchange      (m/s)
    enddo
    if( NSEDS2 > NSEDS )then
      do NS = NSEDS+1, NSEDS2                                                                     ! *** Fast class deposition volumes
        QSBDTOP(L) = QSBDTOP(L) + SSGI(NS)*SEDF(L,0,NS)                                           ! *** Volume of sediment exchange   (m/s)
        QWBDTOP(L) = QWBDTOP(L) + SSGI(NS)*SEDF(L,0,NS)*VDRBED(L,KBT(L))                          ! *** Volume of water exchange      (m/s)
      enddo
    endif
  else
    do NS = 1,NSEDS2
      QSBDTOP(L) = QSBDTOP(L) + SSGI(NS)*SEDF(L,0,NS)                                             ! *** Volume of sediment exchange   (m/s)
      QWBDTOP(L) = QWBDTOP(L) + SSGI(NS)*SEDF(L,0,NS)*VDRBED(L,KBT(L))                            ! *** Volume of water exchange      (m/s)
    enddo
  endif
  
  ! *** Handle layers
  TSUM = HBED(L,1) + HBED(L,2)
  if( TSUM > 0. )then
    ! *** Active and/or Deposition layers exist.  Accumulate active and Deposition layers into one layer
    ! *** and collapse any empty layers between.
    SURFACE = -1
    KBT(L) = KB
    do K = 3,KB
      if( TSED(K,L) > 0. )then   
        SURFACE = K
        KBT(L) = K-1
        HBED(L,KBT(L))        = HBED(L,1)        + HBED(L,2) 
        SEDB(L,KBT(L),1:NSEDS) = SEDB(L,1,1:NSEDS) + SEDB(L,2,1:NSEDS)
        SEDBT(L,KBT(L))       = SEDBT(L,1)       + SEDBT(L,2) 
        EXIT
      endif
    enddo
    if( SURFACE == -1 )then
      ! *** All parent and deep layers are missing.  Only Layers 1 and 2 are active
      KBT(L) = KB
      HBED(L,KBT(L))        = HBED(L,1)        + HBED(L,2) 
      SEDB(L,KBT(L),1:NSEDS) = SEDB(L,1,1:NSEDS) + SEDB(L,2,1:NSEDS)
      SEDBT(L,KBT(L))       = SEDBT(L,1)       + SEDBT(L,2) 
    endif
    
    ! *** CFT linkage to propwash
    if( ISTRAN(5) > 0 )then
      ! *** TOXB(L,K)   - MG/M2
      ! *** TSED(K,L)   - G/CM2
      ! *** QBFLUX      - G/CM2
      ! *** HBED(L,1:KB) = 0.01*TSED(1:KB,L)/BULKDENS(1:KB,L)
      K1 = min(MAX(KBT(L)+1,2),KB)
      TEMP = TSED(K1,L) - INITMASS(K1)
      
      if( ICALC_BL > 0 )then
        do NT = 1,NTOX
          TOXB1(L,KBT(L),NT) = TOXB(L,KBT(L),NT)    ! *** Required for bedload transport
        enddo
      endif
      
      ! *** ADJUST TOXB CONCENTRATIONS FOR BED LAYER CHANGES
      if( SURFACE > -1 .and. TEMP /= 0. .and. INITMASS(K1) > 0 )then
        TEMP1 = TEMP/INITMASS(K1)
        ! *** ALLOCATE TOXIC MASS BY BED EXCHANGE
        if( TEMP1 < 0.0 )then
          ! *** SEDIMENT HAS BEEN MOVED TO THE ACTIVE/BUFFER LAYERS
          ! *** CALCULATE MASS OF TOXICS EXCHANGED BETWEEN KBT AND KBT+1
          do NT = 1,NTOX
            TEMP2 = TEMP1*TOXB(L,KBT(L)+1,NT)
            TOXB(L,KBT(L),NT)   = TOXB(L,KBT(L),NT)   - TEMP2
            TOXB(L,KBT(L)+1,NT) = TOXB(L,KBT(L)+1,NT) + TEMP2
          enddo
        elseif( TEMP > 0. )then
          ! *** SEDIMENT HAS BEEN MOVED FROM THE ACTIVE/BUFFER LAYERS
          do NT = 1,NTOX
            TEMP2 = TEMP1*TOXB(L,KBT(L),NT)
            TOXB(L,KBT(L),NT)   = TOXB(L,KBT(L),NT)   - TEMP2
            TOXB(L,KBT(L)+1,NT) = TOXB(L,KBT(L)+1,NT) + TEMP2
          enddo
        endif
      endif
      
      if( SURFOLD > 2 .and. SURFOLD /= KBT(L) )then
        ! *** OLD PARENT LAYER COMPLETELY ERODED.  ADD MASS OF TOXICS FROM SURFOLD TO KBT
        PRINT '(A,I10,F14.6,I6,2I4,3X,3E12.4)','Eroded through Layer:', N, TIMEDAY, Map2Global(L).LG, SURFOLD, KBT(L), TOXB(L,SURFOLD,1), TOXB(L,KBT(L),1), SEDF(L,0,1)
        do NT = 1,NTOX
          TOXB(L,KBT(L),NT)   = TOXB(L,KBT(L),NT) + TOXB(L,SURFOLD,NT)
          TOXB(L,SURFOLD,NT)  = 0.
        enddo
        SEDB1(L,KBT(L),1:NSEDS) = SEDB1(L,KBT(L),1:NSEDS) + SEDB1(L,SURFOLD,1:NSEDS)
      endif
    endif
    
    ! *** Zero any layers above KBT
    K1 = max(KBT(L)-1,1)
    do K = K1,1,-1
      HBED(L,K)         = 0.
      SEDB(L,K,1:NSEDS)  = 0.
      SEDB1(L,K,1:NSEDS) = 0.
      SEDBT(L,K)        = 0.
    enddo

    ! *** Optionally use maximum layer thickness
    if( HBEDMAX > 0.0 )then
      ! *** Check if Deposition layer > Max Layer Thickness
      TEMP1 = TACT/BULKDENS(2,L)*0.01                                       ! *** Thickness of active layer (m)
      if( KBT(L) > 2 .and. HBED(L,KBT(L)) > HBEDMAX + TEMP1 + 0.01 )then    ! *** Exclude the active layer thickness from max thickness
        ! *** The KBT layer is thicker than the max layer thickness and an empty layer exists below
        TEMP = HBEDMAX/HBED(L,KBT(L))         ! *** Ratio of thickness reduction
        K = KBT(L) - 1                        ! *** New KBT after splitting
        K1 = KBT(L)                           ! *** New deep layer after splitting
        
        HBED(L,K)  = HBED(L,K1) - HBEDMAX
        HBED(L,K1) = HBEDMAX
        
        SMASS(1:NSEDS) = SEDB(L,K1,1:NSEDS)
        SEDB(L,K,1:NSEDS)  = SMASS(1:NSEDS)*(ONE-TEMP)
        SEDB(L,K1,1:NSEDS) = SMASS(1:NSEDS) - SEDB(L,K,1:NSEDS)
        SEDBT(L,K)  = SUM(SEDB(L,K,1:NSEDS))
        SEDBT(L,K1) = SUM(SEDB(L,K1,1:NSEDS))
        
        TEMP1 = TSED(1,L) + TSED(2,L)
        TEMP2 = HBEDMAX*100.*BULKDENS(2,L)    ! *** Mass of sediment of HBEDMAX thickness   (g/cm^2)
        if( (TSED(2,L)-TEMP2) > 0.0 )then
          ! *** Deposition layer has sufficient sediment mass          
          PERSED(1:NSEDS,K1,L) = PERSED(1:NSEDS,2,L)
          TSED(2,L)  = TSED(2,L) - TEMP2
          TSED(K1,L) = TEMP2
          BULKDENS(K1,L) = BULKDENS(2,L)
        else
          ! *** Deposition layer has insufficient sediment mass   
          PERSED(1:NSEDS,1,L) = (PERSED(1:NSEDS,1,L)*TSED(1,L) + PERSED(1:NSEDS,2,L)*(TSED(2,L)))/(TSED(1,L) + TSED(2,L))
          TSED(1,L)  = TEMP1 - TEMP2
          
          PERSED(1:NSEDS,K1,L) = PERSED(1:NSEDS,1,L) 
          TSED(K1,L) = TEMP2

          PERSED(1:NSEDS,2,L) = 0.0
          TSED(2,L) = 0.0
        endif
        
        if( ISTRAN(5) > 0 )then
          do NT = 1,NTOX
            TEMP1 = TOXB(L,K1,NT)
            TOXB(L,K,NT)  = TEMP1*(ONE - TEMP)
            TOXB(L,K1,NT) = TEMP1 - TOXB(L,K,NT)
          enddo        
        endif
        
        KBT(L) = K 

        ! *** SET ACTIVE LAYER FLAG
        do K = 1,KB
          if( TSED(K,L) > 0.0 )then
            if( LAYERACTIVE(K,L) < 2 ) LAYERACTIVE(K,L) = 1            ! *** Preserve in-place sediment layer flag LAYERACTIVE(K,L) = 2
          else
            LAYERACTIVE(K,L) = 0
            TSED(K,L) = 0.0
          endif
        enddo
        
      endif
    endif
    
  endif   ! *** End of Active Layer Exchange 

  ! *** Setup chemical processes in the sediment bed
  if( ISTRAN(5) > 0 )then
    if( ISGWIT > 0. )then  
      !K = KBT(L)
      !VOIDCON1 = VDRBED(L,K)                                          ! *** VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
      !HBEDTMP = (1. + VOIDCON1)*HBED(L,K)/(1. + VDRBED(L,K))  
      !HBEDTMP = HBED(L,K)                                              
      !TMPVALO = VDRBED(L,K)*HBED(L,K)/(1. + VDRBED(L,K))  
      !TMPVALN = VOIDCON1*HBEDTMP/(1. + VOIDCON1)  
      !QWBDTOP(L) = DELTI*(TMPVALO-TMPVALN)  
      !HBED(L,K) = HBEDTMP  
      !QWTRBED(L,K) = QWBDTOP(L) + QGW(L)/DXYP(L) 
      !VDRBED(L,K) = VOIDCON1
    
      ! *** Assumes void ratio does not change. (See commented out code above)
      do K = 0,KBT(L)
        QWTRBED(L,K) = QGW(L)/DXYP(L)                                  ! *** QWTRBED is the seepage velocity in m/s
      enddo  
    endif
  endif
  
  return
  
END SUBROUTINE SEDZLJ

