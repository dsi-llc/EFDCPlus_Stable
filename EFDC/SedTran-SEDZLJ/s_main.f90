! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!**********************************************************************!
SUBROUTINE SEDZLJ_MAIN

  ! *** *******************************************************************!
  !
  ! ***  SUBROUTINE SEDZLJ_MAIN CALCULATES COHESIVE SEDIMENT SETTLING,
  ! ***  DEPOSITION AND RESUSPENSION ACCORDING TO SEDZLJ MODEL
  !
  ! ORIGINAL DATE :  May 24, 2006
  !  Craig Jones and Scott James

  ! 2016-11-07  Paul M. Craig  Updated for SIGMA-ZED AND OMP
  ! 2017-01-04  Paul M. Craig  Toxics and bedload mass balance updates
  ! 2021-06-07  Paul M. Craig  Updated for propwash
  ! 2022-01-07  Paul M. Craig  Added mass erosion classes for propwash

  ! *** *******************************************************************!

  use GLOBAL
  use Variables_Propwash
  implicit none

  real(RKD) :: CRNUM, SEDAVG, GRADSED, CLEFT, CRIGHT, WVEL, SMASSD, SMASSU, SMASS
  real(RKD),save :: LASTTIME = -9999.
  real(RKD),dimension(NSEDS) :: BLFLUXU, BLFLUXD
  integer :: L, K, NS, ND, LL, LF, LP, LUTMP, LDTMP, NSB
 
  ! *** ENFORCE STARTUP FOR ENTIRE SEDIMENT PROCESS
  if( LASTTIME == -9999. ) LASTTIME = TIMEDAY
  
  SSGI(1:NSEDS2) = 1.0/(1.0E6*SSG(1:NSEDS2)) 
  
  ! *** INITIALIZE CURRENT TIME STEP VARIABLES FOR SEDZLJ
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = MIN(LF+LDMSED-1,LASED)
    do LP = LF,LL
      L = LSED(LP)
      HPCM(L) = 100.0*HP(L)
    enddo
  enddo 
  !$OMP END PARALLEL DO

  ! *** Calculates the shear stresses from the velocity output of the hydraulic model.
  call SEDZLJ_SHEAR
  
  if( ISSLOPE ) CALL SEDZLJ_SLOPE
  ! *** Calling bed load calculation.  BEDLOADJ subroutine provides CBL, the sediment concentration load.
  if( ICALC_BL > 0 )then
    call BEDLOADJ
  endif
  
  ! *** Calls the sedflume transport. NOTE: ISEDTIME is the time difference between when sediment
  ! *** transport begins versus the beginning of hydraulic transport.
  if( NSEDFLUME > 0 )then

    ! *** **********************************************************************************!
    !   SEDZLJ Sediment Transport
    !
    !$OMP PARALLEL DEFAULT(SHARED)

    !$OMP DO PRIVATE(ND, NS, K, LP, L, WVEL, CLEFT, CRIGHT)
    do ND = 1,NDM  
      ! *** ZERO LAYER FLUXES, INCLUDING BOTTOM LAYER
      do NS = 1,NSEDS2
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            SEDF(L,K,NS) = 0.0
            WSETA(L,K,NS) = WSEDO(NS)   ! *** PMC - NOT NEEDED IF ISTOPT(6) = 0 
          enddo
          
          ! *** ZERO BED/WATER INTERFACE FLUX
          if( K == KS )then
            do LP = 1,LLWET(KC,ND)
              L = LKWET(LP,KC,ND)  
              SEDF(L,0,NS) = 0.0
            enddo
          endif
        enddo
      enddo
    
      !-----------------------------------------------------------------------------------!
      !
      ! *** HORIZONTAL LOOPS
      ! *** These loops account for sediment flux between the water layers due to settling.
      if( KC /= 1 )then
        do NS = 1,NSEDS2
          do LP = 1,LLWET(KC,ND)
            L = LKWET(LP,KC,ND) 
            WVEL = DTSEDJ*HPKI(L,KC)
            CLEFT = 1.0 + WSETA(L,KC-1,NS)*WVEL
            CRIGHT = MAX(SED(L,KC,NS),0.0)
            SED(L,KC,NS) = CRIGHT/CLEFT
            SEDF(L,KC-1,NS) = -WSETA(L,KC-1,NS)*SED(L,KC,NS)
          enddo
           
          if( KC /= 2 )then
            do K = KS,2,-1
              do LP = 1,LLWET(K-1,ND)
                L = LKWET(LP,K-1,ND) 
                WVEL = DTSEDJ*HPKI(L,K)
                CLEFT = 1.0+WSETA(L,K-1,NS)*WVEL
                CRIGHT = MAX(SED(L,K,NS),0.0)-SEDF(L,K,NS)*WVEL
                SED(L,K,NS) = CRIGHT/CLEFT
                SEDF(L,K-1,NS) = -WSETA(L,K-1,NS)*SED(L,K,NS)
              enddo
            enddo
          endif
        enddo
        
        ! *** Handle hard bottom for bottom water layer
        if( ISBEDMAP > 0 )then
          do LP = 1,LLWET(KC,ND)
            L = LKWET(LP,KC,ND)  
            if( LBED(L) .or. HP(L) < HDRY )then         ! *** Disable bed processes for very shallow depths (delme)
              ! *** ADD SETTLING FROM THE LAYER ABOVE
              K = KSZ(L)
              WVEL = DTSEDJ*HPKI(L,K)  
              do NS = 1,NSEDS2
                SED(L,K,NS) = MAX(SED(L,K,NS),0.) - SEDF(L,K,NS)*WVEL 
              enddo
            endif
          enddo
        endif
    
      endif
    enddo
    !$OMP END DO
  
    ! **********************************************************************************
    ! *** Erosion/Deposition Loop
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)
      do LP = LF,LL
        L = LSED(LP)
        ! *** Calculate total flux from bed and determine if there is enough sediment 
        ! *** if shear stesses are sufficient to cause erosion
        call SEDZLJ(L)
        
      enddo
    enddo
    !$OMP END DO

    ! **********************************************************************************
    ! *** Compute cell centered bedload fluxes for outflow or recirculation boundary
    !$OMP SINGLE
    if( NSBDLDBC > 0 )then
      BLFLUXU = 0.0
      BLFLUXD = 0.0
      do NSB = 1,NSBDLDBC
        LUTMP = LSBLBCU(NSB)
        LDTMP = LSBLBCD(NSB)

        SMASSU = TSED(1,LUTMP)
        do NS = 1,NSEDS
          if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
          BLFLUXU(NS) = ( DBL(LUTMP,NS) - EBL(LUTMP,NS) )/10000.             ! *** Flux in/out of bedload CBL                   (g/cm^2)
          TSED(1,LUTMP) = TSED(1,LUTMP) - BLFLUXU(NS)                        ! *** Back out erosion/depostion due to bedload    (g/cm^2)
          QBLFLUX(LUTMP,NS)  = -DBL(LUTMP,NS)/10000.                         ! *** Flux out of bed for bedload                  (g/cm^2)
          QSBDLDOT(LUTMP,NS) =  DBL(LUTMP,NS)*DXYP(LUTMP)/DTSEDJ             ! *** Sediment flux from upstream cell             (g/s)
          EBL(LUTMP,NS) = 0.0                            
        enddo
        
        ! *** UPDATE PERSED WITH MACHINE PRECISION A CONSIDERATION: UPSTREAM
        SMASS = 0.0
        do NS = 1,NSEDS
          SMASS = SMASS + MAX(PERSED(NS,1,LUTMP)*SMASSU - BLFLUXU(NS),0.0)
        enddo
        if( SMASS > 0.0 )then 
          do NS = 1,NSEDS
            PERSED(NS,1,LUTMP) = MAX(PERSED(NS,1,LUTMP)*SMASSU - BLFLUXU(NS),0.0)/SMASS
          enddo
          TSED(1,LUTMP) = SMASS
        else
          TSED(1,LUTMP) = 0.0
          PERSED(1:NSEDS,1,LUTMP) = 0.0
        endif

        ! *** Update EFDC variables
        SEDDIA50(LUTMP,1) = SUM(PERSED(1:NSEDS,1,LUTMP)*D50(1:NSEDS))
        HBED(LUTMP,1)  = 0.01_8*TSED(1,LUTMP)/BULKDENS(1,LUTMP)  
        SEDBT(LUTMP,1) = TSED(1,LUTMP)*10000.
        SEDB(LUTMP,1,1:NSEDS)  = SEDBT(LUTMP,1)*PERSED(1:NSEDS,1,LUTMP)

        if( LDTMP > 0 )then
          ! *** The flow is returning.  Add to downstream cell
          SMASSD = TSED(1,LDTMP)
          do NS = 1,NSEDS
            if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
            BLFLUXD(NS)    = DBL(LUTMP,NS)/10000.*DXYP(LUTMP)*DXYIP(LDTMP)     ! *** Depostion due to bedload from upstream cell  (g/cm^2)
            TSED(1,LDTMP) = TSED(1,LDTMP) + BLFLUXD(NS)                        ! *** Add mass to downstream cell                  (g/cm^2)
          enddo
          
          ! *** UPDATE PERSED WITH MACHINE PRECISION A CONSIDERATION: DOWNSTREAM
          SMASS = 0.0
          do NS = 1,NSEDS
            SMASS = SMASS + MAX(PERSED(NS,1,LDTMP)*SMASSD + BLFLUXD(NS),0.0)
          enddo
          if( SMASS > 0.0 )then 
            do NS = 1,NSEDS
              PERSED(NS,1,LDTMP) = MAX(PERSED(NS,1,LDTMP)*SMASSD + BLFLUXD(NS),0.0)/SMASS
            enddo
            TSED(1,LDTMP) = SMASS
          else
            TSED(1,LUTMP) = 0.0
            PERSED(1:NSEDS,1,LUTMP) = 0.0
          endif
          
          ! *** Update EFDC variables
          SEDDIA50(LDTMP,1) = SUM(PERSED(1:NSEDS,1,LDTMP)*D50(1:NSEDS))
          HBED(LDTMP,1)  = 0.01_8*TSED(1,LDTMP)/BULKDENS(1,LDTMP)  
          SEDBT(LDTMP,1) = TSED(1,LDTMP)*10000.
          SEDB(LDTMP,1,1:NSEDS)  = SEDBT(LDTMP,1)*PERSED(1:NSEDS,1,LDTMP)
        endif
          
      enddo
    endif
    !$OMP END SINGLE
    
    ! **********************************************************************************
    ! ***  ANTI-DIFFUSION AND FINAL BED/WATER FLUX FOR SEDIMENTS WHEN KC > 1
    if( KC > 1 .and. ISTOPT(6) > 0 )then
      !$OMP DO PRIVATE(ND,LF,LL,NS,K,LP,L,CRNUM,GRADSED,SEDAVG)
      do ND = 1,NDM  
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)
      
        do NS = 1,NSEDS2
          ! ***  ANTI-DIFFUSION OF SEDIMENT FOR KC > 1
          do K = 1,KS
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              CRNUM = 1.+DTSEDJ*WSETA(L,K,NS)*HPKI(L,K+1)
              GRADSED = (SED(L,K+1,NS)-SED(L,K,NS))/(DZC(L,K+1)+DZC(L,K))
              SEDAVG = 0.5*(SED(L,K+1,NS)+SED(L,K,NS)+1.E-16)
              WSETA(L,K,NS) = -CRNUM*DZC(L,K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
            enddo
          enddo

          ! *** TVAR1S = LOWER DIAGONAL
          do LP = LF,LL
            L = LWET(LP)  
            TVAR1S(L,KSZ(L)) = 0.
          enddo
          do K = 2,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR1S(L,K) = MIN(WSETA(L,K-1,NS),0.)
            enddo
          enddo

          ! *** TVAR1N = UPPER DIAGONAL
          do LP = LF,LL
            L = LWET(LP)  
            TVAR1N(L,KC) = 0.
          enddo
          do K = 1,KS
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR1N(L,K) = -MAX(WSETA(L,K,NS),0.)
            enddo
          enddo
        
          ! *** TVAR1W = MAIN DIAGONAL
          do LP = LF,LL
            L = LWET(LP)  
            TVAR1W(L,KSZ(L)) = DELTI*HPK(L,KSZ(L)) - MIN(WSETA(L,1,NS),0.)
            TVAR1W(L,KC)     = DELTI*HPK(L,KC)     + MAX(WSETA(L,KC-1,NS),0.)
          enddo
          do K = 2,KS
            do LP = 1,LLWET(K-1,ND)
              L = LKWET(LP,K-1,ND) 
              TVAR1W(L,K) = DELTI*HPK(L,K)+MAX(WSETA(L,K-1,NS),0.)-MIN(WSETA(L,K,NS),0.)
            enddo
          enddo

          ! *** TVAR1E = RIGHT HAND SIDE
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR1E(L,K) = DELTI*HPK(L,K)*SED(L,K,NS)
            enddo
          enddo

          ! *** TVAR3S = BET,TVAR2N = U,TVAR2S = GAM ARE WORKING ARRAYS
          do LP = LF,LL
            L = LWET(LP)  
            TVAR3S(L) = TVAR1W(L,KSZ(L))
          enddo
          do LP = LF,LL
            L = LWET(LP)  
            TVAR2N(L,KSZ(L)) = TVAR1E(L,KSZ(L))/TVAR3S(L)
          enddo
          do K = 2,KC
            do LP = 1,LLWET(K-1,ND)
              L = LKWET(LP,K-1,ND) 
              TVAR2S(L,K) = TVAR1N(L,K-1)/TVAR3S(L)
              TVAR3S(L) = TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
              TVAR2N(L,K) = (TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/TVAR3S(L)
            enddo
          enddo
          do K = KS,1,-1
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR2N(L,K) = TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
            enddo
          enddo
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              SED(L,K,NS) = TVAR2N(L,K)
            enddo
          enddo
          ! *** END OF ANTI-DIFFUSION CALCULATIONS
    
          ! *** RECOMPUTE FLUXES FROM FINAL SED: KC-1 LAYER
          do LP = LF,LL
            L = LWET(LP)  
            SEDF(L,KS,NS) = DELTI*HPK(L,KC)*( SED(L,KC,NS)-SED2(L,KC,NS) )
          enddo 

          ! *** RECOMPUTE FLUXES FROM FINAL SED: MIDDLE LAYERS
          do K = KS-1,1,-1
            do LP = 1,LLWET(K+1,ND)
              L = LKWET(LP,K+1,ND)  
              SEDF(L,K,NS) = DELTI*HPK(L,K+1)*( SED(L,K+1,NS)-SED2(L,K+1,NS) ) + SEDF(L,K+1,NS)
            enddo  
          enddo  
        enddo   ! *** End do of NS = 1,NSEDS loop
        
      enddo  ! *** End of Domain Loop
      !$OMP END DO
    endif    ! *** END OF KC > 1 ANTI-DIFFUSION BLOCK    
    
    !$OMP END PARALLEL
     
  endif
  ! *** End of condition  NSEDFLUME > 0 .and. N >= ISEDTIME
  ! ************************************************************************************

  ! *** PMC deleted the MORPHJ routine as it was redundant and not as complete as the morphology update in SSEDTOX  (2016-12)

  NDYCOUNT = NDYCOUNT+1

  return

END SUBROUTINE SEDZLJ_MAIN
