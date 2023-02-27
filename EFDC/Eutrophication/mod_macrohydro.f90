! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE WQ_MACROPHYTE_FEEDBACK
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Module: Module name WQ_MACROPHYTE_FEEDBACK
!
!> @details Handle the option feedback of macrophyte growth and density on the 
!           hydrodynamics.  Part of the eutrophication module.
!
!> @author Paul M. Craig
!> @date 06/2021
!---------------------------------------------------------------------------!
  
  Use GLOBAL    
  Use Variables_WQ
  Use Variables_MPI
  
  IMPLICIT NONE

  REAL    :: BETAMAC_P = 1.0, BETAMAC_D = 0.39, CE4MAC = 1.05 ! *** James & O'Donncha (2019)
  
CONTAINS

!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine Feedback
!
!> @details Set parameters for vegetative drag to be used in CALEXP and CALTBXY
!---------------------------------------------------------------------------!
!> @author P.M. Craig
!---------------------------------------------------------------------------!
SUBROUTINE feedback
  
  integer :: NAL

  DO NAL = 1,NALGAE
    !IF( .NOT. ALGAES(NAL).ISMOBILE .AND. K == KSZ(L) )THEN
    !  ! *** Macrophytes and periphyton
    !  WQA6A(NAL) = ALGAES(NAL).WQFCDB + (1.0 - ALGAES(NAL).WQFCDB)*(ALGAES(NAL).WQKHRA(IZ)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)) 
    !  WQA6 = WQA6 + ( WQA6A(NAL)*WQBM(L,NAL) + ALGAES(NAL).WQFCDP*WQPR(L,NAL) )*WQVO(L,K,19+NAL) 
    !ENDIF
  ENDDO
    
END SUBROUTINE feedback
  
SUBROUTINE Macro_Veg(NAL)
  
  INTEGER, INTENT(IN) :: NAL
    
  INTEGER :: L, K, K1, K2, KTOP
  REAL :: LAYER_THRESHOLD, THICK, EXCESS_MASS, TVAL1, TVALHEI
    
  DO L = 2,LA
    ! *** Set active layers               delme - TODO - handle growth downward from suspended base
    DO K = KSZ(L),KC
      IF( HP(L)*Z(L,K-1) >= (HP(L)-ALGAES(NAL).BASEDEPTH+HEIGHT_MAC(L,NAL)) )THEN
        LAYERTOP(NAL,L) = K                                                    ! *** Top active layer
        EXIT                                                                   ! *** Jump out of the layer loop
      ENDIF
    ENDDO
    
    ! *** Set vegetation height/diameter
    MAC_CELL(L) = .FALSE.                                                      ! *** Reset total vegetation height flag
    HEIGHT_MAC(L,NAL) = 0.0
    LayerRatio_MAC(L,:,NAL) = 0.0
    DIAMETER_MAC(L,NAL) = ALGAES(NAL).MINDIAMETER
    
    if( timeday > 180.5 )then
      k = 1  ! delme
    endif
    
    KTOP = 0
    K2 = MIN(LAYERTOP(NAL,L)+1,KC)
    ! *** Loop over layers to redistribute vegetation mass
    DO K = LAYERBOT(NAL,L), K2
      IF( WQV(L,K,19+NAL) < 1E-8 ) CYCLE
      
      if( k == kc )then
        k1 = 0  ! delme
      endif
      
      ! *** Allow for growth between layers, if allowed
      IF( WQV(L,K,19+NAL) > ALGAES(NAL).THRESHOLD .AND. K < KC )THEN
        ! *** Redistribute mass to other layers, i.e. "Grow" any excess mass into the layer above
        !LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K), 0.1)               ! *** Handle thin layers (g C/m3)
        TVAL1 = WQV(L,K,19+NAL) - ALGAES(NAL).THRESHOLD
        TVAL1 = TVAL1*HPK(L,K)                                                 ! *** Convert to mass of macrophyte carbon to move/grow from layer
        TVAL1 = TVAL1*HPKI(L,K+1)                                              ! *** Convert back to concentration based on layer thickness above.
        WQV(L,K,19+NAL)   = ALGAES(NAL).THRESHOLD
        WQV(L,K+1,19+NAL) = WQV(L,K+1,19+NAL) + TVAL1

        WQVO(L,K,19+NAL)   = ALGAES(NAL).THRESHOLD
        WQVO(L,K+1,19+NAL) = WQVO(L,K+1,19+NAL) + TVAL1
        !ELSE
          !! *** Macrophyte is at its full depth but mass is larger than threshold.
          !! *** Therefore, increase stem diameter to account for mass
          !THICK = Z(L,K) - Z(L,LAYERBOT(NAL,L)-1)
          !IF( THICK <= 0. ) THICK = 1.0                                             ! *** Handle single layer models
          !EXCESS_MASS = WQV(L,K,19+NAL) - ALGAES(NAL).THRESHOLD
          !
          !! *** Distribute excess mass
          !DO K1 = LAYERBOT(NAL,L),K
          !  LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K1), 0.1)             ! *** Handle thin layers
          !  IF( K1 < LAYERTOP(NAL,L) )THEN
          !    WQV(L,K1,19+NAL)  = WQV(L,K1,19+NAL) + EXCESS_MASS*DZC(L,K1)/THICK    ! *** Update layer concentrations
          !    WQVO(L,K1,19+NAL) = WQV(L,K1,19+NAL)                                  ! *** Update layer concentrations
          !  ENDIF
          !  DIAMETER_MAC(L,NAL) = DIAMETER_MAC(L,NAL)*EXCESS_MASS/LAYER_THRESHOLD/ALGAES(NAL).STEMDENSITY   ! *** Update stem diameter
          !  !IF( DIAMETER_MAC(L,NAL) < ALGAES(NAL).MINDIAMETER )THEN
          !  !  PRINT '(I6,2I5,F8.3,4F10.3)', NITER, NAL, L, HP(L), ALGAES(NAL).THRESHOLD, LAYER_THRESHOLD, DIAMETER_MAC(L,NAL), ALGAES(NAL).MINDIAMETER  ! DELME
          !  !ENDIF
          !  DIAMETER_MAC(L,NAL) = MAX(DIAMETER_MAC(L,NAL),ALGAES(NAL).MINDIAMETER)
          !ENDDO
          !EXIT
        !ENDIF
      ENDIF
      KTOP = K
    ENDDO
    K1 = KTOP
    
    ! *** Update vegetation heights
    DO K = K1, LAYERBOT(NAL,L), -1
      IF( WQV(L,K,19+NAL) > 1E-8 )THEN
        !LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K), 0.1)             ! *** Handle thin layers (g C/m3)
        IF( K >= LAYERBOT(NAL,L)+1 )THEN
          IF( WQV(L,K-1,19+NAL) < ALGAES(NAL).THRESHOLD )THEN
            ! *** Move mass to layer below
            TVAL1 = ALGAES(NAL).THRESHOLD - WQV(L,K-1,19+NAL)
            TVAL1 = MAX(TVAL1,0.0)
            TVAL1 = TVAL1*HPK(L,K-1)                                           ! *** Convert to mass of macrophyte carbon to move/grow from layer
            TVAL1 = TVAL1*HPKI(L,K)                                            ! *** Convert back to concentration based on layer thickness above.
            WQV(L,K-1,19+NAL) = ALGAES(NAL).THRESHOLD
            WQV(L,K,19+NAL)   = WQV(L,K,19+NAL) - TVAL1
            IF( WQV(L,K,19+NAL) < 1E-8 )THEN
              WQV(L,K,19+NAL) = 0.0
              KTOP = MAX(KTOP-1, LAYERBOT(NAL,L))
            ENDIF
            
            WQVO(L,K-1,19+NAL) = ALGAES(NAL).THRESHOLD
            WQVO(L,K,19+NAL)   = WQVO(L,K,19+NAL) - TVAL1
            WQVO(L,K,19+NAL)   = MAX(WQVO(L,K,19+NAL), 0.0)
          ENDIF
        ENDIF
        
        ! *** Complete or partial penetration of current layer
        TVALHEI = MIN(WQV(L,K,19+NAL)/ALGAES(NAL).THRESHOLD, 1.0)              ! *** Height of macrophyte in current layer   DELME - /ALGAES(NAL).STEMDENSITY HERE AND EVERYWHERE ?
        HEIGHT_MAC(L,NAL) = HEIGHT_MAC(L,NAL) + TVALHEI*HPK(L,K)               ! *** Total height of macrophyte
        LayerRatio_MAC(L,K,NAL) = TVALHEI                                      ! *** Fraction of layer containing vegetation  
        IF( LayerRatio_MAC(L,K,NAL) < 0.01 ) LayerRatio_MAC(L,K,NAL) = 0.0     ! *** Ignore very small heights for hydrodynamics
      ENDIF
    ENDDO
  
    ! *** 
    IF( KTOP == 0 )THEN
      KTOP = LAYERBOT(NAL,L)
      LAYERTOP(NAL,L) = KTOP
      K = KTOP
      !LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD/MAX(HPK(L,K), 0.1)               ! *** Handle thin layers (g C/m3)
      
      ! *** Partial penetration of current layer
      TVALHEI = MIN(WQV(L,K,19+NAL)/ALGAES(NAL).THRESHOLD, 1.0)                ! *** Height of macrophyte in current layer   DELME - /ALGAES(NAL).STEMDENSITY HERE AND EVERYWHERE ?
      HEIGHT_MAC(L,NAL) = HEIGHT_MAC(L,NAL) + TVALHEI*HPK(L,K)                 ! *** Total height of macrophyte
      LayerRatio_MAC(L,K,NAL) = TVALHEI                                        ! *** Fraction of layer containing vegetation  
      IF( LayerRatio_MAC(L,K,NAL) < 0.01 ) LayerRatio_MAC(L,K,NAL) = 0.0       ! *** Ignore very small heights for hydrodynamics
    ELSE
      LAYERTOP(NAL,L) = KTOP      
    ENDIF
          
    IF( HEIGHT_MAC(L,NAL) > 0.0 ) MAC_CELL(L) = .TRUE.                         ! *** Total vegetation height is used as a flag to compute turbulence
    
          ! *** ACCUMULATE ACTIVE VEGETATION LAYERS
          !VEGK(L) = VEGK(L) + HVGTC*HPKI(L,K)

  ENDDO        ! *** End of L index loop  
    
  TVAL1 = SUM( WQV(2,:,20)*HPK(2,:) ) ! DELME
  
  L = 2  ! DELME
  IF( KTOP == KC )THEN
    K = 0  ! DELME
  ENDIF
  IF( WQV(2,KC,19+NAL) > ALGAES(NAL).THRESHOLD .AND. Ktop == KC )THEN
    K = 0  ! DELME
  ENDIF
  
  
END SUBROUTINE Macro_Veg
  
END MODULE
  
