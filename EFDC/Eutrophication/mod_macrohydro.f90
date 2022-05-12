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
    
  INTEGER :: L, K, K1
  REAL :: LAYER_THRESHOLD, THICK, EXCESS_MASS, TVAL1
    
  DO L = 2,LA
    ! *** Set active layers
    DO K = KSZ(L),KC
      IF( HP(L)*Z(L,K-1) >= ALGAES(NAL).HEIGHTMIN )THEN
        LAYERBOT(NAL,L) = K                ! *** Bottom active layer
      ENDIF
      IF( HP(L)*Z(L,K) >= ALGAES(NAL).HEIGHTMAX .OR. K == KC )THEN
        LAYERTOP(NAL,L) = K                ! *** Top active layer
        EXIT                               ! *** Jump out of the layer loop
      ENDIF
    ENDDO
    
    ! *** Set vegetation height/diameter
    MAC_CELL(L) = .FALSE.                                                   ! *** Reset total vegetation height flag
    HEIGHT_MAC(L,NAL) = 0.0
    LayerRatio_MAC(L,:,NAL) = 0.0
    DIAMETER_MAC(L,NAL) = ALGAES(NAL).MINDIAMETER
      
    DO K = LAYERBOT(NAL,L),LAYERTOP(NAL,L)
      LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD*MAX(HPK(L,K), 0.1)            ! *** Handle thin layers (g C/m3)
            
      ! *** Allow for growth between layers, if allowed
      IF( WQV(L,K,19+NAL) > LAYER_THRESHOLD )THEN
        ! *** Height of macrophyte
        HEIGHT_MAC(L,NAL) = HEIGHT_MAC(L,NAL) + HPK(L,K)
        LayerRatio_MAC(L,K,NAL) = 1.0                                       ! *** Fraction of layer containing vegetation
        
        ! *** Redistribute mass to other layers
        IF( K < LAYERTOP(NAL,L) )THEN
          ! *** "Grow" any excess mass into the layer above
          TVAL1 = WQV(L,K,19+NAL) - LAYER_THRESHOLD
          TVAL1 = TVAL1*HPK(L,K)                                            ! *** Mass of macrophyte carbon to move/grow from layer
          TVAL1 = TVAL1*HPKI(L,K+1)                                         ! *** Concentration of macrophyte carbon to add/grow from layer below
          WQV(L,K,19+NAL)   = LAYER_THRESHOLD
          WQV(L,K+1,19+NAL) = WQV(L,K+1,19+NAL) + TVAL1

          WQVO(L,K,19+NAL)   = LAYER_THRESHOLD
          WQVO(L,K+1,19+NAL) = WQVO(L,K+1,19+NAL) + TVAL1
        ELSE
          ! *** Macrophyte is at its full depth but mass is larger than threshold.
          ! *** Therefore, increase stem diameter to account for mass
          THICK = Z(L,K) - Z(L,LAYERBOT(NAL,L)-1)
          IF( THICK <= 0. ) THICK = 1.0                                             ! *** Handle single layer models
          EXCESS_MASS = WQV(L,K,19+NAL) - LAYER_THRESHOLD
          
          ! *** Distribute excess mass
          DO K1 = LAYERBOT(NAL,L),K
            LAYER_THRESHOLD = ALGAES(NAL).THRESHOLD*MAX(HPK(L,K1), 0.1)             ! *** Handle thin layers
            IF( K1 < LAYERTOP(NAL,L) )THEN
              WQV(L,K1,19+NAL)  = WQV(L,K1,19+NAL) + EXCESS_MASS*DZC(L,K1)/THICK    ! *** Update layer concentrations
              WQVO(L,K1,19+NAL) = WQV(L,K1,19+NAL)                                  ! *** Update layer concentrations
            ENDIF
            DIAMETER_MAC(L,NAL) = DIAMETER_MAC(L,NAL)*EXCESS_MASS/LAYER_THRESHOLD/ALGAES(NAL).STEMDENSITY   ! *** Update stem diameter
            !IF( DIAMETER_MAC(L,NAL) < ALGAES(NAL).MINDIAMETER )THEN
            !  PRINT '(I6,2I5,F8.3,4F10.3)', NITER, NAL, L, HP(L), ALGAES(NAL).THRESHOLD, LAYER_THRESHOLD, DIAMETER_MAC(L,NAL), ALGAES(NAL).MINDIAMETER  ! DELME
            !ENDIF
            DIAMETER_MAC(L,NAL) = MAX(DIAMETER_MAC(L,NAL),ALGAES(NAL).MINDIAMETER)
          ENDDO
          EXIT
        ENDIF
      
      ELSE
        ! *** Partial penetration of current layer
        HEIGHT_MAC(L,NAL) = HEIGHT_MAC(L,NAL) + WQV(L,K,19+NAL)/LAYER_THRESHOLD/ALGAES(NAL).STEMDENSITY    ! *** Height of macrophyte
        LayerRatio_MAC(L,K,NAL) = HEIGHT_MAC(L,NAL)*HPKI(L,K)                                              ! *** Fraction of layer containing vegetation  
        IF( LayerRatio_MAC(L,K,NAL) < 0.01 ) LayerRatio_MAC(L,K,NAL) = 0.0                                 ! *** Ignore very small heights for hydrodynamics
      ENDIF

    ENDDO
          
    ! *** Adjust maximum layer height based on actual stem height
    IF( HP(L)*(Z(L,LAYERTOP(NAL,L)) - Z(L,LAYERBOT(NAL,L)-1)) > HEIGHT_MAC(L,NAL) )THEN
      THICK = 0.0
      DO K = LAYERBOT(NAL,L),LAYERTOP(NAL,L)
        THICK = THICK + HPK(L,K)
        IF( HEIGHT_MAC(L,NAL) < THICK )THEN
          LAYERTOP(NAL,L) = K
          EXIT
        ENDIF
      ENDDO
    ENDIF
          
    IF( HEIGHT_MAC(L,NAL) > 0.0 ) MAC_CELL(L) = .TRUE.                                                     ! *** Total vegetation height is used as a flag to compute turbulence
    
          ! *** ACCUMULATE ACTIVE VEGETATION LAYERS
          !VEGK(L) = VEGK(L) + HVGTC*HPKI(L,K)

  ENDDO        ! *** End of L index loop  
    
END SUBROUTINE Macro_Veg
  
END MODULE
  
