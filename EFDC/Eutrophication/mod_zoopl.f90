! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE WQ_ZOOPLANKTON

!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Module: Module name WQ_ZOOPLANKTON
!
!> @details Including zooplankton in the eutrophication model
!
!> @author Duc Kien TRAN
!> @date 02/2020
!---------------------------------------------------------------------------!
  
  Use GLOBAL    
  Use INFOMOD
  Use JULIANMOD
  Use Variables_MPI
  Use Broadcast_Routines
  Use Variables_WQ
  
  IMPLICIT NONE
  
  CONTAINS

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine ZOOPL_JNP
  !
  !> @details Reading JSON input file for zooplankton
  !---------------------------------------------------------------------------!
  !> @author D.K. TRAN
  !---------------------------------------------------------------------------!
  SUBROUTINE ZOOPL_CONTROL
    USE fson
    USE fson_value_m, ONLY: fson_value_count, fson_value_get
    
    TYPE(fson_value), POINTER :: json_data, zoogroups, pointsource, item
    
    INTEGER :: M,NAL,L,K,NS
    INTEGER :: I,J,ITMP,KK, NT
    REAL :: TEMPOZ,WQTT 
    
    IF( process_id == master_id )then
      ! *** Parse the json file
      json_data => fson_parse("wq_zooplankton.jnp")
      Write(*,'(A)') 'WQ: READING WQ_ZOOPLANKTON.JNP - ZOOPLANKTON CONTROL FILE'     
      ! *** Get the wq_zooplankton group's parameters as an array
      CALL fson_get(json_data, "groups", zoogroups)
      
      ! *** Loop through each array item
      DO NZO = 1, fson_value_count(zoogroups)
        ! *** Get the array item 
        item => fson_value_get(zoogroups, NZO)
        
        ! *** Lookup the values from the array
        CALL fson_get(item, "index", ZOOPL(NZO).IDZ)
        CALL fson_get(item, "predator_flag", ZOOPL(NZO).ISPREDATOR)
        CALL fson_get(item, "prey_flag", ZOOPL(NZO).ISPREY)
        
        ! *** General parameters
        CALL fson_get(item, "stoichiometry.nitrogen_to_carbon_ratio", ZOOPL(NZO).ANCZ)
        CALL fson_get(item, "stoichiometry.phosphorus_to_carbon_ratio", ZOOPL(NZO).APCZ)
        CALL fson_get(item, "stoichiometry.silica_to_carbon_ratio", ZOOPL(NZO).ASCZ)
        
        ! *** Grazing parameters
        CALL fson_get(item, "grazing.carbon_threshold", ZOOPL(NZO).CTZ)
        CALL fson_get(item, "grazing.halved_grazing_prey_density", ZOOPL(NZO).KHCZ)
        CALL fson_get(item, "grazing.maximum_ration", ZOOPL(NZO).RMAXZ(1))
        CALL fson_get(item, "grazing.utilization_of_LPOC", ZOOPL(NZO).ULZ)
        CALL fson_get(item, "grazing.utilization_of_RPOC", ZOOPL(NZO).URZ)
        CALL fson_get(item, "grazing.utilization_of_DOC", ZOOPL(NZO).UDZ)

        CALL fson_get(item, "grazing.optimal_grazing.lower_temperature", ZOOPL(NZO).TMZG1)
        CALL fson_get(item, "grazing.optimal_grazing.upper_temperature", ZOOPL(NZO).TMZG2)
        CALL fson_get(item, "grazing.optimal_grazing.lower_coefficient", ZOOPL(NZO).KTGZ1)
        CALL fson_get(item, "grazing.optimal_grazing.upper_coefficient", ZOOPL(NZO).KTGZ2)
        CALL fson_get(item, "grazing.utilization_of_ZOOPL", ZOOPL(NZO).UZPL)
        
        ! *** Basal Metabolism parameters
        CALL fson_get(item, "metabolism.reference_temperature", ZOOPL(NZO).TRZB)
        CALL fson_get(item, "metabolism.reference_rate", ZOOPL(NZO).BMRZ(1))
        CALL fson_get(item, "metabolism.effect_of_temperature", ZOOPL(NZO).KTBZ)
        CALL fson_get(item, "metabolism.fraction_of_phosphorus.DOP", ZOOPL(NZO).FPDBZ)   
        CALL fson_get(item, "metabolism.fraction_of_phosphorus.PO4", ZOOPL(NZO).FPIBZ)
        CALL fson_get(item, "metabolism.fraction_of_nitrogen.DON", ZOOPL(NZO).FNDBZ)
        CALL fson_get(item, "metabolism.fraction_of_nitrogen.NH4", ZOOPL(NZO).FNIBZ) 
        
        ! *** Predation parameters  
        CALL fson_get(item, "predation.reference_temperature", ZOOPL(NZO).TRZP)
        CALL fson_get(item, "predation.reference_rate", ZOOPL(NZO).PRRZ(1))
        CALL fson_get(item, "predation.effect_of_temperature", ZOOPL(NZO).KTPZ)
        CALL fson_get(item, "predation.fraction_of_carbon.LPOC", ZOOPL(NZO).FCLPZ)
        CALL fson_get(item, "predation.fraction_of_carbon.RPOC", ZOOPL(NZO).FCRPZ)
        CALL fson_get(item, "predation.fraction_of_carbon.DOC", ZOOPL(NZO).FCDPZ)
        CALL fson_get(item, "predation.fraction_of_phosphorus.LPOP", ZOOPL(NZO).FPLPZ)
        CALL fson_get(item, "predation.fraction_of_phosphorus.RPOP", ZOOPL(NZO).FPRPZ)
        CALL fson_get(item, "predation.fraction_of_phosphorus.DOP", ZOOPL(NZO).FPDPZ)
        CALL fson_get(item, "predation.fraction_of_phosphorus.PO4", ZOOPL(NZO).FPIPZ)
        CALL fson_get(item, "predation.fraction_of_nitrogen.LPON", ZOOPL(NZO).FNLPZ)
        CALL fson_get(item, "predation.fraction_of_nitrogen.RPON", ZOOPL(NZO).FNRPZ)
        CALL fson_get(item, "predation.fraction_of_nitrogen.DON", ZOOPL(NZO).FNDPZ)
        CALL fson_get(item, "predation.fraction_of_nitrogen.NH4", ZOOPL(NZO).FNIPZ)
        CALL fson_get(item, "predation.fraction_of_silica.SU", ZOOPL(NZO).FSPPZ)
        CALL fson_get(item, "predation.fraction_of_silica.SU", ZOOPL(NZO).FSAPZ)
      
        ! *** Death parameters     
        CALL fson_get(item, "death.critical_DO", ZOOPL(NZO).DOCRIT)
        CALL fson_get(item, "death.zero_DO_rate", ZOOPL(NZO).DZEROZ)
        CALL fson_get(item, "death.fraction_of_carbon.LPOC", ZOOPL(NZO).FCLDZ)
        CALL fson_get(item, "death.fraction_of_carbon.RPOC", ZOOPL(NZO).FCRDZ)
        CALL fson_get(item, "death.fraction_of_carbon.DOC", ZOOPL(NZO).FCDDZ)
        CALL fson_get(item, "death.fraction_of_phosphorus.LPOP", ZOOPL(NZO).FPLDZ)
        CALL fson_get(item, "death.fraction_of_phosphorus.RPOP", ZOOPL(NZO).FPRDZ)
        CALL fson_get(item, "death.fraction_of_phosphorus.DOP", ZOOPL(NZO).FPDDZ)
        CALL fson_get(item, "death.fraction_of_phosphorus.PO4", ZOOPL(NZO).FPIDZ)
        CALL fson_get(item, "death.fraction_of_nitrogen.LPON", ZOOPL(NZO).FNLDZ)
        CALL fson_get(item, "death.fraction_of_nitrogen.RPON", ZOOPL(NZO).FNRDZ)
        CALL fson_get(item, "death.fraction_of_nitrogen.DON", ZOOPL(NZO).FNDDZ)
        CALL fson_get(item, "death.fraction_of_nitrogen.NH4", ZOOPL(NZO).FNIDZ)
        CALL fson_get(item, "death.fraction_of_silica.SU", ZOOPL(NZO).FSPDZ)
        CALL fson_get(item, "death.fraction_of_silica.SA", ZOOPL(NZO).FSADZ)
        
        If (NALGAE > 0) CALL fson_get(item, "grazing.utilization_of_algae", ZOOPL(NZO).UBZ) 
      ENDDO
      
      ! *** Reading the kinetics file for spatial zones Zooplankton parameters 
      IF( NWQZ > 1 )THEN
        WRITE(*,'(A)')' WQ: WQ_BIO_ZOO_ZN.JNP'
        json_data => fson_parse("wq_bio_zoo_zn.jnp")  
        CALL fson_get(json_data, "groups", zoogroups)
    
        ! *** Loop through each array item
        DO NZO = 1, fson_value_count(zoogroups)
          ! *** Get the array item 
          item => fson_value_get(zoogroups, NZO)
          
          Call fson_get(item, "zooplankton_dynamic.max_grazing_ration", ZOOPL(NZO).RMAXZ)
          Call fson_get(item, "zooplankton_dynamic.metabolism_rate", ZOOPL(NZO).BMRZ)
          Call fson_get(item, "zooplankton_dynamic.predation_rate", ZOOPL(NZO).PRRZ)
        ENDDO
      END IF
      
    END IF ! *** End of master_id block
    
    ! *** Broadcast zooplankton global settings, not dependent on grid.
 
    Do NZO = 1, NZOOPL
     Call Broadcast_Scalar(ZOOPL(NZO).IDZ,        master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).ISPREDATOR, master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).ISPREY,     master_id)
     
     Call Broadcast_Scalar(ZOOPL(NZO).ANCZ,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).APCZ,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).ASCZ,       master_id)
     
     Call Broadcast_Scalar(ZOOPL(NZO).CTZ,        master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).KHCZ,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).RMAXZ(1),   master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).ULZ,        master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).URZ,        master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).UDZ,        master_id)

     Call Broadcast_Scalar(ZOOPL(NZO).TMZG1,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).TMZG2,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).KTGZ1,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).KTGZ2,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).UZPL,       master_id)
          
     Call Broadcast_Scalar(ZOOPL(NZO).TRZB,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).BMRZ(1),    master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).KTBZ,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPDBZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPIBZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNDBZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNIBZ,      master_id)
     
     Call Broadcast_Scalar(ZOOPL(NZO).TRZP,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).PRRZ(1),    master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).KTPZ,       master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FCLPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FCRPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FCDPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPLPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPRPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPDPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPIPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNLPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNRPZ,      master_id)    
     Call Broadcast_Scalar(ZOOPL(NZO).FNDPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNIPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FSPPZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FSAPZ,      master_id)
     
     Call Broadcast_Scalar(ZOOPL(NZO).DOCRIT,     master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).DZEROZ,     master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FCLDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FCRDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FCDDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPLDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPRDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPDDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FPIDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNLDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNRDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNDDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FNIDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FSPDZ,      master_id)
     Call Broadcast_Scalar(ZOOPL(NZO).FSADZ,      master_id)
     
     Call Broadcast_Array (ZOOPL(NZO).UBZ,        master_id)
     
    Enddo
    
    ! *** Zonal parameters
    IF( NWQZ > 1 )THEN
      DO NZO = 1, NZOOPL
        Call Broadcast_Array(ZOOPL(NZO).RMAXZ,    master_id)
        Call Broadcast_Array(ZOOPL(NZO).BMRZ,     master_id)
        Call Broadcast_Array(ZOOPL(NZO).PRRZ,     master_id)
      ENDDO
    END IF
    
    ! Set up look-up table for temperature dependency over -1O°C to 50°C
    WQTDZMIN =-10
    WQTDZMAX = 50
    WTEMPZ = WQTDZMIN
    WQTDZINC = (WQTDZMAX-WQTDZMIN)/NWQTD
    
    DO M = 1,NWQTD
      ! *** Loop for all zooplankton groups
      DO NZO = 1,NZOOPL
        ! *** Reference temperature for zooplankton grazing
        WQTDGZ(M,NZO) = 1.
        IF( WTEMPZ < ZOOPL(NZO).TMZG1 )THEN
          WQTDGZ(M,NZO) = EXP(-ZOOPL(NZO).KTGZ1*(WTEMPZ - ZOOPL(NZO).TMZG1)*(WTEMPZ - ZOOPL(NZO).TMZG1))
        ENDIF
        IF( WTEMPZ > ZOOPL(NZO).TMZG2 )THEN
          WQTDGZ(M,NZO) = EXP(-ZOOPL(NZO).KTGZ2*(WTEMPZ - ZOOPL(NZO).TMZG2)*(WTEMPZ - ZOOPL(NZO).TMZG2))
        ENDIF
      
        ! *** Reference temperature for zooplankton metabolism
        WQTDBZ(M,NZO) = EXP(-ZOOPL(NZO).KTBZ*(WTEMPZ - ZOOPL(NZO).TRZB))
        ! *** Reference temperature for zooplankton predation
        WQTDPZ(M,NZO) = EXP(-ZOOPL(NZO).KTPZ*(WTEMPZ - ZOOPL(NZO).TRZP))
      
      ENDDO ! *** End loop for zooplankton
      WTEMPZ = WTEMPZ + WQTDZINC
    ENDDO    
    
  END SUBROUTINE ZOOPL_CONTROL
    
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine ZOOPL_KINETIC
  !
  !> @details Computing the kinectics of zooplankton
  !> @details Interaction between zooplankton and eutrophication model
  !---------------------------------------------------------------------------!
  !> @author D.K. TRAN
  !---------------------------------------------------------------------------!
  SUBROUTINE ZOOPL_KINETIC
    
    INTEGER :: L, K, ND, LF, LL, LP, NAL, IZ
    REAL :: DOREF
    REAL :: WQACZ, WQRCZ
    REAL :: BZPAL
    REAL :: UBZT        !< Total utilization of algaes by zooplankton
    REAL :: ZPLPREY     !< Total prey available of zooplankton prey
    REAL :: LPOCZ       !< Effect of zooplankton on LPOC
    REAL :: RPOCZ       !< Effect of zooplankton on RPOC
    REAL :: DOCZ        !< Effect of zooplankton on DOC
    REAL :: LPONZ       !< Effect of zooplankton on LPON
    REAL :: RPONZ       !< Effect of zooplankton on RPON
    REAL :: DONZ        !< Effect of zooplankton on DON
    REAL :: NH4Z        !< Effect of zooplankton on NH4
    REAL :: LPOPZ       !< Effect of zooplankton on LPOP
    REAL :: RPOPZ       !< Effect of zooplankton on RPOP
    REAL :: DOPZ        !< Effect of zooplankton on DOP
    REAL :: PO4Z        !< Effect of zooplankton on PO4
    REAL :: SUZ         !< Effect of zooplankton on SU
    REAL :: SAZ         !< Effect of zooplankton on SA
    REAL :: DOZ         !< Effect of zooplankton on DO
    
    ! *** Zooplankton time step
    DTWQZO2 = 0.5*DTWQ
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, K, NAL, NZO)  &
    !$OMP             PRIVATE(WQACZ, WQRCZ, ZPLPREY, UBZT, DOREF, BZPAL)       &
    !$OMP             PRIVATE(LPOCZ, LPOPZ, LPONZ, RPOCZ, RPOPZ, RPONZ)  &
    !$OMP             PRIVATE(DOCZ,  DOPZ,  PO4Z,  DONZ, NH4Z, SUZ, SAZ, DOZ)
    DO ND = 1,NDM
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      DO K = KC,1,-1         
        DO LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)           
          ! *** Find an index for look-up table for temperature dependency
          IWQZT(L) = NINT((TWQ(L) - WQTDZMIN)/WQTDZINC) + 1  
          IF( IWQZT(L) < 1 .OR. IWQZT(L) > NWQTD  )THEN
            IWQZT(L) = MAX(IWQZT(L),1)  
            IWQZT(L) = MIN(IWQZT(L),NWQTD) 
          ENDIF
        ENDDO
        
        ! *** Begin horizontal loop for zooplankton parameters
        DO LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)
          
          IF( LMASKDRY(L) )THEN
            
            ! *** Computation of prey available from prey zooplankton (gC/m3)
            ZPLPREY = 0.
            DO NZO = 1, NZOOPL
              IF( ZOOPL(NZO).ISPREY == 1 )THEN
                ZPLPREY = ZPLPREY + MAX(0., WQV(L,K,NWQVZ+NZO) - ZOOPL(NZO).CTZ)
              END IF
            END DO
            
            DO NZO = 1, NZOOPL         
              ! *** Computation of prey available from carbon pools to zooplankton (gC/m3)
              UBZT = 0.
              DO NAL = 1, NALGAE  ! *** Loop for all algal groups except macroalgal
                IF( ALGAES(NAL).ISMOBILE )THEN
                  BAZ(L,NZO,NAL) = MAX(0., WQV(L,K,19+NAL) - ZOOPL(NZO).CTZ)
                  UBZT = UBZT + ZOOPL(NZO).UBZ(NAL)*BAZ(L,NZO,NAL)
                ENDIF
              END DO
            
              DOCAZ (L,NZO) = MAX(0., WQV(L,K,IDOC) - ZOOPL(NZO).CTZ)
              RPOCAZ(L,NZO) = MAX(0., WQV(L,K,IROC) - ZOOPL(NZO).CTZ)
              LPOCAZ(L,NZO) = MAX(0., WQV(L,K,ILOC) - ZOOPL(NZO).CTZ)
              ! *** Prey available
              PRAZ(L,NZO) = UBZT + ZOOPL(NZO).URZ*RPOCAZ(L,NZO) &
                          + ZOOPL(NZO).ULZ*LPOCAZ(L,NZO) + ZOOPL(NZO).UDZ*DOCAZ(L,NZO)
                            
              IF( ZOOPL(NZO).ISPREDATOR == 1 )THEN
                PRAZ(L,NZO) = PRAZ(L,NZO) + ZOOPL(NZO).UZPL * ZPLPREY
              END IF
            
              ! *** Computation of zooplankton growth rate (1/day)
              ZOOPL(NZO).WQGZ(L) = ZOOPL(NZO).RMAXZ(IMWQZT(L))*WQTDGZ(IWQZT(L),NZO)*PRAZ(L,NZO)/(PRAZ(L,NZO) + ZOOPL(NZO).KHCZ + 1.E-18)
              
              ! *** Computation of zooplankton metabolism (1/day)
              ZOOPL(NZO).WQBZ(L) = ZOOPL(NZO).BMRZ(IMWQZT(L))*WQTDBZ(IWQZT(L),NZO)
              
              ! *** Computation of zooplankton predation (1/day)
              ZOOPL(NZO).WQPZ(L) = ZOOPL(NZO).PRRZ(IMWQZT(L))*WQTDPZ(IWQZT(L),NZO)
              
              ! *** Computation of zooplankton death (1/day)
              DOREF = MIN(ZOOPL(NZO).DOCRIT, WQV(L,K,IDOX))
              ZOOPL(NZO).WQDZ(L) = ZOOPL(NZO).DZEROZ*(1.0 - DOREF/ZOOPL(NZO).DOCRIT)
              
            END DO ! *** End of loop for zooplankton groups
          ELSE ! *** Dry cell bypass
            DO NZO = 1,NZOOPL
              RPOCAZ(L,NZO) = 0.
              LPOCAZ(L,NZO) = 0.
              DOCAZ (L,NZO) = 0.
              PRAZ  (L,NZO) = 0.
              ZOOPL(NZO).WQGZ(L) = 0.
              ZOOPL(NZO).WQBZ(L) = 0.
              ZOOPL(NZO).WQPZ(L) = 0.
              ZOOPL(NZO).WQDZ(L) = 0.
            ENDDO ! *** End of loop for zooplankton groups
          ENDIF
        ENDDO ! *** End horizontal loop for zooplankton parameters
        
        ! *** Computation of kinetics for each zooplankton group
        DO NZO = 1, NZOOPL
          DO LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)          
            ! ***              grazing          basal metab               death             predation
            WQACZ = (ZOOPL(NZO).WQGZ(L) - ZOOPL(NZO).WQBZ(L) - ZOOPL(NZO).WQDZ(L) - ZOOPL(NZO).WQPZ(L))*DTWQZO2
            WQZKK(L) = 1.0 / (1.0 - WQACZ)
            ! ***   point source    volume
            WQRCZ = WQWPSZ(L,K,NZO)*VOLWQ(L)
            WQZRR(L) = WQV(L,K,NWQVZ+NZO) + WQACZ*WQV(L,K,NWQVZ+NZO) + DTWQ*WQRCZ
            
            WQV(L,K,NWQVZ+NZO) = WQZRR(L)*WQZKK(L)
          ENDDO        
        ENDDO
        
        ! *** Effect of zooplankton on Phytoplankton
        DO NAL = 1,NALGAE
          IF( ALGAES(NAL).ISMOBILE )THEN 
            DO LP = 1, LLWET(K,ND)
              L = LKWET(LP,K,ND)
              SBZPAL(L,K,NAL) = 0.
              DO NZO = 1, NZOOPL 
                BZPAL = ZOOPL(NZO).WQGZ(L)*ZOOPL(NZO).UBZ(NAL)*BAZ(L,NZO,NAL)/(PRAZ(L,NZO) + 1.E-18)*WQV(L,K,NWQVZ+NZO)
                SBZPAL(L,K,NAL) = SBZPAL(L,K,NAL) + BZPAL
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        
        ! *** Effect of zooplankton on Carbon
        DO NZO = 1,NZOOPL
          DO LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            FRLP(L,NZO)  = ZOOPL(NZO).ULZ*LPOCAZ(L,NZO)/(PRAZ(L,NZO) + 1.E-18)
            FRRP(L,NZO)  = ZOOPL(NZO).URZ*RPOCAZ(L,NZO)/(PRAZ(L,NZO) + 1.E-18)
          ENDDO          
        ENDDO
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SRPOCZ(L,K) = 0.
          SLPOCZ(L,K) = 0.
          SDOCZ (L,K) = 0.
          DO NZO = 1,NZOOPL
            ! *** RPOC
            IF( ISKINETICS(IROC) == 1 )THEN 
              RPOCZ = (ZOOPL(NZO).FCRDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FCRPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRRP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)
              SRPOCZ(L,K) = SRPOCZ(L,K) + RPOCZ
            ENDIF
            ! *** LPOC
            IF( ISKINETICS(ILOC) == 1 )THEN 
              LPOCZ = (ZOOPL(NZO).FCLDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FCLPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRLP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)
              SLPOCZ(L,K) = SLPOCZ(L,K) + LPOCZ
            ENDIF
            ! *** DOC
            IF( ISKINETICS(IDOC) == 1 )THEN 
              DOCZ = (ZOOPL(NZO).FCDDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FCDPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)
              SDOCZ (L,K) = SDOCZ (L,K) + DOCZ
            ENDIF
          ENDDO
        ENDDO
        
        ! *** Effect of zooplankton on Phosphorus
        DO LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SRPOPZ(L,K) = 0.
          SLPOPZ(L,K) = 0.
          SDOPZ (L,K) = 0.
          SPO4Z (L,K) = 0.
          DO NZO = 1,NZOOPL
            ! *** RPOP
            IF( ISKINETICS(IROP) == 1 )THEN 
              RPOPZ = (ZOOPL(NZO).FPRDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPRPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRRP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SRPOPZ(L,K) = SRPOPZ(L,K) + RPOPZ
            ENDIF
            ! *** LPOP
            IF( ISKINETICS(ILOP) == 1 )THEN 
              LPOPZ = (ZOOPL(NZO).FPLDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPLPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRLP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SLPOPZ(L,K) = SLPOPZ(L,K) + LPOPZ
            ENDIF
            ! *** DOP
            IF( ISKINETICS(IDOP) == 1 )THEN
              DOPZ = (ZOOPL(NZO).FPDDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPDPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FPDBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SDOPZ(L,K) = SDOPZ (L,K) + DOPZ
            ENDIF
            ! *** PO4
            IF( ISKINETICS(IP4D) == 1 )THEN
              PO4Z = (ZOOPL(NZO).FPIDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPIPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FPIBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SPO4Z(L,K) = SPO4Z (L,K) + PO4Z
            ENDIF
          ENDDO
        ENDDO 
        ! *** Effect of zooplankton on Nitrogen
        DO LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SRPONZ(L,K) = 0.
          SLPONZ(L,K) = 0.
          SDONZ (L,K) = 0.
          SNH4Z (L,K) = 0.
          DO NZO = 1,NZOOPL
            !*** RPON
            IF( ISKINETICS(IRON) == 1 )THEN
              RPONZ = (ZOOPL(NZO).FNRDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNRPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRRP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SRPONZ(L,K) = SRPONZ(L,K) + RPONZ
            ENDIF
            !*** LPON
            IF( ISKINETICS(ILON) == 1 )THEN
              LPONZ = (ZOOPL(NZO).FNLDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNLPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRLP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SLPONZ(L,K) = SLPONZ(L,K) + LPONZ
            ENDIF
            ! *** DON
            IF( ISKINETICS(IDON) == 1 )THEN
              DONZ = (ZOOPL(NZO).FNDDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNDPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FNDBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SDONZ(L,K) =  SDONZ(L,K) + DONZ
            ENDIF
            ! *** NH4
            IF( ISKINETICS(INHX) == 1 )THEN
              NH4Z = (ZOOPL(NZO).FNIDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNIPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FNIBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SNH4Z(L,K) = SNH4Z(L,K) + NH4Z
            ENDIF
          ENDDO
        ENDDO         
        
        ! *** Effect of zooplankton on Silica
        IF( IWQSI == 1  )THEN
          DO NZO = 1,NZOOPL
            DO LP = 1, LLWET(K,ND)
              L = LKWET(LP,K,ND)
              FRSI(L,NZO) = 0.0
              DO NAL = 1,NALGAE
                IF( ALGAES(NAL).ISILICA > 0 ) FRSI(L,NZO) = FRSI(L,NZO) + ZOOPL(NZO).UBZ(NAL)*BAZ(L,NZO,NAL)/(PRAZ(L,NZO) + 1.E-18)
              ENDDO
            ENDDO          
          ENDDO
        ENDIF
        ! *** Particulate Biogenic Silica
        IF( ISKINETICS(ISUU) == 1 )THEN
          DO LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SSUZ(L,K) = 0.
            DO NZO = 1,NZOOPL
              SUZ = (ZOOPL(NZO).FSPDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FSPPZ*ZOOPL(NZO).WQPZ(L)) &
                         *WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ASCZ*FRSI(L,NZO)
              SSUZ(L,K) = SSUZ(L,K) + SUZ
            ENDDO
          ENDDO
        ENDIF
        ! *** Available Silica
        IF( ISKINETICS(ISAA) == 1 )THEN      
          DO LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SSAZ(L,K) = 0.
            DO NZO = 1,NZOOPL
              SAZ = (ZOOPL(NZO).WQBZ(L) + ZOOPL(NZO).FSADZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FSAPZ*ZOOPL(NZO).WQPZ(L)) &
                         *WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ASCZ*FRSI(L,NZO)
              SSAZ(L,K) = SSAZ(L,K) + SAZ
            ENDDO
          ENDDO
        ENDIF
        
        ! *** Effect of zooplankton on Dissolved Oxygen
        IF( ISKINETICS(IDOX) == 1 )THEN
          DO LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SDOZ(L,K) = 0.
            DO NZO = 1,NZOOPL
              DOZ = -ZOOPL(NZO).WQBZ(L)*WQV(L,K,NWQVZ+NZO)*WQAOCR   
              SDOZ(L,K) = SDOZ(L,K) + DOZ
            ENDDO
          ENDDO 
        ENDIF
      ENDDO ! *** End of loop for layers
    ENDDO ! *** End of loop for domains
    !$OMP END PARALLEL DO

  END SUBROUTINE ZOOPL_KINETIC
  
END MODULE WQ_ZOOPLANKTON
