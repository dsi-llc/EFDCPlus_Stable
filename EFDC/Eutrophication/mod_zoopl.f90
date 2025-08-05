! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
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
!> @author Duc Kien Tran, Paul M Craig
!> @date 02/2020
!---------------------------------------------------------------------------!
  
  use GLOBAL    
  use INFOMOD
  use JULIANMOD
  use Variables_MPI
  use Broadcast_Routines
  use Variables_WQ
  
  implicit none
  
  contains

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine ZOOPL_JNP
  !
  !> @details Reading JSON input file for zooplankton
  !---------------------------------------------------------------------------!
  !> @author Duc Kien Tran
  !---------------------------------------------------------------------------!
  SUBROUTINE ZOOPL_CONTROL
    use fson
    use mod_fson_value, only: fson_value_count, fson_value_get
    
    type(fson_value), pointer :: json_data, zoogroups, pointsource, item
    
    integer :: M,NAL,L,K,NS
    integer :: I,J,ITMP,KK, NT
    real :: TEMPOZ,WQTT 
    
    if( process_id == master_id )then
      ! *** Parse the json file
      json_data => fson_parse("wq_zooplankton.jnp")
      Write(*,'(A)') 'WQ: READING WQ_ZOOPLANKTON.JNP - ZOOPLANKTON CONTROL FILE'     
      ! *** Get the wq_zooplankton group's parameters as an array
      call fson_get(json_data, "groups", zoogroups)
      
      ! *** Loop through each array item
      do NZO = 1, fson_value_count(zoogroups)
        ! *** Get the array item 
        item => fson_value_get(zoogroups, NZO)
        
        ! *** Lookup the values from the array
        call fson_get(item, "index", ZOOPL(NZO).IDZ)
        call fson_get(item, "predator_flag", ZOOPL(NZO).ISPREDATOR)
        call fson_get(item, "prey_flag", ZOOPL(NZO).ISPREY)
        
        ! *** General parameters
        call fson_get(item, "stoichiometry.nitrogen_to_carbon_ratio", ZOOPL(NZO).ANCZ)
        call fson_get(item, "stoichiometry.phosphorus_to_carbon_ratio", ZOOPL(NZO).APCZ)
        call fson_get(item, "stoichiometry.silica_to_carbon_ratio", ZOOPL(NZO).ASCZ)
        
        ! *** Grazing parameters
        call fson_get(item, "grazing.carbon_threshold", ZOOPL(NZO).CTZ)
        call fson_get(item, "grazing.halved_grazing_prey_density", ZOOPL(NZO).KHCZ)
        call fson_get(item, "grazing.maximum_ration", ZOOPL(NZO).RMAXZ(1))
        call fson_get(item, "grazing.utilization_of_LPOC", ZOOPL(NZO).ULZ)
        call fson_get(item, "grazing.utilization_of_RPOC", ZOOPL(NZO).URZ)
        call fson_get(item, "grazing.utilization_of_DOC", ZOOPL(NZO).UDZ)

        call fson_get(item, "grazing.optimal_grazing.lower_temperature", ZOOPL(NZO).TMZG1)
        call fson_get(item, "grazing.optimal_grazing.upper_temperature", ZOOPL(NZO).TMZG2)
        call fson_get(item, "grazing.optimal_grazing.lower_coefficient", ZOOPL(NZO).KTGZ1)
        call fson_get(item, "grazing.optimal_grazing.upper_coefficient", ZOOPL(NZO).KTGZ2)
        call fson_get(item, "grazing.utilization_of_ZOOPL", ZOOPL(NZO).UZPL)
        
        ! *** Basal Metabolism parameters
        call fson_get(item, "metabolism.reference_temperature", ZOOPL(NZO).TRZB)
        call fson_get(item, "metabolism.reference_rate", ZOOPL(NZO).BMRZ(1))
        call fson_get(item, "metabolism.effect_of_temperature", ZOOPL(NZO).KTBZ)
        call fson_get(item, "metabolism.fraction_of_phosphorus.DOP", ZOOPL(NZO).FPDBZ)   
        call fson_get(item, "metabolism.fraction_of_phosphorus.PO4", ZOOPL(NZO).FPIBZ)
        call fson_get(item, "metabolism.fraction_of_nitrogen.DON", ZOOPL(NZO).FNDBZ)
        call fson_get(item, "metabolism.fraction_of_nitrogen.NH4", ZOOPL(NZO).FNIBZ) 
        
        ! *** Predation parameters  
        call fson_get(item, "predation.reference_temperature", ZOOPL(NZO).TRZP)
        call fson_get(item, "predation.reference_rate", ZOOPL(NZO).PRRZ(1))
        call fson_get(item, "predation.effect_of_temperature", ZOOPL(NZO).KTPZ)
        call fson_get(item, "predation.fraction_of_carbon.LPOC", ZOOPL(NZO).FCLPZ)
        call fson_get(item, "predation.fraction_of_carbon.RPOC", ZOOPL(NZO).FCRPZ)
        call fson_get(item, "predation.fraction_of_carbon.DOC", ZOOPL(NZO).FCDPZ)
        call fson_get(item, "predation.fraction_of_phosphorus.LPOP", ZOOPL(NZO).FPLPZ)
        call fson_get(item, "predation.fraction_of_phosphorus.RPOP", ZOOPL(NZO).FPRPZ)
        call fson_get(item, "predation.fraction_of_phosphorus.DOP", ZOOPL(NZO).FPDPZ)
        call fson_get(item, "predation.fraction_of_phosphorus.PO4", ZOOPL(NZO).FPIPZ)
        call fson_get(item, "predation.fraction_of_nitrogen.LPON", ZOOPL(NZO).FNLPZ)
        call fson_get(item, "predation.fraction_of_nitrogen.RPON", ZOOPL(NZO).FNRPZ)
        call fson_get(item, "predation.fraction_of_nitrogen.DON", ZOOPL(NZO).FNDPZ)
        call fson_get(item, "predation.fraction_of_nitrogen.NH4", ZOOPL(NZO).FNIPZ)
        call fson_get(item, "predation.fraction_of_silica.SU", ZOOPL(NZO).FSPPZ)
        call fson_get(item, "predation.fraction_of_silica.SU", ZOOPL(NZO).FSAPZ)
      
        ! *** Death parameters     
        call fson_get(item, "death.critical_DO", ZOOPL(NZO).DOCRIT)
        call fson_get(item, "death.zero_DO_rate", ZOOPL(NZO).DZEROZ)
        call fson_get(item, "death.fraction_of_carbon.LPOC", ZOOPL(NZO).FCLDZ)
        call fson_get(item, "death.fraction_of_carbon.RPOC", ZOOPL(NZO).FCRDZ)
        call fson_get(item, "death.fraction_of_carbon.DOC", ZOOPL(NZO).FCDDZ)
        call fson_get(item, "death.fraction_of_phosphorus.LPOP", ZOOPL(NZO).FPLDZ)
        call fson_get(item, "death.fraction_of_phosphorus.RPOP", ZOOPL(NZO).FPRDZ)
        call fson_get(item, "death.fraction_of_phosphorus.DOP", ZOOPL(NZO).FPDDZ)
        call fson_get(item, "death.fraction_of_phosphorus.PO4", ZOOPL(NZO).FPIDZ)
        call fson_get(item, "death.fraction_of_nitrogen.LPON", ZOOPL(NZO).FNLDZ)
        call fson_get(item, "death.fraction_of_nitrogen.RPON", ZOOPL(NZO).FNRDZ)
        call fson_get(item, "death.fraction_of_nitrogen.DON", ZOOPL(NZO).FNDDZ)
        call fson_get(item, "death.fraction_of_nitrogen.NH4", ZOOPL(NZO).FNIDZ)
        call fson_get(item, "death.fraction_of_silica.SU", ZOOPL(NZO).FSPDZ)
        call fson_get(item, "death.fraction_of_silica.SA", ZOOPL(NZO).FSADZ)
        
        if( NALGAE > 0) call fson_get(item, "grazing.utilization_of_algae", ZOOPL(NZO).UBZ) 
      enddo
      
      ! *** Reading the kinetics file for spatial zones Zooplankton parameters 
      if( NWQZ > 1 )then
        write(*,'(A)')' WQ: WQ_BIO_ZOO_ZN.JNP'
        json_data => fson_parse("wq_bio_zoo_zn.jnp")  
        call fson_get(json_data, "groups", zoogroups)
    
        ! *** Loop through each array item
        do NZO = 1, fson_value_count(zoogroups)
          ! *** Get the array item 
          item => fson_value_get(zoogroups, NZO)
          
          call fson_get(item, "zooplankton_dynamic.max_grazing_ration", ZOOPL(NZO).RMAXZ)
          call fson_get(item, "zooplankton_dynamic.metabolism_rate", ZOOPL(NZO).BMRZ)
          call fson_get(item, "zooplankton_dynamic.predation_rate", ZOOPL(NZO).PRRZ)
        enddo
      endif
      
    endif ! *** End of master_id block
    
    ! *** Broadcast zooplankton global settings, not dependent on grid.
 
    Do NZO = 1, NZOOPL
     call Broadcast_Scalar(ZOOPL(NZO).IDZ,        master_id)
     call Broadcast_Scalar(ZOOPL(NZO).ISPREDATOR, master_id)
     call Broadcast_Scalar(ZOOPL(NZO).ISPREY,     master_id)
     
     call Broadcast_Scalar(ZOOPL(NZO).ANCZ,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).APCZ,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).ASCZ,       master_id)
     
     call Broadcast_Scalar(ZOOPL(NZO).CTZ,        master_id)
     call Broadcast_Scalar(ZOOPL(NZO).KHCZ,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).RMAXZ(1),   master_id)
     call Broadcast_Scalar(ZOOPL(NZO).ULZ,        master_id)
     call Broadcast_Scalar(ZOOPL(NZO).URZ,        master_id)
     call Broadcast_Scalar(ZOOPL(NZO).UDZ,        master_id)

     call Broadcast_Scalar(ZOOPL(NZO).TMZG1,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).TMZG2,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).KTGZ1,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).KTGZ2,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).UZPL,       master_id)
          
     call Broadcast_Scalar(ZOOPL(NZO).TRZB,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).BMRZ(1),    master_id)
     call Broadcast_Scalar(ZOOPL(NZO).KTBZ,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPDBZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPIBZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNDBZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNIBZ,      master_id)
     
     call Broadcast_Scalar(ZOOPL(NZO).TRZP,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).PRRZ(1),    master_id)
     call Broadcast_Scalar(ZOOPL(NZO).KTPZ,       master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FCLPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FCRPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FCDPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPLPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPRPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPDPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPIPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNLPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNRPZ,      master_id)    
     call Broadcast_Scalar(ZOOPL(NZO).FNDPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNIPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FSPPZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FSAPZ,      master_id)
     
     call Broadcast_Scalar(ZOOPL(NZO).DOCRIT,     master_id)
     call Broadcast_Scalar(ZOOPL(NZO).DZEROZ,     master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FCLDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FCRDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FCDDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPLDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPRDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPDDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FPIDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNLDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNRDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNDDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FNIDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FSPDZ,      master_id)
     call Broadcast_Scalar(ZOOPL(NZO).FSADZ,      master_id)
     
     call Broadcast_Array (ZOOPL(NZO).UBZ,        master_id)
     
    Enddo
    
    ! *** Zonal parameters
    if( NWQZ > 1 )then
      do NZO = 1, NZOOPL
        call Broadcast_Array(ZOOPL(NZO).RMAXZ,    master_id)
        call Broadcast_Array(ZOOPL(NZO).BMRZ,     master_id)
        call Broadcast_Array(ZOOPL(NZO).PRRZ,     master_id)
      enddo
    endif
    
    ! *** Set up look-up table for temperature dependency over -10 °C to 60 °C
    WQTDZMIN = -10
    WQTDZMAX =  60
    WTEMPZ = WQTDZMIN
    WQTDZINC = (WQTDZMAX-WQTDZMIN)/NWQTD
    
    do M = 1,NWQTD
      ! *** Loop for all zooplankton groups
      do NZO = 1,NZOOPL
        ! *** Reference temperature for zooplankton grazing
        WQTDGZ(M,NZO) = 1.
        if( WTEMPZ < ZOOPL(NZO).TMZG1 )then
          WQTDGZ(M,NZO) = EXP(-ZOOPL(NZO).KTGZ1 * (WTEMPZ - ZOOPL(NZO).TMZG1) * (WTEMPZ - ZOOPL(NZO).TMZG1))
        endif
        if( WTEMPZ > ZOOPL(NZO).TMZG2 )then
          WQTDGZ(M,NZO) = EXP(-ZOOPL(NZO).KTGZ2 * (WTEMPZ - ZOOPL(NZO).TMZG2) * (WTEMPZ - ZOOPL(NZO).TMZG2))
        endif
      
        ! *** Reference temperature for zooplankton metabolism
        WQTDBZ(M,NZO) = EXP(ZOOPL(NZO).KTBZ * (WTEMPZ - ZOOPL(NZO).TRZB))
        ! *** Reference temperature for zooplankton predation
        WQTDPZ(M,NZO) = EXP(ZOOPL(NZO).KTPZ * (WTEMPZ - ZOOPL(NZO).TRZP))
      
      enddo ! *** End loop for zooplankton
      WTEMPZ = WTEMPZ + WQTDZINC
    enddo    
    
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
    
    integer :: L, K, ND, LF, LL, LP, NAL, IZ
    real :: DOREF
    real :: WQACZ, WQRCZ
    real :: BZPAL
    real :: UBZT        !< Total utilization of algaes by zooplankton
    real :: ZPLPREY     !< Total prey available of zooplankton prey
    real :: LPOCZ       !< Effect of zooplankton on LPOC
    real :: RPOCZ       !< Effect of zooplankton on RPOC
    real :: DOCZ        !< Effect of zooplankton on DOC
    real :: LPONZ       !< Effect of zooplankton on LPON
    real :: RPONZ       !< Effect of zooplankton on RPON
    real :: DONZ        !< Effect of zooplankton on DON
    real :: NH4Z        !< Effect of zooplankton on NH4
    real :: LPOPZ       !< Effect of zooplankton on LPOP
    real :: RPOPZ       !< Effect of zooplankton on RPOP
    real :: DOPZ        !< Effect of zooplankton on DOP
    real :: PO4Z        !< Effect of zooplankton on PO4
    real :: SUZ         !< Effect of zooplankton on SU
    real :: SAZ         !< Effect of zooplankton on SA
    real :: DOZ         !< Effect of zooplankton on DO
    
    ! *** Zooplankton time step
    DTWQZO2 = 0.5*DTWQ
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, K, NAL, NZO)  &
    !$OMP             PRIVATE(WQACZ, WQRCZ, ZPLPREY, UBZT, DOREF, BZPAL)       &
    !$OMP             PRIVATE(LPOCZ, LPOPZ, LPONZ, RPOCZ, RPOPZ, RPONZ)  &
    !$OMP             PRIVATE(DOCZ,  DOPZ,  PO4Z,  DONZ, NH4Z, SUZ, SAZ, DOZ)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      do K = KC,1,-1         
        do LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)           
          ! *** Find an index for look-up table for temperature dependency
          IWQZT(L) = NINT((TWQ(L) - WQTDZMIN)/WQTDZINC) + 1  
          if( IWQZT(L) < 1 .or. IWQZT(L) > NWQTD  )then
            IWQZT(L) = MAX(IWQZT(L),1)  
            IWQZT(L) = MIN(IWQZT(L),NWQTD) 
          endif
        enddo
        
        ! *** Begin horizontal loop for zooplankton parameters
        do LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)
          
          if( LMASKDRY(L) )then
            
            ! *** Computation of prey available from prey zooplankton (gC/m3)
            ZPLPREY = 0.
            do NZO = 1, NZOOPL
              if( ZOOPL(NZO).ISPREY == 1 )then
                ZPLPREY = ZPLPREY + MAX(0., WQV(L,K,NWQVZ+NZO) - ZOOPL(NZO).CTZ)
              endif
            enddo
            
            do NZO = 1, NZOOPL         
              ! *** Computation of prey available from carbon pools to zooplankton (gC/m3)
              UBZT = 0.
              do NAL = 1, NALGAE  ! *** Loop for all algal groups except macroalgal
                if( ALGAES(NAL).ISMOBILE )then
                  BAZ(L,NZO,NAL) = MAX(0., WQV(L,K,19+NAL) - ZOOPL(NZO).CTZ)
                  UBZT = UBZT + ZOOPL(NZO).UBZ(NAL)*BAZ(L,NZO,NAL)
                endif
              enddo
            
              DOCAZ (L,NZO) = MAX(0., WQV(L,K,IDOC) - ZOOPL(NZO).CTZ)
              RPOCAZ(L,NZO) = MAX(0., WQV(L,K,IROC) - ZOOPL(NZO).CTZ)
              LPOCAZ(L,NZO) = MAX(0., WQV(L,K,ILOC) - ZOOPL(NZO).CTZ)
              ! *** Prey available
              PRAZ(L,NZO) = UBZT + ZOOPL(NZO).URZ*RPOCAZ(L,NZO) &
                          + ZOOPL(NZO).ULZ*LPOCAZ(L,NZO) + ZOOPL(NZO).UDZ*DOCAZ(L,NZO)
                            
              if( ZOOPL(NZO).ISPREDATOR == 1 )then
                PRAZ(L,NZO) = PRAZ(L,NZO) + ZOOPL(NZO).UZPL * ZPLPREY
              endif
            
              ! *** Computation of zooplankton growth rate (1/day)
              ZOOPL(NZO).WQGZ(L) = ZOOPL(NZO).RMAXZ(IMWQZT(L))*WQTDGZ(IWQZT(L),NZO)*PRAZ(L,NZO)/(PRAZ(L,NZO) + ZOOPL(NZO).KHCZ + 1.E-18)
              
              ! *** Computation of zooplankton metabolism (1/day)
              ZOOPL(NZO).WQBZ(L) = ZOOPL(NZO).BMRZ(IMWQZT(L))*WQTDBZ(IWQZT(L),NZO)
              
              ! *** Computation of zooplankton predation (1/day)
              ZOOPL(NZO).WQPZ(L) = ZOOPL(NZO).PRRZ(IMWQZT(L))*WQTDPZ(IWQZT(L),NZO)
              
              ! *** Computation of zooplankton death (1/day)
              DOREF = MIN(ZOOPL(NZO).DOCRIT, WQV(L,K,IDOX))
              ZOOPL(NZO).WQDZ(L) = ZOOPL(NZO).DZEROZ*(1.0 - DOREF/ZOOPL(NZO).DOCRIT)
              
            enddo ! *** End of loop for zooplankton groups
          else ! *** Dry cell bypass
            do NZO = 1,NZOOPL
              RPOCAZ(L,NZO) = 0.
              LPOCAZ(L,NZO) = 0.
              DOCAZ (L,NZO) = 0.
              PRAZ  (L,NZO) = 0.
              ZOOPL(NZO).WQGZ(L) = 0.
              ZOOPL(NZO).WQBZ(L) = 0.
              ZOOPL(NZO).WQPZ(L) = 0.
              ZOOPL(NZO).WQDZ(L) = 0.
            enddo ! *** End of loop for zooplankton groups
          endif
        enddo ! *** End horizontal loop for zooplankton parameters
        
        ! *** Computation of kinetics for each zooplankton group
        do NZO = 1, NZOOPL
          do LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)          
            ! ***              grazing          basal metab               death             predation
            WQACZ = (ZOOPL(NZO).WQGZ(L) - ZOOPL(NZO).WQBZ(L) - ZOOPL(NZO).WQDZ(L) - ZOOPL(NZO).WQPZ(L))*DTWQZO2
            WQZKK(L) = 1.0 / (1.0 - WQACZ)
            ! ***   point source    volume
            WQRCZ = WQWPSZ(L,K,NZO)*VOLWQ(L)
            WQZRR(L) = WQV(L,K,NWQVZ+NZO) + WQACZ*WQV(L,K,NWQVZ+NZO) + DTWQ*WQRCZ
            
            WQV(L,K,NWQVZ+NZO) = WQZRR(L)*WQZKK(L)
          enddo        
        enddo
        
        ! *** Effect of zooplankton on Phytoplankton
        do NAL = 1,NALGAE
          if( ALGAES(NAL).ISMOBILE )then 
            do LP = 1, LLWET(K,ND)
              L = LKWET(LP,K,ND)
              SBZPAL(L,K,NAL) = 0.
              do NZO = 1, NZOOPL 
                BZPAL = ZOOPL(NZO).WQGZ(L)*ZOOPL(NZO).UBZ(NAL)*BAZ(L,NZO,NAL)/(PRAZ(L,NZO) + 1.E-18)*WQV(L,K,NWQVZ+NZO)
                SBZPAL(L,K,NAL) = SBZPAL(L,K,NAL) + BZPAL
              enddo
            enddo
          endif
        enddo
        
        ! *** Effect of zooplankton on Carbon
        do NZO = 1,NZOOPL
          do LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            FRLP(L,NZO)  = ZOOPL(NZO).ULZ*LPOCAZ(L,NZO)/(PRAZ(L,NZO) + 1.E-18)
            FRRP(L,NZO)  = ZOOPL(NZO).URZ*RPOCAZ(L,NZO)/(PRAZ(L,NZO) + 1.E-18)
          enddo          
        enddo
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SRPOCZ(L,K) = 0.
          SLPOCZ(L,K) = 0.
          SDOCZ (L,K) = 0.
          do NZO = 1,NZOOPL
            ! *** RPOC
            if( ISKINETICS(IROC) == 1 )then 
              RPOCZ = (ZOOPL(NZO).FCRDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FCRPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRRP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)
              SRPOCZ(L,K) = SRPOCZ(L,K) + RPOCZ
            endif
            ! *** LPOC
            if( ISKINETICS(ILOC) == 1 )then 
              LPOCZ = (ZOOPL(NZO).FCLDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FCLPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRLP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)
              SLPOCZ(L,K) = SLPOCZ(L,K) + LPOCZ
            endif
            ! *** DOC
            if( ISKINETICS(IDOC) == 1 )then 
              DOCZ = (ZOOPL(NZO).FCDDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FCDPZ*ZOOPL(NZO).WQPZ(L))*WQV(L,K,NWQVZ+NZO)
              SDOCZ (L,K) = SDOCZ (L,K) + DOCZ
            endif
          enddo
        enddo
        
        ! *** Effect of zooplankton on Phosphorus
        do LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SRPOPZ(L,K) = 0.
          SLPOPZ(L,K) = 0.
          SDOPZ (L,K) = 0.
          SPO4Z (L,K) = 0.
          do NZO = 1,NZOOPL
            ! *** RPOP
            if( ISKINETICS(IROP) == 1 )then 
              RPOPZ = (ZOOPL(NZO).FPRDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPRPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRRP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SRPOPZ(L,K) = SRPOPZ(L,K) + RPOPZ
            endif
            ! *** LPOP
            if( ISKINETICS(ILOP) == 1 )then 
              LPOPZ = (ZOOPL(NZO).FPLDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPLPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRLP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SLPOPZ(L,K) = SLPOPZ(L,K) + LPOPZ
            endif
            ! *** DOP
            if( ISKINETICS(IDOP) == 1 )then
              DOPZ = (ZOOPL(NZO).FPDDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPDPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FPDBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SDOPZ(L,K) = SDOPZ (L,K) + DOPZ
            endif
            ! *** PO4
            if( ISKINETICS(IP4D) == 1 )then
              PO4Z = (ZOOPL(NZO).FPIDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FPIPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FPIBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).APCZ
              SPO4Z(L,K) = SPO4Z (L,K) + PO4Z
            endif
          enddo
        enddo 
        ! *** Effect of zooplankton on Nitrogen
        do LP = 1, LLWET(K,ND)
          L = LKWET(LP,K,ND)
          SRPONZ(L,K) = 0.
          SLPONZ(L,K) = 0.
          SDONZ (L,K) = 0.
          SNH4Z (L,K) = 0.
          do NZO = 1,NZOOPL
            !*** RPON
            if( ISKINETICS(IRON) == 1 )then
              RPONZ = (ZOOPL(NZO).FNRDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNRPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRRP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SRPONZ(L,K) = SRPONZ(L,K) + RPONZ
            endif
            !*** LPON
            if( ISKINETICS(ILON) == 1 )then
              LPONZ = (ZOOPL(NZO).FNLDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNLPZ*ZOOPL(NZO).WQPZ(L) &
                           - FRLP(L,NZO)*ZOOPL(NZO).WQGZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SLPONZ(L,K) = SLPONZ(L,K) + LPONZ
            endif
            ! *** DON
            if( ISKINETICS(IDON) == 1 )then
              DONZ = (ZOOPL(NZO).FNDDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNDPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FNDBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SDONZ(L,K) =  SDONZ(L,K) + DONZ
            endif
            ! *** NH4
            if( ISKINETICS(INHX) == 1 )then
              NH4Z = (ZOOPL(NZO).FNIDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FNIPZ*ZOOPL(NZO).WQPZ(L) &
                           + ZOOPL(NZO).FNIBZ*ZOOPL(NZO).WQBZ(L))*WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ANCZ
              SNH4Z(L,K) = SNH4Z(L,K) + NH4Z
            endif
          enddo
        enddo         
        
        ! *** Effect of zooplankton on Silica
        if( IWQSI == 1  )then
          do NZO = 1,NZOOPL
            do LP = 1, LLWET(K,ND)
              L = LKWET(LP,K,ND)
              FRSI(L,NZO) = 0.0
              do NAL = 1,NALGAE
                if( ALGAES(NAL).ISILICA > 0 ) FRSI(L,NZO) = FRSI(L,NZO) + ZOOPL(NZO).UBZ(NAL)*BAZ(L,NZO,NAL)/(PRAZ(L,NZO) + 1.E-18)
              enddo
            enddo          
          enddo
        endif
        ! *** Particulate Biogenic Silica
        if( ISKINETICS(ISUU) == 1 )then
          do LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SSUZ(L,K) = 0.
            do NZO = 1,NZOOPL
              SUZ = (ZOOPL(NZO).FSPDZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FSPPZ*ZOOPL(NZO).WQPZ(L)) &
                         *WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ASCZ*FRSI(L,NZO)
              SSUZ(L,K) = SSUZ(L,K) + SUZ
            enddo
          enddo
        endif
        ! *** Available Silica
        if( ISKINETICS(ISAA) == 1 )then      
          do LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SSAZ(L,K) = 0.
            do NZO = 1,NZOOPL
              SAZ = (ZOOPL(NZO).WQBZ(L) + ZOOPL(NZO).FSADZ*ZOOPL(NZO).WQDZ(L) + ZOOPL(NZO).FSAPZ*ZOOPL(NZO).WQPZ(L)) &
                         *WQV(L,K,NWQVZ+NZO)*ZOOPL(NZO).ASCZ*FRSI(L,NZO)
              SSAZ(L,K) = SSAZ(L,K) + SAZ
            enddo
          enddo
        endif
        
        ! *** Effect of zooplankton on Dissolved Oxygen
        if( ISKINETICS(IDOX) == 1 )then
          do LP = 1, LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SDOZ(L,K) = 0.
            do NZO = 1,NZOOPL
              DOZ = -ZOOPL(NZO).WQBZ(L)*WQV(L,K,NWQVZ+NZO)*WQAOCR   
              SDOZ(L,K) = SDOZ(L,K) + DOZ
            enddo
          enddo 
        endif
      enddo ! *** End of loop for layers
    enddo ! *** End of loop for domains
    !$OMP END PARALLEL DO

  END SUBROUTINE ZOOPL_KINETIC
  
END MODULE WQ_ZOOPLANKTON
