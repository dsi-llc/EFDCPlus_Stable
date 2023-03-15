! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE SHELLFISHMOD

  Use GLOBAL    
  Use Variables_WQ

  IMPLICIT NONE

    TYPE SHELLFISHTYPE
      INTEGER :: IFILT                          !< Option for filtration model    
      INTEGER :: IFTEMP                         !< Option for temperature effect on filtration
      INTEGER :: IFSALT                         !< Option for salinity effect on filtration
      INTEGER :: IFTSS                          !< Option for sediment effect on filtration
      INTEGER :: IFDO                           !< Option for dissolved oxygen effect on filtration
      INTEGER :: IRES                           !< Option for respiration model
      INTEGER :: IRTEMP                         !< Option for temperature effect on respiration
      INTEGER :: IRSALT                         !< Option for salinity effect on respiration
      INTEGER :: IRTSS                          !< Option for sediment effect on respiration
      INTEGER :: IRDO                           !< Option for dissolved oxygen effect on respiration
      INTEGER :: ITDEATH                        !< Temperature dependent mortality flag
      INTEGER :: ISPAWN                         !< Option for spawning
      REAL(RKD) :: IC                           !< Initial biomass [g C]
      REAL(RKD) :: DRATIO                       !< Density ratio [-]
      REAL(RKD) :: CWLEN                        !< Factor for the length - weight relationship [-]
      REAL(RKD) :: PWLEN                        !< Exponent for the length - weight relationship [-]
      REAL(RKD) :: FRMAX                        !< Maximum filtration rate [-]
      REAL(RKD) :: TOPT                         !< Temperature for optimal filtration [deg. C]
      REAL(RKD) :: KTG                          !< Effect of temperature on filtration [-]
      REAL(RKD) :: KHSOY                        !< Salinity at which filtration rate is halved [ppt]
      REAL(RKD) :: DOHX                         !< DO at which filtration rate is one-half [g/m^3]
      REAL(RKD) :: DOQX                         !< DO at which filtration rate is one-fourth [g/m^3]
      REAL(RKD) :: K_ASSIM                      !< Assimilation efficiency [-]
      REAL(RKD) :: K_RESPI                      !< Basal metabolic rate [1/day]
      REAL(RKD) :: K_GROEF                      !< Growth efficiency [-]
      REAL(RKD) :: K_DEATH                      !< Specific mortality rate [1/day]
      REAL(RKD) :: TEMP_RESPI                   !< Reference temperature for respiration [deg. C]
      REAL(RKD) :: THETA_RESPI                  !< Temperature coefficient for respiration [-]
      REAL(RKD) :: TEMP_DEATH                   !< Reference temperature for mortality [deg. C]
      REAL(RKD) :: THETA_DEATH                  !< Temperature coefficient for shellfish mortality [-]
      REAL(RKD) :: DFROC                   !< Fraction of RPOC from dead oyster [-]
      REAL(RKD) :: DFLOC                   !< Fraction of LPOC from dead oyster [-]
      REAL(RKD) :: DFDOC                   !< Fraction of DOC from dead oyster [-]
      REAL(RKD) :: FFROC                   !< Fraction of RPOC in feces [-]
      REAL(RKD) :: FFLOC                   !< Fraction of LPOC in feces [-]
      REAL(RKD) :: FFDOC                   !< Fraction of DOC in feces [-]
      REAL(RKD) :: FFSUU                   !< Fraction of SUU in feces [-]
      REAL(RKD) :: UFROP                     !< Fraction of ROP in urine [-]
      REAL(RKD) :: UFLOP                     !< Fraction of LOP in urine [-]
      REAL(RKD) :: UFDOP                     !< Fraction of DOP in urine [-]
      REAL(RKD) :: UFP4D                     !< Fraction of P4D in urine [-]
      REAL(RKD) :: UFRON                     !< Fraction of RON in urine [-]
      REAL(RKD) :: UFLON                     !< Fraction of LON in urine [-]
      REAL(RKD) :: UFDON                     !< Fraction of DON in urine [-]
      REAL(RKD) :: UFNHX                     !< Fraction of NHX in urine [-]
      REAL(RKD) :: WQACF                        !< WQACF: Fraction of carbon in dry meat weight [-]
      REAL(RKD) :: WQANC                        !< WQANC: N/C ratio of shellfish [g N/g C]
      REAL(RKD) :: KACHC                        !< KACHC: Assimilation efficiency for cyanobacteria
      REAL(RKD) :: KACHD                        !< KACHD: Assimilation efficiency for diatom algae
      REAL(RKD) :: KACHG                        !< KACHG: Assimilation efficiency for green algae
      REAL(RKD) :: KAPOC                        !< KAPOC: Assimilation efficiency for POC      
      REAL(RKD) :: SPAWN_WEIGHT                 !< Min. dry weight for spawning [g C]
      REAL(RKD) :: SPAWN_TEMP                   !< Min. temperature for spawning [deg C]
      REAL(RKD) :: SPAWN_RATIO                  !< Fraction of cumulative reproductive biomass for spawing
      REAL(RKD) :: EGG_WEIGHT                   !< egg dry weight [g C]
      REAL(RKD) :: EGG_CAL                      !< Caloric content of egg [cal/g dry weight]
      REAL(RKD) :: SHF_CAL                      !< Caloric content of shellfish [cal/g dry weight]
      REAL(RKD) :: RESCAL                       !< Energy comsumption by respiration [cal/mL O2]
      REAL(RKD) :: EGG_HATCHED                  !< fraction of egg hatched to larvae
      REAL(RKD) :: EGG_TEMP1                    !< Min. temprature for spawned eggs to survive [deg. C]
      REAL(RKD) :: EGG_TEMP2                    !< Max. temprature for spawned eggs to survive [deg. C]
      REAL(RKD) :: EGG_SALT1                    !< Min. salinity for spawned eggs to survive [ppt]
      REAL(RKD) :: EGG_SALT2                    !< Max. salinity for spawned eggs to survive [ppt]
      REAL(RKD) :: LARVAE_RECRUITED             !< Fraction of larvae recruited per spawn [-]
      REAL(RKD) :: LARVAE_SURVIVAL              !< Fraction of larvae to spat survival [-]
      REAL(RKD) :: LARVAE_SPAN                  !< Larvae life span [days]
      REAL(RKD),ALLOCATABLE :: SPAWN(:,:)       !< Spawn
      REAL(RKD),ALLOCATABLE :: CRB(:,:)         !< Cumulative reproductive biomass [g C]
      REAL(RKD),ALLOCATABLE :: PR(:,:)          !< Reproductive tissue production
      REAL(RKD),ALLOCATABLE :: NP(:,:)          !< Net production
      REAL(RKD),ALLOCATABLE :: CELL_C(:,:)      !< Biomass/carbon weight per grid cell layer [g C]
      REAL(RKD),ALLOCATABLE :: INDI_C(:,:)      !< INDI_C(LCM,KC): carbon weight per individual in cell (Dry meat weight [g C])
      REAL(RKD),ALLOCATABLE :: SHLEN(:,:)       !< Shell length [cm]
      REAL(RKD),ALLOCATABLE :: INDI(:,:)        !< Number of shellfish individuals [#]
      REAL(RKD),ALLOCATABLE :: INDI_D(:,:)      !< Number of death individuals [#]
      REAL(RKD),ALLOCATABLE :: INDI_H(:,:)      !< Number of harvested individuals [#]
      REAL(RKD),ALLOCATABLE :: HARV_C(:,:)      !< Harvested biomass per grid cell layer [g C]
      REAL(RKD),ALLOCATABLE :: FR(:)            !< FR(LCM):  Filtration rate [1 filtered / ind. / hour]  
      REAL(RKD),ALLOCATABLE :: BMG(:)           !< Basal metabolic rate [1/day]
      REAL(RKD),ALLOCATABLE :: B_RESPI(:)       !< Basal metabolic rate [1/day]
      REAL(RKD),ALLOCATABLE :: B_URINE(:) 
      REAL(RKD),ALLOCATABLE :: B_GRAZI(:) 
      REAL(RKD),ALLOCATABLE :: B_FECAL(:) 
      REAL(RKD),ALLOCATABLE :: B_DEATH(:) 
      REAL(RKD),ALLOCATABLE :: B_RPOC(:) 
      REAL(RKD),ALLOCATABLE :: B_LPOC(:) 
      REAL(RKD),ALLOCATABLE :: B_DOC(:) 
      REAL(RKD),ALLOCATABLE :: AQ_RATE(:,:)     !< Harvest rate [-]
      REAL(RKD),ALLOCATABLE :: WQAQ(:)          !< Dry meat weight per grid cell [g C]
      REAL(RKD),ALLOCATABLE :: WQEA(:)          !< Shellfish individuals per grid cell [#]
      REAL(RKD),ALLOCATABLE :: WQO(:)
      REAL(RKD),ALLOCATABLE :: WQVO(:,:)
      REAL(RKD),ALLOCATABLE :: KAALGAE(:)       !< Assimilation efficiency for algae                       
    END TYPE
    
    TYPE SHELLFISHZONE
      INTEGER :: ITYPE                          !< Type of aquaculture farm (not used)    
      REAL(RKD) :: AQ_DEP                       !< Depth of aquaculture farm [m]
      REAL(RKD),ALLOCATABLE :: KDEATH(:)        !< Specific mortality rate [1/day]
    END TYPE
    
    INTEGER :: NSF                              !< NSF: Number of shellfish classes (oyster, scallop, etc)
    INTEGER :: NSFZONES                         !< NSF: Number of shellfish classes (oyster, scallop, etc)
    TYPE(SHELLFISHTYPE),ALLOCATABLE :: SF(:)    !< SF: A shellfish class (oyster, scallop, etc)
    TYPE(SHELLFISHZONE),ALLOCATABLE :: SFZONES(:)
    INTEGER,ALLOCATABLE :: FARMCELL(:)
    INTEGER,ALLOCATABLE :: FARMGR(:)            !< FARMGR(LCM): Type of shellfish farm (Zone ID)
    !REAL(RKD) :: TEMPB                          !< Base temperature [deg C]
    REAL(RKD),ALLOCATABLE :: HARVTIMES(:),HARVDW(:,:)   !< mid-harvest
    LOGICAL,ALLOCATABLE :: HARVDONE(:)
    INTEGER :: ISFFARM,NSFCELLS,NSFHARV,ISHFWQ

  CONTAINS

  SUBROUTINE INIT_SHELLFISH(LCM,NK)
    INTEGER, INTENT(IN) :: LCM,NK
    INTEGER :: I,J,K,L
    
    ALLOCATE(SF(NSF))
    DO I=1,NSF
      ALLOCATE(SF(i).CELL_C(LCM,NK),SF(i).INDI_C(LCM,NK),SF(i).INDI(LCM,NK),SF(i).INDI_D(LCM,NK),SF(i).SHLEN(LCM,NK))
      ALLOCATE(SF(i).SPAWN(LCM,NK),SF(i).CRB(LCM,NK),SF(i).NP(LCM,NK),SF(i).PR(LCM,NK))
      ALLOCATE(SF(i).HARV_C(LCM,NK),SF(i).INDI_H(LCM,NK))
      ALLOCATE(SF(i).FR(LCM))
      ALLOCATE(SF(i).BMG(LCM),SF(i).B_RESPI(LCM),SF(i).B_GRAZI(LCM))      
      ALLOCATE(SF(i).B_FECAL(LCM),SF(i).B_DEATH(LCM),SF(i).B_URINE(LCM)) 
      ALLOCATE(SF(i).B_RPOC(LCM),SF(i).B_LPOC(LCM),SF(i).B_DOC(LCM)) 
      !ALLOCATE(SF(i).AQ_RATE(NSFHARV,LCM))
      ALLOCATE(SF(i).WQAQ(LCM),SF(i).WQEA(LCM))
      ALLOCATE(SF(i).WQO(LCM),SF(i).WQVO(LCM,NK))
      ALLOCATE(SF(i).KAALGAE(NALGAEM))
      
      DO L=1,LCM
        DO K=1,NK
          SF(i).CELL_C(L,K) = 0.
          SF(i).INDI(L,K) = 0.
          SF(i).INDI_C(L,K) = 0.
          SF(i).INDI_D(L,K) = 0.
          SF(i).INDI_H(L,K) = 0.  ! *** DKT
          SF(i).HARV_C(L,K) = 0.  ! *** DKT
          SF(i).NP(L,K) = 0.
          SF(i).PR(L,K) = 0.
          SF(i).CRB(L,K) = 0.
          SF(i).SPAWN(L,K) = 0.
        ENDDO
        SF(i).WQAQ(L) = 0.
        SF(i).WQEA(L) = 0.
        SF(i).FR(L) = 0.
        SF(i).BMG(L) = 0.
        SF(i).B_RESPI(L) = 0.
        SF(i).B_GRAZI(L) = 0.
        SF(i).B_DEATH(L) = 0.
        SF(i).B_FECAL(L) = 0.
        SF(i).B_URINE(L) = 0.
        SF(i).B_RPOC(L) = 0.
        SF(i).B_LPOC(L) = 0.
        SF(i).B_DOC(L) = 0.        
      ENDDO
      SF(i).KAALGAE = 0.
    ENDDO

    
    ALLOCATE(FARMCELL(LCM),FARMGR(LCM))    
    DO L=1,LCM
      FARMCELL(L)=0
      DO I=1,NSF
        SF(i).WQAQ(L)=0.
        SF(i).WQEA(L)=0.
        DO K=1,KC
            SF(i).INDI_C(L,K) = 0.
            SF(i).INDI(L,K) = 0.
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE 
  
  SUBROUTINE FREE_SHELLFISH()
    INTEGER :: I

     DO I=1,NSF
      DEALLOCATE(SF(i).CELL_C,SF(i).INDI_C,SF(i).INDI,SF(i).INDI_D,SF(i).SHLEN)
      DEALLOCATE(SF(i).SPAWN,SF(i).CRB,SF(i).NP,SF(i).PR)
      DEALLOCATE(SF(i).HARV_C,SF(i).INDI_H)
      DEALLOCATE(SF(i).FR)
      DEALLOCATE(SF(i).BMG,SF(i).B_RESPI,SF(i).B_GRAZI)      
      DEALLOCATE(SF(i).B_FECAL,SF(i).B_DEATH,SF(i).B_URINE) 
      DEALLOCATE(SF(i).B_RPOC,SF(i).B_LPOC,SF(i).B_DOC) 
      DEALLOCATE(SF(i).AQ_RATE)
      DEALLOCATE(SF(i).WQAQ,SF(i).WQEA)
      DEALLOCATE(SF(i).WQO,SF(i).WQVO)
      DEALLOCATE(SF(i).KAALGAE)
    ENDDO
    DO I=1,NSFZONES
      DEALLOCATE(SFZONES(i).KDEATH)
    ENDDO
    DEALLOCATE(SF)
    DEALLOCATE(SFZONES)
    DEALLOCATE(HARVTIMES,HARVDW,HARVDONE)
    DEALLOCATE(FARMCELL,FARMGR)
  END SUBROUTINE 

  SUBROUTINE SHELLFISH_REDIST(L)
    REAL(RKD) :: AQDEP,ADEP,BDEP
    INTEGER :: I,K,L,IDX

    IF(FARMCELL(L) > 0 )THEN
      AQDEP = HP(L)
      IDX = FARMGR(L)          
      IF(IDX > 0 .AND. IDX <= NSFZONES )THEN
        AQDEP = SFZONES(IDX).AQ_DEP
      ENDIF
      DO I = 1,NSF
        IF(AQDEP .LE. HP(L) )THEN
          ADEP = AQDEP
          BDEP = AQDEP
        ELSE
          ADEP = HP(L)
          BDEP = HP(L)
        ENDIF
        DO K=KC,1,-1
          IF( ADEP > 0 .AND. ADEP >= HPK(L,K) )THEN
            IF( SF(i).WQEA(L) > 0.0 )THEN
              SF(i).INDI(L,K) = SF(i).WQEA(L)*(HPK(L,K)/BDEP)       !!! WQEA(L) : vertical sum of ea of pre-step
              SF(i).CELL_C(L,K) = SF(i).WQAQ(L)*(HPK(L,K)/BDEP)     !!! WQAQ(L) : vertical sum of WQV(i) of pre-step
              SF(i).INDI_C(L,K) = SF(i).CELL_C(L,K)/SF(i).INDI(L,K)
            ENDIF
            ADEP=ADEP-HPK(L,K)
          ELSE
            IF( SF(i).WQEA(L) > 0.0 )THEN
              SF(i).CELL_C(L,K) = (ADEP/BDEP)*SF(i).WQAQ(L)
              SF(i).INDI(L,K) = (ADEP/BDEP)*SF(i).WQEA(L)
              SF(i).INDI_C(L,K) = SF(i).CELL_C(L,K)/SF(i).INDI(L,K)
            ENDIF
            ADEP=0
            EXIT
          ENDIF
        ENDDO
      ENDDO      
    ENDIF
  END SUBROUTINE 

  SUBROUTINE SHELLFISH_LAYERSUM(L)
    REAL(RKD) :: TYAQ
    INTEGER :: I,L,K

    DO I = 1,NSF
      TYAQ = 0.0
      SF(i).WQAQ(L)=0.
      SF(i).WQEA(L)=0
      DO K=KC,1,-1
        SF(i).WQAQ(L) = SF(i).WQAQ(L) + SF(i).CELL_C(L,K)
        SF(i).WQEA(L) = SF(i).WQEA(L) + SF(i).INDI(L,K)
      ENDDO
!      WQV(L,KC,18)=WQAQ1(L)
!      WQV(L,KC,1)=WQAQ2(L)
!      WQV(L,KC,3)=WQEA1(L)
      TYAQ = TYAQ + SF(i).WQAQ(L)
    ENDDO
  END SUBROUTINE 

  SUBROUTINE SHELLFISH_LENGTH(L)
    INTEGER :: I,K,L
    IF(FARMCELL(L) > 0 )THEN
      DO I = 1,NSF
        DO K=1,KC
          SF(i).SHLEN(L,K) = SF(i).CWLEN*(SF(i).INDI_C(L,K)/SF(i).WQACF)**SF(i).PWLEN
        ENDDO
      ENDDO     
    ENDIF
  END SUBROUTINE
  
  SUBROUTINE SHELLFISH_GROWTH(J,SF,L,K)
    USE GLOBAL
    INTEGER, INTENT(IN) :: J
    TYPE(SHELLFISHTYPE), INTENT(INOUT) :: SF
    REAL(RKD) :: DW,KDEATH,INGEST,SFCHC,SFCHD,SFCHG,SFROC,SFLOC,FAC,WVOL,TSS,DOX
    REAL(RKD) :: SHFLEN,REFF,FRATIO,EGGS,LARVAE,RW,RT,GRAZIND,ASSIMIND,RESPIND,FDEATH
    REAL(RKD) :: SFALGAE(NALGAEM)
    INTEGER :: L,K,IDX,NAL

    IF( SF.CELL_C(L,K) > 0.0 )THEN
      
      DW = SF.INDI_C(L,K)/SF.WQACF          ! Biomass (gram C dry weight)
      
      SHFLEN = SF.CWLEN*DW**SF.PWLEN        ! Shellfish length (Kobayashi et al., 1997)
      IF( SF.SHLEN(L,K) < SHFLEN) SF.SHLEN(L,K) = SHFLEN
      SHFLEN = SF.SHLEN(L,K)
      
      ! Max. filtration rate (liter filtered per ind per hour) 
      IF( SF.IFILT == 0 )THEN
        ! Cerco and Noel (2007)
        SF.FR(L) = SF.FRMAX*(1000./24.)*DW 
      ELSEIF (SF.IFILT == 2 )THEN
        ! Doering and Oviatt (1986)
        SF.FR(L) = (60./1000.)/2.95*(SHFLEN**0.96)*(TWQ(L)**0.95) 
      ELSEIF (SF.IFILT == 3 )THEN
        ! Cloern (1982)
        SF.FR(L) = 7.0*(DW**0.67)
      ELSEIF (SF.IFILT == 4 )THEN
        ! Coughlan & Ansell (1964)
        SF.FR(L) = 2.59*(DW**0.73)
      ELSE ! Default SF.IFILT = 1
        ! Kobayashi et al. (1997)
        IF( DW < 2.0 )THEN
          SF.FR(L) = 0.117*DW**3 - 1.05*DW**2 + 3.09*DW + 0.133
        ELSE
          SF.FR(L) = 2.51*DW**0.279 
        ENDIF
      ENDIF
      
      ! Temperature effect
      IF( SF.IFTEMP == 1 )THEN
        ! Kobayashi et al. (1997)
        IF( TWQ(L) < 7.0 )THEN
          SF.FR(L) = SF.FR(L) * 0.59
        ELSE
          SF.FR(L) = SF.FR(L) * (TWQ(L)**0.5)/4.47
        ENDIF
      ELSEIF (SF.IFTEMP == 2 )THEN
        ! Cerco and Noel (2007)
          SF.FR(L) = SF.FR(L) * EXP(-SF.KTG*(TWQ(L)-SF.TOPT)**2)
      ENDIF
      
      ! Salinity effect
      IF( SF.IFSALT == 1 )THEN
        ! Buzzelli et al. (2015)
        SF.FR(L) = SF.FR(L) * ((-0.0017*SWQ(L)+0.0084)*SWQ(L)-0.1002)
      ELSEIF (SF.IFSALT == 2 )THEN
        ! Fulford (2007)
        SF.FR(L) = SF.FR(L) * (0.0926*SWQ(L)-0.139)
      ELSEIF (SF.IFSALT == 3 )THEN
        ! Cerco and Noel (2005)
        SF.FR(L) = SF.FR(L) * 0.5*(1+TANH(SWQ(L)-SF.KHSOY))
      ELSEIF (SF.IFSALT == 4 )THEN
        ! Quayle (1988) and Mann et al. (1991)
        IF( SWQ(L) <= 10.0 )THEN
    	    SF.FR(L) = 0.0
        ELSEIF (SWQ(L) < 20.0 )THEN
    	    SF.FR(L) = SF.FR(L) * (SWQ(L)-10.0)/10.0
        ENDIF
      ELSEIF (SF.IFSALT == 5 )THEN
        ! Loosanoff (1958)
        IF( SWQ(L) <= 3.5 )THEN
    	    SF.FR(L) = 0.0
        ELSEIF (SWQ(L) < 7.5 )THEN
    	    SF.FR(L) = SF.FR(L) * (SWQ(L)-3.5)/4.0
        ENDIF
      ENDIF

      ! TSS effect
      IF( SF.IFTSS > 1 )THEN
        TSS = 0.
        IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
          TSS = SEDT(L,K)+SNDT(L,K) 
        ENDIF
        IF( TSS > 0 )THEN
          IF( SF.IFTSS == 1 )THEN
            ! Fulford (2007)
            IF( TSS > 25 )THEN
              SF.FR(L) = SF.FR(L) * (10.364*(LOG(TSS))**(-2.0477))   
            ENDIF
          ELSEIF (SF.IFTSS == 2 )THEN
            ! Cerco and Noel (2005)
            IF( TSS > 100.0 )THEN
    	        SF.FR(L) = 0.0
            ELSEIF (TSS > 25. )THEN
    	        SF.FR(L) = SF.FR(L) * 0.2
            ELSEIF (TSS < 5. )THEN
    	        SF.FR(L) = SF.FR(L) * 0.1
            ENDIF
          ELSEIF (SF.IFTSS == 3 )THEN
            ! Hofmann et al. (1992)
            TSS = TSS/1000. ! Convert mg/L to g/L
            SF.FR(L) = SF.FR(L) * (1. - 0.01*(LOG10(TSS) + 3.38)/0.418)   
          ENDIF
        ENDIF
      ENDIF
      
      ! DO effect
      IF( SF.IFDO == 1 .AND. ISKINETICS(IDOX) > 0 )THEN
        DOX = WQV(L,K,IDOX)
        ! Cerco and Noel (2005)
        SF.FR(L) = SF.FR(L) / (1.  + EXP(1.1*(SF.DOHX-DOX)/(SF.DOHX-SF.DOQX)))   
      ENDIF
      
      INGEST = SF.FR(L) * 24.0/1000.0 * DTWQ          ! Ingestion [L]
      IF( ISWQLVL == 0 )THEN
        SFCHC =  INGEST * SF.KACHC*WQV(L,K,1) 
        SFCHD =  INGEST * SF.KACHD*WQV(L,K,2) 
        SFCHG =  INGEST * SF.KACHG*WQV(L,K,3) 
      ELSE
        DO NAL = 1,NALGAE
          SFALGAE(NAL) = INGEST * SF.KAALGAE(NAL)*WQV(L,K,19+NAL)
        ENDDO
      ENDIF
      
      SFROC =  INGEST * SF.KAPOC*WQV(L,K,IROC) 
      SFLOC =  INGEST * SF.KAPOC*WQV(L,K,ILOC) 
      IF(ISWQLVL == 0 )THEN
        GRAZIND = SFCHC + SFCHD + SFCHG + SFROC + SFLOC ! Grazing [g C]
      ELSE
        GRAZIND = SFROC + SFLOC
        DO NAL = 1,NALGAE
          GRAZIND = GRAZIND +  SFALGAE(NAL)
        ENDDO
      ENDIF
      ASSIMIND = SF.K_ASSIM * GRAZIND                 ! Assimilation [g C]
      
      ! Respiration
      IF( SF.IRES == 1 )THEN
        IF( TWQ(L) > SF.TEMP_RESPI )THEN
          SF.BMG(L) = SF.K_RESPI * SF.THETA_RESPI**(TWQ(L)-SF.TEMP_RESPI)
        ELSE
          SF.BMG(L) = 0.
        ENDIF
        RESPIND = SF.INDI_C(L,K) * SF.BMG(L)*DTWQ
      ELSEIF (SF.IRES == 2 )THEN
        ! Cerco and Noel (2005)
        SF.BMG(L) = SF.K_RESPI * EXP(SF.THETA_RESPI*(TWQ(L)-SF.TEMP_RESPI))
        RESPIND = SF.INDI_C(L,K) * SF.BMG(L)*DTWQ
      ELSE  ! Default SF.IRES == 0
        RW = (12.6*TWQ(L) + 69.7)*DW**-0.25       ! Dame (1972)
        !RW = (31.0*TWQ(L) - 22.0)*DW**-0.3
        RT = 1.0
        IF( SF.IRTEMP == 1 )THEN
          ! Shumway and Koehn (1992)
          IF( TWQ(L) < 20.0 )THEN
              RT = (0.007*TWQ(L) + 2.099)
          ELSE
              RT = (0.0915*TWQ(L) + 1.324)
          ENDIF
        ENDIF
        IF( SF.IRSALT == 1 )THEN
          ! Shumway and Koehn (1982)
          IF( SWQ(L) <= 15.0 )THEN
            RW = RW*RT
          ELSEIF (SWQ(L) < 20.0 )THEN
    	      RW = RW*(1.0+(((RT-1.0)/5.0)*(20.0-SWQ(L))))
          ENDIF
        ELSEIF (SF.IRSALT == 2 )THEN
          ! Hoffman et al. (1992)
          IF( SWQ(L) <= 10.0 )THEN
            RW = RW*RT
          ELSEIF (SWQ(L) < 15.0 )THEN
    	      RW = RW*(1.0+(((RT-1.0)/5.0)*(15.0-SWQ(L))))
          ENDIF
        ENDIF
        RESPIND = RW * (SF.RESCAL/SF.SHF_CAL) * 24.0/1000.0 * DTWQ
      ENDIF

      SF.NP(L,K) = ASSIMIND - RESPIND           ! Net production

      SF.PR(L,K) = 0.0
      SF.SPAWN(L,K) = 0.0
      IF( DW > SF.SPAWN_WEIGHT )THEN            ! Adult weight
        IF(  SF.NP(L,K) > 0.0 )THEN
          ! Kusaka et al. (1991)
          IF( TWQ(L) >= 27.0  )THEN
            REFF = 0.8
          ELSEIF (TWQ(L) > 23.0 )THEN
            REFF = 0.2*TWQ(L) - 4.6
          ELSE
            REFF = 0.0
          ENDIF
          SF.PR(L,K) = REFF * SF.NP(L,K)            ! Reproductive tissue production
          SF.CRB(L,K) = SF.CRB(L,K) + SF.PR(L,K)    ! Cumulative reproductive biomass
        ELSEIF( SF.NP(L,K) < 0.0 )THEN
          IF( SF.CRB(L,K) > ABS(SF.NP(L,K)) )THEN
            SF.CRB(L,K) = SF.CRB(L,K) + SF.NP(L,K)
            SF.NP(L,K) = 0.0
          ELSE
            SF.NP(L,K) = SF.NP(L,K) + SF.CRB(L,K)
            SF.CRB(L,K) = 0.0
          ENDIF
        ENDIF
               
        ! Reproduction
        IF( SF.ISPAWN > 0 )THEN
          IF( (TWQ(L) >= SF.SPAWN_TEMP) .AND. (SF.CRB(L,K) >= SF.SPAWN_RATIO*DW) )THEN
        	FRATIO = 0.021*10*SHFLEN - 0.62          ! Kennedy (1982)
            SF.SPAWN(L,K) = FRATIO/(FRATIO+1) * SF.CRB(L,K)
            EGGS = SF.SPAWN(L,K)*SF.SHF_CAL/(SF.EGG_CAL*SF.EGG_WEIGHT)
            IF( TWQ(L) >= SF.EGG_TEMP1 .AND. TWQ(L) <= SF.EGG_TEMP2 .AND. & 
                SWQ(L) >= SF.EGG_SALT1 .AND. SWQ(L) <= SF.EGG_SALT2 )THEN
              !LARVAE(L,K) = SF.EGG_HATCHED*EGGS  ! Number of shellfish lavae
              !LARVAE_LIFE(L,K) = 0.0             ! A larvae life starts
            ELSE
              !LARVAE(L,K) = 0.0
            ENDIF
            DW = DW - SF.SPAWN(L,K)
            SF.CRB(L,K) = 0.0
          ENDIF
        ENDIF
      ENDIF      
      
      DW = DW + SF.NP(L,K)
      SF.INDI_C(L,K) = DW * SF.WQACF
            
      IDX = FARMGR(L)          
      IF(IDX > 0 .AND. IDX <= NSFZONES )THEN
        KDEATH = SFZONES(IDX).KDEATH(J)
      ELSE
        KDEATH = SF.K_DEATH
      ENDIF
      SF.B_DEATH(L) = SF.CELL_C(L,K)*KDEATH*DTWQ                          
      SF.INDI_D(L,K) = SF.INDI(L,K)*KDEATH*DTWQ
      IF( SF.ITDEATH > 0 )THEN
        FDEATH = SF.THETA_DEATH**(TWQ(L)-SF.TEMP_DEATH)
        SF.B_DEATH(L) = SF.B_DEATH(L) * FDEATH
        SF.INDI_D(L,K) = SF.INDI_D(L,K) * FDEATH
      ENDIF

      ! Cell - layer values 
      SF.B_GRAZI(L) = SF.INDI(L,K) * GRAZIND                                               ! [g C]
      SF.B_RESPI(L) = SF.INDI(L,K) * RESPIND                                               ! [g C]
      SF.B_FECAL(L) = SF.B_GRAZI(L)*(1 - SF.K_ASSIM)                                       ! [g C]                
      SF.B_URINE(L) = (SF.B_GRAZI(L)-SF.B_FECAL(L))*(SF.K_ASSIM-SF.K_GROEF)                ! [g C]
                
      SF.B_RPOC(L) = SF.B_FECAL(L) * SF.FFROC + SF.B_DEATH(L) * SF.DFROC
      SF.B_LPOC(L) = SF.B_FECAL(L) * SF.FFLOC + SF.B_DEATH(L) * SF.DFLOC
      SF.B_DOC(L)  = SF.B_FECAL(L) * SF.FFDOC + SF.B_DEATH(L) * SF.DFDOC
      
      SF.CELL_C(L,K) = SF.CELL_C(L,K) + SF.B_GRAZI(L) - SF.B_FECAL(L) - SF.B_RESPI(L) - SF.B_DEATH(L) - SF.B_URINE(L)
      SF.CELL_C(L,K) = SF.CELL_C(L,K) - SF.INDI(L,K) * SF.SPAWN(L,K)
      SF.INDI(L,K) = SF.INDI(L,K) - SF.INDI_D(L,K)
     
      !WQV(L,K,20) = SF.INDI_C(L,K)
      !WQO(L,20) = WQVO(L,K,20) + WQV(L,K,20)
      SF.WQO(L) = SF.WQVO(L,K) + SF.INDI_C(L,K)
      IF( ISHFWQ > 0 .AND. HPK(L,K) > 0 )THEN
        WVOL = DXYP(L) * HPK(L,K)                 ! water volume [m3]
        FAC = 1.0  / WVOL                         ! Converts [g C] to [mg/L]
        IF(ISWQLVL == 0 )THEN
          IF( ISKINETICS(1) > 0) WQV(L,K,1) = WQV(L,K,1) - FAC*SFCHC*SF.INDI(L,K)
          IF( ISKINETICS(2) > 0) WQV(L,K,2) = WQV(L,K,2) - FAC*SFCHD*SF.INDI(L,K)
          IF( ISKINETICS(3) > 0) WQV(L,K,3) = WQV(L,K,3) - FAC*SFCHG*SF.INDI(L,K)
        ELSE
          DO NAL = 1,NALGAE
            WQV(L,K,19+NAL) = WQV(L,K,19+NAL) - FAC*SFALGAE(NAL)*SF.INDI(L,K)
          ENDDO
        ENDIF
        
        IF( ISKINETICS(IROC) > 0) WQV(L,K,IROC) = WQV(L,K,IROC) + FAC*(SF.B_RPOC(L) - SFROC*SF.INDI(L,K))
        IF( ISKINETICS(ILOC) > 0) WQV(L,K,ILOC) = WQV(L,K,ILOC) + FAC*(SF.B_LPOC(L) - SFLOC*SF.INDI(L,K))
        IF( ISKINETICS(IDOC) > 0) WQV(L,K,IDOC) = WQV(L,K,IDOC) + FAC*(SF.B_DOC(L))
      
        IF( ISKINETICS(IROP) > 0) WQV(L,K,IROP) = WQV(L,K,IROP) + FAC*(SF.UFROP * SF.B_URINE(L))
        IF( ISKINETICS(ILOP) > 0) WQV(L,K,ILOP) = WQV(L,K,ILOP) + FAC*(SF.UFLOP * SF.B_URINE(L))
        IF( ISKINETICS(IDOP) > 0) WQV(L,K,IDOP) = WQV(L,K,IDOP) + FAC*(SF.UFDOP * SF.B_URINE(L))
        IF( ISKINETICS(IP4D) > 0) WQV(L,K,IP4D) = WQV(L,K,IP4D) + FAC*(SF.UFP4D * SF.B_URINE(L))
        IF( ISKINETICS(IRON) > 0) WQV(L,K,IRON) = WQV(L,K,IRON) + FAC*(SF.UFRON * SF.B_URINE(L))
        IF( ISKINETICS(ILON) > 0) WQV(L,K,ILON) = WQV(L,K,ILON) + FAC*(SF.UFLON * SF.B_URINE(L))
        IF( ISKINETICS(IDON) > 0) WQV(L,K,IDON) = WQV(L,K,IDON) + FAC*(SF.UFDON * SF.B_URINE(L))
        IF( ISKINETICS(INHX) > 0) WQV(L,K,INHX) = WQV(L,K,INHX) + FAC*(SF.UFNHX * SF.B_URINE(L))
      
        IF( ISKINETICS(ISUU) > 0) WQV(L,K,ISUU) = WQV(L,K,ISUU) + FAC*(SF.FFSUU * SF.B_FECAL(L))
        IF( ISKINETICS(IDOX) > 0) WQV(L,K,IDOX) = WQV(L,K,IDOX) - FAC*RESPIND*SF.INDI(L,K)
      ENDIF
    ENDIF
  END SUBROUTINE
  
  SUBROUTINE SHELLFISH_HARVEST(LF,LL)
    INTEGER, INTENT(IN) :: LF,LL
    INTEGER :: I,J,K,L
    REAL :: TIME
    
    TIME = TIMESEC/TCON  
    
    DO J=1,NSFHARV
      IF((TIME >= HARVTIMES(J)) .AND. (.NOT. HARVDONE(J)) )THEN
        DO L=LF,LL
          DO I = 1,NSF
            SF(i).WQAQ(L)=0.
            SF(i).WQEA(L)=0
            DO K=KC,1,-1
              IF( SF(i).CELL_C(L,K) >= 0.0  )THEN
                SF(i).INDI_C(L,K) = SF(i).CELL_C(L,K)/(SF(i).INDI(L,K) + 1.E-18)         ! Individual weight
                IF( SF(i).INDI_C(L,K) >= HARVDW(I,J) )THEN                    ! Is greater than or equal to harvest size?
                  SF(i).INDI_H(L,K) = SF(i).AQ_RATE(J,L)*SF(i).INDI(L,K)      ! Number of harvested individuals
                  SF(i).HARV_C(L,K) = SF(i).INDI_H(L,K)*SF(i).INDI_C(L,K)     ! Harvested biomass
                  SF(i).INDI(L,K) = SF(i).INDI(L,K) - SF(i).INDI_H(L,K)       ! Remaining individuals
                  SF(i).INDI_C(L,K) = SF(i).INDI_C(L,K) - SF(i).HARV_C(L,K)   ! Remaining biomass
                ENDIF
                SF(i).WQAQ(L) = SF(i).WQAQ(L) + SF(i).CELL_C(L,K)
                SF(i).WQEA(L) = SF(i).WQEA(L) + SF(i).INDI(L,K)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        HARVDONE(J) = .TRUE.
      ENDIF
    ENDDO
  END SUBROUTINE
  
  SUBROUTINE READ_SHELLFISH_JSON(LCM,NK)
    USE fson
    USE mod_fson_value, ONLY: fson_value_count, fson_value_get
      
    TYPE(fson_value), POINTER :: json_data, items, item
    CHARACTER*80 :: STR*200
    INTEGER :: I,J,K,L,LCM,NK,IX,JX,IDX,NAL
    REAL :: ADEP,AQDEP
    REAL,ALLOCATABLE :: ARR1D(:),ARR2D(:,:),SFWEIGHT(:)

    json_data => fson_parse("shellfish.jnp")   
    CALL fson_get(json_data, "number_of_species", NSF)
    CALL fson_get(json_data, "wq_interact_flag", ISHFWQ)    
    
    CALL INIT_SHELLFISH(LCM,NK)
    
    CALL fson_get(json_data, "species", items)
    DO i = 1, fson_value_count(items) !NSF
      item => fson_value_get(items, i)
      !CALL fson_get(item, "name", SF(i).NAME)
      !CALL fson_get(item, "index", SF(i).IDX)
        
      ! *** General parameters
      CALL fson_get(item, "stoichiometry.nitrogen_to_carbon_ratio", SF(i).WQANC)
      CALL fson_get(item, "dry_meat_carbon", SF(i).WQACF)
      CALL fson_get(item, "dry_meat_caloric", SF(i).SHF_CAL)        
      CALL fson_get(item, "length_weight_factor", SF(i).CWLEN)
      CALL fson_get(item, "length_weight_exponent", SF(i).PWLEN)        

      CALL fson_get(item, "filtration.max_filtration_opt", SF(i).IFILT)
      CALL fson_get(item, "filtration.temperature_opt", SF(i).IFTEMP)
      CALL fson_get(item, "filtration.salinity_opt", SF(i).IFSALT)
      CALL fson_get(item, "filtration.sediment_opt", SF(i).IFTSS)
      CALL fson_get(item, "filtration.dissolved_oxygen_opt", SF(i).IFDO)
      CALL fson_get(item, "filtration.max_filtration_rate", SF(i).FRMAX)
      CALL fson_get(item, "filtration.temperature_optimal", SF(i).TOPT)
      CALL fson_get(item, "filtration.temperature_effect", SF(i).KTG)
      CALL fson_get(item, "filtration.halved_rate_salinity", SF(i).KHSOY)
      CALL fson_get(item, "filtration.one_half_rate_DO", SF(i).DOHX)
      CALL fson_get(item, "filtration.one_fourth_rate_DO", SF(i).DOQX)

      CALL fson_get(item, "respiration.respiration_opt", SF(i).IRES)
      CALL fson_get(item, "respiration.temperature_opt", SF(i).IRTEMP)
      CALL fson_get(item, "respiration.salinity_opt", SF(i).IRSALT)
      CALL fson_get(item, "respiration.sediment_opt", SF(i).IRTSS)
      CALL fson_get(item, "respiration.dissolved_oxygen_opt", SF(i).IRDO)
      CALL fson_get(item, "respiration.basal_metabolic_rate", SF(i).K_RESPI)
      CALL fson_get(item, "respiration.ref_temperature", SF(i).TEMP_RESPI)
      CALL fson_get(item, "respiration.temperature_coeff", SF(i).THETA_RESPI)
      CALL fson_get(item, "respiration.energy_consumption", SF(i).RESCAL)
      
      ALLOCATE(ARR1D(3))
      CALL fson_get(item, "growth.assimilation_of_algaes", ARR1D)
      CALL fson_get(item, "growth.assimilation_of_POC", SF(i).KAPOC)
      CALL fson_get(item, "growth.assimilation_efficiency", SF(i).K_ASSIM)
      CALL fson_get(item, "growth.growth_efficiency", SF(i).K_GROEF)
      IF(ISWQLVL == 0 )THEN
        SF(i).KACHC = ARR1D(1)
        SF(i).KACHD = ARR1D(2)
        SF(i).KACHG = ARR1D(3)
      ELSE
        DO NAL = 1,NALGAE
          SF(i).KAALGAE(NAL) = ARR1D(NAL)
        ENDDO
      ENDIF
      DEALLOCATE(ARR1D)
      
      CALL fson_get(item, "excretion.fraction_of_carbon.RPOC", SF(i).FFROC)
      CALL fson_get(item, "excretion.fraction_of_carbon.LPOC", SF(i).FFLOC)
      CALL fson_get(item, "excretion.fraction_of_carbon.DOC", SF(i).FFDOC)
      CALL fson_get(item, "excretion.fraction_of_phosphorus.RPOP", SF(i).UFROP)
      CALL fson_get(item, "excretion.fraction_of_phosphorus.LPOP", SF(i).UFLOP)
      CALL fson_get(item, "excretion.fraction_of_phosphorus.DOP", SF(i).UFDOP)
      CALL fson_get(item, "excretion.fraction_of_phosphorus.PO4", SF(i).UFP4D)
      CALL fson_get(item, "excretion.fraction_of_nitrogen.RPON", SF(i).UFRON)
      CALL fson_get(item, "excretion.fraction_of_nitrogen.LPON", SF(i).UFLON)
      CALL fson_get(item, "excretion.fraction_of_nitrogen.DON", SF(i).UFDON)
      CALL fson_get(item, "excretion.fraction_of_nitrogen.NH4", SF(i).UFNHX)
      CALL fson_get(item, "excretion.fraction_of_silica.SU", SF(i).FFSUU)
      
      CALL fson_get(item, "mortality.specific_rate", SF(i).K_DEATH)
      CALL fson_get(item, "mortality.temperature_dependent_flag", SF(i).ITDEATH)
      CALL fson_get(item, "mortality.ref_temperature", SF(i).TEMP_DEATH)      
      CALL fson_get(item, "mortality.temperature_coeff", SF(i).THETA_DEATH)      
      CALL fson_get(item, "mortality.dead_shellfish.RPOC", SF(i).DFROC)
      CALL fson_get(item, "mortality.dead_shellfish.LPOC", SF(i).DFLOC)
      CALL fson_get(item, "mortality.dead_shellfish.DOC", SF(i).DFDOC)
     
      CALL fson_get(item, "spawning.spawning_flag", SF(i).ISPAWN)
      CALL fson_get(item, "spawning.min_dry_weight", SF(i).SPAWN_WEIGHT)
      CALL fson_get(item, "spawning.min_temperature", SF(i).SPAWN_TEMP)
      CALL fson_get(item, "spawning.min_cumulative_reproductive_biomass", SF(i).SPAWN_RATIO)
     
      CALL fson_get(item, "eggs.dry_weight", SF(i).EGG_WEIGHT)
      CALL fson_get(item, "eggs.caloric_content", SF(i).EGG_CAL)
      CALL fson_get(item, "eggs.min_temp_to_survive", SF(i).EGG_TEMP1)
      CALL fson_get(item, "eggs.max_temp_to_survive", SF(i).EGG_TEMP2)
      CALL fson_get(item, "eggs.min_salt_to_survive", SF(i).EGG_SALT1)
      CALL fson_get(item, "eggs.max_salt_to_survive", SF(i).EGG_SALT2)

      CALL fson_get(item, "larvae.fraction_hatched", SF(i).EGG_HATCHED)
      CALL fson_get(item, "larvae.fraction_to_spat_survival", SF(i).LARVAE_SURVIVAL)
      CALL fson_get(item, "larvae.fraction_recruited_per_spawn", SF(i).LARVAE_RECRUITED)
      CALL fson_get(item, "larvae.life_span", SF(i).LARVAE_SPAN)

      ! *** Initial condition
      CALL fson_get(item, "initial_conditions",  SF(i).IC)    
    ENDDO
    
    json_data => fson_parse("shffarm.jnp")   
    CALL fson_get(json_data, "number_of_zones", NSFZONES)
    CALL fson_get(json_data, "number_of_cells", NSFCELLS)    
    CALL fson_get(json_data, "number_of_harvests", NSFHARV)   
    
    ALLOCATE(SFZONES(NSFZONES))
    DO J=1,NSFZONES
      ALLOCATE(SFZONES(J).KDEATH(NSF))
    ENDDO
    ALLOCATE(HARVTIMES(NSFHARV),HARVDONE(NSFHARV))
    ALLOCATE(HARVDW(NSF,NSFHARV))
    DO J=1,NSFHARV
      HARVDONE(J) = .FALSE.
    ENDDO
    DO i = 1,NSF
      ALLOCATE(SF(i).AQ_RATE(NSFHARV,LCM))
    ENDDO
    ALLOCATE(ARR1D(NSF), SFWEIGHT(NSF))    
    
    CALL fson_get(json_data, "shellfish_harvest.mid_days", HARVTIMES)   
    CALL fson_get(json_data, "shellfish_harvest.min_sizes", ARR2D)
    IF( .NOT. ALLOCATED(ARR2D) )THEN
      ALLOCATE(ARR2D(NSF,NSFHARV))
      ARR2D = 0.
    ENDIF
    
    DO i = 1, NSF
        DO j=1, NSFHARV
            HARVDW(i, j) = ARR2D(i, j)
        ENDDO
    ENDDO

    CALL fson_get(json_data, "zones", items)
    DO i = 1, fson_value_count(items) !NSFZONES
      item => fson_value_get(items, i)
      CALL fson_get(item, "Type", SFZONES(i).ITYPE)
      CALL fson_get(item, "installed_depth", SFZONES(i).AQ_DEP)
      CALL fson_get(item, "death_rate", SFZONES(i).KDEATH)
      !DO j=1,NSF
      !  SFZONES(i).KDEATH(j) = VALUES(j)
      !ENDDO
    ENDDO
    
    CALL fson_get(json_data, "cells", items)
    DO k = 1, fson_value_count(items) !NSFCELLS
      item => fson_value_get(items, k)
      CALL fson_get(item, "address.L", L)
      !CALL fson_get(item, "address.I", IX)
      !CALL fson_get(item, "address.J", JX)
      CALL fson_get(item, "zone", FARMGR(L))
      CALL fson_get(item, "initial_conditions.individuals", ARR1D)
      CALL fson_get(item, "initial_conditions.weights", SFWEIGHT)
      CALL fson_get(item, "harvest_rates", ARR2D)
      IF( .NOT. ALLOCATED(ARR2D) )THEN
        ALLOCATE(ARR2D(NSFHARV,LCM))
        ARR2D = 0.
      ENDIF
      FARMCELL(L)=1
      
      ADEP=HMP(L)
      IDX = FARMGR(L)          
      IF(IDX > 0 .AND. IDX <= NSFZONES )THEN
        AQDEP = SFZONES(IDX).AQ_DEP
      ELSE
        WRITE(*,*) 'INVALID SHELLFISH ZONE NUMBER FOR CELL',L
      ENDIF
       DO i = 1,NSF
        DO j = 1,NSFHARV
          SF(i).AQ_RATE(J,L) = ARR2D(i, j)
        ENDDO
        
        IF(AQDEP <= HMP(L) )THEN
          ADEP=AQDEP
        ELSE
          ADEP=HMP(L)
        ENDIF
        
        SF(i).IC = SFWEIGHT(i)
        SF(i).WQEA(L) = ARR1D(i)  !*DXYP(L)*ADEP*SF(i).DRATIO
        IF( ICONTINUE /= 1 )THEN      
          SF(i).WQAQ(L) = SF(i).IC*SF(i).WQEA(L)
        ELSE
          !WQAQ1(L)=WQV(L,3,18)
          !WQAQ2(L)=WQV(L,3,1)
          !WQEA1(L)=WQV(L,3,3)
          !WQEA2(L)=WQEA2(L)
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(ARR1D,SFWEIGHT)
  END SUBROUTINE 
  
END MODULE
