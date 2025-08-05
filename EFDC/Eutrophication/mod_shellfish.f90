! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE SHELLFISHMOD

  use GLOBAL    
  use Variables_WQ

  implicit none

    type SHELLFISHTYPE
      integer :: IFILT                          !< Option for filtration model    
      integer :: IFTEMP                         !< Option for temperature effect on filtration
      integer :: IFSALT                         !< Option for salinity effect on filtration
      integer :: IFTSS                          !< Option for sediment effect on filtration
      integer :: IFDO                           !< Option for dissolved oxygen effect on filtration
      integer :: IRES                           !< Option for respiration model
      integer :: IRTEMP                         !< Option for temperature effect on respiration
      integer :: IRSALT                         !< Option for salinity effect on respiration
      integer :: IRTSS                          !< Option for sediment effect on respiration
      integer :: IRDO                           !< Option for dissolved oxygen effect on respiration
      integer :: ITDEATH                        !< Temperature dependent mortality flag
      integer :: ISPAWN                         !< Option for spawning
      real(RKD) :: IC                           !< Initial biomass [g C]
      real(RKD) :: DRATIO                       !< Density ratio [-]
      real(RKD) :: CWLEN                        !< Factor for the length - weight relationship [-]
      real(RKD) :: PWLEN                        !< Exponent for the length - weight relationship [-]
      real(RKD) :: FRMAX                        !< Maximum filtration rate [-]
      real(RKD) :: TOPT                         !< Temperature for optimal filtration [deg. C]
      real(RKD) :: KTG                          !< Effect of temperature on filtration [-]
      real(RKD) :: KHSOY                        !< Salinity at which filtration rate is halved [ppt]
      real(RKD) :: DOHX                         !< DO at which filtration rate is one-half [g/m^3]
      real(RKD) :: DOQX                         !< DO at which filtration rate is one-fourth [g/m^3]
      real(RKD) :: K_ASSIM                      !< Assimilation efficiency [-]
      real(RKD) :: K_RESPI                      !< Basal metabolic rate [1/day]
      real(RKD) :: K_GROEF                      !< Growth efficiency [-]
      real(RKD) :: K_DEATH                      !< Specific mortality rate [1/day]
      real(RKD) :: TEMP_RESPI                   !< Reference temperature for respiration [deg. C]
      real(RKD) :: THETA_RESPI                  !< Temperature coefficient for respiration [-]
      real(RKD) :: TEMP_DEATH                   !< Reference temperature for mortality [deg. C]
      real(RKD) :: THETA_DEATH                  !< Temperature coefficient for shellfish mortality [-]
      real(RKD) :: DFROC                   !< Fraction of RPOC from dead oyster [-]
      real(RKD) :: DFLOC                   !< Fraction of LPOC from dead oyster [-]
      real(RKD) :: DFDOC                   !< Fraction of DOC from dead oyster [-]
      real(RKD) :: FFROC                   !< Fraction of RPOC in feces [-]
      real(RKD) :: FFLOC                   !< Fraction of LPOC in feces [-]
      real(RKD) :: FFDOC                   !< Fraction of DOC in feces [-]
      real(RKD) :: FFSUU                   !< Fraction of SUU in feces [-]
      real(RKD) :: UFROP                     !< Fraction of ROP in urine [-]
      real(RKD) :: UFLOP                     !< Fraction of LOP in urine [-]
      real(RKD) :: UFDOP                     !< Fraction of DOP in urine [-]
      real(RKD) :: UFP4D                     !< Fraction of P4D in urine [-]
      real(RKD) :: UFRON                     !< Fraction of RON in urine [-]
      real(RKD) :: UFLON                     !< Fraction of LON in urine [-]
      real(RKD) :: UFDON                     !< Fraction of DON in urine [-]
      real(RKD) :: UFNHX                     !< Fraction of NHX in urine [-]
      real(RKD) :: WQACF                        !< WQACF: Fraction of carbon in dry meat weight [-]
      real(RKD) :: WQANC                        !< WQANC: N/C ratio of shellfish [g N/g C]
      real(RKD) :: KACHC                        !< KACHC: Assimilation efficiency for cyanobacteria
      real(RKD) :: KACHD                        !< KACHD: Assimilation efficiency for diatom algae
      real(RKD) :: KACHG                        !< KACHG: Assimilation efficiency for green algae
      real(RKD) :: KAPOC                        !< KAPOC: Assimilation efficiency for POC      
      real(RKD) :: SPAWN_WEIGHT                 !< Min. dry weight for spawning [g C]
      real(RKD) :: SPAWN_TEMP                   !< Min. temperature for spawning [deg C]
      real(RKD) :: SPAWN_RATIO                  !< Fraction of cumulative reproductive biomass for spawing
      real(RKD) :: EGG_WEIGHT                   !< egg dry weight [g C]
      real(RKD) :: EGG_CAL                      !< Caloric content of egg [cal/g dry weight]
      real(RKD) :: SHF_CAL                      !< Caloric content of shellfish [cal/g dry weight]
      real(RKD) :: RESCAL                       !< Energy comsumption by respiration [cal/mL O2]
      real(RKD) :: EGG_HATCHED                  !< fraction of egg hatched to larvae
      real(RKD) :: EGG_TEMP1                    !< Min. temprature for spawned eggs to survive [deg. C]
      real(RKD) :: EGG_TEMP2                    !< Max. temprature for spawned eggs to survive [deg. C]
      real(RKD) :: EGG_SALT1                    !< Min. salinity for spawned eggs to survive [ppt]
      real(RKD) :: EGG_SALT2                    !< Max. salinity for spawned eggs to survive [ppt]
      real(RKD) :: LARVAE_RECRUITED             !< Fraction of larvae recruited per spawn [-]
      real(RKD) :: LARVAE_SURVIVAL              !< Fraction of larvae to spat survival [-]
      real(RKD) :: LARVAE_SPAN                  !< Larvae life span [days]
      real(RKD),allocatable :: SPAWN(:,:)       !< Spawn
      real(RKD),allocatable :: CRB(:,:)         !< Cumulative reproductive biomass [g C]
      real(RKD),allocatable :: PR(:,:)          !< Reproductive tissue production
      real(RKD),allocatable :: NP(:,:)          !< Net production
      real(RKD),allocatable :: CELL_C(:,:)      !< Biomass/carbon weight per grid cell layer [g C]
      real(RKD),allocatable :: INDI_C(:,:)      !< INDI_C(LCM,KC): carbon weight per individual in cell (Dry meat weight [g C])
      real(RKD),allocatable :: SHLEN(:,:)       !< Shell length [cm]
      real(RKD),allocatable :: INDI(:,:)        !< Number of shellfish individuals [#]
      real(RKD),allocatable :: INDI_D(:,:)      !< Number of death individuals [#]
      real(RKD),allocatable :: INDI_H(:,:)      !< Number of harvested individuals [#]
      real(RKD),allocatable :: HARV_C(:,:)      !< Harvested biomass per grid cell layer [g C]
      real(RKD),allocatable :: FR(:)            !< FR(LCM):  Filtration rate [1 filtered / ind. / hour]  
      real(RKD),allocatable :: BMG(:)           !< Basal metabolic rate [1/day]
      real(RKD),allocatable :: B_RESPI(:)       !< Basal metabolic rate [1/day]
      real(RKD),allocatable :: B_URINE(:) 
      real(RKD),allocatable :: B_GRAZI(:) 
      real(RKD),allocatable :: B_FECAL(:) 
      real(RKD),allocatable :: B_DEATH(:) 
      real(RKD),allocatable :: B_RPOC(:) 
      real(RKD),allocatable :: B_LPOC(:) 
      real(RKD),allocatable :: B_DOC(:) 
      real(RKD),allocatable :: AQ_RATE(:,:)     !< Harvest rate [-]
      real(RKD),allocatable :: WQAQ(:)          !< Dry meat weight per grid cell [g C]
      real(RKD),allocatable :: WQEA(:)          !< Shellfish individuals per grid cell [#]
      real(RKD),allocatable :: WQO(:)
      real(RKD),allocatable :: WQVO(:,:)
      real(RKD),allocatable :: KAALGAE(:)       !< Assimilation efficiency for algae                       
    end type 
    
    type SHELLFISHZONE
      integer :: ITYPE                          !< Type of aquaculture farm (not used)    
      real(RKD) :: AQ_DEP                       !< Depth of aquaculture farm [m]
      real(RKD),allocatable :: KDEATH(:)        !< Specific mortality rate [1/day]
    end type 
    
    integer :: NSF                              !< NSF: Number of shellfish classes (oyster, scallop, etc)
    integer :: NSFZONES                         !< NSF: Number of shellfish classes (oyster, scallop, etc)
    type(SHELLFISHTYPE),allocatable :: SF(:)    !< SF: A shellfish class (oyster, scallop, etc)
    type(SHELLFISHZONE),allocatable :: SFZONES(:)
    integer,allocatable :: FARMCELL(:)
    integer,allocatable :: FARMGR(:)            !< FARMGR(LCM): Type of shellfish farm (Zone ID)
    !REAL(RKD) :: TEMPB                          !< Base temperature [deg C]
    real(RKD),allocatable :: HARVTIMES(:),HARVDW(:,:)   !< mid-harvest
    logical,allocatable :: HARVDONE(:)
    integer :: ISFFARM,NSFCELLS,NSFHARV,ISHFWQ

  contains

  SUBROUTINE INIT_SHELLFISH(LCM,NK)
    integer, intent(IN) :: LCM,NK
    integer :: I,J,K,L
    
    allocate(SF(NSF))
    do I = 1,NSF
      allocate(SF(i).CELL_C(LCM,NK),SF(i).INDI_C(LCM,NK),SF(i).INDI(LCM,NK),SF(i).INDI_D(LCM,NK),SF(i).SHLEN(LCM,NK))
      allocate(SF(i).SPAWN(LCM,NK),SF(i).CRB(LCM,NK),SF(i).NP(LCM,NK),SF(i).PR(LCM,NK))
      allocate(SF(i).HARV_C(LCM,NK),SF(i).INDI_H(LCM,NK))
      allocate(SF(i).FR(LCM))
      allocate(SF(i).BMG(LCM),SF(i).B_RESPI(LCM),SF(i).B_GRAZI(LCM))      
      allocate(SF(i).B_FECAL(LCM),SF(i).B_DEATH(LCM),SF(i).B_URINE(LCM)) 
      allocate(SF(i).B_RPOC(LCM),SF(i).B_LPOC(LCM),SF(i).B_DOC(LCM)) 
      !ALLOCATE(SF(i).AQ_RATE(NSFHARV,LCM))
      allocate(SF(i).WQAQ(LCM),SF(i).WQEA(LCM))
      allocate(SF(i).WQO(LCM),SF(i).WQVO(LCM,NK))
      allocate(SF(i).KAALGAE(NALGAEM))
      
      do L = 1,LCM
        do K = 1,NK
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
        enddo
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
      enddo
      SF(i).KAALGAE = 0.
    enddo

    
    allocate(FARMCELL(LCM),FARMGR(LCM))    
    do L = 1,LCM
      FARMCELL(L) = 0
      do I = 1,NSF
        SF(i).WQAQ(L) = 0.
        SF(i).WQEA(L) = 0.
        do K = 1,KC
            SF(i).INDI_C(L,K) = 0.
            SF(i).INDI(L,K) = 0.
        enddo
      enddo
    enddo
  END SUBROUTINE 
  
  SUBROUTINE FREE_SHELLFISH()
    integer :: I

     do I = 1,NSF
      deallocate(SF(i).CELL_C,SF(i).INDI_C,SF(i).INDI,SF(i).INDI_D,SF(i).SHLEN)
      deallocate(SF(i).SPAWN,SF(i).CRB,SF(i).NP,SF(i).PR)
      deallocate(SF(i).HARV_C,SF(i).INDI_H)
      deallocate(SF(i).FR)
      deallocate(SF(i).BMG,SF(i).B_RESPI,SF(i).B_GRAZI)      
      deallocate(SF(i).B_FECAL,SF(i).B_DEATH,SF(i).B_URINE) 
      deallocate(SF(i).B_RPOC,SF(i).B_LPOC,SF(i).B_DOC) 
      deallocate(SF(i).AQ_RATE)
      deallocate(SF(i).WQAQ,SF(i).WQEA)
      deallocate(SF(i).WQO,SF(i).WQVO)
      deallocate(SF(i).KAALGAE)
    enddo
    do I = 1,NSFZONES
      deallocate(SFZONES(i).KDEATH)
    enddo
    deallocate(SF)
    deallocate(SFZONES)
    deallocate(HARVTIMES,HARVDW,HARVDONE)
    deallocate(FARMCELL,FARMGR)
  END SUBROUTINE 

  SUBROUTINE SHELLFISH_REDIST(L)
    real(RKD) :: AQDEP,ADEP,BDEP
    integer :: I,K,L,IDX

    if(FARMCELL(L) > 0 )then
      AQDEP = HP(L)
      IDX = FARMGR(L)          
      if(IDX > 0 .and. IDX <= NSFZONES )then
        AQDEP = SFZONES(IDX).AQ_DEP
      endif
      do I = 1,NSF
        if(AQDEP .LE. HP(L) )then
          ADEP = AQDEP
          BDEP = AQDEP
        else
          ADEP = HP(L)
          BDEP = HP(L)
        endif
        do K = KC,1,-1
          if( ADEP > 0 .and. ADEP >= HPK(L,K) )then
            if( SF(i).WQEA(L) > 0.0 )then
              SF(i).INDI(L,K) = SF(i).WQEA(L)*(HPK(L,K)/BDEP)       !!! WQEA(L) : vertical sum of ea of pre-step
              SF(i).CELL_C(L,K) = SF(i).WQAQ(L)*(HPK(L,K)/BDEP)     !!! WQAQ(L) : vertical sum of WQV(i) of pre-step
              SF(i).INDI_C(L,K) = SF(i).CELL_C(L,K)/SF(i).INDI(L,K)
            endif
            ADEP = ADEP-HPK(L,K)
          else
            if( SF(i).WQEA(L) > 0.0 )then
              SF(i).CELL_C(L,K) = (ADEP/BDEP)*SF(i).WQAQ(L)
              SF(i).INDI(L,K) = (ADEP/BDEP)*SF(i).WQEA(L)
              SF(i).INDI_C(L,K) = SF(i).CELL_C(L,K)/SF(i).INDI(L,K)
            endif
            ADEP = 0
            EXIT
          endif
        enddo
      enddo      
    endif
  END SUBROUTINE 

  SUBROUTINE SHELLFISH_LAYERSUM(L)
    real(RKD) :: TYAQ
    integer :: I,L,K

    do I = 1,NSF
      TYAQ = 0.0
      SF(i).WQAQ(L) = 0.
      SF(i).WQEA(L) = 0
      do K = KC,1,-1
        SF(i).WQAQ(L) = SF(i).WQAQ(L) + SF(i).CELL_C(L,K)
        SF(i).WQEA(L) = SF(i).WQEA(L) + SF(i).INDI(L,K)
      enddo
!      WQV(L,KC,18) = WQAQ1(L)
!      WQV(L,KC,1) = WQAQ2(L)
!      WQV(L,KC,3) = WQEA1(L)
      TYAQ = TYAQ + SF(i).WQAQ(L)
    enddo
  END SUBROUTINE 

  SUBROUTINE SHELLFISH_LENGTH(L)
    integer :: I,K,L
    if(FARMCELL(L) > 0 )then
      do I = 1,NSF
        do K = 1,KC
          SF(i).SHLEN(L,K) = SF(i).CWLEN*(SF(i).INDI_C(L,K)/SF(i).WQACF)**SF(i).PWLEN
        enddo
      enddo     
    endif
  END SUBROUTINE
  
  SUBROUTINE SHELLFISH_GROWTH(J,SF,L,K)
    use GLOBAL
    integer, intent(IN) :: J
    type(SHELLFISHTYPE), intent(INOUT) :: SF
    real(RKD) :: DW,KDEATH,INGEST,SFCHC,SFCHD,SFCHG,SFROC,SFLOC,FAC,WVOL,TSS,DOX
    real(RKD) :: SHFLEN,REFF,FRATIO,EGGS,LARVAE,RW,RT,GRAZIND,ASSIMIND,RESPIND,FDEATH
    real(RKD) :: SFALGAE(NALGAEM)
    integer :: L,K,IDX,NAL

    if( SF.CELL_C(L,K) > 0.0 )then
      
      DW = SF.INDI_C(L,K)/SF.WQACF          ! Biomass (gram C dry weight)
      
      SHFLEN = SF.CWLEN*DW**SF.PWLEN        ! Shellfish length (Kobayashi et al., 1997)
      if( SF.SHLEN(L,K) < SHFLEN) SF.SHLEN(L,K) = SHFLEN
      SHFLEN = SF.SHLEN(L,K)
      
      ! Max. filtration rate (liter filtered per ind per hour) 
      if( SF.IFILT == 0 )then
        ! Cerco and Noel (2007)
        SF.FR(L) = SF.FRMAX*(1000./24.)*DW 
      elseif( SF.IFILT == 2 )then
        ! Doering and Oviatt (1986)
        SF.FR(L) = (60./1000.)/2.95*(SHFLEN**0.96)*(TWQ(L)**0.95) 
      elseif( SF.IFILT == 3 )then
        ! Cloern (1982)
        SF.FR(L) = 7.0*(DW**0.67)
      elseif( SF.IFILT == 4 )then
        ! Coughlan & Ansell (1964)
        SF.FR(L) = 2.59*(DW**0.73)
      else ! Default SF.IFILT = 1
        ! Kobayashi et al. (1997)
        if( DW < 2.0 )then
          SF.FR(L) = 0.117*DW**3 - 1.05*DW**2 + 3.09*DW + 0.133
        else
          SF.FR(L) = 2.51*DW**0.279 
        endif
      endif
      
      ! Temperature effect
      if( SF.IFTEMP == 1 )then
        ! Kobayashi et al. (1997)
        if( TWQ(L) < 7.0 )then
          SF.FR(L) = SF.FR(L) * 0.59
        else
          SF.FR(L) = SF.FR(L) * (TWQ(L)**0.5)/4.47
        endif
      elseif( SF.IFTEMP == 2 )then
        ! Cerco and Noel (2007)
          SF.FR(L) = SF.FR(L) * EXP(-SF.KTG*(TWQ(L)-SF.TOPT)**2)
      endif
      
      ! Salinity effect
      if( SF.IFSALT == 1 )then
        ! Buzzelli et al. (2015)
        SF.FR(L) = SF.FR(L) * ((-0.0017*SWQ(L)+0.0084)*SWQ(L)-0.1002)
      elseif( SF.IFSALT == 2 )then
        ! Fulford (2007)
        SF.FR(L) = SF.FR(L) * (0.0926*SWQ(L)-0.139)
      elseif( SF.IFSALT == 3 )then
        ! Cerco and Noel (2005)
        SF.FR(L) = SF.FR(L) * 0.5*(1+TANH(SWQ(L)-SF.KHSOY))
      elseif( SF.IFSALT == 4 )then
        ! Quayle (1988) and Mann et al. (1991)
        if( SWQ(L) <= 10.0 )then
    	    SF.FR(L) = 0.0
        elseif( SWQ(L) < 20.0 )then
    	    SF.FR(L) = SF.FR(L) * (SWQ(L)-10.0)/10.0
        endif
      elseif( SF.IFSALT == 5 )then
        ! Loosanoff (1958)
        if( SWQ(L) <= 3.5 )then
    	    SF.FR(L) = 0.0
        elseif( SWQ(L) < 7.5 )then
    	    SF.FR(L) = SF.FR(L) * (SWQ(L)-3.5)/4.0
        endif
      endif

      ! TSS effect
      if( SF.IFTSS > 1 )then
        TSS = 0.
        if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
          TSS = SEDT(L,K)+SNDT(L,K) 
        endif
        if( TSS > 0 )then
          if( SF.IFTSS == 1 )then
            ! Fulford (2007)
            if( TSS > 25 )then
              SF.FR(L) = SF.FR(L) * (10.364*(LOG(TSS))**(-2.0477))   
            endif
          elseif( SF.IFTSS == 2 )then
            ! Cerco and Noel (2005)
            if( TSS > 100.0 )then
    	        SF.FR(L) = 0.0
            elseif( TSS > 25. )then
    	        SF.FR(L) = SF.FR(L) * 0.2
            elseif( TSS < 5. )then
    	        SF.FR(L) = SF.FR(L) * 0.1
            endif
          elseif( SF.IFTSS == 3 )then
            ! Hofmann et al. (1992)
            TSS = TSS/1000. ! Convert mg/L to g/L
            SF.FR(L) = SF.FR(L) * (1. - 0.01*(LOG10(TSS) + 3.38)/0.418)   
          endif
        endif
      endif
      
      ! DO effect
      if( SF.IFDO == 1 .and. ISKINETICS(IDOX) > 0 )then
        DOX = WQV(L,K,IDOX)
        ! Cerco and Noel (2005)
        SF.FR(L) = SF.FR(L) / (1.  + EXP(1.1*(SF.DOHX-DOX)/(SF.DOHX-SF.DOQX)))   
      endif
      
      INGEST = SF.FR(L) * 24.0/1000.0 * DTWQ          ! Ingestion [L]
      if( ISWQLVL == 0 )then
        SFCHC =  INGEST * SF.KACHC*WQV(L,K,1) 
        SFCHD =  INGEST * SF.KACHD*WQV(L,K,2) 
        SFCHG =  INGEST * SF.KACHG*WQV(L,K,3) 
      else
        do NAL = 1,NALGAE
          SFALGAE(NAL) = INGEST * SF.KAALGAE(NAL)*WQV(L,K,19+NAL)
        enddo
      endif
      
      SFROC =  INGEST * SF.KAPOC*WQV(L,K,IROC) 
      SFLOC =  INGEST * SF.KAPOC*WQV(L,K,ILOC) 
      if(ISWQLVL == 0 )then
        GRAZIND = SFCHC + SFCHD + SFCHG + SFROC + SFLOC ! Grazing [g C]
      else
        GRAZIND = SFROC + SFLOC
        do NAL = 1,NALGAE
          GRAZIND = GRAZIND +  SFALGAE(NAL)
        enddo
      endif
      ASSIMIND = SF.K_ASSIM * GRAZIND                 ! Assimilation [g C]
      
      ! Respiration
      if( SF.IRES == 1 )then
        if( TWQ(L) > SF.TEMP_RESPI )then
          SF.BMG(L) = SF.K_RESPI * SF.THETA_RESPI**(TWQ(L)-SF.TEMP_RESPI)
        else
          SF.BMG(L) = 0.
        endif
        RESPIND = SF.INDI_C(L,K) * SF.BMG(L)*DTWQ
      elseif( SF.IRES == 2 )then
        ! Cerco and Noel (2005)
        SF.BMG(L) = SF.K_RESPI * EXP(SF.THETA_RESPI*(TWQ(L)-SF.TEMP_RESPI))
        RESPIND = SF.INDI_C(L,K) * SF.BMG(L)*DTWQ
      else  ! Default SF.IRES == 0
        RW = (12.6*TWQ(L) + 69.7)*DW**-0.25       ! Dame (1972)
        !RW = (31.0*TWQ(L) - 22.0)*DW**-0.3
        RT = 1.0
        if( SF.IRTEMP == 1 )then
          ! Shumway and Koehn (1992)
          if( TWQ(L) < 20.0 )then
              RT = (0.007*TWQ(L) + 2.099)
          else
              RT = (0.0915*TWQ(L) + 1.324)
          endif
        endif
        if( SF.IRSALT == 1 )then
          ! Shumway and Koehn (1982)
          if( SWQ(L) <= 15.0 )then
            RW = RW*RT
          elseif( SWQ(L) < 20.0 )then
    	      RW = RW*(1.0+(((RT-1.0)/5.0)*(20.0-SWQ(L))))
          endif
        elseif( SF.IRSALT == 2 )then
          ! Hoffman et al. (1992)
          if( SWQ(L) <= 10.0 )then
            RW = RW*RT
          elseif( SWQ(L) < 15.0 )then
    	      RW = RW*(1.0+(((RT-1.0)/5.0)*(15.0-SWQ(L))))
          endif
        endif
        RESPIND = RW * (SF.RESCAL/SF.SHF_CAL) * 24.0/1000.0 * DTWQ
      endif

      SF.NP(L,K) = ASSIMIND - RESPIND           ! Net production

      SF.PR(L,K) = 0.0
      SF.SPAWN(L,K) = 0.0
      if( DW > SF.SPAWN_WEIGHT )then            ! Adult weight
        if( SF.NP(L,K) > 0.0 )then
          ! Kusaka et al. (1991)
          if( TWQ(L) >= 27.0  )then
            REFF = 0.8
          elseif( TWQ(L) > 23.0 )then
            REFF = 0.2*TWQ(L) - 4.6
          else
            REFF = 0.0
          endif
          SF.PR(L,K) = REFF * SF.NP(L,K)            ! Reproductive tissue production
          SF.CRB(L,K) = SF.CRB(L,K) + SF.PR(L,K)    ! Cumulative reproductive biomass
        elseif( SF.NP(L,K) < 0.0 )then
          if( SF.CRB(L,K) > ABS(SF.NP(L,K)) )then
            SF.CRB(L,K) = SF.CRB(L,K) + SF.NP(L,K)
            SF.NP(L,K) = 0.0
          else
            SF.NP(L,K) = SF.NP(L,K) + SF.CRB(L,K)
            SF.CRB(L,K) = 0.0
          endif
        endif
               
        ! Reproduction
        if( SF.ISPAWN > 0 )then
          if( (TWQ(L) >= SF.SPAWN_TEMP) .and. (SF.CRB(L,K) >= SF.SPAWN_RATIO*DW) )then
        	FRATIO = 0.021*10*SHFLEN - 0.62          ! Kennedy (1982)
            SF.SPAWN(L,K) = FRATIO/(FRATIO+1) * SF.CRB(L,K)
            EGGS = SF.SPAWN(L,K)*SF.SHF_CAL/(SF.EGG_CAL*SF.EGG_WEIGHT)
            if( TWQ(L) >= SF.EGG_TEMP1 .and. TWQ(L) <= SF.EGG_TEMP2 .and. & 
                SWQ(L) >= SF.EGG_SALT1 .and. SWQ(L) <= SF.EGG_SALT2 )then
              !LARVAE(L,K) = SF.EGG_HATCHED*EGGS  ! Number of shellfish lavae
              !LARVAE_LIFE(L,K) = 0.0             ! A larvae life starts
            else
              !LARVAE(L,K) = 0.0
            endif
            DW = DW - SF.SPAWN(L,K)
            SF.CRB(L,K) = 0.0
          endif
        endif
      endif      
      
      DW = DW + SF.NP(L,K)
      SF.INDI_C(L,K) = DW * SF.WQACF
            
      IDX = FARMGR(L)          
      if(IDX > 0 .and. IDX <= NSFZONES )then
        KDEATH = SFZONES(IDX).KDEATH(J)
      else
        KDEATH = SF.K_DEATH
      endif
      SF.B_DEATH(L) = SF.CELL_C(L,K)*KDEATH*DTWQ                          
      SF.INDI_D(L,K) = SF.INDI(L,K)*KDEATH*DTWQ
      if( SF.ITDEATH > 0 )then
        FDEATH = SF.THETA_DEATH**(TWQ(L)-SF.TEMP_DEATH)
        SF.B_DEATH(L) = SF.B_DEATH(L) * FDEATH
        SF.INDI_D(L,K) = SF.INDI_D(L,K) * FDEATH
      endif

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
      if( ISHFWQ > 0 .and. HPK(L,K) > 0 )then
        WVOL = DXYP(L) * HPK(L,K)                 ! water volume [m3]
        FAC = 1.0  / WVOL                         ! Converts [g C] to [mg/L]
        if(ISWQLVL == 0 )then
          if( ISKINETICS(1) > 0) WQV(L,K,1) = WQV(L,K,1) - FAC*SFCHC*SF.INDI(L,K)
          if( ISKINETICS(2) > 0) WQV(L,K,2) = WQV(L,K,2) - FAC*SFCHD*SF.INDI(L,K)
          if( ISKINETICS(3) > 0) WQV(L,K,3) = WQV(L,K,3) - FAC*SFCHG*SF.INDI(L,K)
        else
          do NAL = 1,NALGAE
            WQV(L,K,19+NAL) = WQV(L,K,19+NAL) - FAC*SFALGAE(NAL)*SF.INDI(L,K)
          enddo
        endif
        
        if( ISKINETICS(IROC) > 0) WQV(L,K,IROC) = WQV(L,K,IROC) + FAC*(SF.B_RPOC(L) - SFROC*SF.INDI(L,K))
        if( ISKINETICS(ILOC) > 0) WQV(L,K,ILOC) = WQV(L,K,ILOC) + FAC*(SF.B_LPOC(L) - SFLOC*SF.INDI(L,K))
        if( ISKINETICS(IDOC) > 0) WQV(L,K,IDOC) = WQV(L,K,IDOC) + FAC*(SF.B_DOC(L))
      
        if( ISKINETICS(IROP) > 0) WQV(L,K,IROP) = WQV(L,K,IROP) + FAC*(SF.UFROP * SF.B_URINE(L))
        if( ISKINETICS(ILOP) > 0) WQV(L,K,ILOP) = WQV(L,K,ILOP) + FAC*(SF.UFLOP * SF.B_URINE(L))
        if( ISKINETICS(IDOP) > 0) WQV(L,K,IDOP) = WQV(L,K,IDOP) + FAC*(SF.UFDOP * SF.B_URINE(L))
        if( ISKINETICS(IP4D) > 0) WQV(L,K,IP4D) = WQV(L,K,IP4D) + FAC*(SF.UFP4D * SF.B_URINE(L))
        if( ISKINETICS(IRON) > 0) WQV(L,K,IRON) = WQV(L,K,IRON) + FAC*(SF.UFRON * SF.B_URINE(L))
        if( ISKINETICS(ILON) > 0) WQV(L,K,ILON) = WQV(L,K,ILON) + FAC*(SF.UFLON * SF.B_URINE(L))
        if( ISKINETICS(IDON) > 0) WQV(L,K,IDON) = WQV(L,K,IDON) + FAC*(SF.UFDON * SF.B_URINE(L))
        if( ISKINETICS(INHX) > 0) WQV(L,K,INHX) = WQV(L,K,INHX) + FAC*(SF.UFNHX * SF.B_URINE(L))
      
        if( ISKINETICS(ISUU) > 0) WQV(L,K,ISUU) = WQV(L,K,ISUU) + FAC*(SF.FFSUU * SF.B_FECAL(L))
        if( ISKINETICS(IDOX) > 0) WQV(L,K,IDOX) = WQV(L,K,IDOX) - FAC*RESPIND*SF.INDI(L,K)
      endif
    endif
  END SUBROUTINE
  
  SUBROUTINE SHELLFISH_HARVEST(LF,LL)
    integer, intent(IN) :: LF,LL
    integer :: I,J,K,L
    real :: TIME
    
    TIME = TIMESEC/TCON  
    
    do J = 1,NSFHARV
      if((TIME >= HARVTIMES(J)) .and. (.not. HARVDONE(J)) )then
        do L = LF,LL
          do I = 1,NSF
            SF(i).WQAQ(L) = 0.
            SF(i).WQEA(L) = 0
            do K = KC,1,-1
              if( SF(i).CELL_C(L,K) >= 0.0  )then
                SF(i).INDI_C(L,K) = SF(i).CELL_C(L,K)/(SF(i).INDI(L,K) + 1.E-18)         ! Individual weight
                if( SF(i).INDI_C(L,K) >= HARVDW(I,J) )then                    ! Is greater than or equal to harvest size?
                  SF(i).INDI_H(L,K) = SF(i).AQ_RATE(J,L)*SF(i).INDI(L,K)      ! Number of harvested individuals
                  SF(i).HARV_C(L,K) = SF(i).INDI_H(L,K)*SF(i).INDI_C(L,K)     ! Harvested biomass
                  SF(i).INDI(L,K) = SF(i).INDI(L,K) - SF(i).INDI_H(L,K)       ! Remaining individuals
                  SF(i).INDI_C(L,K) = SF(i).INDI_C(L,K) - SF(i).HARV_C(L,K)   ! Remaining biomass
                endif
                SF(i).WQAQ(L) = SF(i).WQAQ(L) + SF(i).CELL_C(L,K)
                SF(i).WQEA(L) = SF(i).WQEA(L) + SF(i).INDI(L,K)
              endif
            enddo
          enddo
        enddo
        HARVDONE(J) = .TRUE.
      endif
    enddo
  END SUBROUTINE
  
  SUBROUTINE READ_SHELLFISH_JSON(LCM,NK)
    use fson
    use mod_fson_value, only: fson_value_count, fson_value_get
      
    type(fson_value), pointer :: json_data, items, item
    character*80 :: STR*200
    integer :: I,J,K,L,LCM,NK,IX,JX,IDX,NAL
    real :: ADEP,AQDEP
    real,allocatable :: ARR1D(:),ARR2D(:,:),SFWEIGHT(:)

    json_data => fson_parse("shellfish.jnp")   
    call fson_get(json_data, "number_of_species", NSF)
    call fson_get(json_data, "wq_interact_flag", ISHFWQ)    
    
    call INIT_SHELLFISH(LCM,NK)
    
    call fson_get(json_data, "species", items)
    do i = 1, fson_value_count(items) !NSF
      item => fson_value_get(items, i)
      !CALL fson_get(item, "name", SF(i).NAME)
      !CALL fson_get(item, "index", SF(i).IDX)
        
      ! *** General parameters
      call fson_get(item, "stoichiometry.nitrogen_to_carbon_ratio", SF(i).WQANC)
      call fson_get(item, "dry_meat_carbon", SF(i).WQACF)
      call fson_get(item, "dry_meat_caloric", SF(i).SHF_CAL)        
      call fson_get(item, "length_weight_factor", SF(i).CWLEN)
      call fson_get(item, "length_weight_exponent", SF(i).PWLEN)        

      call fson_get(item, "filtration.max_filtration_opt", SF(i).IFILT)
      call fson_get(item, "filtration.temperature_opt", SF(i).IFTEMP)
      call fson_get(item, "filtration.salinity_opt", SF(i).IFSALT)
      call fson_get(item, "filtration.sediment_opt", SF(i).IFTSS)
      call fson_get(item, "filtration.dissolved_oxygen_opt", SF(i).IFDO)
      call fson_get(item, "filtration.max_filtration_rate", SF(i).FRMAX)
      call fson_get(item, "filtration.temperature_optimal", SF(i).TOPT)
      call fson_get(item, "filtration.temperature_effect", SF(i).KTG)
      call fson_get(item, "filtration.halved_rate_salinity", SF(i).KHSOY)
      call fson_get(item, "filtration.one_half_rate_DO", SF(i).DOHX)
      call fson_get(item, "filtration.one_fourth_rate_DO", SF(i).DOQX)

      call fson_get(item, "respiration.respiration_opt", SF(i).IRES)
      call fson_get(item, "respiration.temperature_opt", SF(i).IRTEMP)
      call fson_get(item, "respiration.salinity_opt", SF(i).IRSALT)
      call fson_get(item, "respiration.sediment_opt", SF(i).IRTSS)
      call fson_get(item, "respiration.dissolved_oxygen_opt", SF(i).IRDO)
      call fson_get(item, "respiration.basal_metabolic_rate", SF(i).K_RESPI)
      call fson_get(item, "respiration.ref_temperature", SF(i).TEMP_RESPI)
      call fson_get(item, "respiration.temperature_coeff", SF(i).THETA_RESPI)
      call fson_get(item, "respiration.energy_consumption", SF(i).RESCAL)
      
      allocate(ARR1D(3))
      call fson_get(item, "growth.assimilation_of_algaes", ARR1D)
      call fson_get(item, "growth.assimilation_of_POC", SF(i).KAPOC)
      call fson_get(item, "growth.assimilation_efficiency", SF(i).K_ASSIM)
      call fson_get(item, "growth.growth_efficiency", SF(i).K_GROEF)
      if(ISWQLVL == 0 )then
        SF(i).KACHC = ARR1D(1)
        SF(i).KACHD = ARR1D(2)
        SF(i).KACHG = ARR1D(3)
      else
        do NAL = 1,NALGAE
          SF(i).KAALGAE(NAL) = ARR1D(NAL)
        enddo
      endif
      deallocate(ARR1D)
      
      call fson_get(item, "excretion.fraction_of_carbon.RPOC", SF(i).FFROC)
      call fson_get(item, "excretion.fraction_of_carbon.LPOC", SF(i).FFLOC)
      call fson_get(item, "excretion.fraction_of_carbon.DOC", SF(i).FFDOC)
      call fson_get(item, "excretion.fraction_of_phosphorus.RPOP", SF(i).UFROP)
      call fson_get(item, "excretion.fraction_of_phosphorus.LPOP", SF(i).UFLOP)
      call fson_get(item, "excretion.fraction_of_phosphorus.DOP", SF(i).UFDOP)
      call fson_get(item, "excretion.fraction_of_phosphorus.PO4", SF(i).UFP4D)
      call fson_get(item, "excretion.fraction_of_nitrogen.RPON", SF(i).UFRON)
      call fson_get(item, "excretion.fraction_of_nitrogen.LPON", SF(i).UFLON)
      call fson_get(item, "excretion.fraction_of_nitrogen.DON", SF(i).UFDON)
      call fson_get(item, "excretion.fraction_of_nitrogen.NH4", SF(i).UFNHX)
      call fson_get(item, "excretion.fraction_of_silica.SU", SF(i).FFSUU)
      
      call fson_get(item, "mortality.specific_rate", SF(i).K_DEATH)
      call fson_get(item, "mortality.temperature_dependent_flag", SF(i).ITDEATH)
      call fson_get(item, "mortality.ref_temperature", SF(i).TEMP_DEATH)      
      call fson_get(item, "mortality.temperature_coeff", SF(i).THETA_DEATH)      
      call fson_get(item, "mortality.dead_shellfish.RPOC", SF(i).DFROC)
      call fson_get(item, "mortality.dead_shellfish.LPOC", SF(i).DFLOC)
      call fson_get(item, "mortality.dead_shellfish.DOC", SF(i).DFDOC)
     
      call fson_get(item, "spawning.spawning_flag", SF(i).ISPAWN)
      call fson_get(item, "spawning.min_dry_weight", SF(i).SPAWN_WEIGHT)
      call fson_get(item, "spawning.min_temperature", SF(i).SPAWN_TEMP)
      call fson_get(item, "spawning.min_cumulative_reproductive_biomass", SF(i).SPAWN_RATIO)
     
      call fson_get(item, "eggs.dry_weight", SF(i).EGG_WEIGHT)
      call fson_get(item, "eggs.caloric_content", SF(i).EGG_CAL)
      call fson_get(item, "eggs.min_temp_to_survive", SF(i).EGG_TEMP1)
      call fson_get(item, "eggs.max_temp_to_survive", SF(i).EGG_TEMP2)
      call fson_get(item, "eggs.min_salt_to_survive", SF(i).EGG_SALT1)
      call fson_get(item, "eggs.max_salt_to_survive", SF(i).EGG_SALT2)

      call fson_get(item, "larvae.fraction_hatched", SF(i).EGG_HATCHED)
      call fson_get(item, "larvae.fraction_to_spat_survival", SF(i).LARVAE_SURVIVAL)
      call fson_get(item, "larvae.fraction_recruited_per_spawn", SF(i).LARVAE_RECRUITED)
      call fson_get(item, "larvae.life_span", SF(i).LARVAE_SPAN)

      ! *** Initial condition
      call fson_get(item, "initial_conditions",  SF(i).IC)    
    enddo
    
    json_data => fson_parse("shffarm.jnp")   
    call fson_get(json_data, "number_of_zones", NSFZONES)
    call fson_get(json_data, "number_of_cells", NSFCELLS)    
    call fson_get(json_data, "number_of_harvests", NSFHARV)   
    
    allocate(SFZONES(NSFZONES))
    do J = 1,NSFZONES
      allocate(SFZONES(J).KDEATH(NSF))
    enddo
    allocate(HARVTIMES(NSFHARV),HARVDONE(NSFHARV))
    allocate(HARVDW(NSF,NSFHARV))
    do J = 1,NSFHARV
      HARVDONE(J) = .FALSE.
    enddo
    do i = 1,NSF
      allocate(SF(i).AQ_RATE(NSFHARV,LCM))
    enddo
    allocate(ARR1D(NSF), SFWEIGHT(NSF))    
    
    call fson_get(json_data, "shellfish_harvest.mid_days", HARVTIMES)   
    call fson_get(json_data, "shellfish_harvest.min_sizes", ARR2D)
    if( .not. allocated(ARR2D) )then
      allocate(ARR2D(NSF,NSFHARV))
      ARR2D = 0.
    endif
    
    do i = 1, NSF
        do j = 1, NSFHARV
            HARVDW(i, j) = ARR2D(i, j)
        enddo
    enddo

    call fson_get(json_data, "zones", items)
    do i = 1, fson_value_count(items) !NSFZONES
      item => fson_value_get(items, i)
      call fson_get(item, "Type", SFZONES(i).ITYPE)
      call fson_get(item, "installed_depth", SFZONES(i).AQ_DEP)
      call fson_get(item, "death_rate", SFZONES(i).KDEATH)
      !DO j = 1,NSF
      !  SFZONES(i).KDEATH(j) = VALUES(j)
      !ENDDO
    enddo
    
    call fson_get(json_data, "cells", items)
    do k = 1, fson_value_count(items) !NSFCELLS
      item => fson_value_get(items, k)
      call fson_get(item, "address.L", L)
      !CALL fson_get(item, "address.I", IX)
      !CALL fson_get(item, "address.J", JX)
      call fson_get(item, "zone", FARMGR(L))
      call fson_get(item, "initial_conditions.individuals", ARR1D)
      call fson_get(item, "initial_conditions.weights", SFWEIGHT)
      call fson_get(item, "harvest_rates", ARR2D)
      if( .not. allocated(ARR2D) )then
        allocate(ARR2D(NSFHARV,LCM))
        ARR2D = 0.
      endif
      FARMCELL(L) = 1
      
      ADEP = HP(L)
      IDX = FARMGR(L)          
      if(IDX > 0 .and. IDX <= NSFZONES )then
        AQDEP = SFZONES(IDX).AQ_DEP
      else
        write(*,*) 'INVALID SHELLFISH ZONE NUMBER FOR CELL',L
      endif
       do i = 1,NSF
        do j = 1,NSFHARV
          SF(i).AQ_RATE(J,L) = ARR2D(i, j)
        enddo
        
        if(AQDEP <= HP(L) )then
          ADEP = AQDEP
        else
          ADEP = HP(L)
        endif
        
        SF(i).IC = SFWEIGHT(i)
        SF(i).WQEA(L) = ARR1D(i)  !*DXYP(L)*ADEP*SF(i).DRATIO
        if( ICONTINUE /= 1 )then      
          SF(i).WQAQ(L) = SF(i).IC*SF(i).WQEA(L)
        else
          !WQAQ1(L) = WQV(L,3,18)
          !WQAQ2(L) = WQV(L,3,1)
          !WQEA1(L) = WQV(L,3,3)
          !WQEA2(L) = WQEA2(L)
        endif
      enddo
    enddo
    deallocate(ARR1D,SFWEIGHT)
  END SUBROUTINE 
  
END MODULE
