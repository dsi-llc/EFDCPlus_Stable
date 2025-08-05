! ----------------------------------------------------------------------
!   This file is a part of EFDC + 
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE BEDINIT
  !
  ! *** SUBROUTINE BEDINIT INITIALIZES SEDIMENT AND TOXIC VARIABLES
  ! *** IN SEDIMENT BED FOR HOT AND COLD START CONDITIONS
  !     CHANGE RECORD
  !     ADDED ADDITIONAL DIAGNOSTIC OUTPUT
  !     MOVED TOXIC INITIALIZATIONS FROM SSEDTOX
  !

  use GLOBAL
  use Variables_MPI
  use MPI
  use Broadcast_Routines

  use INFOMOD,only:SKIPCOM,READSTR

  implicit none

  integer :: K, L, LP, NS, NX, NT, KTOPP1, IVAL, KTOPTP, IHOTSTRT
  integer :: ISSTYPE, LD, ID, JD, LCORE, IDUM, NL, KT1, KT2
  integer :: ISO, NSITES, NCOHSEDS, NCOHSEDL
  character :: STR*200

  integer,allocatable,dimension(:) :: LSSCOHSED


  real :: CSEDTAUS, CSEDRESS, CSEDTAUB, CSEDRESB, TMPCVT
  real :: SURF, FVOLSSD, FVOLSED, FVOLSND, TMP1, SXD
  real :: RMULADJ, RUMLADJC, FRACT1, FRACT2, FRACT2C, HBEDP, TTHICK
  real(RKD) :: DTOTAL
  
  real,allocatable,dimension(:) :: FRACACT
  real,allocatable,dimension(:) :: FRACPAR
  real,allocatable,dimension(:) :: RADJCOHSEDS
  real,allocatable,dimension(:) :: SUMSED
  real,allocatable,dimension(:) :: TAUDSS
  real,allocatable,dimension(:,:) :: TAURSS
  real,allocatable,dimension(:,:) :: TAUNSS
  real,allocatable,dimension(:,:) :: TEXPSS
  real,allocatable,dimension(:,:) :: DEPBBSS
  real,allocatable,dimension(:,:) :: WRSPOSS
  real,allocatable,dimension(:,:) :: SEDBALL

  ! *** MPI
  integer,allocatable,dimension(:) :: LSSCOHSED_Global
  real,allocatable,dimension(:) :: RADJCOHSEDS_Global
  integer :: l_local
  
  allocate(LSSCOHSED(LCM))
  allocate(FRACACT(LCM))
  allocate(FRACPAR(LCM))
  allocate(RADJCOHSEDS(LCM))
  allocate(SEDBALL(LCM,KBM))
  allocate(SUMSED(KBM))

  ! *** For MPI
  allocate(LSSCOHSED_Global(LCM_Global))
  allocate(RADJCOHSEDS_Global(LCM_Global))
  
  LSSCOHSED = 0
  FRACACT = 0.0
  FRACPAR = 0.0
  RADJCOHSEDS = 0.0
  SEDBALL = 0.0
  SUMSED  = 0.0

  ! *** SEDIMENT PROPERTIES ARE IN MASS PER UNIT AREA
  ! ***
  ! *** SEDB - COHESIVE SEDIMENTS (G/M2)
  ! *** SNDB - NONCOHESIVE SEDIMENTS (G/M2)
  ! *** TOXB - TOXICS (MG/M2)

  ! *** ZERO LOCAL ARRAYS
  FRACACT = 0.0
  FRACPAR = 0.0

  ! SET BOTTOM LAYER NUMBER: ORIGINAL:KBB = 1, SEDZLJ:KBB = KB
  if( LSEDZLJ )then
    KBB = KB
  else
    KBB = 1
  endif

  ! *** DETERMINE START UP MODE
  IHOTSTRT = 0
  if( ISRESTI /= 0 )then
    if( ISCI(6) /= 0 .or. ISCI(7) /= 0 )then
      IHOTSTRT = 1
    endif
  endif
  
  ! ***********************************************************************************
  ! *** HOT START INITIALIZATION.  SEDZLJ HOTSTART IS HANDLED IN S_SEDIC
  if( IHOTSTRT /= 0 .and. .not. LSEDZLJ )then

    ! *** SET POROSITY
    do K = 1,KB
      do L = 2,LA
        PORBED(L,K) = VDRBED(L,K)/(1.0 + VDRBED(L,K))
        PORBED1(L,K) = VDRBED1(L,K)/(1.0 + VDRBED1(L,K))
      enddo
    enddo

    ! *** SET BULK DENSITY - CURRENT TIME
    do K = 1,KB
      do L = 2,LA
        SEDBT(L,K) = 0.0
        SNDBT(L,K) = 0.0
      enddo
    enddo
    if( ISTRAN(6) >= 1 )then
      do NS = 1,NSED
        do K = 1,KB
          do L = 2,LA
            SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
          enddo
        enddo
      enddo
    endif
    if( ISTRAN(7) >= 1 )then
      do NX = 1,NSND
        do K = 1,KB
          do L = 2,LA
            SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NX)
          enddo
        enddo
      enddo
    endif

    ! *** COMPUTE TOTAL/WET DENSITY
    do K = 1,KB
      do L = 2,LA
        if( HBED(L,K) > 0.0 )then
          BDENBED(L,K) = 1000.0*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)
        else
          BDENBED(L,K) = 0.0
        endif
      enddo
    enddo

    ! *** SET BULK DENSITY - PREVIOUS TIME
    do K = 1,KB
      do L = 2,LA
        SEDBT(L,K) = 0.0
        SNDBT(L,K) = 0.0
      enddo
    enddo
    if( ISTRAN(6) >= 1 )then
      do NS = 1,NSED
        do K = 1,KB
          do L = 2,LA
            SEDBT(L,K) = SEDBT(L,K) + SEDB1(L,K,NS)
          enddo
        enddo
      enddo
    endif
    if( ISTRAN(7) >= 1 )then
      do NX = 1,NSND
        do K = 1,KB
          do L = 2,LA
            SNDBT(L,K) = SNDBT(L,K) + SNDB1(L,K,NX)
          enddo
        enddo
      enddo
    endif
    do K = 1,KB
      do L = 2,LA
        if( HBED1(L,K) > 0.0 )then
          BDENBED1(L,K) = 1000.0*PORBED1(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED1(L,K)
        else
          BDENBED1(L,K) = 0.0
        endif
      enddo
    enddo

    ! *** SET TOP BED LAYER
    do L = 2,LA
      KBT(L) = 1
      do K = 1,KB
        if( HBED(L,K) <= 0.0 )then
          KBT(L) = max(1,K - 1)
          exit
        endif
      enddo
    enddo
    
  endif
  ! *** END HOT START INITIALIZATION

  ! ***********************************************************************************
  ! ***  COLD START INITIALIZATION: SPATIALLY UNIFORM BED
  if( IHOTSTRT == 0 .and. ISEDINT <= 1 .and. .not. LSEDZLJ )then
    ! ***  SET POROSITY AND VOID RATIO

    ! *** ORIGINAL SEDIMENT MODULE
    HBED = 0.0
   
    ! ***  UNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS  
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY  
    if( ISTRAN(6) >= 1 )then
      do NS = 1,NSED
        do K = 1,KB
          do L = 2,LA
            HBED(L,K) = HBED(L,K) + SDEN(NS)*SEDB(L,K,NS)
          enddo
        enddo
      enddo
    endif
    if( ISTRAN(7) >= 1 )then
      do NX = 1,NSND
        NS = NSED+NX
        do K = 1,KB
          do L = 2,LA
            HBED(L,K) = HBED(L,K) + SDEN(NS)*SNDB(L,K,NX)
          enddo
        enddo
      enddo
    endif
  
    ! *** SET TOP LAYER
    do L = 2,LA
      KBT(L) = 1
      do K = 1,KB
        if( HBED(L,K) <= 0.0 )then
          KBT(L) = max(1,K - 1)
          exit
        endif
      enddo
    enddo

    ! *** HANDLE ACTIVE ARMORING LAYER IF NOT PRESENT IN INITIAL OR RESTART
    if( ISNDAL == 2 .and. IALSTUP > 0 .and. KB > 1 )then
      do L = 2,LA
        KBT(L)  = max(1,KBT(L) - 1)
        SEDB(L,KBT(L)+1,:) = 0.0
        SNDB(L,KBT(L)+1,:) = 0.0
        HBED(L,KBT(L)+1) = 0.0
      enddo
    endif

  endif
  ! *** END COLD START INITIALIZATION: SPATIALLY UNIFORM BED

  ! ***********************************************************************************
  ! *** COLD START INITIALIZATION: SPATIALLY VARYING BED
  if( IHOTSTRT == 0 .and. (ISEDINT > 1 .or. NSEDFLUME == 3) )then

    ! *** INITIALIZE VOID RATIO AND POROSITY
    if( IBMECH > 0 .or. LSEDZLJ )then
      do K = 1,KB
        do L = 2,LA
          VDRBED(L,K) = BEDDINIT(L,K)
          PORBED(L,K) = VDRBED(L,K)/(1. + VDRBED(L,K))
        enddo
      enddo
    endif
            
    ! *** INITIAL CONDITIONS IN SEDB.INP 
    if( ISEDBINT == 0 )then
      ! *** SEDIMENT INITIALIZATION: MASS per UNIT AREA
      HBED = 0.0
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do L = 2,LA
              HBED(L,K) = HBED(L,K) + SDEN(NS)*SEDB(L,K,NS)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do L = 2,LA
              HBED(L,K) = HBED(L,K) + SDEN(NS)*SNDB(L,K,NX)
            enddo
          enddo
        enddo
      endif
      
      ! *** ADJUST LAYER THICKNESSES FOR VOOID RATIO
      do K = 1,KB
        do L = 2,LA
          HBED(L,K) = (1. + VDRBED(L,K))*HBED(L,K)
        enddo
      enddo
      
      !  do K = 1,KB
      !    do L = 2,LA
      !      if( HBED(L,K) > 0.0 )then
      !        ! *** COMPUTE TOTAL/WET DENSITY
      !        BDENBED(L,K) = 1000.*PORBED(L,K)  +0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)
      !      else
      !        BDENBED(L,K) = 0.0
      !      endif
      !    enddo
      !  enddo
      !  do K = 1,KB
      !    do L = 2,LA
      !      HBED1(L,K) = HBED(L,K)
      !      BDENBED1(L,K) = BDENBED(L,K)
      !    enddo
      !  enddo
    endif     ! *** END OF SEDIMENT INITIALIZATION: MASS per UNIT AREA
    
    ! *** INITIAL CONDITIONS IN SEDB.INP 
    if( ISEDBINT == 1 )then
      ! *** SEDIMENT INITIALIZATION: MASS FRACTION
      ! *** THIS OPTION REQUIRES INITIAL LAYER THICKNESSES


      ! *** INITIALIZE BED LAYER THICKNESSES
      do K = 1,KB
        do L = 2,LA
          HBED(L,K) = BEDLINIT(L,K)
        enddo
      enddo
    
      ! *** FIRST COMPUTE DRY DENSITY (G/M2)
      BDENBED = 0.0
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do L = 2,LA
              BDENBED(L,K) = BDENBED(L,K) + 1000.0*SSG(NS)*SEDB(L,K,NS)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do L = 2,LA
              BDENBED(L,K) = BDENBED(L,K) + 1000.0*SSG(NS)*SNDB(L,K,NX)
            enddo
          enddo
        enddo
      endif

      ! *** CONVERT DRY DENSITY TO TOTAL/WET DENSITY
      do K = 1,KB
        do L = 2,LA
          BDENBED(L,K) = 1000.0*PORBED(L,K) + (1.0 - PORBED(L,K))*BDENBED(L,K)
        enddo
      enddo

      ! *** TOTAL SEDIMENT FOR ALL CLASSES (G/M2)
      do K = 1,KB
        do L = 2,LA
          SEDBALL(L,K) = 1000.0*HBED(L,K)*(BDENBED(L,K) - 1000.0*PORBED(L,K))
        enddo
      enddo

      ! *** SPLIT TOTAL SEDIMENT INTO CLASSES
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do L = 2,LA
              SEDB(L,K,NS) = SEDB(L,K,NS)*SEDBALL(L,K)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          do K = 1,KB
            do L = 2,LA
              SNDB(L,K,NX) = SNDB(L,K,NX)*SEDBALL(L,K)
            enddo
          enddo
        enddo
      endif
    endif     ! *** END OF SEDIMENT INITIALIZATION: MASS FRACTION
    
    ! *** SET TOP LAYER
    if( LSEDZLJ )then
      do L = 2,LA
        KBT(L) = KBB
        do K = KBB,1,-1
          if( HBED(L,K) <= 0.0 )then
            KBT(L) = min(KB,K + 1)
            exit
          endif
        enddo
      enddo
    else
      do L = 2,LA
        KBT(L) = 1
        do K = 1,KB
          if( HBED(L,K) <= 0.0 )then
            KBT(L) = max(1,K - 1)
            exit
          endif
        enddo
      enddo
    endif
  endif       ! *** END OF COLD START ISEDINT > 1
  ! *** END OF OPTIONAL HOT/COLD START LOGIC
  
  ! *** SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
  if( ISTRAN(6) >= 1 .and. .not. LSEDZLJ )then
    if( IWRSP(1) == 0 )then
      do K = 1,KB
        do L = 2,LA
          TAURS(L,K) = TAUR(1)
          WRSPS(L,K) = WRSPO(1)
        enddo
      enddo
    endif
    if( IWRSPB(1) == 0 )then
      do K = 1,KB
        do L = 2,LA
          TAURB(L,K) = 1.E6
          WRSPB(L,K) = 0.0
        enddo
      enddo
    endif
    if( IWRSP(1) >= 1 .and. IWRSP(1) < 99 )then
      do K = 1,KB
        do L = 2,LA
          TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1), VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
          WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1))
        enddo
      enddo
    endif

    if( IWRSPB(1) >= 1 )then
      do K = 1,KB
        do L = 2,LA
          TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))
          WRSPB(L,K) = CSEDRESB
        enddo
      enddo
    endif
  endif

  ! *** SET SEDIMENT VOLUME FRACTIONS
  do K = 1,KB
    do L = 2,LA
      BEDLINIT(L,K) = 0.0
      BEDDINIT(L,K) = 0.0
    enddo
  enddo
  if( ISTRAN(6) >= 1 )then
    do NS = 1,NSED
      do K = 1,KB
        do L = 2,LA
          VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)
          VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)
        enddo
      enddo
    enddo
  endif
  if( ISTRAN(7) >= 1 )then
    do NX = 1,NSND
      NS = NSED + NX
      do K = 1,KB
        do L = 2,LA
          VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)
          VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)
        enddo
      enddo
    enddo
  endif
  if( ISTRAN(6) >= 1 )then
    do NS = 1,NSED
      do K = 1,KB
        do L = 2,LA
          BEDLINIT(L,K) = BEDLINIT(L,K) + VFRBED(L,K,NS)
          BEDDINIT(L,K) = BEDDINIT(L,K) + VFRBED1(L,K,NS)
        enddo
      enddo
    enddo
  endif
  if( ISTRAN(7) >= 1 )then
    do NX = 1,NSND
      NS = NSED + NX
      do K = 1,KB
        do L = 2,LA
          BEDLINIT(L,K) = BEDLINIT(L,K) + VFRBED(L,K,NS)
          BEDDINIT(L,K) = BEDDINIT(L,K) + VFRBED1(L,K,NS)
        enddo
      enddo
    enddo
  endif
  if( ISTRAN(6) >= 1 )then
    do NS = 1,NSED
      do K = 1,KB
        do L = 2,LA
          if( BEDLINIT(L,K) > 0.0 )then
            VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
          else
            VFRBED(L,K,NS) = 0.0
            BEDLINIT(L,K)  = 0.0
          endif
          if( BEDDINIT(L,K) > 0.0 )then
            VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
          else
            VFRBED1(L,K,NS) = 0.0
            BEDDINIT(L,K)  = 0.0
          endif
        enddo
      enddo
    enddo
  endif
  if( ISTRAN(7) >= 1 )then
    do NX = 1,NSND
      NS = NSED + NX
      do K = 1,KB
        do L = 2,LA
          if( BEDLINIT(L,K) > 0.0 )then
            VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
          else
            VFRBED(L,K,NS) = 0.0
            BEDLINIT(L,K)  = 0.0
          endif
          if( BEDDINIT(L,K) > 0.0 )then
            VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
          else
            VFRBED1(L,K,NS) = 0.0
            BEDDINIT(L,K)  = 0.0
          endif
        enddo
      enddo
    enddo
  endif

  ! *** INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA
  do K = 1,KB
    do L = 2,LA
      SEDBT(L,K) = 0.0
      SNDBT(L,K) = 0.0
    enddo
  enddo
  if( ISTRAN(6) >= 1 )then
    do NS = 1,NSED
      do K = 1,KB
        do L = 2,LA
          SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
        enddo
      enddo
    enddo
  endif
  if( ISTRAN(7) >= 1 )then
    do NX = 1,NSND
      do K = 1,KB
        do L = 2,LA
          SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NX)
        enddo
      enddo
    enddo
  endif

  if( LSEDZLJ )then
    do L = 2,LA
      do K = 1,KB
        SEDDIA50(L,K) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))
      enddo
    enddo
  endif

  ! ***********************************************************************************
  ! *** INITIALIZE BED BOTTOM ELEVATION
  if( NSEDFLUME /= 1 .and. NSEDFLUME /= 2 )then
    ! *** INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA
    SEDBT = 0.0
    SNDBT = 0.0
    if( ISTRAN(6) >= 1 )then
      do NS = 1,NSED
        do K = 1,KB
          do L = 2,LA
            if( SEDB(L,K,NS) > 0.0 )then
              SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
            else
              SEDB(L,K,NS) = 0.0
            endif
          enddo
        enddo
      enddo
    endif
    if( ISTRAN(7) >= 1 )then
      do NX = 1,NSND
        do K = 1,KB
          do L = 2,LA
            if( SNDB(L,K,NX) > 0.0 )then
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NX)
            else
              SNDB(L,K,NX) = 0.0
            endif
          enddo
        enddo
      enddo
    endif

    ! *** GET TOTAL THICKNESSES AND BEDROCK ELEVATION
    do L = 2,LA
      HBEDA(L) = 0.0
    enddo
    do L = 2,LA
      if( LBED(L) ) CYCLE            ! *** HARDBOTTOM BYPASS
      do K = 1,KB
        if( SEDBT(L,K) <= 1E-6 )then
          SEDB(L,K,1:NSED) = 0.0
          SEDBT(L,K) = 0.0
        endif
        if( SNDBT(L,K) <= 1E-6 )then
          SNDB(L,K,1:NSND) = 0.0
          SNDBT(L,K) = 0.0
        endif
        if( (SEDBT(L,K) + SNDBT(L,K)) <= 1E-6 )then
          HBED(L,K) = 0.0
          SEDB(L,K,1:NSED) = 0.0
          SNDB(L,K,1:NSND) = 0.0
        endif
        HBEDA(L) = HBEDA(L) + HBED(L,K)

        ! *** QC
        if( HBED(L,K) > 0. .and. PORBED(L,K) > 0. )then
          TTHICK = 0.
          do NS = 1,NSED
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SEDB(L,K,NS)*DSEDGMM
          enddo
          do NX = 1,NSND
            NS = NSED + NX
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SNDB(L,K,NX)*DSEDGMM
          enddo
          TTHICK = TTHICK/(1.0 - PORBED(L,K))
          if( ABS(HBED(L,K) - TTHICK) > 1E-4 )then
            PRINT '(" *** WARNING - BED LAYER THICKNESS IS NOT CONSISTENT.  ADJUSTING THICKNESS (L, K, HOLD, HNEW): ",2I5,2F12.5)', Map2Global(L).LG,K,HBED(L,K),TTHICK
            write(mpi_efdc_out_unit,'(" *** WARNING - BED LAYER THICKNESS IS NOT CONSISTENT.  ADJUSTING THICKNESS (L, K, HOLD, HNEW): ",2I5,2F12.5)') Map2Global(L).LG,K,HBED(L,K),TTHICK
            HBED(L,K) = TTHICK
          endif

        else
          HBED(L,K)   = 0.0
          PORBED(L,K) = BEDPORC
          VDRBED(L,K) = PORBED(L,K)/(1.0 - PORBED(L,K))
          SEDB(L,K,:) = 0.
          SNDB(L,K,:) = 0.
        endif
      enddo

      ! *** Check for empty layers and fill as needed
      if( NSEDFLUME == 3 .and. HBED(L,KB) > 0.0 )then
        do K = KB - 1,1,-1
          if( HBED(L,K) == 0.0 )then
            PORBED(L,K) = PORBED(L,K+1)
            VDRBED(L,K) = VDRBED(L,K+1)
          endif
        enddo
      endif
    enddo

    do L = 2,LA
      ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** Hard bottom elevation
    enddo

    ! *** SET TOP LAYER
    if( LSEDZLJ )then
      do L = 2,LA
        KBT(L) = KBB
        do K = KBB,1,-1
          if( HBED(L,K) <= 0.0 )then
            KBT(L) = min(KB,K + 1)
            exit
          endif
        enddo
      enddo
    else
      do L = 2,LA
        KBT(L) = 1
        do K = 1,KB
          if( HBED(L,K) <= 0.0 )then
            KBT(L) = max(1,K - 1)
            exit
          endif
        enddo
      enddo
    endif

  endif

  ! ***********************************************************************************
  ! *** IF N = 1 AND ISTRAN(5) = 1 CHECK INITIAL TOXIC CONCENTRATIONS IN
  ! *** BED AND REINITILIZE IF NECESSARY
  if( ISTRAN(5) >= 1 )then
    ! *** GET EACH SOLIDS CLASS FOR EACH TOXIC: NSP2
    do NT = 1,NTOX
      NSP2(NT) = NSED + NSND                        ! *** Kd  Approach ISTOC(NT) = 0
      if( ISTOC(NT) == 1 ) NSP2(NT) = NSP2(NT) + 2  ! *** DOC AND POC (POC NON-SEDIMENT RELATED)  (3 Phase)
      if( ISTOC(NT) == 2 ) NSP2(NT) = NSP2(NT) + 1  ! *** DOC AND POC FRACTIONALLY DISTRIBUTED    (3 Phase)
      ! *** ISTOC(NT) = 0 and ISTOC(NT) = 3               POC fOC*SED/SND BASED ONLY              (2 Phase)
    enddo

    if( ISRESTI == 0 .or. ISCI(5) == 0 )then

      ! *** CALCULATE TOTAL SEDIMENT IN THE BED
      do K = 1,KB
        do L = 1,LC
          SEDBALL(L,K) = SEDBT(L,K) + SNDBT(L,K)
        enddo
      enddo

      ! *** *******************************************************************C
      !
      ! *** CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
      !
      ! *** TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED  (DIMENSIONLESS)
      ! *** TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED  (DIMENSIONLESS)
      call CALTOXB_FRACTIONS(1,LASED)

      ! *** COMPUTE FINAL FRACTION (DIMENSIONLESS) ADSORBED ON ALL SOLIDS AND DOC (I.E. NOT DISSOLVED)
      do NT = 1,NTOX
        do K = 1,KB
          do L = 2,LA
            if( SEDBALL(L,K) > 0.0 )then
              ! ***                M           (   POREWATER (M)        SOLIDS (M)   )
              TOXPFTB(L,K,NT) = TOXPFTB(L,K,NT)/(PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT))
            else
              TOXPFTB(L,K,NT) = 1.0
            endif
          enddo
        enddo
      enddo

      ! *** CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC
      ! *** CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG
      ! *** TO TOXB UNITS OF OF MG/M**2
      do NT = 1,NTOX
        if( ITXBDUT(NT) == 0 )then
          ! *** INPUT UNITS OF UG/L ( = MG/M3)
          do K = 1,KB
            do L = 2,LA
              ! *** MG/M2   =   M          MG/M3
              TOXB(L,K,NT)  = HBED(L,K)*TOXB(L,K,NT)
              TOXB1(L,K,NT) = TOXB(L,K,NT)
            enddo
          enddo
        endif
        if( ITXBDUT(NT) == 1 )then
          ! *** PMC - BEFORE 2017-06 EFDC ASSUMED ONLY PARTICLE MASS IN FILE.  TYPICAL LABRATORY ANALYSIS ARE TOTAL MASS TOXIC/MASS OF SEDIMENT.
          ! ***      ITXBDUT(NT) = 1 NOW REPRESENTS THE INPUT TOTAL TOXB AS MG/KG
          ! *** INPUT UNITS OF MG/KG
          do K = 1,KB
            do L = 2,LA
              !TOXB(L,K,NT)  = 0.001*TOXB(L,K,NT)*(SEDBT(L,K) + SNDBT(L,K))/TOXPFTB(L,K,NT)   DEPRECATED
              ! *** MG/M2   =  KG/G    MG/KG     (          G/M2         )
              TOXB(L,K,NT)  = 0.001*TOXB(L,K,NT)*(SEDBT(L,K) + SNDBT(L,K))
              TOXB1(L,K,NT) = TOXB(L,K,NT)
            enddo
          enddo
        endif
      enddo

      ! *** DIAGNOSTICS OF INITIALIZATION
      call writebreak(mpi_efdc_out_unit)
      NS = 0
      write(mpi_efdc_out_unit,*) ' Sediment Bed Check: TOXBED.DIA'
      do NT = 1,NTOX
        write(mpi_efdc_out_unit,*) 'TOXIC BED FOR CLASS: ',NT
        do LP = 1,LASED
          L = LSED(LP)
          do K = 1,KB
            if( HBED(L,K) > 0.0 )then
              TMP1 = TOXB(L,K,NT)/HBED(L,K)   ! *** MG/M^3
            else
              if( TOXB(L,K,NT) > 0.0 )then
                write(mpi_efdc_out_unit,*) 'WARNING - TOXIC IC HAS TOXB > 0 FOR ZERO THICKNESS LAYER: ', Map2Global(L).LG, K, NT, TOXB(L,K,NT)
                TOXB(L,K,NT) = 0.0
                NS = 1
              endif
              TMP1 = 0.0
            endif
            write(mpi_efdc_out_unit,2222)  Map2Global(L).IG, Map2Global(L).JG, K, TOXPFTB(L,K,NT), TOXB(L,K,NT), TMP1, TOX(L,KSZ(L),NT)
          enddo
        enddo
      enddo
      if( NS == 1 ) PRINT *, 'BEDINIT:  Found one or more cells with bad TOX/SED settings.  Check mpi_efdc_out_file for details. '        
    endif
  endif

2222 FORMAT(3I5,7E13.4)

  ! *** INITIALIZE FRACTION OF PARTICULATE ORGANIC CARBON IN BED BY FUNCTION
  if( ISTPOCB == 4 )then
    IVAL = 0
    do NT = 1,NTOX
      if( ISTOC(NT) >= 2)IVAL = 1
    enddo
    if( IVAL == 1 )CALL SETFPOCB(0)
  endif

  ! *** CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS
  if( .not. LSEDZLJ )then
    do K = 1,KB
      do L = 2,LA
        if( K <= KBT(L) )then
          FVOLSSD = 1.0/(1.0 + VDRBED(L,K))
          FVOLSED = 0.0
          do NS = 1,NSED
            FVOLSED = FVOLSED + VFRBED(L,K,NS)
          enddo
          FVOLSND = 0.0
          do NX = 1,NSND
            NS = NSED + NX
            FVOLSND = FVOLSND + VFRBED(L,K,NS)
          enddo
          FVOLSED = FVOLSSD*FVOLSED
          FVOLSND = FVOLSSD*FVOLSND
          VDRBEDSND(L,K) = SNDVDRD
          if( FVOLSED > 1.0E-18 )then
            VDRBEDSED(L,K) = ((FVOLSED + FVOLSND)*VDRBED(L,K)-  FVOLSND*SNDVDRD)/FVOLSED
          else
            VDRBEDSED(L,K) = 0.0
          endif
        else
          VDRBEDSND(L,K) = 0.0
          VDRBEDSED(L,K) = 0.0
        endif
      enddo
    enddo

    ! *** ADD ACTIVE ARMORING LAYER IF NOT PRESENT IN INITIAL OR RESTART
    if( IHOTSTRT == 0 .and. ISNDAL == 2 .and. IALSTUP > 0 .and. KB > 1 )then

      do L = 2,LA
        KTOPTP = KBT(L)       ! *** NEW PARENT LAYER
        KTOPP1 = KBT(L) + 1   ! *** NEW ACTIVE LAYER

        TTHICK = 0.
        do K = 1,KB
          TTHICK = TTHICK + HBED(L,K)
        enddo

        ! *** MAKE SURE THERE IS SEDIMENT IN THE CELL
        if( TTHICK > 0. )then
          if( KTOPTP == KB )then
            ! *** CASE KBT = KB

            if( HBED(L,KTOPTP) < HBEDAL )then
              ! *** ACTIVE LAYER IS NOT THICK ENOUGH. REPARTITIION KBT AND KBT-1

              KT1 = KBT(L)       ! *** ACTIVE LAYER
              KT2 = KBT(L) - 1   ! *** PARENT LAYER
              if( (HBED(L,KT1) + HBED(L,KT2)) < HBEDAL )then
                write(STR,'(A,I5)')'BEDINIT: ISNDAL = 2 AND IALSTUP>0: BED LAYER IS NOT THICK ENOUGH.  L = ', Map2Global(L).LG
                call STOPP(STR)
              endif
              HBEDP = HBEDAL

              FRACT1  = HBEDP/HBED(L,KT1)
              FRACT2  = (HBED(L,KT2) - (HBEDP - HBED(L,KT1)))/HBED(L,KT2)
              FRACT2C = 1. - FRACT2

              ! *** INCREASE TOP LAYER TO HBEDP
              HBED(L,KT1 ) = HBEDP
              SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)
              SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)
              do NS = 1,NSED
                SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)
                SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)
              enddo
              do NX = 1,NSND
                SNDB(L,KT1,NX)  = FRACT1*SNDB(L,KT1,NX)
                SNDB1(L,KT1,NX) = FRACT1*SNDB1(L,KT1,NX)
              enddo
              do NT = 1,NTOX
                TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)
                TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)
              enddo

              ! *** REDUCE LAYER BELOW
              HBED(L,KT2)  = FRACT2*HBED(L,KT2)
              SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)
              SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)
              do NS = 1,NSED
                SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)
                SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)
              enddo
              do NX = 1,NSND
                SNDB(L,KT2,NX)  = FRACT2*SNDB(L,KT2,NX)
                SNDB1(L,KT2,NX) = FRACT2*SNDB1(L,KT2,NX)
              enddo
              do NT = 1,NTOX
                TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)
                TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)
              enddo
            else
              ! *** ACTIVE LAYER IS TOO THICK. REPARTITIION KBT AND KBT-1
              KT1 = KBT(L)     ! *** ACTIVE LAYER
              KT2 = KBT(L) - 1   ! *** PARENT LAYER

              HBEDP = HBEDAL

              FRACT1  = HBEDP/HBED(L,KT1)
              FRACT2  = (HBED(L,KT2) + (HBED(L,KT1) - HBEDP))/HBED(L,KT2)
              FRACT2C = 1. - FRACT2

              ! *** DECREASE TOP LAYER TO HBEDP
              HBED(L,KT1)  = HBEDP
              SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)
              SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)
              do NS = 1,NSED
                SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)
                SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)
              enddo
              do NX = 1,NSND
                SNDB(L,KT1,NX)  = FRACT1*SNDB(L,KT1,NX)
                SNDB1(L,KT1,NX) = FRACT1*SNDB1(L,KT1,NX)
              enddo
              do NT = 1,NTOX
                TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)
                TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)
              enddo

              ! *** INCREASE THE LAYER BELOW
              HBED(L,KT2)  = FRACT2*HBED(L,KT2)
              SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)
              SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)
              do NS = 1,NSED
                SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)
                SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)
              enddo
              do NX = 1,NSND
                SNDB(L,KT2,NX)  = FRACT2*SNDB(L,KT2,NX)
                SNDB1(L,KT2,NX) = FRACT2*SNDB1(L,KT2,NX)
              enddo
              do NT = 1,NTOX
                TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)
                TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)
              enddo

            endif
          else

            ! *** STANDARD CASE OF KBT < KB
            if( HBED(L,KTOPTP) < HBEDAL )then
              ! *** PARENT LAYER IS NOT THICK ENOUGH. REPARTITIION KBT AND KBT-1

              KT1 = KBT(L)     ! *** PARENT LAYER
              KT2 = KBT(L) - 1   ! *** LAYER BELOW PARENT
              if( KT2 > 0 )then
                HBEDP = min(HBEDAL*2.,HBED(L,KT2))

                FRACT1  = HBEDP/HBED(L,KT1)
                FRACT2  = (HBED(L,KT2) - (HBEDP - HBED(L,KT1)))/HBED(L,KT2)
                FRACT2C = 1. - FRACT2

                ! *** INCREASE TOP LAYER TO HBEDP
                HBED(L,KT1) = HBEDP
                if( ISTRAN(6) > 0 )then
                  SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)
                  do NS = 1,NSED
                    SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)
                    SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)
                  enddo
                endif
                if( ISTRAN(7) > 0 )then
                  SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)
                  do NX = 1,NSND
                    SNDB(L,KT1,NX)  = FRACT1*SNDB(L,KT1,NX)
                    SNDB1(L,KT1,NX) = FRACT1*SNDB1(L,KT1,NX)
                  enddo
                endif
                if( ISTRAN(5) > 0 )then
                  do NT = 1,NTOX
                    TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)
                    TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)
                  enddo
                endif

                ! *** REDUCE LAYER BELOW
                HBED(L,KT2) = FRACT2*HBED(L,KT2)
                if( ISTRAN(6) > 0 )then
                  SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)
                  do NS = 1,NSED
                    SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)
                    SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)
                  enddo
                endif
                if( ISTRAN(7) > 0 )then
                  SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)
                  do NX = 1,NSND
                    SNDB(L,KT2,NX)  = FRACT2*SNDB(L,KT2,NX)
                    SNDB1(L,KT2,NX) = FRACT2*SNDB1(L,KT2,NX)
                  enddo
                endif
                if( ISTRAN(5) > 0 )then
                  do NT = 1,NTOX
                    TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)
                    TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)
                  enddo
                endif
              endif
            endif

            ! *** SPLIT THE PARENT LAYER INTO THE PARENT (KTOPTP) AND THE ACTIVE (KTOPP1)
            FRACACT(L) = HBEDAL/HBED(L,KTOPTP)
            FRACPAR(L) = (HBED(L,KTOPTP) - HBEDAL)/HBED(L,KTOPTP)
            HBED(L,KTOPP1) = FRACACT(L)*HBED(L,KTOPTP)  ! ACTIVE
            HBED(L,KTOPTP) = FRACPAR(L)*HBED(L,KTOPTP)  ! PARENT
            PORBED(L,KTOPP1) = PORBED(L,KTOPTP)
            PORBED1(L,KTOPP1) = PORBED1(L,KTOPTP)
            VDRBED(L,KTOPP1) = VDRBED(L,KTOPTP)
            VDRBED1(L,KTOPP1) = VDRBED1(L,KTOPTP)
            BDENBED(L,KTOPP1) = BDENBED(L,KTOPTP)
            BDENBED1(L,KTOPP1) = BDENBED1(L,KTOPTP)
            SEDBT(L,KTOPP1) = FRACACT(L)*SEDBT(L,KTOPTP)
            SEDBT(L,KTOPTP) = FRACPAR(L)*SEDBT(L,KTOPTP)
            SNDBT(L,KTOPP1) = FRACACT(L)*SNDBT(L,KTOPTP)
            SNDBT(L,KTOPTP) = FRACPAR(L)*SNDBT(L,KTOPTP)
            STDOCB(L,KTOPP1) = STDOCB(L,KTOPTP)
            STPOCB(L,KTOPP1) = STPOCB(L,KTOPTP)

            if( ISTRAN(6) > 0 )then
              do NS = 1,NSED
                SEDB(L,KTOPP1,NS)    = FRACACT(L)*SEDB(L,KTOPTP,NS)
                SEDB1(L,KTOPP1,NS)   = FRACACT(L)*SEDB1(L,KTOPTP,NS)
                SEDB(L,KTOPTP,NS)    = FRACPAR(L)*SEDB(L,KTOPTP,NS)
                SEDB1(L,KTOPTP,NS)   = FRACPAR(L)*SEDB1(L,KTOPTP,NS)
                STFPOCB(L,KTOPP1,NS) = STFPOCB(L,KTOPTP,NS)
                VFRBED(L,KTOPP1,NS)  = SDEN(NS)*SEDB(L,KTOPP1,NS)
                VFRBED1(L,KTOPP1,NS) = SDEN(NS)*SEDB1(L,KTOPP1,NS)
                VFRBED(L,KTOPTP,NS)  = SDEN(NS)*SEDB(L,KTOPTP,NS)
                VFRBED1(L,KTOPTP,NS) = SDEN(NS)*SEDB1(L,KTOPTP,NS)
              enddo
            endif
            if( ISTRAN(7) > 0 )then
              do NX = 1,NSND
                NS = NSED + NX
                SNDB(L,KTOPP1,NX)    = FRACACT(L)*SNDB(L,KTOPTP,NX)
                SNDB1(L,KTOPP1,NX)   = FRACACT(L)*SNDB1(L,KTOPTP,NX)
                SNDB(L,KTOPTP,NX)    = FRACPAR(L)*SNDB(L,KTOPTP,NX)
                SNDB1(L,KTOPTP,NX)   = FRACPAR(L)*SNDB1(L,KTOPTP,NX)
                STFPOCB(L,KTOPP1,NS) = STFPOCB(L,KTOPTP,NS)
                VFRBED(L,KTOPP1,NS)  = SDEN(NS)*SNDB(L,KTOPP1,NX)
                VFRBED1(L,KTOPP1,NS) = SDEN(NS)*SNDB1(L,KTOPP1,NX)
                VFRBED(L,KTOPTP,NS)  = SDEN(NS)*SNDB(L,KTOPTP,NX)
                VFRBED1(L,KTOPTP,NS) = SDEN(NS)*SNDB1(L,KTOPTP,NX)
              enddo
            endif
            
            ! *** Compute fractions from totals
            if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
              SUMSED = 0.0
              do NS = 1,NSEDS
                do K = 1,KB
                  SUMSED(K)  = SUMSED(K) + VFRBED(L,K,NS)
                enddo
              enddo
              do NS = 1,NSEDS
                do K = 1,KB
                  if( SUMSED(K) > 0.0 )then
                    VFRBED(L,K,NS)  = VFRBED(L,K,NS)/SUMSED(K)
                  else
                    VFRBED(L,K,NS) = 0.0
                  endif
                enddo
              enddo            
              SUMSED = 0.0
              do NS = 1,NSEDS
                do K = 1,KB
                  SUMSED(K)  = SUMSED(K) + VFRBED1(L,K,NS)
                enddo
              enddo
              do NS = 1,NSEDS
                do K = 1,KB
                  if( SUMSED(K) > 0.0 )then
                    VFRBED1(L,K,NS)  = VFRBED1(L,K,NS)/SUMSED(K)
                  else
                    VFRBED1(L,K,NS) = 0.0
                  endif
                enddo
              enddo            
            endif
            
            if( ISTRAN(5) > 0 )then
              do NT = 1,NTOX
                TOXB(L,KTOPP1,NT) = FRACACT(L)*TOXB(L,KTOPTP,NT)
                TOXB1(L,KTOPP1,NT) = FRACACT(L)*TOXB1(L,KTOPTP,NT)
                TOXB(L,KTOPTP,NT) = FRACPAR(L)*TOXB(L,KTOPTP,NT)
                TOXB1(L,KTOPTP,NT) = FRACPAR(L)*TOXB1(L,KTOPTP,NT)
                TOXPFTB(L,KTOPP1,NT) = TOXPFTB(L,KTOPTP,NT)
              enddo
              do NT = 1,NTOX
                do NS = 1,NSED + NSND + 2
                  TOXPFB(L,KTOPP1,NS,NT) = TOXPFB(L,KTOPTP,NS,NT)
                enddo
              enddo
            endif

            ! *** UPDATE TOP LAYER
            KBT(L) = KBT(L) + 1

          endif  ! *** END OF KBT<KB CHECK
        endif    ! *** END OF TOTAL THICKNESS CHECK: TTHICK>0
      enddo  ! *** END OF LA LOOP

    endif    ! *** END OF REPARTITIONING ACTIVE ARMORING LAYER AT STARTUP IF REQUESTED

  endif   ! *** END OF LSEDZLJ BYPASS

  if( .not. LSEDZLJ )then
    ! *** *******************************************************************C
    ! *** SET MEAN D50 AND D90 (meters)
    if( ISTRAN(7) >= 1 )then
      do K = 1,KB
        do L = 2,LA
          SEDDIA50(L,K) = 0.
          SEDDIA90(L,K) = 0.
          SNDBT(L,K) = 0.
        enddo
      enddo
      do NX = 1,NSND
        NS = NSED + NX
        do K = 1,KB
          do L = 2,LA
            SEDDIA50(L,K) = SEDDIA50(L,K) + SNDB(L,K,NX)*LOG(SEDDIA(NS))   ! *** SNDDIA is in meters
            SNDBT(L,K)    = SNDBT(L,K) + SNDB(L,K,NX)
          enddo
        enddo
      enddo
      do K = 1,KB
        do L = 2,LA
          if( SNDBT(L,K) > 0. )then
            SEDDIA50(L,K) = SEDDIA50(L,K)/SNDBT(L,K)
          endif
        enddo
      enddo

      do K = 1,KB
        do L = 2,LA
          if( SNDBT(L,K) > 0.1 )then  ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
            SEDDIA50(L,K) = EXP(SEDDIA50(L,K))
          else
            SEDDIA50(L,K) = SEDDIA(NSED + 1)
          endif
        enddo
      enddo

    endif  ! *** END OF NON-COHESIVE GRAIN SIZE CALCS

  elseif( NSEDFLUME == 3 .and. IHTSTRT == 0 )then
    ! *** CLEAN UP SPATIALLY VARYING INPUTS TO WORK WITH SEDZLJ
    do LP = 1,LASED
      L = LSED(LP)
      do K = KB,1,-1
        if( HBED(L,K) > 0.0 )then
          BULKDENS(K,L) = SEDBT(L,K)/HBED(L,K)*1.E-6   ! *** Convert g/m^3 to g/cm^3
          LAYERACTIVE(K,L) = 2                         ! *** Flag original in-place sediment layers
        else
          if( K < KB )then
            TMP1 = BULKDENS(K + 1,L)
          else
            TMP1 = 2.65*(1. - BEDPORC)
          endif
          BULKDENS(K,L) = TMP1
          LAYERACTIVE(K,L) = 0
        endif

        do NS = 1,NSEDS
          if( SEDBT(L,K) > 0.0 )then
            PERSED(NS,K,L) = SEDB(L,K,NS)/SEDBT(L,K)
          else
            PERSED(NS,K,L) = 0.0
          endif
        enddo
        DTOTAL = SUM(PERSED(:,K,L))
        if( DTOTAL > 0. )then
          PERSED(1:NSEDS,K,L) = PERSED(1:NSEDS,K,L)/DTOTAL   ! *** Ensure precision for mass balance
        else
          PERSED(1:NSEDS,K,L) = 0.0
        endif
        TSED(K,L)  = SEDBT(L,K)*1.E-4                      ! *** SEDBT is in units of g/m^2 and TSED in units of g/cm^2.
        TSED0(K,L) = TSED(K,L)
      enddo
    enddo

    ! *** Initialize KBT to the first layer with mass
    do LP = 1,LASED
      L = LSED(LP)
      KBT(L) = -1
      do K = 1,KB
        if( TSED0(K,L) > 0. )then
          KBT(L) = K
          exit
        endif
      enddo
      if( KBT(L) == -1 ) KBT(L) = KB
    enddo

  endif    ! *** END OF LSEDZLJ BYPASS

  ! ***********************************************************************************************************
1000 continue     ! JUMP FROM THE HOTSTART (SEDZLJ INITIALIZATION CONTINUES FROM THE HOTSTART SECTION)

  do K = 1,KB
    do L = 2,LA
      if( (K <= KBT(L) .and. .not. LSEDZLJ) .or. (K >= KBT(L) .and. LSEDZLJ) )then
        SDENAVG(L,K) = (BDENBED(L,K) - 1000.0*PORBED(L,K))/(1.0 - PORBED(L,K))
      else
        SDENAVG(L,K) = (BDENBED(L,KBT(L)) - 1000.0*PORBED(L,KBT(L)))/(1.0 - PORBED(L,KBT(L)))
      endif
      if( SDENAVG(L,K) <= 0. ) SDENAVG(L,K) = 0.0
    enddo
  enddo

  if(master_id == process_id )then
    ! ***********************************************************************************************************
    ! *** WRITE DIAGNOSTIC FILES FOR BED INITIALIZATION
    write(6,'(A)') 'WRITE DIAGNOSTIC FILES FOR BED INITIALIZATION INTO mpi_efdc_out_file'
  endif
  
  call WriteBreak(mpi_efdc_out_unit)
  write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.ELV'
  write(mpi_efdc_out_unit,111)
  do L = 2,LA
    SURF = HP(L) + BELV(L)
    write(mpi_efdc_out_unit,101) Map2Global(L).IG, Map2Global(L).JG, ZELBEDA(L), HBEDA(L), BELV(L), HP(L), SURF
  enddo

  call WriteBreak(mpi_efdc_out_unit)
  write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.ZHB'
  write(mpi_efdc_out_unit,112) (K,K = 1,KB)
  do L = 2,LA
    write(mpi_efdc_out_unit,101) Map2Global(L).IG,  Map2Global(L).JG,  ZELBEDA(L), HBEDA(L),(HBED(L,K),K = 1,KB)
  enddo

  if( ISTRAN(6) > 0 )then
    call WriteBreak(mpi_efdc_out_unit)
    write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.SED'
    write(mpi_efdc_out_unit,113) (K,K = 1,KB)
    do L = 2,LA
      write(mpi_efdc_out_unit,101)  Map2Global(L).IG, Map2Global(L).JG, (SEDBT(L,K),K = 1,KB), SUM(SEDBT(L,:))
    enddo
  endif
    
  if( ISTRAN(7) > 0 )then
    call WriteBreak(mpi_efdc_out_unit)
    write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.SND'
    write(mpi_efdc_out_unit,114) (K,K = 1,KB)
    do L = 2,LA
      write(mpi_efdc_out_unit,101)  Map2Global(L).IG, Map2Global(L).JG, (SNDBT(L,K),K = 1,KB), SUM(SNDBT(L,:))
    enddo
  endif
    
  if( ISTRAN(5) > 0 )then
    call WriteBreak(mpi_efdc_out_unit)
    write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.TOX'
    do NT = 1,NTOX
      write(mpi_efdc_out_unit,115) NT, (K,K = 1,KB)
      do L = 2,LA
        write(mpi_efdc_out_unit,101) Map2Global(L).IG, Map2Global(L).JG, (TOXB(L,K,NT),K = 1,KB)
      enddo
    enddo
  endif
    
  call WriteBreak(mpi_efdc_out_unit)
  write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.VDR'
  write(mpi_efdc_out_unit,116) (K,K = 1,KB)
  do L = 2,LA
    write(mpi_efdc_out_unit,101) Map2Global(L).IG, Map2Global(L).JG, (VDRBED(L,K),K = 1,KB)
  enddo

  call WriteBreak(mpi_efdc_out_UNIT)
  write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.POR'
  write(mpi_efdc_out_unit,117) (K,K = 1,KB)
  do L = 2,LA
    write(mpi_efdc_out_unit,101) Map2Global(L).IG, Map2Global(L).JG, (PORBED(L,K),K = 1,KB)
  enddo

  call WriteBreak(mpi_efdc_out_unit)
  write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.BDN'
  write(mpi_efdc_out_unit,118) (K,K = 1,KB)
  do L = 2,LA
    write(mpi_efdc_out_unit,101) Map2Global(L).IG, Map2Global(L).JG, (BDENBED(L,K),K = 1,KB)
  enddo

  call WriteBreak(mpi_efdc_out_unit)
  write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINIT.VR'
  write(mpi_efdc_out_unit,119) (K,K = 1,KB)
  do L = 2,LA
    write(mpi_efdc_out_unit,191) Map2Global(L).IG, Map2Global(L).JG, (VDRBED(L,K),K = 1,KB)
  enddo

  if( ISTRAN(6) > 0 )then
    call WriteBreak(mpi_efdc_out_unit)
    write(mpi_efdc_out_unit,*) ' Sediment Bed Check: BEDINITC.VV'
    write(mpi_efdc_out_unit,120)
    do L = 2,LA
      K = KBT(L)
      write(mpi_efdc_out_unit,191,IOSTAT = ISO) Map2Global(L).LG, KBT(L), VFRBED(L,K,1), PORBED(L,K), VDRBED(L,K), VDRBEDSED(L,K)
    enddo
  endif
  ! *** END DIAGNOSTIC FILE DUMP

  ! *** Update previous time step and save startup bed parameters for later
  do K = 1,KB
    do L = 2,LA
      if( LSEDZLJ )then
        if( K < KBT(L) )then
          VDRBED(L,K)  = VDRBED(L,KBT(L))
          PORBED(L,K)  = PORBED(L,KBT(L))
          VDRBED0(L,K) = VDRBED(L,KBT(L))
        endif
      else
        if( K > KBT(L) )then
          VDRBED(L,K)  = VDRBED(L,KBT(L))
          PORBED(L,K)  = PORBED(L,KBT(L))
          VDRBED0(L,K) = VDRBED(L,KBT(L))
        endif
      endif
      HBED1(L,K)    = HBED(L,K)
      VDRBED1(L,K)  = VDRBED(L,K)
      PORBED1(L,K)  = PORBED(L,K)
      BDENBED1(L,K) = BDENBED(L,K)
      VDRBED0(L,K)  = VDRBED(L,K)
    enddo
  enddo

101 FORMAT(2I5,100E13.5)
102 FORMAT(10X,100E13.5)
    
111 FORMAT('   IL   JL        ZBEDB        HBEDT         BELV        HWCOL         SELV')
112 FORMAT(' ZBEDB HBEDT HBED(K = 1,KB)',/,'   IL   JL        ZBEDB        HBEDT',100I13)
113 FORMAT(' SEDBT(K = 1,KB)'  ,/,'   IL   JL',100I13)
114 FORMAT(' SNDBT(K = 1,KB)',  /,'   IL   JL',100I13)
115 FORMAT(' TOXB(K = 1,KB,NT)  NT = ',I5,/,'   IL   JL',100I10)
116 FORMAT(' VRDBED(K = 1,KB)', /,'   IL   JL',100I13)
117 FORMAT(' PORBED(K = 1,KB)', /,'   IL   JL',100I13)
118 FORMAT(' BDENBED(K = 1,KB)',/,'   IL   JL',20X,100I13)
119 FORMAT(' VDRBED(K = 1,KB)', /,'   IL   JL',100I10)
120 FORMAT('    L  KBT    VFRBED    PORBED    VDRBED VDRBEDSED')

191 FORMAT(2I5,100F10.3)
192 FORMAT(10X,100F10.3)

  return

  END

