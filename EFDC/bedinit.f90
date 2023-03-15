! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE BEDINIT
  !
  ! ***  SUBROUTINE BEDINIT INITIALIZES SEDIMENT AND TOXIC VARIABLES
  ! ***  IN SEDIMENT BED FOR HOT AND COLD START CONDITIONS
  !     CHANGE RECORD
  !     ADDED ADDITIONAL DIAGNOSTIC OUTPUT
  !     MOVED TOXIC INITIALIZATIONS FROM SSEDTOX
  !

  USE GLOBAL
  USE Variables_MPI
#ifdef _MPI
  USE MPI
  Use Broadcast_Routines
#endif

  USE INFOMOD,ONLY:SKIPCOM,READSTR

  IMPLICIT NONE

  INTEGER :: K,L,LP,NS,NX,NT,KTOPP1,IVAL,KTOPTP,IHOTSTRT
  INTEGER :: NSKIP,ISSTYPE,LD,ID,JD,LCORE,IDUM,NL,KT1,KT2
  INTEGER :: ISO,NSITES,NCOHSEDS,NCOHSEDL
  CHARACTER :: STR*200

  INTEGER,ALLOCATABLE,DIMENSION(:) :: LSSCOHSED


  REAL :: CSEDTAUS, CSEDRESS, CSEDTAUB, CSEDRESB, TMPCVT
  REAL :: SURF, FVOLSSD, FVOLSED, FVOLSND, TMP1, SXD
  REAL :: RMULADJ, RUMLADJC, FRACT1, FRACT2, FRACT2C, HBEDP, TTHICK
  REAL(RKD) :: DTOTAL
  
  REAL,ALLOCATABLE,DIMENSION(:) :: FRACACT
  REAL,ALLOCATABLE,DIMENSION(:) :: FRACPAR
  REAL,ALLOCATABLE,DIMENSION(:) :: RADJCOHSEDS
  REAL,ALLOCATABLE,DIMENSION(:) :: TAUDSS
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TAURSS
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TAUNSS
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TEXPSS
  REAL,ALLOCATABLE,DIMENSION(:,:) :: DEPBBSS
  REAL,ALLOCATABLE,DIMENSION(:,:) :: WRSPOSS
  REAL,ALLOCATABLE,DIMENSION(:,:) :: SEDBALL

  ! *** MPI
  INTEGER,ALLOCATABLE,DIMENSION(:) :: LSSCOHSED_Global
  REAL,ALLOCATABLE,DIMENSION(:) :: RADJCOHSEDS_Global
  Integer :: l_local
  
  ALLOCATE(LSSCOHSED(LCM))
  ALLOCATE(FRACACT(LCM))
  ALLOCATE(FRACPAR(LCM))
  ALLOCATE(RADJCOHSEDS(LCM))
  ALLOCATE(SEDBALL(LCM,KBM))

  ! *** For MPI
  Allocate(LSSCOHSED_Global(LCM_Global))
  Allocate(RADJCOHSEDS_Global(LCM_Global))
  
  LSSCOHSED=0
  FRACACT=0.0
  FRACPAR=0.0
  RADJCOHSEDS=0.0
  SEDBALL=0.0

  ! *** SEDIMENT PROPERTIES ARE IN MASS PER UNIT AREA
  ! ***
  ! *** SEDB - COHESIVE SEDIMENTS (G/M2)
  ! *** SNDB - NONCOHESIVE SEDIMENTS (G/M2)
  ! *** TOXB - TOXICS (MG/M2)

  ! *** ZERO LOCAL ARRAYS
  FRACACT=0.0
  FRACPAR=0.0

  ! *** SITE SPECIFIC RESUSPENSION INFORMATION BASED ON HAMRICK'S
  ! *** ANALYSIS OF SEDFLUME CORES JANUARY 2008
  IF( NSED > 0 .AND. IWRSP(1) >= 99 .AND.  .NOT. LSEDZLJ )THEN
    if( process_id == master_id )THEN
      WRITE(*,'(A)')'READING SSCOHSEDPROP.INP'
      OPEN(99,FILE='sscohsedprop.inp')
      STR=READSTR(99)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(99,*,IOSTAT=ISO)NCOHSEDS,NCOHSEDL,RMULADJ
    end if

    ALLOCATE(TAUDSS(NCOHSEDS))
    ALLOCATE(WRSPOSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(TAURSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(TAUNSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(TEXPSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(DEPBBSS(NCOHSEDL,NCOHSEDS))

    if( process_id == master_id )THEN
      DO NSITES=1,NCOHSEDS
        READ(99,*,IOSTAT=ISO)IDUM,RUMLADJC,TAUDSS(NSITES)
        READ(99,*,IOSTAT=ISO)IDUM,(WRSPOSS(NL,NSITES),NL=1,NCOHSEDL)
        DO NL=1,NCOHSEDL
          WRSPOSS(NL,NSITES) = RUMLADJC*RMULADJ*WRSPOSS(NL,NSITES)
        ENDDO
        READ(99,*,IOSTAT=ISO)IDUM,(TAURSS(NL,NSITES),NL=1,NCOHSEDL)
        READ(99,*,IOSTAT=ISO)IDUM,(TAUNSS(NL,NSITES),NL=1,NCOHSEDL)
        READ(99,*,IOSTAT=ISO)IDUM,(TEXPSS(NL,NSITES),NL=1,NCOHSEDL)
        READ(99,*,IOSTAT=ISO)IDUM,(DEPBBSS(NL,NSITES),NL=1,NCOHSEDL)
      ENDDO
      CLOSE(99)
      IF( ISO > 0 )THEN
        CALL STOPP('ERROR READING SSCOHSEDPROP.INP FILE')
      ENDIF
    end if
#   ifdef _MPI
    ! ****************************************************************************
    Call Broadcast_Array(TAUDSS, master_id)
    Call Broadcast_Array(WRSPOSS, master_id)
    Call Broadcast_Array(TAURSS, master_id)
    Call Broadcast_Array(TAUNSS, master_id)
    Call Broadcast_Array(TEXPSS, master_id)
    Call Broadcast_Array(DEPBBSS, master_id)
    ! ****************************************************************************
#   endif

  ENDIF

  ! SET BOTTOM LAYER NUMBER: ORIGINAL:KBB=1, SEDZLJ:KBB=KB
  IF( LSEDZLJ )THEN
    KBB = KB
  ELSE
    KBB = 1
  ENDIF

  ! ***  DETERMINE START UP MODE
  IHOTSTRT = 0
  IF( ISRESTI /= 0 )THEN
    IF( ISCI(6) /= 0 .OR. ISCI(7) /= 0 )THEN
      IHOTSTRT = 1
    ENDIF
  ENDIF

  ! ***********************************************************************************
  ! ***  HOT START INITIALIZATION.  SEDZLJ HOTSTART IS HANDLED IN S_SEDIC
  IF( IHOTSTRT /= 0 .AND. .NOT. LSEDZLJ )THEN

    ! ***  SET POROSITY
    DO K=1,KB
      DO L=2,LA
        PORBED(L,K) = VDRBED(L,K)/(1.0+VDRBED(L,K))
        PORBED1(L,K) = VDRBED1(L,K)/(1.0+VDRBED1(L,K))
      ENDDO
    ENDDO

    ! ***  SET BULK DENSITY
    DO K=1,KB
      DO L=2,LA
        SEDBT(L,K)=0.0
        SNDBT(L,K)=0.0
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NS=1,NSND
        DO K=1,KB
          DO L=2,LA
            SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DO K=1,KB
      DO L=2,LA
        IF( HBED(L,K) > 0.0 )THEN
          ! *** COMPUTE TOTAL/WET DENSITY
          BDENBED(L,K) = 1000.0*PORBED(L,K) + 0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)
        ELSE
          BDENBED(L,K) = 0.0
        ENDIF
      ENDDO
    ENDDO
    DO K=1,KB
      DO L=2,LA
        SEDBT(L,K)=0.0
        SNDBT(L,K)=0.0
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K) = SEDBT(L,K) + SEDB1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NS=1,NSND
        DO K=1,KB
          DO L=2,LA
            SNDBT(L,K) = SNDBT(L,K) + SNDB1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    DO K=1,KB
      DO L=2,LA
        IF( HBED1(L,K) > 0.0 )THEN
          BDENBED1(L,K) = 1000.0*PORBED1(L,K) + 0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED1(L,K)
        ELSE
          BDENBED1(L,K) = 0.0
        ENDIF
      ENDDO
    ENDDO

    ! ***  SET TOP BED LAYER
    DO L=2,LA
      KBT(L) = 1
    ENDDO
    DO K=1,KB
      DO L=2,LA
        IF( HBED(L,K) > 0.) KBT(L) = K
      ENDDO
    ENDDO

    ! ***  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
    IF( ISTRAN(6) >= 1 )THEN
      IF( IWRSP(1) == 0 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURS(L,K)=TAUR(1)
            WRSPS(L,K)=WRSPO(1)
          ENDDO
        ENDDO
      ENDIF
      IF( IWRSPB(1) == 0 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURB(L,K)=1.E6
            WRSPB(L,K)=0.0
          ENDDO
        ENDDO
      ENDIF
      IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1))
          ENDDO
        ENDDO
      ENDIF

      ! *** READ IN THE SPATIALLY ASSIGNED COHESIVE CRITICAL BED SHEAR STRESS AND SURFACE EROSION RATE
      IF( IWRSP(1) >= 99 .AND. ISHOUSATONIC == 0 )THEN
        ! *** Read on master
        if(process_id == master_id )then 
            WRITE(*,'(A)') 'READING SSCOHSEDPMAP.INP'
            OPEN(1,FILE='sscohsedpmap.inp')
            OPEN(2,FILE=OUTDIR//'SSCOHSEDPMAP.OUT')
            STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
            READ(1,*)ISSTYPE
            IF( ISSTYPE == 0 )THEN
              DO L=2,LA_Global
                READ(1,*)LD,ID,JD,LSSCOHSED_Global(L)
                RADJCOHSEDS_Global(L) = 1.
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*)LD,ID,JD,LSSCOHSED_Global(L),RADJCOHSEDS_Global(L)
              ENDDO
            ENDIF
        end if ! *** end on master
        
        ! *** send to all processes
        Call Broadcast_Array(LSSCOHSED_Global, master_id)
        Call Broadcast_Array(RADJCOHSEDS_Global, master_id)
        
        ! *** Map to local
        DO L=2,LA_Global
            l_local = map2local(l).ll
            IF( l_local > 1 )THEN
                LCORE = LSSCOHSED_Global(L)
                TAUDS(l_local) = TAUDSS(LCORE)
            end if
        ENDDO
        
        IF( NCOHSEDL == 1 )THEN
          DO K=1,KB
            ! *** Map to local
            DO L=2,LA_Global
                l_local = map2local(l).ll
                IF( l_local > 1 )THEN
                    LCORE=LSSCOHSED_Global(L)
                    TAURS(l_local,K) = TAURSS(1,LCORE)
                    TAUNS(l_local,K) = TAUNSS(1,LCORE)
                    WRSPS(l_local,K) = RADJCOHSEDS_Global(L)*WRSPOSS(1,LCORE)
                    TEXPS(l_local,K) = TEXPSS(1,LCORE)
                end if
            ENDDO
            
          ENDDO
        ELSE
          DO K=1,KB
            ! *** Map to local
            DO L=2,LA_Global
                l_local = map2local(l).ll
                IF( l_local > 1 )THEN
                    LCORE = LSSCOHSED_Global(L)
                    TAURS(l_local,K) = TAURSS(K,LCORE)
                    !TAUNS(L,K) = TAUNSS(K,LCORE)
                    WRSPS(l_local,K) = RADJCOHSEDS_Global(L)*WRSPOSS(K,LCORE)
                    !TEXPS(L,K) = TEXPSS(K,LCORE)
                end if
            ENDDO
          ENDDO
        ENDIF
        
        if(process_id == master_id )then
            DO L=2,LA
              K = KBT(L)
              WRITE(2,*) L,Map2Global(L).IG, Map2Global(L).JG, TAURS(L,K), WRSPS(L,K)
            ENDDO
            CLOSE(1)
            CLOSE(2)
        end if
        
      ENDIF

      IF( IWRSPB(1) >= 1 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))
            WRSPB(L,K) = CSEDRESB
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    ! ***  SET SEDIMENT VOLUME FRACTIONS
    DO K=1,KB
      DO L=2,LA
        BEDLINIT(L,K)=0.0
        BEDDINIT(L,K)=0.0
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)
            VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)
            VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            ! *** BEGIN DSI
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
            ELSE
              VFRBED(L,K,NS)=0.0
              BEDLINIT(L,K) =0.0
            ENDIF
            IF( BEDDINIT(L,K) > 0.0 )THEN
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
            ELSE
              VFRBED1(L,K,NS)=0.0
              BEDDINIT(L,K) =0.0
            ENDIF
            ! *** END DSI
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
            ELSE
              VFRBED(L,K,NS)=0.0
              BEDLINIT(L,K) =0.0
            ENDIF
            IF( BEDDINIT(L,K) > 0.0 )THEN
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
            ELSE
              VFRBED1(L,K,NS)=0.0
              BEDDINIT(L,K) =0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! ***  INITIALIZE BED BOTTOM ELEVATION
    DO L=2,LA
      HBEDA(L)=0.0
    ENDDO
    DO L=2,LA
      DO K=1,KBT(L)
        HBEDA(L) = HBEDA(L)+HBED(L,K)
      END DO
    ENDDO
    DO L=2,LA
      ZELBEDA(L) = BELV(L)-HBEDA(L)
    ENDDO

    ! ***  INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA
    DO K=1,KB
      DO L=2,LA
        SEDBT(L,K)=0.0
        SNDBT(L,K)=0.0
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NS=1,NSND
        DO K=1,KB
          DO L=2,LA
            SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF( LSEDZLJ )THEN
      DO L=2,LA
        FORALL(K=1:KB) SEDDIA50(L,K) = SUM(PERSED(1:NSCM,K,L)*D50(1:NSCM))
      ENDDO
    ENDIF

    GOTO 1000    ! *** JUMP TO FINAL DIAGNOSTICS PRINTOUT

  ENDIF
  ! ***  END HOT START INITIALIZATION

  ! ***********************************************************************************
  ! ***  COLD START INITIALIZATION: IBMECH=0  OR NSEDFLUME OPTION 3
  IF( ( (IBMECH == 0 .OR. ISEDINT <= 1) .AND. .NOT. LSEDZLJ ) .OR. (IHTSTRT == 0 .AND. NSEDFLUME == 3) )THEN
    ! ***  SET POROSITY AND VOID RATIO
    IF( NSEDFLUME == 3 )THEN
      ! ***  CONVERT AND DRY DENSITY OF BED
      ! ***  IBEDDDNU=0, 1 ACTUAL DRY DENSITY, = 2 POROSITY, = 3 VOID RATIO
      IF( IBEDDDNU == 1 )THEN
        DO K=1,KB
          DO L=2,LA
            BEDDINIT(L,K) = 1000.*BEDDINIT(L,K)
          ENDDO
        ENDDO
      ENDIF

      ! ***  CALCULATE POROSITY AND VOID RATIO
      IF( IBEDDDNU <= 1 )THEN
        DO K=1,KB
          DO L=2,LA
            PORBED(L,K) = 0.001*(BEDBINIT(L,K)-BEDDINIT(L,K))
            VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
          ENDDO
        ENDDO
      ENDIF
      IF( IBEDDDNU == 2 )THEN
        DO K=1,KB
          DO L=2,LA
            PORBED(L,K) = BEDDINIT(L,K)
            VDRBED(L,K) = PORBED(L,K)/(1.-PORBED(L,K))
          ENDDO
        ENDDO
      ENDIF
      IF( IBEDDDNU == 3 )THEN
        DO K=1,KB
          DO L=2,LA
            VDRBED(L,K) = BEDDINIT(L,K)
            PORBED(L,K) = VDRBED(L,K)/(1.+VDRBED(L,K))
          ENDDO
        ENDDO
      ENDIF
      DO K=1,KB
        DO L=2,LA
          VDRBED1(L,K) = VDRBED(L,K)
          PORBED1(L,K) = PORBED(L,K)
          HBED(L,K) = 0.0
          HBED1(L,K) = 0.0
          KBT(L) = 1
        ENDDO
      ENDDO

    ELSE
      ! *** ORIGINAL SEDIMENT MODULE
      DO K=1,KB
        DO L=2,LA
          PORBED(L,K) = BEDPORC
          PORBED1(L,K) = BEDPORC
          VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
          VDRBED1(L,K) = PORBED1(L,K)/(1.0-PORBED1(L,K))
          HBED(L,K) = 0.0
          HBED1(L,K) = 0.0
          KBT(L) = 1
        ENDDO
      ENDDO
    ENDIF
    
    ! ***  UNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS  
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY  
    IF( ISEDINT <= 1 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          SEDBT(L,K)=0.0
          SNDBT(L,K)=0.0
        ENDDO
      ENDDO
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO L=2,LA
              HBED(L,K) = HBED(L,K) + SDEN(NS)*SEDB(L,K,NS)
              SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO L=2,LA
              HBED(L,K) = HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
              SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      DO K=1,KB
        DO L=2,LA
          HBED(L,K) = (1.0 + VDRBED(L,K))*HBED(L,K)
          IF( HBED(L,K) > 0.0) KBT(L) = K
        ENDDO
      ENDDO
      DO K=1,KB
        DO L=2,LA
          IF( HBED(L,K) > 0.0 )THEN
            ! *** COMPUTE TOTAL/WET DENSITY
            BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)
          ELSE
            BDENBED(L,K) = 0.0
          ENDIF
        ENDDO
      ENDDO
      DO K=1,KB
        DO L=2,LA
          HBED1(L,K) = HBED(L,K)
          BDENBED1(L,K) = BDENBED(L,K)
        ENDDO
      ENDDO

      ! ***  HANDLE ACTIVE ARMORING LAYER IF NOT PRESENT IN INITIAL OR RESTART
      IF( ISNDAL == 2 .AND. IALSTUP > 0 .AND. KB > 1 )THEN
        DO L=2,LA
          KBT(L) = KBT(L) - 1
          IF( KBT(L) < 1) KBT(L) = 1
          HBED(L,KBT(L)+1)  = 0.
          HBED1(L,KBT(L)+1) = 0.
        ENDDO
      ENDIF

    ENDIF

    ! ***  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS  
    ! ***  AND INITIAL CONDITIONS IN SEDB.INP ARE IN MASS PER UNIT AREA  
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY  
    IF( ISEDINT >= 2 )THEN  
      IF( ISEDBINT == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            SEDBT(L,K) = 0.0
            SNDBT(L,K) = 0.0
          ENDDO
        ENDDO
        IF( ISTRAN(6) >= 1 )THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                HBED(L,K) = HBED(L,K) + SDEN(NS)*SEDB(L,K,NS)
                SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF( ISTRAN(7) >= 1 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                HBED(L,K) = HBED(L,K) + SDEN(NS)*SNDB(L,K,NX)
                SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NX)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        DO K=1,KB
          DO L=2,LA
            HBED(L,K) = (1. + VDRBED(L,K))*HBED(L,K)
            IF( HBED(L,K) > 0.0) KBT(L)=K
          ENDDO
        ENDDO
        DO K=1,KB
          DO L=2,LA
            IF( HBED(L,K) > 0.0 )THEN
              ! *** COMPUTE TOTAL/WET DENSITY
              BDENBED(L,K) = 1000.*PORBED(L,K)  +0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)
            ELSE
              BDENBED(L,K) = 0.0
            ENDIF
          ENDDO
        ENDDO
        DO K=1,KB
          DO L=2,LA
            HBED1(L,K) = HBED(L,K)
            BDENBED1(L,K) = BDENBED(L,K)
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    ! ***  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS
    ! ***  AND INITIAL CONDITIONS  IN SEDB.INP ARE IN MASS FRACTION
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY
    ! ***  THIS OPTION REQUIRES INITIAL LAYER THICKNESSES
    IF( ISEDINT >= 2 )THEN
      IF( ISEDBINT == 1 )THEN
        IF( IBEDLAYU == 1 )THEN
          DO K=1,KB
            DO L=2,LA
              BEDLINIT(L,K) = 0.1*BEDLINIT(L,K)
            ENDDO
          ENDDO
        ENDIF
        DO K=1,KB
          DO L=2,LA
            HBED(L,K) = BEDLINIT(L,K)
            HBED1(L,K) = BEDLINIT(L,K)
            IF( HBED(L,K) > 0.0) KBT(L) = K
          ENDDO
        ENDDO
        DO K=1,KB
          DO L=2,LA
            BDENBED(L,K) = 0.0
          ENDDO
        ENDDO

        ! *** FIRST COMPUTE DRY DENSITY (G/M2)
        IF( ISTRAN(6) >= 1 )THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                BDENBED(L,K) = BDENBED(L,K) + 1000.0*SSG(NS)*SEDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF( ISTRAN(7) >= 1 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO K=1,KB
              DO L=2,LA
                BDENBED(L,K) = BDENBED(L,K) + 1000.0*SSG(NS)*SNDB(L,K,NX)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        ! *** CONVERT DRY DENSITY TO TOTAL/WET DENSITY
        DO K=1,KB
          DO L=2,LA
            BDENBED(L,K) = 1000.0*PORBED(L,K) + (1.0-PORBED(L,K))*BDENBED(L,K)
            BDENBED1(L,K) = BDENBED(L,K)
          ENDDO
        ENDDO

        ! *** TOTAL SEDIMENT FOR ALL CLASSES (G/M2)
        DO K=1,KB
          DO L=2,LA
            SEDBT(L,K) = 1000.0*HBED(L,K)*(BDENBED(L,K) - 1000.0*PORBED(L,K))
          ENDDO
        ENDDO

        ! *** SPLIT TOTAL SEDIMENT INTO CLASSES
        IF( ISTRAN(6) >= 1 )THEN
          DO NS=1,NSED
            DO K=1,KB
              DO L=2,LA
                SEDB(L,K,NS) = SEDB(L,K,NS)*SEDBT(L,K)
                SEDB1(L,K,NS) = SEDB(L,K,NS)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF( ISTRAN(7) >= 1 )THEN
          DO NX=1,NSND
            DO K=1,KB
              DO L=2,LA
                SNDB(L,K,NX) = SNDB(L,K,NX)*SEDBT(L,K)
                SNDB1(L,K,NX) = SNDB(L,K,NX)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    ! ***  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
    IF( ISTRAN(6) >= 1 .AND. .NOT. LSEDZLJ )THEN
      IF( IWRSP(1) == 0 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURS(L,K) = TAUR(1)
            WRSPS(L,K) = WRSPO(1)
          ENDDO
        ENDDO
      ENDIF
      IF( IWRSPB(1) == 0 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURB(L,K) = 1.E6
            WRSPB(L,K) = 0.0
          ENDDO
        ENDDO
      ENDIF
      IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1), VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1))
          ENDDO
        ENDDO
      ENDIF

      ! *** READ IN THE SPATIALLY ASSIGNED COHESIVE CRITICAL BED SHEAR STRESS AND SURFACE EROSION RATE
      IF( IWRSP(1) >= 99 .AND. ISHOUSATONIC == 0 )THEN
        if(process_id == master_id )then
            WRITE(*,'(A)')'READING SSCOHSEDPMAP.INP'
            OPEN(1,FILE='sscohsedpmap.inp')
            OPEN(2,FILE=OUTDIR//'SSCOHSEDPMAP.OUT')
            STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
            READ(1,*)ISSTYPE
            IF( ISSTYPE == 0 )THEN
              DO L=2,LA_Global
                READ(1,*)LD,ID,JD,LSSCOHSED_Global(L)
                RADJCOHSEDS(L) = 1.
              ENDDO
            ELSE
              DO L=2,LA_Global
                READ(1,*)LD,ID,JD,LSSCOHSED_Global(L),RADJCOHSEDS_Global(L)
              ENDDO
            ENDIF
        end if
        ! *** send to all processes
        Call Broadcast_Array(LSSCOHSED_Global, master_id) 
        Call Broadcast_Array(RADJCOHSEDS_Global,master_id)
        
        
        DO L=2,LA_Global
            l_local = Map2Local(l).ll
            if(l_local > 1 )then
               LCORE = LSSCOHSED_Global(L)
               TAUDS(l_local) = TAUDSS(LCORE)
            end if
        ENDDO
        
        IF( NCOHSEDL == 1 )THEN
            DO K=1,KB
                DO L=2,LA_Global
                    l_local = Map2Local(l).ll
                    if(l_local > 1 )then
                        LCORE = LSSCOHSED_Global(L)
                        TAURS(l_local,K) = TAURSS(1,LCORE)
                        TAUNS(l_local,K) = TAUNSS(1,LCORE)
                        WRSPS(l_local,K) = RADJCOHSEDS_Global(L)*WRSPOSS(1,LCORE)
                        TEXPS(l_local,K) = TEXPSS(1,LCORE)
                    end if
                ENDDO
            ENDDO
        ELSE
          DO K=1,KB
              DO L=2,LA_Global
                  l_local = Map2Local(l).ll
                  if(l_local > 1 )then
                      LCORE = LSSCOHSED_Global(L)
                      TAURS(l_local,K) = TAURSS(K,LCORE)
                      TAUNS(l_local,K) = TAUNSS(K,LCORE)
                      WRSPS(l_local,K) = RADJCOHSEDS_Global(L)*WRSPOSS(K,LCORE)
                      TEXPS(l_local,K) = TEXPSS(K,LCORE)
                  end if
              ENDDO
          ENDDO
        ENDIF
        
        if(process_id == master_id )then
            DO L=2,LA
              K = KBT(L)
              WRITE(2,*)L, Map2Global(L).IG, Map2Global(L).JG, TAURS(L,K),WRSPS(L,K)
            ENDDO
            CLOSE(1)
            CLOSE(2)
        end if ! *** end on master
        
      ENDIF
      
      IF( IWRSPB(1) >= 1 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))
            WRSPB(L,K) = CSEDRESB
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    ! ***  SET SEDIMENT VOLUME FRACTIONS
    DO K=1,KB
      DO L=2,LA
        BEDLINIT(L,K) = 0.0
        BEDDINIT(L,K) = 0.0
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)
            VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)
            VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
            ELSE
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0
            ENDIF
            IF( BEDDINIT(L,K) > 0.0 )THEN
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
            ELSE
              VFRBED1(L,K,NS) = 0.0
              BEDDINIT(L,K)  = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
            ELSE
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0
            ENDIF
            IF( BEDDINIT(L,K) > 0.0 )THEN
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
            ELSE
              VFRBED1(L,K,NS) = 0.0
              BEDDINIT(L,K)  = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  ! ***  END COLD START INITIALIZATION: IBMECH=0

  ! ***********************************************************************************
  ! ***  COLD START INITIALIZATION: IBMECH >= 1
  IF( IBMECH >= 1 .AND. ISEDINT > 1 .AND. .NOT. LSEDZLJ )THEN

    ! ***  CONVERT AND INITIALIZE BED LAYER THICKNESS AND DEFINE
    ! ***  INITIAL TOP LAYER
    IF( IBEDLAYU == 1 ) TMPCVT=0.001    ! *** THICKNESS IN MM
    IF( IBEDLAYU == 2 ) TMPCVT=0.01     ! *** THICKNESS IN CM
    IF( IBEDLAYU == 3 ) TMPCVT=1.0      ! *** THICKNESS IN M
    IF( IBEDLAYU >= 1 )THEN
      DO K=1,KB
        DO L=2,LA
          BEDLINIT(L,K) = TMPCVT*BEDLINIT(L,K)
        ENDDO
      ENDDO
    ENDIF
    DO L=2,LA
      KBT(L) = 1
    ENDDO
    DO K=1,KB
      DO L=2,LA
        HBED(L,K) = BEDLINIT(L,K)
        HBED1(L,K) = BEDLINIT(L,K)
        IF( HBED(L,K) > 0.0) KBT(L) = K
      ENDDO
    ENDDO

    ! ***  CONVERT AND INITIALIZE BED BULK DENSITY
    ! ***   IBEDBDNU=0 BEDBINIT IS NOT BULK DENSITY
    ! ***   IBEDBDNU=1 BEDBINIT BULK DENSITY IN KG/M**3
    ! ***   IBEDBDNU=3 BEDBINIT BULK DENSITY IN GM/CM**3
    IF( IBEDBDNU >= 1 )THEN
      IF( IBEDBDNU == 2 )THEN
        DO K=1,KB
          DO L=2,LA
            BEDBINIT(L,K) = 1000.*BEDBINIT(L,K)
          ENDDO
        ENDDO
      ENDIF
      DO K=1,KB
        DO L=2,LA
          BDENBED(L,K) = BEDBINIT(L,K)
          BDENBED1(L,K) = BEDBINIT(L,K)
        ENDDO
      ENDDO
    ENDIF

    ! ***  CONVERT AND DRY DENSITY OF BED
    ! ***  IBEDDDNU=0, 1 ACTUAL DRY DENSITY, = 2 POROSITY, = 3 VOID RATIO
    IF( IBEDDDNU == 1 )THEN
      DO K=1,KB
        DO L=2,LA
          BEDDINIT(L,K) = 1000.*BEDDINIT(L,K)
        ENDDO
      ENDDO
    ENDIF

    ! ***  CALCULATE POROSITY AND VOID RATIO
    IF( IBEDDDNU <= 1 )THEN
      DO K=1,KB
        DO L=2,LA
          PORBED(L,K) = 0.001*(BEDBINIT(L,K)-BEDDINIT(L,K))
          VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
        ENDDO
      ENDDO
    ENDIF
    IF( IBEDDDNU == 2 )THEN
      DO K=1,KB
        DO L=2,LA
          PORBED(L,K) = BEDDINIT(L,K)
          VDRBED(L,K) = PORBED(L,K)/(1.-PORBED(L,K))
        ENDDO
      ENDDO
    ENDIF
    IF( IBEDDDNU == 3 )THEN
      DO K=1,KB
        DO L=2,LA
          VDRBED(L,K) = BEDDINIT(L,K)
          PORBED(L,K) = VDRBED(L,K)/(1.+VDRBED(L,K))
        ENDDO
      ENDDO
    ENDIF
    DO K=1,KB
      DO L=2,LA
        VDRBED1(L,K) = VDRBED(L,K)
        PORBED1(L,K) = PORBED(L,K)
      ENDDO
    ENDDO

    ! ***  INITIALIZE BED SEDIMENT FOR MASS FRACTION INPUT BY CACLULATING
    ! ***  AND STORING TOTAL MASS OF SED/AREA IN BEDDINIT(L,K)
    DO K=1,KB
      DO L=2,LA
        ! ***              M         KG/M3             KG/M3
        BEDDINIT(L,K) = HBED(L,K)*(BDENBED(L,K) - 1000.*PORBED(L,K))     ! *** BEDDINIT is Dry Density in KG/M2
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        IF( ISEDBU(NS) == 1 )THEN
          DO K=1,KB
            DO L=2,LA
              ! *** G/M2     G/KG    FRAC (NO DIM)     KG/M2
              SEDB(L,K,NS) = 1000.*SEDBINIT(L,K,NS)*BEDDINIT(L,K)
              SEDB1(L,K,NS) = SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NS=1,NSND
        IF( ISNDBU(NS) == 1 )THEN
          DO K=1,KB
            DO L=2,LA
              SNDB(L,K,NS) = 1000.*SNDBINIT(L,K,NS)*BEDDINIT(L,K)
              SNDB1(L,K,NS) = SNDB(L,K,NS)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    ! ***  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
    IF( ISTRAN(6) >= 1 )THEN
      IF( IWRSP(1) == 0 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURS(L,K) = TAUR(1)
            WRSPS(L,K) = WRSPO(1)
          ENDDO
        ENDDO
      ENDIF
      IF( IWRSPB(1) == 0 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURB(L,K) = 1.E6
            WRSPB(L,K) = 0.0
          ENDDO
        ENDDO
      ENDIF
      IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1))
          ENDDO
        ENDDO
      ENDIF

      ! *** READ IN THE SPATIALLY ASSIGNED COHESIVE CRITICAL BED SHEAR STRESS AND SURFACE EROSION RATE
      IF( IWRSP(1) >= 99 .AND. ISHOUSATONIC == 0 )THEN
        ! *** Read on master
        if(process_id == master_id )then
          WRITE(*,'(A)')'***WARNING*** SSCOHSEDPMAP input file is not working with MPI'
          WRITE(*,'(A)')'READING SSCOHSEDPMAP.INP'
          OPEN(1, FILE='sscohsedpmap.inp')
          OPEN(2, FILE=OUTDIR//'SSCOHSEDPMAP.OUT')
          STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
          READ(1,*)ISSTYPE
          IF( ISSTYPE == 0 )THEN
            DO L=2,LA_Global
              READ(1,*)LD,ID,JD,LSSCOHSED_Global(L)
              RADJCOHSEDS_Global(L)=1.
            ENDDO
          ELSE
            DO L=2,LA_Global
              READ(1,*)LD,ID,JD,LSSCOHSED_Global(L),RADJCOHSEDS_Global(L)
            ENDDO
          ENDIF
        end if !*** end on master
        
        ! *** send to all processes
        Call Broadcast_Array(LSSCOHSED_Global, master_id)
        Call Broadcast_Array(RADJCOHSEDS_Global, master_id)

        ! *** Map to local
        DO L=2, LA_Global
            l_local = Map2Local(l).ll
            if(l_local > 1 )then
                LCORE = LSSCOHSED_Global(L)
                TAUDS(l_local) = TAUDSS(LCORE)
            end if
        ENDDO
        
        IF( NCOHSEDL == 1 )THEN
          DO K=1,KB
            ! *** Map to local
            DO L=2,LA_Global
              l_local = Map2Local(l).ll
              if(l_local > 1 )then 
                  LCORE = LSSCOHSED_Global(L)
                  TAURS(l_local,K) = TAURSS(1,LCORE)
                  TAUNS(l_local,K) = TAUNSS(1,LCORE)
                  WRSPS(l_local,K) = RADJCOHSEDS_Global(L)*WRSPOSS(1,LCORE)
                  TEXPS(l_local,K) = TEXPSS(1,LCORE)
              end if
              
            ENDDO
          ENDDO
        ELSE
          DO K=1,KB
            ! *** Map to local
            DO L=2,LA_Global
                l_local = Map2Local(l).ll
                if(l_local > 1 )then 
                    LCORE = LSSCOHSED_Global(L)
                    TAURS(l_local,K) = TAURSS(K,LCORE)
                    TAUNS(l_local,K) = TAUNSS(K,LCORE)
                    WRSPS(l_local,K) = RADJCOHSEDS_Global(L)*WRSPOSS(K,LCORE)
                    TEXPS(l_local,K) = TEXPSS(K,LCORE)
                end if
            ENDDO
          ENDDO
        ENDIF
        
        if(process_id == master_id )then
            DO L=2,LA
              K = KBT(L)
              WRITE(2,*)L, Map2Global(L).IG, Map2Global(L).JG, TAURS(L,K), TAUNS(L,K),WRSPS(L,K), TEXPS(L,K)
            ENDDO
            CLOSE(1)
            CLOSE(2)   
        end if ! *** end on master
        
      ENDIF

      IF( IWRSPB(1) >= 1 )THEN
        DO K=1,KB
          DO L=2,LA
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))
            WRSPB(L,K) = CSEDRESB
          ENDDO
        ENDDO
      ENDIF
    ENDIF    ! *** END OF SETTING RESUSPENSION RATES FOR COHESIVES

    ! ***  SET SEDIMENT VOLUME FRACTIONS
    DO K=1,KB
      DO L=2,LA
        BEDLINIT(L,K) = 0.
        BEDDINIT(L,K) = 0.
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)
            VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)
            VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            ! *** BEGIN DSI
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
            ELSE
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0
            ENDIF
            IF( BEDDINIT(L,K) > 0.0 )THEN
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
            ELSE
              VFRBED1(L,K,NS) = 0.0
              BEDDINIT(L,K)  = 0.0
            ENDIF
            ! *** END DSI
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            ! *** BEGIN DSI
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
            ELSE
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0
            ENDIF
            IF( BEDDINIT(L,K) > 0.0 )THEN
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
            ELSE
              VFRBED1(L,K,NS) = 0.0
              BEDDINIT(L,K)  = 0.0
            ENDIF
            ! *** END DSI
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  ! ***  END COLD START INITIALIZATION: IBMECH >= 1

  ! ***********************************************************************************
  ! ***  INITIALIZE BED BOTTOM ELEVATION
  IF( NSEDFLUME /= 1 .AND. NSEDFLUME /= 2 )THEN
    ! ***  INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA
    DO K=1,KB
      DO L=2,LA
        SEDBT(L,K) = 0.0
        SNDBT(L,K) = 0.0
      ENDDO
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            IF( SEDB(L,K,NS) > 0.0 )THEN
              SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
            ELSE
              SEDB(L,K,NS) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NS=1,NSND
        DO K=1,KB
          DO L=2,LA
            IF( SNDB(L,K,NS) > 0.0 )THEN
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NS)
            ELSE
              SNDB(L,K,NS) = 0.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! *** GET TOTAL THICKNESSES AND BEDROCK ELEVATION
    DO L=2,LA
      HBEDA(L)=0.0
    ENDDO
    DO L=2,LA
      IF( LBED(L) ) CYCLE            ! *** HARDBOTTOM BYPASS
      DO K=1,KB
        IF( SEDBT(L,K) <= 1E-6 )THEN
          SEDB(L,K,1:NSED) = 0.0
          SEDBT(L,K) = 0.0
        ENDIF
        IF( SNDBT(L,K) <= 1E-6 )THEN
          SNDB(L,K,1:NSND) = 0.0
          SNDBT(L,K) = 0.0
        ENDIF
        IF( (SEDBT(L,K) + SNDBT(L,K)) <= 1E-6 )THEN
          HBED(L,K) = 0.0
          SEDB(L,K,1:NSED) = 0.0
          SNDB(L,K,1:NSND) = 0.0
        ENDIF
        HBEDA(L) = HBEDA(L) + HBED(L,K)

        ! *** QC
        IF( HBED(L,K) > 0. .AND. PORBED(L,K) > 0. )THEN
          TTHICK = 0.
          DO NS=1,NSED
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SEDB(L,K,NS)*DSEDGMM
          ENDDO
          DO NX=1,NSND
            NS=NSED+NX
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SNDB(L,K,NX)*DSEDGMM
          ENDDO
          TTHICK = TTHICK/(1.0 - PORBED(L,K))
          IF( ABS(HBED(L,K) - TTHICK) > 1E-4 )THEN
            PRINT '(" *** WARNING - BED LAYER THICKNESS IS NOT CONSISTENT.  ADJUSTING THICKNESS (L, K, HOLD, HNEW): ",2I5,2F12.5)', Map2Global(L).LG,K,HBED(L,K),TTHICK
            WRITE(mpi_log_unit,'(" *** WARNING - BED LAYER THICKNESS IS NOT CONSISTENT.  ADJUSTING THICKNESS (L, K, HOLD, HNEW): ",2I5,2F12.5)') Map2Global(L).LG,K,HBED(L,K),TTHICK
            HBED(L,K) = TTHICK
          ENDIF

        ELSE
          HBED(L,K) = 0.0
          PORBED(L,K) = BEDPORC
          VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
          SEDB(L,K,:) = 0.
          SNDB(L,K,:) = 0.
        ENDIF
      END DO

      ! *** Check for empty layers and fill as needed
      IF( NSEDFLUME == 3 .AND. HBED(L,KB) > 0.0 )THEN
        DO K = KB-1,1,-1
          IF( HBED(L,K) == 0.0 )THEN
            PORBED(L,K) = PORBED(L,K+1)
            VDRBED(L,K) = VDRBED(L,K+1)
          ENDIF
        ENDDO
      ENDIF
      
      IF( HBEDA(L) <= 1E-6 ) KBT(L) = 1
    ENDDO
    DO L=2,LA
      ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** HARD BOTTOM ELEVATION
    ENDDO

  ENDIF

  ! ***********************************************************************************
  ! ***  IF N=1 AND ISTRAN(5)=1 CHECK INITIAL TOXIC CONCENTRATIONS IN
  ! ***  BED AND REINITILIZE IF NECESSARY
  IF( ISTRAN(5) >= 1 )THEN
    ! *** GET EACH SOLIDS CLASS FOR EACH TOXIC: NSP2
    DO NT=1,NTOX
      NSP2(NT) = NSED + NSND                      ! *** Kd  Approach ISTOC(NT)=0
      IF( ISTOC(NT) == 1 ) NSP2(NT) = NSP2(NT)+2  ! *** DOC AND POC (POC NON-SEDIMENT RELATED)  (3 Phase)
      IF( ISTOC(NT) == 2 ) NSP2(NT) = NSP2(NT)+1  ! *** DOC AND POC FRACTIONALLY DISTRIBUTED    (3 Phase)
      ! *** ISTOC(NT)=0 and ISTOC(NT)=3                 POC fOC*SED/SND BASED ONLY              (2 Phase)
    END DO

    IF( ISRESTI == 0 .OR. ISCI(5) == 0 )THEN

      ! ***  CALCULATE TOTAL SEDIMENT IN THE BED
      DO K=1,KB
        DO L=1,LC
          SEDBALL(L,K) = 0.
        ENDDO
      ENDDO
      DO K=1,KB
        DO L=1,LC
          SEDBALL(L,K) = SEDBT(L,K)+SNDBT(L,K)
        ENDDO
      ENDDO

      ! *** *******************************************************************C
      !
      ! **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
      !
      ! **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED  (DIMENSIONLESS)
      ! **  TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED  (DIMENSIONLESS)
      CALL CALTOXB_FRACTIONS(1,LASED)

      ! *** COMPUTE FINAL FRACTION (DIMENSIONLESS) ADSORBED ON ALL SOLIDS AND DOC (I.E. NOT DISSOLVED)
      DO NT=1,NTOX
        DO K=1,KB
          DO L=2,LA
            IF( SEDBALL(L,K) > 0.0 )THEN
              ! ***                 M           (   POREWATER (M)        SOLIDS (M)   )
              TOXPFTB(L,K,NT) = TOXPFTB(L,K,NT)/(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
            ELSE
              TOXPFTB(L,K,NT) = 1.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! ***  CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC
      ! ***  CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG
      ! ***  TO TOXB UNITS OF OF MG/M**2
      DO NT=1,NTOX
        IF( ITXBDUT(NT) == 0 )THEN
          ! *** INPUT UNITS OF UG/L ( = MG/M3)
          DO K=1,KB
            DO L=2,LA
              ! *** MG/M2   =   M          MG/M3
              TOXB(L,K,NT)  = HBED(L,K)*TOXB(L,K,NT)
              TOXB1(L,K,NT) = TOXB(L,K,NT)
            ENDDO
          ENDDO
        ENDIF
        IF( ITXBDUT(NT) == 1 )THEN
          ! *** PMC - BEFORE 2017-06 EFDC ASSUMED ONLY PARTICLE MASS IN FILE.  TYPICAL LABRATORY ANALYSIS ARE TOTAL MASS TOXIC/MASS OF SEDIMENT.
          ! ***       ITXBDUT(NT)=1 NOW REPRESENTS THE INPUT TOTAL TOXB AS MG/KG
          ! *** INPUT UNITS OF MG/KG
          DO K=1,KB
            DO L=2,LA
              !TOXB(L,K,NT)  = 0.001*TOXB(L,K,NT)*(SEDBT(L,K) + SNDBT(L,K))/TOXPFTB(L,K,NT)   DEPRECATED
              ! *** MG/M2   =  KG/G    MG/KG     (          G/M2         )
              TOXB(L,K,NT)  = 0.001*TOXB(L,K,NT)*(SEDBT(L,K) + SNDBT(L,K))
              TOXB1(L,K,NT) = TOXB(L,K,NT)
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      ! *** DIAGNOSTICS OF INITIALIZATION
      Call writebreak(mpi_log_unit)
      NS = 0
      WRITE(mpi_log_unit,*) ' Sediment Bed Check: TOXBED.DIA'
      DO NT=1,NTOX
        WRITE(mpi_log_unit,*) 'TOXIC BED FOR CLASS: ',NT
        DO LP=1,LASED
          L=LSED(LP)
          DO K=1,KB
            IF( HBED(L,K) > 0.0 )THEN
              TMP1 = TOXB(L,K,NT)/HBED(L,K)   ! *** MG/M^3
            ELSE
              IF( TOXB(L,K,NT) > 0.0 )THEN
                WRITE(mpi_log_unit,*) 'WARNING - TOXIC IC HAS TOXB > 0 FOR ZERO THICKNESS LAYER: ', Map2Global(L).LG, K, NT, TOXB(L,K,NT)
                TOXB(L,K,NT) = 0.0
                NS = 1
              ENDIF
              TMP1 = 0.0
            ENDIF
            WRITE(mpi_log_unit,2222)  Map2Global(L).IG, Map2Global(L).JG, K, TOXPFTB(L,K,NT), TOXB(L,K,NT), TMP1, TOX(L,KSZ(L),NT)
          ENDDO
        ENDDO
      ENDDO
      IF( NS == 1 ) PRINT *, 'BEDINIT:  Found one or more cells with bad TOX/SED settings.  Check mpi_log_file for details. '        
    ENDIF
  ENDIF

2222 FORMAT(3I5,7E13.4)

  ! ***  INITIALIZE FRACTION OF PARTICULATE ORGANIC CARBON IN BED BY FUNCTION
  IF( ISTPOCB == 4 )THEN
    IVAL=0
    DO NT=1,NTOX
      IF( ISTOC(NT) >= 2)IVAL=1
    ENDDO
    IF( IVAL == 1 )CALL SETFPOCB(0)
  ENDIF

  ! ***  CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS
  IF( .NOT. LSEDZLJ )THEN
    DO K=1,KB
      DO L=2,LA
        IF( K <= KBT(L) )THEN
          FVOLSSD = 1.0/(1.0+VDRBED(L,K))
          FVOLSED = 0.0
          DO NS=1,NSED
            FVOLSED = FVOLSED+VFRBED(L,K,NS)
          ENDDO
          FVOLSND = 0.0
          DO NX=1,NSND
            NS = NSED+NX
            FVOLSND = FVOLSND+VFRBED(L,K,NS)
          ENDDO
          FVOLSED = FVOLSSD*FVOLSED
          FVOLSND = FVOLSSD*FVOLSND
          VDRBEDSND(L,K) = SNDVDRD
          IF( FVOLSED > 1.0E-18 )THEN
            VDRBEDSED(L,K) = ((FVOLSED+FVOLSND)*VDRBED(L,K)-  FVOLSND*SNDVDRD)/FVOLSED
          ELSE
            VDRBEDSED(L,K) = 0.0
          ENDIF
        ELSE
          VDRBEDSND(L,K) = 0.0
          VDRBEDSED(L,K) = 0.0
        ENDIF
      ENDDO
    ENDDO

    ! ***  ADD ACTIVE ARMORING LAYER IF NOT PRESENT IN INITIAL OR RESTART
    IF( ISNDAL == 2 .AND. IALSTUP > 0 .AND. KB > 1 )THEN

      DO L=2,LA
        KTOPTP = KBT(L)     ! *** NEW PARENT LAYER
        KTOPP1 = KBT(L)+1   ! *** NEW ACTIVE LAYER

        TTHICK=0.
        DO K=1,KB
          TTHICK = TTHICK+HBED(L,K)
        ENDDO

        ! *** MAKE SURE THERE IS SEDIMENT IN THE CELL
        IF( TTHICK > 0. )THEN
          IF( KTOPTP == KB )THEN
            ! *** CASE KBT = KB

            IF( HBED(L,KTOPTP) < HBEDAL )THEN
              ! *** ACTIVE LAYER IS NOT THICK ENOUGH. REPARTITIION KBT AND KBT-1

              KT1 = KBT(L)     ! *** ACTIVE LAYER
              KT2 = KBT(L)-1   ! *** PARENT LAYER
              IF( (HBED(L,KT1)+HBED(L,KT2))<HBEDAL )THEN
                WRITE(STR,'(A,I5)')'BEDINIT: ISNDAL=2 AND IALSTUP>0: BED LAYER IS NOT THICK ENOUGH.  L = ', Map2Global(L).LG
                CALL STOPP(STR)
              ENDIF
              HBEDP = HBEDAL

              FRACT1  = HBEDP/HBED(L,KT1)
              FRACT2  = (HBED(L,KT2)-(HBEDP-HBED(L,KT1)))/HBED(L,KT2)
              FRACT2C = 1.-FRACT2

              ! *** INCREASE TOP LAYER TO HBEDP
              HBED(L,KT1 ) = HBEDP
              SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)
              SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)
              DO NS=1,NSED
                SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)
                SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)
              ENDDO
              DO NS=1,NSND
                SNDB(L,KT1,NS)  = FRACT1*SNDB(L,KT1,NS)
                SNDB1(L,KT1,NS) = FRACT1*SNDB1(L,KT1,NS)
              ENDDO
              DO NT=1,NTOX
                TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)
                TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)
              ENDDO

              ! *** REDUCE LAYER BELOW
              HBED(L,KT2)  = FRACT2*HBED(L,KT2)
              SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)
              SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)
              DO NS=1,NSED
                SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)
                SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)
              ENDDO
              DO NS=1,NSND
                SNDB(L,KT2,NS)  = FRACT2*SNDB(L,KT2,NS)
                SNDB1(L,KT2,NS) = FRACT2*SNDB1(L,KT2,NS)
              ENDDO
              DO NT=1,NTOX
                TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)
                TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)
              ENDDO
            ELSE
              ! *** ACTIVE LAYER IS TOO THICK. REPARTITIION KBT AND KBT-1
              KT1 = KBT(L)     ! *** ACTIVE LAYER
              KT2 = KBT(L)-1   ! *** PARENT LAYER

              HBEDP = HBEDAL

              FRACT1  = HBEDP/HBED(L,KT1)
              FRACT2  = (HBED(L,KT2)+(HBED(L,KT1)-HBEDP))/HBED(L,KT2)
              FRACT2C = 1.-FRACT2

              ! *** DECREASE TOP LAYER TO HBEDP
              HBED(L,KT1)  = HBEDP
              SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)
              SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)
              DO NS=1,NSED
                SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)
                SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)
              ENDDO
              DO NS=1,NSND
                SNDB(L,KT1,NS)  = FRACT1*SNDB(L,KT1,NS)
                SNDB1(L,KT1,NS) = FRACT1*SNDB1(L,KT1,NS)
              ENDDO
              DO NT=1,NTOX
                TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)
                TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)
              ENDDO

              ! *** INCREASE THE LAYER BELOW
              HBED(L,KT2)  = FRACT2*HBED(L,KT2)
              SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)
              SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)
              DO NS=1,NSED
                SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)
                SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)
              ENDDO
              DO NS=1,NSND
                SNDB(L,KT2,NS)  = FRACT2*SNDB(L,KT2,NS)
                SNDB1(L,KT2,NS) = FRACT2*SNDB1(L,KT2,NS)
              ENDDO
              DO NT=1,NTOX
                TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)
                TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)
              ENDDO

            ENDIF

          ELSE

            ! *** STANDARD CASE OF KBT<KB
            IF( HBED(L,KTOPTP) < HBEDAL )THEN
              ! *** PARENT LAYER IS NOT THICK ENOUGH. REPARTITIION KBT AND KBT-1

              KT1 = KBT(L)     ! *** PARENT LAYER
              KT2 = KBT(L)-1   ! *** LAYER BELOW PARENT
              IF( KT2 > 0 )THEN
                HBEDP = MIN(HBEDAL*2.,HBED(L,KT2))

                FRACT1  = HBEDP/HBED(L,KT1)
                FRACT2  = (HBED(L,KT2)-(HBEDP-HBED(L,KT1)))/HBED(L,KT2)
                FRACT2C = 1.-FRACT2

                ! *** INCREASE TOP LAYER TO HBEDP
                HBED(L,KT1) = HBEDP
                IF( ISTRAN(6) > 0 )THEN
                  SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)
                  DO NS=1,NSED
                    SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)
                    SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)
                  ENDDO
                ENDIF
                IF( ISTRAN(7) > 0 )THEN
                  SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)
                  DO NS=1,NSND
                    SNDB(L,KT1,NS)  = FRACT1*SNDB(L,KT1,NS)
                    SNDB1(L,KT1,NS) = FRACT1*SNDB1(L,KT1,NS)
                  ENDDO
                ENDIF
                IF( ISTRAN(5) > 0 )THEN
                  DO NT=1,NTOX
                    TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)
                    TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)
                  ENDDO
                ENDIF

                ! *** REDUCE LAYER BELOW
                HBED(L,KT2) = FRACT2*HBED(L,KT2)
                IF( ISTRAN(6) > 0 )THEN
                  SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)
                  DO NS=1,NSED
                    SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)
                    SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)
                  ENDDO
                ENDIF
                IF( ISTRAN(7) > 0 )THEN
                  SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)
                  DO NS=1,NSND
                    SNDB(L,KT2,NS)  = FRACT2*SNDB(L,KT2,NS)
                    SNDB1(L,KT2,NS) = FRACT2*SNDB1(L,KT2,NS)
                  ENDDO
                ENDIF
                IF( ISTRAN(5) > 0 )THEN
                  DO NT=1,NTOX
                    TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)
                    TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF

            ! *** SPLIT THE PARENT LAYER INTO THE PARENT (KTOPTP) AND THE ACTIVE (KTOPP1)
            FRACACT(L) = HBEDAL/HBED(L,KTOPTP)
            FRACPAR(L) = (HBED(L,KTOPTP)-HBEDAL)/HBED(L,KTOPTP)
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

            IF( ISTRAN(6) > 0 )THEN
              DO NS=1,NSED
                SEDB(L,KTOPP1,NS) = FRACACT(L)*SEDB(L,KTOPTP,NS)
                SEDB1(L,KTOPP1,NS) = FRACACT(L)*SEDB1(L,KTOPTP,NS)
                SEDB(L,KTOPTP,NS) = FRACPAR(L)*SEDB(L,KTOPTP,NS)
                SEDB1(L,KTOPTP,NS) = FRACPAR(L)*SEDB1(L,KTOPTP,NS)
                STFPOCB(L,KTOPP1,NS) = STFPOCB(L,KTOPTP,NS)
              ENDDO
            ENDIF
            IF( ISTRAN(7) > 0 )THEN
              DO NS=1,NSND
                NX=NSED+NS
                SNDB(L,KTOPP1,NS) = FRACACT(L)*SNDB(L,KTOPTP,NS)
                SNDB1(L,KTOPP1,NS) = FRACACT(L)*SNDB1(L,KTOPTP,NS)
                SNDB(L,KTOPTP,NS) = FRACPAR(L)*SNDB(L,KTOPTP,NS)
                SNDB1(L,KTOPTP,NS) = FRACPAR(L)*SNDB1(L,KTOPTP,NS)
                STFPOCB(L,KTOPP1,NX) = STFPOCB(L,KTOPTP,NX)
              ENDDO
            ENDIF
            IF( ISTRAN(5) > 0 )THEN
              DO NT=1,NTOX
                TOXB(L,KTOPP1,NT) = FRACACT(L)*TOXB(L,KTOPTP,NT)
                TOXB1(L,KTOPP1,NT) = FRACACT(L)*TOXB1(L,KTOPTP,NT)
                TOXB(L,KTOPTP,NT) = FRACPAR(L)*TOXB(L,KTOPTP,NT)
                TOXB1(L,KTOPTP,NT) = FRACPAR(L)*TOXB1(L,KTOPTP,NT)
                TOXPFTB(L,KTOPP1,NT) = TOXPFTB(L,KTOPTP,NT)
              ENDDO
              DO NT=1,NTOX
                DO NS=1,NSED+NSND+2
                  TOXPFB(L,KTOPP1,NS,NT) = TOXPFB(L,KTOPTP,NS,NT)
                ENDDO
              ENDDO
            ENDIF

            ! *** UPDATE TOP LAYER
            KBT(L) = KBT(L)+1

          ENDIF  ! *** END OF KBT<KB CHECK
        ENDIF    ! *** END OF TOTAL THICKNESS CHECK: TTHICK>0
      ENDDO  ! *** END OF LA LOOP

    ENDIF    ! *** END OF REPARTITIONING ACTIVE ARMORING LAYER AT STARTUP IF REQUESTED

    ! ***********************************************************************!
    ! ***  ADJUST POROSITY AND VOID RATIO FOR IBMECH == 99
    IF( IBMECH == 99 .AND. .NOT. LSEDZLJ )THEN
      DO K=1,KB
        DO L=2,LA
          FRACCOH(L,K) = 0.0
          FRACNON(L,K) = 0.0
        ENDDO
      ENDDO

      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            IF( K <= KBT(L) )THEN
              FRACCOH(L,K) = FRACCOH(L,K)+VFRBED(L,K,NS)
              FRACNON(L,K) = FRACNON(L,K)+VFRBED1(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO K=1,KB
        DO L=2,LA
          IF( K <= KBT(L) )THEN
            PORBED(L,K) = BMECH1*(FRACCOH(L,K)**BMECH2)+BMECH3
            PORBED1(L,K) = BMECH1*(FRACNON(L,K)**BMECH2)+BMECH3
            VDRBED(L,K) = PORBED(L,K)/(1.-PORBED(L,K))
            VDRBED1(L,K) = PORBED1(L,K)/(1.-PORBED1(L,K))
          ENDIF
        ENDDO
      ENDDO

    ENDIF

  ENDIF   ! *** END OF LSEDZLJ BYPASS

  IF( .NOT. LSEDZLJ )THEN
    ! *** *******************************************************************C
    ! *** SET MEAN D50 AND D90 (meters)
    IF( ISTRAN(7) >= 1 )THEN
      DO K=1,KB
        DO L=2,LA
          SEDDIA50(L,K)=0.
          SEDDIA90(L,K)=0.
          SNDBT(L,K)=0.
        ENDDO
      ENDDO
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            SEDDIA50(L,K) = SEDDIA50(L,K) + SNDB(L,K,NX)*LOG(SEDDIA(NS))   ! *** SNDDIA is in meters
            SNDBT(L,K)    = SNDBT(L,K) + SNDB(L,K,NX)
          ENDDO
        ENDDO
      ENDDO
      DO K=1,KB
        DO L=2,LA
          IF( SNDBT(L,K) > 0. )THEN
            SEDDIA50(L,K) = SEDDIA50(L,K)/SNDBT(L,K)
          ENDIF
        ENDDO
      ENDDO

      DO K=1,KB
        DO L=2,LA
          IF( SNDBT(L,K) > 0.1 )THEN  ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
            SEDDIA50(L,K) = EXP(SEDDIA50(L,K))
          ELSE
            SEDDIA50(L,K) = SEDDIA(NSED+1)
          ENDIF
        ENDDO
      ENDDO

    ENDIF  ! *** END OF NON-COHESIVE GRAIN SIZE CALCS

  ELSEIF( NSEDFLUME == 3 .AND. IHTSTRT == 0 )THEN
    ! *** CLEAN UP SPATIALLY VARYING INPUTS TO WORK WITH SEDZLJ
    DO LP=1,LASED
      L=LSED(LP)
      DO K=KB,1,-1
        IF( HBED(L,K) > 0.0 )THEN
          BULKDENS(K,L) = SEDBT(L,K)/HBED(L,K)*1.E-6   ! *** Convert g/m^3 to g/cm^3
          LAYERACTIVE(K,L) = 2                         ! *** Flag original in-place sediment layers
        ELSE
          IF( K < KB )THEN
            TMP1 = BULKDENS(K+1,L)
          ELSE
            TMP1 = 2.65*(1.-BEDPORC)
          ENDIF
          BULKDENS(K,L) = TMP1
          LAYERACTIVE(K,L) = 0
        ENDIF

        DO NS=1,NSCM
          IF( SEDBT(L,K) > 0.0 )THEN
            PERSED(NS,K,L) = SEDB(L,K,NS)/SEDBT(L,K)
          ELSE
            PERSED(NS,K,L) = 0.0
          ENDIF
        ENDDO
        DTOTAL = SUM(PERSED(:,K,L))
        IF( DTOTAL > 0. )THEN
          PERSED(1:NSCM,K,L) = PERSED(1:NSCM,K,L)/DTOTAL   ! *** Ensure precision for mass balance
        ELSE
          PERSED(1:NSCM,K,L) = 0.0
        ENDIF
        TSED(K,L)  = SEDBT(L,K)*1.E-4                      ! *** SEDBT is in units of g/m^2 and TSED in units of g/cm^2.
        TSED0(K,L) = TSED(K,L)
      ENDDO
    ENDDO

    ! *** Initialize KBT to the first layer with mass
    DO LP=1,LASED
      L = LSED(LP)
      KBT(L) = -1
      DO K=1,KB
        IF( TSED0(K,L) > 0. )THEN
          KBT(L) = K
          EXIT
        ENDIF
      ENDDO
      IF( KBT(L) == -1 ) KBT(L) = KB
    ENDDO

  ENDIF    ! *** END OF LSEDZLJ BYPASS

  ! ***********************************************************************************************************
1000 CONTINUE     ! JUMP FROM THE HOTSTART (SEDZLJ INITIALIZATION CONTINUES FROM THE HOTSTART SECTION)

  DO K=1,KB
    DO L=2,LA
      IF( (K <= KBT(L) .AND. .NOT. LSEDZLJ) .OR. (K >= KBT(L) .AND. LSEDZLJ) )THEN
        SDENAVG(L,K) = (BDENBED(L,K)-1000.0*PORBED(L,K))/(1.0-PORBED(L,K))
      ELSE
        SDENAVG(L,K) = (BDENBED(L,KBT(L))-1000.0*PORBED(L,KBT(L)))/(1.0-PORBED(L,KBT(L)))
      ENDIF
      IF( SDENAVG(L,K) <= 0. ) SDENAVG(L,K) = 0.0
    ENDDO
  ENDDO

  if(master_id == process_id )THEN
    ! ***********************************************************************************************************
    ! ***  WRITE DIAGNOSTIC FILES FOR BED INITIALIZATION
    WRITE(6,'(A)') 'WRITE DIAGNOSTIC FILES FOR BED INITIALIZATION INTO mpi_log_file'
  endif
  
  Call WriteBreak(mpi_log_unit)
  WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.ELV'
  WRITE(mpi_log_unit,111)
  DO L=2,LA
    SURF = HP(L) + BELV(L)
    WRITE(mpi_log_unit,101) Map2Global(L).IG, Map2Global(L).JG, ZELBEDA(L),HBEDA(L),BELV(L),HP(L),SURF
  ENDDO

  Call WriteBreak(mpi_log_unit)
  WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.ZHB'
  WRITE(mpi_log_unit,112) (K,K=1,KB)
  DO L=2,LA
    WRITE(mpi_log_unit,101) Map2Global(L).IG, Map2Global(L).JG, ZELBEDA(L),HBEDA(L),(HBED(L,K),K=1,KB)
  ENDDO

  IF( ISTRAN(6) > 0 )THEN
    Call WriteBreak(mpi_log_unit)
    WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.SED'
    WRITE(mpi_log_unit,113) (K,K=1,KB)
    DO L=2,LA
      WRITE(mpi_log_unit,101)  Map2Global(L).IG, Map2Global(L).JG, (SEDBT(L,K),K=1,KB), SUM(SEDBT(L,:))
    ENDDO
  ENDIF
    
  IF( ISTRAN(7) > 0 )THEN
    Call WriteBreak(mpi_log_unit)
    WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.SND'
    WRITE(mpi_log_unit,114) (K,K=1,KB)
    DO L=2,LA
      WRITE(mpi_log_unit,101)  Map2Global(L).IG, Map2Global(L).JG, (SNDBT(L,K),K=1,KB), SUM(SNDBT(L,:))
    ENDDO
  ENDIF
    
  IF( ISTRAN(5) > 0 )THEN
    Call WriteBreak(mpi_log_unit)
    WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.TOX'
    DO NT=1,NTOX
      WRITE(mpi_log_unit,115) NT, (K,K=1,KB)
      DO L=2,LA
        WRITE(mpi_log_unit,101) Map2Global(L).IG, Map2Global(L).JG, (TOXB(L,K,NT),K=1,KB)
      ENDDO
    ENDDO
  ENDIF
    
  Call WriteBreak(mpi_log_unit)
  WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.VDR'
  WRITE(mpi_log_unit,116) (K,K=1,KB)
  DO L=2,LA
    WRITE(mpi_log_unit,101) Map2Global(L).IG, Map2Global(L).JG, (VDRBED(L,K),K=1,KB)
  ENDDO

  Call WriteBreak(mpi_log_UNIT)
  WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.POR'
  WRITE(mpi_log_unit,117) (K,K=1,KB)
  DO L=2,LA
    WRITE(mpi_log_unit,101) Map2Global(L).IG, Map2Global(L).JG, (PORBED(L,K),K=1,KB)
  ENDDO

  Call WriteBreak(mpi_log_unit)
  WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.BDN'
  WRITE(mpi_log_unit,118) (K,K=1,KB)
  DO L=2,LA
    WRITE(mpi_log_unit,101) Map2Global(L).IG, Map2Global(L).JG, (BDENBED(L,K),K=1,KB)
  ENDDO

  Call WriteBreak(mpi_log_unit)
  WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINIT.VR'
  WRITE(mpi_log_unit,119) (K,K=1,KB)
  DO L=2,LA
    WRITE(mpi_log_unit,191) Map2Global(L).IG, Map2Global(L).JG, (VDRBED(L,K),K=1,KB)
  ENDDO

  IF( ISTRAN(6) > 0 )THEN
    Call WriteBreak(mpi_log_unit)
    WRITE(mpi_log_unit,*) ' Sediment Bed Check: BEDINITC.VV'
    WRITE(mpi_log_unit,120)
    DO L=2,LA
      K=KBT(L)
      WRITE(mpi_log_unit,191,IOSTAT=ISO) Map2Global(L).IG, Map2Global(L).JG,  VFRBED(L,K,1), PORBED(L,K), VDRBED(L,K),  VDRBEDSED(L,K)
    ENDDO
  ENDIF
  ! *** END DIAGNOSTIC FILE DUMP

  ! *** SAVE STARTUP BED PARAMETERS FOR LATER
  DO K=1,KB
    DO L=2,LA
      IF( K <= KBT(L) )THEN
        VDRBED2(L,K) = VDRBED(L,K)
      ELSE
        VDRBED2(L,K) = VDRBED(L,KBT(L))
      ENDIF
    ENDDO
  ENDDO

  DO L=2,LA
    DO K=1,KB
      BEDTHKSV(L,K) = HBED(L,K)
      BEDPORSV(L,K) = PORBED(L,K)
      BEDVDRSV(L,K) = VDRBED(L,K)
      BEDBKDSV(L,K) = BDENBED(L,K)
    ENDDO
  ENDDO

101 FORMAT(2I5,100E13.5)
102 FORMAT(10X,100E13.5)
    
111 FORMAT('   IL   JL        ZBEDB        HBEDT         BELV        HWCOL         SELV')
112 FORMAT(' ZBEDB HBEDT HBED(K=1,KB)',/,'   IL   JL        ZBEDB        HBEDT',100I13)
113 FORMAT(' SEDBT(K=1,KB)'  ,/,'   IL   JL',100I13)
114 FORMAT(' SNDBT(K=1,KB)',  /,'   IL   JL',100I13)
115 FORMAT(' TOXB(K=1,KB,NT)  NT = ',I5,/,'   IL   JL',100I10)
116 FORMAT(' VRDBED(K=1,KB)', /,'   IL   JL',100I13)
117 FORMAT(' PORBED(K=1,KB)', /,'   IL   JL',100I13)
118 FORMAT(' BDENBED(K=1,KB)',/,'   IL   JL',20X,100I13)
119 FORMAT(' VDRBED(K=1,KB)', /,'   IL   JL',100I10)
120 FORMAT('   IL   JL    VFRBED    PORBED      BELV    VDRBED VDRBEDSED')

191 FORMAT(2I5,100F10.3)
192 FORMAT(10X,100F10.3)

RETURN

  END

