! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEDIC
  
  USE GLOBAL
  USE INFOMOD,ONLY:SKIPCOM,READSTR
  Use Variables_MPI
#ifdef _MPI
  Use Broadcast_Routines
  Use Variables_MPI_Write_Out
#endif

  IMPLICIT NONE
  
  INTEGER :: CORE, INCORE, I, J, K, KT, L, LG, LL, M, NS, NT, VAR_BED, FDIR, NWV, NSC, ILocal, JLocal, SURFACE 
  INTEGER :: IWV, JWV, NSKIP, KTOP
  CHARACTER(LEN=80)  :: STR_LINE
  CHARACTER(LEN=120) :: STR_120
  
  !PT- real values are written in DOUBLE PRECISION. 7/16/08
  REAL(RKD) :: DTOTAL, STWVHTMP, STWVTTMP, STWVDTMP, TSUM
  REAL(RKD),ALLOCATABLE,DIMENSION(:,:)   :: BDEN          ! *** (INCORE,KB)
  REAL(RKD),ALLOCATABLE,DIMENSION(:,:)   :: TAUTEMP       ! *** (KB)
  REAL(RKD),ALLOCATABLE,DIMENSION(:,:,:) :: PNEW          ! *** (INCORE,KB,NSCM)
  REAL(RKD),ALLOCATABLE,DIMENSION(:,:)   :: TSED0S        ! *** (KB,INCORE)
  REAL(RKD),ALLOCATABLE,DIMENSION(:,:,:) :: ERATETEMP     ! *** (INCORE,KB,ITBM)

  ! Reads in Initial Erosion data.
  ! Reads in Erosion Data for Newly deposited 
  ! Sediments.
  ! Calculates bed parameters.
  ! REVISION DATE :  May 24, 2006
  ! Craig Jones and Scott James
  ! *** ***********************************************************************
  ! Set Maximum Number of Size Clases to be read in
  ! for newly deposited sediment erosion rates NSICM 
  ! Open Input Files
  ! Craig Jones
  ! *************************************************************************
  ! 2016-12
  ! Rearranged SEDZLJ initialization and added parameters for better toxics 
  ! similations.  Deprecated NSEDFLUME = 99 (i.e. SEDZLJ toxics) 
  ! Toxics are handled by ISTRAN(5)>0.  NEQUIL no longer used
  ! Paul M. Craig
 
  if( process_id == master_id )THEN
    WRITE(*,'(A)')'READING SEDZLJ FILES'
  
    !CALL SEDDATA !calls routine to convert SEDflume data into a form useable by this code
    OPEN(UNIT=10,FILE='erate.sdf')
    OPEN(UNIT=30,FILE='bed.sdf')
  
    ! Read in Sediment Transport Parameters
    ! VAR_BED  = 1 for variable sediment bed
    ! ICALC_BL > 0 for bedload  calculation
    READ (30,'(A80)') STR_LINE
    READ (30,*) VAR_BED, KB, ICALC_BL, SEDSTEP, SEDSTART, IHTSTRT, IMORPH, ISWNWAVE, MAXDEPLIMIT, HPMIN
    IF( HPMIN < 0.003 .OR. HPMIN >= 1.0 ) HPMIN = 0.25
  
    IF( SEDSTEP < TIDALP/REAL(NTSPTC) ) SEDSTEP = TIDALP/REAL(NTSPTC)
    DTSED = SEDSTEP
    DTSEDJ = REAL(DTSED,8)
  
    IF( SEDSTART <= TBEGIN ) SEDSTART=TBEGIN
  
    ! *** NSCM  = Maximum number of grainsize classes defined.  NSCM=NSED most times.
    ! *** ITBM  = Number of Sedflume Shear Categories
    ! *** NSICM = Maximum number of redeposited grainsize classes defined.  NSICM=NSCM most times.
    READ (30,'(A80)') STR_LINE 
    READ (30,'(A80)') STR_LINE    ! DATA READ BY SCANSEDZLJ; ITBM and NSICM.
  
    READ (30,'(A80)') STR_LINE
    READ (30,*)  ZBSKIN,TAUCONST,ISSLOPE,BEDLOAD_CUTOFF

    ! *** Median grain size for each size class
    READ (30,'(A80)') STR_LINE
    READ (30,*) (D50(K),K=1,NSCM)  

    ! *** Erosion
    READ (30,'(A80)') STR_LINE
    READ (30,*) (TCRE(K),K=1,NSCM) 

    ! *** Suspension
    READ (30,'(A80)') STR_LINE
    READ (30,*) (TCRSUS(K),K=1,NSCM)

    ! *** Settling Velocities
    READ (30,'(A80)') STR_LINE
    READ (30,*) (DWSIN(K),K=1,NSCM)

    ! Hydrophobic Contaminant Information (PMC - 2016-11-09 - Deprecated - Handled by CALTOX/CALTOXB)
    !IF( NSEDFLUME == 2 )THEN
      !READ (30,'(A80)') STR_LINE
      !READ (30,*) (KPART(K),K=1,NSCM)
      !READ (30,'(A80)') STR_LINE
      !READ (30,*) (DIFFCOFF(K),K=1,NSCM)
      !READ (30,'(A80)') STR_LINE
      !DO LL=1,KB
      !   READ(30,*) (PCONTEMP(K,LL),K=1,NSCM)
      !ENDDO
    !ENDIF
  
    ! *** ***********************************************************************
    !Reading in Erate.sdf starting with the layer's thickness.
    READ (10,'(A80)') STR_LINE
    READ (10,*) TACTM !read in active layer multiplier
    ! Read in Initial Erosion Data
  end if !***End calculation on master process

  Call Broadcast_Scalar(VAR_BED,  master_id)
  Call Broadcast_Scalar(KB,       master_id)
  Call Broadcast_Scalar(ICALC_BL, master_id)
  Call Broadcast_Scalar(SEDSTEP,  master_id)
  Call Broadcast_Scalar(SEDSTART, master_id)
  Call Broadcast_Scalar(IHTSTRT,  master_id)
  Call Broadcast_Scalar(IMORPH,   master_id)
  Call Broadcast_Scalar(ISWNWAVE, master_id)
  Call Broadcast_Scalar(MAXDEPLIMIT, master_id)
  Call Broadcast_Scalar(HPMIN,    master_id)
  Call Broadcast_Scalar(DTSED,    master_id)
  Call Broadcast_Scalar(DTSEDJ,   master_id)

  Call Broadcast_Scalar(ZBSKIN,   master_id)
  Call Broadcast_Scalar(TAUCONST, master_id)
  Call Broadcast_Scalar(ISSLOPE,  master_id)
  Call Broadcast_Scalar(BEDLOAD_CUTOFF, master_id)
  
  Call Broadcast_Scalar(TACTM, master_id)

  Call Broadcast_Array(D50, master_id)
  Call Broadcast_Array(TCRE, master_id)
  Call Broadcast_Array(TCRSUS, master_id)
  Call Broadcast_Array(DWSIN, master_id)

  IF( VAR_BED >= 1 )THEN    
    ! Variable Bed *************************************************
    ALLOCATE(I2D_Global(IC_Global,JC_Global))
    I2D_Global = 0
    
    if( process_id == master_id )then
      OPEN(UNIT=20,FILE='core_field.sdf')

      ! Determine File Format
      READ(20,'(A120)')STR_120
      I = 1
      DO WHILE ( i < 117 )
        IF( STR_120(I:I+2) == 'DSI' .OR. STR_120(I:I+2) == 'dsi' )EXIT
        I = I + 1
      END DO
     
      CLOSE(20)
      
      ! *** READ THE DATA
      OPEN(UNIT=20,FILE='core_field.sdf')
      IF( I < 117 )THEN
        ! *** DSI STANDARD
        READ(20,'(A80)') STR_LINE
        READ(20,'(A80)') STR_LINE
        READ(20,'(A80)') STR_LINE
        READ(20,*) INCORE !read the number of cores    
        READ(20,'(A80)') STR_LINE
        READ(20,'(A80)') STR_LINE
        DO LG=2,LA_Global
          READ(20,*) I, J, CORE
          I2D_Global(I,J) = CORE                                ! *** NCORENO
        ENDDO
      ELSE
        ! *** SNL STANDARD
        READ(20,*)INCORE !read the number of cores    
        DO J=JC,1,-1
          READ(20,'(120(I1,1X))')(I2D_Global(I,J),I=1,IC)       ! *** NCORENO
        ENDDO
      ENDIF
      CLOSE(20)
    end if !***End calculation on master process
    call Broadcast_Scalar(INCORE,    master_id)
    call Broadcast_Array(I2D_Global, master_id)
    
    ! *** Map to local domain
    DO J=3,JC_GLOBAL-2
      DO I=3,IC_GLOBAL-2
        LG = LIJ_Global(I,J)
        IF( LG > 0 )THEN
          L = Map2Local(LG).LL
          IF( L > 0 )THEN
            ILocal = Map2Local(LG).IL
            JLocal = Map2Local(LG).JL
            NCORENO(ILocal,JLocal) = I2D_Global(I,J)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(I2D_Global)
  ELSE
    ! Constant Erosion in Horizontal ********************************  
    INCORE = 1

    DO L=2,LA
      I = IL(L)
      J = JL(L)
      NCORENO(I,J) = 1
    ENDDO
  ENDIF
 
  ! *** ALLOCATIONS AND INITIALIZATIONS
  IF( NSEDFLUME == 1 )THEN
    ALLOCATE(ERATETEMP(INCORE,KB,ITBM))
    ERATETEMP = 0.0
  ELSE
    ITBM = 2
    ALLOCATE(EA(INCORE,KBM))
    ALLOCATE(EN(INCORE,KBM))
    ALLOCATE(MAXRATE(INCORE,KBM))
    EA = 0.0
    EN = 0.0
    MAXRATE = 10000.0
  ENDIF

  ALLOCATE(ACTDEPA(NSICM))
  ALLOCATE(ACTDEPN(NSICM))
  ALLOCATE(ACTDEPMAX(NSICM))
  ALLOCATE(BDEN(INCORE,KBM))
  ALLOCATE(ERATEND(NSICM,ITBM))
  ALLOCATE(ERATE(KB,LCM,ITBM))
  ALLOCATE(PNEW(INCORE,KBM,NSCM+1))
  ALLOCATE(SEDDENS(INCORE))
  ALLOCATE(TAUTEMP(INCORE,KBM))
  ALLOCATE(TSED0S(KB,INCORE))
  ALLOCATE(TAULOC(ITBM))

  ACTDEPA=0.0
  ACTDEPN=0.0
  ACTDEPMAX=10000.0
  BDEN=0.0   
  ERATEND=0.0 
  ERATE=0.0   
  PNEW=0.0
  SEDDENS=0.0 
  TAUTEMP=0.0
  TSED0S = 0.0
  TAULOC=0.0    

  if( process_id == master_id )THEN
    ! *** **********************************************************   
    DO CORE=1,INCORE  ! for each core of data      
      READ (10,'(A80)') STR_LINE
      READ(10,*)(TAUTEMP(CORE,K),K=1,KB)           ! *** read the critical shear stresses of the core
      READ (10,'(A80)') STR_LINE
      READ (10,*) (TSED0S(K,CORE),K=1,KB)          ! *** read in layer thickness 
      READ (10,'(A80)') STR_LINE      
      READ(10,*) (BDEN(CORE,K),K=1,KB)             ! *** read in the bulk density of the core 
      READ (10,'(A80)') STR_LINE      
      READ(10,*) WATERDENS, SEDDENS(CORE)          ! *** read in the water density and sediment solid's density
      READ (10,'(A80)') STR_LINE  
      DO K=1,KB
        READ(10,*)(PNEW(CORE,K,NS),NS=1,NSCM)
      ENDDO       
      READ (10,'(A80)') STR_LINE
      IF( NSEDFLUME == 1 )THEN
        DO M=1,ITBM
          READ(10,*)TAULOC(M)                      ! *** shear stress used to erode a portion of the core
          READ(10,*)(ERATETEMP(CORE,K,M),K=1,KB)   ! *** erosion rate for each layer subject to shear stress TAULOC
        ENDDO
      ELSE
        TAULOC(1) = 0.0
        TAULOC(2) = 1000.                          ! *** Set to a very large maximum shear stress
        DO K=1,KB
          READ(10,*) EA(CORE,K), EN(CORE,K), MAXRATE(CORE,K)
        ENDDO
      ENDIF
    ENDDO     
  end if !***End calculation on master process

  Call Broadcast_Scalar(WATERDENS,  master_id)

  Call Broadcast_Array(TAUTEMP,     master_id)
  Call Broadcast_Array(TSED0S,      master_id)
  Call Broadcast_Array(BDEN,        master_id)
  Call Broadcast_Array(SEDDENS,     master_id)
  Call Broadcast_Array(PNEW,        master_id)
  IF( NSEDFLUME == 1 )THEN
    Call Broadcast_Array(TAULOC,    master_id)
    Call Broadcast_Array(ERATETEMP, master_id)
  ELSE
    Call Broadcast_Array(TAULOC,    master_id)
    Call Broadcast_Array(EA,        master_id)
    Call Broadcast_Array(EN,        master_id)
    Call Broadcast_Array(MAXRATE,   master_id)
      
    DO K=1,KB
      DO CORE=1,INCORE  ! for each core of data      
        MAXRATE(CORE,K) = MAXRATE(CORE,K)*BDEN(CORE,K)                             ! *** Convert from cm/s to g/m2/s using dry bulk density
      ENDDO
    ENDDO
  ENDIF
  
  
  ! ***
  DO L=2,LA
    I=IL(L)   ! *** I location as a function of L
    J=JL(L)   ! *** J location as a function of L
    
    IF( NCORENO(I,J) > 0 )THEN
      CORE = NCORENO(I,J)
      DO K=1,KB
        TAUCOR(K,L) = TAUTEMP(CORE,K)                    ! *** Critical shear stresses from cores
        IF( NSEDFLUME == 1 )THEN
          DO M=1,ITBM
            ERATE(K,L,M) = ERATETEMP(CORE,K,M)           ! *** Set erosion rate to measured value
          ENDDO
        ENDIF
        
        IF( NSEDFLUME == 3 )THEN
          
        ELSE
          ! *** NSEDFLUME = 1 OR 2
          DO NS=1,NSCM
            PERSED(NS,K,L) = PNEW(CORE,K,NS)/100.0_8         ! *** Set mass fraction to measured value
          ENDDO
          DTOTAL = SUM(PERSED(:,K,L))
          IF( DTOTAL > 0. )THEN
            PERSED(1:NSCM,K,L) = PERSED(1:NSCM,K,L)/DTOTAL   ! *** Ensure precision for mass balance
          ELSE
            PERSED(1:NSCM,K,L) = 0.0
          ENDIF

          ! *** Compute bed porosity and dry bulk density, depending on whether the SEDZLJ input files have wet or dry density
          IF( .true. )THEN
            ! *** BDEN is dry bulk density in g/cm3
            PORBED(L,K) = 1.-BDEN(CORE,K)/SEDDENS(CORE)
            BULKDENS(K,L) = BDEN(CORE,K)                                             ! *** Dry Bulk Density (BULKDENS)
          ELSE
            ! *** BDEN is wet bulk density in g/cm3
            PORBED(L,K) = ( BDEN(1,K)-SEDDENS(CORE) ) / ( 1.-SEDDENS(CORE) )
            BULKDENS(K,L) = (1.-PORBED(L,K))*BDEN(CORE,K)                            ! *** Dry Bulk Density (BULKDENS)
          ENDIF
          IF( BULKDENS(K,L) <= 0. )THEN
            PRINT '(" INVALID BULK DENSITY FOR CORE",I5," AT LAYER = ",I3)',NCORENO(I,J),K
            CALL STOPP('.')
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO  
    
  ! *** Disable Non-Cohesives.  Handled by SEDZLJ and CALTRAN
  SNDBT = 0.0
  ISTRAN(7) = 0
  
  ! *** ***********************************************************************
  ! *** Set Initial Layer Mass (TSED) and Active Layer Flag
  IF( NSEDFLUME /= 3 )THEN
    FORALL(L=2:LA)
      WHERE(TSED0S(1:KB,NCORENO(IL(L),JL(L))) > 0.0 )
        LAYERACTIVE(1:KB,L) = 2                                    ! *** Flag original in-place sediment layers
      ELSEWHERE
        LAYERACTIVE(1:KB,L) = 0
      ENDWHERE
      FORALL(K=1:KB)
        TSED(K,L)  = TSED0S(K,NCORENO(IL(L),JL(L)))*BULKDENS(K,L)  ! *** TSED  in units of (g/cm^2).
        TSED0(K,L) = TSED0S(K,NCORENO(IL(L),JL(L)))*BULKDENS(K,L)  ! *** TSED0 in units of (g/cm^2).
      ENDFORALL
    ENDFORALL
  
    ! *** POST-PROCESS THE DATA TO ENSURE VALID ASSIGNMENTS
    DO L=2,LA
      DO K=1,KB
        IF( TSED0(K,L)/BULKDENS(K,L) < 1.E-8 .OR. K <= 2 )THEN
          HBED(L,K) = 0.0
          TSED(K,L)  = 0.0
          TSED0(K,L) = 0.0
          TAUCOR(K,L) = 1000.
          IF( NSEDFLUME == 1 )THEN
            DO M=1,ITBM
              ERATE(K,L,M) = 0.
            ENDDO
          ENDIF
          DO NS=1,NSCM
            PERSED(NS,K,L) = 1./FLOAT(NSCM)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    ! *** GET HARD BOTTOM ELEVATIONS AND TOTAL SEDIMENT THICKNESS
    DO L=2,LA
      TSET0T(L)=0.0
      HBEDA(L)=0.0
      KBT(L)=-1                          ! *** Initialize KBT to the first layer with mass
      DO K=1,KB                          ! *** Topdown layer loop
        IF( TSED0(K,L) > 0. .AND. KBT(L) == -1 .AND. K > 2 ) KBT(L) = K
        TSET0T(L) = TSET0T(L) + TSED0(K,L)/BULKDENS(K,L)
        HBED(L,K) = 0.01*TSED(K,L)/BULKDENS(K,L)
        HBEDA(L)  = HBEDA(L) + HBED(L,K)
      ENDDO
      IF( KBT(L) == -1 ) KBT(L) = KB
      ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** Hard Bottom Elevation
    ENDDO
  ENDIF
  
  ! *** Back to BED.SDF
  ! *** Read in Newly Deposited Sediments
  ! ***    Erosion Rates (ERATEND) and Critical Shear Stress for Erosion (TAUCRITE).
  if( process_id == master_id )then
    READ (30,'(A80)') STR_LINE
    READ (30,*)  (SCND(NSC),NSC=1,NSICM)
    READ (30,'(A80)') STR_LINE
    READ (30,*)  (TAUCRITE(NSC),NSC=1,NSICM)
    READ (30,'(A80)') STR_LINE
    IF( NSEDFLUME == 1 )THEN
      DO NSC=1,NSICM
          READ(30,*)(ERATEND(NSC,M),M=1,ITBM)
      ENDDO
    ELSE
      DO NS=1,NSICM
        READ(30,*) ACTDEPA(NS), ACTDEPN(NS), ACTDEPMAX(NS)
      ENDDO
    ENDIF
  end if !***End calculation on master process

  Call Broadcast_Array(SCND,        master_id)
  Call Broadcast_Array(TAUCRITE,    master_id)
  IF( NSEDFLUME == 1 )THEN
    Call Broadcast_Array(ERATEND,   master_id)
  ELSE
    Call Broadcast_Array(ACTDEPA,   master_id)
    Call Broadcast_Array(ACTDEPN,   master_id)
    Call Broadcast_Array(ACTDEPMAX, master_id)
  ENDIF
  
  DO K=1,NSCM
    DISTAR(K) = D50(K)/10000.0*(((SEDDENS(1)/WATERDENS)-1.0)*980.0/0.01**2)**(1.0/3.0)

    ! Settling speed (DWS) calculated from D50 (micron) if DWS is < 0
    IF( DWSIN(K) == -1 )THEN
      ! input diameter using Cheng's model (1998).  Settling speed in cm/s
      DWS(K) = 0.01/(D50(K)*0.0001)*(SQRT(25.0+1.2*DISTAR(K)**2)-5.0)**1.5
    !ELSEIF( DWSIN(K) == -2 )THEN
      ! TBD
    ELSE
       DWS(K) = DWSIN(K)
    ENDIF
  ENDDO
  
  IF( IHTSTRT > 0  )THEN
    if( process_id == master_id )then
      WRITE(*,'(A)')'READING SEDBED_HOT.SDF'
      OPEN (514,FILE='SEDBED_HOT.SDF',FORM='FORMATTED',STATUS='old')
      
      READ (514,34569) ((LAYERACTIVE_Global(K,LG),K=1,KB),LG=2,LA_Global)                ! *** LAYERACTIVE
      READ (514,34569) (KBT_Global(LG),LG=2,LA_Global)                                   ! *** KBT
      READ (514,34567) (D50AVG_Global(LG),LG=2,LA_Global)                                ! *** D50AVG
      READ (514,34568) ((BULKDENS_Global(K,LG),K=1,KB),LG=2,LA_Global)                   ! *** BULKDENS
      READ (514,34568) ((TSED_Global(K,LG),K=1,KB),LG=2,LA_Global)                       ! *** TSED
      READ (514,34568) (((PERSED_Global(NS,K,LG),NS=1,NSCM),K=1,KB),LG=2,LA_Global)      ! *** PERSED
    end if !***End calculation on master process
    
    Call Broadcast_Array(LAYERACTIVE_Global, master_id)
    Call Broadcast_Array(KBT_Global, master_id)
    Call Broadcast_Array(D50AVG_Global, master_id)
    Call Broadcast_Array(BULKDENS_Global, master_id)
    Call Broadcast_Array(TSED_Global, master_id)
    Call Broadcast_Array(PERSED_Global, master_id)

    ! *** Map to local domain
    DO LG=2,LA_Global
      LL = Map2Local(LG).LL
      IF( LL > 0 )THEN
        LAYERACTIVE(:,LL) = LAYERACTIVE_Global(:,LG)
        KBT(LL)           = KBT_Global(LG)
        D50AVG(LL)        = D50AVG_Global(LG)
        BULKDENS(:,LL)    = BULKDENS_Global(:,LG)
        TSED(:,LL)        = TSED_Global(:,LG)
        PERSED(:,:,LL)    = PERSED_Global(:,:,LG)
        DO K=1,KB
          DTOTAL = SUM(PERSED(:,K,LL))
          IF( DTOTAL > 0. )THEN
            PERSED(1:NSCM,K,LL) = PERSED(1:NSCM,K,LL)/DTOTAL                             ! *** Ensure precision for mass balance
          ELSE
            PERSED(1:NSCM,K,LL) = 0.0
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    
    IF( ICALC_BL > 0 .AND. Restart_In_Ver > 1000 )THEN
      if( process_id == master_id )then
        READ (514,34568) ((CBL_Global(L,NS),NS=1,NSCM),L=2,LA_Global)                    ! *** CBL
      end if !***End calculation on master process

      Call Broadcast_Array(CBL_Global, master_id)
      DO LG=2,LA_Global
        LL = Map2Local(LG).LL
        IF( LL > 0 )THEN
          CBL(LL,1:NSCM) = CBL_Global(LG,1:NSCM)
        ENDIF
      ENDDO

      IF( ISTRAN(5) > 0 )THEN
        if( process_id == master_id )then
          READ (514,34568) ((CBLTOX_Global(L,NT),NT=1,NTOX),L=2,LA_Global)               ! *** CBLTOX
        end if !***End calculation on master process
        
        Call Broadcast_Array(CBLTOX_Global, master_id)
        DO LG=2,LA_Global
          LL = Map2Local(LG).LL
          IF( LL > 0 )THEN
            CBLTOX(LL,1:NTOX) = CBLTOX_Global(LG,1:NTOX)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    
    DO L=2,LA
      IF( LBED(L) )THEN
        LAYERACTIVE(:,L) = 0
        TSED(:,L)        = 0.0
        TSED0(:,L)       = 0.0
        PORBED(L,:)      = BEDPORC
      ELSE
        CORE = NCORENO(IL(L),JL(L))
        IF( CORE < 1 ) CORE = 1
        DO K=1,KB
          TSED0(K,L)  = TSED(K,L)
          PORBED(L,K) = 1. - BULKDENS(K,L)/SEDDENS(CORE)
        END DO
      ENDIF
    END DO    
  
    if( process_id == master_id )then
      CLOSE(514)
    end if !***End calculation on master process

    
34567 FORMAT(E17.9)
34568 FORMAT(6E17.9)
34569 FORMAT(8I8)
  ENDIF
  
  ! *** ***********************************************************************
  
  ! Read in Wave Fetch or STWAVE Data if Used (DELME - TODO FOR MPI)
  IF( ISWNWAVE == 1 )THEN
    WRITE(*,'(A)')'READING SEDZLJ: FETCH.INP'
    OPEN(UNIT=50,FILE='fetch.inp')
    DO L=2,LA
      READ (50,*) I,J,(FWDIR(LIJ(I,J),FDIR),FDIR=1,8)
    ENDDO
    CLOSE(50)
     
  ELSEIF (ISWNWAVE == 2 )THEN
    WRITE(*,'(A)')'READING SEDZLJ: STWAVE.INP'
    OPEN(UNIT=51,FILE='stwave.inp')
    CALL SKIPCOM(51,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES 
     
    READ(51,*) STWVNUM,STWVTIM
     
    STWVTIM=STWVTIM*DT/3600.0
     
    DO NWV=1,STWVNUM  
      READ (51,'(A80)') STR_LINE
      !READ(51,*)IWV,JWV,STWVHT(2,NWV),STWVTP(2,NWV),STWVDR(2,NWV)
        
      DO L=2,LA
        READ(51,*)IWV,JWV,STWVHTMP,STWVTTMP,STWVDTMP
        STWVHT(LIJ(IWV,JWV),NWV)=STWVHTMP
        STWVTP(LIJ(IWV,JWV),NWV)=STWVTTMP
        STWVDR(LIJ(IWV,JWV),NWV)=STWVDTMP
        !STWVHT(L,NWV)=STWVHT(2,NWV)
        !STWVTP(L,NWV)=STWVTP(2,NWV)
        !STWVDR(L,NWV)=STWVDR(2,NWV)     
      ENDDO
      ! Incremental counter for which wave data set we are on
      STINC=0
      NWVCOUNT=STWVTIM-1
    ENDDO
     
    CLOSE(51)
  ENDIF
  
  ! *** ***********************************************************************
  FORALL(NS=1:NSCM) SSGI(NS) = 1.0/(1.0E6*SSG(NS))  !initialize SSGI.  Specific Volume (m**3/g)
  
  if( process_id == master_id )then
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Input Sizes (micron)'
    WRITE(6,*) (SCND(K),K=1,NSICM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Input Critical Shear Ero dynes/cm^2'
    WRITE(6,*) (TAUCRITE(K),K=1,NSICM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Sediment Sizes (micron)'
    WRITE(6,*) (D50(K),K=1,NSCM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'DISTAR for each size class'
    WRITE(6,*) (DISTAR(K),K=1,NSCM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Critical Shear Sus dynes/cm^2 '
    WRITE(6,*) (TCRSUS(K),K=1,NSCM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Critical Shear Ero dynes/cm^2 '
    WRITE(6,*) (TCRE(K),K=1,NSCM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Settling Speeds in cm/s'
    WRITE(6,*) (DWS(K),K=1,NSCM) 
    WRITE(6,*) '**************************************'
    WRITE(6,*) 'Surface Sediment Density  (g/cm^3) '
    KTOP = KB
    DO K=1,KB
      IF( TSED(K,LSED(1)) > 0.0 )THEN
        KTOP = K
        EXIT
      ENDIF
    ENDDO
    WRITE(6,*) BULKDENS(KTOP,LSED(1))
    WRITE(6,*) '**************************************'
  
    ! *** EFDCLOG.OUT
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Input Sizes (micron)'
    WRITE(8,*) (SCND(K),K=1,NSICM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Input Critical Shear Ero dynes/cm^2'
    WRITE(8,*) (TAUCRITE(K),K=1,NSICM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Sediment Sizes (micron)'
    WRITE(8,*) (D50(K),K=1,NSCM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'DISTAR for each size class'
    WRITE(8,*) (DISTAR(K),K=1,NSCM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Critical Shear Sus dynes/cm^2 '
    WRITE(8,*) (TCRSUS(K),K=1,NSCM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Critical Shear Ero dynes/cm^2 '
    WRITE(8,*) (TCRE(K),K=1,NSCM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Settling Speeds in cm/s'
    WRITE(8,*) (DWS(K),K=1,NSCM) 
    WRITE(8,*) '**************************************'
    WRITE(8,*) 'Surface Sediment Density  (g/cm^3) '
    WRITE(8,*) BULKDENS(KTOP,LSED(1))
    WRITE(8,*) '**************************************'

    CLOSE(10)
    CLOSE(30)
    CLOSE(40)
  end if !***End calculation on master process
  
  ! *** ADD SMALL OFFSET FOR PRECISION ISSUES
  SCND(1)     = SCND(1)     - 1E-6
  SCND(NSICM) = SCND(NSICM) + 1E-6

  DO L=2,LA
    SH_SCALE(L) = 1.0
    IF( IHTSTRT == 0  )THEN
      DO K=1,KB
        IF( HBED(L,K) > 0.0 )THEN
          D50AVG(L) = SUM(PERSED(1:NSCM,K,L)*D50(1:NSCM))     ! *** Calculate local d50 at sediment bed surface
          EXIT
        ENDIF
      ENDDO
      D50AVG(L) = MAX(D50AVG(L),D50(1))
    ENDIF
  ENDDO

  ! *** WSEDO - Fixed/Specified settling velocity.  Convert from cm/s (DWS) to m/s (WSEDO)
  WSEDO(1:NSCM) = DWS(1:NSCM)/100.0_8
  
  ! *** Bedload size cutoff for bedload (e.g. > 64) and approach for probability of deposition
  ! *** approach (Gessler > or Krone <)
  IF( BEDLOAD_CUTOFF < 10. ) BEDLOAD_CUTOFF = 64.
  
  ! *** INITIALIZE BED ARRAYS FOR USE BY SSEDTOX AND TOXICS SUBROUTINES CALTOX AND CALTOXB
  ! *** SEDZLJ DOES NOT HAVE CONSOLIDATION.  SO PORBED, PORBED1, VDRBED ARE TEMPORALLY CONSTANT
  DO L=2,LA
    DO K=1,KB
      IF( BULKDENS(K,L) > 0.0 )THEN
        HBED(L,K) = 0.01*TSED(K,L)/BULKDENS(K,L)                                      ! *** Sediment bed layer thickness (m)
      ELSEIF( TSED(K,L) > 1E-6 .AND. BDEN(1,KB) > 0 )THEN
        HBED(L,K) = 0.01*TSED(K,L)/BDEN(1,KB)                                         ! *** Sediment bed layer thickness (m)
      ELSE
        HBED(L,K) = 0.0                                                               ! *** Sediment bed layer thickness (m)
      ENDIF
    ENDDO
    FORALL(K=1:KB) VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))                        ! *** Sediment bed void ratio. (dimensionless)  

    FORALL(K=1:KB) SEDBT(L,K) = TSED(K,L)*10000.                                      ! *** Total sediment mass (g/m^2) in a layer, TSED-sediment layer unit mass (g/cm^2)
    FORALL(K=1:KB) SEDB(L,K,1:NSCM)  = SEDBT(L,K)*PERSED(1:NSCM,K,L)                  ! *** Sediment mass (g/m^2) by class in each layer
    FORALL(K=1:KB) SEDDIA50(L,K) = SUM(PERSED(1:NSCM,K,L)*D50(1:NSCM))                ! *** D50 for sediment layer.  

    FORALL(K=1:KB) HBED1(L,K)   = HBED(L,K)
    FORALL(K=1:KB) PORBED1(L,K) = PORBED(L,K)
    FORALL(K=1:KB) VDRBED1(L,K) = VDRBED(L,K)
    FORALL(K=1:KB) SEDB1(L,K,1:NSCM) = SEDB(L,K,1:NSCM)
    
    ! *** Collapse empty layers to top of existing sediment bed
    IF( IHTSTRT > 0 )THEN
      TSUM = HBED(L,1) + HBED(L,2)
      IF( TSUM > 0. )THEN
        ! *** Active and/or Deposition layers exist.  Accumulate active and Deposition layers into one layer
        ! *** and collapse any empty layers between.
        SURFACE = -1
        DO K=3,KB
          IF( TSED(K,L) > 0. )THEN   
            SURFACE = K
            KT = K-1
            SEDB(L,KBT(L),1:NSCM) = SEDB(L,1,1:NSCM) + SEDB(L,2,1:NSCM)
            SEDBT(L,KBT(L))       = SEDBT(L,1)       + SEDBT(L,2) 
            HBED(L,KBT(L))        = HBED(L,1)        + HBED(L,2) 
            EXIT
          ENDIF
        ENDDO
        IF( SURFACE == -1 )THEN
          ! *** All parent and deep layers are missing.  Only Layers 1 and 2 are active
          KBT(L) = KB
          HBED(L,KBT(L))        = HBED(L,1)        + HBED(L,2) 
          SEDB(L,KBT(L),1:NSCM) = SEDB(L,1,1:NSCM) + SEDB(L,2,1:NSCM)
          SEDBT(L,KBT(L))       = SEDBT(L,1)       + SEDBT(L,2) 
        ENDIF
    
        ! *** Zero any layers above KBT
        KT = MIN(KBT(L)-1,1)
        DO K=KT,1,-1
          HBED(L,K)         = 0.
          SEDB(L,K,1:NSCM)  = 0.
          SEDB1(L,K,1:NSCM) = 0.
          SEDBT(L,K)        = 0.
        ENDDO
      ENDIF
    ENDIF
  
  ENDDO
  SNDVDRD = BEDPORC/(1.-BEDPORC)                                                     ! *** Non-Cohesive settling void ratio (not used in SEDZLJ)
  
  DEALLOCATE(TAUTEMP,BDEN,PNEW)

  ! *** SET THE WC TRANSPORT CUTOFF FOR NON-COHESIVES
  DO I = 1,NACTIVEWC
    IF( IACTIVEWC1(I) == 6 )THEN
      ! *** FOUND THE SEDIMENT SECTION
      DO NSC=1,NSCM
        IF( D50(NSC) > BEDLOAD_CUTOFF )THEN
          IF( PROCESS_ID == 0 ) PRINT *,'WC LIMIT SET FOR ', WCV(I+NSC-1).ID
          WCV(I+NSC-1).WCLIMIT = 1.0E-3
        ENDIF
      ENDDO
      EXIT
    ENDIF
  ENDDO
  
  RETURN

  END SUBROUTINE SEDIC

