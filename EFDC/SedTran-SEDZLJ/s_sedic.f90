! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEDIC
  
  use GLOBAL
  use INFOMOD,only: SKIPCOM, READSTR, NUMCOL

  use Variables_MPI
  use Broadcast_Routines
  use Variables_MPI_Write_Out

  implicit none
  
  integer :: CORE, INCORE, I, J, K, KT, L, LG, LL, M, NS, NT, VAR_BED, FDIR, NWV, NSC, ILocal, JLocal, SURFACE 
  integer :: IWV, JWV, KTOP, NSKIP, NCOL, VER
  character(LEN = 80)  :: STR_LINE
  character(LEN = 120) :: STR_120
  
  real(RKD) :: TBEGINSEDZLJ
  real(RKD) :: DTOTAL, STWVHTMP, STWVTTMP, STWVDTMP, TSUM
  real(RKD),allocatable,dimension(:,:)   :: BDEN          ! *** (INCORE,KB)
  real(RKD),allocatable,dimension(:,:)   :: TAUTEMP       ! *** (KB)
  real(RKD),allocatable,dimension(:,:,:) :: PNEW          ! *** (INCORE,KB,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: TSED0S        ! *** (KB,INCORE)
  real(RKD),allocatable,dimension(:,:,:) :: ERATETEMP     ! *** (INCORE,KB,ITBM)

  character(100) :: STR                    !< String buffer for input

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
 
  if( process_id == master_id )then
    write(*,'(A)')'READING SEDZLJ FILES'
  
    !CALL SEDDATA !calls routine to convert SEDflume data into a form useable by this code
    open(UNIT = 10,FILE = 'erate.sdf')
    open(UNIT = 30,FILE = 'bed.sdf')
  
    ! Read in Sediment Transport Parameters
    ! VAR_BED  = 1 for variable sediment bed
    ! ICALC_BL > 0 for bedload  calculation
    read(30,'(A80)') STR_LINE
    read(30,*) VAR_BED, KB, ICALC_BL, SEDSTEP, SEDSTART, IHTSTRT, IMORPH, ISWNWAVE, MAXDEPLIMIT, HPMIN
    if( HPMIN < 0.003 .or. HPMIN >= 1.0 ) HPMIN = 0.25
  
    if( SEDSTEP < TIDALP/REAL(NTSPTC) ) SEDSTEP = TIDALP/REAL(NTSPTC)
    DTSED = SEDSTEP
    DTSEDJ = REAL(DTSED,8)
  
    if( SEDSTART <= TBEGIN ) SEDSTART = TBEGIN
  
    ! *** NSEDS = Maximum number of grainsize classes defined. 
    ! *** ITBM  = Number of Sedflume Shear Categories
    ! *** NSICM = Maximum number of redeposited grainsize classes defined.  NSICM = NSEDS most times.
    read(30,'(A80)') STR_LINE 
    read(30,'(A80)') STR_LINE    ! DATA READ BY SCANSEDZLJ; ITBM and NSICM.
  
    read(30,'(A80)') STR_LINE
    read(30,*)  ZBSKIN,TAUCONST,ISSLOPE,BEDLOAD_CUTOFF

    ! *** Median grain size for each size class
    read(30,'(A80)') STR_LINE
    read(30,*) (D50(K),K = 1,NSEDS)  

    ! *** Erosion
    read(30,'(A80)') STR_LINE
    read(30,*) (TCRE(K),K = 1,NSEDS) 

    ! *** Suspension
    read(30,'(A80)') STR_LINE
    read(30,*) (TCRSUS(K),K = 1,NSEDS)

    ! *** Settling Velocities
    read(30,'(A80)') STR_LINE
    read(30,*) (DWSIN(K),K = 1,NSEDS)

    ! Hydrophobic Contaminant Information (PMC - 2016-11-09 - Deprecated - Handled by CALTOX/CALTOXB)
    !IF( NSEDFLUME == 2 )then
      !READ (30,'(A80)') STR_LINE
      !READ (30,*) (KPART(K),K = 1,NSEDS)
      !READ (30,'(A80)') STR_LINE
      !READ (30,*) (DIFFCOFF(K),K = 1,NSEDS)
      !READ (30,'(A80)') STR_LINE
      !DO LL = 1,KB
      !   read(30,*) (PCONTEMP(K,LL),K = 1,NSEDS)
      !ENDDO
    !ENDIF
  
    ! *** ***********************************************************************
    !Reading in Erate.sdf starting with the layer's thickness.
    read(10,'(A80)') STR_LINE
    read(10,*) TACTM !read in active layer multiplier
    ! Read in Initial Erosion Data
  endif !***End calculation on master process

  call Broadcast_Scalar(VAR_BED,  master_id)
  call Broadcast_Scalar(KB,       master_id)
  call Broadcast_Scalar(ICALC_BL, master_id)
  call Broadcast_Scalar(SEDSTEP,  master_id)
  call Broadcast_Scalar(SEDSTART, master_id)
  call Broadcast_Scalar(IHTSTRT,  master_id)
  call Broadcast_Scalar(IMORPH,   master_id)
  call Broadcast_Scalar(ISWNWAVE, master_id)
  call Broadcast_Scalar(MAXDEPLIMIT, master_id)
  call Broadcast_Scalar(HPMIN,    master_id)
  call Broadcast_Scalar(DTSED,    master_id)
  call Broadcast_Scalar(DTSEDJ,   master_id)

  call Broadcast_Scalar(ZBSKIN,   master_id)
  call Broadcast_Scalar(TAUCONST, master_id)
  call Broadcast_Scalar(ISSLOPE,  master_id)
  call Broadcast_Scalar(BEDLOAD_CUTOFF, master_id)
  
  call Broadcast_Scalar(TACTM, master_id)

  call Broadcast_Array(D50, master_id)
  call Broadcast_Array(TCRE, master_id)
  call Broadcast_Array(TCRSUS, master_id)
  call Broadcast_Array(DWSIN, master_id)

  if( VAR_BED >= 1 )then    
    ! Variable Bed *************************************************
    allocate(I2D_Global(IC_Global,JC_Global))
    I2D_Global = 0
    
    if( process_id == master_id )then
      open(UNIT = 20,FILE = 'core_field.sdf')

      ! Determine File Format
      read(20,'(A120)')STR_120
      I = 1
      do while ( i < 117 )
        if( STR_120(I:I+2) == 'DSI' .or. STR_120(I:I+2) == 'dsi' )EXIT
        I = I + 1
      enddo
     
      close(20)
      
      ! *** READ THE DATA
      open(UNIT = 20,FILE = 'core_field.sdf')
      if( I < 117 )then
        ! *** DSI STANDARD
        read(20,'(A80)') STR_LINE
        read(20,'(A80)') STR_LINE
        read(20,'(A80)') STR_LINE
        read(20,*) INCORE !read the number of cores    
        read(20,'(A80)') STR_LINE
        read(20,'(A80)') STR_LINE
        do LG = 2,LA_Global
          read(20,*) I, J, CORE
          I2D_Global(I,J) = CORE                                ! *** NCORENO
        enddo
      else
        ! *** SNL STANDARD
        read(20,*)INCORE !read the number of cores    
        do J = JC,1,-1
          read(20,'(120(I1,1X))')(I2D_Global(I,J),I = 1,IC)       ! *** NCORENO
        enddo
      endif
      close(20)
    endif !***End calculation on master process
    call Broadcast_Scalar(INCORE,    master_id)
    call Broadcast_Array(I2D_Global, master_id)
    
    ! *** Map to local domain
    do J = 3,JC_GLOBAL-2
      do I = 3,IC_GLOBAL-2
        LG = LIJ_Global(I,J)
        if( LG > 0 )then
          L = Map2Local(LG).LL
          if( L > 0 )then
            ILocal = Map2Local(LG).IL
            JLocal = Map2Local(LG).JL
            NCORENO(ILocal,JLocal) = I2D_Global(I,J)
          endif
        endif
      enddo
    enddo
    deallocate(I2D_Global)
  else
    ! Constant Erosion in Horizontal ********************************  
    INCORE = 1

    do L = 2,LA
      I = IL(L)
      J = JL(L)
      NCORENO(I,J) = 1
    enddo
  endif
 
  ! *** ALLOCATIONS AND INITIALIZATIONS
  if( NSEDFLUME == 1 )then
    allocate(ERATETEMP(INCORE,KB,ITBM))
    ERATETEMP = 0.0
  else
    ITBM = 2
    allocate(EA(INCORE,KBM))
    allocate(EN(INCORE,KBM))
    allocate(MAXRATE(INCORE,KBM))
    EA = 0.0
    EN = 0.0
    MAXRATE = 10000.0
  endif

  allocate(ACTDEPA(NSICM))
  allocate(ACTDEPN(NSICM))
  allocate(ACTDEPMAX(NSICM))
  allocate(BDEN(INCORE,KBM))
  allocate(ERATEND(NSICM,ITBM))
  allocate(ERATE(KB,LCM,ITBM))
  allocate(PNEW(INCORE,KBM,NSEDS+1))
  allocate(SEDDENS(INCORE))
  allocate(TAUTEMP(INCORE,KBM))
  allocate(TSED0S(KB,INCORE))
  allocate(TAULOC(ITBM))

  ACTDEPA = 0.0
  ACTDEPN = 0.0
  ACTDEPMAX = 10000.0
  BDEN = 0.0   
  ERATEND = 0.0 
  ERATE = 0.0   
  PNEW = 0.0
  SEDDENS = 0.0 
  TAUTEMP = 0.0
  TSED0S = 0.0
  TAULOC = 0.0    

  if( process_id == master_id )then
    ! *** **********************************************************   
    do CORE = 1,INCORE  ! for each core of data      
      read(10,'(A80)') STR_LINE
      read(10,*)(TAUTEMP(CORE,K),K = 1,KB)           ! *** read the critical shear stresses of the core
      read(10,'(A80)') STR_LINE
      read(10,*) (TSED0S(K,CORE),K = 1,KB)          ! *** read in layer thickness 
      read(10,'(A80)') STR_LINE      
      read(10,*) (BDEN(CORE,K),K = 1,KB)             ! *** read in the bulk density of the core 
      read(10,'(A80)') STR_LINE      
      read(10,*) WATERDENS, SEDDENS(CORE)          ! *** read in the water density and sediment solid's density
      read(10,'(A80)') STR_LINE  
      do K = 1,KB
        read(10,*)(PNEW(CORE,K,NS),NS = 1,NSEDS)
      enddo       
      read(10,'(A80)') STR_LINE
      if( NSEDFLUME == 1 )then
        do M = 1,ITBM
          read(10,*)TAULOC(M)                      ! *** shear stress used to erode a portion of the core
          read(10,*)(ERATETEMP(CORE,K,M),K = 1,KB)   ! *** erosion rate for each layer subject to shear stress TAULOC
        enddo
      else
        TAULOC(1) = 0.0
        TAULOC(2) = 1000.                          ! *** Set to a very large maximum shear stress
        do K = 1,KB
          read(10,*) EA(CORE,K), EN(CORE,K), MAXRATE(CORE,K)
        enddo
      endif
    enddo     
  endif !***End calculation on master process

  call Broadcast_Scalar(WATERDENS,  master_id)

  call Broadcast_Array(TAUTEMP,     master_id)
  call Broadcast_Array(TSED0S,      master_id)
  call Broadcast_Array(BDEN,        master_id)
  call Broadcast_Array(SEDDENS,     master_id)
  call Broadcast_Array(PNEW,        master_id)
  if( NSEDFLUME == 1 )then
    call Broadcast_Array(TAULOC,    master_id)
    call Broadcast_Array(ERATETEMP, master_id)
  else
    call Broadcast_Array(TAULOC,    master_id)
    call Broadcast_Array(EA,        master_id)
    call Broadcast_Array(EN,        master_id)
    call Broadcast_Array(MAXRATE,   master_id)
      
    do K = 1,KB
      do CORE = 1,INCORE  ! for each core of data      
        MAXRATE(CORE,K) = MAXRATE(CORE,K)*BDEN(CORE,K)                             ! *** Convert from cm/s to g/m2/s using dry bulk density
      enddo
    enddo
  endif
  
  
  ! ***
  do L = 2,LA
    I = IL(L)   ! *** I location as a function of L
    J = JL(L)   ! *** J location as a function of L
    
    if( NCORENO(I,J) > 0 )then
      CORE = NCORENO(I,J)
      do K = 1,KB
        TAUCOR(K,L) = TAUTEMP(CORE,K)                    ! *** Critical shear stresses from cores
        if( NSEDFLUME == 1 )then
          do M = 1,ITBM
            ERATE(K,L,M) = ERATETEMP(CORE,K,M)           ! *** Set erosion rate to measured value
          enddo
        endif
        
        if( NSEDFLUME == 3 )then
          
        else
          ! *** NSEDFLUME = 1 OR 2
          do NS = 1,NSEDS
            PERSED(NS,K,L) = PNEW(CORE,K,NS)/100.0_8         ! *** Set mass fraction to measured value
          enddo
          DTOTAL = SUM(PERSED(:,K,L))
          if( DTOTAL > 0. )then
            PERSED(1:NSEDS,K,L) = PERSED(1:NSEDS,K,L)/DTOTAL   ! *** Ensure precision for mass balance
          else
            PERSED(1:NSEDS,K,L) = 0.0
          endif

          ! *** Compute bed porosity and dry bulk density, depending on whether the SEDZLJ input files have wet or dry density
          if( .true. )then
            ! *** BDEN is dry bulk density in g/cm3
            PORBED(L,K) = 1.-BDEN(CORE,K)/SEDDENS(CORE)
            BULKDENS(K,L) = BDEN(CORE,K)                                             ! *** Dry Bulk Density (BULKDENS)
          else
            ! *** BDEN is wet bulk density in g/cm3
            PORBED(L,K) = ( BDEN(1,K)-SEDDENS(CORE) ) / ( 1.-SEDDENS(CORE) )
            BULKDENS(K,L) = (1.-PORBED(L,K))*BDEN(CORE,K)                            ! *** Dry Bulk Density (BULKDENS)
          endif
          if( BULKDENS(K,L) <= 0. )then
            PRINT '(" INVALID BULK DENSITY FOR CORE",I5," AT LAYER = ",I3)',NCORENO(I,J),K
            call STOPP('', 1)
          endif
        endif
      enddo
    endif
  enddo  
    
  ! *** Disable Non-Cohesives.  Handled by SEDZLJ and CALTRAN
  SNDBT = 0.0
  ISTRAN(7) = 0
  
  ! *** ***********************************************************************
  ! *** Set Initial Layer Mass (TSED) and Active Layer Flag
  if( NSEDFLUME /= 3 )then
    FORALL(L = 2:LA)
      WHERE(TSED0S(1:KB,NCORENO(IL(L),JL(L))) > 0.0 )
        LAYERACTIVE(1:KB,L) = 2                                    ! *** Flag original in-place sediment layers
      elseWHERE
        LAYERACTIVE(1:KB,L) = 0
      ENDWHERE
      FORALL(K = 1:KB)
        TSED(K,L)  = TSED0S(K,NCORENO(IL(L),JL(L)))*BULKDENS(K,L)  ! *** TSED  in units of (g/cm^2).
        TSED0(K,L) = TSED0S(K,NCORENO(IL(L),JL(L)))*BULKDENS(K,L)  ! *** TSED0 in units of (g/cm^2).
      ENDFORALL
    ENDFORALL
  
    ! *** POST-PROCESS THE DATA TO ENSURE VALID ASSIGNMENTS
    do L = 2,LA
      do K = 1,KB
        if( TSED0(K,L)/BULKDENS(K,L) < 1.E-8 .or. K <= 2 )then
          HBED(L,K) = 0.0
          TSED(K,L)  = 0.0
          TSED0(K,L) = 0.0
          TAUCOR(K,L) = 1000.
          if( NSEDFLUME == 1 )then
            do M = 1,ITBM
              ERATE(K,L,M) = 0.
            enddo
          endif
          do NS = 1,NSEDS
            PERSED(NS,K,L) = 1./FLOAT(NSEDS)
          enddo
        endif
      enddo
    enddo

    ! *** GET HARD BOTTOM ELEVATIONS AND TOTAL SEDIMENT THICKNESS
    do L = 2,LA
      TSET0T(L) = 0.0
      HBEDA(L) = 0.0
      KBT(L) = -1                          ! *** Initialize KBT to the first layer with mass
      do K = 1,KB                          ! *** Topdown layer loop
        if( TSED0(K,L) > 0. .and. KBT(L) == -1 .and. K > 2 ) KBT(L) = K
        TSET0T(L) = TSET0T(L) + TSED0(K,L)/BULKDENS(K,L)
        HBED(L,K) = 0.01*TSED(K,L)/BULKDENS(K,L)
        HBEDA(L)  = HBEDA(L) + HBED(L,K)
      enddo
      if( KBT(L) == -1 ) KBT(L) = KB
      ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** Hard Bottom Elevation
    enddo
  endif
  
  ! *** Back to BED.SDF
  ! *** Read in Newly Deposited Sediments
  ! ***    Erosion Rates (ERATEND) and Critical Shear Stress for Erosion (TAUCRITE).
  if( process_id == master_id )then
    read(30,'(A80)') STR_LINE
    read(30,*)  (SCND(NSC),NSC = 1,NSICM)
    read(30,'(A80)') STR_LINE
    read(30,*)  (TAUCRITE(NSC),NSC = 1,NSICM)
    read(30,'(A80)') STR_LINE
    if( NSEDFLUME == 1 )then
      do NSC = 1,NSICM
          read(30,*)(ERATEND(NSC,M),M = 1,ITBM)
      enddo
    else
      do NS = 1,NSICM
        read(30,*) ACTDEPA(NS), ACTDEPN(NS), ACTDEPMAX(NS)
      enddo
    endif
  endif !***End calculation on master process

  call Broadcast_Array(SCND,        master_id)
  call Broadcast_Array(TAUCRITE,    master_id)
  if( NSEDFLUME == 1 )then
    call Broadcast_Array(ERATEND,   master_id)
  else
    call Broadcast_Array(ACTDEPA,   master_id)
    call Broadcast_Array(ACTDEPN,   master_id)
    call Broadcast_Array(ACTDEPMAX, master_id)
  endif
  
  do K = 1,NSEDS
    DISTAR(K) = D50(K)/10000.0*(((SEDDENS(1)/WATERDENS)-1.0)*980.0/0.01**2)**(1.0/3.0)

    ! Settling speed (DWS) calculated from D50 (micron) if DWS is < 0
    if( DWSIN(K) == -1 )then
      ! input diameter using Cheng's model (1998).  Settling speed in cm/s
      DWS(K) = 0.01/(D50(K)*0.0001)*(SQRT(25.0+1.2*DISTAR(K)**2)-5.0)**1.5
    !ELSEIF( DWSIN(K) == -2 )then
      ! TBD
    else
       DWS(K) = DWSIN(K)
    endif
  enddo
  
  if( IHTSTRT > 0  )then
    if( process_id == master_id )then
      write(*,'(A)')'READING SEDBED_HOT.SDF'
      OPEN (514,FILE = 'SEDBED_HOT.SDF',FORM = 'FORMATTED',STATUS = 'old')
      
      ! *** Check restart file version
      read(514,'(A)') STR
      NCOL = NUMCOL(STR)
      VER = -1
      if( NCOL == 2 )then
        read(STR,*) TBEGINSEDZLJ, VER                   !< Restart time and version
      endif

      if( VER < 1240 )then
        ! *** Reset file
        close(514)
        OPEN (514,FILE = 'SEDBED_HOT.SDF',FORM = 'FORMATTED',STATUS = 'old')
      endif
      
      read(514,34569) ((LAYERACTIVE_Global(K,LG),K = 1,KB),LG = 2,LA_Global)                   ! *** LAYERACTIVE
      read(514,34569) (KBT_Global(LG),LG = 2,LA_Global)                                        ! *** KBT
      read(514,34567) (D50AVG_Global(LG),LG = 2,LA_Global)                                     ! *** D50AVG
      read(514,34568) ((BULKDENS_Global(K,LG),K = 1,KB),LG = 2,LA_Global)                      ! *** BULKDENS
      read(514,34568) ((TSED_Global(K,LG),K = 1,KB),LG = 2,LA_Global)                          ! *** TSED
      if( VER >= 1240 )then
        read(514,34568) ((TSED0_Global(K,LG),K = 1,KB),LG = 2,LA_Global)                       ! *** TSED
      endif
      read(514,34568) (((PERSED_Global(NS,K,LG),NS = 1,NSEDS),K = 1,KB),LG = 2,LA_Global)      ! *** PERSED
    endif !***End calculation on master process
    
    call Broadcast_Array(LAYERACTIVE_Global, master_id)
    call Broadcast_Array(KBT_Global,         master_id)
    call Broadcast_Array(D50AVG_Global,      master_id)
    call Broadcast_Array(BULKDENS_Global,    master_id)
    call Broadcast_Array(TSED_Global,        master_id)
    call Broadcast_Array(TSED0_Global,       master_id)
    call Broadcast_Array(PERSED_Global,      master_id)

    ! *** Map to local domain
    do LG = 2,LA_Global
      LL = Map2Local(LG).LL
      if( LL > 0 )then
        LAYERACTIVE(:,LL) = LAYERACTIVE_Global(:,LG)
        KBT(LL)           = KBT_Global(LG)
        D50AVG(LL)        = D50AVG_Global(LG)
        BULKDENS(:,LL)    = BULKDENS_Global(:,LG)
        TSED(:,LL)        = TSED_Global(:,LG)
        if( VER >= 12400 )then
          TSED0(:,LL)        = TSED0_Global(:,LG)
        endif
        PERSED(:,:,LL)    = PERSED_Global(:,:,LG)
        do K = 1,KB
          DTOTAL = SUM(PERSED(:,K,LL))
          if( DTOTAL > 0. )then
            PERSED(1:NSEDS,K,LL) = PERSED(1:NSEDS,K,LL)/DTOTAL                             ! *** Ensure precision for mass balance
          else
            PERSED(1:NSEDS,K,LL) = 0.0
          endif
        enddo
      endif
    enddo
    
    if( ICALC_BL > 0 .and. Restart_In_Ver > 1000 )then
      if( process_id == master_id )then
        read(514,34568) ((CBL_Global(L,NS),NS = 1,NSEDS),L = 2,LA_Global)                    ! *** CBL
      endif !***End calculation on master process

      call Broadcast_Array(CBL_Global, master_id)
      do LG = 2,LA_Global
        LL = Map2Local(LG).LL
        if( LL > 0 )then
          CBL(LL,1:NSEDS) = CBL_Global(LG,1:NSEDS)
        endif
      enddo

      if( ISTRAN(5) > 0 )then
        if( process_id == master_id )then
          read(514,34568) ((CBLTOX_Global(L,NT),NT = 1,NTOX),L = 2,LA_Global)               ! *** CBLTOX
        endif !***End calculation on master process
        
        call Broadcast_Array(CBLTOX_Global, master_id)
        do LG = 2,LA_Global
          LL = Map2Local(LG).LL
          if( LL > 0 )then
            CBLTOX(LL,1:NTOX) = CBLTOX_Global(LG,1:NTOX)
          endif
        enddo
      endif
    endif
    
    do L = 2,LA
      if( LBED(L) )then
        LAYERACTIVE(:,L) = 0
        TSED(:,L)        = 0.0
        TSED0(:,L)       = 0.0
        PORBED(L,:)      = BEDPORC
      else
        CORE = NCORENO(IL(L),JL(L))
        if( CORE < 1 ) CORE = 1
        do K = 1,KB
          if( VER < 12400 ) TSED0(K,L)  = TSED(K,L)
          PORBED(L,K) = 1. - BULKDENS(K,L)/SEDDENS(CORE)
        enddo
      endif
    enddo    
  
    if( process_id == master_id )then
      close(514)
    endif !***End calculation on master process

    
34567 FORMAT(E17.9)
34568 FORMAT(6E17.9)
34569 FORMAT(8I8)
  endif
  
  ! *** ***********************************************************************
  
  ! Read in Wave Fetch or STWAVE Data if Used
  if( ISWNWAVE == 1 )then
    write(*,'(A)')'READING SEDZLJ: FETCH.INP'
    open(UNIT = 50,FILE = 'fetch.inp')
    do L = 2,LA
      read(50,*) I,J,(FWDIR(LIJ(I,J),FDIR),FDIR = 1,8)
    enddo
    close(50)
     
  elseif( ISWNWAVE == 2 )then
    write(*,'(A)')'READING SEDZLJ: STWAVE.INP'
    open(UNIT = 51,FILE = 'stwave.inp')
    call SKIPCOM(51,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES 
     
    read(51,*) STWVNUM,STWVTIM
     
    STWVTIM = STWVTIM*DT/3600.0
     
    do NWV = 1,STWVNUM  
      read(51,'(A80)') STR_LINE
      !READ(51,*)IWV,JWV,STWVHT(2,NWV),STWVTP(2,NWV),STWVDR(2,NWV)
        
      do L = 2,LA
        read(51,*)IWV,JWV,STWVHTMP,STWVTTMP,STWVDTMP
        STWVHT(LIJ(IWV,JWV),NWV) = STWVHTMP
        STWVTP(LIJ(IWV,JWV),NWV) = STWVTTMP
        STWVDR(LIJ(IWV,JWV),NWV) = STWVDTMP
        !STWVHT(L,NWV) = STWVHT(2,NWV)
        !STWVTP(L,NWV) = STWVTP(2,NWV)
        !STWVDR(L,NWV) = STWVDR(2,NWV)     
      enddo
      ! Incremental counter for which wave data set we are on
      STINC = 0
      NWVCOUNT = STWVTIM-1
    enddo
     
    close(51)
  endif
  
  ! *** ***********************************************************************
  FORALL(NS = 1:NSEDS2) SSGI(NS) = 1.0/(1.0E6*SSG(NS))  !initialize SSGI.  Specific Volume (m**3/g)
  
  if( process_id == master_id )then
    write(6,*) '**************************************'
    write(6,*) 'Input Sizes (micron)'
    write(6,*) (SCND(K),K = 1,NSICM) 
    write(6,*) '**************************************'
    write(6,*) 'Input Critical Shear Ero dynes/cm^2'
    write(6,*) (TAUCRITE(K),K = 1,NSICM) 
    write(6,*) '**************************************'
    write(6,*) 'Sediment Sizes (micron)'
    write(6,*) (D50(K),K = 1,NSEDS) 
    write(6,*) '**************************************'
    write(6,*) 'DISTAR for each size class'
    write(6,*) (DISTAR(K),K = 1,NSEDS) 
    write(6,*) '**************************************'
    write(6,*) 'Critical Shear Sus dynes/cm^2 '
    write(6,*) (TCRSUS(K),K = 1,NSEDS) 
    write(6,*) '**************************************'
    write(6,*) 'Critical Shear Ero dynes/cm^2 '
    write(6,*) (TCRE(K),K = 1,NSEDS) 
    write(6,*) '**************************************'
    write(6,*) 'Settling Speeds in cm/s'
    write(6,*) (DWS(K),K = 1,NSEDS) 
    write(6,*) '**************************************'
    write(6,*) 'Surface Sediment Density  (g/cm^3) '
    KTOP = KB
    do K = 1,KB
      if( TSED(K,LSED(1)) > 0.0 )then
        KTOP = K
        EXIT
      endif
    enddo
    write(6,*) BULKDENS(KTOP,LSED(1))
    write(6,*) '**************************************'
  
    ! *** EFDCLOG.OUT
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Input Sizes (micron)'
    write(mpi_efdc_out_unit,*) (SCND(K),K = 1,NSICM) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Input Critical Shear Ero dynes/cm^2'
    write(mpi_efdc_out_unit,*) (TAUCRITE(K),K = 1,NSICM) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Sediment Sizes (micron)'
    write(mpi_efdc_out_unit,*) (D50(K),K = 1,NSEDS) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'DISTAR for each size class'
    write(mpi_efdc_out_unit,*) (DISTAR(K),K = 1,NSEDS) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Critical Shear Sus dynes/cm^2 '
    write(mpi_efdc_out_unit,*) (TCRSUS(K),K = 1,NSEDS) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Critical Shear Ero dynes/cm^2 '
    write(mpi_efdc_out_unit,*) (TCRE(K),K = 1,NSEDS) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Settling Speeds in cm/s'
    write(mpi_efdc_out_unit,*) (DWS(K),K = 1,NSEDS) 
    write(mpi_efdc_out_unit,*) '**************************************'
    write(mpi_efdc_out_unit,*) 'Surface Sediment Density  (g/cm^3) '
    write(mpi_efdc_out_unit,*) BULKDENS(KTOP,LSED(1))
    write(mpi_efdc_out_unit,*) '**************************************'

    close(10)
    close(30)
    close(40)
  endif !***End calculation on master process
  
  ! *** ADD SMALL OFFSET FOR PRECISION ISSUES
  SCND(1)     = SCND(1)     - 1E-6
  SCND(NSICM) = SCND(NSICM) + 1E-6

  do L = 2,LA
    SH_SCALE(L) = 1.0
    if( IHTSTRT == 0  )then
      do K = 1,KB
        if( HBED(L,K) > 0.0 )then
          D50AVG(L) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))     ! *** Calculate local d50 at sediment bed surface
          EXIT
        endif
      enddo
      D50AVG(L) = max(D50AVG(L),D50(1))
    endif
  enddo

  ! *** Initialize for use outside of SEDZLJ specific routines
  WSEDO(1:NSEDS) = DWS(1:NSEDS)/100.0_8      ! *** WSEDO  - Fixed/Specified settling velocity.  Convert from cm/s (DWS) to m/s (WSEDO)
  SEDDIA(1:NSEDS) = D50(1:NSEDS)/1.d6        ! *** SEDDIA - Nominal diameter for each sediment class, even cohesives.  Microns to m.
  
  ! *** Bedload size cutoff for bedload (e.g. > 64) and approach for probability of deposition
  ! *** approach (Gessler > or Krone <)
  if( BEDLOAD_CUTOFF < 10. ) BEDLOAD_CUTOFF = 64.
  
  ! *** INITIALIZE BED ARRAYS FOR use BY SSEDTOX AND TOXICS SUBROUTINES CALTOX AND CALTOXB
  ! *** SEDZLJ DOES NOT HAVE CONSOLIDATION.  SO PORBED, PORBED1, VDRBED ARE TEMPORALLY CONSTANT
  do L = 2,LA
    do K = 1,KB
      if( BULKDENS(K,L) > 0.0 )then
        HBED(L,K) = 0.01*TSED(K,L)/BULKDENS(K,L)                                      ! *** Sediment bed layer thickness (m)
      elseif( TSED(K,L) > 1E-6 .and. BDEN(1,KB) > 0 )then
        HBED(L,K) = 0.01*TSED(K,L)/BDEN(1,KB)                                         ! *** Sediment bed layer thickness (m)
      else
        HBED(L,K) = 0.0                                                               ! *** Sediment bed layer thickness (m)
      endif
    enddo
    FORALL(K = 1:KB) VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))                        ! *** Sediment bed void ratio. (dimensionless)  

    FORALL(K = 1:KB) SEDBT(L,K) = TSED(K,L)*10000.                                      ! *** Total sediment mass (g/m^2) in a layer, TSED-sediment layer unit mass (g/cm^2)
    FORALL(K = 1:KB) SEDB(L,K,1:NSEDS)  = SEDBT(L,K)*PERSED(1:NSEDS,K,L)                  ! *** Sediment mass (g/m^2) by class in each layer
    FORALL(K = 1:KB) SEDDIA50(L,K) = SUM(PERSED(1:NSEDS,K,L)*D50(1:NSEDS))                ! *** D50 for sediment layer.  

    FORALL(K = 1:KB) HBED1(L,K)   = HBED(L,K)
    FORALL(K = 1:KB) PORBED1(L,K) = PORBED(L,K)
    FORALL(K = 1:KB) VDRBED1(L,K) = VDRBED(L,K)
    FORALL(K = 1:KB) SEDB1(L,K,1:NSEDS) = SEDB(L,K,1:NSEDS)
    
    ! *** Collapse empty layers to top of existing sediment bed
    if( IHTSTRT > 0 )then
      TSUM = HBED(L,1) + HBED(L,2)
      if( TSUM > 0. )then
        ! *** Active and/or Deposition layers exist.  Accumulate active and Deposition layers into one layer
        ! *** and collapse any empty layers between.
        SURFACE = -1
        do K = 3,KB
          if( TSED(K,L) > 0. )then   
            SURFACE = K
            KT = K-1
            SEDB(L,KBT(L),1:NSEDS) = SEDB(L,1,1:NSEDS) + SEDB(L,2,1:NSEDS)
            SEDBT(L,KBT(L))        = SEDBT(L,1)        + SEDBT(L,2) 
            HBED(L,KBT(L))         = HBED(L,1)         + HBED(L,2) 
            EXIT
          endif
        enddo
        if( SURFACE == -1 )then
          ! *** All parent and deep layers are missing.  Only Layers 1 and 2 are active
          KBT(L) = KB
          HBED(L,KBT(L))         = HBED(L,1)         + HBED(L,2) 
          SEDB(L,KBT(L),1:NSEDS) = SEDB(L,1,1:NSEDS) + SEDB(L,2,1:NSEDS)
          SEDBT(L,KBT(L))        = SEDBT(L,1)        + SEDBT(L,2) 
        endif
    
        ! *** Zero any layers above KBT
        KT = max(KBT(L)-1,1)
        do K = KT,1,-1
          HBED(L,K)          = 0.
          SEDB(L,K,1:NSEDS)  = 0.
          SEDB1(L,K,1:NSEDS) = 0.
          SEDBT(L,K)         = 0.
        enddo
      endif
    endif
  
  enddo
  SNDVDRD = BEDPORC/(1.-BEDPORC)                                                     ! *** Non-Cohesive settling void ratio (not used in SEDZLJ)
  
  deallocate(TAUTEMP,BDEN,PNEW)
  
  return
  END SUBROUTINE SEDIC

