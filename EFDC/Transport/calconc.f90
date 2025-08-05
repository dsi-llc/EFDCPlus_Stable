! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALCONC

  ! *** SUBROUTINE CALCULATES THE CONCENTRATION OF DISSOLVED AND
  ! *** SUSPENDED CONSTITUENTS, INCLUDING SALINITY, TEMPERATURE, DYE,
  ! *** AND SUSPENDED SEDIMENT AT TIME LEVEL (N+1). THE VALUE OF ISTL
  ! *** INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2015-01       PAUL M. CRAIG     Added fully coupled Ice Sub-model with Frazil Ice Transport
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2014-08       D H CHUNG         SET EXPLICIT PRECISIONS OF INTEGER & REAL
  !    2011-03       PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP
  !    2002-05       John Hamrick      Modified calls to calbal and budget subroutines
  !                                     added calls to bal2t2, bal2t3
  !------------------------------------------------------------------------------------------------!

  use GLOBAL
  use Allocate_Initialize      
  use Variables_Propwash
  use turbulence, only: eps_min,k_min
  use OMP_LIB
  use HEAT_MODULE, only:CALHEAT
  
  use MPI
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out
  use Communicate_Ghost_Routines

  implicit none

  integer :: K, L, LP, NS, ND, IT, IW, I, J
  integer :: NTMP, LF, LL, LE, LN, LG, IERR
  integer, save :: NICE, IFIRST, NANTIDIFF
  integer,save,allocatable,dimension(:) :: ISKIP

  real      :: CDTMP, RCDZKMK, RCDZKK
  real,save :: SEDTIME

  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, TTDS1, TTDS2, TMPCOMM, TWAIT   ! *** Model timing temporary variables
  real(RKD)           :: CCUBTMP, CCMBTMP                     ! *** Vertical diffusion temporary variables

  ! *** VERTICAL DIFFUSION VARIABLES
  real(RKD),save,allocatable,dimension(:,:) :: CCLBTMP
  real(RKD),save,allocatable,dimension(:,:) :: EEB
  real(RKD),save,allocatable,dimension(:,:) :: VCU

  !REAL,save,allocatable,dimension(:) :: TOXASM
  !REAL,save,allocatable,dimension(:) :: SEDASM
  !REAL,save,allocatable,dimension(:) :: SNDASM
  
  ! *** GOTM arrays 
  ! *** The data type must be consistent to global TKE and EPS - DKT
  real,save,allocatable,dimension(:,:) :: TKE_TEMP
  real,save,allocatable,dimension(:,:) :: TKE_TEMP1
  real,save,allocatable,dimension(:,:) :: EPS_TEMP
  real,save,allocatable,dimension(:,:) :: EPS_TEMP1
  
  logical,save :: BYPASS
  
  if( .not. allocated(EEB) )then
    call AllocateDSI(CCLBTMP,   LCM,  KCM,  0.0)
    call AllocateDSI(EEB,       LCM,  KCM,  0.0)
    call AllocateDSI(VCU,       LCM,  KCM,  0.0)
    call AllocateDSI(TKE_TEMP,  LCM,  KCM,  0.0)
    call AllocateDSI(TKE_TEMP1, LCM,  KCM,  0.0)
    call AllocateDSI(EPS_TEMP,  LCM,  KCM,  0.0)
    call AllocateDSI(EPS_TEMP1, LCM,  KCM,  0.0)
    
    !Call AllocateDSI(TOXASM,  NTXM,     0.0)
    !Call AllocateDSI(SEDASM,  NSCM2,    0.0)
    !Call AllocateDSI(SNDASM,  NSNM2,    0.0)
    if( ISGOTM > 0 )then 
      call AllocateDSI(ISKIP, NACTIVEWC+2,    0)
    else
      call AllocateDSI(ISKIP, NACTIVEWC,    0)
    endif

    NICE = 0
    SEDTIME = 0.0
    IFIRST = 0
    
    NANTIDIFF = 0
    do IW = 1,NACTIVEWC
      if( ISADAC(IACTIVEWC1(IW)) > 0 )then 
        NANTIDIFF = NANTIDIFF + 1
      endif
    enddo
  
    ! *** Transport bypass
    BYPASS = .FALSE.
    if( ANY(WCV(:).WCLIMIT > 0.0) ) BYPASS = .TRUE.
  endif

  IT = 1

  ! ***************************************************************************************
  ! *** 3D ADVECTI0N TRANSPORT CALCULATION, STANDARD TRANSPORT
  TTDS1 = DSTIME(0)
  TWAIT = 0.
  TMPITMP = 0.

  ! *** PRESPECIFY THE UPWIND CELLS FOR 3D ADVECTION
  !$OMP PARALLEL DO PRIVATE(ND,K,LP,L,LE,LN)
  do ND = 1,NDM
    do K = 1,KC
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        if( UHDY2(L,K) >= 0.0 )then
          LUPU(L,K) = LWC(L)
        else
          LUPU(L,K) = L
        endif
        if( VHDX2(L,K) >= 0.0 )then
          LUPV(L,K) = LSC(L)
        else
          LUPV(L,K) = L
        endif
      enddo
    enddo

    if( KC > 1 )then
      do K = 1,KS
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          if( W2(L,K) >= 0. )then
            KUPW(L,K) = K
          else
            KUPW(L,K) = K + 1
          endif
        enddo
      enddo
    endif
  enddo
  !$OMP END PARALLEL DO
  
  ! *** EE7.2 - MUST ZERO ALL THREADS FOR EVERY INSTANCE DO TO POTENTIAL OF DIFFERENT THREADS MAY NOT BE ZEROED FOR CURRENT TIME
  ! *** ZERO DRY CELL FLUXES
  if( LADRY > 0 )then
    do LP = 1,LADRY
      L = LDRY(LP)
      FUHUD(L,:,:) = 0.
      FVHUD(L,:,:) = 0.
      FWUU(L,:,:) = 0.0

      LN = LNC(L)
      LE = LEC(L)
      FUHUD(LEC(L),:,:) = 0.
      FVHUD(LNC(L),:,:) = 0.
    enddo
    
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LADRY
          L = LDRY(LP)
          FUHVD(L,K,ND) = 0.0
          FVHVD(L,K,ND) = 0.0
          UUUU(L,K,ND)  = 0.0
          VVVV(L,K,ND)  = 0.0
          DUU(L,K,ND)   = 0.0
          DVV(L,K,ND)   = 0.0
          POS(L,K,ND)   = 0.0
          WWWW(L,K,ND)  = 0.0

          LN = LNC(L)
          LE = LEC(L)
          FUHVD(LE,K,ND) = 0.0
          FVHVD(LN,K,ND) = 0.0
          UUUU(LE,K,ND)  = 0.0
          VVVV(LN,K,ND)  = 0.0
          DUU(LE,K,ND)   = 0.0
          DVV(LN,K,ND)   = 0.0
        enddo
      enddo
    enddo
  endif

  ! *** 3D ADVECTI0N-DIFFUSION TRANSPORT FOR ALL WATER COLUMN CONSITUENTS (SEE VARINIT FOR WCV INITIALIZATION)
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(IW,IT,L) SCHEDULE(STATIC,1)
  do IW = 1,NACTIVEWC
    !$  IT = OMP_GET_THREAD_NUM() + 1
    call CALTRAN( IACTIVEWC1(IW), IACTIVEWC2(IW), WCV(IW).VAL0, WCV(IW).VAL1, IW, IT, WCV(IW).WCLIMIT, ISKIP(IW) )
    
    if( ISICE > 2 .and. IACTIVEWC1(IW) == 8 .and. IACTIVEWC2(IW) == MSVDOX .and. ISKIP(IW) == 0 )then
      ! *** ZERO SURFACE MELT FLUX
      do L = 1,LC
        FQC(L,KC,IT) = 0.
      enddo
    endif
    
    ! *** Transport bypass reporting
    if( BYPASS )then
      if( NITER < 2 .or. TIMEDAY > HOUR12NEXT )then
        if( IS2TL > 0 .or. ISTL == 3 )then
          PRINT '(A,2I5,A10,E12.4,L5)', ' ISTL, THREAD, ID, WCLIMIT, SKIP: ', ISTL, IT, WCV(IW).ID, WCV(IW).WCLIMIT, (ISKIP(IW) == 1 )
        endif
      endif
    endif
  enddo
  !$OMP END DO

  ! *** APPLY ANTI-DIFFUSION AND FLUX CORRECTOR
  !$OMP DO PRIVATE(IW,IT) SCHEDULE(STATIC,1)
  do IW = 1,NACTIVEWC
    if( ISKIP(IW) == 0 )then
      !$  IT = OMP_GET_THREAD_NUM() + 1
      call CALTRAN_AD( IACTIVEWC1(IW), IACTIVEWC2(IW), WCV(IW).VAL0, WCV(IW).VAL1, IW, IT )
    endif
  enddo
  !$OMP END DO

  ! ****************************************************************************
  ! *** MPI communication for FUHUD, FVHUD & FWUU
  !$OMP SINGLE
  call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)

  call Communicate_CON2
  
  TMPITMP = DSTIME(0) - TTDS
  DSITIMING(6) = DSITIMING(6) + TMPITMP
  !$OMP END SINGLE
  ! ****************************************************************************

  !$OMP DO PRIVATE(IW,K,LP,L,CDTMP) SCHEDULE(STATIC,1)
  do IW = 1,NACTIVEWC
    if( ISADAC(IACTIVEWC1(IW)) == 1 )then 
      ! *** APPLY THE ANTI-DIFFUSIVE ADVECTION CALCULATION TO STANDARD DONOR CELL SCHEME
      if( ISKIP(IW) == 0 )then
        do K = 1,KC  
          do LP = 1,LLWET(K,0)
            L = LKWET(LP,K,0)
            CDTMP             = WCV(IW).VAL0(L,K)*HPK(L,K) + DELT*( ( FUHUD(L,K,IW)-FUHUD(LEC(L),K,IW) + FVHUD(L,K,IW)-FVHUD(LNC(L),K,IW) ) * DXYIP(L)   &
                                                           +        ( FWUU(L,K-1,IW)-FWUU(L,K,IW) )  )  
            WCV(IW).VAL0(L,K) = CDTMP*HPKI(L,K)
          enddo  
        enddo    

        ! *** RESET OPEN BC CONCENTRATIONS
        do K = 1,KC
          do LP = 1,NBCSOP
            L = LOBCS(LP)
            WCV(IW).VAL0(L,K) = WQBCCON(LP,K,IW)
          enddo
        enddo
      endif
    endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL  
  
  ! *** Advection of the turbulence quantities
  if( ISGOTM > 0 )then
    
    !$OMP SECTIONS PRIVATE(IT)  
    !$OMP SECTION
    !$  IT = OMP_GET_THREAD_NUM() + 1
    ! *** TKE transport section
    TKE_TEMP  = tke3d (1:LCM,1:KCM)
    TKE_TEMP1 = tke3d1(1:LCM,1:KCM)
    call CALTRAN(0, 0, tke_temp, tke_temp1, NACTIVEWC+1, IT, 0.0, ISKIP(NACTIVEWC+1))
    
    tke3d (1:LCM,1:KCM) = TKE_TEMP
    if( IS2TIM == 0 .and. ISTL /= 2 )then
      tke3d1 (1:LCM,1:KCM) = TKE_TEMP1
    endif
    tke3d = max(k_min,tke3d)
    
    ! *** EPS transport section
    !$OMP SECTION
    !$  IT = OMP_GET_THREAD_NUM() + 1
    EPS_TEMP  = eps3d (1:LCM,1:KCM)
    EPS_TEMP1 = eps3d1(1:LCM,1:KCM)
    call CALTRAN(0, 0, eps_temp, eps_temp1, NACTIVEWC+2, IT, 0.0, ISKIP(NACTIVEWC+2))
    
    eps3d (1:LCM,1:KCM) = EPS_TEMP
    if( IS2TIM == 0 .and. ISTL /= 2 )then
      eps3d1 (1:LCM,1:KCM) = EPS_TEMP1
    endif
    eps3d = max(eps_min,eps3d)

    !$OMP END SECTIONS
  endif
  
  if( ISICE == 4 .and. (IS2TL > 0 .or. (IS2TL == 0 .and. NCTBC /= NTSTBC)) )then
    if( LFRAZIL )then
      ! *** FRAZIL ICE TRANSPORT
      call CALTRANICE(FRAZILICE,FRAZILICE1,1)
      NICE = NICE+1

    elseif( NICE > 0 )then
      ! *** ADVANCE THE FRAZIL ICE VARIABLE
      FRAZILICE1 = FRAZILICE
      NICE = 0
    endif
  endif

  TTDS2 = DSTIME(0)
  call MPI_barrier(MPI_Comm_World, ierr)
  TWAIT = TWAIT + (DSTIME(0)- TTDS2)

  TSADV = TSADV + (DSTIME(0)-TTDS1) - TMPITMP - TWAIT

  ! ******************************************************************************************
  ! *** VERTICAL DIFFUSION IMPLICIT HALF STEP CALCULATION
  if( KC > 1 .and. NACTIVEWC > 0 )then

    TTDS1 = DSTIME(0)

    if( LADRY > 0 )then
      do LP = 1,LADRY
        L = LDRY(LP)
        CCLBTMP(L,:) = 0.
        EEB(L,:) = 1.
        VCU(L,:) = 0.
      enddo
    endif
  
    !$OMP PARALLEL DO PRIVATE(ND,LF,LL,LP,L,K,RCDZKK,CCUBTMP,CCMBTMP,RCDZKMK,IW) SCHEDULE(STATIC,1)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)
    
      ! -------------------------------------------------------------------------------
      ! *** COMPUTE THE DIFFUSIVE FLUXES
    
      ! *** BOTTOM LAYER
      do LP = 1,LLWET(KS,ND)
        L = LKWET(LP,KS,ND)
        RCDZKK  = -DELT*CDZKK(L,KSZ(L))
        CCUBTMP = RCDZKK*HPI(L)*AB(L,KSZ(L))
        CCMBTMP = 1._8-CCUBTMP
        EEB(L,KSZ(L)) = 1._8/CCMBTMP
        VCU(L,KSZ(L)) = CCUBTMP*EEB(L,KSZ(L))
      enddo
    
      ! *** MIDDLE LAYERS
      do K = 2,KS
        do LP = 1,LLWET(K-1,ND)
          L = LKWET(LP,K-1,ND)
          RCDZKMK      = -DELT*CDZKMK(L,K)
          RCDZKK       = -DELT*CDZKK(L,K)
          CCLBTMP(L,K) = RCDZKMK*HPI(L)*AB(L,K-1)
          CCUBTMP      = RCDZKK*HPI(L)*AB(L,K)
          CCMBTMP      = 1._8-CCLBTMP(L,K)-CCUBTMP
          EEB(L,K)     = 1._8/( CCMBTMP - CCLBTMP(L,K)*VCU(L,K-1) )
          VCU(L,K)     = CCUBTMP*EEB(L,K)
        enddo
      enddo
    
      ! *** TOP LAYER
      K = KC
      do LP = 1,LLWET(KS,ND)
        L = LKWET(LP,KS,ND)
        RCDZKMK      = -DELT*CDZKMK(L,K)
        CCLBTMP(L,K) = RCDZKMK*HPI(L)*AB(L,K-1)
        CCMBTMP      = 1._8-CCLBTMP(L,K)
        EEB(L,K)     = 1._8/( CCMBTMP - CCLBTMP(L,K)*VCU(L,K-1) )
      enddo
    
      ! -------------------------------------------------------------------------------
      ! *** APPLY THE DIFFUSION
    
      ! *** BOTTOM LAYER
      do IW = 1,NACTIVEWC
        if( ISKIP(IW) == 1 ) CYCLE
        do LP = 1,LLWET(KS,ND)
          L = LKWET(LP,KS,ND)
          K = KSZ(L)
          WCV(IW).VAL0(L,K) = WCV(IW).VAL0(L,K)*EEB(L,K)
        enddo
    
        ! *** MIDDLE AND TOP LAYERS
        do K = 2,KC
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND)
            WCV(IW).VAL0(L,K) = ( WCV(IW).VAL0(L,K) - CCLBTMP(L,K)*WCV(IW).VAL0(L,K-1) )*EEB(L,K)
          enddo
        enddo
    
        ! *** FINAL PASS
        do K = KS,1,-1
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            WCV(IW).VAL0(L,K) = WCV(IW).VAL0(L,K) - VCU(L,K)*WCV(IW).VAL0(L,K+1)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    TVDIF = TVDIF + (DSTIME(0)-TTDS1)
    
  endif    ! *** END OF VERTICAL DIFFUSION STEP
  
  ! ***  COMPUTE TOTALS FOR SEDIMENT TRANSPORT (AFTER ADVECTION/DIFFUSION)
  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    !$OMP PARALLEL DEFAULT(SHARED)
    ! *** ZERO SED SEDIMENT ACCUMULATION ARRAYS
    if( ISTRAN(6) > 0 )then
      !$OMP DO PRIVATE(ND,L,K,LF,LL,LP,NS) SCHEDULE(STATIC,1)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)

        ! *** WATER COLUMN
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SEDT(L,K) = 0.
          enddo
        enddo

        do NS = 1,NSED2
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              SEDT(L,K) = SEDT(L,K) + SED(L,K,NS)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
    endif

    ! *** ZERO SND SEDIMENT ACCUMULATION ARRAYS
    if( ISTRAN(7) > 0 )then
      !$OMP DO PRIVATE(ND,L,K,LF,LL,LP,NS) SCHEDULE(STATIC,1)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = MIN(LF+LDMWET-1,LAWET)

        ! *** WATER COLUMN
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            SNDT(L,K) = 0.
          enddo
        enddo

        do NS = 1,NSND
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              SNDT(L,K) = SNDT(L,K) + SND(L,K,NS)
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
    endif
    !$OMP END PARALLEL
  endif
        
  ! ***************************************************************************************
  ! *** SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION
  if( ISTRAN(2) >= 1 )then
    TTDS1 = DSTIME(0)
    call CALHEAT

    THEAT = THEAT + (DSTIME(0)-TTDS1)
  endif

  ! *** APPLY DYE PROCESSES
  if( ISTRAN(3) >= 1 ) CALL CALDYE

  ! **************************************************************************************************
  ! *** BOTTOM AND INTERNAL SEDIMENT AND TOXIC CONTAMINATION SOURCE-SINK CALCULATION

  if( ISPROPWASH == 2 )then
    ! *** If using propeller momentum, compute propeller velocities at every timestep, but skip if computing below
    if( .not. ((SEDTIME+DELT) >= SEDSTEP .and. (ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 ) .and. TIMEDAY >= SEDSTART) )then
      call Propwash_Calc_Sequence(1)
    endif
  endif

  ! *** SEDIMENT AND TOXICS SETTLING,DEPOSITION,RESUSPENSION,ETC
  if( ( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 ) .and. TIMEDAY >= SEDSTART )then
    if( IS2TIM >= 1 )then
      ! *** FOR TWO TIME LEVEL SIMULATION
      SEDTIME = SEDTIME + DELT
      if( SEDTIME >= SEDSTEP )then
        DTSED = SEDTIME
        
        ! *** Call Propwash module if that option is selected
        if( propwash_on .and. IFIRST > 0 )then
          call Propwash_Calc_Sequence(0)
        endif
        
        call SSEDTOX
        SEDTIME = 0.0
      endif
    else
      ! *** FOR THREE TIME LEVEL SIMULATION
      if( NCTBC == 1 )then
        DTSED = FLOAT(NTSTBC)*DT
        
        ! *** Call Propwash module if that option is selected
        if( propwash_on .and. IFIRST > 0 )then
          call Propwash_Calc_Sequence(0)
        endif
        
        call SSEDTOX
        SEDTIME = 0.0
      endif
    endif
  endif
  IFIRST = 1

  ! *** 2014-09 - REMOVED DATA ASSIMILATION CODE (PMC)

  return

END

