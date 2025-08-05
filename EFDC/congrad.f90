! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CONGRAD

  ! SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE GRADIENT SCHEME 

  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP

  use GLOBAL  
  use EFDCOUT
  use Allocate_Initialize      
  use Variables_MPI
  use Mod_Map_Write_EE_Binary
  
  implicit none
  
  integer            :: ND, LF, LL, L, LE, LS, LN, LW
  
  real               :: ALPHA, BETA, RSQ, RSQ0
  real               :: RPCGN, RPCG, PAPCG, PMC, PMC2
  
  real(RKD), external :: DSTIME 
  real(RKD)           :: TTDS         ! MODEL TIMING TEMPORARY VARIABLE
  
  real,save,allocatable,dimension(:) :: PCG
  real,save,allocatable,dimension(:) :: RCG
  real,save,allocatable,dimension(:) :: TMPCG  
  
  ! *** Reporting
  real,save,allocatable,dimension(:) :: ALPHAs
  real,save,allocatable,dimension(:) :: BETAs
  real,save,allocatable,dimension(:) :: RPCGs  
  real,save,allocatable,dimension(:) :: RSQs  

  if( .not. allocated(PCG) )then
    call AllocateDSI(PCG,   LCM, 0.0)  
    call AllocateDSI(RCG,   LCM, 0.0)  
    call AllocateDSI(TMPCG, LCM, 0.0)
    
    ! *** Reporting
    call AllocateDSI(ALPHAs, ITERM+1,  0.0)  
    call AllocateDSI(BETAs,  ITERM+1,  0.0)  
    call AllocateDSI(RPCGs,  ITERM+1,  0.0)  
    call AllocateDSI(RSQs,   ITERM+1, -1.0)  
  endif

  ! *** START THE TIMING
  TTDS = DSTIME(0)  

  !$OMP PARALLEL DO PRIVATE(ND,L,LF,LL,LN,LS,LE,LW) NUM_THREADS(NOPTIMAL(0))
  do ND = 1,NOPTIMAL(0)  
    LF = 2+(ND-1)*LDMOPT(0)  
    LL = min(LF+LDMOPT(0)-1,LA)
    
    do L = LF,LL  
      LN = LNC(L)
      LS = LSC(L)
      LE = LEC(L)
      LW = LWC(L)
      RCG(L) = FPTMP(L) - CCC(L)*P(L) - CCN(L)*P(LN) - CCS(L)*P(LS) - CCW(L)*P(LW) - CCE(L)*P(LE)  
    enddo  

    do L = LF,LL 
      PCG(L) = RCG(L)*CCCI(L)
    enddo  
  enddo
  !$OMP END PARALLEL DO
    
  ! *** Initialize for reporting
  ALPHAs = 0.
  BETAS = 0.
  RPCGs = 0.
  RSQs = -1.
  
  RPCG = 0.  
  do L = 2,LA 
    RPCG = RPCG + RCG(L)*PCG(L)  
  enddo

  if( RPCG == 0.0 ) return   ! *** DSI SINGLE LINE

  ! *** BEGIN THE ITERATIVE SOLUTION LOOP
  do ITER = 1,ITERM

    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(PAPCG) NUM_THREADS(NOPTIMAL(0))   &
    !$OMP    SHARED(NOPTIMAL, LDMOPT, LA, ITER, LNC, LSC, LEC, LWC, APCG, PCG, CCC, CCS, CCN, CCW, CCE, CCCI)   &
    !$OMP    SHARED(P, ALPHA, ALPHAS, RCG, RPCG, RPCGN, RSQ, RSQ0, BETA, BETAS, RPCGS, RSQS, RSQM, TMPCG)
   
    !$OMP DO PRIVATE(ND,L,LF,LL,LN,LS,LE,LW)
    do ND = 1,NOPTIMAL(0)  
      LF = 2+(ND-1)*LDMOPT(0)  
      LL = min(LF+LDMOPT(0)-1,LA)

      do L = LF,LL 
        LN = LNC(L)
        LS = LSC(L)
        LE = LEC(L)
        LW = LWC(L)
        APCG(L) = CCC(L)*PCG(L) + CCS(L)*PCG(LS) + CCN(L)*PCG(LN) + CCW(L)*PCG(LW) + CCE(L)*PCG(LE)  
      enddo  
    enddo  ! *** END OF DOMAIN
    !$OMP END DO
    
    ! *** PULLED OUT OF DOMAIN LOOP TO PREVENT ROUNDOFF ERROR
    !$OMP SINGLE
    PAPCG = 0.
    do L = 2,LA
      PAPCG = PAPCG + APCG(L)*PCG(L)  
    enddo  
    ALPHA = RPCG/PAPCG  
    ALPHAS(ITER) = ALPHA
    !$OMP END SINGLE
    
    !$OMP DO PRIVATE(ND,L,LF,LL)   
    do ND = 1,NOPTIMAL(0)
      LF = 2+(ND-1)*LDMOPT(0)  
      LL = min(LF+LDMOPT(0)-1,LA)
      
      do L = LF,LL 
        P(L)      = P(L) + ALPHA*PCG(L)
        RCG(L)    = RCG(L) - ALPHA*APCG(L)  
        TMPCG(L)  = CCCI(L)*RCG(L)  
      enddo  
    enddo
    !$OMP END DO

    ! *** PULLED OUT OF DOMAIN LOOP TO PREVENT ROUNDOFF ERROR
    !$OMP SINGLE
    RPCGN = 0.
    RSQ = 0.
    do L = 2,LA
      RPCGN = RPCGN + RCG(L)*TMPCG(L)  
      RSQ   = RSQ   + RCG(L)*RCG(L)  
    enddo
    BETA = RPCGN/RPCG  
    RPCG = RPCGN  
    if( ITER == 1 ) RSQ0 = RSQ*1.E-6
    
    ! *** Reporting
    BETAS(ITER) = BETA
    RPCGS(ITER) = RPCG
    RSQs(ITER) = RSQ

    !$OMP END SINGLE
    
    if( RSQ > RSQM )then
      ! *** PREPARE THE NEXT ITERATION
      !$OMP DO PRIVATE(ND,L,LF,LL)
      do ND = 1,NOPTIMAL(0)  
        LF = 2+(ND-1)*LDMOPT(0)  
        LL = min(LF+LDMOPT(0)-1,LA)
        
        do L = LF,LL 
          PCG(L) = TMPCG(L)+BETA*PCG(L)  
        enddo  
      enddo
      !$OMP END DO
    endif
    !$OMP END PARALLEL
      
    if( RSQ <= RSQM )then
      exit
    endif
  
  enddo  ! *** END OF ITERATION LOOP

  if( RSQ > RSQM )then
    ! *** WET/DRY MAXIMUM ITERATIONS EXCEEDED   
    L = MAXLOC(RCG,DIM = 1)
    write(6,600) Map2Global(L).LG, process_id
 
    print *, ' See CONGRAD.ERR for last iteration history'
    ITER = min(ITER,ITERM)
    
    Open(1000, FILE = OUTDIR//'CONGRAD.ERR' ,STATUS = 'UNKNOWN')
    write(1000,'(a,i5,f12.4,i5,2e14.6)') ' *** Error. Too many conjugate solution iterations. PROCESS_ID TIMEDAY ITER RSQ RSQM  = ', PROCESS_ID, TIMEDAY, ITER, RSQ, RSQM
    do LF = 1,ITER
      write(1000,'(I5,4E14.6)') LF, ALPHAS(LF), BETAS(LF), RPCGS(LF), RSQS(LF)
    enddo 
    close(1000)
    
    ! *** Write to unique log per process 
    Open(mpi_error_unit, FILE = OUTDIR//mpi_error_file,STATUS = 'OLD')
    write(mpi_error_unit,600) Map2Global(L).LG, process_id
    
    write(mpi_error_unit,610)
    do L = 2,LA  
      write(mpi_error_unit,800) Map2Global(L).LG, L, IL(L), JL(L), CCS(L), CCW(L), CCC(L), CCE(L), CCN(L), FPTMP(L), RCG(L)
    enddo
    close(mpi_error_unit)

    ! *** save A SNAPSHOT FOR EE
    if( ISPPH == 1 )then
      write(6,*) ' The latest model results have been saved to the EE linkage files.'
      call Map_Write_EE_Binary
      if( process_id == master_id )then
        call EE_LINKAGE(-1)  
      endif
    endif
        
    call STOPP('')  
    
600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION.  MAX ERROR AT LG = ',I6,', PROCESS ID = ',I4)
610 FORMAT('    LG     L     I     J          CCS          CCW          CCC          CCE          CCN       FPTMP'    &
                             //     '       RCG^2            P          RCG        TMPCG          PCG')
800 FORMAT(4I6,12E13.4)

  endif

  TCONG = TCONG + (DSTIME(0)-TTDS)  
  
  return 

END  

