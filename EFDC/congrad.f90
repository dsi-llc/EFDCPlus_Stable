! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CONGRAD(NOPTIMAL, LDMOPT)  

  ! SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE GRADIENT SCHEME 

  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP

  USE GLOBAL  
  USE EFDCOUT
  Use Variables_MPI
  Use Mod_Map_Write_EE_Binary
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: NOPTIMAL,LDMOPT
  
  INTEGER            :: ND, LF, LL, L, LE, LS, LN, LW
  
  REAL               :: ALPHA, BETA, RSQ, RSQ0
  REAL               :: RPCGN, RPCG, PAPCG, PMC, PMC2
  
  REAL(RKD), EXTERNAL   :: DSTIME 
  REAL(RKD)             :: TTDS         ! MODEL TIMING TEMPORARY VARIABLE
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PCG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RCG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TMPCG  
  
  ! *** delme
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: ALPHAs
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: BETAs
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RPCGs  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RSQs  
  ! *** delme

  IF( .NOT. ALLOCATED(PCG) )THEN
    ALLOCATE(PCG(LCM))  
    ALLOCATE(RCG(LCM))  
    ALLOCATE(TMPCG(LCM))
    PCG=0.0
    RCG=0.0
    TMPCG=0.0 
    
    ! *** delme
    ALLOCATE(ALPHAs(ITERM+1))  
    ALLOCATE(BETAs(ITERM+1))  
    ALLOCATE(RPCGs(ITERM+1))  
    ALLOCATE(RSQs(ITERM+1))  
    ALPHAs = 0.
    BETAS = 0.
    RPCGs = 0.
    RSQs = -1.
    ! *** delme
  ENDIF

  ! *** START THE TIMING
  TTDS = DSTIME(0)  

  !$OMP PARALLEL DO PRIVATE(ND,L,LF,LL,LN,LS,LE,LW)
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    
    DO L=LF,LL  
      LN=LNC(L)
      LS=LSC(L)
      LE=LEC(L)
      LW=LWC(L)
      RCG(L) = FPTMP(L) - CCC(L)*P(L) - CCN(L)*P(LN) - CCS(L)*P(LS) - CCW(L)*P(LW) - CCE(L)*P(LE)  
    ENDDO  

    DO L=LF,LL 
      PCG(L) = RCG(L)*CCCI(L)
    ENDDO  
  ENDDO
  !$OMP END PARALLEL DO
    
  ! *** delme
  ALPHAs = 0.
  BETAS = 0.
  RPCGs = 0.
  RSQs = -1.
  ! delme
  
  RPCG=0.  
  DO L=2,LA 
    RPCG = RPCG + RCG(L)*PCG(L)  
  ENDDO

  IF( RPCG == 0.0 ) RETURN   ! *** DSI SINGLE LINE

  ! *** BEGIN THE ITERATIVE SOLUTION LOOP
  DO ITER=1,ITERM

    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(ND,L,LF,LL,LN,LS,LE,LW)
    DO ND=1,NOPTIMAL  
      LF=2+(ND-1)*LDMOPT  
      LL=MIN(LF+LDMOPT-1,LA)

      DO L=LF,LL 
        LN=LNC(L)
        LS=LSC(L)
        LE=LEC(L)
        LW=LWC(L)
        APCG(L) = CCC(L)*PCG(L) + CCS(L)*PCG(LS) + CCN(L)*PCG(LN) + CCW(L)*PCG(LW) + CCE(L)*PCG(LE)  
      ENDDO  
    ENDDO  ! *** END OF DOMAIN
    !$OMP END DO
    
    ! *** PULLED OUT OF DOMAIN LOOP TO PREVENT ROUNDOFF ERROR
    !$OMP SINGLE
    PAPCG=0.
    DO L=2,LA
      PAPCG = PAPCG + APCG(L)*PCG(L)  
    ENDDO  
    ALPHA = RPCG/PAPCG  
    ALPHAS(ITER) = ALPHA    ! *** delme
    !$OMP END SINGLE
    
    !$OMP DO PRIVATE(ND,L,LF,LL)
    DO ND=1,NOPTIMAL  
      LF=2+(ND-1)*LDMOPT  
      LL=MIN(LF+LDMOPT-1,LA)
      
      DO L=LF,LL 
        P(L)      = P(L) + ALPHA*PCG(L)
        RCG(L)    = RCG(L) - ALPHA*APCG(L)  
        TMPCG(L)  = CCCI(L)*RCG(L)  
      ENDDO  
    ENDDO
    !$OMP END DO

    ! *** PULLED OUT OF DOMAIN LOOP TO PREVENT ROUNDOFF ERROR
    !$OMP SINGLE
    RPCGN = 0.
    RSQ = 0.
    DO L=2,LA
      RPCGN = RPCGN + RCG(L)*TMPCG(L)  
      RSQ   = RSQ   + RCG(L)*RCG(L)  
    ENDDO
    BETA = RPCGN/RPCG  
    RPCG = RPCGN  
    IF( ITER == 1 ) RSQ0 = RSQ*1.E-6
    
    ! *** delme
    BETAS(ITER) = BETA
    RPCGS(ITER) = RPCG
    RSQs(ITER) = RSQ
    ! *** delme
    !$OMP END SINGLE
    
    IF( RSQ > RSQM )THEN
      ! *** PREPARE THE NEXT ITERATION
      !$OMP DO PRIVATE(ND,L,LF,LL)
      DO ND=1,NOPTIMAL  
        LF=2+(ND-1)*LDMOPT  
        LL=MIN(LF+LDMOPT-1,LA)
        
        DO L=LF,LL 
          PCG(L)=TMPCG(L)+BETA*PCG(L)  
        ENDDO  
      ENDDO
      !$OMP END DO
    ENDIF
    !$OMP END PARALLEL
      
    IF( RSQ <= RSQM )THEN
      EXIT
    ENDIF
  
  ENDDO  ! *** END OF ITERATION LOOP
  
  ! *** delme
  if( ISINWV == 2 )then
    DO L = 2,LA
      PMCTESTX(1,L) = P(L)
      PMCTESTX(2,L) = PCG(L)
      PMCTESTX(3,L) = RCG(L)
      PMCTESTX(4,L) = TMPCG(L)
    ENDDO
  endif
  ! *** delme

  IF( RSQ > RSQM )THEN
    ! *** WET/DRY MAXIMUM ITERATIONS EXCEEDED   
    L = MAXLOC(RCG,DIM=1)
    WRITE(6,600) Map2Global(L).LG, process_id
 
    print *, ' See CONGRAD.ERR for last iteration history'
    ITER = min(ITER,ITERM)
    
    Open(1000, FILE=OUTDIR//'CONGRAD.ERR' ,STATUS='UNKNOWN')
    WRITE(1000,'(a,i5,f12.4,i5,2e14.6)') ' *** Error. Too many conjugate solution iterations. PROCESS_ID TIMEDAY ITER RSQ RSQM =', PROCESS_ID, TIMEDAY, ITER, RSQ, RSQM
    DO LF = 1,ITER
      WRITE(1000,'(I5,4E14.6)') LF, ALPHAS(LF), BETAS(LF), RPCGS(LF), RSQS(LF)
    ENDDO 
    CLOSE(1000)
    
    ! *** Write to unique log per process 
    Open(unit_efdc_out, FILE=OUTDIR//filename_out ,STATUS='OLD')
    WRITE(unit_efdc_out,610)

    WRITE(unit_efdc_out,600) Map2Global(L).LG, process_id
    DO L=2,LA  
      WRITE(unit_efdc_out,800) Map2Global(L).LG, L, IL(L), JL(L), CCS(L), CCW(L), CCC(L), CCE(L), CCN(L), &
                               FPTMP(L), RCG(L)*RCG(L), P(L), RCG(L), TMPCG(L), PCG(L)
    ENDDO
    CLOSE(unit_efdc_out)

    ! *** SAVE A SNAPSHOT FOR EE
    IF( ISPPH == 1 )THEN
      WRITE(6,*) ' The latest model results have been saved to the EE linkage files.'
      Call Map_Write_EE_Binary
      if( process_id == master_id )THEN
        CALL EE_LINKAGE(-1)  
      endif
    ENDIF
        
    CALL STOPP('')  
    
600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION.  MAX ERROR AT LG = ',I6,', PROCESS ID = ',I4)
610 FORMAT('    LG     L     I     J          CCS          CCW          CCC          CCE          CCN       FPTMP'    &
                             //     '       RCG^2            P          RCG        TMPCG          PCG')
800 FORMAT(4I6,12E13.4)

  ENDIF

  TCONG = TCONG + (DSTIME(0)-TTDS)  
  
  RETURN  

END  

