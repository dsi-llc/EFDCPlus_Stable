! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CONGRAD(NOPTIMAL,LDMOPT)  

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
  
  REAL               :: ALPHA, BETA, RSQ
  REAL               :: RPCGN, RPCG, PAPCG
  
  REAL(RKD), EXTERNAL   :: DSTIME 
  REAL(RKD)             :: TTDS         ! MODEL TIMING TEMPORARY VARIABLE
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PCG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RCG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TMPCG  
  
  IF(  .NOT. ALLOCATED(PCG) )THEN
    ALLOCATE(PCG(LCM))  
    ALLOCATE(RCG(LCM))  
    ALLOCATE(TMPCG(LCM))
    PCG=0.0
    RCG=0.0
    TMPCG=0.0 
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
    
  RPCG=0.  
  !DIR$ SIMD
  DO L=2,LA 
    RPCG = RPCG + RCG(L)*PCG(L)  
  ENDDO

  IF( RPCG == 0.0 )RETURN   ! *** DSI SINGLE LINE

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
    !DIR$ SIMD
    DO L=2,LA
      PAPCG = PAPCG + APCG(L)*PCG(L)  
    ENDDO  
    ALPHA=RPCG/PAPCG  
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
    !DIR$ SIMD
    DO L=2,LA
      RPCGN = RPCGN+RCG(L)*TMPCG(L)  
      RSQ   = RSQ  +RCG(L)*RCG(L)  
    ENDDO
    BETA=RPCGN/RPCG  
    RPCG=RPCGN  
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
  
  IF( RSQ > RSQM )THEN
    ! *** WET/DRY MAXIMUM ITERATIONS EXCEEDED   
    WRITE(6,*) '  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION' 

    ! *** SAVE A SNAPSHOT FOR EE
    IF( ISPPH == 1 )THEN
      WRITE(6,*) '  THE LATEST MODEL RESULTS HAVE BEEN SAVED TO THE EE LINKAGE'
      Call Map_Write_EE_Binary
      if( process_id == master_id )THEN
        CALL EE_LINKAGE(-1)  
      endif
    ENDIF
        
    CALL STOPP('')  
  ENDIF

  TCONG = TCONG + (DSTIME(0)-TTDS)  
  IDRYTBP = ITER
  
  RETURN  

END  

