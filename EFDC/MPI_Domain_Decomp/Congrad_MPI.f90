! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! @details Modified conjugate gradient method borrow from O'Donncha's code
! @author Zander Mausolff
! @date 9/4/2018
!---------------------------------------------------------------------------!

SUBROUTINE CONGRAD_MPI(NOPTIMAL,LDMOPT)

  USE GLOBAL
  USE EFDCOUT
  
  USE MPI
  Use Variables_MPI
  Use Communicate_Ghost_Routines
  Use MPI_All_Reduce
  Use Variables_MPI_Mapping
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NOPTIMAL, LDMOPT

  ! *** Local Variables
  INTEGER       :: ND, LL, LF, L, LE, LS, LN, LW, I, J
  INTEGER       :: IERR

  Real  :: ALPHA, BETA
  REAL  :: RPCGN, RPCG, PAPCG, RSQ
  REAL(8) :: RPCGN_OUT, RPCG_OUT, PAPCG_OUT, RSQ_OUT, TWAIT
  
  Real(RKD), EXTERNAL   :: DSTIME
  REAL(RKD)             :: TTDS, TTDS1, TTDS2, TMPCOMM, TMPALLREDUCE         ! MODEL TIMING TEMPORARY VARIABLES

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PCG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RCG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TMPCG

  IF( .NOT. ALLOCATED(PCG) )THEN
    ALLOCATE(PCG(LCM))
    ALLOCATE(RCG(LCM))
    ALLOCATE(TMPCG(LCM))

    PCG   = 0.0
    RCG   = 0.0
    TMPCG = 0.0 
  ENDIF

  ! *** START THE TIMING
  TTDS1 = DSTIME(0)
  TMPCOMM = 0.
  TMPALLREDUCE = 0.
  
  !$OMP SIMD
  DO L=2,LA
    RCG(L) = FPTMP(L) - CCC(L)*P(L) - CCN(L)*P(LNC(L)) - CCS(L)*P(LSC(L)) - CCW(L)*P(LWC(L)) - CCE(L)*P(LEC(L))
  ENDDO

  !$OMP SIMD
  DO L=2,LA
    PCG(L) = RCG(L)*CCCI(L)
  ENDDO

  ! *** Exchange computed RCG ghost cells
  TTDS2 = DSTIME(0)
  Call Communicate_1D2(RCG, PCG)
  TMPCOMM = TMPCOMM + (DSTIME(0)-TTDS2)
  
  RPCG = 0.0
  !$OMP SIMD REDUCTION(+:RPCG)
  DO L=2,LA 
    RPCG = RPCG + GhostMask(L)*RCG(L)*PCG(L)  
  ENDDO

  ! *** Sum PAPCG values across all processes
  CALL DSI_All_Reduce(RPCG, RPCG_OUT, MPI_SUM, TTDS, 0, TWAIT)
  TMPALLREDUCE = TMPALLREDUCE + TTDS
  RPCG = RPCG_OUT
  
  IF( RPCG == 0.0 )RETURN   ! *** DSI SINGLE LINE
  
  !--------------------------------------------------------------------------------
  ! *** BEGIN THE ITERATIVE SOLUTION LOOP
  DO ITER=1,ITERM

    PAPCG = 0.
    RPCGN = 0.
    RSQ   = 0.
    
    !$OMP SIMD
    DO L=2,LA
      APCG(L) = CCC(L)*PCG(L) + CCS(L)*PCG(LSC(L)) + CCN(L)*PCG(LNC(L)) + CCW(L)*PCG(LWC(L)) + CCE(L)*PCG(LEC(L))
    ENDDO

    !$OMP SIMD REDUCTION(+:PAPCG)
    DO L=2,LA
      PAPCG = PAPCG + GhostMask(L)*APCG(L)*PCG(L)  
    ENDDO    
    
    ! *** Sum PAPCG values across all proceses
    CALL DSI_All_Reduce(PAPCG, PAPCG_OUT, MPI_SUM, TTDS, 0, TWAIT)
    TMPALLREDUCE = TMPALLREDUCE + TTDS
    PAPCG = PAPCG_OUT
    
    ALPHA = RPCG/PAPCG

    ! *** RPCGN and RSQ
    !$OMP SIMD
    DO L=2,LA
      P(L)     = P(L) + ALPHA*PCG(L)
      RCG(L)   = RCG(L) - ALPHA*APCG(L)
      TMPCG(L) = CCCI(L)*RCG(L)
    ENDDO

    !$OMP SIMD REDUCTION(+:RPCGN)
    DO L=2,LA
      RPCGN = RPCGN + GhostMask(L)*RCG(L)*TMPCG(L)  
    ENDDO  
    
    !$OMP SIMD REDUCTION(+:RSQ)
    DO L=2,LA
      RSQ   = RSQ   + GhostMask(L)*RCG(L)*RCG(L)  
    ENDDO  
    
    ! *** Sum all values across RPCGN and RSQ_DBL
    CALL DSI_All_Reduce(RPCGN, RPCGN_OUT, MPI_SUM, TTDS, 0, TWAIT)
    TMPALLREDUCE = TMPALLREDUCE + TTDS
    RPCGN = RPCGN_OUT

    CALL DSI_All_Reduce(RSQ, RSQ_OUT, MPI_SUM, TTDS, 0, TWAIT)
    TMPALLREDUCE = TMPALLREDUCE + TTDS
    RSQ  = RSQ_OUT

    BETA = RPCGN/RPCG
    RPCG = RPCGN
    
    IF( RSQ > RSQM )THEN
      ! *** PREPARE THE NEXT ITERATION
      !$OMP SIMD
      DO L=2,LA
        PCG(L) = TMPCG(L) + BETA*PCG(L)
      ENDDO
    ENDIF

    ! *** IF the tolerance on the residual has been met then exit the iteration loop
    IF( RSQ <= RSQM )THEN
      EXIT
    ENDIF

    Call MPI_barrier(MPI_Comm_World, IERR)
    TTDS2 = DSTIME(0)
    CALL Communicate_Ghost_Cells(PCG, 'PCG')
    TMPCOMM = TMPCOMM + (DSTIME(0)-TTDS2)

  ENDDO  ! *** END OF ITERATION LOOP

  IF( RSQ > RSQM )THEN
    WRITE(6,600) 
 
    ! *** Write to unique log per process 
    Open(unit_efdc_out, FILE=OUTDIR//filename_out,STATUS='OLD')
    WRITE(unit_efdc_out,610)

    DO L=2,LA  
      WRITE(unit_efdc_out,800)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),FPTMP(L)  
    ENDDO
    CLOSE(unit_efdc_out)

    CALL STOPP('')

600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
610 FORMAT('   I    J     CCS     CCW      CCC      CCE      CCN     FPTMP')
800 FORMAT(2I6,6E13.4)

  ENDIF

  TCONG = TCONG + (DSTIME(0)-TTDS1) - TMPCOMM - TMPALLREDUCE
  DSITIMING(1)  = DSITIMING(1)  + TMPCOMM
  DSITIMING(12) = DSITIMING(12) + TMPALLREDUCE

  IDRYTBP = ITER

  RETURN

End Subroutine Congrad_MPI
