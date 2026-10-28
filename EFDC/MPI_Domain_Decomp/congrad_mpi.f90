! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! @details Modified conjugate gradient method borrow from O'Donncha's code
! @author Zander Mausolff
! @date 9/4/2018
!---------------------------------------------------------------------------!

SUBROUTINE CONGRAD_MPI()  !NOPTIMAL(0),LDMOPT(0))

  use GLOBAL
  use EFDCOUT
  
  use MPI
  use Variables_MPI
  use Communicate_Ghost_Routines
  use MPI_All_Reduce
  use Variables_MPI_Mapping
  use Mod_Map_Write_EE_Binary
  
  implicit none

  !INTEGER, intent(IN) :: NOPTIMAL, LDMOPT

  ! *** Local Variables
  integer       :: ND, LL, LF, L, LE, LS, LN, LW, I, J
  integer       :: IERR

  real  :: ALPHA, BETA
  real  :: RPCGN, RPCG, PAPCG, RSQ
  real(8) :: RPCGN_OUT, RPCG_OUT, PAPCG_OUT, RSQ_OUT, TWAIT
  
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, TTDS1, TTDS2, TMPCOMM, TMPALLREDUCE         ! MODEL TIMING TEMPORARY VARIABLES

  real,save,allocatable,dimension(:) :: PCG
  real,save,allocatable,dimension(:) :: RCG
  real,save,allocatable,dimension(:) :: TMPCG

  if( .not. allocated(PCG) )then
    allocate(PCG(LCM))
    allocate(RCG(LCM))
    allocate(TMPCG(LCM))

    PCG   = 0.0
    RCG   = 0.0
    TMPCG = 0.0 
  endif

  ! *** START THE TIMING
  TTDS1 = DSTIME(0)
  TMPCOMM = 0.
  TMPALLREDUCE = 0.
  
  do L = 2,LA
    RCG(L) = FPTMP(L) - CCC(L)*P(L) - CCN(L)*P(LNC(L)) - CCS(L)*P(LSC(L)) - CCW(L)*P(LWC(L)) - CCE(L)*P(LEC(L))
  enddo

  do L = 2,LA
    PCG(L) = RCG(L)*CCCI(L)
  enddo

  ! *** Exchange computed RCG ghost cells
  TTDS2 = DSTIME(0)
  call Communicate_1D2(RCG, PCG)
  TMPCOMM = TMPCOMM + (DSTIME(0)-TTDS2)
  
  RPCG = 0.0
  do L = 2,LA 
    RPCG = RPCG + GhostMask(L)*RCG(L)*PCG(L)  
  enddo

  ! *** Sum PAPCG values across all processes
  call DSI_All_Reduce(RPCG, RPCG_OUT, MPI_SUM, TTDS, 0, TWAIT)
  TMPALLREDUCE = TMPALLREDUCE + TTDS
  RPCG = RPCG_OUT
  
  if( RPCG == 0.0 ) return   ! *** DSI SINGLE LINE
  
  !--------------------------------------------------------------------------------
  ! *** BEGIN THE ITERATIVE SOLUTION LOOP
  do ITER = 1,ITERM

    PAPCG = 0.
    RPCGN = 0.
    RSQ   = 0.
    
    do L = 2,LA
      APCG(L) = CCC(L)*PCG(L) + CCS(L)*PCG(LSC(L)) + CCN(L)*PCG(LNC(L)) + CCW(L)*PCG(LWC(L)) + CCE(L)*PCG(LEC(L))
    enddo

    do L = 2,LA
      PAPCG = PAPCG + GhostMask(L)*APCG(L)*PCG(L)  
    enddo    
    
    ! *** Sum PAPCG values across all proceses
    call DSI_All_Reduce(PAPCG, PAPCG_OUT, MPI_SUM, TTDS, 0, TWAIT)
    TMPALLREDUCE = TMPALLREDUCE + TTDS
    PAPCG = PAPCG_OUT
    
    ALPHA = RPCG/PAPCG

    ! *** RPCGN and RSQ
    do L = 2,LA
      P(L)     = P(L) + ALPHA*PCG(L)
      RCG(L)   = RCG(L) - ALPHA*APCG(L)
      TMPCG(L) = CCCI(L)*RCG(L)
    enddo

    do L = 2,LA
      RPCGN = RPCGN + GhostMask(L)*RCG(L)*TMPCG(L)  
    enddo  
    
    do L = 2,LA
      RSQ   = RSQ   + GhostMask(L)*RCG(L)*RCG(L)  
    enddo  
    
    ! *** Sum all values across RPCGN and RSQ_DBL
    call DSI_All_Reduce(RPCGN, RPCGN_OUT, MPI_SUM, TTDS, 0, TWAIT)
    TMPALLREDUCE = TMPALLREDUCE + TTDS
    RPCGN = RPCGN_OUT

    call DSI_All_Reduce(RSQ, RSQ_OUT, MPI_SUM, TTDS, 0, TWAIT)
    TMPALLREDUCE = TMPALLREDUCE + TTDS
    RSQ  = RSQ_OUT

    BETA = RPCGN/RPCG
    RPCG = RPCGN
    
    if( RSQ > RSQM )then
      ! *** PREPARE THE NEXT ITERATION
      do L = 2,LA
        PCG(L) = TMPCG(L) + BETA*PCG(L)
      enddo
    endif

    ! *** IF the tolerance on the residual has been met then exit the iteration loop
    if( RSQ <= RSQM )then
      exit
    endif

    call MPI_barrier(DSIcomm, IERR)
    TTDS2 = DSTIME(0)
    call Communicate_Ghost_Cells(PCG, 'PCG')
    TMPCOMM = TMPCOMM + (DSTIME(0)-TTDS2)

  enddo  ! *** END OF ITERATION LOOP

  if( RSQ > RSQM )then
    L = MAXLOC(RCG,DIM = 1 )
    write(6,600) Map2Global(L).LG, process_id
 
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

600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION.  MAX ERROR AT L = ',I6,', PROCESS ID = ',I4)
610 FORMAT('    LG    L     I     J          CCS          CCW          CCC          CCE          CCN       FPTMP          RCG')
800 FORMAT(4I6,8E13.4)

  endif

  TCONG = TCONG + (DSTIME(0)-TTDS1) - TMPCOMM - TMPALLREDUCE
  DSITIMING(1)  = DSITIMING(1)  + TMPCOMM
  DSITIMING(12) = DSITIMING(12) + TMPALLREDUCE

  return

End Subroutine Congrad_MPI
