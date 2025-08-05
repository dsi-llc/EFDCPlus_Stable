! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CONGRADC()

  ! ***  SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE
  ! ***  GRADIENT SCHEME

  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2014-08           D H CHUNG          SET EXPLICIT PRECISIONS OF INTEGER & REAL
  ! 2011-03           Paul M. Craig      Rewritten to F90

  use GLOBAL
  use EFDCOUT
  use Variables_MPI
  use Mod_Map_Write_EE_Binary

  ! *** DSI
  implicit none
  
  integer    :: L, NMD
  integer    :: LHOST,LCHNU,LCHNV
  
  real       :: RPCGN, RPCG, PAPCG, ALPHA, BETA, RSQ
  real(RKD)  :: TTDS         ! MODEL TIMING TEMPORARY VARIABLE
  real(RKD), external :: DSTIME

  real,save,allocatable,dimension(:) :: PNORTH
  real,save,allocatable,dimension(:) :: PSOUTH
  real,save,allocatable,dimension(:) :: TMPCG
  real,save,allocatable,dimension(:) :: PCG
  real,save,allocatable,dimension(:) :: RCG
  real,save,allocatable,dimension(:) :: SRPCGN
  real,save,allocatable,dimension(:) :: SRSQ

  if( .not. allocated(PNORTH) )then
    allocate(PNORTH(LCM))
    allocate(PSOUTH(LCM))
    allocate(PCG(LCM))
    allocate(RCG(LCM))
    allocate(SRPCGN(NDM))
    allocate(SRSQ(NDM))
    allocate(TMPCG(LCM))
    PNORTH = 0.0
    PSOUTH = 0.0
    PCG = 0.0
    RCG = 0.0
    SRSQ = 0.0
    TMPCG = 0.0
  endif

  TTDS = DSTIME(0)
  do L = 2,LA
    PNORTH(L) = P(LNC(L))
    PSOUTH(L) = P(LSC(L))
  enddo
  do L = 2,LA
    RCG(L) = FPTMP(L) - CCC(L)*P(L) - CCN(L)*PNORTH(L) - CCS(L)*PSOUTH(L) &
                                    - CCW(L)*P(LWC(L)) - CCE(L)*P(LEC(L))
  enddo
  if( MDCHH >= 1 )then
    do NMD = 1,MDCHH
      LHOST = LMDCHH(NMD)
      LCHNU = LMDCHU(NMD)
      LCHNV = LMDCHV(NMD)

      !         X-DIRECTION CHANNEL
      if( MDCHTYP(NMD) == 1 )then
        RCG(LCHNU) = RCG(LCHNU)+CCCCHH(NMD)*P(LHOST)
        RCG(LHOST) = RCG(LHOST)+CCCCHH(NMD)*P(LCHNU)
      endif

      !         Y-DIRECTION CHANNEL
      if( MDCHTYP(NMD) == 2 )then
        RCG(LCHNV) = RCG(LCHNV)+CCCCHH(NMD)*P(LHOST)
        RCG(LHOST) = RCG(LHOST)+CCCCHH(NMD)*P(LCHNV)
      endif
    enddo
  endif

  do L = 2,LA
    PCG(L) = RCG(L)*CCCI(L)
  enddo
  RPCG = 0.
  do L = 2,LA
    RPCG = RPCG+RCG(L)*PCG(L)
  enddo
  if( RPCG == 0.0) return 
  ITER = 0
  
100 continue
  ITER = ITER+1
  do L = 2,LA
    PNORTH(L) = PCG(LNC(L))
    PSOUTH(L) = PCG(LSC(L))
  enddo
  do L = 2,LA
    APCG(L) = CCC(L)*PCG(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L) &
        +CCW(L)*PCG(LWC(L))+CCE(L)*PCG(LEC(L))
  enddo
  if( MDCHH >= 1 )then
    do NMD = 1,MDCHH
      LHOST = LMDCHH(NMD)
      LCHNU = LMDCHU(NMD)
      LCHNV = LMDCHV(NMD)

      !         X-DIRECTION CHANNEL
      if( MDCHTYP(NMD) == 1 )then
        APCG(LCHNU) = APCG(LCHNU)+CCCCHH(NMD)*PCG(LHOST)
        APCG(LHOST) = APCG(LHOST)+CCCCHH(NMD)*PCG(LCHNU)
      endif

      !         Y-DIRECTION CHANNEL
      if( MDCHTYP(NMD) == 2 )then
        APCG(LCHNV) = APCG(LCHNV)+CCCCHH(NMD)*PCG(LHOST)
        APCG(LHOST) = APCG(LHOST)+CCCCHH(NMD)*PCG(LCHNV)
      endif
    enddo
  endif

  PAPCG = 0.
  do L = 2,LA
    PAPCG = PAPCG+APCG(L)*PCG(L)
  enddo
  ALPHA = RPCG/PAPCG
  do L = 2,LA
    P(L) = P(L)+ALPHA*PCG(L)
  enddo
  do L = 2,LA
    RCG(L) = RCG(L)-ALPHA*APCG(L)
  enddo
  do L = 2,LA
    TMPCG(L) = CCCI(L)*RCG(L)
  enddo
  RPCGN = 0.
  RSQ = 0.
  do L = 2,LA
    RPCGN = RPCGN+RCG(L)*TMPCG(L)
    RSQ = RSQ+RCG(L)*RCG(L)
  enddo
  if( RSQ  <=  RSQM) GOTO 200
  if( ITER  >=  ITERM )then
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
    
    STOP
  endif
  
  BETA = RPCGN/RPCG
  RPCG = RPCGN
  do L = 2,LA
    PCG(L) = TMPCG(L)+BETA*PCG(L)
  enddo
  
  GOTO 100
600 FORMAT('  MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION.  MAX ERROR AT LG = ',I6,', PROCESS ID = ',I4)
610 FORMAT('    LG     L     I     J          CCS          CCW          CCC          CCE          CCN       FPTMP'    &
                             //     '       RCG^2            P          RCG        TMPCG          PCG')
800 FORMAT(4I6,12E13.4)

  ! *** CALCULATE FINAL RESIDUAL
  200 continue

  if( ISLOG >= 1 )then
    do L = 2,LA
      PNORTH(L) = P(LNC(L))
      PSOUTH(L) = P(LSC(L))
    enddo
    RSQ = 0.
    do L = 2,LA
      RCG(L) = CCC(L)*P(L)+CCS(L)*PSOUTH(L)+CCN(L)*PNORTH(L) &
        +CCW(L)*P(LWC(L))+CCE(L)*P(LEC(L))-FPTMP(L)
    enddo
    if( MDCHH >= 1 )then
      do NMD = 1,MDCHH
        LHOST = LMDCHH(NMD)
        LCHNU = LMDCHU(NMD)
        LCHNV = LMDCHV(NMD)

        !         X-DIRECTION CHANNEL
        if( MDCHTYP(NMD) == 1 )then
          RCG(LCHNU) = RCG(LCHNU)-CCCCHH(NMD)*P(LHOST)
          RCG(LHOST) = RCG(LHOST)-CCCCHH(NMD)*P(LCHNU)
        endif

        !         Y-DIRECTION CHANNEL
        if( MDCHTYP(NMD) == 2 )then
          RCG(LCHNV) = RCG(LCHNV)-CCCCHH(NMD)*P(LHOST)
          RCG(LHOST) = RCG(LHOST)-CCCCHH(NMD)*P(LCHNV)
        endif
      enddo
    endif
    do L = 2,LA
      RCG(L) = RCG(L)*CCCI(L)
    enddo
    do L = 2,LA
      RSQ = RSQ+RCG(L)*RCG(L)
    enddo
  endif

  TCONG = TCONG+(DSTIME(0)-TTDS)

  return
END

