! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CALSTEPD

  ! CHANGE RECORD
  ! *** SUBROUTINE CALSTEP ESTIMATE THE CURRENT MAXIMUM TIME STEP SIZE
  ! *** FORM LINEAR STABILITY CRITERIA AND A FACTOR OF SAFETY
  !
  !----------------------------------------------------------------------C
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP

  use GLOBAL
  use Allocate_Initialize
  
  use Variables_MPI
  use MPI
  use Variables_MPI_Mapping
  use MPI_All_Reduce
  use Broadcast_Routines

  implicit none

  integer :: IERR, ND, L, K, LF, LL, LE, LN, LP, KM, LLOC, ITRNTMP, NX, NINCRMTLoc
  integer :: NMD, LMDCHHT, LMDCHUT, LMDCHVT, ITMPR, LLOCOLD,  MINTYPE
  integer,save :: NUP
  !INTEGER,save,allocatable,dimension(:) :: LDYN
  
  real(RKD)      :: QUKTMP, QVKTMP, RTMPR, TESTTEMP, DTMAXX
  real(RKD)      :: TMPUUU, DTTMP, TMPVVV, TMPVAL, THP1, THP2
  real(RKD)      :: TOP, QXPLUS, QYPLUS, QZPLUS, QXMINS, QYMINS, QZMINS, QTOTAL, QSRC, BOT
  real(RKD)      :: DTCOMP, DTWARN, DTDYNP
  real(RKD),save :: HPLIM, HPLIMOLD
  real(RKD),save,allocatable,dimension(:) :: DTL1
  real(RKD),save,allocatable,dimension(:) :: DTL2
  real(RKD),save,allocatable,dimension(:) :: DTL3
  real(RKD),save,allocatable,dimension(:) :: DTL4
  real(RKD),save,allocatable,dimension(:,:) :: QSUBINN
  real(RKD),save,allocatable,dimension(:,:) :: QSUBOUT
  real(RKD),save,allocatable,dimension(:) :: DTHISTORY

  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, TWAIT                 ! MODEL TIMING TEMPORARY VARIABLE
  
  if( .not. allocated(DTL1) )then
    call AllocateDSI(DTL1,    LCM,  2.*TIDALP)
    call AllocateDSI(DTL2,    LCM,  2.*TIDALP)
    call AllocateDSI(DTL3,    LCM,  2.*TIDALP)
    call AllocateDSI(DTL4,    LCM,  2.*TIDALP)
    call AllocateDSI(QSUBINN, LCM,  KCM,  0.0)
    call AllocateDSI(QSUBOUT, LCM,  KCM,  0.0)

    allocate(DTHISTORY(10))
   
    DTHISTORY = DT
    NUP = 0
  endif

  ! ********************************************************************************
  
  ! *** SET WATER COLUMN CONSTITUENT TRANSPORT FLAG
  ITRNTMP = SUM(ISTRAN)
  
  if( NITER <= 1 ) DTDYN = DT

  DTMIN = DT
  DTMAXX = TIDALP

  ! *** DETERMINE SOURCE/SINKS FOR SUBGRID SCALE CHANNEL EXCHANGES
  if( MDCHH >= 1 )then
    do K = 1,KC
      do L = 2,LA
        QSUBOUT(L,K) = 0.0
        QSUBINN(L,K) = 0.0
      enddo
    enddo

    do K = 1,KC
      do NMD = 1,MDCHH
        if( LKSZ(L,K) ) CYCLE
        LMDCHHT = LMDCHH(NMD)
        LMDCHUT = LMDCHU(NMD)
        LMDCHVT = LMDCHV(NMD)
        if( MDCHTYP(NMD) == 1 )then
          QUKTMP = QCHANU(NMD)*DZC(LMDCHUT,K)
          QVKTMP = 0.
        endif
        if( MDCHTYP(NMD) == 2 )then
          QVKTMP = QCHANV(NMD)*DZC(LMDCHVT,K)
          QUKTMP = 0.
        endif
        if( MDCHTYP(NMD) == 3 )then
          QUKTMP = QCHANU(NMD)*DZC(LMDCHUT,K)
          QVKTMP = QCHANV(NMD)*DZC(LMDCHVT,K)
        endif
        QSUBOUT(LMDCHHT,K) = QSUBOUT(LMDCHHT,K) + min(QUKTMP,0.) + min(QVKTMP,0.)
        QSUBINN(LMDCHHT,K) = QSUBINN(LMDCHHT,K) + max(QUKTMP,0.) + max(QVKTMP,0.)
        QSUBOUT(LMDCHUT,K) = QSUBOUT(LMDCHUT,K) - max(QUKTMP,0.)
        QSUBINN(LMDCHUT,K) = QSUBINN(LMDCHUT,K) - min(QUKTMP,0.)
        QSUBOUT(LMDCHVT,K) = QSUBOUT(LMDCHVT,K) - max(QVKTMP,0.)
        QSUBINN(LMDCHVT,K) = QSUBINN(LMDCHVT,K) - min(QVKTMP,0.)
      enddo
    enddo
  endif

  ! *** Initialize
  DTL1 = DTMAXX
  if( ITRNTMP > 0 )   DTL2 = DTMAXX
  if( DTSSDHDT > 0. ) DTL4 = DTMAXX
  
  !$OMP PARALLEL DO DEFAULT(NONE)                                                               &
  !$OMP  SHARED(NDM, LDM, LA, LAWET, LWET, KC, ITRNTMP, LLWET, LKWET, LEC, LNC, IsGhost)        &
  !$OMP  SHARED(DTSSDHDT, DTMAXX, DXYP, DZC, UHDY2, VHDX2, W2, QSUM, DTL1, DTL2, DTL3, DTL4)    &
  !$OMP  SHARED(HDRY, HP, H1P, UHE, VHE, HUI, HVI, DXIU, DYIV, QSUBINN, QSUBOUT, DTDYN)         &
  !$OMP  PRIVATE(ND, L, K, LF, LL, LP, LE, LN, KM)                                              &
  !$OMP  PRIVATE(TMPUUU, DTTMP, TMPVVV, TMPVAL, TESTTEMP, TOP, QXPLUS, QYPLUS, QZPLUS)          &
  !$OMP  PRIVATE(QXMINS, QYMINS, QZMINS, QTOTAL, QSRC, BOT, THP1, THP2)
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = min(LF+LDM-1,LA)

    ! *** Method 1: COURANT–FRIEDRICHS–LEWY
    do LP = 1,LAWET
      L = LWET(LP)   
      TMPVAL = 1.E-16
      TMPUUU = ABS(UHE(L))*HUI(L)*DXIU(L)
      TMPVVV = ABS(VHE(L))*HVI(L)*DYIV(L)
      TMPVAL = TMPUUU + TMPVVV
      if( TMPVAL > 0.  ) DTL1(L)  = 1./TMPVAL
    enddo

    ! *** Method 2: Positivity of advected material, DTL2
    if( ITRNTMP >= 1 )then
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          if( IsGhost(L) ) CYCLE     ! *** Do not use ghost cells to set min DT
          
          LE = LEC(L)
          LN = LNC(L)
          KM = K - 1
          
          ! *** Volume of layer
          TOP    = H1P(L)*DXYP(L)*DZC(L,K)
          
          ! *** Positive Fluxes
          QXPLUS = UHDY2(LE,K)
          QXPLUS = max(QXPLUS,0.0)
          QYPLUS = VHDX2(LN,K)
          QYPLUS = max(QYPLUS,0.0)
          QZPLUS = W2(L,K)*DXYP(L)
          QZPLUS = max(QZPLUS,0.0)
          
          ! *** Negative fluxes
          QXMINS = UHDY2(L,K)
          QXMINS = -MIN(QXMINS,0.0)
          QYMINS = VHDX2(L,K)
          QYMINS = -MIN(QYMINS,0.0)
          QZMINS = W2(L,KM)*DXYP(L)
          QZMINS = -MIN(QZMINS,0.0)
          
          QTOTAL = QSUM(L,K) + QSUBOUT(L,K) + QSUBINN(L,K)
          QSRC   = -MIN(QTOTAL,0.0)
          BOT = QXPLUS + QYPLUS + QZPLUS + QXMINS + QYMINS + QZMINS + QSRC
          
          if( BOT > 1.E-12 )then
            DTTMP = TOP/BOT
            DTL2(L) = min(DTL2(L),DTTMP)
          endif
        enddo
      enddo
    endif

    ! *** METHOD 3: IMPLICIT BOTTOM FRICTION AND ROTATIONAL ACCELERATION DAMPING
    !DO L = LF,LL
    !  if( LACTIVE(L) )then
    !    TMPVAL = SUB(L)+SUB(LEC(L))+SVB(L)+SVB(LNC(L))
    !    if( TMPVAL > 0.5 )then
    !      LN = LNC(L)
    !      TAUBC = QQ(L,0)/CTURB2
    !      UCTR = 0.5*(U(L,1)+U(LEC(L),1))
    !      VCTR = 0.5*(V(L,1)+V(LN,1))
    !      UHMAG = HP(L)*SQRT(UCTR*UCTR+VCTR*VCTR)
    !      if( UHMAG > 0.0 )then
    !        FRIFRE = TAUBC/UHMAG
    !        FRIFRE2 = FRIFRE*FRIFRE
    !        ACACTMP = (CAC(L,KC)*HPI(L)*DXYIP(L))**2
    !        if( ACACTMP > FRIFRE2 )then
    !          DTTMP = 2.*FRIFRE/(ACACTMP-FRIFRE2)
    !          DTL3(L) = min(DTL3(L),DTTMP)
    !        endif
    !      endif
    !    endif
    !  endif
    !ENDDO

    ! ***  Method 4: Limit rate of depth change
    if( DTSSDHDT > 0. )then
      do LP = 1,LLWET(KC,ND)
        L = LKWET(LP,KC,ND)  
        if( IsGhost(L) ) CYCLE     ! *** Do not use ghost cells to set min DT
        
        THP1 = HP(L) 
        if( THP1 < HDRY/10.) THP1 = HDRY/10.
        THP2 = H1P(L) 
        if( THP2 < HDRY/10.) THP2 = HDRY/10.
        TESTTEMP = max(ABS(THP1-THP2),1.E-06)
        TMPVAL   = DTDYN*HP(L)/TESTTEMP
        DTL4(L)  = DTSSDHDT*TMPVAL
      enddo
    endif
  enddo   ! *** END OF DOMAIN
  !$OMP END PARALLEL DO

  ! *** Choose the minimum of the three methods
  MINTYPE = 0

  L1LOC = MINLOC(DTL1,DIM = 1)
  DTL1MN = DTL1(L1LOC)

  ! *** Initialize Method 2, if advective transport is activated
  if( ITRNTMP >= 1 )then
    L2LOC = MINLOC(DTL2,DIM = 1)
    DTL2MN = DTL2(L2LOC)
  endif

  !DTL3MN = 2.*DTMAXX
  !DO L = 2,LA
  !  if( DTL3MN > DTL3(L) )then
  !    DTL3MN = DTL3(L)
  !    L3LOC = L
  !  endif
  !ENDDO

  DTL4MN = 2.*DTMAXX
  if( DTSSDHDT > 0. )then
    L4LOC = MINLOC(DTL4,DIM = 1)
    DTL4MN = DTL4(L4LOC)
  endif

  ! *** Find minimum & apply a safety factor
  DTL1MN = DTL1MN*DTSSFAC
  DTTMP  = 2.*DTMAXX               ! *** DTTMP - Minimum time step without safety factor
  if( DTTMP > DTL1MN )then
    DTTMP  = DTL1MN
    DTCOMP = 1./DTSSFAC
    LLOC   = L1LOC
    MINTYPE = 1
  endif

  if( ITRNTMP >= 1 )then
    DTL2MN = DTL2MN*0.5             ! *** Fixed safety factor for advection of 0.5
    if( DTTMP > DTL2MN )then
      DTTMP  = DTL2MN
      DTCOMP = 2.
      LLOC   = L2LOC
      MINTYPE = 2
    endif
  endif
  
  if( DTSSDHDT > 0. )then
    if( DTTMP > DTL4MN )then
      DTTMP  = DTL4MN
      DTCOMP = 1.0
      LLOC   = L4LOC
      MINTYPE = 4
    endif
  endif
  
  LLOCOLD = LMINSTEP              ! *** Global Old LMINSTEP

  ! *** Synchronize all processes with minimum timestep and corresponding global LMINSTEP
  TMPVAL = min(DTTMP,DTMAX)
  LL = MAP2GLOBAL(LLOC).LG        ! *** Global LMINSTEP for current process

  call MPI_barrier(DSIcomm, IERR)
  call DSI_All_Reduce(TMPVAL, LL, DTTMP, LMINSTEP, MPI_MIN, TTDS, 1, TWAIT)

  call Broadcast_Scalar(DTTMP, master_id)             ! *** Timestep with safety factor
  call Broadcast_Scalar(LMINSTEP, master_id)

  DSITIMING(11) = DSITIMING(11) + TTDS


  DTCOMP = DTTMP*DTCOMP           ! *** Backout the safety factor to have the minimum delta T without any safety factors
  if( DTCOMP < DTMIN )then
    write(6,800) TIMEDAY, DTTMP, DTMIN, Map2Global(LLOC).IG, Map2Global(LLOC).JG, HP(LLOC)
    write(6,801) Map2Global(L1LOC).IG, Map2Global(L1LOC).JG, DTL1MN, HP(L1LOC)
    write(6,802) Map2Global(L2LOC).IG, Map2Global(L2LOC).JG, DTL2MN, HP(L2LOC)
    if( DTSSDHDT > 0. )then
      write(6,804) Map2Global(L4LOC).IG, Map2Global(L4LOC).JG, DTL4MN
    endif

    open(mpi_log_unit,FILE = OUTDIR//mpi_log_file,POSITION = 'APPEND')
    write(mpi_log_unit,800) TIMEDAY, DTTMP, DTMIN, Map2Global(LLOC).IG, Map2Global(LLOC).JG, HP(LLOC)
    write(mpi_log_unit,801) Map2Global(L1LOC).IG, Map2Global(L1LOC).JG, DTL1MN, HP(L1LOC)
    write(mpi_log_unit,802) Map2Global(L2LOC).IG, Map2Global(L2LOC).JG, DTL2MN, HP(L2LOC)

    if( DTSSDHDT > 0. )then
      write(mpi_log_unit,804) Map2Global(L4LOC).IG, Map2Global(L4LOC).JG, DTL4MN
    endif
    close(mpi_log_unit)
    DTTMP = DTMIN

  elseif( DTTMP < DTMIN )then
    DTWARN = DTTMP
    DTTMP  = DTMIN
  else
    TMPVAL = DTTMP/DTMIN
    ITMPR = NINT(TMPVAL)
    RTMPR = DBLE(ITMPR)
    if( RTMPR < TMPVAL )then
      DTTMP = RTMPR*DTMIN
    else
      DTTMP = (RTMPR-1.)*DTMIN
    endif
  endif

  ! *** FORCE TO MINIMUM TIME STEP ON STARTUP
  if( NITER < 2 ) DTTMP = DTMIN

  ! *** RESTRICT INCREASE IN TIME STEP TO DTMIN
  !DTDYNP = DTHISTORY(10) + DTMIN
  DTDYNP = DTDYN + DTMIN               ! *** AT THIS POINT DTDYN IS THE PREVIOUS ITERATIONS TIME STEP
  if( DTTMP > DTDYNP )then
    NUP = NUP + 1
    if( NUP >= NUPSTEP )then
      DTTMP = DTDYNP
      NUP = 0
    else
      DTTMP = DTDYN
    endif
  elseif( DTTMP < DTDYN )then
    if( (NUPSTEP >= 25 .and. TIMEDAY > TBEGIN+1.) .or. 1.2*DTTMP < DTDYN )then
      if( process_id == master_id )then
        ! *** REPORT WHEN THERE IS A DECREASE
        if( 2.*DTTMP < DTDYN )then
          write(6,'(A,I10,F14.5,2(A,F10.4,I7))') 'DTDYN REDUCED',NITER,TIMEDAY,'     [DT,L]   OLD: ',DTDYN,LLOCOLD,'     NEW: ',DTTMP,LMINSTEP
        endif
        
        open(9,FILE = OUTDIR//'TIME.LOG',POSITION = 'APPEND')
        write(9,'(A,I10,F14.5,2(A,F10.4,I7,F8.3))') 'DYNAMIC STEP DOWN',NITER,TIMEDAY,'     [DT,L,HP]   OLD: ',DTDYN,LLOCOLD,-99.99,'     NEW: ',DTTMP,LMINSTEP,-99.99
        close(9)
      endif
    endif

    NUP = 0
  else
    DTTMP = DTDYN
  endif
  DTDYN = min(DTTMP,DTMAX)

  ! *** Set incremental increase in output counter
  NINCRMTLoc = NINT(DTDYN/DTMIN)

  call MPI_barrier(DSIcomm, IERR)
  call DSI_All_Reduce(NINCRMTLoc, NINCRMT, MPI_MIN, TTDS, 1, TWAIT)

  DTDYN   = DBLE(NINCRMT)*DTMIN

100 FORMAT(5I5,5F12.5,E13.5)
101 FORMAT(3I5,E13.5)
800 FORMAT('  TIME,DTDYN,DTMIN,I,J,HP = ',F12.5,2E12.4,2I7,F10.3)
801 FORMAT('  MOM  ADV,I,J,DT,HP = ',2I7,E13.4,F10.3)
802 FORMAT('  MASS ADV,I,J,DT,HP = ',2I7,E13.4,F10.3)
803 FORMAT('  CURV ACC,I,J,DT,HP = ',2I7,E13.4,F10.3)
804 FORMAT('  LIM DHDT,I,J,DT,HP = ',2I7,E13.4,F10.3)
880 FORMAT(3I5,8E13.4)
8899 FORMAT(' DT3 ERROR ',2I5,6E13.5)

  !DO K = 10,2,-1
  !  DTHISTORY(K) = DTHISTORY(K-1)
  !ENDDO
  !DTHISTORY(1) = DTDYN
     
  ! *** *******************************************************************C
  return
END
