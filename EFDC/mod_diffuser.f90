! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE DIFFUSER_MODULE

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  ! - ---------------------------------------------------------------------!
  !    2018 - 05       PAUL M. CRAIG     RESTRUCTRED EFDC+ DIFFUSER CODE FOR BETTER MANAGEMENT

  USE GLOBAL
  Use Variables_WQ
  
  Use Variables_MPI
  
  IMPLICIT NONE

  REAL :: SALA
  REAL :: SFLA
  REAL :: TEMA

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: UAGD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: VAGD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WAGD

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DYEA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DYEAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SALAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SEDA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SEDAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SFLAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SNDA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SNDAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TEMAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TOXA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WQVA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: WQVAD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: ZAD
  LOGICAL :: LDEBUG
  
  CONTAINS
  
SUBROUTINE JPEFDC

  ! *** PROGRAM JPEFDC IS STAND ALONE VERSION OF EFDC JET - PLUME MODEL
  ! *** BASED ON LEE AND CHEUNG'S LAGRANGIAN BUOYANT JET MODEL
  ! *** AND EXTENDED FOR THREE-DIMENSIONAL AMBIENT CURRENTS
  ! ***   REF: LEE, J.H.W., AND V. CHEUNG, J. ENVIRON. ENGR., 116, 1085 - 
  ! ***   1106, 1990.

  Use Allocate_Initialize
  
  IMPLICIT NONE

  INTEGER,PARAMETER :: NJELM = 2, NATDM = 1
  INTEGER           :: NJP, K, LJP, KJP, KFLAG, ISAMB, ISJPD, ISTOP, ISOUT, NPRTE, NIMAX, ITMP, JTMPVAL, KL
  INTEGER           :: NT, NV, NX, NS, L, LU, KU, LE, LW, LN, LS, N2, NJE, N1, NI, ISHEAR, IFORCE
  INTEGER           :: KTMP, KZERO, IFILE, MD, NW
  INTEGER           :: IDUMP
  INTEGER,SAVE      :: NPRT(10)
  
  REAL :: TIME, RADDEG, ZTMP, ALPHA, ALPH2, TUJP, XERRM, YERRM, ZERRM, RRERM, RKC
  REAL :: RLERM, DMRERM, XJET, YJET, RJET, AJET, QVJET, QSERTAVG, QVJET0, VJET, WTC
  REAL :: WTV, TMPVAL, SALJET, TEMJET, SFLJET, ZBOT, ZSUR, DZPRT, ZJPRT
  REAL :: DTJP, UJ0, VJ0, WJ0, SEDJETT, UTMP, VTMP, TMPEXP, VA, VAH, SEDATT
  REAL :: FRD2I, EBOT, ETOP, ENTS, QJTOT, QJTOTO, RLSCL, SALJETI, TEMJETI, SFLJETI
  REAL :: THAG, SINPHJ, DVEL, DRHO, FTOP, DJETI, XJGNE, YJGNE, ZJGNE, RADJNE, SIGNE
  REAL :: RLEJNE, SALJNE, TEMJNE, SFLJNE, DRMAJ, COSPHJ, COSPHJM, SINTHJL, COSTHJL
  REAL :: COSTHJLM, DRMAJSO, DRMAJSA, DRMAJFO, ENTF1
  REAL :: PTINI, RHO_INI, RHO_ELE , DELTEM, DELRHO, EXITTEM, EXITRHO, RHO_DELTA, RHO_DELTA2
  ! REAL :: XJOLD, YJOLD, ZJOLD   ! DELME
  REAL :: DELSIG, ENTF2, ENTF3, ENTF, DRMAJFA, RMAJI, DXTMP, DYTMP, DZTMP, DS, SEDAA, TOTTIME
  REAL :: DRMAJSE, DRMAJFE, ZJGTOP, ZJGBOT, DRHOT, ZLOWER, ZUPPER, RDUM, QUJ0, QVJ0, QWJ0
  REAL :: QENTTMP, QUAG, QVAG, QWAG, ZVAL, RADLAST, VJET0
  REAL :: UAG, VAG, WAG
  REAL,SAVE :: CURDAY
  
  !REAL,EXTERNAL :: FUNDEN
  
  CHARACTER :: FNNUM*3, PROCNUM*3, FNJPLOG*40

  LOGICAL*4,SAVE :: LOUTJET
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DRHONS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DRMAJF
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DRMAJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DRMAJS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DYEJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DYEJET
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DYEJETI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DYEJNE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DYEJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: PHJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: VJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: VJH
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: QJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RADJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RADJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RADJOLD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RDQA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RDQANS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RHOA
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RHOJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RLEJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RLEJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RMAJP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RMAJPNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SALJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SALJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SEDJET
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SEDJETI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SEDJNE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SFLJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SIG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SIGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SNDJET
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SNDJETI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SNDJNE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TEMJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TEMJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: THJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: THJL
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TIMJP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TOXJET
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TOXJETI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TOXJNE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: UJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: UJGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: VJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: VJGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WJGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: XJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: XJGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: YJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: YJGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: ZJG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: ZJGNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SEDJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SEDJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SNDJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SNDJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:) :: TOXPFJP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXPFTJP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXPFTNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: UJPAVG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: VJPAVG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: WJPAVG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: WQVJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WQVJET
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WQVJETI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: WQVJNE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: WQVJNS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RPORTS
  
  LDEBUG = DEBUG
  !LDEBUG = .TRUE.   ! delme
  
  ! **************************************************************************************************
  ! **************************************************************************************************
  ! *** Allocation Block
  IF(  .NOT. ALLOCATED(DRHONS) )THEN
    WRITE(*,'(A,I8)')'JET/PLUME COMPUTATIONS STARTED.  NQJPIJ = ',NQJPIJ
    
    ! *** MODULE LEVEL VARIABLES FIRST
    Call AllocateDSI( DYEA,     NDYM,   0.0)
    Call AllocateDSI( DYEAD,    KCM,    NDYM,   0.0)
    Call AllocateDSI( SALAD,    KCM,    0.0)    
    Call AllocateDSI( SEDA,     NSED2,  0.0)      
    Call AllocateDSI( SEDAD,    KCM,    NSCM2,  0.0)
    Call AllocateDSI( SFLAD,    KCM,    0.0)    
    Call AllocateDSI( SNDA,     NSND,   0.0)    
    Call AllocateDSI( SNDAD,    KCM,    NSND,   0.0)
    Call AllocateDSI( TEMAD,    KCM,    0.0)    
    Call AllocateDSI( TAD,      KCM,    0.0)    
    Call AllocateDSI( TOXA,     NTXM,   0.0)    
    Call AllocateDSI( TOXAD,    KCM,    NTXM,   0.0)
    Call AllocateDSI( UAGD,     KCM,    0.0)
    Call AllocateDSI( VAGD,     KCM,    0.0)
    Call AllocateDSI( WAGD,     KCM,    0.0)
    Call AllocateDSI( ZAD,      KCM,    0.0)
                      
    Call AllocateDSI( DRHONS,   NJPSM,  0.0)      
    Call AllocateDSI( DRMAJF,   NJELM,  0.0)      
    Call AllocateDSI( DRMAJNS,  NJPSM,  0.0)      
    Call AllocateDSI( DRMAJS,   NJELM,  0.0)      
    Call AllocateDSI( DYEJ,     NJELM,  NDYM,   0.0)
    Call AllocateDSI( DYEJET,   NDYM,   0.0)      
    Call AllocateDSI( DYEJETI,  NDYM,   0.0)      
    Call AllocateDSI( DYEJNE,   NDYM,   0.0)      
    Call AllocateDSI( DYEJNS,   NDYM,   NJPSM,  0.0)
    Call AllocateDSI( PHJ,      NJELM,  0.0)      
    Call AllocateDSI( VJ,       NJELM,  0.0)      
    Call AllocateDSI( VJH,      NJELM,  0.0)      
    Call AllocateDSI( QJNS,     NJPSM,  0.0)      
    Call AllocateDSI( RADJ,     NJELM,  0.0)   
    Call AllocateDSI( RADJNS,   NJPSM,  0.0)      
    Call AllocateDSI( RADJOLD,  NJELM,  0.0)   
    Call AllocateDSI( RDQA,     NJELM,  0.0)      
    Call AllocateDSI( RDQANS,   NJPSM,  0.0)      
    Call AllocateDSI( RHOA,     NJELM,  0.0)      
    Call AllocateDSI( RHOJ,     NJELM,  0.0)      
    Call AllocateDSI( RLEJ,     NJELM,  0.0)      
    Call AllocateDSI( RLEJNS,   NJPSM,  0.0)      
    Call AllocateDSI( RMAJP,   -NJELM,  0.0)      
    Call AllocateDSI( RMAJPNS,  NJPSM,  0.0)      
    Call AllocateDSI( SALJ,     NJELM,  0.0)      
    Call AllocateDSI( SALJNS,   NJPSM,  0.0)      
    Call AllocateDSI( SEDJ,     NJELM,  NSCM,   0.0)
    Call AllocateDSI( SEDJET,   NSCM2,  0.0)      
    Call AllocateDSI( SEDJETI,  NSCM2,  0.0)      
    Call AllocateDSI( SEDJNE,   NSCM2,  0.0)      
    Call AllocateDSI( SEDJNS,   NSCM2,  NJPSM,  0.0)
    Call AllocateDSI( SFLJ,     NJELM,  0.0)      
    Call AllocateDSI( SIG,      NJELM,  0.0)      
    Call AllocateDSI( SIGNS,    NJPSM,  0.0)      
    Call AllocateDSI( SNDJ,     NJELM,  NSNM,   0.0)
    Call AllocateDSI( SNDJET,   NSNM,   0.0)      
    Call AllocateDSI( SNDJETI,  NSNM,   0.0)      
    Call AllocateDSI( SNDJNE,   NSNM,   0.0)      
    Call AllocateDSI( SNDJNS,   NSNM,   NJPSM,  0.0)
    Call AllocateDSI( TEMJ,     NJELM,  0.0)      
    Call AllocateDSI( TEMJNS,   NJPSM,  0.0)      
    Call AllocateDSI( THJG,     NJELM,  0.0)      
    Call AllocateDSI( THJL,     NJELM,  0.0)      
    Call AllocateDSI( TIMJP,    NJPSM,  0.0)      
    Call AllocateDSI( TOXJ,     NJELM,  NTXM,   0.0)
    Call AllocateDSI( TOXJET,   NTXM,   0.0)      
    Call AllocateDSI( TOXJETI,  NTXM,   0.0)      
    Call AllocateDSI( TOXJNE,   NTXM,   0.0)      
    Call AllocateDSI( TOXJNS,   NTXM,   NJPSM,  0.0)
    Call AllocateDSI( TOXPFJP,  NJELM,  NSTM,   NTXM, 0.0)
    Call AllocateDSI( TOXPFTJP, NJELM,  NTXM,   0.0)
    Call AllocateDSI( TOXPFTNS, NTXM,   NJPSM,  0.0)
    Call AllocateDSI( UJG,      NJELM,  0.0)
    Call AllocateDSI( UJGNS,    NJPSM,  0.0)
    Call AllocateDSI( UJPAVG,   KCM,    NJPSM,  0.0)
    Call AllocateDSI( VJG,      NJELM,  0.0)    
    Call AllocateDSI( VJGNS,    NJPSM,  0.0)    
    Call AllocateDSI( VJPAVG,   KCM,    NJPSM,  0.0)
    Call AllocateDSI( WJG,      NJELM,  0.0)    
    Call AllocateDSI( WJGNS,    NJPSM,  0.0)    
    Call AllocateDSI( WJPAVG,   KCM,    NJPSM,  0.0)
    Call AllocateDSI( XJG,      NJELM,  0.0)
    Call AllocateDSI( XJGNS,    NJPSM,  0.0)
    Call AllocateDSI( YJG,      NJELM,  0.0)
    Call AllocateDSI( YJGNS,    NJPSM,  0.0)
    Call AllocateDSI( ZJG,      NJELM,  0.0)
    Call AllocateDSI( ZJGNS,    NJPSM,  0.0)
    Call AllocateDSI( RPORTS,   NQJPIJ, 0.0)

    IF( ISTRAN(8) > 0 )THEN
      Call AllocateDSI( WQVJ,     NJELM, -NWQV,   0.0)
      Call AllocateDSI( WQVA,    -NWQV,   0.0)
      Call AllocateDSI( WQVAD,    KCM,   -NWQV,   0.0)
      Call AllocateDSI( WQVJET,  -NWQV,   0.0)      
      Call AllocateDSI( WQVJETI, -NWQV,   0.0)      
      Call AllocateDSI( WQVJNE,  -NWQV,   0.0)      
      Call AllocateDSI( WQVJNS,  -NWQV,   NJPSM,  0.0)
    ENDIF
    
    ! *** INITIALIZE LOCAL ARRAYS
    NPRT = 0
    CURDAY = TIMEDAY + 1./24.
    RMAJP = 0.0

    LOUTJET = .FALSE.
    DO NJP = 1,NQJPIJ
      IF( IOUTJP(NJP) > 0 ) LOUTJET = .TRUE.
      RPORTS(NJP) = FLOAT(NPORTJP(NJP))  
    ENDDO
    
    IF( LOUTJET )THEN
      OPEN(88,FILE = OUTDIR//'JPBUG.DIA', STATUS = 'UNKNOWN')
      CLOSE(88,STATUS = 'DELETE')
    
      DO NJP = 1,NQJPIJ
        IF( IOUTJP(NJP) > 0 )THEN
          WRITE(FNNUM,'(I3.3)') NJP
          WRITE(PROCNUM,'(I3.3)') process_id
          FNJPLOG = OUTDIR//'JPLOG' // FNNUM // '_' // PROCNUM // '.OUT'
          OPEN(10,FILE = FNJPLOG, STATUS = 'UNKNOWN')
          CLOSE(10,STATUS = 'DELETE')
          
          IF( IOUTJP(NJP) > 1 )THEN
            OPEN(10,FILE = FNJPLOG, STATUS = 'UNKNOWN')
            WRITE(10,630) 'NITER', 'TIMEDAY', 'NJE', 'NI', 'TOTTIME', 'LENGTH',                                                          &
                          'XJG(', 'YJG', 'ZJG', 'DIA', 'QVJET', 'SUM QJPENT', 'QJTOT', 'RHOA', 'RHOJ', 'TEMJ', 'VJ', 'VJH', 'WJG',       &
                          'WAG', 'RHOJ-RHOA', 'RHO_DELTA', 'RHO_DELTA2', 'TEMA', 'VA', 'VAH', 'RMAJP', 'DRMAJS', 'DRMAJF'
            CLOSE(10)
          ENDIF
        ENDIF
      ENDDO
630   FORMAT(A10,A15,A8,A5,2A9,2A12,A10,A8,3A12,8A10,10A12)  
      
    ENDIF

    IF( LDEBUG )THEN
      OPEN(88,FILE = OUTDIR//'JPBUG.DIA', STATUS = 'UNKNOWN')
      CLOSE(88,STATUS = 'DELETE')
    ENDIF
    
  ENDIF   ! *** END OF INITIALIZATION

  ! *** SET CONSTANTS
  IF( ISDYNSTP == 0 )THEN
    TIME = DT*FLOAT(N) + TCON*TBEGIN
    TIME = TIME/TCON
  ELSE
    TIME = TIMESEC/TCON
  ENDIF
  
  IF( CURDAY < TIMEDAY )THEN
    CURDAY = CURDAY + 1./24.
    NPRT = 0
  ENDIF

  ! *** G   = EFDC+ global constant for acceleration due to gravity
  ! *** RPI = 3.14159
  RADDEG = PI/180.
  
  IFILE = -1
  IF( LDEBUG )THEN
    OPEN(88,FILE = OUTDIR//'JPBUG.DIA',STATUS = 'UNKNOWN',POSITION = 'APPEND')
    CLOSE(88,STATUS = 'DELETE')
    WRITE(88,*) TIMEDAY
  ENDIF
  
  ! **************************************************************************************************
  ! **************************************************************************************************
  ! *** Prepare for plume calculations

  ! *** LOOP OVER ALL JET/PLUME LOCATIONS
  DO NJP = 1,NQJPIJ
    ! *** Inititalize and prepare for the plume calculatopn for the current jet NJP
    DO K = 1,KC
      QJPENT(K,NJP) = 0.0
      UJPAVG(K,NJP) = 0.0
      VJPAVG(K,NJP) = 0.0
      WJPAVG(K,NJP) = 0.0
    ENDDO

    IF( NITER == 1 )THEN
      KEFFJP(NJP) = KQJP(NJP)               ! *** Initialize the plume discharge elevation to the nominal elevation
    ENDIF
    
    IF( IOUTJP(NJP) > 0 )THEN
      WRITE(FNNUM,'(I3.3)') NJP
      WRITE(PROCNUM,'(I3.3)') process_id
      FNJPLOG = OUTDIR//'JPLOG' // FNNUM // '_' // PROCNUM // '.OUT'
      OPEN(10,FILE = FNJPLOG,STATUS = 'UNKNOWN',POSITION = 'APPEND')
      IF( IOUTJP(NJP) > 1 ) WRITE(10,134) NJP, TIMEDAY
    ENDIF
    
    IF( LDEBUG )THEN
      OPEN(88,FILE = OUTDIR//'JPBUG.DIA', STATUS = 'UNKNOWN',POSITION = 'APPEND')      
      WRITE(88,134) NJP, TIMEDAY
    ENDIF
    
    LJP   = LIJ(IQJP(NJP),JQJP(NJP))         ! *** L index of the receiving/returning water
    KL    = KSZ(LJP)                         ! *** Current layer
    RKC   = FLOAT(KC - KL + 1)               ! *** Number of layers from surface
    KJP   = KQJP(NJP)                        ! *** Default/initial discharge layer
    KFLAG = 0
    
    ZTMP = (ZJET(NJP) - BELV(LJP))/HP(LJP)   ! *** Fraction of the water column from the bottom to the port elevation
    KJP  = NINT(RKC*ZTMP)
    IF( KJP < KSZ(LJP) ) KJP = KSZ(LJP)
    IF( KJP > KC ) KJP = KC
    
    ! *** ICAL > 0 - JET/PLUME BC ACTIVATED, OTHERWISE IGNORED
    IF( ICALJP(NJP) > 0 )THEN
      !    NJEL     MAXIMUM NUMBER OF ELEMENTS ALONG JET/PLUME LENGTH
      !    ISAMB  0 FOR SPATIALLY AND TEMPORALLY CONSTANT AMBIENT CONDITIONS
      !           1 FOR SPATIALLY VARYING AND TEMPORALLY CONSTANT CONDITIONS
      !           2 FOR SPATIALLY AND TEMPORALLY VARYING AMBIENT CONDITIONS
      !    ISJPD  0 FOR TEMPORALLY CONSTANT JET/PLUME DISCHARGE
      !           1 FOR TEMPORALLY VARYING JET/PLUME DISCHARGE
      !    ISENT  0 USE MAXIMUM OF SHEAR AND FORCED ENTRAINMENT
      !           1 USE SUM OF SHEAR AND FORCED ENTRAINMENT
      !    ISTOP  0 STOP AT SPECIFIED NUMBER OF ELEMENTS
      !           1 STOP WHEN CENTERLINE PENETRATES BOTTOM OR SURFACE
      !           2 STOP WITH BOUNDARY PENETRATES BOTTOM OR SURFACE
      !    ISOUT  0 DIMENSIONAL OUTPUT,
      !           1 NONDIM OUTPUT LENGTH SCALED BY DJET
      !           2 NONDIM OUTPUT LENGTH SCALED BY SQRT(PI)*RJET
      ISAMB = 1         
      ISJPD = 0         
      ISTOP = ISTJP(NJP)
      ISOUT = 0         
      NPRTE = 0              ! *** NPRTE    Element output print frequency
      ALPHA = CFRD(NJP)      ! *** CFRD     Adjustment factor for froude number
      ALPH2 = ALPHA*ALPHA
      TUJP = 0               ! *** TUJP     Temporal frequency for updating jet/plume (sec) if ISJPD = 1
      
      NIMAX = NJPMX(NJP)     ! *** NIMAX    Maximum numeber of iteration
      XERRM = 1000.          ! *** XYERR    Horizontal trajectory error critera (m)
      YERRM = 1000.           
      ZERRM = 1000.          ! *** ZERR     Vertical trajectory error critera (m)
      RRERM = 1000.          ! *** RRER     Horizontal trajectory error critera (m)
      RLERM = 1000.          ! *** RLER     Vertical trajectory error critera (m)
      DMRERM = DJPER(NJP)    ! *** DMRERR   Entrainment error criteria
      
      ! ***    TDCY   Contaminant decay rate (1./sec)
      ! ***    WSET,TPAR,TDCY not used in embedded version.  Appropriate efdc variables are used instead
      ! ***    XJET   East coordinate of discharge (m)   (not used)
      ! ***    YJET   North coordinate of discharge (m)  (not used)
      ! ***    ZJET   Elevation of discharge (m)
      ! ***    PHJET  Vertical jet angle positive from horizontal (degrees)
      ! ***    THJET  Horizontal jet angle pos counter clockwise from east (degrees)
      ! ***    AJET   Area of the jet (m2)
      ! ***    VJET   Constant port (jet at port) velocity (m/s)   (negative value indicates withdrawal)
      ! ***    QQCJP  Constant discharge rate (cms)                (negative value indicates withdrawal)
      ! ***    QVJET  Total flow (constant + time variable)
      ! *** QSERTAVG  Total time variable flow
      
      XJET = 0.               ! *** WSET   Sediment settling velocity (m/s)
      YJET = 0.               ! *** TPAR   Partition coefficient (l/mg or m**3/gm)
      RJET = 0.5*DJET(NJP)    ! *** RJET   Radius of discharge port (m),   
                              ! *** DJET   Diameter of discharge port (m)
      AJET = PI*RJET*RJET     ! *** This is area of 1 port but the flows used below are total flow into the cell

      ! *** Total Flow over all layers (QSER Type, ICAL = 1 )
      IF( ICALJP(NJP) == 1 )THEN
        QVJET = QQCJP(NJP)      
        QSERTAVG = 0.
        DO K = 1,KC
          QVJET =    QVJET    + QSERT(K,NQSERJP(NJP))
          QSERTAVG = QSERTAVG + QSERT(K,NQSERJP(NJP))
        ENDDO
      ENDIF

      ! *** Total Flow over all layers (Withdrawal/Return type, ICAL = 2 )
      IF( ICALJP(NJP) == 2 )THEN
        QVJET = QWRCJP(NJP) + QWRSERT(NQWRSERJP(NJP))
      ENDIF
      IF( QVJET <= 1.E-16 ) GOTO 9000   ! *** NO FLOW, SKIP CALCS
      
      ! *** THE QSER OR NQWRSERJP FLOWS ARE ON A PER PORT BASIS.  CALQVS ADJUSTS PORT FLOWS FOR NUMBER OF PORTS IN CELL.  [ ALTERNATE: QVJET = QVJET/RPORTS(NJP) ]
      QVJET0 = QVJET      
      VJET = QVJET/AJET   ! *** Velocity of jet in port (m/s)
      
      ! *** Get jet/plume mass loading from QSER or W/R 
      !    SALJET   SALINITY (PPT)
      !    SEDJET   SEDIMENT (MG/L OR G/M^3)
      !    TEMJET   TEMPERATURE
      !    TOXJET   CONTAMINANT CONCENTRATION (UG/L OR MG/M**3)
      !    WQVJET   WATER QUALITY CONCENTRATION (G/M3)
      
      ! *** STANDARD FLOW (QSER) TYPE (ICAL = 1)
      IF( ICALJP(NJP) == 1 )THEN
        WTC = QQCJP(NJP)              ! *** Constant flow component
        WTV = QSERTAVG                ! *** Variable flow component
        TMPVAL = WTC + WTV
        WTC = WTC/TMPVAL              ! *** Fraction of total flow for constant flows
        WTV = WTV/TMPVAL              ! *** Fraction of total flow for variable flows

        ! ***      Constant                   Series
        SALJET = WTC*CQCJP(1,NJP,1) + WTV*CSERT(1,NCSERJP(NJP,1),1)
        TEMJET = WTC*CQCJP(1,NJP,2) + WTV*CSERT(1,NCSERJP(NJP,2),2)

        DO MD = 1,NDYE
          NV = MSVDYE(MD)
          DYEJET(MD) = WTC*CQCJP(1,NJP,NV) + WTV*CSERT(1,NCSERJP(NJP,3),NV)
        ENDDO
        
        NV = 3 + NDYM
        SFLJET = WTC*CQCJP(1,NJP,4) + WTV*CSERT(1,NCSERJP(NJP,4),NV)

        DO NT = 1,NTOX
          NV = MSVTOX(NT)
          TOXJET(NT) = WTC*CQCJP(1,NJP,NV) + WTV*CSERT(1,NCSERJP(NJP,5),NV)
        ENDDO
        DO NX = 1,NSED     
          NV = MSVSED(NX)
          SEDJET(NX) = WTC*CQCJP(1,NJP,NV) + WTV*CSERT(1,NCSERJP(NJP,6),NV)
        ENDDO
        DO NX = 1,NSND
          NV = MSVSND(NX)
          SNDJET(NX) = WTC*CQCJP(1,NJP,NV) + WTV*CSERT(1,NCSERJP(NJP,7),NV)
        ENDDO
        DO NX = 1,NWQV
          NV = MSVWQV(NX)
          WQVJET(NX) = WTC*CQCJP(1,NJP,NV) + WTV*CSERT(1,NCSERJP(NJP,8),NV)
        ENDDO
      ENDIF
      
      ! *** WITHDRAWAL/RETURN TYPE (ICAL = 2)
      IF( ICALJP(NJP) == 2 )THEN
        WTC = QWRCJP(NJP)                  ! *** Current withdrawal/return flow
        WTV = QWRSERT(NQWRSERJP(NJP))      ! *** Times series number
        TMPVAL = WTC + WTV
        WTC = WTC/TMPVAL                   ! *** Fraction of total flow for constant flows
        WTV = WTV/TMPVAL                   ! *** Fraction of total flow for variable flows
        NS = NQWRSERJP(NJP)
        LU = LIJ(IUPCJP(NJP),JUPCJP(NJP))  ! *** L Index at the withdrawal cell
        KU = KUPCJP(NJP)                   ! *** Layer at withdrawal cell
        
        ! ***      Constant                   Series
        SALJET = WTC*CWRCJP(NJP,1) + WTV*CQWRSERT(NS,1) + SAL1(LU,KU)
        TEMJET = WTC*CWRCJP(NJP,2) + WTV*CQWRSERT(NS,2) + TEM1(LU,KU)
        
        DO MD = 1,NDYE
          NV = MSVDYE(MD)
          DYEJET(MD) = WTC*CWRCJP(NJP,NV) + WTV*CQWRSERT(NS,NV) + DYE1(LU,KU,MD)
        ENDDO

        NV = 3 + NDYM
        SFLJET = WTC*CWRCJP(NJP,NV) + WTV*CQWRSERT(NS,NV) + SFL2(LU,KU)

        DO NT = 1,NTOX
          NV = MSVTOX(NT)
          TOXJET(NT) = WTC*CWRCJP(NJP,NV) + WTV*CQWRSERT(NS,NV) + TOX1(LU,KU,NT)
        ENDDO
        DO NX = 1,NSED
          NV = MSVSED(NX)
          SEDJET(NX) = WTC*CWRCJP(NJP,NV) + WTV*CQWRSERT(NS,NV) + SED1(LU,KU,NX)
        ENDDO
        DO NX = 1,NSND
          NV = MSVSND(NX)
          SNDJET(NX) = WTC*CWRCJP(NJP,NV) + WTV*CQWRSERT(NS,NV) + SND1(LU,KU,NX)
        ENDDO
        DO NW = 1,NWQV
          NV = MSVWQV(NW)
          WQVJET(NW) = WTC*CWRCJP(NJP,NV) + WTV*CQWRSERT(NS,NV) + WQV(LU,KU,NW)
        ENDDO
      ENDIF
      
      ZBOT = BELV(LJP)              ! *** ZBOT   BOTTOM ELEVATION (M)
      ZSUR = BELV(LJP) + HP(LJP)    ! *** ZSUR   SURFACE ELEVATION (M)
                                    ! *** KC     NUMBER OF AMBIENT DATA IN VERTICAL
      NATD = 1                      ! *** NATD   NUMBER OF AMBIENT DATA IN TIME

      !    SPATIAL PRINT/SAVE INTERVAL IN VERTICAL
      DZPRT = HP(LJP)/FLOAT(NZPRJP(NJP) - 2)
      ZJPRT = ZJET(NJP) + DZPRT

      ! *** ZA     Elevation of ambient data (m)
      ! *** TA     Time of ambient data (sec, hr, or day for output reference)
      ! *** UAG    East velocity (m/s)
      ! *** VAG    North velocity (m/s)
      ! *** WAG    Vertical velocity (m/s)
      ! *** SALA   Ambient salinity (ppt)
      ! *** SEDA   Ambient sediment concentration (mg/l or gm/m**3)
      ! *** TEMA   Ambient temperature (deg c)
      ! *** TOXA   Ambient toxic contaminant concentration (ug/l or mg/m**3)
      
      ! *** ZAD is the layer midpoint elevation for each layer
      ZAD(1) = BELV(LJP) + 0.5*HPK(LJP,KSZ(LJP))
      DO K = KSZ(LJP)+1,KC
        ZAD(K) = ZAD(K-1) + 0.5*(HPK(LJP,K-1) + HPK(LJP,K))
      ENDDO
      
      ! *** Get layer specific cell centroid velocities
      DO K = 1,KC
        L = LJP
        LE = LEC(L)
        LW = LWC(L)
        LN = LNC(L)
        LS = LSC(L)
        UTMP = 0.5*(U(LE,K) + U(L,K))                   ! *** Cell centroid U velocity
        VTMP = 0.5*(V(LN,K) + V(L,K))                   ! *** Cell centroid V velocity
        UAGD(K) = CUE(L)*UTMP + CVE(L)*VTMP             ! *** Rotating to world coordinates
        VAGD(K) = CUN(L)*UTMP + CVN(L)*VTMP             ! *** Rotating to world coordinates
        
        ! *** Vertical velocities must be computed because W includes artificial velocities to account for plume entrainment.
        ! *** For a rising or falling plume, the entrainment flows come from the sides, not vertically         
        !WAGD(K) = 0.5*(W(L,K) + W(L,K - 1)) +    *** Revised                                                         &
        WAGD(K) =  GI*ZZ(L,K)*( DTI*(P(L) - P1(L)) +                                                     &
                               0.5*( U(LE,K)*(P(LE)-P(L))*DXIU(LE) + U(L,K)*(P(L)-P(LW))*DXIU(L) +       &
                                     V(LN,K)*(P(LN)-P(L))*DYIV(LN) + V(L,K)*(P(L)-P(LS))*DYIV(L) ) ) +   &
                  0.5*(1.-ZZ(L,K))*( U(LE,K)*(BELV(LE)-BELV(L)) *DXIU(LE) +                              &
                                     U(L,K) *(BELV(L) -BELV(LW))*DXIU(L)  +                              &
                                     V(LN,K)*(BELV(LN)-BELV(L)) *DYIV(LN) +                              &
                                     V(L,K)* (BELV(L) -BELV(LS))*DYIV(L) )
      ENDDO
889   FORMAT(5I6,3E14.5)
      
      ! *** Initialize cell water column concentrations (xxxAD) for receiving/returning water cell
      DO K = 1,KC
        SALAD(K) = SAL(LJP,K)
        TEMAD(K) = TEM(LJP,K)
        SFLAD(K) = SFL(LJP,K)
      ENDDO
      DO MD = 1,NDYE
        DO K = 1,KC
          DYEAD(K,NT) = DYE(LJP,K,MD)
        ENDDO
      ENDDO
      DO NT = 1,NTOX
        DO K = 1,KC
          TOXAD(K,NT) = TOX(LJP,K,NT)
        ENDDO
      ENDDO
      DO NS = 1,NSED2
        DO K = 1,KC
          SEDAD(K,NS) = SED(LJP,K,NS)
        ENDDO
      ENDDO
      DO NX = 1,NSND
        DO K = 1,KC
          SNDAD(K,NX) = SND(LJP,K,NX)
        ENDDO
      ENDDO
      DO NW = 1,NWQV
        DO K = 1,KC
          WQVAD(K,NT) = WQV(LJP,K,NW)
        ENDDO
      ENDDO  

      ! *** Set initial conditions for current plume
      DTJP = 0.1*RJET/VJET    ! *** 1/S

      ! ***  JET/PLUME DISCHARGE, GLOBAL COORDINATES
      XJG(1) = XJET          ! *** X coordinate of starting location   - Element 1
      YJG(1) = YJET          ! *** Y coordinate of starting location   - Element 1
      ZJG(1) = ZJET(NJP)     ! *** Elevation of starting location      - Element 1
      XJG(2) = XJET          ! *** X coordinate of starting location   - Element 2
      YJG(2) = YJET          ! *** Y coordinate of starting location   - Element 2
      ZJG(2) = ZJET(NJP)     ! *** Elevation of starting location      - Element 2
      SIG(1) = 0.            ! *** Total plume centerline displacement - Element 1
      SIG(2) = 0.            ! *** Total plume centerline displacement - Element 2

      ! ***  INITIALIZE JET DISCHARGE GEOMENTRY FOR THE COMPUTATIONAL ELEMENTS
      RADJ(1) = RJET         ! *** Radius of the jet - Element 1
      RADJ(2) = RJET         ! *** Radius of the jet - Element 2
      RADJOLD(1) = RJET      ! *** Radius of the jet of previous iteration - Element 1
      RADJOLD(2) = RJET      ! *** Radius of the jet of previous iteration - Element 2
      RLEJ(1) = 0.1*RJET     ! *** 0.1*Radius = Element length of the plume - Element 1
      RLEJ(2) = 0.1*RJET     ! *** 0.1*Radius = Element length of the plume - Element 2
      RADLAST = RJET         ! *** Initializes write tolerance
      
      ! *** JET DISCHARGE VELOCITY MAGNITUDE AND COMPONENTS
      ! *** PHJET IS THE VERTICAL ANGLE FROM HORIZONTAL
      ! *** THJET IS THE HORIZONTAL ANGLE
      AJET = PI*RJET*RJET

      ! *** INITIALIZE ELEMENT 1 VELOCITIES
      VJ(1)  = QVJET/AJET                       ! *** Jet centerline velocity - Total
      VJH(1) = VJ(1)*COS(0.0175*PHJET(NJP))     ! *** Jet centerline velocity - Horizontal only
      UJG(1) = VJ(1)*COS(0.0175*PHJET(NJP))*COS(0.0175*THJET(NJP))
      VJG(1) = VJ(1)*COS(0.0175*PHJET(NJP))*SIN(0.0175*THJET(NJP))
      WJG(1) = VJ(1)*SIN(0.0175*PHJET(NJP))
      UJ0   = UJG(1)
      VJ0   = VJG(1)
      WJ0   = WJG(1)
      VJET0 = VJ(1)
            
      ! *** INITIALIZE ELEMENT 2 VELOCITIES
      VJ(2) = QVJET/AJET
      VJH(2) = VJ(2)*COS(0.0175*PHJET(NJP))
      UJG(2) = VJ(2)*COS(0.0175*PHJET(NJP))*COS(0.0175*THJET(NJP))
      VJG(2) = VJ(2)*COS(0.0175*PHJET(NJP))*SIN(0.0175*THJET(NJP))
      WJG(2) = VJ(2)*SIN(0.0175*PHJET(NJP))
      
      DTJP = RLEJ(1)/VJ(1)                    ! *** Initialize time step (s)
      TOTTIME = 0.                            ! *** Initialize total elapsed time of plume development

      ! ***  INITIAL JET DENSITY AND MASS
      SEDJETT = 0.
      DO NS = 1,NSED
        SEDJETT = SEDJETT + SEDJET(NS)
      ENDDO
      DO NX = 1,NSND
        SEDJETT = SEDJETT + SNDJET(NX)
      ENDDO
      RHOJ(1) = FUNDEN(SALJET,SEDJETT,TEMJET)         ! *** Jet density
      RMAJP(1) = PI*RHOJ(1)*RADJ(1)*RADJ(1)*RLEJ(1)   ! *** Jet mass

      ! ***  INITIAL JET CONCENTRATIONS
      SALJ(1) = SALJET
      TEMJ(1) = TEMJET
      DO MD = 1,NDYE
        DYEJ(1,MD) = DYEJET(MD)
      ENDDO
      SFLJ(1) = SFLJET
      DO NT = 1,NTOX
        TOXJ(1,NT) = TOXJET(NT)
      ENDDO
      DO NS = 1,NSED
        SEDJ(1,NS) = SEDJET(NS)
      ENDDO
      DO NX = 1,NSND
        SNDJ(1,NX) = SNDJET(NX)
      ENDDO
      DO NW = 1,NWQV
        WQVJ(1,NW) = WQVJET(NW)
      ENDDO

      ! ***  INITIAL JET TOXIC CONTAMINANT PARTICULATE FRACTIONS
      DO NT = 1,NTOX
        IF( ISTRAN(6) >= 1 )THEN
          DO NS = 1,NSED
            TMPEXP = CONPARW(NS,NT)
            IF( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
            IF( ITXPARW(NS,NT) == 1 )THEN
              IF( SEDJET(NS) > 0.) TMPVAL = SEDJET(NX)**TMPEXP
            ENDIF
            TOXPFJP(1,NS,NT) = TMPVAL*SEDJET(NS)*TOXPARW(LJP,NS,NT)
          ENDDO
        ENDIF
        IF( ISTRAN(7) >= 1 )THEN
          DO NX = 1,NSND
            NS = NX + NSED
            TMPEXP = CONPARW(NS,NT)
            IF( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
            IF( ITXPARW(NS,NT) == 1 )THEN
              IF( SNDJET(NX) > 0.) TMPVAL = SNDJET(NX)**TMPEXP
            ENDIF
            TOXPFJP(1,NS,NT) = TMPVAL*SNDJET(NX)*TOXPARW(LJP,NS,NT)
          ENDDO
        ENDIF
      ENDDO
      DO NT = 1,NTOX
        TOXPFTJP(1,NT) = 0.
      ENDDO
      DO NT = 1,NTOX
        IF( ISTRAN(6) >= 1 )THEN
          DO NS = 1,NSED
            TOXPFTJP(1,NT) = TOXPFTJP(1,NT) + TOXPFJP(1,NS,NT)
          ENDDO
        ENDIF
        IF( ISTRAN(7) >= 1 )THEN
          DO NX = 1,NSND
            NS = NX + NSED
            TOXPFTJP(1,NT) = TOXPFTJP(1,NT) + TOXPFJP(1,NS,NT)
          ENDDO
        ENDIF
      ENDDO
      DO NT = 1,NTOX
        DO NS = 1,NSED + NSND
          TOXPFJP(1,NS,NT) = TOXPFJP(1,NS,NT)/(1. + TOXPFTJP(1,NT))
        ENDDO
      ENDDO
      DO NT = 1,NTOX
        TOXPFTJP(1,NT) = TOXPFTJP(1,NT)/(1. + TOXPFTJP(1,NT))
      ENDDO

      ! ***  INITIALIZE JET ELEMENT 2 VALUE TO ELEMENT 1 VALUES
      RHOJ(2) = FUNDEN(SALJET,SEDJETT,TEMJET)            ! *** Density 
      RMAJP(2) = PI*RHOJ(2)*RADJ(2)*RADJ(2)*RLEJ(2)      ! *** Mass (density * volume)
      SALJ(2) = SALJ(1)
      TEMJ(2) = TEMJ(1)
      RHO_INI = RHOJ(1)                                  ! *** Initial plume density
      PTINI = TEMJ(1)                                    ! *** Initial plume temperature
      DO MD = 1,NDYE
        DYEJ(2,MD) = DYEJ(1,MD)
      ENDDO
      SFLJ(2) = SFLJ(1)
      DO NT = 1,NTOX
        TOXJ(2,NT) = TOXJ(1,NT)
      ENDDO
      DO NT = 1,NTOX
        DO NS = 1,NSED + NSND
          TOXPFJP(2,NS,NT) = TOXPFJP(1,NS,NT)
        ENDDO
      ENDDO
      DO NT = 1,NTOX
        TOXPFTJP(2,NT) = TOXPFTJP(1,NT)
      ENDDO
      DO NS = 1,NSED
        SEDJ(2,NS) = SEDJ(1,NS)
      ENDDO
      DO NX = 1,NSND
        SNDJ(2,NX) = SNDJ(1,NX)
      ENDDO
      DO NW = 1,NWQV
        WQVJ(2,NW) = WQVJ(1,NW)
      ENDDO

      ! ***  INITIALIZE AMBIENT CONDITIONS
      IF( ISAMB == 0 )THEN
        ! *** CONSTANT AMBIENT CONDITIONS (LAYER = KQJP(NJP))
        K = KQJP(NJP)
        UAG = UAGD(K)
        VAG = VAGD(K)
        WAG = WAGD(K)
        SALA = SALAD(K)
        TEMA = TEMAD(K)
        DO MD = 1,NDYE
          DYEA(MD) = DYEAD(K,MD)
        ENDDO
        SFLA = SFLAD(K)
        DO NT = 1,NTOX
          TOXA(NT) = TOXAD(K,NT)
        ENDDO
        DO NS = 1,NSED2
          SEDA(NS) = SEDAD(K,NS)
        ENDDO
        DO NX = 1,NSND
          SNDA(NX) = SNDAD(K,NX)
        ENDDO
        DO NW = 1,NWQV
          WQVA(NW) = WQVAD(K,NW)
        ENDDO
      ELSEIF( ISAMB >= 1 )THEN
        ! *** SPATIALLY VARIABLE AMBIENT CONDITIONS
        ZVAL = ZJG(1)
        CALL JPACON(KSZ(LJP),ZVAL,UAG,VAG,WAG)
      ENDIF

      ! ***  AMBIENT VELOCITY MAGNITUDES
      VA = SQRT(UAG*UAG + VAG*VAG + WAG*WAG)
      VAH = SQRT(UAG*UAG + VAG*VAG)

      ! ***  AMBIENT DENSITY
      SEDATT = 0.
      DO NS = 1,NSED2
        SEDATT = SEDATT + SEDA(NS)
      ENDDO
      DO NX = 1,NSND
        SEDATT = SEDATT + SNDA(NX)
      ENDDO
      RHOA(1) = FUNDEN(SALA,SEDATT,TEMA)
      RHOA(2) = RHOA(1)
      RHO_ELE = RHOA(1)
      DELTEM = PTINI - TEMA
      DELRHO = RHO_INI - RHO_ELE
      EXITTEM = 0.05*ABS(DELTEM)          ! *** Exit criteria based on initial temperature delta
      EXITRHO = 0.05*ABS(DELRHO)          ! *** Exit criteria based on initial density delta

      ! ***  Global Ambient and Jet orientations
      PHJ(1) = PHJET(NJP)
      THJG(1) = THJET(NJP)
      PHJ(2) = PHJET(NJP)
      THJG(2) = THJET(NJP)
      THAG = 57.2958*ATAN2(VAG,UAG)

      ! *** Local jet horizontal orientation
      THJL(1) = THJG(1) - THAG            ! *** Jet - ambient
      THJL(2) = THJG(2) - THAG            ! *** Jet - ambient

      ! *** Projection of ambient current to jet direction
      RDQA(1) = COS(0.0175*PHJ(1))*COS(0.0175*THJL(1))*VAH       ! *** VAH is only horizontal flow
      RDQA(2) = RDQA(1)                                          ! *** VAH is only horizontal flow

      ! *** Initial shear entrainment
      DVEL = ABS(VJ(1) - RDQA(1))                                ! *** Difference in velocities (m/s)
      DRHO = (RHOA(1) - RHOJ(1))/RHOA(1)
      FTOP = G*ABS(DRHO)*RADJ(1)
      FRD2I = 0.
      EBOT = 1.
      ETOP = 0.057
      ENTS = 0.
      IF( DVEL > 0. )THEN
        SINPHJ = SIN(0.0175*PHJ(1))
        FRD2I = FTOP/(ALPH2*DVEL*DVEL)                         ! *** Relative density Froude Number
        EBOT = 1. + 5.*RDQA(1)/DVEL
        ETOP = 0.057 + 0.554*SINPHJ*FRD2I
        ENTS = 1.414*ETOP/EBOT
      ENDIF
      
      ! ***           s      kg/m3     -   m/s   m       m
      DRMAJS(1) = 2.*DTJP*PI*RHOA(1)*ENTS*DVEL*RADJ(1)*RLEJ(1) ! *** Initialize shear entranment - Element 1
      DRMAJS(2) = DRMAJS(1)                                    ! *** Initialize shear entranment - Element 2
      
      DRMAJF(1) = 0.                                           ! *** Initialize forced entranment - Element 1
      DRMAJF(2) = 0.                                           ! *** Initialize forced entranment - Element 2
      DRMAJFO = 0.5*(DRMAJF(2) + DRMAJF(1))

      ! *** INITIALIZE VOLUMETRIC ENTRAINMENT
      QJTOT = QVJET
      QJTOTO = QVJET

      ! *** OUTPUT INITIAL CONDITIONS
      RLSCL = 1.
      IF( ISOUT == 1 ) RLSCL = DJET(NJP)
      IF( ISOUT == 2 ) RLSCL = RJET*SQRT(PI)
      SALJETI = 1.
      TEMJETI = 1.
      DO MD = 1,NDYE
        DYEJETI(MD) = 1.
      ENDDO
      SFLJETI = 1.
      DO NT = 1,NTOX
        TOXJETI(NT) = 1.
      ENDDO
      DO NS = 1,NSED
        SEDJETI(NS) = 1.
      ENDDO
      DO NX = 1,NSND
        SNDJETI(NX) = 1.
      ENDDO
      DO NW = 1,NWQV
        WQVJETI(NW) = 1.
      ENDDO
      
      ! *** APPLY SCALING FACTOR 
      IF( ISOUT >= 1 )THEN
        IF( SALJET > 0.) SALJETI = 1./SALJET
        IF( TEMJET > 0.) TEMJETI = 1./TEMJET
        DO MD = 1,NDYE
          IF( DYEJET(MD) > 0.) DYEJETI(MD) = 1./DYEJET(MD)
        ENDDO
        IF( SFLJET > 0.) SFLJETI = 1./SFLJET
        DO NT = 1,NTOX
          IF( TOXJET(NT) > 0.) TOXJETI(NT) = 1./TOXJET(NT)
        ENDDO
        DO NS = 1,NSED
          IF( SEDJET(NS) > 0.) SEDJETI(NS) = 1./SEDJET(NS)
        ENDDO
        DO NX = 1,NSND
          IF( SNDJET(NX) > 0.) SNDJETI(NX) = 1./SNDJET(NX)
        ENDDO
        DO NW = 1,NWQV
          IF( WQVJET(NW) > 0.) WQVJETI(NW) = 1./WQVJET(NW)
        ENDDO
      ENDIF
      N2 = 1
      NJE = 1
      DJETI = 1./RLSCL
      RHO_DELTA = 999.                  ! *** Initialize exit criteria
      RHO_DELTA2 = -999.                ! *** Initialize exit criteria
      
      XJGNE = DJETI*XJG(N2)  ! ZERO
      YJGNE = DJETI*YJG(N2)  ! ZERO
      ZJGNE = DJETI*ZJG(N2)
      SIGNE = DJETI*SIG(N2)  ! ZERO
      RADJNE = DJETI*RADJ(N2)
      RLEJNE = DJETI*RLEJ(N2)
      SALJNE = SALJETI*SALJ(N2)
      TEMJNE = TEMJETI*TEMJ(N2)
      DO MD = 1,NDYE
        DYEJNE(MD) = DYEJETI(MD)*DYEJ(N2,MD)
      ENDDO
      SFLJNE = TEMJETI*SFLJ(N2)
      DO NT = 1,NTOX
        TOXJNE(NT) = TOXJETI(NT)*TOXJ(N2,NT)
      ENDDO
      DO NS = 1,NSED
        SEDJNE(NS) = SEDJETI(NS)*SEDJ(N2,NS)
      ENDDO
      DO NX = 1,NSND
        SNDJNE(NX) = SNDJETI(NX)*SNDJ(N2,NX)
      ENDDO
      DO NW = 1,NWQV
        WQVJNE(NW) = WQVJETI(NW)*WQVJ(N2,NW)
      ENDDO
      
      DRMAJ = RMAJP(2) - RMAJP(1)                ! *** Mass of newly entrained water at the current iteration
      DRHO = (RHOA(1) - RHOJ(N2))/RHOA(1)        ! *** Relative density diffference between ambient conditions and plume

      ! **************************************************************************************************
      ! **************************************************************************************************
      ! *** Plume entrainment loop for NJEL elements until either max elements or plume hits end point
      DO NJE = 2,NJEL(NJP)

        N2 = 2   ! *** Time element at end of current element
        N1 = 1   ! *** Time element at start of current element
        
        ! *** Setup plume conditions and trajectory
        SINPHJ = SIN(0.0175*PHJ(N1))        ! *** 
        COSPHJ = COS(0.0175*PHJ(N1))        ! *** 
        COSPHJM = COS(0.0175*PHJ(N1))       ! *** 
        SINTHJL = SIN(0.0175*THJL(N1))      ! *** 
        COSTHJL = COS(0.0175*THJL(N1))      ! *** 
        COSTHJLM = COS(0.0175*THJL(N1))     ! *** 
        DTJP = RLEJ(N1)/VJ(N1)              ! *** Current time step (element length/velocity) (s)
        
        TOTTIME = TOTTIME + DTJP
        
        ! *** Loop over mass change due to entrainment for the current element
        NI = 1
1000    CONTINUE

        ! *** Calculate shear entrainment
        DRMAJSO = 0.5*(DRMAJS(N2) + DRMAJS(N1))      ! *** Average shear between N2 and N1
        DVEL = ABS(VJ(N2) - RDQA(N2))                ! *** Velocity shear (m/s):  VJ (3D centerline velocity)  RDQA (difference between jet and ambient velocities)
        DRHO = (RHOA(N2) - RHOJ(N2))/RHOA(N2)        ! *** Relative density between the jet and the ambient waters
        FTOP = G*ABS(DRHO)*RADJ(N2)                  ! *** Density momentum
        FRD2I = 0.
        EBOT = 1.
        ETOP = 0.057
        ENTS = 0.
        IF( DVEL > 0. )THEN
          FRD2I = FTOP/(ALPH2*DVEL*DVEL)
          EBOT = 1. + 5.*RDQA(1)/DVEL 
          ETOP = 0.057 + 0.554*SINPHJ*FRD2I
          ENTS = 1.414*ETOP/EBOT
        ENDIF
        
        ! ***             s     kg/m3    (-)   m/s    m       m
        DRMAJS(N2) = 2.*DTJP*PI*RHOA(N2)*ENTS*DVEL*RADJ(N2)*RLEJ(N2)     ! *** Shear at the end of the current element time step (DTJP) (kg)
        DRMAJSA = 0.5*(DRMAJS(N2) + DRMAJS(N1))                          ! *** New average shear between N2 and N1

        ! ***  Calculate forced entrainment (hardwired)
        DRMAJFO = 0.5*(DRMAJF(N2) + DRMAJF(N1))

        ENTF1 = 2.*SQRT( 1. - COSPHJ*COSPHJ*COSTHJL*COSTHJL )
        DELSIG = SIG(N2) - SIG(N1)
        ENTF2 = 0.
        ENTF3 = 0.
        IF( DELSIG > 0. )THEN
          ENTF2 = PI*COSPHJ*COSTHJL*(RADJ(N2) - RADJ(N1))/(DELSIG)
          ENTF3 = 0.5*PI*RADJ(N2)*(COSPHJ*COSTHJL - COSPHJM*COSTHJLM)/(DELSIG)
        ENDIF
        ENTF = ENTF1 + ENTF2 + ENTF3
        ENTF = MAX(ENTF,0.)
        DRMAJF(N2) = DTJP*RHOA(N2)*RADJ(N2)*RLEJ(N2)*VAH*ENTF
        IF( NJE == 2 .AND. NI == 1 )DRMAJF(N2) = 0.
        DRMAJFA = 0.5*(DRMAJF(N2) + DRMAJF(N1))
        
        ISHEAR = 0
        IFORCE = 0
        IF( ISENT(NJP) == 0 )THEN
          ! ***  TAKE MAX OF SHEAR AND FORCED
          DRMAJ = MAX(DRMAJSA,DRMAJFA)
        ELSE
          ! ***  TAKE SUM OF SHEAR AND FORCED
          DRMAJ = DRMAJSA + DRMAJFA
        ENDIF
        IF( DRMAJSA > DRMAJFA) ISHEAR = 1
        IF( DRMAJFA > DRMAJSA) IFORCE = 1

        ! *** Preparations for mass entrainment
        RMAJP(N2) = RMAJP(N1) + DRMAJ       ! *** Total volume at the end of the current time step
        RMAJI = 1./RMAJP(N2)

        ! *** Entrain mass based on the total shear (DRMAJ) and compute new concentrations
        SALJ(N2)      = RMAJI*( RMAJP(N1)*SALJ(N1)    + DRMAJ*SALA )
        TEMJ(N2)      = RMAJI*( RMAJP(N1)*TEMJ(N1)    + DRMAJ*TEMA )
        DO MD = 1,NDYE
          DYEJ(N2,MD) = RMAJI*( RMAJP(N1)*DYEJ(N1,MD) + DRMAJ*DYEA(MD) )
        ENDDO           
        SFLJ(N2)      = RMAJI*( RMAJP(N1)*SFLJ(N1)    + DRMAJ*TEMA )
        DO NT = 1,NTOX    
          TOXJ(N2,NT) = RMAJI*( RMAJP(N1)*TOXJ(N1,NT) + DRMAJ*TOXA(NT) )
        ENDDO           
        DO NS = 1,NSED   
          SEDJ(N2,NS) = RMAJI*( RMAJP(N1)*SEDJ(N1,NS) + DRMAJ*SEDA(NS) )
        ENDDO           
        DO NX = 1,NSND    
          SNDJ(N2,NX) = RMAJI*( RMAJP(N1)*SNDJ(N1,NX) + DRMAJ*SNDA(NX))
        ENDDO           
        DO NW = 1,NWQV    
          WQVJ(N2,NW) = RMAJI*( RMAJP(N1)*WQVJ(N1,NW) + DRMAJ*WQVA(NW) )
        ENDDO

        ! *** ADVANCE TOXIC PARTICULATE FRACTION
        IF( ISTRAN(5) >= 1 )THEN
          DO NT = 1,NTOX
            IF( ISTRAN(6) >= 1 )THEN
              DO NS = 1,NSED
                TMPEXP = CONPARW(NS,NT)
                IF( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
                IF( ITXPARW(NS,NT) == 1 )THEN
                  IF( SEDJ(N2,NS) > 0.) TMPVAL = SEDJ(N2,NS)**TMPEXP
                ENDIF
                TOXPFJP(N2,NS,NT) = TMPVAL*SEDJ(N2,NS)*TOXPARW(LJP,NS,NT)
              ENDDO
            ENDIF
            IF( ISTRAN(7) >= 1 )THEN
              DO NX = 1,NSND
                NS = NX + NSED
                TMPEXP = CONPARW(NS,NT)
                IF( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
                IF( ITXPARW(NS,NT) == 1 )THEN
                  IF( SNDJ(N2,NX) > 0.) TMPVAL = SNDJ(N2,NX)**TMPEXP
                ENDIF
                TOXPFJP(N2,NS,NT) = TMPVAL*SNDJ(N2,NX)*TOXPARW(LJP,NS,NT)
              ENDDO
            ENDIF
          ENDDO
          DO NT = 1,NTOX
            TOXPFTJP(N2,NT) = 0.
          ENDDO
          DO NT = 1,NTOX
            IF( ISTRAN(6) >= 1 )THEN
              DO NS = 1,NSED
                TOXPFTJP(N2,NT) = TOXPFTJP(N2,NT) + TOXPFJP(N2,NS,NT)
              ENDDO
            ENDIF
            IF( ISTRAN(7) >= 1 )THEN
              DO NX = 1,NSND
                NS = NX + NSED
                TOXPFTJP(N2,NT) = TOXPFTJP(N2,NT) + TOXPFJP(N2,NS,NT)
              ENDDO
            ENDIF
          ENDDO
          DO NT = 1,NTOX
            DO NS = 1,NSED + NSND
            TOXPFJP(N2,NS,NT) = TOXPFJP(N2,NS,NT)/(1. + TOXPFTJP(N2,NT))
            ENDDO
          ENDDO
          DO NT = 1,NTOX
            TOXPFTJP(N2,NT) = TOXPFTJP(1,NT)/(1. + TOXPFTJP(N2,NT))
          ENDDO
        ENDIF
        
        ! *** ADVANCE DENSITY
        SEDJETT = 0.
        DO NS = 1,NSED
          SEDJETT = SEDJETT + SEDJ(N2,NS)
        ENDDO
        DO NX = 1,NSND
          SEDJETT = SEDJETT + SNDJ(N2,NX)
        ENDDO
        RHOJ(N2) = FUNDEN(SALJ(N2),SEDJETT,TEMJ(N2))
        DRHO = (RHOA(N2) - RHOJ(N2))/RHOA(N2)

        ! ***  ADVANCE VELOCITY COMPONENTS
        VJG(N2) = RMAJI*( RMAJP(N1)*VJG(N1) + DRMAJ*VAG )                   ! *** 
        UJG(N2) = RMAJI*( RMAJP(N1)*UJG(N1) + DRMAJ*UAG )                   ! ***
        WJG(N2) = RMAJI*( RMAJP(N1)*WJG(N1) + DRMAJ*WAG ) + G*DRHO*DTJP     ! ***

        ! ***  NEW JET COORDINATES
        DXTMP = DTJP*UJG(N2)
        DYTMP = DTJP*VJG(N2)
        DZTMP = DTJP*WJG(N2)
        XJG(N2) = XJG(N1) + DXTMP                             ! *** New X location  delme - can be used to update cell location
        YJG(N2) = YJG(N1) + DYTMP                             ! *** New Y location
        ZJG(N2) = ZJG(N1) + DZTMP                             ! *** New Elevation
        DS = SQRT( DXTMP*DXTMP + DYTMP*DYTMP + DZTMP*DZTMP )  ! *** 3D movement
        SIG(N2) = SIG(N1) + DS                                ! *** Cumulative plume movement
        
        ! ***  Reset ambient conditions
        SEDAA = 0.0
        IF( ISAMB == 0 )THEN
          K = KQJP(NJP)
          UAG = UAGD(K)
          VAG = VAGD(K)
          WAG = WAGD(K)
          SALA = SALAD(K)
          TEMA = TEMAD(K)
          DO MD = 1,NDYE
            DYEA(MD) = DYEAD(K,MD)
          ENDDO
          SFLA = SFLAD(K)
          DO NT = 1,NTOX
            TOXA(NT) = TOXAD(K,NT)
          ENDDO
          DO NS = 1,NSED2
            SEDA(NS) = SEDAD(K,NS)
            SEDAA = SEDAA + SEDA(NS)
          ENDDO
          DO NX = 1,NSND
            SNDA(NX) = SNDAD(K,NX)
            SEDAA = SEDAA + SEDA(NX)
          ENDDO
          DO NW = 1,NWQV
            WQVA(NW) = WQVAD(K,NW)
          ENDDO
        ENDIF
        IF( ISAMB >= 1 )THEN
          ! *** Get ambient conditions from the current plume conditions
          ZVAL = ZJG(N2)
          CALL JPACON(KSZ(LJP),ZVAL,UAG,VAG,WAG)
        ENDIF
        
        VA  = SQRT(UAG*UAG + VAG*VAG + WAG*WAG)     ! *** 3D water velocity
        VAH = SQRT(UAG*UAG + VAG*VAG)               ! *** 2DH water velocity
        RHOA(N2) = FUNDEN(SALA,SEDAA,TEMA)          ! *** Ambient density

        VJ(N2)  = SQRT(UJG(N2)*UJG(N2) + VJG(N2)*VJG(N2) + WJG(N2)*WJG(N2))   ! *** 3D jet velocity
        VJH(N2) = SQRT(UJG(N2)*UJG(N2) + VJG(N2)*VJG(N2))                     ! *** 2DH jet velocity

        ! ***  Calculate global jet discharge orientation
        PHJ(N2)  = 57.2958*ATAN2(WJG(N2),VJH(N2))                             ! *** Azimuth or vertical angle (PHI) (degrees)
        THJG(N2) = 57.2958*ATAN2(VJG(N2),UJG(N2))                             ! *** Horizontal angle of the jet
        THAG     = 57.2958*ATAN2(VAG,UAG)                                     ! *** Horizontal angle of the ambient velocities
        THJL(N2) = THJG(N2) - THAG                                            ! *** Difference in ambient and jet angles
        SINPHJ   = SIN(0.0175*PHJ(N2))                                        ! *** Convert back to radians and get SIN of vertical angle
        COSPHJ   = COS(0.0175*PHJ(N2))                                        ! *** Convert back to radians and get COS of vertical angle 
        SINTHJL  = SIN(0.0175*THJL(N2))                                       ! *** Convert back to radians and get SIN of jet/ambient difference
        COSTHJL  = COS(0.0175*THJL(N2))                                       ! *** Convert back to radians and get COS of jet/ambient difference

        RDQA(N2) = COS(0.0175*PHJ(N2))*COS(0.0175*THJL(N2))*VAH               ! *** Projection of ambient on jet

        ! ***  CALCULATE NEW RADIUS AND ELEMENT LENGTH
        RLEJ(N2) = VJ(N2)*RLEJ(N1)/VJ(N1)                                     ! *** New element length
        RADJ(N2) = SQRT( RMAJP(N2)/(PI*RHOJ(N2)*RLEJ(N2)) )                   ! *** New plume radius (RMAJP - Total mass, RHOJ - Plume density, RLEJ - Element length)

        ! *** Check for convergence  
        ! ***
        ! *** DRMAJSE and DRMAJFE are relative differences in entrainment due to shear and forced, respecitvely, at the beginning and end of the NI loop
        DRMAJSE = ABS(DRMAJSA - DRMAJSO)
        DRMAJFE = ABS(DRMAJFA - DRMAJFO)
        IF( DRMAJSO > 0. ) DRMAJSE = DRMAJSE/DRMAJSO
        IF( DRMAJFO > 0. ) DRMAJFE = DRMAJFE/DRMAJFO
        ITMP = 0
        IF( DRMAJFE > DMRERM ) ITMP = 1
        IF( DRMAJSE > DMRERM ) ITMP = 1
        
        IF( LDEBUG )THEN
          IF( NJE == 2 ) WRITE(88,620)NJP, NJE, NI, ITMP, DRMAJSA, DRMAJSO, DRMAJFA, DRMAJFO
          IF( ISDJP(NJP) == 1 )THEN
            ! *** ISDJP > 0 write diagnostics
            JTMPVAL = MOD(NJE,100)
            IF( JTMPVAL == 0 )THEN
              WRITE(88,620)NJP, NJE, NI, ITMP, DRMAJSA, DRMAJSO, DRMAJFA, DRMAJFO
            ENDIF
          ENDIF
        ENDIF

        ! ***  STOP IF MAXIMUM ITERATIONS EXCEEDED
        IF( NI > NIMAX )THEN
          KFLAG = 1
          IF( LDEBUG )WRITE(88,620)NJP, NJE, NI, ITMP, DRMAJSA, DRMAJSO, DRMAJFA, DRMAJFO
          IF( LDEBUG )WRITE(88,601)TIMEDAY, NJE, NI, LJP, HP(LJP)
          
          WRITE(6,'(A,3I6,F8.3)' ) ' JET/PLUME ITERATIONS EXCEEDED N2,NI @ L,HP = ',NJE,NI,LJP,HP(LJP)
          
          IF( LOUTJET )THEN
            WRITE(10,'(A,I10,F15.5,3I6,F8.3)' ) ' JET/PLUME ITERATIONS EXCEEDED N2,NI @ L,HP = ',NITER,TIMEDAY,NJE,NI,LJP,HP(LJP)
          ENDIF

          GOTO 1500
        ENDIF
        
        ! *** Add criteria to exit entraiment calculations based on stagnant plume growth
        IF( NI > 5 .AND. RADJ(N2) <= RADJOLD(N2) .AND. RADJ(N1) <= RADJOLD(N1) )THEN
          ITMP = 0          ! *** Flag causes loop to exit below
        ENDIF

        ! *** Check for exit conditions based on current centerline or boundary
        ! ***  STOP if jet centerline penetrates surface
        IF( ISTOP == 1 )THEN
          ZJGTOP = ZJG(N2)
          IF( ZJGTOP > ZSUR )THEN
            NPRT(5) = NPRT(5) + 1
            IF( NPRT(5) < 10 )THEN
              WRITE(6,605) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            IF( LOUTJET )THEN
              WRITE(10,615) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF

            EXIT
          ENDIF
        ENDIF

        ! ***  STOP if jet centerline penetrates bottom
        IF( ISTOP == 1 )THEN
          ZJGBOT = ZJG(N2)
          IF( ZJGBOT < ZBOT )THEN
            NPRT(6) = NPRT(6) + 1
            IF( NPRT(6) < 10 )THEN
              WRITE(6,606) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            IF( LOUTJET )THEN
              WRITE(10,616) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            EXIT
          ENDIF
        ENDIF

        ! ***  STOP if jet boundary penetrates surface
        IF( ISTOP == 2 )THEN
          ZJGTOP = ZJG(N2) + RADJ(N2)*COS(0.0175*PHJ(N2))
          IF( ZJGTOP > ZSUR )THEN
            NPRT(6) = NPRT(6) + 1
            IF( NPRT(6) < 10 )THEN
              WRITE(6,602) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            IF( LOUTJET )THEN
              WRITE(10,612) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            EXIT
          ENDIF
        ENDIF

        ! ***  STOP if jet boundary penetrates bottom
        IF( ISTOP == 2 )THEN
          ZJGBOT = ZJG(N2) - RADJ(N2)*COS(0.0175*PHJ(N2))
          IF( ZJGBOT < ZBOT )THEN
            NPRT(7) = NPRT(7) + 1
            IF( NPRT(7) < 10 )THEN
              WRITE(6,603) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            IF( LOUTJET )THEN
              WRITE(10,613) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            EXIT
          ENDIF
        ENDIF

        ! ***  STOP if neutral level is reached
        ! ***  Rising plume
        IF( NJE > 4 .AND. RHOJ(N2) >= RHO_INI )THEN
          DRHOT = (RHOA(N2) - RHOJ(N2))/RHOA(N2)
          IF( DRHOT < 0. )THEN
            NPRT(1) = NPRT(1) + 1
            IF( NPRT(1) < 10 )THEN
              WRITE(6,604) "R", NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            IF( LOUTJET )THEN
              WRITE(10,614)  "R", TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            EXIT
          ENDIF
        ENDIF

        ! ***  Falling plume
        IF( NJE > 4 .AND. RHOJ(N2) < RHO_INI )THEN
          DRHOT = (RHOA(N2) - RHOJ(N2))/RHOA(N2)
          IF( DRHOT > 0. )THEN
            NPRT(2) = NPRT(2) + 1
            IF( NPRT(2) < 10 )THEN
              WRITE(6,604) "R", NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            IF( LOUTJET )THEN
              WRITE(10,614)  "R", TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            ENDIF
            
            EXIT
          ENDIF
        ENDIF

        ! ***  RETURN FOR ANOTHER ITERATION
        IF( ITMP == 1 .OR. NI == 1 )THEN
          ! XJOLD = XJG(N2)
          ! YJOLD = YJG(N2)
          ! ZJOLD = ZJG(N2)
          ! RLOLD = RLEJ(N2)
          RADJOLD = RADJ
          NI = NI + 1
          GOTO 1000
        ENDIF
        
        ! *************************************
        ! *** Entrainment iterations complete
1500    CONTINUE        
        
        ! ***  WRITE OUTPUT AND PROCEED TO NEXT JET ELEMENT
        IF( ZJG(N2) >= ZJPRT )THEN
          NPRTE = 1
          ZJPRT = ZJPRT + DZPRT
        ELSE
          NPRTE = 0
        ENDIF
        DJETI = 1./RLSCL
        XJGNE = DJETI*XJG(N2)
        YJGNE = DJETI*YJG(N2)
        ZJGNE = DJETI*ZJG(N2)
        SIGNE = DJETI*SIG(N2)
        RADJNE = DJETI*RADJ(N2)
        RLEJNE = DJETI*RLEJ(N2)
        SALJNE = SALJETI*SALJ(N2)
        TEMJNE = TEMJETI*TEMJ(N2)
        DO MD = 1,NDYE
          DYEJNE(MD) = DYEJETI(MD)*DYEJ(N2,MD)
        ENDDO
        SFLJNE = SFLJETI*SFLJ(N2)
        DO NT = 1,NTOX
          TOXJNE(NT) = TOXJETI(NT)*TOXJ(N2,NT)
        ENDDO
        DO NS = 1,NSED2
          SEDJNE(NS) = SEDJETI(NS)*SEDJ(N2,NS)
        ENDDO
        DO NX = 1,NSND
          SNDJNE(NX) = SNDJETI(NX)*SNDJ(N2,NX)
        ENDDO
        DO NW = 1,NWQV
          WQVJNE(NW) = WQVJETI(NW)*WQVJ(N2,NW)
        ENDDO

        ! *** CALCULATE ENTRAINMENT
        QJTOTO = QJTOT
        QJTOT = QJTOT*RMAJP(N2)/RMAJP(N1)
        DO K = KSZ(L),KC
          ZLOWER = Z(L,K-1)*HP(LJP) + BELV(LJP)
          ZUPPER = ZLOWER + DZC(L,K)*HP(LJP)
          IF( ZJG(N2) >= ZLOWER .AND. ZJG(N2) < ZUPPER )THEN
            QJPENT(K,NJP) = QJPENT(K,NJP) + (QJTOT - QJTOTO)
            IF( RMAJI > 0.0 )THEN
              UJPAVG(K,NJP) = UJG(N2)/RMAJI
              VJPAVG(K,NJP) = VJG(N2)/RMAJI
              WJPAVG(K,NJP) = WJG(N2)/RMAJI
            ENDIF
            EXIT
          ENDIF
        ENDDO

        ! *** CALCULATE TOTAL ENTRAINMENT
        QJPENTT(NJP) = 0.0
        DO K = KSZ(L),KC
          QJPENTT(NJP) = QJPENTT(NJP) + QJPENT(K,NJP)
        ENDDO
        
        ! *** Added additional plume migration exit criteria
        !IF( ABS(RHOJ(N2)-RHOA(N2)) < EXITRHO )THEN
        ! 
        !  IF( MOD(NJE,20) == 0 )THEN
        !    RHO_DELTA = ABS(RHOJ(N2) - RHOA(N1))
        !    IF( RHO_DELTA2 == -999. )THEN
        !      RHO_DELTA2 = RHO_DELTA
        !    ENDIF
        !  
        !    IF( RHO_INI < RHO_ELE )THEN
        !      ! *** Rising plume
        !      IF( RHO_DELTA > RHO_DELTA2 )THEN
        !        NPRT(3) = NPRT(3) + 1
        !        IF( NPRT(3) < 10 )THEN
        !          PRINT '(A,I10,F15.4,I8,I5,4F8.3,A)', 'Rising Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), '  New Exit'
        !        ENDIF
        !        
        !        IF( LOUTJET )THEN
        !          !WRITE(10, '(A,I10,F15.4,I8,I5,4F8.3,3f10.3,A)') 'Rising Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), RADJ(N2), SIG(N2), ZJG(N2), '  New Exit'
        !        ENDIF
        !    
        !        !EXIT                         ! *** Jet calculations diverging from convergence
        !      ENDIF
        !    ELSE
        !      ! *** Falling plume
        !      IF( RHOJ(N2) < RHOA(N1) )THEN
        !        NPRT(4) = NPRT(4) + 1
        !        IF( NPRT(4) < 10 )THEN
        !          PRINT '(A,I10,F15.4,I8,I5,4F8.3,A)', 'Falling Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), '  New Exit'
        !        ENDIF
        !        
        !        IF( LOUTJET )THEN
        !          WRITE(10, '(A,I10,F15.4,I8,I5,4F8.3,A)') 'Rising Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), ZJG(N2), SIG(N2), '  New Exit'
        !        ENDIF
        !        
        !        EXIT                         ! *** Jet calculations diverging from convergence
        !      ENDIF              
        !    ENDIF
        !    RHO_DELTA2 = RHO_DELTA
        !    RHOA(N1)   = RHOA(N2)               ! *** Density of ambient water at face
        !    RHOJ(N1)   = RHOJ(N2)               ! *** Density of jet
        !  ENDIF        
        !else
        !  nx = 0 ! DELME
        !ENDIF

        ! *** Write to file when diameter changes by 0.1 m
        IDUMP = 0
        IF( IOUTJP(NJP) > 1 )THEN
          IF(  NJE == 1 .OR. RADJ(N2)-RADLAST > 0.05 )THEN
            WRITE(10,631) NITER, TIMEDAY, NJE, NI, TOTTIME, SIG(N2),                                                                          &
                          XJG(N2), YJG(N2), ZJG(N2), 2.*RADJ(N2), QVJET, SUM(QJPENT(:,NJP)), QJTOT, RHOA(N2), -999., RHOJ(N2), TEMJ(N2),      &
                          VJ(N2), VJH(N2), WJG(N2), WAG, RHOJ(N2)-RHOA(N2), RHO_DELTA, RHO_DELTA2, TEMA, VA, VAH, -999., RMAJP(N2), DRMAJS(N2), DRMAJF(N2)
            RADLAST = RADJ(N2)
            IDUMP = 1
          ENDIF
        ENDIF
        
        ! *** Advance the variables
        XJG(N1)    = XJG(N2)                  ! *** X coordinate of plume face
        YJG(N1)    = YJG(N2)                  ! *** Y coordinate of plume face 
        ZJG(N1)    = ZJG(N2)                  ! *** Z coordinate of plume face 
        RADJ(N1)   = RADJ(N2)                 ! *** Radius of jet
        RLEJ(N1)   = RLEJ(N2)                 ! *** Element length
        SIG(N1)    = SIG(N2)                  ! *** Cumulative distance traveled along plume centerline
        UJG(N1)    = UJG(N2)                  ! *** Jet velocity - u component
        VJG(N1)    = VJG(N2)                  ! *** Jet velocity - v component
        WJG(N1)    = WJG(N2)                  ! *** Jet velocity - w component
        VJ(N1)     = VJ(N2)                   ! *** 
        VJH(N1)    = VJH(N2)                  ! *** 
        RMAJP(N1)  = RMAJP(N2)                ! *** Jet mass
        DRMAJS(N1) = DRMAJS(N2)               ! *** 
        DRMAJF(N1) = DRMAJF(N2)               ! *** 
        PHJ(N1)    = PHJ(N2)                  ! *** Vertical angle of jet
        THJG(N1)   = THJG(N2)                 ! *** Angle of jet
        THJL(N1)   = THJL(N2)                 ! *** Difference in angle of jet - ambient
        RDQA(N1)   = RDQA(N2)                 ! *** 
        SALJ(N1)   = SALJ(N2)                 ! *** 
        TEMJ(N1)   = TEMJ(N2)                 ! *** 
        DO MD = 1,NDYE
          DYEJ(N1,MD) = DYEJ(N2,MD)
        ENDDO
        SFLJ(N1) = SFLJ(N2)
        DO NT = 1,NTOX
          TOXJ(N1,NT) = TOXJ(N2,NT)
        ENDDO
        DO NS = 1,NSED2
          SEDJ(N1,NS) = SEDJ(N2,NS)
        ENDDO
        DO NX = 1,NSND
          SNDJ(N1,NX) = SNDJ(N2,NX)
        ENDDO
        DO NW = 1,NWQV
          WQVJ(N1,NW) = WQVJ(N2,NW)
        ENDDO
       
      ENDDO  ! *** End NJEL element integration

      ! *** Output results for last element
      IF( IDUMP == 1 )THEN
            WRITE(10,631) NITER, TIMEDAY, NJE, NI, TOTTIME, SIG(N2),                                                                          &
                          XJG(N2), YJG(N2), ZJG(N2), 2.*RADJ(N2), QVJET, SUM(QJPENT(:,NJP)), QJTOT, RHOA(N2), -999., RHOJ(N2), TEMJ(N2),      &
                          VJ(N2), VJH(N2), WJG(N2), WAG, RHOJ(N2)-RHOA(N2), RHO_DELTA, RHO_DELTA2, TEMA, VA, VAH, -999., RMAJP(N2), DRMAJS(N2), DRMAJF(N2)
      ENDIF
631   FORMAT(I10,F15.4,I8,I5,2F9.3,2F12.4,F10.3,F8.2,3E12.4,8f10.3,10E12.4)  
    
      ! *** COMPUTE THE LAYER SPECIFIC JET/PLUME ENTRAINMENT
      QJTOTO = QJTOT
      QJTOT = QJTOT*RMAJP(N2)/RMAJP(N1)
      DO K = KSZ(L),KC
        ZLOWER = Z(L,K - 1)*HP(LJP) + BELV(LJP)
        ZUPPER = ZLOWER + DZC(L,K)*HP(LJP)
        IF( ZJG(N2) >= ZLOWER .AND. ZJG(N2) < ZUPPER )THEN
          ! *** Found Layer
          QJPENT(K,NJP) = QJPENT(K,NJP) + (QJTOT - QJTOTO)
          EXIT
        ENDIF
      ENDDO
      
      ! *** GET THE TOTAL ENTRAINMENT VOLUMES FOR CURRENT NJP
      QJPENTT(NJP) = 0.0
      DO K = KSZ(L),KC
        QJPENTT(NJP) = QJPENTT(NJP) + QJPENT(K,NJP)
      ENDDO

      IF( LDEBUG )THEN
        ! *** OPEN LOF FILE
        IF( IFILE == -1 )THEN
          IFILE = 8
          OPEN(8,FILE = OUTDIR//'EFDCLOG.OUT',POSITION = 'APPEND')
        ENDIF
        WRITE(8,898)NJP,TIME,(QJPENT(K,NJP),K = 1,KC),QJPENTT(NJP)
        WRITE(88,898)NJP,TIME,(QJPENT(K,NJP),K = 1,KC),QJPENTT(NJP)
      ENDIF
      DJETI = 1./RLSCL
      XJGNE = DJETI*XJG(N2)
      YJGNE = DJETI*YJG(N2)
      ZJGNE = DJETI*ZJG(N2)
      SIGNE = DJETI*SIG(N2)
      RADJNE = DJETI*RADJ(N2)
      RLEJNE = DJETI*RLEJ(N2)
      SALJNE = SALJETI*SALJ(N2)
      TEMJNE = TEMJETI*TEMJ(N2)
      DO MD = 1,NDYE
        DYEJNE(MD) = DYEJETI(MD)*DYEJ(N2,MD)
      ENDDO
      SFLJNE = SFLJETI*SFLJ(N2)
      DO NT = 1,NTOX
        TOXJNE(NT) = TOXJETI(NT)*TOXJ(N2,NT)
      ENDDO
      DO NS = 1,NSED2
        SEDJNE(NS) = SEDJETI(NS)*SEDJ(N2,NS)
      ENDDO
      DO NX = 1,NSND
        SNDJNE(NX) = SNDJETI(NX)*SNDJ(N2,NX)
      ENDDO
      DO NW = 1,NWQV
        WQVJNE(NW) = WQVJETI(NW)*WQVJ(N2,NW)
      ENDDO
      DRMAJ = RMAJP(N2) - RMAJP(N2-1)
      DRHO = (RHOA(N2) - RHOJ(N2))/RHOA(N2)

      ! *** DIAGNOSTIC OUPUT
      IF( LOUTJET .AND. (IOUTJP(NJP) == 1 .OR. IOUTJP(NJP) == 3 ) )THEN

      ENDIF

      ! *** DIAGNOSTIC OUPUT

      ! *** RELOCATE VOLUME SOURCE
      ZTMP = ( ZJGNE-BELV(LJP) )/HP(LJP)
      ZTMP = RKC*ZTMP
      KTMP = INT(ZTMP) + KL 
      KTMP = MAX(KL,KTMP)
      KTMP = MIN(KC,KTMP)
      IF( KFLAG == 0 ) KEFFJP(NJP) = KTMP
      GOTO 9001

      ! *** JET/PLUME COMPUTATIONS BYPASSED
      9000 CONTINUE
      KEFFJP(NJP) = KQJP(NJP)
      
9001  CONTINUE
      
      IF( LDEBUG )THEN
        ! *** OPEN LOF FILE
        IF( IFILE == -1 )THEN
          IFILE = 8
          OPEN(8,FILE = OUTDIR//'EFDCLOG.OUT',POSITION = 'APPEND')
        ENDIF
        WRITE(8 ,899)NJP,TIME,(QJPENT(K,NJP),K = 1,KC)
        WRITE(8 ,135)NJP,TIME,KFLAG,KEFFJP(NJP),KQJP(NJP),QVJET,QJTOT
        WRITE(88,899)NJP,TIME,(QJPENT(K,NJP),K = 1,KC)
        WRITE(88,135)NJP,TIME,KFLAG,KEFFJP(NJP),KQJP(NJP),QVJET,QJTOT
      ENDIF

      ! *** CALCULATION MOMENT INTERFACE QUANTITIES
      RDUM = 0.
      IF( LDEBUG )THEN
        WRITE(11,1110)NJP,N
        KZERO = 0
        QUJ0 = UJ0*QVJET0
        QVJ0 = VJ0*QVJET0
        QWJ0 = WJ0*QVJET0
        WRITE(11,1111)KZERO,QUJ0,QVJ0,QWJ0,QVJET0,RDUM,RDUM,RDUM,RDUM
        QENTTMP = QVJET0
        DO K = 1,KEFFJP(NJP)
          QENTTMP = QENTTMP + QJPENT(K,NJP)
          QUAG = UAG*QJPENT(K,NJP)
          QVAG = VAG*QJPENT(K,NJP)
          QWAG = WAG*QJPENT(K,NJP)
          WRITE(11,1111)K,UJPAVG(K,NJP),VJPAVG(K,NJP),WJPAVG(K,NJP), &
                        QENTTMP,QUAG,QVAG,QWAG,QJPENT(K,NJP)
        ENDDO
      ENDIF

      ! *** END LOOP OVER ALL JET/PLUME LOCATIONS
    ENDIF   ! *** End of active jet/plume BC flag
  ENDDO     ! *** End of NQJPIJ loop

  ! *** CLOSE LOG FILE
  IF( IFILE == 8 ) CLOSE(8)
  
  IF( LOUTJET )THEN
    CLOSE(10)
  ENDIF
  IF( LDEBUG )THEN
    CLOSE(88)
  ENDIF
    
1110 FORMAT(/,'NJP,N = ',I5,I10,/)
1111 FORMAT(I5,50E14.5)
 899 FORMAT(' JPENT ',I5,F12.6,50E12.4)
 898 FORMAT(' FINAL JPENT ',I5,F12.6,50E12.4)
 100 FORMAT(120X)

 134 FORMAT(' BEGIN JET/PLUME NJP,TIME = ',I6,F12.5)
 135 FORMAT(' END JET/PLUME NJP,TIME,KFLAG,KEFFJP,KQJP,QVJET,QVJTOT',' = ',I6,F13.5,3I4,2E12.4)
       
 601 FORMAT(' MAXIMUM ITERATIONS EXCEEDED TIMEDAY N2,NI @ L,HP = ',F10.3,3I6,F8.3,' !!')
 602 FORMAT(' JET/PLUME BOUNDARY PENETRATES SURFACE @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
 603 FORMAT(' JET/PLUME BOUNDARY PENETRATES BOTTOM @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
 604 FORMAT(' JET/PLUME AT NEUTRAL LEVEL @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',A,I5,I8,I5,3F10.3)
 605 FORMAT(' JET/PLUME CENTERLINE PENETRATES SURFACE @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
 606 FORMAT(' JET/PLUME CENTERLINE PENETRATES BOTTOM @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
       
 ! *** Write to log file formats
 612 FORMAT(' JET/PLUME BOUNDARY PENETRATES SURFACE @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
 613 FORMAT(' JET/PLUME BOUNDARY PENETRATES BOTTOM @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
 614 FORMAT(' JET/PLUME AT NEUTRAL LEVEL @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',A,F15.5,I5,I8,I5,3F10.3)
 615 FORMAT(' JET/PLUME CENTERLINE PENETRATES SURFACE @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
 616 FORMAT(' JET/PLUME CENTERLINE PENETRATES BOTTOM @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
       
 888 FORMAT(A80,/)
 620 FORMAT('NJ,N2,NI,IT,DS,DSO,DF,DFO = ',4I6,6E13.4)
  
  RETURN
  
END

FUNCTION FUNDEN(SALIN, SEDIN, TEMIN)

  ! **  FUNDEN CALCULATED DENSITY AS A FUNCTION OF SAL,TEM,AND SED
  ! CHANGE RECORD
  !  2022-04 Szu-Ting Lee & Paul M. Craig 
  !             Updated funciton to Millero & Huang (2008) to allow for broader salinity and temperature ranges

  IMPLICIT NONE

  REAL, INTENT(IN) :: SALIN, SEDIN, TEMIN
  REAL             :: FUNDEN, TTMP, SSG, SDEN, SSTMP, RHTMP, RHO, RHOCA, RHOCB, RHOCC 

  SSG = 2.5
  SDEN = 1./2500000.

  ! **  DENSITY AT GIVEN VALUES OF SAL AND TEM
  SSTMP = MIN(MAX(SALIN,0.),70.)
  TTMP  = MIN(MAX(TEMIN,0.),90.)
  RHTMP = 999.842594
  
  RHTMP = RHTMP + + 6.793952E-2*TTMP - 9.095290E-3*TTMP*TTMP +1.001685E-4*TTMP*TTMP*TTMP - 1.120083E-6*TTMP*TTMP*TTMP*TTMP + 6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
  
  RHOCA = 0.8174451 - 3.638577E-3*TTMP + 6.480811E-5*TTMP*TTMP - 7.312404E-7*TTMP*TTMP*TTMP +5.330431E-9*TTMP*TTMP*TTMP*TTMP - 1.657628E-11*TTMP*TTMP*TTMP*TTMP*TTMP
  RHOCB = -5.481436E-3 + 3.486075E-5*TTMP - 3.049727E-7*TTMP*TTMP
  RHOCC = 5.346196E-4
  
  RHO = RHTMP + SSTMP*RHOCA + SQRT(SSTMP)*SSTMP*RHOCB + SSTMP*SSTMP*RHOCC
  
  ! **  CORRECTION FOR SEDIMENT
  RHO = RHO*( (1. - SDEN*SEDIN) + (SSG - 1.)*SDEN*SEDIN )

  ! **  RETURN DENSITY
  FUNDEN = RHO

END FUNCTION
  
  
SUBROUTINE JPACON(KBOT, ZVAL, UAG, VAG, WAG)

  ! *** JET/PLUME SUB - MODEL SUBROUTINE
  ! *** RETURNS THE AMBIENT CONDITIONS FOR VARIABLE WQ CONDITIONS BASED 

  ! CHANGE RECORD
  
  USE GLOBAL
  IMPLICIT NONE   
  
  INTEGER,INTENT(IN) :: KBOT                                                                                                     
  INTEGER :: K, NT, NS, NX, NZP, MD, NW                                                                        
  REAL    :: WTNZP, WTNZ, ZVAL                                                                                                
  REAL    :: UAG, VAG, WAG, DZK

  K = KBOT

  IF( ZVAL < ZAD(K) )THEN
    ! *** Jet/Plume is below the bottom layer
    UAG = UAGD(K)
    VAG = VAGD(K)
    WAG = WAGD(K)
    SALA = SALAD(K)
    TEMA = TEMAD(K)
    DO MD = 1,NDYE
      DYEA(MD) = DYEAD(K,MD)
    ENDDO
    SFLA = SFLAD(K)
    DO NT = 1,NTOX
      TOXA(NT) = TOXAD(K,NT)
    ENDDO
    DO NS = 1,NSED2
      SEDA(NS) = SEDAD(K,NS)
    ENDDO
    DO NX = 1,NSND
      SNDA(NX) = SNDAD(K,NX)
    ENDDO
    DO NW = 1,NWQV
      WQVA(NW) = WQVAD(K,NW)
    ENDDO
    RETURN
  ENDIF

  ! *** GET CONCENTRATIONS BASED ON PLUME ELEVATION (ZVAL)
  IF( ZVAL >= ZAD(KC) )THEN
    ! Jet/plume elevation is above the top layer midpoint
    UAG = UAGD(KC)
    VAG = VAGD(KC)
    WAG = WAGD(KC)
    SALA = SALAD(KC)
    TEMA = TEMAD(KC)
    DO MD = 1,NDYE
      DYEA(MD) = DYEAD(KC,MD)
    ENDDO
    SFLA = SFLAD(KC)
    DO NT = 1,NTOX
      TOXA(NT) = TOXAD(KC,NT)
    ENDDO
    DO NS = 1,NSED2
      SEDA(NS) = SEDAD(KC,NS)
    ENDDO
    DO NX = 1,NSND
      SNDA(NX) = SNDAD(KC,NX)
    ENDDO
    DO NW = 1,NWQV
      WQVA(NW) = WQVAD(KC,NW)
    ENDDO
    RETURN
  ENDIF

  ! *** Determine which layer to extract the ambient conditions from
  
  ! *** Loop over layers
   1000 CONTINUE

  ! *** NZP = LAYER ABOVE
  NZP = K + 1
  IF( ZVAL >= ZAD(K) .AND. ZVAL < ZAD(NZP) )THEN
    ! *** INTERPOLATE AMBIENT CONDITIONS FROM LAYER INTERVAL
    DZK =1./(ZAD(NZP) - ZAD(K))
    WTNZ = DZK*(ZAD(NZP) - ZVAL)
    WTNZP = DZK*(ZVAL - ZAD(K))
    UAG = WTNZ*UAGD(K) + WTNZP*UAGD(NZP)
    VAG = WTNZ*VAGD(K) + WTNZP*VAGD(NZP)
    WAG = WTNZ*WAGD(K) + WTNZP*WAGD(NZP)
    SALA = WTNZ*SALAD(K) + WTNZP*SALAD(NZP)
    TEMA = WTNZ*TEMAD(K) + WTNZP*TEMAD(NZP)
    DO MD = 1,NDYE
      DYEA(MD) = WTNZ*DYEAD(K,MD) + WTNZP*DYEAD(NZP,MD)
    ENDDO
    SFLA = WTNZ*SFLAD(K) + WTNZP*SFLAD(NZP)
    DO NT = 1,NTOX
      TOXA(NT) = WTNZ*TOXAD(K,NT) + WTNZP*TOXAD(NZP,NT)
    ENDDO
    DO NS = 1,NSED2
      SEDA(NS) = WTNZ*SEDAD(K,NS) + WTNZP*SEDAD(NZP,NS)
    ENDDO
    DO NX = 1,NSND
      SNDA(NX) = WTNZ*SNDAD(K,NX) + WTNZP*SNDAD(NZP,NX)
    ENDDO
    DO NW = 1,NWQV
      WQVA(NW) = WTNZ*WQVAD(K,NW) + WTNZP*WQVAD(NZP,NW)
    ENDDO
    RETURN
  ELSE
    K = K + 1
    GOTO 1000
  ENDIF
  RETURN
END

END MODULE
