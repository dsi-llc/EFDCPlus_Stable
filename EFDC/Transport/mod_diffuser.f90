! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE DIFFUSER_MODULE

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  ! - ---------------------------------------------------------------------!
  !    2018 - 05       PAUL M. CRAIG     RESTRUCTRED EFDC+ DIFFUSER CODE FOR BETTER MANAGEMENT

  use GLOBAL
  use Variables_WQ
  
  use Variables_MPI
  
  implicit none

  real :: SALA
  real :: SFLA
  real :: TEMA

  real,save,allocatable,dimension(:)   :: UAGD
  real,save,allocatable,dimension(:)   :: VAGD
  real,save,allocatable,dimension(:)   :: WAGD

  real,save,allocatable,dimension(:)   :: DYEA
  real,save,allocatable,dimension(:,:) :: DYEAD
  real,save,allocatable,dimension(:)   :: SALAD
  real,save,allocatable,dimension(:)   :: SEDA
  real,save,allocatable,dimension(:,:) :: SEDAD
  real,save,allocatable,dimension(:)   :: SFLAD
  real,save,allocatable,dimension(:)   :: SNDA
  real,save,allocatable,dimension(:,:) :: SNDAD
  real,save,allocatable,dimension(:)   :: TAD
  real,save,allocatable,dimension(:)   :: TEMAD
  real,save,allocatable,dimension(:)   :: TOXA
  real,save,allocatable,dimension(:,:) :: TOXAD
  real,save,allocatable,dimension(:)   :: WQVA
  real,save,allocatable,dimension(:,:) :: WQVAD
  real,save,allocatable,dimension(:)   :: ZAD
  logical :: LDEBUG
  
  contains
  
SUBROUTINE JPEFDC

  ! *** PROGRAM JPEFDC IS STAND ALONE VERSION OF EFDC JET - PLUME MODEL
  ! *** BASED ON LEE AND CHEUNG'S LAGRANGIAN BUOYANT JET MODEL
  ! *** AND EXTENDED FOR THREE-DIMENSIONAL AMBIENT CURRENTS
  ! ***   REF: LEE, J.H.W., AND V. CHEUNG, J. ENVIRON. ENGR., 116, 1085 - 
  ! ***   1106, 1990.

  use Allocate_Initialize
  
  implicit none

  integer,parameter :: NJELM = 2, NATDM = 1
  integer           :: NJP, K, LJP, KJP, KFLAG, ISAMB, ISJPD, ISTOP, ISOUT, NPRTE, NIMAX, ITMP, JTMPVAL, KL
  integer           :: NT, NV, NX, NS, L, LU, KU, LE, LW, LN, LS, N2, NJE, N1, NI, ISHEAR, IFORCE
  integer           :: KTMP, KZERO, IFILE, MD, NW
  integer           :: IDUMP
  integer,save      :: NPRT(10)
  
  real :: TIME, RADDEG, ZTMP, ALPHA, ALPH2, TUJP, XERRM, YERRM, ZERRM, RRERM, RKC
  real :: RLERM, DMRERM, XJET, YJET, RJET, AJET, QVJET, QSERTAVG, QVJET0, VJET, WTC
  real :: WTV, TMPVAL, SALJET, TEMJET, SFLJET, ZBOT, ZSUR, DZPRT, ZJPRT
  real :: DTJP, UJ0, VJ0, WJ0, SEDJETT, UTMP, VTMP, TMPEXP, VA, VAH, SEDATT
  real :: FRD2I, EBOT, ETOP, ENTS, QJTOT, QJTOTO, RLSCL, SALJETI, TEMJETI, SFLJETI
  real :: THAG, SINPHJ, DVEL, DRHO, FTOP, DJETI, XJGNE, YJGNE, ZJGNE, RADJNE, SIGNE
  real :: RLEJNE, SALJNE, TEMJNE, SFLJNE, DRMAJ, COSPHJ, COSPHJM, SINTHJL, COSTHJL
  real :: COSTHJLM, DRMAJSO, DRMAJSA, DRMAJFO, ENTF1
  real :: PTINI, RHO_INI, RHO_ELE , DELTEM, DELRHO, EXITTEM, EXITRHO, RHO_DELTA, RHO_DELTA2
  ! REAL :: XJOLD, YJOLD, ZJOLD   ! DELME
  real :: DELSIG, ENTF2, ENTF3, ENTF, DRMAJFA, RMAJI, DXTMP, DYTMP, DZTMP, DS, SEDAA, TOTTIME
  real :: DRMAJSE, DRMAJFE, ZJGTOP, ZJGBOT, DRHOT, ZLOWER, ZUPPER, RDUM, QUJ0, QVJ0, QWJ0
  real :: QENTTMP, QUAG, QVAG, QWAG, ZVAL, RADLAST, VJET0
  real :: UAG, VAG, WAG
  real,save :: CURDAY
  
  !REAL,external :: FUNDEN
  
  character :: FNNUM*3, PROCNUM*3, FNJPLOG*40

  logical*4,save :: LOUTJET
  
  real,save,allocatable,dimension(:)   :: DRHONS
  real,save,allocatable,dimension(:)   :: DRMAJF
  real,save,allocatable,dimension(:)   :: DRMAJNS
  real,save,allocatable,dimension(:)   :: DRMAJS
  real,save,allocatable,dimension(:,:) :: DYEJ
  real,save,allocatable,dimension(:)   :: DYEJET
  real,save,allocatable,dimension(:)   :: DYEJETI
  real,save,allocatable,dimension(:)   :: DYEJNE
  real,save,allocatable,dimension(:,:) :: DYEJNS
  real,save,allocatable,dimension(:)   :: PHJ
  real,save,allocatable,dimension(:)   :: VJ
  real,save,allocatable,dimension(:)   :: VJH
  real,save,allocatable,dimension(:)   :: QJNS
  real,save,allocatable,dimension(:)   :: RADJ
  real,save,allocatable,dimension(:)   :: RADJNS
  real,save,allocatable,dimension(:)   :: RADJOLD
  real,save,allocatable,dimension(:)   :: RDQA
  real,save,allocatable,dimension(:)   :: RDQANS
  real,save,allocatable,dimension(:)   :: RHOA
  real,save,allocatable,dimension(:)   :: RHOJ
  real,save,allocatable,dimension(:)   :: RLEJ
  real,save,allocatable,dimension(:)   :: RLEJNS
  real,save,allocatable,dimension(:)   :: RMAJP
  real,save,allocatable,dimension(:)   :: RMAJPNS
  real,save,allocatable,dimension(:)   :: SALJ
  real,save,allocatable,dimension(:)   :: SALJNS
  real,save,allocatable,dimension(:)   :: SEDJET
  real,save,allocatable,dimension(:)   :: SEDJETI
  real,save,allocatable,dimension(:)   :: SEDJNE
  real,save,allocatable,dimension(:)   :: SFLJ
  real,save,allocatable,dimension(:)   :: SIG
  real,save,allocatable,dimension(:)   :: SIGNS
  real,save,allocatable,dimension(:)   :: SNDJET
  real,save,allocatable,dimension(:)   :: SNDJETI
  real,save,allocatable,dimension(:)   :: SNDJNE
  real,save,allocatable,dimension(:)   :: TEMJ
  real,save,allocatable,dimension(:)   :: TEMJNS
  real,save,allocatable,dimension(:)   :: THJG
  real,save,allocatable,dimension(:)   :: THJL
  real,save,allocatable,dimension(:)   :: TIMJP
  real,save,allocatable,dimension(:)   :: TOXJET
  real,save,allocatable,dimension(:)   :: TOXJETI
  real,save,allocatable,dimension(:)   :: TOXJNE
  real,save,allocatable,dimension(:)   :: UJG
  real,save,allocatable,dimension(:)   :: UJGNS
  real,save,allocatable,dimension(:)   :: VJG
  real,save,allocatable,dimension(:)   :: VJGNS
  real,save,allocatable,dimension(:)   :: WJG
  real,save,allocatable,dimension(:)   :: WJGNS
  real,save,allocatable,dimension(:)   :: XJG
  real,save,allocatable,dimension(:)   :: XJGNS
  real,save,allocatable,dimension(:)   :: YJG
  real,save,allocatable,dimension(:)   :: YJGNS
  real,save,allocatable,dimension(:)   :: ZJG
  real,save,allocatable,dimension(:)   :: ZJGNS
  real,save,allocatable,dimension(:,:) :: SEDJ
  real,save,allocatable,dimension(:,:) :: SEDJNS
  real,save,allocatable,dimension(:,:) :: SNDJ
  real,save,allocatable,dimension(:,:) :: SNDJNS
  real,save,allocatable,dimension(:,:) :: TOXJ
  real,save,allocatable,dimension(:,:) :: TOXJNS
  real,save,allocatable,dimension(:,:,:) :: TOXPFJP
  real,save,allocatable,dimension(:,:) :: TOXPFTJP
  real,save,allocatable,dimension(:,:) :: TOXPFTNS
  real,save,allocatable,dimension(:,:) :: UJPAVG
  real,save,allocatable,dimension(:,:) :: VJPAVG
  real,save,allocatable,dimension(:,:) :: WJPAVG
  real,save,allocatable,dimension(:,:) :: WQVJ
  real,save,allocatable,dimension(:)   :: WQVJET
  real,save,allocatable,dimension(:)   :: WQVJETI
  real,save,allocatable,dimension(:)   :: WQVJNE
  real,save,allocatable,dimension(:,:) :: WQVJNS
  real,save,allocatable,dimension(:)   :: RPORTS
  
  LDEBUG = DEBUG
  !LDEBUG = .TRUE.   ! delme
  
  ! **************************************************************************************************
  ! **************************************************************************************************
  ! *** Allocation Block
  if( .not. allocated(DRHONS) )then
    write(*,'(A,I8)')'JET-PLUME COMPUTATIONS STARTED.  NQJPIJ = ',NQJPIJ
    
    ! *** MODULE LEVEL VARIABLES FIRST
    call AllocateDSI( DYEA,     NDYM,   0.0)
    call AllocateDSI( DYEAD,    KCM,    NDYM,   0.0)
    call AllocateDSI( SALAD,    KCM,    0.0)    
    call AllocateDSI( SEDA,     NSED2,  0.0)      
    call AllocateDSI( SEDAD,    KCM,    NSEDS2,  0.0)
    call AllocateDSI( SFLAD,    KCM,    0.0)    
    call AllocateDSI( SNDA,     NSND,   0.0)    
    call AllocateDSI( SNDAD,    KCM,    NSND,   0.0)
    call AllocateDSI( TEMAD,    KCM,    0.0)    
    call AllocateDSI( TAD,      KCM,    0.0)    
    call AllocateDSI( TOXA,     NTXM,   0.0)    
    call AllocateDSI( TOXAD,    KCM,    NTXM,   0.0)
    call AllocateDSI( UAGD,     KCM,    0.0)
    call AllocateDSI( VAGD,     KCM,    0.0)
    call AllocateDSI( WAGD,     KCM,    0.0)
    call AllocateDSI( ZAD,      KCM,    0.0)
                      
    call AllocateDSI( DRHONS,   NJPSM,  0.0)      
    call AllocateDSI( DRMAJF,   NJELM,  0.0)      
    call AllocateDSI( DRMAJNS,  NJPSM,  0.0)      
    call AllocateDSI( DRMAJS,   NJELM,  0.0)      
    call AllocateDSI( DYEJ,     NJELM,  NDYM,   0.0)
    call AllocateDSI( DYEJET,   NDYM,   0.0)      
    call AllocateDSI( DYEJETI,  NDYM,   0.0)      
    call AllocateDSI( DYEJNE,   NDYM,   0.0)      
    call AllocateDSI( DYEJNS,   NDYM,   NJPSM,  0.0)
    call AllocateDSI( PHJ,      NJELM,  0.0)      
    call AllocateDSI( VJ,       NJELM,  0.0)      
    call AllocateDSI( VJH,      NJELM,  0.0)      
    call AllocateDSI( QJNS,     NJPSM,  0.0)      
    call AllocateDSI( RADJ,     NJELM,  0.0)   
    call AllocateDSI( RADJNS,   NJPSM,  0.0)      
    call AllocateDSI( RADJOLD,  NJELM,  0.0)   
    call AllocateDSI( RDQA,     NJELM,  0.0)      
    call AllocateDSI( RDQANS,   NJPSM,  0.0)      
    call AllocateDSI( RHOA,     NJELM,  0.0)      
    call AllocateDSI( RHOJ,     NJELM,  0.0)      
    call AllocateDSI( RLEJ,     NJELM,  0.0)      
    call AllocateDSI( RLEJNS,   NJPSM,  0.0)      
    call AllocateDSI( RMAJP,   -NJELM,  0.0)      
    call AllocateDSI( RMAJPNS,  NJPSM,  0.0)      
    call AllocateDSI( SALJ,     NJELM,  0.0)      
    call AllocateDSI( SALJNS,   NJPSM,  0.0)      
    call AllocateDSI( SEDJ,     NJELM,  NSEDS,   0.0)
    call AllocateDSI( SEDJET,   NSEDS2, 0.0)      
    call AllocateDSI( SEDJETI,  NSEDS2, 0.0)      
    call AllocateDSI( SEDJNE,   NSEDS2, 0.0)      
    call AllocateDSI( SEDJNS,   NSEDS2, NJPSM,  0.0)
    call AllocateDSI( SFLJ,     NJELM,  0.0)      
    call AllocateDSI( SIG,      NJELM,  0.0)      
    call AllocateDSI( SIGNS,    NJPSM,  0.0)      
    call AllocateDSI( SNDJ,     NJELM,  NSNM,   0.0)
    call AllocateDSI( SNDJET,   NSNM,   0.0)      
    call AllocateDSI( SNDJETI,  NSNM,   0.0)      
    call AllocateDSI( SNDJNE,   NSNM,   0.0)      
    call AllocateDSI( SNDJNS,   NSNM,   NJPSM,  0.0)
    call AllocateDSI( TEMJ,     NJELM,  0.0)      
    call AllocateDSI( TEMJNS,   NJPSM,  0.0)      
    call AllocateDSI( THJG,     NJELM,  0.0)      
    call AllocateDSI( THJL,     NJELM,  0.0)      
    call AllocateDSI( TIMJP,    NJPSM,  0.0)      
    call AllocateDSI( TOXJ,     NJELM,  NTXM,   0.0)
    call AllocateDSI( TOXJET,   NTXM,   0.0)      
    call AllocateDSI( TOXJETI,  NTXM,   0.0)      
    call AllocateDSI( TOXJNE,   NTXM,   0.0)      
    call AllocateDSI( TOXJNS,   NTXM,   NJPSM,  0.0)
    call AllocateDSI( TOXPFJP,  NJELM,  NSTM,   NTXM, 0.0)
    call AllocateDSI( TOXPFTJP, NJELM,  NTXM,   0.0)
    call AllocateDSI( TOXPFTNS, NTXM,   NJPSM,  0.0)
    call AllocateDSI( UJG,      NJELM,  0.0)
    call AllocateDSI( UJGNS,    NJPSM,  0.0)
    call AllocateDSI( UJPAVG,   KCM,    NJPSM,  0.0)
    call AllocateDSI( VJG,      NJELM,  0.0)    
    call AllocateDSI( VJGNS,    NJPSM,  0.0)    
    call AllocateDSI( VJPAVG,   KCM,    NJPSM,  0.0)
    call AllocateDSI( WJG,      NJELM,  0.0)    
    call AllocateDSI( WJGNS,    NJPSM,  0.0)    
    call AllocateDSI( WJPAVG,   KCM,    NJPSM,  0.0)
    call AllocateDSI( XJG,      NJELM,  0.0)
    call AllocateDSI( XJGNS,    NJPSM,  0.0)
    call AllocateDSI( YJG,      NJELM,  0.0)
    call AllocateDSI( YJGNS,    NJPSM,  0.0)
    call AllocateDSI( ZJG,      NJELM,  0.0)
    call AllocateDSI( ZJGNS,    NJPSM,  0.0)
    call AllocateDSI( RPORTS,   NQJPIJ, 0.0)

    if( ISTRAN(8) > 0 )then
      call AllocateDSI( WQVJ,     NJELM, -NWQV,   0.0)
      call AllocateDSI( WQVA,    -NWQV,   0.0)
      call AllocateDSI( WQVAD,    KCM,   -NWQV,   0.0)
      call AllocateDSI( WQVJET,  -NWQV,   0.0)      
      call AllocateDSI( WQVJETI, -NWQV,   0.0)      
      call AllocateDSI( WQVJNE,  -NWQV,   0.0)      
      call AllocateDSI( WQVJNS,  -NWQV,   NJPSM,  0.0)
    endif
    
    ! *** INITIALIZE LOCAL ARRAYS
    NPRT = 0
    CURDAY = TIMEDAY + 1./24.
    RMAJP = 0.0

    LOUTJET = .FALSE.
    do NJP = 1,NQJPIJ
      if( JET_PLM(NJP).IOUTJP > 0 ) LOUTJET = .TRUE.
      RPORTS(NJP) = FLOAT(JET_PLM(NJP).NPORTJP)  
    enddo
    
    if( LOUTJET )then
      open(88,FILE = OUTDIR//'JPBUG.DIA', STATUS = 'UNKNOWN')
      close(88,STATUS = 'DELETE')
    
      do NJP = 1,NQJPIJ
        if( JET_PLM(NJP).IOUTJP > 0 )then
          write(FNNUM,'(I3.3)') NJP
          write(PROCNUM,'(I3.3)') process_id
          FNJPLOG = OUTDIR//'JPLOG' // FNNUM // '_' // PROCNUM // '.OUT'
          open(10,FILE = FNJPLOG, STATUS = 'UNKNOWN')
          close(10,STATUS = 'DELETE')
          
          if( JET_PLM(NJP).IOUTJP > 1 )then
            open(10,FILE = FNJPLOG, STATUS = 'UNKNOWN')
            write(10,630) 'NITER', 'TIMEDAY', 'NJE', 'NI', 'TOTTIME', 'LENGTH',                                                          &
                          'XJG(', 'YJG', 'ZJG', 'DIA', 'QVJET', 'SUM QJPENT', 'QJTOT', 'RHOA', 'RHOJ', 'TEMJ', 'VJ', 'VJH', 'WJG',       &
                          'WAG', 'RHOJ-RHOA', 'RHO_DELTA', 'RHO_DELTA2', 'TEMA', 'VA', 'VAH', 'RMAJP', 'DRMAJS', 'DRMAJF'
            close(10)
          endif
        endif
      enddo
630   FORMAT(A10,A15,A8,A5,2A9,2A12,A10,A8,3A12,8A10,10A12)  
      
    endif

    if( LDEBUG )then
      open(88,FILE = OUTDIR//'JPBUG.DIA', STATUS = 'UNKNOWN')
      close(88,STATUS = 'DELETE')
    endif
    
  endif   ! *** END OF INITIALIZATION

  ! *** SET CONSTANTS
  TIME = TIMESEC/TCON
  
  if( CURDAY < TIMEDAY )then
    CURDAY = CURDAY + 1./24.
    NPRT = 0
  endif

  ! *** G   = EFDC+ global constant for acceleration due to gravity
  ! *** RPI = 3.14159
  RADDEG = PI/180.
  
  IFILE = -1
  if( LDEBUG )then
    open(88,FILE = OUTDIR//'JPBUG.DIA',STATUS = 'UNKNOWN',POSITION = 'APPEND')
    close(88,STATUS = 'DELETE')
    write(88,*) TIMEDAY
  endif
  
  ! **************************************************************************************************
  ! **************************************************************************************************
  ! *** Prepare for plume calculations

  ! *** LOOP OVER ALL JET-PLUME LOCATIONS
  do NJP = 1,NQJPIJ
    ! *** Inititalize and prepare for the plume calculatopn for the current jet NJP
    do K = 1,KC
      QJPENT(K,NJP) = 0.0
      UJPAVG(K,NJP) = 0.0
      VJPAVG(K,NJP) = 0.0
      WJPAVG(K,NJP) = 0.0
    enddo

    if( NITER == 1 )then
      KEFFJP(NJP) = JET_PLM(NJP).KQJP               ! *** Initialize the plume discharge elevation to the nominal elevation
    endif
    
    if( JET_PLM(NJP).IOUTJP > 0 )then
      write(FNNUM,'(I3.3)') NJP
      write(PROCNUM,'(I3.3)') process_id
      FNJPLOG = OUTDIR//'JPLOG' // FNNUM // '_' // PROCNUM // '.OUT'
      open(10,FILE = FNJPLOG,STATUS = 'UNKNOWN',POSITION = 'APPEND')
      if( JET_PLM(NJP).IOUTJP > 1 ) WRITE(10,134) NJP, TIMEDAY
    endif
    
    if( LDEBUG )then
      open(88,FILE = OUTDIR//'JPBUG.DIA', STATUS = 'UNKNOWN',POSITION = 'APPEND')      
      write(88,134) NJP, TIMEDAY
    endif
    
    LJP   = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)         ! *** L index of the receiving/returning water
    KL    = KSZ(LJP)                         ! *** Current layer
    RKC   = FLOAT(KC - KL + 1)               ! *** Number of layers from surface
    KJP   = JET_PLM(NJP).KQJP                        ! *** Default/initial discharge layer
    KFLAG = 0
    
    ZTMP = (JET_PLM(NJP).ZJET - BELV(LJP))/HP(LJP)   ! *** Fraction of the water column from the bottom to the port elevation
    KJP  = NINT(RKC*ZTMP)
    if( KJP < KSZ(LJP) ) KJP = KSZ(LJP)
    if( KJP > KC ) KJP = KC
    
    ! *** ICAL > 0 - JET-PLUME BC ACTIVATED, OTHERWISE IGNORED
    if( JET_PLM(NJP).ICALJP > 0 )then
      !    NJEL     MAXIMUM NUMBER OF ELEMENTS ALONG JET-PLUME LENGTH
      !    ISAMB  0 FOR SPATIALLY AND TEMPORALLY CONSTANT AMBIENT CONDITIONS
      !           1 FOR SPATIALLY VARYING AND TEMPORALLY CONSTANT CONDITIONS
      !           2 FOR SPATIALLY AND TEMPORALLY VARYING AMBIENT CONDITIONS
      !    ISJPD  0 FOR TEMPORALLY CONSTANT JET-PLUME DISCHARGE
      !           1 FOR TEMPORALLY VARYING JET-PLUME DISCHARGE
      !    ISENT  0 use MAXIMUM OF SHEAR AND FORCED ENTRAINMENT
      !           1 use SUM OF SHEAR AND FORCED ENTRAINMENT
      !    ISTOP  0 STOP AT SPECIFIED NUMBER OF ELEMENTS
      !           1 STOP WHEN CENTERLINE PENETRATES BOTTOM OR SURFACE
      !           2 STOP WITH BOUNDARY PENETRATES BOTTOM OR SURFACE
      !    ISOUT  0 DIMENSIONAL OUTPUT,
      !           1 NONDIM OUTPUT LENGTH SCALED BY DJET
      !           2 NONDIM OUTPUT LENGTH SCALED BY SQRT(PI)*RJET
      ISAMB = 1         
      ISJPD = 0         
      ISTOP = JET_PLM(NJP).ISTJP
      ISOUT = 0         
      NPRTE = 0              ! *** NPRTE    Element output print frequency
      ALPHA = JET_PLM(NJP).CFRD      ! *** CFRD     Adjustment factor for froude number
      ALPH2 = ALPHA*ALPHA
      TUJP = 0               ! *** TUJP     Temporal frequency for updating jet-plume (sec) if ISJPD = 1
      
      NIMAX = JET_PLM(NJP).NJPMX     ! *** NIMAX    Maximum number of iterations
      XERRM = 1000.          ! *** XYERR    Horizontal trajectory error critera (m)
      YERRM = 1000.           
      ZERRM = 1000.          ! *** ZERR     Vertical trajectory error critera (m)
      RRERM = 1000.          ! *** RRER     Horizontal trajectory error critera (m)
      RLERM = 1000.          ! *** RLER     Vertical trajectory error critera (m)
      DMRERM = JET_PLM(NJP).DJPER    ! *** DMRERR   Entrainment error criteria
      
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
      RJET = 0.5*JET_PLM(NJP).DJET    ! *** RJET   Radius of discharge port (m),   
                              ! *** DJET   Diameter of discharge port (m)
      AJET = PI*RJET*RJET     ! *** This is area of 1 port but the flows used below are total flow into the cell

      ! *** Total Flow over all layers (QSER Type, ICAL = 1 )
      if( JET_PLM(NJP).ICALJP == 1 )then
        QVJET = JET_PLM(NJP).QQCJP      
        QSERTAVG = 0.
        do K = 1,KC
          QVJET =    QVJET    + QSERT(K,JET_PLM(NJP).NQSERJP)
          QSERTAVG = QSERTAVG + QSERT(K,JET_PLM(NJP).NQSERJP)
        enddo
      endif

      ! *** Total Flow over all layers (Withdrawal/Return type, ICAL = 2 )
      if( JET_PLM(NJP).ICALJP == 2 )then
        QVJET = JET_PLM(NJP).QWRCJP + QWRSERT(JET_PLM(NJP).NQWRSERJP)
      endif
      if( QVJET <= 1.E-16 ) GOTO 9000   ! *** NO FLOW, SKIP CALCS
      
      ! *** THE QSER OR NQWRSERJP FLOWS ARE ON A PER PORT BASIS.  CALQVS ADJUSTS PORT FLOWS FOR NUMBER OF PORTS IN CELL.  [ ALTERNATE: QVJET = QVJET/RPORTS(NJP) ]
      QVJET0 = QVJET      
      VJET = QVJET/AJET   ! *** Velocity of jet in port (m/s)
      
      ! *** Get jet-plume mass loading from QSER or W/R 
      !    SALJET   SALINITY (PPT)
      !    SEDJET   SEDIMENT (MG/L OR G/M^3)
      !    TEMJET   TEMPERATURE
      !    TOXJET   CONTAMINANT CONCENTRATION (UG/L OR MG/M**3)
      !    WQVJET   WATER QUALITY CONCENTRATION (G/M3)
      
      ! *** STANDARD FLOW (QSER) TYPE (ICAL = 1)
      if( JET_PLM(NJP).ICALJP == 1 )then
        WTC = JET_PLM(NJP).QQCJP      ! *** Constant flow component
        WTV = QSERTAVG                ! *** Variable flow component
        TMPVAL = WTC + WTV
        WTC = WTC/TMPVAL              ! *** Fraction of total flow for constant flows
        WTV = WTV/TMPVAL              ! *** Fraction of total flow for variable flows

        ! ***            Constant                              Series
        SALJET = WTC*JET_PLM(NJP).CQCJP(1,1) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(1),1)
        TEMJET = WTC*JET_PLM(NJP).CQCJP(1,2) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(2),2)

        do MD = 1,NDYE
          NV = MSVDYE(MD)
          DYEJET(MD) = WTC*JET_PLM(NJP).CQCJP(1,NV) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(3),NV)
        enddo
        
        NV = 3 + NDYM
        SFLJET = WTC*JET_PLM(NJP).CQCJP(1,4) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(4),NV)

        do NT = 1,NTOX
          NV = MSVTOX(NT)
          TOXJET(NT) = WTC*JET_PLM(NJP).CQCJP(1,NV) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(5),NV)
        enddo
        do NX = 1,NSED     
          NV = MSVSED(NX)
          SEDJET(NX) = WTC*JET_PLM(NJP).CQCJP(1,NV) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(6),NV)
        enddo
        do NX = 1,NSND
          NV = MSVSND(NX)
          SNDJET(NX) = WTC*JET_PLM(NJP).CQCJP(1,NV) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(7),NV)
        enddo
        do NX = 1,NWQV
          NV = MSVWQV(NX)
          WQVJET(NX) = WTC*JET_PLM(NJP).CQCJP(1,NV) + WTV*CSERT(1,JET_PLM(NJP).NCSERJP(8),NV)
        enddo
      endif
      
      ! *** WITHDRAWAL/RETURN TYPE (ICAL = 2)
      if( JET_PLM(NJP).ICALJP == 2 )then
        WTC = JET_PLM(NJP).QWRCJP                          ! *** Current withdrawal/return flow
        WTV = QWRSERT(JET_PLM(NJP).NQWRSERJP)              ! *** Times series number
        TMPVAL = WTC + WTV                                 
        WTC = WTC/TMPVAL                                   ! *** Fraction of total flow for constant flows
        WTV = WTV/TMPVAL                                   ! *** Fraction of total flow for variable flows
        NS = JET_PLM(NJP).NQWRSERJP                        ! *** Withdrawal/Return series ID
        LU = LIJ(JET_PLM(NJP).IUPCJP,JET_PLM(NJP).JUPCJP)  ! *** L Index at the withdrawal cell
        KU = JET_PLM(NJP).KUPCJP                           ! *** Layer at withdrawal cell
        
        ! ***            Constant                   Series
        SALJET = WTC*JET_PLM(NJP).CWRCJP(1) + WTV*CQWRSERT(NS,1) + SAL1(LU,KU)
        TEMJET = WTC*JET_PLM(NJP).CWRCJP(2) + WTV*CQWRSERT(NS,2) + TEM1(LU,KU)
        
        do MD = 1,NDYE
          NV = MSVDYE(MD)
          DYEJET(MD) = WTC*JET_PLM(NJP).CWRCJP(NV) + WTV*CQWRSERT(NS,NV) + DYE1(LU,KU,MD)
        enddo

        NV = 3 + NDYM
        SFLJET = WTC*JET_PLM(NJP).CWRCJP(NV) + WTV*CQWRSERT(NS,NV) + SFL2(LU,KU)

        do NT = 1,NTOX
          NV = MSVTOX(NT)
          TOXJET(NT) = WTC*JET_PLM(NJP).CWRCJP(NV) + WTV*CQWRSERT(NS,NV) + TOX1(LU,KU,NT)
        enddo
        do NX = 1,NSED
          NV = MSVSED(NX)
          SEDJET(NX) = WTC*JET_PLM(NJP).CWRCJP(NV) + WTV*CQWRSERT(NS,NV) + SED1(LU,KU,NX)
        enddo
        do NX = 1,NSND
          NV = MSVSND(NX)
          SNDJET(NX) = WTC*JET_PLM(NJP).CWRCJP(NV) + WTV*CQWRSERT(NS,NV) + SND1(LU,KU,NX)
        enddo
        do NW = 1,NWQV
          NV = MSVWQV(NW)
          WQVJET(NW) = WTC*JET_PLM(NJP).CWRCJP(NV) + WTV*CQWRSERT(NS,NV) + WQV(LU,KU,NW)
        enddo
      endif
      
      ZBOT = BELV(LJP)              ! *** ZBOT   BOTTOM ELEVATION (M)
      ZSUR = BELV(LJP) + HP(LJP)    ! *** ZSUR   SURFACE ELEVATION (M)
                                    ! *** KC     NUMBER OF AMBIENT DATA IN VERTICAL

      ! *** Spatial print/save interval in vertical
      DZPRT = HP(LJP)/FLOAT(JET_PLM(NJP).NZPRJP - 2)
      ZJPRT = JET_PLM(NJP).ZJET + DZPRT

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
      do K = KSZ(LJP)+1,KC
        ZAD(K) = ZAD(K-1) + 0.5*(HPK(LJP,K-1) + HPK(LJP,K))
      enddo
      
      ! *** Get layer specific cell centroid velocities
      do K = 1,KC
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
      enddo
889   FORMAT(5I6,3E14.5)
      
      ! *** Initialize cell water column concentrations (xxxAD) for receiving/returning water cell
      do K = 1,KC
        SALAD(K) = SAL(LJP,K)
        TEMAD(K) = TEM(LJP,K)
        SFLAD(K) = SFL(LJP,K)
      enddo
      do MD = 1,NDYE
        do K = 1,KC
          DYEAD(K,NT) = DYE(LJP,K,MD)
        enddo
      enddo
      do NT = 1,NTOX
        do K = 1,KC
          TOXAD(K,NT) = TOX(LJP,K,NT)
        enddo
      enddo
      do NS = 1,NSED2
        do K = 1,KC
          SEDAD(K,NS) = SED(LJP,K,NS)
        enddo
      enddo
      do NX = 1,NSND
        do K = 1,KC
          SNDAD(K,NX) = SND(LJP,K,NX)
        enddo
      enddo
      do NW = 1,NWQV
        do K = 1,KC
          WQVAD(K,NT) = WQV(LJP,K,NW)
        enddo
      enddo  

      ! *** Set initial conditions for current plume
      DTJP = 0.1*RJET/VJET    ! *** 1/S

      ! ***  JET-PLUME DISCHARGE, GLOBAL COORDINATES
      XJG(1) = XJET                  ! *** X coordinate of starting location   - Element 1
      YJG(1) = YJET                  ! *** Y coordinate of starting location   - Element 1
      ZJG(1) = JET_PLM(NJP).ZJET     ! *** Elevation of starting location      - Element 1
      XJG(2) = XJET                  ! *** X coordinate of starting location   - Element 2
      YJG(2) = YJET                  ! *** Y coordinate of starting location   - Element 2
      ZJG(2) = JET_PLM(NJP).ZJET     ! *** Elevation of starting location      - Element 2
      SIG(1) = 0.                    ! *** Total plume centerline displacement - Element 1
      SIG(2) = 0.                    ! *** Total plume centerline displacement - Element 2

      ! ***  INITIALIZE JET DISCHARGE GEOMETRY FOR THE COMPUTATIONAL ELEMENTS
      RADJ(1) = RJET                 ! *** Radius of the jet                        - Element 1
      RADJ(2) = RJET                 ! *** Radius of the jet                        - Element 2
      RADJOLD(1) = RJET              ! *** Radius of the jet of previous iteration  - Element 1
      RADJOLD(2) = RJET              ! *** Radius of the jet of previous iteration  - Element 2
      RLEJ(1) = 0.1*RJET             ! *** 0.1*Radius = Element length of the plume - Element 1
      RLEJ(2) = 0.1*RJET             ! *** 0.1*Radius = Element length of the plume - Element 2
      RADLAST = RJET                 ! *** Initializes write tolerance
      
      ! *** JET DISCHARGE VELOCITY MAGNITUDE AND COMPONENTS
      ! *** PHJET IS THE VERTICAL ANGLE FROM HORIZONTAL
      ! *** THJET IS THE HORIZONTAL ANGLE
      AJET = PI*RJET*RJET

      ! *** INITIALIZE ELEMENT 1 VELOCITIES
      VJ(1)  = QVJET/AJET                               ! *** Jet centerline velocity - Total
      VJH(1) = VJ(1)*COS(0.0175*JET_PLM(NJP).PHJET)     ! *** Jet centerline velocity - Horizontal only
      UJG(1) = VJ(1)*COS(0.0175*JET_PLM(NJP).PHJET)*COS(0.0175*JET_PLM(NJP).THJET)
      VJG(1) = VJ(1)*COS(0.0175*JET_PLM(NJP).PHJET)*SIN(0.0175*JET_PLM(NJP).THJET)
      WJG(1) = VJ(1)*SIN(0.0175*JET_PLM(NJP).PHJET)
      UJ0   = UJG(1)
      VJ0   = VJG(1)
      WJ0   = WJG(1)
      VJET0 = VJ(1)
            
      ! *** INITIALIZE ELEMENT 2 VELOCITIES
      VJ(2) = QVJET/AJET
      VJH(2) = VJ(2)*COS(0.0175*JET_PLM(NJP).PHJET)
      UJG(2) = VJ(2)*COS(0.0175*JET_PLM(NJP).PHJET)*COS(0.0175*JET_PLM(NJP).THJET)
      VJG(2) = VJ(2)*COS(0.0175*JET_PLM(NJP).PHJET)*SIN(0.0175*JET_PLM(NJP).THJET)
      WJG(2) = VJ(2)*SIN(0.0175*JET_PLM(NJP).PHJET)
      
      DTJP = RLEJ(1)/VJ(1)                    ! *** Initialize time step (s)
      TOTTIME = 0.                            ! *** Initialize total elapsed time of plume development

      ! ***  INITIAL JET DENSITY AND MASS
      SEDJETT = 0.
      do NS = 1,NSED
        SEDJETT = SEDJETT + SEDJET(NS)
      enddo
      do NX = 1,NSND
        SEDJETT = SEDJETT + SNDJET(NX)
      enddo
      RHOJ(1) = FUNDEN(SALJET,SEDJETT,TEMJET)         ! *** Jet density
      RMAJP(1) = PI*RHOJ(1)*RADJ(1)*RADJ(1)*RLEJ(1)   ! *** Jet mass

      ! ***  INITIAL JET CONCENTRATIONS
      SALJ(1) = SALJET
      TEMJ(1) = TEMJET
      do MD = 1,NDYE
        DYEJ(1,MD) = DYEJET(MD)
      enddo
      SFLJ(1) = SFLJET
      do NT = 1,NTOX
        TOXJ(1,NT) = TOXJET(NT)
      enddo
      do NS = 1,NSED
        SEDJ(1,NS) = SEDJET(NS)
      enddo
      do NX = 1,NSND
        SNDJ(1,NX) = SNDJET(NX)
      enddo
      do NW = 1,NWQV
        WQVJ(1,NW) = WQVJET(NW)
      enddo

      ! ***  INITIAL JET TOXIC CONTAMINANT PARTICULATE FRACTIONS
      do NT = 1,NTOX
        if( ISTRAN(6) >= 1 )then
          do NS = 1,NSED
            TMPEXP = CONPARW(NS,NT)
            if( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
            if( ITXPARW(NS,NT) == 1 )then
              if( SEDJET(NS) > 0.) TMPVAL = SEDJET(NX)**TMPEXP
            endif
            TOXPFJP(1,NS,NT) = TMPVAL*SEDJET(NS)*TOXPARW(LJP,NS,NT)
          enddo
        endif
        if( ISTRAN(7) >= 1 )then
          do NX = 1,NSND
            NS = NX + NSED
            TMPEXP = CONPARW(NS,NT)
            if( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
            if( ITXPARW(NS,NT) == 1 )then
              if( SNDJET(NX) > 0.) TMPVAL = SNDJET(NX)**TMPEXP
            endif
            TOXPFJP(1,NS,NT) = TMPVAL*SNDJET(NX)*TOXPARW(LJP,NS,NT)
          enddo
        endif
      enddo
      do NT = 1,NTOX
        TOXPFTJP(1,NT) = 0.
      enddo
      do NT = 1,NTOX
        if( ISTRAN(6) >= 1 )then
          do NS = 1,NSED
            TOXPFTJP(1,NT) = TOXPFTJP(1,NT) + TOXPFJP(1,NS,NT)
          enddo
        endif
        if( ISTRAN(7) >= 1 )then
          do NX = 1,NSND
            NS = NX + NSED
            TOXPFTJP(1,NT) = TOXPFTJP(1,NT) + TOXPFJP(1,NS,NT)
          enddo
        endif
      enddo
      do NT = 1,NTOX
        do NS = 1,NSED + NSND
          TOXPFJP(1,NS,NT) = TOXPFJP(1,NS,NT)/(1. + TOXPFTJP(1,NT))
        enddo
      enddo
      do NT = 1,NTOX
        TOXPFTJP(1,NT) = TOXPFTJP(1,NT)/(1. + TOXPFTJP(1,NT))
      enddo

      ! ***  INITIALIZE JET ELEMENT 2 VALUE TO ELEMENT 1 VALUES
      RHOJ(2) = FUNDEN(SALJET,SEDJETT,TEMJET)            ! *** Density 
      RMAJP(2) = PI*RHOJ(2)*RADJ(2)*RADJ(2)*RLEJ(2)      ! *** Mass (density * volume)
      SALJ(2) = SALJ(1)
      TEMJ(2) = TEMJ(1)
      RHO_INI = RHOJ(1)                                  ! *** Initial plume density
      PTINI = TEMJ(1)                                    ! *** Initial plume temperature
      do MD = 1,NDYE
        DYEJ(2,MD) = DYEJ(1,MD)
      enddo
      SFLJ(2) = SFLJ(1)
      do NT = 1,NTOX
        TOXJ(2,NT) = TOXJ(1,NT)
      enddo
      do NT = 1,NTOX
        do NS = 1,NSED + NSND
          TOXPFJP(2,NS,NT) = TOXPFJP(1,NS,NT)
        enddo
      enddo
      do NT = 1,NTOX
        TOXPFTJP(2,NT) = TOXPFTJP(1,NT)
      enddo
      do NS = 1,NSED
        SEDJ(2,NS) = SEDJ(1,NS)
      enddo
      do NX = 1,NSND
        SNDJ(2,NX) = SNDJ(1,NX)
      enddo
      do NW = 1,NWQV
        WQVJ(2,NW) = WQVJ(1,NW)
      enddo

      ! ***  INITIALIZE AMBIENT CONDITIONS
      if( ISAMB == 0 )then
        ! *** CONSTANT AMBIENT CONDITIONS (LAYER = JET_PLM(NJP).KQJP)
        K = JET_PLM(NJP).KQJP
        UAG = UAGD(K)
        VAG = VAGD(K)
        WAG = WAGD(K)
        SALA = SALAD(K)
        TEMA = TEMAD(K)
        do MD = 1,NDYE
          DYEA(MD) = DYEAD(K,MD)
        enddo
        SFLA = SFLAD(K)
        do NT = 1,NTOX
          TOXA(NT) = TOXAD(K,NT)
        enddo
        do NS = 1,NSED2
          SEDA(NS) = SEDAD(K,NS)
        enddo
        do NX = 1,NSND
          SNDA(NX) = SNDAD(K,NX)
        enddo
        do NW = 1,NWQV
          WQVA(NW) = WQVAD(K,NW)
        enddo
      elseif( ISAMB >= 1 )then
        ! *** SPATIALLY VARIABLE AMBIENT CONDITIONS
        ZVAL = ZJG(1)
        call JPACON(KSZ(LJP),ZVAL,UAG,VAG,WAG)
      endif

      ! ***  AMBIENT VELOCITY MAGNITUDES
      VA = SQRT(UAG*UAG + VAG*VAG + WAG*WAG)
      VAH = SQRT(UAG*UAG + VAG*VAG)

      ! ***  AMBIENT DENSITY
      SEDATT = 0.
      do NS = 1,NSED2
        SEDATT = SEDATT + SEDA(NS)
      enddo
      do NX = 1,NSND
        SEDATT = SEDATT + SNDA(NX)
      enddo
      RHOA(1) = FUNDEN(SALA,SEDATT,TEMA)
      RHOA(2) = RHOA(1)
      RHO_ELE = RHOA(1)
      DELTEM = PTINI - TEMA
      DELRHO = RHO_INI - RHO_ELE
      EXITTEM = 0.05*ABS(DELTEM)          ! *** Exit criteria based on initial temperature delta
      EXITRHO = 0.05*ABS(DELRHO)          ! *** Exit criteria based on initial density delta

      ! ***  Global Ambient and Jet orientations
      PHJ(1) = JET_PLM(NJP).PHJET
      THJG(1) = JET_PLM(NJP).THJET
      PHJ(2) = JET_PLM(NJP).PHJET
      THJG(2) = JET_PLM(NJP).THJET
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
      if( DVEL > 0. )then
        SINPHJ = SIN(0.0175*PHJ(1))
        FRD2I = FTOP/(ALPH2*DVEL*DVEL)                         ! *** Relative density Froude Number
        EBOT = 1. + 5.*RDQA(1)/DVEL
        ETOP = 0.057 + 0.554*SINPHJ*FRD2I
        ENTS = 1.414*ETOP/EBOT
      endif
      
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
      if( ISOUT == 1 ) RLSCL = JET_PLM(NJP).DJET
      if( ISOUT == 2 ) RLSCL = RJET*SQRT(PI)
      SALJETI = 1.
      TEMJETI = 1.
      do MD = 1,NDYE
        DYEJETI(MD) = 1.
      enddo
      SFLJETI = 1.
      do NT = 1,NTOX
        TOXJETI(NT) = 1.
      enddo
      do NS = 1,NSED
        SEDJETI(NS) = 1.
      enddo
      do NX = 1,NSND
        SNDJETI(NX) = 1.
      enddo
      do NW = 1,NWQV
        WQVJETI(NW) = 1.
      enddo
      
      ! *** APPLY SCALING FACTOR 
      if( ISOUT >= 1 )then
        if( SALJET > 0.) SALJETI = 1./SALJET
        if( TEMJET > 0.) TEMJETI = 1./TEMJET
        do MD = 1,NDYE
          if( DYEJET(MD) > 0.) DYEJETI(MD) = 1./DYEJET(MD)
        enddo
        if( SFLJET > 0.) SFLJETI = 1./SFLJET
        do NT = 1,NTOX
          if( TOXJET(NT) > 0.) TOXJETI(NT) = 1./TOXJET(NT)
        enddo
        do NS = 1,NSED
          if( SEDJET(NS) > 0.) SEDJETI(NS) = 1./SEDJET(NS)
        enddo
        do NX = 1,NSND
          if( SNDJET(NX) > 0.) SNDJETI(NX) = 1./SNDJET(NX)
        enddo
        do NW = 1,NWQV
          if( WQVJET(NW) > 0.) WQVJETI(NW) = 1./WQVJET(NW)
        enddo
      endif
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
      do MD = 1,NDYE
        DYEJNE(MD) = DYEJETI(MD)*DYEJ(N2,MD)
      enddo
      SFLJNE = TEMJETI*SFLJ(N2)
      do NT = 1,NTOX
        TOXJNE(NT) = TOXJETI(NT)*TOXJ(N2,NT)
      enddo
      do NS = 1,NSED
        SEDJNE(NS) = SEDJETI(NS)*SEDJ(N2,NS)
      enddo
      do NX = 1,NSND
        SNDJNE(NX) = SNDJETI(NX)*SNDJ(N2,NX)
      enddo
      do NW = 1,NWQV
        WQVJNE(NW) = WQVJETI(NW)*WQVJ(N2,NW)
      enddo
      
      DRMAJ = RMAJP(2) - RMAJP(1)                ! *** Mass of newly entrained water at the current iteration
      DRHO = (RHOA(1) - RHOJ(N2))/RHOA(1)        ! *** Relative density diffference between ambient conditions and plume

      ! **************************************************************************************************
      ! **************************************************************************************************
      ! *** Plume entrainment loop for NJEL elements until either max elements or plume hits end point
      do NJE = 2,JET_PLM(NJP).NJEL

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
1000    continue

        ! *** Calculate shear entrainment
        DRMAJSO = 0.5*(DRMAJS(N2) + DRMAJS(N1))      ! *** Average shear between N2 and N1
        DVEL = ABS(VJ(N2) - RDQA(N2))                ! *** Velocity shear (m/s):  VJ (3D centerline velocity)  RDQA (difference between jet and ambient velocities)
        DRHO = (RHOA(N2) - RHOJ(N2))/RHOA(N2)        ! *** Relative density between the jet and the ambient waters
        FTOP = G*ABS(DRHO)*RADJ(N2)                  ! *** Density momentum
        FRD2I = 0.
        EBOT = 1.
        ETOP = 0.057
        ENTS = 0.
        if( DVEL > 0. )then
          FRD2I = FTOP/(ALPH2*DVEL*DVEL)
          EBOT = 1. + 5.*RDQA(1)/DVEL 
          ETOP = 0.057 + 0.554*SINPHJ*FRD2I
          ENTS = 1.414*ETOP/EBOT
        endif
        
        ! ***             s     kg/m3    (-)   m/s    m       m
        DRMAJS(N2) = 2.*DTJP*PI*RHOA(N2)*ENTS*DVEL*RADJ(N2)*RLEJ(N2)     ! *** Shear at the end of the current element time step (DTJP) (kg)
        DRMAJSA = 0.5*(DRMAJS(N2) + DRMAJS(N1))                          ! *** New average shear between N2 and N1

        ! ***  Calculate forced entrainment (hardwired)
        DRMAJFO = 0.5*(DRMAJF(N2) + DRMAJF(N1))

        ENTF1 = 2.*SQRT( 1. - COSPHJ*COSPHJ*COSTHJL*COSTHJL )
        DELSIG = SIG(N2) - SIG(N1)
        ENTF2 = 0.
        ENTF3 = 0.
        if( DELSIG > 0. )then
          ENTF2 = PI*COSPHJ*COSTHJL*(RADJ(N2) - RADJ(N1))/(DELSIG)
          ENTF3 = 0.5*PI*RADJ(N2)*(COSPHJ*COSTHJL - COSPHJM*COSTHJLM)/(DELSIG)
        endif
        ENTF = ENTF1 + ENTF2 + ENTF3
        ENTF = MAX(ENTF,0.)
        DRMAJF(N2) = DTJP*RHOA(N2)*RADJ(N2)*RLEJ(N2)*VAH*ENTF
        if( NJE == 2 .and. NI == 1 )DRMAJF(N2) = 0.
        DRMAJFA = 0.5*(DRMAJF(N2) + DRMAJF(N1))
        
        ISHEAR = 0
        IFORCE = 0
        if( JET_PLM(NJP).ISENT == 0 )then
          ! ***  TAKE MAX OF SHEAR AND FORCED
          DRMAJ = MAX(DRMAJSA,DRMAJFA)
        else
          ! ***  TAKE SUM OF SHEAR AND FORCED
          DRMAJ = DRMAJSA + DRMAJFA
        endif
        if( DRMAJSA > DRMAJFA) ISHEAR = 1
        if( DRMAJFA > DRMAJSA) IFORCE = 1

        ! *** Preparations for mass entrainment
        RMAJP(N2) = RMAJP(N1) + DRMAJ       ! *** Total volume at the end of the current time step
        RMAJI = 1./RMAJP(N2)

        ! *** Entrain mass based on the total shear (DRMAJ) and compute new concentrations
        SALJ(N2)      = RMAJI*( RMAJP(N1)*SALJ(N1)    + DRMAJ*SALA )
        TEMJ(N2)      = RMAJI*( RMAJP(N1)*TEMJ(N1)    + DRMAJ*TEMA )
        do MD = 1,NDYE
          DYEJ(N2,MD) = RMAJI*( RMAJP(N1)*DYEJ(N1,MD) + DRMAJ*DYEA(MD) )
        enddo           
        SFLJ(N2)      = RMAJI*( RMAJP(N1)*SFLJ(N1)    + DRMAJ*TEMA )
        do NT = 1,NTOX    
          TOXJ(N2,NT) = RMAJI*( RMAJP(N1)*TOXJ(N1,NT) + DRMAJ*TOXA(NT) )
        enddo           
        do NS = 1,NSED   
          SEDJ(N2,NS) = RMAJI*( RMAJP(N1)*SEDJ(N1,NS) + DRMAJ*SEDA(NS) )
        enddo           
        do NX = 1,NSND    
          SNDJ(N2,NX) = RMAJI*( RMAJP(N1)*SNDJ(N1,NX) + DRMAJ*SNDA(NX))
        enddo           
        do NW = 1,NWQV    
          WQVJ(N2,NW) = RMAJI*( RMAJP(N1)*WQVJ(N1,NW) + DRMAJ*WQVA(NW) )
        enddo

        ! *** ADVANCE TOXIC PARTICULATE FRACTION
        if( ISTRAN(5) >= 1 )then
          do NT = 1,NTOX
            if( ISTRAN(6) >= 1 )then
              do NS = 1,NSED
                TMPEXP = CONPARW(NS,NT)
                if( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
                if( ITXPARW(NS,NT) == 1 )then
                  if( SEDJ(N2,NS) > 0.) TMPVAL = SEDJ(N2,NS)**TMPEXP
                endif
                TOXPFJP(N2,NS,NT) = TMPVAL*SEDJ(N2,NS)*TOXPARW(LJP,NS,NT)
              enddo
            endif
            if( ISTRAN(7) >= 1 )then
              do NX = 1,NSND
                NS = NX + NSED
                TMPEXP = CONPARW(NS,NT)
                if( ITXPARW(NS,NT) == 0 ) TMPVAL = 1.
                if( ITXPARW(NS,NT) == 1 )then
                  if( SNDJ(N2,NX) > 0.) TMPVAL = SNDJ(N2,NX)**TMPEXP
                endif
                TOXPFJP(N2,NS,NT) = TMPVAL*SNDJ(N2,NX)*TOXPARW(LJP,NS,NT)
              enddo
            endif
          enddo
          do NT = 1,NTOX
            TOXPFTJP(N2,NT) = 0.
          enddo
          do NT = 1,NTOX
            if( ISTRAN(6) >= 1 )then
              do NS = 1,NSED
                TOXPFTJP(N2,NT) = TOXPFTJP(N2,NT) + TOXPFJP(N2,NS,NT)
              enddo
            endif
            if( ISTRAN(7) >= 1 )then
              do NX = 1,NSND
                NS = NX + NSED
                TOXPFTJP(N2,NT) = TOXPFTJP(N2,NT) + TOXPFJP(N2,NS,NT)
              enddo
            endif
          enddo
          do NT = 1,NTOX
            do NS = 1,NSED + NSND
            TOXPFJP(N2,NS,NT) = TOXPFJP(N2,NS,NT)/(1. + TOXPFTJP(N2,NT))
            enddo
          enddo
          do NT = 1,NTOX
            TOXPFTJP(N2,NT) = TOXPFTJP(1,NT)/(1. + TOXPFTJP(N2,NT))
          enddo
        endif
        
        ! *** ADVANCE DENSITY
        SEDJETT = 0.
        do NS = 1,NSED
          SEDJETT = SEDJETT + SEDJ(N2,NS)
        enddo
        do NX = 1,NSND
          SEDJETT = SEDJETT + SNDJ(N2,NX)
        enddo
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
        if( ISAMB == 0 )then
          K = JET_PLM(NJP).KQJP
          UAG = UAGD(K)
          VAG = VAGD(K)
          WAG = WAGD(K)
          SALA = SALAD(K)
          TEMA = TEMAD(K)
          do MD = 1,NDYE
            DYEA(MD) = DYEAD(K,MD)
          enddo
          SFLA = SFLAD(K)
          do NT = 1,NTOX
            TOXA(NT) = TOXAD(K,NT)
          enddo
          do NS = 1,NSED2
            SEDA(NS) = SEDAD(K,NS)
            SEDAA = SEDAA + SEDA(NS)
          enddo
          do NX = 1,NSND
            SNDA(NX) = SNDAD(K,NX)
            SEDAA = SEDAA + SEDA(NX)
          enddo
          do NW = 1,NWQV
            WQVA(NW) = WQVAD(K,NW)
          enddo
        endif
        if( ISAMB >= 1 )then
          ! *** Get ambient conditions from the current plume conditions
          ZVAL = ZJG(N2)
          call JPACON(KSZ(LJP),ZVAL,UAG,VAG,WAG)
        endif
        
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
        if( DRMAJSO > 0. ) DRMAJSE = DRMAJSE/DRMAJSO
        if( DRMAJFO > 0. ) DRMAJFE = DRMAJFE/DRMAJFO
        ITMP = 0
        if( DRMAJFE > DMRERM ) ITMP = 1
        if( DRMAJSE > DMRERM ) ITMP = 1
        
        if( LDEBUG )then
          if( NJE == 2 ) WRITE(88,620)NJP, NJE, NI, ITMP, DRMAJSA, DRMAJSO, DRMAJFA, DRMAJFO
          if( JET_PLM(NJP).ISDJP == 1 )then
            ! *** ISDJP > 0 write diagnostics
            JTMPVAL = MOD(NJE,100)
            if( JTMPVAL == 0 )then
              write(88,620)NJP, NJE, NI, ITMP, DRMAJSA, DRMAJSO, DRMAJFA, DRMAJFO
            endif
          endif
        endif

        ! ***  STOP IF MAXIMUM ITERATIONS EXCEEDED
        if( NI > NIMAX )then
          KFLAG = 1
          if( LDEBUG )WRITE(88,620)NJP, NJE, NI, ITMP, DRMAJSA, DRMAJSO, DRMAJFA, DRMAJFO
          if( LDEBUG )WRITE(88,601)TIMEDAY, NJE, NI, LJP, HP(LJP)
          
          write(6,'(A,3I6,F8.3)' ) ' JET-PLUME ITERATIONS EXCEEDED N2,NI @ L,HP = ',NJE,NI,LJP,HP(LJP)
          
          if( LOUTJET )then
            write(10,'(A,I10,F15.5,3I6,F8.3)' ) ' JET-PLUME ITERATIONS EXCEEDED N2,NI @ L,HP = ',NITER,TIMEDAY,NJE,NI,LJP,HP(LJP)
          endif

          GOTO 1500
        endif
        
        ! *** Add criteria to exit entraiment calculations based on stagnant plume growth
        if( NI > 5 .and. RADJ(N2) <= RADJOLD(N2) .and. RADJ(N1) <= RADJOLD(N1) )then
          ITMP = 0          ! *** Flag causes loop to exit below
        endif

        ! *** Check for exit conditions based on current centerline or boundary
        ! ***  STOP if jet centerline penetrates surface
        if( ISTOP == 1 )then
          ZJGTOP = ZJG(N2)
          if( ZJGTOP > ZSUR )then
            NPRT(5) = NPRT(5) + 1
            if( NPRT(5) < 10 )then
              write(6,605) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            if( LOUTJET )then
              write(10,615) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif

            EXIT
          endif
        endif

        ! ***  STOP if jet centerline penetrates bottom
        if( ISTOP == 1 )then
          ZJGBOT = ZJG(N2)
          if( ZJGBOT < ZBOT )then
            NPRT(6) = NPRT(6) + 1
            if( NPRT(6) < 10 )then
              write(6,606) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            if( LOUTJET )then
              write(10,616) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            EXIT
          endif
        endif

        ! ***  STOP if jet boundary penetrates surface
        if( ISTOP == 2 )then
          ZJGTOP = ZJG(N2) + RADJ(N2)*COS(0.0175*PHJ(N2))
          if( ZJGTOP > ZSUR )then
            NPRT(6) = NPRT(6) + 1
            if( NPRT(6) < 10 )then
              write(6,602) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            if( LOUTJET )then
              write(10,612) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            EXIT
          endif
        endif

        ! ***  STOP if jet boundary penetrates bottom
        if( ISTOP == 2 )then
          ZJGBOT = ZJG(N2) - RADJ(N2)*COS(0.0175*PHJ(N2))
          if( ZJGBOT < ZBOT )then
            NPRT(7) = NPRT(7) + 1
            if( NPRT(7) < 10 )then
              write(6,603) NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            if( LOUTJET )then
              write(10,613) TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            EXIT
          endif
        endif

        ! ***  STOP if neutral level is reached
        ! ***  Rising plume
        if( NJE > 4 .and. RHOJ(N2) >= RHO_INI )then
          DRHOT = (RHOA(N2) - RHOJ(N2))/RHOA(N2)
          if( DRHOT < 0. )then
            NPRT(1) = NPRT(1) + 1
            if( NPRT(1) < 10 )then
              write(6,604) "R", NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            if( LOUTJET )then
              write(10,614)  "R", TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            EXIT
          endif
        endif

        ! ***  Falling plume
        if( NJE > 4 .and. RHOJ(N2) < RHO_INI )then
          DRHOT = (RHOA(N2) - RHOJ(N2))/RHOA(N2)
          if( DRHOT > 0. )then
            NPRT(2) = NPRT(2) + 1
            if( NPRT(2) < 10 )then
              write(6,604) "F", NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            if( LOUTJET )then
              write(10,614)  "F", TIMEDAY, NJP, NJE, NI, BELV(LJP), ZJG(N2), ZSUR
            endif
            
            EXIT
          endif
        endif

        ! ***  RETURN FOR ANOTHER ITERATION
        if( ITMP == 1 .or. NI == 1 )then
          ! XJOLD = XJG(N2)
          ! YJOLD = YJG(N2)
          ! ZJOLD = ZJG(N2)
          ! RLOLD = RLEJ(N2)
          RADJOLD = RADJ
          NI = NI + 1
          GOTO 1000
        endif
        
        ! *************************************
        ! *** Entrainment iterations complete
1500    continue        
        
        ! ***  WRITE OUTPUT AND PROCEED TO NEXT JET ELEMENT
        if( ZJG(N2) >= ZJPRT )then
          NPRTE = 1
          ZJPRT = ZJPRT + DZPRT
        else
          NPRTE = 0
        endif
        DJETI = 1./RLSCL
        XJGNE = DJETI*XJG(N2)
        YJGNE = DJETI*YJG(N2)
        ZJGNE = DJETI*ZJG(N2)
        SIGNE = DJETI*SIG(N2)
        RADJNE = DJETI*RADJ(N2)
        RLEJNE = DJETI*RLEJ(N2)
        SALJNE = SALJETI*SALJ(N2)
        TEMJNE = TEMJETI*TEMJ(N2)
        do MD = 1,NDYE
          DYEJNE(MD) = DYEJETI(MD)*DYEJ(N2,MD)
        enddo
        SFLJNE = SFLJETI*SFLJ(N2)
        do NT = 1,NTOX
          TOXJNE(NT) = TOXJETI(NT)*TOXJ(N2,NT)
        enddo
        do NS = 1,NSED2
          SEDJNE(NS) = SEDJETI(NS)*SEDJ(N2,NS)
        enddo
        do NX = 1,NSND
          SNDJNE(NX) = SNDJETI(NX)*SNDJ(N2,NX)
        enddo
        do NW = 1,NWQV
          WQVJNE(NW) = WQVJETI(NW)*WQVJ(N2,NW)
        enddo

        ! *** Compute the layer specific JET-PLUME entrainment
        QJTOTO = QJTOT
        QJTOT = QJTOT*RMAJP(N2)/RMAJP(N1)
        do K = KSZ(L),KC
          ZLOWER = Z(L,K-1)*HP(LJP) + BELV(LJP)
          ZUPPER = ZLOWER + DZC(L,K)*HP(LJP)
          if( ZJG(N2) >= ZLOWER .and. ZJG(N2) < ZUPPER )then
            QJPENT(K,NJP) = QJPENT(K,NJP) + (QJTOT - QJTOTO)
            if( RMAJI > 0.0 )then
              UJPAVG(K,NJP) = UJG(N2)/RMAJI
              VJPAVG(K,NJP) = VJG(N2)/RMAJI
              WJPAVG(K,NJP) = WJG(N2)/RMAJI
            endif
            EXIT
          endif
        enddo

        ! *** Added additional plume migration exit criteria
        !IF( ABS(RHOJ(N2)-RHOA(N2)) < EXITRHO )then
        ! 
        !  if( MOD(NJE,20) == 0 )then
        !    RHO_DELTA = ABS(RHOJ(N2) - RHOA(N1))
        !    if( RHO_DELTA2 == -999. )then
        !      RHO_DELTA2 = RHO_DELTA
        !    endif
        !  
        !    if( RHO_INI < RHO_ELE )then
        !      ! *** Rising plume
        !      if( RHO_DELTA > RHO_DELTA2 )then
        !        NPRT(3) = NPRT(3) + 1
        !        if( NPRT(3) < 10 )then
        !          PRINT '(A,I10,F15.4,I8,I5,4F8.3,A)', 'Rising Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), '  New Exit'
        !        endif
        !        
        !        if( LOUTJET )then
        !          !WRITE(10, '(A,I10,F15.4,I8,I5,4F8.3,3f10.3,A)') 'Rising Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), RADJ(N2), SIG(N2), ZJG(N2), '  New Exit'
        !        endif
        !    
        !        !EXIT                         ! *** Jet calculations diverging from convergence
        !      endif
        !    else
        !      ! *** Falling plume
        !      if( RHOJ(N2) < RHOA(N1) )then
        !        NPRT(4) = NPRT(4) + 1
        !        if( NPRT(4) < 10 )then
        !          PRINT '(A,I10,F15.4,I8,I5,4F8.3,A)', 'Falling Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), '  New Exit'
        !        endif
        !        
        !        if( LOUTJET )then
        !          write(10, '(A,I10,F15.4,I8,I5,4F8.3,A)') 'Rising Exit ',NITER, TIMEDAY, NJE, NI, ABS(RHOJ(N2)-RHOA(N2))/EXITRHO, abs(RHO_DELTA)/exitrho, TEMA, TEMJ(N2), ZJG(N2), SIG(N2), '  New Exit'
        !        endif
        !        
        !        EXIT                         ! *** Jet calculations diverging from convergence
        !      endif              
        !    endif
        !    RHO_DELTA2 = RHO_DELTA
        !    RHOA(N1)   = RHOA(N2)               ! *** Density of ambient water at face
        !    RHOJ(N1)   = RHOJ(N2)               ! *** Density of jet
        !  endif        
        !else
        !  nx = 0 ! DELME
        !ENDIF

        ! *** Write to file when diameter changes by 0.1 m
        IDUMP = 0
        if( JET_PLM(NJP).IOUTJP > 1 )then
          if( NJE == 1 .or. RADJ(N2)-RADLAST > 0.05 )then
            write(10,631) NITER, TIMEDAY, NJE, NI, TOTTIME, SIG(N2),                                                                          &
                          XJG(N2), YJG(N2), ZJG(N2), 2.*RADJ(N2), QVJET, SUM(QJPENT(:,NJP)), QJTOT, RHOA(N2), -999., RHOJ(N2), TEMJ(N2),      &
                          VJ(N2), VJH(N2), WJG(N2), WAG, RHOJ(N2)-RHOA(N2), RHO_DELTA, RHO_DELTA2, TEMA, VA, VAH, -999., RMAJP(N2), DRMAJS(N2), DRMAJF(N2)
            RADLAST = RADJ(N2)
            IDUMP = 1
          endif
        endif
        
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
        do MD = 1,NDYE
          DYEJ(N1,MD) = DYEJ(N2,MD)
        enddo
        SFLJ(N1) = SFLJ(N2)
        do NT = 1,NTOX
          TOXJ(N1,NT) = TOXJ(N2,NT)
        enddo
        do NS = 1,NSED2
          SEDJ(N1,NS) = SEDJ(N2,NS)
        enddo
        do NX = 1,NSND
          SNDJ(N1,NX) = SNDJ(N2,NX)
        enddo
        do NW = 1,NWQV
          WQVJ(N1,NW) = WQVJ(N2,NW)
        enddo
       
      enddo  ! *** End NJEL element integration

      ! *** Output results for last element
      if( IDUMP == 1 )then
            write(10,631) NITER, TIMEDAY, NJE, NI, TOTTIME, SIG(N2),                                                                          &
                          XJG(N2), YJG(N2), ZJG(N2), 2.*RADJ(N2), QVJET, SUM(QJPENT(:,NJP)), QJTOT, RHOA(N2), -999., RHOJ(N2), TEMJ(N2),      &
                          VJ(N2), VJH(N2), WJG(N2), WAG, RHOJ(N2)-RHOA(N2), RHO_DELTA, RHO_DELTA2, TEMA, VA, VAH, -999., RMAJP(N2), DRMAJS(N2), DRMAJF(N2)
      endif
631   FORMAT(I10,F15.4,I8,I5,2F9.3,2F12.4,F10.3,F8.2,3E12.4,8f10.3,10E12.4)  
    
      ! *** Get the total entrainment volumes for current NJP
      QJPENTT(NJP) = 0.0
      do K = KSZ(L),KC
        QJPENTT(NJP) = QJPENTT(NJP) + QJPENT(K,NJP)
      enddo

      if( LDEBUG )then
        ! *** OPEN LOF FILE
        if( IFILE == -1 )then
          IFILE = mpi_efdc_out_unit
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        endif
        write(mpi_efdc_out_unit,898)NJP,TIME,(QJPENT(K,NJP),K = 1,KC),QJPENTT(NJP)
      endif
      DJETI = 1./RLSCL
      XJGNE = DJETI*XJG(N2)
      YJGNE = DJETI*YJG(N2)
      ZJGNE = DJETI*ZJG(N2)
      SIGNE = DJETI*SIG(N2)
      RADJNE = DJETI*RADJ(N2)
      RLEJNE = DJETI*RLEJ(N2)
      SALJNE = SALJETI*SALJ(N2)
      TEMJNE = TEMJETI*TEMJ(N2)
      do MD = 1,NDYE
        DYEJNE(MD) = DYEJETI(MD)*DYEJ(N2,MD)
      enddo
      SFLJNE = SFLJETI*SFLJ(N2)
      do NT = 1,NTOX
        TOXJNE(NT) = TOXJETI(NT)*TOXJ(N2,NT)
      enddo
      do NS = 1,NSED2
        SEDJNE(NS) = SEDJETI(NS)*SEDJ(N2,NS)
      enddo
      do NX = 1,NSND
        SNDJNE(NX) = SNDJETI(NX)*SNDJ(N2,NX)
      enddo
      do NW = 1,NWQV
        WQVJNE(NW) = WQVJETI(NW)*WQVJ(N2,NW)
      enddo
      DRMAJ = RMAJP(N2) - RMAJP(N2-1)
      DRHO = (RHOA(N2) - RHOJ(N2))/RHOA(N2)

      ! *** DIAGNOSTIC OUPUT
      if( LOUTJET .and. (JET_PLM(NJP).IOUTJP == 1 .or. JET_PLM(NJP).IOUTJP == 3 ) )then

      endif

      ! *** DIAGNOSTIC OUPUT

      ! *** RELOCATE VOLUME SOURCE
      ZTMP = ( ZJGNE-BELV(LJP) )/HP(LJP)
      ZTMP = RKC*ZTMP
      KTMP = INT(ZTMP) + KL 
      KTMP = MAX(KL,KTMP)
      KTMP = MIN(KC,KTMP)
      if( KFLAG == 0 ) KEFFJP(NJP) = KTMP
      GOTO 9001

      ! *** JET-PLUME COMPUTATIONS BYPASSED
      9000 continue
      KEFFJP(NJP) = JET_PLM(NJP).KQJP
      
9001  continue
      
      if( LDEBUG )then
        ! *** OPEN LOF FILE
        if( IFILE == -1 )then
          IFILE = mpi_efdc_out_unit
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        endif
        write(mpi_efdc_out_unit ,899)NJP,TIME,(QJPENT(K,NJP),K = 1,KC)
        write(mpi_efdc_out_unit ,135)NJP,TIME,KFLAG,KEFFJP(NJP),JET_PLM(NJP).KQJP,QVJET,QJTOT
      endif

      ! *** CALCULATION MOMENT INTERFACE QUANTITIES
      RDUM = 0.
      if( LDEBUG )then
        write(11,1110)NJP,N
        KZERO = 0
        QUJ0 = UJ0*QVJET0
        QVJ0 = VJ0*QVJET0
        QWJ0 = WJ0*QVJET0
        write(11,1111)KZERO,QUJ0,QVJ0,QWJ0,QVJET0,RDUM,RDUM,RDUM,RDUM
        QENTTMP = QVJET0
        do K = 1,KEFFJP(NJP)
          QENTTMP = QENTTMP + QJPENT(K,NJP)
          QUAG = UAG*QJPENT(K,NJP)
          QVAG = VAG*QJPENT(K,NJP)
          QWAG = WAG*QJPENT(K,NJP)
          write(11,1111)K,UJPAVG(K,NJP),VJPAVG(K,NJP),WJPAVG(K,NJP), &
                        QENTTMP,QUAG,QVAG,QWAG,QJPENT(K,NJP)
        enddo
      endif

      ! *** END LOOP OVER ALL JET-PLUME LOCATIONS
    endif   ! *** End of active jet-plume BC flag
  enddo     ! *** End of NQJPIJ loop

  ! *** CLOSE LOG FILE
  if( IFILE == mpi_efdc_out_unit ) close(mpi_efdc_out_unit)
  
  if( LOUTJET )then
    close(10)
  endif
  if( LDEBUG )then
    close(88)
  endif
    
1110 FORMAT(/,'NJP,N = ',I5,I10,/)
1111 FORMAT(I5,50E14.5)
 899 FORMAT(' JPENT ',I5,F12.6,50E12.4)
 898 FORMAT(' FINAL JPENT ',I5,F12.6,50E12.4)
 100 FORMAT(120X)

 134 FORMAT(' BEGIN JET-PLUME NJP,TIME = ',I6,F12.5)
 135 FORMAT(' END JET-PLUME NJP,TIME,KFLAG,KEFFJP,KQJP,QVJET,QVJTOT',' = ',I6,F13.5,3I4,2E12.4)
       
 601 FORMAT(' MAXIMUM ITERATIONS EXCEEDED TIMEDAY N2,NI @ L,HP = ',F10.3,3I6,F8.3,' !!')
 602 FORMAT(' JET-PLUME BOUNDARY PENETRATES SURFACE @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
 603 FORMAT(' JET-PLUME BOUNDARY PENETRATES BOTTOM @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
 604 FORMAT(' JET-PLUME AT NEUTRAL LEVEL @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',A,I5,I8,I5,3F10.3)
 605 FORMAT(' JET-PLUME CENTERLINE PENETRATES SURFACE @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
 606 FORMAT(' JET-PLUME CENTERLINE PENETRATES BOTTOM @ IJ,NE,NI,ZBOT,ZJET,ZSURF = ',I5,I8,I5,3F10.3)
       
 ! *** Write to log file formats
 612 FORMAT(' JET-PLUME BOUNDARY PENETRATES SURFACE @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
 613 FORMAT(' JET-PLUME BOUNDARY PENETRATES BOTTOM @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
 614 FORMAT(' JET-PLUME AT NEUTRAL LEVEL @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',A,F15.5,I5,I8,I5,3F10.3)
 615 FORMAT(' JET-PLUME CENTERLINE PENETRATES SURFACE @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
 616 FORMAT(' JET-PLUME CENTERLINE PENETRATES BOTTOM @ DAY,IJ,NE,NI,ZBOT,ZJET,ZSURF = ',F15.5,I5,I8,I5,3F10.3)
       
 888 FORMAT(A80,/)
 620 FORMAT('NJ,N2,NI,IT,DS,DSO,DF,DFO = ',4I6,6E13.4)
  
  return  
END

FUNCTION FUNDEN(SALIN, SEDIN, TEMIN)

  ! ***  FUNDEN CALCULATED DENSITY AS A FUNCTION OF SAL,TEM,AND SED
  ! CHANGE RECORD
  !  2022-04 Szu-Ting Lee & Paul M. Craig 
  !             Updated funciton to Millero & Huang (2008) to allow for broader salinity and temperature ranges

  implicit none

  real, intent(IN) :: SALIN, SEDIN, TEMIN
  real             :: FUNDEN, TTMP, SSG, SDEN, SSTMP, RHTMP, RHO, RHOCA, RHOCB, RHOCC 

  SSG = 2.5
  SDEN = 1./2500000.

  ! ***  DENSITY AT GIVEN VALUES OF SAL AND TEM
  SSTMP = MIN(MAX(SALIN,0.),70.)
  TTMP  = MIN(MAX(TEMIN,0.),90.)
  RHTMP = 999.842594
  
  RHTMP = RHTMP + 6.793952E-2*TTMP - 9.095290E-3*TTMP*TTMP +1.001685E-4*TTMP*TTMP*TTMP - 1.120083E-6*TTMP*TTMP*TTMP*TTMP + 6.536332E-9*TTMP*TTMP*TTMP*TTMP*TTMP
  
  RHOCA = 0.8174451 - 3.638577E-3*TTMP + 6.480811E-5*TTMP*TTMP - 7.312404E-7*TTMP*TTMP*TTMP +5.330431E-9*TTMP*TTMP*TTMP*TTMP - 1.657628E-11*TTMP*TTMP*TTMP*TTMP*TTMP
  RHOCB = -5.481436E-3 + 3.486075E-5*TTMP - 3.049727E-7*TTMP*TTMP
  RHOCC = 5.346196E-4
  
  RHO = RHTMP + SSTMP*RHOCA + SQRT(SSTMP)*SSTMP*RHOCB + SSTMP*SSTMP*RHOCC
  
  ! ***  CORRECTION FOR SEDIMENT
  RHO = RHO*( (1. - SDEN*SEDIN) + (SSG - 1.)*SDEN*SEDIN )

  ! ***  RETURN DENSITY
  FUNDEN = RHO

END FUNCTION
  
  
SUBROUTINE JPACON(KBOT, ZVAL, UAG, VAG, WAG)

  ! *** JET-PLUME SUB - MODEL SUBROUTINE
  ! *** RETURNS THE AMBIENT CONDITIONS FOR VARIABLE WQ CONDITIONS BASED 

  ! CHANGE RECORD
  
  use GLOBAL
  implicit none   
  
  integer,intent(IN) :: KBOT                                                                                                     
  integer :: K, NT, NS, NX, NZP, MD, NW                                                                        
  real    :: WTNZP, WTNZ, ZVAL                                                                                                
  real    :: UAG, VAG, WAG, DZK

  K = KBOT

  if( ZVAL < ZAD(K) )then
    ! *** Jet/Plume is below the bottom layer
    UAG = UAGD(K)
    VAG = VAGD(K)
    WAG = WAGD(K)
    SALA = SALAD(K)
    TEMA = TEMAD(K)
    do MD = 1,NDYE
      DYEA(MD) = DYEAD(K,MD)
    enddo
    SFLA = SFLAD(K)
    do NT = 1,NTOX
      TOXA(NT) = TOXAD(K,NT)
    enddo
    do NS = 1,NSED2
      SEDA(NS) = SEDAD(K,NS)
    enddo
    do NX = 1,NSND
      SNDA(NX) = SNDAD(K,NX)
    enddo
    do NW = 1,NWQV
      WQVA(NW) = WQVAD(K,NW)
    enddo
    return
  endif

  ! *** GET CONCENTRATIONS BASED ON PLUME ELEVATION (ZVAL)
  if( ZVAL >= ZAD(KC) )then
    ! Jet/plume elevation is above the top layer midpoint
    UAG = UAGD(KC)
    VAG = VAGD(KC)
    WAG = WAGD(KC)
    SALA = SALAD(KC)
    TEMA = TEMAD(KC)
    do MD = 1,NDYE
      DYEA(MD) = DYEAD(KC,MD)
    enddo
    SFLA = SFLAD(KC)
    do NT = 1,NTOX
      TOXA(NT) = TOXAD(KC,NT)
    enddo
    do NS = 1,NSED2
      SEDA(NS) = SEDAD(KC,NS)
    enddo
    do NX = 1,NSND
      SNDA(NX) = SNDAD(KC,NX)
    enddo
    do NW = 1,NWQV
      WQVA(NW) = WQVAD(KC,NW)
    enddo
    return
  endif

  ! *** Determine which layer to extract the ambient conditions from
  
  ! *** Loop over layers
   1000 continue

  ! *** NZP = LAYER ABOVE
  NZP = K + 1
  if( ZVAL >= ZAD(K) .and. ZVAL < ZAD(NZP) )then
    ! *** INTERPOLATE AMBIENT CONDITIONS FROM LAYER INTERVAL
    DZK  = 1./(ZAD(NZP) - ZAD(K))
    WTNZ = DZK*(ZAD(NZP) - ZVAL)
    WTNZP = DZK*(ZVAL - ZAD(K))
    UAG = WTNZ*UAGD(K) + WTNZP*UAGD(NZP)
    VAG = WTNZ*VAGD(K) + WTNZP*VAGD(NZP)
    WAG = WTNZ*WAGD(K) + WTNZP*WAGD(NZP)
    SALA = WTNZ*SALAD(K) + WTNZP*SALAD(NZP)
    TEMA = WTNZ*TEMAD(K) + WTNZP*TEMAD(NZP)
    do MD = 1,NDYE
      DYEA(MD) = WTNZ*DYEAD(K,MD) + WTNZP*DYEAD(NZP,MD)
    enddo
    SFLA = WTNZ*SFLAD(K) + WTNZP*SFLAD(NZP)
    do NT = 1,NTOX
      TOXA(NT) = WTNZ*TOXAD(K,NT) + WTNZP*TOXAD(NZP,NT)
    enddo
    do NS = 1,NSED2
      SEDA(NS) = WTNZ*SEDAD(K,NS) + WTNZP*SEDAD(NZP,NS)
    enddo
    do NX = 1,NSND
      SNDA(NX) = WTNZ*SNDAD(K,NX) + WTNZP*SNDAD(NZP,NX)
    enddo
    do NW = 1,NWQV
      WQVA(NW) = WTNZ*WQVAD(K,NW) + WTNZP*WQVAD(NZP,NW)
    enddo
    return
    
  else
    K = K + 1
    GOTO 1000
  endif
  return
  
END

END MODULE
