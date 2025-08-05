! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details This module creates variables that modified globally by various subroutines
! @date 2015-06  IMPLEMENTED SIGMA-Z (SGZ) IN EE
! @author Paul M. Craig
   
MODULE GLOBAL

  implicit none
  
  !***The following 'dp' variable has 15 significant digits of precision 
  !   and an exponent range of at least 307 (Lemmon and Schafer, 2005, p. 20�22) 
  integer, parameter :: dp = selected_real_kind(15, 307)

  integer(4),parameter :: RKD  = 8    !< REAL(RKD):    THIS ALWAYS DOUBLE PRECISION
  integer(4),parameter :: RK4  = 4    !< REAL(RK4):    THIS ALWAYS REAL*4
  integer(4),parameter :: IK8  = 8    !< INTEGER(IK8): THIS ALWAYS INTEGER*8
  integer(4),parameter :: IKV  = 4    !< INTEGER: DEPRECATED - use COMPILER DATA DEFAULTS TO SPECIFY DEFAULT PRECISION
  integer(4),parameter :: IK4  = 4    !< INTEGER(4): THIS ALWAYS INTEGER*4
  
  integer(4),parameter :: ULOC = 301  !< OPEN UNIT FOR DRIFTER.INP
  integer(4),parameter :: ULGR = 302  !< OPEN UNIT FOR EE_DRIFTER.OUT
  integer(4),parameter :: UCOR = 303  !< OPEN UNIT FOR CORNERS.INP
  integer(4),parameter :: WUNI = 304
  integer(4),parameter :: UGRP = 305
  integer(4),parameter :: UTBL = 306
  integer(4),parameter :: UINP = 307
  integer(4),parameter :: UOUT = 308
  
  character(10) :: EFDC_VER
  character(40) :: EFDC_EXE

  ! *****************************************************************************************************************************
  ! *** EFDC+ Main Controls and Hydrodynamics
  
  character :: projectId*80, title*100, RunID*24
  
  ! *** Global Control Variables
  integer :: ISADAC(0:10)
  integer :: ISFCT(0:10)
  integer :: ISPHXY(10)
  integer :: ISRSPH(10)
  integer :: ISSPH(10)
  integer :: ISTOPT(0:10)    !< Module option flags, see code for options
  integer :: ISTRAN(0:10)    !< Module activation/deactivation flag
  integer :: IWRSP(20)       !< SEDZLJ
  
  ! *** Optimal Number of Threads
  integer :: LDMOPT(0:50)                   !< Number of subdomains for optimal threads
  integer :: NOPTIMAL(0:50)                 !< Optimal number of threads for each hydrodynamic process
  real,allocatable,dimension(:,:) :: TIMES  !< Store the computational time for each number of threads

  logical DEBUG                                         !< Global flag to turn on additional output statements for debugging and model evaluation
  logical,allocatable,dimension(:) :: LMASKDRY          !< LMASKDRY = TRUE CELL IS WET.  CONTROLS ACTIVE CELL OPERATIONS
  logical,allocatable,dimension(:) :: OLDMASK  
  
  character*8  :: OUTDIR
  character*8  :: CARDNO
  character*10 :: RUNTITLE

  integer(IK8) :: N                                     !< Global counter of timesteps based on DT 
  integer(IK8) :: NRESTART                              !< Global counter of timesteps at the time of writing the restart file
  integer(IK8) :: NITER                                 !< Global counter of actual timesteps
  integer(IK8) :: NTS                                   !< Global counter of maximum timesteps based on DT

  integer :: ISTL        !< Time step level: Normal: 2TL - Always 2, 3TL - 3 for standard timestep, 2 for trapezoidal corrector step
  integer :: NTC         !< Number of reference time periods in run
  integer :: NTSPTC      !< Number of time steps per reference time period
  integer :: NTSPTC2     !< 
  integer :: IS2TIM      !< Time stepping option: 0 - 3TL, 1 - 2TL
  integer :: IS2TL       !< Time stepping option: 0 - 3TL, 1 - 2TL
  
  real :: TBEGIN         !< Starting Julian date of the run
  real :: TCON           !< Conversion factor to convert TBEGIN to seconds
  real :: TIDALP         !< Reference time period in seconds (e.g. 3600.0 or 86400.0) 

  integer :: NTSTBC      !< 3tl hydrodynamic step counters
  integer :: NCTBC       !< 3tl hydrodynamic step counters
  
  ! *** Model time counters and intervals
 !REAL(RKD) :: DAYNEXT      !< The next whole day,     i.e.  INT(TIMEDAY) + 1
  real(RKD) :: HOURNEXT     !< The next whole hour,    i.e.  INT(TIMEDAY) + 1/24
  real(RKD) :: HOUR06NEXT   !< The next 6 hour block,  i.e.  INT(TIMEDAY) + 6/24
  real(RKD) :: HOUR12NEXT   !< The next 12 hour block, i.e.  INT(TIMEDAY) + 12/24
  real(RKD) :: TIMEDAY      !< Current simulation time in days
  real(RKD) :: TIMESEC      !< Current simulation time in seconds
  real(RKD) :: TIMEEND      !< Simulation end time in days
  real(RKD) :: TIMERST      !< Restart interval in days
  real(RKD) :: VARYINC      !< Minimum timestep increment for 3TL DT adjustment
  real(RKD) :: VARYMAX      !< Maximum 3TL varying DT
  real(RKD) :: VARYMIN      !< Minimum 3TL varying DT
  real(RKD) :: VARYTHRESH   !< Inflow metric threshold to adjust 3TL DT up or down 

  !< Time and time stepping variables
  real(RKD) :: DT           !< Base delta T (s)
  real(RKD) :: DT2          !< Fpr 3TL DT*2
  real(RKD) :: DTD          !< Delta T (days)
  real(RKD) :: DTDYN        !< Current dynamic timestep (2TL only)
  real(RKD) :: DTDYN1       !< Current dynamic timestep (2TL only) - Previous timesstep
  real(RKD) :: DTI          !< 
  real(RKD) :: DTL1MN       !< Dynamic timestepping: Maximum delta T for method 1
  real(RKD) :: DTL2MN       !< Dynamic timestepping: Maximum delta T for method 2
  real(RKD) :: DTL3MN       !< Dynamic timestepping: Maximum delta T for method 3
  real(RKD) :: DTL4MN       !< Dynamic timestepping: Maximum delta T for method 4
  real(RKD) :: DTMAX        !< 
  real(RKD) :: DTMIN        !< 
  real(RKD) :: DTSSFAC      !< Dynamic timestepping: Safety factor (> 0.0 to enable dynamic timestepping)
  real(RKD) :: DTSSDHDT     !< Dynamic timestepping: Maximum allowable change in HP (method 4)

  real :: DELT              !< Time step
  real :: DELTI             !< 1/DELT

  !< Runtime timing variables
  real(RKD) :: TAVB
  real(RKD) :: TTBXY
  real(RKD) :: TCEXP
  real(RKD) :: TCONG
  real(RKD) :: THDMT
  real(RKD) :: THEAT
  real(RKD) :: THMDF
  real(RKD) :: TLRPD
  real(RKD) :: TMISC
  real(RKD) :: TPUV
  real(RKD) :: TQCTL
  real(RKD) :: TQQQ
  real(RKD) :: TSADV
  real(RKD) :: TTSED
  real(RKD) :: TSSTX
  real(RKD) :: TTWAIT
  real(RKD) :: TUVW
  real(RKD) :: TVDIF
  real(RKD) :: TMPIEE
  real(RKD) :: TMPIGH
  real(RKD) :: TMPITMP
  real(RKD) :: TWQKIN
  real(RKD) :: TWQRPEM
  real(RKD) :: TWQSED
  real(RKD) :: DSITIMING(20)
  
  character*8 :: DSITIME(20)
  
  ! *****************************************************************************************************************************
  ! *** Restart Variables
  character(120) :: RESFILE
  character(120) :: RSTFILE
  character(10)  :: RSTDATE
  character(2)   :: RSTTIME

  integer :: Restart_In_Ver
  
  integer :: ISRESTI     !< Global restart input flag
  integer :: ISRESTIOPT
  integer :: ISRESTO     !< Global restart write flag
  integer :: ISRESTR     !< Global restart
  integer :: ISCI(0:10)      !< RESTART FLAG - INPUT
  integer :: ISCO(0:10)      !< RESTART FLAG - OUTPUT

  integer(IK4)    :: ISGREGOR
  integer(IK4)    :: ICONTINUE
  integer(IK4)    :: RSTFIRST_WC
  integer(IK4)    :: RSTFIRST_WS
  integer(IK4)    :: RSTFIRST_VEL
  integer(IK4)    :: RSTFIRST_WQ
  integer(IK4)    :: RSTFIRST_SEDZLJ
  integer(IK4)    :: RSTFIRST_BED
  integer(IK4)    :: RSTFIRST_SD
  integer(IK4)    :: RSTFIRST_BC
  integer(IK4)    :: RSTFIRST_RPEM
  
  ! *****************************************************************************************************************************
  ! *** OMP
  integer(4)  :: NTHREADS
  real(RKD)   :: TIME_END
  real(RKD)   :: TIME_START
  real(RKD)   :: TCGRS
  
  ! *****************************************************************************************************************************
  ! *** EFDC+ Structured Grid
  integer :: IC          !< Number of columns
  integer :: ICM         !< Maximum number of columns
  integer :: JC          !< Number of rows
  integer :: JCM         !< Maximum number of rows
  integer :: ITRICELL    !< Flag to enable the use of triangular cells

  integer,allocatable,dimension(:,:)   :: IJCT       !< Cell type array
  integer,allocatable,dimension(:,:)   :: IJCTLT     !< Cell type array for wq model interface
  integer,allocatable,dimension(:)     :: IL         !< I index of cell L
  integer,allocatable,dimension(:)     :: JL         !< J index as function of L index
  
  integer,allocatable,dimension(:,:)   :: LADJ       !< List of L indices surrounding the current cell, LADJ(5,L) = L
  integer,allocatable,dimension(:)     :: LCT
  integer,allocatable,dimension(:)     :: LCTLT
  integer,allocatable,dimension(:,:)   :: LIJ        !< L index referenced by I and J
  integer,allocatable,dimension(:,:)   :: LIJLT
     
  integer,target,allocatable,dimension(:) :: LNC     !< Index of the cell to the NORTH of the current L
  integer,allocatable,dimension(:)     :: LNEC       !< Index of the cell to the North-East of the current L
  integer,allocatable,dimension(:)     :: LNWC       !< Index of the cell to the North-West of the current L
  integer,target,allocatable,dimension(:) :: LEC     !< Index of the cell to the EAST of the current L
  integer,allocatable,dimension(:)     :: LSC        !< Index of the cell to the SOUTH of the current L
  integer,allocatable,dimension(:)     :: LSEC       !< Index of the cell to the SOUTH-EAST of the current L
  integer,allocatable,dimension(:)     :: LSWC       !< Index of the cell to the SOUTH-WEST of the current L
  integer,allocatable,dimension(:)     :: LWC        !< Index of the cell to the SOUTH of the current L
  
  real :: DX
  real :: DY
  real :: DZI           
  real(8) :: center_x, center_y                      !< Domain centroid, used to reduce coordinates from double precision to single precision

  real,target,allocatable,dimension(:) :: CUE        !< Cell rotation matrix
  real,target,allocatable,dimension(:) :: CUN        !< Cell rotation matrix
  real,target,allocatable,dimension(:) :: CVE        !< Cell rotation matrix
  real,target,allocatable,dimension(:) :: CVN        !< Cell rotation matrix
  
  integer :: KC          !< Number of water column layers, uniform for Sigma stretch grid, variable for Sigma-Zed grid
  integer :: KCM         !< Maximum KC
  integer :: KS          !< KS = KC-1
  integer :: KSM         !< Maximum KS

  integer :: LA          !< Number of active cells minus one.  L = 2 is first active cell.  Active cells range from 2 to LA
  integer :: LC          !< LC = LA + 1
  integer :: LCM         !< Maximum number of cells, LCM = LA + 2
  integer :: LVC
  
  integer :: NDM         !< Number of OpemMP threads used to divide the domain
  integer :: LDM         !< Number of active cells in each multi-thread region for OpenMP implementation

  integer :: IGRIDH      !< Horizontal grid option (currently unused)
  integer :: KMINV       !< Minimum allowable vertical layers for Sigma-ZED grids
  integer :: IGRIDV      !< Vertical gridding option: 0 - Sigma stretch, > 0 for Sigma-ZED grids
  real    :: SGZHPDELTA  !< Water level offset for Sigma-ZED grid initialization

  integer :: ISMASK      !< Flag to enable masks (thin barriers) at cell faces 
  
  real :: BELADJ         !< Bottom elevation conversion to meters: Offset
  real :: BELCVRT        !< Bottom elevation conversion to meters: Multiplier

  real,target,allocatable,dimension(:) :: BELV       !< Sediment bed-water column interface elevation
  real,allocatable,dimension(:)     :: BELV0         !< Initial condition BELV
  real,allocatable,dimension(:)     :: BELV1         !< BELV for non-hydrostatic solution
  real,allocatable,dimension(:)     :: BELAGW        !< BELV with groundwater offset
  
  real,allocatable,dimension(:)     :: DLAT          !< Cell center lat or northing coordinate
  real,allocatable,dimension(:)     :: DLON          !< Cell center lon or easting coordinate

  ! *** Depths
  real,allocatable,dimension(:)     :: HRU           !< SUB*HMU*DYU*DXIU (m)
  real,allocatable,dimension(:)     :: HRV           !< SVB*HMV*DXV*DYIV (m)
  real,allocatable,dimension(:)     :: H1C           !< Corner or grid vertex value of depth - Previous timestep
  real,allocatable,dimension(:)     :: H1P           !< HP - Previous timestep
  real,allocatable,dimension(:)     :: H1U           !< HU - Previous timestep
  real,allocatable,dimension(:)     :: H1UI          !< 1/H1U
  real,allocatable,dimension(:)     :: H1V           !< HV - Previous timestep
  real,allocatable,dimension(:)     :: H1VI          !< 1/H1V
  real,allocatable,dimension(:)     :: H2P           !< HP - N-2 timestep
  real,allocatable,dimension(:)     :: HMP           !< Still water depth at cell center
  real,allocatable,dimension(:)     :: HMPW          
  real,allocatable,dimension(:)     :: HMU           !< Still water depth at U velocity face
  real,allocatable,dimension(:)     :: HMUW          
  real,allocatable,dimension(:)     :: HMV           !< Still water depth at V velcoity face
  real,allocatable,dimension(:)     :: HMVW
  real,target,allocatable,dimension(:) :: HP         !< Water depth at cell center
  real,allocatable,dimension(:)     :: HPI           !< 1/HP
  real,allocatable,dimension(:)     :: HPTMP         !< Temporary value of HP
  
  real,allocatable,dimension(:)     :: HU            !< Water depth at U velocity face
  real,allocatable,dimension(:)     :: HUI           !< 1/HU
  real,allocatable,dimension(:)     :: HUTMP         !< Temporary value of HU
  real,allocatable,dimension(:)     :: HV            !< Water depth at V velocity face
  real,allocatable,dimension(:)     :: HVI           !< 1/HV
  real,allocatable,dimension(:)     :: HVTMP         !< Temporary value of HV

  ! *** Computational flags/metrics
  real,allocatable,dimension(:)     :: DXIU          !< 1./DXU
  real,allocatable,dimension(:)     :: DXIV          !< 1./DXV
  real,allocatable,dimension(:)     :: DXP           !< Cell dimension in X or I direction at cell center L
  real,allocatable,dimension(:)     :: DXU           !< Cell dimension in X or I direction at U velocity face L
  real,allocatable,dimension(:)     :: DXV           !<  0.5*(DXP+DXP(LS))  (m)
  real,allocatable,dimension(:)     :: DXYIP         !<  1./(STCAP*DXP*DYP)   (1/m2)
  real,allocatable,dimension(:)     :: DXYIU         !<  1./(DXU*DYU)   (1/m2)
  real,allocatable,dimension(:)     :: DXYIV         !<  1./(DXV*DYV)   (1/m2)
  real,allocatable,dimension(:)     :: DXYP          !<  STCAP*DXP*DYP  (m2)
  real,allocatable,dimension(:)     :: DXYU          !<  DXU*DYU        (m2)
  real,allocatable,dimension(:)     :: DXYV          !<  DXV*DYV        (m2)
  
  real,allocatable,dimension(:)     :: DYIU          !< 1./DYU
  real,allocatable,dimension(:)     :: DYIV          !< 1./DYV
  real,allocatable,dimension(:)     :: DYP           !< CELL dimension IN Y OR J DIRECTION AT CELL CENTER L
  real,allocatable,dimension(:)     :: DYU           !< CELL dimension IN Y OR J DIRECTION AT U VELOCITY FACE L
  real,allocatable,dimension(:)     :: DYV           !< CELL dimension IN Y OR J DIRECTION AT V VELOCITY FACE 
  real,allocatable,dimension(:)     :: HRUO          !< SUBO*DYU*DXIU   (dimensionless)
  real,allocatable,dimension(:)     :: HRVO          !< SVBO*DXV*DYIV   (dimensionless)
  real,allocatable,dimension(:)     :: SBX           !<  0.5*SBX*DYU      (m)
  real,allocatable,dimension(:)     :: SBXO          !<  SBX
  real,allocatable,dimension(:)     :: SBY           !<  0.5*SBY*DXV      (m)
  real,allocatable,dimension(:)     :: SBYO          !<  SBY

  real,allocatable,dimension(:)     :: STCAP         !< Area adjustment for diagonal cells
  real,allocatable,dimension(:)     :: SUB           !< Flag to turn on/off cell boundary flows: U face
  real,allocatable,dimension(:)     :: SUBO          !< Flag to turn on/off cell boundary flows: U face - initial condition
  real,allocatable,dimension(:)     :: SVB           !< Flag to turn on/off cell boundary flows: V face
  real,allocatable,dimension(:)     :: SVBO          !< Flag to turn on/off cell boundary flows: V face - initial condition

  integer :: ISCONNECT   !< Flag to enable user specified cell face connections
  integer :: NPEWBP      !< Number of E-W connectors
  integer :: NPNSBP      !< Number of N-S connectors

  integer,allocatable,dimension(:)     :: ISPNS      !< North-South cell connector, I index of South cell
  integer,allocatable,dimension(:)     :: INPNS      !< North-South cell connector, I index of North cell
  integer,allocatable,dimension(:)     :: JSPNS      !< North-South cell connector, J index of South cell
  integer,allocatable,dimension(:)     :: JNPNS      !< North-South cell connector, J index of North cell
  
  integer,allocatable,dimension(:)     :: IEPEW      !< East-West cell connector, I index of East cell 
  integer,allocatable,dimension(:)     :: IWPEW      !< East-West cell connector, I index of West cell 
  integer,allocatable,dimension(:)     :: JEPEW      !< East-West cell connector, J index of East cell     
  integer,allocatable,dimension(:)     :: JWPEW      !< East-West cell connector, J index of West cell    

  real :: QCHERR   
  real,allocatable,dimension(:) :: CHANFRIC      !< Subrgid scale channel friction coefficient dimensionless
  real,allocatable,dimension(:) :: CHANLEN       !< Subgrid scale channel connector length
  real,allocatable,dimension(:) :: PMDCH         !< SURFACE ELEVATION POTENTIAL IN CHANNEL HOST CELL
  real,allocatable,dimension(:) :: QCHANU        !< U DOMINANT SUBGRID SCALE CHANNEL CONNETION FLOW L*L*L/T
  real,allocatable,dimension(:) :: QCHANUN       !< U DOMINANT SUBGRID SCALE CHANNEL CONNETION FLOW L*L*L/T
  real,allocatable,dimension(:) :: QCHANV        !< V DOMINANT SUBGRID SCALE CHANNEL CONNETION FLOW L*L*L/T
  real,allocatable,dimension(:) :: QCHANVN       !< V DOMINANT SUBGRID SCALE CHANNEL CONNETION FLOW L*L*L/T

  ! *** SIGMA-Z SGZ
  integer :: LAWET                                !< NUMBER OF ACTIVE/WET CELLS FOR CURRENT ITERATION
  integer :: LADRY                                !< NUMBER OF FIRST TIME INACTIVE/DRY CELLS FOR CURRENT ITERATION
  integer :: LASED                                !< NUMBER OF ACTIVE/WET SEDIMENT CELLS FOR CURRENT ITERATION
  integer :: LASEDO                               !< NUMBER OF ACTIVE/WET SEDIMENT CELLS FOR CURRENT ITERATION
  integer :: LASGZ1                               !< NUMBER OF ACTIVE/WET CELLS FOR CURRENT ITERATION WHERE KC-KSZ+1 = 1 (EQUIVALENT KC = 1)
  integer :: LDMWET                               !< NUMBER OF ACTIVE/WET CELLS PER DOMAIN
  integer :: LDMDRY                               !< NUMBER OF ACTIVE/WET CELLS PER DOMAIN 
  integer :: LDMSED                               !< NUMBER OF ACTIVE/WET SEDIMENT CELLS PER DOMAIN
  integer :: LDMSEDO                              !< NUMBER OF ACTIVE/WET SEDIMENT CELLS PER DOMAIN
  integer :: LDMSGZ1                              !< NUMBER OF ACTIVE/WET CELLS PER DOMAIN  (EQUIVALENT KC = 1)
  integer :: NHDMF                                !< NUMBER OF ACTIVE/WET HDMF CELLS
  
  integer,allocatable,dimension(:,:)   :: LLWET   !< NUMBER OF ACTIVE CELLS FOR EACH LAYER BY DOMAIN.  ND = 0 FOR ENTIRE LAYER
  integer,allocatable,dimension(:,:,:) :: LKWET   !< L INDEX FOR THE SPECIFIED LAYER AND DOMAIN
  integer,allocatable,dimension(:,:)   :: LLWETZ  !< NUMBER OF ACTIVE CELLS WHERE KSZ/=KC FOR EACH LAYER BY DOMAIN
  integer,allocatable,dimension(:,:,:) :: LKWETZ  !< L INDEX FOR THE SPECIFIED LAYER AND DOMAIN
  integer,allocatable,dimension(:,:)   :: LLHDMF  !< NUMBER OF HMDF ACTIVE CELLS FOR EACH LAYER BY DOMAIN.  ND = 0 FOR ENTIRE LAYER
  integer,allocatable,dimension(:,:,:) :: LKHDMF  !< L INDEX FOR THE HDMF CELLS
  integer,allocatable,dimension(:,:)   :: LLVEG   !< NUMBER OF HMDF ACTIVE CELLS FOR EACH LAYER BY DOMAIN.  ND = 0 FOR ENTIRE LAYER
  integer,allocatable,dimension(:,:,:) :: LKVEG   !< L INDEX FOR THE HDMF CELLS
  
  integer,allocatable,dimension(:)   :: LWET      !< L INDEX FOR THE WET CELLS
  integer,allocatable,dimension(:)   :: LDRY      !< L INDEX FOR THE DRY CELLS.  DRY CELLS ONLY INITIALIZED ON FIRST DRY.
  integer,allocatable,dimension(:)   :: LSED      !< L INDEX FOR THE SEDIMENT CELLS
  integer,allocatable,dimension(:)   :: LSEDO     !< L INDEX FOR THE SEDIMENT CELLS
  integer,allocatable,dimension(:)   :: LSGZ1     !< L INDEX FOR THE WET CELLS WHOSE LAYER COUNT IS 1
  integer,allocatable,dimension(:)   :: KSZ       !< LAYER NUMBER OF BOTTOM MOST ACTIVE LAYER FOR THE SPECIFIED CELL
  integer,allocatable,dimension(:)   :: KSZU      !< LAYER NUMBER OF ACTIVE LAYER MATCHING THE BOTTOM LAYER OF THE WEST CELL
  integer,allocatable,dimension(:)   :: KSZV      !< LAYER NUMBER OF ACTIVE LAYER MATCHING THE BOTTOM LAYER OF THE SOUTH CELL
  
  integer,allocatable,dimension(:)   :: KSZE      !< LAYER NUMBER OF ACTIVE LAYER MATCHING THE BOTTOM LAYER OF THE EAST CELL
  integer,allocatable,dimension(:)   :: KSZW      !< LAYER NUMBER OF ACTIVE LAYER MATCHING THE BOTTOM LAYER OF THE WEST CELL
  integer,allocatable,dimension(:)   :: KSZN      !< LAYER NUMBER OF ACTIVE LAYER MATCHING THE BOTTOM LAYER OF THE NORTH CELL
  integer,allocatable,dimension(:)   :: KSZS      !< LAYER NUMBER OF ACTIVE LAYER MATCHING THE BOTTOM LAYER OF THE SOUTH CELL
  
  logical,allocatable,dimension(:,:) :: LKSZ      !< LOGICAL FLAG IF THE SPECIFIED CELL AND LAYER IS NOT ACTIVE (TRUE = NOT ACTIVE)
  logical,allocatable,dimension(:,:) :: LHDMF     !< LOGICAL FLAG IF THE SPECIFIED HDMF CELL AND LAYER IS SELECTED
  logical,allocatable,dimension(:)   :: LBED      !< LOGICAL FLAG IF THE SPECIFIED VEGETATION CELL IS SELECTED
  logical,allocatable,dimension(:)   :: LVEG      !< LOGICAL FLAG IF THE SPECIFIED VEGETATION CELL IS SELECTED

  real,allocatable,dimension(:)   :: BELVE
  real,allocatable,dimension(:)   :: BELVW
  real,allocatable,dimension(:)   :: BELVN
  real,allocatable,dimension(:)   :: BELVS
  real,allocatable,dimension(:,:) :: BE
  real,allocatable,dimension(:,:) :: BW
  real,allocatable,dimension(:,:) :: BN
  real,allocatable,dimension(:,:) :: BS
  real,allocatable,dimension(:)   :: BEE
  real,allocatable,dimension(:)   :: BEW
  real,allocatable,dimension(:)   :: BEN
  real,allocatable,dimension(:)   :: BES
  real,allocatable,dimension(:)   :: BI1E
  real,allocatable,dimension(:)   :: BI1W
  real,allocatable,dimension(:)   :: BI1N
  real,allocatable,dimension(:)   :: BI1S
  real,allocatable,dimension(:)   :: BI2E
  real,allocatable,dimension(:)   :: BI2W
  real,allocatable,dimension(:)   :: BI2N
  real,allocatable,dimension(:)   :: BI2S
  real,allocatable,dimension(:,:) :: DZIGSD4U
  real,allocatable,dimension(:,:) :: DZIGSD4V
  real,allocatable,dimension(:)   :: HPE
  real,allocatable,dimension(:)   :: HPW
  real,allocatable,dimension(:)   :: HPN
  real,allocatable,dimension(:)   :: HPS
  real,allocatable,dimension(:,:) :: CDZDU
  real,allocatable,dimension(:,:) :: CDZDV
  real,allocatable,dimension(:,:) :: CDZFU
  real,allocatable,dimension(:,:) :: CDZFV
  real,allocatable,dimension(:,:) :: CDZMU
  real,allocatable,dimension(:,:) :: CDZLU
  real,allocatable,dimension(:,:) :: CDZLV
  real,allocatable,dimension(:,:) :: CDZMV
  real,allocatable,dimension(:,:) :: CDZRU
  real,allocatable,dimension(:,:) :: CDZRV
  real,allocatable,dimension(:,:) :: CDZUU
  real,allocatable,dimension(:,:) :: CDZUV
  real,allocatable,dimension(:,:) :: ZE
  real,allocatable,dimension(:,:) :: ZW
  real,allocatable,dimension(:,:) :: ZN
  real,allocatable,dimension(:,:) :: ZS
  real,allocatable,dimension(:,:) :: ZZE
  real,allocatable,dimension(:,:) :: ZZW
  real,allocatable,dimension(:,:) :: ZZN
  real,allocatable,dimension(:,:) :: ZZS

  real,allocatable,dimension(:)   :: DZCK         !< DEFAULT/GLOBAL SIGMA LAYER SPLITS
  real,target,allocatable,dimension(:,:) :: DZC   !< CELL BY CELL SIGMA LAYER SPLITS
  real,allocatable,dimension(:,:) :: DZIC         !< 1./DZC
  real,allocatable,dimension(:,:) :: DZCU
  real,allocatable,dimension(:,:) :: DZCV
  real,allocatable,dimension(:,:) :: DZG          !< VERTICAL LAYER THICKNESS AT LAYER INTERFACE
  real,allocatable,dimension(:,:) :: DZIG         !< 1./DZG
  real,allocatable,dimension(:,:) :: DZGU
  real,allocatable,dimension(:,:) :: DZGV
  real,allocatable,dimension(:,:) :: SUB3D
  real,allocatable,dimension(:,:) :: SUB3DO
  real,allocatable,dimension(:,:) :: SVB3D
  real,allocatable,dimension(:,:) :: SVB3DO
  real,allocatable,dimension(:,:) :: SGZKE
  real,allocatable,dimension(:,:) :: SGZKW
  real,allocatable,dimension(:,:) :: SGZKN
  real,allocatable,dimension(:,:) :: SGZKS
  real,allocatable,dimension(:,:) :: SGZE
  real,allocatable,dimension(:,:) :: SGZW
  real,allocatable,dimension(:,:) :: SGZN
  real,allocatable,dimension(:,:) :: SGZS
  real,allocatable,dimension(:,:) :: SGZU
  real,allocatable,dimension(:,:) :: SGZV
  real,allocatable,dimension(:,:) :: FSGZU
  real,allocatable,dimension(:,:) :: FSGZV
  
  logical,allocatable,dimension(:,:) :: LSGZU
  logical,allocatable,dimension(:,:) :: LSGZV

  real,allocatable,dimension(:,:) :: HPK           !< Water layer thickness (m)
  real,allocatable,dimension(:,:) :: HPKI          !< 1/HPK
  real,allocatable,dimension(:,:) :: H1PK          !< Water layer thickness (m) - Previous timestep
  real,allocatable,dimension(:,:) :: H2PK          !< Water layer thickness (m) - N-2 timestep
  real,allocatable,dimension(:,:) :: UHDY          !< Layer specific U flows = U(L,K)*HU(L)*DYU(L)
  real,allocatable,dimension(:,:) :: VHDX          !< Layer specific V flows = V(L,K)*HV(L)*DXV(L)
  real,allocatable,dimension(:,:) :: UHDY1         !< UHDY - Previous timestep
  real,allocatable,dimension(:,:) :: VHDX1         !< VHDX - Previous timestep
  real,target,allocatable,dimension(:,:) :: UHDY2  !< 3TL uses previous 3TL step flows for W calc in CALUVW. 2TL-3TL computes average discharge between the two most recent time steps,  
  real,target,allocatable,dimension(:,:) :: VHDX2  !< 3TL uses previous 3TL step flows for W calc in CALUVW. 2TL-3TL computes average discharge between the two most recent time steps, 

  real,allocatable,dimension(:,:) :: UHDYF         !< Total discharge on the U face by layer with shear (m3/s)
  real,allocatable,dimension(:,:) :: UHDYF1        !< Total discharge on the U face by layer with shear (m3/s) - Previous timestep
  real,allocatable,dimension(:,:) :: UHDYF2        !< Total discharge on the U face by layer with shear (m3/s) - N-2 timestep
  real,allocatable,dimension(:,:) :: VHDXF         !< Total discharge on the V face by layer with shear (m3/s)
  real,allocatable,dimension(:,:) :: VHDXF1        !< Total discharge on the V face by layer with shear (m3/s) - Previous timestep
  real,allocatable,dimension(:,:) :: VHDXF2        !< Total discharge on the V face by layer with shear (m3/s) - N-2 timestep
 
  ! *** Blocked layer face option (thin barriers similar to masks)
  integer :: ISBLOCKED                             !< Flag if cells are using partially blocked u and/or v faces 
  integer :: NBLOCKED                              !< Number of cells whose u and/or v faces are partially blocked
  integer,allocatable,dimension(:) :: KBBU         !< Bottom active layer accounting for the bottom layer of the west cell
  integer,allocatable,dimension(:) :: KBBV         !< Bottom active layer accounting for the bottom layer of the south cell
  integer,allocatable,dimension(:) :: KTBU         !< Top active layer accounting for the top layer of the west cell
  integer,allocatable,dimension(:) :: KTBV         !< Top active layer accounting for the top layer of the south cell
  integer,allocatable,dimension(:) :: LBLOCKED     !< List of nblocked cells whose u and/or v faces are blocked 
                                                   
  real,allocatable,dimension(:) :: BLANCHORU       !< Elevation of the anchor point for the U draft, 0.0 to allow draftu to float [m]
  real,allocatable,dimension(:) :: BLANCHORV       !< Elevation of the anchor point for the V draft, 0.0 to allow draftv to float [m]
  real,allocatable,dimension(:) :: BLTSILLU        !< Height above bottom to block flow on U face [m]
  real,allocatable,dimension(:) :: BLTSILLV        !< Height above bottom to block flow on V face [m]
  real,allocatable,dimension(:) :: BLDRAFTUO       !< Depth of blocked flow if blanchoru = 0 or fixed height below anchor flow at u interface [m]
  real,allocatable,dimension(:) :: BLDRAFTVO       !< Depth of blocked flow if blanchorv = 0 or fixed height below anchor flow at v interface [m]
  real,allocatable,dimension(:) :: BLSILLU         !< Height from bottom of blocked flow at U interface, e.g. sill  [m]
  real,allocatable,dimension(:) :: BLSILLV         !< Height from bottom of blocked flow at V interface, e.g. sill  [m]
  real,allocatable,dimension(:) :: BLDRAFTU        !< Depth from surface of blocked flow at U interface, e.g. vessel draft  [m]
  real,allocatable,dimension(:) :: BLDRAFTV        !< Depth from surface of blocked flow at V interface, e.g. vessel draft  [m]

  ! *** Channel Modifiers/Pipes
  integer :: ISCHAN      !< Channel modifier flag (i.e. "pipes"), 0 - Do not use, 
  integer :: NCHANM      !< Number of channel modifiers
  integer :: MDCHH
  integer :: MDCHHD
  integer :: MDCHHD2
  
  integer,allocatable,dimension(:) :: IMDCHH       !< I index of channel host cell
  integer,allocatable,dimension(:) :: JMDCHH       !< J index of channel host cell
  integer,allocatable,dimension(:) :: IMDCHU       !< I index of U channel
  integer,allocatable,dimension(:) :: JMDCHU       !< J Index of channel host cell
  integer,allocatable,dimension(:) :: IMDCHV       !< I index of V channel
  integer,allocatable,dimension(:) :: JMDCHV       !< J Index of U channel cell
  integer,allocatable,dimension(:) :: MDCHTYP      !< Channel modifier direction
                                                   
  integer,allocatable,dimension(:) :: LMDCHH       !< List of cells for channel modifiers
  integer,allocatable,dimension(:) :: LMDCHU       !< List of cells for channel modifiers
  integer,allocatable,dimension(:) :: LMDCHV       !< List of cells for channel modifiers
                                                   
  real,allocatable,dimension(:) :: CCCCHH          !< External mode channel interaction coef
  real,allocatable,dimension(:) :: CCCCHU          !< External mode channel interaction coef
  real,allocatable,dimension(:) :: CCCCHV          !< External mode channel interaction coef
  
  ! *****************************************************************************************************************************
  ! *** Waves
  
  type WAVE
    real :: TWX                !<  BED SHEAR STRESS BY WAVE X     (M2/S2)
    real :: TWY                !<  BED SHEAR STRESS BY WAVE Y     (M2/S2)
    real(RKD) :: K             !<  ANGULAR WAVE NUMBER            (RAD/M)
    real(RKD) :: KHP           !<  ANGULAR WAVE NUMBER*DEPTH        (RAD)
    real(RKD) :: UDEL          !<  ORBITAL VELOCITY                 (M/S)
    real(RKD) :: LENGTH        !<  WAVE LENGTH                        (M)
    real(RKD) :: HEIGHT        !<  WAVE HEIGHT, APPLIED               (M)
    real(RKD) :: HEISIG        !<  WAVE HEIGHT, SIGNIFICANT           (M)
    real(RKD) :: DIR           !<  WAVE DIRECTION (RADIANS) COUNTER-CLOCKWISE (CELL-EAST AXIS,WAVE)
    real(RKD) :: FREQ          !<  WAVE FREQUENCY                 (1/SEC) 
    real(RKD) :: PERIOD        !<  WAVE PERIOD                      (SEC)
    real(RKD), allocatable :: DISSIPA(:)    !< WAVE DISSIPATION     (SEC)
  end type  WAVE
  
  logical,allocatable,dimension(:) :: LWVMASK           !< MASK TO DETERMINE IF WAVE CALCUATIONS ARE ON/OFF FOR EACH CELL

  integer :: ISWRSR            !< 
  integer :: ISWVSD            !< 
  integer :: IUSEWVCELLS       !< Spatially varying wave option: 0 - Use all cells, Use varying wave by cell
  integer :: NWVCELLS    !< Number of wave cells

  integer,allocatable,dimension(:,:,:) :: LWDIR      !< List of cells along the fetch for windwave module
  integer,allocatable,dimension(:)     :: LWVCELL    !< L index of wave cells  

  real :: RSWRSI         !< Waves - Irrotational
  real :: RSWRSR         !< Waves - Rotaional
  real :: WVDISH         !< Fraction Of Wave Dissipation As Source In Horiz Smagorinky's Subgrid Closure
  real :: WVDISV         !< Fraction Of Wave Dissipation As Source In Vertical TKE Closure

  real(RKD), parameter :: WLMIN = 0.01
  real(RKD), parameter :: WLMAX = 300
  real(RKD), parameter :: EPS = 0.001
  real(RKD), parameter :: QQMAX = 0.2       !<  MAXIMUM TURBULENT INTENSITY 20% (I = u'/U)
  real(RKD), parameter :: WHMI = 1D-3       !<  MINIMUM WAVE HEIGHT
  real(RKD)            :: SHLIM             !<  MAXIMUM ANGULAR WAVE NUMBER*DEPTH

  real,allocatable,dimension(:)     :: UWVSQ       !< Waves - Square of the orbital velocity
  real,allocatable,dimension(:)     :: WVDTKEM     !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVDTKEP     !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVENEP      !< Wave current interaction variable
  real,allocatable,dimension(:,:)   :: WVHUU       !< Wave current interaction variable: Radiation Stresses SXX  (kg/m^3)(m/s^2)(m^2)
  real,allocatable,dimension(:,:)   :: WVHUV       !< Wave current interaction variable: Radiation Stresses SXY
  real,allocatable,dimension(:,:)   :: WVHVV       !< Wave current interaction variable: Radiation Stresses SYY
  real,allocatable,dimension(:)     :: WVKHC       !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVKHU       !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVKHV       !< Wave current interaction variable
  real,allocatable,dimension(:,:)   :: WVPP        !< Wave current interaction variable
  real,allocatable,dimension(:,:)   :: WVPU        !< Wave
  real,allocatable,dimension(:,:)   :: WVPV        !< Wave
  real,allocatable,dimension(:)     :: WVTMP1      !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVTMP2      !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVTMP3      !< Wave current interaction variable
  real,allocatable,dimension(:)     :: WVTMP4      !< Wave current interaction variable

  type(WAVE), target, allocatable :: WV(:)
  
  ! *** External Waves
  real :: TAUCMIN        !< Minimum bottom shear for wave calculations with radiation diffusion (m2/s2)
  real :: WVPRD          !<  Wave Period (sec)
  
  ! *** SWAN
  type CELL
    integer, allocatable :: ICEL(:),JCEL(:)
    real(RKD),    allocatable :: XCEL(:)                 
    real(RKD),    allocatable :: YCEL(:)  
    character(20),allocatable :: NAME(:)       !< Used to write the location name in HF netCDF output
  end type 
  
  integer(4) :: SWANGRP, ISSTEAD
  integer    :: NLOC, NWVTIM, WVLCAL
  real(RKD), allocatable :: WAVEDAY(:)
  type(CELL) :: SWNLOC

  ! *****************************************************************************************************************************
  ! *** Boundary Conditions - External Forcing
  integer :: NASER       !< Number of atmospheric time series 
  integer :: NASERM      !< Number of atmospheric time series - Maximum
  integer :: NDASER      !< Maximum number of points in a ASER series

  integer :: NWSER       !< Number of wind time series 
  integer :: NWSERM      !< Number of wind time series - Maximum

  integer :: NQSER       !< Number of flow series
  integer :: NQSERM      !< Number of flow series - Maximum
  integer :: NDQSER      !< Maximum number of points in a flow series

  integer :: NDCSER      !< Number of points in a concentration series
  integer :: NCSERM      !< Number of points in a concentration series - Maximum
  integer,allocatable,dimension(:,:) :: MCSER      !< Number points in concentration series
  
  integer :: NQWRSR      !< Number of withdrawal/return series
  integer :: NQWRSRM     !< Number of withdrawal/return series - Maximum
  integer :: NDQWRSR     !< Maximum number of points in a withdrawal/return series
  
  integer :: NGWSER      !< Number of groundwater series
  integer :: NGWSERM     !< Number of groundwater series - Maximum
  integer :: NDGWSER     !< Maximum number of points in a groundwater series
  
  integer :: NPSER       !< Number of pressure/elevation series
  integer :: NPSERM      !< Number of pressure/elevation series - <aximum
  real :: PDGINIT        !< PSER series elevation offset (times G) (m2/s2)
  integer,allocatable,dimension(:)   :: INTPSER    !< Open BC pressure series option

  integer :: NQCTLT      !< Number of lookup tables for hydraulic structures
  integer :: NQCTTM      !< Number of lookup tables for hydraulic structures - Maximum
  integer :: NDQCLT      !< Number of upstream/downstream rating curve lookup pairs
  integer,allocatable,dimension(:)   :: MQCTL      !< Number of points in hydraulic structure

  real :: DAGWZ          !< GW Calc: Groundwater offset to initialize GW elevation (m)
  real :: RIFTRM         !< GW Calc: Maximum infiltration rate
  real :: RNPOR          !< GW Calc: Soil porosity (dimensionless)
  
  real :: TSFSER(100)    !< SFL
  real :: TCSFSER        !< SFL

  real :: RAINCVT        !< ASER rainfall conversion to m/s
  real :: EVAPCVT        !< ASER series - conversion factor for evaporation input units to m/s
  integer,allocatable,dimension(:)  :: IRELH         !< ASER relative humidity/wet bulb temp switch

  real,allocatable,dimension(:)     :: WINDST        !< WSER wind speed at current time (m/s)
  real,allocatable,dimension(:)     :: WINDH         !< Wind anemometer height

  integer,allocatable,dimension(:)  :: IGWSER        !< Groundwater series: Flag to set type of gw flow term.  0-M^3/S,  1-M/S
  integer,allocatable,dimension(:)  :: MGWSER        !< Groundwater series: Number of points in groundwater series
  real,allocatable,dimension(:)     :: TAGWSER       !< Groundwater series: time offset
  real,allocatable,dimension(:)     :: TCGWSER       !< Groundwater series: time multiplier
  real,allocatable,dimension(:,:)   :: TGWSER        !< Groundwater series: Time array
  real,allocatable,dimension(:,:)   :: GWSER         !< Groundwater series:
  real,allocatable,dimension(:)     :: GWSERT        !< Groundwater series:

  real,allocatable,dimension(:,:)   :: TCCSER        !< Constituent concentration series: time multiplier
    
  real,allocatable,dimension(:,:,:) :: GWCSER        !< Ground water concentration series data
  real,allocatable,dimension(:,:)   :: GWCSERT       !< Current groundwater inflow concentration
  real,allocatable,dimension(:)     :: GWFAC         !< Groundwater inflow factor for specific cell
  real,allocatable,dimension(:,:)   :: PSERAVG       !< 
  real,allocatable,dimension(:)     :: PSERT         !< PSER at the current time
  real,allocatable,dimension(:)     :: PSERST        !< PSERS at the current time
  real,allocatable,dimension(:,:,:,:) :: QCTL        !< Hydraulic structure table flows
  real,allocatable,dimension(:,:)   :: QCTLST        !< 
  real,allocatable,dimension(:,:)   :: QCTLSTO       !< 
  real,allocatable,dimension(:,:,:) :: QCTLT         !< Hydraulic structure flow at the current time
  real,allocatable,dimension(:,:)   :: QCTLTLP       !< 
  real,allocatable,dimension(:,:)   :: QCTLTO        !< 
  real,allocatable,dimension(:)     :: AQCTL         !< Hydraulic structure - minimum acceleration depth (m)

  real,allocatable,dimension(:,:)   :: QSERT         !< Current source sink flow
  real,allocatable,dimension(:,:)   :: QSERCELL      !< 
  
  real,allocatable,dimension(:)     :: QWRSERT       !< W/R flow at the current time
  real,allocatable,dimension(:)     :: QWRSERTLP     !< W/R
  real,allocatable,dimension(:,:)   :: CQWRSERT      !< W/R

  integer,allocatable,dimension(:)   :: NCSER        !< Number of concentration series
  integer,allocatable,dimension(:,:) :: NCSERQ       !< Concentration series ID for flow boundary

  real,allocatable,dimension(:)     :: WCOREST       !< 
  real,allocatable,dimension(:)     :: WCORWST       !< 
  real,allocatable,dimension(:)     :: WCORNTH       !< 
  real,allocatable,dimension(:)     :: WCORSTH       !< 
  real,allocatable,dimension(:)     :: WINDCD10      !< 
  real,allocatable,dimension(:)     :: WINDSTKA      !< Wind sheltering
  real,allocatable,dimension(:)     :: WINDSTKA_SAVE !<
  real,target,allocatable,dimension(:) :: WNDVELE    !< Wind speed component: X direction
  real,target,allocatable,dimension(:) :: WNDVELN    !< Wind speed component: Y direction
  
  integer,allocatable,dimension(:)   :: MTSALAST     !< Current index for time series: Atmospheric - ASER
  integer,allocatable,dimension(:,:) :: MTSCLAST     !< Current index for time series: Constituent - SSER, TSER, DSER, SFSER, TXSER, SDSER, SNSER, WQ
  integer,allocatable,dimension(:)   :: MTSGWLAST    !< Current index for time series: Groundwater - GWSER
  integer,allocatable,dimension(:)   :: MTSPLAST     !< Current index for time series: Pressure - PSER
  integer,allocatable,dimension(:)   :: MTSQLAST     !< Current index for time series: Flows - QSER
  integer,allocatable,dimension(:)   :: MTSWLAST     !< Current index for time series: Winds - WSER
  integer,allocatable,dimension(:)   :: MTSWRLAST    !< Current index for time series: Withdrawal/Return - QWR
  integer,allocatable,dimension(:)   :: MTSCUR       !< Time series output counter 
  integer,allocatable,dimension(:)   :: MVEGTLAST    !< Current index for time series
  integer :: MSFTLST     !< SFL
  
  integer,allocatable,dimension(:)     :: MSVDYE     !< List of constituent indices associated with Dye
  integer,allocatable,dimension(:)     :: MSVSED     !< List of constituent indices associated with Cohesive Sediments
  integer,allocatable,dimension(:)     :: MSVSND     !< List of constituent indices associated with Non-cohesive Sediments
  integer,allocatable,dimension(:)     :: MSVTOX     !< List of constituent indices associated with Toxics
  integer,allocatable,dimension(:)     :: MSVWQV     !< List of constituent indices associated with WQ variables
  integer :: NMAXBC, MSVDOX

  integer :: MTIDE       !< Number of tidal constituents to use for open boundary harmonic series
  integer :: MTM         !< Number of tidal constituents to use for open boundary harmonic series - Maximum
  real,allocatable,dimension(:)     :: PSERZDF       !< 
  real,allocatable,dimension(:)     :: PSERZDS       !< 
  real,allocatable,dimension(:,:)   :: PFAM          !< 
  real,allocatable,dimension(:,:)   :: PFPH          !< 
  real,allocatable,dimension(:)     :: TCP           !< Tidal period for each component
                                                     !< 
  integer :: NDPSER      !< Maximum number of points in a pressure series        delme
  type TSR                                           
    integer :: NREC                                  
    real    :: TMULT                                 !< Time series multiplier
    real    :: TOFFSET                               !< Time series offset
    real,allocatable :: VAL(:,:)                     !< Value of time series
    real,allocatable :: TIM(:)                       !< Time array
    real,allocatable :: CONC(:,:)                    !< Constituent value (optional)
  end type 

  ! *** Restructured approach for time series
  integer :: NSTA(8)              
  type(TSR), allocatable :: TSFL(:)                  !< Time series of Flows (QSER) 
  type(TSR), allocatable :: TSPS(:)                  !< Time series of Pressure (PSER) 
  type(TSR), allocatable :: TSWR(:)                  !< Time series of Withdrawal/Return (QWRS)
  type(TSR), allocatable :: TSATM(:)                 !< Time series of Atmosphere (ASER)
  type(TSR), allocatable :: TSWND(:)                 !< Time series of Wind (WSER) 
                                                     
  type(TSR), allocatable :: TSSAL(:)                 !< Time series of Salinity
  type(TSR), allocatable :: TSTEM(:)                 !< Time series of Temperature
  type(TSR), allocatable :: TSDYE(:,:)               !< Time series of Dye
  type(TSR), allocatable :: TSSFL(:)                 !< Time series of Shellfish
  type(TSR), allocatable :: TSICE(:)                 !< Time series of Frazil ice
  type(TSR), allocatable :: TSTOX(:,:)               !< Time series of ChemFate toxics (TXSER)
  type(TSR), allocatable :: TSSED(:,:)               !< Time series of Original-Cohesives, SEDZLJ-All sediments
  type(TSR), allocatable :: TSSND(:,:)               !< Time series of Original-Noncohesives
  type(TSR), allocatable :: TSWQ(:,:)                !< Time series of Eutrophication
  
  ! *** Volatilization Variables
  integer :: IVOLTEMP               !< Flag to indicate source of water temperatures for volatilization calculations: 0 - if ISTRAN(2) > 0, 1 - Constant using TEMO, > 1 - TSER series ID
  real :: VOL_DEP_MIN               !< Volatilization limit: Minimum depth for quiescent conditions    - Depth (m)
  real :: VOL_VEL_MAX               !< Volatilization limit: Maximum velocity for quiescent conditions - Velocity (m/s)
  
  ! *****************************************************************************************************************************
  ! *** Transport
  ! ***  Water column pointer assignment
  real,target,allocatable,dimension(:,:)   :: SAL    !< Salinity                            (PPT) 
  real,target,allocatable,dimension(:,:)   :: SAL1   !< SAL - Previous timestep
  real,target,allocatable,dimension(:,:)   :: TEM    !< Temperature                         (degC) 
  real,target,allocatable,dimension(:,:)   :: TEM1   !< TEM - Previous timestep             
  real,target,allocatable,dimension(:,:,:) :: DYE    !< Dye concentrations, by class        (g/m3)
  real,target,allocatable,dimension(:,:,:) :: DYE1   !< DYE - Previous timestep             
  real,target,allocatable,dimension(:,:)   :: SFL    !< Shellfish                           (g/m3)
  real,target,allocatable,dimension(:,:)   :: SFL2   !< SFL - Previous timestep
  real,target,allocatable,dimension(:,:,:) :: TOX    !< Toxic contaminant concentration     (mg/m3) or (Mass/L**3)
  real,target,allocatable,dimension(:,:,:) :: TOX1   !< TOX - Previous timestep
  real,target,allocatable,dimension(:,:,:) :: SED    !< Cohesive sediment concentration     (g/m3), Also used for all sediments classes when NSEDFLUME > 0
  real,target,allocatable,dimension(:,:,:) :: SED1   !< SED - Previous timestep
  real,target,allocatable,dimension(:,:,:) :: SND    !< Non-Cohesive sediment concentration (g/m3)
  real,target,allocatable,dimension(:,:,:) :: SND1   !< SND - Previous timestep
  real,target,allocatable,dimension(:,:,:) :: WQV    !< Water quality model variable        (g/m3) other than bacteria
  real,target,allocatable,dimension(:,:,:) :: WQVO   !< WQV - Previous timestep

  real ::  DYESTEPW      !< TIME STEP FOR DYE KINETIC UPDATE            (SECONDS)
  real ::  DYESTEPB      !< TIME STEP FOR DYE BED DIFFUSION AND MIXING  (SECONDS)
  
  integer :: ISTRANACTIVE  !< Flag indicating if ANY(ISTRAN(:)) > 0

  type VOLTYPE
    ! *** Volatilization variables
    integer :: KL_OPT = 0           !< Liquid film coefficient option for quiescent (i.e. lake) conditions, 1 - O'Connor, 2 - MacKay & Yeun");
    real :: MW        = 0.0         !< Molecular Weight  (g/mole)
    real :: HE        = 0.0         !< Henry's constant, ( PA(atm)*m^3(water)/mole(gas) )
    real :: TCOEFF    = 0.0         !< Mass Transfer Temperature (dimensionless)
    real :: AIRCON    = 0.0         !< Atmospheric concentration  (ug/l)
    real :: MULT      = 0.0         !< Volatilization rate adjustment factor  (dimensionless)
  end type  VOLTYPE
  
  type DYECLASS
    integer :: ITYPE                !< Dye class type, 0 - conservative, 1 - non-conservative, 2 - water age
    integer :: ICFLAG               !< IC options: 0 - constant (use IC below). 1 - use spatially varying ic from dye.inp 
    integer :: IVOL                 !< Volatilization option: 0 - Do not use, 1 - calculate volatilization losses
    real    :: KRATE0               !< 0th order decay/growth rate at reference temperature (TREF) degc (1/day)
    real    :: KRATE1               !< First order decay/growth rate at reference temperature (TREF) degc (1/day)
    real    :: TADJ                 !< Not used
    real    :: TREF                 !< Reference temperature (degc)
    real    :: SETTLE               !< Settling rate (m/day)
    real    :: IC                   !< Initial condition (mg/l or days)
    real    :: WCLIMIT              !< Maximum water column concentration allowed to bypass transport calculations.  0.0 to disable
    type(VOLTYPE) :: VOL            !< Volatilization parameters for dye class
  end type  DYECLASS

  type TOXCLASS
    real :: BIO_KB   = 0.0          !< Biodegredation rate in the sediment bed   (1/day)
    real :: BIO_KW   = 0.0          !< Biodegredation rate in the sediment bed   (1/day)
    real :: BIO_MXD  = 0.0          !< Maximum depth in sediment bed where biodegredation is active (m)
    real :: BIO_Q10B = 0.0          !< Reference temperature - Bed    (C)
    real :: BIO_Q10W = 0.0          !< Reference temperature - Water  (C)
    real :: BIO_TB   = 0.0          !< Temperature - Bed   (C)
    real :: BIO_TW   = 0.0          !< Temperature - Water (C)
    real :: BLK_KB   = 0.0          !< Bulk decay coefficient - Bed   (1/day)
    real :: BLK_KW   = 0.0          !< Bulk decay coefficient - Water (1/day)
    real :: BLK_MXD  = 0.0          !< Maximum depth in sediment bed where bulk degradation is active (m)
    real :: WCLIMIT  = 0.0          !< Maximum water column concentration allowed to bypass transport calculations.  0.0 to disable
    type(VOLTYPE) :: VOL            !< Volatilization parameters for TOX class
  end type  TOXCLASS
  
  type SEDCLASS
    real    :: WCLIMIT              !< Maximum water column concentration allowed to bypass transport calculations.  0.0 to disable
  end type  SEDCLASS
  
  type TOX_DEPOSITION
    real    :: TXDRY                !< Constant dry deposition flux (mg/m2/day)
    real    :: TXWET                !< Constant wet deposition concentration (mg/m^3)
    integer :: ITXDRY               !< Dry deposition flux usage: 0-only constant, 1-temporally varying, 2-both 
    integer :: ITXWET               !< Wet deposition flux usage: 0-only constant, 1-temporally varying, 2-both
    integer :: ITDRY                !< Last time step index for dry series
    integer :: ITWET                !< Last time step index for wet series
    real    :: TXDRYCUR             !< Current time dry deposition flux (mg/m2/day)
    real    :: TXWETCUR             !< Current time wet deposition concentration (mg/m^3)
    real    :: WCLIMIT              !< Maximum water column concentration allowed to bypass transport calculations.  0.0 to disable
  end type  TOX_DEPOSITION
  
  type(TSR) :: TXDRYSER(1)          !< Time series of dry deposition fluxes  (mg/m2/day)
  type(TSR) :: TXWETSER(1)          !< Time series of wet deposition concentrations (mg/m^3)
  
  type WCVPOINTER
    character     :: ID*8 = '        '
    real, pointer :: VAL0(:,:) => NULL()
    real, pointer :: VAL1(:,:) => NULL()
    real          :: WCLIMIT        !< Maximum water column concentration allowed to bypass transport calculations.  0.0 to disable
  end type 

  integer :: NDYM
  type(DYECLASS),allocatable,dimension(:) :: DYES     !< DYE CLASSES
  
  ! *****************************************************************************************************************************
  ! *** Hydrodynamics
  logical ISCURVATURE

  integer :: IDRYTBP
  integer :: IRVEC       !< Conjugate gradient option
  integer :: ISDRY
  integer :: ISHELTERVARY
  integer :: ISITB       !< Implicit vegetation drag option
  integer :: ISNEGH      !< Flag to check for negative depths
  integer :: ISRLID      !< Rigid lid flag
  integer :: ITER        !< Iteration count for conjugate gradient solution
  integer :: ITERM       !< Maximum iterations allowed for conjugate gradient solution
  integer :: ITERHPM     !< Wet/Dry option for resetting SUB/SVB for each iteration: 0 - Always reset (default), 1 - Reset only when needed
  integer :: NASPECT     !< Number of large aspect ratio cells

  integer :: NDRYSTP     !< 
  
  integer :: IBSC        !< handle buoyancy temp calcs
  real :: BSC            !< Water density flag (multiplier)
  
  character*20,allocatable,dimension(:) :: CLSL
  character* 5,allocatable,dimension(:) :: SYMBOL
  logical,allocatable,dimension(:) :: LASPECT      !< LARGE ASPECT RATIO LIST
  
  integer :: ICALTB      !< Bed Shear calculation option
  integer :: ICK2COR     !<
  integer :: ICSHOW      !<
  integer :: IINTPG      !<
  integer :: IOSWD       !< Flag for wind shear effect on oil spill movement
  integer :: IS2LMC      !< 
  integer :: ISAVBMX     !< 
  integer :: ISAVCOMP    !< 
  integer :: ISSQL       !< 
  integer :: ISBEXP      !< ee linkage bed save flag
  integer :: ISBLFUC     !< 
  integer :: ISBODYF     !< 
  integer :: ISBSDFUF    !< 
  integer :: ISCDMA      !< 
  integer :: ISCFL       !< 
  integer :: ISCFLM      !< 
  integer :: ISCORTBC    !< Corner correction flag for bed shear stress
  integer :: ISCORTBCD   !< 
  integer :: ISCORV      !< 
  integer :: ISDCCA      !< 
  integer :: ISDIQ       !< 
  integer :: ISDIVEX     !< Compute divergence on a cell by cell basis
  integer :: ISDSOLV     !< 
  integer :: ISDYNSTP    !< 
  integer :: ISDZBR      !< waves

  integer :: ISHDMF      !< HMD
  integer :: ISHDMFILTER !< HMD
  integer :: IHMDSUB     !< HMD flag: 0 - Use all cells, 1 - Use spatially subset of cells
  integer :: ISFAVB      !< 
  integer :: ISGWIE      !< 
  integer :: ISGWIT      !< 
  integer :: ISHOW       !< 
  integer :: ISHPRT      !< 
  integer :: ISLOG       !< 
  integer :: ISLTMT      !< Flag to activate decoupled solution (not available)
  integer :: ISLTMTS     !< 
  integer :: ISMMC       !< 
  integer :: ISPD        !< 
  integer :: ISPERC      !< 
  integer :: ISPNHYDS    !< 
  integer :: ISQQ        !< 

  integer :: ISPPH       !< Primary EE linkage flag
  integer :: NPPPH       !< EE Linkage output frequency per reference period
  integer :: ISINWV      !< 1 - Activate CFL diagnostics, 2 - EE_Arrays linkage

  integer :: ISWAVE      !< Waves - 
  integer :: ISWCBL      !< Waves - 
  integer :: ISWRSI      !<  1 Activates Inclusion Of Irrotational Component Of Rad Stress
  integer :: ISWQFLUX    !< Flag to indicate if any WQ / Shellfish modules actived
  
  integer :: IWDRAG      !< Wind drag calculation option
  integer :: IWVCOUNT    !< Waves
  
  integer :: JCSHOW      !< SHOWVAL J index
  integer :: JSWAVE      !< Initialization flag - Waves
  
  integer :: L1LOC       !< 2TL dynamic time steps - limiting location - Type 1
  integer :: L2LOC       !< 2TL dynamic time steps - limiting location - Type 2
  integer :: L3LOC       !< 2TL dynamic time steps - limiting location - Type 3
  integer :: L4LOC       !< 2TL dynamic time steps - limiting location - Type 4
  integer :: LMINSTEP    !< L index of the cell that has the lowest time step when using dynamic time stepping, i.e. the controlling cell

  integer :: NFLTMT ! delme
  integer :: NINCRMT     !< 
  integer :: NITERAT     !< 
  integer :: NLTC        !< Ramp-up
  integer :: NLTS        !< Ramp-up
  integer :: NDYE        !< Number of dye classes
  integer :: NPD         !< Number of drifters
  integer :: NPDM        !< Number of drifters - Maximum
  integer :: NPFORM      !< Open BC harmonics
  integer :: NPFORT      !< Open BC harmonics
  
  integer :: NRAMPUP     !< 2TL number of rampup loops
  integer :: NUPSTEP     !< 2TL dynamic time stepping
  integer :: NSHOWC      !< Runtime reporting to screen
  integer :: NSHOWR      !< Runtime reporting to screen
  integer :: NSHTYPE     !< Runtime reporting to screen
  integer :: NSNAPSHOTS  !< 
  integer :: NSTVM       !< Total number of all WC constituents
  integer :: NSTVM2      !< 2*NSTVM, Used for fast settling of propwash resuspended materials
  integer :: NSVSFP      !< 
  integer :: NTCVB       !< buoyancy rampup
  integer :: NTIMER      !< 
          
  integer :: NTSSTSP     !< 
  integer :: NTSSTSPM    !< 
  integer :: NTSVB       !< 
  integer :: NTSWV       !< Transition - Number Of Time Steps For Gradual Introduction Of Wave Forcing - Rampup
  integer :: NTTC        !< Transition - Non-Linear
  integer :: NTTS        !< Transition - Non-Linear
  integer :: NTXM

  integer :: ISWGS84     !< Horizontal grid coordinate system: 0 = UTM (meters), 1 = Geographic (WGS84) (Lon/Lat) (degrees)
  integer :: ISHOUSATONIC
  
  integer :: ISLLIM
  integer :: IFPROX
  integer :: ISVTURB
  real :: BC_EDGEFACTOR

  ! *** General variables
  integer,allocatable,dimension(:)     :: ICFLMP     !< CFL diagnostics
  
  integer,allocatable,dimension(:)     :: ISCDRY     !< Dry cell flag
  
  
  integer,allocatable,dimension(:)     :: ISSBCP     !< Corner corrections - 2TL Only
  
  integer,allocatable,dimension(:,:)   :: KUPW       !< K index of the upwind direction for 
  integer,allocatable,dimension(:)     :: LBERC      !< List of edge cells for inflow momentum adjustments
  integer,allocatable,dimension(:)     :: LBNRC      !< List of edge cells for inflow momentum adjustments
  integer,allocatable,dimension(:)     :: LBSRC      !< List of edge cells for inflow momentum adjustments
  integer,allocatable,dimension(:)     :: LBWRC      !< List of edge cells for inflow momentum adjustments
  
  integer,allocatable,dimension(:,:)   :: LQSPATH    ! *** List if cells followed from original L to final L.  Only used if HDRYMOVE > 0.0
  integer,allocatable,dimension(:)     :: LQSSAVE    ! *** Point source flow boundary cell list - Original cell list.  Only used if HDRYMOVE > 0.0
  integer,allocatable,dimension(:)     :: LQSSAVE0   ! *** Point source flow boundary cell list - Original cell list.  Only used if HDRYMOVE > 0.0
  
  integer,allocatable,dimension(:,:)   :: LUPU       !< L index of cell upwind (upstream) in the U direction
  integer,allocatable,dimension(:,:)   :: LUPV       !< L index of cell upwind (upstream) in the V direction
  
  integer,allocatable,dimension(:)     :: MVPSL      !< Point source series pointer for point source cell
  integer,allocatable,dimension(:)     :: NATDRY     !< WET/DRY time step counter

  integer,allocatable,dimension(:)     :: NSERWQ     !< Temporary holder of global wq series number
  
  integer,allocatable,dimension(:)     :: NWET       !< Number of iterations since cell was wet
  
  real :: AHD
  real :: AHO
  real :: XYRATIO

  real :: ABMAX
  real :: ABMIN
  real :: ABMX
  real :: ABO
  real :: AVCON
  real :: AVCON1
  real :: AVMAX
  real :: AVMIN
  real :: AVMX
  real :: AVO

  real :: CDRAG1         !< Lower drag coefficient for user defined wind drag relationship
  real :: CDRAG2         !< Upper drag coefficient for user defined wind drag relationship
  real :: CE4SUP         !< Vegetation induced turbulence coefficient
  real :: CE4VEG         !< Vegetation induced turbulence coefficient
  real :: CF             !< Coriolis factor
  real :: CFMAX          !< Coriolis - Maximum
  real :: CK2FCX         !< 2TL curvature correction
  real :: CK2FCY         !< 2TL curvature correction
  real :: CK2UUC         !< 2TL curvature correction
  real :: CK2UUM         !< 2TL curvature correction
  real :: CK2UVC         !< 2TL curvature correction
  real :: CK2UVM         !< 2TL curvature correction
  real :: CK2VVC         !< 2TL curvature correction
  real :: CK2VVM         !< 2TL curvature correction
  
  real :: CTE1           !< Vertical turbulence factors
  real :: CTE2           !< Vertical turbulence factors
  real :: CTE3           !< Vertical turbulence factors
  real :: CTE4           !< Vertical turbulence factors
  real :: CTE5           !< Vertical turbulence factors
  real :: CTURB          !< Vertical turbulence factors
  real :: CTURB2         !< Vertical turbulence factors
  real :: CTURB2B        !< Vertical turbulence factors
  real :: CTURB3         !< Vertical turbulence factors
  
  real :: DMLMIN         !< Vertical turbulence
  real :: DS_LAT         !< Model domain centroid for use in computing clear sky solar radiation
  real :: DS_LONG        !< Model domain centroid for use in computing clear sky solar radiation
  
  real :: ERR            !< Global read error code
  
  real :: FOURDPI        !< 4./PI
  real :: FSCORTBC
  real :: FSWRATF        !< Fraction of solar radiation attenuated fast
  real :: FSOLRADMIN     !< Minimum fraction of solar radiation absorbed in the top layer of the water column (dimensionless)
  
  real :: G              !< Gravity default = 9.81 m/s2
  real :: GI             !< Gravity inverse
  real :: GP             !< Gravity with rampup factor - varies from 0.0 to 1.0*G
  real :: GPO            !< Gravity with buoyancy multiplier
  
  real :: HCVRT          !< Initial condition depth conversion factor
  real :: HDRY           !< If ISDRY > 0, this is the depth to control wetting and drying (m)
  real :: HDRYMOVE       !< If ISDRY > 0, if this depth is > 0.0, then this is the minimum depth allowed for QSER inflows into a cell (m)
  real :: HDRYWAV
  real :: HDRYICE
  real :: HMIN           !< Minimum water depth for initial conditions
  real :: HWET           !< Minimum water depth to allow for QSER withdrawals
  
  real :: OSWDA          !< Wind speed dependent coefficient for oil spill wind shear (s/m)
  real :: OSWDB          !< Fixed coefficient for oil spill wind shear (dimensionless fraction)
  
  real :: PI
  real :: PI2
  
  real :: QQLMIN
  real :: QQMIN
  
  real :: RITB           !< Implicit factor
  real :: RITB1          !< Implicit factor
  real :: RNEW           !< 
  real :: ROLD           !< 
  real :: RP             !< 
  real :: RIQMAX         !< 
  real :: RSQM           !< Square of the maximum conjugate gradient error
  
  real :: S2TL           !< Time integration level real switch
  real :: S3TL           !< Time integration level real switch
  real :: SNLT           !< 
  real :: SWRATNF        !< Solar radiation extinction coefficient - Fast (1/m)
  real :: SWRATNS        !< Solar radiation extinction coefficient - Slow (1/m)
  
  real :: VKC            !< Von Karman's coefficient

  real :: WDRAG1         !< Lower wind speed for user defined wind drag relationship
  real :: WDRAG2         !< Upper wind speed for user defined wind drag relationship
  real :: WVLSH          !< Weight For Depth As The Horiz SSG Eddy Viscosity Length Scale 
  real :: WVLSX          !< Weight For Sqrt(Dxdy) As The Horiz SSG Eddy Viscosity Length Scale
  
  real :: ZBRADJ         !< Bottom roughness - Offset
  real :: ZBRCVRT        !< Bottom roughness - Multiplier
  real :: ZBRWALL        !< 
  real :: ZSSMAX         !< SHOWVAL
  real :: ZSSMIN         !< SHOWVAL
  
  real,allocatable,dimension(:,:)   :: AB            !<  Vertical Molecular Diffusiviy, depth normalized  (m/s)
  real,allocatable,dimension(:,:)   :: ACOEF         !< TEMPORARY CONSOLIDATION VARIABLE
  real,target,allocatable,dimension(:) :: AGWELV     !< Groundwater elevation (m)
  real,allocatable,dimension(:)     :: AGWELV1       !< Groundwater elevation - Previous timestep(m)
  real,allocatable,dimension(:)     :: AGWELV2       !< Groundwater elevation - 2nd previous (m)
  real,allocatable,dimension(:,:)   :: AH            !< Horizontal Turbulent Viscosity, depth normalized (m/s)
  real,allocatable,dimension(:,:)   :: AHC           !< Corner value of AH for 3TL otherwize not used
  real,allocatable,dimension(:)     :: AHDXY         !< AHD * DX * DY
  real,allocatable,dimension(:)     :: AHOXY         !< AHO * DX * DY
  real,allocatable,dimension(:,:,:) :: ALOW          !< Lower diagonal of tridiagonal matrix
  
  real,allocatable,dimension(:)     :: APCG          !< Conjugate gradient          
  real,allocatable,dimension(:,:)   :: AQ            !< Vertical turbulence
                                                     
  real,allocatable,dimension(:,:)   :: AV            !< Vertical Turbulent Viscosity, depth normalized   (m/s)  (m2/s [std SI units] / m = m/s)
  real,allocatable,dimension(:)     :: AVOXY         !< Spatially varying AVO
  real,allocatable,dimension(:)     :: AVBXY         !< Spatially varying ABO
  real,allocatable,dimension(:,:)   :: AVUI          !< Inverse of av at U point
  real,allocatable,dimension(:,:)   :: AVVI          !< Inverse of av at V point

  real,allocatable,dimension(:,:)   :: B             !< Buoyancy (dimensionless)
  real,allocatable,dimension(:,:)   :: B1            !< Buoyancy one time step back
  real,allocatable,dimension(:,:)   :: NNTEM         !< Buoyancy (dimensionless) - temperature 
  real,allocatable,dimension(:,:)   :: NNSAL         !< Buoyancy (dimensionless) - salinity

  real,allocatable,dimension(:,:)   :: CAC           !< Coriolis and curvature parameter l*l/t
  real,allocatable,dimension(:)     :: CCC           !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CCCI          !< External mode linear equation coeff - Inverse
  real,allocatable,dimension(:)     :: CCE           !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CCN           !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CCS           !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CCW           !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CCNHTT        !< Conductive (sensible) heat transfer coefficient
  real,allocatable,dimension(:,:)   :: CDZKK         !< Hydrodynamics - layer metrics
  real,allocatable,dimension(:,:)   :: CDZKKP        !< Hydrodynamics - layer metrics
  real,allocatable,dimension(:,:)   :: CDZKMK        !< Hydrodynamics - layer metrics
  real,allocatable,dimension(:)     :: CC            !< External mode diagonal in calpuv
  real,allocatable,dimension(:)     :: CE            !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CN            !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CS            !< External mode linear equation coeff
  real,allocatable,dimension(:)     :: CW            !< External mode linear equation coeff
  real,allocatable,dimension(:,:)   :: CFLCAC        !< Coriolis-curvature courant number
  real,allocatable,dimension(:,:)   :: CFLUUU        !< U velocity courant number
  real,allocatable,dimension(:,:)   :: CFLVVV        !< V velocity cournat number
  real,allocatable,dimension(:,:)   :: CFLWWW        !< W velocity courant number
  real,allocatable,dimension(:,:)   :: COEFK         !< Temporary consolidation variable
  real,allocatable,dimension(:,:)   :: COEFSK        !< Temporary consolidation variable
  real,allocatable,dimension(:,:)   :: CONGW         !< Concentration of toxics in GW flux
  real,allocatable,dimension(:,:)   :: CONPARW       !< Exponent for concentration dependent sed-toxic partitioning in water column
  real,allocatable,dimension(:,:)   :: CONPARWC      !< Exponent for concentration dependent doc-toxic partitioning in water column
  real,allocatable,dimension(:,:)   :: CONT          
  real,allocatable,dimension(:,:)   :: CQBEDLOADX    !< time series
  real,allocatable,dimension(:,:)   :: CQBEDLOADY    !< time series
  real,allocatable,dimension(:,:,:) :: CSERT         !< Input concentration time series data
  real,allocatable,dimension(:,:)   :: CTURBB1       !< Dimensionless turbulence closure coefficient
  real,allocatable,dimension(:,:)   :: CTURBB2       !< Dimensionless turbulence closure coefficient
  real,allocatable,dimension(:,:)   :: CU1           
  real,allocatable,dimension(:,:)   :: CU2           
  real,allocatable,dimension(:,:,:) :: CUPP          !< Upper diagonal of tridiagonal matrix
  real,allocatable,dimension(:)     :: CUU           
  real,allocatable,dimension(:)     :: CVV           
                                                     
  real,allocatable,dimension(:)     :: DIFTOX        !< Toxic contaminant pore water diffusion coef l*l/t
  real,allocatable,dimension(:)     :: DIFTOXS       !< Diffusion rate between the sediment/water column by toxic
  real,allocatable,dimension(:,:)   :: DIFTOXBW      !< Spatially varying diffusion rate between the sediment/water column by toxic
  real,allocatable,dimension(:,:)   :: DML           !< Ratio of QQL/QQ, Nondimensional Length Scale
  real,allocatable,dimension(:)     :: DPDIFTOX      
  real,allocatable,dimension(:,:)   :: DSTRSE        !< Derivative of effective stress with respect to void ratio (L/T)**2
  real,allocatable,dimension(:,:)   :: DU            !< Temporary delta U array
  real,allocatable,dimension(:,:)   :: DV            !< Temporary delta U array
  real,allocatable,dimension(:)     :: DXDJ          !< Change in dx in J direction L
  real,allocatable,dimension(:)     :: DYDI          !< Change in DY IN I direction L
  real,allocatable,dimension(:,:,:) :: DYEINIT       !< Initial dye concentration
  real,allocatable,dimension(:,:)   :: DZBTR         !< TRANSFORMED BED LAYER THICKNESS L
  real,allocatable,dimension(:,:)   :: DZBTR1        !< DZBTR - Previous timestep
                                                     
  real,allocatable,dimension(:,:)   :: FBBX          !< Bouyancy forcing by layer: X
  real,allocatable,dimension(:,:)   :: FBBY          !< Bouyancy forcing by layer: Y
  real,allocatable,dimension(:,:)   :: FBODYFX       !< 
  real,allocatable,dimension(:,:)   :: FBODYFY       !< 
  real,allocatable,dimension(:,:)   :: FCAX          !< 
  real,allocatable,dimension(:)     :: FCAXE         !< 
  real,allocatable,dimension(:,:)   :: FCAY          !< 
  real,allocatable,dimension(:,:)   :: FCAY1         !< 
  real,allocatable,dimension(:)     :: FCAY1E        !< 
  real,allocatable,dimension(:)     :: FCAYE         !< 
  real,allocatable,dimension(:)     :: FCORC         !< 
                                          
  real,allocatable,dimension(:,:)   :: FMDUX         !< Horizontal momentum diffusion
  real,allocatable,dimension(:,:)   :: FMDUY         !< Horizontal momentum diffusion
  real,allocatable,dimension(:,:)   :: FMDVX         !< Horizontal momentum diffusion
  real,allocatable,dimension(:,:)   :: FMDVY         !< Horizontal momentum diffusion
  
  real,allocatable,dimension(:)     :: FP            !< EXTERNAL MODE
  real,allocatable,dimension(:)     :: FP1           !< EXTERNAL MODE
  real,allocatable,dimension(:,:)   :: FPOCB         !< 
  real,allocatable,dimension(:)     :: FPGYE         !< 
  real,allocatable,dimension(:)     :: FPGXE         !<  
  real,allocatable,dimension(:,:)   :: FPOCBST       !< 
  real,allocatable,dimension(:,:)   :: FPOCWST       !< 
  real,allocatable,dimension(:,:)   :: FPROX         !< Dimensionless
  real,allocatable,dimension(:)     :: FPTMP         !< 
  real,allocatable,dimension(:)     :: FSCORTBCV     !< 
  real,allocatable,dimension(:)     :: FUHDYE        !< 
  real,allocatable,dimension(:,:)   :: FUHU          !<  (m5/s3) Temporary variable in calqq
  real,allocatable,dimension(:,:)   :: FUHV          !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:)     :: FVHDXE  
  
  real,allocatable,dimension(:,:)   :: FVHU          !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:,:)   :: FVHV          !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:,:)   :: FWQQ          !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:,:)   :: FWQQL         !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:,:)   :: FWU           !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:,:)   :: FWV           !< TEMPORARY FLUX VARIABLE
  real,allocatable,dimension(:,:)   :: FX            !< External solution - X component
  real,allocatable,dimension(:,:)   :: FX1           !< External solution - X component - previous timestep
  real,allocatable,dimension(:)     :: FXE           !< TEMPORARY FORCING VARIABLE
  real,allocatable,dimension(:,:)   :: FXVEG         !< TEMPORARY FORCING VARIABLE
  real,allocatable,dimension(:)     :: FXVEGE        
  real,target,allocatable,dimension(:,:) :: FXWAVE   !< TEMPORARY FORCING VARIABLE
  real,allocatable,dimension(:,:)   :: FY            !< External solution - Y component
  real,allocatable,dimension(:,:)   :: FY1           !< External solution - Y component - previous timestep
  real,allocatable,dimension(:)     :: FYE           !< TEMPORARY FORCING VARIABLE
  real,allocatable,dimension(:,:)   :: FYVEG         !< TEMPORARY FORCING VARIABLE
  real,allocatable,dimension(:)     :: FYVEGE        !< TEMPORARY FORCING VARIABLE
  real,target,allocatable,dimension(:,:) :: FYWAVE   !< TEMPORARY FORCING VARIABLE
  
  real,allocatable,dimension(:)     :: HCTLDA        !< FLOW CONTROL TABLE VARIABLE
  real,allocatable,dimension(:)     :: HCTLDM        !< FLOW CONTROL TABLE VARIABLE
  real,allocatable,dimension(:)     :: HCTLUA        !< FLOW CONTROL TABLE VARIABLE
  real,allocatable,dimension(:)     :: HCTLUM        !< FLOW CONTROL TABLE VARIABLE
  real,allocatable,dimension(:)     :: HDFUFX        
  real,allocatable,dimension(:)     :: HDFUFY        
  real,allocatable,dimension(:)     :: HDFUF         
  real,allocatable,dimension(:,:)   :: HDIFCTD       !< FLOW CONTROL TABLE VARIABLE
  real,allocatable,dimension(:,:)   :: HDIFCTL       !< FLOW CONTROL TABLE VARIABLE
  real,allocatable,dimension(:)     :: HGDH

  real,allocatable,dimension(:,:)   :: HYDCN         !< BED HYDRAULIC CONDUCTIVITY L/T
  real,allocatable,dimension(:)     :: P             !< G TIMES WATER SURFACE ELEVATION
  real,allocatable,dimension(:)     :: P1            !< P - Previous timestep
  real,allocatable,dimension(:,:)   :: PCBE          !< Open BC - Tides
  real,allocatable,dimension(:,:)   :: PCBN          !< Open BC - Tides
  real,allocatable,dimension(:,:)   :: PCBS          !< Open BC - Tides
  real,allocatable,dimension(:,:)   :: PCBW          !< Open BC - Tides
  real,allocatable,dimension(:)     :: PDIFTOX       !< 
  
  real,allocatable,dimension(:,:)   :: PNHYDS
  real,allocatable,dimension(:)     :: PPH
  real,allocatable,dimension(:,:)   :: PSBE          !< Open BC - Tides
  real,allocatable,dimension(:,:)   :: PSBN          !< Open BC - Tides
  real,allocatable,dimension(:,:)   :: PSBS          !< Open BC - Tides
  real,allocatable,dimension(:,:)   :: PSBW          !< Open BC - Tides
  real,allocatable,dimension(:)     :: PSHADE        !< Solar radiation shading factor (dimensionless)

  real,allocatable,dimension(:)     :: QCELLCTR      !< HORIZONTAL VELOCITY MAGNITUDE AT CELL CENTER
  real,allocatable,dimension(:,:)   :: QCOEF         !< TEMPORARY CONSOLIDATION VARIABLE
  real,allocatable,dimension(:)     :: QDNEG
  real,allocatable,dimension(:)     :: QDWASTE
  real,target,allocatable,dimension(:) :: QGW
  real,allocatable,dimension(:,:)   :: QJPENT
  real,allocatable,dimension(:)     :: QJPENTT
  real,allocatable,dimension(:)     :: QMORPH
  real,allocatable,dimension(:,:)   :: QQ            !< TURBULENT INTENSITY
  real,allocatable,dimension(:,:)   :: QQ1           !< QQ - Previous timestep
  real,allocatable,dimension(:,:)   :: QQ2           !< QQ - N-2 timestep
  real,allocatable,dimension(:,:)   :: QQSQR         
  real,allocatable,dimension(:,:)   :: QQL           !< TURBULENT INTENSITY X TURBULENT LENGTH SCALE
  real,allocatable,dimension(:,:)   :: QQL1          !< QQL - Previous timestep
  real,allocatable,dimension(:,:)   :: QQL2          
  real,allocatable,dimension(:)     :: QQWC          !<  
  real,allocatable,dimension(:)     :: QQWCR         
  real,allocatable,dimension(:)     :: QQWV1         !<  (m2/s2) Bed Turbulent Intensity Due To Waves Only 
  real,allocatable,dimension(:)     :: QQWV2         !<  (m2/s2) Water Column Turbulent Intensity Due To Waves  
  real,allocatable,dimension(:)     :: QQWV3
  real,allocatable,dimension(:)     :: QRAIN
  real,allocatable,dimension(:,:)   :: QSRTLPN       
  real,allocatable,dimension(:,:)   :: QSRTLPP       
  real,target,allocatable,dimension(:,:) :: QSUM     !< NET INFLOW INTO CELL FROM ALL SOURCES AND SINKS
  real,target,allocatable,dimension(:)   :: QSUME    !< SUM OF QSUM OVER ALL LAYERS
  real,allocatable,dimension(:)     :: QSUM1E        
  real,allocatable,dimension(:)     :: QWATPA        
  
  real,allocatable,dimension(:,:)   :: RADBOT        !< Solar radiation at the bottom of the layer
  real,allocatable,dimension(:,:)   :: RADNET        !< Solar radiation absorbed in the current layer                 
  real,allocatable,dimension(:,:)   :: RADTOP        !< Solar radiation at the top of the layer
  real,allocatable,dimension(:,:)   :: RADKE         !< Water column light extinction coefficient         
  real,allocatable,dimension(:)     :: RCX           !< IMPLICIT DRAG RESISTANCE IN X DIRECTION
  real,allocatable,dimension(:)     :: RCY           !< IMPLICIT DRAG RESISTANCE IN Y DIRECITON
  real,allocatable,dimension(:)     :: RHS           
  real,allocatable,dimension(:)     :: RKTOXP        !< PHOTOLOSIS BASE RATE (NOT IMPLEMENTED)  DELME
  
  real,allocatable,dimension(:)     :: RSSBCE        !< Hard boundary flag for cell centered velocity calculations
  real,allocatable,dimension(:)     :: RSSBCN        !< Hard boundary flag for cell centered velocity calculations
  real,allocatable,dimension(:)     :: RSSBCS        !< Hard boundary flag for cell centered velocity calculations
  real,allocatable,dimension(:)     :: RSSBCW        !< Hard boundary flag for cell centered velocity calculations
  
  real,allocatable,dimension(:)     :: SAAX          !< BC switch for E/W open boundaries (dimensionless)
  real,allocatable,dimension(:)     :: SAAY          !< BC switch for N/S open boundaries (dimensionless)
  real,allocatable,dimension(:,:)   :: SALINIT       !< INITIAL SALINITY
  real,allocatable,dimension(:)     :: SCAX          !< BC switch for E/W open boundaries (dimensionless), Radiation option dependent
  real,allocatable,dimension(:)     :: SCAY          !< BC switch for N/S open boundaries (dimensionless), Radiation option dependent
  real,allocatable,dimension(:)     :: SCB           !< BC switch for kinetic processes 
  real,allocatable,dimension(:)     :: SDX           !< BC switch for Horizontal Eddy Viscosity
  real,allocatable,dimension(:)     :: SDY           !< BC switch for Horizontal Eddy Viscosity
  real,allocatable,dimension(:,:)   :: RHOW          !< Water density (kg/m3)
                                                     
  real,allocatable,dimension(:)     :: SKTOXP        !< BASE SOLAR RADIATION FOR TOXIC PHOTOLYSIS
                                                     
  real,allocatable,dimension(:)     :: SPB           !< BOUNDARY PRESSURE FLAG

  real,allocatable,dimension(:)     :: SSSS          
  real,allocatable,dimension(:)     :: STBX          !< Drag coefficient in the U direction
  real,allocatable,dimension(:)     :: STBXO         !< Drag coefficient in the U direction - INITIAL
  real,allocatable,dimension(:)     :: STBY          !< Drag coefficient in the V direction
  real,allocatable,dimension(:)     :: STBYO         !< Drag coefficient in the V direction - INITIAL
  real,allocatable,dimension(:)     :: STCUV         !< Flag to zero shear for diagonal cells
  real,allocatable,dimension(:)     :: SVKEBACK      !< 
  real,allocatable,dimension(:)     :: SVREVC        !< Full heat balance - surface heat exchange - Evaporative (Latent)    
  real,allocatable,dimension(:)     :: SVRCHC        !< Full heat balance - surface heat exchange - Conductive (Sensible)
  real,allocatable,dimension(:)     :: TBX           !< Bed stress in x direction (l/t)**2
  real,allocatable,dimension(:)     :: TBX1          !< TBX - Previous timestep
  real,allocatable,dimension(:)     :: TBY           !< Bed stress in y direction (l/t)**2
  real,allocatable,dimension(:)     :: TBY1          !< TBY - Previous timestep
  
  real,allocatable,dimension(:)     :: TSX          !< SURFACE SHEAR IN THE U DIRECTION
  real,allocatable,dimension(:)     :: TSX1         !< TSX - Previous timestep
  real,allocatable,dimension(:)     :: TSY          !< SURFACE SHEAR IN THE V DIRECTION
  real,allocatable,dimension(:)     :: TSY1         !< TSY - Previous timestep
  
  real,allocatable,dimension(:,:)   :: TVAR1E       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR1N       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR1S       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR1W       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR2C       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR2E       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR2N       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR2S       !< Temporary variable
  real,allocatable,dimension(:,:)   :: TVAR2W       !< Temporary variable
  real,allocatable,dimension(:)     :: TVAR3C       !< Temporary variable
  real,allocatable,dimension(:)     :: TVAR3E       !< Temporary variable
  real,allocatable,dimension(:)     :: TVAR3N       !< Temporary variable
  real,allocatable,dimension(:)     :: TVAR3S       !< Temporary variable
  real,allocatable,dimension(:)     :: TVAR3W       !< Temporary variable
  real,allocatable,dimension(:)     :: TWATER       !< Temporary variable used for Tecplot (mostly unused)
  
  real,target,allocatable,dimension(:,:) :: U       !< U face layer velocity (m/s)
  real,allocatable,dimension(:,:)   :: U1           !< U face layer velocity - Previous timestep
  real,allocatable,dimension(:,:)   :: U2           !< U face layer velocity - Midpoint or two timesteps back (3TL)
  real,allocatable,dimension(:,:)   :: UCTR         !< Cell center U velocity (m/s)
  real,allocatable,dimension(:,:)   :: UCTR1        !< Cell center U velocity - Previous timestep
  real,allocatable,dimension(:,:)   :: UCTR2        !< Cell center U velocity - Two timesteps back
  real,allocatable,dimension(:)     :: UHDYE        !< Total discharge at the U face (m3/s)
  real,allocatable,dimension(:)     :: UHDY1E       !< Total discharge at the U face (m3/s) - Previous timestep
  real,allocatable,dimension(:)     :: UHDY2E       !< Total discharge at the U face (m3/s) - Midpoint
  real,allocatable,dimension(:,:)   :: UHDYEK       !< Total discharge at the U face split into layers (m3/s)
  real,allocatable,dimension(:,:)   :: UHDY1EK      !< Total discharge at the U face split into layers - Previous timestep
  real,allocatable,dimension(:,:)   :: UHDY2EK      !< Total discharge at the U face split into layers - Midpoint 
  real,allocatable,dimension(:)     :: UHE          !< Depth averaged unit width discharge at the U face (2/s)
  real,allocatable,dimension(:,:)   :: UUU          !< Temporary variable in calqq (m3/s2) 
  real,allocatable,dimension(:)     :: UV           !< Corner velocity of U at the southeast corner (m/s)
  real,allocatable,dimension(:)     :: U1V          !< Corner velocity of U at the southeast corner - Previous timestep
  
  real,target,allocatable,dimension(:,:) :: V       !< V face layer velocity (m/s)
  real,allocatable,dimension(:,:)   :: V1           !< V face layer velocity - Previous timestep
  real,allocatable,dimension(:,:)   :: V2           !< V face layer velocity - Midpoint or two timesteps back (3TL)
  real,allocatable,dimension(:,:)   :: VCTR         !< Cell center V velocity 
  real,allocatable,dimension(:,:)   :: VCTR1        !< Cell center V velocity - Previous timestep
  real,allocatable,dimension(:,:)   :: VCTR2        !< Cell center V velocity - Two timesteps back
  real,allocatable,dimension(:)     :: VHDXE        !< Total discharge at the V face (m3/s)
  real,allocatable,dimension(:)     :: VHDX1E       !< Total discharge at the V face (m3/s) - Previous timestep
  real,allocatable,dimension(:)     :: VHDX2E       !< Total discharge at the V face (m3/s) - Two timesteps back
  real,allocatable,dimension(:,:)   :: VHDXEK       !< Total discharge at the V face split into layers (m3/s)
  real,allocatable,dimension(:,:)   :: VHDX1EK      !< Total discharge at the V face split into layers - Previous timestep
  real,allocatable,dimension(:,:)   :: VHDX2EK      !< Total discharge at the V face split into layers - Midpoint 
  real,allocatable,dimension(:)     :: VHE          !< Depth averaged unit width discharge at the V face (m2/s)
  real,allocatable,dimension(:,:)   :: VVV          !<  (m4/s2) Temporary variable in calqq
  real,allocatable,dimension(:)     :: VU           !< Corner velocity of V at the northwest corner (m/s)
  real,allocatable,dimension(:)     :: V1U          !< Corner velocity of V at the northwest corner - Previous timestep
  
  real,allocatable,dimension(:)     :: VDWASTE      !< Cummulative volume of "wasted" water (m3/s)
  
  real,target,allocatable,dimension(:,:) :: W       !< Vertical velocity at layer top (m/s)
  real,allocatable,dimension(:,:)   :: W1           !< W - Previous timestep
  real,allocatable,dimension(:,:)   :: W2           !< W - Midpoint or two timesteps back (3TL)
  
  real,allocatable,dimension(:)     :: WKQ          !< Temporary variable used for inputting layer values
  real,allocatable,dimension(:,:,:) :: WNDWHT      !< Wind station weigh by cell
  real,allocatable,dimension(:)     :: WNDWHT_TEMP !< Temporary value used in domain remapping

  real,allocatable,dimension(:,:)   :: WTCI        !< Open BC constant concentration depth interpolation
  real,allocatable,dimension(:,:)   :: WWW         !< Temporary array

  real,allocatable,dimension(:,:)   :: Z           !< COORDINATE AT TOP OF LAYER INTERFACE DIMENSIONLESS
  real,allocatable,dimension(:,:)   :: ZBEDC       !< TEMPORARY BED LAYER CELL CENTER VERTICAL COORDINATE L
  real,allocatable,dimension(:,:)   :: ZBEDG       !< TEMPROARY BED LAYER CELL INTERFACE (GRID) COORDINATE .
  real,allocatable,dimension(:)     :: ZBEDGT      !< TEMPORARY BED SURFACE COORDINATE L
  real,allocatable,dimension(:)     :: ZBR         !< LOG LAW ROUGHNESS HEIGHT
  real,allocatable,dimension(:)     :: ZBRE        !< EFFECTIVE VALUE OF ZBR
  real,allocatable,dimension(:)     :: ZSRE        !< Suface roughness computed by COARE 3.6
  real,allocatable,dimension(:)     :: ZBRSED      
  real,allocatable,dimension(:,:)   :: ZELBED      !< BED LAYER ELEVATION L
  real,allocatable,dimension(:,:)   :: ZELBED1     !< ZELBED - Previous timestep
  real,allocatable,dimension(:)     :: ZELBEDA     !< AVERAGE OF ZELBED L
  real,allocatable,dimension(:)     :: ZELBEDA1    !< ZELBEDA - Previous timestep
  real,allocatable,dimension(:)     :: ZEQ         !< SEDIMENT CONCENTRATION REFERENCE LEVEL
  real,allocatable,dimension(:)     :: ZEQD        !< ZEQ/(BOTTOM WATER COLUMN LAYER FRACTIONAL THICKNESS
  real,allocatable,dimension(:)     :: ZEQDI       !< 1/ZEQD
  real,allocatable,dimension(:)     :: ZEQI        !< 1/ZEQ

  real,allocatable,dimension(:,:)   :: ZZ          !< Relative height of the mid-point of the layer (LCM,KCM)
  real,allocatable,dimension(:,:)   :: ZZC         !< Relative height of the mid-point of the layer (KCM,LCM)

  ! *** Vegetation
  integer :: ISVEG       !< Main vegetation activation flag
  integer :: ISVEGL      !< Veg - Laminar flow option
  integer :: MVEGOW      !< Vegetation - Open water class
  integer :: MVEGTYP     !< Vegetation
  integer :: NDVEGSER    !< 
  integer :: NVEGCELLS
  integer :: NVEGSER
  integer :: NVEGSERM
  integer :: NVEGTPM

  real :: UVEGSCL        !< Maximum allowable velocity for vegetation drag
  
  integer,allocatable,dimension(:)     :: MVEGL      !< Vegetation type for cell L
  integer,allocatable,dimension(:)     :: MVEGSER    !< Number of points in each vegetation series: VEGSER.INP
  integer,allocatable,dimension(:)     :: NVEGSERV   !< Vegetation series ID
  
  real,allocatable,dimension(:)     :: ALPVEG        !< Vegetation resistance parameter
  real,allocatable,dimension(:)     :: BDLPSQ        !< Vegetation - Projected plant area per unit volume  (1/m)
  real,allocatable,dimension(:)     :: BPVEG         !< Vegetation - Stem diameter
  real,allocatable,dimension(:)     :: HPVEG         !< HEIGHT OF PLANT IN VEGE RESISTANCE .
  real,allocatable,dimension(:)     :: PVEGZ
  real,allocatable,dimension(:)     :: RDLPSQ        !< Vegetation 
  real,allocatable,dimension(:)     :: SCVEG         

  real,allocatable,dimension(:)     :: TAVEGSER      
  real,allocatable,dimension(:)     :: TCVEGSER      
  real,allocatable,dimension(:,:)   :: TVEGSER      
    
  real,allocatable,dimension(:)     :: VEGK
  real,allocatable,dimension(:,:)   :: VEGSERB
  real,allocatable,dimension(:)     :: VEGSERBT
  real,allocatable,dimension(:,:)   :: VEGSERH
  real,allocatable,dimension(:)     :: VEGSERHT
  real,allocatable,dimension(:,:)   :: VEGSERR
  real,allocatable,dimension(:)     :: VEGSERRT
  
  ! *****************************************************************************************************************************
  ! *** Boundary Conditions
  
  integer :: NBCS        !< Number of boundary cells  
  integer,allocatable,dimension(:)     :: LBCS       !< List of all boundary cells

  ! *** Hydraulic Structures
  integer :: NQCTL       !< Number of hydraulic structures
  integer :: NQCTLM      !< Number of hydraulic structures - Maximum

  integer :: NQJPIJ      !< Number of jet-plume cells
  integer :: NJPSM       !< Number of Jet-Plume cells - Maximum
  integer,allocatable,dimension(:)     :: KEFFJP     !< Jet-plume

  ! *** Withdrawal/Return
  integer :: NQWR        !< Number of withdrawal/return cells
  integer :: NQWRM       !< Number of withdrawal/return cells - Maximum
  real,allocatable,dimension(:,:)   :: CQWR          !< Constant rise (+) or fall (-) for each of the constituents being modeled 
    
  ! *** Open Boundaries
  integer :: NCBE        !< Number of open boundary cells - East
  integer :: NCBN        !< Number of open boundary cells - North
  integer :: NCBS        !< Number of open boundary cells - South 
  integer :: NCBW        !< Number of open boundary cells - West

  integer :: NBCSOP      !< Total number of open boundary cells  
  integer :: NBCSOP2     !< Total number of open boundary cells
  integer :: NBCSOP3     !< Total number of open boundary cells
  
  integer :: NBBEM       !< Maximum number of cells in the East  direction
  integer :: NBBNM       !< Maximum number of cells in the North direction
  integer :: NBBSM       !< Maximum number of cells in the South direction
  integer :: NBBWM       !< Maximum number of cells in the West  direction

  integer :: NPBE        !< Number of open boundary cells: East
  integer :: NPBEM       !< Number of open boundary cells: East (max)
  integer :: NPBN        !< Number of open boundary cells: North
  integer :: NPBNM       !< Number of open boundary cells: North (max)
  integer :: NPBS        !< Number of open boundary cells: South
  integer :: NPBSM       !< Number of open boundary cells: South (max)
  integer :: NPBW        !< Number of open boundary cells: West
  integer :: NPBWM       !< Number of open boundary cells: West (max)

  logical,allocatable,dimension(:)     :: LOPENBCDRY !< FLAG TO DETERMINE IF PSER IS BELOW CELL BOTTOM
  
  integer,allocatable,dimension(:,:)   :: NCSERE     !< Concentration series ID for open boundary: EAST
  integer,allocatable,dimension(:,:)   :: NCSERN     !< Concentration series ID for open boundary: NORTH
  integer,allocatable,dimension(:,:)   :: NCSERS     !< Concentration series ID for open boundary: SOUTH
  integer,allocatable,dimension(:,:)   :: NCSERW     !< Concentration series ID for open boundary: WEST
  
  integer,allocatable,dimension(:)     :: NPSERE     !< Open boundary pressure series ID: EAST
  integer,allocatable,dimension(:)     :: NPSERN     !< Open boundary pressure series ID: NORTH
  integer,allocatable,dimension(:)     :: NPSERS     !< Open boundary pressure series ID: SOUTH
  integer,allocatable,dimension(:)     :: NPSERW     !< Open boundary pressure series ID: WEST
  
  integer,allocatable,dimension(:,:,:) :: NLOE       !< Open Boundary last iteration of flows leaving the domain      
  integer,allocatable,dimension(:,:,:) :: NLON       !< Open Boundary last iteration of flows leaving the domain
  integer,allocatable,dimension(:,:,:) :: NLOS       !< Open Boundary last iteration of flows leaving the domain
  integer,allocatable,dimension(:,:,:) :: NLOW       !< Open Boundary last iteration of flows leaving the domain
  
  integer,allocatable,dimension(:)     :: NPSERE1    !< Second open boundary pressure series (NPFORT > 1) ID: EAST
  integer,allocatable,dimension(:)     :: NPSERN1    !< Second open boundary pressure series (NPFORT > 1) ID: NORTH
  integer,allocatable,dimension(:)     :: NPSERS1    !< Second open boundary pressure series (NPFORT > 1) ID: SOUTH
  integer,allocatable,dimension(:)     :: NPSERW1    !< Second open boundary pressure series (NPFORT > 1) ID: WEST
  
  integer,allocatable,dimension(:)     :: ICBE       !< I index of EAST  concentration open boundary
  integer,allocatable,dimension(:)     :: ICBN       !< I index of NORTH concentration open boundary
  integer,allocatable,dimension(:)     :: ICBS       !< I index of SOUTH concentration open boundary
  integer,allocatable,dimension(:)     :: ICBW       !< I index of WEST  concentration open boundary
  
  integer,allocatable,dimension(:)     :: IPBE       !< I index of EAST  open pressure boundary cell
  integer,allocatable,dimension(:)     :: IPBN       !< I index of NORTH open pressure boundary cell
  integer,allocatable,dimension(:)     :: IPBS       !< I index of SOUTH open pressure boundary cell
  integer,allocatable,dimension(:)     :: IPBW       !< I index of WEST  open pressure boundary cell
  integer,allocatable,dimension(:)     :: ISPBE      !< Open boundary pressure specification method - EAST
  integer,allocatable,dimension(:)     :: ISPBN      !< Open boundary pressure specification method - NORTH
  integer,allocatable,dimension(:)     :: ISPBS      !< Open boundary pressure specification method - SOUTH
  integer,allocatable,dimension(:)     :: ISPBW      !< Open boundary pressure specification method - WEST
  integer,allocatable,dimension(:)     :: JCBE       !< J index of concentration EAST  open boundary
  integer,allocatable,dimension(:)     :: JCBN       !< J index of concentration NORTH open boundary
  integer,allocatable,dimension(:)     :: JCBS       !< J index of concentration SOUTH open boundary
  integer,allocatable,dimension(:)     :: JCBW       !< J index of concentration WEST  open boundary
  
  integer,allocatable,dimension(:)     :: JPBE       !< J Index of East open boundary cell
  integer,allocatable,dimension(:)     :: JPBN       !< J Index of North open boundary cell
  integer,allocatable,dimension(:)     :: JPBS       !< J Index of South open boundary cell
  integer,allocatable,dimension(:)     :: JPBW       !< J Index of West open boundary cell
 
  integer,allocatable,dimension(:)     :: LCBE       !< L index of the concentration open boundary for East
  integer,allocatable,dimension(:)     :: LCBN       !< L index of the concentration open boundary for North
  integer,allocatable,dimension(:)     :: LCBS       !< L index of the concentration open boundary for South
  integer,allocatable,dimension(:)     :: LCBW       !< L index of the concentration open boundary for West
  
  integer,allocatable,dimension(:)     :: LPBE       !< L index of the open boundary for East
  integer,allocatable,dimension(:)     :: LPBN       !< L index of the open boundary for North
  integer,allocatable,dimension(:)     :: LPBS       !< L index of the open boundary for South
  integer,allocatable,dimension(:)     :: LPBW       !< L index of the open boundary for West

  integer,allocatable,dimension(:)     :: LOBCS      !< List of open boundary cells
  integer,allocatable,dimension(:)     :: LOBCS2     !< List of cells adjacent to the main open BC cell
  
  integer,allocatable,dimension(:)     :: NGWSL      !< Groundwater series ID
  
  integer,allocatable,dimension(:)     :: NTSCRE     !< Number of time steps to recover changing from outflow to inflow
  integer,allocatable,dimension(:)     :: NTSCRN     !< Number of time steps to recover changing from outflow to inflow
  integer,allocatable,dimension(:)     :: NTSCRS     !< Number of time steps to recover changing from outflow to inflow
  integer,allocatable,dimension(:)     :: NTSCRW     !< Number of time steps to recover changing from outflow to inflow
  
  real,allocatable,dimension(:,:,:) :: CLOE          !< Last outflowing concentration on EAST open boundary
  real,allocatable,dimension(:,:,:) :: CLON          !< Last outflowing concentration on NORTH open boundary
  real,allocatable,dimension(:,:,:) :: CLOS          !< Last outflowing concentration on SOUTH open boundary
  real,allocatable,dimension(:,:,:) :: CLOW          !< Last outflowing concentration on WEST open boundary
  real,allocatable,dimension(:,:,:) :: CBE           !< EAST  open boundary constant concentration
  real,allocatable,dimension(:,:,:) :: CBN           !< NORTH open boundary constant concentration
  real,allocatable,dimension(:,:,:) :: CBS           !< SOUTH open boundary constant concentration
  real,allocatable,dimension(:,:,:) :: CBW           !< WEST  open boundary constant concentration

  ! *** Hydraulic Structures
  real,   allocatable :: HS_TIMES(:,:)          !< TIME (DAYS) OF WHEN TO TRANSITION TO NEW HYDRAULIC STRUCTURE DEFINITION (BY CELL)
  integer,allocatable :: HS_SERIESQ(:,:)        !< HYDRAULIC STRUCTURE LOOKUP TABLE OR EQUATION NUMBER (MUST MATCH NQCTYP) (BY CELL)

  integer             :: NHYDST                 !< Number of equation defined hydraulic structure groups
  integer,allocatable :: HS_XSTYPE(:)           !< Cross-section type: 1:Circle, 2:Ellip, 3:Rectangle, 4:Parabola, 5:V Notch, 6:Trapezoid
  integer,allocatable :: HS_REVERSE(:)          !< Reverse flow flag:  0-Don't allow reverse flow,  1-Allow reverse flow

  real,   allocatable :: HS_WIDTH(:)            !< Width  (diameter/2*ra of ellip) of hydraulic sturcture       [m]
  real,   allocatable :: HS_HEIGHT(:)           !< Height (2*rb of ellip) of hydraulic sturcture                [m]  
  real,   allocatable :: HS_LENGTH(:)           !< Length of culvert                                            [m]  
  real,   allocatable :: HS_MANN(:)             !< Mannings roughness factor  (metric)
  real,   allocatable :: HS_ANGLE(:)            !< V Notch - angle                                              [deg]
  real,   allocatable :: HS_USELEV(:)           !< Upstream invert elevation                                    [m]
  real,   allocatable :: HS_DSELEV(:)           !< Downstream invert elevation                                  [m]
  real,   allocatable :: HS_COEFF(:,:)          !< Hydraulic coefficient (depends on equation type, i.e. NQCTYP)
    
  ! *** Flow
  integer :: NQSIJ       !< Number of flow boundary cells
  integer :: NQSIJM      !< Number of flow boundary cells - Maximum
  integer :: NGRPID      !< Number of unique flow boundary groups
  real,allocatable,dimension(:,:,:)  :: CQS           !< Concentration associated in volume sources
  real,allocatable,dimension(:,:)    :: QSS           !< Constant flows by layer

  type FLOWBC
    integer :: GRPID                 !< Boundary group number
    integer :: I       = 0           !< Cell - I index
    integer :: J       = 0           !< Cell - J index
    integer :: L       = 0           !< Cell - L index
    integer :: NQSMF   = 0           !< Momentum flux option
    integer :: NQSMUL  = 0           !< Flow BC - Type of specified flow in QSER
    integer :: NQSERQ  = 0           !< Index of the QSER time series
    real    :: QFACTOR = 1.0         !< Flow split to convert QSER flow series to the cell flow (dimensionless) 
    real    :: QSSE    = 1.0         !< Constant flows for the cell
    real    :: QWIDTH  = 1.0         !< Flow width for momentum option 
    real    :: RQSMUL  = 1.0         !< Flow factor multilier to convert input values to m3/s
    integer,allocatable :: NCSERQ(:) !< Index of the constituent time series
    real,allocatable    :: CQSE(:)   !< Concentration associated in volume sources
  end type FLOWBC
  
  type HYDRAULIC_STRUCTURES
    integer :: IQCTLU                !< Upstream cell - I index
    integer :: JQCTLU                !< Upstream cell - j index
    integer :: IQCTLD                !< Downstream cell - I index
    integer :: JQCTLD                !< Downstream cell - j index
    integer :: ISTATE                !< Flag indicating the status of the lock 0 - Closed, 1 - Opening soon-Rising, 2 - Opening soon-Lowering, 3 - Open
    integer :: IUSBC                 !< Location within the list of hydraulic structures of the upstream   navigation lock, 0 for most upstream   lock
    integer :: IDSBC                 !< Location within the list of hydraulic structures of the downstream navigation lock, 0 for most downstream lock
    integer :: ICLOSEDUS             !< Flag to indicate if there is a closed gate upstream   of lock
    integer :: ICLOSEDDS             !< Flag to indicate if there is a closed gate downstream of lock
    integer :: LDR                   !< L index for the most downstream closed gate
    integer :: LUR                   !< L index for the most upstream closed gate
    integer :: NQCTYP                !< Type of control structure
    integer :: NQCTLQ                !< Hydraulic structure table id
    integer :: MMASKS                !< Flag to bypass default masks usage in SETBCS
    integer :: NQCMUL                !< Hydraulic structure multiplier option
    integer :: NQCMINS               !< Minimum number of steps required above lower chord
    integer :: QCTLGRP               !< Number identifier to associate physically based flow groups               [only used if NQCTYP = -2]
    integer :: NWRGRP                !< Withdrawal/Return boundary condition group to provide flows for filling/draining operation
    real :: QCTLMU                   !< Multiplier to split the total QCTL rating table into cell specific flows  [only used if NQCTYP = -2]
    real :: BQCLCE                   !< Lower chord elevation (m)
    real :: HQCTLU = 0.0             !< Offset for upstream head (m)
    real :: HQCTLD = 0.0             !< Offset for upstream head (m)
    real :: HS_FACTOR = 1.0          !< Discharge distribution factor
    real :: TRANSIT                  !< Number of seconds to transition from time (t) to time (t+1) for lock filling/emptying
    real :: RQCMUL = 1.0
    real :: CURHEI = 0.0             !< Current status of the structure height (m)
    real :: CURWID = 0.0             !< Current status of the structure width  (m)
    real :: CURSIL = 0.0             !< Current status of the structure sill height (m)
    real :: CURQ   = 0.0             !< Current computed flow rate through the structure (m3/s)
    real :: QNWR   = 0.0             !< Maximum flows used to fill  the chamber (m3/s)
    real :: ZHDR   = -999.           !< Downstream reference water surface elevation (m) 
    real :: ZHUR   = -999.           !< Upstream  reference water surface elevation (m) 
  end type HYDRAULIC_STRUCTURES
  
  type WITHDRAWAL_RETURN
    integer :: GROUPID              !< Boundary group number
    integer :: IQWRU                !< Upstream cell - I index
    integer :: JQWRU                !< Upstream cell - j index
    integer :: KQWRU                !< Upstream cell - k index
    integer :: IQWRD                !< Downstream cell - I index
    integer :: JQWRD                !< Downstream cell - j index
    integer :: KQWRD                !< Downstream cell - k index
    integer :: NQWRSERQ             !< Number of associated volume withdrawal-return flow and concentration rise time series
    integer :: NQWRMFU              !< Flag to account for withdrawal flow momentum flux
    integer :: NQWRMFD              !< Flag to account for return flow momentum flux
    real :: QWR                     !< Constant volume flow rate from withdrawal to return
    real :: BQWRMFU                 !< Upstream momentum flux width (m)
    real :: BQWRMFD                 !< Downstream momentum flux width (m)
    real :: ANGWRMFD                !< Angle for horizontal for return flow momentum flux
  end type WITHDRAWAL_RETURN
  
  integer :: NUDJPC = 1      !< Iteration update count for Jet-Plume
  integer :: NUDJP           !< Iteration update frequency for Jet-Plume
  type JETPLUME  
    integer :: ICALJP        !< Jet-Plume calculation approach
    integer :: IQJP          !< I index of the jet-plume discharge
    integer :: JQJP          !< J index of the jet-plume discharge
    integer :: KQJP          !< K index of the jet-plume discharge   Note - The ZJET overwrites this value
    integer :: IUPCJP        !< I index of the withdrawal cell (W/R option)
    integer :: JUPCJP        !< J index of the withdrawal cell (W/R option)
    integer :: KUPCJP        !< K index of the withdrawal cell (W/R option)
    integer :: ISENT         !< Shear entrainment option
    integer :: ISTJP         !< Plume grouwth stop option
    integer :: IOUTJP        !< Frequency of diagnostic output
    integer :: ISDJP         !< Plume diagnostics option
    integer :: NJEL          !< Maximum number of elements along jet-plume length
    integer :: NJPMX         !< Maximum number of iterations
    integer :: NPORTJP       !< Number of ports in this jet-plume cell
    integer :: NQSERJP       !< Flow series
    integer :: NQWRSERJP     !< W/R series
    integer :: NZPRJP        !< Number of spatial print/save point in vertical

    real    :: QQCJP         !< Constant flow (m3/s)
    real    :: QWRCJP        !< Constant W/R flow (m3/s)
    real    :: XJETL         !< X Coordinate  (m)
    real    :: YJETL         !< Y Coordinate  (m)
    real    :: CFRD          !< Froude Number adjustment factor
    real    :: ZJET          !< Z Coordinate  (m)
    real    :: PHJET         !< Horizontal angle (rad)
    real    :: THJET         !< Vertical angle   (rad)
    real    :: DJET          !< Diameter of jet (m)
    real    :: DJPER         !< Entrainment error criteria
  
    integer,allocatable :: NCSERJP(:)   !< Concentration series
    real,allocatable    :: CQCJP(:,:)   !< Constant rise/fall for each constituent (Constant Flow)
    real,allocatable    :: CWRCJP(:)    !< Constant rise/fall for each constituent (W/R Flow)
  end type JETPLUME
  
  type(FlowBC),save,allocatable,dimension(:) :: BCFL                      !< Flow boundary cell (Local)
  type(FlowBC),save,allocatable,dimension(:) :: BCFL_GL                   !< Flow boundary cell (Global)

  type(HYDRAULIC_STRUCTURES),save,allocatable,dimension(:) :: HYD_STR     !< Hydraulic Structures (Local)
  type(HYDRAULIC_STRUCTURES),save,allocatable,dimension(:) :: HYD_STR_GL  !< Hydraulic Structures (Global)
  
  type(WITHDRAWAL_RETURN),save,allocatable,dimension(:) :: WITH_RET       !< Withdrawal/return (Local)  
  type(WITHDRAWAL_RETURN),save,allocatable,dimension(:) :: WITH_RET_GL    !< Withdrawal/return (Global)
  
  type(JETPLUME),save,allocatable,dimension(:) :: JET_PLM                 !< Jet-Plume (Local)
  type(JETPLUME),save,allocatable,dimension(:) :: JET_PLM_GL              !< Jet-Plume (Global)
  
  ! *****************************************************************************************************************************
  !*** Temperature/Heat
  logical :: LFRAZIL     !< FLAG IS TRUE OF FRAZIL ICE EXISTS
  logical :: LCHECKICE   !< FLAG TO DETERMINE IF ICE CONDITIONS MAY EXIST

  integer :: IASWRAD     !< Solar ratiation adsorbtion option
  integer :: IATMP       !< 
  integer :: ISVHEAT     !< Flag to use Variable surface heat coefficients
  
  real :: HTBED1         !< Bed heat - convective heat transfer coefficient (no dim)
  real :: HTBED2         !< Bed heat - heat transfer coefficient (m/s)

  real :: RAINTSUM1      !< Sum of previous rainfall
  real :: RAINTSUM0      !< Sum of current rainfall
  real :: RCHC           !< 1000*Convective heat transfer coefficient (sensible heat)
  real :: REVC           !< 1000*Evaprative heat transfer coefficient (latent heat)

  real :: TEMO           !< 
  real :: TEMBO          !< 
  real :: TEMTHKO        !< 

  logical COMPUTESOLRAD
  logical USESHADE
  logical LDAYLIGHT                                     !< Flag to indicate if there solar radiation > 0 anywhere in the domain
  logical LRAIN                                         !< Flag to determine if rainfall is available

  logical,allocatable,dimension(:) :: LSVHTWINDE        !< FLAG TO DETERMINE WIND HEAT EXCHANGE COEFFICIENTS ARE BEING USED - EVAPORATION
  logical,allocatable,dimension(:) :: LSVHTWINDC        !< FLAG TO DETERMINE WIND HEAT EXCHANGE COEFFICIENTS ARE BEING USED - CONDUCTION 

  ! *** Atmospheric Conditions
  real,allocatable,dimension(:)        :: PATMT   !< Currant value of ATM pressure (mb)
  real,allocatable,dimension(:)        :: TATMT   !< Currant value of dry bulb temperature (C)
  real,allocatable,dimension(:)        :: TDEWT   !< Currant value of dew point temperature (C)
  real,target,allocatable,dimension(:) :: RAINT   !< Currant value of rainfall (m/s)
  real,allocatable,dimension(:)        :: SOLSWRT !< Currant value of solar radiation (W/m2)
  real,allocatable,dimension(:)        :: CLOUDT  !< Currant value of cloud cover 
  real,allocatable,dimension(:)        :: RHAT    !< Currant value of relative humidity
  real,target,allocatable,dimension(:) :: EVAPT   !< Currant value of evaporation (m/s)
  real,allocatable,dimension(:)        :: VPAT    !< Currant value of vapor pressure
  real,allocatable,dimension(:)        :: SVPAT   !< Currant value of saturated vapor pressure

  integer                              :: IEVAP    !< Evaporation option used for water balance
  real                                 :: WINDFA   !< WIND FACTOR COEFFICIENT A
  real                                 :: WINDFB   !< WIND FACTOR COEFFICIENT B
  real                                 :: WINDFC   !< WIND FACTOR COEFFICIENT C
  real,allocatable,dimension(:)        :: CLEVAP   !< EVAPORATIVE TRANSFER COEFFICIENT DIMENSIONLESS
  real,target,allocatable,dimension(:) :: EVAPGW   !< CURRENT EVAPORATION TAKEN FROM GROUNDWATER LAYER L/T
  real,target,allocatable,dimension(:) :: EVAPSW   !< CURRENT EVAPORATION TAKEN FRON TOP LAYER OF WATER COLUMN L/T
  real,allocatable,dimension(:)        :: EVACOARE !< CURRENT TOTAL EVAPORATION L/T
  real,allocatable,dimension(:)        :: SVPW
                                       
  

  real,allocatable,dimension(:)     :: ATMP          !< Atmospheric pressure at each cell
  real,allocatable,dimension(:,:,:) :: ATMWHT        !< ASER series weighting for each cell
  real,allocatable,dimension(:)     :: ATMWHT_TEMP   !< ASER series temporary variable

  real,allocatable,dimension(:)     :: TEMB          !< Bed temperature deg c
  real,allocatable,dimension(:)     :: TEMB1         !< TEMB - previous timestep
  real,allocatable,dimension(:,:)   :: TEMINIT       !< Initial value of temperature

  ! *** Ice
  integer :: ISICE
  integer :: ISICECOV    !< ISICE = 1  ICE FUNCTION (NOT USED)
  integer :: NISER       !< ISICE = 1  NUMBER OF ICE INPUT SERIES
  integer :: ICETHKFUN   !< ISICE = 3  ICE FUNCTION (NOT USED)
  integer(4) :: ISRHEVAP !< FLAG TO use RYAN HARLEMAN WIND FUNCTION
  
  real ::  AFW           !< WIND FUNCTION A   AFW + BFW*WIND2M**CFW (WATER)
  real ::  BFW           !< WIND FUNCTION B
  real ::  CFW           !< WIND FUNCTION C
  real ::  AFWI          !< WIND FUNCTION A   AFW + BFW*WIND2M**CFW (ICE)
  real ::  BFWI          !< WIND FUNCTION B
  real ::  CFWI          !< WIND FUNCTION C
  real ::  DYEBEG
  real ::  DYEBEG2T
  real ::  DYICEBEG      !< ISICE = 2  
  real ::  DYICEEND      !< ISICE = 2  
  real ::  DYICEM1       !< ISICE = 2  
  real ::  DYICEM2       !< ISICE = 2  
  real ::  RICETHK0      !< ISICE = 2  ICE THICKNESS
  real ::  TEMPICE       !< ISICE = 1  TEMPERTURE OF THE ICE
  real ::  ICETHMX       !< MAXIMUM ICE THICKNESS (M)
  real ::  MINICETHICK   !< ICE GROWTH MINIMUM TO ALLOW ICE COVER FORMATION (m)
  real ::  HWI           !< COEFFICIENT OF WATER-ICE HEAT EXCHANGE, W/(M^2 C)
  real ::  ICEK          !< ICE CONDUCTIVITY, ~2.12 W/M/DEGC
  real ::  ALBEDOI       !< RATIO OF REFLECTION TO INCIDENT RADIATION (ICE)
  real ::  BETAI         !< FRACTION OF SOLAR RADIATION ABSORBED IN THE ICE SURFACE
  real ::  GAMMAI        !< SOLAR RADIATION EXTINCTION COEFFICIENT, 1/m
  real ::  CDICE         !< DRAG COEFFICIENT FOR ICE
  real ::  RHOI          !< DENSITY OF ICE (KG/M3)  (DEFAULT = 916.0)
  real ::  RISEVEL       !< FRAZIL ICE RISE VELOCITY (M/S), TYPICAL VALUES RANGE FROM .003 TO .022 M/S
  real ::  MELTFACTOR    !< ICE MELT FACTOR.  < 1.0 TO REDUCE THE SPEED OF ICE MELT

  real, parameter :: REFL = 0.06    !< SHORTWAVE RADIATION RELFECTED OFF SURFACE  
  

  integer,allocatable :: MITLAST(:)
  
  real,   allocatable :: FRAZILICE(:,:)
  real,   allocatable :: FRAZILICE1(:,:)
  real,   allocatable :: ICECOVER(:)
  real,   allocatable :: ICERATE(:)
  real,target,   allocatable :: ICETHICK(:)
  real,   allocatable :: ICETHICK1(:)
  real,target,   allocatable :: ICETEMP(:)
  real,   allocatable :: ICEVOL(:)
  real,   allocatable :: RICEWHT(:,:,:)
  real,   allocatable :: RICECOVT(:)  
  real,   allocatable :: RICETHKT(:)  
  
  logical(IK4), allocatable :: ICECELL(:)
  
  ! *****************************************************************************************************************************
  ! *** Shellfish Larvae
  integer :: ISSFLDN     !< 
  integer :: ISSFLFE     !< 
  integer :: MSFSER      !< 

  real :: DSFLMNT        !< 
  real :: DSFLMXT        !< 
  real :: RKDSFLT        !< 
  real :: SFATBTT        !< 
  real :: SFLKILL        !< 
  real :: SFNTBET        !< 
  real :: TASFSER        !< 
  real :: TSRSF          !< 
  real :: TSSSF          !< 
  real :: WSFLSMT        !< 
  real :: WSFLSTT        !< 

  real :: DSFLMN(100)    !< 
  real :: DSFLMX(100)    !< 
  real :: RKDSFL(100)    !< 
  real :: SFATBT(100)    !< 
  real :: SFNTBE(100)    !< 
  real :: WSFLSM(100)    !< 
  real :: WSFLST(100)    !< 
  
  real,allocatable,dimension(:,:)   :: ACCWFLD       !< SFL
  real,allocatable,dimension(:,:)   :: SFLINIT       
  real,allocatable,dimension(:)     :: SFLSBOT      
  
  
  ! *****************************************************************************************************************************
  ! *** Chemical fate and transport (ChemFate)  
  integer :: NTOX        !< Number of ChemFate/Toxic constituents
  integer :: ISDTXBUG    !< Debug print flag
  integer :: ISTDOCB     !< ChemFate
  integer :: ISTDOCW     !< ChemFate
  integer :: ISTPOCB     !< ChemFate
  integer :: ISTPOCW     !< ChemFate
  integer :: NPMXPTSM    !< ChemFate - Particle Mixing
  integer :: NPMXZM      !< ChemFate - Particle Mixing
  
  real :: STDOCBC        !< ChemFate
  real :: STDOCWC        !< ChemFate
  real :: STPOCBC        !< ChemFate
  real :: STPOCWC        !< ChemFate

  real ::  TOXSTEPW      !< TIME STEP FOR TOXIC KINETIC UPDATE            (SECONDS)
  real ::  TOXSTEPB      !< TIME STEP FOR TOXIC BED DIFFUSION AND MIXING  (SECONDS)
  
  integer,allocatable,dimension(:)     :: ISTOC      !< Switch to activate total organic carbon/toxics effects
  integer,allocatable,dimension(:)     :: ITOXBU     !< Toxic initialization switch
  integer,allocatable,dimension(:)     :: ISDIFBW    !< ChemFate
  integer,allocatable,dimension(:,:)   :: ITOXKIN    !< ChemFate loss options
  integer,allocatable,dimension(:)     :: ITOXWU     !< Toxic initialization switch
  integer,allocatable,dimension(:)     :: ITXBDUT    !< Toxic initialization switch
  integer,allocatable,dimension(:)     :: ITXINT     !< Toxic initialization switch
  integer,allocatable,dimension(:,:)   :: ITXPARW    !< Water column sediment-toxics partitioning option flag
  integer,allocatable,dimension(:,:)   :: ITXPARWC   !< Water column oc-toxics partitioning option flag
  integer,allocatable,dimension(:)     :: NSP2       !< ChemFate - Total number of sediment classes plus 2 (for POC and DOC)

  integer,allocatable,dimension(:)     :: ISPMXZ
  integer,allocatable,dimension(:)     :: LPMXZ

  real,allocatable,dimension(:,:)   :: PARTMIXZ
  real,allocatable,dimension(:,:)   :: PFPOCB        !< 
  real,allocatable,dimension(:,:)   :: PMXDEPTH
  real,allocatable,dimension(:,:)   :: PMXCOEF
  
  real,allocatable,dimension(:,:,:) :: RRHS          !< ChemFate - RIGHT HAND SIDE OF TRIDIAGONAL LINEAR SYSTEM
  real,allocatable,dimension(:,:)   :: STDOCB        
  real,allocatable,dimension(:,:,:) :: STFPOCB       
  real,allocatable,dimension(:,:,:) :: STFPOCW       
  real,allocatable,dimension(:,:)   :: STPOCB        
  real,allocatable,dimension(:,:)   :: STPOCW        
  real,allocatable,dimension(:,:)   :: STDOCW        
  real,target,allocatable,dimension(:,:,:) :: TOXB   !< Bed TOX by class and layer (mg/m2)
  real,allocatable,dimension(:,:)   :: TADFLUX       !< Total flux by toxic class from/to the sediment bed (mg/m2/s)
  real,allocatable,dimension(:,:,:) :: TOXB1         !< Bed TOX - Previous timestep
  real,allocatable,dimension(:,:,:) :: TOXBINIT      !< Input values of TOXB
  real,allocatable,dimension(:,:,:) :: TOXCDFB
  real,allocatable,dimension(:,:,:) :: TOXCDFW
  real,allocatable,dimension(:,:,:) :: TOXF         !< Bed toxic contaminant flux M/L*L*T
  real,allocatable,dimension(:,:)   :: TOXFB        !< Erosional flux of toxics from the bed to the water column
  real,allocatable,dimension(:,:)   :: TOXFBL
  real,allocatable,dimension(:,:,:) :: TOXFDFB
  real,allocatable,dimension(:,:,:) :: TOXFDFW
  real,allocatable,dimension(:,:,:) :: TOXINIT      !< Low pass filter of TOX
  real,allocatable,dimension(:)     :: TOXINTB      !< Initial toxic contamiant concentration in bed            mg/m2
  real,allocatable,dimension(:)     :: TOXINTW      !< Initial toxic contamiant concentration in water column   mg/m3
  real,allocatable,dimension(:,:,:) :: TOXPARB      !< SED-TOX partition coefficient in bed  L/mg
  real,allocatable,dimension(:,:)   :: TOXPARBC     !< OC-TOX partition coefficient in bed   L/mg
  real,allocatable,dimension(:,:,:) :: TOXPARW      !< SED-TOXIC partition coefficient in water column  L/mg
  real,allocatable,dimension(:,:)   :: TOXPARWC     !< OC-TOXIC Partition coefficient in water column   L/mg
  real,allocatable,dimension(:,:,:,:) :: TOXPFB     !< Toxic particulate phase fraction in bed                     dimensionless
  real,allocatable,dimension(:,:,:)   :: TOXPFTB    !< Total toxic particulate phase fraction in bed               dimensionless
  real,allocatable,dimension(:,:,:,:) :: TOXPFW     !< Toxic particulate phase fraction in water column            dimensionless
  real,allocatable,dimension(:,:,:)   :: TOXPFTW    !< Total toxic particulate phase fraction in water column      dimensionless
  real,allocatable,dimension(:,:,:) :: TOXTMP       !< Temporary variable

  real(RKD),target,allocatable,dimension(:,:) :: CBLTOX       !(LCM,NTXM) 
  real(RKD),allocatable,dimension(:,:,:) :: CBLTXCON          !(LCM,NSEDS,NTXM) 

  ! *** Toxic volatilization variables
  type(VOLTYPE),allocatable,dimension(:) :: TOXVOL

  type(TOX_DEPOSITION), allocatable :: TOXDEP(:)
  type(TOXCLASS),allocatable,dimension(:) :: TOXS    !< All toxic/ChemFate classes
                                                              
  
  ! *****************************************************************************************************************************
  ! *** Sediment Transport (SedTran)
  integer :: IALSTUP     !< During initialization, split the top sediment layer into armor and parent layer 
  integer :: IALTYP      !< Armor layer option: 0 - Constant thickness, 1 - Constant mass
  integer :: IBMECH      !< Bed consolidation option
  integer :: IMORPH      !< Morphologic feedback to hydrodynamics: 0 - No feedback, 1 - Use feedback
  integer :: ICALC_BL    !< Bedload flag:  0 - No bedload, 1 - Bedload
  integer :: ISBEDMAP    !< Flag for erosion/deposition calculations:  0 - All cells, 1 - Only cells listed is BEDMAP.INP
  integer :: ISBEDSTR    !< Bed shear stress otion
  integer :: ISEDBED     !<
  integer :: ISEDBINT    !<
  integer :: ISEDEFF     !< Cohesive hiding option
  integer :: ISEDINT     !<
  integer :: ISEDVW      !< sediment settling vel option
  integer :: ISMUD       !< Flag for cohesive viscous effects for SedTran
  integer :: ISNDAL      !< Activate armoring
  integer :: ISNDVW      !< 
  integer :: ISBSDIAM    !< Sediment diameter option used for bed shear
  
  integer :: KB          !< Maximum number of sediment bed layers.  KB at top for original sedtran module and at the bottom for SEDZLJ
  integer :: KBB         !< BOTTOM LAYER (ORIGINAL:KBB = 1, SEDZLJ:KBB = KB)
  integer :: KBH         !< NOT USED, BUT MAINTAINED FOR GVC
  integer :: KBM         !< Maximum KB

  integer :: NSBDLDBC    !< SedTran - Number of cells to allow bedload to exit model
  integer :: NSCM        !< Number of cohesive sediment classes
  integer :: NSED        !< Number of cohesive sediment classes
  integer :: NSED2       !< 2*NSED, Used for fast settling of propwash resuspended materials
  integer :: NSEDS       !< Total number of sediment classes, NSEDS = NSED + NSND
  integer :: NSEDS2      !< 2*NSCM, Used for fast settling of propwash resuspended materials
  integer :: NSND        !< Number of non-cohesive sediment classes
  integer :: NSND2       !< 2*NSND, Used for fast settling of propwash resuspended materials
  integer :: NSNM        !< Number of non-cohesive sediment classes
  integer :: NSNM2       !< 2*NSND, Used for fast settling of propwash resuspended materials
  
  integer :: NSTM        !< Number of cohesive + non-cohesive sediment classes
  integer :: NSTM2       !< 2*NSTM, Used for fast settling of propwash resuspended materials

  integer,allocatable,dimension(:)     :: IBLTAUC    !< SedTran - Bedload
  integer,allocatable,dimension(:)     :: IROUSE     !< SedTran
  integer,allocatable,dimension(:)     :: ISDBLDIR   !< SedTran
  integer,allocatable,dimension(:)     :: ISBDLD     !< Bed load transport option switch

  integer,allocatable,dimension(:)     :: ISEDBU     !< SedTran
  integer,allocatable,dimension(:)     :: ISEDSCOR   !< SedTran
  integer,allocatable,dimension(:)     :: ISEDWU     !< SedTran
  integer,allocatable,dimension(:)     :: ISLTAUC    !< SedTran
  integer,allocatable,dimension(:)     :: ISPROBDEP  !< SedTran
  integer,allocatable,dimension(:)     :: ISNDBU     !< SedTran
  integer,allocatable,dimension(:)     :: ISNDEQ     !< Equilibrium sediment concentration option
  integer,allocatable,dimension(:)     :: ISNDM1     !< SedTran    
  integer,allocatable,dimension(:)     :: ISNDM2     !< SedTran     
  integer,allocatable,dimension(:)     :: IWRSPB     !< SedTran   
  
  integer,target,allocatable,dimension(:) :: KBT     !< Current K index of top layer of sediment bed
  integer,allocatable,dimension(:)     :: LSBLBCD    !< SedTran - Bedload exit
  integer,allocatable,dimension(:)     :: LSBLBCU    !< SedTran - Bedload exit
  
  real(RKD) :: DTSED

  real :: BEDPORC
  real :: BLBSNT         !< SedTran - Bedload adverse slope
  real :: BMECH1
  real :: BMECH2
  real :: BMECH3
  real :: BMECH4
  real :: BMECH5
  real :: BMECH6
  
  real :: COEFTSBL       !< SedTran
  real :: COEHEFF        !< SedTran
  real :: COEHEFF2       !< SedTran
  real :: DIASED         !< Temporary sediment diameter
  real :: DSEDGMM        !< 1/(sediment density in grams)
  real :: GPDIASED       !< Gravity*(SSG-1)*(sediment diameter)
  
  real :: HBEDAL         !< Sedtran - sediment bed layer thickness
  real :: HBEDALMN       !< Sedtran - sediment bed layer thickness
  real :: HBEDMAX        !< Sedtran - sediment bed layer thickness
  real :: HBEDMIN        !< Minimum bed layer thickness

  real :: SEDMDGM        !< Geometric mean
  real :: SEDMDMN        !< 
  real :: SEDMDMX        !< 
  real :: SEDSTART       !< 
  real :: SEDSTEP        !< 
  real :: SEDVDRD        !< 
  real :: SEDVDRM        !< 
  real :: SEDVRDT        !< 
  real :: SNDDMX         !< 
  real :: SNDVDRD        !< 

  real :: VISMUD         !< 
  real :: VISMUDST       !< 

  real,target,allocatable,dimension(:,:) :: BDENBED  !< Sediment bed bulk density m/l*l*l
  real,allocatable,dimension(:,:)   :: BDENBED1      !< Bdenbed one time step back
  real,allocatable,dimension(:)     :: BDENBEDA      !< Average of bdenbed of depth of bed
  real,allocatable,dimension(:)     :: BDENBEDA1     !< Bdenbeda one time step back
  real,allocatable,dimension(:,:)   :: BEDBINIT      !< Temporary bed initialization array - BULK DENSITY
  real,allocatable,dimension(:,:)   :: BEDDINIT      !< Temporary bed initialization array - VOID RATIO
  real,allocatable,dimension(:,:)   :: BEDLINIT      !< Temporary bed initialization array - LAYER THICKNESS
  real,allocatable,dimension(:)     :: CBEDTOTAL     !< delme - move to ssedtox
  real,allocatable,dimension(:)     :: COSEDHID      
  real,allocatable,dimension(:)     :: CSHIELDS50    !< Shields parameter based on nonchohesive d50
  real,allocatable,dimension(:)     :: CTAUC         

  real,allocatable,dimension(:)     :: FACBEDL       !< Fraction of sediment transported as bed load
  real,allocatable,dimension(:)     :: FACSUSL       !< Fraction of sediment transported as suspended load
  real,allocatable,dimension(:,:)   :: FRACCOH       
  real,allocatable,dimension(:,:)   :: FRACNON       

  real,target,allocatable,dimension(:,:) :: HBED     !< Sediment bed layer thickness
  real,allocatable,dimension(:,:)   :: HBED1         !< HBED - Previous timestep
  real,allocatable,dimension(:)     :: HBEDA         !< Total bed thickness (m)
  real,allocatable,dimension(:)     :: HBEDA1        !< Total bed thicknes - previous timestep

  real,allocatable,dimension(:,:)   :: PEXP          !< Noncohesive sediment armoring exposure function
  real,allocatable,dimension(:,:)   :: PHID          !< Noncohesive sediment armoring hidding function
  real,target,allocatable,dimension(:,:) :: PORBED   !< Porosity of sediment bed
  real,allocatable,dimension(:,:)   :: PORBED1       !< Probed - previous timestep

  real,allocatable,dimension(:,:)   :: QSBDLDIN
  real,allocatable,dimension(:,:)   :: QSBDLDOT
  real,allocatable,dimension(:)     :: QSBDLDP       !< Cell center bed load transport rate l*l/t
  real,target,allocatable,dimension(:,:) :: QSBDLDX  !< Bed load volume flux in X direction l*l*l/t
  real,target,allocatable,dimension(:,:) :: QSBDLDY  !< Bed load volume flux in Y direction l*l*l/t
  real,allocatable,dimension(:)     :: QSBDTOP       !< Volume of sediment exchange at bed surface due to erosion/deposition (m/s)
  real,allocatable,dimension(:,:,:) :: QSEDBED       
  real,allocatable,dimension(:,:,:) :: QSEDBED1      !< Source sink time series data
  real,allocatable,dimension(:,:)   :: QSEDBEDA
  real,allocatable,dimension(:,:)   :: QSEDBEDA1
  real,allocatable,dimension(:)     :: QSSDPA        
  real,allocatable,dimension(:)     :: QWBDTOP       !< Volume of water exchange due to deposition or erosion (Porewater)     (m/s)
  real,allocatable,dimension(:,:)   :: QWTRBED       !< Volume of water exchange (specific discharge) at bed layer interfaces (m/s)

  real,allocatable,dimension(:)     :: RBPSBL        !< Bed load transport mask
  real,allocatable,dimension(:)     :: ROUSE         !<
  real,allocatable,dimension(:)     :: RSNDM         !<
                                                    
  real,allocatable,dimension(:)     :: SBDLDA        !<
  real,allocatable,dimension(:)     :: SBDLDB        !<
  real,allocatable,dimension(:)     :: SBDLDG1       !<
  real,allocatable,dimension(:)     :: SBDLDG2       !<
  real,allocatable,dimension(:)     :: SBDLDG3       !<
  real,allocatable,dimension(:)     :: SBDLDG4       !<
  real,allocatable,dimension(:)     :: SBDLDP        !<
  real,allocatable,dimension(:)     :: SDEN          !< 1./sediment density
  real,allocatable,dimension(:,:)   :: SDENAVG       !< 
  
  real,allocatable,dimension(:)     :: SEXP          !< Sediment processes variable
  real,allocatable,dimension(:,:)   :: SGSM1         !< Specific gravity of sediment minus one dimensionless
  real,allocatable,dimension(:,:)   :: SIGPHI        
  real,allocatable,dimension(:)     :: SIGPHIA       
  real,allocatable,dimension(:,:)   :: STRSE         !< Effective stress in sediment bed (l/t)**2
  real,allocatable,dimension(:)     :: TCSHIELDS     !< Shield's coefficient
  real,allocatable,dimension(:)     :: TEXP          !< Sediment processes coefficient

  
  ! *** COHESIVE CLASS SEDIMENT VARIABLES            
  real,allocatable,dimension(:)     :: SED3DMAX      !< Cohesive sediment concentration scale factor
  real,allocatable,dimension(:)     :: SED3DMIN      !< Cohesive sediment concentration scale factor
  real,target,allocatable,dimension(:,:,:) :: SEDB   !< Sedb - Previous timestep
  real,allocatable,dimension(:,:,:) :: SEDB1         !< Initial cohesive sediment bed concentration, sedb UNITS
  real,allocatable,dimension(:,:)   :: SEDBA         !< Total sediment in bed
  real,allocatable,dimension(:,:,:) :: SEDBINIT      !< Low pass filter sediment bed concentration
  real,allocatable,dimension(:)     :: SEDBO         !< Initial constant bed sediment concentration
  real,allocatable,dimension(:,:)   :: SEDBT         !< Total sediment in each bed layer
  real,allocatable,dimension(:)     :: SEDDIA        !< Sediment diameter l
  real,allocatable,dimension(:,:)   :: SEDDIA50      !< D50 0f noncohesive sediment l
  real,allocatable,dimension(:,:)   :: SEDDIA90      !< 
  real,allocatable,dimension(:,:,:) :: SEDF          !< Vertical mass flux rate for cohesive class sedimenTS
  real,allocatable,dimension(:,:)   :: SEDFPA        !< 
  real,allocatable,dimension(:,:,:) :: SEDINIT       !< 
  real,allocatable,dimension(:)     :: SEDN          !< Normalizing sediment concentration
  real,allocatable,dimension(:)     :: SEDO          !< Initial constant water column sediment concentratiON
  real,allocatable,dimension(:,:,:) :: SED2          !< Temporary cohesive sediment concentration water coLUMN  M/L*L*L
  real,allocatable,dimension(:,:)   :: SEDT          !< Total cohesive sediment concentration in water colUMN
                                                     
  ! *** NON-COHESIVE CLASS SEDIMENT VARIABLES        
  real,allocatable,dimension(:,:,:) :: SNDB          !< SND - Previous timestep
  real,allocatable,dimension(:,:,:) :: SNDB1         !< Initial nocohesiver sediment concentration in bed M/L*L
  real,allocatable,dimension(:,:)   :: SNDBA         !<
  real,allocatable,dimension(:,:,:) :: SNDBINIT      !< Temporary variable used for bed initialization
  real,allocatable,dimension(:,:)   :: SNDBT         !< Total noncohesive sediment bed concentration m/l*l
  real,allocatable,dimension(:)     :: SNDEQ         !< 
  real,allocatable,dimension(:,:)   :: SNDEQSAV      !<
  real,allocatable,dimension(:)     :: SNDEQB        !< Temporary variable
  real,allocatable,dimension(:,:,:) :: SNDF          !< Vertical mass flux rate for non-cohesive class sedIMENTS
  real,allocatable,dimension(:,:)   :: SNDFBL        !< Non-cohesive sediment flux from bed due to bed loaD TRANSPORT  M/L*L*L
  real,allocatable,dimension(:,:)   :: SNDFPA        !< 
  real,allocatable,dimension(:,:,:) :: SNDINIT       !< 
  real,allocatable,dimension(:,:,:) :: SNDS          !< Temporary noncohesive sediment concentration water COLUMN  M/L*L*L
  real,allocatable,dimension(:,:)   :: SNDT          !< Total water column noncohesive sediment concentratION M/L*L*L
  real,allocatable,dimension(:)     :: SSG           !< Sediment specific gravity
  
  real,target,allocatable,dimension(:) :: TAUB       !< Bed shear stress
  real,target,allocatable,dimension(:) :: TAUBSED    
  real,target,allocatable,dimension(:) :: TAUBSND    
  real,allocatable,dimension(:,:)   :: TAUCRCOH      
  real,allocatable,dimension(:)     :: TAUD          !< Critical stress for sediment deposition
  real,allocatable,dimension(:)     :: TAUN          !< Normalizing stress for sediment processes
  real,allocatable,dimension(:)     :: TAUR          !< Critical stress for sediment resuspension
  real,allocatable,dimension(:,:)   :: TAURB         !< Shear strength of sediment in bed layers for bulk RESUSPENSION
  real,allocatable,dimension(:)     :: TAUDS         
  real,allocatable,dimension(:,:)   :: TAUNS         
  real,allocatable,dimension(:,:)   :: TAURS         !< Shear strength of sediment in bed layers for surfaCE RESUSPENSION
  real,allocatable,dimension(:,:)   :: TEXPS         
  real,allocatable,dimension(:)     :: USTAR         !< Shear velocity
  real,allocatable,dimension(:)     :: USTARSED      
  real,allocatable,dimension(:)     :: USTARSND      
  real,allocatable,dimension(:)     :: UCELLCTR      !< Cell center u velocity
  real,allocatable,dimension(:)     :: VCELLCTR      
                                                     
  real,allocatable,dimension(:,:)   :: VDRBED0       !< Initial conditions of vdrbed
  real,allocatable,dimension(:,:)   :: VDRBED        !< Void ratio of the bed
  real,allocatable,dimension(:,:)   :: VDRBED1       !< Vdrbed - Previous timestep
  real,allocatable,dimension(:)     :: VDRBEDA       
  real,allocatable,dimension(:,:)   :: VDRBEDSED     
  real,allocatable,dimension(:,:)   :: VDRBEDSND     
  real,allocatable,dimension(:)     :: VDRDEPO       !< Void ratio - Constant for cohesives and constant for non-cohesives
  real,allocatable,dimension(:)     :: VDRRSPO       !< Void ratio - Input by class
  real,target,allocatable,dimension(:,:,:) :: VFRBED !< Volume fraction of the bed
  real,allocatable,dimension(:,:,:) :: VFRBED1       !< VFRBED - Previous timestep

  real,allocatable,dimension(:,:)   :: WRSPB         !< 
  real,allocatable,dimension(:)     :: WRSPBA        !< 
  real,allocatable,dimension(:)     :: WRSPO         !< Base cohesive sediment resuspension rate
  real,allocatable,dimension(:,:)   :: WRSPS         !< Cohesive sediment resuspension rate with depth in bed
  real,allocatable,dimension(:)     :: WS            !< Settling velocity
  real,allocatable,dimension(:)     :: WSEDO         !< Sediment settling velocity class (m/s)
  real,allocatable,dimension(:,:,:) :: WSETA         !< Temporary cohesive and noncohesive sediment settling velocity
  
  ! *** Hard bottom boundary list
  type BEDEDGE
    integer :: NEDGE
    integer, allocatable,dimension(:) :: LEDGE
  end type 
  
  type(BEDEDGE) :: BEDEDGEE 
  type(BEDEDGE) :: BEDEDGEW    
  type(BEDEDGE) :: BEDEDGEN    
  type(BEDEDGE) :: BEDEDGES    
  
  ! *** Bank erosion
  integer :: ISBKERO     !< Bank erosion activation flag
  integer :: NBEPAIR,NBEPAIRM
  integer :: NBESER,NBESERM
  integer :: NDBESER
  integer,allocatable,dimension(:)     :: IBANKBE
  integer,allocatable,dimension(:)     :: JBANKBE
  integer,allocatable,dimension(:)     :: ICHANBE
  integer,allocatable,dimension(:)     :: JCHANBE
  integer,allocatable,dimension(:)     :: NBESERN
  integer,allocatable,dimension(:)     :: MBESER
  integer,allocatable,dimension(:)     :: MBETLAST

  real,allocatable,dimension(:)     :: FBESER
  real,allocatable,dimension(:)     :: BESERT
  real,allocatable,dimension(:)     :: FWCBESERT
  real,allocatable,dimension(:)     :: TCBESER
  real,allocatable,dimension(:)     :: TABESER
  real,allocatable,dimension(:,:)   :: TBESER
  real,allocatable,dimension(:,:)   :: BESER
  real,allocatable,dimension(:,:)   :: FWCBESER
  real,allocatable,dimension(:)     :: QSBDTOPBEBKB
  real,allocatable,dimension(:)     :: QSBDTOPBECHB
  real,allocatable,dimension(:)     :: QSBDTOPBECHW
  real,allocatable,dimension(:)     :: QWBDTOPBEBKB
  real,allocatable,dimension(:)     :: QWBDTOPBECHB
  real,allocatable,dimension(:)     :: QWBDTOPBECHW
  real,allocatable,dimension(:,:)   :: SEDFBEBKB
  real,allocatable,dimension(:,:)   :: SEDFBECHB
  real,allocatable,dimension(:,:)   :: SEDFBECHW
  real,allocatable,dimension(:,:)   :: SNDFBEBKB
  real,allocatable,dimension(:,:)   :: SNDFBECHB
  real,allocatable,dimension(:,:)   :: SNDFBECHW
  real,allocatable,dimension(:,:)   :: TOXFBEBKB
  real,allocatable,dimension(:,:)   :: TOXFBECHB
  real,allocatable,dimension(:,:)   :: TOXFBECHW

  ! *** Propwash
  real,target,allocatable,dimension(:,:,:) :: SDF    !< Cohesive sediment concentration without fast settling classes     (g/m3), Also used for all sediments classes when NSEDFLUME > 0

  integer              :: NACTIVEWC
  integer, allocatable :: IACTIVEWC1(:)
  integer, allocatable :: IACTIVEWC2(:)

  type(WCVPOINTER), dimension(100) :: WCV
  
  ! *** Cross section fluxes
  real,allocatable,dimension(:,:,:,:) :: WC_UP        ! *** Mass of water column parameter moving upstream   (g)
  real,allocatable,dimension(:,:,:,:) :: WC_DN        ! *** Mass of water column parameter moving downstream (g)
  real,allocatable,dimension(:,:)     :: WC_QU        ! *** Discharge upstream   (m3)
  real,allocatable,dimension(:,:)     :: WC_QD        ! *** Discharge downstream (m3)
  
  real(RKD),allocatable,dimension(:) :: SNAPSHOTS

  ! *****************************************************************************************************************************
  ! Begin SEDZLJ variables
  ! PT- All REAL(RKD) :: values are explicitly specified as DOUBLE PRECISION for accuracy. 7/16/08
  logical :: LSEDZLJ
  integer :: IFWAVE
  integer :: IHTSTRT
  integer :: ISSLOPE
  integer :: ISWNWAVE
  integer :: ITBM
  integer :: KZ      
  integer :: NDYCOUNT
  integer :: NEQUIL
  integer :: NSEDFLUME
  integer :: NSICM
  integer :: NVAR_BED
  integer :: NWARNING
  integer :: NWVCOUN
  integer :: NWVCOUNT
  integer :: STWVNUM
  integer :: STWVTIM
  integer :: STINC

  real(RKD) :: DTSEDJ       !< SEDZLJ timestep in REAL(RKD) for better mass balance
  real(RKD) :: HPMIN        !< Minimum depth to compute shears for SEDZLJ
  real(RKD) :: MAXDEPLIMIT  !< The maximum limit of mass from 1st layer deposited on active bed layer.
  real(RKD) :: RHO          !< EFDC variable in units of MKS 1000. kg/m^3
  real(RKD) :: TACTM        !< Non-dimensional active layer thickness multiplier
  real(RKD) :: TAUCONST     
  real(RKD) :: WATERDENS    !< Add new parameter for water density rho_w
  real(RKD) :: ZBSKIN

  ! *** NEW parameterS FOR UPDATED SEDZLJ
  real(RKD) :: ALPHA_DEP
  real(RKD) :: ALAY_EXP
  real(RKD) :: ALAY_EXP_D50
  real(RKD) :: ALAY_CON
  real(RKD) :: ALPHA_D90
  real(RKD) :: BEDLOAD_CUTOFF
  real(RKD) :: CF_MIN
  real(RKD) :: EXP_HFACT_COH
  real(RKD) :: SIGMA_PDEP_COH
  real(RKD) :: SIGMA_PDEP_NON
  real(RKD) :: TAUB_MIN
  real(RKD) :: Z0_BED

  integer,target,allocatable,dimension(:,:) :: LAYERACTIVE   !(KB,LCM)
  integer,allocatable,dimension(:,:)     :: BLFLAG           !(LCM,NSEDS)
  integer,allocatable,dimension(:)       :: BEDMAP           !(LCM)
  integer,allocatable,dimension(:,:)     :: NCORENO          !(ICM,JCM)
                                                             
  real(RKD),allocatable,dimension(:)     :: ALPHA_PX          !(LCM)
  real(RKD),allocatable,dimension(:)     :: ALPHA_PY          !(LCM)
  real(RKD),allocatable,dimension(:,:)   :: ALPHA_RX          !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: ALPHA_RY          !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:)     :: ACTDEPA           !(KB)          2019-03-18   ACTIVE AND DEPOSITED EROSION RATE MULTIPLER FOR  E = A*tau**N
  real(RKD),allocatable,dimension(:)     :: ACTDEPMAX         !(KB)          2019-04-21   ACTIVE AND DEPOSITED MAX EROSION RATE
  real(RKD),allocatable,dimension(:)     :: ACTDEPN           !(KB)          2019-03-18   ACTIVE AND DEPOSITED EROSION RATE EXPONENT FOR   E = A*tau**N
  real(RKD),allocatable,dimension(:,:)   :: BLVEL             !(LCM,NSEDS)
  real(RKD),target,allocatable,dimension(:,:) :: BULKDENS     !(KB,LCM)
  real(RKD),target,allocatable,dimension(:,:) :: CBL          !(LCM,NSEDS) 
  real(RKD),allocatable,dimension(:)     :: D50               !(NSEDS)
  real(RKD),target,allocatable,dimension(:) :: D50AVG         !(LCM)
  real(RKD),allocatable,dimension(:,:)   :: DBL               !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:)     :: DISTAR            !(NSEDS)
  real(RKD),allocatable,dimension(:)     :: DWS               !(NSEDS)
  real(RKD),allocatable,dimension(:)     :: DWSIN             !(NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: DZBL              !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: DZBL_LAST         !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: EA                !(INCORE,KB)  2019-03-18   CORE INPLACE EROSION RATE MULTIPLER FOR  E = A*tau**N
  real(RKD),allocatable,dimension(:,:)   :: EN                !(INCORE,KB)  2019-03-18   CORE INPLACE EROSION RATE EXPONENT FOR   E = A*tau**N
  real(RKD),allocatable,dimension(:,:)   :: EBL               !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: ERATEND           !(NSICM,ITBM)
  real(RKD),allocatable,dimension(:,:,:) :: ERATE             !(KB,LCM,ITBM)
  real(RKD),allocatable,dimension(:)     :: ESUS              !(LCM)
  real(RKD),allocatable,dimension(:,:)   :: FWDIR             !(LCM,8)
  real(RKD),allocatable,dimension(:)     :: FWVTP             !(LCM)   
  real(RKD),allocatable,dimension(:)     :: HPCM              !(LCM)
  real(RKD),allocatable,dimension(:)     :: KPART             !(NTXM)
  real(RKD),allocatable,dimension(:,:)   :: MAXRATE           !(INCORE,KB)  2019-04-21   CORE INPLACE MAX EROSION RATE (CM/S)
  real(RKD),allocatable,dimension(:,:)   :: PSUS              !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: QBLFLUX           !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:)     :: SCND              !(NSICM)
  real(RKD),allocatable,dimension(:)     :: SEDDENS           !(INCORE)
  real(RKD),allocatable,dimension(:)     :: SH_SCALE          !(LCM)
  real(RKD),allocatable,dimension(:)     :: SSGI              !(NSEDS)   DSEDGMM array
  real(RKD),allocatable,dimension(:,:)   :: STWVHT            !(LCM,200)
  real(RKD),allocatable,dimension(:,:)   :: STWVTP            !(LCM,200)
  real(RKD),allocatable,dimension(:,:)   :: STWVDR            !(LCM,200)
  real(RKD),allocatable,dimension(:)     :: TAU               !(LCM)
  real(RKD),allocatable,dimension(:,:)   :: TAUCOR            !(KB,LCM)
  real(RKD),allocatable,dimension(:)     :: TAUCRITE          !(NSICM)
  real(RKD),allocatable,dimension(:)     :: TCRE              !(NSEDS)
  real(RKD),allocatable,dimension(:)     :: TAULOC            !(ITBM)
  real(RKD),allocatable,dimension(:)     :: TCRSUS            !(NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: TRANS             !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: TSED0             !(KB,LCM)
  real(RKD),allocatable,dimension(:)     :: TSEDT             !(LCM),TOXFTMP,TOXXTMP
  real(RKD),allocatable,dimension(:)     :: TSET0T            !(LCM)
  real(RKD),allocatable,dimension(:,:)   :: UBL               !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: USW               !(LCM,NSEDS)
  real(RKD),allocatable,dimension(:,:)   :: VBL               !(LCM,NSEDS)
  real(RKD),target,allocatable,dimension(:,:) :: DEP_SED_FLX  !(LCM,NSEDS)
  real(RKD),target,allocatable,dimension(:,:) :: ERO_SED_FLX  !(LCM,NSEDS)
  real(RKD),target,allocatable,dimension(:,:,:) :: PERSED     !(NSEDS,KB,LCM)
  real(RKD),target,allocatable,dimension(:,:)   :: TSED       !(KB,LCM)
  ! *** End SEDZLJ variables

  ! *** SEDIMENT CLASS VARIABLES                              
  type(SEDCLASS),allocatable,dimension(:) :: SEDS             !< All sediment classes (cohesive + non-cohesive)

  ! *****************************************************************************************************************************
  ! *** LAGRANGIAN TRAJECTORIES OF DRIFTERS
  integer,allocatable :: NACTT(:)
  integer,allocatable :: JSPD(:)
  
  real,      allocatable :: MOC(:)            !< MASS CONCENTRATION OF OIL [MG/L]

  real(RKD) ,allocatable :: XLA(:)            !< XLA(NPD) DRIFTER COORDINATES
  real(RKD) ,allocatable :: YLA(:)            !< YLA(NPD)
  real(RKD) ,allocatable :: ZLA(:)            !< ZLA(NPD)
  real(RKD) ,allocatable :: VELB(:)           !< VELB(1:KC+1) VETICAL DISTRIBUTION FOR U,V,W

  real(RKD) ,allocatable :: DLA(:)            !< DLA(NPD) INITIAL DEPTH OF DRIFTER
  real(RKD),allocatable :: XCOR(:,:)          !< XCOR(L,1:4) X POLYGON VERTICES, X(L,5): CENTROID X
  real(RKD),allocatable :: YCOR(:,:)          !< YCOR(L,1:4) Y POLYGON VERTICES, Y(L,5): CENTROID Y
  real(RKD),allocatable :: AREA(:)            !< AREA OF POLYGON LIJ
  real(RKD),allocatable :: HPLA(:)            !< TOTAL WATER DEPTH OF DRIFTERS
  real(RKD),allocatable :: BELVLA(:)          !< BOTTOM ELEVATION OF DRIFTERS
  real(RKD),allocatable :: DARE(:)            !< (1:NGRP) AREA OF OIL DRIFTER OF EACH GROUP   [M2]
  real(RKD),allocatable :: DDEN(:)            !< (1:NGRP) DENSITY OF OIL SPILL OF EACH GROUP[KG/M3]
  real(RKD),allocatable :: DVOL(:)            !< (1:NPD ) VOLUME OF EACH OIL DRIFTER [M3]
  real(RKD),allocatable :: DVOL0(:)           !< (1:NPD ) VOLUME OF EACH OIL DRIFTER [M3]
  real(RKD),allocatable :: GVOL(:)            !< (1:NGRP) TOTAL VOLUME OF OIL FOR EACH GROUP [M3]
  real(RKD),allocatable :: DRAT(:)            !< (1:NGRP) BIOGRADATION RATE OF EACH OIL DRIFTER OF EACH GROUP
  real(RKD),allocatable :: DTEM(:)            !< (1:NGRP) TEMPERATURE (C) OF EACH OIL DRIFTER OF EACH GROUP
  real(RKD),allocatable :: DVAP(:)            !< (1:NGRP) VAPOUR PRESSURE OF THE BULK LIQUID [PA]
  real(RKD),allocatable :: DVMO(:)            !< (1:NGRP) LIQUID'S MOLAR VOLUME [M3/MOL]
  
  integer,allocatable :: LLA(:)               !< LIJ(NP) INDEX OF POLYGON
  integer,allocatable :: KLA(:)               !< K OF DRIFTER NP
  integer,allocatable :: ISOILSPI(:)          !< (1:NGRP), 0: NORMAL DRIFTER/1:OIL DRIFTER
  integer(4)            :: LA_ZCAL            !< OPTION TO CALCULATE Z
  integer(4)            :: LA_PRAN            !< OPTION TO ADD A RANDOM MOVEMENT
  integer(4)            :: LA_DIFOP           !< OPTION FOR DIFFUSION COEFFICIENT 
  integer(4)            :: LA_FIRSTCALL
  integer(4)            :: DEPOP              !< OPTION FOR READING DEPTH  COLMUN IN DRIFTER.INP 
  integer(4)            :: XYZSCL

  real(8)               :: LA_BEGTI0          !< TIME BEGINNING FOR LAGRA
  real(8)               :: LA_ENDTI0          !< TIME ENDING FOR LAGRA
  real(RKD)             :: LA_FREQ            !< LA_FREQ IN HOURS FOR LAGRA OUTPUT
  real(RKD)             :: LA_HORDIF          !< HORIZONTAL DIFFUSION
  real(RKD)             :: LA_VERDIF          !< VERTICAL DIFFUSION
  real(RKD)             :: SLIPFACTOR         !< WALL SLIP FACTOR FOR PARTICLES NEAR WALLS (DIMENSIONLESS)  
  
  ! ISWAVE = 3
  real(RKD)             :: KSW                !< NIKURADSE ROUGHNESS (APPROXIMATE 2.5*D50)
  integer, allocatable  :: UMASK(:)
  integer, allocatable  :: VMASK(:)

  ! *****************************************************************************************************************************
  ! *** Begin MHK variables
  logical :: LMHK                                    !< Flag to indicate MHK is active
  integer :: MHKTYP                                  !< number of different types of MHK devices
  integer :: TCOUNT                                  !< number of cells with MHK devices
  integer :: UPSTREAM                                !< flag to decide if upstream or cell velocities are used
  
  real :: BETAMHK_P      !< Adjustable parameters that affect wake structure through K-e terms
  real :: BETAMHK_D      !< Adjustable parameters that affect wake structure through K-e terms
  real :: BETASUP_P      !< 
  real :: BETASUP_D      !< 
  real :: BETAVEG_P      !<    
  real :: BETAVEG_D      !< 
  real :: CE4MHK         !< Adjustable parameters that affect wake structure through K-e terms
  real :: PB_COEF        !< Adjustable parameters that affect wake structure through K-e terms
  
  integer :: OUTPUTFLAG                              !< this tells the Tecplot routine what to output
  integer,allocatable,dimension(:,:) :: IJLTURB      !< (TCOUNT,3) # of cells with turbines:1 = I-loc, 2 = J-loc, 3 = L-loc
  
  real,allocatable,dimension(:)      :: CTMHK        !< (MHKTYPE) # of turbine types
  real,allocatable,dimension(:)      :: CDSUP        !< (MHKTYPE)
  real,allocatable,dimension(:)      :: DENMHK       !< density of MHK devices (#/cell) (MHKTYPE)
  real,allocatable,dimension(:,:)    :: ESUP         !< energy dissipated from MHK support (LCM,TCOUNT)
  real,allocatable,dimension(:,:)    :: EMHK         !< energy dissipated from MHK device (LCM,TCOUNT)
  real,allocatable,dimension(:,:)    :: FXMHK        !< (LCM,KCM)
  real,allocatable,dimension(:)      :: FXMHKE       !< (LCM)
  real,allocatable,dimension(:,:)    :: FXSUP        !< (LCM,KCM)
  real,allocatable,dimension(:)      :: FXSUPE       !< (LCM)
  real,allocatable,dimension(:,:)    :: FYMHK        !< (LCM,KCM)
  real,allocatable,dimension(:)      :: FYMHKE       !< (LCM)
  real,allocatable,dimension(:,:)    :: FYSUP        !< (LCM,KCM)
  real,allocatable,dimension(:)      :: FYSUPE       !< (LCM)
  real,allocatable,dimension(:)      :: HEIGHTMHK    !< (MHKTYPE) # of turbine types
  real,allocatable,dimension(:)      :: HEIGHTSUP    !< (MHKTYPE) # of turbine types
  real,allocatable,dimension(:)      :: REFELEV      !< (MHKTYPE) # of turbine types
  real,allocatable,dimension(:,:)    :: PMHK         !< array that accumulates MHK power (LCM,KCM)
  real,allocatable,dimension(:,:)    :: PSUP         !< array that accumulates vegetative power (LCM,KCM)
  real,allocatable,dimension(:)      :: VMAXCUT      !< (MHKTYPE)
  real,allocatable,dimension(:)      :: VMINCUT      !< (MHKTYPE)
  real,allocatable,dimension(:)      :: WIDTHMHK     !< (MHKTYPE)
  real,allocatable,dimension(:)      :: WIDTHSUP     !< (MHKTYPE)
  real,allocatable,dimension(:,:)    :: ZMAXMHK      !< (MHKTYPE,LCM)
  real,allocatable,dimension(:,:)    :: ZMAXSUP      !< (MHKTYPE,LCM)
  real,allocatable,dimension(:,:)    :: ZMINMHK      !< (MHKTYPE,LCM)
  real,allocatable,dimension(:,:)    :: ZMINSUP      !< (MHKTYPE,LCM)
  ! *** End MHK variables

  ! *****************************************************************************************************************************
  ! *** OMP TEMPORARY VARIABLES,  LCM,KCM,NTHREADS
  real,allocatable,dimension(:,:,:) :: CD  
  real,allocatable,dimension(:,:,:) :: FQC !< EXTERNAL MODE FLOW BOUNDARY FORCING
  real,allocatable,dimension(:,:,:) :: FUHUD  
  real,allocatable,dimension(:,:,:) :: FVHUD  
  real,allocatable,dimension(:,:,:) :: FUHVD  
  real,allocatable,dimension(:,:,:) :: FVHVD  
  real,allocatable,dimension(:,:,:) :: FWUU  
  real,allocatable,dimension(:,:,:) :: FWVV  
  real,allocatable,dimension(:,:,:) :: FQCPAD  
  real,allocatable,dimension(:,:,:) :: QSUMNAD  
  real,allocatable,dimension(:,:,:) :: QSUMPAD  
  real,allocatable,dimension(:,:,:) :: POS  
  real,allocatable,dimension(:,:,:) :: UUUU  
  real,allocatable,dimension(:,:,:) :: VVVV  
  real,allocatable,dimension(:,:,:) :: WWWW  
  real,allocatable,dimension(:,:,:) :: DUU  
  real,allocatable,dimension(:,:,:) :: DVV  
  real,allocatable,dimension(:,:,:) :: CONQ
  real,allocatable,dimension(:,:,:) :: CON2
  real,allocatable,dimension(:,:,:) :: CMAX  
  real,allocatable,dimension(:,:,:) :: CMIN  
  real,allocatable,dimension(:,:,:) :: CONTD  
  real,allocatable,dimension(:,:,:) :: CONTMN  
  real,allocatable,dimension(:,:,:) :: CONTMX  
  real,allocatable,dimension(:,:,:) :: WQBCCON  
  real,allocatable,dimension(:,:,:) :: WQBCCON1  

  ! *****************************************************************************************************************************
  ! *** NetCDF Variables
  character(5)   :: PROJ
  character(10)  :: BASEDATE
  character(10)  :: SDATE
  character(8)   :: BASETIME
  character(120) :: NCFILE
  character(120) :: TEMFILE

  integer(IK4) :: NCDFOUT           
  integer(IK4) :: DEFLEV            
  integer(IK4) :: BLK
  integer(IK4) :: IMN
  integer(IK4) :: JMN
  integer(IK4) :: IMX
  integer(IK4) :: JMX
  integer(IK4) :: HEMI
  integer(IK4) :: UTMZ
  integer(IK4) :: ROTA
  integer(IK4) :: ISSGLFIL
  integer(4)   :: IS_NC_OUT(50)           ! *** Options for NetCDF output: 0 = NO, 1 = YES
  
  real      :: HREST
  real(RKD) :: TBEGINC
  real(RKD) :: TBEGNCDF
  real(RKD) :: TENDNCDF
  real(RKD) :: NCFREQ                  ! *** NetCDF output interval in minutes
  real(RKD) :: NCDFSHOT = 1.
    
  real,   allocatable :: XCR(:,:,:)
  real,   allocatable :: YCR(:,:,:)
  
  ! *** HIGH FREQUENCY OUTPUT
  integer :: HFREOUT,NSUBSET
  real,       allocatable :: UK(:),VK(:)           
  type(CELL), allocatable :: HFREGRP(:)
  integer,    allocatable :: IJHFRE(:),NPNT(:)  
  real(RKD),  allocatable :: HFREDAYBG(:),HFREDAYEN(:)
  real(RKD),  allocatable :: HFREDUR(:),HFREDAY(:),HFREMIN(:)

  ! *** NEW VARIABLES FOR QCTL NQCTYP = 3 & 4 
  integer,allocatable,dimension(:) :: NLOWCHORD
  real,allocatable,dimension(:)    :: LOWCHORDU
  real,allocatable,dimension(:)    :: LOWCHORDV
  real,allocatable,dimension(:,:)  :: SAVESUB
  real,allocatable,dimension(:,:)  :: SAVESVB

  real,allocatable,dimension(:) :: WINDSXY         !< Wind sheltering adjusted cell rotation matrix
  real,allocatable,dimension(:) :: WINDSYX         !< Wind sheltering adjusted cell rotation matrix
  real,allocatable,dimension(:) :: WINDSXX         !< Wind sheltering adjusted cell rotation matrix
  real,allocatable,dimension(:) :: WINDSYY         !< Wind sheltering adjusted cell rotation matrix
  
  ! *** TIME BLOCKS
  integer :: NWNDMAP
  integer :: NATMMAP
  integer :: NICEMAP
  real,allocatable,dimension(:) :: TWNDMAPBEG
  real,allocatable,dimension(:) :: TWNDMAPEND
  real,allocatable,dimension(:) :: TATMMAPBEG
  real,allocatable,dimension(:) :: TATMMAPEND
  real,allocatable,dimension(:) :: TICEMAPBEG
  real,allocatable,dimension(:) :: TICEMAPEND
    
  ! *** Food chain linkage
  integer :: ISFDCH
  integer :: NFDCHZ 

  real :: HBFDCH    
  real :: TFCAVG    
  
  integer,allocatable,dimension(:)     :: MFDCHZ

  ! *** GOTM turbulence variables
  integer :: ISGOTM                                  !< Flag for application of GOTM
  Real, allocatable :: tke3d(:,:)                    !< Turbulence kinetic energy (m2/s2)
  Real, allocatable :: tke3d1(:,:)                   !< Turbulence kinetic energy (m2/s2) one level time step
  Real, allocatable :: eps3d(:,:)                    !< Dissipation rate
  Real, allocatable :: eps3d1(:,:)                   !< Dissipation rate one level time step
  Real, allocatable :: gl3d (:,:)                    !< Turbulence length scale (m)
  integer :: ICALSS                                  !< Calculate Shear Frequency Squared (SS) Methods             
  integer :: ICALNN                                  !< Calculate Buoyancy Frequency Squared (NN) Methods   
  integer :: IFRICTION                               !< Calculate Bottom Friction Velocity Methods  
  integer :: charnock                                !< Calculate surface roughness based on charnock formulation 
                                                     
  integer :: iGOTM_Test                              !< Flag for application of surface forcing (KPP GOTM test-case)
  real    :: SurfaceShearX, SurfaceShearY, HeatFlux
  

  ! *** Temporary use only for debugging and other testing
  integer :: LB, IIB, JJB, KKB
  
  integer,allocatable,dimension(:,:) :: DOMAIN_LIST  !< Debugging cell list of each sub-domain
  integer,allocatable,dimension(:)   :: ITEST        !< Created and allocated for debug assignments - Not used for routine runs
  real,allocatable,dimension(:)      :: TEST1D       !< Created and allocated for debug assignments - Not used for routine runs
  real,allocatable,dimension(:,:)    :: TEST1        !< Created and allocated for debug assignments - Not used for routine runs
  real,allocatable,dimension(:,:)    :: TEST2        !< Created and allocated for debug assignments - Not used for routine runs
  
  ! *** COARE36 Variables
  real :: PBLZ                                       !< Planetary boundary layer height (m)
  real :: TEM_HRZ                                    !< Temperature/RH sensor height (m)
  real :: COARE_NITS                                 !< Number of iterations of COARE algorithm
  real :: TCTMSR                                     !< Time conversion from seconds to user units 
  real,allocatable,dimension(:) :: CDCOARE           !< COARE 3.6 Wind drag coefficient
  real,allocatable,dimension(:) :: HSCOARE           !< Sensible heat flux
  real,allocatable,dimension(:) :: HLCOARE           !< Latent heat flux
  
  ! *** Legacy Output - Time Series
  integer :: ISTMSR          !< Activate flag
  integer :: JSTMSR          !< 
  integer :: NSMTS           !< 
  integer :: NSMTSM          !< 
  integer :: MLTMSR          !< 
  integer :: MLTMSRM         !<  
  integer :: MTSSTSPM        !<  
  integer :: NBTMSR          !< Iteration to begin output
  integer :: NSTMSR          !< Iteration to stop output
  integer :: NCTMSR          !< 
  integer :: NWTMSR          !< Write ferquency for time series


  character*20,allocatable,dimension(:) :: CLTMSR

  integer,allocatable,dimension(:)     :: ILTMSR     !< 
  integer,allocatable,dimension(:)     :: JLTMSR     !< 
  integer,allocatable,dimension(:)     :: NTSSSS     !<      
  integer,allocatable,dimension(:)     :: MTMSRA     !< 
  integer,allocatable,dimension(:)     :: MTMSRC     !< 
  integer,allocatable,dimension(:)     :: MTMSRP     !< 
  integer,allocatable,dimension(:)     :: MTMSRQ     !< 
  integer,allocatable,dimension(:)     :: MTMSRQE    !< 
  integer,allocatable,dimension(:)     :: MTMSRU     !< 
  integer,allocatable,dimension(:)     :: MTMSRUE    !< 
  integer,allocatable,dimension(:)     :: MTMSRUT    !< 

  integer,allocatable,dimension(:)     :: MTSSTSP    !< Legacy output - time series
  real,allocatable,dimension(:,:)   :: TSSTOP   ! time series
  real,allocatable,dimension(:,:)   :: TSSTRT   ! time series
 
  !*** Residuals/CALMMT
  integer :: JSRSPH(10)      !< 
  integer :: JSRESTR         !< 
  integer :: ISSSMMT         !< 
  
  real :: RESSTEP        !< Model linkage/mean mass transport linkage timestep (s)

  integer,allocatable,dimension(:) :: ISPRE      !< Residual current method for ISPE - EAST
  integer,allocatable,dimension(:) :: ISPRN      !< Residual current method for ISPN - NORTH
  integer,allocatable,dimension(:) :: ISPRS      !< Residual current method for ISPS - SOUTH
  integer,allocatable,dimension(:) :: ISPRW      !< Residual current method for ISPW - WEST 
  
  real,allocatable,dimension(:,:)     :: ABEFF       
  real,allocatable,dimension(:)       :: GWLPF         !< Low pass filter      
  real,allocatable,dimension(:,:)     :: VTLPF         !< Low pass filter
  real,allocatable,dimension(:,:)     :: WIRT          !< Low pass filter
                                                
  real,allocatable,dimension(:,:)     :: ABLPF         !< Low pass filter of AB
  real,allocatable,dimension(:,:)     :: AHULPF        !< Low pass filter of AH in the U direction
  real,allocatable,dimension(:,:)     :: AHVLPF        !< Low pass filter of AH in the V direction
  real,allocatable,dimension(:,:,:)   :: DYELPF        !< Low pass filter of DYE
  real,allocatable,dimension(:)       :: EVPGLPF       !< Low pass filter of EVAPGW
  real,allocatable,dimension(:)       :: EVPSLPF       !< Low pass filter of EVAPSW
  real,allocatable,dimension(:)       :: QSUMELPF      
  real,allocatable,dimension(:,:)     :: QSUMLPF       
  real,allocatable,dimension(:)       :: QCHNULP       !< LOW PASS FILTER OF U DOMINANT SUBGRID SCALE CHANNEL CONNETION FLOW L*L*L/T
  real,allocatable,dimension(:)       :: QCHNVLP       !< LOW PASS FILTER OF V DOMINANT SUBGRID SCALE CHANNEL CONNETION FLOW L*L*L/T
  real,allocatable,dimension(:)       :: RAINLPF       !< ACCUMULATED RAINFALL (M3/S)
  real,allocatable,dimension(:)       :: RINFLPF       
  real,allocatable,dimension(:,:,:)   :: SEDBLPF       !< WATER COLUMN COHESIVE SEDIMENT FLUX M/L*L*T
  real,allocatable,dimension(:,:)     :: SEDBTLPF      !< LOW PASS FILTER OF SEDBO
  real,allocatable,dimension(:,:,:)   :: SEDLPF        !< NONCOHESIVE SEDIMENT CONCENTRATION IN WATER COLUMN
  real,allocatable,dimension(:,:)     :: SEDTLPF       !< LOW PASS FILTERED VALUE OF SEDT
  real,allocatable,dimension(:,:)     :: SFLLPF        
  real,allocatable,dimension(:,:,:)   :: SNDBLPF       !< NONCOHESIVE SEDIMENT VERTICAL FLUX IN WATER COLUMN M/L*L*T
  real,allocatable,dimension(:,:)     :: SNDBTLPF      !< BED MASS 
  real,allocatable,dimension(:,:,:)   :: SNDLPF        !< NON-COHESIVE CLASS WATER COLUMN CONCENTRATION ACCUMULATOR
  real,allocatable,dimension(:,:)     :: SNDTLPF       !< LOW PASS FILTER OF SNDT
  real,allocatable,dimension(:,:)     :: TEMLPF        
  real,allocatable,dimension(:,:,:)   :: TOXBLPF       !< WATER COLUMN TOXIC CONTAMIANT FLUX M/L*L*T
  real,allocatable,dimension(:,:,:)   :: TOXLPF        !< TOTAL TOXIC PARTICULATE PHASE FRACTION IN BED
  real,allocatable,dimension(:,:,:,:) :: TXPFLPF       !< LOW PASS FILTER OF TOXIC PARTICULATE FRACTION
  real,allocatable,dimension(:)       :: UELPF         !< Rotated depth average U
  real,allocatable,dimension(:,:)     :: UHLPF         !< Layer specific U
  real,allocatable,dimension(:,:)     :: UTLPF         !< 
  real,allocatable,dimension(:)       :: VELPF         !< Rotated depth average U 
  real,allocatable,dimension(:,:)     :: VHLPF         !< 
  real,allocatable,dimension(:,:)     :: UIRT  
  real,allocatable,dimension(:,:)     :: VIRT  
  real,allocatable,dimension(:,:)     :: WTLPF         !<        
  
  
  ! *** WASP/ICM/RCA linkage variables
  character*20 :: HYDFIL
  integer(4)   :: IHL_HANDLE

  integer :: ISICM       !< icm
  integer :: IAUXICM     !< ICM
  integer :: ISTICM      !< ICM
  integer :: NCICM       !< ICM 
  integer :: NCRCA1      !< ICM
  integer :: NINTFL      !< ICM
  integer :: NQINFLM     !< ICM
  integer :: NTSMMT      !< ICM
  
  integer,allocatable,dimension(:)     :: KCEFDC     !< ICM
  integer,allocatable,dimension(:)     :: LCEFDC     !< ICM
  integer,allocatable,dimension(:)     :: KFEFDC     !< ICM
  integer,allocatable,dimension(:)     :: LFEFDC     !< ICM
  
  real,allocatable,dimension(:)     :: HTMP   

  integer :: NMMT
  integer :: NWASPOUT
  integer :: NTIME = 1
  
  integer :: ISWASP      !< WASP
  integer :: ISWASPD     !< WASP
  integer :: IBEDV       !< WASP
  integer :: IQOPT       !< WASP
  integer :: ISDHD       !< WASP
  integer :: ISDICM      !< ICM
  integer :: ISNKH       !< WASP
  integer :: ISRCA       !< RCA 
  integer :: IVOPT       !< WASP
  integer :: JSWASP      !< Initialization flag - WASP
  
  integer :: LWASP
  integer :: LALT        !< Number of WASP cells
  integer :: LCLT        !< LCLT = LALT + 1

  integer :: NFICM       !< ICM
  integer :: NFIELD      !< WASP
  integer :: NRFLD       !< WASP

  real :: CONVQ          !< WASP
  real :: CONVR          !< WASP
  real :: CONVV          !< WASP
  real :: DEPSED         !< WASP
  real :: DEXP           !< WASP
  real :: DMULT          !< WASP
  real :: SCALQ          !< WASP
  real :: SCALV          !< WASP
  real :: SEDIFF         !< WASP
  real :: TDINTS         !< WASP
  real :: VEXP           !< WASP
  real :: VMULT          !< WASP

  ! *** Connections WASP/LT Models
  integer,allocatable,dimension(:)     :: ILLT       !< 
  integer,allocatable,dimension(:)     :: JLLT       !< 
  integer,allocatable,dimension(:,:)   :: LCHNC      !< 
  integer,allocatable,dimension(:)     :: NCHNC      !< 

  integer,allocatable,dimension(:)     :: LNCLT      !< List of cells used for WASP/LT models - Cell to the North
  integer,allocatable,dimension(:)     :: LSCLT      !< List of cells used for WASP/LT models - Cell to the South
  integer,allocatable,dimension(:)     :: LECLT      !< List of cells used for WASP/LT models - Cell to the East
  integer,allocatable,dimension(:)     :: LWCLT      !< List of cells used for WASP/LT models - Cell to the West

  real,allocatable,dimension(:)  :: TPCOORDS
  real,allocatable,dimension(:)  :: TPCOORDW
  real,allocatable,dimension(:)  :: TPCOORDE
  real,allocatable,dimension(:)  :: TPCOORDN
 
  real,allocatable,dimension(:)     :: HLPF          !< LOW PASS FILTER OF DEPTH - wasp
  real,allocatable,dimension(:,:)   :: SALLPF        !< LOW PASS FILTER SALINITY  wasp
  real,allocatable,dimension(:)     :: SWB           !< WASP FLAG FOR OPEN BOUNDARIES
  real,allocatable,dimension(:,:)   :: ULPF          !< 
  real,allocatable,dimension(:,:)   :: UVPT  
  real,allocatable,dimension(:,:)   :: VLPF          !< 
  real,allocatable,dimension(:,:)   :: VVPT  
  real,allocatable,dimension(:,:)   :: WVPT 
  real,allocatable,dimension(:,:)   :: WLPF          !< 
  
  real(rkd),allocatable,dimension(:) :: WASPTIME     !< 
  
#ifdef GNU  
  character(11) :: FMT_BINARY                        !< 'UNFORMATTED' file format
#else
  character(6) :: FMT_BINARY                         !< 'BINARY' file format
#endif

  END MODULE
