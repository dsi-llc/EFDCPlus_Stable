! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE HYDSTRUCMOD

  ! *** PURPOSE:  calculate flow at the hydraulic structure BCS
  ! ***           including culverts, weirs and sluice gates
  ! *** Based on the paper by NL DILL: Estuarine and Coastal Modeling (2011)
  ! *** Modeling Hydraulic Control Structures in Estuarine Environments with EFDC
  ! *** DATE: August 2015
  ! ***
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-12       P M CRAIG         Improved hydraulic structures flow calculation
  !    2015-11       D H CHUNG         Implemented new hydraulic structures based on DILL
  
  use GLOBAL 
  use INFOMOD
  use Variables_MPI
  
  implicit none
  
  type OPERATIONDATA                !< Trigger for operation of gates/pumps
    integer :: STATE                !< Gate position
                                    !    0 = Completely Closed/Off
                                    !    1 = Fully Open/On
                                    !    2 = Opening
                                    !    3 = Closing
    integer :: ID                   !< Index of rating curve/mask for control parameter
    real :: LEVEL                   !< Trigger level
    real :: HEIGHT                  !< Opening height (m), for upward opening
    real :: WIDTH                   !< Opening width (m), for sideward opening
    real :: SILL                    !< Sill level change (m), for downward opening
    real :: FLOW                    !< Flow rate (m3/s), for pumps
    real :: RATE                    !< Rate of opening/closing
    !INTEGER :: UNITS                !< Number of gates, pump units  (NOT USED)
  end type 

  type CONTROLRULES                 !< Operation rules for of gates/pumps
    integer :: PARAM                !< Mask for control parameters
    integer :: STATE                !< Structure State, 0 = OFF, 1 = ON
    integer :: NTRGON, NTROFF  
    type(OPERATIONDATA),allocatable :: TRGON(:), TROFF(:)
  end type   
  
  type CONTROLSERIES                !< Gate opening time-series
    integer :: ITYPE                !< Structure control type: 0 = UNCONTROLED, 1 = TS, 2 = US, 3 = HDIFF
    integer :: COUNT                !< # Points
    integer :: IPOS = 1             !< Index of current time
    integer :: PARAM                !< Mask for control parameters
    integer,allocatable :: ID(:)
    
    real, allocatable :: HEIGHT(:), WIDTH(:), SILL(:), FLOW(:), RATE(:) 
    real(RKD) :: TMUL, TADD
    real(RKD),allocatable :: TIME(:)
    !INTEGER,allocatable :: NUM(:)  ! NUMBER OF GATES/PUMP UNITS (NOT USED)
  end type 
  
  type STRCONTROL                    !< Info on the control of hydraulic structure
    integer :: ITYPE                 !< Structure control type:
                                     !<   0 - Uncontroled
                                     !<   1 - Time series
                                     !<   2 - Upstream
                                     !<   3 - Head difference
    integer :: ID                    !< ID of associated control data
    integer :: IREFUP                !< Upstream   reference cell. Allows reference cell be different from US cell
                                     !<                            For navigation locks, this is the upstream lock structure I
                                     !<                            For the most upstream gate use I = 0
    integer :: JREFUP                !< Upstream   reference cell. Allows reference cell be different from US cell
                                     !<                            For navigation locks, this is the upstream lock structure J
                                     !<                            For the most upstream gate use J = 0
    integer :: IREFDN                !< Downstream reference cell. Allows reference cell be different from DS cell
    integer :: JREFDN                !< Downstream reference cell. Allows reference cell be different from DS cell
    integer :: LUR                   !< Upstream   reference L Index
    integer :: LDR                   !< Downstream reference L Index
    integer :: SUBID                 !< Sub index of associated control data
    real(RKD) :: TSTART, TOPEN
    type(OPERATIONDATA) :: CUR, LIM  !< Initial condition (current value) & limit
    real :: PREV                     !< Previous control level
  end type   

  type(STRCONTROL),save,allocatable,dimension(:)    :: HSCTL_GL  !< Hydraulic sructure control
  type(STRCONTROL),save,allocatable,dimension(:)    :: HSCTL     !< Hydraulic sructure control
  type(CONTROLSERIES),save,allocatable,dimension(:) :: QCTLSER   !< Control time-series
  type(CONTROLRULES),save,allocatable,dimension(:)  :: QCTRULES  !< Control rules

  integer :: NQCTLSER, NQCTRULES

  type(STRCONTROL),save,allocatable,dimension(:)        :: WITH_RET_CTL_GL    !< Withdrawal/return control (Global)
  type(STRCONTROL),save,allocatable,dimension(:)        :: WITH_RET_CTL       !< Withdrawal/return control (Local)

  integer :: HSOPTION
  integer,save,allocatable,dimension(:) :: USCELL    ! *** Save the original upstream cell
  integer,save,allocatable,dimension(:) :: DSCELL    ! *** Save the original downstream cell
    
  contains 

  ! *** Determine the flows for the current time steps based on all approaches for flow contro
  SUBROUTINE COMPUTE_HSFLOW(NCTL)
  
    integer,save,allocatable,dimension(:) :: NDIRCULV   ! *** Flow direction counter
    integer,save,allocatable,dimension(:) :: NLIMITQ    ! *** Flow limitation counter
    real,save,allocatable,dimension(:)    :: CULVQ      ! *** Last culvert flow rate

    integer, intent(IN) :: NCTL
    integer             :: LU, LD, LUR, LDR, K, KM, IHYD, IBLOCK, ID, IFLAG, M1, M2, NWR, STATE
    
    real :: ZHU, ZHD, ZHUR, ZHDR, HUD, HDD, S, SCR, HCR, VCR, VN, QHS, QHS1, QHS2, USINV, DSINV, DIA, CDD
    real :: QCR, YN, TOPZ, TOLZ, DELZ, CUMZ, HB, W, RATIO, Q1, Q2, QMAX
    real :: FAREA, WETPER, HRAD, KCON, TWID, HVAL
    real :: HOPEN, HOPEN2, WOPEN, WOPEN2, SOPEN, SOPEN2, TOPEN1, TOPEN2
    real :: CVAL, LVAL, RVAL, FOPEN
    real :: HSZ(KC)
    
    if( .not. allocated(NDIRCULV) )then
      allocate(NDIRCULV(NQCTL))
      allocate(NLIMITQ(NQCTL))
      allocate(CULVQ(NQCTL))
      NDIRCULV = 100
      NLIMITQ = 0
      CULVQ = 0.
    endif

    if( ISDYNSTP == 0 )then
      DELT = DT
    else
      DELT = DTDYN
    endif

    HSOPTION = 1                         ! *** Option for sluice gate, change later
                                         
    QCTLT(1:KC,NCTL,1:2) = 0             ! *** Reset for every time step
    
    IHYD = HYD_STR(NCTL).NQCTLQ          ! *** Hydraulic structure defintion parameters 
    if( HYD_STR(NCTL).NQCTYP /= 9 )then
      HOPEN = HS_HEIGHT(IHYD)            ! *** Default opening height is the structure dimension
      WOPEN = HS_WIDTH(IHYD)             ! *** Default opening width is the structure dimension
      SOPEN = 0.                         ! *** Default is no sill level change
    endif
    
    ! *** Get settings for strutures that are controlled for the current time step
    if( HSCTL(NCTL).ITYPE == 1 .or. HSCTL(NCTL).ITYPE == 4 )then    
      ! *** Structure is controlled by time-series
      ID = HSCTL(NCTL).ID
      IBLOCK = 0
      if( HSCTL(NCTL).ITYPE == 4 )then
          IBLOCK = 1
          IHYD = 0
      endif
      call GATE_OPENING_INTERP(QCTLSER(ID), TIMESEC, HOPEN, WOPEN, SOPEN, IBLOCK)
      
    elseif( HSCTL(NCTL).ITYPE == 2 .or. HSCTL(NCTL).ITYPE == 3 )then
      ! *** Structure is controlled by operation rules
      LU = LIJ(HSCTL(NCTL).IREFUP, HSCTL(NCTL).JREFUP)
      ZHU = BELV(LU) + HP(LU)
      
      if( HSCTL(NCTL).ITYPE == 3 )then
        ! *** Structure is controlled by operation rules depending on difference of water surface elevations
        LD = LIJ(HSCTL(NCTL).IREFDN, HSCTL(NCTL).JREFDN)      
        ZHD = BELV(LD) + HP(LD)
        HVAL = ZHU - ZHD
      else
        ! *** Structure is controlled by operation rules depending on specific water surface elevation
        HVAL = ZHU
      endif
      
      ID = HSCTL(NCTL).ID
      if( (QCTRULES(ID).PARAM .and. 16) == 16 )then
        ! *** Pump flows
        CVAL = HSCTL(NCTL).CUR.FLOW

        call PUMP_OPERATION_RULES(HSCTL(NCTL), HVAL, CVAL, RVAL)
        QHS = RVAL
        HSCTL(NCTL).CUR.FLOW = RVAL
        LD = LC
        GOTO 1000
      else  
        ! *** For gates, lookup table should be handled in CALQVS.F90
        if( (QCTRULES(ID).PARAM .and. 1) == 1 )then
          ! *** Gate opening upward (raising/lowering sill)
          CVAL = HSCTL(NCTL).CUR.HEIGHT
          LVAL = HSCTL(NCTL).LIM.HEIGHT
          call GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
          HOPEN = RVAL
          HSCTL(NCTL).CUR.HEIGHT = RVAL
          
        elseif( (QCTRULES(ID).PARAM .and. 2) == 2 )then
          ! *** Gate opening sideward
          CVAL = HSCTL(NCTL).CUR.WIDTH
          LVAL = HSCTL(NCTL).LIM.WIDTH
          call GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
          WOPEN = RVAL
          HSCTL(NCTL).CUR.WIDTH = RVAL
          
        elseif( (QCTRULES(ID).PARAM .and. 4) == 4 )then
          ! *** Gate opening downward (raising/lowering gate)
          CVAL = HSCTL(NCTL).CUR.SILL
          LVAL = HSCTL(NCTL).LIM.SILL
          call GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
          SOPEN = RVAL
          HSCTL(NCTL).CUR.SILL = RVAL
        endif
      endif

    else
      ! *** Structure is uncontrolled, do nothing
    endif
    
    ! *** Now set appropriate head for structure based on hydraulic structure cells
    LU = LIJ(HYD_STR(NCTL).IQCTLU,HYD_STR(NCTL).JQCTLU)
    LD = LIJ(HYD_STR(NCTL).IQCTLD,HYD_STR(NCTL).JQCTLD)
    ZHU = BELV(LU) + HP(LU)
    USINV = HS_USELEV(IHYD) + SOPEN         ! 2017-11-06, NTL: ACCOUNTS FOR SILL LEVEL CHANGE
    
    if( LD == 0 )then
      ! *** Downstream cell is outside domain
      LD = LC
      ZHD = -9999
    else
      ! *** Downstream cell is active, allow reverse flow
      ZHD = BELV(LD) + HP(LD)
      DSINV = HS_DSELEV(IHYD) + SOPEN     ! 2017-11-06, NTL: ACCOUNTS FOR SILL LEVEL CHANGE
      if( ZHU < ZHD .and. HS_REVERSE(IHYD) == 1 )then
        if( NDIRCULV(NCTL) < 6 .and. HYD_STR(NCTL).NQCTYP == 5 )then
          NDIRCULV(NCTL) = NDIRCULV(NCTL)+1
          QHS = 0.5*CULVQ(NCTL)
          CULVQ(NCTL) = QHS
          HDD = ZHD - DSINV                
          call CROSS_SECTION(IHYD, HDD, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID) 
          GOTO 1000
        endif
        NDIRCULV(NCTL) = 0
        call SWAP_USDS(NCTL,LU,LD,ZHU,ZHD,USINV,DSINV)
      endif
    endif
       
    ! *** Determine if elevations/depth allow any flow
    HUD = ZHU - USINV
    if( HYD_STR(NCTL).NQCTYP /= 9 .and. (HUD <= 0. .or. HP(LU) < HWET .or. (HS_REVERSE(IHYD) == 0 .and. ZHD > ZHU)) )then
      return
    endif    
    
    if( HYD_STR(NCTL).NQCTYP == 5  )then
      ! *** CULVERT 
      ! *** FORMULA FROM N L DILL PAPER:
      
      call CROSS_SECTION(IHYD, HUD, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID, KCON)
      
      HCR =  HUD - 0.5*HRAD               ! *** Critical depth
      HCR = MIN(HCR,DIA)
      VCR = SQRT(9.81*HRAD)
      QCR = VCR*FAREA
      SCR = (QCR/KCON)**2   
      
      ! *** Downstream head
      if( LD == LC )then
        ! *** Flow leaves domain (assume critical depth at outlet)
        HDD = HCR
      else
        HDD = ZHD - DSINV
      endif      

      ! *** Pipe slope
      S = (HS_USELEV(IHYD)-HS_DSELEV(IHYD))/HS_LENGTH(IHYD)
      if( S < 1.E-6 )then
        ! *** Pipe zero or adverse slope
        if( LD == LC )then
          ! *** Flow leaves domain (assume critical depth at outlet)
          S = HCR/HS_LENGTH(IHYD)
        else
          S =  (ZHU-ZHD)/HS_LENGTH(IHYD)
        endif      
      else
        if( LD == LC )then
          ! *** Flow leaves domain
          S = ABS(S)
        else
          S = MIN(S,(ZHU-ZHD)/HS_LENGTH(IHYD))  ! *** Minimum of energy grade or pipe slope
        endif      
      endif

      if( HDD > DIA .or. HUD > 1.5*DIA  )then
        ! *** Full flow         
        S = (ZHU-ZHD)/HS_LENGTH(IHYD)
        QHS = KCON*SQRT(S)
        
      elseif( S >= SCR  )then
        ! *** Supercritical flow control at inlet LD
        QHS = QCR
          
      elseif( S < SCR .and. HDD >= HCR  )then
        ! *** Subcritical flow tailwater control LD
        DELZ = MIN(HDD,HUD)
        call CROSS_SECTION(IHYD, DELZ, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID)
        VN = HRAD**0.6667 * S**0.5 / HS_MANN(IHYD)
        QHS = VN*FAREA     ! *** Q AT HDD
          
      elseif( HDD < HCR .and. S < SCR )then
        ! *** Subcritical flow control at outleT LU
        YN = HUD - HRAD**1.3333*S/(19.62*HS_MANN(IHYD)**2)
        call CROSS_SECTION(IHYD, YN, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID)
        VN = (FAREA/WETPER)**0.6667*S**0.5/HS_MANN(IHYD)
        QHS = VN*FAREA
      endif
      
      ! *** limit the rate of growth/reduction (5%)
      if( CULVQ(NCTL) > 0.01 )then
        if( QHS < 0.95*CULVQ(NCTL) )then
          QHS = 0.95*CULVQ(NCTL)
        elseif( QHS > 1.05*CULVQ(NCTL) )then
          QHS = 1.05*CULVQ(NCTL)
        endif
      elseif( QHS > 0.015 )then
        QHS = 0.0125
      endif
      
      ! *** Limit outflow from culvert cell
      QMAX = 0.01*DXYP(LU)/DELT
      if( HUD < 2.*HDRY ) QMAX = 0.1*DXYP(LU)*HDRY/DELT
      !IF( 10.*DIA > DXYP(LU) ) QMAX = 0.1*QMAX  ! *** LIMIT CONDITION WITH SMALL CELLS RELATIVE TO CULVERT SIZES
      if( (HP(LU)-H1P(LU)) < -0.01 )then
        QMAX = CULVQ(NCTL)/(H1P(LU)-HP(LU))*0.01
      endif
      if( QHS > QMAX )then
        NLIMITQ(NCTL) = NLIMITQ(NCTL) + 1
        if( NLIMITQ(NCTL) == 1 .or. MOD(NLIMITQ(NCTL),1000) == 0 )then
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit, "(A,I10,F12.5,I5,2F10.4)" )'CULVERT FLOW LIMITATION (NITER,TIMEDAY,NCTL,QHS,QMAX): ',NITER,TIMEDAY,NCTL,QHS,QMAX
          close(mpi_efdc_out_unit)
        endif
        QHS = QMAX
      endif

      NDIRCULV(NCTL) = NDIRCULV(NCTL)+1
      CULVQ(NCTL) = QHS

    elseif( HYD_STR(NCTL).NQCTYP == 6  )then 
      ! *** SLUICE GATE
      HDD = ZHD - USINV
      W  = WOPEN              ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      HB = HOPEN              ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      DIA = MIN(HB,HUD)
      
      if( HB <= 1E-6 .or. W <= 1E-6 ) return
      
      ! *** APPROACH FROM N L DILL PAPER:
      ! *** QHS = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5   C = 1.45  ! *** Super-critical weir flow
      ! *** QHS = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))   C = 0.61  ! *** Sub-critical weir flow
      ! *** QHS = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)          C = 0.60  ! *** Free sluice flow
      ! *** QHS = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))    C = 0.50  ! *** Submerged orifice flow
          
      ! *** Comparison of HUD & HDD
      if( HDD <= 0.64*HUD  )then
        ! *** Supercritical flow
        if( HUD >= 1.25*HB )then
          QHS = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)           ! *** Free sluice flow
          
        elseif( HUD > HB .and. HUD < 1.25*HB )then
          ! *** Weighted average weir/orifice critical flow
          RATIO = (HUD-HB)/(0.25*HB)
          Q1 = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)            ! *** Free sluice flow
          Q2 = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5     ! *** Super-critical weir flow
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        else
          ! *** Broad crested weir flow
          QHS = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5    ! *** Super-critical weir flow
          
        endif
        
      elseif( 0.64*HUD < HDD .and. HDD < 0.68*HUD )then
        ! *** Weighted average of supercritical & subcritical
        RATIO = (HDD-0.64*HUD)/(0.04*HUD)
        if( HUD >= 1.25*HB )then
          ! *** Weighted average super/sub critical orifice flow
          Q1 = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)            ! *** Free sluice flow
          Q2 = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))      ! *** Submerged orifice flow
          QHS = RATIO*Q2 + (1.-RATIO)*Q1
          
        elseif( HUD > HB .and. HUD < 1.25*HB )then
          ! *** Weighted average weir/orifice 
          Q1 = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)            ! *** Free sluice flow
          Q2 = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))     ! *** Sub-critical weir flow
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        else
          ! *** Weighted average weir/orifice 
          Q1 = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))     ! *** Sub-critical weir flow
          Q2 = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5     ! *** Super-critical weir flow
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        endif
          
      elseif( HDD >= 0.68*HUD )then
        ! *** Subcritical flows
        if( HUD >= 1.25*HB )then
          ! *** Weighted average super/sub critical orifice flow
          QHS = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))     ! *** Submerged orifice flow
          
        elseif( HUD > HB .and. HUD < 1.25*HB )then
          ! *** Weighted average weir/orifice 
          RATIO = (HUD-HB)/(0.25*HB)
          Q1 = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))      ! *** Submerged orifice flow
          Q2 = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))     ! *** Sub-critical weir flow
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        else
          ! *** Sub-critical weir flow
          QHS = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))    ! *** Sub-critical weir flow
          
        endif
             
      endif
      
    elseif( HYD_STR(NCTL).NQCTYP == 7  )then 
      ! *** WEIR
      ! *** FORMULA FROM http://www.codecogs.com/library/engineering/fluid_mechanics/weirs/discharge.php
      HUD = ZHU - USINV
      HDD = ZHD - USINV
      DIA = HUD
      CDD = HS_COEFF(IHYD,2)
      if( HS_XSTYPE(IHYD) == 5  )then
        ! *** Rectangular weir       
        QHS = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5   ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        
      elseif( HS_XSTYPE(IHYD) == 6  )then
        ! *** Triangle notch
        QHS = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
        
      elseif( HS_XSTYPE(IHYD) == 7  )then
        ! *** Trapezoidal 
        QHS1 = 2.0D0/3.0D0*CDD*HS_WIDTH(IHYD)*SQRT(2.0*G)*HUD**1.5   ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        QHS2 = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
        QHS  = QHS1 + QHS2

      elseif( HS_XSTYPE(IHYD) == 8  )then
        ! *** Broad crested weir (CLH**1.5)
        QHS = CDD*WOPEN*HUD**1.5        ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        
      else
        PRINT '(A29,2I4)', ' *** BAD CROSS-SECTION TYPE:',IHYD, HS_XSTYPE(IHYD)
        call STOPP('.')
        
      endif
      
      ! *** Account for submergence (Villemonte_1947_Submerged_weir_discharge_studies)
      if( ZHD > USINV )then
        QHS = QHS*(1.-HDD/HUD)**0.385
      endif
      
    elseif( HYD_STR(NCTL).NQCTYP == 8  )then 
      ! *** ORIFICE
      HDD = ZHD - USINV
      W  = WOPEN               ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      HB = HOPEN                ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      CDD = HS_COEFF(IHYD,2)
      
      call CROSS_SECTION(IHYD,HUD,HOPEN,WOPEN,SOPEN,DIA,FAREA,WETPER,HRAD,TWID)       

      if( DIA <= 1.E-6 ) return

      ! *** Comparison of HUD & HDD
      if( HDD >= DIA  )then
        ! *** Submerged orifice
        QHS = CDD*FAREA*(2.*G*(HUD-HDD))**0.5
        
      elseif( HUD > DIA .and. HDD < DIA )then
        ! *** Free flowing orifice jet
        HUD = ZHU - USINV + 0.5*HB  ! *** HUD based on centerline
        QHS = CDD*FAREA*(2.*G*HUD)**0.5
        
        ! *** Account for submergence (Villemonte_1947_Submerged_weir_discharge_studies)
        if( ZHD > USINV )then
          QHS = QHS*(1.-HDD/HUD)**0.385
        endif

      else
        ! *** HUD < DIA.  TREAT AS WEIR FLOW
        CDD = HS_COEFF(IHYD,1)
        if( HS_XSTYPE(IHYD) == 5  )then
          ! *** Rectangular weir       
          QHS = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5                      ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        
        elseif( HS_XSTYPE(IHYD) == 6  )then
          ! *** Triangle notch
          QHS = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
        
        elseif( HS_XSTYPE(IHYD) == 7  )then
          ! *** Trapezoidal 
          QHS1 = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5                     ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
          QHS2 = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
          QHS  = QHS1 + QHS2

        else
          ! *** Catch all weir shapes (circles, half-circles, etc).  base on area.
          ! ***  CDD must include all geometric effects
          QHS = CDD*FAREA*SQRT(2.0*G)*HUD**1.5
        endif
      
        ! *** Account for submergence (Villemonte_1947_Submerged_weir_discharge_studies)
        if( ZHD > USINV )then
          QHS = QHS*(1.-HDD/HUD)**0.385
        endif
      endif

    elseif( HYD_STR(NCTL).NQCTYP == 9  )then 
      ! *** Navigation Lock.  No hydraulic structure flow, just a toggle for masks

      TOPEN1 = HOPEN  + WOPEN  + SOPEN                ! *** Current gate positions
        
      ! *** Lock filling/draining flows and timing (only apply after first iteration due to initializations)
      if( HYD_STR(NCTL).TRANSIT > 0.0 .and. NITER > 1 )then                   
        
        ! *** Get gate positions at the end of the transit period
        call GATE_OPENING_INTERP(QCTLSER(ID), TIMESEC+HYD_STR(NCTL).TRANSIT, HOPEN2, WOPEN2, SOPEN2, IBLOCK)
        
        TOPEN2 = HOPEN2 + WOPEN2 + SOPEN2             ! *** Total gate positions at the end of the transit period
        
        ! *** Check if gate openings will open in the next transit time
        if( TOPEN2 > TOPEN1+0.001 )then
          
          ! *** Find the current upstream reference elevation
          HYD_STR(NCTL).ICLOSEDUS = 0    ! *** Default to all open gates
          M1 = NCTL
          M2 = NCTL
          do while (M2 > 0)
            M2 = HYD_STR(M2).IUSBC     ! *** Next upstream cell
            if( M2 == 0 .or. HYD_STR(M1).ISTATE == 0 )then
              HYD_STR(NCTL).LUR  = LIJ(HYD_STR(M1).IQCTLU,HYD_STR(M1).JQCTLU)
              HYD_STR(NCTL).ZHUR = BELV(HYD_STR(NCTL).LUR) + HP(HYD_STR(NCTL).LUR)
              if( M2 > 0 ) HYD_STR(NCTL).ICLOSEDUS = 1    ! *** Found a closed gate
              EXIT
            endif
            M1 = M2                     ! *** Previous upstream cell (always valid)
          enddo        
          
          ! *** Find the current downstream reference elevation
          HYD_STR(NCTL).ICLOSEDDS = 0    ! *** Default to all open gates
          M1 = NCTL
          M2 = NCTL
          do while (M2 > 0)
            M2 = HYD_STR(M2).IDSBC     ! *** Next downstream cell
            if( M2 == 0 .or. HYD_STR(M1).ISTATE == 0 )then
              HYD_STR(NCTL).LDR  = LIJ(HYD_STR(M1).IQCTLD,HYD_STR(M1).JQCTLD)
              HYD_STR(NCTL).ZHDR = BELV(HYD_STR(NCTL).LDR) + HP(HYD_STR(NCTL).LDR)
              if( M2 > 0 ) HYD_STR(NCTL).ICLOSEDDS = 1    ! *** Found a closed gate
              EXIT
            endif
            M1 = M2                     ! *** Previous downstream cell (always valid)
          enddo  
          
          ! *** Turn on flows 
          do NWR = 1,NQWR
            if( HYD_STR(NCTL).NWRGRP == WITH_RET(NWR).GROUPID )then
              WITH_RET(NWR).QWR = WITH_RET(NWR).QWR + HYD_STR(NCTL).QNWR*0.1   ! *** Rampup the flows
              if( WITH_RET(NWR).QWR >= HYD_STR(NCTL).QNWR )then
                WITH_RET(NWR).QWR = HYD_STR(NCTL).QNWR
                if( HYD_STR(NCTL).ICLOSEDDS == 1 )then
                  HYD_STR(NCTL).ISTATE = 1        ! *** Rising/Filling
                else   !IF( ZHU <= (ZHD+0.001) )then
                  HYD_STR(NCTL).ISTATE = 2        ! *** Lowering/Draining
                endif
              endif
            endif
          enddo

          ! *** Check if need to close valves
          IFLAG = 0
          if( HYD_STR(NCTL).ISTATE == 1 )then         ! *** Rising/Filling
            ZHUR = HYD_STR(NCTL).ZHUR
            if( ZHU >= (ZHUR-0.001) )then
              IFLAG = 1
            endif
          endif
          if( HYD_STR(NCTL).ISTATE == 2 )then         ! *** Lowering/Draining
            if( ZHU <= (ZHD+0.001) )then
              IFLAG = 1
            endif
          endif
          
          ! *** Zero flows for all upstream lock chambers
          if( IFLAG == 1 )then
            M1 = NCTL
            do while (M1 > 0)
              do NWR = 1,NQWR
                if( HYD_STR(M1).NWRGRP == WITH_RET(NWR).GROUPID )then
                  WITH_RET(NWR).QWR = 0.0             ! *** Turn off flows
                endif
              enddo
              M1 = HYD_STR(M1).IUSBC
            enddo
          endif
        
        endif
       
      endif
      
      QHS = 0.0
      HSCTL(NCTL).CUR.FLOW = QHS
      
      ! *** Check if navigation lock gates are fully open or closed
      TOPEN1 = HYD_STR(NCTL).CURHEI + HYD_STR(NCTL).CURWID + HYD_STR(NCTL).CURSIL   ! *** Total gate positions: Previous time
      TOPEN2 = HOPEN + WOPEN + SOPEN                                                ! *** Total gate positions: Current time
      
      if( ABS(TOPEN2 - TOPEN1) < 0.001 )then
        return   ! *** Lock status has not changed.  Do nothing.  
      elseif( TOPEN2 > TOPEN1 )then
        ! *** Lock gates are CLOSED, now OPEN the gates and activate flow between LU and LD
        FOPEN   = 1.0    ! *** Remove mask and activate flow
        PRINT '(" Lock Gate: ",I3, " @ ", F12.4," opened ",2F8.2,2F10.3)', NCTL, TIMEDAY, TOPEN1, TOPEN2, ZHU, ZHD
        HYD_STR(NCTL).ISTATE = 3        ! *** Fully open
        
        ! *** Ensure flows are turned off
        M1 = NCTL
        do while (M1 > 0)
          do NWR = 1,NQWR
            if( HYD_STR(M1).NWRGRP == WITH_RET(NWR).GROUPID )then
              WITH_RET(NWR).QWR = 0.0             ! *** Turn off flows
            endif
          enddo
          M1 = HYD_STR(M1).IUSBC
        enddo
      else
        ! *** Lock gates are OPEN, now CLOSE the gates and deactivate any flow between LU and LD
        FOPEN   = 0.0    ! *** Add mask and deactivate flow
        PRINT '(" Lock Gate: ",I3, " @ ", F12.4," closed ",2F8.2,2F10.3)', NCTL, TIMEDAY, TOPEN1, TOPEN2, ZHU, ZHD
        HYD_STR(NCTL).ISTATE = 0        ! *** Fully Closed
      endif
      
      ! *** Find the appropriate face and set the corresponding computational flags
      if( LD == LSC(LU) )then
        SVB(LU)  = FOPEN
        SVBO(LU) = FOPEN
        SAAY(LU) = FOPEN
        SVB3D(LU,KSZ(LU):KC)  = FOPEN
        SVB3DO(LU,KSZ(LU):KC) = FOPEN
      elseif( LU == LSC(LD) )then
        SVB(LD)  = FOPEN
        SVBO(LD) = FOPEN
        SAAY(LD) = FOPEN
        SVB3D(LD,KSZ(LD):KC)  = FOPEN
        SVB3DO(LD,KSZ(LD):KC) = FOPEN
      elseif( LD == LWC(LU) )then
        SUB(LU)  = FOPEN
        SUBO(LU) = FOPEN
        SAAX(LU) = FOPEN
        SUB3D(LU,KSZ(LU):KC)  = FOPEN
        SUB3DO(LU,KSZ(LU):KC) = FOPEN
      elseif( LU == LWC(LD) )then
        SUB(LD)  = FOPEN
        SUBO(LD) = FOPEN
        SAAX(LD) = FOPEN
        SUB3D(LD,KSZ(LD):KC)  = FOPEN
        SUB3DO(LD,KSZ(LD):KC) = FOPEN
      endif
      HYD_STR(NCTL).CURHEI = HOPEN     ! *** Set the curent time for the height (m)
      HYD_STR(NCTL).CURWID = WOPEN     ! *** Set the curent time for the width (m)
      HYD_STR(NCTL).CURSIL = SOPEN     ! *** Set the curent time for the sill height (m)

      return   ! *** No need to compute the remaining flow dependent variables
      
    elseif( HYD_STR(NCTL).NQCTYP == 10  )then 
      ! *** Floating skimmer wall
      call STOPP('*** HY-STRUCTURE IS NOT SUPPORTED: 9')
      
    elseif( HYD_STR(NCTL).NQCTYP == 11  )then 
      ! *** Submerged weir
      QHS = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5                        ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING

    endif 
    HSCTL(NCTL).CUR.FLOW = QHS
    
1000 continue    
    ! *** Horizontal distribution
    QHS = QHS*HYD_STR(NCTL).HS_FACTOR
    
    ! *** SET THE CURRENT HYDRAULIC STRUCTURE FLOWS BY LAYER
    ! *** Use LAYER ELEVATIONS AND INVERTS TO SELECT WHICH LAYERS TO WITHDRAW FROM
    if( KC == 1 )then
      HSZ(KC) = 1.0
      CUMZ = 1.
    else
      HSZ = 0.
      CUMZ = 0.
      TOPZ = MIN(DIA+USINV,ZHU)
      do K = 1,KC
        TOLZ = BELV(LU) + Z(LU,K)*HP(LU)
        if( TOLZ > USINV )then
          DELZ = MIN( (TOLZ-USINV)/DIA, 1.0) - CUMZ
          HSZ(K) = DELZ
          CUMZ = CUMZ + DELZ
          if( CUMZ >= 1.0 .or. TOLZ > TOPZ )EXIT
        endif
      enddo
    endif
    
    ! *** ADD FLOWS TO GLOBAL VARIABLE QSUM
    do K = 1,KC
      QCTLT(K,NCTL,1) = QHS*HSZ(K)/CUMZ
      QSUM(LU,K) = QSUM(LU,K) - QCTLT(K,NCTL,1)      ! *** UPSTREAM
    enddo
    
    ! *** Use LAYER ELEVATIONS AND INVERTS TO SELECT WHICH LAYERS TO RETURN THE FLOW TO
    if( LD /= LC )then
      ! *** DOWNSTREAM ACTIVE
      HSZ = 0.
      CUMZ = 0.
      TOPZ = MIN(DIA+DSINV,ZHD)
      if( (DSINV+0.0001) >= TOPZ .or. KC == 1 )then
        ! *** PLACE ALL OF THE FLOW IN THE TOP LAYER
        HSZ(KC) = 1.0
        CUMZ = 1. 
      else
        do K = 1,KC
          TOLZ = BELV(LD) + Z(LD,K)*HP(LD)
          if( TOLZ > DSINV )then
            DELZ = MIN( (TOLZ-DSINV)/DIA, 1.0) - CUMZ
            HSZ(K) = DELZ
            CUMZ = CUMZ + DELZ
            if( CUMZ >= 1.0 .or. TOLZ > TOPZ )EXIT
          endif
        enddo
      endif
      if( ABS(CUMZ - 1.0) > 0.01 )then
        !CUMZ = 1.                        ! DELME
      endif

      ! *** ADD RETURN FLOWS TO GLOBAL VARIABLE QSUM
      do K = 1,KC
        QCTLT(K,NCTL,2) = QHS*HSZ(K)/CUMZ
        QSUM(LD,K) = QSUM(LD,K) + QCTLT(K,NCTL,2)  ! *** DOWNSTREAM
      enddo
    endif
    
    
  END SUBROUTINE COMPUTE_HSFLOW
  
  SUBROUTINE CROSS_SECTION(IHYD,HUD,HOPEN,WOPEN,SOPEN,DIA,FAREA,WETPER,HRAD,TWID,KCON)
    ! *** CALCULATE HYD. STRUCTURE parameterS:
    ! *** FAREA  : FLOW AREA            [M2]
    ! *** WETPER : WET PERIMETER        [M]
    ! *** HRAD   : HYDRAULIC RADIUS     [M]
    ! *** KCON   : CONVEYANCE           [M3/S]
    ! *** TWID   : TOP WIDTH            [M]
    ! *** DIA    : DIAMETER OR HEIGHT   [M]
    ! *** HOPEN  : OPENING HEIGHT       [M]
    ! *** WOPEN  : OPENING WIDTH        [M]
    ! *** SOPEN  : CHANGE IN SILL LEVEL [M]
    integer(4),parameter   :: NMAX = 100
    integer(4),intent(IN ) :: IHYD
    real,      intent(IN ) :: HUD,HOPEN,WOPEN,SOPEN
    real,      intent(OUT) :: FAREA,WETPER,HRAD,TWID,DIA
    real,OPTIONAL,intent(OUT) :: KCON
    real       :: RA,RB,XH,DEL,ALP,WETPER4,YH,A
    integer(4) :: NN
    real       :: X(NMAX),Y(NMAX)
  
    if( HS_XSTYPE(IHYD) == 1 .or.  HS_XSTYPE(IHYD) == 2  )then
      ! *** 1 - CIRCLE: RA = HS_WIDTH, RB IS NOT USED
      ! *** 2 - HALF CIRCLE: RA = HS_WIDTH, RB IS NOT USED
      RA = HS_WIDTH(IHYD)/2.
      DIA = 2.*RA
      if( HS_XSTYPE(IHYD) == 2 )then
        YH = MIN(HUD+RA,DIA)
      else
        YH = MIN(HUD,DIA)
      endif
      if( YH >= DIA  )then
        ! *** PIPE IS FULL
        FAREA = PI*RA**2
        WETPER = 2*PI*RA
        TWID = 0
      else
        ! *** PARTIALLY FILLED PIPE
        DEL = ABS(RA-YH)
        ALP = ACOS(DEL/RA)
        TWID = 2*SQRT(RA**2-DEL**2)
        FAREA = ALP*RA**2 - 0.5*DEL*TWID
        WETPER = 2*ALP*RA
        if( YH > RA )then
          FAREA = PI*RA**2 - FAREA
          WETPER = 2*PI*RA-WETPER
        endif
      endif
      if( HS_XSTYPE(IHYD) == 2 )then
        ! *** HALF PIPE CORRECTION: YH > RA
        FAREA  = FAREA - 0.5*PI*RA**2
        WETPER = WETPER - PI*RA + 2.*RA
        DIA = RA
      endif
      
    elseif( HS_XSTYPE(IHYD) == 3 .or.  HS_XSTYPE(IHYD) == 4  )then
      ! *** ELLIP: SEMI-AXIS RA = HS_WIDTH/2, RB = HS_HEIGHT/2
      RA = HS_WIDTH(IHYD)/2.
      DIA = HS_HEIGHT(IHYD)

      if( HS_XSTYPE(IHYD) == 4 )then
        RB = HS_HEIGHT(IHYD)
        YH = MIN(HUD+RB,DIA)
      else
        RB = HS_HEIGHT(IHYD)/2.
        YH = MIN(HUD,DIA)
      endif

      X = (/ ((NN-1)*RA/(NMAX-1),NN = 1,NMAX) /) 
      Y  = RB*SQRT(1.0-X**2/RA**2)
      WETPER4 = 0
      do NN = 2,NMAX
        WETPER4 = WETPER4+SQRT((X(NN)-X(NN-1))**2+(Y(NN)-Y(NN-1))**2) 
      enddo
      
      if( YH >= 2.*RB )then
        ! *** FULL FLOW
        FAREA = PI*RA*RB
        WETPER = 4.*WETPER4
        TWID = 0
      else
        ! *** PARTIAL FLOW
        DEL = ABS(RB-YH)
        XH = RA*SQRT(1.0-DEL**2/RB**2) 
        TWID = 2.*XH

        ! *** YH < RB        
        WETPER = 0
        do NN = 2,NMAX
          if( X(NN) <= XH )then
            WETPER = WETPER+SQRT((X(NN)-X(NN-1))**2+(Y(NN)-Y(NN-1))**2)
          else
            EXIT
          endif
        enddo
        WETPER = 2*WETPER
    
        FAREA = 0
        do NN = 2,NMAX
          if( X(NN) <= XH )then
            FAREA = FAREA + 0.5*(X(NN)-X(NN-1))*(Y(NN)+Y(NN-1))
          else
            EXIT
          endif
        enddo
        FAREA = 2.*( FAREA - DEL*XH )   
    
        if( YH > RB )then
          WETPER = 4.*WETPER4 - WETPER
          FAREA  = PI*RA*RB - FAREA
        endif
      endif
    
      if( HS_XSTYPE(IHYD) == 4 )then
        ! *** HALF PIPE CORRECTION     
        FAREA  = FAREA - 0.5*PI*RA*RB
        WETPER = WETPER - 2.*WETPER4 + 2.*RA
      endif
      
    elseif( HS_XSTYPE(IHYD) == 5  )then
      ! *** RECTANGLE: RA = HS_WIDTH, RB = HS_HEIGHT
      RA = WOPEN                                    ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      RB = HOPEN                                    ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      DIA = RB
      YH = MIN(HUD,RB)
      TWID = RA
      FAREA  = TWID*YH
      WETPER = TWID + 2.*YH   
      if( HUD > RB ) WETPER = WETPER + RB
      
    elseif( HS_XSTYPE(IHYD) == 6  )then
      ! *** V NOTCH: RB = HS_HEIGHT, HS_ANGLE
      DIA   = HS_HEIGHT(IHYD)
      RA    = DIA*TAN(0.5*HS_ANGLE(IHYD))
      TWID  = 2.*RA
      FAREA = DIA*RA
      WETPER= 2.*SQRT(RA**2+DIA**2)
      
    elseif( HS_XSTYPE(IHYD) == 7  )then
      ! *** TRAPEZOID: RA = HS_WIDTH, RB = HS_HEIGHT
      DIA   = HS_HEIGHT(IHYD)      
      RA    = DIA*TAN(0.5*HS_ANGLE(IHYD))
      TWID  = HS_WIDTH(IHYD) + 2.*RA    
      FAREA = DIA*HS_WIDTH(IHYD) + DIA*RA
      WETPER= 2.*SQRT(RA**2+DIA**2) + HS_WIDTH(IHYD)
      
    elseif( HS_XSTYPE(IHYD) == 8  )then
      ! *** PARABOLA: RA = HS_WIDTH, RB = HS_HEIGHT
      RA = HS_WIDTH(IHYD)/2.
      RB = HS_HEIGHT(IHYD)/2.
      DIA = 2.*RB
      YH = MIN(HUD,2.*RB)
      A  = 2*RB/RA**2
      X = (/ ((NN-1)*RA/(NMAX-1),NN = 1,NMAX) /) 
      Y  = A*X**2
      XH = SQRT(YH/A) 
      WETPER = 0
      do NN = 2,NMAX
        if( X(NN) <= XH )then
          WETPER = WETPER+SQRT((X(NN)-X(NN-1))**2+(Y(NN)-Y(NN-1))**2)
        else
          EXIT
        endif
      enddo
      TWID = 2*XH
      WETPER = 2*WETPER  
      FAREA = 4.0*XH*YH/3.0    
    endif
    
    HRAD = FAREA/WETPER                          ! *** HYDRAULIC RADIUS
    if( PRESENT(KCON) )then
      if( HS_MANN(IHYD) > 0 )then
        KCON = FAREA/HS_MANN(IHYD)*HRAD**0.6667  ! *** CONVEYANCE
      else
        KCON = 0.
      endif
    endif
  
  END SUBROUTINE CROSS_SECTION
  
  SUBROUTINE SWAP_USDS(NCTL,LU,LD,ZHU,ZHD,USINV,DSINV)
    ! *** BASED ON WATER SURFACE ELEVATIONS, SWAP UPSTREAM AND DOWNSTREAM CELLS FOR FLOW AND TRANSPORT
    integer, intent(IN ) :: NCTL
    integer, intent(OUT) :: LU,LD
    real, intent(OUT)   :: ZHU,ZHD
    real, intent(INOUT) :: USINV,DSINV
    integer(4) :: ITMP,JTMP
    real       :: TMP
    
    ! *** SWAP
    ITMP = HYD_STR(NCTL).IQCTLU
    JTMP = HYD_STR(NCTL).JQCTLU
    HYD_STR(NCTL).IQCTLU = HYD_STR(NCTL).IQCTLD
    HYD_STR(NCTL).JQCTLU = HYD_STR(NCTL).JQCTLD
    HYD_STR(NCTL).IQCTLD = ITMP
    HYD_STR(NCTL).JQCTLD = JTMP
    
    TMP   = USINV
    USINV = DSINV
    DSINV = TMP
    
    LU = LIJ(HYD_STR(NCTL).IQCTLU,HYD_STR(NCTL).JQCTLU)
    LD = LIJ(HYD_STR(NCTL).IQCTLD,HYD_STR(NCTL).JQCTLD)
    ZHU = BELV(LU) + HP(LU)
    ZHD = BELV(LD) + HP(LD)
  
  END SUBROUTINE SWAP_USDS

  SUBROUTINE HYDSTRUCT_CHECK(NX,L)

    integer, intent(IN)::L,NX
    integer :: NBAD
    
    if( HS_XSTYPE(L) < 0 .or. HS_XSTYPE(L) > 8 )then
      PRINT*,' *** BAD CROSS-SECTION TYPE, EQ#, HS_XSTYPE = ',L,HS_XSTYPE(L)
      call STOPP('.')
    endif   
    NBAD = 0
    if( NX == 5 .or. NX == 8 )then
      SELECT CASE(HS_XSTYPE(L))
      CASE (1,2)
        if(HS_WIDTH(L) <= 1E-6 )then
          ! *** 1 - CIRCLE: RA = HS_WIDTH, RB IS NOT USED
          ! *** 2 - HALF CIRCLE: RA = HS_WIDTH, RB IS NOT USED
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          call STOPP('.') 
        endif
        
      CASE (3:5)
        ! *** ELLIP/RECTANGLE/PARABOL: SEMI-AXIS RA = HS_WIDTH/2, RB = HS_HEIGHT/2
        if( HS_WIDTH(L) <= 1E-6  )then
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          NBAD = 1
        endif
        if( HS_HEIGHT(L) <= 1E-6 )then
          PRINT*,' *** BAD HS_HEIGHT: EQ#, HS_HEIGHT = ',L,HS_HEIGHT(L)
          NBAD = 1
        endif 
        if( NBAD > 0 ) CALL STOPP('.')
        
      CASE (6)
        ! *** V NOTCH: RB = HS_HEIGHT, HS_ANGLE
        if( HS_HEIGHT(L) <= 1E-6  )then
          PRINT*,' *** BAD HS_HEIGHT: EQ#, HS_HEIGHT = ',L,HS_HEIGHT(L)
          NBAD = 1
        endif
        if( HS_ANGLE(L) <= 1E-6 )then
          PRINT*,' *** BAD HS_ANGLE: EQ#, HS_ANGLE = ',L,HS_ANGLE(L)
          NBAD = 1
        endif   
        if( NBAD > 0 ) CALL STOPP('.')
          
      CASE (7)
        ! *** TRAPEZOID: RA = HS_WIDTH, RB = HS_HEIGHT, HS_ANGLE
        if( HS_WIDTH(L) <= 1E-6  )then
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          NBAD = 1
        endif
        if( HS_HEIGHT(L) <= 1E-6  )then
          PRINT*,' *** BAD HS_HEIGHT: EQ#, HS_HEIGHT = ',L,HS_HEIGHT(L)
          NBAD = 1
        endif
        if( HS_ANGLE(L) <= 1E-6 )then
          PRINT*,' *** BAD HS_ANGLE: EQ#, HS_ANGLE = ',L,HS_ANGLE(L)
          NBAD = 1
        endif     
        if( NBAD > 0 ) CALL STOPP('.')

      CASE (8)
        ! *** BROAD CRESTED WEIR: RA = HS_WIDTH
        if( HS_WIDTH(L) <= 1E-6  )then
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          call STOPP('.')
        endif
      END SELECT
    endif
        
  END SUBROUTINE HYDSTRUCT_CHECK
  
  SUBROUTINE PUMP_OPERATION_RULES(CTL, HVAL, CVAL, RVAL)
    ! *** CTL - STRUCTURE FOR PUMP
    ! *** HVAL - WATER LEVEL (HEAD) OF TRIGGER
    ! *** CVAL - OLD FLOWS
    ! *** RVAL - NEW FLOWS
  
    type(STRCONTROL), intent(INOUT) :: CTL
    real, intent(IN) :: HVAL, CVAL
    real, intent(OUT) :: RVAL
    real(RKD) :: TSUM
    integer :: K, ID

    ID = CTL.ID
    CTL.SUBID = 0
    TSUM = CTL.TOPEN + DELT
    
    if( CTL.CUR.STATE == 2 )then
      ! *** GATE OPENING
      RVAL = CVAL + CTL.CUR.RATE * DELT/60.
      CTL.TOPEN = TSUM
      if( RVAL >= CTL.LIM.FLOW )then
        RVAL = CTL.LIM.FLOW
        CTL.CUR.STATE = 1                   ! *** FULLY OPEN OR MAXIMUM FLOW
        CTL.CUR.RATE = 0.
      endif
    elseif( CTL.CUR.STATE == 3 )then
      ! *** GATE CLOSING
      RVAL = CVAL - CTL.CUR.RATE * DELT/60.
      CTL.TOPEN = TSUM
      if( RVAL <= CTL.LIM.FLOW )then
        RVAL = CTL.LIM.FLOW
        CTL.CUR.STATE = 0                   ! *** FULLY CLOSED OR MINIMUM FLOW
        CTL.CUR.RATE = 0.
        CTL.TOPEN = 0.
      endif
    else
      ! *** FULLY OPEN OR CLOSED
      RVAL = CVAL
      
      ! *** CHECK FOR START OF OPENING
      do K = 1,QCTRULES(ID).NTRGON
        if( (HVAL >= QCTRULES(ID).TRGON(K).LEVEL) .and. &
          ((CTL.CUR.STATE == 0) .or. (QCTRULES(ID).TRGON(K).LEVEL > CTL.CUR.LEVEL)) )then
          ! *** WATER LEVELS TRIGGER OPENING
          CTL.CUR.STATE = 2
          CTL.CUR.LEVEL = QCTRULES(ID).TRGON(K).LEVEL
          CTL.CUR.RATE = QCTRULES(ID).TRGON(K).RATE
          CTL.LIM.FLOW = QCTRULES(ID).TRGON(K).FLOW
          CTL.LIM.SILL = QCTRULES(ID).TRGON(K).SILL
          CTL.LIM.WIDTH = QCTRULES(ID).TRGON(K).WIDTH
          CTL.LIM.HEIGHT = QCTRULES(ID).TRGON(K).HEIGHT
          CTL.TSTART = TIMESEC
          TSUM = 0.
          CTL.SUBID = K
          EXIT
        endif
      enddo

      ! *** CHECK FOR START OF CLOSING
      do K = QCTRULES(ID).NTROFF,1,-1
        if( (HVAL <= QCTRULES(ID).TROFF(K).LEVEL) .and. &
          ! *** WATER LEVELS TRIGGER CLOSING
          ((CTL.CUR.STATE == 1) .or. (QCTRULES(ID).TROFF(K).LEVEL < CTL.CUR.LEVEL)) )then
          CTL.CUR.STATE = 3
          CTL.CUR.LEVEL = QCTRULES(ID).TROFF(K).LEVEL
          CTL.CUR.RATE = QCTRULES(ID).TROFF(K).RATE
          CTL.LIM.FLOW = QCTRULES(ID).TROFF(K).FLOW
          CTL.LIM.SILL = QCTRULES(ID).TROFF(K).SILL
          CTL.LIM.WIDTH = QCTRULES(ID).TROFF(K).WIDTH
          CTL.LIM.HEIGHT = QCTRULES(ID).TROFF(K).HEIGHT
          CTL.TSTART = TIMESEC
          TSUM = 0.
          CTL.SUBID = K
          EXIT
        endif
      enddo
      CTL.TOPEN = TSUM   
    endif
  END SUBROUTINE PUMP_OPERATION_RULES
  
  SUBROUTINE GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
    ! *** NCTL - INDEX OF STRUCTURE FOR GATE
    ! *** HVAL - WATER LEVEL (HEAD) OF TRIGGER
    ! *** CVAL - OLD LEVEL
    ! *** LVAL - TRIGGER LEVEL
    ! *** RVAL - NEW LEVEL

    integer, intent(IN) :: NCTL
    real, intent(IN) :: HVAL, CVAL, LVAL
    real, intent(OUT) :: RVAL
    real(RKD) :: TSUM
    integer :: K, ID

    ID = HSCTL(NCTL).ID
    HSCTL(NCTL).SUBID = 0
    
    TSUM = HSCTL(NCTL).TOPEN + DELT
    if( HSCTL(NCTL).CUR.STATE == 2 )then
      ! *** GATE OPENING
      RVAL = CVAL + HSCTL(NCTL).CUR.RATE * DELT/60.
      HSCTL(NCTL).TOPEN = TSUM
      if( RVAL >= LVAL )then
        RVAL = LVAL
        HSCTL(NCTL).CUR.STATE = 1           ! *** FULLY OPEN OR MAXIMUM OPENING
        HSCTL(NCTL).CUR.RATE = 0.
        
        if( .TRUE. ) WRITE(6, '("GATE ",I3," IS FULLY OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        write(mpi_efdc_out_unit, '("GATE ",I3," IS FULLY OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        close(mpi_efdc_out_unit)
      endif
      
    elseif( HSCTL(NCTL).CUR.STATE == 3 )then
      ! *** GATE CLOSING
      RVAL = CVAL - HSCTL(NCTL).CUR.RATE * DELT/60.
      HSCTL(NCTL).TOPEN = TSUM
      if( RVAL <= LVAL )then
        RVAL = LVAL
        HSCTL(NCTL).CUR.STATE = 0           ! *** FULLY CLOSED OR MINIMUM CLOSURE
        HSCTL(NCTL).CUR.RATE = 0.
        HSCTL(NCTL).TOPEN = 0.

        if( .TRUE. ) WRITE(6, '("GATE ",I3," IS FULLY CLOSED (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        write(mpi_efdc_out_unit, '("GATE ",I3," IS FULLY CLOSED (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        close(mpi_efdc_out_unit)
      endif
      
    else
      ! *** FULLY OPEN OR CLOSED
      RVAL = CVAL                           ! *** SET INITIALLY OPENING TO FINAL, I.E. NO CHANGE IN OPENING
      
      ! *** IF CLOSED, CHECK FOR START OF OPENING
      do K = 1,QCTRULES(ID).NTRGON
        if( HVAL >= QCTRULES(ID).TRGON(K).LEVEL )then
          if( HSCTL(NCTL).CUR.STATE == 0 .or. QCTRULES(ID).TRGON(K).LEVEL > HSCTL(NCTL).CUR.LEVEL )then
            ! *** WATER LEVELS TRIGGER OPENING
            HSCTL(NCTL).CUR.STATE  = 2
            HSCTL(NCTL).CUR.LEVEL  = QCTRULES(ID).TRGON(K).LEVEL
            HSCTL(NCTL).CUR.RATE   = QCTRULES(ID).TRGON(K).RATE
            HSCTL(NCTL).LIM.FLOW   = QCTRULES(ID).TRGON(K).FLOW
            HSCTL(NCTL).LIM.SILL   = QCTRULES(ID).TRGON(K).SILL
            HSCTL(NCTL).LIM.WIDTH  = QCTRULES(ID).TRGON(K).WIDTH
            HSCTL(NCTL).LIM.HEIGHT = QCTRULES(ID).TRGON(K).HEIGHT
            HSCTL(NCTL).TSTART     = TIMESEC
            HSCTL(NCTL).SUBID      = K
            TSUM = 0.

            if( .TRUE. ) WRITE(6,'("GATE ",I3," STARTING TO OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            write(mpi_efdc_out_unit,'("GATE ",I3," STARTING TO OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            close(mpi_efdc_out_unit)
            
            EXIT
          endif
        endif
      enddo

      ! *** IF FULLY OPEN, CHECK FOR START OF CLOSING
      do K = QCTRULES(ID).NTROFF,1,-1
        if( HVAL <= QCTRULES(ID).TROFF(K).LEVEL )then
          if( HSCTL(NCTL).CUR.STATE == 1 .or. QCTRULES(ID).TROFF(K).LEVEL < HSCTL(NCTL).CUR.LEVEL )then
            ! *** WATER LEVELS TRIGGER CLOSING
            HSCTL(NCTL).CUR.STATE = 3
            HSCTL(NCTL).CUR.LEVEL  = QCTRULES(ID).TROFF(K).LEVEL
            HSCTL(NCTL).CUR.RATE   = QCTRULES(ID).TROFF(K).RATE
            HSCTL(NCTL).LIM.FLOW   = QCTRULES(ID).TROFF(K).FLOW
            HSCTL(NCTL).LIM.SILL   = QCTRULES(ID).TROFF(K).SILL
            HSCTL(NCTL).LIM.WIDTH  = QCTRULES(ID).TROFF(K).WIDTH
            HSCTL(NCTL).LIM.HEIGHT = QCTRULES(ID).TROFF(K).HEIGHT
            HSCTL(NCTL).TSTART     = TIMESEC
            HSCTL(NCTL).SUBID      = K
            TSUM = 0.

            if( .TRUE. ) WRITE(6, '("GATE ",I3," STARTING TO CLOSE (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            write(mpi_efdc_out_unit, '("GATE ",I3," STARTING TO CLOSE (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            close(mpi_efdc_out_unit)
            !PRINT '(20X,3F10.3)', HVAL, LVAL, HSCTL(NCTL).CUR.LEVEL  ! DELME
            EXIT
          endif
        endif
      enddo
      HSCTL(NCTL).TOPEN = TSUM
    endif
    
  END SUBROUTINE GATE_OPERATION_RULES
  
  SUBROUTINE GATE_OPENING_INTERP(TS, T, H, W, S, IBLOCK)
  
    type(CONTROLSERIES), intent(INOUT) :: TS
    integer, intent(IN) :: IBLOCK
    real(RKD), intent(IN) :: T
    real, intent(OUT) :: H, W, S
    real :: FACTOR
    integer :: K

    if( T < TS.TIME(1) )then
      K = 1
    elseif( T > TS.TIME(TS.COUNT) )then
      K = TS.COUNT - 1
    else
      do K = TS.IPOS,TS.COUNT-1
        if( T >= TS.TIME(K) .and. T <= TS.TIME(K+1)) EXIT
      enddo
    endif
    FACTOR = TS.TIME(K+1) - TS.TIME(K)
    
    if( IBLOCK == 1 )then
      FACTOR = 0.0
    elseif( ABS(FACTOR) < 1.0E-9 )then
      FACTOR = 0.5
    else
      FACTOR = (T - TS.TIME(K))/FACTOR
    endif
    
    H = TS.HEIGHT(K)+ FACTOR*(TS.HEIGHT(K+1) - TS.HEIGHT(K))
    W = TS.WIDTH(K) + FACTOR*(TS.WIDTH(K+1)  - TS.WIDTH(K))
    S = TS.SILL(K)  + FACTOR*(TS.SILL(K+1)   - TS.SILL(K))
    if( ABS(T-TIMESEC) < 1.0 ) TS.IPOS = K                              ! *** Updated current index
    
  END SUBROUTINE GATE_OPENING_INTERP
  
  SUBROUTINE CHANGE_RATING_CURVE_BY_RULES(NCTL, IRET)
  
    integer, intent(IN) :: NCTL
    integer, intent(OUT) :: IRET
    real(RKD) :: TSUM, ZHU, ZHD, HVAL
    integer :: LU, LD, K, ID, IVAL

    ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
    LU = LIJ(HSCTL(NCTL).IREFUP, HSCTL(NCTL).JREFUP)
    ZHU = BELV(LU) + HP(LU)
      
    if( HSCTL(NCTL).ITYPE == 3 )then
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
        LD = LIJ(HSCTL(NCTL).IREFDN, HSCTL(NCTL).JREFDN)      
        ZHD = BELV(LD) + HP(LD)
        HVAL = ZHU - ZHD
    else
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
        HVAL = ZHU
    endif
    
    ID = HSCTL(NCTL).ID
    HSCTL(NCTL).SUBID = 0
    TSUM = HSCTL(NCTL).TOPEN + DELT
    if( HVAL > HSCTL(NCTL).PREV )then  ! ASCENDING
      do K = 1,QCTRULES(ID).NTRGON
        if( HVAL >= QCTRULES(ID).TRGON(K).LEVEL )then
          HSCTL(NCTL).CUR.STATE = 2      ! START OPENING
          IRET = QCTRULES(ID).TRGON(K).ID
          HSCTL(NCTL).CUR.ID = IRET
          HSCTL(NCTL).TSTART = TIMESEC
          HSCTL(NCTL).TOPEN = 0.
          HSCTL(NCTL).SUBID = K
          return
        endif
      enddo
    else  ! DESCENDING
      HSCTL(NCTL).TOPEN = TSUM
      do K = QCTRULES(ID).NTROFF,1,-1
        if( HVAL <= QCTRULES(ID).TROFF(K).LEVEL )then
          HSCTL(NCTL).CUR.STATE = 3      ! STOP OPERATIONAL
          IRET = QCTRULES(ID).TROFF(K).ID
          HSCTL(NCTL).CUR.ID = IRET
          HSCTL(NCTL).TSTART = TIMESEC
          HSCTL(NCTL).TOPEN = 0.
          HSCTL(NCTL).SUBID = K
          return
        endif
      enddo
    endif
    HSCTL(NCTL).PREV = HVAL
    
  END SUBROUTINE CHANGE_RATING_CURVE_BY_RULES
  
  SUBROUTINE CHANGE_RATING_CURVE_BY_TIME(TS,T,IRET)
  
    type(CONTROLSERIES), intent(IN) :: TS
    real(RKD), intent(IN) :: T
    integer, intent(OUT) :: IRET
    integer :: K

    if( T < TS.TIME(1) )then
      IRET = TS.ID(1)
    elseif( T > TS.TIME(TS.COUNT) )then
      K = TS.COUNT
      IRET = TS.ID(K)
    else
      do K = 1,TS.COUNT-1
        if( T >= TS.TIME(K) .and. T <= TS.TIME(K+1) )then
          IRET = TS.ID(K)
          EXIT
        endif
      enddo
    endif
  END SUBROUTINE CHANGE_RATING_CURVE_BY_TIME
  
END MODULE
