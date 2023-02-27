! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE HYDSTRUCMOD

  ! *** PURPOSE:  CALCULATE FLOW AT THE HYDRAULIC STRUCTURE BCS
  ! ***           INCLUDING CULVERTS, WEIRS AND SLUICE GATES
  ! *** BASED ON THE PAPER BY NL DILL: ESTUARINE AND COASTAL MODELING (2011)
  ! *** Modeling Hydraulic Control Structures in Estuarine Environments with EFDC
  ! *** DATE: AUGUST 2015
  ! ***
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-12       P M CRAIG         IMPROVED HYDRAULIC STRUCTURES FLOW CALCULATION
  !    2015-11       D H CHUNG         IMPLEMENTED NEW HYDRAULIC STRUCTURES BASED ON DILL
  
  USE GLOBAL 
  USE INFOMOD
  USE Variables_MPI
  
  IMPLICIT NONE
  
  TYPE OPERATIONDATA                !< Trigger for operation of gates/pumps
    REAL :: LEVEL                   !< Trigger level
    INTEGER :: STATE                !< Gate position
                                    !    0 = Completely Closed/Off
                                    !    1 = Fully Open/On
                                    !    2 = Opening
                                    !    3 = Closing
    INTEGER :: ID                   !< Index of rating curve/mask for control parameter
    REAL :: HEIGHT                  !< Opening height (m), for upward opening
    REAL :: WIDTH                   !< Opening width (m), for sideward opening
    REAL :: SILL                    !< Sill level change (m), for downward opening
    REAL :: FLOW                    !< Flow rate (m3/s), for pumps
    !INTEGER :: UNITS                !< Number of gates, pump units  (NOT USED)
    REAL :: RATE                    !< Rate of opening/closing
  END TYPE

  TYPE CONTROLRULES                 !< Operation rules for of gates/pumps
    INTEGER :: PARAM                !< Mask for control parameters
    INTEGER :: STATE                !< Structure State, 0 = OFF, 1 = ON
    INTEGER :: NTRGON,NTROFF  
    TYPE(OPERATIONDATA),ALLOCATABLE :: TRGON(:), TROFF(:)
  END TYPE  
  
  TYPE CONTROLSERIES                !< GATE OPENING TIME-SERIES
    INTEGER :: ITYPE                !< STRUCTURE CONTROL TYPE: 0 = UNCONTROLED, 1=TS, 2=US, 3=HDIFF
    INTEGER :: COUNT                !< # POINTS
    INTEGER :: PARAM                !< MASK FOR CONTROL PARAMETERS
    REAL :: TMUL,TADD
    REAL,ALLOCATABLE :: TIME(:),HEIGHT(:),WIDTH(:),SILL(:),FLOW(:),RATE(:) 
    INTEGER,ALLOCATABLE :: ID(:)
    !INTEGER,ALLOCATABLE :: NUM(:)  ! NUMBER OF GATES/PUMP UNITS (NOT USED)
  END TYPE
  
  TYPE STRCONTROL                   !< INFO ON THE CONTROL OF HYDRAULIC STRUCTURE
    INTEGER :: ITYPE                !< STRUCTURE CONTROL TYPE: 0 = UNCONTROLED, 1=TS, 2=US, 3=HDIFF
    INTEGER :: ID                   !< ID OF ASSOCIATED CONTROL DATA
    INTEGER :: IREFUP,JREFUP        !< UPSTREAM   REFERENCE CELL. ALLOWS REFERENCE CELL BE DIFFERENT FROM US CELL
    INTEGER :: IREFDN,JREFDN        !< DOWNSTREAM REFERENCE CELL. ALLOWS REFERENCE CELL BE DIFFERENT FROM DS CELL
    INTEGER :: SUBID                !< SUB INDEX OF ASSOCIATED CONTROL DATA
    REAL(RKD) :: TSTART, TOPEN
    TYPE(OPERATIONDATA) :: CUR,LIM  !< INITIAL CONDITION (CURRENT VALUE) & LIMIT
    REAL :: PREV                    !< PREVIOUS CONTROL LEVEL
  END TYPE  

  TYPE(STRCONTROL),SAVE,ALLOCATABLE,DIMENSION(:)    :: HSCTL_GL  !< Hydraulic sructure control
  TYPE(STRCONTROL),SAVE,ALLOCATABLE,DIMENSION(:)    :: HSCTL     !< Hydraulic sructure control
  TYPE(CONTROLSERIES),SAVE,ALLOCATABLE,DIMENSION(:) :: QCTLSER   !< Control time-series
  TYPE(CONTROLRULES),SAVE,ALLOCATABLE,DIMENSION(:)  :: QCTRULES  !< Control rules

  INTEGER :: NQCTLSER,NQCTRULES

  TYPE(STRCONTROL),SAVE,ALLOCATABLE,DIMENSION(:)        :: WRCTL_GL       !< Withdrawal/return control (Global)
  TYPE(STRCONTROL),SAVE,ALLOCATABLE,DIMENSION(:)        :: WRCTL          !< Withdrawal/return control (Local)

  INTEGER :: HSOPTION
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: USCELL    ! *** SAVE THE ORIGINAL UPSTREAM CELL
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: DSCELL    ! *** SAVE THE ORIGINAL DOWNSTREAM CELL
    
  CONTAINS 

  SUBROUTINE COMPUTE_HSFLOW(NCTL)
  
    INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: NDIRCULV   ! *** FLOW DIRECTION COUNTER
    INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: NLIMITQ    ! *** FLOW LIMITATION COUNTER
    REAL,SAVE,ALLOCATABLE,DIMENSION(:)    :: CULVQ      ! *** LAST CULVERT FLOW RATE

    INTEGER, INTENT(IN) :: NCTL
    INTEGER             :: LU,LD,K,IHYD,ID,STATE
    
    REAL :: ZHU,ZHD,HUD,HDD,S,SCR,HCR,VCR,VN,QHS,QHS1,QHS2,USINV,DSINV,DIA,CDD
    REAL :: QCR,YN,TOPZ,TOLZ,DELZ,CUMZ,HB,W,RATIO,Q1,Q2,QMAX,DELT
    REAL :: FAREA,WETPER,HRAD,KCON,TWID, HVAL, HOPEN, WOPEN, SOPEN
    REAL :: CVAL, LVAL, RVAL, TSTART, TOPEN
    REAL :: HSZ(KC)
    
    IF( .NOT. ALLOCATED(NDIRCULV) )THEN
      ALLOCATE(NDIRCULV(NQCTL))
      ALLOCATE(NLIMITQ(NQCTL))
      ALLOCATE(CULVQ(NQCTL))
      NDIRCULV = 100
      NLIMITQ = 0
      CULVQ = 0.
    ENDIF

    IF( ISDYNSTP == 0 )THEN
      DELT=DT
    ELSE
      DELT=DTDYN
    END IF

    HSOPTION = 1     ! *** OPTION FOR SLUICE GATE, CHANGE LATER
  
    QCTLT(1:KC,NCTL,1:2) = 0                 ! ** RESET FOR EVERY TIME STEP
    
    IHYD = HYD_STR(NCTL).NQCTLQ ! *** HYDRAULIC STRUCTURE DEFINTION PARAMETERS 

    HOPEN = HS_HEIGHT(IHYD) ! *** DEFAULT OPENING HEIGHT IS THE STRUCTURE DIMENSION
    WOPEN = HS_WIDTH(IHYD)  ! *** DEFAULT OPENING WIDTH IS THE STRUCTURE DIMENSION
    SOPEN = 0.              ! *** DEFAULT IS NO SILL LEVEL CHANGE
    IF( HSCTL(NCTL).ITYPE == 1 )THEN    
      ! *** STRUCTURE IS CONTROLLED BY TIME-SERIES
      ID = HSCTL(NCTL).ID
      CALL GATE_OPENING_INTERP(QCTLSER(ID),TIMESEC,HOPEN, WOPEN, SOPEN)
      
    ELSEIF (HSCTL(NCTL).ITYPE == 2 .OR. HSCTL(NCTL).ITYPE == 3 )THEN
      ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
      LU = LIJ(HSCTL(NCTL).IREFUP, HSCTL(NCTL).JREFUP)
      ZHU = BELV(LU) + HP(LU)
      
      IF( HSCTL(NCTL).ITYPE == 3 )THEN
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
        LD = LIJ(HSCTL(NCTL).IREFDN, HSCTL(NCTL).JREFDN)      
        ZHD = BELV(LD) + HP(LD)
        HVAL = ZHU - ZHD
      ELSE
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
        HVAL = ZHU
      ENDIF
      
      ID = HSCTL(NCTL).ID
      IF( (QCTRULES(ID).PARAM .AND. 16) == 16 )THEN
        ! *** PUMP FLOWS
        CVAL = HSCTL(NCTL).CUR.FLOW
        !LVAL = HSCTL(NCTL).LIM.FLOW
        CALL PUMP_OPERATION_RULES(HSCTL(NCTL), HVAL, CVAL, RVAL)
        QHS = RVAL
        HSCTL(NCTL).CUR.FLOW = RVAL
        LD = LC
        GOTO 1000
      ELSE  
        ! *** FOR GATES, LOOKUP TABLE SHOULD BE HANDLED IN CALQVS.F90
        IF( (QCTRULES(ID).PARAM .AND. 1) == 1 )THEN
          ! *** GATE OPENING UPWARD (RAISING/LOWERING SILL)
          CVAL = HSCTL(NCTL).CUR.HEIGHT
          LVAL = HSCTL(NCTL).LIM.HEIGHT
          CALL GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
          HOPEN = RVAL
          HSCTL(NCTL).CUR.HEIGHT = RVAL
          
        ELSEIF( (QCTRULES(ID).PARAM .AND. 2) == 2 )THEN
          ! *** GATE OPENING SIDEWARD
          CVAL = HSCTL(NCTL).CUR.WIDTH
          LVAL = HSCTL(NCTL).LIM.WIDTH
          CALL GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
          WOPEN = RVAL
          HSCTL(NCTL).CUR.WIDTH = RVAL
          
        ELSEIF( (QCTRULES(ID).PARAM .AND. 4) == 4 )THEN
          ! *** GATE OPENING DOWNWARD (RAISING/LOWERING GATE)
          CVAL = HSCTL(NCTL).CUR.SILL
          LVAL = HSCTL(NCTL).LIM.SILL
          CALL GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
          SOPEN = RVAL
          HSCTL(NCTL).CUR.SILL = RVAL
        ENDIF
      ENDIF

    ELSE
      ! *** STRUCTURE IS UNCONTROLLED, DO NOTHING
    ENDIF
    
    ! *** NOW SET APPROPRIATE HEAD FOR STRUCTURE BASED ON HYDRAULIC STRUCTURE CELLS
    LU = LIJ(HYD_STR(NCTL).IQCTLU,HYD_STR(NCTL).JQCTLU)
    LD = LIJ(HYD_STR(NCTL).IQCTLD,HYD_STR(NCTL).JQCTLD)
    ZHU = BELV(LU) + HP(LU)
    USINV = HS_USELEV(IHYD) + SOPEN         ! 2017-11-06, NTL: ACCOUNTS FOR SILL LEVEL CHANGE
    
    IF( LD == 0 )THEN
      ! *** DOWNSTREAM CELL IS OUTSIDE DOMAIN
      LD = LC
      ZHD = -9999
    ELSE
      ! *** DOWNSTREAM CELL IS ACTIVE, ALLOW REVERSE FLOW
      ZHD = BELV(LD) + HP(LD)
      DSINV = HS_DSELEV(IHYD) + SOPEN     ! 2017-11-06, NTL: ACCOUNTS FOR SILL LEVEL CHANGE
      IF( ZHU < ZHD .AND. HS_REVERSE(IHYD) == 1 )THEN
        IF( NDIRCULV(NCTL) < 6 .AND. HYD_STR(NCTL).NQCTYP == 5 )THEN
          NDIRCULV(NCTL) = NDIRCULV(NCTL)+1
          QHS = 0.5*CULVQ(NCTL)
          CULVQ(NCTL) = QHS
          HDD = ZHD - DSINV                
          CALL CROSS_SECTION(IHYD, HDD, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID) 
          GOTO 1000
        ENDIF
        NDIRCULV(NCTL) = 0
        CALL SWAP_USDS(NCTL,LU,LD,ZHU,ZHD,USINV,DSINV)
      ENDIF
    ENDIF
       
    ! *** DETERMINE IF ELEVATIONS/DEPTH ALLOW ANY FLOW
    HUD = ZHU - USINV
    IF( HUD <= 0. .OR. HP(LU) < HWET .OR. (HS_REVERSE(IHYD) == 0 .AND. ZHD > ZHU) )THEN
      RETURN
    ENDIF    
    
    IF( HYD_STR(NCTL).NQCTYP == 5  )THEN
      ! ** CULVERT 
      ! ** FORMULA FROM N L DILL PAPER:
      
      CALL CROSS_SECTION(IHYD, HUD, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID, KCON)
      
      HCR =  HUD - 0.5*HRAD               !CRITICAL DEPTH
      HCR = MIN(HCR,DIA)
      VCR = SQRT(9.81*HRAD)
      QCR = VCR*FAREA
      SCR = (QCR/KCON)**2   
      
      ! *** DOWNSTREAM HEAD
      IF( LD == LC )THEN
        ! *** FLOW LEAVES DOMAIN (ASSUME CRITICAL DEPTH AT OUTLET)
        HDD = HCR
      ELSE
        HDD = ZHD - DSINV
      ENDIF      

      ! *** PIPE SLOPE
      S = (HS_USELEV(IHYD)-HS_DSELEV(IHYD))/HS_LENGTH(IHYD)
      IF( S < 1.E-6 )THEN
        ! *** PIPE ZERO OR ADVERSE SLOPE
        IF( LD == LC )THEN
          ! *** FLOW LEAVES DOMAIN (ASSUME CRITICAL DEPTH AT OUTLET)
          S = HCR/HS_LENGTH(IHYD)
        ELSE
          S =  (ZHU-ZHD)/HS_LENGTH(IHYD)
        ENDIF      
      ELSE
        IF( LD == LC )THEN
          ! *** FLOW LEAVES DOMAIN
          S = ABS(S)
        ELSE
          S = MIN(S,(ZHU-ZHD)/HS_LENGTH(IHYD))  ! *** MINIMUM OF ENERGY GRADE OR PIPE SLOPE
        ENDIF      
      ENDIF

      IF( HDD > DIA .OR. HUD > 1.5*DIA  )THEN
        ! ** FULL FLOW         
        S = (ZHU-ZHD)/HS_LENGTH(IHYD)
        QHS = KCON*SQRT(S)
        
      ELSEIF( S >= SCR  )THEN
        ! ** SUPERCRITICAL FLOW CONTROL AT INLET LD
        QHS = QCR
          
      ELSEIF( S < SCR .AND. HDD >= HCR  )THEN
        ! ** SUBCRITICAL FLOW TAILWATER CONTROL LD
        DELZ = MIN(HDD,HUD)
        CALL CROSS_SECTION(IHYD, DELZ, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID)
        VN = HRAD**0.6667 * S**0.5 / HS_MANN(IHYD)
        QHS = VN*FAREA     ! *** Q AT HDD
          
      ELSEIF( HDD < HCR .AND. S < SCR )THEN
        ! ** SUBCRITICAL FLOW CONTROL AT OUTLET LU
        YN = HUD - HRAD**1.3333*S/(19.62*HS_MANN(IHYD)**2)
        CALL CROSS_SECTION(IHYD, YN, HOPEN, WOPEN, SOPEN, DIA, FAREA, WETPER, HRAD, TWID)
        VN = (FAREA/WETPER)**0.6667*S**0.5/HS_MANN(IHYD)
        QHS = VN*FAREA
      ENDIF
      
      ! *** LIMIT THE RATE OF GROWTH/REDUCTION (5%)
      IF( CULVQ(NCTL) > 0.01 )THEN
        IF( QHS < 0.95*CULVQ(NCTL) )THEN
          QHS = 0.95*CULVQ(NCTL)
        ELSEIF( QHS > 1.05*CULVQ(NCTL) )THEN
          QHS = 1.05*CULVQ(NCTL)
        ENDIF
      ELSEIF( QHS > 0.015 )THEN
        QHS = 0.0125
      ENDIF
      
      ! *** LIMIT OUTFLOW FROM CULVERT CELL
      QMAX = 0.01*DXYP(LU)/DELT
      IF( HUD < 2.*HDRY ) QMAX = 0.1*DXYP(LU)*HDRY/DELT
      !IF( 10.*DIA > DXYP(LU) ) QMAX = 0.1*QMAX  ! *** LIMIT CONDITION WITH SMALL CELLS RELATIVE TO CULVERT SIZES
      IF( (HP(LU)-H1P(LU)) < -0.01 )THEN
        QMAX = CULVQ(NCTL)/(H1P(LU)-HP(LU))*0.01
      ENDIF
      IF( QHS > QMAX )THEN
        NLIMITQ(NCTL) = NLIMITQ(NCTL) + 1
        IF( NLIMITQ(NCTL) == 1 .OR. MOD(NLIMITQ(NCTL),1000) == 0 )THEN
          OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
          WRITE(mpi_log_unit, "(A,I10,F12.5,I5,2F10.4)" )'CULVERT FLOW LIMITATION (NITER,TIMEDAY,NCTL,QHS,QMAX): ',NITER,TIMEDAY,NCTL,QHS,QMAX
          CLOSE(mpi_log_unit)
        ENDIF
        QHS = QMAX
      ENDIF

      NDIRCULV(NCTL) = NDIRCULV(NCTL)+1
      CULVQ(NCTL) = QHS

    ELSEIF( HYD_STR(NCTL).NQCTYP == 6  )THEN 
      ! *** SLUICE GATE
      HDD = ZHD - USINV
      W  = WOPEN              ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      HB = HOPEN              ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      DIA = MIN(HB,HUD)
      
      IF( HB <= 1E-6 .OR. W <= 1E-6 ) RETURN
      
      ! *** APPROACH FROM N L DILL PAPER:
      ! *** QHS = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5   C=1.45  ! *** SUPER-CRITICAL WEIR FLOW
      ! *** QHS = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))   C=0.61  ! *** SUB-CRITICAL WEIR FLOW
      ! *** QHS = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)          C=0.60  ! *** FREE SLUICE FLOW
      ! *** QHS = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))    C=0.50  ! *** SUBMERGED ORIFICE FLOW
          
      ! ** COMPARISON OF HUD & HDD
      IF( HDD <= 0.64*HUD  )THEN
        ! *** SUPERCRITICAL FLOW
        IF( HUD >= 1.25*HB )THEN
          QHS = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)           ! *** FREE SLUICE FLOW
          
        ELSEIF( HUD > HB .AND. HUD < 1.25*HB )THEN
          ! *** WEIGHTED AVERAGE WEIR/ORIFICE CRITICAL FLOW
          RATIO = (HUD-HB)/(0.25*HB)
          Q1 = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)            ! *** FREE SLUICE FLOW
          Q2 = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5     ! *** SUPER-CRITICAL WEIR FLOW
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        ELSE
          ! *** BROAD CRESTED WEIR FLOW
          QHS = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5    ! *** SUPER-CRITICAL WEIR FLOW
          
        ENDIF
        
      ELSEIF( 0.64*HUD < HDD .AND. HDD < 0.68*HUD )THEN
        ! *** WEIGHTED AVERAGE OF SUPERCRITICAL & SUBCRITICAL
        RATIO = (HDD-0.64*HUD)/(0.04*HUD)
        IF( HUD >= 1.25*HB )THEN
          ! *** WEIGHTED AVERAGE SUPER/SUB CRITICAL ORIFICE FLOW
          Q1 = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)            ! *** FREE SLUICE FLOW
          Q2 = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))      ! *** SUBMERGED ORIFICE FLOW
          QHS = RATIO*Q2 + (1.-RATIO)*Q1
          
        ELSEIF( HUD > HB .AND. HUD < 1.25*HB )THEN
          ! *** WEIGHTED AVERAGE WEIR/ORIFICE 
          Q1 = HS_COEFF(IHYD,3)*W*HB*SQRT(2.*G*HUD)            ! *** FREE SLUICE FLOW
          Q2 = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))     ! *** SUB-CRITICAL WEIR FLOW
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        ELSE
          ! *** WEIGHTED AVERAGE WEIR/ORIFICE 
          Q1 = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))     ! *** SUB-CRITICAL WEIR FLOW
          Q2 = HS_COEFF(IHYD,1)*W*SQRT(G)*(2./3.*HUD)**1.5     ! *** SUPER-CRITICAL WEIR FLOW
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        ENDIF
          
      ELSEIF( HDD >= 0.68*HUD )THEN
        ! *** SUBCRITICAL FLOWS
        IF( HUD >= 1.25*HB )THEN
          ! *** WEIGHTED AVERAGE SUPER/SUB CRITICAL ORIFICE FLOW
          QHS = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))     ! *** SUBMERGED ORIFICE FLOW
          
        ELSEIF( HUD > HB .AND. HUD < 1.25*HB )THEN
          ! *** WEIGHTED AVERAGE WEIR/ORIFICE 
          RATIO = (HUD-HB)/(0.25*HB)
          Q1 = HS_COEFF(IHYD,4)*W*HB*SQRT(2.*G*(HUD-HDD))      ! *** SUBMERGED ORIFICE FLOW
          Q2 = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))     ! *** SUB-CRITICAL WEIR FLOW
          QHS = RATIO*Q1 + (1.-RATIO)*Q2
          
        ELSE
          ! *** SUB-CRITICAL WEIR FLOW
          QHS = HS_COEFF(IHYD,2)*W*HDD*SQRT(2.*G*(HUD-HDD))    ! *** SUB-CRITICAL WEIR FLOW
          
        ENDIF
             
      ENDIF
      
    ELSEIF( HYD_STR(NCTL).NQCTYP == 7  )THEN 
      ! ** WEIR
      ! ** FORMULA FROM http://www.codecogs.com/library/engineering/fluid_mechanics/weirs/discharge.php
      HUD = ZHU - USINV
      HDD = ZHD - USINV
      DIA = HUD
      CDD = HS_COEFF(IHYD,2)
      IF( HS_XSTYPE(IHYD) == 5  )THEN
        ! ** RECTANGULAR WEIR       
        QHS = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5   ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        
      ELSEIF( HS_XSTYPE(IHYD) == 6  )THEN
        ! ** TRIANGLE NOTCH
        QHS = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
        
      ELSEIF( HS_XSTYPE(IHYD) == 7  )THEN
        ! ** TRAPEZOIDAL 
        QHS1 = 2.0D0/3.0D0*CDD*HS_WIDTH(IHYD)*SQRT(2.0*G)*HUD**1.5   ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        QHS2 = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
        QHS  = QHS1 + QHS2

      ELSEIF( HS_XSTYPE(IHYD) == 8  )THEN
        ! ** BROAD CRESTED WEIR (CLH**1.5)
        QHS = CDD*WOPEN*HUD**1.5        ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        
      ELSE
        PRINT '(A29,2I4)', ' *** BAD CROSS-SECTION TYPE:',IHYD, HS_XSTYPE(IHYD)
        CALL STOPP('.')
        
      ENDIF
      
      ! *** ACCOUNT FOR SUBMERGENCE (Villemonte_1947_Submerged_weir_discharge_studies)
      IF( ZHD > USINV )THEN
        QHS = QHS*(1.-HDD/HUD)**0.385
      ENDIF
      
    ELSEIF( HYD_STR(NCTL).NQCTYP == 8  )THEN 
      ! *** ORIFICE
      HDD = ZHD - USINV
      W  = WOPEN               ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      HB = HOPEN                ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      CDD = HS_COEFF(IHYD,2)
      
      CALL CROSS_SECTION(IHYD,HUD,HOPEN,WOPEN,SOPEN,DIA,FAREA,WETPER,HRAD,TWID)       

      IF( DIA <= 1.E-6 )RETURN

      ! ** COMPARISON OF HUD & HDD
      IF( HDD >= DIA  )THEN
        ! *** SUBMERGED ORIFICE
        QHS = CDD*FAREA*(2.*G*(HUD-HDD))**0.5
        
      ELSEIF( HUD > DIA .AND. HDD < DIA )THEN
        ! *** FREE FLOWING ORIFICE JET
        HUD = ZHU - USINV + 0.5*HB  ! *** HUD BASED ON CENTERLINE
        QHS = CDD*FAREA*(2.*G*HUD)**0.5
        
        ! *** ACCOUNT FOR SUBMERGENCE (Villemonte_1947_Submerged_weir_discharge_studies)
        IF( ZHD > USINV )THEN
          QHS = QHS*(1.-HDD/HUD)**0.385
        ENDIF

      ELSE
        ! *** HUD < DIA.  TREAT AS WEIR FLOW
        CDD = HS_COEFF(IHYD,1)
        IF( HS_XSTYPE(IHYD) == 5  )THEN
          ! ** RECTANGULAR WEIR       
          QHS = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5                      ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
        
        ELSEIF( HS_XSTYPE(IHYD) == 6  )THEN
          ! ** TRIANGLE NOTCH
          QHS = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
        
        ELSEIF( HS_XSTYPE(IHYD) == 7  )THEN
          ! ** TRAPEZOIDAL 
          QHS1 = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5                     ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
          QHS2 = 8.0D0/15.0D0*CDD*SQRT(2.0*G)*TAN(0.5*HS_ANGLE(IHYD))*HUD**2.5
          QHS  = QHS1 + QHS2

        ELSE
          ! *** CATCH ALL WEIR SHAPES (CIRCLES, HALF-CIRCLES, ETC).  BASE ON AREA.
          ! **  CDD MUST INCLUDE ALL GEOMETRIC EFFECTS
          QHS = CDD*FAREA*SQRT(2.0*G)*HUD**1.5
        ENDIF
      
        ! *** ACCOUNT FOR SUBMERGENCE (Villemonte_1947_Submerged_weir_discharge_studies)
        IF( ZHD > USINV )THEN
          QHS = QHS*(1.-HDD/HUD)**0.385
        ENDIF
      ENDIF

    ELSEIF( HYD_STR(NCTL).NQCTYP == 9  )THEN 
      ! *** FLOATING SKIMMER WALL
      CALL STOPP('*** HY-STRUCTURE IS NOT SUPPORTED: 9')
      
    ELSEIF( HYD_STR(NCTL).NQCTYP == 10  )THEN 
      ! *** SUBMERGED WEIR
      QHS = 2.0D0/3.0D0*CDD*WOPEN*SQRT(2.0*G)*HUD**1.5                        ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING

    ENDIF 
    HSCTL(NCTL).CUR.FLOW = QHS
    
1000 CONTINUE    
    ! ** HORIZONTAL DISTRIBUTION
    QHS = QHS*HYD_STR(NCTL).HS_FACTOR
    
    ! *** SET THE CURRENT HYDRAULIC STRUCTURE FLOWS BY LAYER
    ! *** USE LAYER ELEVATIONS AND INVERTS TO SELECT WHICH LAYERS TO WITHDRAW FROM
    IF( KC == 1 )THEN
      HSZ(KC) = 1.0
      CUMZ = 1.
    ELSE
      HSZ = 0.
      CUMZ = 0.
      TOPZ = MIN(DIA+USINV,ZHU)
      DO K=1,KC
        TOLZ = BELV(LU) + Z(LU,K)*HP(LU)
        IF( TOLZ > USINV )THEN
          DELZ = MIN( (TOLZ-USINV)/DIA, 1.0) - CUMZ
          HSZ(K) = DELZ
          CUMZ = CUMZ + DELZ
          IF( CUMZ >= 1.0 .OR. TOLZ > TOPZ )EXIT
        ENDIF
      ENDDO
    ENDIF
    
    ! *** ADD FLOWS TO GLOBAL VARIABLE QSUM
    DO K=1,KC
      QCTLT(K,NCTL,1) = QHS*HSZ(K)/CUMZ
      QSUM(LU,K) = QSUM(LU,K) - QCTLT(K,NCTL,1)      ! *** UPSTREAM
    ENDDO
    
    ! *** USE LAYER ELEVATIONS AND INVERTS TO SELECT WHICH LAYERS TO RETURN THE FLOW TO
    IF( LD /= LC )THEN
      ! *** DOWNSTREAM ACTIVE
      HSZ = 0.
      CUMZ = 0.
      TOPZ = MIN(DIA+DSINV,ZHD)
      IF( (DSINV+0.0001) >= TOPZ .OR. KC == 1 )THEN
        ! *** PLACE ALL OF THE FLOW IN THE TOP LAYER
        HSZ(KC) = 1.0
        CUMZ = 1. 
      ELSE
        DO K=1,KC
          TOLZ = BELV(LD) + Z(LD,K)*HP(LD)
          IF( TOLZ > DSINV )THEN
            DELZ = MIN( (TOLZ-DSINV)/DIA, 1.0) - CUMZ
            HSZ(K) = DELZ
            CUMZ = CUMZ + DELZ
            IF( CUMZ >= 1.0 .OR. TOLZ > TOPZ )EXIT
          ENDIF
        ENDDO
      ENDIF
      IF( ABS(CUMZ - 1.0) > 0.01 )THEN
        !CUMZ = 1.                        ! DELME
      ENDIF

      ! *** ADD RETURN FLOWS TO GLOBAL VARIABLE QSUM
      DO K=1,KC
        QCTLT(K,NCTL,2) = QHS*HSZ(K)/CUMZ
        QSUM(LD,K) = QSUM(LD,K) + QCTLT(K,NCTL,2)  ! *** DOWNSTREAM
      ENDDO
    ENDIF
    
    
  END SUBROUTINE COMPUTE_HSFLOW
  
  SUBROUTINE CROSS_SECTION(IHYD,HUD,HOPEN,WOPEN,SOPEN,DIA,FAREA,WETPER,HRAD,TWID,KCON)
    ! ** CALCULATE HYD. STRUCTURE PARAMETERS:
    ! ** FAREA  : FLOW AREA            [M2]
    ! ** WETPER : WET PERIMETER        [M]
    ! ** HRAD   : HYDRAULIC RADIUS     [M]
    ! ** KCON   : CONVEYANCE           [M3/S]
    ! ** TWID   : TOP WIDTH            [M]
    ! ** DIA    : DIAMETER OR HEIGHT   [M]
    ! ** HOPEN  : OPENING HEIGHT       [M]
    ! ** WOPEN  : OPENING WIDTH        [M]
    ! ** SOPEN  : CHANGE IN SILL LEVEL [M]
    INTEGER(4),PARAMETER   :: NMAX=100
    INTEGER(4),INTENT(IN ) :: IHYD
    REAL,      INTENT(IN ) :: HUD,HOPEN,WOPEN,SOPEN
    REAL,      INTENT(OUT) :: FAREA,WETPER,HRAD,TWID,DIA
    REAL,OPTIONAL,INTENT(OUT) :: KCON
    REAL       :: RA,RB,XH,DEL,ALP,WETPER4,YH,A
    INTEGER(4) :: NN
    REAL       :: X(NMAX),Y(NMAX)
  
    IF( HS_XSTYPE(IHYD) == 1 .OR.  HS_XSTYPE(IHYD) == 2  )THEN
      ! ** 1 - CIRCLE: RA = HS_WIDTH, RB IS NOT USED
      ! ** 2 - HALF CIRCLE: RA = HS_WIDTH, RB IS NOT USED
      RA = HS_WIDTH(IHYD)/2.
      DIA = 2.*RA
      IF( HS_XSTYPE(IHYD) == 2 )THEN
        YH = MIN(HUD+RA,DIA)
      ELSE
        YH = MIN(HUD,DIA)
      ENDIF
      IF( YH >= DIA  )THEN
        ! *** PIPE IS FULL
        FAREA = PI*RA**2
        WETPER = 2*PI*RA
        TWID = 0
      ELSE
        ! *** PARTIALLY FILLED PIPE
        DEL = ABS(RA-YH)
        ALP = ACOS(DEL/RA)
        TWID = 2*SQRT(RA**2-DEL**2)
        FAREA = ALP*RA**2 - 0.5*DEL*TWID
        WETPER = 2*ALP*RA
        IF( YH > RA )THEN
          FAREA = PI*RA**2 - FAREA
          WETPER = 2*PI*RA-WETPER
        ENDIF
      ENDIF
      IF( HS_XSTYPE(IHYD) == 2 )THEN
        ! *** HALF PIPE CORRECTION: YH > RA
        FAREA  = FAREA - 0.5*PI*RA**2
        WETPER = WETPER - PI*RA + 2.*RA
        DIA = RA
      ENDIF
      
    ELSEIF( HS_XSTYPE(IHYD) == 3 .OR.  HS_XSTYPE(IHYD) == 4  )THEN
      ! ** ELLIP: SEMI-AXIS RA = HS_WIDTH/2, RB = HS_HEIGHT/2
      RA = HS_WIDTH(IHYD)/2.
      DIA = HS_HEIGHT(IHYD)

      IF( HS_XSTYPE(IHYD) == 4 )THEN
        RB = HS_HEIGHT(IHYD)
        YH = MIN(HUD+RB,DIA)
      ELSE
        RB = HS_HEIGHT(IHYD)/2.
        YH = MIN(HUD,DIA)
      ENDIF

      X = (/ ((NN-1)*RA/(NMAX-1),NN=1,NMAX) /) 
      Y  = RB*SQRT(1.0-X**2/RA**2)
      WETPER4 = 0
      DO NN=2,NMAX
        WETPER4=WETPER4+SQRT((X(NN)-X(NN-1))**2+(Y(NN)-Y(NN-1))**2) 
      ENDDO
      
      IF( YH >= 2.*RB )THEN
        ! ** FULL FLOW
        FAREA = PI*RA*RB
        WETPER = 4.*WETPER4
        TWID = 0
      ELSE
        ! ** PARTIAL FLOW
        DEL = ABS(RB-YH)
        XH = RA*SQRT(1.0-DEL**2/RB**2) 
        TWID = 2.*XH

        ! ** YH < RB        
        WETPER = 0
        DO NN=2,NMAX
          IF( X(NN) <= XH )THEN
            WETPER=WETPER+SQRT((X(NN)-X(NN-1))**2+(Y(NN)-Y(NN-1))**2)
          ELSE
            EXIT
          ENDIF
        ENDDO
        WETPER = 2*WETPER
    
        FAREA = 0
        DO NN=2,NMAX
          IF( X(NN) <= XH )THEN
            FAREA=FAREA + 0.5*(X(NN)-X(NN-1))*(Y(NN)+Y(NN-1))
          ELSE
            EXIT
          ENDIF
        ENDDO
        FAREA = 2.*( FAREA - DEL*XH )   
    
        IF( YH > RB )THEN
          WETPER = 4.*WETPER4 - WETPER
          FAREA  = PI*RA*RB - FAREA
        ENDIF
      ENDIF
    
      IF( HS_XSTYPE(IHYD) == 4 )THEN
        ! *** HALF PIPE CORRECTION     
        FAREA  = FAREA - 0.5*PI*RA*RB
        WETPER = WETPER - 2.*WETPER4 + 2.*RA
      ENDIF
      
    ELSEIF( HS_XSTYPE(IHYD) == 5  )THEN
      ! ** RECTANGLE: RA = HS_WIDTH, RB = HS_HEIGHT
      RA = WOPEN                                    ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      RB = HOPEN                                    ! 2017-11-06, NTL: ACCOUNTS FOR GATE OPENING
      DIA = RB
      YH = MIN(HUD,RB)
      TWID = RA
      FAREA  = TWID*YH
      WETPER = TWID + 2.*YH   
      IF( HUD > RB ) WETPER = WETPER + RB
      
    ELSEIF( HS_XSTYPE(IHYD) == 6  )THEN
      ! ** V NOTCH: RB = HS_HEIGHT, HS_ANGLE
      DIA   = HS_HEIGHT(IHYD)
      RA    = DIA*TAN(0.5*HS_ANGLE(IHYD))
      TWID  = 2.*RA
      FAREA = DIA*RA
      WETPER= 2.*SQRT(RA**2+DIA**2)
      
    ELSEIF( HS_XSTYPE(IHYD) == 7  )THEN
      ! ** TRAPEZOID: RA = HS_WIDTH, RB = HS_HEIGHT
      DIA   = HS_HEIGHT(IHYD)      
      RA    = DIA*TAN(0.5*HS_ANGLE(IHYD))
      TWID  = HS_WIDTH(IHYD) + 2.*RA    
      FAREA = DIA*HS_WIDTH(IHYD) + DIA*RA
      WETPER= 2.*SQRT(RA**2+DIA**2) + HS_WIDTH(IHYD)
      
    ELSEIF( HS_XSTYPE(IHYD) == 8  )THEN
      ! ** PARABOLA: RA = HS_WIDTH, RB = HS_HEIGHT
      RA = HS_WIDTH(IHYD)/2.
      RB = HS_HEIGHT(IHYD)/2.
      DIA = 2.*RB
      YH = MIN(HUD,2.*RB)
      A  = 2*RB/RA**2
      X = (/ ((NN-1)*RA/(NMAX-1),NN=1,NMAX) /) 
      Y  = A*X**2
      XH = SQRT(YH/A) 
      WETPER = 0
      DO NN=2,NMAX
        IF( X(NN) <= XH )THEN
          WETPER=WETPER+SQRT((X(NN)-X(NN-1))**2+(Y(NN)-Y(NN-1))**2)
        ELSE
          EXIT
        ENDIF
      ENDDO
      TWID = 2*XH
      WETPER = 2*WETPER  
      FAREA = 4.0*XH*YH/3.0    
    ENDIF
    
    HRAD = FAREA/WETPER                          ! *** HYDRAULIC RADIUS
    IF( PRESENT(KCON) )THEN
      IF( HS_MANN(IHYD) > 0 )THEN
        KCON = FAREA/HS_MANN(IHYD)*HRAD**0.6667  ! *** CONVEYANCE
      ELSE
        KCON = 0.
      ENDIF
    ENDIF
  
  END SUBROUTINE CROSS_SECTION
  
  SUBROUTINE SWAP_USDS(NCTL,LU,LD,ZHU,ZHD,USINV,DSINV)
    ! *** BASED ON WATER SURFACE ELEVATIONS, SWAP UPSTREAM AND DOWNSTREAM CELLS FOR FLOW AND TRANSPORT
    INTEGER, INTENT(IN ) :: NCTL
    INTEGER, INTENT(OUT) :: LU,LD
    REAL, INTENT(OUT)   :: ZHU,ZHD
    REAL, INTENT(INOUT) :: USINV,DSINV
    INTEGER(4) :: ITMP,JTMP
    REAL       :: TMP
    
    ! ** SWAP
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

    INTEGER, INTENT(IN)::L,NX
    INTEGER :: NBAD
    
    IF( HS_XSTYPE(L) < 0 .OR. HS_XSTYPE(L) > 8 )THEN
      PRINT*,' *** BAD CROSS-SECTION TYPE, EQ#, HS_XSTYPE = ',L,HS_XSTYPE(L)
      CALL STOPP('.')
    ENDIF   
    NBAD = 0
    IF( NX == 5 .OR. NX == 8 )THEN
      SELECT CASE(HS_XSTYPE(L))
      CASE (1,2)
        IF(HS_WIDTH(L) <= 1E-6 )THEN
          ! ** 1 - CIRCLE: RA = HS_WIDTH, RB IS NOT USED
          ! ** 2 - HALF CIRCLE: RA = HS_WIDTH, RB IS NOT USED
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          CALL STOPP('.') 
        ENDIF
        
      CASE (3:5)
        ! ** ELLIP/RECTANGLE/PARABOL: SEMI-AXIS RA = HS_WIDTH/2, RB = HS_HEIGHT/2
        IF( HS_WIDTH(L) <= 1E-6  )THEN
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          NBAD = 1
        ENDIF
        IF( HS_HEIGHT(L) <= 1E-6 )THEN
          PRINT*,' *** BAD HS_HEIGHT: EQ#, HS_HEIGHT = ',L,HS_HEIGHT(L)
          NBAD = 1
        ENDIF 
        IF( NBAD > 0 ) CALL STOPP('.')
        
      CASE (6)
        ! ** V NOTCH: RB = HS_HEIGHT, HS_ANGLE
        IF( HS_HEIGHT(L) <= 1E-6  )THEN
          PRINT*,' *** BAD HS_HEIGHT: EQ#, HS_HEIGHT = ',L,HS_HEIGHT(L)
          NBAD = 1
        ENDIF
        IF( HS_ANGLE(L) <= 1E-6 )THEN
          PRINT*,' *** BAD HS_ANGLE: EQ#, HS_ANGLE = ',L,HS_ANGLE(L)
          NBAD = 1
        ENDIF   
        IF( NBAD > 0 ) CALL STOPP('.')
          
      CASE (7)
        ! ** TRAPEZOID: RA = HS_WIDTH, RB = HS_HEIGHT, HS_ANGLE
        IF( HS_WIDTH(L) <= 1E-6  )THEN
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          NBAD = 1
        ENDIF
        IF( HS_HEIGHT(L) <= 1E-6  )THEN
          PRINT*,' *** BAD HS_HEIGHT: EQ#, HS_HEIGHT = ',L,HS_HEIGHT(L)
          NBAD = 1
        ENDIF
        IF( HS_ANGLE(L) <= 1E-6 )THEN
          PRINT*,' *** BAD HS_ANGLE: EQ#, HS_ANGLE = ',L,HS_ANGLE(L)
          NBAD = 1
        ENDIF     
        IF( NBAD > 0 ) CALL STOPP('.')

      CASE (8)
        ! ** BROAD CRESTED WEIR: RA = HS_WIDTH
        IF( HS_WIDTH(L) <= 1E-6  )THEN
          PRINT*,' *** BAD HS_WIDTH: EQ#, HS_WIDTH = ',L,HS_WIDTH(L)
          CALL STOPP('.')
        ENDIF
      END SELECT
    ENDIF
        
  END SUBROUTINE HYDSTRUCT_CHECK
  
  SUBROUTINE PUMP_OPERATION_RULES(CTL, HVAL, CVAL, RVAL)
    ! *** CTL - STRUCTURE FOR PUMP
    ! *** HVAL - WATER LEVEL (HEAD) OF TRIGGER
    ! *** CVAL - OLD FLOWS
    ! *** RVAL - NEW FLOWS
  
    TYPE(STRCONTROL), INTENT(INOUT) :: CTL
    REAL, INTENT(IN) :: HVAL, CVAL
    REAL, INTENT(OUT) :: RVAL
    REAL(RKD) :: TSUM
    INTEGER :: K, ID

    ID = CTL.ID
    CTL.SUBID = 0
    TSUM = CTL.TOPEN + DELT
    
    IF( CTL.CUR.STATE == 2 )THEN
      ! *** GATE OPENING
      RVAL = CVAL + CTL.CUR.RATE * DELT/60.
      CTL.TOPEN = TSUM
      IF( RVAL >= CTL.LIM.FLOW )THEN
        RVAL = CTL.LIM.FLOW
        CTL.CUR.STATE = 1                   ! *** FULLY OPEN OR MAXIMUM FLOW
        CTL.CUR.RATE = 0.
      ENDIF
    ELSEIF (CTL.CUR.STATE == 3 )THEN
      ! *** GATE CLOSING
      RVAL = CVAL - CTL.CUR.RATE * DELT/60.
      CTL.TOPEN = TSUM
      IF( RVAL <= CTL.LIM.FLOW )THEN
        RVAL = CTL.LIM.FLOW
        CTL.CUR.STATE = 0                   ! *** FULLY CLOSED OR MINIMUM FLOW
        CTL.CUR.RATE = 0.
        CTL.TOPEN = 0.
      ENDIF
    ELSE
      ! *** FULLY OPEN OR CLOSED
      RVAL = CVAL
      
      ! *** CHECK FOR START OF OPENING
      DO K=1,QCTRULES(ID).NTRGON
        IF( (HVAL >= QCTRULES(ID).TRGON(K).LEVEL) .AND. &
          ((CTL.CUR.STATE == 0) .OR. (QCTRULES(ID).TRGON(K).LEVEL > CTL.CUR.LEVEL)) )THEN
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
        ENDIF
      ENDDO

      ! *** CHECK FOR START OF CLOSING
      DO K=QCTRULES(ID).NTROFF,1,-1
        IF( (HVAL <= QCTRULES(ID).TROFF(K).LEVEL) .AND. &
          ! *** WATER LEVELS TRIGGER CLOSING
          ((CTL.CUR.STATE == 1) .OR. (QCTRULES(ID).TROFF(K).LEVEL < CTL.CUR.LEVEL)) )THEN
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
        ENDIF
      ENDDO
      CTL.TOPEN = TSUM   
    ENDIF
  END SUBROUTINE PUMP_OPERATION_RULES
  
  SUBROUTINE GATE_OPERATION_RULES(NCTL, HVAL, CVAL, LVAL, RVAL)
    ! *** NCTL - INDEX OF STRUCTURE FOR GATE
    ! *** HVAL - WATER LEVEL (HEAD) OF TRIGGER
    ! *** CVAL - OLD LEVEL
    ! *** LVAL - TRIGGER LEVEL
    ! *** RVAL - NEW LEVEL

    INTEGER, INTENT(IN) :: NCTL
    REAL, INTENT(IN) :: HVAL, CVAL, LVAL
    REAL, INTENT(OUT) :: RVAL
    REAL(RKD) :: TSUM
    INTEGER :: K, ID

    ID = HSCTL(NCTL).ID
    HSCTL(NCTL).SUBID = 0
    
    TSUM = HSCTL(NCTL).TOPEN + DELT
    IF( HSCTL(NCTL).CUR.STATE == 2 )THEN
      ! *** GATE OPENING
      RVAL = CVAL + HSCTL(NCTL).CUR.RATE * DELT/60.
      HSCTL(NCTL).TOPEN = TSUM
      IF( RVAL >= LVAL )THEN
        RVAL = LVAL
        HSCTL(NCTL).CUR.STATE = 1           ! *** FULLY OPEN OR MAXIMUM OPENING
        HSCTL(NCTL).CUR.RATE = 0.
        
        IF( .TRUE. ) WRITE(6, '("GATE ",I3," IS FULLY OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
        WRITE(mpi_log_unit, '("GATE ",I3," IS FULLY OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        CLOSE(mpi_log_unit)
      ENDIF
      
    ELSEIF (HSCTL(NCTL).CUR.STATE == 3 )THEN
      ! *** GATE CLOSING
      RVAL = CVAL - HSCTL(NCTL).CUR.RATE * DELT/60.
      HSCTL(NCTL).TOPEN = TSUM
      IF( RVAL <= LVAL )THEN
        RVAL = LVAL
        HSCTL(NCTL).CUR.STATE = 0           ! *** FULLY CLOSED OR MINIMUM CLOSURE
        HSCTL(NCTL).CUR.RATE = 0.
        HSCTL(NCTL).TOPEN = 0.

        IF( .TRUE. ) WRITE(6, '("GATE ",I3," IS FULLY CLOSED (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
        WRITE(mpi_log_unit, '("GATE ",I3," IS FULLY CLOSED (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
        CLOSE(mpi_log_unit)
      ENDIF
      
    ELSE
      ! *** FULLY OPEN OR CLOSED
      RVAL = CVAL                           ! *** SET INITIALLY OPENING TO FINAL, I.E. NO CHANGE IN OPENING
      
      ! *** IF CLOSED, CHECK FOR START OF OPENING
      DO K=1,QCTRULES(ID).NTRGON
        IF( HVAL >= QCTRULES(ID).TRGON(K).LEVEL )THEN
          IF( HSCTL(NCTL).CUR.STATE == 0 .OR. QCTRULES(ID).TRGON(K).LEVEL > HSCTL(NCTL).CUR.LEVEL )THEN
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

            IF( .TRUE. ) WRITE(6,'("GATE ",I3," STARTING TO OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
            WRITE(mpi_log_unit,'("GATE ",I3," STARTING TO OPEN (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            CLOSE(mpi_log_unit)
            
            EXIT
          ENDIF
        ENDIF
      ENDDO

      ! *** IF FULLY OPEN, CHECK FOR START OF CLOSING
      DO K=QCTRULES(ID).NTROFF,1,-1
        IF( HVAL <= QCTRULES(ID).TROFF(K).LEVEL )THEN
          IF( HSCTL(NCTL).CUR.STATE == 1 .OR. QCTRULES(ID).TROFF(K).LEVEL < HSCTL(NCTL).CUR.LEVEL )THEN
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

            IF( .TRUE. ) WRITE(6, '("GATE ",I3," STARTING TO CLOSE (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
            WRITE(mpi_log_unit, '("GATE ",I3," STARTING TO CLOSE (",F8.3," M) AT:",I10,F12.5,F10.3)') NCTL, RVAL, NITER, TIMEDAY, HVAL
            CLOSE(mpi_log_unit)
            !PRINT '(20X,3F10.3)', HVAL, LVAL, HSCTL(NCTL).CUR.LEVEL  ! DELME
            EXIT
          ENDIF
        ENDIF
      ENDDO
      HSCTL(NCTL).TOPEN = TSUM
    ENDIF
    
  END SUBROUTINE GATE_OPERATION_RULES
  
  SUBROUTINE GATE_OPENING_INTERP(TS,T,H,W,S)
    TYPE(CONTROLSERIES), INTENT(IN) :: TS
    REAL(RKD), INTENT(IN) :: T
    REAL, INTENT(OUT) :: H,W,S
    REAL :: FACTOR
    INTEGER :: K

    IF( T < TS.TIME(1) )THEN
      K = 1
    ELSEIF (T > TS.TIME(TS.COUNT) )THEN
      K = TS.COUNT - 1
    ELSE
      DO K=1,TS.COUNT-1
        IF( T >= TS.TIME(K) .AND. T <= TS.TIME(K+1)) EXIT
      ENDDO
    ENDIF
    FACTOR = TS.TIME(K+1) - TS.TIME(K)
    IF( ABS(FACTOR) < 1.0E-9 )THEN
      FACTOR = 0.5
    ELSE
      FACTOR = (T - TS.TIME(K))/FACTOR
    ENDIF
    H = TS.HEIGHT(K) + (TS.HEIGHT(K+1) - TS.HEIGHT(K))*FACTOR
    W = TS.WIDTH(K) + (TS.WIDTH(K+1) - TS.WIDTH(K))*FACTOR
    S = TS.SILL(K) + (TS.SILL(K+1) - TS.SILL(K))*FACTOR
  END SUBROUTINE GATE_OPENING_INTERP
  
  SUBROUTINE CHANGE_RATING_CURVE_BY_RULES(NCTL, IRET)
    INTEGER, INTENT(IN) :: NCTL
    INTEGER, INTENT(OUT) :: IRET
    REAL(RKD) :: TSUM, ZHU, ZHD, HVAL
    INTEGER :: LU, LD, K, ID, IVAL

    ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
    LU = LIJ(HSCTL(NCTL).IREFUP, HSCTL(NCTL).JREFUP)
    ZHU = BELV(LU) + HP(LU)
      
    IF( HSCTL(NCTL).ITYPE == 3 )THEN
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
        LD = LIJ(HSCTL(NCTL).IREFDN, HSCTL(NCTL).JREFDN)      
        ZHD = BELV(LD) + HP(LD)
        HVAL = ZHU - ZHD
    ELSE
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
        HVAL = ZHU
    ENDIF
    
    ID = HSCTL(NCTL).ID
    HSCTL(NCTL).SUBID = 0
    TSUM = HSCTL(NCTL).TOPEN + DELT
    IF( HVAL > HSCTL(NCTL).PREV )THEN  ! ASCENDING
      DO K=1,QCTRULES(ID).NTRGON
        IF( HVAL >= QCTRULES(ID).TRGON(K).LEVEL )THEN
          HSCTL(NCTL).CUR.STATE = 2      ! START OPENING
          IRET = QCTRULES(ID).TRGON(K).ID
          HSCTL(NCTL).CUR.ID = IRET
          HSCTL(NCTL).TSTART = TIMESEC
          HSCTL(NCTL).TOPEN = 0.
          HSCTL(NCTL).SUBID = K
          RETURN
        ENDIF
      ENDDO
    ELSE  ! DESCENDING
      HSCTL(NCTL).TOPEN = TSUM
      DO K=QCTRULES(ID).NTROFF,1,-1
        IF( HVAL <= QCTRULES(ID).TROFF(K).LEVEL )THEN
          HSCTL(NCTL).CUR.STATE = 3      ! STOP OPERATIONAL
          IRET = QCTRULES(ID).TROFF(K).ID
          HSCTL(NCTL).CUR.ID = IRET
          HSCTL(NCTL).TSTART = TIMESEC
          HSCTL(NCTL).TOPEN = 0.
          HSCTL(NCTL).SUBID = K
          RETURN
        ENDIF
      ENDDO
    ENDIF
    HSCTL(NCTL).PREV = HVAL
  END SUBROUTINE CHANGE_RATING_CURVE_BY_RULES
  
  SUBROUTINE CHANGE_RATING_CURVE_BY_TIME(TS,T,IRET)
    TYPE(CONTROLSERIES), INTENT(IN) :: TS
    REAL(RKD), INTENT(IN) :: T
    INTEGER, INTENT(OUT) :: IRET
    INTEGER :: K

    IF( T < TS.TIME(1) )THEN
      IRET = TS.ID(1)
    ELSEIF (T > TS.TIME(TS.COUNT) )THEN
      K = TS.COUNT
      IRET = TS.ID(K)
    ELSE
      DO K=1,TS.COUNT-1
        IF( T >= TS.TIME(K) .AND. T <= TS.TIME(K+1) )THEN
          IRET = TS.ID(K)
          EXIT
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE CHANGE_RATING_CURVE_BY_TIME
  
END MODULE
