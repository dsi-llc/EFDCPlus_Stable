! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALPUV2C 

  ! *** *******************************************************************C
  !
  ! ** SUBROUTINE CALPUV2C CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
  ! ** AND VHDXE, FOR FREE SURFACE FLOWS FOR THE 2TL SOLUTION.
  ! ** WITH PROVISIONS FOR WETTING AND DRYING OF CELLS
  !
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP
  !    2002-02-28    John Hamrick      Modified drying and wetting scheme.
  !                                      the old formulation remains  
  !       See (isdry > 0 .and. isdry < 98). The new formulation is activated  
  !       by ( isdry == 99 ). Also added option to waste water from essentially  
  !       dry cells having water depths greater than HDRY, i.e. the high and  
  !       wet cells blocked by dry cells. This is actived by a negative value  
  !       of ndrystp parameter is the efdc.inp file.  
  !       Added save of old values of horizontal flow face switches SUB1 & SVB1  
  !       and transport bypass mask, LMASKDRY for dry cells.
  !       added QDWASTE(L) to save source equivalent of volume loss rate  
  !       for reducing depth of high/dry cells.  Also added concentration  
  !       adjustments. 

  USE GLOBAL 
  USE EFDCOUT
  Use Variables_MPI

#ifdef _MPI
  Use MPI
  Use MPI_All_Reduce
  Use Communicate_Ghost_Routines
  Use Mod_Map_Write_EE_Binary
#endif

  IMPLICIT NONE

  INTEGER :: NMD, ITERHP, NCORDRY, ICORDRY, ND, LF, LL, L, LP, LS, LN, LW, LE, LHOST, LCHNU, LCHNV, LG
  INTEGER :: IUW, IUE, IVS, IVN, IFACE, LMAX, LMIN, NTMP, IMAX, IMIN, JMAX, JMIN, I, J, K, NNEGFLG, KM
  INTEGER, SAVE :: INOTICE, NCORDRYMAX, NCORDRYAVG, ITERMAX, ITERAVG, NITERAVG, NCOUNT
  INTEGER, SAVE :: NOPTIMAL
  INTEGER, SAVE :: LDMOPT
  REAL :: DELTD2, RLAMN, RLAMO, TMPX, TMPY, C1, TMPVAL, HDRY2, BELVAVG, SVPW1
  REAL :: SUBW, SUBE, SVBS, SVBN, DHPDT, DHPDT2, RDRY, CCMNM, CCMNMI, HDRY90
  REAL :: RVAL, HOLDTMP, RNPORI, DIVEXMX, DIVEXMN, DIVEX, ETGWTMP, ETGWAVL, EL, EW, ES
  REAL(RKD), SAVE :: DAYOLD, DAYOLD30
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: IACTIVE  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ICORDRYD
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: NNATDRY
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:) :: LNATDRY

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: QCHANUT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: QCHANVT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: QSUMTMP 
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SUB1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SVB1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: CCMNMD

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: FSGZUDXYPI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: FSGZVDXYPI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: HPOLD
  
  REAL(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS, TWAIT                 ! MODEL TIMING TEMPORARY VARIABLE

  ! *** New variables for MPI
  Real    :: CCMNM_Local
  Integer :: ICORDRY_GLobal !< Dry couunter for the entire domain, takes the sum from all processes
  Integer :: ierr, IUPDATE, IUPDATE_Global
  
  IF( .NOT. ALLOCATED(IACTIVE) )THEN
    ! *** SET THE OPTIMAL NUMBER OF THREADS.  USE 100 CELLS PER THREAD AS GENERAL RULE
    NOPTIMAL = MIN(NTHREADS,8,MAX(LA/1000,2))
    LDMOPT   = INT(FLOAT(LA-1)/FLOAT(NOPTIMAL)) + 1
    
    if( process_id == master_id )then
      WRITE(*,'(A,I5)')'FIRST CALL TO 2TL PRESSURE SOLUTION.  CALPUV THREADS:',NOPTIMAL
    endif

    ALLOCATE(IACTIVE(NCHANM))  
    ALLOCATE(QCHANUT(NCHANM))  
    ALLOCATE(QCHANVT(NCHANM))  
    ALLOCATE(QSUMTMP(LCM))  
    ALLOCATE(SUB1(LCM))  
    ALLOCATE(SVB1(LCM))
    ALLOCATE(CCMNMD(NOPTIMAL))
    ALLOCATE(ICORDRYD(NOPTIMAL))
    ALLOCATE(FSGZUDXYPI(LCM))
    ALLOCATE(FSGZVDXYPI(LCM))
    ALLOCATE(HPOLD(LCM))
    ALLOCATE(LNATDRY(NOPTIMAL,(INT(LCM/NOPTIMAL)+1)))
    ALLOCATE(NNATDRY(NOPTIMAL))
    
    IACTIVE = 0
    QCHANUT = 0.
    QCHANVT = 0.
    QSUMTMP = 0.
    SUB1 = SUB
    SVB1 = SVB
    ISCDRY = 0
    CCMNMD = 0
    ICORDRYD = 0
    RCX = 1.0
    RCY = 1.0
    RCX(1) = 0.  
    RCY(1) = 0.  
    RCX(LC) = 0.  
    RCY(LC) = 0.  
    LNATDRY = 0
    NNATDRY = 0
    INOTICE = -999
    NCORDRYMAX = 0
    NCORDRYAVG = 0
    ITERMAX = 0
    ITERAVG = 0
    NITERAVG = 0
    NCOUNT = 0
    DAYOLD = INT(TIMEDAY)
    DAYOLD30 = DAYOLD + 30.
    
    ! INITIALIZE DIAGONAL
    CC=1.0
    
    FSGZUDXYPI=1.
    FSGZVDXYPI=1.
    DO L=2,LA
      FSGZUDXYPI(L) =  0.5/DXYU(L)
      FSGZVDXYPI(L) =  0.5/DXYV(L)
    ENDDO
    DO L=2,LA
      HPOLD(L) = HP(L)
    ENDDO
    
  ENDIF !***End allocation and initialization of various arrays
  TTWAIT = 0.

  ! *** Select time step depending on whether a static or dynamic time step is used
  IF( ISDYNSTP == 0 )THEN  
    DELT = DT  
    DELTD2 = 0.5*DT  
    DELTI = 1./DELT  
  ELSE  
    DELT = DTDYN  
    DELTD2 = 0.5*DTDYN  
    DELTI = 1./DELT  
  ENDIF  
  ISTL=2  
  
  ! *** INITIALIZE SUBGRID SCALE CHANNEL INTERACTIONS  
  IF( MDCHH >= 1 )THEN  
    DO NMD=1,MDCHH  
      QCHANUT(NMD) = QCHANU(NMD)  
      QCHANVT(NMD) = QCHANV(NMD)  
    ENDDO  
  ENDIF  

  ! *** SET SWITCHES FOR DRYING AND WETTING
  HDRY90 = 0.9*HDRY
  ITERHP = 0  
  NCORDRY = 0  
  NNEGFLG = 0
  ICORDRY_GLobal = 0

  IF( ISDRY > 0 )THEN
    DO LL=1,NBCSOP
      L = LOBCS(LL)
      LOPENBCDRY(L) = .FALSE.
    ENDDO
  ENDIF

  ! ***************************************************************************
  ! *** CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N)
  IF( BSC > 1.E-6 )Then
    CALL CALEBI
  End if

  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LE,LS,LN,LW,LG,K,TMPX,TMPY,C1,RDRY,EL,EW,ES,NTMP)
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    
    ! *** SET SWITCHES FOR DRYING AND WETTING  
    IF( ISDRY > 0 )THEN
      DO L=LF,LL  
        ISCDRY(L) = 0
        SUB1(L) = SUB(L)
        SVB1(L) = SVB(L)
        OLDMASK(L) = LMASKDRY(L)
      ENDDO
  
      ! *** ZERO VOLUME WASTING ARRAY COUNTERS
      IF( NDRYSTP > 0 )THEN
        NNATDRY(ND)=0
        IF( ISBAL > 0 )THEN
          DO L=LF,LL  
            QDWASTE(L)=0.
          ENDDO
        ENDIF
      ENDIF

    ENDIF

    IF( BSC > 1.E-6 )THEN
      ! *** CALCULATE EXPLICIT EXTERNAL DENSITY GRADIENTS  
      IF( IGRIDV == 0 )THEN
        !$OMP SIMD PRIVATE(LW,LS)
        DO L = LF,LL  
          LW=LWC(L) 
          FPGXE(L) = -SBX(L)*HU(L)*GP*((BI2W(L)+BI2W(LW))*(HP(L)-HP(LW))+2.0*HU(L)*(BI1W(L)-BI1W(LW))+(BEW(L)+BEW(LW))*(BELV(L)-BELV(LW)))  

          LS = LSC(L)  
          FPGYE(L) = -SBY(L)*HV(L)*GP*((BI2S(L)+BI2S(LS))*(HP(L)-HP(LS))+2.0*HV(L)*(BI1S(L)-BI1S(LS))+(BES(L)+BES(LS))*(BELV(L)-BELV(LS)))  
        ENDDO  
      ELSE
        !$OMP SIMD PRIVATE(LW,LS)
        DO L=LF,LL  
          LW = LWC(L) 
          FPGXE(L) = -SBX(L)*HU(L)*GP*( (BI2W(L)+BI2E(LW))*(HPW(L)-HPE(LW)) + 2.0*HU(L)*(BI1W(L)-BI1E(LW)) + (BEW(L)+BEE(LW))*(BELVW(L)-BELVE(LW)) )  
          
          LS = LSC(L)  
          FPGYE(L) = -SBY(L)*HV(L)*GP*( (BI2S(L)+BI2N(LS))*(HPS(L)-HPN(LS)) + 2.0*HV(L)*(BI1S(L)-BI1N(LS)) + (BES(L)+BEN(LS))*(BELVS(L)-BELVN(LS)) )  
        ENDDO  
      ENDIF
    ENDIF
    
    ! *** CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS  
    !$OMP SIMD PRIVATE(LW,LS)
    DO L=LF,LL
      LS = LSC(L)
      LW = LWC(L) 
      FUHDYE(L) = UHDYE(L) - DELTD2*SUB(L)*HRUO(L)*HU(L)*(P(L)-P(LW)) + SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L)) + FCAXE(L) + FPGXE(L) - SNLT*FXE(L))   ! *** m3/s  

      FVHDXE(L) = VHDXE(L) - DELTD2*SVB(L)*HRVO(L)*HV(L)*(P(L)-P(LS)) + SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L)) - FCAYE(L) + FPGYE(L) - SNLT*FYE(L))   ! *** m3/s 
    ENDDO  

    ! *** SET IMPLICIT BOTTOM AND VEGETATION DRAG AS APPROPRIATE  
    IF( ISITB >= 1 )THEN  
      ! *** IMPLICIT BOTTOM DRAG WITH VEGETATION  
      DO L=LF,LL
        TMPX = 1.0
        TMPY = 1.0
        IF( UHE(L) /= 0.0) TMPX = U(L,KSZU(L))*HU(L)/UHE(L)
        IF( VHE(L) /= 0.0) TMPY=V(L,KSZV(L))*HV(L)/VHE(L)
        RCX(L) = 1./( 1. + TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L) + U(L,KSZU(L))*U(L,KSZU(L))) + DELT*FXVEGE(L) )
        RCY(L) = 1./( 1. + TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L) + V(L,KSZV(L))*V(L,KSZV(L))) + DELT*FYVEGE(L) )
        FUHDYE(L) = FUHDYE(L)*RCX(L)
        FVHDXE(L) = FVHDXE(L)*RCY(L)
      ENDDO
    ELSEIF( ISVEG > 0 )THEN
      ! *** IMPLICIT VEGETATION DRAG ONLY.  REDUCE BY THE TOTAL OF ENERGY
      DO L=LF,LL
        !                    S    1/S
        RCX(L) = 1./( 1. + DELT*FXVEGE(L) )  ! *** dimensionless
        RCY(L) = 1./( 1. + DELT*FYVEGE(L) )  ! *** dimensionless
        FUHDYE(L) = FUHDYE(L)*RCX(L)
        FVHDXE(L) = FVHDXE(L)*RCY(L)
      ENDDO
    ENDIF

    ! *** RESET BOUNDARY CONDITIONS SWITCHES
    IF( ISDRY > 0 )THEN  
#ifdef _MPI 
      IF( NITER < 100 .OR. MOD(NITER,50) == 0 .OR. ITERHPM == 0 )THEN
        ! *** FORCE A COMPLETE UPDATE EVERY 100 INTERATIONS
        DO L=LF,LL
          SUB(L) = SUBO(L)  
          SVB(L) = SVBO(L)  
          SBX(L) = SBXO(L)  
          SBY(L) = SBYO(L)  
        ENDDO  
      ELSE
        ! *** UPDATED TO REDUCE MPI COMMUNICATIONS
        DO L=LF,LL
          LE = LEC(L)
          LN = LNC(L)
          RDRY = SUB(L) + SUB(LE) + SVB(L) + SVB(LN)
          
          IF( RDRY < 0.5 )THEN
            IF( HP(L) > HDRY )THEN
              ! *** ISOLATED CELL, CHECK ADJACENT CELLS
              EL = BELV(L)  + HP(L)                         ! *** WSEL OF CURRENT CELL
              
              IF( EL > BELV(LWC(L)) )THEN                   ! *** WEST FACE
                SUB(LWC(L)) = SUBO(LWC(L))            
                SBX(LWC(L)) = SBXO(LWC(L))
              ENDIF
              IF( EL > BELV(LEC(L)) )THEN                   ! *** EAST FACE
                SUB(LEC(L)) = SUBO(LEC(L))            
                SBX(LEC(L)) = SBXO(LEC(L))
              ENDIF
              IF( EL > BELV(LSC(L)) )THEN                   ! *** SOUTH FACE
                SVB(L) = SVBO(L)            
                SBY(L) = SBYO(L)
              ENDIF
              IF( EL > BELV(LNC(L)) )THEN                   ! *** NORTH FACE
                SVB(LNC(L)) = SVBO(LNC(L))            
                SBY(LNC(L)) = SBYO(LNC(L))
              ENDIF   
            ENDIF
            
          ELSEIF( RDRY < 3.5 )THEN
            IF( RDRY == (SUBO(L) + SUBO(LE) + SVBO(L) + SVBO(LN)) ) CYCLE  ! *** ALL FACES ACTIVE
            
            ! *** AT LEAST ONE FACE IS INACTIVE.  CHECK DEPTHS
            
            ! *** CHECK EAST/WEST FACES
            IF( SUB(L) < 0.5 .AND. SUBO(L) > 0.5 )THEN
              LW = LWC(L)
              IF( LW > 1 )THEN
                EL = BELV(L)  + HP(L)                       ! *** WSEL OF CURRENT CELL
                EW = BELV(LW) + HP(LW)                      ! *** WSEL OF WEST CELL
                IF( HP(L) < HDRY .AND. HP(LW) < HDRY )THEN
                  ! *** DO NOTHING
                ELSEIF( HP(L) >= HDRY .AND. HP(LW) >= HDRY )THEN
                  ! *** BOTH CELLS ARE WET, RESET CURRENT CELL
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)
                  
                ELSEIF( HP(LW) < HDRY .AND. EL > EW )THEN
                  ! *** RESET WEST CELL
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)
                  
                  ! *** OPEN OTHER FACES IN WEST CELL
                  IF( EW > BELV(LWC(LW)) )THEN               ! *** WEST FACE
                    SUB(LW) = SUBO(LW)            
                    SBX(LW) = SBXO(LW)
                  ENDIF
                  IF( EW > BELV(LSC(LW)) )THEN               ! *** SOUTH FACE
                    SVB(LW) = SVBO(LW)            
                    SBY(LW) = SBYO(LW)
                  ENDIF
                  IF( EW > BELV(LNC(LW)) )THEN               ! *** NORTH FACE
                    SVB(LNC(LW)) = SVBO(LNC(LW))            
                    SBY(LNC(LW)) = SBYO(LNC(LW))
                  ENDIF
                  
                ELSEIF( HP(L) < HDRY .AND. EW >= EL )THEN
                  ! *** RESET CURRENT CELL
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)

                  ! *** OPEN OTHER FACES IN CURRENT CELL
                  IF( EL > BELV(LEC(L)) )THEN                ! *** EAST FACE
                    SUB(LEC(L)) = SUBO(LEC(L))            
                    SBX(LEC(L)) = SBXO(LEC(L))
                  ENDIF
                  IF( EL > BELV(LSC(L)) )THEN                ! *** SOUTH FACE
                    SVB(L) = SVBO(L)            
                    SBY(L) = SBYO(L)
                  ENDIF
                  IF( EL > BELV(LNC(L)) )THEN                ! *** NORTH FACE
                    SVB(LNC(L)) = SVBO(LNC(L))            
                    SBY(LNC(L)) = SBYO(LNC(L))
                  ENDIF

                ENDIF
              ENDIF
            ENDIF

            ! *** CHECK NORTH/SOUTH FACES
            IF( SVB(L) < 0.5 .AND. SVBO(L) > 0.5 )THEN
              LS = LSC(L)
              IF( LS > 1 )THEN
                EL = BELV(L)  + HP(L)
                ES = BELV(LS) + HP(LS)
                IF( HP(L) < HDRY .AND. HP(LS) < HDRY )THEN
                  ! *** DO NOTHING
                ELSEIF( HP(L) >= HDRY .AND. HP(LS) >= HDRY )THEN
                  ! *** BOTH CELLS ARE WET, RESET CURRENT CELL
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)
                  
                ELSEIF( HP(LS) < HDRY .AND. EL > ES )THEN
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)
                  
                  ! *** OPEN OTHER FACES IN SOUTH CELL
                  IF( ES > BELV(LWC(LS)) )THEN               ! *** WEST FACE
                    SUB(LS) = SUBO(LS)            
                    SBX(LS) = SBXO(LS)
                  ENDIF
                  IF( ES > BELV(LEC(LS)) )THEN               ! *** EAST FACE
                    SUB(LEC(LS)) = SUBO(LEC(LS))            
                    SBX(LEC(LS)) = SBXO(LEC(LS))
                  ENDIF
                  IF( ES > BELV(LSC(LS)) )THEN               ! *** SOUTH FACE
                    SVB(LSC(LS)) = SVBO(LSC(LS))            
                    SBY(LSC(LS)) = SBYO(LSC(LS))
                  ENDIF
                  
                ELSEIF( HP(L) < HDRY .AND. ES >= EL )THEN
                  ! *** RESET CURRENT CELL
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)

                  ! *** OPEN OTHER FACES IN CURRENT CELL
                  IF( EL > BELV(LEC(L)) )THEN                ! *** EAST FACE
                    SUB(LEC(L)) = SUBO(LEC(L))            
                    SBX(LEC(L)) = SBXO(LEC(L))
                  ENDIF
                  IF( EL > BELV(LWC(L)) )THEN                ! *** WEST FACE
                    SUB(L) = SUBO(L)            
                    SBX(L) = SBXO(L)
                  ENDIF
                  IF( EL > BELV(LNC(L)) )THEN                ! *** NORTH FACE
                    SVB(LNC(L)) = SVBO(LNC(L))            
                    SBY(LNC(L)) = SBYO(LNC(L))
                  ENDIF

                ENDIF
              ENDIF
            ENDIF
          ELSEIF( QSUME(L) > 0.0 )THEN
            SUB(L) = SUBO(L)
            SVB(L) = SVBO(L)
            SBX(L) = SBXO(L)
            SBY(L) = SBYO(L)
          ENDIF
        ENDDO
        
        ! *** ENSURE OPEN BC'S ARE "ON"
        IF( ND == 1 )THEN
          DO K = 1,NBCSOP
            L = LOBCS(K)
            SUB(L) = SUBO(L)
            SVB(L) = SVBO(L)
            SBX(L) = SBXO(L)
            SBY(L) = SBYO(L)
          ENDDO
          DO K=1,NBCSOP2
            L = LOBCS2(K)
            SUB(L) = SUBO(L)
            SVB(L) = SVBO(L)
            SBX(L) = SBXO(L)
            SBY(L) = SBYO(L)
          ENDDO 
        ENDIF
      ENDIF
#else     
      DO L=LF,LL
        SUB(L) = SUBO(L)  
        SVB(L) = SVBO(L)  
        SBX(L) = SBXO(L)  
        SBY(L) = SBYO(L)  
      ENDDO  
#endif
    ENDIF

    ! *** ADVANCE EXTERNAL VARIABLES  
    DO L=LF,LL  
      UHDY1E(L) = UHDYE(L)  
      VHDX1E(L) = VHDXE(L)  
      P1(L)   = P(L)  
      H1U(L)  = HU(L)  
      H1V(L)  = HV(L)  
      H1UI(L) = HUI(L)  
      H1VI(L) = HVI(L)  
      H2P(L)  = H1P(L)
      H1P(L)  = HP(L)
    ENDDO

    DO K=1,KC  
      DO L=LF,LL
        UHDY1EK(L,K) = UHDYEK(L,K)
        VHDX1EK(L,K) = VHDXEK(L,K) 
      ENDDO
    ENDDO
    DO K=1,KC  
      DO L=LF,LL
        H1PK(L,K) = HPK(L,K)   
      ENDDO
    ENDDO

    IF( ISGWIE >= 1 )THEN
      DO L=LF,LL  
        AGWELV2(L)=AGWELV1(L)  
        AGWELV1(L)=AGWELV(L)  
      ENDDO
    ENDIF  

    ! *** SET OLD TIME LEVEL TERMS IN CONTINUITY EQUATION FOR NON BOUNDARY POINTS  
    ! *** DXYP = DXP*DYP (M^2),  G-9.82 (M/S^2), P (M2/S2), FP1 (M4/S3)
    C1=0.5*G
    !$OMP SIMD PRIVATE(LE,LN)
    DO L=LF,LL
      LE = LEC(L)
      LN = LNC(L)
      FP1(L) = DELTI*DXYP(L)*P(L) - C1*( UHDYE(LE)-UHDYE(L)+VHDXE(LN)-VHDXE(L) )
    ENDDO  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO

  ! *** LAYER FACE BLOCKING
  IF( BSC > 0.0 .AND. NBLOCKED > 0 .AND. N > 1 )THEN
    ND = 1
    C1=0.5*G
    DO LP=1,NBLOCKED
      L = LBLOCKED(LP)
      IF( KSZ(L) == KC ) CYCLE
      
      ! *** LAYER U BLOCKING
      IF( SUB(L) > 0. .AND. (BLDRAFTU(LP)+BLSILLU(LP)) > 0. )THEN
        LW=LWC(L) 
        FPGXE(L) = -SBX(L)*HU(L)*GP*( (BI2W(L)+BI2E(LW))*(HPW(L)-HPE(LW)) + 2.0*HU(L)*(BI1W(L)-BI1E(LW)) + (BEW(L)+BEE(LW))*(BELVW(L)-BELVE(LW)) )                  ! *** m3/s

        FUHDYE(L) = UHDYE(L) - DELTD2*SUB(L)*HRUO(L)*HU(L)*(P(L)-P(LW)) + SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L)) + FCAXE(L) + FPGXE(L) - SNLT*FXE(L))   ! *** m3/s
      ENDIF

      ! *** LAYER V BLOCKING
      IF( SVB(L) > 0. .AND. (BLDRAFTV(LP)+BLSILLV(LP)) > 0. )THEN
        LS=LSC(L)  
        FPGYE(L) = -SBY(L)*HV(L)*GP*( (BI2S(L)+BI2N(LS))*(HPS(L)-HPN(LS)) + 2.0*HV(L)*(BI1S(L)-BI1N(LS)) + (BES(L)+BEN(LS))*(BELVS(L)-BELVN(LS)) )  
        
        FVHDXE(L) = VHDXE(L) - DELTD2*SVB(L)*HRVO(L)*HV(L)*(P(L)-P(LS)) + SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L)) - FCAYE(L) + FPGYE(L) - SNLT*FYE(L)) 
      ENDIF

      LE=LEC(L)
      LN=LNC(L)
      FP1(L) = DELTI*DXYP(L)*P(L) - C1*( UHDYE(LE)-UHDYE(L)+VHDXE(LN)-VHDXE(L) )
    ENDDO
  ENDIF
  
  ! *************************************************************************
  ! ***  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
  ! ***  HOST-GUEST CHANNEL INTERACTION FOR NON BOUNDARY POINTS 

  !----------------------------------------------------------------------C
  ! ***  REENTER AT 1000 FOR WETTING-DRYING CORRECTION AND CHANNEL 
  ! ***  INTERACTION
  1000 CONTINUE  
  !----------------------------------------------------------------------C
  
  ! *** SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM  
  CCMNMD = 1.E+18 
  
  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LS,LN,C1) 
  DO ND=1,NOPTIMAL  
    LF = 2+(ND-1)*LDMOPT  
    LL = MIN(LF+LDMOPT-1,LA)

    C1 = 0.5*G
    !$OMP SIMD PRIVATE(LE,LN)
    DO L=LF,LL
      LE = LEC(L)
      LN = LNC(L)
      ! ***  USE THE SUB & SVB SWITCHES FOR MOMENTUM FOR THE CURRENT WET/DRY ITERATION
      ! ***   m4/s3   m/s2            m3/s                 m3/s              m3/s               m3/s            m3/s
      FP(L) = FP1(L) - C1*( SUB(LE)*FUHDYE(LE) - SUB(L)*FUHDYE(L) + SVB(LN)*FVHDXE(LN) - SVB(L)*FVHDXE(L) - 2.0*QSUME(L) )    ! *** m4/s3
    ENDDO  

    IF( ISGWIE >= 1 )THEN  
      DO L=LF,LL  
        !       m4/s3   m/s2             m3/s
        FP(L) = FP(L) - G*SPB(L)*(EVAPSW(L)-QGW(L))   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
      ENDDO  
    ENDIF  

    C1 = -0.5*DELTD2*G
    DO L=LF,LL
      CS(L) = C1*SVB(L)*HRVO(L)*RCY(L)*HV(L) 
      CW(L) = C1*SUB(L)*HRUO(L)*RCX(L)*HU(L)  
    ENDDO
    !DIR$ NOFUSION
    !$OMP SIMD PRIVATE(LE,LN)
    DO L=LF,LL
      LE = LEC(L)
      LN = LNC(L) 
      CE(L) = C1*SUB(LE)*HRUO(LE)*RCX(LE)*HU(LE)  
      CN(L) = C1*SVB(LN)*HRVO(LN)*RCY(LN)*HV(LN) 
    ENDDO
    !DIR$ NOFUSION
    DO L=LF,LL
      ! *** SET THE CENTROID
      CC(L) = DELTI*DXYP(L) - CS(L)-CW(L)-CE(L)-CN(L)  
    ENDDO  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP SINGLE
  
  ! *** APPLY THE OPEN BOUNDARY CONDITIONS
  IF( NBCSOP > 0 ) CALL SETOPENBC(DELTD2, HU, HV, NCORDRY)    

  ! *** INSERT IMPLICT SUB-GRID SCALE CHANNEL INTERACTIONS  
  IF( MDCHH >= 1 )THEN
    RLAMN=QCHERR  
    RLAMO=1.-RLAMN  
    CALL SUBCHAN(QCHANUT,QCHANVT,IACTIVE,DELT)  
  ENDIF
  !$OMP END SINGLE

  ! *** SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM
  !$OMP DO PRIVATE(ND,LF,LL,L) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    DO L=LF,LL
      CCMNMD(ND) = MIN(CCMNMD(ND),CC(L))
      FPTMP(L)   = FP(L)  
    ENDDO  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! $OMP BARRIER
  !$OMP SINGLE
  CCMNM_Local =1.E+18  ! ***  Modified to have a local copy
  DO ND=1,NOPTIMAL
    CCMNM_Local = MIN(CCMNM_Local, CCMNMD(ND))
  ENDDO

#ifdef _MPI 
  Call DSI_All_Reduce(CCMNM_Local, CCMNM, MPI_Min, TTDS, 1, TWAIT)
  DSITIMING(11) = DSITIMING(11) + TTDS
  TTWAIT = TTWAIT + TWAIT
#else
  CCMNM = CCMNM_Local
#endif
  CCMNMI = 1./CCMNM

  !$OMP END SINGLE  
  ! $OMP BARRIER

  ! *** SCALE BY MINIMUM DIAGONAL  
  !$OMP DO PRIVATE(ND,LF,LL,L) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    DO L=LF,LL  
      CCS(L)   = CS(L)*CCMNMI  
      CCW(L)   = CW(L)*CCMNMI  
      CCE(L)   = CE(L)*CCMNMI  
      CCN(L)   = CN(L)*CCMNMI  
      CCC(L)   = CC(L)*CCMNMI  
      FPTMP(L) = FPTMP(L)*CCMNMI  
      CCCI(L)  = 1./CCC(L)  
    ENDDO  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL

  IF( MDCHH >= 1 )THEN
    DO NMD=1,MDCHH
      CCCCHH(NMD) = CCCCHH(NMD)*CCMNMI
    ENDDO
  ENDIF
  
  ! *********************************************************************************************
  ! *** CALL THE PRECONDITIONED CONJUGATE GRADIENT SOLVER
#ifdef _MPI
  IF( num_processors > 1 )THEN
    IF( MDCHH == 0 ) Call Congrad_MPI(NOPTIMAL,LDMOPT)   ! *** MSCHH>=1 not parallelized yet @todo

    TTDS = DSTIME(0)
    Call MPI_barrier(MPI_Comm_World, ierr)
    TTWAIT = TTWAIT + (DSTIME(0)- TTDS)

    TTDS = DSTIME(0)
    Call communicate_ghost_cells(P, 'P')
    TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  ELSE
    IF( MDCHH == 0 ) CALL CONGRAD(NOPTIMAL,LDMOPT)
    IF( MDCHH >= 1 ) CALL CONGRADC
  ENDIF
#else
  IF( MDCHH == 0 ) CALL CONGRAD(NOPTIMAL,LDMOPT)
  IF( MDCHH >= 1 ) CALL CONGRADC
#endif  
  ! *********************************************************************************************
  
  ITERMAX   = MAX(ITERMAX,ITER)
  ITERAVG   = ITERAVG + ITER
  NITERAVG  = NITERAVG + 1
  
  !if( process_id == 0 ) write(6,'(I10,3i5)') niter, process_id, NCORDRY, ITER   ! delme
  
  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND,LF,LL,L,LS,LW) 
  DO ND=1,NOPTIMAL  
    LF = 2+(ND-1)*LDMOPT  
    LL = MIN(LF+LDMOPT-1,LA)
    
    ! *** CELL FACE DISCHARGE (M3/S)
    !$OMP SIMD PRIVATE(LW,LS)
    DO L=LF,LL
      LS = LSC(L)  
      LW = LWC(L) 
      UHDYE(L) = SUB(L)*( FUHDYE(L) - DELTD2*HRUO(L)*RCX(L)*HU(L)*(P(L)-P(LW)) )
      VHDXE(L) = SVB(L)*( FVHDXE(L) - DELTD2*HRVO(L)*RCY(L)*HV(L)*(P(L)-P(LS)) )  
    ENDDO  
    
    ! *** UNIT WIDTH DISCHARGE AT CELL FACE (M2/S)
    !DIR$ NOFUSION
    DO L=LF,LL  
      UHE(L) = UHDYE(L)*DYIU(L)  
      VHE(L) = VHDXE(L)*DXIV(L)  
    ENDDO  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  ! *** CALCULATE NEW SUB-GRID SCALE CHANNEL EXCHANGE FLOWS  
  !$OMP SINGLE
  IF( MDCHH >= 1 )THEN  
    DO NMD=1,MDCHH  
      IF( IACTIVE(NMD) > 0 )THEN  
        LHOST=LMDCHH(NMD)  
        LCHNU=LMDCHU(NMD)  
        LCHNV=LMDCHV(NMD)  
        IF( MDCHTYP(NMD) == 1 )THEN  
          QCHANU(NMD)=CCCCHU(NMD)*QCHANUT(NMD)-RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNU))-RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNU))  
          QCHANV(NMD)=0.  
        ENDIF 
        IF( MDCHTYP(NMD) == 2 )THEN  
          QCHANU(NMD)=0.  
          QCHANV(NMD)=CCCCHU(NMD)*QCHANVT(NMD)-RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNV))-RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNV))  
        ENDIF  
      ELSE  
        QCHANV(NMD)=0.  
        QCHANVN(NMD)=0.  
        QCHANU(NMD)=0.  
        QCHANUN(NMD)=0.  
      ENDIF  
    ENDDO  
  ENDIF  
  !$OMP END SINGLE
  
  ! **  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL TRANSPORTS AT (N+1)
  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LN) 
  DO ND=1,NOPTIMAL  
    LF = 2+(ND-1)*LDMOPT  
    LL = MIN(LF+LDMOPT-1,LA)

    ! *** CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL TRANSPORTS AT (N+1)  
    !$OMP SIMD PRIVATE(LE,LN)
    DO L=LF,LL  
      LE = LEC(L)
      LN = LNC(L)  
      HP(L) = H1P(L) + DELTD2*DXYIP(L)*( 2.*QSUME(L) - ( UHDYE(LE)+UHDY1E(LE)-UHDYE(L)-UHDY1E(L)  &
                                                       + VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L) ) )  
    ENDDO 

    IF( ISGWIE >= 1 )THEN
      DO L=LF,LL  
        HP(L) = HP(L) - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
      ENDDO  
    ENDIF  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL
  
  ! ***  APPLY OPEN BOUNDARYS
  DO LL=1,NBCSOP
    L = LOBCS(LL)
    HP(L) = GI*P(L) - BELV(L)  

    ! *** CHECK FOR CONDITION OF PSERT<BELV
    IF( ISDRY > 0 )THEN

      ! ***  HP CHECKING IS FOR RADIATION BC'S **
      IF( (LOPENBCDRY(L) .OR. (H1P(L) < 0.9*HDRY .AND. HP(L) < H1P(L)) .OR. HP(L) < 0. ) .AND. ISCDRY(L) == 0 )THEN
        IF( HP(L) <= 0. .AND. LOPENBCDRY(L) )THEN
          PRINT '(A,I5,3F10.4,F12.4)',' WARNING!  NEG DEPTH AT OPEN BC: L,HP,H1P,H2P,TIMEDAY', Map2Global(L).LG, HP(L), H1P(L), H2P(L), TIMEDAY
          OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
          WRITE(mpi_log_unit,'(A,I5,3F10.4,F12.4)')' WARNING!  NEG DEPTH AT OPEN BC: L,HP,H1P,H2P,TIMEDAY', Map2Global(L).LG, HP(L), H1P(L), H2P(L), TIMEDAY
          CLOSE(mpi_log_unit)

          HP(L) = 0.2*HDRY
        ELSE
          HP(L) = MIN(MAX(H1P(L),0.1*HDRY),0.9*HDRY)
        ENDIF
        
        ISCDRY(L) = 1  

        LE = LEC(L)
        LN = LNC(L)
        
        SUB(L) = 0.
        SVB(L) = 0.
        SUB(LE) = 0.
        SVB(LN) = 0.
        SBX(L) = 0.  
        SBY(L) = 0.  
        SBX(LE) = 0.  
        SBY(LN) = 0.  
        UHDYE(L) = 0.
        UHDYE(LE) = 0.
        VHDXE(L)  = 0.
        VHDXE(LN) = 0.
        
        LOPENBCDRY(L) = .TRUE.
        CC(L)  = DELTI*DXYP(L)
        P(L)   = (HP(L)+BELV(L))*G
        FP1(L) = DELTI*DXYP(L)*P(L)
      ENDIF
    ENDIF
  ENDDO

  ! *** Exchange HP ghost values
#ifdef _MPI
  TTDS = DSTIME(0)
  Call MPI_barrier(MPI_Comm_World, ierr)
  TTWAIT = TTWAIT + (DSTIME(0)- TTDS)
  
  TTDS = DSTIME(0)
  Call communicate_ghost_cells(HP, 'HP')
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
#endif

  ! *** ADD CHANNEL INTERACTION EXCHANGES
  IF( MDCHH >= 1 )THEN
    DO NMD=1,MDCHH
      IF( IACTIVE(NMD) > 0 )THEN
        LHOST=LMDCHH(NMD)
        LCHNU=LMDCHU(NMD)
        LCHNV=LMDCHV(NMD)
        IF( MDCHTYP(NMD) == 1 )THEN
          TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUT(NMD))
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)
        ENDIF
        IF( MDCHTYP(NMD) == 2 )THEN
          TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVT(NMD))
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  ! *** CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY
  IF( ISDRY > 0 .AND. ISDRY < 98 )THEN
    ICORDRYD = 0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LN,LE) 
    DO ND=1,NOPTIMAL  
      LF = 2+(ND-1)*LDMOPT  
      LL = MIN(LF+LDMOPT-1,LA)

      DO L=LF,LL
        IF( HP(L) <= HDRY )THEN  
          IF( ISCDRY(L) == 0 )THEN  
            ISCDRY(L) = 1  
            ICORDRYD(ND) = 1  
          ENDIF  
          LE = LEC(L) 
          LN = LNC(L)  
          SUB(L)  = 0.  
          SVB(L)  = 0.  
          SUB(LE) = 0.  
          SVB(LN) = 0.  
          SBX(L)  = 0.  
          SBY(L)  = 0.  
          SBX(LE) = 0.  
          SBY(LN) = 0.  
        ENDIF  
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF

  ! *** CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY  
  IF( ISDRY == 99 )THEN  
    HDRY2  = 2.*HDRY  
    ICORDRYD = 0  
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LW,LE)  &
    !$OMP             PRIVATE(SUBW,SUBE,SVBS,SVBN,DHPDT,DHPDT2,RDRY,TMPVAL)
    DO ND=1,NOPTIMAL  
      LF = 2+(ND-1)*LDMOPT  
      LL = MIN(LF+LDMOPT-1,LA)

      DO L=LF,LL
        IF( HP(L) <= HDRY )THEN
          LE = LEC(L) 
          LN = LNC(L)  
          SUBW = SUB(L)  
          SUBE = SUB(LE)  
          SVBS = SVB(L)  
          SVBN = SVB(LN)  
          DHPDT = (HP(L) - H1P(L))
          
          ! *** ALLOW RE-WETTING
          IF( DHPDT > 0.0 )THEN 
            LW = LWC(L) 
            LS = LSC(L)  
            SUB(L)  = 0.0  
            SUB(LE) = 0.0  
            SVB(L)  = 0.0  
            SVB(LN) = 0.0  
            SBX(L)  = 0.0  
            SBX(LE) = 0.0  
            SBY(L)  = 0.0  
            SBY(LN) = 0.0  
            
            ! *** RAISING WATER, IS IT FAST ENOUGH TO STAY WET
            DHPDT2 = (HDRY - H1P(L))*0.01
            IF( DHPDT > DHPDT2 )THEN
              IF( UHDYE(L) /= 0.0 )THEN
                IF( HP(LW) > HDRY )THEN
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)
                ENDIF
              ENDIF  
              IF( UHDYE(LE) /= 0.0 )THEN
                IF( HP(LE) > HDRY )THEN
                  SUB(LE) = SUBO(LE)
                  SBX(LE) = SBXO(LE)  
                ENDIF  
              ENDIF  
              IF( VHDXE(L) /= 0.0 )THEN
                IF( HP(LS) > HDRY )THEN
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)
                ENDIF  
              ENDIF  
              IF( VHDXE(LN) /= 0.0 )THEN
                IF( HP(LN) > HDRY )THEN
                  SVB(LN) = SVBO(LN)
                  SBY(LN) = SBYO(LN)
                ENDIF  
              ENDIF  
              
              RDRY = SUB(L) + SUB(LE) + SVB(L) + SVB(LN)  
              IF( RDRY < 0.5 )THEN  
                ISCDRY(L) = 1  
              ELSE  
                ISCDRY(L) = 0  
              ENDIF  
              
              TMPVAL = ABS(SUB(L) - SUBW)  
              IF( TMPVAL > 0.5 )THEN
                ICORDRYD(ND) = 1  
              ELSE
                TMPVAL = ABS(SUB(LE) - SUBE)  
                IF( TMPVAL > 0.5 )THEN
                  ICORDRYD(ND) = 1  
                ELSE
                  TMPVAL = ABS(SVB(L) - SVBS)  
                  IF( TMPVAL > 0.5 )THEN
                    ICORDRYD(ND) = 1  
                  ELSE
                    TMPVAL = ABS(SVB(LN) - SVBN)  
                    IF( TMPVAL > 0.5)THEN
                      ICORDRYD(ND) = 1  
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ELSE
              ! *** CASE: HP < HDRY BUT RISING, JUST NOT FAST ENOUGH
              IF( ISCDRY(L) == 0 )THEN
                ISCDRY(L) = 1  
                ICORDRYD(ND) = 1  
              ENDIF  
            ENDIF     ! *** END OF REWETTING SECTION, DHPDT > DHPDT2
            
          ELSEIF( HP(L) < HDRY90 .OR. H1P(L) < HDRY )THEN  
            ! *** HP < HDRY.  SET SWITCHES TO DRY
            SUB(L)  = 0.0
            SUB(LE) = 0.0  
            SVB(L)  = 0.0  
            SVB(LN) = 0.0  
            SBX(L)  = 0.0  
            SBX(LE) = 0.0  
            SBY(L)  = 0.0  
            SBY(LN) = 0.0  
            IF( ISCDRY(L) == 0 )THEN  
              ISCDRY(L)=1  
              ICORDRYD(ND) = 1  
            ENDIF  
          ENDIF  
        ENDIF  
        
      ENDDO  
    ENDDO    ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
  End if !***End if on IDRY == 99

  ! *** PERFORM UPDATE OF P
  DO L=2,LA
    P(L) = G*(HP(L) + BELV(L))
  ENDDO

  ICORDRY = SUM(ICORDRYD(:))   ! *** OMP flag
  
#ifdef _MPI
  ! *** Now gather sum of ICCORDRY
  TTDS = DSTIME(0)
  Call MPI_barrier(MPI_Comm_World, ierr)
  TTWAIT = TTWAIT + (DSTIME(0)- TTDS)
  
  TTDS = DSTIME(0)
  Call communicate_ghost_cells(P, 'P')
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)

  Call DSI_All_Reduce(ICORDRY, ICORDRY_Global, MPI_SUM, TTDS, 1, TWAIT)
  DSITIMING(11) = DSITIMING(11) + TTDS
  TTWAIT = TTWAIT + TWAIT

#else
  ICORDRY_Global = ICORDRY
#endif

  ! *** IF ANY NODE OR SUB-DOMAIN WET/DRY STATUS CHANGED THEN UPDATE VARIABLES AND REPEAT PRESSURE SOLUTION
  IF( ICORDRY_Global > 0 )THEN

    NCORDRY = NCORDRY + 1
    
    IF( NCORDRY < 500  )THEN
      GOTO 1000
    ELSE
      ! *** WET/DRY MAXIMUM ITERATIONS EXCEEDED   
      WRITE (6,*)'THE LATEST MODEL RESULTS HAVE BEEN SAVED TO THE EE LINKAGE'

      ! *** SAVE A SNAPSHOT FOR EE
      IF( ISPPH == 1 )THEN
        Call Map_Write_EE_Binary
        if( process_id == master_id )THEN
          CALL EE_LINKAGE(-1)  
        endif
      ENDIF
        
      CALL STOPP('*** NCORDRY > 500 WHEN ISDRY > 0')
      
    ENDIF
    
  ELSEIF( ISDRY > 0 )THEN
    
#ifdef _MPI
    ! ***  IF SUB/SVB CHANGED THEN COMMUNICATE TO OTHER NODES
    TTDS = DSTIME(0)
    Call MPI_barrier(MPI_Comm_World, ierr)
    TTWAIT = TTWAIT + (DSTIME(0)- TTDS)
    
    TTDS = DSTIME(0)
    Call Communicate_PUV1
    TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
#endif
  ENDIF    ! *** END OF CHECKING FOR ISDRY ITERATIONS

  NCORDRYMAX = MAX(NCORDRY,NCORDRYMAX)
  NCORDRYAVG = NCORDRYAVG + NCORDRY
  NCOUNT = NCOUNT + 1
  
  ! *** *******************************************************************C
  ! *** FINISHED WITH WETTING/DRYING ITERATIONS
  IF( INT(TIMEDAY) /= DAYOLD .OR. TIMEDAY >= (TIMEEND-DELT/86400.) )THEN
  
    ! *** Write CALPUV.LOG to every process log
    OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
    
    IF( INT(TBEGIN) == DAYOLD .OR. INT(TIMEDAY) >= DAYOLD30 )THEN
      WRITE(mpi_log_unit,'(A)') '              N     NITER     TIMEDAY             NPUVAVG   NPUVMAX            NCONGAVG  NCONGMAX' 
      DAYOLD30 = DAYOLD + 30.
    ENDIF
    WRITE(mpi_log_unit, '(I15,I10,F12.3,3(10X,2I10))') N, NITER, TIMEDAY, NINT(FLOAT(NCORDRYAVG)/FLOAT(NCOUNT)), NCORDRYMAX, &
                                                                          NINT(FLOAT(ITERAVG)/FLOAT(NITERAVG)),  ITERMAX
    CLOSE(mpi_log_unit)
    
    NCORDRYAVG = 0
    ITERMAX = 0
    ITERAVG = 0
    NITERAVG = 0
    NCOUNT = 0
    DAYOLD = INT(TIMEDAY)
  ENDIF
  
  ! *** REPORT CELLS THAT HAVE JUST REWETTED
  IF( process_id == master_id .AND. ISDRY > 0 .AND. DEBUG )THEN
    OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
    DO L=2,LA
      IF( .NOT. LMASKDRY(L) )THEN 
        IF( HP(L) >= HDRY )THEN
          ! *** PREVIOUSLY CELL WAS DRY, NOW IT IS WET
          WRITE(8,'(A,2I5,F12.5,I5,L5,3F10.5)') 'REWETTED CELL: ',N,2,TIMEDAY,Map2Global(L).LG,LMASKDRY(L),H1P(L),HP(L)
        ENDIF
      ENDIF
    ENDDO
    CLOSE(8)
  ENDIF

  ! *** CHECK IF ISOLATED CELL WATER NEEDS TO BE REMOVED
  IUPDATE = 0
  IF( ISDRY > 0 .AND. NDRYSTP > 0 )THEN
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LE,LW,NTMP,IUW,IUE,IVS,IVN,IFACE,K)  &
    !$OMP                             PRIVATE(RDRY,BELVAVG,RVAL,HOLDTMP,TMPVAL,SVPW1,RNPORI,ETGWTMP,ETGWAVL) REDUCTION(MAX:IUPDATE)
    DO ND=1,NOPTIMAL
      LF = 2+(ND-1)*LDMOPT
      LL = MIN(LF+LDMOPT-1,LA)

      ! *** COUNT THE NUMBER TO TIME STEPS A CELL IS ISOLATED, AND IF IT HAS BEEN
      ! *** ISOLATED FOR MORE THAN NDRYSTP, AND ITS BOTTOM ELEVATION IS HIGHER
      ! *** THAN THE SURROUNDING DRY CELLS, THEN REDUCE ITS DEPTH BELOW THE
      ! *** DRYING DEPTH IF NECESSARY.  SAVE VOLUME REDUCTION RATE AS QDWASTE
      ! *** DEFINED AS POSITIVE OUT.
      NTMP = NDRYSTP
      DO L=LF,LL
        IF( HP(L) >= HDRY )THEN
          ! *** WET CELL, DETERMINE IF ISOLATED
          LE = LEC(L)
          LN = LNC(L)
          RDRY = SUB(L) + SUB(LE) + SVB(L) + SVB(LN)
          IF( RDRY > 0.5 )THEN
            NATDRY(L) = 0                         ! *** CELL IS NOT ISOLATED
          ELSE
            NATDRY(L) = NATDRY(L) + 1             ! *** CELL IS ISOLATED

            IF( NATDRY(L) > NTMP )THEN
              ! *** EXCEEDED THE NUMBER OF ISOLATED STEPS SO DETERMINE IF NEED TO DRY OUT THE CELL
              LW = LWC(L)
              LS = LSC(L)
              BELVAVG = 0.0
              RVAL = 0.0
              IF( SUBO(LE) > 0.5 .AND. BELV(LE) < BELV(L) )THEN
                RVAL = RVAL + 1.
              ENDIF
              IF( SUBO(L)  > 0.5 .AND. BELV(LW) < BELV(L) )THEN
                RVAL = RVAL + 1.
              ENDIF
              IF( SVBO(LN) > 0.5 .AND. BELV(LN) < BELV(L) )THEN
                RVAL = RVAL + 1.
              ENDIF
              IF( SVBO(L)  > 0.5 .AND. BELV(LS) < BELV(L) )THEN
                RVAL = RVAL + 1.
              ENDIF
              IF( RVAL > 0. .OR. HP(L) < HDRY*2. )THEN
                
                ! *** CHECK FOR OPEN BOUNDARY CELLS
                IF( ANY(LOBCS == L) )THEN
                  NATDRY(L) = 0
                  CYCLE
                ENDIF
                
                ! *** CHECK FOR FLOW BOUNDARY CELLS
                IF( QSUME(L) > 0.0 )THEN
                  NATDRY(L) = 0
                  CYCLE
                ENDIF
                
                ! *** SUM VOLUME OF "WASTED" WATER
                HOLDTMP = HP(L)
                HP(L)   = 0.90*HDRY
                P(L)    = G*(HP(L) + BELV(L))

                ! *** UPDATE H1P FOR THE DRY/WASTED CELL
                H1P(L)  = HP(L)
                DO K=KSZ(L),KC 
                  H1PK(L,K) = H1P(L)*DZC(L,K)
                ENDDO

                NATDRY(L) = 0
                QDWASTE(L) = DELTI*DXYP(L)*(HOLDTMP-HP(L))
                VDWASTE(L) = VDWASTE(L) + DXYP(L)*(HOLDTMP-HP(L))
                NNATDRY(ND) = NNATDRY(ND)+1
                LNATDRY(ND,NNATDRY(ND)) = L
                IUPDATE = 1
              ENDIF
            ENDIF
          ENDIF
        ENDIF  ! *** END OF WET CELL TEST
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  END IF

#ifdef _MPI 
  Call DSI_All_Reduce(IUPDATE, IUPDATE_Global, MPI_SUM, TTDS, 1, TWAIT)
  DSITIMING(11) = DSITIMING(11) + TTDS
  TTWAIT = TTWAIT + TWAIT
  
  IF( IUPDATE_Global > 0 )THEN
    TTDS = DSTIME(0)
    Call communicate_1D2(P, HP)
    TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  ENDIF
#endif

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LE,LW,NTMP,IUW,IUE,IVS,IVN,IFACE,K)  &
  !$OMP                             PRIVATE(RDRY,BELVAVG,RVAL,HOLDTMP,TMPVAL,SVPW1,RNPORI,ETGWTMP,ETGWAVL)
  DO ND=1,NOPTIMAL
    LF = 2+(ND-1)*LDMOPT
    LL = MIN(LF+LDMOPT-1,LA)
    
    ! *** *******************************************************************C
    ! *** PERFORM FINAL UPDATES OF HU, AND HV FOR SIGMA STRETCH GRID
    IF( IGRIDV == 0 )THEN
      !$OMP SIMD PRIVATE(LS,LW)
      DO L=LF,LL  
        LS = LSC(L)
        LW = LWC(L)
        HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )*FSGZUDXYPI(L)
        HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )*FSGZVDXYPI(L)
      ENDDO  

      !DIR$ NOFUSION
      DO L=LF,LL  
        HPI(L) = 1./HP(L)  
        HUI(L) = 1./HU(L)  
        HVI(L) = 1./HV(L)  
      ENDDO  
    ENDIF

    ! *** SET TRANSPORT MASK FOR DRY CELLS  
    IF( ISDRY > 0 )THEN  
      DO L=LF,LL  
        LMASKDRY(L)=.TRUE.  
      END DO  
      DO L=LF,LL  
        ! *** Bypass dry cells unless they are actively have boundary flows
        IF( HP(L) < HDRY )THEN
          IF( QSUME(L) /= 0.0 ) CYCLE        ! *** Check if cell is an active boundary
          LE=LEC(L)
          LN=LNC(L) 
          IUW=0  
          IUE=0  
          IVS=0  
          IVN=0  
          ! *** THIS REQUIRES THE CELL HAVE HP<HDRY FOR 2 ITERATIONS
          IF( SUB1(L)   < 0.5 .AND. SUB(L)   < 0.5 ) IUE = 1  
          IF( SUB1(LE)  < 0.5 .AND. SUB(LE)  < 0.5 ) IUW = 1  
          IF( SVB1(L)   < 0.5 .AND. SVB(L)   < 0.5 ) IVS = 1  
          IF( SVB1(LN)  < 0.5 .AND. SVB(LN)  < 0.5 ) IVN = 1  
          IFACE = IUW + IUE + IVS + IVN  
          IF( IFACE == 4 )THEN  
            LMASKDRY(L) = .FALSE.  
          END IF  
        ENDIF
      END DO  
    END IF  

    ! *** PERFORM UPDATE ON GROUNDWATER ELEVATION  
    IF( ISGWIE >= 1 )THEN
      !$OMP SIMD
      DO L=LF,LL  
        QSUM(L,KC)     = QSUM(L,KC)     - EVAPSW(L)  
        QSUM(L,KSZ(L)) = QSUM(L,KSZ(L)) + QGW(L)   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
      ENDDO  

      ! *** INFILTRATION STEP  
      RNPORI=1./RNPOR  
      IF( ISTL == 3 )THEN  
        DO L=LF,LL  
          AGWELV(L) = AGWELV2(L) - RNPORI*DELT*DXYIP(L)*QGW(L)  
        ENDDO  
      ELSE  
        DO L=LF,LL  
          AGWELV(L) = AGWELV1(L) - RNPORI*DELT*DXYIP(L)*QGW(L)  
        ENDDO  
      ENDIF  
      DO L=LF,LL  
        AGWELV(L)=MIN(AGWELV(L),BELV(L))  
      ENDDO  

      ! *** ET STEP  
      DO L=LF,LL  
        IF( IEVAP > 1 )THEN  
          SVPW1 = (10.**((0.7859+0.03477*TEM(L,KC))/(1.+0.00412*TEM(L,KC))))  
          EVAPT(L) = CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW1-VPAT(L))/PATMT(L)  
        ENDIF  
        ETGWTMP = EVAPT(L)-EVAPSW(L)*DXYIP(L)        ! *** EXCESS EVAPORATION
        ETGWTMP = MAX(ETGWTMP,0.0)  
        ETGWAVL = RNPOR*DELTI*(AGWELV(L)-BELAGW(L))  
        ETGWAVL = MAX(ETGWAVL,0.0)  
        ETGWTMP = MIN(ETGWTMP,ETGWAVL)  
        EVAPGW(L) = ETGWTMP*DXYP(L)                  ! *** TRANSPIRATION
      ENDDO  
      DO L=LF,LL  
        AGWELV(L)=AGWELV(L)-RNPORI*DELT*DXYIP(L)*EVAPGW(L)  
      ENDDO  
      DO L=LF,LL  
        AGWELV(L)=MAX(AGWELV(L),BELAGW(L))  
      ENDDO  
    ENDIF  

    IF( ISDRY > 0 )THEN
      DO K=1,KC
        !$OMP SIMD
        DO L=LF,LL
          SUB3D(L,K) = SUB(L)*SUB3DO(L,K)
          SVB3D(L,K) = SVB(L)*SVB3DO(L,K)
        ENDDO
      ENDDO
    ENDIF

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO
  
  ! *** IF ANY CELLS WASTED WATER THEN REPORT
  NTMP = SUM(NNATDRY)
  IF( NTMP > 0 )THEN
    IF( INOTICE /= INT(TIMEDAY) )THEN    ! *** ONLY DISPLAY WARNING ONCE PER DAY
      PRINT '(A,F14.5,I6,I4)',' QDWASTED CELLS. SEE log_mpi_proc_xxx.log FOR DETAILS. [TIMEDAY, # OF CELLS, Process ID]: ',TIMEDAY,NTMP,process_id
      INOTICE = INT(TIMEDAY)
    ENDIF

    OPEN(mpi_log_unit,FILE=OUTDIR//mpi_log_file,POSITION='APPEND')
    DO ND=1,NOPTIMAL
      DO I=1,NNATDRY(ND)
        L = LNATDRY(ND,I)
        
        TMPVAL = QDWASTE(L)/DXYP(L)
        WRITE(mpi_log_unit,8888) TIMEDAY, ND, Map2Global(L).IG, Map2Global(L).JG, TIMEDAY, H1P(L), HP(L), QDWASTE(L), TMPVAL
        
        ! *** ZERO QDWASTE IF NOT CONDUCTING MASS BALANCE
        IF( ISBAL == 0 ) QDWASTE(L) = 0.0
      ENDDO
    ENDDO
    CLOSE(mpi_log_unit)
    8888 FORMAT(' QDWASTE ',F12.5,3I6,F12.4,2F10.4,E14.6,F10.4)
  ENDIF

  ! *** CHECK FOR NEGATIVE DEPTHS  
  CALL NEGDEP(NOPTIMAL, LDMOPT, QCHANUT, QCHANVT, 2, SUB1, SVB1, HPOLD, NNEGFLG)  
  
  ! *** CALCULATE THE EXTERNAL DIVERGENCE  
  IF( ISDIVEX > 0 )THEN  
    DIVEXMX=0.  
    DIVEXMN=1000000.  
    DO L=2,LA  
      IF( SPB(L) /= 0 )THEN  
        LN=LNC(L)  
        DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI+0.5*(UHDYE(LEC(L))+UHDY1E(LEC(L))-UHDYE(L)-UHDY1E(L)+VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L)-QGW(L)+EVAPSW(L))  
        IF( ISDIVEX == 2 ) DIVEX = DIVEX*DXYIP(L)*HPI(L)  !  *** RELATIVE DIVERGENCE
        IF( DIVEX > DIVEXMX )THEN  
          DIVEXMX=DIVEX  
          LMAX=L  
        ENDIF  
        IF( DIVEX < DIVEXMN )THEN  
          DIVEXMN=DIVEX  
          LMIN=L  
        ENDIF  
      ENDIF  
    ENDDO  
    IMAX=IL(LMAX)  
    JMAX=JL(LMAX)  
    IMIN=IL(LMIN)  
    JMIN=JL(LMIN)  
    WRITE(6,6628)DIVEXMX,IMAX,JMAX  
    WRITE(6,6629)DIVEXMN,IMIN,JMIN  
    6628 FORMAT('  DIVEXMX=',E13.5,5X,2I10)  
    6629 FORMAT('  DIVEXMN=',E13.5,5X,2I10)  
  ENDIF  

  ! *** DETERMINE THE WET AND DRY CELL LIST
  IF( ISDRY == 0 )THEN
    LAWET = LA-1
    LADRY = 0
  ELSE
    LAWET = 0
    LADRY = 0
    DO L=2,LA
      IF( LMASKDRY(L) )THEN
        LAWET = LAWET+1
        LWET(LAWET) = L
      ELSEIF( OLDMASK(L) .OR. N < NTSTBC+1 )THEN
        ! *** ONLY FLAG NEWLY DRY CELLS
        LADRY = LADRY + 1
        LDRY(LADRY) = L
        
        ! *** UPDATE DEPTHS FOR DRY CELLS
        DO K=KSZ(L),KC
          HPK(L,K)  = HP(L)*DZC(L,K)
          HPKI(L,K) = 1./HPK(L,K)
        ENDDO    
        HU(L) = ( DXYP(L)*HP(L) + DXYP(LWC(L))*HP(LWC(L)) )*FSGZUDXYPI(L)
        HV(L) = ( DXYP(L)*HP(L) + DXYP(LSC(L))*HP(LSC(L)) )*FSGZVDXYPI(L)
      ENDIF
    ENDDO
  ENDIF
    
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    IF( ISBEDMAP > 0 )THEN
      LASED = 0
      LSED = 0
      DO L=1,LAWET
        IF( BEDMAP(LWET(L)) > 0 )THEN
          LASED = LASED + 1
          LSED(LASED) = LWET(L)
        ENDIF
      ENDDO
    ELSE
      LASED = LAWET
      IF( ISDRY > 0 .OR. NITER < 5 ) LSED = LWET
    ENDIF
  ENDIF

  LDMWET = INT(FLOAT(LAWET)/FLOAT(NDM))+1
  LDMDRY = INT(FLOAT(LADRY)/FLOAT(NDM))+1
  LDMSED = INT(FLOAT(LASED)/FLOAT(NDM))+1

  ! *** GET CELL LIST FOR ACTIVE LAYERS = 1
  IF( IGRIDV > 0 .AND. KMINV == 1 )THEN
    LASGZ1 = 0
    DO L=2,LA
      IF( KSZ(L) == KC )THEN
        IF( LMASKDRY(L) )THEN
          LASGZ1 = LASGZ1+1
          LSGZ1(LASGZ1)=L
        ENDIF
      ENDIF
    ENDDO
    LDMSGZ1 = INT(LASGZ1/NDM)+1
  ENDIF
  
  ! *** COMPUTATIONAL CELL LIST FOR EACH SUB-DOMAIN
  IF( ISDRY > 0 )THEN
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LDMWET,LAWET,LWET,LKSZ,LLWET,LKWET,KS,KC) PRIVATE(ND,LF,LL,K,LN,LP,L)
    DO ND=1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)
    
      DO K=1,KC
        LN = 0
        DO LP=LF,LL
          L = LWET(LP)  
          IF( LKSZ(L,K) )CYCLE
          LN = LN + 1
          LKWET(LN,K,ND)=L
        ENDDO
        LLWET(K,ND)=LN    ! *** NUMBER OF WET CELLS FOR THE CURRENT LAYER
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LLWET,LLWETZ,LKWET,LKWETZ,KS,KC) PRIVATE(ND,K,LP,L)
    DO ND=1,NDM
      DO K=1,KS
        LLWETZ(K,ND) = LLWET(K,ND)
        DO LP=1,LLWET(K,ND)
          LKWETZ(LP,K,ND) = LKWET(LP,K,ND)  
        ENDDO 
      ENDDO
      
      LLWETZ(KC,ND) = LLWET(KS,ND)
      DO LP=1,LLWET(KS,ND)
        LKWETZ(LP,KC,ND) = LKWET(LP,KS,ND)  
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  ENDIF

  ! *** THIRD PASS CELL CONSTANTS
  IF( IGRIDV > 0 )THEN
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NOPTIMAL,LDMOPT,LA,KC,LWC,LSC,KSZ,DZC,HP,HPK,HPKI,HU,HV,HPI,HUI,HVI)  &
    !$OMP                           SHARED(DXYP,FSGZUDXYPI,FSGZVDXYPI,LLWET,LKWET)                               &
    !$OMP                           PRIVATE(ND,LF,LL,K,LP,L,LW,LS) 
    DO ND=1,NOPTIMAL
      LF=2+(ND-1)*LDMOPT
      LL=MIN(LF+LDMOPT-1,LA)

      ! *** GLOBAL UPDATE OF LAYER THICKNESSES
      DO K=1,KC
        !$OMP SIMD PRIVATE(L)
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          HPK(L,K)  = HP(L)*DZC(L,K)
          HPKI(L,K) = 1./HPK(L,K)
        ENDDO
      ENDDO

      ! *** UPDATE HU & HV FOR ALL CELLS
      DO L=LF,LL
        LW = LWC(L)
        LS = LSC(L)

        IF( KSZ(LW) > KSZ(L) )THEN
          HU(L) = MAX( 0.5*HPK(L,KSZ(LW)), HP(LW)*(1.+DZC(L,KSZ(LW))*0.1) )
        ELSEIF( KSZ(LW) < KSZ(L) )THEN
          HU(L) = MAX( 0.5*HPK(LW,KSZ(L)), HP(L)*(1.+DZC(LW,KSZ(L))*0.1) )
        ELSE
          HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )*FSGZUDXYPI(L)
        ENDIF

        IF( KSZ(LS) > KSZ(L) )THEN
          HV(L) = MAX( 0.5*HPK(L,KSZ(LS)), HP(LS)*(1.+DZC(L,KSZ(LS))*0.1) )
        ELSEIF( KSZ(LS) < KSZ(L) )THEN
          HV(L) = MAX( 0.5*HPK(LS,KSZ(L)), HP(L)*(1.+DZC(LS,KSZ(L))*0.1) )
        ELSE
          HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )*FSGZVDXYPI(L)
        ENDIF
      ENDDO
  
      !DIR$ NOFUSION
      DO L=LF,LL
        HPI(L) = 1./HP(L)  
        HUI(L) = 1./HU(L)  
        HVI(L) = 1./HV(L)  
      ENDDO  

    ENDDO
    !$OMP END PARALLEL DO
  ENDIF

  ! *** BLOCKED LAYER FACE OPTION
  IF( NBLOCKED > 0 )THEN
    DO LP=1,NBLOCKED
      L = LBLOCKED(LP)
      IF( KSZ(L) == KC ) CYCLE
      
      ! *** U FACE BLOCKING
      IF( SUBO(L) > 0. .AND. (BLDRAFTUO(LP)+BLSILLU(LP)) > 0.0 )THEN
        IF( BLANCHORU(LP) /= 0. )THEN
          BLDRAFTU(LP) = BLDRAFTUO(LP) - (BLANCHORU(LP) - (BELV0(L) + HP(L)))
          BLDRAFTU(LP) = MAX(BLDRAFTU(LP),0.0)
        ELSE
          BLDRAFTU(LP) = BLDRAFTUO(LP) 
        ENDIF

        HU(L) = HU(L) - BLDRAFTU(LP) - BLSILLU(LP)
        HU(L) = MAX(HU(L),HWET)
        HUI(L) = 1./HU(L)
      
        ! *** RESET TO DEFAULT LAYER
        KM=MAX(KSZ(LWC(L)), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR U FACE
        DO K=1,KC
          SUB3D(L,K) = 0.0
          IF( K >= KM )THEN
            SUB3D(L,K) = 1.0
            IF( IGRIDV > 0 )THEN
              IF( KSZ(LWC(L)) > KSZ(L) )THEN
                SGZU(L,K)  = DZC(LSC(L),K)
              ELSE
                SGZU(L,K)  = DZC(L,K)
              ENDIF
            ELSE
              SGZU(L,K)  = MAX(DZC(LWC(L),K),DZC(L,K))
            ENDIF
          ENDIF
        ENDDO
      
        ! *** BLOCK CELL FACES FROM TOP
        IF( BLDRAFTU(LP) > 0.0 )THEN
          TMPX = 0.
          DO K=KC,KSZU(L)+1,-1
            TMPX = TMPX + SGZU(L,K)*HP(L)
            KTBU(LP) = K - 1
            SUB3D(L,K) = 0.0
            SGZU(L,K) = 0.0
            IF( TMPX > BLDRAFTU(LP) ) EXIT
          ENDDO
          KTBU(LP) = MAX(KTBU(LP),KSZU(L))
        ELSE
          KTBU(LP) = KC
        ENDIF
        
        ! *** BLOCK CELL FACES FROM BOTTOM
        IF( BLSILLU(LP) > 0.0 )THEN
          TMPX = 0.
          DO K=KSZU(L),KC-1
            TMPX = TMPX + SGZU(L,K)*HP(L)
            KBBU(LP) = K + 1
            SUB3D(L,K) = 0.0
            SGZU(L,K) = 0.0
            IF( TMPX > BLSILLU(LP) ) EXIT
          ENDDO
          KBBU(LP) = MIN(KBBU(LP),KTBU(LP))
        ELSE
          KBBU(LP) = KSZU(L)
        ENDIF
        
        TMPX = SUM(SGZU(L,1:KC))
        DO K=1,KC
          SGZU(L,K) = SGZU(L,K) / TMPX
          DZGU(L,K) = 0.0
          CDZFU(L,K) = 0.0          ! *** USED FOR GLOBAL DU SHEAR
          CDZUU(L,K) = 0.0          ! *** USED FOR SURFACE SHEAR DUE TO WIND
          CDZLU(L,K) = 0.0
          CDZMU(L,K) = 0.0
          CDZRU(L,K) = 0.0
          CDZDU(L,K) = 0.0
        ENDDO
      
        DO K=KBBU(LP),KTBU(LP)-1
          DZGU(L,K) = 0.5*(SGZU(L,K)+SGZU(L,K+1))
          CDZFU(L,K) = SGZU(L,K)*SGZU(L,K+1)/(SGZU(L,K) + SGZU(L,K+1))
          CDZUU(L,K) = -SGZU(L,K)  /(SGZU(L,K) + SGZU(L,K+1))
          CDZLU(L,K) = -SGZU(L,K+1)/(SGZU(L,K) + SGZU(L,K+1))
        ENDDO
        
        CDZRU(L,KBBU(LP)) = SGZU(L,KBBU(LP)) - 1.
        CDZDU(L,KBBU(LP)) = SGZU(L,KBBU(LP))
        DO K=KBBU(LP)+1,KTBU(LP)-1
          CDZRU(L,K) = CDZRU(L,K-1) + SGZU(L,K)
          CDZDU(L,K) = CDZDU(L,K-1) + SGZU(L,K)
        ENDDO
        
        DO K=KBBU(LP),KTBU(LP)-1
          CDZRU(L,K) = CDZRU(L,K)*DZGU(L,K)*CDZLU(L,KSZU(L))
          CDZMU(L,K) = 0.5*SGZU(L,K)*SGZU(L,K+1)
        ENDDO
      ENDIF
      
      ! *** V FACE BLOCKING
      IF( SVBO(L) > 0. .AND. (BLDRAFTVO(LP)+BLSILLV(LP)) > 0.0 )THEN
        IF( BLANCHORV(LP) /= 0. )THEN
          BLDRAFTV(LP) = BLDRAFTVO(LP) - (BLANCHORV(LP) - (BELV0(L) + HP(L)))
          BLDRAFTV(LP) = MAX(BLDRAFTV(LP),0.0)
        ELSE
          BLDRAFTV(LP) = BLDRAFTVO(LP) 
        ENDIF

        HV(L) = HV(L) - BLDRAFTV(LP) - BLSILLV(LP)
        HV(L) = MAX(HV(L),HWET)
        HVI(L) = 1./HV(L)
      
        ! *** RESET TO DEFAULT LAYER
        KM=MAX(KSZ(LSC(L)), KSZ(L))      ! *** MINUMUM ACTIVE LAYERS FOR V FACE
        DO K=1,KC
          SVB3D(L,K) = 0.0
          IF( K >= KM )THEN
            SVB3D(L,K) = 1.0
            IF( IGRIDV > 0 )THEN
              IF( KSZ(LSC(L)) > KSZ(L) )THEN
                SGZV(L,K)  = DZC(LSC(L),K)
              ELSE
                SGZV(L,K)  = DZC(L,K)
              ENDIF
            ELSE
              SGZV(L,K)  = MAX(DZC(LSC(L),K),DZC(L,K))
            ENDIF
          ENDIF
        ENDDO
      
        ! *** BLOCK CELL FACES
        IF( BLDRAFTV(LP) > 0.0 )THEN
          TMPY = 0.
          DO K=KC,KSZV(L)+1,-1
            TMPY = TMPY + SGZV(L,K)*HP(L)
            KTBV(LP) = K - 1
            SVB3D(L,K) = 0.0
            SGZV(L,K) = 0.0
            IF( TMPY > BLDRAFTV(LP) ) EXIT
          ENDDO
          KTBV(LP) = MAX(KTBV(LP),KSZV(L))
        ELSE
          KTBV(LP) = KC
        ENDIF
        
        ! *** BLOCK CELL FACES FROM BOTTOM
        IF( BLSILLV(LP) > 0.0 )THEN
          TMPX = 0.
          DO K=KSZV(L),KC-1
            TMPX = TMPX + SGZV(L,K)*HP(L)
            KBBV(LP) = K + 1
            SVB3D(L,K) = 0.0
            SGZV(L,K) = 0.0
            IF( TMPX > BLSILLV(LP) ) EXIT
          ENDDO
          KBBV(LP) = MIN(KBBV(LP),KTBV(LP))
        ELSE
          KBBV(LP) = KSZV(L)
        ENDIF
        
        TMPY = SUM(SGZV(L,1:KC))
        DO K=1,KC
          SGZV(L,K) = SGZV(L,K) / TMPY
          DZGV(L,K) = 0.0
          CDZFV(L,K) = 0.0
          CDZUV(L,K) = 0.0
          CDZLV(L,K) = 0.0
          CDZRV(L,K) = 0.0
          CDZDV(L,K) = 0.0
          CDZMV(L,K) = 0.0
        ENDDO
      
        DO K=KBBV(LP),KTBV(LP)-1
          DZGV(L,K) = 0.5*(SGZV(L,K)+SGZV(L,K+1))
          CDZFV(L,K) = SGZV(L,K)*SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
          CDZUV(L,K) = -SGZV(L,K)  /(SGZV(L,K)+SGZV(L,K+1))
          CDZLV(L,K) = -SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
        ENDDO

        CDZRV(L,KBBV(LP)) = SGZV(L,KBBV(LP))-1.
        CDZDV(L,KBBV(LP)) = SGZV(L,KBBV(LP))
        DO K=KBBV(LP)+1,KTBV(LP)-1
          CDZRV(L,K) = CDZRV(L,K-1)+SGZV(L,K)
          CDZDV(L,K) = CDZDV(L,K-1)+SGZV(L,K)
        ENDDO

        DO K=KBBV(LP),KTBV(LP)-1
          CDZRV(L,K) = CDZRV(L,K)*DZGV(L,K)*CDZLV(L,KSZV(L))
          CDZMV(L,K) = 0.5*SGZV(L,K)*SGZV(L,K+1)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  ! *** COMPUTATIONAL CELL LIST FOR ENTIRE DOMAIN
  IF( ISDRY > 0 )THEN
    DO K=1,KC
      LN=0
      DO L=2,LA
        IF( LKSZ(L,K) )CYCLE
        IF( LMASKDRY(L) )THEN
          LN = LN+1
          LKWET(LN,K,0) = L   ! *** Wet Cell for Layer K
        ENDIF
      ENDDO
      LLWET(K,0) = LN         ! *** Total Wet Cells for Layer K
    ENDDO
    
  ENDIF
  
#ifdef _MPI
  TTDS = DSTIME(0)
  Call MPI_barrier(MPI_Comm_World, ierr)
  TTWAIT = TTWAIT + (DSTIME(0)- TTDS)
  
  TTDS = DSTIME(0)
  Call Communicate_PUV3
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
#endif

  RETURN

END  

