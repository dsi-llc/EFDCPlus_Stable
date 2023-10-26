! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE WATERQUALITY
    
  ! *** WATER QUALITY MODEL KINETIC AND SEDIMENT PROCESSES

  ! *** ------------------------------------------------------------------------------
  ! *** CHANGE RECORD  
  ! *** DATE MODIFIED     BY               DESCRIPTION
  ! *** ------------------------------------------------------------------------------
  ! *** 2019-12           D.K. TRAN          Derived data type for algae classes
  ! *** 2019-11           D.K. TRAN          Converted to Module
  ! *** 2018-04           Paul M. Craig      Rearranged entire water quality code
  ! ***                                      to be better integrated with EFDC+
  ! *** 2011-03           Paul M. Craig      Rewritten to F90 
  ! *** 2008                                 Merged SNL and DSI

  Use GLOBAL    
  Use INFOMOD
  Use JULIANMOD
  
  USE Variables_MPI
  Use Broadcast_Routines
  Use Variables_MPI_Mapping
  Use Mod_Map_Global_to_Local
  Use Variables_MPI_Write_Out
  
  Use Variables_WQ
  USE WQ_DIAGENESIS
  USE WQ_RPEM_MODULE
  USE SHELLFISHMOD
  USE WQ_ZOOPLANKTON
  USE WQ_MACROPHYTE_FEEDBACK

  ! *** LB is a globally declared integer that can be hardwired for debugging (defined in AAEFDC)
  
  CONTAINS
  
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ3D
  !
  !> @details WATER QUALITY MODEL KINETIC AND SEDIMENT PROCESSES
  !---------------------------------------------------------------------------!
  ! 2018-04           Paul M. Craig      Rearranged entire water quality code
  !                                      to be better integrated with EFDC+
  ! 2011-03           Paul M. Craig      Rewritten to F90 
  ! 2008  
  !---------------------------------------------------------------------------! 
  SUBROUTINE WQ3D

  USE CALCSERMOD,ONLY: CALCSER

  IMPLICIT NONE

  REAL(RKD), STATIC :: DAYNEXT
  REAL(RKD), STATIC :: SUNDAY1, SUNDAY2
  REAL,      STATIC :: SUNSOL1, SUNSOL2
  REAL,      STATIC :: SUNFRC1, SUNFRC2
  REAL,      STATIC :: WQKCNT=0.
  
  REAL       :: TIMTMP, RATIO, SOLARAVG, WTEMP, WQTT, TT20, HP2I, WQ2
  INTEGER    :: IACTION ! not used
  INTEGER    :: IWQTAGR, IWQTSTL, ISMTICI ! not used
  INTEGER    :: M1, M2, L, K, NMALG, NW
  INTEGER    :: LF, LL, LP, NAL, ND
  INTEGER, STATIC :: M
  
  REAL(RKD), EXTERNAL :: DSTIME 
  REAL(RKD)           :: TTDS       ! MODEL TIMING TEMPORARY VARIABLE

  DATA IWQTAGR,IWQTSTL,ISMTICI/3*0/  
  
  ! *** SET THE HYDRODYNAMIC TIMESTEP
  IF( ISDYNSTP == 0 )THEN  
    DELT=DT  
  ELSE  
    DELT=DTDYN  
  ENDIF  

  ! *** INITIALIZE PARAMETERS ON FIRST CALL
  IF( ITNWQ == 0 )THEN
    ! *** INITIALIZE DAYNEXT VALUE
    DAYNEXT = DBLE(INT(TIMEDAY)) + 1.

    ! *** PMC - NEW IMPLEMENTATION TO USE DAILY (FROM HOURLY) SOLAR RADIATION FOR ALGAL GROWTH
    IF( NASER > 0 )THEN
      IF( IWQSUN == 3 )THEN
        ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
        SUNDAY1 = DAYNEXT - 1.
        SUNDAY2 = DAYNEXT

        ! *** FIND 1ST POINT
        M = 1
        DO WHILE (TSATM(1).TIM(M) < SUNDAY1)
          M = M+1
        END DO
    
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL1 = 0.0
        DO WHILE (TSATM(1).TIM(M) < SUNDAY1)
          M1 = M1+1
          IF( TSATM(1).VAL(M,6) > 0. )THEN
            M2 = M2+1
            SUNSOL1=SUNSOL1+TSATM(1).VAL(M,6)
          ENDIF
          M = M+1
        END DO
        IF( M1 > 0 )THEN
          SUNFRC1=FLOAT(M2)/FLOAT(M1)
          SUNSOL1=SUNSOL1/FLOAT(M1)
        ELSE
          SUNFRC1=1.0
        ENDIF
    
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL2 = 0.
        DO WHILE (TSATM(1).TIM(M) < SUNDAY2)
          M1 = M1+1
          IF( TSATM(1).VAL(M,6) > 0. )THEN
            M2 = M2+1
            SUNSOL2=SUNSOL2+TSATM(1).VAL(M,6)
          ENDIF
          M = M+1
        END DO
        IF( M1 > 0 )THEN
          SUNFRC2=FLOAT(M2)/FLOAT(M1)
          SUNSOL2=SUNSOL2/FLOAT(M1)
          IF( SUNSOL1 == 0.0 )THEN
            SUNFRC1=SUNFRC2
            SUNSOL1=SUNSOL2
          ENDIF            
        ELSE
          SUNFRC2=1.
        ENDIF
      ENDIF
    ENDIF    ! *** NASER > 0
  ENDIF      ! *** End of initialization

  ! *** WQI1 = SOLAR RADIATION ON PREVIOUS DAY  
  ! *** WQI2 = SOLAR RADIATION TWO DAYS AGO  
  ! *** WQI3 = SOLAR RADIATION THREE DAYS AGO  
  ! *** Update occurs only when the simulation day changes.  
  IF( TIMEDAY > DAYNEXT )THEN  
    WQI3 = WQI2  
    WQI2 = WQI1  
    
    ! *** Check if no solar radiation the previous day. Only update WQI1 if SUNFRC2 > 0
    IF( IWQSUN == 2 )THEN
      IF( SUNFRC2 > 0.0 )THEN
        WQI1 = WQI0OPT/SUNFRC2  
      ENDIF
      WQI0OPT = 0.0            ! *** Reset daily average
      SUNFRC2 = 0.             ! *** Reset sun fraction
    ELSE
      WQI1 = WQI0
    ENDIF
    DAYNEXT = DAYNEXT + 1.
  ENDIF

  ! *** INCREMENT THE WATER QUALITY TIMESTEP
  WQKCNT = WQKCNT + DELT/86400.
  
  ! *** APPLY WATER DEPTH UPDATE TO WQ CONSTITUENTS THAT DO NOT TRANSPORT, I.E. MACROPHYTES
  IF( NFIXED > 0 )THEN
    DO NAL = 1, NALGAE
      IF( .NOT. ALGAES(NAL).ISMOBILE )THEN
        NW = 19 + NAL
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,HP2I,WQ2) 
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
        
          DO LP=LF,LL  
            L = LWET(LP)
            K = KSZ(L)
            HP2I = H1PK(L,K)
            IF( ISTL == 3 ) HP2I = H2PK(L,K)
            WQ2 = WQV(L,K,NW)*HP2I              ! *** Mass from previous timestep
            WQ2 = WQ2*HPKI(L,K)                 ! *** New concentration
            WQV(L,K,NW) = WQ2
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
    ENDDO
  ENDIF
   
  ! *** UPDATE WATER COLUMN KINETICS AND SEDIMENT MODEL  
  ! *** OVER LONGER TIME INTERVALS THAN PHYSICAL TRANSPORT  
  IF( ITNWQ == 0 .OR. WQKCNT >= WQKINUPT )THEN  
    DTWQ   = WQKCNT         ! *** KINETIC TIME STEP USED (DAYS)
    DTWQO2 = DTWQ*0.5  
    WQKCNT = 0.
    
    ! *** INTIALIZE THE TIMER FOR WQ KINETICS, ALL COMPONENTS
    TTDS = DSTIME(0) 

    ! *** GET SOLAR RADIATION INTENSITY AND DAYLIGHT LENGTH  
    ! *** NOTE: IWQSUN=1 CALLS SUBROUTINE WQSUN WHICH READS THE DAILY  
    ! ***                SOLAR RADIATION DATA FROM FILE SUNDAY.INP WHICH  
    ! ***                ARE IN UNITS OF LANGLEYS/DAY.  
    ! ***       IWQSUN=2 USES THE HOURLY SOLAR RADIATION DATA FROM ASER.INP  
    ! ***                COUPLED WITH THE COMPUTED OPTIMAL DAILY LIGHT TO
    ! ***                LIMIT ALGAL GROWTH.
    ! ***       IWQSUN=3 USES THE DAILY AVERAGE SOLAR RADIATION DATA COMPUTED 
    ! ***                FROM THE HOURLY ASER.INP AND THE COMPUTED OPTIMAL DAILY
    ! ***                LIGHT TO LIMIT ALGAL GROWTH.
    ! ***       IWQSUN>1 USES THE DAILY AVERAGE SOLAR RADIATION DATA COMPUTED 
    ! ***                FROM THE HOURLY ASER.INP DATA.  CONVERTS WATTS/M**2 TO
    ! ***                LANGLEYS/DAY USING 2.065.  COMPUTES THE FRACTION OF
    ! ***                DAYLIGHT AND ADJUSTS FOR PHOTOSYNTHETIC ACTIVE RADIATION BY 
    ! ***                PARADJ (~0.43) 

    IF( IWQSUN == 0 )THEN
      WQI1 = WQI0         ! *** Constant
      
    ELSEIF( IWQSUN == 1 )THEN  
      ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
      CALL WQSUN             ! *** Read SUNDAY.INP
      WQI0 = SOLSRDT  
      WQFD = SOLFRDT  
      
    ELSEIF( NASER > 0 )THEN
      ! *** SOLAR RADIAION COMES FROM ASER FILE.  IWQSUN: 2-USE TIMING FROM ASER, 3-DAILY AVERAGE COMPUTED FROM ASER
      IF( IWQSUN == 3 )THEN
        ! *** Check if new day
        IF( TIMEDAY > SUNDAY2 )THEN
          ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
          SUNDAY1 = SUNDAY2
          SUNSOL1 = SUNSOL2
          SUNFRC1 = SUNFRC2
        
          ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
          M1 = 0
          M2 = 0
          SUNSOL2 = 0.
          SUNDAY2 = SUNDAY2 + 1.
          DO WHILE (TSATM(1).TIM(M) < SUNDAY2)
            M1 = M1+1
            IF( TSATM(1).VAL(M,6) > 0. )THEN
              M2 = M2+1
              SUNSOL2=SUNSOL2+TSATM(1).VAL(M,6)
            ENDIF
            M = M+1
            IF( M > TSATM(1).NREC )THEN
              M = TSATM(1).NREC
              EXIT
            ENDIF
          END DO
          IF( M1 > 0 )THEN
            SUNFRC2=FLOAT(M2)/FLOAT(M1)
            SUNSOL2=SUNSOL2/FLOAT(M1)
          ELSE
            SUNFRC2=1.
          ENDIF
        ENDIF

        RATIO = (TIMEDAY-SUNDAY1)
        SOLARAVG = RATIO*(SUNSOL2-SUNSOL1)+SUNSOL1

        ! *** SOLAR RADIATION IN LANGLEYS/DAY
        WQI0 = PARADJ*2.065*SOLARAVG  
        WQFD = RATIO*(SUNFRC2-SUNFRC1)+SUNFRC1
    
      ELSEIF( IWQSUN == 2 )THEN
        ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
        IF( LDAYLIGHT .AND. (NASER > 1 .OR. USESHADE) )THEN  
          SOLARAVG = 0.                           ! *** SOLARAVG is the domain average solar radiation
          DO L=2,LA  
            SOLARAVG = SOLARAVG + SOLSWRT(L)      ! *** SOLSWRT already include surface albedo and shading
          ENDDO  
          SOLARAVG = SOLARAVG/FLOAT(LA-1)
        ELSE
          ! *** Spatially Constant Atmospheric Parameters
          SOLARAVG = SOLSWRT(2)
        ENDIF  
        ! *** SOLAR RADIATION IN LANGLEYS/DAY
        WQI0 = PARADJ*2.065*SOLARAVG               ! *** Current light for growth
        WQFD = 1.                                  ! *** Set fraction of day to 1 for 
        
        WQI0OPT = WQI0OPT + SOLARAVG*DTWQ          ! *** Sum current light for daily average
        IF( SOLARAVG > 0.0 ) SUNFRC2 = SUNFRC2 + DTWQ
      ENDIF
    ENDIF  

    ! *** MASS LOADING BC'S.  WQWPSL IS ONLY USED IN THE KINETIC ROUTINES.  
    ! *** IF CONCENTRATION BASED LOADING, CALCSER AND CALFQC ALREADY HANDLED LOADING
    IF( IWQPSL == 1 ) CALL WQPSL  

    ! *** CALL SPATIALLY AND TIME VARYING BENTHIC FLUX HERE, IF REQUIRED  
    ! *** IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.  
    IF( IWQBEN  ==  2 )THEN  
      IF( ISDYNSTP == 0 )THEN  
        TIMTMP=(DT*FLOAT(N)+TCON*TBEGIN)/86400.  
      ELSE  
        TIMTMP=TIMEDAY
      ENDIF  
      IF( TIMTMP  >=  BENDAY )THEN  
        CALL WQBENTHIC(TIMTMP)  
      ENDIF  
    ENDIF  

    ! *** UPDATE WET DEPOSITION (LRAIN ENSURES THAT THERE IS VALID RAINFALL)
    IF( LRAIN ) CALL WQWET  
    
    DO NW = 1, NWQV
      IF( ISKINETICS(NW) > 0 )THEN 
        DO K=1,KC  
          DO L=2,LA  
            WQVO(L,K,NW) = WQV(L,K,NW)
          ENDDO  
        ENDDO
      ENDIF  
    ENDDO  

    ! *** SET UP LOOK-UP TABLE FOR BACTERIA (FCB) TEMPERATURE DEPENDENCY OVER -10 degC TO 50 degC  
    IF( ISKINETICS(IFCB) > 0 )THEN
      WTEMP = WQTDMIN
      DO M1=1,NWQTD  
        TT20 = WTEMP - 20.0  
        WQTT = WQKFCB * WQTFCB**TT20 * DTWQO2  
        WQTD1FCB(M1) = 1.0 - WQTT  
        WQTD2FCB(M1) = 1.0 / (1.0 + WQTT)  
        WTEMP = WTEMP + WQTDINC
      ENDDO  
    ENDIF
    
    
    ! ***   CALCULATE KINETIC SOURCES AND SINKS  
    IF( ISWQLVL == 0 ) CALL WQSKE0  
    IF( ISWQLVL == 1 ) CALL WQSKE1 ! *** Extension of CEQUAL-ICM for unlimited algae + zooplankton
    IF( ISWQLVL == 2 ) CALL WQSKE2  
    IF( ISWQLVL == 3 ) CALL WQSKE3  
    IF( ISWQLVL == 4 ) CALL WQSKE4  
    
    ! *** ZOOPLANKTON
    ! *** It is called after WQSKE1 since the kinetic zones map is defined in WQSKE1
    IF( IWQZPL > 0 )THEN
      CALL ZOOPL_KINETIC
    ENDIF
    
    TWQKIN = TWQKIN + (DSTIME(0)-TTDS) 

    ! ***   DIAGNOSE NEGATIVE CONCENTRATIONS  
    IF( IWQNC > 0 ) CALL WWQNC  

    ! *** CALL SEDIMENT DIAGENSIS MODEL AFTER THE FIRST WQ ITERATION
    IF( IWQBEN == 1 .AND. ITNWQ > 0 )THEN
      TTDS = DSTIME(0) 
      CALL SMMBE  
      TWQSED = TWQSED + (DSTIME(0)-TTDS) 
    ENDIF  

    ! *** RPEM
    IF( ISRPEM > 0 .AND. ITNWQ > 0  )THEN
      TTDS = DSTIME(0) 
      CALL CAL_RPEM
      TWQRPEM = TWQRPEM + (DSTIME(0)-TTDS) 
    ENDIF

  ENDIF    ! *** ENDIF ON KINETIC AND SEDIMENT UPDATE  
  
  ! *** UPDATE WATER QUALITY TIMESTEP
  ITNWQ = ITNWQ + 1

  RETURN  
  
  END SUBROUTINE WQ3D
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ3DINP
  !
  !> @details  READ WATER QUALITY SUBMODEL INPUT FILES
  !---------------------------------------------------------------------------!
  !  ORGINALLY CODED BY K.-Y. PARK
  !  OPTIMIZED AND MODIFIED BY J. M. HAMRICK  
  !---------------------------------------------------------------------------! 
  SUBROUTINE WQ3DINP
  
  USE RESTART_MODULE, ONLY:WQ_WCRST_IN, WQSDRST_IN
  IMPLICIT NONE

  INTEGER :: J, L, LG, K, IWQTICI, IWQTAGR, IWQTSTL, IWQTSUN, IWQTBEN, IWQTPSL
  INTEGER :: IWQTNPL, ISMTICI, NWQVOUT, NS, NW, NAL
  REAL    :: O2WQ_, WQTAMD 

  CHARACTER*3 CWQHDR(NWQVM)
  DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
  DATA ISMTICI/0/

  IWQTICI=IWQTICI
  IWQTAGR=IWQTAGR
  IWQTSTL=IWQTSTL
  IWQTSUN=IWQTSUN
  IWQTBEN=IWQTBEN
  IWQTPSL=IWQTPSL
  IWQTNPL=IWQTNPL
  ISMTICI=ISMTICI
  
  if( process_id == master_id )then
    OPEN(2,FILE=OUTDIR//'WQ3D.OUT',STATUS='UNKNOWN')
    CLOSE(2,STATUS='DELETE')
    OPEN(2,FILE=OUTDIR//'WQ3D.OUT',STATUS='UNKNOWN')
  endif
  
  WQKINUPT = DT/86400.  ! WQ VARIABLE DT
  UHEQ(1)=0.0
  UHEQ(LC)=0.0
  DO L=2,LA
    UHEQ(L)=1.0
  ENDDO
  
  ITNWQ = 0
  ! *** Fractional of Depth at the Top of the Layer
  RKCWQ = 1.0/REAL(KC)
  DO K=1,KC
    WQHT(K)=REAL(KC-K)*RKCWQ
  ENDDO

  DO K=1,KC
    IWQPSC(1,K)=0
    IWQPSC(LC,K)=0
  ENDDO
  DO K=1,KC
    DO L=2,LA
      IWQPSC(L,K)=0
      IWQPSV(L,K)=0
    ENDDO
  ENDDO
  DO J=1,NWQV
    DO K=1,KC
      DO L=1,LC
        WQWDSL(L,K,J)=0.0
        WQWPSL(L,K,J)=0.0
      ENDDO
    ENDDO
  ENDDO
  
  ! *** Read main water quality control file
  CALL WQ3DCONTROL
  CALL ALGAECONTROL

  NWQVOUT=0
  if( process_id == master_id )then
    OPEN(1,FILE=OUTDIR//'WQWCTS.OUT',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'WQWCTS.OUT',STATUS='UNKNOWN')
    DO NW=1,NWQV
      IF( ISKINETICS(NW) == 1 )THEN
        NWQVOUT=NWQVOUT+1
        CWQHDR(NWQVOUT) = WQCONSTIT(NW)
      ENDIF
    ENDDO
    WRITE(1,1969)(CWQHDR(NW),NW=1,NWQVOUT)
     1969 FORMAT('C   I    J    K    TIME',7X,A3,8X,A3,8X,A3, &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3, &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3, &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3)
    CLOSE(1)
  endif
  
  ! **  INITIALIZE DIURNAL DO ANALYSIS
  IF( WQHRAVG > 0.0)THEN
    if( process_id == master_id )then
      OPEN(1,FILE=OUTDIR//'DIURNDO.OUT')
      CLOSE(1,STATUS='DELETE')
    endif
    DO K=1,KC
      DO L=2,LA
        DDOMAX(L,K)=-1.E6
        DDOMIN(L,K)=1.E6
      ENDDO
    ENDDO
  ENDIF

  ! **  INITIALIZE LIGHT EXTINCTION ANALYSIS
  NDLTCNT=0
  IF( NDLTAVG >= 1 )THEN
    if( process_id == master_id )then
      OPEN(1,FILE=OUTDIR//'LIGHT.OUT')
      CLOSE(1,STATUS='DELETE')
    endif
    DO K=1,KC
      DO L=2,LA
        RLIGHTT(L,K)=0.
        RLIGHTC(L,K)=0.
      ENDDO
    ENDDO
  ENDIF

  ! *** READ INITIAL CONDITIONS
  IF( IWQICI > 0 )THEN
    IF( process_id == master_id )THEN
      IF( IWQICI == 1 ) CALL WQICI 
      IF( IWQICI == 2 ) CALL WQ_WCRST_IN 
    ENDIF
    
    IF( IWQICI == 2 )THEN
      Call Broadcast_Scalar(WQI1,    master_id)
      Call Broadcast_Scalar(WQI2,    master_id)
      Call Broadcast_Scalar(WQI3,    master_id)
    ENDIF
    Call Broadcast_Array(WQV_Global, master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        DO NW=1,NWQV
          DO K=1,KC
            WQV(L,K,NW) = WQV_Global(LG,K,NW)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  ! *** WQCHLX=1/WQCHLX
  DO L=2,LA
    DO K=1,KC
      WQCHL(L,K) = 0
      DO NAL = 1, NALGAE
        IF( ALGAES(NAL).ISMOBILE )THEN
          WQCHL(L,K) = WQCHL(L,K) + WQV(L,K,19+NAL)*ALGAES(NAL).WQCHLA
        ENDIF
      ENDDO
      IF( IWQSRP == 1 )THEN
        O2WQ_ = MAX(WQV(L,K,IDOX), 0.0)
        WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ_), WQV(L,K,ITAM) )
        WQTAMP(L,K) = WQV(L,K,ITAM) - WQTAMD
        WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*WQTAMP(L,K))
        WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*WQTAMP(L,K))
      ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN
        WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*SEDT(L,K))
        WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*SEDT(L,K))
      ELSE
        WQPO4D(L,K) = WQV(L,K,IP4D)
        WQSAD(L,K)  = WQV(L,K,ISAA)
      ENDIF
    ENDDO
  ENDDO
  
  ! *** Set vegetative growth and drag 
  DO NAL = 1, NALGAE
    IF( .NOT. ALGAES(NAL).ISMOBILE )THEN
      ! *** Set base layer for all cells
      DO L=2,LA
        ! *** Set base layer               delme - TODO - handle growth downward from suspended base
        DO K = KSZ(L),KC
          IF( HP(L)*Z(L,K-1) >= (HP(L)-ALGAES(NAL).BASEDEPTH) )THEN
            LAYERBOT(NAL,L) = K                ! *** Bottom active layer
            EXIT                               ! *** Jump out of the layer loop
          ENDIF
        ENDDO
      ENDDO
      
      IF( ALGAES(NAL).ISDRAG > 0 )THEN
        ISVEG = 2                                               ! *** Turn on the vegetation drag/turbulence model
      ENDIF      
      IF( ALGAES(NAL).THRESHOLD /= 0 )THEN
        Call Macro_Veg(NAL)
      ELSE
        ! *** Use specified heights to set vegetation height
        HEIGHT_MAC(L,NAL) = ALGAES(NAL).MAXLENGTH
        
        ! *** Set base layer for all cells
        DO L=2,LA
          ! *** Set top layer               delme - TODO - handle growth downward from suspended base
          DO K = KSZ(L),KC
            IF( HP(L)*Z(L,K-1) >= (HP(L) - ALGAES(NAL).BASEDEPTH + HEIGHT_MAC(L,NAL)) )THEN
              LAYERTOP(NAL,L) = K                ! *** Top active layer
              EXIT                               ! *** Jump out of the layer loop
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  
  ! *** ZOOPLANKTON
  IF(IWQZPL > 0 )THEN
    CALL ZOOPL_CONTROL
  ENDIF
  
  ! *** SEDIMENT DIAGENESIS
  IF( IWQBEN == 1 )THEN
    ! *** Both Sediment Diagenesis initializations use R3D_Global in place of WQV
    ALLOCATE(R3D_Global(LCM_Global,KCM,NWQVM))
    R3D_Global = 0
    
    DO L=2,LA
      SMHYST(L)=.FALSE.
    ENDDO

    ! *** Read main sediment diagenesis control file  
    CALL SMRIN1_JNP
    !CALL SMRIN1
    
    IF( process_id == master_id )then
      ! *** READ SEDIMENT MODEL INITIAL CONDITION  
      IF( ISMICI == 1 )THEN
        CALL WQSDICI
      ENDIF
      
      IF( ISMICI == 2 )THEN
        CALL WQSDRST_IN
      ENDIF
    endif
    
    Call Broadcast_Array(SMPON_Global,  master_id)
    Call Broadcast_Array(SMPOP_Global,  master_id)
    Call Broadcast_Array(SMPOC_Global,  master_id)
  
    Call Broadcast_Array(SM1NH4_Global, master_id)
    Call Broadcast_Array(SM2NH4_Global, master_id)
    Call Broadcast_Array(SM2NO3_Global, master_id)
    Call Broadcast_Array(SM2PO4_Global, master_id)
    Call Broadcast_Array(SM2H2S_Global, master_id)
    Call Broadcast_Array(SMPSI_Global,  master_id)
    Call Broadcast_Array(SM2SI_Global,  master_id)
    Call Broadcast_Array(SMBST_Global,  master_id)
    Call Broadcast_Array(SMT_Global,    master_id)
    
    IF( ISMICI == 1 .OR. ISMICI == 2)THEN
      ! *** Map to Local Domain
      DO LG=2,LA_GLOBAL
        L = Map2Local(LG).LL
        IF( L > 1 )THEN
          SMPON(L,:) = SMPON_Global(LG,:)
          SMPOP(L,:) = SMPOP_Global(LG,:)
          SMPOC(L,:) = SMPOC_Global(LG,:)
        
          SM1NH4(L)  = SM1NH4_Global(LG)
          SM2NH4(L)  = SM2NH4_Global(LG)
          SM2NO3(L)  = SM2NO3_Global(LG)
          SM2PO4(L)  = SM2PO4_Global(LG)
          SM2H2S(L)  = SM2H2S_Global(LG)
          SMPSI(L)   = SMPSI_Global(LG)
          SM2SI(L)   = SM2SI_Global(LG)
          SMBST(L)   = SMBST_Global(LG)
          SMT(L)     = SMT_Global(LG)
        ENDIF
      ENDDO
    ENDIF
  
  ENDIF

  ! *** RPEM
  IF( ISRPEM > 0 )THEN
    CALL INIT_RPEMVARS
    CALL RPEMINP_JNP
  ENDIF
  
  ! *** READ WQ TIMESERIES
  IF( IWQPSL /= 0 )THEN        ! ** SCJ skip reading WQCSR if there are only constant point source loads
    if( process_id == master_id )then
      IF( ISWQLVL == 0 )then
        CALL WQCSR
      else 
        CALL WQCSR2
      endif
    endif

    Call Broadcast_Array(MCSER,   master_id)
    Call Broadcast_Array(TCCSER,  master_id)
    Call Broadcast_Array(TACSER,  master_id)

    DO NW=1,NWQV
      IF( NWQCSR(NW) >= 1 )THEN
        DO NS=1,NWQCSR(NW)
          Call Broadcast_Array(TSWQ(NS,NW).TIM,   master_id)
          Call Broadcast_Array(TSWQ(NS,NW).VAL,   master_id)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  CLOSE(2)

  RETURN
  
  END SUBROUTINE WQ3DINP

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQ3DCONTROL
  !
  !> @details  Read the Json formatted control file for water quality
  !---------------------------------------------------------------------------!
  !    Merged SNL and DSINTL codes
  !    READ IN FROM THE UNIT #8
  !    I/O CONTROL VARIABLES
  !    SPATIALLY AND TEMPORALLY CONSTANT REAL PARAMETERS
  !---------------------------------------------------------------------------!
  SUBROUTINE WQ3DCONTROL
  
  Use Allocate_Initialize      
  USE Variables_MPI
  Use Broadcast_Routines
  Use fson
  Use mod_fson_value, Only: fson_value_count, fson_value_get
  Use INFOMOD,ONLY:SKIPCOM,READSTR
  
  IMPLICIT NONE

  Type(fson_value), Pointer :: json_data, item, cell_group, pointsource, zones, phytogroups
  Character(len=79), allocatable :: TITLE(:)
  Character(len=1) :: CCMRM, SNUM
  Character*80 STR*200
  Integer :: I, IBIO, IZ, J, K, L, M, NW, NT, LL, N1, IM, KK, LG, NWQV0, NAL, ICOUNT
  Integer :: IWQDT, IJKC, IWQZX, ITMP,  II, JJ
  Real    :: XMRM1,  XMRM2,  XMRM3, XMRM4, XMRMA, XMRMB, XMRMC, XMRMD, XMRME  ! MACROALGAE
  Real    :: XPSQ, WQTT, XDSQ, O2WQ_, WQTAMD
  Real    :: XWQCHL, XWQTAMP, XWQPO4D, XWQSAD, TVARWQ, WTEMP, TT20
  Real    :: XMUD, IZMUD, IZSAND, WQKHRA_TEMP, WQDOP_TEMP, WQKDCALM_TEMP
  Real, parameter :: CONV1 = 1.0
  Real, parameter :: CONV2 = 8.64E4
  Real, allocatable :: WQCI(:), WQCPPRM(:), XPSL(:), XDSL(:)

  ! *** Nutrient components name
  WQCONSTIT(IROC) = 'ROC'
  WQCONSTIT(ILOC) = 'LOC'
  WQCONSTIT(IDOC) = 'DOC'
  WQCONSTIT(IROP) = 'ROP'
  WQCONSTIT(ILOP) = 'LOP'
  WQCONSTIT(IDOP) = 'DOP'
  WQCONSTIT(IP4D) = 'P4D'
  WQCONSTIT(IRON) = 'RON'
  WQCONSTIT(ILON) = 'LON'
  WQCONSTIT(IDON) = 'DON'
  WQCONSTIT(INHX) = 'NHX'
  WQCONSTIT(INOX) = 'NOX'
  WQCONSTIT(ISUU) = 'SUU'
  WQCONSTIT(ISAA) = 'SAA'
  WQCONSTIT(ICOD) = 'COD'
  WQCONSTIT(IDOX) = 'DOX'
  WQCONSTIT(ITAM) = 'TAM'
  WQCONSTIT(IFCB) = 'FCB'
  WQCONSTIT(ICO2) = 'CO2'
  
  ! *** Algae name if WQSKE1
  IF( ISWQLVL == 0 )THEN
    WQCONSTIT(ICHC) = 'CHC'
    WQCONSTIT(ICHD) = 'CHD'
    WQCONSTIT(ICHG) = 'CHG'
  ELSE
  ! *** Algae name if WQSKE2
    DO NW = 1,NALGAE
      WRITE(SNUM,'(I1)') NW
      WQCONSTIT(19+NW) = 'ALG'//SNUM
    ENDDO
  ENDIF
  
  ! *** Zooplankton name
  IF( IWQZPL > 0 )THEN
    DO NW = 1,NZOOPL
      WRITE(SNUM,'(I1)') NW
      WQCONSTIT(19+NALGAE+NW) = 'ZOO'//SNUM
    ENDDO
  ENDIF
  
  IF( .NOT. ALLOCATED(XDSL) )THEN
    ALLOCATE(XDSL(NWQVM))
    ALLOCATE(XPSL(NWQVM))
    XDSL = 0.0
    XPSL = 0.0
  ENDIF
  
  if( process_id == master_id )then
    json_data => fson_parse("wq_3dwc.jnp")
    Write(*,'(A)') 'WQ: READING WQ_3DWC.JNP - MAIN WATER QUALITY CONTROL FILE'
    Write(2,'(/,A)') 'READING WQ_3DWC.JNP - MAIN WATER QUALITY CONTROL FILE'
         
    Call fson_get(json_data, "title",                                           TITLE)
    Call fson_get(json_data, "kinetics_option",                                 ISWQLVL)
    Call fson_get(json_data, "number_of_variables",                             NWQV)
    Call fson_get(json_data, "use_kinetic_zones",                               IWQZONES)
    Call fson_get(json_data, "number_of_kinetic_zones",                         NWQZ)
    Call fson_get(json_data, "temperature_lookup_table_size",                   NWQTD)
    Call fson_get(json_data, "number_of_time_series_output_locations",          NWQTS)
    Call fson_get(json_data, "number_of_time_series_output_variables",          NTSWQV)
    Call fson_get(json_data, "number_of_sediment_zones",                        NSMZ)

    Call fson_get(json_data, "number_of_sediment_time_series_output_variables", NTSSMV)
    Call fson_get(json_data, "max_number_of_time_series_output_locations",      NSMTS)
    Call fson_get(json_data, "kinetic_update_time_step",                        WQKINUPT)
    
    If( ISWQLVL < 0 .OR. ISWQLVL > 4 )Call STOPP('BAD KINETICS OPTION')
    DTD = DT/86400.0
    WQKINUPT = WQKINUPT/86400.
    DTWQ = WQKINUPT
    DTWQO2 = DTWQ*0.5
    
    Call fson_get(json_data, "active_constituents.ROC",   ISKINETICS(IROC))
    Call fson_get(json_data, "active_constituents.LOC",   ISKINETICS(ILOC))
    Call fson_get(json_data, "active_constituents.DOC",   ISKINETICS(IDOC))
    Call fson_get(json_data, "active_constituents.ROP",   ISKINETICS(IROP))
    Call fson_get(json_data, "active_constituents.LOP",   ISKINETICS(ILOP))
    Call fson_get(json_data, "active_constituents.DOP",   ISKINETICS(IDOP))
    Call fson_get(json_data, "active_constituents.P4D",   ISKINETICS(IP4D))
    Call fson_get(json_data, "active_constituents.RON",   ISKINETICS(IRON))
    Call fson_get(json_data, "active_constituents.LON",   ISKINETICS(ILON))
    Call fson_get(json_data, "active_constituents.DON",   ISKINETICS(IDON))
    Call fson_get(json_data, "active_constituents.NHX",   ISKINETICS(INHX))
    Call fson_get(json_data, "active_constituents.NOX",   ISKINETICS(INOX))
    Call fson_get(json_data, "active_constituents.SUU",   ISKINETICS(ISUU))
    Call fson_get(json_data, "active_constituents.SAA",   ISKINETICS(ISAA))
    Call fson_get(json_data, "active_constituents.COD",   ISKINETICS(ICOD))
    Call fson_get(json_data, "active_constituents.DOX",   ISKINETICS(IDOX))
    Call fson_get(json_data, "active_constituents.TAM",   ISKINETICS(ITAM))
    Call fson_get(json_data, "active_constituents.FCB",   ISKINETICS(IFCB))
    Call fson_get(json_data, "active_constituents.CO2",   ISKINETICS(ICO2))
    If( ISWQLVL == 0  )then
      Call fson_get(json_data, "active_constituents.CHC", ISKINETICS(ICHC))
      Call fson_get(json_data, "active_constituents.CHD", ISKINETICS(ICHD))
      Call fson_get(json_data, "active_constituents.CHG", ISKINETICS(ICHG))
    End if

    If( ISWQLVL == 0  )then
      ISKINETICS = 0      
      ISKINETICS(IDOC) = 1
      ISKINETICS(INHX) = 1
      ISKINETICS(IDOX) = 1
    End if
       
    Call fson_get(json_data, "number_of_algae_groups", NALGAE)
    IF( ISWQLVL == 0 )then
      ALG_COUNT = 3
    Else
      ALG_COUNT = NALGAE
    End if
    
    ! *** Global Activation Options
    Call fson_get(json_data, "silica_activate",                  IWQSI)
    Call fson_get(json_data, "cyanobacteria_salinity_toxicity",  IWQSTOX)
    Call fson_get(json_data, "shellfish_farm_activate",          ISFFARM)
    Call fson_get(json_data, "number_of_shellfish_species",      NSF)
    Call fson_get(json_data, "number_of_shellfish_cells",        NSFCELLS)
    Call fson_get(json_data, "zooplankton_activate",             IWQZPL)
    Call fson_get(json_data, "number_of_zooplankton_groups",     NZOOPL)
    Call fson_get(json_data, "rpem_activate",                    ISRPEM)
    
    If( IWQZPL > 0 )then
      NWQVZ = NWQV - NZOOPL
    End if

    Call fson_get(json_data, "po4_sorption_option",                                 IWQSRP)
    Call fson_get(json_data, "log_negative_concentrations",                         IWQNC)
    Call fson_get(json_data, "write_restart",                                       IWQRST)
    
    Call fson_get(json_data, "formulation_for_DO_saturation",                       IDOSFRM)
    Call fson_get(json_data, "elevation_adjustment_for_DO_saturation",              IDOSELE)
    Call fson_get(json_data, "elevation_offset_for_DO_saturation",                  DOELEV)
    
    Call fson_get(json_data, "number_of_hours_averaging_DO",                        WQHRAVG)
                                                                                    
    Call fson_get(json_data, "initial_condition_option",                            IWQICI)
    Call fson_get(json_data, "point_source_load_option",                            IWQPSL)
    Call fson_get(json_data, "use_atmospheric_dry_deposition",                      IWQNPL)
    Call fson_get(json_data, "solar_radiation.source_option",                       IWQSUN)
    Call fson_get(json_data, "solar_radiation.initial_optimal_sr",                  WQI0)
    Call fson_get(json_data, "solar_radiation.minimum_optimal_sr",                  WQISMIN)
    Call fson_get(json_data, "solar_radiation.fraction_of_daylight",                WQFD)
    Call fson_get(json_data, "solar_radiation.photoactive_radiation_fraction",      PARADJ)
    Call fson_get(json_data, "solar_radiation.daily_weighting_factors",             WQCI)
    WQCIA = WQCI(1)
    WQCIB = WQCI(2)
    WQCIC = WQCI(3)
    WQCIM = WQCI(4)
    
    WQI0 = PARADJ*WQI0        ! *** Apply conversion to photosynthesis active light fraction
    WQI1 = WQI0
    WQI2 = WQI0
    WQI3 = WQI0
    WQI0OPT = 0.0
    
    Call fson_get(json_data, "light_extinction.light_extinction_diagnostics",      NDLTAVG)
    Call fson_get(json_data, "light_extinction.chlorophyll_coefficient",           WQKECHL)
    Call fson_get(json_data, "light_extinction.chlorophyll_exponent",              WQKECHLE)
    Call fson_get(json_data, "light_extinction.particular_organic_matter_coeff",   WQKEPOC)
    Call fson_get(json_data, "light_extinction.dissolved_organic_carbon_coeff",    WQKEDOM)
    Call fson_get(json_data, "light_extinction.background_coeff",                  WQKEB(1))
    Call fson_get(json_data, "light_extinction.total_suspended_solids_coeff",      WQKETSS)
    If( ISTRAN(6) == 0 .AND. ISTRAN(7) == 0  )then
      WQKETSS = 0.0
    End if
    
    Call fson_get(json_data, "reaeration.reaeration_option",                               IWQKA(1))
    Call fson_get(json_data, "reaeration.reaeration_constant",                             WQKRO(1))
    Call fson_get(json_data, "reaeration.temperature_rate_const",                          WQKTR(1))
    Call fson_get(json_data, "reaeration.adjustment_factor",                               REAC(1))
    Call fson_get(json_data, "nutrient_sorption.partition_coeff_for_sorbed_dissolved_PO4", WQKPO4P)
    Call fson_get(json_data, "nutrient_sorption.partition_coeff_for_sorbed_dissolved_SA",  WQKSAP)
    If( IWQSRP /= 1 .AND. IWQSRP /= 2  )then
      WQKPO4P = 0.0
      WQKSAP = 0.0
    End if
          
    Call fson_get(json_data, "hydrolysis.reference_temperature",                       WQTRHDR)
    Call fson_get(json_data, "hydrolysis.effect_of_temperature",                       WQKTHDR)
    Call fson_get(json_data, "hydrolysis.carbon.minimum_rate.RPOC",                    WQKRC)
    Call fson_get(json_data, "hydrolysis.carbon.minimum_rate.LPOC",                    WQKLC)
    Call fson_get(json_data, "hydrolysis.carbon.constant_relating_to_algae.RPOC",      WQKRCALG)
    Call fson_get(json_data, "hydrolysis.carbon.constant_relating_to_algae.LPOC",      WQKLCALG)
    Call fson_get(json_data, "hydrolysis.phosphorus.minimum_rate.RPOP",                WQKRP)
    Call fson_get(json_data, "hydrolysis.phosphorus.minimum_rate.LPOP",                WQKLP)
    Call fson_get(json_data, "hydrolysis.phosphorus.constant_relating_to_algae.RPOP",  WQKRPALG)
    Call fson_get(json_data, "hydrolysis.phosphorus.constant_relating_to_algae.LPOP",  WQKLPALG)
    Call fson_get(json_data, "hydrolysis.phosphorus.carbon_to_phosphorus_ratio",       WQCPPRM)
    Call fson_get(json_data, "hydrolysis.nitrogen.minimum_rate.RPON",                  WQKRN)
    Call fson_get(json_data, "hydrolysis.nitrogen.minimum_rate.LPON",                  WQKLN)
    Call fson_get(json_data, "hydrolysis.nitrogen.constant_relating_to_algae.RPON",    WQKRNALG)
    Call fson_get(json_data, "hydrolysis.nitrogen.constant_relating_to_algae.LPON",    WQKLNALG)

    WQCP1PRM = WQCPPRM(1)
    WQCP2PRM = WQCPPRM(2)
    WQCP3PRM = WQCPPRM(3)

    Call fson_get(json_data, "mineralization.reference_temperature",                      WQTRMNL)
    Call fson_get(json_data, "mineralization.effect_of_temperature",                      WQKTMNL)
    Call fson_get(json_data, "mineralization.carbon.minimum_rate.DOC",                    WQKDC(1))
    Call fson_get(json_data, "mineralization.carbon.constant_relating_to_algae.DOC",      WQKDCALG)
    Call fson_get(json_data, "mineralization.carbon.constant_relating_to_macroalgae.DOC", WQKDCALM_TEMP)
    
    Call fson_get(json_data, "mineralization.phosphorus.minimum_rate.DOP",                WQKDP)
    Call fson_get(json_data, "mineralization.phosphorus.constant_relating_to_algae.DOP",  WQKDPALG)

    Call fson_get(json_data, "mineralization.nitrogen.minimum_rate.DON",                  WQKDN)
    Call fson_get(json_data, "mineralization.nitrogen.constant_relating_to_algae.DON",    WQKDNALG)
    
                                                                                       
    Call fson_get(json_data, "nitrification.mass_NO3_reduces_per_DOC_oxidized",        WQANDC)
    Call fson_get(json_data, "nitrification.max_rate",                                 WQNITM)
    Call fson_get(json_data, "nitrification.half_sat_const_for_DO",                    WQKHNDO)
    Call fson_get(json_data, "nitrification.half_sat_const_for_NH4",                   WQKHNN)
    Call fson_get(json_data, "nitrification.reference_temperature",                    WQTNIT)
    Call fson_get(json_data, "nitrification.suboptimal_temperature_effect_const",      WQKN1)
    Call fson_get(json_data, "nitrification.superoptimal_temperature_effect_const",    WQKN2)
    
    Call fson_get(json_data, "denitrification.oxic_respiration_half_sat_const_for_DO", WQKHORDO)
    Call fson_get(json_data, "denitrification.half_sat_const",                         WQKHDNN)
    Call fson_get(json_data, "denitrification.ratio_to_oxic_DOC_respiration",          WQAANOX)
    WQAANOX = WQAANOX*WQKHORDO
    
    Call fson_get(json_data, "silica_dissolution.dissolution_rate",           WQKSU)
    Call fson_get(json_data, "silica_dissolution.reference_temperature",      WQTRSUA)
    Call fson_get(json_data, "silica_dissolution.effect_of_temperature",      WQKTSUA)
                                                                              
    Call fson_get(json_data, "TAM_release.half_anoxic_rate_DO",               WQKHBMF)
    Call fson_get(json_data, "TAM_release.anoxic_release_rate",               WQBFTAM)
    Call fson_get(json_data, "TAM_release.reference_temperature",             WQTTAM)
    Call fson_get(json_data, "TAM_release.effect_of_temperature",             WQKTAM)
    Call fson_get(json_data, "TAM_release.solubility_at_anoxic_conditions",   WQTAMDMX)
    Call fson_get(json_data, "TAM_release.solubility_to_DO_const",            WQKDOTAM)
                                                                              
    Call fson_get(json_data, "coliform_decay_rate.first_order_decay_rate",    WQKFCB)
    Call fson_get(json_data, "coliform_decay_rate.temperature_effect_const",  WQTFCB)
    
    Call fson_get(json_data, "COD_decay.oxygen_half_sat_const_for_COD_decay", WQKHCOD(1))
    Call fson_get(json_data, "COD_decay.COD_decay_rate",                      WQKCD(1))
    Call fson_get(json_data, "COD_decay.reference_temperature",               WQTRCOD)
    Call fson_get(json_data, "COD_decay.effect_of_temperature",               WQKTCOD)
    
    Call fson_get(json_data, "settling_velocity.refractory_POM",              WQWSRP(1))
    Call fson_get(json_data, "settling_velocity.labile_POM",                  WQWSLP(1))
    Call fson_get(json_data, "settling_velocity.particles_sorbed_to_TAM",     WQWSS(1))

    !---------------------------------------------------------------------------
    ! C47: Spatially/Temporally constant benthic fluxes
    !---------------------------------------------------------------------------
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_option",                   IWQBEN)
    Call fson_get(json_data, "sediment_diagenesis.number_of_reactive_classes",            NSMG)
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.phosphate",          WQBFPO4D(1))
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.ammonia",            WQBFNH4(1))
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.nitrate",            WQBFNO3(1))
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.silica",             WQBFSAD(1))
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.COD",                WQBFCOD(1))
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.SOD",                WQBFO2(1))
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.temperature_factor", STEMFAC)     

    ! *** Initialize all zones
    Do IZ = 2, NWQZ
      WQKEB(IZ)   = WQKEB(1)
      IWQKA(IZ)   = IWQKA(1)
      WQKRO(IZ)   = WQKRO(1)
      WQKTR(IZ)   = WQKTR(1)
      REAC(IZ)    = REAC(1)
      WQKDC(IZ)   = WQKDC(1)
      WQKHCOD(IZ) = WQKHCOD(1)
      WQKCD(IZ)   = WQKCD(1)
      WQWSRP(IZ)  = WQWSRP(1)
      WQWSLP(IZ)  = WQWSLP(1)
      WQWSS(IZ)   = WQWSS(1)
    End do

    !---------------------------------------------------------------------------
    ! C30: Concentration time series data for open boundaries
    ! Number of time series for each state variables of nutrient
    !---------------------------------------------------------------------------
    DO NW = 1,NWQV
      Call fson_get(json_data,'number_of_time_series.'//WQCONSTIT(NW), NWQCSR(NW))
    ENDDO
    
    NT = 0
    DO NW = 1,NWQV
      NT = MAX(NT,NWQCSR(NW))
    ENDDO
    NCSER(8) = NT
    
    !---------------------------------------------------------------------------
    ! C31 + C32 + C33 + C34: South Boundary
    !---------------------------------------------------------------------------
    Call fson_get(json_data, "open_boundaries.south.number_of_cells", NWQOBS)
    Call fson_get(json_data, "open_boundaries.south.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      Call fson_get(item, "I", IWQCBS(M))
      Call fson_get(item, "J", JWQCBS(M))
      Do NW = 1, NWQV
        Call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBS(M,NW))
        Call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCS(M,1,NW))
        Call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCS(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      IF( IWQCBS(M) == ICBS(M) .AND. JWQCBS(M) == JCBS(M) )THEN
        NCSERS(M,8) = IWQOBS(M,IDOX)   ! *** ALL CONSTITUENTS USE THE SAME SERIES
        Do NW = 1, NWQV
          IF( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom
            CBS(M,1,NT) = WQOBCS(M,1,NW)
            ! *** Top
            CBS(M,2,NT) = WQOBCS(M,2,NW)
          End if
        End do
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  CALL STOPP('WQ: SOUTH OBC: MISS MATCH BETWEEN NCBS & NWQOBS')
      ENDIF
    Enddo
    
    !---------------------------------------------------------------------------  
    ! C31 + C35 + C36 + C37: West Boundary  
    !---------------------------------------------------------------------------
    Call fson_get(json_data, "open_boundaries.west.number_of_cells", NWQOBW) 
    Call fson_get(json_data, "open_boundaries.west.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      Call fson_get(item, "I", IWQCBW(M))
      Call fson_get(item, "J", JWQCBW(M))
      Do NW = 1, NWQV
        Call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBW(M,NW))
        Call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCW(M,1,NW))
        Call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCW(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      IF( IWQCBW(M) == ICBW(M) .AND. JWQCBW(M) == JCBW(M) )THEN
        NCSERW(M,8) = IWQOBW(M,IDOX)   ! *** ALL CONSTITUENTS USE THE SAME SERIES
        Do NW = 1, NWQV
          IF( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom
            CBW(M,1,NT) = WQOBCW(M,1,NW)
            ! *** Top
            CBW(M,2,NT) = WQOBCW(M,2,NW)
          End if
        End do
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  CALL STOPP('WQ: WST OBC: MISS MATCH BETWEEN NCBW & NWQOBW')
      ENDIF
    End do
    !---------------------------------------------------------------------------  
    ! C30 + C38 + C39 + C40: East Boundary  
    !---------------------------------------------------------------------------
    Call fson_get(json_data, "open_boundaries.east.number_of_cells", NWQOBE) 
    Call fson_get(json_data, "open_boundaries.east.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      Call fson_get(item, "I", IWQCBE(M))
      Call fson_get(item, "J", JWQCBE(M))
      Do NW = 1, NWQV
        Call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBE(M,NW))
        Call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCE(M,1,NW))
        Call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCE(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      IF( IWQCBE(M) == ICBE(M) .AND. JWQCBE(M) == JCBE(M) )THEN
        NCSERE(M,8) = IWQOBE(M,IDOX)   ! *** ALL CONSTITUENTS USE THE SAME SERIES
        Do NW = 1, NWQV
          IF( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom  
            CBE(M,1,NT) = WQOBCE(M,1,NW)
            ! *** Top
             CBE(M,2,NT) = WQOBCE(M,2,NW)
          End if
        End do
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  CALL STOPP('WQ: EAST OBC: MISS MATCH BETWEEN NCBE & NWQOBE')
      ENDIF
    Enddo
    !---------------------------------------------------------------------------  
    ! C30 + C41 + C42 + C43: North Boundary  
    !---------------------------------------------------------------------------
    Call fson_get(json_data, "open_boundaries.north.number_of_cells", NWQOBN)    
    Call fson_get(json_data, "open_boundaries.north.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      Call fson_get(item, "I", IWQCBN(M))
      Call fson_get(item, "J", JWQCBN(M))
      Do NW = 1, NWQV
        Call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBN(M,NW))
        Call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCN(M,1,NW))
        Call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCN(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      IF( IWQCBN(M) == ICBN(M) .AND. JWQCBN(M) == JCBN(M) )THEN
        NCSERN(M,8) = IWQOBN(M,IDOX)   ! *** ALL CONSTITUENTS USE THE SAME SERIES
        Do NW = 1, NWQV
          IF( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom  
            CBN(M,1,NT) = WQOBCN(M,1,NW)
            ! *** Top
            CBN(M,2,NT) = WQOBCN(M,2,NW)
          End if
        End do
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  CALL STOPP('WQ: NORTH OBC: MISS MATCH BETWEEN NCBN & NWQOBN')
      ENDIF
    Enddo
    !---------------------------------------------------------------------------
    ! C48: TEMPORALLY-CONSTANT VALUES FOR POINT SOURCE CONCENTRATIONS IN MG/L
    !---------------------------------------------------------------------------
    Call fson_get(json_data, "mass_loading_point_sources.number_of_point_sources",  NWQPS) 
    Call fson_get(json_data, "mass_loading_point_sources.number_of_time_series",    NPSTMSR) 
    Call fson_get(json_data, "mass_loading_point_sources.constant_point_sources",   pointsource)
    ICOUNT = fson_value_count(pointsource)
    DO N1 = 1, ICOUNT   !  fson_value_count(pointsource)
      item => fson_value_get(pointsource, N1)
      Call fson_get(item, "I", I)
      Call fson_get(item, "J", J)
      Call fson_get(item, "K", K)
      Call fson_get(item, "NSR", ITMP)
      Call fson_get(item, "PSQ", XPSQ)
      Do NW = 1, NWQV
        Call fson_get(item,'concentrations.'//WQCONSTIT(NW), XPSL(NW))
      Enddo
      
      IF( IJCT_Global(I,J) < 1 .OR. IJCT_Global(I,J) > 8 )THEN
        CALL STOPP('ERROR!! INVALID (I,J) IN FILE WQ_3DWC.INP FOR PSL')
      ENDIF 
    
      ! *** HANDLE CONCENTRATION BASED POINT SOURCE
      IF( IWQPSL == 2 )THEN
        IF( LIJ_Global(I,J) /= LQS_GL(N1) )THEN
          CALL STOPP('MISMATCH NQSIJ BETWEEN EFDC.INP & WQ_3DWC.INP')
        ENDIF

        ! *** ASSIGN GLOBAL CONCENTRATION TIME SERIES INDEX
        NSERWQ(N1) = ITMP  ! *** ALL WQ VARIABLES USE SAME TIME SERIES (NSERWQ IS TEMPORARY STORAGE FOR LATER MAPPING TO LOCAL NCSERQ(N1,8))
        DO NW = 1,NWQV
          IF( ISKINETICS(NW) > 0 )THEN
            WQWPSLC(N1,NW) = XPSL(NW)
          ENDIF
        ENDDO
      ELSE            
        ICPSL(N1) = I
        JCPSL(N1) = J
        KCPSL(N1) = K
        MVPSL(N1) = ITMP
        
        ! *** SAVE FOR NOW, CONVERT IN SUB-DOMAIN MAPPING
        WQTT = XPSQ*CONV2                   ! *** CONVERT FROM M3/S TO M3/DAY
        DO NW=1,NWQV
          WQWPSLC(N1,NW) = XPSL(NW)*WQTT
        ENDDO  
      ENDIF
    ENDDO
    
    !---------------------------------------------------------------------------
    ! C44: Spatially/Temporally Constant Initial Concentration
    !---------------------------------------------------------------------------
    Do NW = 1,NWQV
      Call fson_get(json_data,'const_initial_conditions.'//WQCONSTIT(NW), WQV(1,1,NW))
    End do
    !---------------------------------------------------------------------------
    ! C49: Constant Dry Atmospheric Deposition(g / m2 / day; MPN / m2 / day)
    !---------------------------------------------------------------------------
    Call fson_get(json_data,"dry_atmospheric_deposition.DSQ", XDSQ)
    Do NW = 1,NWQV
      Call fson_get(json_data,'dry_atmospheric_deposition.'//WQCONSTIT(NW), XDSL(NW))
    End do
    
    !---------------------------------------------------------------------------
    ! C50: Wet Atmospheric Deposition concentrations (mg / L, TAM - MOLES / L, FCB - MPN / 100ml)
    !---------------------------------------------------------------------------
    Do NW = 1,NWQV
      Call fson_get(json_data,'wet_atmospheric_deposition.'//WQCONSTIT(NW), WQATM(NW,1))
    Enddo        
 
    ! *** DEFAULT ZONAL KINETICS,   IWQZONES = 0 OR 1
    DO I = 2,NWQZ
      REAC(I)    = REAC(1)
      IWQKA(I)   = IWQKA(1)
      WQKRO(I)   = WQKRO(1)
      WQKTR(I)   = WQKTR(1)
      WQKDC(I)   = WQKDC(1)
      WQKHCOD(I) = WQKHCOD(1)
      WQKCD(I)   = WQKCD(1)
      WQWSRP(I)  = WQWSRP(1)
      WQWSLP(I)  = WQWSLP(1)
      WQWSS(I)   = WQWSS(1)
      WQKEB(I)   = WQKEB(1)
    ENDDO
    
    DO I = 1,NWQZ
      DO NAL = 1,NALGAE
        ALGAES(NAL).WQKDCALM(I) = WQKDCALM_TEMP
      ENDDO
    ENDDO
           
    ! *** Reading the kinetics file for spatial zones WQ parameters 
    IF( IWQZONES  > 0 )THEN
      WRITE(*,'(A)')' WQ: WQ_KIN_ZONES.JNP'
      WRITE(2,'(/,A)')'READING WQ_KIN_ZONES.JNP - Zonal kinetics'
      json_data => fson_parse("wq_kin_zones.jnp")
       
      Call fson_get(json_data, "reaeration.reaeration_option",                  IWQKA)    ! EE WQZoneP[0]
      Call fson_get(json_data, "reaeration.reaeration_constant",                WQKRO)    ! EE WQZoneP[1]
      Call fson_get(json_data, "reaeration.temperature_rate_const",             WQKTR)    ! EE WQZoneP[2]
      Call fson_get(json_data, "reaeration.adjustment_factor",                  REAC)     ! EE WQZoneP[3]
      Call fson_get(json_data, "minimum_DOC_hydrolysis_rate",                   WQKDC)    ! EE WQZoneP[4]
      Call fson_get(json_data, "COD_decay.COD_decay_rate",                      WQKCD)    ! EE WQZoneP[5]
      Call fson_get(json_data, "COD_decay.oxygen_half_sat_const_for_COD_decay", WQKHCOD)  ! EE WQZoneP[6]
      CALL fson_get(json_data, "background_light_extinction_coeff",             WQKEB)    ! EE WQZoneP[7]
      CALL fson_get(json_data, "settling_velocity.RPOM",                        WQWSRP)   ! EE WQZoneP[8]
      CALL fson_get(json_data, "settling_velocity.LPOM",                        WQWSLP)   ! EE WQZoneP[9]
      CALL fson_get(json_data, "settling_velocity.TAM",                         WQWSS)    ! EE WQZoneP[10]
    ENDIF
    
  End if  ! *** End of master_id block
  
  ! *** Broadcast water quality global settings, not dependent on grid. 

  ! *** Broadcast scalars
  Call Broadcast_Scalar(ISWQLVL,     master_id)
  Call Broadcast_Scalar(IWQZONES,    master_id)        
  Call Broadcast_Scalar(NWQZ,        master_id)
  Call Broadcast_Scalar(NWQTD,       master_id)
  Call Broadcast_Scalar(NWQTS,       master_id)
  Call Broadcast_Scalar(NTSWQV,      master_id)
  Call Broadcast_Scalar(NSMZ,        master_id)
                                     
  Call Broadcast_Scalar(NTSSMV,      master_id)
  Call Broadcast_Scalar(NSMTS,       master_id)
  Call Broadcast_Scalar(WQKINUPT,    master_id)                                   
  Call Broadcast_Scalar(DTD,         master_id)

  Call Broadcast_Scalar(DTWQ,        master_id)
  Call Broadcast_Scalar(DTWQO2,      master_id)
                                     
  Call Broadcast_Scalar(NALGAE,      master_id)
  Call Broadcast_Scalar(ALG_COUNT,   master_id)
                                     
  Call Broadcast_Scalar(IWQSI,       master_id)
  Call Broadcast_Scalar(IWQSTOX,     master_id)
  Call Broadcast_Scalar(ISFFARM,     master_id)
  Call Broadcast_Scalar(NSF,         master_id)
  Call Broadcast_Scalar(NSFCELLS,    master_id)
  Call Broadcast_Scalar(IWQZPL,      master_id)
  Call Broadcast_Scalar(NZOOPL,      master_id)
  Call Broadcast_Scalar(ISRPEM,      master_id)
  Call Broadcast_Scalar(NWQVZ,       master_id)
                                     
  Call Broadcast_Scalar(IWQSRP,      master_id)
  Call Broadcast_Scalar(IWQNC,       master_id)
  Call Broadcast_Scalar(IWQRST,      master_id)
    
  Call Broadcast_Scalar(IDOSFRM,     master_id)
  Call Broadcast_Scalar(IDOSELE,     master_id)
  Call Broadcast_Scalar(DOELEV,      master_id)
                                     
  Call Broadcast_Scalar(WQHRAVG,     master_id)
                                     
  Call Broadcast_Scalar(IWQICI,      master_id)
  Call Broadcast_Scalar(IWQPSL,      master_id)
  Call Broadcast_Scalar(IWQNPL,      master_id)
  Call Broadcast_Scalar(IWQSUN,      master_id)
  Call Broadcast_Scalar(WQI0,        master_id)
  Call Broadcast_Scalar(WQISMIN,     master_id)
  Call Broadcast_Scalar(WQFD,        master_id)
  Call Broadcast_Scalar(PARADJ,      master_id)
                                     
  Call Broadcast_Scalar(WQCIA,       master_id)
  Call Broadcast_Scalar(WQCIB,       master_id)
  Call Broadcast_Scalar(WQCIC,       master_id)
  Call Broadcast_Scalar(WQCIM,       master_id)
  Call Broadcast_Scalar(WQI1,        master_id)
  Call Broadcast_Scalar(WQI2,        master_id)
  Call Broadcast_Scalar(WQI3,        master_id)
  Call Broadcast_Scalar(WQI0OPT,     master_id)
                                     
  Call Broadcast_Scalar(NDLTAVG,     master_id)
  Call Broadcast_Scalar(WQKECHL,     master_id)
  Call Broadcast_Scalar(WQKECHLE,    master_id)
  Call Broadcast_Scalar(WQKEPOC,     master_id)
  Call Broadcast_Scalar(WQKEDOM,     master_id)
  Call Broadcast_Scalar(WQKEB(1),    master_id)
  Call Broadcast_Scalar(WQKETSS,     master_id)
                                     
  Call Broadcast_Scalar(IWQKA(1),    master_id)
  Call Broadcast_Scalar(WQKRO(1),    master_id)
  Call Broadcast_Scalar(WQKTR(1),    master_id)
  Call Broadcast_Scalar(REAC(1),     master_id)
  Call Broadcast_Scalar(WQKPO4P,     master_id)
  Call Broadcast_Scalar(WQKSAP,      master_id)
                                     
  Call Broadcast_Scalar(WQTRHDR,     master_id)
  Call Broadcast_Scalar(WQKTHDR,     master_id)
  Call Broadcast_Scalar(WQKRC,       master_id)
  Call Broadcast_Scalar(WQKLC,       master_id)
  Call Broadcast_Scalar(WQKDC(1),    master_id)
  Call Broadcast_Scalar(WQKRCALG,    master_id)
  Call Broadcast_Scalar(WQKLCALG,    master_id)
  Call Broadcast_Scalar(WQKDCALG,    master_id)
  Call Broadcast_Scalar(WQKRP,       master_id)
  Call Broadcast_Scalar(WQKLP,       master_id)
  Call Broadcast_Scalar(WQKDP,       master_id)
  Call Broadcast_Scalar(WQKRPALG,    master_id)
  Call Broadcast_Scalar(WQKLPALG,    master_id)
  Call Broadcast_Scalar(WQKDPALG,    master_id)
                                     
  Call Broadcast_Scalar(WQCP1PRM,    master_id)
  Call Broadcast_Scalar(WQCP2PRM,    master_id)
  Call Broadcast_Scalar(WQCP3PRM,    master_id)
                                     
  Call Broadcast_Scalar(WQKRN,       master_id)
  Call Broadcast_Scalar(WQKLN,       master_id)
  Call Broadcast_Scalar(WQKDN,       master_id)
  Call Broadcast_Scalar(WQKRNALG,    master_id)
  Call Broadcast_Scalar(WQKLNALG,    master_id)
  Call Broadcast_Scalar(WQKDNALG,    master_id)
                                     
  Call Broadcast_Scalar(WQTRMNL,     master_id)
  Call Broadcast_Scalar(WQKTMNL,     master_id)
                                     
  Call Broadcast_Scalar(WQANDC,      master_id)
  Call Broadcast_Scalar(WQNITM,      master_id)
  Call Broadcast_Scalar(WQKHNDO,     master_id)
  Call Broadcast_Scalar(WQKHNN,      master_id)
  Call Broadcast_Scalar(WQTNIT,      master_id)
  Call Broadcast_Scalar(WQKN1,       master_id)
  Call Broadcast_Scalar(WQKN2,       master_id)
                                     
  Call Broadcast_Scalar(WQKHORDO,    master_id)
  Call Broadcast_Scalar(WQKHDNN,     master_id)
  Call Broadcast_Scalar(WQAANOX,     master_id)
                                     
  Call Broadcast_Scalar(WQKSU,       master_id)
  Call Broadcast_Scalar(WQTRSUA,     master_id)
  Call Broadcast_Scalar(WQKTSUA,     master_id)
                                     
  Call Broadcast_Scalar(WQKHBMF,     master_id)
  Call Broadcast_Scalar(WQBFTAM,     master_id)
  Call Broadcast_Scalar(WQTTAM,      master_id)
  Call Broadcast_Scalar(WQKTAM,      master_id)
  Call Broadcast_Scalar(WQTAMDMX,    master_id)
  Call Broadcast_Scalar(WQKDOTAM,    master_id) 
                                     
  Call Broadcast_Scalar(WQKFCB,      master_id)
  Call Broadcast_Scalar(WQTFCB,      master_id)
  
  Call Broadcast_Scalar(WQKHCOD(1),  master_id)
  Call Broadcast_Scalar(WQKCD(1),    master_id)
  Call Broadcast_Scalar(WQTRCOD,     master_id)
  Call Broadcast_Scalar(WQKTCOD,     master_id)
                                     
  Call Broadcast_Scalar(WQWSRP(1),   master_id)
  Call Broadcast_Scalar(WQWSLP(1),   master_id)
  Call Broadcast_Scalar(WQWSS(1),    master_id)
                                     
  Call Broadcast_Scalar(IWQBEN,      master_id)
  Call Broadcast_Scalar(NSMG,        master_id) 
  Call Broadcast_Scalar(WQBFPO4D(1), master_id)
  Call Broadcast_Scalar(WQBFNH4(1),  master_id)
  Call Broadcast_Scalar(WQBFNO3(1),  master_id)
  Call Broadcast_Scalar(WQBFSAD(1),  master_id)
  Call Broadcast_Scalar(WQBFCOD(1),  master_id)
  Call Broadcast_Scalar(WQBFO2(1),   master_id)
  Call Broadcast_Scalar(STEMFAC,     master_id)
  
  Call Broadcast_Scalar(NCSER(8),    master_id)
                                     
  ! *** Arrays                       
  Call broadcast_array(ISKINETICS,   master_id)
  Call Broadcast_Array(NWQCSR,       master_id)

  ! *** Zonal Kinetics

  Call Broadcast_Array(IWQKA,      master_id)
  Call Broadcast_Array(WQKRO,      master_id)
  Call Broadcast_Array(WQKTR,      master_id)
  Call Broadcast_Array(REAC,       master_id)
  Call Broadcast_Array(WQKDC,      master_id)
  Call Broadcast_Array(WQKCD,      master_id)
  Call Broadcast_Array(WQKHCOD,    master_id)
  Call Broadcast_Array(WQKEB,      master_id)
  Call Broadcast_Array(WQWSRP,     master_id)
  Call Broadcast_Array(WQWSLP,     master_id)
  Call Broadcast_Array(WQWSS,      master_id)
  
  ! *** Dry Deposition
  Call Broadcast_Scalar(XDSQ,      master_id)
  Call Broadcast_Array(XDSL,       master_id)

  ! *** Wet Deposition
  Call Broadcast_Array(WQATM,      master_id)
  ! *** End of global water quality settings   
  
  ! **************************************************************************
  ! *** Open Boundaries
  Call Broadcast_Scalar(NWQOBS,    master_id)
  Call Broadcast_Scalar(NWQOBW,    master_id)
  Call Broadcast_Scalar(NWQOBE,    master_id)
  Call Broadcast_Scalar(NWQOBN,    master_id)

  Call Broadcast_Array(IWQCBS,     master_id)
  Call Broadcast_Array(JWQCBS,     master_id)
  Call Broadcast_Array(IWQOBS,     master_id)
  Call Broadcast_Array(WQOBCS,     master_id)
  Call Broadcast_Array(CBS,        master_id)

  Call Broadcast_Array(IWQCBW,     master_id)
  Call Broadcast_Array(JWQCBW,     master_id)
  Call Broadcast_Array(IWQOBW,     master_id)
  Call Broadcast_Array(WQOBCW,     master_id)
  Call Broadcast_Array(CBW,        master_id)

  Call Broadcast_Array(IWQCBE,     master_id)
  Call Broadcast_Array(JWQCBE,     master_id)
  Call Broadcast_Array(IWQOBE,     master_id)
  Call Broadcast_Array(WQOBCE,     master_id)
  Call Broadcast_Array(CBE,        master_id)

  Call Broadcast_Array(IWQCBN,     master_id)
  Call Broadcast_Array(JWQCBN,     master_id)
  Call Broadcast_Array(IWQOBN,     master_id)
  Call Broadcast_Array(WQOBCN,     master_id)
  Call Broadcast_Array(CBN,        master_id)
  
  ! *** Map to Local
  Call Map_OpenBC_Eutrophication

  ! **************************************************************************
  ! *** Point Sources
  Call Broadcast_Scalar(NWQPS,   master_id)
  Call Broadcast_Scalar(NPSTMSR, master_id)

  Call Broadcast_Array(ICPSL,    master_id)
  Call Broadcast_Array(JCPSL,    master_id)
  Call Broadcast_Array(KCPSL,    master_id)
  Call Broadcast_Array(MVPSL,    master_id)
  Call Broadcast_Array(NSERWQ,   master_id)
  Call Broadcast_Array(WQWPSLC,  master_id)
  Call Broadcast_Array(IWQPSC,   master_id)
  Call Broadcast_Array(IWQPSV,   master_id)  
  Call Broadcast_Array(CQS,      master_id)
  
  ! *** Map to Local
  Call Map_WQ_PointSource
  
  ! **************************************************************************
  ! *** Set up look-up table for temperature dependency over -10 °C to 50 °C
  WQTDMIN =-10
  WQTDMAX = 50
  WTEMP = WQTDMIN
  WQTDINC = (WQTDMAX-WQTDMIN)/NWQTD
  ALLOCATE (WQTDTEMP(NWQTD))
  
  DO M = 1,NWQTD
    WQTDTEMP(M) = WTEMP

    WQTDHDR(M) = EXP( WQKTHDR*(WTEMP-WQTRHDR) )   ! *** TEMPERATURE EFFECT ON HYDROLYSIS
    WQTDMNL(M) = EXP( WQKTMNL*(WTEMP-WQTRMNL) )   ! *** TEMPERATURE EFFECT ON MINERALIZATION
    
    ! *** TEMPERATURE ADJUSTED RATES FOR NITRIFICATION
    IF( WTEMP > WQTNIT )THEN
      WQTDNIT(M) = WQNITM*EXP(-WQKN2*(WQTNIT-WTEMP)**2)
    ELSE
      WQTDNIT(M) = WQNITM*EXP(-WQKN1*(WTEMP-WQTNIT)**2)
    ENDIF
    
    WQKSUA(M)    = WQKSU * EXP( WQKTSUA*(WTEMP-WQTRSUA) )        ! *** TEMPERATURE EFFECT ON PSI DISSOLUTION
    WQKCOD(M,1)  = WQKCD(1) * EXP( WQKTCOD*(WTEMP-WQTRCOD) )     ! *** TEMPERATURE EFFECT ON COD OXIDATION
    
    TT20 = WTEMP - 20.0
    WQTDKR(M,1) = WQKTR(1)**TT20                                 ! *** TEMPERATURE EFFECT ON REAERATION
    
    WQTDTAM(M) = WQKHBMF * WQBFTAM * EXP( WQKTAM*(WTEMP-WQTTAM) )
    
    WQTT = WQKFCB * WQTFCB**TT20 * DTWQO2
    WQTD1FCB(M) = 1.0 - WQTT
    WQTD2FCB(M) = 1.0 / (1.0 + WQTT)
    
    WTEMP = WTEMP + WQTDINC
  ENDDO
  
  if( process_id == master_id )then
    WRITE(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for Water Column Temperature Dependency",                               &
                                             "WQTDMIN = ", WQTDMIN, "WQTDMAX = ", WQTDMAX, "WQTDINC = ", WQTDINC, "NWQTD = ", NWQTD
    WRITE(2,'(/,A5,A10,A15,20A10)') "IT", "TEMP", "WQTDMNL", "WQTDHDR", "WQTDNIT", "WQKSUA", "WQKCOD", "WQTDKR", "WQTDTAM", "WQTD1FCB", "WQTD2FCB"
                                            
    DO M = 1,NWQTD
      WRITE(2,'(I5,F10.3,F15.3,20F10.5)') M, WQTDTEMP(M), WQTDMNL(M), WQTDHDR(M), WQTDNIT(M)  , WQKSUA(M),  WQKCOD(M,1),  WQTDKR(M,1),  WQTDTAM(M), WQTD1FCB(M), WQTD2FCB(M)
    ENDDO
  endif
  
  ! *** Set up look-up table for temperature dependency over WQTDMIN to WQTDMAX
  DO M = 1,NWQTD
    WTEMP = WQTDTEMP(M)
    TT20 = WTEMP - 20.0
    DO I = 1,NWQZ
      WQKCOD(M,I) = WQKCD(I) * EXP( WQKTCOD*(WTEMP-WQTRCOD) )
      WQTDKR(M,I) = WQKTR(I)**TT20
    ENDDO
  ENDDO
  
  ! **************************************************************************
  ! *** Dry Deposition, compute mass loadings, g/day, moles/day & mpn/day
  IF( IWQNPL /= 1 )THEN
    DO L=2,LA
      DO NW=1,NWQV
        ! M. MORTON MODIFIED THE LINE BELOW SO THAT CONSTANT ATMOSPHERIC DEPOSIT
        ! CAN BE ADDED VIA THIS ROUTINE INSTEAD OF CONSTANT NPS INPUT WHICH THE
        ! ORIGINAL CODE CALLED FOR AND WHICH WAS NOT PARTICULARLY USEFUL.
        ! INPUT DATA (XDSL) ARE IN G/M2/DAY AND ARE MULTIPLIED BY THE CELL SURFA
        ! AREA (DXYP) TO GET G/DAY.  ATMOSPHERIC DEPOSITION ONLY ENTERS THRU SUR
        ! LAYER (KC):
        WQWDSL(L,KC,NW) = XDSL(NW) * DXYP(L)
      ENDDO
      WQWDSL(L,KC,ITAM) = XDSL(ITAM) * CONV1
    ENDDO
  ENDIF
  
  ! **************************************************************************
  ! *** Initial conditions or spatial varing parameters
  Call Broadcast_Array(WQV, master_id)

  ! *** APPLY SPATIALLY CONSTAND WQ IC'S AT 1 AND LC
  DO NW = 1,NWQV
    TVARWQ = WQV(1,1,NW)
    DO K = 1,KC
      WQV(LC,K,NW)  = TVARWQ
      WQV(1 ,K,NW)  = TVARWQ
      WQVO(LC,K,NW) = TVARWQ
      WQVO(1 ,K,NW) = TVARWQ
    ENDDO
  ENDDO

  ! *** APPLY SPATIALLY CONSTAND WQ IC'S
  ! *** Initialize nutrient + biota (g C/m^3)
  DO NW = 1,NWQV
    TVARWQ = WQV(1,1,NW)
    DO K = 1,KC
      DO L = 2,LA
        WQV(L,K,NW)  = TVARWQ
        WQVO(L,K,NW) = TVARWQ
      ENDDO
    ENDDO
  ENDDO
  
  DO NAL = 1,NALGAE
    IF( ISMOB(NAL) == 0 )THEN  ! *** ALGAES(NAL).ISMOBILE is only available after calling ALGAECONTROL - DKT
      ! *** Zero concentration if fixed biota
      DO K = 1,KC
        DO L = 2,LA
          WQV(L,K,19+NAL) = 0
        ENDDO
      ENDDO
      ! *** Updated fixed biota
      ! *** Fixed biota can be converted to biomass density by dividing by bottom layer thickness
      TVARWQ = WQV(1,1,19+NAL)  
      DO L = 1,LCM
        WQV(L,KSZ(L),19+NAL)  = TVARWQ
        WQVO(L,KSZ(L),19+NAL) = TVARWQ
      ENDDO
    ENDIF
  ENDDO
  
  ! *** Setup PO4 sorption options
  IF( IWQSRP == 1 )THEN
    ! *** TAM sorption
    O2WQ_ = MAX(WQV(1,1,IDOX), 0.0)
    WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ_), WQV(1,1,ITAM) )
    WQTAMP(1,1) = WQV(1,1,ITAM) - WQTAMD
    WQPO4D(1,1) = WQV(1,1,IP4D) / (1.0 + WQKPO4P*WQTAMP(1,1))
    WQSAD(1,1)  = WQV(1,1,ISAA) / (1.0 + WQKSAP*WQTAMP(1,1))
  ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN
    ! *** Cohesive sorption
    WQPO4D(1,1) = WQV(1,1,IP4D) / (1.0 + WQKPO4P*SEDT(1,1))
    WQSAD(1,1)  = WQV(1,1,ISAA) / (1.0 + WQKSAP*SEDT(1,1))
  ELSE
    ! *** No sortpion
    WQPO4D(1,1) = WQV(1,1,IP4D)
    WQSAD(1,1)  = WQV(1,1,ISAA)
  ENDIF
  
  XWQTAMP = WQTAMP(1,1)
  XWQPO4D = WQPO4D(1,1)
  XWQSAD  = WQSAD(1,1)
  DO K = 1,KC
    WQTAMP(LC,K) = XWQTAMP
    WQPO4D(LC,K) = XWQPO4D
    WQSAD(LC,K)  = XWQSAD
    WQTAMP(1,K)  = XWQTAMP
    WQPO4D(1,K)  = XWQPO4D
    WQSAD(1,K)   = XWQSAD
  ENDDO

  DO K = 1,KC
    DO L = 2,LA
      WQTAMP(L,K) = XWQTAMP
      WQPO4D(L,K) = XWQPO4D
      WQSAD(L,K) = XWQSAD
    ENDDO
  ENDDO
  
  ! *** Initialize biota.  Units are in g C/m^3
  ! *** Fixed biota can be converted to biomass density by dividing by bottom layer thickness
  
  ! *** C47 SPATIALLY/TEMPORALLY CONSTANT BENTHIC FLUXES
  Call Broadcast_Array(WQBFPO4D,   master_id)
  Call Broadcast_Array(WQBFNH4,    master_id)
  Call Broadcast_Array(WQBFNO3,    master_id)
  Call Broadcast_Array(WQBFSAD,    master_id)
  Call Broadcast_Array(WQBFCOD,    master_id)
  Call Broadcast_Array(WQBFO2,     master_id)

  IF( IWQBEN == 0 )THEN
    ! *** Initialize domain
    DO L=2,LA
      WQBFPO4D(L)= WQBFPO4D(1)
      WQBFNH4(L) = WQBFNH4(1)
      WQBFNO3(L) = WQBFNO3(1)
      WQBFSAD(L) = WQBFSAD(1)
      WQBFCOD(L) = WQBFCOD(1)
      WQBFO2(L)  = WQBFO2(1)
    ENDDO
  ENDIF
   
  
  ! *** READ IN MAPPING INFORMATION FOR SPATIALLY-VARYING PARAMETERS (UNIT #7)
  DO K = 1,KC
    DO L = 2,LA
      IWQZMAP(L,K) = 1
    ENDDO
  ENDDO
        
  IF( NWQZ > 1 )THEN
    Allocate(I2D_Global(LCM_Global,KCM))
    I2D_Global = 0
    
    if( process_id == master_id )then
      WRITE(*,'(A)')' WQ: WQWCMAP.INP'
      WRITE(2,'(/,A)')'Reading WQWCMAP.INP - Water quality zone map'
      OPEN(1,FILE='wqwcmap.inp',STATUS='UNKNOWN')
      CALL SKIPCOM(1,'*',2)  ! *** SKIP OVER TITLE AND AND HEADER LINES

      IM = 0
      IJKC = IC_Global*JC_Global*KC

      DO M = 1,IJKC
        READ(1,*,END=1111) I, J, K, IWQZX
        IM = IM + 1
        IF( IJCT_Global(I,J) < 1 .OR. IJCT_Global(I,J) > 8 )THEN
          PRINT*, 'I, J, K, IJCT(I,J) = ', I,J,K,IJCT_Global(I,J)
          CALL STOPP('ERROR!! INVALID (I,J) IN FILE WQWCMAP.INP')
        ENDIF
        L = LIJ_Global(I,J)
        
        I2D_Global(L,K) = IWQZX
      ENDDO
1111  CONTINUE
      IF( IM /= (LA_Global-1)*KC )THEN
        PRINT *, 'WARNING: ALL ACTIVE WATER CELLS SHOULD BE MAPPED FOR WQ PAR.'
        PRINT *, '         NUMBER OF LINES IN FILE WQWCMAP.INP =\ (LA-1)'
      ENDIF
      CLOSE(1)
    endif   ! *** End of master_id block

    Call Broadcast_Array(I2D_Global, master_id)

    ! *** Map to Local Domain
    DO LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        DO K = 1,KC
          IWQZMAP(L,K) = I2D_Global(LG,K)
        ENDDO
      ENDIF
    ENDDO
    DEALLOCATE(I2D_Global)   

  ENDIF
  
  ! READ IN MAPPING INFORMATION FOR SPATIALLY-VARYING BENTHIC FLUXES.
  ! FORMULATED FOR PECONIC BAY DATA WHICH INCLUDES %MUD FOR EACH CELL AS
  ! WELL AS MAPPING TO BOTH MUD AND SAND FLUXES.  SUBROUTINE WQBENTHIC
  ! CONTAINS THE CODE TO INTERPOLATE THE FINAL FLUX FOR THE CELL BASED
  ! ON THE PERCENT MUD AND THE MUD/SAND FLUXES.
    
  IF( IWQBEN == 2 )THEN
    Call AllocateDSI(R2D_Global, LCM_Global, 3, 0.0)
    
    DO K=1,2
      DO L=2,LA
        IBENMAP(L,K)=1
        XBENMUD(L) = 0.50
      ENDDO
    ENDDO
    if( process_id == master_id )then
      WRITE(*,'(A)')' WQ: WQBENMAP.INP'
      WRITE(2,'(/,A)')'Reading WQBENMAP.INP - Benthic flux rate map for user specified benthic fluxes'
      OPEN(1,FILE='wqbenmap.inp',STATUS='UNKNOWN')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES

      REWIND(1)
      CCMRM = '#'
      CALL SKIPCOM(1, CCMRM,2)
      
      R2D_Global(:,1) = -999.
      IM = 0
      IJKC=IC*JC
      DO M=1,IJKC
        READ(1,*,END=1112) I, J, XMUD, IZMUD, IZSAND
        IM = IM + 1
        
        IF( IJCT_Global(I,J) < 1 .OR. IJCT_Global(I,J) > 8 )THEN
          PRINT *, 'I, J, K, IJCT(I,J) = ', I,J,IJCT_Global(I,J)
          CALL STOPP('ERROR!! INVALID (I,J) IN FILE WQBENMAP.INP')
        ENDIF
        L = LIJ_Global(I,J)
        
        R2D_Global(L,1) = XMUD
        R2D_Global(L,2) = IZMUD
        R2D_Global(L,3) = IZSAND
        
      ENDDO
1112  CONTINUE
      IF( IM  /=  (LA_Global-1) )THEN
        PRINT *, 'WARNING: ALL ACTIVE WATER CELLS SHOULD BE MAPPED FOR BENTHIC FLUXES.'
        PRINT *, '         NUMBER OF LINES IN FILE WQBENMAP.INP <> (LA-1)'
        PRINT *, '         LA-1     = ', LA_Global-1
        PRINT *, '         WQBENMAP = ', IM
        PRINT *, '         Details can be found in log_mpi_proc_000'
        
        
        Call WriteBreak(mpi_log_unit)
        write(mpi_log_unit,'(a)') 'Warning.  The following cells have not been mapped in WQBENMAP.INP'
        DO L = 2,LA_Global
          IF( R2D_Global(L,1) == -999. )THEN
            write(mpi_log_unit,'(a,i10,2(a,i8))') 'L = ', L, '  I = ', Map2Global(L).iG,  '  J = ', Map2Global(L).jG
          endif
        ENDDO
        
      ENDIF
      CLOSE(1)
    endif   ! *** End of master_id block

    Call Broadcast_Array(R2D_Global,     master_id)

    ! *** Map to Local Domain
    DO LG=2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        IBENMAP(L,1) = R2D_Global(LG,2)        ! *** IZMUD  - Zone to be used for mud
        IBENMAP(L,2) = R2D_Global(LG,3)        ! *** IZSAND - Zone to be used for sand
        XBENMUD(L)   = R2D_Global(LG,1)/100.0  ! *** XMUD   - Percent mud for interpolation between IZMUD and IZSAND
      ENDIF
    ENDDO
    DEALLOCATE(R2D_Global)   
  ENDIF

  ! *** These are required to initialize WQ variables along boundaries so that concentrations discontinuities can be smoothed.
  DO NW=1,NWQV
    M = MSVWQV(NW)
    IF( M > 0 )THEN
      DO K=1,KC
        DO LL=1,NCBS
          CLOS(LL,K,M) = WQV(1,1,NW)
          NLOS(LL,K,M) = NITER
        ENDDO
        DO LL=1,NCBW
          CLOW(LL,K,M) = WQV(1,1,NW)
          NLOW(LL,K,M) = NITER
        ENDDO
        DO LL=1,NCBE
          CLOE(LL,K,M) = WQV(1,1,NW)
          NLOE(LL,K,M) = NITER
        ENDDO
        DO LL=1,NCBN
          CLON(LL,K,M) = WQV(1,1,NW)
          NLON(LL,K,M) = NITER
        ENDDO
      ENDDO
    ENDIF
  ENDDO
    

  END SUBROUTINE WQ3DCONTROL
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine ALGAECONTROL
  !
  !> @details  READ THE ALGAE CONTROL JNP FILE FOR WATER QUALITY
  !---------------------------------------------------------------------------!
  !    I/O CONTROL VARIABLES
  !    SPATIALLY AND TEMPORALLY CONSTANT PARAMETERS
  !---------------------------------------------------------------------------! 
  SUBROUTINE ALGAECONTROL
  
  USE fson
  USE mod_fson_value, ONLY: fson_value_count, fson_value_get
  
  IMPLICIT NONE
  
  INTEGER :: IZ, M, ITMP, NS, NAL, L, K, ISMOB
  REAL :: WTEMP
  TYPE(fson_value), POINTER :: json_data, phytogroups, item
  
  if( process_id == master_id  .and. NALGAE > 0 )then
    json_data => fson_parse("wq_biota.jnp")
    
    ! *** Get the phytoplankton group's parameters as an array
    CALL fson_get(json_data, "groups", phytogroups)
    
    ! *** Loop through each array item
    DO NAL = 1, fson_value_count(phytogroups)
      ! *** Get the array item 
      item => fson_value_get(phytogroups, NAL)
      
      ! *** Lookup the values from the array
      CALL fson_get(item, "index",                                          ALGAES(NAL).IDN)            ! *** 
      CALL fson_get(item, "mobility_flag",                                  ISMOB)                      ! *** 
      CALL fson_get(item, "salinity_toxicity_flag",                         ALGAES(NAL).ISTOX)          ! *** 
      CALL fson_get(item, "winter_bloom_flag",                              ALGAES(NAL).ISBLOOM)        ! *** 
      CALL fson_get(item, "silica_active_flag",                             ALGAES(NAL).ISILICA)        ! *** 
      CALL fson_get(item, "settling_velocity",                              ALGAES(NAL).WQWS(1))        ! *** 
        
      IF( ISMOB == 0 )THEN
        ALGAES(NAL).ISMOBILE = .FALSE.     
      ELSE
        ALGAES(NAL).ISMOBILE = .TRUE.
      ENDIF
      
      ! *** General parameters
      CALL fson_get(item, "carbon_to_chla_ratio",                           ALGAES(NAL).WQCHLA)         ! *** 
      Call fson_get(item, "light_extinction_due_to_shading",                ALGAES(NAL).WQKEMAC)        ! *** Light_extinction_due_to_shading
      CALL fson_get(item, "O2_to_C_ratio",                                  ALGAES(NAL).WQALGOCR)       ! *** 
      CALL fson_get(item, "O2_to_N_ratio",                                  ALGAES(NAL).WQALGONT)       ! *** 
      CALL fson_get(item, "stoichiometry.nitrogen_to_carbon_ratio",         ALGAES(NAL).WQANCA)         ! *** 
      CALL fson_get(item, "stoichiometry.silica_to_carbon_ratio",           ALGAES(NAL).WQASC)          ! *** 
      CALL fson_get(item, "stoichiometry.factor_to_modify_C_to_P_ratio",    ALGAES(NAL).WQAPCM)         ! *** 
      
      ! *** Growth parameters
      CALL fson_get(item, "growth.max_growth_rate",                         ALGAES(NAL).WQPMA(1))       ! *** 
      CALL fson_get(item, "growth.photosynthesis_O2_to_C_ratio",            ALGAES(NAL).WQAOCRP)        ! *** 
      CALL fson_get(item, "growth.phosphorus_half_saturation",              ALGAES(NAL).WQKHPA)         ! *** 
      CALL fson_get(item, "growth.nitrogen_half_saturation",                ALGAES(NAL).WQKHNA)         ! *** 
      CALL fson_get(item, "growth.silica_half_saturation",                  ALGAES(NAL).WQKHS)          ! *** 
      CALL fson_get(item, "growth.CO2_half_saturation",                     ALGAES(NAL).WQKHCO2)        ! *** 
      CALL fson_get(item, "growth.halved_microsystis_growth_salinity",      ALGAES(NAL).WQSTOX)         ! *** 
      CALL fson_get(item, "growth.optimal_depth",                           ALGAES(NAL).WQDOP(1))       ! *** 
      CALL fson_get(item, "growth.optimal_growth.lower_temperature",        ALGAES(NAL).WQTM1)          ! *** 
      CALL fson_get(item, "growth.optimal_growth.upper_temperature",        ALGAES(NAL).WQTM2)          ! *** 
      CALL fson_get(item, "growth.optimal_growth.lower_coefficient",        ALGAES(NAL).WQKG1)          ! *** 
      CALL fson_get(item, "growth.optimal_growth.upper_coefficient",        ALGAES(NAL).WQKG2)          ! *** 
      CALL fson_get(item, "growth.plant_density_half_saturation_factor",    ALGAES(NAL).WQKBP(1))       ! *** 
      Call fson_get(item, "growth.minimum_biomass",                         ALGAES(NAL).WQBMIN(1))      ! *** 
      Call fson_get(item, "growth.maximum_biomass",                         ALGAES(NAL).WQBMAX)         ! *** 

      ! *** Metabolism parameters
      CALL fson_get(item, "metabolism.reference_temperature",               ALGAES(NAL).WQTR)           ! *** 
      CALL fson_get(item, "metabolism.reference_rate",                      ALGAES(NAL).WQBMRA(1))      ! *** 
      CALL fson_get(item, "metabolism.effect_of_temperature",               ALGAES(NAL).WQKTB)          ! *** 
      CALL fson_get(item, "metabolism.respiration_O2_to_C_ratio",           ALGAES(NAL).WQAOCRR)        ! *** 
      CALL fson_get(item, "metabolism.half_sat_DOC_excretion",              ALGAES(NAL).WQKHRA(1))      ! *** 
      CALL fson_get(item, "metabolism.fraction_of_carbon.DOC",              ALGAES(NAL).WQFCDB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_phosphorus.RPOP",         ALGAES(NAL).WQFPRB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_phosphorus.LPOP",         ALGAES(NAL).WQFPLB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_phosphorus.DOP",          ALGAES(NAL).WQFPDB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_phosphorus.PO4",          ALGAES(NAL).WQFPIB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_nitrogen.RPON",           ALGAES(NAL).WQFNRB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_nitrogen.LPON",           ALGAES(NAL).WQFNLB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_nitrogen.DON",            ALGAES(NAL).WQFNDB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_nitrogen.NH4",            ALGAES(NAL).WQFNIB)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_silica.SU",               ALGAES(NAL).WQFSPD)         ! *** 
      CALL fson_get(item, "metabolism.fraction_of_silica.SA",               ALGAES(NAL).WQFSID)         ! *** 
      
      ! *** Predation parameters
      CALL fson_get(item, "predation.predation_rate", ALGAES(NAL).WQPRRA(1))
      ! *** Set algal death rate if zooplankton is simulated
      IF( IWQZPL == 1  ) ALGAES(NAL).WQDRA(1) = ALGAES(NAL).WQPRRA(1)
        
      CALL fson_get(item, "predation.optimal_predation.lower_temperature",  ALGAES(NAL).WQTP1)          ! *** 
      CALL fson_get(item, "predation.optimal_predation.upper_temperature",  ALGAES(NAL).WQTP2)          ! *** 
      CALL fson_get(item, "predation.optimal_predation.lower_coefficient",  ALGAES(NAL).WQKP1)          ! *** 
      CALL fson_get(item, "predation.optimal_predation.upper_coefficient",  ALGAES(NAL).WQKP2)          ! *** 
      CALL fson_get(item, "predation.fraction_of_carbon.RPOC",              ALGAES(NAL).WQFCRP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_carbon.LPOC",              ALGAES(NAL).WQFCLP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_carbon.DOC",               ALGAES(NAL).WQFCDP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_phosphorus.RPOP",          ALGAES(NAL).WQFPRP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_phosphorus.LPOP",          ALGAES(NAL).WQFPLP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_phosphorus.DOP",           ALGAES(NAL).WQFPDP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_phosphorus.PO4",           ALGAES(NAL).WQFPIP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_nitrogen.RPON",            ALGAES(NAL).WQFNRP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_nitrogen.LPON",            ALGAES(NAL).WQFNLP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_nitrogen.DON",             ALGAES(NAL).WQFNDP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_nitrogen.NH4",             ALGAES(NAL).WQFNIP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_silica.SU",                ALGAES(NAL).WQFSPP)         ! *** 
      CALL fson_get(item, "predation.fraction_of_silica.SA",                ALGAES(NAL).WQFSID)         ! *** 
           
      Call fson_get(item, "growth_velocity_limitation_option",              ALGAES(NAL).IWQVLIM)      
			CALL fson_get(item, "vel_lim.half_saturation_velocity",               ALGAES(NAL).WQKMV(1))       ! *** 
      CALL fson_get(item, "vel_lim.mininum_velocity",                       ALGAES(NAL).WQKMVMIN(1))    ! *** 
			CALL fson_get(item, "vel_lim.lf_param_a",                             ALGAES(NAL).WQKMVA(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_b",                             ALGAES(NAL).WQKMVB(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_c",                             ALGAES(NAL).WQKMVC(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_d",                             ALGAES(NAL).WQKMVD(1))      ! *** 
			CALL fson_get(item, "vel_lim.lf_param_e",                             ALGAES(NAL).WQKMVE(1))      ! *** 

      IF( .NOT. ALGAES(NAL).ISMOBILE )THEN
        Call fson_get(item, "layerthreshold",                               ALGAES(NAL).THRESHOLD)      ! *** Biomass concentration limit to force macrophyte mass into the next layer (mg C/L/m)
        Call fson_get(item, "base_depth",                                   ALGAES(NAL).BASEDEPTH)      ! *** Distance below the surface (i.e. depth) of the "base" of the macrophyte/periphyton growth (m)
        Call fson_get(item, "maximum_length",                               ALGAES(NAL).MAXLENGTH)      ! *** Maximum length from the "base" to allow macrophyte/periphyton growth (m)
        Call fson_get(item, "use_hydro_feedback",                           ALGAES(NAL).ISDRAG)      

        Call fson_get(item, "hydro_feedback.dragcoefficient",               ALGAES(NAL).DRAGCOEFF)      ! *** Vegetagive drag coefficient.  Includes stems and leaves
        Call fson_get(item, "hydro_feedback.stemdiameter",                  ALGAES(NAL).MINDIAMETER)    ! *** Minimum stem diameter. As growth occurs, stemdiameter is updated based on % dry matter, stemdensity and concentration (m)
        Call fson_get(item, "hydro_feedback.stemdensity",                   ALGAES(NAL).STEMDENSITY)    ! *** As growth occurs, stemdensity is constant (# stems/m2)

        ! *** Apply QC and special conditions

      ENDIF      
    ENDDO   
    
  Endif    ! *** End of master_id block
    
  ! *** Broadcast water quality global settings, not dependent on grid.
  
  ! *** Algae's parameters
  Do NAL = 1, NALGAE
     Call Broadcast_Scalar(ALGAES(NAL).IDN,         master_id)
     Call Broadcast_Scalar(ALGAES(NAL).ISMOBILE,    master_id)
     Call Broadcast_Scalar(ALGAES(NAL).ISTOX,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).ISBLOOM,     master_id)
     Call Broadcast_Scalar(ALGAES(NAL).ISILICA,     master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQWS(1),     master_id)
          
     Call Broadcast_Scalar(ALGAES(NAL).WQCHLA,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKEMAC,     master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQALGOCR,    master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQALGONT,    master_id)

     Call Broadcast_Scalar(ALGAES(NAL).WQANCA,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQASC,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQPMA(1),    master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQAOCRP,     master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKHPA,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKHNA,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKHS,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKHCO2,     master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQSTOX,      master_id)  
     Call Broadcast_Scalar(ALGAES(NAL).WQDOP(1),    master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQTM1,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQTM2,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKG1,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKG2,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKBP(1),    master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQBMIN(1),   master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQBMAX,      master_id) 
     Call Broadcast_Scalar(ALGAES(NAL).WQTR,        master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQBMRA(1),   master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKTB,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQAOCRR,     master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKHRA(1),   master_id) 
     Call Broadcast_Scalar(ALGAES(NAL).WQFCDB,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFPRB,      master_id)     
     Call Broadcast_Scalar(ALGAES(NAL).WQFPLB,      master_id)         
     Call Broadcast_Scalar(ALGAES(NAL).WQFPDB,      master_id)          
     Call Broadcast_Scalar(ALGAES(NAL).WQFPIB,      master_id)        
     Call Broadcast_Scalar(ALGAES(NAL).WQFNRB,      master_id)   
     Call Broadcast_Scalar(ALGAES(NAL).WQFNLB,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFNDB,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFNIB,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFSPD,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFSID,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQPRRA(1),   master_id)
     IF( IWQZPL == 1 ) Call Broadcast_Scalar(ALGAES(NAL).WQDRA(1),   master_id)
     
     Call Broadcast_Scalar(ALGAES(NAL).WQTP1,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQTP2,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKP1,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQKP2,       master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFCRP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFCLP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFCDP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFPRP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFPLP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFPDP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFPIP,      master_id)

     Call Broadcast_Scalar(ALGAES(NAL).WQFNRP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFNLP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFNDP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFNIP,      master_id)                                            
     Call Broadcast_Scalar(ALGAES(NAL).WQFSPP,      master_id)
     Call Broadcast_Scalar(ALGAES(NAL).WQFSIP,      master_id)  
     
     Call Broadcast_Scalar(ALGAES(NAL).IWQVLIM,     master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).WQKMV(1),    master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQKMVMIN(1), master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQKMVA(1),   master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQKMVB(1),   master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQKMVC(1),   master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQKMVD(1),   master_id)                                              
     Call Broadcast_Scalar(ALGAES(NAL).WQKMVE(1),   master_id)                                              

     ! *** Macrophyte/Periphyton settings
     Call Broadcast_Scalar(ALGAES(NAL).THRESHOLD,   master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).BASEDEPTH,    master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).MAXLENGTH,   master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).ISDRAG,      master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).DRAGCOEFF,   master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).MINDIAMETER, master_id)                                               
     Call Broadcast_Scalar(ALGAES(NAL).STEMDENSITY, master_id)                                               
  End do
        
  Do NAL = 1, NALGAE
    ALGAES(NAL).WQCHLA = 1.0/(ALGAES(NAL).WQCHLA + 1.E-12)
    ALGAES(NAL).WQSTOX = ALGAES(NAL).WQSTOX*ALGAES(NAL).WQSTOX
    ! *** Zeroing transport bypass flag of macro algae
    IF( .NOT. ALGAES(NAL).ISMOBILE )THEN 
      ISTRWQ(19+NAL) = 0
    END IF   
    
    Call Broadcast_Scalar(ALGAES(NAL).WQCHLA,   master_id)
    Call Broadcast_Scalar(ALGAES(NAL).WQSTOX,   master_id)
  End do
  
  Call Broadcast_Array(ISTRWQ,   master_id)
  
  WQAOCR = ALGAES(1).WQALGOCR
  WQAONT = ALGAES(1).WQALGONT
  Call Broadcast_Scalar(WQAOCR,   master_id)
  Call Broadcast_Scalar(WQAONT,   master_id)
  
  !IF(  IWQBEN > 0  )then         ! *** Deprecated - Multiple alagal classes can interact with silica and no need to duplicate variables
  !  Do NAL = 1, NALGAE
  !    SMDWQANC(NAL) = ALGAES(NAL).WQANCA     
  !    IF( ALGAES(NAL).ISILICA == 1 )THEN    
  !      SMWQASC = ALGAES(NAL).WQASC
  !      DIATOM  = ALGAES(NAL).IDN
  !    END IF
  !  Enddo
  !  Call Broadcast_Array (SMDWQANC,   master_id)
  !  Call Broadcast_Scalar(SMWQASC,    master_id)
  !  Call Broadcast_Scalar(DIATOM,     master_id)
  !End if

  ! *** SET UP LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY OVER -10 OC TO 50 OC
  WQTDMIN =-10
  WQTDMAX = 50
  WTEMP = WQTDMIN
  WQTDINC = (WQTDMAX - WQTDMIN)/NWQTD
  
  DO M = 1, NWQTD
    ! *** Reference temperature for algae growth
    DO NAL = 1, NALGAE
      WQTDG(M,NAL) = 1.
      IF( WTEMP < ALGAES(NAL).WQTM1 )THEN
        WQTDG(M,NAL) = EXP(-ALGAES(NAL).WQKG1*(WTEMP - ALGAES(NAL).WQTM1)*(WTEMP - ALGAES(NAL).WQTM1) )
      ENDIF
      IF( WTEMP > ALGAES(NAL).WQTM2 )THEN
        WQTDG(M,NAL) = EXP(-ALGAES(NAL).WQKG2*(WTEMP - ALGAES(NAL).WQTM2)*(WTEMP - ALGAES(NAL).WQTM2) )
      ENDIF
    ENDDO
    
    ! *** Temperature adjustment for algae metabolism & Predation
    DO NAL = 1, NALGAE 
      WQTDR(M,NAL) = EXP( ALGAES(NAL).WQKTB*(WTEMP - ALGAES(NAL).WQTR) ) 
    ENDDO
    
    ! *** Temperature adjustment for winter bloom.  Configure WQTP1/2 & WQKP1.1 so WQTDP is low/zero at low temps
    DO NAL = 1,NALGAE
      WQTDP(M,NAL) = 1.
      IF( ALGAES(NAL).ISBLOOM > 0 )THEN
        IF( WTEMP < ALGAES(NAL).WQTP1 )THEN
          WQTDP(M,NAL) = EXP(-ALGAES(NAL).WQKP1*(WTEMP - ALGAES(NAL).WQTP1)*(WTEMP - ALGAES(NAL).WQTP1))
        ENDIF
        IF( WTEMP > ALGAES(NAL).WQTP2 )THEN
          WQTDP(M,NAL) = EXP(-ALGAES(NAL).WQKP2*(WTEMP - ALGAES(NAL).WQTP2)*(WTEMP - ALGAES(NAL).WQTP2))
        ENDIF
      ENDIF
    ENDDO
    
    WTEMP = WTEMP + WQTDINC
  ENDDO
  
  if( process_id == master_id )then
    WRITE(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for Biota Temperature Dependency",                               &
                                             "WQTDMIN = ", WQTDMIN, "WQTDMAX = ", WQTDMAX, "WQTDINC = ", WQTDINC, "NWQTD = ", NWQTD
    WRITE(2,'(/,A5,A10,30(" |",I2,3A10))') "IT", "TEMP", (NAL, "WQTDG", "WQTDR", "WQTDP", NAL=1,NALGAE)
                                            
    DO M = 1,NWQTD
      WRITE(2,'(I5,F10.3,30(4X,3F10.5))') M, WQTDTEMP(M), (WQTDG(M,NAL), WQTDR(M,NAL), WQTDP(M,NAL), NAL=1,NALGAE)
    ENDDO
  endif
  
  ! *** C44
  WQCHL(1,1) = 0.0
  IF( ISWQLVL == 0 )THEN
    WQCHL(1,1) = WQV(1,1,1)*ALGAES(1).WQCHLA + WQV(1,1,2)*ALGAES(2).WQCHLA + WQV(1,1,3)*ALGAES(3).WQCHLA
  ELSE
    DO NAL = 1,NALGAE
      IF( ALGAES(NAL).ISMOBILE )THEN
        ! *** Phytoplankton
        WQCHL(1,1) = WQCHL(1,1) + WQV(1,1,19+NAL)*ALGAES(NAL).WQCHLA
      ENDIF
    ENDDO
  ENDIF

  DO K=1,KC
    WQCHL(LC,K) = WQCHL(1,1)
  ENDDO
  DO K = 1,KC
    DO L = 2,LA
      WQCHL(L,K) = WQCHL(1,1)
    ENDDO
  ENDDO
  
  ! *** Read Zonal biota class values
  IF( IWQZONES > 0 )THEN
    ! *** Defaults
    ! *** WQKMV    = 0.25     ! *** KMV velocity half-saturation
    ! *** WQKMVMIN = 0.15     ! *** Minimum velocity
    ! *** WQKBP    = 6.5      ! *** M density half saturation
    ! *** WQKMVA   = 1.0      ! *** Five-parameter logistic - A
    ! *** WQKMVB   = 12.0     ! *** Five-parameter logistic - B
    ! *** WQKMVC   = 0.3      ! *** Five-parameter logistic - C
    ! *** WQKMVD   = 0.25     ! *** Five-parameter logistic - D
    ! *** WQKMVE   = 2.0      ! *** Five-parameter logistic - E
    
    if( process_id == master_id )then
      WRITE(*,'(A)')' WQ: WQ_BIO_ZONES.JNP'
      WRITE(2,'(1/,A)')'Reading WQ_BIO_ZONES.JNP - Zonal kinetics for Biota parameters'
      json_data => fson_parse("wq_bio_zones.jnp")  
      CALL fson_get(json_data, "groups", phytogroups)  
      
      ! *** Macrophytes, periphyton and phytoplankton
      DO NAL = 1, fson_value_count(phytogroups)
          item => fson_value_get(phytogroups, NAL)
          
          ! *** BIO = "index"
          CALL fson_get(item, "algae_dynamic.max_growth_rate",                 ALGAES(NAL).WQPMA)      ! *** 
          CALL fson_get(item, "algae_dynamic.predation(death)_rate",           ALGAES(NAL).WQPRRA)     ! *** 
          CALL fson_get(item, "algae_dynamic.metabolism_rate",                 ALGAES(NAL).WQBMRA)     ! *** 
          CALL fson_get(item, "growth_optimal_depth",                          ALGAES(NAL).WQDOP)      ! *** 
          CALL fson_get(item, "plant_density_half_saturation_factor",          ALGAES(NAL).WQKBP)      ! *** 
          CALL fson_get(item, "settling_velocity",                             ALGAES(NAL).WQWS)       ! *** 
          Call fson_get(item, "minimum_biomass_zone",                          ALGAES(NAL).WQBMIN)

          Call fson_get(item, "hydrolysis_constant_to_macrophytes",            ALGAES(NAL).WQKDCALM)   ! *** 
          CALL fson_get(item, "half_saturation_DOC_excretion",                 ALGAES(NAL).WQKHRA)     ! *** 
          CALL fson_get(item, "velocity_limitation.half_saturation_velocity",  ALGAES(NAL).WQKMV)      ! *** 
          CALL fson_get(item, "velocity_limitation.mininum_velocity",          ALGAES(NAL).WQKMVMIN)   ! *** 

          CALL fson_get(item, "velocity_limitation.logistic_function.param_a", ALGAES(NAL).WQKMVA)     ! *** 
          CALL fson_get(item, "velocity_limitation.logistic_function.param_b", ALGAES(NAL).WQKMVB)     ! *** 
          CALL fson_get(item, "velocity_limitation.logistic_function.param_c", ALGAES(NAL).WQKMVC)     ! *** 
          CALL fson_get(item, "velocity_limitation.logistic_function.param_d", ALGAES(NAL).WQKMVD)     ! *** 
          CALL fson_get(item, "velocity_limitation.logistic_function.param_e", ALGAES(NAL).WQKMVE)     ! *** 
      ENDDO
    endif   ! *** End of master_id block
    
    ! *** Broadcast
    DO NAL = 1, NALGAE
      Call Broadcast_Array(ALGAES(NAL).WQPMA,    master_id)
      Call Broadcast_Array(ALGAES(NAL).WQPRRA,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQBMRA,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQDOP,    master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKBP,    master_id)
      Call Broadcast_Array(ALGAES(NAL).WQWS,     master_id)
      Call Broadcast_Array(ALGAES(NAL).WQBMIN,   master_id)
      
      Call Broadcast_Array(ALGAES(NAL).WQKDCALM, master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKHRA,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKMV,    master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKMVMIN, master_id)
      
      Call Broadcast_Array(ALGAES(NAL).WQKMVA,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKMVB,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKMVC,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKMVD,   master_id)
      Call Broadcast_Array(ALGAES(NAL).WQKMVE,   master_id)
    ENDDO
    
  ENDIF  
  
  END SUBROUTINE ALGAECONTROL
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQICI
  !
  !> @details  READ IN SPATIALLY AND/OR TEMPORALLY VARYING ICS (UNIT INWQICI)
  !---------------------------------------------------------------------------!
  !
  !---------------------------------------------------------------------------! 
  SUBROUTINE WQICI
  
  IMPLICIT NONE

  INTEGER :: M,I,J,NW,L,K
  REAL, SAVE,ALLOCATABLE,DIMENSION(:) :: XWQV
  CHARACTER HEADER(3)*300
  
  IF(  .NOT. ALLOCATED(XWQV) )THEN
    ALLOCATE(XWQV(NWQVM))
    XWQV=0.0
  ENDIF

  WRITE(*,'(A)')' WQ: WQICI.INP'
  WRITE(2,'(/,A)')'Reading WQICI.INP - Reading initial conditions'
  OPEN(1,FILE="WQICI.INP",STATUS='UNKNOWN')

  READ(1,50) (HEADER(M),M=1,3)
  READ(1,999)
  READ(1,50) HEADER(1)
  WRITE(2,50) HEADER(1)
  DO M=2,LA_Global
    READ(1,*) I, J, (XWQV(NW),NW=1,NWQV)
    IF( IJCT_Global(I,J) < 1 .OR. IJCT_Global(I,J) > 8 )THEN
      PRINT *, 'I, J, LINE # = ', I,J,M-1
      CALL STOPP('ERROR!! INVALID (I,J) IN WQICI.INP')
    ENDIF
    
    L = LIJ_Global(I,J)
    DO K=1,KC
      DO NW=1,NWQV
        WQV_Global(L,K,NW) = XWQV(NW)                   ! *** WQV
      ENDDO
    ENDDO
    WRITE(2,84) I, J, (WQV_Global(L,1,NW),NW=1,NWQV)    ! *** WQV
  ENDDO
 
  CLOSE(1)

    999 FORMAT(1X)
     50 FORMAT(A300)
     52 FORMAT(I7, 1X, A3)
     60 FORMAT(/, A24, I5, A24)
     84 FORMAT(2I5, 25E12.4)
  RETURN
  
  END SUBROUTINE WQICI
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQCSR
  !
  !> @details  READ TIME SERIES FILES FOR WQ CONCENTRATIONS
  !
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQCSR

  USE INFOMOD,ONLY:SKIPCOM,READSTR

  IMPLICIT NONE

  CHARACTER*11 FNWQSR(40)
  CHARACTER*2  SNUM
  CHARACTER*80 STR*200
  INTEGER NC,NW,IS,NS,ISO,ISTYP,K,M
  REAL    RMULADJ,ADDADJ,CSERTMP

  WRITE(*,'(A)') ' WQ: READING WQCSRxx.INP - WQ CONCENTRATION TIME SERIES'
  WRITE(2,'(/,A)') 'Reading WQCSRxx.INP - WQ CONCENTRATION TIME SERIES'

  ! *** DEFINE THE INPUT FILE NAMES
  DO NW = 1,NWQV
    WRITE(SNUM,'(I2.2)')NW
    FNWQSR(NW)='wqcsr'//SNUM//'.inp'
  ENDDO
  
  IF( IWQZPL > 0 )THEN
    DO NW = 1,NZOOPL
      WRITE(SNUM,'(I2.2)')NW
      FNWQSR(22+NW)='zoosr'//SNUM//'.inp'
    ENDDO
  ENDIF
  
  ! **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
  ! **  TIME SERIES FROM THE FILES WQCSRNN.INP
  NC = 8
  DO NW=1,NWQV
    IF( NWQCSR(NW) >= 1 )THEN
      OPEN(1,FILE=FNWQSR(NW),STATUS='UNKNOWN')
      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      DO NS=1,NWQCSR(NW)
        MTSCLAST(NS,NC) = 2
        READ(1,*,IOSTAT=ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TACSER(NS,NC), RMULADJ, ADDADJ
        IF( ISO > 0 ) GOTO 900
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) GOTO 900
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSWQ(NS,NW).TIM(M),CSERTMP
            IF( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSWQ(NS,NW).VAL(M,K) = ( RMULADJ*(CSERTMP + ADDADJ) )*WKQ(K)
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSWQ(NS,NW).TIM(M),(TSWQ(NS,NW).VAL(M,K), K=1,KC)
            IF( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSWQ(NS,NW).VAL(M,K) = RMULADJ*( TSWQ(NS,NW).VAL(M,K) + ADDADJ )
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDDO

  GOTO 901

  900 CONTINUE
  WRITE(6,601) NW, NS, M
  CALL STOPP('.')

  901 CONTINUE
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQ TIME SERIES, NWQ,NSER,MDATA = ',3I5)
  602 FORMAT(' READ OF FILES WQCSRNN.INP SUCCESSFUL'/)

  RETURN

  END SUBROUTINE WQCSR

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQCSR
  !
  !> @details  READ TIME SERIES FILES FOR WQ CONCENTRATIONS
  !
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQCSR2

  USE INFOMOD,ONLY:SKIPCOM,READSTR

  IMPLICIT NONE

  CHARACTER*15 FNWQSR(40)
  CHARACTER*2  SNUM
  CHARACTER*80 STR*200
  INTEGER :: NC,NW,IS,NS,ISO,ISTYP,K,M
  REAL ::   RMULADJ,ADDADJ,CSERTMP
  
  WRITE(*,'(A)') ' WQ: READING WQ CONCENTRATION TIME SERIES'
  
  ! *** DEFINE THE INPUT FILE NAMES
  ! *** Nutrient
  DO NW = 1,19
    WRITE(SNUM,'(I2.2)') NW
    FNWQSR(NW) = 'wqcsr'//SNUM//'.inp'
  ENDDO
  ! *** Algae
  DO NW = 1,NALGAE
    WRITE(SNUM,'(I2.2)')NW
    FNWQSR(NW+19)='wqalgsr'//SNUM//'.inp'
  ENDDO
  ! *** Zooplankton
  DO NW = 1,NZOOPL
     WRITE(SNUM,'(I2.2)')NW
     FNWQSR(NW+19+NALGAE)='wqzoosr'//SNUM//'.inp'
  ENDDO

  ! **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
  ! **  TIME SERIES FROM THE FILES WQCSRNN.INP
  NC = 8
  DO NW=1,NWQV
    IF( NWQCSR(NW) >= 1 )THEN
      OPEN(1,FILE=FNWQSR(NW),STATUS='UNKNOWN')
      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      DO NS=1,NWQCSR(NW)
        MTSCLAST(NS,NC) = 2
        READ(1,*,IOSTAT=ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TACSER(NS,NC), RMULADJ, ADDADJ
        IF( ISO > 0 ) GOTO 900
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) GOTO 900
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSWQ(NS,NW).TIM(M),CSERTMP
            IF( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSWQ(NS,NW).VAL(M,K) = ( RMULADJ*(CSERTMP + ADDADJ) )*WKQ(K)
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSWQ(NS,NW).TIM(M),(TSWQ(NS,NW).VAL(M,K), K=1,KC)
            IF( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TACSER(NS,NC)
            DO K=1,KC
              TSWQ(NS,NW).VAL(M,K) = RMULADJ*( TSWQ(NS,NW).VAL(M,K) + ADDADJ )
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDDO

  GOTO 901

  900 CONTINUE
  WRITE(6,601) NW, NS, M
  CALL STOPP('.')

  901 CONTINUE
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQ TIME SERIES, NWQ,NSER,MDATA = ',3I5)
  602 FORMAT(' READ OF FILES WQCSRNN.INP SUCCESSFUL'/)

  RETURN

  END SUBROUTINE WQCSR2
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSUN
  !
  !> @details READ IN TEMPORALLY VARYING PARAMETERS FOR DAILY SOLAR RADIATION (WQI0)
  ! AND FRACTIONAL DAYLENGTH (WQFD) (UNIT INWQSUN).
  !
  !---------------------------------------------------------------------------!
  ! NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
  ! READS AND INTERPOLATES DAILY AVERAGE SOLAR RADIATION AND DAYLIGHT FRACTION
  !
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSUN
  
  USE INFOMOD,ONLY:SKIPCOM,READSTR

  IMPLICIT NONE

  INTEGER :: IS, M, ISO, NSUNDAY, M1, M2
  REAL    :: TCSUNDAY, TASUNDAY, RMULADJ, ADDADJ, TDIFF, WTM1, WTM2, TIME
  CHARACTER*80 STR*200

  IF( ITNWQ > 0 ) GOTO 1000
  !
  ! **  READ IN DAILY AVERAGE SOLAR SW RAD SERIES FROM FILE 'SUNDAY.INP'
  !
  WRITE(*,'(A)')' WQ: SUNDAY.INP'
  OPEN(1,FILE='sunday.inp',STATUS='UNKNOWN')

  STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
  M=0
  ISPAR=1
  !
  !      MCSUNDAY=1 TEMP USE ISPAR FOR MCSUNDAY
  !
  READ(1,*,IOSTAT=ISO)NSUNDAY,TCSUNDAY,TASUNDAY,RMULADJ,ADDADJ
  IF( ISO > 0 ) GOTO 900
  DO M=1,NSUNDAY
    READ(1,*,IOSTAT=ISO)TSSRD(M),SOLSRD(M),SOLFRD(M)
    IF( ISO > 0 ) GOTO 900
    TSSRD(M)=TCSUNDAY*( TSSRD(M)+TASUNDAY )
    SOLSRD(M)=RMULADJ*(SOLSRD(M)+ADDADJ) * PARADJ
  ENDDO
  CLOSE(1)
  GOTO 901
    900 CONTINUE
  WRITE(6,601)M
  CALL STOPP('.')
    901 CONTINUE
      1 FORMAT(120X)
    601 FORMAT(' READ ERROR FILE SUNDAY.INP ')
   1000 CONTINUE
  !
  ! **  DAILY AVERAGE SOLAR SW RADIATION INTERPOLTATION FOR WATER QUALITY
  !
  TIME=TIMESEC/86400.

  M1=ISPAR
  !
  !      TEMP USE ISPAR FOR MCSUNDAY
  !
    100 CONTINUE
  M2=M1+1
  IF( TIME > TSSRD(M2) )THEN
    M1=M2
    GOTO 100
  ELSE
    ISPAR=M1
  !
  !      TEMP USE ISPAR FOR MCSUNDAY
  !
  ENDIF
  TDIFF=TSSRD(M2)-TSSRD(M1)
  WTM1=(TSSRD(M2)-TIME)/TDIFF
  WTM2=(TIME-TSSRD(M1))/TDIFF
  SOLSRDT=WTM1*SOLSRD(M1)+WTM2*SOLSRD(M2)
  SOLFRDT=WTM1*SOLFRD(M1)+WTM2*SOLFRD(M2)
  RETURN
  
  END SUBROUTINE WQSUN
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQPSL
  !
  !> @details READ IN TEMPORALLY VARYING POINT SOURCE INPUT (UNIT INWQPSL)
  !
  !---------------------------------------------------------------------------!
  ! Merged SNL & DSI
  !---------------------------------------------------------------------------!

  SUBROUTINE WQPSL

  ! *** INPUT UNITS (KG/D) EXCEPT:  TAM(KMOL/D), FCB(MPN/D).
  ! *** COMPUTATIONAL UNITS, WQ CONSTITUENT LOADS ARE IN G/DAY,
  ! ***                      EXCEPT TAM IN (MOL/D) & FCB IN (MPN/D).

  USE INFOMOD,ONLY:SKIPCOM,READSTR
  
  IMPLICIT NONE

  INTEGER :: IS,NS,ISO,M,NW,M2,M1,K,L,ITMP,KK
  REAL    :: TAWQPSR,RMULADJ,ADDADJ,TIME,TDIFF,WTM1,WTM2
  REAL    :: RLDTMP(NTSWQVM)
  CHARACTER*80 STR*200

  IF( ITNWQ > 0 ) GOTO 1000

  ! **  READ IN LOADING SERIES FROM FILE 'WQPSL.INP'
  IF( NPSTMSR >= 1 )THEN
      
    ! *** Read only on master
    if(process_id == master_id )then

      WRITE(*,'(A)')' WQ: WQPSL.INP'
      OPEN(1,FILE="WQPSL.INP",STATUS='UNKNOWN')
      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      DO NS=1,NPSTMSR
        MWQPTLT(NS)=2 
        READ(1,*,IOSTAT=ISO)MWQPSR(NS),TCWQPSR(NS),TAWQPSR,RMULADJ,ADDADJ
        IF( ISO > 0 ) GOTO 900

        ! *** CONVERT WQ VAR 1-19, 22 FROM KG/D TO G/D  !VB
        ! *** CONVERT WQ VAR 20 (TAM) FROM KMOLS/D TO MOLES/D
        ! *** CONVERT FECAL COLIFORM FROM MPN/DAY TO MPN/D FOR FCM
        RMULADJ=1000.*RMULADJ
        !ADDADJ=ADDADJ

        DO M=1,MWQPSR(NS)
          READ(1,*,IOSTAT=ISO)TWQPSER(M,NS),(RLDTMP(NW),NW=1,7)
          IF( ISO > 0 ) GOTO 900
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=8,14)
          IF( ISO > 0 ) GOTO 900
          READ(1,*,IOSTAT=ISO)(RLDTMP(NW),NW=15,NWQV)  ! PMC HARDWIRED FOR TENKILLER
          IF( ISO > 0 ) GOTO 900

          ! *** STANDARD CONVERSIONS
          TWQPSER(M,NS)=TWQPSER(M,NS)+TAWQPSR
          DO NW=1,NWQV
            WQPSSER(M,NW,NS)=RMULADJ*RLDTMP(NW)
          ENDDO
          WQPSSER(M,IFCB,NS) = WQPSSER(M,IFCB,NS)/1000.
        ENDDO
      ENDDO
      CLOSE(1)
                
    end if ! *** end read on master
    
    ! *** broadcast to all processes
    Call Broadcast_Array(MWQPSR , master_id)
    Call Broadcast_Array(TCWQPSR, master_id)
    Call Broadcast_Scalar(TAWQPSR, master_id)
    Call Broadcast_Array(TWQPSER, master_id)
    Call Broadcast_Array(WQPSSER, master_id)
    Call Broadcast_Scalar(NPSTMSR, master_id)
    call Broadcast_Array(MWQPTLT, master_id)
  ENDIF
  
  GOTO 901


  900 CONTINUE
  WRITE(6,601)NS,M
  CALL STOPP('.')

  901 CONTINUE
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQPS TIME SERIES, NSER,MDATA = ',2I5)
  602 FORMAT(' READ OF FILE WQPSL.INP SUCCESSFUL'/)
      
      
  ! *** BEGINNING OF UPDATING MASS LOADING TO THE CURRENT TIME
  1000 CONTINUE

  ! **  INITIALIZE NULL SERIES LOADING TO ZERO
  DO NW=1,NWQV
    WQPSSRT(NW,0)=0.
  ENDDO

  ! **  LOADING SERIES INTERPOLTATION
  TIME=TIMESEC/86400.

  DO NS=1,NPSTMSR
    TIME=TIMESEC/TCWQPSR(NS)

    M2=MWQPTLT(NS)
    DO WHILE (TIME > TWQPSER(M2,NS))
      M2=M2+1
      IF( M2 > MWQPSR(NS) )THEN
        M2=MWQPSR(NS)
        EXIT
      ENDIF
    END DO
    MWQPTLT(NS)=M2  
    M1 = M2-1
    TDIFF=TWQPSER(M2,NS)-TWQPSER(M1,NS)
    WTM1=(TWQPSER(M2,NS)-TIME)/TDIFF
    WTM2=(TIME-TWQPSER(M1,NS))/TDIFF
    DO NW=1,NWQV
      WQPSSRT(NW,NS) = WTM1*WQPSSER(M1,NW,NS) + WTM2*WQPSSER(M2,NW,NS)
    ENDDO
  ENDDO

  if(process_id == master_id )then
      IF( ITNWQ == 0 .AND. DEBUG )THEN
        OPEN(1,FILE=OUTDIR//'WQPSLT.DIA',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=OUTDIR//'WQPSLT.DIA',STATUS='UNKNOWN')
        WRITE(1,112)NITER,TIME
        DO NS=1,NPSTMSR
          WRITE(1,111)NS,(WQPSSRT(NW,NS),NW=1,NWQV)
        ENDDO
        CLOSE(1)
      ENDIF
  end if
  
  ! **  COMBINE CONSTANT AND TIME VARIABLE PS LOADS
  ! M.R. MORTON 02/20/1999
  ! MODIFIED SO MULTIPLE POINT SOURCES CAN BE ADDED TO ANY GRID CELL
  ! AND ANY LAYER (HAD TO CHANGE WQWPSL ARRAY FROM 2D TO 3D).
  
  IF( ITNWQ == 0 )THEN
    DO NW=1,NWQV
      DO K=1,KC
        DO L=2,LA
          WQWPSL(L,K,NW) = 0.0
        ENDDO
      ENDDO
    ENDDO
  
    if(process_id == master_id )then
        OPEN(1,FILE=OUTDIR//'WQPSL.DIA',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=OUTDIR//'WQPSL.DIA',STATUS='UNKNOWN')
        WRITE(1,112)NITER,TIME
    end if
    
  ENDIF

  ! *** ZERO THE ACTIVE BOUNDARY CELLS
  DO NS=1,NWQPS
    L = LIJ(ICPSL(NS), JCPSL(NS))
    K = KCPSL(NS)
    ! *** DSI BEGIN BLOCK
    IF( K >= 1 )THEN
      DO NW=1,NWQV
        WQWPSL(L,K,NW) = 0.0
      ENDDO
    ELSE
      DO K=1,KC
        DO NW=1,NWQV
          WQWPSL(L,K,NW) = 0.0
        ENDDO
      ENDDO
    ENDIF 
    ! *** DSI END BLOCK
  ENDDO

  ! *** LOOP OVER THE WQ BOUNDARY CELLS
  DO NS=1,NWQPS
    L = LIJ(ICPSL(NS), JCPSL(NS))
    K = KCPSL(NS)
    ITMP = MVPSL(NS)
    ! *** write out some debug code?
    if(process_id == master_id )then
      IF( ITNWQ == 0 ) WRITE(1,121) NS, L, ICPSL(NS), JCPSL(NS), K, ITMP
    end if
    
    IF( K > 0 )THEN
      ! *** K>0, ASSIGN A SPECIFIC LAYER
      IF( K < KSZ(L) )K=KSZ(L)  ! *** FORCE TO A VALID LAYER
      DO NW=1,NWQV
        WQWPSL(L,K,NW) = WQWPSL(L,K,NW) + WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP)
      ENDDO
    ELSE
      ! *** K=0, DISTRIBUTE OVER ALL THE LAYERS
      DO KK=KSZ(L),KC
        DO NW=1,NWQV
          WQWPSL(L,KK,NW) = WQWPSL(L,KK,NW) + DZC(L,KK)*( WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP) )
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  if(process_id == master_id )then
    IF( ITNWQ == 0 )THEN
      DO L=2,LA
        ITMP=IWQPSC(L,1)
        IF( ITMP > 0 )THEN
          DO K=1,KC
            WRITE(1,110) ITMP, IL(L), JL(L), K, (WQWPSL(L,K,NW),NW=1,NWQV)
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  end if
  
  110 FORMAT(1X,4I4,2X,7E12.4,/,19X,7E12.4,/,19X,20E12.4)
  111 FORMAT(1X,I4,2X,7E12.4,/,7X,7E12.4,/,7X,20E12.4)
  112 FORMAT(' N, TIME = ', I10, F12.5/)
  121 FORMAT(' NS,L,I,J,K,ITMP = ', 6I5/)

  RETURN
  
  END SUBROUTINE WQPSL

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQBENTHIC(TIMTMP)
  !
  !> @details  READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR BENTHIC
  !            FLUXES OF PO4D, NH4, NO3, SAD, COD, O2
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQBENTHIC(TIMTMP)

  ! FORMAT OF BENFN FILE IS:
  !    TITLE 1
  !    TITLE 2
  !    TITLE 3
  !  270.00000  <-- DAY AT WHICH FOLLOWING FLUXES BECOME ACTIVE
  !  350.00000  <-- DAY AT WHICH FOLLOWING FLUXES BECOME ACTIVE
  ! 9999.99999 <-- ENTER LARGE DAY AT END OF FILE

  USE INFOMOD,ONLY:SKIPCOM
 
  IMPLICIT NONE
  
  REAL,INTENT(IN) :: TIMTMP
  INTEGER :: M,MM,IBENZ,I,L,IZM,IZS
  REAL    :: XBSFAD,BDAY,XM
  
  CHARACTER TITLE(3)*79, CCMRM*1
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: IZONE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: XBFCOD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: XBFNH4
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: XBFNO3
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: XBFO2
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: XBFPO4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: XBFSAD

  IF( .NOT. ALLOCATED(IZONE) )THEN
    ALLOCATE(IZONE(0:NSMZM))
    ALLOCATE(XBFCOD(0:NSMZM))
    ALLOCATE(XBFNH4(0:NSMZM))
    ALLOCATE(XBFNO3(0:NSMZM))
    ALLOCATE(XBFO2(0:NSMZM))
    ALLOCATE(XBFPO4D(0:NSMZM))
    ALLOCATE(XBFSAD(0:NSMZM))
    IZONE=0
    XBFCOD=0.0
    XBFNH4=0.0
    XBFNO3=0.0
    XBFO2=0.0
    XBFPO4D=0.0
    XBSFAD=0.0
  ENDIF

  if( process_id == master_id )then
    WRITE(*,'(A)')' WQ: WQBENFLX.INP'
    OPEN(1,FILE='WQBENFLX.INP',STATUS='UNKNOWN')

    ! SKIP OVER THREE HEADER RECORDS:
    READ(1,50) (TITLE(M),M=1,3)
    WRITE(2,999)
    WRITE(2,50) (TITLE(M),M=1,3)

    ! SKIP OVER ALL COMMENT CARDS AT BEGINNING OF FILE:
    REWIND(1)
    CCMRM = '#'
    CALL SKIPCOM(1, CCMRM,2)
    READ(1, *) IBENZ
    
    WRITE(2, 65) TIMTMP, IBENZ
    65 FORMAT(' * BENTHIC FLUXES AT     ', F10.5,' DAYS OF MODEL RUN',/, & 
              '   NUMBER OF BENTHIC FLUX ZONES = ', I4)

    ! SEQUENTIALLY READ THROUGH BENTHIC FLUX FILE UNTIL THE APPROPRIATE
    ! TIME IS FOUND:
    !   BDAY   = CURRENT DAY AT WHICH BENTHIC FLUX IS IN EFFECT
    !   BENDAY = NEXT DAY AT WHICH BENTHIC FLUX CHANGES (PASSED TO MAIN PROG

    10 READ(1, *, END=15) BENDAY

    IF( BENDAY  >  TIMTMP ) GOTO 20
    BDAY = BENDAY
    DO I=1,IBENZ
      READ(1,*,END=15) MM, XBFPO4D(MM), XBFNH4(MM), XBFNO3(MM), XBFSAD(MM), XBFCOD(MM), XBFO2(MM)
      IZONE(I) = MM
    ENDDO
    
    ! *** Loop back to read the next day
    GOTO 10 

    ! *** UNEXPECTED END-OF-FILE ENCOUNTERED:
    15 WRITE(2,16) 'WQBENFLX.INP'
    16 FORMAT(//,' ************* WARNING *************',/, &
                 ' END-OF-FILE ENCOUNTERED IN FILE: ', A20,/,/ &
                 ' BENTHIC FLUXES SET TO VALUES CORRESPONDING TO LAST DAY IN FILE.',/)
    BENDAY = (TCON*TBEGIN + NTC*TIDALP)/86400.0  ! *** PMC SINGLE LINE

  endif   ! *** End of master_id block

  ! *** Prepare the model for this data block
20 CONTINUE
   
  Call Broadcast_Scalar(IBENZ,    master_id)
  Call Broadcast_Scalar(BENDAY,   master_id)
  
  Call Broadcast_Array(IZONE,     master_id)
  Call Broadcast_Array(XBFPO4D,   master_id)
  Call Broadcast_Array(XBFNH4,    master_id)
  Call Broadcast_Array(XBFNO3,    master_id)
  Call Broadcast_Array(XBFSAD,    master_id)
  Call Broadcast_Array(XBFCOD,    master_id)
  Call Broadcast_Array(XBFO2,     master_id)
 
  if( process_id == master_id )then
    WRITE(2, 48) BDAY
    48 FORMAT(/,' DAY IN BENTHIC FLUX FILE: ',F10.5,/, &
                '    ZONE    FPO4    FNH4    FNO3    FSAD    FCOD    FSOD')
    DO I=1,IBENZ
      MM = IZONE(I)
      WRITE(2,51) MM, XBFPO4D(MM), XBFNH4(MM), XBFNO3(MM), XBFSAD(MM), XBFCOD(MM), XBFO2(MM)
    ENDDO

    CLOSE(1)
  endif
  
  ! DETERMINE BENTHIC FLUX FOR EACH CELL (L) BY INTERPOLATING BETWEEN
  ! THE MUD AND SAND FLUXES.  XBENMUD(L) IS THE PERCENT MUD FOR EACH CELL.
  DO L=2,LA
    IZM = IBENMAP(L,1)
    IZS = IBENMAP(L,2)
    XM = XBENMUD(L)
    WQBFPO4D(L) = XM*XBFPO4D(IZM) + (1.0-XM)*XBFPO4D(IZS)
    WQBFNH4(L)  = XM*XBFNH4(IZM)  + (1.0-XM)*XBFNH4(IZS)
    WQBFNO3(L)  = XM*XBFNO3(IZM)  + (1.0-XM)*XBFNO3(IZS)
    WQBFSAD(L)  = XM*XBFSAD(IZM)  + (1.0-XM)*XBFSAD(IZS)
    WQBFCOD(L)  = XM*XBFCOD(IZM)  + (1.0-XM)*XBFCOD(IZS)
    WQBFO2(L)   = XM*XBFO2(IZM)   + (1.0-XM)*XBFO2(IZS)
  ENDDO
  
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)

  RETURN

  END SUBROUTINE WQBENTHIC

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQBENTHIC(TIMTMP)
  !
  !> @details  COMPUTES WET ATMOSPHERIC DEPOSITION USING CONSTANT CONCENTRATIONS
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQWET

  ! ** FOR THE 22 STATE VARIABLES MULTIPLIED BY THE RAINFALL FLOW RATE      !VB CHANGED 21 TO 22
  ! ** ENTERING EACH GRID CELL.  COMPUTED LOADS ARE IN G/DAY.

  IMPLICIT NONE

  INTEGER :: L,NW
  REAL    :: CV2,TIME

  !  CV2 = CONVERSION TO GET UNITS OF G/DAY
  !  WQATM(NW) HAS UNITS OF MG/L
  !  RAINT(L) HAS UNITS OF M/SEC
  !  DXYP(L) HAS UNITS OF M2
  !  WQATML(L,KC,NW) HAS UNITS OF G/DAY

  CV2 = 86400.0
  DO NW=1,NWQV
    DO L=2,LA
      WQATML(L,KC,NW) = WQATM(NW,1)*RAINT(L)*DXYP(L)*CV2
    ENDDO
  ENDDO
  
  IF( ITNWQ == 0 .AND. DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'WQATM.DIA',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'WQATM.DIA',STATUS='UNKNOWN')
    TIME=TIMESEC/86400.

    WRITE(1,112) NITER,TIME
    DO L=2,LA
      WRITE(1,110) IL(L),JL(L),(WQATML(L,KC,NW),NW=1,NWQV)
    ENDDO
    CLOSE(1)
  ENDIF
  
    110 FORMAT(1X,2I4,2X,1P,7E11.3,/,15X,7E11.3,/,15X,7E11.3)
    112 FORMAT('# WET ATMOSPHERIC DEPOSITION DIAGNOSTIC FILE',/, &
      ' N, TIME = ', I10, F12.5/)
  RETURN
  
  END SUBROUTINE WQWET
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSKE0
  !
  !> @details  Solve Kinetic Eq from K=KC (surface layer) to K=1 (bottom).
  !            Simplified version that only updates:
  !            IPARAM: 09 Dissolved Organic Phosphorus,
  !                    14 Ammonia Nitrogen
  !                    19 Dissolved Oxygen
  !           After computing new values, store WQVO+WQV into WQVO(L,K,NWQV)
  !           NWQV=15,19,21.
  !---------------------------------------------------------------------------!
  !  ORGINALLY CODED BY K.-Y. PARK
  !---------------------------------------------------------------------------!
  !
  !  1) CHC - cyanobacteria
  !  2) CHD - diatom algae
  !  3) CHG - green algae
  !  4) ROC - refractory particulate organic carbon
  !  5) LOC - labile particulate organic carbon
  !  6) DOC - dissolved organic carbon
  !  7) ROP - refractory particulate organic phosphorus
  !  8) LOP - labile particulate organic phosphorus
  !  9) DOP - dissolved organic phosphorus
  ! 10) P4D - total phosphate
  ! 11) RON - refractory particulate organic nitrogen
  ! 12) LON - labile particulate organic nitrogen
  ! 13) DON - dissolved organic nitrogen
  ! 14) NHX - ammonia nitrogen
  ! 15) NOX - nitrate nitrogen
  ! 16) SUU - particulate biogenic silica
  ! 17) SAA - dissolved available silica
  ! 18) COD - chemical oxygen demand
  ! 19) DOX - dissolved oxygen
  ! 20) TAM - total active metal
  ! 21) FCB - fecal coliform bacteria
  ! 22) CO2 - DISSOLVED CO2
  ! 22) macroalgae

  SUBROUTINE WQSKE0
    
  INTEGER :: L,K,IZ,NS,NAL
  REAL    :: CNS1,TIMTMP,WQOBT0T,XMRM,WQKD0C,O2WQ_,WQTT1,WQKHN,RNH4NO3,WQTTA
  REAL    :: WINDREA,WQWREA,UMRM,VMRM,YMRM,WQD6,WQR6,WQVREA

  
  
  Call STOPP('BAD KINETICS OPTION: ISWQLVL = 0 Does not work!  Contact DSI')   ! delme
  
  
  
  CNS1 = 2.718
  NS = 1
  DO L = 2,LA
    WQI0BOT(L) = WQI0
  ENDDO
  !
  DO K = KC,1,-1
  !
    DO L = 2,LA
      TWQ(L) = MAX(TEM(L,K), 0.0)
      SWQ(L) = MAX(SAL(L,K), 0.0)
      !DZWQ(L) = 1.0 / (HPK(L,K))   deprecated 10.4
      VOLWQ(L) = HPKI(L,K) / DXYP(L)
      IMWQZT(L) = IWQZMAP(L,K)
    ENDDO
  !
  ! FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY
  !
    DO L = 2,LA
      IWQT(L) = NINT((TWQ(L)-WQTDMIN)/WQTDINC)  ! *** DSI SINGLE L!INE
      IF( IWQT(L) < 1 .OR. IWQT(L) > NWQTD )THEN
        if( process_id == master_id )THEN  
          TIMTMP = TIMESEC/86400.
          
          OPEN(3,FILE=OUTDIR//'ERROR.LOG',POSITION='APPEND' &
                ,STATUS='UNKNOWN')
          WRITE(3,*)' *** ERROR IN WATER QUALITY'
          WRITE(3,911) TIMTMP, L, IL(L), JL(L), K, TWQ(L)
          CLOSE(3)
          WRITE(6,600)IL(L),JL(L),K,TWQ(L)
          IWQT(L) = MAX(IWQT(L),1)
          IWQT(L) = MIN(IWQT(L),NWQTD)
  !             CALL STOPP('ERROR!! INVALID WATER TEMPERATURE')
        end if
      ENDIF
    ENDDO
    600 FORMAT(' I,J,K,TEM = ',3I5,E13.4)
    911 FORMAT(/,'ERROR: TIME, L, I, J, K, TWQ(L) = ', F10.5, 4I4, F10.4)
  !
    DO L = 2,LA
      IZ = IWQZMAP(L,K)
  !
  ! UPDATE SOLAR RADIATION AT BOTTOM OF THIS LAYER
  !
      !WQF2IM = WQF2IM * PSHADE(L)      PMC
  !
  ! ALGAL BASAL METABOLISM & PREDATION
  !
      WQBM(L,1) = ALGAES(1).WQBMRA(IMWQZT(L)) * WQTDR(IWQT(L),1)            ! *** Basal metabolism temperature adjustment
      IF( IWQZPL == 0  )THEN
        WQPR(L,1) = ALGAES(1).WQPRRA(IMWQZT(L)) * WQTDR(IWQT(L),1)          ! *** Basal metabolism temperature adjustment - Zooplankton not simulated
      ELSE
        WQPR(L,1) = ALGAES(1).WQDRA(1)                                      ! *** Death rate when zooplankton is used - Zooplankton simulated
      ENDIF
  
      ! *** THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A
      ! *** LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO
      ! *** BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.
      WQOBT0T = WQVO(L,K,20) + WQVO(L,K,21) + WQVO(L,K,22)
      WQKRPC(L) = (WQKRC + WQKRCALG*WQOBT0T) * WQTDHDR(IWQT(L))   ! *** Hydrolysis --> DOC
      WQKLPC(L) = (WQKLC + WQKLCALG*WQOBT0T) * WQTDHDR(IWQT(L))   ! *** Hydrolysis --> DOC
      XMRM = 0.0
      DO NAL = 1, NALGAE
        IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
          ! *** Macrophytes and periphyton
          XMRM = XMRM + ALGAES(NAL).WQKDCALM(1) * WQVO(L,K,19+NAL)
        ENDIF
      ENDDO

      ! M. MORTON 08/28/99: ADDED SPATIALLY VARIABLE DOC HYDROLYSIS RATE WQKDC
      !    TO ACHIEVE BETTER CONTROL IN SYSTEMS WITH A COMBINATION OF FRESHWAT
      !    STREAMS AND TIDAL RIVERS WITH DIFFERENT CHARACTERISTICS.
      WQKD0C = (WQKDC(1) + WQKDCALG*WQOBT0T + XMRM)*WQTDMNL(IWQT(L))
      O2WQ_ = MAX(WQVO(L,K,16), 0.0)
      WQTT1 = WQKD0C / (WQKHORDO + O2WQ_+ 1.E-18)
      WQKHR(L) = WQTT1 * O2WQ_
      WQDENIT(L) = 0.0
  !
  ! 7-10 PHOSPHORUS
  ! 11-15 NITROGEN
  !
      WQKHN = (ALGAES(1).WQKHNA + ALGAES(2).WQKHNA + ALGAES(3).WQKHNA) / 3.0
      RNH4WQ_ = MAX (WQVO(L,K,11), 0.0)
      RNO3WQ_ = MAX (WQVO(L,K,12), 0.0)
      RNH4NO3_ = RNH4WQ_ + RNO3WQ_
      WQTT1 = WQKHN / (WQKHN+RNH4NO3_+ 1.E-18) * WQOBT0T
      
      IF( RNH4NO3_ == 0.0 )THEN
        DO NAL = 1,NALGAE
          WQPN(L,NAL) = 0.0      
        ENDDO
      ELSE
        DO NAL = 1,NALGAE
          WQTTA = RNH4WQ_/(ALGAES(NAL).WQKHNA+RNO3WQ_+ 1.E-18)
          WQPN(L,NAL) = (RNO3WQ_/(ALGAES(NAL).WQKHNA+RNH4WQ_+ 1.E-18) + ALGAES(NAL).WQKHNA/(RNH4NO3_+ 1.E-18)) * WQTTA
        ENDDO
      ENDIF
      WQNIT(L) = WQTDNIT(IWQT(L)) * (O2WQ(L) / (WQKHNDO + O2WQ(L) + 1.E-18)) * (RNH4WQ_ / (WQKHNN + RNH4WQ_ + 1.E-18))
        
      WQDOS(L) = DO_SAT(L)
      XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*HPK(L,K)
      IF( K == KC )THEN
  !
  ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION:
  !
        WINDREA = WINDST(L)
        WQWREA = 0.728*SQRT(WINDREA) + (0.0372*WINDREA-0.317)*WINDREA
  !
  !        WQWREA = 0.728*SQRT(WINDST(L))
  !
        IF( IWQKA(IZ)  ==  0 )THEN
          WQVREA = WQKRO(IZ)
          WQWREA = 0.0
        ENDIF
  !
  !                 WIND VELOCITY COMPUTED ABOVE:
  !
        IF( IWQKA(IZ)  ==  1 )THEN
          WQVREA = WQKRO(IZ)
        ENDIF
  !
  !    WQKRO = 3.933 TYPICALLY
  !
        IF( IWQKA(IZ)  ==  2 )THEN
          UMRM = 0.5*(U(L,K)+U(LEC(L),K))
          VMRM = 0.5*(V(L,K)+V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**0.5
        ENDIF
  !
  !    WQKRO = 5.32 TYPICALLY
  !
        IF( IWQKA(IZ)  ==  3 )THEN
          UMRM = MAX(U(L,K), U(LEC(L),K))
          VMRM = MAX(V(L,K), V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85
        ENDIF
  !
  ! MODIFIED OWENS AND GIBBS REAERATION EQUATION:
  ! NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE
  !       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER
  !       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.
  ! WQKRO = 5.32 TYPICALLY
  !
        IF( IWQKA(IZ)  ==  4 )THEN
          UMRM = MAX(U(L,K), U(LEC(L),K))
          VMRM = MAX(V(L,K), V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L) + 0.1524))
          WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85
        ENDIF
        IF( IWQKA(IZ)  ==  5 )THEN
          UMRM = MAX(U(L,K), U(LEC(L),K))
          VMRM = MAX(V(L,K), V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          WQVREA = 3.7*XMRM
        ENDIF
  !
  ! NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS:
  !
        WQVREA = WQVREA * REAC(IZ)
        WQWREA = WQWREA * REAC(IZ)
        WQP19(L) = - (WQVREA + WQWREA) * HPKI(L,K)* WQTDKR(IWQT(L),IZ)
        WQKRDOS(L) = - WQP19(L)*WQDOS(L)
      ELSE
        WQP19(L) = 0.0
      ENDIF
    666 FORMAT(' K,IWQ,IZ,WQTDKR = ',3I5,E12.4)
  ENDDO
  !
  ! TRAPEZOIDAL SOLUTION OF KINETIC EQS: AFTER COMPUTING NEW VALUES, STORE
  ! WQVO+WQV INTO WQVO(L,K,NWQV)
  !
    DO L = 2,LA
      IZ = IWQZMAP(L,K)
      WQD6 = - WQKHR(L)
      WQKK(L) = 1.0 / (1.0 - DTWQO2*WQD6)
      WQR6 = (WQWDSL(L,K,3) + WQWPSL(L,K,3)) * VOLWQ(L)
      WQRR(L) = WQVO(L,K,3) + DTWQ*WQR6 +  DTWQO2*WQD6*WQVO(L,K,3)
      WQV(L,K,3) = SCB(L)*(WQRR(L)*WQKK(L)) + (1. - SCB(L))*WQVO(L,K,3)
      WQVO(L,K,3) = WQVO(L,K,3) + WQV(L,K,3)
    ENDDO
    DO L = 2,LA
      WQRR(L) = (WQWDSL(L,K,11)+WQWPSL(L,K,11)) * VOLWQ(L)
    ENDDO
    DO L = 2,LA
      WQKK(L) = 1.0 / (1.0 + DTWQO2*WQNIT(L))
      WQRR(L) = WQVO(L,K,11) + DTWQ*WQRR(L) &
          - DTWQO2*( WQNIT(L)*WQVO(L,K,11) )
      WQV(L,K,11) = SCB(L)*(WQRR(L)*WQKK(L)) + (1.-SCB(L))*WQVO(L,K,11)
      WQVO(L,K,11) = WQVO(L,K,11)+WQV(L,K,11)
    ENDDO
    DO L = 2,LA
      WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L))
      WQRR(L) = (WQWDSL(L,K,16) + WQWPSL(L,K,16)) * VOLWQ(L)
    ENDDO
    IF( K == KC )THEN
      DO L = 2,LA
        WQRR(L) = WQRR(L) + WQKRDOS(L)
      ENDDO
    ENDIF
    IF( K == 1 )THEN
      DO L = 2,LA
        WQRR(L) = WQRR(L) + WQBFO2(L)*HPKI(L,K)
      ENDDO
    ENDIF
    DO L = 2,LA
  !
  ! MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED
  !   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE:
  !
      WQRR(L) = WQVO(L,K,16) + DTWQ*WQRR(L) + DTWQO2*( &
          - WQAOCR*WQKHR(L)*WQVO(L,K,3) - WQAONT*WQNIT(L)*WQVO(L,K,11) &
          + WQP19(L)*WQVO(L,K,16) )
      
      WQV(L,K,16) = SCB(L)*(WQRR(L)*WQKK(L)) + (1. - SCB(L))*WQVO(L,K,16)
      WQV(L,K,16) = MAX (WQV(L,K,16), 0.0)
      WQVO(L,K,16) = WQVO(L,K,16) + WQV(L,K,16)
    ENDDO
  ENDDO
  !
  ! INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:
  !
  TIMTMP = TIMESEC/TCTMSR
  TIMESUM3 = TIMESUM3 + TIMTMP
  
  ! COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX=1/WQCHLX
  ! COUPLING TO SEDIMENT MODEL
  !: EVALUATE DEP. FLUX USING NEW VALUES CAUSE IMPLICIT SCHEME IS USED IN SPM

  ! DIURNAL DO ANALYSIS
  ! LIGHT EXTINCTION ANALYSIS

   1111 FORMAT(I12,F10.4)
   1112 FORMAT(2I5,12F7.2)
   1113 FORMAT(2I5,12E12.4)
   1414 FORMAT(I12,11E12.4)
  RETURN
  
  END SUBROUTINE WQSKE0  
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSKE1
  !
  !> @details  Modeling of unlimited algae groups
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSKE1
      
  ! *** ORGINALLY CODED BY K.-Y. PARK 
  ! *** WQSKE1 - CE-QUAL-ICM KINETICS WITH UPDATES
  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2020-04           D.K. TRAN        Adapted to unlimited phytoplankton groups
  ! 2013-03           Paul M. Craig    Added OMP
  ! 2012-08           Paul M. Craig    Added DOM component to light extinction WQKEDOM
  ! 2012-06           Paul M. Craig    Removed thickness from concentration based WQKETSS & WQKEPOC
  ! 2011-09           Paul M. Craig    Fixed Reaeration
  ! 2011-07           Paul M. Craig    Rewritten to F90
  ! 2008              SCOTT JAMES      ADDED CARBON DIOXIDE
  ! 2006-01-12        PAUL M. CRAIG    MAJOR REWRITE
 
  USE SHELLFISHMOD
  
  IMPLICIT NONE

  INTEGER NQ, NW, NS, IZ, IMWQZ, NSTPTMP, IOBC
  INTEGER ND, LF, LL, LP, L, K, IFLAG, NAL, ILAST, IPMC, i
  INTEGER, SAVE :: NWQITER                   ! *** Number of kinetic iterations
  INTEGER, SAVE :: IRPOM(3), ILPOM(3)        ! *** Refractory and Labile POM indicies
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LUSED
  
  !< In oder to adapt to the general algal group, these arrays 
  !< replace the temporary variables using for computation
  !< of the kinetic processes of water quality state variables --- DKT
  
  REAL :: WQF1N (NALGAEM)  !< Nitrogen limiation factor for current cell's biota, by class.  Silica limitation also included
  REAL :: WQF2I (NALGAEM)  !< Light limiation factor for current cell's biota, by class
  REAL :: WQTTAA(NALGAEM)
  REAL :: WQGN  (NALGAEM)  !< Concentration of available inorganic nitrogen - DKT
  REAL :: WQGP  (NALGAEM)  !< Concentration of available orthophosphate  - DKT
  REAL :: WQGCO2(NALGAEM)  !< CO2 Limitation Consts added by AA
  REAL :: WQA6A (NALGAEM)
  REAL :: WQA7A (NALGAEM)
  REAL :: WQA8A (NALGAEM)
  REAL :: WQA9A (NALGAEM)
  REAL :: WQA10A(NALGAEM)
  REAL :: WQA11A(NALGAEM)
  REAL :: WQA12A(NALGAEM)
  REAL :: WQA13A(NALGAEM)
  REAL :: WQA14A(NALGAEM)
  REAL :: WQA15A(NALGAEM)
  REAL :: WQA19A(NALGAEM)
  REAL :: WQA22A(NALGAEM)
  REAL :: RPOMF(3), LPOMF(3)

  REAL TIME, RLIGHT1, RLIGHT2, CNS1, RMULTMP
  REAL DTWQxH, DTWQxH2, TEMFAC
  REAL WQAVGIO, WQSROPT
  REAL XMRM, YMRM, WQTT1, WQKHN
  REAL WQFDI0, WQHTT, WQTTT, WQTTA
  REAL SADWQ, WQGSD, WQTTB, WQFDM
  REAL UMRM, VMRM, WQVEL, WQLVF, WQF4SC, WQKDOC, WQKHP, WQTTS
  REAL XNUMER, XDENOM, WQLDF, WQTTM
  REAL WINDREA, WQWREA, WQVREA, WQAC, WQVA1C, WQRC
  REAL WQB4, WQA4, WQR4
  REAL WQC5, WQA5, WQR5
  REAL WQD6, WQA6, WQR6
  REAL WQE7, WQA7, WQR7
  REAL WQF8, WQA8, WQR8
  REAL WQF9, WQA9, WQR9
  REAL WQR10, WQKKL
  REAL WQI11, WQA11, WQR11
  REAL WQJ12, WQA12, WQR12
  REAL WQF13, WQA13, WQR13
  REAL WQR14, WQF14, WQA14
  REAL WQR15, WQA15, WQB15
  REAL WQM16, WQA16D, WQR16, WQR17, WQR18
  REAL WQA19, WQSUM, WQRea, WQPOC, WQDOC, WQNH3, WQCOD, WQRes
  REAL WQT20, WQR21, TIMTMP, WQTAMD
  REAL PPCDO, TMP22, WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC
  REAL WQCDREA, WQCDSUM, ALGCOUNT
  REAL WQKESS, EXPA0, EXPA1, WQISM                         ! VARIABLES FOR LIGHT EXTINCTION
  REAL PSMLMULTIPLIER
  
  REAL(RKD) :: TVAL1
  REAL(RKD), STATIC :: AVGNEXT, AVGLAST
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: WQIS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQISC
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQISD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQISG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQIBOT
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQI0TOP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: WQO
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:) :: WQOLD
  
  ! ***  1) ROC - Refractory particulate organic carbon
  ! ***  2) LOC - Labile particulate organic carbon
  ! ***  3) DOC - Dissolved organic carbon
  ! ***  4) ROP - Refractory particulate organic phosphorus
  ! ***  5) LOP - Labile particulate organic phosphorus
  ! ***  6) DOP - Dissolved organic phosphorus
  ! ***  7) P4D - Total phosphate
  ! ***  8) RON - Refractory particulate organic nitrogen
  ! ***  9) LON - Labile particulate organic nitrogen
  ! *** 10) DON - Dissolved organic nitrogen
  ! *** 11) NHX - Ammonia nitrogen
  ! *** 12) NOX - Nitrate nitrogen
  ! *** 13) SUU - Particulate biogenic silica
  ! *** 14) SAA - Dissolved available silica
  ! *** 15) COD - Chemical oxygen demand
  ! *** 16) DOX - Dissolved oxygen
  ! *** 17) TAM - Total active metal
  ! *** 18) FCB - Fecal coliform bacteria
  ! *** 19) CO2 - Dissolved carbon dioxide
  ! *** 20) BIO - Any number of biota classes (NALGAE), phytoplankton, macrophytes and zooplankton
  ! ***           Can include green, blue-green, cyanobacteria, etc.  
  IF(.NOT. ALLOCATED(WQIS))THEN
    ALLOCATE(LUSED(LCM))                 ! *** Flag to indicate whether the L index has already been used for dry BC conditions
    
    ALLOCATE(WQIS(LCM,NALGAEM))
    ALLOCATE(WQISC(LCM))  
    ALLOCATE(WQISD(LCM))  
    ALLOCATE(WQISG(LCM))  
    
    ALLOCATE(WQIBOT(LCM))                ! *** Solar Radiation at the bottom  of the current layer, accounting for shade and converted PAR
    ALLOCATE(WQI0TOP(LCM))               ! *** Solar Radiation at the surface of the current layer, accounting for shade and converted PAR
    ALLOCATE(WQO(LCM,NWQVM))
    ALLOCATE(WQOLD(NBCSOP,KCM,0:NWQVM))  ! *** ALLOCATE MEMORY FOR VARIABLE TO STORE CONCENTRATIONS AT OPEN BOUNDARIES
    
    WQIS = 0.0
    WQISC = 0.0
    WQISD = 0.0
    WQISG = 0.0       
    WQIBOT = 0.0
    WQI0TOP = 0.0
    WQO = 0.0
    WQOLD = 0.0
    AVGNEXT = TBEGIN + WQHRAVG
    AVGLAST = TBEGIN
    NWQITER = 0
    
    LUSED = 0
    IRPOM(1) = IROC
    IRPOM(2) = IROP
    IRPOM(3) = IRON
    ILPOM(1) = ILOC
    ILPOM(2) = ILOP
    ILPOM(3) = ILON 
    IPMC = 0
    
  ENDIF
  NWQITER = NWQITER + 1
  PSMLMULTIPLIER = 1.0
  IF( ISTL == 3 )THEN
    PSMLMULTIPLIER = (DTWQ*86400. + DT/FLOAT(NTSTBC))/DT   ! 2.0
  ENDIF
  
  ! *** Compute the optimal average light intensity over the last three days
  IF( IWQSUN  == 2 )THEN  
    WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2
  ELSE
    WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3
  ENDIF 
  
  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO NQ = 0,NWQVM
    DO K = 1,KC
      DO IOBC = 1,NBCSOP  
        L = LOBCS(IOBC)
        WQOLD(IOBC,K,NQ) = WQV(L,K,NQ)
      ENDDO
    ENDDO  
  ENDDO  
  
  ! *** ZERO RATES
  DO L = 2,LA
    ! *** DRY CELL BYPASS
    IF( .NOT. LMASKDRY(L) )THEN
      DO NAL = 1,NALGAE
        WQPA(L,NAL) = 0.
        WQBM(L,NAL) = 0.
        WQPR(L,NAL) = 0.
        
        ! *** Apply "death rate" to fixed biota if cell is dry
        IF( .NOT. ALGAES(NAL).ISMOBILE )THEN
          DO K = KSZ(L),KC
            IF( K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L) )THEN
              WQAC = -ALGAES(NAL).WQPRRA(IWQZMAP(L,K))*DTWQO2
              WQVA1C = 1.0 / (1.0 - WQAC)
              WQV(L,K,19+NAL) = (WQV(L,K,19+NAL) + WQAC*WQV(L,K,19+NAL))*WQVA1C
              WQV(L,K,19+NAL) = MAX(WQV(L,K,19+NAL), ALGAES(NAL).WQBMIN(IWQZMAP(L,K)))
              WQO(L,19+NAL) = WQVO(L,K,19+NAL) + WQV(L,K,19+NAL)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  
  ! *** Set vegetative growth and drag 
  DO NAL = 1, NALGAE
    IF( .NOT. ALGAES(NAL).ISMOBILE )THEN
      !IF( ALGAES(NAL).ISDRAG > 0 )THEN
      IF( ALGAES(NAL).THRESHOLD /= 0 )THEN
        Call Macro_Veg(NAL)
      ENDIF
    ENDIF
  ENDDO
  
  CNS1 = 2.718  
  NS = 1  
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, ILAST, K, NQ, NW, IZ, IMWQZ, NAL)  &
  !$OMP             PRIVATE(WQF1N, WQF2I, WQTTA, WQTTAA, WQGN, WQGP, WQGCO2, WQA6A)               &
  !$OMP             PRIVATE(WQA7A, WQA8A, WQA9A, WQA10A, WQA11A, WQA12A, WQA13A)                  &
  !$OMP             PRIVATE(WQA14A, WQA15A, WQA19A, WQA22A)                                       &
  !$OMP             PRIVATE(DTWQxH, DTWQxH2, TEMFAC, ALGCOUNT)                                    &
  !$OMP             PRIVATE(XMRM, YMRM, WQTT1, WQKHN)                                             &
  !$OMP             PRIVATE(WQFDI0, WQHTT, WQTTT)                                                 &
  !$OMP             PRIVATE(SADWQ, WQGSD, WQTTB, WQFDM)                                           &
  !$OMP             PRIVATE(UMRM, VMRM, WQVEL, WQLVF, WQF4SC, WQKDOC, WQKHP, WQTTS)               &
  !$OMP             PRIVATE(XNUMER, XDENOM, WQLDF, WQTTM)                                         &
  !$OMP             PRIVATE(WINDREA, WQWREA, WQVREA, WQAC, WQVA1C, WQRC)                          &
  !$OMP             PRIVATE(WQB4,  WQA4, WQR4)                                                    &
  !$OMP             PRIVATE(WQC5,  WQA5, WQR5)                                                    &
  !$OMP             PRIVATE(WQD6,  WQA6, WQR6)                                                    &
  !$OMP             PRIVATE(WQE7,  WQA7, WQR7)                                                    &
  !$OMP             PRIVATE(WQF8,  WQA8, WQR8)                                                    &
  !$OMP             PRIVATE(WQF9,  WQA9, WQR9)                                                    &
  !$OMP             PRIVATE(WQR10, WQKKL)                                                         &
  !$OMP             PRIVATE(WQI11, WQA11, WQR11)                                                  &
  !$OMP             PRIVATE(WQJ12, WQA12, WQR12)                                                  &
  !$OMP             PRIVATE(WQF13, WQA13, WQR13)                                                  &
  !$OMP             PRIVATE(WQR14, WQF14, WQA14)                                                  &
  !$OMP             PRIVATE(WQR15, WQA15, WQB15)                                                  &
  !$OMP             PRIVATE(WQM16, WQA16D, WQR16)                                                 &
  !$OMP             PRIVATE(WQR17, WQR18)                                                         &
  !$OMP             PRIVATE(WQA19, WQSUM, WQRea, WQPOC, WQDOC, WQNH3, WQCOD, WQRes)               &
  !$OMP             PRIVATE(WQT20, WQR21, WQTAMD)                                                 &
  !$OMP             PRIVATE(PPCDO, TMP22, WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC)                 &
  !$OMP             PRIVATE(WQCDREA, WQCDSUM)                                                     &
  !$OMP             PRIVATE(WQKESS, EXPA0, EXPA1, WQISM)                                          &
  !$OMP             PRIVATE(WQSROPT)
  DO ND = 1,NDM
    LF = (ND - 1)*LDMWET + 1  
    LL = MIN(LF+LDMWET-1,LAWET)
  
    ! COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX=1/WQCHLX  

    ! ***************************************************************************
    ! *** INITIALIZE SOLAR RADIATION AND OPTIMAL LIGHT
    IF( LDAYLIGHT .AND. IWQSUN == 2 )THEN
      ! *** INITIAL SOLAR RADIATION AT TOP OF SURFACE LAYER (SHADING AND ICECOVER ALREADY ACCOUNTED FOR)
      DO LP = LF,LL
        L = LWET(LP)
        WQI0TOP(L) = PARADJ*2.065*RADTOP(L,KC)   ! *** Solar radiation in Langleys/day, PAR adjusted
      ENDDO
    ELSE
      ! *** Initial solar radiation at top of surface layer (shading and icecover already accounted for)
      DO LP = LF,LL
        L = LWET(LP)
        WQI0TOP(L) = WQI0                        ! *** Solar radiation in Langleys/day 
      ENDDO
    ENDIF

    ! ***************************************************************************
    ! *** REDISTRIBUTION OF SHELLFISH TO EACH LAYER *** shellfish 1
    IF( ISFFARM > 0 .AND. NSF > 0 )THEN
      DO L = LF,LL
        CALL SHELLFISH_REDIST(L)
      ENDDO
    ENDIF

    ! ***************************************************************************
    ! *** DZWQ=1/H (for a layer), VOLWQ=1/VOL (m^-3)
    DO K = KC,1,-1  
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        ! *** HPK(L,K)                    ! *** Layer thickness of a cell in meters
        ! *** HPKI(L,K)                   ! *** Inverse layer thickness
        TWQ(L)    = MAX(TEM(L,K), 0.0)    ! *** Layer temperature for WQ calcs (ice formation allows for small negative temperatures)
        SWQ(L)    = MAX(SAL(L,K), 0.0)    ! *** Layer salinity for WQ calcs              (KG/M3)
        VOLWQ(L)  = HPKI(L,K)*DXYIP(L)    ! *** Inverse volume of each cell in a layer   (1/M^3)
        IMWQZT(L) = IWQZMAP(L,K)          ! *** WQ Zone Map for current layer
      ENDDO  
  
      ! *** ZERO WQWPSL IF FLOWS ARE NEGATIVE.  THESE ARE HANDLED IN CALFQC (PMC)
      IF( IWQPSL /= 2 )THEN
        DO NQ = 1,NQSIJ  
          IF( (QSERCELL(K,NQ) + QSS(K,NQ)) < 0.0 )THEN
            ! *** ZERO THE FLUX
            L = LQS(NQ)  
            DO NW = 1,NWQV
              WQWPSL(L,K,NW) = 0.0
            ENDDO
          ENDIF
        ENDDO
      ENDIF
        
      ! *** DETERMINE THE RATE OF ALGAE LEAVING THE CELL THROUGH SETTLING OR FLOATING
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
    
        DO NAL = 1,NALGAE   !< Loop for all alage groups  - DKT
          IF( ALGAES(NAL).ISMOBILE )THEN
            ! *** Phytoplankton
            IF( ALGAES(NAL).WQWS(IMWQZT(L)) < 0  )THEN ! *** PERMITS ALGAES TO FLOAT AND/OR SETTLE
              IF( K == KC )THEN      
                WQBSET(L,1,NAL) = 0.0                         ! *** ALGAES AT THE WATER SURFACE CAN'T LEAVE CELL
              ELSE
                WQBSET(L,1,NAL) = -ALGAES(NAL).WQWS(IMWQZT(L))*HPKI(L,K)   ! *** CYANOBACTERIA NEEDS TO BE A POSITIVE QTY  
              ENDIF
            ELSE
              WQBSET(L,1,NAL) = ALGAES(NAL).WQWS(IMWQZT(L))*HPKI(L,K)   
            ENDIF
          ENDIF
        ENDDO           

        ! *** ZONE SPECIFIC SETTING VELOCITIES for POM, (m/day)   
        WQRPSET(L,1) = WQWSRP(IMWQZT(L))*HPKI(L,K)  ! *** Refractory POM  (1/day)
        WQLPSET(L,1) = WQWSLP(IMWQZT(L))*HPKI(L,K)  ! *** Labile POM      (1/day)
      ENDDO

      ! *** SET SETTLING FOR TAM SORPTION: CURRENT LAYER  
      IF( IWQSRP == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQWSSET(L,1) = WQWSS(IMWQZT(L))*HPKI(L,K)  
        ENDDO  
      ENDIF  
 
      ! *** SET CURRENT LAYER WQ ZONE FOR LAYERS ABOVE AND BELOW                                Note - IMWQZT1 and IMWQZT2 should be a 2D array so can set once and use every iteration...
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IF( K /= KC )THEN 
          IMWQZT1(L) = IWQZMAP(L,K+1)   ! *** IMWQZT1 - LAYER ABOVE CURRENT LAYER
        ENDIF
        IF( K /= KSZ(L) )THEN
          IMWQZT2(L) = IWQZMAP(L,K-1)   ! *** IMWQZT2 - LAYER BELOW CURRENT LAYER
        ENDIF
      ENDDO 

      ! *** COMPUTE THE MATERIAL COMING INTO THE CURRENT LAYER FROM LAYER ABOVE
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        DO NAL = 1,NALGAE   !< Loop for all alage groups  - DKT
          ! *** Phytoplankton
          ! *** All layers except the top layer
          IF( K /= KC .AND. ALGAES(NAL).ISMOBILE )THEN  
            IF( ALGAES(NAL).WQWS(IMWQZT1(L)) < 0 )THEN
              WQBSET(L,2,NAL) = 0.0  
            ELSE
              WQBSET(L,2,NAL) = ALGAES(NAL).WQWS(IMWQZT1(L))*HPKI(L,K)  
            ENDIF
          ENDIF
          ! *** ALL LAYERS EXCEPT THE BOTTOM LAYER
          IF( K /= KSZ(L) .AND. ALGAES(NAL).ISMOBILE )THEN
              IF( ALGAES(NAL).WQWS(IMWQZT2(L)) < 0 )THEN
                WQBSET(L,2,NAL) = WQBSET(L,2,NAL) - ALGAES(NAL).WQWS(IMWQZT2(L))*HPKI(L,K)  
              ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF( K /= KC )THEN
        ! *** Flux of particulate into the layer from the layer above
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRPSET(L,2) = WQWSRP(IMWQZT1(L))*HPKI(L,K)    ! *** (1/day)
          WQLPSET(L,2) = WQWSLP(IMWQZT1(L))*HPKI(L,K)    ! *** (1/day)  
        ENDDO  
      ENDIF

      ! *** Set settling for tam sorption: One layer up
      IF( IWQSRP == 1 .AND. K /= KC )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQWSSET(L,2) = WQWSS(IMWQZT1(L))*HPKI(L,K)  
        ENDDO  
      ENDIF

      ! *** FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY  
      DO LP = 1,LLWET(K,ND) 
        L = LKWET(LP,K,ND)  
        IWQT(L) = NINT((TWQ(L) - WQTDMIN)/WQTDINC) + 1  
        
        IF( IWQT(L) < 1 .OR. IWQT(L) > NWQTD )THEN
          WRITE(6,600) Map2Global(L).IG, Map2Global(L).JG, K, TWQ(L), HP(L)
          if( process_id == master_id )THEN
            OPEN(1,FILE=OUTDIR//'ERROR.LOG',POSITION='APPEND',STATUS='UNKNOWN')
            WRITE(1,*)' *** ERROR IN WQSKE1:TEMPERATURE LOOKUP TABLE'
            WRITE(1,911) TIMEDAY, Map2Global(L).LG, Map2Global(L).IG, Map2Global(L).JG, K, TWQ(L),TEM1(L,K), HP(L), H1P(L)
            WRITE(1,'(A)')'SURROUNDING DEPTHS'
            WRITE(1,'(2X,A14,I5,4F14.4)')'HP  ',L,HP(LWC(L)), HP(LEC(L)), HP(LSC(L)), HP(LNC(L))
            WRITE(1,'(2X,A14,I5,4F14.4)')'H1P ',L,H1P(LWC(L)),H1P(LEC(L)),H1P(LSC(L)),H1P(LNC(L))
            WRITE(1,'(A)')'FLUX TERMS'
            WRITE(1,'(2X,A14,I5,4E14.6)')'UHDYE/VHDXE'  ,L,UHDYE(L), UHDYE(LEC(L)), VHDXE(L), VHDXE(LNC(L))
            WRITE(1,'(2X,A14,I5,4E14.6)')'UHDY1E/VHDX1E',L,UHDY1E(L),UHDY1E(LEC(L)),VHDX1E(L),VHDX1E(LNC(L))
            CLOSE(1,STATUS='KEEP')
          end if

          IWQT(L) = MAX(IWQT(L),1)
          IWQT(L) = MIN(IWQT(L),NWQTD)

        ENDIF  
      ENDDO  

      600 FORMAT(' WQ TEM LOOKUP TABLE ERROR:  I,J,K,TEM,HP = ',3I5,4E12.4)  
      911 FORMAT('ERROR: TIME, L, I, J, K, TWQ, TEM, HP, H1P = ',F10.5, I7, 3I4, 4E12.4,/)  

      ! *** BEGIN HORIZONTAL LOOP FOR NUTRIENT AND LIGHT LIMITATION OF ALGAL GROWTH
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IZ = IWQZMAP(L,K)
        
        RNH4WQ(L) = MAX (WQV(L,K,INHX), 0.0)                                          ! *** Ammonia
        RNO3WQ(L) = MAX (WQV(L,K,INOX), 0.0)                                          ! *** Nitrate
        PO4DWQ(L) = MAX (WQPO4D(L,K), 0.0)                                            ! *** Phosphate
        RNH4NO3(L) = RNH4WQ(L) + RNO3WQ(L)                                            ! *** Total Inorganic Nitrogen
          
        ! *** Loop over biota for N/P/Si limiations
        DO NAL = 1,NALGAE
          IF( ALGAES(NAL).ISMOBILE  )THEN
            ! *** Phytoplankton
            WQGN(NAL) = RNH4NO3(L) / (ALGAES(NAL).WQKHNA + RNH4NO3(L)+ 1.E-18)
            WQGP(NAL) = PO4DWQ(L)  / (ALGAES(NAL).WQKHPA + PO4DWQ(L) + 1.E-18)
            WQF1N(NAL) = MIN(WQGN(NAL), WQGP(NAL))                                  ! *** Minimum of the N/P 
            IF( ISKINETICS(ICO2) > 0  )THEN
              CO2WQ(L)  = MAX (WQV(L,K,ICO2), 0.0)                                  ! *** CO2 
              WQGCO2(NAL) = CO2WQ(L) / (ALGAES(NAL).WQKHCO2 + CO2WQ(L) + 1.E-18) 
              WQF1N(NAL) = MIN(WQGN(NAL), WQGP(NAL), WQGCO2(NAL))                   ! *** Minimum of the N/P/CO2     
            ENDIF
          ENDIF
            
          IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
            ! *** Calculate nutrient limitations for Macrophytes and periphyton
            WQGN(NAL) = RNH4NO3(L) / (ALGAES(NAL).WQKHNA + RNH4NO3(L) + 1.E-18)
            WQGP(NAL) = PO4DWQ(L)  / (ALGAES(NAL).WQKHPA + PO4DWQ(L)  + 1.E-18)
            WQF1N(NAL) = MIN(WQGN(NAL), WQGP(NAL))                                  ! ***  Minimum of the N/P
            IF( ISKINETICS(ICO2) > 0  )THEN
              WQGCO2(NAL) = CO2WQ(L) / (ALGAES(NAL).WQKHCO2 + CO2WQ(L) + 1.E-18)
              WQF1N(NAL) = MIN(WQGN(NAL), WQGP(NAL), WQGCO2(NAL))                   ! ***  Minimum of the N/P/CO2
            ENDIF
          ENDIF
            
          IF( IWQSI == 1 )THEN  
            ! *** Silica limitation
            SADWQ = MAX (WQSAD(L,K), 0.0)
            WQGSD = SADWQ / (ALGAES(NAL).WQKHS + SADWQ + 1.E-18)
            IF( ALGAES(NAL).ISILICA /= 0 )THEN
              WQF1N(NAL) = MIN(WQF1N(NAL),WQGSD)                                    ! *** Limit: Diatoms - Minimum of the N/P/S
            ENDIF
          ENDIF                
        ENDDO   ! *** End of biota loop
                    
        ! *** LIGHT EXTINCTION (THIS WILL ALWAYS BE TRUE EXCEPT FOR IWQSUN=2)
        IF( WQI0 > 0.1 )THEN
          ! *** GET THE EXTINCTION COEFFICIENT
          WQKESS = RADKE(L,K)
        
          ! *** OPTIMAL LIGHT INTENSITY AT OPTIMAL DEPTH
          IF( K == KC )THEN
            ! *** Only compute optimal light once, even if optimal depth is below layer KC
            DO NAL = 1,NALGAE
              IF( ALGAES(NAL).ISMOBILE )THEN
                ! *** Phytoplankton
                WQIS(L,NAL) = MAX( WQAVGIO*EXP(-WQKESS*ALGAES(NAL).WQDOP(IZ)), WQISMIN )  
                
                ! *** HARDWIRE TO SET OPTIMAL GROWTH TO A FIXED VALUE USED BY CHAPRA (Algal Growth Example (Chapra 33.2 with Fixed IS250).
                !WQIS(L,NAL) = 250./WQFD         ! *** Is
              ENDIF
            ENDDO
          ENDIF

          ! *** Current light growth limiting factor. 
          IF( K == KC )THEN    
            WQITOP(L,K) = WQI0TOP(L)                            ! *** WQITOP is solar radiation at the TOP of layer K
            WQIBOT(L)   = WQI0TOP(L)*EXP(-WQKESS*HPK(L,K))      ! *** WQIBOT is solar radiation at the BOTTOM of layer K
          ELSE            
            WQITOP(L,K) = WQIBOT(L)
            WQIBOT(L)   = WQITOP(L,K)*EXP(-WQKESS*HPK(L,K)) 
          ENDIF !SEE DiTORO ET AL (1971, EQNS. (11)&(12)) 

          WQTT1 = WQFD*2.718                                     ! *** EXP(1.0) = 2.718
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE )THEN
              ! *** Phytoplankton
              EXPA0 = EXP(-WQITOP(L,K)/(WQFD*WQIS(L,NAL)))              ! *** WQITOP(L,K)/WQIS(L,NAL) is the fraction of actual SR to average SR at the TOP 
              EXPA1 = EXP(-WQIBOT(L)  /(WQFD*WQIS(L,NAL)))              ! *** WQITOP(L,K)/WQIS(L,NAL) is the fraction of actual SR to average SR at the BOTTOM
              WQF2I(NAL) = WQTT1/(HPK(L,K)*WQKESS)*(EXPA1 - EXPA0)
            ELSE
              ! *** Macrophytes and periphyton - Light Limitation at top of growth layer
              WQF2I(NAL) = 0.0
              IF( WQITOP(L,K) > 1.0E-18 )THEN  
                WQFDI0 = - WQIBOT(L)/(WQFD + 1.E-18)
                WQISM  = MAX(WQAVGIO*EXP(-WQKESS*ALGAES(NAL).WQDOP(IZ)), WQISMIN)  ! *** Optimal Light
                WQFDM  = WQFDI0/(WQISM + 1.E-18)                                   ! *** Ratio of Actual to Optimal Light
                WQHTT  = WQHT(K) * HP(L)  
                WQTTB  = EXP(-WQKESS * (WQHTT + 1.0/HPKI(L,K)))  
                WQTTT  = EXP(-WQKESS * WQHTT)  
                WQTT1  = (CNS1 * WQFD * HPKI(L,K))/WQKESS  
                WQF2I(NAL) = WQTT1 * (EXP(WQFDM*WQTTB) - EXP(WQFDM*WQTTT))         ! *** Light Based Macroalgae Growth Limiting Factor   
                WQIS(L,NAL) = WQISM
              ENDIF
            ENDIF
          ENDDO
        ELSE
          ! *** No Light Case
          WQIBOT(L) = 0.
          WQITOP(L,K) = 0.
          WQKESS = 0.
          DO NAL = 1,NALGAE
            WQF2I(NAL) = 0.0
          ENDDO
        ENDIF
            
        DO NAL = 1,NALGAE
          IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
            ! *** Velocity Limitation for macrophytes and periphyton
            WQLVF = 1.0
            IF( ALGAES(NAL).IWQVLIM > 0 )THEN
              UMRM  = 0.5*(U(L,K) + U(LEC(L),K))  
              VMRM  = 0.5*(V(L,K) + V(LNC(L),K))  
              WQVEL = SQRT(UMRM*UMRM + VMRM*VMRM)  
              
              ! *** OPTION 1 FOR VELOCITY LIMITATION ASSUMES MACROALGAE GROWTH  
              ! *** IS LIMITED AT LOW VELOCITIES DUE TO REDUCED AVAILABILITY OF  
              ! *** NUTRIENTS REACHING THE ALGAE BIOMASS.    
              ! *** USES A MICHAELIS-MENTON or MONOD TYPE OF EQUATION.  
              IF( ALGAES(NAL).IWQVLIM  ==  1 )THEN  
                IF( WQVEL  > ALGAES(NAL).WQKMVMIN(IZ) )THEN  
                  WQLVF = WQVEL / (ALGAES(NAL).WQKMV(IZ) + WQVEL)  
                ELSE  
                  WQLVF = ALGAES(NAL).WQKMVMIN(IZ) / (ALGAES(NAL).WQKMV(IZ) + ALGAES(NAL).WQKMVMIN(IZ))  
                ENDIF  
              ENDIF       
              
              ! *** OPTION 2 FOR VELOCITY LIMITATION APPLIES A FIVE-PARAMETER LOGISTIC  
              ! *** FUNCTION THAT CAN BE ADJUSTED TO LIMIT MACROALGAE GROWTH FOR  
              ! *** EITHER LOW OR HIGH (SCOUR) VELOCITIES.  IN STREAMS WITH LOW NUTRIENTS,  
              ! *** THE LOW VELOCITY WILL LIKELY BE LIMITING SINCE AMPLE NUTRIENTS MAY  
              ! *** NOT REACH THE ALGAE BIOMASS DUE TO REDUCED FLOW.  IN STREAMS WITH  
              ! *** ABUNDANT NUTRIENTS, LOW VELOCITIES WILL NOT LIMIT MACROALGAE GROWTH,  
              ! *** INSTEAD, HIGH VELOCITIES WILL LIKELY SCOUR THE MACROALGAE AND DETACH  
              ! *** IT FROM THE SUBSTRATE.  
              IF( ALGAES(NAL).IWQVLIM  == 2 )THEN  
                XNUMER = ALGAES(NAL).WQKMVA(IZ) - ALGAES(NAL).WQKMVD(IZ)  
                XDENOM = 1.0 + (WQVEL/ALGAES(NAL).WQKMVC(IZ))**ALGAES(NAL).WQKMVB(IZ)  
                WQLVF = ALGAES(NAL).WQKMVD(IZ) + ( XNUMER / (XDENOM**ALGAES(NAL).WQKMVE(IZ)) )  
              ENDIF  
            ENDIF
              
            ! *** USE THE MORE SEVERELY LIMITING OF VELOCITY OR NUTRIENT FACTORS:  
            WQF1N(NAL) = MIN(WQLVF, WQF1N(NAL)) 

            ! *** Crowding limitation factor based on WQKBP optimal plant density
            XMRM = WQV(L,K,19+NAL)*HPK(L,K)          ! *** Convert WQV (gC/M3) to a density: XMRM (gC/M2)  
            WQLDF = ALGAES(NAL).WQKBP(IZ) / (ALGAES(NAL).WQKBP(IZ) + XMRM)
              
            ! ***                 Max growth rate       Nutr/Vel    Light     Temperature       Crowding
            WQPA(L,NAL) = ALGAES(NAL).WQPMA(IMWQZT(L))*WQF1N(NAL)*WQF2I(NAL)*WQTDG(IWQT(L),NAL)*WQLDF   ! *** Macroalgae growth rate
              
            WQBM(L,NAL) = ALGAES(NAL).WQBMRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted macroalgae metabolism rate
            WQPR(L,NAL) = ALGAES(NAL).WQPRRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted macroalgae predation rate 
            
          ENDIF
        ENDDO
          
        ! *** Compute the Growth Rate based on Maximums & Limiting Factors             
        DO NAL = 1,NALGAE
          IF( ALGAES(NAL).ISMOBILE )THEN
            ! *** Phytoplankton
            WQPA(L,NAL) = ALGAES(NAL).WQPMA (IMWQZT(L)) * WQF1N(NAL)*WQF2I(NAL)*WQTDG(IWQT(L),NAL)      ! *** Nutrient & Temperature Adjusted
            WQBM(L,NAL) = ALGAES(NAL).WQBMRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted
            WQPR(L,NAL) = ALGAES(NAL).WQPRRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted
          ENDIF
            
          IF( ALGAES(NAL).ISTOX /= 0 )THEN 
            WQF4SC  = ALGAES(NAL).WQSTOX / (ALGAES(NAL).WQSTOX + SWQ(L)*SWQ(L) + 1.E-12)
            WQPA(L,NAL) = WQPA(L,NAL) * WQF4SC
          ENDIF
            
          ! THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A  
          ! LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO  
          ! BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.  
          IF( ALGAES(NAL).ISBLOOM /= 0 )THEN
            WQBM(L,NAL) = WQBM(L,NAL)*WQTDP(IWQT(L),NAL)                                                ! *** Temperature Adjusted for winter bloom
            WQPR(L,NAL) = WQPR(L,NAL)*WQTDP(IWQT(L),NAL)                                                ! *** Temperature Adjusted for winter bloom
          ENDIF
            
          ! The term WQPR is replaced by a constant algae death when zooplankton is actived.
          ! Also the fraction of nutrient produced by algae predation are still used for death process,
          ! So that the solving process for kinetic equation of algae does not change.
          IF( IWQZPL > 0 .AND. ALGAES(NAL).ISMOBILE )THEN
            WQPR(L,NAL) = ALGAES(NAL).WQDRA(IZ)
          ENDIF
        ENDDO
      ENDDO      ! *** END ACTIVE CELL LOOP FOR ALGAE PARAMETERS  
        
      ! *** SET UP OTHER NUTRIENT VARIABLES FOR MAIN KINETICS SECTION
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IZ = IWQZMAP(L,K)
        WQOBTOT(L) = 0.0
        XMRM = 0.0
        
        DO NAL = 1,NALGAE
          IF( ALGAES(NAL).ISMOBILE )THEN           ! *** Phytoplankton
            WQOBTOT(L) = WQOBTOT(L) + WQV(L,K,19+NAL)           
          ELSEIF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN    ! *** Fixed biota
            XMRM = XMRM + ALGAES(NAL).WQKDCALM(IZ) * WQV(L,K,19+NAL)  
          ENDIF
        ENDDO
        
        WQKRPC(L) = (WQKRC + WQKRCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))        ! *** Hydrolysis:     RPOP--> DOC
        WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))        ! *** Hydrolysis:     RPOP--> DOC
        
        ! *** M. MORTON added spatially variable doc minerization rate WQKDC  
        ! *** to achieve better control in systems with a combination of freshwater
        ! *** streams and tidal rivers with different characteristics.  

        ! ***   Min Mineral     Factor   MobC   Factor*MacC   T Correct
        WQKDOC = (WQKDC(IZ) + WQKDCALG*WQOBTOT(L) + XMRM)*WQTDMNL(IWQT(L))     ! *** Heterotrophic respiration rate of DOC at infinite DO (1/day)
        
        O2WQ(L) = MAX(WQV(L,K,IDOX), 0.0)  
        WQTT1 = WQKDOC / (WQKHORDO + O2WQ(L) + 1.E-18)  
        WQKHR(L)   = WQTT1*O2WQ(L)  
        WQDENIT(L) = WQTT1*WQAANOX*RNO3WQ(L)/(WQKHDNN + RNO3WQ(L) + 1.E-18)    ! *** Denitrification Rate  
      ENDDO  

      ! ***********************************
      ! 7-10 PHOSPHORUS  
      ! *** HYDROLYSIS
      
      WQKHP = 0.0
      ALGCOUNT = 0.
      DO NAL = 1,NALGAE
        IF( ALGAES(NAL).ISMOBILE )THEN
          WQKHP = WQKHP + ALGAES(NAL).WQKHPA                            ! *** Mean phosphorus half-saturation for algae
          ALGCOUNT =  ALGCOUNT + 1.
        ENDIF
      ENDDO
      IF( ALGCOUNT > 0. ) WQKHP = WQKHP/ALGCOUNT   ! *** Note - better if mass weighted
      
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        WQAPC(L) = 1.0/(WQCP1PRM + WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ(L)))   ! *** A mean Phosphorus-to-Carbon ration for all algal groups
        WQTT1 = WQKHP/(WQKHP + PO4DWQ(L) + 1.E-18) * WQOBTOT(L)  
        WQKRPP(L) = (WQKRP + WQKRPALG*WQTT1) * WQTDHDR(IWQT(L))         ! *** Hydrolysis:     RPOP--> DOP
        WQKLPP(L) = (WQKLP + WQKLPALG*WQTT1) * WQTDHDR(IWQT(L))         ! *** Hydrolysis:     LPOP--> DOP
        WQKDOP(L) = (WQKDP + WQKDPALG*WQTT1) * WQTDMNL(IWQT(L))         ! *** Mineralization: DOP --> PO4
      ENDDO
  
      ! *** PHOSPHATE SETTLING   Note - THIS COULD BE SPED UP BY CHECKING OPTIONS OUTSIDE OF THE L LOOP
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IF( IWQSRP == 1 )THEN  
          WQTTM = WQKPO4P*WQTAMP(L,K)                                   ! *** Sorbed mass fraction to TAM
          WQH10(L) = - WQWSSET(L,1) * WQTTM / (1.0 + WQTTM)             ! *** Loss from the layer out bottom
          IF( K /= KC )THEN  
            WQTTM = WQKPO4P*WQTAMP(L,K+1)  
            WQT10(L) = WQWSSET(L,2) * WQTTM / (1.0 + WQTTM)             ! *** Gain to the layer from above
          ENDIF  
        ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN  
          WQTTS = WQKPO4P*SEDT(L,K)                                     ! *** Sorbed mass fraction to cohesive sediments  
          WQH10(L) = - WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)    ! *** Loss from the layer 
          IF( K /= KC )THEN  
            WQTTS = WQKPO4P*SEDT(L,K+1)  
            WQT10(L) = WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)    ! *** Gain to the layer from above  
          ENDIF  
        ELSE  
          WQH10(L) = 0.0  
          WQT10(L) = 0.0  
        ENDIF 
        WQH10(L) = WQH10(L)*DTWQO2 
      ENDDO  

      ! ***********************************
      ! 11-15 NITROGEN  
      ! *** HYDROLYSIS
      ! *** Mean nitrogen half-saturation for algae
      WQKHN = 0.0
      ALGCOUNT = 0.
      DO NAL = 1,NALGAE
        IF( ALGAES(NAL).ISMOBILE )THEN
          WQKHN = WQKHN + ALGAES(NAL).WQKHNA
          ALGCOUNT =  ALGCOUNT + 1.
        ENDIF
      ENDDO
      IF( ALGCOUNT > 0. ) WQKHN = WQKHN/ALGCOUNT   ! *** Note - better if mass weighted
      
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        WQTT1 = WQKHN / (WQKHN + RNH4NO3(L)+ 1.E-18) * WQOBTOT(L)  
        WQKRPN(L) = (WQKRN + WQKRNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** Hydrolysis:     RPON-->DON
        WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** Hydrolysis:     LPON-->DON
        WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** Mineralization: DON -->NH4
      ENDDO
      
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        
        ! *** Ammonium Preference
        DO NAL = 1, NALGAE
          WQTTA = RNH4WQ(L)/(ALGAES(NAL).WQKHNA + RNO3WQ(L) + 1.E-18)
          WQPN(L,NAL) = ( RNO3WQ(L)/(ALGAES(NAL).WQKHNA + RNH4WQ(L) + 1.E-18) + ALGAES(NAL).WQKHNA/(RNH4NO3(L) + 1.E-18) ) * WQTTA   ! *** Ammonium preference
        ENDDO
        WQNIT(L) = WQTDNIT(IWQT(L)) * O2WQ(L) / (WQKHNDO + O2WQ(L) + 1.E-18) * RNH4WQ(L) / (WQKHNN + RNH4WQ(L) + 1.E-18)
      ENDDO
      
      IF( IWQSI == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IF( IWQSRP == 1 )THEN  
            WQTTM = WQKSAP*WQTAMP(L,K)  
            WQN17(L) = - WQWSSET(L,1) * WQTTM / (1.0 + WQTTM)  
            IF( K /= KC )THEN  
              WQTTM = WQKSAP*WQTAMP(L,K+1)  
              WQT17(L) = WQWSSET(L,2) * WQTTM / (1.0 + WQTTM)  
            ENDIF  
          ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN  
            WQTTS = WQKSAP*SEDT(L,K)  
            WQN17(L) = - WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)  
            IF( K /= KC )THEN  
              WQTTS = WQKSAP*SEDT(L,K+1)  
              WQT17(L) = WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)  
            ENDIF  
          ELSE  
            WQN17(L) = 0.0  
            WQT17(L) = 0.0  
          ENDIF  
        ENDDO  
        WQN17(L) = WQN17(L)*DTWQO2 
      ENDIF  

      ! ***********************************
      ! *** DISSOLVED OXYGEN
      PPCDO = -3.45  !PARTIAL PRES OF CO2 IN 10^ppcdo ATM; TEMPORARILY DECLARED HERE. SHLD BE READ IN FROM INPUT FILE
      DO LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IZ = IWQZMAP(L,K)  
        WQO18(L) = -DTWQO2*WQKCOD(IWQT(L),IZ)*O2WQ(L)/(WQKHCOD(IZ) + O2WQ(L) + 1.E-18)  

        WQDOS(L) = DO_SAT(L)

        XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*HPK(L,K)
  
        !************* CO2 parameters
        !VB COMPUTING THE pK FOR SAT CONC OF CO2; K - HENRY'S CONST
        CDOSATIDX(L) = -2385.73/(TWQ(L) + 273.15) -  0.0152642 * (TWQ(L) + 273.15) + 14.0184
        !          K * MOL WT OF CO2 * PARTAL PRES OF CO2 IN ATM
        WQCDOS(L) = 10.**(-CDOSATIDX(L) + PPCDO) * (44.* 1000.) !VB EVALUATING CONC OF CO2 IN G/M^3 
        !************* CO2 parameters

        ! *** Compute Reaeration
        IF( K == KC )THEN       
          ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION
          WINDREA = MIN(WINDST(L),11.)  
          WQWREA = 0.728*SQRT(WINDREA) + (0.0372*WINDREA - 0.317)*WINDREA  
          WQWREA = WQWREA * HPKI(L,K)
          
          IF( IWQKA(IZ) == 0 )THEN
            ! *** Constant  
            WQVREA = WQKRO(IZ) * HPKI(L,K) 
            WQWREA = 0.0  
          ELSEIF( IWQKA(IZ) == 1 )THEN  
            ! *** Constant plus Wind
            WQVREA = WQKRO(IZ) * HPKI(L,K) 
          ELSEIF( IWQKA(IZ) == 2 )THEN
            ! *** OCONNOR-DOBBINS REAERATION FORMULA
            UMRM = 0.5*(U(L,K) + U(LEC(L),K))  
            VMRM = 0.5*(V(L,K) + V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            ! *** WQKRO = 3.933 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**0.5         ! *** Multi-layer implementation  
            WQVREA = WQVREA * HPKI(L,K)                         ! *** Multi-layer implementation
          ELSEIF( IWQKA(IZ) == 3 )THEN
            ! *** OWENS & GIBBS (1964) REAERATION FORMULA
            UMRM = MAX(U(L,K), U(LEC(L),K))  
            VMRM = MAX(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            ! *** WQKRO = 5.32 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85  
          ELSEIF( IWQKA(IZ) == 4 )THEN  
            ! *** MODIFIED OWENS AND GIBBS REAERATION EQUATION:  
            ! *** NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE  
            ! ***       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER  
            ! ***       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.  
            UMRM = MAX(U(L,K), U(LEC(L),K))  
            VMRM = MAX(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L) + 0.1524))  
            ! *** WQKRO = 5.32 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85  
          ELSEIF( IWQKA(IZ) == 5 )THEN  
            UMRM = MAX(U(L,K), U(LEC(L),K))  
            VMRM = MAX(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            WQVREA = 3.7*XMRM * HPKI(L,K)  
          ENDIF  

          ! *** NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS
          WQVREA = WQVREA * REAC(IZ)  ! *** Diffusive Flux
          WQWREA = WQWREA * REAC(IZ)  ! *** Wind Component
          WQP19(L)    = -(WQVREA + WQWREA) * WQTDKR(IWQT(L),IZ)   ! *** Reaeration rate (1/day)  
          WQKRDOS(L)  = -WQP19(L)*WQDOS(L)                        ! *** O2 Flux due to DO saturation (offset by current DO below)
          WQP22(L)    = WQP19(L)*((32./44.)**0.25)                ! *** Kr FOR CO2 ANALOGOUS TO WQP19 ; 44 = MOL WT OF CO2
          WQKRCDOS(L) = -WQP22(L) * WQCDOS(L)                     ! *** EVALUATING Kr*SAT CONC OF CO2
        ELSE
          WQKRDOS(L)  = 0.0
          WQKRCDOS(L) = 0.0  
          WQP19(L)    = 0.0  
          WQP22(L)    = 0.0              !VB Kr FOR CO2 IS ZERO FOR CELLS NOT AT THE SURFACE
        ENDIF  
      ENDDO  
  
      ! *** ICE
      IF( ISICE > 0 .AND. K == KC )THEN
        ! *** REDUCE REAERATION AND SURFACE EXCHANGE DUE TO ICE
        DO LP = LF,LL
          L = LWET(LP)
          WQP19(L)   = (1.-ICECOVER(L))*WQP19(L)
          WQKRDOS(L) = (1.-ICECOVER(L))*WQKRDOS(L)
        ENDDO
      ENDIF
    
      ! *** Trapezoidal solution of kinetic eqs: After computing new values, store WQVO + WQV into WQO

      !*******************************************************************************************************************************
      ! *** MACROPHYTES AND PERIPHYTON
      DO NAL = 1,NALGAE
        IF( .NOT. ALGAES(NAL).ISMOBILE )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            IF( K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L) )THEN
              IMWQZ = IWQZMAP(L,KSZ(L))
            
              ! *** Apply maximum allowable concentration
              IF( WQV(L,K,19+NAL) > ALGAES(NAL).WQBMAX )THEN
                WQPA(L,NAL) = 0.0                                                                    ! *** Zero the growth
              ENDIF
            
              ! ***    GROWTH        METABOLISM  PREDATION/DEATH           SETTLING(NEEDED?)
              WQAC = (WQPA(L,NAL) - WQBM(L,NAL) - WQPR(L,NAL) - ALGAES(NAL).WQWS(IMWQZ)*HPKI(L,K))*DTWQO2  

              WQVA1C = 1.0 / (1.0 - WQAC)
              WQV(L,K,19+NAL) = (WQV(L,K,19+NAL) + WQAC*WQV(L,K,19+NAL))*WQVA1C*SMAC(L)  
              WQV(L,K,19+NAL) = MAX(WQV(L,K,19+NAL), ALGAES(NAL).WQBMIN(IMWQZ))*SMAC(L)  
              WQO(L,19+NAL) = WQVO(L,K,19+NAL) + WQV(L,K,19+NAL)
            
              ! *** Hydrodynamic feedback parameters
              IF( ALGAES(NAL).ISDRAG > 0 )THEN
                ! *** Compute vegetation equivalents for drag calculations
              
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
  
      ! *** Phytoplankton       
      DO NAL = 1,NALGAE
        IF( ALGAES(NAL).ISMOBILE )THEN
          ! *** Phytoplankton
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            ! ***   GROWTH        BASAL_METAB   PREDATION     SETTLING         TIME STEP  
            WQAC = (WQPA(L,NAL) - WQBM(L,NAL) - WQPR(L,NAL) - WQBSET(L,1,NAL))*DTWQO2       ! *** Production per unit time multiplied by half time step
            WQKK(L) = 1.0 / (1.0 - WQAC) 
      
            ! ***   PT_SRC_LOADS    VOLUME  
            WQRC = WQWPSL(L,K,19+NAL) * VOLWQ(L) * PSMLMULTIPLIER            ! *** Point source load rate multiplied by inverse cell volume  g/m^3/t
            WQRR(L) = WQV(L,K,19+NAL) + DTWQ*WQRC + WQAC*WQV(L,K,19+NAL)     ! *** Transported biomass conc. (CALWQC) + point source load rate X time step + growth rate X previous biomass conc.
            
            ! *** Coupling to zooplankton - DKT
            IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) - DTWQ*SBZPAL(L,K,NAL)
          ENDDO
          
          ! *** Add in settling
          IF( K /= KC )THEN
            IF( ALGAES(NAL).WQWS(IMWQZT2(2)) > 0.0 )THEN
              ! *** From above
              DO LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)  
                WQRR(L) = WQRR(L) + DTWQO2*WQBSET(L,2,NAL)*WQVO(L,K+1,19+NAL)  ! *** Biomass conc. + DtX(1/t)* biomass conc.
              ENDDO
            ELSE
              ! *** From below
              DO LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)  
                WQRR(L) = WQRR(L) + DTWQO2*WQBSET(L,2,NAL)*WQVO(L,K-1,19+NAL)  ! *** Biomass conc. + DtX(1/t)* biomass conc.
              ENDDO
            ENDIF
          ELSE  ! K == KC
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***    ATM DRY DEP           ATM WET DEP         VOLUME  
              WQRC = (WQWDSL(L,KC,19+NAL) + WQATML(L,KC,19+NAL))*VOLWQ(L)    ! *** Atmospheric loading mass per time / cell volume
              WQRR(L) = WQRR(L) + DTWQ*WQRC                                  ! *** Biomass conc. + Dt*loading rate per unit volume
            ENDDO
          ENDIF  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQV(L,K,19+NAL) = WQRR(L)*WQKK(L)
            WQV(L,K,19+NAL) = MAX(WQV(L,K,19+NAL), ALGAES(NAL).WQBMIN(1))    ! *** Apply a user specified minimum concentration
            WQO(L,19+NAL)   = WQVO(L,K,19+NAL) + WQV(L,K,19+NAL)             ! *** Depth totaled biomass conc = old biomass conc in cell + biomass conc from this iteration
          ENDDO
        ENDIF
      ENDDO

      !*******************************************************************************************************************************
      ! *** NOW COMPUTE KINETICS FOR EACH CONSTITUENT
      
      ! ****  PARAM 01  ROC - refractory particulate organic carbon
      IF( ISKINETICS(1) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS    SETTLING  
          WQB4 = -( WQKRPC(L) + WQRPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQB4)  

          ! ***  ALGAE PREDATION SOURCE OF RPOC
          WQA4 = 0.0
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE )THEN
              ! *** Phytoplankton
              WQA4 = WQA4 + ALGAES(NAL).WQFCRP*WQPR(L,NAL)*WQO(L,19+NAL)
            ELSEIF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA4 = WQA4 + ALGAES(NAL).WQFCRP*WQPR(L,NAL)*WQO(L,19+NAL)
            ENDIF
          ENDDO 
          
          ! ***  PT_SRC_LOADS    VOLUME  
          WQR4 = WQWPSL(L,K,IROC) * VOLWQ(L) * PSMLMULTIPLIER
          
          WQRR(L) = WQV(L,K,IROC) + DTWQ*WQR4 + DTWQO2*WQA4 + WQB4*WQVO(L,K,IROC)
          
          !*** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SRPOCZ(L,K)
        ENDDO

        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,IROC)  
          ENDDO  
        ELSE  ! K == KC
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP            ATM WET DEP    VOLUME  
            WQR4 = (WQWDSL(L,KC,IROC) + WQATML(L,KC,IROC))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR4  
          ENDDO
        ENDIF  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IROC) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,IROC)   = WQVO(L,K,IROC) + WQV(L,K,IROC)
        ENDDO  
      ENDIF  

      ! ****  PARAM 02  LOC - labile particulate organic carbon
      IF( ISKINETICS(2) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS    SETTLING  
          WQC5 = - (WQKLPC(L) + WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQC5)
          
          ! *** Algae predation source                         
          WQA5 = 0.0
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE )THEN
              ! *** Phytoplankton
              WQA5 = WQA5 + ALGAES(NAL).WQFCLP*WQPR(L,NAL)*WQO(L,19+NAL)
            ELSE IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN  
              WQA5 = WQA5 + ALGAES(NAL).WQFCLP*WQPR(L,NAL)*WQO(L,19+NAL)  
            ENDIF
          ENDDO
          
          ! ***  PT_SRC_LOADS       VOLUME  
          WQR5 = WQWPSL(L,K,ILOC) * VOLWQ(L) * PSMLMULTIPLIER

          WQRR(L) = WQV(L,K,ILOC) + DTWQ*WQR5 + DTWQO2*WQA5 + WQC5*WQVO(L,K,ILOC)
          
          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SLPOCZ(L,K)
        ENDDO  
        
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,ILOC)  
          ENDDO  
        ELSE  ! K == KC
          ! *** Add surface loads
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP       ATM WET DEP        VOLUME  
            WQR5 = (WQWDSL(L,K,ILOC) + WQATML(L,KC,ILOC))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR5  
          ENDDO
        ENDIF  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,ILOC) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,ILOC)   = WQVO(L,K,ILOC) + WQV(L,K,ILOC)
        ENDDO  
      ENDIF  

      ! ****  PARAM 03  DOC - dissolved organic carbon
      IF( ISKINETICS(3) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IZ = IWQZMAP(L,K) 
          
          ! ***    RESPIRATION  DENITRIFICATION
          WQD6 = -(WQKHR(L) + WQDENIT(L)) *DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQD6)
          
          ! *** Algae metabolism and predation source 
          WQA6 = 0.0
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA6A(NAL) = ALGAES(NAL).WQFCDB + (1.0 - ALGAES(NAL).WQFCDB)*(ALGAES(NAL).WQKHRA(IZ)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18))
              WQA6 = WQA6 + ( WQA6A(NAL)*WQBM(L,NAL) + ALGAES(NAL).WQFCDP*WQPR(L,NAL) )*WQO(L,19+NAL)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA6A(NAL) = ALGAES(NAL).WQFCDB + (1.0 - ALGAES(NAL).WQFCDB)*(ALGAES(NAL).WQKHRA(IZ)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)) 
              WQA6 = WQA6 + ( WQA6A(NAL)*WQBM(L,NAL) + ALGAES(NAL).WQFCDP*WQPR(L,NAL) )*WQO(L,19+NAL) 
            ENDIF
          ENDDO

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR6 = WQWPSL(L,K,IDOC) * VOLWQ(L) * PSMLMULTIPLIER

          WQRR(L) = WQV(L,K,IDOC) + DTWQ*WQR6 + WQD6*WQVO(L,K,IDOC) + DTWQO2*(WQA6 + WQKRPC(L)*WQO(L,IROC) + WQKLPC(L)*WQO(L,ILOC))
          
          ! *** Coupling to zooplankton - DKT              
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDOCZ(L,K)
        ENDDO

        IF( K == KC )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR6 = (WQWDSL(L,K,IDOC) + WQATML(L,KC,IDOC))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR6  
          ENDDO
        ENDIF
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IDOC) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,IDOC)   = WQVO(L,K,IDOC) + WQV(L,K,IDOC)
        ENDDO  
      ENDIF  

      ! ****  PARAM 04  ROP - refractory particulate organic phosphorus
      IF( ISKINETICS(4) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQE7 = -(WQKRPP(L) + WQRPSET(L,1))*DTWQO2 
          WQKK(L) = 1.0 / (1.0 - WQE7)
          
          ! *** Algal metabolism and predation source
          WQA7 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA7A(NAL) = (ALGAES(NAL).WQFPRB*WQBM(L,NAL) + ALGAES(NAL).WQFPRP*WQPR(L,NAL)) * WQO(L,19+NAL)
              WQA7 = WQA7 + WQA7A(NAL) * WQAPC(L)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L))  )THEN
              ! *** Macrophytes and periphyton
              WQA7A(NAL) = (ALGAES(NAL).WQFPRB*WQBM(L,NAL) + ALGAES(NAL).WQFPRP*WQPR(L,NAL)) * WQO(L,19+NAL)
              WQA7 = WQA7 + WQA7A(NAL) * WQAPC(L) * ALGAES(NAL).WQAPCM 
            ENDIF
          ENDDO
          
          ! ***  PT_SRC_LOADS    VOLUME  
          WQR7 = WQWPSL(L,K,IROP) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,IROP) + DTWQ*WQR7 + DTWQO2*WQA7 + WQE7*WQVO(L,K,IROP) 
         
          ! *** Coupling to zooplankton - DKT              
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SRPOPZ(L,K)
        ENDDO  
        
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,IROP)
          ENDDO  
        ELSE
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP    ATM WET DEP    VOLUME  
            WQR7 = (WQWDSL(L,K,IROP) + WQATML(L,KC,IROP))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR7  
          ENDDO
        ENDIF  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IROP) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,IROP)   = WQVO(L,K,IROP) + WQV(L,K,IROP)
        ENDDO  
      ENDIF  

      ! ****  PARAM 05  LOP - labile particulate organic phosphorus
      IF( ISKINETICS(5) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***    HYDROLYSIS  SETTLING
          WQF8 = - (WQKLPP(L) + WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQF8)
          
          ! *** Algae metabolism and predation source
          WQA8 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA8A(NAL) = (ALGAES(NAL).WQFPLB*WQBM(L,NAL) + ALGAES(NAL).WQFPLP*WQPR(L,NAL)) * WQO(L,19+NAL)
              WQA8 = WQA8 + WQA8A(NAL) * WQAPC(L) 
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA8A(NAL) = (ALGAES(NAL).WQFPLB*WQBM(L,NAL) + ALGAES(NAL).WQFPLP*WQPR(L,NAL)) * WQO(L,19+NAL)
              WQA8 = WQA8 + WQA8A(NAL) * WQAPC(L) * ALGAES(NAL).WQAPCM 
            ENDIF
          ENDDO
          
          ! ***  PT_SRC_LOADS    VOLUME  
          WQR8 = WQWPSL(L,K,ILOP) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,ILOP) + DTWQ*WQR8 + DTWQO2*WQA8 + WQF8*WQVO(L,K,ILOP)  
          
          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SLPOPZ(L,K)
        ENDDO 
        
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,ILOP)
          ENDDO  
        ELSE
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP       ATM WET DEP        VOLUME  
            WQR8 = (WQWDSL(L,K,ILOP) + WQATML(L,KC,ILOP))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR8  
          ENDDO
        ENDIF  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,ILOP) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,ILOP)   = WQVO(L,K,ILOP) + WQV(L,K,ILOP)
        ENDDO  
      ENDIF  

      ! ****  PARAM 06  DOP - dissolved organic phosphorus
      IF( ISKINETICS(6) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQF9 = - DTWQO2*WQKDOP(L)  
          WQKK(L) = 1.0 / (1.0 - WQF9)
          
          ! *** Algae metabolism and predation source
          WQA9 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE )THEN
              ! *** Phytoplankton
              WQA9A(NAL) = (ALGAES(NAL).WQFPDB*WQBM(L,NAL) + ALGAES(NAL).WQFPDP*WQPR(L,NAL)) * WQO(L,19+NAL)
              WQA9 = WQA9 + WQA9A(NAL) * WQAPC(L)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA9A(NAL) = (ALGAES(NAL).WQFPDB*WQBM(L,NAL) + ALGAES(NAL).WQFPDP*WQPR(L,NAL)) * WQO(L,19+NAL)
              WQA9 = WQA9 + WQA9A(NAL) * WQAPC(L) * ALGAES(NAL).WQAPCM 
            ENDIF
          ENDDO

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR9 = WQWPSL(L,K,IDOP) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,IDOP) + DTWQ*WQR9 + WQF9*WQVO(L,K,IDOP) + DTWQO2*(WQA9 + WQKRPP(L)*WQO(L,IROP) + WQKLPP(L)*WQO(L,ILOP))
           
          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDOPZ(L,K)
        ENDDO

        IF( K == KC )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP         ATM WET DEP       VOLUME  
            WQR9 = (WQWDSL(L,KC,IDOP) + WQATML(L,KC,IDOP))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR9  
          ENDDO
        ENDIF

        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IDOP) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,IDOP)   = WQVO(L,K,IDOP) + WQV(L,K,IDOP)
        ENDDO  
      ENDIF  

      ! ****  PARAM 07  P4D - total phosphate (PO4t)
      IF( ISKINETICS(7) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)
          
          ! *** Algal metabolism and predation source
          WQKK(L) = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE /= 0  )THEN
              ! *** Phytoplankton
              WQA10A(NAL) = (ALGAES(NAL).WQFPIB*WQBM(L,NAL) + ALGAES(NAL).WQFPIP*WQPR(L,NAL) - WQPA(L,NAL))*WQO(L,19+NAL)
              WQKK(L) = WQKK(L) + WQA10A(NAL)*WQAPC(L)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA10A(NAL) = (ALGAES(NAL).WQFPIB*WQBM(L,NAL) + ALGAES(NAL).WQFPIP*WQPR(L,NAL) - WQPA(L,NAL))*WQO(L,19+NAL)
              WQKK(L) = WQKK(L) + WQA10A(NAL)*WQAPC(L)*ALGAES(NAL).WQAPCM  
            ENDIF
          ENDDO 

          ! ***      PT_SRC_LOADS      VOLUME  
          WQRR(L) = WQWPSL(L,K,IP4D) * VOLWQ(L) * PSMLMULTIPLIER
           
          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + SPO4Z(L,K)
        ENDDO  

        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN                       ! *** Note - CAN ADD A PER CELL ARRAY OF 1 TO KC WHERE K=KSZ IS 1.0 AND ALL OTHER LAYERS = 0.0
            WQRR(L) = WQRR(L) + WQBFPO4D(L)*HPKI(L,K) ! *** Add in Benthic Flux
          ENDIF
        ENDDO  

        IF( K == KC )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***      ATM DRY DEP       ATM WET DEP        VOLUME  
            WQR10 = (WQWDSL(L,KC,IP4D) + WQATML(L,KC,IP4D))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR10  
          ENDDO
        ENDIF  

        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          WQRR(L) = WQVO(L,K,IP4D) + DTWQ*WQRR(L) + WQH10(L)*WQVO(L,K,IP4D) + DTWQO2*(WQKK(L) + WQKDOP(L)*WQO(L,IDOP))
        ENDDO  
  
        IF( K /= KC )THEN                              ! *** Note - no settling if IWQSRP = 0
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQT10(L)*WQVO(L,K+1,IP4D)
          ENDDO  
        ENDIF 
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQKKL = 1.0 / (1.0 - WQH10(L)) 
          WQV(L,K,IP4D) = MAX(WQRR(L)*WQKKL,0.0)
          WQO(L,IP4D)   = WQVO(L,K,IP4D) + WQV(L,K,IP4D)
        ENDDO  
      ENDIF  

      ! ****  PARAM 08  RON - refractory particulate organic nitrogen
      IF( ISKINETICS(8) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS     SETTLING
          WQI11 = - (WQKRPN(L) + WQRPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQI11)
          
          ! *** Algae metabolism and predation
          WQA11 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA11A(NAL) = (ALGAES(NAL).WQFNRB*WQBM(L,NAL) + ALGAES(NAL).WQFNRP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,19+NAL)  
              WQA11 = WQA11 + WQA11A(NAL)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA11A(NAL) = (ALGAES(NAL).WQFNRB*WQBM(L,NAL) + ALGAES(NAL).WQFNRP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,19+NAL)   
              WQA11 = WQA11 + WQA11A(NAL)
            ENDIF
          ENDDO

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR11 = WQWPSL(L,K,IRON) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,IRON) + DTWQ*WQR11 + DTWQO2*WQA11 + WQI11*WQVO(L,K,IRON)
          
          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SRPONZ(L,K)
        ENDDO 

        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,IRON)
          ENDDO  
        ELSE   ! K == KC
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
            WQR11 = (WQWDSL(L,KC,IRON) + WQATML(L,KC,IRON))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR11  
          ENDDO
        ENDIF
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IRON) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,IRON)   = WQVO(L,K,IRON) + WQV(L,K,IRON)
        ENDDO  
      ENDIF  

      ! ****  PARAM 09  LON - labile particulate organic nitrogen
      IF( ISKINETICS(9) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS     SETTLING
          WQJ12 = - (WQKLPN(L) + WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQJ12)
          
          ! *** Algae metabolism and predation
          WQA12 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA12A(NAL) = (ALGAES(NAL).WQFNLB*WQBM(L,NAL) + ALGAES(NAL).WQFNLP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,19+NAL)  
              WQA12 = WQA12 + WQA12A(NAL)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA12A(NAL) = (ALGAES(NAL).WQFNLB*WQBM(L,NAL) + ALGAES(NAL).WQFNLP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,19+NAL)  
              WQA12 = WQA12 + WQA12A(NAL)
            ENDIF
          ENDDO

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR12 = WQWPSL(L,K,ILON) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,ILON) + DTWQ*WQR12 + DTWQO2*WQA12 + WQJ12*WQVO(L,K,ILON)
                     

          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SLPONZ(L,K)
        ENDDO  
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,ILON)
          ENDDO  
        ELSE   ! K == KC
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP       ATM WET DEP     VOLUME  
            WQR12 = (WQWDSL(L,KC,ILON) + WQATML(L,KC,ILON))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR12  
          ENDDO
        ENDIF  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,ILON) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,ILON)   = WQVO(L,K,ILON) + WQV(L,K,ILON)
        ENDDO  
      ENDIF  

      ! ****  PARAM 10  DON - dissolved organic nitrogen
      IF( ISKINETICS(10) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQF13 = - DTWQO2*WQKDON(L)  
          WQKK(L) = 1.0 / (1.0 - WQF13)
          
          ! *** Algal metabolism and predation source
          WQA13 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA13A(NAL) = (ALGAES(NAL).WQFNDB*WQBM(L,NAL) + ALGAES(NAL).WQFNDP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,19+NAL)  
              WQA13 = WQA13 + WQA13A(NAL) 
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA13A(NAL) = (ALGAES(NAL).WQFNDB*WQBM(L,NAL) + ALGAES(NAL).WQFNDP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,19+NAL)  
              WQA13 = WQA13 + WQA13A(NAL) 
            ENDIF
          ENDDO
          
          ! ***    PT_SRC_LOADS    VOLUME  
          WQR13 = WQWPSL(L,K,IDON) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQVO(L,K,IDON) + DTWQ*WQR13 + WQF13*WQVO(L,K,IDON) + DTWQO2*(WQA13 + WQKRPN(L)*WQO(L,IRON)  + WQKLPN(L)*WQO(L,ILON))
           
          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDONZ(L,K)
        ENDDO
        
        IF( K == KC )THEN                          ! *** Note - add a check to see of the deposition concentrations are > 0.  Do this for all constituents
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
            WQR13 = (WQWDSL(L,KC,IDON) + WQATML(L,KC,IDON))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR13  
          ENDDO
        ENDIF  

        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IDON) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,IDON)   = WQVO(L,K,IDON) + WQV(L,K,IDON)
        ENDDO  
      ENDIF  
      
      ! ****  PARAM 11  NHX - ammonia nitrogen
      IF( ISKINETICS(11) == 1 )THEN 
        IF( IWQPSL < 2 )THEN
          DO LP = 1,LLWET(K,ND)                      ! *** Note - SKIP IF IWQPSL = 2
            L = LKWET(LP,K,ND)  
            ! ***      PT_SRC_LOADS     VOLUME  
            WQRR(L) = WQWPSL(L,K,INHX) * VOLWQ(L) * PSMLMULTIPLIER
          ENDDO
        ELSE
          ! *** Zero the loadings array
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = 0.0
          ENDDO   
        ENDIF

        DO LP = 1,LLWET(K,ND)                        ! *** Note - CAN ADD A PER CELL ARRAY OF 1 TO KC WHERE K=KSZ IS 1.0 AND ALL OTHER LAYERS = 0.0
          L = LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQRR(L) = WQRR(L) + WQBFNH4(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          ENDIF  
        ENDDO  

        IF( K == KC )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
            WQR14 = (WQWDSL(L,KC,INHX) + WQATML(L,KC,INHX))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR14  
          ENDDO
        ENDIF  

        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQF14 = - DTWQO2*WQNIT(L)  
          WQKK(L) = 1.0 / (1.0 - WQF14)
          
          ! *** Algal metabolism and predation source
          WQA14 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA14A(NAL) = ALGAES(NAL).WQFNIB*WQBM(L,NAL) + ALGAES(NAL).WQFNIP*WQPR(L,NAL) - WQPN(L,1)*WQPA(L,NAL)  
              WQA14 = WQA14 + WQA14A(NAL)*ALGAES(NAL).WQANCA*WQO(L,19+NAL)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA14A(NAL) = ALGAES(NAL).WQFNIB*WQBM(L,NAL) + ALGAES(NAL).WQFNIP*WQPR(L,NAL) - WQPN(L,1)*WQPA(L,NAL)  
              WQA14 = WQA14 + WQA14A(NAL)*ALGAES(NAL).WQANCA*WQO(L,19+NAL)
            ENDIF
          ENDDO
 
          WQRR(L) = WQV(L,K,INHX) + DTWQ*WQRR(L) + WQF14*WQVO(L,K,INHX) + DTWQO2*(WQA14 + WQKDON(L)*WQO(L,IDON))

          ! *** Coupling to zooplankton - DKT
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SNH4Z(L,K)
          
          WQV(L,K,INHX) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,INHX)   = WQVO(L,K,INHX) + WQV(L,K,INHX)
        ENDDO  
      ENDIF  

      ! ****  PARAM 12  NOX - nitrate nitrogen
      IF( ISKINETICS(12) == 1 )THEN 
        IF( IWQPSL < 2 )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***      PT_SRC_LOADS      VOLUME  
            WQRR(L) = WQWPSL(L,K,INOX) * VOLWQ(L) * PSMLMULTIPLIER
          ENDDO  
        ELSE
          ! *** Zero the loadings array
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = 0.0
          ENDDO   
        ENDIF
    
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQRR(L) = WQRR(L) + WQBFNO3(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          ENDIF  
        ENDDO  
    
        IF( K == KC )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP         ATM WET DEP       VOLUME  
            WQR15 = (WQWDSL(L,KC,INOX) + WQATML(L,KC,INOX))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR15  
          ENDDO
        ENDIF  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          ! *** Algal predation source
          WQA15 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE  )THEN
              ! *** Phytoplankton
              WQA15A(NAL) = (WQPN(L,NAL) - 1.0)*WQPA(L,NAL) * ALGAES(NAL).WQANCA * WQO(L,19+NAL)  
              WQA15 = WQA15 + WQA15A(NAL)
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN
              ! *** Macrophytes and periphyton
              WQA15A(NAL) = (WQPN(L,NAL) - 1.0)*WQPA(L,NAL) * ALGAES(NAL).WQANCA * WQO(L,19+NAL)  
              WQA15 = WQA15 + WQA15A(NAL)
            ENDIF
          ENDDO
  
          WQB15 = WQV(L,K,INOX) + DTWQ*WQRR(L) + DTWQO2*( WQA15 - WQANDC*WQDENIT(L)*WQO(L,IDOC) + WQNIT(L)*WQO(L,INHX))

          WQV(L,K,INOX) = MAX(WQB15,0.0)
          WQO(L,INOX)   = WQVO(L,K,INOX) + WQV(L,K,INOX)
        ENDDO  
      ENDIF  

      ! ****  PARAM 13  SUU - particulate biogenic silica
      IF( ISKINETICS(13) == 1 )THEN  
        IF( IWQSI == 1 )THEN  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            
            ! *** Algal metabolism and predation source
            WQA16D = 0.0
            WQM16 = 0.0
            ILAST = 1
            DO NAL = 1,NALGAE
              IF( ALGAES(NAL).ISILICA )THEN
                WQA16D = WQA16D + (ALGAES(NAL).WQFSPD*WQBM(L,NAL) + ALGAES(NAL).WQFSPP*WQPR(L,NAL)) * ALGAES(NAL).WQASC * WQO(L,19+NAL)
                WQM16 =  WQM16  - (WQKSUA(IWQT(L)) + WQBSET(L,1,NAL)) * DTWQO2
                ILAST = NAL
              ENDIF
            ENDDO
            WQKK(L) = 1.0 / (1.0 - WQM16)
            
            ! ***    PT_SRC_LOADS      VOLUME  
            WQR16 = WQWPSL(L,K,ISUU) * VOLWQ(L) * PSMLMULTIPLIER
            WQRR(L) = WQV(L,K,ISUU) + DTWQ*WQR16 + DTWQO2*WQA16D + WQM16*WQVO(L,K,ISUU)
             
            ! *** Coupling to zooplankton - DKT
            IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SSUZ(L,K)
          ENDDO  
          
          IF( K /= KC )THEN  
            ! *** Add in settling from above
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQBSET(L,2,ILAST)*WQVO(L,K+1,ISUU)
            ENDDO  
          ELSE
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***     ATM DRY DEP         ATM WET DEP       VOLUME  
              WQR16 = (WQWDSL(L,KC,ISUU) + WQATML(L,KC,ISUU))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR16  
            ENDDO
          ENDIF  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQV(L,K,ISUU) = MAX(WQRR(L)*WQKK(L),0.0)
            WQO(L,ISUU)   = WQVO(L,K,ISUU) + WQV(L,K,ISUU)
          ENDDO  
        ENDIF  
      ENDIF  

      ! ****  PARAM 14  SAA - dissolved available silica
      IF( ISKINETICS(14) == 1 )THEN  
        IF( IWQSI == 1 )THEN  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            WQKK(L) = 0.
            DO NAL = 1,NALGAE
              IF( ALGAES(NAL).ISILICA )THEN
                WQKK(L) = WQKK(L) + (ALGAES(NAL).WQFSID*WQBM(L,NAL) + ALGAES(NAL).WQFSIP*WQPR(L,NAL) - WQPA(L,NAL)) * ALGAES(NAL).WQASC * WQO(L,19+NAL)
              ENDIF
            ENDDO
            ! ***      PT_SRC_LOADS      VOLUME  
            WQRR(L) = WQWPSL(L,K,ISAA) * VOLWQ(L) * PSMLMULTIPLIER
          ENDDO  
      
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            IF( K == KSZ(L) )THEN
              WQRR(L) = WQRR(L) + WQBFSAD(L)*HPKI(L,K)   ! *** Add in Benthic Flux
            ENDIF  
          ENDDO  
          
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQV(L,K,ISAA) + DTWQ*WQRR(L) + WQN17(L)*WQVO(L,K,ISAA) + DTWQO2*( WQKK(L) + WQKSUA(IWQT(L))*WQO(L,ISUU))
             
            ! *** Coupling to zooplankton - DKT
            IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SSAZ(L,K)
          ENDDO  

          IF( K /= KC )THEN  
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQT17(L)*WQVO(L,K+1,ISAA)  
            ENDDO  
          ELSE
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
              WQR17 = (WQWDSL(L,KC,ISAA) + WQATML(L,KC,ISAA))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR17  
            ENDDO
          ENDIF  
          
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - WQN17(L))
            WQV(L,K,ISAA) = MAX(WQRR(L)*WQKK(L),0.0)
            WQO(L,ISAA)   = WQVO(L,K,ISAA) + WQV(L,K,ISAA)
          ENDDO  
        ENDIF  
      ENDIF  

      ! ****  PARAM 15  COD - chemical oxygen demand
      IF( ISKINETICS(15) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQKK(L) = 1.0 / (1.0 - WQO18(L))  
          
            ! ***    PT_SRC_LOADS      VOLUME  
          WQRR(L) = WQWPSL(L,K,ICOD) * VOLWQ(L) * PSMLMULTIPLIER
        ENDDO  
        
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQRR(L) = WQRR(L) + WQBFCOD(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          ENDIF  
        ENDDO  
        
        IF( K == KC )THEN
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP         ATM WET DEP       VOLUME  
            WQR18 = (WQWDSL(L,KC,ICOD) + WQATML(L,KC,ICOD))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR18  
          ENDDO
        ENDIF  
         
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRR(L) = WQV(L,K,ICOD) + DTWQ*WQRR(L) + WQO18(L)*WQV(L,K,ICOD)  
          WQV(L,K,ICOD) = MAX(WQRR(L)*WQKK(L),0.0)
          WQO(L,ICOD)   = WQVO(L,K,ICOD) + WQV(L,K,ICOD)  
        ENDDO  
      ENDIF

      ! ****  PARAM 16  DOX - dissolved oxygen

      ! ***  1) ROC - Refractory particulate organic carbon
      ! ***  2) LOC - Labile particulate organic carbon
      ! ***  3) DOC - Dissolved organic carbon
      ! ***  4) ROP - Refractory particulate organic phosphorus
      ! ***  5) LOP - Labile particulate organic phosphorus
      ! ***  6) DOP - Dissolved organic phosphorus
      ! ***  7) P4D - Total phosphate
      ! ***  8) RON - Refractory particulate organic nitrogen
      ! ***  9) LON - Labile particulate organic nitrogen
      ! *** 10) DON - Dissolved organic nitrogen
      ! *** 11) NHX - Ammonia nitrogen
      ! *** 12) NOX - Nitrate nitrogen
      ! *** 13) SUU - Particulate biogenic silica
      ! *** 14) SAA - Dissolved available silica
      ! *** 15) COD - Chemical oxygen demand
      ! *** 16) DOX - Dissolved oxygen
      ! *** 17) TAM - Total active metal
      ! *** 18) FCB - Fecal coliform bacteria
      ! *** 19) CO2 - Dissolved carbon dioxide
      ! *** 20) BIO - Any number of biota classes (NALGAE), phytoplankton, macrophytes and zooplankton
      ! ***           Can include green, blue-green, cyanobacteria, etc.  

      ! *** EE7.2 REMOVED THE DO COMPONENT ANALYSIS

      IF( ISKINETICS(16) == 1 )THEN 
        
        ! *** Point Source Loading 
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRR(L) = WQWPSL(L,K,IDOX) * VOLWQ(L) * PSMLMULTIPLIER
        ENDDO
  
        ! *** Handle Surface Processes
        IF( K == KC )THEN  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L)) 
            ! ***                  ATM DRY DEP        ATM WET DEP       VOLUME  
            WQRR(L) = WQRR(L) + (WQWDSL(L,KC,IDOX) + WQATML(L,KC,IDOX))*VOLWQ(L)
            ! *** Reaeration - Atm to water flux  (Offset later in O2 Mass Balance by WQRea term)
            WQRR(L) = WQRR(L) + WQKRDOS(L)
          ENDDO
        ELSE  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0
          ENDDO
        ENDIF 

        ! *** Bottom Processes 
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            TEMFAC = STEMFAC**(TEM(L,KSZ(L)) - 20.)
            WQRR(L) = WQRR(L) + TEMFAC*WQBFO2(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          ENDIF
        ENDDO
  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          DTWQxH = DTWQ*HPK(L,K)  
          DTWQxH2= DTWQO2*HPK(L,K)

          ! *** Photosynthesis
          IF( WQI0  <=  0.001 )THEN
            DO NAL = 1, NALGAE
              WQTTAA(NAL) = 0.0
            ENDDO
          ELSE
            DO NAL = 1,NALGAE
              IF( ALGAES(NAL).ISMOBILE )THEN
                ! *** Phytoplankton
                WQTTAA(NAL) = (1.3 - 0.3*WQPN(L,NAL)) * WQPA(L,NAL) 
              ENDIF
              IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN  
                ! *** Macrophytes and periphyton
                WQTTAA(NAL) = (1.3 - 0.3*WQPN(L,NAL)) * WQPA(L,NAL)
              ENDIF
            ENDDO               
          ENDIF
          
          ! *** Respiration and final net O2 rate
          WQA19 = 0.
          IZ = IWQZMAP(L,K)  
          DO NAL = 1,NALGAE
            ! *** Respiration 
            IF( ALGAES(NAL).ISMOBILE )THEN
              ! *** Phytoplankton
              WQRes = (1.0 - ALGAES(NAL).WQFCDB)*O2WQ(L)*WQBM(L,NAL)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)       
            ENDIF
            IF( .NOT. ALGAES(NAL).ISMOBILE .AND. (K >= LAYERBOT(NAL,L) .AND. K <= LAYERTOP(NAL,L)) )THEN  
              ! *** Macrophytes and periphyton
              WQRes = (1.0 - ALGAES(NAL).WQFCDB)*O2WQ(L)*WQBM(L,NAL)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)
            ENDIF
            
            ! *** Net algal O2 production due to photosynthesis and respiration
            WQA19A(NAL) = WQTTAA(NAL)*ALGAES(NAL).WQAOCRP - WQRes*ALGAES(NAL).WQAOCRR
            WQA19A(NAL) = WQA19A(NAL)*WQO(L,19+NAL)                                      ! *** Net O2 production by group
            WQA19 = WQA19 + WQA19A(NAL)                                                  ! *** Total O2 net production
          ENDDO

          ! WQA19                                  ! *** Total Net Respiration/Photosynthesis
          WQSUM = DTWQ*WQRR(L)                     ! *** Sum of Loadings/Demands
          WQRea = DTWQO2*WQP19(L)*WQV(L,K,IDOX)    ! *** Reaeration - Offsetting Flux (WQP19 is negative)
          WQDOC = WQAOCR*WQKHR(L)*WQO(L,IDOC)      ! *** DOC
          WQNH3 = WQAONT*WQNIT(L)*WQO(L,INHX)      ! *** Ammonia
          WQCOD = WQO18(L)*WQO(L,ICOD)             ! *** COD
          WQRR(L) = WQV(L,K,IDOX) + WQSUM  + WQCOD  + WQRea + DTWQO2*(WQA19 - WQDOC - WQNH3)
          
          IF( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDOZ(L,K)    ! *** Coupling to zooplankton

          WQV(L,K,IDOX) = MAX(WQRR(L)*WQKK(L),0.0)

        ENDDO
        
      ENDIF  ! *** END OF ISKINETICS(16)

      ! ****  PARAM 17  TAM - total active metal
      IF( ISKINETICS(17) == 1 )THEN  
        IF( IWQSRP == 1 )THEN  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQR20(L) = WQWPSL(L,K,17)*VOLWQ(L) * PSMLMULTIPLIER + (WQV(L,K,ITAM) - WQTAMP(L,K))*WQWSSET(L,1)  
          ENDDO

          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            IF(K == KSZ(L)) WQR20(L) = WQR20(L) + WQTDTAM(IWQT(L))*HPKI(L,K)/(WQKHBMF + O2WQ(L) + 1.E-18)
          ENDDO

          IF( K /= KC )THEN
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQR20(L) = WQR20(L) + (WQV(L,K+1,ITAM) - WQTAMP(L,K+1)) * WQWSSET(L,2)  
            ENDDO 
          ELSE     ! K == KC
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQR20(L) = WQR20(L) + (WQWDSL(L,KC,ITAM) + WQATML(L,KC,ITAM))*VOLWQ(L)
            ENDDO
          ENDIF 

          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQT20 = - DTWQ*WQWSSET(L,1)    ! *** DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQT20)  
            WQRR(L) = WQV(L,K,ITAM) + DTWQ*WQR20(L) + WQT20*WQV(L,K,ITAM)  
          ENDDO  
          
          IF( K /= KC )THEN  
            ! *** Add in settling from above
            DO LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQWSSET(L,2)*WQVO(L,K+1,ITAM)
            ENDDO  
          ENDIF  
          
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQV(L,K,ITAM) = MAX(WQRR(L)*WQKK(L),0.0)
            WQO(L,ITAM)   = WQVO(L,K,ITAM) + WQV(L,K,ITAM)
          ENDDO  
        ENDIF  
      ENDIF  

      ! ****  PARAM 18  FCB - fecal coliform bacteria
      IF( ISKINETICS(18) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQKK(L) = WQTD2FCB(IWQT(L))  
          
          ! ***   PT SRC LOADS     VOLUME  
          WQR21 = WQWPSL(L,K,IFCB)*VOLWQ(L) * PSMLMULTIPLIER
          
          IF( K == KC )THEN
            ! ***            ATM DRY DEP        ATM WET DEP       VOLUME  
            WQR21 = WQR21 + (WQWDSL(L,K,IFCB) + WQATML(L,KC,IFCB))*VOLWQ(L)
          ENDIF
      
          WQRR(L) = WQV(L,K,IFCB)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21
          WQV(L,K,IFCB) = MAX(WQRR(L)*WQKK(L),0.0)
        ENDDO  
      ENDIF

      ! ***************************************************************************
      ! *** SHELLFISH KINETICS *** shellfish 2
      IF(  ISFFARM > 0 .AND. NSF > 0 )THEN
        DO L=LF,LL
          DO NQ = 1,NSF
            CALL SHELLFISH_GROWTH(NQ, SF(NQ), L, K)
          ENDDO
        ENDDO
        CALL SHELLFISH_HARVEST(LF,LL)
      ENDIF

      ! ****  PARAM 19 DISSOLVED CARBON DIOXIDE

      !C THE FOLLOWING ARRAYS WERE ADDED TO KEEP TRACK OF THE VARIOUS COMPONENT  
      !C OF DISSOLVED CARBON DIOXIDE.  
      !C THE ARRAY DESCRIPTIONS ARE:  
      !C  XCDOKAR(L,K) = CDO. COMPONENT FOR REAERATION  
      !C  XCDODOC(L,K) = CDO. COMPONENT FOR DISS. ORG. CARBON DECAY  
      !C  XCDOPPB(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
      !C  XCDORRB(L,K) = CDO. COMPONENT FOR RESPIRATION OF TOTAL CHLOROPHYLL  
      !C  XCDOPPM(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF MACROALGAE  
      !C  XCDORRM(L,K) = CDO. COMPONENT FOR RESPIRATION OF MACROALGAE  
      !C  XCDOALL(L,K) = SUM OF THE ABOVE 6 CDO. COMPONENTS  
      
      ! *** CO2
      IF( ISKINETICS(19) == 1 )THEN  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRR(L) = WQWPSL(L,K,ICO2) * VOLWQ(L) * PSMLMULTIPLIER  
          TMP22 = WQRR(L)*DTWQ*HPK(L,K) 
        ENDDO  
  
        ! *** Handle Surface Processes
        IF( K == KC )THEN  
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP22(L))
            
            ! ***                  ATM DRY DEP     ATM WET DEP     VOLUME  
            WQRR(L) = WQRR(L) + (WQWDSL(L,KC,ICO2) + WQATML(L,KC,ICO2)*VOLWQ(L))  
            
            ! *** Reaeration
            WQRR(L) = WQRR(L) + WQKRCDOS(L) 
          ENDDO
        ELSE
          DO LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0
          ENDDO
        ENDIF   
  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          DTWQxH  = DTWQ*HPK(L,K)  
          DTWQxH2 = DTWQO2*HPK(L,K)
          IF( WQI0  <=  0.001 )THEN
            DO NAL = 1,NALGAE
              IF( ALGAES(NAL).ISMOBILE )THEN
                WQTTAA(NAL) = 0.0
              ENDIF
            ENDDO
          ELSE
            DO NAL = 1,NALGAE
              ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL 
              IF( ALGAES(NAL).ISMOBILE )THEN
                WQTTAA(NAL) = (1.3 - 0.3*WQPN(L,NAL)) * WQPA(L,NAL) 
              ENDIF
            ENDDO
          ENDIF
          
          ! *** RESPIRATION OF TOTAL CHLOROPHYLL - general algal group
          WQA22 = 0.
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISMOBILE )THEN
              ! *** Phytoplankton
              XMRM = (1.0 - ALGAES(NAL).WQFCDB)*O2WQ(L)*WQBM(L,NAL)/(ALGAES(NAL).WQKHRA(1) + O2WQ(L) + 1.E-18)  
              WQA22A(NAL) = WQTTAA(NAL) - XMRM
              WQA22 = WQA22 + 3.67*WQA22A(NAL)*WQO(L,19+NAL)  !VB 3.67 CONVERTS g CARBON TO g CO2
            ENDIF
          ENDDO

          ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
          ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
          ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
          ! *** CO2 Mass Balance
          ! WQA22                                               ! *** Total Net Respiration/Photosynthesis
          WQCDSUM = DTWQ*WQRR(L)                                ! *** Sum of Loadings/Demands
          WQCDRea = WQP22(L)*WQV(L,K,ICO2)                      ! *** Reaeration
          WQCDDOC = (WQKHR(L) + WQDENIT(L))*WQO(L,IDOC)*3.67    ! *** DOC FROM HYDROLYSIS AND DENITRIFICATION 3.67 CONVERTS G CARBON TO G CO2    

          WQRR(L) = WQV(L,K,ICO2) + WQCDSUM + DTWQO2*(-WQA22 + WQCDDOC + WQCDRea)
          WQV(L,K,ICO2) = MAX(WQRR(L)*WQKK(L),0.0)
        ENDDO  
      ENDIF  
    ENDDO  ! *** END OF THE KC LOOP

    ! ***************************************************************************
    ! *** SUM OF SHELLFISH FOR ALL LAYERS *** shellfish 4
    IF(  ISFFARM > 0 .AND. NSF > 0 )THEN
      DO L=LF,LL
        CALL SHELLFISH_LAYERSUM(L)
      ENDDO
    ENDIF

    ! *** PHOSPHATE AND SILICA SORPTION OPTIONS  
    IF( IWQSRP == 1 )THEN  
      ! *** Sorption Option: TAM
      DO K = 1,KC  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          O2WQ(L) = MAX(WQV(L,K,IDOX), 0.0)  
          WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ(L)), WQV(L,K,ITAM) )  
          WQTAMP(L,K) = WQV(L,K,ITAM) - WQTAMD  
          WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*WQTAMP(L,K))  
          WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*WQTAMP(L,K))  
        ENDDO  
      ENDDO  
    ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN
      ! *** Sorption Option: Sediments
      DO K = 1,KC  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*SEDT(L,K))  
          WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*SEDT(L,K))  
        ENDDO  
      ENDDO  
    ELSE  
      DO K = 1,KC  
        DO LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQPO4D(L,K) = WQV(L,K,IP4D)  
          WQSAD(L,K)  = WQV(L,K,ISAA)  
        ENDDO  
      ENDDO  
    ENDIF  

    ! ***************************************************************************
    ! COUPLING TO SEDIMENT MODEL  
    ! EVALUATE DEP. FLUX USING NEW VALUES CAUSE IMPLICIT SCHEME IS USED IN SPM
    IF( IWQBEN == 1 )THEN  
      DO LP = LF,LL
        L = LWET(LP)
        IMWQZ = IWQZMAP(L,KSZ(L))
        DO NAL = 1,NALGAE
          IF( ALGAES(NAL).ISMOBILE )THEN
            ! *** Phytoplankton
            WQDFB(L,NAL) = SCB(L)*ALGAES(NAL).WQWS(IMWQZ)*WQV(L,KSZ(L),19+NAL) 
          ELSE
            WQDFB(L,NAL) = ALGAES(NAL).WQWS(IMWQZ)*HPKI(L,K)*WQV(L,KSZ(L),19+NAL)
          ENDIF
        ENDDO       
        WQDFRC(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,KSZ(L),IROC)   !< Depositional flux of RPOC 
        WQDFLC(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,KSZ(L),ILOC)   !< Depositional flux of LPOC   
        WQDFRP(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,KSZ(L),IROP)   !< Depositional flux of RPOP
        WQDFLP(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,KSZ(L),ILOP)   !< Depositional flux of LPOP   
        WQDFRN(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,KSZ(L),IRON)   !< Depositional flux of RPON
        WQDFLN(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,KSZ(L),ILON)   !< Depositional flux of LPON   
        
        WQDFSI(L) = 0.0
        IF( IWQSI == 1 )THEN
          DO NAL = 1,NALGAE
            IF( ALGAES(NAL).ISILICA > 0 )THEN
               WQDFSI(L) = WQDFSI(L) + SCB(L)*ALGAES(NAL).WQWS(IMWQZ)*WQV(L,KSZ(L),ISUU)  
            ENDIF
          ENDDO       
        ENDIF      
      ENDDO  

      IF( IWQSRP == 1 )THEN  
        DO LP = LF,LL
          L = LWET(LP)
          IMWQZ = IWQZMAP(L,KSZ(L))  
          WQDFLP(L) = SCB(L)*(WQDFLP(L) + WQWSS(IMWQZ)*(WQV(L,KSZ(L),IP4D) - WQPO4D(L,1)))  
          IF( IWQSI == 1 ) WQDFSI(L) = SCB(L)*(WQDFSI(L) + WQWSS(IMWQZ)*(WQV(L,KSZ(L),ISAA) - WQSAD(L,1)))  
        ENDDO  
      ELSE IF( IWQSRP == 2 )THEN  
        DO LP = LF,LL
          L = LWET(LP)
          WQDFLP(L) = SCB(L)*(WQDFLP(L)+WSEDO(NS)*( WQV(L,KSZ(L),IP4D)-WQPO4D(L,1)))  
          IF( IWQSI == 1 ) WQDFSI(L) = SCB(L)*(WQDFSI(L) + WSEDO(NS)*(WQV(L,KSZ(L),ISAA) - WQSAD(L,1)))  
        ENDDO  
      ENDIF  
    ENDIF  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO
  
  ! *** ZERO THIN LAYERS WHEN DELTA T WQ IS LARGE
  IF( KC > 1 .AND. ISDRY > 0 )THEN
    IFLAG = 0
    !$OMP PARALLEL DO PRIVATE(ND,LF,LL,LP,L,NW,K) REDUCTION(+:IFLAG)
    DO ND=1,NDM 
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      DO LP = LF,LL
        L = LWET(LP)
        IF( HP(L) < 2.*HDRY )THEN
          DO NW = 1,NWQV
            IF( ISKINETICS(NW) > 0 )THEN
              DO K = KSZ(L),KC
                IF( WQV(L,K,NW) < 0.0 )THEN
                  WQV(L,K,NW) = 0.0
                  IFLAG = IFLAG + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    if( process_id == master_id )THEN 
      IF( IFLAG > 0 )THEN
        OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
        WRITE(8,'(A,F12.4,I10)')' WARNING!  NEGATIVE WQ CONCENTRATIONS: TIMEDAY, # OF VALUES < 0.0:',TIMEDAY,IFLAG
        CLOSE(8)
      ENDIF
    endif
  ENDIF
  
  ! *** RESTORE PREVIOUS WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO NQ = 0,NWQVM
    DO K = 1,KC
      DO IOBC = 1,NBCSOP  
        L = LOBCS(IOBC)
        WQV(L,K,NQ) = WQOLD(IOBC,K,NQ)
      ENDDO
    ENDDO  
  ENDDO  
  
  ! ***************************************************************************
  ! *** DIURNAL DO ANALYSIS  
  IF( WQHRAVG >  0.0 )THEN 
    NDDOCNT = NDDOCNT + 1
    DO K=1,KC  
      DO L=2,LA
        DDOAVG(L,K) = DDOAVG(L,K) + WQV(L,K,IDOX)*DTWQ
        DDOMAX(L,K) = MAX(DDOMAX(L,K),WQV(L,K,IDOX))  
        DDOMIN(L,K) = MIN(DDOMIN(L,K),WQV(L,K,IDOX))  
      ENDDO  
    ENDDO  

    IF( TIMEDAY >= WQHRAVG )THEN  
      NDDOCNT = 0  

      write(mpi_log_unit, '(a)') '*************************************************************'
      write(mpi_log_unit, '(a,f10.2)') 'Water Quality D.O. Analysis : ',TIMEDAY
      Call WriteBreak(mpi_log_unit)

      TVAL1 = 1./((TIMEDAY-AVGLAST)*86400.)
      DO L=2,LA
        DDOAVG(L,K) = DDOAVG(L,K)*TVAL1
        WRITE(mpi_log_unit,1112)  Map2Global(L).IG, Map2Global(L).JG, ((K, DDOMIN(L,K), DDOAVG(L,K), DDOMAX(L,K)),K=1,KC)  
      ENDDO  
      
      DO K=1,KC  
        DO L=2,LA  
          DDOAVG(L,K) = 0.0
          DDOMAX(L,K) = -1.E6  
          DDOMIN(L,K) = 1.E6  
        ENDDO  
      ENDDO  
      AVGLAST = TIMEDAY
      AVGNEXT = TIMEDAY + WQHRAVG
    ENDIF  
  ENDIF  

  ! ***************************************************************************
  ! *** LIGHT EXTINCTION ANALYSIS  
  IF( NDLTAVG >= 1 )THEN 
    IF( process_id == master_id )THEN
      OPEN(1,FILE=OUTDIR//'LIGHT.OUT',POSITION='APPEND')  
    endif
    NDLTCNT=NDLTCNT+1  
    NSTPTMP=NDLTAVG*NTSPTC/2  
    RMULTMP=1./FLOAT(NSTPTMP)  
    DO K=1,KC  
      DO L=2,LA 
        IF(  ISTRAN(6) > 0  )THEN
          RLIGHT1=WQKEB(IWQZMAP(L,K))+WQKETSS*SEDT(L,K) 
        ELSE
          RLIGHT1=WQKEB(IWQZMAP(L,K))
        ENDIF
        XMRM = WQKECHL*WQCHL(L,K)  
        IF( WQKECHL  < 0.0 )THEN  
          XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)  
        ENDIF  
        RLIGHT2 = XMRM  
        RLIGHTT(L,K)=RLIGHTT(L,K)+RLIGHT1  
        RLIGHTC(L,K)=RLIGHTC(L,K)+RLIGHT1+RLIGHT2  
      ENDDO  
    ENDDO  
    IF( NDLTCNT == NSTPTMP )THEN  
      NDLTCNT=0  
      TIME=TIMESEC/TCON  

      DO K=1,KC  
        DO L=2,LA
          RLIGHTT(L,K)=RMULTMP*RLIGHTT(L,K)  
          RLIGHTC(L,K)=RMULTMP*RLIGHTC(L,K)  
        ENDDO  
      ENDDO 
      
      IF( process_id == master_id )THEN
        WRITE(1,1111)NITER,TIME  
        DO L=2,LA  
          WRITE(1,1113)IL(L),JL(L),(RLIGHTT(L,K),K=1,KC),(RLIGHTC(L,K),K=1,KC)  
        ENDDO  
      endif
      DO K=1,KC  
        DO L=2,LA  
          RLIGHTT(L,K)=0.  
          RLIGHTC(L,K)=0.  
        ENDDO  
      ENDDO  
    ENDIF  
    CLOSE(1)  
  ENDIF  

  ! ***************************************************************************
  ! INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:  
  TIMTMP = TIMESEC/TCON  
  TIMESUM3 = TIMESUM3 + TIMTMP  
  
  1111 FORMAT(I12,F10.4)  
  1112 FORMAT(2I7,200(I4,3F7.2))  
  1113 FORMAT(2I7,12E12.4)  
  1414 FORMAT(I12,11E12.4)  

  RETURN    
  END SUBROUTINE WQSKE1
    
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSKE2
  !
  !> @details  
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSKE2

  PRINT *,'WQSKE2 HAS BEEN REMOVED.  NEEDS CODE IF IWQSKE=2 IMPLEMENTED'
  
  RETURN
  
  END  SUBROUTINE WQSKE2
      
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSKE3
  !
  !> @details  
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSKE3

  PRINT *,'WQSKE3 HAS BEEN REMOVED.  NEEDS CODE IF IWQSKE=3 IMPLEMENTED'
  
  RETURN
  
  END  SUBROUTINE WQSKE3
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSKE4
  !
  !> @details  
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSKE4

    PRINT *,'WQSKE4 HAS BEEN REMOVED.  NEEDS CODE IF IWQSKE=4 IMPLEMENTED'
  
    RETURN
  
    END SUBROUTINE WQSKE4
  
    !---------------------------------------------------------------------------!
    ! EFDC+ Developed by DSI, LLC.
    !---------------------------------------------------------------------------!
    ! Subroutine: Subroutine WQ_QC
    !
    !> @details  Write out a summary of the WQ settings
    !---------------------------------------------------------------------------!
    !---------------------------------------------------------------------------!
    SUBROUTINE WQ_QC

      IMPLICIT NONE

  !  WRITE(2,83)': FREQUENCY OF DIURNAL DO OUTPUT (IN DT UNIT) =', IWQDIUDT
  !  WRITE(2,83)'* IWQDT (DTWQ(D) = DT(S)*IWQDT/86400)        = ', IWQDT
  !
  !  WRITE(2,80)'* FULL VERSION WITH 21 VARIABLES IS ACTIVATED     '
  !  IF( IWQBEN == 1 )THEN
  !    WRITE(2,80)'* SEDIMENT PROCESS MODEL IS ACTIVATED             '
  !  ELSE IF( IWQBEN == 0 )THEN
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY CONSTANT BENTHIC FLUX IS SPECIFIED   '
  !  ELSE IF( IWQBEN == 2 )THEN
  !    WRITE(2,80)'* SPATIALLY AND/OR TEMPORALLY-VARYING BENTHIC FLUX SPECIFIED'
  !  ELSE
  !    CALL STOPP('** ERROR!!! INVALID IWQBEN VALUE **')
  !  ENDIF
  !
  !  IF( IWQSI == 1 )THEN
  !    WRITE(2,80)'* SILICA STATE VARIABLES (SU & SA) ARE MODELED    '
  !  ELSE
  !    WRITE(2,80)'* NO SILICA (SU & SU) LIMITATION                  '
  !  ENDIF
  !
  !  IF( IWQFCB == 1 )THEN
  !    WRITE(2,80)'* FCB (FECAL COLIFORM BACTERIA) IS MODELED        '
  !  ELSE
  !    WRITE(2,80)'* FCB (FECAL COLIFORM BACTERIA) IS NOT MODELED    '
  !  ENDIF
  !
  !  IF( IWQSRP == 1 )THEN
  !    WRITE(2,80)'* TAM IS USED FOR SORPTION OF PO4T/SA: MODEL TAM  '
  !  ELSE IF( IWQSRP == 2 )THEN
  !    WRITE(2,80)'* TSS IS USED FOR SORPTION OF PO4T/SA: MODEL TSS  '
  !    IF( ISTRAN(6) /= 1) CALL STOPP('ERROR! INCOMPATIBLE ISTRAN(6)/IWQSRP')
  !  ELSE
  !    WRITE(2,80)'* NO SORPTION OF PO4T/SA: MAY MODEL TSS & NO TAM  '
  !  ENDIF
  !
  !  IF( IWQSTOX == 1 )THEN
  !    WRITE(2,80)'* SALINITY TOXICITY IS APPLIED TO CYANOBACTERIA   '
  !  ELSE
  !    WRITE(2,80)'* NO SALINITY TOXICITY: SALTWATER CYANOBACTERIA   '
  !  ENDIF
  !
  !  IF( IWQKA(1) == 0 )THEN
  !    WRITE(2,80)'* USER-SPECIFIED CONSTANT REAERATION SET TO WQKRO '
  !    WRITE(2,80)'*   REAERATION DUE TO WIND SET TO ZERO            '
  !  ENDIF
  !
  !  IF( IWQKA(1) == 1 )THEN
  !    WRITE(2,80)'* USER-SPECIFIED CONSTANT REAERATION SET TO WQKRO '
  !    WRITE(2,80)'*   REAERATION DUE TO WIND ADDED TO WQKRO         '
  !  ENDIF
  !  IF( IWQKA(1) == 2 )THEN
  !    WRITE(2,80)'* OCONNOR-DOBBINS REAERATION FORMULA IS USED      '
  !  ENDIF
  !  IF( IWQKA(1) == 3 )THEN
  !    WRITE(2,80)'* OWENS & GIBBS (1964) REAERATION FORMULA IS USED '
  !  ENDIF
  !  IF( IWQKA(1) == 4 )THEN
  !    WRITE(2,80)'* MODIFIED OWENS & GIBBS REAERATION IS USED       '
  !  ENDIF
  !
  !  IF( IWQVLIM == 0 )THEN
  !    WRITE(2,80)'* MACROALGAE GROWTH IS NOT LIMITED BY VELOCITY    '
  !  ENDIF
  !  IF( IWQVLIM == 1 )THEN
  !    WRITE(2,80)'* MACROALGAE VELOCITY LIMIT, MICHAELIS-MENTON EQU.'
  !  ENDIF
  !  IF( IWQVLIM == 2 )THEN
  !    WRITE(2,80)'*MACROALGAE VEL. LIMIT, 5-PARAM LOGISTIC FUNCTION'
  !  ENDIF
  !
  !  WRITE(2,83)'* # OF ZONES FOR SPATIALLY VARYING PARAMETERS =',IWQZ
  !  IF( IWQZ > NWQZ) CALL STOPP('ERROR!! IWQZ SHOULD BE <= NWQZ')
  !
  !  IF( IWQNC > 0 )THEN
  !    WRITE(2,80)'* WRITE NEGATIVE CONC. INFORMATION TO NEG-CONC.LOG'
  !  ELSE
  !    WRITE(2,80)'* NO WRTING OF NEGATIVE CONCENTRATION INFORMATION '
  !  ENDIF
  !
  !  IF( IWQRST == 1 )THEN
  !    WRITE(2,80)'* WRITE SPATIAL DISTRIBUTIONS TO IWQORST          '
  !  ELSE
  !    WRITE(2,80)'* NO WRITING TO IWQORST                           '
  !  ENDIF
  !  WRITE(2,999)
  !
  !  IF( IWQICI == 1 )THEN
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY-VARYING ICS FROM INWQICI   '
  !  ELSE IF( IWQICI == 2 )THEN
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY-VARYING ICS FROM INWQRST   '
  !  ELSE
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY CONSTANT INITIAL CONDITIONS'
  !  ENDIF
  !
  !  IF( IWQAGR == 1 )THEN
  !    WRITE(2,80)'* SPATIALLY A/O TEMPORALLY-VARYING ALGAL KINETICS '
  !  ELSE
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY CONSTANT ALGAL KINETICS    '
  !  ENDIF
  !
  !  IF( IWQSTL == 1 )THEN
  !    WRITE(2,80)'* SPATIALLY AND/OR TEMPORALLY-VARYING SETTLING VEL'
  !  ELSE
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY CONSTANT SETTLING VELOCITY '
  !  ENDIF
  !
  !  IF( IWQSUN >= 1 )THEN
  !    WRITE(2,80)'* TEMPORALLY-VARYING IO & FD                      '
  !  ELSE
  !    WRITE(2,80)'* TEMPORALLY CONSTANT IO & FD                     '
  !  ENDIF
  !
  !  IF( IWQNPL == 1 )THEN
  !    WRITE(2,80)'* SPATIALLY AND/OR TEMPORALLY-VARYING NPS INPUT   '
  !  ELSE
  !    WRITE(2,80)'* SPATIALLY/TEMPORALLY CONSTANT NPS INPUT         '
  !  ENDIF
  !
  !  IF( IWQKIN == 1 )THEN
  !    WRITE(2,80)'* SPATIALLY VARYING KINETICS FROM KINETICS.INP    '
  !  ELSE
  !    WRITE(2,80)'* FILE KINETICS.INP NOT USED                      '
  !  ENDIF
  !  WRITE(2,999)
  !
  !  WRITE(2,84) '* TIME-SERIES OUTPUT FROM ', WQTSB, ' DAY ', &
  !    '                       TO ', WQTSE, ' DAY ', &
  !    '                    EVERY ', WQTSDT, ' HOUR', &
  !    '                       AT ', IWQTS,  ' LOCATIONS', &
  !    ' BIN FILE SWITCH ISWQAVG =', ISWQAVG,' (0=OFF)  ', &
  !    ' BIN FILE SWITCH ISWQMIN =', ISWQMIN,' (0=OFF)  ', &
  !    ' BIN FILE SWITCH ISWQMAX =', ISWQMAX,' (0=OFF)  ', &
  !    ' BIN FILE SWITCH ISCOMP  =', ISCOMP, ' (0=OFF)  '
  !  WRITE(2,999)
  !
  !  IF( IWQTS >= 1 )THEN
  !    WRITE(2,80)': ICWQTS(I)=1, TIME-SERIES OUTPUT FOR VARIABLE I  '
  !    WRITE(2,80)': ICWQTS(I)\=1, NO TIME-SERIES OUTPUT FOR VAR. I  '
  !  endif
  !999 FORMAT(1X)
  !90 FORMAT(A79)
  !91 FORMAT(10I8)
  !92 FORMAT(10F8.4)
  !93 FORMAT(I8,3F8.4)
  !94 FORMAT(2I5, 13I5, /, 10X, 9I5)
  !95 FORMAT(A254)
  !80 FORMAT(A50)
  !81 FORMAT(A27, 4(F8.4,2X))
  !82 FORMAT((A45, F10.4))
  !83 FORMAT(A47, I10)
  !84 FORMAT(3(A26,F10.4,A5,/), 5(A26,I8,A10,/))
  !86 FORMAT(' I,J,M = ',3I10)
  !
  !
  !  IF( ISTRWQ(22)  > 0 )THEN
  !    !AA! Added read of CO2 half-saturation consts
  !    WRITE(2,80)'* HALF-SAT. CONSTANT (G/M^3) FOR NUTRIENT UPTAKE '
  !    WRITE(2,81)' : (KHNC, KHPC, KHCO2C)       = ', WQKHNC,WQKHPC,WQKHCO2C
  !    WRITE(2,81)' : (KHND, KHPD, KHS, KHCO2D)  = ', WQKHND,WQKHPD,WQKHS,WQKHCO2D
  !    WRITE(2,81)' : (KHNG, KHPG, KHCO2G)       = ', WQKHND,WQKHPG,WQKHCO2G
  !    WRITE(2,81)' : (KHNM, KHPM, KHCO2M)       = ', WQKHNM,WQKHPM,WQKHCO2M
  !  ELSE
  !    WRITE(2,80)'* HALF-SAT. CONSTANT (G/M^3) FOR NUTRIENT UPTAKE '
  !    WRITE(2,81)' : (KHNC, KHPC)       = ', WQKHNC,WQKHPC
  !    WRITE(2,81)' : (KHND, KHPD, KHS)  = ', WQKHND,WQKHPD,WQKHS
  !    WRITE(2,81)' : (KHNG, KHPG)       = ', WQKHND,WQKHPG
  !    WRITE(2,81)' : (KHNM, KHPM)       = ', WQKHNM,WQKHPM
  !  ENDIF
  !  WRITE(2,82)'* SAL. WHERE MICROSYSTIS GROWTH IS HALVED  = ', WQSTOX
  !
  !  WRITE(2,80)'* LIGHT EXTINC. COEFF. DUE TO TSS, CHL & POM      '
  !  WRITE(2,81)' : KETSS (/M PER G/M^3)  = ', WQKETSS
  !  WRITE(2,81)' : KECHL (/M PER MG/M^3) = ', WQKECHL
  !  WRITE(2,81)' : KECHLE (CHL EXPONENT) = ', WQKECHLE
  !  IF( WQKECHL  < 0.0 )THEN
  !    WRITE(2,80) '* USE RILEY (1956) EQUATION FOR WQKECHL          '
  !    WRITE(2,80) ' : KECHL = 0.054*CHL**0.667 + 0.0088*CHL         '
  !  ENDIF
  !  WRITE(2,81)' : KEPOM (/M PER G/M^3)  = ', WQKEPOC
  !  WRITE(2,80)'* CARBON-TO-CHL RATIO (G C PER MG CHL)            '
  !  WRITE(2,81)' : (CCHLC, CCHLD, CCHLG) = ', WQCHLC,WQCHLD,WQCHLG
  !  WRITE(2,80)'* DEPTH (M) OF MAXIMUM ALGAL GROWTH               '
  !  WRITE(2,81)' : (DOPTC, DOPTD, DOPTG) = ', WQDOPC,WQDOPD,WQDOPG
  !
  !  WRITE(2,82)'*INITIAL IO (LY/D) AT WATER SURFACE       = ',WQI0 &
  !    ,' MINIMUM OPTIMUM SOLAR RADIATION (LY/D)   = ',WQISMIN &
  !    ,' FRACTIONAL DAYLENGTH                     = ',WQFD &
  !    ,' WEIGHTING FACTOR FOR RAD. AT CURRENT DAY = ',WQCIA &
  !    ,' WEIGHTING FACTOR FOR RAD. AT (-1) DAY    = ',WQCIB &
  !    ,' WEIGHTING FACTOR FOR RAD. AT (-2) DAYS   = ',WQCIC &
  !    ,' FRACTION OF SOLAR RADIATION THAT IS PAR  = ',PARADJ
  !
  !  WRITE(2,80)'* LOWER OPTIMUM TEMP FOR ALGAL GROWTH (DEGC)     '
  !  WRITE(2,81)' : (TMC1, TMD1, TMG1   ) = ', WQTMC1,WQTMD1,WQTMG1
  !  WRITE(2,80)'* UPPER OPTIMUM TEMP FOR ALGAL GROWTH (DEGC)     '
  !  WRITE(2,81)' : (TMC2, TMD2, TMG2   ) = ', WQTMC2,WQTMD2,WQTMG2
  !
  !  WRITE(2,80)'* REFERENCE TEMPERATURE FOR ALGAL METABOLISM (OC) '
  !  WRITE(2,81)' : (TRC, TRD, TRG)       = ', WQTRC,WQTRD,WQTRG
  !  WRITE(2,80)'* TEMPERATURE EFFECT FOR ALGAL METABOLISM         '
  !  WRITE(2,81)' : (KTBC, KTBD, KTBG)    = ', WQKTBC,WQKTBD,WQKTBG
  !
  !  WRITE(2,80)'* CARBON DISTRIBUTION COEFF FOR ALGAL PREDATION   '
  !  WRITE(2,81)' : (FCRP, FCLP, FCDP)    = ', WQFCRP,WQFCLP,WQFCDP
  !  WRITE(2,80)'* CARBON DISTRIBUTION COEFF FOR ALGAL METABOLISM  '
  !  WRITE(2,81)' : (FCDC, FCDD, FCDG)    = ', WQFCDC,WQFCDD,WQFCDG
  !  WRITE(2,80)'* HALF-SAT. CONSTANT (GO/M*3) FOR ALGAL DOC EXCRET'
  !  WRITE(2,81)' : (KHRC, KHRD, KHRG)    = ', WQKHRC,WQKHRD,WQKHRG
  !
  !  WRITE(2,80)'* MINIMUM DISSOLUTION RATE (/DAY) OF ORGANIC C    '
  !  WRITE(2,81)' : (KRC, KLC, KDC)       = ', WQKRC,WQKLC,WQKDC(1)
  !  WRITE(2,80)'* CONSTANT RELATING DISSOLUTION RATE TO ALGAE     '
  !  WRITE(2,81)' : (KRCALG,KLCALG,KDCALG)= ', WQKRCALG,WQKLCALG,WQKDCALG
  !
  !  WRITE(2,80)'* REFERENCE TEMP FOR HYDROLYSIS/MINERALIZATION(OC)'
  !  WRITE(2,81)' : (TRHDR, TRMNL)        = ', WQTRHDR,WQTRMNL
  !  WRITE(2,80)'* TEMPERATURE EFFECT ON HYDROLYSIS/MINERALIZATION '
  !  WRITE(2,81)' : (KTHDR, KTMNL)        = ', WQKTHDR,WQKTMNL
  !  WRITE(2,80)'* HALF-SAT. CONSTANT FOR OXIC RESP/DENITRIFICATION'
  !  WRITE(2,81)' : (KHORDO, KHDNN)       = ', WQKHORDO,WQKHDNN
  !  WRITE(2,80)'* RATION OF DENITRIFICATION TO OXIC DOC RESP      '
  !  WRITE(2,81)' : (AANOX)               = ', WQAANOX
  !
  !  WRITE(2,80)'* PHOSPHORUS DISTRIBUTION COEF FOR ALGAL PREDATION'
  !  WRITE(2,81)' : (FPRP,FPLP,FPDP,FPIP) = ', WQFPRP,WQFPLP,WQFPDP,WQFPIP
  !  WRITE(2,80)'* PHOSPHORUS DIST COEF OF RPOP FOR ALGAL METABOLIS'
  !  WRITE(2,81)' : (FPRC, FPRD, FPRG)    = ', WQFPRC,WQFPRD,WQFPRG
  !  WRITE(2,80)'* PHOSPHORUS DIST COEF OF LPOP FOR ALGAL METABOLIS'
  !  WRITE(2,81)' : (FPLC, FPLD, FPLG)    = ', WQFPLC,WQFPLD,WQFPLG
  !
  !
  !  IF( IWQSRP /= 1 .AND. IWQSRP /= 2 )THEN
  !    WQKPO4P = 0.0
  !    WRITE(2,80)': NO SORPTION OF PO4T/SA, SO KPO4P IS FORCED TO 0 '
  !  ENDIF
  !  WRITE(2,80)'* PHOSPHORUS DIST COEF OF DOP FOR ALGAL METABOLISM'
  !  WRITE(2,81)' : (FPDC, FPDD, FPDG)    = ', WQFPDC,WQFPDD,WQFPDG
  !  WRITE(2,80)'* PHOSPHORUS DIST COEF OF NH4 FOR ALGAL METABOLISM'
  !  WRITE(2,81)' : (FPIC, FPID, FPIG)    = ', WQFPIC,WQFPID,WQFPIG
  !  WRITE(2,82)'* PARITITION COEFF FOR SORBED/DISSOLVED PO4 =', WQKPO4P
  !
  !  WRITE(2,80)'* MINIMUM HYDROLYSIS RATE (/DAY) OF ORGANIC P     '
  !  WRITE(2,81)' : (KRP, KLP, KDP)       = ', WQKRP,WQKLP,WQKDP
  !  WRITE(2,80)'* CONSTANT RELATING HYDROLYSIS RATE TO ALGAE      '
  !  WRITE(2,81)' : (KRPALG,KLPALG,KDPALG)= ', WQKRPALG,WQKLPALG,WQKDPALG
  !  WRITE(2,80)'* CONSTANT USED IN DETERMINING P-TO-C RATIO       '
  !  WRITE(2,81)' : (CPPRM1,CPPRM2,CPPRM3)= ', WQCP1PRM,WQCP2PRM,WQCP3PRM
  !
  !  WRITE(2,80)'* NITROGEN DISTRIBUTION COEFF FOR ALGAL PREDATION '
  !  WRITE(2,81)' : (FNRP,FNLP,FNDP,FNIP) = ', WQFNRP,WQFNLP,WQFNDP,WQFNIP
  !  WRITE(2,80)'* NITROGEN DIST COEF OF RPON FOR ALGAL METABOLISM '
  !  WRITE(2,81)' : (FNRC, FNRD, FNRG)    = ', WQFNRC,WQFNRD,WQFNRG
  !  WRITE(2,80)'* NITROGEN DIST COEF OF LPON FOR ALGAL METABOLISM '
  !  WRITE(2,81)' : (FNLC, FNLD, FNLG)    = ', WQFNLC,WQFNLD,WQFNLG
  !
  !  WRITE(2,80)'* NITROGEN DIST COEF OF DON FOR ALGAL METABOLISM  '
  !  WRITE(2,81)' : (FNDC, FNDD, FNDG)    = ', WQFNDC,WQFNDD,WQFNDG
  !  WRITE(2,80)'* NITROGEN DIST COEF OF NH4 FOR ALGAL METABOLISM  '
  !  WRITE(2,81)' : (FNIC, FNID, FNIG)    = ', WQFNIC,WQFNID,WQFNIG
  !  WRITE(2,80)'* NITROGEN-TO-CARBON RATIO IN ALGAE               '
  !  WRITE(2,81)' : (ANCC, ANCD, ANCG)    = ', WQANCC,WQANCD,WQANCG
  !
  !  WRITE(2,82)'* MASS NO3 REDUCED PER DOC OXIDIZED (GN/GC)= ',WQANDC &
  !    ,'  MAXIMUM NITRIFICATION RATE (/DAY)        = ',WQNITM &
  !    ,'  REFERENCE TEMP FOR NITRIFICATION (DEGC)  = ',WQTNIT
  !  WRITE(2,80)'* NITRIFICATION HALF-SAT CONSTANT FOR DO & NH4    '
  !  WRITE(2,81)' : (KHNITDO, KHNITN)     = ', WQKHNDO,WQKHNN
  !  WRITE(2,80)'* SUB & SUPER-OPTIMUM TEMP EFFECT ON NITRIFICATION'
  !  WRITE(2,81)' : (KNIT1, KNIT2)        = ', WQKN1,WQKN2
  !
  !  WRITE(2,80)'* MINIMUM HYDROLYSIS RATE (/DAY) OF ORGANIC N     '
  !  WRITE(2,81)' : (KRN, KLN, KDN)       = ', WQKRN,WQKLN,WQKDN
  !  WRITE(2,80)'* CONSTANT RELATING HYDROLYSIS RATE TO ALGAE      '
  !  WRITE(2,81)' : (KRNALG,KLNALG,KDNALG)= ', WQKRNALG,WQKLNALG,WQKDNALG
  !
  !  IF( IWQSRP /= 1 .AND. IWQSRP /= 2 )THEN
  !    WQKSAP = 0.0
  !    WRITE(2,80)': NO SORPTION OF PO4T/SA, SO KSAP IS FORCED TO 0  '
  !  ENDIF
  !  WRITE(2,80)'* SILICA DISTRIBUTION COEFF FOR DIATOM PREDATION  '
  !  WRITE(2,81)' : (FSPP, FSIP)          = ', WQFSPP,WQFSIP
  !  WRITE(2,80)'* SILICA DISTRIBUTION COEFF FOR DIATOM METABOLISM '
  !  WRITE(2,81)' : (FSPD, FSID)          = ', WQFSPD,WQFSID
  !  WRITE(2,82)'*SILICA-TO-CARBON RATIO IN DIATOMS        = ',WQASCD &
  !    ,'*PARITITION COEFF FOR SORBED/DISSOLVED SA = ',WQKSAP &
  !    ,'*DISSOLUTION RATE (/D) OF PSI             = ',WQKSU &
  !    ,' REFERENCE TEMP FOR PSI DISSOLUTION (OC)  = ',WQTRSUA &
  !    ,' TEMPERATURE EFFECT ON PSI DISSOLUTION    = ',WQKTSUA
  !
  !  WRITE(2,82)'* DO-TO-CARBON RATIO IN RESPIRATION        = ',WQAOCR &
  !    ,':MASS DO CONSUMED PER NH4-N NITRIFIED     = ',WQAONT &
  !    ,':PROPORN. CONSTANT FOR DO-REAERATION (MKS)= ',WQKRO(1) &
  !    ,' TEMPERATURE EFFECT ON DO-REAERATION      = ',WQKTR(1) &
  !    ,'*HALF-SAT CONSTANT OF DO FOR COD (GO2/M^3)= ',WQKHCOD(1) &
  !    ,':OXIDATION RATE OF COD (/DAY)             = ',WQKCD(1) &
  !    ,'  REFERENCE TEMP FOR COD OXIDATION (OC)    = ',WQTRCOD &
  !    ,'  TEMPERATURE EFFECT ON COD OXIDATION      = ',WQKTCOD &
  !    ,': DO-TO-CARBON RATIO MACROALGAE PHOTOSYNTH = ',WQAOCRPM &
  !    ,': DO-TO-CARBON RATIO MACROALGAE RESPIRATION= ',WQAOCRRM
  !
  !  WRITE(2,82) &
  !    '* DO WHERE TAM RELEASE IS HALF ANOXIC ONE  = ',WQKHBMF &
  !    ,'  ANOXIC RELEASE OF TAM (MOL/M^2/D)        = ',WQBFTAM &
  !    ,'  REFERENCE TEMP FOR TAM RELEASE (OC)      = ',WQTTAM &
  !    ,'  TEMPERATURE EFFECT ON TAM RELEASE        = ',WQKTAM &
  !    ,': TAM SOLUBILITY AT ANOXIC COND. (MOL/M^3) = ',WQTAMDMX &
  !    ,'  CONSTANT RELATING TAM SOLUBILITY TO DO   = ',WQKDOTAM &
  !    ,'* FIRST-ORDER DIE-OFF RATE AT 20OC (/D)    = ',WQKFCB &
  !    ,'  TEMPERATURE EFFECT ON BACTERIA DIE-OFF   = ',WQTFCB
  !
  !  !WRITE(2,23)'* # OF OPEN BOUNDARY CELLS ON SOUTH BOUNDARY       = ',NWQOBS
  !  !WRITE(2,23)'* # OF OPEN BOUNDARY CELLS ON WEST BOUNDARY        = ',NWQOBW
  !  !WRITE(2,23)'* # OF OPEN BOUNDARY CELLS ON EAST BOUNDARY        = ',NWQOBE
  !  !WRITE(2,23)'* # OF OPEN BOUNDARY CELLS ON NORTH BOUNDARY       = ',NWQOBN
  !
  !  IF( NWQOBS > NBBSM) CALL STOPP('ERROR!! NWQOBS SHOULD <= NBBSM')
  !  IF( NWQOBW > NBBWM) CALL STOPP('ERROR!! NWQOBW SHOULD <= NBBWM')
  !  IF( NWQOBE > NBBEM) CALL STOPP('ERROR!! NWQOBE SHOULD <= NBBEM')
  !  IF( NWQOBN > NBBNM) CALL STOPP('ERROR!! NWQOBN SHOULD <= NBBNM')
  !
  !  WRITE(2,999)
  !  WRITE(2,80)'* CONSTANT OBC AT (ICBX(M),JCBX(M)) IF IWQOBX(M)=0'
  !  WRITE(2,80)': READ TIME-SERIES OBCS IWQOBX TIMES IF IWQOBX > 0'
  !
  !  !  *** C44
  !  ! SPATIALLY/TEMPORALLY CONSTANT INITIAL CONDITIONS: WQCHLX=1/WQCHLX
  !  ! READ DATA POINTS & DO INTERNAL INTERPOLATION?
  !
  !  !IF(IWQICI == 0 )THEN    PMC - ALLOW INITIALIZATION, EVEN IF NOT USING CONSTANT IC.  THESE WILL BE OVERWRITTEN LATER, IF NEEDED
  !  WRITE(2,999)
  !  !WRITE(2,90) TITLE(1)
  !  !WRITE(2,21)' : (BC, BD, BG)         = ', (WQV(1,1,NW),NW=1,3)
  !  !WRITE(2,21)' : (RPOC, LPOC, DOC)    = ', (WQV(1,1,NW),NW=4,6)
  !  !WRITE(2,21)' : (RPOP,LPOP,DOP,PO4T) = ', (WQV(1,1,NW),NW=7,10)
  !  !WRITE(2,21)' : (RPON, LPON, DON)    = ', (WQV(1,1,NW),NW=11,13)
  !  !WRITE(2,21)' : (NH4, NO3)           = ', (WQV(1,1,NW),NW=14,15)
  !  !WRITE(2,21)' : (SU, SA, COD, DO)    = ', (WQV(1,1,NW),NW=16,19)
  !  !WRITE(2,981)' : (TAM, FCB, CO2, MAC)= ', (WQV(1,1,NW),NW=20,NWQV),WQV(1,1,IDNOTRVA)    !VB! CHANGED NWQV+1 TO NWQV   !VB added CO2
  !
  !  DO NAL = 1,NALAGE
  !    IF( ALGAES(NAL).ISMOBILE == 0 )THEN
  !    ENDIF
  !  ENDDO
  !  
  !  DO NAL = 1,NALAGE
  !    IF( ALGAES(NAL).ISMOBILE == 0 )THEN
  !    WRITE(2,9003)
  !9003 FORMAT(/,' MACALGMP.INP - MACROALGAE MAP FILE',/, &
  !      ' PSHADE = SHADE FACTOR FOR TREE CANOPY (1.0=NO CANOPY)',/, &
  !      ' KMV    = MACROALGAE HALF-SATURATION VELOCITY (M/SEC)',/, &
  !      ' KMVMIN = MACROALGAE VELOCITY LIMITATION MINIMUM (M/SEC)',/, &
  !      ' KBP    = MACROALGAE HALF-SATURATION DENSITY (G C/M2)',/, &
  !      ' KMVA   = MACROALGAE VEL. LIMIT LOGISTIC FUNC. PARAM. A',/, &
  !      ' KMVB   = MACROALGAE VEL. LIMIT LOGISTIC FUNC. PARAM. B',/, &
  !      ' KMVC   = MACROALGAE VEL. LIMIT LOGISTIC FUNC. PARAM. C',/, &
  !      ' KMVD   = MACROALGAE VEL. LIMIT LOGISTIC FUNC. PARAM. D',/, &
  !      ' KMVE   = MACROALGAE VEL. LIMIT LOGISTIC FUNC. PARAM. E',/, &
  !      '   I   J   L PSHADE    KMV KMVMIN    KBP   KMVA   KMVB', &
  !      '   KMVC   KMVD   KMVE')
  !
  !    ! *** Loop over the IJ map for each cell
  !    WRITE(2,9004) II, JJ, LL, PSHADE(LL), WQKMV(LL), WQKMVMIN(LL), WQKBP(LL), WQKMVA(LL), WQKMVB(LL), WQKMVC(LL), WQKMVD(LL), WQKMVE(LL)
  !9004 FORMAT(' ',I4,' ',I4,' ',I5, 9F7.3)
  !    endif
  !  ENDDO
  !
  !  WRITE(2,80)'* ALGAL GROWTH RATE (/DAY)                        '
  !  !WRITE(2,21)' : (PMC, PMD, PMG)       = ', WQPMC(1),WQPMD(1), WQPMG(1)
  !  WRITE(2,80)'* ALGAL BASAL METABOLISM RATE (/DAY)              '
  !  !WRITE(2,21)' : (BMRC, BMRD, BMRG)    = ', WQBMRC(1),WQBMRD(1), WQBMRG(1)
  !  WRITE(2,80)'* ALGAL PREDATION RATE (/DAY)                     '
  !  !WRITE(2,21)' : (PRRC, PRRD, PRRG)    = ', WQPRRC(1),WQPRRD(1), WQPRRG(1)
  !
  !  WRITE(2,82) '* BASE LIGHT EXTINCTION COEFFICIENT (/M)   = ',WQKEB(1)
  !
  !  IF( IWQSTL /= 1 )THEN
  !    WRITE(2,80)'* ALGAL SETTLING RATE (M/DAY)                     '
  !    !WRITE(2,21)' : (WSC, WSD, WSG)       = ', WQWSC(1),WQWSD(1),WQWSG(1)
  !    WRITE(2,80)'* POM SETTLING RATE (M/DAY)                       '
  !    !WRITE(2,21)' : (WSRP, WSLP)          = ', WQWSRP(1),WQWSLP(1)
  !    WRITE(2,80)'* SETTLING RATE OF PARTICULATE METAL (M/DAY)      '
  !    !WRITE(2,21)' : (WSS)                 = ', WQWSS(1)
  !  ENDIF
  !
  !  ! *** C47 SPATIALLY/TEMPORALLY CONSTANT BENTHIC FLUXES
  !  IF( IWQBEN == 0 )THEN
  !    WRITE(2,999)
  !    WRITE(2,90) TITLE(1)
  !    !WRITE(2,21)' : (PO4D, NH4, NO3)     = ',WQBFPO4D(1),WQBFNH4(1),WQBFNO3(1)
  !    !WRITE(2,21)' : (SAD, COD, DO)       = ',WQBFSAD(1),WQBFCOD(1),WQBFO2(1)
  !  ENDIF
  !
  !  !WRITE(2,23)'* NUMBER OF CELLS FOR POINT SOURCE INPUT  = ',IWQPS
  !  !WRITE(2,23)'* NUMBER WITH VARIABLE POINT SOURCE INPUT = ',NPSTMSR
  !  IF( IWQPS > NWQPS) CALL STOPP('ERROR!! IWQPS SHOULD BE <= NWQPS')
  !
  !  IF( IWQNPL /= 1 )THEN
  !    WRITE(2,999)
  !    WRITE(2,90) TITLE(1)
  !    !WRITE(2,21)' : (DSQ, CHC, CHD, CHG)  = ',XDSQ,(XDSL(NW),NW=1,3)
  !    !WRITE(2,21)' : (ROC, LOC, DOC)       = ',(XDSL(NW),NW=4,6)
  !    !WRITE(2,21)' : (ROP, LOP, DOP, P4D)  = ',(XDSL(NW),NW=7,10)
  !    !WRITE(2,21)' : (RON, LON, DON)       = ',(XDSL(NW),NW=11,13)
  !    !WRITE(2,21)' : (NHX, NOX)            = ',(XDSL(NW),NW=14,15)
  !    !WRITE(2,21)' : (SUU, SAA, COD, DOX)  = ',(XDSL(NW),NW=16,19)
  !    !WRITE(2,981)' : (TAM, FCB, CO2)  = ',(XDSL(NW),NW=20,NWQV)    !V!B ADDED co2
  !
  !    ! *** CONVERT FROM M3/S TO M3/DAY
  !    WQTT = XDSQ*CONV2
  !
  !    ! *** C50 WET DEPOSTION (MULTIPLIED BY RAINFALL VOLUME IN WQWET)
  !    !WRITE(2, 21)' : (CHC, CHD, CHG)       = ',(WQATM(NW),NW=1,3)
  !    !WRITE(2, 21)' : (ROC, LOC, DOC)       = ',(WQATM(NW),NW=4,6)
  !    !WRITE(2, 21)' : (ROP, LOP, DOP, P4D)  = ',(WQATM(NW),NW=7,10)
  !    !WRITE(2, 21)' : (RON, LON, DON)       = ',(WQATM(NW),NW=11,13)
  !    !WRITE(2, 21)' : (NHX, NOX)            = ',(WQATM(NW),NW=14,15)
  !    !WRITE(2, 21)' : (SUU, SAA, COD, DOX)  = ',(WQATM(NW),NW=16,19)
  !    !WRITE(2,981)' : (TAM, FCB, CO2)       = ',(WQATM(NW),NW=20,NWQV)  !!VB ADDED CO2
  !endif
  !
  !  !READ(1,295) RSTOFN
  !  !WRITE(2,85)'* OUTPUT FILE FOR RESTART WRITING         = ', RSTOFN
  !  IF( IWQRST == 1 )THEN
  !  ELSE
  !    !IF( RSTOFN(1:4) /= 'NONE' .AND. RSTOFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQORST/RSTOFN')
  !  ENDIF
  !
  !  !READ(1,295) ICIFN
  !  !WRITE(2,85)'* FILE FOR INITIAL CONDITIONS             = ', ICIFN
  !  IF( IWQICI == 1 )THEN
  !  ELSE IF( IWQICI == 2 )THEN
  !  ELSE
  !    !IF( ICIFN(1:4) /= 'NONE' .AND. ICIFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQICI/ICIFN')
  !  ENDIF
  !
  !  !READ(1,295) AGRFN
  !  !WRITE(2,85)'* FILE FOR ALGAL GROWTH, RESP., PREDATAT. = ', AGRFN
  !  IF( IWQAGR == 1 )THEN
  !  ELSE
  !    !IF( AGRFN(1:4) /= 'NONE' .AND. AGRFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQAGR/AGRFN')
  !  ENDIF
  !
  !  !READ(1,295) STLFN
  ! ! WRITE(2,85)'* FILE FOR SETTLING RATES OF ALGAE, PART. = ', STLFN
  !  IF( IWQSTL == 1 )THEN
  !  ELSE
  !    !IF( STLFN(1:4) /= 'NONE' .AND. STLFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQSTL/STLFN')
  !  ENDIF
  !
  !  !READ(1,295) SUNFN
  !  !WRITE(2,85)'* FILE FOR IO, FD, TE, KT                 = ', SUNFN
  !  IF( IWQSUN == 1 )THEN
  !  ELSE
  !    !   IF( SUNFN(1:4) /= 'NONE' .AND. SUNFN(1:4) /= 'none')
  !    !&     CALL STOPP('ERROR!! INVALID IWQSUN/SUNFN')
  !  ENDIF
  !
  !  !READ(1,295) BENFN
  !  !WRITE(2,85)'* FILE FOR BENTHIC FLUX                   = ', BENFN
  !  IF( IWQBEN == 2 )THEN
  !  ELSE
  !    !IF( BENFN(1:4) /= 'NONE' .AND. BENFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQBEN/BENFN')
  !  ENDIF
  !
  !  !READ(1,295) PSLFN
  !  !WRITE(2,85)'* FILE FOR POINT SOURCE INPUT             = ', PSLFN
  !
  !  !READ(1,295) NPLFN
  !  !WRITE(2,85)'* FILE FOR NPS INPUT INCLUDING ATM. INPUT = ', NPLFN
  !  IF( IWQNPL == 1 )THEN
  !  ELSE
  !    !IF( NPLFN(1:4) /= 'NONE' .AND. NPLFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQNPL/NPLFN')
  !  ENDIF
  !
  !  !READ(1,295) NCOFN
  !  !WRITE(2,85)'* DIAGNOSTIC FILE FOR NEGATIVE CONCENTRAT = ', NCOFN
  !  CLOSE(1)
  !
  !  IF( IWQZ  > 1 .AND. IWQKIN  > 0 )THEN
  !    !WRITE(*,'(A)')' WQ: KINETICS.INP'
  !    DO I=1,IWQZ
  !      WRITE(2,9112) IZ, IWQKA(IZ), WQKRO(IZ), WQKTR(IZ), REAC(IZ),WQKDC(IZ),WQKDCALM(IZ),WQKHRM(IZ),WQDOPM(IZ),WQKCD(IZ),WQKHCOD(IZ)
  !    ENDDO
  !  ENDIF
  !9111 FORMAT(/,'ZONE IWQKA   KRO   KTR  REAC   KDC KDCALGM  KHRM DOPTM   KCD KHCOD')
  !9112 FORMAT(I4, I6, 4F6.3, F8.3, 4F6.3)
  !
  !  IF( IWQBEN == 2 )THEN
  !    WRITE(*,'(A)')' WQ: WQBENMAP.INP'
  !    WRITE(2,999)
  !    !WRITE(2,34) L, I, J, XBENMUD(L), IBENMAP(L,1), IBENMAP(L,2)
  !  ENDIF


  END SUBROUTINE WQ_QC

  REAL(RKD) FUNCTION DO_SAT(L)
  
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: L
    REAL(RKD) :: ELEV, TVAL, RLNSAT1, RLNSAT2, RLNSAT3

    ! Elevation adjustment factor for D.O. saturation
    RLNSAT3 = 1.
    IF( IDOSELE > 0 )THEN
      ELEV = BELV(L) + HP(L) + DOELEV
      IF(IDOSELE == 1) THEN
        ! *** Chapra (1997) pg. 3
        ELEV = 1.0e-3 * ELEV                      ! Convert elevation to km
        RLNSAT3 = (((-1.60747e-4 * ELEV + 6.10834e-3)*ELEV - 0.11988)*ELEV + 1.)
      ELSEIF( IDOSELE == 2 )THEN
        ! *** Zison et al. (1978)
        RLNSAT3 = (1.0 - 0.1148e-3 * ELEV)
      ENDIF
    ENDIF

    ! THE D.O. SATURATION CALCULATION
    IF( IDOSFRM == 2 )THEN
      ! *** Genet et al. (1974)
      RLNSAT1 = ((+0.0054258 * TWQ(L) -0.38217)*TWQ(L) + 14.5532)
      DO_SAT = RLNSAT1 * RLNSAT3
    ELSEIF( IDOSFRM == 1 )THEN
        ! *** Chapra (1997) pg. 3
      TVAL = 1./(TWQ(L) + 273.15)
      RLNSAT1 = ((((-8.621949e11 * TVAL + 1.2438e10)*TVAL -6.642308e7)*TVAL + 1.575701e5)*TVAL -139.34411)
      RLNSAT2 = -SWQ(L)*((-2140.7 * TVAL + 10.754)*TVAL + 1.7674e-2)
      DO_SAT = RLNSAT3 * EXP(RLNSAT1 + RLNSAT2)
    ELSE
      ! *** DO Saturation, Modified by SCJ, see Garcia and Gordon, Limnology and Oceanography 37(6), 1992, Eqn. 8 and Table 1
      TVAL = LOG((298.15 - TWQ(L))/(273.15 + TWQ(L)))
      RLNSAT1 = (((((1.41575*TVAL + 1.01567)*TVAL + 4.93845)*TVAL + 4.11890)*TVAL + 3.20684)*TVAL + 5.80818)
      RLNSAT2 = -SWQ(L)*( 1.32412e-7*SWQ(L) + (((5.54491e-3*TVAL + 7.93334E-3)*TVAL + 7.25958e-3)*TVAL +7.01211e-3) )
      DO_SAT = 32.0e-3 * RLNSAT3 * EXP(RLNSAT1 + RLNSAT2)         ! *** 32E-3 approximately converts micromol/L to mg/L or g/m^3
    ENDIF
  END FUNCTION DO_SAT

END MODULE WATERQUALITY
