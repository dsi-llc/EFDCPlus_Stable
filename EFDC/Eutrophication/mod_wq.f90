! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
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

  use GLOBAL    
  use INFOMOD
  use JULIANMOD
  
  use Variables_MPI
  use Broadcast_Routines
  use Variables_MPI_Mapping
  use Mod_Map_Global_to_Local
  use Variables_MPI_Write_Out
  
  use Variables_WQ
  use WQ_DIAGENESIS
  use WQ_RPEM_MODULE
  use SHELLFISHMOD
  use WQ_ZOOPLANKTON
  use WQ_BIOTA

  ! *** LB is a globally declared integer that can be hardwired for debugging (defined in AAEFDC)
  
  contains
  
  
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

  use CALCSERMOD,only: CALCSER

  implicit none

  real(RKD), STATIC :: DAYNEXT
  real(RKD), STATIC :: SUNDAY1, SUNDAY2
  real,      STATIC :: SUNSOL1, SUNSOL2
  real,      STATIC :: SUNFRC1, SUNFRC2
  real,      STATIC :: WQKCNT = 0.
  
  real       :: TIMTMP, RATIO, SOLARAVG, WTEMP, WQTT, TT20, HP2I, WQ2
  integer    :: IACTION ! not used
  integer    :: IWQTAGR, IWQTSTL, ISMTICI ! not used
  integer    :: M1, M2, L, K, NMALG, NW
  integer    :: LF, LL, LP, NAL, ND
  integer, STATIC :: M
  
  real(RKD), external :: DSTIME 
  real(RKD)           :: TTDS       ! MODEL TIMING TEMPORARY VARIABLE

  DATA IWQTAGR,IWQTSTL,ISMTICI/3*0/  
  
  ! *** SET THE HYDRODYNAMIC TIMESTEP
  if( ISDYNSTP == 0 )then  
    DELT = DT  
  else  
    DELT = DTDYN  
  endif  
  DAYINT = DBLE(INT(TIMEDAY))
  
  ! *** INITIALIZE parameterS ON FIRST CALL
  if( ITNWQ == 0 )then
    ! *** INITIALIZE DAYNEXT VALUE
    DAYNEXT = DBLE(INT(TIMEDAY)) + 1.

    ! *** PMC - NEW IMPLEMENTATION TO USE DAILY (FROM HOURLY) SOLAR RADIATION FOR ALGAL GROWTH
    if( NASER > 0 )then
      if( IWQSUN == 3 )then
        ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
        SUNDAY1 = DAYNEXT - 1.
        SUNDAY2 = DAYNEXT

        ! *** FIND 1ST POINT
        M = 1
        do while (TSATM(1).TIM(M) < SUNDAY1)
          M = M+1
        enddo
    
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL1 = 0.0
        do while (TSATM(1).TIM(M) < SUNDAY1)
          M1 = M1+1
          if( TSATM(1).VAL(M,6) > 0. )then
            M2 = M2+1
            SUNSOL1 = SUNSOL1+TSATM(1).VAL(M,6)
          endif
          M = M+1
        enddo
        if( M1 > 0 )then
          SUNFRC1 = FLOAT(M2)/FLOAT(M1)
          SUNSOL1 = SUNSOL1/FLOAT(M1)
        else
          SUNFRC1 = 1.0
        endif
    
        ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
        M1 = 0
        M2 = 0
        SUNSOL2 = 0.
        do while (TSATM(1).TIM(M) < SUNDAY2)
          M1 = M1+1
          if( TSATM(1).VAL(M,6) > 0. )then
            M2 = M2+1
            SUNSOL2 = SUNSOL2+TSATM(1).VAL(M,6)
          endif
          M = M+1
        enddo
        if( M1 > 0 )then
          SUNFRC2 = FLOAT(M2)/FLOAT(M1)
          SUNSOL2 = SUNSOL2/FLOAT(M1)
          if( SUNSOL1 == 0.0 )then
            SUNFRC1 = SUNFRC2
            SUNSOL1 = SUNSOL2
          endif            
        else
          SUNFRC2 = 1.
        endif
      endif
    endif    ! *** NASER > 0
  endif      ! *** End of initialization

  ! *** WQI1 = SOLAR RADIATION ON PREVIOUS DAY  
  ! *** WQI2 = SOLAR RADIATION TWO DAYS AGO  
  ! *** WQI3 = SOLAR RADIATION THREE DAYS AGO  
  ! *** Update occurs only when the simulation day changes.  
  if( TIMEDAY > DAYNEXT )then  
    WQI3 = WQI2  
    WQI2 = WQI1  
    
    ! *** Check if no solar radiation the previous day. Only update WQI1 if SUNFRC2 > 0
    if( IWQSUN == 2 )then
      if( SUNFRC2 > 0.0 )then
        WQI1 = WQI0OPT/SUNFRC2  
      endif
      WQI0OPT = 0.0            ! *** Reset daily average
      SUNFRC2 = 0.             ! *** Reset sun fraction
    else
      WQI1 = WQI0
    endif
    DAYNEXT = DAYNEXT + 1.
  endif

  ! *** INCREMENT THE WATER QUALITY TIMESTEP
  WQKCNT = WQKCNT + DTWQ
  
  ! *** APPLY WATER DEPTH UPDATE TO WQ CONSTITUENTS THAT DO NOT TRANSPORT, I.E. MACROPHYTES
  if( NFIXED > 0 )then
    do NAL = 1, NALGAE
      if( .not. ALGAES(NAL).ISMOBILE )then
        NW = 19 + NAL
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,HP2I,WQ2) 
        do ND = 1,NDM  
          LF = (ND-1)*LDMWET+1  
          LL = min(LF+LDMWET-1,LAWET)
        
          do LP = LF,LL  
            L = LWET(LP)
            K = KSZ(L)
            HP2I = H1PK(L,K)
            if( ISTL == 3 ) HP2I = H2PK(L,K)
            WQ2 = WQV(L,K,NW)*HP2I              ! *** Mass from previous timestep
            WQ2 = WQ2*HPKI(L,K)                 ! *** New concentration
            WQV(L,K,NW) = WQ2
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
    enddo
  endif
   
  ! *** UPDATE WATER COLUMN KINETICS AND SEDIMENT MODEL  
  ! *** OVER LONGER TIME INTERVALS THAN PHYSICAL TRANSPORT  
  if( ITNWQ == 0 .or. WQKCNT >= WQKINUPT )then
    DTWQ   = WQKCNT/86400.         ! *** KINETIC TIME STEP USED (DAYS)
    DTWQO2 = DTWQ*0.5  
    WQKCNT = 0.
    
    ! *** INTIALIZE THE TIMER FOR WQ KINETICS, ALL COMPONENTS
    TTDS = DSTIME(0) 

    ! *** GET SOLAR RADIATION INTENSITY AND DAYLIGHT LENGTH  
    ! *** NOTE: IWQSUN = 1 callS SUBROUTINE WQSUN WHICH READS THE DAILY  
    ! ***                  SOLAR RADIATION DATA FROM FILE SUNDAY.INP WHICH  
    ! ***                  ARE IN UNITS OF LANGLEYS/DAY.  
    ! ***       IWQSUN = 2 USES THE HOURLY SOLAR RADIATION DATA FROM ASER.INP  
    ! ***                  COUPLED WITH THE COMPUTED OPTIMAL DAILY LIGHT TO
    ! ***                  LIMIT ALGAL GROWTH.
    ! ***       IWQSUN = 3 USES THE DAILY AVERAGE SOLAR RADIATION DATA COMPUTED 
    ! ***                  FROM THE HOURLY ASER.INP AND THE COMPUTED OPTIMAL DAILY
    ! ***                  LIGHT TO LIMIT ALGAL GROWTH.
    ! ***                  DAYLIGHT AND ADJUSTS FOR PHOTOSYNTHETIC ACTIVE RADIATION BY 
    ! ***                  PARADJ (~0.43) 

    if( IWQSUN == 0 )then
      WQI1 = WQI0         ! *** Constant
      
    elseif( IWQSUN == 1 )then  
      ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
      call WQSUN             ! *** Read SUNDAY.INP
      WQI0 = SOLSRDT  
      WQFD = SOLFRDT  
      
    elseif( NASER > 0 )then
      ! *** SOLAR RADIAION COMES FROM ASER FILE.  IWQSUN: 2-USE TIMING FROM ASER, 3-DAILY AVERAGE COMPUTED FROM ASER
      if( IWQSUN == 3 )then
        ! *** Check if new day
        if( TIMEDAY > SUNDAY2 )then
          ! *** BUILD THE DAILY AVERAGE SOLAR RADIATION FROM THE ASER DATA
          SUNDAY1 = SUNDAY2
          SUNSOL1 = SUNSOL2
          SUNFRC1 = SUNFRC2
        
          ! *** BUILD THE AVERAGE DAILY SOLAR RADIATION        
          M1 = 0
          M2 = 0
          SUNSOL2 = 0.
          SUNDAY2 = SUNDAY2 + 1.
          do while (TSATM(1).TIM(M) < SUNDAY2)
            M1 = M1+1
            if( TSATM(1).VAL(M,6) > 0. )then
              M2 = M2+1
              SUNSOL2 = SUNSOL2+TSATM(1).VAL(M,6)
            endif
            M = M+1
            if( M > TSATM(1).NREC )then
              M = TSATM(1).NREC
              exit
            endif
          enddo
          if( M1 > 0 )then
            SUNFRC2 = FLOAT(M2)/FLOAT(M1)
            SUNSOL2 = SUNSOL2/FLOAT(M1)
          else
            SUNFRC2 = 1.
          endif
        endif

        RATIO = (TIMEDAY-SUNDAY1)
        SOLARAVG = RATIO*(SUNSOL2-SUNSOL1)+SUNSOL1

        ! *** SOLAR RADIATION IN LANGLEYS/DAY
        WQI0 = PARADJ*2.065*SOLARAVG  
        WQFD = RATIO*(SUNFRC2-SUNFRC1)+SUNFRC1
    
      elseif( IWQSUN == 2 )then
        ! *** OPTIMAL SOLAR RADIATION IS ALWAYS UPDATED BASED ON DAY AVERAGED
        if( LDAYLIGHT .and. (NASER > 1 .or. USESHADE) )then  
          SOLARAVG = 0.                           ! *** SOLARAVG is the domain average solar radiation
          do L = 2,LA  
            SOLARAVG = SOLARAVG + SOLSWRT(L)      ! *** SOLSWRT already include surface albedo and shading
          enddo  
          SOLARAVG = SOLARAVG/FLOAT(LA-1)
        else
          ! *** Spatially Constant Atmospheric Parameters
          SOLARAVG = SOLSWRT(2)
        endif  
        ! *** SOLAR RADIATION IN LANGLEYS/DAY
        WQI0 = PARADJ*2.065*SOLARAVG               ! *** Current light for growth
        WQFD = 1.                                  ! *** Set fraction of day to 1 for 
        
        WQI0OPT = WQI0OPT + SOLARAVG*DTWQ          ! *** Sum current light for daily average
        if( SOLARAVG > 0.0 ) SUNFRC2 = SUNFRC2 + DTWQ
      endif
    endif  

    ! *** MASS LOADING BC'S.  WQWPSL IS ONLY USED IN THE KINETIC ROUTINES.  
    ! *** IF CONCENTRATION BASED LOADING, CALCSER AND CALFQC ALREADY HANDLED LOADING
    if( IWQPSL == 1 ) CALL WQPSL  

    ! *** CALL SPATIALLY AND TIME VARYING BENTHIC FLUX HERE, IF REQUIRED  
    ! *** IF SIMULATION TIME IS >= THE NEXT TIME IN THE BENTHIC FILE.  
    if( IWQBEN  ==  2 )then  
      TIMTMP = TIMEDAY
      
      if( TIMTMP  >=  BENDAY )then  
        call WQBENTHIC(TIMTMP)  
      endif  
    endif  

    ! *** UPDATE WET DEPOSITION (LRAIN ENSURES THAT THERE IS VALID RAINFALL)
    if( LRAIN ) CALL WQWET  
    
    do NW = 1, NWQV
      if( ISKINETICS(NW) > 0 )then 
        do K = 1,KC  
          do L = 2,LA  
            WQVO(L,K,NW) = WQV(L,K,NW)
          enddo  
        enddo
      endif  
    enddo  

    ! *** SET UP LOOK-UP TABLE FOR BACTERIA (FCB) TEMPERATURE DEPENDENCY OVER -10 degC TO 50 degC  
    if( ISKINETICS(IFCB) > 0 )then
      WTEMP = WQTDMIN
      do M1 = 1,NWQTD  
        TT20 = WTEMP - 20.0  
        WQTT = WQKFCB * WQTFCB**TT20 * DTWQO2  
        WQTD1FCB(M1) = 1.0 - WQTT  
        WQTD2FCB(M1) = 1.0 / (1.0 + WQTT)  
        WTEMP = WTEMP + WQTDINC
      enddo  
    endif
    
    
    ! ***   CALCULATE KINETIC SOURCES AND SINKS  
    if( ISWQLVL == 0 ) CALL WQSKE0  
    if( ISWQLVL == 1 ) CALL WQSKE1 ! *** Extension of CEQUAL-ICM for unlimited algae + zooplankton
    if( ISWQLVL == 2 ) CALL WQSKE2  
    if( ISWQLVL == 3 ) CALL WQSKE3  
    if( ISWQLVL == 4 ) CALL WQSKE4  
    
    ! *** ZOOPLANKTON
    ! *** It is called after WQSKE1 since the kinetic zones map is defined in WQSKE1
    if( IWQZPL > 0 )then
      call ZOOPL_KINETIC
    endif
    
    TWQKIN = TWQKIN + (DSTIME(0)-TTDS) 

    ! ***   DIAGNOSE NEGATIVE CONCENTRATIONS  
    if( IWQNC > 0 ) CALL WWQNC  

    ! *** CALL SEDIMENT DIAGENSIS MODEL AFTER THE FIRST WQ ITERATION
    if( IWQBEN == 1 .and. ITNWQ > 0 )then
      TTDS = DSTIME(0) 
      call SMMBE  
      TWQSED = TWQSED + (DSTIME(0)-TTDS) 
    endif  

    ! *** RPEM
    if( ISRPEM > 0 .and. ITNWQ > 0  )then
      TTDS = DSTIME(0) 
      call CAL_RPEM
      TWQRPEM = TWQRPEM + (DSTIME(0)-TTDS) 
    endif

  endif    ! *** ENDIF ON KINETIC AND SEDIMENT UPDATE  
  
  ! *** UPDATE WATER QUALITY TIMESTEP
  ITNWQ = ITNWQ + 1

  return
  
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
  
  use RESTART_MODULE, only:WQ_WCRST_IN, WQSDRST_IN
  implicit none

  integer :: J, L, LG, K, IWQTICI, IWQTAGR, IWQTSTL, IWQTSUN, IWQTBEN, IWQTPSL
  integer :: IWQTNPL, ISMTICI, NWQVOUT, NS, NW, NAL
  real    :: O2WQ_, WQTAMD 
  real    :: RKCWQ         !< 

  character*3 CWQHDR(NWQVM)
  DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
  DATA ISMTICI/0/

  IWQTICI = IWQTICI
  IWQTAGR = IWQTAGR
  IWQTSTL = IWQTSTL
  IWQTSUN = IWQTSUN
  IWQTBEN = IWQTBEN
  IWQTPSL = IWQTPSL
  IWQTNPL = IWQTNPL
  ISMTICI = ISMTICI
  
  if( process_id == master_id )then
    open(2,FILE = OUTDIR//'WQ3D.OUT',STATUS = 'UNKNOWN')
    close(2,STATUS = 'DELETE')
    open(2,FILE = OUTDIR//'WQ3D.OUT',STATUS = 'UNKNOWN')
  endif
  
  WQKINUPT = DT/86400.  ! WQ VARIABLE DT
  UHEQ(1) = 0.0
  UHEQ(LC) = 0.0
  do L = 2,LA
    UHEQ(L) = 1.0
  enddo
  
  ITNWQ = 0
  ! *** Fractional of Depth at the Top of the Layer
  RKCWQ = 1.0/REAL(KC)
  do K = 1,KC
    WQHT(K) = REAL(KC-K)*RKCWQ
  enddo

  do K = 1,KC
    IWQPSC(1,K) = 0
    IWQPSC(LC,K) = 0
  enddo
  do K = 1,KC
    do L = 2,LA
      IWQPSC(L,K) = 0
      IWQPSV(L,K) = 0
    enddo
  enddo
  do J = 1,NWQV
    do K = 1,KC
      do L = 1,LC
        WQWDSL(L,K,J) = 0.0
        WQWPSL(L,K,J) = 0.0
      enddo
    enddo
  enddo
  
  ! *** Read main water quality control file
  call WQ3DCONTROL
  call ALGAECONTROL

  NWQVOUT = 0
  if( process_id == master_id )then
    open(1,FILE = OUTDIR//'WQWCTS.OUT',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'WQWCTS.OUT',STATUS = 'UNKNOWN')
    do NW = 1,NWQV
      if( ISKINETICS(NW) == 1 )then
        NWQVOUT = NWQVOUT+1
        CWQHDR(NWQVOUT) = WQCONSTIT(NW)
      endif
    enddo
    write(1,1969)(CWQHDR(NW),NW = 1,NWQVOUT)
     1969 FORMAT('C   I    J    K    TIME',7X,A3,8X,A3,8X,A3, &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3, &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3, &
        8X,A3,8X,A3,8X,A3,8X,A3,8X,A3,8X,A3)
    close(1)
  endif
  
  ! ***  INITIALIZE DIURNAL DO ANALYSIS
  if( WQHRAVG > 0.0 )then
    if( process_id == master_id )then
      open(1,FILE = OUTDIR//'DIURNDO.OUT')
      close(1,STATUS = 'DELETE')
    endif
    do K = 1,KC
      do L = 2,LA
        DDOMAX(L,K) = -1.E6
        DDOMIN(L,K) = 1.E6
      enddo
    enddo
  endif

  ! ***  INITIALIZE LIGHT EXTINCTION ANALYSIS
  NDLTCNT = 0
  if( NDLTAVG >= 1 )then
    if( process_id == master_id )then
      open(1,FILE = OUTDIR//'LIGHT.OUT')
      close(1,STATUS = 'DELETE')
    endif
    do K = 1,KC
      do L = 2,LA
        RLIGHTT(L,K) = 0.
        RLIGHTC(L,K) = 0.
      enddo
    enddo
  endif

  ! *** READ INITIAL CONDITIONS
  if( IWQICI > 0 )then
    if( process_id == master_id )then
      if( IWQICI == 1 )then
        call WQICI
      else
        call WQ_WCRST_IN
      endif
    endif
    
    if( IWQICI == 2 )then
      call Broadcast_Scalar(WQI1,    master_id)
      call Broadcast_Scalar(WQI2,    master_id)
      call Broadcast_Scalar(WQI3,    master_id)
    endif
    call Broadcast_Array(WQV_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        do NW = 1,NWQV
          do K = 1,KC
            WQV(L,K,NW) = WQV_Global(LG,K,NW)
          enddo
        enddo
      endif
    enddo
  endif
  
  ! *** WQCHLX = 1/WQCHLX
  do L = 2,LA
    do K = 1,KC
      WQCHL(L,K) = 0
      do NAL = 1, NALGAE
        if( ALGAES(NAL).ISMOBILE )then
          WQCHL(L,K) = WQCHL(L,K) + WQV(L,K,19+NAL)*ALGAES(NAL).WQCHLA
        endif
      enddo
      if( IWQSRP == 1 )then
        O2WQ_ = max(WQV(L,K,IDOX), 0.0)
        WQTAMD = min( WQTAMDMX*EXP(-WQKDOTAM*O2WQ_), WQV(L,K,ITAM) )
        WQTAMP(L,K) = WQV(L,K,ITAM) - WQTAMD
        WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*WQTAMP(L,K))
        WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*WQTAMP(L,K))
      elseif( IWQSRP == 2 .and. ISTRAN(6) > 0 )then
        WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*SEDT(L,K))
        WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*SEDT(L,K))
      else
        WQPO4D(L,K) = WQV(L,K,IP4D)
        WQSAD(L,K)  = WQV(L,K,ISAA)
      endif
    enddo
  enddo
  
  ! *** Set vegetative growth and drag 
  do NAL = 1, NALGAE
    if( .not. ALGAES(NAL).ISMOBILE )then
      ! *** Set base layer for all cells
      do L = 2,LA
        ! *** Set base layer               delme - TODO - handle growth downward from suspended base
        do K = KSZ(L),KC
          if( HP(L)*Z(L,K-1) >= (HP(L)-ALGAES(NAL).BASEDEPTH) )then
            LAYERBOT(NAL,L) = K                ! *** Bottom active layer
            exit                               ! *** Jump out of the layer loop
          endif
        enddo
      enddo
            
      if( ALGAES(NAL).THRESHOLD /= 0 )then
        call Macro_Veg(NAL)
      else
        ! *** Set base layer for all cells
        do L = 2,LA
          if( SUM(WQV(L,:,19 + NAL)) >= 0.001 )then
          ! *** Use specified heights to set vegetation height
            HEIGHT_MAC(L,NAL) = ALGAES(NAL).MAXLENGTH
            
            ! *** Set top layer               delme - TODO - handle growth downward from suspended base
            do K = KSZ(L),KC
              if( HP(L)*Z(L,K-1) >= (HP(L) - ALGAES(NAL).BASEDEPTH + HEIGHT_MAC(L,NAL)) )then
                LAYERTOP(NAL,L) = K                ! *** Top active layer
                exit                               ! *** Jump out of the layer loop
              endif
            enddo
          else
              HEIGHT_MAC(L,NAL) = 0.
          endif
        enddo
      endif
    endif
  enddo
  
  ! *** ZOOPLANKTON
  if(IWQZPL > 0 )then
    call ZOOPL_CONTROL
  endif
  
  ! *** SEDIMENT DIAGENESIS
  if( IWQBEN == 1 )then
    ! *** Both Sediment Diagenesis initializations use R3D_Global in place of WQV
    allocate(R3D_Global(LCM_Global,KCM,NWQVM))
    R3D_Global = 0
    
    do L = 2,LA
      SMHYST(L) = .FALSE.
    enddo

    ! *** Read main sediment diagenesis control file  
    call SMRIN1_JNP
    !CALL SMRIN1
    
    if( process_id == master_id )then
      ! *** READ SEDIMENT MODEL INITIAL CONDITION  
      if( ISMICI == 1 )then
        call WQSDICI
      endif
      
      if( ISMICI == 2 )then
        call WQSDRST_IN
      endif
    endif
    
    call Broadcast_Array(SMPON_Global,  master_id)
    call Broadcast_Array(SMPOP_Global,  master_id)
    call Broadcast_Array(SMPOC_Global,  master_id)
  
    call Broadcast_Array(SM1NH4_Global, master_id)
    call Broadcast_Array(SM2NH4_Global, master_id)
    call Broadcast_Array(SM2NO3_Global, master_id)
    call Broadcast_Array(SM2PO4_Global, master_id)
    call Broadcast_Array(SM2H2S_Global, master_id)
    call Broadcast_Array(SMPSI_Global,  master_id)
    call Broadcast_Array(SM2SI_Global,  master_id)
    call Broadcast_Array(SMBST_Global,  master_id)
    call Broadcast_Array(SMT_Global,    master_id)
    
    if( ISMICI == 1 .or. ISMICI == 2 )then
      ! *** Map to Local Domain
      do LG = 2,LA_GLOBAL
        L = Map2Local(LG).LL
        if( L > 1 )then
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
        endif
      enddo
    endif
  
  endif

  ! *** RPEM
  if( ISRPEM > 0 )then
    call INIT_RPEMVARS
    call RPEMINP_JNP
  endif
  
  ! *** READ WQ TIMESERIES
  if( IWQPSL /= 0 )then        ! *** SCJ skip reading WQCSR if there are only constant point source loads
    if( process_id == master_id )then
      if( ISWQLVL == 0 )then
        call WQCSR
      else 
        call WQCSR2
      endif
    endif

    call Broadcast_Array(MCSER,   master_id)
    call Broadcast_Array(TCCSER,  master_id)

    do NW = 1,NWQV
      if( NWQCSR(NW) >= 1 )then
        do NS = 1,NWQCSR(NW)
          call Broadcast_Array(TSWQ(NS,NW).TIM,   master_id)
          call Broadcast_Array(TSWQ(NS,NW).VAL,   master_id)
        enddo
      endif
    enddo
  endif
  
  close(2)

  return  
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
  !    SPATIALLY AND TEMPORALLY CONSTANT REAL parameterS
  !---------------------------------------------------------------------------!
  SUBROUTINE WQ3DCONTROL
  
  use Allocate_Initialize      
  use Variables_MPI
  use Broadcast_Routines
  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  use INFOMOD,only:SKIPCOM,READSTR
  
  implicit none

  type(fson_value), Pointer :: json_data, item, cell_group, pointsource, zones, phytogroups
  Character(len = 79), allocatable :: TITLE_(:)
  Character(len = 1) :: CCMRM
  Character(len = 2) :: SNUM
  Character*80 STR*200
  integer :: I, IBIO, IZ, J, K, L, M, NW, NT, LL, N1, IM, KK, LG, NWQV0, NAL, ICOUNT
  integer :: IWQDT, IJKC, IWQZX, ITMP,  II, JJ
  real    :: XMRM1,  XMRM2,  XMRM3, XMRM4, XMRMA, XMRMB, XMRMC, XMRMD, XMRME  ! MACROALGAE
  real    :: XPSQ, WQTT, XDSQ, O2WQ_, WQTAMD
  real    :: XWQCHL, XWQTAMP, XWQPO4D, XWQSAD, TVARWQ, WTEMP, TT20
  real    :: XMUD, IZMUD, IZSAND, WQKHRA_TEMP, WQDOP_TEMP, WQKDCALM_TEMP
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
  if( ISWQLVL == 0 )then
    WQCONSTIT(ICHC) = 'CHC'
    WQCONSTIT(ICHD) = 'CHD'
    WQCONSTIT(ICHG) = 'CHG'
  else
  ! *** Algae name if WQSKE2
    do NW = 1,NALGAE
      if( NW < 10 )then
        write(SNUM,'(I1)') NW
        WQCONSTIT(19+NW) = 'ALG'//SNUM
      else
        write(SNUM,'(I2)') NW
        WQCONSTIT(19+NW) = 'ALG'//SNUM
      endif
    enddo
  endif
  
  ! *** Zooplankton name
  if( IWQZPL > 0 )then
    do NW = 1,NZOOPL
      if( NW < 10 )then
        write(SNUM,'(I1)') NW
        WQCONSTIT(19+NALGAE+NW) = 'ZOO'//SNUM
      else
        write(SNUM,'(I2)') NW
        WQCONSTIT(19+NALGAE+NW) = 'ZOO'//SNUM
      endif
    enddo
  endif
  
  if( .not. allocated(XDSL) )then
    allocate(XDSL(NWQVM))
    allocate(XPSL(NWQVM))
    XDSL = 0.0
    XPSL = 0.0
  endif
  
  if( process_id == master_id )then
    json_data => fson_parse("wq_3dwc.jnp")
    Write(*,'(A)') 'WQ: READING WQ_3DWC.JNP - MAIN WATER QUALITY CONTROL FILE'
    Write(2,'(/,A)') 'READING WQ_3DWC.JNP - MAIN WATER QUALITY CONTROL FILE'
         
    call fson_get(json_data, "title",                                           TITLE_)
    call fson_get(json_data, "kinetics_option",                                 ISWQLVL)
    call fson_get(json_data, "number_of_variables",                             NWQV)
    call fson_get(json_data, "use_kinetic_zones",                               IWQZONES)
    call fson_get(json_data, "number_of_kinetic_zones",                         NWQZ)
    call fson_get(json_data, "temperature_lookup_table_size",                   NWQTD)
    call fson_get(json_data, "number_of_time_series_output_locations",          NWQTS)
    call fson_get(json_data, "number_of_time_series_output_variables",          NTSWQV)
    call fson_get(json_data, "number_of_sediment_zones",                        NSMZ)

    call fson_get(json_data, "max_number_of_time_series_output_locations",      NSMTS)
    call fson_get(json_data, "kinetic_update_time_step",                        WQKINUPT)
    
    If( ISWQLVL < 0 .or. ISWQLVL > 4 )Call STOPP('BAD KINETICS OPTION')
    DTD = DT/86400.0
    if( IS2TIM == 0 .and. WQKINUPT < 2.*DT )then
      PRINT *,' WARNING:  Kinetic time step too small.  Setting WQKINUPT = 2*DT.'
      WQKINUPT = 2.*DT
    endif
    
    ! *** DTWQ is used later to initialize some variables in sediment diagenesis
    ! *** DTWQO2 is used later to set up look-up table for temperature dependency of FCB
    DTWQ = WQKINUPT/86400.
    DTWQO2 = DTWQ*0.5
    
    call fson_get(json_data, "active_constituents.ROC",   ISKINETICS(IROC))
    call fson_get(json_data, "active_constituents.LOC",   ISKINETICS(ILOC))
    call fson_get(json_data, "active_constituents.DOC",   ISKINETICS(IDOC))
    call fson_get(json_data, "active_constituents.ROP",   ISKINETICS(IROP))
    call fson_get(json_data, "active_constituents.LOP",   ISKINETICS(ILOP))
    call fson_get(json_data, "active_constituents.DOP",   ISKINETICS(IDOP))
    call fson_get(json_data, "active_constituents.P4D",   ISKINETICS(IP4D))
    call fson_get(json_data, "active_constituents.RON",   ISKINETICS(IRON))
    call fson_get(json_data, "active_constituents.LON",   ISKINETICS(ILON))
    call fson_get(json_data, "active_constituents.DON",   ISKINETICS(IDON))
    call fson_get(json_data, "active_constituents.NHX",   ISKINETICS(INHX))
    call fson_get(json_data, "active_constituents.NOX",   ISKINETICS(INOX))
    call fson_get(json_data, "active_constituents.SUU",   ISKINETICS(ISUU))
    call fson_get(json_data, "active_constituents.SAA",   ISKINETICS(ISAA))
    call fson_get(json_data, "active_constituents.COD",   ISKINETICS(ICOD))
    call fson_get(json_data, "active_constituents.DOX",   ISKINETICS(IDOX))
    call fson_get(json_data, "active_constituents.TAM",   ISKINETICS(ITAM))
    call fson_get(json_data, "active_constituents.FCB",   ISKINETICS(IFCB))
    call fson_get(json_data, "active_constituents.CO2",   ISKINETICS(ICO2))
    If( ISWQLVL == 0  )then
      call fson_get(json_data, "active_constituents.CHC", ISKINETICS(ICHC))
      call fson_get(json_data, "active_constituents.CHD", ISKINETICS(ICHD))
      call fson_get(json_data, "active_constituents.CHG", ISKINETICS(ICHG))
    endif

    If( ISWQLVL == 0  )then
      ISKINETICS = 0      
      ISKINETICS(IDOC) = 1
      ISKINETICS(INHX) = 1
      ISKINETICS(IDOX) = 1
    endif
       
    call fson_get(json_data, "number_of_algae_groups", NALGAE)
    if( ISWQLVL == 0 )then
      ALG_COUNT = 3
    Else
      ALG_COUNT = NALGAE
    endif
    
    ! *** Global Activation Options
    call fson_get(json_data, "silica_activate",                  IWQSI)
    call fson_get(json_data, "cyanobacteria_salinity_toxicity",  IWQSTOX)
    call fson_get(json_data, "shellfish_farm_activate",          ISFFARM)
    call fson_get(json_data, "number_of_shellfish_species",      NSF)
    call fson_get(json_data, "number_of_shellfish_cells",        NSFCELLS)
    call fson_get(json_data, "zooplankton_activate",             IWQZPL)
    call fson_get(json_data, "number_of_zooplankton_groups",     NZOOPL)
    call fson_get(json_data, "rpem_activate",                    ISRPEM)
    
    If( IWQZPL > 0 )then
      NWQVZ = NWQV - NZOOPL
    endif

    call fson_get(json_data, "po4_sorption_option",                                 IWQSRP)
    call fson_get(json_data, "log_negative_concentrations",                         IWQNC)
    ! call fson_get(json_data, "write_restart",                                       IWQRST)
    
    call fson_get(json_data, "formulation_for_DO_saturation",                       IDOSFRM)
    call fson_get(json_data, "elevation_adjustment_for_DO_saturation",              IDOSELE)
    call fson_get(json_data, "elevation_offset_for_DO_saturation",                  DOELEV)
    
    call fson_get(json_data, "number_of_hours_averaging_DO",                        WQHRAVG)
                                                                                    
    call fson_get(json_data, "initial_condition_option",                            IWQICI)
    call fson_get(json_data, "point_source_load_option",                            IWQPSL)
    call fson_get(json_data, "use_atmospheric_dry_deposition",                      IWQNPL)
    call fson_get(json_data, "solar_radiation.source_option",                       IWQSUN)
    call fson_get(json_data, "solar_radiation.initial_optimal_sr",                  WQI0)
    call fson_get(json_data, "solar_radiation.minimum_optimal_sr",                  WQISMIN)
    call fson_get(json_data, "solar_radiation.fraction_of_daylight",                WQFD)
    call fson_get(json_data, "solar_radiation.photoactive_radiation_fraction",      PARADJ)
    call fson_get(json_data, "solar_radiation.daily_weighting_factors",             WQCI)
    WQCIA = WQCI(1)
    WQCIB = WQCI(2)
    WQCIC = WQCI(3)
    WQCIM = WQCI(4)
    
    WQI0 = PARADJ*WQI0        ! *** Apply conversion to photosynthesis active light fraction
    WQI1 = WQI0
    WQI2 = WQI0
    WQI3 = WQI0
    WQI0OPT = 0.0
    
    call fson_get(json_data, "light_extinction.light_extinction_diagnostics",      NDLTAVG)
    call fson_get(json_data, "light_extinction.chlorophyll_coefficient",           WQKECHL)
    call fson_get(json_data, "light_extinction.chlorophyll_exponent",              WQKECHLE)
    call fson_get(json_data, "light_extinction.particular_organic_matter_coeff",   WQKEPOC)
    call fson_get(json_data, "light_extinction.dissolved_organic_carbon_coeff",    WQKEDOM)
    !Call fson_get(json_data, "light_extinction.background_coeff",                  WQKEB(1))
    !Call fson_get(json_data, "light_extinction.total_suspended_solids_coeff",      WQKETSS)
    
    call fson_get(json_data, "reaeration.reaeration_option",                               IWQKA(1))
    call fson_get(json_data, "reaeration.reaeration_constant",                             WQKRO(1))
    call fson_get(json_data, "reaeration.temperature_rate_const",                          WQKTR(1))
    call fson_get(json_data, "reaeration.adjustment_factor",                               REAC(1))
    call fson_get(json_data, "nutrient_sorption.partition_coeff_for_sorbed_dissolved_PO4", WQKPO4P)
    call fson_get(json_data, "nutrient_sorption.partition_coeff_for_sorbed_dissolved_SA",  WQKSAP)
    If( IWQSRP /= 1 .and. IWQSRP /= 2  )then
      WQKPO4P = 0.0
      WQKSAP = 0.0
    endif
          
    call fson_get(json_data, "hydrolysis.reference_temperature",                       WQTRHDR)
    call fson_get(json_data, "hydrolysis.effect_of_temperature",                       WQKTHDR)
    call fson_get(json_data, "hydrolysis.carbon.minimum_rate.RPOC",                    WQKRC)
    call fson_get(json_data, "hydrolysis.carbon.minimum_rate.LPOC",                    WQKLC)
    call fson_get(json_data, "hydrolysis.carbon.constant_relating_to_algae.RPOC",      WQKRCALG)
    call fson_get(json_data, "hydrolysis.carbon.constant_relating_to_algae.LPOC",      WQKLCALG)
    call fson_get(json_data, "hydrolysis.phosphorus.minimum_rate.RPOP",                WQKRP)
    call fson_get(json_data, "hydrolysis.phosphorus.minimum_rate.LPOP",                WQKLP)
    call fson_get(json_data, "hydrolysis.phosphorus.constant_relating_to_algae.RPOP",  WQKRPALG)
    call fson_get(json_data, "hydrolysis.phosphorus.constant_relating_to_algae.LPOP",  WQKLPALG)
    call fson_get(json_data, "hydrolysis.phosphorus.carbon_to_phosphorus_ratio",       WQCPPRM)
    call fson_get(json_data, "hydrolysis.nitrogen.minimum_rate.RPON",                  WQKRN)
    call fson_get(json_data, "hydrolysis.nitrogen.minimum_rate.LPON",                  WQKLN)
    call fson_get(json_data, "hydrolysis.nitrogen.constant_relating_to_algae.RPON",    WQKRNALG)
    call fson_get(json_data, "hydrolysis.nitrogen.constant_relating_to_algae.LPON",    WQKLNALG)

    WQCP1PRM = WQCPPRM(1)
    WQCP2PRM = WQCPPRM(2)
    WQCP3PRM = WQCPPRM(3)

    call fson_get(json_data, "mineralization.reference_temperature",                      WQTRMNL)
    call fson_get(json_data, "mineralization.effect_of_temperature",                      WQKTMNL)
    call fson_get(json_data, "mineralization.carbon.minimum_rate.DOC",                    WQKDC(1))
    call fson_get(json_data, "mineralization.carbon.constant_relating_to_algae.DOC",      WQKDCALG)
    call fson_get(json_data, "mineralization.carbon.constant_relating_to_macroalgae.DOC", WQKDCALM_TEMP)
    
    call fson_get(json_data, "mineralization.phosphorus.minimum_rate.DOP",                WQKDP)
    call fson_get(json_data, "mineralization.phosphorus.constant_relating_to_algae.DOP",  WQKDPALG)

    call fson_get(json_data, "mineralization.nitrogen.minimum_rate.DON",                  WQKDN)
    call fson_get(json_data, "mineralization.nitrogen.constant_relating_to_algae.DON",    WQKDNALG)
    
                                                                                       
    call fson_get(json_data, "nitrification.mass_NO3_reduces_per_DOC_oxidized",        WQANDC)
    call fson_get(json_data, "nitrification.max_rate",                                 WQNITM)
    call fson_get(json_data, "nitrification.half_sat_const_for_DO",                    WQKHNDO)
    call fson_get(json_data, "nitrification.half_sat_const_for_NH4",                   WQKHNN)
    call fson_get(json_data, "nitrification.reference_temperature",                    WQTNIT)
    call fson_get(json_data, "nitrification.suboptimal_temperature_effect_const",      WQKN1)
    call fson_get(json_data, "nitrification.superoptimal_temperature_effect_const",    WQKN2)
    
    call fson_get(json_data, "denitrification.oxic_respiration_half_sat_const_for_DO", WQKHORDO)
    call fson_get(json_data, "denitrification.half_sat_const",                         WQKHDNN)
    call fson_get(json_data, "denitrification.ratio_to_oxic_DOC_respiration",          WQAANOX)
    WQAANOX = WQAANOX*WQKHORDO
    
    call fson_get(json_data, "silica_dissolution.dissolution_rate",           WQKSU)
    call fson_get(json_data, "silica_dissolution.reference_temperature",      WQTRSUA)
    call fson_get(json_data, "silica_dissolution.effect_of_temperature",      WQKTSUA)
                                                                              
    call fson_get(json_data, "TAM_release.half_anoxic_rate_DO",               WQKHBMF)
    call fson_get(json_data, "TAM_release.anoxic_release_rate",               WQBFTAM)
    call fson_get(json_data, "TAM_release.reference_temperature",             WQTTAM)
    call fson_get(json_data, "TAM_release.effect_of_temperature",             WQKTAM)
    call fson_get(json_data, "TAM_release.solubility_at_anoxic_conditions",   WQTAMDMX)
    call fson_get(json_data, "TAM_release.solubility_to_DO_const",            WQKDOTAM)
                                                                              
    call fson_get(json_data, "coliform_decay_rate.first_order_decay_rate",    WQKFCB)
    call fson_get(json_data, "coliform_decay_rate.temperature_effect_const",  WQTFCB)
    
    call fson_get(json_data, "COD_decay.oxygen_half_sat_const_for_COD_decay", WQKHCOD(1))
    call fson_get(json_data, "COD_decay.COD_decay_rate",                      WQKCD(1))
    call fson_get(json_data, "COD_decay.reference_temperature",               WQTRCOD)
    call fson_get(json_data, "COD_decay.effect_of_temperature",               WQKTCOD)
    
    call fson_get(json_data, "settling_velocity.refractory_POM",              WQWSRP(1))
    call fson_get(json_data, "settling_velocity.labile_POM",                  WQWSLP(1))
    call fson_get(json_data, "settling_velocity.particles_sorbed_to_TAM",     WQWSS(1))

    !---------------------------------------------------------------------------
    ! C47: Spatially/Temporally constant benthic fluxes
    !---------------------------------------------------------------------------
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_option",                   IWQBEN)
    call fson_get(json_data, "sediment_diagenesis.number_of_reactive_classes",            NSMG)
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.phosphate",          WQBFPO4D(1))
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.ammonia",            WQBFNH4(1))
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.nitrate",            WQBFNO3(1))
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.silica",             WQBFSAD(1))
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.COD",                WQBFCOD(1))
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.SOD",                WQBFO2(1))
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_rates.temperature_factor", STEMFAC)     

    ! *** Initialize all zones
    Do IZ = 2, NWQZ
      WQKEB(IZ)   = SWRATNF
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
    enddo

    !---------------------------------------------------------------------------
    ! C30: Concentration time series data for open boundaries
    ! Number of time series for each state variables of nutrient
    !---------------------------------------------------------------------------
    do NW = 1,NWQV
      call fson_get(json_data,'number_of_time_series.'//WQCONSTIT(NW), NWQCSR(NW))
    enddo
    
    NT = 0
    do NW = 1,NWQV
      NT = max(NT,NWQCSR(NW))
    enddo
    NCSER(8) = NT
    
    !---------------------------------------------------------------------------
    ! C31 + C32 + C33 + C34: South Boundary
    !---------------------------------------------------------------------------
    call fson_get(json_data, "open_boundaries.south.number_of_cells", NWQOBS)
    call fson_get(json_data, "open_boundaries.south.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      call fson_get(item, "I", IWQCBS(M))
      call fson_get(item, "J", JWQCBS(M))
      Do NW = 1, NWQV
        call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBS(M,NW))
        call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCS(M,1,NW))
        call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCS(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      if( IWQCBS(M) == ICBS(M) .and. JWQCBS(M) == JCBS(M) )then
        NCSERS(M,8) = IWQOBS(M,IDOX)   ! *** ALL CONSTITUENTS USE THE SAME SERIES
        Do NW = 1, NWQV
          if( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom
            CBS(M,1,NT) = WQOBCS(M,1,NW)
            ! *** Top
            CBS(M,2,NT) = WQOBCS(M,2,NW)
          endif
        enddo
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  call STOPP('WQ: SOUTH OBC: MISS MATCH BETWEEN NCBS & NWQOBS')
      endif
    Enddo
    
    !---------------------------------------------------------------------------  
    ! C31 + C35 + C36 + C37: West Boundary  
    !---------------------------------------------------------------------------
    call fson_get(json_data, "open_boundaries.west.number_of_cells", NWQOBW) 
    call fson_get(json_data, "open_boundaries.west.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      call fson_get(item, "I", IWQCBW(M))
      call fson_get(item, "J", JWQCBW(M))
      Do NW = 1, NWQV
        call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBW(M,NW))
        call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCW(M,1,NW))
        call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCW(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      if( IWQCBW(M) == ICBW(M) .and. JWQCBW(M) == JCBW(M) )then
        NCSERW(M,8) = IWQOBW(M,IDOX)   ! *** ALL CONSTITUENTS use THE SAME SERIES
        Do NW = 1, NWQV
          if( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom
            CBW(M,1,NT) = WQOBCW(M,1,NW)
            ! *** Top
            CBW(M,2,NT) = WQOBCW(M,2,NW)
          endif
        enddo
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  call STOPP('WQ: WST OBC: MISS MATCH BETWEEN NCBW & NWQOBW')
      endif
    enddo
    !---------------------------------------------------------------------------  
    ! C30 + C38 + C39 + C40: East Boundary  
    !---------------------------------------------------------------------------
    call fson_get(json_data, "open_boundaries.east.number_of_cells", NWQOBE) 
    call fson_get(json_data, "open_boundaries.east.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      call fson_get(item, "I", IWQCBE(M))
      call fson_get(item, "J", JWQCBE(M))
      Do NW = 1, NWQV
        call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBE(M,NW))
        call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCE(M,1,NW))
        call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCE(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      if( IWQCBE(M) == ICBE(M) .and. JWQCBE(M) == JCBE(M) )then
        NCSERE(M,8) = IWQOBE(M,IDOX)   ! *** ALL CONSTITUENTS use THE SAME SERIES
        Do NW = 1, NWQV
          if( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom  
            CBE(M,1,NT) = WQOBCE(M,1,NW)
            ! *** Top
             CBE(M,2,NT) = WQOBCE(M,2,NW)
          endif
        enddo
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  call STOPP('WQ: EAST OBC: MISS MATCH BETWEEN NCBE & NWQOBE')
      endif
    Enddo
    !---------------------------------------------------------------------------  
    ! C30 + C41 + C42 + C43: North Boundary  
    !---------------------------------------------------------------------------
    call fson_get(json_data, "open_boundaries.north.number_of_cells", NWQOBN)    
    call fson_get(json_data, "open_boundaries.north.cells",           cell_group)
    Do M = 1, fson_value_count(cell_group)
      item => fson_value_get(cell_group, M)
      call fson_get(item, "I", IWQCBN(M))
      call fson_get(item, "J", JWQCBN(M))
      Do NW = 1, NWQV
        call fson_get(item,'time_series.'//WQCONSTIT(NW), IWQOBN(M,NW))
        call fson_get(item,'const_bottom_concentration.'//WQCONSTIT(NW),  WQOBCN(M,1,NW))
        call fson_get(item,'const_surface_concentration.'//WQCONSTIT(NW), WQOBCN(M,2,NW))
      Enddo
      ! *** CONCENTRATION ASSIGNMENTS
      if( IWQCBN(M) == ICBN(M) .and. JWQCBN(M) == JCBN(M) )then
        NCSERN(M,8) = IWQOBN(M,IDOX)   ! *** ALL CONSTITUENTS use THE SAME SERIES
        Do NW = 1, NWQV
          if( ISKINETICS(NW) > 0 )then
            NT = MSVWQV(NW)
            ! *** Bottom  
            CBN(M,1,NT) = WQOBCN(M,1,NW)
            ! *** Top
            CBN(M,2,NT) = WQOBCN(M,2,NW)
          endif
        enddo
      !ELSE
      ! *** This condition will be checked lately by calling Map_OpenBC_Eutrophication - DKT
      !  call STOPP('WQ: NORTH OBC: MISS MATCH BETWEEN NCBN & NWQOBN')
      endif
    Enddo
    !---------------------------------------------------------------------------
    ! C48: TEMPORALLY-CONSTANT VALUES FOR POINT SOURCE CONCENTRATIONS IN MG/L
    !---------------------------------------------------------------------------
    call fson_get(json_data, "mass_loading_point_sources.number_of_point_sources",  NWQPS) 
    call fson_get(json_data, "mass_loading_point_sources.number_of_time_series",    NPSTMSR) 
    call fson_get(json_data, "mass_loading_point_sources.constant_point_sources",   pointsource)
    ICOUNT = fson_value_count(pointsource)
    do N1 = 1, ICOUNT   !  fson_value_count(pointsource)
      item => fson_value_get(pointsource, N1)
      call fson_get(item, "I", I)
      call fson_get(item, "J", J)
      call fson_get(item, "K", K)
      call fson_get(item, "NSR", ITMP)
      call fson_get(item, "PSQ", XPSQ)
      Do NW = 1, NWQV
        call fson_get(item,'concentrations.'//WQCONSTIT(NW), XPSL(NW))
      Enddo
      
      if( IJCT_Global(I,J) < 1 .or. IJCT_Global(I,J) > 8 )then
        call STOPP('ERROR!! INVALID (I,J) IN FILE WQ_3DWC.INP FOR PSL')
      endif 
    
      ! *** HANDLE CONCENTRATION BASED POINT SOURCE
      if( IWQPSL == 2 )then
        if( LIJ_Global(I,J) /= BCPS_GL(N1).L )then
          call STOPP('MISMATCH NQSIJ BETWEEN EFDC.INP & WQ_3DWC.INP')
        endif
        ICPSL(N1) = I
        JCPSL(N1) = J

        ! *** ASSIGN GLOBAL CONCENTRATION TIME SERIES INDEX
        NSERWQ(N1) = ITMP  ! *** ALL WQ VARIABLES use SAME TIME SERIES (NSERWQ IS TEMPORARY STORAGE FOR LATER MAPPING TO LOCAL NCSERQ(N1,8))
        do NW = 1,NWQV
          if( ISKINETICS(NW) > 0 )then
            WQWPSLC(N1,NW) = XPSL(NW)
          endif
        enddo
      else            
        ICPSL(N1) = I
        JCPSL(N1) = J
        KCPSL(N1) = K
        MVPSL(N1) = ITMP
        
        ! *** save FOR NOW, CONVERT IN SUB-DOMAIN MAPPING
        WQTT = XPSQ*CONV2                   ! *** CONVERT FROM M3/S TO M3/DAY
        do NW = 1,NWQV
          WQWPSLC(N1,NW) = XPSL(NW)*WQTT
        enddo  
      endif
    enddo
    
    !---------------------------------------------------------------------------
    ! C44: Spatially/Temporally Constant Initial Concentration
    !---------------------------------------------------------------------------
    Do NW = 1,NWQV
      call fson_get(json_data,'const_initial_conditions.'//WQCONSTIT(NW), WQV(1,1,NW))
    enddo
    !---------------------------------------------------------------------------
    ! C49: Constant Dry Atmospheric Deposition(g / m2 / day; MPN / m2 / day)
    !---------------------------------------------------------------------------
    call fson_get(json_data,"dry_atmospheric_deposition.DSQ", XDSQ)
    Do NW = 1,NWQV
      call fson_get(json_data,'dry_atmospheric_deposition.'//WQCONSTIT(NW), XDSL(NW))
    enddo
    
    !---------------------------------------------------------------------------
    ! C50: Wet Atmospheric Deposition concentrations (mg / L, TAM - MOLES / L, FCB - MPN / 100ml)
    !---------------------------------------------------------------------------
    Do NW = 1,NWQV
      call fson_get(json_data,'wet_atmospheric_deposition.'//WQCONSTIT(NW), WQATM(NW,1))
    Enddo        
 
    call fson_get(json_data, "netcdf_output.water_quality.ROC",   IS_NC_OUT(25+IROC))
    call fson_get(json_data, "netcdf_output.water_quality.LOC",   IS_NC_OUT(25+ILOC))
    call fson_get(json_data, "netcdf_output.water_quality.DOC",   IS_NC_OUT(25+IDOC))
    call fson_get(json_data, "netcdf_output.water_quality.ROP",   IS_NC_OUT(25+IROP))
    call fson_get(json_data, "netcdf_output.water_quality.LOP",   IS_NC_OUT(25+ILOP))
    call fson_get(json_data, "netcdf_output.water_quality.DOP",   IS_NC_OUT(25+IDOP))
    call fson_get(json_data, "netcdf_output.water_quality.P4D",   IS_NC_OUT(25+IP4D))
    call fson_get(json_data, "netcdf_output.water_quality.RON",   IS_NC_OUT(25+IRON))
    call fson_get(json_data, "netcdf_output.water_quality.LON",   IS_NC_OUT(25+ILON))
    call fson_get(json_data, "netcdf_output.water_quality.DON",   IS_NC_OUT(25+IDON))
    call fson_get(json_data, "netcdf_output.water_quality.NHX",   IS_NC_OUT(25+INHX))
    call fson_get(json_data, "netcdf_output.water_quality.NOX",   IS_NC_OUT(25+INOX))
    call fson_get(json_data, "netcdf_output.water_quality.SUU",   IS_NC_OUT(25+ISUU))
    call fson_get(json_data, "netcdf_output.water_quality.SAA",   IS_NC_OUT(25+ISAA))
    call fson_get(json_data, "netcdf_output.water_quality.COD",   IS_NC_OUT(25+ICOD))
    call fson_get(json_data, "netcdf_output.water_quality.DOX",   IS_NC_OUT(25+IDOX))
    call fson_get(json_data, "netcdf_output.water_quality.TAM",   IS_NC_OUT(25+ITAM))
    call fson_get(json_data, "netcdf_output.water_quality.FCB",   IS_NC_OUT(25+IFCB))
    call fson_get(json_data, "netcdf_output.water_quality.CO2",   IS_NC_OUT(25+ICO2))
    
    call fson_get(json_data, "netcdf_output.phytoplankton",     IS_NC_OUT(45))
    call fson_get(json_data, "netcdf_output.zooplankton",       IS_NC_OUT(46))
    call fson_get(json_data, "netcdf_output.rpem",              IS_NC_OUT(47))
    call fson_get(json_data, "netcdf_output.shellfish",         IS_NC_OUT(48))
    call fson_get(json_data, "netcdf_output.sediment_fluxes",   IS_NC_OUT(50))
    
    ! *** DEFAULT ZONAL KINETICS,   IWQZONES = 0 OR 1
    do I = 2,NWQZ
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
      WQKEB(I)   = SWRATNF
    enddo
    
    do I = 1,NWQZ
      do NAL = 1,NALGAE
        ALGAES(NAL).WQKDCALM(I) = WQKDCALM_TEMP
      enddo
    enddo
           
    ! *** Reading the kinetics file for spatial zones WQ parameters 
    if( IWQZONES  > 0 )then
      write(*,'(A)')' WQ: WQ_KIN_ZONES.JNP'
      write(2,'(/,A)')'READING WQ_KIN_ZONES.JNP - Zonal kinetics'
      json_data => fson_parse("wq_kin_zones.jnp")
       
      call fson_get(json_data, "reaeration.reaeration_option",                  IWQKA)    ! EE WQZoneP[0]
      call fson_get(json_data, "reaeration.reaeration_constant",                WQKRO)    ! EE WQZoneP[1]
      call fson_get(json_data, "reaeration.temperature_rate_const",             WQKTR)    ! EE WQZoneP[2]
      call fson_get(json_data, "reaeration.adjustment_factor",                  REAC)     ! EE WQZoneP[3]
      call fson_get(json_data, "minimum_DOC_hydrolysis_rate",                   WQKDC)    ! EE WQZoneP[4]
      call fson_get(json_data, "COD_decay.COD_decay_rate",                      WQKCD)    ! EE WQZoneP[5]
      call fson_get(json_data, "COD_decay.oxygen_half_sat_const_for_COD_decay", WQKHCOD)  ! EE WQZoneP[6]
      call fson_get(json_data, "background_light_extinction_coeff",             WQKEB)    ! EE WQZoneP[7]
      call fson_get(json_data, "settling_velocity.RPOM",                        WQWSRP)   ! EE WQZoneP[8]
      call fson_get(json_data, "settling_velocity.LPOM",                        WQWSLP)   ! EE WQZoneP[9]
      call fson_get(json_data, "settling_velocity.TAM",                         WQWSS)    ! EE WQZoneP[10]
    endif
    
  endif  ! *** End of master_id block
  
  ! *** Broadcast water quality global settings, not dependent on grid. 

  ! *** Broadcast scalars
  call Broadcast_Scalar(ISWQLVL,     master_id)
  call Broadcast_Scalar(IWQZONES,    master_id)        
  call Broadcast_Scalar(NWQZ,        master_id)
  call Broadcast_Scalar(NWQTD,       master_id)
  call Broadcast_Scalar(NWQTS,       master_id)
  call Broadcast_Scalar(NTSWQV,      master_id)
  call Broadcast_Scalar(NSMZ,        master_id)
                                     
  call Broadcast_Scalar(NSMTS,       master_id)
  call Broadcast_Scalar(WQKINUPT,    master_id)                                   
  call Broadcast_Scalar(DTD,         master_id)

  call Broadcast_Scalar(DTWQ,        master_id)
  call Broadcast_Scalar(DTWQO2,      master_id)
                                     
  call Broadcast_Scalar(NALGAE,      master_id)
  call Broadcast_Scalar(ALG_COUNT,   master_id)
                                     
  call Broadcast_Scalar(IWQSI,       master_id)
  call Broadcast_Scalar(IWQSTOX,     master_id)
  call Broadcast_Scalar(ISFFARM,     master_id)
  call Broadcast_Scalar(NSF,         master_id)
  call Broadcast_Scalar(NSFCELLS,    master_id)
  call Broadcast_Scalar(IWQZPL,      master_id)
  call Broadcast_Scalar(NZOOPL,      master_id)
  call Broadcast_Scalar(ISRPEM,      master_id)
  call Broadcast_Scalar(NWQVZ,       master_id)
                                     
  call Broadcast_Scalar(IWQSRP,      master_id)
  call Broadcast_Scalar(IWQNC,       master_id)
  ! call Broadcast_Scalar(IWQRST,      master_id)
    
  call Broadcast_Scalar(IDOSFRM,     master_id)
  call Broadcast_Scalar(IDOSELE,     master_id)
  call Broadcast_Scalar(DOELEV,      master_id)
                                     
  call Broadcast_Scalar(WQHRAVG,     master_id)
                                     
  call Broadcast_Scalar(IWQICI,      master_id)
  call Broadcast_Scalar(IWQPSL,      master_id)
  call Broadcast_Scalar(IWQNPL,      master_id)
  call Broadcast_Scalar(IWQSUN,      master_id)
  call Broadcast_Scalar(WQI0,        master_id)
  call Broadcast_Scalar(WQISMIN,     master_id)
  call Broadcast_Scalar(WQFD,        master_id)
  call Broadcast_Scalar(PARADJ,      master_id)
                                     
  call Broadcast_Scalar(WQCIA,       master_id)
  call Broadcast_Scalar(WQCIB,       master_id)
  call Broadcast_Scalar(WQCIC,       master_id)
  call Broadcast_Scalar(WQCIM,       master_id)
  call Broadcast_Scalar(WQI1,        master_id)
  call Broadcast_Scalar(WQI2,        master_id)
  call Broadcast_Scalar(WQI3,        master_id)
  call Broadcast_Scalar(WQI0OPT,     master_id)
                                     
  call Broadcast_Scalar(NDLTAVG,     master_id)
  call Broadcast_Scalar(WQKECHL,     master_id)
  call Broadcast_Scalar(WQKECHLE,    master_id)
  call Broadcast_Scalar(WQKEPOC,     master_id)
  call Broadcast_Scalar(WQKEDOM,     master_id)
                                     
  call Broadcast_Scalar(IWQKA(1),    master_id)
  call Broadcast_Scalar(WQKRO(1),    master_id)
  call Broadcast_Scalar(WQKTR(1),    master_id)
  call Broadcast_Scalar(REAC(1),     master_id)
  call Broadcast_Scalar(WQKPO4P,     master_id)
  call Broadcast_Scalar(WQKSAP,      master_id)
                                     
  call Broadcast_Scalar(WQTRHDR,     master_id)
  call Broadcast_Scalar(WQKTHDR,     master_id)
  call Broadcast_Scalar(WQKRC,       master_id)
  call Broadcast_Scalar(WQKLC,       master_id)
  call Broadcast_Scalar(WQKDC(1),    master_id)
  call Broadcast_Scalar(WQKRCALG,    master_id)
  call Broadcast_Scalar(WQKLCALG,    master_id)
  call Broadcast_Scalar(WQKDCALG,    master_id)
  call Broadcast_Scalar(WQKRP,       master_id)
  call Broadcast_Scalar(WQKLP,       master_id)
  call Broadcast_Scalar(WQKDP,       master_id)
  call Broadcast_Scalar(WQKRPALG,    master_id)
  call Broadcast_Scalar(WQKLPALG,    master_id)
  call Broadcast_Scalar(WQKDPALG,    master_id)
                                     
  call Broadcast_Scalar(WQCP1PRM,    master_id)
  call Broadcast_Scalar(WQCP2PRM,    master_id)
  call Broadcast_Scalar(WQCP3PRM,    master_id)
                                     
  call Broadcast_Scalar(WQKRN,       master_id)
  call Broadcast_Scalar(WQKLN,       master_id)
  call Broadcast_Scalar(WQKDN,       master_id)
  call Broadcast_Scalar(WQKRNALG,    master_id)
  call Broadcast_Scalar(WQKLNALG,    master_id)
  call Broadcast_Scalar(WQKDNALG,    master_id)
                                     
  call Broadcast_Scalar(WQTRMNL,     master_id)
  call Broadcast_Scalar(WQKTMNL,     master_id)
                                     
  call Broadcast_Scalar(WQANDC,      master_id)
  call Broadcast_Scalar(WQNITM,      master_id)
  call Broadcast_Scalar(WQKHNDO,     master_id)
  call Broadcast_Scalar(WQKHNN,      master_id)
  call Broadcast_Scalar(WQTNIT,      master_id)
  call Broadcast_Scalar(WQKN1,       master_id)
  call Broadcast_Scalar(WQKN2,       master_id)
                                     
  call Broadcast_Scalar(WQKHORDO,    master_id)
  call Broadcast_Scalar(WQKHDNN,     master_id)
  call Broadcast_Scalar(WQAANOX,     master_id)
                                     
  call Broadcast_Scalar(WQKSU,       master_id)
  call Broadcast_Scalar(WQTRSUA,     master_id)
  call Broadcast_Scalar(WQKTSUA,     master_id)
                                     
  call Broadcast_Scalar(WQKHBMF,     master_id)
  call Broadcast_Scalar(WQBFTAM,     master_id)
  call Broadcast_Scalar(WQTTAM,      master_id)
  call Broadcast_Scalar(WQKTAM,      master_id)
  call Broadcast_Scalar(WQTAMDMX,    master_id)
  call Broadcast_Scalar(WQKDOTAM,    master_id) 
                                     
  call Broadcast_Scalar(WQKFCB,      master_id)
  call Broadcast_Scalar(WQTFCB,      master_id)
  
  call Broadcast_Scalar(WQKHCOD(1),  master_id)
  call Broadcast_Scalar(WQKCD(1),    master_id)
  call Broadcast_Scalar(WQTRCOD,     master_id)
  call Broadcast_Scalar(WQKTCOD,     master_id)
                                     
  call Broadcast_Scalar(WQWSRP(1),   master_id)
  call Broadcast_Scalar(WQWSLP(1),   master_id)
  call Broadcast_Scalar(WQWSS(1),    master_id)
                                     
  call Broadcast_Scalar(IWQBEN,      master_id)
  call Broadcast_Scalar(NSMG,        master_id) 
  call Broadcast_Scalar(WQBFPO4D(1), master_id)
  call Broadcast_Scalar(WQBFNH4(1),  master_id)
  call Broadcast_Scalar(WQBFNO3(1),  master_id)
  call Broadcast_Scalar(WQBFSAD(1),  master_id)
  call Broadcast_Scalar(WQBFCOD(1),  master_id)
  call Broadcast_Scalar(WQBFO2(1),   master_id)
  call Broadcast_Scalar(STEMFAC,     master_id)
  
  call Broadcast_Scalar(NCSER(8),    master_id)
                                     
  ! *** Arrays                       
  call broadcast_array(ISKINETICS,   master_id)
  call Broadcast_Array(NWQCSR,       master_id)

  ! *** Zonal Kinetics

  call Broadcast_Array(IWQKA,      master_id)
  call Broadcast_Array(WQKRO,      master_id)
  call Broadcast_Array(WQKTR,      master_id)
  call Broadcast_Array(REAC,       master_id)
  call Broadcast_Array(WQKDC,      master_id)
  call Broadcast_Array(WQKCD,      master_id)
  call Broadcast_Array(WQKHCOD,    master_id)
  call Broadcast_Array(WQKEB,      master_id)
  call Broadcast_Array(WQWSRP,     master_id)
  call Broadcast_Array(WQWSLP,     master_id)
  call Broadcast_Array(WQWSS,      master_id)
  
  ! *** Dry Deposition
  call Broadcast_Scalar(XDSQ,      master_id)
  call Broadcast_Array(XDSL,       master_id)

  ! *** Wet Deposition
  call Broadcast_Array(WQATM,      master_id)
  ! *** End of global water quality settings   
  
  ! **************************************************************************
  ! *** Open Boundaries
  call Broadcast_Scalar(NWQOBS,    master_id)
  call Broadcast_Scalar(NWQOBW,    master_id)
  call Broadcast_Scalar(NWQOBE,    master_id)
  call Broadcast_Scalar(NWQOBN,    master_id)

  call Broadcast_Array(IWQCBS,     master_id)
  call Broadcast_Array(JWQCBS,     master_id)
  call Broadcast_Array(IWQOBS,     master_id)
  call Broadcast_Array(WQOBCS,     master_id)
  call Broadcast_Array(CBS,        master_id)

  call Broadcast_Array(IWQCBW,     master_id)
  call Broadcast_Array(JWQCBW,     master_id)
  call Broadcast_Array(IWQOBW,     master_id)
  call Broadcast_Array(WQOBCW,     master_id)
  call Broadcast_Array(CBW,        master_id)

  call Broadcast_Array(IWQCBE,     master_id)
  call Broadcast_Array(JWQCBE,     master_id)
  call Broadcast_Array(IWQOBE,     master_id)
  call Broadcast_Array(WQOBCE,     master_id)
  call Broadcast_Array(CBE,        master_id)

  call Broadcast_Array(IWQCBN,     master_id)
  call Broadcast_Array(JWQCBN,     master_id)
  call Broadcast_Array(IWQOBN,     master_id)
  call Broadcast_Array(WQOBCN,     master_id)
  call Broadcast_Array(CBN,        master_id)
  
  ! *** Map to Local
  call Map_OpenBC_Eutrophication

  ! **************************************************************************
  ! *** Point Sources
  call Broadcast_Scalar(NWQPS,   master_id)
  call Broadcast_Scalar(NPSTMSR, master_id)

  call Broadcast_Array(ICPSL,    master_id)
  call Broadcast_Array(JCPSL,    master_id)
  call Broadcast_Array(KCPSL,    master_id)
  call Broadcast_Array(MVPSL,    master_id)
  call Broadcast_Array(NSERWQ,   master_id)
  call Broadcast_Array(WQWPSLC,  master_id)
  call Broadcast_Array(IWQPSC,   master_id)
  call Broadcast_Array(IWQPSV,   master_id)  
  call Broadcast_Array(CQS,      master_id)

  call Broadcast_Array(IS_NC_OUT,master_id)

  ! *** Map to Local
  call Map_WQ_PointSource
  
  ! **************************************************************************
  ! *** Set up look-up table for temperature dependency over -10 �C to 60 �C
  WQTDMIN = -10
  WQTDMAX =  60
  WTEMP = WQTDMIN
  WQTDINC = (WQTDMAX-WQTDMIN)/NWQTD
  ALLOCATE (WQTDTEMP(NWQTD))
  
  do M = 1,NWQTD
    WQTDTEMP(M) = WTEMP

    WQTDHDR(M) = EXP( WQKTHDR*(WTEMP-WQTRHDR) )   ! *** TEMPERATURE EFFECT ON HYDROLYSIS
    WQTDMNL(M) = EXP( WQKTMNL*(WTEMP-WQTRMNL) )   ! *** TEMPERATURE EFFECT ON MINERALIZATION
    
    ! *** TEMPERATURE ADJUSTED RATES FOR NITRIFICATION
    if( WTEMP > WQTNIT )then
      WQTDNIT(M) = WQNITM*EXP(-WQKN2*(WQTNIT-WTEMP)**2)
    else
      WQTDNIT(M) = WQNITM*EXP(-WQKN1*(WTEMP-WQTNIT)**2)
    endif
    
    WQKSUA(M)    = WQKSU * EXP( WQKTSUA*(WTEMP-WQTRSUA) )        ! *** TEMPERATURE EFFECT ON PSI DISSOLUTION
    WQKCOD(M,1)  = WQKCD(1) * EXP( WQKTCOD*(WTEMP-WQTRCOD) )     ! *** TEMPERATURE EFFECT ON COD OXIDATION
    
    TT20 = WTEMP - 20.0
    WQTDKR(M,1) = WQKTR(1)**TT20                                 ! *** TEMPERATURE EFFECT ON REAERATION
    
    WQTDTAM(M) = WQKHBMF * WQBFTAM * EXP( WQKTAM*(WTEMP-WQTTAM) )
    
    WQTT = WQKFCB * WQTFCB**TT20 * DTWQO2
    WQTD1FCB(M) = 1.0 - WQTT
    WQTD2FCB(M) = 1.0 / (1.0 + WQTT)
    
    WTEMP = WTEMP + WQTDINC
  enddo
  
  if( process_id == master_id )then
    write(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for Water Column Temperature Dependency",                               &
                                             "WQTDMIN = ", WQTDMIN, "WQTDMAX = ", WQTDMAX, "WQTDINC = ", WQTDINC, "NWQTD = ", NWQTD
    write(2,'(/,A5,A10,A15,20A10)') "IT", "TEMP", "WQTDMNL", "WQTDHDR", "WQTDNIT", "WQKSUA", "WQKCOD", "WQTDKR", "WQTDTAM", "WQTD1FCB", "WQTD2FCB"
                                            
    do M = 1,NWQTD
      write(2,'(I5,F10.3,F15.3,20F10.5)') M, WQTDTEMP(M), WQTDMNL(M), WQTDHDR(M), WQTDNIT(M)  , WQKSUA(M),  WQKCOD(M,1),  WQTDKR(M,1),  WQTDTAM(M), WQTD1FCB(M), WQTD2FCB(M)
    enddo
  endif
  
  ! *** Set up look-up table for temperature dependency over WQTDMIN to WQTDMAX
  do M = 1,NWQTD
    WTEMP = WQTDTEMP(M)
    TT20 = WTEMP - 20.0
    do I = 1,NWQZ
      WQKCOD(M,I) = WQKCD(I) * EXP( WQKTCOD*(WTEMP-WQTRCOD) )
      WQTDKR(M,I) = WQKTR(I)**TT20
    enddo
  enddo
  
  ! **************************************************************************
  ! *** Dry Deposition, compute mass loadings, g/day, moles/day & mpn/day
  if( IWQNPL /= 1 )then
    do L = 2,LA
      do NW = 1,NWQV
        ! M. MORTON MODIFIED THE LINE BELOW SO THAT CONSTANT ATMOSPHERIC DEPOSIT
        ! CAN BE ADDED VIA THIS ROUTINE INSTEAD OF CONSTANT NPS INPUT WHICH THE
        ! ORIGINAL CODE CALLED FOR AND WHICH WAS NOT PARTICULARLY USEFUL.
        ! INPUT DATA (XDSL) ARE IN G/M2/DAY AND ARE MULTIPLIED BY THE CELL SURFA
        ! AREA (DXYP) TO GET G/DAY.  ATMOSPHERIC DEPOSITION ONLY ENTERS THRU SUR
        ! LAYER (KC):
        WQWDSL(L,KC,NW) = XDSL(NW) * DXYP(L)
      enddo
      WQWDSL(L,KC,ITAM) = XDSL(ITAM) * CONV1
    enddo
  endif
  
  ! **************************************************************************
  ! *** Initial conditions or spatial varing parameters
  call Broadcast_Array(WQV, master_id)

  ! *** APPLY SPATIALLY CONSTAND WQ IC'S AT 1 AND LC
  do NW = 1,NWQV
    TVARWQ = WQV(1,1,NW)
    do K = 1,KC
      WQV(LC,K,NW)  = TVARWQ
      WQV(1 ,K,NW)  = TVARWQ
      WQVO(LC,K,NW) = TVARWQ
      WQVO(1 ,K,NW) = TVARWQ
    enddo
  enddo

  ! *** APPLY SPATIALLY CONSTAND WQ IC'S
  ! *** Initialize nutrient + biota (g C/m^3)
  do NW = 1,NWQV
    TVARWQ = WQV(1,1,NW)
    do K = 1,KC
      do L = 2,LA
        WQV(L,K,NW)  = TVARWQ
        WQVO(L,K,NW) = TVARWQ
      enddo
    enddo
  enddo
  
  do NAL = 1,NALGAE
    if( ISMOB(NAL) == 0 )then  ! *** ALGAES(NAL).ISMOBILE is only available after calling ALGAECONTROL - DKT
      ! *** Zero concentration if fixed biota
      do K = 1,KC
        do L = 2,LA
          WQV(L,K,19+NAL) = 0
        enddo
      enddo
      ! *** Updated fixed biota
      
      ! *** Fixed biota can be converted to biomass density by dividing by bottom layer thickness
      TVARWQ = WQV(1,1,19+NAL)  
      do L = 1,LCM
        WQV(L,KSZ(L),19+NAL)  = TVARWQ
        WQVO(L,KSZ(L),19+NAL) = TVARWQ
      enddo
    endif
  enddo
  
  ! *** Setup PO4 sorption options
  if( IWQSRP == 1 )then
    ! *** TAM sorption
    O2WQ_ = max(WQV(1,1,IDOX), 0.0)
    WQTAMD = min( WQTAMDMX*EXP(-WQKDOTAM*O2WQ_), WQV(1,1,ITAM) )
    WQTAMP(1,1) = WQV(1,1,ITAM) - WQTAMD
    WQPO4D(1,1) = WQV(1,1,IP4D) / (1.0 + WQKPO4P*WQTAMP(1,1))
    WQSAD(1,1)  = WQV(1,1,ISAA) / (1.0 + WQKSAP*WQTAMP(1,1))
  elseif( IWQSRP == 2 .and. ISTRAN(6) > 0 )then
    ! *** Cohesive sorption
    WQPO4D(1,1) = WQV(1,1,IP4D) / (1.0 + WQKPO4P*SEDT(1,1))
    WQSAD(1,1)  = WQV(1,1,ISAA) / (1.0 + WQKSAP*SEDT(1,1))
  else
    ! *** No sortpion
    WQPO4D(1,1) = WQV(1,1,IP4D)
    WQSAD(1,1)  = WQV(1,1,ISAA)
  endif
  
  XWQTAMP = WQTAMP(1,1)
  XWQPO4D = WQPO4D(1,1)
  XWQSAD  = WQSAD(1,1)
  do K = 1,KC
    WQTAMP(LC,K) = XWQTAMP
    WQPO4D(LC,K) = XWQPO4D
    WQSAD(LC,K)  = XWQSAD
    WQTAMP(1,K)  = XWQTAMP
    WQPO4D(1,K)  = XWQPO4D
    WQSAD(1,K)   = XWQSAD
  enddo

  do K = 1,KC
    do L = 2,LA
      WQTAMP(L,K) = XWQTAMP
      WQPO4D(L,K) = XWQPO4D
      WQSAD(L,K) = XWQSAD
    enddo
  enddo
  
  ! *** Initialize biota.  Units are in g C/m^3
  ! *** Fixed biota can be converted to biomass density by dividing by bottom layer thickness
  
  ! *** C47 SPATIALLY/TEMPORALLY CONSTANT BENTHIC FLUXES
  call Broadcast_Array(WQBFPO4D,   master_id)
  call Broadcast_Array(WQBFNH4,    master_id)
  call Broadcast_Array(WQBFNO3,    master_id)
  call Broadcast_Array(WQBFSAD,    master_id)
  call Broadcast_Array(WQBFCOD,    master_id)
  call Broadcast_Array(WQBFO2,     master_id)

  if( IWQBEN == 0 )then
    ! *** Initialize domain
    do L = 2,LA
      WQBFPO4D(L)= WQBFPO4D(1)
      WQBFNH4(L) = WQBFNH4(1)
      WQBFNO3(L) = WQBFNO3(1)
      WQBFSAD(L) = WQBFSAD(1)
      WQBFCOD(L) = WQBFCOD(1)
      WQBFO2(L)  = WQBFO2(1)
    enddo
  endif
   
  
  ! *** READ IN MAPPING INFORMATION FOR SPATIALLY-VARYING parameterS (UNIT #7)
  do K = 1,KC
    do L = 2,LA
      IWQZMAP(L,K) = 1
    enddo
  enddo
        
  if( NWQZ > 1 )then
    allocate(I2D_Global(LCM_Global,KCM))
    I2D_Global = 0
    
    if( process_id == master_id )then
      write(*,'(A)')' WQ: WQWCMAP.INP'
      write(2,'(/,A)')'Reading WQWCMAP.INP - Water quality zone map'
      open(1,FILE = 'wqwcmap.inp',STATUS = 'UNKNOWN')
      call SKIPCOM(1,'*',2)  ! *** SKIP OVER TITLE AND AND HEADER LINES

      IM = 0
      IJKC = IC_Global*JC_Global*KC

      do M = 1,IJKC
        read(1,*,end = 1111) I, J, K, IWQZX
        IM = IM + 1
        if( IJCT_Global(I,J) < 1 .or. IJCT_Global(I,J) > 8 )then
          PRINT*, 'I, J, K, IJCT(I,J) = ', I,J,K,IJCT_Global(I,J)
          call STOPP('ERROR!! INVALID (I,J) IN FILE WQWCMAP.INP')
        endif
        L = LIJ_Global(I,J)
        
        I2D_Global(L,K) = IWQZX
      enddo
1111  continue
      if( IM /= (LA_Global-1)*KC )then
        PRINT *, 'WARNING: ALL ACTIVE WATER CELLS SHOULD BE MAPPED FOR WQ PAR.'
        PRINT *, '         NUMBER OF LINES IN FILE WQWCMAP.INP  = \ (LA-1)'
      endif
      close(1)
    endif   ! *** End of master_id block

    call Broadcast_Array(I2D_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        do K = 1,KC
          IWQZMAP(L,K) = I2D_Global(LG,K)
        enddo
      endif
    enddo
    deallocate(I2D_Global)   

  endif
  
  ! *** Read in mapping information for spatially-varying benthic fluxes.
  ! *** formulated for Peconic Bay data which includes %mud for each cell as
  ! *** well as mapping to both mud and sand fluxes.  Subroutine WQBENTHIC
  ! *** contains the code to interpolate the final flux for the cell based
  ! *** on the percent mud and the mud/sand fluxes.
    
  if( IWQBEN == 2 )then
    call AllocateDSI(R2D_Global, LCM_Global, 3, 0.0)
    
    do K = 1,2
      do L = 2,LA
        IBENMAP(L,K) = 1
        XBENMUD(L) = 0.50
      enddo
    enddo
    if( process_id == master_id )then
      write(*,'(A)')' WQ: WQBENMAP.INP'
      write(2,'(/,A)')'Reading WQBENMAP.INP - Benthic flux rate map for user specified benthic fluxes'
      open(1,FILE = 'wqbenmap.inp',STATUS = 'UNKNOWN')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES

      REWIND(1)
      CCMRM = '#'
      call SKIPCOM(1, CCMRM,2)
      
      R2D_Global(:,1) = -999.
      IM = 0
      IJKC = IC*JC
      do M = 1,IJKC
        read(1,*,end = 1112) I, J, XMUD, IZMUD, IZSAND
        IM = IM + 1
        
        if( IJCT_Global(I,J) < 1 .or. IJCT_Global(I,J) > 8 )then
          PRINT *, 'I, J, K, IJCT(I,J) = ', I,J,IJCT_Global(I,J)
          call STOPP('ERROR!! INVALID (I,J) IN FILE WQBENMAP.INP')
        endif
        L = LIJ_Global(I,J)
        
        R2D_Global(L,1) = XMUD
        R2D_Global(L,2) = IZMUD
        R2D_Global(L,3) = IZSAND
        
      enddo
1112  continue
      if( IM  /=  (LA_Global-1) )then
        PRINT *, 'WARNING: ALL ACTIVE WATER CELLS SHOULD BE MAPPED FOR BENTHIC FLUXES.'
        PRINT *, '         NUMBER OF LINES IN FILE WQBENMAP.INP <> (LA-1)'
        PRINT *, '         LA-1     = ', LA_Global-1
        PRINT *, '         WQBENMAP = ', IM
        PRINT *, '         Details can be found in log_mpi_proc_000'
        
        
        call WriteBreak(mpi_efdc_out_unit)
        write(mpi_efdc_out_unit,'(a)') 'Warning.  The following cells have not been mapped in WQBENMAP.INP'
        do L = 2,LA_Global
          if( R2D_Global(L,1) == -999. )then
            write(mpi_efdc_out_unit,'(a,i10,2(a,i8))') 'L = ', L, '  I = ', Map2Global(L).iG,  '  J = ', Map2Global(L).jG
          endif
        enddo
        
      endif
      close(1)
    endif   ! *** End of master_id block

    call Broadcast_Array(R2D_Global,     master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        IBENMAP(L,1) = R2D_Global(LG,2)        ! *** IZMUD  - Zone to be used for mud
        IBENMAP(L,2) = R2D_Global(LG,3)        ! *** IZSAND - Zone to be used for sand
        XBENMUD(L)   = R2D_Global(LG,1)/100.0  ! *** XMUD   - Percent mud for interpolation between IZMUD and IZSAND
      endif
    enddo
    deallocate(R2D_Global)   
  endif

  ! *** These are required to initialize WQ variables along boundaries so that concentrations discontinuities can be smoothed.
  do NW = 1,NWQV
    M = MSVWQV(NW)
    if( M > 0 )then
      do K = 1,KC
        do LL = 1,NCBS
          CLOS(LL,K,M) = WQV(1,1,NW)
          NLOS(LL,K,M) = NITER
        enddo
        do LL = 1,NCBW
          CLOW(LL,K,M) = WQV(1,1,NW)
          NLOW(LL,K,M) = NITER
        enddo
        do LL = 1,NCBE
          CLOE(LL,K,M) = WQV(1,1,NW)
          NLOE(LL,K,M) = NITER
        enddo
        do LL = 1,NCBN
          CLON(LL,K,M) = WQV(1,1,NW)
          NLON(LL,K,M) = NITER
        enddo
      enddo
    endif
  enddo
    

  END SUBROUTINE WQ3DCONTROL
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQICI
  !
  !> @details  Read in spatially varying - depth averaged ICs
  !---------------------------------------------------------------------------!
  !
  !---------------------------------------------------------------------------! 
  SUBROUTINE WQICI
  
  implicit none

  integer  :: I, J, L, LG, K, KG, NW, NCOL
  real :: WQVDA(100)
  character*80 STR*200
  
  write(*,'(A)')' WQ: WQICI.INP'
  write(2,'(/,A)')'Reading WQICI.INP - Reading Spatially Variable - Depth Averaged Initial Conditions'
  open(1,FILE = "WQICI.INP",STATUS = 'UNKNOWN')

  NCOL  = 0

  ! *** Skip the first 5 lines
  do k = 1,5
    read(1,*)
  enddo

  ! *** Loop over cells and assign layers
  do LG = 2,LA_Global
    read(1,*) I, J, (WQVDA(NW),NW = 1,NWQV)
    L = LIJ(I,J)
    do NW = 1,NWQV
      do KG = KSZ(L),KC
        WQV_Global(L,KG,NW) = WQVDA(NW)
      enddo
    enddo
  enddo
    
  close(1)
  
  return
  
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

  use INFOMOD,only:SKIPCOM,READSTR

  implicit none

  character*11 FNWQSR(40)
  character*2  SNUM
  character*80 STR*200
  integer NC,NW,IS,NS,ISO,ISTYP,K,M
  real    RMULADJ,ADDADJ,CSERTMP,TOFFSET

  write(*,'(A)') ' WQ: READING WQCSRxx.INP - WQ CONCENTRATION TIME SERIES'
  write(2,'(/,A)') 'Reading WQCSRxx.INP - WQ CONCENTRATION TIME SERIES'

  ! *** DEFINE THE INPUT FILE NAMES
  do NW = 1,NWQV
    write(SNUM,'(I2.2)')NW
    FNWQSR(NW) = 'wqcsr'//SNUM//'.inp'
  enddo
  
  if( IWQZPL > 0 )then
    do NW = 1,NZOOPL
      write(SNUM,'(I2.2)')NW
      FNWQSR(22+NW) = 'zoosr'//SNUM//'.inp'
    enddo
  endif
  
  ! ***  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
  ! ***  TIME SERIES FROM THE FILES WQCSRNN.INP
  NC = 8
  do NW = 1,NWQV
    if( NWQCSR(NW) >= 1 )then
      open(1,FILE = FNWQSR(NW),STATUS = 'UNKNOWN')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      do NS = 1,NWQCSR(NW)
        MTSCLAST(NS,NC) = 2
        read(1,*,IOSTAT = ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TOFFSET, RMULADJ, ADDADJ
        if( ISO > 0 ) GOTO 900
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) GOTO 900
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSWQ(NS,NW).TIM(M),CSERTMP
            if( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TOFFSET
            do K = 1,KC
              TSWQ(NS,NW).VAL(M,K) = ( RMULADJ*(CSERTMP + ADDADJ) )*WKQ(K)
            enddo
          enddo
        else
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSWQ(NS,NW).TIM(M),(TSWQ(NS,NW).VAL(M,K), K = 1,KC)
            if( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TOFFSET
            do K = 1,KC
              TSWQ(NS,NW).VAL(M,K) = RMULADJ*( TSWQ(NS,NW).VAL(M,K) + ADDADJ )
            enddo
          enddo
        endif
      enddo
      close(1)
    endif
  enddo

  GOTO 901

  900 continue
  write(6,601) NW, NS, M
  call STOPP('.')

  901 continue
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQ TIME SERIES, NWQ,NSER,MDATA = ',3I5)
  602 FORMAT(' READ OF FILES WQCSRNN.INP SUCCESSFUL'/)

  return
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

  use INFOMOD,only:SKIPCOM,READSTR

  implicit none

  character*15 FNWQSR(40)
  character*2  SNUM
  character*80 STR*200
  integer :: NC,NW,IS,NS,ISO,ISTYP,K,M
  real ::   RMULADJ,ADDADJ,CSERTMP,TOFFSET
  
  write(*,'(A)') ' WQ: READING WQ CONCENTRATION TIME SERIES'
  
  ! *** DEFINE THE INPUT FILE NAMES
  ! *** Nutrient
  do NW = 1,19
    write(SNUM,'(I2.2)') NW
    FNWQSR(NW) = 'wqcsr'//SNUM//'.inp'
  enddo
  ! *** Algae
  do NW = 1,NALGAE
    write(SNUM,'(I2.2)')NW
    FNWQSR(NW+19) = 'wqalgsr'//SNUM//'.inp'
  enddo
  ! *** Zooplankton
  do NW = 1,NZOOPL
     write(SNUM,'(I2.2)')NW
     FNWQSR(NW+19+NALGAE) = 'wqzoosr'//SNUM//'.inp'
  enddo

  ! ***  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
  ! ***  TIME SERIES FROM THE FILES WQCSRNN.INP
  NC = 8
  do NW = 1,NWQV
    if( NWQCSR(NW) >= 1 )then
      open(1,FILE = FNWQSR(NW),STATUS = 'UNKNOWN')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      do NS = 1,NWQCSR(NW)
        MTSCLAST(NS,NC) = 2
        read(1,*,IOSTAT = ISO) ISTYP, MCSER(NS,NC), TCCSER(NS,NC), TOFFSET, RMULADJ, ADDADJ
        if( ISO > 0 ) GOTO 900
        if( ISTYP == 1 )then
          read(1,*,IOSTAT = ISO) (WKQ(K),K = 1,KC)
          if( ISO > 0 ) GOTO 900
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSWQ(NS,NW).TIM(M),CSERTMP
            if( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TOFFSET
            do K = 1,KC
              TSWQ(NS,NW).VAL(M,K) = ( RMULADJ*(CSERTMP + ADDADJ) )*WKQ(K)
            enddo
          enddo
        else
          do M = 1,MCSER(NS,NC)
            read(1,*,IOSTAT = ISO) TSWQ(NS,NW).TIM(M),(TSWQ(NS,NW).VAL(M,K), K = 1,KC)
            if( ISO > 0 ) GOTO 900
            TSWQ(NS,NW).TIM(M) = TSWQ(NS,NW).TIM(M) + TOFFSET
            do K = 1,KC
              TSWQ(NS,NW).VAL(M,K) = RMULADJ*( TSWQ(NS,NW).VAL(M,K) + ADDADJ )
            enddo
          enddo
        endif
      enddo
      close(1)
    endif
  enddo

  GOTO 901

  900 continue
  write(6,601) NW, NS, M
  call STOPP('.')

  901 continue
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQ TIME SERIES, NWQ,NSER,MDATA = ',3I5)
  602 FORMAT(' READ OF FILES WQCSRNN.INP SUCCESSFUL'/)

  return
  END SUBROUTINE WQCSR2
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSUN
  !
  !> @details READ IN TEMPORALLY VARYING parameterS FOR DAILY SOLAR RADIATION (WQI0)
  ! AND FRACTIONAL DAYLENGTH (WQFD) (UNIT INWQSUN).
  !
  !---------------------------------------------------------------------------!
  ! NEW VERSION BY J. M. HAMRICK  7 APRIL 1997
  ! READS AND INTERPOLATES DAILY AVERAGE SOLAR RADIATION AND DAYLIGHT FRACTION
  !
  !---------------------------------------------------------------------------!
  SUBROUTINE WQSUN
  
  use INFOMOD,only:SKIPCOM,READSTR

  implicit none

  integer :: IS, M, ISO, NSUNDAY, M1, M2
  
  real    :: TCSUNDAY, TASUNDAY, RMULADJ, ADDADJ, TDIFF, WTM1, WTM2, TIME
  integer :: ISPAR

  character*80 STR*200

  if( ITNWQ > 0 ) GOTO 1000
  !
  ! ***  READ IN DAILY AVERAGE SOLAR SW RAD SERIES FROM FILE 'SUNDAY.INP'
  !
  write(*,'(A)')' WQ: SUNDAY.INP'
  open(1,FILE = 'sunday.inp',STATUS = 'UNKNOWN')

  STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
  M = 0
  ISPAR = 1
  !
  !      MCSUNDAY = 1 TEMP use ISPAR FOR MCSUNDAY
  !
  read(1,*,IOSTAT = ISO)NSUNDAY,TCSUNDAY,TASUNDAY,RMULADJ,ADDADJ
  if( ISO > 0 ) GOTO 900
  do M = 1,NSUNDAY
    read(1,*,IOSTAT = ISO)TSSRD(M),SOLSRD(M),SOLFRD(M)
    if( ISO > 0 ) GOTO 900
    TSSRD(M) = TCSUNDAY*( TSSRD(M)+TASUNDAY )
    SOLSRD(M) = RMULADJ*(SOLSRD(M)+ADDADJ) * PARADJ
  enddo
  close(1)
  GOTO 901
    900 continue
  write(6,601)M
  call STOPP('.')
    901 continue
      1 FORMAT(120X)
    601 FORMAT(' READ ERROR FILE SUNDAY.INP ')
   1000 continue
  !
  ! ***  DAILY AVERAGE SOLAR SW RADIATION INTERPOLTATION FOR WATER QUALITY
  !
  TIME = TIMESEC/86400.

  M1 = ISPAR
  !
  !      TEMP use ISPAR FOR MCSUNDAY
  !
    100 continue
  M2 = M1+1
  if( TIME > TSSRD(M2) )then
    M1 = M2
    GOTO 100
  else
    ISPAR = M1
  !
  !      TEMP use ISPAR FOR MCSUNDAY
  !
  endif
  TDIFF = TSSRD(M2)-TSSRD(M1)
  WTM1 = (TSSRD(M2)-TIME)/TDIFF
  WTM2 = (TIME-TSSRD(M1))/TDIFF
  SOLSRDT = WTM1*SOLSRD(M1)+WTM2*SOLSRD(M2)
  SOLFRDT = WTM1*SOLFRD(M1)+WTM2*SOLFRD(M2)
  return  
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

  use INFOMOD,only:SKIPCOM,READSTR
  
  implicit none

  integer :: IS,NS,ISO,M,NW,M2,M1,K,L,ITMP,KK
  real    :: TAWQPSR,RMULADJ,ADDADJ,TIME,TDIFF,WTM1,WTM2
  real    :: RLDTMP(NTSWQVM)
  
  character*80 STR*200

  if( ITNWQ > 0 ) GOTO 1000

  ! ***  READ IN LOADING SERIES FROM FILE 'WQPSL.INP'
  if( NPSTMSR >= 1 )then
      
    ! *** Read only on master
    if(process_id == master_id )then

      write(*,'(A)')' WQ: WQPSL.INP'
      open(1,FILE = "WQPSL.INP",STATUS = 'UNKNOWN')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      do NS = 1,NPSTMSR
        MWQPTLT(NS) = 2 
        read(1,*,IOSTAT = ISO)NMLSER(NS),TCWQPSR(NS),TAWQPSR,RMULADJ,ADDADJ
        if( ISO > 0 ) GOTO 900

        ! *** CONVERT WQ VAR 1-19, 22 FROM KG/D TO G/D  !VB
        ! *** CONVERT WQ VAR 20 (TAM) FROM KMOLS/D TO MOLES/D
        ! *** CONVERT FECAL COLIFORM FROM MPN/DAY TO MPN/D FOR FCM
        RMULADJ = 1000.*RMULADJ
        !ADDADJ = ADDADJ

        do M = 1,NMLSER(NS)
          read(1,*,IOSTAT = ISO) T_MLSER(M,NS),(RLDTMP(NW),NW = 1,7)
          if( ISO > 0 ) GOTO 900
          read(1,*,IOSTAT = ISO)(RLDTMP(NW),NW = 8,14)
          if( ISO > 0 ) GOTO 900
          read(1,*,IOSTAT = ISO)(RLDTMP(NW),NW = 15,NWQV)  ! PMC HARDWIRED FOR TENKILLER
          if( ISO > 0 ) GOTO 900

          ! *** STANDARD CONVERSIONS
          T_MLSER(M,NS) = T_MLSER(M,NS) + TAWQPSR
          do NW = 1,NWQV
            V_MLSER(M,NW,NS) = RMULADJ*RLDTMP(NW)
          enddo
          V_MLSER(M,IFCB,NS) = V_MLSER(M,IFCB,NS)/1000.
        enddo
      enddo
      close(1)
                
    endif ! *** end read on master
    
    ! *** Broadcast to all processes
    call Broadcast_Array(NMLSER , master_id)
    call Broadcast_Array(TCWQPSR, master_id)
    call Broadcast_Scalar(TAWQPSR, master_id)
    call Broadcast_Array(T_MLSER, master_id)
    call Broadcast_Array(V_MLSER, master_id)
    call Broadcast_Scalar(NPSTMSR, master_id)
    call Broadcast_Array(MWQPTLT, master_id)
  endif
  
  GOTO 901


  900 continue
  write(6,601)NS,M
  call STOPP('.')

  901 continue
    1 FORMAT(120X)
  601 FORMAT(' READ ERROR WQPS TIME SERIES, NSER,MDATA = ',2I5)
  602 FORMAT(' READ OF FILE WQPSL.INP SUCCESSFUL'/)
      
      
  ! *** BEGINNING OF UPDATING MASS LOADING TO THE CURRENT TIME
  1000 continue

  ! ***  INITIALIZE NULL SERIES LOADING TO ZERO
  do NW = 1,NWQV
    WQPSSRT(NW,0) = 0.
  enddo

  ! ***  LOADING SERIES INTERPOLTATION
  TIME = TIMESEC/86400.

  do NS = 1,NPSTMSR
    TIME = TIMESEC/TCWQPSR(NS)

    M2 = MWQPTLT(NS)
    do while (TIME > T_MLSER(M2,NS))
      M2 = M2+1
      if( M2 > NMLSER(NS) )then
        M2 = NMLSER(NS)
        exit
      endif
    enddo
    MWQPTLT(NS) = M2  
    M1 = M2-1
    TDIFF = T_MLSER(M2,NS)-T_MLSER(M1,NS)
    WTM1 = (T_MLSER(M2,NS)-TIME)/TDIFF
    WTM2 = (TIME-T_MLSER(M1,NS))/TDIFF
    do NW = 1,NWQV
      WQPSSRT(NW,NS) = WTM1*V_MLSER(M1,NW,NS) + WTM2*V_MLSER(M2,NW,NS)
    enddo
  enddo

  if(process_id == master_id )then
    if( ITNWQ == 0 .and. DEBUG )then
      open(1,FILE = OUTDIR//'WQPSLT.DIA',STATUS = 'UNKNOWN')
      close(1,STATUS = 'DELETE')
      open(1,FILE = OUTDIR//'WQPSLT.DIA',STATUS = 'UNKNOWN')
      write(1,112)NITER,TIME
      do NS = 1,NPSTMSR
        write(1,111)NS,(WQPSSRT(NW,NS),NW = 1,NWQV)
      enddo
      close(1)
    endif
  endif

  ! ***  COMBINE CONSTANT AND TIME VARIABLE PS LOADS
  ! M.R. MORTON 02/20/1999
  ! MODIFIED SO MULTIPLE POINT SOURCES CAN BE ADDED TO ANY GRID CELL
  ! AND ANY LAYER (HAD TO CHANGE WQWPSL ARRAY FROM 2D TO 3D).
  
  if( ITNWQ == 0 )then
    do NW = 1,NWQV
      do K = 1,KC
        do L = 2,LA
          WQWPSL(L,K,NW) = 0.0
        enddo
      enddo
    enddo
  
    if(process_id == master_id )then
        open(1,FILE = OUTDIR//'WQPSL.DIA',STATUS = 'UNKNOWN')
        close(1,STATUS = 'DELETE')
        open(1,FILE = OUTDIR//'WQPSL.DIA',STATUS = 'UNKNOWN')
        write(1,112)NITER,TIME
    endif
    
  endif

  ! *** ZERO THE ACTIVE BOUNDARY CELLS
  do NS = 1,NWQPS
    L = LIJ(ICPSL(NS), JCPSL(NS))
    WQWPSL(L,:,:) = 0.0
  enddo

  ! *** LOOP OVER THE WQ BOUNDARY CELLS
  do NS = 1,NWQPS
    L = LIJ(ICPSL(NS), JCPSL(NS))
    K = KCPSL(NS)
    ITMP = MVPSL(NS)
    ! *** write out some debug code?
    if(process_id == master_id )then
      if( ITNWQ == 0 ) WRITE(1,121) NS, L, ICPSL(NS), JCPSL(NS), K, ITMP
    endif
    
    if( K > 0 )then
      ! *** K>0, ASSIGN A SPECIFIC LAYER
      if( K < KSZ(L) )K = KSZ(L)  ! *** FORCE TO A VALID LAYER
      do NW = 1,NWQV
        WQWPSL(L,K,NW) = WQWPSL(L,K,NW) + WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP)
      enddo
    else
      ! *** K = 0, DISTRIBUTE OVER ALL THE LAYERS
      do KK = KSZ(L),KC
        do NW = 1,NWQV
          WQWPSL(L,KK,NW) = WQWPSL(L,KK,NW) + DZC(L,KK)*( WQWPSLC(NS,NW) + WQPSSRT(NW,ITMP) )
        enddo
      enddo
    endif
  enddo

  if(process_id == master_id )then
    if( ITNWQ == 0 )then
      do L = 2,LA
        ITMP = IWQPSC(L,1)
        if( ITMP > 0 )then
          do K = 1,KC
            write(1,110) ITMP, IL(L), JL(L), K, (WQWPSL(L,K,NW),NW = 1,NWQV)
          enddo
        endif
      enddo
      close(1)
    endif
  endif
  
  110 FORMAT(1X,4I4,2X,7E12.4,/,19X,7E12.4,/,19X,20E12.4)
  111 FORMAT(1X,I4,2X,7E12.4,/,7X,7E12.4,/,7X,20E12.4)
  112 FORMAT(' N, TIME = ', I10, F12.5/)
  121 FORMAT(' NS,L,I,J,K,ITMP = ', 6I5/)

  return  
  END SUBROUTINE WQPSL

  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQBENTHIC(TIMTMP)
  !
  !> @details  READ IN SPATIALLY AND/OR TEMPORALLY VARYING parameterS FOR BENTHIC
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

  use INFOMOD,only:SKIPCOM
 
  implicit none
  
  real,intent(IN) :: TIMTMP
  integer :: M,MM,IBENZ,I,L,IZM,IZS
  real    :: XBSFAD,BDAY,XM
  
  character TITLE_(3)*79, CCMRM*1
  integer,save,allocatable,dimension(:) :: IZONE
  real,save,allocatable,dimension(:) :: XBFCOD
  real,save,allocatable,dimension(:) :: XBFNH4
  real,save,allocatable,dimension(:) :: XBFNO3
  real,save,allocatable,dimension(:) :: XBFO2
  real,save,allocatable,dimension(:) :: XBFPO4D
  real,save,allocatable,dimension(:) :: XBFSAD

  if( .not. allocated(IZONE) )then
    allocate(IZONE(0:NSMZM))
    allocate(XBFCOD(0:NSMZM))
    allocate(XBFNH4(0:NSMZM))
    allocate(XBFNO3(0:NSMZM))
    allocate(XBFO2(0:NSMZM))
    allocate(XBFPO4D(0:NSMZM))
    allocate(XBFSAD(0:NSMZM))
    IZONE = 0
    XBFCOD = 0.0
    XBFNH4 = 0.0
    XBFNO3 = 0.0
    XBFO2 = 0.0
    XBFPO4D = 0.0
    XBSFAD = 0.0
  endif

  if( process_id == master_id )then
    write(*,'(A)')' WQ: WQBENFLX.INP'
    open(1,FILE = 'WQBENFLX.INP',STATUS = 'UNKNOWN')

    ! SKIP OVER THREE HEADER RECORDS:
    read(1,50) (TITLE_(M),M = 1,3)
    write(2,999)
    write(2,50) (TITLE_(M),M = 1,3)

    ! SKIP OVER ALL COMMENT CARDS AT BEGINNING OF FILE:
    REWIND(1)
    CCMRM = '#'
    call SKIPCOM(1, CCMRM,2)
    read(1, *) IBENZ
    
    write(2, 65) TIMTMP, IBENZ
    65 FORMAT(' * BENTHIC FLUXES AT     ', F10.5,' DAYS OF MODEL RUN',/, & 
              '   NUMBER OF BENTHIC FLUX ZONES = ', I4)

    ! SEQUENTIALLY READ THROUGH BENTHIC FLUX FILE UNTIL THE APPROPRIATE
    ! TIME IS FOUND:
    !   BDAY   = CURRENT DAY AT WHICH BENTHIC FLUX IS IN EFFECT
    !   BENDAY = NEXT DAY AT WHICH BENTHIC FLUX CHANGES (PASSED TO MAIN PROG

    10 READ(1, *, end = 15) BENDAY

    if( BENDAY  >  TIMTMP ) GOTO 20
    BDAY = BENDAY
    do I = 1,IBENZ
      read(1,*,end = 15) MM, XBFPO4D(MM), XBFNH4(MM), XBFNO3(MM), XBFSAD(MM), XBFCOD(MM), XBFO2(MM)
      IZONE(I) = MM
    enddo
    
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
20 continue
   
  call Broadcast_Scalar(IBENZ,    master_id)
  call Broadcast_Scalar(BENDAY,   master_id)
  
  call Broadcast_Array(IZONE,     master_id)
  call Broadcast_Array(XBFPO4D,   master_id)
  call Broadcast_Array(XBFNH4,    master_id)
  call Broadcast_Array(XBFNO3,    master_id)
  call Broadcast_Array(XBFSAD,    master_id)
  call Broadcast_Array(XBFCOD,    master_id)
  call Broadcast_Array(XBFO2,     master_id)
 
  if( process_id == master_id )then
    write(2, 48) BDAY
    48 FORMAT(/,' DAY IN BENTHIC FLUX FILE: ',F10.5,/, &
                '    ZONE    FPO4    FNH4    FNO3    FSAD    FCOD    FSOD')
    do I = 1,IBENZ
      MM = IZONE(I)
      write(2,51) MM, XBFPO4D(MM), XBFNH4(MM), XBFNO3(MM), XBFSAD(MM), XBFCOD(MM), XBFO2(MM)
    enddo

    close(1)
  endif
  
  ! DETERMINE BENTHIC FLUX FOR EACH CELL (L) BY INTERPOLATING BETWEEN
  ! THE MUD AND SAND FLUXES.  XBENMUD(L) IS THE PERCENT MUD FOR EACH CELL.
  do L = 2,LA
    IZM = IBENMAP(L,1)
    IZS = IBENMAP(L,2)
    XM = XBENMUD(L)
    WQBFPO4D(L) = XM*XBFPO4D(IZM) + (1.0-XM)*XBFPO4D(IZS)
    WQBFNH4(L)  = XM*XBFNH4(IZM)  + (1.0-XM)*XBFNH4(IZS)
    WQBFNO3(L)  = XM*XBFNO3(IZM)  + (1.0-XM)*XBFNO3(IZS)
    WQBFSAD(L)  = XM*XBFSAD(IZM)  + (1.0-XM)*XBFSAD(IZS)
    WQBFCOD(L)  = XM*XBFCOD(IZM)  + (1.0-XM)*XBFCOD(IZS)
    WQBFO2(L)   = XM*XBFO2(IZM)   + (1.0-XM)*XBFO2(IZS)
  enddo
  
  999 FORMAT(1X)
   50 FORMAT(A79)
   51 FORMAT(I8, 10F8.3)
   52 FORMAT(I7, 1X, A3)
   60 FORMAT(/, A24, I5, A24)

  return
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

  ! *** FOR THE 22 STATE VARIABLES MULTIPLIED BY THE RAINFALL FLOW RATE      !VB CHANGED 21 TO 22
  ! *** ENTERING EACH GRID CELL.  COMPUTED LOADS ARE IN G/DAY.

  implicit none

  integer :: L,NW
  real    :: CV2,TIME

  !  CV2 = CONVERSION TO GET UNITS OF G/DAY
  !  WQATM(NW) HAS UNITS OF MG/L
  !  RAINT(L) HAS UNITS OF M/SEC
  !  DXYP(L) HAS UNITS OF M2
  !  WQATML(L,KC,NW) HAS UNITS OF G/DAY

  CV2 = 86400.0
  do NW = 1,NWQV
    do L = 2,LA
      WQATML(L,KC,NW) = WQATM(NW,1)*RAINT(L)*DXYP(L)*CV2
    enddo
  enddo
  
  if( ITNWQ == 0 .and. DEBUG )then
    open(1,FILE = OUTDIR//'WQATM.DIA',STATUS = 'UNKNOWN')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'WQATM.DIA',STATUS = 'UNKNOWN')
    TIME = TIMESEC/86400.

    write(1,112) NITER,TIME
    do L = 2,LA
      write(1,110) IL(L),JL(L),(WQATML(L,KC,NW),NW = 1,NWQV)
    enddo
    close(1)
  endif
  
    110 FORMAT(1X,2I4,2X,1P,7E11.3,/,15X,7E11.3,/,15X,7E11.3)
    112 FORMAT('# WET ATMOSPHERIC DEPOSITION DIAGNOSTIC FILE',/, &
      ' N, TIME = ', I10, F12.5/)
  return  
  END SUBROUTINE WQWET
  
  !---------------------------------------------------------------------------!
  ! EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! Subroutine: Subroutine WQSKE0
  !
  !> @details  Solve Kinetic Eq from K = KC (surface layer) to K = 1 (bottom).
  !            Simplified version that only updates:
  !            IPARAM: 09 Dissolved Organic Phosphorus,
  !                    14 Ammonia Nitrogen
  !                    19 Dissolved Oxygen
  !           After computing new values, store WQVO+WQV into WQVO(L,K,NWQV)
  !           NWQV = 15,19,21.
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
    
  integer :: L,K,IZ,NS,NAL
  real    :: CNS1,TIMTMP,WQOBT0T,XMRM,WQKD0C,O2WQ_,WQTT1,WQKHN,RNH4NO3,WQTTA
  real    :: WINDREA,WQWREA,UMRM,VMRM,YMRM,WQD6,WQR6,WQVREA

  
  
  call STOPP('BAD KINETICS OPTION: ISWQLVL = 0 Does not work!  Contact DSI')   ! delme
  
  
  
  CNS1 = 2.718
  NS = 1
  do L = 2,LA
    WQI0BOT(L) = WQI0
  enddo
  !
  do K = KC,1,-1
  !
    do L = 2,LA
      TWQ(L) = max(TEM(L,K), 0.0)
      SWQ(L) = max(SAL(L,K), 0.0)
      !DZWQ(L) = 1.0 / (HPK(L,K))   deprecated 10.4
      VOLWQ(L) = HPKI(L,K) / DXYP(L)
      IMWQZT(L) = IWQZMAP(L,K)
    enddo
  !
  ! FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY
  !
    do L = 2,LA
      IWQT(L) = NINT((TWQ(L)-WQTDMIN)/WQTDINC)  ! *** DSI SINGLE L!INE
      if( IWQT(L) < 1 .or. IWQT(L) > NWQTD )then
        if( process_id == master_id )then  
          TIMTMP = TIMESEC/86400.
          
          open(3,FILE = OUTDIR//'ERROR.LOG',POSITION = 'APPEND' &
                ,STATUS = 'UNKNOWN')
          write(3,*)' *** ERROR IN WATER QUALITY'
          write(3,911) TIMTMP, L, IL(L), JL(L), K, TWQ(L)
          close(3)
          write(6,600)IL(L),JL(L),K,TWQ(L)
          IWQT(L) = max(IWQT(L),1)
          IWQT(L) = min(IWQT(L),NWQTD)
  !             call STOPP('ERROR!! INVALID WATER TEMPERATURE')
        endif
      endif
    enddo
    600 FORMAT(' I,J,K,TEM = ',3I5,E13.4)
    911 FORMAT(/,'ERROR: TIME, L, I, J, K, TWQ(L) = ', F10.5, 4I4, F10.4)
  !
    do L = 2,LA
      IZ = IWQZMAP(L,K)
  !
  ! UPDATE SOLAR RADIATION AT BOTTOM OF THIS LAYER
  !
      !WQF2IM = WQF2IM * PSHADE(L)      PMC
  !
  ! ALGAL BASAL METABOLISM & PREDATION
  !
      WQBM(L,1) = ALGAES(1).WQBMRA(IMWQZT(L)) * WQTDR(IWQT(L),1)          ! *** Basal metabolism temperature adjustment
      
      WQPR(L,1) = ALGAES(1).WQPRRA(IMWQZT(L)) * WQTDR(IWQT(L),1)          ! *** Predation temperature adjustment
      
      !*** When zooplankton is actived, predation variables of algaes inluding the predation rate, 
      !*** the fractions of producing nutrient will represent the death process of algaes.
      !*** The death rate also include the effect of temperature adjusted predation
  
      ! *** THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A
      ! *** LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO
      ! *** BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.
      WQOBT0T = WQVO(L,K,20) + WQVO(L,K,21) + WQVO(L,K,22)
      WQKRPC(L) = (WQKRC + WQKRCALG*WQOBT0T) * WQTDHDR(IWQT(L))   ! *** Hydrolysis --> DOC
      WQKLPC(L) = (WQKLC + WQKLCALG*WQOBT0T) * WQTDHDR(IWQT(L))   ! *** Hydrolysis --> DOC
      XMRM = 0.0
      do NAL = 1, NALGAE
        if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
          ! *** Macrophytes and periphyton
          XMRM = XMRM + ALGAES(NAL).WQKDCALM(1) * WQVO(L,K,19+NAL)
        endif
      enddo

      ! M. MORTON 08/28/99: ADDED SPATIALLY VARIABLE DOC HYDROLYSIS RATE WQKDC
      !    TO ACHIEVE BETTER CONTROL IN SYSTEMS WITH A COMBINATION OF FRESHWAT
      !    STREAMS AND TIDAL RIVERS WITH DIFFERENT CHARACTERISTICS.
      WQKD0C = (WQKDC(1) + WQKDCALG*WQOBT0T + XMRM)*WQTDMNL(IWQT(L))
      O2WQ_ = max(WQVO(L,K,16), 0.0)
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
      
      if( RNH4NO3_ == 0.0 )then
        do NAL = 1,NALGAE
          WQPN(L,NAL) = 0.0      
        enddo
      else
        do NAL = 1,NALGAE
          WQTTA = RNH4WQ_/(ALGAES(NAL).WQKHNA+RNO3WQ_+ 1.E-18)
          WQPN(L,NAL) = (RNO3WQ_/(ALGAES(NAL).WQKHNA+RNH4WQ_+ 1.E-18) + ALGAES(NAL).WQKHNA/(RNH4NO3_+ 1.E-18)) * WQTTA
        enddo
      endif
      WQNIT(L) = WQTDNIT(IWQT(L)) * (O2WQ(L) / (WQKHNDO + O2WQ(L) + 1.E-18)) * (RNH4WQ_ / (WQKHNN + RNH4WQ_ + 1.E-18))
        
      WQDOS(L) = DO_SAT(L)
      XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*HPK(L,K)
      if( K == KC )then
  !
  ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION:
  !
        WINDREA = WINDST(L)
        WQWREA = 0.728*SQRT(WINDREA) + (0.0372*WINDREA-0.317)*WINDREA
  !
  !        WQWREA = 0.728*SQRT(WINDST(L))
  !
        if( IWQKA(IZ)  ==  0 )then
          WQVREA = WQKRO(IZ)
          WQWREA = 0.0
        endif
  !
  !                 WIND VELOCITY COMPUTED ABOVE:
  !
        if( IWQKA(IZ)  ==  1 )then
          WQVREA = WQKRO(IZ)
        endif
  !
  !    WQKRO = 3.933 TYPICALLY
  !
        if( IWQKA(IZ)  ==  2 )then
          UMRM = 0.5*(U(L,K)+U(LEC(L),K))
          VMRM = 0.5*(V(L,K)+V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**0.5
        endif
  !
  !    WQKRO = 5.32 TYPICALLY
  !
        if( IWQKA(IZ)  ==  3 )then
          UMRM = max(U(L,K), U(LEC(L),K))
          VMRM = max(V(L,K), V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85
        endif
  !
  ! MODIFIED OWENS AND GIBBS REAERATION EQUATION:
  ! NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE
  !       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER
  !       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.
  ! WQKRO = 5.32 TYPICALLY
  !
        if( IWQKA(IZ)  ==  4 )then
          UMRM = max(U(L,K), U(LEC(L),K))
          VMRM = max(V(L,K), V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L) + 0.1524))
          WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85
        endif
        if( IWQKA(IZ)  ==  5 )then
          UMRM = max(U(L,K), U(LEC(L),K))
          VMRM = max(V(L,K), V(LNC(L),K))
          XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)
          WQVREA = 3.7*XMRM
        endif
  !
  ! NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS:
  !
        WQVREA = WQVREA * REAC(IZ)
        WQWREA = WQWREA * REAC(IZ)
        WQP19(L) = - (WQVREA + WQWREA) * HPKI(L,K)* WQTDKR(IWQT(L),IZ)
        WQKRDOS(L) = - WQP19(L)*WQDOS(L)
      else
        WQP19(L) = 0.0
      endif
    666 FORMAT(' K,IWQ,IZ,WQTDKR = ',3I5,E12.4)
  enddo
  !
  ! TRAPEZOIDAL SOLUTION OF KINETIC EQS: AFTER COMPUTING NEW VALUES, STORE
  ! WQVO+WQV INTO WQVO(L,K,NWQV)
  !
    do L = 2,LA
      IZ = IWQZMAP(L,K)
      WQD6 = - WQKHR(L)
      WQKK(L) = 1.0 / (1.0 - DTWQO2*WQD6)
      WQR6 = (WQWDSL(L,K,3) + WQWPSL(L,K,3)) * VOLWQ(L)
      WQRR(L) = WQVO(L,K,3) + DTWQ*WQR6 +  DTWQO2*WQD6*WQVO(L,K,3)
      WQV(L,K,3) = SCB(L)*(WQRR(L)*WQKK(L)) + (1. - SCB(L))*WQVO(L,K,3)
      WQVO(L,K,3) = WQVO(L,K,3) + WQV(L,K,3)
    enddo
    do L = 2,LA
      WQRR(L) = (WQWDSL(L,K,11)+WQWPSL(L,K,11)) * VOLWQ(L)
    enddo
    do L = 2,LA
      WQKK(L) = 1.0 / (1.0 + DTWQO2*WQNIT(L))
      WQRR(L) = WQVO(L,K,11) + DTWQ*WQRR(L) &
          - DTWQO2*( WQNIT(L)*WQVO(L,K,11) )
      WQV(L,K,11) = SCB(L)*(WQRR(L)*WQKK(L)) + (1.-SCB(L))*WQVO(L,K,11)
      WQVO(L,K,11) = WQVO(L,K,11)+WQV(L,K,11)
    enddo
    do L = 2,LA
      WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L))
      WQRR(L) = (WQWDSL(L,K,16) + WQWPSL(L,K,16)) * VOLWQ(L)
    enddo
    if( K == KC )then
      do L = 2,LA
        WQRR(L) = WQRR(L) + WQKRDOS(L)
      enddo
    endif
    if( K == 1 )then
      do L = 2,LA
        WQRR(L) = WQRR(L) + WQBFO2(L)*HPKI(L,K)
      enddo
    endif
    do L = 2,LA
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
    enddo
  enddo
  !
  ! INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:
  !
  TIMTMP = TIMESEC/TCTMSR
  
  ! COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX = 1/WQCHLX
  ! COUPLING TO SEDIMENT MODEL
  !: EVALUATE DEP. FLUX USING NEW VALUES CAUSE IMPLICIT SCHEME IS USED IN SPM

  ! DIURNAL DO ANALYSIS
  ! LIGHT EXTINCTION ANALYSIS

   1111 FORMAT(I12,F10.4)
   1112 FORMAT(2I5,12F7.2)
   1113 FORMAT(2I5,12E12.4)
   1414 FORMAT(I12,11E12.4)
  return  
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
 
  use SHELLFISHMOD
  
  implicit none

  integer NQ, NW, NS, IZ, IMWQZ, NSTPTMP, IOBC
  integer ND, LF, LL, LP, L, K, IFLAG, NAL, IPMC, i
  integer, save :: NWQITER                   ! *** Number of kinetic iterations
  integer, save :: IRPOM(3), ILPOM(3)        ! *** Refractory and Labile POM indicies
  integer,save,allocatable,dimension(:) :: LUSED
  
  !< In oder to adapt to the general algal group, these arrays 
  !< replace the temporary variables using for computation
  !< of the kinetic processes of water quality state variables --- DKT
  
  real :: WQF1N (NALGAEM)  !< Nitrogen limiation factor for current cell's biota, by class.  Silica limitation also included
  real :: WQF2I (NALGAEM)  !< Light limiation factor for current cell's biota, by class
  real :: WQTTAA(NALGAEM)
  real :: WQGN  (NALGAEM)  !< Concentration of available inorganic nitrogen - DKT
  real :: WQGP  (NALGAEM)  !< Concentration of available orthophosphate  - DKT
  real :: WQGCO2(NALGAEM)  !< CO2 Limitation Consts added by AA
  real :: WQA6A (NALGAEM)
  real :: WQA7A (NALGAEM)
  real :: WQA8A (NALGAEM)
  real :: WQA9A (NALGAEM)
  real :: WQA10A(NALGAEM)
  real :: WQA11A(NALGAEM)
  real :: WQA12A(NALGAEM)
  real :: WQA13A(NALGAEM)
  real :: WQA14A(NALGAEM)
  real :: WQA15A(NALGAEM)
  real :: WQA19A(NALGAEM)
  real :: WQA22A(NALGAEM)
  real :: RPOMF(3), LPOMF(3)

  real TIME, RLIGHT1, RLIGHT2, CNS1, RMULTMP
  real DTWQxH, DTWQxH2, TEMFAC
  real WQAVGIO, WQSROPT
  real XMRM, YMRM, WQTT1, WQKHN
  real WQFDI0, WQHTT, WQTTT, WQTTA
  real SADWQ, WQGSD, WQTTB, WQFDM
  real UMRM, VMRM, WQVEL, WQLVF, WQF4SC, WQKDOC, WQKHP, WQTTS
  real XNUMER, XDENOM, WQLDF, WQTTM
  real WINDREA, WQWREA, WQVREA, WQAC, WQVA1C, WQRC
  real WQB4, WQA4, WQR4
  real WQC5, WQA5, WQR5
  real WQD6, WQA6, WQR6
  real WQE7, WQA7, WQR7
  real WQF8, WQA8, WQR8
  real WQF9, WQA9, WQR9
  real WQR10, WQKKL
  real WQI11, WQA11, WQR11
  real WQJ12, WQA12, WQR12
  real WQF13, WQA13, WQR13
  real WQR14, WQF14, WQA14
  real WQR15, WQA15, WQB15
  real WQM16, WQA16D, WQR16, WQR17, WQR18
  real WQA19, WQSUM, WQRea, WQPOC, WQDOC, WQNH3, WQCOD, WQRes
  real WQT20, WQR21, TIMTMP, WQTAMD
  real PPCDO, TMP22, WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC
  real WQCDREA, WQCDSUM, ALGCOUNT
  real WQKESS, EXPA0, EXPA1, WQISM                         ! VARIABLES FOR LIGHT EXTINCTION
  real PSMLMULTIPLIER
  
  real(RKD) :: TVAL1
  real(RKD), STATIC :: AVGNEXT, AVGLAST
  
  real,save,allocatable,dimension(:,:) :: WQIS
  real,save,allocatable,dimension(:)   :: WQISC
  real,save,allocatable,dimension(:)   :: WQISD
  real,save,allocatable,dimension(:)   :: WQISG
  real,save,allocatable,dimension(:)   :: WQIBOT
  real,save,allocatable,dimension(:)   :: WQI0TOP
  real,save,allocatable,dimension(:,:,:) :: WQO      !< Current plus previous time step WQV for use in trapezoidal time stepping
  real,save,allocatable,dimension(:,:,:) :: WQOLD
  
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
  if(.not. allocated(WQIS) )then
    allocate(LUSED(LCM))                 ! *** Flag to indicate whether the L index has already been used for dry BC conditions
    
    allocate(WQIS(LCM,NALGAEM))
    allocate(WQISC(LCM))  
    allocate(WQISD(LCM))  
    allocate(WQISG(LCM))  
    
    allocate(WQIBOT(LCM))                ! *** Solar Radiation at the bottom  of the current layer, accounting for shade and converted PAR
    allocate(WQI0TOP(LCM))               ! *** Solar Radiation at the surface of the current layer, accounting for shade and converted PAR
    allocate(WQO(LCM,KCM,NWQVM))         ! *** Sum of current and previous concentrations
    allocate(WQOLD(NBCSOP,KCM,0:NWQVM))  ! *** Allocate memory for variable to store concentrations at open boundaries
    
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
    
  endif
  NWQITER = NWQITER + 1
  PSMLMULTIPLIER = 1.0
  if( ISTL == 3 )then
    PSMLMULTIPLIER = (DTWQ*86400. + DT/FLOAT(NTSTBC))/DT   ! 2.0
  endif
  
  ! *** Compute the optimal average light intensity over the last three days
  if( IWQSUN  == 2 )then  
    WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2
  else
    WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3
  endif 
  
  ! *** save OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  do NQ = 0,NWQVM
    do K = 1,KC
      do IOBC = 1,NBCSOP  
        L = LOBCS(IOBC)
        WQOLD(IOBC,K,NQ) = WQV(L,K,NQ)
      enddo
    enddo  
  enddo  
  
  ! *** ZERO RATES
  do L = 2,LA
    ! *** DRY CELL BYPASS
    if( .not. LMASKDRY(L) )then
      do NAL = 1,NALGAE
        WQPA(L,NAL) = 0.
        WQBM(L,NAL) = 0.
        WQPR(L,NAL) = 0.
        
        ! *** Apply "death rate" to fixed biota if cell is dry
        if( .not. ALGAES(NAL).ISMOBILE )then
          do K = KSZ(L),KC
            if( K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L) )then
              WQAC = -ALGAES(NAL).WQPRRA(IWQZMAP(L,K))*DTWQO2
              WQVA1C = 1.0 / (1.0 - WQAC)
              WQV(L,K,19+NAL) = (WQV(L,K,19+NAL) + WQAC*WQV(L,K,19+NAL))*WQVA1C
              WQV(L,K,19+NAL) = max(WQV(L,K,19+NAL), ALGAES(NAL).WQBMIN(IWQZMAP(L,K)))
              WQO(L,K,19+NAL) = WQVO(L,K,19+NAL) + WQV(L,K,19+NAL)
            endif
          enddo
        endif
      enddo
    endif
  enddo
  
  ! *** Set vegetative growth and drag 
  do NAL = 1, NALGAE
    if( .not. ALGAES(NAL).ISMOBILE )then
      !IF( ALGAES(NAL).ISDRAG > 0 )then
      if( ALGAES(NAL).THRESHOLD /= 0 )then
        call Macro_Veg(NAL)
      endif
    endif
  enddo
  
  CNS1 = 2.718  
  NS = 1  
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, K, NQ, NW, IZ, IMWQZ, NAL)         &
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
  do ND = 1,NDM
    LF = (ND - 1)*LDMWET + 1  
    LL = min(LF+LDMWET-1,LAWET)
  
    ! COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX = 1/WQCHLX  

    ! ***************************************************************************
    ! *** INITIALIZE SOLAR RADIATION AND OPTIMAL LIGHT
    if( LDAYLIGHT .and. IWQSUN == 2 )then
      ! *** INITIAL SOLAR RADIATION AT TOP OF SURFACE LAYER (SHADING AND ICECOVER ALREADY ACCOUNTED FOR)
      do LP = LF,LL
        L = LWET(LP)
        WQI0TOP(L) = PARADJ*2.065*RADTOP(L,KC)   ! *** Solar radiation in Langleys/day, PAR adjusted
      enddo
    else
      ! *** Initial solar radiation at top of surface layer (shading and icecover already accounted for)
      do LP = LF,LL
        L = LWET(LP)
        WQI0TOP(L) = WQI0                        ! *** Solar radiation in Langleys/day 
      enddo
    endif

    ! ***************************************************************************
    ! *** REDISTRIBUTION OF SHELLFISH TO EACH LAYER *** shellfish 1
    if( ISFFARM > 0 .and. NSF > 0 )then
      do L = LF,LL
        call SHELLFISH_REDIST(L)
      enddo
    endif

        
    ! ***************************************************************************
    ! *** Set the cyanobacteria temporally/spatially varying settling rates
    do NAL = 1,NALGAE
      if( ALGAES(NAL).ISVARSETTLE > 0 )then
        do LP = LF,LL  
          L = LWET(LP)
          if( HP(L) > ALGAES(NAL).HVM ) call ALGAE_SETTLING(NAL, L)
        enddo
      endif
    enddo   
    
    ! ***************************************************************************
    ! *** DZWQ = 1/H (for a layer), VOLWQ = 1/VOL (m^-3)
    do K = KC,1,-1  
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        ! *** HPK(L,K)                    ! *** Layer thickness of a cell in meters
        ! *** HPKI(L,K)                   ! *** Inverse layer thickness
        TWQ(L)    = max(TEM(L,K), 0.0)    ! *** Layer temperature for WQ calcs (ice formation allows for small negative temperatures)
        SWQ(L)    = max(SAL(L,K), 0.0)    ! *** Layer salinity for WQ calcs              (KG/M3)
        VOLWQ(L)  = HPKI(L,K)*DXYIP(L)    ! *** Inverse volume of each cell in a layer   (1/M^3)
        IMWQZT(L) = IWQZMAP(L,K)          ! *** WQ Zone Map for current layer
      enddo  
  
      ! *** ZERO WQWPSL IF FLOWS ARE NEGATIVE.  THESE ARE HANDLED IN CALFQC (PMC)
      if( IWQPSL /= 2 )then
        do NQ = 1,NQSIJ  
          if( (QSERCELL(K,NQ) + QSS(K,NQ)) < 0.0 )then
            ! *** ZERO THE FLUX
            L = BCPS(NQ).L
            do NW = 1,NWQV
              WQWPSL(L,K,NW) = 0.0
            enddo
            if( HDRYMOVE > 0.0 )then
              L = LQSSAVED(NQ)  
              do NW = 1,NWQV
                WQWPSL(L,K,NW) = 0.0
              enddo
            endif
          endif
        enddo
      endif
    
      ! *** Determine the rate of algae leaving the cell through settling or rising
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
    
        do NAL = 1,NALGAE
          if( ALGAES(NAL).ISMOBILE )then
            ! *** Phytoplankton
            WQBSETL(L,1,NAL) = 0.0      ! *** Phytoplankton migration into current layer from layer above
            WQBSETL(L,2,NAL) = 0.0      ! *** Phytoplankton migration from current layer to the layer below
            WQBSETL(L,3,NAL) = 0.0      ! *** Phytoplankton migration from current layer to the layer above
            WQBSETL(L,4,NAL) = 0.0      ! *** Phytoplankton migration into current layer from layer below
            
            if( K == KC )then
              if( ALGAES(NAL).SETTLING(L,K) > 0  )then
                WQBSETL(L,2,NAL) = ALGAES(NAL).SETTLING(L,K)*HPKI(L,K)      ! *** Flux out the bottom
              endif
              if( ALGAES(NAL).SETTLING(L,K-1) < 0  )then
                WQBSETL(L,4,NAL) = ALGAES(NAL).SETTLING(L,K-1)*HPKI(L,K)    ! *** Flux in from layer below
              endif
              
            elseif( K == KSZ(L) )then
              if( ALGAES(NAL).SETTLING(L,K+1) > 0  )then
                WQBSETL(L,1,NAL) = ALGAES(NAL).SETTLING(L,K+1)*HPKI(L,K)    ! *** Flux in from layer above
              endif
              if( ALGAES(NAL).SETTLING(L,K) < 0  )then
                WQBSETL(L,3,NAL) = ALGAES(NAL).SETTLING(L,K)*HPKI(L,K)      ! *** Flux out to the layer above
              else
                WQBSETL(L,2,NAL) = ALGAES(NAL).SETTLING(L,K)*HPKI(L,K)      ! *** Flux to the sediment
              endif
            
            else
              ! *** Middle layers
              if( ALGAES(NAL).SETTLING(L,K+1) > 0  )then
                WQBSETL(L,1,NAL) = ALGAES(NAL).SETTLING(L,K+1)*HPKI(L,K)    ! *** Flux in from layer above
              endif
              if( ALGAES(NAL).SETTLING(L,K) < 0  )then
                WQBSETL(L,3,NAL) = ALGAES(NAL).SETTLING(L,K)*HPKI(L,K)      ! *** Flux out to the layer above
              else
                WQBSETL(L,2,NAL) = ALGAES(NAL).SETTLING(L,K)*HPKI(L,K)      ! *** Flux out the bottom
              endif
              if( ALGAES(NAL).SETTLING(L,K-1) < 0  )then
                WQBSETL(L,4,NAL) = ALGAES(NAL).SETTLING(L,K-1)*HPKI(L,K)    ! *** Flux in from layer below
              endif
            endif
          endif
        enddo           

        ! *** ZONE SPECIFIC SETTING VELOCITIES for POM, (m/day)   
        WQRPSET(L,1) = WQWSRP(IMWQZT(L))*HPKI(L,K)  ! *** Refractory POM  (1/day)
        WQLPSET(L,1) = WQWSLP(IMWQZT(L))*HPKI(L,K)  ! *** Labile POM      (1/day)
      enddo

      ! *** SET SETTLING FOR TAM SORPTION: CURRENT LAYER  
      if( IWQSRP == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQWSSET(L,1) = WQWSS(IMWQZT(L))*HPKI(L,K)  
        enddo  
      endif  

      if( K /= KC )then
        ! *** Flux of particulate into the layer from the layer above
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRPSET(L,2) = WQWSRP(IWQZMAP(L,K+1))*HPKI(L,K)    ! *** (1/day)
          WQLPSET(L,2) = WQWSLP(IWQZMAP(L,K+1))*HPKI(L,K)    ! *** (1/day)  
        enddo  
      endif

      ! *** Set settling for tam sorption: One layer up
      if( IWQSRP == 1 .and. K /= KC )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQWSSET(L,2) = WQWSS(IWQZMAP(L,K+1))*HPKI(L,K)  
        enddo  
      endif

      ! *** FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY  
      do LP = 1,LLWET(K,ND) 
        L = LKWET(LP,K,ND)  
        IWQT(L) = NINT((TWQ(L) - WQTDMIN)/WQTDINC) + 1  
        
        if( IWQT(L) < 1 .or. IWQT(L) > NWQTD )then
          write(6,600) Map2Global(L).IG, Map2Global(L).JG, K, TWQ(L), HP(L)
          if( process_id == master_id )then
            open(1,FILE = OUTDIR//'ERROR.LOG',POSITION = 'APPEND',STATUS = 'UNKNOWN')
            write(1,*)' *** ERROR IN WQSKE1:TEMPERATURE LOOKUP TABLE'
            write(1,911) TIMEDAY, Map2Global(L).LG, Map2Global(L).IG, Map2Global(L).JG, K, TWQ(L),TEM1(L,K), HP(L), H1P(L)
            write(1,'(A)')'SURROUNDING DEPTHS'
            write(1,'(2X,A14,I5,4F14.4)')'HP  ',L,HP(LWC(L)), HP(LEC(L)), HP(LSC(L)), HP(LNC(L))
            write(1,'(2X,A14,I5,4F14.4)')'H1P ',L,H1P(LWC(L)),H1P(LEC(L)),H1P(LSC(L)),H1P(LNC(L))
            write(1,'(A)')'FLUX TERMS'
            write(1,'(2X,A14,I5,4E14.6)')'UHDYE/VHDXE'  ,L,UHDYE(L), UHDYE(LEC(L)), VHDXE(L), VHDXE(LNC(L))
            write(1,'(2X,A14,I5,4E14.6)')'UHDY1E/VHDX1E',L,UHDY1E(L),UHDY1E(LEC(L)),VHDX1E(L),VHDX1E(LNC(L))
            close(1,STATUS = 'KEEP')
          endif

          IWQT(L) = max(IWQT(L),1)
          IWQT(L) = min(IWQT(L),NWQTD)

        endif  
      enddo  

      600 FORMAT(' WQ TEM LOOKUP TABLE ERROR:  I,J,K,TEM,HP = ',3I5,4E12.4)  
      911 FORMAT('ERROR: TIME, L, I, J, K, TWQ, TEM, HP, H1P = ',F10.5, I7, 3I4, 4E12.4,/)  

      ! *** BEGIN HORIZONTAL LOOP FOR NUTRIENT AND LIGHT LIMITATION OF ALGAL GROWTH
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IZ = IWQZMAP(L,K)
        
        RNH4WQ(L) = MAX (WQV(L,K,INHX), 0.0)                                          ! *** Ammonia
        RNO3WQ(L) = MAX (WQV(L,K,INOX), 0.0)                                          ! *** Nitrate
        PO4DWQ(L) = MAX (WQPO4D(L,K), 0.0)                                            ! *** Phosphate
        RNH4NO3(L) = RNH4WQ(L) + RNO3WQ(L)                                            ! *** Total Inorganic Nitrogen
          
        ! *** Loop over biota for N/P/Si limiations
        do NAL = 1,NALGAE
          if( ALGAES(NAL).ISMOBILE  )then
            ! *** Phytoplankton
            WQGN(NAL) = RNH4NO3(L) / (ALGAES(NAL).WQKHNA + RNH4NO3(L)+ 1.E-18)
            WQGP(NAL) = PO4DWQ(L)  / (ALGAES(NAL).WQKHPA + PO4DWQ(L) + 1.E-18)
            WQF1N(NAL) = min(WQGN(NAL), WQGP(NAL))                                  ! *** Minimum of the N/P 
            if( ISKINETICS(ICO2) > 0  )then
              CO2WQ(L)  = MAX (WQV(L,K,ICO2), 0.0)                                  ! *** CO2 
              WQGCO2(NAL) = CO2WQ(L) / (ALGAES(NAL).WQKHCO2 + CO2WQ(L) + 1.E-18) 
              WQF1N(NAL) = min(WQGN(NAL), WQGP(NAL), WQGCO2(NAL))                   ! *** Minimum of the N/P/CO2     
            endif
          endif
            
          if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
            ! *** Calculate nutrient limitations for Macrophytes and periphyton
            WQGN(NAL) = RNH4NO3(L) / (ALGAES(NAL).WQKHNA + RNH4NO3(L) + 1.E-18)
            WQGP(NAL) = PO4DWQ(L)  / (ALGAES(NAL).WQKHPA + PO4DWQ(L)  + 1.E-18)
            WQF1N(NAL) = min(WQGN(NAL), WQGP(NAL))                                  ! ***  Minimum of the N/P
            if( ISKINETICS(ICO2) > 0  )then
              WQGCO2(NAL) = CO2WQ(L) / (ALGAES(NAL).WQKHCO2 + CO2WQ(L) + 1.E-18)
              WQF1N(NAL) = min(WQGN(NAL), WQGP(NAL), WQGCO2(NAL))                   ! ***  Minimum of the N/P/CO2
            endif
          endif
            
          if( IWQSI == 1 )then  
            ! *** Silica limitation
            SADWQ = MAX (WQSAD(L,K), 0.0)
            WQGSD = SADWQ / (ALGAES(NAL).WQKHS + SADWQ + 1.E-18)
            if( ALGAES(NAL).ISILICA /= 0 )then
              WQF1N(NAL) = min(WQF1N(NAL),WQGSD)                                    ! *** Limit: Diatoms - Minimum of the N/P/S
            endif
          endif                
        enddo   ! *** End of biota loop
                    
        ! *** LIGHT EXTINCTION (THIS WILL ALWAYS BE TRUE EXCEPT FOR IWQSUN = 2)
        if( WQI0 > 0.1 )then
          ! *** GET THE EXTINCTION COEFFICIENT
          WQKESS = RADKE(L,K)
        
          ! *** OPTIMAL LIGHT INTENSITY AT OPTIMAL DEPTH
          if( K == KC )then
            ! *** Only compute optimal light once, even if optimal depth is below layer KC
            do NAL = 1,NALGAE
              if( ALGAES(NAL).ISMOBILE )then
                ! *** Phytoplankton
                WQIS(L,NAL) = max( WQAVGIO*EXP(-WQKESS*ALGAES(NAL).WQDOP(IZ)), WQISMIN )  
                
                ! *** HARDWIRE TO SET OPTIMAL GROWTH TO A FIXED VALUE USED BY CHAPRA (Algal Growth Example (Chapra 33.2 with Fixed IS250).
                !WQIS(L,NAL) = 250./WQFD         ! *** Is
              endif
            enddo
          endif

          ! *** Current light growth limiting factor. 
          if( K == KC )then    
            WQITOP(L,K) = WQI0TOP(L)                            ! *** WQITOP is solar radiation at the TOP of layer K
            WQIBOT(L)   = WQI0TOP(L)*EXP(-WQKESS*HPK(L,K))      ! *** WQIBOT is solar radiation at the BOTTOM of layer K
          else            
            WQITOP(L,K) = WQIBOT(L)
            WQIBOT(L)   = WQITOP(L,K)*EXP(-WQKESS*HPK(L,K)) 
          endif !SEE DiTORO ET AL (1971, EQNS. (11)&(12)) 

          WQTT1 = WQFD*2.718                                     ! *** EXP(1.0) = 2.718
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              EXPA0 = EXP(-WQITOP(L,K)/(WQFD*WQIS(L,NAL)))              ! *** WQITOP(L,K)/WQIS(L,NAL) is the fraction of actual SR to average SR at the TOP 
              EXPA1 = EXP(-WQIBOT(L)  /(WQFD*WQIS(L,NAL)))              ! *** WQITOP(L,K)/WQIS(L,NAL) is the fraction of actual SR to average SR at the BOTTOM
              WQF2I(NAL) = WQTT1/(HPK(L,K)*WQKESS)*(EXPA1 - EXPA0)
            else
              ! *** Macrophytes and periphyton - Light Limitation at top of growth layer
              WQF2I(NAL) = 0.0
              if( WQITOP(L,K) > 1.0E-18 )then  
                WQFDI0 = - WQIBOT(L)/(WQFD + 1.E-18)
                WQISM  = max(WQAVGIO*EXP(-WQKESS*ALGAES(NAL).WQDOP(IZ)), WQISMIN)  ! *** Optimal Light
                WQFDM  = WQFDI0/(WQISM + 1.E-18)                                   ! *** Ratio of Actual to Optimal Light
                WQHTT  = WQHT(K) * HP(L)  
                WQTTB  = EXP(-WQKESS * (WQHTT + 1.0/HPKI(L,K)))  
                WQTTT  = EXP(-WQKESS * WQHTT)  
                WQTT1  = (CNS1 * WQFD * HPKI(L,K))/WQKESS  
                WQF2I(NAL) = WQTT1 * (EXP(WQFDM*WQTTB) - EXP(WQFDM*WQTTT))         ! *** Light Based Macroalgae Growth Limiting Factor   
                WQIS(L,NAL) = WQISM
              endif
            endif
          enddo
        else
          ! *** No Light Case
          WQIBOT(L) = 0.
          WQITOP(L,K) = 0.
          WQKESS = 0.
          do NAL = 1,NALGAE
            WQF2I(NAL) = 0.0
          enddo
        endif
            
        do NAL = 1,NALGAE
          if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
            ! *** Velocity Limitation for macrophytes and periphyton
            WQLVF = 1.0
            if( ALGAES(NAL).IWQVLIM > 0 )then
              UMRM  = 0.5*(U(L,K) + U(LEC(L),K))  
              VMRM  = 0.5*(V(L,K) + V(LNC(L),K))  
              WQVEL = SQRT(UMRM*UMRM + VMRM*VMRM)  
              
              ! *** OPTION 1 FOR VELOCITY LIMITATION ASSUMES MACROALGAE GROWTH  
              ! *** IS LIMITED AT LOW VELOCITIES DUE TO REDUCED AVAILABILITY OF  
              ! *** NUTRIENTS REACHING THE ALGAE BIOMASS.    
              ! *** USES A MICHAELIS-MENTON or MONOD TYPE OF EQUATION.  
              if( ALGAES(NAL).IWQVLIM  ==  1 )then  
                if( WQVEL  > ALGAES(NAL).WQKMVMIN(IZ) )then  
                  WQLVF = WQVEL / (ALGAES(NAL).WQKMV(IZ) + WQVEL)  
                else  
                  WQLVF = ALGAES(NAL).WQKMVMIN(IZ) / (ALGAES(NAL).WQKMV(IZ) + ALGAES(NAL).WQKMVMIN(IZ))  
                endif  
              endif       
              
              ! *** OPTION 2 FOR VELOCITY LIMITATION APPLIES A FIVE-parameter LOGISTIC  
              ! *** FUNCTION THAT CAN BE ADJUSTED TO LIMIT MACROALGAE GROWTH FOR  
              ! *** EITHER LOW OR HIGH (SCOUR) VELOCITIES.  IN STREAMS WITH LOW NUTRIENTS,  
              ! *** THE LOW VELOCITY WILL LIKELY BE LIMITING SINCE AMPLE NUTRIENTS MAY  
              ! *** NOT REACH THE ALGAE BIOMASS DUE TO REDUCED FLOW.  IN STREAMS WITH  
              ! *** ABUNDANT NUTRIENTS, LOW VELOCITIES WILL NOT LIMIT MACROALGAE GROWTH,  
              ! *** INSTEAD, HIGH VELOCITIES WILL LIKELY SCOUR THE MACROALGAE AND DETACH  
              ! *** IT FROM THE SUBSTRATE.  
              if( ALGAES(NAL).IWQVLIM  == 2 )then  
                XNUMER = ALGAES(NAL).WQKMVA(IZ) - ALGAES(NAL).WQKMVD(IZ)  
                XDENOM = 1.0 + (WQVEL/ALGAES(NAL).WQKMVC(IZ))**ALGAES(NAL).WQKMVB(IZ)  
                WQLVF = ALGAES(NAL).WQKMVD(IZ) + ( XNUMER / (XDENOM**ALGAES(NAL).WQKMVE(IZ)) )  
              endif  
            endif
              
            ! *** Use THE MORE SEVERELY LIMITING OF VELOCITY OR NUTRIENT FACTORS:  
            WQF1N(NAL) = min(WQLVF, WQF1N(NAL)) 

            ! *** Crowding limitation factor based on WQKBP optimal plant density
            XMRM = WQV(L,K,19+NAL)*HPK(L,K)          ! *** Convert WQV (gC/M3) to a density: XMRM (gC/M2)  
            WQLDF = ALGAES(NAL).WQKBP(IZ) / (ALGAES(NAL).WQKBP(IZ) + XMRM)
              
            ! ***                 Max growth rate       Nutr/Vel    Light     Temperature       Crowding
            WQPA(L,NAL) = ALGAES(NAL).WQPMA(IMWQZT(L))*WQF1N(NAL)*WQF2I(NAL)*WQTDG(IWQT(L),NAL)*WQLDF   ! *** Macroalgae growth rate
              
            WQBM(L,NAL) = ALGAES(NAL).WQBMRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted macroalgae metabolism rate
            WQPR(L,NAL) = ALGAES(NAL).WQPRRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted macroalgae predation rate 
            
          endif
        enddo
          
        ! *** Compute the Growth Rate based on Maximums & Limiting Factors             
        do NAL = 1,NALGAE
          if( ALGAES(NAL).ISMOBILE )then
            ! *** Phytoplankton
            WQPA(L,NAL) = ALGAES(NAL).WQPMA (IMWQZT(L)) * WQF1N(NAL)*WQF2I(NAL)*WQTDG(IWQT(L),NAL)      ! *** Nutrient & Temperature Adjusted
            WQBM(L,NAL) = ALGAES(NAL).WQBMRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted
            WQPR(L,NAL) = ALGAES(NAL).WQPRRA(IMWQZT(L)) * WQTDR(IWQT(L),NAL)                            ! *** Temperature Adjusted
          endif
            
          if( ALGAES(NAL).ISTOX /= 0 )then 
            WQF4SC  = ALGAES(NAL).WQSTOX / (ALGAES(NAL).WQSTOX + SWQ(L)*SWQ(L) + 1.E-12)
            WQPA(L,NAL) = WQPA(L,NAL) * WQF4SC
          endif
            
          ! THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A  
          ! LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO  
          ! BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.  
          if( ALGAES(NAL).ISBLOOM /= 0 )then
            WQBM(L,NAL) = WQBM(L,NAL)*WQTDP(IWQT(L),NAL)                                                ! *** Temperature Adjusted for winter bloom
            WQPR(L,NAL) = WQPR(L,NAL)*WQTDP(IWQT(L),NAL)                                                ! *** Temperature Adjusted for winter bloom
          endif
          !*** When zooplankton is actived, predation variables of algaes inluding the predation rate, 
          !*** the fractions of producing nutrient will represent the death process of algaes.
          !*** The death rate also include the effect of temperature adjusted predation
        enddo
      enddo      ! *** END ACTIVE CELL LOOP FOR ALGAE parameterS  
        
      ! *** SET UP OTHER NUTRIENT VARIABLES FOR MAIN KINETICS SECTION
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        IZ = IWQZMAP(L,K)
        WQOBTOT(L) = 0.0
        XMRM = 0.0
        
        do NAL = 1,NALGAE
          if( ALGAES(NAL).ISMOBILE )then           ! *** Phytoplankton
            WQOBTOT(L) = WQOBTOT(L) + WQV(L,K,19+NAL)           
          elseif( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then    ! *** Fixed biota
            XMRM = XMRM + ALGAES(NAL).WQKDCALM(IZ) * WQV(L,K,19+NAL)  
          endif
        enddo
        
        WQKRPC(L) = (WQKRC + WQKRCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))        ! *** Hydrolysis:     RPOP--> DOC
        WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))        ! *** Hydrolysis:     RPOP--> DOC
        
        ! *** M. MORTON added spatially variable doc minerization rate WQKDC  
        ! *** to achieve better control in systems with a combination of freshwater
        ! *** streams and tidal rivers with different characteristics.  

        ! ***   Min Mineral     Factor   MobC   Factor*MacC   T Correct
        WQKDOC = (WQKDC(IZ) + WQKDCALG*WQOBTOT(L) + XMRM)*WQTDMNL(IWQT(L))     ! *** Heterotrophic respiration rate of DOC at infinite DO (1/day)
        
        O2WQ(L) = max(WQV(L,K,IDOX), 0.0)  
        WQTT1 = WQKDOC / (WQKHORDO + O2WQ(L) + 1.E-18)  
        WQKHR(L)   = WQTT1*O2WQ(L)  
        WQDENIT(L) = WQTT1*WQAANOX*RNO3WQ(L)/(WQKHDNN + RNO3WQ(L) + 1.E-18)    ! *** Denitrification Rate  
      enddo  

      ! ***********************************
      ! 7-10 PHOSPHORUS  
      ! *** HYDROLYSIS
      
      WQKHP = 0.0
      ALGCOUNT = 0.
      do NAL = 1,NALGAE
        if( ALGAES(NAL).ISMOBILE )then
          WQKHP = WQKHP + ALGAES(NAL).WQKHPA                            ! *** Mean phosphorus half-saturation for algae
          ALGCOUNT =  ALGCOUNT + 1.
        endif
      enddo
      if( ALGCOUNT > 0. ) WQKHP = WQKHP/ALGCOUNT   ! *** Note - better if mass weighted
      
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        WQAPC(L) = 1.0/(WQCP1PRM + WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ(L)))   ! *** A mean Phosphorus-to-Carbon ration for all algal groups
        WQTT1 = WQKHP/(WQKHP + PO4DWQ(L) + 1.E-18) * WQOBTOT(L)  
        WQKRPP(L) = (WQKRP + WQKRPALG*WQTT1) * WQTDHDR(IWQT(L))         ! *** Hydrolysis:     RPOP--> DOP
        WQKLPP(L) = (WQKLP + WQKLPALG*WQTT1) * WQTDHDR(IWQT(L))         ! *** Hydrolysis:     LPOP--> DOP
        WQKDOP(L) = (WQKDP + WQKDPALG*WQTT1) * WQTDMNL(IWQT(L))         ! *** Mineralization: DOP --> PO4
      enddo
  
      ! *** PHOSPHATE SETTLING   Note - THIS COULD BE SPED UP BY CHECKING OPTIONS OUTSIDE OF THE L LOOP
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        if( IWQSRP == 1 )then  
          WQTTM = WQKPO4P*WQTAMP(L,K)                                   ! *** Sorbed mass fraction to TAM
          WQH10(L) = - WQWSSET(L,1) * WQTTM / (1.0 + WQTTM)             ! *** Loss from the layer out bottom
          if( K /= KC )then  
            WQTTM = WQKPO4P*WQTAMP(L,K+1)  
            WQT10(L) = WQWSSET(L,2) * WQTTM / (1.0 + WQTTM)             ! *** Gain to the layer from above
          endif  
        elseif( IWQSRP == 2 .and. ISTRAN(6) > 0 )then  
          WQTTS = WQKPO4P*SEDT(L,K)                                     ! *** Sorbed mass fraction to cohesive sediments  
          WQH10(L) = - WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)    ! *** Loss from the layer 
          if( K /= KC )then  
            WQTTS = WQKPO4P*SEDT(L,K+1)  
            WQT10(L) = WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)    ! *** Gain to the layer from above  
          endif  
        else  
          WQH10(L) = 0.0  
          WQT10(L) = 0.0  
        endif 
        WQH10(L) = WQH10(L)*DTWQO2 
      enddo  

      ! ***********************************
      ! 11-15 NITROGEN  
      ! *** HYDROLYSIS
      ! *** Mean nitrogen half-saturation for algae
      WQKHN = 0.0
      ALGCOUNT = 0.
      do NAL = 1,NALGAE
        if( ALGAES(NAL).ISMOBILE )then
          WQKHN = WQKHN + ALGAES(NAL).WQKHNA
          ALGCOUNT =  ALGCOUNT + 1.
        endif
      enddo
      if( ALGCOUNT > 0. ) WQKHN = WQKHN/ALGCOUNT   ! *** Note - better if mass weighted
      
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        WQTT1 = WQKHN / (WQKHN + RNH4NO3(L)+ 1.E-18) * WQOBTOT(L)  
        WQKRPN(L) = (WQKRN + WQKRNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** Hydrolysis:     RPON-->DON
        WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** Hydrolysis:     LPON-->DON
        WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** Mineralization: DON -->NH4
      enddo
      
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        
        ! *** Ammonium Preference
        do NAL = 1, NALGAE
          WQTTA = RNH4WQ(L)/(ALGAES(NAL).WQKHNA + RNO3WQ(L) + 1.E-18)
          WQPN(L,NAL) = ( RNO3WQ(L)/(ALGAES(NAL).WQKHNA + RNH4WQ(L) + 1.E-18) + ALGAES(NAL).WQKHNA/(RNH4NO3(L) + 1.E-18) ) * WQTTA   ! *** Ammonium preference
        enddo
        WQNIT(L) = WQTDNIT(IWQT(L)) * O2WQ(L) / (WQKHNDO + O2WQ(L) + 1.E-18) * RNH4WQ(L) / (WQKHNN + RNH4WQ(L) + 1.E-18)
      enddo
      
      if( IWQSI == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          if( IWQSRP == 1 )then  
            WQTTM = WQKSAP*WQTAMP(L,K)  
            WQN17(L) = - WQWSSET(L,1) * WQTTM / (1.0 + WQTTM)  
            if( K /= KC )then  
              WQTTM = WQKSAP*WQTAMP(L,K+1)  
              WQT17(L) = WQWSSET(L,2) * WQTTM / (1.0 + WQTTM)  
            endif  
          elseif( IWQSRP == 2 .and. ISTRAN(6) > 0 )then  
            WQTTS = WQKSAP*SEDT(L,K)  
            WQN17(L) = - WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)  
            if( K /= KC )then  
              WQTTS = WQKSAP*SEDT(L,K+1)  
              WQT17(L) = WSEDO(NS) * WQTTS * HPKI(L,K) / (1.0 + WQTTS)  
            endif  
          else  
            WQN17(L) = 0.0  
            WQT17(L) = 0.0  
          endif  
        enddo  
        WQN17(L) = WQN17(L)*DTWQO2 
      endif  

      ! ***********************************
      ! *** DISSOLVED OXYGEN
      PPCDO = -3.45  !PARTIAL PRES OF CO2 IN 10^ppcdo ATM; TEMPORARILY DECLARED HERE. SHLD BE READ IN FROM INPUT FILE
      do LP = 1,LLWET(K,ND)
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
        if( K == KC )then       
          ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION
          WINDREA = min(WINDST(L),11.)  
          WQWREA = 0.728*SQRT(WINDREA) + (0.0372*WINDREA - 0.317)*WINDREA  
          WQWREA = WQWREA * HPKI(L,K)
          
          if( IWQKA(IZ) == 0 )then
            ! *** Constant  
            WQVREA = WQKRO(IZ) * HPKI(L,K) 
            WQWREA = 0.0  
          elseif( IWQKA(IZ) == 1 )then  
            ! *** Constant plus Wind
            WQVREA = WQKRO(IZ) * HPKI(L,K) 
          elseif( IWQKA(IZ) == 2 )then
            ! *** OCONNOR-DOBBINS REAERATION FORMULA
            UMRM = 0.5*(U(L,K) + U(LEC(L),K))  
            VMRM = 0.5*(V(L,K) + V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            ! *** WQKRO = 3.933 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**0.5         ! *** Multi-layer implementation  
            WQVREA = WQVREA * HPKI(L,K)                         ! *** Multi-layer implementation
          elseif( IWQKA(IZ) == 3 )then
            ! *** OWENS & GIBBS (1964) REAERATION FORMULA
            UMRM = max(U(L,K), U(LEC(L),K))  
            VMRM = max(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            ! *** WQKRO = 5.32 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85  
          elseif( IWQKA(IZ) == 4 )then  
            ! *** MODIFIED OWENS AND GIBBS REAERATION EQUATION:  
            ! *** NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE  
            ! ***       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER  
            ! ***       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.  
            UMRM = max(U(L,K), U(LEC(L),K))  
            VMRM = max(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L) + 0.1524))  
            ! *** WQKRO = 5.32 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85  
          elseif( IWQKA(IZ) == 5 )then  
            UMRM = max(U(L,K), U(LEC(L),K))  
            VMRM = max(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            WQVREA = 3.7*XMRM * HPKI(L,K)  
          endif  

          ! *** NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS
          WQVREA = WQVREA * REAC(IZ)  ! *** Diffusive Flux
          WQWREA = WQWREA * REAC(IZ)  ! *** Wind Component
          WQP19(L)    = -(WQVREA + WQWREA) * WQTDKR(IWQT(L),IZ)   ! *** Reaeration rate (1/day)  
          WQKRDOS(L)  = -WQP19(L)*WQDOS(L)                        ! *** O2 Flux due to DO saturation (offset by current DO below)
          WQP22(L)    = WQP19(L)*((32./44.)**0.25)                ! *** Kr FOR CO2 ANALOGOUS TO WQP19 ; 44 = MOL WT OF CO2
          WQKRCDOS(L) = -WQP22(L) * WQCDOS(L)                     ! *** EVALUATING Kr*SAT CONC OF CO2
          WQRREA(L)   = WQP19(L)                                  ! *** Store reaeration rate for array out writing   
        else
          WQKRDOS(L)  = 0.0
          WQKRCDOS(L) = 0.0  
          WQP19(L)    = 0.0  
          WQP22(L)    = 0.0              !VB Kr FOR CO2 IS ZERO FOR CELLS NOT AT THE SURFACE
        endif  
      enddo  
  
      ! *** ICE
      if( ISICE > 0 .and. K == KC )then
        ! *** REDUCE REAERATION AND SURFACE EXCHANGE DUE TO ICE
        do LP = LF,LL
          L = LWET(LP)
          WQP19(L)   = (1.-ICECOVER(L))*WQP19(L)
          WQKRDOS(L) = (1.-ICECOVER(L))*WQKRDOS(L)
          WQRREA(L)  = WQP19(L)        ! *** Store reaeration rate for array out writing   
        enddo
      endif
    
      ! *** Trapezoidal solution of kinetic eqs: After computing new values, store WQVO + WQV into WQO

      !*******************************************************************************************************************************
      ! *** MACROPHYTES AND PERIPHYTON
      do NAL = 1,NALGAE
        if( .not. ALGAES(NAL).ISMOBILE )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            if( K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L) )then
              IMWQZ = IWQZMAP(L,KSZ(L))
            
              ! *** Apply maximum allowable concentration
              if( WQV(L,K,19+NAL) > ALGAES(NAL).WQBMAX )then
                WQPA(L,NAL) = 0.0                                                                    ! *** Zero the growth
              endif
            
              ! ***    GROWTH        METABOLISM  PREDATION/DEATH           SETTLING(NEEDED?)
              WQAC = (WQPA(L,NAL) - WQBM(L,NAL) - WQPR(L,NAL) - ALGAES(NAL).SETTLING(L,K)*HPKI(L,K))*DTWQO2  

              WQVA1C = 1.0 / (1.0 - WQAC)
              WQV(L,K,19+NAL) = (WQV(L,K,19+NAL) + WQAC*WQV(L,K,19+NAL))*WQVA1C*SMAC(L)  
              WQV(L,K,19+NAL) = max(WQV(L,K,19+NAL), ALGAES(NAL).WQBMIN(IMWQZ))*SMAC(L)  
              WQO(L,K,19+NAL) = WQVO(L,K,19+NAL) + WQV(L,K,19+NAL)
            
              ! *** Hydrodynamic feedback parameters
              if( ALGAES(NAL).ISDRAG > 0 )then
                ! *** Compute vegetation equivalents for drag calculations
              
              endif
            endif
          enddo
        endif
      enddo
  
      ! *** Phytoplankton       
      do NAL = 1,NALGAE
        if( ALGAES(NAL).ISMOBILE )then
          ! *** Phytoplankton
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            ! ***   GROWTH        BASAL_METAB   PREDATION     SETTLING           TIME STEP  
            WQAC = (WQPA(L,NAL) - WQBM(L,NAL) - WQPR(L,NAL) - WQBSETL(L,2,NAL) )*DTWQO2           ! *** Production per unit time multiplied by half time step
            WQKK(L) = 1.0 / (1.0 - WQAC) 
      
            ! ***   PT_SRC_LOADS        VOLUME  
            WQRC = WQWPSL(L,K,19+NAL) * VOLWQ(L) * PSMLMULTIPLIER              ! *** Point source load rate multiplied by inverse cell volume  g/m^3/t
            WQRR(L) = WQV(L,K,19+NAL) + DTWQ*WQRC + WQAC*WQV(L,K,19+NAL)       ! *** Transported biomass conc. (CALWQC) + point source load rate X time step + growth rate X previous biomass conc.
          enddo
          
          ! *** Coupling to zooplankton 
          if( IWQZPL > 0 )then
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              WQRR(L) = WQRR(L) - DTWQ*SBZPAL(L,K,NAL)
            enddo
          endif
          
          ! *** Vertical fluxes from various sources
          if( K /= KC )then
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQBSETL(L,1,NAL)*WQO(L,K+1,19+NAL)    ! *** Vertical migration into the current layer from layer above
            enddo
          else  ! K == KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***    ATM DRY DEP           ATM WET DEP         VOLUME  
              WQRC = (WQWDSL(L,KC,19+NAL) + WQATML(L,KC,19+NAL))*VOLWQ(L)      ! *** Atmospheric loading mass per time / cell volume
              WQRR(L) = WQRR(L) + DTWQ*WQRC                                    ! *** Biomass conc. + Dt*loading rate per unit volume
            enddo
          endif 
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) - DTWQ*( WQBSETL(L,4,NAL)*WQVO(L,K-1,19+NAL)   & ! *** Vertical migration into the current layer from layer below
                                     - WQBSETL(L,3,NAL)*WQVO(L,K,19+NAL) )
            
          enddo

          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQV(L,K,19+NAL) = WQRR(L)*WQKK(L)
            WQV(L,K,19+NAL) = max(WQV(L,K,19+NAL), ALGAES(NAL).WQBMIN(1))      ! *** Apply a user specified minimum concentration
            WQO(L,K,19+NAL) = WQVO(L,K,19+NAL) + WQV(L,K,19+NAL)               ! *** Depth totaled biomass conc = old biomass conc in cell + biomass conc from this iteration
          enddo
        endif
      enddo

      !*******************************************************************************************************************************
      ! *** NOW COMPUTE KINETICS FOR EACH CONSTITUENT
      
      ! ****  PARAM 01  ROC - refractory particulate organic carbon
      if( ISKINETICS(1) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS    SETTLING  
          WQB4 = -( WQKRPC(L) + WQRPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQB4)  

          ! ***  ALGAE PREDATION SOURCE OF RPOC
          WQA4 = 0.0
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              WQA4 = WQA4 + ALGAES(NAL).WQFCRP*WQPR(L,NAL)*WQO(L,K,19+NAL)
            elseif( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA4 = WQA4 + ALGAES(NAL).WQFCRP*WQPR(L,NAL)*WQO(L,K,19+NAL)
            endif
          enddo 
          
          ! ***  PT_SRC_LOADS    VOLUME  
          WQR4 = WQWPSL(L,K,IROC) * VOLWQ(L) * PSMLMULTIPLIER
          
          WQRR(L) = WQV(L,K,IROC) + DTWQ*WQR4 + DTWQO2*WQA4 + WQB4*WQVO(L,K,IROC)
          
          !*** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SRPOCZ(L,K)
        enddo

        if( K /= KC )then  
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,IROC)  
          enddo  
        else  ! K == KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP            ATM WET DEP    VOLUME  
            WQR4 = (WQWDSL(L,KC,IROC) + WQATML(L,KC,IROC))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR4  
          enddo
        endif  
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IROC) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,IROC) = WQVO(L,K,IROC) + WQV(L,K,IROC)
        enddo  
      endif  

      ! ****  PARAM 02  LOC - labile particulate organic carbon
      if( ISKINETICS(2) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS    SETTLING  
          WQC5 = - (WQKLPC(L) + WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQC5)
          
          ! *** Algae predation source                         
          WQA5 = 0.0
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              WQA5 = WQA5 + ALGAES(NAL).WQFCLP*WQPR(L,NAL)*WQO(L,K,19+NAL)
            elseif( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then  
              WQA5 = WQA5 + ALGAES(NAL).WQFCLP*WQPR(L,NAL)*WQO(L,K,19+NAL)  
            endif
          enddo
          
          ! ***  PT_SRC_LOADS       VOLUME  
          WQR5 = WQWPSL(L,K,ILOC) * VOLWQ(L) * PSMLMULTIPLIER

          WQRR(L) = WQV(L,K,ILOC) + DTWQ*WQR5 + DTWQO2*WQA5 + WQC5*WQVO(L,K,ILOC)
          
          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SLPOCZ(L,K)
        enddo  
        
        if( K /= KC )then  
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,ILOC)  
          enddo  
        else  ! K == KC
          ! *** Add surface loads
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP       ATM WET DEP        VOLUME  
            WQR5 = (WQWDSL(L,K,ILOC) + WQATML(L,KC,ILOC))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR5  
          enddo
        endif  
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,ILOC) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,ILOC) = WQVO(L,K,ILOC) + WQV(L,K,ILOC)
        enddo  
      endif  

      ! ****  PARAM 03  DOC - dissolved organic carbon
      if( ISKINETICS(3) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          IZ = IWQZMAP(L,K) 
          
          ! ***    RESPIRATION  DENITRIFICATION
          WQD6 = -(WQKHR(L) + WQDENIT(L)) *DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQD6)
          
          ! *** Algae metabolism and predation source 
          WQA6 = 0.0
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA6A(NAL) = ALGAES(NAL).WQFCDB + (1.0 - ALGAES(NAL).WQFCDB)*(ALGAES(NAL).WQKHRA(IZ)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18))
              WQA6 = WQA6 + ( WQA6A(NAL)*WQBM(L,NAL) + ALGAES(NAL).WQFCDP*WQPR(L,NAL) )*WQO(L,K,19+NAL)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA6A(NAL) = ALGAES(NAL).WQFCDB + (1.0 - ALGAES(NAL).WQFCDB)*(ALGAES(NAL).WQKHRA(IZ)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)) 
              WQA6 = WQA6 + ( WQA6A(NAL)*WQBM(L,NAL) + ALGAES(NAL).WQFCDP*WQPR(L,NAL) )*WQO(L,K,19+NAL) 
            endif
          enddo

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR6 = WQWPSL(L,K,IDOC) * VOLWQ(L) * PSMLMULTIPLIER

          WQRR(L) = WQV(L,K,IDOC) + DTWQ*WQR6 + WQD6*WQVO(L,K,IDOC) + DTWQO2*(WQA6 + WQKRPC(L)*WQO(L,K,IROC) + WQKLPC(L)*WQO(L,K,ILOC))
          
          ! *** Coupling to zooplankton - DKT              
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDOCZ(L,K)
        enddo

        if( K == KC )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR6 = (WQWDSL(L,K,IDOC) + WQATML(L,KC,IDOC))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR6  
          enddo
        endif
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IDOC) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,IDOC) = WQVO(L,K,IDOC) + WQV(L,K,IDOC)
        enddo  
      endif  

      ! ****  PARAM 04  ROP - refractory particulate organic phosphorus
      if( ISKINETICS(4) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQE7 = -(WQKRPP(L) + WQRPSET(L,1))*DTWQO2 
          WQKK(L) = 1.0 / (1.0 - WQE7)
          
          ! *** Algal metabolism and predation source
          WQA7 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA7A(NAL) = (ALGAES(NAL).WQFPRB*WQBM(L,NAL) + ALGAES(NAL).WQFPRP*WQPR(L,NAL)) * WQO(L,K,19+NAL)
              WQA7 = WQA7 + WQA7A(NAL) * WQAPC(L)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L))  )then
              ! *** Macrophytes and periphyton
              WQA7A(NAL) = (ALGAES(NAL).WQFPRB*WQBM(L,NAL) + ALGAES(NAL).WQFPRP*WQPR(L,NAL)) * WQO(L,K,19+NAL)
              WQA7 = WQA7 + WQA7A(NAL) * WQAPC(L) * ALGAES(NAL).WQAPCM 
            endif
          enddo
          
          ! ***  PT_SRC_LOADS    VOLUME  
          WQR7 = WQWPSL(L,K,IROP) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,IROP) + DTWQ*WQR7 + DTWQO2*WQA7 + WQE7*WQVO(L,K,IROP) 
         
          ! *** Coupling to zooplankton - DKT              
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SRPOPZ(L,K)
        enddo  
        
        if( K /= KC )then  
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,IROP)
          enddo  
        else
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP    ATM WET DEP    VOLUME  
            WQR7 = (WQWDSL(L,K,IROP) + WQATML(L,KC,IROP))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR7  
          enddo
        endif  
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IROP) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,IROP) = WQVO(L,K,IROP) + WQV(L,K,IROP)
        enddo  
      endif  

      ! ****  PARAM 05  LOP - labile particulate organic phosphorus
      if( ISKINETICS(5) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***    HYDROLYSIS  SETTLING
          WQF8 = - (WQKLPP(L) + WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQF8)
          
          ! *** Algae metabolism and predation source
          WQA8 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA8A(NAL) = (ALGAES(NAL).WQFPLB*WQBM(L,NAL) + ALGAES(NAL).WQFPLP*WQPR(L,NAL)) * WQO(L,K,19+NAL)
              WQA8 = WQA8 + WQA8A(NAL) * WQAPC(L) 
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA8A(NAL) = (ALGAES(NAL).WQFPLB*WQBM(L,NAL) + ALGAES(NAL).WQFPLP*WQPR(L,NAL)) * WQO(L,K,19+NAL)
              WQA8 = WQA8 + WQA8A(NAL) * WQAPC(L) * ALGAES(NAL).WQAPCM 
            endif
          enddo
          
          ! ***  PT_SRC_LOADS    VOLUME  
          WQR8 = WQWPSL(L,K,ILOP) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,ILOP) + DTWQ*WQR8 + DTWQO2*WQA8 + WQF8*WQVO(L,K,ILOP)  
          
          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SLPOPZ(L,K)
        enddo 
        
        if( K /= KC )then  
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,ILOP)
          enddo  
        else
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP       ATM WET DEP        VOLUME  
            WQR8 = (WQWDSL(L,K,ILOP) + WQATML(L,KC,ILOP))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR8  
          enddo
        endif  
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,ILOP) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,ILOP) = WQVO(L,K,ILOP) + WQV(L,K,ILOP)
        enddo  
      endif  

      ! ****  PARAM 06  DOP - dissolved organic phosphorus
      if( ISKINETICS(6) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQF9 = - DTWQO2*WQKDOP(L)  
          WQKK(L) = 1.0 / (1.0 - WQF9)
          
          ! *** Algae metabolism and predation source
          WQA9 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              WQA9A(NAL) = (ALGAES(NAL).WQFPDB*WQBM(L,NAL) + ALGAES(NAL).WQFPDP*WQPR(L,NAL)) * WQO(L,K,19+NAL)
              WQA9 = WQA9 + WQA9A(NAL) * WQAPC(L)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA9A(NAL) = (ALGAES(NAL).WQFPDB*WQBM(L,NAL) + ALGAES(NAL).WQFPDP*WQPR(L,NAL)) * WQO(L,K,19+NAL)
              WQA9 = WQA9 + WQA9A(NAL) * WQAPC(L) * ALGAES(NAL).WQAPCM 
            endif
          enddo

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR9 = WQWPSL(L,K,IDOP) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,IDOP) + DTWQ*WQR9 + WQF9*WQVO(L,K,IDOP) + DTWQO2*(WQA9 + WQKRPP(L)*WQO(L,K,IROP) + WQKLPP(L)*WQO(L,K,ILOP))
           
          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDOPZ(L,K)
        enddo

        if( K == KC )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP         ATM WET DEP       VOLUME  
            WQR9 = (WQWDSL(L,KC,IDOP) + WQATML(L,KC,IDOP))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR9  
          enddo
        endif

        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IDOP) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,IDOP) = WQVO(L,K,IDOP) + WQV(L,K,IDOP)
        enddo  
      endif  

      ! ****  PARAM 07  P4D - total phosphate (PO4t)
      if( ISKINETICS(7) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          
          ! *** Algal metabolism and predation source
          WQKK(L) = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              WQA10A(NAL) = (ALGAES(NAL).WQFPIB*WQBM(L,NAL) + ALGAES(NAL).WQFPIP*WQPR(L,NAL) - WQPA(L,NAL))*WQO(L,K,19+NAL)
              WQKK(L) = WQKK(L) + WQA10A(NAL)*WQAPC(L)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA10A(NAL) = (ALGAES(NAL).WQFPIB*WQBM(L,NAL) + ALGAES(NAL).WQFPIP*WQPR(L,NAL) - WQPA(L,NAL))*WQO(L,K,19+NAL)
              WQKK(L) = WQKK(L) + WQA10A(NAL)*WQAPC(L)*ALGAES(NAL).WQAPCM  
            endif
          enddo 

          ! ***      PT_SRC_LOADS      VOLUME  
          WQRR(L) = WQWPSL(L,K,IP4D) * VOLWQ(L) * PSMLMULTIPLIER
           
          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + SPO4Z(L,K)
        enddo  

        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          if( K == KSZ(L) )then                       ! *** Note - CAN ADD A PER CELL ARRAY OF 1 TO KC WHERE K = KSZ IS 1.0 AND ALL OTHER LAYERS = 0.0
            WQRR(L) = WQRR(L) + WQBFPO4D(L)*HPKI(L,K) ! *** Add in Benthic Flux
          endif
        enddo  

        if( K == KC )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***      ATM DRY DEP       ATM WET DEP        VOLUME  
            WQR10 = (WQWDSL(L,KC,IP4D) + WQATML(L,KC,IP4D))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR10  
          enddo
        endif  

        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          WQRR(L) = WQVO(L,K,IP4D) + DTWQ*WQRR(L) + WQH10(L)*WQVO(L,K,IP4D) + DTWQO2*(WQKK(L) + WQKDOP(L)*WQO(L,K,IDOP))
        enddo  
  
        if( K /= KC )then                              ! *** Note - no settling if IWQSRP = 0
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQT10(L)*WQVO(L,K+1,IP4D)
          enddo  
        endif 
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQKKL = 1.0 / (1.0 - WQH10(L)) 
          WQV(L,K,IP4D) = max(WQRR(L)*WQKKL,0.0)
          WQO(L,K,IP4D) = WQVO(L,K,IP4D) + WQV(L,K,IP4D)
        enddo  
      endif  

      ! ****  PARAM 08  RON - refractory particulate organic nitrogen
      if( ISKINETICS(8) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS     SETTLING
          WQI11 = - (WQKRPN(L) + WQRPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQI11)
          
          ! *** Algae metabolism and predation
          WQA11 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA11A(NAL) = (ALGAES(NAL).WQFNRB*WQBM(L,NAL) + ALGAES(NAL).WQFNRP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)  
              WQA11 = WQA11 + WQA11A(NAL)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA11A(NAL) = (ALGAES(NAL).WQFNRB*WQBM(L,NAL) + ALGAES(NAL).WQFNRP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)   
              WQA11 = WQA11 + WQA11A(NAL)
            endif
          enddo

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR11 = WQWPSL(L,K,IRON) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,IRON) + DTWQ*WQR11 + DTWQO2*WQA11 + WQI11*WQVO(L,K,IRON)
          
          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SRPONZ(L,K)
        enddo 

        if( K /= KC )then  
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQVO(L,K+1,IRON)
          enddo  
        else   ! K == KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
            WQR11 = (WQWDSL(L,KC,IRON) + WQATML(L,KC,IRON))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR11  
          enddo
        endif
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IRON) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,IRON) = WQVO(L,K,IRON) + WQV(L,K,IRON)
        enddo  
      endif  

      ! ****  PARAM 09  LON - labile particulate organic nitrogen
      if( ISKINETICS(9) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS     SETTLING
          WQJ12 = - (WQKLPN(L) + WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQJ12)
          
          ! *** Algae metabolism and predation
          WQA12 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA12A(NAL) = (ALGAES(NAL).WQFNLB*WQBM(L,NAL) + ALGAES(NAL).WQFNLP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)  
              WQA12 = WQA12 + WQA12A(NAL)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA12A(NAL) = (ALGAES(NAL).WQFNLB*WQBM(L,NAL) + ALGAES(NAL).WQFNLP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)  
              WQA12 = WQA12 + WQA12A(NAL)
            endif
          enddo

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR12 = WQWPSL(L,K,ILON) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQV(L,K,ILON) + DTWQ*WQR12 + DTWQO2*WQA12 + WQJ12*WQVO(L,K,ILON)
                     

          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SLPONZ(L,K)
        enddo  
        if( K /= KC )then  
          ! *** Add in settling from above
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQVO(L,K+1,ILON)
          enddo  
        else   ! K == KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP       ATM WET DEP     VOLUME  
            WQR12 = (WQWDSL(L,KC,ILON) + WQATML(L,KC,ILON))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR12  
          enddo
        endif  
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,ILON) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,ILON) = WQVO(L,K,ILON) + WQV(L,K,ILON)
        enddo  
      endif  

      ! ****  PARAM 10  DON - dissolved organic nitrogen
      if( ISKINETICS(10) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQF13 = - DTWQO2*WQKDON(L)  
          WQKK(L) = 1.0 / (1.0 - WQF13)
          
          ! *** Algal metabolism and predation source
          WQA13 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA13A(NAL) = (ALGAES(NAL).WQFNDB*WQBM(L,NAL) + ALGAES(NAL).WQFNDP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)  
              WQA13 = WQA13 + WQA13A(NAL) 
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA13A(NAL) = (ALGAES(NAL).WQFNDB*WQBM(L,NAL) + ALGAES(NAL).WQFNDP*WQPR(L,NAL))*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)  
              WQA13 = WQA13 + WQA13A(NAL) 
            endif
          enddo
          
          ! ***    PT_SRC_LOADS    VOLUME  
          WQR13 = WQWPSL(L,K,IDON) * VOLWQ(L) * PSMLMULTIPLIER
          WQRR(L) = WQVO(L,K,IDON) + DTWQ*WQR13 + WQF13*WQVO(L,K,IDON) + DTWQO2*(WQA13 + WQKRPN(L)*WQO(L,K,IRON)  + WQKLPN(L)*WQO(L,K,ILON))
           
          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDONZ(L,K)
        enddo
        
        if( K == KC )then                          ! *** Note - add a check to see of the deposition concentrations are > 0.  Do this for all constituents
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
            WQR13 = (WQWDSL(L,KC,IDON) + WQATML(L,KC,IDON))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR13  
          enddo
        endif  

        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQV(L,K,IDON) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,IDON) = WQVO(L,K,IDON) + WQV(L,K,IDON)
        enddo  
      endif  
      
      ! ****  PARAM 11  NHX - ammonia nitrogen
      if( ISKINETICS(11) == 1 )then 
        if( IWQPSL < 2 )then
          do LP = 1,LLWET(K,ND)                      ! *** Note - SKIP IF IWQPSL = 2
            L = LKWET(LP,K,ND)  
            ! ***      PT_SRC_LOADS     VOLUME  
            WQRR(L) = WQWPSL(L,K,INHX) * VOLWQ(L) * PSMLMULTIPLIER
          enddo
        else
          ! *** Zero the loadings array
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = 0.0
          enddo   
        endif

        do LP = 1,LLWET(K,ND)                        ! *** Note - CAN ADD A PER CELL ARRAY OF 1 TO KC WHERE K = KSZ IS 1.0 AND ALL OTHER LAYERS = 0.0
          L = LKWET(LP,K,ND)  
          if( K == KSZ(L) )then
            WQRR(L) = WQRR(L) + WQBFNH4(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          endif  
        enddo  

        if( K == KC )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
            WQR14 = (WQWDSL(L,KC,INHX) + WQATML(L,KC,INHX))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR14  
          enddo
        endif  

        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQF14 = - DTWQO2*WQNIT(L)  
          WQKK(L) = 1.0 / (1.0 - WQF14)
          
          ! *** Algal metabolism and predation source
          WQA14 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA14A(NAL) = ALGAES(NAL).WQFNIB*WQBM(L,NAL) + ALGAES(NAL).WQFNIP*WQPR(L,NAL) - WQPN(L,1)*WQPA(L,NAL)  
              WQA14 = WQA14 + WQA14A(NAL)*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA14A(NAL) = ALGAES(NAL).WQFNIB*WQBM(L,NAL) + ALGAES(NAL).WQFNIP*WQPR(L,NAL) - WQPN(L,1)*WQPA(L,NAL)  
              WQA14 = WQA14 + WQA14A(NAL)*ALGAES(NAL).WQANCA*WQO(L,K,19+NAL)
            endif
          enddo
 
          WQRR(L) = WQV(L,K,INHX) + DTWQ*WQRR(L) + WQF14*WQVO(L,K,INHX) + DTWQO2*(WQA14 + WQKDON(L)*WQO(L,K,IDON))

          ! *** Coupling to zooplankton - DKT
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SNH4Z(L,K)
          
          WQV(L,K,INHX) = max(WQRR(L)*WQKK(L),0.0)
          WQO(L,K,INHX) = WQVO(L,K,INHX) + WQV(L,K,INHX)
        enddo  
      endif  

      ! ****  PARAM 12  NOX - nitrate nitrogen
      if( ISKINETICS(12) == 1 )then 
        if( IWQPSL < 2 )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***      PT_SRC_LOADS      VOLUME  
            WQRR(L) = WQWPSL(L,K,INOX) * VOLWQ(L) * PSMLMULTIPLIER
          enddo  
        else
          ! *** Zero the loadings array
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = 0.0
          enddo   
        endif
    
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          if( K == KSZ(L) )then
            WQRR(L) = WQRR(L) + WQBFNO3(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          endif  
        enddo  
    
        if( K == KC )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP         ATM WET DEP       VOLUME  
            WQR15 = (WQWDSL(L,KC,INOX) + WQATML(L,KC,INOX))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR15  
          enddo
        endif  
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          ! *** Algal predation source
          WQA15 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE  )then
              ! *** Phytoplankton
              WQA15A(NAL) = (WQPN(L,NAL) - 1.0)*WQPA(L,NAL) * ALGAES(NAL).WQANCA * WQO(L,K,19+NAL)  
              WQA15 = WQA15 + WQA15A(NAL)
            endif
            if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then
              ! *** Macrophytes and periphyton
              WQA15A(NAL) = (WQPN(L,NAL) - 1.0)*WQPA(L,NAL) * ALGAES(NAL).WQANCA * WQO(L,K,19+NAL)  
              WQA15 = WQA15 + WQA15A(NAL)
            endif
          enddo
  
          WQB15 = WQV(L,K,INOX) + DTWQ*WQRR(L) + DTWQO2*( WQA15 - WQANDC*WQDENIT(L)*WQO(L,K,IDOC) + WQNIT(L)*WQO(L,K,INHX))

          WQV(L,K,INOX) = max(WQB15,0.0)
          WQO(L,K,INOX) = WQVO(L,K,INOX) + WQV(L,K,INOX)
        enddo  
      endif  

      ! ****  PARAM 13  SUU - particulate biogenic silica
      if( ISKINETICS(13) == 1 )then  
        if( IWQSI == 1 )then  
          ! *** Algal metabolism and predation source
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            
            WQA16D = 0.0
            WQM16 = 0.0
            do NAL = 1,NALGAE
              if( ALGAES(NAL).ISILICA > 0 )then
                WQA16D = WQA16D + (ALGAES(NAL).WQFSPD*WQBM(L,NAL) + ALGAES(NAL).WQFSPP*WQPR(L,NAL)) * ALGAES(NAL).WQASC * WQO(L,K,19+NAL)
                WQM16 =  WQM16  - (WQKSUA(IWQT(L)) + WQBSETL(L,2,NAL) + WQBSETL(L,3,NAL)) * DTWQO2  
              endif
            enddo
            WQKK(L) = 1.0 / (1.0 - WQM16)
                        
            ! ***    PT_SRC_LOADS      VOLUME  
            WQR16 = WQWPSL(L,K,ISUU) * VOLWQ(L) * PSMLMULTIPLIER
            WQRR(L) = WQV(L,K,ISUU) + DTWQ*WQR16 + DTWQO2*WQA16D + WQM16*WQO(L,K,ISUU)
          enddo 
             
          ! *** Coupling to zooplankton 
          if( IWQZPL > 0 )then
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              WQRR(L) = WQRR(L) + DTWQ*SSUZ(L,K)
            enddo
          endif
          
          ! *** Vertical fluxes from various sources
          if( K /= KC )then  
            do NAL = 1,NALGAE
              if( ALGAES(NAL).ISILICA > 0 )then
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  WQRR(L) = WQRR(L) + DTWQO2*WQBSETL(L,1,NAL)*WQO(L,K+1,ISUU)    ! *** Vertical migration into the current layer from layer above
                enddo
              endif
            enddo  
          else
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***     ATM DRY DEP         ATM WET DEP       VOLUME  
              WQR16 = (WQWDSL(L,KC,ISUU) + WQATML(L,KC,ISUU))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR16  
            enddo
          endif 

          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQV(L,K,ISUU) = max(WQRR(L)*WQKK(L),0.0)
            WQO(L,K,ISUU) = WQVO(L,K,ISUU) + WQV(L,K,ISUU)
          enddo  
        endif  
      endif  

      ! ****  PARAM 14  SAA - dissolved available silica
      if( ISKINETICS(14) == 1 )then  
        if( IWQSI == 1 )then  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            WQKK(L) = 0.0
          enddo
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISILICA > 0 )then
              do LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)
                WQKK(L) = WQKK(L) + (ALGAES(NAL).WQFSID*WQBM(L,NAL) + ALGAES(NAL).WQFSIP*WQPR(L,NAL) - WQPA(L,NAL)) * ALGAES(NAL).WQASC * WQO(L,K,19+NAL)
              enddo
            endif
          enddo
          
          ! *** Add mass loading BCs
          if( IWQPSL /= 2 )then
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              ! ***      PT_SRC_LOADS      VOLUME  
              WQRR(L) = WQWPSL(L,K,ISAA) * VOLWQ(L) * PSMLMULTIPLIER
            enddo
          else
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              WQRR(L) = 0.0
            enddo
          endif
      
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            if( K == KSZ(L) )then
              WQRR(L) = WQRR(L) + WQBFSAD(L)*HPKI(L,K)   ! *** Add in Benthic Flux
            endif  
          enddo  
          
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = WQV(L,K,ISAA) + DTWQ*WQRR(L) + WQN17(L)*WQVO(L,K,ISAA) + DTWQO2*( WQKK(L) + WQKSUA(IWQT(L))*WQO(L,K,ISUU))
          enddo  
  
          ! *** Coupling to zooplankton
          if( IWQZPL > 0 )then
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQ*SSAZ(L,K)
            enddo  
          endif
           
           if( K /= KC )then  
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQT17(L)*WQVO(L,K+1,ISAA)  
            enddo  
          else
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***     ATM DRY DEP        ATM WET DEP        VOLUME  
              WQR17 = (WQWDSL(L,KC,ISAA) + WQATML(L,KC,ISAA))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR17  
            enddo
          endif  
          
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQR17 = 1.0 / (1.0 - WQN17(L))
            WQV(L,K,ISAA) = max(WQRR(L)*WQR17,0.0)
            WQO(L,K,ISAA) = WQVO(L,K,ISAA) + WQV(L,K,ISAA)
          enddo  
        endif  
      endif  

      ! ****  PARAM 15  COD - chemical oxygen demand
      if( ISKINETICS(15) == 1 )then  
        if( IWQPSL /= 2 )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***      PT_SRC_LOADS      VOLUME  
            WQRR(L) = WQWPSL(L,K,ICOD) * VOLWQ(L) * PSMLMULTIPLIER
          enddo  
        else
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQRR(L) = 0.0
          enddo  
        endif
        
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          if( K == KSZ(L) )then
            WQRR(L) = WQRR(L) + WQBFCOD(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          endif  
        enddo  
        
        if( K == KC )then
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP         ATM WET DEP       VOLUME  
            WQR18 = (WQWDSL(L,KC,ICOD) + WQATML(L,KC,ICOD))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR18  
          enddo
        endif  
         
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQR18 = 1.0 / (1.0 - WQO18(L))  
          WQRR(L) = WQV(L,K,ICOD) + DTWQ*WQRR(L) + WQO18(L)*WQV(L,K,ICOD)  
          WQV(L,K,ICOD) = max(WQRR(L)*WQR18,0.0)
          WQO(L,K,ICOD) = WQVO(L,K,ICOD) + WQV(L,K,ICOD)  
        enddo  
      endif

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

      if( ISKINETICS(16) == 1 )then 
        
        ! *** Point Source Loading 
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRR(L) = WQWPSL(L,K,IDOX) * VOLWQ(L) * PSMLMULTIPLIER
        enddo
  
        ! *** Handle Surface Processes
        if( K == KC )then  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L)) 
            ! ***                  ATM DRY DEP        ATM WET DEP       VOLUME  
            WQRR(L) = WQRR(L) + (WQWDSL(L,KC,IDOX) + WQATML(L,KC,IDOX))*VOLWQ(L)
            ! *** Reaeration - Atm to water flux  (Offset later in O2 Mass Balance by WQRea term)
            WQRR(L) = WQRR(L) + WQKRDOS(L)
          enddo
        else  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0
          enddo
        endif 

        ! *** Bottom Processes 
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          if( K == KSZ(L) )then
            TEMFAC = STEMFAC**(TEM(L,KSZ(L)) - 20.)
            WQRR(L) = WQRR(L) + TEMFAC*WQBFO2(L)*HPKI(L,K)   ! *** Add in Benthic Flux
          endif
        enddo
  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          DTWQxH = DTWQ*HPK(L,K)  
          DTWQxH2= DTWQO2*HPK(L,K)

          ! *** Photosynthesis
          if( WQI0  <=  0.001 )then
            do NAL = 1, NALGAE
              WQTTAA(NAL) = 0.0
            enddo
          else
            do NAL = 1,NALGAE
              if( ALGAES(NAL).ISMOBILE )then
                ! *** Phytoplankton
                WQTTAA(NAL) = (1.3 - 0.3*WQPN(L,NAL)) * WQPA(L,NAL) 
              endif
              if( .not. ALGAES(NAL).ISMOBILE .and. (K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L)) )then  
                ! *** Macrophytes and periphyton
                WQTTAA(NAL) = (1.3 - 0.3*WQPN(L,NAL)) * WQPA(L,NAL)
              endif
            enddo               
          endif
          
          ! *** Respiration and final net O2 rate
          WQA19 = 0.
          IZ = IWQZMAP(L,K)  
          do NAL = 1,NALGAE
            ! *** Respiration 
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              WQRes = (1.0 - ALGAES(NAL).WQFCDB)*O2WQ(L)*WQBM(L,NAL)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)       
            endif
            if( .not. ALGAES(NAL).ISMOBILE )then
              if( K >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L) )then  
                ! *** Macrophytes and periphyton
                WQRes = (1.0 - ALGAES(NAL).WQFCDB)*O2WQ(L)*WQBM(L,NAL)/(ALGAES(NAL).WQKHRA(IZ) + O2WQ(L) + 1.E-18)
              else
                WQRes = 0.
              endif
            endif
            
            ! *** Net algal O2 production due to photosynthesis and respiration
            WQA19A(NAL) = WQTTAA(NAL)*ALGAES(NAL).WQAOCRP - WQRes*ALGAES(NAL).WQAOCRR
            WQA19A(NAL) = WQA19A(NAL)*WQO(L,K,19+NAL)                                      ! *** Net O2 production by group
            WQA19 = WQA19 + WQA19A(NAL)                                                    ! *** Total O2 net production
          enddo

          ! WQA19                                    ! *** Total Net Respiration/Photosynthesis
          WQSUM = DTWQ*WQRR(L)                       ! *** Sum of Loadings/Demands
          WQRea = DTWQO2*WQP19(L)*WQV(L,K,IDOX)      ! *** Reaeration - Offsetting Flux (WQP19 is negative)
          WQDOC = WQAOCR*WQKHR(L)*WQO(L,K,IDOC)      ! *** DOC
          WQNH3 = WQAONT*WQNIT(L)*WQO(L,K,INHX)      ! *** Ammonia
          WQCOD = WQO18(L)*WQO(L,K,ICOD)             ! *** COD
          WQRR(L) = WQV(L,K,IDOX) + WQSUM  + WQCOD  + WQRea + DTWQO2*(WQA19 - WQDOC - WQNH3)
          
          if( IWQZPL > 0 ) WQRR(L) = WQRR(L) + DTWQ*SDOZ(L,K)    ! *** Coupling to zooplankton

          WQV(L,K,IDOX) = max(WQRR(L)*WQKK(L),0.0)

        enddo
        
      endif  ! *** END OF ISKINETICS(16)

      ! ****  PARAM 17  TAM - total active metal
      if( ISKINETICS(17) == 1 )then  
        if( IWQSRP == 1 )then  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQR20(L) = WQWPSL(L,K,17)*VOLWQ(L) * PSMLMULTIPLIER + (WQV(L,K,ITAM) - WQTAMP(L,K))*WQWSSET(L,1)  
          enddo

          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            if(K == KSZ(L)) WQR20(L) = WQR20(L) + WQTDTAM(IWQT(L))*HPKI(L,K)/(WQKHBMF + O2WQ(L) + 1.E-18)
          enddo

          if( K /= KC )then
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQR20(L) = WQR20(L) + (WQV(L,K+1,ITAM) - WQTAMP(L,K+1)) * WQWSSET(L,2)  
            enddo 
          else     ! K == KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQR20(L) = WQR20(L) + (WQWDSL(L,KC,ITAM) + WQATML(L,KC,ITAM))*VOLWQ(L)
            enddo
          endif 

          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQT20 = - DTWQ*WQWSSET(L,1)    ! *** DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQT20)  
            WQRR(L) = WQV(L,K,ITAM) + DTWQ*WQR20(L) + WQT20*WQV(L,K,ITAM)  
          enddo  
          
          if( K /= KC )then  
            ! *** Add in settling from above
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQWSSET(L,2)*WQVO(L,K+1,ITAM)
            enddo  
          endif  
          
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQV(L,K,ITAM) = max(WQRR(L)*WQKK(L),0.0)
            WQO(L,K,ITAM) = WQVO(L,K,ITAM) + WQV(L,K,ITAM)
          enddo  
        endif  
      endif  

      ! ****  PARAM 18  FCB - fecal coliform bacteria
      if( ISKINETICS(18) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQKK(L) = WQTD2FCB(IWQT(L))  
          
          ! ***   PT SRC LOADS     VOLUME  
          WQR21 = WQWPSL(L,K,IFCB)*VOLWQ(L) * PSMLMULTIPLIER
          
          if( K == KC )then
            ! ***            ATM DRY DEP        ATM WET DEP       VOLUME  
            WQR21 = WQR21 + (WQWDSL(L,K,IFCB) + WQATML(L,KC,IFCB))*VOLWQ(L)
          endif
      
          WQRR(L) = WQV(L,K,IFCB)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21
          WQV(L,K,IFCB) = max(WQRR(L)*WQKK(L),0.0)
        enddo  
      endif

      ! ***************************************************************************
      ! *** SHELLFISH KINETICS *** shellfish 2
      if( ISFFARM > 0 .and. NSF > 0 )then
        do L = LF,LL
          do NQ = 1,NSF
            call SHELLFISH_GROWTH(NQ, SF(NQ), L, K)
          enddo
        enddo
        call SHELLFISH_HARVEST(LF,LL)
      endif

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
      if( ISKINETICS(19) == 1 )then  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQRR(L) = WQWPSL(L,K,ICO2) * VOLWQ(L) * PSMLMULTIPLIER  
          TMP22 = WQRR(L)*DTWQ*HPK(L,K) 
        enddo  
  
        ! *** Handle Surface Processes
        if( K == KC )then  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP22(L))
            
            ! ***                  ATM DRY DEP     ATM WET DEP     VOLUME  
            WQRR(L) = WQRR(L) + (WQWDSL(L,KC,ICO2) + WQATML(L,KC,ICO2)*VOLWQ(L))  
            
            ! *** Reaeration
            WQRR(L) = WQRR(L) + WQKRCDOS(L) 
          enddo
        else
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            WQKK(L) = 1.0
          enddo
        endif   
  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          DTWQxH  = DTWQ*HPK(L,K)  
          DTWQxH2 = DTWQO2*HPK(L,K)
          if( WQI0  <=  0.001 )then
            do NAL = 1,NALGAE
              if( ALGAES(NAL).ISMOBILE )then
                WQTTAA(NAL) = 0.0
              endif
            enddo
          else
            do NAL = 1,NALGAE
              ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL 
              if( ALGAES(NAL).ISMOBILE )then
                WQTTAA(NAL) = (1.3 - 0.3*WQPN(L,NAL)) * WQPA(L,NAL) 
              endif
            enddo
          endif
          
          ! *** RESPIRATION OF TOTAL CHLOROPHYLL - general algal group
          WQA22 = 0.
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISMOBILE )then
              ! *** Phytoplankton
              XMRM = (1.0 - ALGAES(NAL).WQFCDB)*O2WQ(L)*WQBM(L,NAL)/(ALGAES(NAL).WQKHRA(1) + O2WQ(L) + 1.E-18)  
              WQA22A(NAL) = WQTTAA(NAL) - XMRM
              WQA22 = WQA22 + 3.67*WQA22A(NAL)*WQO(L,K,19+NAL)  !VB 3.67 CONVERTS g CARBON TO g CO2
            endif
          enddo

          ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
          ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
          ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
          ! *** CO2 Mass Balance
          ! WQA22                                                 ! *** Total Net Respiration/Photosynthesis
          WQCDSUM = DTWQ*WQRR(L)                                  ! *** Sum of Loadings/Demands
          WQCDRea = WQP22(L)*WQV(L,K,ICO2)                        ! *** Reaeration
          WQCDDOC = (WQKHR(L) + WQDENIT(L))*WQO(L,K,IDOC)*3.67    ! *** DOC FROM HYDROLYSIS AND DENITRIFICATION 3.67 CONVERTS G CARBON TO G CO2    

          WQRR(L) = WQV(L,K,ICO2) + WQCDSUM + DTWQO2*(-WQA22 + WQCDDOC + WQCDRea)
          WQV(L,K,ICO2) = max(WQRR(L)*WQKK(L),0.0)
        enddo  
      endif  
    enddo  ! *** END OF THE KC LOOP

    ! ***************************************************************************
    ! *** SUM OF SHELLFISH FOR ALL LAYERS *** shellfish 4
    if( ISFFARM > 0 .and. NSF > 0 )then
      do L = LF,LL
        call SHELLFISH_LAYERSUM(L)
      enddo
    endif

    ! *** PHOSPHATE AND SILICA SORPTION OPTIONS  
    if( IWQSRP == 1 )then  
      ! *** Sorption Option: TAM
      do K = 1,KC  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          O2WQ(L) = max(WQV(L,K,IDOX), 0.0)  
          WQTAMD = min( WQTAMDMX*EXP(-WQKDOTAM*O2WQ(L)), WQV(L,K,ITAM) )  
          WQTAMP(L,K) = WQV(L,K,ITAM) - WQTAMD  
          WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*WQTAMP(L,K))  
          WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*WQTAMP(L,K))  
        enddo  
      enddo  
    elseif( IWQSRP == 2 .and. ISTRAN(6) > 0 )then
      ! *** Sorption Option: Sediments
      do K = 1,KC  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQPO4D(L,K) = WQV(L,K,IP4D) / (1.0 + WQKPO4P*SEDT(L,K))  
          WQSAD(L,K)  = WQV(L,K,ISAA) / (1.0 + WQKSAP*SEDT(L,K))  
        enddo  
      enddo  
    else  
      do K = 1,KC  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          WQPO4D(L,K) = WQV(L,K,IP4D)  
          WQSAD(L,K)  = WQV(L,K,ISAA)  
        enddo  
      enddo  
    endif  

    ! ***************************************************************************
    ! *** Coupling to sediment diagenesis model  
    if( IWQBEN == 1 )then
    
      ! *** Depositional flux of biota and particulate organic nutrients   ToDo - eliminate SCB and zero fluxes wden resetting OBC's
      do LP = LF,LL
        L = LWET(LP)
        K = KSZ(L)
        
        do NAL = 1,NALGAE
          if( ALGAES(NAL).ISMOBILE )then
            if( ALGAES(NAL).SETTLING(L,K) > 0.0 )then
              ! *** Depositional flux of phytoplankton
              WQDFB(L,NAL) = SCB(L)*ALGAES(NAL).SETTLING(L,K)*WQV(L,K,19+NAL) 
            else
              WQDFB(L,NAL) = 0.0
            endif
          else
            ! *** Depositional flux of macrophytes
            WQDFB(L,NAL) = SCB(L)*ALGAES(NAL).SETTLING(L,K)*WQV(L,K,19+NAL)*HPKI(L,K)
          endif
        enddo       
        
        IMWQZ = IWQZMAP(L,KSZ(L))
        WQDFRC(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,K,IROC)   !< Depositional flux of RPOC 
        WQDFLC(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,K,ILOC)   !< Depositional flux of LPOC   
        WQDFRP(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,K,IROP)   !< Depositional flux of RPOP
        WQDFLP(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,K,ILOP)   !< Depositional flux of LPOP   
        WQDFRN(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,K,IRON)   !< Depositional flux of RPON
        WQDFLN(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,K,ILON)   !< Depositional flux of LPON   
      enddo  
        
      ! *** Depositional flux of particulate silica
      if( IWQSI == 1 )then
        do LP = LF,LL
          L = LWET(LP)
          K = KSZ(L)
          
          WQDFSI(L) = 0.0
          do NAL = 1,NALGAE
            if( ALGAES(NAL).ISILICA > 0 .and. ALGAES(NAL).SETTLING(L,K) > 0.0 )then
              WQDFSI(L) = WQDFSI(L) + SCB(L)*ALGAES(NAL).SETTLING(L,K)*WQV(L,K,ISUU)  
            endif
          enddo    
        enddo
      endif      

      ! *** Depositional flux of sorbed phosphate and dissolved silica
      if( IWQSRP == 1 )then  
        ! *** Sorbed to TAM
        do LP = LF,LL
          L = LWET(LP)
          IMWQZ = IWQZMAP(L,KSZ(L))  
          WQDFLP(L) = SCB(L)*(WQDFLP(L) + WQWSS(IMWQZ)*( WQV(L,KSZ(L),IP4D) - WQPO4D(L,1)) )       ! *** Add PO4 flux to LOP flux
          if( IWQSI == 1 ) WQDFSI(L) = (WQDFSI(L) + WQWSS(IMWQZ)*(WQV(L,KSZ(L),ISAA) - WQSAD(L,1)))  
        enddo  
      elseif( IWQSRP == 2 .and. ISTRAN(6) > 0 )then  
        ! *** Sorbed to cohesive class 1
        do LP = LF,LL
          L = LWET(LP)
          WQDFLP(L) = SCB(L)*(WQDFLP(L) + WSEDO(NS)*( WQV(L,KSZ(L),IP4D) - WQPO4D(L,1)) )          ! *** Add PO4 flux to LOP flux
          if( IWQSI == 1 ) WQDFSI(L) = (WQDFSI(L) + WSEDO(NS)*(WQV(L,KSZ(L),ISAA) - WQSAD(L,1)))  
        enddo  
      endif  
    endif  
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO
  
  ! *** ZERO THIN LAYERS WHEN DELTA T WQ IS LARGE
  if( KC > 1 .and. ISDRY > 0 )then
    IFLAG = 0
    !$OMP PARALLEL DO PRIVATE(ND,LF,LL,LP,L,NW,K) REDUCTION(+:IFLAG)
    do ND = 1,NDM 
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)
      do LP = LF,LL
        L = LWET(LP)
        if( HP(L) < 2.*HDRY )then
          do NW = 1,NWQV
            if( ISKINETICS(NW) > 0 )then
              do K = KSZ(L),KC
                if( WQV(L,K,NW) < 0.0 )then
                  WQV(L,K,NW) = 0.0
                  IFLAG = IFLAG + 1
                endif
              enddo
            endif
          enddo
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO
    if( IFLAG > 0 )then
      open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
      write(mpi_efdc_out_unit,'(A,F12.4,I10)')' WARNING!  NEGATIVE WQ CONCENTRATIONS: TIMEDAY, # OF VALUES < 0.0:',TIMEDAY,IFLAG
      close(mpi_efdc_out_unit)
    endif
  endif
  
  ! *** RESTORE PREVIOUS WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  do NQ = 0,NWQVM
    do K = 1,KC
      do IOBC = 1,NBCSOP  
        L = LOBCS(IOBC)
        WQV(L,K,NQ) = WQOLD(IOBC,K,NQ)
      enddo
    enddo  
  enddo  
  
  ! ***************************************************************************
  ! *** DIURNAL DO ANALYSIS  
  if( WQHRAVG >  0.0 )then 
    NDDOCNT = NDDOCNT + 1
    do K = 1,KC  
      do L = 2,LA
        DDOAVG(L,K) = DDOAVG(L,K) + WQV(L,K,IDOX)*DTWQ
        DDOMAX(L,K) = max(DDOMAX(L,K),WQV(L,K,IDOX))  
        DDOMIN(L,K) = min(DDOMIN(L,K),WQV(L,K,IDOX))  
      enddo  
    enddo  

    if( TIMEDAY >= WQHRAVG )then  
      NDDOCNT = 0  

      write(mpi_efdc_out_unit, '(a)') '*************************************************************'
      write(mpi_efdc_out_unit, '(a,f10.2)') 'Water Quality D.O. Analysis : ',TIMEDAY
      call WriteBreak(mpi_efdc_out_unit)

      TVAL1 = 1./((TIMEDAY-AVGLAST)*86400.)
      do L = 2,LA
        DDOAVG(L,K) = DDOAVG(L,K)*TVAL1
        write(mpi_efdc_out_unit,1112)  Map2Global(L).IG, Map2Global(L).JG, (K, DDOMIN(L,K), DDOAVG(L,K), DDOMAX(L,K),K = 1,KC)  
      enddo  
      
      do K = 1,KC  
        do L = 2,LA  
          DDOAVG(L,K) = 0.0
          DDOMAX(L,K) = -1.E6  
          DDOMIN(L,K) = 1.E6  
        enddo  
      enddo  
      AVGLAST = TIMEDAY
      AVGNEXT = TIMEDAY + WQHRAVG
    endif  
  endif  

  ! ***************************************************************************
  ! *** LIGHT EXTINCTION ANALYSIS  
  if( NDLTAVG >= 1 )then 
    if( process_id == master_id )then
      open(1,FILE = OUTDIR//'LIGHT.OUT',POSITION = 'APPEND')  
    endif
    NDLTCNT = NDLTCNT+1  
    NSTPTMP = NDLTAVG*NTSPTC/2  
    RMULTMP = 1./FLOAT(NSTPTMP)  
    do K = 1,KC  
      do L = 2,LA 
        if( ISTRAN(6) > 0  )then
          RLIGHT1 = WQKEB(IWQZMAP(L,K))+WQKETSS*SEDT(L,K) 
        else
          RLIGHT1 = WQKEB(IWQZMAP(L,K))
        endif
        XMRM = WQKECHL*WQCHL(L,K)  
        if( WQKECHL  < 0.0 )then  
          XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)  
        endif  
        RLIGHT2 = XMRM  
        RLIGHTT(L,K) = RLIGHTT(L,K)+RLIGHT1  
        RLIGHTC(L,K) = RLIGHTC(L,K)+RLIGHT1+RLIGHT2  
      enddo  
    enddo  
    if( NDLTCNT == NSTPTMP )then  
      NDLTCNT = 0  
      TIME = TIMESEC/TCON  

      do K = 1,KC  
        do L = 2,LA
          RLIGHTT(L,K) = RMULTMP*RLIGHTT(L,K)  
          RLIGHTC(L,K) = RMULTMP*RLIGHTC(L,K)  
        enddo  
      enddo 
      
      if( process_id == master_id )then
        write(1,1111)NITER,TIME  
        do L = 2,LA  
          write(1,1113)IL(L),JL(L),(RLIGHTT(L,K),K = 1,KC),(RLIGHTC(L,K),K = 1,KC)  
        enddo  
      endif
      do K = 1,KC  
        do L = 2,LA  
          RLIGHTT(L,K) = 0.  
          RLIGHTC(L,K) = 0.  
        enddo  
      enddo  
    endif  
    close(1)  
  endif  

  ! ***************************************************************************
  ! INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:  
  TIMTMP = TIMESEC/TCON  
  
  1111 FORMAT(I12,F10.4)  
  1112 FORMAT(2I7,200(I4,3F7.2))  
  1113 FORMAT(2I7,12E12.4)  
  1414 FORMAT(I12,11E12.4)  

  return  
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

  PRINT *,'WQSKE2 HAS BEEN REMOVED.  NEEDS CODE IF IWQSKE = 2 IMPLEMENTED'
  
  return  
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

  PRINT *,'WQSKE3 HAS BEEN REMOVED.  NEEDS CODE IF IWQSKE = 3 IMPLEMENTED'
  
  return  
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

    PRINT *,'WQSKE4 HAS BEEN REMOVED.  NEEDS CODE IF IWQSKE = 4 IMPLEMENTED'
  
    return  
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

      implicit none

  !  write(2,83)': FREQUENCY OF DIURNAL DO OUTPUT (IN DT UNIT)  = ', IWQDIUDT
  !  write(2,83)'* IWQDT (DTWQ(D) = DT(S)*IWQDT/86400)        = ', IWQDT
  !
  !  write(2,80)'* FULL VERSION WITH 21 VARIABLES IS ACTIVATED     '
  !  if( IWQBEN == 1 )then
  !    write(2,80)'* SEDIMENT PROCESS MODEL IS ACTIVATED             '
  !  elseif( IWQBEN == 0 )then
  !    write(2,80)'* SPATIALLY/TEMPORALLY CONSTANT BENTHIC FLUX IS SPECIFIED   '
  !  elseif( IWQBEN == 2 )then
  !    write(2,80)'* SPATIALLY AND/OR TEMPORALLY-VARYING BENTHIC FLUX SPECIFIED'
  !  else
  !    call STOPP('** ERROR!!! INVALID IWQBEN VALUE **')
  !  endif
  !
  !  if( IWQSI == 1 )then
  !    write(2,80)'* SILICA STATE VARIABLES (SU & SA) ARE MODELED    '
  !  else
  !    write(2,80)'* NO SILICA (SU & SU) LIMITATION                  '
  !  endif
  !
  !  if( IWQFCB == 1 )then
  !    write(2,80)'* FCB (FECAL COLIFORM BACTERIA) IS MODELED        '
  !  else
  !    write(2,80)'* FCB (FECAL COLIFORM BACTERIA) IS NOT MODELED    '
  !  endif
  !
  !  if( IWQSRP == 1 )then
  !    write(2,80)'* TAM IS USED FOR SORPTION OF PO4T/SA: MODEL TAM  '
  !  elseif( IWQSRP == 2 )then
  !    write(2,80)'* TSS IS USED FOR SORPTION OF PO4T/SA: MODEL TSS  '
  !    if( ISTRAN(6) /= 1) CALL STOPP('ERROR! INCOMPATIBLE ISTRAN(6)/IWQSRP')
  !  else
  !    write(2,80)'* NO SORPTION OF PO4T/SA: MAY MODEL TSS & NO TAM  '
  !  endif
  !
  !  if( IWQSTOX == 1 )then
  !    write(2,80)'* SALINITY TOXICITY IS APPLIED TO CYANOBACTERIA   '
  !  else
  !    write(2,80)'* NO SALINITY TOXICITY: SALTWATER CYANOBACTERIA   '
  !  endif
  !
  !  if( IWQKA(1) == 0 )then
  !    write(2,80)'* USER-SPECIFIED CONSTANT REAERATION SET TO WQKRO '
  !    write(2,80)'*   REAERATION DUE TO WIND SET TO ZERO            '
  !  endif
  !
  !  if( IWQKA(1) == 1 )then
  !    write(2,80)'* USER-SPECIFIED CONSTANT REAERATION SET TO WQKRO '
  !    write(2,80)'*   REAERATION DUE TO WIND ADDED TO WQKRO         '
  !  endif
  !  if( IWQKA(1) == 2 )then
  !    write(2,80)'* OCONNOR-DOBBINS REAERATION FORMULA IS USED      '
  !  endif
  !  if( IWQKA(1) == 3 )then
  !    write(2,80)'* OWENS & GIBBS (1964) REAERATION FORMULA IS USED '
  !  endif
  !  if( IWQKA(1) == 4 )then
  !    write(2,80)'* MODIFIED OWENS & GIBBS REAERATION IS USED       '
  !  endif
  !
  !  if( IWQVLIM == 0 )then
  !    write(2,80)'* MACROALGAE GROWTH IS NOT LIMITED BY VELOCITY    '
  !  endif
  !  if( IWQVLIM == 1 )then
  !    write(2,80)'* MACROALGAE VELOCITY LIMIT, MICHAELIS-MENTON EQU.'
  !  endif
  !  if( IWQVLIM == 2 )then
  !    write(2,80)'*MACROALGAE VEL. LIMIT, 5-PARAM LOGISTIC FUNCTION'
  !  endif
  !
  !  write(2,83)'* # OF ZONES FOR SPATIALLY VARYING parameterS  = ',IWQZ
  !  if( IWQZ > NWQZ) CALL STOPP('ERROR!! IWQZ SHOULD BE <= NWQZ')
  !
  !  if( IWQNC > 0 )then
  !    write(2,80)'* WRITE NEGATIVE CONC. INFORMATION TO NEG-CONC.LOG'
  !  else
  !    write(2,80)'* NO WRTING OF NEGATIVE CONCENTRATION INFORMATION '
  !  endif
  !
  !  if( IWQRST == 1 )then
  !    write(2,80)'* WRITE SPATIAL DISTRIBUTIONS TO IWQORST          '
  !  else
  !    write(2,80)'* NO WRITING TO IWQORST                           '
  !  endif
  !  write(2,999)
  !
  !  if( IWQICI == 1 )then
  !    write(2,80)'* SPATIALLY/TEMPORALLY-VARYING ICS FROM INWQICI   '
  !  elseif( IWQICI == 2 )then
  !    write(2,80)'* SPATIALLY/TEMPORALLY-VARYING ICS FROM INWQRST   '
  !  else
  !    write(2,80)'* SPATIALLY/TEMPORALLY CONSTANT INITIAL CONDITIONS'
  !  endif
  !
  !  if( IWQAGR == 1 )then
  !    write(2,80)'* SPATIALLY A/O TEMPORALLY-VARYING ALGAL KINETICS '
  !  else
  !    write(2,80)'* SPATIALLY/TEMPORALLY CONSTANT ALGAL KINETICS    '
  !  endif
  !
  !  if( IWQSTL == 1 )then
  !    write(2,80)'* SPATIALLY AND/OR TEMPORALLY-VARYING SETTLING VEL'
  !  else
  !    write(2,80)'* SPATIALLY/TEMPORALLY CONSTANT SETTLING VELOCITY '
  !  endif
  !
  !  if( IWQSUN >= 1 )then
  !    write(2,80)'* TEMPORALLY-VARYING IO & FD                      '
  !  else
  !    write(2,80)'* TEMPORALLY CONSTANT IO & FD                     '
  !  endif
  !
  !  if( IWQNPL == 1 )then
  !    write(2,80)'* SPATIALLY AND/OR TEMPORALLY-VARYING NPS INPUT   '
  !  else
  !    write(2,80)'* SPATIALLY/TEMPORALLY CONSTANT NPS INPUT         '
  !  endif
  !
  !  if( IWQKIN == 1 )then
  !    write(2,80)'* SPATIALLY VARYING KINETICS FROM KINETICS.INP    '
  !  else
  !    write(2,80)'* FILE KINETICS.INP NOT USED                      '
  !  endif
  !  write(2,999)
  !
  !  write(2,84) '* TIME-SERIES OUTPUT FROM ', WQTSB, ' DAY ', &
  !    '                       TO ', WQTSE, ' DAY ', &
  !    '                    EVERY ', WQTSDT, ' HOUR', &
  !    '                       AT ', IWQTS,  ' LOCATIONS', &
  !    ' BIN FILE SWITCH ISWQAVG  = ', ISWQAVG,' (0 = OFF)  ', &
  !    ' BIN FILE SWITCH ISWQMIN  = ', ISWQMIN,' (0 = OFF)  ', &
  !    ' BIN FILE SWITCH ISWQMAX  = ', ISWQMAX,' (0 = OFF)  ', &
  !    ' BIN FILE SWITCH ISCOMP   = ', ISCOMP, ' (0 = OFF)  '
  !  write(2,999)
  !
  !  if( IWQTS >= 1 )then
  !    write(2,80)': ICWQTS(I) = 1, TIME-SERIES OUTPUT FOR VARIABLE I  '
  !    write(2,80)': ICWQTS(I)\ = 1, NO TIME-SERIES OUTPUT FOR VAR. I  '
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
  !  if( ISTRWQ(22)  > 0 )then
  !    !AA! Added read of CO2 half-saturation consts
  !    write(2,80)'* HALF-SAT. CONSTANT (G/M^3) FOR NUTRIENT UPTAKE '
  !    write(2,81)' : (KHNC, KHPC, KHCO2C)       = ', WQKHNC,WQKHPC,WQKHCO2C
  !    write(2,81)' : (KHND, KHPD, KHS, KHCO2D)  = ', WQKHND,WQKHPD,WQKHS,WQKHCO2D
  !    write(2,81)' : (KHNG, KHPG, KHCO2G)       = ', WQKHND,WQKHPG,WQKHCO2G
  !    write(2,81)' : (KHNM, KHPM, KHCO2M)       = ', WQKHNM,WQKHPM,WQKHCO2M
  !  else
  !    write(2,80)'* HALF-SAT. CONSTANT (G/M^3) FOR NUTRIENT UPTAKE '
  !    write(2,81)' : (KHNC, KHPC)       = ', WQKHNC,WQKHPC
  !    write(2,81)' : (KHND, KHPD, KHS)  = ', WQKHND,WQKHPD,WQKHS
  !    write(2,81)' : (KHNG, KHPG)       = ', WQKHND,WQKHPG
  !    write(2,81)' : (KHNM, KHPM)       = ', WQKHNM,WQKHPM
  !  endif
  !  write(2,82)'* SAL. WHERE MICROSYSTIS GROWTH IS HALVED  = ', WQSTOX
  !
  !  write(2,80)'* LIGHT EXTINC. COEFF. DUE TO TSS, CHL & POM      '
  !  write(2,81)' : KETSS (/M PER G/M^3)  = ', WQKETSS
  !  write(2,81)' : KECHL (/M PER MG/M^3) = ', WQKECHL
  !  write(2,81)' : KECHLE (CHL EXPONENT) = ', WQKECHLE
  !  if( WQKECHL  < 0.0 )then
  !    write(2,80) '* Use RILEY (1956) EQUATION FOR WQKECHL          '
  !    write(2,80) ' : KECHL = 0.054*CHL**0.667 + 0.0088*CHL         '
  !  endif
  !  write(2,81)' : KEPOM (/M PER G/M^3)  = ', WQKEPOC
  !  write(2,80)'* CARBON-TO-CHL RATIO (G C PER MG CHL)            '
  !  write(2,81)' : (CCHLC, CCHLD, CCHLG) = ', WQCHLC,WQCHLD,WQCHLG
  !  write(2,80)'* DEPTH (M) OF MAXIMUM ALGAL GROWTH               '
  !  write(2,81)' : (DOPTC, DOPTD, DOPTG) = ', WQDOPC,WQDOPD,WQDOPG
  !
  !  write(2,82)'*INITIAL IO (LY/D) AT WATER SURFACE       = ',WQI0 &
  !    ,' MINIMUM OPTIMUM SOLAR RADIATION (LY/D)   = ',WQISMIN &
  !    ,' FRACTIONAL DAYLENGTH                     = ',WQFD &
  !    ,' WEIGHTING FACTOR FOR RAD. AT CURRENT DAY = ',WQCIA &
  !    ,' WEIGHTING FACTOR FOR RAD. AT (-1) DAY    = ',WQCIB &
  !    ,' WEIGHTING FACTOR FOR RAD. AT (-2) DAYS   = ',WQCIC &
  !    ,' FRACTION OF SOLAR RADIATION THAT IS PAR  = ',PARADJ
  !
  !  write(2,80)'* LOWER OPTIMUM TEMP FOR ALGAL GROWTH (DEGC)     '
  !  write(2,81)' : (TMC1, TMD1, TMG1   ) = ', WQTMC1,WQTMD1,WQTMG1
  !  write(2,80)'* UPPER OPTIMUM TEMP FOR ALGAL GROWTH (DEGC)     '
  !  write(2,81)' : (TMC2, TMD2, TMG2   ) = ', WQTMC2,WQTMD2,WQTMG2
  !
  !  write(2,80)'* REFERENCE TEMPERATURE FOR ALGAL METABOLISM (OC) '
  !  write(2,81)' : (TRC, TRD, TRG)       = ', WQTRC,WQTRD,WQTRG
  !  write(2,80)'* TEMPERATURE EFFECT FOR ALGAL METABOLISM         '
  !  write(2,81)' : (KTBC, KTBD, KTBG)    = ', WQKTBC,WQKTBD,WQKTBG
  !
  !  write(2,80)'* CARBON DISTRIBUTION COEFF FOR ALGAL PREDATION   '
  !  write(2,81)' : (FCRP, FCLP, FCDP)    = ', WQFCRP,WQFCLP,WQFCDP
  !  write(2,80)'* CARBON DISTRIBUTION COEFF FOR ALGAL METABOLISM  '
  !  write(2,81)' : (FCDC, FCDD, FCDG)    = ', WQFCDC,WQFCDD,WQFCDG
  !  write(2,80)'* HALF-SAT. CONSTANT (GO/M*3) FOR ALGAL DOC EXCRET'
  !  write(2,81)' : (KHRC, KHRD, KHRG)    = ', WQKHRC,WQKHRD,WQKHRG
  !
  !  write(2,80)'* MINIMUM DISSOLUTION RATE (/DAY) OF ORGANIC C    '
  !  write(2,81)' : (KRC, KLC, KDC)       = ', WQKRC,WQKLC,WQKDC(1)
  !  write(2,80)'* CONSTANT RELATING DISSOLUTION RATE TO ALGAE     '
  !  write(2,81)' : (KRCALG,KLCALG,KDCALG)= ', WQKRCALG,WQKLCALG,WQKDCALG
  !
  !  write(2,80)'* REFERENCE TEMP FOR HYDROLYSIS/MINERALIZATION(OC)'
  !  write(2,81)' : (TRHDR, TRMNL)        = ', WQTRHDR,WQTRMNL
  !  write(2,80)'* TEMPERATURE EFFECT ON HYDROLYSIS/MINERALIZATION '
  !  write(2,81)' : (KTHDR, KTMNL)        = ', WQKTHDR,WQKTMNL
  !  write(2,80)'* HALF-SAT. CONSTANT FOR OXIC RESP/DENITRIFICATION'
  !  write(2,81)' : (KHORDO, KHDNN)       = ', WQKHORDO,WQKHDNN
  !  write(2,80)'* RATION OF DENITRIFICATION TO OXIC DOC RESP      '
  !  write(2,81)' : (AANOX)               = ', WQAANOX
  !
  !  write(2,80)'* PHOSPHORUS DISTRIBUTION COEF FOR ALGAL PREDATION'
  !  write(2,81)' : (FPRP,FPLP,FPDP,FPIP) = ', WQFPRP,WQFPLP,WQFPDP,WQFPIP
  !  write(2,80)'* PHOSPHORUS DIST COEF OF RPOP FOR ALGAL METABOLIS'
  !  write(2,81)' : (FPRC, FPRD, FPRG)    = ', WQFPRC,WQFPRD,WQFPRG
  !  write(2,80)'* PHOSPHORUS DIST COEF OF LPOP FOR ALGAL METABOLIS'
  !  write(2,81)' : (FPLC, FPLD, FPLG)    = ', WQFPLC,WQFPLD,WQFPLG
  !
  !
  !  if( IWQSRP /= 1 .and. IWQSRP /= 2 )then
  !    WQKPO4P = 0.0
  !    write(2,80)': NO SORPTION OF PO4T/SA, SO KPO4P IS FORCED TO 0 '
  !  endif
  !  write(2,80)'* PHOSPHORUS DIST COEF OF DOP FOR ALGAL METABOLISM'
  !  write(2,81)' : (FPDC, FPDD, FPDG)    = ', WQFPDC,WQFPDD,WQFPDG
  !  write(2,80)'* PHOSPHORUS DIST COEF OF NH4 FOR ALGAL METABOLISM'
  !  write(2,81)' : (FPIC, FPID, FPIG)    = ', WQFPIC,WQFPID,WQFPIG
  !  write(2,82)'* PARITITION COEFF FOR SORBED/DISSOLVED PO4  = ', WQKPO4P
  !
  !  write(2,80)'* MINIMUM HYDROLYSIS RATE (/DAY) OF ORGANIC P     '
  !  write(2,81)' : (KRP, KLP, KDP)       = ', WQKRP,WQKLP,WQKDP
  !  write(2,80)'* CONSTANT RELATING HYDROLYSIS RATE TO ALGAE      '
  !  write(2,81)' : (KRPALG,KLPALG,KDPALG)= ', WQKRPALG,WQKLPALG,WQKDPALG
  !  write(2,80)'* CONSTANT USED IN DETERMINING P-TO-C RATIO       '
  !  write(2,81)' : (CPPRM1,CPPRM2,CPPRM3)= ', WQCP1PRM,WQCP2PRM,WQCP3PRM
  !
  !  write(2,80)'* NITROGEN DISTRIBUTION COEFF FOR ALGAL PREDATION '
  !  write(2,81)' : (FNRP,FNLP,FNDP,FNIP) = ', WQFNRP,WQFNLP,WQFNDP,WQFNIP
  !  write(2,80)'* NITROGEN DIST COEF OF RPON FOR ALGAL METABOLISM '
  !  write(2,81)' : (FNRC, FNRD, FNRG)    = ', WQFNRC,WQFNRD,WQFNRG
  !  write(2,80)'* NITROGEN DIST COEF OF LPON FOR ALGAL METABOLISM '
  !  write(2,81)' : (FNLC, FNLD, FNLG)    = ', WQFNLC,WQFNLD,WQFNLG
  !
  !  write(2,80)'* NITROGEN DIST COEF OF DON FOR ALGAL METABOLISM  '
  !  write(2,81)' : (FNDC, FNDD, FNDG)    = ', WQFNDC,WQFNDD,WQFNDG
  !  write(2,80)'* NITROGEN DIST COEF OF NH4 FOR ALGAL METABOLISM  '
  !  write(2,81)' : (FNIC, FNID, FNIG)    = ', WQFNIC,WQFNID,WQFNIG
  !  write(2,80)'* NITROGEN-TO-CARBON RATIO IN ALGAE               '
  !  write(2,81)' : (ANCC, ANCD, ANCG)    = ', WQANCC,WQANCD,WQANCG
  !
  !  write(2,82)'* MASS NO3 REDUCED PER DOC OXIDIZED (GN/GC)= ',WQANDC &
  !    ,'  MAXIMUM NITRIFICATION RATE (/DAY)        = ',WQNITM &
  !    ,'  REFERENCE TEMP FOR NITRIFICATION (DEGC)  = ',WQTNIT
  !  write(2,80)'* NITRIFICATION HALF-SAT CONSTANT FOR DO & NH4    '
  !  write(2,81)' : (KHNITDO, KHNITN)     = ', WQKHNDO,WQKHNN
  !  write(2,80)'* SUB & SUPER-OPTIMUM TEMP EFFECT ON NITRIFICATION'
  !  write(2,81)' : (KNIT1, KNIT2)        = ', WQKN1,WQKN2
  !
  !  write(2,80)'* MINIMUM HYDROLYSIS RATE (/DAY) OF ORGANIC N     '
  !  write(2,81)' : (KRN, KLN, KDN)       = ', WQKRN,WQKLN,WQKDN
  !  write(2,80)'* CONSTANT RELATING HYDROLYSIS RATE TO ALGAE      '
  !  write(2,81)' : (KRNALG,KLNALG,KDNALG)= ', WQKRNALG,WQKLNALG,WQKDNALG
  !
  !  if( IWQSRP /= 1 .and. IWQSRP /= 2 )then
  !    WQKSAP = 0.0
  !    write(2,80)': NO SORPTION OF PO4T/SA, SO KSAP IS FORCED TO 0  '
  !  endif
  !  write(2,80)'* SILICA DISTRIBUTION COEFF FOR DIATOM PREDATION  '
  !  write(2,81)' : (FSPP, FSIP)          = ', WQFSPP,WQFSIP
  !  write(2,80)'* SILICA DISTRIBUTION COEFF FOR DIATOM METABOLISM '
  !  write(2,81)' : (FSPD, FSID)          = ', WQFSPD,WQFSID
  !  write(2,82)'*SILICA-TO-CARBON RATIO IN DIATOMS        = ',WQASCD &
  !    ,'*PARITITION COEFF FOR SORBED/DISSOLVED SA = ',WQKSAP &
  !    ,'*DISSOLUTION RATE (/D) OF PSI             = ',WQKSU &
  !    ,' REFERENCE TEMP FOR PSI DISSOLUTION (OC)  = ',WQTRSUA &
  !    ,' TEMPERATURE EFFECT ON PSI DISSOLUTION    = ',WQKTSUA
  !
  !  write(2,82)'* DO-TO-CARBON RATIO IN RESPIRATION        = ',WQAOCR &
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
  !  write(2,82) &
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
  !  if( NWQOBS > NBBSM) CALL STOPP('ERROR!! NWQOBS SHOULD <= NBBSM')
  !  if( NWQOBW > NBBWM) CALL STOPP('ERROR!! NWQOBW SHOULD <= NBBWM')
  !  if( NWQOBE > NBBEM) CALL STOPP('ERROR!! NWQOBE SHOULD <= NBBEM')
  !  if( NWQOBN > NBBNM) CALL STOPP('ERROR!! NWQOBN SHOULD <= NBBNM')
  !
  !  write(2,999)
  !  write(2,80)'* CONSTANT OBC AT (ICBX(M),JCBX(M)) IF IWQOBX(M) = 0'
  !  write(2,80)': READ TIME-SERIES OBCS IWQOBX TIMES IF IWQOBX > 0'
  !
  !  !  *** C44
  !  ! SPATIALLY/TEMPORALLY CONSTANT INITIAL CONDITIONS: WQCHLX = 1/WQCHLX
  !  ! READ DATA POINTS & DO INTERNAL INTERPOLATION?
  !
  !  !IF(IWQICI == 0 )then    PMC - ALLOW INITIALIZATION, EVEN IF NOT USING CONSTANT IC.  THESE WILL BE OVERWRITTEN LATER, IF NEEDED
  !  write(2,999)
  !  !WRITE(2,90) TITLE_(1)
  !  !WRITE(2,21)' : (BC, BD, BG)         = ', (WQV(1,1,NW),NW = 1,3)
  !  !WRITE(2,21)' : (RPOC, LPOC, DOC)    = ', (WQV(1,1,NW),NW = 4,6)
  !  !WRITE(2,21)' : (RPOP,LPOP,DOP,PO4T) = ', (WQV(1,1,NW),NW = 7,10)
  !  !WRITE(2,21)' : (RPON, LPON, DON)    = ', (WQV(1,1,NW),NW = 11,13)
  !  !WRITE(2,21)' : (NH4, NO3)           = ', (WQV(1,1,NW),NW = 14,15)
  !  !WRITE(2,21)' : (SU, SA, COD, DO)    = ', (WQV(1,1,NW),NW = 16,19)
  !  !WRITE(2,981)' : (TAM, FCB, CO2, MAC)= ', (WQV(1,1,NW),NW = 20,NWQV),WQV(1,1,IDNOTRVA)    !VB! CHANGED NWQV+1 TO NWQV   !VB added CO2
  !
  !  do NAL = 1,NALAGE
  !    if( ALGAES(NAL).ISMOBILE == 0 )then
  !    endif
  !  enddo
  !  
  !  do NAL = 1,NALAGE
  !    if( ALGAES(NAL).ISMOBILE == 0 )then
  !    write(2,9003)
  !9003 FORMAT(/,' MACALGMP.INP - MACROALGAE MAP FILE',/, &
  !      ' PSHADE = SHADE FACTOR FOR TREE CANOPY (1.0 = NO CANOPY)',/, &
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
  !    write(2,9004) II, JJ, LL, PSHADE(LL), WQKMV(LL), WQKMVMIN(LL), WQKBP(LL), WQKMVA(LL), WQKMVB(LL), WQKMVC(LL), WQKMVD(LL), WQKMVE(LL)
  !9004 FORMAT(' ',I4,' ',I4,' ',I5, 9F7.3)
  !    endif
  !  enddo
  !
  !  write(2,80)'* ALGAL GROWTH RATE (/DAY)                        '
  !  !WRITE(2,21)' : (PMC, PMD, PMG)       = ', WQPMC(1),WQPMD(1), WQPMG(1)
  !  write(2,80)'* ALGAL BASAL METABOLISM RATE (/DAY)              '
  !  !WRITE(2,21)' : (BMRC, BMRD, BMRG)    = ', WQBMRC(1),WQBMRD(1), WQBMRG(1)
  !  write(2,80)'* ALGAL PREDATION RATE (/DAY)                     '
  !  !WRITE(2,21)' : (PRRC, PRRD, PRRG)    = ', WQPRRC(1),WQPRRD(1), WQPRRG(1)
  !
  !  write(2,82) '* BASE LIGHT EXTINCTION COEFFICIENT (/M)   = ',WQKEB(1)
  !
  !  if( IWQSTL /= 1 )then
  !    write(2,80)'* ALGAL SETTLING RATE (M/DAY)                     '
  !    !WRITE(2,21)' : (WSC, WSD, WSG)       = ', WQWSC(1),WQWSD(1),WQWSG(1)
  !    write(2,80)'* POM SETTLING RATE (M/DAY)                       '
  !    !WRITE(2,21)' : (WSRP, WSLP)          = ', WQWSRP(1),WQWSLP(1)
  !    write(2,80)'* SETTLING RATE OF PARTICULATE METAL (M/DAY)      '
  !    !WRITE(2,21)' : (WSS)                 = ', WQWSS(1)
  !  endif
  !
  !  ! *** C47 SPATIALLY/TEMPORALLY CONSTANT BENTHIC FLUXES
  !  if( IWQBEN == 0 )then
  !    write(2,999)
  !    write(2,90) TITLE_(1)
  !    !WRITE(2,21)' : (PO4D, NH4, NO3)     = ',WQBFPO4D(1),WQBFNH4(1),WQBFNO3(1)
  !    !WRITE(2,21)' : (SAD, COD, DO)       = ',WQBFSAD(1),WQBFCOD(1),WQBFO2(1)
  !  endif
  !
  !  !WRITE(2,23)'* NUMBER OF CELLS FOR POINT SOURCE INPUT  = ',IWQPS
  !  !WRITE(2,23)'* NUMBER WITH VARIABLE POINT SOURCE INPUT = ',NPSTMSR
  !  if( IWQPS > NWQPS) CALL STOPP('ERROR!! IWQPS SHOULD BE <= NWQPS')
  !
  !  if( IWQNPL /= 1 )then
  !    write(2,999)
  !    write(2,90) TITLE_(1)
  !    !WRITE(2,21)' : (DSQ, CHC, CHD, CHG)  = ',XDSQ,(XDSL(NW),NW = 1,3)
  !    !WRITE(2,21)' : (ROC, LOC, DOC)       = ',(XDSL(NW),NW = 4,6)
  !    !WRITE(2,21)' : (ROP, LOP, DOP, P4D)  = ',(XDSL(NW),NW = 7,10)
  !    !WRITE(2,21)' : (RON, LON, DON)       = ',(XDSL(NW),NW = 11,13)
  !    !WRITE(2,21)' : (NHX, NOX)            = ',(XDSL(NW),NW = 14,15)
  !    !WRITE(2,21)' : (SUU, SAA, COD, DOX)  = ',(XDSL(NW),NW = 16,19)
  !    !WRITE(2,981)' : (TAM, FCB, CO2)  = ',(XDSL(NW),NW = 20,NWQV)    !V!B ADDED co2
  !
  !    ! *** CONVERT FROM M3/S TO M3/DAY
  !    WQTT = XDSQ*CONV2
  !
  !    ! *** C50 WET DEPOSTION (MULTIPLIED BY RAINFALL VOLUME IN WQWET)
  !    !WRITE(2, 21)' : (CHC, CHD, CHG)       = ',(WQATM(NW),NW = 1,3)
  !    !WRITE(2, 21)' : (ROC, LOC, DOC)       = ',(WQATM(NW),NW = 4,6)
  !    !WRITE(2, 21)' : (ROP, LOP, DOP, P4D)  = ',(WQATM(NW),NW = 7,10)
  !    !WRITE(2, 21)' : (RON, LON, DON)       = ',(WQATM(NW),NW = 11,13)
  !    !WRITE(2, 21)' : (NHX, NOX)            = ',(WQATM(NW),NW = 14,15)
  !    !WRITE(2, 21)' : (SUU, SAA, COD, DOX)  = ',(WQATM(NW),NW = 16,19)
  !    !WRITE(2,981)' : (TAM, FCB, CO2)       = ',(WQATM(NW),NW = 20,NWQV)  !!VB ADDED CO2
  !endif
  !
  !  !READ(1,295) RSTOFN
  !  !WRITE(2,85)'* OUTPUT FILE FOR RESTART WRITING         = ', RSTOFN
  !  if( IWQRST == 1 )then
  !  else
  !    !IF( RSTOFN(1:4) /= 'NONE' .and. RSTOFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQORST/RSTOFN')
  !  endif
  !
  !  !READ(1,295) ICIFN
  !  !WRITE(2,85)'* FILE FOR INITIAL CONDITIONS             = ', ICIFN
  !  if( IWQICI == 1 )then
  !  elseif( IWQICI == 2 )then
  !  else
  !    !IF( ICIFN(1:4) /= 'NONE' .and. ICIFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQICI/ICIFN')
  !  endif
  !
  !  !READ(1,295) AGRFN
  !  !WRITE(2,85)'* FILE FOR ALGAL GROWTH, RESP., PREDATAT. = ', AGRFN
  !  if( IWQAGR == 1 )then
  !  else
  !    !IF( AGRFN(1:4) /= 'NONE' .and. AGRFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQAGR/AGRFN')
  !  endif
  !
  !  !READ(1,295) STLFN
  ! ! WRITE(2,85)'* FILE FOR SETTLING RATES OF ALGAE, PART. = ', STLFN
  !  if( IWQSTL == 1 )then
  !  else
  !    !IF( STLFN(1:4) /= 'NONE' .and. STLFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQSTL/STLFN')
  !  endif
  !
  !  !READ(1,295) SUNFN
  !  !WRITE(2,85)'* FILE FOR IO, FD, TE, KT                 = ', SUNFN
  !  if( IWQSUN == 1 )then
  !  else
  !    !   if( SUNFN(1:4) /= 'NONE' .and. SUNFN(1:4) /= 'none')
  !    !&     call STOPP('ERROR!! INVALID IWQSUN/SUNFN')
  !  endif
  !
  !  !READ(1,295) BENFN
  !  !WRITE(2,85)'* FILE FOR BENTHIC FLUX                   = ', BENFN
  !  if( IWQBEN == 2 )then
  !  else
  !    !IF( BENFN(1:4) /= 'NONE' .and. BENFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQBEN/BENFN')
  !  endif
  !
  !  !READ(1,295) PSLFN
  !  !WRITE(2,85)'* FILE FOR POINT SOURCE INPUT             = ', PSLFN
  !
  !  !READ(1,295) NPLFN
  !  !WRITE(2,85)'* FILE FOR NPS INPUT INCLUDING ATM. INPUT = ', NPLFN
  !  if( IWQNPL == 1 )then
  !  else
  !    !IF( NPLFN(1:4) /= 'NONE' .and. NPLFN(1:4) /= 'none') CALL STOPP('ERROR!! INVALID IWQNPL/NPLFN')
  !  endif
  !
  !  !READ(1,295) NCOFN
  !  !WRITE(2,85)'* DIAGNOSTIC FILE FOR NEGATIVE CONCENTRAT = ', NCOFN
  !  close(1)
  !
  !  if( IWQZ  > 1 .and. IWQKIN  > 0 )then
  !    !WRITE(*,'(A)')' WQ: KINETICS.INP'
  !    do I = 1,IWQZ
  !      write(2,9112) IZ, IWQKA(IZ), WQKRO(IZ), WQKTR(IZ), REAC(IZ),WQKDC(IZ),WQKDCALM(IZ),WQKHRM(IZ),WQDOPM(IZ),WQKCD(IZ),WQKHCOD(IZ)
  !    enddo
  !  endif
  !9111 FORMAT(/,'ZONE IWQKA   KRO   KTR  REAC   KDC KDCALGM  KHRM DOPTM   KCD KHCOD')
  !9112 FORMAT(I4, I6, 4F6.3, F8.3, 4F6.3)
  !
  !  if( IWQBEN == 2 )then
  !    write(*,'(A)')' WQ: WQBENMAP.INP'
  !    write(2,999)
  !    !WRITE(2,34) L, I, J, XBENMUD(L), IBENMAP(L,1), IBENMAP(L,2)
  !  endif


  END SUBROUTINE WQ_QC

  real(RKD) FUNCTION DO_SAT(L)
  
    implicit none
  
    integer, intent(IN) :: L
    real(RKD) :: ELEV, TVAL, RLNSAT1, RLNSAT2, RLNSAT3

    ! Elevation adjustment factor for D.O. saturation
    RLNSAT3 = 1.
    if( IDOSELE > 0 )then
      ELEV = BELV(L) + HP(L) + DOELEV
      if( IDOSELE == 1 )then
        ! *** Chapra (1997) pg. 3
        ELEV = 1.0e-3 * ELEV                      ! Convert elevation to km
        RLNSAT3 = (((-1.60747e-4 * ELEV + 6.10834e-3)*ELEV - 0.11988)*ELEV + 1.)
      elseif( IDOSELE == 2 )then
        ! *** Zison et al. (1978)
        RLNSAT3 = (1.0 - 0.1148e-3 * ELEV)
      endif
    endif

    ! THE D.O. SATURATION CALCULATION
    if( IDOSFRM == 2 )then
      ! *** Genet et al. (1974)
      RLNSAT1 = ((+0.0054258 * TWQ(L) -0.38217)*TWQ(L) + 14.5532)
      DO_SAT = RLNSAT1 * RLNSAT3
    elseif( IDOSFRM == 1 )then
        ! *** Chapra (1997) pg. 3
      TVAL = 1./(TWQ(L) + 273.15)
      RLNSAT1 = ((((-8.621949e11 * TVAL + 1.2438e10)*TVAL -6.642308e7)*TVAL + 1.575701e5)*TVAL -139.34411)
      RLNSAT2 = -SWQ(L)*((-2140.7 * TVAL + 10.754)*TVAL + 1.7674e-2)
      DO_SAT = RLNSAT3 * EXP(RLNSAT1 + RLNSAT2)
    else
      ! *** DO Saturation, Modified by SCJ, see Garcia and Gordon, Limnology and Oceanography 37(6), 1992, Eqn. 8 and Table 1
      TVAL = LOG((298.15 - TWQ(L))/(273.15 + TWQ(L)))
      RLNSAT1 = (((((1.41575*TVAL + 1.01567)*TVAL + 4.93845)*TVAL + 4.11890)*TVAL + 3.20684)*TVAL + 5.80818)
      RLNSAT2 = -SWQ(L)*( 1.32412e-7*SWQ(L) + (((5.54491e-3*TVAL + 7.93334E-3)*TVAL + 7.25958e-3)*TVAL +7.01211e-3) )
      DO_SAT = 32.0e-3 * RLNSAT3 * EXP(RLNSAT1 + RLNSAT2)         ! *** 32E-3 approximately converts micromol/L to mg/L or g/m^3
    endif
  END FUNCTION DO_SAT

END MODULE WATERQUALITY
