! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CALQVS

  ! *** SUBROUTINE CALQVS UPDATES TIME VARIABLE VOLUME SOURCES

  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !--------------------------------------------------------------------------------------!
  !    2018-03       PAUL M. CRAIG     IMPLEMENTED WHOLE CHANNEL ELEVATION RATING CURVE (NQCTYP = 2)
  !    2018-02       N T LAM           IMPLEMENTED TIME VARIABLE GATE STRUCTURES
  !    2015-11       D H CHUNG  &      IMPLEMENTED NEW HYDRAULIC STRUCTURES BASED ON DILL
  !                  PAUL M. CRAIG
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2015-01       PAUL M. CRAIG     ADDED FULLY COUPLED ICE SUB-MODEL
  !                  DANG H CHUNG
  !    2014-09       PAUL M. CRAIG     IMPLEMENTED NEW LOW CHORD BC TYPE
  !    2014-08       D H CHUNG         SET EXPLICIT PRECISIONS of INTEGER & REAL

  use GLOBAL
  use Allocate_Initialize
  
  use Variables_WQ
  use HEAT_MODULE, only:ICECOMP
  use HYDSTRUCMOD
  use DIFFUSER_MODULE, only:JPEFDC
  use FIELDS

  use Variables_MPI
  use Variables_MPI_Mapping

  implicit none

  integer :: L, K, LL, LDR, LUR, NS, NCTMP, M1, M2, ITYP, NCTL, NCTLT, I, NC, IG
  integer :: IU, JU, ID, JD, LU, LD, ND, LF, IPMC, MU1, MU2, MD1, MD2
  integer :: NTMP, NWR, KU, KD, NJP, LJP, KTMP, ITMPD, NTT, LEVELFLAG, MWET
  integer :: IBLOCK
  integer, save, allocatable, dimension(:)   :: FLOWTYPE,  NUNIQQSERQ,  NQWET
  integer, save, allocatable, dimension(:,:) :: UNIQQSERQ
  integer, save :: NLIMITICE = 0
  
  real :: QWRABS, QSERTCELL, CTIM, TDIFF, WTM1, WTM2, HUP, HDW, DELH, HDRY9, TF, HPICE
  real :: TMPVAL, QCTLMAX, HTMPD, TDIFFU, WTM1U, WTM2U, TDIFFD, WTM1D, WTM2D, CLEVAPTMP
  real :: TVW, TVA, DTV, DTVL, LAMBDA, TM, VPTG                                           !< RYAN-HARLEMAN
  real :: RPORTS, QVJPTMP, RAVAIL, QSUMIET, QEAVAIL, DIFQVOL, DTAGW, RIFTRL, QSUMTMP
  real :: QSUM1, QSUM2                                                                    !< Flow accumulators to handle specified flows below active bottom layer (SGZ)
  real :: ZHU, ZHD, HVAL, CVAL, RVAL, DENAIR, FWUP, QSW, QAA
  real :: HOPEN, WOPEN, SOPEN
  
  real(RKD) :: SPLIT, QLAYER1, QLAYER2
  
  real,save :: ICESTEP,DELTICE1,DELTICE2,ICEDAY

  real,save,allocatable,dimension(:,:) :: QLAYER
  real,save,allocatable,dimension(:)   :: QSFACTTOT
  real,save,allocatable,dimension(:)   :: QSFACT
  real(RKD),save,allocatable,dimension(:,:) :: QSUM3                                !< Flow accumulator for in/out flows for a dry cell
  real(RKD),save,allocatable,dimension(:)   :: QSUM4                                !< Flow accumulator for in/out flows for a dry cell
  real(RKD) :: TIMESEC2                                                             !< Current simulation time in seconds, adjusted for time stepping option

  logical,save :: LGROUPS

  if( .not. allocated(QLAYER) )then
    ! *** FIRST CALL  FOR TYP = 3/4
    call AllocateDSI( USCELL,    NQCTL,   0)
    call AllocateDSI( DSCELL,    NQCTL,   0)

    call AllocateDSI( QSFACT,    NQCTL, 0.0)
    call AllocateDSI( QSFACTTOT, NQCTL, 0.0)
    call AllocateDSI( QLAYER,    NQCTL, KCM, 0.0)
    
    if( NGRPID > 0 .and. NQSIJ > 0 )then
      call AllocateDSI( FLOWTYPE,  NGRPID,   0)
      call AllocateDSI( NUNIQQSERQ, NGRPID,   0)
      call AllocateDSI( NQWET,      NGRPID,   0)
      call AllocateDSI( QSUM4,     NGRPID, 0.0)
      call AllocateDSI( UNIQQSERQ,  NQSIJ, NGRPID,   0)
      call AllocateDSI( QSUM3,     2,     NGRPID, 0.0)

      do LL = 1,NQSIJ
        FLOWTYPE(BCPS(LL).GRPID) = BCPS(LL).NQSMUL
      enddo

      ! *** LIST OF UNIQUE QSER SERIES FOR EACH GROUP
      do LL = 1,NQSIJ
        do I = 1,NUNIQQSERQ(BCPS(LL).GRPID)
          if( BCPS(LL).NQSERQ == UNIQQSERQ(I,BCPS(LL).GRPID) ) exit
        enddo
        if( I > NUNIQQSERQ(BCPS(LL).GRPID) )then
          NUNIQQSERQ(BCPS(LL).GRPID) = NUNIQQSERQ(BCPS(LL).GRPID) + 1
          UNIQQSERQ(NUNIQQSERQ(BCPS(LL).GRPID),BCPS(LL).GRPID) = BCPS(LL).NQSERQ
        endif
      enddo
    endif

    ! *** Save the initialize location for each NQSIJ
    do LL = 1,NQSIJ
      LQSSAVED(LL) = BCPS(LL).L
    enddo
    
    LEVELFLAG = 0
    do NCTL = 1,NQCTL
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      USCELL(NCTL) = LU
      if( ISRESTI == 0 )then
        SAVESUB(1,NCTL) = SUBO(LU)
        SAVESVB(1,NCTL) = SVBO(LU)
      endif

      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      DSCELL(NCTL) = LC
      if( ID > 0 .and. JD > 0 )then
        LD = LIJ(ID,JD)
        DSCELL(NCTL) = LD
        if( ISRESTI == 0 )then
          SAVESUB(2,NCTL) = SUBO(LD)
          SAVESVB(2,NCTL) = SVBO(LD)
        endif
      else
        LD = LC
      endif

      ! *** RESET IF USING RESTART
      if( ISRESTI /= 0 )then
        HRUO(LU) = SAVESUB(1,NCTL)*DYU(LU)*DXIU(LU)
        HRVO(LU) = SAVESVB(1,NCTL)*DXV(LU)*DYIV(LU)
        HRUO(LD) = SAVESUB(2,NCTL)*DYU(LD)*DXIU(LD)
        HRVO(LD) = SAVESVB(2,NCTL)*DXV(LD)*DYIV(LD)
      endif

      if( HYD_STR(NCTL).NQCTYP == 3 .or. HYD_STR(NCTL).NQCTYP == 4 )then
        ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE DIFFERENCE (4) DEPENDANT FLOWS
        NCTLT = HYD_STR(NCTL).NQCTLQ  ! *** SERIES INDEX

        ! *** ADJUST THE FLOWS TO ZERO AT LOWER CHORD DEPTHS
        if( HYD_STR(NCTL).NQCTYP == 3 )then
          DELH = HCTLUM(NCTLT)*(HYD_STR(NCTL).BQCLCE-BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU)

          ! *** LOOKUP THE FLOWS
          M1 = 0
          M2 = 1
        70 M1 = M1+1
          M2 = M2+1
          if( M2 > MQCTL(NCTLT) )then
            write(6,*)' BAD QCTL BOUNDARY SPECIFICATION.  LOW CHORD OUT OF RANGE OF RATING CURVE.'
            write(6,6667)NCTL,NCTLT,IU,JU,ID,JD
            write(6,*)' LOW CHORD:  ',HYD_STR(NCTL).BQCLCE,DELH
            open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            write(mpi_efdc_out_unit,*)' BAD QCTL BOUNDARY SPECIFICATION.  LOW CHORD OUT OF RANGE OF RATING CURVE.'
            write(mpi_efdc_out_unit,6667)NCTL,NCTLT,IU,JU,ID,JD
            write(mpi_efdc_out_unit,*)' LOW CHORD:  ',HYD_STR(NCTL).BQCLCE,DELH
            close(mpi_efdc_out_unit)
            call STOPP('.')
          endif
          if( DELH >= HDIFCTL(M1,NCTLT) .and. DELH <= HDIFCTL(M2,NCTLT) )then
            ! *** FOUND LOWER CHORD ELEVATION, DETERMINE THE FLOW OFFSET
            TDIFF = HDIFCTL(M2,NCTLT)-HDIFCTL(M1,NCTLT)
            WTM1 = (HDIFCTL(M2,NCTLT)-DELH)/TDIFF
            WTM2 = (DELH-HDIFCTL(M1,NCTLT))/TDIFF
            do K = 1,KC
              QLAYER(NCTL,K) = WTM1*QCTL(M1,1,K,NCTLT)+WTM2*QCTL(M2,1,K,NCTLT)
            enddo

            CYCLE
          else
            GOTO 70
          endif

          ! *** INITIALIZE IF WSEL>LOWCHORD
          if( ISRESTI == 0 )then
            HUP = BELV(LU)+HP(LU)
            if( HUP <= HYD_STR(NCTL).BQCLCE )NLOWCHORD(NCTL) = -(HYD_STR(NCTL).NQCMINS+1)
          endif

        endif
        
      elseif( HYD_STR(NCTL).NQCTYP == 9 )then
        ! *** Navigation locks.  Initialize current gate positions
        if( HSCTL(NCTL).ITYPE == 4 )then  
          ID = HSCTL(NCTL).ID
          IBLOCK = 0
          if( HSCTL(NCTL).ITYPE == 4 ) IBLOCK = 1
          call GATE_OPENING_INTERP(QCTLSER(ID), TIMESEC, HOPEN, WOPEN, SOPEN, IBLOCK)
          
          HYD_STR(NCTL).CURHEI = HOPEN     ! *** Set the curent time for the height (m)
          HYD_STR(NCTL).CURWID = WOPEN     ! *** Set the curent time for the width (m)
          HYD_STR(NCTL).CURSIL = SOPEN     ! *** Set the curent time for the sill height (m)
          
          ! *** Set gates open/closed status based on initial WSEL
          ZHU = BELV(LU) + HP(LU)
          ZHD = BELV(LD) + HP(LD)
          if( ABS(ZHU - ZHD) > 0.01 )then
            HYD_STR(NCTL).CURHEI = HOPEN + 1.0     ! *** Initialize gates to properly open/close on first call to COMPUTE_HSFLOW
            HYD_STR(NCTL).ISTATE = 0               ! *** Fully closed
          else
            HYD_STR(NCTL).ISTATE = 3               ! *** Fully open
          endif
        endif   
        
        ! *** Set navigation lock filling/draining flows
        if( HYD_STR(NCTL).TRANSIT > 0.0 )then
          ! *** Flows
          if( HYD_STR(NCTL).NWRGRP > 0 )then
            do NWR = 1,NQWR
              if( HYD_STR(NCTL).NWRGRP == WITH_RET(NWR).GROUPID )then
                HYD_STR(NCTL).QNWR = WITH_RET(NWR).QWR
              endif
            enddo
            HYD_STR(NCTL).ZHUR = -999.                         ! *** Flag uninitialized
          endif
          
          ! *** Find upstream Lock index
          HYD_STR(NCTL).IUSBC = 0
          LUR = LIJ(HSCTL(NCTL).IREFUP, HSCTL(NCTL).JREFUP)    ! *** Upstream gate cell
          if( LUR > 0 )then
            ! *** Find the BC index for the upstream cell
            do M1 = 1,NQCTL
              if( LUR == LIJ(HYD_STR(M1).IQCTLU,HYD_STR(M1).JQCTLU) )then
                HYD_STR(NCTL).IUSBC = M1
                exit
              endif
            enddo
          endif

          ! *** Find downstream Lock index
          HYD_STR(NCTL).IDSBC = 0
          LDR = LIJ(HSCTL(NCTL).IREFDN, HSCTL(NCTL).JREFDN)    ! *** Downstream gate cell
          if( LDR > 0 )then
            ! *** Find the BC index for the upstream cell
            do M1 = 1,NQCTL
              if( LDR == LIJ(HYD_STR(M1).IQCTLD,HYD_STR(M1).JQCTLD) )then
                HYD_STR(NCTL).IDSBC = M1
                exit
              endif
            enddo
          endif
        endif
      endif

    enddo

    ! *** Lock filling/draining - Reset associated W/R flow settings
    do NCTL = 1,NQCTL
      if( HYD_STR(NCTL).NQCTYP == 9 )then
        if( HYD_STR(NCTL).TRANSIT > 0.0 )then
          if( HYD_STR(NCTL).NWRGRP > 0 )then
            do NWR = 1,NQWR
              if( HYD_STR(NCTL).NWRGRP == WITH_RET(NWR).GROUPID )then
                WITH_RET(NWR).QWR = 0.0       ! *** Zero the constant flow and use lockage to control
                WITH_RET(NWR).NQWRSERQ = 0    ! *** Ensure that no flow series has been assigned
              endif
            enddo
          endif
        endif
        
        ! *** Set the current most upstream reference cell and WSEL
        M1 = NCTL
        M2 = NCTL
        do while (M2 > 0)
          M2 = HYD_STR(M2).IUSBC     ! *** Next upstream cell
          if( M2 == 0 .or. HYD_STR(M1).ISTATE == 0 )then
            HYD_STR(NCTL).LUR  = LIJ(HYD_STR(M1).IQCTLU,HYD_STR(M1).JQCTLU)
            HYD_STR(NCTL).ZHUR = BELV(HYD_STR(NCTL).LUR) + HP(HYD_STR(NCTL).LUR)
            exit
          endif
          M1 = M2                     ! *** Previous upstream cell (always valid)
        enddo        
      endif
    enddo
        
    LGROUPS = .FALSE.
    do NCTL = 1,NQCTL
      if( HYD_STR(NCTL).NQCTYP == -2 )then
        LGROUPS = .TRUE.
        if( HYD_STR(NCTL).QCTLGRP < 1  )then
          write(6,'(" BAD GROUPID AT NCTL: ",I5,I10)') NCTL,HYD_STR(NCTL).QCTLGRP
          call STOPP('.')
        endif
      endif
    enddo

    ! *** ICE SETTINGS
    ICESTEP = 60.*15.  ! *** DURATION TO ACHEIVE "QUASI-STEADY STATE" CONDITIONS FOR ICE SURFACE TEMPERATURE
    DELTICE1  = 0.0
    DELTICE2  = 0.0
    ICEDAY = INT(TIMEDAY) + 1.
  endif      ! *** End of first time step initializations

  ! *** Adjust current time based on time stepping option
  if( ISDYNSTP == 0 )then
    if( ISTL == 2 )then
      TIMESEC2 = TIMESEC - DT/2.    ! *** 3TL when ISTL = 2 and 2TL without dynamic timestepping
      DELT = DT
    else
      TIMESEC2 = TIMESEC - DT       ! *** 3TL when ISTL = 3
      DELT = DT2
    endif
  else
    TIMESEC2 = TIMESEC              ! *** 2TL dynamic time stepping
    DELT = DTDYN
  endif 

  HDRY9 = 0.9*HDRY

  ! ***  INITIALIZE NULL (0) FLOW SERIES
  GWSERT(0) = 0.
  QWRSERT(0) = 0.
  QSERTCELL = 0.0
  do K = 1,KC
    QSERT(K,0) = 0.
    QCTLT(K,0,1) = 0.
    QCTLT(K,0,2) = 0.
    QCTLTO(K,0) = 0.
  enddo

  ! *** 2018-10-08, NTL: ADD TIME VARIABLE RATING CURVE
  ! IMPLEMENTATION: RATING CURVES ARE SUDDENLY CHANGED BETWEEN TWO TIME STEPS
  do NCTL = 1,NQCTL
    if( HYD_STR(NCTL).NQCTYP < 5 )then
      if( HSCTL(NCTL).ITYPE == 1 )then
        ! *** STRUCTURE IS CONTROLLED BY TIME-SERIES
        ID = HSCTL(NCTL).ID
        call CHANGE_RATING_CURVE_BY_TIME(QCTLSER(ID), TIMESEC, HYD_STR(NCTL).NQCTLQ)
      elseif( HSCTL(NCTL).ITYPE == 2 .or. HSCTL(NCTL).ITYPE == 3 )then
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES
        call CHANGE_RATING_CURVE_BY_RULES(NCTL, HYD_STR(NCTL).NQCTLQ)
      else
        ! *** STRUCTURE IS UNCONTROLLED, DO NOTHING
      endif
    endif
  enddo

  !$OMP PARALLEL DEFAULT(SHARED)

  if( NGWSER >= 1 )then
    !$OMP SINGLE
    NCTMP = 3 + NDYM + NSED + NSND + NTOX
    do NC = 1,NCTMP
      GWCSERT(0,NC) = 0.
    enddo
    !$OMP END SINGLE

    if( ISTRAN(5) > 0 )then
      !$OMP DO PRIVATE(ND,LF,LL,L,NC)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do NC = 1,NCTMP
          do L = LF,LL
            CONGW(L,NC) = 0.0
          enddo
        enddo
      enddo
      !$OMP END DO
    endif
  endif

  ! ***  INITIALIZE TOTAL FLOW SERIES
  !$OMP DO PRIVATE(ND,LF,LL,L)
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = min(LF+LDM-1,LA)
    do L = LF,LL
      QSUM1E(L) = QSUME(L)
      QSUME(L)  = 0.
    enddo
  enddo
  !$OMP END DO

  ! *** SELECTIVE ZEROING
  if( KC > 1 )then
    if( NGWSER > 0 .or. ISGWIT /= 0 )then
      !$OMP DO PRIVATE(ND,LF,LL,L,K)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          QSUM(L,KSZ(L)) = 0.
        enddo
      enddo
      !$OMP END DO
    endif

    ! *** ZERO EVAP/RAINFALL/ICE
    !$OMP DO PRIVATE(ND,LF,LL,L)
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do L = LF,LL
        QSUM(L,KC) = 0.
      enddo
    enddo
    !$OMP END DO

    ! *** ZERO ALL DEFINED BC'S
    !$OMP SINGLE
    do NS = 1,NBCS
      L = LBCS(NS)
      do K = 1,KC
        QSUM(L,K) = 0.
      enddo
    enddo
    !$OMP END SINGLE

  else
    ! *** SINGLE LAYER
    !$OMP DO PRIVATE(ND,LF,LL,L)
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do L = LF,LL
        QSUM(L,KC) = 0.
      enddo
    enddo
    !$OMP END DO
  endif
  !$OMP END PARALLEL

  ! *************************************************************************************
  ! *** VOLUME SOURCE/SINK INTERPOLATION
  do NS = 1,NQSER
    CTIM = TIMESEC2/TSFL(NS).TMULT
    
    M2 = MTSQLAST(NS)
    do while (CTIM > TSFL(NS).TIM(M2))
      M2 = M2+1
      if( M2 > TSFL(NS).NREC )then
        M2 = TSFL(NS).NREC
        exit
      endif
    enddo
    MTSQLAST(NS) = M2
    M1 = M2-1
    TDIFF = TSFL(NS).TIM(M2)-TSFL(NS).TIM(M1)
    WTM1 = (TSFL(NS).TIM(M2)-CTIM)/TDIFF
    WTM2 = (CTIM-TSFL(NS).TIM(M1))/TDIFF
    do K = 1,KC
      QSERT(K,NS) = WTM1*TSFL(NS).VAL(M1,K) + WTM2*TSFL(NS).VAL(M2,K)
    enddo
  enddo
  
  if( N == 1 )then
    open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
    do LL = 1,NQSIJ
      L = BCPS(LL).L
      ITYP = LCT(L)
      if( ITYP <= 0 .or. ITYP >= 8 )then
        write(6,6111) LL,BCPS(LL).I,BCPS(LL).J
        write(mpi_efdc_out_unit,6111) LL,BCPS(LL).I,BCPS(LL).J
      endif
    enddo
    close(mpi_efdc_out_unit)
  endif

  ! *** If HDRYMOVE > 0.0 then make sure there is sufficient depth to receive inflows
  
  if( HDRYMOVE > 0.0 )then
    do LL = 1,NQSIJ
      LQSMOVED(LL) = BCPS(LL).L          ! *** Save the previous location for each NQSIJ
    enddo
    do LL = 1,NQSIJ
      if( BCPS(LL).NQSMF > 0 ) CYCLE     ! *** Only allow total flow rates at this time
      L = LQSSAVED(LL)                   ! *** Reset cell to original
      BCPS(LL).L = L                     ! *** Reset cell to original
      LBCS(LL) = L                       ! *** Reset cell to original
      LU = L
      LD = L
      HDW = BELV(LD) 
      M1 = 0
      do while( HP(LU) < HDRYMOVE )
        ! *** Search the surrounding cells to determine if the inflows can migrate
        ! *** Find the lowest adjacent cell
        do NS = 1,9
          LF = LADJ(NS,LU)
          if( LF == LU ) CYCLE
          if( LF == 1 )  CYCLE
          if( BELV(LF) < HDW )then
            LD = LF
            HDW = BELV(LD)
          endif                    
        enddo
        if( LD == LU ) exit
        LU = LD
        M1 = M1 + 1
        if( HP(LU) >= HDRYMOVE ) exit
      enddo
      BCPS(LL).L  = LU      ! *** Update receiving cell in the QSER list
      LBCS(LL) = LU         ! *** Update receiving cell in the global BC list
      QSUM(LU,:) = 0.       ! *** Zero new BC's cell
    enddo
  endif
  
  do LL = 1,NQSIJ
    NS = BCPS(LL).NQSERQ
    L = BCPS(LL).L
    QSUM1 = 0.
    QSUM2 = 0.
    do K = 1,KC
      ! *** APPLY MULTIPLIERS HERE TO CORRECT MASS BALANCE PROBLEMS
      if( BCPS(LL).NQSMUL < 5 .or. ISDRY == 0 )then
        ! *** FIXED BC GROUP DISTRIBUTION BY QFACTOR
        QSS(K,LL)      = BCPS(LL).QSSE   *BCPS(LL).RQSMUL*DZC(L,K)
        QSERCELL(K,LL) = QSERT(K,NS)*BCPS(LL).RQSMUL*BCPS(LL).QFACTOR

      elseif( BCPS(LL).NQSMUL == 5 )then
        ! *** REDISTRIBUTION OF DRY CELL FLOWS
        QSS(K,LL)      = BCPS(LL).QSSE*DZC(L,K)
        QSERCELL(K,LL) = QSERT(K,NS)*BCPS(LL).QFACTOR

        if( HP(L) < HDRY )then
          QSUM3(1,BCPS(LL).GRPID) = QSUM3(1,BCPS(LL).GRPID) + QSS(K,LL)
          QSUM3(2,BCPS(LL).GRPID) = QSUM3(2,BCPS(LL).GRPID) + QSERCELL(K,LL)
          QSS(K,LL) = 0.
          QSERCELL(K,LL) = 0.
        endif

      elseif( BCPS(LL).NQSMUL == 6 )then
        QSS(K,LL) = 0.
        QSERCELL(K,LL) = 0.
        CYCLE
      endif
      QSUM(L,K) = QSUM(L,K) + QSS(K,LL) + QSERCELL(K,LL)

      if( LKSZ(L,K) )then
        ! *** LAYER IS INACTIVE
        QSUM1 = QSUM1 + QSS(K,LL)
        QSUM2 = QSUM2 + QSERCELL(K,LL)
        QSS(K,LL)      = 0.
        QSERCELL(K,LL) = 0.
        QSUM(L,K)      = 0.
      endif
    enddo
    
    ! *** ADD FLOWS BELOW BOTTOM ACTIVE LAYER BACK TO ACTIVE LAYERS
    if( QSUM1 /= 0. )then
      ! *** CONSTANT FLOW
      do K = KSZ(L),KC
        QSS(K,LL) = QSS(K,LL) + QSUM1*DZC(L,K)
        QSUM(L,K) = QSUM(L,K) + QSUM1*DZC(L,K)
      enddo
    endif

    if( QSUM2 /= 0. )then
      ! *** TIME VARIABLE
      do K = KSZ(L),KC
        QSERCELL(K,LL) = QSERCELL(K,LL) + QSUM2*DZC(L,K)
        QSUM(L,K)      = QSUM(L,K)      + QSUM2*DZC(L,K)
      enddo
    endif
  enddo

  ! *** REDISTRIBUTE FLOWS INTO OTHER WET CELLS IN THE SAME GROUP
  if( NQSIJ > 0 )then
    do IG = 1,NGRPID
      if( FLOWTYPE(IG) == 5 )then
        if( QSUM3(1,IG) /= 0. .or. QSUM3(2,IG) /= 0. )then
          ! *** FOUND SOME FLOW, SO SPLIT FLOWS BY CELL VOLUME
          QSUM1 = 0.
          QSUM2 = 0.
          do LL = 1,NQSIJ
            if( BCPS(LL).GRPID == IG )then
              ! *** FOUND A CELL IN THE GROUP.  GET VOLUME
              L = BCPS(LL).L
              if( HP(L) > HDRY .or. ISDRY == 0 )then
                QSUM1 = QSUM1 + HP(L)*DXYP(L)            ! *** TOTAL WET VOLUME
              else
                QSUM2 = QSUM2 + HP(L)*DXYP(L)            ! *** TOTAL DRY VOLUME
              endif
            endif
          enddo

          ! *** NOW REDISTRIBUTE IF WET CELLS AVAILABLE
          do LL = 1,NQSIJ
            if( BCPS(LL).GRPID == IG )then
              ! *** TEST IF THERE ARE ANY WET CELLS IN GROUP
              L = BCPS(LL).L
              if( QSUM1 > 0. )then
                ! *** WET CELLS ARE AVAILABLE.  DISTRIBUTE FLOWS INTO GROUP
                if( HP(L) > HDRY .or. ISDRY == 0 )then
                  SPLIT = HP(L)*DXYP(L)/QSUM1                 ! *** FLOW SPLIT
                  do K = KSZ(L),KC
                    QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                    QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                    QSS(K,LL)      = QSS(K,LL)      + QLAYER1
                    QSERCELL(K,LL) = QSERCELL(K,LL) + QLAYER2
                    QSUM(L,K)      = QSUM(L,K)      + QLAYER1 + QLAYER2
                  enddo
                endif
              elseif( QSUM2 > 0. )then
                ! *** NO WET CELLS AVAILABLE.  JUST USE THE ORIGINAL FLOWS
                SPLIT = HP(L)*DXYP(L)/QSUM2                 ! *** FLOW SPLIT
                do K = KSZ(L),KC
                  QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                  QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                  QSS(K,LL)      = QSS(K,LL)      + QLAYER1
                  QSERCELL(K,LL) = QSERCELL(K,LL) + QLAYER2
                  QSUM(L,K)      = QSUM(L,K)      + QLAYER1 + QLAYER2
                enddo
              endif
            endif
          enddo
        endif

      elseif( FLOWTYPE(IG) == 6 )then
        ! *** REDISTRIBUTE TOTAL FLOW INTO WET CELLS BY VOLUME

        ! *** GET VOLUMES
        QSUM1 = 0.
        QSUM2 = 0.
        QSUM3(:,IG) = 0.
        MWET = 0

        ! *** GET TOTAL CONSTANT FLOW (QSUM3(1,IG)) INTO THIS GROUP 
        do LL = 1,NQSIJ
          if( BCPS(LL).GRPID == IG )then
            ! *** FOUND A CELL IN THE GROUP.  GET VOLUME
            NS = BCPS(LL).NQSERQ
            L = BCPS(LL).L
            if( HP(L) >= HDRY .or. ISDRY == 0 )then
              QSUM1 = QSUM1 + HP(L)*DXYP(L)            ! *** TOTAL WET VOLUME
              MWET = MWET + 1
            else
              QSUM2 = QSUM2 + HP(L)*DXYP(L)            ! *** TOTAL DRY VOLUME
            endif
            QSUM3(1,IG) = QSUM3(1,IG) + BCPS(LL).QSSE       ! *** SUM CONSTANT FLOW
            do K = 1,KC
              QSS(K,LL)      = 0.
              QSERCELL(K,LL) = 0.
            enddo
          endif
        enddo

        ! *** GET TOTAL TIME SERIES FLOW (QSUM3(2,IG)) INTO THIS GROUP 
        do I = 1,NUNIQQSERQ(IG)
          NS = UNIQQSERQ(I,IG)
          do K = 1,KC
            QSUM3(2,IG) = QSUM3(2,IG) + QSERT(K,NS)  ! *** SUM TIME SERIES FLOWS
          enddo
        enddo

        ! *** NOW REDISTRIBUTE IF WET CELLS AVAILABLE
        IU = 0
        do LL = 1,NQSIJ
          if( BCPS(LL).GRPID == IG )then
            ! *** TEST IF THERE ARE ANY WET CELLS IN GROUP
            L = BCPS(LL).L

            ! *** TEST IF UPDATE IS NEEDED
            if( QSUM1 > 0. )then
              if( MWET == NQWET(IG) .and. ABS(QSUM1-QSUM4(IG))/QSUM1 < 0.10 )then
                ! *** Use PREVIOUS TIME STEPS SPLITS IF LITTLE HAS CHANGED SINCE THE LAST UPDATED SPLITS
                SPLIT = BCPS(LL).QFACTOR
                do K = KSZ(L),KC
                  QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                  QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                  QSS(K,LL)      = QLAYER1
                  QSERCELL(K,LL) = QLAYER2
                  QSUM(L,K)      = QSUM(L,K) + QLAYER1 + QLAYER2
                enddo
              else
                ! *** WET CELLS ARE AVAILABLE.  DISTRIBUTE FLOWS INTO GROUP
                if( HP(L) >= HDRY .or. ISDRY == 0 )then
                  SPLIT = HP(L)*DXYP(L)/QSUM1                 ! *** FLOW SPLIT
                  do K = KSZ(L),KC
                    QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                    QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                    QSS(K,LL)      = QLAYER1
                    QSERCELL(K,LL) = QLAYER2
                    QSUM(L,K)      = QSUM(L,K) + QLAYER1 + QLAYER2
                  enddo
                  BCPS(LL).QFACTOR = SPLIT
                else
                  BCPS(LL).QFACTOR = 0.0
                endif
                IU = 1
              endif

            elseif( QSUM2 > 0. )then
              ! *** NO WET CELLS AVAILABLE.  JUST USE THE DRY VOLUMES
              SPLIT = HP(L)*DXYP(L)/QSUM2                     ! *** FLOW SPLIT
              do K = KSZ(L),KC
                QLAYER1 = SPLIT*QSUM3(1,IG)*DBLE(DZC(L,K))
                QLAYER2 = SPLIT*QSUM3(2,IG)*DBLE(DZC(L,K))
                QSS(K,LL)      = QLAYER1
                QSERCELL(K,LL) = QLAYER2
                QSUM(L,K)      = QSUM(L,K) + QLAYER1 + QLAYER2
              enddo
              BCPS(LL).QFACTOR = SPLIT
            endif
          endif
        enddo

        ! *** REPORT WHEN CHANGING FLOW SPLITS
        if( NITER == 1 .or. IU == 1 )then
          PRINT '(A,F10.4,I8,3I5,2E14.6)', 'Change in flow splits: ',TIMEDAY, NITER, IG, NQWET(IG), MWET, QSUM4(IG), QSUM1
        endif

        ! *** save CURRENT WET CELL CONDITIONS
        NQWET(IG) = MWET
        if( IU == 1 ) QSUM4(IG) = QSUM1

      endif
    enddo
  endif

  ! *************************************************************************************
  ! ***  GROUNDWATER SOURCE/SINK INTERPOLATION
  if( NGWSER >= 1 .and. GWSP.IFLAG == 0 )then
    NCTMP = 3 + NDYM + NTOX + NSED + NSND
    do NS = 1,NGWSER
      CTIM = TIMESEC2/TCGWSER(NS)
      
      M2 = MTSGWLAST(NS)
      do while ( CTIM > TGWSER(M2,NS) )
        M2 = M2+1
        if( M2 > MGWSER(NS) )then
          M2 = MGWSER(NS)
          exit
        endif
      enddo
      MTSGWLAST(NS) = M2
      M1 = M2-1

      TDIFF = TGWSER(M2,NS) - TGWSER(M1,NS)
      WTM1 = (TGWSER(M2,NS)-CTIM)/TDIFF
      WTM2 = (CTIM-TGWSER(M1,NS))/TDIFF
      GWSERT(NS) = WTM1*GWSER(M1,NS) + WTM2*GWSER(M2,NS)
      do NC = 1,NCTMP
        GWCSERT(NS,NC) = WTM1*GWCSER(M1,NS,NC) + WTM2*GWCSER(M2,NS,NC)
        if( GWCSERT(NS,NC) < 0.0 )then
          PRINT '(A,I10,F10.4,2I6)','GWCSERT(NS,NC) is < 0.0!', niter, timeday, ns, nc  
          PAUSE
        endif
      enddo
    enddo

  endif

  ! *************************************************************************************
  ! *** CONTROL STRUCTURES AND TIDAL INLETS

  ! -------------------------------------------------------------------------------------
  ! *** SET CURRENT TIME QSFACTOR BASED ON ELEVATIONS
  if( LGROUPS )then
    QSFACTTOT = 0.
    do NCTL = 1,NQCTL
      if( HYD_STR(NCTL).NQCTYP == -2 )then
        IU = HYD_STR(NCTL).IQCTLU
        JU = HYD_STR(NCTL).JQCTLU
        LU = LIJ(IU,JU)
        if( HP(LU) >= HDRY )then
          QSFACTTOT(HYD_STR(NCTL).QCTLGRP) = QSFACTTOT(HYD_STR(NCTL).QCTLGRP) + HYD_STR(NCTL).QCTLMU
        endif
      endif
    enddo

    ! *** NOW COMPUTE CURRENT TIME STEP FLOW SPLITS
    do NCTL = 1,NQCTL
      if( HYD_STR(NCTL).NQCTYP == -2 )then
        if( QSFACTTOT(HYD_STR(NCTL).QCTLGRP) > 0. )then
          QSFACT(NCTL) = HYD_STR(NCTL).QCTLMU/QSFACTTOT(HYD_STR(NCTL).QCTLGRP)
        else
          QSFACT(NCTL) = 0.
        endif
      endif
    enddo

    ! *** UPSTREAM ELEVATION DEPENDANT FLOWS
    do NCTL = 1,NQCTL
      if( HYD_STR(NCTL).NQCTYP == -2 )then
        NCTLT = HYD_STR(NCTL).NQCTLQ  ! *** SERIES INDEX
        IU = HYD_STR(NCTL).IQCTLU
        JU = HYD_STR(NCTL).JQCTLU
        LU = LIJ(IU,JU)
        HUP = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
        HUP = HUP*HCTLUM(NCTLT)

        if( HP(LU) < HWET )then
          do K = 1,KC
            QCTLT(K,NCTL,1) = 0.
          enddo
        else
          ! *** SUFFICIENT DEPTH, GET FLOWS
          M1 = 0
          M2 = 1
      500 M1 = M1+1
          M2 = M2+1
          if( M2 > MQCTL(NCTLT) )then
            write(6,*) ' US ELEVATION CONTROL STRUCTURE TABLE OUT OF BOUNDS'
            write(6,'(" NCTL,NCTLT,IU,JU = ",4I5)') NCTL,NCTLT,IU,JU
            write(6,'(" WSEL,HU = ",2F10.3)') HUP,HP(LU)
            open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            write(mpi_efdc_out_unit,*) ' US ELEVATION CONTROL STRUCTURE TABLE OUT OF BOUNDS'
            write(mpi_efdc_out_unit,'(" NCTL,NCTLT,IU,JU = ",4I5)') NCTL,NCTLT,IU,JU
            write(mpi_efdc_out_unit,'(" WSEL,HU = ",2F10.3)') HUP,HP(LU)
            close(mpi_efdc_out_unit)
            call STOPP('.')
          endif
          if( HUP >= HDIFCTL(M1,NCTLT) .and. HUP <= HDIFCTL(M2,NCTLT) )then
            TDIFF = HDIFCTL(M2,NCTLT) - HDIFCTL(M1,NCTLT)
            WTM1 = (HDIFCTL(M2,NCTLT) - HUP)/TDIFF
            WTM2 = (HUP - HDIFCTL(M1,NCTLT))/TDIFF
            do K = 1,KC
              QCTLT(K,NCTL,1) = WTM1*QCTL(M1,1,K,NCTLT) + WTM2*QCTL(M2,1,K,NCTLT)
              QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*QSFACT(NCTL)
            enddo
          elseif( HUP < HDIFCTL(1,NCTLT) )then
            ! *** Water levels below rating curve.  Assume zero flow
            do K = 1,KC
              QCTLT(K,NCTL,1) = 0.
            enddo
          else
            GOTO 500
          endif
        endif

        ! *** APPLY RQCMUL FLOW SCALING
        do K = 1,KC
          QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*HYD_STR(NCTL).RQCMUL
        enddo
      endif   ! **** HYD_STR(NCTL).NQCTYP == -2
    enddo     ! *** END OF NQCTL LOOP
  endif

  ! -------------------------------------------------------------------------------------
  ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE
  do NCTL = 1,NQCTL
    if( HYD_STR(NCTL).NQCTYP >= -1 .and. HYD_STR(NCTL).NQCTYP <= 1 )then
      ! *** UPSTREAM ELEVATION/STAGE (0 & 1) OR DEPTH (-1) DEPENDANT FLOWS
      NCTLT = HYD_STR(NCTL).NQCTLQ  ! *** SERIES INDEX

      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      HUP = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
      if( HYD_STR(NCTL).NQCTYP == -1 ) HUP = HP(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU

      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      if( ID <= 1 .or. JD <= 1 )then
        LD = 1
        HDW = 0.
      else
        LD = LIJ(ID,JD)
        HDW = HP(LD) + BELV(LD) + HCTLDA(NCTLT) + HYD_STR(NCTL).HQCTLD
        if( HYD_STR(NCTL).NQCTYP == -1 )then
          ! *** UPSTREAM DEPTH ONLY
          HDW = 0.
        endif
      endif

      DELH = HCTLUM(NCTLT)*HUP - HCTLDM(NCTLT)*HDW
      if( DELH <= 0. .or. HP(LU) < HWET )then
        ! *** PREVENT REVERSE FLOW OR INSUFFICENT DEPTH
        do K = 1,KC
          QCTLT(K,NCTL,1) = 0.
        enddo
      else
        ! *** SUFFICIENT DEPTH
        if( HYD_STR(NCTL).NQCTYP == 1 ) DELH = SQRT(DELH)
        M1 = 0
        M2 = 1
    600 M1 = M1+1
        M2 = M2+1
        if( M2 > MQCTL(NCTLT) )then
          write(6,6666)
          write(6,6667)NCTL,NCTLT,IU,JU,ID,JD
          write(6,6668)HUP,HP(LU),HDW,HP(LD)
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit,6666)
          write(mpi_efdc_out_unit,6667)NCTL,NCTLT,IU,JU,ID,JD
          write(mpi_efdc_out_unit,6668)HUP,HP(LU),HDW,HP(LD)
          close(mpi_efdc_out_unit)
          call STOPP('.')
        endif
        if( DELH >= HDIFCTL(M1,NCTLT) .and. DELH <= HDIFCTL(M2,NCTLT) )then
          TDIFF = HDIFCTL(M2,NCTLT) - HDIFCTL(M1,NCTLT)
          WTM1 = (HDIFCTL(M2,NCTLT) - DELH)/TDIFF
          WTM2 = (DELH - HDIFCTL(M1,NCTLT))/TDIFF
          do K = 1,KC
            QCTLT(K,NCTL,1) = WTM1*QCTL(M1,1,K,NCTLT) + WTM2*QCTL(M2,1,K,NCTLT)
          enddo
        else
          GOTO 600
        endif
      endif

      if( HYD_STR(NCTL).NQCTYP == 1 )then
        ! *** ACCELERATION TYPE
        if( AQCTL(NCTLT) > 0.0 )then
          if( ISTL == 3 )then
            do K = 1,KC
              QCTLST(K,NCTL) = QCTLT(K,NCTL,1)
              TMPVAL = QCTLTO(K,NCTL) + DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
              QCTLT(K,NCTL,1) = TMPVAL/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
              QCTLTO(K,NCTL) = QCTLT(K,NCTL,1)
              QCTLSTO(K,NCTL) = QCTLST(K,NCTL)
            enddo
          else
            do K = 1,KC
              QCTLST(K,NCTL) = QCTLT(K,NCTL,1)
              TMPVAL = QCTLTO(K,NCTL) + DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
              QCTLT(K,NCTL,1) = TMPVAL/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
              QCTLT(K,NCTL,1) = 0.5*(QCTLT(K,NCTL,1)+QCTLTO(K,NCTL))
            enddo
          endif
        endif
      endif

      ! *** APPLY RQCMUL FLOW SCALING
      do K = 1,KC
        QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*HYD_STR(NCTL).RQCMUL
      enddo

    endif
  enddo

  ! -------------------------------------------------------------------------------------
  ! *** UPSTREAM/DOWNSTREAM ELEVATION RATING CURVE  HYD_STR(NCTL).NQCTYP == 2
  do NCTL = 1,NQCTL
    if( HYD_STR(NCTL).NQCTYP == 2 )then
      NCTLT = HYD_STR(NCTL).NQCTLQ

      ! *** UPSTREAM ELEVATION
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      HUP = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU

      ! *** CHECK FOR ACTIVE FLOW CONDITIONS
      if( HUP < HDIFCTL(1,NCTLT) .or. HP(LU) < HWET )then
        ! *** INSUFFICIENT DEPTHS
        do K = 1,KC
          QCTLT(K,NCTL,1) = 0.
        enddo
        LD = 1
      else
        ! *** SUFFICIENT DEPTHS
        ! *** DOWNSTREAM ELEVATION
        ID = HYD_STR(NCTL).IQCTLD
        JD = HYD_STR(NCTL).JQCTLD
        LD = LIJ(ID,JD)
        HDW = HP(LD) + BELV(LD) + HCTLDA(NCTLT) + HYD_STR(NCTL).HQCTLD
        HTMPD = HDIFCTD(1,NCTLT)+0.001
        HDW = max(HDW,HTMPD)

        MU1 = 0
        MU2 = 1
        MD1 = 0
        MD2 = 1

        ! *** FIND UPSTREAM ELEVATION BRACKET
    700 MU1 = MU1+1
        MU2 = MU1+1
        if( MU2 > MQCTL(NCTLT) )then
          write(6,6676)
          write(6,6677)NCTL,NCTLT,IU,JU,ID,JD
          write(6,6678)HUP,HP(LU),HDW,HP(LD)
          write(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit,6676)
          write(mpi_efdc_out_unit,6677)NCTL,NCTLT,IU,JU,ID,JD
          write(mpi_efdc_out_unit,6678)HUP,HP(LU),HDW,HP(LD)
          write(mpi_efdc_out_unit,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          close(mpi_efdc_out_unit)
          call STOPP('.')
        endif
        if( HUP >= HDIFCTL(MU1,NCTLT) .and. HUP <= HDIFCTL(MU2,NCTLT) )then
          ! *** FOUND VALID UPSTREAM ELEVATION RANGE
          TDIFFU = HDIFCTL(MU2,NCTLT)-HDIFCTL(MU1,NCTLT)
          WTM1U = (HDIFCTL(MU2,NCTLT)-HUP)/TDIFFU
          WTM2U = (HUP-HDIFCTL(MU1,NCTLT))/TDIFFU
        else
          GOTO 700
        endif

        ! *** FIND DOWNSTREAM ELEVATION BRACKET
    750 MD1 = MD1+1
        MD2 = MD1+1
        if( MD2 > MQCTL(NCTLT) )then
          write(6,6686)
          write(6,6687)NCTL,NCTLT,IU,JU,ID,JD
          write(6,6688)HUP,HP(LU),HDW,HP(LD)
          write(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit,6686)
          write(mpi_efdc_out_unit,6687)NCTL,NCTLT,IU,JU,ID,JD
          write(mpi_efdc_out_unit,6688)HUP,HP(LU),HDW,HP(LD)
          write(mpi_efdc_out_unit,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
          close(mpi_efdc_out_unit)
          call STOPP('.')
        endif
        if( HDW >= HDIFCTD(MD1,NCTLT) .and. HDW <= HDIFCTD(MD2,NCTLT) )then
          ! *** FOUND VALID DOWNSTREAM ELEVATION RANGE
          TDIFFD = HDIFCTD(MD2,NCTLT)-HDIFCTD(MD1,NCTLT)
          WTM1D = (HDIFCTD(MD2,NCTLT)-HDW)/TDIFFD
          WTM2D = (HDW-HDIFCTD(MD1,NCTLT))/TDIFFD
        else
          GOTO 750
        endif

        ! *** DETERMINE FLOWS BASED ON UPSTREAM AND DOWNSTREAM ELEVATIONS
        do K = 1,KC
          QCTLT(K,NCTL,1) = WTM1U*( WTM1D*QCTL(MU1,MD1,K,NCTLT) + WTM2D*QCTL(MU1,MD2,K,NCTLT) ) &
            + WTM2U*( WTM1D*QCTL(MU2,MD1,K,NCTLT) + WTM2D*QCTL(MU2,MD2,K,NCTLT) )
        enddo
      endif

      ! *** APPLY RQCMUL FLOW SCALING
      do K = 1,KC
        QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1)*HYD_STR(NCTL).RQCMUL
      enddo
    endif
  enddo

  ! -------------------------------------------------------------------------------------
  ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE WITH A LOWER CHORD ACTIVATION TOGGLE
  do NCTL = 1,NQCTL
    if( HYD_STR(NCTL).NQCTYP == 3 .or. HYD_STR(NCTL).NQCTYP == 4 )then
      ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE DIFFERENCE (4) DEPENDANT FLOWS
      NCTLT = HYD_STR(NCTL).NQCTLQ

      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)

      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      if( ID == 0 .and. JD == 0 )then
        ! *** INVALID SPECIFICATION
        CYCLE
      endif
      LD = LIJ(ID,JD)

      ! *** SKIP WHEN INSUFFICENT WATER
      if( HP(LU) < HWET ) CYCLE

      HUP = HP(LU) + BELV(LU) ! *** WSEL
      if( HUP <= HYD_STR(NCTL).BQCLCE .and. (N < HYD_STR(NCTL).NQCMINS .or. NLOWCHORD(NCTL) < 0 .or. NLOWCHORD(NCTL) > HYD_STR(NCTL).NQCMINS/2) )then
        ! *** WATER BELOW LOW CHORD

        ! *** SET THE CELL FACE FLAGS - ALLOW FULL HYDRODYNAMICS
        if( LOWCHORDU(NCTL) /= -9999. )then
          write(*,*)' *** LOWER CHORD BC: OFF   NCTL = ',NCTL

          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit,*)' *** LOWER CHORD BC: OFF   NCTL = ',NCTL
          write(mpi_efdc_out_unit,*)'                           TIMEDAY = ',TIMEDAY
          write(mpi_efdc_out_unit,*)'                           LU,LD = ',LU,LD
          write(mpi_efdc_out_unit,*)'                           HUP,HDP = ',HUP,HP(LD)+BELV(LD)
          write(mpi_efdc_out_unit,*)' ***'
          close(mpi_efdc_out_unit)

          do K = 1,KC
            QCTLT(K,NCTL,1) = 0.
          enddo
          if( ID > IU )then
            if( SAVESUB(2,NCTL) > 0.5 )then
              SUB(LD)   = 1.0
              SUBO(LD)  = 1.0
              UHDYE(LD) = UHDYE(LU)
              TBX(LD)   = TBX(LU)
              TSX(LD)   = TSX(LU)
            endif
          else
            if( SAVESUB(1,NCTL) > 0.5 )then
              SUB(LU)   = 1.0
              SUBO(LU)  = 1.0
              UHDYE(LU) = UHDYE(LD)
              TBX(LU)   = TBX(LD)
              TSX(LU)   = TSX(LD)
            endif
          endif
          if( JD > JU )then
            if( SAVESVB(2,NCTL) > 0.5 )then
              SVB(LD)   = 1.0
              SVBO(LD)  = 1.0
              VHDXE(LD) = VHDXE(LU)
              TBY(LD)   = TBY(LU)
              TSY(LD)   = TSY(LU)
            endif
          else
            if( SAVESVB(1,NCTL) > 0.5 )then
              SVB(LU) = 1.0
              SVBO(LU) = 1.0
              VHDXE(LU) = VHDXE(LD)
              TBY(LU)   = TBY(LD)
              TSY(LU)   = TSY(LD)
            endif
          endif
          LOWCHORDU(NCTL) = -9999.
          LOWCHORDV(NCTL) = -9999.
          NLOWCHORD(NCTL) = 0
          IPMC = 0
        endif
        NLOWCHORD(NCTL) = NLOWCHORD(NCTL)-1
        IPMC = 0
        CYCLE

      else

        ! *** ELEVATION ABOVE LOW CHORD
        if( LOWCHORDU(NCTL) == -9999. )then
          if( NLOWCHORD(NCTL) < -HYD_STR(NCTL).NQCMINS )then
            ! *** SET THE CELL FACE FLAGS - BLOCK FULL HYDRODYNAMICS

            write(*,*)' *** LOWER CHORD BC: ON    NCTL = ',NCTL

            if( ID > IU )then
              SUB(LD) = 0.0
              SUBO(LD) = 0.0
              ! *** save THE FLOWS
              LOWCHORDU(NCTL) = ABS(UHDYE(LD))
              UHDYE(LD) = 0.0
            elseif( ID < IU )then
              SUB(LU) = 0.0
              SUBO(LU) = 0.0
              ! *** save THE FLOWS
              LOWCHORDU(NCTL) = ABS(UHDYE(LU))
              UHDYE(LU) = 0.0
            else
              LOWCHORDU(NCTL) = 0.
            endif
            if( JD > JU )then
              SVB(LD) = 0.0
              SVBO(LD) = 0.0
              ! *** save THE FLOWS
              LOWCHORDV(NCTL) = ABS(VHDXE(LD))
              VHDXE(LD) = 0.0
            elseif( JD < JU )then
              SVB(LU) = 0.0
              SVBO(LU) = 0.0
              ! *** save THE FLOWS
              LOWCHORDV(NCTL) = ABS(VHDXE(LU))
              VHDXE(LU) = 0.0
            else
              LOWCHORDV(NCTL) = 0.
            endif
            open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            write(mpi_efdc_out_unit,*)' *** LOWER CHORD BC: ON    NCTL = ',NCTL
            write(mpi_efdc_out_unit,*)'                           TIMEDAY = ',TIMEDAY
            write(mpi_efdc_out_unit,*)'                           LU,LD = ',LU,LD
            write(mpi_efdc_out_unit,*)'                           HUP,HDP = ',HUP,HP(LD)+BELV(LD)
            write(mpi_efdc_out_unit,*)'                           QU,QV = ',LOWCHORDU(NCTL),LOWCHORDV(NCTL)
            write(mpi_efdc_out_unit,*)' ***'
            close(mpi_efdc_out_unit)

            NLOWCHORD(NCTL) = 0
            IPMC = 0
          else
            NLOWCHORD(NCTL) = NLOWCHORD(NCTL)-1
            CYCLE
          endif
        endif

        if( HYD_STR(NCTL).NQCTYP == 3 )then
          ! *** UPSTREAM DEPTH ONLY
          HUP  = HP(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
          DELH = HCTLUM(NCTLT)*HUP
        else
          ! *** ELEVATION DIFFERENCE
          HUP  = HP(LU) + BELV(LU) + HCTLUA(NCTLT) + HYD_STR(NCTL).HQCTLU
          HDW  = HP(LD) + BELV(LD) + HCTLDA(NCTLT) + HYD_STR(NCTL).HQCTLD
          DELH = HCTLUM(NCTLT)*HUP - HCTLDM(NCTLT)*HDW
        endif
        DELH = max( DELH,HDIFCTL(1,NCTLT) )

        ! *** LOOKUP THE FLOWS
        M1 = 0
        M2 = 1
    800 M1 = M1+1
        M2 = M2+1
        if( M2 > MQCTL(NCTLT) )then
          write(6,6666)
          write(6,6667)NCTL,NCTLT,IU,JU,ID,JD
          write(6,6668)HUP,HP(LU),HDW,HP(LD)
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit,6666)
          write(mpi_efdc_out_unit,6667)NCTL,NCTLT,IU,JU,ID,JD
          write(mpi_efdc_out_unit,6668)HUP,HP(LU),HDW,HP(LD)
          close(mpi_efdc_out_unit)
          call STOPP('.')
        endif
        if( DELH >= HDIFCTL(M1,NCTLT) .and. DELH <= HDIFCTL(M2,NCTLT) )then
          TDIFF = HDIFCTL(M2,NCTLT) - HDIFCTL(M1,NCTLT)
          WTM1 = (HDIFCTL(M2,NCTLT) - DELH)/TDIFF
          WTM2 = (DELH - HDIFCTL(M1,NCTLT))/TDIFF
          do K = 1,KC
            QCTLT(K,NCTL,1) = WTM1*QCTL(M1,1,K,NCTLT) + WTM2*QCTL(M2,1,K,NCTLT)
            ! *** SUBTRACT LOW CHORD FLOWS
            QCTLT(K,NCTL,1) = max(QCTLT(K,NCTL,1) - QLAYER(NCTL,K),0.)
          enddo
          do K = 1,KC
            ! ***                                 ** TOTAL FLOW OPEN CHANNEL FLOW  **
            QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1) + DZC(LU,K)*(LOWCHORDU(NCTL)+LOWCHORDV(NCTL))
          enddo
        else
          GOTO 800
        endif

        NLOWCHORD(NCTL) = NLOWCHORD(NCTL)+1
      endif
    endif
  enddo   ! *** END OF LOW/HI CHORD CHECK

  ! -------------------------------------------------------------------------------------
  ! *** COMPUTED FLOWS USING EQUATIONS OR DO FINAL UPDATES FOR OTHER TYPES
  do NCTL = 1,NQCTL
    if( HYD_STR(NCTL).NQCTYP > 4 )then
      ! *** COMPUTED FLOWS USING HYDRAULIC STRUCTURE EQUATIONS:  HYD_STR(NCTL).NQCTYP > 4
      call COMPUTE_HSFLOW(NCTL)

    else
      ! *** ADJUST FOR INACTIVE LAYERS AND ADD FINAL LAYERS FLOWS TO QSUM
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)

      ! *** CHECK VALID LAYERS
      QSUM1 = 0.
      do K = 1,KSZ(LU)-1
        QSUM1 = QSUM1 + QCTLT(K,NCTL,1)
        QCTLT(K,NCTL,1) = 0.
      enddo

      ! *** ADD FLOWS BELOW BOTTOM ACTIVE LAYER BACK TO ACTIVE LAYERS
      if( QSUM1 > 0. )then
        do K = KSZ(LU),KC
          QCTLT(K,NCTL,1) = QCTLT(K,NCTL,1) + QSUM1*DZC(LU,K)
        enddo
      endif

      ! *** LIMIT OUTFLOWS TO AVAILABLE WATER
      QCTLMAX = (HP(LU) - HDRY)*DXYP(LU)/(DELT*DZI)
      if( QCTLMAX < 0. ) QCTLMAX = 0.
      do K = 1,KC
        QCTLT(K,NCTL,1) = min(QCTLT(K,NCTL,1),QCTLMAX)
      enddo

      ! *** FINAL ADJUSTED LAYER FLOWS ADDED TO TOTAL FLOWS
      do K = KSZ(LU),KC
        QSUM(LU,K) = QSUM(LU,K) - QCTLT(K,NCTL,1)
      enddo

      ! *** DOWNSTREAM
      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      if( ID /= 0 .and. JD /= 0 )then
        LD = LIJ(ID,JD)
        QSUM2 = 0.
        do K = 1,KSZ(LD)-1
          QSUM2 = QSUM2 + QCTLT(K,NCTL,1)
        enddo

        ! *** SET THE RETURN FLOW VARIABLE FOR LOOKUP TABLE CONTROL TYPE STRUCTURES
        do K = KSZ(LD),KC
          QCTLT(K,NCTL,2) = QCTLT(K,NCTL,1) + QSUM2*DZC(LD,K)
          QSUM(LD,K) = QSUM(LD,K) + QCTLT(K,NCTL,2)
        enddo
      endif
    endif
  enddo
  ! *************************************************************************************

  NTMP = 3 + NDYM + NTOX + NSED + NSND
  if( ISTRAN(8) > 0 ) NTMP = NTMP + NWQV

  ! *************************************************************************************
  ! ***  FLOW WITHDRAWAL AND RETURN
  do NC = 1,NTMP
    CQWRSERT(0,NC) = 0.0
  enddo
  do NS = 1,NQWRSR
    CTIM = TIMESEC2/TSWR(NS).TMULT
    
    M2 = MTSWRLAST(NS)
    do while( CTIM > TSWR(NS).TIM(M2))
      M2 = M2+1
      if( M2 > TSWR(NS).NREC )then
        M2 = TSWR(NS).NREC
        exit
      endif
    enddo
    MTSWRLAST(NS) = M2
    M1 = M2-1
    TDIFF = TSWR(NS).TIM(M2) - TSWR(NS).TIM(M1)
    WTM1 = (TSWR(NS).TIM(M2) - CTIM)/TDIFF
    WTM2 = (CTIM - TSWR(NS).TIM(M1))/TDIFF
    
    QWRSERT(NS) = WTM1*TSWR(NS).VAL(M1,0) + WTM2*TSWR(NS).VAL(M2,0) 
    do NC = 1,NTMP
      CQWRSERT(NS,NC) = WTM1*TSWR(NS).VAL(M1,NC) + WTM2*TSWR(NS).VAL(M2,NC)
    enddo
  enddo

  if( NQWR > 0 )then
    do NWR = 1,NQWR
      NS = WITH_RET(NWR).NQWRSERQ
      if( WITH_RET_CTL(NWR).ITYPE > 0 .and. WITH_RET_CTL(NWR).ITYPE < 3 )then
        ! *** Withdrawal/Return is controlled by operation rules
        LU = LIJ(WITH_RET_CTL(NWR).IREFUP, WITH_RET_CTL(NWR).JREFUP)
        ZHU = BELV(LU) + HP(LU)

        if( WITH_RET_CTL(NWR).ITYPE == 2 )then
          ! *** W/R CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
          LD = LIJ(WITH_RET_CTL(NWR).IREFDN, WITH_RET_CTL(NWR).JREFDN)
          ZHD = BELV(LD) + HP(LD)
          HVAL = ZHU - ZHD
        else
          ! *** W/R IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
          HVAL = ZHU
        endif
        CVAL = WITH_RET_CTL(NWR).CUR.FLOW
        call PUMP_OPERATION_RULES(WITH_RET_CTL(NWR), HVAL, CVAL, RVAL)
        QWRSERT(NS) = RVAL
      endif

      ! *** Handle +/- Flows for Withdrawal/Return Structures
      if( QWRSERT(NS) >= 0. )then
        ! *** Original Withdrawal/Return
        IU = WITH_RET(NWR).IQWRU
        JU = WITH_RET(NWR).JQWRU
        KU = WITH_RET(NWR).KQWRU
        ID = WITH_RET(NWR).IQWRD
        JD = WITH_RET(NWR).JQWRD
        KD = WITH_RET(NWR).KQWRD
      else
        ! *** Reverse Flow Withdrawal/Return
        ID = WITH_RET(NWR).IQWRU
        JD = WITH_RET(NWR).JQWRU
        KD = WITH_RET(NWR).KQWRU
        IU = WITH_RET(NWR).IQWRD
        JU = WITH_RET(NWR).JQWRD
        KU = WITH_RET(NWR).KQWRD
        WITH_RET(NWR).QWR = 0.  ! *** Only allow time variable flows when using! -W/R
      endif
      LU = LIJ(IU,JU)
      LD = LIJ(ID,JD)
      QWRABS = ABS(QWRSERT(NS))

      QSUM(LU,KU) = QSUM(LU,KU) - WITH_RET(NWR).QWR - QWRABS
      QSUM(LD,KD) = QSUM(LD,KD) + WITH_RET(NWR).QWR + QWRABS
      WITH_RET_CTL(NWR).CUR.FLOW = QWRSERT(NS)
    enddo
  endif

  ! ***  call JPEFDC AND PLACE JET-PLUME VOLUMES SOURCES
  if( NQJPIJ > 0 .and. N == 1 ) CALL JPEFDC
  if( NQJPIJ > 0 .and. ISTL == 3 )then
    if( NUDJPC >= NUDJP )then
      call JPEFDC
      NUDJPC = 1
    else
      NUDJPC = NUDJPC + 1
    endif
  endif
  if( NQJPIJ > 0 .and. IS2TIM >= 1 )then
    if( NUDJPC >= NUDJP )then
      call JPEFDC
      NUDJPC = 1
    else
      NUDJPC = NUDJPC + 1   ! *** delme - need to convert to time (seconds) based, not N
    endif
  endif

  ! *** Place JET-PLUME volume sources
  if( NQJPIJ > 0 )then
    do NJP = 1,NQJPIJ
      ! *** QVJPTMP = JET-PLUME discharge per port
      ! *** QJPENT  = Ambient water entrainment volumetric rate by layer (negative for inflow, positive for outflow)
      ! *** QJPENTT = Ambient water entrainment volumetric rate water column total.  Used for ensuring mass balance.
      if( JET_PLM(NJP).ICALJP == 1 )then
        ! *** Discharge only using qser and standard concentration series
        RPORTS = FLOAT(JET_PLM(NJP).NPORTJP)
        LJP = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)
        KTMP = KEFFJP(NJP)
        
        QVJPTMP  = JET_PLM(NJP).QQCJP         ! *** QQCJP   = Constant flow per cell
        do K = 1,KC
          QVJPTMP = QVJPTMP + QSERT(K,JET_PLM(NJP).NQSERJP)
        enddo

        ! *** Remove the entrainment from each layer
        do K = KSZ(LJP),KC
          QSUM(LJP,K) = QSUM(LJP,K) - RPORTS*QJPENT(K,NJP)
        enddo
        
        ! *** Place discharge and total entrainment at effective layer
        QSUM(LJP,KTMP) = QSUM(LJP,KTMP) + RPORTS*(QVJPTMP + QJPENTT(NJP))
      endif
      if( JET_PLM(NJP).ICALJP == 2 )then
        ! *** JET-PLUME type using withdrawal/return series
      
        ! *** Downstream/JET-PLUME discharge
        RPORTS = FLOAT(JET_PLM(NJP).NPORTJP)
        LJP = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)
        KTMP = KEFFJP(NJP)

        ! *** QVJPTMP = JET-PLUME discharge per port.  QWRCJP = Constant flow per cell
        QVJPTMP = JET_PLM(NJP).QWRCJP + QWRSERT(JET_PLM(NJP).NQWRSERJP)

        ! *** Remove entrainment from each layer
        do K = KSZ(LJP),KC
          QSUM(LJP,K) = QSUM(LJP,K) - RPORTS*QJPENT(K,NJP)
        enddo

        ! *** Place discharge and total entrainment at effective layer
        QSUM(LJP,KTMP) = QSUM(LJP,KTMP) + RPORTS*(QVJPTMP + QJPENTT(NJP))

        ! *** Remove discharge from upstream/intake cell
        LU = LIJ(JET_PLM(NJP).IUPCJP,JET_PLM(NJP).JUPCJP)
        KU = JET_PLM(NJP).KUPCJP
        QSUM(LU,KU) = QSUM(LU,KU) - RPORTS*QVJPTMP
      endif
    enddo
  endif

  ! *** Compute ice met/growth volumetric rates (m3/s)
  if( ISICE > 2 .and. N > 0 .and. (IS2TIM > 0 .or. ISTL == 3) )then
    if( LCHECKICE .or. LFRAZIL )then
      if( ISTL == 3 )then
        DELTICE1 = DELTICE1 + 0.5*DELT   !*** DKT - Use the same time increment as in 2TL
        DELTICE2 = DELTICE2 + 0.5*DELT
      else                  
        DELTICE1 = DELTICE1 + DELT
        DELTICE2 = DELTICE2 + DELT
      endif

      call ICECOMP(DELTICE2, ICESTEP)
      if( DELTICE2 >= ICESTEP ) DELTICE2 = 0.0  ! *** Reset ice cover timer and volume accumulator
    endif
  endif
  
  ! *** Update time variable groundwater field (must include all cells)
  if( GWSP.IFLAG > 0 )then
    call UPDATEFIELD(GWSP,TIMEDAY,1,QGW)
  endif
  if( ISTL == 3 )then
    TMPVAL = 0.5*999.8426/RHOI*DELT
  else
    TMPVAL = 999.8426/RHOI*DELT
  endif
  
  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** DETERMIME GROUNDWATER FLUXES
  if( ISGWIT == 2 .and. GWSP.IFLAG == 0 )then
    if( NGWSER == 1 )then
      if( IGWSER(1) == 0 )then
        ! *** FLOW SERIES
        !$OMP DO PRIVATE(ND,LF,LL,L)
        do ND = 1,NDM
          LF = 2+(ND-1)*LDM
          LL = min(LF+LDM-1,LA)
          do L = LF,LL
            QGW(L) = GWFAC(L)*GWSERT(NGWSL(L))
          enddo
        enddo
        !$OMP END DO
      else
        ! *** SEEPAGE VELOCITY SERIES
        !$OMP DO PRIVATE(ND,LF,LL,L)
        do ND = 1,NDM
          LF = 2+(ND-1)*LDM
          LL = min(LF+LDM-1,LA)
          do L = LF,LL
            QGW(L) = GWFAC(L)*DXYP(L)*GWSERT(NGWSL(L))
          enddo
        enddo
        !$OMP END DO
      endif
    else
      !$OMP DO PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          if( IGWSER(NGWSL(L)) == 0 )then
            ! *** FLOW SERIES
            QGW(L) = GWFAC(L)*GWSERT(NGWSL(L))
          else
            ! *** SEEPAGE VELOCITY SERIES
            QGW(L) = GWFAC(L)*DXYP(L)*GWSERT(NGWSL(L))
          endif
        enddo
      enddo
      !$OMP END DO
    endif
  endif

  ! *** ACCUMULATE GROUNDWATER FLUXES
  if( ISGWIT > 0 )then
    !$OMP DO PRIVATE(ND,LF,LL,L)
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do L = LF,LL
        QSUM(L,KSZ(L)) = QSUM(L,KSZ(L)) + QGW(L)
      enddo
    enddo
    !$OMP END DO
  endif

  ! *** ADD ICE MET/GROWTH VOLUMETRIC RATES (M3/S)
  if( ISICE > 2 .and. ( LCHECKICE .or. LFRAZIL ) .and. DELTICE1 >= 60. )then
    !$OMP DO PRIVATE(ND,LF,LL,L,HPICE) REDUCTION(+:NLIMITICE)
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do L = LF,LL

        ! *** DUE TO MACHINE PRECISION ISSUES, ONLY APPLY ICE MET/FREEZE FLOW RATES WHEN LARGE ENOUGH
        HPICE = ICEVOL(L)*DXYIP(L)
        if( ABS(HPICE) > HP(L)*1.E-6 )then
          ! *** LIMIT ICE GROWTH TO AVAILABLE WATER
          if( HPICE < 0. .and. ISDRY > 0 )then
            if( HP(L)+HPICE < HDRYICE )then
              NLIMITICE = NLIMITICE + 1
              ICEVOL(L) = 0.0
            endif
          endif

          ICERATE(L) = ICEVOL(L)/DELT
          QSUM(L,KC) = QSUM(L,KC) + ICERATE(L)
          ICEVOL(L)  = 0.0
        elseif( ICETHICK(L) < HPICE )then
          ICERATE(L) = ICEVOL(L)/DELT
          QSUM(L,KC) = QSUM(L,KC) + ICERATE(L)
          ICEVOL(L)  = 0.0
        endif
      enddo
    enddo
    !$OMP END DO
    !$OMP SINGLE
    if( TIMEDAY >= ICEDAY .and. NLIMITICE > 0 )then  ! *** ONLY DISPLAY WARNING ONCE PER DAY
      write(6,'(F10.1,A,I5,A)') TIMEDAY, ' : LIMITED ICE PRODUCTION IN CELLS ', NLIMITICE, ' TIMES'
      ICEDAY = INT(TIMEDAY) + 1.
      NLIMITICE = 0
    endif
    DELTICE1 = 0.
    !$OMP END SINGLE
  endif

  ! *** EVAPORATION AND RAINFALL.  BOTH RAINT AND EVAPT ARE IN M/S
  if( NASER > 0  .or. EVAP.IFLAG > 0 .or. RAIN.IFLAG > 0 )then
    if( IEVAP > 1 .or. ISTOPT(2) == 1 )then
      ! *** COMPUTE PARTIAL PRESSURE OF WATER VAPOR AT WATER TEMPERATURE (mb)
      !$OMP DO PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          SVPW(L) = (10.**((0.7859+0.03477*TEM(L,KC))/(1.+0.00412*TEM(L,KC))))
        enddo
      enddo
      !$OMP END DO
    endif

    ! *** PRECIPITATION, TREAT AS ICE, IF BELOW FREEZING AND HAS EXISTING ICE COVER
    IF( ISTRAN(2) > 0 )THEN
      IF( ISICE > 2 .AND. ( ISTL == 3 .OR. IS2TL > 0 ))THEN
        
        !$OMP DO PRIVATE(ND,LF,LL,L)
        do ND = 1,NDM
          LF = 2+(ND-1)*LDM
          LL = min(LF+LDM-1,LA)

          if( LFRAZIL .and. IEVAP == 1 )then
            do L = LF,LL
              if( ICETHICK(L) > 0.0 )then
                ICETHICK(L) = ICETHICK(L) - EVAPT(L)*TMPVAL          ! *** SUBLIMATION
                if( ICETHICK(L) < 0.0 )then
                  ICETHICK(L) = 0.0
                endif
                if( ICETHICK(L) < MINICETHICK ) ICECELL(L) = .FALSE.
                EVAPT(L) = 0.0
              endif
            enddo
          endif

          do L = LF,LL
            if( TATMT(L) < 0.0 )then
              if( ICECELL(L) .or. .not. LMASKDRY(L) )then
                ICETHICK(L) = ICETHICK(L) + RAINT(L)*TMPVAL
                if( ICETHICK(L) < MINICETHICK )then
                  ICECELL(L) = .FALSE.
                else
                  ICECELL(L) = .TRUE.
                endif
                RAINT(L) = 0.0
              endif
            endif
          enddo
        enddo
        !$OMP END DO
      endif
    endif

    if( IEVAP == 0 )then
      ! *** IGNORE EVAPORATION AND RAINFALL
    
    elseif( IEVAP == 1 .or. EVAP.IFLAG > 0 .or. RAIN.IFLAG > 0 )then
      ! *** Use INPUT EVAPORATION
      !$OMP DO PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        enddo
      enddo
      !$OMP END DO

    elseif( IEVAP == 2 )then
      ! *** COMPUTE EVAPORATION (ORIGINAL EFDC APPROACH)
      !$OMP DO PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          EVAPT(L) = CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW(L) - VPAT(L))/PATMT(L)
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        enddo
      enddo
      !$OMP END DO

    elseif( IEVAP > 2 .and. IEVAP < 11 )then
      ! *** COMPUTE EVAPORATION USING THE WIND FUNCTION APPROACH  f(W) = A+B*W+C*W**2
      ! *** THE f(W) UNITS HAVE BEEN CONVERTED TO M/S/millibar
      !$OMP DO PRIVATE(ND,LF,LL,L,CLEVAPTMP)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          CLEVAPTMP  = WINDFA + WINDFB*WINDST(L) + WINDFC*WINDST(L)**2
          EVAPT(L)   = CLEVAPTMP*(SVPW(L)-VPAT(L))
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        enddo
      enddo
      !$OMP END DO

    elseif( IEVAP == 11 )then
      ! *** COMPUTE EVAPORATION: RYAN-HARLEMAN (FROM CE-QUAL-W2, COLE & WELLS, 2011)
      !$OMP DO PRIVATE(ND,LF,LL,L,TVW,TVA,DTV,DTVL,LAMBDA,TM,VPTG)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          TVA  = (TATMT(L) +273.0)/(1.0-0.378*VPAT(L)/PATMT(L))
          TVW  = (TEM(L,KC)+273.0)/(1.0-0.378*SVPW(L)/PATMT(L))
          DTV  = TVW-TVA
          DTVL =  0.0084*WINDST(L)**3
          if( DTV < DTVL) DTV = DTVL
          LAMBDA = 3.59*DTV**0.3333333
          LAMBDA = LAMBDA+4.26*WINDST(L)

          TM       = (TEM(L,KC)+TDEWT(L))*0.5
          VPTG     =  0.35 + 0.015*TM + 0.0012*TM*TM
          EVAPT(L) = VPTG*(TEM(L,KC) - TDEWT(L))*LAMBDA/2.45E9
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        enddo
      enddo
      !$OMP END DO

    elseif( IEVAP == 12 )then
      ! *** Use EVAPORATION BY COARE 3.6
      !$OMP DO PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          EVAPT(L) = EVACOARE(L)
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        enddo
      enddo
      !$OMP END DO

    elseif( IEVAP == 13 )then
      ! *** Compute evaporation using the approach from Arifin et al., 2016
      ! *** This approach combines the adjusment for wind speed from Croley 1989 with the algorithsm in Quinn 1979
      !$OMP DO PRIVATE(ND,LF,LL,L,DENAIR,CLEVAPTMP,FWUP,QSW,QAA)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          DENAIR = PATMT(L)*100./287.05 /(TATMT(L) + 273.15)                          ! *** Density of air (kg/m^3)
          CLEVAPTMP = 8.*10E-4*DENAIR/RHOW(L,KC)
          FWUP = 1.607 + 0.92*1.2475*WINDST(L) - 0.28*(TATMT(L) - TEM(L,KC))          ! *** Adjusted wind speed factor
          QSW = 0.622*SVPW(L) /(PATMT(L) - 0.378*SVPW(L))                             ! *** The specific humidity of water
          QAA = 0.622*VPAT(L) /(PATMT(L) - 0.378*VPAT(L))                             ! *** The specific humidity of air
          
          EVAPT(L)   = CLEVAPTMP*FWUP*(QSW - QAA)
          QSUM(L,KC) = QSUM(L,KC) + DXYP(L)*(RAINT(L)-EVAPT(L))
        enddo
      enddo
      !$OMP END DO

    endif

    ! *** REMOVE EVAP & RAIN FROM CELLS WITH ICE COVER
    if( ISICE > 0 .and. IEVAP > 1 )then
      !$OMP DO PRIVATE(ND,LF,LL,L,TVW,TVA,DTV,DTVL,LAMBDA,TM,VPTG)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          if( ICECELL(L) )then
            ! *** BACK OUT ANY ADDED FLOWS
            QSUM(L,KC) = QSUM(L,KC) - DXYP(L)*( RAINT(L)-EVAPT(L) )
            EVAPT(L)   = 0.
            RAINT(L)   = 0.
          endif
        enddo
      enddo
      !$OMP END DO
    endif

  endif

  ! ***  DETERMINE NET EXTERNAL VOLUME SOURCE/SINK
  !$OMP DO PRIVATE(ND,LF,LL,L,K)
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = min(LF+LDM-1,LA)
    do K = 1,KC
      do L = LF,LL
        QSUME(L) = QSUME(L) + QSUM(L,K)
      enddo
    enddo
  enddo
  !$OMP END DO

  ! *** ADJUST VOLUME SOURCE AND SINKS
  if( ISGWIE == 0 )then
    !$OMP DO PRIVATE(ND,LF,LL,L,K,QEAVAIL,RAVAIL,TMPVAL,TF,I)
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do L = LF,LL
        if( QSUME(L) < 0. )then
          ! *** LIMIT WITHDRAWAL VOLUMES TO AVAILABLE WATER IN CELL
          if( HP(L) < HWET )then
            if( HP(L) <= HDRY9 )then
              ! *** TURN OFF ANY OUTFLOWS
              QGW(L)   = 0.
              EVAPT(L) = 0.
              RAINT(L) = 0.

              ! *** ADJUST LAYER BY LAYER FLOWS
              do K = 1,KC
                QSUM(L,K) = 0.0
              enddo
              QSUME(L) = 0.

              ! *** REDUCE ICE THICKNESS
              if( ISICE > 2 .and. ( ISTL == 3 .or. IS2TL > 0 ) )then
                if( ICEVOL(L) < 0. )then
                  TMPVAL = ICEVOL(L)*999.8426/RHOI/DXYP(L)
                  ICETHICK(L) = ICETHICK(L) + TMPVAL
                  ICETHICK(L) = max(ICETHICK(L),0.)
                  if( ICETHICK(L) < MINICETHICK )ICECELL(L) = .FALSE.
                  ICEVOL(L) = 0.
                endif
              endif
            else
              ! *** CHECK AVAILBILITY OF WATER IN CELL
              QEAVAIL = DXYP(L)*(HP(L)-HDRY9)*DELTI
              if( QEAVAIL < ABS(QSUME(L)) )then
                ! *** LIMIT WITHDRAWAL VOLUMES TO AVAILABLE WATER IN CELL
                QEAVAIL  = -QEAVAIL
                RAVAIL   = QEAVAIL/QSUME(L)

                QSUME(L) = QEAVAIL
                ! *** ADJUST LAYER BY LAYER FLOWS
                do K = 1,KC
                  QSUM(L,K) = QSUM(L,K)*RAVAIL
                enddo

                QGW(L)   = QGW(L)*RAVAIL
                EVAPT(L) = EVAPT(L)*RAVAIL
                RAINT(L) = RAINT(L)*RAVAIL

                ! *** REDUCE ICE THICKNESS
                if( ISICE > 2 .and. ( ISTL == 3 .or. IS2TL > 0 ) )then
                  if( ICEVOL(L) < 0. )then
                    ! *** FREEZING TEMPERATURE OF WATER
                    if( ISTRAN(1) > 0 )then
                      if( SAL(L,KC) < 35. )then
                        TF = -0.0545*SAL(L,KC)
                      else
                        TF = -0.31462-0.04177*SAL(L,KC)-0.000166*SAL(L,KC)*SAL(L,KC)
                      endif
                    else
                      TF = 0.0
                    endif

                    ICEVOL(L) = ICEVOL(L)*RAVAIL
                    TMPVAL = ICEVOL(L)*999.8426/RHOI/DXYP(L)
                    if( ISICE == 3 )then
                      ! *** REDUCE ICE COVER THICKNESS
                      ICETHICK(L) = ICETHICK(L) + TMPVAL
                      ICETHICK(L) = max(ICETHICK(L),0.)
                    else
                      ! *** RESTORE THE FRAZIL ICE
                      if( TMPVAL+ICETHICK(L) > 0. )then
                        ICETHICK(L) = ICETHICK(L) + TMPVAL
                      else
                        FRAZILICE(L,KC) = FRAZILICE(L,KC) + TMPVAL
                        FRAZILICE(L,KC) = max(FRAZILICE(L,KC),0.)
                      endif
                    endif
                    if( ICETHICK(L) < MINICETHICK )then
                      ICECELL(L) = .FALSE.
                    else
                      ICECELL(L) = .TRUE.
                    endif

                    ! *** RESTORE THE TEMPERATURE (HEAT)
                    TEM(L,KC) = TEM(L,KC)*RAVAIL + (TEM1(L,KC)-TF)*(1.-RAVAIL)
                  endif
                endif  ! *** End of ISICE > 2
              endif
            endif
          endif
        endif
      enddo
    enddo
    !$OMP END DO

  else    ! *** if( ISGWIE >= 1 )then

    ! *** ADJUST SOURCES AND SINKS ESTIMATING SURFACE AND GROUNDWATER
    ! *** AVAILABLE FOR EVAPOTRANSPIRATON AND INFILTRATION
    !$OMP DO PRIVATE(ND,LF,LL,L,K,DTAGW,RIFTRL,RAVAIL,QSUMIET,QEAVAIL,QSUMTMP,DIFQVOL)
    do ND = 1,NDM
      LF = 2+(ND-1)*LDM
      LL = min(LF+LDM-1,LA)
      do L = LF,LL
        EVAPSW(L) = 0.
        EVAPGW(L) = 0.

        if( HP(L) > HDRY9 )then
          ! *** APPLY MAXIMUM ET
          EVAPSW(L) = EVAPT(L)*DXYP(L)      ! *** m3/s
          QGW(L) = 0.

          ! *** CALCULATE DEPTH OF ACTIVE GROUNDWATER ELEV BELOW SURFACE
          DTAGW = BELV(L)-AGWELV(L)
          if( DTAGW > 0.0 )then

            ! *** INFLITRATION CAN OCCUR, CALCULATE LIMITING RATE TO BRING
            ! *** GW ELEV TO SOIL SURFACE
            RIFTRL = RNPOR*DTAGW*DELTI

            ! *** SET RIFTRL TO MIN OF LIMITING RATE OR ACTUAL RATE
            RIFTRL = min(RIFTRM,RIFTRL)

            ! *** ESTIMATE RATE BASED ON AVAILABLE SURFACE WATER
            RAVAIL = (H1P(L)-HDRY)*DELTI-EVAPT(L)

            ! *** SET RIFTRL TO MIN OF AVAILABLE RATE OR LIMITING RATE
            RIFTRL = min(RAVAIL,RIFTRL)

            ! *** CONVERT TO VOLUME FLOW UNITS
            QGW(L) = -RIFTRL*DXYP(L)           ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT)
          endif

          ! *** ADJUST VOLUME OUTFLOWS OF WET CELLS
          if( QSUME(L) < 0.0 .and. HP(L) < HWET )then
            QSUMIET = EVAPSW(L)-QGW(L)
            QEAVAIL = DXYP(L)*(H1P(L)-HDRY)*DELTI-QSUMIET
            QEAVAIL = max(QEAVAIL,0.0)
            QEAVAIL = -QEAVAIL
            QSUMTMP = max(QSUME(L),QEAVAIL)

            ! *** UPDATE LAYER FLOWS
            DIFQVOL = QSUME(L)-QSUMTMP
            do K = KSZ(L),KC
              QSUM(L,K) = QSUM(L,K) - DIFQVOL*DZC(L,K)
            enddo
            QSUME(L) = QSUMTMP
          endif
        else
          ! *** CELL IS 'DRY'
          QGW(L)  = 0.
          EVAPSW(L) = 0.
          QSUME(L) = max(QSUME(L),0.0)      ! *** ONLY ALLOW INFLOWS
          do K = 1,KC
            QSUM(L,K) = max(QSUM(L,K),0.0)
          enddo
        endif
      enddo
    enddo
    !$OMP END DO
  endif
  !$OMP END PARALLEL

  ! ***  WRITE DIAGNOSTIC FILE FOR VOLUME SOURCES,SINKS, ETC
  ITMPD = 0
  if( ISDIQ == 2 .and. ISTL == 2 ) ITMPD = 1
  if( ISDIQ == 1 ) ITMPD = 1
  NTT = 3 + NDYM + NTOX + NSED + NSND
  
6665 FORMAT(' US ELEVATION CONTROL STRUCTURE TABLE OUT OF BOUNDS ')
6666 FORMAT(' SINGLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS ')
6667 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
6668 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
6676 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, UP ')
6677 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
6678 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
6679 FORMAT(' HUF,HUL,HDF,HDL = ',4(2X,E12.4))
6686 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, DW ')
6687 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
6688 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
6111 FORMAT(' INVALID NQSIJ LOCATION, NQSIJ,I,J = ',3I5)
  return

  END SUBROUTINE

