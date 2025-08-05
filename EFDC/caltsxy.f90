! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALTSXY(INITFLAG)

  ! ***SUBROUTINE CALTSXY UPDATES TIME VARIABLE SURFACE WIND STRESS  
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------!
  !    2015-01       Paul M. Craig     Added fully coupled Ice Sub-model
  !                  Dang H Chung        
  !    2014-12       Paul M. Craig     Added the new Evaporation Approach
  !    2014-10       D H CHUNG         ADD ICE THICKNESS CALCULATION USING CE-QUAL-W2
  !    2011-03       PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP

  use GLOBAL
  use Allocate_Initialize
  use HEAT_MODULE, only:EQUILIBRIUM_TEMPERATURE
  use FIELDS
  use CYCLONE
  use Variables_MPI
  use MPI
  use Communicate_Ghost_Routines
  use CALCSERMOD,only:LIN_INTER_COEF
  
  implicit none
  
  integer, intent(IN) :: INITFLAG
  
  integer :: ND, K, LF, LL, LP, L, LS, LN, I, J, NA, NS, NW, M1, M2, NI, ICECOVL, IWMP, ITMP, IMAP, IICE, IERR
  integer, save :: NUMSTEPS
  
  real    :: TIME, TDIFF, WTM1, WTM2, DEGM1, DEGM2, WINDS1, WINDS2, WINDE1, WINDE2, WINDN1, WINDN2, TAUICE !, CONVRT2
  real    :: WNDFAC, C2, TSEAST, TSNORT, WINDXX, WINDYY, TMPVAL, SVPWET, CLEVAPTMP, CCNHTTTMP, TMPVL1, U10, CD10
  real    :: VAPORP, WSPD, ET, CSHE
  
  character(1) :: VE1, VC1
  character(6) :: WS4(20), AS4(20)
  
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS                 ! MODEL TIMING TEMPORARY VARIABLE
  real(RKD), save     :: DAYNEXT              ! *** Whole day counter
  real, save          :: SUMT                 ! *** Temperature sum

  real,save,allocatable,dimension(:) :: CLOUDTT  
  real,save,allocatable,dimension(:) :: EVAPTT  
  real,save,allocatable,dimension(:) :: PATMTT  
  real,save,allocatable,dimension(:) :: RAINTT  
  real,save,allocatable,dimension(:) :: RHATT  
  real,save,allocatable,dimension(:) :: SOLSWRTT  
  real,save,allocatable,dimension(:) :: SVPATT  
  real,save,allocatable,dimension(:) :: TATMTT  
  real,save,allocatable,dimension(:) :: TDEWTT  
  real,save,allocatable,dimension(:) :: TWETTT  
  real,save,allocatable,dimension(:) :: VPATT  
  real,save,allocatable,dimension(:) :: WINDE  
  real,save,allocatable,dimension(:) :: WINDN  
  real,save,allocatable,dimension(:) :: RAINSUM
  real,save,allocatable,dimension(:) :: TLAYER  
  
  if( .not. allocated(CLOUDTT) )then
    call AllocateDSI(CLOUDTT,  NASERM,  0.0)
    call AllocateDSI(EVAPTT,   NASERM,  0.0)
    call AllocateDSI(RAINTT,   NASERM,  0.0)
    call AllocateDSI(PATMTT,   NASERM,  0.0)
    call AllocateDSI(RHATT,    NASERM,  0.0)
    call AllocateDSI(SOLSWRTT, NASERM,  0.0)
    call AllocateDSI(SVPATT,   NASERM,  0.0)
    call AllocateDSI(TATMTT,   NASERM,  0.0)
    call AllocateDSI(TDEWTT,   NASERM,  0.0)
    call AllocateDSI(TWETTT,   NASERM,  0.0)
    call AllocateDSI(VPATT,    NASERM,  0.0)
    
    call AllocateDSI(WINDE,    NWSERM,  0.0)
    call AllocateDSI(WINDN,    NWSERM,  0.0)
    
    call AllocateDSI(WINDSXX,  LCM,     0.0)
    call AllocateDSI(WINDSXY,  LCM,     0.0)
    call AllocateDSI(WINDSYX,  LCM,     0.0)
    call AllocateDSI(WINDSYY,  LCM,     0.0)
    call AllocateDSI(RAINSUM,  LCM,     0.0)
    call AllocateDSI(TLAYER,   KCM,     0.0)

    if( ISVHEAT > 0 )then
      do L = 2,LA  
        CLEVAP(L) = SVREVC(L)       ! *** Already converted by 0.001 when read
        CCNHTT(L) = SVRCHC(L)       ! *** Already converted by 0.001 when read 
      enddo 
    else
      do L = 2,LA  
        CLEVAP(L) = 0.001*ABS(REVC)
        CCNHTT(L) = 0.001*ABS(RCHC)  
      enddo 
    endif
    
    ! *** Setup daily temperature average
    NUMSTEPS = 0
    SUMT = 0.0
    DAYNEXT = DBLE(INT(TIMEDAY)) 
    
    if( NASER == 0 ) TATMT(2:LA) = TEM(2:LA,KC)  ! *** ASSIGN INITIAL WATER TEMPERATURE TO DRY BULB
    
    ! *** WRITE TO LOG  FILE
    if( ISTRAN(2) > 0 .and. ISVHEAT > 0 )then
      open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
      write(mpi_efdc_out_unit,'(//A)') 'SURFACE HEAT parameterS'
      FORALL(I = 1:NWSER) WS4(I) = 'WSER'//CHAR(48+I)
      FORALL(I = 1:NASER) AS4(I) = 'ASER'//CHAR(48+I)
      write(mpi_efdc_out_unit,'(A5,3A10,"  | ",20A8)') 'L','REVC','RCHC','BACKKE', WS4(1:NWSER), AS4(1:NASER)
      do L = 2,LA
        if( LSVHTWINDE(L) )then
          VE1 = '+'
        else
          VE1 = ' '
        endif
        if( LSVHTWINDC(L) )then
          VC1 = '+'
        else
          VC1 = ' '
        endif
        if( NWSER > 1 .and. NASER > 1 )then
          write(mpi_efdc_out_unit,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), WNDWHT(1:NWSER,L,1), ATMWHT(1:NASER,L,1)
        elseif( NWSER > 1 )then
          write(mpi_efdc_out_unit,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), WNDWHT(1:NWSER,L,1), 1.
        elseif( NASER > 1 )then
          write(mpi_efdc_out_unit,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), 1., ATMWHT(1:NASER,L,1)
        else
          write(mpi_efdc_out_unit,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), 1., 1.
        endif
      enddo
    endif
    
    if( ISTRAN(2) == 0 )then
      ! *** SET CONSTANT TEMPERATURE
      do L = 2,LA
        do K = 1,KC
          TEM(L,K) = TEMO
        enddo
      enddo
      
      ! *** DELME - PROVIDE USER INPUTS
      do L = 2,LA
        TEMB(L) = TEMO
        TATMT(L) = TEMO + 5.0
        PATMT(L) = 1000.0
      enddo
      
    endif
    close(mpi_efdc_out_unit)
  endif

  ! *** *******************************************************************!
  ! Initialize wind sheltered surface gas transfer
  if( NITER == 1 .and. (NWSER > 0 .or. WIND.IFLAG > 0) )then
    if( DEBUG )then
      open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
      write(mpi_efdc_out_unit,'(a)'), ' *** Writing WINDSHELT.OUT'
    endif
    
    ! *** Compute wind sheltering adjusted wind curvature coefficients
    do L = 2,LA  
      I = IL(L)  
      J = JL(L)  
      if( WINDSTKA(L) > 0.0 )then  
        ! *** IF WINDSTKA > 0 BOTH X AND Y COMPONENTS ARE APPLIED  
        WINDSXX(L) = CVN(L)  
        WINDSXY(L) = -CVE(L)  
        WINDSYX(L) = -CUN(L)  
        WINDSYY(L) = CUE(L)  
      else  
        ! *** IF WINDSTKA < 0 SLECTIVELY APPLY X AND Y COMPONENTS  
        ! *** FIRST CASE IS FULLY OPEN WATER  
        WINDSXX(L) = CVN(L)  
        WINDSXY(L) = -CVE(L)  
        WINDSYX(L) = -CUN(L)  
        WINDSYY(L) = CUE(L)  
        LS = LSC(L)  
        LN = LNC(L)  
        ! *** SECOND CASE IS 1D CHANNEL IN COMP X DIRECTION  
        if( SVB(L) < 0.5 .and. IJCT(I,J-1) /= 5 )then  
          if( SVB(LN) < 0.5 .and. IJCT(I,J+1) /= 5 )then  
            WINDSXX(L) = CVN(L)  
            WINDSXY(L) = -CVE(L)  
            WINDSYX(L) = -1000.  
            WINDSYY(L) = 0.  
          endif  
        endif  
        ! *** THIRD CASE IS 1D CHANNEL IN COMP Y DIRECTION  
        if( SUB(L) < 0.5 .and. IJCT(I-1,J) /= 5 )then  
          if( SUB(LEC(L)) < 0.5 .and. IJCT(I+1,J) /= 5 )then  
            WINDSXX(L) = 0.  
            WINDSXY(L) = -1000.  
            WINDSYX(L) = -CUN(L)  
            WINDSYY(L) = CUE(L)  
          endif  
        endif  
      endif
      
      if( DEBUG )then
        open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        write(mpi_efdc_out_unit,1111) IL(L), JL(L), WINDSTKA(L), WINDSXX(L), WINDSXY(L), WINDSYX(L), WINDSYY(L)  
        close(mpi_efdc_out_unit)  ! *** EFDC+ log file File
      endif
      
    enddo  
1111 FORMAT(2I5,10F10.6)  
    
  endif  

  LDAYLIGHT = .FALSE.
  LCHECKICE = .FALSE.
  LRAIN     = .FALSE.
  
  ! *** *******************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** *******************************************************************
  ! *** INTERPOLATE THE WINDS FOR EACH SERIES TO THE CURRENT TIME
  if( NWSER > 0 .or. WIND.IFLAG > 0 )then  !2018-10-12, NTL: UPDATE TIME VARIABLE WIND FIELD
    ! *** UPDATE THE FORCING WIND DATA TO THE CURRENT TIME
    !$OMP SINGLE
    do NW = 1,NWSER  
      TIME = TIMESEC/TSWND(NW).TMULT
      M2 = MTSWLAST(NW)
      do while (TIME > TSWND(NW).TIM(M2))
        M2 = M2 + 1
        if( M2 > TSWND(NW).NREC )then
          M2 = TSWND(NW).NREC   !** THIS ALLOWS USING EXTRAPOLATION !!!
          exit
        endif
      enddo
      MTSWLAST(NW) = M2
      M1 = M2 - 1
            
      TDIFF     = TSWND(NW).TIM(M2)-TSWND(NW).TIM(M1)  
      WTM1      = (TSWND(NW).TIM(M2)-TIME)/TDIFF  
      WTM2      = (TIME-TSWND(NW).TIM(M1))/TDIFF  
      DEGM1     = 90.-TSWND(NW).VAL(M1,2)
      DEGM2     = 90.-TSWND(NW).VAL(M2,2)  
      WINDS1    = WTM1*TSWND(NW).VAL(M1,1)+WTM2*TSWND(NW).VAL(M2,1)  
      WINDS2    = WTM1*TSWND(NW).VAL(M1,1)+WTM2*TSWND(NW).VAL(M2,1)  
      WINDE1    = TSWND(NW).VAL(M1,1)*COS(DEGM1/57.29578)  
      WINDN1    = TSWND(NW).VAL(M1,1)*SIN(DEGM1/57.29578)  
      WINDE2    = TSWND(NW).VAL(M2,1)*COS(DEGM2/57.29578)  
      WINDN2    = TSWND(NW).VAL(M2,1)*SIN(DEGM2/57.29578)  
      WINDE(NW) = WTM1*WINDE1+WTM2*WINDE2  
      WINDN(NW) = WTM1*WINDN1+WTM2*WINDN2
      
      ! *** EE7.3 CONVERT TO 2 METERS FOR ALL CALCULATIONS (0.003 IS OPEN GRASSLAND Z0)
      CONVRT2   = LOG(2.0/0.003)/LOG(WINDH(NW)/0.003)
      WINDE(NW) = CONVRT2*WINDE(NW)
      WINDN(NW) = CONVRT2*WINDN(NW)
        
    enddo  
    
    ! --------------------------------------------------------------------------
    ! *** CALCULATE THE WIND STRESS
    ! *** HANDLE TIME VARIABLE WIND MAPS
    if( NWSER > 1 )then      
      IWMP = 1
      if( NWNDMAP > 1 )then
        do ITMP = 1,NWNDMAP
          if( TIMEDAY >= TWNDMAPBEG(ITMP) .and. TIMEDAY < TWNDMAPEND(ITMP) )then
            IWMP = ITMP
            exit
          endif
        enddo
      endif
    endif

    ! *** CONVERSION FACTOR FROM 2 METERS TO 10 METERS (0.003 IS OPEN GRASSLAND Z0)
    CONVRT2  = LOG(10.0/0.003)/LOG(2./0.003)
    
    !2018-10-12, NTL: UPDATE TIME VARIABLE WIND FIELD
    if(WIND.IFLAG > 0 )then
      call UPDATEFIELD(WIND,TIMEDAY,1,WNDVELE)
      call UPDATEFIELD(WIND,TIMEDAY,2,WNDVELN)
    endif
    !$OMP END SINGLE
    
    !$OMP DO PRIVATE(ND,LF,LL,L,WNDFAC,C2,TSEAST,TSNORT,WSPD,WINDXX,WINDYY,U10,CD10)
    do ND = 1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = min(LF+LDM-1,LA)

      if( WIND.IFLAG == 0 )then !2018-10-12, NTL: UPDATE TIME VARIABLE WIND FIELD
        ! *** Convert 2m wind vectors to 10m wind vectors for each cell
        if( NWSER > 1 )then
          ! *** MULTIPLE STATIONS
          do L = LF,LL 
            WNDVELE(L) = SUM(WNDWHT(1:NWSER,L,IWMP)*WINDE(1:NWSER))*CONVRT2
            WNDVELN(L) = SUM(WNDWHT(1:NWSER,L,IWMP)*WINDN(1:NWSER))*CONVRT2
          enddo              
        elseif( NWSER == 1 )then
          ! *** SINGLE STATION
          do L = LF,LL
            WNDVELE(L) = WINDE(1)*CONVRT2
            WNDVELN(L) = WINDN(1)*CONVRT2
          enddo 
        endif
      endif
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO
    
    if( INITFLAG == 0 .or. ISRESTI == 0 .or. Restart_In_Ver < 1210 )then
      !$OMP DO PRIVATE(ND,LF,LL,L)
      do ND = 1,NDM  
        LF = 2+(ND-1)*LDM  
        LL = min(LF+LDM-1,LA)
      
        do L = LF,LL
          call WINDSTRESS(L)
        enddo  
      
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    endif
    
  endif    ! *** END OF NWSER>0

  ! *************************************************************************
  ! *** UPDATE ASER DATA AND COMPUTE SURFACE EXCHANGE parameterS
  if( NASER > 0 )then  
    ! *** UPDATE ATMOSPHERIC DATA TO THE CURRENT TIME
    !$OMP SINGLE
    do NA = 1,NASER  
      TIME = TIMESEC/TSATM(NA).TMULT 
      M2 = MTSALAST(NA)   
      do while (TIME > TSATM(NA).TIM(M2))
        M2 = M2+1
        if( M2 > TSATM(NA).NREC )then
          M2 = TSATM(NA).NREC
          exit
        endif
      enddo
      MTSALAST(NA) = M2 
      M1 = M2-1
                  
      TDIFF = TSATM(NA).TIM(M2)-TSATM(NA).TIM(M1)  
      WTM1 = (TSATM(NA).TIM(M2)-TIME)/TDIFF  
      WTM2 = (TIME-TSATM(NA).TIM(M1))/TDIFF  
      PATMTT(NA) = WTM1*TSATM(NA).VAL(M1,1)+WTM2*TSATM(NA).VAL(M2,1)  
      TATMTT(NA) = WTM1*TSATM(NA).VAL(M1,2)+WTM2*TSATM(NA).VAL(M2,2)  
      TWETTT(NA) = WTM1*TSATM(NA).VAL(M1,3)+WTM2*TSATM(NA).VAL(M2,3)  
      RAINTT(NA) = WTM1*TSATM(NA).VAL(M1,4)+WTM2*TSATM(NA).VAL(M2,4)
      if( RAINTT(NA) > 0.0 ) LRAIN = .TRUE.
      EVAPTT(NA) = WTM1*TSATM(NA).VAL(M1,5)+WTM2*TSATM(NA).VAL(M2,5)  

      SOLSWRTT(NA) = WTM1*TSATM(NA).VAL(M1,6)+WTM2*TSATM(NA).VAL(M2,6)  
      SOLSWRTT(NA) = SOLSWRTT(NA)*(1.-REFL)     ! *** ADJUST INCIDENT LIGHT FOR REFLECTANCE
      if( SOLSWRTT(NA) > 0.2 ) LDAYLIGHT = .TRUE.
    
      CLOUDTT(NA) = WTM1*TSATM(NA).VAL(M1,7)+WTM2*TSATM(NA).VAL(M2,7)  
      SVPATT(NA) = 10.**((0.7859+0.03477*TATMTT(NA))/(1.+0.00412*TATMTT(NA)))  
      if( IRELH(NA) == 0 )then  
        ! *** COMPUTE RHA FROM WET BULB
        TMPVAL = 0.00066*(1.0 + 0.00115*TWETTT(NA))
        SVPWET = 10.**((0.7859 + 0.03477*TWETTT(NA))/(1. + 0.00412*TWETTT(NA)))  
        TMPVL1 = SVPWET - TMPVAL*PATMTT(NA)*(TATMTT(NA) - TWETTT(NA))
        RHATT(NA) = max(TMPVL1/SVPATT(NA),0.01)
      else
        ! *** DIRECT ENTRY OF RELATIVE HUMIDITY
        RHATT(NA) = TWETTT(NA)  
      endif 

      ! *** PREVENT LOG OF ZERO
      RHATT(NA) = max(RHATT(NA),0.0001) 
      
      ! *** Compute Dew Point Temp in C.  From Jensen et al. (1990) ASCE Manual No. 70
      ! *** Ambient vapor pressure in kPa 
      VAPORP = RHATT(NA)* 0.611*EXP(17.27*TATMTT(NA)/(TATMTT(NA)+237.3)) 
      TDEWTT(NA) = (116.9+237.3*LOG(VAPORP))/(16.78-LOG(VAPORP))
      
      ! *** ACTUAL VAPOR PRESSURE (MB)
      VPATT(NA) = RHATT(NA)*SVPATT(NA)  
      
      ! *** Accumulate average air temperature
      NUMSTEPS = NUMSTEPS + 1
      SUMT = SUMT + TATMTT(NA)
    enddo
    
    ! *** Update average air temperature
    if( TIMEDAY >= DAYNEXT .and. NUMSTEPS > 0 )then
      TATMT(1) = SUMT/FLOAT(NUMSTEPS)
      TATMT(1) = max(TATMT(1), 0.0)     ! *** TATMT(1) is only used for QC bed temperatures
      NUMSTEPS = 0
      SUMT = 0.0
      DAYNEXT = DAYNEXT + 1.
    endif
    
    ! *** HANDLE TIME VARIABLE ATMOSPHERIC MAPS
    if( NASER > 1 )then
      IMAP = 1
      if( NATMMAP > 1 )then
        do ITMP = 1,NATMMAP
          if( TIMEDAY >= TATMMAPBEG(ITMP) .and. TIMEDAY < TATMMAPEND(ITMP) )then
            IMAP = ITMP
            exit
          endif
        enddo
      endif
      
    endif
    !$OMP END SINGLE
    
    ! *** UPDATE CELL ATMOSPHERIC parameterS FOR ALL CELLS NOT JUST WET
    !$OMP DO PRIVATE(ND,LF,LL,L)
    do ND = 1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = min(LF+LDM-1,LA)
      RAINSUM(ND) = 0.0
      
      ! *************************************************************************
      ! *** ADJUST parameterS BASED ON ATM WEIGHTING, IF NEEDED
      if( NASER > 1 )then  
        do L = LF,LL
          PATMT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*PATMTT(1:NASER))
          TATMT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*TATMTT(1:NASER)) 
          RAINT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*RAINTT(1:NASER))  
          EVAPT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*EVAPTT(1:NASER))  
          SOLSWRT(L) = SUM(ATMWHT(1:NASER,L,IMAP)*SOLSWRTT(1:NASER)) 
          CLOUDT(L)  = SUM(ATMWHT(1:NASER,L,IMAP)*CLOUDTT(1:NASER))  
          SVPAT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*SVPATT(1:NASER)) 
          RHAT(L)    = SUM(ATMWHT(1:NASER,L,IMAP)*RHATT(1:NASER))  
          TDEWT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*TDEWTT(1:NASER))  
          VPAT(L)    = SUM(ATMWHT(1:NASER,L,IMAP)*VPATT(1:NASER))  
          RAINSUM(ND) = RAINSUM(ND) + RAINT(L)
        enddo  
        
        if( IATMP > 0 )then
          do L = LF,LL
            ! *** CONVERT ATM PRESSURE (MB) TO METERS OF WATER * G FOR UNITS OF M2/S2 
            ATMP(L) = PATMT(L)*0.0101974*G 
          enddo  
        endif
        
      elseif( NASER == 1 )then
        do L = LF,LL
          PATMT(L)   = PATMTT(1)  
          TATMT(L)   = TATMTT(1)  
          RAINT(L)   = RAINTT(1)  
          EVAPT(L)   = EVAPTT(1)  
          SOLSWRT(L) = SOLSWRTT(1)  
          CLOUDT(L)  = CLOUDTT(1)  
          SVPAT(L)   = SVPATT(1)  
          RHAT(L)    = RHATT(1)  
          TDEWT(L)   = TDEWTT(1)  
          VPAT(L)    = VPATT(1)    
          RAINSUM(ND) = RAINSUM(ND) + RAINT(L)
        enddo
      endif
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO

    ! *** UPDATE TIME VARIABLE PRESSURE FIELD
    !$OMP SINGLE
    RAINTSUM1 = RAINTSUM0
    RAINTSUM0 = SUM(RAINSUM(:))
    
    if( PRESS.IFLAG > 0 )then
      call UPDATEFIELD(PRESS,TIMEDAY,1,PATMT)
      do L = 2,LA
        ! *** CONVERT ATM PRESSURE (MB) TO METERS OF WATER * G FOR UNITS OF M2/S2 
        ATMP(L) = PATMT(L)*0.0101974*G 
      enddo
    endif
    
    ! *** UPDATE TIME VARIABLE RAINFALL FIELD
    if( RAIN.IFLAG > 0 )then
      call UPDATEFIELD(RAIN,TIMEDAY,1,RAINT)
      RAINTSUM0 = SUM(RAINT(:))
    endif
    
    ! *** UPDATE TIME VARIABLE EVAPORATION FIELD
    if( EVAP.IFLAG > 0 )then
      call UPDATEFIELD(EVAP,TIMEDAY,1,EVAPT)
    endif
    !$OMP END SINGLE
    
    ! CYCLONE PRESSURE AND WIND FIELDS OVERLAID ON TOP OF AMBIENT ONES
    if( ICYCLONE > 0 )then
        call CycloneFields(TIMEDAY)
        if( INITFLAG == 0 .or. ISRESTI == 0 .or. Restart_In_Ver < 1210 )then
          !$OMP DO PRIVATE(ND,LF,LL,L)    !,WNDFAC,C2,TSEAST,TSNORT,WSPD,WINDXX,WINDYY,U10,CD10)
          do ND = 1,NDM  
            LF = 2+(ND-1)*LDM  
            LL = min(LF+LDM-1,LA)
        
            do L = LF,LL
                call WINDSTRESS(L)
            enddo  
          enddo  
          !$OMP END DO
        endif
    endif
    
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,NA,CLEVAPTMP,CCNHTTTMP)
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)
        
      ! *** EVAPORATION WIND FUNCTION (FW)
      if( IEVAP == 2 .or. ISTOPT(2) == 1 )then
        if( ISVHEAT > 0 )then
          do LP = LF,LL
            L = LWET(LP)  
            if( LSVHTWINDE(L) )then
              CLEVAPTMP = 1.E-3*(0.8 + 0.065*WINDST(L))  
              CLEVAP(L) = max(SVREVC(L),CLEVAPTMP)  
            endif
          enddo
        else
          if( REVC < 0. )then  
            CLEVAPTMP = 0.001*ABS(REVC)
            do LP = LF,LL
              L = LWET(LP)  
              CLEVAP(L) = 1.E-3*(0.8 + 0.065*WINDST(L))  
              CLEVAP(L) = max(CLEVAP(L),CLEVAPTMP)  
            enddo  
          endif
        endif
        
      endif

      ! *** CONDUCTIVE HEAT TRANSFER WIND FUNCTION
      if( ISTOPT(2) == 1 )then
        if( ISVHEAT > 0 )then
          do LP = LF,LL
            L = LWET(LP)  
            if( LSVHTWINDC(L) )then
              CCNHTTTMP = 1.E-3*(0.8+0.065*WINDST(L))  
              CCNHTT(L) = max(SVRCHC(L),CCNHTTTMP)  
            endif
          enddo
        else
          if( RCHC < 0. )then  
            CCNHTTTMP = 0.001*ABS(RCHC)
            do LP = LF,LL
              L = LWET(LP)  
              CCNHTT(L) = 1.E-3*(0.8+0.065*WINDST(L))  
              CCNHTT(L) = max(CCNHTT(L),CCNHTTTMP)  
            enddo  
          endif
        endif
      endif
      
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO
    !$OMP SINGLE
  
    ! *** CHECK FOR POSSIBLE ICE CONDITIONS
    if( ISICE > 2 .and. NWSER > 0 .and. NASER > 0 .and. (IS2TIM > 0 .or. ISTL == 3) )then
      do NA = 1,NASER
        ! *** COMPUTE EQUILIBRIUM TEMPERATURE
        WSPD = SQRT(WINDE(1)*WINDE(1) + WINDN(1)*WINDN(1))
        call EQUILIBRIUM_TEMPERATURE(WSPD, TDEWTT(NA), TATMTT(NA), SOLSWRTT(NA), ET, CSHE) 
        if( ET < 1.0 ) LCHECKICE = .TRUE.
      enddo
    endif
    !$OMP END SINGLE
  
  elseif( EVAP.IFLAG > 0 .or. RAIN.IFLAG > 0 )then

    ! *** UPDATE TIME VARIABLE RAINFALL FIELD
    !$OMP SINGLE
    if( RAIN.IFLAG > 0 )then
      call UPDATEFIELD(RAIN,TIMEDAY,1,RAINT)
      RAINTSUM0 = SUM(RAINT(:))
    endif
    
    ! *** UPDATE TIME VARIABLE EVAPORATION FIELD
    if( EVAP.IFLAG > 0 )then
      call UPDATEFIELD(EVAP,TIMEDAY,1,EVAPT)
    endif
    !$OMP END SINGLE
    
  endif    ! *** END OF NASER>0 SECTION

  ! *** ASSIGN TEMPERATURES, IF NOT COMPUTED
  if( ISTRAN(2) == 0 )then
    if( IVOLTEMP > 1 )then
      ! *** ASSIGN TIME VARYING TEMPERATURES
      !$OMP SINGLE
      NS = IVOLTEMP - 1
      TIME = TIMEDAY
      call LIN_INTER_COEF(TIME, TSTEM(NS).TIM, NS, 2, M1, M2, WTM1, WTM2)
      do K = 1,KC
        TLAYER(K) = WTM1*TSTEM(NS).VAL(M1,K) + WTM2*TSTEM(NS).VAL(M2,K)
      enddo 
      !$OMP END SINGLE
      
      !$OMP DO PRIVATE(ND,LP,L,K) 
      do ND = 1,NDM  
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TEM(L,K) = TLAYER(K)
          enddo
        enddo
      enddo
      !$OMP END DO
    endif
  endif
  
  ! *** *******************************************************************C
  if( ISICE > 0 )then
    if( (ISICE == 1 .or. ISICE == 2) .and. NISER > 0 )then
      !$OMP SINGLE
      do NI = 1,NISER    
        TIME = TIMESEC/TSICE(NI).TMULT         ! *** THE SAME TIME UNIT IN ISER.INP
        M2 = MITLAST(NI)
        do while( TIME > TSICE(NI).TIM(M2) )
          M2 = M2 + 1
          if( M2 > TSICE(NI).NREC )then 
            M2 = TSICE(NI).NREC
            exit
          endif
        enddo
        MITLAST(NI) = M2           
        M1 = M2 - 1
        
        if( ISICE == 1 )then
          TDIFF = TSICE(NI).TIM(M2) - TSICE(NI).TIM(M1)
          WTM1  = (TSICE(NI).TIM(M2)-TIME)/TDIFF
          WTM2  = (TIME-TSICE(NI).TIM(M1))/TDIFF
          RICECOVT(NI) = WTM1*TSICE(NI).VAL(M1,1) + WTM2*TSICE(NI).VAL(M2,1)
          
        elseif( ISICE == 2 )then
          RICECOVT(NI) = TSICE(NI).VAL(M1,1)
        endif
        
      enddo
      !$OMP END SINGLE
                   
      if( NISER > 1 )then
        ! *** ONLY FOR ISICE = 1 & NISER>1
        ! *** HANDLE SERIES WEIGHTING ICE MAPS
        !$OMP SINGLE
        IICE = 1
        if( NICEMAP > 1 )then
          do ITMP = 1,NICEMAP
            if( TIMEDAY >= TICEMAPBEG(ITMP) .and. TIMEDAY < TICEMAPEND(ITMP) )then
              IICE = ITMP
              exit
            endif
          enddo
        endif      
        !$OMP END SINGLE

        !$OMP DO PRIVATE(ND,LF,LL,LP,L,CLEVAPTMP,CCNHTTTMP)
        do ND = 1,NDM  
          LF = (ND-1)*LDMWET+1  
          LL = min(LF+LDMWET-1,LAWET)
        
          do LP = LF,LL
            L = LWET(LP)  
            ICECOVER(L) = SUM(RICEWHT(IICE,L,1:NISER)*RICECOVT(1:NISER))
            ICETHICK(L) = ICECOVER(L)      ! *** IF ISICE == 1 or 2 ICETHICK is not used.  
          enddo
        enddo   ! *** END OF DOMAIN LOOP
        !$OMP END DO
        
      elseif( NISER == 1 )then
        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM  
          LF = (ND-1)*LDMWET+1  
          LL = min(LF+LDMWET-1,LAWET)
        
          do LP = LF,LL
            L = LWET(LP)  
            ICECOVER(L) = RICECOVT(1) 
            ICETHICK(L) = ICECOVER(L)      ! *** IF ISICE == 1 or 2 ICETHICK is not used.
          enddo
        enddo   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      endif

      ! *** Set ICECELL flag
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,ICECOVL)
      do ND = 1,NDM  
        LF = (ND-1)*LDMWET+1  
        LL = min(LF+LDMWET-1,LAWET)
        
        do LP = LF,LL
          L = LWET(LP)  
          ICECOVL = NINT(ICECOVER(L))
          ICECOVER(L) = min(FLOAT(ICECOVL),1.0)
          if( ICECOVL == 0 )then
            ICETHICK(L) = 0.0
            ICECELL(L) = .FALSE.
          else
            ICECELL(L) = .TRUE.
          endif
        enddo
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    
    endif     ! *** ISICE = 1 OR ISICE = 2

    ! *** *******************************************************************
    if( INITFLAG == 0 .or. ISRESTI == 0 .or. Restart_In_Ver < 1210 )then
      !$OMP DO PRIVATE(ND,LF,LL,L,TAUICE)
      do ND = 1,NDM  
        LF = 2+(ND-1)*LDM  
        LL = min(LF+LDM-1,LA)

        do L = LF,LL
          if( ICECELL(L) )then
            ! *** SHEAR ON TOP WATER LAYER DUE TO DRAG ON BOTTOM OF ICE (CDICE = 0.001)
            TAUICE = -CDICE*SQRT( U(L,KC)*U(L,KC) + V(L,KC)*V(L,KC) )
            TAUICE = TAUICE*(1. + 9.17*ICETHICK(L)/(HP(L) + ICETHICK(L)) )
            if( TAUICE < -0.075 )then
              TAUICE = -0.075
              if( HP(L) < 0.5 ) TAUICE = HP(L)*TAUICE
            endif
            TSX(L) = TAUICE*U(L,KC)
            TSY(L) = TAUICE*V(L,KC)
          endif
        enddo
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    endif
  endif     ! *** END OF ISICE>0
  
  ! *** *******************************************************************C
  !$OMP END PARALLEL
    
  ! *** GOTM testing
  if( iGOTM_Test > 0 )then
    do L = 2,LA
      TSX(L) = SurfaceShearX/1000.
      TSY(L) = SurfaceShearY/1000.
    enddo
  endif

  ! ****************************************************************************
  if( NWSER > 1 .or. iGOTM_Test )then
    call MPI_barrier(DSIcomm, ierr)
    TTDS = DSTIME(0)
    call Communicate_1D2(TSX, TSY)
    TMPITMP = DSTIME(0) - TTDS
  endif
  ! ****************************************************************************

  return 
  
END  

SUBROUTINE WINDSTRESS(L)
  use GLOBAL
  use Variables_MPI
  use MPI

  implicit none
  real    :: WSPD, WNDFAC, TSEAST, TSNORT, C2, U10, CD10, WINDXX, WINDYY, CONVRT2
  integer, intent(IN):: L

  ! *** CONVERSION FACTOR FROM 2 METERS TO 10 METERS (0.003 IS OPEN GRASSLAND Z0)
  CONVRT2  = LOG(10.0/0.003)/LOG(2./0.003)
  ! *** CASE 0 MAGNITUDE SHELTERING AND NO DIRECTIONAL SHELTERING
  if( WINDSTKA(L) > 0.0 )then
    WNDFAC = WINDSTKA(L)
    WNDVELE(L) = WNDFAC*WNDVELE(L)
    WNDVELN(L) = WNDFAC*WNDVELN(L)
    WSPD = SQRT( WNDVELE(L)*WNDVELE(L) + WNDVELN(L)*WNDVELN(L) )    ! *** Compute magnitude of 10m wind components
    U10 = max(WSPD,0.25)
    
    if( IWDRAG < 2 )then
      ! *** ORIGINAL EFDC WIND DRAG
      if( U10 < 5. )then
        C2 = 1./U10
        CD10 = 3.83111E-005*(C2**3) - 0.000308715*(C2**2) + 0.00116012*C2 + 0.000899602
      elseif( U10 >= 5.0 .and. U10 <= 7. )then
        CD10 = -5.37642E-006*(U10**3) + 0.000112556*(U10**2) - 0.000721203*U10 + 0.00259657
      elseif( U10 >= 7. )then
        CD10 = -3.99677E-007*(U10**2) + 7.32937E-005*U10 + 0.000726716
      endif
    elseif( IWDRAG == 2 )then
      ! *** Use HERSBACH 2011, EUROPEAN CENTRE FOR MEDIUM-RANGE WEATHER FORECASTS (ECMWF)
      CD10 = (0.00103 + 0.00004 * U10**1.48) / U10**0.21
    elseif( IWDRAG == 3 )then
      ! *** Use COARE3.6 (EDSON, ET.AL. 2013) SIMPLIFICATION AT NEUTRAL BOUYANCY
      if( U10 > 20. )then
        CD10 = 8.394E-05*U10 + 8.285E-04
      else
        CD10 = 1.376E-08*U10**4 - 9.841E-07*U10**3 + 2.324E-05*U10**2 - 1.012E-04*U10 + 8.971E-04
      endif
    elseif( IWDRAG == 4 )then
      ! *** USER DEFINED WIND DRAG
      if( U10 <= WDRAG1 )then
        CD10 = CDRAG1
      elseif( U10 >= WDRAG2 )then
        CD10 = CDRAG2
      else
        CD10 = CDRAG1 + (CDRAG2 - CDRAG1)/(WDRAG2 - WDRAG1)*(U10 - WDRAG1)
      endif
    endif
    
    if( IWDRAG > 0 )then
      ! *** COMPUTE WIND SPEED RELATIVE TO WATER
      TSEAST = WINDSXX(L)*WNDVELE(L) + WINDSXY(L)*WNDVELN(L) - U(L,KC)
      TSNORT = WINDSYX(L)*WNDVELE(L) + WINDSYY(L)*WNDVELN(L) - V(L,KC)
      ! *** 1.225E-3 = RHOa/RHOw (DIMENSIONLESS)
      TSX(L) = 1.225E-3*CD10*U10*TSEAST             ! *** TSX IS THE WIND SHEAR IN THE U DIRECTION (M2/S2)
      TSY(L) = 1.225E-3*CD10*U10*TSNORT             ! *** TSY IS THE WIND SHEAR IN THE V DIRECTION (M2/S2)
     if( IWDRAG == 3 .and. ISTOPT(2) == 2 )then
      TSX(L) = 1.225E-3*CDCOARE(L)*U10*TSEAST 
      TSY(L) = 1.225E-3*CDCOARE(L)*U10*TSNORT
     endif
    else
      TSEAST = 1.225E-3*CD10*U10*WNDVELE(L)
      TSNORT = 1.225E-3*CD10*U10*WNDVELN(L)
      TSX(L) = WINDSXX(L)*TSEAST + WINDSXY(L)*TSNORT
      TSY(L) = WINDSYX(L)*TSEAST + WINDSYY(L)*TSNORT
    endif
    WINDCD10(L) = CD10
    U10 = WSPD                                      ! *** Reset to 10 meters
    
  elseif( WINDSTKA(L) < 0.0 )then
    ! *** CASE 1 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, OPEN WATER
    if( WINDSYX(L) > -99.0 .and. WINDSXY(L) > -99.0 )then
      WNDFAC = ABS(WINDSTKA(L))
      WNDVELE(L) = WNDFAC*WNDVELE(L)
      WNDVELN(L) = WNDFAC*WNDVELN(L)
      U10 = SQRT( WNDVELE(L)*WNDVELE(L) + WNDVELN(L)*WNDVELN(L) )    ! *** Compute magnitude of 10m wind components

      C2 = 1.2E-6*(0.8 + 0.065*U10)
      TSEAST = C2*U10*WNDVELE(L)
      TSNORT = C2*U10*WNDVELN(L)
      TSX(L) = WINDSXX(L)*TSEAST+WINDSXY(L)*TSNORT
      TSY(L) = WINDSYX(L)*TSEAST+WINDSYY(L)*TSNORT
    endif

    ! *** CASE 2 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, X CHANNEL
    if( WINDSYX(L) < -99.0 )then
      WINDXX = WINDSXX(L)*WNDVELE(L)+WINDSXY(L)*WNDVELN(L)
      WNDFAC = ABS(WINDSTKA(L))
      WINDXX = WNDFAC*WNDVELE(L)
      U10 = ABS(WINDXX)
      TSX(L) = 1.2E-6*(0.8 + 0.065*U10)*U10*WINDXX
      TSY(L) = 0.
    endif

    ! *** CASE 3 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, Y CHANNEL
    if( WINDSXY(L) < -99.0 )then
      WINDYY = WINDSYX(L)*WNDVELE(L)+WINDSYY(L)*WNDVELN(L)
      WNDFAC = ABS(WINDSTKA(L))
      WINDYY = WNDFAC*WINDYY
      U10 = ABS(WINDYY)
      TSX(L) = 0
      TSY(L) = 1.2E-6*(0.8 + 0.065*U10)*U10*WINDYY
    endif
  endif
  
  ! *** WINDST is wind magnitude at 2 meters (heat exchange)
  WINDST(L) = U10/CONVRT2

END SUBROUTINE

