! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE WAVEBL 
  ! *** WAVEBL subroutine is used to specify information for wave boundary layer only        
  ! *** NTSWV  = Number of time steps for gradual introduction of wave forcing
  ! *** WVLCAL = 1: To calculate wave length / 0: wave length is from wave.inp
  ! *** ISDZBR = Write diagnostics for effective wave current boundary layer roughness
  ! *** IFWAVE = 0 To use wave.inp
  ! ***          1 To directly input SWAN formatted files 
  ! *** SWANGRP= 1 To use SWAN output for whole grid / 0 SWAN output for locations (x,y)
      
  ! CHANGE RECORD 
  ! DATE MODIFIED     BY               DESCRIPTION        
  !-- ------------------------------------------------------------------
  ! 2017-05-23        Dang Chung       IMPROVED WAVE CALCULATION
  ! 2012-08-04        Dang Chung       IMPROVED WAVE INPUT & CALCULATION 
  ! 2011-08-08        PAUL M. CRAIG &  RESTRUCTURED TO F90 AND CORRECTED CODE
  !                   Dang Chung      
  

  use GLOBAL
  use GETSWANMOD     ! *** DIRECTLY READ SWAN OUTPUT
  use WAVELENGTH
  use Variables_MPI

  implicit none
  integer :: NWV,IWV,JWV,NW,L,K
  integer :: IOS,ND,LF,LL
  
  real(RKD) :: WDEP1,RRLS,WPRD

  real   :: WWVH,WANGLE,WWPRDP,WVLEN,DISPTMP
  real   :: AEXTMP,UWORBIT,REYWAVE
  real   :: RA,FCW,CDTMP,VISMUDD
  real, external :: CSEDVIS

  
  if( JSWAVE == 0 )then
    ! *** INITIALIZE AND INPUT WAVE INFORMATION       
    do L = 1,LC
      HMPW(L) = 0.
      HMUW(L) = 0.
      HMVW(L) = 0.
      WVKHC(L) = 0.
      WVKHU(L) = 0.
      WVKHV(L) = 0.
      WVTMP1(L) = 0.
      WVTMP2(L) = 0.
      WVTMP3(L) = 0.
      WVTMP4(L) = 0.
      WVENEP(L) = 0.
      UWVSQ(L) = 0.
      QQWC(L) = 1.E-12
      QQWCR(L) = 1.E-12
      QQWV1(L) = 1.E-12   ! *** BOUNDARY LAYER STRESS DUE TO WAVES
      QQWV2(L) = 1.E-12   ! *** WATER COLUMN STRESS DUE TO WAVES
      QQWV3(L) = 1.E-12
    enddo
    do K = 1,KC
      do L = 1,LC
        WVHUU(L,K) = 0.
        WVHVV(L,K) = 0.
        WVHUV(L,K) = 0.
        WVPP(L,K) = 0.
        WVPU(L,K) = 0.
        WVPV(L,K) = 0.
        FXWAVE(L,K) = 0.
        FYWAVE(L,K) = 0.
      enddo
    enddo

    if( ISSTEAD == 0 )then
      ! *** UNSTEADY WAVE
      if( process_id == master_id ) write(*,'(A)')'WAVE: READING WAVETIME.INP'     
      call GETWAVEDAY
      if( WAVEDAY(1) > TIMEDAY .or. WAVEDAY(NWVTIM) < TIMEDAY) STOP 'TIMEDAY IS OUT OF THE WAVE-TIME RANGE!'
      do NW = 1,NWVTIM-1  
        if( WAVEDAY(NW+1) > TIMEDAY .and. WAVEDAY(NW) <= TIMEDAY )then
          IWVCOUNT = NW
          EXIT
        endif
      enddo
    else
      ! *** STEADY WAVE
      IWVCOUNT = 1
      NWVTIM   = 1
      allocate(WAVEDAY(2))
      WAVEDAY(1:2) = TIMEDAY
    endif
    
    ! *** READ FIRST WAVE PERIOD FIELD
    if( IFWAVE == 0 )then
      if( process_id == master_id )then 
        write(*,'(A)')'READING WAVE.INP'
        write(*,'(A36,F12.4)')'WAVE: READING WAVE.INP AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
      endif
      call GETWAVEINP
      
    elseif( IFWAVE == 1 )then
      if( process_id == master_id ) write(*,'(A38,F12.4)')'WAVE: READING SWAN OUPUT AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
      if( SWANGRP == 1 )then
        call GETSWAN_GRP
      else
        call GETSWAN_LOC
      endif
    endif    
    JSWAVE = 1
    
    elseif( JSWAVE == 1 .and. IWVCOUNT < NWVTIM .and. TIMEDAY >= WAVEDAY(IWVCOUNT+1) )then
      ! *** UPDATE WAVE parameterS
      IWVCOUNT = IWVCOUNT + 1   
      if( IFWAVE == 0 )then
        if( process_id == master_id ) write(*,'(A36,F12.4)')'WAVE: READING WAVE.INP AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
        call GETWAVEINP

      elseif( IFWAVE == 1 )then
        if( process_id == master_id ) write(*,'(A38,F12.4)')'WAVE: READING SWAN OUPUT AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
        if( SWANGRP == 1 )then
          call GETSWAN_GRP
        else
          call GETSWAN_LOC
        endif
      endif
    endif
    
  ! ****************************************************************************
  ! *** FINISHED READING DATA, NOW SETUP THE WAVE TABLE AND COMPUTE WAVE parameterS
  ! ***  GENERATE WAVE TABLE
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(ND,L,LF,LL,WDEP1,WPRD,RRLS)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = min(LF+LDM-1,LA)
        
    do L = LF,LL
      WV(L).HEIGHT = min(0.75*HP(L),WV(L).HEISIG)                     ! *** INCLUDING BREAKING WAVE
      if( WV(L).HEIGHT >= WHMI .and. HP(L) > HDRYWAV )then
        WDEP1 = HP(L)
        WPRD  = 2.*PI/WV(L).FREQ
        if( WVLCAL == 1 )then
          call BISEC(DISRELATION,WLMIN,WLMAX,EPS,WPRD,WDEP1,0._8,0._8,RRLS)
          WV(L).LENGTH = RRLS
        endif
        WV(L).K   = max( 2.*PI/WV(L).LENGTH, 0.01 )                   ! *** ANGULAR WAVE NUMBER (RAD/M)
        WV(L).KHP = min(WV(L).K*HP(L),SHLIM)
      endif
    enddo
  enddo
  !$OMP END DO
  
  ! ***  HDRYWAV IS USED TO LIMIT WAVE HEIGHT INSTEAD OF HP(L) = 0.55
  ! ***  THE WAVE TURBULENT INTENSITY, QQWV
  ! ***  AND SQUARED HORIZONTAL WAVE OBRITAL VELOCITY MAGNITUDE

  if( process_id == master_id )then
      ! *** WRITE DIAGNOSTICS
      !$OMP SINGLE
      if( ISDZBR > 0 )then
        open(1,FILE = OUTDIR//'WAVEBL_DIA.OUT')
        close(1,STATUS = 'DELETE')
        open(1,FILE = OUTDIR//'WAVEBL_DIA.OUT')
      endif
      !$OMP END SINGLE
   endif

  !$OMP DO PRIVATE(ND,L,LF,LL,AEXTMP,UWORBIT,VISMUDD,REYWAVE,RA,FCW,CDTMP)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = min(LF+LDM-1,LA)
    do L= LF,LL  !2,LA
      ! *** SET ZBRE AS N IKURADSE ROUGHNESS
      if( ISTRAN(7) > 0 )then
        ! *** BASE ON NIKURADSE ROUGHNESS (APPROXIMATE 2.5*D50)
        ZBRE(L) = max(SEDDIA50(L,KBT(L)),1E-6)*2.5
      else
        ZBRE(L) = KSW    !Z0 = KSW/30
      endif

      if( WV(L).HEIGHT >= WHMI .and. HP(L) > HDRYWAV )then
        AEXTMP   = 0.5*WV(L).HEIGHT/SINH(WV(L).KHP)
        UWORBIT  = AEXTMP*WV(L).FREQ
        UWVSQ(L) = UWORBIT*UWORBIT
        if( UWORBIT < 1.E-6 )then
          UWVSQ(L) = 0.    
          QQWV1(L) = 0.
          CYCLE
        endif
        VISMUDD = 1.36E-6      !DHC: 2010-05-06
        if( ISMUD >= 1 ) VISMUDD = CSEDVIS(SED(L,KSZ(L),1))
        REYWAVE = UWORBIT*AEXTMP/VISMUDD
        RA= AEXTMP/ZBRE(L)       !Relative roughness = A/Ksw

        ! *** COMPUTE FRICTION FACTOR DUE TO WAVE
        if( REYWAVE <= 5D5 )then
          !LAMINAR
          FCW  = 2*REYWAVE**(-0.5)

        elseif( REYWAVE>5D5 .and. RA>1.57 )then
          !** TURBULENT SMOOTH WAVE BOUNDARY LAYER
          FCW = 0.09*REYWAVE**(-0.2)

        elseif( REYWAVE>5D5 .and. RA <= 1.57 )then
          !** TURBULENT ROUGH WAVE BOUNDARY LAYER
          FCW = EXP(5.2*RA**(-0.19)-6)    ! *** Baird's paper
          FCW = min(FCW,0.3)
        endif
        CDTMP = 0.5*FCW
        QQWV1(L) = min(CDTMP*UWORBIT*UWORBIT,QQMAX)

      else
        QQWV1(L) = 0.
        UWVSQ(L) = 0.
      endif
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  if( process_id == master_id )then
      if( ISDZBR > 0 )then
        write(1,'(A)') '**  L    I    J    WVWHA(L)   WVFRQL(L)    WVLEN(L)    WVKHP(L)    UWVSQ(L)    QQWV1(L)     ZBRE(L)'
        do L = 2,LA
            write(1,'(3I5,7f12.5)') L,IL(L),JL(L),WV(L).HEIGHT ,WV(L).FREQ ,WV(L).LENGTH ,WV(L).K ,UWVSQ(L) ,QQWV1(L) ,ZBRE(L)
        enddo
        close(1)
      endif
  endif
  
END SUBROUTINE
