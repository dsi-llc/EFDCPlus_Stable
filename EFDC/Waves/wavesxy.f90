! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE WAVESXY 
  ! *** WAVESXY subroutine links and processes externally generated wave data for
  ! *** * First Line: Initial Data Settings (1 line)
  ! *** NWVCELLS  - Number Of Cells With Wave Data                                                                      
  ! *** WVPRD   - Wave Period (sec)                                                                                   
  ! *** ISWCBL  - 1 Activates Wave-Current BL Model                                                                   
  ! *** ISWRSR  - 1 Activates Inclusion Of Rotational Component Of Rad Stress                                        
  ! *** ISWRSI  - 1 Activates Inclusion Of Irrotational Component Of Rad Stress                                      
  ! *** NTSWV   - Number Of Time Steps For Gradual Introduction Of Wave Forcing                                      
  ! *** WVDISV  - Fraction Of Wave Dissipation As Source In Vertical TKE Closure                                     
  ! *** WVDISH  - Fraction Of Wave Dissipation As Source In Horiz Smagorinky's Subgrid Closure  (NOT USED)                       
  ! *** WVLSH   - Weight For Depth As The Horiz SSG Eddy Viscosity Length Scale,       ISHMDF>0                                      
  ! *** WVLSX   - Weight For Sqrt(Dxdy) As The Horiz SSG Eddy Viscosity Length Scale,  ISHMDF>0
  ! *** ISWVSD  - 1 Include Nondiverg Wave Stokes Drift In Mass Transport                                             
  ! *** ISDZBR  - 1 Write Diagnostics For Effect Wave Current Bndry Layer Roughness                                  
  ! ***
  ! *** * Second Line: (repeated NWVCELLS times)                                                                        
  ! *** IWV     - I Of The Wave Cell                                                                                  
  ! *** JWV     - J Of The Wave Cell                                                                                  
  ! *** ENETMP  - Wave Energy, 0.5*G*Abs(Amplitude)**2  (M3/S2)                                                       
  ! *** SXXTMP  - Radiation Stresses- XX Component (kg/s^2)                                                           
  ! *** SYYTMP  - Radiation Stresses- YY Component (kg/s^2)                                                           
  ! *** SXYTMP  - Radiation Stresses- XY Component (kg/s^2)                                                           
  ! *** WVDISP  - Energy Dissipation  (M3/S3)                        [WV.DISSIPA]                                                                       
  ! *** WANGLE  - Wave Angle Measured from East positive CCW (deg)   
                                                                                   
  ! CHANGE RECORD                                                                                                     
  ! DATE MODIFIED     BY               DESCRIPTION                                                                    
  !-- -------------------------------------------------------------------!   
  ! 2017-05-23        Dang Chung       IMPROVED WAVE CALCULATION
  ! 2011-04-26        Dang Chung       Corrected Wave Formulation & Orbital Velocity                                 
  ! 2011-03-02        PAUL M. CRAIG    RESTRUCTURED AND CORRECTED CODE     
  ! 2012-08-04        Dang Chung       IMPROVED WAVE INPUT & CALCULATION 

  use GLOBAL
  use GETSWANMOD
  use WAVELENGTH
  
  implicit none
  
  real :: UWORBIT, AEXTMP
  real :: WWVH,  WANGLE,  WVLEN, DISPTMP
  real :: TAUTMP, CORZBR, CDRGTMP
  real :: WVWHAUT, WVWHAVT, DZITMP
  real(RKD) :: SNHBOTU, SNHTOPV, SNHBOTV, TMPPU, TMPPV, SNHTOPU
  real(RKD) :: WG, WN, WK, WE, WDEP1, RRLS, WPRD
  real(RKD) :: SXXTMP, SXYTMP, SYYTMP, RATIO3, RKHM1, RKHM2
  real(RKD) :: ZTOP, ZBOT, SINHTOP, SINHBOT, SINH3, COSH3
  real(RKD) :: COSHTOP, COSHBOT, TMPVAL3, TMPP1, TMPP2, TMP3
  real, external :: CSEDVIS

 
  integer :: NW,L,K,IOS,NWV,IWV,JWV,LS,LW,LSW,LE
  integer :: LN,LNW,LSE,ND,LF,LL
  
  ! ***  INPUT WAVE INFORMATION                                                                                            
  if( JSWAVE == 0 )then
    RSWRSR = FLOAT(ISWRSR)
    RSWRSI = FLOAT(ISWRSI)
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
      QQWV1(L) = 1.E-12  ! *** BED TURBULENT INTENSITY DUE TO WAVES ONLY
      QQWV2(L) = 1.E-12  ! *** WATER COLUMN TURBULENT INTENSITY DUE TO WAVES
      QQWV3(L) = 1.E-12  ! *** WATER COLUMN TURBULENT INTENSITY DUE TO WAVES MODIFIED FOR NON-COHESIVE MOVING BED
    enddo
    do K = 1,KC
      do L = 1,LC
        WVHUU(L,K)  = 0.
        WVHVV(L,K)  = 0.
        WVHUV(L,K)  = 0.
        WVPP(L,K)   = 0.
        WVPU(L,K)   = 0.
        WVPV(L,K)   = 0.
        FXWAVE(L,K) = 0.
        FYWAVE(L,K) = 0.
      enddo
    enddo
    ! ***  INITIALIZE VERTICAL DISTRIBUTION OF WAVE DISSIPATION AS SOURCE                                                    
    ! ***  TO VERTICAL TKE CLOSURE                                                                                           
    if( KC == 2 )then
      WVDTKEM(1) = WVDISV
      WVDTKEP(1) = WVDISV
    endif
    if( KC == 3 )then
      WVDTKEM(1) = WVDISV
      WVDTKEP(1) = 0.5*WVDISV
      WVDTKEM(2) = 0.5*WVDISV
      WVDTKEP(2) = WVDISV
    endif
    if( KC >= 4 )then
      WVDTKEM(1) = WVDISV
      WVDTKEP(1) = 0.5*WVDISV
      WVDTKEM(KS) = 0.5*WVDISV
      WVDTKEP(KS) = WVDISV
      do K = 2,KS-1
        WVDTKEM(K) = 0.5*WVDISV  ! *** BOTTOM
        WVDTKEP(K) = 0.5*WVDISV  ! *** TOP
      enddo
    endif
    
    if( ISSTEAD == 0 )then
      ! *** UNSTEADY WAVE
      if( process_id == master_id ) write(*,'(A)')'WAVE: READING WAVETIME.INP'     
      call GETWAVEDAY
      if( WAVEDAY(1) > TIMEDAY .or. WAVEDAY(NWVTIM) < TIMEDAY ) CALL STOPP('TIMEDAY IS OUT OF THE WAVE-TIME RANGE!')
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
       
    ! *** INPUT THE WAVE FILE: FIRST READ 
    if( IFWAVE == 0 )then
      ! *** Read wave parameters from WAVES.INP
      ! ***   WWVH    - Wave height (m)
      ! ***   WANGLE  - Wave Angle Measured from East positive CCW (deg)
      ! ***   WVPRD   - Wave period (sec)
      ! ***   WVLEN   - Wave length (m)
      ! ***   DISPTMP - Wave energy dissipation  (M3/S3)
      if( process_id == master_id )then 
        write(*,'(A)')'WAVE: READING WAVE.INP'
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
    call HDEP0_PLSWAVE
    
  elseif( JSWAVE == 1 .and. IWVCOUNT<NWVTIM .and. TIMEDAY >= WAVEDAY(IWVCOUNT+1) )then
    ! *** WAVE UPDATE
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
    call HDEP0_PLSWAVE
  endif
  ! *** END OF WAVE UPDATE *************************************************
  
  ! ***  DISTRIBUTE WVHUU, WVHVV, AND WV.DISSIPA OVER DEPTH            
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(ND,L,LF,LL,WDEP1,WPRD,RRLS,WK,WE,WG,WN,SXXTMP,SXYTMP,SYYTMP)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)
    
    do L = LF,LL
      WV(L).HEIGHT = MIN(0.75*HP(L),WV(L).HEISIG)            ! *** INCLUDING BREAKING WAVE
      if( WV(L).HEIGHT >= WHMI .and. HP(L) > HDRYWAV )then
        WDEP1 = HP(L)
        WPRD  = 2.*PI/WV(L).FREQ
        if( WVLCAL == 1 )then
          call BISEC(DISRELATION,WLMIN,WLMAX,EPS,WPRD,WDEP1,0._8,0._8,RRLS)
          WV(L).LENGTH = RRLS
        endif
        WV(L).K   = MAX( 2.*PI/WV(L).LENGTH, 0.01 ) 
        WV(L).KHP = MIN(WV(L).K*HP(L),SHLIM)
        WK = 2.0*WV(L).KHP                             ! *** 4*PI/WV(L).LENGTH*HP(L) = 2KH
        WE = 9.81*1000.*WV(L).HEIGHT**2 /16.           ! *** TOTAL WAVE ENERGY : KG/S2  
        WVENEP(L) = WE/1000.                           ! *** WAVE ENERGY/RHO: M3/S2
        WG = WK/SINH(WK)                               ! *** 2KH/SINH(2KH)
        WN = (WG+1)/2                               
        SXXTMP = WE*(WG/2+WN*COS(WV(L).DIR)**2)        ! *** RADIATION SHEAR STRESS [Kg/S2] 
        SXYTMP = WE*WN/2*SIN(2*WV(L).DIR)              ! *** RADIATION SHEAR STRESS [Kg/S2] 
        SYYTMP = WE*(WG/2+WN*SIN(WV(L).DIR)**2)        ! *** RADIATION SHEAR STRESS [Kg/S2] 
        WVHUU(L,KC) = SXXTMP/1000.
        WVHVV(L,KC) = SYYTMP/1000.
        WVHUV(L,KC) = SXYTMP/1000.
      else
        WVHUU(L,KC) = 0.
        WVHVV(L,KC) = 0.
        WVHUV(L,KC) = 0.    
      endif
    enddo
  enddo
  !$OMP END DO
  
  ! *** DISTRIBUTE VALUES ACROSS KC    
  !$OMP DO PRIVATE(ND,L,LF,LL,K,RKHM1,RKHM2,SINH3,COSH3,RATIO3,ZTOP,ZBOT) &
  !$OMP    PRIVATE(SINHTOP,SINHBOT,COSHTOP,COSHBOT,TMPVAL3,TMPP1,TMP3,TMPP2)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)
    do L = LF,LL 
      if( WV(L).HEIGHT >= WHMI .and. HP(L) > HDRYWAV )then
        RKHM1 = WV(L).KHP
        RKHM2 = 2.*RKHM1
        SINH3 = SINH(RKHM2)
        COSH3 = COSH(RKHM1)
        RATIO3 = 0.5+(RKHM1/SINH3)
        do K = 1,KC
          ZTOP = Z(L,K)
          ZBOT = Z(L,K-1)
          SINHTOP = SINH(RKHM2*ZTOP)
          SINHBOT = SINH(RKHM2*ZBOT)
          COSHTOP = COSH(RKHM1*ZTOP)
          COSHBOT = COSH(RKHM1*ZBOT)
          TMPVAL3 = (RKHM2*(ZTOP-ZBOT)+SINHTOP-SINHBOT)/(RKHM2+SINH3)

          ! *** APPLY FACTOR
          WVHUU(L,K) = TMPVAL3*WVHUU(L,KC)
          WVHVV(L,K) = TMPVAL3*WVHVV(L,KC)
          WVHUV(L,K) = TMPVAL3*WVHUV(L,KC)
          WV(L).DISSIPA(K) = TMPVAL3*WV(L).DISSIPA(KC)  ! FROM GETSWANMOD
          
          TMPP1 = -0.5*(ZTOP-ZBOT)+(ZTOP*COSHTOP-ZBOT*COSHBOT)/COSH3
          TMP3 = SINH3-2.
          if( ABS(TMP3)<1D-3) TMP3 = SIGN(1D-3,TMP3)
          TMPP2 = (RATIO3-1.)*(SINHTOP-SINHBOT-2.*(ZTOP-ZBOT))/TMP3

          ! *** LIMIT RANGE WHEN WV.K~0.72
          if( ABS(TMPP1)>0.5 )then
            TMPP1 = SIGN(0.5,TMPP1)
          endif
          if( ABS(TMPP2)>0.5 )then
            TMPP2 = SIGN(0.5,TMPP2)
          endif
          WVPP(L,K) = WVENEP(L)*(TMPP1+TMPP2)
        enddo
      else
        WVHUU(L,1:KC) = 1D-6
        WVHVV(L,1:KC) = 1D-6
        WVHUV(L,1:KC) = 1D-6
        WV(L).DISSIPA(1:KC) = 1D-6
        WVPP(L,1:KC)  = 1D-6
      endif
    enddo
  enddo
  !$OMP END DO
  
  ! ***  INITIALIZE WAVE-CURRENT BOUNDARY LAYER MODEL CALCULATING                                                          
  ! ***  THE WAVE TURBULENT INTENSITY, QQWV                                                                                
  ! ***  AND SQUARED HORIZONTAL WAVE ORBITAL VELOCITY MAGNITUDE            
  !$OMP DO PRIVATE(ND,L,LF,LL,AEXTMP,UWORBIT) &
  !$OMP    PRIVATE(TMPVAL3,TAUTMP,CORZBR,CDRGTMP)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)
    
    do L = LF,LL
      ! *** SET ZBRE AS GRAIN/SKIN ROUGHNESS (M)
      if( ISTRAN(7) > 0 )then
        ! *** BASE ON NIKURADSE ROUGHNESS (APPROXIMATE 2.5*D50)
        !ZBRE(L) = SEDDIA50(L,KBT(L)) !*2.5/30.
        ZBRE(L) = MAX(SEDDIA50(L,KBT(L)),1E-6)*2.5
      else
        ZBRE(L) = KSW
      endif

      if( WV(L).HEIGHT >= WHMI .and. HP(L) > HDRYWAV )then

        ! *** WAVE HEIGHT / HYPERBOLIC SINE OF WAVE NUMBER
        AEXTMP = 0.5*WV(L).HEIGHT/SINH(WV(L).KHP)

        ! *** SQUARED HORIZONTAL WAVE ORBITAL VELOCITY MAGNITUDE
        UWORBIT  = AEXTMP*WV(L).FREQ  !WVFRQ
        UWVSQ(L) = UWORBIT*UWORBIT

        if( UWORBIT >= 1.E-6 )then
 
          ! *** CHECK FOR NON-COHESIVE TRANSPORT
          if( ISTRAN(7) > 0 )then
            TMPVAL3 = UWORBIT*SQRT( AEXTMP/(30.*ZBRE(L)) )  
            TAUTMP = TMPVAL3/TAUR(NSED+1)
            CORZBR = 1.+1.2*TAUTMP/(1.+0.2*TAUTMP)
            ZBRE(L) = ZBRE(L)*CORZBR
          endif
          CDRGTMP = (30.*ZBRE(L)/AEXTMP)**0.2
          CDRGTMP = 5.57*CDRGTMP-6.13
          CDRGTMP = EXP(CDRGTMP)
          CDRGTMP = MIN(CDRGTMP,0.22)
          TMPVAL3 = 0.5*CDRGTMP*UWVSQ(L)
          QQWV2(L) = MIN(CTURB2*TMPVAL3,QQMAX)
        else
          QQWV2(L) = QQLMIN
        endif
      else
        UWVSQ(L)   = 0.0
        WV(L).FREQ = 1.0
        QQWV2(L) = QQLMIN
      endif
    enddo
  enddo
  !$OMP END DO
  
  ! ***  COMPUTE CELL FACE QUANTITIES WVPU,WVPV         
  !$OMP DO PRIVATE(ND,L,LF,LL)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)
    do L = LF,LL
      if( WV(L).HEIGHT >= WHMI ) WVKHU(L)= MIN(HMUW(L)*WV(L).K,SHLIM)     !VALKH(HFFDG)
      if( WV(L).HEIGHT >= WHMI ) WVKHV(L)= MIN(HMVW(L)*WV(L).K,SHLIM)     !VALKH(HFFDG)
    enddo
  enddo
  !$OMP END DO
  
  !$OMP DO PRIVATE(ND,L,LF,LL,LS,TMPVAL3,WVWHAUT,WVWHAVT)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)
    do L = LF,LL
      LS = LSC(L)
      TMPVAL3   = 0.5*WV(L).FREQ*WV(L).FREQ     !WVFRQ*WVFRQ
      WVTMP1(L) = MAX(SINH(WVKHU(L)),1D-6)
      WVWHAUT   = (WV(L).HEIGHT + SUB(L)*WV(LWC(L)).HEIGHT)/(1.+SUB(L))
      WVTMP2(L) = TMPVAL3*WVWHAUT*WVWHAUT /(WVTMP1(L)*WVTMP1(L))
      WVWHAVT   = (WV(L).HEIGHT + SVB(L)*WV(LSC(L)).HEIGHT)/(1.+SVB(L))
      WVTMP3(L) = MAX(SINH(WVKHV(L)),1D-6)
      WVTMP4(L) = TMPVAL3*WVWHAVT*WVWHAVT /(WVTMP3(L)*WVTMP3(L))
    enddo
  enddo
  !$OMP END DO
  
  !$OMP DO PRIVATE(ND,LF,LL,K,L,ZTOP,ZBOT,SNHTOPU,SNHBOTU,SNHTOPV,SNHBOTV,TMPPU,TMPPV)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)
    do K = 1,KC
      do L = LF,LL
        ZTOP = Z(L,K)
        ZBOT = Z(L,K-1)
        SNHTOPU = SINH(WVKHU(L)*ZTOP)
        SNHBOTU = SINH(WVKHU(L)*ZBOT)
        SNHTOPV = SINH(WVKHV(L)*ZTOP)
        SNHBOTV = SINH(WVKHV(L)*ZBOT)
        TMPPU = (1.-ZTOP)*SNHTOPU*(ZTOP*WVTMP1(L)-SNHTOPU) - (1.-ZBOT)*SNHBOTU*(ZBOT*WVTMP1(L)-SNHBOTU)
        TMPPV = (1.-ZTOP)*SNHTOPV*(ZTOP*WVTMP3(L)-SNHTOPV) - (1.-ZBOT)*SNHBOTV*(ZBOT*WVTMP3(L)-SNHBOTV)
        WVPU(L,K) = WVTMP2(L)*TMPPU
        WVPV(L,K) = WVTMP4(L)*TMPPV
      enddo
    enddo
  enddo
  !$OMP END DO

  ! ***  CALCULATE THE INTERNAL MODE NET X AND Y WAVE REYNOLDS STRESS FORCINGS                                             
  !$OMP DO PRIVATE(ND,LF,LL,K,L,LS,LN,LW,LE,LNW,LSE,DZITMP)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)

    do K = 1,KC
      do L = LF,LL
        if( LKSZ(L,K) ) CYCLE
        DZITMP = 1.*DZIC(L,K)
        LS = LSC(L)
        LN = LNC(L)
        LW = LWC(L)
        LE = LEC(L)
        LNW = LNWC(L)
        LSE = LSEC(L)
        FXWAVE(L,K) = DZITMP*SUB(L)*SPB(L)*( RSWRSI*(DYU(L)*(WVPP(L,K) - WVPP(LW,K)) + DYU(L)*WVPU(L,K)*(HMPW(L)-HMPW(LW)))    &
                                            +RSWRSR*(DYP(L)*WVHUU(L,K) - DYP(LW)*WVHUU(LW,K)                                   &
                                                    + 0.5*(DXV(LN)+DXV(LNW))*WVHUV(LN,K) - 0.5*(DXV(L)+DXV(LW))*WVHUV(L,K)) )
      
        FYWAVE(L,K) = DZITMP*SVB(L)*SPB(L)*( RSWRSI*(DXV(L)*(WVPP(L,K) - WVPP(LS,K)) + DXV(L)*WVPV(L,K)*(HMPW(L)-HMPW(LS )))   &
                                            +RSWRSR*(DXP(L)*WVHVV(L,K) - DXP(LS)*WVHVV(LS,K)                                   &
                                                    + 0.5*(DYU(LE)+DYU(LSE))*WVHUV(LE,K) - 0.5*(DYU(L)+DYU(LS))*WVHUV(L,K)) )
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

contains

! *** 
SUBROUTINE HDEP0_PLSWAVE
  ! *** COMPUTE THE DERIVED WATER DEPTHS AT U/V CELL FACES
  do L = 2,LA
    LS = LSC(L)
    LW = LWC(L)
    LSW = LSWC(L)
    HMPW(L) = HP(L) + WV(L).HEIGHT
    
    ! HEIGHT U AND V FACE
    HMUW(L) = 0.5*(DXYP(L)*HMPW(L)+DXYP(LW)*HMPW(LW))/(DXU(L)*DYU(L))
    HMVW(L) = 0.5*(DXYP(L)*HMPW(L)+DXYP(LS)*HMPW(LS))/(DXV(L)*DYV(L))
  enddo
 END SUBROUTINE

END SUBROUTINE

SUBROUTINE GETWAVEINP
  
  use GLOBAL
  use INFOMOD
  use Variables_MPI
  use Variables_MPI_Write_Out
  use Broadcast_Routines
  
  ! *** READ WAVE PARAMS FROM WAVE.INP
  integer :: L,IOS,NWV,IWV,JWV, LG
  real   :: WVCX,WVCY,WVDX,WVDY,DISPTMP
  real   :: WWVH,WANGLE,WVLEN
  real,allocatable,dimension(:) :: WVDISSIPA
  allocate(WVDISSIPA(LCM_Global))
  
  if( .not. allocated(LWVCELL_Global) )then
    allocate(LWVCELL_Global(LCM_Global))
    allocate(LWVMASK_Global(LCM_Global))
    LWVMASK_Global = .false.
  endif
  
  if( process_id == master_id )then
    write(*,'(A)')'WAVE: READING WAVE.INP'
    open(WUNI,FILE = 'wave.inp',STATUS = 'UNKNOWN')
    call SKIPCOM(WUNI,'*')
  
    ! *** READ BUFFER DATA 
    do NW = 1,IWVCOUNT-1
      do NWV = 2,LA_Global
        read(WUNI,*,IOSTAT = IOS) IWV, JWV, WWVH, WANGLE, WVPRD, WVLEN, DISPTMP
        if( IOS > 0 )then
          write(*,'(A45,I5,A12,I5)') '***  READ ERROR ON FILE WAVE.INP AT CELL L = ',NWV,'BLOCK NO = ',NW
          call STOPP('.')
        endif
      enddo
    enddo
  
    NWVCELLS = 0
    do NWV = 2,LA_Global
      read(WUNI,*,IOSTAT = IOS)IWV,JWV,WWVH,WANGLE,WVPRD,WVLEN,DISPTMP
      if( IOS > 0 )then
        write(*,'(A45,I5)') '***  READ ERROR ON FILE WAVE.INP AT CELL L = ',NWV
        call STOPP('.')
      endif
      LG = LIJ_GLOBAL(IWV,JWV)
      if( LG > 1 .and. LG <= LA_Global )then
        NWVCELLS = NWVCELLS + 1
        LWVCELL_Global(NWVCELLS) = LG   ! *** PMC - 2015-12 - This was pointing to L-1 but LWVCELL is used as an L array reference.  This has to be L not L-1.
        LWVMASK_Global(LG) = .TRUE.
        
        ! *** WAVE HEIGHT (M )then ZERO FOR VEGETATION EFFECT
        WV_Global(LG).HEISIG = WWVH 
        if( ISVEG > 0 )then
          if( MVEG_Global(LG) /= MVEGOW  )then
            WV_Global(LG).HEIGHT = 0.
          else
            if( (MVEG_Global(LWC_Global(LG)) /= MVEGOW) .or. &
                (MVEG_Global(LEC_Global(LG)) /= MVEGOW) .or. &
                (MVEG_Global(LSC_Global(LG)) /= MVEGOW) .or. &
                (MVEG_Global(LNC_Global(LG)) /= MVEGOW)      )  WV_Global(LG).HEIGHT = 0.5*WV_Global(LG).HEIGHT    ! delme - this needs to be updated
          endif
        endif
        
        if( WV_Global(LG).HEISIG >= WHMI .and. WVPRD > 0. )then      
          WVDX= COS(WANGLE*PI/180)
          WVDY= SIN(WANGLE*PI/180)       
          WVCX =  CVN_Global(LG)*WVDX - CVE_Global(LG)*WVDY       
          WVCY = -CUN_Global(LG)*WVDX + CUE_Global(LG)*WVDY
          WV_Global(LG).DIR    = ATAN2(WVCY,WVCX)        
          WV_Global(LG).LENGTH = WVLEN     
          WV_Global(LG).FREQ   = 2.*PI/WVPRD
          WVDISSIPA(LG) = DISPTMP                        ! *** ENERGY DISSIPATION DUE TO BREAKING WAVE (M3/S3)
          WV_Global(LG).K = MAX( 2.*PI/WV_Global(LG).LENGTH, 0.01 )          ! *** ANGULAR WAVE NUMBER (RAD/M)
        else
          WV_Global(LG).DIR    = 0.
          WV_Global(LG).LENGTH = 0.
          WV_Global(LG).FREQ   = 1.
          WV_Global(LG).DISSIPA(KC) = 0.
          WV_Global(LG).HEISIG = 0.
          WV_Global(LG).K      = 1.
        endif          
      endif
    enddo
    WV_Global(2:LA_Global).HEIGHT = WV_Global(2:LA_Global).HEISIG  ! *** ASSIGN INCIDENT SIGNIFICANT WAVE HEIGHT TO APPLIED WAVE HEIGHT
  endif
  
  ! *** Broadcast wave cells properties
  Call Broadcast_Array(LWVCELL_Global,   master_id)
  Call Broadcast_Array(LWVMASK_Global,   master_id)
  Call Broadcast_Array(WV_Global.HEIGHT, master_id)
  Call Broadcast_Array(WV_Global.DIR,    master_id)
  Call Broadcast_Array(WV_Global.LENGTH, master_id)
  Call Broadcast_Array(WV_Global.FREQ,   master_id)
  Call Broadcast_Array(WV_Global.K,      master_id)
  Call Broadcast_Array(WV_Global.HEISIG, master_id)
  Call Broadcast_Array(WVDISSIPA,        master_id)
  
  ! *** Map wave cells properties to each process
  NWVCELLS = 0
  do LG = 2, LA_Global
    ! *** Get local L
    L = Map2Local(LG).LL
    ! *** Valid if greater than zero
    if(L > 0 )then
      ! *** set local copy of this
      NWVCELLS = NWVCELLS + 1
  
      LWVCELL(NWVCELLS) = L
      LWVMASK(L)   = LWVMASK_Global(LG)
      WV(L).DIR    = WV_Global(LG).DIR
      WV(L).LENGTH = WV_Global(LG).LENGTH
      WV(L).FREQ   = WV_Global(LG).FREQ  
      WV(L).HEIGHT = WV_Global(LG).HEIGHT
      WV(L).K      = WV_Global(LG).K
      WV(L).HEISIG = WV_Global(LG).HEISIG
      WV(L).DISSIPA(KC) = WVDISSIPA(LG)
    endif
  enddo
      
  if( IWVCOUNT == NWVTIM )then
    close(WUNI)
    if( process_id == master_id ) write(*,'(A)') '** WAVE: END OF WAVE RECORD'
  endif 
  
END SUBROUTINE 
  
