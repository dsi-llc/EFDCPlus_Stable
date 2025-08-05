! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALTOX_KINETICS

  ! ***  SUBROUTINE CALTOX_KINETICS CALCULATES TOXICS DECAY FOR THE WATER COLUMN AND
  ! ***  THE SEDIMENT BED AND IS CALLED FROM SSEDTOX
  !
  !--------------------------------------------------------------------------------
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !--------------------------------------------------------------------------------
  ! 2012-12-05        PAUL M. CRAIG    RESTRUCTURED AND ADDED OMP
  ! 2017-09-01        PAUL M. CRAIG    REWRITTEN TO UPDATE FRAMEWORK TO ADD
  !                                    MULTIPLE KINETIC PROCESSES.
  !                                    ADDED BIODEGRADATAION AND VOLATILIZATION 
  !********************************************************************************

  ! *** ITOXKIN(1,NT) = BULK DECAY,         0-DO NOT USE, 1-USE
  ! *** ITOXKIN(2,NT) = BIODEGRADATION,     0-DO NOT USE, 1-USE
  ! *** ITOXKIN(3,NT) = VOLATILIZATION,     0-DO NOT USE, 1-SIMPLE, 2-COMPUTED
  ! *** ITOXKIN(4,NT) = PHOTOLYSIS,         0-DO NOT USE, 1-SIMPLE  (NOT IMPLEMENTED)
  ! *** ITOXKIN(5,NT) = HYDROLYSIS,         0-DO NOT USE, 1-SIMPLE  (NOT IMPLEMENTED)
  ! *** ITOXKIN(6,NT) = DAUGHTER PRODUCTS,  0-DO NOT USE, 1-SIMPLE  (NOT IMPLEMENTED)
  
  use GLOBAL

  implicit none
  
  integer :: K, L, NT, ND, LF, LL, LP, IDECAYB, IDECAYW, KBOT, KINC
  real :: COEFF, DENA, DENW, DA, DW, CDRAG, HE, SCA, SCW, TKA, TKW, TOXMW, VISCA, VISCW, USTR
  real :: U_AVG, V_AVG, VEL_AVG, K_L, K_G, KV, HPU, TOXDIS, VOLTERM, DEPTH, TIME
  real,save :: TOXTIME, RA, RCONST
  
  real :: PFTWC
  
  real,save,allocatable,dimension(:,:) :: CDECAYB
  real,save,allocatable,dimension(:,:) :: CDECAYW
  real,save,allocatable,dimension(:)   :: TLAYER
  
  real, external :: CALVOLTERM

  if( .not. allocated(CDECAYB) )then
    allocate(CDECAYB(LCM,KBM))  
    allocate(CDECAYW(LCM,KCM))  
    allocate(TLAYER(KCM))
    !ALLOCATE(KINTERM(LCM,KCM,NTXM))
  
    CDECAYB = 0.0
    CDECAYW = 0.0
    TLAYER = 0.0
    
    TOXTIME = 0.0
    RCONST = 8.206E-5                                                ! *** Universal gas constant (atm-m3/mole K)
    TOXSTEPW = TOXSTEPW - DELT/10.
  endif
  
  TOXTIME = TOXTIME + DTSED
  if( TOXTIME < TOXSTEPW ) return

  ! *** SET LAYER ORIENTATION
  if( LSEDZLJ )then
    KINC  = 1
    KBOT  = KB
  else
    KINC  = -1
    KBOT  = 1
  endif

  !KINTERM = 0.0    
  !RCONST = 8.3144
  IDECAYW = 0
  IDECAYB = 0

  ! *** LOOP OVER EACH TOXIC AND APPLY LOSS TERMS
  do NT = 1,NTOX
    ! *** INITIALIZE TOXIC CONSTANTS
    if( TOXS(NT).VOL.MW > 0. ) TOXMW = 1./TOXS(NT).VOL.MW**0.66667

    !$OMP PARALLEL DEFAULT(SHARED)
    
    ! *** BULK DECAY
    if( ITOXKIN(1,NT) > 0 )then
      ! *** COMPUTE DECAY COEFFICIENTS FOR EVERY TIMESTEP IF USING DYNAMIC TIME STEPPING
  
      ! *** BULK DECAY COEFFICIENT: WATER
      if( TOXS(NT).BLK_KW > 0. )then
        IDECAYW = 1
        !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,COEFF) 
        do ND = 1,NDM  
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              CDECAYW(L,K) = TOXS(NT).BLK_KW
            enddo
          enddo
        enddo   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      endif
      
      ! *** BULK DECAY COEFFICIENT: SEDIMENT BED
      if( TOXS(NT).BLK_KB > 0. )then
        IDECAYB = 1
        !$OMP DO PRIVATE(ND,LF,LL,L,K,COEFF,DEPTH) 
        do ND = 1,NDM  
          LF = 2+(ND-1)*LDM
          LL = min(LF+LDM-1,LA)
          do L = LF,LL
            DEPTH = 0.0
            do K = KBT(L),KBOT,KINC
              DEPTH = DEPTH + HBED(L,K)
              if( DEPTH > TOXS(NT).BLK_MXD ) CYCLE
              CDECAYB(L,K) = TOXS(NT).BLK_KB
            enddo
          enddo
        enddo   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      endif
    endif
  
    ! *** BIODEGRADATION
    if( ITOXKIN(2,NT) > 0 .and. ( ISTRAN(2) > 0 .or. IVOLTEMP > 0 ) )then
      ! *** DECAY COEFFICIENT: WATER
      if( TOXS(NT).BIO_KW > 0. )then
        IDECAYW = 1
        !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,COEFF) 
        do ND = 1,NDM  
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              COEFF = TOXS(NT).BIO_KW*TOXS(NT).BIO_Q10W**( 0.1*(TEM(L,K) - TOXS(NT).BIO_TW) )
              CDECAYW(L,K) = CDECAYW(L,K) + COEFF
            enddo
          enddo
        enddo   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      endif
      
      ! *** DECAY COEFFICIENT: SEDIMENT BED
      if( TOXS(NT).BIO_KB > 0. )then
        IDECAYB = 1
        if( TEMBO > 0.0 .and. (HTBED1 + HTBED2) > 0.0 .and. IVOLTEMP == 0 )then
          ! *** BED TEMPERATURES ARE SIMULATED, SO USE
          !$OMP DO PRIVATE(ND,LF,LL,L,K,COEFF,DEPTH) 
          do ND = 1,NDM  
            LF = 2+(ND-1)*LDM
            LL = min(LF+LDM-1,LA)
            do L = LF,LL
              DEPTH = 0.0
              do K = KBT(L),KBOT,KINC
                DEPTH = DEPTH + HBED(L,K)
                if( DEPTH > TOXS(NT).BIO_MXD ) CYCLE
                COEFF = TOXS(NT).BIO_KB * TOXS(NT).BIO_Q10B**( 0.1*(TEMB(L) - TOXS(NT).BIO_TB) )
                CDECAYB(L,K) = CDECAYB(L,K) + COEFF
              enddo
            enddo
          enddo   ! *** END OF DOMAIN LOOP
          !$OMP END DO
        else
          ! *** BED TEMPERATURES ARE NOT SIMULATED, USE BOTTOM LAYER WATER TEMPERATURES
          !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,COEFF,DEPTH) 
          do ND = 1,NDM  
            LF = 2+(ND-1)*LDM
            LL = min(LF+LDM-1,LA)
            do L = LF,LL
              DEPTH = 0.0
              do K = KBT(L),KBOT,KINC
                DEPTH = DEPTH + HBED(L,K)
                if( DEPTH > TOXS(NT).BIO_MXD ) CYCLE
                COEFF = TOXS(NT).BIO_KB*TOXS(NT).BIO_Q10B**( 0.1*(TEM(L,KSZ(L)) - TOXS(NT).BIO_TB) )
                CDECAYB(L,K) = CDECAYB(L,K) + COEFF
              enddo
            enddo
          enddo   ! *** END OF DOMAIN LOOP
          !$OMP END DO
        endif
      endif  ! *** END OF BED DEGRADATION BLOCK  
    endif    ! *** END OF BIODEGRADATION BLOCK

    if( IDECAYW == 1 )then
      ! *** APPLY DEGRADATION TO THE WATER COLUMN
      !$OMP DO PRIVATE(ND,LP,L,K) 
      do ND = 1,NDM  
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOX(L,K,NT) = TOX(L,K,NT)*(1.0 - TOXTIME*CDECAYW(L,K))
            CDECAYW(L,K) = 0.0
          enddo
        enddo
      enddo
      !$OMP END DO
    endif

    if( IDECAYB == 1 )then
      ! *** APPLY DEGRADATION TO THE SEDIMENT BED
      !$OMP DO PRIVATE(ND,LL,LF,L,K) 
      do ND = 1,NDM  
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          do K = KBOT,KBT(L),-KINC
            TOXB(L,K,NT) = TOXB(L,K,NT)*(1.0 - TOXTIME*CDECAYB(L,K))
            CDECAYB(L,K) = 0.0
          enddo
        enddo
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    endif
    
    if( TOXS(NT).VOL.KL_OPT > 0 )then
      ! ***    VOLATILIZATION
      !$OMP DO PRIVATE(ND, LP, L, PFTWC, VOLTERM)
      do ND = 1,NDM  
        do LP = 1,LLWET(KC,ND)
          L = LKWET(LP,KC,ND) 
          
          ! *** Use CALVOLTERM function to calculate volatilization term
          PFTWC = TOXPFTW(L,KC,NT)              ! *** Solids partitioning fraction
          VOLTERM = CALVOLTERM(HP(L), HPKI(L,KC), STCUV(L), UHE(L), VHE(L), HU(L), HV(L), TEM(L,KC), SAL(L,KC),     &
             TATMT(L), PATMT(L), WINDST(L), VOL_VEL_MAX, VOL_DEP_MIN, TOXMW, TOXS(NT).VOL.HE, TOXS(NT).VOL.AIRCON,  &
             TOXS(NT).VOL.TCOEFF, TOXS(NT).VOL.MULT, TOX(L,KC,NT), TOXS(NT).VOL.KL_OPT, PFTWC)
                         
          TOX(L,KC,NT) = TOX(L,KC,NT) - VOLTERM*TOXTIME
          TOX(L,KC,NT) = max(TOX(L,KC,NT), 0.0)
          
        enddo
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    endif     ! *** END OF VOLATILIZATION
  
    if( ITOXKIN(4,NT) > 0 )then
      ! ***    PHOTOLOSIS
      ! ***      RKTOXP(NT) = BASE RATE
      ! ***      SKTOXP(NT) = SOLAR RADIATION AT BASE RATE
    endif
  
    !$OMP END PARALLEL
  
  enddo   ! *** END OF NTOX LOOP
    
  TOXTIME = 0.0
  
  return

  END

!********************************************************************************
!  Subroutine for calculating water dynamic viscosity from temperature
!********************************************************************************
REAL FUNCTION ViscosityW(TEMP)

  ! *** TEMP = Temperature at which to calculate water dynamic viscosity
  ! *** VISC = Calculated dynamic (absolute) viscosity of water  (kg/m/s)

  real, intent(IN) :: TEMP
  
  real    :: VISC, TK

  ! *** Units conversions
  ! *** 1.0 kg/m/s  =  1.0 Pa*s 
  ! *** 1.0 kg/m/s  =  1.0 N*s/m^2
  ! *** 1.0 kg/m/s  =  10.0 Poise (g/cm/s)
  ! *** 1.0 kg/m/s  =  1000.0 cP

  ! *** LINEARLY INTERPOLATE BETWEEN TABLE VALUES
  if( TEMP < 0. )then
    ! *** CALHEAT handles this case, so just apply minimum
    VISC = 0.0017914     ! *** Viscosity at 0C

  elseif( TEMP >= 100. )then
    write(6,*) 'Water temperature abovE 100 deg.C. encountered!'
    write(6,*) 'EFDC+ stopped in ViscosityW'
    call STOPP('.')
  
  else
    TK = TEMP + 273.15                       ! *** Convert to Kelvin
    VISC = 2.414E-5*10**(247.8/(TK-140))     ! *** Dynamic viscosity (kg/m/s) (Pa*s)
  endif
  ViscosityW = VISC
  
END FUNCTION ViscosityW
