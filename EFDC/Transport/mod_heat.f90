! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE HEAT_MODULE

! CHANGE RECORD
! DATE MODIFIED     BY               DESCRIPTION
!----------------------------------------------------------------------!
!    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
!    2015-01       Paul M. Craig     Added fully coupled Ice Sub-model
!                  Dang H Chung
! 2014-05-08       Paul M. Craig     Converted to Module and fixed OMP issue with equilibrium temperature

use GLOBAL
use INFOMOD, only:READSTR
use DIFFUSER_MODULE, only:FUNDEN

use Variables_WQ
use WATERQUALITY, only:ALGAES, LAYERBOT, LAYERTOP

use Variables_MPI
use Variables_MPI_Write_Out
use Broadcast_Routines

implicit none

! ******* parameter DECLARATIONS
real, save, PRIVATE :: MPS_TO_MPH          = 2.23714     ! *** W2 equilibrium temperature method (English units)
real, save, PRIVATE :: W_M2_TO_BTU_FT2_DAY = 7.60796     ! *** W2 equilibrium temperature method (English units)
real, save, PRIVATE :: FLUX_BR_TO_FLUX_SI  = 0.23659     ! *** W2 equilibrium temperature method (English units)
real, save, PRIVATE :: BTU_FT2_DAY_TO_W_M2 = 0.1314      ! *** W2 equilibrium temperature method (English units)
real, save, PRIVATE :: BCONV               = 1.520411    ! *** W2 equilibrium temperature method (English units)

real, save, PRIVATE :: LHF     = 333507.0                ! *** LATENT HEAT OF FUSION FOR ICE (J/KG)
real, save, PRIVATE :: CP      = 4179.0                  ! *** SPECIFIC HEAT (J/KG/degC)    (Previously EFDC used 4179.0)   4184

real, save, PRIVATE :: REFICE                            ! ***
real, save, PRIVATE :: RHOWCP                            ! *** 1393725753 (W s)/(m^3 degC)
real, save, PRIVATE :: RHOWCPI                           ! *** 0.2393E-6 (m^3 degC)/(W s)
real, save, PRIVATE :: RHOICP                            ! ***
real, save, PRIVATE :: RHOICPI                           ! ***
real, save, PRIVATE :: CREFLI                            ! ***
real, save, PRIVATE :: RHOILHFI                          ! ***

contains

SUBROUTINE CALHEAT

  !   Subroutine CALHEAT takes the information from the atmospheric boundary
  !   file and the wind forcing file and calculates the net heat flux across
  !   the water surface boundary.  The heat flux is then used to update the
  !   water temperature either in the surface cells, or distributed across
  !   the cells in the vertical and into the bottom.
  !
  !   The heat flux terms are derived from a paper by Rosati
  !   and Miyakoda (1988) entitled "A General Circulation Model for Upper Ocean
  !   Simulation".  The heat flux is prescribed by term for the following
  !   influxes and outfluxes:
  !
  !     - Short Wave Incoming Radiation (+)
  !     - Net Long Wave Radiation (+/-)
  !     - Sensible Heat Flux (convection -)
  !     - Latent Heat Flux (evaporation +/-)
  !
  !   Two formulations of the Latent Heat Flux are provided.  The first is from
  !   the Rosati and Miyakoda paper, the second is an alternate formulation by
  !   Brady, Graves, and Geyer (1969).  The second formulation was taken from
  !   "Hydrodynamics and Transport for Water Quality Modeling" (Martin and
  !   McCutcheon, 1999).  The Rosati and Miyakoda formulation will have zero
  !   evaporative cooling or heating if wind speed goes to zero.  The Brady,
  !   Graves, and Geyer formulation provides for a minimum evaporative cooling
  !   under zero wind speed.
  !
  !
  ! VARIABLE LIST:
  !
  !   CLOUDT  = Cloud Cover (0 to 10)
  !   HCON    = Sensible Heat Flux (W/m2)
  !   HLAT    = Latent Heat Flux (W/m2)
  !   HLWBR   = Net Longwave Radiation (Atmospheric Long Wave plus Back Radiation, W/m2)
  !   SOLSWRT = Total Short Wave Incoming Radiation (W/m2)
  !   SVPW    = Saturation Vapor Pressure in mb Based on the Water Surface Temperature (mb)
  !   TATMT   = Temperature of Air Above Water Surface (deg C)
  !   TEM     = Water Temperature in Cell (deg C)
  !   VPAT    = Vapor Pressure of Air at Near Surface Air Temperature (mb)
  !   WINDST  = Wind Speed at 2 Meters Over Cell Surface (m/s)
  !
  ! MODIFICATION HISTORY:
  !
  !   Date       Author             Comments*6
  !   ---------- ------------------ -----------------------------------------------------------------
  !   2015-01    Paul M. Craig      Added Fully Heat Coupled Ice Sub-model
  !              Dang Chung
  !   2014-12    Paul M. Craig      Added the new Evaporation Approach
  !   2012-10    Paul M. Craig      Added DOC to the light extinction and
  !                                 standardized for all other routines
  !   2011-03    Paul M. Craig      Converted to F90, added OMP
  !   11/01/2005 Paul M. Craig      Added Option 4, to use the Equilibrium Temperature
  !                                 algorithym from CE-QUAL-W2.  Also added the sub-option
  !                                 under this option to couple or decouple the bottom temperature
  !                                 to the water column temperatures.
  !                                 Added the ability to input spatially variable bed temps and
  !                                 thermally active bed thicknesses.
  !                                 Also cleaned up the code and added structure.
  !   11/07/2000 Steven Peene       Cleaned code, provided more detailed
  !                                 descriptions, added alternate formulation
  !                                 for latent heat flux, separated out
  !                                 individual heat flux terms
  !   06/01/1992 John M. Hamrick    Original author

  use INFOMOD,only:SKIPCOM
  use Variables_MPI
  use Broadcast_Routines
  use Allocate_Initialize

  integer :: L, LF, LG, LL, LP, K, I, J, L1, ND, IOS, NAL

  integer,save,allocatable,dimension(:) :: ICOUNT

  real       :: T1, T2, SVPW1, DTHEQT, WCKESS, TFLUX, BOT, FRACLAYER
  real       :: TFAST, TFAST1, TSLOW, TSLOW1, RSN, C1, C2, UBED, VBED
  real       :: USPD, TMPKC, SRO, SRON
  real       :: THICK, TMIN, FLUXQB
  real       :: ET, CSHE, TEMP, RICETHKL0, ET0, CSHE0, WSHLTR0, PSHADE0
  real       :: HBLW, HBCD, HBCV, HBEV, HBNT, HBLWW
  real,save  :: TIMEDRY

  real,save,allocatable,dimension(:) :: TBEDTHK                ! *** Thermal Thickness of the bed (m)
  real,save,allocatable,dimension(:) :: FLUXTB                 ! *** Accumulator for conductive and convective heat exchange with the bed  (C*m)
  real,save,allocatable,dimension(:) :: RADBOTT                ! *** Accumulator for solar radiation at the bed surface                    (W/m2)
  real,save,allocatable,dimension(:) :: PSHADE_OLD
  real,save,allocatable,dimension(:) :: WSHLTR_OLD
  real,save,allocatable,dimension(:) :: TMIND
  real,save,allocatable,dimension(:) :: ALGMOBILE
  character(200) :: STR

  ! *** COARE variable
  real,save :: zuu   !< height of wind measurement [m]
  real,save :: ztt   !< height of temperature sensor [m]
  real,save :: zqq   !< height of RH sensor [m]
  real,save :: zii   !< PBL height [m]
  real :: tau_, hsb, hlb, Le_, sw_net, lw_dn_abs, lw_up, rainn, rhh
  real :: zoo, Ss, vcp, sigH,dz_skin, cd_

  if( .not. allocated(TBEDTHK) )then
    call AllocateDSI(TBEDTHK,    LCM, 0.0)
    call AllocateDSI(FLUXTB,     LCM, 0.0)
    call AllocateDSI(RADBOTT,    LCM, 0.0)

    call AllocateDSI(PSHADE_OLD, NDM, 0.0)
    call AllocateDSI(WSHLTR_OLD, NDM, 1.0)
    call AllocateDSI(TMIND,      NDM, 0.0)
    call AllocateDSI(ICOUNT,     NDM, 1)

    if( ISTRAN(8) > 0 .and. NALGAE > 0 )then
      allocate(ALGMOBILE(NALGAEM))
      ALGMOBILE = 0.0

      do NAL = 1,NALGAE
        if( ALGAES(NAL).ISMOBILE )then
          ALGMOBILE(NAL) = ALGAES(NAL).WQCHLA
        endif
      enddo
    endif

    if( process_id == master_id )then
      write(*,'(A)')'CALHEAT: INITIALIZING'
    endif

    TIMEDRY = TBEGIN

    ! *** Initialize Heat Exchange Terms
    RADNET  = 0.
    TBEDTHK = ABS(TEMTHKO)
    if( .not. ( ISRESTI > 0 .and. ISCI(2) == 1 )  )then
      TEMB        = ABS(TEMBO)
      TEMB_Global = ABS(TEMBO)
    endif

    if( TEMTHKO > 0. .and. .not. ( ISRESTI > 0 .and. ISCI(2) == 1 )  )then
      call AllocateDSI(R1D_Global, LCM_Global, 0.0)

      if( process_id == master_id )then
        ! *** READ IN THE SPATIALLY VARYING INIT T AND BED THICKNESS (TEMTHKO)
        write(*,'(A)')' CALHEAT: READ IN THE SPATIALLY VARYING INIT T AND BED THICKNESS: TEMB.INP'
        open(1001,FILE = 'temb.inp',ERR = 1000,STATUS = 'OLD')

        STR = READSTR(1001)
        do L1 = 2,LA_Global
          read(1001,*,END = 1000) I, J, T1, T2
          LG = LIJ_Global(I,J)
          T1 = max(MIN(T1,50.),0.)                    ! *** Bed emperature, Prevent bad IC settings (EE should catch this)
          TEMB_Global(LG) = T1
          T2 = max(T2,0.1)                            ! *** Thermal thickness, Prevent bad IC settings (EE should catch this)
          R1D_Global(LG)  = T2
        enddo
1000    close(1001)
      endif

      call Broadcast_Array(TEMB_Global, master_id)
      call Broadcast_Array(R1D_Global, master_id)

      do LG = 2,LA_Global
        L = Map2Local(LG).LL
        if( L > 0 )then
          TEMB(L)    = TEMB_Global(LG)
          TBEDTHK(L) = R1D_Global(LG)
        endif
      enddo

      deallocate(R1D_Global)
    endif

    PSHADE = 1.0          ! *** DEFAULT SHADE FACTOR

    if( USESHADE )then
      if( process_id == master_id )then
        ! *** READ IN SPATIALLY VARYING SHADING FACTORS
        write(*,'(A)')' CALHEAT: READ IN SPATIALLY VARYING SHADE: PSHADE.INP'
        open(1001,FILE = 'pshade.inp',ERR = 1010,STATUS = 'OLD')

        STR = READSTR(1001)

        do L1 = 2,LA_Global
          read(1001,*,END = 1010) I, J, T1
          LG = LIJ_Global(I,J)
          SHAD_Global(LG) = T1
        enddo
1010    close(1001)
      endif

      call Broadcast_Array(SHAD_Global, master_id)

      do LG = 2,LA_Global
        L = Map2Local(LG).LL
        if( L > 0 )then
          PSHADE(L) = SHAD_Global(LG)
        endif
      enddo
    else
      if( process_id == master_id )then
        write(*,'(A)')' CALHEAT: SETTING CONSTANT SHADE TO: 1.0 (NO SHADE)'
      endif
    endif

    ! *** READING ICEINIT.INP FILE FOR THE INITIAL ICE THICKNESS
    if( ISICE > 2 )then
      if( ISRESTI == 0 )then
        if( process_id == master_id )then
          write(*,'(A)')' CALHEAT: READ IN SPATIALLY VARYING ICE CONDITIONS: ICE.INP'
          open(1001,FILE = 'ice.inp',STATUS = 'OLD')
          call SKIPCOM(1001,'*')
          do LL = 2,LA_Global
            read(1001,*,IOSTAT = IOS) I, J, RICETHKL0
            if( IOS > 0 ) STOP 'ICE.INP: READING ERROR'
            LG = LIJ_Global(I,J)
            ICETHICK_Global(LG) = RICETHKL0        ! INITIAL ICE THICKNESS IN M
          enddo
          close(1001)
        endif
        call Broadcast_Array(ICETHICK_Global, master_id)

        ! *** Map to Local Domain
        do LG = 2,LA_GLOBAL
          L = Map2Local(LG).LL
          if( L > 1 )then
            ICETHICK(L) = ICETHICK_Global(LG)
          endif
        enddo
      endif

      do L = 2,LA
        if( ICETHICK(L)  >= MINICETHICK )then
          ICECELL(L) = .TRUE.
          ICECOVER(L) = 1.0            ! *** Fraction of ice cover over cell: 0 to 1.0
          LFRAZIL = .TRUE.
        else
          ICECELL(L) = .FALSE.
          ICEVOL(L) = 0.0
        endif
      enddo
    endif

    ICOUNT  = 1  ! *** FORCE INITIAL ZEROING

    ! *** INITIALIZE EXTINCTION COEFFICIENTS
    if( IASWRAD == 3 )then
      do K = 1,KC
        do L = 2,LA
          WCKESS = GET_EXT_COEFF(L,K)
        enddo
      enddo
    endif
    RHOWCP  = 1000.0*CP             ! *** CP-SPECIFIC HEAT (J/KG/degC) * RHO (KG/M3)
    RHOWCPI = 1./RHOWCP
    if( ISICE < 3 ) ALBEDOI = 1.0   ! *** Reflect 100% of solar radiation for ice cells
    
    if( DEBUG )then
      write(6,'(A,F8.1)')' CALHEAT: Bed Temp(L = 2):', TEMB(2), process_id

      !OPEN(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
      !WRITE(mpi_efdc_out_unit,'(a)') 'CALHEAT Diagnostics'
      !WRITE(mpi_efdc_out_unit,'(A11,8A9)')'TIMEDAY','SRON','ET','TD_C','TA_C','TDEW_F','TAIR_F','FW'
      !CLOSE(mpi_efdc_out_unit)
    endif

    ! *** Initialize COARE variables
    zuu = 2.   !< Wind speed is at 2m
    ztt = TEM_HRZ
    zqq = TEM_HRZ
    zii = PBLZ
  endif    ! *** END OF CALHEAT INITIALIZATION BLOCK

  if( ISTL == 2 )then
    if( ISDYNSTP == 0 )then
      DELT = DT
    else
      DELT = DTDYN
    endif
  else
    DELT = DT2
  endif
  TMIND = 9999.

  ! *** OVERWRITE THE INPUT SOLAR RAD WITH A COMPUTED ONE
  if( COMPUTESOLRAD )then
    call SHORT_WAVE_RADIATION(CLOUDT(2), SRO, SRON)
    if( SRON > 1.0E-4 )then
      LDAYLIGHT = .TRUE.
    else
      LDAYLIGHT = .FALSE.
    endif
  endif

  ! *************************************************************************************************************
  ! *************************************************************************************************************
  ! *** PREPARE SOLAR RADIATION AND SURFACE PROPERTIES BEFORE
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, K, NAL)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET+1
    LL = min(LF+LDMWET-1,LAWET)

    if( COMPUTESOLRAD )then
      ! *** OVERWRITE THE INPUT SOLAR RAD WITH A COMPUTED ONE
      do LP = LF,LL
        L = LWET(LP)
        SOLSWRT(L) = SRON
      enddo
    endif

    if( LDAYLIGHT )then
      if( (IASWRAD == 3 .or. IASWRAD == 0) .and. ISTRAN(8) > 0 )then
        ! *** If using WQ then use the WQ Coefficients for light extinction
        ! *** Water Quality is Active so account for Chlorophyll
        ! *** Compute WQCHL (Chlorophyll) Using Algal Biomass & factors
        ! *** IASWRAD == 0 is included for Chapra test-case
        if( ISWQLVL == 0 )then
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)
              WQCHL(L,K) = WQV(L,K,1)*ALGAES(1).WQCHLA + WQV(L,K,2)*ALGAES(2).WQCHLA + WQV(L,K,3)*ALGAES(3).WQCHLA
            enddo
          enddo
        else
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)

              WQCHL(L,K) = 0.0
              do NAL = 1,NALGAE
                WQCHL(L,K) = WQCHL(L,K) + WQV(L,K,NAL+19)*ALGMOBILE(NAL)
              enddo
            enddo
          enddo
        endif
      endif
      
      ! *** SET SOLAR RADIATION AT SURFACE
      ICOUNT(ND) = ICOUNT(ND) + 1
      
      if( IASWRAD > 0 .or. ISTRAN(8) > 0 )then
        if( USESHADE )then
          do LP = LF,LL
            L = LWET(LP)
            ! *** APPLY PSHADE FACTORS
            SOLSWRT(L) = SOLSWRT(L)*PSHADE(L)
          enddo
        endif
        
        do LP = LF,LL
          L = LWET(LP)
          call SET_LIGHT(L)
        enddo
      endif
      
    elseif( ICOUNT(ND) > 0 )then
      ! *** NO LIGHT CASE (ONLY ZERO IF LAST TIMESTEP WAS DAYLIGHT)
      ICOUNT(ND) = 0
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          RADBOT(L,K) = 0.
          RADTOP(L,K) = 0.
          RADNET(L,K) = 0.
        enddo
      enddo
      
      ! *** ZERO SURFACE RADIATION (ONLY ONCE PER ITERATION)
      if( ND == 1 )then
        do LP = LF,LL
          L = LWET(LP)
          RADTOP(L,0) = 0.
        enddo
      endif
    else
      ! *** ALWAYS ZERO TOP NET RADIATION
      K = KC
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        RADBOT(L,K) = 0.
        RADTOP(L,K) = 0.
        RADNET(L,K) = 0.
      enddo
    endif   ! *** END OF LDAYLIGHT BLOCK

  enddo   ! *** END OF DOMAIN
  !$OMP END PARALLEL DO
  
  ! *************************************************************************************************************
  ! *************************************************************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  
  ! *** Computation of Solar Radiation
  if( IASWRAD == 0 )then                         ! *** Ignore solar radiation
  
    if( ISTOPT(2) == 0 .and. ISTRAN(8) > 0 )then
      ! *** Computing the light extinction for Chapra test cases
      !$OMP DO PRIVATE(ND,K,LP,L,C1,RSN)
      do ND = 1,NDM
        do K = KC,1,-1
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            
            ! *** Update Chl-a regardless of LDAYLIGHT
            WQCHL(L,K) = 0.0
            do NAL = 1,NALGAE
              WQCHL(L,K) = WQCHL(L,K) + WQV(L,K,NAL+19)*ALGMOBILE(NAL)
            enddo
            WCKESS = GET_EXT_COEFF(L,K)
              
            ! *** Computing fraction of light at bottom of the current layer for Visser, 1997 test-cases
            BOT=MAX(-WCKESS*HPK(L,K),-40.0)
            FRACLAYER=EXP(BOT)
            RADBOT(L,K) = RADTOP(L,K)*FRACLAYER
            
            ! *** Compute Net Energy (W/M2)
            RADNET(L,K) = RADTOP(L,K) - RADBOT(L,K)
            
            ! *** Update layer variables
            RADTOP(L,K-1) = RADBOT(L,K)
          enddo
        enddo
      enddo
      !$OMP END DO
    
    else
      ! ***  
      !$OMP DO PRIVATE(ND,K,LP,L,C1,RSN)
      do ND = 1,NDM
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            RADBOT(L,K) = 0.
            RADTOP(L,K) = 0.
            RADNET(L,K) = 0.
          enddo
        enddo
      enddo
      !$OMP END DO
      
    endif
    
  elseif( IASWRAD == 1 )then                     ! *** Absorb 100% solar radiation in top layer
  
    if( LDAYLIGHT )then
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,C1,RSN)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)
          do LP = LF,LL
            L = LWET(LP)
            C1 = DELT*DZIC(L,KC)*RHOWCPI
            RSN = RADTOP(L,KC)
            RADNET(L,KC) = RSN*C1
            RADBOT(L,KC) = 0.
          enddo
      enddo   ! *** END OF DOMAIN
      !$OMP END DO
    endif

  elseif( IASWRAD == 2 )then                     ! *** Use fast/slow extinction coefficients
  
    if( LDAYLIGHT )then
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,TFAST,TFAST1,TSLOW,TSLOW1,C1,RSN)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)

        ! *** SURFACE LAYER
        do LP = LF,LL
          L = LWET(LP)
          TFAST  = SWRATNF*(Z(L,KC)  -1.)
          TFAST1 = SWRATNF*(Z(L,KC-1)-1.)
          TSLOW  = SWRATNS*(Z(L,KC)  -1.)
          TSLOW1 = SWRATNS*(Z(L,KC-1)-1.)
          C1 = DELT*DZIC(L,KC)*RHOWCPI
          RSN = RADTOP(L,KC)*( FSWRATF*EXP(TFAST*HP(L))  + (1.-FSWRATF)*EXP(TSLOW*HP(L))        &
                             - FSWRATF*EXP(TFAST1*HP(L)) - (1.-FSWRATF)*EXP(TSLOW1*HP(L)) )
                        
          RADNET(L,KC) = RSN*C1                                                  ! *** m*degC
          RADBOT(L,KC) = RADTOP(L,KC) - RSN                                      ! *** Radiation at the bottom of the layer (W/m2)
        enddo
      enddo   ! *** END OF DOMAIN
      !$OMP END DO

      ! *** ALL REMAINING LAYERS
      if( KC > 1 )then
        !$OMP DO PRIVATE(ND,K,LP,L,TFAST,TFAST1,TSLOW,TSLOW1,C1,RSN)
        do ND = 1,NDM
          do K = KS,1,-1
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)

              TFAST  = SWRATNF*(Z(L,K)-1.)
              TFAST1 = SWRATNF*(Z(L,K-1)-1.)
              TSLOW  = SWRATNS*(Z(L,K)-1.)
              TSLOW1 = SWRATNS*(Z(L,K-1)-1.)
              C1 = DELT*DZIC(L,K)*RHOWCPI
              RSN = RADTOP(L,KC)*( FSWRATF*EXP(TFAST*HP(L))  + (1.-FSWRATF)*EXP(TSLOW*HP(L))      &
                                 - FSWRATF*EXP(TFAST1*HP(L)) - (1.-FSWRATF)*EXP(TSLOW1*HP(L)) )

              RADNET(L,K) = RSN*C1                                                                ! *** m*degC
              RADTOP(L,K) = RADBOT(L,K+1)                                                         ! *** W/m2
              RADBOT(L,K) = RADTOP(L,K) - RSN                                                     ! *** W/m2
            enddo
          enddo
        enddo   ! *** END OF DOMAIN
        !$OMP END DO
      endif
    endif   ! *** END OF DAYLIGHT BLOCK

  elseif( IASWRAD == 3 )then                     ! *** Spatially & temporally varying extinction coefficients
  
    if( LDAYLIGHT )then
      !$OMP DO PRIVATE(ND,K,LF,LL,LP,L,WCKESS,BOT,FRACLAYER,RSN,C1)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)

        ! *** SURFACE LAYER
        K = KC
        do LP = LF,LL
          L = LWET(LP)
          C1 = DELT*DZIC(L,K)*RHOWCPI
          ! *** Extinction Coefficient
          WCKESS = GET_EXT_COEFF(L,K)

          ! *** FRACTION OF LIGHT AT THE BOTTOM OF THE CURRENT LAYER
          BOT = max(-WCKESS*HPK(L,K),-40.0)
          FRACLAYER = EXP(BOT)
          if( FSOLRADMIN > 0.0 )then
            BOT = 1.0 - FRACLAYER
            BOT = max(BOT, FSOLRADMIN)
            FRACLAYER = 1.0 - BOT
          endif
          
          RADBOT(L,K) = RADTOP(L,K)*FRACLAYER

          ! *** Compute Net Energy (W/M2)
          RSN = RADTOP(L,K) - RADBOT(L,K)
          RADNET(L,K) = RSN*C1

          ! *** Update layer variables
          RADTOP(L,K-1) = RADBOT(L,K)
        enddo

        if( KC > 1 )then
          do K = KS,1,-1
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)

              C1 = DELT*DZIC(L,K)*RHOWCPI
              ! *** Extinction Coefficient
              WCKESS = GET_EXT_COEFF(L,K)

              ! *** FRACTION OF LIGHT AT THE BOTTOM OF THE CURRENT LAYER
              BOT = max(-WCKESS*HPK(L,K),-40.0)
              FRACLAYER = EXP(BOT)

              RADBOT(L,K) = RADTOP(L,K)*FRACLAYER

              ! *** Compute Net Energy (W/M2)
              RSN = RADTOP(L,K) - RADBOT(L,K)
              RADNET(L,K) = RSN*C1

              ! *** Update layer variables
              RADTOP(L,K-1) = RADBOT(L,K)

            enddo
          enddo
        endif
      enddo
      !$OMP END DO
    endif
    
  endif      ! *** End of solar radiation options
  !$OMP END PARALLEL

  ! *************************************************************************************************************
  ! *************************************************************************************************************
  
  ! *** Conduct heat transfer at surface, bottom and solar radiation extinction, based on selected method
  
  if( ISTOPT(2) == 0 )then                           ! *** No heat transfer at surface or bed
  
    return
  
  elseif( ISTOPT(2) == 1 )then                       ! *** Full heat balance
  
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,C2,SVPW1,HBLW,HBCV,HBEV,TEMP)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)

      ! *** SURFACE HEAT FLUX WITHOUT SOLAR RADIATION (m * degC)
      do LP = LF,LL
        L = LWET(LP)
        C2 = -DELT*DZIC(L,KC)
        if( .not. ICECELL(L) )then
          ! *** CELL WITHOUT ICE
          SVPW1 = SVPW(L)
          HBLW = 1.312E-14*((TEM(L,KC)+273.)**4)*(0.39-0.05*SQRT(VPAT(L)))*(1.-.8*CLOUDT(L)) + &
            5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))
          HBCV = CCNHTT(L)*0.288E-3*WINDST(L)*(TEM(L,KC)-TATMT(L)) 
          HBEV = CLEVAP(L)*0.445*WINDST(L)*(SVPW1-VPAT(L))/PATMT(L)
          
          RADNET(L,KC) = RADNET(L,KC) + C2*(HBLW + HBCV + HBEV)
                    
          ! *** Writing in array out for diagnostic
          if( ISINWV == 2 )then
            HS_OUT(L) = HBCV    ! *** m*degC/s
            HL_OUT(L) = HBEV    ! *** m*degC/s
            HW_OUT(L) = HBLW    ! *** m*degC/s
          endif
                  
          ! *** FINALIZE TEMPERATURE
          TEMP = TEM(L,KC) + HPI(L)*RADNET(L,KC)
          TEM(L,KC) = TEMP
        else
          ! *** FINALIZE TEMPERATURE
          TEMP = TEM(L,KC) + HPI(L)*RADNET(L,KC)
          TEM(L,KC) = TEMP
        endif
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO

    ! *** NOW FINALIZE THE TEMPERATURE
    if( LDAYLIGHT .and. KC > 1 )then
      !$OMP DO PRIVATE(ND,K,LP,L,TEMP)
      do ND = 1,NDM
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            TEMP = TEM(L,K) + HPI(L)*RADNET(L,K)
            TEM(L,K) = TEMP
          enddo
        enddo
      enddo   ! *** END OF DOMAIN
      !$OMP END DO
    endif

    !$OMP END PARALLEL

  elseif( ISTOPT(2) == 2 )then                   ! *** COARE 3.6
  
    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,C2,HBLW,HBLWW,rainn,rhh,Ss,vcp,sigH,HBCV,HBEV, TEMP) &
    !$OMP    PRIVATE(tau_,hsb,hlb,Le_,sw_net,lw_dn_abs,lw_up,zoo,cd_,dz_skin)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)
      ! *** SURFACE HEAT FLUX WITHOUT SOLAR RADIATION (m * degC)
      do LP = LF,LL
        L = LWET(LP)
        C2 = -DELT*DZIC(L,KC)
        if( .not. ICECELL(L) )then
          ! *** Preparing input arguments for COARE subroutine
          HBLW = 1.312E-14*((TEM(L,KC)+273.)**4)*(0.39-0.05*SQRT(VPAT(L)))*(1.-.8*CLOUDT(L)) + &
            5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))

          ! *** convert to COARE unit
          HBLWW = RHOWCP*HBLW
          rainn = RAINT(L)*3600000.
          rhh = RHAT(L)*100.

          if( ISTOPT(1)  >= 1 )then
            Ss = SAL(L,KC)
          else
            Ss = 0.
          endif

          if( ISWAVE  >= 3 .and. WV(L).PERIOD > 0 )then
            vcp   = WV(L).LENGTH/WV(L).PERIOD
            sigH  = WV(L).HEIGHT
          else
            vcp  = 0.
            sigH = 0.
          endif

          ! *** call subroutine coare 36
          call coare36flux_coolskin(WINDST(L),zuu,TATMT(L),ztt,rhh,zqq,PATMT(L),TEM(L,KC),SOLSWRT(L),HBLWW,G,zii, &
            rainn,Ss,vcp,sigH,tau_,hsb,hlb,Le_,sw_net,lw_dn_abs,lw_up,zoo,cd_,dz_skin)

          ! *** convert back to EFDC unit
          HBCV = RHOWCPI*hsb    ! *** m*degC/s
          HBEV = RHOWCPI*hlb    ! *** m*degC/s
          RADNET(L,KC) = RADNET(L,KC) + C2*(HBLW + HBCV + HBEV)
          CDCOARE(L)  = cd_
          ZSRE(L)     = zoo
          EVACOARE(L) = hlb/Le_/1000. !*** evap rate m/s

          ! *** Writing in array out for diagnostic
          if( ISINWV == 2 )then
            HS_OUT(L) = HBCV    ! *** m*degC/s
            HL_OUT(L) = HBEV    ! *** m*degC/s
            HW_OUT(L) = HBLW    ! *** m*degC/s
          endif
          
          ! *** FINALIZE TEMPERATURE
          TEMP = TEM(L,KC) + HPI(L)*RADNET(L,KC)
          TEM(L,KC) = TEMP
        else
          ! *** FINALIZE TEMPERATURE
          TEMP = TEM(L,KC) + HPI(L)*RADNET(L,KC)
          TEM(L,KC) = TEMP
        endif
      enddo
    enddo
    !$OMP END DO
    
    ! *** NOW FINALIZE THE TEMPERATURE
    if( LDAYLIGHT .and. KC > 1 )then
      !$OMP DO PRIVATE(ND,K,LP,L,TEMP)
      do ND = 1,NDM
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)
            TEMP = TEM(L,K) + HPI(L)*RADNET(L,K)
            TEM(L,K) = TEMP
          enddo
        enddo
      enddo   ! *** END OF DOMAIN
      !$OMP END DO
    endif
    
    !$OMP END PARALLEL

  elseif( ISTOPT(2) == 3 )then                   ! *** W2 Equilibrium temperature formulation

    if( .not. COMPUTESOLRAD )then
      ! *** MUST MAKE AT LEAST ONE CALL TO THIS TO INITIALIZE VARIABLES
      call SHORT_WAVE_RADIATION(CLOUDT(2),SRO,SRON)
    endif

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K)  &
    !$OMP                             PRIVATE(ET,CSHE,THICK,TFLUX,TEMP,ET0,CSHE0,WSHLTR0,PSHADE0)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)

      ! *** INITIALIZE DOMAIN EVAPOTRANSPIRATION (ET)
      L = LWET(LF)
      call EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)
      PSHADE_OLD(ND) = PSHADE(L)
      PSHADE0 = PSHADE(L)
      WSHLTR_OLD(ND) = WINDSTKA(L)
      WSHLTR0 = WINDSTKA(L)
      ET0 = ET
      CSHE0 = CSHE

      ! *** SURFACE HEAT EXCHANGE PROCESSES
      do LP = LF,LL
        L = LWET(LP)
        if( .not. ICECELL(L) )then
          if( NASER > 1 )then
            call EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)

          elseif( PSHADE_OLD(ND) /= PSHADE(L) .or. WSHLTR_OLD(ND) /= WINDSTKA(L) )then
            if( WINDSTKA(L) == WSHLTR0 .and. (PSHADE(L) == PSHADE0 .or. SOLSWRT(L) == 0.0) )then
              ! *** Use PREVIOUS CONDITIONS
              call SWAP(ET,ET0)
              call SWAP(CSHE,CSHE0)
              call SWAP(PSHADE0,PSHADE_OLD(ND))
              call SWAP(WSHLTR0,WSHLTR_OLD(ND))
            else
              ET0 = ET
              CSHE0 = CSHE
              PSHADE0 = PSHADE_OLD(ND)
              WSHLTR0 = WSHLTR_OLD(ND)
              call EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)
            endif
            PSHADE_OLD(ND) = PSHADE(L)
            WSHLTR_OLD(ND) = WINDSTKA(L)
          endif

          ! *** SURFACE HEAT FLUX
          THICK = HPK(L,KC)
          TFLUX = CSHE*(ET - TEM(L,KC))/THICK*DELT
          TEM(L,KC) = TEM(L,KC) + TFLUX
        endif  ! *** END OF ICE COVER CHECK FOR CELL
      enddo

      ! *** Distribute Solar Radiation Across Water Column
      if( LDAYLIGHT .and. KC > 1 )then
        ! *** NOW FINALIZE THE TEMPERATURE
        do LP = LF,LL
          L = LWET(LP)
          do K = KS,KSZ(L),-1
            ! *** RHO = 1000.0  Density (kg / m^3)
            ! *** CP  = 4179.0  Specific Heat (J / kg / degC)
            ! *** RHOWCPI = 1/RHO/CP = 0.2393E-6 (m^3 degC)/(W s)
            !C1 = DELT*DZIC(L,K)*RHOWCPI
            TEMP = TEM(L,K) + HPI(L)*RADNET(L,K)
            TEM(L,K) = TEMP
          enddo
        enddo

      endif    ! *** END OF ACTIVE SOLAR RADIATION

    enddo   ! *** END OF DOMAIN
    !$OMP END PARALLEL DO

  elseif( ISTOPT(2) == 4 )then                   ! *** Externally specified equilibrium temperature formulation
  
    ! *** This option uses "cloudt" in a specialized manor to contain the surface exchange coefficient
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)
      do LP = LF,LL
        L = LWET(LP)
        TMPKC     = DELT*DZIC(L,KC)
        TEM(L,KC) = TEM(L,KC) - TMPKC*CLOUDT(L)*HPI(L)*(TEM(L,KC) - TATMT(L)) ! *** CLOUDT(L) IS NOT FRACTION OF CLOUD COVER
      enddo
    enddo
    
  endif    ! *** End of surface heat calculation options

  ! *** Hardwire heat fluxes for GOTM testing
  if( iGOTM_Test > 0 )then
    do L = 2,LA
      C1 = DELT*RHOWCPI
      RADNET(L,KC) = heatflux
      TEMP = TEM(L,KC) + HPKI(L,KC)*C1*RADNET(L,KC)
      TEM(L,KC) = TEMP
    enddo
  endif

  ! *** Apply dry cell corrections
  if( ISDRY > 0 .and. LADRY > 0 )then
    do K = 1,KC
      do LP = 1,LADRY
        L = LDRY(LP)

        ! *** ZERO HEAT ARRAYS
        RADBOT(L,K) = 0.
        RADTOP(L,K) = 0.
        RADNET(L,K) = 0.
      enddo
    enddo
  endif

  ! *** CHECK NON-ICE "DRY" CELLS
  if( LAWET < LA-1 .and. TIMEDAY > TIMEDRY )then
    ! *** Only check every 6 hrs, i.e. TIMEDRY increments 0.25 of a day.  Incremented below.
    if( ISICE > 2 )then
      do L = 2,LA
        if( .not. LMASKDRY(L) .or. HP(L) < HDRY )then
          if( ICETHICK(L) == 0.0 )then
            if( TEM(L,KC) < -5 .or. TEM(L,KC) > 30. )then
              TEM(L,:) = max(TATMT(L),0.0)
            endif
          endif
        endif
      enddo
    endif
  endif

  ! *************************************************************************************************************
  ! *************************************************************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** Bed/water column interface heat flux
  if( NASER > 0 .and. ISTOPT(2) /= 0 )then

    ! *** Update bottom water column layer temperatures
    if( (HTBED1 + HTBED2) > 0.0 )then
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,UBED,VBED,USPD,TFLUX,THICK,TEMP)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)
        do LP = LF,LL
          L = LWET(LP)
          UBED = 0.5*( U(L,KSZU(L)) + U(LEC(L),KSZU(LEC(L))) )
          VBED = 0.5*( V(L,KSZV(L)) + V(LNC(L),KSZV(LNC(L))) )
          USPD = SQRT( UBED*UBED + VBED*VBED )
          TEMP = max(TEM(L,KSZ(L)), 0.1 )

          !        (nodim) (m/s)   (m/s)      (C)       (C)    (s)
          TFLUX = ( HTBED1*USPD + HTBED2 )*( TEMB(L) - TEMP )*DELT   ! *** C*m
                    
          ! *** Accumulate bottom heat FLUXTB (to avoid truncation error, apply minimum FLUXTB)
          FLUXTB(L) = FLUXTB(L) + TFLUX
          THICK = HPK(L,KSZ(L))

          if( ABS(FLUXTB(L))/THICK > 0.001 )then
            TEM(L,KSZ(L)) = TEM(L,KSZ(L)) + FLUXTB(L)/THICK       ! *** Update water column for bed heat exchange
            if( FLUXTB(L) < 0.0 )then
              TEM(L,KSZ(L)) = max(TEM(L,KSZ(L)),TEMB(L))
            else
              TEM(L,KSZ(L)) = min(TEM(L,KSZ(L)),TEMB(L))
            endif
            if( TEMBO <= 0.0 ) FLUXTB(L) = 0.0                    ! *** Zero accumulator here, if not updating the bed T's
          endif
        enddo
      enddo   ! *** END OF DOMAIN
      !$OMP END DO
    endif

    ! *** Update bed temperatures, if needed
    if( TEMBO > 0.0 )then
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,FLUXQB,TEMP,THICK)
      do ND = 1,NDM
        LF = (ND-1)*LDMWET+1
        LL = min(LF+LDMWET-1,LAWET)

        ! *** Accumulate solar radiation on sediment surface
        if( LDAYLIGHT )then
          do LP = LF,LL
            L = LWET(LP)
            RADBOTT(L) = RADBOTT(L) + RADBOT(L,KSZ(L))*DELT
          enddo
        endif

        ! *** Updated bed temperatures, if required
        do LP = LF,LL
          L = LWET(LP)

          ! *** Calculate longwave heat emission from the bed
          if( HP(L) < 0.5 .or. (HP(L) < 3.0 .and. RADKE(L,KSZ(L)) < 1.0) )then
            if( .not. ICECELL(L) )then
              ! *** 4.43E-14 = 1/rhob * 1/cpb * 5.67e-8, where rhob = 1600 kg/m3 and cpb = 800 J/kg/C
              TEMP = max(TEM(L,KSZ(L)), 0.1 )
              FLUXQB = 4.43E-14*((TEMB(L) + 273.)**4 - (TEMP + 273.)**4)*DELT

              ! *** Update bed temperature
              TEMB(L) = TEMB(L) - FLUXQB/TBEDTHK(L)
            endif
          endif
          
          THICK = HPK(L,KSZ(L))
          if( ABS(FLUXTB(L))/THICK > 0.001 )then
            ! ***       C         (     C*m        )     C*m          1/m
            TEMB(L) = TEMB(L) + ( RHOWCPI*RADBOTT(L) - FLUXTB(L) )/TBEDTHK(L)        ! *** Update bed temperature
            FLUXTB(L)  = 0.0
            RADBOTT(L) = 0.0
            TEMB(L) = max(MIN(TEMB(L), TATMT(1)+20.), 0.)                            ! *** Prevent shallow cells from excessive heat
          endif
        enddo
      enddo   ! *** END OF DOMAIN
      !$OMP END DO
            
      ! *** Update bed temperature for dried cells to get reasonable value when they are wet again
      ! *** Here we only account for the conductive heat loss and longwave heat emission to the air
      ! *** Heat conductivity of sediment and air is taken as a half of heat conductivity of sediment and water
      ! *** 4.43E-14 = 1/rhob * 1/cpb * 5.67e-8, where rhob = 1600 kg/m3 and cpb = 800 J/kg/C
      !$OMP SINGLE
      do L = 2,LA
        if( HP(L) < HDRY )then
          TEMB(L) = TEMB(L) - 0.5*HTBED2*( TEMB(L) - TATMT(L))*DELT/TBEDTHK(L)  &
                   - 4.43E-14*((TEMB(L) + 273.)**4 - (TATMT(L) + 273.)**4)*DELT/TBEDTHK(L)
          TEMB(L) = max(MIN(TEMB(L), TATMT(1)+20.), 0.) 
        endif
      enddo
      !$OMP END SINGLE
    endif
  endif

  ! *** Limit surface temperatures based on ISICE options **********************************

  ! *** Get the minimum surface temperature
  if( ISICE < 3 )then
    !$OMP DO PRIVATE(ND,LP,L)
    do ND = 1,NDM
      do LP = 1,LLWET(KC,ND)
        L = LKWET(LP,KC,ND)
        TMIND(ND) = min(TEM(L,KC),TMIND(ND))
      enddo
    enddo   ! *** END OF DOMAIN
    !$OMP END DO

    !$OMP BARRIER
    !$OMP SINGLE
    TMIN = MINVAL(TMIND)
    !$OMP END SINGLE

  endif

  ! *** LIMIT WATER TEMPERATURES WHEN ICE FREEZE/MELT IS NOT COUPLED TO THE HEAT SUB-MODEL
  if( ISICE < 3 )then
    if( TMIN < 1.0 )then
      if( ISICE == 0 )then
        ! *** LIMIT WATER TEMPERATURES TO > FREEZING TEMPERATURE IF ICE IS NOT SIMULATED
        K = KC
        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM
          LF = (ND-1)*LDMWET+1
          LL = min(LF+LDMWET-1,LAWET)

          do LP = LF,LL
            L = LWET(LP)
            ! *** LIMIT THE WATER TEMPERATURE
            if( ISTRAN(1) > 0 )then
              if( TEM(L,KC) < -1.3 )then
                TEM(L,KC) = -1.3*(SAL(L,KC)/35.)
              endif
            else
              if( TEM(L,KC) < 0.001 )then
                ! *** LIMIT THE WATER TEMPERATURE
                TEM(L,KC) = 0.001
              endif
            endif
          enddo
        enddo   ! *** END OF DOMAIN
        !$OMP END DO

      elseif( ISICE == 1 .or. ISICE == 2 )then
        ! *** FORCE A LIMITATION OF TEMPERATURE TO ICE TEMP FOR USER SPECIFIED ICE COVER
        K = KC
        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        do ND = 1,NDM
          LF = (ND-1)*LDMWET+1
          LL = min(LF+LDMWET-1,LAWET)

          do LP = LF,LL
            L = LWET(LP)

            ! *** LIMIT THE WATER TEMPERATURE
            if( TEM(L,K) < TEMPICE )then
              TEM(L,K) = TEMPICE
            else
              TEM(L,K) = ICECOVER(L)*TEMPICE + (1.-ICECOVER(L))*TEM(L,K)
            endif
          enddo
        enddo   ! *** END OF DOMAIN
        !$OMP END DO
      endif
    endif
  endif   ! *** END OF TMIN < 1.0 BLOCK
  !$OMP END PARALLEL

  ! *** Check shallow cells temperature for excessive heating
  if( LAWET < LA-1 .and. TIMEDAY > TIMEDRY )then
    TIMEDRY = TIMEDRY + 0.25   ! *** Only check every 6 hrs
    if( TEMBO  > 0.0 )then
      do L = 2,LA
        if( .not. LMASKDRY(L) .or. HP(L) < HDRY )then
          TEMB(L) = max(MIN(TEMB(L), TATMT(1)+20.), 0.)                            ! *** Prevent shallow cells from excessive heat
          !TEMB(L) =  max(TATMT(L),0.0)   ! *** Limit bed temperature to air temp for dry cells
        endif
      enddo
    endif
  endif

  return

END SUBROUTINE CALHEAT

SUBROUTINE ICECOMP(DELTD2, ICESTEP)

  ! *************************************************************************************************
  ! *** HEAT COUPLED ICE COMPUTATIONS WITH MASS BALANCE IN FREEZE/MELT
  ! ***
  ! *** ISICE = 0 - ICE NOT INCLUDED
  ! *** ISICE = 1 - USER SPECIFIED TEMPORALLY/SPATIALLY VARYING ICE COVER
  ! *** ISICE = 2 - USER SPECIFIED CTIME VARYING COMPLETE ICE COVER USING BINARY ON/OFF
  ! *** ISICE = 3 - FULLY HEAT COUPLED
  ! *** ISICE = 4 - FULLY HEAT COUPLED WITH FRAZIL ICE TRANSPORT

  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------------------!
  ! 2015-05           PAUL M. CRAIG    ADDED FULLY HEAT COUPLED ICE FORMATION/MELT
  !                   DANG H CHUNG

  real, intent(IN)    :: DELTD2,ICESTEP

  integer :: L, K, ND, LF, LL, IFILE, ICEITER, ITERI
  integer,save :: ITERMAX, ITERMAX2

  real       :: THICKMIN, MINICE, TMPVAL, ET, CSHE
  real       :: WTEMP, RHOA
  real       :: DEL, RANLW, RN, TF, TFS, RB, RC, RE, RT, TICEBOT, DELTICE
  real       :: RATIODENS, DICETHI, DICETHW, HICE

  real,save,allocatable,dimension(:) :: ICETEM_OLD

  !REAL,external :: FUNDEN

  character*9 :: FMODE

  if( .not. allocated(ICETEM_OLD) )then
    allocate(ICETEM_OLD(LCM))
    ICETEM_OLD = -9999.

    ! *** RHOI = Density of Ice (kg / m^3)
    ! *** CP  = 4179.0  Specific Heat (J / kg / degC)
    ! *** RHOWICPI = 1/RHOI/CP
    RHOICP   = RHOI*CP
    RHOICPI  = 1./RHOICP
    RHOILHFI = 1./RHOI/LHF
    CREFLI   = 1./(1.-REFL)*(1.0-ALBEDOI)

    ITERMAX = 500
    ITERMAX2= 0
  endif

  IFILE = -1

  if( ISTL == 2 )then
    if( ISDYNSTP == 0 )then
      DELTICE = DT
    else
      DELTICE = DTDYN
    endif
  else
    DELTICE = DT  !*** DKT - Use the same time increment as in 2TL
  endif

  RATIODENS = RHOI/999.8426                                 ! *** 999.82 kg/m3 is water density at 20 degC
  LFRAZIL = .FALSE.

  ! *** NEED TO LOOP OVER ALL THE CELLS BECAUSE SOME ICE COVERED CELLS MAY BE "DRY"
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, L, K, ICEITER, ITERI)                     &
  !$OMP                             PRIVATE(DICETHW, DICETHI, MINICE, TF, TFS, TMPVAL)  &
  !$OMP                             PRIVATE(WTEMP, RHOA, THICKMIN, ET, CSHE, TICEBOT)            &
  !$OMP                             PRIVATE(RANLW, RT, RB, RC, RE, RN, DEL, HICE)
  do ND = 1,NDM
    LF = 2+(ND-1)*LDM
    LL = min(LF+LDM-1,LA)
    do L = LF,LL
      DICETHW   = 0.0   ! *** CHANGE IN ICE THICKNESS DUE TO WATER TEMPERATURES FOR OPEN WATER (NO ICE)/ICE WATER INTERFACE (WITH ICE COVER)  (M)
      DICETHI   = 0.0   ! *** CHANGE IN ICE THICKNESS DUE TO ICE TEMPERATURES, I.E. BOTTOM GROWTH (M)

      ICETHICK1(L) = ICETHICK(L)
      ! *** HANDLE SHALLOW CELLS
      MINICE = MINICETHICK
      if( HP(L) <= HDRY*1.2 .and. ISDRY > 0 )then
        MINICE = MINICETHICK*DZC(L,KC)
      endif

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
      TFS = TF - 0.01  ! *** SUPERCOOLED STATE

      ! *** DETERMINE ICE FORMATION AND HEAT FLUX DUE TO ICE FORMATION.  ALLOW ONLY FOR "WET" CELLS
      if( ISICE == 4 .and. .not. ICECELL(L) )then
        ! *** COUPLED HEAT/ICE MDOEL WITH FRAZIL ICE TRANSPORT
        do K = KSZ(L),KC
          if( TEM(L,K) <= TFS )then
            if( HP(L)  >= HDRYICE )then
              ! *** FRAZIL ICE GROWTH
              TMPVAL = -(TF-TEM(L,K))*RHOW(L,K)*CP*HP(L)*DZC(L,K)
              TMPVAL = -TMPVAL*RHOILHFI
              FRAZILICE(L,K) = FRAZILICE(L,K) + TMPVAL

              ! *** REMOVE THE VOLUME OF WATER FROM THE WATER COLUMN
              ICEVOL(L) = ICEVOL(L) - TMPVAL*RATIODENS*DXYP(L)
              LFRAZIL = .TRUE.
            endif

            ! *** ICE FORMATION LATENT HEAT OF FUSION ADDS HEAT
            if( IS2TL == 0 )then
              if( ISTL == 3 )then
                TEM1(L,K) = TF
              endif
            else
              TEM1(L,K) = TF
              TEM(L,K)  = TF
            endif

            ! ***  AMBIENT DENSITY
            if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
              RHOW(L,K) = FUNDEN(SAL(L,KC),SEDT(L,KC),TF)
            else
              RHOW(L,K) = FUNDEN(SAL(L,KC),0.,TF)
            endif

          elseif( TEM(L,K) > 0.0 )then
            ! *** CHECK FOR MELTING FRAZIL ICE
            if( FRAZILICE(L,K)  > 0.0 )then
              TMPVAL     = -DELTICE*HWI*(TEM(L,K)-TF)*RHOILHFI
              FRAZILICE(L,K) = FRAZILICE(L,K) + TMPVAL
              FRAZILICE(L,K) = max(FRAZILICE(L,K),0.)
              if( FRAZILICE(L,K)  > 0.0 )LFRAZIL = .TRUE.

              ! *** ICE MELT ADDS WATER AT T = 0.0 DEGC
              ICEVOL(L) = ICEVOL(L) - TMPVAL*RATIODENS*DXYP(L)
            endif

            ! *** ACCOUNT FOR ICE "COVER" MELT
            if( ICETHICK(L)  > 0.0 .and. K == KC )then
              TMPVAL     = -DELTICE*HWI*(TEM(L,K)-TF)*RHOILHFI
              ICETHICK(L) = ICETHICK(L) + TMPVAL
              ICETHICK(L) = max(ICETHICK(L),0.)
              if( ICETHICK(L)  > 0.0 )LFRAZIL = .TRUE.

              ! *** ICE MELT ADDS WATER AT T = 0.0 DEGC
              ICEVOL(L) = ICEVOL(L) - TMPVAL*RATIODENS*DXYP(L)
            endif

          endif

        enddo

      elseif( ICETHICK(L) < 1E-6 )then   ! *** ISICE = 3

        ! *** COUPLED HEAT MODEL WITHOUT FRAZIL ICE TRANSPORT
        if( TEM(L,KC) < TFS )then
          ! *** GET ANY "FREEZING" WATER THAT MAY HAVE SETTLED/TRANSPORTED IN PRIOR TIMESTEPS WHERE ICE THICKNESS < MINICETHICK
          WTEMP  = 0.
          do K = KSZ(L),KC
            if( TEM(L,K) <= TF )then
              WTEMP  = WTEMP  + (TF-TEM(L,K))*DZC(L,K)
            endif
          enddo
          WTEMP = WTEMP*DZIC(L,KC)

          ! ***  AMBIENT DENSITY
          if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
            RHOA = FUNDEN(SAL(L,KC),SEDT(L,KC),TF)
          else
            RHOA = FUNDEN(SAL(L,KC),0.,TF)
          endif

          ! *** INITIAL ICE FORMATION FOR OPEN WATER
          THICKMIN = HP(L)*DZC(L,KC)
          if( THICKMIN < HDRY )then
            THICKMIN = 0.
            do K = KSZ(L),KC
              if( TEM(L,K) <= TF )then
                THICKMIN  = THICKMIN  + HP(L)*DZC(L,K)
              endif
            enddo
          endif
          TMPVAL  = WTEMP*RHOA*CP*THICKMIN
          DICETHW = TMPVAL*RHOILHFI

          ! *** CHECK IF SUFFICIENT WATER EXISTS
          if( DICETHW > 0.1*HP(L) .and. HP(L) < HDRYICE )then
            DICETHW = 0.0
          endif
          
          ! *** UPDATE ALL THE TEMPERATURES
          do K = KSZ(L),KC
            if( TEM(L,K) <= TF )then
              TEM1(L,K) = TF
              TEM(L,K)  = TF                            ! *** ICE FORMATION LATENT HEAT OF FUSION ADDS HEAT

              ! ***  AMBIENT DENSITY
              if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
                RHOW(L,K) = FUNDEN(SAL(L,K), SEDT(L,K),TF)
              else
                RHOW(L,K) = FUNDEN(SAL(L,K),0.,TF)
              endif
            endif
          enddo

          ICEVOL(L) = ICEVOL(L) - DICETHW*RATIODENS*DXYP(L)    ! *** REMOVE THE VOLUME OF WATER FROM THE WATER COLUMN

        endif       ! *** END BLOCK TEM(L,KC) < TF
      endif         ! *** ISICE OPTION

      ! *** ICE BALANCE, IF ICE COVER ALREADY EXISTS FOR THE CELL. ALLOW ICE MELT EVEN FOR "DRY" CELLS
      if( ICECELL(L) )then
        ! ******************************************************************************
        ! *** DETERMINE ICE TEMPERATURE AT TOP OF ICE COVER (Ts) AFTER REACHING
        ! *** QUASI STEADY STATE CONDITIONS
        if( DELTD2 >= ICESTEP )then
          if( .not. LMASKDRY(L) ) CALL SET_LIGHT(L)   ! *** UPDATE DRY CELL SOLAR RADIATION CONDITIONS

          if( ICETEM_OLD(L) == -9999. )then
            ICETEMP(L) = TATMT(L)
          else
            ICETEMP(L) = ICETEM_OLD(L) + SIGN(0.005,TATMT(L))
          endif

          ! *** INCIDENT LONG WAVE RADIATION
          if( TATMT(L)  >= 5.0 )then
            RANLW = 5.31E-13*(273.15+TATMT(L))**6*(1.0+0.17*CLOUDT(L)**2)*0.97
          else
            RANLW = 5.62E-8*(273.15+TATMT(L))**4*(1.-0.261*EXP(-7.77E-4*TATMT(L)**2))*(1.0+0.17*CLOUDT(L)**2)*0.97
          endif

          ! *** RT = TOTAL RADIATION AT TOP OF ICE, SHORTWAVE PLUS LONGWAVE
          RT = SOLSWRT(L)*CREFLI + RANLW

          ICEITER = 0
          do ITERI = 1, ITERMAX
            ICEITER = ICEITER+1
            call SURFACE_TERMS (ICETEMP(L),L,RB,RC,RE)          ! *** OUT: RB,RC,RE
            RN = RT-RB-RE-RC                                    ! *** W/M2
            DEL = RN + ICEK*(TF-ICETEMP(L))/ICETHICK(L)         ! *** W/M2
            TMPVAL = DEL*ICETHICK(L)/10.
            if( ABS(TMPVAL) > 1. )then
              TMPVAL = SIGN(1.0,TMPVAL)
            endif
            ICETEMP(L) = ICETEMP(L) + TMPVAL
            if( ABS(TMPVAL) < 0.01 ) EXIT
          enddo

          if( ICEITER > ITERMAX2 )then
            ! PROVIDE SOME INFORMATION ON THE ITERATIONS

            if( IFILE == -1 )then
              IFILE = mpi_efdc_out_unit
              open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
              open(mpi_error_unit,FILE = OUTDIR//mpi_error_file,POSITION = 'APPEND')
            endif
            ITERMAX2 = ICEITER

            write(mpi_efdc_out_unit, "(A,2I5,2F8.3,F10.4)" )'ICE TEMPERATURE: NEW MAX ITERATION (MAX L ICETHK ICET TIME): ',ICEITER,L,ICETHICK(L),ICETEMP(L),TIMEDAY

            if( ICEITER > ITERMAX )then
              PRINT '(A,2I5)','ERROR: ICE TEMPERATURE ITERATIONS EXCEEDED MAX  ALLOWABLE: (MAX L ICETHK ICET TIME): ',ICEITER,L
              write(mpi_error_unit,"(A,2I5,2F8.3,F10.4)")'ERROR: ICE TEMPERATURE ITERATIONS EXCEEDED MAX  ALLOWABLE: (MAX L ICETHK ICET TIME): ',ICEITER,L,ICETHICK(L),ICETEMP(L),TIMEDAY

              if( ICETEM_OLD(L) == -9999. )then
                ICETEMP(L) = min(0.5*TATMT(L),0.0)
              else
                ICETEMP(L) = ICETEM_OLD(L)  ! *** Use LAST GOOD ICE TEMPERATURE
              endif
            endif
          endif
          ICETEM_OLD(L) = ICETEMP(L)

          ! ******************************************************************************
          ! *** SOLAR RADIATION ATTENUATION THROUGH ICE IS HANDLED IN SET_LIGHT FUNCTION

          ! *** GET INTERACTION AT THE ICE INTERFACE DUE TO ICE TEMPERATURE
          if( ICETEMP(L) > 0.0 )then
            ! *** ICE MELT DUE TO ICE WARMING FROM AIR/ICE INTERFACE.  ICE IS PURE WATER SO, TF = 0.0
            ICETEMP(L) = min(ICETEMP(L),0.001/ICETHICK(L))
            DICETHI = -CP*ICETEMP(L)*ICETHICK(L)*MELTFACTOR/LHF*ICESTEP/DELTICE

            ! *** ICE MELT ADDS WATER AT T = 0.0 DEGC
            ICEVOL(L) = ICEVOL(L) - DICETHI*RATIODENS*DXYP(L)
          endif
        endif

        if( ICETEMP(L) <= 0.0 )then
          ! *** ICE GROWTH AT ICE/WATER INTERFACE (ICE BOTTOM)
          if( LMASKDRY(L) .and. ICETHICK(L) <= ICETHMX )then
            TICEBOT = min(ICETEMP(L) + TEM(L,KC),0.0)       ! *** ADJUST ICE TEMPERATURE FOR WATER CONTACT
            HICE    = ICEK*(TF-TICEBOT)/ICETHICK(L)         ! *** W/M2
            DICETHI = DELTICE*HICE*RHOILHFI

            ! *** ICE FORMATION LATENT HEAT OF FUSION CONSUMES THE HEAT
            ICEVOL(L) = ICEVOL(L) - DICETHI*RATIODENS*DXYP(L)

            if( TEM(L,KC) > 1.0 )then
              TEM(L,KC)  = TEM(L,KC)*MAX(1.-DICETHI*RATIODENS*HPKI(L,KC),0.0)
              TEM1(L,KC) = TEM(L,KC)
            endif
          endif
        else
          ICETEMP(L) =  0.0
        endif
      endif

      IF( ICETHICK(L)  >= 1.E-6 )THEN
        ! *** ICE MELT/FREEZE FROM WATER-ICE INTERFACE
        if( TEM(L,KC) > 0.0 )then
          ! *** MELT
          HICE   = -HWI*TEM(L,KC)*RHOILHFI
          DICETHW = DELTICE*HICE

          ! *** ICE MELT ADDS WATER AT T = 0.0 DEGC
          ICEVOL(L) = ICEVOL(L) - DICETHW*RATIODENS*DXYP(L)
        elseif( TEM(L,KC) < TFS )then
          ! *** FREEZE
          WTEMP = TF - TEM(L,KC)
          if( HP(L) < HDRY )then
            THICKMIN = HP(L)
          else
            THICKMIN = HP(L)*DZC(L,KC)
          endif
          TMPVAL  = WTEMP*CP*THICKMIN*999.8426
          DICETHW = TMPVAL*RHOILHFI

          ICEVOL(L) = ICEVOL(L) - DICETHW*RATIODENS*DXYP(L)    ! *** REMOVE THE VOLUME OF WATER FROM THE WATER COLUMN
          TEM1(L,KC) = TF
          TEM(L,KC)  = TF
        endif
      endif   ! *** END BLOCK ICECELL(L)

      ! *** TOTAL ICE THICKNESS AT THE TOP OF THE WATER COLUMN (M)
      ICETHICK(L) = ICETHICK(L) + (DICETHI + DICETHW)

      if( ICETHICK(L) >= MINICE )then
        ICECELL(L) = .TRUE.
        ICECOVER(L) = 1.
        LFRAZIL = .TRUE.
      else
        ICETEM_OLD(L) = -9999.
        if( ISICE == 3 )then
          if( ICETHICK(L)  >  0.0 .and. ICETHICK(L) < 1.E-6 )then
            call EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)
            if( ET  >  1.0 )then
              if( DICETHW  > 0.0 ) ICETHICK(L) = max(ICETHICK(L) - DICETHW,0.0)      ! *** PREVENT A FREEZE AND A THAW IN THE SAME TIME STEP FOR THE SAME CELL
              ! *** COMPLETE MELT
              TMPVAL = ICETHICK(L)*RATIODENS*DXYP(L)
              ICEVOL(L) = ICEVOL(L) + TMPVAL
              ICETHICK(L) = 0.
            endif
          endif
              
          ICECOVER(L) = 0.
          ! *** ACCUMULATE ANY REMAINING ICE VOLUME AND ADD IT BACK TO THE WATER COLUMN
          if( ICEVOL(L)  /= 0.0 ) LFRAZIL = .TRUE.   ! ACCOUNT FOR FINAL "MELT"
        else   ! *** ISICE = 4
          ICECOVER(L) = ICETHICK(L)/MINICETHICK
          if( FRAZILICE(L,KC)  > 0.0 .or. ICETHICK(L) > 0. ) LFRAZIL = .TRUE.   ! ALLOW FRAZIL ICE TRANSPORT EVEN IF NO ICE COVER OR NEW FORMATION
        endif

        ! *** CHECK FOR SPECIAL CASES
        if( ICECELL(L) )then
          ! *** COMPLETE MELT
          if( ISDRY > 0 .and. HP(L) < 3.*HDRY )then
            ! *** COMPUTE EQUILIBRIUM TEMPERATURE
            call EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)
            ET = max(ET,0.1)
            if( ET+2. < TEM(L,KC) )then
              do K = KSZ(L),KC
                TEM(L,K) = ET
              enddo
            endif
            ICETEMP(L) =  0.0       ! *** Ensure ice temperature is zeroed
          endif
        endif

        ICECELL(L) = .FALSE.
        ICETEMP(L) =  0.0       ! *** Ensure ice temperature is zeroed
      endif
            
      ! *** Update water depth for 3TL 
      if( IS2TL == 0 )then
         TMPVAL = ICEVOL(L)/DXYP(L)
         if( TMPVAL > 0.0 .or. HP(L) + TMPVAL > 0.1*HDRY )then
           HP(L) = HP(L) + TMPVAL
           H1P(L) = H1P(L) + TMPVAL
           do K = KSZ(L),KC
             H1PK(L,K) = H1P(L)*DZC(L,K)
           enddo
         endif
         ICEVOL(L) = 0.
      endif
    enddo
  enddo   ! *** END OF DOMAIN
  !$OMP END PARALLEL DO

  if( IFILE == mpi_efdc_out_unit )then
    close(mpi_efdc_out_unit)
    close(mpi_error_unit)
  endif

  return

END SUBROUTINE ICECOMP

SUBROUTINE SHORT_WAVE_RADIATION(CLD,SRO,SRON)

  ! ************************************************************************
  ! ***                 S H O R T  W A V E  R A D I A T I O N              **
  ! ***                      FROM CE-QUAL-W2 (VER 3.1)                     **
  ! ************************************************************************

  ! ******* Type declaration
  real, intent(IN)    ::  CLD
  real, intent(OUT)   :: SRO, SRON
  integer   :: IDAY
  real      :: JDAY
  real      :: STANDARD, THOUR, PMC1, EQTNEW, H, DECL, SINAL, A0, CLD10

  ! ******* Input Conversions
  CLD10 = CLD*10.                ! *** CONVERT FROM 0-1 (EFDC) TO 0-10 (W2)

  ! ******* Shortwave Radiation
  STANDARD = 15.0*INT(DS_LONG/15.0)

  ! *** Day of the Year
  THOUR    = (TIMEDAY-INT(TIMEDAY))*24.0
  IDAY     =  TIMEDAY-INT(TIMEDAY/365.)*365.
  IDAY     =  IDAY+INT(INT(TIMEDAY/365.)/4.)
  JDAY     = real(IDAY)
  PMC1     = (2.*PI*(JDAY-1.))/365.
  EQTNEW   =  0.170*SIN(4.*PI*(JDAY-80.)/373.)-0.129*SIN(2.*PI*(JDAY-8.)/355.)
  H     = 0.2618*(THOUR-(DS_LONG-STANDARD)*0.066667+EQTNEW-12.0)
  DECL =  0.006918-0.399912*COS(PMC1)  +0.070257*SIN(PMC1)-0.006758*COS(2*PMC1)+0.000907*SIN(2*PMC1)-0.002697*COS(3*PMC1)+0.001480*SIN(3*PMC1)
  SINAL = SIN(DS_LAT*.017453)*SIN(DECL)+COS(DS_LAT*.017453)*COS(DECL)*COS(H)
  A0    = 57.2958*ASIN(SINAL)

  if( A0 > 0.0 )then
    SRO  = 2.044*A0+0.1296*A0**2-1.941E-3*A0**3+7.591E-6*A0**4
    SRO  = (1.0-0.0065*CLD10**2)*SRO*24.0
    SRO  = SRO*BTU_FT2_DAY_TO_W_M2
    SRON = SRO*(1.-REFL)                    ! *** Adjust for surface reflection
  else
    SRO  = 0.0
    SRON = 0.0
  endif

  return

END SUBROUTINE SHORT_WAVE_RADIATION

SUBROUTINE EQUILIBRIUM_TEMPERATURE(WSPD, TD, TAIR, SRON, ET, CSHE)
  ! *** *********************************************************************
  !**             E Q U I L I B R I U M  T E M P E R A T U R E           **
  !**                      FROM CE-QUAL-W2 (VER 3.1)                     **
  ! *** *********************************************************************

  ! ******* Type declaration
  real, intent(IN)    :: WSPD, TD, TAIR, SRON
  real, intent(OUT)   :: ET, CSHE
  real      :: TDEW_F, TAIR_F, WIND_2M
  real      :: WIND_MPH
  real      :: SRO_BR, TSTAR, BETA, FW, RA, ETP

  integer :: J, NAL

  ! ******* Input Conversions

  ! ******* British units
  TDEW_F   = TD*1.8+32.0
  TAIR_F   = TAIR*1.8+32.0
  WIND_MPH = WSPD*MPS_TO_MPH
  WIND_2M  = WIND_MPH         ! *** EE7.3 AND LATER APPLY CORRECTION TO ALL WINDS
  WIND_2M = max(WIND_2M, 1.)  ! PMC - APPLY A MINUMUM TO ALLOW SURFACE HEAT EXCHANGE EVEN WITH CALM CONDITIONS

  ! *** SRON Should already be adjusted for Shading & Reflection
  SRO_BR   = SRON*W_M2_TO_BTU_FT2_DAY

  ! ******* Equilibrium temperature and heat exchange coefficient

  ET    = TDEW_F
  TSTAR = (ET+TDEW_F)*0.5
  BETA  = 0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
  FW    = W_M2_TO_BTU_FT2_DAY*AFW+BCONV*BFW*WIND_2M**CFW
  CSHE  = 15.7+(0.26+BETA)*FW
  RA    = 3.1872E-08*(TAIR_F+459.67)**4
  ETP   = (SRO_BR+RA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE*(0.26+BETA))
  J     = 0
  do while (ABS(ETP-ET) > 0.05 .and. J < 10)
    ET    = ETP
    TSTAR = (ET+TDEW_F)*0.5
    BETA  = 0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
    CSHE  = 15.7+(0.26+BETA)*FW
    ETP   = (SRO_BR+RA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE*(0.26+BETA))
    J     = J+1
  enddo

  ! ******* SI units

  ! *** RHOWCPI = 1/RHO/CP = 0.2393E-6
  ET   = (ET-32.0)*5.0/9.0                 ! *** DEG C
  CSHE = CSHE*FLUX_BR_TO_FLUX_SI*RHOWCPI   ! *** M/S

END SUBROUTINE EQUILIBRIUM_TEMPERATURE

REAL FUNCTION GET_EXT_COEFF(L,K)
  ! *** COMPUTE THE LIGHT EXTINCTION COEFFICIENT
  
  ! *** RADKE  - EXTINCTION COEFFICIENT AT THE MIDPOINT OF THE LAYER

  integer, intent(IN) :: L, K
  real    :: WCKESS, CHLKE, TSS
  integer :: NAL

  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    TSS = SEDT(L,K) + SNDT(L,K)
  else
    TSS = 0.
  endif

  if( ISTRAN(8) > 0 )then
    ! *** If using WQ then use the WQ Coefficients
    WCKESS = WQKEB(IWQZMAP(L,K))                                     ! *** BACKGROUND EXTINCTION (BY WQ ZONE)
    WCKESS = WCKESS + WQKETSS*TSS                                    ! *** INORGANIC SOLIDS COMPONENT
    WCKESS = WCKESS + WQKEPOC*(WQV(L,K,IROC) + WQV(L,K,ILOC))        ! *** PARTICULATE ORGANIC MATTER
    WCKESS = WCKESS + WQKEDOM*WQV(L,K,IDOC)                          ! *** DISSOLVED ORGANIC MATTER

    ! *** Light extinction due to shading by stems and leaves (EEMS10.4)
    do NAL = 1,NALGAE
      if( .not. ALGAES(NAL).ISMOBILE )then
        if( K  >= LAYERBOT(NAL,L) .and. K <= LAYERTOP(NAL,L) )then
          WCKESS = WCKESS + ALGAES(NAL).WQKEMAC*WQV(L,K,19+NAL)      ! *** MACROPHYTE MASS G/M3 TO 1/M
        endif
      endif
    enddo

    if( WQKECHL < 0.0 )then                                          ! *** CHLOROPHYLL
      ! *** Compute Extinction Factor as a fn(Chla)
      CHLKE = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
    else
      CHLKE = WQKECHL*WQCHL(L,K)**WQKECHLE
    endif
    WCKESS = WCKESS + CHLKE

  else
    WCKESS = SVKEBACK(L) + WQKETSS*TSS
  endif
  RADKE(L,K) = WCKESS
  GET_EXT_COEFF = WCKESS

END FUNCTION GET_EXT_COEFF

SUBROUTINE SET_LIGHT(L)
  ! *** SETS THE SOLAR RADIATION FOR THE TOP OF THE WATER COLUMN

  ! *** RADTOP - SOLAR RADIATION AT THE TOP OF THE LAYER (W/M2)
  ! *** RADNET - NET SOLAR RADIATION ABSORBED IN THE LAYER (W/M2)
  ! *** RADBOT - SOLAR RADIATION AT THE BOTTOM OF THE LAYER (W/M2)

  integer, intent(IN) :: L
  integer :: K
  real    :: WCKESS, BOT, FRACLAYER, SOLBOT

  ! *** INITIALIZE CELL

  ! *** SOLSWRT has been adjusted before here to include reflectance and ice (if any)
  if( ICECELL(L) )then
    ! *** Adjust incident solar radiation at the top of the water due to ice
    SOLBOT = EXP(-GAMMAI*ICETHICK(L))
    BETAI = 1.0 - SOLBOT
    REFICE  = (1. - ALBEDOI)/(1.-REFL)*SOLBOT
    RADTOP(L,KC) = SOLSWRT(L)*REFICE
  else
    RADTOP(L,KC) = SOLSWRT(L)
  endif

  !IF( ISTRAN(8) > 0 .and. IASWRAD == 3 )then
  !  ! *** Set extinction for water quality for cases without surface heat exchange (i.e. test cases)
  !  do K = KC,KSZ(L),-1
  !    WCKESS = GET_EXT_COEFF(L,K)
  !  enddo
  !ENDIF

END SUBROUTINE

SUBROUTINE SURFACE_TERMS(TSUR,L,RB,RC,RE)
  ! *** CALCULATION OF HEAT EXCHANGE TERMS ON WATER SURFACE
  ! *** OUTPUT:
  ! *** FW,RE,RB,RC

  integer,intent(IN)    :: L
  real,   intent(IN)    :: TSUR
  real,   intent(INOUT) :: RB,RC,RE
  real :: VPA, VPS, TAIRV, DTV, DTVL, BOWEN_CONSTANT, FW
  DATA BOWEN_CONSTANT /0.47/

  ! *** ATMOSPHERIC CONDITIONS: VAPOR PRESSURE (mmHg)
  VPA = EXP(2.3026*(7.5*TDEWT(L)/(TDEWT(L)+237.3)+0.6609))

  ! *** SATURATION VAPOR PRESSURE AT WATER TEMPERATURE
  if( TSUR < 0.0 )then
    VPS = EXP(2.3026*(9.5*TSUR/(TSUR+265.5)+0.6609))
  else
    VPS = EXP(2.3026*(7.5*TSUR/(TSUR+237.3)+0.6609))
  endif

  ! *** EVAPORATIVE WIND SPEED FUNCTION, W/M2/mmHG
  if( ISRHEVAP == 1 )then
    TAIRV = (TATMT(L)+273.0)/(1.0-0.378*VPA/760.0)    ! *** Tav
    DTV   = (TSUR+273.0)/(1.0-0.378*VPS/760.0)-TAIRV  ! *** Tsv-Tav
    DTVL  =  0.0084*WINDST(L)**3
    if( DTV < DTVL) DTV = DTVL
    FW = (3.59*DTV**0.3333+4.26*WINDST(L))           ! *** Rayan-Harleman,1974
  else
    FW = AFWI+BFWI*WINDST(L)**CFWI
  endif

  ! *** EVAPORATIVE HEAT LOSS, (W/M2)
  RE = FW*(VPS-VPA)

  ! *** HEAT CONDUCTION, (W/M2)
  RC = FW*BOWEN_CONSTANT*(TSUR-TATMT(L))

  ! *** BACK RADIATION FROM WATER SURFACE, (W/M2)
  RB = 5.51E-8*(TSUR+273.15)**4

END SUBROUTINE

SUBROUTINE SWAP (X1,X2)
  ! *** FUNCTION TO SWAP THE VALUES OF TWO VARIABLES
  real, intent(INOUT) :: X1,X2
  real                :: X
  X = X1
  X1 = X2
  X2 = X
  return
END SUBROUTINE

END MODULE
