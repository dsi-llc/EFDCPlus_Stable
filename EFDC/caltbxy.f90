! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details SUBROUTINE CALTBXY CALCULATES BOTTOM FRICTION OR DRAG  
!! COEFFICIENTS IN QUADRATIC LAW FORM REFERENCED TO NEAR  
!! BOTTOM OR DEPTH AVERAGED HORIZONTAL VELOCITIES  
!! FOR VEGETATION RESISTANCE IN DEPTH INTEGRATED FLOW  
!! THE COEFFICIENT REPRESENTS BOTTOM AND WATER COLUMN VEGETATION  
!! RESISTANCE 

! @date  2019-06       PAUL M. CRAIG     CHANGED DRAG COEFFICIENT APPROACH FOR CONSISTENT TREATMENT OF LAYERS
!!                                      AND THIN BOTTOM LAYERS RELATIVE TO Z0
!!       2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
!!       2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
!!       2013-08-05    DANG CHUNG        Correction of drag coefficient for wave conditions
!!       2011-04-26    DANG CHUNG        Corrected Wave Formulation & Orbital Velocity
!!       2011-03-02    PAUL M. CRAIG     RESTRUCTURED AND CORRECTED CODE, ADDED OMP
!!       2011-01-15    PAUL M. CRAIG     ADDED THE FRACTIONAL SCALING FOR PARTIAL PENETRATION OF VEGETATION
!!       2010-XX-XX    SCOTT JAMES       ADDED MHK
!!       11/08/2001    john hamrick      REMOVED DRAG COEFFICIENT CONSTRAINT FOR MULIPLE LAYER ROUGHNESS
!!                                         BOUNDARIES WHEN DYNAMIC TIME STEPPING IS ACTIVE
!!       01/28/2002    john hamrick      FIXED POSSIBLE DIVIDE BY ZERO FOR SUB GRID CHANNEL FRICTION IN 
!!                                      ABSENCE OF VEGETATION RESISTANCE

SUBROUTINE CALTBXY

  use GLOBAL  
  use Variables_WQ
  
  use MPI
  use Variables_MPI 
  use Communicate_Ghost_Routines
  
  implicit none

  integer :: L, K, LS, M, LW, LE, LN, LNW, LSE, MW, MS, LF, LL, ND, LP, LG, IOBC
  integer :: NMD, LHOST, LCHNU, LCHNV, MH, MU, MV, NTMP, LDMW, LWAVE, NAL, IERR

  real :: CDTOTUM, CDTOTVM, CDMAXUM, CDMAXVM
  real :: UMAGTMP, VMAGTMP, CDMAXU, CDMAXV, ZBRATUC, ZBRATVC
  real :: HURTMP, HVRTMP, HUDZBR, HVDZBR, VTMPATU, UTMPATV, CPVEGU
  real :: CPVEGV, HVGTC, HVGTW, HVGTS, VISEXP, VISFAC, VISMUDU
  real :: VISMUDV, SEDTMP, CSEDVIS, VISDHU, VISDHV, DZHUDZBR, DZHVDZBR
  real :: FRACLAY, FHLAYC, FHLAYW, FHLAYS, WCHAN, RLCHN, HCHAN, STBXCH
  real :: FXVEGCH, STBYCH, FYVEGCH, TMPVALW, WVFACT, QQWCTMP, TWCTMP
  real :: AEXTMP, TMPVAL, USTARC, CDRGTMP, TAUBTMP, TAUE, RIPAMP
  real :: RIPSTP, RIPFAC, ZBREU
  real :: WVDTMP, RKZTURB, UTMP, VTMP, DWVDZ, DWUDZ, DWVD2Z
  real :: DWUD2Z, HZRVDZ, HZRUDZ, ZDHZRV, ZDHZRU, ZBREV, HZREFV, HZREFU
  real :: QWDQCV, QWDQCU, QCTMPV, QCTMPU, CDTMPVY
  real :: BOTTMP, DWVDHR, DWUDHR, QWCTMPV, QWCTMPU
  real :: CDTMPV, CDTMPU, COSWC, CURANG, CDTMPUX
  real :: WVDELV, WVDELU, TAUTMP, QQWVTMP
  real(rkd) :: ztemp, rr, z0b_gotm, GTAUB
  integer :: itr, itz0b = 10
  
  logical, save, allocatable :: LOCALWAVEMASK(:)
  real, save :: CDLIMIT
  real, save, allocatable :: ZBRATU(:)
  real, save, allocatable :: ZBRATV(:)
  real, save, allocatable :: SGZUU(:,:)
  real, save, allocatable :: SGZVV(:,:)

  DELT = DT2  
  if( ISTL /= 3 )then  
    DELT = DT  
  endif  
  if( IS2TL == 1 )then  
    if( ISDYNSTP == 0 )then  
      DELT = DT  
    else  
      DELT = DTDYN  
    endif  
  endif  
  DELTI = 1./DELT  
  
  ! ***  INITIALIZE IMPLICIT BOTTOM FRICTION AND SET DIAGNOSTIC FILES ON FIRST CALL  
  if( .not. allocated(ZBRATU) )then
    ! *** FIRST CALL
    allocate(LOCALWAVEMASK(LCM))
    do L = 2,LA
      LOCALWAVEMASK(L) = .not. LWVMASK(L)
    enddo
    
    ! *** ZBRATU & ZBRATV - AVERAGE WEIGHTED ROUGHNESS HEIGHTS
    allocate(ZBRATU(LCM))
    allocate(ZBRATV(LCM))
    allocate(SGZUU(LCM,KCM))
    allocate(SGZVV(LCM,KCM))
    
    do L = 2,LA  
      LW = LWC(L)
      LS = LSC(L)
      ZBRATU(L) = 0.5*(DXP(LW)*ZBR(LW) + DXP(L)*ZBR(L))*DXIU(L)
      ZBRATV(L) = 0.5*(DYP(LS)*ZBR(LS) + DYP(L)*ZBR(L))*DYIV(L)  
      do K = 1,KC
        SGZUU(L,K) = 0.5*MAX(DZC(L,K),SGZU(L,K))
        SGZVV(L,K) = 0.5*MAX(DZC(L,K),SGZV(L,K))
      enddo
    enddo

    ! *** Set drag coefficient limits and weighting fractions
    if( ISITB >= 1 )then  
      ! *** ISVEG > 0 case and implicit drag is used
      if( ISITB == 1 )then  
        RITB1 = 0.45      ! *** WEIGHTED EXPLICIT FRACTION
        RITB = 0.55       ! *** WEIGHTED IMPLICIT FRACTION
        CDLIMIT = 10.  
      else  
        RITB1 = 0.0       ! *** ZERO EXPLICIT FRACTION  
        RITB = 1.0        ! *** FULL IMPLICIT FRACTION
        CDLIMIT = 100.  
      endif  
    else  
      RITB1 = 1.0         ! *** FILL EXPLICIT FRACTION  
      RITB = 0.0          ! *** ZERO IMPLICIT FRACTION 
      CDLIMIT = 0.5  
    endif  
    
    do L = 2,LA  
      STBXO(L) = STBX(L)  
      STBYO(L) = STBY(L)  
    enddo  
    do L = 1,LC  
      STBX(L) = 0.        ! *** STBX - DRAG COEFFICIENT IN THE U DIRECTION
      STBY(L) = 0.        ! *** STBY - DRAG COEFFICIENT IN THE V DIRECTION
      ZBRE(L) = KSW/30.   ! *** ZBRE - BED ROUGHNESS USED FOR WAVE CALCULATIONS
    enddo
    
    if( ISVEG > 0 )then
      do K = 1,KC  
        do L = 1,LC  
          FXVEG(L,K) = 0.  
          FYVEG(L,K) = 0.  
        enddo  
      enddo  

      ! *** GET THE VEGETATION CELL FLAG
      LVEG = .FALSE.
      do L = 2,LA
        M = MVEGL(L)  
        if( M /= MVEGOW .and. M /= 0 )then
          LVEG(L) = .TRUE.                                ! *** Not open water or uninitialized devices
        endif
      enddo
      
      ! *** GET THE VEGETATION CELL FLAG FOR MACROPHYTES
      if( ISTRAN(8) > 0 )then
        do L = 2,LA
          do NAL = 1,NALGAE
            if( .not. ALGAES(NAL).ISMOBILE .and. ALGAES(NAL).ISDRAG > 0 )then
              if( SUM(WQV(L,:,19 + NAL)) >= 0.001 )then
                LVEG(L) = .TRUE.                          ! *** This is a macrophyte cell that includes drag
                MAC_CELL(L) = .TRUE.
              endif
            endif
          enddo
        enddo
        if( UVEGSCL <= 0. ) UVEGSCL = 1.E-12
      endif
      
      ! *** GET THE VEGETATION CELL COUNT AND LIST BY LAYER
      allocate(LLVEG(KCM,NDM))
      allocate(LKVEG(LCM,KCM,NDM))
      LLVEG = 0
      LKVEG = 0
      do ND = 1,NDM  
        do K = 1,KC  
          LN = 0
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            if( LVEG(L) )then
              LN = LN + 1
              LKVEG(LN,K,ND) = L
            endif
          enddo
          LLVEG(K,ND) = LN    ! *** NUMBER OF WET VEG CELLS FOR THE CURRENT LAYER
        enddo
      enddo

      ! *** DEACTIVATE VEGETATION FOR OPEN BC CELLS
      if( NBCSOP > 0 )then
        do IOBC = 1,NBCSOP  
          L = LOBCS(IOBC)  
          MVEGL(L) = 0.
        enddo  
        do IOBC = 1,NBCSOP2 
          L = LOBCS2(IOBC)  
          MVEGL(L) = 0.
        enddo
      endif
    endif  ! *** END OF ISVEG > 0 
    N = -2  

    ! *** CHECK FOR SMOOTH DRAG FORMULATION
    do L = 2,LA
      if( ZBR(L) <= 1.E-6 )then
        ICALTB = 1   ! *** SET BOTTOM DRAG APPROACH TO EFDCPLUS APPROACH
        EXIT
      endif
    enddo
    
    ! ****************************************************************************
    ! *** MPI communication
    call MPI_barrier(MPI_Comm_World, ierr)
    call communicate_ghost_cells(ZBRATU, 'ZBRATU')
    call communicate_ghost_cells(ZBRATV, 'ZBRATV')
    call communicate_ghost_cells(SGZUU)
    call communicate_ghost_cells(SGZVV)
    ! ****************************************************************************

  endif   ! *** End of first call initializations

  if( ISWAVE > 0 )then
    ! *** COMPUTE RAMPUP FACTOR
    NTMP = MAX(N,1)  
    if( NTMP < NTSWV )then  
      TMPVALW = FLOAT(NTMP)/FLOAT(NTSWV)  
      WVFACT = 0.5-0.5*COS(PI*TMPVALW)  
    else  
      WVFACT = 1.0  
    endif  
    RKZTURB = 0.4/CTURB3
    LDMW = INT(FLOAT(NWVCELLS)/FLOAT(NDM)) + 1
  endif
  
  ! ***  INITIALIZED DIAGNOSTICS FOR STANDARD AND VEGE RESISTANCE CALCULATION  
  CDTOTUM = 0.  
  CDTOTVM = 0.  
  CDMAXUM = 0.  
  CDMAXVM = 0.  
  if( ISVEG == 0 ) UVEGSCL = 1.E-12

  ! *********************************************************************************
  ! *** NORMAL ENTRY INTO STANDARD AND VEGE RESISTANCE CALCULATION 
  ! ***  FOR SINGLE AND MULTIPLE LAYERS
  VISEXP = 2./7.
  VISFAC = 0.0258*(COEFTSBL**VISEXP)
  
  ! *** ZERO NEWLY DRY CELLS
  if( LADRY > 0 )then
    do LP = 1,LADRY
      L = LDRY(LP)  
      VEGK(L) = 0.
      QQ(L,0) = 0.
      QQSQR(L,0) = 0.
      TBX(L) = 0.
      TBY(L) = 0.
      TBX(LEC(L)) = 0.
      TBY(LNC(L)) = 0.
      TSX(L) = 0.
      TSY(L) = 0.
      TSX(LEC(L)) = 0.
      TSY(LNC(L)) = 0.
    enddo
 
    if( ISVEG > 0 )then
      do K = 1,KC
        do LP = 1,LADRY
          L = LDRY(LP)  
          FXVEG(L,K) = 0.
          FYVEG(L,K) = 0.
        enddo
      enddo
    endif
  endif
  
  !$OMP PARALLEL DEFAULT(SHARED) 

  ! *********************************************************************************
  ! *** WAVE-CURRENT BOUNDARY LAYER  
  if( (ISWAVE == 2 .or. ISWAVE == 4) .and. .not. LSEDZLJ )then
    !$OMP DO PRIVATE(ND,L,LF,LL,LWAVE)  &
    !$OMP    PRIVATE(QQWCTMP,TWCTMP,AEXTMP,TAUTMP,TMPVAL,USTARC,CDRGTMP) &
    !$OMP    PRIVATE(TAUBTMP,TAUE,RIPAMP,RIPSTP,RIPFAC,QQWVTMP)
    do ND = 1,NDM  
      LF = 1 + (ND-1)*LDMW  
      LL = MIN(LF + LDMW-1,NWVCELLS)

      ! *** Update water column turbulent intensity by wave action
      do LWAVE = LF,LL  
        L = LWVCELL(LWAVE)

        ! *** SET ZBRE AS GRAIN/SKIN ROUGHNESS (M)  (NIKURADSE ROUGHNESS)
        ZBRE(L) = KSW/30.

        if( UWVSQ(L) > 1.E-6 .and. LMASKDRY(L) .and. HP(L) > HDRYWAV )then  
          
          ! *** QQ(L,0) - (m2/s2) N-1 Total Bed Shear Turbulent Intensity
          ! *** QQWV2   - (m2/s2) N-1 Water Column Turbulent Intensity due to waves only
          QQWCTMP = SQRT( QQWV2(L)*QQWV2(L) + QQ(L,0)*QQ(L,0) )          

          TWCTMP = QQWCTMP/CTURB2
          AEXTMP = 0.5*WV(L).HEIGHT/SINH(WV(L).KHP)                 
          if( ISTRAN(7) > 0 )then  
            TAUTMP = TWCTMP/TAUCMIN
            TMPVAL = 1. + 1.2*TAUTMP/(1. + 0.2*TAUTMP)
            ZBRE(L) = ZBRE(L)*TMPVAL
          elseif( QQ(L,0) > 0. )then
            USTARC = SQRT(QQ(L,0)/CTURB2)  
            TMPVAL = TWCTMP/USTARC  
            ZBRE(L) = ZBRE(L)*(1. + 0.19*TMPVAL)  
          endif
          CDRGTMP = (30.*ZBRE(L)/AEXTMP)**0.2  
          CDRGTMP = 5.57*CDRGTMP-6.13  
          CDRGTMP = EXP(CDRGTMP)  
          CDRGTMP = MIN(CDRGTMP,0.22)  
          TAUTMP = 0.5*CDRGTMP*UWVSQ(L)   

          ! *** Recalculate Turbulent Intensity due to waves
          
          QQWVTMP = CTURB2*TAUTMP*WVFACT
          
          ! *** COMPUTE MOVING BED EFFECTS
          if( ISTRAN(7) > 0 .and. ISWCBL == 2 )then
            QQWC(L) = SQRT( QQWVTMP*QQWVTMP + QQ(L,0)*QQ(L,0) )   ! *** Combined wave/current for stationary bed (m2/s2)
            TWCTMP = QQWC(L)/CTURB2
            TAUBTMP = QQWVTMP/CTURB2
            TAUE = TWCTMP/TAUN(NSED + 1)
            RIPAMP = 0.
            RIPSTP = 0.
            if( TAUBTMP>TAUN(NSED + 1) .and. TAUBTMP <= TAUD(NSED + 1) )then
             RIPAMP = 0.22/(TAUE**0.16)
             RIPSTP = 0.16/(TAUE**0.04)
            endif
            if( TAUBTMP>TAUD(NSED + 1) )then
             RIPAMP = 0.78/(TAUE**1.5)
             RIPSTP = 0.41/TAUE
            endif
            RIPAMP = RIPAMP*WV(L).HEIGHT/SINH(WV(L).KHP)
            TMPVAL = 0.
            if( RIPAMP>0.) TMPVAL = LOG(RIPAMP/ZBRE(L))-1.
            TMPVAL = MAX(TMPVAL,0.)
            RIPFAC = 1. + 3.125*TMPVAL*TMPVAL*RIPSTP
            QQWV3(L) = RIPFAC*QQWVTMP
          else
            QQWV3(L) = QQWVTMP  
          endif

          ! *** Maximum magnitude without accounting for wave/current interactions (m2/s2)
          QQWCR(L) = SQRT( QQWV3(L)*QQWV3(L) + QQ(L,0)*QQ(L,0) )

        else  
          QQWV3(L) = 1.E-12
          QQWCR(L) = 1.E-12
        endif  
      enddo   ! *** END OF LWAVE LOOP
    enddo     ! *** END OF DOMAIN
    !$OMP END DO
  
    ! *** Bed shear stress by current & wave
    !$OMP DO PRIVATE(ND,L,LF,LL,LE,LS,LN,LW,LWAVE)  &
    !$OMP    PRIVATE(UTMP,VTMP,CURANG,COSWC,UMAGTMP,VMAGTMP,CDMAXU,CDMAXV,CDTMPU,CDTMPV,WVDTMP,WVDELU,WVDELV)  &
    !$OMP    PRIVATE(QWCTMPU,QWCTMPV,QCTMPU,QCTMPV,QWDQCU,QWDQCV,HZREFU,HZREFV,ZBREU,ZBREV,BOTTMP )  &
    !$OMP    PRIVATE(ZDHZRU,ZDHZRV,HZRUDZ,HZRVDZ,DWUD2Z,DWVD2Z,DWUDZ,DWVDZ,DWUDHR,DWVDHR,CDTMPUX,CDTMPVY)
    do ND = 1,NDM  
      LF = 1 + (ND-1)*LDMW  
      LL = MIN(LF + LDMW-1,NWVCELLS)

      ! *** Update bed drag coefficients STBX and STBY based on wave and current angles
      do LWAVE = LF,LL  
        L = LWVCELL(LWAVE)

        if( UWVSQ(L) > 1.E-6 .and. LMASKDRY(L) )then
          LE = LEC(L)
          LS = LSC(L)  
          LN = LNC(L)
          LW = LWC(L)
              
          ! *** Avg Velocities at U & V Faces
          UTMP = 0.5*STCUV(L)*(U(LE,KSZ(LE)) + U(L,KSZ(L))) + 1.E-12  
          VTMP = 0.5*STCUV(L)*(V(LN,KSZ(LN)) + V(L,KSZ(L)))  
          CURANG = ATAN2(VTMP,UTMP)
            
          ! *** Cosine of the Current - Wave  
          COSWC = COS(CURANG-WV(L).DIR)
            
          ! *** Velocites at Corners  
          UMAGTMP = SQRT( U1(L,KSZ(L))*U1(L,KSZ(L)) + V1U(L)      *V1U(L)      +1.E-12 )  
          VMAGTMP = SQRT( U1V(L)      *U1V(L)       + V1(L,KSZ(L))*V1(L,KSZ(L)) + 1.E-12 )
            
          ! *** Set Initial Drag Coefficients  
          CDMAXU = STBXO(L)*H1U(L)/( 4.*DELT*UMAGTMP )  
          CDMAXV = STBYO(L)*H1V(L)/( 4.*DELT*VMAGTMP )  
          CDTMPU = -1.  
          CDTMPV = -1.
            
          ! *** Avg Wave Turbulence at U & V Faces  
          QWCTMPU = 0.5*( QQWV3(L) + QQWV3(LE) )  
          QWCTMPV = 0.5*( QQWV3(L) + QQWV3(LS) )  
            
          if( ISWCBL == 2 )then  
            QWCTMPU = 0.5*( QQWC(L) + QQWC(LE) )  
            QWCTMPV = 0.5*( QQWC(L) + QQWC(LS) )  
          endif
             
          if( WV(L).FREQ > 1E-6 )then
            WVDTMP = 0.4/(WV(L).FREQ*CTURB3)
          else
            WVDTMP = 0.
          endif  
          WVDELU = WVDTMP*SQRT(QWCTMPU)  
          WVDELV = WVDTMP*SQRT(QWCTMPV)
            
          ! *** Avg Wave & Current Stress  
          QWCTMPU = 0.5*( QQWCR(L) + QQWCR(LE) )  
          QWCTMPV = 0.5*( QQWCR(L) + QQWCR(LS) )
            
          ! *** Avg Shear Velocities  
          QWCTMPU = SQRT(QWCTMPU)  
          QWCTMPV = SQRT(QWCTMPV)
            
          ! *** Avg Total Bottom Shear  
          QCTMPU = 0.5*( QQ(L,0) + QQ(LE,0) )  
          QCTMPV = 0.5*( QQ(L,0) + QQ(LS,0) ) 
            
          if( QCTMPU /= 0 )then
            QWDQCU = QWCTMPU/SQRT(QCTMPU)  
          else
            QWDQCU = 0
          endif
          if( QCTMPV /= 0 )then
            QWDQCV = QWCTMPV/SQRT(QCTMPV)
          else
            QWDQCV = 0
          endif
            
          ! *** Thickness of Bottom Layer  
          HZREFU = DZC(L,KSZ(L))*H1U(L)  
          HZREFV = DZC(L,KSZ(L))*H1V(L)
            
          ! *** Avg Bottom roughness  
          ZBREU = 0.5*(ZBRE(L) + ZBRE(LE))  
          ZBREV = 0.5*(ZBRE(L) + ZBRE(LS))
            
          ! *** Ratio of Roughness to layer thickness  
          ZDHZRU = ZBREU/HZREFU  
          ZDHZRV = ZBREV/HZREFV  
          HZRUDZ = 1./ZDHZRU  
          HZRVDZ = 1./ZDHZRV  
          DWUD2Z = 0.5*WVDELU/ZBREU  
          DWVD2Z = 0.5*WVDELV/ZBREV  
          DWUDZ = 2.*DWUD2Z  
          DWVDZ = 2.*DWVD2Z  
          DWUDHR = WVDELU/HZREFU  
          DWVDHR = WVDELV/HZREFV  
          CDTMPUX = RKZTURB*QWCTMPU  
          CDTMPVY = RKZTURB*QWCTMPV  
          if( HZRUDZ <= DWUD2Z )then  
            CDTMPU = CDTMPUX/( (1. + ZDHZRU)*LOG(1. + HZRUDZ)-1. )  
          endif  
          if( HZRVDZ <= DWVD2Z )then  
            CDTMPV = CDTMPVY/( (1. + ZDHZRV)*LOG(1. + HZRVDZ)-1. )  
          endif  
          if( HZRUDZ > DWUD2Z .and. HZRUDZ <= DWUDZ )then  
            BOTTMP = (1. + ZDHZRU)*LOG(1. + DWUD2Z)-0.5*DWUDHR + 0.5*HZRUDZ*(1.-0.5*DWUDHR)*(1.-0.5*DWUDHR)/(1. + DWUD2Z)  
            CDTMPU = CDTMPUX/BOTTMP  
          endif  
          if( HZRVDZ > DWVD2Z .and. HZRVDZ <= DWVDZ )then  
            BOTTMP = (1. + ZDHZRV)*LOG(1. + DWVD2Z)-0.5*DWVDHR + 0.5*HZRVDZ*(1.-0.5*DWVDHR)*(1.-0.5*DWVDHR)/(1. + DWVD2Z)  
            CDTMPV = CDTMPVY/BOTTMP  
          endif  
          if( HZRUDZ > DWUDZ )then  
            BOTTMP = QWDQCU*( (1. + ZDHZRU)*(LOG(1. + HZRUDZ)-LOG(1. + DWUDZ)) + DWUDHR-1. )  
            BOTTMP = BOTTMP + (1. + ZDHZRU)*LOG(1. + DWUD2Z) + DWUD2Z*(1.-1.25*DWUDHR-ZDHZRU)/(1. + DWUD2Z)
            CDTMPU = CDTMPUX/BOTTMP  
          endif  
          if( HZRVDZ > DWVDZ )then  
            BOTTMP = QWDQCV*( (1. + ZDHZRV)*(LOG(1. + HZRVDZ)-LOG(1. + DWVDZ)) + DWVDHR-1. )  
            BOTTMP = BOTTMP + (1. + ZDHZRV)*LOG(1. + DWVD2Z) + DWVD2Z*(1.-1.25*DWVDHR-ZDHZRV)/(1. + DWVD2Z)  
            CDTMPV = CDTMPVY/BOTTMP  
          endif  
          CDTMPU = CDTMPU/UMAGTMP  
          CDTMPV = CDTMPV/VMAGTMP
            
          if( CDTMPU <= 0.) CDTMPU = CDMAXU  
          if( CDTMPV <= 0.) CDTMPV = CDMAXV  
            
          STBX(L) = STBXO(L)*CDTMPU
          STBY(L) = STBYO(L)*CDTMPV  
          STBX(L) = MIN(CDMAXU,STBX(L))
          STBY(L) = MIN(CDMAXV,STBY(L))
 
        endif
      enddo  ! *** END OF LWAVE LOOP
    enddo    ! *** END OF DOMAIN
    !$OMP END DO
  
  elseif( (ISWAVE == 1 .or. ISWAVE == 3) .and. .not. LSEDZLJ )then
    !$OMP SINGLE
    QQWV3 = QQWV1
    !$OMP END SINGLE
  endif      ! *** END OF BLOCK FOR WAVE-BOUNDARY LAYER APPROACH

  ! ******************************************************************************
  ! *** BED STRESS FOR NON-WAVE CELLS
  
  !$OMP DO PRIVATE(ND,LF,LL,LP,L)  &
  !$OMP    PRIVATE(UMAGTMP,VMAGTMP,CDMAXU,CDMAXV,VISMUDU,VISMUDV,SEDTMP,VISDHU,VISDHV)  &
  !$OMP    PRIVATE(HURTMP,HVRTMP,HUDZBR,HVDZBR,DZHUDZBR,DZHVDZBR)                       &
  !$OMP    PRIVATE(itr,rr,ztemp,z0b_gotm,GTAUB)
  do ND = 1,NDM  
    LF = (ND-1)*LDMWET + 1  
    LL = MIN(LF + LDMWET-1,LAWET)
    
    if( ISGOTM > 0 .and. IFRICTION == 1 )then  
      do LP = LF,LL
        L = LWET(LP)
        ztemp = HPK(L,KSZ(L))
        GTAUB = 0.
        do itr = 1, itz0b    
          z0b_gotm = 0.1*AVO/max(AVO,GTAUB) + 0.03*ZBR(L)    
          rr = VKC/(log((z0b_gotm + ztemp)/z0b_gotm))          
          GTAUB = rr*sqrt(UCTR(L,KSZ(L))*UCTR(L,KSZ(L)) + VCTR(L,KSZ(L))*VCTR(L,KSZ(L)))
        enddo
         STBX(L) = rr
         STBY(L) = rr
      enddo
    elseif( ICALTB == 0 )then
      ! *** DSI DEFAULT.  CONSISTENT WITH THEORY AND HANDLES THIN BOTTOM LAYERS
      do LP = LF,LL
        L = LWET(LP)  
        if( LOCALWAVEMASK(L) .or. QQWV2(L) <= 1.E-12 )  THEN 
          ! *** SET DRAG COEFFICIENT LIMITS
          UMAGTMP = SQRT( VU(L)*VU(L) + U(L,KSZU(L))*U(L,KSZU(L)) + 1.E-12 )  
          VMAGTMP = SQRT( UV(L)*UV(L) + V(L,KSZV(L))*V(L,KSZV(L)) + 1.E-12 )  
          CDMAXU = CDLIMIT*HU(L)/( DELT*UMAGTMP )  
          CDMAXV = CDLIMIT*HV(L)/( DELT*VMAGTMP )  
      
          ! *** HU & HV - FLOW DEPTHS AT CELL INTERFACE
          HUDZBR = HU(L)/ZBRATU(L)
          HVDZBR = HV(L)/ZBRATV(L)  
          HUDZBR = MAX(HUDZBR,7.5)
          HVDZBR = MAX(HVDZBR,7.5)
            
          ! *** NEZU & NAKAGAWA (1993) WHERE THE WAKE parameter IS 0.2  FOR TURBULENT CONDITIONS.  -0.8 = WAKE - 1.
          STBX(L) = (VKC/(LOG( HUDZBR ) - 0.8))**2 
          STBY(L) = (VKC/(LOG( HVDZBR ) - 0.8))**2  
              
          STBX(L) = MIN(CDMAXU, STBX(L) )
          STBY(L) = MIN(CDMAXV, STBY(L) )
        endif
      enddo
      
    else
      ! *** LEGACY APPROACHES
      do LP = LF,LL
        L = LWET(LP)  
        if( LOCALWAVEMASK(L) .or. QQWV2(L) <= 1.E-12 )  THEN 
          if( ZBR(L) <= 1.E-6 )then  
            ! ***  BEGIN SMOOTH DRAG FORMULATION  
            UMAGTMP = SQRT( U(L,KSZ(L))*U(L,KSZ(L)) + VU(L)      *VU(L)       + 1.E-12 )  
            VMAGTMP = SQRT( UV(L)      *UV(L)       + V(L,KSZ(L))*V(L,KSZ(L)) + 1.E-12 )  
            CDMAXU = CDLIMIT*STBXO(L)*HU(L)/( DELT*UMAGTMP )  
            CDMAXV = CDLIMIT*STBYO(L)*HV(L)/( DELT*VMAGTMP )
            
            VISMUDU = VISMUD
            VISMUDV = VISMUD
            if( ISMUD >= 1 )then  
              SEDTMP = 0.5*(SED(L,KSZ(L),1) + SED(LWC(L),KSZ(LWC(L)),1))  
              VISMUDU = CSEDVIS(SEDTMP)  
              SEDTMP = 0.5*(SED(L,KSZ(L),1) + SED(LSC(L),KSZ(LSC(L)),1))  
              VISMUDV = CSEDVIS(SEDTMP)  
            endif  
            VISDHU = 0.0
            VISDHV = 0.0
            if( UMAGTMP > 0.0) VISDHU = (VISMUDU*HUI(L)/UMAGTMP)*VISEXP
            if( VMAGTMP > 0.0) VISDHV = (VISMUDV*HVI(L)/VMAGTMP)*VISEXP
            STBX(L) = VISFAC*STBXO(L)*VISDHU
            STBY(L) = VISFAC*STBYO(L)*VISDHV
            STBX(L) = MIN(CDMAXU,STBX(L))  
            STBY(L) = MIN(CDMAXV,STBY(L))  
        
          elseif( ISAVCOMP > 0 )then  
            ! ***  BEGIN ROUGH DRAG FORMULATION  
            if( ICALTB == 1 )then
              ! *** PRE DSI VERSION 10.0 APPROACH
              UMAGTMP = SQRT( VU(L)*VU(L) + U(L,KSZU(L))*U(L,KSZU(L)) + 1.E-12 )  
              VMAGTMP = SQRT( UV(L)*UV(L) + V(L,KSZV(L))*V(L,KSZV(L)) + 1.E-12 )  
              CDMAXU = CDLIMIT*STBXO(L)*HU(L)/( DELT*UMAGTMP )  
              CDMAXV = CDLIMIT*STBYO(L)*HV(L)/( DELT*VMAGTMP )  
    
              ! *** HURTMP & HVRTMP - FLOW DEPTHS AT FACE
              HURTMP = MAX(ZBRATU(L), H1U(L))
              HVRTMP = MAX(ZBRATV(L), H1V(L))
    
              ! *** HANDLE LAYER ISSUES              
              if( KSZ(L) == KC )then    ! *** Alberta
                HUDZBR = 0.5*HURTMP/ZBRATU(L)  
                HUDZBR = MAX(HUDZBR,7.5)
                HVDZBR = 0.5*HVRTMP/ZBRATV(L)  
                HVDZBR = MAX(HVDZBR,7.5)
                STBX(L) = STBXO(L)*0.16/( (LOG( HUDZBR ) -1.)**2)  
                STBY(L) = STBYO(L)*0.16/( (LOG( HVDZBR ) -1.)**2)  
              else                
                DZHUDZBR = 1. + SGZUU(L,KSZU(L))*HURTMP/ZBRATU(L)
                DZHVDZBR = 1. + SGZVV(L,KSZV(L))*HVRTMP/ZBRATV(L)
                DZHUDZBR = MAX(DZHUDZBR,7.5)  ! *** APPLY THE SAME LIMIT AS KC = 1
                DZHVDZBR = MAX(DZHVDZBR,7.5)  ! *** APPLY THE SAME LIMIT AS KC = 1
                STBX(L) = STBXO(L)*0.16/((LOG(DZHUDZBR))**2)  
                STBY(L) = STBYO(L)*0.16/((LOG(DZHVDZBR))**2) 
              endif
              STBX(L) = MIN(CDMAXU, STBX(L) )
              STBY(L) = MIN(CDMAXV, STBY(L) )
    
            elseif( ICALTB == 2 )then
              ! *** LEGACY VERSION
              UMAGTMP = SQRT( V1U(L)*V1U(L) + U1(L,KSZU(L))*U1(L,KSZU(L)) + 1.E-12 )  
              VMAGTMP = SQRT( U1V(L)*U1V(L) + V1(L,KSZV(L))*V1(L,KSZV(L)) + 1.E-12 )  
              CDMAXU = CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )  
              CDMAXV = CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP )  
    
              ! *** HURTMP & HVRTMP - FLOW DEPTHS AT FACE
              HURTMP = MAX(ZBRATU(L), H1U(L))
              HVRTMP = MAX(ZBRATV(L), H1V(L))
    
              DZHUDZBR = 1. + SGZUU(L,KSZU(L))*HURTMP/ZBRATU(L)
              DZHVDZBR = 1. + SGZVV(L,KSZV(L))*HVRTMP/ZBRATV(L)
            
              STBX(L) = STBXO(L)*.16/((LOG(DZHUDZBR))**2)  
              STBY(L) = STBYO(L)*.16/((LOG(DZHVDZBR))**2) 
              STBX(L) = MIN(CDMAXU, STBX(L) )
              STBY(L) = MIN(CDMAXV, STBY(L) )
            endif
          endif  
        endif  
      enddo  
    endif
  enddo   ! *** END OF DOMAIN
  !$OMP END DO

  ! *********************************************************************************
  ! *** BEGIN INTERNAL MODE VEGETATION DRAG
  if( ISVEG >= 1 )then  
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,LS,LW,LE,LN,LNW,LSE,MW,M,MS,NAL)  &
    !$OMP    PRIVATE(VTMPATU,UTMPATV,UMAGTMP,VMAGTMP,CPVEGU,CPVEGV )  &
    !$OMP    PRIVATE(HVGTC,HVGTW,HVGTS,FRACLAY,FHLAYC,FHLAYW,FHLAYS,CDMAXU,CDMAXV)  
    do ND = 1,NDM  

      do LP = 1,LLVEG(KC,ND)
        ! *** ZERO ACTIVE VEGETATION LAYERS FOR WET CELLS
        L = LKVEG(LP,KC,ND) 
        VEGK(L) = 0.
      enddo

      do K = 1,KC  
        do LP = 1,LLVEG(K,ND)
          L = LKVEG(LP,K,ND)
          
          FXVEG(L,K) = 0.  
          FYVEG(L,K) = 0.  

          LW = LWC(L)
          LE = LEC(L)
          LS = LSC(L)  
          LN = LNC(L)  
          LNW = LNWC(L)  
          LSE = LSEC(L)  
          
          VTMPATU = 0.25*(V(L,K) + V(LW,K) + V(LN,K) + V(LNW,K))         ! *** Average U velocity at the V face
          UTMPATV = 0.25*(U(L,K) + U(LE,K) + U(LS,K) + U(LSE,K))         ! *** Average V velocity at the U face  
          UMAGTMP = SQRT( U(L,K)*U(L,K)   + VTMPATU*VTMPATU + 1.E-12 )   ! *** Magnitude of velocity at U face
          VMAGTMP = SQRT( UTMPATV*UTMPATV + V(L,K)*V(L,K)   + 1.E-12 )   ! *** Magnitude of velocity at V face
                  
          UMAGTMP = MAX(UMAGTMP,UVEGSCL)  
          VMAGTMP = MAX(VMAGTMP,UVEGSCL)  
          CDMAXU  = CDLIMIT*STBXO(L)*H1U(L)/( DELT*UMAGTMP )
          CDMAXV  = CDLIMIT*STBYO(L)*H1V(L)/( DELT*VMAGTMP ) 
          
          ! *** Handle macrophyte feedback in WQ model
          if( MAC_CELL(L) )then
            do NAL = 1, NALGAE
              if( .not. ALGAES(NAL).ISMOBILE .and. ALGAES(NAL).ISDRAG > 0 )then
                ! *** DRAG COEFFICIENT: U COMPONENT 
                CPVEGU = 1.0
                if( ALGAES(NAL).ISMACL == 1 ) CPVEGU = CPVEGU + 10.E-6/((DIAMETER_MAC(LW,NAL) + DIAMETER_MAC(L,NAL))*UMAGTMP )  
                if( CPVEGU > 1.0 )then  
                  ! *** CALCULATE R FOR LAMINAR FLOW  
                  CPVEGU = CPVEGU-0.5  
                endif
                CPVEGU = ALGAES(NAL).DRAGCOEFF*CPVEGU  
                
                ! *** DRAG COEFFICIENT: V COMPONENT 
                CPVEGV = 1.0
                if( ALGAES(NAL).ISMACL == 1 ) CPVEGV = CPVEGV + 10.E-6/((DIAMETER_MAC(LS,NAL) + DIAMETER_MAC(L,NAL))*VMAGTMP )  
                if( CPVEGV > 1.0 )then  
                  ! *** CALCULATE R FOR LAMINAR FLOW  
                  CPVEGV = CPVEGV-0.5  
                endif  
                CPVEGV = ALGAES(NAL).DRAGCOEFF*CPVEGV
                            
                ! *** HANDLE LAYER ISSUES
                if( KSZ(L) == KC )then
                  HVGTC = MIN(HEIGHT_MAC(L ,NAL),HP(L ))  
                  HVGTW = MIN(HEIGHT_MAC(LW,NAL),HP(LW))  
                  HVGTS = MIN(HEIGHT_MAC(LS,NAL),HP(LS))  
                else   !IF( MVEGL(L)<91 )then
                  FRACLAY = Z(L,K)
                  FHLAYC = FRACLAY*HP(L)  
                  FHLAYW = FRACLAY*HP(LW)  
                  FHLAYS = FRACLAY*HP(LS)  
                  HVGTC = HPK(L,K)  
                  HVGTW = HPK(LW,K)
                  HVGTS = HPK(LS,K)  
                  if( HEIGHT_MAC(L,NAL) < FHLAYC )then         
                    ! *** GRADUALLY DECREASE MACROPHYTE EFFECTS  (PMC)
                    HVGTC = HEIGHT_MAC(L,1)-HP(L)*Z(L,K-1)
                    HVGTC = MAX(HVGTC,0.)
                  endif
                  if( HEIGHT_MAC(LW,NAL) < FHLAYW )then 
                    ! *** GRADUALLY DECREASE MACROPHYTE EFFECTS  (PMC)
                    HVGTW = HEIGHT_MAC(LW,NAL)-HP(LW)*Z(L,K-1)
                    HVGTW = MAX(HVGTW,0.)
                  endif
                  if( HEIGHT_MAC(LS,NAL) < FHLAYS )then
                    ! *** GRADUALLY DECREASE MACROPHYTE EFFECTS  (PMC)
                    HVGTS = HEIGHT_MAC(LS,NAL)-HP(LS)*Z(L,K-1)
                    HVGTS = MAX(HVGTS,0.)
                  endif
                endif
                
                ! *** 0.25 factor is 0.5(drag force equation)*0.5(averaging face quantity)
                FXVEG(L,K) = FXVEG(L,K) + 0.25*CPVEGU*(DXP(L) *(MACAD(L,NAL) *HVGTC/MVEGZ(L,NAL)) +   &          ! *** dimensionless
                                          DXP(LW)*(MACAD(LW,NAL)*HVGTW/MVEGZ(LW,NAL)))*DXIU(L)  
                FYVEG(L,K) = FYVEG(L,K) + 0.25*CPVEGV*(DYP(L) *(MACAD(L,NAL) *HVGTC/MVEGZ(L,NAL)) +   &          ! *** dimensionless
                                        DYP(LS)*(MACAD(LS,NAL)*HVGTS/MVEGZ(LS,NAL)))*DYIV(L)
              endif
            enddo
            
            FXVEG(L,K) = MIN(FXVEG(L,K),CDMAXU)  
            FYVEG(L,K) = MIN(FYVEG(L,K),CDMAXV)
            
            VEGK(L) = VEGK(L) + HVGTC*HPKI(L,K)
                    
          else
            ! *** Standard rigid 
            M = MVEGL(L)  
            if( M > 90 ) M = M-90  !SCJ vegetation below MHK device, which is OK
          
            MW = MVEGL(LW)  
            if( MW > 90 ) MW = MW-90  !SCJ, vegetation below MHK device, which is OK
            MS = MVEGL(LS)  
            if( MS > 90 ) MS = MS-90  !SCJ, vegetation below MHK device, which is OK
              
            ! *** DRAG COEFFICIENT: U COMPONENT 
            CPVEGU = 1.0
            if( ISVEGL == 1 ) CPVEGU = CPVEGU + 10.E-6/((BPVEG(MW) + BPVEG(M))*UMAGTMP )  
            if( CPVEGU > 1.0 )then  
              ! *** CALCULATE R FOR LAMINAR FLOW  
              CPVEGU = CPVEGU-0.5  
            endif
            CPVEGU = SCVEG(M)*CPVEGU  

            ! *** DRAG COEFFICIENT: V COMPONENT 
            CPVEGV = 1.0
            if( ISVEGL == 1 ) CPVEGV = CPVEGV + 10.E-6/((BPVEG(MS) + BPVEG(M))*VMAGTMP )  
            if( CPVEGV > 1.0 )then  
              ! *** CALCULATE R FOR LAMINAR FLOW  
              CPVEGV = CPVEGV-0.5  
            endif  
            CPVEGV = SCVEG(M)*CPVEGV
              
            ! *** HANDLE LAYER ISSUES
            if( KSZ(L) == KC )then    ! *** Alberta
              HVGTC = MIN(HPVEG(M),HP(L))  
              HVGTW = MIN(HPVEG(MW),HP(LW))  
              HVGTS = MIN(HPVEG(MS),HP(LS))  
            else   !IF( MVEGL(L)<91 )then
              FRACLAY = Z(L,K) 
              FHLAYC = FRACLAY*HP(L)  
              FHLAYW = FRACLAY*HP(LW)  
              FHLAYS = FRACLAY*HP(LS)  
              HVGTC = HPK(L,K)  
              HVGTW = HPK(LW,K)
              HVGTS = HPK(LS,K)  
              if( HPVEG(M) < FHLAYC )then         
                ! *** GRADUALLY DECREASE VEGETATION EFFECTS  (PMC)
                HVGTC = HPVEG(M)-HP(L)*Z(L,K-1)
                HVGTC = MAX(HVGTC,0.)
              endif
              if( HPVEG(MW) < FHLAYW )then 
                ! *** GRADUALLY DECREASE VEGETATION EFFECTS  (PMC)
                HVGTW = HPVEG(MW)-HP(LW)*Z(L,K-1)
                HVGTW = MAX(HVGTW,0.)
              endif
              if( HPVEG(MS) < FHLAYS )then
                ! *** GRADUALLY DECREASE VEGETATION EFFECTS  (PMC)
                HVGTS = HPVEG(MS)-HP(LS)*Z(L,K-1)
                HVGTS = MAX(HVGTS,0.)
              endif
            endif
            
            ! *** 0.25 factor is 0.5(drag force equation)*0.5(averaging face quantity)
            FXVEG(L,K) = 0.25*CPVEGU*(DXP(L) *(BDLPSQ(M) *HVGTC/PVEGZ(M)) +   &          ! *** dimensionless
                                      DXP(LW)*(BDLPSQ(MW)*HVGTW/PVEGZ(MW)))*DXIU(L)  
            FYVEG(L,K) = 0.25*CPVEGV*(DYP(L) *(BDLPSQ(M) *HVGTC/PVEGZ(M)) +   &          ! *** dimensionless
                                      DYP(LS)*(BDLPSQ(MS)*HVGTS/PVEGZ(MS)))*DYIV(L)  
            FXVEG(L,K) = MIN(FXVEG(L,K),CDMAXU)  
            FYVEG(L,K) = MIN(FYVEG(L,K),CDMAXV)

            ! *** ACCUMULATE ACTIVE VEGETATION LAYERS
            VEGK(L) = VEGK(L) + HVGTC*HPKI(L,K)
          endif
        enddo   ! *** END OF LWET LOOP
      enddo     ! *** END OF KC LOOP  
    enddo       ! *** END OF DOMAIN
    !$OMP END DO
    
  endif    ! *** END OF ISVEG>0
  ! *** END OF VEGETATION DRAG
  !$OMP END PARALLEL

  ! *********************************************************************************
  ! *** SUBGRID SCALE CHANNEL FRICTION  
  if( MDCHH >= 1 )then  
    do NMD = 1,MDCHH  
      LHOST = LMDCHH(NMD)  
      LCHNU = LMDCHU(NMD)  
      LCHNV = LMDCHV(NMD)  
      MH = MVEGL(LHOST)  

      ! *** X-DIRECTION CHANNEL  
      if( MDCHTYP(NMD) == 1 )then  
        MU = 0  
        if( ISVEG >= 1 ) MU = MVEGL(LCHNU)  
        WCHAN = DXP(LCHNU)  
        RLCHN = 0.5*DYP(LCHNU) + CHANLEN(NMD)  
        HCHAN = 0.5*DYP(LCHNU)*H1P(LCHNU) + CHANLEN(NMD)*H1P(LHOST)  
        HCHAN = HCHAN/RLCHN  
        ZBRATUC = 0.5*DYP(LCHNU)*ZBR(LCHNU) + CHANLEN(NMD)*ZBR(LHOST)  
        ZBRATUC = ZBRATUC/RLCHN  
        HURTMP = MAX(ZBRATUC,HCHAN)  
        HUDZBR = HURTMP/ZBRATUC  
        if( HUDZBR < 7.5) HUDZBR = 7.5  
        STBXCH = 0.16/( (LOG( HUDZBR ) -1.)**2)  
        CDMAXU = HCHAN*HCHAN*WCHAN/( DELT*(QCHANU(NMD) + 1.E-12) )  
        STBXCH = MAX(STBXCH,CDMAXU)  
        STBXCH = MAX(STBXCH,0.1)  
        FXVEGCH = 0.0  
        if( MU > 0 ) FXVEGCH = 0.5*( 0.5*DYP(LCHNU)*(BDLPSQ(MU)*H1P(LCHNU)/PVEGZ(MU))  &
                                   + CHANLEN(NMD)*(BDLPSQ(MH)*H1P(LHOST)/PVEGZ(MH)) )/RLCHN  
        CHANFRIC(NMD) = FXVEGCH + STBXCH  
      endif  

      ! *** Y-DIRECTION CHANNEL  
      if( MDCHTYP(NMD) == 2 )then  
        MV = 0  
        if( ISVEG >= 1 ) MV = MVEGL(LCHNV)  
        WCHAN = DYP(LCHNV)  
        RLCHN = 0.5*DXP(LCHNV) + CHANLEN(NMD)  
        HCHAN = 0.5*DXP(LCHNV)*H1P(LCHNV) + CHANLEN(NMD)*H1P(LHOST)  
        HCHAN = HCHAN/RLCHN  
        ZBRATVC = 0.5*DXP(LCHNV)*ZBR(LCHNV) + CHANLEN(NMD)*ZBR(LHOST)  
        ZBRATVC = ZBRATVC/RLCHN  
        HVRTMP = MAX(ZBRATVC,HCHAN)  
        HVDZBR = HVRTMP/ZBRATVC  
        if( HVDZBR < 7.5) HVDZBR = 7.5  
        STBYCH = 0.16/( (LOG( HVDZBR ) -1.)**2)  
        CDMAXV = HCHAN*HCHAN*WCHAN/( DELT*(QCHANV(NMD) + 1.E-12) )  
        STBYCH = MAX(STBYCH,CDMAXV)  
        STBYCH = MAX(STBYCH,0.1)  
        FYVEGCH = 0.0  
        if( MV > 0 ) FYVEGCH = 0.5*(0.5*DXP(LCHNV)*(BDLPSQ(MV)*H1P(LCHNV)/PVEGZ(MV))  &
                                  + CHANLEN(NMD)*(BDLPSQ(MH)*H1P(LHOST)/PVEGZ(MH)) )/RLCHN  
        CHANFRIC(NMD) = FYVEGCH + STBYCH  
      endif  
    enddo  
  endif  
  
  ! *********************************************************************************
  return 

END  

