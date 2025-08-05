! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE MHKPWRDIS  
  ! *** *******************************************************************C  
  ! ***  SUBROUTINE MHKPWRDIS CALCULATES POWER DISSIPATION FROM MARINE HYDROKINETIC  
  ! ***  DEVICES AS A FUNCTION OF WATER VELOCITY  
  !
  !     04/2010  Bill Arnold and Scott James
  !              Sandia National Laboratories
  !  
  ! *** *******************************************************************C  
  use GLOBAL
  implicit none
  
  real :: THRSTCOEF = 0.0,ZTOP = 0.0,ZBOTTOM = 0.0
  real :: FS_WF,FS_NF,FS_EF,FS_SF,MAXSPD,SPDN,SPDS,SPDE,SPDW,DXMHK
  real :: VELUP,UVEC,VVEC,FXMHKA,FYMHKA,FXSUPA,FYSUPA
  real :: FMHK,FSUP,UATVFACE,VATUFACE,UATVFACEN,VATUFACEE
  real :: BLOCKAGE_RATIO = 1.0
  real :: LAYFRACM(KC),LAYFRACS(KC),FLOWSPEED(KC)
  real :: AWEIGHTXW,AWEIGHTXE,AWEIGHTYS,AWEIGHTYN,SUMLAYM,SUMLAYS

  real, save, allocatable :: FXTEMP(:,:)
  real, save, allocatable :: FYTEMP(:,:)
  
  integer :: MHKCOUNT,M,L,K,LW,LE,LS,LN,LNW,LSE,LNE,LP
  
  logical :: STATUS
  logical,save :: FIRSTTIME = .FALSE.

  if( .not. allocated(FXTEMP) )then
    allocate(FXTEMP(LCM,LCM))  
    allocate(FYTEMP(LCM,LCM))  
    FXTEMP = 0.0
    FYTEMP = 0.0
  endif
  
  ! *** *******************************************************************C 
  INQUIRE(FILE = OUTDIR//'POWERBUG.DAT',EXIST = STATUS)
  if( STATUS .and. DEBUG .and.  .not. FIRSTTIME )then
    open(357,FILE = OUTDIR//'POWERBUG.DAT') !this file is for debugging purposes
    FIRSTTIME = .TRUE.
  endif
  
  ! *** Initializations
  MHKCOUNT = 0        ! Running count of the number of cells with MHKdevices (will equal TCOUNT if no cells dry)
  VELUP = 0.0         ! Flow speeds

  do LP = 1,TCOUNT  !loop over the MHK cells
    L = IJLTURB(LP,3)
    LE = LEC(L)         ! east cell, I+1,J
    LN = LNC(L)         ! north cell, I,J+1
    FLOWSPEED(:) = 0.0  ! Flow speeds
    FXMHKE(L) = 0.0     ! x-force from MHK device for the external solution
    FXMHKE(LE) = 0.0    ! east face
    FYMHKE(L) = 0.0     ! y-force from MHK device for the external solution
    FYMHKE(LN) = 0.0    ! north face
    FXSUPE(L) = 0.0     ! x-force from MHK support for the external solution
    FXSUPE(LE) = 0.0    ! east face
    FYSUPE(L) = 0.0     ! y-force from MHK support for the external solution
    FYSUPE(LN) = 0.0    ! north face
    FXMHK(L,KSZ(L):KC) = 0.0    ! x-force from the MHK device for the internal solution
    FXMHK(LE,KSZ(LE):KC) = 0.0  ! east face
    FYMHK(L,KSZ(L):KC) = 0.0    ! y-force from the MHK device for the internal solution
    FYMHK(LN,KSZ(LN):KC) = 0.0  ! north face
    FXSUP(L,KSZ(L):KC) = 0.0    ! x-force from the MHK support for the internal solution
    FXSUP(LE,KSZ(LE):KC) = 0.0  ! east face
    FYSUP(L,KSZ(L):KC) = 0.0    ! y-force from the MHK support for the internal solution
    FYSUP(LN,KSZ(LN):KC) = 0.0  ! north face
    PMHK(L,KSZ(L):KC) = 0.0     ! Power extracted from the MHK device from flow
    PSUP(L,KSZ(L):KC) = 0.0     ! Power extracted from the MHK support from flow
  enddo
  
  do LP = 1,TCOUNT  !loop over the MHK cells
    L = IJLTURB(LP,3)
    if( .not. LMASKDRY(L) ) CYCLE  ! if the cell is dry skip this cell
    MHKCOUNT = MHKCOUNT+1
    LW = LWC(L)         ! west cell, I-1,J
    LE = LEC(L)         ! east cell, I+1,J
    LS = LSC(L)         ! south cell, I,J-1
    LN = LNC(L)         ! north cell, I,J+1
    LNW = LNWC(L)       ! northwest cell, I-1,J+1
    LNE = LNEC(L)       ! northeast cell, I+1,J+1
    LSE = LSEC(L)       ! southeast cell, I+1,J-1
    M = MVEGL(L)-90     ! This was put in to have MHKs be vegetative inputs > 90; this is the mhktype
    LAYFRACM(:) = 0.0   ! initialize the layer fraction variable MHK
    LAYFRACS(:) = 0.0   ! initialize the layer fraction variable support
    
    do K = KSZ(L),KC       ! MHK device layer filler - which layers does the device occupy and at what fraction
      ZTOP = HP(L)*Z(L,K)+BELV(L)     ! layer top elevation
      if( ZTOP < ZMINMHK(M,L) ) CYCLE  ! layer is below device
      ZBOTTOM = HP(L)*Z(L,K-1)+BELV(L)  ! layer bottom elevation
      if( ZBOTTOM > ZMAXMHK(M,L) ) CYCLE !layer is above device
      if( ZTOP >= ZMAXMHK(M,L) .and. ZBOTTOM <= ZMINMHK(M,L) )then  ! device is wholly contained in this layer (special case)
        LAYFRACM(K) = (ZMAXMHK(M,L)-ZMINMHK(M,L))/(HP(L)*DZC(L,K))    ! calculate fraction of layer that is occupied
        EXIT
      endif
      if( ZMAXMHK(M,L) >= ZTOP .and. ZMINMHK(M,L) <= ZBOTTOM )then  ! this layer is fully occupied by the device
        LAYFRACM(K) = 1.0
        CYCLE
      endif
      if( ZBOTTOM<ZMINMHK(M,L) .and. ZMAXMHK(M,L) >= ZTOP )then     ! this layer is partially occupied by the device (bottom)
        LAYFRACM(K) = (ZTOP-ZMINMHK(M,L))/(HP(L)*DZC(L,K))            ! calculate the fraction of layer that is occupied
        CYCLE
      endif
      if( ZTOP >= ZMAXMHK(M,L) .and. ZMINMHK(M,L)<ZBOTTOM )then    ! this layer is partially occupied by the device (top)
        LAYFRACM(K) = (ZMAXMHK(M,L)-ZBOTTOM)/(HP(L)*DZC(L,K))        ! calculate the fraction of layer that is occupied
        CYCLE
      endif
    enddo
    
    SUMLAYM = SUM(LAYFRACM(KSZ(L):KC));SUMLAYS = SUM(LAYFRACS(KSZ(L):KC))        ! Sum of MHK/SUP layer fractions
    if( SUMLAYM == 0.0 .and. CTMHK(M) /= 0 )PRINT*,"Check MHK turbine parameters, looks like the turbine is missing."
    
    do K = KSZ(L),KC                                                 ! MHK support layer filler - which layers does the support occupy and at what fraction
      ZTOP = HP(L)*Z(L,K)+BELV(L)                                    ! layer top elevation
      if( ZTOP<ZMINSUP(M,L)) CYCLE                                  ! layer is below support
      ZBOTTOM = HP(L)*Z(L,K-1)+BELV(L)                               ! layer bottom elevation
      if( ZBOTTOM>ZMAXSUP(M,L)) CYCLE                               ! layer is above support
      if( ZTOP >= ZMAXSUP(M,L) .and. ZBOTTOM <= ZMINSUP(M,L) )then ! support is wholly contained in this layer (special case)
        LAYFRACS(K) = (ZMAXSUP(M,L)-ZMINSUP(M,L))/(HP(L)*DZC(L,K))   ! calculate fraction of layer that is occupied
        EXIT
      endif
      if( ZMAXSUP(M,L) >= ZTOP .and. ZMINSUP(M,L) <= ZBOTTOM )then ! this layer is fully occupied by the support
        LAYFRACS(K) = 1.0
        CYCLE
      endif
      if( ZBOTTOM<ZMINSUP(M,L) .and. ZMAXSUP(M,L) >= ZTOP )then    ! this layer is partially occupied by the support (bottom)
        LAYFRACS(K) = (ZTOP-ZMINSUP(M,L))/(HP(L)*DZC(L,K))           ! calculate the fraction of layer that is occupied
        CYCLE
      endif
      if( ZTOP >= ZMAXSUP(M,L) .and. ZMINSUP(M,L)<ZBOTTOM )then    ! this layer is partially occupied by the support (top)
        LAYFRACS(K) = (ZMAXSUP(M,L)-ZBOTTOM)/(HP(L)*DZC(L,K))        ! calculate the fraction of layer that is occupied
        CYCLE
      endif
    enddo
    
    if( ZMAXMHK(M,L) > (HP(L)+BELV(L)) )then                     ! Check if protruding from water
      PRINT*,'MHK DEVICE PROTRUDING FROM WATER (M, HP, ZMAXMHK, L, I ,J): ',M,HP(L),ZMAXMHK(M,L),L,IL(L),JL(L)
      STOP
    elseif( ZMINMHK(M,L) < BELV(L) )then                           ! Check if below bed elevation?
      PRINT*,'MKH DEVICE IN SEDIMENT (M, BELV, ZMINMHK, L, I ,J): ',M,BELV(L),ZMINMHK(M,L),L,IL(L),JL(L)
      STOP
    endif

    do K = KSZ(L),KC
      UVEC = 0.5*(U(L,K)+U(LE,K))                                    ! I,J cell center u-velocity
      VVEC = 0.5*(V(L,K)+V(LN,K))                                    ! I,J cell center v-velocity
      FLOWSPEED(K) = SQRT(UVEC*UVEC+VVEC*VVEC)                       ! I,J cell center speed
      if( (LAYFRACM(K) == 0.0 .and. LAYFRACS(K) == 0.0) .or. FLOWSPEED(K)<1.0E-03) CYCLE  ! no MHK or support or velocity in this layer
      FMHK = 0.0;FSUP = 0.0 !initialize variables

      UATVFACE  = 0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))   ! u-velocity at south face (the v-face)
      VATUFACE  = 0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))   ! v-velocity at west face  (the u-face)
      UATVFACEN = 0.25*(U(L,K)+U(LE,K)+U(LN,K)+U(LNE,K))   ! u-velocity at north face (the u-north-face)
      VATUFACEE = 0.25*(V(L,K)+V(LE,K)+V(LN,K)+V(LNE,K))   ! v-velocity at east face  (the v-east-face)
      FS_WF = U(L,K); FS_EF = U(LE,K)                          ! velocities on the west/east   faces (u velocities into the cell)
      FS_SF = V(L,K); FS_NF = V(LN,K)                          ! velocities on the south/north faces (v velocities into the cell) !Bug found! Had FS_SF = V(LS,K)
      SPDN = SQRT(UATVFACEN*UATVFACEN+V(LN,K)*V(LN,K))     ! speed at north face
      SPDS = SQRT(UATVFACE*UATVFACE+V(L,K)*V(L,K))         ! speed at south face
      SPDE = SQRT(U(LE,K)*U(LE,K)+VATUFACEE*VATUFACEE)     ! speed at east face
      SPDW = SQRT(U(L,K)*U(L,K)+VATUFACE*VATUFACE)         ! speed at west face
      if( FS_NF > -0.01 )SPDN = 0.0                          ! flow is OUT of north face
      if( FS_SF <  0.01 )SPDS = 0.0                          ! flow is OUT of south face
      if( FS_WF <  0.01 )SPDW = 0.0                          ! flow is OUT of west face
      if( FS_EF > -0.01 )SPDE = 0.0                          ! flow is OUT of east face
      MAXSPD = MAX(SPDN,SPDS,SPDE,SPDW)                      ! identify maximum speed

      if( UPSTREAM == 1 )then         ! use the upstream flowspeed to assess power extraction (not typical)
        ! what face is it on?
        if( MAXSPD == SPDN )then      ! North      
          VELUP = SQRT((0.25*(U(LN,K)+U(LNE,K)+U(LNC(LN),K)+U(LNC(LN)+1,K)))**2+V(LNC(LN),K)**2) 
          DXMHK = DXP(L)
        elseif( MAXSPD == SPDS )then  ! South
          VELUP = SQRT((0.25*(U(LS,K)+U(LSE,K)+U(LSC(LS),K)+U(LSC(LS)+1,K)))**2+V(LS,K)**2)
          DXMHK = DXP(L)
        elseif( MAXSPD == SPDE )then  ! East
          VELUP = SQRT(U(LE+1,K)**2+(0.25*(V(LE,K)+V(LE+1,K)+V(LNE,K)+V(MIN(LC,LN+2),K)))**2)
          DXMHK = DYP(L)
        else                          ! West
          VELUP = SQRT(U(LW  ,K)**2+(0.25*(V(LW,K)+V(LW-1,K)+V(LNW,K)+V(LN-2,K)))**2)
          DXMHK = DYP(L)
        endif
      elseif( UPSTREAM == 0 )then                ! use the local cell's flowspeed to assess power extraction
        VELUP = FLOWSPEED(K)
        ! what face is it on?
        if(MAXSPD==SPDN.or.MAXSPD==SPDS )then     ! North-South
          DXMHK = DYP(L)
        elseif(MAXSPD==SPDE.or.MAXSPD==SPDW )then ! East-West
          DXMHK = DXP(L)
        endif
      endif
      if( LAYFRACM(K) > 0.0 )then     ! MHK device exists in this layer 
        if( VELUP<VMINCUT(M) )then    ! no power generation
          THRSTCOEF = 0.0               ! no need for these calcs
        elseif(VELUP<VMAXCUT(M) )then  ! optimal power generation
          THRSTCOEF = 4.0*(1.0-SQRT(1.0-CTMHK(M)))/(1.0+SQRT(1.0-CTMHK(M))) !From Roc paper
        else                          ! superoptimal flow speed limits power generation to VMAXCUT
          THRSTCOEF = 4.0*(1.0-SQRT(1.0-CTMHK(M)))/(1.0+SQRT(1.0-CTMHK(M)))
          VELUP = VMAXCUT(M)
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   VELUP = U(LWC(L)0,K)  !Special case for calibration for flow from west to east
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! FMHK = 0.5*ThrustCoef*Area*(U_inf)^2 where U_inf is the upstream velocity, VELUP [m^4/s^2]
        FMHK = (B(L,K)+1.0)/DXMHK*0.5*THRSTCOEF*VELUP*VELUP*HP(L)*DZC(L,K)*WIDTHMHK(M)*LAYFRACM(K)*DENMHK(M)  ! area is ASSUMED square
        ! PMHK = FMHK*U where U is the local flowspeed
        PMHK(L,K) = FMHK*FLOWSPEED(K) !ThrustCoef*|u|u^2*area [m^5/s^3] (will yield different power outputs depending on UPSTREAM)
        AWEIGHTXW = DYU(L)*HU(L)/(DYU(L)*HU(L)+DYU(LE)*HU(LE));AWEIGHTXE = 1.0-AWEIGHTXW   ! area-weight for west/east faces
        AWEIGHTYS = DXV(L)*HV(L)/(DXV(L)*HV(L)+DXV(LN)*HV(LN));AWEIGHTYN = 1.0-AWEIGHTYS   ! area-weight for south/north faces
        
        ! To get the x and y components, multiply by a velocity vector divided by the local flow speed FXMHK = FMHK*UVEC/FLOWSPEED(K)
        FXMHK(L ,K) = FXMHK(L ,K) + AWEIGHTXW*SUB(L )*FMHK*UVEC/FLOWSPEED(K)               ! SUB(L)*FMHK(L,K)*(Uvel/q) [m^4/s^2]
        FXMHK(LE,K) = FXMHK(LE,K) + AWEIGHTXE*SUB(LE)*FMHK*UVEC/FLOWSPEED(K)               ! distribute forces on each U-face of the cell
        FYMHK(L ,K) = FYMHK(L ,K) + AWEIGHTYS*SVB(L )*FMHK*VVEC/FLOWSPEED(K)               ! y components of "forces" [m^4/s^2]
        FYMHK(LN,K) = FYMHK(LN,K) + AWEIGHTYN*SVB(LN)*FMHK*VVEC/FLOWSPEED(K)               ! distribute forces on each V-face of the cell
        PMHK(L,K) = PMHK(L,K)*(B(L,K)+1.0)*1000.0
        
        if( DEBUG .and. .not. UPSTREAM .and. MOD(N,100) == 0 )WRITE(357,'(I5,1X,3(I3,1X),5(E12.5,1X))')N,IL(L),JL(L),K,PMHK(L,K),VELUP,UVEC,VVEC,HP(L)
        !IF( DEBUG .and.       UPSTREAM .and. MOD(N,100) == 0 )WRITE(357,'(I5,1X,3(I3,1X),5(E12.5,1X))')N,IL(L),JL(L),K,PMHK(L,K),VELUP,UVEC,VVEC,HP(L)
        !IF( DEBUG .and. M == 1 .and. K == 3 .and. MOD(N,1) == 0 )WRITE(357,'(I5,1X,3(I3,1X),5(E12.5,1X))')N,IL(L),JL(L),K,PMHK(L,K),VELUP,FXMHK(L ,K),FYMHK(L ,K),HP(L)
      endif
      if( LAYFRACS(K) > 0.0 )then                                                     ! MHK support exists in this layer
        FSUP = (B(L,K)+1.0)/DXMHK*0.5*LAYFRACS(K)*CDSUP(M)*FLOWSPEED(K)*FLOWSPEED(K)*HP(L)*DZC(L,K)*WIDTHSUP(M)*DENMHK(M)  ! calculate the force on a cell [m^4/s^2]
        PSUP(L,K) = FSUP*FLOWSPEED(K) !0.5*C_d*Asup*u^3 [m^5/s^3]
        AWEIGHTXW = DYU(L)*HU(L)/(DYU(L)*HU(L)+DYU(LE)*HU(LE));AWEIGHTXE = 1.0-AWEIGHTXW  ! area-weight for west/east face
        AWEIGHTYS = DXV(L)*HV(L)/(DXV(L)*HV(L)+DXV(LN)*HV(LN));AWEIGHTYN = 1.0-AWEIGHTYS  ! area-weight for south/north face
        FXSUP(L ,K) = FXSUP(L ,K)+AWEIGHTXW*SUB(L )*FSUP*UVEC/FLOWSPEED(K)              ! x-component [m^4/s^2] Note that FLOWSPEED(K) is multiplied in but divided back out again when normalizing by UVEC/FLOWSPEED(K)
        FXSUP(LE,K) = FXSUP(LE,K)+AWEIGHTXE*SUB(LE)*FSUP*UVEC/FLOWSPEED(K)              ! distribute forces on both U-faces of the cell 
        FYSUP(L ,K) = FYSUP(L ,K)+AWEIGHTYS*SVB(L )*FSUP*VVEC/FLOWSPEED(K)              ! y-component [m^4/s^2] Note that FLOWSPEED(K) is multiplied in but divided back out again when normalizing by VVEC/FLOWSPEED(K)
        FYSUP(LN,K) = FYSUP(LN,K)+AWEIGHTYN*SVB(LN)*FSUP*VVEC/FLOWSPEED(K)              ! distribute forces on both V-faces of the cell
        PSUP(L,K) = PSUP(L,K)*(B(L,K)+1.0)*1000.0
      endif
      if( VELUP < 1.0E-3 ) CYCLE         ! avoid divide by zero errors
      FXMHKE(L) = FXMHKE(L)+ABS(FXMHK(L,K))/VELUP;FXMHKE(LE) = FXMHKE(LE)+ABS(FXMHK(LE,K))/VELUP  ! Sum layer force magnitudes for external mode solution (need absolute value of forces)
      FYMHKE(L) = FYMHKE(L)+ABS(FYMHK(L,K))/VELUP;FYMHKE(LN) = FYMHKE(LN)+ABS(FYMHK(LN,K))/VELUP  ! Sum layer force magnitudes for external mode solution (these forces are later multiplied by a directional velocity so they need to be absolute values here)
      if(FLOWSPEED(K)<1.0E-3) CYCLE  ! avoid divide by zero errors
      FXSUPE(L) = FXSUPE(L)+ABS(FXSUP(L,K))/FLOWSPEED(K);FXSUPE(LE) = FXSUPE(LE)+ABS(FXSUP(LE,K))/FLOWSPEED(K)  ! Sum layer force magnitudes for external mode solution (absolute values because these are later multiplied by the local velocity to apply a direction)
      FYSUPE(L) = FYSUPE(L)+ABS(FYSUP(L,K))/FLOWSPEED(K);FYSUPE(LN) = FYSUPE(LN)+ABS(FYSUP(LN,K))/FLOWSPEED(K)  ! Sum layer force magnitudes for external mode solution
    enddo
    
    EMHK(MHKCOUNT,L) = EMHK(MHKCOUNT,L) + DT*SUM(PMHK(L,KSZ(L):KC))*2.7778E-10  ! factor converts to MW-hr
    ESUP(MHKCOUNT,L) = ESUP(MHKCOUNT,L) + DT*SUM(PSUP(L,KSZ(L):KC))*2.7778E-10  ! factor converts to MW-hr
  enddo

  ! CALEXP2T is expecting units of [1/s] for FXMHKE, which is the sum of absolute values of FXMHK divided by the average speed in the cell
  ! CALEXP2T divides by water-column volume before passing this "force" onto FUHDYE (in units of [1/s]), which is used for momentum conservation in CALPUV
  ! Units of FXMHKE (etc) were the same as FX and FXMHK (etc) [m^4/s^2] so they must be divided by the average speed in this water column
  
  do LP = 1,TCOUNT  !loop over the MHK cells
    L = IJLTURB(LP,3)
    M = MVEGL(L)-90
    if( PB_COEF == 0.0 )then
        BLOCKAGE_RATIO = 1.0 + (ZMAXMHK(M,L)-ZMINMHK(M,L))/HP(L)  ! This acts as a force multiplier to take the blockage ratio into account. If it is specified as nonzero in MHK.INP, then PB_COEF is used, otherwise it is approximated as the fraction of the water column that is occupied by the turbine
    else
        BLOCKAGE_RATIO = PB_COEF
    endif

    !SUM(FXMHK(L,1:KC)*DZC(1:KC))calculate the AVGERAGE x-force applied by the MHK/SUP against the flow, multiplying it by KC below approximates TOTAL force
    !SUM(FYMHK(L,1:KC)*DZC(1:KC))calculate the AVGERAGE y-force applied by the MHK/SUP against the flow, multiplying it by KC below approximates TOTAL force
    FXMHKA = SUM(FXMHK(L,KSZ(L):KC)*DZC(L,KSZ(L):KC)) 
    FYMHKA = SUM(FYMHK(L,KSZ(L):KC)*DZC(L,KSZ(L):KC))
    FXSUPA = SUM(FXSUP(L,KSZ(L):KC)*DZC(L,KSZ(L):KC))
    FYSUPA = SUM(FYSUP(L,KSZ(L):KC)*DZC(L,KSZ(L):KC))
    do K = KSZ(L),KC
       FXTEMP(L,K) = BLOCKAGE_RATIO*( FXMHK(L,K)-FXMHKA + FXSUP(L,K)-FXSUPA )*FLOAT(KC-KSZ(L)+1)  ! Pull x-force out of MHK/support layer for internal mode - push forces in other layers
       FYTEMP(L,K) = BLOCKAGE_RATIO*( FYMHK(L,K)-FYMHKA + FYSUP(L,K)-FYSUPA )*FLOAT(KC-KSZ(L)+1)  ! Pull y-force out of MHK/support layer for internal mode - push forces in other layers
    enddo       

    ! *** Treat opposite face
    LE = LEC(L)
    if( MVEGL(LE) < 90 )then
      FXMHKA = SUM(FXMHK(LE,KSZ(LE):KC)*DZC(LE,KSZ(LE):KC))
      FXSUPA = SUM(FXSUP(LE,KSZ(LE):KC)*DZC(LE,KSZ(LE):KC))
      do K = KSZ(LE),KC
         FXTEMP(LE,K) = BLOCKAGE_RATIO*( FXMHK(LE,K)-FXMHKA + FXSUP(LE,K)-FXSUPA )*FLOAT(KC-KSZ(LE)+1)  ! Pull x-force out of MHK/support layer for internal mode - push forces in other layers
      enddo       
    endif
    
    LN = LNC(L)
    if( MVEGL(LN) < 90 )then
      FYMHKA = SUM(FYMHK(LN,KSZ(LN):KC)*DZC(LN,KSZ(LN):KC))
      FYSUPA = SUM(FYSUP(LN,KSZ(LN):KC)*DZC(LN,KSZ(LN):KC))
      do K = KSZ(LN),KC
         FYTEMP(LN,K) = BLOCKAGE_RATIO*( FYMHK(LN,K)-FYMHKA + FYSUP(LN,K)-FYSUPA )*FLOAT(KC-KSZ(LN)+1)  ! Pull y-force out of MHK/support layer for internal mode - push forces in other layers
      enddo       
    endif
  enddo

  do LP = 1,TCOUNT  !loop over the MHK cells
    L = IJLTURB(LP,3)
    FX(L,KSZ(L):KC) = FX(L,KSZ(L):KC) + FXTEMP(L,KSZ(L):KC)
    FY(L,KSZ(L):KC) = FY(L,KSZ(L):KC) + FYTEMP(L,KSZ(L):KC)
    
    FXMHKE(L) = FXMHKE(L)*DXYIU(L)*HUI(L)  !external mode solution units of [1/s] divide by water-column volume: Turbine X base calculations on U face so DXYIP --> DXYIU and HPI --> HUI
    FYMHKE(L) = FYMHKE(L)*DXYIV(L)*HVI(L)  !external mode solution units of [1/s] divide by water-column volume: Turbine Y base calculations on U face so DXYIP --> DXYIV and HPI --> HVI
    FXSUPE(L) = FXSUPE(L)*DXYIU(L)*HUI(L)  !external mode solution units of [1/s] divide by water-column volume: Support X
    FYSUPE(L) = FYSUPE(L)*DXYIV(L)*HVI(L)  !external mode solution units of [1/s] divide by water-column volume: Support Y

    ! *** Treat opposite face
    LE = LEC(L)
    if( MVEGL(LE) < 90 )then
      FX(LE,KSZ(LE):KC) = FX(LE,KSZ(LE):KC) + FXTEMP(LE,KSZ(LE):KC)
      FXMHKE(LE) = FXMHKE(LE)*DXYIU(LE)*HUI(LE)  !external mode solution units of [1/s] divide by water-column volume: Turbine X
      FXSUPE(LE) = FXSUPE(LE)*DXYIU(LE)*HUI(LE)  !external mode solution units of [1/s] divide by water-column volume: Support X
    endif
    
    LN = LNC(L)
    if( MVEGL(LN) < 90 )then
      FY(LN,KSZ(LN):KC) = FY(LN,KSZ(LN):KC) + FYTEMP(LN,KSZ(LN):KC)
      FYMHKE(LN) = FYMHKE(LN)*DXYIV(LN)*HVI(LN)  !external mode solution units of [1/s] divide by water-column volume: Turbine Y
      FYSUPE(LN) = FYSUPE(LN)*DXYIV(LN)*HVI(LN)  !external mode solution units of [1/s] divide by water-column volume: Support Y
    endif
  enddo
  
  return
  
END SUBROUTINE MHKPWRDIS
