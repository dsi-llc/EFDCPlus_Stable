! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE BEDLOADJ
  !  Bedload transport subroutine based on Van Rijn's transport
  !  Equations.
  !
  !  University of California, Santa Barbara
  !  Craig Jones and Wilbert Lick
  
  ! ORIGINAL:  May 24, 2006
  !  Craig Jones and Scott James
  ! REVISED: Added toxics linkage and bedload mass balance updates - 2017-01-04
  !  Paul M. Craig
  ! REVISED: SIGMA-ZED AND OMP - 2016-11-07
  !  Paul M. Craig
  ! REVISED: Added bybass by class - 2020-01-12
  !  Paul M. Craig
  ! REVISED: Eliminated use of NNONCO.  Always check D50 - 2020-09-23
  !  Paul M. Craig
  
  use GLOBAL    
  !use MPI
  use Variables_MPI
  use Communicate_Ghost_Routines

  implicit none 
  
  integer :: L, NS, NT, LW, LS, LE, LN, ND, LF, LL, LP, IERR
  integer :: ISKIP(NSEDS)
  
  real(RKD) :: UTMP, VTMP 

  real(RKD) ,save,allocatable,dimension(:) :: DXUCM
  real(RKD) ,save,allocatable,dimension(:) :: DYVCM
  real(RKD) ,save,allocatable,dimension(:) :: DXYIPCM

  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, TWAIT          ! *** MODEL TIMING TEMPORARY VARIABLES

  if( .not. allocated(DXUCM) )then
    allocate(DXUCM(LCM))
    allocate(DYVCM(LCM))
    allocate(DXYIPCM(LCM))
    DXUCM(2:LA) = DBLE(DXYP(2:LA))/DBLE(DXP(2:LA))*100._8
    DYVCM(2:LA) = DBLE(DXYP(2:LA))/DBLE(DYP(2:LA))*100._8
    DXYIPCM(2:LA) = DBLE(DXYIP(2:LA))/10000._8
    
    PSUS = 1.0   ! *** Default to any eroded cohesive sediments into the water column.  Logic below for non-cohesives
  endif
  ISKIP = 1
  
  !$OMP PARALLEL DEFAULT(SHARED)
  
  ! *********************************************************************
  ! *** Setup for Toxics Transport using previous timestep CBL and CBLTOX
  if( ISTRAN(5) > 0 )then
    !$OMP DO PRIVATE(ND,LF,LL,NT,NS,LP,L)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      do NT = 1,NTOX
        do NS = 1,NSEDS
          if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
          do LP = LF,LL
            L = LSED(LP)
            if( CBL(L,NS) > 1.E-6 )then
              ! *** MG/G             MG/M2      CM2/G    M2/CM2 
              CBLTXCON(L,NS,NT) = CBLTOX(L,NT)/CBL(L,NS)/10000.
            else
              CBLTXCON(L,NS,NT) = 0.0
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
  endif
  
  ! *** ************************************************************
  ! *** Calculate Percentage of erosion into suspension PSUS
  ! *** and whether the cell has bedload or not BLFLAG
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,NS)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)
    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
      do LP = LF,LL
        L = LSED(LP)
        USW(L,NS) = SQRT(TAU(L))/DWS(NS)   ! *** DWS is settling speed.  USW is the shear velocity
      enddo
    enddo
  enddo
  !$OMP END DO
  
  ! *** Loop below determines if bedload exists or not.  There are three regimes of tranport in this loop.
  ! *** the first conditional in the where check for a large enough particle's diameter and small enough
  ! *** shear velocity.  If the particle is too small or the shear velocity is too large, then all the sediment 
  ! *** transport is in the suspended load, specified by suspended probability (PSUS = 1).  If the particle
  ! *** is large enough and the shear velocity is small enough then we have two situations. In the first case, 
  ! *** shear stress tau is smaller than the critical shear velocity or if shear velocity is negative or zero
  ! *** then there is neither bedload transport or suspended load transport.  Otherwise, both bedload and suspended
  ! *** load transport exists.  Also calculated is the probability of suspension for suspended load PSUS (eqn. 8).
  !$OMP DO PRIVATE(ND,LF,LL,NS,LP,L)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)
    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
      do LP = LF,LL
        L = LSED(LP)
        ! ***                     
        if( USW(L,NS) < 4.0 )then                                 ! *** "USW(L,NS) < 4.0" is an out of range check
          if( TAU(L) <= TCRE(NS) .or. USW(L,NS) <= 0.0 )then
            ! *** Shear is too small to erode anything for current class
            PSUS(L,NS) = 0.0
          else      
            ISKIP(NS) = 0
            PSUS(L,NS) = max((LOG(USW(L,NS))-LOG(SQRT(TCRSUS(NS))/DWS(NS)))/(LOG(4.0)-LOG(SQRT(TCRSUS(NS))/DWS(NS))),0.0)
          endif
        else
          ! *** Shear is high enough to move all eroded material to suspension
          PSUS(L,NS) = 1.0
        endif 
      enddo
    enddo  
  enddo
  !$OMP END DO

  ! *** Compute the bedload velocities (cm/s) in the U and V directions
  !$OMP DO PRIVATE(ND,LF,LL,NS,LP,L)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)
    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
      if( ISKIP(NS) == 0 )then
        do LP = LF,LL
          L = LSED(LP)
          TRANS(L,NS) = max((TAU(L)-TCRE(NS))/TCRE(NS),0.0)                     ! *** eqn. 21
          DZBL(L,NS)  = D50(NS)/10000.0*0.3*DISTAR(NS)**0.7*SQRT(TRANS(L,NS))   ! *** eqn. 20b
          DZBL(L,NS)  = min(DZBL(L,NS), HPCM(L))                                ! *** Don't allow bedload height to exceed water column depth
          BLVEL(L,NS) = 1.5*TRANS(L,NS)**0.6*SQRT(((SEDDENS(NCORENO(IL(L),JL(L)))/WATERDENS) -1.0)*980.0*D50(NS)/10000.0)    ! *** eqn. 20a  (cm/s)
        enddo
      endif
    enddo
  enddo
  !$OMP END DO
      
  !$OMP DO PRIVATE(ND,LF,LL,NS,LP,L,LS,LW)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)
    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
      if( ISKIP(NS) == 0 )then
        do LP = LF,LL
          ! *** Interpolate BLVEL onto faces
          L = LSED(LP)
          LS = LSC(L)
          LW = LWC(L)
          UBL(L,NS) = 0.5*SUB(L)*(BLVEL(L,NS)*UCELLCTR(L) + BLVEL(LW,NS)*UCELLCTR(LW))
          VBL(L,NS) = 0.5*SVB(L)*(BLVEL(L,NS)*VCELLCTR(L) + BLVEL(LS,NS)*VCELLCTR(LS))
        enddo
      endif
    enddo
  enddo
  !$OMP END DO
 
  if( ISSLOPE )then                                            ! *** if bedslope is calculated
    !$OMP DO PRIVATE(ND,LF,LL,NS,LP,L,UTMP,VTMP)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      do NS = 1,NSEDS
        if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
        if( ISKIP(NS) == 0 )then
          do LP = LF,LL
            L = LSED(LP)
            UTMP = UBL(L,NS)                                     ! *** save original x-bedload velocity
            VTMP = VBL(L,NS)                                     ! *** save original x-bedload velocity
            UBL(L,NS) = ALPHA_PX(L)*UBL(L,NS)                    ! *** modify by pitch angle
            VBL(L,NS) = ALPHA_PY(L)*VBL(L,NS)                    ! *** modify by roll angle
            if( UBL(L,NS)>VBL(L,NS) )then                        ! *** find dominant velocity direction
              VBL(L,NS) = VBL(L,NS) + ALPHA_RX(L,NS)*UTMP        ! *** Bedload velocity (x is dominant roll rirection)
              UBL(L,NS) = UBL(L,NS) - ALPHA_RX(L,NS)*VBL(L,NS)   ! *** as impacted by bedslope
              UBL(L,NS) = UBL(L,NS) + ALPHA_RY(L,NS)*VTMP        ! *** (secondary roll due to y)
              VBL(L,NS) = VBL(L,NS) - ALPHA_RY(L,NS)*UBL(L,NS)   ! *** see Lesser (2004) Ikeda (1982)
            else
              UBL(L,NS) = UBL(L,NS) + ALPHA_RY(L,NS)*VTMP        ! *** Bedload velocity (y is dominant roll direction)
              VBL(L,NS) = VBL(L,NS) - ALPHA_RY(L,NS)*UBL(L,NS)   ! *** as impacted by bedslope
              VBL(L,NS) = VBL(L,NS) + ALPHA_RX(L,NS)*UTMP        ! *** (secondary roll due to x)
              UBL(L,NS) = UBL(L,NS) - ALPHA_RX(L,NS)*VBL(L,NS)   ! *** see Lesser (2004) Ikeda (1982)
            endif
          enddo
        endif
      enddo
    enddo
    !$OMP END DO
  endif
     
  ! *** *******************************************************************!
  ! All the equations below are solving the pde in eqn.18.
  !$OMP DO PRIVATE(ND,LF,LL,NS,LP,L,LS,LW)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)

    ! *** CBL     - Bedload concentration (g/cm^2)     (Original SNL was in g/cm^3)
    ! *** QSBDLDX - Bedload flux in X direction (g/s)  (Original SNL was in g/cm^2)
    ! *** QSBDLDY - Bedload flux in Y direction (g/s)  (Original SNL was in g/cm^2)
    ! *** DZBL    - Bedload (i.e. saltation) height (cm)
    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
      if( ISKIP(NS) == 0 )then
        do LP = LF,LL
          L = LSED(LP)
          LS = LSC(L)
          LW = LWC(L)
          
          !  *** X Bedload flux at I-1/2 interface
          if( UBL(L,NS) == 0. )then
            QSBDLDX(L,NS) = 0.
          elseif( UBL(L,NS) > 0. )then
            ! *** g/s                cm           g/cm2    cm/s
            QSBDLDX(L,NS) = SUB(L)*DXUCM(LW)*CBL(LW,NS)*UBL(L,NS)    ! *** FLOWING EAST
          else
            QSBDLDX(L,NS) = SUB(L)*DXUCM(L) *CBL(L,NS) *UBL(L,NS)    ! *** FLOWING WEST
          endif
        
          ! *** Y Bedload flux at J-1/2 interface
          if( VBL(L,NS) == 0. )then
            QSBDLDY(L,NS) = 0.
          elseif( VBL(L,NS) > 0. )then
            QSBDLDY(L,NS) = SVB(L)*DYVCM(LS)*CBL(LS,NS)*VBL(L,NS)    ! *** FLOWING NORTH
          else
            QSBDLDY(L,NS) = SVB(L)*DYVCM(L) *CBL(L,NS) *VBL(L,NS)    ! *** FLOWING SOUTH
          endif
        enddo
      else
        ! *** ZERO THE FLUXES
        do LP = LF,LL
          L = LSED(LP)
          QSBDLDX(L,NS) = 0.0
          QSBDLDY(L,NS) = 0.0
        enddo
      endif
    enddo
  enddo
  !$OMP END DO
     
  !$OMP SINGLE
  TTDS = DSTIME(0)
  call MPI_barrier(DSIcomm, ierr)
  TWAIT = DSTIME(0) - TTDS
  TTSED = TTSED - TWAIT

  TTDS = DSTIME(0)
  call Communicate_BEDLOAD(1,NSEDS)
  DSITIMING(8) = DSITIMING(8) + (DSTIME(0) - TTDS)
  !$OMP END SINGLE

  ! *** *************************************************************************
  ! *** Transport Equation for bedload concentration (g/cm2)
  !$OMP DO PRIVATE(ND,LF,LL,NS,LP,L,LE,LN)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)

    do NS = 1,NSEDS
      if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
      do LP = LF,LL
        L = LSED(LP)
        LE = LEC(L)
        LN = LNC(L)
        CBL(L,NS) = CBL(L,NS) + ( DXYIPCM(L)*DTSEDJ*( QSBDLDX(L,NS)-QSBDLDX(LE,NS) + QSBDLDY(L,NS)-QSBDLDY(LN,NS) ) + QBLFLUX(L,NS) )
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  return

END SUBROUTINE BEDLOADJ
