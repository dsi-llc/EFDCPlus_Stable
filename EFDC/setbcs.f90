! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
!> @details SETS BOUNDARY CONDITION SWITCHES
!
!---------------------------------------------------------------------------!
!> @author
!---------------------------------------------------------------------------!
! CHANGE RECORD
! DATE MODIFIED     BY               DESCRIPTION
!---------------------------------------------------------------------------!
!   2019-10      Zander     Adding calls for mapping global RSSBC* arrays
!   2019-09      Zander     Added mapping routine calls for domain decomp
!   Unknown      Unknown    MODIFIED BOUNDARY CONDITION FLAGS FOR TYPE 2
!                            OPEN BOUNDARIES
!   Unknown      Unknown    ADDED REAL FLAGS RSSBCE(L),RSSBCW(L),
!                           RSSBCN(L),RSSBCS(L)
!   Unknown      Unknown    TO MODIFIED CALCULATION OF CELL CENTER BED STRESS
!                           (STORED AS QQ(L,0)) AND THE OUTPUTTED CELL
!                           CENTER VELOCITY FOR CELLS HAVE SOURCE/SINKS
!
!---------------------------------------------------------------------------!

SUBROUTINE SETBCS

  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Write_Out

  implicit none

  ! *** Local variables
  integer :: L, I, J, NPN, LD, LE, LW, LS, LN, LL, LU, NCTL, IU, JU, IQ, L1, ID, JD, ierr
  integer :: LTMP, NWR, NJP, NMD, NDRYTMP, NC
  real    :: DDYDDDX,DDXDDDY,RQDW
  real,save,allocatable,dimension(:) :: SUBEW
  real,save,allocatable,dimension(:) :: SVBNS

  if( .not. allocated(SUBEW) )then
    allocate(SUBEW(LCM))
    allocate(SVBNS(LCM))
    SUBEW = 0.0
    SVBNS = 0.0
  endif


  ! *** SET LAND-WATER BOUNDARY SWITCHES
  ITRICELL = 0

  do L = 2,LA
    I = IL(L)
    J = JL(L)
    if( LCT(L) == 1 )then
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
    if( LCT(L) == 2 )then
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
    if( LCT(L) == 3 )then
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
    if( LCT(L) == 4 )then
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
    if( LCT(L) == 5 )then
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
    if( LCT(L) == 6 )then
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 6 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 7 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 6 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 7 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
    if( LCT(L) == 7 )then
      if( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      if( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 6 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 7 ) SUB(L) = 1.
      if( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      if( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      if( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 6 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 7 ) SVB(L) = 1.
      if( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    endif
  enddo
  SUB(1) = 0.
  SVB(1) = 0.
  SUB(LC) = 0.
  SVB(LC) = 0.

  ! *** Modify land-water boundary conditions for connectors in E-W direction
  if( ISCONNECT >=  2 )then
    do NPN = 1,NPEWBP
      L = LIJ(IWPEW(NPN),JWPEW(NPN))
      SUB(L) = 1.
      SUBO(L) = 1.
    enddo
  endif

  ! *** Modify land-water boundary conditions for connectors in N-S direction
  if( ISCONNECT == 1 .or. ISCONNECT == 3 )then
    do NPN = 1,NPNSBP
      LS = LIJ(ISPNS(NPN),JSPNS(NPN))
      SVB(LS)  = 1.
      SVBO(LS) = 1.
    enddo
  endif
  
  ! *** New Mapping routines for MPI Domain Decomp
  call Map_River
  call Map_OpenBC_Pressure
  call Map_OpenBC_Conc
  call Map_Hydraulic_Structures
  call Map_Withdrawal_Return
  call Map_Jet_Plume
  
  ! *** SET WATER-WATER (P OR SURFACE ELEVATION) BOUNDARY SWITCHES
  do LL = 1,NPBW
    I = IPBW(LL)
    J = JPBW(LL)
    L = LIJ(I,J)
    LPBW(LL) = L
    LE = LEC(L)
    SWB(L) = 0.               ! *** Used for On/Off of Vertical Velocities
    SPB(L) = 0.               ! *** Used for On/Off of various BC forcings at open boundaries
    
    SUB(L) = 0.
    SVB(L) = 0.               ! *** Bottom of cell
    SVB(LNC(L)) = 0.          ! *** Top of cell
    SAAX(LE) = 0.             ! *** Used for On/Off of Explicit Solution Horizontal Momentum
    SCAX(LE) = 0.             ! *** Used for On/Off of Coriolis Acceleration
    SDX(LE) = 0.              ! *** Used for On/Off of Horizontal Eddy Viscosity
    if( ISPBW(LL) == 0 .or. ISPBW(LL) == 1 .or. ISPBW(LL) == 4 )then
      ! *** Zero Tangiential (East of West BC)
      SWB(LE) = 0.
      SVB(LE) = 0.     
      SVB(LNC(LE)) = 0.
    endif
  enddo
  
  do LL = 1,NPBE
    I = IPBE(LL)
    J = JPBE(LL)
    L = LIJ(I,J)
    LPBE(LL) = L
    LW = LWC(L)
    SWB(L) = 0.               ! *** Used for On/Off of Vertical Velocities
    SPB(L) = 0.               ! *** Used for On/Off of various BC forcings at open boundaries
    
    SVB(L) = 0.               ! *** Bottom of cell
    SVB(LNC(L)) = 0.          ! *** Top of cell
    SAAX(L) = 0.              ! *** Used for On/Off of Explicit Solution Horizontal Momentum
    SCAX(L) = 0.              ! *** Used for On/Off of Coriolis Acceleration
    SDX(L) = 0.               ! *** Used for On/Off of Horizontal Eddy Viscosity
    if( ISPBE(LL) == 0 .or. ISPBE(LL) == 1 .or. ISPBE(LL) == 4 )then
      ! *** Zero Tangiential (West of East BC)
      SWB(LW) = 0.
      SVB(LW) = 0.
      SVB(LNC(LW)) = 0.
    endif
  enddo
  
  do LL = 1,NPBS
    I = IPBS(LL)
    J = JPBS(LL)
    L = LIJ(I,J)
    LPBS(LL) = L
    LN = LNC(L)
    SWB(L) = 0.               ! *** Used for On/Off of Vertical Velocities
    SPB(L) = 0.               ! *** Used for On/Off of various BC forcings at open boundaries
    
    SVB(L) = 0.     
    SUB(L) = 0.               ! *** West of cell  
    SUB(LEC(L)) = 0.          ! *** East of cell
    SAAY(LN) = 0.             ! *** Used for On/Off of Explicit Solution Horizontal Momentum
    SCAY(LN) = 0.             ! *** Used for On/Off of Coriolis Acceleration
    SDY(LN) = 0.              ! *** Used for On/Off of Horizontal Eddy Viscosity
    if( ISPBS(LL) == 0 .or. ISPBS(LL) == 1 .or. ISPBS(LL) == 4 )then
      ! *** Zero tangiential (North of South BC)
      SWB(LN) = 0.
      SUB(LN) = 0.
      SUB(LEC(LN)) = 0.
    endif
  enddo
  
  do LL = 1,NPBN
    I = IPBN(LL)
    J = JPBN(LL)
    L = LIJ(I,J)
    LPBN(LL) = L
    LS = LSC(L)
    SWB(L) = 0.               ! *** Used for On/Off of Vertical Velocities
    SPB(L) = 0.               ! *** Used for On/Off of various BC forcings at open boundaries
    
    SUB(L) = 0.               ! *** West of cell  
    SUB(LEC(L)) = 0.          ! *** East of cell
    SAAY(L) = 0.              ! *** Used for On/Off of Explicit Solution Horizontal Momentum
    SCAY(L) = 0.              ! *** Used for On/Off of Coriolis Acceleration
    SDY(L) = 0.               ! *** Used for On/Off of Horizontal Eddy Viscosity
    if( ISPBN(LL) == 0 .or. ISPBN(LL) == 1 .or. ISPBN(LL) == 4 )then
      ! *** Zero Tangiential (South of North BC)
      SWB(LS) = 0.
      SUB(LS) = 0.
      SUB(LEC(LS)) = 0.
    endif
  enddo

  ! *********************************************************************
  ! *** SET THE CELL FACES SWITCHES FOR HEAD CONTROL STRUCTURES
  ! *** UPSTREAM CONTROL
  do NCTL = 1,NQCTL
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    L = LIJ(IU,JU)
    SAVESUB(1,NCTL) = SUB(L)
    SAVESVB(1,NCTL) = SVB(L)

    ! *** SET U FACE
    LW = LWC(L)
    do NC = 1,NQCTL
      if( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
      if( NC /= NCTL .and. HYD_STR(NC).NQCTYP /=  3 .and. HYD_STR(NC).NQCTYP /=  4 )then
        I = HYD_STR(NC).IQCTLU
        J = HYD_STR(NC).JQCTLU
        L1 = LIJ(I,J)
        if( L1 == LW )then
          SUB(L) = 0.0
          EXIT
        endif
      endif
    enddo

    ! *** SET V FACE
    LS = LSC(L)
    do NC = 1,NQCTL
      if( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
      if( NC /=  NCTL .and. HYD_STR(NC).NQCTYP /=  3 .and. HYD_STR(NC).NQCTYP /=  4 )then
        I = HYD_STR(NC).IQCTLU
        J = HYD_STR(NC).JQCTLU
        L1 = LIJ(I,J)
        if( L1 == LS )then
          SVB(L) = 0.0
          EXIT
        endif
      endif
    enddo
  enddo

  ! *** DOWNSTREAM CONTROL
  do NCTL = 1,NQCTL
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD
    if( ID /=  0 .and. JD /=  0 )then
      L = LIJ(ID,JD)
      SAVESUB(2,NCTL) = SUB(L)
      SAVESVB(2,NCTL) = SVB(L)

      ! *** SET U FACE
      LW = LWC(L)
      do NC = 1,NQCTL
        if( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
        if( NC /=  NCTL .and. HYD_STR(NC).NQCTYP /=  3 .and. HYD_STR(NC).NQCTYP /=  4 )then
          I = HYD_STR(NC).IQCTLD
          J = HYD_STR(NC).JQCTLD
          if( I > 0 .and. J > 0 )then
            L1 = LIJ(I,J)
            if( L1 == LW )then
              SUB(L) = 0.0
              EXIT
            endif
          endif
        endif
      enddo

      ! *** SET V FACE
      LS = LSC(L)
      do NC = 1,NQCTL
        if( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
        if( NC /=  NCTL .and. HYD_STR(NC).NQCTYP /=  3 .and. HYD_STR(NC).NQCTYP /=  4 )then
          I = HYD_STR(NC).IQCTLD
          J = HYD_STR(NC).JQCTLD
          if( I > 0 .and. J > 0 )then
            L1 = LIJ(I,J)
            if( L1 == LS )then
              SVB(L) = 0.0
              EXIT
            endif
          endif
        endif
      enddo
    endif
  enddo

  ! *** RESET DXU,DYU,DXV,DYV BASED ON BOUNDARY CONDITION SWITCHES
  do L = 2,LA
    if( SUB(L) > 0.5 )then
      DXU(L) = 0.5*(DXP(L)+DXP(LWC(L)))
      DYU(L) = 0.5*(DYP(L)+DYP(LWC(L)))
    endif
    if( SUB(L) < 0.5 .and. SUB(LEC(L)) > 0.5 )then
      DXU(L) = DXP(L)
      DDYDDDX = 2.*(DYP(LEC(L))-DYP(L))/(DXP(L)+DXP(LEC(L)))
      DYU(L) = DYP(L)-0.5*DXP(L)*DDYDDDX
    endif
    if( SUB(L) < 0.5 .and. SUB(LEC(L)) < 0.5 )then
      DXU(L) = DXP(L)
      DYU(L) = DYP(L)
    endif
  enddo

  do L = 2,LA
    LN = LNC(L)
    LS = LSC(L)
    if( SVB(L) > 0.5 )then
      DXV(L) = 0.5*(DXP(L)+DXP(LS))
      DYV(L) = 0.5*(DYP(L)+DYP(LS))
    endif
    if( SVB(L) < 0.5 .and. SVB(LN) > 0.5 )then
      DDXDDDY = 2.*(DXP(LN)-DXP(L))/(DYP(L)+DYP(LN))
      DXV(L) = DXP(L)-0.5*DYP(L)*DDXDDDY
      DYV(L) = DYP(L)
    endif
    if( SVB(L) < 0.5 .and. SVB(LN) < 0.5 )then
      DXV(L) = DXP(L)
      DYV(L) = DYP(L)
    endif
  enddo

  ! *** Set thin barriers by calling CELLMASK
  if( ISMASK == 1 )  call CELLMASK

  ! *** Navigation locks - Set initial masks to off
  do NCTL = 1,NQCTL
    LU = LIJ(HYD_STR(NCTL).IQCTLU,HYD_STR(NCTL).JQCTLU)
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD

    if( ID /=  0 .and. JD /=  0 .and. HYD_STR(NCTL).NQCTYP == 9 )then
      ! *** Find the appropriate face and set the masks
      LD = LIJ(ID,JD)
      if( LD == LSC(LU) )then
        SVB(LU)  = 1.0
      elseif( LU == LSC(LD) )then
        SVB(LD)  = 1.0
      elseif( LD == LWC(LU) )then
        SUB(LU)  = 1.0
      elseif( LU == LWC(LD) )then
        SUB(LD)  = 1.0
      endif
    endif
  enddo 
    
  ! *** Set cell-face flag of external density gradient in external/internal mode solution 
  ! *** for two cell wide channels
  IU = 0
  JU = 0
  do L = 2,LA
    I = IL(L)
    J = JL(L)
    LE = LEC(L)
    LW = LWC(L)
    LN = LNC(L)
    LS = LSC(L)
    
    ! *** Two cell flow channel - West face
    if( SUB(L) > 0.5 .and. SUB(LW) < 0.5 )then
      if( SUB(LE) < 0.5 )then
        SUBD(L) = 0.0
        IU = IU + 1
      endif
    endif
    
    ! *** Two cell flow channel - South face
    if( SVB(L) > 0.5 .and. SVB(LS) < 0.5 )then
      if( SVB(LN) < 0.5 )then
        SVBD(L) = 0.0
        JU = JU + 1
      endif
    endif

  enddo
  PRINT "(' EDG Test: Number of West Faces: ',I6,', Number of South Faces:',I6 )", IU, JU
  
  ! *** SET VOLUMETRIC & CONCENTRATION SOURCE LOCATIONS AND BED STRESS
  ! *** AND CELL CENTER BED STRESS AND VELOCITY MODIFERS
  do LL = 1,NQSIJ
    I = BCPS(LL).I
    J = BCPS(LL).J
    LTMP = LIJ(I,J)
    BCPS(LL).L = LTMP
    BCPS(LL).RQSMUL = 1.0                                             ! *** DEFAULT - Flows are already in m3/s
    if( BCPS(LL).NQSMUL == 1 ) BCPS(LL).RQSMUL = DYP(LTMP)            ! *** Unit discharge along DY
    if( BCPS(LL).NQSMUL == 2 ) BCPS(LL).RQSMUL = DXP(LTMP)            ! *** Unit discharge along DX
    if( BCPS(LL).NQSMUL == 3 ) BCPS(LL).RQSMUL = DXP(LTMP)+DYP(LTMP)  ! *** Unit discharge along DX & DY
    if( BCPS(LL).NQSMUL == 4 ) BCPS(LL).RQSMUL = DXP(LTMP)*DYP(LTMP)  ! *** Seepage velocity over the cell area
  enddo

  do NCTL = 1,NQCTL
    RQDW = 1.
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    LTMP = LIJ(IU,JU)
    HYD_STR(NCTL).RQCMUL = 1.                                                   ! *** DEFAULT
    if( HYD_STR(NCTL).NQCMUL == 1 ) HYD_STR(NCTL).RQCMUL = DYP(LTMP)
    if( HYD_STR(NCTL).NQCMUL == 2 ) HYD_STR(NCTL).RQCMUL = DXP(LTMP)
    if( HYD_STR(NCTL).NQCMUL == 3 ) HYD_STR(NCTL).RQCMUL = DXP(LTMP)+DYP(LTMP)
    if( HYD_STR(NCTL).NQCMUL == 4 ) HYD_STR(NCTL).RQCMUL = DXP(LTMP)*DYP(LTMP)
  enddo

  ! *********************************************************************
  ! *** SET THE VELOCITY AND FLUX BOUNDARY CONDITIONS MULTIPLIERS

  ! *** DEFAULT CONDITION
  do L = 2,LA
    RSSBCE(L) = 1.0
    RSSBCW(L) = 1.0
    RSSBCN(L) = 1.0
    RSSBCS(L) = 1.0
    SUBEW(L) = SUB(L) + SUB(LEC(L))
    SVBNS(L) = SVB(L) + SVB(LNC(L))
  enddo

  ! *** STANDARD BORDER CELLS
  do L = 2,LA
    LE = LEC(L)
    LN = LNC(L)
    if( SUBEW(L) < 1.5 )then
      if( SUB(L) < 0.5 .and. SUB(LE) > 0.5 )then
        RSSBCW(L) = 0.0
        RSSBCE(L) = 1.0
      elseif( SUB(L) > 0.5 .and. SUB(LE) < 0.5 )then
        RSSBCW(L) = 1.0
        RSSBCE(L) = 0.0
      else
        RSSBCW(L) = 0.0
        RSSBCE(L) = 0.0
      endif
    endif
    if( SVBNS(L) < 1.5 )then
      if( SVB(L) < 0.5 .and. SVB(LN) > 0.5 )then
        RSSBCS(L) = 0.0
        RSSBCN(L) = 1.0
      elseif( SVB(L) > 0.5 .and. SVB(LN) < 0.5 )then
        RSSBCS(L) = 1.0
        RSSBCN(L) = 0.0
      else
        RSSBCN(L) = 0.0
        RSSBCS(L) = 0.0
      endif
    endif
  enddo

  ! *** FLOW BOUNDARY CONDITIONS
  do LL = 1,NQSIJ
    L = BCPS(LL).L
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    if( SUBEW(L) < 1.5 .and. DXP(L) > DYP(L) .and. SVB(L) == 0. .and. SVB(LNC(L)) == 0. )then
      if( SUB(L) < 0.5 .and. SUB(LE) >  0.5 )then
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      endif
      if( SUB(L) > 0.5 .and. SUB(LE) <  0.5 )then
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      endif
    endif
    if( SVBNS(L) < 1.5 .and. DYP(L) > DXP(L) .and. SUB(L) == 0. .and. SUB(LEC(L)) == 0. )then
      if( SVB(L) < 0.5 .and. SVB(LN) >  0.5 )then
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      endif
      if( SVB(L) >  0.5 .and. SVB(LN) <  0.5 )then
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      endif
    endif
  enddo

  ! *** HYDRAULIC STRUCTURE: UPSTREAM
  do NCTL = 1,NQCTL
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    L = LIJ(IU,JU)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    if( SUBEW(L) < 1.5 .and. DXP(L) > DYP(L) .and. SVB(L) == 0. .and. SVB(LNC(L)) == 0. )then
      if( SUB(L) < 0.5 .and. SUB(LE) > 0.5 )then
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      endif
      if( SUB(L) >  0.5 .and. SUB(LE) < 0.5 )then
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      endif
    endif
    if( SVBNS(L) < 1.5 .and. DYP(L) > DXP(L) .and. SUB(L) == 0. .and. SUB(LEC(L)) == 0. )then
      if( SVB(L) < 0.5 .and. SVB(LN) > 0.5 )then
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      endif
      if( SVB(L) > 0.5 .and. SVB(LN) < 0.5 )then
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      endif
    endif
  enddo

  ! *** HYDRAULIC STRUCTURE: DOWNSTREAM
  do NCTL = 1,NQCTL
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD
    if( ID /=  0 .and. JD /=  0 )then
      L = LIJ(ID,JD)
      LE = LEC(L)
      LN = LNC(L)
      ! *** ACCOUNT FOR DEAD END CHANNELS
      if( SUBEW(L) < 1.5 .and. DXP(L) > DYP(L) .and. SVB(L) == 0. .and. SVB(LNC(L)) == 0. )then
        if( SUB(L) < 0.5 .and. SUB(LE) > 0.5 )then
          RSSBCW(L) = 0.0
          RSSBCE(L) = 2.0
        endif
        if( SUB(L) > 0.5 .and. SUB(LE) < 0.5 )then
          RSSBCW(L) = 2.0
          RSSBCE(L) = 0.0
        endif
      endif
      if( SVBNS(L) < 1.5 .and. DYP(L) > DXP(L) .and. SUB(L) == 0. .and. SUB(LEC(L)) == 0. )then
        if( SVB(L) < 0.5 .and. SVB(LN) > 0.5 )then
          RSSBCS(L) = 0.0
          RSSBCN(L) = 2.0
        endif
        if( SVB(L) > 0.5 .and. SVB(LN) < 0.5 )then
          RSSBCS(L) = 2.0
          RSSBCN(L) = 0.0
        endif
      endif
    endif
  enddo

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: UPSTREAM
  do NWR = 1,NQWR
    IU = WITH_RET(NWR).IQWRU
    JU = WITH_RET(NWR).JQWRU
    L = LIJ(IU,JU)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    if( SUBEW(L) < 1.5 .and. DXP(L) > DYP(L) .and. SVB(L) == 0. .and. SVB(LNC(L)) == 0. )then
      if( SUB(L) < 0.5 .and. SUB(LE) > 0.5 )then
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      endif
      if( SUB(L) > 0.5 .and. SUB(LE) < 0.5 )then
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      endif
    endif
    if( SVBNS(L) < 1.5 .and. DYP(L) > DXP(L) .and. SUB(L) == 0. .and. SUB(LEC(L)) == 0. )then
      if( SVB(L) < 0.5 .and. SVB(LN) > 0.5 )then
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      endif
      if( SVB(L) > 0.5 .and. SVB(LN) < 0.5 )then
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      endif
    endif
  enddo

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: DOWNSTREAM
  do NWR = 1,NQWR
    ID = WITH_RET(NWR).IQWRD
    JD = WITH_RET(NWR).JQWRD
    L = LIJ(ID,JD)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    if( SUBEW(L) < 1.5 .and. DXP(L) > DYP(L) .and. SVB(L) == 0. .and. SVB(LNC(L)) == 0. )then
      if( SUB(L) <  0.5 .and. SUB(LE) >  0.5 )then
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      endif
      if( SUB(L) >  0.5 .and. SUB(LE) <  0.5 )then
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      endif
    endif
    if( SVBNS(L) <  1.5 .and. DYP(L) > DXP(L) .and. SUB(L) == 0. .and. SUB(LEC(L)) == 0. )then
      if( SVB(L) <  0.5 .and. SVB(LN) >  0.5 )then
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      endif
      if( SVB(L) >  0.5 .and. SVB(LN) <  0.5 )then
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      endif
    endif
  enddo

  ! *** OPEN BOUNDARIES
  do LL = 1,NPBS
    I = IPBS(LL)
    J = JPBS(LL)
    L = LIJ(I,J)
    RSSBCS(L) = 0.0
    RSSBCN(L) = 2.0
  enddo
  do LL = 1,NPBW
    I = IPBW(LL)
    J = JPBW(LL)
    L = LIJ(I,J)
    RSSBCW(L) = 0.0
    RSSBCE(L) = 2.0
  enddo
  do LL = 1,NPBE
    I = IPBE(LL)
    J = JPBE(LL)
    L = LIJ(I,J)
    RSSBCW(L) = 2.0
    RSSBCE(L) = 0.0
  enddo
  do LL = 1,NPBN
    I = IPBN(LL)
    J = JPBN(LL)
    L = LIJ(I,J)
    RSSBCS(L) = 2.0
    RSSBCN(L) = 0.0
  enddo

  ! *********************************************************************
  ! *** SET BOUNDARY MOMENTUM SWITCHES FOR FLOW & HEAD CONTROL

  ! *** GLOBAL BOUNDARY CELL LIST
  NBCS = 0

  ! *** FLOW BC'S
  do LL = 1,NQSIJ
    I = BCPS(LL).I
    J = BCPS(LL).J
    L = LIJ(I,J)
    
    ! *** Reset SUBD and SVBD at flow BC's
    SUBD(L) = 1.
    SVBD(L) = 1.
    
    NBCS = NBCS + 1
    LBCS(NBCS) = L
    
    ! *** Set up momentum at the closed face of a cell for boundary flows
    if( ABS(BCPS(LL).NQSMF) > 0  .and. BCPS(LL).NQSMF /= 5 )then
      if( BCPS(LL).QWIDTH <=  0.0 )then
        if( ABS(BCPS(LL).NQSMF) == 1 .or. ABS(BCPS(LL).NQSMF) == 3 ) BCPS(LL).QWIDTH = DYP(L)
        if( ABS(BCPS(LL).NQSMF) == 2 .or. ABS(BCPS(LL).NQSMF) == 4 ) BCPS(LL).QWIDTH = DXP(L)
      endif
    endif
    
    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    if( SUB(L) < 0.5 )then
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    endif
    if( SVB(L) < 0.5 )then
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    endif

    if( BCPS(LL).NQSMF == 5 )then
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
      SAAX(LEC(L)) = 0.        ! *** EAST/WEST MOMENTUM
      SAAY(LNC(L)) = 0.        ! *** NORTH/SOUTH MOMENTUM
      LBERC(NBCS) = LEC(L)
      LBNRC(NBCS) = LNC(L)
    endif

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    if( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS) = L
    if( L < LA-1 )then
      if( SUB(L) < 0.5 .and. (SUB(LEC(L)) > 0.5 .and. SUB(LEC(LEC(L))) > 0.5) )then
        LBERC(NBCS) = LEC(L)
        SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
        SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
      endif
    endif
    if( L > 2 )then
      if( (SUB(L) > 0.5 .and. SUB(LEC(L)) < 0.5) .and. (SUB(LWC(L)) > 0.5 .and. SUB(LWC(LWC(L))) > 0.5) )then
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      endif
    endif

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    if( SVB(L) < 0.5 .and. (SVB(LNC(L)) > 0.5 .and. SVB(LNC(LNC(L))) > 0.5) )then
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    endif
    if( (SVB(L) > 0.5 .and. SVB(LNC(L)) < 0.5) .and. (SVB(LSC(L)) > 0.5 .and. SVB(LSC(LSC(L))) > 0.5) )then
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    endif

  enddo

  ! *** HEAD CONTROL: UPSTREAM
  do NCTL = 1,NQCTL
    RQDW = 1.
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    L = LIJ(IU,JU)
    NBCS = NBCS + 1
    LBCS(NBCS) = L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    if( SUB(L) < 0.5 )then
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    endif
    if( SVB(L) < 0.5 )then
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    endif

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    if( BC_EDGEFACTOR <= 0.0 ) CYCLE

    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS) = L
    if( SUB(L) < 0.5 .and. (SUB(LEC(L)) > 0.5 .and. SUB(LEC(LEC(L))) > 0.5) )then
      LBERC(NBCS) = LEC(L)
      SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
    endif
    if( L > 2 )then
      if( (SUB(L) > 0.5 .and. SUB(LEC(L)) < 0.5) .and. (SUB(LWC(L)) > 0.5 .and. SUB(LWC(LWC(L))) > 0.5) )then
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      endif
    endif

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    if( SVB(L) < 0.5 .and. (SVB(LNC(L)) > 0.5 .and. SVB(LNC(LNC(L))) > 0.5) )then
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    endif
    if( (SVB(L) > 0.5 .and. SVB(LNC(L)) < 0.5) .and. (SVB(LSC(L)) > 0.5 .and. SVB(LSC(LSC(L))) > 0.5) )then
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    endif

  enddo

  ! *** HEAD CONTROL: DOWNSTREAM
  do NCTL = 1,NQCTL
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD
    if( ID /=  0 .and. JD /=  0 )then
      L = LIJ(ID,JD)
      NBCS = NBCS + 1
      LBCS(NBCS) = L

      ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
      if( SUB(L) < 0.5 )then
        SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
      endif
      if( SVB(L) < 0.5 )then
        SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
      endif

      ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
      if( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
      ! *** EAST/WEST MOMENTUM
      LBERC(NBCS) = L
      if( L < LA-2 )then
        if( SUB(L) < 0.5 .and. (SUB(LEC(L)) > 0.5 .and. SUB(LEC(LEC(L))) > 0.5) )then
          LBERC(NBCS) = LEC(L)
          SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
          SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
        endif
      endif
      if( L > 2 )then
        if( (SUB(L  ) > 0.5 .and. SUB(LEC(L)) < 0.5) .and.  &
          (SUB(LWC(L)) > 0.5 .and. SUB(LWC(LWC(L))) > 0.5) )then
          SAAX(L) = BC_EDGEFACTOR
          SAAY(L) = BC_EDGEFACTOR
        endif
      endif

      ! *** NORTH/SOUTH MOMENTUM
      LBNRC(NBCS) = L
      if( SVB(L) < 0.5 .and. (SVB(LNC(L)) > 0.5 .and.  &
        SVB(LNC(LNC(L))) > 0.5) )then
        LBNRC(NBCS) = LNC(L)
        SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
        SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
      endif
      if( (SVB(L     ) > 0.5 .and. SVB(LNC(L)) < 0.5) .and.  &
        (SVB(LSC(L)) > 0.5 .and. SVB(LSC(LSC(L))) > 0.5) )then
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      endif
    endif
  enddo

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: UPSTREAM
  do NWR = 1,NQWR
    IU = WITH_RET(NWR).IQWRU
    JU = WITH_RET(NWR).JQWRU
    L = LIJ(IU,JU)
    NBCS = NBCS + 1
    LBCS(NBCS) = L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    if( SUB(L) < 0.5 )then
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    endif
    if( SVB(L) < 0.5 )then
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    endif

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    if( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS) = L
    if( SUB(L) < 0.5 .and. (SUB(LEC(L)) > 0.5 .and. SUB(LEC(LEC(L))) > 0.5) )then
      LBERC(NBCS) = LEC(L)
      SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
    endif
    if( L > 2 )then
      if( (SUB(L  ) > 0.5 .and. SUB(LEC(L)) < 0.5) .and.  &
        (SUB(LWC(L)) > 0.5 .and. SUB(LWC(LWC(L))) > 0.5) )then
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      endif
    endif

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    if( SVB(L) < 0.5 .and. (SVB(LNC(L)) > 0.5 .and.  &
      SVB(LNC(LNC(L))) > 0.5) )then
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    endif
    if( (SVB(L     ) > 0.5 .and. SVB(LNC(L)) < 0.5) .and. (SVB(LSC(L)) > 0.5 .and. SVB(LSC(LSC(L))) > 0.5) )then
      !LBNRC(NBCS) = LSC(L)
      !SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      !SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    endif

  enddo

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: DOWNSTREAM
  do NWR = 1,NQWR
    ID = WITH_RET(NWR).IQWRD
    JD = WITH_RET(NWR).JQWRD
    L = LIJ(ID,JD)
    NBCS = NBCS + 1
    LBCS(NBCS) = L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    if( SUB(L) < 0.5 )then
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    endif
    if( SVB(L) < 0.5 )then
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    endif

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    if( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
    if( SUB(L) < 0.5 .and. (SUB(LEC(L)) > 0.5 .and. SUB(LEC(LEC(L))) > 0.5) )then
      LBERC(NBCS) = LEC(L)
      SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
    endif
    if( L > 2 )then
      if( (SUB(L  ) > 0.5 .and. SUB(LEC(L)) < 0.5) .and.  &
        (SUB(LWC(L)) > 0.5 .and. SUB(LWC(LWC(L))) > 0.5) )then
        !LBERC(NBCS) = LWC(L)
        !SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
        !SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      endif
    endif

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    if( SVB(L) < 0.5 .and. (SVB(LNC(L)) > 0.5 .and. SVB(LNC(LNC(L))) > 0.5) )then
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    endif
    if( (SVB(L     ) > 0.5 .and. SVB(LNC(L)) < 0.5) .and.  &
      (SVB(LSC(L)) > 0.5 .and. SVB(LSC(LSC(L))) > 0.5) )then
      !LBNRC(NBCS) = LSC(L)
      !SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      !SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    endif
  enddo

  ! *** OPEN BOUNDARIES
  NBCSOP = 0
  NBCSOP2 = 0
  NBCSOP3 = 0
  do LL = 1,NPBS
    I = IPBS(LL)
    J = JPBS(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP + 1
    LOBCS(NBCSOP) = L
    NBCS = NBCS + 1
    LBCS(NBCS) = LNC(L)
    if( ISPBS(LL) <= 1 .or. ISPBS(LL) >= 4 )then
      NBCSOP2 = NBCSOP2 + 1
      LOBCS2(NBCSOP2) = LNC(L)
    endif
  enddo

  do LL = 1,NPBW
    I = IPBW(LL)
    J = JPBW(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP + 1
    LOBCS(NBCSOP) = L
    NBCS = NBCS + 1
    LBCS(NBCS) = LEC(L)
    if( ISPBW(LL) <= 1 .or. ISPBW(LL) >= 4 )then
      NBCSOP2 = NBCSOP2 + 1
      LOBCS2(NBCSOP2) = LEC(L)
    endif
  enddo
  do LL = 1,NPBE
    I = IPBE(LL)
    J = JPBE(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP + 1
    LOBCS(NBCSOP) = L
    NBCS = NBCS + 1
    LBCS(NBCS) = L
    if( ISPBE(LL) <= 1 .or. ISPBE(LL) >= 4 )then
      NBCSOP2 = NBCSOP2 + 1
      LOBCS2(NBCSOP2) = LWC(L)
    endif
  enddo
  do LL = 1,NPBN
    I = IPBN(LL)
    J = JPBN(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP + 1
    LOBCS(NBCSOP) = L
    NBCS = NBCS + 1
    LBCS(NBCS) = L
    if( ISPBN(LL) <= 1 .or. ISPBN(LL) >= 4 )then
      NBCSOP2 = NBCSOP2 + 1
      LOBCS2(NBCSOP2) = LSC(L)
    endif
  enddo

  ! *********************************************************************
  ! *** SET OPEN BOUNDARY FLAGS FOR CONSTITUENTS
  do LL = 1,NCBS
    I = ICBS(LL)
    J = JCBS(LL)
    LCBS(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  enddo
  do LL = 1,NCBW
    I = ICBW(LL)
    J = JCBW(LL)
    LCBW(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  enddo
  do LL = 1,NCBE
    I = ICBE(LL)
    J = JCBE(LL)
    LCBE(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  enddo
  do LL = 1,NCBN
    I = ICBN(LL)
    J = JCBN(LL)
    LCBN(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  enddo

  ! *********************************************************************
  ! *** SET JET-PLUME VOLUMES SOURCES
  do NJP = 1,NQJPIJ
    L = LIJ(JET_PLM(NJP).IQJP,JET_PLM(NJP).JQJP)
    NBCS = NBCS + 1
    LBCS(NBCS) = L
    LBERC(NBCS) = 1
    LBNRC(NBCS) = 1

    if( JET_PLM(NJP).ICALJP == 2 )then
      ! *** WITHDRAWAL CELL
      L = LIJ(JET_PLM(NJP).IUPCJP,JET_PLM(NJP).JUPCJP)
      NBCS = NBCS + 1
      LBCS(NBCS) = L
      LBERC(NBCS) = 1
      LBNRC(NBCS) = 1
    endif
  enddo

  ! *** SET CHANNEL HOST AND GUEST LOCATION MAPPINGS
  if( MDCHH >=  1 )then
    do NMD = 1,MDCHH
      L = LIJ(IMDCHH(NMD),JMDCHH(NMD))
      ! *** HOST
      LMDCHH(NMD) = L
      NBCS = NBCS + 1
      LBCS(NBCS) = L

      ! *** DOWNSTREAM
      if( IMDCHU(NMD) == 1 .and. JMDCHU(NMD) == 1 )then
        LMDCHU(NMD) = 1
      else
        L = LIJ(IMDCHU(NMD),JMDCHU(NMD))
        LMDCHU(NMD) = L
      endif
      if( IMDCHV(NMD) == 1 .and. JMDCHV(NMD) == 1 )then
        LMDCHV(NMD) = 1
      else
        L = LIJ(IMDCHV(NMD),JMDCHV(NMD))
        LMDCHV(NMD) = L
      endif
      NBCS = NBCS + 1
      LBCS(NBCS) = L
    enddo
  endif

  ! *** SET PERMANENT FACE SWITCHES
  do L = 1,LC
    SUBO(L) = SUB(L)
    SVBO(L) = SVB(L)
  enddo

  ! *** Set partial masks
  if( NBLOCKED > 0 ) CALL BLOCKING

  if( process_id == master_id )then
    ! *** DIAGNOSTIC OUTPUT
    if( DEBUG )then
      open(1,FILE = OUTDIR//'SETBC.DIA',STATUS = 'UNKNOWN')
      close(1,STATUS = 'DELETE')
      open(1,FILE = OUTDIR//'SETBC.DIA')
      do L = 2,LA
        write(1,1001)IL(L),JL(L),SUB(L),SUB(LEC(L)),SVB(L),SVB(LNC(L)),SPB(L)
      enddo
      close(1)
    endif
1001 FORMAT(2I5,8E13.4)
  endif ! End of reading of files from master

  ! *** OBTAIN ACTIVE CELLS SURROUNDING EACH CELL
  LADJ = 1
  do L = 2,LA
    ! *** ORDER OF CELLS
    ! ***   1  2  3
    ! ***   4  5  6
    ! ***   7  8  9
    LADJ(5,L) = L
    if( SUBO(L)      + SVBO(LNWC(L)) > 1.5 .or. SVBO(LNC(L)) + SUBO(LNC(L))  > 1.5 ) LADJ(1,L) = LNWC(L)
    if( SVBO(LNC(L)) + SUBO(LNEC(L)) > 1.5 .or. SUBO(LEC(L)) + SVBO(LNEC(L)) > 1.5 ) LADJ(3,L) = LNEC(L)
    if( SUBO(L)      + SVBO(LWC(L))  > 1.5 .or. SVBO(L)      + SUBO(LSC(L))  > 1.5 ) LADJ(7,L) = LSWC(L)
    if( SVBO(L)      + SUBO(LSEC(L)) > 1.5 .or. SUBO(LEC(L)) + SVBO(LEC(L))  > 1.5 ) LADJ(9,L) = LSEC(L)

    if( SVBO(LNC(L)) > 0.5 ) LADJ(2,L) = LNC(L)
    if( SUBO(L) > 0.5 )      LADJ(4,L) = LWC(L)
    if( SUBO(LEC(L)) > 0.5 ) LADJ(6,L) = LEC(L)
    if( SVBO(L) > 0.5 )      LADJ(8,L) = LSC(L)
  enddo
  
  ! *** Add mapping of RSSBCE,RSSBCW, RSSBCN, RSSBCS to global
  call Get_LA_Local_No_Ghost

  ! *** Calls all of the remapping routines to get the correct global spatial values from RSSBC*
  call ReMap_RSSBC

 RETURN
  
END SUBROUTINE SETBCS
    