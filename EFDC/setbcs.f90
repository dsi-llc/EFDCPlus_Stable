! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

  USE GLOBAL
  USE Variables_MPI
  Use Variables_MPI_Write_Out

  IMPLICIT NONE

  ! *** Local variables
  INTEGER :: L, I, J, NPN, LE, LW, LS, LN, LL, NCTL, IU, JU, IQ, L1, ID, JD, ierr
  INTEGER :: LTMP, NWR, NJP, NMD, NDRYTMP, NC
  REAL    :: DDYDDDX,DDXDDDY,RQDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SUBEW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SVBNS

  IF(  .NOT. ALLOCATED(SUBEW) )THEN
    ALLOCATE(SUBEW(LCM))
    ALLOCATE(SVBNS(LCM))
    SUBEW = 0.0
    SVBNS = 0.0
  ENDIF


  ! *** SET LAND-WATER BOUNDARY SWITCHES
  ITRICELL = 0

  DO L=2,LA
    I = IL(L)
    J = JL(L)
    IF( LCT(L) == 1 )THEN
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
    IF( LCT(L) == 2 )THEN
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
    IF( LCT(L) == 3 )THEN
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
    IF( LCT(L) == 4 )THEN
      STCUV(L) = 0.
      ITRICELL = 1
      STCAP(L) = 0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
    IF( LCT(L) == 5 )THEN
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
    IF( LCT(L) == 6 )THEN
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 6 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 7 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 6 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 7 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
    IF( LCT(L) == 7 )THEN
      IF( IJCT(I-1,J) == 1 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 2 ) SUB(L) = 0.
      IF( IJCT(I-1,J) == 3 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 4 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 5 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 6 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 7 ) SUB(L) = 1.
      IF( IJCT(I-1,J) == 9 ) SUB(L) = 0.
      IF( IJCT(I,J-1) == 1 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 2 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 3 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 4 ) SVB(L) = 0.
      IF( IJCT(I,J-1) == 5 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 6 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 7 ) SVB(L) = 1.
      IF( IJCT(I,J-1) == 9 ) SVB(L) = 0.
    ENDIF
  ENDDO
  SUB(1) = 0.
  SVB(1) = 0.
  SUB(LC) = 0.
  SVB(LC) = 0.

  ! *** Modify land-water boundary conditions for connectors in E-W direction
  IF( ISCONNECT >=  2 )THEN
    DO NPN = 1,NPEWBP
      L = LIJ(IWPEW(NPN),JWPEW(NPN))
      SUB(L) = 1.
      SUBO(L) = 1.
    ENDDO
  ENDIF

  ! *** Modify land-water boundary conditions for connectors in N-S direction
  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    DO NPN = 1,NPNSBP
      LS = LIJ(ISPNS(NPN),JSPNS(NPN))
      SVB(LS)  = 1.
      SVBO(LS) = 1.
    ENDDO
  ENDIF

  ! *** New Mapping routines for MPI Domain Decomp
  Call Map_River
  Call Map_OpenBC_Pressure
  Call Map_OpenBC_Conc
  Call Map_Hydraulic_Structures
  Call Map_Withdrawal_Return
  Call Map_Jet_Plume
  
  ! *** SET WATER-WATER (P OR SURFACE ELEVATION) BOUNDARY SWITCHES
  DO LL = 1,NPBW
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
    IF( ISPBW(LL) == 0 .OR. ISPBW(LL) == 1 .OR. ISPBW(LL) == 4 )THEN
      ! *** Zero Tangiential (East of West BC)
      SWB(LE) = 0.
      SVB(LE) = 0.     
      SVB(LNC(LE)) = 0.
    END IF
  ENDDO
  
  DO LL = 1,NPBE
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
    IF( ISPBE(LL) == 0 .OR. ISPBE(LL) == 1 .OR. ISPBE(LL) == 4 )THEN
      ! *** Zero Tangiential (West of East BC)
      SWB(LW) = 0.
      SVB(LW) = 0.
      SVB(LNC(LW)) = 0.
    END IF
  ENDDO
  
  DO LL = 1,NPBS
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
    IF( ISPBS(LL) == 0 .OR. ISPBS(LL) == 1 .OR. ISPBS(LL) == 4 )THEN
      ! *** Zero tangiential (North of South BC)
      SWB(LN) = 0.
      SUB(LN) = 0.
      SUB(LEC(LN)) = 0.
    END IF
  ENDDO
  
  DO LL = 1,NPBN
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
    IF( ISPBN(LL) == 0 .OR. ISPBN(LL) == 1 .OR. ISPBN(LL) == 4 )THEN
      ! *** Zero Tangiential (South of North BC)
      SWB(LS) = 0.
      SUB(LS) = 0.
      SUB(LEC(LS)) = 0.
    END IF
  ENDDO

  ! *********************************************************************
  ! *** SET THE CELL FACES SWITCHES FOR HEAD CONTROL STRUCTURES
  ! *** UPSTREAM CONTROL
  DO NCTL = 1,NQCTL
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    L = LIJ(IU,JU)
    SAVESUB(1,NCTL) = SUB(L)
    SAVESVB(1,NCTL) = SVB(L)

    ! *** SET U FACE
    LW = LWC(L)
    DO NC = 1,NQCTL
      IF( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
      IF( NC /= NCTL .AND. HYD_STR(NC).NQCTYP /=  3 .AND. HYD_STR(NC).NQCTYP /=  4 )THEN
        I = HYD_STR(NC).IQCTLU
        J = HYD_STR(NC).JQCTLU
        L1 = LIJ(I,J)
        IF( L1 == LW )THEN
          SUB(L) = 0.0
          EXIT
        ENDIF
      ENDIF
    ENDDO

    ! *** SET V FACE
    LS = LSC(L)
    DO NC = 1,NQCTL
      IF( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
      IF( NC /=  NCTL .AND. HYD_STR(NC).NQCTYP /=  3 .AND. HYD_STR(NC).NQCTYP /=  4 )THEN
        I = HYD_STR(NC).IQCTLU
        J = HYD_STR(NC).JQCTLU
        L1 = LIJ(I,J)
        IF( L1 == LS )THEN
          SVB(L) = 0.0
          EXIT
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! *** DOWNSTREAM CONTROL
  DO NCTL = 1,NQCTL
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD
    IF( ID /=  0 .AND. JD /=  0 )THEN
      L = LIJ(ID,JD)
      SAVESUB(2,NCTL) = SUB(L)
      SAVESVB(2,NCTL) = SVB(L)

      ! *** SET U FACE
      LW = LWC(L)
      DO NC = 1,NQCTL
        IF( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
        IF( NC /=  NCTL .AND. HYD_STR(NC).NQCTYP /=  3 .AND. HYD_STR(NC).NQCTYP /=  4 )THEN
          I = HYD_STR(NC).IQCTLD
          J = HYD_STR(NC).JQCTLD
          IF( I > 0 .AND. J > 0 )THEN
            L1 = LIJ(I,J)
            IF( L1 == LW )THEN
              SUB(L) = 0.0
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      ! *** SET V FACE
      LS = LSC(L)
      DO NC = 1,NQCTL
        IF( HYD_STR(NC).MMASKS == 1 ) CYCLE                  ! *** User must manually set masks
        IF( NC /=  NCTL .AND. HYD_STR(NC).NQCTYP /=  3 .AND. HYD_STR(NC).NQCTYP /=  4 )THEN
          I = HYD_STR(NC).IQCTLD
          J = HYD_STR(NC).JQCTLD
          IF( I > 0 .AND. J > 0 )THEN
            L1 = LIJ(I,J)
            IF( L1 == LS )THEN
              SVB(L) = 0.0
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  
  ! *** RESET DXU,DYU,DXV,DYV BASED ON BOUNDARY CONDITION SWITCHES
  DO L=2,LA
    IF( SUB(L) > 0.5 )THEN
      DXU(L) = 0.5*(DXP(L)+DXP(LWC(L)))
      DYU(L) = 0.5*(DYP(L)+DYP(LWC(L)))
    ENDIF
    IF( SUB(L) < 0.5 .AND. SUB(LEC(L)) > 0.5 )THEN
      DXU(L) = DXP(L)
      DDYDDDX = 2.*(DYP(LEC(L))-DYP(L))/(DXP(L)+DXP(LEC(L)))
      DYU(L) = DYP(L)-0.5*DXP(L)*DDYDDDX
    ENDIF
    IF( SUB(L) < 0.5 .AND. SUB(LEC(L)) < 0.5 )THEN
      DXU(L) = DXP(L)
      DYU(L) = DYP(L)
    ENDIF
  ENDDO

  DO L=2,LA
    LN = LNC(L)
    LS = LSC(L)
    IF( SVB(L) > 0.5 )THEN
      DXV(L) = 0.5*(DXP(L)+DXP(LS))
      DYV(L) = 0.5*(DYP(L)+DYP(LS))
    ENDIF
    IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
      DDXDDDY = 2.*(DXP(LN)-DXP(L))/(DYP(L)+DYP(LN))
      DXV(L) = DXP(L)-0.5*DYP(L)*DDXDDDY
      DYV(L) = DYP(L)
    ENDIF
    IF( SVB(L) < 0.5 .AND. SVB(LN) < 0.5 )THEN
      DXV(L) = DXP(L)
      DYV(L) = DYP(L)
    ENDIF
  ENDDO

  ! *** SET THIN BARRIERS BY CALLING CELLMASK
  IF( ISMASK == 1 )  CALL CELLMASK
  IF( NBLOCKED > 0 ) CALL BLOCKING

  ! *** SET VOLUMETRIC & CONCENTRATION SOURCE LOCATIONS AND BED STRESS
  ! *** AND CELL CENTER BED STRESS AND VELOCITY MODIFERS
  DO LL = 1,NQSIJ
    I = IQS(LL)
    J = JQS(LL)
    LTMP = LIJ(I,J)
    LQS(LL) = LTMP
    RQSMUL(LL) = 1.                                 ! *** DEFAULT
    IF( NQSMUL(LL) == 1 ) RQSMUL(LL) = DYP(LTMP)
    IF( NQSMUL(LL) == 2 ) RQSMUL(LL) = DXP(LTMP)
    IF( NQSMUL(LL) == 3 ) RQSMUL(LL) = DXP(LTMP)+DYP(LTMP)
    IF( NQSMUL(LL) == 4 ) RQSMUL(LL) = DXP(LTMP)*DYP(LTMP)
  ENDDO

  DO NCTL = 1,NQCTL
    RQDW = 1.
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    LTMP = LIJ(IU,JU)
    HYD_STR(NCTL).RQCMUL = 1.                                                   ! *** DEFAULT
    IF( HYD_STR(NCTL).NQCMUL == 1 ) HYD_STR(NCTL).RQCMUL = DYP(LTMP)
    IF( HYD_STR(NCTL).NQCMUL == 2 ) HYD_STR(NCTL).RQCMUL = DXP(LTMP)
    IF( HYD_STR(NCTL).NQCMUL == 3 ) HYD_STR(NCTL).RQCMUL = DXP(LTMP)+DYP(LTMP)
    IF( HYD_STR(NCTL).NQCMUL == 4 ) HYD_STR(NCTL).RQCMUL = DXP(LTMP)*DYP(LTMP)
  ENDDO

  ! *********************************************************************
  ! *** SET THE VELOCITY AND FLUX BOUNDARY CONDITIONS MULTIPLIERS

  ! *** DEFAULT CONDITION
  DO L=2,LA
    RSSBCE(L) = 1.0
    RSSBCW(L) = 1.0
    RSSBCN(L) = 1.0
    RSSBCS(L) = 1.0
    SUBEW(L) = SUB(L) + SUB(LEC(L))
    SVBNS(L) = SVB(L) + SVB(LNC(L))
  ENDDO

  ! *** STANDARD BORDER CELLS
  DO L=2,LA
    LE = LEC(L)
    LN = LNC(L)
    IF( SUBEW(L) < 1.5 )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
        RSSBCW(L) = 0.0
        RSSBCE(L) = 1.0
      ELSEIF( SUB(L) > 0.5 .AND. SUB(LE) < 0.5 )THEN
        RSSBCW(L) = 1.0
        RSSBCE(L) = 0.0
      ELSE
        RSSBCW(L) = 0.0
        RSSBCE(L) = 0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
        RSSBCS(L) = 0.0
        RSSBCN(L) = 1.0
      ELSEIF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
        RSSBCS(L) = 1.0
        RSSBCN(L) = 0.0
      ELSE
        RSSBCN(L) = 0.0
        RSSBCS(L) = 0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** FLOW BOUNDARY CONDITIONS
  DO LL = 1,NQSIJ
    L = LQS(LL)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) >  0.5 )THEN
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      ENDIF
      IF( SUB(L) > 0.5 .AND. SUB(LE) <  0.5 )THEN
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) >  0.5 )THEN
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      ENDIF
      IF( SVB(L) >  0.5 .AND. SVB(LN) <  0.5 )THEN
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** HYDRAULIC STRUCTURE: UPSTREAM
  DO NCTL = 1,NQCTL
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    L = LIJ(IU,JU)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      ENDIF
      IF( SUB(L) >  0.5 .AND. SUB(LE) < 0.5 )THEN
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      ENDIF
      IF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** HYDRAULIC STRUCTURE: DOWNSTREAM
  DO NCTL = 1,NQCTL
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD
    IF( ID /=  0 .AND. JD /=  0 )THEN
      L = LIJ(ID,JD)
      LE = LEC(L)
      LN = LNC(L)
      ! *** ACCOUNT FOR DEAD END CHANNELS
      IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
        IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
          RSSBCW(L) = 0.0
          RSSBCE(L) = 2.0
        ENDIF
        IF( SUB(L) > 0.5 .AND. SUB(LE) < 0.5 )THEN
          RSSBCW(L) = 2.0
          RSSBCE(L) = 0.0
        ENDIF
      ENDIF
      IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
        IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
          RSSBCS(L) = 0.0
          RSSBCN(L) = 2.0
        ENDIF
        IF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
          RSSBCS(L) = 2.0
          RSSBCN(L) = 0.0
        ENDIF
      ENDIF
    END IF
  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: UPSTREAM
  DO NWR = 1,NQWR
    IU = WITH_RET(NWR).IQWRU
    JU = WITH_RET(NWR).JQWRU
    L = LIJ(IU,JU)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      ENDIF
      IF( SUB(L) > 0.5 .AND. SUB(LE) < 0.5 )THEN
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      ENDIF
      IF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: DOWNSTREAM
  DO NWR = 1,NQWR
    ID = WITH_RET(NWR).IQWRD
    JD = WITH_RET(NWR).JQWRD
    L = LIJ(ID,JD)
    LE = LEC(L)
    LN = LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) <  0.5 .AND. SUB(LE) >  0.5 )THEN
        RSSBCW(L) = 0.0
        RSSBCE(L) = 2.0
      ENDIF
      IF( SUB(L) >  0.5 .AND. SUB(LE) <  0.5 )THEN
        RSSBCW(L) = 2.0
        RSSBCE(L) = 0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) <  1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) <  0.5 .AND. SVB(LN) >  0.5 )THEN
        RSSBCS(L) = 0.0
        RSSBCN(L) = 2.0
      ENDIF
      IF( SVB(L) >  0.5 .AND. SVB(LN) <  0.5 )THEN
        RSSBCS(L) = 2.0
        RSSBCN(L) = 0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** OPEN BOUNDARIES
  DO LL = 1,NPBS
    I = IPBS(LL)
    J = JPBS(LL)
    L = LIJ(I,J)
    RSSBCS(L) = 0.0
    RSSBCN(L) = 2.0
  ENDDO
  DO LL = 1,NPBW
    I = IPBW(LL)
    J = JPBW(LL)
    L = LIJ(I,J)
    RSSBCW(L) = 0.0
    RSSBCE(L) = 2.0
  ENDDO
  DO LL = 1,NPBE
    I = IPBE(LL)
    J = JPBE(LL)
    L = LIJ(I,J)
    RSSBCW(L) = 2.0
    RSSBCE(L) = 0.0
  ENDDO
  DO LL = 1,NPBN
    I = IPBN(LL)
    J = JPBN(LL)
    L = LIJ(I,J)
    RSSBCS(L) = 2.0
    RSSBCN(L) = 0.0
  ENDDO

  ! *********************************************************************
  ! *** SET BOUNDARY MOMENTUM SWITCHES FOR FLOW & HEAD CONTROL

  ! *** GLOBAL BOUNDARY CELL LIST
  NBCS = 0

  ! *** FLOW BC'S
  DO LL = 1,NQSIJ
    I = IQS(LL)
    J = JQS(LL)
    L = LIJ(I,J)
    NBCS = NBCS+1
    LBCS(NBCS) = L
    
    ! *** Set up momentum at the closed face of a cell for boundary flows
    IF( ABS(NQSMF(LL)) > 0  .AND. NQSMF(LL) /= 5 )THEN
      IF( QWIDTH(LL) <=  0.0 )THEN
        IF( ABS(NQSMF(LL)) == 1 .OR. ABS(NQSMF(LL)) == 3 ) QWIDTH(LL) = DYP(L)
        IF( ABS(NQSMF(LL)) == 2 .OR. ABS(NQSMF(LL)) == 4 ) QWIDTH(LL) = DXP(L)
      ENDIF
    ENDIF
    
    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    IF( SUB(L) < 0.5 )THEN
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    ENDIF

    IF( NQSMF(LL) == 5 )THEN
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
      SAAX(LEC(L)) = 0.        ! *** EAST/WEST MOMENTUM
      SAAY(LNC(L)) = 0.        ! *** NORTH/SOUTH MOMENTUM
      LBERC(NBCS) = LEC(L)
      LBNRC(NBCS) = LNC(L)
    ENDIF

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    IF( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS) = L
    IF( L < LA-1 )THEN
      IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
        LBERC(NBCS) = LEC(L)
        SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
        SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
      ENDIF
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND. (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND. SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND. (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    ENDIF

  ENDDO

  ! *** HEAD CONTROL: UPSTREAM
  DO NCTL = 1,NQCTL
    RQDW = 1.
    IU = HYD_STR(NCTL).IQCTLU
    JU = HYD_STR(NCTL).JQCTLU
    L = LIJ(IU,JU)
    NBCS = NBCS + 1
    LBCS(NBCS) = L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    IF( SUB(L) < 0.5 )THEN
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    ENDIF

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    IF( BC_EDGEFACTOR <= 0.0 ) CYCLE

    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS) = L
    IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
      LBERC(NBCS) = LEC(L)
      SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND. (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND. SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND. (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    ENDIF

  ENDDO

  ! *** HEAD CONTROL: DOWNSTREAM
  DO NCTL = 1,NQCTL
    ID = HYD_STR(NCTL).IQCTLD
    JD = HYD_STR(NCTL).JQCTLD
    IF( ID /=  0 .AND. JD /=  0 )THEN
      L = LIJ(ID,JD)
      NBCS = NBCS + 1
      LBCS(NBCS) = L

      ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
      IF( SUB(L) < 0.5 )THEN
        SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
      ENDIF
      IF( SVB(L) < 0.5 )THEN
        SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
      ENDIF

      ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
      IF( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
      ! *** EAST/WEST MOMENTUM
      LBERC(NBCS) = L
      IF( L < LA-2 )THEN
        IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
          LBERC(NBCS) = LEC(L)
          SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
          SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
        ENDIF
      ENDIF
      IF( L > 2 )THEN
        IF( (SUB(L  ) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND.  &
          (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
          SAAX(L) = BC_EDGEFACTOR
          SAAY(L) = BC_EDGEFACTOR
        ENDIF
      ENDIF

      ! *** NORTH/SOUTH MOMENTUM
      LBNRC(NBCS) = L
      IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND.  &
        SVB(LNC(LNC(L))) > 0.5) )THEN
        LBNRC(NBCS) = LNC(L)
        SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
        SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
      ENDIF
      IF( (SVB(L     ) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND.  &
        (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      ENDIF
    ENDIF
  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: UPSTREAM
  DO NWR = 1,NQWR
    IU = WITH_RET(NWR).IQWRU
    JU = WITH_RET(NWR).JQWRU
    L = LIJ(IU,JU)
    NBCS = NBCS + 1
    LBCS(NBCS) = L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    IF( SUB(L) < 0.5 )THEN
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    ENDIF

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    IF( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS) = L
    IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
      LBERC(NBCS) = LEC(L)
      SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L  ) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND.  &
        (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND.  &
      SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L     ) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND. (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      !LBNRC(NBCS) = LSC(L)
      !SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      !SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    ENDIF

  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: DOWNSTREAM
  DO NWR = 1,NQWR
    ID = WITH_RET(NWR).IQWRD
    JD = WITH_RET(NWR).JQWRD
    L = LIJ(ID,JD)
    NBCS = NBCS + 1
    LBCS(NBCS) = L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    IF( SUB(L) < 0.5 )THEN
      SAAX(L) = 0.             ! *** EAST/WEST MOMENTUM
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L) = 0.             ! *** NORTH/SOUTH MOMENTUM
    ENDIF

    ! *** BC_EDGEFACTOR IS A USER DEFINED MULTIPLIER TO ADJUST BC CELL MOMENTUM FOR SIMPLE CHANNEL US OR DS CELLS
    IF( BC_EDGEFACTOR <= 0.0 ) CYCLE
    
    IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
      LBERC(NBCS) = LEC(L)
      SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L  ) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND.  &
        (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        !LBERC(NBCS) = LWC(L)
        !SAAX(LBERC(NBCS)) = BC_EDGEFACTOR
        !SAAY(LBERC(NBCS)) = BC_EDGEFACTOR
        SAAX(L) = BC_EDGEFACTOR
        SAAY(L) = BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS) = L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND. SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS) = LNC(L)
      SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L     ) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND.  &
      (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      !LBNRC(NBCS) = LSC(L)
      !SAAX(LBNRC(NBCS)) = BC_EDGEFACTOR
      !SAAY(LBNRC(NBCS)) = BC_EDGEFACTOR
      SAAX(L) = BC_EDGEFACTOR
      SAAY(L) = BC_EDGEFACTOR
    ENDIF
  ENDDO

  ! *** OPEN BOUNDARIES
  NBCSOP = 0
  NBCSOP2 = 0
  DO LL = 1,NPBS
    I = IPBS(LL)
    J = JPBS(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP + 1
    LOBCS(NBCSOP) = L
    NBCS = NBCS + 1
    LBCS(NBCS) = LNC(L)
    !LBERC(NBCS) = L
    !LBNRC(NBCS) = LNC(L)
    IF( ISPBS(LL) <= 1 .OR. ISPBS(LL) == 4 )THEN
      NBCSOP2 = NBCSOP2+1
      LOBCS2(NBCSOP2) = LNC(L)
    ENDIF
  ENDDO

  DO LL = 1,NPBW
    I = IPBW(LL)
    J = JPBW(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP+1
    LOBCS(NBCSOP) = L
    NBCS = NBCS+1
    LBCS(NBCS) = LEC(L)
    !LBERC(NBCS) = LEC(L)
    !LBNRC(NBCS) = L
    IF( ISPBW(LL) <= 1 .OR. ISPBW(LL) == 4 )THEN
      NBCSOP2 = NBCSOP2+1
      LOBCS2(NBCSOP2) = LEC(L)
    ENDIF
  ENDDO
  DO LL = 1,NPBE
    I = IPBE(LL)
    J = JPBE(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP+1
    LOBCS(NBCSOP) = L
    NBCS = NBCS+1
    LBCS(NBCS) = L
    !LBERC(NBCS) = LWC(L)
    !LBNRC(NBCS) = L
    IF( ISPBE(LL) <= 1 .OR. ISPBE(LL) == 4 )THEN
      NBCSOP2 = NBCSOP2+1
      LOBCS2(NBCSOP2) = LWC(L)
    ENDIF
  ENDDO
  DO LL = 1,NPBN
    I = IPBN(LL)
    J = JPBN(LL)
    L = LIJ(I,J)
    NBCSOP = NBCSOP+1
    LOBCS(NBCSOP) = L
    NBCS = NBCS+1
    LBCS(NBCS) = L
    !LBERC(NBCS) = L
    !LBNRC(NBCS) = LSC(L)
    IF( ISPBN(LL) <= 1 .OR. ISPBN(LL) == 4 )THEN
      NBCSOP2 = NBCSOP2+1
      LOBCS2(NBCSOP2) = LSC(L)
    ENDIF
  ENDDO

  ! *********************************************************************
  ! ***  SET OPEN BOUNDARY FLAGS FOR CONSTITUENTS
  DO LL = 1,NCBS
    I = ICBS(LL)
    J = JCBS(LL)
    LCBS(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  ENDDO
  DO LL = 1,NCBW
    I = ICBW(LL)
    J = JCBW(LL)
    LCBW(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  ENDDO
  DO LL = 1,NCBE
    I = ICBE(LL)
    J = JCBE(LL)
    LCBE(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  ENDDO
  DO LL = 1,NCBN
    I = ICBN(LL)
    J = JCBN(LL)
    LCBN(LL) = LIJ(I,J)
    L = LIJ(I,J)
    SCB(L) = 0.
  ENDDO

  ! *********************************************************************
  ! ***  SET JET-PLUME VOLUMES SOURCES
  DO NJP = 1,NQJPIJ
    L = LIJ(IQJP(NJP),JQJP(NJP))
    NBCS = NBCS+1
    LBCS(NBCS) = L
    LBERC(NBCS) = 1
    LBNRC(NBCS) = 1

    IF( ICALJP(NJP) == 2 )THEN
      ! *** WITHDRAWAL CELL
      L = LIJ(IUPCJP(NJP),JUPCJP(NJP))
      NBCS = NBCS+1
      LBCS(NBCS) = L
      LBERC(NBCS) = 1
      LBNRC(NBCS) = 1
    ENDIF
  ENDDO

  ! *** SET CHANNEL HOST AND GUEST LOCATION MAPPINGS
  IF( MDCHH >=  1 )THEN
    DO NMD = 1,MDCHH
      L = LIJ(IMDCHH(NMD),JMDCHH(NMD))
      ! *** HOST
      LMDCHH(NMD) = L
      NBCS = NBCS+1
      LBCS(NBCS) = L

      ! *** DOWNSTREAM
      IF( IMDCHU(NMD) == 1 .AND. JMDCHU(NMD) == 1 )THEN
        LMDCHU(NMD) = 1
      ELSE
        L = LIJ(IMDCHU(NMD),JMDCHU(NMD))
        LMDCHU(NMD) = L
      ENDIF
      IF( IMDCHV(NMD) == 1 .AND. JMDCHV(NMD) == 1 )THEN
        LMDCHV(NMD) = 1
      ELSE
        L = LIJ(IMDCHV(NMD),JMDCHV(NMD))
        LMDCHV(NMD) = L
      ENDIF
      NBCS = NBCS+1
      LBCS(NBCS) = L
    ENDDO
  ENDIF

  ! *** SET CELL FACE WET DEPTHS
  HUWET(1) = HWET
  HUWET(LC) = HWET
  HVWET(1) = HWET
  HVWET(LC) = HWET
  HUDRY(1) = HDRY
  HUDRY(LC) = HDRY
  HVDRY(1) = HDRY
  HVDRY(LC) = HDRY
  DO L=2,LA
    LS = LSC(L)
    HUDRY(L) = HDRY+0.5*ABS(BELV(L)-BELV(LWC(L)))
    HVDRY(L) = HDRY+0.5*ABS(BELV(L)-BELV(LS))
    HUWET(L) = HWET+0.5*ABS(BELV(L)-BELV(LWC(L)))
    HVWET(L) = HWET+0.5*ABS(BELV(L)-BELV(LS))
  ENDDO

  IF( ISDRY > 0 )THEN
    NDRYTMP = MOD(ISDRY,2)
    IF( NDRYTMP /=  0 )THEN
      DO L=2,LA
        HUWET(L) = HWET
        HVWET(L) = HWET
        HUDRY(L) = HDRY
        HVDRY(L) = HDRY
      ENDDO
    ENDIF
  ENDIF

  ! *** SET PERMANENT FACE SWITCHES
  DO L = 1,LC
    SUBO(L) = SUB(L)
    SVBO(L) = SVB(L)
  ENDDO

  if( process_id == master_id )THEN
    ! *** DIAGNOSTIC OUTPUT
    IF( DEBUG )THEN
      OPEN(1,FILE = OUTDIR//'SETBC.DIA',STATUS = 'UNKNOWN')
      CLOSE(1,STATUS = 'DELETE')
      OPEN(1,FILE = OUTDIR//'SETBC.DIA')
      DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),SUB(L),SUB(LEC(L)),SVB(L),SVB(LNC(L)),SPB(L)
      ENDDO
      CLOSE(1)
    ENDIF
1001 FORMAT(2I5,8E13.4)
  End if ! End of reading of files from master

  ! *** Add mapping of RSSBCE,RSSBCW, RSSBCN, RSSBCS to global
  Call Get_LA_Local_No_Ghost

  ! *** Calls all of the remapping routines to get the correct global spatial values from RSSBC*
  Call ReMap_RSSBC

 RETURN
  
END SUBROUTINE SETBCS

