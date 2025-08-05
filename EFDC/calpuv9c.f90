  ! ----------------------------------------------------------------------
  !   This file is a part of EFDC+
  !   Website:  https://eemodelingsystem.com/
  !   Repository: https://github.com/dsi-llc/EFDC_Plus.git
  ! ----------------------------------------------------------------------
  ! Copyright 2021-2024 DSI, LLC
  ! Distributed under the GNU GPLv2 License.
  ! ----------------------------------------------------------------------
  SUBROUTINE CALPUV9C

  ! *** *******************************************************************C
  !
  ! *** SUBROUTINE CALPUV9C CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
  ! *** AND VHDXE, FOR FREE SURFACE FLOWS FOR THE 3TL SOLUTION.
  ! *** WITH PROVISIONS FOR WETTING AND DRYING OF CELLS

  !----------------------------------------------------------------------C
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
  !    2015-02       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH WITH DRY BYPASS
  !    2012-09       Chung Dang        Added OMP
  !    2011-03       Paul M. Craig     Rewritten to F90

  use GLOBAL
  use EFDCOUT
  use Mod_Map_Write_EE_Binary
  use Allocate_Initialize

  use MPI
  use MPI_All_Reduce
  use Communicate_Ghost_Routines
  use Mod_Map_Write_EE_Binary

  implicit none

  integer :: NMD, ITERHP, NCORDRY, ICORDRY, ND, LF, LL, L, LP, LS, LN, LW, LE, LHOST, LCHNU, LCHNV, NTMP, I, KM
  integer :: IUW, IUE, IVS, IVN, IFACE, LMAX, LMIN, IFLAG, IMAX, IMIN, JMAX, JMIN, K, NEGFLAG, NNEG, NEGCOUNT
  integer, save :: INOTICE, NCORDRYMAX, NCORDRYAVG, ITERMAX, ITERAVG, NITERAVG, NCOUNT
  real    :: DELTD2, RLAMN, RLAMO, TMPX, TMPY, C1, TMPVAL, HOLDTMP, BELVAVG, RVAL
  real    :: SVPW1, HPPMC, HDRY10, HDRY90
  real    :: SUBW, SUBE, SVBS, SVBN, DHPDT, DHPDT2, RDRY, CCMNM, CCMNMI
  real    :: RNPORI, DIVEXMX, DIVEXMN, DIVEX, ETGWTMP, ETGWAVL, EL, EW, ES
  real(RKD), save :: DAYOLD, DAYOLD30

  integer,save,allocatable,dimension(:) :: IACTIVE
  integer,save,allocatable,dimension(:) :: ICORDRYD
  integer,save,allocatable,dimension(:) :: NNATDRY
  integer,save,allocatable,dimension(:,:) :: LNATDRY

  real,save,allocatable,dimension(:) :: QCHANUT
  real,save,allocatable,dimension(:) :: QCHANVT
  real,save,allocatable,dimension(:) :: QSUMTMP
  real,save,allocatable,dimension(:) :: SUB1
  real,save,allocatable,dimension(:) :: SVB1

  real,save,allocatable,dimension(:) :: CCMNMD

  real,save,allocatable,dimension(:) :: FSGZUDXYPI
  real,save,allocatable,dimension(:) :: FSGZVDXYPI

  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, TWAIT, TTCON, TTPUV                 ! MODEL TIMING TEMPORARY VARIABLE

  ! *** New variables for MPI
  real    :: CCMNM_Local
  integer :: ICORDRY_GLobal !< Dry couunter for the entire domain, takes the sum from all processes
  integer :: ierr, IUPDATE, IUPDATE_Global

  ! *** Initialize timing
  TTWAIT = 0.
  TTPUV = DSTIME(0)
  TTCON = TCONG
  TMPITMP = 0.

  if( .not. allocated(IACTIVE) )then
    call AllocateDSI(IACTIVE, NCHANM,   0)
    call AllocateDSI(QCHANUT, NCHANM, 0.0)
    call AllocateDSI(QCHANVT, NCHANM, 0.0)
    call AllocateDSI(FSGZUDXYPI, LCM, 0.0)
    call AllocateDSI(FSGZVDXYPI, LCM, 0.0)
    call AllocateDSI(QSUMTMP,    LCM, 0.0)
    call AllocateDSI(SUB1,       LCM, 0.0)
    call AllocateDSI(SVB1,       LCM, 0.0)

    call AllocateDSI(LNATDRY,  NTHREADS, LCM,  0)
    call AllocateDSI(CCMNMD,   NTHREADS, 0.0)
    call AllocateDSI(ICORDRYD, NTHREADS,   0)
    call AllocateDSI(NNATDRY,  NTHREADS,   0)

    IACTIVE = 0
    QCHANUT = 0.
    QCHANVT = 0.
    QSUMTMP = 0.
    SUB1 = SUB
    SVB1 = SVB
    ISCDRY = 0
    CCMNMD = 0
    ICORDRYD = 0
    RCX = 1.0
    RCY = 1.0
    RCX(1) = 0.
    RCY(1) = 0.
    RCX(LC) = 0.
    RCY(LC) = 0.
    LNATDRY = 0
    NNATDRY = 0
    INOTICE = -999
    NCORDRYMAX = 0
    NCORDRYAVG = 0
    ITERMAX = 0
    ITERAVG = 0
    NITERAVG = 0
    NCOUNT = 0
    DAYOLD = INT(TIMEDAY)
    DAYOLD30 = DAYOLD + 30.

    ! INITIALIZE DIAGONAL
    CC = 1.0

    FSGZUDXYPI = 1.0
    FSGZVDXYPI = 1.0
    do L = 2,LA
      FSGZUDXYPI(L) =  0.5/DXYU(L)
      FSGZVDXYPI(L) =  0.5/DXYV(L)
    enddo
  endif   ! *** End allocation and initialization of various arrays

  DELT = DT2
  DELTD2 = DT
  if( ISTL == 2 )then
    DELT = DT
    DELTD2 = 0.5*DT
  endif
  DELTI = 1./DELT

  ! ***  INITIALIZE SUBGRID SCALE CHANNEL INTERACTIONS
  if( MDCHH >= 1 )then
    RLAMN = QCHERR
    RLAMO = 1.-RLAMN
    do NMD = 1,MDCHH
      QCHANUT(NMD) = QCHANUN(NMD)
      QCHANVT(NMD) = QCHANVN(NMD)
    enddo
    if( ISTL == 3 )then
      do NMD = 1,MDCHH
        QCHANUN(NMD) = QCHANU(NMD)
        QCHANVN(NMD) = QCHANV(NMD)
      enddo
    endif
  endif

  ! *** SET SWITCHES FOR DRYING AND WETTING
  ITERHP  = 0
  NCORDRY = 0
  NEGFLAG = 0
  NEGCOUNT = 0
  ICORDRY_GLobal = 0

  if( ISDRY > 0 )then
    do LL = 1,NBCSOP
      L = LOBCS(LL)
      LOPENBCDRY(L) = .FALSE.
    enddo
  endif

  ! ***************************************************************************
  ! *** CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N)
  if( BSC > 1.E-6 )then
    call CALEBI
  endif

  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NOPTIMAL(1))

  if( ISDRY > 0 )then
    !$OMP DO PRIVATE(ND,LF,LL,L)
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)

      ! *** SET SWITCHES FOR DRYING AND WETTING
      do L = LF,LL
        ISCDRY(L) = 0
        SUB1(L) = SUB(L)
        SVB1(L) = SVB(L)
        OLDMASK(L) = LMASKDRY(L)
      enddo

      ! *** ZERO VOLUME WASTING ARRAY COUNTERS
      if( NDRYSTP > 0 )then
        NNATDRY(ND) = 0
        !IF( ISBAL > 0 )then
        !  do L = LF,LL
        !    QDWASTE(L) = 0.0
        !  enddo
        !ENDIF
      endif
    enddo
    !$OMP END DO
  endif

  ! *** NOTE:  H1P ARE THE DEPTHS AND UHDY1E/VHDX1E ARE THE FLOWS FROM THE N-2 ITERATION
  ! ***        UNTIL THE VARIABLES ARE ADVANCED WHEN ISTL = 3 BELOW

  !$OMP DO PRIVATE(ND,LF,LL,L,LW,LS,LN)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    if( BSC > 1.E-6 )then
      ! *** CALCULATE EXPLICIT EXTERNAL DENSITY GRADIENTS
      if( IGRIDV == 0 )then
        do L = LF,LL
          LW = LWC(L)
          FPGXE(L) = ROLD*FPGXE(L) + RNEW*(-SBX(L)*SUBD(L)*HU(L)*GP*( (BI2W(L)+BI2W(LW))*(HP(L)-HP(LW)) + 2.0*HU(L)*(BI1W(L)-BI1W(LW)) + (BEW(L)+BEW(LW))*(BELV(L)-BELV(LW)) ) )

          LS = LSC(L)
          FPGYE(L) = ROLD*FPGYE(L) + RNEW*(-SBY(L)*SVBD(L)*HV(L)*GP*( (BI2S(L)+BI2S(LS))*(HP(L)-HP(LS)) + 2.0*HV(L)*(BI1S(L)-BI1S(LS)) + (BES(L)+BES(LS))*(BELV(L)-BELV(LS)) ) )
        enddo
      else
        do L = LF,LL
          LW = LWC(L)
          FPGXE(L) = ROLD*FPGXE(L) + RNEW*(-SBX(L)*SUBD(L)*HU(L)*GP*( (BI2W(L)+BI2E(LW))*(HPW(L)-HPE(LW)) + 2.0*HU(L)*(BI1W(L)-BI1E(LW)) + (BEW(L)+BEE(LW))*(BELVW(L)-BELVE(LW)) ) )

          LS = LSC(L)
          FPGYE(L) = ROLD*FPGYE(L) + RNEW*(-SBY(L)*SVBD(L)*HV(L)*GP*( (BI2S(L)+BI2N(LS))*(HPS(L)-HPN(LS)) + 2.0*HV(L)*(BI1S(L)-BI1N(LS)) + (BES(L)+BEN(LS))*(BELVS(L)-BELVN(LS)) ) )
        enddo
      endif
    endif

    ! *** SET THE CURRENT FACE DEPTHS INTO HUTMP AND HVTMP
    if( ISTL == 2 )then
      do L = LF,LL
        HUTMP(L) = 0.5*(HU(L)+H1U(L))
        HVTMP(L) = 0.5*(HV(L)+H1V(L))
      enddo
    else
      do L = LF,LL
        HUTMP(L) = HU(L)
        HVTMP(L) = HV(L)
      enddo
    endif
  enddo   ! *** END OF DOMAIN
  !$OMP END DO

  !$OMP DO PRIVATE(ND,LF,LL,L,LW,LS,TMPX,TMPY,k,LN,LE,IFLAG)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)
    ! *** CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
    do L = LF,LL
      LW = LWC(L)
      LS = LSC(L)
      FUHDYE(L) = UHDY1E(L) - DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P1(L)-P1(LW)) + SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX1(L)-RITB1*TBX1(L)) + FCAXE(L) + FPGXE(L) - SNLT*FXE(L))
      FVHDXE(L) = VHDX1E(L) - DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P1(L)-P1(LS)) + SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY1(L)-RITB1*TBY1(L)) - FCAYE(L) + FPGYE(L) - SNLT*FYE(L))
    enddo

    ! *** SET IMPLICIT BOTTOM AND VEGETATION DRAG AS APPROPRIATE
    IFLAG = 0
    if( ISITB >= 1 )then
      ! *** IMPLICIT BOTTOM DRAG WITH VEGETATION
      do L = LF,LL
        TMPX = 1.0
        TMPY = 1.0
        if( UHE(L) /= 0.0) TMPX = U(L,KSZU(L))*HU(L)/UHE(L)
        if( VHE(L) /= 0.0) TMPY = V(L,KSZV(L))*HV(L)/VHE(L)
        RCX(L) = 1./( 1. + TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L) + U(L,KSZU(L))*U(L,KSZU(L))) + DELT*FXVEGE(L) )
        RCY(L) = 1./( 1. + TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L) + V(L,KSZV(L))*V(L,KSZV(L))) + DELT*FYVEGE(L) )
      enddo
      IFLAG = 1
    elseif( ISVEG == 1 )then
      ! *** IMPLICIT VEGETATION DRAG ONLY.  REDUCE BY THE TOTAL OF ENERGY
      do L = LF,LL
        !                    S    1/S
        RCX(L) = 1./( 1. + DELT*FXVEGE(L) )  ! *** dimensionless
        RCY(L) = 1./( 1. + DELT*FYVEGE(L) )  ! *** dimensionless
      enddo
      IFLAG = 1
    endif

    ! *** Apply implicit momentum adjustment
    if( IFLAG > 0 )then

      do L = LF,LL
        FUHDYE(L) = FUHDYE(L)*RCX(L)
        FVHDXE(L) = FVHDXE(L)*RCY(L)
      enddo
    endif
  enddo   ! *** END OF DOMAIN
  !$OMP END DO

  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LS,LN,LW,RDRY,EL,EW,ES,K)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    ! *** RESET BOUNDARY CONDITIONS SWITCHES
    if( ISDRY > 0 )then
      ! ****************************************************************************
      if( NITER < 100 .or. MOD(NITER,50) == 0 .or. ITERHPM == 0 )then
        ! *** FORCE A COMPLETE UPDATE EVERY 100 INTERATIONS
        do L = LF,LL
          SUB(L) = SUBO(L)
          SVB(L) = SVBO(L)
          SBX(L) = SBXO(L)
          SBY(L) = SBYO(L)
        enddo
      else
        ! *** UPDATED TO REDUCE MPI COMMUNICATIONS
        do L = LF,LL
          LE = LEC(L)
          LN = LNC(L)
          RDRY = SUB1(L) + SUB1(LE) + SVB1(L) + SVB1(LN)

          if( RDRY < 0.5 )then
            if( HP(L) > HDRY )then
              ! *** ISOLATED CELL, CHECK ADJACENT CELLS
              EL = BELV(L)  + HP(L)                         ! *** WSEL OF CURRENT CELL

              if( EL > BELV(LWC(L)) )then                   ! *** WEST FACE
                SUB(LWC(L)) = SUBO(LWC(L))
                SBX(LWC(L)) = SBXO(LWC(L))
              endif
              if( EL > BELV(LEC(L)) )then                   ! *** EAST FACE
                SUB(LEC(L)) = SUBO(LEC(L))
                SBX(LEC(L)) = SBXO(LEC(L))
              endif
              if( EL > BELV(LSC(L)) )then                   ! *** SOUTH FACE
                SVB(L) = SVBO(L)
                SBY(L) = SBYO(L)
              endif
              if( EL > BELV(LNC(L)) )then                   ! *** NORTH FACE
                SVB(LNC(L)) = SVBO(LNC(L))
                SBY(LNC(L)) = SBYO(LNC(L))
              endif
            endif

          elseif( RDRY < 3.5 )then
            if( RDRY == (SUBO(L) + SUBO(LE) + SVBO(L) + SVBO(LN)) ) CYCLE  ! *** ALL FACES ACTIVE

            ! *** AT LEAST ONE FACE IS INACTIVE.  CHECK DEPTHS

            ! *** CHECK EAST/WEST FACES
            if( SUB1(L) < 0.5 .and. SUBO(L) > 0.5 )then
              LW = LWC(L)
              if( LW > 1 )then
                EL = BELV(L)  + HP(L)                       ! *** WSEL OF CURRENT CELL
                EW = BELV(LW) + HP(LW)                      ! *** WSEL OF WEST CELL
                if( HP(L) < HDRY .and. HP(LW) < HDRY )then
                  ! *** DO NOTHING
                elseif( HP(L) >= HDRY .and. HP(LW) >= HDRY )then
                  ! *** BOTH CELLS ARE WET, RESET CURRENT CELL
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)

                elseif( HP(LW) < HDRY .and. EL > EW )then
                  ! *** RESET WEST CELL
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)

                  ! *** OPEN OTHER FACES IN WEST CELL
                  if( EW > BELV(LWC(LW)) )then               ! *** WEST FACE
                    SUB(LW) = SUBO(LW)
                    SBX(LW) = SBXO(LW)
                  endif
                  if( EW > BELV(LSC(LW)) )then               ! *** SOUTH FACE
                    SVB(LW) = SVBO(LW)
                    SBY(LW) = SBYO(LW)
                  endif
                  if( EW > BELV(LNC(LW)) )then               ! *** NORTH FACE
                    SVB(LNC(LW)) = SVBO(LNC(LW))
                    SBY(LNC(LW)) = SBYO(LNC(LW))
                  endif

                elseif( HP(L) < HDRY .and. EW >= EL )then
                  ! *** RESET CURRENT CELL
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)

                  ! *** OPEN OTHER FACES IN CURRENT CELL
                  if( EL > BELV(LEC(L)) )then                ! *** EAST FACE
                    SUB(LEC(L)) = SUBO(LEC(L))
                    SBX(LEC(L)) = SBXO(LEC(L))
                  endif
                  if( EL > BELV(LSC(L)) )then                ! *** SOUTH FACE
                    SVB(L) = SVBO(L)
                    SBY(L) = SBYO(L)
                  endif
                  if( EL > BELV(LNC(L)) )then                ! *** NORTH FACE
                    SVB(LNC(L)) = SVBO(LNC(L))
                    SBY(LNC(L)) = SBYO(LNC(L))
                  endif

                endif
              endif
            endif

            ! *** CHECK NORTH/SOUTH FACES
            if( SVB1(L) < 0.5 .and. SVBO(L) > 0.5 )then
              LS = LSC(L)
              if( LS > 1 )then
                EL = BELV(L)  + HP(L)
                ES = BELV(LS) + HP(LS)
                if( HP(L) < HDRY .and. HP(LS) < HDRY )then
                  ! *** DO NOTHING
                elseif( HP(L) >= HDRY .and. HP(LS) >= HDRY )then
                  ! *** BOTH CELLS ARE WET, RESET CURRENT CELL
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)

                elseif( HP(LS) < HDRY .and. EL > ES )then
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)

                  ! *** OPEN OTHER FACES IN SOUTH CELL
                  if( ES > BELV(LWC(LS)) )then               ! *** WEST FACE
                    SUB(LS) = SUBO(LS)
                    SBX(LS) = SBXO(LS)
                  endif
                  if( ES > BELV(LEC(LS)) )then               ! *** EAST FACE
                    SUB(LEC(LS)) = SUBO(LEC(LS))
                    SBX(LEC(LS)) = SBXO(LEC(LS))
                  endif
                  if( ES > BELV(LSC(LS)) )then               ! *** SOUTH FACE
                    SVB(LSC(LS)) = SVBO(LSC(LS))
                    SBY(LSC(LS)) = SBYO(LSC(LS))
                  endif

                elseif( HP(L) < HDRY .and. ES >= EL )then
                  ! *** RESET CURRENT CELL
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)

                  ! *** OPEN OTHER FACES IN CURRENT CELL
                  if( EL > BELV(LEC(L)) )then                ! *** EAST FACE
                    SUB(LEC(L)) = SUBO(LEC(L))
                    SBX(LEC(L)) = SBXO(LEC(L))
                  endif
                  if( EL > BELV(LWC(L)) )then                ! *** WEST FACE
                    SUB(L) = SUBO(L)
                    SBX(L) = SBXO(L)
                  endif
                  if( EL > BELV(LNC(L)) )then                ! *** NORTH FACE
                    SVB(LNC(L)) = SVBO(LNC(L))
                    SBY(LNC(L)) = SBYO(LNC(L))
                  endif

                endif
              endif
            endif
          elseif( QSUME(L) /= 0.0 )then
            SUB(L) = SUBO(L)
            SVB(L) = SVBO(L)
            SBX(L) = SBXO(L)
            SBY(L) = SBYO(L)
          endif
        enddo

        ! *** ENSURE OPEN BC'S ARE "ON"
        if( ND == 1 )then
          do K = 1,NBCSOP
            L = LOBCS(K)
            SUB(L) = SUBO(L)
            SVB(L) = SVBO(L)
            SBX(L) = SBXO(L)
            SBY(L) = SBYO(L)
          enddo
          do K = 1,NBCSOP2
            L = LOBCS2(K)
            SUB(L) = SUBO(L)
            SVB(L) = SVBO(L)
            SBX(L) = SBXO(L)
            SBY(L) = SBYO(L)
          enddo
        endif
      endif
    else
      ! ****************************************************************************
    endif
  enddo   ! *** END OF DOMAIN
  !$OMP END DO

  ! *** SET OLD TIME LEVEL TERMS IN CONTINUITY EQUATION FOR NON BOUNDARY POINTS
  ! *** DXYP = DXP*DYP (M^2),  G-9.82 (M/S^2), P (M2/S2), FP1 (M4/S3)
  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LN,C1)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    C1 = 0.5*G
    do L = LF,LL
      LE = LEC(L)
      LN = LNC(L)
      FP1(L) = DELTI*DXYP(L)*P1(L) - C1*( UHDY1E(LE)-UHDY1E(L) + VHDX1E(LN)-VHDX1E(L) )
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! *** ADVANCE EXTERNAL VARIABLES FOR THREE TIME LEVEL STEP
  if( ISTL == 3 )then
    !$OMP DO PRIVATE(ND,LF,LL,L,K)
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)

      ! *** ADVANCE EXTERNAL VARIABLES
      do L = LF,LL
        UHDY2E(L) = UHDY1E(L)
        VHDX2E(L) = VHDX1E(L)
        UHDY1E(L) = UHDYE(L)
        VHDX1E(L) = VHDXE(L)
        U1V(L) = UV(L)
        V1U(L) = VU(L)
        P1(L)  = P(L)
        H1U(L) = HU(L)
        H1V(L) = HV(L)
        H1UI(L) = HUI(L)
        H1VI(L) = HVI(L)
        H2P(L)  = H1P(L)
        H1P(L)  = HP(L)
      enddo

      do K = 1,KC
        do L = LF,LL
          UHDY2EK(L,K) = UHDY1EK(L,K)
          VHDX2EK(L,K) = VHDX1EK(L,K)
          UHDY1EK(L,K) = UHDYEK(L,K)
          VHDX1EK(L,K) = VHDXEK(L,K)
        enddo
      enddo

      do K = 1,KC
        do L = LF,LL
          H2PK(L,K) = H1PK(L,K)
          H1PK(L,K) = HPK(L,K)
        enddo
      enddo

      if( ISGWIE >= 1 )then
        do L = LF,LL
          AGWELV2(L) = AGWELV1(L)
          AGWELV1(L) = AGWELV(L)
        enddo
      endif
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO
  endif
  !$OMP END PARALLEL

  ! *** LAYER FACE BLOCKING
  if( BSC > 0.0 .and. NBLOCKED > 0 .and. N > 1 )then
    ND = 1
    C1 = 0.5*G
    do LP = 1,NBLOCKED
      L = LBLOCKED(LP)
      if( KSZ(L) == KC ) CYCLE

      ! *** LAYER U BLOCKING
      if( SUB(L) > 0.0 .and. (BLDRAFTU(LP)+BLSILLU(LP)) > 0.0 )then
        LW = LWC(L)
        FPGXE(L) = -SBX(L)*SUBD(L)*HU(L)*GP*( (BI2W(L)+BI2E(LW))*(HPW(L)-HPE(LW)) + 2.0*HU(L)*(BI1W(L)-BI1E(LW)) + (BEW(L)+BEE(LW))*(BELVW(L)-BELVE(LW)) )          ! *** m3/s

        FUHDYE(L) = UHDYE(L) - DELTD2*SUB(L)*HRUO(L)*HU(L)*(P(L)-P(LW)) + SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L)) + FCAXE(L) + FPGXE(L) - SNLT*FXE(L))   ! *** m3/s
      endif

      ! *** LAYER V BLOCKING
      if( SVB(L) > 0.0 .and. (BLDRAFTV(LP)+BLSILLV(LP)) > 0.0 )then
        LS = LSC(L)
        FPGYE(L) = -SBY(L)*SVBD(L)*HV(L)*GP*( (BI2S(L)+BI2N(LS))*(HPS(L)-HPN(LS)) + 2.0*HV(L)*(BI1S(L)-BI1N(LS)) + (BES(L)+BEN(LS))*(BELVS(L)-BELVN(LS)) )

        FVHDXE(L) = VHDXE(L) - DELTD2*SVB(L)*HRVO(L)*HV(L)*(P(L)-P(LS)) + SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L)) - FCAYE(L) + FPGYE(L) - SNLT*FYE(L))
      endif

      LE = LEC(L)
      LN = LNC(L)

      FP1(L) = DELTI*DXYP(L)*P1(L) - C1*( UHDY1E(LE)-UHDY1E(L) + VHDX1E(LN)-VHDX1E(L) )
    enddo
  endif

  ! *************************************************************************
  ! ***  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
  ! ***  HOST-GUEST CHANNEL INTERACTION FOR NON BOUNDARY POINTS

  !----------------------------------------------------------------------C
  ! ***  REENTER AT 1000 FOR WETTING-DRYING CORRECTION AND CHANNEL
  ! ***  INTERACTION
1000 continue
  !----------------------------------------------------------------------C

  ! *** SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM
  CCMNMD = 1.E+18

  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NOPTIMAL(1))
  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LS,LN,C1)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    C1 = 0.5*G
    do L = LF,LL
      LE = LEC(L)
      LN = LNC(L)
      ! ***  USE THE SUB & SVB SWITCHES FOR MOMENTUM FOR THE CURRENT WET/DRY ITERATION
      ! ***   m4/s3   m/s2            m3/s                 m3/s              m3/s               m3/s            m3/s
      FP(L) = FP1(L) - C1*( SUB(LE)*FUHDYE(LE) - SUB(L)*FUHDYE(L) + SVB(LN)*FVHDXE(LN) - SVB(L)*FVHDXE(L) - 2.0*QSUME(L) )    ! *** m4/s3
    enddo

    if( ISGWIE >= 1 )then
      do L = LF,LL
        !       m4/s3   m/s2             m3/s
        FP(L) = FP(L)-G*SPB(L)*(EVAPSW(L)-QGW(L))    ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT)
      enddo
    endif

    C1 = -0.5*DELTD2*G
    do L = LF,LL
      LE = LEC(L)
      LN = LNC(L)
      CS(L) = C1*SVB(L) *HRVO(L) *RCY(L) *HVTMP(L)
      CW(L) = C1*SUB(L) *HRUO(L) *RCX(L) *HUTMP(L)
      CE(L) = C1*SUB(LE)*HRUO(LE)*RCX(LE)*HUTMP(LE)
      CN(L) = C1*SVB(LN)*HRVO(LN)*RCY(LN)*HVTMP(LN)
    enddo

    ! *** SET THE CENTROID
    do L = LF,LL
      CC(L) = DELTI*DXYP(L) - CS(L) - CW(L) - CE(L) - CN(L)
    enddo

  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP SINGLE

  ! *** APPLY THE OPEN BOUNDARY CONDITIONS
  if( NBCSOP > 0 ) CALL SETOPENBC(DELTD2, HUTMP, HVTMP, NCORDRY)

  ! *** INSERT IMPLICT SUB-GRID SCALE CHANNEL INTERACTIONS
  if( MDCHH >= 1 )CALL SUBCHAN(QCHANUT,QCHANVT,IACTIVE,DELT)
  !$OMP END SINGLE

  ! *** SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM
  !$OMP DO PRIVATE(ND,LF,LL,L)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)
    do L = LF,LL
      CCMNMD(ND) = min(CCMNMD(ND),CC(L))
      FPTMP(L) = FP(L)
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  !$OMP SINGLE
  CCMNM  = 1.E+18
  do ND = 1,NOPTIMAL(1)
    CCMNM = min(CCMNM,CCMNMD(ND))
  enddo

  CCMNMI = 1./CCMNM

  !$OMP END SINGLE

  ! *** SCALE BY MINIMUM DIAGONAL  (IRVEC == 9 IS THE ONLY OPTION NOW)
  ! *** BEGIN DOMAIN LOOP
  !$OMP DO PRIVATE(ND,LF,LL,L)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    do L = LF,LL
      CCS(L) = CS(L)*CCMNMI
      CCW(L) = CW(L)*CCMNMI
      CCE(L) = CE(L)*CCMNMI
      CCN(L) = CN(L)*CCMNMI
      CCC(L) = CC(L)*CCMNMI
      FPTMP(L) = FPTMP(L)*CCMNMI
      CCCI(L) = 1./CCC(L)
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL

  if( MDCHH >= 1 )then
    do NMD = 1,MDCHH
      CCCCHH(NMD) = CCCCHH(NMD)*CCMNMI
    enddo
  endif

  ! *********************************************************************************************
  ! *** CALL THE PRECONDITIONED CONJUGATE GRADIENT SOLVER
  if( num_processors > 1 )then
    if( MDCHH == 0 ) CALL Congrad_MPI    ! *** MSCHH> = 1 not parallelized yet @todo
    TTDS = DSTIME(0)        ! delme - 2TL does not have this communicaiton...
    call MPI_barrier(DSIcomm, ierr)
    TTWAIT = TTWAIT + (DSTIME(0)- TTDS)

    TTDS = DSTIME(0)
    call communicate_ghost_cells(P, 'P')
    TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  else
    if( MDCHH == 0 ) CALL CONGRAD
    if( MDCHH >= 1 ) CALL CONGRADC
  endif
  ! *********************************************************************************************

  ITERMAX = max(ITERMAX,ITER)
  ITERAVG = ITERAVG + ITER
  NITERAVG = NITERAVG + 1

  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NOPTIMAL(1))
  !$OMP DO PRIVATE(ND,LF,LL,L,LS,LW)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    ! *** CELL FACE DISCHARGE (M3/S)
    do L = LF,LL
      LS = LSC(L)
      LW = LWC(L)
      UHDYE(L) = SUB(L)*( FUHDYE(L) - DELTD2*HRUO(L)*RCX(L)*HUTMP(L)*(P(L)-P(LW)) )
      VHDXE(L) = SVB(L)*( FVHDXE(L) - DELTD2*HRVO(L)*RCY(L)*HVTMP(L)*(P(L)-P(LS)) )
    enddo

    ! *** UNIT WIDTH DISCHARGE AT CELL FACE (M2/S)
    do L = LF,LL
      UHE(L) = UHDYE(L)*DYIU(L)
      VHE(L) = VHDXE(L)*DXIV(L)
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! *** CALCULATE NEW SUB-GRID SCALE CHANNEL EXCHANGE FLOWS
  !$OMP SINGLE
  if( MDCHH >= 1 )then
    do NMD = 1,MDCHH
      if( IACTIVE(NMD) > 0 )then
        LHOST = LMDCHH(NMD)
        LCHNU = LMDCHU(NMD)
        LCHNV = LMDCHV(NMD)
        if( MDCHTYP(NMD) == 1 )then
          QCHANU(NMD) = CCCCHU(NMD)*QCHANUT(NMD)-RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNU))-RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNU))
          QCHANV(NMD) = 0.
        endif
        if( MDCHTYP(NMD) == 2 )then
          QCHANU(NMD) = 0.
          QCHANV(NMD) = CCCCHU(NMD)*QCHANVT(NMD)-RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNV))-RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNV))
        endif
      else
        QCHANV(NMD) = 0.
        QCHANVN(NMD) = 0.
        QCHANU(NMD) = 0.
        QCHANUN(NMD) = 0.
      endif
    enddo
  endif
  !$OMP END SINGLE

  ! ***  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL TRANSPORTS AT (N+1)
  if( ISTL == 3 )then
    !$OMP DO PRIVATE(ND,LF,LL,L,LE,LN)
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)
      do L = LF,LL
        LE = LEC(L)
        LN = LNC(L)
        HP(L) = H2P(L) + DELT*DXYIP(L)*( QSUME(L) - 0.5*( UHDYE(LE) + UHDY2E(LE) - UHDYE(L) - UHDY2E(L) &
          + VHDXE(LN) + VHDX2E(LN) - VHDXE(L) - VHDX2E(L)) )
      enddo
    enddo
    !$OMP END DO
  else
    !$OMP DO PRIVATE(ND,LF,LL,L,LE,LN)
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)
      do L = LF,LL
        LE = LEC(L)
        LN = LNC(L)
        HP(L) = H1P(L) + DELT*DXYIP(L)*( QSUME(L) - 0.5*( UHDYE(LE) + UHDY1E(LE) - UHDYE(L) - UHDY1E(L) &
          + VHDXE(LN) + VHDX1E(LN) - VHDXE(L) - VHDX1E(L)) )
      enddo
    enddo
    !$OMP END DO
  endif
  !$OMP END PARALLEL

  if( ISGWIE >= 1 )then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LN) NUM_THREADS(NOPTIMAL(1))
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)
      do L = LF,LL
        HP(L) = HP(L) - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT)
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
  endif

  ! *** APPLY OPEN BOUNDARYS
  do LL = 1,NBCSOP
    L = LOBCS(LL)
    HP(L) = GI*P(L) - BELV(L)

    ! *** CHECK FOR CONDITION OF PSERT<BELV
    if( ISDRY > 0 )then

      ! *** HP CHECKING IS FOR RADIATION BC'S **
      if( (LOPENBCDRY(L) .or. (H1P(L) < 0.9*HDRY .and. HP(L) < H1P(L)) .or. HP(L) < 0. ) .and. ISCDRY(L) == 0 )then
        if( HP(L) <= 0. .and. LOPENBCDRY(L) )then
          PRINT '(A,I5,3F10.4,F12.4)',' WARNING!  NEG DEPTH AT OPEN BC: L,HP,H1P,H2P,TIMEDAY',Map2Global(L).LG,HP(L),H1P(L),H2P(L),TIMEDAY
          open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
          write(mpi_efdc_out_unit,'(A,I5,3F10.4,F12.4)')' WARNING!  NEG DEPTH AT OPEN BC: L,HP,H1P,H2P,TIMEDAY',Map2Global(L).LG,HP(L),H1P(L),H2P(L),TIMEDAY
          close(mpi_efdc_out_unit)
          HP(L) = 0.2*HDRY
        else
          HP(L) = min(MAX(H1P(L),0.1*HDRY),0.9*HDRY)
        endif

        ISCDRY(L) = 1
        ICORDRY = -999

        LE = LEC(L)
        LN = LNC(L)
        SUB(L) = 0.
        SVB(L) = 0.
        SUB(LE) = 0.
        SVB(LN) = 0.
        SBX(L) = 0.
        SBY(L) = 0.
        SBX(LE) = 0.
        SBY(LN) = 0.
        UHDYE(L) = 0.
        UHDYE(LE) = 0.
        VHDXE(L) = 0.
        VHDXE(LN) = 0.

        LOPENBCDRY(L) = .TRUE.
        CC(L)  = DELTI*DXYP(L)
        P(L)   = (HP(L)+BELV(L))*G
        FP1(L) = DELTI*DXYP(L)*P(L)
      endif
    endif
  enddo

  ! *** Exchange HP ghost values
  ! ****************************************************************************
  TTDS = DSTIME(0)
  call MPI_barrier(DSIcomm, ierr)
  TTWAIT = TTWAIT + (DSTIME(0)- TTDS)

  TTDS = DSTIME(0)
  call communicate_ghost_cells(HP, 'HP')
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  ! ****************************************************************************

  ! *** ADD CHANNEL INTERACTION EXCHANGES
  if( MDCHH >= 1 )then
    do NMD = 1,MDCHH
      if( IACTIVE(NMD) > 0 )then
        LHOST = LMDCHH(NMD)
        LCHNU = LMDCHU(NMD)
        LCHNV = LMDCHV(NMD)
        if( MDCHTYP(NMD) == 1 )then
          TMPVAL = DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUT(NMD))
          HP(LHOST) = HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNU) = HP(LCHNU)-TMPVAL*DXYIP(LCHNU)
        endif
        if( MDCHTYP(NMD) == 2 )then
          TMPVAL = DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVT(NMD))
          HP(LHOST) = HP(LHOST)+TMPVAL*DXYIP(LHOST)
          HP(LCHNV) = HP(LCHNV)-TMPVAL*DXYIP(LCHNV)
        endif
      endif
    enddo
  endif

  ! *** Check for wet and dry cell changes
  ICORDRY = 0
  ICORDRYD = 0
  !IF( ISDRY > 0 )then
  HDRY90 = 0.9*HDRY
  HDRY10 = 0.1*HDRY

  !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NOPTIMAL(1))

  ! *** CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY
  if( ISDRY > 0 .and. ISDRY < 98 )then
    !$OMP DO PRIVATE(ND,LF,LL,L,LN,LE,K,RVAL,HOLDTMP,HPPMC)
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)

      do L = LF,LL
        if( HP(L) < HDRY )then
          if( ISTL == 3 )then
            HPPMC = H2P(L)
          else
            HPPMC = H1P(L)
          endif
          if( ( HP(L) < HDRY90 .or. HPPMC < HDRY ) .and. HP(L) < (HPPMC+HDRY10) )then
            if( ISCDRY(L) == 0 )then
              ISCDRY(L) = ISCDRY(L)+1
              ICORDRYD(ND) = 1
            endif
            LE = LEC(L)
            LN = LNC(L)
            SUB(L) = 0.
            SVB(L) = 0.
            SUB(LE) = 0.
            SVB(LN) = 0.
            SBX(L) = 0.
            SBY(L) = 0.
            SBX(LE) = 0.
            SBY(LN) = 0.
          endif
        endif
      enddo
    enddo    ! *** END OF DOMAIN LOOP
    !$OMP END DO
  endif

  if( ISDRY == 99 )then
    !$OMP DO PRIVATE(ND,LF,LL,L,LS,LN,LW,LE)  &
    !$OMP    PRIVATE(SUBW,SUBE,SVBS,SVBN,HPPMC,DHPDT,DHPDT2,RDRY,TMPVAL)
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)

      do L = LF,LL
        if( HP(L) <= HDRY )then
          LE = LEC(L)
          LN = LNC(L)

          SUBW = SUB(L)
          SUBE = SUB(LE)
          SVBS = SVB(L)
          SVBN = SVB(LN)
          if( ISTL == 3 )then
            HPPMC = H2P(L)
          else
            HPPMC = H1P(L)
          endif
          DHPDT = HP(L) - HPPMC

          ! *** ALLOW RE-WETTING
          if( DHPDT > 0.0 )then
            ! *** Rising water levels
            LW = LWC(L)
            LS = LSC(L)
            SUB(L)  = 0.0
            SUB(LE) = 0.0
            SVB(L)  = 0.0
            SVB(LN) = 0.0
            SBX(L)  = 0.0
            SBX(LE) = 0.0
            SBY(L)  = 0.0
            SBY(LN) = 0.0

            ! *** RAISING WATER, IS IT FAST ENOUGH TO STAY WET
            DHPDT2 = (HDRY - HPPMC)*0.01
            if( DHPDT > DHPDT2 )then
              if( UHDYE(L) /= 0.0 )then
                if( HP(LW) > HDRY )then
                  SUB(L) = SUBO(L)
                  SBX(L) = SBXO(L)
                endif
              endif
              if( UHDYE(LE) /= 0.0 )then
                if( HP(LE) > HDRY )then
                  SUB(LE) = SUBO(LE)
                  SBX(LE) = SBXO(LE)
                endif
              endif
              if( VHDXE(L) /= 0.0 )then
                if( HP(LS) > HDRY )then
                  SVB(L) = SVBO(L)
                  SBY(L) = SBYO(L)
                endif
              endif
              if( VHDXE(LN) /= 0.0 )then
                if( HP(LN) > HDRY )then
                  SVB(LN) = SVBO(LN)
                  SBY(LN) = SBYO(LN)
                endif
              endif
              RDRY = SUB(L) + SUB(LE) + SVB(L) + SVB(LN)
              if( RDRY < 0.5 )then
                ISCDRY(L) = 1
              else
                ISCDRY(L) = 0
              endif
              TMPVAL = ABS(SUB(L)-SUBW)
              if( TMPVAL > 0.5 )then
                ICORDRYD(ND) = 1
              else
                TMPVAL = ABS(SUB(LE)-SUBE)
                if( TMPVAL > 0.5 )then
                  ICORDRYD(ND) = 1
                else
                  TMPVAL = ABS(SVB(L)-SVBS)
                  if( TMPVAL > 0.5 )then
                    ICORDRYD(ND) = 1
                  else
                    TMPVAL = ABS(SVB(LN)-SVBN)
                    if( TMPVAL > 0.5 )then
                      ICORDRYD(ND) = 1
                    endif
                  endif
                endif
              endif
            else
              ! *** CASE: HP < HDRY BUT RISING, JUST NOT FAST ENOUGH
              if( ISCDRY(L) == 0 )then
                ISCDRY(L) = 1
                ICORDRYD(ND) = 1
              endif
            endif     ! *** END OF REWETTING SECTION, DHPDT > DHPDT2

          elseif( HP(L) < HDRY90 .or. HPPMC < HDRY )then
            ! *** Falling water levels
            ! *** HP < HDRY.  Set switches to dry
            SUB(L)  = 0.0
            SUB(LE) = 0.0
            SVB(L)  = 0.0
            SVB(LN) = 0.0
            SBX(L)  = 0.0
            SBX(LE) = 0.0
            SBY(L)  = 0.0
            SBY(LN) = 0.0
            if( ISCDRY(L) == 0 )then
              ISCDRY(L) = 1
              ICORDRYD(ND) = 1
            endif
          endif
        endif
      enddo  ! *** END OF LOOP OVER LA
    enddo    ! *** END OF DOMAIN LOOP
    !$OMP END DO
  endif ! *** End if on IDRY = 99

  ! *** PERFORM UPDATE OF P
  !$OMP DO PRIVATE(ND,LF,LL,L)
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    do L = LF,LL
      P(L) = G*(HP(L) + BELV(L))
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  ICORDRY = SUM(ICORDRYD(:))   ! *** OMP flag
  ICORDRY_Global = ICORDRY

  ! *** Now gather sum of ICORDRY
  TTDS = DSTIME(0)
  call MPI_barrier(DSIcomm, ierr)
  TTWAIT = TTWAIT + (DSTIME(0)- TTDS)

  TTDS = DSTIME(0)
  call communicate_ghost_cells(P, 'P')
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)

  call DSI_All_Reduce(ICORDRY, ICORDRY_Global, MPI_SUM, TTDS, 1, TWAIT)
  DSITIMING(11) = DSITIMING(11) + TTDS
  TTWAIT = TTWAIT + TWAIT

  ! *** IF ANY NODE OR SUB-DOMAIN WET/DRY STATUS CHANGED THEN UPDATE VARIABLES AND REPEAT PRESSURE SOLUTION
  if( ICORDRY_Global > 0 )then

    NCORDRY = NCORDRY + 1

    if( NCORDRY < 500 )then
      GOTO 1000
    else
      ! *** WET/DRY MAXIMUM ITERATIONS EXCEEDED

      ! *** save A SNAPSHOT FOR EE
      if( ISPPH == 1 )then
        WRITE (6,*)'THE LATEST MODEL RESULTS HAVE BEEN SAVED TO THE EE LINKAGE'
        call Map_Write_EE_Binary
        if( process_id == master_id )then
          call EE_LINKAGE(-1)
        endif
      endif

      call STOPP('*** NCORDRY > 500 WHEN ISDRY > 0')
    endif

  elseif( ISDRY > 0 )then
    ! ****************************************************************************
    ! ***  IF SUB/SVB CHANGED THEN COMMUNICATE TO OTHER NODES
    TTDS = DSTIME(0)
    call MPI_barrier(DSIcomm, ierr)
    TTWAIT = TTWAIT + (DSTIME(0)- TTDS)

    TTDS = DSTIME(0)
    call Communicate_PUV1
    TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
    ! ****************************************************************************
  endif

  NCORDRYMAX = max(NCORDRY,NCORDRYMAX)
  NCORDRYAVG = NCORDRYAVG + NCORDRY + 1
  NCOUNT = NCOUNT + 1

  ! *** *******************************************************************C
  ! *** FINISHED WITH WETTING/DRYING ITERATIONS
  if( INT(TIMEDAY) /= DAYOLD .or. TIMEDAY >= (TIMEEND-DELT/86400.) )then

    ! *** Write CALPUV.LOG to every process log
    open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')

    if( INT(TBEGIN) == DAYOLD .or. INT(TIMEDAY) >= DAYOLD30 )then
      write(mpi_efdc_out_unit,'(A//)') ' *** CALPUV Interation summary'
      write(mpi_efdc_out_unit,'(A)') '              N     NITER     TIMEDAY             NPUVAVG   NPUVMAX            NCONGAVG  NCONGMAX'
      DAYOLD30 = DAYOLD + 30.
    endif
    write(mpi_efdc_out_unit, '(I15,I10,F12.3,3(10X,2I10))') N, NITER, TIMEDAY, NINT(FLOAT(NCORDRYAVG)/FLOAT(NCOUNT)), NCORDRYMAX, &
      NINT(FLOAT(ITERAVG)/FLOAT(NITERAVG)),  ITERMAX
    close(mpi_efdc_out_unit)

    NCORDRYAVG = 0
    ITERMAX = 0
    ITERAVG = 0
    NITERAVG = 0
    NCOUNT = 0
    DAYOLD = INT(TIMEDAY)
  endif

  ! *** REPORT CELLS THAT HAVE JUST REWETTED
  if( ISDRY > 0 .and. DEBUG )then
    open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
    do L = 2,LA
      if( .not. LMASKDRY(L) )then
        if( HP(L) >= HDRY )then
          ! *** PREVIOUSLY CELL WAS DRY, NOW IT IS WET
          if( ISTL == 3 )then
            write(mpi_efdc_out_unit,'(A,2I5,F12.5,I5,L5,3F10.5)') 'REWETTED CELL: ',N,ISTL,TIMEDAY,Map2Global(L),LMASKDRY(L),H2P(L),HP(L)
          else
            write(mpi_efdc_out_unit,'(A,2I5,F12.5,I5,L5,3F10.5)') 'REWETTED CELL: ',N,ISTL,TIMEDAY,Map2Global(L),LMASKDRY(L),H1P(L),HP(L)
          endif
        endif
      endif
    enddo
    close(mpi_efdc_out_unit)
  endif

  ! *** CHECK IF ISOLATED CELL WATER NEEDS TO BE REMOVED
  if( ISDRY > 0 .and. NDRYSTP > 0 )then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LE,LW,NTMP,IUW,IUE,IVS,IVN,IFACE,K)  &
    !$OMP          PRIVATE(RDRY,BELVAVG,RVAL,HOLDTMP,TMPVAL,SVPW1,RNPORI,ETGWTMP,ETGWAVL) NUM_THREADS(NOPTIMAL(1))
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)

      ! *** Count the number to time steps a cell is isolated, and if it has been
      ! *** isolated for more than ndrystp, and its bottom elevation is higher
      ! *** than the surrounding dry cells, then reduce its depth below the
      ! *** drying depth if necessary.  Save volume reduction rate as QDWASTE
      ! *** Defined as positive out.
      NTMP = NDRYSTP
      do L = LF,LL
        if( HP(L) >= HDRY )then
          ! *** WET CELL, DETERMINE IF ISOLATED
          LE = LEC(L)
          LN = LNC(L)
          RDRY = SUB(L) + SUB(LE) + SVB(L) + SVB(LN)
          if( RDRY > 0.5 )then
            NATDRY(L) = 0                 ! *** CELL IS NOT ISOLATED
          else
            NATDRY(L) = NATDRY(L) + 1     ! *** CELL IS ISOLATED

            if( NATDRY(L) > NTMP )then
              ! *** EXCEEDED THE NUMBER OF ISOLATED STEPS SO DETERMINE IF NEED TO DRY OUT THE CELL
              LW = LWC(L)
              LS = LSC(L)
              BELVAVG = 0.0
              RVAL = 0.0
              if( SUBO(LE) > 0.5 .and. BELV(LE) < BELV(L) )then
                RVAL = RVAL+1.
              endif
              if( SUBO(L)  > 0.5 .and. BELV(LW)<BELV(L) )then
                RVAL = RVAL+1.
              endif
              if( SVBO(LN) > 0.5 .and. BELV(LN)<BELV(L) )then
                RVAL = RVAL+1.
              endif
              if( SVBO(L)  > 0.5 .and. BELV(LS)<BELV(L) )then
                RVAL = RVAL+1.
              endif
              if( RVAL > 0. .or. HP(L) < HDRY*2. )then
                ! *** CHECK FOR OPEN BOUNDARY CELLS

                if( ANY(LOBCS == L) )then
                  NATDRY(L) = 0
                  CYCLE
                endif

                ! *** CHECK FOR FLOW BOUNDARY CELLS
                if( QSUME(L) > 0.0 )then
                  NATDRY(L) = 0
                  CYCLE
                endif

                ! *** SUM VOLUME OF "WASTED" WATER
                HOLDTMP = HP(L)
                HP(L)   = 0.90*HDRY
                P(L)    = G*(HP(L) + BELV(L))

                ! *** UPDATE H1P FOR THE DRY/WASTED CELL
                H1P(L) = HP(L)
                H2P(L) = HP(L)
                do K = KSZ(L),KC
                  H1PK(L,K) = H1P(L)*DZC(L,K)
                  H2PK(L,K) = H1P(L)
                enddo

                NATDRY(L) = 0
                QDWASTE(L) = DELTI*DXYP(L)*(HOLDTMP-HP(L))
                VDWASTE(L) = VDWASTE(L) + DXYP(L)*(HOLDTMP-HP(L))
                NNATDRY(ND) = NNATDRY(ND)+1
                LNATDRY(ND,NNATDRY(ND)) = L
              endif
            endif
          endif
        endif  ! *** END OF WET CELL TEST
      enddo
    enddo      ! *** End of Domain loop
    !$OMP END PARALLEL DO
  endif

  ! ****************************************************************************
  call DSI_All_Reduce(IUPDATE, IUPDATE_Global, MPI_SUM, TTDS, 1, TWAIT)
  DSITIMING(11) = DSITIMING(11) + TTDS
  TTWAIT = TTWAIT + TWAIT

  if( IUPDATE_Global > 0 )then
    TTDS = DSTIME(0)
    call communicate_1D2(P, HP)
    TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  endif
  ! ****************************************************************************


  ! *** *******************************************************************C
  ! *** PERFORM FINAL UPDATES OF HU, AND HV
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LE,LW,NTMP,IUW,IUE,IVS,IVN,IFACE,K)            &
  !$OMP                             PRIVATE(RDRY,BELVAVG,RVAL,HOLDTMP,TMPVAL,SVPW1,RNPORI,ETGWTMP,ETGWAVL)  &
  !$OMP                             NUM_THREADS(NOPTIMAL(1))
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    !IF( ISDRY == 0 )then
    !  do L = LF,LL
    !    P(L) = G*(HP(L)+BELV(L))
    !  enddo
    !ENDIF

    if( IGRIDV == 0 )then
      do L = LF,LL
        LS = LSC(L)
        LW = LWC(L)
        HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )*FSGZUDXYPI(L)
        HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )*FSGZVDXYPI(L)
      enddo
    endif

    ! *** SET TRANSPORT MASK FOR DRY CELLS
    if( ISDRY > 0 )then
      do L = LF,LL
        LMASKDRY(L) = .TRUE.
      enddo
      do L = LF,LL
        ! *** Bypass dry cells unless they are actively have boundary flows
        if( HP(L) < HDRY )then
          if( QSUME(L) /= 0.0 ) CYCLE        ! *** Check if cell is an active boundary
          LE = LEC(L)
          LN = LNC(L)
          IUW = 0
          IUE = 0
          IVS = 0
          IVN = 0
          ! *** THIS REQUIRES THE CELL HAVE HP<HDRY FOR 2 ITERATIONS
          if( SUB1(L)  < 0.5 .and. SUB(L)  < 0.5 ) IUE = 1
          if( SUB1(LE) < 0.5 .and. SUB(LE) < 0.5 ) IUW = 1
          if( SVB1(L)  < 0.5 .and. SVB(L)  < 0.5 ) IVS = 1
          if( SVB1(LN) < 0.5 .and. SVB(LN) < 0.5 ) IVN = 1

          IFACE = IUW + IUE + IVS + IVN
          if( IFACE == 4 )then
            LMASKDRY(L) = .FALSE.
          endif
        endif
      enddo
    endif

    ! *** PERFORM UPDATE ON GROUNDWATER ELEVATION
    if( ISGWIE >= 1 )then
      do L = LF,LL
        QSUM(L,KSZ(L)) = QSUM(L,KSZ(L)) + QGW(L)     ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT)
      enddo

      ! *** INFILTRATION STEP
      RNPORI = 1./RNPOR
      if( ISTL == 3 )then
        do L = LF,LL
          AGWELV(L) = AGWELV2(L) - RNPORI*DELT*DXYIP(L)*QGW(L)
        enddo
      else
        do L = LF,LL
          AGWELV(L) = AGWELV1(L) - RNPORI*DELT*DXYIP(L)*QGW(L)
        enddo
      endif

      do L = LF,LL
        AGWELV(L) = min(AGWELV(L),BELV(L))
      enddo

      ! *** ET STEP
      do L = LF,LL
        if( IEVAP > 1 )then
          SVPW1 = (10.**((0.7859 + 0.03477*TEM(L,KC))/(1. + 0.00412*TEM(L,KC))))
          EVAPT(L) = CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW1 - VPAT(L))/PATMT(L)
        endif
        ETGWTMP = EVAPT(L)-EVAPSW(L)*DXYIP(L)        ! *** EXCESS EVAPORATION
        ETGWTMP = max(ETGWTMP,0.0)
        ETGWAVL = RNPOR*DELTI*(AGWELV(L)-BELAGW(L))
        ETGWAVL = max(ETGWAVL,0.0)
        ETGWTMP = min(ETGWTMP,ETGWAVL)
        EVAPGW(L) = ETGWTMP*DXYP(L)                  ! *** TRANSPIRATION
      enddo

      do L = LF,LL
        AGWELV(L) = AGWELV(L)-RNPORI*DELT*DXYIP(L)*EVAPGW(L)
      enddo

      do L = LF,LL
        AGWELV(L) = max(AGWELV(L),BELAGW(L))
      enddo
    endif

    if( ISDRY > 0 )then
      do K = 1,KC
        do L = LF,LL
          SUB3D(L,K) = SUB(L)*SUB3DO(L,K)
          SVB3D(L,K) = SVB(L)*SVB3DO(L,K)
        enddo
      enddo
    endif

  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO

  ! *** IF ANY CELLS WASTED WATER THEN REPORT
  NTMP = SUM(NNATDRY)
  if( NTMP > 0 )then
    if( INOTICE /= INT(TIMEDAY) )then    ! *** ONLY DISPLAY WARNING ONCE PER DAY
      PRINT '(A,F14.5,I6,i4)',' QDWASTED CELLS. SEE mpi_qdwaste_proc_xxx.log FOR DETAILS. [TIMEDAY, # OF CELLS, Process ID]: ', TIMEDAY, NTMP, process_id
      INOTICE = INT(TIMEDAY)
    endif

    open(mpi_qdwaste_unit,FILE = OUTDIR//mpi_qdwaste,POSITION = 'APPEND')
    do ND = 1,NOPTIMAL(1)
      do I = 1,NNATDRY(ND)
        L = LNATDRY(ND,I)

        TMPVAL = QDWASTE(L)/DXYP(L)
        write(mpi_qdwaste_unit,8888) TIMEDAY, ND, Map2Global(L).IG, Map2Global(L).JG, TIMEDAY, H1P(L), HP(L), QDWASTE(L), TMPVAL

        ! *** ZERO QDWASTE IF NOT CONDUCTING MASS BALANCE
        QDWASTE(L) = 0.0
      enddo
    enddo
    close(mpi_qdwaste_unit)
8888 FORMAT(' QDWASTE ',F12.5,3I6,F12.4,2F10.4,E14.6,F10.4)
  endif

  ! *** CHECK FOR NEGATIVE DEPTHS
  call NEGDEP(NOPTIMAL(1), LDMOPT(1), QCHANUT, QCHANVT, 3, SUB1, SVB1)

  ! ***  CALCULATE THE EXTERNAL DIVERGENCE
  if( ISDIVEX > 0 )then
    DIVEXMX = 0.
    DIVEXMN = 1000000.
    if( ISTL == 3 )then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LN,LE,DIVEX,DIVEXMX,LMAX,LMIN,DIVEXMN) NUM_THREADS(NDM)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          if( SPB(L) /= 0 )then
            LN = LNC(L)
            LE = LEC(L)
            DIVEX = SPB(L)*(DXYP(L)*(HP(L)-H2P(L))*DELTI + 0.5*( UHDYE(LE)+UHDY2E(LE)-UHDYE(L)-UHDY2E(L) &
              +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L) ) - QSUME(L) - QGW(L) + EVAPSW(L))
            if( ISDIVEX == 2 ) DIVEX = DIVEX*DXYIP(L)*HPI(L)  !  *** RELATIVE DIVERGENCE
            if( DIVEX > DIVEXMX )then
              DIVEXMX = DIVEX
              LMAX = L
            endif
            if( DIVEX < DIVEXMN )then
              DIVEXMN = DIVEX
              LMIN = L
            endif
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LN,LE,DIVEX,DIVEXMX,LMAX,LMIN,DIVEXMN) NUM_THREADS(NDM)
      do ND = 1,NDM
        LF = 2+(ND-1)*LDM
        LL = min(LF+LDM-1,LA)
        do L = LF,LL
          if( SPB(L) /= 0 )then
            LN = LNC(L)
            LE = LEC(L)
            DIVEX = SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI + 0.5*( UHDYE(LE)+UHDY1E(LE)-UHDYE(L)-UHDY1E(L) &
              +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L) ) - QSUME(L) - QGW(L) + EVAPSW(L))
            if( ISDIVEX == 2 ) DIVEX = DIVEX*DXYIP(L)*HPI(L)  !  *** RELATIVE DIVERGENCE
            if( DIVEX > DIVEXMX )then
              DIVEXMX = DIVEX
              LMAX = L
            endif
            if( DIVEX < DIVEXMN )then
              DIVEXMN = DIVEX
              LMIN = L
            endif
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif
    IMAX = IL(LMAX)
    JMAX = JL(LMAX)
    IMIN = IL(LMIN)
    JMIN = JL(LMIN)
    write(6,6628)DIVEXMX,IMAX,JMAX
    write(6,6629)DIVEXMN,IMIN,JMIN

6628 FORMAT('  DIVEXMX = ',E13.5,5X,2I10)
6629 FORMAT('  DIVEXMN = ',E13.5,5X,2I10)
  endif  ! *** End DIVEX

  ! *** DETERMINE THE WET AND DRY CELL LIST
  if( ISDRY == 0 )then
    LAWET = LA-1
    LADRY = 0
  else
    LAWET = 0
    LADRY = 0
    do L = 2,LA
      if( LMASKDRY(L) )then
        LAWET = LAWET + 1
        LWET(LAWET) = L
        NWET(L) = NWET(L) + 1
      elseif( OLDMASK(L) .or. N < NTSTBC+1 )then
        ! *** ONLY FLAG NEWLY DRY CELLS
        LADRY = LADRY + 1
        LDRY(LADRY) = L

        ! *** Reset wet counter for adjacent cells
        do LL = 1,9
          NWET(LADJ(LL,L)) = 0
        enddo

        ! *** UPDATE DEPTHS FOR DRY CELLS
        do K = KSZ(L),KC
          HPK(L,K)  = HP(L)*DZC(L,K)
        enddo
        HU(L) = ( DXYP(L)*HP(L) + DXYP(LWC(L))*HP(LWC(L)) )*FSGZUDXYPI(L)
        HV(L) = ( DXYP(L)*HP(L) + DXYP(LSC(L))*HP(LSC(L)) )*FSGZVDXYPI(L)
      endif
    enddo
  endif

  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then
    if( ISBEDMAP > 0 )then
      LASED = 0
      LSED = 0
      do L = 1,LAWET
        if( BEDMAP(LWET(L)) > 0 )then
          LASED = LASED + 1
          LSED(LASED) = LWET(L)
        endif
      enddo
    else
      LASED = LAWET
      if( ISDRY > 0 .or. NITER < 5 ) LSED = LWET
    endif
  endif

  LDMWET = INT(FLOAT(LAWET)/FLOAT(NDM))+1
  LDMDRY = INT(FLOAT(LADRY)/FLOAT(NDM))+1
  LDMSED = INT(FLOAT(LASED)/FLOAT(NDM))+1

  ! *** GET CELL LIST FOR ACTIVE LAYERS = 1
  if( IGRIDV > 0 .and. KMINV == 1 )then
    LASGZ1 = 0
    do L = 2,LA
      if( KSZ(L) == KC )then
        if( LMASKDRY(L) )then
          LASGZ1 = LASGZ1+1
          LSGZ1(LASGZ1) = L
        endif
      endif
    enddo
    LDMSGZ1 = INT(LASGZ1/NDM)+1
  endif

  ! *** COMPUTATIONAL CELL LIST FOR EACH SUB-DOMAIN
  if( ISDRY > 0 )then
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LDMWET,LAWET,LWET,LKSZ,LLWET,LKWET,KS,KC)  &
    !$OMP                           PRIVATE(ND,LF,LL,K,LN,LP,L) NUM_THREADS(NDM)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)

      do K = 1,KC
        LN = 0
        do LP = LF,LL
          L = LWET(LP)
          if( LKSZ(L,K) ) CYCLE
          LN = LN+1
          LKWET(LN,K,ND) = L
        enddo
        LLWET(K,ND) = LN    ! *** NUMBER OF WET CELLS FOR THE CURRENT LAYER
      enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LLWET,LLWETZ,LKWET,LKWETZ,KS,KC) PRIVATE(ND,K,LP,L) NUM_THREADS(NDM)
    do ND = 1,NDM
      do K = 1,KS
        LLWETZ(K,ND) = LLWET(K,ND)
        do LP = 1,LLWET(K,ND)
          LKWETZ(LP,K,ND) = LKWET(LP,K,ND)
        enddo
      enddo

      LLWETZ(KC,ND) = LLWET(KS,ND)
      do LP = 1,LLWET(KS,ND)
        LKWETZ(LP,KC,ND) = LKWET(LP,KS,ND)
      enddo
    enddo
    !$OMP END PARALLEL DO

  endif

  ! *** THIRD PASS CELL CONSTANTS
  if( IGRIDV > 0 )then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,K,LP,L) NUM_THREADS(NDM)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = min(LF+LDMWET-1,LAWET)

      ! *** GLOBAL UPDATE OF LAYER THICKNESSES
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          HPK(L,K)  = HP(L)*DZC(L,K)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! *** UPDATE HU & HV FOR ALL CELLS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,K,LP,L,LW,LS) NUM_THREADS(NOPTIMAL(1))
    do ND = 1,NOPTIMAL(1)
      LF = 2+(ND-1)*LDMOPT(1)
      LL = min(LF+LDMOPT(1)-1,LA)
      do L = LF,LL
        LW = LWC(L)
        LS = LSC(L)

        if( KSZ(LW) > KSZ(L) )then
          HU(L) = max( 0.5*HPK(L,KSZ(LW)), HP(LW)*(1.+DZC(L,KSZ(LW))*0.1) )
        elseif( KSZ(LW) < KSZ(L) )then
          HU(L) = max( 0.5*HPK(LW,KSZ(L)), HP(L)*(1.+DZC(LW,KSZ(L))*0.1) )
        else
          HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )*FSGZUDXYPI(L)
        endif

        if( KSZ(LS) > KSZ(L) )then
          HV(L) = max( 0.5*HPK(L,KSZ(LS)), HP(LS)*(1.+DZC(L,KSZ(LS))*0.1) )
        elseif( KSZ(LS) < KSZ(L) )then
          HV(L) = max( 0.5*HPK(LS,KSZ(L)), HP(L)*(1.+DZC(LS,KSZ(L))*0.1) )
        else
          HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )*FSGZVDXYPI(L)
        endif
      enddo

    enddo
    !$OMP END PARALLEL DO
  endif

  ! *** BLOCKED LAYER FACE OPTION
  if( NBLOCKED > 0 )then
    do LP = 1,NBLOCKED
      L = LBLOCKED(LP)
      if( KSZ(L) == KC ) CYCLE

      ! *** U FACE BLOCKING
      if( BLDRAFTUO(LP)+BLSILLU(LP) > 0.0 )then
        if( BLANCHORU(LP) /= 0. )then
          BLDRAFTU(LP) = BLDRAFTUO(LP) - (BLANCHORU(LP) - (BELV0(L) + HP(L)))
          BLDRAFTU(LP) = max(BLDRAFTU(LP),0.0)
        else
          BLDRAFTU(LP) = BLDRAFTUO(LP)
        endif

        ! *** RESET TO DEFAULT LAYER
        KSZU(L) = max(KSZ(L), KSZ(LWC(L)))      ! *** MINUMUM ACTIVE LAYERS FOR U FACE
        do K = 1,KC
          SUB3D(L,K) = 0.0
          if( K >= KSZU(L) )then
            SUB3D(L,K)  = 1.0
            SUB3DO(L,K) = 1.0
            if( IGRIDV > 0 )then
              if( KSZ(LWC(L)) > KSZ(L) )then
                SGZU(L,K)  = DZC(LSC(L),K)
              else
                SGZU(L,K)  = DZC(L,K)
              endif
            else
              SGZU(L,K)  = max(DZC(LWC(L),K),DZC(L,K))
            endif
          endif
        enddo

        HU(L) = HU(L) - BLDRAFTU(LP) - BLSILLU(LP)

        if( HU(L) < 0.0 )then
          ! *** Flows fully blocked
          SUB(L)  = 0.0
          SUBO(L) = 0.0
          SAAX(L) = 0.0
          SUB3D(L,KSZ(L):KC)  = 0.0
          SUB3DO(L,KSZ(L):KC) = 0.0
          HU(L) = HDRY
        else
          ! *** Flows partially blocked
          SUB(L)  = 1.0
          SUBO(L) = 1.0
          SAAX(L) = 1.0

          HU(L) = max(HU(L),HWET)
          HUI(L) = 1./HU(L)

          ! *** BLOCK CELL FACES FROM TOP
          if( BLDRAFTU(LP) > 0.0 )then
            TMPX = 0.
            do K = KC,KSZU(L)+1,-1
              TMPX = TMPX + SGZU(L,K)*HP(L)
              KTBU(LP) = K - 1
              SUB3D(L,K) = 0.0
              SGZU(L,K) = 0.0
              if( TMPX > BLDRAFTU(LP) ) exit
            enddo
            KTBU(LP) = max(KTBU(LP),KSZU(L))
          else
            KTBU(LP) = KC
          endif

          ! *** BLOCK CELL FACES FROM BOTTOM
          if( BLSILLU(LP) > 0.0 )then
            TMPX = 0.
            do K = KSZU(L),KC-1
              TMPX = TMPX + SGZU(L,K)*HP(L)
              KBBU(LP) = K + 1
              SUB3D(L,K) = 0.0
              SGZU(L,K) = 0.0
              if( TMPX > BLSILLU(LP) ) exit
            enddo
            KBBU(LP) = min(KBBU(LP),KTBU(LP))
          else
            KBBU(LP) = KSZU(L)
          endif
          KSZU(L) = KBBU(LP)

          TMPX = SUM(SGZU(L,1:KC))
          do K = 1,KC
            SGZU(L,K) = SGZU(L,K) / TMPX
            DZGU(L,K) = 0.0
            CDZFU(L,K) = 0.0          ! *** USED FOR GLOBAL DU SHEAR
            CDZUU(L,K) = 0.0          ! *** USED FOR SURFACE SHEAR DUE TO WIND
            CDZLU(L,K) = 0.0
            CDZMU(L,K) = 0.0
            CDZRU(L,K) = 0.0
            CDZDU(L,K) = 0.0
          enddo

          do K = KBBU(LP),KTBU(LP)-1
            DZGU(L,K) = 0.5*(SGZU(L,K)+SGZU(L,K+1))
            CDZFU(L,K) = SGZU(L,K)*SGZU(L,K+1)/(SGZU(L,K) + SGZU(L,K+1))
            CDZUU(L,K) = -SGZU(L,K)  /(SGZU(L,K) + SGZU(L,K+1))
            CDZLU(L,K) = -SGZU(L,K+1)/(SGZU(L,K) + SGZU(L,K+1))
          enddo

          CDZRU(L,KBBU(LP)) = SGZU(L,KBBU(LP)) - 1.
          CDZDU(L,KBBU(LP)) = SGZU(L,KBBU(LP))
          do K = KBBU(LP)+1,KTBU(LP)-1
            CDZRU(L,K) = CDZRU(L,K-1) + SGZU(L,K)
            CDZDU(L,K) = CDZDU(L,K-1) + SGZU(L,K)
          enddo

          do K = KBBU(LP),KTBU(LP)-1
            CDZRU(L,K) = CDZRU(L,K)*DZGU(L,K)*CDZLU(L,KSZU(L))
            CDZMU(L,K) = 0.5*SGZU(L,K)*SGZU(L,K+1)
          enddo
        endif
      endif

      ! *** V FACE BLOCKING
      if( BLDRAFTVO(LP)+BLSILLV(LP) > 0.0 )then
        if( BLANCHORV(LP) /= 0. )then
          BLDRAFTV(LP) = BLDRAFTVO(LP) - (BLANCHORV(LP) - (BELV0(L) + HP(L)))
          BLDRAFTV(LP) = max(BLDRAFTV(LP),0.0)
        else
          BLDRAFTV(LP) = BLDRAFTVO(LP)
        endif

        ! *** RESET TO DEFAULT LAYER
        KSZV(L) = max(KSZ(L), KSZ(LSC(L)))      ! *** MINUMUM ACTIVE LAYERS FOR V FACE
        do K = 1,KC
          SVB3D(L,K)  = 0.0
          SVB3DO(L,K) = 0.0
          if( K >= KSZV(L) )then
            SVB3D(L,K)  = 1.0
            SVB3DO(L,K) = 1.0
            if( IGRIDV > 0 )then
              if( KSZ(LSC(L)) > KSZ(L) )then
                SGZV(L,K)  = DZC(LSC(L),K)
              else
                SGZV(L,K)  = DZC(L,K)
              endif
            else
              SGZV(L,K)  = max(DZC(LSC(L),K),DZC(L,K))
            endif
          endif
        enddo

        HV(L) = HV(L) - BLDRAFTV(LP) - BLSILLV(LP)

        if( HV(L) < 0.0 )then
          ! *** Flows fully blocked
          SVB(L)  = 0.0
          SVBO(L) = 0.0
          SAAY(L) = 0.0
          SVB3D(L,KSZ(L):KC)  = 0.0
          SVB3DO(L,KSZ(L):KC) = 0.0
          HV(L) = HDRY
        else
          ! *** Flows partially blocked
          SVB(L) = 1.0
          SVBO(L) = 1.0
          SAAY(L) = 1.0

          HV(L) = max(HV(L),HWET)
          HVI(L) = 1./HV(L)

          ! *** BLOCK CELL FACES FROM TOP
          if( BLDRAFTV(LP) > 0.0 )then
            TMPY = 0.
            do K = KC,KSZV(L)+1,-1
              TMPY = TMPY + SGZV(L,K)*HP(L)
              KTBV(LP) = K - 1
              SVB3D(L,K) = 0.0
              SGZV(L,K) = 0.0
              if( TMPY > BLDRAFTV(LP) ) exit
            enddo
            KTBV(LP) = max(KTBV(LP),KSZV(L))
          else
            KTBV(LP) = KC
          endif

          ! *** BLOCK CELL FACES FROM BOTTOM
          if( BLSILLV(LP) > 0.0 )then
            TMPX = 0.
            do K = KSZV(L),KC-1
              TMPX = TMPX + SGZV(L,K)*HP(L)
              KBBV(LP) = K + 1
              SVB3D(L,K) = 0.0
              SGZV(L,K) = 0.0
              if( TMPX > BLSILLV(LP) ) exit
            enddo
            KBBV(LP) = min(KBBV(LP),KTBV(LP))
          else
            KBBV(LP) = KSZV(L)
          endif
          KSZV(L) = KBBV(LP)

          TMPY = SUM(SGZV(L,1:KC))
          do K = 1,KC
            SGZV(L,K) = SGZV(L,K) / TMPY
            DZGV(L,K) = 0.0
            CDZFV(L,K) = 0.0
            CDZUV(L,K) = 0.0
            CDZLV(L,K) = 0.0
            CDZRV(L,K) = 0.0
            CDZDV(L,K) = 0.0
            CDZMV(L,K) = 0.0
          enddo

          do K = KBBV(LP),KTBV(LP)-1
            DZGV(L,K) = 0.5*(SGZV(L,K)+SGZV(L,K+1))
            CDZFV(L,K) = SGZV(L,K)*SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
            CDZUV(L,K) = -SGZV(L,K)  /(SGZV(L,K)+SGZV(L,K+1))
            CDZLV(L,K) = -SGZV(L,K+1)/(SGZV(L,K)+SGZV(L,K+1))
          enddo

          CDZRV(L,KBBV(LP)) = SGZV(L,KBBV(LP))-1.
          CDZDV(L,KBBV(LP)) = SGZV(L,KBBV(LP))
          do K = KBBV(LP)+1,KTBV(LP)-1
            CDZRV(L,K) = CDZRV(L,K-1)+SGZV(L,K)
            CDZDV(L,K) = CDZDV(L,K-1)+SGZV(L,K)
          enddo

          do K = KBBV(LP),KTBV(LP)-1
            CDZRV(L,K) = CDZRV(L,K)*DZGV(L,K)*CDZLV(L,KSZV(L))
            CDZMV(L,K) = 0.5*SGZV(L,K)*SGZV(L,K+1)
          enddo
        endif
      endif
    enddo
  endif

  ! *** COMPUTATIONAL CELL LIST FOR ENTIRE DOMAIN
  if( ISDRY > 0 )then
    do K = 1,KC
      LN = 0
      do L = 2,LA
        if( LKSZ(L,K) ) CYCLE
        if( LMASKDRY(L) )then
          LN = LN+1
          LKWET(LN,K,0) = L   ! *** Wet Cell for Layer K
        endif
      enddo
      LLWET(K,0) = LN        ! *** Total Wet Cells for Layer K
    enddo
  endif

  ! ****************************************************************************
  TTDS = DSTIME(0)
  call MPI_barrier(DSIcomm, ierr)
  TTWAIT = TTWAIT + (DSTIME(0)- TTDS)

  TTDS = DSTIME(0)
  call Communicate_PUV3
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  ! ****************************************************************************

  ! *** Update variables that are dependent on MPI communicated values
  !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NOPTIMAL, LDMOPT, LA, HP, HPI, HU, HUI, HV, HVI) PRIVATE(ND,LF,LL,L) NUM_THREADS(NOPTIMAL(1))
  do ND = 1,NOPTIMAL(1)
    LF = 2+(ND-1)*LDMOPT(1)
    LL = min(LF+LDMOPT(1)-1,LA)

    do L = LF,LL
      HPI(L) = 1./HP(L)
      HUI(L) = 1./HU(L)
      HVI(L) = 1./HV(L)
    enddo
  enddo
  !$OMP END PARALLEL DO

  ! *** Find optimal number of threads
  TTCON = TCONG - TTCON
  call HYD_THREADS(0, TTCON)

  TTDS = DSTIME(0)
  TTPUV = TTDS - TTPUV - TTCON
  call HYD_THREADS(1, TTPUV)

  return

  END

