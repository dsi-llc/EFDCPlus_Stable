! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Calculates internal solution at time level (N+1)
! @author Paul Craig
! @date    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3
!!         2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
!!         2011-03       Paul M. Craig     Rewritten to F90 and added OMP
! @parameter ISTL ISTL indicates the # of time levels in the step

SUBROUTINE CALUVW

  use GLOBAL
  use Variables_MPI
  use Allocate_Initialize      
  use MPI 
  use Communicate_Ghost_Routines

  implicit none

  integer :: NMD, ND, LF, LL, L, LS, LN, K, LE, LW, LNN, LP, LHOST, LCHNU, LCHNV, LG
  integer :: ICFL, JCFL, KCFL, IVAL, IDTCFL, LTMP, L1P, IOBC, IFILE
  integer, save :: NKCE, NKCN, NSTEP, IUPDATE_UVHE

  real      :: Q1, Q2, CMU, CMV, CRU, CRV, EU, EV, RCDZM, RCDZU, RCDZL, HPPTMP, TMPVALN, DTMAXX, HPDEL
  real      :: RLAMN, RLAMO, TMPVAL, CFLUUUT, CFLVVVT, CFLWWWT, CFLCACT, DTCFL, UWTMP, UETMP, VSTMP, VNTMP, WBTMP, WTTMP, DTMAXI
  real      :: STOKESU, STOKESX, STOKESY, STOKESW, C1, C2, THETA, AMPL, SINHH

  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS                 ! MODEL TIMING TEMPORARY VARIABLE

  integer,save,allocatable,dimension(:) :: LKCE
  integer,save,allocatable,dimension(:) :: LKCN
  integer,save,allocatable,dimension(:) :: LHOLE
  integer,save,allocatable,dimension(:) :: LSTEP

  real,save,allocatable,dimension(:)   :: AAU           
  real,save,allocatable,dimension(:)   :: AAV           
  real,save,allocatable,dimension(:)   :: BBU           
  real,save,allocatable,dimension(:)   :: BBV           
  real,save,allocatable,dimension(:)   :: UHDYEE
  real,save,allocatable,dimension(:)   :: VHDXEE
  real,save,allocatable,dimension(:)   :: RCXX
  real,save,allocatable,dimension(:)   :: RCYY
  real,save,allocatable,dimension(:)   :: TVARE
  real,save,allocatable,dimension(:)   :: TVARN
  real,save,allocatable,dimension(:,:) :: CERRU
  real,save,allocatable,dimension(:,:) :: CERRV
  real,save,allocatable,dimension(:)   :: SINH2
  real,save,allocatable,dimension(:)   :: SINH4

  integer :: ierr

  if( .not. allocated(RCXX) )then
    call AllocateDSI( AAU,    LCM, 0.0)     
    call AllocateDSI( AAV,    LCM, 0.0)     
    call AllocateDSI( BBU,    LCM, 0.0)     
    call AllocateDSI( BBV,    LCM, 0.0)     
    call AllocateDSI( LKCE,   LCM,   0)
    call AllocateDSI( LKCN,   LCM,   0)
    call AllocateDSI( LHOLE,  LCM,   0)
    call AllocateDSI( LSTEP,  LCM,   0)
    call AllocateDSI( UHDYEE, LCM, 0.0)
    call AllocateDSI( VHDXEE, LCM, 0.0)
    call AllocateDSI( RCXX,   LCM, 0.0)
    call AllocateDSI( RCYY,   LCM, 0.0)
    call AllocateDSI( TVARE,  LCM, 0.0)
    call AllocateDSI( TVARN,  LCM, 0.0)
    call AllocateDSI( CERRU,  LCM, KCM, 0.0)
    call AllocateDSI( CERRV,  LCM, KCM, 0.0)

    if( (ISWAVE == 2 .or. ISWAVE == 4) .and. ISWVSD >= 1 )then
      call AllocateDSI(SINH2, LCM, 0.0)
      call AllocateDSI(SINH4, LCM, 0.0)
    endif
    
    do L = 2,LA
      ! *** CERRU
      if( SGZU(L,KSZU(L)) == 1.0 )then
        CERRU(L,KSZU(L)) = 1.  ! *** PREVENTS DEVISION BY ZERO FOR BLOCKING OTPION
      elseif( KSZU(L) < KC )then
        Q1 = 0.
        do K = KSZU(L),KC
          Q1 = Q1 + 1./(1.-SGZU(L,K))
        enddo
        if( Q1 > 0. )then
          do K = KSZU(L),KC
            CERRU(L,K) = (1./(1.-SGZU(L,K)))/Q1
            CERRU(L,K) = CERRU(L,K)*REAL(KC-KSZU(L)+1)
          enddo
        endif
      else
        CERRU(L,KC) = 1.
      endif

      ! *** CERRV
      if( SGZV(L,KSZV(L)) == 1.0 )then
        CERRV(L,KSZV(L)) = 1.    ! *** PREVENTS DEVISION BY ZERO FOR BLOCKING OTPION
      elseif( KSZV(L) < KC )then
        Q1 = 0.
        do K = KSZV(L),KC
          Q1 = Q1 + 1./(1.-SGZV(L,K))
        enddo
        if( Q1 > 0. )then
          do K = KSZV(L),KC
            CERRV(L,K) = (1./(1.-SGZV(L,K)))/Q1
            CERRV(L,K) = CERRV(L,K)*REAL(KC-KSZV(L)+1)
          enddo
        endif
      else
        CERRV(L,KC) = 1.
      endif
    enddo

    ! *** FACE BLOCKING OPTION
    if( NBLOCKED > 0 )then
      do LP = 1,NBLOCKED
        L = LBLOCKED(LP)
        if( KSZ(L) == KC ) CYCLE

        ! *** CERRU
        if( BLDRAFTU(LP)+BLSILLU(LP) > 0.0 )then
          CERRU(L,:) = 0.0
          if( SGZU(L,KBBU(LP)) == 1.0 )then
            CERRU(L,KBBU(LP)) = 1.  ! *** PREVENTS DEVISION BY ZERO FOR BLOCKING OTPION
          elseif( KBBU(LP) < KC )then
            Q1 = 0.
            do K = KBBU(LP),KTBU(LP)
              Q1 = Q1 + 1./(1.-SGZU(L,K))
            enddo
            if( Q1 > 0. )then
              do K = KBBU(LP),KTBU(LP)
                CERRU(L,K) = (1./(1.-SGZU(L,K)))/Q1
                CERRU(L,K) = CERRU(L,K)*REAL(KTBU(LP)-KBBU(LP)+1)
              enddo
            endif
          else
            CERRU(L,KTBU(LP)) = 1.
          endif
        endif

        ! *** CERRV
        if( BLDRAFTV(LP)+BLSILLV(LP) > 0.0 )then
          CERRV(L,:) = 0.0
          if( SGZV(L,KBBV(LP)) == 1.0 )then
            CERRV(L,KBBV(LP)) = 1.  ! *** PREVENTS DEVISION BY ZERO FOR BLOCKING OTPION
          elseif( KBBV(LP) < KC )then
            Q1 = 0.
            do K = KBBV(LP),KTBV(LP)
              Q1 = Q1 + 1./(1.-SGZV(L,K))
            enddo
            if( Q1 > 0. )then
              do K = KBBV(LP),KTBV(LP)
                CERRV(L,K) = (1./(1.-SGZV(L,K)))/Q1
                CERRV(L,K) = CERRV(L,K)*REAL(KTBV(LP)-KBBV(LP)+1)
              enddo
            endif
          else
            CERRV(L,KTBV(LP)) = 1.
          endif
        endif

      enddo
    endif

    ! *** HANDLE KSZ E&N = KC CASES
    NKCE = 0
    NKCN = 0
    do L = 2,LA
      if( KSZ(LWC(L)) < KC .and. SUBO(L) > 0.5 .and. KSZ(L) == KC )then
        NKCE = NKCE+1
        LKCE(NKCE) = L
      endif
      if( KSZ(LSC(L)) < KC .and. SVBO(L) > 0.5 .and. KSZ(L) == KC  )then
        NKCN = NKCN+1
        LKCN(NKCN) = L
      endif
    enddo

    ! *** BUILD A LIST OF ALL THE STEPS BETWEEN KC AND KS BOTTOM LAYERS
    NSTEP = 0
    if( IGRIDV > 0 .and. KMINV == 1 )then
      do L = 2,LA
        if( KSZ(L) == KS )then
          K = KS
          Q1 = SGZU(L,K) + SGZU(LEC(L),K) + SGZV(L,K) + SGZV(LNC(L),K)
          Q2 = 0.

          ! *** EXCLUDE HOLES
          if( Q1 /= 0. )then

            ! *** ONE OR MORE U FACES ARE ACTIVE
            if( ( SUBO(LEC(L)) > 0. .and. KSZ(LEC(L)) == KC ) .or. ( SUBO(L) > 0. .and. KSZ(LWC(L)) == KC ) )then
              Q2 = Q2 + 1.
            endif

            ! *** ONE OR MORE V FACES ARE ACTIVE
            if( ( SVBO(LNC(L)) > 0. .and. KSZ(LNC(L)) == KC ) .or. ( SVBO(L) > 0. .and. KSZ(LSC(L)) == KC ) )then
              Q2 = Q2 + 1.
            endif
            if( Q2 > 0. )then
              NSTEP = NSTEP+1
              LSTEP(NSTEP) = L
            endif
          endif
        endif
      enddo
    endif
    
    ! *** Set skip flag for updating UHE and VHE
    IUPDATE_UVHE = 0
    if( ISDYNSTP > 0 .or. ISITB >= 1 .or. LSEDZLJ .or. (ISTRAN(5) > 0 .and. ISTRAN(2) > 0 .and. ANY(ITOXKIN(3,:) > 0)) ) IUPDATE_UVHE = 1
      
    ! ****************************************************************************
    call MPI_barrier(MPI_Comm_World, ierr)
    call communicate_ghost_cells(CERRU)
    call communicate_ghost_cells(CERRV)
    ! ****************************************************************************

  endif

  if( ISDYNSTP == 0 )then
    DELT = DT2
    if( ISTL  ==  2 )then
      DELT = DT
    endif
    DELTI = 1./DELT
  else
    DELT = DTDYN
    DELTI = 1./DELT
  endif
  IFILE = -1

  ! *** ZERO NEWLY DRY CELLS
  if( LADRY > 0 )then
    do LP = 1,LADRY
      L = LDRY(LP)
      AAU(L) = 0.0
      AAV(L) = 0.0
      BBU(L) = 1.0
      BBV(L) = 1.0
      RCXX(L) = 0.0
      RCYY(L) = 0.0
      TVARE(L) = 0.0
      TVARN(L) = 0.0
    enddo

    do K = 1,KC
      do LP = 1,LADRY
        L = LDRY(LP)
        DU(L,K) = 0.0
        DV(L,K) = 0.0
        U(L,K)  = 0.0
        V(L,K)  = 0.0
        W(L,K)  = 0.0
        UUU(L,K) = 0.0
        VVV(L,K) = 0.0
        UHDY(L,K) = 0.0
        VHDX(L,K) = 0.0
        UHDYF(L,K) = 0.0
        VHDXF(L,K) = 0.0
        UHDYEK(L,K) = 0.0
        VHDXEK(L,K) = 0.0
      enddo
    enddo
  endif

  if( KC == 1 ) GOTO 30

  if( IGRIDV > 0 .and. KMINV == 1 )then
    ! *** HANDLE SINGLE LAYER CELLS
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LASGZ1,LDMSGZ1,LSGZ1,KSZ,UHDYE,UHDYF,UHDY,HUI,DYIU,U,VHDXE,VHDXF,VHDX,HVI,DXIV,V,W) PRIVATE(ND,LF,LL,LP,L)
    do ND = 1,NDM
      LF = (ND-1)*LDMSGZ1+1
      LL = MIN(LF+LDMSGZ1-1,LASGZ1)

      do LP = LF,LL
        L = LSGZ1(LP)
        UHDYF(L,KSZ(L)) = UHDYE(L)
        UHDY(L,KSZ(L))  = UHDYE(L)
        U(L,KSZ(L))     = UHDYE(L)*HUI(L)*DYIU(L)
        VHDXF(L,KSZ(L)) = VHDXE(L)
        VHDX(L,KSZ(L))  = VHDXE(L)
        V(L,KSZ(L))     = VHDXE(L)*HVI(L)*DXIV(L)
        W(L,KSZ(L))     = 0.0
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif

  !$OMP PARALLEL DEFAULT(SHARED)

  ! ***************************************************************************
  ! *** CALCULATE BOTTOM FRICTION COEFFICIENT
  if( ISTL == 3 )then
    !$OMP DO PRIVATE(ND,LP,L)
    do ND = 1,NDM
      do LP = 1,LLWETZ(KC,ND)
        L = LKWETZ(LP,KC,ND)
        RCXX(L) = STBX(L)*SQRT( V1U(L)*V1U(L) + U1(L,KSZU(L))*U1(L,KSZU(L)) )
        RCYY(L) = STBY(L)*SQRT( U1V(L)*U1V(L) + V1(L,KSZV(L))*V1(L,KSZV(L)) )
      enddo
    enddo
    !$OMP END DO
  else
    if( AVCON1 < 0.00001 .or. ICALTB == 0 )then
      ! *** FOR 2TL U1 & U AND V1 & V ARE THE SAME
      ! *** THESE ARE ONLY DIFFERENCE FOR 3TL ISTL = 2 TRAP CORRECTION STEP

      !$OMP DO PRIVATE(ND,LP,L,Q1,Q2)
      do ND = 1,NDM
        do LP = 1,LLWETZ(KC,ND)
          L = LKWETZ(LP,KC,ND)
          Q1      = SQRT( U1(L,KSZU(L))*U1(L,KSZU(L)) + V1U(L)*V1U(L) )
          Q2      = SQRT( U(L,KSZU(L)) *U(L,KSZU(L))  + VU(L) *VU(L) )
          RCXX(L) = STBX(L)*SQRT(Q1*Q2)
          Q1      = SQRT( V1(L,KSZV(L))*V1(L,KSZV(L)) + U1V(L)*U1V(L) )
          Q2      = SQRT( V(L,KSZV(L)) *V(L,KSZV(L))  + UV(L) *UV(L) )
          RCYY(L) = STBY(L)*SQRT(Q1*Q2)
        enddo
      enddo
      !$OMP END DO
    else
      ! *** CONSTANT AVO
      !$OMP DO PRIVATE(ND,LP,L,Q1,Q2)
      do ND = 1,NDM
        do LP = 1,LLWETZ(KC,ND)
          L = LKWETZ(LP,KC,ND)
          RCXX(L) = AVCON1/SQRT(H1U(L)*HU(L))
          RCYY(L) = AVCON1/SQRT(H1V(L)*HV(L))
        enddo
      enddo
      !$OMP END DO
    endif

  endif

  ! ***************************************************************************
  ! *** SPLIT THE DEPTH AVERAGED FLOWS INTO THE LAYER SPECIFIC FLOWS
  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN)
  do ND = 1,NDM
    do K = 1,KS
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        LN = LNC(L)
        LE = LEC(L)
        UHDYEK(L,K) = DZC(L,K)*( UHDYE(L) - UHDYE(LE) )
        VHDXEK(L,K) = DZC(L,K)*( VHDXE(L) - VHDXE(LN) )
      enddo
    enddo
  enddo
  !$OMP END DO

  ! *** MAKE ADJUSTMENTS FOR SIGMA-ZED STEPS AND HOLES
  if( IGRIDV > 0 )then
    ! *** HANDLE KS-KC STEPS
    !$OMP SINGLE
    if( KMINV == 1 )then
      do LP = 1,NSTEP
        L = LSTEP(LP)
        LN = LNC(L)
        LE = LEC(L)

        K = KS
        HPDEL = HP(L)-H1P(L)
        Q1 = DZC(L,K)*HPDEL*DXYP(L)/DELT - DZC(L,K)*QSUME(L)                           ! *** FLUX REQUIRED DUE TO HPK CHANGE
        Q2 = -0.5*DZC(L,K)*( UHDYE(LE)+UHDY1E(LE) - UHDYE(L)-UHDY1E(L) + &
                             VHDXE(LN)+VHDX1E(LN) - VHDXE(L)-VHDX1E(L) - QSUME(L) )    ! *** FLUX DUE TO HORIZONTAL IN/OUT FLOWS
        if( HPDEL > 0. )then
          UHDYEK(L,K) = MAX(Q1,Q2)
        else
          UHDYEK(L,K) = MIN(Q1,Q2)
        endif
        VHDXEK(L,K) = 0.0
      enddo
    endif
    !$OMP END SINGLE
  endif

  ! ***************************************************************************
  ! *** CALCULATE THE U AND V SHEARS
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LE,LW,LS,LN,RCDZL,RCDZM,RCDZU,CMU,EU,CMV,EV,CRU,CRV)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET+1
    LL = MIN(LF+LDMWET-1,LAWET)

    do K = 1,KS
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        DU(L,K) = SUB3D(L,K)*DU(L,K)
        if( K > KSZU(L) )then
          RCDZL = CDZLU(L,K)
          RCDZM = CDZMU(L,K)*DELTI
          RCDZU = CDZUU(L,K)
          CMU = 1.+RCDZM*HU(L)*AVUI(L,K)
          EU  = 1./(CMU-RCDZL*CU1(L,K-1))
          CU1(L,K) = RCDZU*EU
          DU(L,K)  = (DU(L,K)-RCDZL*DU(L,K-1))*EU
          UUU(L,K) = -RCDZL*UUU(L,K-1)*EU
        elseif( K == KSZU(L) )then
          ! *** K = KSZU(L)
          RCDZL = CDZLU(L,K)
          RCDZM = CDZMU(L,K)*DELTI
          RCDZU = CDZUU(L,K)
          CMU = 1.+RCDZM*HU(L)*AVUI(L,K)
          EU  = 1./CMU
          CU1(L,K) = RCDZU*EU
          DU(L,K)  = (DU(L,K)-RCDZL*RCXX(L)*UHE(L)*HUI(L))*EU
          UUU(L,K) = EU
        else
          CU1(L,K) = 0.
          DU(L,K) = 0.
          UUU(L,K) = 0.
        endif

        DV(L,K) = SVB3D(L,K)*DV(L,K)
        if( K > KSZV(L) )then
          RCDZL = CDZLV(L,K)
          RCDZM = CDZMV(L,K)*DELTI
          RCDZU = CDZUV(L,K)
          CMV = 1.+RCDZM*HV(L)*AVVI(L,K)
          EV  = 1./(CMV-RCDZL*CU2(L,K-1))
          CU2(L,K) = RCDZU*EV
          DV(L,K)  = (DV(L,K)-RCDZL*DV(L,K-1))*EV
          VVV(L,K) = -RCDZL*VVV(L,K-1)*EV
        elseif( K == KSZV(L) )then
          ! *** K = KSZV(L)
          RCDZL = CDZLV(L,K)
          RCDZM = CDZMV(L,K)*DELTI
          RCDZU = CDZUV(L,K)
          CMV = 1.+RCDZM*HV(L)*AVVI(L,K)
          EV  = 1./CMV
          CU2(L,K) = RCDZU*EV
          DV(L,K)  = (DV(L,K)-RCDZL*RCYY(L)*VHE(L)*HVI(L))*EV
          VVV(L,K) = EV
        else
          CU2(L,K) = 0.
          DV(L,K) = 0.
          VVV(L,K) = 0.
        endif
      enddo
    enddo

    ! ***  BACK SUBSTITUTION
    do K = KS-1,1,-1
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        DU(L,K) = DU(L,K)-CU1(L,K)*DU(L,K+1)
        DV(L,K) = DV(L,K)-CU2(L,K)*DV(L,K+1)
        UUU(L,K) = UUU(L,K)-CU1(L,K)*UUU(L,K+1)
        VVV(L,K) = VVV(L,K)-CU2(L,K)*VVV(L,K+1)
      enddo
    enddo

    ! *** SHERMAN-MORRISON BACK SUBSTITUTION
    do LP = 1,LLWETZ(KC,ND)
      L = LKWETZ(LP,KC,ND)
      AAU(L) = 0.0
      AAV(L) = 0.0
      BBU(L) = 1.0
      BBV(L) = 1.0
    enddo

    do K = 1,KS
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        CRU = CDZRU(L,K)*RCXX(L)*AVUI(L,K)
        AAU(L) = AAU(L) + CRU*DU(L,K)
        BBU(L) = BBU(L) + CRU*UUU(L,K)
        CRV = CDZRV(L,K)*RCYY(L)*AVVI(L,K)
        AAV(L) = AAV(L) + CRV*DV(L,K)
        BBV(L) = BBV(L) + CRV*VVV(L,K)
      enddo
    enddo
    do LP = 1,LLWETZ(KC,ND)
      L = LKWETZ(LP,KC,ND)
      AAU(L) = AAU(L)/BBU(L)
      AAV(L) = AAV(L)/BBV(L)
    enddo

    do K = 1,KS
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        DU(L,K) = SUB3D(L,K)*DZGU(L,K)*HU(L)*AVUI(L,K)*( DU(L,K) - AAU(L)*UUU(L,K) )
        DV(L,K) = SVB3D(L,K)*DZGV(L,K)*HV(L)*AVVI(L,K)*( DV(L,K) - AAV(L)*VVV(L,K) )
      enddo
    enddo

  enddo    ! ***  END OF DOMAIN LOOP
  !$OMP END DO

  ! ***************************************************************************
  ! *** CALCULATE U AND V USING THE INTERNAL SHEARS DU AND DV (DU/DV ARE IN M2/S)
  ! *** DUSUM+UHE = UHE, DVSUM+VHE = VHE
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET+1
    LL = MIN(LF+LDMWET-1,LAWET)

    ! *** ADJUST FLOWS BASED ON INTERNAL SHEARS

    ! *** INTERIM: COMPUTE SUM OF UNIT WIDTH Q PLUS UNIT FLOWS BY LAYER
    do K = 1,KS
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        UHE(L) = UHE(L) + CDZDU(L,K)*DU(L,K)
        VHE(L) = VHE(L) + CDZDV(L,K)*DV(L,K)
      enddo
    enddo

    ! *** LAYER SPECIFIC DEPTH AVERAGED FLOWS (M3/S)
    do LP = 1,LLWETZ(KC,ND)
      L = LKWETZ(LP,KC,ND)
      UHDYF(L,KC) = UHE(L)*SUB(L)
      VHDXF(L,KC) = VHE(L)*SVB(L)
    enddo
    do K = KS,1,-1
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        UHDYF(L,K) = SUB3D(L,K)*(UHDYF(L,K+1) - DU(L,K))
        VHDXF(L,K) = SVB3D(L,K)*(VHDXF(L,K+1) - DV(L,K))
      enddo
    enddo

    ! *** COMPUTE FLOWS IN M3/S FROM UNIT WIDTH Q
    do K = 1,KC
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        UHDYF(L,K) = UHDYF(L,K)*DYU(L)
        VHDXF(L,K) = VHDXF(L,K)*DXV(L)
      enddo
    enddo
  enddo    ! ***  END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP SINGLE

  ! *** BLOCKED LAYER FACE OPTION
  if( NBLOCKED > 0 )then
    do LP = 1,NBLOCKED
      L = LBLOCKED(LP)

      ! *** LAYER U BLOCKING
      if( BLDRAFTU(LP)+BLSILLU(LP) > 0.0 )then
        ! *** LAYER SPECIFIC DEPTH AVERAGED FLOWS (M3/S)
        UHDYF(L,:) = 0.0
        UHDYF(L,KTBU(LP)) = UHE(L)*SUB(L)

        do K = KTBU(LP)-1,KBBU(LP),-1
          UHDYF(L,K) = UHDYF(L,K+1) - DU(L,K)
        enddo

        ! *** COMPUTE FLOWS IN M3/S FROM UNIT WIDTH Q
        do K = KBBU(LP),KTBU(LP)
          UHDYF(L,K) = UHDYF(L,K)*DYU(L)
        enddo
      endif

      ! *** LAYER V BLOCKING
      if( BLDRAFTV(LP)+BLSILLV(LP) > 0.0 )then
        ! *** LAYER SPECIFIC DEPTH AVERAGED FLOWS (M3/S)
        VHDXF(L,:) = 0.0
        VHDXF(L,KTBV(LP)) = VHE(L)*SVB(L)

        do K = KTBV(LP)-1,KBBV(LP),-1
          VHDXF(L,K) = VHDXF(L,K+1) - DV(L,K)
        enddo

        ! *** COMPUTE FLOWS IN M3/S FROM UNIT WIDTH Q
        do K = KBBV(LP),KTBV(LP)
          VHDXF(L,K) = VHDXF(L,K)*DXV(L)
        enddo
      endif
    enddo
  endif
  !$OMP END SINGLE

  ! ***************************************************************************
  ! *** ADD ADJUSTMENT TO 3D HORIZONTAL TRANSPORT
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET+1
    LL = MIN(LF+LDMWET-1,LAWET)

    do LP = 1,LLWETZ(KC,ND)
      L = LKWETZ(LP,KC,ND)
      TVARE(L) = 0.0
      TVARN(L) = 0.0
    enddo

    ! *** SCALE TO DEPTH AVERAGE FLOW
    do K = 1,KC
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        TVARE(L) = TVARE(L) + UHDYF(L,K)*SGZU(L,K)
        TVARN(L) = TVARN(L) + VHDXF(L,K)*SGZV(L,K)
      enddo
    enddo

    ! *** COMPUTE DIFFERENCE FROM EXTERNAL SOLUTION
    do LP = 1,LLWETZ(KC,ND)
      L = LKWETZ(LP,KC,ND)
      TVARE(L) = TVARE(L) - UHDYE(L)
      TVARN(L) = TVARN(L) - VHDXE(L)
    enddo

    ! *** CORRECT INTERNAL SOLUTION LAYER SPECIFIC FLOWS
    do K = 1,KC
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        UHDYF(L,K) = UHDYF(L,K) - TVARE(L)*CERRU(L,K)
        UHDY(L,K)  = UHDYF(L,K)*SGZU(L,K)               ! *** LAYER SPECIFIC FLOWS
        VHDXF(L,K) = VHDXF(L,K) - TVARN(L)*CERRV(L,K)
        VHDX(L,K)  = VHDXF(L,K)*SGZV(L,K)               ! *** LAYER SPECIFIC FLOWS
      enddo
    enddo
  enddo    ! ***  END OF DOMAIN LOOP
  !$OMP END DO

  ! ***************************************************************************
  ! *** RESET VELOCITIES
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET+1
    LL = MIN(LF+LDMWET-1,LAWET)

    do K = 1,KC
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        U(L,K) = UHDYF(L,K)*HUI(L)
        V(L,K) = VHDXF(L,K)*HVI(L)
      enddo
    enddo
    do K = 1,KC
      do LP = 1,LLWETZ(K,ND)
        L = LKWETZ(LP,K,ND)
        U(L,K) = U(L,K)*DYIU(L)
        V(L,K) = V(L,K)*DXIV(L)
      enddo
    enddo
    
    ! *** Only needed for specific computational options
    if( IUPDATE_UVHE > 0 )then 
      do LP = 1,LLWETZ(KC,ND)
        L = LKWETZ(LP,KC,ND)
        UHE(L) = 0.0
        VHE(L) = 0.0
      enddo

      do K = 1,KC
        do LP = 1,LLWETZ(K,ND)
          L = LKWETZ(LP,K,ND)
          UHE(L) = UHE(L) + UHDYF(L,K)*SGZU(L,K)
          VHE(L) = VHE(L) + VHDXF(L,K)*SGZV(L,K)
        enddo
      enddo

      do LP = 1,LLWETZ(KC,ND)
        L = LKWETZ(LP,KC,ND)
        UHE(L) = UHE(L)*DYIU(L)
        VHE(L) = VHE(L)*DXIV(L)
      enddo
    endif
  enddo    ! ***  END OF DOMAIN LOOP
  !$OMP END DO

  !$OMP SINGLE
  call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)
  call Communicate_UVW1   ! *** Communicate UHDY and VHDX before W calculations
  
  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  !$OMP END SINGLE

  ! ***************************************************************************
  ! *** CALCULATE W
  if( IS2TL == 0 )then
    if( ISTL == 3 )then
    ! *** 
      !$OMP DO PRIVATE(ND,K,LP,L,LN,LE)
      do ND = 1,NDM
        do K = 1,KS
          do LP = 1,LLWETZ(K,ND)
            L = LKWETZ(LP,K,ND)
            LN = LNC(L)
            LE = LEC(L)
            W(L,K) = W(L,K-1) - W2(L,K) + W2(L,K-1)  &
                   - DXYIP(L)*( UHDY(LE,K)-UHDY(L,K) + UHDYEK(L,K) + UHDY2(LE,K)-UHDY2(L,K) + UHDY2EK(L,K)   &
                              + VHDX(LN,K)-VHDX(L,K) + VHDXEK(L,K) + VHDX2(LN,K)-VHDX2(L,K) + VHDX2EK(L,K) ) &
                   + 2.*( QSUM(L,K)-DZC(L,K)*QSUME(L) )*DXYIP(L)
          enddo
        enddo
      enddo    ! ***  END OF DOMAIN LOOP
      !$OMP END DO

    else
      ! *** 3TL CORRECTOR STEP
      !$OMP DO PRIVATE(ND,K,LP,L,LN,LE)
      do ND = 1,NDM
        do K = 1,KS
          do LP = 1,LLWETZ(K,ND)
            L = LKWETZ(LP,K,ND)
            LN = LNC(L)
            LE = LEC(L)
            W(L,K) = W(L,K-1) - W1(L,K) + W1(L,K-1)  &
                   - DXYIP(L)*( UHDY(LE,K)-UHDY(L,K) + UHDYEK(L,K) + UHDY1(LE,K)-UHDY1(L,K) + UHDY1EK(L,K)   &
                              + VHDX(LN,K)-VHDX(L,K) + VHDXEK(L,K) + VHDX1(LN,K)-VHDX1(L,K) + VHDX1EK(L,K) ) &
                   + 2.*( QSUM(L,K)-DZC(L,K)*QSUME(L) )*DXYIP(L)
          enddo
        enddo
      enddo    ! ***  END OF DOMAIN LOOP
      !$OMP END DO
    endif
  else
    ! *** TWO TIME LEVEL SOLUTION
    !$OMP DO PRIVATE(ND,K,LP,L,LN,LE)
    do ND = 1,NDM
      do K = 1,KS
        do LP = 1,LLWETZ(K,ND)
          L = LKWETZ(LP,K,ND)
          LN = LNC(L)
          LE = LEC(L)
          W(L,K) = W(L,K-1) - 0.5*DXYIP(L)  &
                   *( UHDY(LE,K) - UHDY(L,K) + UHDYEK(L,K) + UHDY1(LE,K)-UHDY1(L,K) + UHDY1EK(L,K)   &
                    + VHDX(LN,K) - VHDX(L,K) + VHDXEK(L,K) + VHDX1(LN,K)-VHDX1(L,K) + VHDX1EK(L,K) ) &
                    + ( QSUM(L,K) - DZC(L,K)*QSUME(L) )*DXYIP(L)
        enddo
      enddo
    enddo    ! ***  END OF DOMAIN LOOP
    !$OMP END DO  
  endif
  
  !$OMP END PARALLEL
  
  ! *** APPLY OPEN BOUNDARYS
  do LL = 1,NBCSOP
    L = LOBCS(LL)
    do K = 1,KS
      W(L,K) = 0.0
    enddo
  enddo

  ! ***************************************************************************
  ! *** JUMP TO POINT FOR KC = 1
30 continue

  ! ***************************************************************************
  ! *** CALCULATE U AND V ON OPEN BOUNDARIES
  do K = 1,KC
    do LL = 1,NCBS
      if( ISPBS(LL) /= 2 )then
        L = LCBS(LL)
        LN = LNC(L)
        LNN = LNC(LN)
        if( LN /= LC .and. K >= KSZ(L) )then
          VHDXF(LN,K) = VHDXF(LNN,K) - VHDXE(LNN) + VHDXE(LN)
          VHDX(LN,K)  = VHDXF(LN,K)*SGZV(LN,K)
          V(LN,K)     = VHDXF(LN,K)/(HV(LN)*DXV(LN))
          W(LN,K)     = 0.
        else
          W(LN,K)     = 0.
        endif
      endif
    enddo
  enddo

  do K = 1,KC
    do LL = 1,NCBW
      if( ISPBW(LL) /= 2 )then
        L = LCBW(LL)
        LE = LEC(L)
        L1P = LEC(LE)
        if( LE /= LC .and. K >= KSZ(L) )then
          UHDYF(LE,K) = UHDYF(L1P,K) - UHDYE(L1P) + UHDYE(LE)
          UHDY(LE,K)  = UHDYF(LE,K)*SGZU(LE,K)
          U(LE,K)     = UHDYF(LE,K)/(HU(LE)*DYU(LE))
          W(LE,K)     = 0.
        else
          W(LE,K)     = 0.
        endif
      endif
    enddo
  enddo

  do K = 1,KC
    do LL = 1,NCBE
      L = LCBE(LL)
      LW = LWC(L)
      if( ISPBE(LL) /= 2 .and. K >= KSZ(L) )then
        UHDYF(L,K) = UHDYF(LW,K) - UHDYE(LW) + UHDYE(L)
        UHDY(L,K)  = UHDYF(L,K)*SGZU(L,K)
        U(L,K)     = UHDYF(L,K)/(HU(L)*DYU(L))
        W(L,K)     = 0.
      else
        W(L,K)     = 0.
      endif
    enddo
  enddo

  do K = 1,KC
    do LL = 1,NCBN
      L = LCBN(LL)
      if( ISPBN(LL) /= 2 .and. K >= KSZ(L) )then
        LS = LSC(L)
        VHDXF(L,K) = VHDXF(LS,K) - VHDXE(LS) + VHDXE(L)
        VHDX(L,K)  = VHDXF(L,K)*SGZV(L,K)
        V(L,K)     = VHDXF(L,K)/(HV(L)*DXV(L))
        W(L,K)     = 0.
      else
        W(L,K)     = 0.
      endif
    enddo
  enddo

  ! ***************************************************************************
  ! *** MPI Communication
  call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)
  call Communicate_UVW3   ! *** Communicate U, V, and W

  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  ! ***************************************************************************
               
  !$OMP PARALLEL DEFAULT(SHARED)

  ! ***************************************************************************
  ! *** CALCULATE AVERAGE CELL FACE TRANSPORTS FOR SALT, TEMPERATURE AND
  ! *** SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2
  if( ISTL == 2 )then
    ! *** 3TL ISTL = 2 OR IS2TIM>0
    !$OMP DO PRIVATE(ND,K,LP,L,LN)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          UHDYF2(L,K) = 0.5*(UHDYF(L,K) + UHDYF1(L,K))
          VHDXF2(L,K) = 0.5*(VHDXF(L,K) + VHDXF1(L,K))
          UHDY2(L,K) = 0.5*(UHDY(L,K) + UHDY1(L,K))
          VHDX2(L,K) = 0.5*(VHDX(L,K) + VHDX1(L,K))
          U2(L,K) = 0.5*(U(L,K) + U1(L,K))
          V2(L,K) = 0.5*(V(L,K) + V1(L,K))
          W2(L,K) = 0.5*(W(L,K) + W1(L,K))
        enddo
      enddo
    enddo
    !$OMP END DO

  else

    ! *** 3TL ISTL = 3
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          UHDYF2(L,K) = 0.5*(UHDYF(L,K) + UHDYF2(L,K))
          VHDXF2(L,K) = 0.5*(VHDXF(L,K) + VHDXF2(L,K))
          UHDY2(L,K) = 0.5*(UHDY(L,K) + UHDY2(L,K))
          VHDX2(L,K) = 0.5*(VHDX(L,K) + VHDX2(L,K))
          U2(L,K) = 0.5*(U(L,K) + U2(L,K))
          V2(L,K) = 0.5*(V(L,K) + V2(L,K))
          W2(L,K) = 0.5*(W(L,K) + W2(L,K))
        enddo
      enddo
    enddo
    !$OMP END DO
  endif

  ! ***************************************************************************
  ! *** ADDITIONAL 3D CONTINUITY ADJUSTED ADDED BELOW
  if( KC > 1 )then
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K)
    do ND = 1,NDM
      do LP = 1,LLWETZ(KC,ND)
        L = LKWETZ(LP,KC,ND)
        TVARE(L) = 0.
        TVARN(L) = 0.
      enddo
      do K = 1,KC
        do LP = 1,LLWETZ(K,ND)
          L = LKWETZ(LP,K,ND)
          TVARE(L) = TVARE(L) + UHDY2(L,K)
          TVARN(L) = TVARN(L) + VHDX2(L,K)
        enddo
      enddo
    enddo
    !$OMP END DO

    ! *** HANDLE KSZ EAST = KC CASES
    !$OMP SINGLE
    do LP = 1,NKCE
      L = LKCE(LP)
      TVARE(L) = UHDY2(L,KC)
    enddo

    ! *** HANDLE KSZ NORTH = KC CASES
    do LP = 1,NKCN
      L = LKCN(LP)
      TVARN(L) = VHDX2(L,KC)
    enddo
    !$OMP END SINGLE

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LE,LW,LS,LN,K,IOBC,HPPTMP)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)

      if( ISTL == 3 )then
        ! *** 3TL AND ISTL = 3
        do LP = 1,LLWETZ(KC,ND)
          L = LKWETZ(LP,KC,ND)
          LE = LEC(L)
          LN = LNC(L)
          HPPTMP = H2P(L) + DELT*DXYIP(L)*( QSUME(L) - TVARE(LE)+TVARE(L)-TVARN(LN)+TVARN(L) )
          if( ISGWIE >= 1 ) HPPTMP = HPPTMP - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))

          ! *** NEGATIVE DEPTH CHECK
          if( HPPTMP <= 0. .and. HP(L) > 0. )then
            ! *** CHECK OPEN BC
            LW = 0
            do IOBC = 1,NBCSOP
              LW = LOBCS(IOBC)
              if( L == LW )EXIT
            enddo
            if( L == LW ) CYCLE

            LW = LWC(L)
            LS = LSC(L)
            if( HP(L) <= HDRY .and. ISDRY > 0 )then
              PRINT '(A,I8,I5,3F10.4,F12.4)',' WARNING! NEG DEPTH IN CONTINUITY CHECK: NITER,L,H2P,HP,HPNEW,TIMEDAY',NITER,map2global(L).LG,H2P(L),HP(L),HPPTMP,TIMEDAY
              if( IFILE == -1 )then
                IFILE = mpi_efdc_out_unit
                open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
              endif
              write(mpi_efdc_out_unit,'(A,I5,3F10.4,F12.4)')    ' WARNING!  NEG DEPTH IN CONTINUITY CHECK:'
              write(mpi_efdc_out_unit,'(A,I10,F12.4,I5,3F10.4)')'           NITER,TIMEDAY,L,H2P,HP,HPNEW',NITER,TIMEDAY,map2global(L).LG,H2P(L),HP(L),HPPTMP

              HPPTMP = HP(L)
              UHDYE(L) = 0.
              VHDXE(L) = 0.
              do K = 1,KC
                UHDYF(L,K) = 0.
                VHDXF(L,K) = 0.
                UHDY(L,K) = 0.
                VHDX(L,K) = 0.
                U(L,K) = 0.
                V(L,K) = 0.
                W(L,K) = 0.

                UHDYF(LE,K) = 0.
                VHDXF(LN,K) = 0.
                UHDY(LE,K) = 0.
                VHDX(LN,K) = 0.
                U(LE,K) = 0.
                V(LN,K) = 0.

                UHDYF2(L,K) = 0.
                VHDXF2(L,K) = 0.
                UHDY2(L,K) = 0.
                VHDX2(L,K) = 0.
                U2(L,K) = 0.
                V2(L,K) = 0.
                W2(L,K) = 0.
                UHDY2(LE,K) = 0.
                VHDX2(LN,K) = 0.
                U2(LE,K) = 0.
                V2(LN,K) = 0.
              enddo
            else
              PRINT '(A,I8,I5,3F10.4,F12.4)',' ERROR! NEG DEPTH IN CONTINUITY CHECK. N,L,H2P,HP,HPNEW,TIMEDAY', N, map2global(L).LG, H2P(L), HP(L), HPPTMP, TIMEDAY
                if( IFILE == -1 )then
                  IFILE = mpi_error_unit
                  open(mpi_error_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
                endif
                write(mpi_error_unit,'(A)')              ' ERROR!  NEGATIVE DEPTH IN CONTINUITY CHECK.'
                write(mpi_error_unit,'(A,I15,3I5,F15.5)')'         N L ISTL NCTBC TIMEDAY   ',N,L,ISTL,NCTBC,TIMEDAY
                write(mpi_error_unit,'(A,4F10.4)')       '         H1P H2P HP HPNEW         ',H1P(L),H2P(L),HP(L),HPPTMP
                write(mpi_error_unit,'(A,4F10.4)')       '         HP  WESN                 ',HP(LW),   HP(LE),   HP(LS),    HP(LN)
                write(mpi_error_unit,'(A,4F10.4)')       '         H1P WESN                 ',H1P(LW),  H1P(LE),  H1P(LS),   H1P(LN)
                write(mpi_error_unit,'(A,4F10.4)')       '         H2P WESN                 ',H2P(LW),  H2P(LE),  H2P(LS),   H2P(LN)
                write(mpi_error_unit,'(A,4F10.4)')       '         SUB/SVB                  ',SUB(L),   SUB(LE),   SVB(L),   SVB(LN)
                write(mpi_error_unit,'(A,4E12.4)')       '         UHDYE/VHDXE              ',UHDYE(L), UHDYE(LE), VHDXE(L), VHDXE(LN)
                write(mpi_error_unit,'(A,4E12.4)')       '         UHDYE1/VHDXE1            ',UHDY1E(L),UHDY1E(LE),VHDX1E(L),VHDX1E(LN)
                write(mpi_error_unit,'(A,4E12.4)')       '         UHDYE2/VHDXE2            ',UHDY2E(L),UHDY2E(LE),VHDX2E(L),VHDX2E(LN)
                write(mpi_error_unit,'(A,4E12.4,L5)')    '         WESN FLOWS SUM LAYERS    ',TVARE(L), TVARE(LE), TVARN(L), TVARN(LN), LMASKDRY(L)
                do K = KC,KSZ(L),-1
                  write(mpi_error_unit,'(A,I5,4E12.4)')  '         LAYER FLOWS/CURRENT      ',K,UHDY(L,K),UHDY(LE,K),VHDX(L,K),VHDX(LN,K)
                enddo
                do K = KC,KSZ(L),-1
                  write(mpi_error_unit,'(A,I5,4E12.4)')  '         LAYER FLOWS/OLD          ',K,UHDY1(L,K),UHDY1(LE,K),VHDX1(L,K),VHDX1(LN,K)
                enddo
                do K = KC,KSZ(L),-1
                  write(mpi_error_unit,'(A,I5,4E12.4)')  '         SUB3D/SVB3D              ',K,SUB3D(L,K),SUB3D(LE,K),SVB3D(L,K),SVB3D(LN,K)
                enddo
            endif
          endif

          HP(L)  = HPPTMP
          HPI(L) = 1./HP(L)
        enddo

      else
        ! *** 2TL / ISTL = 2
        do LP = 1,LLWETZ(KC,ND)
          L = LKWETZ(LP,KC,ND)
          LE = LEC(L)
          LN = LNC(L)
          HPPTMP = H1P(L) + DELT*DXYIP(L)*( QSUME(L) + TVARE(L)-TVARE(LE) + TVARN(L)-TVARN(LN) )
          if( ISGWIE >= 1 ) HPPTMP = HPPTMP - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))

          ! *** NEGATIVE DEPTH CHECK
          if( HPPTMP <= 0. .and. HP(L) > 0. )then
            ! *** CHECK OPEN BC
            LW = 0
            do IOBC = 1,NBCSOP
              LW = LOBCS(IOBC)
              if( L == LW )EXIT
            enddo
            if( L == LW ) CYCLE

            LW = LWC(L)
            LS = LSC(L)
            if( HP(L) <= HDRY .and. ISDRY > 0 )then
              PRINT '(A,I8,I5,3F10.4,F12.4)',' WARNING!  NEG DEPTH IN CONTINUITY CHECK. N,L,H1P,HP,HPNEW,TIMEDAY',N,map2global(L).LG,H1P(L),HP(L),HPPTMP,TIMEDAY
              if( IFILE == -1 )then
                IFILE = mpi_efdc_out_unit
                open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
              endif
              write(mpi_efdc_out_unit,'(A,I5,3F10.4,F12.4)')    ' WARNING!  NEG DEPTH IN CONTINUITY CHECK:'
              write(mpi_efdc_out_unit,'(A,I10,F12.4,I5,3F10.4)')'                N,TIMEDAY,L,H1P,HP,HPNEW',N,TIMEDAY,map2global(L).LG,H1P(L),HP(L),HPPTMP
              write(mpi_efdc_out_unit,'(A,4E12.4,L5)')          '                WESN FLOWS              ',TVARE(L),TVARE(LE),TVARN(L),TVARN(LN),LMASKDRY(L)

              HPPTMP = HP(L)
              UHDYE(L) = 0.
              VHDXE(L) = 0.
              do K = 1,KC
                LN = LNC(L)
                LE = LEC(L)

                UHDYF(L,K) = 0.
                VHDXF(L,K) = 0.
                UHDY(L,K) = 0.
                VHDX(L,K) = 0.
                U(L,K) = 0.
                V(L,K) = 0.
                W(L,K) = 0.

                UHDYF(LE,K) = 0.
                VHDXF(LN,K) = 0.
                UHDY(LE,K) = 0.
                VHDX(LN,K) = 0.
                U(LE,K) = 0.
                V(LN,K) = 0.

                UHDY2(L,K) = 0.
                VHDX2(L,K) = 0.
                U2(L,K) = 0.
                V2(L,K) = 0.
                W2(L,K) = 0.
                UHDY2(LE,K) = 0.
                VHDX2(LN,K) = 0.
                U2(LE,K) = 0.
                V2(LN,K) = 0.
              enddo
            else
              PRINT '(A,I8,I5,3F10.4,F12.4)',' ERROR! NEG DEPTH IN CONTINUITY CHECK. N,L,H1P,HP,HPNEW,TIMEDAY',N,map2global(L).LG,H1P(L),HP(L),HPPTMP,TIMEDAY
                if( IFILE == -1 )then
                  IFILE = mpi_error_unit
                  open(mpi_error_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
                endif
                write(mpi_error_unit,'(A)')              ' ERROR!  NEGATIVE DEPTH IN CONTINUITY CHECK.'
                write(mpi_error_unit,'(A,I15,3I5,F15.5)')'         N L ISTL NCTBC TIMEDAY   ',N,L,ISTL,NCTBC,TIMEDAY
                write(mpi_error_unit,'(A,4F10.4)')       '         H1P H2P HP HPNEW         ',H1P(L),H2P(L),HP(L),HPPTMP
                write(mpi_error_unit,'(A,4F10.4)')       '         HP  WESN                 ',HP(LW),   HP(LE),   HP(LS),    HP(LN)
                write(mpi_error_unit,'(A,4F10.4)')       '         H1P WESN                 ',H1P(LW),  H1P(LE),  H1P(LS),   H1P(LN)
                if( IS2TIM /= 0 )WRITE(mpi_error_unit,'(A,4F10.4)')       '         H2P WESN                 ',H2P(LW),  H2P(LE),  H2P(LS),   H2P(LN)
                write(mpi_error_unit,'(A,4F12.5)')       '         SUB/SVB                  ',SUB(L),   SUB(LE),   SVB(L),   SVB(LN)
                write(mpi_error_unit,'(A,4E12.4)')       '         UHDYE/VHDXE              ',UHDYE(L), UHDYE(LE), VHDXE(L), VHDXE(LN)
                write(mpi_error_unit,'(A,4E12.4)')       '         UHDYE1/VHDXE1            ',UHDY1E(L),UHDY1E(LE),VHDX1E(L),VHDX1E(LN)
                write(mpi_error_unit,'(A,4E12.4)')       '         UHDYE2/VHDXE2            ',UHDY2E(L),UHDY2E(LE),VHDX2E(L),VHDX2E(LN)
                write(mpi_error_unit,'(A,4E12.4,L5)')    '         WESN FLOWS SUM LAYERS    ',TVARE(L), TVARE(LE), TVARN(L), TVARN(LN), LMASKDRY(L)
                do K = KC,KSZ(L),-1
                  write(mpi_error_unit,'(A,I5,4E12.4)')  '         LAYER FLOWS/CURRENT      ',K,UHDY(L,K),UHDY(LE,K),VHDX(L,K),VHDX(LN,K)
                enddo
                do K = KC,KSZ(L),-1
                  write(mpi_error_unit,'(A,I5,4E12.4)')  '         LAYER FLOWS/OLD          ',K,UHDY1(L,K),UHDY1(LE,K),VHDX1(L,K),VHDX1(LN,K)
                enddo
                do K = KC,KSZ(L),-1
                  write(mpi_error_unit,'(A,I5,4E12.4)')  '         SUB3D/SVB3D              ',K,SUB3D(L,K),SUB3D(LE,K),SVB3D(L,K),SVB3D(LN,K)
                enddo
            endif
          endif

          HP(L)  = HPPTMP
          HPI(L) = 1./HP(L)
        enddo
      endif

    enddo
    !$OMP END DO
  endif

  ! ***************************************************************************
  ! *** Include Nondiverg Wave Stokes Drift In Mass Transport
  if( (ISWAVE == 2 .or. ISWAVE == 4) .and. ISWVSD >= 1 )then
    ! *** WV(L).HEIGHT  - WAVE HEIGHT (M)
    ! *** WV(L).DIR     - WAVE DIRECTION (RADIANS) COUNTER-CLOCKWISE (CELL-EAST AXIS,WAVE)
    ! *** WV(L).FREQ    - WAVE FREQENCY (SEC)
    ! *** WV(L).PERIOD  - WAVE PERIOD (SEC)
    ! *** WV(L).K       - WAVE NUMBER
    ! *** WV(L).LENGTH  - WAVE LENGTH (M)

    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,AMPL,SINHH,C1,C2,STOKESU,STOKESW,THETA,STOKESX,STOKESY)
    do ND = 1,NDM
      LF = (ND-1)*LDMWET+1
      LL = MIN(LF+LDMWET-1,LAWET)

      if( ISSSMMT == 2 )then
        do LP = LF,LL
          L = LWET(LP)
          if( LWVMASK(L) )then
            C1 = WV(L).KHP
            if( C1 < 20. .and. C1 > 1E-8 )then
              AMPL = 0.5*WV(L).HEIGHT
              AMPL = AMPL**2
              SINHH = SINH(C1)**2
              SINHH = AMPL*WV(L).FREQ*WV(L).K/SINHH
              SINH2(L) = 0.50*SINHH
              if( ABS(UHDYE(L)) > ABS(VHDXE(L)) )then
                SINH4(L) = 0.25*SINHH/WV(L).FREQ*WV(L).K * SUB(L)*(HP(LWC(L))-HP(L))*DXIU(L)
              else
                SINH4(L) = 0.25*SINHH/WV(L).FREQ*WV(L).K * SUB(L)*(HP(LSC(L))-HP(L))*DYIV(L)
              endif
            else
              SINH2(L) = 0.
              SINH4(L) = 0.
            endif

          endif
        enddo

        do K = 1,KC
          do LP = LF,LL
            L = LWET(LP)
            if( LWVMASK(L) )then
              ! *** COMPUTE STOKES DRIFT U/V
              C1 = 2.*WV(L).KHP*ZZ(L,K)
              STOKESU = SINH2(L)*COSH(C1)
              STOKESW = SINH4(L)*SINH(C1)
              THETA = WV(L).DIR
              STOKESX = SUB(L)*STOKESU*COS(THETA)
              STOKESY = SVB(L)*STOKESU*SIN(THETA)
              UHDY2(L,K) = UHDY2(L,K) + SUB(L)*STOKESX*DYU(L)*HU(L)
              VHDX2(L,K) = VHDX2(L,K) + SVB(L)*STOKESY*DXV(L)*HV(L)
              U2(L,K)    = U2(L,K) + STOKESX
              V2(L,K)    = V2(L,K) + STOKESY
              W2(L,K)    = W2(L,K) + STOKESW
            endif
          enddo
        enddo

      elseif( NTSMMT > 1 )then
        do K = 1,KC
          do LP = LF,LL
            L = LWET(LP)
            if( LWVMASK(L) )then
              UHDY2(L,K) = UHDY2(L,K) + UVPT(L,K)*DYU(L)
              VHDX2(L,K) = VHDX2(L,K) + VVPT(L,K)*DXV(L)
              U2(L,K)    = U2(L,K) + UVPT(L,K)*HUI(L)
              V2(L,K)    = V2(L,K) + VVPT(L,K)*HVI(L)
              W2(L,K)    = W2(L,K) + WVPT(L,K)
            endif
          enddo
        enddo
      endif
    enddo
    !$OMP END DO
  endif

  !$OMP SINGLE
  ! *** RESET OPEN BC CONCENTRATIONS
  do IOBC = 1,NBCSOP
    L = LOBCS(IOBC)
    HP(L)  = GI*P(L)-BELV(L)
    HPI(L) = 1./HP(L)
  enddo

  if( MDCHH >= 1 )then
    RLAMN = QCHERR
    RLAMO = 1.-RLAMN
    do NMD = 1,MDCHH
      LHOST = LMDCHH(NMD)
      LCHNU = LMDCHU(NMD)
      LCHNV = LMDCHV(NMD)
      if( MDCHTYP(NMD) == 1 )then
        TMPVAL = DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUN(NMD))
        HP(LHOST) = HP(LHOST)+TMPVAL*DXYIP(LHOST)
        HP(LCHNU) = HP(LCHNU)-TMPVAL*DXYIP(LCHNU)
        HPI(LHOST) = 1./HP(LHOST)
        HPI(LCHNU) = 1./HP(LCHNU)
      endif
      if( MDCHTYP(NMD) == 2 )then
        TMPVAL = DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVN(NMD))
        HP(LHOST) = HP(LHOST)+TMPVAL*DXYIP(LHOST)
        HP(LCHNV) = HP(LCHNV)-TMPVAL*DXYIP(LCHNV)
        HPI(LHOST) = 1./HP(LHOST)
        HPI(LCHNV) = 1./HP(LCHNV)
      endif
    enddo
  endif
  
  call MPI_barrier(MPI_Comm_World, ierr)
  TTDS = DSTIME(0)
  call Communicate_1D2(HP, HPI)   ! *** Communicate HP and HPI after continuity checks
  !Call communicate_ghost_cells(HP, 'HP')   ! *** Communicate HP after continuity checks

  TMPITMP = TMPITMP + (DSTIME(0)-TTDS)
  !$OMP END SINGLE

  ! *** COMPUTE FACE DEPTHS AND LAYER THICKNESSES
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET+1
    LL = MIN(LF+LDMWET-1,LAWET)
    do K = 1,KC
      do LP = 1,LLWET(K,ND)
        L = LKWET(LP,K,ND)
        HPK(L,K) = HP(L)*DZC(L,K)
        HPKI(L,K) = 1./HPK(L,K)
      enddo
    enddo

    if( IGRIDV > 0 )then
      do LP = LF,LL
        L = LWET(LP)
        ! *** SET DIRECTIONAL DEPTHS
        HPW(L) = HP(L) + BELV(L) - BELVW(L)
        HPE(L) = HP(L) + BELV(L) - BELVE(L)
        HPS(L) = HP(L) + BELV(L) - BELVS(L)
        HPN(L) = HP(L) + BELV(L) - BELVN(L)
      enddo
    endif
  enddo
  !$OMP END DO

  !$OMP END PARALLEL

  ! ***************************************************************************
  ! *** ACCUMULTATE MAX COURANT NUMBERS
  if( ISINWV == 1 )then
    do K = 1,KC
      do L = 2,LA
        CFLUUUT     = DELT*ABS(DXIU(L)*U(L,K))
        CFLUUU(L,K) = MAX(CFLUUUT,CFLUUU(L,K))
        CFLVVVT     = DELT*ABS(DYIV(L)*V(L,K))
        CFLVVV(L,K) = MAX(CFLVVVT,CFLVVV(L,K))
        CFLWWWT     = DELT*ABS(HPI(L)*DZIG(L,K)*W(L,K))
        CFLWWW(L,K) = MAX(CFLWWWT,CFLWWW(L,K))
        CFLCACT     = DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))
        CFLCAC(L,K) = MAX(CFLCACT,CFLCAC(L,K))
      enddo
    enddo
  endif

  ! ***************************************************************************
  ! ***CALCULATE NONHYDROSTATIC PRESSURE
  if( KC > 1 .and. ISPNHYDS >= 1 ) CALL CALPNHS

  ! ***************************************************************************
  ! *** WRITE TO DIAGNOSTIC FILE CFL.OUT WITH DIAGNOSTICS OF MAXIMUM TIME STEP
  ! *** SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2
  if( ISCFL >= 1 .and. ISTL == 3 .and. DEBUG )then

    if( IFILE == -1 )then
      IFILE = mpi_efdc_out_unit
      open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
    endif
    if( NITER == 1 .and. (ISCFL >= 1 .or. ISCFLM >= 1) )then
      write(mpi_efdc_out_unit,'(a)'), ' *** Writing CFL.OUT'

      if( ISCFLM >= 1 )then
        do L = 1,LC
          ICFLMP(L) = 0
        enddo
      endif
    endif
    
    DTCFL = 1.E+18
    K = 1
    do L = 2,LA
      LE = LEC(L)
      LN = LNC(L)
      UWTMP = ABS(DXIU(L) *U2(L,K))
      UETMP = ABS(DXIU(LE)*U2(LE,K))
      VSTMP = ABS(DYIV(L) *V2(L,K))
      VNTMP = ABS(DYIV(LN)*U2(LN,K))
      WBTMP = 0.
      WTTMP = ABS(HPKI(L,K)*W2(L,K))
      DTMAXI = MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  +1.0E-12
      DTMAXX = 0.5/DTMAXI
      if( DTMAXX < DTCFL )then
        DTCFL = DTMAXX
        ICFL = IL(L)
        JCFL = JL(L)
        KCFL = K
      endif
    enddo
    if( KC > 1 )then
      K = KC
      do L = 2,LA
        LN = LNC(L)
        UWTMP = ABS(DXIU(L  )*U2(L  ,K))
        UETMP = ABS(DXIU(LE)*U2(LE,K))
        VSTMP = ABS(DYIV(L  )*V2(L  ,K))
        VNTMP = ABS(DYIV(LN )*U2(LN ,K))
        WTTMP = 0.
        WBTMP = ABS(HPKI(L,K)*W2(L,K-1))
        DTMAXI = MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  +1.0E-12
        DTMAXX = 0.5/DTMAXI
        if( DTMAXX < DTCFL )then
          DTCFL = DTMAXX
          ICFL = IL(L)
          JCFL = JL(L)
          KCFL = K
        endif
      enddo
    endif
    if( KC > 2 )then
      do K = 2,KS
        do L = 2,LA
          LN = LNC(L)
          UWTMP = ABS(DXIU(L) *U2(L,K))
          UETMP = ABS(DXIU(LE)*U2(LE,K))
          VSTMP = ABS(DYIV(L) *V2(L,K))
          VNTMP = ABS(DYIV(LN)*U2(LN,K))
          WBTMP = ABS(HPKI(L,K)*W2(L,K-1))
          WTTMP = ABS(HPKI(L,K)*W2(L,K))
          DTMAXI = MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  +1.0E-12
          DTMAXX = 0.5/DTMAXI
          if( DTMAXX < DTCFL )then
            DTCFL = DTMAXX
            ICFL = IL(L)
            JCFL = JL(L)
            KCFL = K
          endif
        enddo
      enddo
    endif
    IVAL = MOD(N,ISCFL)
    IDTCFL = NINT(DTCFL)

    if( ISCFL == 1 ) WRITE(mpi_efdc_out_unit,1212) DTCFL,N,ICFL,JCFL,KCFL
    if( ISCFL >= 2 .and. IVAL == 0 ) WRITE(mpi_efdc_out_unit,1213) IDTCFL
    if( ISCFLM >= 1 )then
      LTMP = LIJ(ICFL,JCFL)
      ICFLMP(LTMP) = ICFLMP(LTMP)+1
    endif
    if( ISCFLM >= 1 .and. N == NTS )then
      TMPVALN = 1./FLOAT(NTS)
      do L = 2,LA
        TMPVAL = TMPVALN*FLOAT(ICFLMP(L))
        write(mpi_efdc_out_unit,1214) IL(L),JL(L),ICFLMP(L),TMPVAL
      enddo
    endif
1212  FORMAT(' MAX TIME STEP  = ',F10.2,' SEC FOR N,I,J,K  = ',I8,3I5)
1213  FORMAT(I4)
1214  FORMAT(2I5,I12,F10.2)
  endif

  if( IFILE == mpi_efdc_out_unit ) close(mpi_efdc_out_unit)
  if( IFILE == mpi_error_unit ) close(mpi_error_unit)
  ! ***************************************************************************

  return

END

