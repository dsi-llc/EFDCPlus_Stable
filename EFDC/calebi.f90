! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALEBI

  ! ***  CALEBI CALCULATES THE EXTERNAL BUOYANCY INTEGRALS  
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2016-02       PAUL M. CRAIG     UPDATED SIGMA-Z (SGZ) FOR EE8.0 
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP

  use GLOBAL  
  use Allocate_Initialize
  
  implicit none

  integer :: K, L, LP, ND, LF, LL, LW, LS

  real :: DZCBK, DBK, BEC, BI1C, BI2C, SUMT
  
  real,save,allocatable,dimension(:,:) :: DZCB
  real,save,allocatable,dimension(:,:) :: BK
  real,save,allocatable,dimension(:,:) :: DZCC
  
  logical,save,allocatable,dimension(:) :: KSZFLAT
  
  if( .not. allocated(DZCB) )then
    call AllocateDSI( DZCB,     KCM,  NDM,  0.0 )
    call AllocateDSI( BK,       KCM,  NDM,  0.0 )
    call AllocateDSI( DZCC,     KCM,  LCM,  0.0 )
    call AllocateDSI( KSZFLAT,  LCM,  .FALSE. ) 

    ! *** SET UP K ORDERED ARRAYS
    do L = 1,LA
      do K = 1,KC
        DZCC(K,L) = DZC(L,K)
      enddo
    enddo
    do L = 1,LA
      do K = 0,KC
        ZZC(K,L)  = ZZ(L,K)
      enddo
    enddo

    ! ***
    if( IGRIDV == 1 )then
      do L = 2,LA
        if( ABS(SUBO(L)*KSZ(L)-SUBO(L)*KSZ(LWC(L))) < 1. .and. ABS(SUBO(LEC(L))*KSZ(L)-SUBO(LEC(L))*KSZ(LEC(L))) < 1. .and. &    ! *** Addresses real4/8 precision
            ABS(SVBO(L)*KSZ(L)-SVBO(L)*KSZ(LSC(L))) < 1. .and. ABS(SVBO(LNC(L))*KSZ(L)-SVBO(LNC(L))*KSZ(LNC(L))) < 1. )then
          KSZFLAT(L) = .TRUE.
        endif
      enddo
    endif
  endif

  if( IGRIDV > 0 )then               
    ! *** INTERPOLATE B ARRAYS MIDPOINT FOR EACH FACE
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, K)
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)
      
      do LP = LF,LL
        L = LWET(LP)
        if( IGRIDV > 1 )then
          CALL INTERPB(L)
        else
          ! *** Using B at cells layer center in SGZ-IVGRID 1 for consistency
          ! *** in the implementation of external density gradient and buoyancy shears - DKT
          do K = KSZ(L),KC
            BW(L,K) = B(L,K)
            BE(L,K) = B(L,K)
            BS(L,K) = B(L,K)
            BN(L,K) = B(L,K)
          enddo
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif
  
  
  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,K,DBK,DZCBK,BEC,BI1C,BI2C)            &
  !$OMP                           SHARED(NDM,LDMWET,LAWET,LWET,LEC,LNC,LSC,LWC,KSZ,KC,IGRIDV) &
  !$OMP                           SHARED(B,SVB,SUB,DZCC,ZZC,DZCB,BK,BW,BE,BS,BN,KSZFLAT)      &
  !$OMP                           SHARED(BI1W,BI2W,BEW,KSZW,ZZW,SGZKW)    &
  !$OMP                           SHARED(BI1E,BI2E,BEE,KSZE,ZZE,SGZKE)    &
  !$OMP                           SHARED(BI1S,BI2S,BES,KSZS,ZZS,SGZKS)    &
  !$OMP                           SHARED(BI1N,BI2N,BEN,KSZN,ZZN,SGZKN)
  do ND = 1,NDM  
    LF = (ND-1)*LDMWET+1  
    LL = min(LF+LDMWET-1,LAWET)
      
    do LP = LF,LL
      L = LWET(LP)
      
      ! *** CENTROID INTEGRAL 
      do K = KSZ(L),KC
        DZCB(K,ND) = DZCC(K,L)*B(L,K)
      enddo

      DBK = 0.  
      do K = KC,KSZ(L),-1
        DBK = DBK + DZCB(K,ND)       
        BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
      enddo

      BEC  = 0.
      BI1C = 0.
      BI2C = 0.
      do K = KC,KSZ(L),-1
        BEC   = BEC  + DZCB(K,ND)
        DZCBK = DZCC(K,L)*BK(K,ND)
        BI1C  = BI1C + DZCBK  
        BI2C  = BI2C + DZCBK + ZZC(K,L)*DZCB(K,ND) 
      enddo

      if( IGRIDV == 0 )then
        ! *** Standard SIGMA code
        BI1W(L) = BI1C
        BI2W(L) = BI2C
        BEW(L)  = BEC
        BI1S(L) = BI1C
        BI2S(L) = BI2C
        BES(L)  = BEC
      else   ! *** SGZ FACE INTEGRALS
        if( KSZFLAT(L) )then
          ! *** Initialize to standard SIGMA code
          BI1W(L) = BI1C
          BI2W(L) = BI2C
          BEW(L)  = BEC

          BI1E(L) = BI1C
          BI2E(L) = BI2C
          BEE(L)  = BEC

          BI1S(L) = BI1C
          BI2S(L) = BI2C
          BES(L)  = BEC

          BI1N(L) = BI1C
          BI2N(L) = BI2C
          BEN(L)  = BEC
          
          CYCLE
        endif
        ! *** FACE INTEGRALS: WEST
        if( SUB(L) > 0. )then
          BI1W(L) = 0.
          BI2W(L) = 0.
          BEW(L)  = 0.
          do K = KSZW(L),KC
            DZCB(K,ND) = SGZKW(K,L)*BW(L,K)
          enddo

          DBK = 0.  
          do K = KC,KSZW(L),-1
            DBK = DBK+DZCB(K,ND)       
            BK(K,ND) = DBK-0.5*DZCB(K,ND) 
          enddo

          do K = KC,KSZW(L),-1
            BEW(L)  = BEW(L)  + DZCB(K,ND)
            DZCBK   = SGZKW(K,L)*BK(K,ND)
            BI1W(L) = BI1W(L) + DZCBK  
            BI2W(L) = BI2W(L) + DZCBK + ZZW(K,L)*DZCB(K,ND) 
          enddo
        else
          BI1W(L) = BI1C
          BI2W(L) = BI2C
          BEW(L)  = BEC
        endif

        ! *** FACE INTEGRALS: EAST
        if( SUB(LEC(L)) > 0. )then
          BI1E(L) = 0.
          BI2E(L) = 0.
          BEE(L)  = 0.
          do K = KSZE(L),KC
            DZCB(K,ND) = SGZKE(K,L)*BE(L,K)
          enddo

          DBK = 0.  
          do K = KC,KSZE(L),-1
            DBK = DBK+DZCB(K,ND)       
            BK(K,ND) = DBK-0.5*DZCB(K,ND) 
          enddo

          do K = KC,KSZE(L),-1
            BEE(L)  = BEE(L)  + DZCB(K,ND)
            DZCBK   = SGZKE(K,L)*BK(K,ND)
            BI1E(L) = BI1E(L) + DZCBK  
            BI2E(L) = BI2E(L) + DZCBK + ZZE(K,L)*DZCB(K,ND) 
          enddo
        else
          BI1E(L) = BI1C
          BI2E(L) = BI2C
          BEE(L)  = BEC
        endif

        ! *** FACE INTEGRALS: SOUTH
        if( SVB(L) > 0. )then
          BI1S(L) = 0.
          BI2S(L) = 0.
          BES(L)  = 0.
          do K = KSZS(L),KC
            DZCB(K,ND) = SGZKS(K,L)*BS(L,K)
          enddo

          DBK = 0.  
          do K = KC,KSZS(L),-1
            DBK = DBK+DZCB(K,ND)       
            BK(K,ND) = DBK-0.5*DZCB(K,ND) 
          enddo

          do K = KC,KSZS(L),-1
            BES(L)  = BES(L)  + DZCB(K,ND)
            DZCBK   = SGZKS(K,L)*BK(K,ND)
            BI1S(L) = BI1S(L) + DZCBK  
            BI2S(L) = BI2S(L) + DZCBK + ZZS(K,L)*DZCB(K,ND) 
          enddo
        else
          BI1S(L) = BI1C
          BI2S(L) = BI2C
          BES(L)  = BEC
        endif

        ! *** FACE INTEGRALS: NORTH
        if( SVB(LNC(L)) > 0. )then
          BI1N(L) = 0.
          BI2N(L) = 0.
          BEN(L)  = 0.
          do K = KSZN(L),KC
            DZCB(K,ND) = SGZKN(K,L)*BN(L,K)
          enddo

          DBK = 0.  
          do K = KC,KSZN(L),-1
            DBK = DBK+DZCB(K,ND)       
            BK(K,ND) = DBK-0.5*DZCB(K,ND) 
          enddo

          do K = KC,KSZN(L),-1
            BEN(L)  = BEN(L)  + DZCB(K,ND)
            DZCBK   = SGZKN(K,L)*BK(K,ND)
            BI1N(L) = BI1N(L) + DZCBK  
            BI2N(L) = BI2N(L) + DZCBK + ZZN(K,L)*DZCB(K,ND) 
          enddo
        else
          BI1N(L) = BI1C
          BI2N(L) = BI2C
          BEN(L)  = BEC
        endif
      endif
      
    enddo
  enddo
  !$OMP END PARALLEL DO

  ! *** ADJUST INTEGRALS FOR FACE BLOCKING
  if( NBLOCKED > 0 .and. N > 1 )then
    ND = 1
    do LP = 1,NBLOCKED
      L = LBLOCKED(LP)
      if( KSZ(L) == KC ) CYCLE
      
      ! *** LAYER U BLOCKING
      if( SUB(L) > 0. .and. (BLDRAFTU(LP)+BLSILLU(LP)) > 0. )then
        ! *** WEST FACE INTEGRAL OF L

        ! *** METRICS
        SUMT = 0.
        DZCC(:,L) = 0.
        do K = KBBU(LP),KC
          DZCC(K,L) = DZC(L,K)
          SUMT = SUMT + DZCC(K,L)
        enddo
        DZCC(:,L) = DZCC(:,L)/SUMT

        ZW(L,KBBU(LP)-1) = 0.
        do K = KBBU(LP),KTBU(LP)
          ZW(L,K)  = ZW(L,K-1) + SGZU(L,K)          ! *** TOP OF LAYER Z FOR OPENING
        enddo

        SUMT = 0.
        ZZW(:,L) = 0.
        do K = KBBU(LP),KC
          SUMT = SUMT + DZCC(K,L)
          ZZW(K,L) = SUMT - 0.5*DZCC(K,L)           ! *** MID LAYER Z TO SURFACE
        enddo
        
        ! *** BUOYANCY INTEGRAL
        do K = KBBU(LP),KC  
          DZCB(K,ND) = DZCC(K,L)*B(L,K)
        enddo

        DBK = 0.  
        do K = KC,KBBU(LP),-1
          DBK = DBK + DZCB(K,ND)       
          BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
        enddo
            
        BEC = 0.
        BI1C = 0.
        BI2C = 0.
        do K = KC,KBBU(LP),-1
          BEC   = BEC  + DZCB(K,ND)
          DZCBK = DZCC(K,L)*BK(K,ND)
          BI1C  = BI1C + DZCBK  
          BI2C  = BI2C + DZCBK + ZZW(K,L)*DZCB(K,ND) 
        enddo
        BI1W(L) = BI1C
        BI2W(L) = BI2C
        BEW(L)  = BEC
      
        ! *** EAST FACE INTEGRAL OF LW
        if( SUBO(L) > 0.0 )then
          LW = LWC(L)

          ZE(LW,KBBU(LP)-1) = 0.
          do K = KBBU(LP),KTBU(LP)
            ZE(LW,K)  = ZE(LW,K-1) + SGZU(L,K)         ! *** TOP OF LAYER Z FOR OPENING
          enddo

          SUMT = 0.
          ZZE(:,LW) = 0.
          do K = KBBU(LP),KC
            SUMT = SUMT + DZCC(K,L)
            ZZE(K,LW) = SUMT - 0.5*DZCC(K,L)           ! *** MID LAYER Z TO SURFACE
          enddo

          do K = KBBU(LP),KC
            DZCB(K,ND) = DZCC(K,L)*B(LW,K)
          enddo

          DBK = 0.  
          do K = KC,KBBU(LP),-1
            DBK = DBK + DZCB(K,ND)       
            BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
          enddo

          BEC = 0.
          BI1C = 0.
          BI2C = 0.
          do K = KC,KBBU(LP),-1
            BEC   = BEC  + DZCB(K,ND)
            DZCBK = DZCC(K,L)*BK(K,ND)
            BI1C  = BI1C + DZCBK  
            BI2C  = BI2C + DZCBK + ZZE(K,LW)*DZCB(K,ND) 
          enddo
          BI1E(LW) = BI1C
          BI2E(LW) = BI2C
          BEE(LW)  = BEC
        endif
      endif

      ! *** LAYER V BLOCKING
      if( SVB(L) > 0. .and. (BLDRAFTV(LP)+BLSILLV(LP)) > 0. )then
        ! *** V FACE INTEGRAL 

        ! *** METRICS
        SUMT = 0.
        DZCC(:,L) = 0.
        do K = KBBV(LP),KC
          DZCC(K,L) = DZC(L,K)
          SUMT = SUMT + DZCC(K,L)
        enddo
        DZCC(:,L) = DZCC(:,L)/SUMT

        ZS(L,:) = 0.
        do K = KBBV(LP),KTBV(LP)
          ZS(L,K)  = ZS(L,K-1) + SGZV(L,K)          ! *** TOP OF LAYER Z FOR OPENING
        enddo
        
        SUMT = 0.
        ZZS(:,L) = 0.
        do K = KBBV(LP),KC
          SUMT = SUMT + DZCC(K,L)
          ZZS(K,L) = SUMT - 0.5*DZCC(K,L)           ! *** MID LAYER Z TO SURFACE
        enddo
        
        ! *** BUOYANCY INTEGRAL
        do K = KBBV(LP),KC
          DZCB(K,ND) = DZCC(K,L)*B(L,K)
        enddo

        DBK = 0.  
        do K = KC,KBBV(LP),-1
          DBK = DBK + DZCB(K,ND)       
          BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
        enddo

        BEC = 0.
        BI1C = 0.
        BI2C = 0.
        do K = KC,KBBV(LP),-1
          BEC   = BEC  + DZCB(K,ND)
          DZCBK = DZCC(K,L)*BK(K,ND)
          BI1C  = BI1C + DZCBK  
          BI2C  = BI2C + DZCBK + ZZS(K,L)*DZCB(K,ND) 
        enddo
        BI1S(L) = BI1C
        BI2S(L) = BI2C
        BES(L)  = BEC

        ! *** NORTH FACE INTEGRAL OF LS
        if( SVBO(L) > 0.0 )then
          LS = LSC(L)
          ZN(LS,KBBV(LP)-1) = 0.
          do K = KBBV(LP),KTBV(LP)
            ZN(LS,K)  = ZN(LS,K-1) +     SGZV(L,K)    ! *** TOP OF LAYER Z
            ZZN(K,LS) = ZN(LS,K)   - 0.5*SGZV(L,K)    ! *** MID LAYER Z
          enddo
          
          do K = KBBV(LP),KTBV(LP)
            DZCB(K,ND) = SGZV(L,K)*B(LS,K)
          enddo

          DBK = 0.  
          do K = KTBV(LP),KBBV(LP),-1
            DBK = DBK + DZCB(K,ND)       
            BK(K,ND) = DBK-0.5*DZCB(K,ND) 
          enddo

          BEC = 0.
          BI1C = 0.
          BI2C = 0.
          do K = KTBV(LP),KBBV(LP),-1
            BEC   = BEC  + DZCB(K,ND)
            DZCBK = SGZV(L,K)*BK(K,ND)
            BI1C  = BI1C + DZCBK  
            BI2C  = BI2C + DZCBK + ZZN(K,LS)*DZCB(K,ND) 
          enddo
          BI1N(L) = BI1C
          BI2N(L) = BI2C
          BEN(L)  = BEC
        endif
      endif
    enddo
  endif
  
  return 

END  

SUBROUTINE INTERPB(L)

  ! *** RETURNS THE B ARRAY INTERPOLATED ONTO EACH FACE AT THE LAYER MIDPOINT
  use GLOBAL
  
  integer, intent(IN) :: L

  integer :: K,KF,LW,LE,LS,LN,NMAX,KK
  real :: ZF, XDAT(1:KC),YDAT(1:KC)
  
  ! *** CELL CENTROID VALUES APPLY FOR ALL LAYERS BUT 
  NMAX = KC-KSZ(L)+1
  KK = 0
  do K = KSZ(L),KC
    KK = KK+1
    BW(L,K) = B(L,K)
    BE(L,K) = B(L,K)
    BS(L,K) = B(L,K)
    BN(L,K) = B(L,K)
    XDAT(KK)= BELV(L)  + HP(L)*ZZC(K,L)
    YDAT(KK)= B(L,K)
  enddo

  LW = LWC(L)
  LE = LEC(L)
  LS = LSC(L)
  LN = LNC(L)
  
  ! *** INTERPOLATE VARIN FOR THE WEST FACE
  if( KSZ(LW) > KSZ(L) .OR. ( KSZ(LW) == KSZ(L) .AND. BELV(L) < BELV(LW) ) )then    
    do KF = KSZW(L),KC
      ZF = BELVE(LW) + HPE(LW)*ZZE(KF,LW)
      call INTERPOL(NMAX,XDAT,YDAT,ZF,BW(L,KF))
    enddo
  endif

  ! *** INTERPOLATE VARIN FOR THE EAST FACE
  if( KSZ(LE) > KSZ(L) .OR. ( KSZ(LE) == KSZ(L) .AND. BELV(L) < BELV(LE) ) )then
    do KF = KSZE(L),KC
      ZF = BELVW(LE) + HPW(LE)*ZZW(KF,LE)
      call INTERPOL(NMAX,XDAT,YDAT,ZF,BE(L,KF))
    enddo
  endif

  ! *** INTERPOLATE VARIN FOR THE SOUTH FACE
  if( KSZ(LS) > KSZ(L) .OR. ( KSZ(LS) == KSZ(L) .AND. BELV(L) < BELV(LS) ) )then   
    do KF = KSZS(L),KC
      ZF = BELVN(LS) + HPN(LS)*ZZN(KF,LS)
      call INTERPOL(NMAX,XDAT,YDAT,ZF,BS(L,KF))
    enddo
  endif

  ! *** INTERPOLATE VARIN FOR THE NORTH FACE
  if( KSZ(LN) > KSZ(L) .OR. ( KSZ(LN) == KSZ(L) .AND. BELV(L) < BELV(LN) ) )then
    do KF = KSZN(L),KC
      ZF = BELVS(LN) + HPS(LN)*ZZS(KF,LN)
      call INTERPOL(NMAX,XDAT,YDAT,ZF,BN(L,KF))
    enddo
  endif
  
  contains
  
  SUBROUTINE INTERPOL(NMAX,XDAT,YDAT,XVAL,YVAL)
    ! *** INTERPOLATION FOR YVAL OF XVAL
    integer, intent(IN ) :: NMAX
    real,    intent(IN ) :: XVAL,XDAT(NMAX),YDAT(NMAX)
    real,    intent(OUT) :: YVAL
    integer       :: N

    if( XVAL < XDAT(1) )then 
      YVAL = YDAT(1)
    elseif(  XVAL > XDAT(NMAX) )then 
      YVAL = YDAT(NMAX)
    else 
      do N = 1,NMAX-1
        if( XVAL >= XDAT(N) .and. XVAL <= XDAT(N+1) )then 
          YVAL  =  (YDAT(N+1)-YDAT(N))*(XVAL-XDAT(N)) /(XDAT(N+1)-XDAT(N))+YDAT(N)
          exit
        endif 
      enddo
    endif
  END SUBROUTINE
  
 END SUBROUTINE
