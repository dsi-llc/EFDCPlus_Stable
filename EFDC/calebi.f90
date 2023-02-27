! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALEBI

  ! **  CALEBI CALCULATES THE EXTERNAL BUOYANCY INTEGRALS  
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2016-02       PAUL M. CRAIG     UPDATED SIGMA-Z (SGZ) FOR EE8.0 
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP

  USE GLOBAL  
  IMPLICIT NONE

  INTEGER :: K,L,LP,ND,LF,LL,LW,LS

  REAL :: DZCBK,DBK,BEC,BI1C,BI2C,SUMT
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DZCB
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: BK
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DZCC
  
  LOGICAL,SAVE,ALLOCATABLE,DIMENSION(:) :: KSZFLAT
  
  IF( .NOT. ALLOCATED(DZCB) )THEN
    ALLOCATE(DZCB(KCM,NDM))
    ALLOCATE(BK(KCM,NDM))
    ALLOCATE(DZCC(KCM,LCM))
    ALLOCATE(KSZFLAT(LCM))

    DZCB = 0.0
    BK   = 0.0

    ! *** SET UP K ORDERED ARRAYS
    DO L=1,LA
      DO K=1,KC
        DZCC(K,L) = DZC(L,K)
      ENDDO
    ENDDO
    DO L=1,LA
      DO K=0,KC
        ZZC(K,L)  = ZZ(L,K)
      ENDDO
    ENDDO

    ! ***
    KSZFLAT = .FALSE.
    IF( IGRIDV == 1 )THEN
      DO L=2,LA
        IF( ABS(SUBO(L)*KSZ(L)-SUBO(L)*KSZ(LWC(L))) < 1. .AND. ABS(SUBO(LEC(L))*KSZ(L)-SUBO(LEC(L))*KSZ(LEC(L))) < 1. .AND. &    ! *** Addresses real4/8 precision
            ABS(SVBO(L)*KSZ(L)-SVBO(L)*KSZ(LSC(L))) < 1. .AND. ABS(SVBO(LNC(L))*KSZ(L)-SVBO(LNC(L))*KSZ(LNC(L))) < 1. )THEN
          KSZFLAT(L) = .TRUE.
        ENDIF
      ENDDO
    ENDIF
  ENDIF


  IF( IGRIDV > 1 )THEN
    ! *** INTERPOLATE B ARRAYS MIDPOINT FOR EACH FACE
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      DO LP=LF,LL
        L=LWET(LP)
        IF( .NOT. KSZFLAT(L) ) CALL INTERPB(L)
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
  
  
  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,K,DBK,DZCBK,BEC,BI1C,BI2C)            &
  !$OMP                           SHARED(NDM,LDMWET,LAWET,LWET,LEC,LNC,LSC,LWC,KSZ,KC,IGRIDV) &
  !$OMP                           SHARED(B,SVB,SUB,DZCC,ZZC,DZCB,BK,BW,BE,BS,BN,KSZFLAT)      &
  !$OMP                           SHARED(BI1W,BI2W,BEW,KSZW,ZZW,SGZKW)    &
  !$OMP                           SHARED(BI1E,BI2E,BEE,KSZE,ZZE,SGZKE)    &
  !$OMP                           SHARED(BI1S,BI2S,BES,KSZS,ZZS,SGZKS)    &
  !$OMP                           SHARED(BI1N,BI2N,BEN,KSZN,ZZN,SGZKN)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
      
    DO LP=LF,LL
      L=LWET(LP)
      
      ! *** CENTROID INTEGRAL 
      DO K=KSZ(L),KC
        DZCB(K,ND) = DZCC(K,L)*B(L,K)
      ENDDO

      DBK=0.  
      DO K=KC,KSZ(L),-1
        DBK = DBK + DZCB(K,ND)       
        BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
      ENDDO

      BEC  = 0.
      BI1C = 0.
      BI2C = 0.
      DO K=KC,KSZ(L),-1
        BEC   = BEC  + DZCB(K,ND)
        DZCBK = DZCC(K,L)*BK(K,ND)
        BI1C  = BI1C + DZCBK  
        BI2C  = BI2C + DZCBK + ZZC(K,L)*DZCB(K,ND) 
      ENDDO

      IF( IGRIDV == 0 )THEN
        ! *** STANDARD SIGMA CODE
        BI1W(L) = BI1C
        BI2W(L) = BI2C
        BEW(L)  = BEC
        BI1S(L) = BI1C
        BI2S(L) = BI2C
        BES(L)  = BEC
      
      ELSE   ! *** SGZ FACE INTEGRALS
        IF( KSZFLAT(L) )THEN
          ! *** STANDARD SIGMA CODE
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
        ENDIF
                                           
        ! *** FACE INTEGRALS: WEST
        IF( SUB(L) > 0. )THEN
          BI1W(L) = 0.
          BI2W(L) = 0.
          BEW(L)  = 0.
          DO K=KSZW(L),KC
            DZCB(K,ND)=SGZKW(K,L)*BW(L,K)
          ENDDO

          DBK=0.  
          DO K=KC,KSZW(L),-1
            DBK=DBK+DZCB(K,ND)       
            BK(K,ND)=DBK-0.5*DZCB(K,ND) 
          ENDDO

          DO K=KC,KSZW(L),-1
            BEW(L)  = BEW(L)  + DZCB(K,ND)
            DZCBK   = SGZKW(K,L)*BK(K,ND)
            BI1W(L) = BI1W(L) + DZCBK  
            BI2W(L) = BI2W(L) + DZCBK + ZZW(K,L)*DZCB(K,ND) 
          ENDDO
        ELSE
          BI1W(L) = BI1C
          BI2W(L) = BI2C
          BEW(L)  = BEC
        ENDIF

        ! *** FACE INTEGRALS: EAST
        IF( SUB(LEC(L)) > 0. )THEN
          BI1E(L) = 0.
          BI2E(L) = 0.
          BEE(L)  = 0.
          DO K=KSZE(L),KC
            DZCB(K,ND)=SGZKE(K,L)*BE(L,K)
          ENDDO

          DBK=0.  
          DO K=KC,KSZE(L),-1
            DBK=DBK+DZCB(K,ND)       
            BK(K,ND)=DBK-0.5*DZCB(K,ND) 
          ENDDO

          DO K=KC,KSZE(L),-1
            BEE(L)  = BEE(L)  + DZCB(K,ND)
            DZCBK   = SGZKE(K,L)*BK(K,ND)
            BI1E(L) = BI1E(L) + DZCBK  
            BI2E(L) = BI2E(L) + DZCBK + ZZE(K,L)*DZCB(K,ND) 
          ENDDO
        ELSE
          BI1E(L) = BI1C
          BI2E(L) = BI2C
          BEE(L)  = BEC
        ENDIF

        ! *** FACE INTEGRALS: SOUTH
        IF( SVB(L) > 0. )THEN
          BI1S(L) = 0.
          BI2S(L) = 0.
          BES(L)  = 0.
          DO K=KSZS(L),KC
            DZCB(K,ND)=SGZKS(K,L)*BS(L,K)
          ENDDO

          DBK=0.  
          DO K=KC,KSZS(L),-1
            DBK=DBK+DZCB(K,ND)       
            BK(K,ND)=DBK-0.5*DZCB(K,ND) 
          ENDDO

          DO K=KC,KSZS(L),-1
            BES(L)  = BES(L)  + DZCB(K,ND)
            DZCBK   = SGZKS(K,L)*BK(K,ND)
            BI1S(L) = BI1S(L) + DZCBK  
            BI2S(L) = BI2S(L) + DZCBK + ZZS(K,L)*DZCB(K,ND) 
          ENDDO
        ELSE
          BI1S(L) = BI1C
          BI2S(L) = BI2C
          BES(L)  = BEC
        ENDIF

        ! *** FACE INTEGRALS: NORTH
        IF( SVB(LNC(L)) > 0. )THEN
          BI1N(L) = 0.
          BI2N(L) = 0.
          BEN(L)  = 0.
          DO K=KSZN(L),KC
            DZCB(K,ND)=SGZKN(K,L)*BN(L,K)
          ENDDO

          DBK=0.  
          DO K=KC,KSZN(L),-1
            DBK=DBK+DZCB(K,ND)       
            BK(K,ND)=DBK-0.5*DZCB(K,ND) 
          ENDDO

          DO K=KC,KSZN(L),-1
            BEN(L)  = BEN(L)  + DZCB(K,ND)
            DZCBK   = SGZKN(K,L)*BK(K,ND)
            BI1N(L) = BI1N(L) + DZCBK  
            BI2N(L) = BI2N(L) + DZCBK + ZZN(K,L)*DZCB(K,ND) 
          ENDDO
        ELSE
          BI1N(L) = BI1C
          BI2N(L) = BI2C
          BEN(L)  = BEC
        ENDIF
      ENDIF
      
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  ! *** ADJUST INTEGRALS FOR FACE BLOCKING
  IF( NBLOCKED > 0 .AND. N > 1 )THEN
    ND = 1
    DO LP=1,NBLOCKED
      L = LBLOCKED(LP)
      IF( KSZ(L) == KC ) CYCLE
      
      ! *** LAYER U BLOCKING
      IF( SUB(L) > 0. .AND. (BLDRAFTU(LP)+BLSILLU(LP)) > 0. )THEN
        ! *** WEST FACE INTEGRAL OF L

        ! *** METRICS
        SUMT = 0.
        DZCC(:,L) = 0.
        DO K=KBBU(LP),KC
          DZCC(K,L) = DZC(L,K)
          SUMT = SUMT + DZCC(K,L)
        ENDDO
        DZCC(:,L) = DZCC(:,L)/SUMT

        ZW(L,KBBU(LP)-1) = 0.
        DO K=KBBU(LP),KTBU(LP)
          ZW(L,K)  = ZW(L,K-1) + SGZU(L,K)          ! *** TOP OF LAYER Z FOR OPENING
        ENDDO

        SUMT = 0.
        ZZW(:,L) = 0.
        DO K=KBBU(LP),KC
          SUMT = SUMT + DZCC(K,L)
          ZZW(K,L) = SUMT - 0.5*DZCC(K,L)           ! *** MID LAYER Z TO SURFACE
        ENDDO
        
        ! *** BUOYANCY INTEGRAL
        DO K=KBBU(LP),KC  
          DZCB(K,ND) = DZCC(K,L)*B(L,K)
        ENDDO

        DBK = 0.  
        DO K=KC,KBBU(LP),-1
          DBK = DBK + DZCB(K,ND)       
          BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
        ENDDO
            
        BEC=0.
        BI1C=0.
        BI2C=0.
        DO K=KC,KBBU(LP),-1
          BEC   = BEC  + DZCB(K,ND)
          DZCBK = DZCC(K,L)*BK(K,ND)
          BI1C  = BI1C + DZCBK  
          BI2C  = BI2C + DZCBK + ZZW(K,L)*DZCB(K,ND) 
        ENDDO
        BI1W(L) = BI1C
        BI2W(L) = BI2C
        BEW(L)  = BEC
      
        ! *** EAST FACE INTEGRAL OF LW
        IF( SUBO(L) > 0.0 )THEN
          LW = LWC(L)

          ZE(LW,KBBU(LP)-1) = 0.
          DO K=KBBU(LP),KTBU(LP)
            ZE(LW,K)  = ZE(LW,K-1) + SGZU(L,K)         ! *** TOP OF LAYER Z FOR OPENING
          ENDDO

          SUMT = 0.
          ZZE(:,LW) = 0.
          DO K=KBBU(LP),KC
            SUMT = SUMT + DZCC(K,L)
            ZZE(K,LW) = SUMT - 0.5*DZCC(K,L)           ! *** MID LAYER Z TO SURFACE
          ENDDO

          DO K=KBBU(LP),KC
            DZCB(K,ND) = DZCC(K,L)*B(LW,K)
          ENDDO

          DBK=0.  
          DO K=KC,KBBU(LP),-1
            DBK = DBK + DZCB(K,ND)       
            BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
          ENDDO

          BEC=0.
          BI1C=0.
          BI2C=0.
          DO K=KC,KBBU(LP),-1
            BEC   = BEC  + DZCB(K,ND)
            DZCBK = DZCC(K,L)*BK(K,ND)
            BI1C  = BI1C + DZCBK  
            BI2C  = BI2C + DZCBK + ZZE(K,LW)*DZCB(K,ND) 
          ENDDO
          BI1E(LW) = BI1C
          BI2E(LW) = BI2C
          BEE(LW)  = BEC
        ENDIF
      ENDIF

      ! *** LAYER V BLOCKING
      IF( SVB(L) > 0. .AND. (BLDRAFTV(LP)+BLSILLV(LP)) > 0. )THEN
        ! *** V FACE INTEGRAL 

        ! *** METRICS
        SUMT = 0.
        DZCC(:,L) = 0.
        DO K=KBBV(LP),KC
          DZCC(K,L) = DZC(L,K)
          SUMT = SUMT + DZCC(K,L)
        ENDDO
        DZCC(:,L) = DZCC(:,L)/SUMT

        ZS(L,:)=0.
        DO K=KBBV(LP),KTBV(LP)
          ZS(L,K)  = ZS(L,K-1) + SGZV(L,K)          ! *** TOP OF LAYER Z FOR OPENING
        ENDDO
        
        SUMT = 0.
        ZZS(:,L) = 0.
        DO K=KBBV(LP),KC
          SUMT = SUMT + DZCC(K,L)
          ZZS(K,L) = SUMT - 0.5*DZCC(K,L)           ! *** MID LAYER Z TO SURFACE
        ENDDO
        
        ! *** BUOYANCY INTEGRAL
        DO K=KBBV(LP),KC
          DZCB(K,ND) = DZCC(K,L)*B(L,K)
        ENDDO

        DBK=0.  
        DO K=KC,KBBV(LP),-1
          DBK = DBK + DZCB(K,ND)       
          BK(K,ND) = DBK - 0.5*DZCB(K,ND) 
        ENDDO

        BEC=0.
        BI1C=0.
        BI2C=0.
        DO K=KC,KBBV(LP),-1
          BEC   = BEC  + DZCB(K,ND)
          DZCBK = DZCC(K,L)*BK(K,ND)
          BI1C  = BI1C + DZCBK  
          BI2C  = BI2C + DZCBK + ZZS(K,L)*DZCB(K,ND) 
        ENDDO
        BI1S(L) = BI1C
        BI2S(L) = BI2C
        BES(L)  = BEC

        ! *** NORTH FACE INTEGRAL OF LS
        IF( SVBO(L) > 0.0 )THEN
          LS = LSC(L)
          ZN(LS,KBBV(LP)-1)=0.
          DO K=KBBV(LP),KTBV(LP)
            ZN(LS,K)  = ZN(LS,K-1) +     SGZV(L,K)    ! *** TOP OF LAYER Z
            ZZN(K,LS) = ZN(LS,K)   - 0.5*SGZV(L,K)    ! *** MID LAYER Z
          ENDDO
          
          DO K=KBBV(LP),KTBV(LP)
            DZCB(K,ND) = SGZV(L,K)*B(LS,K)
          ENDDO

          DBK=0.  
          DO K=KTBV(LP),KBBV(LP),-1
            DBK = DBK + DZCB(K,ND)       
            BK(K,ND) = DBK-0.5*DZCB(K,ND) 
          ENDDO

          BEC=0.
          BI1C=0.
          BI2C=0.
          DO K=KTBV(LP),KBBV(LP),-1
            BEC   = BEC  + DZCB(K,ND)
            DZCBK = SGZV(L,K)*BK(K,ND)
            BI1C  = BI1C + DZCBK  
            BI2C  = BI2C + DZCBK + ZZN(K,LS)*DZCB(K,ND) 
          ENDDO
          BI1N(L) = BI1C
          BI2N(L) = BI2C
          BEN(L)  = BEC
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  
  RETURN  

END  

SUBROUTINE INTERPB(L)

  ! *** RETURNS THE B ARRAY INTERPOLATED ONTO EACH FACE AT THE LAYER MIDPOINT
  USE GLOBAL
  
  INTEGER, INTENT(IN) :: L
  INTEGER :: K,KF,LW,LE,LS,LN,NMAX,KK
  REAL :: ZF, XDAT(1:KC),YDAT(1:KC)
  
  ! *** CELL CENTROID VALUES APPLY FOR ALL LAYERS BUT 
  NMAX = KC-KSZ(L)+1
  KK = 0
  DO K=KSZ(L),KC
    KK = KK+1
    BW(L,K) = B(L,K)
    BE(L,K) = B(L,K)
    BS(L,K) = B(L,K)
    BN(L,K) = B(L,K)
    XDAT(KK)= BELV(L)  + HP(L)*ZZC(K,L)
    YDAT(KK)= B(L,K)
  ENDDO
  IF( IGRIDV < 1 ) RETURN

  LW = LWC(L)
  LE = LEC(L)
  LS = LSC(L)
  LN = LNC(L)
  
  ! *** INTERPOLATE VARIN FOR THE WEST FACE
  IF( KSZ(LW) > KSZ(L) )THEN    
    DO KF=KSZW(L),KC
      ZF = BELVW(L) + HPW(L)*ZZW(KF,L)
      CALL INTERPOL(NMAX,XDAT,YDAT,ZF,BW(L,KF))
    ENDDO
  ENDIF

  ! *** INTERPOLATE VARIN FOR THE EAST FACE
  IF( KSZ(LE) > KSZ(L) )THEN
    DO KF=KSZE(L),KC
      ZF = BELVE(L) + HPE(L)*ZZE(KF,L)
      CALL INTERPOL(NMAX,XDAT,YDAT,ZF,BE(L,KF))
    ENDDO
  ENDIF

  ! *** INTERPOLATE VARIN FOR THE SOUTH FACE
  IF( KSZ(LS) > KSZ(L) )THEN   
    DO KF=KSZS(L),KC
      ZF = BELVS(L) + HPS(L)*ZZS(KF,L)
      CALL INTERPOL(NMAX,XDAT,YDAT,ZF,BS(L,KF))
    ENDDO
  ENDIF

  ! *** INTERPOLATE VARIN FOR THE NORTH FACE
  IF( KSZ(LN) > KSZ(L) )THEN
    DO KF=KSZN(L),KC
      ZF = BELVN(L) + HPN(L)*ZZN(KF,L)
      CALL INTERPOL(NMAX,XDAT,YDAT,ZF,BN(L,KF))
    ENDDO
  ENDIF
  
  CONTAINS
  
  SUBROUTINE INTERPOL(NMAX,XDAT,YDAT,XVAL,YVAL)
    ! ** INTERPOLATION FOR YVAL OF XVAL
    INTEGER, INTENT(IN ) :: NMAX
    REAL,    INTENT(IN ) :: XVAL,XDAT(NMAX),YDAT(NMAX)
    REAL,    INTENT(OUT) :: YVAL
    INTEGER       :: N

    IF(  XVAL < XDAT(1) )THEN 
      YVAL = YDAT(1)
    ELSEIF ( XVAL > XDAT(NMAX) )THEN 
      YVAL = YDAT(NMAX)
    ELSE 
      DO N=1,NMAX-1
        IF( XVAL >= XDAT(N) .AND. XVAL <= XDAT(N+1) )THEN 
          YVAL  =  (YDAT(N+1)-YDAT(N))*(XVAL-XDAT(N)) /(XDAT(N+1)-XDAT(N))+YDAT(N)
          EXIT
        ENDIF 
      ENDDO
    ENDIF
  END SUBROUTINE
  
 END SUBROUTINE
