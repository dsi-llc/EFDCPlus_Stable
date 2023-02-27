! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALPNHS

  ! CHANGE RECORD
  ! **  SUBROUTINE CALPNHS CALCULATES QUASI-NONHYDROSTATIC PRESSURE

  ! *** PNHYDS HAS UNITS OF M2/S2
  
  USE GLOBAL

  IMPLICIT NONE
  
  INTEGER :: K,L,LE,LW,LN,LS,NWR,NS,IU,JU,KU,LU,ID,JD,KD,LD,LP
  
  REAL    :: DELTD2,UHUW,VHVW,WB,QMF,QUMF,ADIFF,TMPANG,TMPVAL
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: PNHYDSS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: FWJET
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: WZ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: WZ1
  REAL :: QWRABS

  IF(  .NOT. ALLOCATED(PNHYDSS) )THEN
    ALLOCATE(PNHYDSS(LCM,KCM))
    ALLOCATE(FWJET(LCM,KCM))
    ALLOCATE(WZ(LCM,0:KCM))
    ALLOCATE(WZ1(LCM,0:KCM))
    PNHYDSS=0.0
    FWJET=0.0
    WZ = 0.0
    WZ1 = 0.0
  ENDIF

  IF( ISDYNSTP == 0 )THEN
    DELT=DT
    DELTD2=0.5*DT
    DELTI=1./DELT
  ELSE
    DELT=DTDYN
    DELTD2=0.5*DTDYN
    DELTI=1./DELT
  END IF

  ! **  CALCULATE THE PHYSICAL VERTICAL VELOCIY
  DO L=2,LA
    LE=LEC(L)
    LW=LWC(L)
    LN=LNC(L)
    LS=LSC(L)
    WZ(L,0)=DELTI*(BELV(L)-BELV1(L))
    WZ(L,KC)=GI*( DELTI*(P(L)-P1(L)) + 0.5*U(LE,KC)*(P(LE)-P(L))*DXIU(LE) + 0.5*U(L,KC)*(P(L)-P(LW))*DXIU(L) &
                                     + 0.5*V(LN,KC)*(P(LN)-P(L))*DYIV(LN) + 0.5*V(L,KC)*(P(L)-P(LS))*DYIV(L) )
  ENDDO
  
  IF( KC > 2 )THEN
    DO K=1,KS
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)
        LE=LEC(L)
        LW=LWC(L)
        LN=LNC(L)
        LS=LSC(L)
        WZ(L,K)=W(L,K)+GI*ZZ(L,K)*( DELTI*(P(L)-P1(L)) &
                                  + 0.5*U(LE,K)*(P(LE)-P(L))*DXIU(LE) + 0.5*U(L,K)*(P(L)-P(LW))*DXIU(L) &
                                  + 0.5*V(LN,K)*(P(LN)-P(L))*DYIV(LN) + 0.5*V(L,K)*(P(L)-P(LS))*DYIV(L) ) &
                    + (1.-ZZ(L,K))*( DELTI*(BELV(L)-BELV1(L)) &
                                  + 0.5*U(LE,K)*(BELV(LE)-BELV(L))*DXIU(LE) + 0.5*U(L,K)*(BELV(L)-BELV(LW))*DXIU(L) &
                                  + 0.5*V(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN) + 0.5*V(L,K)*(BELV(L)-BELV(LS))*DYIV(L) )
      ENDDO
    ENDDO
  ENDIF

  ! **  CALCULATE FLUXES
  DO K=1,KC
    DO L=1,LC
      PNHYDSS(L,K) = PNHYDS(L,K)
      FUHU(L,K) = 0.
      FVHU(L,K) = 0.
      FWQQ(L,KC) = 0.
    ENDDO
  ENDDO
  DO K=1,KS
    DO L=2,LA
      LW=LWC(L)
      LS=LSC(L)
      UHUW = 0.5*(UHDY(L,K)+UHDY(L,K+1))
      VHVW = 0.5*(VHDX(L,K)+VHDX(L,K+1))
      FUHU(L,K) = MAX(UHUW,0.)*WZ(LW,K) + MIN(UHUW,0.)*WZ(L,K)
      FVHU(L,K) = MAX(VHVW,0.)*WZ(LS,K)  + MIN(VHVW,0.)*WZ(L,K)
    ENDDO
  ENDDO
  DO K=1,KC
    DO L=2,LA
      WB = 0.5*DXYP(L)*(W(L,K-1)+W(L,K))
      FWQQ(L,K) = MAX(WB,0.)*WZ(L,K-1) + MIN(WB,0.)*WZ(L,K)
      FWJET(L,K) = 0.
    ENDDO
  ENDDO

  ! **  ADD RETURN FLOW MOMENTUM FLUX
  DO NWR=1,NQWR
    IF( WITH_RET(NWR).NQWRMFU > 0 )THEN
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS=WITH_RET(NWR).NQWRSERQ
      IF( QWRSERT(NS) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=WITH_RET(NWR).IQWRU
        JU=WITH_RET(NWR).JQWRU
        KU=WITH_RET(NWR).KQWRU
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        IU=WITH_RET(NWR).IQWRD
        JU=WITH_RET(NWR).JQWRD
        KU=WITH_RET(NWR).KQWRD
      ENDIF
      QWRABS = ABS(QWRSERT(NS))
      LU=LIJ(IU,JU)

      QMF=WITH_RET(NWR).QWR+QWRABS
      QUMF=QMF*QMF/(H1P(LU)*DZC(LU,KU)*WITH_RET(NWR).BQWRMFU)
      IF( WITH_RET(NWR).NQWRMFU == 1 )  FWJET(LU     ,KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == 2 )  FWJET(LU     ,KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == 3 )  FWJET(LU+1   ,KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == 4 )  FWJET(LNC(LU),KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == -1 ) FWJET(LU     ,KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == -2 ) FWJET(LU     ,KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == -3 ) FWJET(LU+1   ,KU)=-QUMF
      IF( WITH_RET(NWR).NQWRMFU == -4 ) FWJET(LNC(LU),KU)=-QUMF
    ENDIF
    IF( WITH_RET(NWR).NQWRMFD > 0 )THEN
      ID=WITH_RET(NWR).IQWRD
      JD=WITH_RET(NWR).JQWRD
      KD=WITH_RET(NWR).KQWRD
      LD=LIJ(ID,JD)
      ADIFF=ABS(WITH_RET(NWR).ANGWRMFD-90.)
      IF( ADIFF < 1.0 )THEN
        TMPANG=1.
      ELSE
        TMPANG=0.017453*WITH_RET(NWR).ANGWRMFD
        TMPANG=SIN(TMPANG)
      ENDIF
      NS=WITH_RET(NWR).NQWRSERQ
      QMF=WITH_RET(NWR).QWR+QWRSERT(NS)
      QUMF=TMPANG*QMF*QMF/(H1P(LD)*DZC(LD,KD)*WITH_RET(NWR).BQWRMFD)
      IF( WITH_RET(NWR).NQWRMFD == 1 )  FWJET(LD     ,KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == 2 )  FWJET(LD     ,KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == 3 )  FWJET(LD+1   ,KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == 4 )  FWJET(LNC(LD),KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == -1) FWJET(LD     ,KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == -2) FWJET(LD     ,KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == -3) FWJET(LD+1   ,KD)=QUMF
      IF( WITH_RET(NWR).NQWRMFD == -4) FWJET(LNC(LD),KD)=QUMF
    ENDIF
  ENDDO

  ! **  CALCULATE QUASI-NONHYDROSTATIC PRESSURE
  DO L=2,LA
    LE=LEC(L)
    LN=LNC(L)
    TMPVAL=0.5*DZC(L,KC)/DXYP(L)
    PNHYDS(L,KC)= 0.75*TMPVAL*( DELTI*DXYP(L)*(HP(L)*WZ(L,KC)-H1P(L)*WZ1(L,KC)) + FUHU(LE,KC)-FUHU(L,KC)+FVHU(LN,KC)-FVHU(L,KC) )            &
                + 0.25*TMPVAL*( DELTI*DXYP(L)*(HP(L)*WZ(L,KS)-H1P(L)*WZ1(L,KS)) + FUHU(LE,KS)-FUHU(L,KS)+FVHU(LN,KS)-FVHU(L,KS) ) -FWQQ(L,KC)
  ENDDO
  
  DO K=KS,1,-1
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)
      LE=LEC(L)
      LN=LNC(L)
      TMPVAL=0.5*(DZC(L,K+1)+DZC(L,K))/DXYP(L)
      PNHYDS(L,K) = PNHYDS(L,K+1) + FWQQ(L,K+1)-FWQQ(L,K) - FWJET(L,K)   &
                                  + TMPVAL*( DELTI*DXYP(L)*(HP(L)*WZ(L,K)-H1P(L)*WZ1(L,K)) + FUHU(LE,K)-FUHU(L,K)+FVHU(LN,K)-FVHU(L,K) )
    ENDDO
  ENDDO
  
  DO K=0,KC
    DO L=2,LA
      WZ1(L,K)=WZ(L,K)
    ENDDO
  ENDDO
  
  DO K=1,KC
    DO L=1,LC
      PNHYDS(L,K) = 0.5*( PNHYDSS(L,K)+PNHYDS(L,K) )
    ENDDO
  ENDDO
  
  RETURN

END

