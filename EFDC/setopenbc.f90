! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SETOPENBC(DELTD2, HUT, HVT)
 
  ! *** SUBROUTINE SETOPENBC SETS OPEN BOUNDARY CONDITIONS FOR CALPUV2C AND CALPUV9C

  ! *** Modified by Paul M. Craig to include outgoing wave (anti-reflection) BC type
  ! *** Incoming wave radiation type BC (Bennett and McIntosh, 1982) adjusted for
  ! *** for tidal offsets and datum shifts to only include the displacement heights.

  ! *** User Specified Elevation Note
  ! ***   For older versions of EFDC, including the 2020 GVC code, SETOPENBC uses only the 
  ! ***   celerity based term CxT for the open BC face to assign the adjacent cell's FP. 
  ! ***   Later versions of EFDC add the current FP of the adjacent cell plus the CxT component.  
  ! ***   The newer approach allows for non-zero elevation offsets to be used.
  ! ***   
  ! ***   For example, for West open boundary cells EFDC+ uses
  ! ***      FP(LEC(L)) = FP(LEC(L)) + CET*FP(L)                            ! *** m4/s3
  
  ! *** DELTD2 is DELT/2     (s)
  ! *** DELTI  is 1/DELT,  Passed via the GLOBAL module     (1/s)        
  ! *** HUT(LCM), HVT(LCM) are the depths at the cell interfaces of U and V respectively   (m)
  ! *** 
  
  USE GLOBAL
  
  IMPLICIT NONE

  ! *** Passed in
  REAL, INTENT(IN) :: DELTD2
  REAL, INTENT(IN) :: HUT(LCM), HVT(LCM)
  
  ! *** Local Variables
  INTEGER :: L, LE, LN, M, LL
  REAL    :: C1, HDRY2, TM, TMPVAL, TC, TS, FP1G, TMP, CET, CWT, CNT, CST
  REAL    :: FPE, FPW, FPN, FPS, C, R, UVEL, VVEL

  C1 = 0.5*G
  HDRY2 = 0.2*HDRY
  
  ! **  SET OPEN BOUNDARY SURFACE ELEVATIONS
  
  IF( ISDYNSTP == 0 )THEN
    TN=DT*FLOAT(N)+TCON*TBEGIN
  ELSE
    TN=TIMESEC
  ENDIF
  DO M=1,MTIDE
    TM=MOD(TN,TCP(M))
    TM=PI2*TM/TCP(M)
    CCCOS(M)=COS(TM)
    SSSIN(M)=SIN(TM)
  ENDDO

  ! *** WEST OPEN BOUNDARY
  DO LL=1,NPBW
    L = LPBW(LL)
    LE = LEC(L)
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    IF( LOPENBCDRY(L) )CYCLE

    FP(L) = PSERT(NPSERW(LL)) + 0.5*PSERZDF(NPSERW(LL)) + PSERST(NPSERW(LL)) + 0.5*PSERZDS(NPSERW(LL))         ! *** m2/s2
    IF( NPFORT >= 1 .AND. NPSERW1(LL) > 0 )THEN
      TMPVAL = PSERT(NPSERW1(LL)) + 0.5*PSERZDF(NPSERW1(LL)) + PSERST(NPSERW1(LL)) + 0.5*PSERZDS(NPSERW1(LL))
      FP(L) = FP(L) + TPCOORDW(LL)*(TMPVAL-FP(L))
    ENDIF
    DO M = 1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
    ENDDO

    IF( ISPBW(LL) == 1 .OR. ISPBW(LL) == 2 )THEN
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERW(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation

      ! ***     (s)   (m/s2) (-)  (-)     (m)
      CET = 0.5*DELTD2*G*HRUO(LE)*RCX(LE)*HUT(LE)                    ! *** m2/s

      ! ***   (s)             (m/s)          (1/m)  
      TMP = DELTD2*SQRT(G*HUT(LE))*DXIU(LE)                          ! *** dimensionless
      CC(L) = CET*(1.+TMP)/TMP                                       ! *** m2/s
      CE(L) = -CET                                                   ! *** m2/s

      ! ***   m2/s    m2/s2          m/s         m3/s     1/m     1/m     (-)
      FP(L) = CET*(2.*FP(L) - SQRT(G*HUT(LE))*FUHDYE(LE)*DYIU(LE)/HUT(LE))/TMP    ! *** m4/s3

      FP(L) = FP(L) + CET*PSERAVG(NPSERW(LL),1)/TMP                  ! *** Add back offset

    ELSEIF( ISPBW(LL) == 4 .OR. ISPBW(LL) == 5 )THEN                 
      ! *** Outgoing wave (anti-reflection) type                     
                                                                     
      ! ***     (s)   (m/s2) (-)   (-)     (m)
      CET = 0.5*DELTD2*G*HRUO(LE)*RCX(LE)*HUT(LE)                    ! *** m2/s

      ! ***   (s)             (m/s)          (1/m)  
      TMP = DELTD2*SQRT(G*HUT(LE))*DXIU(LE)                          ! *** dimensionless
      CC(L) = CET*(1.+TMP)/TMP                                       ! *** m2/s
      CE(L) = -CET
      
      ! *** Setup adjacent cell
      ! ***        m4/s3        m2/s m2/s2
      FP(LE) = FP(LE) + CET*FP(L)                                    ! *** m4/s3
      CW(LE) = 0.

      ! **   m2/s         m/s       m3/s     1/m    1/m     (-)
      FPW = -CET*(SQRT(G*HUT(LE))*FUHDYE(LE)*DYIU(L)/HUT(LE))/TMP    ! *** m4/s3

      FP(L) = FPW + CET*FP(L)/TMP                                    ! *** Add specified elevation

    ELSE
      ! *** INACTIVATE BC'S WHEN ELEVATIONS DROP BELOW BOTTOM+HDRY
      FP1G = FP(L)/G-HDRY2
      IF( FP1G < BELV(L) .OR. FP1G < BELV(LE) )THEN
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CET = 0.
        LOPENBCDRY(L) = .TRUE.
      ELSE
        CET = 0.5*DELTD2*G*HRUO(LE)*RCX(LE)*HUT(LE)
      ENDIF
      
      ! *** Setup adjacent cell
      ! ***        m4/s3        m2/s m2/s2
      FP(LE) = FP(LE) + CET*FP(L)                                    ! *** m4/s3
      CW(LE) = 0.

      ! ***    1/s   m2    m2/s2
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
    ENDIF
  ENDDO

  ! *** EAST OPEN BOUNDARY
  DO LL=1,NPBE
    L = LPBE(LL)
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    
    IF( LOPENBCDRY(L) )CYCLE
    
    FP(L) = PSERT(NPSERE(LL)) + 0.5*PSERZDF(NPSERE(LL)) + PSERST(NPSERE(LL)) + 0.5*PSERZDS(NPSERE(LL))         ! *** m2/s2
    IF( NPFORT >= 1 .AND. NPSERE1(LL) > 0 )THEN
      TMPVAL = PSERT(NPSERE1(LL)) +0.5*PSERZDF(NPSERE1(LL)) + PSERST(NPSERE1(LL)) + 0.5*PSERZDS(NPSERE1(LL))
      FP(L) = FP(L)+TPCOORDE(LL)*(TMPVAL-FP(L))
    ENDIF
    
    DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
    ENDDO

    IF( ISPBE(LL) == 1 .OR. ISPBE(LL) == 2 )THEN
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERE(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation
                                                                     
      ! ***        (s)  (m/s2) (-)     (-)    (m)                    
      CWT = 0.5*DELTD2* G*     HRUO(L)*RCX(L)*HUT(L)                 ! *** m2/s
                                                                     
      ! ***   (s)  (m/s)          (1/m)                              
      TMP = DELTD2*SQRT(G*HUT(L))*DXIU(L)                            ! *** dimensionless
      CC(L) = CWT*(1.+TMP)/TMP                                       ! *** m2/s
      CW(L) = -CWT                                                   ! *** m2/s
      
      ! ***   m2/s    m2/s2          m/s       m3/s     1/m    1/m     (-)
      FP(L) = CWT*(2.*FP(L) + SQRT(G*HUT(L))*FUHDYE(L)*DYIU(L)/HUT(L))/TMP       ! *** m4/s3

      FP(L) = FP(L) + CWT*PSERAVG(NPSERE(LL),1)/TMP                  ! *** Add back offset
                                                                     
    ELSEIF( ISPBE(LL) == 4 .OR. ISPBE(LL) == 5 )THEN                 
      ! *** Outgoing wave (anti-reflection) type                     
                                                                     
      ! ***        (s)  (m/s2) (-)     (-)    (m)                    
      CWT = 0.5*DELTD2* G*     HRUO(L)*RCX(L)*HUT(L)                 ! *** m2/s
                                                                     
      ! ***   (s)  (m/s)          (1/m)                              
      TMP = DELTD2*SQRT(G*HUT(L))*DXIU(L)                            ! *** dimensionless
      CC(L) = CWT*(1.+TMP)/TMP                                       ! *** m2/s
      CW(L) = -CWT                                                   ! *** m2/s
      
      ! *** Setup adjacent cell
      ! ***        m4/s3        m2/s m2/s2                           
      FP(LWC(L)) = FP(LWC(L)) + CWT*FP(L)                            ! *** m4/s3
      CE(LWC(L)) = 0.0
      
      ! **  m2/s         m/s       m3/s     1/m    1/m     (-)
      FPE = CWT*( SQRT(G*HUT(L))*FUHDYE(L)*DYIU(L)/HUT(L))/TMP       ! *** m4/s3

      FP(L) = FPE + CWT*FP(L)/TMP                                    ! *** Add specified elevation

    ELSE
      ! *** INACTIVATE BC'S WHEN ELEVATIONS DROP BELOW BOTTOM+HDRY
      FP1G = FP(L)/G-HDRY2
      IF( FP1G < BELV(L) .OR. FP1G < BELV(LWC(L)) )THEN
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CWT = 0.                                                     
        LOPENBCDRY(L) = .TRUE.                                       
      ELSE                                                           
        ! ***        (s)  (m/s2) (-)     (-)    (m)                  
        CWT = 0.5*DELTD2* G*     HRUO(L)*RCX(L)*HUT(L)               ! *** m2/s
      ENDIF
      
      ! *** Setup adjacent cell
      ! ***        m4/s3        m2/s m2/s2                           
      FP(LWC(L)) = FP(LWC(L)) + CWT*FP(L)                            ! *** m4/s3
      CE(LWC(L)) = 0.0
                                                                     
      ! ***    1/s   m2    m2/s2                                     
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
    ENDIF
  ENDDO

  ! *** SOUTH OPEN BOUNDARY
  DO LL=1,NPBS
    L = LPBS(LL)
    LN = LNC(L)
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    IF( LOPENBCDRY(L) )CYCLE

    FP(L) = PSERT(NPSERS(LL)) + 0.5*PSERZDF(NPSERS(LL)) + PSERST(NPSERS(LL)) + 0.5*PSERZDS(NPSERS(LL))         ! *** m2/s2
    IF( NPFORT >= 1 .AND. NPSERS1(LL) > 0 )THEN
      TMPVAL = PSERT(NPSERS1(LL)) + 0.5*PSERZDF(NPSERS1(LL)) + PSERST(NPSERS1(LL)) + 0.5*PSERZDS(NPSERS1(LL))
      FP(L) = FP(L)+TPCOORDS(LL)*(TMPVAL-FP(L))
    ENDIF
    DO M=1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
    ENDDO
    
    IF( ISPBS(LL) == 1 .OR. ISPBS(LL) == 2 )THEN
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERS(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation
                                                                     
      ! ***      (s)  (m/s2) (-)     (-)    (m)                      
      CNT = 0.5*DELTD2*G   *HRVO(LN)*RCY(LN)*HVT(LN)                 ! *** m2/s
                                                                     
      ! ***   (s)         (m/s)     (1/m)                            
      TMP = DELTD2*SQRT(G*HVT(LN))*DYIV(LN)                          ! *** dimensionless
      CC(L) = CNT*(1.+TMP)/TMP                                       ! *** m2/s
      CN(L) = -CNT                                                   ! *** m2/s
      
      ! ***   m2/s      m2/s2          m/s       m3/s      1/m     1/m       (-)
      FP(L) = CNT*( 2.*FP(L) - SQRT(G*HVT(LN))*FVHDXE(LN)*DXIV(LN)/HVT(LN) )/TMP       ! *** m4/s3
      
      FP(L) = FP(L) + CNT*PSERAVG(NPSERS(LL),1)/TMP                  ! *** Add back offset
      
    ELSEIF( ISPBS(LL) == 4 .OR. ISPBS(LL) == 5 )THEN
      ! *** Outgoing wave (anti-reflection) type                     
      
      ! ***      (s)  (m/s2) (-)     (-)    (m)                      
      CNT = 0.5*DELTD2*G   *HRVO(LN)*RCY(LN)*HVT(LN)                 ! *** m2/s
                                                                     
      ! ***   (s)         (m/s)     (1/m)                            
      TMP = DELTD2*SQRT(G*HVT(LN))*DYIV(LN)                          ! *** dimensionless
      CC(L) = CNT*(1.+TMP)/TMP                                       ! *** m2/s
      CN(L) = -CNT                                                   ! *** m2/s
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LN) = FP(LN) + CNT*FP(L)                                    ! *** m4/s3
      CS(LN) = 0.


      ! ***   m2/s         m/s       m3/s      1/m     1/m       (-)
      FPS = -CNT*( SQRT(G*HVT(LN))*FVHDXE(LN)*DXIV(LN)/HVT(LN) )/TMP       ! *** m4/s3
      
      FP(L) = FPS + CNT*FP(L)/TMP                                    ! *** Add specified elevation
      
    ELSE
      ! *** INACTIVATE BC'S WHEN ELEVATIONS DROP BELOW BOTTOM+HDRY
      FP1G = FP(L)/G-HDRY2
      IF( FP1G < BELV(L) .OR. FP1G < BELV(LN) )THEN
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CNT = 0.
        LOPENBCDRY(L) = .TRUE.
      ELSE
        ! ***        (s)  (m/s2) (-)     (-)    (m)
        CNT = 0.5*DELTD2*G    *HRVO(LN)*RCY(LN)*HVT(LN)              ! *** m2/s
      ENDIF

      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LN) = FP(LN) + CNT*FP(L)                                    ! *** m4/s3
      CS(LN) = 0.

      ! ***    1/s   m2    m2/s2
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
    ENDIF
  ENDDO

  ! *** NORTH OPEN BOUNDARY
  DO LL=1,NPBN
    L = LPBN(LL)
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    IF( LOPENBCDRY(L) )CYCLE
    
    FP(L) = PSERT(NPSERN(LL)) + 0.5*PSERZDF(NPSERN(LL)) + PSERST(NPSERN(LL)) + 0.5*PSERZDS(NPSERN(LL))         ! *** m2/s2
    IF( NPFORT >= 1 .AND. NPSERN1(LL) > 0 )THEN
      TMPVAL = PSERT(NPSERN1(LL)) + 0.5*PSERZDF(NPSERN1(LL)) + PSERST(NPSERN1(LL)) + 0.5*PSERZDS(NPSERN1(LL))
      FP(L) = FP(L)+TPCOORDN(LL)*(TMPVAL-FP(L))
    ENDIF
    DO M=1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
    ENDDO

    IF( ISPBN(LL) == 1 .OR. ISPBN(LL) == 2 )THEN
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERN(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation
                                                                     
      ! ***      (s)  (m/s2) (-)     (-)    (m)                      
      CST = 0.5*DELTD2*G   *HRVO(L)*RCY(L)*HVT(L)                    ! *** m2/s

      ! ***   (s)         (m/s)     (1/m)                            
      TMP = DELTD2*SQRT(G*HVT(L))*DYIV(L)                            ! *** dimensionless
      CC(L) = CST*(1.+TMP)/TMP                                       ! *** m2/s
      CS(L) = -CST                                                   ! *** m2/s
      
      ! ***   m2/s      m2/s2         m/s       m3/s    1/m     1/m     (-)
      FP(L) = CST*( 2.*FP(L) + SQRT(G*HVT(L))*FVHDXE(L)*DXIV(L)/HVT(L))/TMP     ! *** m4/s3
      
      FP(L) = FP(L) + CNT*PSERAVG(NPSERN(LL),1)/TMP                  ! *** Add back offset
      
    ELSEIF( ISPBN(LL) == 4 .OR. ISPBN(LL) == 5 )THEN
      ! *** Outgoing wave (anti-reflection) type                     
      
      ! ***      (s)  (m/s2) (-)     (-)    (m)                      
      CST = 0.5*DELTD2*G   *HRVO(L)*RCY(L)*HVT(L)                    ! *** m2/s

      ! ***   (s)         (m/s)     (1/m)                            
      TMP = DELTD2*SQRT(G*HVT(L))*DYIV(L)                            ! *** dimensionless
      CC(L) = CST*(1.+TMP)/TMP                                       ! *** m2/s
      CS(L) = -CST                                                   ! *** m2/s
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LSC(L)) = FP(LSC(L)) + CST*FP(L)                            ! *** m4/s3
      CN(LSC(L)) = 0.

      ! *** m2/s         m/s       m3/s    1/m     1/m     (-)
      FPN = CST*( SQRT(G*HVT(L))*FVHDXE(L)*DXIV(L)/HVT(L))/TMP       ! *** m4/s3
      
      FP(L) = FPN + CST*FP(L)/TMP                                    ! *** Add specified elevation
      
    ELSE
      ! *** INACTIVATE BC'S WHEN ELEVATIONS DROP BELOW BOTTOM+HDRY
      FP1G = FP(L)/G-HDRY2
      IF( FP1G < BELV(L) .OR. FP1G < BELV(LSC(L)) )THEN
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CST = 0.
        LOPENBCDRY(L) = .TRUE.
      ELSE
        ! ***        (s)  (m/s2) (-)     (-)    (m)
        CST = 0.5*DELTD2*G     *HRVO(L)*RCY(L)*HVT(L)                ! *** m2/s
      ENDIF

      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LSC(L)) = FP(LSC(L)) + CST*FP(L)                            ! *** m4/s3
      CN(LSC(L)) = 0.

      ! ***    1/s   m2    m2/s2
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
    ENDIF
  ENDDO

  RETURN
END

