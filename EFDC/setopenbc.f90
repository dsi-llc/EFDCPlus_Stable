! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SETOPENBC(DELTD2, HUT, HVT, NCORDRY)
 
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
  
  use GLOBAL
  use Allocate_Initialize
  
  implicit none

  ! *** Passed in
  integer, intent(IN) :: NCORDRY
  real, intent(IN)    :: DELTD2
  real, intent(IN)    :: HUT(LCM), HVT(LCM)
  
  ! *** Local Variables
  integer :: IBC, L, LE, LN, LS, LW, LL, M
  integer, save :: ISRUNNING, ISSUMQ
  integer, save,allocatable,dimension(:) :: RESTYPE          ! *** Type of Residual flow at open BC cells
  
  real    :: C1, HDRY2, HDRY5, HDRY10, TM, TMPVAL, TC, TS, FP1G, FACTOR, TMP, CET, CWT, CNT, CST
  real    :: ELEV, FPE, FPW, FPN, FPS, C, R, UVEL, VVEL
  real(RKD) :: VOLSUM, VOLINC, SUMQ, RINC
  real :: TN             !< 
  
  real(RKD), save :: DAYOLD, SUMDAY

  real(RKD),save,allocatable,dimension(:) :: OBCRESIDUALQ    ! *** Residual flow at open BC cells (m3/s)
  real(RKD),save,allocatable,dimension(:) :: OBCRESIDUALQ1   ! *** Volume (m3)
  real(RKD),save,allocatable,dimension(:) :: OBCRESIDUALQ2   ! *** Incremental flow at previous time step (m3/s)
  real,save,allocatable,dimension(:) :: CCCOS                     ! 
  real,save,allocatable,dimension(:) :: SSSIN         
  
  if( .not. allocated(OBCRESIDUALQ) )then
    call AllocateDSI( OBCRESIDUALQ,  NBCSOP, 0.0 )
    call AllocateDSI( OBCRESIDUALQ1, NBCSOP, 0.0 )
    call AllocateDSI( OBCRESIDUALQ2, NBCSOP, 0.0 )
    call AllocateDSI( RESTYPE,       NBCSOP,   0 )
    call AllocateDSI( CCCOS,         MTM,    0.0 )
    call AllocateDSI( SSSIN,         MTM,    0.0)
                             
    SUMDAY = 0.0

    ! *** Open BC's in order of NBCSOP assignment
    IBC = 0
    do LL = 1,NPBS
      IBC = IBC + 1
      if( ISPBS(LL) > 3 .and. ISPRS(LL) < 2 ) RESTYPE(IBC) = ISPRS(LL) + 1
    enddo                                                                    
    do LL = 1,NPBW                                                           
      IBC = IBC + 1
      if( ISPBW(LL) > 3 .and. ISPRW(LL) < 2 ) RESTYPE(IBC) = ISPRW(LL) + 1
    enddo                                                                    
    do LL = 1,NPBE                                                           
      IBC = IBC + 1
      if( ISPBE(LL) > 3 .and. ISPRE(LL) < 2 ) RESTYPE(IBC) = ISPRE(LL) + 1
    enddo                                                                    
    do LL = 1,NPBN                                                           
      IBC = IBC + 1
      if( ISPBN(LL) > 3 .and. ISPRN(LL) < 2 ) RESTYPE(IBC) = ISPRN(LL) + 1
    enddo   
        
    ! *** Check for outgoing wave and residual flow method
    ISRUNNING = 0
    if( ANY(ISPBS > 3) .and. ANY(ISPRS == 2) ) ISRUNNING = 1
    if( ANY(ISPBW > 3) .and. ANY(ISPRW == 2) ) ISRUNNING = 1
    if( ANY(ISPBE > 3) .and. ANY(ISPRE == 2) ) ISRUNNING = 1
    if( ANY(ISPBN > 3) .and. ANY(ISPRN == 2) ) ISRUNNING = 1
    
    ISSUMQ = 0
    if( ANY(ISPBS > 3) .and. ANY(ISPRS == 1) ) ISSUMQ = 1
    if( ANY(ISPBW > 3) .and. ANY(ISPRW == 1) ) ISSUMQ = 1
    if( ANY(ISPBE > 3) .and. ANY(ISPRE == 1) ) ISSUMQ = 1
    if( ANY(ISPBN > 3) .and. ANY(ISPRN == 1) ) ISSUMQ = 1
  endif
  
  ! *** Check for radiation boundary types.  Set residual flows, if needed
  if( NCORDRY == 0 )then
    if( ISSUMQ > 0 )then
      ! *** Sum all inflow/outflows from domain.  Used for small domains and test cases.
      SUMQ = 0.0
      do IBC = 1,NBCS
        L = LBCS(IBC)
        if( ANY(LOBCS == L) ) CYCLE                                  ! *** Exclude OBC's
        SUMQ = SUMQ + QSUME(L)
      enddo
    
      ! *** Total volume of OBC's
      VOLSUM = 0.0
      do IBC = 1,NBCSOP
        if( RESTYPE(IBC) == 2 )then
          L = LOBCS(IBC)
          VOLINC = DXYP(L)*HP(L)
          VOLSUM = VOLSUM + VOLINC         
        endif
      enddo
    
      ! *** Distribute the total to the cells
      OBCRESIDUALQ1 = OBCRESIDUALQ                                   ! *** Previous time step OBC flows
      OBCRESIDUALQ = 0.0
      if( VOLSUM > 0.0 )then
        do IBC = 1,NBCSOP
          if( RESTYPE(IBC) == 2 )then          
            L = LOBCS(IBC)
            VOLINC = DXYP(L)*HP(L)
            OBCRESIDUALQ(IBC) = VOLINC/VOLSUM*SUMQ
            OBCRESIDUALQ(IBC) = (OBCRESIDUALQ1(IBC) + 2.*OBCRESIDUALQ(IBC))/3.
          endif
        enddo
      endif
    endif
    
    if( ISRUNNING > 0 )then
      ! *** Sum flows of current period for application to the next tidal period.
      RINC = 0.0
      if( IS2TIM == 0 )then
        if( ISTL == 3 )then
          RINC = 2.*DELTD2
        endif
      else
        RINC = 2.*DELTD2
      endif
    
      OBCRESIDUALQ2 = OBCRESIDUALQ                                   ! *** Previous time step OBC flows
      if( RINC > 0.0 )then
      
        ! *** Convert flows to volumes
        do IBC = 1,NBCSOP
          OBCRESIDUALQ1(IBC) = OBCRESIDUALQ(IBC)*SUMDAY              ! *** m3    Total volume at beginning
        enddo
        
        ! *** Remove one increment
        if( SUMDAY >= TIDALP )then      
          do IBC = 1,NBCSOP
            if( LOPENBCDRY(L) ) CYCLE
            VOLSUM  = OBCRESIDUALQ1(IBC)                             ! *** m3    Total volume at beginning
            OBCRESIDUALQ1(IBC) = (VOLSUM - RINC*OBCRESIDUALQ(IBC))   ! *** m3    Intermediate volume
          enddo
          SUMDAY = TIDALP - RINC
        endif
      
        ! *** Open BC's in order of NBCSOP assignment
        IBC = 0
        do LL = 1,NPBS
          IBC = IBC + 1
          if( ISPBS(LL) > 3 .and. ISPRS(LL) == 0 )then
            L = LIJ(IPBS(LL),JPBS(LL))
            if( LOPENBCDRY(L) ) CYCLE
            VOLINC = DELTD2*(FVHDXE(LNC(L)) + OBCRESIDUALQ(IBC))     ! *** m3    Incremental volume for current time step
            OBCRESIDUALQ(IBC) = VOLINC + OBCRESIDUALQ1(IBC)          ! *** m3    Total volume at end    
          endif
        enddo                                                                    
        do LL = 1,NPBW                                                           
          IBC = IBC + 1                                                          
          if( ISPBW(LL) > 3 .and. ISPRW(LL) == 0 )then
            L = LIJ(IPBW(LL),JPBW(LL))                                             
            if( LOPENBCDRY(L) ) CYCLE
            VOLINC = DELTD2*(FUHDYE(LEC(L)) + OBCRESIDUALQ(IBC))     ! *** m3    Incremental volume for current time step
            OBCRESIDUALQ(IBC) = VOLINC + OBCRESIDUALQ1(IBC)          ! *** m3    Total volume at end          
          endif
        enddo                                                                    
        do LL = 1,NPBE                                                           
          IBC = IBC + 1                                                          
          if( ISPBE(LL) > 3 .and. ISPRE(LL) == 0 )then
            L = LIJ(IPBE(LL),JPBE(LL))                                             
            if( LOPENBCDRY(L) ) CYCLE
            VOLINC = DELTD2*(FUHDYE(L) + OBCRESIDUALQ(IBC))          ! *** m3    Incremental volume for current time step
            OBCRESIDUALQ(IBC) = VOLINC + OBCRESIDUALQ1(IBC)          ! *** m3    Total volume at end          
          endif
        enddo   
        do LL = 1,NPBN                                                           
          IBC = IBC + 1                                                          
          if( ISPBN(LL) > 3 .and. ISPRN(LL) == 0 )then
            L = LIJ(IPBN(LL),JPBN(LL)) 
            if( LOPENBCDRY(L) ) CYCLE
            VOLINC = DELTD2*(FVHDXE(L) + OBCRESIDUALQ(IBC))          ! *** m3    Incremental volume for current time step
            OBCRESIDUALQ(IBC) = VOLINC + OBCRESIDUALQ1(IBC)          ! *** m3    Total volume at end          
          endif
        enddo                                                        
                                                                   
        SUMDAY = SUMDAY + RINC                                       ! *** s     Increment the timing counter
                                                                   
        do IBC = 1,NBCSOP                                            
          L = LOBCS(IBC)
          if( LOPENBCDRY(L) ) CYCLE
          OBCRESIDUALQ(IBC) = OBCRESIDUALQ(IBC)/SUMDAY               ! *** m3/s  Residual flow for current time step
          !OBCRESIDUALQ(IBC) = (OBCRESIDUALQ2(IBC) + 2.*OBCRESIDUALQ(IBC))/3.
        enddo

      endif
    endif
  endif           ! *** End of Residual Flow Update
  
  C1 = 0.5*G
  HDRY2  = 0.2*HDRY
  if( ISDRY > 0 )then
    HDRY5  = 5.0*HDRY
    HDRY10 = 10.*HDRY
  else
    HDRY5  = 0.0
    HDRY10 = 1E-12
  endif
  IBC = 0
  
  ! *******************************************************************************************************
  ! *** Set open boundary surface elevations
  TN = TIMESEC
  do M = 1,MTIDE
    TM = MOD(TN,TCP(M))
    TM = PI2*TM/TCP(M)
    CCCOS(M) = COS(TM)
    SSSIN(M) = SIN(TM)
  enddo

  ! *** SOUTH OPEN BOUNDARY
  do LL = 1,NPBS
    L = LPBS(LL)
    LN = LNC(L)
    IBC = IBC+ 1
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    if( LOPENBCDRY(L) ) CYCLE

    FP(L) = PSERT(NPSERS(LL)) + 0.5*PSERZDF(NPSERS(LL)) + PSERST(NPSERS(LL)) + 0.5*PSERZDS(NPSERS(LL))         ! *** m2/s2
    if( NPFORT >= 1 .and. NPSERS1(LL) > 0 )then
      TMPVAL = PSERT(NPSERS1(LL)) + 0.5*PSERZDF(NPSERS1(LL)) + PSERST(NPSERS1(LL)) + 0.5*PSERZDS(NPSERS1(LL))
      FP(L) = FP(L)+TPCOORDS(LL)*(TMPVAL-FP(L))
    endif
    do M = 1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
    enddo
    
    ELEV = FP(L)/G
    FP1G = ELEV - HDRY2                                              ! *** Prepare for boundary below minimum depths

    if( ISPBS(LL) == 1 .or. ISPBS(LL) == 2 )then
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERS(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation
                                                                     
      ! ***        s   m/s2   (-)     (-)     m                      
      CNT = 0.5*DELTD2*G   *HRVO(LN)*RCY(LN)*HVT(LN)                 ! *** m2/s
                                                                     
      ! ***    s        m/s        1/m                              
      TMP = DELTD2*SQRT(G*HVT(LN))*DYIV(LN)                          ! *** dimensionless
      CC(L) = CNT*(1.+TMP)/TMP                                       ! *** m2/s
      CN(L) = -CNT                                                   ! *** m2/s
      
      ! ***   m2/s      m2/s2          m/s       m3/s      1/m     1/m       (-)
      FP(L) = CNT*( 2.*FP(L) - SQRT(G*HVT(LN))*FVHDXE(LN)*DXIV(LN)/HVT(LN) )/TMP       ! *** m4/s3
      
      FP(L) = FP(L) + CNT*PSERAVG(NPSERS(LL),1)/TMP                  ! *** Add back offset
      
    elseif( ISPBS(LL) == 4 .or. ISPBS(LL) == 5 )then
      ! *** Outgoing wave (anti-reflection) type                     
      
      ! ***        s   m/s2   (-)     (-)     m                      
      CNT = 0.5*DELTD2*G   *HRVO(LN)*RCY(LN)*HVT(LN)                 ! *** m2/s
                                                                     
      ! ***    s        m/s        1/m                              
      TMP = DELTD2*SQRT(G*HVT(LN))*DYIV(LN)                          ! *** dimensionless
      CC(L) = CNT*(1.+TMP)/TMP                                       ! *** m2/s
      CN(L) = -CNT                                                   ! *** m2/s
      
      CST = FVHDXE(LN) - OBCRESIDUALQ(IBC)                           ! *** m3/s   Excess flow 
      if( FP1G < BELV(L) .or. FP1G < BELV(LN) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        FP(L) = DELTI*DXYP(L)*FP(L)                                  ! *** m4/s3
        CC(L) = DELTI*DXYP(L)                                        ! *** m2/s
        CS(L) = 0.
        CW(L) = 0.
        CE(L) = 0.
        CN(L) = 0.                                                   ! *** m2/s
        CS(LN) = 0.0
        LOPENBCDRY(L) = .TRUE.
        CYCLE                                                        ! *** Skip to next cell
      elseif( ELEV < (BELV(L) + HDRY10) )then
        ! *** Gradually reduce flows
        FACTOR = max((ELEV - BELV(L) - HDRY5)/HDRY10, 0.0)
        !CST = FACTOR*CST
      endif
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LN) = FP(LN) + CNT*FP(L)                                    ! *** m4/s3

      ! ***   m2/s         m/s       m3/s   1/m     1/m       (-)
      FPS = -CNT*( SQRT(G*HVT(LN))* CST   *DXIV(LN)/HVT(LN) )/TMP    ! *** m4/s3
      
      FP(L) = FPS + CNT*FP(L)/TMP                                    ! *** Add specified elevation
      
    else
      ! *** Inactivate BC'S when elevations drop below BELV + HDRY
      if( FP1G < BELV(L) .or. FP1G < BELV(LN) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CNT = 0.
        LOPENBCDRY(L) = .TRUE.
      else
        ! ***        s   m/s2   (-)     (-)      m
        CNT = 0.5*DELTD2*G    *HRVO(LN)*RCY(LN)*HVT(LN)              ! *** m2/s
      endif

      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LN) = FP(LN) + CNT*FP(L)                                    ! *** m4/s3
      CS(LN) = 0.

      ! ***    1/s   m2     m2/s2
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
      
    endif

  enddo

  ! *** WEST OPEN BOUNDARY
  do LL = 1,NPBW
    L = LPBW(LL)
    LE = LEC(L)
    IBC = IBC+ 1
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    if( LOPENBCDRY(L) ) CYCLE

    FP(L) = PSERT(NPSERW(LL)) + 0.5*PSERZDF(NPSERW(LL)) + PSERST(NPSERW(LL)) + 0.5*PSERZDS(NPSERW(LL))         ! *** m2/s2
    if( NPFORT >= 1 .and. NPSERW1(LL) > 0 )then
      TMPVAL = PSERT(NPSERW1(LL)) + 0.5*PSERZDF(NPSERW1(LL)) + PSERST(NPSERW1(LL)) + 0.5*PSERZDS(NPSERW1(LL))
      FP(L) = FP(L) + TPCOORDW(LL)*(TMPVAL-FP(L))
    endif
    do M = 1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
    enddo

    ELEV = FP(L)/G
    FP1G = ELEV - HDRY2                                           ! *** Prepare for boundary below minimum depths

    if( ISPBW(LL) == 1 .or. ISPBW(LL) == 2 )then
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERW(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation

      ! ***        s   m/s2   (-)  (-)     m
      CET = 0.5*DELTD2*G*HRUO(LE)*RCX(LE)*HUT(LE)                    ! *** m2/s

      ! ***    s        m/s        1/m 
      TMP = DELTD2*SQRT(G*HUT(LE))*DXIU(LE)                          ! *** dimensionless
      CC(L) = CET*(1.+TMP)/TMP                                       ! *** m2/s
      CE(L) = -CET                                                   ! *** m2/s

      ! ***   m2/s    m2/s2          m/s         m3/s     1/m     1/m     (-)
      FP(L) = CET*(2.*FP(L) - SQRT(G*HUT(LE))*FUHDYE(LE)*DYIU(LE)/HUT(LE))/TMP    ! *** m4/s3

      FP(L) = FP(L) + CET*PSERAVG(NPSERW(LL),1)/TMP                  ! *** Add back offset

    elseif( ISPBW(LL) == 4 .or. ISPBW(LL) == 5 )then                 
      ! *** Outgoing wave (anti-reflection) type                     
                                                                     
      ! ***        s   m/s2 (-)   (-)      m
      CET = 0.5*DELTD2*G*HRUO(LE)*RCX(LE)*HUT(LE)                    ! *** m2/s

      ! ***    s        m/s        1/m 
      TMP = DELTD2*SQRT(G*HUT(LE))*DXIU(LE)                          ! *** dimensionless
      CC(L) = CET*(1.+TMP)/TMP                                       ! *** m2/s
      CE(L) = -CET
      
      CWT = FUHDYE(LE) - OBCRESIDUALQ(IBC)                           ! *** m3/s   Excess flow 
      if( FP1G < BELV(L) .or. FP1G < BELV(LE) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        FP(L) = DELTI*DXYP(L)*FP(L)                                  ! *** m4/s3
        CC(L) = DELTI*DXYP(L)                                        ! *** m2/s
        CS(L) = 0.
        CW(L) = 0.
        CE(L) = 0.
        CN(L) = 0.                                                   ! *** m2/s
        CW(LE) = 0.0
        LOPENBCDRY(L) = .TRUE.
        CYCLE                                                        ! *** Skip to next cell
      elseif( ELEV < (BELV(L) + HDRY10) )then
        ! *** Gradually reduce flows
        FACTOR = max((ELEV - BELV(L) - HDRY5)/HDRY10, 0.0)
        !CWT = FACTOR*CWT
      endif
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LE) = FP(LE) + CET*FP(L)                                    ! *** m4/s3

      ! ***   m2/s         m/s       m3/s       1/m    1/m     (-)
      FPW = -CET*(SQRT(G*HUT(LE))* CWT*      DYIU(L)/HUT(LE))/TMP    ! *** m4/s3

      FP(L) = FPW + CET*FP(L)/TMP                                    ! *** Add specified elevation

    else
      ! *** Inactivate BC'S when elevations drop below BELV + HDRY
      if( FP1G < BELV(L) .or. FP1G < BELV(LE) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CET = 0.
        LOPENBCDRY(L) = .TRUE.
      else
        CET = 0.5*DELTD2*G*HRUO(LE)*RCX(LE)*HUT(LE)
      endif
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LE) = FP(LE) + CET*FP(L)                                    ! *** m4/s3
      CW(LE) = 0.

      ! ***    1/s   m2     m2/s2
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
    endif
  enddo

  ! *** EAST OPEN BOUNDARY
  do LL = 1,NPBE
    L = LPBE(LL)
    LW = LWC(L)
    IBC = IBC+ 1
    CC(L) = DELTI*DXYP(L)                                            ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    
    if( LOPENBCDRY(L) ) CYCLE
    
    FP(L) = PSERT(NPSERE(LL)) + 0.5*PSERZDF(NPSERE(LL)) + PSERST(NPSERE(LL)) + 0.5*PSERZDS(NPSERE(LL))         ! *** m2/s2
    if( NPFORT >= 1 .and. NPSERE1(LL) > 0 )then
      TMPVAL = PSERT(NPSERE1(LL)) +0.5*PSERZDF(NPSERE1(LL)) + PSERST(NPSERE1(LL)) + 0.5*PSERZDS(NPSERE1(LL))
      FP(L) = FP(L)+TPCOORDE(LL)*(TMPVAL-FP(L))
    endif
    
    do M = 1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
    enddo

    ELEV = FP(L)/G
    FP1G = ELEV - HDRY2                                           ! *** Prepare for boundary below minimum depths

    if( ISPBE(LL) == 1 .or. ISPBE(LL) == 2 )then
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERE(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation
                                                                     
      ! ***        s    m/s2   (-)     (-)     m                
      CWT = 0.5*DELTD2* G*     HRUO(L)*RCX(L)*HUT(L)                 ! *** m2/s
                                                                     
      ! ***    s        m/s        1/m                             
      TMP = DELTD2*SQRT(G*HUT(L))*DXIU(L)                            ! *** dimensionless
      CC(L) = CWT*(1.+TMP)/TMP                                       ! *** m2/s
      CW(L) = -CWT                                                   ! *** m2/s
      
      ! ***   m2/s    m2/s2          m/s       m3/s     1/m    1/m     (-)
      FP(L) = CWT*(2.*FP(L) + SQRT(G*HUT(L))*FUHDYE(L)*DYIU(L)/HUT(L))/TMP       ! *** m4/s3

      FP(L) = FP(L) + CWT*PSERAVG(NPSERE(LL),1)/TMP                  ! *** Add back offset
                                                                     
    elseif( ISPBE(LL) == 4 .or. ISPBE(LL) == 5 )then                 
      ! *** Outgoing wave (anti-reflection) type       
      
      ! ***        s    m/s2    (-)     (-)    m                
      CWT = 0.5*DELTD2* G*     HRUO(L)*RCX(L)*HUT(L)                 ! *** m2/s
                                                                     
      ! ***    s        m/s        1/m                             
      TMP = DELTD2*SQRT(G*HUT(L))*DXIU(L)                            ! *** dimensionless
      CC(L) = CWT*(1.+TMP)/TMP                                       ! *** m2/s
      CW(L) = -CWT                                                   ! *** m2/s
      
      CET = FUHDYE(L) - OBCRESIDUALQ(IBC)                            ! *** m3/s   Excess flow 
      if( FP1G < BELV(L) .or. FP1G < BELV(LW) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        FP(L) = DELTI*DXYP(L)*FP(L)                                  ! *** m4/s3
        CC(L) = DELTI*DXYP(L)                                        ! *** m2/s
        CS(L) = 0.
        CW(L) = 0.
        CE(L) = 0.
        CN(L) = 0.                                                   ! *** m2/s
        CE(LW) = 0.0
        LOPENBCDRY(L) = .TRUE.
        CYCLE                                                        ! *** Skip to next cell
      elseif( ELEV < (BELV(L) + HDRY10) )then
        ! *** Gradually reduce flows
        FACTOR = max((ELEV - BELV(L) - HDRY5)/HDRY10, 0.0)
        !CET = FACTOR*CET
      endif   
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2                           
      FP(LW) = FP(LW) + CWT*FP(L)                                    ! *** m4/s3
      
      ! ***  m2/s         m/s       m3/s     1/m    1/m     (-)
      FPE = CWT*( SQRT(G*HUT(L))* CET*     DYIU(L)/HUT(L))/TMP       ! *** m4/s3

      FP(L) = FPE + CWT*FP(L)/TMP                                    ! *** Add specified elevation
      
    else
      ! *** Inactivate BC'S when elevations drop below BELV + HDRY
      if( FP1G < BELV(L) .or. FP1G < BELV(LW) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CWT = 0.                                                     
        LOPENBCDRY(L) = .TRUE.                                       
      else                                                           
        ! ***        s    m/s2    (-)     (-)    m              
        CWT = 0.5*DELTD2* G*     HRUO(L)*RCX(L)*HUT(L)               ! *** m2/s
      endif
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2                           
      FP(LW) = FP(LW) + CWT*FP(L)                                    ! *** m4/s3
      CE(LW) = 0.0
                                                                     
      ! ***    1/s   m2    m2/s2                                     
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3

    endif
    
  enddo

  ! *** NORTH OPEN BOUNDARY
  do LL = 1,NPBN
    L = LPBN(LL)
    LS = LSC(L)
    IBC = IBC+ 1
    CC(L) = DELTI*DXYP(L)   ! *** m2/s
    CS(L) = 0.
    CW(L) = 0.
    CE(L) = 0.
    CN(L) = 0.
    if( LOPENBCDRY(L) ) CYCLE
    
    FP(L) = PSERT(NPSERN(LL)) + 0.5*PSERZDF(NPSERN(LL)) + PSERST(NPSERN(LL)) + 0.5*PSERZDS(NPSERN(LL))         ! *** m2/s2
    if( NPFORT >= 1 .and. NPSERN1(LL) > 0 )then
      TMPVAL = PSERT(NPSERN1(LL)) + 0.5*PSERZDF(NPSERN1(LL)) + PSERST(NPSERN1(LL)) + 0.5*PSERZDS(NPSERN1(LL))
      FP(L) = FP(L)+TPCOORDN(LL)*(TMPVAL-FP(L))
    endif
    do M = 1,MTIDE
      TC = CCCOS(M)
      TS = SSSIN(M)
      FP(L) = FP(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
    enddo

    ELEV = FP(L)/G
    FP1G = ELEV - HDRY2                                           ! *** Prepare for boundary below minimum depths

    if( ISPBN(LL) == 1 .or. ISPBN(LL) == 2 )then
      ! *** RADIATION BC TYPES
      
      FP(L) = FP(L) - PSERAVG(NPSERN(LL),1)                          ! *** Set radiation calcs in terms of displacement not elevation
                                                                     
      ! ***        s   m/s2   (-)     (-)     m                      
      CST = 0.5*DELTD2*G   *HRVO(L)*RCY(L)*HVT(L)                    ! *** m2/s

      ! ***    s        m/s        1/m                              
      TMP = DELTD2*SQRT(G*HVT(L))*DYIV(L)                            ! *** dimensionless
      CC(L) = CST*(1.+TMP)/TMP                                       ! *** m2/s
      CS(L) = -CST                                                   ! *** m2/s
      
      ! ***   m2/s      m2/s2         m/s       m3/s    1/m     1/m     (-)
      FP(L) = CST*( 2.*FP(L) + SQRT(G*HVT(L))*FVHDXE(L)*DXIV(L)/HVT(L))/TMP     ! *** m4/s3
      
      FP(L) = FP(L) + CNT*PSERAVG(NPSERN(LL),1)/TMP                  ! *** Add back offset
      
    elseif( ISPBN(LL) == 4 .or. ISPBN(LL) == 5 )then
      ! *** Outgoing wave (anti-reflection) type                     
      
      ! ***        s   m/s2   (-)     (-)     m                      
      CST = 0.5*DELTD2*G   *HRVO(L)*RCY(L)*HVT(L)                    ! *** m2/s

      ! ***    s        m/s        1/m                              
      TMP = DELTD2*SQRT(G*HVT(L))*DYIV(L)                            ! *** dimensionless
      CC(L) = CST*(1.+TMP)/TMP                                       ! *** m2/s
      CS(L) = -CST                                                   ! *** m2/s
      
      CNT = FVHDXE(L) - OBCRESIDUALQ(IBC)                            ! *** m3/s   Excess flow 
      if( FP1G < BELV(L) .or. FP1G < BELV(LS) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        FP(L) = DELTI*DXYP(L)*FP(L)                                  ! *** m4/s3
        CC(L) = DELTI*DXYP(L)                                        ! *** m2/s
        CS(L) = 0.
        CW(L) = 0.
        CE(L) = 0.
        CN(L) = 0.                                                   ! *** m2/s
        CN(LS) = 0.0
        LOPENBCDRY(L) = .TRUE.
        CYCLE                                                        ! *** Skip to next cell
      elseif( ELEV < (BELV(L) + HDRY10) )then
        ! *** Gradually reduce flows
        FACTOR = max((ELEV - BELV(L) - HDRY5)/HDRY10, 0.0)
        !CNT = FACTOR*CNT
      endif
      
      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LS) = FP(LS) + CST*FP(L)                                    ! *** m4/s3
      
      ! *** m2/s          m/s     m3/s    1/m     1/m    (-)
      FPN = CST*( SQRT(G*HVT(L)) *CNT   *DXIV(L)/HVT(L))/TMP         ! *** m4/s3
       
      FP(L) = FPN + CST*FP(L)/TMP                                    ! *** Add specified elevation
      
    else
      ! *** Inactivate BC'S when elevations drop below BELV + HDRY
      if( FP1G < BELV(L) .or. FP1G < BELV(LS) )then
        FP(L) = (BELV(L) + HDRY2)*G                                  ! *** m2/s2
        CST = 0.
        LOPENBCDRY(L) = .TRUE.
      else
        ! ***        s   m/s2    (-)     (-)     m
        CST = 0.5*DELTD2*G     *HRVO(L)*RCY(L)*HVT(L)                ! *** m2/s
      endif

      ! *** Setup adjacent cell
      ! ***    m4/s3    m2/s m2/s2
      FP(LS) = FP(LS) + CST*FP(L)                                    ! *** m4/s3
      CN(LS) = 0.

      ! ***    1/s   m2    m2/s2
      FP(L) = DELTI*DXYP(L)*FP(L)                                    ! *** m4/s3
      
    endif
    
  enddo

  return
END

