! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE BEDLOAD(NX,NS)

  !****************************************************************************
  ! *** SUBROUTINE CALSND CALCULATES NON-COHESIVE SEDIMENT SETTLING,                                                      
  ! *** DEPOSITION AND RESUSPENSION AND IS CALLED FROM SSEDTOX                                                             
  ! *** NOT USED BY SEDZLJ
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !   2011-03       Paul M. Craig      Converted to F90, added OMP

  use GLOBAL
  use Variables_Propwash
  
  implicit none                                                                                                          
  
  integer, intent(IN) :: NX,NS
  integer :: L, NSB, LUTMP, LDTMP, LE, LS, LN, LW, ND, LF, LL, LP
  real    :: SLOPE, ASNDFBL, SNDBTMP, SNDFBLM,  FLUXFAC
  real    :: UCELLCTRM, VCELLCTRM, SHIELDS, BDLDTMPB, CSHIELDSC                                                        
  real    :: BDLDTMPA, CSHIELDS, TMPVAL, BDLDTMPP, BDLDTMP, XOUT, YOUT                                                           
  
  real, external :: FSEDMODE
  real, external :: FSBDLD
  
  ! *** ZERO BOUNDARY BEDLOAD FLUXES
  if( NSBDLDBC > 0 )then
    do NSB = 1,NSBDLDBC
      LUTMP = LSBLBCU(NSB)
      QSBDLDOT(LUTMP,NX) = 0.       ! *** BC OUTFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
      QSBDLDIN(LUTMP,NX) = 0.       ! *** BC INFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
      
      LDTMP = LSBLBCD(NSB)
      if( LDTMP > 0 )then
        QSBDLDOT(LDTMP,NX) = 0.     ! *** BC OUTFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
        QSBDLDIN(LDTMP,NX) = 0.     ! *** BC INFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
      endif
    enddo
  endif
  
  ! *** BED LOAD TRANSPORT HORIZONTAL LOOP                                                                                

  !$OMP PARALLEL DEFAULT(SHARED)
  
  !$OMP DO PRIVATE(ND,LF,LL,L)
  do ND = 1,NDM  
    LF = 2+(ND-1)*LDM  
    LL = MIN(LF+LDM-1,LA)

    ! *** ZERO BED LOAD TRANSPORTS                                                                                    
    do L = LF,LL
      QSBDLDP(L) = 0.         ! *** CELL CENTER BED LOAD TRANSPORT RATE (G/M/S)
      QSBDLDX(L,NX) = 0.      ! *** U FACE SND FLUX DUE TO BEDLOAD (G/S)
      QSBDLDY(L,NX) = 0.      ! *** V FACE SND FLUX DUE TO BEDLOAD (G/S) 
      SNDFBL(L,NX) = 0.       ! *** BED/WATER INTERFACE SND FLUX DUE TO BEDLOAD (G/M2/S)
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  !*********************************************************************************
  !*** COMPUTE THE CELL CENTERED BEDLOAD TRANSPORT RATES USING THE SPECIFIED OPTION
      
  if( ISBDLD(NS) == 0 )then
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING GENERIC BED LOAD EQUATION

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS,LN) &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA,BDLDTMPB) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        if( LMASKDRY(L) )then
          CSHIELDS = TCSHIELDS(NS)
          if( ISEDEFF == 2 )then
            TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
            CSHIELDS = TMPVAL*CSHIELDS
          endif
          BDLDTMPP = SBDLDP(NX)
          BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
          if( BDLDTMPP > 0.0 )then
            FACBEDL(L) = FSEDMODE(WSETA(L,0,NS), USTAR(L), USTARSND(L), RSNDM(NX), ISNDM1(NX), ISNDM2(NX),1)
            SHIELDS = TAUBSND(L)/GPDIASED
            if( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDS )then
              if( SBDLDA(NX) > 0.0 )then
                BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
              else
                BDLDTMPA = 1.0
              endif
              if( SBDLDB(NX) > 0.0 )then
                BDLDTMPB = (SBDLDG3(NX)*SQRT(SHIELDS)-SBDLDG4(NX)*SQRT(CSHIELDS))**SBDLDB(NX)
              else
                BDLDTMPB = 1.0
              endif
              QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA*BDLDTMPB     ! *** (G/M/S)
              if( ISEDEFF == 1 )then
                TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
                QSBDLDP(L) = TMPVAL*QSBDLDP(L)
              endif
            endif
          endif
        endif
      enddo
    
    enddo   ! *** END OF DOMAIN LOOP FOR GENERIC BED LOAD EQUATION
    !$OMP END DO

  elseif( ISBDLD(NS) == 1 )then
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING VAN RIJN BED LOAD EQUATION

    !$OMP DO PRIVATE(ND,LF,LL,LP,L)   &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,CSHIELDSC,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        CSHIELDS = TCSHIELDS(NS)
        if( ISEDEFF == 2 )then
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        endif
        if( ISNDAL == 1 )then
          TMPVAL = LOG10(19.*DIASED/SEDDIA50(L,KBT(L)))
          TMPVAL = 1.66667/(TMPVAL**2)
          CSHIELDSC = CSHIELDS50(L)*TMPVAL
        else
          CSHIELDSC = TCSHIELDS(NS)
        endif
        if( ISEDEFF == 2 )then
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDSC = TMPVAL*CSHIELDS
        endif
        BDLDTMPP = FSBDLD(DIASED, GPDIASED, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
        if( ISNDAL == 1 )then
          BDLDTMPP = ((DIASED/SEDDIA50(L,KBT(L)))**0.3)*BDLDTMPP
        endif
        ! ***            M/S       M     G/M3
        BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM                                                              ! *** G/M/S
        FACBEDL(L) = FSEDMODE(WSETA(L,0,NS), USTAR(L), USTARSND(L), RSNDM(NX), ISNDM1(NX), ISNDM2(NX),1)     ! *** Dimensionless
        SHIELDS = TAUBSND(L)/GPDIASED
        if( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDSC )then
          BDLDTMPA = (SBDLDG1(NX)*SHIELDS -SBDLDG2(NX)*CSHIELDSC)**SBDLDA(NX)
          ! ***           (-)       (-)                G/M/S    (-)       (-) 
          QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA                               ! *** G/M/S
          if( ISEDEFF == 1 )then
            TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
            QSBDLDP(L) = TMPVAL*QSBDLDP(L)
          endif
        endif
      enddo
    enddo   ! *** END OF DOMAIN LOOP FOR VAN RIJN BED LOAD EQUATION
    !$OMP END DO

  elseif( ISBDLD(NS) == 2 )then
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING ENGELUND-HANSEN

    !$OMP DO PRIVATE(ND,LF,LL,LP,L)   &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      CSHIELDS = 0.
      if( IBLTAUC(NS) == 1 ) CSHIELDS = TCSHIELDS(NS)
      do LP = LF,LL
        L = LSED(LP)  

        if( IBLTAUC(NS) == 2 ) CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
        if( IBLTAUC(NS) == 3 ) CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
        if( ISEDEFF == 2 )then
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        endif
        SHIELDS = TAUBSND(L)/GPDIASED
        if( SHIELDS >= CSHIELDS )then
          BDLDTMPP = FSBDLD(DIASED, GPDIASED, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
          if( HGDH(L) > 0.0) BDLDTMPP = BDLDTMPP/(HGDH(L)**0.333)
          BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
          FACBEDL(L) = FSEDMODE(WSETA(L,0,NS), USTAR(L), USTARSND(L), RSNDM(NX), ISNDM1(NX), ISNDM2(NX),1)
          BDLDTMPA = (SBDLDG1(NX)*SHIELDS)**SBDLDA(NX)
          QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA
          if( ISEDEFF == 1 )then
            TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
            QSBDLDP(L) = TMPVAL*QSBDLDP(L)
          endif
        endif
      enddo
    enddo   ! *** END OF DOMAIN LOOP FOR ENGELUND-HANSEN SECTION
    !$OMP END DO

  elseif( ISBDLD(NS) == 3 )then
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING WU, WANG, AND JIA

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS,LN) &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        CSHIELDS = TCSHIELDS(NS)
        BDLDTMPP = FSBDLD(DIASED, GPDIASED, SEDDIA50(L,KBT(L)), HP(L), PEXP(L,NX), PHID(L,NX), CSHIELDS, SBDLDP(NX), ISBDLD(NS))
        CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
        if( ISEDEFF == 2 )then
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        endif
        BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
        FACBEDL(L) = FSEDMODE(WSETA(L,0,NS), USTAR(L), USTARSND(L), RSNDM(NX), ISNDM1(NX), ISNDM2(NX), 1)
        SHIELDS = TAUBSND(L)/GPDIASED
        if( SBDLDG1(NX)*SHIELDS > SBDLDG2(NX)*CSHIELDS )then
          BDLDTMPA = (SBDLDG1(NX)*SHIELDS - SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
          QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA
          if( ISEDEFF == 1 )then
            TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
            QSBDLDP(L) = TMPVAL*QSBDLDP(L)
          endif
        endif
      enddo
    enddo   ! *** END OF DOMAIN LOOP FOR WU, WANG, AND JIA SECTION
    !$OMP END DO
    
  endif  ! *** END OF CELL CENTERED BEDLOAD TRANSPORT CALCULATIONS
  
  ! ********************************************************************************
  ! *** INCORPORATE PROPWASH BEDLOAD
  if( ISPROPWASH > 0 )then
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)
        if( PROP_BLD(L,NS) > 0.0 )then
          !print '(i8,2i5,f10.4)', niter, l, ns, PROP_BLD(L,NS)
          QSBDLDP(L) = PROP_BLD(L,NS)              ! *** Overwrite ambient current based bedload (g/m/s)
        endif
      enddo
    enddo
    !$OMP END DO
  
  endif
  
  ! ********************************************************************************
  ! *** CALCULATE CELL FACE TRANSPORT RATES (G/M/S)
  if( ISBLFUC == 0 )then
    ! *** CALCULATE CELL FACE TRANSPORT RATES BY DOWN WIND PROJECTION

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LW,LS)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        LW = LWC(L)
        LS = LSC(L)
        
        if( UCELLCTR(LW) > 0.0 )then
          QSBDLDX(L,NX) = SUB(L)*QSBDLDP(LW)*UCELLCTR(LW)
        elseif( UCELLCTR(L) < 0.0 )then
          QSBDLDX(L,NX) = SUB(L)*QSBDLDP(L)*UCELLCTR(L)
        endif
        
        if( VCELLCTR(LS) > 0.0 )then
          QSBDLDY(L,NX) = SVB(L)*QSBDLDP(LS)*VCELLCTR(LS)
        elseif( VCELLCTR(L) < 0.0 )then
          QSBDLDY(L,NX) = SVB(L)*QSBDLDP(L)*VCELLCTR(L)
        endif
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  elseif( ISBLFUC == 1 )then
    ! *** CALCULATE CELL FACE TRANSPORT RATES BY DOWN WIND PROJECTION WITH CORNER EFFECTS CORRECTION                                                                                    

    !$OMP DO PRIVATE(ND,LF,LL,L,LE,LN)       &
    !$OMP    PRIVATE(UCELLCTRM,VCELLCTRM)
    do ND = 1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = MIN(LF+LDM-1,LA)
      
      do L = LF,LL
        if( LMASKDRY(L) )then
          LE = LEC(L)
          LN = LNC(L)
          if( UCELLCTR(L) >= 0.0 .and. VCELLCTR(L) >= 0.0 )then
            UCELLCTRM = SUB(LE)*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(LN)*ABS(VCELLCTR(L))
            QSBDLDX(LE,NX) = SUB(LE)*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(LN,NX) = SVB(LN)*QSBDLDP(L)*VCELLCTR(L)
            if( UCELLCTRM < 1.0E-9 )then
              QSBDLDY(LN ,NX) = SVB(LN )*SIGN(QSBDLDP(L),VCELLCTR(L))
            endif
            if( VCELLCTRM < 1.0E-9 )then
              QSBDLDX(LE,NX) = SUB(LE)*SIGN(QSBDLDP(L),UCELLCTR(L))
            endif
          endif
          if( UCELLCTR(L) >= 0.0 .and. VCELLCTR(L) < 0.0 )then
            UCELLCTRM = SUB(LE)*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(L )*ABS(VCELLCTR(L))
            QSBDLDX(LE,NX) = SUB(LE)*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(L ,NX) = SVB(L )*QSBDLDP(L)*VCELLCTR(L)
            if( UCELLCTRM < 1.0E-9 )then
              QSBDLDY(L  ,NX) = SVB(L  )*SIGN(QSBDLDP(L),VCELLCTR(L))
            endif
            if( VCELLCTRM < 1.0E-9 )then
              QSBDLDX(LE,NX) = SUB(LE)*SIGN(QSBDLDP(L),UCELLCTR(L))
            endif
          endif
          if( UCELLCTR(L) < 0.0 .and. VCELLCTR(L) >= 0.0 )then
            UCELLCTRM = SUB(L  )*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(LN )*ABS(VCELLCTR(L))
            QSBDLDX(L  ,NX) = SUB(L  )*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(LN ,NX) = SVB(LN )*QSBDLDP(L)*VCELLCTR(L)
            if( UCELLCTRM < 1.0E-9 )then
              QSBDLDY(LN ,NX) = SVB(LN )*SIGN(QSBDLDP(L),VCELLCTR(L))
            endif
            if( VCELLCTRM < 1.0E-9 )then
              QSBDLDX(L  ,NX) = SUB(L  )*SIGN(QSBDLDP(L),UCELLCTR(L))
            endif
          endif
          if( UCELLCTR(L) < 0.0 .and. VCELLCTR(L) < 0.0 )then
            UCELLCTRM = SUB(L)*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(L)*ABS(VCELLCTR(L))
            QSBDLDX(L,NX) = SUB(L)*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(L,NX) = SVB(L)*QSBDLDP(L)*VCELLCTR(L)
            if( UCELLCTRM < 1.0E-9 )then
              QSBDLDY(L,NX) = SVB(L)*SIGN(QSBDLDP(L),VCELLCTR(L))
            endif
            if( VCELLCTRM < 1.0E-9 )then
              QSBDLDX(L,NX) = SUB(L)*SIGN(QSBDLDP(L),UCELLCTR(L))
            endif
          endif
        endif
      enddo

    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  elseif( ISBLFUC == 2 )then
    ! *** CALCULATE CELL FACE TRANSPORT RATES BY AVERAGING VECTOR COMPONENTS FROM CELL CENTERS TO FACES                                                                      
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS,LW)   
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        LS = LSC(L)
        LW = LWC(L)
        QSBDLDX(L,NX) = 0.5*SUB(L)*(QSBDLDP(L)*UCELLCTR(L) + QSBDLDP(LW)*UCELLCTR(LW))
        QSBDLDY(L,NX) = 0.5*SVB(L)*(QSBDLDP(L)*VCELLCTR(L) + QSBDLDP(LS)*UCELLCTR(LS))
      enddo

    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  endif    ! *** END OF CELL FACE TRANSPORT RATES SECTION
    
    
  ! ********************************************************************************
  ! *** CONVERT TRANSPORT VECTORS TO FACE VECTORS  (G/S)                                                                        
  !$OMP DO PRIVATE(ND,LF,LL,LP,L)   
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = MIN(LF+LDMSED-1,LASED)

    do LP = LF,LL
      L = LSED(LP)
      QSBDLDX(L,NX) = SUB(L)*DYU(L)*QSBDLDX(L,NX)
      QSBDLDY(L,NX) = SVB(L)*DXV(L)*QSBDLDY(L,NX)
    enddo

  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! ********************************************************************************
  ! *** ELIMINATE BEDLOAD TRANSPORT UP ADVERSE SLOPES IN DIRECTION OF FLOW                                                
  if( BLBSNT > 0.0 )then

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,SLOPE)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        if( QSBDLDX(L,NX) > 0.0 )then
          SLOPE = (BELV(L)-BELV(LWC(L)))*DXIU(L)
          if( SLOPE > BLBSNT ) QSBDLDX(L,NX) = 0.0
        endif
        if( QSBDLDX(L,NX) < 0.0 )then
          SLOPE = (BELV(LWC(L))-BELV(L))*DXIU(L)
          if( SLOPE > BLBSNT ) QSBDLDX(L,NX) = 0.0
        endif
        if( QSBDLDY(L,NX) > 0.0 )then
          SLOPE = (BELV(L)-BELV(LSC(L)))*DYIV(L)
          if( SLOPE > BLBSNT ) QSBDLDY(L,NX) = 0.0
        endif
        if( QSBDLDY(L,NX) < 0.0 )then
          SLOPE = (BELV(LSC(L))-BELV(L))*DYIV(L)
          if( SLOPE > BLBSNT ) QSBDLDY(L,NX) = 0.0
        endif
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  endif
  
  ! *** Provide a rampup time if SEDSTART > TBEGIN
  if( TIMEDAY < SEDSTART + 1.0 )then   ! DELME - allow user defined rampup period
    FLUXFAC = TIMEDAY - SEDSTART
    if( FLUXFAC > 0.0 )then
      FLUXFAC = FLUXFAC / 1.0    ! DELME - 1.0 IS RAMPUP TIME
      
      !$OMP DO PRIVATE(ND,LF,LL,LP,L)   
      do ND = 1,NDM  
        LF = (ND-1)*LDMSED+1  
        LL = MIN(LF+LDMSED-1,LASED)
  
        do LP = LF,LL
          L = LSED(LP)
          QSBDLDX(L,NX) = FLUXFAC*QSBDLDX(L,NX)
          QSBDLDY(L,NX) = FLUXFAC*QSBDLDY(L,NX)
          QSBDLDP(L)    = FLUXFAC*QSBDLDP(L)   
        enddo
  
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    endif
  endif
  
  !$OMP END PARALLEL

  ! ********************************************************************************
  ! *** LIMIT OUTGOING FLUXES IN EACH CELL                                                                                
  ! *** PMC - DELETED FOR 8.4.4.  SUFFICIENT SEDIMENT CHECK IN CALSND
  
  return

END

