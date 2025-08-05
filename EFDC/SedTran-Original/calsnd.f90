! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALSND  

  !**********************************************************************!
  ! *** SUBROUTINE CALSND CALCULATES NONCOHESIVER SEDIMENT SETTLING,  
  ! *** DEPOSITION AND RESUSPENSION AND IS CALLED FROM SSEDTOX  
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2019-01           PAUL M. CRAIG    ADDED HARD BOTTOM BYPASS
  ! 2011-03-02        PAUL M. CRAIG    REWRITTEN TO CODE F90
  !                                     REMOVED KC DEPENDENT DUPLCIATE CODE
  !                                     ADDED OMP
  ! 2010-XX-XX        SCOTT JAMES      ADDED MHK
  !
  !**********************************************************************!
  
  use GLOBAL  
  use Allocate_Initialize
  use MPI
  use Variables_MPI
  use Variables_Propwash

  implicit none

  integer :: IERR, NX, NS, K, L, LP, ISGPFLAG, ISEHFLAG, KTOP, NXX, NSS
  integer :: IFLAG, LE, LW, LN, LS, ND, LF, LL, L1, NSB, LUTMP, LDTMP
  real    :: TIME, GRADSED, SIGP, CRNUM, DUM1, DUM3, DUM4
  real    :: ZEQMAX, CSHIELDS, TMPVAL
  real    :: SHIELDS, TOP, BOT, WSFAC, WESE, FLUXFAC
  real    :: WSETMP, WVEL, CLEFT, CRIGHT, SNDBTMP, SEDAVG
  real, external :: FSEDMODE, CSNDZEQ, CSNDSET, CSNDEQC
  
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS, TWAIT                     ! MODEL TIMING TEMPORARY VARIABLES

  real,save,allocatable,dimension(:) :: SEDPHI

  ! *** FIRST CALL ALLOCATIONS AND ASSIGNMENTS
  if( .not. allocated(SEDPHI) )then
    call AllocateDSI(SEDPHI, NSTM2, 0.0)      

    ! *** GP - CONVERT SEDIMENT DIAMETERS IN M TO MM AND SET PHI SIZE  
    do NX = 1,NSND  
      NS = NSED + NX  
      SEDPHI(NS) = -LOG(1000.*SEDDIA(NS))/LOG(2.)  
    enddo  

    ! *** IF ISNDAL == 0 THEN NO ARMORING, CONSTANTS FOR ENTIRE SIMULATION
    do K = 1,KB  
      do L = 2,LA  
        SIGPHI(L,K) = 0.  
      enddo  
    enddo  

    do NX = 1,NSND  
      do L = 2,LA  
        PEXP(L,NX) = 1.  
        PHID(L,NX) = 1.  
      enddo  
    enddo  

    ! *** SET MAXIMUM NONCOHESIVE SEDIMENT DIAMETER  
    SNDDMX = 0.0
    do NX = 1,NSND  
      NS = NSED + NX  
      SNDDMX = max(SNDDMX,SEDDIA(NS))  
    enddo  

    if( ISNDVW == 0 )then 
      ! *** CONSTANT (ONLY ASSIGN AT START OF RUN)
      do NX = 1,NSND  
        NS = NSED + NX
        do K = 0,KS  
          ! *** USING 2,LA TO ENSURE ASSIGNMENT EVEN IF DRY CELLS EXIST
          do L = 2,LA                       
            WSETA(L,K,NS) = WSEDO(NS)  
          enddo  
        enddo  
      enddo
    endif
  endif
  
  !**********************************************************************!  
  TIME = TIMESEC/TCON  

  !**********************************************************************!  
  ! *** SET TIME/SPATIALLY VARYING NONCOHESIVE ARMORING parameterS  
  if( ISNDAL >= 1 )then  
    ! *** 1 ACTIVATE NON-COHESIVE ARMORING EFFECTS (GARCIA & PARKER)
    ! *** 2 SAME AS 1 WITH ACTIVE-PARENT LAYER FORMULATION
    ISGPFLAG = 0  
    ISEHFLAG = 0  
    do NX = 1,NSND  
      NS = NSED + NX 
      if( ISNDEQ(NS) == 1 ) ISGPFLAG = 1  
      if( ISBDLD(NS) >= 2 ) ISEHFLAG = 1  
    enddo  

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,KTOP,NXX,NSS)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      ! *** SET SIGPHI FOR GARCIA AND PARKER (1991) EQS 46 AND 47  
      if( ISGPFLAG == 1 )then  

        ! *** GP - SET MEAN PHI FOR TOP LAYER OF BED  
        do LP = LF,LL
          L = LSED(LP)
          TVAR3W(L) = 0.  
          TVAR3E(L) = 0.  
          SIGPHI(L,KBT(L)) = 0.  
        enddo  

        do NX = 1,NSND  
          NS = NSED + NX  
          do LP = LF,LL
            L = LSED(LP)
            KTOP = KBT(L)  
            TVAR3W(L) = TVAR3W(L) + SEDPHI(NS)*VFRBED(L,KTOP,NS)  
            TVAR3E(L) = TVAR3E(L) + VFRBED(L,KTOP,NS)  
          enddo  
        enddo  

        do LP = LF,LL
          L = LSED(LP)
          if( TVAR3E(L) <= 0. ) TVAR3E(L) = 1.  
          TVAR3W(L) = TVAR3W(L)/TVAR3E(L)  
        enddo  

        do NX = 1,NSND  
          NS = NSED + NX  
          do LP = LF,LL
            L = LSED(LP)
            KTOP = KBT(L)  
            SIGPHI(L,KTOP) = SIGPHI(L,KTOP) + ((SEDPHI(NS)-TVAR3W(L))**2)*VFRBED(L,KTOP,NS)/TVAR3E(L)
          enddo  
        enddo  

        do LP = LF,LL
          L = LSED(LP)
          KTOP = KBT(L)
          if( SIGPHI(L,KTOP) < 0. )then  
            SIGPHI(L,KTOP) = -SQRT(ABS(SIGPHI(L,KTOP)))   ! *** PMC TEMP PATCH  
          else
            SIGPHI(L,KTOP) = SQRT(SIGPHI(L,KTOP))  
          endif
        enddo  

      endif    ! *** END CALCULATION OF SIGPHI FOR GARCIA AND PARKER (1991)  

      ! *** SET EXPOSURE AND HIDING FUNCTIONS FOR ENGULAND-HANSEN AND WU,WANG,  
      if( ISEHFLAG == 1 .or. ISGPFLAG == 1 )then  

        do NX = 1,NSND  
          do LP = LF,LL
            L = LSED(LP)
            PEXP(L,NX) = 0.0  
            PHID(L,NX) = 0.0  
          enddo  
        enddo  

        do NX = 1,NSND  
          NS = NSED + NX  
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)  
            do NXX = 1,NSND  
              NSS = NSED + NXX  
              PEXP(L,NX) = PEXP(L,NX) + SNDB(L,K,NXX)*SEDDIA(NS) /(SEDDIA(NS)+SEDDIA(NSS))  
              PHID(L,NX) = PHID(L,NX) + SNDB(L,K,NXX)*SEDDIA(NSS)/(SEDDIA(NS)+SEDDIA(NSS))  
            enddo  
          enddo  
        enddo  

        do NX = 1,NSND  
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)  
            if( SNDBT(L,K) > 1E-12 )then  
              PEXP(L,NX) = PEXP(L,NX)/SNDBT(L,K)  
              PHID(L,NX) = PHID(L,NX)/SNDBT(L,K)  
            else  
              PEXP(L,NX) = 1.0  
              PHID(L,NX) = 1.0  
            endif  
          enddo  
        enddo  
      endif   ! *** END SET EXPOSURE AND HIDING FUNCTIONS FOR ENGULAND-HANSEN  

    enddo  ! *** END OF DOMAIN
    !$OMP END PARALLEL DO 

  endif     ! *** END OF ARMORING SECTION

  !**********************************************************************!  
  ! ***  SET CRITICAL SHIELD'S parameter FOR D50  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,DUM1,DUM3,DUM4)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)
    do LP = LF,LL
      L = LSED(LP)
      call SETSHLD(DUM1,CSHIELDS50(L),SEDDIA50(L,KBT(L)),SSG(NSED+1),DUM3,DUM4)  
    enddo  
  enddo
  !$OMP END PARALLEL DO 
  
  !**********************************************************************!  
  ! *** NONCOHESIVE SEDIMENT SCOUR, DEPOSITION AND VERTICAL PROCESSES 
  do NX = 1,NSND  
    NS = NSED + NX 
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,K,LP,L) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)

      !----------------------------------------------------------------------!  
      ! *** SET SETTLING VELOCITIES  
      if( ISNDVW == 0 .and. ISTOPT(7) == 1 .and. KC > 1 )then 
        ! *** Reset settling velocities back to user specified if using sediment anti-diffusion (ISTOPT(7) == 1)
        do K = 0,KS  
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSEDO(NS)  
          enddo  
        enddo  
      endif  

      if( ISNDVW >= 1 )then  
        do K = 0,KS  
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSEDO(NS)*CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)  
          enddo  
        enddo  
      endif  

      ! *** HANDLE LAYER 0 FOR SIGMA-ZED GRIDS
      if( IGRIDV > 0 )then
        ! *** ASSIGN LAYER 0 WHEN KSZ(L) > 1
        do LP = 1,LLWET(KC,ND)
          L = LKWET(LP,KC,ND)  
          if( KSZ(L) > 1 )then
            WSETA(L,0,NS) = WSETA(L,KSZ(L)-1,NS) 
          endif
        enddo
      endif
      
      ! *** ZERO SETTLING VELOCITIES AT OPEN BOUNDARIES
      do K = 0,KS
        do LP = 1,LLWET(K+1,ND)
          L = LKWET(LP,K+1,ND)  
          WSETA(L,K,NS) = SPB(L)*WSETA(L,K,NS)
        enddo
      enddo
    enddo  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO 

    !----------------------------------------------------------------------!  
    ! *** COMPUTE THE BEDLOAD COMPONENT FOR THE CURRENT SEDIMENT CLASS
    if( ICALC_BL > 0 )then
      DSEDGMM  = 1./(1.E6*SSG(NS))          ! *** SPECIFIC VOLUME (M**3/G)
      DIASED   = SEDDIA(NS)                 ! *** NOMINAL SEDIMENT PARTICLE DIAMTER (M)
      GPDIASED = G*(SSG(NS)-1.)*DIASED      ! *** "EXCESS" DENSITY (M2/S2)

      call BEDLOAD(NX,NS)
    endif
  enddo  ! *** END OF NON-COHESIVE CLASS LOOP
      
  ! *****************************************************************************************
  ! *** Communicate Bedload, if needed
  if( ICALC_BL > 0 )then
    ! ****************************************************************************
    TTDS = DSTIME(0)
    call MPI_barrier(DSIcomm, ierr)
    TWAIT = DSTIME(0) - TTDS
    TTSED = TTSED - TWAIT
  
    TTDS = DSTIME(0)
    call Communicate_BEDLOAD(1,NSND)
    DSITIMING(8) = DSITIMING(8) + (DSTIME(0) - TTDS)
    ! ****************************************************************************

    ! *****************************************************************************************
    ! *** FINISH BEDLOAD CALCULATES AFTER GHOST CELL COMMUNICATION
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,NX,LP,L,LE,LN)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      do NX = 1,NSND  
        ! *** CALCULATE MASS PER UNIT AREA CHANGE IN BED CONCENTRATION DUE TO TO NET BED LOAD (G/M2/S)                                                                                                  
        do LP = LF,LL
          L = LSED(LP)  
          LE = LEC(L)
          LN = LNC(L)
          ! ***                              U COMPONENT                    V COMPONENT
          SNDFBL(L,NX) = DXYIP(L)*( QSBDLDX(LE,NX)-QSBDLDX(L,NX) + QSBDLDY(LN,NX)-QSBDLDY(L,NX) )
        enddo
      enddo  ! *** END OF NON-COHESIVE CLASS LOOP
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
  
    ! *****************************************************************************************
    ! *** COMPUTE CELL CENTERED BEDLOAD FLUXES FOR OUTFLOW OR RECIRCULATION BOUNDARY  
    if( NSBDLDBC > 0 )then
      do NX = 1,NSND
        do NSB = 1,NSBDLDBC
          LUTMP = LSBLBCU(NSB)
          LDTMP = LSBLBCD(NSB)
        
          QSBDLDOT(LUTMP,NX) = SNDFBL(LUTMP,NX)*DXYP(LUTMP)
          if( LDTMP > 0 )then
            QSBDLDIN(LDTMP,NX) = QSBDLDIN(LDTMP,NX) + QSBDLDOT(LUTMP,NX)     
            SNDFBL(LDTMP,NX) = SNDFBL(LDTMP,NX) + QSBDLDOT(LUTMP,NX)*DXYIP(LDTMP)
          endif
          
          SNDFBL(LUTMP,NX) = 0.0
        enddo
      enddo
    endif
  endif

  ! *****************************************************************************************
  ! *** NONCOHESIVE SEDIMENT SCOUR, DEPOSITION AND VERTICAL PROCESSES 
  do NX = 1,NSND  
    NS = NSED + NX
    
    DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
    DIASED   = SEDDIA(NS)  
    GPDIASED = G*(SSG(NS)-1.)*DIASED  

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND, LF, LL, LP, L, K, KTOP, LE, LW, LN, LS, L1) &
    !$OMP             PRIVATE(WVEL, CLEFT, CRIGHT, WSFAC, WESE, ZEQMAX, SIGP, CSHIELDS, SHIELDS, TMPVAL, SNDBTMP, WSETMP)  &
    !$OMP             PRIVATE(TOP, BOT, FLUXFAC, CRNUM, GRADSED, SEDAVG, NSB, LUTMP, LDTMP) &
    !$OMP  SHARED(NDM, LDMWET, LAWET, LLWET, LKWET, LWET, KC, SNDF, SND, SND1, SNDS, WSETA, USTAR, USTARSND, RSNDM, ISNDM1, ISNDM2, NX, NS, KS) &
    !$OMP  SHARED(DIASED, DSEDGMM, GPDIASED, DELTI, SNDB, SNDB1, KBT, SNDFBL, TAUR, TAUBSND, SEDDIA50, HP, ISNDEQ, SSG, FACSUSL, CSHIELDS50) &
    !$OMP  SHARED(DZC, ZEQ, ZEQD, ZEQDI, SIGPHI, VDRBED, ISNDAL, ISLTAUC, TCSHIELDS, ISEDEFF, PHID, PEXP, VFRBED, VDRDEPO, HPK, HPKI, DZIC)  &
    !$OMP  SHARED(ROUSE, COEHEFF2, COEHEFF, FRACCOH, SNDEQB, SNDEQ, SNDEQSAV, DTSED, LEC, LNC, S2TL, S3TL, N, IBMECH, SEDVRDT, VKC) &
    !$OMP  SHARED(QSBDLDX, QSBDLDY, QSBDLDOT, DXYIP, QSBDLDIN, ISTOPT, LMASKDRY, LSC, LWC, ISBEDMAP, IROUSE, IFLAG, ISPROPWASH) &
    !$OMP  SHARED(TVAR1S, TVAR2S, TVAR3S, TVAR1N, TVAR2N, TVAR1W, TVAR1E, KSZ, ICALC_BL, NSBDLDBC, LSBLBCU, LSBLBCD, LBED, PROP_ERO) 
    do ND = 1, NDM  
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)
      
      ! *** SET THE ROUSE parameter
      if( IROUSE(NX) == 0 )then  
        do LP = LF,LL
          L = LWET(LP)
          if( LBED(L) ) CYCLE
          if( USTAR(L) > 1E-12 )then  
            ROUSE(L) = WSETA(L,0,NS)/(VKC*USTAR(L))  
          else  
            ROUSE(L) = 250000.*WSETA(L,0,NS)  
          endif  
        enddo  
      else  
        do LP = LF,LL
          L = LWET(LP)
          if( LBED(L) ) CYCLE
          if( USTARSND(L) > 1E-12 )then  
            ROUSE(L) = WSETA(L,0,NS)/(VKC*USTARSND(L))  
          else  
            ROUSE(L) = 250000.*WSETA(L,0,NS)  
          endif  
        enddo  
      endif  

      ! ----------------------------------------------------------------------
      ! *** CALCULATE WATER COLUMN SETTLING
      
      ! *** SET FLUX FOR THE BOTTOM OF THE TOP LAYER
      if( KC > 1 )then
        K = KC  
        do LP = LF,LL
          L = LWET(LP)
          SNDF(L,K,NX) = 0.                              !  (G/M2/S)
          WVEL = DTSED*HPKI(L,K)  
          CLEFT = 1.+WSETA(L,K-1,NS)*WVEL  
          CRIGHT = max(SND(L,K,NX),0.)  
          SND(L,K,NX) = CRIGHT/CLEFT  
          SNDF(L,K-1,NX) = -WSETA(L,K-1,NS)*SND(L,K,NX)  
        enddo  

        ! *** SET FLUX FOR MIDDLE LAYERS
        do K = KS,2,-1  
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            WVEL = DTSED*HPKI(L,K)  
            CLEFT = 1.+WSETA(L,K-1,NS)*WVEL  
            CRIGHT = max(SND(L,K,NX),0.)-SNDF(L,K,NX)*WVEL  
            SND(L,K,NX) = CRIGHT/CLEFT  
            SNDF(L,K-1,NX) = -WSETA(L,K-1,NS)*SND(L,K,NX)  
          enddo  
        enddo  

      endif
    
      ! *** Handle hard bottom for bottom water layer
      if( ISBEDMAP > 0 )then
        do LP = LF,LL
          L = LWET(LP)
          if( LBED(L) )then
            ! *** ADD SETTLING FROM THE LAYER ABOVE
            K = KSZ(L)
            WVEL = DTSED*HPKI(L,K)  
            SND(L,K,NX) = max(SND(L,K,NX),0.) - SNDF(L,K,NX)*WVEL  
          endif
        enddo
      endif
    
      ! *** FSEDMODE SETS BEDLOAD (IMODE = 1) AND SUSPENDED LOAD (IMODE = 2) TRANSPORT FRACTIONS
      do LP = LF,LL
        L = LWET(LP)
        FACSUSL(L) = FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)  
      enddo  

      ! ----------------------------------------------------------------------
      ! *** UPDATE SEDIMENT BED MASS & BOTTOM LAYER CONCENTRATION
      IFLAG = 0
      WSFAC = 1.0  ! *** FOR KC>1
      do LP = LF,LL
        L = LWET(LP)
        if( LBED(L) ) CYCLE

        WESE = 0.0
        if( ISNDEQ(NS) == 0 )then
          ! *** USER SPECIFIED EQUILIBRIUM CONCENTRATION
          ZEQ(L)    = 0.05         ! *** GARCIA & PARKER
          ZEQD(L)   = ZEQ(L)*DZIC(L,KSZ(L))  
          ZEQDI(L)  = 1./ZEQD(L)  
          SNDEQB(L) = TAUR(NS)
        else
          ! *** SET EQUILIBRUIM CONCENTRATION REFERENCE HEIGHT (DIMENSIONLESS)  
          ZEQ(L) = CSNDZEQ(ISNDEQ(NS),DIASED,GPDIASED,TAUR(NS),TAUBSND(L),SEDDIA50(L,KBT(L)),HP(L),SSG(NS),WSETA(L,0,NS))
          
          ZEQMAX   = 0.5*DZC(L,KSZ(L))  
          ZEQ(L)   = min(ZEQ(L),ZEQMAX)  
          ZEQD(L)  = ZEQ(L)*DZIC(L,KSZ(L))  
          ZEQDI(L) = 1./ZEQD(L)  
          SIGP = SIGPHI(L,KBT(L))  

          ! *** SET EQUILIBRUIM CONCENTRATION  
          SNDEQB(L) = CSNDEQC(ISNDEQ(NS),DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUBSND(L),SEDDIA50(L,KBT(L)),SIGP,ZEQ(L),VDRBED(L,KBT(L)),ISNDAL)      !  

          ! *** APPLIED LIMITOR TO GARCIA AND PARKER  
          if( ISNDEQ(NS) == 1 )then  
            if( ISLTAUC(NS) == 1 )then  
              CSHIELDS = TCSHIELDS(NS)  
              if( ISEDEFF == 2 )then  
                TMPVAL = 1. + (COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
                CSHIELDS = TMPVAL*CSHIELDS  
              endif  
              SHIELDS = TAUBSND(L)/GPDIASED  
              if( SHIELDS < CSHIELDS) SNDEQB(L) = 0.  
            endif  
            if( ISLTAUC(NS) == 2 )then  
              CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED  
              if( ISEDEFF == 2 )then  
                TMPVAL = 1. + (COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
                CSHIELDS = TMPVAL*CSHIELDS  
              endif  
              SHIELDS = TAUBSND(L)/GPDIASED  
              if( SHIELDS < CSHIELDS) SNDEQB(L) = 0.  
            endif  
            if( ISLTAUC(NS) == 3 )then  
              CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)  
              if( ISEDEFF == 2 )then  
                TMPVAL = 1. + (COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
                CSHIELDS = TMPVAL*CSHIELDS  
              endif  
              SHIELDS = TAUBSND(L)/GPDIASED  
              if( SHIELDS < CSHIELDS) SNDEQB(L) = 0.  
            endif  
          endif  

          if( ISEDEFF == 1 ) SNDEQB(L) = SNDEQB(L)*EXP(-COEHEFF*FRACCOH(L,KBT(L)))  
        endif    ! *** ISNDEQ(NS) == 0
        
        if( ROUSE(L) < 0.999 .or. ROUSE(L) > 1.001 )then  
          TOP = (ZEQD(L)**(ROUSE(L)-1.))-1.  
          BOT = (1.-ROUSE(L))*(ZEQDI(L)-1.)  
          SNDEQ(L) = SNDEQB(L)*TOP/BOT  
          SNDEQ(L) = FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)  
          SNDEQSAV(L,NX) = SNDEQ(L)  
        else  
          TOP = LOG(ZEQDI(L))  
          BOT = (ZEQDI(L)-1.)  
          SNDEQ(L) = SNDEQB(L)*TOP/BOT  
          SNDEQ(L) = FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)  
          SNDEQSAV(L,NX) = SNDEQ(L)  
        endif  

        ! *** SET RESUSPENSION FLUX
        if( KSZ(L) == KC )then    ! *** Alberta
          WSFAC = 2.*(1.+ROUSE(L))/(2.+ROUSE(L)*(1.-ZEQ(L)))  
        endif
        WESE = WSFAC*WSETA(L,0,NS)*SNDEQ(L)  

        ! *** SET DEPOSITION VELOCITY  
        WSETMP = WSFAC*WSETA(L,0,NS)  
        WVEL   = DTSED*HPKI(L,KSZ(L))
        
        ! *** Handle Propwash
        if( PROP_ERO(L,0) > 0.0 )then
          PROP_ERO(L,NS) = PROP_ERO(L,NS)*DXYIP(L)*DELTI                        ! *** Convert mass from g to g/m**2/s
          WSETMP = 0.0                                                          ! *** Disable settling for active propwash cell
          !IF( ISPROPWASH == 2 )then
          !  WESE = PROP_ERO(L,NS)                                              ! *** Only allow propwash erosion to avoid double counting
          !ELSE
            WESE = WESE + PROP_ERO(L,NS)                                        ! *** Allow ambient current erosion for all propwash options (2021-11-10)
          !ENDIF
        endif
        
        ! *** SET BOTTOM LAYER SND CONCENTRATION AND FLUX TO SUSPENDED LOAD
        CLEFT  = 1. + WSETMP*WVEL  
        CRIGHT = max(SND(L,KSZ(L),NX),0.) + ( WESE-SNDF(L,KSZ(L),NX) )*WVEL  
        SND(L,KSZ(L),NX) = CRIGHT/CLEFT  
        SNDF(L,0,NX) = -WSETMP*SND(L,KSZ(L),NX) + WESE  

        ! *** ADD BED LOAD FLUX TO SUSPENDED LOAD FLUX
        SNDBTMP = SNDB(L,KBT(L),NX) - DTSED*SNDF(L,0,NX) - DTSED*SNDFBL(L,NX)  

        ! *** HANDLE CASE IF INSUFFICIENT BED MATERIAL
        if( SNDBTMP < 0.0 )then  
          ! *** ADJUST BEDLOAD FLUX, TRYING TO KEEP SUSPENDED LOAD FLUX CONSTANT
          SNDBTMP = SNDB(L,KBT(L),NX) - DTSED*SNDF(L,0,NX) 
          LE = LEC(L)
          LN = LNC(L)
          if( SNDBTMP < 0.0 )then
            if( ICALC_BL > 0 )then
              IFLAG = 2
              ! *** ZERO BEDLOAD OUTFLUX
              QSBDLDX(L,NX) = max(QSBDLDX(L,NX),0.0)
              QSBDLDY(L,NX) = max(QSBDLDY(L,NX),0.0)
              QSBDLDX(LE,NX) = min(QSBDLDX(LE,NX),0.0)
              QSBDLDY(LN,NX) = min(QSBDLDY(LN,NX),0.0)
              if( NSBDLDBC > 0 )then
                do NSB = 1,NSBDLDBC
                  LUTMP = LSBLBCU(NSB)
                  LDTMP = LSBLBCD(NSB)
                  if( L == LUTMP )then
                    if( LDTMP > 0 ) QSBDLDIN(LDTMP,NX) = QSBDLDIN(LDTMP,NX) - QSBDLDOT(LUTMP,NX)
                    QSBDLDOT(LUTMP,NX) = 0.
                    exit
                  endif
                enddo
              endif
            endif
            
            ! *** REDUCE SUSPENDED LOAD FLUX
            SNDF(L,0,NX)     = SNDB(L,KBT(L),NX)/DTSED
            SND(L,KSZ(L),NX) = SNDS(L,KSZ(L),NX) + (SNDF(L,0,NX) - SNDF(L,KSZ(L),NX))*WVEL  
          else
            ! *** REDUCE BEDLOAD FLUX
            IFLAG = 1
            SNDBTMP = SNDBTMP/DTSED
            FLUXFAC = SNDBTMP/SNDFBL(L,NX)

            if( QSBDLDX(L,NX) < 0. ) QSBDLDX(L,NX) = FLUXFAC*QSBDLDX(L,NX)  
            if( QSBDLDY(L,NX) < 0. ) QSBDLDY(L,NX) = FLUXFAC*QSBDLDY(L,NX)  
            if( QSBDLDX(LE,NX) > 0. ) QSBDLDX(LE,NX) = FLUXFAC*QSBDLDX(LE,NX)  
            if( QSBDLDY(LN,NX) > 0. ) QSBDLDY(LN,NX) = FLUXFAC*QSBDLDY(LN,NX)  
            if( NSBDLDBC > 0 )then
              do NSB = 1,NSBDLDBC
                LUTMP = LSBLBCU(NSB)
                LDTMP = LSBLBCD(NSB)
                if( L == LUTMP )then
                  if( LDTMP > 0 ) QSBDLDIN(LDTMP,NX) = QSBDLDIN(LDTMP,NX) - QSBDLDOT(LUTMP,NX)
                  QSBDLDOT(L,NX) = FLUXFAC*QSBDLDOT(L,NX) 
                  exit
                endif
              enddo
            endif
          endif
        endif  
      enddo 
        
      !----------------------------------------------------------------------!  
      ! ***  ANTI-DIFFUSION OF NON-COHESIVE SEDIMENT  KC > 1
      if( ISTOPT(7) == 1 .and. KC > 1 )then
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            CRNUM = 1.+DTSED*WSETA(L,K,NS)*HPKI(L,K+1)
            GRADSED = (SND(L,K+1,NX)-SND(L,K,NX))/(DZC(L,K+1)+DZC(L,K))  
            SEDAVG = 0.5*(SND(L,K+1,NX)-SND(L,K,NX)+1.E-16)  
            WSETA(L,K,NS) = -CRNUM*DZC(L,K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG  
          enddo  
        enddo  

        ! *** TVAR1S = LOWER DIAGONAL  
        do LP = LF,LL
          L = LWET(LP)
          TVAR1S(L,KSZ(L)) = 0  
        enddo  
        do K = 2,KC  
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            TVAR1S(L,K) = min(WSETA(L,K-1,NS),0.)  
          enddo  
        enddo  

        ! *** TVAR1N = UPPER DIAGONAL  
        do LP = LF,LL
          L = LWET(LP)
          TVAR1N(L,KC) = 0  
        enddo  
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TVAR1N(L,K) = -MAX(WSETA(L,K,NS),0.)  
          enddo  
        enddo  

        ! *** TVAR1W = MAIN DIAGONAL  
        do LP = LF,LL
          L = LWET(LP)
          TVAR1W(L,KSZ(L)) = DELTI*HPK(L,KSZ(L))-MIN(WSETA(L,KSZ(L),NS),0.)
          TVAR1W(L,KC) = DELTI*DZC(L,KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)  
        enddo  
        do K = 2,KS  
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            TVAR1W(L,K) = DELTI*DZC(L,KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)-MIN(WSETA(L,K,NS),0.)  
          enddo  
        enddo  

        ! *** TVAR1E = RIGHT HAND SIDE  
        do K = 1,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TVAR1E(L,K) = DELTI*DZC(L,KC)*HP(L)*SND(L,K,NX)  
          enddo  
        enddo  

        ! *** TVAR3S = BET,TVAR2N = U,TVAR2S = GAM ARE WORKING ARRAYS  
        do LP = LF,LL
          L = LWET(LP)
          TVAR3S(L) = TVAR1W(L,KSZ(L))  
        enddo  
        do LP = LF,LL
          L = LWET(LP)
          TVAR2N(L,KSZ(L)) = TVAR1E(L,KSZ(L))/TVAR3S(L)  
        enddo  
        do K = 2,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TVAR2S(L,K) = TVAR1N(L,K-1)/TVAR3S(L)  
            TVAR3S(L) = TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)  
            TVAR2N(L,K) = (TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/TVAR3S(L)  
          enddo  
        enddo  
        do K = KS,1,-1  
          do L = LF,LL  
            TVAR2N(L,K) = TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)  
          enddo  
        enddo  
        do K = 1,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            SND(L,K,NX) = TVAR2N(L,K)  
          enddo  
        enddo  

      endif   ! *** END OF ANTI DIFFUSION
      
      !----------------------------------------------------------------------!  
      ! *** FINAL FLUX  
      if( KC > 1 )then
        do LP = 1,LLWET(KS,ND)
          L = LKWET(LP,KS,ND)  
          SNDF(L,KS,NX) = DELTI*DZC(L,KC)*HP(L)*(SND(L,KC,NX)-SNDS(L,KC,NX))  
        enddo  

        do K = KS-1,1,-1  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            SNDF(L,K,NX) = DELTI*DZC(L,K+1)*HP(L)*(SND(L,K+1,NX)-SNDS(L,K+1,NX)) + SNDF(L,K+1,NX)  
          enddo  
        enddo  
      endif
            
    enddo  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO 
    
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,KTOP,LE,LW,LS,LN,SNDBTMP,WVEL)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)
      
      if( ICALC_BL > 0 .and. IFLAG > 0 )then
        ! *** RECALCULATE THE BEDLOAD FLUXES CONSISTENT WITH THE NEW FACE FLUXES
        do LP = LF,LL
          L = LSED(LP)
          LE = LEC(L)
          LN = LNC(L)
        
          ! *** UPDATE BEDLOAD FLUX FOR CURRENT CELL
          ! ***                              U COMPONENT                    V COMPONENT
          SNDFBL(L,NX) = DXYIP(L)*( QSBDLDX(LE,NX)-QSBDLDX(L,NX) + QSBDLDY(LN,NX)-QSBDLDY(L,NX) + QSBDLDOT(L,NX)-QSBDLDIN(L,NX) )  
          
          if( NSBDLDBC > 0 )then
            do NSB = 1,NSBDLDBC
              LUTMP = LSBLBCU(NSB)
              if( L == LUTMP )then
                LDTMP = LSBLBCD(NSB)
                if( LDTMP > 0 )then
                  SNDFBL(LDTMP,NX) = SNDFBL(LDTMP,NX) + SNDFBL(LUTMP,NX)
                endif
                SNDFBL(LUTMP,NX) = 0.0
                exit
              endif
            enddo
          endif
        enddo
      endif
      
      ! *** UPDATE BED MASS AND VOLUME FLUXES
      do LP = LF,LL
        L = LSED(LP)
        KTOP = KBT(L)
        
        ! *** ADD BED LOAD FLUX TO SUSPENDED LOAD FLUX
        SNDBTMP = SNDB(L,KTOP,NX) - DTSED*SNDF(L,0,NX) - DTSED*SNDFBL(L,NX)  
        
        ! *** CHECK ONE LAST TIME FOR MASS SUFFICIENCY
        if( SNDBTMP < 0.0 )then
          if( SNDFBL(L,NX) > SNDF(L,0,NX) )then
            ! *** REDUCE BEDLOAD FLUX AND ZERO SUSPENDED LOAD FLUX
            SNDFBL(L,NX) = SNDB(L,KTOP,NX)/DTSED
            SNDF(L,0,NX) = 0.0
          else
            ! *** REDUCE SUSPENDED LOAD FLUX AND ZERO BEDLOAD FLUX
            SNDFBL(L,NX) = 0.0
            SNDF(L,0,NX) = SNDB(L,KTOP,NX)/DTSED
          endif
          WVEL   = DTSED*HPKI(L,KSZ(L))
          SND(L,KSZ(L),NX) = SNDS(L,KSZ(L),NX) + (SNDF(L,0,NX) - SNDF(L,KSZ(L),NX))*WVEL    ! *** REDUCE SUSPENDED LOAD BY THE REDUCTION IN SNDF
          SNDBTMP = 0.0                                                                     ! *** ZERO BED SEDIMENT FOR CLASS NX
        endif
        
        SNDB1(L,KTOP,NX) = S3TL*SNDB(L,KTOP,NX) + S2TL*SNDB1(L,KTOP,NX)  
        SNDB(L,KTOP,NX)  = SNDBTMP  
        SNDF(L,0,NX) = SNDF(L,0,NX) + SNDFBL(L,NX)       ! *** BED/WATER INTEFACE SND FLUX COMBINES BEDLOAD AND SUSPENSION

        ! *** EROSION/DEPOSITION RATE DUE TO SND CLASS NX SOLIDS (QSBDTOP) AND VOID FRACTION (QWBDTOP)  (M/S)
        QSBDTOP(L) = QSBDTOP(L) + DSEDGMM*SNDF(L,0,NX)   
        if( IBMECH == 0 .or. SEDVRDT < 0.0 )then
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*VDRBED(L,KTOP)*SNDF(L,0,NX)  ! *** IF EITHER CONSTANT OR INSTANTLY CONSOLIDATING, USE BED VR  
        else
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*( VDRBED(L,KTOP)*MAX(SNDF(L,0,NX),0.) + VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
        endif
      enddo

      ! ******************************************************************************
      ! *** MOVE ANY EXITING BEDLOAD FROM AN ACTIVE CELL TO HARD BOTTOM SUSPENDED LOAD                                                                                               
      if( ND == NDM .and. ICALC_BL > 0 .and. ISBEDMAP > 0 )then
        ! *** EAST CELL IS A HARD BOTTOM CELL
        do LP = 1,BEDEDGEE.NEDGE
          L = BEDEDGEE.LEDGE(LP)
          LE = LEC(L)
          if( QSBDLDX(LE,NX) > 0. )then 
            SNDBTMP = QSBDLDX(LE,NX)*DXYIP(LE)*DTSED*HPKI(LE,KSZ(LE))   ! EQUIVALENT MG/L
            SND(LE,KSZ(LE),NX) = SND(LE,KSZ(LE),NX) + SNDBTMP
          endif
        enddo

        ! *** WEST CELL IS A HARD BOTTOM CELL
        do LP = 1,BEDEDGEW.NEDGE
          L = BEDEDGEW.LEDGE(LP)
          LW = LWC(L)
          if( QSBDLDX(L,NX) < 0. )then 
            SNDBTMP = -QSBDLDX(L,NX)*DXYIP(LW)*DTSED*HPKI(LW,KSZ(LW))   ! EQUIVALENT MG/L
            SND(LW,KSZ(LW),NX) = SND(LW,KSZ(LW),NX) + SNDBTMP
          endif
        enddo    
    
        ! *** SOUTH CELL IS A HARD BOTTOM CELL
        do LP = 1,BEDEDGES.NEDGE
          L = BEDEDGES.LEDGE(LP)
          LS = LSC(L)
          if( QSBDLDY(L,NX) < 0. )then 
            SNDBTMP = -QSBDLDY(L,NX)*DXYIP(LS)*DTSED*HPKI(LS,KSZ(LS))   ! EQUIVALENT MG/L
            SND(LS,KSZ(LS),NX) = SND(LS,KSZ(LS),NX) + SNDBTMP
          endif
        enddo    

        ! *** NORTH CELL IS A HARD BOTTOM CELL
        do LP = 1,BEDEDGEN.NEDGE
          L = BEDEDGEN.LEDGE(LP)
          LN = LNC(L)
          if( QSBDLDY(LN,NX) > 0. )then 
            SNDBTMP = QSBDLDY(LN,NX)*DXYIP(LN)*DTSED*HPKI(LN,KSZ(LN))   ! EQUIVALENT MG/L
            SND(LN,KSZ(LN),NX) = SND(LN,KSZ(LN),NX) + SNDBTMP
          endif
        enddo    
      endif
  
    enddo  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO 

  enddo  ! *** END OF NON-COHESIVE CLASS LOOP

  ! ********************************************************************************
  ! *** DO A CHECK FOR SND AND SNDB < 0.0 AND LOG IT
  if( ISDTXBUG > 0 )then
    IFLAG = 0  
    do NS = 1,NSND  
      do L = 2,LA  
        K = KSZ(L)
        if( SND(L,K,NS) < 0. )then  
          if( IFLAG == 0 )then  
            open(1,FILE = OUTDIR//'NEGSEDSND.OUT',POSITION = 'APPEND')  
            IFLAG = 1  
          endif  
          write(1,107) TIME,NS,IL(L),JL(L),K,SND(L,K,NS)  
          SND(L,K,NS) = 0.0    ! *** Continue with warning
        endif  
      enddo  
    enddo  

    do NS = 1,NSND  
      do L = 2,LA
        K = KBT(L)
        if( SNDB(L,K,NS) < 0. .or. HBED(L,K) < 0. .or. SEDDIA50(L,K) < 0.0 )then
          if( IFLAG == 0 )then  
            open(1,FILE = OUTDIR//'NEGSEDSND.OUT',POSITION = 'APPEND')  
            IFLAG = 1  
          endif  
          write(1,108) TIME,NS,IL(L),JL(L),K,SNDB(L,K,NS),SNDF(L,0,NS),HBED(L,K),SEDDIA50(L,K)*1000000.
          if( SNDB(L,K,NS) < 0. )  SNDB(L,K,NS) = 0.0    ! *** Continue with warning
          if( HBED(L,K) < 0. )     HBED(L,K) = 0.0       ! *** Continue with warning
          if( SEDDIA50(L,K) < 0. ) SEDDIA50(L,K) = 0.0   ! *** Continue with warning
        endif  
      enddo  
    enddo  

    if( IFLAG == 1 ) close(1) 
  endif
  
  107 FORMAT(' Warning: WC  SND < 0: TIME, NS, I, J, K, NEGSND = ',F12.4,4I5,4E13.4)       
  108 FORMAT(' Warning: BED SND < 0: TIME, NS, I, J, K, NEGSNDB, SNDF, HBED, D50 = ',F12.4,4I5,4E13.4)       

  !**********************************************************************!  
  return 

END  
