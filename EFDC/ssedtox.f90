! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SSEDTOX

  !----------------------------------------------------------------------C                                                
                          
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !---------------------------------------------------------------------------------
  ! 2019-01        Paul M. Craig       Added Hard Bottom Bypass
  ! 2016-12        Paul M. Craig       Implemented SEDZLJ for Toxics  
  ! 2016-12        Paul M. Craig       Changed bed layering order approach to swtich for each model
  !                                      SEDZLJ   -  1 (top) to KB (bottom)
  !                                      Original - KB (top) to  1 (bottom)
  ! 2014-08        D H CHUNG           SET EXPLICIT PRECISIONS OF INTEGER & REAL
  !  2013          Dang H  Chung       Adjusted Bed shear stress by wave
  !  2011          Paul M. Craig       Merged latest SNL, TT GVC and DSI codes
  !                                     Restructed code to F90 and OMP
  !  2008          Paul M. Craig       Merged SNL with DSI                                                                                         
  !----------------------------------------------------------------------C                                                
  !  2002          John Hamrick       Many updates to work with bedload, morphology and toxics
                          
  ! ***********************************************************************
                          
  ! *** SUBROUTINE SSEDTOX CALCULATES SETTLING AND WATER COLUMN-BED                                                       
  ! *** EXCHANGE OF SEDIMENT AND SORBED TOXIC CONTAMINANTS                                                                
                          
  ! ***********************************************************************CC                                               

  USE GLOBAL
  USE RESTART_MODULE
  Use Variables_Propwash
  
  IMPLICIT NONE

  INTEGER :: IFILE
  INTEGER :: NVAL, LP, L, K, NS, NX, KK, ND, LF, LL, NP, ipmc
  INTEGER :: KTOPTP, KTOPM1, NT, LE, LN, ITMP, LBANK, LCHAN
  INTEGER,SAVE :: NKINETICS

  real pmc0, pmc1, pmc2, pmc3, pmc4
  REAL :: HDZBR, DSEDGMMI, MINTHICK
  REAL :: UTMP, VTMP, CURANG, SHEAR, RAT1O
  REAL :: TAUBC, TAUB2, CSEDRESB, HDFUFXX, HDFUFYY, QSWNEG
  REAL :: TMPEXP, QCELLCTRA, TMPVAL, QSWPOS, TMP
  REAL :: ZOGRAINCOH, TMPTOP, TMPBOT, RVAL, SXD1, SXD, TTHICK
  REAL, EXTERNAL :: CSEDTAUS, CSEDRESS, CSEDTAUB

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DELBED
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: QCELLAD1SQ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: QCELLAD1ZZ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TAUBSEDS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SEDGEOSTD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SEDDIAGS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: ZOGRAIN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: ZOTOTAL
  
  REAL(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS,TTDSX        ! MODEL TIMING TEMPORARY VARIABLE
  
  INTEGER(IK4),SAVE,ALLOCATABLE,DIMENSION(:) :: LSEDPLUS
  INTEGER(IK4),SAVE                          :: NSEDPLUS
  
  IF( .NOT. ALLOCATED(DELBED) )THEN
    ALLOCATE(DELBED(LCM))
    ALLOCATE(QCELLAD1SQ(LCM))
    ALLOCATE(QCELLAD1ZZ(LCM))
    ALLOCATE(TAUBSEDS(LCM))
    ALLOCATE(SEDGEOSTD(LCM,KBM))
    ALLOCATE(SEDDIAGS(LCM,KBM))
    ALLOCATE(ZOGRAIN(LCM))
    ALLOCATE(ZOTOTAL(LCM))
    ALLOCATE(LSEDPLUS(LCM))
    
    DELBED=0.0
    QCELLAD1SQ=1.0  ! IF KC=1 ONLY NEED TO INITIALIZE ONCE
    QCELLAD1ZZ=1.0  ! IF KC=1 ONLY NEED TO INITIALIZE ONCE
    TAUBSEDS=0.0
    SEDGEOSTD=0.0
    SEDDIAGS=0.0
    ZOGRAIN=0.0
    ZOTOTAL=0.0
    NKINETICS=0
    DO NT=1,NTOX
      NKINETICS = NKINETICS+SUM(ITOXKIN(1:3,NT))
    ENDDO
    
    LSEDPLUS = 0
    NSEDPLUS = 0
    IF( ISBEDMAP > 0 )THEN
      ! *** EAST CELL IS A HARD BOTTOM CELL
      DO LP=1,BEDEDGEE.NEDGE
        L = BEDEDGEE.LEDGE(LP)
        IF( SUBO(LEC(L)) > 0. )THEN
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LEC(L)
        ENDIF
      ENDDO

      ! *** WEST CELL IS A HARD BOTTOM CELL
      DO LP=1,BEDEDGEW.NEDGE
        L = BEDEDGEW.LEDGE(LP)
        IF( SUBO(L) > 0. )THEN
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LWC(L)
        ENDIF
      ENDDO    
    
      ! *** SOUTH CELL IS A HARD BOTTOM CELL
      DO LP=1,BEDEDGES.NEDGE
        L = BEDEDGES.LEDGE(LP)
        IF( SVBO(L) > 0. )THEN
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LSC(L)
        ENDIF
      ENDDO    

      ! *** NORTH CELL IS A HARD BOTTOM CELL
      DO LP=1,BEDEDGEN.NEDGE
        L = BEDEDGEN.LEDGE(LP)
        IF( SVBO(LNC(L)) > 0. )THEN
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LNC(L)
        ENDIF
      ENDDO    
    ENDIF
  ENDIF

  TTDS = DSTIME(0)

  IFILE = -1

  ! *** S2TL AND S3TL ARE GLOBAL VARIABLES USED IN CALTOX AND CALSND
  IF( IS2TL == 1 )THEN
    S3TL=1.0
    S2TL=0.0
  ELSE
    IF( ISTL /= 3 )THEN
      S3TL=0.0
      S2TL=1.0
    ELSE
      S3TL=1.0
      S2TL=0.0
    ENDIF
  ENDIF

  ! *** ALL TIMING IS BASED ON THE SEDTIME
  DTSEDJ = REAL(DTSED,8)
  DELTI = 1._8/DTSEDJ

  SEDMDGM = SQRT(SEDMDMX*SEDMDMN)
  BEDEX = 1.
  NVAL = MOD(N,2)
  FOURDPI = 4./PI
  MINTHICK = MIN(1E-6,DTSED*1E-6)           ! *** Minimum bed thickness
  
  ! ***********************************************************************C                                                
  ! *** IF N=1 CALCULATE INITIAL SEDIMENT BED THICKNESS
                          
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,LN,LE,KK,NT) &
  !$OMP             PRIVATE(TAUBC,UTMP,VTMP,CURANG,TAUB2,HDFUFXX,HDFUFYY,TMPEXP,TMPVAL)  &
  !$OMP             PRIVATE(HDZBR,TMP,TMPTOP,TMPBOT,RVAL,ZOGRAINCOH,SHEAR,RAT1O,QCELLCTRA)
  DO ND=1,NDM  
    LF=(ND-1)*LDMSED+1  
    LL=MIN(LF+LDMSED-1,LASED)
                            
    ! *** SET SEDIMENT VOLUME FRACTIONS                                                                                     
    IF( .NOT. LSEDZLJ )THEN
      DO K=1,KB
        DO LP=LF,LL
          L = LSED(LP)
          BEDLINIT(L,K) = 0.
          BEDDINIT(L,K) = 0.
        ENDDO
      ENDDO
      DO NX=1,NSED+NSND
        DO K=1,KB
          DO LP=LF,LL
            L = LSED(LP)
            VFRBED(L,K,NX) = 0.
            VFRBED1(L,K,NX) = 0.
          ENDDO
        ENDDO
      ENDDO
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS)  = SDEN(NS)*SEDB(L,K,NS)
              VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS)  = SDEN(NS)*SNDB(L,K,NX)
              VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
              BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
              BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              IF( BEDLINIT(L,K) > 0.0 )THEN
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              ELSE
                VFRBED(L,K,NS) = 0.0
              ENDIF
              IF( BEDDINIT(L,K) > 0.0 )THEN
                VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
              ELSE
                VFRBED1(L,K,NS) = 0.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              IF( BEDLINIT(L,K) > 0.0 )THEN
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              ELSE
                VFRBED(L,K,NS) = 0.0
              ENDIF
              IF( BEDDINIT(L,K) > 0.0 )THEN
                VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
              ELSE
                VFRBED1(L,K,NS) = 0.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! *** INITIALIZE FRACTIONS FOR THE CURRENT TIME STEP
      DO K=1,KB
        DO LP=LF,LL
          L = LSED(LP)
          FRACCOH(L,K) = 0.0
          FRACNON(L,K) = 0.0
        ENDDO
      ENDDO
      DO NS=1,NSED
        DO K=1,KB
          DO LP=LF,LL
            L = LSED(LP)
            IF( K <= KBT(L) )THEN
              FRACCOH(L,K) = FRACCOH(L,K)+VFRBED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DO NX=1,NSND
        NS=NX+NSED
        DO K=1,KB
          DO LP=LF,LL
            L = LSED(LP)
            IF( K <= KBT(L) )THEN
              FRACNON(L,K) = FRACNON(L,K)+VFRBED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ! *** INITIALIZE SEDIMENT BED/WATER COLUMN INTERFACE FLUXES FOR THE CURRENT TIME STEP
      DO LP=LF,LL
        L = LSED(LP)
        QWBDTOP(L) = 0.   ! *** VOLUME OF WATER EXCHANGE DUE TO ENTRAPMENT/SCOUR     (M/S)
        QSBDTOP(L) = 0.   ! *** VOLUME OF SEDIMENT EXCHANGE                          (M/S)
      ENDDO

      ! *** SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES                                                         
      IF( ISTRAN(6) >= 1 )THEN 
        IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN
          ! *** COMPUTE CRITICAL SHEAR STRESS FOR RESUSPENSION
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)  ! PMC - REMOVED LOOP AND ONLY COMPUTE THE TOP LAYER
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1), VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1))
          ENDDO
        ENDIF
        IF( IWRSPB(1) >= 1 )THEN
          ! *** COMPUTE CRITICAL SHEAR STRESS FOR BULK/MASS EROSION
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)  ! PMC - REMOVED LOOP AND ONLY COMPUTE THE TOP LAYER
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))
            WRSPB(L,K) = CSEDRESB
          ENDDO
        ENDIF
      ENDIF
    ENDIF    ! *** END OF LSEDZLJ BYPASS
                              
    ! *** CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC                                                        
    ! *** CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG                                                        
    ! *** TO TOXB UNITS OF OF MG/M**2                                                                                       
                            
    ! *** SAVE OLD VALUES                                                                                                   
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED2
        DO K=1,KC
          DO LP=LF,LL
            L = LSED(LP)
            SEDS(L,K,NS) = SED(L,K,NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    
    IF( .NOT. LSEDZLJ )THEN
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          DO K=1,KC
            DO LP=LF,LL
              L = LSED(LP)
              SNDS(L,K,NX) = SND(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! ***********************************************************************C                                                
      ! *** SET MEAN D50 AND D90                                                                                              
      IF( ISTRAN(7) >= 1 )THEN
        ! *** SEDIMENT MASS
        DO K=1,KB
          DO LP=LF,LL
            L = LSED(LP)
            SNDBT(L,K) = 0.
          ENDDO
        ENDDO
        
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NX)
            ENDDO
          ENDDO
        ENDDO

        IF( ISBSDIAM > 1 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            SEDDIA90(L,KBT(L)) = 0.
            SEDGEOSTD(L,KBT(L)) = 0.
          ENDDO
        ENDIF
        
        ! *** D50
        DO LP=LF,LL
          L = LSED(LP)
          SEDDIA50(L,KBT(L)) = 0.
        ENDDO
        DO NX=1,NSND
          NS=NSED+NX
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)
            SEDDIA50(L,K) = SEDDIA50(L,K) + SNDB(L,K,NX)*LOG(SEDDIA(NS))
          ENDDO
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          K = KBT(L)
          IF( SNDBT(L,K) > 0. )THEN
            SEDDIA50(L,K) = SEDDIA50(L,K)/SNDBT(L,K)
          ENDIF
        ENDDO
                            
        IF( ISBSDIAM > 1 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            DO LP=LF,LL
              L = LSED(LP)
              K = KBT(L)
              SEDGEOSTD(L,K) = SEDGEOSTD(L,K) + SNDB(L,K,NX)*(( LOG(SEDDIA(NS))-SEDDIA50(L,K) )**2)
            ENDDO
          ENDDO
                            
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)
            IF( SNDBT(L,K) > 0. )THEN
              SEDGEOSTD(L,K) = SEDGEOSTD(L,K)/SNDBT(L,K)
            ENDIF
          ENDDO
        ENDIF
        
        DO LP=LF,LL
          L = LSED(LP)
          K = KBT(L)
          ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
          IF( SNDBT(L,K) > 5.0 )THEN  
            SEDDIA50(L,K)  = EXP(SEDDIA50(L,K))
          ELSEIF( SNDBT(L,K) > 1E-12 )THEN
            ! *** MASS WEIGHTED D50
            SEDDIA50(L,K) = 0.
            DO NX=1,NSND
              NS=NSED+NX
              SEDDIA50(L,K) = SEDDIA50(L,K) + SNDB(L,K,NX)*SEDDIA(NS)
            ENDDO
            SEDDIA50(L,K)  = SEDDIA50(L,K)/SNDBT(L,K)
          ELSE
            SEDDIA50(L,K)  = SEDDIA(NSED+1)
          ENDIF
        ENDDO
                            
        IF( ISBSDIAM > 1 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)
            ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
            IF( SNDBT(L,K) > 5.0 )THEN  
              SEDGEOSTD(L,K) = EXP(SEDGEOSTD(L,K))
            ELSE
              SEDGEOSTD(L,K) = 1.0
            ENDIF
          ENDDO
                            
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)
            SEDDIA90(L,K) = (SEDGEOSTD(L,K)**1.28)*SEDDIA50(L,K)
          ENDDO
        ENDIF
        
        ! *** SET SEDIMENT GRAINSIZE TO USE FOR SHEAR STRESS CALCS  
        IF( ISBSDIAM == 0 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = SEDDIA50(L,KBT(L))
          ENDDO
        ELSEIF( ISBSDIAM == 1 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = 2.*SEDDIA50(L,KBT(L))
          ENDDO
        ELSEIF( ISBSDIAM == 2 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = SEDDIA90(L,KBT(L))
          ENDDO
        ELSEIF( ISBSDIAM == 3 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = 2.*SEDDIA90(L,KBT(L))
          ENDDO
        ENDIF                                                                                                                                        
      ENDIF    ! *** END OF NON-COHESIVE GRAIN SIZE CALCS
    ENDIF    ! *** END OF LSEDZLJ BYPASS
                            
    ! ***********************************************************************C                                                
    ! *** SET CELL CENTER BED STRESS FOR SEDIMENT RESUSPENSION AND DEPOSITION                                               
    IF( .NOT. LSEDZLJ )THEN
      IF( ISWAVE > 0 )THEN
        ! *** WITH WAVE INTERACTIONS
        DO LP=LF,LL
          L = LSED(LP)
          IF( LWVMASK(L) )THEN
            TAUBC = QQ(L,0)/CTURB2  ! *** CURRENT PLUS BOTTOM SHEAR FROM QQWV1
            UTMP = 0.5*STCUV(L)*(U(LEC(L),KSZ(L)) + U(L,KSZ(L)))+1.E-12
            VTMP = 0.5*STCUV(L)*(V(LNC(L),KSZ(L)) + V(L,KSZ(L)))
            CURANG = ATAN2(VTMP,UTMP)
            TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
            TAUB2 = MAX(TAUB2,0.)
            TAUB(L)  = SQRT(TAUB2)
            USTAR(L) = SQRT(TAUB(L))
          ELSE
            TAUBC = QQ(L,0)/CTURB2
            TAUB(L)  = TAUBC
            USTAR(L) = SQRT(TAUB(L))
          ENDIF
        ENDDO
      ELSE
        ! *** WITHOUT WAVE INTERACTIONS
        DO LP=LF,LL
          L = LSED(LP)
          TAUBC = QQ(L,0)/CTURB2
          TAUB(L)  = TAUBC
          USTAR(L) = SQRT(TAUB(L))
        ENDDO
      ENDIF
    ENDIF    ! *** END OF LSEDZLJ BYPASS
    
    ! *** COMPUTE CELL CENTERED VELOCITY, WITH BOUNDARY AND CORNER CORRECTIONS
    IF( (ICALC_BL > 0 .AND. .NOT. LSEDZLJ) .OR. (ICALC_BL > 0 .AND. LSEDZLJ) )THEN
      DO LP=LF,LL
        L = LSED(LP)
        LN = LNC(L)
        LE = LEC(L)
        UCELLCTR(L) = 0.5*( RSSBCW(L)*WCORWST(L)*U(L,KSZ(L)) + RSSBCE(L)*WCOREST(L)*U(LE,KSZ(L)) )
        VCELLCTR(L) = 0.5*( RSSBCS(L)*WCORSTH(L)*V(L,KSZ(L)) + RSSBCN(L)*WCORNTH(L)*V(LN,KSZ(L)) )
        QCELLCTR(L) = SQRT(UCELLCTR(L)*UCELLCTR(L)+VCELLCTR(L)*VCELLCTR(L))
      
        ! *** CALCULATE UNIT VECTOR
        IF( QCELLCTR(L) > 0.0 )THEN
          UCELLCTR(L)  = UCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, X COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          VCELLCTR(L)  = VCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, Y COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          CBEDTOTAL(L) = TAUB(L)/(QCELLCTR(L)*QCELLCTR(L))
        ELSE
          UCELLCTR(L)  = 0.
          VCELLCTR(L)  = 0.
          CBEDTOTAL(L) = 0.0
        ENDIF
      ENDDO
    ENDIF
    
    ! *** INITIALIZE GRAIN SHEARS WITH TOTAL
    DO LP=LF,LL
      L = LSED(LP)
      TAUBSED(L)  = TAUB(L)
      TAUBSND(L)  = TAUB(L)
      USTARSED(L) = USTAR(L)
      USTARSND(L) = USTAR(L)
    ENDDO

    ! *** CALCULATE CELL CENTER DEPTH DEVIATION FROM UNIFORM FLOW                                                           
    IF( .NOT. LSEDZLJ )THEN
      IF( ISBSDFUF >= 1 )THEN
        DO LP=LF,LL
          L = LSED(LP)
          LN = LNC(L)
          HDFUFXX = 0.5*(RSSBCW(L)*WCORWST(L)*HDFUFX(L)+RSSBCE(L)*WCOREST(L)*HDFUFX(LEC(L)))
          IF( HDFUFXX <= 0.0) HDFUFXX=1.0
          HDFUFYY = 0.5*(RSSBCS(L)*WCORSTH(L)*HDFUFY(L)+RSSBCN(L)*WCORNTH(L)*HDFUFY(LN ))
          IF( HDFUFYY <= 0.0) HDFUFYY=1.0
          HDFUF(L) = SQRT(HDFUFXX*HDFUFXX+HDFUFYY*HDFUFYY)
        ENDDO
      ENDIF
                            
      IF( KC > 1 )THEN
        ! *** CALCULATE CELL CENTER FLOW CORRECTION, DEPTH AVG/BOTTOM
        TMPEXP=16./7.
        DO LP=LF,LL
          L = LSED(LP)
          IF( QCELLCTR(L) > 0.0 )THEN
            LN = LNC(L)
            UTMP = 0.5*(RSSBCW(L)*WCORWST(L)*UHDYE(L)+RSSBCE(L)*WCOREST(L)*UHDYE(LEC(L)))/(DYP(L)*HP(L))
            VTMP = 0.5*(RSSBCS(L)*WCORSTH(L)*VHDXE(L)+RSSBCN(L)*WCORNTH(L)*VHDXE(LN ))/(DXP(L)*HP(L))
            QCELLCTRA = SQRT(UTMP*UTMP+VTMP*VTMP)
            QCELLAD1SQ(L) = (QCELLCTRA/QCELLCTR(L))**2
            QCELLAD1ZZ(L) = (QCELLCTRA/QCELLCTR(L))**TMPEXP
          ELSE
            QCELLAD1SQ(L) = 1.
            QCELLAD1ZZ(L) = 1.
          ENDIF
        ENDDO
      ENDIF

      ! *****************************************************************************************
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      IF( ISBEDSTR == 1 .OR. ISBEDSTR == 2 )THEN
        TMPEXP=2./7.
        TMPVAL=1./(COEFTSBL*VISMUDST)
        DO LP=LF,LL
          L = LSED(LP)
          TVAR3S(L) = TMPVAL*HP(L)*QCELLCTR(L)
          TVAR3N(L) = SEDDIAGS(L,KBT(L))*HPI(L) 
          HGDH(L) = 0.0
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          IF( TVAR3S(L) > 0.0 )THEN
            TVAR3S(L) = 1./TVAR3S(L)
          ENDIF
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          ! *** PREVIOUS TIME STEP'S SEDIMENT BED SHEAR
          IF( TAUBSEDS(L) > 0.0 )THEN
            TVAR3E(L) = QCELLCTR(L)/SQRT(TAUBSEDS(L))
          ELSE
            TVAR3E(L) = 0.0
          ENDIF
        ENDDO
                              
        ! ***  TVAR3E(L) = VELOCITY_MAGNITUDE/COHESIVE_BEDSTRESS                                                                  
        ! ***  TVAR3N(L) = D50/DEPTH                                                                                              
        DO LP=LF,LL
          L = LSED(LP)
          TVAR3W(L) = TVAR3S(L)**0.333333
          TVAR3E(L) = TVAR3E(L)**0.333333
          TVAR3N(L) = TVAR3N(L)**0.333333
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          TVAR3S(L) = TVAR3S(L)**TMPEXP
        ENDDO
                              
        ! *** 
        DO LP=LF,LL
          L = LSED(LP)
          IF( CBEDTOTAL(L) > 0.0 )THEN
            K = KBT(L)
            HGDH(L) = 0.014*HDFUF(L)*QCELLAD1SQ(L)*(FRACCOH(L,K)*TVAR3W(L)*TVAR3E(L) &
                                                    +FRACNON(L,K)*TVAR3N(L))/CBEDTOTAL(L)
          ENDIF
            
          HGDH(L) = HGDH(L)**0.75
          HGDH(L) = MIN(HGDH(L),1.0)
            
          ! *** CONVERT HGDH TO 1/HGDH  IE HDHG                                                                                     
          IF( HGDH(L) > 0.) HGDH(L) = 1./HGDH(L)
        ENDDO

        DO LP=LF,LL
          L = LSED(LP)
          IF( TAUB(L) > 0.0 )THEN
            TAUBSED(L) = HGDH(L)**TMPEXP
            TAUBSND(L) = HGDH(L)**0.333333
            TVAR3W(L)  = QCELLCTR(L)*QCELLCTR(L)
          ENDIF
        ENDDO

        DO LP=LF,LL
          L = LSED(LP)
          IF( TAUB(L) > 0.0 )THEN
            TAUBSED(L) = 0.026*TVAR3S(L)*TAUBSED(L)*TVAR3W(L)*QCELLAD1ZZ(L)
            TAUBSND(L) = 0.014*TVAR3N(L)*TAUBSND(L)*TVAR3W(L)*QCELLAD1SQ(L)
            ! *** SAVE FOR NEXT ITERATION
            TAUBSEDS(L) = TAUBSED(L)
          ENDIF
        ENDDO
                              
        !----------------------------------------------------------------------C                                                
        ! *** IF ISBEDSTR=2, APPLY WEIGHTED AVERAGE TO BOTH SED AND SND                                                         
        IF( ISBEDSTR == 2 .OR. N == 1 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            K = KBT(L)
            TVAR3E(L) = FRACCOH(L,K)*TAUBSED(L)+FRACNON(L,K)*TAUBSND(L)
          ENDDO
          DO LP=LF,LL
            L = LSED(LP)
            TAUBSED(L) = TVAR3E(L)
            TAUBSND(L) = TVAR3E(L)
          ENDDO
        ENDIF
        
        ! *** ADJUST WAVE CELLS
        IF( ISWAVE > 0 )THEN
          DO LP=LF,LL
            L = LSED(LP)
            IF( LWVMASK(L) )THEN
              ! *** Shear due to Current Only
              SHEAR = CTAUC(L)
              IF( SHEAR > 0. )THEN
                ! *** COHESIVE FRACTION
                RAT1O = MIN(TAUBSED(L)/SHEAR,2.)
                RAT1O = MAX(RAT1O,0.2)
                TAUBSED(L) = TAUB(L)*RAT1O
                !IF(L == 16459)PRINT *,L,RAT1O,TAUB(L),TAUBSED(L)

                ! *** COHESIVE FRACTION
                RAT1O = MIN(TAUBSND(L)/SHEAR,2.)
                RAT1O = MAX(RAT1O,0.2)
                TAUBSND(L) = TAUB(L)*RAT1O
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      
        ! *** COMPUTE SHEAR VELOCITIES
        DO LP=LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))
          USTARSND(L) = SQRT(TAUBSND(L))
        ENDDO

      ENDIF  ! *** ENDIF ON GRAIN STRESS PARTITIONING FOR ISBEDSTR == 1 .OR. ISBEDSTR == 2
                            
      !----------------------------------------------------------------------C                                                
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      ! *** INDEPENDENTLY SET GRAIN STRESS                                                                                    
      IF( ISBEDSTR == 3 )THEN
        DO LP=LF,LL
          L = LSED(LP)
          HDZBR=HP(L)/ZBRSED(L)
          TVAR3E(L) = 0.16/( (LOG(HDZBR)-1.)**2 )
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          TAUBSED(L) = TVAR3E(L)*QCELLCTR(L)*QCELLCTR(L)
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          TAUBSND(L) = TAUBSED(L)
        ENDDO
        DO LP=LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))
          USTARSND(L) = SQRT(TAUBSND(L))
        ENDDO
                              
      ENDIF
    
      !----------------------------------------------------------------------C                                                
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      !     USING WEIGHTED ROUGHNESS AND LOG RESISTANCE LAW                                                                   
      IF( ISBEDSTR == 4 )THEN
        !     RECALCUATE ZOTOTAL AND CALCULATE ZOGRAIN                                                                          
        DO LP=LF,LL
          L = LSED(LP)
          IF( CBEDTOTAL(L) > 0.0 )THEN                                                                                          
            TMP=EXP(1.+0.4/SQRT(CBEDTOTAL(L)))                                                                                 
            ZOTOTAL(L) = HP(L)/TMP
          ELSE
            ZOTOTAL(L) = ZBR(L)                                                                                                  
          ENDIF                                                                                                                
          ZOGRAIN(L) = FRACNON(L,KBT(L))*SEDDIAGS(L,KBT(L))/30.
        ENDDO                                                                                                                  

        DO LP=LF,LL
          L = LSED(LP)
          IF( TAUBSED(L) > 0.0 )THEN                                                                                            
            ZOGRAINCOH = 0.041*VISMUDST/SQRT(TAUBSED(L))
            ZOGRAIN(L) = ZOGRAIN(L)+FRACCOH(L,KBT(L))*ZOGRAINCOH
            ZOGRAIN(L) = MAX(ZOGRAIN(L),ZOGRAINCOH)
          ENDIF
        ENDDO                                                                                                                  

        !     ITERATE RVAL = SQRT(HG/H)                                                                                         
        DO LP=LF,LL
          L = LSED(LP)
          TMPTOP=LOG(HP(L)/ZOTOTAL(L))-1.
          TMPBOT=LOG(HP(L)/ZOGRAIN(L))-1.
          RVAL=1.                                                                                                              
          DO KK=1,100                                                                                                          
            RVAL = TMPTOP/(LOG(RVAL*RVAL)+TMPBOT)                                                                                
          ENDDO
          RVAL=MIN(RVAL,1.0)
          TAUBSED(L) = RVAL*RVAL*TAUB(L)
          TAUBSND(L) = RVAL*RVAL*TAUB(L)                                                                                         
        ENDDO                                                                                                                  

        DO LP=LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))                                                                                         
          USTARSND(L) = SQRT(TAUBSND(L))                                                                                         
        ENDDO                                                                                                                  

      ENDIF
                              
      !----------------------------------------------------------------------C                                                
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      !     USING WEIGHTED ROUGHNESS AND POWER RESISTANCE LAW                                                                 
      IF( ISBEDSTR == 5 )THEN
                              
        ! *** RECALCUATE ZOTOTAL AND CALCULATE ZOGRAIN                                                                          
        DO LP=LF,LL
          L = LSED(LP)
          IF( CBEDTOTAL(L) > 0.0 )THEN                                                                                          
            TMP=CBEDTOTAL(L)/0.04736                                                                                           
            ZOTOTAL(L) = HP(L)*(TMP**3)
          ELSE
            ZOTOTAL(L) = ZBR(L)                                                                                                  
          ENDIF                                                                                                                
          ZOGRAIN(L) = FRACNON(L,KBT(L))*SEDDIAGS(L,KBT(L))/30.
        ENDDO                                                                                                                  
                              
        DO LP=LF,LL
          L = LSED(LP)
          IF( TAUBSED(L) > 0.0 )THEN                                                                                            
            ZOGRAINCOH = 0.041*VISMUDST/SQRT(TAUBSED(L))
            ZOGRAIN(L) = ZOGRAIN(L)+FRACCOH(L,KBT(L))*ZOGRAINCOH
            ZOGRAIN(L) = MAX(ZOGRAIN(L),ZOGRAINCOH)
          ENDIF                                                                                                                
        ENDDO                                                                                                                  
                              
        !     CALCULATE GRAIN STRESS DIRECTLY                                                                                   
        DO LP=LF,LL
          L = LSED(LP)
          TMP = (ZOGRAIN(L)/ZOTOTAL(L))**0.25
          TMP = MIN(TMP,1.0)
          TAUBSED(L) = TMP*TAUB(L)
          TAUBSND(L) = TMP*TAUB(L)                                                                                               
        ENDDO                                                                                                                  
                              
        DO LP=LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))                                                                                         
          USTARSND(L) = SQRT(TAUBSND(L))                                                                                         
        ENDDO                                                                                                                  
                              
      ENDIF
                              
      ! ***********************************************************************C                                                
      DO LP=LF,LL
        L = LSED(LP)
        HBEDA(L) = 0.0
        DO K=1,KBT(L)
          HBEDA(L) = HBEDA(L)+HBED(L,K)
        END DO
      ENDDO
    ENDIF  ! *** END OF SEDZLJ BYPASS
  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO

  ! *** COMPUTE CELL CENTERED VELOCITY, WITH BOUNDARY AND CORNER CORRECTIONS
  IF( (ICALC_BL > 0 .AND. .NOT. LSEDZLJ) .OR. (ICALC_BL > 0 .AND. LSEDZLJ) )THEN
    ! *** COMPUTE CELL CENTERED TRANSPORT VECTORS FOR ONE ROW OF CELLS AROUND THE ACTIVE REGION
    IF( ISBEDMAP > 0 )THEN
      DO LP=1,NSEDPLUS
        L = LSEDPLUS(LP)
        LN = LNC(L)
        LE = LEC(L)
        UCELLCTR(L) = 0.5*( RSSBCW(L)*WCORWST(L)*U(L,KSZ(L)) + RSSBCE(L)*WCOREST(L)*U(LE,KSZ(L)) )
        VCELLCTR(L) = 0.5*( RSSBCS(L)*WCORSTH(L)*V(L,KSZ(L)) + RSSBCN(L)*WCORNTH(L)*V(LN,KSZ(L)) )
        QCELLCTR(L) = SQRT(UCELLCTR(L)*UCELLCTR(L)+VCELLCTR(L)*VCELLCTR(L))
      
        ! *** CALCULATE UNIT VECTOR
        IF( QCELLCTR(L) > 0.0 )THEN
          UCELLCTR(L)  = UCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, X COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          VCELLCTR(L)  = VCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, Y COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          CBEDTOTAL(L) = TAUB(L)/(QCELLCTR(L)*QCELLCTR(L))
        ELSE
          UCELLCTR(L)  = 0.
          VCELLCTR(L)  = 0.
          CBEDTOTAL(L) = 0.
        ENDIF
      ENDDO
    ENDIF
  ENDIF
    
  ! ***********************************************************************C                                                
  ! *** CALCULATE COHESIVE SEDIMENT SETTLING, DEPOSITION AND RESUSPENSION                                                 
  IF( ISTRAN(6) >= 1 )THEN
    IF( LSEDZLJ )THEN
      CALL SEDZLJ_MAIN
      GOTO 1000  ! *** BYPASS ORIGINAL BED/WATER INTERFACE CALCULATION
    ELSE
      CALL CALSED
    ENDIF
  ENDIF
                            
  ! ***********************************************************************C                                                
  ! *** CALCULATE NONCOHESIVE SEDIMENT BEDLOAD TRANSPORT, SETTLING,                                                       
  ! *** DEPOSITION AND RESUSPENSION                                                                                       
  IF( ISTRAN(7) >= 1 ) CALL CALSND
                          
  ! ***********************************************************************C                                                
  ! *** CALCULATE BANK EROSION AND ADJUST SEDIMENT AND WATER VOLUME FLUXES                                                                                                            
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    IF( ISBKERO >= 1 )THEN
      CALL BANKEROSED
      DO NP=1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        QSBDTOP(LBANK) = QSBDTOP(LBANK)+QSBDTOPBEBKB(NP)
        QWBDTOP(LBANK) = QWBDTOP(LBANK)+QWBDTOPBEBKB(NP)
        QSBDTOP(LCHAN) = QSBDTOP(LCHAN)+QSBDTOPBECHB(NP)
        QWBDTOP(LCHAN) = QWBDTOP(LCHAN)+QWBDTOPBECHB(NP)
      ENDDO
    ENDIF
  ENDIF
                            
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,LN,KTOPTP,KTOPM1) &
  !$OMP             PRIVATE(DSEDGMM,DSEDGMMI,QSWPOS,QSWNEG,TMPVAL,TTHICK,SXD,SXD1)
  DO ND=1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = MIN(LF+LDMSED-1,LASED)

    ! ***********************************************************************C                                                
    ! *** CALCULATE PARENT TO ACTIVE LAYER SEDIMENT FLUX                                                                    
    DO LP=LF,LL
      L = LSED(LP)
      QWATPA(L) = 0.0
      QSSDPA(L) = 0.0
    ENDDO

    IF( ISNDAL == 2 )THEN
      IF( IALTYP == 0 )THEN
        ! *** CONSTANT ACTIVE ARMOR LAYER THICKNESS
        DO NS=1,NSED
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          DO LP=LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            IF( KTOPM1 > 0 )THEN
              QSWPOS = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPM1) )
              QSWNEG = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPTP) )
              SEDFPA(L,NS) = DSEDGMMI*( VFRBED(L,KTOPM1,NS)*MAX(QSWPOS,0.) + VFRBED(L,KTOPTP,NS)*MIN(QSWNEG,0.) )     ! *** FLUX OF SND BETWEEN ACTIVE/PARENT LAYERS (G/M2/S)
              QSSDPA(L) = QSSDPA(L) + DSEDGMM*SEDFPA(L,NS)
              QWATPA(L) = QWATPA(L) + DSEDGMM*( VDRBED(L,KTOPM1)*MAX(SEDFPA(L,NS),0.) + VDRBED(L,KTOPTP)*MIN(SEDFPA(L,NS),0.))
              SEDB(L,KTOPTP,NS) = SEDB(L,KTOPTP,NS) + DTSED*SEDFPA(L,NS)
              SEDB(L,KTOPM1,NS) = SEDB(L,KTOPM1,NS) - DTSED*SEDFPA(L,NS)
            ENDIF
          ENDDO
        ENDDO
        DO NX=1,NSND
          NS=NX+NSED
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          DO LP=LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            IF( KTOPM1 > 0 )THEN
              QSWPOS = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPM1) )
              QSWNEG = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPTP) )
              SNDFPA(L,NX) = DSEDGMMI*( VFRBED(L,KTOPM1,NS)*MAX(QSWPOS,0.) + VFRBED(L,KTOPTP,NS)*MIN(QSWNEG,0.) )
              QSSDPA(L) = QSSDPA(L) + DSEDGMM*SNDFPA(L,NX)
              QWATPA(L) = QWATPA(L) + DSEDGMM*( VDRBED(L,KTOPM1)*MAX(SNDFPA(L,NX),0.) + VDRBED(L,KTOPTP)*MIN(SNDFPA(L,NX),0.) )
              SNDB(L,KTOPTP,NX) = SNDB(L,KTOPTP,NX) + DTSED*SNDFPA(L,NX)
              SNDB(L,KTOPM1,NX) = SNDB(L,KTOPM1,NX) - DTSED*SNDFPA(L,NX)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        ! *** CONSTANT ACTIVE ARMOR LAYER TOTAL SEDIMENT MASS                                                                   
        DO NS=1,NSED
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          DO LP=LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            IF( KTOPM1 > 0 )THEN
              SEDFPA(L,NS) = VFRBED(L,KTOPM1,NS)    *MAX(QSBDTOP(L),0.)  +VFRBED(L,KTOPTP,NS)*MIN(QSBDTOP(L),0.)
              QSSDPA(L) = QSSDPA(L)+SEDFPA(L,NS)
              QWATPA(L) = QWATPA(L)+VDRBED(L,KTOPM1)*MAX(SEDFPA(L,NS),0.)+VDRBED(L,KTOPTP)*MIN(SEDFPA(L,NS),0.)
              SEDFPA(L,NS) = DSEDGMMI*SEDFPA(L,NS)
              SEDB(L,KTOPTP,NS) = SEDB(L,KTOPTP,NS) + DTSED*SEDFPA(L,NS)
              SEDB(L,KTOPM1,NS) = SEDB(L,KTOPM1,NS) - DTSED*SEDFPA(L,NS)
            ENDIF
          ENDDO
        ENDDO
        DO NX=1,NSND
          NS=NX+NSED
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          DO LP=LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            IF( KTOPM1 > 0 )THEN
              SNDFPA(L,NX) = VFRBED(L,KTOPM1,NS)    *MAX(QSBDTOP(L),0.)  +VFRBED(L,KTOPTP,NS)*MIN(QSBDTOP(L),0.)
              QSSDPA(L) = QSSDPA(L)+SNDFPA(L,NX)
              QWATPA(L) = QWATPA(L)+VDRBED(L,KTOPM1)*MAX(SNDFPA(L,NX),0.)+VDRBED(L,KTOPTP)*MIN(SNDFPA(L,NX),0.)
              SNDFPA(L,NX) = DSEDGMMI*SNDFPA(L,NX)
              SNDB(L,KTOPTP,NX) = SNDB(L,KTOPTP,NX) + DTSED*SNDFPA(L,NX)
              SNDB(L,KTOPM1,NX) = SNDB(L,KTOPM1,NX) - DTSED*SNDFPA(L,NX)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
                            
    ! *****************************************************************************                                           
    ! *** UPDATE TOP BED LAYER THICKNESS AND VOID RATIO                                                                     
    ! *** FOR DEPOSITION-RESUSPENSION STEP                                                                                  
    DO LP=LF,LL
      L = LSED(LP)
      K = KBT(L)
      HBED1(L,K) = HBED(L,K)
      VDRBED1(L,K) = VDRBED(L,K)

      IF( DEBUG )THEN
        ! *** PMC - FOR DEBUGGING
        SXD  = SUM(SEDB(L,K,:))  + SUM(SNDB(L,K,:))
        SXD1 = SUM(SEDB1(L,K,:)) + SUM(SNDB1(L,K,:))
        HBED(L,K) = HBED(L,K) - DTSED*( QSBDTOP(L) + QWBDTOP(L) ) + DTSED*( QSSDPA(L) + QWATPA(L) )
      
        IF( SIGN(1.0,(HBED(L,K)-HBED1(L,K))) /= SIGN(1.0,(SXD-SXD1)) )THEN 
          TTHICK = 0.
          DO NS=1,NSED  
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SEDB(L,K,NS)*DSEDGMM
          ENDDO  
          DO NX=1,NSND  
            NS = NSED+NX  
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SNDB(L,K,NX)*DSEDGMM
          ENDDO 
          TTHICK = TTHICK/(1.0 - PORBED(L,K))
          IF( ABS(TTHICK-HBED(L,K)) > 1E-3 )THEN
             PRINT '(F10.4,2I5,2F10.5,4E14.6,F10.5)', TIMEDAY,L,K,HBED(L,K),HBED1(L,K),SXD,SXD1, DTSED*QSBDTOP(L), DTSED*QWBDTOP(L), TTHICK
          ENDIF   
        ENDIF
      ENDIF

      IF( K == 1 )THEN
        IF( HBED(L,K) < MINTHICK )THEN
          ! *** ZERO NEGATIVE THICKNESSES
          HBED(L,K) = 0.0
          VDRBED(L,K) = SNDVDRD
          PORBED(L,K) = BEDPORC
          STDOCB(L,K) = 0.0
          STPOCB(L,K) = 0.0
          
          SEDB(L,K,:) = 0.0
          SNDB(L,K,:) = 0.0
        ENDIF
      ENDIF
      
      ! *** Update void ratio with the depositing sediment mass, using the last layer thickness.
      ! *** HBED then gets updated in CALBED
      TMPVAL = HBED1(L,K)/(1. + VDRBED1(L,K))
      TMPVAL = TMPVAL - DTSED*( QSBDTOP(L) - QSSDPA(L) )
      IF( TMPVAL > 0.0 )THEN
        IF( HBED(L,K) > 0.0 )THEN        
          VDRBED(L,K) = (HBED(L,K)/TMPVAL)-1.
          IF( K == 1 )THEN
            ! *** LIMIT VOID RATIOS TO 0.01 >= N <= 99
            IF( VDRBED(L,K) < 0.01 .OR. VDRBED(L,K) > 99. )THEN
              VDRBED(L,K) = SNDVDRD  
            ENDIF
          ENDIF
        ELSE
          ! *** Bed thickness is zero.  Initialize to depositing sediment mass
          HBED(L,K) = -DTSED*( QSBDTOP(L) - QSSDPA(L) )/(1. - PORBED(L,K))
          VDRBED(L,K) = PORBED(L,K)/(1. - PORBED(L,K))
        ENDIF
      ELSE
        VDRBED(L,K) = SNDVDRD          ! *** Assign default VR if layer empty or completely eroded away
      ENDIF
    ENDDO
                            
    ! *** UPDATE PARENT LAYER BED LAYER THICKNESS AND VOID RATIO                                                            
    ! *** FOR DEPOSITION-RESUSPENSION STEP WHEN USING ACTIVE-PARENT SCHEME
    IF( ISNDAL == 2 )THEN
      DO LP=LF,LL
        L = LSED(LP)
        K = KBT(L)-1
        IF( K > 0 )THEN
          HBED1(L,K) = HBED(L,K)
          VDRBED1(L,K) = VDRBED(L,K)
          HBED(L,K) = HBED(L,K) - DTSED*(QSSDPA(L)+QWATPA(L))

          IF( K == 1 )THEN
            IF( HBED(L,K) < 0.0 )THEN
              ! *** ZERO NEGATIVE THICKNESSES
              HBED(L,K) = 0.0
              VDRBED(L,K) = SNDVDRD
              PORBED(L,K) = BEDPORC
              STDOCB(L,K) = 0.0
              STPOCB(L,K) = 0.0
            ENDIF
          ENDIF
          TMPVAL = HBED1(L,K)/(1.+VDRBED1(L,K))
          TMPVAL = TMPVAL - DTSED*QSSDPA(L)
          IF( TMPVAL > 0.0 )THEN
            VDRBED(L,K) = (HBED(L,K)/TMPVAL)-1.
          ELSE
            VDRBED(L,K) = SNDVDRD
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO
    
  ! **********************************************************************
  1000 CONTINUE  ! *** JUMP DIRECTLY HERE FROM SEDZLJ_MAIN CALL
                           
  ! ***********************************************************************C                                                
  ! *** CALCULATE TOXIC SETTLING, DEPOSITION AND RESUSPENSION                                                             
  IF( ISTRAN(5) > 0 )THEN
    TTDSX = DSTIME(0)
    CALL CALTOX
    TSSTX = TSSTX + (DSTIME(0)-TTDSX)
  ENDIF                      
                          
  IF( .NOT. LSEDZLJ )THEN
    ! ***********************************************************************C                                                
    ! *** UPDATE TOP BED LAYER THICKNESS AND VOID RATIO                                                                     
    ! *** FOR DEPOSITION-RESUSPENSION STEP                                                                                  
    ! *** PRESENTLY ACTIVE BEFORE THE WATER COLUMN-BED TOXICS EXCHANGE                                                      
    ! *** CHECK PLACEMENT THERE AND HERE FOR                                                                                
    !        VDRTMP=(TMPVAL/HBED(L,K))-1.                                                                                   
    !        VDRBED(L,K) = VDRTMP                                                                                             

    ! ***********************************************************************C                                                
    ! *** UPDATE SEDIMENT BED LAYERING                                                                                      
    IF( KB > 1 ) CALL CALBLAY
                          
    ! ***********************************************************************C                                                
    ! *** RESET SEDIMENT VOLUME FRACTIONS                                                                                   
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX)
    DO ND=1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      DO K=1,KB
        DO LP=LF,LL
          L = LSED(LP)
          BEDLINIT(L,K) = 0.
        ENDDO
      ENDDO
                            
      DO NX=1,NSED+NSND
        DO K=1,KB
          DO LP=LF,LL
            L = LSED(LP)
            VFRBED(L,K,NX) = 0.
          ENDDO
        ENDDO
      ENDDO                                                                                                                  

      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
                            
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
                            
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
                            
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
                            
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              IF( BEDLINIT(L,K) > 0.0 )THEN
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              ELSE
                VFRBED(L,K,NS) = 0.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
                            
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NSED+NX
          DO K=1,KB
            DO LP=LF,LL
              L = LSED(LP)
              IF( BEDLINIT(L,K) > 0.0 )THEN
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              ELSE
                VFRBED(L,K,NS) = 0.0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
                            
      ! *** UPDATE VOLUME OF COHESIVES AND NON-COHESIVES
      DO K=1,KB                                                                                                              
        DO LP=LF,LL
          L = LSED(LP)
          FRACCOH(L,K) = 0.0                                                                                                     
          FRACNON(L,K) = 0.0
        ENDDO
      ENDDO
                            
      DO NS=1,NSED
        DO K=1,KB                                                                                                              
          DO LP=LF,LL
            L = LSED(LP)
            IF( K <= KBT(L) )THEN                                                                                                  
              FRACCOH(L,K) = FRACCOH(L,K)+VFRBED(L,K,NS)                                                                           
            ENDIF
          ENDDO
        ENDDO
      ENDDO
                            
      DO NX=1,NSND
        NS=NX+NSED                                                                                                             
        DO K=1,KB                                                                                                              
          DO LP=LF,LL
            L = LSED(LP)
            IF( K <= KBT(L) )THEN                                                                                                  
              FRACNON(L,K) = FRACNON(L,K)+VFRBED(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
                            
      ! *** CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS                                                                    
      DO K=1,KB
        DO LP=LF,LL
          L = LSED(LP)
          IF( K <= KBT(L) )THEN                                                                                                  
            VDRBEDSND(L,K) = SNDVDRD                                                                                             
            VDRBEDSED(L,K) = 0.0
            IF( FRACCOH(L,K) > 0.0 )THEN
              VDRBEDSED(L,K) = ( (FRACCOH(L,K)+FRACNON(L,K))*VDRBED(L,K)-FRACNON(L,K)*SNDVDRD )/FRACCOH(L,K)
            ENDIF
          ELSE
            VDRBEDSND(L,K) = 0.0                                                                                                 
            VDRBEDSED(L,K) = 0.0                                                                                                 
          ENDIF                                                                                                                
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
  ENDIF  ! *** END OF SEDZLJ BYPASS
  
  ! ***********************************************************************C                                                
  ! *** UPDATE SEDIMENT BED PHYSICAL PROPERTIES  
  IF( LSEDZLJ )THEN
    ! *** ALL SEDZLJ VARIABLES NEEDED BY EFDC ARE UPDATED IN S_SEDZLJ.F90
  ELSE
    IF( IBMECH == 9 )THEN
       CALL CALBED9
    ELSE
       CALL CALBED
    ENDIF
  ENDIF
  
  ! ***********************************************************************C                                                
  ! *** CHANGE BED MORPHOLOGY.  NOW INCLUDING SEDZLJ OPTION (2017-01)
  IF( IMORPH > 0 )THEN
    ITMP=0

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K)
    DO ND=1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      DO LP=LF,LL
        L = LSED(LP)
        BELV1(L) = BELV(L)
        P1(L)    = P(L)
        HBEDA(L) = 0.0
      ENDDO
      DO K=1,KB
        DO LP=LF,LL
          L = LSED(LP)
          HBEDA(L) = HBEDA(L) + HBED(L,K)
        END DO
      ENDDO 
      DO LP=LF,LL
        
        L = LSED(LP)
        BELV(L) = ZELBEDA(L) + HBEDA(L)
        P(L) =  ( HP(L) + BELV(L) )*G
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
    
    ! *** ADJUST CONCENTRATIONS OF TRANSPORT VARIABLES IN RESPONSE TO CHANGE IN BED MORPHOLOGY
    ! *** 2017-01 - PMC REMOVED CONCENTRATION ADJUSTMENT SINCE NEW APPROACH MAINTAINS WATER BALANCE
  ENDIF
  
  ! ***********************************************************************C                                                
  ! *** TOXICS CALCULATIONS
  IF( ISTRAN(5) > 0 )THEN
    ! *** POREWATER ADVECTION AND DIFFUSION OF TOXICS                                                                       
    TTDSX=DSTIME(0)
    CALL CALTOXB
                          
    ! *** TOXIC CONTAMINANT REACTIONS
    IF( NKINETICS > 0 ) CALL CALTOX_KINETICS
                      
    ! *** TOXIC CONTAMINANT OUTPUT TO FOOD CHAIN MODEL                                                                      
    IF( ISFDCH >= 1 ) CALL FOODCHAIN(0)

    TSSTX = TSSTX + (DSTIME(0)-TTDSX)
  ENDIF

  if( nactiveships > 0 )then
    prop_ero(:,:) = 0.0                       ! *** Zero erosion due to propwash for next iteration
    if( icalc_bl > 0 ) prop_bld(:,:) = 0.0    ! *** Zero erosion due to propwash for next iteration
  endif
  
  ! ***********************************************************************C                                                
  
  TTSED = TTSED + (DSTIME(0) - TTDS) - TTWAIT

  RETURN

END

