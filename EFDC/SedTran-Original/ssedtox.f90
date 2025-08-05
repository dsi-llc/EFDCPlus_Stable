! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
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

  use GLOBAL
  use Allocate_Initialize      
  use RESTART_MODULE
  use Variables_Propwash
  
  implicit none

  integer :: IFILE
  integer :: NVAL, LP, L, K, NS, NX, KK, ND, LF, LL, NP, ipmc
  integer :: KTOPTP, KTOPM1, NT, LE, LN, ITMP, LBANK, LCHAN
  integer,save :: NKINETICS

  real pmc0, pmc1, pmc2, pmc3, pmc4
  real :: HDZBR, DSEDGMMI, MINTHICK
  real :: UTMP, VTMP, CURANG, SHEAR, RAT1O
  real :: TAUBC, TAUB2, CSEDRESB, HDFUFXX, HDFUFYY, QSWNEG
  real :: TMPEXP, QCELLCTRA, TMPVAL, QSWPOS, TMP
  real :: ZOGRAINCOH, TMPTOP, TMPBOT, RVAL, SXD1, SXD, TTHICK
  real, external :: CSEDTAUS, CSEDRESS, CSEDTAUB

  real,save,allocatable,dimension(:)   :: DELBED
  real,save,allocatable,dimension(:)   :: QCELLAD1SQ
  real,save,allocatable,dimension(:)   :: QCELLAD1ZZ
  real,save,allocatable,dimension(:,:) :: SEDGEOSTD
  real,save,allocatable,dimension(:,:) :: SEDDIAGS
  real,save,allocatable,dimension(:)   :: TAUBSEDS
  real,save,allocatable,dimension(:)   :: ZOGRAIN
  real,save,allocatable,dimension(:)   :: ZOTOTAL
  
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS,TTDSX        ! MODEL TIMING TEMPORARY VARIABLE
  
  integer(IK4),save,allocatable,dimension(:) :: LSEDPLUS
  integer(IK4),save                          :: NSEDPLUS
  
  if( .not. allocated(DELBED) )then
    call AllocateDSI( DELBED,     KBM,  0.0)
    call AllocateDSI( QCELLAD1SQ, LCM,  1.0)
    call AllocateDSI( QCELLAD1ZZ, LCM,  1.0)
    call AllocateDSI( SEDGEOSTD,  LCM,  KBM,  0.0)
    call AllocateDSI( SEDDIAGS,   LCM,  KBM,  0.0)
    call AllocateDSI( TAUBSEDS,   LCM,  0.0)
    call AllocateDSI( ZOGRAIN,    LCM,  0.0)
    call AllocateDSI( ZOTOTAL,    LCM,  0.0)
                                        
    call AllocateDSI( LSEDPLUS,   LCM,  0)
    
    NKINETICS = 0
    do NT = 1,NTOX
      NKINETICS = NKINETICS+SUM(ITOXKIN(1:3,NT))
    enddo
    
    NSEDPLUS = 0
    if( ISBEDMAP > 0 )then
      ! *** EAST CELL IS A HARD BOTTOM CELL
      do LP = 1,BEDEDGEE.NEDGE
        L = BEDEDGEE.LEDGE(LP)
        if( SUBO(LEC(L)) > 0. )then
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LEC(L)
        endif
      enddo

      ! *** WEST CELL IS A HARD BOTTOM CELL
      do LP = 1,BEDEDGEW.NEDGE
        L = BEDEDGEW.LEDGE(LP)
        if( SUBO(L) > 0. )then
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LWC(L)
        endif
      enddo    
    
      ! *** SOUTH CELL IS A HARD BOTTOM CELL
      do LP = 1,BEDEDGES.NEDGE
        L = BEDEDGES.LEDGE(LP)
        if( SVBO(L) > 0. )then
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LSC(L)
        endif
      enddo    

      ! *** NORTH CELL IS A HARD BOTTOM CELL
      do LP = 1,BEDEDGEN.NEDGE
        L = BEDEDGEN.LEDGE(LP)
        if( SVBO(LNC(L)) > 0. )then
          NSEDPLUS = NSEDPLUS + 1
          LSEDPLUS(NSEDPLUS) = LNC(L)
        endif
      enddo    
    endif
  endif

  TTDS = DSTIME(0)

  IFILE = -1

  ! *** S2TL AND S3TL ARE GLOBAL VARIABLES USED IN CALTOX AND CALSND
  if( IS2TL == 1 )then
    S3TL = 1.0
    S2TL = 0.0
  else
    if( ISTL /= 3 )then
      S3TL = 0.0
      S2TL = 1.0
    else
      S3TL = 1.0
      S2TL = 0.0
    endif
  endif

  ! *** ALL TIMING IS BASED ON THE SEDTIME
  DTSEDJ = REAL(DTSED,8)
  DELTI = 1._8/DTSEDJ

  SEDMDGM = SQRT(SEDMDMX*SEDMDMN)
  NVAL = MOD(N,2)
  FOURDPI = 4./PI
  MINTHICK = min(1E-6,DTSED*1E-6)           ! *** Minimum bed thickness
  
  ! ***********************************************************************C                                                
  ! *** IF N = 1 CALCULATE INITIAL SEDIMENT BED THICKNESS
                          
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,LN,LE,KK,NT) &
  !$OMP             PRIVATE(TAUBC,UTMP,VTMP,CURANG,TAUB2,HDFUFXX,HDFUFYY,TMPEXP,TMPVAL)  &
  !$OMP             PRIVATE(HDZBR,TMP,TMPTOP,TMPBOT,RVAL,ZOGRAINCOH,SHEAR,RAT1O,QCELLCTRA)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)
                            
    ! *** SET SEDIMENT VOLUME FRACTIONS                                                                                     
    if( .not. LSEDZLJ )then
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          BEDLINIT(L,K) = 0.
          BEDDINIT(L,K) = 0.
        enddo
      enddo
      do NX = 1,NSED+NSND
        do K = 1,KB
          do LP = LF,LL
            L = LSED(LP)
            VFRBED(L,K,NX) = 0.
            VFRBED1(L,K,NX) = 0.
          enddo
        enddo
      enddo
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS)  = SDEN(NS)*SEDB(L,K,NS)
              VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS)  = SDEN(NS)*SNDB(L,K,NX)
              VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)
            enddo
          enddo
        enddo
      endif
    
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
              BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
              BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              if( BEDLINIT(L,K) > 0.0 )then
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              else
                VFRBED(L,K,NS) = 0.0
              endif
              if( BEDDINIT(L,K) > 0.0 )then
                VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
              else
                VFRBED1(L,K,NS) = 0.0
              endif
            enddo
          enddo
        enddo
      endif
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              if( BEDLINIT(L,K) > 0.0 )then
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              else
                VFRBED(L,K,NS) = 0.0
              endif
              if( BEDDINIT(L,K) > 0.0 )then
                VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)
              else
                VFRBED1(L,K,NS) = 0.0
              endif
            enddo
          enddo
        enddo
      endif

      ! *** INITIALIZE FRACTIONS FOR THE CURRENT TIME STEP
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          FRACCOH(L,K) = 0.0
          FRACNON(L,K) = 0.0
        enddo
      enddo
      do NS = 1,NSED
        do K = 1,KB
          do LP = LF,LL
            L = LSED(LP)
            if( K <= KBT(L) )then
              FRACCOH(L,K) = FRACCOH(L,K)+VFRBED(L,K,NS)
            endif
          enddo
        enddo
      enddo
      do NX = 1,NSND
        NS = NSED + NX
        do K = 1,KB
          do LP = LF,LL
            L = LSED(LP)
            if( K <= KBT(L) )then
              FRACNON(L,K) = FRACNON(L,K)+VFRBED(L,K,NS)
            endif
          enddo
        enddo
      enddo

      ! *** INITIALIZE SEDIMENT BED/WATER COLUMN INTERFACE FLUXES FOR THE CURRENT TIME STEP
      do LP = LF,LL
        L = LSED(LP)
        QWBDTOP(L) = 0.   ! *** VOLUME OF WATER EXCHANGE DUE TO ENTRAPMENT/SCOUR     (M/S)
        QSBDTOP(L) = 0.   ! *** VOLUME OF SEDIMENT EXCHANGE                          (M/S)
      enddo

      ! *** SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES                                                         
      if( ISTRAN(6) >= 1 )then 
        if( IWRSP(1) >= 1 .and. IWRSP(1) < 99 )then
          ! *** COMPUTE CRITICAL SHEAR STRESS FOR RESUSPENSION
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)  ! PMC - REMOVED LOOP AND ONLY COMPUTE THE TOP LAYER
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1), VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),VDRBED(L,K),VDRBED(L,K),IWRSP(1))
          enddo
        endif
        if( IWRSPB(1) >= 1 )then
          ! *** COMPUTE CRITICAL SHEAR STRESS FOR BULK/MASS EROSION
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)  ! PMC - REMOVED LOOP AND ONLY COMPUTE THE TOP LAYER
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))
            WRSPB(L,K) = CSEDRESB
          enddo
        endif
      endif
    endif    ! *** END OF LSEDZLJ BYPASS
                              
    ! *** CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC                                                        
    ! *** CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG                                                        
    ! *** TO TOXB UNITS OF OF MG/M**2                                                                                       
                            
    ! *** save OLD VALUES                                                                                                   
    if( ISTRAN(6) >= 1 )then
      do NS = 1,NSED2
        do K = 1,KC
          do LP = LF,LL
            L = LSED(LP)
            SED2(L,K,NS) = SED(L,K,NS)
          enddo
        enddo
      enddo
    endif
    
    if( .not. LSEDZLJ )then
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          do K = 1,KC
            do LP = LF,LL
              L = LSED(LP)
              SNDS(L,K,NX) = SND(L,K,NX)
            enddo
          enddo
        enddo
      endif

      ! ***********************************************************************C                                                
      ! *** SET MEAN D50 AND D90                                                                                              
      if( ISTRAN(7) >= 1 )then
        ! *** SEDIMENT MASS
        do K = 1,KB
          do LP = LF,LL
            L = LSED(LP)
            SNDBT(L,K) = 0.
          enddo
        enddo
        
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NX)
            enddo
          enddo
        enddo

        if( ISBSDIAM > 1 )then
          do LP = LF,LL
            L = LSED(LP)
            SEDDIA90(L,KBT(L)) = 0.
            SEDGEOSTD(L,KBT(L)) = 0.
          enddo
        endif
        
        ! *** D50
        do LP = LF,LL
          L = LSED(LP)
          SEDDIA50(L,KBT(L)) = 0.
        enddo
        do NX = 1,NSND
          NS = NSED + NX
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)
            SEDDIA50(L,K) = SEDDIA50(L,K) + SNDB(L,K,NX)*LOG(SEDDIA(NS))
          enddo
        enddo
        do LP = LF,LL
          L = LSED(LP)
          K = KBT(L)
          if( SNDBT(L,K) > 0. )then
            SEDDIA50(L,K) = SEDDIA50(L,K)/SNDBT(L,K)
          endif
        enddo
                            
        if( ISBSDIAM > 1 )then
          do NX = 1,NSND
            NS = NSED + NX
            do LP = LF,LL
              L = LSED(LP)
              K = KBT(L)
              SEDGEOSTD(L,K) = SEDGEOSTD(L,K) + SNDB(L,K,NX)*(( LOG(SEDDIA(NS))-SEDDIA50(L,K) )**2)
            enddo
          enddo
                            
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)
            if( SNDBT(L,K) > 0. )then
              SEDGEOSTD(L,K) = SEDGEOSTD(L,K)/SNDBT(L,K)
            endif
          enddo
        endif
        
        do LP = LF,LL
          L = LSED(LP)
          K = KBT(L)
          ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
          if( SNDBT(L,K) > 5.0 )then  
            SEDDIA50(L,K)  = EXP(SEDDIA50(L,K))
          elseif( SNDBT(L,K) > 1E-12 )then
            ! *** MASS WEIGHTED D50
            SEDDIA50(L,K) = 0.
            do NX = 1,NSND
              NS = NSED + NX
              SEDDIA50(L,K) = SEDDIA50(L,K) + SNDB(L,K,NX)*SEDDIA(NS)
            enddo
            SEDDIA50(L,K)  = SEDDIA50(L,K)/SNDBT(L,K)
          else
            SEDDIA50(L,K)  = SEDDIA(NSED+1)
          endif
        enddo
                            
        if( ISBSDIAM > 1 )then
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)
            ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
            if( SNDBT(L,K) > 5.0 )then  
              SEDGEOSTD(L,K) = EXP(SEDGEOSTD(L,K))
            else
              SEDGEOSTD(L,K) = 1.0
            endif
          enddo
                            
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)
            SEDDIA90(L,K) = (SEDGEOSTD(L,K)**1.28)*SEDDIA50(L,K)
          enddo
        endif
        
        ! *** SET SEDIMENT GRAINSIZE TO use FOR SHEAR STRESS CALCS  
        if( ISBSDIAM == 0 )then
          do LP = LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = SEDDIA50(L,KBT(L))
          enddo
        elseif( ISBSDIAM == 1 )then
          do LP = LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = 2.*SEDDIA50(L,KBT(L))
          enddo
        elseif( ISBSDIAM == 2 )then
          do LP = LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = SEDDIA90(L,KBT(L))
          enddo
        elseif( ISBSDIAM == 3 )then
          do LP = LF,LL
            L = LSED(LP)
            SEDDIAGS(L,KBT(L)) = 2.*SEDDIA90(L,KBT(L))
          enddo
        endif                                                                                                                                        
      endif    ! *** END OF NON-COHESIVE GRAIN SIZE CALCS
    endif    ! *** END OF LSEDZLJ BYPASS
                            
    ! ***********************************************************************C                                                
    ! *** SET CELL CENTER BED STRESS FOR SEDIMENT RESUSPENSION AND DEPOSITION                                               
    if( .not. LSEDZLJ )then
      if( ISWAVE > 0 )then
        ! *** WITH WAVE INTERACTIONS
        do LP = LF,LL
          L = LSED(LP)
          if( LWVMASK(L) )then
            TAUBC = QQ(L,0)/CTURB2  ! *** CURRENT PLUS BOTTOM SHEAR FROM QQWV1
            UTMP = 0.5*STCUV(L)*(U(LEC(L),KSZ(L)) + U(L,KSZ(L)))+1.E-12
            VTMP = 0.5*STCUV(L)*(V(LNC(L),KSZ(L)) + V(L,KSZ(L)))
            CURANG = ATAN2(VTMP,UTMP)
            TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
            TAUB2 = max(TAUB2,0.)
            TAUB(L)  = SQRT(TAUB2)
            USTAR(L) = SQRT(TAUB(L))
          else
            TAUBC = QQ(L,0)/CTURB2
            TAUB(L)  = TAUBC
            USTAR(L) = SQRT(TAUB(L))
          endif
        enddo
      else
        ! *** WITHOUT WAVE INTERACTIONS
        do LP = LF,LL
          L = LSED(LP)
          TAUBC = QQ(L,0)/CTURB2
          TAUB(L)  = TAUBC
          USTAR(L) = SQRT(TAUB(L))
        enddo
      endif
    endif    ! *** END OF LSEDZLJ BYPASS
    
    ! *** COMPUTE CELL CENTERED VELOCITY, WITH BOUNDARY AND CORNER CORRECTIONS
    if( (ICALC_BL > 0 .and. .not. LSEDZLJ) .or. (ICALC_BL > 0 .and. LSEDZLJ) )then
      do LP = LF,LL
        L = LSED(LP)
        LN = LNC(L)
        LE = LEC(L)
        UCELLCTR(L) = 0.5*( RSSBCW(L)*WCORWST(L)*U(L,KSZ(L)) + RSSBCE(L)*WCOREST(L)*U(LE,KSZ(L)) )
        VCELLCTR(L) = 0.5*( RSSBCS(L)*WCORSTH(L)*V(L,KSZ(L)) + RSSBCN(L)*WCORNTH(L)*V(LN,KSZ(L)) )
        QCELLCTR(L) = SQRT(UCELLCTR(L)*UCELLCTR(L)+VCELLCTR(L)*VCELLCTR(L))
      
        ! *** CALCULATE UNIT VECTOR
        if( QCELLCTR(L) > 0.0 )then
          UCELLCTR(L)  = UCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, X COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          VCELLCTR(L)  = VCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, Y COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          CBEDTOTAL(L) = TAUB(L)/(QCELLCTR(L)*QCELLCTR(L))
        else
          UCELLCTR(L)  = 0.
          VCELLCTR(L)  = 0.
          CBEDTOTAL(L) = 0.0
        endif
      enddo
    endif
    
    ! *** INITIALIZE GRAIN SHEARS WITH TOTAL
    if( .not. LSEDZLJ )then
      do LP = LF,LL
        L = LSED(LP)
        TAUBSED(L)  = TAUB(L)
        TAUBSND(L)  = TAUB(L)
        USTARSED(L) = USTAR(L)
        USTARSND(L) = USTAR(L)
      enddo
    endif

    ! *** CALCULATE CELL CENTER DEPTH DEVIATION FROM UNIFORM FLOW                                                           
    if( .not. LSEDZLJ )then
      if( ISBSDFUF >= 1 )then
        do LP = LF,LL
          L = LSED(LP)
          LN = LNC(L)
          HDFUFXX = 0.5*(RSSBCW(L)*WCORWST(L)*HDFUFX(L)+RSSBCE(L)*WCOREST(L)*HDFUFX(LEC(L)))
          if( HDFUFXX <= 0.0) HDFUFXX = 1.0
          HDFUFYY = 0.5*(RSSBCS(L)*WCORSTH(L)*HDFUFY(L)+RSSBCN(L)*WCORNTH(L)*HDFUFY(LN ))
          if( HDFUFYY <= 0.0) HDFUFYY = 1.0
          HDFUF(L) = SQRT(HDFUFXX*HDFUFXX+HDFUFYY*HDFUFYY)
        enddo
      endif
                            
      if( KC > 1 )then
        ! *** CALCULATE CELL CENTER FLOW CORRECTION, DEPTH AVG/BOTTOM
        TMPEXP = 16./7.
        do LP = LF,LL
          L = LSED(LP)
          if( QCELLCTR(L) > 0.0 )then
            LN = LNC(L)
            UTMP = 0.5*(RSSBCW(L)*WCORWST(L)*UHDYE(L)+RSSBCE(L)*WCOREST(L)*UHDYE(LEC(L)))/(DYP(L)*HP(L))
            VTMP = 0.5*(RSSBCS(L)*WCORSTH(L)*VHDXE(L)+RSSBCN(L)*WCORNTH(L)*VHDXE(LN ))/(DXP(L)*HP(L))
            QCELLCTRA = SQRT(UTMP*UTMP+VTMP*VTMP)
            QCELLAD1SQ(L) = (QCELLCTRA/QCELLCTR(L))**2
            QCELLAD1ZZ(L) = (QCELLCTRA/QCELLCTR(L))**TMPEXP
          else
            QCELLAD1SQ(L) = 1.
            QCELLAD1ZZ(L) = 1.
          endif
        enddo
      endif

      ! *****************************************************************************************
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      if( ISBEDSTR == 1 .or. ISBEDSTR == 2 )then
        TMPEXP = 2./7.
        TMPVAL = 1./(COEFTSBL*VISMUDST)
        do LP = LF,LL
          L = LSED(LP)
          TVAR3S(L) = TMPVAL*HP(L)*QCELLCTR(L)
          TVAR3N(L) = SEDDIAGS(L,KBT(L))*HPI(L) 
          HGDH(L) = 0.0
        enddo
        do LP = LF,LL
          L = LSED(LP)
          if( TVAR3S(L) > 0.0 )then
            TVAR3S(L) = 1./TVAR3S(L)
          endif
        enddo
        do LP = LF,LL
          L = LSED(LP)
          ! *** PREVIOUS TIME STEP'S SEDIMENT BED SHEAR
          if( TAUBSEDS(L) > 0.0 )then
            TVAR3E(L) = QCELLCTR(L)/SQRT(TAUBSEDS(L))
          else
            TVAR3E(L) = 0.0
          endif
        enddo
                              
        ! ***  TVAR3E(L) = VELOCITY_MAGNITUDE/COHESIVE_BEDSTRESS                                                                  
        ! ***  TVAR3N(L) = D50/DEPTH                                                                                              
        do LP = LF,LL
          L = LSED(LP)
          TVAR3W(L) = TVAR3S(L)**0.333333
          TVAR3E(L) = TVAR3E(L)**0.333333
          TVAR3N(L) = TVAR3N(L)**0.333333
        enddo
        do LP = LF,LL
          L = LSED(LP)
          TVAR3S(L) = TVAR3S(L)**TMPEXP
        enddo
                              
        ! *** 
        do LP = LF,LL
          L = LSED(LP)
          if( CBEDTOTAL(L) > 0.0 )then
            K = KBT(L)
            HGDH(L) = 0.014*HDFUF(L)*QCELLAD1SQ(L)*(FRACCOH(L,K)*TVAR3W(L)*TVAR3E(L) &
                                                    +FRACNON(L,K)*TVAR3N(L))/CBEDTOTAL(L)
          endif
            
          HGDH(L) = HGDH(L)**0.75
          HGDH(L) = min(HGDH(L),1.0)
            
          ! *** CONVERT HGDH TO 1/HGDH  IE HDHG                                                                                     
          if( HGDH(L) > 0.) HGDH(L) = 1./HGDH(L)
        enddo

        do LP = LF,LL
          L = LSED(LP)
          if( TAUB(L) > 0.0 )then
            TAUBSED(L) = HGDH(L)**TMPEXP
            TAUBSND(L) = HGDH(L)**0.333333
            TVAR3W(L)  = QCELLCTR(L)*QCELLCTR(L)
          endif
        enddo

        do LP = LF,LL
          L = LSED(LP)
          if( TAUB(L) > 0.0 )then
            TAUBSED(L) = 0.026*TVAR3S(L)*TAUBSED(L)*TVAR3W(L)*QCELLAD1ZZ(L)
            TAUBSND(L) = 0.014*TVAR3N(L)*TAUBSND(L)*TVAR3W(L)*QCELLAD1SQ(L)
            ! *** save FOR NEXT ITERATION
            TAUBSEDS(L) = TAUBSED(L)
          endif
        enddo
                              
        !----------------------------------------------------------------------C                                                
        ! *** IF ISBEDSTR = 2, APPLY WEIGHTED AVERAGE TO BOTH SED AND SND                                                         
        if( ISBEDSTR == 2 .or. N == 1 )then
          do LP = LF,LL
            L = LSED(LP)
            K = KBT(L)
            TVAR3E(L) = FRACCOH(L,K)*TAUBSED(L)+FRACNON(L,K)*TAUBSND(L)
          enddo
          do LP = LF,LL
            L = LSED(LP)
            TAUBSED(L) = TVAR3E(L)
            TAUBSND(L) = TVAR3E(L)
          enddo
        endif
        
        ! *** ADJUST WAVE CELLS
        if( ISWAVE > 0 )then
          do LP = LF,LL
            L = LSED(LP)
            if( LWVMASK(L) )then
              ! *** Shear due to Current Only
              SHEAR = CTAUC(L)
              if( SHEAR > 0. )then
                ! *** COHESIVE FRACTION
                RAT1O = min(TAUBSED(L)/SHEAR,2.)
                RAT1O = max(RAT1O,0.2)
                TAUBSED(L) = TAUB(L)*RAT1O
                !IF(L == 16459)PRINT *,L,RAT1O,TAUB(L),TAUBSED(L)

                ! *** COHESIVE FRACTION
                RAT1O = min(TAUBSND(L)/SHEAR,2.)
                RAT1O = max(RAT1O,0.2)
                TAUBSND(L) = TAUB(L)*RAT1O
              endif
            endif
          enddo
        endif
      
        ! *** COMPUTE SHEAR VELOCITIES
        do LP = LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))
          USTARSND(L) = SQRT(TAUBSND(L))
        enddo

      endif  ! *** ENDIF ON GRAIN STRESS PARTITIONING FOR ISBEDSTR == 1 .or. ISBEDSTR == 2
                            
      !----------------------------------------------------------------------C                                                
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      ! *** INDEPENDENTLY SET GRAIN STRESS                                                                                    
      if( ISBEDSTR == 3 )then
        do LP = LF,LL
          L = LSED(LP)
          HDZBR = HP(L)/ZBRSED(L)
          TVAR3E(L) = 0.16/( (LOG(HDZBR)-1.)**2 )
        enddo
        do LP = LF,LL
          L = LSED(LP)
          TAUBSED(L) = TVAR3E(L)*QCELLCTR(L)*QCELLCTR(L)
        enddo
        do LP = LF,LL
          L = LSED(LP)
          TAUBSND(L) = TAUBSED(L)
        enddo
        do LP = LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))
          USTARSND(L) = SQRT(TAUBSND(L))
        enddo
                              
      endif
    
      !----------------------------------------------------------------------C                                                
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      !     USING WEIGHTED ROUGHNESS AND LOG RESISTANCE LAW                                                                   
      if( ISBEDSTR == 4 )then
        !     RECALCUATE ZOTOTAL AND CALCULATE ZOGRAIN                                                                          
        do LP = LF,LL
          L = LSED(LP)
          if( CBEDTOTAL(L) > 0.0 )then                                                                                          
            TMP = EXP(1.+0.4/SQRT(CBEDTOTAL(L)))                                                                                 
            ZOTOTAL(L) = HP(L)/TMP
          else
            ZOTOTAL(L) = ZBR(L)                                                                                                  
          endif                                                                                                                
          ZOGRAIN(L) = FRACNON(L,KBT(L))*SEDDIAGS(L,KBT(L))/30.
        enddo                                                                                                                  

        do LP = LF,LL
          L = LSED(LP)
          if( TAUBSED(L) > 0.0 )then                                                                                            
            ZOGRAINCOH = 0.041*VISMUDST/SQRT(TAUBSED(L))
            ZOGRAIN(L) = ZOGRAIN(L)+FRACCOH(L,KBT(L))*ZOGRAINCOH
            ZOGRAIN(L) = max(ZOGRAIN(L),ZOGRAINCOH)
          endif
        enddo                                                                                                                  

        !     ITERATE RVAL = SQRT(HG/H)                                                                                         
        do LP = LF,LL
          L = LSED(LP)
          TMPTOP = LOG(HP(L)/ZOTOTAL(L))-1.
          TMPBOT = LOG(HP(L)/ZOGRAIN(L))-1.
          RVAL = 1.                                                                                                              
          do KK = 1,100                                                                                                          
            RVAL = TMPTOP/(LOG(RVAL*RVAL)+TMPBOT)                                                                                
          enddo
          RVAL = min(RVAL,1.0)
          TAUBSED(L) = RVAL*RVAL*TAUB(L)
          TAUBSND(L) = RVAL*RVAL*TAUB(L)                                                                                         
        enddo                                                                                                                  

        do LP = LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))                                                                                         
          USTARSND(L) = SQRT(TAUBSND(L))                                                                                         
        enddo                                                                                                                  

      endif
                              
      !----------------------------------------------------------------------C                                                
      ! *** PARTITION BED STRESS BETWEEN TOTAL AND GRAIN STRESS                                                               
      !     USING WEIGHTED ROUGHNESS AND POWER RESISTANCE LAW                                                                 
      if( ISBEDSTR == 5 )then
                              
        ! *** RECALCUATE ZOTOTAL AND CALCULATE ZOGRAIN                                                                          
        do LP = LF,LL
          L = LSED(LP)
          if( CBEDTOTAL(L) > 0.0 )then                                                                                          
            TMP = CBEDTOTAL(L)/0.04736                                                                                           
            ZOTOTAL(L) = HP(L)*(TMP**3)
          else
            ZOTOTAL(L) = ZBR(L)                                                                                                  
          endif                                                                                                                
          ZOGRAIN(L) = FRACNON(L,KBT(L))*SEDDIAGS(L,KBT(L))/30.
        enddo                                                                                                                  
                              
        do LP = LF,LL
          L = LSED(LP)
          if( TAUBSED(L) > 0.0 )then                                                                                            
            ZOGRAINCOH = 0.041*VISMUDST/SQRT(TAUBSED(L))
            ZOGRAIN(L) = ZOGRAIN(L)+FRACCOH(L,KBT(L))*ZOGRAINCOH
            ZOGRAIN(L) = max(ZOGRAIN(L),ZOGRAINCOH)
          endif                                                                                                                
        enddo                                                                                                                  
                              
        !     CALCULATE GRAIN STRESS DIRECTLY                                                                                   
        do LP = LF,LL
          L = LSED(LP)
          TMP = (ZOGRAIN(L)/ZOTOTAL(L))**0.25
          TMP = min(TMP,1.0)
          TAUBSED(L) = TMP*TAUB(L)
          TAUBSND(L) = TMP*TAUB(L)                                                                                               
        enddo                                                                                                                  
                              
        do LP = LF,LL
          L = LSED(LP)
          USTARSED(L) = SQRT(TAUBSED(L))                                                                                         
          USTARSND(L) = SQRT(TAUBSND(L))                                                                                         
        enddo                                                                                                                  
                              
      endif
                              
      ! ***********************************************************************C                                                
      do LP = LF,LL
        L = LSED(LP)
        HBEDA(L) = 0.0
        do K = 1,KBT(L)
          HBEDA(L) = HBEDA(L)+HBED(L,K)
        enddo
      enddo
    endif  ! *** END OF SEDZLJ BYPASS
  
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO

  ! *** COMPUTE CELL CENTERED VELOCITY, WITH BOUNDARY AND CORNER CORRECTIONS
  if( (ICALC_BL > 0 .and. .not. LSEDZLJ) .or. (ICALC_BL > 0 .and. LSEDZLJ) )then
    ! *** COMPUTE CELL CENTERED TRANSPORT VECTORS FOR ONE ROW OF CELLS AROUND THE ACTIVE REGION
    if( ISBEDMAP > 0 )then
      do LP = 1,NSEDPLUS
        L = LSEDPLUS(LP)
        LN = LNC(L)
        LE = LEC(L)
        UCELLCTR(L) = 0.5*( RSSBCW(L)*WCORWST(L)*U(L,KSZ(L)) + RSSBCE(L)*WCOREST(L)*U(LE,KSZ(L)) )
        VCELLCTR(L) = 0.5*( RSSBCS(L)*WCORSTH(L)*V(L,KSZ(L)) + RSSBCN(L)*WCORNTH(L)*V(LN,KSZ(L)) )
        QCELLCTR(L) = SQRT(UCELLCTR(L)*UCELLCTR(L)+VCELLCTR(L)*VCELLCTR(L))
      
        ! *** CALCULATE UNIT VECTOR
        if( QCELLCTR(L) > 0.0 )then
          UCELLCTR(L)  = UCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, X COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          VCELLCTR(L)  = VCELLCTR(L)/QCELLCTR(L)             ! *** UNIT VECTOR, Y COMPONENT TO SPLIT BEDLOAD (DIMENSIONLESS)
          CBEDTOTAL(L) = TAUB(L)/(QCELLCTR(L)*QCELLCTR(L))
        else
          UCELLCTR(L)  = 0.
          VCELLCTR(L)  = 0.
          CBEDTOTAL(L) = 0.
        endif
      enddo
    endif
  endif
    
  ! ***********************************************************************C                                                
  ! *** CALCULATE COHESIVE SEDIMENT SETTLING, DEPOSITION AND RESUSPENSION                                                 
  if( ISTRAN(6) >= 1 )then
    if( LSEDZLJ )then
      call SEDZLJ_MAIN
      GOTO 1000  ! *** BYPASS ORIGINAL BED/WATER INTERFACE CALCULATION
    else
      call CALSED
    endif
  endif
                            
  ! ***********************************************************************C                                                
  ! *** CALCULATE NONCOHESIVE SEDIMENT BEDLOAD TRANSPORT, SETTLING,                                                       
  ! *** DEPOSITION AND RESUSPENSION                                                                                       
  if( ISTRAN(7) >= 1 ) CALL CALSND
                          
  ! ***********************************************************************C                                                
  ! *** CALCULATE BANK EROSION AND ADJUST SEDIMENT AND WATER VOLUME FLUXES                                                                                                            
  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
    if( ISBKERO >= 1 )then
      call BANKEROSED
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        QSBDTOP(LBANK) = QSBDTOP(LBANK)+QSBDTOPBEBKB(NP)
        QWBDTOP(LBANK) = QWBDTOP(LBANK)+QWBDTOPBEBKB(NP)
        QSBDTOP(LCHAN) = QSBDTOP(LCHAN)+QSBDTOPBECHB(NP)
        QWBDTOP(LCHAN) = QWBDTOP(LCHAN)+QWBDTOPBECHB(NP)
      enddo
    endif
  endif
                            
  !  $OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,LN,KTOPTP,KTOPM1) &
  !  $OMP             PRIVATE(DSEDGMM,DSEDGMMI,QSWPOS,QSWNEG,TMPVAL,TTHICK,SXD,SXD1)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = min(LF+LDMSED-1,LASED)

    ! ***********************************************************************C                                                
    ! *** CALCULATE PARENT TO ACTIVE LAYER SEDIMENT FLUX                                                                    
    do LP = LF,LL
      L = LSED(LP)
      QWATPA(L) = 0.0
      QSSDPA(L) = 0.0
    enddo
    
    if( ISNDAL == 2 )then
      if( IALTYP == 0 )then
        ! *** CONSTANT ACTIVE ARMOR LAYER THICKNESS
      
        ! *** Adjust cohesive sediment masses
        do NS = 1,NSED
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          do LP = LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            if( KTOPM1 > 0 )then
              QSWPOS = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPM1) )
              QSWNEG = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPTP) )
              SEDFPA(L,NS) = DSEDGMMI*( VFRBED(L,KTOPM1,NS)*MAX(QSWPOS,0.) + VFRBED(L,KTOPTP,NS)*MIN(QSWNEG,0.) )     ! *** FLUX OF SND BETWEEN ACTIVE/PARENT LAYERS (G/M2/S)
              QSSDPA(L) = QSSDPA(L) + DSEDGMM*SEDFPA(L,NS)
              QWATPA(L) = QWATPA(L) + DSEDGMM*( VDRBED(L,KTOPM1)*MAX(SEDFPA(L,NS),0.) + VDRBED(L,KTOPTP)*MIN(SEDFPA(L,NS),0.))
              SEDB(L,KTOPTP,NS) = SEDB(L,KTOPTP,NS) + DTSED*SEDFPA(L,NS)
              SEDB(L,KTOPM1,NS) = SEDB(L,KTOPM1,NS) - DTSED*SEDFPA(L,NS)
            endif
          enddo
        enddo

        ! *** Adjust non-cohesive sediment masses
        do NX = 1,NSND
          NS = NSED + NX
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          do LP = LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            if( KTOPM1 > 0 )then
              QSWPOS = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPM1) )
              QSWNEG = ( QSBDTOP(L) + QWBDTOP(L) )/( 1. + VDRBED1(L,KTOPTP) )
              SNDFPA(L,NX) = DSEDGMMI*( VFRBED(L,KTOPM1,NS)*MAX(QSWPOS,0.) + VFRBED(L,KTOPTP,NS)*MIN(QSWNEG,0.) )
              QSSDPA(L) = QSSDPA(L) + DSEDGMM*SNDFPA(L,NX)
              QWATPA(L) = QWATPA(L) + DSEDGMM*( VDRBED(L,KTOPM1)*MAX(SNDFPA(L,NX),0.) + VDRBED(L,KTOPTP)*MIN(SNDFPA(L,NX),0.) )
              SNDB(L,KTOPTP,NX) = SNDB(L,KTOPTP,NX) + DTSED*SNDFPA(L,NX)
              SNDB(L,KTOPM1,NX) = SNDB(L,KTOPM1,NX) - DTSED*SNDFPA(L,NX)
            endif
          enddo
        enddo
        
      else
        ! *** CONSTANT ACTIVE ARMOR LAYER TOTAL SEDIMENT MASS    
      
        ! *** Adjust cohesive sediment masses
        do NS = 1,NSED
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          do LP = LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            if( KTOPM1 > 0 )then
              SEDFPA(L,NS) = VFRBED(L,KTOPM1,NS)    *MAX(QSBDTOP(L),0.)   + VFRBED(L,KTOPTP,NS)*MIN(QSBDTOP(L),0.)
              QSSDPA(L) = QSSDPA(L)+SEDFPA(L,NS)                          
              QWATPA(L) = QWATPA(L)+VDRBED(L,KTOPM1)*MAX(SEDFPA(L,NS),0.) + VDRBED(L,KTOPTP)*MIN(SEDFPA(L,NS),0.)
              SEDFPA(L,NS) = DSEDGMMI*SEDFPA(L,NS)
              SEDB(L,KTOPTP,NS) = SEDB(L,KTOPTP,NS) + DTSED*SEDFPA(L,NS)
              SEDB(L,KTOPM1,NS) = SEDB(L,KTOPM1,NS) - DTSED*SEDFPA(L,NS)
            endif
          enddo
        enddo

        ! *** Adjust non-cohesive sediment masses
        do NX = 1,NSND
          NS = NSED + NX
          DSEDGMM  = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)
          DSEDGMMI = 1.E6*SSG(NS)         ! *** PURE SEDIMENT DENSITY (G/M**3)
          do LP = LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            if( KTOPM1 > 0 )then
              SNDFPA(L,NX) = VFRBED(L,KTOPM1,NS)    *MAX(QSBDTOP(L),0.)   + VFRBED(L,KTOPTP,NS)*MIN(QSBDTOP(L),0.)
              QSSDPA(L) = QSSDPA(L)+SNDFPA(L,NX)                          
              QWATPA(L) = QWATPA(L)+VDRBED(L,KTOPM1)*MAX(SNDFPA(L,NX),0.) + VDRBED(L,KTOPTP)*MIN(SNDFPA(L,NX),0.)
              SNDFPA(L,NX) = DSEDGMMI*SNDFPA(L,NX)
              SNDB(L,KTOPTP,NX) = SNDB(L,KTOPTP,NX) + DTSED*SNDFPA(L,NX)
              SNDB(L,KTOPM1,NX) = SNDB(L,KTOPM1,NX) - DTSED*SNDFPA(L,NX)
            endif
          enddo
        enddo
      endif
    endif
    
    ! *****************************************************************************                                           
    ! *** UPDATE TOP BED LAYER THICKNESS AND VOID RATIO                                                                     
    ! *** FOR DEPOSITION-RESUSPENSION STEP                                                                                  
    do LP = LF,LL
      L = LSED(LP)
      K = KBT(L)
      HBED1(L,K) = HBED(L,K)
      VDRBED1(L,K) = VDRBED(L,K)

      if( DEBUG )then
        ! *** PMC - FOR DEBUGGING
        SXD  = SUM(SEDB(L,K,:))  + SUM(SNDB(L,K,:))
        SXD1 = SUM(SEDB1(L,K,:)) + SUM(SNDB1(L,K,:))
        HBED(L,K) = HBED(L,K) - DTSED*( QSBDTOP(L) + QWBDTOP(L) ) + DTSED*( QSSDPA(L) + QWATPA(L) )
      
        if( SIGN(1.0,(HBED(L,K)-HBED1(L,K))) /= SIGN(1.0,(SXD-SXD1)) )then 
          TTHICK = 0.
          do NS = 1,NSED  
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SEDB(L,K,NS)*DSEDGMM
          enddo  
          do NX = 1,NSND  
            NS = NSED + NX  
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SNDB(L,K,NX)*DSEDGMM
          enddo 
          TTHICK = TTHICK/(1.0 - PORBED(L,K))
          if( ABS(TTHICK-HBED(L,K)) > 1E-3 )then
             PRINT '(F10.4,2I5,2F10.5,4E14.6,F10.5)', TIMEDAY ,L, K, HBED(L,K), HBED1(L,K), SXD, SXD1, DTSED*QSBDTOP(L), DTSED*QWBDTOP(L), TTHICK
          endif   
        endif
      endif

      if( K == 1 )then
        if( HBED(L,K) < MINTHICK )then
          ! *** ZERO NEGATIVE THICKNESSES
          HBED(L,K) = 0.0
          VDRBED(L,K) = SNDVDRD
          PORBED(L,K) = BEDPORC
          STDOCB(L,K) = 0.0
          STPOCB(L,K) = 0.0
          
          SEDB(L,K,:) = 0.0
          SNDB(L,K,:) = 0.0
        endif
      endif
      
      ! *** Update void ratio with the depositing sediment mass, using the last layer thickness.
      ! *** HBED then gets updated in CALBED
      TMPVAL = HBED1(L,K)/(1. + VDRBED1(L,K))
      TMPVAL = TMPVAL - DTSED*( QSBDTOP(L) - QSSDPA(L) )
      if( TMPVAL > 0.0 )then
        if( HBED(L,K) > 0.0 )then        
          VDRBED(L,K) = (HBED(L,K)/TMPVAL)-1.
          if( K == 1 )then
            ! *** LIMIT VOID RATIOS TO 0.01 >= N <= 99
            if( VDRBED(L,K) < 0.01 .or. VDRBED(L,K) > 99. )then
              VDRBED(L,K) = SNDVDRD  
            endif
          endif
        else
          ! *** Bed thickness is zero.  Initialize to depositing sediment mass
          HBED(L,K) = -DTSED*( QSBDTOP(L) - QSSDPA(L) )/(1. - PORBED(L,K))
          VDRBED(L,K) = PORBED(L,K)/(1. - PORBED(L,K))
        endif
      else
        VDRBED(L,K) = SNDVDRD          ! *** Assign default VR if layer empty or completely eroded away
      endif
    enddo
                            
    ! *** UPDATE PARENT LAYER BED LAYER THICKNESS AND VOID RATIO                                                            
    ! *** FOR DEPOSITION-RESUSPENSION STEP WHEN USING ACTIVE-PARENT SCHEME
    if( ISNDAL == 2 )then
      do LP = LF,LL
        L = LSED(LP)
        K = KBT(L)-1
        if( K > 0 )then
          HBED1(L,K) = HBED(L,K)
          VDRBED1(L,K) = VDRBED(L,K)
          HBED(L,K) = HBED(L,K) - DTSED*(QSSDPA(L)+QWATPA(L))

          if( K == 1 )then
            if( HBED(L,K) < 0.0 )then
              ! *** ZERO NEGATIVE THICKNESSES
              HBED(L,K) = 0.0
              VDRBED(L,K) = SNDVDRD
              PORBED(L,K) = BEDPORC
              STDOCB(L,K) = 0.0
              STPOCB(L,K) = 0.0
            endif
          endif
          TMPVAL = HBED1(L,K)/(1. + VDRBED1(L,K))
          TMPVAL = TMPVAL - DTSED*QSSDPA(L)
          if( TMPVAL > 0.0 )then
            VDRBED(L,K) = (HBED(L,K)/TMPVAL) - 1.
          else
            VDRBED(L,K) = SNDVDRD
          endif
        endif
      enddo
    endif
    
  enddo   ! *** END OF DOMAIN LOOP
  !  $OMP END PARALLEL DO
    
  ! **********************************************************************
  1000 continue  ! *** JUMP DIRECTLY HERE FROM SEDZLJ_MAIN CALL
                           
  ! ***********************************************************************C                                                
  ! *** CALCULATE TOXIC SETTLING, DEPOSITION AND RESUSPENSION                                                             
  if( ISTRAN(5) > 0 )then
    TTDSX = DSTIME(0)
    call CALTOX
    TSSTX = TSSTX + (DSTIME(0)-TTDSX)
  endif                      
                          
  if( .not. LSEDZLJ )then
    ! ***********************************************************************C                                                
    ! *** UPDATE TOP BED LAYER THICKNESS AND VOID RATIO                                                                     
    ! *** FOR DEPOSITION-RESUSPENSION STEP                                                                                  
    ! *** PRESENTLY ACTIVE BEFORE THE WATER COLUMN-BED TOXICS EXCHANGE                                                      
    ! *** CHECK PLACEMENT THERE AND HERE FOR                                                                                
    !        VDRTMP = (TMPVAL/HBED(L,K))-1.                                                                                   
    !        VDRBED(L,K) = VDRTMP                                                                                             

    ! ***********************************************************************C                                                
    ! *** UPDATE SEDIMENT BED LAYERING                                                                                      
    if( KB > 1 ) CALL CALBLAY
                          
    ! ***********************************************************************C                                                
    ! *** RESET SEDIMENT VOLUME FRACTIONS                                                                                   
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          BEDLINIT(L,K) = 0.
        enddo
      enddo
                            
      do NX = 1,NSED+NSND
        do K = 1,KB
          do LP = LF,LL
            L = LSED(LP)
            VFRBED(L,K,NX) = 0.
          enddo
        enddo
      enddo                                                                                                                  

      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)
            enddo
          enddo
        enddo
      endif
                            
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)
            enddo
          enddo
        enddo
      endif
                            
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            enddo
          enddo
        enddo
      endif
                            
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)
            enddo
          enddo
        enddo
      endif
                            
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              if( BEDLINIT(L,K) > 0.0 )then
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              else
                VFRBED(L,K,NS) = 0.0
              endif
            enddo
          enddo
        enddo
      endif
                            
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NSED + NX
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              if( BEDLINIT(L,K) > 0.0 )then
                VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K)
              else
                VFRBED(L,K,NS) = 0.0
              endif
            enddo
          enddo
        enddo
      endif
                            
      ! *** UPDATE VOLUME OF COHESIVES AND NON-COHESIVES
      do K = 1,KB                                                                                                              
        do LP = LF,LL
          L = LSED(LP)
          FRACCOH(L,K) = 0.0                                                                                                     
          FRACNON(L,K) = 0.0
        enddo
      enddo
                            
      do NS = 1,NSED
        do K = 1,KB                                                                                                              
          do LP = LF,LL
            L = LSED(LP)
            if( K <= KBT(L) )then                                                                                                  
              FRACCOH(L,K) = FRACCOH(L,K)+VFRBED(L,K,NS)                                                                           
            endif
          enddo
        enddo
      enddo
                            
      do NX = 1,NSND
        NS = NSED + NX                                                                                                             
        do K = 1,KB                                                                                                              
          do LP = LF,LL
            L = LSED(LP)
            if( K <= KBT(L) )then                                                                                                  
              FRACNON(L,K) = FRACNON(L,K)+VFRBED(L,K,NS)
            endif
          enddo
        enddo
      enddo
                            
      ! *** CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS                                                                    
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          if( K <= KBT(L) )then                                                                                                  
            VDRBEDSND(L,K) = SNDVDRD                                                                                             
            VDRBEDSED(L,K) = 0.0
            if( FRACCOH(L,K) > 0.0 )then
              VDRBEDSED(L,K) = ( (FRACCOH(L,K)+FRACNON(L,K))*VDRBED(L,K)-FRACNON(L,K)*SNDVDRD )/FRACCOH(L,K)
            endif
          else
            VDRBEDSND(L,K) = 0.0                                                                                                 
            VDRBEDSED(L,K) = 0.0                                                                                                 
          endif                                                                                                                
        enddo
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
  endif  ! *** END OF SEDZLJ BYPASS
  
  ! ***********************************************************************C                                                
  ! *** UPDATE SEDIMENT BED PHYSICAL PROPERTIES 
  ! *** FOR SEDZLJ, THE VARIABLES ARE UPDATED IN S_SEDZLJ.F90
  if( .not. LSEDZLJ )then
    call CALBED
  endif
  
  ! ***********************************************************************C                                                
  ! *** CHANGE BED MORPHOLOGY.  NOW INCLUDING SEDZLJ OPTION (2017-01)
  if( IMORPH > 0 )then
    ITMP = 0

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = min(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)
        P1(L)    = P(L)
        HBEDA(L) = 0.0
      enddo
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          HBEDA(L) = HBEDA(L) + HBED(L,K)
        enddo
      enddo 
      do LP = LF,LL
        L = LSED(LP)
        BELV(L) = ZELBEDA(L) + HBEDA(L)
        P(L) =  ( HP(L) + BELV(L) )*G
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
    
    ! *** ADJUST CONCENTRATIONS OF TRANSPORT VARIABLES IN RESPONSE TO CHANGE IN BED MORPHOLOGY
    ! *** 2017-01 - PMC REMOVED CONCENTRATION ADJUSTMENT SINCE NEW APPROACH MAINTAINS WATER BALANCE
  endif
  
  ! ***********************************************************************C                                                
  ! *** TOXICS CALCULATIONS
  if( ISTRAN(5) > 0 )then
    ! *** POREWATER ADVECTION AND DIFFUSION OF TOXICS                                                                       
    TTDSX = DSTIME(0)
      
    call CALTOXB
                          
    ! *** TOXIC CONTAMINANT REACTIONS
    if( NKINETICS > 0 ) CALL CALTOX_KINETICS
                      
    ! *** TOXIC CONTAMINANT OUTPUT TO FOOD CHAIN MODEL                                                                      
    if( ISFDCH >= 1 ) CALL FOODCHAIN(0)

    TSSTX = TSSTX + (DSTIME(0)-TTDSX)
  endif

  if( nactiveships > 0 )then
    prop_ero(:,:) = 0.0                       ! *** Zero erosion due to propwash for next iteration
    if( icalc_bl > 0 ) prop_bld(:,:) = 0.0    ! *** Zero erosion due to propwash for next iteration
  endif
  
  ! ***********************************************************************C                                                
  
  TTSED = TTSED + (DSTIME(0) - TTDS) - TTWAIT

  return

END

