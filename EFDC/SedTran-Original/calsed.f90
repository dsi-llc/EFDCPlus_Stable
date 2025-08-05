! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALSED  

  !**********************************************************************!
  ! ***  SUBROUTINE CALSED CALCULATES COHESIVE SEDIMENT SETTLING,  
  ! ***  DEPOSITION AND RESUSPENSION AND IS CALLED FROM SSEDTOX  
  !  
  ! *** STANDARD EFDC COHESIVE SEDIMENT TRANSPORT
  ! *** NOT USED FOR SEDZLJ
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2019-01       PAUL M. CRAIG     ADDED HARD BOTTOM BYPASS
  !    2015-07       PAUL M. CRAIG     CHANGED COHESIVE CONCENTRATIONS FOR SETTLING
  !                                      TO TOTAL CONCENTRATION (SEDT) RATHER THAN
  !                                      BY CLASS (SED)
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! 2011-03-02       PAUL M. CRAIG     RESTRUCTURED AND CORRECTED CODE
  !                                      REMOVED KC DEPENDENT DUPLCIATE CODE
  !                                      EXPANDED CAPABILITIES AND ADDED OMP
  !**********************************************************************!

  use GLOBAL
  use Variables_Propwash
  
  implicit none
  
  integer :: K, LP, L, NS, ND, LF, LL, LN, K1, IOBC, IFLAG
  real    :: TIME, STRESS, SHEAR, TAUBC, UTMP, VTMP, CURANG, TAUB2, UUSTARTMP, TAUDSS
  real    :: WVEL, CLEFT, CRIGHT, WESE, TAUE, TMPSTR, TMPSEDHID, TAUBHYDRO, PROBDEP, WSETMP
  real    :: SEDBTMP, CRNUM, GRADSED, SEDAVG, WESEMX, TAURTMP
  real,  external :: CSEDSET, FPROBDEP 

  integer,save,allocatable,dimension(:) :: IFIRST

  if( .not. allocated(IFIRST) )then
    allocate(IFIRST(2))
    IFIRST = 0

    if( ISEDVW == 0 )then
      ! *** CONSTANT (ONLY ASSIGN AT START OF RUN)
      do NS = 1,NSED2
        do K = 0,KS
          ! *** USING 2,LA INTENTIONAL TO ENSURE ASSIGNMENT EVEN IF DRY CELLS EXIST
          do L = 2,LA
            WSETA(L,K,NS) = WSEDO(NS)
          enddo
        enddo
      enddo
    endif
  endif
  
  !**********************************************************************!  
  if( ISDYNSTP == 0 )then  
    TIME = (DT*FLOAT(N)+TCON*TBEGIN)/TCON  
  else  
    TIME = TIMESEC/TCON  
  endif  

  !**********************************************************************!
  ! **
  ! ***  COHESIVE SEDIMENT FLUX

  if( LADRY > 0 )then
    do K = 1,KC
      do LP = 1,LADRY
        L = LDRY(LP)  
        TVAR2E(L,K) = 0.
        TVAR2N(L,K) = 0.
      enddo
    enddo
  endif
  
  do NS = 1,NSED2
    DSEDGMM = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,K1,LN,IOBC)  &
    !$OMP             PRIVATE(STRESS,SHEAR,TAUBC,UTMP,VTMP,CURANG,TAUB2,UUSTARTMP,TAUDSS) &
    !$OMP             PRIVATE(WVEL,CLEFT,CRIGHT,WESE,TAUE,TMPSTR,TMPSEDHID,TAUBHYDRO,PROBDEP,WSETMP) &
    !$OMP             PRIVATE(SEDBTMP,CRNUM,GRADSED,SEDAVG,WESEMX,TAURTMP)
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      
      !----------------------------------------------------------------------!
      ! ***  SET SETTLING VELOCITIES  ( ISEDVW == 0 SET ABOVE )
      if( ISEDVW == 0 .and. ISTOPT(6) > 0 .and. KC > 1 )then
        ! *** Reset settling velocities back to user specified if using sediment anti-diffusion (ISTOPT(6) == 1)
        do K = 0,KS
          ! *** USING 2,LA INTENTIONAL TO ENSURE ASSIGNMENT EVEN IF DRY CELLS EXIST
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSEDO(NS)
          enddo
        enddo
      endif

      if( ISEDVW == 1 )then
        ! *** SIMPLE CONCENTRATION DEPENDENT (HUANG AND METHA)
        do K = 0,KS
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),0.0,ISEDVW)
          enddo
        enddo
      endif

      if( ISEDVW == 2 )then
        ! *** CONCENTRATION AND SHEAR/TURBULENCE DEPENDENT COHESIVE SETTLING  (SHRESTA AND ORLOB)
        if( ISWAVE == 0 )then
          do K = 0,KS
            do LP = 1,LLWET(K+1,ND)
              L = LKWET(LP,K+1,ND)  
              STRESS = QQ(L,0)/CTURB2
              SHEAR = 2.*HPKI(L,K+1)*SQRT(STRESS)/VKC
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),SHEAR,ISEDVW)
            enddo
          enddo
        else
          K = 0
          do LP = LF,LL
            L = LWET(LP)  
            K1 = KSZ(L)
            if( LWVMASK(L) )then
              TAUBC = QQ(L,0)/CTURB2
              UTMP = 0.5*STCUV(L)*(U(LEC(L),   K1) + U(L,K1))+1.E-12
              VTMP = 0.5*STCUV(L)*(V(LNC(L),K1) + V(L,K1))
              CURANG = ATAN2(VTMP,UTMP)
              TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
              TAUB2 = MAX(TAUB2,0.)
              STRESS = SQRT(TAUB2)
              SHEAR = 2.*HPKI(L,K1)*SQRT(STRESS)/VKC
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),SHEAR,ISEDVW)
            else
              STRESS = QQ(L,0)/CTURB2
              SHEAR = 2.*HPKI(L,K1)*SQRT(STRESS)/VKC
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),SHEAR,ISEDVW)
            endif
          enddo
          if( KC > 1 )then
            do K = 1,KS
              do LP = 1,LLWET(K+1,ND)
                L = LKWET(LP,K+1,ND)  
                if( LWVMASK(L) )then
                  SHEAR = HPI(L)*SQRT( DZIGSD4U(L,K) )*SQRT( (U(LEC(L),K+1)-U(LEC(L),K)+U(L,K+1)-U(L,K))**2  &
                                                           + (V(LNC(L),K+1)-V(LNC(L),K)+V(L,K+1)-V(L,K))**2 )
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),SHEAR,ISEDVW)
                else
                  STRESS = QQ(L,0)/CTURB2
                  SHEAR = 2.*HPKI(L,K+1)*SQRT(STRESS)/VKC
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),SHEAR,ISEDVW)
                endif
              enddo
            enddo
          endif
        endif  
      endif

      if( ISEDVW >= 3 .and. ISEDVW <= 4 )then
        ! *** CONCENTRATION AND SHEAR/TURBULENCE DEPENDENT COHESIVE SETTLING  (ZIEGLER AND NESBIT)
        if( ISWAVE == 0 )then
          do K = 0,KS
            do LP = 1,LLWET(K+1,ND)
              L = LKWET(LP,K+1,ND)  
              STRESS = 0.5*QQ(L,0)/CTURB2
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
            enddo
          enddo
        else
          K = 0
          do LP = LF,LL
            L = LWET(LP)  
            K1 = KSZ(L)
            if( LWVMASK(L) )then
              TAUBC = QQ(L,0)/CTURB2
              UTMP = 0.5*STCUV(L)*(U(LEC(L),   K1) + U(L,K1))+1.E-12
              VTMP = 0.5*STCUV(L)*(V(LNC(L),K1) + V(L,K1))
              CURANG = ATAN2(VTMP,UTMP)
              TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
              TAUB2 = MAX(TAUB2,0.)
              STRESS = SQRT(TAUB2)
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),STRESS,ISEDVW)
            else
              STRESS = 0.5*QQ(L,0)/CTURB2
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),STRESS,ISEDVW)
            endif
          enddo

          if( KC > 1 )then
            do K = 1,KS
              do LP = 1,LLWET(K+1,ND)
                L = LKWET(LP,K+1,ND)  
                if( LWVMASK(L) )then
                  LN = LNC(L)
                  STRESS = AV(L,K)*SQRT( DZIGSD4U(L,K)*(U(LEC(L),K+1)-U(LEC(L),K) + U(L,K+1)-U(L,K))**2  &
                                    +  DZIGSD4V(L,K)*(V(LN ,K+1)-V(LN ,K) + V(L,K+1)-V(L,K))**2 )
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
                else
                  STRESS = 0.5*QQ(L,0)/CTURB2
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
                endif
              enddo
            enddo
          endif
        endif
      endif

      if( ISEDVW == 5 )then
        ! *** CONCENTRATION AND SHEAR/TURBULENCE DEPENDENT COHESIVE SETTLING  (HOUSATONIC)
        do K = 0,KS
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            UUSTARTMP = QCELLCTR(L)*SQRT(QQ(L,0)/CTURB2)
            STRESS = SQRT(HPI(L)*HPI(L)*UUSTARTMP)
            WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
          enddo
        enddo
      endif

      ! *** Treat "fast" settling classes for propwash resuspended sediments
      if( NS > NSED )then
        do K = 0,KS
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSETA(L,K,NS)*fast_multiplier
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
        do IOBC = 1,NBCSOP  
          L = LOBCS(IOBC)  
          WSETA(L,K,NS) = 0.0
        enddo  
      enddo  

      !----------------------------------------------------------------------!
      !
      ! ***  HORIZONTAL LOOPS

      ! *** WATER COLUMN
      if( KC > 1 )then
        ! *** SET FLUX FOR THE BOTTOM OF THE TOP LAYER
        K = KC
        do LP = LF,LL
          L = LWET(LP)  
          SEDF(L,K,NS) = 0.
          WVEL = DTSED*HPKI(L,K)
          CLEFT = 1. + WSETA(L,K-1,NS)*WVEL
          CRIGHT = MAX(SED(L,K,NS),0.)
          SED(L,K,NS) = CRIGHT/CLEFT
          SEDF(L,K-1,NS) = -WSETA(L,K-1,NS)*SED(L,K,NS)
        enddo

        ! *** SET FLUX FOR MIDDLE LAYERS
        do K = KS,2,-1
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND)
            WVEL = DTSED*HPKI(L,K)
            CLEFT = 1. + WSETA(L,K-1,NS)*WVEL
            CRIGHT = MAX(SED(L,K,NS),0.) - SEDF(L,K,NS)*WVEL
            SED(L,K,NS) = CRIGHT/CLEFT
            SEDF(L,K-1,NS) = -WSETA(L,K-1,NS)*SED(L,K,NS)
          enddo
        enddo
        
        ! *** HANDLE HARD BOTTOM BOTTOM LAYER SETTLING
        if( ISBEDMAP > 0 )then
          do LP = LF,LL
            L = LWET(LP)  
            if( LBED(L) )then
              ! *** ADD SETTLING FROM THE LAYER ABOVE
              K = KSZ(L)
              WVEL = DTSED*HPKI(L,K)  
              SED(L,K,NS) = MAX(SED(L,K,NS),0.) - SEDF(L,K,NS)*WVEL  
            endif
          enddo
        endif
      endif
      
      ! *** Collapse water column sediments to the for propwash resuspended sediments
      if( NS > NSED )then
        do K = 0,KS
          do LP = 1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSETA(L,K,NS)*fast_multiplier
          enddo
        enddo
      endif
      
      ! *** SET CONSTANTS
      TMPSEDHID = 1.0
      if( IWRSP(1) == 4 )then
        TMPSTR = 0.0
      else
        TMPSTR = 1.0
      endif

      ! *** UPDATE SEDIMENT BED MASS & BOTTOM LAYER CONCENTRATION
      do LP = LF,LL
        L = LWET(LP)
        if( LBED(L) ) CYCLE     
        
        ! *** SET MAXIMUM EROSION RATE
        WESEMX = DELTI*SEDB(L,KBT(L),NS)

        if( TAUBSED(L) > TAURB(L,KBT(L)) )then
          ! *** MASS/BULK EROSION
          WESE = WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
          WESE = MIN(WESE,WESEMX)
        else
          TAUE = 0.
          if( TAUBSED(L) > TAURS(L,KBT(L)) )then
            ! *** SURFACE EROSION
            WESE = WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)

            ! *** SET NORMALIZING SHEAR BASED ON INPUT OPTION
            if( IWRSP(1) >= 99 )then
              ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU FROM SSCOHSEDPMAP
              TAURTMP = TAUNS(L,KBT(L))
            elseif( (IWRSP(1) >= 2) .and. (IWRSP(1) < 4) )then
              ! *** 1 HWANG AND METHA - LAKE OKEECHOBEE
              ! *** 2 HAMRICK'S MODIFICATION OF SANFORD AND MAA
              ! *** 3 SAME AS 2 EXCEPT VOID RATIO OF COHESIVE SEDIMENT FRACTION IS USED
              TAURTMP = TAUR(1)
            else
              ! *** Use DIRECTLY INPUT PROPERTIES
              TAURTMP = TAURS(L,KBT(L))
            endif

            TAUE = (TAUBSED(L)-TMPSTR*TAURS(L,KBT(L)))/TAURTMP
            TAUE = MAX(TAUE,0.0)

            ! *** SET NON-COHESIVE HIDING FACTOR
            if( ISTRAN(7) >= 1 .and. COSEDHID(1) /= 0.0 )then
              TMPSEDHID = (FRACCOH(L,KBT(L)))**COSEDHID(1)
            endif
              
            ! *** FINALIZE SURFACE EROSION RATE
            if( IWRSP(1) < 99 )then
              WESE = TMPSEDHID*WESE*( TAUE**TEXP(NS) )
            else
              ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU EXPONENT FROM SSCOHSEDPMAP
              WESE = TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
            endif
            WESE = MIN(WESE,WESEMX)
          else
            ! *** NO EROSION 
            WESE = 0.0
          endif
        endif

        ! *** SET PROBABILITY OF DEPOSITION 
        if( IWRSP(1) < 99 )then
          TAUDSS = TAUD(NS)
        else
          TAUDSS = TAUDS(L)
        endif
        TAUBHYDRO = QQ(L,0)/CTURB2
        PROBDEP = 0.
        if( ISPROBDEP(NS) == 0 )then
          ! *** KRONE PROBABILITY OF DEPOSITION USING COHESIVE GRAIN STRESS
          if( TAUBSED(L) < TAUDSS) PROBDEP = (TAUDSS-TAUBSED(L))/TAUDSS
        elseif( ISPROBDEP(NS) == 1 )then
          ! *** KRONE PROBABILITY OF DEPOSITION USING TOTAL BED STRESS
          if( TAUBHYDRO < TAUDSS) PROBDEP = (TAUDSS-TAUBHYDRO)/TAUDSS
        elseif( ISPROBDEP(NS) == 2 )then
          ! *** PARTHENIADES PROBABILITY OF DEPOSITION USING COHESIVE GRAIN STRESS
          if( TAUBSED(L) <= TAUDSS )then
            PROBDEP = 1.0
          else
            PROBDEP = FPROBDEP(TAUDSS,TAUBSED(L))
          endif
        elseif( ISPROBDEP(NS) == 3 )then
          ! *** PARTHENIADES PROBABILITY OF DEPOSITION USING TOTAL BED STRESS
          if( TAUBHYDRO <= TAUDSS )then
            PROBDEP = 1.0
          else
            PROBDEP = FPROBDEP(TAUDSS,TAUBHYDRO)
          endif
        endif
        if( SED(L,KSZ(L),NS) > SEDMDGM ) PROBDEP = 1.
        
        ! *** Handle Propwash
        if( PROP_ERO(L,0) > 0.0 )then
          PROP_ERO(L,NS) = PROP_ERO(L,NS)*DXYIP(L)*DELTI                        ! *** Convert mass from g to g/m**2/s
          WSETMP = 0.0                                                          ! *** Disable settling for active propwash cell
          !IF( ISPROPWASH == 2 )then
          !  WESE = PROP_ERO(L,NS)                                              ! *** Only allow propwash erosion to avoid double counting
          !ELSE
            WESE = WESE + PROP_ERO(L,NS)                                        ! *** Allow ambient current erosion for all propwash options (2021-11-10)
          !ENDIF
        else
          WSETMP = PROBDEP*WSETA(L,0,NS)
        endif
        
        WVEL   = DTSED*HPKI(L,KSZ(L))
        CLEFT  = 1. + WSETMP*WVEL
        CRIGHT = MAX(SED(L,KSZ(L),NS),0.) + (WESE-SEDF(L,KSZ(L),NS))*WVEL
        SED(L,KSZ(L),NS) = CRIGHT/CLEFT
        SEDF(L,0,NS)     = -WSETMP*SED(L,KSZ(L),NS) + WESE
        SEDBTMP          = SEDB(L,KBT(L),NS) - DTSED*SEDF(L,0,NS)

        ! *** Limit Erosion to Available Mass
        if( SEDBTMP < 0.0 )then
          SEDF(L,0,NS) = DELTI*SEDB(L,KBT(L),NS)
          SEDBTMP = 0.0
          SED(L,KSZ(L),NS) = SED2(L,KSZ(L),NS) + (SEDF(L,0,NS)-SEDF(L,KSZ(L),NS))*WVEL
        endif
        SEDB1(L,KBT(L),NS) = SEDB(L,KBT(L),NS)
        SEDB(L,KBT(L),NS) = SEDBTMP
        QSBDTOP(L) = QSBDTOP(L) + DSEDGMM*SEDF(L,0,NS)
        
        if( IBMECH == 0 .or. SEDVRDT < 0.0 )then
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*VDRBED(L,KBT(L))*SEDF(L,0,NS)  ! *** IF EITHER CONSTANT OR INSTANTLY CONSOLIDATING, USE BED VR  
        else
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.) + VDRDEPO(NS)*MIN(SEDF(L,0,NS),0.) )
        endif
      enddo      ! *** END OF UPDATE SEDIMENT BED MASS & BOTTOM LAYER LOOP

      !----------------------------------------------------------------------!
      ! ***  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC > 1
      if( ISTOPT(6) == 1 .and. KC > 1 )then

        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            CRNUM = 1.+DTSED*WSETA(L,K,NS)*HPKI(L,K+1)
            GRADSED = (SED(L,K+1,NS)-SED(L,K,NS))/(DZC(L,K+1)+DZC(L,K))
            SEDAVG = 0.5*(SED(L,K+1,NS)+SED(L,K,NS)+1.E-16)
            WSETA(L,K,NS) = MIN(-CRNUM*DZC(L,K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG, 1.0)
          enddo
        enddo

        ! *** TVAR1S = LOWER DIAGONAL
        do LP = LF,LL
          L = LWET(LP)  
          TVAR1S(L,KSZ(L)) = 0.
        enddo
        do K = 2,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TVAR1S(L,K) = MIN(WSETA(L,K-1,NS),0.)
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
          TVAR1W(L,KSZ(L)) = DELTI*HPK(L,KSZ(L)) - MIN(WSETA(L,KSZ(L),NS),0.)
          TVAR1W(L,KC)     = DELTI*HPK(L,KC)     + MAX(WSETA(L,KC-1,NS),0.)
        enddo
        do K = 2,KS
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            TVAR1W(L,K) = DELTI*HPK(L,K)+MAX(WSETA(L,K-1,NS),0.)-MIN(WSETA(L,K,NS),0.)
          enddo
        enddo

        ! *** TVAR1E = RIGHT HAND SIDE
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TVAR1E(L,K) = DELTI*HPK(L,K)*SED(L,K,NS)
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
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TVAR2N(L,K) = TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
          enddo
        enddo
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            SED(L,K,NS) = TVAR2N(L,K)
          enddo
        enddo
      endif   ! *** END ANTI-DIFFUSION

      !----------------------------------------------------------------------!
      !
      ! ***  FINAL FLUX
      if( KC > 1 .and. ( ISTRAN(5) > 0 )  )then
        ! *** KC-1 LAYER
        do LP = 1,LLWET(KS,ND)
          L = LKWET(LP,KS,ND)  
          SEDF(L,KS,NS) = DELTI*DZC(L,KC)*HP(L)*(SED(L,KC,NS)-SED2(L,KC,NS))
        enddo 

        ! *** MIDDLE LAYERS
        do K = KS-1,1,-1
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            SEDF(L,K,NS) = DELTI*DZC(L,K+1)*HP(L)*(SED(L,K+1,NS)-SED2(L,K+1,NS))+SEDF(L,K+1,NS)
          enddo  
        enddo  
      endif
      
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
    
  enddo     ! *** END OF SEDIMENT CLASS LOOP

  !**********************************************************************!
  if( ISDTXBUG > 0 )then
    IFLAG = 0  
    do NS = 1,NSED2  
      do L = 2,LA
        K = KSZ(L)
        if( SED(L,K,NS) < 0. )then
          if( IFLAG == 0 )then  
            open(1,FILE = OUTDIR//'NEGSEDSND.OUT',POSITION = 'APPEND')  
            IFLAG = 1  
          endif  
          write(1,107) TIME,NS,IL(L),JL(L),K,SED(L,K,NS)  
          SED(L,K,NS) = 0.0    ! *** Continue with warning
        endif  
      enddo  
    enddo  

    do NS = 1,NSED  
      do L = 2,LA
        K = KBT(L)
        if( SEDB(L,K,NS) < 0. .or. HBED(L,K) < 0. )then
          if( IFLAG == 0 )then  
            open(1,FILE = OUTDIR//'NEGSEDSND.OUT',POSITION = 'APPEND')  
            IFLAG = 1  
          endif  
          write(1,108) TIME,NS,IL(L),JL(L),K,SEDB(L,K,NS),SEDF(L,0,NS),HBED(L,K)
          if( SEDB(L,K,NS) < 0. ) SEDB(L,K,NS) = 0.0    ! *** Continue with warning
          if( HBED(L,K) < 0. )    HBED(L,K) = 0.0       ! *** Continue with warning
        endif  
      enddo  
    enddo  

    if( IFLAG == 1 ) close(1) 
  endif
  
  107 FORMAT(' Warning: WC  SED < 0: TIME, NS, I, J, K, NEGSED = ',F12.4,4I5,4E13.4)       
  108 FORMAT(' Warning: BED SED < 0: TIME, NS, I, J, K, NEGSEDB, SEDF, HBED = ',F12.4,4I5,4E13.4)       

  !**********************************************************************!

  return

END

