! ----------------------------------------------------------------------
!   This file is a part of EFDC + 
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALTOX

  ! ***  SUBROUTINE CALTOX CALCULATES CONTAMINANT TRANSPORT.
  ! ***  IT IS CALLED FROM SSEDTOX
  ! ***  USED BY BOTH THE STANDARD EFDC SEDTRAN MODULE AND SEDZLJ
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2014-08-26        PAUL M. CRAIG    Updated to allow for any mix of ISTOC
  ! 2012-12-05        PAUL M. CRAIG    Updated OMP
  ! 2012-10-02        Dang H Chung     Added OMP
  !**********************************************************************!

  use GLOBAL
  use Variables_MPI
  
  implicit none

  integer :: ND, LF, LL, LP, L, IVAL, NT, NS, NX, K, NFD, LUTMP, LDTMP, KTOPM1, NP
  integer :: LBANK, LCHAN, KBANK, LE, LW, LS, LN, NSB, KCHAN, KTOPTP
  integer,  save :: NDUMP
  
  real :: TMPEXP, TMPVAL, TMPTOXB, TMPTOXC, TMPTOXE, TMPTOXW, TMPTOXN, TMPTOXS, TOXFLUX, CBLTOXTMP
  real :: AA11, BB11, BB22, SNDFEFF, CLEFT, CRIGHT, WVEL
  real :: FTPOS, FTNEG, HBEDMIN0, SORBMIN
  real,save,allocatable,dimension(:) :: TOXFPA
  
  if( .not. allocated(TOXFPA) )then
    allocate(TOXFPA(LCM))
    TOXFPA = 0.0
    NDUMP = 0
    
    ! *** VALIDATE INITIAL CONDITIONS
    do L  = 2,LA
      do K = 1,KB
        if( HBED(L,K) <= 1E-9 )then
          TOXB(L,K,:) = 0.0
          TOXB1(L,K,:) = 0.0
        endif
      enddo
    enddo
    
  endif

  if( LSEDZLJ )then
      HBEDMIN0 = 1E-9
  else
      HBEDMIN0 = 1E-4
  endif

  !**********************************************************************CC
  ! ***  UPDATE FRACTION OF PARTICULATE ORGANIC CARBON IN BED
  if( ISTPOCB == 4 )then
    IVAL = 0
    do NT = 1,NTOX
      if( ISTOC(NT) >= 2) IVAL = 1
    enddo
  
    if( IVAL == 1 )then
      ! *** HOUSATONIC
      call SETFPOCB(1)
    endif
  endif
  
  !**********************************************************************C
  ! ***  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN WATER COLUMN
  ! **
  ! ***  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
  ! ***  TOXPARW(L,NS,NT) = PARTITION COEFFICIENT IN WATER COLUMN
  ! ***  TOXPFTW(L,K,NT)  = TOTAL PARTICULATE FRACTION IN WATER COLUMN
  !                       USED AS TEMPORARY VARIBLE IN THIS AND
  !                       FOLLOWING CODE BLOCK
  NFD = NSED2 + NSND + 1

  NP = PRECISION(SORBMIN)
  SORBMIN = 0.99  ! 1. - 1./10**NP
  
  ! *** OMP LOOP OVER THE TOXICS COMPUTATIONS (SETTLING, DEPOSTION AND EROSION ONLY)
  !$OMP PARALLEL DEFAULT(SHARED)
 
  !$OMP DO PRIVATE(ND, LP, L, K, NT, NS, NX, TMPEXP, TMPVAL)
  do ND = 1,NDM
    !**********************************************************************C
    ! ***  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN WATER COLUMN
    ! **
    ! ***  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
    ! ***  TOXPARW(L,NS,NT)  = PARTITION COEFFICIENT IN WATER COLUMN
    ! ***  TOXPFTW(L,K,NT)   = TOTAL PARTICULATE FRACTION IN WATER COLUMN

    do NT = 1,NTOX
      do NS = 1,NSP2(NT)
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOXPFW(L,K,NS,NT) = 0.
          enddo
        enddo
      enddo
    enddo

    do NT = 1,NTOX
      ! *** PARTITION TO COHESIVE
      if( ISTRAN(6) >= 1 )then
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do NS = 1,NSED2
            if( ITXPARW(NS,NT) == 0 )then   ! *** NON-SOLIDS BASED PARTITIONING
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  ! ***            =      G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SED(L,K,NS)*STFPOCW(L,K,NS)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif
            if( ITXPARW(NS,NT) == 1 )then   ! *** SOLIDS BASED PARTITIONING
              TMPEXP = CONPARW(NS,NT)
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  TMPVAL = SED(L,K,NS)**TMPEXP
                  ! ***             =           G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SED(L,K,NS)*STFPOCW(L,K,NS)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif
          enddo
          
          ! *** Propwash fast settling
          if( NSED2 > NSED )then
          endif

        elseif( ISTOC(NT) == 0 )then
          ! *** Kd APPROACH
          do NS = 1,NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  ! ***             =     G/M3      L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SED(L,K,NS)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif
            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  TMPVAL = SED(L,K,NS)**TMPEXP
                  ! ***             =          G/M3       L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SED(L,K,NS)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif
          enddo

          ! *** Propwash fast settling
          if( NSED2 > NSED )then
          endif
        endif
      endif

      ! *** PARTITION TO NONCOHESIVE
      if( ISTRAN(7) >= 1 )then
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do NX = 1,NSND
            NS = NX + NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  ! ***             =     G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SND(L,K,NX)*STFPOCW(L,K,NS)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif

            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  TMPVAL = SND(L,K,NX)**TMPEXP
                  ! ***             =            G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SND(L,K,NX) *STFPOCW(L,K,NS)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif
          enddo
        elseif( ISTOC(NT) == 0 )then
          ! *** Kd APPROACH
          do NX = 1,NSND
            NS = NX + NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  ! ***             =     G/M3     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SND(L,K,NX)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif

            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do K = 1,KC
                do LP = 1,LLWET(K,ND)
                  L = LKWET(LP,K,ND)  
                  TMPVAL = SND(L,K,NX)**TMPEXP
                  ! ***             =            G/M3      L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SND(L,K,NX)*TOXPARW(L,NS,NT)
                enddo
              enddo
            endif
          enddo
        endif
      endif

      ! *** PARTITION (COMPLEX TO DISSOLVED ORGANIC CARBON)
      if( ISTOC(NT) == 1 .or. ISTOC(NT) == 2 )then
        if( ITXPARWC(1,NT) == 0 )then
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***                 DOC(G/M3)    L/MG (M3/G)
              TOXPFW(L,K,NFD,NT) = STDOCW(L,K)*TOXPARWC(1,NT)
            enddo
          enddo
        endif

        if( ITXPARWC(1,NT) == 1 )then
          TMPEXP = CONPARWC(1,NT)
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TMPVAL = 1.
              if( STDOCW(L,K) > 0.) TMPVAL = STDOCW(L,K)**TMPEXP
              ! ***                       DOC(G/M3)    L/MG (M3/G)
              TOXPFW(L,K,NFD,NT) = TMPVAL*STDOCW(L,K)*TOXPARWC(1,NT)
            enddo
          enddo
        endif
      endif

      ! *** PARTITION TO PARTICULATE ORGANIC CARBON
      if( ISTOC(NT) == 1 )then
        if( ITXPARWC(2,NT) == 0 )then
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              ! ***                  POC(G/M3)    L/MG (M3/G)
              TOXPFW(L,K,NFD+1,NT) = STPOCW(L,K)*TOXPARWC(2,NT)
            enddo
          enddo
        endif

        if( ITXPARW(NS,NT) == 1 )then
          TMPEXP = CONPARW(NS,NT)
          do K = 1,KC
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TMPVAL = STPOCW(L,K)**TMPEXP
              ! ***                         POC(G/M3)    L/MG (M3/G)
              TOXPFW(L,K,NFD+1,NT) = TMPVAL*STPOCW(L,K)*TOXPARWC(2,NT)
            enddo
          enddo
        endif
      endif
    enddo   ! *** NTOX

    ! *** TOXPFTW IS TEMPORARILY USED TO STORE TOTAL SORBED (DIMENSIONLESS)
    do NT = 1,NTOX
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          TOXPFTW(L,K,NT) = 0.
        enddo
      enddo
      do NS = 1,NSP2(NT)
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOXPFTW(L,K,NT) = TOXPFTW(L,K,NT) + TOXPFW(L,K,NS,NT)
          enddo
        enddo
      enddo
    enddo   ! *** NTOX
  enddo
  !$OMP END DO   
  
  !$OMP DO PRIVATE(ND, K, LP, L, NT, NS)
  do ND = 1,NDM  
    ! *** COMPUTE THE SORBED FRACTION (TOXPFW) (DIMENSIONLESS)
    do NT = 1,NTOX
      do NS = 1,NSP2(NT)
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOXPFW(L,K,NS,NT) = TOXPFW(L,K,NS,NT)/(1. + TOXPFTW(L,K,NT))
          enddo
        enddo
      enddo
    enddo   ! *** NTOX

    if( ISTMSR >= 1 .or. ISFDCH > 0 )then
      ! *** ONLY NEEDED FOR TIME SERIES USING TMSR OR FOODCHAIN
      do NT = 1,NTOX
        do K = 1,KC
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOXFDFW(L,K,NT) = 1./(1. + TOXPFTW(L,K,NT))
            TOXCDFW(L,K,NT) = TOXPFW(L,K,NFD,NT)
          enddo
        enddo
      enddo
    endif
  enddo
  !$OMP END DO
  
  !**********************************************************************C
  !
  ! ***  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
  !
  ! ***  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED 
  ! ***  TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED (INTERIM UNITS: METERS)
  !$OMP DO PRIVATE(ND, LF, LL)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)
    call CALTOXB_FRACTIONS(LF,LL)
  enddo
  !$OMP END DO 
  
  ! ******************************************************************************
  ! ******************************************************************************
  !
  ! *** CALCULATE PARTICULATE TOXIC CONTAMINANT SETTLING AND BED EXCHANGE FLUXES

  ! *** TOXF(L,1:KS,NT) = TOXIC CONTAMINANT SETTLING AND BED EXCHANGE FLUX 
  ! ***                    ( + ) UPWARD FLUX,  (-) DOWNWARD FLUX   (M/S)
  ! *** TOXF(L,0,NT)    = TOXIC CONTAMINANT DEPOSITIONAL FLUX      (M/S)
  ! *** TOXFB(L,NT)     = TOXIC CONTAMINANT EROSIONAL FLUX         (FINAL 1/S, INTERIM M/S)
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT, NS, NX, TMPEXP, TMPVAL)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)
    
    do NT = 1,NTOX
      ! *** Water column
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          TOXF(L,K,NT) = 0.
        enddo
      enddo
      ! *** Bed/water column interface
      do LP = LF,LL
        L = LSED(LP)
        TOXF(L,0,NT) = 0.
      enddo
    enddo

    ! ******************************************************************************
    ! *** TOXICS VERTICAL FLUX (TOXF) IN THE WATER COLUMN (K = 1,KS)
    ! *** SEDF (G/M2/S) IS ALWAYS NEGATIVE IN THE WATER COLUMN
    if( KC >= 2 )then
      do NT = 1,NTOX
        ! *** PARTICLE COHESIVE FLUX
        if( ISTRAN(6) >= 1 )then
          if( ISTOC(NT) > 1 )then
            ! *** fPOC BASED
            do NS = 1,NSED2
              if( ITXPARW(NS,NT) == 0 )then
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    ! *** M/S          M/S           G/M2/S                         M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SEDF(L,K,NS)*STFPOCW(L,K+1,NS)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              elseif( ITXPARW(NS,NT) == 1 )then
                TMPEXP = CONPARW(NS,NT)
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    TMPVAL = 1.
                    if( SED(L,K+1,NS) > 0.) TMPVAL = SED(L,K+1,NS)**TMPEXP
                    ! *** M/S          M/S                  G/M2/S                         M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SEDF(L,K,NS)*STFPOCW(L,K+1,NS)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              endif
            enddo
          elseif( ISTOC(NT) == 0 )then
            ! *** Kd APPROACH
            do NS = 1,NSED2
              if( ITXPARW(NS,NT) == 0 )then
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND) 
                    ! *** M/S          M/S            G/M2/S       M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SEDF(L,K,NS)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              elseif( ITXPARW(NS,NT) == 1 )then
                TMPEXP = CONPARW(NS,NT)
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    TMPVAL = 1.
                    if( SED(L,K+1,NS) > 0.) TMPVAL = SED(L,K+1,NS)**TMPEXP
                    ! *** M/S          M/S                   G/M2/S       M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SEDF(L,K,NS)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              endif
            enddo
          endif            
        endif

        ! *** PARTICLE NONCOHESIVE FLUX
        if( ISTRAN(7) >= 1 )then
          if( ISTOC(NT) > 1 )then
            ! *** fPOC BASED
            do NX = 1,NSND
              NS = NX + NSED2
              if( ITXPARW(NS,NT) == 0 )then
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SNDF(L,K,NX)*STFPOCW(L,K+1,NS)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              endif
              if( ITXPARW(NS,NT) == 1 )then
                TMPEXP = CONPARW(NS,NT)
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    TMPVAL = SND(L,K+1,NX)**TMPEXP
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SNDF(L,K,NX)*STFPOCW(L,K+1,NS)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              endif
            enddo
          elseif( ISTOC(NT) == 0 )then
            ! *** Kd APPROACH
            do NX = 1,NSND
              NS = NX + NSED2
              if( ITXPARW(NS,NT) == 0 )then
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SNDF(L,K,NX)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              endif
              if( ITXPARW(NS,NT) == 1 )then
                TMPEXP = CONPARW(NS,NT)
                do K = 1,KS
                  do LP = 1,LLWET(K,ND)
                    L = LKWET(LP,K,ND)  
                    TMPVAL = SND(L,K+1,NX)**TMPEXP
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SNDF(L,K,NX)*TOXPARW(L,NS,NT)
                  enddo
                enddo
              endif
            enddo
          endif
        endif
      enddo   ! *** NTOX
    endif    ! *** KC >= 2
  enddo
  !$OMP END DO
  
  if( KC >= 2 )then
    !$OMP DO PRIVATE(ND, LP, L, K, NT)
    do ND = 1,NDM
      do NT = 1,NTOX
        do K = 1,KS
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOXF(L,K,NT) = TOXF(L,K,NT)/(1. + TOXPFTW(L,K+1,NT))
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
  endif    ! *** KC >= 2

  ! ******************************************************************************
  ! *** DEPOSITIONAL TOXICS FLUX (TOXF) AT THE SEDIMENT BED/WATER COLUMN INTERFACE 
  ! *** FROM THE WATER COLUMMN ONLY
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT, NS, NX, TMPEXP, TMPVAL, SNDFEFF)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)
    
    do NT = 1,NTOX
      ! *** PARTICLE COHESIVE DEPOSITIONAL FLUX,  (SEDF) < 0 
      if( ISTRAN(6) >= 1 )then
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do NS = 1,NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do LP = LF,LL
                L = LSED(LP)
                ! *** M/S         M/S                G/M2/S                                M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SEDF(L,0,NS),0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(L,NS,NT)
              enddo
            endif
                
            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do LP = LF,LL
                L = LSED(LP)
                TMPVAL = 1.
                if( SED(L,KSZ(L),NS) > 0.) TMPVAL = SED(L,KSZ(L),NS)**TMPEXP
                ! *** M/S         M/S                     G/M2/S                               M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SEDF(L,0,NS),0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(L,NS,NT)
              enddo
            endif
          enddo
        elseif( ISTOC(NT) == 0 )then
          ! *** Kd APPROACH
          do NS = 1,NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do LP = LF,LL
                L = LSED(LP)
                ! *** M/S         M/S                G/M2/S           M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SEDF(L,0,NS),0.)*TOXPARW(L,NS,NT)
              enddo
            endif
            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do LP = LF,LL
                L = LSED(LP)
                TMPVAL = 1.
                if( SED(L,KSZ(L),NS) > 0.) TMPVAL = SED(L,KSZ(L),NS)**TMPEXP
                ! *** M/S         M/S                     G/M2/S           M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SEDF(L,0,NS),0.)*TOXPARW(L,NS,NT)
              enddo
            endif
          enddo
        endif
      endif

      ! *** PARTICLE NONCOHESIVE DEPOSITIONAL FLUX,  (SNDF-SNDFBL) < 0 
      if( ISTRAN(7) >= 1 )then
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do NX = 1,NSND
            NS = NX + NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do LP = LF,LL
                L = LSED(LP)
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SNDFEFF,0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(L,NS,NT)
              enddo
            endif
            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do LP = LF,LL
                L = LSED(LP)
                TMPVAL = SND(L,KSZ(L),NX)**TMPEXP
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SNDFEFF,0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(L,NS,NT)
              enddo
            endif
          enddo
        elseif( ISTOC(NT) == 0 )then
          ! *** Kd APPROACH
          do NX = 1,NSND
            NS = NX + NSED2
            if( ITXPARW(NS,NT) == 0 )then
              do LP = LF,LL
                L = LSED(LP)
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SNDFEFF,0.)*TOXPARW(L,NS,NT)
              enddo
            endif
            if( ITXPARW(NS,NT) == 1 )then
              TMPEXP = CONPARW(NS,NT)
              do LP = LF,LL
                L = LSED(LP)
                TMPVAL = SND(L,KSZ(L),NX)**TMPEXP
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SNDFEFF,0.)*TOXPARW(L,NS,NT)
              enddo
            endif
          enddo

        endif
      endif
    enddo   ! *** NTOX
  enddo
  !$OMP END DO 
  
  ! *** FINALIZE THE DEPOSITIONAL TOXICS FLUX (M/S)
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, NT)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)
    
    do NT = 1,NTOX
      do LP = LF,LL
        L = LSED(LP)
        TOXF(L,0,NT) = TOXF(L,0,NT)/(1. + TOXPFTW(L,KSZ(L),NT))
      enddo
    enddo
  enddo
  !$OMP END DO 
  
  ! ******************************************************************************
  ! *** EROSIONAL TOXICS FLUX (TOXFB) AT THE SEDIMENT BED/WATER COLUMN INTERFACE
  ! *** WHEN SOLIDS FLUX IS > 0
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, NT, NS, NX, TMPVAL, SNDFEFF)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)
    
    do NT = 1,NTOX
      do LP = LF,LL
        L = LSED(LP)
        TOXFB(L,NT) = 0.   ! *** INTERIM UNITS: M/S, FINAL UNITS: 1/S
      enddo
    enddo

    do NT = 1,NTOX
      ! *** PARTICLE COHESIVE EROSIONAL FLUX,  (SEDF) > 0
      if( ISTRAN(6) >= 1 )then
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do NS = 1,NSED2
            do LP = LF,LL
              L = LSED(LP)
              ! *** M/S          M/S             G/M2/S                               M3/G
              TOXFB(L,NT) = TOXFB(L,NT) + MAX(SEDF(L,0,NS),0.)*STFPOCB(L,KBT(L),NS)*TOXPARB(L,NS,NT)
            enddo
          enddo
        elseif( ISTOC(NT) == 0 )then
          ! *** Kd APPROACH
          do NS = 1,NSED2
            do LP = LF,LL
              L = LSED(LP)
              ! ***  M/S         M/S                  G/M2/S          M3/G
              TOXFB(L,NT) = TOXFB(L,NT) + MAX(SEDF(L,0,NS),0.)*TOXPARB(L,NS,NT)
            enddo
          enddo
        endif
      endif
      
      ! *** PARTICLE NONCOHESIVE EROSIONAL FLUX,  (SNDF-SNDFBL) > 0
      if( ISTRAN(7) >= 1 )then
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do NX = 1,NSND
            NS = NX + NSED2
            do LP = LF,LL
              L = LSED(LP)
              SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
              TOXFB(L,NT) = TOXFB(L,NT) + MAX(SNDFEFF,0.)*STFPOCB(L,KBT(L),NS)*TOXPARB(L,NS,NT)
            enddo
          enddo
        elseif( ISTOC(NT) == 0 )then
          ! *** Kd APPROACH
          do NX = 1,NSND
            NS = NX + NSED2
            do LP = LF,LL
              L = LSED(LP)
              SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
              TOXFB(L,NT) = TOXFB(L,NT) + MAX(SNDFEFF,0.)*TOXPARB(L,NS,NT)
            enddo
          enddo
        endif      
      endif
    enddo

    ! *** FINALIZE EROSIONAL FLUX TOXFB (1/S)
    do NT = 1,NTOX
      do LP = LF,LL
        L = LSED(LP)
        if( HBED(L,KBT(L)) > 1E-7 )then
          ! ***                   M                           M
          TMPVAL = PORBED(L,KBT(L))*HBED(L,KBT(L)) + TOXPFTB(L,KBT(L),NT)    ! *** UNITS: M
          ! *** 1/S     M/S           M
          TOXFB(L,NT) = TOXFB(L,NT)/TMPVAL
        else
          TOXFB(L,NT) = 0.
        endif
      enddo
    enddo

  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! ******************************************************************************
  ! *** BANK EROSION
  if( ISBKERO >= 1 )then
    !$OMP SINGLE
    do NT = 1,NTOX
      do L = 2,LA
        TOXFBEBKB(L,NT) = 0.
        TOXFBECHB(L,NT) = 0.
        TOXFBECHW(L,NT) = 0.
      enddo
    enddo

    do NT = 1,NTOX
      ! PARTICLE cohesive flux
      if( ISTRAN(6) >= 1 )then
        do NS = 1,NSED
          do NP = 1,NBEPAIR
            LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
            KBANK = KBT(LBANK)
            TOXFBECHB(NP,NT) = TOXFBECHB(NP,NT) + SEDFBECHB(NP,NS)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(LBANK,NS,NT)
            TOXFBECHW(NP,NT) = TOXFBECHW(NP,NT) + SEDFBECHW(NP,NS)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(LBANK,NS,NT)
            TOXFBEBKB(NP,NT) = TOXFBEBKB(NP,NT) + SEDFBEBKB(NP,NS)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(LBANK,NS,NT)
          enddo
        enddo
      endif

      ! PARTICLE noncohesive flux
      if( ISTRAN(7) >= 1 )then
        do NX = 1,NSND
          NS = NX + NSED
          do NP = 1,NBEPAIR
            LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
            KBANK = KBT(LBANK)
            TOXFBEBKB(NP,NT) = TOXFBEBKB(NP,NT) + SNDFBEBKB(NP,NX)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(LBANK,NS,NT)
            TOXFBECHB(NP,NT) = TOXFBECHB(NP,NT) + SNDFBECHB(NP,NX)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(LBANK,NS,NT)
            TOXFBECHW(NP,NT) = TOXFBECHW(NP,NT) + SNDFBECHW(NP,NX)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(LBANK,NS,NT)
          enddo
        enddo
      endif
    enddo

    do NT = 1,NTOX
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        KBANK = KBT(LBANK)
        if( HBED(LBANK,KBANK) > 0.0 )then
          TOXFBEBKB(NP,NT) = TOXFBEBKB(NP,NT)/(PORBED(LBANK,KBANK)*HBED(LBANK,KBANK) + TOXPFTB(LBANK,KBANK,NT))
          TOXFBECHB(NP,NT) = TOXFBECHB(NP,NT)/(PORBED(LBANK,KBANK)*HBED(LBANK,KBANK) + TOXPFTB(LBANK,KBANK,NT))
          TOXFBECHW(NP,NT) = TOXFBECHW(NP,NT)/(PORBED(LBANK,KBANK)*HBED(LBANK,KBANK) + TOXPFTB(LBANK,KBANK,NT))
        else
          TOXFBEBKB(NP,NT) = 0.
          TOXFBECHB(NP,NT) = 0.
          TOXFBECHW(NP,NT) = 0.
        endif
      enddo
    enddo
    !$OMP END SINGLE
  endif 
  ! *** END BANK EROSION

  !$OMP DO PRIVATE(ND,LF,LL,NT,K,LP,L)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)

    !**********************************************************************CC
    ! ***  CALCULATE TOTAL PARTICULATE FRACTIONS IN WATER COLUMN AND BED
    ! ***  NOTING THAT TO THIS POINT TOXPFTW AND TOXPFTB HAVE BEEN USED
    ! ***  TO TEMPORARILY STORE THE SORBED PORTION
    do NT = 1,NTOX
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***                                                          SORBED TO DOC
          TOXPFTW(L,K,NT) = ( TOXPFTW(L,K,NT)/(1. + TOXPFTW(L,K,NT)) ) - TOXPFW(L,K,NFD,NT)
          TOXPFTW(L,K,NT) = MIN(TOXPFTW(L,K,NT),SORBMIN)
        enddo
      enddo
    enddo
    
    ! *** COMPUTE MASS FRACTION (TOXPFTB) OF SORBED TOXICS TO PARTICULATES (DIMENSIONLESS)
    do NT = 1,NTOX
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          if( HBED(L,K) >= HBEDMIN0 .and. PORBED(L,K) > 0. )then
            ! ***                                                                                 SORBED TO DOC
            TOXPFTB(L,K,NT) = ( TOXPFTB(L,K,NT)/( PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT) ) ) - TOXPFB(L,K,NFD,NT)
            TOXPFTB(L,K,NT) = MAX(TOXPFTB(L,K,NT),SORBMIN)
          else
            TOXPFTB(L,K,NT) = SORBMIN
          endif
        enddo
      enddo
    enddo
  enddo       ! *** END OF DOMAIN LOOP
  !$OMP END DO

  !**********************************************************************C
  ! ***  DETERMINE TOXIC FLUX FROM BED LOAD SORBED MATERIAL
  if( ICALC_BL > 0 .and. .not. LSEDZLJ )then  
    !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT, NS, NX, LE, LW, LN, LS) &
    !$OMP    PRIVATE(TMPTOXC, TMPTOXE, TMPTOXW, TMPTOXN, TMPTOXS)
    do ND = 1,NDM
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED - 1,LASED)

      do NT = 1,NTOX
        do LP = LF,LL
          L = LSED(LP)
          TOXFBL(L,NT) = 0.0     ! *** UNITS:  MG/M2/S
        enddo
      enddo

      do NT = 1,NTOX
        do NX = 1,NSND
          NS = NX + NSED2
          do LP = LF,LL
            L = LSED(LP)
            LE = LEC(L)
            LW = LWC(L)
            LN = LNC(L)
            LS = LSC(L)
            TMPTOXC = 0.0
            TMPTOXE = 0.0
            TMPTOXW = 0.0
            TMPTOXN = 0.0
            TMPTOXS = 0.0
            K = KBT(L)
            if( SNDB1(L,K,NX) > 0.0 )then
              ! ***         DIMENSIONLESS    MG/M2        M          M2/G      (ORIGINAL CODE)
              !TMPTOXC = TOXPFB(L,K,NS,NT)*TOXB(L,K,NT)*HBED(L,K)/SNDB(L,K,NX)
              ! MG/G         DIMENSIONLESS     MG/M2        M2/G 
              TMPTOXC = TOXPFB(L,K,NS,NT)*TOXB(L,K,NT)/SNDB1(L,K,NX) 
            endif
            K = KBT(LE)
            if( SUB(LE) > 0.5 .and. SNDB1(LE,K,NX) > 0.0 )then
              TMPTOXE = TOXPFB(LE,K,NS,NT)*TOXB(LE,K,NT)/SNDB1(LE,K,NX)
            endif
            K = KBT(LW)
            if( SUB(L) > 0.5 .and. SNDB1(LW,K,NX) > 0.0 )then
              TMPTOXW = TOXPFB(LW,K,NS,NT)*TOXB(LW,K,NT)/SNDB1(LW,K,NX)
            endif
            K = KBT(LN)
            if( SVB(LN) > 0.5 .and. SNDB1(LN,K,NX) > 0.0 )then
              TMPTOXN = TOXPFB(LN,K,NS,NT)*TOXB(LN,K,NT)/SNDB1(LN,K,NX)
            endif
            K = KBT(LS)
            if( SVB(L) > 0.5 .and. SNDB1(LS,K,NX) > 0.0 )then
              TMPTOXS = TOXPFB(LS,K,NS,NT)*TOXB(LS,K,NT)/SNDB1(LS,K,NX)
            endif

            ! ***             MG/M2/S       1/M2                        G/S           MG/G
            TOXFBL(L,NT) = TOXFBL(L,NT) + DXYIP(L)*( SUB(LE)*MAX(QSBDLDX(LE,NX),0.)*TMPTOXC + SUB(LE)*MIN(QSBDLDX(LE,NX),0.)*TMPTOXE  &
                                                   - SUB(L) *MAX(QSBDLDX(L,NX),0.) *TMPTOXW - SUB(L) *MIN(QSBDLDX(L,NX),0.) *TMPTOXC  &
                                                   + SVB(LN)*MAX(QSBDLDY(LN,NX),0.)*TMPTOXC + SVB(LN)*MIN(QSBDLDY(LN,NX),0.)*TMPTOXN  &
                                                   - SVB(L) *MAX(QSBDLDY(L,NX),0.) *TMPTOXS - SVB(L) *MIN(QSBDLDY(L,NX),0.) *TMPTOXC )
          enddo
        enddo
      enddo   ! *** NTOX
    enddo       ! *** END OF DOMAIN LOOP
    !$OMP END DO

    !$OMP SINGLE
    ! *** ADJUST FOR BEDLOAD TRANSPORT AT BOUNDARY CELLS
    if( NSBDLDBC > 0 )then
      ! *** UPSTREAM AND DOWNSTREAM 
      do NSB = 1,NSBDLDBC
        LUTMP = LSBLBCU(NSB)
        LDTMP = LSBLBCD(NSB)
        do NT = 1,NTOX
          do NX = 1,NSND
            NS = NX + NSED2
            TMPTOXC = 0.0
            K = KBT(LUTMP)
            if( SNDB1(LUTMP,K,NX) > 0.0 )then
              TMPTOXC = TOXPFB(LUTMP,K,NS,NT)*TOXB(LUTMP,K,NT)/SNDB1(LUTMP,K,NX)
            endif  

            ! *** MG/M2/S            MG/M2/S          1/M2          G/S            MG/G          
            TOXFBL(LUTMP,NT) = TOXFBL(LUTMP,NT) - DXYIP(LUTMP)*QSBDLDOT(LUTMP,NX)*TMPTOXC
            if( LDTMP /= 0 ) TOXFBL(LDTMP,NT) = TOXFBL(LDTMP,NT) + DXYIP(LDTMP)*QSBDLDOT(LUTMP,NX)*TMPTOXC
          enddo
        enddo
      enddo
    endif   ! *** END OF NSBDLDBC > 0 BLOCK
    !$OMP END SINGLE
  endif
  
  ! *** HANDLE BEDLOAD FOR SEDZLJ
  if( ICALC_BL > 0 .and. LSEDZLJ )then
    !$OMP DO PRIVATE(ND,LF,LL,NT,LP,L,NS,LE,LW,LN,LS) &
    !$OMP    PRIVATE(TMPTOXB,TMPTOXC,TMPTOXE,TMPTOXW,TMPTOXN,TMPTOXS,TMPVAL,TOXFLUX,CBLTOXTMP)
    do ND = 1,NDM
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED-1,LASED)

      do NT = 1,NTOX
        do LP = LF,LL
          L = LSED(LP)
          TOXFBL(L,NT) = 0.0     ! *** UNITS:  MG/M2/S  (+ is Erosional, - is Depositional)
        enddo
      enddo

      ! *** Compute bedload toxics concentration
      do NT = 1,NTOX
        do NS = 1,NSEDS2
          if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
          do LP = LF,LL
            L = LSED(LP)
            LE = LEC(L)
            LW = LWC(L)
            LN = LNC(L)
            LS = LSC(L)
            if( SEDB1(L,KBT(L),NS) > 1.E-6 )then
              ! MG/G          DIMENSIONLESS          MG/M2            M2/G
              TMPTOXB = TOXPFB(L,KBT(L),NS,NT)*TOXB1(L,KBT(L),NT)/SEDB1(L,KBT(L),NS) 
            else
              TMPTOXB = 0.0
            endif
            ! MG/G                 DIMENSIONLESS           MG/G
            TMPTOXC =         TOXPFB(L,KBT(L),NS,NT)  *CBLTXCON(L,NS,NT)
            TMPTOXE = SUB(LE)*TOXPFB(LE,KBT(LE),NS,NT)*CBLTXCON(LE,NS,NT)
            TMPTOXW = SUB(L)* TOXPFB(LW,KBT(LW),NS,NT)*CBLTXCON(LW,NS,NT)
            TMPTOXN = SVB(LN)*TOXPFB(LN,KBT(LN),NS,NT)*CBLTXCON(LN,NS,NT)
            TMPTOXS = SVB(L)* TOXPFB(LS,KBT(LS),NS,NT)*CBLTXCON(LS,NS,NT)

            ! *** NET TOXIC FLUX IN/OUT OF THE BED (MG/M2)
            TOXFLUX = EBL(L,NS)*TMPTOXB - DBL(L,NS)*TMPTOXC

            ! *** CONCENTRATION OF TOXIC IN BEDLOAD LAYER (MG/M2)
            TMPVAL = DTSED*DXYIP(L)
            CBLTOXTMP = CBLTOX(L,NT) + TMPVAL*(  MAX(QSBDLDX(L,NS),0.)*TMPTOXW - MIN(QSBDLDX(LE,NS),0.)*TMPTOXE   &   ! *** EAST/WEST: IN 
                                               + MIN(QSBDLDX(L,NS),0.)*TMPTOXC - MAX(QSBDLDX(LE,NS),0.)*TMPTOXC   &   ! *** EAST/WEST: OUT 
                                               + MAX(QSBDLDY(L,NS),0.)*TMPTOXS - MIN(QSBDLDY(LN,NS),0.)*TMPTOXN   &   ! *** NORTH/SOUTH: IN
                                               + MIN(QSBDLDY(L,NS),0.)*TMPTOXC - MAX(QSBDLDY(LN,NS),0.)*TMPTOXC ) &   ! *** NORTH/SOUTH: OUT
                                               + TOXFLUX
            
            if( CBLTOXTMP >= 0. )then
              CBLTOX(L,NT) = CBLTOXTMP
            else
              ! *** LIMIT FLUX TO AVAILABLE MASS
              CBLTOX(L,NT) = 0.0
              TOXFBL(L,NT) = TOXFBL(L,NT) - CBLTOXTMP*DELTI
              EXIT
            endif
            
            TOXFBL(L,NT) = TOXFBL(L,NT) + TOXFLUX*DELTI                                                                   ! *** TOXFBL FLUX RATE   (MG/M2/S)
          enddo
        enddo
      enddo   ! *** NTOX

      ! *** CHECK FOR DEPLETED BEDLOAD MASS
      do LP = LF,LL
        L = LSED(LP)
        if( SUM(CBL(L,1:NSEDS)) <= 1.E-8 )then
          ! *** MOVE ANY REMAINING MASS INTO BED
          TOXB(L,KBT(L),1:NTOX) = TOXB(L,KBT(L),1:NTOX) + CBLTOX(L,1:NTOX)
          CBLTOX(L,1:NTOX) = 0.0
        endif
      enddo
      
    enddo       ! *** END OF DOMAIN LOOP
    !$OMP END DO

    ! *** ADJUST FOR BEDLOAD TRANSPORT AT BOUNDARY CELLS
    if( NSBDLDBC > 0 )then
      !$OMP SINGLE
      ! *** UPSTREAM AND DOWNSTREAM 
      do NSB = 1,NSBDLDBC
        LUTMP = LSBLBCU(NSB)
        LDTMP = LSBLBCD(NSB)
        do NT = 1,NTOX
          TOXFBL(LUTMP,NT) = 0.0
          do NS = 1,NSEDS
            if( D50(NS) < BEDLOAD_CUTOFF ) CYCLE
            TMPTOXC = TOXPFB(LUTMP,KBT(LUTMP),NS,NT)*CBLTXCON(LUTMP,NS,NT)

            ! *** MG/M2/S            MG/M2/S          1/M2          G/S            MG/G
            !TOXFBL(LUTMP,NT) = TOXFBL(LUTMP,NT) - DXYIP(LUTMP)*QSBDLDOT(LUTMP,NS)*TMPTOXC
            
            ! *** MG/M2              MG/M2/S          1/M2          G/S            MG/G    S
            !CBLTOX(LUTMP,NT) = CBLTOX(LUTMP,NT) - DXYIP(LUTMP)*QSBDLDOT(LUTMP,NS)*TMPTOXC*DELT
            if( LDTMP /= 0 )then
              ! *** MG/M2/S            MG/M2/S          1/M2          G/S            MG/G
              TOXFBL(LDTMP,NT) = TOXFBL(LDTMP,NT) - DXYIP(LDTMP)*QSBDLDOT(LUTMP,NS)*TMPTOXC
            endif
          enddo
        enddo
      enddo
      !$OMP END SINGLE
    endif   ! *** END OF NSBDLDBC > 0 BLOCK
  endif     ! *** END OF ICALC_BL > 0 BLOCK

  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,NT,NS,NX,TMPVAL)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED - 1,LASED)

    !**********************************************************************C
    ! ***  ADJUST TOXIC FLUXES ACROSS WATER COLUMN - BED INTERFACE TO
    ! ***  INCLUDE WATER ENTRAINMENT AND EXPULSION ASSOCIATED WITH
    ! ***  DEPOSITION AND RESUSPENSION
    
    ! *** DEPSITIONAL FLUX ADJUSTMENT DUE TO ENTRAPMENT (M/S)
    do NT = 1,NTOX
      do LP = LF,LL
        L = LSED(LP)
        TMPVAL = ( MIN(QSBDTOP(L),0.0) + MIN(QWBDTOP(L),0.0) )    ! *** TOTAL DEPOSITIONAL VOLUMETRIC RATE (M/S)
        TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*(1.-TOXPFTW(L,KSZ(L),NT))
      enddo
    enddo
  
    ! *** EROSONAL FLUX ADJUSTMENT DUE TO POREWATER EXPULSION (1/S)
    do NT = 1,NTOX
      do LP = LF,LL
        L = LSED(LP)
        K = KBT(L)
        if( HBED(L,K) > 1.E-6 )then
          TMPVAL = ( MAX(QSBDTOP(L),0.0) + MAX(QWBDTOP(L),0.0) )/HBED(L,K)      ! *** TOTAL EROSIONAL VOLUMETRIC RATE (M/S)/BED THICKNESS (M)
        else
          TMPVAL = 0.
        endif       
        TOXFB(L,NT) = TOXFB(L,NT) + TMPVAL*(1.-TOXPFTB(L,K,NT))
      enddo
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! *** START BANK EROSION
  if( ISBKERO >= 1 )then
    !$OMP SINGLE
    do NT = 1,NTOX
      do NP = 1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        KBANK = KBT(LBANK)
        TMPVAL = (1.-TOXPFTB(LBANK,KBANK,NT))/HBED(LBANK,KBANK)
        TOXFBEBKB(NP,NT) = TOXFBEBKB(NP,NT) + TMPVAL*(QSBDTOPBEBKB(NP) + QWBDTOPBEBKB(NP))
        TOXFBECHB(NP,NT) = TOXFBECHB(NP,NT) + TMPVAL*(QSBDTOPBECHB(NP) + QWBDTOPBECHB(NP))
        TOXFBECHW(NP,NT) = TOXFBECHW(NP,NT) + TMPVAL*(QSBDTOPBECHW(NP) + QWBDTOPBECHW(NP))
      enddo
    enddo
    !$OMP END SINGLE
  endif
  !---END BANK EROSION-------------------------

  ! ***********************************************************************************
  ! ***  UPDATE WATER COLUMN BOTTOM LAYER AND TOP SEDIMENT LAYER FOR TOXIC CONTAMINANT
  ! ***  PMC - NEW APPROACH FOR BETTER MASS BALANCE 2017-07
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT, NS, NX, LE, LW, LN, LS)  &
  !$OMP    PRIVATE(AA11, BB11, BB22, TMPVAL, CLEFT, CRIGHT, WVEL)
  do ND = 1,NDM
    LF = (ND-1)*LDMWET + 1  
    LL = MIN(LF + LDMWET-1,LAWET)

    ! ----------------------------------------------------------------
    ! *** KC = 1 (SINGLE LAYER IN VERTICAL)
    if( KC == 1 )then
      do NT = 1,NTOX
        ! *** WATER COLUMN/SEDIMENT BED INTERFACE
        do LP = LF,LL
          L = LWET(LP)
          
          ! *** BB11 = 0.0                                  ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          ! *** BB22 = 0.0                                  ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S
          if( LBED(L) ) CYCLE        
          
          ! *** DETERMINE FLUXES
          BB11 = TOXFB(L,NT)*TOXB(L,KBT(L),NT)              ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          BB22 = TOXF(L,0,NT)*TOX(L,KSZ(L),NT)              ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S

          ! *** LIMIT DEPOSITING MASS TO AVAILABLE MASS
          TMPVAL = TOX(L,KSZ(L),NT) + DTSED*BB22*HPKI(L,KSZ(L))
          if( TMPVAL < 0.0 )then
            BB22 = -TOX(L,KSZ(L),NT)*DELTI*HPK(L,KSZ(L))
          endif
            
          ! *** CHECK FOR AVAILABLE MASS IN THE BED
          TMPVAL = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          if( TMPVAL < 0.0 )then
            BB11 = TOXB(L,KBT(L),NT)*DELTI - BB22 - TOXFBL(L,NT)
            BB11 = MAX(BB11,0.0)
          endif
            
          ! *** MASS BALANCE AROUND WATER COLUMN BOTTOM LAYER
          TOX(L,KSZ(L),NT) = TOX(L,KSZ(L),NT) + DTSED*( BB11 + BB22 )*HPKI(L,KSZ(L))
          if( TOX(L,KSZ(L),NT) < 0.0 )then
            open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            write(mpi_efdc_out_unit, '(A,F15.6,I7,I5,E14.6,I5,2E14.6)' ) 'Bad TOX:  TIMEDAY, L, KBT, HBED, NT, TOX1,  TOX:  ', TIMEDAY, MAP2GLOBAL(L).LG, KBT(L), HBED(L,KBT(L)), NT, TOX1(L,KSZ(L),NT) , TOX(L,KSZ(L),NT)
            close(mpi_efdc_out_unit)
            PRINT '(A,F15.6,I7,I5,E14.6,I5,2E14.6)', 'Bad TOX', TIMEDAY, MAP2GLOBAL(L).LG, KBT(L), HBED(L,KBT(L)), NT, BB22, TOX(L,KSZ(L),NT) 
            TOX(L,KSZ(L),NT) = 0.0
          endif
          
          ! *** MASS BALANCE AROUND TOP SEDIMENT LAYER
          TOXB1(L,KBT(L),NT) = TOXB(L,KBT(L),NT)
          TOXB(L,KBT(L),NT) = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          if( TOXB(L,KBT(L),NT) < 0.0 )then
            !OPEN(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            !WRITE(mpi_efdc_out_unit, '(A,F15.6,I7,I5,E14.6,I5,2E14.6)' ) 'Bad TOXB: TIMEDAY, L, KBT, HBED, NT, TOXB1, TOXB: ', TIMEDAY, MAP2GLOBAL(L).LG, KBT(L), HBED(L,KBT(L)), NT, TOXB1(L,KBT(L),NT), TOXB(L,KBT(L),NT)
            !CLOSE(mpi_efdc_out_unit)
            !PRINT '(A,F15.6,I7,I5,E14.6,I5,2E14.6)', 'Bad TOXB', TIMEDAY, MAP2GLOBAL(L).LG, KBT(L), HBED(L,KBT(L)), NT, TOXB1(L,KBT(L),NT), TOXB(L,KBT(L),NT) 
            !IF( ABS(TOXB(L,KBT(L),NT)) > 1E-5*TOXB1(L,KBT(L),NT) )then
            !  PAUSE  
            !ENDIF
            TOXB(L,KBT(L),NT) = 0.0 
          endif          
        enddo
      enddo   ! *** NTOX
    endif     ! *** KC = 1

    ! ----------------------------------------------------------------
    ! *** KC = 2 (TWO LAYERS IN VERTICAL)
    if( KC == 2 )then
      do NT = 1,NTOX
        ! *** WATER COLUMN
        K = KC
        do LP = LF,LL
          L = LWET(LP)
          WVEL = DELTI*HPK(L,K)
          CLEFT = WVEL-TOXF(L,K-1,NT)
          CRIGHT = WVEL*TOX(L,K,NT)
          TOX(L,K,NT) = CRIGHT/CLEFT
        enddo
       
        ! *** WATER COLUMN/SEDIMENT BED INTERFACE
        do LP = LF,LL
          L = LWET(LP)

          ! *** BB11 = 0.0                                  ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          ! *** BB22 = 0.0                                  ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S
          if( LBED(L) ) CYCLE        
          
          ! *** DETERMINE FLUXES
          AA11 = TOXF(L,KSZ(L),NT)*TOX(L,KSZ(L)+1,NT)       ! *** SETTLING FROM TOP OF BOTTOM LAYER,        MG/M2/S
          BB11 = TOXFB(L,NT)*TOXB(L,KBT(L),NT)              ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          BB22 = TOXF(L,0,NT)*TOX(L,KSZ(L),NT)              ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S

          ! *** LIMIT DEPOSITING MASS TO AVAILABLE MASS
          TMPVAL = TOX(L,KSZ(L),NT) + DTSED*( BB22 - AA11 )*HPKI(L,KSZ(L))
          if( TMPVAL < 0.0 )then
            BB22 = -TOX(L,KSZ(L),NT)*DELTI*HPK(L,KSZ(L)) - AA11
          endif
            
          ! *** CHECK FOR AVAILABLE MASS IN THE BED
          TMPVAL = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          if( TMPVAL < 0.0 )then
            BB11 = TOXB(L,KBT(L),NT)*DELTI - BB22 - TOXFBL(L,NT)
            if( BB11 < 0.0 )then
              !PRINT '(A,F10.6,I7,I5,2E14.6)', 'Bad TOXB correction', TIMEDAY, L, KBT(L), BB11, TOXB(L,KBT(L),NT)   ! DELME
              !pause
              BB11 = 0.0   ! delme
            endif
          endif
            
          ! *** MASS BALANCE AROUND WATER COLUMN BOTTOM LAYER
          TOX(L,KSZ(L),NT) = TOX(L,KSZ(L),NT) + DTSED*( BB11 + BB22 - AA11 )*HPKI(L,KSZ(L))
          if( TOX(L,KSZ(L),NT) < 0.0 )then
            !PRINT '(A,F10.6,I7,I5,2E14.6)', 'Bad TOX correction', TIMEDAY, L, KSZ(L), BB22, TOX(L,KSZ(L),NT)   ! DELME
            !pause
            TOX(L,KSZ(L),NT) = 0.0
          endif
          
          ! *** MASS BALANCE AROUND TOP SEDIMENT LAYER
          TOXB1(L,KBT(L),NT) = TOXB(L,KBT(L),NT)
          TOXB(L,KBT(L),NT) = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          if( TOXB(L,KBT(L),NT) < 0.0 )then
            !OPEN(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            !WRITE(mpi_efdc_out_unit, '(A,F10.6,I7,I5,E14.6,I5,2E14.6)' ) 'Bad TOXB: TIMEDAY, L, KBT, HBED, NT, TOXB1, TOXB: ', TIMEDAY, L, KBT(L), HBED(L,KBT(L)), NT, TOXB1(L,KBT(L),NT), TOXB(L,KBT(L),NT)
            !CLOSE(mpi_efdc_out_unit)
            !PRINT '(A,F10.6,I7,I5,E14.6,I5,2E14.6)', 'Bad TOXB', TIMEDAY, L, KBT(L), HBED(L,KBT(L)), NT, TOXB1(L,KBT(L),NT), TOXB(L,KBT(L),NT)   ! DELME
            !IF( ABS(TOXB(L,KBT(L),NT)) > 1E-5*TOXB1(L,KBT(L),NT) )then
            !  pause   ! delme
            !ENDIF
            TOXB(L,KBT(L),NT) = 0.0    !  delme
          endif   
        enddo
      enddo   ! *** NTOX
    endif     ! *** KC = 2

    ! ----------------------------------------------------------------
    ! *** KC >= 3 (THREE OR MORE LAYERS IN VERTICAL)
    if( KC >= 3 )then
      do NT = 1,NTOX
        ! *** WATER COLUMN
        K = KC
        do LP = LF,LL
          L = LWET(LP)
          WVEL = DELTI*HPK(L,K)
          CLEFT = WVEL-TOXF(L,K-1,NT)
          CRIGHT = WVEL*TOX(L,K,NT)
          TOX(L,K,NT) = CRIGHT/CLEFT
        enddo
       
        do K = KS,2,-1 
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND)  
            WVEL = DELTI*HPK(L,K)
            CLEFT = WVEL - TOXF(L,K-1,NT)                                 ! *** M/S
            CRIGHT = WVEL*TOX(L,K,NT) - TOXF(L,K,NT)*TOX(L,K+1,NT)        ! *** MG/M2/S
            TOX(L,K,NT) = CRIGHT/CLEFT                                    ! *** MG/M3
          enddo
        enddo

        ! *** HANDLE HARD BOTTOM BOTTOM LAYER SETTLING
        if( ISBEDMAP > 0 )then
          do LP = LF,LL
            L = LWET(LP)  
            if( LBED(L) )then
              ! *** ADD SETTLING FROM THE LAYER ABOVE
              K = KSZ(L)
              WVEL = DELTI*HPK(L,K)
              CLEFT = WVEL                                                ! *** M/S
              CRIGHT = WVEL*TOX(L,K,NT) - TOXF(L,K,NT)*TOX(L,K+1,NT)      ! *** MG/M2/S
              TOX(L,K,NT) = CRIGHT/CLEFT                                  ! *** MG/M3
            endif
          enddo
        endif
      
        ! *** WATER COLUMN/SEDIMENT BED INTERFACE
        do LP = LF,LL
          L = LWET(LP)

          ! *** BB11 = 0.0                                  ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          ! *** BB22 = 0.0                                  ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S
          if( LBED(L) ) CYCLE        

          ! *** DETERMINE FLUXES
          AA11 = TOXF(L,KSZ(L),NT)*TOX(L,KSZ(L)+1,NT)       ! *** SETTLING FROM TOP OF BOTTOM LAYER,        MG/M2/S
          BB11 = TOXFB(L,NT)*TOXB(L,KBT(L),NT)              ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          BB22 = TOXF(L,0,NT)*TOX(L,KSZ(L),NT)              ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S

          ! *** LIMIT DEPOSITING MASS TO AVAILABLE MASS
          TMPVAL = TOX(L,KSZ(L),NT) + DTSED*( BB22 - AA11 )*HPKI(L,KSZ(L))
          if( TMPVAL < 0.0 )then
            !BB22 = (TOX(L,KSZ(L),NT) - TMPVAL)*DELTI*HPK(L,KSZ(L)) - AA11
            BB22 = -TOX(L,KSZ(L),NT)*DELTI*HPK(L,KSZ(L)) - AA11
          endif
            
          ! *** CHECK FOR AVAILABLE MASS IN THE BED
          TMPVAL = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          if( TMPVAL < 0.0 )then
            !BB11 = (TOXB(L,KBT(L),NT) - TMPVAL)*DELTI - BB22 - TOXFBL(L,NT)
            BB11 = TOXB(L,KBT(L),NT)*DELTI - BB22 - TOXFBL(L,NT)
            if( BB11 < 0.0 )then
              !PRINT '(A,F10.6,I7,I5,2E14.6)', 'Bad TOXB correction', TIMEDAY, L, KBT(L), BB11, TOXB(L,KBT(L),NT)   ! DELME
              !pause
              BB11 = 0.0   ! delme
            endif
          endif
            
          ! *** MASS BALANCE AROUND WATER COLUMN BOTTOM LAYER
          TOX(L,KSZ(L),NT) = TOX(L,KSZ(L),NT) + DTSED*( BB11 + BB22 - AA11 )*HPKI(L,KSZ(L))
          if( TOX(L,KSZ(L),NT) < 0.0 )then
            !PRINT '(A,F10.6,I7,I5,2E14.6)', 'Bad TOX correction', TIMEDAY, L, KSZ(L), BB22, TOX(L,KSZ(L),NT)   ! DELME
            !pause
            TOX(L,KSZ(L),NT) = 0.0
          endif
          
          ! *** MASS BALANCE AROUND TOP SEDIMENT LAYER
          TOXB1(L,KBT(L),NT) = TOXB(L,KBT(L),NT)
          TOXB(L,KBT(L),NT) = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          if( TOXB(L,KBT(L),NT) < 0.0 )then
            !OPEN(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
            !WRITE(mpi_efdc_out_unit, '(A,2F10.6,I7,I5,E14.6,I5,2E14.6)' ) 'Bad TOXB: TIMEDAY, DELT, L, KBT, HBED, NT, TOXB1, TOXB: ', TIMEDAY, L, KBT(L), HBED(L,KBT(L)), NT, TOXB1(L,KBT(L),NT), TOXB(L,KBT(L),NT)
            !CLOSE(mpi_efdc_out_unit)
            !PRINT '(A,F10.6,I7,I5,E14.6,I5,2E14.6)', 'Bad TOXB', TIMEDAY, L, KBT(L), HBED(L,KBT(L)), NT, TOXB1(L,KBT(L),NT), TOXB(L,KBT(L),NT)   ! DELME
            !IF( ABS(TOXB(L,KBT(L),NT)) > 1E-5*TOXB1(L,KBT(L),NT) )then
            !  pause   ! delme
            !ENDIF
            TOXB(L,KBT(L),NT) = 0.0
          endif   
        enddo
      enddo   ! *** NTOX
    endif     ! *** KC >= 3
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  ! *** START BANK EROSION
  if( ISBKERO >= 1 )then
    !$OMP SINGLE
    do NP = 1,NBEPAIR
      LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
      LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
      KBANK = KBT(LBANK)
      KCHAN = KBT(LCHAN)
      TOXFBEBKB(NP,NT) = TOXFBEBKB(NP,NT)*TOXB(LBANK,KBANK,NT)
      TOXFBECHB(NP,NT) = TOXFBECHB(NP,NT)*TOXB(LBANK,KBANK,NT)
      TOXFBECHW(NP,NT) = TOXFBECHW(NP,NT)*TOXB(LBANK,KBANK,NT)
      TOXB(LBANK,KBANK,NT) = TOXB(LBANK,KBANK,NT)-DTSED*TOXFBEBKB(NP,NT)
      TOXB(LCHAN,KCHAN,NT) = TOXB(LCHAN,KCHAN,NT)-DTSED*TOXFBECHB(NP,NT)
      TOX(LCHAN,KSZ(LCHAN),NT) = TOX(LCHAN,KSZ(LCHAN),NT) + DTSED*HPKI(L,KSZ(LCHAN))*TOXFBECHW(NP,NT)
    enddo
    !$OMP END SINGLE
  endif

  if( .not. LSEDZLJ )then
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,NT,NS,NX,LE,LW,LN,LS,KTOPTP,KTOPM1)  &
    !$OMP    PRIVATE(FTPOS,FTNEG)
    do ND = 1,NDM
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED-1,LASED)
    
      !**********************************************************************C
      !
      ! ***  ADD PARENT TO ACTIVE LAYER TOXIC TRANSPORT
      !
      if( ISNDAL >= 2 )then
        do NT = 1,NTOX
          do LP = LF,LL
            L = LSED(LP)
            TOXFPA(L) = 0.0
          enddo
          if( ISTRAN(6) > 0 )then
            do NS = 1,NSED
              do LP = LF,LL
                L = LSED(LP)
                KTOPTP = KBT(L)
                KTOPM1 = KBT(L)-1
                if( KTOPM1 < 1 ) CYCLE
                FTPOS = 0.0
                FTNEG = 0.0
                if( SEDFPA(L,NS) > 0.0 .and. SEDB(L,KTOPM1,NS) > 0.0 ) FTPOS = SEDFPA(L,NS)*TOXPFB(L,KTOPM1,NS,NT)/SEDB(L,KTOPM1,NS)
                if( SEDFPA(L,NS) < 0.0 .and. SEDB(L,KTOPTP,NS) > 0.0 ) FTNEG = SEDFPA(L,NS)*TOXPFB(L,KTOPTP,NS,NT)/SEDB(L,KTOPTP,NS)
                TOXFPA(L) = TOXFPA(L) + FTPOS*TOXB(L,KTOPM1,NT) + FTNEG*TOXB(L,KTOPTP,NT)
              enddo
            enddo
          endif
          
          if( ISTRAN(7) > 0 )then
            do NX = 1,NSND
              NS = NX + NSED2
              do LP = LF,LL
                L = LSED(LP)
                KTOPTP = KBT(L)
                KTOPM1 = KBT(L)-1
                if( KTOPM1 < 1 ) CYCLE
                FTPOS = 0.0
                FTNEG = 0.0
                if( SNDFPA(L,NX) > 0.0 .and. SNDB(L,KTOPM1,NX) > 0.0 ) FTPOS = SNDFPA(L,NX)*TOXPFB(L,KTOPM1,NS,NT)/SNDB(L,KTOPM1,NX)
                if( SNDFPA(L,NX) < 0.0 .and. SNDB(L,KTOPTP,NX) > 0.0 ) FTNEG = SNDFPA(L,NX)*TOXPFB(L,KTOPTP,NS,NT)/SNDB(L,KTOPTP,NX)
                TOXFPA(L) = TOXFPA(L) + FTPOS*TOXB(L,KTOPM1,NT) + FTNEG*TOXB(L,KTOPTP,NT)
              enddo
            enddo
          endif
          
          do LP = LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            if( KTOPM1 < 1 ) CYCLE
            FTPOS = 0.0
            FTNEG = 0.0
            if( QWATPA(L) > 0.0 ) FTPOS = QSSDPA(L)*(1.-TOXPFTB(L,KTOPM1,NT))
            if( QWATPA(L) < 0.0 ) FTNEG = QSSDPA(L)*(1.-TOXPFTB(L,KTOPTP,NT))
            TOXFPA(L) = TOXFPA(L) + FTPOS*TOXB(L,KTOPM1,NT) + FTNEG*TOXB(L,KTOPTP,NT)
          enddo
          do LP = LF,LL
            L = LSED(LP)
            KTOPTP = KBT(L)
            KTOPM1 = KBT(L)-1
            if( KTOPM1 < 1 ) CYCLE
            TOXB(L,KTOPTP,NT) = TOXB(L,KTOPTP,NT) + DTSED*TOXFPA(L)
            TOXB(L,KTOPM1,NT) = TOXB(L,KTOPM1,NT)-DTSED*TOXFPA(L)
          enddo
        enddo
      endif   ! *** ISNDAL >= 2
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO
  endif

  !$OMP END PARALLEL
    
  return

  
END

