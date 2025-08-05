! ----------------------------------------------------------------------
!   This file is a part of EFDC + 
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALTOXB

  ! *** SUBROUTINE CALTOXB CALCULATES CONTAMINANT TRANSPORT WITHIN THE SEDIMENT BED
  ! *** AND THE FLUX OF CONTAMINANTS AT THE SEDIMENT/WATER INTERFACE
  ! *** 
  ! *** USED FOR BOTH STANDARD EFDC SEDIMENT TRANSPORT AND SEDZLJ 
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2017-05-15        PAUL M. CRAIG    Updated to work for Original and SEDZLJ
  ! 2014-08-26        PAUL M. CRAIG    Updated to allow for any mix of ISTOC
  ! 2012-12-05        PAUL M. CRAIG    Updated OMP
  ! 2012-10-02        Dang H Chung     Added OMP
  ! *** *******************************************************************!

  use GLOBAL
  use Allocate_Initialize

  implicit none

  integer   :: NC, ND, LF, LL, LP, L, NT, K, KBTM1, NS
  integer   :: NFD, KM, KK, KBTP1, KBOT, KINC, KBOT2
  real      :: DEPINBED, DIFBWFAC, BETTMP, BSMALL
  real      :: HBEDMIN0, SORBMIN, TOXTIMEI, GAMTOX(KBM)
  real(RKD) :: CELLMASS, ERRT, ERRB, ERRW
  real,save :: TOXTIME
  
  real(RKD),save,allocatable,dimension(:)   :: DERRB
  real,save,allocatable,dimension(:,:,:)    :: PARTDIF     ! *** Sediment particle mixing and diffusion between layers, by toxic class
  real(RKD),save,allocatable,dimension(:,:) :: TOXBBALN    ! *** Total toxic mass in bed   after  advection/diffusion calculations
  real(RKD),save,allocatable,dimension(:,:) :: TOXBBALO    ! *** Total toxic mass in bed   before advection/diffusion calculations
  real(RKD),save,allocatable,dimension(:,:) :: TOXWBALN    ! *** Total toxic mass in water after  advection/diffusion calculations
  real(RKD),save,allocatable,dimension(:,:) :: TOXWBALO    ! *** Total toxic mass in water before advection/diffusion calculations
  real,save,allocatable,dimension(:,:,:)    :: BMNN        ! *** Main diagonal of tridiagonal matrix

  if( .not. allocated(DERRB) )then
    call AllocateDSI( DERRB,    KBM, 0.0)
    call AllocateDSI( PARTDIF,  LCM, KBM,  NTXM, 0.0)
    call AllocateDSI( TOXBBALN, LCM, NTXM, 0.0)
    call AllocateDSI( TOXBBALO, LCM, NTXM, 0.0)
    call AllocateDSI( TOXWBALN, LCM, NTXM, 0.0)
    call AllocateDSI( TOXWBALO, LCM, NTXM, 0.0)
    call AllocateDSI( BMNN,     LCM, -(KBM+1), NTXM, 0.0)
    
    ! *** Set up spatial surface diffusion/surface flux array (delme - spatially varying not implemented)
    do ND = 1,NDM
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED-1,LASED)
      do LP = LF,LL
        L = LSED(LP)
        DIFTOXBW(L,:) = DIFTOXS(:)    ! *** Implicit loop over NTOX
      enddo    
    enddo
  
    TOXTIME = 0.0
    TOXSTEPB = TOXSTEPB - DTSED/10.         ! *** Add a tolerance to prevent machine precision impacting update frequencies
  endif
  NFD = NSED2 + NSND + 1
  
  BSMALL = 1.0E-12
  TOXTIME = TOXTIME + DTSED
  
  ! *** ONLY COMPUTE FLUX TERMS AND BED MIXING IF TIME INTERVAL EXCEEDS TOXSTEPB
  if( TOXTIME < TOXSTEPB ) return         
  
  TOXTIMEI = 1./TOXTIME
  
  ! *** SET LAYER ORIENTATION
  if( LSEDZLJ )then
    KINC  = 1
    KBOT  = KB
    KBOT2 = KB-1
    HBEDMIN0 = BSMALL  ! *** ACCOUNT FOR ACTIVE LAYERS
  else
    KINC  = -1
    KBOT  = 1
    KBOT2 = 2
    HBEDMIN0 = BSMALL
  endif
  
  SORBMIN = 0.99
  
  ! *** *******************************************************************C
  !
  ! *** UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
  !
  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED-1,LASED)

    ! *** *******************************************************************C
    !
    ! *** CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
    !
    ! *** TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED 
    ! *** TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED
    call CALTOXB_FRACTIONS(LF,LL)
  enddo
  !$OMP END DO

  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED-1,LASED)
    ! *** COMPUTE MASS FRACTION (TOXPFTB) OF SORBED TOXICS TO PARTICULATES (DIMENSIONLESS)
    do NT = 1,NTOX
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          if( HBED(L,K) >= HBEDMIN0 .and. PORBED(L,K) > 0. )then
            ! ***                                                                             SORBED TO DOC
            TOXPFTB(L,K,NT) = ( TOXPFTB(L,K,NT)/(PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT)) ) - TOXPFB(L,K,NFD,NT)
            TOXPFTB(L,K,NT) = MAX(TOXPFTB(L,K,NT),SORBMIN)
          else
            TOXPFTB(L,K,NT) = SORBMIN
          endif
        enddo
      enddo
    enddo
  enddo   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, NC, NT, K, KM, KK, KBTP1, KBTM1)                     &
  !$OMP    PRIVATE(DEPINBED, DIFBWFAC, BETTMP, ERRT, CELLMASS, ERRB, ERRW, GAMTOX, DERRB)
  do ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED-1,LASED)

    ! *** *******************************************************************C
    !
    ! *** ADVECT AND DIFFUSE TOXICS IN BED AND INTO BOTTOM WATER
    ! *** COLUMN LAYER

    ! *** ADD PARTICLE MIXING AND SCALE PARTICLE DIFFUSION FOR SOLUTION
    do NT = 1,NTOX
      NC = 3 + NDYM + NT
      do LP = LF,LL
        L = LSED(LP)
          
        DEPINBED = 0.
        
        ! *** Particle Mixing
        PARTDIF(L,KBT(L),NT) = 0.0
        if( KBT(L) /= KBOT2 )then
          do K = KBT(L),KBOT2,KINC
            KM = K + KINC
            DEPINBED = DEPINBED + HBED(L,K)
            PARTDIF(L,KM,NT) = 0.0
            if( DEPINBED < DPDIFTOX(NT) )then
              PARTDIF(L,KM,NT) = 2.*PDIFTOX(NT)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
            endif
          enddo
        else
          K  = KBOT2
          KM = K + KINC
          DEPINBED = DEPINBED + HBED(L,K)
          PARTDIF(L,KM,NT) = 0.0
          if( DEPINBED < DPDIFTOX(NT) )then
            PARTDIF(L,KM,NT) = 2.*PDIFTOX(NT)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
          endif
        endif
      enddo
    
      do LP = LF,LL
        L = LSED(LP)
        if( HBED(L,KBT(L)) <= HBEDMIN0 ) CYCLE

        DIFBWFAC = 2./HBED(L,KBT(L))
        if( ISDIFBW(NT) == 1 ) DIFBWFAC = 1.0   ! *** Surface flux rate used
        
        TOXBBALO(L,NT) = 0.
        KBTP1 = KBT(L) - KINC   ! *** Layer above KBT
        KBTM1 = KBT(L) + KINC   ! *** Layer below KBT
        ALOW(L,KBOT,NT) = 0. 
        CUPP(L,KBTP1,NT) = 0.
        
        do K = KBOT,KBTM1,-KINC
          CUPP(L,K,NT) = MIN(QWTRBED(L,K),0.) - (DIFTOX(NT) + PARTDIF(L,K,NT))*(PORBED(L,K) + PORBED(L,K-KINC))/(HBED(L,K) + HBED(L,K-KINC))
        enddo
        CUPP(L,KBT(L),NT) = MIN(QWTRBED(L,KBT(L)),0.) - DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))       ! *** Surface diffusion/flux
        
        do K = KBOT2,KBT(L),-KINC
          ALOW(L,K,NT) = -MAX(QWTRBED(L,K+KINC),0.)   - (DIFTOX(NT) + PARTDIF(L,K+KINC,NT))*(PORBED(L,K+KINC) + PORBED(L,K))/(HBED(L,K + KINC) + HBED(L,K))
        enddo
        ALOW(L,KBTP1,NT) = -MAX(QWTRBED(L,KBT(L)),0.) - DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))       ! *** Surface diffusion/flux
        
        do K = KBOT,KBT(L),-KINC
          BMNN(L,K,NT) = TOXTIMEI*HBED(L,K)*PORBED(L,K)/(1.-TOXPFTB(L,K,NT))
        enddo
        BMNN(L,KBTP1,NT) = TOXTIMEI*HPK(L,KSZ(L))/(1.-TOXPFTW(L,KSZ(L),NT))                         ! *** Bottom layer of water column
        
        ! *** BOTTOM LAYER
        BMNN(L,KBOT,NT)  = BMNN(L,KBOT,NT) - MIN(QWTRBED(L,0),0.)
        BMNN(L,KBOT,NT)  = BMNN(L,KBOT,NT) + MAX(QWTRBED(L,KBOT),0.) + (DIFTOX(NT) + PARTDIF(L,KBOT,NT))*(PORBED(L,KBOT2) + PORBED(L,KBOT))/(HBED(L,KBOT2) + HBED(L,KBOT))
        
        ! *** MIDDLE LAYERS
        do K = KBOT2,KBTM1,-KINC
          BMNN(L,K,NT) = BMNN(L,K,NT) + MAX(QWTRBED(L,K),0.)      + (DIFTOX(NT) + PARTDIF(L,K,NT))     *(PORBED(L,K-KINC) + PORBED(L,K))/(HBED(L,K-KINC) + HBED(L,K)) &
                                      - MIN(QWTRBED(L,K+KINC),0.) + (DIFTOX(NT) + PARTDIF(L,K+KINC,NT))*(PORBED(L,K+KINC) + PORBED(L,K))/(HBED(L,K+KINC) + HBED(L,K))
        enddo
        
        ! *** Top layer of sediment column
        K = KBT(L)
        if( K == KBOT )then
          BMNN(L,K,NT) = BMNN(L,K,NT) + MAX(QWTRBED(L,K),0.)      + DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))
        else
          BMNN(L,K,NT) = BMNN(L,K,NT) + MAX(QWTRBED(L,K),0.)      + DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L)) &
                                      - MIN(QWTRBED(L,K+KINC),0.) + (DIFTOX(NT) + PARTDIF(L,K+KINC,NT))*(PORBED(L,K+KINC) + PORBED(L,K))/(HBED(L,K+KINC) + HBED(L,K))  
        endif
      
        ! *** ABOVE TOP LAYER
        BMNN(L,KBTP1,NT) = BMNN(L,KBTP1,NT) - MIN(QWTRBED(L,KBT(L)),0.) + DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))       ! *** Surface diffusion/flux

        ! *** COMPUTE RRHS KG/S
        do K = KBOT,KBT(L),-KINC
          RRHS(L,K,NT) =  TOXTIMEI*TOXB(L,K,NT)
          TOXBBALO(L,NT) = TOXBBALO(L,NT) + TOXB(L,K,NT)                            ! *** Total mass of toxics in bed before diffusion steps (mg/m2)
        enddo                                                       
        TOXWBALO(L,NT) = HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)                             ! *** Total mass of toxics in bottom layer before diffusion steps (mg/m2)
        
        ! *** ADD FLUX OF GROUNDWATER INTO THE BOTTOM SEDIMENT LAYER
        RRHS(L,KBOT,NT) = RRHS(L,KBOT,NT) + MAX(QWTRBED(L,0),0.)*CONGW(L,NC)
        RRHS(L,KBTP1,NT) = TOXTIMEI*TOXWBALO(L,NT)                                  ! *** RRHS - Mass flux rate (mg/m2/s) 
      enddo

      ! *** TRI-DIAGONAL SOLVER
      do LP = LF,LL
        L = LSED(LP)
        
        KBTP1 = KBT(L)-KINC           ! *** Layer above the top layer
        BETTMP = BMNN(L,KBOT,NT) + BSMALL
        TOXTMP(L,KBOT,NT) = RRHS(L,KBOT,NT)/BETTMP
        do KK = KBOT2,KBTP1,-KINC
          GAMTOX(KK) = CUPP(L,KK+KINC,NT)/BETTMP
          BETTMP = BMNN(L,KK,NT) - ALOW(L,KK,NT)*GAMTOX(KK) + BSMALL
          TOXTMP(L,KK,NT) = (RRHS(L,KK,NT) - ALOW(L,KK,NT)*TOXTMP(L,KK+KINC,NT))/BETTMP
        enddo
        do KK = KBT(L),KBOT,KINC
          TOXTMP(L,KK,NT) = TOXTMP(L,KK,NT) - GAMTOX(KK-KINC)*TOXTMP(L,KK-KINC,NT)
        enddo
      enddo

      ! *** CONVERT SCALED SOLUTION VARIABLES AND CALCULATE FINAL MASS
      do LP = LF,LL
        L = LSED(LP)

        TOXBBALN(L,NT) = 0.0
        do K = KBOT,KBT(L),-KINC
          TOXB(L,K,NT) = HBED(L,K)*PORBED(L,K)*TOXTMP(L,K,NT)/(1.0-TOXPFTB(L,K,NT))
          TOXBBALN(L,NT) = TOXBBALN(L,NT) + TOXB(L,K,NT)                            ! *** Total mass of toxics in bed after diffusion steps (mg/m2)
        enddo
        
        KBTP1 = KBT(L)-KINC
        TOX(L,KSZ(L),NT) = TOXTMP(L,KBTP1,NT)/(1.-TOXPFTW(L,KSZ(L),NT))
        TOXWBALN(L,NT) = HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)                             ! *** Total mass of toxics in bottom layer after diffusion steps (mg/m2)
      enddo
  
      ! *** Correct mass error and determine net flux from bed to water column
      do LP = LF,LL
        L = LSED(LP)
        
        CELLMASS = TOXBBALN(L,NT) + TOXWBALN(L,NT)
        ERRT     = CELLMASS - (TOXBBALO(L,NT) + TOXWBALO(L,NT))                     ! *** Difference in bed mass at the start and end of calculations
        ERRT     = ERRT - MAX(QWTRBED(L,0),0.)*CONGW(L,NC)*TOXTIME                  ! *** Remove the external mass loading
        
        ! *** HANDLE ZERO SEDIMENT AND/OR WATER CONCENTRATIONS
        if( TOXBBALN(L,NT) > 1.E-12 .and. ERRT /= 0. )then
          ERRB = ERRT*TOXBBALN(L,NT)/CELLMASS
          ERRW = ERRT - ERRB
          TOX(L,KSZ(L),NT) = TOX(L,KSZ(L),NT) - ERRW/HPK(L,KSZ(L))
          do K = KBOT,KBT(L),-KINC
            DERRB(K) = TOXB(L,K,NT)/TOXBBALN(L,NT)
          enddo
          do K = KBOT,KBT(L),-KINC
            TOXB(L,K,NT) = TOXB(L,K,NT) - DERRB(K)*ERRB
          enddo
          TADFLUX(L,NT) = (HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)-TOXWBALO(L,NT))/TOXTIME
        else
          TADFLUX(L,NT) = 0.
        endif


      enddo
      
    enddo ! *** End of contaminant class
  enddo   ! *** End of domain loop
  !$OMP END DO
  !$OMP END PARALLEL

  TOXTIME = 0.0
  
  return

END

SUBROUTINE CALTOXB_FRACTIONS(LF,LL)

  ! *** *******************************************************************C
  !
  ! ***  UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
  !

  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2014-08-26        PAUL M. CRAIG    Updated to allow for any mix of ISTOC

  ! ***  TOXPFB(L,NS,NT)  = PARTICULATE FRACTION IN SEDIMENT BED (DIMENSIONLESS, INTERMEDIATE UNITS METERS) 
  ! ***  TOXPARB(L,NS,NT) = PARTITION COEFFICIENT IN SEDIMENT BED (L/MG = M3/G)
  ! ***  TOXPFTB(L,NT)    = TOTAL PARTICULATE FRACTION IN SEDIMENT BED (DIMENSIONLESS)
  ! ***  STFPOCB(L,K,NS)  = FRACTION OF OC ON SEDIMENT CLASS NS (DIMENSIONLESS)
  ! ***  SEDB(L,K)        = COHESIVE SEDIMENTS (G/M2)
  ! ***  SNDB(L,K)        = NON-COHESIVE SEDIMENTS (G/M2)
  ! ***  STDOCB(L,K)      = DOC CONCENTRATION IN SEDIMENTS (G/M3)
  ! ***  STPOCB(L,K)      = POC CONCENTRATION IN SEDIMENTS, NON SEDIMENT COMPONENT (E.G. ALGAE) (G/M3)

  use GLOBAL
  
  implicit none

  integer, intent(IN) :: LL, LF
  integer             :: LP, L, K, NT, NS, NX, NFD

  real                :: HBEDMIN0, TMPVAL,  PMC(10)

  if( LSEDZLJ )then
      HBEDMIN0 = 1E-12
  else
      HBEDMIN0 = 1E-12
  endif
  NFD = NSED2 + NSND + 1

  ! *** ZERO THE TOXIC SORBED FRACTIONS
  do NT = 1,NTOX
    do NS = 1,NSP2(NT)
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          TOXPFB(L,K,NS,NT) = 0.
        enddo
      enddo
    enddo
  enddo

  do NT = 1,NTOX

    ! *** PARTITION TO COHESIVES
    if( ISTRAN(6) >= 1 )then
      if( ISTOC(NT) > 1 )then
        ! *** fPOC BASED
        do NS = 1,NSED2
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2                       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SEDB(L,K,NS)*STFPOCB(L,K,NS)*TOXPARB(L,NS,NT)
            enddo
          enddo
        enddo
      elseif( ISTOC(NT) == 0 )then
        ! *** Kd APPROACH
        do NS = 1,NSED2
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SEDB(L,K,NS)*TOXPARB(L,NS,NT)
            enddo
          enddo
        enddo
      endif
    endif
    
    ! *** PARTITION TO NONCOHESIVES
    if( ISTRAN(7) >= 1 )then
      do NX = 1,NSND
        NS = NX + NSED2
        if( ISTOC(NT) > 1 )then
          ! *** fPOC BASED
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2                       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SNDB(L,K,NX)*STFPOCB(L,K,NS)*TOXPARB(L,NS,NT)
            enddo
          enddo
        elseif( ISTOC(NT) == 0 )then 
          ! *** Kd APPROACH
          do K = 1,KB
            do LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2      L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SNDB(L,K,NX)*TOXPARB(L,NS,NT)
            enddo
          enddo
        endif
      enddo
    endif
    
    ! *** PARTITION (COMPLEX) TO DOC
    if( ISTOC(NT) == 1 .or. ISTOC(NT) == 2 )then
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          ! ***   M                      M           DOC(G/M3)    L/MG (M3/G)
          TOXPFB(L,K,NFD,NT) = PORBED(L,K)*HBED(L,K)*STDOCB(L,K)*TOXPARBC(1,NT)
        enddo
      enddo
    endif
      
    ! *** POC SORPTION (NON-SEDIMENT RELATED) COMPONENT
    if( ISTOC(NT) == 1 )then
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          ! ***    M        =   M        POC(G/M3)    L/MG (M3/G)
          TOXPFB(L,K,NFD+1,NT) = HBED(L,K)*STPOCB(L,K)*TOXPARBC(2,NT)
        enddo
      enddo
    endif
  enddo
  
  ! *** TOXPFTB (M) IS TEMPORARILY USED TO STORE TOTAL SORBED PER TOXIC
  do NT = 1,NTOX
    do K = 1,KB
      do LP = LF,LL
        L = LSED(LP)
        TOXPFTB(L,K,NT) = 0.
      enddo
    enddo
    do NS = 1,NSP2(NT)
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          ! *** M
          TOXPFTB(L,K,NT) = TOXPFTB(L,K,NT) + TOXPFB(L,K,NS,NT)
        enddo
      enddo
    enddo
  enddo
  
  ! *** COMPUTE MASS FRACTION (TOXPFB) OF SORBED TOXICS FOR EACH SEDIMENT CLASS (DIMENSIONLESS)
  do NT = 1,NTOX
    do NS = 1,NSP2(NT)
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          TMPVAL = PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT)
          if( TMPVAL > 0. )then
            TOXPFB(L,K,NS,NT) = TOXPFB(L,K,NS,NT)/TMPVAL
          else
            TOXPFB(L,K,NS,NT) = 1.
          endif
        enddo
      enddo
    enddo
  enddo
  
  if( ISTMSR >= 1 .or. ISFDCH > 0 )then
    ! *** ONLY NEEDED FOR TIME SERIES USING TMSR OR FOODCHAIN
    do NT = 1,NTOX
      do K = 1,KB
        do LP = LF,LL
          L = LSED(LP)
          if( HBED(L,K) > HBEDMIN0 )then
            TOXFDFB(L,K,NT) = PORBED(L,K)*HBED(L,K) /(PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT))
            TOXCDFB(L,K,NT) = TOXPFB(L,K,NFD,NT)
          else
            TOXFDFB(L,K,NT) = 0.
            TOXCDFB(L,K,NT) = 0.
          endif
        enddo
      enddo
    enddo
  endif

END 
