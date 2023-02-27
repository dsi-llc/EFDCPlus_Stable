! ----------------------------------------------------------------------
!   This file is a part of EFDC + 
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

  USE GLOBAL
  Use Allocate_Initialize

  IMPLICIT NONE

  INTEGER   :: NC, ND, LF, LL, LP, L, NT, K, KBTM1, NS
  INTEGER   :: NFD, KM, KK, KBTP1, KBOT, KINC, KBOT2
  REAL      :: DEPINBED, DIFBWFAC, BETTMP
  REAL      :: HBEDMIN0, SORBMIN, TOXTIMEI, GAMTOX(KBM)
  REAL(RKD) :: CELLMASS, ERRT, ERRB, ERRW
  REAL,SAVE :: TOXTIME
  
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:)   :: DERRB
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:)    :: PARTDIF     ! *** Sediment particle mixing and diffusion between layers, by toxic class
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXBBALN    ! *** Total toxic mass in bed   after  advection/diffusion calculations
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXBBALO    ! *** Total toxic mass in bed   before advection/diffusion calculations
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXWBALN    ! *** Total toxic mass in water after  advection/diffusion calculations
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: TOXWBALO    ! *** Total toxic mass in water before advection/diffusion calculations

  IF( .NOT. ALLOCATED(DERRB) )THEN
    Call AllocateDSI( DERRB,    KBM, 0.0)
    Call AllocateDSI( PARTDIF,  LCM, KBM,  NTXM, 0.0)
    Call AllocateDSI( TOXBBALN, LCM, NTXM, 0.0)
    Call AllocateDSI( TOXBBALO, LCM, NTXM, 0.0)
    Call AllocateDSI( TOXWBALN, LCM, NTXM, 0.0)
    Call AllocateDSI( TOXWBALO, LCM, NTXM, 0.0)
    
    ! *** Set up spatial surface diffusion/surface flux array (delme - spatially varying not implemented)
    DO ND = 1,NDM
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED-1,LASED)
      DO LP = LF,LL
        L = LSED(LP)
        DIFTOXBW(L,:) = DIFTOXS(:)    ! *** Implicit loop over NTOX
      ENDDO    
    ENDDO
  
    TOXTIME = 0.0
    TOXSTEPB = TOXSTEPB - DTSED/10.         ! *** Add a tolerance to prevent machine precision impacting update frequencies
  ENDIF
  NFD = NSED2 + NSND + 1
  
  TOXTIME = TOXTIME + DTSED
  
  ! *** ONLY COMPUTE FLUX TERMS AND BED MIXING IF TIME INTERVAL EXCEEDS TOXSTEPB
  IF( TOXTIME < TOXSTEPB ) RETURN         
  
  TOXTIMEI = 1./TOXTIME
  
  ! *** SET LAYER ORIENTATION
  IF( LSEDZLJ )THEN
    KINC  = 1
    KBOT  = KB
    KBOT2 = KB-1
    HBEDMIN0 = 1E-12  ! *** ACCOUNT FOR ACTIVE LAYERS
  ELSE
    KINC  = -1
    KBOT  = 1
    KBOT2 = 2
    HBEDMIN0 = 1E-12
  ENDIF
  
  SORBMIN = 0.99
  
  ! *** *******************************************************************C
  !
  ! *** UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
  !
  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT)
  DO ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED-1,LASED)

    ! *** *******************************************************************C
    !
    ! *** CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
    !
    ! *** TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED 
    ! *** TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED
    CALL CALTOXB_FRACTIONS(LF,LL)
  ENDDO
  !$OMP END DO

  !$OMP DO PRIVATE(ND, LF, LL, LP, L, K, NT)
  DO ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED-1,LASED)
    ! *** COMPUTE MASS FRACTION (TOXPFTB) OF SORBED TOXICS TO PARTICULATES (DIMENSIONLESS)
    DO NT = 1,NTOX
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          IF( HBED(L,K) >= HBEDMIN0 .AND. PORBED(L,K) > 0. )THEN
            ! ***                                                                             SORBED TO DOC
            TOXPFTB(L,K,NT) = ( TOXPFTB(L,K,NT)/(PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT)) ) - TOXPFB(L,K,NFD,NT)
            TOXPFTB(L,K,NT) = MAX(TOXPFTB(L,K,NT),SORBMIN)
          ELSE
            TOXPFTB(L,K,NT) = SORBMIN
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  !$OMP DO PRIVATE(ND, LF, LL, LP, L, NC, NT, K, KM, KK, KBTP1, KBTM1)                     &
  !$OMP    PRIVATE(DEPINBED, DIFBWFAC, BETTMP, ERRT, CELLMASS, ERRB, ERRW, GAMTOX, DERRB)
  DO ND = 1,NDM
    LF = (ND-1)*LDMSED + 1  
    LL = MIN(LF + LDMSED-1,LASED)

    ! *** *******************************************************************C
    !
    ! *** ADVECT AND DIFFUSE TOXICS IN BED AND INTO BOTTOM WATER
    ! *** COLUMN LAYER

    ! *** ADD PARTICLE MIXING AND SCALE PARTICLE DIFFUSION FOR SOLUTION
    DO NT = 1,NTOX
      NC = 3 + NDYM + NT
      DO LP = LF,LL
        L = LSED(LP)
          
        DEPINBED = 0.
        
        ! *** Particle Mixing
        PARTDIF(L,KBT(L),NT) = 0.0
        IF( KBT(L) /= KBOT2 )THEN
          DO K = KBT(L),KBOT2,KINC
            KM = K + KINC
            DEPINBED = DEPINBED + HBED(L,K)
            PARTDIF(L,KM,NT) = 0.0
            IF( DEPINBED < DPDIFTOX(NT) )THEN
              PARTDIF(L,KM,NT) = 2.*PDIFTOX(NT)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
            ENDIF
          ENDDO
        ELSE
          K  = KBOT2
          KM = K + KINC
          DEPINBED = DEPINBED + HBED(L,K)
          PARTDIF(L,KM,NT) = 0.0
          IF( DEPINBED < DPDIFTOX(NT) )THEN
            PARTDIF(L,KM,NT) = 2.*PDIFTOX(NT)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
          ENDIF
        ENDIF
      ENDDO
    
      DO LP = LF,LL
        L = LSED(LP)
        IF( HBED(L,KBT(L)) <= HBEDMIN0 ) CYCLE

        DIFBWFAC = 2./HBED(L,KBT(L))
        IF( ISDIFBW(NT) == 1 ) DIFBWFAC = 1.0   ! *** Surface flux rate used
        
        TOXBBALO(L,NT) = 0.
        KBTP1 = KBT(L) - KINC   ! *** Layer above KBT
        KBTM1 = KBT(L) + KINC   ! *** Layer below KBT
        ALOW(L,KBOT,NT) = 0. 
        CUPP(L,KBTP1,NT) = 0.
        
        DO K = KBOT,KBTM1,-KINC
          CUPP(L,K,NT) = MIN(QWTRBED(L,K),0.) - (DIFTOX(NT) + PARTDIF(L,K,NT))*(PORBED(L,K) + PORBED(L,K-KINC))/(HBED(L,K) + HBED(L,K-KINC))
        ENDDO
        CUPP(L,KBT(L),NT) = MIN(QWTRBED(L,KBT(L)),0.) - DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))       ! *** Surface diffusion/flux
        
        DO K = KBOT2,KBT(L),-KINC
          ALOW(L,K,NT) = -MAX(QWTRBED(L,K+KINC),0.)   - (DIFTOX(NT) + PARTDIF(L,K+KINC,NT))*(PORBED(L,K+KINC) + PORBED(L,K))/(HBED(L,K + KINC) + HBED(L,K))
        ENDDO
        ALOW(L,KBTP1,NT) = -MAX(QWTRBED(L,KBT(L)),0.) - DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))       ! *** Surface diffusion/flux
        
        DO K = KBOT,KBT(L),-KINC
          BMNN(L,K,NT) = TOXTIMEI*HBED(L,K)*PORBED(L,K)/(1.-TOXPFTB(L,K,NT))
        ENDDO
        BMNN(L,KBTP1,NT) = TOXTIMEI*HPK(L,KSZ(L))/(1.-TOXPFTW(L,KSZ(L),NT))                         ! *** Bottom layer of water column
        
        ! *** BOTTOM LAYER
        BMNN(L,KBOT,NT)  = BMNN(L,KBOT,NT) - MIN(QWTRBED(L,0),0.)
        BMNN(L,KBOT,NT)  = BMNN(L,KBOT,NT) + MAX(QWTRBED(L,KBOT),0.) + (DIFTOX(NT) + PARTDIF(L,KBOT,NT))*(PORBED(L,KBOT2) + PORBED(L,KBOT))/(HBED(L,KBOT2) + HBED(L,KBOT))
        
        ! *** MIDDLE LAYERS
        DO K = KBOT2,KBTM1,-KINC
          BMNN(L,K,NT) = BMNN(L,K,NT) + MAX(QWTRBED(L,K),0.)      + (DIFTOX(NT) + PARTDIF(L,K,NT))     *(PORBED(L,K-KINC) + PORBED(L,K))/(HBED(L,K-KINC) + HBED(L,K)) &
                                      - MIN(QWTRBED(L,K+KINC),0.) + (DIFTOX(NT) + PARTDIF(L,K+KINC,NT))*(PORBED(L,K+KINC) + PORBED(L,K))/(HBED(L,K+KINC) + HBED(L,K))
        ENDDO
        
        ! *** Top layer of sediment column
        K = KBT(L)
        IF( K == KBOT )THEN
          BMNN(L,K,NT) = BMNN(L,K,NT) + MAX(QWTRBED(L,K),0.)      + DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))
        ELSE
          BMNN(L,K,NT) = BMNN(L,K,NT) + MAX(QWTRBED(L,K),0.)      + DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L)) &
                                      - MIN(QWTRBED(L,K+KINC),0.) + (DIFTOX(NT) + PARTDIF(L,K+KINC,NT))*(PORBED(L,K+KINC) + PORBED(L,K))/(HBED(L,K+KINC) + HBED(L,K))  
        ENDIF
      
        ! *** ABOVE TOP LAYER
        BMNN(L,KBTP1,NT) = BMNN(L,KBTP1,NT) - MIN(QWTRBED(L,KBT(L)),0.) + DIFBWFAC*DIFTOXBW(L,NT)*PORBED(L,KBT(L))       ! *** Surface diffusion/flux

        ! *** COMPUTE RRHS KG/S
        DO K = KBOT,KBT(L),-KINC
          RRHS(L,K,NT) =  TOXTIMEI*TOXB(L,K,NT)
          TOXBBALO(L,NT) = TOXBBALO(L,NT) + TOXB(L,K,NT)                  ! *** Total mass of toxics in bed before diffusion steps (mg/m2)
        ENDDO                                                       
        TOXWBALO(L,NT) = HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)                ! *** Total mass of toxics in bottom layer before diffusion steps (mg/m2)
        
        ! *** ADD FLUX OF GROUNDWATER INTO THE BOTTOM SEDIMENT LAYER
        RRHS(L,KBOT,NT) = RRHS(L,KBOT,NT) + MAX(QWTRBED(L,0),0.)*CONGW(L,NC)
        RRHS(L,KBTP1,NT) = TOXTIMEI*TOXWBALO(L,NT)                     ! *** RRHS - Mass flux rate (mg/m2/s) 
      ENDDO

      ! *** TRI-DIAGONAL SOLVER
      DO LP = LF,LL
        L = LSED(LP)
        
        KBTP1 = KBT(L)-KINC           ! *** Layer above the top layer
        BETTMP = BMNN(L,KBOT,NT)
        TOXTMP(L,KBOT,NT) = RRHS(L,KBOT,NT)/BETTMP
        DO KK = KBOT2,KBTP1,-KINC
          GAMTOX(KK) = CUPP(L,KK+KINC,NT)/BETTMP
          BETTMP = BMNN(L,KK,NT) - ALOW(L,KK,NT)*GAMTOX(KK) + 1E-12
          TOXTMP(L,KK,NT) = (RRHS(L,KK,NT) - ALOW(L,KK,NT)*TOXTMP(L,KK+KINC,NT))/BETTMP
        ENDDO
        DO KK = KBT(L),KBOT,KINC
          TOXTMP(L,KK,NT) = TOXTMP(L,KK,NT) - GAMTOX(KK-KINC)*TOXTMP(L,KK-KINC,NT)
        ENDDO
      ENDDO

      ! *** CONVERT SCALED SOLUTION VARIABLES AND CALCULATE FINAL MASS
      DO LP = LF,LL
        L = LSED(LP)

        TOXBBALN(L,NT) = 0.0
        DO K = KBOT,KBT(L),-KINC
          TOXB(L,K,NT) = HBED(L,K)*PORBED(L,K)*TOXTMP(L,K,NT)/(1.0-TOXPFTB(L,K,NT))
          TOXBBALN(L,NT) = TOXBBALN(L,NT) + TOXB(L,K,NT)                                 ! *** Total mass of toxics in bed after diffusion steps (mg/m2)
        ENDDO
        
        KBTP1 = KBT(L)-KINC
        TOX(L,KSZ(L),NT) = TOXTMP(L,KBTP1,NT)/(1.-TOXPFTW(L,KSZ(L),NT))
        TOXWBALN(L,NT) = HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)                               ! *** Total mass of toxics in bottom layer after diffusion steps (mg/m2)
      ENDDO
  
      ! *** CORRECT MASS ERROR AND DETERMINE NET FLUX FROM BED TO WATER COLUMN
      DO LP = LF,LL
        L = LSED(LP)
        
        CELLMASS = TOXBBALN(L,NT) + TOXWBALN(L,NT)
        ERRT     = CELLMASS - (TOXBBALO(L,NT) + TOXWBALO(L,NT))                          ! *** Difference in bed mass at the start and end of calculations
        ERRT     = ERRT - MAX(QWTRBED(L,0),0.)*CONGW(L,NC)*TOXTIME                 ! *** Remove the external mass loading
        
        ! *** HANDLE ZERO SEDIMENT AND/OR WATER CONCENTRATIONS
        IF( TOXBBALN(L,NT) > 1.E-12 .AND. ERRT /= 0. )THEN
          ERRB = ERRT*TOXBBALN(L,NT)/CELLMASS
          ERRW = ERRT - ERRB
          TOX(L,KSZ(L),NT) = TOX(L,KSZ(L),NT) - ERRW/HPK(L,KSZ(L))
          DO K = KBOT,KBT(L),-KINC
            DERRB(K) = TOXB(L,K,NT)/TOXBBALN(L,NT)
          ENDDO
          DO K = KBOT,KBT(L),-KINC
            TOXB(L,K,NT) = TOXB(L,K,NT) - DERRB(K)*ERRB
          ENDDO
          TADFLUX(L,NT) = (HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)-TOXWBALO(L,NT))/TOXTIME
        ELSE
          TADFLUX(L,NT) = 0.
        ENDIF


      ENDDO
      
    ENDDO ! *** End of contaminant class
  ENDDO   ! *** End of domain loop
  !$OMP END DO
  !$OMP END PARALLEL

  TOXTIME = 0.0
  
  RETURN

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

  USE GLOBAL
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: LL, LF
  INTEGER             :: LP, L, K, NT, NS, NX, NFD

  REAL                :: HBEDMIN0, TMPVAL,  PMC(10)

  IF( LSEDZLJ )THEN
      HBEDMIN0 = 1E-12
  ELSE
      HBEDMIN0 = 1E-12
  ENDIF
  NFD = NSED2 + NSND + 1

  ! *** ZERO THE TOXIC SORBED FRACTIONS
  DO NT = 1,NTOX
    DO NS = 1,NSP2(NT)
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          TOXPFB(L,K,NS,NT) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO NT = 1,NTOX

    ! *** PARTITION TO COHESIVES
    IF( ISTRAN(6) >= 1 )THEN
      IF( ISTOC(NT) > 1 )THEN
        ! *** fPOC BASED
        DO NS = 1,NSED2
          DO K = 1,KB
            DO LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2                       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SEDB(L,K,NS)*STFPOCB(L,K,NS)*TOXPARB(L,NS,NT)
            ENDDO
          ENDDO
        ENDDO
      ELSEIF( ISTOC(NT) == 0 )THEN
        ! *** Kd APPROACH
        DO NS = 1,NSED2
          DO K = 1,KB
            DO LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SEDB(L,K,NS)*TOXPARB(L,NS,NT)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    ! *** PARTITION TO NONCOHESIVES
    IF( ISTRAN(7) >= 1 )THEN
      DO NX = 1,NSND
        NS = NX + NSED2
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO K = 1,KB
            DO LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2                       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SNDB(L,K,NX)*STFPOCB(L,K,NS)*TOXPARB(L,NS,NT)
            ENDDO
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN 
          ! *** Kd APPROACH
          DO K = 1,KB
            DO LP = LF,LL
              L = LSED(LP)
              ! ***     M       =    G/M2      L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SNDB(L,K,NX)*TOXPARB(L,NS,NT)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    
    ! *** PARTITION (COMPLEX) TO DOC
    IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2 )THEN
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          ! ***   M                      M           DOC(G/M3)    L/MG (M3/G)
          TOXPFB(L,K,NFD,NT) = PORBED(L,K)*HBED(L,K)*STDOCB(L,K)*TOXPARBC(1,NT)
        ENDDO
      ENDDO
    ENDIF
      
    ! *** POC SORPTION (NON-SEDIMENT RELATED) COMPONENT
    IF( ISTOC(NT) == 1 )THEN
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          ! ***    M        =   M        POC(G/M3)    L/MG (M3/G)
          TOXPFB(L,K,NFD+1,NT) = HBED(L,K)*STPOCB(L,K)*TOXPARBC(2,NT)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  
  ! ** TOXPFTB (M) IS TEMPORARILY USED TO STORE TOTAL SORBED PER TOXIC
  DO NT = 1,NTOX
    DO K = 1,KB
      DO LP = LF,LL
        L = LSED(LP)
        TOXPFTB(L,K,NT) = 0.
      ENDDO
    ENDDO
    DO NS = 1,NSP2(NT)
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          ! *** M
          TOXPFTB(L,K,NT) = TOXPFTB(L,K,NT) + TOXPFB(L,K,NS,NT)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  
  ! *** COMPUTE MASS FRACTION (TOXPFB) OF SORBED TOXICS FOR EACH SEDIMENT CLASS (DIMENSIONLESS)
  DO NT = 1,NTOX
    DO NS = 1,NSP2(NT)
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          TMPVAL = PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT)
          IF( TMPVAL > 0. )THEN
            TOXPFB(L,K,NS,NT) = TOXPFB(L,K,NS,NT)/TMPVAL
          ELSE
            TOXPFB(L,K,NS,NT) = 1.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  
  IF( ISTMSR >= 1 .OR. ISFDCH > 0 )THEN
    ! *** ONLY NEEDED FOR TIME SERIES USING TMSR OR FOODCHAIN
    DO NT = 1,NTOX
      DO K = 1,KB
        DO LP = LF,LL
          L = LSED(LP)
          IF( HBED(L,K) > HBEDMIN0 )THEN
            TOXFDFB(L,K,NT) = PORBED(L,K)*HBED(L,K) /(PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT))
            TOXCDFB(L,K,NT) = TOXPFB(L,K,NFD,NT)
          ELSE
            TOXFDFB(L,K,NT) = 0.
            TOXCDFB(L,K,NT) = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

END 
