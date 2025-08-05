! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALAVB ()  

  ! ****************************************************************************
  ! *** SUBROUTINE CALAVB CALCULATES VERTICAL VISCOSITY AND DIFFUSIVITY  
  ! *** USING GLAPERIN ET AL'S MODIFICATION OF THE MELLOR-YAMADA MODEL  
  ! *** (NOTE AV, AB, AND AQ ARE ACTUALLY DIVIDED BY H)  
  ! *** IF ISFAVB = 0 NO TIME FILTER
  ! *** IF ISFAVB = 1 VALUES ARE ARITHMETIC AVERAGES WITH THE PREVIOUS VALUE
  ! *** IF ISFAVB = 2 VALUES ARE GEOMETRIC AVERAGES WITH THE PREVIOUS VALUE
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2018-10       PAUL M. CRAIG     IMPLEMENTED SPATIALLY VARYING AVO AND ABO
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE8.0
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2011-03-11    PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP
  !                  JOHN HAMRICK      ADDED DRYCELL BYPASS AND CONSISTENT INITIALIZATION OF DRY VALUES  

  use GLOBAL  
  use Allocate_Initialize      

  implicit none

  integer :: L, K, ND, LF, LL, LP
  real    :: QQIMAX, RIQMIN, RIQ, SFAV, SFAB, ABTMP, AVTMP
  real    :: SFAV0, SFAV1, SFAV2, SFAV3, SFAB0, SFAB1, RIQMAXX

  real,save,allocatable,dimension(:) :: AVOXYHPI

  ! *** FIRST CALL ALLOCATIONS AND ASSIGNMENTS
  if( .not. allocated(AVOXYHPI) )then
    call AllocateDSI(AVOXYHPI, LCM, 0.0)
  endif
  
  ! ***  GALPERIN (DEFAULT)
  SFAV0= 0.392010
  SFAV1= 7.760050
  SFAV2 = 34.676440
  SFAV3= 6.127200
  SFAB0= 0.493928
  SFAB1 = 34.676440
  RIQMIN = -0.999/SFAB1

  ! *** KANTHA AND CLAYSON (1994)
  if( ISTOPT(0) == 2 )then
    SFAV0= 0.392010
    SFAV1= 8.679790
    SFAV2 = 30.192000
    SFAV3= 6.127200
    SFAB0= 0.493928
    SFAB1 = 30.192000
    RIQMIN = -0.999/SFAB1
  endif

  ! ***  KANTHA (2003)
  if( ISTOPT(0) == 3 )then
    SFAV0= 0.392010
    SFAV1 = 14.509100
    SFAV2 = 24.388300
    SFAV3= 3.236400
    SFAB0= 0.490025
    SFAB1 = 24.388300
    RIQMIN = -0.999/SFAB1
  endif

  RIQMAXX = 1.E32
  if( ISLLIM >= 1 ) RIQMAXX = RIQMAX
  
  QQIMAX = 1./QQMIN  
  ! ****************************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  if( ISAVCOMP == 0 )then
    ! *** CONSTANT AV
    !$OMP DO PRIVATE(ND,K,LP,L) 
    do ND = 1,NDM
      ! *** VERTICAL DIFFUSVITY TIME FILTERING - NONE
      do K = 1,KS  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          AV(L,K) = AVOXY(L)*HPI(L) 
          AB(L,K) = AVBXY(L)*HPI(L)
        enddo  
      enddo  
    enddo
    !$OMP END DO
  
  elseif( BSC <= 1.E-6 )then
    ! *** VERTICAL DIFFUSVITY TIME FILTERING - NONE
    !$OMP DO PRIVATE(ND,K,LP,L,ABTMP,AVTMP) 
    do ND = 1,NDM
      if( ISFAVB == 0 .or. N == 1 )then
        ! *** VERTICAL DIFFUSVITY TIME FILTERING - NONE
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            AB(L,K) = SFAB0*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AV(L,K) = SFAV0*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = AV(L,K)*HPI(L) 
            AB(L,K) = AB(L,K)*HPI(L)
          enddo  
        enddo  

      elseif( ISFAVB == 1 .and. N > 1 )then
        ! *** VERTICAL DIFFUSVITY TIME FILTERING - AVERAGE
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ABTMP = SFAB0*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = SFAV0*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = 0.5*(AV(L,K) + AVTMP*HPI(L))  
            AB(L,K) = 0.5*(AB(L,K) + ABTMP*HPI(L))  
          enddo  
        enddo  

      elseif( ISFAVB == 2 .and. N > 1 )then
        ! *** VERTICAL DIFFUSVITY TIME FILTERING - SQRT
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ABTMP = SFAB0*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = SFAV0*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = SQRT(AV(L,K)*AVTMP*HPI(L))  
            AB(L,K) = SQRT(AB(L,K)*ABTMP*HPI(L))  
          enddo  
        enddo  
      endif
    enddo
    !$OMP END DO
  
  else
    
    !----------------------------------------------------------------------C
    if( ISFAVB == 0 .or. N == 1 )then
      ! *** VERTICAL DIFFUSIVITY TIME FILTERING - NONE
      !$OMP DO PRIVATE(ND,K,LP,L,RIQ,SFAV,SFAB) 
      do ND = 1,NDM  
    
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            RIQ = -GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(L,K)*(B(L,K+1)-B(L,K))/QQ(L,K)  
            RIQ = MAX(RIQ,RIQMIN)  
            RIQ = MIN(RIQ,RIQMAXX)
            SFAV = SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB = SFAB0/(1.+SFAB1*RIQ)
            AB(L,K) = SFAB*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)
            AV(L,K) = SFAV*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)
            AV(L,K) = AV(L,K)*HPI(L)  
            AB(L,K) = AB(L,K)*HPI(L)
          enddo  
        enddo  
      enddo
      !$OMP END DO

    !----------------------------------------------------------------------C
    elseif( ISFAVB == 1 .and. N > 1 )then
      ! *** VERTICAL DIFFUSIVITY TIME FILTERING - AVERAGE
      !$OMP DO PRIVATE(ND,K,LP,L,RIQ,SFAV,SFAB,AVTMP,ABTMP) 
      do ND = 1,NDM  
    
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            RIQ = -GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(L,K)*(B(L,K+1)-B(L,K))/QQ(L,K)  
            RIQ = MAX(RIQ,RIQMIN)  
            RIQ = MIN(RIQ,RIQMAXX)  
            SFAV = SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB = SFAB0/(1.+SFAB1*RIQ)
            ABTMP = SFAB*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = SFAV*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = 0.5*(AV(L,K) + AVTMP*HPI(L))  
            AB(L,K) = 0.5*(AB(L,K) + ABTMP*HPI(L))  
          enddo  
        enddo  
      enddo
      !$OMP END DO

    !----------------------------------------------------------------------C
    elseif( ISFAVB == 2 .and. N > 1 )then
      ! *** VERTICAL DIFFUSIVITY TIME FILTERING - SQRT
      !$OMP DO PRIVATE(ND,K,LP,L,RIQ,SFAV,SFAB,AVTMP,ABTMP) 
      do ND = 1,NDM  
    
        do K = 1,KS  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            RIQ = -GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(L,K)*(B(L,K+1)-B(L,K))/QQ(L,K)  
            RIQ = MAX(RIQ,RIQMIN)  
            RIQ = MIN(RIQ,RIQMAXX)  
            SFAV = SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB = SFAB0/(1.+SFAB1*RIQ)
            ABTMP = SFAB*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = SFAV*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = SQRT(AV(L,K)*AVTMP*HPI(L))  
            AB(L,K) = SQRT(AB(L,K)*ABTMP*HPI(L))  
          enddo  
        enddo  
      enddo
      !$OMP END DO
    endif
    
  endif
  
  ! *** APPLY OPEN BOUNDARYS 
  !$OMP SINGLE
  do LL = 1,NBCSOP
    L = LOBCS(LL)
    do K = 1,KS  
      AB(L,K) = 0.0
    enddo  
  enddo 
  !$OMP END SINGLE

  ! ****************************************************************************
  if( ISLOG >= 1 )then  
    !$OMP SINGLE
    AVMAX = AVO  
    ABMAX = ABO  
    AVMIN = 10.  
    ABMIN = 10.  
    do K = 1,KS  
      do L = 2,LA  
        if( LKSZ(L,K) ) CYCLE 
        AVMAX = MAX(AVMAX,AV(L,K))  
        ABMAX = MAX(ABMAX,AB(L,K))  
        AVMIN = MIN(AVMIN,AV(L,K))  
        ABMIN = MIN(ABMIN,AB(L,K))  
      enddo
    enddo
    !$OMP END SINGLE
  endif

  ! ****************************************************************************
  ! *** CHECK FOR DEPTHS LESS THAN ZBR
  if( ISDRY > 0 )then
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
    
      do LP = LF,LL  
        L = LWET(LP)  
        if( HPK(L,KSZ(L)) < ZBR(L) )then
          do K = KS-1,KSZ(L),-1
            if( HP(L)*Z(L,K-1) > ZBR(L) ) CYCLE
            AV(L,K) = AV(L,K+1)
            AB(L,K) = AB(L,K+1)
          enddo
        endif
      enddo
    enddo
    !$OMP END DO
  endif
        
  ! ****************************************************************************
  ! *** NOW APPLY MAXIMUM, IF REQURIED
  if( ISAVBMX >= 1 )then  
    !$OMP DO PRIVATE(ND,K,LP,L,AVTMP,ABTMP)
    do ND = 1,NDM  
    
      do K = 1,KS  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          AVTMP = AVMX*HPI(L)  
          ABTMP = ABMX*HPI(L)  
          AV(L,K) = MIN(AV(L,K),AVTMP)  
          AB(L,K) = MIN(AB(L,K),ABTMP)  
        enddo  
      enddo  
    enddo
    !$OMP END DO
  endif 

  ! ****************************************************************************
  ! *** COMPUTE THE INVERSE OF THE AVERAGE AV AT THE U AND V INTERFACES
  if( IGRIDV == 0 )then
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM  
      do K = 1,KS  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          AVUI(L,K) = 2.0/( AV(L,K) + AV(LWC(L),K) )
          AVVI(L,K) = 2.0/( AV(L,K) + AV(LSC(L),K) )
        enddo  
      enddo 
    enddo
    !$OMP END DO
  else
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM  
      do K = 1,KS  
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          AVUI(L,K) = SUB(L)*(1. + SUB3D(L,K))/( AV(L,K)*HP(L) + SUB3D(L,K)*AV(LWC(L),K)*HP(LWC(L) ) )*HU(L)
          AVVI(L,K) = SVB(L)*(1. + SVB3D(L,K))/( AV(L,K)*HP(L) + SVB3D(L,K)*AV(LSC(L),K)*HP(LSC(L) ) )*HV(L)
        enddo  
      enddo 
    enddo
    !$OMP END DO
  endif
  
  ! ****************************************************************************
  ! *** ISSQL:   0 SETS QQ AND QQL STABILITY FUNCTIONS PROPORTIONAL TO
  ! ***            MOMENTUM STABILITY FUNCTIONS
  if( ISTOPT(0) <= 2 .or. ISSQL == 0 )then
    !$OMP DO PRIVATE(ND,K,LP,L)
    do ND = 1,NDM  
      do K = 2,KS  
        do LP = 1,LLWET(K-1,ND)
          L = LKWET(LP,K-1,ND) 
          AQ(L,K) = 0.205*(AV(L,K-1)+AV(L,K))  
        enddo  
      enddo  
      do LP = 1,LLWET(KS,ND) 
        L = LKWET(LP,KS,ND)
        AQ(L,KSZ(L)) = 0.205*AV(L,KSZ(L))  
        AQ(L,KC)     = 0.205*AV(L,KS)  
      enddo  
    enddo
    !$OMP END DO
    
  endif

  ! ****************************************************************************
  ! *** ISSQL:     1  SETS QQ AND QQL STABILITY FUNCTIONS TO CONSTANTS
  ! *** ISTOPT(0)  3  KANTHA (2003)
  if( ISTOPT(0) == 3 .or. ISSQL == 1 )then
    if( ISFAVB == 1 .and. N > 1 )then
      ! *** AVERAGE FILTER
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
      do ND = 1,NDM  
        LF = (ND-1)*LDMWET+1  
        LL = MIN(LF+LDMWET-1,LAWET)
    
        do LP = 1,LLWET(KS,ND)
          L = LKWET(LP,KS,ND) 
          AVOXYHPI(L) = AVOXY(L)*HPI(L)
        enddo
        
        do K = 2,KS  
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            AQ(L,K) = 0.314*( DML(L,K-1)*QQSQR(L,K-1)+DML(L,K)*QQSQR(L,K) ) + AVOXYHPI(L)
          enddo
        enddo

        do LP = LF,LL  
          L = LWET(LP)  
          AQ(L,KSZ(L)) = 0.314*DML(L,KSZ(L))**QQSQR(L,KSZ(L)) + AVOXYHPI(L)
          AQ(L,KC)     = 0.314*DML(L,KC    )* QQSQR(L,KC    ) + AVOXYHPI(L)
        enddo
      enddo
      !$OMP END DO

    else
      ! *** SQRT FILTER & NO FILTER
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
      do ND = 1,NDM  
        LF = (ND-1)*LDMWET+1  
        LL = MIN(LF+LDMWET-1,LAWET)
        
        do K = 2,KS
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            AQ(L,K) = 0.5*( 0.314*( DML(L,K-1)*QQSQR(L,K-1)+DML(L,K  )*QQSQR(L,K) ) + AVOXYHPI(L) ) + 0.5*AQ(L,K)
          enddo
        enddo

        do LP = LF,LL  
          L = LWET(LP)  
          AQ(L,KSZ(L)) = 0.5*( 0.314*DML(L,KSZ(L))**QQSQR(L,KSZ(L)) + AVOXYHPI(L) ) + 0.5*AQ(L,KSZ(L))
          AQ(L,KC)     = 0.5*( 0.314*DML(L,KC    )* QQSQR(L,KC    ) + AVOXYHPI(L) ) + 0.5*AQ(L,KC)
        enddo
      enddo
      !$OMP END DO

    endif
  endif

  !$OMP END PARALLEL
  
  ! ****************************************************************************
  return 

END  
