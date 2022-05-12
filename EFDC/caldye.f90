! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALDYE

  USE GLOBAL

  IMPLICIT NONE   
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DYEF
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: TTHICK
  
  INTEGER ND,MD,L,K,LP,IFIRST
  REAL :: CDYETMP, DAGE, CRIGHT, CLEFT
  
  IF( .NOT. ALLOCATED(DYEF) )THEN
    ALLOCATE(DYEF(LCM,0:KCM))
    ALLOCATE(TTHICK(LCM,KCM))
    DYEF = 0.0
    TTHICK = 0.0
  ENDIF
  
  DAGE = DELT/86400.
  IFIRST = 1
  DO MD = 1,NDYE
    IF( DYES(MD).ITYPE == 0  )THEN
      ! *** CONSERVATIVE DYE
      CYCLE
    
    ELSEIF( DYES(MD).ITYPE == 1 )THEN   
      ! *** NON-CONSERVATIVE DYE
      
      IF( DYES(MD).SETTLE /= 0 .AND. KC > 1 )THEN
        IF( IFIRST == 1 )THEN
          DO K=1,KC
            DO LP=1,LLWET(K+1,ND)
              TTHICK(L,K) = DELT*HPKI(L,K)
            ENDDO
          ENDDO
          IFIRST = 0
        ENDIF
        
        DO LP=1,LLWET(KC,ND)
          L=LKWET(LP,KC,ND)  
          CLEFT = 1.0 + DYES(MD).SETTLE*TTHICK(L,K)
          CRIGHT = MAX(DYE(L,KC,MD),0.0)
          DYE(L,KC,MD) = CRIGHT/CLEFT
          DYEF(L,KC-1) = -DYES(MD).SETTLE*DYE(L,KC,MD)
        ENDDO
           
        DO K=KS,1,-1
          DO LP=1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            CLEFT = 1.0 + DYES(MD).SETTLE*TTHICK(L,K)
            CRIGHT = MAX(DYE(L,K,MD),0.0) - DYEF(L,K)*TTHICK(L,K)
            DYE(L,K,MD) = CRIGHT/CLEFT
            DYEF(L,K-1) = -DYES(MD).SETTLE*DYE(L,K,MD)
          ENDDO
        ENDDO

        !DO LP=1,LLWET(KS,ND)
        !  L=LKWET(LP,KS,ND)
        !  K = KSZ(L)
        !  CLEFT = 1.0 + DYES(MD).SETTLE*TTHICK(L,K)
        !  CRIGHT = MAX(DYE(L,KC,MD),0.0)
        !  DYE(L,KC,MD) = CRIGHT/CLEFT
        !  DYEF(L,KC-1) = -DYES(MD).SETTLE*DYE(L,KC,MD)
        !ENDDO
      ENDIF
    
      ! *** TEMPERATURE
      IF( DYES(MD).TREF > 0. .AND. ISTRAN(2) > 0 )THEN
        ! *** TEMPERATURE DEPENDENT DECAY
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            DYE(L,K,MD) = DYE(L,K,MD) - ( DYES(MD).KRATE0**(TEM(L,K)-DYES(MD).TREF) + DYE(L,K,MD)*DYES(MD).KRATE1**(TEM(L,K)-DYES(MD).TREF) )*DAGE
          ENDDO
        ENDDO  
        
      ELSE
        ! *** TEMPERATURE INDEPENDENT GROWTH/DECAY
        IF( DYES(MD).KRATE1 < 0.0 )THEN
          CDYETMP = EXP(-DYES(MD).KRATE1*DELT)     ! *** EXPONENTIAL DECAY
        ELSE
          CDYETMP = 1./(1.+DELT*DYES(MD).KRATE1)   ! *** GROWTH RATE 
        ENDIF

        !$OMP PARALLEL DO PRIVATE(ND,K,LP,L)
        DO ND=1,NDM  
          DO K=1,KC  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              DYE(L,K,MD) = CDYETMP*DYE(L,K,MD)  
            ENDDO
          ENDDO  
        ENDDO 
        !$OMP END PARALLEL DO
      
      ENDIF
      
    ELSEIF( DYES(MD).ITYPE == 2 )THEN 
      ! *** AGE OF WATER (0th ORDER GROWTH)
      !$OMP PARALLEL DO PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            DYE(L,K,MD) = DYE(L,K,MD) + DAGE
          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
  
    ENDIF 
  ENDDO
  
  RETURN
END
