! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALDYE

  use GLOBAL
  use Allocate_Initialize
  
  implicit none   
  
  real,save,allocatable,dimension(:,:) :: DYEF
  real,save,allocatable,dimension(:,:) :: TTHICK
  
  integer ND, MD, L, K, LP, IFIRST
  real :: CDYETMP, DAGE, CRIGHT, CLEFT, DYEMW, VOLTERM
  real, save :: DYESTEP
  real, external :: CALVOLTERM
  
  if( .not. allocated(DYEF) )then
    call AllocateDSI(DYEF,   LCM, -KCM, 0.0)
    call AllocateDSI(TTHICK, LCM, -KCM, 0.0)
    DYESTEP = 0.0
  endif
  
  DYESTEP = DYESTEP + DELT
  
  IFIRST = 1
  do MD = 1,NDYE
    ! *** Select approach for each dye class
    if( DYES(MD).ITYPE == 0  )then
      ! *** Conservative dye
      CYCLE
    
    elseif( DYES(MD).ITYPE == 1 )then   
      ! *** Non-Conservative dye

      ! *** Settling/Rising
      if( DYES(MD).SETTLE /= 0 .and. KC > 1 )then
        if( IFIRST == 1 )then
          do K = 1,KC
            do LP = 1,LLWET(K+1,ND)
              TTHICK(L,K) = DELT*HPKI(L,K)
            enddo
          enddo
          IFIRST = 0
        endif
        
        do LP = 1,LLWET(KC,ND)
          L = LKWET(LP,KC,ND)  
          CLEFT = 1.0 + DYES(MD).SETTLE*TTHICK(L,K)
          CRIGHT = max(DYE(L,KC,MD),0.0)
          DYE(L,KC,MD) = CRIGHT/CLEFT
          DYEF(L,KC-1) = -DYES(MD).SETTLE*DYE(L,KC,MD)
        enddo
           
        do K = KS,1,-1
          do LP = 1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND) 
            CLEFT = 1.0 + DYES(MD).SETTLE*TTHICK(L,K)
            CRIGHT = max(DYE(L,K,MD),0.0) - DYEF(L,K)*TTHICK(L,K)
            DYE(L,K,MD) = CRIGHT/CLEFT
            DYEF(L,K-1) = -DYES(MD).SETTLE*DYE(L,K,MD)
          enddo
        enddo
      endif
    
      ! *** Kinetic processes
      if( DYESTEP >= DYESTEPW )then
        
        ! *** Temperature
        if( DYES(MD).TREF > 0. .and. ISTRAN(2) > 0 .and. DYES(MD).KRATE0 > 0.0 )then
          ! *** Temperature dependent decay
          DAGE = DYESTEP/86400.
          do K = 1,KC  
            do LP = 1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              DYE(L,K,MD) = DYE(L,K,MD) - ( DYES(MD).KRATE0**(TEM(L,K)-DYES(MD).TREF) + DYE(L,K,MD)*DYES(MD).KRATE1**(TEM(L,K)-DYES(MD).TREF) )*DAGE
            enddo
          enddo  
        
        elseif( DYES(MD).KRATE1 /= 0.0 )then
          ! *** Temperature independent growth/decay
          if( DYES(MD).KRATE1 < 0.0 )then
            CDYETMP = EXP(-DYES(MD).KRATE1*DYESTEP)     ! *** Exponential decay
          else
            CDYETMP = 1./(1.+DYESTEP*DYES(MD).KRATE1)   ! *** Growth rate 
          endif

          !$OMP PARALLEL DO PRIVATE(ND,K,LP,L)
          do ND = 1,NDM  
            do K = 1,KC  
              do LP = 1,LLWET(K,ND)
                L = LKWET(LP,K,ND)  
                DYE(L,K,MD) = CDYETMP*DYE(L,K,MD)  
              enddo
            enddo  
          enddo 
          !$OMP END PARALLEL DO
      
        endif
      
        ! *** Compute volatilization
        if( DYES(MD).VOL.KL_OPT > 0 .and. ( ISTRAN(2) > 0 .or. IVOLTEMP > 0 ) )then
          DYEMW = 1./DYES(MD).VOL.MW**0.66667
        
          !$OMP PARALLEL DO PRIVATE(ND, LP, L, VOLTERM)
          do ND = 1,NDM  
            do LP = 1,LLWET(KC,ND)
              L = LKWET(LP,KC,ND)  
              
              VOLTERM = CALVOLTERM(HP(L), HPKI(L,KC), STCUV(L), UHE(L), VHE(L), HU(L), HV(L), TEM(L,KC), SAL(L,KC), &
                    TATMT(L), PATMT(L), WINDST(L), VOL_VEL_MAX, VOL_DEP_MIN, DYEMW, DYES(MD).VOL.HE,                &
                    DYES(MD).VOL.AIRCON, DYES(MD).VOL.TCOEFF, DYES(MD).VOL.MULT, DYE(L,KC,MD), DYES(MD).VOL.KL_OPT, 0.0)   
          
              DYE(L,KC,MD) = DYE(L,KC,MD) - VOLTERM*DYESTEP
              DYE(L,KC,MD) = max(DYE(L,KC,MD), 0.0)
            enddo
          enddo 
          !$OMP END PARALLEL DO
        endif
        
        DYESTEP = 0.0
      endif
      
    elseif( DYES(MD).ITYPE == 2 )then 
      ! *** Age of Water (0th order growth)
      DAGE = DELT/86400.
      
      !$OMP PARALLEL DO PRIVATE(ND,K,LP,L)
      do ND = 1,NDM  
        do K = 1,KC  
          do LP = 1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            DYE(L,K,MD) = DYE(L,K,MD) + DAGE
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif 
  enddo  ! *** END OF DYE LOOP (MD)
  
  return
END
