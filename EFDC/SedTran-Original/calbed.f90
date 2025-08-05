! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALBED  

  ! *** *************************************************************************!
  ! ***  SUBROUTINE CALBED UPDATES THE BED POROSITY, VOID RATIO AND BED THICKNESS
  ! ***  CALLED FROM SSEDTOX  
  ! ***  NOT USED FOR SEDZLJ
  ! *** *************************************************************************!
  !
  !------------------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------!
  ! 2019-01           Paul M. Craig      Added Hard Bottom Bypass
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP
  !------------------------------------------------------------------------------!

  use GLOBAL  
  use Allocate_Initialize      

  implicit none

  integer :: K, L, IFLAG, NS, KK, NX, KBTM1, ND, LF, LL, LP
  real :: TMPVAL, WDENKGM3, WDENGMM3, TMPVALK
  real :: TMPVALKP, BETTMP, VOIDCON1, HBEDTMP, TMPVALO
  real :: TMPVALN, TMPEXP, TMPTOP, TMPBOT, FSTRSE, FDSTRSE
  real :: FHYDCN, DSTRESET, POROTMP

  real,save,allocatable,dimension(:,:) :: SNDHYDCN  

  if( .not. allocated(SNDHYDCN) )then
    call AllocateDSI(SNDHYDCN, LCM, KBM, 0.0)
  endif

  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then  
    HBEDMIN = 1.E-4  
    if( ISTRAN(7) >= 1 )then  
      HBEDMIN = max(HBEDMIN,SNDDMX)  
    endif  

    ! *********************************************************************
    ! *** CONSTANT POROSITY BED  
    if( IBMECH == 0 )then

      ! *** UPDATE TOP LAYER THICKNESS TO MAINTAIN CONSTANT POROSITY  
      VOIDCON1 = BEDPORC/(1.-BEDPORC)  

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,LP,L,LF,LL)              &
      !$OMP             PRIVATE(HBEDTMP,TMPVALO,TMPVALN)
      do ND = 1,NDM  
        LF = (ND-1)*LDMSED + 1  
        LL = min(LF+LDMSED-1,LASED)

        do LP = LF,LL  
          L = LSED(LP)
          K = KBT(L)  
          HBEDTMP      = (1. + VOIDCON1)*HBED(L,K)/(1. + VDRBED(L,K))  
          TMPVALO      = VDRBED(L,K)*HBED(L,K)/(1. + VDRBED(L,K))  
          TMPVALN      = VOIDCON1*HBEDTMP/(1. + VOIDCON1)  
          QWBDTOP(L)   = DELTI*(TMPVALO-TMPVALN)  
          HBED(L,K)    = HBEDTMP  
          QWTRBED(L,K) = QWBDTOP(L) + QGW(L)/DXYP(L) 
          VDRBED(L,K)  = VOIDCON1     ! *** PMC 2010-12-06
        enddo  
        do K = 0,KBT(L)-1  
          do LP = LF,LL  
            L = LSED(LP)
            QWTRBED(L,K) = QGW(L)/DXYP(L)  
          enddo  
        enddo  

      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END PARALLEL DO
      
    endif  

    ! *********************************************************************
    ! *** SIMPLE CONSOLIDATING BED  
    if( IBMECH == 1 )then

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,KK, LP,L,LF,LL)              &
      !$OMP             PRIVATE(TMPEXP,TMPTOP,TMPBOT,TMPVAL)
      do ND = 1,NDM  
        LF = (ND-1)*LDMSED + 1  
        LL = min(LF+LDMSED-1,LASED)

      ! ***  IF SEDVRDT > 0.0 CONSOLIDATE TO SEDVRM (THE MINIMUM VOID RATIO)  
        if( SEDVRDT > 0.0 )then  
          TMPEXP = EXP(-DTSED/SEDVRDT)  
          do K = 1,KB  
            do LP = LF,LL  
              L = LSED(LP)
              if( K <= KBT(L) )then  
                VDRBED1(L,K) = VDRBED(L,K)  
                HBED1(L,K)   = HBED(L,K)  
                VDRBED(L,K)  = SEDVDRM + (VDRBED1(L,K)-SEDVDRM)*TMPEXP  
                TMPTOP = 1. + VDRBED(L,K)  
                TMPBOT = 1. + VDRBED1(L,K)  
                HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
              endif  
            enddo  
          enddo  
        endif  

        ! ***  IF SEDVRDT = 0.0 - ALLOW DEPOSITION AND MIXING TO MODIFY VOID RATIOS  
        if( SEDVRDT == 0.0 )then  
          do K = 1,KB  
            do LP = LF,LL  
              L = LSED(LP)
              if( K <= KBT(L) )then  
                VDRBED1(L,K) = VDRBED(L,K)  
                HBED1(L,K)   = HBED(L,K)  
                VDRBED(L,K)  = SEDVDRM + (VDRBED1(L,K)-SEDVDRM)  
                TMPTOP = 1. + VDRBED(L,K)  
                TMPBOT = 1. + VDRBED1(L,K)  
                HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
              endif  
            enddo  
          enddo  
        endif  

        ! ***  IF SEDVRDT>0.0001 CONSOLIDATE TO SEDVRM INSTANTANEOUSLY  
        !IF( SEDVRDT >= 0.0 .and. SEDVRDT <= 0.0001 )then  
        !  TMPEXP = 0.0  
        !  do K = 1,KB  
        !    do LP = LF,LL  
        !      L = LSED(LP)
        !      if( K <= KBT(L) )then  
        !        VDRBED1(L,K) = VDRBED(L,K)  
        !        HBED1(L,K) = HBED(L,K)  
        !        VDRBED(L,K) = SEDVDRM  
        !        TMPTOP = 1. + VDRBED(L,K)  
        !        TMPBOT = 1. + VDRBED1(L,K)  
        !        HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
        !      endif  
        !    enddo  
        !  enddo  
        !ENDIF  

        ! ***  IF SEDVRDT < 0.0 MAINTAIN INITIAL VOID RATIO (SAVED IN VDRBED0)  
        if( SEDVRDT < 0.0 )then  
          TMPEXP = 1.0  
          do LP = LF,LL  
            L = LSED(LP)
            KK = max(1,KBT(L)-2)
            do K = KK,KBT(L)
              VDRBED1(L,K) = VDRBED(L,K)                        ! *** void ratio has been updated based on sediment mass flux in SSEDTOX
              HBED1(L,K) = HBED(L,K)  
              !VDRBED(L,K) = VDRBED0(L,K)                       ! *** Set the current void ratio to the initial condition VR
              !TMPTOP = 1. + VDRBED(L,K)  
              !TMPBOT = 1. + VDRBED1(L,K)  
              !HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT             ! *** Update bed layer thickness
            enddo
          enddo  
        endif  

        ! *** UPDATE POROSITY  
        do K = 1,KB  
          do LP = LF,LL  
            L = LSED(LP)
            if( K <= KBT(L) )then  
              PORBED(L,K)  = VDRBED(L,K) /(1. + VDRBED(L,K))  
              PORBED1(L,K) = VDRBED1(L,K)/(1. + VDRBED1(L,K))  
            endif  
          enddo  
        enddo  

        ! *** UPDATE PORE WATER FLOWS  
        do LP = LF,LL  
          L = LSED(LP)
          QWTRBED(L,0) = QGW(L)/DXYP(L)  
        enddo  
        do K = 1,KB  
          do LP = LF,LL  
            L = LSED(LP)
            if( K <= KBT(L) )then  
              TMPVAL = HBED(L,K)/(1. + VDRBED(L,K))  
              QWTRBED(L,K) = QWTRBED(L,K-1) - DELTI*TMPVAL*(VDRBED(L,K)-VDRBED1(L,K))  
            endif  
          enddo  
        enddo  
      enddo   ! *** END OF DOMAIN LOOP
      !$OMP END PARALLEL DO

    endif  

    ! *********************************************************************
    ! *** UPDATE BULK DENSITY AND TOTAL SEDIMENT MASS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,K,L,NS) 
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED + 1  
      LL = min(LF+LDMSED-1,LASED)

      if( ISTRAN(6) >= 1 )then
        do K = 1,KB
          do LP = LF,LL  
            L = LSED(LP)
            SEDBT(L,K) = 0.0     
          enddo
        enddo

        do NS = 1,NSED
          do K = 1,KB
            do LP = LF,LL  
              L = LSED(LP)
              SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
            enddo
          enddo
        enddo
      endif

      if( ISTRAN(7) >= 1 )then
        do K = 1,KB
          do LP = LF,LL  
            L = LSED(LP)
            SNDBT(L,K) = 0.0
          enddo
        enddo

        do NS = 1,NSND
          do K = 1,KB
            do LP = LF,LL  
              L = LSED(LP)
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NS)
            enddo
          enddo
        enddo
      endif
      
      do K = 1,KB
        do LP = LF,LL  
          L = LSED(LP)
          if( HBED(L,K) > 0. )then
            ! *** UPDATE TOTAL/WET DENSITY
            BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)
          else
            BDENBED(L,K) = 0.
          endif
        enddo
      enddo
    enddo  ! END OF DOMAIN
    !$OMP END PARALLEL DO

  endif

8011 FORMAT('CBED',I12,6E13.5)    

  ! *** *******************************************************************!

  return 
  END  

