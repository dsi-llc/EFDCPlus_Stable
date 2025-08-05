! ----------------------------------------------------------------------
!   This file is a part of EFDC + 
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE CALBLAY  

  ! *** *******************************************************************!
  ! *** SUBROUTINE CALBLAY REMOVES OR ADDS LAYERS TO THE SEDIMENT BED
  ! *** NOT USED FOR SEDZLJ
  !  
  ! *** *******************************************************************!
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2019-01           Paul M. Craig      Added Hard Bottom Bypass
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP
  !----------------------------------------------------------------------!

  use GLOBAL  

  implicit none

  integer :: K, NS, L, NT, NX, ND, LF, LL, LP, k1
  real   :: TMPBOT2, TMPTOP1, TMPTOP2, TMPVAL, HBEDMXT, HOLDTOP, FKBTP
  real   :: SEDBOLD, TOXBOLD, TMPBOT1, FKBT, SNDBOLD, HTMP1, HTMP2, HBEDMAX1

  ! *** FOR TRANSPORT OF COHESIVE SEDIMENT ONLY SET HBEDMIN TO FRACTION
  ! *** OF HBEDMAX
  if( ISTRAN(7) == 0 )then
    HBEDMIN = MAX(0.01,0.01*HBEDMAX)
  else
    HBEDMIN = 1.0*MAXVAL(SEDDIA(:))
  endif

  ! *** WHEN NONCOHESIVE TRANSPORT IS ACTIVE, WITHOUT ACTIVE-PARENT LAYER
  ! *** FORMULATION, SET HBEDMIN PROPORTIONAL TO MAXIMUM GRAIN DIAMETER
  if( ISTRAN(7) >= 1 .and. ISNDAL <= 1 )then
    ! *** NOTE HQI CHANGED ORIGINAL FORMULATION FROM
    ! *** TMPVAL = 2.*SNDDMX
    ! *** TO
    TMPVAL  = SNDDMX
    HBEDMIN = MAX(HBEDMIN,TMPVAL)
  endif

  ! *** WHEN NONCOHESIVE TRANSPORT IS ACTIVE, WITH ACTIVE-PARENT LAYER
  ! *** FORMULATION, SET HBEDMIN PROPORTIONAL ACTIVE LAYER THICKNESS
  if( ISNDAL == 2 )then
    HBEDMIN = MAX(HBEDMIN,0.1*HBEDAL)
    HBEDMAX1 = 1.1*HBEDAL
  endif

  HBEDMXT = HBEDMAX + HBEDMIN  
  
  ! *** ADD OR REMOVE TOP LAYER (NO ACTIVE LAYER ARMORING OPTION)  
  if( ISNDAL <= 1 )then
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,LP,L,LF,LL,NS,NX,NT)              &
    !$OMP             PRIVATE(TMPBOT1,TMPBOT2,TMPTOP1,TMPTOP2,TMPVAL,HOLDTOP) &
    !$OMP             PRIVATE(SEDBOLD,SNDBOLD,TOXBOLD,FKBT,FKBTP,HTMP1,HTMP2)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED-1,LASED)

      ! *** ADD NEW TOP LAYER  
      do LP = LF,LL  
        L = LSED(LP)
        if( HBED(L,KBT(L)) > HBEDMXT )then  
          if( KBT(L) < KB )then  
            HOLDTOP = HBED(L,KBT(L))
            HTMP1 = HOLDTOP-HBEDMAX
            HTMP2 = HBEDMAX
            HBED(L,KBT(L)+1) = MIN(HTMP1,HTMP2)   ! *** Thinner layer on top
            HBED(L,KBT(L)) = MAX(HTMP1,HTMP2)
            if( IBMECH == 1 .and. SEDVRDT < 0.0 )then     ! *** These control maintianing initial void ratio k profile
              VDRBED0(L,KBT(L)+1) = VDRBED0(L,KBT(L))     ! *** Reset void ratio profile
              SDENAVG(L,KBT(L)+1) = SDENAVG(L,KBT(L))     ! *** Reset avg sediment solids density
            endif                                    
            VDRBED(L,KBT(L)+1) = VDRBED(L,KBT(L))  
            PORBED(L,KBT(L)+1) = PORBED(L,KBT(L))  
            STDOCB(L,KBT(L)+1) = STDOCB(L,KBT(L))  
            STPOCB(L,KBT(L)+1) = STPOCB(L,KBT(L))  
            FKBTP = HBED(L,KBT(L)+1)/HOLDTOP  
            FKBT = HBED(L,KBT(L))/HOLDTOP
              
            if( ISTRAN(6) >= 1 )then  
              SEDBT(L,KBT(L)+1) = 0.  
              SEDBT(L,KBT(L)) = 0.  
              do NS = 1,NSED  
                SEDBOLD = SEDB(L,KBT(L),NS)  
                SEDB(L,KBT(L)+1,NS) = FKBTP*SEDBOLD  
                SEDB(L,KBT(L),NS) = FKBT*SEDBOLD  
                SEDBT(L,KBT(L)+1) = SEDBT(L,KBT(L)+1) + SEDB(L,KBT(L)+1,NS)  
                SEDBT(L,KBT(L)) = SEDBT(L,KBT(L)) + SEDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NS) = STFPOCB(L,KBT(L),NS)  
                VFRBED(L,KBT(L)+1,NS) = VFRBED(L,KBT(L),NS)  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              SNDBT(L,KBT(L)+1) = 0.  
              SNDBT(L,KBT(L)) = 0.  
              do NS = 1,NSND  
                NX = NS + NSED  
                SNDBOLD = SNDB(L,KBT(L),NS)  
                SNDB(L,KBT(L)+1,NS) = FKBTP*SNDBOLD  
                SNDB(L,KBT(L),NS) = FKBT*SNDBOLD  
                SNDBT(L,KBT(L)+1) = SNDBT(L,KBT(L)+1) + SNDB(L,KBT(L)+1,NS)  
                SNDBT(L,KBT(L)) = SNDBT(L,KBT(L)) + SNDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NX) = STFPOCB(L,KBT(L),NX)  
                VFRBED(L,KBT(L)+1,NX) = VFRBED(L,KBT(L),NX)  
              enddo  
            endif  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXBOLD = TOXB(L,KBT(L),NT)  
                TOXB(L,KBT(L)+1,NT) = FKBTP*TOXBOLD  
                TOXB(L,KBT(L),NT) = FKBT*TOXBOLD  
              enddo  
            endif  
            KBT(L) = KBT(L)+1  
          endif  
        endif  
      enddo  

      ! *** REZONE WITH NEW TOP LAYER ADDED NEXT TIME STEP  
      do LP = LF,LL  
        L = LSED(LP)
        if( HBED(L,KBT(L)) > HBEDMXT )then  
          if( KBT(L) == KB .and. KB > 1 )then  
            TMPBOT1 = HBED(L,1)/(1. + VDRBED(L,1))  
            TMPBOT2 = HBED(L,2)/(1. + VDRBED(L,2))  
            TMPTOP1 = VDRBED(L,1)*TMPBOT1  
            TMPTOP2 = VDRBED(L,2)*TMPBOT2  
            VDRBED(L,1) = (TMPTOP1 + TMPTOP2)/(TMPBOT1 + TMPBOT2)  
            PORBED(L,1) = VDRBED(L,1)/(1. + VDRBED(L,1))  
            HBED(L,1) = HBED(L,1) + HBED(L,2)  
            if( KB == 2 )then  
              HBED(L,2)   = 0  
              VDRBED(L,2) = 0.0  
              PORBED(L,2) = 0.0  
              STDOCB(L,2) = 0.0  
              STPOCB(L,2) = 0.0  
            endif  
            if( KB > 2 )then  
              do K = 2,KBT(L)-1  
                HBED(L,K)   = HBED(L,K+1)  
                VDRBED(L,K) = VDRBED(L,K+1)  
                PORBED(L,K) = PORBED(L,K+1)  
                STDOCB(L,K) = STDOCB(L,K+1)  
                STPOCB(L,K) = STPOCB(L,K+1)  
              enddo  
              HBED(L,KBT(L))   = 0  
              VDRBED(L,KBT(L)) = 0.0  
              PORBED(L,KBT(L)) = 0.0  
              STDOCB(L,KBT(L)) = 0.0  
              STPOCB(L,KBT(L)) = 0.0  
            endif  
            if( ISTRAN(6) >= 1 )then  
              do NS = 1,NSED  
                SEDB(L,1,NS) = SEDB(L,1,NS) + SEDB(L,2,NS)  
                if( KB == 2 )then  
                  SEDB(L,2,NS) = 0.0  
                  STFPOCB(L,2,NS) = 0.0  
                  VFRBED(L,2,NS) = 0.0  
                endif  
                if( KB > 2 )then  
                  do K = 2,KBT(L)-1  
                    SEDB(L,K,NS) = SEDB(L,K+1,NS)  
                    STFPOCB(L,K,NS) = STFPOCB(L,K+1,NS)  
                    VFRBED(L,K,NS) = VFRBED(L,K+1,NS)  
                  enddo  
                  SEDB(L,KBT(L),NS) = 0  
                  STFPOCB(L,KBT(L),NS) = 0.0  
                  VFRBED(L,KBT(L),NS) = 0.0  
                endif  
              enddo  
              do K = 1,KB  
                SEDBT(L,K) = 0.0  
              enddo  
              do NS = 1,NSED  
                do K = 1,KB  
                  SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)  
                enddo  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              do NS = 1,NSND  
                NX = NS + NSED  
                SNDB(L,1,NS) = SNDB(L,1,NS) + SNDB(L,2,NS)  
                if( KB == 2 )then  
                  SNDB(L,2,NS) = 0.0  
                  STFPOCB(L,2,NX) = 0.0  
                  VFRBED(L,2,NX) = 0.0  
                endif  
                if( KB > 2 )then  
                  do K = 2,KBT(L)-1  
                    SNDB(L,K,NS) = SNDB(L,K+1,NS)  
                    STFPOCB(L,K,NX) = STFPOCB(L,K+1,NX)  
                    VFRBED(L,K,NX) = VFRBED(L,K+1,NX)  
                  enddo  
                  SNDB(L,KBT(L),NS) = 0  
                  STFPOCB(L,KBT(L),NX) = 0.0  
                  VFRBED(L,KBT(L),NX) = 0.0  
                endif  
              enddo  
              do K = 1,KB  
                SNDBT(L,K) = 0.0  
              enddo  
              do NS = 1,NSND  
                do K = 1,KB  
                  SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NS)  
                enddo  
              enddo  
            endif  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXB(L,1,NT) = TOXB(L,1,NT) + TOXB(L,2,NT)  
                if( KB == 2 )then  
                  TOXB(L,2,NT) = 0  
                endif  
                if( KB > 2 )then  
                  do K = 2,KBT(L)-1  
                    TOXB(L,K,NT) = TOXB(L,K+1,NT)  
                  enddo  
                  TOXB(L,KBT(L),NT) = 0  
                endif  
              enddo  
            endif  
            if( KB == 2 )then  
              KBT(L) = 1  
            endif  
            if( KB > 2 )then  
              KBT(L) = KBT(L)-1  
            endif  
          endif  
        endif  
      enddo  

      ! *** REMOVE TOP LAYER  
      do LP = LF,LL  
        L = LSED(LP)
        if( HBED(L,KBT(L)) < HBEDMIN )then  
          if( KBT(L) > 1 )then  
            TMPBOT1 = HBED(L,KBT(L)-1)/(1. + VDRBED(L,KBT(L)-1))  
            TMPBOT2 = HBED(L,KBT(L))/(1. + VDRBED(L,KBT(L)))  
            TMPTOP1 = VDRBED(L,KBT(L)-1)*TMPBOT1  
            TMPTOP2 = VDRBED(L,KBT(L))*TMPBOT2  
            VDRBED(L,KBT(L)-1) = (TMPTOP1 + TMPTOP2)/(TMPBOT1 + TMPBOT2)  
            PORBED(L,KBT(L)-1) = VDRBED(L,KBT(L)-1)/(1. + VDRBED(L,KBT(L)-1))  
            HBED(L,KBT(L)-1) = HBED(L,KBT(L)-1) + HBED(L,KBT(L))  
            HBED(L,KBT(L)) = 0.0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  
            BDENBED(L,KBT(L)) = 0.  ! PMC
            if( ISTRAN(6) >= 1 )then  
              SEDBT(L,KBT(L)-1) = 0.  
              SEDBT(L,KBT(L)) = 0.  
              do NS = 1,NSED  
                SEDB(L,KBT(L)-1,NS) = SEDB(L,KBT(L)-1,NS) + SEDB(L,KBT(L),NS)  
                SEDBT(L,KBT(L)-1) = SEDBT(L,KBT(L)-1)    +SEDB(L,KBT(L)-1,NS)  
                SEDB(L,KBT(L),NS) = 0.0  
                STFPOCB(L,KBT(L),NS) = 0.0  
                VFRBED(L,KBT(L),NS) = 0.0  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              SNDBT(L,KBT(L)-1) = 0.  
              SNDBT(L,KBT(L)) = 0.  
              do NS = 1,NSND  
                SNDB(L,KBT(L)-1,NS) = SNDB(L,KBT(L)-1,NS) + SNDB(L,KBT(L),NS)  
                SNDBT(L,KBT(L)-1) = SNDBT(L,KBT(L)-1) +SNDB(L,KBT(L)-1,NS)  
                SNDB(L,KBT(L),NS) = 0.0  
                STFPOCB(L,KBT(L),NS + NSED) = 0.0  
                VFRBED(L,KBT(L),NS + NSED) = 0.0  
              enddo  
            endif  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXB(L,KBT(L)-1,NT) = TOXB(L,KBT(L)-1,NT) + TOXB(L,KBT(L),NT)  
                TOXB(L,KBT(L),NT) = 0.0  
              enddo  
            endif  
            KBT(L) = KBT(L)-1  

          ! *** PMC BEGIN BLOCK  
          elseif( HBED(L,KBT(L)) < 0.0 )then  
            ! *** ZERO NEGATIVE THICKNESSES
            HBED(L,KBT(L)) = 0.0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  
          endif  
          ! *** PMC END BLOCK  
        endif  
      enddo  

      ! *** UPDATE BULK DENSITY  
      do LP = LF,LL  
        L = LSED(LP)
        ! ***REMOVED KB LOOP, ONLY COMPUTE THE TOP LAYER.  CONSOLIDATION IS HANDLED LATER IN CALBED
        K = KBT(L)  
        if( HBED(L,K) > 0. )then  
          ! *** UPDATE TOTAL/WET DENSITY
          BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)  
        else  
          BDENBED(L,K) = 0.  
        endif  
      enddo  

    enddo  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
         
    ! ***************************************************************************        
  else  
    ! ***************************************************************************        
    ! *** ADD OR REMOVE PARENT LAYER WHEN ARMORING IS ACTIVE  (ISNDAL == 2 )
  
    !  $OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,NT)        &
    !  $OMP             PRIVATE(TMPBOT1,TMPBOT2,TMPTOP1,TMPTOP2,TMPVAL,HOLDTOP)  &
    !  $OMP             PRIVATE(SEDBOLD,TOXBOLD,FKBT,FKBTP,SNDBOLD,HTMP1,HTMP2,  k1)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED + 1  
      LL = MIN(LF + LDMSED-1,LASED)

      ! *** CHECK TO SEE IF NEED TO ADD NEW LAYER BELOW THE ACTIVE LAYER 
      do LP = LF,LL
        L = LSED(LP)
        K = KBT(L) - 1
        if( K == 0 ) CYCLE             ! *** CYCLE IF NO PARENT LAYER
        
        if( HBED(L,K) > HBEDMXT )then  
          ! *** Layer is thicker than the maximum
          if( KBT(L) < KB )then  
            ! *** MOVE ACTIVE LAYER UP  
            HBED(L,KBT(L)+1) = HBED(L,KBT(L))  
            VDRBED(L,KBT(L)+1) = VDRBED(L,KBT(L))  
            PORBED(L,KBT(L)+1) = PORBED(L,KBT(L))  
            STDOCB(L,KBT(L)+1) = STDOCB(L,KBT(L))  
            STPOCB(L,KBT(L)+1) = STPOCB(L,KBT(L))  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXB(L,KBT(L)+1,NT) = TOXB(L,KBT(L),NT)  
              enddo  
            endif  
            if( ISTRAN(6) >= 1 )then  
              SEDBT(L,KBT(L)+1) = SEDBT(L,KBT(L))  
              do NS = 1,NSED  
                SEDB(L,KBT(L)+1,NS) = SEDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NS) = STFPOCB(L,KBT(L),NS)  
                VFRBED(L,KBT(L)+1,NS) = VFRBED(L,KBT(L),NS)  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              SNDBT(L,KBT(L)+1) = SNDBT(L,KBT(L))  
              do NS = 1,NSND  
                NX = NS + NSED  
                SNDB(L,KBT(L)+1,NS) = SNDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NX) = STFPOCB(L,KBT(L),NX)  
                VFRBED(L,KBT(L)+1,NX) = VFRBED(L,KBT(L),NX)  
              enddo  
            endif  

            ! *** SPLIT PARENT INTO TWO LAYERS  
            HOLDTOP = HBED(L,KBT(L)-1)  
            HBED(L,KBT(L)) = HOLDTOP-HBEDMAX  
            HBED(L,KBT(L)-1) = HBEDMAX  
            VDRBED(L,KBT(L)) = VDRBED(L,KBT(L)-1)  
            PORBED(L,KBT(L)) = PORBED(L,KBT(L)-1)  
            STDOCB(L,KBT(L)) = STDOCB(L,KBT(L)-1)  
            STPOCB(L,KBT(L)) = STPOCB(L,KBT(L)-1)  
            FKBTP = HBED(L,KBT(L))/HOLDTOP  
            FKBT = HBED(L,KBT(L)-1)/HOLDTOP  
            if( ISTRAN(6) >= 1 )then  
              SEDBT(L,KBT(L)) = 0.  
              SEDBT(L,KBT(L)-1) = 0.  
              do NS = 1,NSED  
                SEDBOLD = SEDB(L,KBT(L)-1,NS)  
                SEDB(L,KBT(L),NS) = FKBTP*SEDBOLD  
                SEDB(L,KBT(L)-1,NS) = FKBT*SEDBOLD  
                SEDBT(L,KBT(L)) = SEDBT(L,KBT(L))  +SEDB(L,KBT(L),NS)  
                SEDBT(L,KBT(L)-1) = SEDBT(L,KBT(L)-1)  +SEDB(L,KBT(L)-1,NS)  
                STFPOCB(L,KBT(L),NS) = STFPOCB(L,KBT(L)-1,NS)  
                VFRBED(L,KBT(L),NS) = VFRBED(L,KBT(L)-1,NS)  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              SNDBT(L,KBT(L)) = 0.  
              SNDBT(L,KBT(L)-1) = 0.  
              do NS = 1,NSND  
                NX = NS + NSED  
                SNDBOLD = SNDB(L,KBT(L)-1,NS)  
                SNDB(L,KBT(L),NS)   = FKBTP*SNDBOLD  
                SNDB(L,KBT(L)-1,NS) = FKBT*SNDBOLD  
                SNDBT(L,KBT(L))   = SNDBT(L,KBT(L)) + SNDB(L,KBT(L),NS)  
                SNDBT(L,KBT(L)-1) = SNDBT(L,KBT(L)-1) + SNDB(L,KBT(L)-1,NS)  
                STFPOCB(L,KBT(L),NX) = STFPOCB(L,KBT(L)-1,NX)  
                VFRBED(L,KBT(L),NX)  = VFRBED(L,KBT(L)-1,NX)  
              enddo  
            endif  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXBOLD = TOXB(L,KBT(L)-1,NT)  
                TOXB(L,KBT(L),NT) = FKBTP*TOXBOLD  
                TOXB(L,KBT(L)-1,NT) = FKBT*TOXBOLD  
              enddo  
            endif  
            KBT(L) = KBT(L) + 1

          elseif( KB > 2 )then
            ! *** KBT(L) = KB,  REZONE LAYERS BELOW ACTIVE LAYER
            
            ! *** COMBINE THE BOTTOM TWO LAYERS
            TMPBOT1 = HBED(L,1)/(1. + VDRBED(L,1))  
            TMPBOT2 = HBED(L,2)/(1. + VDRBED(L,2))  
            TMPTOP1 = VDRBED(L,1)*TMPBOT1  
            TMPTOP2 = VDRBED(L,2)*TMPBOT2  
            VDRBED(L,1) = (TMPTOP1 + TMPTOP2)/(TMPBOT1 + TMPBOT2)  
            PORBED(L,1) = VDRBED(L,1)/(1. + VDRBED(L,1))  
            HBED(L,1) = HBED(L,1) + HBED(L,2)
            
            ! *** MOVE ALL SEDIMENT DOWN ONE LAYER
            do K = 2,KBT(L)-1  
              HBED(L,K)   = HBED(L,K+1)  
              VDRBED(L,K) = VDRBED(L,K+1)  
              PORBED(L,K) = PORBED(L,K+1)  
              STDOCB(L,K) = STDOCB(L,K+1)  
              STPOCB(L,K) = STPOCB(L,K+1)  
            enddo  
            
            ! *** CREATE A NEW EMPTY KB LAYER
            HBED(L,KBT(L)) = 0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  

            ! *** REASSIGN MASSES
            if( ISTRAN(6) >= 1 )then  
              do NS = 1,NSED  
                SEDB(L,1,NS) = SEDB(L,1,NS) + SEDB(L,2,NS)  
                do K = 2,KBT(L)-1  
                  SEDB(L,K,NS) = SEDB(L,K+1,NS)  
                  STFPOCB(L,K,NS) = STFPOCB(L,K+1,NS)  
                  VFRBED(L,K,NS) = VFRBED(L,K+1,NS)  
                enddo  
                SEDB(L,KBT(L),NS) = 0.0
                STFPOCB(L,KBT(L),NS) = 0.0  
                VFRBED(L,KBT(L),NS) = 0.0  
              enddo  
              do K = 1,KB  
                SEDBT(L,K) = 0.0  
              enddo  
              do NS = 1,NSED  
                do K = 1,KB  
                  SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)  
                enddo  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              do NS = 1,NSND  
                NX = NS + NSED  
                SNDB(L,1,NS) = SNDB(L,1,NS) + SNDB(L,2,NS)  
                do K = 2,KBT(L)-1  
                  SNDB(L,K,NS) = SNDB(L,K+1,NS)  
                  STFPOCB(L,K,NX) = STFPOCB(L,K+1,NX)  
                  VFRBED(L,K,NX) = VFRBED(L,K+1,NX)  
                enddo  
                SNDB(L,KBT(L),NS) = 0.0
                STFPOCB(L,KBT(L),NX) = 0.0  
                VFRBED(L,KBT(L),NX) = 0.0  
              enddo  
              do K = 1,KB  
                SNDBT(L,K) = 0.0  
              enddo  
              do NS = 1,NSND  
                do K = 1,KB  
                  SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NS)  
                enddo  
              enddo  
            endif  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXB(L,1,NT) = TOXB(L,1,NT) + TOXB(L,2,NT)  
                do K = 2,KBT(L)-1  
                  TOXB(L,K,NT) = TOXB(L,K+1,NT)  
                enddo  
                TOXB(L,KBT(L),NT) = 0  
              enddo  
            endif  
            KBT(L) = KBT(L)-1  
          endif  
        endif   
      enddo  
      
      ! *** LIMIT LAYER THICKNESSES FOR ACTIVE LAYERS AFTER SCOUR EVENTS
      do LP = LF,LL  
        L = LSED(LP)
        K = KBT(L)
        if( K == 1 )then
          if( HBED(L,K) > HBEDMAX1 )then  
            TMPBOT1 = HBEDAL/HBED(L,K)
            TMPBOT2 = 1.0 - TMPBOT1
          
            TMPTOP1 = VDRBED(L,KBT(L)-2)*TMPBOT1  
            TMPTOP2 = VDRBED(L,KBT(L)-1)*TMPBOT2
          
            VDRBED(L,K+1) = VDRBED(L,K)  
            PORBED(L,K+1) = PORBED(L,K)  
            STDOCB(L,K+1) = STDOCB(L,K)   
            STPOCB(L,K+1) = STPOCB(L,K) 
            BDENBED(L,K+1) = BDENBED(L,K)    ! delme - not updated anywhere else
            HBED(L,K+1) = HBEDAL
            HBED(L,K)   = HBED(L,K) - HBEDAL

            ! *** SPLIT MASS
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TMPVAL = TOXB(L,K,NT)
                TOXB(L,K+1,NT) = TMPVAL*TMPBOT1
                TOXB(L,K,NT)   = TMPVAL*TMPBOT2
              enddo  
            endif  
            if( ISTRAN(6) >= 1 )then  
              TMPVAL = SEDBT(L,K)
              SEDBT(L,K+1) = TMPVAL*TMPBOT1
              SEDBT(L,K)   = TMPVAL*TMPBOT2
              do NS = 1,NSED  
                TMPVAL = SEDB(L,K,NS)
                SEDB(L,K+1,NS) = TMPVAL*TMPBOT1
                SEDB(L,K,NS)   = TMPVAL*TMPBOT2
                STFPOCB(L,K+1,NS) = STFPOCB(L,K,NS)  
                VFRBED(L,K+1,NS) = VFRBED(L,K,NS)  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              TMPVAL = SNDBT(L,K)
              SNDBT(L,K+1) = TMPVAL*TMPBOT1
              SNDBT(L,K)   = TMPVAL*TMPBOT2
              do NS = 1,NSND  
                NX = NS + NSED  
                TMPVAL = SNDB(L,K,NS)
                SNDB(L,K+1,NS) = TMPVAL*TMPBOT1
                SNDB(L,K,NS)   = TMPVAL*TMPBOT2
                STFPOCB(L,K+1,NX) = STFPOCB(L,K,NX)  
                VFRBED(L,K+1,NX) = VFRBED(L,K,NX)  
              enddo  
            endif
            KBT(L) = 2
          endif
        endif  
      enddo  
  
      ! *** COMBINE THIN PARENT LAYERS WITH LAYER BELOW TO FORM NEW PARENT  
      do LP = LF,LL  
        L = LSED(LP)
        K = KBT(L) - 1
        if( K == 0 ) CYCLE
        if( HBED(L,K) < HBEDMIN )then  ! *** LAYER BELOW TOP/ACTIVE LAYER
          if( KBT(L) > 2 )then  
            if( HBED(L,K) < 1E-6 )then         ! *** 2021-12-09
              HBED(L,K)   = 0.0       
              VDRBED(L,K) = 0.0  
              PORBED(L,K) = 0.0  
              STDOCB(L,K) = 0.0  
              STPOCB(L,K) = 0.0  
            endif              
            TMPBOT1 = HBED(L,KBT(L)-2)/(1. + VDRBED(L,KBT(L)-2))    ! *** 2 LAYERS BELOW TOP/ACTIVE LAYER
            TMPBOT2 = HBED(L,KBT(L)-1)/(1. + VDRBED(L,KBT(L)-1))    ! *** 1 LAYERS BELOW TOP/ACTIVE LAYER
            TMPTOP1 = VDRBED(L,KBT(L)-2)*TMPBOT1  
            TMPTOP2 = VDRBED(L,KBT(L)-1)*TMPBOT2  
            VDRBED(L,KBT(L)-2) = (TMPTOP1 + TMPTOP2)/(TMPBOT1 + TMPBOT2)  
            PORBED(L,KBT(L)-2) = VDRBED(L,KBT(L)-2)/(1. + VDRBED(L,KBT(L)-2))  
            HBED(L,KBT(L)-2) = HBED(L,KBT(L)-2) + HBED(L,KBT(L)-1)  
            HBED(L,KBT(L)-1) = HBED(L,KBT(L))  
            VDRBED(L,KBT(L)-1) = VDRBED(L,KBT(L))  
            PORBED(L,KBT(L)-1) = PORBED(L,KBT(L))  
            HBED(L,KBT(L)) = 0.  
            VDRBED(L,KBT(L)) = 0.  
            PORBED(L,KBT(L)) = 0.  
            STDOCB(L,KBT(L)) = 0.  
            STPOCB(L,KBT(L)) = 0.  
            BDENBED(L,KBT(L)) = 0.
            if( ISTRAN(6) >= 1 )then  
              SEDBT(L,KBT(L)-2) = 0.  
              SEDBT(L,KBT(L)-1) = 0.  
              do NS = 1,NSED  
                SEDB(L,KBT(L)-2,NS) = SEDB(L,KBT(L)-2,NS) + SEDB(L,KBT(L)-1,NS)  
                SEDBT(L,KBT(L)-2) = SEDBT(L,KBT(L)-2)     + SEDB(L,KBT(L)-2,NS)  
                SEDB(L,KBT(L)-1,NS) = SEDB(L,KBT(L),NS)  
                SEDBT(L,KBT(L)-1) = SEDBT(L,KBT(L)-1)     + SEDB(L,KBT(L)-1,NS)  
                SEDB(L,KBT(L),NS) = 0.0  
                SEDBT(L,KBT(L)) = 0.0  
                STFPOCB(L,KBT(L),NS) = 0.0  
                VFRBED(L,KBT(L),NS) = 0.0  
              enddo  
            endif  
            if( ISTRAN(7) >= 1 )then  
              SNDBT(L,KBT(L)-2) = 0.  
              SNDBT(L,KBT(L)-1) = 0.  
              do NS = 1,NSND  
                SNDB(L,KBT(L)-2,NS) = SNDB(L,KBT(L)-2,NS) + SNDB(L,KBT(L)-1,NS)  
                SNDBT(L,KBT(L)-2)   = SNDBT(L,KBT(L)-2)   + SNDB(L,KBT(L)-2,NS)  
                SNDB(L,KBT(L)-1,NS) = SNDB(L,KBT(L),NS)  
                SNDBT(L,KBT(L)-1)   = SNDBT(L,KBT(L)-1)   + SNDB(L,KBT(L)-1,NS)  
                SNDB(L,KBT(L),NS)   = 0.0  
                SNDBT(L,KBT(L))     = 0.0  
                STFPOCB(L,KBT(L),NS+NSED) = 0.0  
                VFRBED(L,KBT(L),NS+NSED) = 0.0  
              enddo  
            endif  
            if( ISTRAN(5) >= 1 )then  
              do NT = 1,NTOX  
                TOXB(L,KBT(L)-2,NT) = TOXB(L,KBT(L)-2,NT) + TOXB(L,KBT(L)-1,NT)  
                TOXB(L,KBT(L)-1,NT) = TOXB(L,KBT(L),NT)  
                TOXB(L,KBT(L),NT)   = 0.0  
              enddo  
            endif  
            KBT(L) = KBT(L)-1  

          elseif( K == 1 )then  ! *** LAYER BELOW TOP/ACTIVE LAYER IS GETTING THIN
            TMPVAL = HBED(L,2) + HBED(L,1)
            if( TMPVAL < HBEDAL )then
              ! *** COLLAPSE THE ACTIVE LAYER BY REMOVING THE PARENT LAYER
              HBED(L,1) = TMPVAL
              HBED(L,2) = 0.
              VDRBED(L,1) = VDRBED(L,2)  
              PORBED(L,1) = PORBED(L,2)  
              if( ISTRAN(6) >= 1 )then  
                SEDBT(L,1) = 0.  
                SEDBT(L,2) = 0.  
                do NS = 1,NSED  
                  SEDB(L,1,NS) = SEDB(L,1,NS) + SEDB(L,2,NS)
                  SEDBT(L,1)   = SEDBT(L,1)   + SEDB(L,1,NS)  
                  SEDB(L,2,NS) = 0. 
                enddo  
              endif  
              if( ISTRAN(7) >= 1 )then  
                SNDBT(L,1) = 0.  
                SNDBT(L,2) = 0.  
                do NS = 1,NSND  
                  SNDB(L,1,NS) = SNDB(L,1,NS) + SNDB(L,2,NS)  
                  SNDBT(L,1)   = SNDBT(L,1)   + SNDB(L,1,NS)  
                  SNDB(L,2,NS) = 0.
                enddo  
              endif  
              if( ISTRAN(5) >= 1 )then  
                do NT = 1,NTOX  
                  TOXB(L,1,NT) = TOXB(L,1,NT) + TOXB(L,2,NT)  
                  TOXB(L,2,NT) = 0.
                enddo  
              endif  

              ! *** ZERO OLD TOP LAYER
              VDRBED(L,KBT(L)) = 0.0  
              PORBED(L,KBT(L)) = 0.0  
              STDOCB(L,KBT(L)) = 0.0  
              STPOCB(L,KBT(L)) = 0.0  

              KBT(L) = KBT(L)-1  
              
            endif
            
          elseif( HBED(L,KBT(L)) < 0.0 )then  
            ! *** ZERO NEGATIVE THICKNESSES
            HBED(L,KBT(L)) = 0.0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  
          endif  
          ! *** PMC END BLOCK  
        endif  
      enddo  

      ! *** UPDATE BULK DENSITY  
      do LP = LF,LL  
        L = LSED(LP)
        ! ***REMOVED KB LOOP, ONLY COMPUTE THE TOP LAYER.  CONSOLIDATION IS HANDLED LATER IN CALBED
        K = KBT(L)
        if( HBED(L,K) > 0. )then  
          ! *** UPDATE TOTAL/WET DENSITY
          BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)  
        else  
          BDENBED(L,K) = 0.  
        endif  
      enddo  

    enddo  ! *** END OF DOMAIN LOOP
    !  $OMP END PARALLEL DO

  endif    ! *** END OF ARMORING OPTION SECTION

  ! ***************************************************************************
  return 

END  

