! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE PARTMIX(NT)

  ! *** *******************************************************************C
  !
  ! ***  SUBROUTINE PARTMIX CALCULATES PARTICLE MIXING OF TOXICS IN THE 
  ! ***  TOP LAYER(S) [PMXDEPTH] OF THE SEDIMENT BED
  !
  ! *** *******************************************************************C

  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-12           Paul M. Craig      ADDED OMP AND CLEANED UP CODE

  use GLOBAL

  implicit none
  
  integer :: ND,LF,LL,LP,K,L,NT,LZ,KM,NPMXPTS,NDP,KK,NP
  real    :: DEPINBED,WT1,WT2,TMPVAL,DELHBED,TERM1,TERM2,TERM3
  real    :: PARTMIXAVG(KBM)

  ! *** *******************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
  do ND = 1,NDM  
    LF = (ND-1)*LDMSED+1  
    LL = MIN(LF+LDMSED-1,LASED)

    do K = 1,KB
      do LP = LF,LL
        L = LSED(LP)  
        PARTMIXZ(L,K) = 0.0
      enddo
    enddo
  enddo   ! *** END OF DOMAIN 
  !$OMP END DO
  
  !----------------------------------------------------------------------
  !  OLD PARTICLE MIXING IS NOW OPTION ISPMXZ = 2
  if( ISPMXZ(NT) == 2 )then

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LZ,K,KM,NP,NDP,DEPINBED,WT1,WT2,TMPVAL)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        DEPINBED = 0.
        LZ = LPMXZ(L)
        
        if( KBT(L) > 2 )then
          do K = KBT(L),2,-1
            KM = K-1
            DEPINBED = DEPINBED+HBED(L,K)
            do NP = 1,NPMXPTS-1
              NDP = NP+1
              if( DEPINBED >= PMXDEPTH(NP,LZ) .and. DEPINBED < PMXDEPTH(NDP,LZ) )then
                WT1 = DEPINBED-PMXDEPTH(NP,LZ)
                WT2 = PMXDEPTH(NDP,LZ)-DEPINBED
                TMPVAL = PMXDEPTH(NP+1,LZ)-PMXDEPTH(NP,LZ)
                PARTMIXZ(L,KM) = (WT2*PMXCOEF(NP,LZ) + WT1*PMXCOEF(NDP,LZ))/TMPVAL
              endif
            enddo
          enddo
        else
          K  = KBT(L)
          KM = K-1
          DEPINBED = DEPINBED + HBED(L,K)
          do NP = 1,NPMXPTS-1
            NDP = NP+1
            if( DEPINBED >= PMXDEPTH(NP,LZ) .and. DEPINBED < PMXDEPTH(NDP,LZ) )then
              WT1 = DEPINBED-PMXDEPTH(NP,LZ)
              WT2 = PMXDEPTH(NDP,LZ)-DEPINBED
              TMPVAL = PMXDEPTH(NP+1,LZ)-PMXDEPTH(NP,LZ)
              PARTMIXZ(L,KM) = (WT2*PMXCOEF(NP,LZ) + WT1*PMXCOEF(NDP,LZ))/TMPVAL
            endif
          enddo
        endif
      enddo

      do LP = LF,LL
        L = LSED(LP)  
        do K = KBT(L),2,-1
          KM = K-1
          PARTMIXZ(L,KM) = 2.*PARTMIXZ(L,KM)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT)+1.E-12)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO
  
  endif

  !----------------------------------------------------------------------
  !  NEW PARTICLE MIXING IS NOW OPTION ISPMXZ = 1
  if( ISPMXZ(NT) == 1 )then

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LZ,K,KK,KM,NP,NDP,DEPINBED,DELHBED,WT1,WT2,TMPVAL,TERM1,TERM2,TERM3)
    do ND = 1,NDM  
      LF = (ND-1)*LDMSED+1  
      LL = MIN(LF+LDMSED-1,LASED)

      do LP = LF,LL
        L = LSED(LP)  
        LZ = LPMXZ(L)

        do K = 1,KB
          PARTMIXAVG(K) = 0.0
          PARTMIXZ(L,K) = 0.0
        enddo

        DEPINBED = 0.
        do K = KBT(L),1,-1
          DELHBED = HBED(L,K)/10.
          do KK = 1,10
            DEPINBED = DEPINBED+DELHBED
            do NP = 1,NPMXPTS-1
              NDP = NP+1
              if( DEPINBED >= PMXDEPTH(NP,LZ) .and. DEPINBED < PMXDEPTH(NDP,LZ) )then
                WT1 = DEPINBED-PMXDEPTH(NP,LZ)
                WT2 = PMXDEPTH(NDP,LZ)-DEPINBED
                TMPVAL = PMXDEPTH(NP+1,LZ)-PMXDEPTH(NP,LZ)
                PARTMIXAVG(K) = PARTMIXAVG(K) + (WT2*PMXCOEF(NP,LZ)+WT1*PMXCOEF(NDP,LZ))/TMPVAL
              endif
            enddo
          enddo
          PARTMIXAVG(K) = PARTMIXAVG(K)/10.
        enddo

        do K = 1,KBT(L)-1
          TERM1 = HBED(L,K)/PARTMIXAVG(K)
          TERM2 = HBED(L,K+1)/PARTMIXAVG(K+1)
          TERM3 = HBED(L,K)+HBED(L,K+1)
          PARTMIXZ(L,K) = TERM3/(TERM1+TERM2)
        enddo

      enddo

      do LP = LF,LL
        L = LSED(LP)  
        do K = KBT(L),2,-1
          KM = K-1
          PARTMIXZ(L,KM) = 2.*PARTMIXZ(L,KM)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT)+1.E-12)
        enddo
      enddo
    enddo   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  endif

  ! *** *******************************************************************C\
  !$OMP END PARALLEL
  
  return

END

