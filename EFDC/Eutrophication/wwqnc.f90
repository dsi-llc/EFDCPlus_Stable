! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE WWQNC

  ! CHANGE RECORD
  ! WRITE INFORMATION OF NEGATIVE WQ STATE VARIABLES (UNIT IWQONC).

  use GLOBAL
  use Variables_WQ
  
  implicit none

  integer     :: L, K, NW, LP, ND, IFLAG
  character*5 :: WQVN(23)

  DATA WQVN/ &
      'BC ','BD ','BG ','RPOC','LPOC','DOC ','RPOP','LPOP','DOP ','PO4T','RPON','LPON', &
      'DON ','NH4 ','NO3 ','SU ','SA   ','COD  ','O2   ','TAM  ','FCB  ', 'CO2 ','MALG '/

  open(1,FILE = OUTDIR//'WQ3DNC.LOG',STATUS = 'UNKNOWN',POSITION = 'APPEND')

  IFLAG = 0
  do L = 2,LA
    do K = 1,KC
      do NW = 1,NWQV
        if( WQV(L,K,NW) < 0.0 )then
          write(1,90) WQVN(NW),ITNWQ,L,IL(L),JL(L),K,WQV(L,K,NW)
          IFLAG = 1
        endif
      enddo
    enddo
  enddo
  close(1)
  90 FORMAT(A5, I8, 4I5, E11.3)
  
  ! *** ZERO NEGATIVE CONCENTRATIONS
  if( IWQNC > 1 .and. IFLAG == 1 )then
    !$OMP PARALLEL DO PRIVATE(ND,K,LP,L,NW)
    do ND = 1,NDM  
      do K = 1,KC
        do LP = 1,LLWET(K,ND)
          L = LKWET(LP,K,ND) 
          do NW = 1,NWQV
            if( ISKINETICS(NW) > 0 )then
              if( WQV(L,K,NW) < 0.0 ) WQV(L,K,NW) = 0.0
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif
  
  return
  
END

