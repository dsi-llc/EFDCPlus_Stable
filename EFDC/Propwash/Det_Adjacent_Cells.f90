! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!< @details Determines adjacent cells for a given cell
!< @author Zander Mausolff
!< @param[in]    
!< @param[out]  
!---------------------------------------------------------------------------!  
subroutine Det_Adjacent_Cells

  Use GLOBAL
  Use Variables_Propwash
    
  implicit none
    
  ! *** Dummy variables
    
  ! *** Local variables
  Integer :: L
    
  Allocate(adjacent_l(9,0:LA))
    
  adjacent_l = 0
    
  do L=2,LA
    ! *** ORDER OF CELLS
    ! ***   1  2  3
    ! ***   4  5  6
    ! ***   7  8  9
    
    ! *** Set center
    adjacent_l(5,L) = L
      
    IF( SUBO(L)     +SVBO(LNWC(L)) > 1.5 .OR. SVBO(LNC(L))+SUBO(LNC(L))  > 1.5 ) adjacent_l(1,L) = LNWC(L)
    IF( SVBO(LNC(L))+SUBO(LNEC(L)) > 1.5 .OR. SUBO(LEC(L))+SVBO(LNEC(L)) > 1.5 ) adjacent_l(3,L) = LNEC(L)
    IF( SUBO(L)     +SVBO(LWC(L))  > 1.5 .OR. SVBO(L)     +SUBO(LSC(L))  > 1.5 ) adjacent_l(7,L) = LSWC(L)
    IF( SVBO(L)     +SUBO(LSEC(L)) > 1.5 .OR. SUBO(LEC(L))+SVBO(LEC(L))  > 1.5 ) adjacent_l(9,L) = LSEC(L)
      
    IF( SVBO(LNC(L)) > 0.5 ) adjacent_l(2,L) = LNC(L)
    IF( SUBO(L) > 0.5 )      adjacent_l(4,L) = LWC(L)
    IF( SUBO(LEC(L)) > 0.5 ) adjacent_l(6,L) = LEC(L)
    IF( SVBO(L) > 0.5 )      adjacent_l(8,L) = LSC(L)
  end do
    
end subroutine Det_Adjacent_Cells