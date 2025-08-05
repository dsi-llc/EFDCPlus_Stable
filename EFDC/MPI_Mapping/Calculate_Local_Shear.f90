! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! @details Calculate local shear to avoid some calculations during the write
!! out of EE
! @author Zander Mausolff
! @date 1/9/2019
!---------------------------------------------------------------------------!

Subroutine Calculate_Local_Shear

  use GLOBAL
  use MPI
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out

  implicit none

  ! *** Local variables
  integer :: L

  if( ISTRAN(6) > 0 .or. ISTRAN(7) > 0 )then

    if( LSEDZLJ )then
      ! *** TAU(LCM) Shear Stress in dynes/cm^2,  WATERDENS  gm/cm3
      ! *** SHEAR IS SAVED AS N/M^2 (Pascals) normalized to water density  (M2/S2)
      do L = 2,LA
        if( LBED(L) )then
          SHEAR_Local(L) = QQ(L,0)/CTURB2
        else
          SHEAR_Local(L) = TAU(L) * 0.1 / (WATERDENS*1000.)
        endif
      enddo

    elseif( ISBEDSTR >= 1 )then
      do L = 2,LA
        if( LBED(L) )then
          SHEAR_Local(L) = QQ(L,0)/CTURB2
        else
          SHEAR_Local(L) = TAUBSED(L)
        endif
      enddo

      if( ISBEDSTR == 1 )then
        do L = 2,LA
          if( LBED(L) )then
            SHEAR_Local2(L) = QQ(L,0)/CTURB2
          else
            SHEAR_Local2(L) = TAUBSND(L)
          endif
        enddo
      endif
    else
      ! *** TOTAL BED SHEAR STRESS
      do L = 2,LA
        if( LBED(L) )then
          SHEAR_Local(L) = QQ(L,0)/CTURB2
        else
          SHEAR_Local(L) = TAUB(L)
        endif
      enddo
    endif

  else
    ! *** TOTAL BED SHEAR STRESS
    if( ISGOTM > 0 )then
      do L = 2,LA
        SHEAR_Local(L) = MAX(QQ(L,KSZ(L)-1),QQMIN)/CTURB2
      enddo
    else    
      do L = 2,LA
        SHEAR_Local(L) = MAX(QQ(L,0),QQMIN)/CTURB2
      enddo
    endif
  endif

End Subroutine Calculate_Local_Shear
