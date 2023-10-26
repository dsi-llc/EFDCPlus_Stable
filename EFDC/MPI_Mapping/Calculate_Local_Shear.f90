! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

  Use Global
  Use MPI
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out

  Implicit None

  ! *** Local variables
  Integer :: L

  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN

    IF( LSEDZLJ )THEN
      ! *** TAU(LCM) Shear Stress in dynes/cm^2,  WATERDENS  gm/cm3
      ! *** SHEAR IS SAVED AS N/M^2 (Pascals) normalized to water density  (M2/S2)
      DO L=2,LA
        IF( LBED(L) )THEN
          SHEAR_Local(L) = QQ(L,0)/CTURB2
        ELSE
          SHEAR_Local(L) = TAU(L) * 0.1 / (WATERDENS*1000.)
        ENDIF
      ENDDO

    ELSEIF( ISBEDSTR >= 1 )THEN
      DO L=2,LA
        IF( LBED(L) )THEN
          SHEAR_Local(L) = QQ(L,0)/CTURB2
        ELSE
          SHEAR_Local(L) = TAUBSED(L)
        ENDIF
      ENDDO

      IF( ISBEDSTR == 1 )THEN
        DO L=2,LA
          IF( LBED(L) )THEN
            SHEAR_Local2(L) = QQ(L,0)/CTURB2
          ELSE
            SHEAR_Local2(L) = TAUBSND(L)
          ENDIF
        ENDDO
      ENDIF
    ELSE
      ! *** TOTAL BED SHEAR STRESS
      DO L=2,LA
        IF( LBED(L) )THEN
          SHEAR_Local(L) = QQ(L,0)/CTURB2
        ELSE
          SHEAR_Local(L) = TAUB(L)
        ENDIF
      ENDDO
    ENDIF

  ELSE
    ! *** TOTAL BED SHEAR STRESS
    DO L=2,LA
      SHEAR_Local(L) = MAX(QQ(L,0),QQMIN)/CTURB2
    ENDDO
  ENDIF

End Subroutine Calculate_Local_Shear
