! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  ! Module: Mod_Map_Write_EE_Binary
  !
  !> @details
  !
  !> @author Zander Mausolff
  !> @date 1/8/2019
  !> @Updated in 2020 and 2021 by Paul M. Craig and Duc Kien TRAN
  !---------------------------------------------------------------------------!
  Module Mod_Map_Write_EE_Binary

  Use GLOBAL
  USE WQ_DIAGENESIS
  Use WQ_RPEM_MODULE
  Use Variables_WQ
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  Use Variables_MPI_MapGatherSort
  Use Mod_Assign_Loc_Glob_For_Write
  USE FIELDS

  Implicit none

  ! *** Local variables
  Integer :: num_arrays_to_write_out

  ! *** End variable definitions

  Contains
  !---------------------------------------------------------------------------!
  ! Subroutine: Map_Write_EE_Binary
  !
  !> @details
  !
  !> @author Zander Mausolff
  !---------------------------------------------------------------------------!
  ! *** Driver subroutine
  Subroutine Map_Write_EE_Binary

  Implicit None
  ! *** Dummy variables

  ! *** Local variables

  ! *** Always call these for every EFDC+ run
  Call Map_Write_WSOUT
  Call Map_Write_VELOUT
  Call Map_Write_BCOUT

  ! *** Conditionally write out the following - logic from EE_LINKAGE subroutine
  IF( ISSPH(8) >= 1 ) Then
    Call Map_Write_WCOUT
  End if

  IF( LSEDZLJ )THEN
    CALL Map_Write_SEDZLJOUT
  ELSE
    IF( ISBEXP >= 1 .AND. KB > 1 )THEN
      IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) CALL Map_Write_BEDOUT
    ENDIF
  ENDIF

  IF( ISTRAN(8) > 0 )THEN
    CALL Map_Write_WQOUT
    IF( IWQBEN > 0 .AND. ISSDBIN /= 0 ) CALL Map_Write_SDOUT
    IF( ISRPEM > 0) CALL Map_Write_RPEMOUT
    !IF( ISFFARM > 0) CALL SHELLFISHOUT() ! @todo shellfish update
  ENDIF

  End Subroutine Map_Write_EE_Binary

  !---------------------------------------------------------------------------!
  ! Subroutine: Map_Write_WSOUT
  !
  !> @details
  !
  !> @author Zander Mausolff
  !> @date 1/8/2019
  !---------------------------------------------------------------------------!
  ! *** WSOUT - EE_WS.OUT
  Subroutine Map_Write_WSOUT

  Implicit None

  ! *** Local variables
  Integer :: i,j

  j=0

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(HP,1), HP, size(HP_Global,1), HP_Global)

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  End Subroutine Map_Write_WSOUT
  !---------------------------------------------------------------------------!
  ! Subroutine: Map_Write_VELOUT
  !
  !> @details
  !
  !> @author Zander Mausolff
  !> @date 1/8/2019
  !---------------------------------------------------------------------------!
  ! *** VELOUT
  Subroutine Map_Write_VELOUT

  Implicit None

  ! *** Local variables
  Integer :: i, j

  ! *** Clear out the pointer array
  j = 0

  ! *** Point to the first array
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(U,1), size(U,2), U, size(U_Global,1), size(U_Global,2), U_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(V,1), size(V,2), V, size(V_Global,1), size(V_Global,2), V_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(W,1), size(W,2), W, size(W_Global,1), size(W_Global,2), W_Global)

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  End Subroutine Map_Write_VELOUT
  !---------------------------------------------------------------------------!
  ! Subroutine: Map_Write_BCOUT
  !
  !> @details
  !
  !> @author Zander Mausolff
  !> @date 1/8/2019
  !---------------------------------------------------------------------------!
  ! *** BCOUT
  Subroutine Map_Write_BCOUT

  Implicit None

  ! *** Local variables
  Integer :: i, j

  num_arrays_to_write_out = 3
  j = 0

  ! *** Depends on if we have more than one vertical layer
  if( KC > 1 )then

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(QSUM,1), size(QSUM,2), QSUM, &
      size(QSUM_Global,1), size(QSUM_Global,2), QSUM_Global)
  else
    ! *** 1D arrays when KC < 1
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(QSUME,1), QSUME, size(QSUME_Global,1), QSUME_Global)
  endif

  ! *** Always write out these to EE binaries
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(VHDX2,1), size(VHDX2,2), VHDX2, &
    size(VHDX2_Global,1), size(VHDX2_Global,2), VHDX2_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(UHDY2,1), size(UHDY2,2), UHDY2, &
    size(UHDY2_Global,1), size(UHDY2_Global,2), UHDY2_Global)

  ! @todo do HSCTL and WRCTL need to be modified at all?

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)


  End Subroutine Map_Write_BCOUT
  !---------------------------------------------------------------------------!
  ! Subroutine: Map_Write_WCOUT
  !
  !> @details This routine maps out the variables necessary for writing out
  ! EE_WC.OUT Follows the logic in WCOUT for determining what variables to write out
  !
  !> @author Zander Mausolff
  !> @date 1/8/2019
  !---------------------------------------------------------------------------!
  ! *** WCOUT
  Subroutine Map_Write_WCOUT

  Implicit None

  ! *** Local variables
  Integer :: i,j,m

  INTEGER,TARGET,ALLOCATABLE,DIMENSION(:) :: LWC_temp_gl
  INTEGER,TARGET,ALLOCATABLE,DIMENSION(:) :: LEC_temp_gl
  INTEGER,TARGET,ALLOCATABLE,DIMENSION(:) :: LSC_temp_gl
  INTEGER,TARGET,ALLOCATABLE,DIMENSION(:) :: LNC_temp_gl
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WVHEIGHT_temp_gl
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WVFREQ_temp_gl
  REAL,TARGET,ALLOCATABLE,DIMENSION(:) :: WVDIR_temp_gl

  Allocate(LWC_temp_gl(LCM))
  Allocate(LEC_temp_gl(LCM))
  Allocate(LSC_temp_gl(LCM))
  Allocate(LNC_temp_gl(LCM))
  Allocate(WVHEIGHT_temp_gl(LCM))
  Allocate(WVFREQ_temp_gl(LCM))
  Allocate(WVDIR_temp_gl(LCM))

  LEC_temp_gl      = 0
  LNC_temp_gl      = 0
  LSC_temp_gl      = 0
  WVHEIGHT_temp_gl = 0
  WVFREQ_temp_gl   = 0
  WVDIR_temp_gl    = 0

  j = 0

  !! *** WRITE THE TOP LAYER INDEX
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(KBT,1), KBT, size(KBT_Global,1), KBT_Global)
  ENDIF

  ! *** Needed for now to determine Shear_Local.  Avoids some of the issues around accumulation
  Call Calculate_Local_Shear

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Shear_Local,1), Shear_Local, size(Shear_Global,1), Shear_Global)
  IF( ISBEDSTR == 1 .AND. .NOT. LSEDZLJ ) THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(Shear_Local2,1), Shear_Local2, size(Shear_Global2,1), Shear_Global2)
  ENDIF

  IF( ISWAVE >= 1 )THEN
    j = j + 1

    Call Assign_Loc_Glob_For_Write(j, size(QQWV3,1), QQWV3, size(QQWV3_Global,1), QQWV3_Global)

    ! *** Shear due to Current Only
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TBY, 1), TBY, size(TBY_Global, 1), TBY_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TBX, 1), TBX, size(TBX_Global, 1), TBX_Global)

  End if

  !*** Need global versions of these for the shear calculation which occurs during the writing out
  j = j + 1
  ! *** need to convert local L to global L
  do m = 1, LA
    LWC_temp_gl(m) = Map2Global(LWC(m)).LG
  end do

  Call Assign_Loc_Glob_For_Write(j, size(LWC_temp_gl,1), LWC_temp_gl, size(LWC_Global, 1), LWC_Global)

  j = j + 1
  ! *** need to convert local L to global L
  do m = 1, LA
    LEC_temp_gl(m) = Map2Global(LEC(m)).LG
  end do

  Call Assign_Loc_Glob_For_Write(j, size(LEC_temp_gl,1), LEC_temp_gl, size(LEC_Global, 1), LEC_Global)

  j = j + 1
  ! *** need to convert local L to global L
  do m = 1, LA
    LSC_temp_gl(m) = Map2Global(LSC(m)).LG
  end do

  Call Assign_Loc_Glob_For_Write(j, size(LSC_temp_gl,1), LSC_temp_gl, size(LSC_Global, 1), LSC_Global)

  j = j + 1
  ! *** need to convert local L to global L
  do m = 1, LA
    LNC_temp_gl(m) = Map2Global(LNC(m)).LG
  end do

  Call Assign_Loc_Glob_For_Write(j, size(LNC_temp_gl,1), LNC_temp_gl, size(LNC_Global, 1), LNC_Global)

  ! @todo develop interface for handling these types
  IF( ISWAVE >= 3 )THEN
    j = j + 1
    do m = 1, LA
      WVHEIGHT_temp_gl(m) = WV(m).HEIGHT
    end do

    Call Assign_Loc_Glob_For_Write(j, size(WVHEIGHT_temp_gl,1), WVHEIGHT_temp_gl, size(WV_HEIGHT_Global,1), WV_HEIGHT_Global)

    j = j + 1
    do m = 1, LA
      WVFREQ_temp_gl(m) = WV(m).FREQ
    end do
    Call Assign_Loc_Glob_For_Write(j, size(WVFREQ_temp_gl,1), WVFREQ_temp_gl, size(WV_FREQ_Global,1), WV_FREQ_Global)

    j = j + 1
    do m = 1, LA
      WVDIR_temp_gl(m) = WV(m).DIR
    end do
    Call Assign_Loc_Glob_For_Write(j, size(WVDIR_temp_gl,1), WVDIR_temp_gl, size(WV_DIR_Global, 1), WV_DIR_Global)

  ENDIF

  IF( ISWAVE == 4 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(WV.DISSIPA(KC),1), WV.DISSIPA(KC), &
      size(WV_Global.DISSIPA(KC),1), WV_Global.DISSIPA(KC))

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(WVHUU,1), size(WVHUU,2), WVHUU, &
      size(WVHUU_Global,1), size(WVHUU_Global,2), WVHUU_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(WVHVV,1), size(WVHVV,2), WVHVV, &
      size(WVHVV_Global,1), size(WVHVV_Global,2), WVHVV_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(WVHUV,1), size(WVHUV,2), WVHUV, &
      size(WVHUV_Global,1), size(WVHUV_Global,2), WVHUV_Global)

  ENDIF

  ! *** Get KSZ_Global when using SigmaZ Uniform option - DKT
  IF (IGRIDV == 2) THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(KSZ,1), KSZ, size(KSZ_Global), KSZ_Global)
  ENDIF

  IF( ISTRAN(1) >= 1 ) Then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SAL,1), size(SAL,2), SAL, &
      size(SAL_Global,1), size(SAL_Global,2), SAL_Global)
  End if

  IF( ISTRAN(2) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TEM,1), size(TEM,2), TEM, &
      size(TEM_Global,1), size(TEM_Global,2), TEM_Global)

    IF( TEMBO > 0.) Then
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(TEMB,1), TEMB, size(TEMB_Global), TEMB_Global)
    End if

    IF( IEVAP > 1 )THEN
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(EVAPT,1), EVAPT, size(EVAPT_Global), EVAPT_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(RAINT,1), RAINT, size(RAINT_Global), RAINT_Global)
    ENDIF
  ENDIF

  IF( ISTRAN(3) >= 1 ) Then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(DYE,1), size(DYE,2), size(DYE,3), DYE, &
      size(DYE_Global,1), size(DYE_Global,2), size(DYE_Global,3), DYE_Global)
  Endif

  IF( ISTRAN(4) >= 1 ) Then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SFL,1), size(SFL,2), SFL, &
      size(SFL_Global,1), size(SFL_Global,2), SFL_Global)
  Endif

  IF( ISTRAN(5) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOXB,1), size(TOXB,2), size(TOXB,3), TOXB,&
      size(TOXB_Global,1), size(TOXB_Global,2), size(TOXB_Global,3), TOXB_Global)
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOX,1), size(TOX,2), size(TOX,3), TOX,&
      size(TOX_Global,1), size(TOX_Global,2), size(TOX_Global,3), TOX_Global)
  ENDIF

  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(BELV,1), BELV, size(BELV_Global,1), BELV_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(HBED,1), size(HBED,2), HBED, &
      size(HBED_Global,1), size(HBED_GLobal,2), HBED_Global)
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(BDENBED,1), size(BDENBED,2), BDENBED, &
      size(BDENBED_Global,1), size(BDENBED_Global,2), BDENBED_Global)
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(PORBED,1), size(PORBED,2), PORBED, &
      size(PORBED_Global,1), size(PORBED_GLobal,2), PORBED_Global)

    IF( ISTRAN(6) >= 1 ) Then
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SEDB,1), size(SEDB,2), size(SEDB,3), SEDB,&
        size(SEDB_Global,1), size(SEDB_Global,2), size(SEDB_Global,3), SEDB_Global)

      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SED,1), size(SED,2), size(SED,3), SED,&
        size(SED_Global,1), size(SED_Global,2), size(SED_Global,3), SED_Global)
    End if

    IF( ISTRAN(7) >= 1 ) Then
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SNDB,1), size(SNDB,2), size(SNDB,3), SNDB,&
        size(SNDB_Global,1), size(SNDB_Global,2), size(SNDB_Global,3), SNDB_Global)
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(SND,1), size(SND,2), size(SND,3), SND,&
        size(SND_Global,1), size(SND_Global,2), size(SND_Global,3), SND_Global)

      IF( ICALC_BL > 0 .AND. NSND > 0 )THEN
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(QSBDLDX,1), size(QSBDLDX,2), QSBDLDX, &
          size(QSBDLDX_Global,1), size(QSBDLDX_Global,2), QSBDLDX_Global)

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(QSBDLDY,1), size(QSBDLDY,2), QSBDLDY, &
          size(QSBDLDY_Global,1), size(QSBDLDY_Global,2), QSBDLDY_Global)
      ENDIF
    Endif

    !j = j + 1
    !Call Assign_Loc_Glob_For_Write(j, size(VDRBED,1), size(VDRBED,2), VDRBED, &
    !                       size(VDRBED_Global,1), size(VDRBED_Global,2), VDRBED_Global)
  ELSEIF (BATHY.IFLAG > 0 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(BELV,1), BELV, size(BELV_Global,1), BELV_Global)
  ENDIF

  IF( ISGWIE > 0 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(EVAPSW,1), EVAPSW, size(EVAPSW_Global,1), EVAPSW_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(EVAPGW,1), EVAPGW, size(EVAPGW_Global,1), EVAPGW_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(QGW,1), QGW, size(QGW_Global,1), QGW_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(AGWELV,1), AGWELV, size(AGWELV_Global,1), AGWELV_Global)

  ENDIF

  IF( ISTRAN(2) > 0 .AND. ISICE  >= 3 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(ICETHICK,1), ICETHICK, size(ICETHICK_Global,1), ICETHICK_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(ICETEMP,1), ICETEMP, size(ICETEMP_Global,1), ICETEMP_Global)
  ENDIF

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  deallocate(LEC_temp_gl)
  deallocate(LNC_temp_gl)
  deallocate(WVHEIGHT_temp_gl)
  deallocate(WVFREQ_temp_gl)
  deallocate(WVDIR_temp_gl)

  End Subroutine Map_Write_WCOUT

  ! *** SEDZLJOUT
  Subroutine Map_Write_SEDZLJOUT

  Implicit None

  ! *** Local
  Integer :: i, j,l,k,nn
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_TSED
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_BULKDENS
  Real(rkd), Target, Allocatable, Dimension(:,:,:) :: Reverse_Temp_3D

  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_TSED
  Real(rkd), Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_BULKDENS
  Real(rkd), Target, Allocatable, Dimension(:,:,:) :: Gl_Reverse_Temp_3D


  Allocate(Reverse_Temp_2D_TSED(LCM,KB))
  Allocate(Reverse_Temp_2D_BULKDENS(LCM,KB))
  Allocate(Reverse_Temp_3D(LCM,KB,NSCM))

  Allocate(Gl_Reverse_Temp_2D_TSED(LCM_Global,KB))
  Allocate(Gl_Reverse_Temp_2D_BULKDENS(LCM_Global,KB))
  Allocate(GL_Reverse_Temp_3D(LCM_Global,KB,NSCM))


  Reverse_Temp_2D_TSED        = 0.0
  Reverse_Temp_2D_BULKDENS    = 0.0
  Reverse_Temp_3D             = 0.0
  Gl_Reverse_Temp_2D_TSED     = 0.0
  Gl_Reverse_Temp_2D_BULKDENS = 0.0
  Gl_Reverse_Temp_3D          = 0.0

  j = 0

  ! *** REAL*4
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(TAU,1), TAU, size(TAU_Global,1), TAU_Global)

  j = j + 1

  ! *** Reverse Order of Array
  DO L=2,LA
    DO K=1,KB
      Reverse_Temp_2D_TSED(L,K) = TSED(K,L)
    ENDDO
  ENDDO

  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_TSED,1), size(Reverse_Temp_2D_TSED,2), Reverse_Temp_2D_TSED, &
    size(Gl_Reverse_Temp_2D_TSED,1),size(Gl_Reverse_Temp_2D_TSED,2), Gl_Reverse_Temp_2D_TSED)

  ! *** Reverse Order of Array
  DO L=2,LA
    DO K=1,KB
      Reverse_Temp_2D_BULKDENS(L,K) = BULKDENS(K,L)
    ENDDO
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_BULKDENS,1), size(Reverse_Temp_2D_BULKDENS,2), Reverse_Temp_2D_BULKDENS, &
    size(Gl_Reverse_Temp_2D_BULKDENS,1),size(Gl_Reverse_Temp_2D_BULKDENS,2), Gl_Reverse_Temp_2D_BULKDENS)

  ! *** Reverse Order of Array
  DO L=2,LA
    DO K=1,KB
      DO nn=1,NSCM
        Reverse_Temp_3D(L,K,nn) = PERSED(nn,K,L)
      ENDDO
    ENDDO
  ENDDO

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_3D,1), size(Reverse_Temp_3D,2), size(Reverse_Temp_3D,3), Reverse_Temp_3D, &
    size(Gl_Reverse_Temp_3D,1),size(Gl_Reverse_Temp_3D,2),size(Gl_Reverse_Temp_3D,3), Gl_Reverse_Temp_3D)


  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(D50AVG,1), D50AVG, size(D50AVG_Global,1), D50AVG_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(ETOTO,1), ETOTO, size(ETOTO_Global,1), ETOTO_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(DEPO,1), DEPO, size(DEPO_Global,1), DEPO_Global)

  IF( ICALC_BL > 0 )THEN
    !
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(CBL,1), size(CBL,2), CBL, &
      size(CBL_Global,1),size(CBL_Global,2), CBL_Global)
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(QSBDLDX,1), size(QSBDLDX,2), QSBDLDX, &
      size(QSBDLDX_Global,1),size(QSBDLDX_Global,2), QSBDLDX_Global)
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(QSBDLDY,1), size(QSBDLDY,2), QSBDLDY, &
      size(QSBDLDY_Global,1),size(QSBDLDY_Global,2), QSBDLDY_Global)

  ENDIF

  IF( ISTRAN(5) > 0 )THEN
    IF( ICALC_BL > 0 )THEN
      j = j + 1
      Call Assign_Loc_Glob_For_Write(j, size(CBLTOX,1), size(CBLTOX,2), CBLTOX, &
        size(CBLTOX_Global,1),size(CBLTOX_Global,2), CBLTOX_Global)
    ENDIF
  ENDIF

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  ! *** Reverse global array order
  Do k = 1, KB
    DO l=2, LA_Global
      TSED_Global(K,L) = Gl_Reverse_Temp_2D_TSED(L,K)
    End do
  End do

  ! *** Reverse global array order
  Do k = 1, KB
    DO l = 2, LA_Global
      BULKDENS_Global(K,L) = Gl_Reverse_Temp_2D_BULKDENS(L,K)
    End do
  End do

  ! *** Reverse global array order
  DO nn = 1, NSCM
    Do k = 1, KB
      DO l = 2, LA_Global
        PERSED_Global(nn,k,l) = Gl_Reverse_Temp_3D(L,K,nn)
      End do
    End do
  End do

  Deallocate(Reverse_Temp_2D_TSED)
  Deallocate(Reverse_Temp_2D_BULKDENS)
  Deallocate(Reverse_Temp_3D)
  Deallocate(Gl_Reverse_Temp_2D_TSED)
  Deallocate(Gl_Reverse_Temp_2D_BULKDENS)
  Deallocate(Gl_Reverse_Temp_3D)

  End Subroutine Map_Write_SEDZLJOUT

  ! *** BEDOUT
  Subroutine Map_Write_BEDOUT

  Implicit None

  ! *** Local Variables
  Integer :: i, j

  j = 0

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(KBT,1), KBT, size(KBT_Global,1), KBT_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(HBED,1), size(HBED,2), HBED,  &
    size(HBED_Global,1), size(HBED_Global,2), HBED_Global)
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(BDENBED,1), size(BDENBED,2), BDENBED,  &
    size(BDENBED_Global,1), size(BDENBED_Global,2), BDENBED_Global)
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(PORBED,1), size(PORBED,2), PORBED,  &
    size(PORBED_Global,1), size(PORBED_Global,2), PORBED_Global)


  IF( ISTRAN(6) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SEDB,1), size(SEDB,2), size(SEDB,3), SEDB, &
      size(SEDB_Global,1), size(SEDB_Global,2), size(SEDB_Global,3), SEDB_Global)
  ENDIF

  IF( ISTRAN(7) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SNDB,1), size(SNDB,2), size(SNDB,3), SNDB, &
      size(SNDB_Global,1), size(SNDB_Global,2), size(SNDB_Global,3), SNDB_Global)
  ENDIF

  IF( ISTRAN(5) >= 1 )THEN
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOXB,1), size(TOXB,2), size(TOXB,3), TOXB, &
      size(TOXB_Global,1), size(TOXB_Global,2), size(TOXB_Global,3), TOXB_Global)
  ENDIF

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  End Subroutine Map_Write_BEDOUT

  ! *** WQOUT
  Subroutine Map_Write_WQOUT

  Implicit None

  ! *** Local Variables
  Integer :: i, j

  j = 0

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQV,1),        size(WQV,2),        size(WQV,3),        WQV, &
    size(WQV_Global,1), size(WQV_Global,2), size(WQV_Global,3), WQV_Global)

  num_arrays_to_write_out = j

  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  End Subroutine Map_Write_WQOUT

  ! *** SDOUT
  Subroutine Map_Write_SDOUT

  Implicit None

  ! *** Local variables
  Integer :: i, j

  j = 0

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMPON,1), size(SMPON,2), SMPON,  &
    size(SMPON_Global,1), size(SMPON_Global,2), SMPON_Global)
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMPOP,1), size(SMPOP,2), SMPOP,  &
    size(SMPOP_Global,1), size(SMPOP_Global,2), SMPOP_Global)
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMPOC,1), size(SMPOC,2), SMPOC,  &
    size(SMPOC_Global,1), size(SMPOC_Global,2), SMPOC_Global)
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMDFN,1), size(SMDFN,2), SMDFN,  &
    size(SMDFN_Global,1), size(SMDFN_Global,2), SMDFN_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMDFP,1), size(SMDFP,2), SMDFP,  &
    size(SMDFP_Global,1), size(SMDFP_Global,2), SMDFP_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMDFC,1), size(SMDFC,2), SMDFC,  &
    size(SMDFC_Global,1), size(SMDFC_Global,2), SMDFC_Global)
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM1NH4,1), SM1NH4, size(SM1NH4_Global,1), SM1NH4_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM2NH4,1), SM2NH4, size(SM2NH4_Global,1), SM2NH4_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM1NO3,1), SM1NO3, size(SM1NO3_Global,1), SM1NO3_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM2NO3,1), SM2NO3, size(SM2NO3_Global,1), SM2NO3_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM1PO4,1), SM1PO4, size(SM1PO4_Global,1), SM1PO4_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM2PO4,1), SM2PO4, size(SM2PO4_Global,1), SM2PO4_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM1H2S,1), SM1H2S, size(SM1H2S_Global,1), SM1H2S_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM2H2S,1), SM2H2S, size(SM2H2S_Global,1), SM2H2S_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM1SI,1), SM1SI, size(SM1SI_Global,1), SM1SI_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SM2SI,1), SM2SI, size(SM2SI_Global,1), SM2SI_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMPSI,1), SMPSI, size(SMPSI_Global,1), SMPSI_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMBST,1), SMBST, size(SMBST_Global,1), SMBST_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMT,1), SMT, size(SMT_Global,1), SMT_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMCSOD,1), SMCSOD, size(SMCSOD_Global,1), SMCSOD_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(SMNSOD,1), SMNSOD, size(SMNSOD_Global,1), SMNSOD_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQBFNH4,1), WQBFNH4, size(WQBFNH4_Global,1), WQBFNH4_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQBFNO3,1), WQBFNO3, size(WQBFNO3_Global,1), WQBFNO3_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQBFO2,1), WQBFO2, size(WQBFO2_Global,1), WQBFO2_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQBFCOD,1), WQBFCOD, size(WQBFCOD_Global,1), WQBFCOD_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQBFPO4D,1), WQBFPO4D, size(WQBFPO4D_Global,1), WQBFPO4D_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQBFSAD,1), WQBFSAD, size(WQBFSAD_Global,1), WQBFSAD_Global)

  num_arrays_to_write_out = j

  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  End Subroutine Map_Write_SDOUT

  ! *** RPEMOUT
  Subroutine Map_Write_RPEMOUT

  Implicit None

  ! *** Local Variables
  Integer :: i, j

  j = 0

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQRPS,1), WQRPS, size(WQRPS_Global,1), WQRPS_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQRPR,1), WQRPR, size(WQRPR_Global,1), WQRPR_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQRPE,1), WQRPE, size(WQRPE_Global,1), WQRPE_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(WQRPD,1), WQRPD, size(WQRPD_Global,1), WQRPD_Global)

  num_arrays_to_write_out = j

  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

  End Subroutine Map_Write_RPEMOUT

  ! @todo Develop method for handling arrays out
  !---------------------------------------------------------------------------!
  ! Subroutine: Handle_Calls_MapGatherSort
  !
  !> @details makes calls to 1D,2D,3D Map/Gather/Sort subroutines
  !
  !> @param[in] num_arrays_to_write: max number of arrays to write out for the
  ! respective routine
  !
  !> @author Zander Mausolff
  !> @date 1/9/2019
  !---------------------------------------------------------------------------!
  Subroutine Handle_Calls_MapGatherSort(num_arrays_to_write)

  Use Mod_Map_Gather_Sort

  Implicit None

  ! *** Read in Variables
  Integer, Intent(in) :: num_arrays_to_write

  ! *** Local Variables
  Integer :: i

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Do i = 1, num_arrays_to_write

    ! *** Write out 1D integer values
    IF( Dim_Array_Written_Out(i) == 10 )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_1D_Int(i).val, 1), &    ! size of the first dim of the local soln
        Local_Arrays_to_Write_1D_Int(i).val, &           ! local array from each process
        size(Global_Arrays_to_Write_1D_Int(i).val, 1), & ! size of the first dim of the global soln
        Global_Arrays_to_Write_1D_Int(i).val)
    end if

    ! *** Write out 2D integer values
    IF( Dim_Array_Written_Out(i) == 20 )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_2D_Int(i).val, 1), &      ! size of the first dim of the local soln
        size(Local_Arrays_to_Write_2D_Int(i).val, 2), &  ! size of the second dim of the local soln
        Local_Arrays_to_Write_2D_Int(i).val, &           ! local array from each process
        size(Global_Arrays_to_Write_2D_Int(i).val, 1), & ! size of the first dim of the global soln
        size(Global_Arrays_to_Write_2D_Int(i).val, 2), & ! size of the second dim of the global soln
        Global_Arrays_to_Write_2D_Int(i).val)
    end if

    ! *** Write out 3D integer values
    IF( Dim_Array_Written_Out(i) == 30 )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_3D_Int(i).val, 1), &
        size(Local_Arrays_to_Write_3D_Int(i).val, 2), &
        size(Local_Arrays_to_Write_3D_Int(i).val, 3), &
        Local_Arrays_to_Write_3D_Int(i).val, &
        size(Global_Arrays_to_Write_3D_Int(i).val, 1), &
        size(Global_Arrays_to_Write_3D_Int(i).val, 2), &
        size(Global_Arrays_to_Write_3D_Int(i).val, 3), &
        Global_Arrays_to_Write_3D_Int(i).val )
    end if

    ! *** Write out 1D arrays of real(4)
    IF( Dim_Array_Written_Out(i) == 1 )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_1D_Real(i).val, 1), &    ! size of the first dim of the local soln
        Local_Arrays_to_Write_1D_Real(i).val, &           ! local array from each process
        size(Global_Arrays_to_Write_1D_Real(i).val, 1), & ! size of the first dim of the global soln
        Global_Arrays_to_Write_1D_Real(i).val)
    endif

    ! *** Write out 2D arrays of real(4)
    IF( Dim_Array_Written_Out(i) == 2 )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_2D_Real(i).val, 1), &      ! size of the first dim of the local soln
        size(Local_Arrays_to_Write_2D_Real(i).val, 2), &  ! size of the second dim of the local soln
        Local_Arrays_to_Write_2D_Real(i).val, &           ! local array from each process
        size(Global_Arrays_to_Write_2D_Real(i).val, 1), & ! size of the first dim of the global soln
        size(Global_Arrays_to_Write_2D_Real(i).val, 2), & ! size of the second dim of the global soln
        Global_Arrays_to_Write_2D_Real(i).val)
    endif

    ! *** Write out 3D arrays of real(4)
    IF( Dim_Array_Written_Out(i) == 3 )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_3D_Real(i).val, 1), &
        size(Local_Arrays_to_Write_3D_Real(i).val, 2), &
        size(Local_Arrays_to_Write_3D_Real(i).val, 3), &
        Local_Arrays_to_Write_3D_Real(i).val, &
        size(Global_Arrays_to_Write_3D_Real(i).val, 1), &
        size(Global_Arrays_to_Write_3D_Real(i).val, 2), &
        size(Global_Arrays_to_Write_3D_Real(i).val, 3), &
        Global_Arrays_to_Write_3D_Real(i).val )
    endif

    ! **************
    ! *** Write out 1D arrays of real(rkd)
    IF( Dim_Array_Written_Out(i) == 100 + rkd )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_1D_Real_RK8(i).val, 1), & ! size of the first dim of the local soln
        Local_Arrays_to_Write_1D_Real_RK8(i).val, &  ! local array from each process
        size(Global_Arrays_to_Write_1D_Real_RK8(i).val, 1), & ! size of the first dim of the global soln
        Global_Arrays_to_Write_1D_Real_RK8(i).val)
    endif

    ! *** Write out 2D arrays of real (rkd)
    IF( Dim_Array_Written_Out(i) == 200 + rkd )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_2D_Real_RK8(i).val, 1), & ! size of the first dim of the local soln
        size(Local_Arrays_to_Write_2D_Real_RK8(i).val, 2), & ! size of the second dim of the local soln
        Local_Arrays_to_Write_2D_Real_RK8(i).val, &  ! local array from each process
        size(Global_Arrays_to_Write_2D_Real_RK8(i).val, 1), & ! size of the first dim of the global soln
        size(Global_Arrays_to_Write_2D_Real_RK8(i).val, 2), & ! size of the second dim of the global soln
        Global_Arrays_to_Write_2D_Real_RK8(i).val)
    endif

    ! *** Write out 3D arrays of real(4)
    IF( Dim_Array_Written_Out(i) == 300 + rkd )then
      Call Map_Gather_Sort(size(Local_Arrays_to_Write_3D_Real_RK8(i).val, 1), &
        size(Local_Arrays_to_Write_3D_Real_RK8(i).val, 2), &
        size(Local_Arrays_to_Write_3D_Real_RK8(i).val, 3), &
        Local_Arrays_to_Write_3D_Real_RK8(i).val, &
        size(Global_Arrays_to_Write_3D_Real_RK8(i).val, 1), &
        size(Global_Arrays_to_Write_3D_Real_RK8(i).val, 2), &
        size(Global_Arrays_to_Write_3D_Real_RK8(i).val, 3), &
        Global_Arrays_to_Write_3D_Real_RK8(i).val )
    endif
  End do

  ! *** Clear the arrays
  !Local_Arrays_to_Write_1D_Real
  !Global_Arrays_to_Write_1D_Real
  !
  !Local_Arrays_to_Write_2D_Real
  !Global_Arrays_to_Write_2D_Real
  !
  !Local_Arrays_to_Write_3D_Real
  !Global_Arrays_to_Write_3D_Real
  !
  !! *** Real(8)
  !Local_Arrays_to_Write_1D_Real_RK8
  !Global_Arrays_to_Write_1D_Real_RK8
  !
  !Local_Arrays_to_Write_2D_Real_RK8
  !Global_Arrays_to_Write_2D_Real_RK8
  !
  !Local_Arrays_to_Write_3D_Real_RK8
  !Global_Arrays_to_Write_3D_Real_RK8
  !
  !! *** Integer
  !Local_Arrays_to_Write_1D_Int
  !Global_Arrays_to_Write_1D_Int
  !
  !Local_Arrays_to_Write_2D_Int
  !Global_Arrays_to_Write_2D_Int
  !
  !Local_Arrays_to_Write_3D_Int
  !Global_Arrays_to_Write_3D_Int

  Dim_Array_Written_Out(:) = 0

  End Subroutine Handle_Calls_MapGatherSort

  End module Mod_Map_Write_EE_Binary
