! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!  !---------------------------------------------------------------------------!
!  !                     EFDC+ Developed by DSI, LLC.
!  !---------------------------------------------------------------------------!
!  ! @details This routine calls all relevant subroutines that collects the
!  ! solution from each subdomain and remaps to the global in preperation for
!  ! writing out to the EE binary files
!  ! @author Zander Mausolff
!  ! @date 11/01/2019
!  !---------------------------------------------------------------------------!
!
!<<<<<<< HEAD
!  Subroutine Global_Mapping_Routines(LA_Local_no_ghost)
!
!  Use Global
!  Use MPI
!  Use Variables_MPI
!  Use Variables_MPI_Concentration
!  Use Variables_MPI_Mapping
!  Use Variables_MPI_Write_Out
!
!  Implicit None
!
!  ! *** Passed in
!  Integer, intent(in) :: LA_Local_no_ghost
!
!  ! *** Local variables
!  Integer :: i, j, ierr, L, k, md, ii, ns
!  Integer :: num_active_l_local
!
!  ! *** L array variables (real)
!  Real, Allocatable, Dimension(:)    :: Soln_Local_1D_L
!  Real, Allocatable, Dimension(:)    :: Soln_Global_1D_L
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_L
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_L
!
!  ! *** L array variables (integer)
!  Integer, Allocatable, Dimension(:) :: Soln_Local_1D_L_I
!  Integer, Allocatable, Dimension(:) :: Soln_Global_1D_L_I
!
!  ! *** L,KC array variables
!  Real, Allocatable, Dimension(:)    :: Soln_Local_1D_KC
!  Real, Allocatable, Dimension(:)    :: Soln_Global_1D_KC
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_KC
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_KC
!
!  ! *** L,KC,NDYE array variables (by DYE class)
!  Real, Allocatable, Dimension(:) :: Soln_Local_1D_DYE
!  Real, Allocatable, Dimension(:) :: Soln_Global_1D_DYE
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_DYE
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_DYE
!
!  ! *** L,KC,NTOX array variables (by TOX class)
!  Real, Allocatable, Dimension(:) :: Soln_Local_1D_TOX
!  Real, Allocatable, Dimension(:) :: Soln_Global_1D_TOX
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_TOX
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_TOX
!
!  ! *** L,KC,NSED array variables (by SED class)
!  Real, Allocatable, Dimension(:) :: Soln_Local_1D_SED
!  Real, Allocatable, Dimension(:) :: Soln_Global_1D_SED
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_SED
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_SED
!
!  ! *** L,KC,NSND array variables (by SND class)
!  Real, Allocatable, Dimension(:) :: Soln_Local_1D_SND
!  Real, Allocatable, Dimension(:) :: Soln_Global_1D_SND
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_SND
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_SND
!
!  ! *** L,0:KC array variables
!  Real, Allocatable, Dimension(:)    ::Soln_Local_1D_QQ
!  Real, Allocatable, Dimension(:)    ::Soln_Global_1D_QQ
!  Integer, Allocatable, Dimension(:) ::Map_Local_L_to_Global_QQ
!  Integer, Allocatable, Dimension(:) ::Global_Local_L_to_Global_QQ
!
!  ! *** L,KB array variables
!  Real, Allocatable, Dimension(:)    :: Soln_Local_1D_KB
!  Real, Allocatable, Dimension(:)    :: Soln_Global_1D_KB
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_KB
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_KB
!
!  ! *** L,KB,NTOX array variables
!  Real, Allocatable, Dimension(:)    :: Soln_Local_1D_TOXB
!  Real, Allocatable, Dimension(:)    :: Soln_Global_1D_TOXB
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_TOXB
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_TOXB
!
!  ! *** L,KB,NSED array variables
!  Real, Allocatable, Dimension(:)    :: Soln_Local_1D_SEDB
!  Real, Allocatable, Dimension(:)    :: Soln_Global_1D_SEDB
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_SEDB
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_SEDB
!
!  ! *** L,KB,NSND array variables
!  Real, Allocatable, Dimension(:)    :: Soln_Local_1D_SNDB
!  Real, Allocatable, Dimension(:)    :: Soln_Global_1D_SNDB
!  Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_SNDB
!  Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_SNDB
!
!  Real, Allocatable, Dimension(:,:)   :: TMPBED
!  Real, Allocatable, Dimension(:,:,:) :: BEDBYCLASS
!  
!  !---------------------------------------------------------------------------!
!
!  If(.not.Allocated(Soln_Local_1D_KC) )THEN
!
!    Allocate(Soln_Local_1D_L(LA_Local_no_ghost))
!    Allocate(Soln_Global_1D_L(LCM_Global))
!    Allocate(Map_Local_L_to_Global_L(LA_Local_no_ghost ))
!    Allocate(Global_Local_L_to_Global_L(LCM_Global))
!
!    Allocate(Soln_Local_1D_L_I(LA_Local_no_ghost))
!    Allocate(Soln_Global_1D_L_I(LCM_Global))
!
!    !Allocate(Soln_Local_1D_KC(LA_Local_no_ghost*kc))
!    !Allocate(Soln_Global_1D_KC(LCM_Global*kc))
!    !Allocate(Map_Local_L_to_Global_KC(LA_local_no_ghost*KC))
!    !Allocate(Global_Local_L_to_Global_KC(LCM_Global*KC))
!
!    Allocate(Soln_Local_1D_KC(LA_Local_no_ghost*kcm))
!    Allocate(Soln_Global_1D_KC(LCM_Global*kcm))
!    Allocate(Map_Local_L_to_Global_KC(LA_local_no_ghost*kcm))
!    Allocate(Global_Local_L_to_Global_KC(LCM_Global*kcm))
!
!    Allocate(Soln_Local_1D_KB(LA_Local_no_ghost*KB))
!    Allocate(Soln_Global_1D_KB(LCM_Global*KB))
!    Allocate(Map_Local_L_to_Global_KB(LA_local_no_ghost*KB))
!    Allocate(Global_Local_L_to_Global_KB(LCM_Global*KB))
!
!    IF( ISTRAN(3) > 0 )THEN
!      Allocate(Soln_Local_1D_DYE(La_local_no_ghost*KC*NDYM))
!      Allocate(Soln_Global_1D_DYE(LCM_Global*KC*NDYM))
!      Allocate(Map_Local_L_to_Global_DYE(LA_local_no_ghost*KC*NDYM))
!      Allocate(Global_Local_L_to_Global_DYE(LCM_Global*KC*NDYM))
!    ENDIF
!
!    IF( ISTRAN(5) > 0 )THEN
!      Allocate(Soln_Local_1D_TOX(La_local_no_ghost*KC*NTXM))
!      Allocate(Soln_Global_1D_TOX(LCM_Global*KC*NTXM))
!      Allocate(Map_Local_L_to_Global_TOX(LA_local_no_ghost*KC*NTXM))
!      Allocate(Global_Local_L_to_Global_TOX(LCM_Global*KC*NTXM))
!
!      Allocate(Soln_Local_1D_TOXB(La_local_no_ghost*KB*NTXM))
!      Allocate(Soln_Global_1D_TOXB(LCM_Global*KB*NTXM))
!      Allocate(Map_Local_L_to_Global_TOXB(LA_local_no_ghost*KB*NTXM))
!      Allocate(Global_Local_L_to_Global_TOXB(LCM_Global*KB*NTXM))
!    ENDIF
!
!    IF( ISTRAN(6) > 0 )THEN
!      Allocate(Soln_Local_1D_SED(La_local_no_ghost*KC*NSCM))
!      Allocate(Soln_Global_1D_SED(LCM_Global*KC*NSCM))
!      Allocate(Map_Local_L_to_Global_SED(LA_local_no_ghost*KC*NSCM))
!      Allocate(Global_Local_L_to_Global_SED(LCM_Global*KC*NSCM))
!      
!      Allocate(Soln_Local_1D_SEDB(La_local_no_ghost*KB*NSCM))
!      Allocate(Soln_Global_1D_SEDB(LCM_Global*KB*NSCM))
!      Allocate(Map_Local_L_to_Global_SEDB(LA_local_no_ghost*KB*NSCM))
!      Allocate(Global_Local_L_to_Global_SEDB(LCM_Global*KB*NSCM))
!    ENDIF
!
!    IF( ISTRAN(7) > 0 )THEN
!      Allocate(Soln_Local_1D_SND(La_local_no_ghost*KC*NSNM))
!      Allocate(Soln_Global_1D_SND(LCM_Global*KC*NSNM))
!      Allocate(Map_Local_L_to_Global_SND(LA_local_no_ghost*KC*NSNM))
!      Allocate(Global_Local_L_to_Global_SND(LCM_Global*KC*NSNM))
!      
!      Allocate(Soln_Local_1D_SNDB(La_local_no_ghost*KB*NSNM))
!      Allocate(Soln_Global_1D_SNDB(LCM_Global*KB*NSNM))
!      Allocate(Map_Local_L_to_Global_SNDB(LA_local_no_ghost*KB*NSNM))
!      Allocate(Global_Local_L_to_Global_SNDB(LCM_Global*KB*NSNM))
!    ENDIF
!
!    Allocate(Soln_Local_1D_QQ(LA_Local_no_ghost*kc+1)) ! + 1 because QQ goes from 0 to KCM
!    Allocate(Soln_Global_1D_QQ(LCM_Global*kc+1))
!    Allocate(Map_Local_L_to_Global_QQ(LA_local_no_ghost*KC+1))
!    Allocate(Global_Local_L_to_Global_QQ(LCM_Global*KC+1))
!  end if
!
!  IF( LSEDZLJ )THEN
!    ALLOCATE(TMPBED(LCM, KBM))
!    ALLOCATE(BEDBYCLASS(LCM, KBM, NSCM))
!  ENDIF
!  
!  !Soln_Local_1D_KB   = 0.0
!  !Soln_Global_1D_KB  = 0.0
!
!  !Soln_Local_1D_QQ  = 0.0
!  !Soln_Global_1D_QQ = 0.0
!  !Map_Local_L_to_Global_QQ = 0
!  !Global_Local_L_to_Global_QQ = 0
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** HP
!  ! *** Map to 1d array containing l values excluding ghost cells and build up array mapping to global l
!  Soln_Local_1D_L            = 0.0
!  Soln_Global_1D_L           = 0.0
!  Map_Local_L_to_Global_L    = 0
!  Global_Local_L_to_Global_L = 0
!
!  Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, HP, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!  ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                 Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!
!  if( process_id == master_id )THEN
!    ! *** Sort the global array so L indexing is correct
!    Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, HP_Global)
!  endif
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** U
!  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  Soln_Local_1D_KC            = 0.0
!  Soln_Global_1D_KC           = 0.0
!  Map_Local_L_to_Global_KC    = 0
!  Global_Local_L_to_Global_KC = 0
!
!  Call Map3D_to_1D(la_local_no_ghost, num_active_l_local, KC, U, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!
!  ! *** Send local values containing values ommitting the ghost cells to master
!  Call Gather_3D(la_local_no_ghost, num_active_l_local, LCM, KC, Soln_Local_1D_KC, Soln_Global_1D_KC, &
!                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!
!  ! *** Handle mapping of global values to the master
!  if( process_id == master_id  )THEN
!    Call Sort_Global_Soln_3D(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, U_Global)
!  end if
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** V
!  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  Soln_Local_1D_KC            = 0.0
!  Soln_Global_1D_KC           = 0.0
!  Map_Local_L_to_Global_KC    = 0
!  Global_Local_L_to_Global_KC = 0
!
!  Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KC, V, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!
!  ! *** Send local values containing values ommitting the ghost cells to master
!  Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KC, Soln_Local_1D_KC, Soln_Global_1D_KC, &
!                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!
!  ! *** Handle mapping of global values to the master
!  if( process_id == master_id  )THEN
!    Call Sort_Global_Soln_3D(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, V_Global)
!  end if
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** W
!  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  Soln_Local_1D_KC            = 0.0
!  Soln_Global_1D_KC           = 0.0
!  Map_Local_L_to_Global_KC    = 0
!  Global_Local_L_to_Global_KC = 0
!
!  Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KC, W, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!
!  ! *** Send local values containing values ommitting the ghost cells to master
!  Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KC, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!
!  ! *** Handle mapping of global values to the master
!  if( process_id == master_id  )THEN
!    Call Sort_Global_Soln_3D(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, W_Global)
!  end if
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** Bed Shear
!  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  Soln_Local_1D_L            = 0.0
!  Soln_Global_1D_L           = 0.0
!  Map_Local_L_to_Global_L    = 0
!  Global_Local_L_to_Global_L = 0
!
!  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
!    IF( LSEDZLJ )THEN
!      ! *** TAU(LCM) Shear Stress in dynes/cm^2,  WATERDENS  gm/cm3
!      ! *** SHEAR IS SAVED AS N/M^2 (Pascals) normalized to water density  (M2/S2)
!      DO L=2,LA
!        IF( LBED(L) .OR. TIMEDAY < SEDSTART )THEN
!          SHEAR_Local(L) = QQ(L,0)/CTURB2
!        ELSE
!          SHEAR_Local(L) = TAU(L) * 0.1 / (WATERDENS*1000.)
!        ENDIF
!      ENDDO
!    ELSEIF( ISBEDSTR >= 1 )THEN
!      DO L=2,LA
!        IF( LBED(L) )THEN
!          SHEAR_Local(L) = QQ(L,0)/CTURB2
!        ELSE
!          SHEAR_Local(L) = TAUBSED(L)
!        ENDIF
!      ENDDO
!
!      IF( ISBEDSTR == 1 )THEN
!        DO L=2,LA
!          IF( LBED(L) )THEN
!            SHEAR_Local(L) = QQ(L,0)/CTURB2
!          ELSE
!            SHEAR_Local(L) = TAUBSND(L)
!          ENDIF
!        ENDDO
!      ENDIF
!    ELSE
!      ! *** TOTAL BED SHEAR STRESS
!      DO L=2,LA
!        IF( LBED(L) )THEN
!          SHEAR_Local(L) = QQ(L,0)/CTURB2
!        ELSE
!          SHEAR_Local(L) = TAUB(L)
!        ENDIF
!      ENDDO
!    ENDIF
!  ELSE
!    ! *** TOTAL BED SHEAR STRESS
!    DO L=2,LA  
!      SHEAR_Local(L) = MAX(QQ(L,0),QQMIN)/CTURB2
!    ENDDO
!  ENDIF
!  
!  ! *** Now gather the calculated bed shears
!  Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, SHEAR_Local, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!  ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                 Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!  if( process_id == master_id )THEN
!    ! *** Sort the global array so L indexing is correct
!    Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, SHEAR_Global)
!  endif
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** Map Salinity
!  IF(ISTRAN(1).EQ.1 )THEN
!    ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!    Soln_Local_1D_KC            = 0.0
!    Soln_Global_1D_KC           = 0.0
!    Map_Local_L_to_Global_KC    = 0
!    Global_Local_L_to_Global_KC = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KC, SAL, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KC, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!                   Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, SAL_Global)
!    end if
!  End if
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** Map Temperature
!  IF(ISTRAN(2).EQ.1 )THEN
!    ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!    Soln_Local_1D_KC            = 0.0
!    Soln_Global_1D_KC           = 0.0
!    Map_Local_L_to_Global_KC    = 0
!    Global_Local_L_to_Global_KC = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KC, TEM, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KC, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!                   Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, TEM_Global)
!    end if
!
!    IF( TBEDIT > 0.)THEN
!      ! *** Bed Temperature
!      Soln_Local_1D_L  = 0.0
!      Soln_Global_1D_L = 0.0
!      Map_Local_L_to_Global_L    = 0
!      Global_Local_L_to_Global_L = 0
!
!      Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, TEMB, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!      ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!      Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                     Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!      if( process_id == master_id )THEN
!        ! *** Sort the global array so L indexing is correct
!        Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, TEMB_Global)
!      endif
!    ENDIF
!
!  End if
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** If dye is turned on then we should map to global domain so we can write it out
!  IF( ISTRAN(3) >= 1  )THEN
!    Soln_Local_1D_DYE            = 0.0
!    Soln_Global_1D_DYE           = 0.0
!    Map_Local_L_to_Global_DYE    = 0
!    Global_Local_L_to_Global_DYE = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NDYE, NDYM, DYE, Soln_Local_1D_DYE, &
!                     Map_Local_L_to_Global_DYE)
!         
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NDYM, Soln_Local_1D_DYE, Soln_Global_1D_DYE, &
!                   Map_Local_L_to_Global_DYE, Global_Local_L_to_Global_DYE)
!
!    ! *** Sorting step onto entire spatial domain
!    if( process_id == master_id )THEN
!      Call Sort_Global_Soln_4D(KC,KCM, NDYM, Soln_Global_1D_DYE, Global_Local_L_to_Global_DYE, DYE_Global)
!    endif
!
!  endif
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** If SFL is turned on then we should map to global domain so we can write it out
!  IF( ISTRAN(4) >= 1  )THEN
!    ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!    Soln_Local_1D_KC            = 0.0
!    Soln_Global_1D_KC           = 0.0
!    Map_Local_L_to_Global_KC    = 0
!    Global_Local_L_to_Global_KC = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KC, SFL, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KC, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!                   Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, SFL_Global)
!    end if
!
!  endif
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** If TOX is turned on then we should map to global domain so we can write it out
!  IF( ISTRAN(5) >= 1  )THEN
!    Soln_Local_1D_TOX            = 0.0
!    Soln_Global_1D_TOX           = 0.0
!    Map_Local_L_to_Global_TOX    = 0
!    Global_Local_L_to_Global_TOX = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NTOX, NTOX, TOX, Soln_Local_1D_TOX, &
!                     Map_Local_L_to_Global_TOX)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NTOX, Soln_Local_1D_TOX, Soln_Global_1D_TOX, &
!                   Map_Local_L_to_Global_TOX, Global_Local_L_to_Global_TOX)
!
!    ! *** Sorting step onto entire spatial domain
!    if( process_id == master_id )THEN
!      Call Sort_Global_Soln_4D(KC, KCM, NDYM, Soln_Global_1D_TOX, Global_Local_L_to_Global_TOX, TOX_Global)
!    endif
!
!    
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** TOXB - Bed Toxics
!    Soln_Local_1D_TOXB            = 0.0
!    Soln_Global_1D_TOXB           = 0.0
!    Map_Local_L_to_Global_TOXB    = 0
!    Global_Local_L_to_Global_TOXB = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KB, KBM, NTOX, NTXM, TOXB, Soln_Local_1D_TOXB, &
!                     Map_Local_L_to_Global_TOXB)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KB, NTOX, Soln_Local_1D_TOXB, Soln_Global_1D_TOXB, &
!                   Map_Local_L_to_Global_TOXB, Global_Local_L_to_Global_TOXB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_4D(KB, KBM, NTXM, Soln_Global_1D_TOXB, Global_Local_L_to_Global_TOXB, TOXB_Global)
!    end if
!  ENDIF
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** If SED or SND is turned on then we should map to global domain so we can write it out
!  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
!    Soln_Local_1D_L_I          = 0
!    Soln_Global_1D_L_I         = 0
!    Map_Local_L_to_Global_L    = 0
!    Global_Local_L_to_Global_L = 0
!
!    ! *** KBT - Sediment Bed Top
!    Call Map1D_no_ghost_int(la_local_no_ghost, num_active_l_local, KBT, Soln_Local_1D_L_I, Map_Local_L_to_Global_L)
!
!    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!    Call Gather_1D_int(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L_I, Soln_Global_1D_L, &
!                   Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!    if( process_id == master_id )THEN
!      ! *** Sort the global array so L indexing is correct
!      Call Sort_Global_Soln_HP_Int(Soln_Global_1D_L, Global_Local_L_to_Global_L, KBT_Global)
!    endif
!
!    !---------------------------------------------------------------------------------------------------------!
!    IF( IMORPH > 0 )THEN
!      ! *** BELV - Bottom Elevation
!      Soln_Local_1D_L            = 0.0
!      Soln_Global_1D_L           = 0.0
!      Map_Local_L_to_Global_L    = 0
!      Global_Local_L_to_Global_L = 0
!
!      Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, BELV, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!      ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!      Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                     Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!
!      if( process_id == master_id )THEN
!        ! *** Sort the global array so L indexing is correct
!        Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, BELV_Global)
!      endif
!    ENDIF
!    
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** HBED - Bed Thickness
!    Soln_Local_1D_KB            = 0.0
!    Soln_Global_1D_KB           = 0.0
!    Map_Local_L_to_Global_KB    = 0
!    Global_Local_L_to_Global_KB = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KB, HBED, Soln_Local_1D_KB, Map_Local_L_to_Global_KB)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!                   Map_Local_L_to_Global_KB, Global_Local_L_to_Global_KB)
!
!      ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KB, Soln_Global_1D_KB, Global_Local_L_to_Global_KB, HBED_Global)
!    end if
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** PORBED - Bed Porosity
!    Soln_Local_1D_KB            = 0.0
!    Soln_Global_1D_KB           = 0.0
!    Map_Local_L_to_Global_KB    = 0
!    Global_Local_L_to_Global_KB = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KB, PORBED, Soln_Local_1D_KB, Map_Local_L_to_Global_KB)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!                   Map_Local_L_to_Global_KB, Global_Local_L_to_Global_KB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KB, Soln_Global_1D_KB, Global_Local_L_to_Global_KB, PORBED_Global)
!    end if
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** BDENBED - Bed Bulk Density
!    Soln_Local_1D_KB            = 0.0
!    Soln_Global_1D_KB           = 0.0
!    Map_Local_L_to_Global_KB    = 0
!    Global_Local_L_to_Global_KB = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KB, BDENBED, Soln_Local_1D_KB, Map_Local_L_to_Global_KB)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!                   Map_Local_L_to_Global_KB, Global_Local_L_to_Global_KB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KB, Soln_Global_1D_KB, Global_Local_L_to_Global_KB, BDENBED_Global)
!    end if
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** VFRBED - Bed Volume Fraction
!    Soln_Local_1D_KB            = 0.0
!    Soln_Global_1D_KB           = 0.0
!    Map_Local_L_to_Global_KB    = 0
!    Global_Local_L_to_Global_KB = 0
!
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KB, VFRBED, Soln_Local_1D_KB, Map_Local_L_to_Global_KB)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!                   Map_Local_L_to_Global_KB, Global_Local_L_to_Global_KB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KB, Soln_Global_1D_KB, Global_Local_L_to_Global_KB, VFRBED_Global)
!    end if
!
!    IF( ISTRAN(7) > 0 .AND. ISBDLDBC > 0 .AND. NSND > 0 )THEN
!      QSBDLDX_Global = -999   ! delme - todo
!      QSBDLDY_Global = -999   ! delme - todo
!    ENDIF
!  ENDIF
!  
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** If SED is turned on then we should map to global domain so we can write it out
!  IF( ISTRAN(6) >= 1  )THEN
!    Soln_Local_1D_SED            = 0.0
!    Soln_Global_1D_SED           = 0.0
!    Map_Local_L_to_Global_SED    = 0
!    Global_Local_L_to_Global_SED = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NSED, NSCM, SED, Soln_Local_1D_SED, &
!                     Map_Local_L_to_Global_SED)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NSCM, Soln_Local_1D_SED, Soln_Global_1D_SED, &
!                   Map_Local_L_to_Global_SED, Global_Local_L_to_Global_SED)
!
!    ! *** Sorting step onto entire spatial domain
!    if( process_id == master_id )THEN
!      Call Sort_Global_Soln_4D(KC,KCM, NSCM, Soln_Global_1D_SED, Global_Local_L_to_Global_SED, SED_Global)
!    endif
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** SEDB - Bed cohesive sediments
!    Soln_Local_1D_SEDB            = 0.0
!    Soln_Global_1D_SEDB           = 0.0
!    Map_Local_L_to_Global_SEDB    = 0
!    Global_Local_L_to_Global_SEDB = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KB, KBM, NSED, NSCM, SEDB, Soln_Local_1D_SEDB, &
!                     Map_Local_L_to_Global_SEDB)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KB, NSED, Soln_Local_1D_SEDB, Soln_Global_1D_SEDB, &
!                   Map_Local_L_to_Global_SEDB, Global_Local_L_to_Global_SEDB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_4D(KB, KBM, NSCM, Soln_Global_1D_SEDB, Global_Local_L_to_Global_SEDB, SEDB_Global)
!    end if
!  ENDIF
!
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** If SND is turned on then we should map to global domain so we can write it out
!  IF( ISTRAN(7) >= 1  )THEN
!    Soln_Local_1D_SND            = 0.0
!    Soln_Global_1D_SND           = 0.0
!    Map_Local_L_to_Global_SND    = 0
!    Global_Local_L_to_Global_SND = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NSND, NSNM, SND, Soln_Local_1D_SND,  &
!                     Map_Local_L_to_Global_SND)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NSNM, Soln_Local_1D_SND, Soln_Global_1D_SND, &
!                   Map_Local_L_to_Global_SND, Global_Local_L_to_Global_SND)
!
!    ! *** Sorting step onto entire spatial domain
!    if( process_id == master_id )THEN
!      Call Sort_Global_Soln_4D(KC,KCM, NSNM, Soln_Global_1D_SND, Global_Local_L_to_Global_SND, SND_Global)
!    endif
!    
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** SNDB - Bed cohesive sediments
!    Soln_Local_1D_SNDB            = 0.0
!    Soln_Global_1D_SNDB           = 0.0
!    Map_Local_L_to_Global_SNDB    = 0
!    Global_Local_L_to_Global_SNDB = 0
!
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KB, KBM, NSND, NSND, SNDB, Soln_Local_1D_SNDB, &
!                     Map_Local_L_to_Global_SNDB)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KB, NSND, Soln_Local_1D_SNDB, Soln_Global_1D_SNDB, &
!                   Map_Local_L_to_Global_SNDB, Global_Local_L_to_Global_SNDB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_4D(KB, KBM, NSNM, Soln_Global_1D_SNDB, Global_Local_L_to_Global_SNDB, SNDB_Global)
!    end if
!
!  ENDIF
!
!  IF( LSEDZLJ )THEN
!    if( process_id == master_id  )THEN
!      DO L=2,LA_Global
!        TAU_Global(L) = SHEAR_Global(L)*10.*WATERDENS*1000.   ! *** Convert to dynes/cm2
!      ENDDO  
!    endif
!  
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** Gather the calculated D50AVG
!    Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, D50AVG, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!    Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                   Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!    if( process_id == master_id )THEN
!      ! *** Sort the global array so L indexing is correct
!      Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, D50AVG_Global)
!    endif
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** Gather the calculated ETOTO
!    Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, ETOTO, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!    Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                   Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!    if( process_id == master_id )THEN
!      ! *** Sort the global array so L indexing is correct
!      Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, ETOTO_Global)
!    endif
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** Gather the calculated DEPO
!    Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, DEPO, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!
!    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!    Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!                   Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!
!    if( process_id == master_id )THEN
!      ! *** Sort the global array so L indexing is correct
!      Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, DEPO_Global)
!    endif
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** TSED - Total sediment bed mass by layer
!    Soln_Local_1D_KB            = 0.0
!    Soln_Global_1D_KB           = 0.0
!    Map_Local_L_to_Global_KB    = 0
!    Global_Local_L_to_Global_KB = 0
!    TMPBED = 0.
!    
!    ! *** Reverse Order of Array
!    DO L=2,LA
!      DO K=1,KB
!        TMPBED(L,K) = TSED(K,L)
!      ENDDO
!    ENDDO
!    
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KB, TMPBED, Soln_Local_1D_KB, Map_Local_L_to_Global_KB)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!                   Map_Local_L_to_Global_KB, Global_Local_L_to_Global_KB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KB, Soln_Global_1D_KB, Global_Local_L_to_Global_KB, TSED_Global)
!    end if
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** BULKDENS - Sediment denisty by layer
!    Soln_Local_1D_KB            = 0.0
!    Soln_Global_1D_KB           = 0.0
!    Map_Local_L_to_Global_KB    = 0
!    Global_Local_L_to_Global_KB = 0
!    TMPBED = 0.
!    
!    ! *** Reverse Order of Array
!    DO L=2,LA
!      DO K=1,KB
!        TMPBED(L,K) = BULKDENS(K,L)
!      ENDDO
!    ENDDO
!    
!    Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KB, TMPBED, Soln_Local_1D_KB, Map_Local_L_to_Global_KB)
!
!    ! *** Send local values containing values ommitting the ghost cells to master
!    Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!                   Map_Local_L_to_Global_KB, Global_Local_L_to_Global_KB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_3D(KB, Soln_Global_1D_KB, Global_Local_L_to_Global_KB, BULKDENS_Global)
!    end if
!
!    !---------------------------------------------------------------------------------------------------------!
!    ! *** PERSED - Sediment mass fraction by layer
!    Soln_Local_1D_SNDB            = 0.0
!    Soln_Global_1D_SNDB           = 0.0
!    Map_Local_L_to_Global_SNDB    = 0
!    Global_Local_L_to_Global_SNDB = 0
!    BEDBYCLASS = 0.
!    
!    ! *** Reverse Order of Array
!    DO L=2,LA
!      DO K=1,KB
!        DO NS=1,NSCM
!          BEDBYCLASS(L,K,NS) = PERSED(NS,K,L)
!        ENDDO
!      ENDDO
!    ENDDO
!    
!    Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KB, KBM, NSED, NSCM, BEDBYCLASS, Soln_Local_1D_SEDB, &
!                     Map_Local_L_to_Global_SEDB)
!
!    ! *** Gather all 1D variables onto single array on master process
!    Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KB, NSED, Soln_Local_1D_SEDB, Soln_Global_1D_SEDB, &
!                   Map_Local_L_to_Global_SEDB, Global_Local_L_to_Global_SEDB)
!
!    ! *** Handle mapping of global values to the master
!    if( process_id == master_id  )THEN
!      Call Sort_Global_Soln_4D(KB, KBM, NSCM, Soln_Global_1D_SEDB, Global_Local_L_to_Global_SEDB, PERSED_Global)
!    end if
!
!  ENDIF
!  
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** For EE_ARRAYS.OUT
!  ! *** Handling QQ seperately because it stores things in the zero index.... -_-
!
!  !  Call Map3D_to_1D_QQ(LA_Local_no_ghost, num_active_l_local, KC, KCM, QQ, Soln_Local_1D_QQ, Map_Local_L_to_Global_QQ)
!=======
!  !Subroutine Global_Mapping_Routines(LA_Local_no_ghost)
!  !
!  !Use Global
!  !Use MPI
!  !Use Variables_MPI
!  !Use Variables_MPI_Concentration
!  !Use Variables_MPI_Mapping
!  !Use Variables_MPI_Write_Out
!  !
!  !Implicit None
!  !
!  !! *** Passed in
!  !Integer, intent(in) :: LA_Local_no_ghost
!  !
!  !! *** Local variables
!  !Integer :: i, j, ierr, L, k, md, ii
!  !Integer :: num_active_l_local
!  !
!  !! *** L array variables (real)
!  !Real, Allocatable, Dimension(:)    :: Soln_Local_1D_L
!  !Real, Allocatable, Dimension(:)    :: Soln_Global_1D_L
!  !Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_L
!  !Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_L
!  !
!  !! *** L array variables (integer)
!  !Integer, Allocatable, Dimension(:) :: Soln_Local_1D_L_I
!  !Integer, Allocatable, Dimension(:) :: Soln_Global_1D_L_I
!  !
!  !! *** L,KC array variables
!  !Real, Allocatable, Dimension(:)    :: Soln_Local_1D_KC
!  !Real, Allocatable, Dimension(:)    :: Soln_Global_1D_KC
!  !Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_KC
!  !Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_KC
!  !
!  !! *** L,KB array variables
!  !Real, Allocatable, Dimension(:) :: Soln_Local_1D_KB
!  !Real, Allocatable, Dimension(:) :: Soln_Global_1D_KB
!  !
!  !! *** L,KC,NDYE array variables (by DYE class)
!  !Real, Allocatable, Dimension(:) :: Soln_Local_1D_DYE
!  !Real, Allocatable, Dimension(:) :: Soln_Global_1D_DYE
!  !Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_DYE
!  !Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_DYE
!  !
!  !! *** L,KC,NTOX array variables (by TOX class)
!  !Real, Allocatable, Dimension(:) :: Soln_Local_1D_TOX
!  !Real, Allocatable, Dimension(:) :: Soln_Global_1D_TOX
!  !Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_TOX
!  !Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_TOX
!  !
!  !! *** L,KC,NSED array variables (by SED class)
!  !Real, Allocatable, Dimension(:) :: Soln_Local_1D_SED
!  !Real, Allocatable, Dimension(:) :: Soln_Global_1D_SED
!  !Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_SED
!  !Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_SED
!  !
!  !! *** L,KC,NSND array variables (by SND class)
!  !Real, Allocatable, Dimension(:) :: Soln_Local_1D_SND
!  !Real, Allocatable, Dimension(:) :: Soln_Global_1D_SND
!  !Integer, Allocatable, Dimension(:) :: Map_Local_L_to_Global_SND
!  !Integer, Allocatable, Dimension(:) :: Global_Local_L_to_Global_SND
!  !
!  !! *** L,0:KC array variables
!  !Real, Allocatable, Dimension(:)    ::Soln_Local_1D_QQ
!  !Real, Allocatable, Dimension(:)    ::Soln_Global_1D_QQ
!  !Integer, Allocatable, Dimension(:) ::Map_Local_L_to_Global_QQ
!  !Integer, Allocatable, Dimension(:) ::Global_Local_L_to_Global_QQ
!  !!---------------------------------------------------------------------------!
!  !
!  !If(.not.Allocated(Soln_Local_1D_KC) )THEN
!  !
!  !  Allocate(Soln_Local_1D_L(LA_Local_no_ghost))
!  !  Allocate(Soln_Global_1D_L(LCM_Global))
!  !  Allocate(Map_Local_L_to_Global_L(LA_Local_no_ghost ))
!  !  Allocate(Global_Local_L_to_Global_L(LCM_Global))
!  !
!  !  Allocate(Soln_Local_1D_L_I(LA_Local_no_ghost))
!  !  Allocate(Soln_Global_1D_L_I(LCM_Global))
!  !
!  !  !Allocate(Soln_Local_1D_KC(LA_Local_no_ghost*kc))
!  !  !Allocate(Soln_Global_1D_KC(LCM_Global*kc))
!  !  !Allocate(Map_Local_L_to_Global_KC(LA_local_no_ghost*KC))
!  !  !Allocate(Global_Local_L_to_Global_KC(LCM_Global*KC))
!  !
!  !  Allocate(Soln_Local_1D_KC(LA_Local_no_ghost*kcm))
!  !  Allocate(Soln_Global_1D_KC(LCM_Global*kcm))
!  !  Allocate(Map_Local_L_to_Global_KC(LA_local_no_ghost*kcm))
!  !  Allocate(Global_Local_L_to_Global_KC(LCM_Global*kcm))
!  !
!  !  Allocate(Soln_Local_1D_KB(LA_Local_no_ghost*KB))
!  !  Allocate(Soln_Global_1D_KB(LCM_Global*KB))
!  !
!  !  IF( ISTRAN(3) > 0 )THEN
!  !    Allocate(Soln_Local_1D_DYE(La_local_no_ghost*KC*NDYE))
!  !    Allocate(Soln_Global_1D_DYE(LCM_Global*KC*NDYM))
!  !    Allocate(Map_Local_L_to_Global_DYE(LA_local_no_ghost*KC*NDYM))
!  !    Allocate(Global_Local_L_to_Global_DYE(LCM_Global*KC*NDYM))
!  !  ENDIF
!  !
!  !  IF( ISTRAN(5) > 0 )THEN
!  !    Allocate(Soln_Local_1D_TOX(La_local_no_ghost*KC*NTOX))
!  !    Allocate(Soln_Global_1D_TOX(LCM_Global*KC*NTXM))
!  !    Allocate(Map_Local_L_to_Global_TOX(LA_local_no_ghost*KC*NTXM))
!  !    Allocate(Global_Local_L_to_Global_TOX(LCM_Global*KC*NTXM))
!  !  ENDIF
!  !
!  !  IF( ISTRAN(6) > 0 )THEN
!  !    Allocate(Soln_Local_1D_SED(La_local_no_ghost*KC*NSCM))
!  !    Allocate(Soln_Global_1D_SED(LCM_Global*KC*NSCM))
!  !    Allocate(Map_Local_L_to_Global_SED(LA_local_no_ghost*KC*NSCM))
!  !    Allocate(Global_Local_L_to_Global_SED(LCM_Global*KC*NSCM))
!  !  ENDIF
!  !
!  !  IF( ISTRAN(7) > 0 )THEN
!  !    Allocate(Soln_Local_1D_SND(La_local_no_ghost*KC*NSNM))
!  !    Allocate(Soln_Global_1D_SND(LCM_Global*KC*NSNM))
!  !    Allocate(Map_Local_L_to_Global_SND(LA_local_no_ghost*KC*NSNM))
!  !    Allocate(Global_Local_L_to_Global_SND(LCM_Global*KC*NSNM))
!  !  ENDIF
!  !
!  !  Allocate(Soln_Local_1D_QQ(LA_Local_no_ghost*kc+1)) ! + 1 because QQ goes from 0 to KCM
!  !  Allocate(Soln_Global_1D_QQ(LCM_Global*kc+1))
!  !  Allocate(Map_Local_L_to_Global_QQ(LA_local_no_ghost*KC+1))
!  !  Allocate(Global_Local_L_to_Global_QQ(LCM_Global*KC+1))
!  !end if
!  !
!  !!Soln_Local_1D_KB   = 0.0
!  !!Soln_Global_1D_KB  = 0.0
!  !
!  !!Soln_Local_1D_QQ  = 0.0
!  !!Soln_Global_1D_QQ = 0.0
!  !!Map_Local_L_to_Global_QQ = 0
!  !!Global_Local_L_to_Global_QQ = 0
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** HP
!  !! *** Map to 1d array containing l values excluding ghost cells and build up array mapping to global l
!  !Soln_Local_1D_L            = 0.0
!  !Soln_Global_1D_L           = 0.0
!  !Map_Local_L_to_Global_L    = 0
!  !Global_Local_L_to_Global_L = 0
!  !
!  !Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, HP, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!  !
!  !! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  !Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!  !               Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!  !
!  !
!  !if( process_id == master_id )THEN
!  !  ! *** Sort the global array so L indexing is correct
!  !  Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, HP_Global)
!  !endif
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** U
!  !! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !Soln_Local_1D_KC            = 0.0
!  !Soln_Global_1D_KC           = 0.0
!  !Map_Local_L_to_Global_KC    = 0
!  !Global_Local_L_to_Global_KC = 0
!  !
!  !Call Map3D_to_1D(la_local_no_ghost, num_active_l_local, KCM, U, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!  !
!  !! *** Send local values containing values ommitting the ghost cells to master
!  !Call Gather_3D(la_local_no_ghost, num_active_l_local, LCM, KCM, Soln_Local_1D_KC, Soln_Global_1D_KC, &
!  !               Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !! *** Handle mapping of global values to the master
!  !if( process_id == 0  )THEN
!  !  Call Sort_Global_Soln_3D(KCM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, U_Global)
!  !end if
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** V
!  !! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !Soln_Local_1D_KC            = 0.0
!  !Soln_Global_1D_KC           = 0.0
!  !Map_Local_L_to_Global_KC    = 0
!  !Global_Local_L_to_Global_KC = 0
!  !
!  !Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local,KCM, V, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!  !
!  !! *** Send local values containing values ommitting the ghost cells to master
!  !Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KCM, Soln_Local_1D_KC, Soln_Global_1D_KC, &
!  !               Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !! *** Handle mapping of global values to the master
!  !if( process_id == 0  )THEN
!  !  Call Sort_Global_Soln_3D(KCM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, V_Global)
!  !end if
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** W
!  !! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !Soln_Local_1D_KC            = 0.0
!  !Soln_Global_1D_KC           = 0.0
!  !Map_Local_L_to_Global_KC    = 0
!  !Global_Local_L_to_Global_KC = 0
!  !
!  !Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KCM, W, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!  !
!  !! *** Send local values containing values ommitting the ghost cells to master
!  !Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KCM, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!  !               Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !! *** Handle mapping of global values to the master
!  !if( process_id == 0  )THEN
!  !  Call Sort_Global_Soln_3D(KCM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, W_Global)
!  !end if
!  !
!  !---------------------------------------------------------------------------------------------------------!
!  ! *** Bed Shear
!  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !Soln_Local_1D_L            = 0.0
!  !Soln_Global_1D_L           = 0.0
!  !Map_Local_L_to_Global_L    = 0
!  !Global_Local_L_to_Global_L = 0
!  !
!    
! 
!  
!  ! *** Now gather the calculated bed shears
!  !Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, SHEAR_Local, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!  !
!  !! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  !Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!  !               Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!  !
!  !if( process_id == master_id )THEN
!  !  ! *** Sort the global array so L indexing is correct
!  !  Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, SHEAR_Global)
!  !endif
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** Map Salinity
!  !IF(ISTRAN(1).EQ.1 )THEN
!  !  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !  Soln_Local_1D_KC            = 0.0
!  !  Soln_Global_1D_KC           = 0.0
!  !  Map_Local_L_to_Global_KC    = 0
!  !  Global_Local_L_to_Global_KC = 0
!  !
!  !  Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local,KCM, SAL, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!>>>>>>> gather_arrays_ee_binary
!  !
!  !  ! *** Send local values containing values ommitting the ghost cells to master
!  !  Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KCM, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!  !                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !  ! *** Handle mapping of global values to the master
!<<<<<<< HEAD
!  !  if( process_id == master_id  )THEN
!  !      Call Sort_Global_Soln_3D_QQ(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, QQ_Global)
!  !end if
!
!  !***end QQ
!  !---------------------------------------------------------------------------!
!
!  !  ! *** For EE_BC.OUT
!=======
!  !  if( process_id == 0  )THEN
!  !    Call Sort_Global_Soln_3D(KCM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, SAL_Global)
!  !  end if
!  !End if
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** Map Temperature
!  !IF(ISTRAN(2).EQ.1 )THEN
!  !  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !  Soln_Local_1D_KC            = 0.0
!  !  Soln_Global_1D_KC           = 0.0
!  !  Map_Local_L_to_Global_KC    = 0
!  !  Global_Local_L_to_Global_KC = 0
!  !
!  !  Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local,KCM, TEM, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!  !
!  !  ! *** Send local values containing values ommitting the ghost cells to master
!  !  Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KCM, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!  !                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !  ! *** Handle mapping of global values to the master
!  !  if( process_id == 0  )THEN
!  !    Call Sort_Global_Soln_3D(KCM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, TEM_Global)
!  !  end if
!  !
!  !  IF( TBEDIT > 0.)THEN
!  !    ! *** Bed Temperature
!  !    Soln_Local_1D_L  = 0.0
!  !    Soln_Global_1D_L = 0.0
!  !    Map_Local_L_to_Global_L    = 0
!  !    Global_Local_L_to_Global_L = 0
!  !
!  !    Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, TEMB, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!  !
!  !    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  !    Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!  !                   Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!  !
!  !    if( process_id == master_id )THEN
!  !      ! *** Sort the global array so L indexing is correct
!  !      Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, TEMB_Global)
!  !    endif
!  !  ENDIF
!  !
!  !End if
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** If dye is turned on then we should map to global domain so we can write it out
!  !IF( ISTRAN(3) >= 1  )THEN
!  !  Soln_Local_1D_DYE            = 0.0
!  !  Soln_Global_1D_DYE           = 0.0
!  !  Map_Local_L_to_Global_DYE    = 0
!  !  Global_Local_L_to_Global_DYE = 0
!  !
!  !  Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NDYE, NDYM, DYE, Soln_Local_1D_DYE, &
!  !                   Map_Local_L_to_Global_DYE)
!  !       
!  !  ! *** Gather all 1D variables onto single array on master process
!  !  Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NDYM, Soln_Local_1D_DYE, Soln_Global_1D_DYE, &
!  !                 Map_Local_L_to_Global_DYE, Global_Local_L_to_Global_DYE)
!  !
!  !  ! *** Sorting step onto entire spatial domain
!  !  if( process_id == 0 )THEN
!  !    Call Sort_Global_Soln_4D(KC,KCM, NDYM, Soln_Global_1D_DYE, Global_Local_L_to_Global_DYE, DYE_Global)
!  !  endif
!  !
!  !endif
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** If SFL is turned on then we should map to global domain so we can write it out
!  !IF( ISTRAN(4) >= 1  )THEN
!  !  ! *** Setup local arrays containing solution without ghost cells and a corresponding global mapping
!  !  Soln_Local_1D_KC            = 0.0
!  !  Soln_Global_1D_KC           = 0.0
!  !  Map_Local_L_to_Global_KC    = 0
!  !  Global_Local_L_to_Global_KC = 0
!  !
!  !  Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local,KCM, SFL, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!  !
!  !  ! *** Send local values containing values ommitting the ghost cells to master
!  !  Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KCM, Soln_Local_1D_KC,   Soln_Global_1D_KC, &
!  !                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !  ! *** Handle mapping of global values to the master
!  !  if( process_id == 0  )THEN
!  !    Call Sort_Global_Soln_3D(KCM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, SFL_Global)
!  !  end if
!  !
!  !endif
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** If TOX is turned on then we should map to global domain so we can write it out
!  !IF( ISTRAN(5) >= 1  )THEN
!  !  Soln_Local_1D_TOX            = 0.0
!  !  Soln_Global_1D_TOX           = 0.0
!  !  Map_Local_L_to_Global_TOX    = 0
!  !  Global_Local_L_to_Global_TOX = 0
!  !
!  !  Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NTOX, NTOX, TOX, Soln_Local_1D_TOX, &
!  !                   Map_Local_L_to_Global_TOX)
!  !
!  !  ! *** Gather all 1D variables onto single array on master process
!  !  Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NTOX, Soln_Local_1D_TOX, Soln_Global_1D_TOX, &
!  !                 Map_Local_L_to_Global_TOX, Global_Local_L_to_Global_TOX)
!  !
!  !  ! *** Sorting step onto entire spatial domain
!  !  if( process_id == 0 )THEN
!  !    Call Sort_Global_Soln_4D(KC,KCM, NDYM, Soln_Global_1D_TOX, Global_Local_L_to_Global_TOX, TOX_Global)
!  !  endif
!  !
!  !  ! DELME - TODO - BED
!  !ENDIF
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** If SED or SND is turned on then we should map to global domain so we can write it out
!  !IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
!  !  Soln_Local_1D_L_I          = 0
!  !  Soln_Global_1D_L_I         = 0
!  !  Map_Local_L_to_Global_L    = 0
!  !  Global_Local_L_to_Global_L = 0
!  !
!  !  ! *** KBT - Sediment Bed Top
!  !  Call Map1D_no_ghost_int(la_local_no_ghost, num_active_l_local, KBT, Soln_Local_1D_L_I, Map_Local_L_to_Global_L)
!  !
!  !  ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  !  Call Gather_1D_int(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L_I, Soln_Global_1D_L, &
!  !                 Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!  !
!  !  if( process_id == master_id )THEN
!  !    ! *** Sort the global array so L indexing is correct
!  !    Call Sort_Global_Soln_HP_Int(Soln_Global_1D_L, Global_Local_L_to_Global_L, KBT_Global)
!  !  endif
!  !
!  !  !---------------------------------------------------------------------------------------------------------!
!  !  IF( IMORPH > 0 )THEN
!  !    Soln_Local_1D_L            = 0.0
!  !    Soln_Global_1D_L           = 0.0
!  !    Map_Local_L_to_Global_L    = 0
!  !    Global_Local_L_to_Global_L = 0
!  !
!  !    Call Map1D_no_ghost(la_local_no_ghost, num_active_l_local, BELV, Soln_Local_1D_L, Map_Local_L_to_Global_L)
!  !
!  !    ! *** Calling gather 3d but just specify a value of 1 for the second variable so that we
!  !    Call Gather_1D(la_local_no_ghost,  num_active_l_local, Soln_Local_1D_L, Soln_Global_1D_L, &
!  !                   Map_Local_L_to_Global_L, Global_Local_L_to_Global_L)
!  !
!  !
!  !    if( process_id == master_id )THEN
!  !      ! *** Sort the global array so L indexing is correct
!  !      Call Sort_Global_Soln_HP(Soln_Global_1D_L, Global_Local_L_to_Global_L, BELV_Global)
!  !    endif
!  !  ENDIF
!  !  
!  !  !---------------------------------------------------------------------------------------------------------!
!  !  Soln_Local_1D_KB            = 0.0
!  !  Soln_Global_1D_KB           = 0.0
!  !  Map_Local_L_to_Global_KC    = 0
!  !  Global_Local_L_to_Global_KC = 0
!  !
!  !  Call Map3D_to_1D(LA_Local_no_ghost, num_active_l_local, KBM, HBED, Soln_Local_1D_KC, Map_Local_L_to_Global_KC)
!  !
!  !  ! *** Send local values containing values ommitting the ghost cells to master
!  !  Call Gather_3D(LA_Local_no_ghost, num_active_l_local, LCM, KBM, Soln_Local_1D_KB, Soln_Global_1D_KB, &
!  !                 Map_Local_L_to_Global_KC, Global_Local_L_to_Global_KC)
!  !
!  !  ! *** Handle mapping of global values to the master
!  !  if( process_id == 0  )THEN
!  !    Call Sort_Global_Soln_3D(KBM, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, HBED_Global)
!  !  end if
!  !
!  !  BDENBED_Global = -999   ! delme - todo
!  !  PORBED_Global  = -999   ! delme - todo
!  !  VFRBED_Global  = -999   ! delme - todo
!  !  IF( ISTRAN(7) > 0 .AND. ISBDLDBC > 0 .AND. NSND > 0 )THEN
!  !    QSBDLDX_Global = -999   ! delme - todo
!  !    QSBDLDY_Global = -999   ! delme - todo
!  !  ENDIF
!  !ENDIF
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** If SED is turned on then we should map to global domain so we can write it out
!  !IF( ISTRAN(6) >= 1  )THEN
!  !  Soln_Local_1D_SED            = 0.0
!  !  Soln_Global_1D_SED           = 0.0
!  !  Map_Local_L_to_Global_SED    = 0
!  !  Global_Local_L_to_Global_SED = 0
!  !
!  !  Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NSED, NSCM, SED, Soln_Local_1D_SED, &
!  !                   Map_Local_L_to_Global_SED)
!  !
!  !  ! *** Gather all 1D variables onto single array on master process
!  !  Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NSCM, Soln_Local_1D_SED, Soln_Global_1D_SED, &
!  !                 Map_Local_L_to_Global_SED, Global_Local_L_to_Global_SED)
!  !
!  !  ! *** Sorting step onto entire spatial domain
!  !  if( process_id == 0 )THEN
!  !    Call Sort_Global_Soln_4D(KC,KCM, NSCM, Soln_Global_1D_SED, Global_Local_L_to_Global_SED, SED_Global)
!  !  endif
!  !
!  !  ! *** Bed
!  !  !Soln_Local_1D_KB            = 0.0
!  !  !Soln_Global_1D_SED           = 0.0
!  !  !Map_Local_L_to_Global_SED    = 0
!  !  !Global_Local_L_to_Global_SED = 0
!  !  !
!  !  !Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NSED, NSED, SED, Soln_Local_1D_SED, &
!  !  !                 Map_Local_L_to_Global_SED)
!  !  !
!  !  !! *** Gather all 1D variables onto single array on master process
!  !  !Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NSED, Soln_Local_1D_SED, Soln_Global_1D_SED, &
!  !  !               Map_Local_L_to_Global_SED, Global_Local_L_to_Global_SED)
!  !  !
!  !  !! *** Sorting step onto entire spatial domain
!  !  !if( process_id == 0 )THEN
!  !  !  Call Sort_Global_Soln_4D(KC,KCM, NDYM, Soln_Global_1D_SED, Global_Local_L_to_Global_SED, SED_Global)
!  !  !endif  
!  !ENDIF
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** If SND is turned on then we should map to global domain so we can write it out
!  !IF( ISTRAN(7) >= 1  )THEN
!  !  Soln_Local_1D_SND            = 0.0
!  !  Soln_Global_1D_SND           = 0.0
!  !  Map_Local_L_to_Global_SND    = 0
!  !  Global_Local_L_to_Global_SND = 0
!  !
!  !  Call Map4D_to_1D(la_local_no_ghost, num_active_l_local, KC, KCM, NSND, NSNM, SND, Soln_Local_1D_SND,  &
!  !                   Map_Local_L_to_Global_SND)
!  !
!  !  ! *** Gather all 1D variables onto single array on master process
!  !  Call Gather_4D(la_local_no_ghost, num_active_l_local, LCM, KC, NSNM, Soln_Local_1D_SND, Soln_Global_1D_SND, &
!  !                 Map_Local_L_to_Global_SND, Global_Local_L_to_Global_SND)
!  !
!  !  ! *** Sorting step onto entire spatial domain
!  !  if( process_id == 0 )THEN
!  !    Call Sort_Global_Soln_4D(KC,KCM, NSNM, Soln_Global_1D_SND, Global_Local_L_to_Global_SND, SND_Global)
!  !  endif
!  !  
!  !  ! DELME - TODO - BED
!  !
!  !ENDIF
!  !
!  !!---------------------------------------------------------------------------------------------------------!
!  !! *** For EE_ARRAYS.OUT
!  !! *** Handling QQ seperately because it stores things in the zero index.... -_-
!  !
!  !!  Call Map3D_to_1D_QQ(LA_Local_no_ghost, num_active_l_local, KC, KCM, QQ, Soln_Local_1D_QQ, Map_Local_L_to_Global_QQ)
!  !!
!  !!  ! *** Send local values containing values ommitting the ghost cells to master
!  !!  Call Gather_3D_QQ(LA_Local_no_ghost, num_active_l_local, LCM, KC, &
!  !!                 Soln_Local_1D_QQ,   Soln_Global_1D_QQ, &
!  !!                 Map_Local_L_to_Global_QQ, Global_Local_L_to_Global_QQ)
!  !!
!  !!  ! *** Handle mapping of global values to the master
!  !!  if( process_id == 0  )THEN
!  !!      Call Sort_Global_Soln_3D_QQ(KC, Soln_Global_1D_KC, Global_Local_L_to_Global_KC, QQ_Global)
!  !!end if
!  !
!  !!***end QQ
!  !!---------------------------------------------------------------------------!
!  !
!  !!  ! *** For EE_BC.OUT
!>>>>>>> gather_arrays_ee_binary
!  !  If( KC > 1  )THEN
!  !      ! *** Zero out variables used in remapping process
!  !      Soln_Local_1D  = 0.0
!  !      Soln_Global_1D = 0.0
!  !      ! *** Map to 1D array not including ghost cells
!  !      Call Map3D_to_1D(KC, KCM, QSUM, Soln_Local_1D)
!  !      ! *** Gather all 1D arrays onto global array
!  !      Call Gather_3D(LCM, KC, Soln_Local_1D, Soln_Global_1D)
!  !      ! *** Map back to global (L,K) array for writing to EE binary files
!  !      Call Map1D_to_3D(KC, KCM, Soln_Global_1D, QSUM_Global)
!  !  else
!  !      Call MapLocalToGlobal_1D(QSUME, QSUME_Local_to_Global)
!  !      Call Gather_1D(QSUME_Local_to_Global, QSUME_Global)
!  !  endif
!  !
!  !  ! *** Zero out variables used in remapping process
!  !  Soln_Local_1D  = 0.0
!  !  Soln_Global_1D = 0.0
!  !  ! *** Map to 1D array not including ghost cells
!  !  Call Map3D_to_1D(KC, KCM, VHDX2, Soln_Local_1D)
!  !  ! *** Gather all 1D arrays onto global array
!  !  Call Gather_3D(LCM, KC, Soln_Local_1D, Soln_Global_1D)
!  !  ! *** Map back to global (L,K) array for writing to EE binary files
!  !  Call Map1D_to_3D(KC, KCM, Soln_Global_1D, VHDX2_Global)
!  !
!  !  ! *** Zero out variables used in remapping process
!  !  Soln_Local_1D  = 0.0
!  !  Soln_Global_1D = 0.0
!  !  ! *** Map to 1D array not including ghost cells
!  !  Call Map3D_to_1D(KC, KCM, UHDY2, Soln_Local_1D)
!  !  ! *** Gather all 1D arrays onto global array
!  !  Call Gather_3D(LCM, KC, Soln_Local_1D, Soln_Global_1D)
!  !  ! *** Map back to global (L,K) array for writing to EE binary files
!  !  Call Map1D_to_3D(KC, KCM, Soln_Global_1D, UHDY2_Global)
!  !!
!  !!
!  !!  !---------------------------------------------------------------------------!
!  !!  ! *** For EE_BED.OUT
!  !!  Call MapLocalToGlobal_1D(KBT, KBT_Local_to_Global)
!  !!  Call Gather_1D(KBT_Local_to_Global, KBT_Global)
!  !!
!  !!  ! *** Zero out variables used in remapping process
!  !!  Soln_Local_1D_KB   = 0.0
!  !!  Soln_Global_1D_KB  = 0.0
!  !!  ! *** Map to 1D array not including ghost cells
!  !!
!  !!  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 .OR. ISBAL > 0 )THEN
!  !!      Call Map3D_to_1D(KB, KBM, HBED, Soln_Local_1D_KB)
!  !!
!  !!      ! *** Gather all 1D arrays onto global array
!  !!      Call Gather_3D(LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB)
!  !!      ! *** Map back to global (L,K) array for writing to EE binary files
!  !!      Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, HBED_Global)
!  !!
!  !!     ! *** Zero out variables used in remapping process
!  !!      Soln_Local_1D_KB   = 0.0
!  !!      Soln_Global_1D_KB  = 0.0
!  !!     ! *** Map to 1D array not including ghost cells
!  !!     Call Map3D_to_1D(KB, KBM, BDENBED, Soln_Local_1D_KB)
!  !!     ! *** Gather all 1D arrays onto global array
!  !!     Call Gather_3D(LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB)
!  !!     ! *** Map back to global (L,K) array for writing to EE binary files
!  !!     Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, BDENBED_Global)
!  !!     ! *** Zero out variables used in remapping process
!  !!     Soln_Local_1D_KB   = 0.0
!  !!     Soln_Global_1D_KB  = 0.0
!  !!
!  !!     ! *** Map to 1D array not including ghost cells
!  !!     Call Map3D_to_1D(KB, KBM, PORBED, Soln_Local_1D_KB)
!  !!     ! *** Gather all 1D arrays onto global array
!  !!     Call Gather_3D(LCM, KB, Soln_Local_1D, Soln_Global_1D_KB)
!  !!     ! *** Map back to global (L,K) array for writing to EE binary files
!  !!     Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, PORBED_Global)
!  !!     ! *** Zero out variables used in remapping process
!  !!     !Soln_Local_1D_KB   = 0.0
!  !!     !Soln_Global_1D_KB  = 0.0
!  !!
!  !!     ! *** 4D
!  !!
!  !!     !!***Map to 1D array not including ghost cells
!  !!     !Call Map3D_to_1D(KB, KBM, SEDB, Soln_Local_1D_KB)
!  !!     !!***Gather all 1D arrays onto global array
!  !!     !Call Gather_3D(LCM, KB, Soln_Local_1D, Soln_Global_1D_KB)
!  !!     !!***Map back to global (L,K) array for writing to EE binary files
!  !!     !Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, SEDB_Global)
!  !!
!  !!
!  !!     Call MapLocalToGlobal_4D(KB, KBM, NSCM, SEDB, SEDB_Local_to_Global)
!  !!     Call Gather_4D(KBM, NSCM, SEDB_Local_to_Global, SEDB_Global)
!  !!
!  !!     Call MapLocalToGlobal_4D(KB, KBM, NSNM, SNDB, SNDB_Local_to_Global)
!  !!     Call Gather_4D(KBM, NSNM, SNDB_Local_to_Global, SNDB_Global)
!  !!
!  !!     Call MapLocalToGlobal_4D(KB, KBM, NTXM, TOXB, TOXB_Local_to_Global)
!  !!     Call Gather_4D(KBM, NTXM, TOXB_Local_to_Global, TOXB_Global)
!  !! End if
!  !!
!  !
!  !Deallocate(Soln_Local_1D_L)
!  !Deallocate(Soln_Global_1D_L)
!  !Deallocate(Map_Local_L_to_Global_L)
!  !Deallocate(Global_Local_L_to_Global_L)
!  !Deallocate(Soln_Local_1D_KC)
!  !Deallocate(Soln_Global_1D_KC)
!  !Deallocate(Map_Local_L_to_Global_KC)
!  !Deallocate(Global_Local_L_to_Global_KC)
!  !Deallocate(Soln_Local_1D_KB)
!  !Deallocate(Soln_Global_1D_KB)
!  !
!<<<<<<< HEAD
!  !  !---------------------------------------------------------------------------!
!  !  ! *** For EE_BED.OUT
!  !  Call MapLocalToGlobal_1D(KBT, KBT_Local_to_Global)
!  !  Call Gather_1D(KBT_Local_to_Global, KBT_Global)
!  !
!  !  ! *** Zero out variables used in remapping process
!  !  Soln_Local_1D_KB   = 0.0
!  !  Soln_Global_1D_KB  = 0.0
!  !  ! *** Map to 1D array not including ghost cells
!  !
!  !  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 .OR. ISBAL > 0 )THEN
!  !      Call Map3D_to_1D(KB, KBM, HBED, Soln_Local_1D_KB)
!  !
!  !      ! *** Gather all 1D arrays onto global array
!  !      Call Gather_3D(LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB)
!  !      ! *** Map back to global (L,K) array for writing to EE binary files
!  !      Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, HBED_Global)
!  !
!  !     ! *** Zero out variables used in remapping process
!  !      Soln_Local_1D_KB   = 0.0
!  !      Soln_Global_1D_KB  = 0.0
!  !     ! *** Map to 1D array not including ghost cells
!  !     Call Map3D_to_1D(KB, KBM, BDENBED, Soln_Local_1D_KB)
!  !     ! *** Gather all 1D arrays onto global array
!  !     Call Gather_3D(LCM, KB, Soln_Local_1D_KB, Soln_Global_1D_KB)
!  !     ! *** Map back to global (L,K) array for writing to EE binary files
!  !     Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, BDENBED_Global)
!  !     ! *** Zero out variables used in remapping process
!  !     Soln_Local_1D_KB   = 0.0
!  !     Soln_Global_1D_KB  = 0.0
!  !
!  !     ! *** Map to 1D array not including ghost cells
!  !     Call Map3D_to_1D(KB, KBM, PORBED, Soln_Local_1D_KB)
!  !     ! *** Gather all 1D arrays onto global array
!  !     Call Gather_3D(LCM, KB, Soln_Local_1D, Soln_Global_1D_KB)
!  !     ! *** Map back to global (L,K) array for writing to EE binary files
!  !     Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, PORBED_Global)
!  !     ! *** Zero out variables used in remapping process
!  !     !Soln_Local_1D_KB   = 0.0
!  !     !Soln_Global_1D_KB  = 0.0
!  !
!  !     ! *** 4D
!  !
!  !     !!***Map to 1D array not including ghost cells
!  !     !Call Map3D_to_1D(KB, KBM, SEDB, Soln_Local_1D_KB)
!  !     !!***Gather all 1D arrays onto global array
!  !     !Call Gather_3D(LCM, KB, Soln_Local_1D, Soln_Global_1D_KB)
!  !     !!***Map back to global (L,K) array for writing to EE binary files
!  !     !Call Map1D_to_3D(KB, KBM, Soln_Global_1D_KB, SEDB_Global)
!  !
!  !
!  !     Call MapLocalToGlobal_4D(KB, KBM, NSCM, SEDB, SEDB_Local_to_Global)
!  !     Call Gather_4D(KBM, NSCM, SEDB_Local_to_Global, SEDB_Global)
!  !
!  !     Call MapLocalToGlobal_4D(KB, KBM, NSNM, SNDB, SNDB_Local_to_Global)
!  !     Call Gather_4D(KBM, NSNM, SNDB_Local_to_Global, SNDB_Global)
!  !
!  !     Call MapLocalToGlobal_4D(KB, KBM, NTXM, TOXB, TOXB_Local_to_Global)
!  !     Call Gather_4D(KBM, NTXM, TOXB_Local_to_Global, TOXB_Global)
!  ! End if
!  !
!
!  Deallocate(Soln_Local_1D_L)
!  Deallocate(Soln_Global_1D_L)
!  Deallocate(Map_Local_L_to_Global_L)
!  Deallocate(Global_Local_L_to_Global_L)
!  Deallocate(Soln_Local_1D_KC)
!  Deallocate(Soln_Global_1D_KC)
!  Deallocate(Map_Local_L_to_Global_KC)
!  Deallocate(Global_Local_L_to_Global_KC)
!  Deallocate(Soln_Local_1D_KB)
!  Deallocate(Soln_Global_1D_KB)
!
!End subroutine Global_Mapping_Routines
!=======
!  !End subroutine Global_Mapping_Routines
!>>>>>>> gather_arrays_ee_binary
