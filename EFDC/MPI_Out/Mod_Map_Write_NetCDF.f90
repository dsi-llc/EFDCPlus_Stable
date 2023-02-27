! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------

    !---------------------------------------------------------------------------!
    !> @details Calls the mapping routines for the NetCDF output
    !
    !> @author Zander Mausolff
    !---------------------------------------------------------------------------!
    module Mod_Map_Write_NetCDF

    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Mapping
    use Variables_MPI_Write_Out
    use Variables_MPI_MapGatherSort
    use DRIFTER, only : Gather_Drifter_Arrays
    use Mod_Map_Write_EE_Binary
    use Mod_Assign_Loc_Glob_For_Write

    implicit none

    contains

    !---------------------------------------------------------------------------!
    !> @details Maps and writes out additional variables needed for NetCDF output
    !
    !> @author Zander Mausolff
    !---------------------------------------------------------------------------!
    subroutine Map_Write_NetCDF

    implicit none
    ! *** Dummy variables

    ! *** Local variables
    integer :: j, l, k, nn, lg, m
    integer :: num_arrays_to_write_out
    integer,   Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_ActLay
    real(rkd), Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_TSED
    real(rkd), Target, Allocatable, Dimension(:,:)   :: Reverse_Temp_2D_BULKDENS
    real(rkd), Target, Allocatable, Dimension(:,:,:) :: Reverse_Temp_3D

    integer,   Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_ActLay
    real(rkd), Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_TSED
    real(rkd), Target, Allocatable, Dimension(:,:)   :: Gl_Reverse_Temp_2D_BULKDENS
    real(rkd), Target, Allocatable, Dimension(:,:,:) :: Gl_Reverse_Temp_3D


    INTEGER,TARGET,ALLOCATABLE,DIMENSION(:) :: LEC_temp_gl
    INTEGER,TARGET,ALLOCATABLE,DIMENSION(:) :: LNC_temp_gl
    
    Allocate(LEC_temp_gl(LCM))
    Allocate(LNC_temp_gl(LCM))
    

    Allocate(Reverse_Temp_2D_ActLay(LCM,KB))
    Allocate(Reverse_Temp_2D_TSED(LCM,KB))
    Allocate(Reverse_Temp_2D_BULKDENS(LCM,KB))
    Allocate(Reverse_Temp_3D(LCM,KB,NSCM))

    Allocate(GL_Reverse_Temp_2D_ActLay(LCM_Global,KB))
    Allocate(Gl_Reverse_Temp_2D_TSED(LCM_Global,KB))
    Allocate(Gl_Reverse_Temp_2D_BULKDENS(LCM_Global,KB))
    Allocate(GL_Reverse_Temp_3D(LCM_Global,KB,NSCM))
    
    LEC_temp_gl = 0
    LNC_temp_gl = 0
    
    Reverse_Temp_2D_ActLay      = 0
    Reverse_Temp_2D_TSED        = 0.0
    Reverse_Temp_2D_BULKDENS    = 0.0
    Reverse_Temp_3D             = 0.0
    Gl_Reverse_Temp_2D_ActLay   = 0
    Gl_Reverse_Temp_2D_TSED     = 0.0
    Gl_Reverse_Temp_2D_BULKDENS = 0.0
    Gl_Reverse_Temp_3D          = 0.0

    j = 0

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(CUE,1), CUE, size(CUE_Global, 1), CUE_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(CUN,1), CUN, size(CUN_Global, 1), CUN_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(CVN,1), CVN, size(CVN_Global, 1), CVN_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(CVE,1), CVE, size(CVE_Global, 1), CVE_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(DZC, 1),size(DZC,2), DZC, size(DZC_Global,1), size(DZC_Global,2), DZC_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(BELV,1), BELV, size(BELV_Global,1), BELV_Global)

    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(HP,1), HP, size(HP_Global,1), HP_Global)
    
    ! *** Point to the first array
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(U,1), size(U,2), U, size(U_Global,1), size(U_Global,2), U_Global)
    
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(V,1), size(V,2), V, size(V_Global,1), size(V_Global,2), V_Global)
    
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(W,1), size(W,2), W, size(W_Global,1), size(W_Global,2), W_Global)
    
    j = j + 1
    ! *** need to convert local L to global L 
    do m = 1, LA
        LEC_temp_gl(m) = Map2Global(LEC(m)).LG
    end do
    
    Call Assign_Loc_Glob_For_Write(j, size(LEC_temp_gl,1), LEC_temp_gl, size(LEC_Global, 1), LEC_Global) 
    
    j = j + 1
    ! *** need to convert local L to global L 
    do m = 1, LA
        LNC_temp_gl(m) = Map2Global(LNC(m)).LG
    end do
    
    Call Assign_Loc_Glob_For_Write(j, size(LNC_temp_gl,1), LNC_temp_gl, size(LNC_Global, 1), LNC_Global) 
    
    !---------------------------------------------------------
    ! *** Shear NC_WRITE_SHEAR(NC_ID, NTI)
    IF( ISNCDF(10) == 1 )then

        IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
            IF( LSEDZLJ )THEN

                j = j + 1
                Call Assign_Loc_Glob_For_Write(j, size(TAU,1), TAU, size(TAU_Global,1), TAU_Global)

            ELSEIF( ISBEDSTR >= 1 )THEN

                j = j + 1
                Call Assign_Loc_Glob_For_Write(j,size(TAUBSED,1), TAUBSED, size(TAUBSED_Global,1), TAUBSED_Global)

                IF( ISBEDSTR == 1 )THEN
                    j = j + 1
                    Call Assign_Loc_Glob_For_Write(j, size(TAUBSND,1), TAUBSND, size(TAUBSND,1), TAUBSND_Global)
                ENDIF
            ELSE
                j = j + 1
                Call Assign_Loc_Glob_For_Write(j, size(TAUB,1), TAUB, size(TAUB_Global,1), TAUB_Global)
            ENDIF
        ELSE
            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(QQ,1), size(QQ,2), QQ, size(QQ_Global,1), size(QQ_Global,2), QQ_Global)
        ENDIF

    end if
    !---------------------------------------------------------
    ! *** NC_WRITE_WIND(NC_ID, NTI)
    IF( NCDFOUT > 0 .OR. HFREOUT > 0 .OR. ISNCDF(11) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(WNDVELE,1), WNDVELE, size(WNDVELE_Global,1), WNDVELE_Global)

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(WNDVELN,1), WNDVELN, size(WNDVELN_Global,1), WNDVELN_Global)
        
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(PATMT,1), PATMT, size(PATMT_Global,1), PATMT_Global)
    end if

    !---------------------------------------------------------
    ! *** NC_WRITE_WAVE(NC_ID, NTI)
    IF( ISNCDF(12) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(WQV,1), size(WQV,2), size(WQV,3), WQV, &
            size(WQV_Global,1), size(WQV_Global,2), size(WQV_Global,3), WQV_Global)
    end if

    !---------------------------------------------------------
    ! *** NC_WRITE_CONS(NC_ID, NTI, 12, SAL_Global, KC, 1)
    IF( ISNCDF(1) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(SAL,1), size(SAL,2), SAL, &
            size(SAL_Global,1), size(SAL_Global,2), SAL_Global)
    end if

    !---------------------------------------------------------
    ! *** NC_WRITE_CONS(NC_ID, NTI, 13, TEM_Global, KC, 1)
    IF( ISNCDF(2) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(TEM,1), size(TEM,2), TEM, &
            size(TEM_Global,1), size(TEM_Global,2), TEM_Global)
    end if

    !---------------------------------------------------------
    ! *** NC_WRITE_CONM(NC_ID, NTI, 14, DYE_Global, NDYE, KC, 1)
    IF( ISNCDF(3) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(DYE,1), size(DYE,2), size(DYE,3), DYE, &
            size(DYE_Global,1), size(DYE_Global,2), size(DYE_Global,3), DYE_Global)
    end if

    !---------------------------------------------------------
    ! *** NC_WRITE_CONS(NC_ID, NTI, 15, SFL_Global, KC, 1)
    IF( ISNCDF(4) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(SFL,1), size(SFL,2), SFL, &
            size(SFL_Global,1), size(SFL_Global,2), SFL_Global)
    end if

    ! *** STATUS = NC_WRITE_CONM(NC_ID, NTI, 16, TOX_Global, NTOX, KC, 1)
    IF( ISNCDF(5) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(TOX,1), size(TOX,2), size(TOX,3), TOX,&
            size(TOX_Global,1), size(TOX_Global,2), size(TOX_Global,3), TOX_Global)
    end if

    ! *** STATUS = NC_WRITE_CONM(NC_ID, NTI, 17, SED_Global, NSED, KC, 1)
    IF( ISNCDF(6) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(SED,1), size(SED,2), size(SED,3), SED,&
            size(SED_Global,1), size(SED_Global,2), size(SED_Global,3), SED_Global)
    end if

    ! *** NC_WRITE_CONM(NC_ID, NTI, 18, SND_Global, NSND, KC, 1)
    IF( ISNCDF(7) == 1 )then
        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(SND,1), size(SND,2), size(SND,3), SND,&
            size(SND_Global,1), size(SND_Global,2), size(SND_Global,3), SND_Global)
    end if

    ! *** NC_WRITE_CONS(NC_ID, NTI, 18 + M, WQV(:,:,M), KC, 1)
    IF( ISNCDF(8) == 1 )THEN

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(WQV,1), size(WQV,2), size(WQV,3), WQV, &
            size(WQV_Global,1), size(WQV_Global,2), size(WQV_Global,3), WQV_Global)
    ENDIF

    ! ** LPT
    IF(ISNCDF(9) == 1 )THEN

        !> @todo Oil spill not functional with MPI
        ! IF( ANY(ISOILSPI == 1)) STATUS = NC_WRITE_FIELD(NC_ID, NTI, 41, MOC)

        ! *** NC_WRITE_LPT(NC_ID, NTI)
        Call Gather_Drifter_Arrays

    ENDIF

    ! *** NC_WRITE_SZLJ_BEDTOP(NC_ID, NTI, 42, LAYERACTIVE, KB)
    IF( LSEDZLJ )THEN

        ! *** Reverse Order of Array (TSED)
        DO L=2,LA
            DO K=1,KB
                Reverse_Temp_2D_ActLay(L,K) = LAYERACTIVE(K,L)
            ENDDO
        ENDDO

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_ActLay,1), size(Reverse_Temp_2D_ActLay,2), &
            Reverse_Temp_2D_ActLay, size(Gl_Reverse_Temp_2D_ActLay,1),size(Gl_Reverse_Temp_2D_ActLay,2), Gl_Reverse_Temp_2D_ActLay)

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(TAU,1), TAU, size(TAU_Global,1), TAU_Global)

        j = j + 1
        ! *** Reverse Order of Array (TSED)
        DO L=2,LA
            DO K=1,KB
                Reverse_Temp_2D_TSED(L,K) = TSED(K,L)
            ENDDO
        ENDDO

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_TSED,1), size(Reverse_Temp_2D_TSED,2), Reverse_Temp_2D_TSED, &
            size(Gl_Reverse_Temp_2D_TSED,1),size(Gl_Reverse_Temp_2D_TSED,2), Gl_Reverse_Temp_2D_TSED)


        ! *** Reverse Order of Array (BULKDENS)
        DO L=2,LA
            DO K=1,KB
                Reverse_Temp_2D_BULKDENS(L,K) = BULKDENS(K,L)
            ENDDO
        ENDDO

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(Reverse_Temp_2D_BULKDENS,1), size(Reverse_Temp_2D_BULKDENS,2), Reverse_Temp_2D_BULKDENS, &
            size(Gl_Reverse_Temp_2D_BULKDENS,1),size(Gl_Reverse_Temp_2D_BULKDENS,2), Gl_Reverse_Temp_2D_BULKDENS)

        ! *** Reverse Order of Array (PERSED)
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

        IF(ISNCDF(5) == 1 .AND. ISTRAN(5) > 0 )THEN

            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(KBT,1), KBT, size(KBT_Global,1), KBT_Global)

            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(HBED,1), size(HBED,2), HBED, &
                size(HBED_Global,1), size(HBED_GLobal,2), HBED_Global)
            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(TOXB,1), size(TOXB,2), size(TOXB,3), TOXB, &
                size(TOXB_Global,1), size(TOXB_Global,2), size(TOXB_Global,3), TOXB_Global)

            IF( ICALC_BL > 0 )THEN
                j = j + 1
                Call Assign_Loc_Glob_For_Write(j, size(CBLTOX,1), size(CBLTOX,2), CBLTOX, &
                    size(CBLTOX_Global,1),size(CBLTOX_Global,2), CBLTOX_Global)
            ENDIF
        ENDIF

    ELSEIF( ISBEXP >= 1 .AND. KB > 1 .AND. ( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) )THEN

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(KBT,1), KBT, size(KBT_Global,1), KBT_Global)

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(HBED,1), size(HBED,2), HBED, &
            size(HBED_Global,1), size(HBED_GLobal,2), HBED_Global)

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(BDENBED,1), size(BDENBED,2), BDENBED,  &
            size(BDENBED_Global,1), size(BDENBED_Global,2), BDENBED_Global)

        j = j + 1
        Call Assign_Loc_Glob_For_Write(j, size(PORBED,1), size(PORBED,2), PORBED,  &
            size(PORBED_Global,1), size(PORBED_Global,2), PORBED_Global)

        IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN

            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(SEDB,1), size(SEDB,2), size(SEDB,3), SEDB, &
                size(SEDB_Global,1), size(SEDB_Global,2), size(SEDB_Global,3), SEDB_Global)

            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(SNDB,1), size(SNDB,2), size(SNDB,3), SNDB, &
                size(SNDB_Global,1), size(SNDB_Global,2), size(SNDB_Global,3), SNDB_Global)

        ENDIF

        IF(ISNCDF(5) == 1 .AND. ISTRAN(5) >= 1 )THEN
            j = j + 1
            Call Assign_Loc_Glob_For_Write(j, size(TOXB,1), size(TOXB,2), size(TOXB,3), TOXB, &
                size(TOXB_Global,1), size(TOXB_Global,2), size(TOXB_Global,3), TOXB_Global)
        ENDIF
    ENDIF


    num_arrays_to_write_out = j

    ! *** Call routine that maps, gathers, and sorts to produce the final Global value
    Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)

    ! *** NC_WRITE_SZLJ_BEDTOP(NC_ID, NTI, 42, LAYERACTIVE, KB)
    IF( LSEDZLJ )THEN
        ! *** Reverse global array order
        DO LG=2,LA_Global
            Do k = 1, KB
                LAYERACTIVE_Global(K,LG) = Gl_Reverse_Temp_2D_ActLay(LG,K)
            End do
        End do

        ! *** Reverse global array order
        DO LG=2,LA_Global
            Do k = 1, KB
                TSED_Global(K,LG) = Gl_Reverse_Temp_2D_TSED(LG,K)
            End do
        End do

        ! *** Reverse global array order
        DO LG = 2,LA_Global
            Do k = 1, KB
                BULKDENS_Global(K,LG) = Gl_Reverse_Temp_2D_BULKDENS(LG,K)
            End do
        End do

        ! *** Reverse global array order
        DO LG = 2,LA_Global
            DO nn = 1, NSCM
                Do k = 1, KB
                    PERSED_Global(nn,k,LG) = Gl_Reverse_Temp_3D(LG,K,nn)
                End do
            End do
        End do
    end if

    Deallocate(Reverse_Temp_2D_ActLay)
    Deallocate(Reverse_Temp_2D_TSED)
    Deallocate(Reverse_Temp_2D_BULKDENS)
    Deallocate(Reverse_Temp_3D)
    Deallocate(Gl_Reverse_Temp_2D_ActLay)
    Deallocate(Gl_Reverse_Temp_2D_TSED)
    Deallocate(Gl_Reverse_Temp_2D_BULKDENS)
    Deallocate(Gl_Reverse_Temp_3D)

    end subroutine Map_Write_NetCDF

    end module Mod_Map_Write_NetCDF
