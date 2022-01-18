! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CELLMASK

  ! **  SUBROUTINE CELLMASK CONVERTS LAND CELLS TO WATER CELLS BY
  ! **  MASKING VARIABLES.  DEPTHS IN THE MASKED CELLS SHOULD BE INPUT AT
  ! **  THE END OF THE DXDY.INP FILE.

  ! CHANGE RECORD

  USE GLOBAL
  USE INFOMOD,ONLY:READSTR

#ifdef _MPI
  Use Broadcast_Routines
  Use Variables_MPI
  Use MPI
#endif

  IMPLICIT NONE

  INTEGER :: I,J,K,L,LN,LE,M,MMASK,MTYPE
  INTEGER :: I_GL, J_GL
  CHARACTER(200) :: STR

  if( process_id == master_id )THEN
    WRITE(*,'(A)')'READING MASK.INP'
    OPEN(1,FILE='mask.inp',STATUS='UNKNOWN')
    STR = READSTR(1)
    READ(1,*,ERR=1000) MMASK
  end if

  Call Broadcast_Scalar(MMASK, master_id)

  DO M=1,MMASK

    if( process_id == master_id )THEN
      READ(1,*,ERR=1000,END=999) I_GL, J_GL, MTYPE
    end if

    Call Broadcast_Scalar(I_GL  ,master_id)
    Call Broadcast_Scalar(J_GL 	,master_id)
    Call Broadcast_Scalar(MTYPE	,master_id)

    !***Remap to local values
    I = IG2IL(I_GL)
    J = JG2JL(J_GL)

    IF( I > 0 .AND. I <= IC )THEN
      IF( J > 0 .AND. J <= JC )THEN
        L = LIJ(I,J)
        IF( MTYPE == 1 )THEN
          SUB(L) = 0.
          SUBO(L) = 0.
          UHDYE(L) = 0.
          UHDY1E(L) = 0.
          UHDY2E(L) = 0.
          UMASK(L) = 1
          DO K = 1,KC
            U(L,K) = 0.
            U1(L,K) = 0.
            U2(L,K) = 0.
            UHDY(L,K) = 0.
            UHDY1(L,K) = 0.
            UHDY2(L,K) = 0.
            UHDYF(L,K) = 0.
            UHDYF1(L,K) = 0.
            UHDYF2(L,K) = 0.
          ENDDO
        ENDIF

        IF( MTYPE == 2 )THEN
          SVB(L) = 0.
          SVBO(L) = 0.
          VHDXE(L) = 0.
          VHDX1E(L) = 0.
          VHDX2E(L) = 0.
          VMASK(L)  =  1
          DO K = 1,KC
            V(L,K) = 0.
            V1(L,K) = 0.
            V2(L,K) = 0.
            VHDX(L,K) = 0.
            VHDX1(L,K) = 0.
            VHDX2(L,K) = 0.
            VHDXF(L,K) = 0.
            VHDXF1(L,K) = 0.
            VHDXF2(L,K) = 0.
          ENDDO
        ENDIF

        IF( MTYPE == 3 )THEN   ! *** PMC
          ! *** U Face
          SUB(L) = 0.
          SUBO(L) = 0.
          UHDYE(L) = 0.
          UHDY1E(L) = 0.
          UHDY2E(L) = 0.
          UMASK(L) = 1
          DO K = 1,KC
            U(L,K) = 0.
            U1(L,K) = 0.
            U2(L,K) = 0.
            UHDY(L,K) = 0.
            UHDY1(L,K) = 0.
            UHDY2(L,K) = 0.
            UHDYF(L,K) = 0.
            UHDYF1(L,K) = 0.
            UHDYF2(L,K) = 0.
          ENDDO
          ! *** V Face
          SVB(L) = 0.
          SVBO(L) = 0.
          VHDXE(L) = 0.
          VHDX1E(L) = 0.
          VHDX2E(L) = 0.
          VMASK(L)  =  1
          DO K = 1,KC
            V(L,K) = 0.
            V1(L,K) = 0.
            V2(L,K) = 0.
            VHDX(L,K) = 0.
            VHDX1(L,K) = 0.
            VHDX2(L,K) = 0.
            VHDXF(L,K) = 0.
            VHDXF1(L,K) = 0.
            VHDXF2(L,K) = 0.
          ENDDO
        ENDIF

        IF( MTYPE == 4 )THEN   ! *** PMC - Change to MTYPE 4 for isolated !waters
          LN = LNC(L)
          LE = LEC(L)
          SUB(L) = 0.
          SUBO(L) = 0.
          UHDYE(L) = 0.
          UHDY1E(L) = 0.
          UHDY2E(L) = 0.
          SUB(LE) = 0.
          SUBO(LE) = 0.
          UHDYE(LE) = 0.
          UHDY1E(LE) = 0.
          UHDY2E(LE) = 0.
          SVB(L) = 0.
          SVBO(L) = 0.
          VHDXE(L) = 0.
          VHDX1E(L) = 0.
          VHDX2E(L) = 0.
          SVB(LN) = 0.
          SVBO(LN) = 0.
          VHDXE(LN) = 0.
          VHDX1E(LN) = 0.
          VHDX2E(LN) = 0.
          P(L) = 0.
          P1(L) = 0.
          DO K = 1,KC
            B(L,K) = 0.
            B1(L,K) = 0.
            SAL(L,K) = 0.
            SAL1(L,K) = 0.
            TEM(L,K) = TEMO
            TEM1(L,K) = TEMO
            DYE(L,K,1:NDYE) = 0.
            DYE1(L,K,1:NDYE) = 0.
            SED(L,K,1) = 0.
            SED1(L,K,1) = 0.
            QQ(L,K) = 0.
            QQ1(L,K) = 0.
            QQL(L,K) = 0.
            QQL1(L,K) = 0.
            U(L,K) = 0.
            U1(L,K) = 0.
            U2(L,K) = 0.
            UHDY(L,K) = 0.
            UHDY1(L,K) = 0.
            UHDY2(L,K) = 0.
            UHDYF(L,K) = 0.
            UHDYF1(L,K) = 0.
            UHDYF2(L,K) = 0.
            U(LE,K) = 0.
            U1(LE,K) = 0.
            U2(LE,K) = 0.
            UHDY(LE,K) = 0.
            UHDY1(LE,K) = 0.
            UHDY2(LE,K) = 0.
            V(L,K) = 0.
            V1(L,K) = 0.
            V2(L,K) = 0.
            VHDX(L,K) = 0.
            VHDX1(L,K) = 0.
            VHDX2(L,K) = 0.
            VHDXF(L,K) = 0.
            VHDXF1(L,K) = 0.
            VHDXF2(L,K) = 0.
            V(LN,K) = 0.
            V1(LN,K) = 0.
            V2(LN,K) = 0.
            VHDX(LN,K) = 0.
            VHDX1(LN,K) = 0.
            VHDX2(LN,K) = 0.
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDDO

999 CONTINUE
  CLOSE(1)

  ! **  WRITE READ ERRORS ON CELLMASK
  GOTO 1002

1000 CALL STOPP('READ ERROR ON FILE MASK.INP ')
1002 CONTINUE

  RETURN
  
END SUBROUTINE CELLMASK

SUBROUTINE BLOCKING

  ! *** ************************************************************************
  !
  ! ***  SUBROUTINE BLOCKING READS THE DRAFT DEPTHS AND/OR SILL HEIGHTS TO
  ! ***  ALLOW LAYER BY LAYER FACE BLOCKING
  !
  !---------------------------------------------------------------------------
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !---------------------------------------------------------------------------
  !    2019-08       PAUL M. CRAIG     IMPLEMENTED LAYER BY LAYER FACE BLOCKNG

  USE GLOBAL
  Use Variables_MPI
  Use Broadcast_Routines
  USE INFOMOD,ONLY:READSTR

  IMPLICIT NONE

  INTEGER :: I, J, L, M, LG, NTMP
  REAL :: BLANCHORUM, BLDRAFTUOM, BLSILLUM, BLANCHORVM, BLDRAFTVOM, BLSILLVM
  CHARACTER(200) :: STR

  if( process_id == master_id )then
    WRITE(*,'(A)')'READING LAYERMASK.INP'
    OPEN(1,FILE='layermask.inp',STATUS='UNKNOWN')
    STR = READSTR(1)
  endif
  
  NTMP = NBLOCKED
  NBLOCKED = 0
  DO M=1,NTMP
    if( process_id == master_id )then
        READ(1,*,ERR=1000,END=999) I,J, BLANCHORUM, BLDRAFTUOM, BLSILLUM, BLANCHORVM, BLDRAFTVOM, BLSILLVM
    endif
    
    Call Broadcast_Scalar(I, master_id)
    Call Broadcast_Scalar(J, master_id)
    Call Broadcast_Scalar(BLANCHORUM, master_id)
    Call Broadcast_Scalar(BLDRAFTUOM, master_id)
    Call Broadcast_Scalar(BLSILLUM,   master_id)
    Call Broadcast_Scalar(BLANCHORVM, master_id)
    Call Broadcast_Scalar(BLDRAFTVOM, master_id) 
    Call Broadcast_Scalar(BLSILLVM,   master_id)
    
    LG = LIJ_Global(I,J)
    L = Map2Local(LG).LL
    
    IF( L > 0 )THEN
      NBLOCKED = NBLOCKED + 1
      LBLOCKED(NBLOCKED) = L

      BLANCHORU(NBLOCKED) = BLANCHORUM
      BLDRAFTUO(NBLOCKED) = BLDRAFTUOM
      BLSILLU(NBLOCKED)   = BLSILLUM
      BLANCHORV(NBLOCKED) = BLANCHORVM
      BLDRAFTVO(NBLOCKED) = BLDRAFTVOM
      BLSILLV(NBLOCKED)   = BLSILLVM
      
      ! *** SET WAVE FETCH FLAGS
      IF( BLDRAFTUO(NBLOCKED) > 0.0 .OR. BLSILLU(NBLOCKED) > HP(L) ) UMASK(L) = 1
      IF( BLDRAFTVO(NBLOCKED) > 0.0 .OR. BLSILLV(NBLOCKED) > HP(L) ) VMASK(L) = 1
    ENDIF
  ENDDO

  999 CONTINUE
  
  if( process_id == master_id )then
    CLOSE(1)  
  endif
  
  !Call Broadcast_Scalar(NBLOCKED, master_id)
  !Call Broadcast_Array(BLANCHORU, master_id)
  !Call Broadcast_Array(BLDRAFTUO, master_id)
  !Call Broadcast_Array(BLSILLU,   master_id)
  !Call Broadcast_Array(BLANCHORV, master_id)
  !Call Broadcast_Array(BLDRAFTVO, master_id)
  !Call Broadcast_Array(BLSILLV,   master_id)
  !Call Broadcast_Array(UMASK,     master_id)
  !Call Broadcast_Array(VMASK,     master_id)
  
  RETURN

  ! **  WRITE READ ERRORS ON CELLMASK
1000 CALL STOPP('READ ERROR ON FILE LAYERMASK.INP ')

END SUBROUTINE BLOCKING
