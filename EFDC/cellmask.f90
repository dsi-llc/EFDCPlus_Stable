! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  SUBROUTINE CELLMASK

  ! ***  SUBROUTINE CELLMASK CONVERTS LAND CELLS TO WATER CELLS BY
  ! ***  MASKING VARIABLES.  DEPTHS IN THE MASKED CELLS SHOULD BE INPUT AT
  ! ***  THE END OF THE DXDY.INP FILE.

  ! CHANGE RECORD

  use GLOBAL
  use INFOMOD,only:READSTR

  use Broadcast_Routines
  use Variables_MPI
  use MPI

  implicit none

  integer :: I,J,K,L,LN,LE,M,MMASK,MTYPE
  integer :: I_GL, J_GL
  character(200) :: STR

  if( process_id == master_id )then
    write(*,'(A)')'READING MASK.INP'
    open(1,FILE = 'mask.inp',STATUS = 'UNKNOWN')
    STR = READSTR(1)
    read(1,*,err = 1000) MMASK
  endif

  call Broadcast_Scalar(MMASK, master_id)

  do M = 1,MMASK

    if( process_id == master_id )then
      read(1,*,err = 1000,end = 999) I_GL, J_GL, MTYPE
    endif

    call Broadcast_Scalar(I_GL  ,master_id)
    call Broadcast_Scalar(J_GL 	,master_id)
    call Broadcast_Scalar(MTYPE	,master_id)

    !***Remap to local values
    I = IG2IL(I_GL)
    J = JG2JL(J_GL)

    if( I > 0 .and. I <= IC )then
      if( J > 0 .and. J <= JC )then
        L = LIJ(I,J)
        if( MTYPE == 1 )then
          SUB(L) = 0.
          SUBO(L) = 0.
          UHDYE(L) = 0.
          UHDY1E(L) = 0.
          UHDY2E(L) = 0.
          UMASK(L) = 1
          do K = 1,KC
            U(L,K) = 0.
            U1(L,K) = 0.
            U2(L,K) = 0.
            UHDY(L,K) = 0.
            UHDY1(L,K) = 0.
            UHDY2(L,K) = 0.
            UHDYF(L,K) = 0.
            UHDYF1(L,K) = 0.
            UHDYF2(L,K) = 0.
          enddo
        endif

        if( MTYPE == 2 )then
          SVB(L) = 0.
          SVBO(L) = 0.
          VHDXE(L) = 0.
          VHDX1E(L) = 0.
          VHDX2E(L) = 0.
          VMASK(L)  =  1
          do K = 1,KC
            V(L,K) = 0.
            V1(L,K) = 0.
            V2(L,K) = 0.
            VHDX(L,K) = 0.
            VHDX1(L,K) = 0.
            VHDX2(L,K) = 0.
            VHDXF(L,K) = 0.
            VHDXF1(L,K) = 0.
            VHDXF2(L,K) = 0.
          enddo
        endif

        if( MTYPE == 3 )then   ! *** PMC
          ! *** U Face
          SUB(L) = 0.
          SUBO(L) = 0.
          UHDYE(L) = 0.
          UHDY1E(L) = 0.
          UHDY2E(L) = 0.
          UMASK(L) = 1
          do K = 1,KC
            U(L,K) = 0.
            U1(L,K) = 0.
            U2(L,K) = 0.
            UHDY(L,K) = 0.
            UHDY1(L,K) = 0.
            UHDY2(L,K) = 0.
            UHDYF(L,K) = 0.
            UHDYF1(L,K) = 0.
            UHDYF2(L,K) = 0.
          enddo
          ! *** V Face
          SVB(L) = 0.
          SVBO(L) = 0.
          VHDXE(L) = 0.
          VHDX1E(L) = 0.
          VHDX2E(L) = 0.
          VMASK(L)  =  1
          do K = 1,KC
            V(L,K) = 0.
            V1(L,K) = 0.
            V2(L,K) = 0.
            VHDX(L,K) = 0.
            VHDX1(L,K) = 0.
            VHDX2(L,K) = 0.
            VHDXF(L,K) = 0.
            VHDXF1(L,K) = 0.
            VHDXF2(L,K) = 0.
          enddo
        endif

        if( MTYPE == 4 )then   ! *** PMC - Change to MTYPE 4 for isolated !waters
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
          do K = 1,KC
            B(L,K) = 0.
            B1(L,K) = 0.
            SAL(L,K) = 0.
            SAL1(L,K) = 0.
            TEM(L,K) = TEMO
            TEM1(L,K) = TEMO
            DYE(L,K,1:NDYE) = 0.
            DYE1(L,K,1:NDYE) = 0.
            SED(L,K,1:NSED2) = 0.
            SED1(L,K,1:NSED2) = 0.
            SND(L,K,1:NSND) = 0.
            SND1(L,K,1:NSND) = 0.
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
          enddo
        endif
      endif
    endif
  enddo

999 continue
  close(1)

  ! ***  WRITE READ ERRORS ON CELLMASK
  GOTO 1002

1000 call STOPP('READ ERROR ON FILE MASK.INP ')
1002 continue

  return
  
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

  use GLOBAL
  use Variables_MPI
  use Broadcast_Routines
  use INFOMOD,only:READSTR

  implicit none

  integer :: I, J, L, M, LG, NTMP
  real :: BLANCHORUM, BLDRAFTUOM, BLSILLUM, BLANCHORVM, BLDRAFTVOM, BLSILLVM
  character(200) :: STR

  if( process_id == master_id )then
    write(*,'(A)')'READING LAYERMASK.INP'
    open(1,FILE = 'layermask.inp',STATUS = 'UNKNOWN')
    STR = READSTR(1)
  endif
  
  NTMP = NBLOCKED
  NBLOCKED = 0
  do M = 1,NTMP
    if( process_id == master_id )then
      read(1,*,err = 1000,end = 999) I,J, BLANCHORUM, BLDRAFTUOM, BLSILLUM, BLANCHORVM, BLDRAFTVOM, BLSILLVM
    endif
    
    call Broadcast_Scalar(I, master_id)
    call Broadcast_Scalar(J, master_id)
    call Broadcast_Scalar(BLANCHORUM, master_id)
    call Broadcast_Scalar(BLDRAFTUOM, master_id)
    call Broadcast_Scalar(BLSILLUM,   master_id)
    call Broadcast_Scalar(BLANCHORVM, master_id)
    call Broadcast_Scalar(BLDRAFTVOM, master_id) 
    call Broadcast_Scalar(BLSILLVM,   master_id)
    
    LG = LIJ_Global(I,J)
    L = Map2Local(LG).LL
    
    if( L > 0 )then
      NBLOCKED = NBLOCKED + 1
      LBLOCKED(NBLOCKED) = L

      if( SUBO(L) > 0.0 )then
        BLANCHORU(NBLOCKED) = BLANCHORUM
        BLDRAFTUO(NBLOCKED) = BLDRAFTUOM
        BLSILLU(NBLOCKED)   = BLSILLUM
      endif
      if( SVBO(L) > 0.0 )then
        BLANCHORV(NBLOCKED) = BLANCHORVM
        BLDRAFTVO(NBLOCKED) = BLDRAFTVOM
        BLSILLV(NBLOCKED)   = BLSILLVM
      endif
      
      ! *** SET WAVE FETCH FLAGS
      if( BLDRAFTUO(NBLOCKED) > 0.0 .or. BLSILLU(NBLOCKED) > HP(L) ) UMASK(L) = 1
      if( BLDRAFTVO(NBLOCKED) > 0.0 .or. BLSILLV(NBLOCKED) > HP(L) ) VMASK(L) = 1
    endif
  enddo

  999 continue
  
  close(1)  
  
  return

  ! ***  WRITE READ ERRORS ON CELLMASK
1000 call STOPP('READ ERROR ON FILE LAYERMASK.INP ')

END SUBROUTINE BLOCKING
