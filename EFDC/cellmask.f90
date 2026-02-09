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
  use Allocate_Initialize

  use MPI
  use Variables_MPI
  use Broadcast_Routines
  use Variables_MPI_Write_Out

  use INFOMOD,only:READSTR

  implicit none

  integer :: I, J, K, L, LN, LE, M, MMASK, MTYPE
  integer :: I_GL, J_GL, LG
  character(200) :: STR

  integer,allocatable,dimension(:) :: MASKLIST

  call AllocateDSI( MASKLIST, LCM_Global, 0 )

  if( process_id == master_id )then
    write(*,'(A)')'READING MASK.INP'
    open(1,FILE = 'mask.inp',STATUS = 'UNKNOWN')
    STR = READSTR(1)
    read(1,*,err = 1000) MMASK

    do M = 1,MMASK
      read(1,*,err = 1000,end = 999) I_GL, J_GL, MTYPE

      LG = LIJ_Global(I_GL, J_GL)
      MASKLIST(LG) = MTYPE

      ! *** Assign wave fetch flags based on masks
      if( MTYPE == 1 )then
        UMASK_Global(LG) = 1

      elseif( MTYPE == 2 )then
        VMASK_Global(LG) = 1

      elseif( MTYPE == 3 )then
        UMASK_Global(LG) = 1
        VMASK_Global(LG) = 1

      elseif( MTYPE == 4 )then
        UMASK_Global(LG) = 1
        VMASK_Global(LG) = 1

        UMASK_Global(LEC_Global(LG)) = 1
        VMASK_Global(LNC_Global(LG)) = 1
      endif
    enddo
  endif

  call Broadcast_Scalar(MMASK,   master_id)
  call Broadcast_Array(MASKLIST, master_id)
  call Broadcast_Array(UMASK_Global, master_id)
  call Broadcast_Array(VMASK_Global, master_id)

  ! *** Map to local values
  do LG = 2,LA_GLOBAL
    L = Map2Local(LG).LL
    if( L > 1 .and. MASKLIST(LG) > 0 )then
      MTYPE = MASKLIST(LG)

      if( MTYPE == 1 )then
        ! *** U face only
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
        ! *** V face only
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

      if( MTYPE == 3 )then
        ! *** U and V faces

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

      if( MTYPE == 4 )then   ! *** Change to MTYPE 4 for isolated !waters
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
  use Allocate_Initialize

  use Variables_MPI
  use Broadcast_Routines
  use Variables_MPI_Write_Out

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

  call AllocateDSI( I1D_Global, NBLOCKED, 0 )
  call AllocateDSI( R2D_Global, NBLOCKED, 6, 0.0 )

  do M = 1,NBLOCKED
    if( process_id == master_id )then
      ! ***                                BLANCHORUM, BLDRAFTUOM, BLSILLUM, BLANCHORVM, BLDRAFTVOM, BLSILLVM
      read(1,*,err = 1000,end = 999) I, J, R2D_Global(M,1:6)

      LG = LIJ_Global(I,J)
      I1D_Global(M) = LG

      ! *** Set wave fetch flags
      if( R2D_Global(M,2) > 0.0 .or. R2D_Global(M,3) > HP_Global(LG) )then
        UMASK_Global(LG) = 1
      endif
      if( R2D_Global(M,5) > 0.0 .or. R2D_Global(M,6) > HP_Global(LG) )then
        VMASK_Global(LG) = 1
      endif

    endif
  enddo

  call Broadcast_Array(I1D_Global, master_id)
  call Broadcast_Array(R2D_Global, master_id)

  ! *** Map masks to local
  NTMP = NBLOCKED
  NBLOCKED = 0
  do M = 1,NTMP
    LG = I1D_Global(M)
    L = Map2Local(LG).LL

    if( L > 0 )then
      NBLOCKED = NBLOCKED + 1
      LBLOCKED(NBLOCKED) = L

      if( SUBO(L) > 0.0 )then
        BLANCHORU(NBLOCKED) = R2D_Global(M,1)
        BLDRAFTUO(NBLOCKED) = R2D_Global(M,2)
        BLSILLU(NBLOCKED)   = R2D_Global(M,3)
      endif
      if( SVBO(L) > 0.0 )then
        BLANCHORV(NBLOCKED) = R2D_Global(M,4)
        BLDRAFTVO(NBLOCKED) = R2D_Global(M,5)
        BLSILLV(NBLOCKED)   = R2D_Global(M,6)
      endif

      ! *** Set path intersect flags
      if( BLDRAFTUO(NBLOCKED) > 0.0 .or. BLSILLU(NBLOCKED) > HP(L) )then
        UMASK(L) = 1
      endif
      if( BLDRAFTVO(NBLOCKED) > 0.0 .or. BLSILLV(NBLOCKED) > HP(L) )then
        VMASK(L) = 1
      endif
    endif
  enddo

999 continue

  close(1)

  deallocate(I1D_Global)   
  deallocate(R2D_Global)   
     
  return

  ! ***  WRITE READ ERRORS ON CELLMASK
1000 call STOPP('READ ERROR ON FILE LAYERMASK.INP ')

  END SUBROUTINE BLOCKING
