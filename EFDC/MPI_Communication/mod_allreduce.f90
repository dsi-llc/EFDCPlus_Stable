! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Performs All_Reduce operations with timing
! @author Paul M. Craig
! @date 2020-05-15
Module MPI_All_Reduce

  ! *** iOP - MPI_MAX = 0x58000001,
  ! *** iOP - MPI_MIN = 0x58000003,
  ! *** iOP - MPI_SUM = 0x58000003,
  
  use GLOBAL
  use MPI
  use Variables_MPI
  
  implicit none
  
  save

  Public :: DSI_All_Reduce

  !***Set up the geneic interface for commuicating ghost cells
  !***This will select the the proper array size depending on what is passed in
  interface DSI_All_Reduce

  Module Procedure DSI_All_Reduce_Integer,   &
                   DSI_All_Reduce_Real4,     &
                   DSI_All_Reduce_Real8,     &
                   DSI_All_Reduce_Real48,    &
                   DSI_All_Reduce2_8,        &
                   DSI_All_Reduce2
  End interface

  ! *** Technical Note:  The STATUS variable returned by MPI_SEND and MPI_RECV is the number of BYTES sent or received.
  
  ! *** Next are each of the subroutines for the procedure
  Contains
 
SUBROUTINE DSI_All_Reduce_Integer(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  integer, Intent(in)  :: ValIn, iOp, iWait
  integer, Intent(out) :: ValOut
  real(8), Intent(out) :: ElapsedTime, WaitTime
  
  Integer(4)          :: IERR
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS

  if( num_Processors > 1 )then
    if( iWait > 0 )then
      TTDS = DSTIME(0)
      call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    elseif( iWait < 0 )then
      WaitTime = 0.
      call MPI_barrier(MPI_Comm_World, ierr)
    endif
    TTDS = DSTIME(0)
    
    call MPI_ALLREDUCE(ValIn, ValOut, 1, MPI_Integer, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  else
    ValOut = ValIn
    ElapsedTime = 0.
    WaitTime = 0.
  endif
  
  
END SUBROUTINE DSI_All_Reduce_Integer

SUBROUTINE DSI_All_Reduce_Real4(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  integer, Intent(in)    :: iOp, iWait
  real(RK4), Intent(in)       :: ValIn
  real(RK4), Intent(out)      :: ValOut
  real(RKD), Intent(out) :: ElapsedTime, WaitTime
  
  integer             :: IERR
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS

  if( num_Processors > 1 )then
    if( iWait > 0 )then
      TTDS = DSTIME(0)
      call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    elseif( iWait < 0 )then
      WaitTime = 0.
      call MPI_barrier(MPI_Comm_World, ierr)
    endif
    TTDS = DSTIME(0)
  
    call MPI_ALLREDUCE(ValIn, ValOut, 1, MPI_Real4, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  else
    ValOut = ValIn
    ElapsedTime = 0.
    WaitTime = 0.
  endif
  
END SUBROUTINE DSI_All_Reduce_Real4

SUBROUTINE DSI_All_Reduce_Real8(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  integer, Intent(in)    :: iOp, iWait
  real(8), Intent(in)    :: ValIn
  real(8), Intent(out)   :: ValOut
  real(RKD), Intent(out) :: ElapsedTime, WaitTime
  
  Integer(4)          :: IERR
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS

  if( num_Processors > 1 )then
    if( iWait > 0 )then
      TTDS = DSTIME(0)
      call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    elseif( iWait < 0 )then
      WaitTime = 0.
      call MPI_barrier(MPI_Comm_World, ierr)
    endif
    TTDS = DSTIME(0)
  
    call MPI_ALLREDUCE(ValIn, ValOut, 1, MPI_Real8, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  else
    ValOut = ValIn
    ElapsedTime = 0.
    WaitTime = 0.
  endif
  
END SUBROUTINE DSI_All_Reduce_Real8

SUBROUTINE DSI_All_Reduce_Real48(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  integer,   Intent(in)  :: iOp, iWait
  real(RK4), Intent(in)  :: ValIn
  real(RKD), Intent(out) :: ValOut
  real(RKD), Intent(out) :: ElapsedTime, WaitTime
  real(RKD)              :: VAL8
  
  Integer(IK4)        :: IERR
  real(RKD), external :: DSTIME
  real(RKD)           :: TTDS

  if( num_Processors > 1 )then
    if( iWait > 0 )then
      TTDS = DSTIME(0)
      call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    elseif( iWait < 0 )then
      WaitTime = 0.
      call MPI_barrier(MPI_Comm_World, ierr)
    endif
    TTDS = DSTIME(0)
  
    VAL8 = DBLE(ValIn)
    call MPI_ALLREDUCE(VAL8, ValOut, 1, MPI_Real8, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  else
    ValOut = DBLE(ValIn)
    ElapsedTime = 0.
    WaitTime = 0.
  endif
  
END SUBROUTINE DSI_All_Reduce_Real48

subroutine DSI_All_Reduce2(VarIn1, VarIn2, VarOut1, VarOut2, iOp, ElapsedTime, iWait, WaitTime)

  use Variables_MPI
  use MPI
  use Broadcast_Routines

  implicit none

  ! *** Passed in variables
  integer,Intent(in)    :: iOp, iWait
  real(RK4),Intent(in)  :: VarIn1
  integer,Intent(in)    :: VarIn2
  real(RK4),Intent(out) :: VarOut1
  integer,Intent(out)   :: VarOut2
  real(RKD),Intent(out) :: ElapsedTime, WaitTime
  real(RKD),external    :: DSTIME
  real(RKD)             :: TTDS
  
  !***local variables
  integer :: i, j, k, II, L, iMin, iMax
  integer :: ierr
  integer :: displ(num_Processors), IRECV(num_Processors)
    
  real(RK4) :: VarOp
  real(RK4) :: VarSend(2), VarRec(num_Processors*2)

  if( num_Processors > 1 )then
    if( iWait > 0 )then
      TTDS = DSTIME(0)
      call MPI_barrier(comm_2d, ierr)
      WaitTime = DSTIME(0) - TTDS
    elseif( iWait < 0 )then
      WaitTime = 0.
      call MPI_barrier(comm_2d, ierr)
    endif
    TTDS = DSTIME(0)
  
    do i = 1,num_Processors
      displ(i) = (i-1)*2
    enddo
    VarRec(:) = 0.
    ierr = 0
    
    VarSend(1) = VarIn1
    VarSend(2) = REAL(VarIn2)
    IRECV = 2
    call MPI_GatherV(VarSend, 2, MPI_Real4, VarRec, IRECV, displ, MPI_Real4, master_id, comm_2d, ierr)
  
    VarOut2 = 0
    if( process_id == master_id )then
      ! *** Perform the operation
      if( iOp == MPI_SUM )then
        ! *** Sum
        VarOp = 0.0
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp + VarRec(II)
        enddo
      
      elseif( iOp == MPI_MAX )then
        ! *** Maximum
        iMax = 0
        VarOp = -1.e32
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          if( VarOp < VarRec(II) )then
            VarOp = VarRec(II)
            iMax = i
            VarOut2 = VarRec(II+1)
          endif
        enddo

      elseif( iOp == MPI_MIN )then
        ! *** Minimum
        iMin = 0
        VarOp = 1.e32
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          if( VarOp > VarRec(II) )then
            VarOp = VarRec(II)
            iMin = i
            VarOut2 = VarRec(II+1)
          endif
        enddo
      
      elseif( iOp == MPI_PROD )then
        ! *** Product
        VarOp = 1._8
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp * VarRec(II)
        enddo

      endif
    endif
    
    call Broadcast_Scalar(VarOp, master_id )
  
    VarOut1 = VarOp

    ElapsedTime = DSTIME(0) - TTDS
  else
    VarOut1 = VarIn1
    VarOut2 = VarIn2
    ElapsedTime = 0.
    WaitTime = 0.
  endif

END SUBROUTINE DSI_All_Reduce2

subroutine DSI_All_Reduce2_8(VarIn1, VarIn2, VarOut1, VarOut2, iOp, ElapsedTime, iWait, WaitTime)

  use Variables_MPI
  use MPI
  use Broadcast_Routines

  implicit none

  ! *** Passed in variables
  integer,Intent(in)    :: iOp, iWait
  real(RKD),Intent(in)  :: VarIn1
  integer,Intent(in)    :: VarIn2
  real(RKD),Intent(out) :: VarOut1
  integer,Intent(out)   :: VarOut2
  real(RKD),Intent(out) :: ElapsedTime, WaitTime
  real(RKD),external    :: DSTIME
  real(RKD)             :: TTDS
  
  !***local variables
  integer :: i, j, k, II, L, iMin, iMax
  integer :: ierr
  integer :: displ(num_Processors), IRECV(num_Processors)
    
  real(RKD) :: VarOp
  real(RKD) :: VarSend(2), VarRec(num_Processors*2)

  if( num_Processors > 1 )then
    if( iWait > 0 )then
      TTDS = DSTIME(0)
      call MPI_barrier(comm_2d, ierr)
      WaitTime = DSTIME(0) - TTDS
    elseif( iWait < 0 )then
      WaitTime = 0.
      call MPI_barrier(comm_2d, ierr)
    endif
    TTDS = DSTIME(0)
  
    do i = 1,num_Processors
      displ(i) = (i-1)*2
    enddo
    VarRec(:) = 0.
    ierr = 0
    
    VarSend(1) = VarIn1
    VarSend(2) = REAL(VarIn2,8)
    IRECV = 2
    call MPI_GatherV(VarSend, 2, MPI_Real8, VarRec, IRECV, displ, MPI_Real8, master_id, comm_2d, ierr)
  
    VarOut2 = 0
    if( process_id == master_id )then
      ! *** Perform the operation
      if( iOp == MPI_SUM )then
        ! *** Sum
        VarOp = 0.0
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp + VarRec(II)
        enddo
      
      elseif( iOp == MPI_MAX )then
        ! *** Maximum
        iMax = 0
        VarOp = -1.e32
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          if( VarOp < VarRec(II) )then
            VarOp = VarRec(II)
            iMax = i
            VarOut2 = VarRec(II+1)
          endif
        enddo

      elseif( iOp == MPI_MIN )then
        ! *** Minimum
        iMin = 0
        VarOp = 1.e32
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          if( VarOp > VarRec(II) )then
            VarOp = VarRec(II)
            iMin = i
            VarOut2 = VarRec(II+1)
          endif
        enddo
      
      elseif( iOp == MPI_PROD )then
        ! *** Product
        VarOp = 1._8
        do i = 1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp * VarRec(II)
        enddo

      endif
    endif
    
    call Broadcast_Scalar(VarOp, master_id )
  
    VarOut1 = VarOp

    ElapsedTime = DSTIME(0) - TTDS
  else
    VarOut1 = VarIn1
    VarOut2 = VarIn2
    ElapsedTime = 0.
    WaitTime = 0.
  endif

END SUBROUTINE DSI_All_Reduce2_8

END MODULE MPI_All_Reduce
