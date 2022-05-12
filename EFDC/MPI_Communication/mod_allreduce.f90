! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Performs All_Reduce operations with timing
! @author Paul M. Craig
! @date 2020-05-15
Module MPI_All_Reduce

  ! *** iOP - MPI_MAX = 0x58000001,
  ! *** iOP - MPI_MIN = 0x58000003,
  ! *** iOP - MPI_SUM = 0x58000003,
  
  USE GLOBAL
  USE MPI
  Use Variables_MPI
  
  IMPLICIT NONE
  
  SAVE

  Public :: DSI_All_Reduce

  !***Set up the geneic interface for commuicating ghost cells
  !***This will select the the proper array size depending on what is passed in
  Interface DSI_All_Reduce

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

  Integer, Intent(in)  :: ValIn, iOp, iWait
  Integer, Intent(out) :: ValOut
  Real(8), Intent(out) :: ElapsedTime, WaitTime
  
  Integer(4)          :: IERR
  Real(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS

  IF( num_Processors > 1 )THEN
    IF( iWait > 0 )THEN
      TTDS = DSTIME(0)
      Call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    ELSEIF( iWait < 0 )THEN
      WaitTime = 0.
      Call MPI_barrier(MPI_Comm_World, ierr)
    ENDIF
    TTDS = DSTIME(0)
    
    CALL MPI_ALLREDUCE(ValIn, ValOut, 1, MPI_Integer, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  ELSE
    ValOut = ValIn
    ElapsedTime = 0.
    WaitTime = 0.
  ENDIF
  
  
END SUBROUTINE DSI_All_Reduce_Integer

SUBROUTINE DSI_All_Reduce_Real4(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  Integer, Intent(in)    :: iOp, iWait
  Real(RK4), Intent(in)       :: ValIn
  Real(RK4), Intent(out)      :: ValOut
  Real(RKD), Intent(out) :: ElapsedTime, WaitTime
  
  Integer             :: IERR
  Real(RKD), EXTERNAL :: DSTIME
  Real(RKD)           :: TTDS

  IF( num_Processors > 1 )THEN
    IF( iWait > 0 )THEN
      TTDS = DSTIME(0)
      Call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    ELSEIF( iWait < 0 )THEN
      WaitTime = 0.
      Call MPI_barrier(MPI_Comm_World, ierr)
    ENDIF
    TTDS = DSTIME(0)
  
    CALL MPI_ALLREDUCE(ValIn, ValOut, 1, MPI_Real4, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  ELSE
    ValOut = ValIn
    ElapsedTime = 0.
    WaitTime = 0.
  ENDIF
  
END SUBROUTINE DSI_All_Reduce_Real4

SUBROUTINE DSI_All_Reduce_Real8(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  Integer, Intent(in)    :: iOp, iWait
  Real(8), Intent(in)    :: ValIn
  Real(8), Intent(out)   :: ValOut
  Real(RKD), Intent(out) :: ElapsedTime, WaitTime
  
  Integer(4)          :: IERR
  Real(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS

  IF( num_Processors > 1 )THEN
    IF( iWait > 0 )THEN
      TTDS = DSTIME(0)
      Call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    ELSEIF( iWait < 0 )THEN
      WaitTime = 0.
      Call MPI_barrier(MPI_Comm_World, ierr)
    ENDIF
    TTDS = DSTIME(0)
  
    CALL MPI_ALLREDUCE(ValIn, ValOut, 1, MPI_Real8, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  ELSE
    ValOut = ValIn
    ElapsedTime = 0.
    WaitTime = 0.
  ENDIF
  
END SUBROUTINE DSI_All_Reduce_Real8

SUBROUTINE DSI_All_Reduce_Real48(ValIn, ValOut, iOp, ElapsedTime, iWait, WaitTime)

  Integer,   Intent(in)  :: iOp, iWait
  Real(RK4), Intent(in)  :: ValIn
  Real(RKD), Intent(out) :: ValOut
  Real(RKD), Intent(out) :: ElapsedTime, WaitTime
  Real(RKD)              :: VAL8
  
  Integer(IK4)        :: IERR
  Real(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS

  IF( num_Processors > 1 )THEN
    IF( iWait > 0 )THEN
      TTDS = DSTIME(0)
      Call MPI_barrier(MPI_Comm_World, ierr)
      WaitTime = DSTIME(0) - TTDS
    ELSEIF( iWait < 0 )THEN
      WaitTime = 0.
      Call MPI_barrier(MPI_Comm_World, ierr)
    ENDIF
    TTDS = DSTIME(0)
  
    VAL8 = DBLE(ValIn)
    CALL MPI_ALLREDUCE(VAL8, ValOut, 1, MPI_Real8, iOp, comm_2d, IERR)

    ElapsedTime = DSTIME(0) - TTDS
  ELSE
    ValOut = DBLE(ValIn)
    ElapsedTime = 0.
    WaitTime = 0.
  ENDIF
  
END SUBROUTINE DSI_All_Reduce_Real48

subroutine DSI_All_Reduce2(VarIn1, VarIn2, VarOut1, VarOut2, iOp, ElapsedTime, iWait, WaitTime)

  USE Variables_MPI
  USE MPI
  USE Broadcast_Routines

  Implicit none

  ! *** Passed in variables
  Integer,Intent(in)    :: iOp, iWait
  Real(RK4),Intent(in)  :: VarIn1
  Integer,Intent(in)    :: VarIn2
  Real(RK4),Intent(out) :: VarOut1
  Integer,Intent(out)   :: VarOut2
  Real(RKD),Intent(out) :: ElapsedTime, WaitTime
  Real(RKD),EXTERNAL    :: DSTIME
  REAL(RKD)             :: TTDS
  
  !***local variables
  Integer :: i, j, k, II, L, iMin, iMax
  Integer :: ierr
  Integer :: displ(num_Processors), IRECV(num_Processors)
    
  Real(RK4) :: VarOp
  Real(RK4) :: VarSend(2), VarRec(num_Processors*2)

  IF( num_Processors > 1 )THEN
    IF( iWait > 0 )THEN
      TTDS = DSTIME(0)
      Call MPI_barrier(comm_2d, ierr)
      WaitTime = DSTIME(0) - TTDS
    ELSEIF( iWait < 0 )THEN
      WaitTime = 0.
      Call MPI_barrier(comm_2d, ierr)
    ENDIF
    TTDS = DSTIME(0)
  
    DO i = 1,num_Processors
      displ(i) = (i-1)*2
    ENDDO
    VarRec(:) = 0.
    ierr = 0
    
    VarSend(1) = VarIn1
    VarSend(2) = REAL(VarIn2)
    IRECV = 2
    Call MPI_GatherV(VarSend, 2, MPI_Real4, VarRec, IRECV, displ, MPI_Real4, master_id, comm_2d, ierr)
  
    VarOut2 = 0
    IF( process_id == master_id )THEN
      ! *** Perform the operation
      IF( iOp == MPI_SUM )THEN
        ! *** Sum
        VarOp = 0.0
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp + VarRec(II)
        ENDDO
      
      ELSEIF( iOp == MPI_MAX )THEN
        ! *** Maximum
        iMax = 0
        VarOp = -1.e32
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          IF( VarOp < VarRec(II) )THEN
            VarOp = VarRec(II)
            iMax = i
            VarOut2 = VarRec(II+1)
          ENDIF
        ENDDO

      ELSEIF( iOp == MPI_MIN )THEN
        ! *** Minimum
        iMin = 0
        VarOp = 1.e32
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          IF( VarOp > VarRec(II) )THEN
            VarOp = VarRec(II)
            iMin = i
            VarOut2 = VarRec(II+1)
          ENDIF
        ENDDO
      
      ELSEIF( iOp == MPI_PROD )THEN
        ! *** Product
        VarOp = 1._8
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp * VarRec(II)
        ENDDO

      ENDIF
    ENDIF
    
    Call Broadcast_Scalar(VarOp, master_id )
  
    VarOut1 = VarOp

    ElapsedTime = DSTIME(0) - TTDS
  ELSE
    VarOut1 = VarIn1
    VarOut2 = VarIn2
    ElapsedTime = 0.
    WaitTime = 0.
  ENDIF

END SUBROUTINE DSI_All_Reduce2

subroutine DSI_All_Reduce2_8(VarIn1, VarIn2, VarOut1, VarOut2, iOp, ElapsedTime, iWait, WaitTime)

  USE Variables_MPI
  USE MPI
  USE Broadcast_Routines

  Implicit none

  ! *** Passed in variables
  Integer,Intent(in)    :: iOp, iWait
  Real(RKD),Intent(in)  :: VarIn1
  Integer,Intent(in)    :: VarIn2
  Real(RKD),Intent(out) :: VarOut1
  Integer,Intent(out)   :: VarOut2
  Real(RKD),Intent(out) :: ElapsedTime, WaitTime
  Real(RKD),EXTERNAL    :: DSTIME
  REAL(RKD)             :: TTDS
  
  !***local variables
  Integer :: i, j, k, II, L, iMin, iMax
  Integer :: ierr
  Integer :: displ(num_Processors), IRECV(num_Processors)
    
  Real(RKD) :: VarOp
  Real(RKD) :: VarSend(2), VarRec(num_Processors*2)

  IF( num_Processors > 1 )THEN
    IF( iWait > 0 )THEN
      TTDS = DSTIME(0)
      Call MPI_barrier(comm_2d, ierr)
      WaitTime = DSTIME(0) - TTDS
    ELSEIF( iWait < 0 )THEN
      WaitTime = 0.
      Call MPI_barrier(comm_2d, ierr)
    ENDIF
    TTDS = DSTIME(0)
  
    DO i = 1,num_Processors
      displ(i) = (i-1)*2
    ENDDO
    VarRec(:) = 0.
    ierr = 0
    
    VarSend(1) = VarIn1
    VarSend(2) = REAL(VarIn2,8)
    IRECV = 2
    Call MPI_GatherV(VarSend, 2, MPI_Real8, VarRec, IRECV, displ, MPI_Real8, master_id, comm_2d, ierr)
  
    VarOut2 = 0
    IF( process_id == master_id )THEN
      ! *** Perform the operation
      IF( iOp == MPI_SUM )THEN
        ! *** Sum
        VarOp = 0.0
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp + VarRec(II)
        ENDDO
      
      ELSEIF( iOp == MPI_MAX )THEN
        ! *** Maximum
        iMax = 0
        VarOp = -1.e32
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          IF( VarOp < VarRec(II) )THEN
            VarOp = VarRec(II)
            iMax = i
            VarOut2 = VarRec(II+1)
          ENDIF
        ENDDO

      ELSEIF( iOp == MPI_MIN )THEN
        ! *** Minimum
        iMin = 0
        VarOp = 1.e32
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          IF( VarOp > VarRec(II) )THEN
            VarOp = VarRec(II)
            iMin = i
            VarOut2 = VarRec(II+1)
          ENDIF
        ENDDO
      
      ELSEIF( iOp == MPI_PROD )THEN
        ! *** Product
        VarOp = 1._8
        DO i=1,num_Processors
          II = (i-1)*2 + 1
          VarOp = VarOp * VarRec(II)
        ENDDO

      ENDIF
    ENDIF
    
    Call Broadcast_Scalar(VarOp, master_id )
  
    VarOut1 = VarOp

    ElapsedTime = DSTIME(0) - TTDS
  ELSE
    VarOut1 = VarIn1
    VarOut2 = VarIn2
    ElapsedTime = 0.
    WaitTime = 0.
  ENDIF

END SUBROUTINE DSI_All_Reduce2_8

END MODULE MPI_All_Reduce
