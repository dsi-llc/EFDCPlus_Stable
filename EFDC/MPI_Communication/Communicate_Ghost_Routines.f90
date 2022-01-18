! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! Generic interface for the communication of ghost cells

Module Communicate_Ghost_Routines

  !Use Mod_Communicate_Functions

  IMPLICIT NONE

  SAVE

  Public :: communicate_ghost_cells

  !***Set up the geneic interface for commuicating ghost cells
  !***This will select the the proper array size depending on what is passed in
  Interface Communicate_ghost_cells

  Module Procedure Communicate_1D_Real,    &
                   Communicate_1D_Real8,   &
                   Communicate_1D_Integer, &
                   Communicate_3D_Real,    &
                   Communicate_3D_Real8,   &
                   Communicate_3D_Integer, &
                   Communicate_4D_Real,    &
                   Communicate_4D_Real8,   &
                   Communicate_4D_Integer, &
                   Communicate_4D0_Real,   &
                   Communicate_4D0_Real8
  End interface

  ! *** Technical Note:  The STATUS variable returned by MPI_SEND and MPI_RECV is the number of BYTES sent or received.
  
  !***Next are each of the subroutines for the procedure
  Contains

  

  ! @details Communicates pressure boundary values amongs processors
  ! @author Zander Mausolff
  ! @date 9/4/2019
  ! @parm pardata(LCM)

  ! ********************************************************************************************
  ! *** Sending values should always be from the active cells
  ! *** Receiving values should always be placed into the ghost cells
  ! ********************************************************************************************
  
  Subroutine Communicate_1D_Real(pardata, var_name)

  Use Variables_MPI
  Use MPI

  Implicit None

  !***Read in variables
  Real(4), Intent(inout), DIMENSION(LCM) ::  pardata
  Character(len=*), Intent(in) :: var_name
  !***Local variables

  Integer :: I, J, II, L
  Integer :: length_arg
  INTEGER :: IERR, status_message(MPI_Status_Size)

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVS

  IF(.NOT.ALLOCATED(PSENDW))THEN
    ALLOCATE(PSENDW(max_width_y*2))
    ALLOCATE(PSENDE(max_width_y*2))
    ALLOCATE(PRECVE(max_width_y*2))
    ALLOCATE(PRECVW(max_width_y*2))
    
    ALLOCATE(PSENDN(max_width_x*2))
    ALLOCATE(PSENDS(max_width_x*2))
    ALLOCATE(PRECVN(max_width_x*2))
    ALLOCATE(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  ENDIF
    
  IF( MPI_DEBUG_FLAG == .TRUE. )THEN
    Call WriteBreak(mpi_comm_unit)
    write(mpi_comm_unit, '(a)') var_name
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = pardata(L)
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(PSENDW, length_arg, mpi_real4, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVE, length_arg, mpi_real4, nbr_east, nbr_east, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      pardata(L) = PRECVE(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDE, length_arg, mpi_real4, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVW, length_arg, mpi_real4, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = PRECVW(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDN, length_arg, mpi_real4, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVS, length_arg, mpi_real4, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = PRECVS(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDS, length_arg, mpi_real4, nbr_south, process_id,  comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVN, length_arg, mpi_real4, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = PRECVN(II)
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_1D_Real
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  ! ********************************************************************************************
  ! *** Sending values should always be from the active cells
  ! *** Receiving values should always be placed into the ghost cells
  ! ********************************************************************************************
  
  Subroutine Communicate_1D_Real8(pardata, var_name)

  Use Variables_MPI
  Use MPI

  Implicit None

  !***Read in variables
  Real(8), Intent(inout), DIMENSION(LCM)::  pardata
  Character(len=*), Intent(in) :: var_name
  !***Local variables

  Integer :: I, J, II, L
  Integer :: length_arg
  INTEGER :: IERR, status_message(MPI_Status_Size)

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVS

  IF(.NOT.ALLOCATED(PSENDW))THEN
    ALLOCATE(PSENDW(max_width_y*2))
    ALLOCATE(PSENDE(max_width_y*2))
    ALLOCATE(PRECVE(max_width_y*2))
    ALLOCATE(PRECVW(max_width_y*2))
    
    ALLOCATE(PSENDN(max_width_x*2))
    ALLOCATE(PSENDS(max_width_x*2))
    ALLOCATE(PRECVN(max_width_x*2))
    ALLOCATE(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  ENDIF
    
  IF( MPI_DEBUG_FLAG == .TRUE. )THEN
    Call WriteBreak(mpi_comm_unit)
    write(mpi_comm_unit, '(a)') var_name
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = pardata(L)
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(PSENDW, length_arg, mpi_real8, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVE, length_arg, mpi_real8, nbr_east, nbr_east, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      pardata(L) = PRECVE(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDE, length_arg, mpi_real8, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVW, length_arg, mpi_real8, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = PRECVW(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDN, length_arg, mpi_real8, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVS, length_arg, mpi_real8, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = PRECVS(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDS, length_arg, mpi_real8, nbr_south, process_id,  comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVN, length_arg, mpi_real8, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = PRECVN(II)
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_1D_Real8
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  ! @details Communicates pressure boundary values amongs processors
  ! @author Zander Mausolff
  ! @date 9/4/2019
  ! @parameters pardata(LCM)

  SUBROUTINE Communicate_1D_Integer(pardata, var_name)

  Use Variables_MPI
  Use MPI

  Implicit None

  !***Read in variables
  Integer, Intent(inout), DIMENSION(LCM)::  pardata
  Character(len=*), Intent(in) :: var_name
  !***Local variables

  Integer :: I, J, II, L
  Integer :: length_arg
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)

  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDW
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDE
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVE
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVW
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDN
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PSENDS
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVN
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::PRECVS

  IF(.NOT.ALLOCATED(PSENDW))THEN
    ALLOCATE(PSENDW(max_width_y*2))
    ALLOCATE(PSENDE(max_width_y*2))
    ALLOCATE(PRECVE(max_width_y*2))
    ALLOCATE(PRECVW(max_width_y*2))
    ALLOCATE(PSENDN(max_width_x*2))
    ALLOCATE(PSENDS(max_width_x*2))
    ALLOCATE(PRECVN(max_width_x*2))
    ALLOCATE(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  ENDIF

  IF( MPI_DEBUG_FLAG == .TRUE. )THEN
    Call WriteBreak(mpi_comm_unit)
    write(mpi_comm_unit, '(a)') var_name
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = pardata(L)
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(PSENDW, length_arg, mpi_integer, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVE, length_arg, mpi_integer, nbr_east, nbr_east, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      pardata(L) = PRECVE(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDE, length_arg, mpi_integer, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVW, length_arg, mpi_integer, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = PRECVW(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDN, length_arg, mpi_integer, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVS, length_arg, mpi_integer, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = PRECVS(II)
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = pardata(L)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(PSENDS, length_arg, mpi_integer, nbr_south, process_id,  comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    ENDDO
    length_arg = II
    CALL MPI_RECV(PRECVN, length_arg, mpi_integer, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = PRECVN(II)
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_1D_Integer
  
  !-----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! @details This routine exchanges data for values for arrays with
  !!         the indexing: (L, K) i.e. ==> (cell index, layer index)
  ! @author zander mausolff - adapted from o'donncha's
  ! @date 8/29/2019
  !---------------------------------------------------------------------------!

  subroutine Communicate_3d_Real(partem)

  use variables_mpi
  use mpi

  Implicit none

  !***passed in variables
  !      dimension  partem(LCM,kcm)
  Real(4),intent(inout) ::partem(LCM, kcm)

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)
    
  Real,save,allocatable,dimension(:)::dsendw
  Real,save,allocatable,dimension(:)::dsende
  Real,save,allocatable,dimension(:)::drecve
  Real,save,allocatable,dimension(:)::drecvw
  Real,save,allocatable,dimension(:)::dsendn
  Real,save,allocatable,dimension(:)::dsends
  Real,save,allocatable,dimension(:)::drecvn
  Real,save,allocatable,dimension(:)::drecvs

  !east_west_size   = (jc-4)*kc*2 !*** Original
  east_west_size   = max_width_y*kcm*4  !*** New
  north_south_size = max_width_x*kcm*4
  
  IF(.not.allocated(dsendw))THEN
    allocate(dsendw(east_west_size))
    allocate(dsende(east_west_size))
    allocate(drecve(east_west_size))
    allocate(drecvw(east_west_size))

    allocate(dsendn(north_south_size))
    allocate(dsends(north_south_size))
    allocate(drecvn(north_south_size))
    allocate(drecvs(north_south_size))

    dsendw = 0.0
    dsende = 0.0
    drecve = 0.0
    drecvw = 0.0
    dsendn = 0.0
    dsends = 0.0
    drecvn = 0.0
    drecvs = 0.0
  ENDIF

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsendw, length_arg, mpi_real, nbr_west, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(DRECVE, length_arg, mpi_real, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    ii = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecve(ii)
      ENDDO
    ENDDO
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsende, length_arg, mpi_real, nbr_east, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(DRECVW, length_arg, mpi_real, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVW(II)
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsendn, length_arg, mpi_real, nbr_north, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south  /=  -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(drecvs, length_arg, mpi_real, nbr_south, nbr_south, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    ii = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO K = KSZ(L),KC
        ii = II + 1
        partem(l,k) = drecvs(ii) 
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO K = KSZ(L),KC
        ii = II + 1
        dsends(ii) = partem(l,k)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsends, length_arg, mpi_real, nbr_south, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north  /=  -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL MPI_RECV(drecvn, length_arg, mpi_real, nbr_north, nbr_north, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    ii = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecvn(ii)
      ENDDO
    ENDDO
  ENDIF


  End Subroutine Communicate_3d_Real

  !-----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! @details This routine exchanges data for values for arrays with
  !!         the indexing: (L, K) i.e. ==> (cell index, layer index)
  ! @author zander mausolff - adapted from o'donncha's
  ! @date 8/29/2019
  !---------------------------------------------------------------------------!

  subroutine Communicate_3d_Real8(partem)

  use variables_mpi
  use mpi

  Implicit none

  !***passed in variables
  !      dimension  partem(LCM,kcm)
  Real(8),intent(inout) ::partem(LCM, kcm)

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)
    
  Real,save,allocatable,dimension(:)::dsendw
  Real,save,allocatable,dimension(:)::dsende
  Real,save,allocatable,dimension(:)::drecve
  Real,save,allocatable,dimension(:)::drecvw
  Real,save,allocatable,dimension(:)::dsendn
  Real,save,allocatable,dimension(:)::dsends
  Real,save,allocatable,dimension(:)::drecvn
  Real,save,allocatable,dimension(:)::drecvs

  !east_west_size   = (jc-4)*kc*2 !*** Original
  east_west_size   = max_width_y*kcm*4  !*** New
  north_south_size = max_width_x*kcm*4
  
  IF(.not.allocated(dsendw))THEN
    allocate(dsendw(east_west_size))
    allocate(dsende(east_west_size))
    allocate(drecve(east_west_size))
    allocate(drecvw(east_west_size))

    allocate(dsendn(north_south_size))
    allocate(dsends(north_south_size))
    allocate(drecvn(north_south_size))
    allocate(drecvs(north_south_size))

    dsendw = 0.0
    dsende = 0.0
    drecve = 0.0
    drecvw = 0.0
    dsendn = 0.0
    dsends = 0.0
    drecvn = 0.0
    drecvs = 0.0

  ENDIF

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsendw, length_arg, mpi_real8, nbr_west, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(DRECVE, length_arg, mpi_real8, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    ii = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecve(ii)
      ENDDO
    ENDDO
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsende, length_arg, mpi_real8, nbr_east, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(DRECVW, length_arg, mpi_real8, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVW(II)
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsendn, length_arg, mpi_real8, nbr_north, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south  /=  -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(drecvs, length_arg, mpi_real8, nbr_south, nbr_south, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    ii = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO K = KSZ(L),KC
        ii = II + 1
        partem(l,k) = drecvs(ii) 
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO K = KSZ(L),KC
        ii = II + 1
        dsends(ii) = partem(l,k)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsends, length_arg, mpi_real8, nbr_south, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north  /=  -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL MPI_RECV(drecvn, length_arg, mpi_real8, nbr_north, nbr_north, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    ii = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecvn(ii)
      ENDDO
    ENDDO
  ENDIF

  End Subroutine Communicate_3d_Real8

  !-----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! @details This routine exchanges data for values for arrays with
  !!         the indexing: (L, K) i.e. ==> (cell index, layer index)
  ! @author zander mausolff - adapted from o'donncha's
  ! @date 8/29/2019
  !---------------------------------------------------------------------------!

  subroutine Communicate_3d_Integer(partem)

  use variables_mpi
  use mpi

  Implicit none

  !***passed in variables
  !      dimension  partem(LCM,kcm)
  Integer,intent(inout) ::partem(LCM, kcm)

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)

  Integer,save,allocatable,dimension(:)::dsendw
  Integer,save,allocatable,dimension(:)::dsende
  Integer,save,allocatable,dimension(:)::drecve
  Integer,save,allocatable,dimension(:)::drecvw
  Integer,save,allocatable,dimension(:)::dsendn
  Integer,save,allocatable,dimension(:)::dsends
  Integer,save,allocatable,dimension(:)::drecvn
  Integer,save,allocatable,dimension(:)::drecvs

  IF(.not.allocated(dsendw))THEN
    allocate(dsendw((JC)*kc*2))
    allocate(dsende((JC)*kc*2))
    allocate(drecve((JC)*kc*2))
    allocate(drecvw((JC)*kc*2))

    allocate(dsendn(IC*kc*2))
    allocate(dsends(IC*kc*2))
    allocate(drecvn(IC*kc*2))
    allocate(drecvs(IC*kc*2))

    dsendw = 0
    dsende = 0
    drecve = 0
    drecvw = 0
    dsendn = 0
    dsends = 0
    drecvn = 0
    drecvs = 0
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    ierr = 0
    CALL MPI_SEND(dsendw, length_arg, MPI_Integer, nbr_west, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II 

    IERR = 0
    CALL MPI_RECV(drecve, length_arg, MPI_Integer, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVE(II)
      ENDDO
    ENDDO
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsende, length_arg, MPI_Integer, nbr_east, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II   !  (jc-4)*2*KC
    
    IERR = 0
    CALL MPI_RECV(DRECVW, length_arg, MPI_Integer, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVW(II)
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_SEND(dsendn, length_arg, MPI_Integer, nbr_north, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south  /=  -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_RECV(drecvs, length_arg, MPI_Integer, nbr_south, nbr_south, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVS(II) 
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO K = KSZ(L),KC
        II = II + 1
        DSENDS(II) = PARTEM(L,K)
      ENDDO
    ENDDO
    LENGTH_ARG = II
    
    IERR = 0
    CALL MPI_SEND(dsends, length_arg, MPI_Integer, nbr_south, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north  /=  -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    LENGTH_ARG = II
    
    IERR = 0
    CALL MPI_RECV(drecvn, length_arg, MPI_Integer, nbr_north, nbr_north, comm_2d, status_message, IERR)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVN(II)
      ENDDO
    ENDDO
  ENDIF

  end subroutine Communicate_3d_Integer


  SUBROUTINE Communicate_4D_Real(partem, n_dim)

  Use GLOBAL
  Use MPI
  Use Variables_MPI

  Implicit None

  !***Passed in variables
  REAL(4),intent(inout)  :: PARTEM(LCM, KCM, n_dim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS=23)
  Integer, intent(in) :: n_dim !< Dimension of the third array index

  !***Local variables
  Integer :: i, j, k, II, L, NW
  Integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)
  INTEGER,SAVE :: MAXLEN = -1
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS_4D

  IF(.NOT.ALLOCATED(DSENDW_4D) .OR. n_dim > MAXLEN )THEN
    NWE = JC*KC*2*n_dim
    IF( MAXLEN > -1 )THEN
      DEALLOCATE(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      DEALLOCATE(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    ENDIF
    ALLOCATE(DSENDW_4D(NWE))
    ALLOCATE(DSENDE_4D(NWE))
    ALLOCATE(DRECVE_4D(NWE))
    ALLOCATE(DRECVW_4D(NWE))
    NNS = IC*KC*2*n_dim
    ALLOCATE(DSENDN_4D(NNS))
    ALLOCATE(DSENDS_4D(NNS))
    ALLOCATE(DRECVN_4D(NNS))
    ALLOCATE(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = n_dim
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(DSENDW_4D, length_arg, mpi_real, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVE_4D, length_arg, mpi_real, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    IERR = 0
    CALL MPI_SEND(DSENDE_4D, length_arg, mpi_real, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVW_4D, length_arg, mpi_real, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN 
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(DSENDN_4D, length_arg, mpi_real, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVS_4D, length_arg, mpi_real, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        ENDDO
      ENDDO
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(DSENDS_4D, length_arg, mpi_real, nbr_south, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVN_4D, length_arg, mpi_real, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_4D_Real

  
  
  SUBROUTINE Communicate_4D_Real8(partem, n_dim)

  Use GLOBAL
  Use MPI
  Use Variables_MPI

  Implicit None

  !***Passed in variables
  REAL(8), intent(inout) :: PARTEM(LCM, KCM, n_dim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS=23)
  Integer, intent(in)    :: n_dim !< Dimension of the third array index

  !***Local variables
  Integer :: i, j, k, II, L, NW
  Integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)
  INTEGER,SAVE :: MAXLEN = -1
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS_4D

  IF(.NOT.ALLOCATED(DSENDW_4D) .OR. n_dim > MAXLEN )THEN
    NWE = JC*KC*2*n_dim
    IF( MAXLEN > -1 )THEN
      DEALLOCATE(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      DEALLOCATE(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    ENDIF
    ALLOCATE(DSENDW_4D(NWE))
    ALLOCATE(DSENDE_4D(NWE))
    ALLOCATE(DRECVE_4D(NWE))
    ALLOCATE(DRECVW_4D(NWE))
    NNS = IC*KC*2*n_dim
    ALLOCATE(DSENDN_4D(NNS))
    ALLOCATE(DSENDS_4D(NNS))
    ALLOCATE(DRECVN_4D(NNS))
    ALLOCATE(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = n_dim
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(DSENDW_4D, length_arg, mpi_real8, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVE_4D, length_arg, mpi_real8, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    IERR = 0
    CALL MPI_SEND(DSENDE_4D, length_arg, mpi_real8, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVW_4D, length_arg, mpi_real8, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN 
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II
    
    IERR = 0
    CALL MPI_SEND(DSENDN_4D, length_arg, mpi_real8, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVS_4D, length_arg, mpi_real8, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        ENDDO
      ENDDO
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(DSENDS_4D, length_arg, mpi_real8, nbr_south, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVN_4D, length_arg, mpi_real8, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_4D_Real8

  ! @details geneic communication routine

  SUBROUTINE Communicate_4D_Integer(partem, n_dim)

  Use GLOBAL
  Use MPI
  Use Variables_MPI

  Implicit None

  !***Passed in variables
  Integer,intent(inout)  :: PARTEM(LCM, KCM, n_dim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS=23)
  Integer, intent(in)    :: n_dim !< Dimension of the third array index

  !***Local variables
  Integer :: i, j, k, II, L, NW, NWE, NNS
  Integer :: length_arg  !< LENGTH OF VECTOR TO BE COMMUNICATED
  INTEGER,SAVE :: MAXLEN = -1
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)

  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN_4D
  Integer,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS_4D

  IF(.NOT.ALLOCATED(DSENDW_4D) .OR. n_dim > MAXLEN )THEN
    NWE = JC*KC*2*n_dim
    IF( MAXLEN > -1 )THEN
      DEALLOCATE(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      DEALLOCATE(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    ENDIF
    ALLOCATE(DSENDW_4D(NWE))
    ALLOCATE(DSENDE_4D(NWE))
    ALLOCATE(DRECVE_4D(NWE))
    ALLOCATE(DRECVW_4D(NWE))
    NNS = IC*KC*2*n_dim
    ALLOCATE(DSENDN_4D(NNS))
    ALLOCATE(DSENDS_4D(NNS))
    ALLOCATE(DRECVN_4D(NNS))
    ALLOCATE(DRECVS_4D(NNS))

    DSENDW_4D = 0
    DSENDE_4D = 0
    DRECVE_4D = 0
    DRECVW_4D = 0
    DSENDN_4D = 0
    DSENDS_4D = 0
    DRECVN_4D = 0
    DRECVS_4D = 0
    MAXLEN = n_dim
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II     ! (jc-4)*2*KC*n_dim    ! *** Length of vector for East-West Communication
    
    IERR = 0
    CALL MPI_SEND(DSENDW_4D, length_arg, MPI_Integer, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II     ! (jc-4)*2*KC*n_dim    ! *** Length of vector for East-West Communication
    ierr = 0
    CALL MPI_RECV(DRECVE_4D, length_arg, MPI_Integer, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    IERR = 0
    CALL MPI_SEND(DSENDE_4D, length_arg, MPI_Integer, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    ierr = 0
    CALL MPI_RECV(DRECVW_4D, length_arg, MPI_Integer, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN 
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II     ! IC*2*KC*n_dim    ! *** Length of vector for North-South Communication
    
    IERR = 0
    CALL MPI_SEND(DSENDN_4D, length_arg, MPI_Integer, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVS_4D, length_arg, MPI_Integer, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II) ! unpack data from north
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(DSENDS_4D, length_arg, MPI_Integer, nbr_south, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVN_4D, length_arg, MPI_Integer, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO NW = 1,n_dim
        DO K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_4D_Integer

  SUBROUTINE Communicate_4D0_Real(partem, n_dim, istart)

  Use GLOBAL
  Use MPI
  Use Variables_MPI

  Implicit None

  !***Passed in variables
  REAL(4),intent(inout)  :: PARTEM(LCM, 0:KCM, n_dim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS=23)
  Integer, intent(in) :: n_dim, istart !< Dimension of the third array index

  !***Local variables
  Integer :: i, j, k, II, L, NW
  Integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)
  INTEGER,SAVE :: MAXLEN = -1
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS_4D

  IF(.NOT.ALLOCATED(DSENDW_4D) .OR. n_dim > MAXLEN )THEN
    NWE = JC*(KC+1)*2*n_dim
    IF( MAXLEN > -1 )THEN
      DEALLOCATE(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      DEALLOCATE(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    ENDIF
    ALLOCATE(DSENDW_4D(NWE))
    ALLOCATE(DSENDE_4D(NWE))
    ALLOCATE(DRECVE_4D(NWE))
    ALLOCATE(DRECVW_4D(NWE))
    NNS = IC*(KC+1)*2*n_dim
    ALLOCATE(DSENDN_4D(NNS))
    ALLOCATE(DSENDS_4D(NNS))
    ALLOCATE(DRECVN_4D(NNS))
    ALLOCATE(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = n_dim
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II    ! (jc-4)*2*(KC+1)*n_dim

    IERR = 0
    CALL MPI_SEND(DSENDW_4D, length_arg, mpi_real, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II     !  (jc-4)*2*(KC+1)*n_dim    ! *** Length of vector for East-West Communication
    
    ierr = 0
    CALL MPI_RECV(DRECVE_4D, length_arg, mpi_real, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    IERR = 0
    CALL MPI_SEND(DSENDE_4D, length_arg, mpi_real, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVW_4D, length_arg, mpi_real, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN 
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II    ! IC*2*(KC+1)*n_dim    ! *** Length of vector for North-South Communication
    
    IERR = 0
    CALL MPI_SEND(DSENDN_4D, length_arg, mpi_real, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II    ! IC*2*(KC+1)*n_dim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    CALL MPI_RECV(DRECVS_4D, length_arg, mpi_real, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        ENDDO
      ENDDO
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(DSENDS_4D, length_arg, mpi_real, nbr_south, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVN_4D, length_arg, mpi_real, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_4D0_Real

  SUBROUTINE Communicate_4D0_Real8(partem, n_dim, istart)

  Use GLOBAL
  Use MPI
  Use Variables_MPI

  Implicit None

  !***Passed in variables
  REAL(8), intent(inout) :: PARTEM(LCM, 0:KCM, n_dim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS=23)
  Integer, intent(in)    :: n_dim, istart !< Dimension of the third array index

  !***Local variables
  Integer :: i, j, k, II, L, NW
  Integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  INTEGER(4) :: IERR, status_message(MPI_Status_Size)
  INTEGER,SAVE :: MAXLEN = -1
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVE_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVW_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DSENDS_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVN_4D
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)::DRECVS_4D

  IF(.NOT.ALLOCATED(DSENDW_4D) .OR. n_dim > MAXLEN )THEN
    NWE = JC*(KC+1)*2*n_dim
    IF( MAXLEN > -1 )THEN
      DEALLOCATE(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      DEALLOCATE(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    ENDIF
    ALLOCATE(DSENDW_4D(NWE))
    ALLOCATE(DSENDE_4D(NWE))
    ALLOCATE(DRECVE_4D(NWE))
    ALLOCATE(DRECVW_4D(NWE))
    NNS = IC*(KC+1)*2*n_dim
    ALLOCATE(DSENDN_4D(NNS))
    ALLOCATE(DSENDS_4D(NNS))
    ALLOCATE(DRECVN_4D(NNS))
    ALLOCATE(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = n_dim
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II    ! (jc-4)*2*(KC+1)*n_dim

    IERR = 0
    CALL MPI_SEND(DSENDW_4D, length_arg, mpi_real8, nbr_west, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II     !  (jc-4)*2*(KC+1)*n_dim    ! *** Length of vector for East-West Communication
    
    ierr = 0
    CALL MPI_RECV(DRECVE_4D, length_arg, mpi_real8, nbr_east, nbr_east, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    IERR = 0
    CALL MPI_SEND(DSENDE_4D, length_arg, mpi_real8, nbr_east, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVW_4D, length_arg, mpi_real8, nbr_west, nbr_west, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN 
    II = 0
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II    ! IC*2*(KC+1)*n_dim    ! *** Length of vector for North-South Communication
    
    IERR = 0
    CALL MPI_SEND(DSENDN_4D, length_arg, mpi_real8, nbr_north, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II    ! IC*2*(KC+1)*n_dim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    CALL MPI_RECV(DRECVS_4D, length_arg, mpi_real8, nbr_south, nbr_south, comm_2d, status_message, IERR)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        ENDDO
      ENDDO
    ENDDO
    length_arg = II

    IERR = 0
    CALL MPI_SEND(DSENDS_4D, length_arg, mpi_real8, nbr_south, process_id, comm_2d, IERR)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC+1)*n_dim
    ENDDO
    length_arg = II
    
    ierr = 0
    CALL MPI_RECV(DRECVN_4D, length_arg, mpi_real8, nbr_north, nbr_north, comm_2d, status_message, IERR)
    
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO NW = 1,n_dim
        DO K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  END SUBROUTINE Communicate_4D0_Real8

End module Communicate_Ghost_Routines
