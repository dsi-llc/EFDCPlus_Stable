! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! Generic interface for the communication of ghost cells

Module Communicate_Ghost_Routines

  !Use Mod_Communicate_Functions

  implicit none

  save

  Public :: communicate_ghost_cells

  !***Set up the geneic interface for commuicating ghost cells
  !***This will select the the proper array size depending on what is passed in
  interface Communicate_ghost_cells

  Module Procedure Communicate_1D_Real,    &
                   Communicate_1D_Real8,   &
                   Communicate_1D_Integer, &
                   Communicate_1D_Logical, &
                   Communicate_3D_Real,    &
                   Communicate_3D_Real8,   &
                   Communicate_3D_Integer, &
                   Communicate_3D_Logical, &
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
  ! @parm pardata(lcm)

  ! ********************************************************************************************
  ! *** Sending values should always be from the active cells
  ! *** Receiving values should always be placed into the ghost cells
  ! ********************************************************************************************
  
  Subroutine Communicate_1D_Real(pardata, var_name)

  use Variables_MPI
  use MPI

  implicit none

  !***Read in variables
  real(4), Intent(inout), dimension(lcm) ::  pardata
  Character(len = *), Intent(in) :: var_name
  !***Local variables

  integer :: I, J, II, L
  integer :: length_arg
  integer :: ierr, status_message(MPI_Status_Size)

  real,save,allocatable,dimension(:)::PSENDW
  real,save,allocatable,dimension(:)::PSENDE
  real,save,allocatable,dimension(:)::PRECVE
  real,save,allocatable,dimension(:)::PRECVW
  real,save,allocatable,dimension(:)::PSENDN
  real,save,allocatable,dimension(:)::PSENDS
  real,save,allocatable,dimension(:)::PRECVN
  real,save,allocatable,dimension(:)::PRECVS

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(PSENDW) )then
    allocate(PSENDW(max_width_y*2))
    allocate(PSENDE(max_width_y*2))
    allocate(PRECVE(max_width_y*2))
    allocate(PRECVW(max_width_y*2))
    
    allocate(PSENDN(max_width_x*2))
    allocate(PSENDS(max_width_x*2))
    allocate(PRECVN(max_width_x*2))
    allocate(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif
    
  if( MPI_DEBUG_FLAG )then
    call WriteBreak(mpi_comm_unit)
    write(mpi_comm_unit, '(a)') var_name
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = pardata(L)
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(PSENDW, length_arg, mpi_real4, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVE, length_arg, mpi_real4, nbr_east, nbr_east, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      pardata(L) = PRECVE(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDE, length_arg, mpi_real4, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVW, length_arg, mpi_real4, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = PRECVW(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDN, length_arg, mpi_real4, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVS, length_arg, mpi_real4, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = PRECVS(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDS, length_arg, mpi_real4, nbr_south, process_id,  DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVN, length_arg, mpi_real4, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = PRECVN(II)
    enddo
  endif

  END SUBROUTINE Communicate_1D_Real
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  ! ********************************************************************************************
  ! *** Sending values should always be from the active cells
  ! *** Receiving values should always be placed into the ghost cells
  ! ********************************************************************************************
  
  Subroutine Communicate_1D_Real8(pardata, var_name)

  use Variables_MPI
  use MPI

  implicit none

  !***Read in variables
  real(8), Intent(inout), dimension(lcm)::  pardata
  Character(len = *), Intent(in) :: var_name
  !***Local variables

  integer :: I, J, II, L
  integer :: length_arg
  integer :: ierr, status_message(MPI_Status_Size)

  real,save,allocatable,dimension(:)::PSENDW
  real,save,allocatable,dimension(:)::PSENDE
  real,save,allocatable,dimension(:)::PRECVE
  real,save,allocatable,dimension(:)::PRECVW
  real,save,allocatable,dimension(:)::PSENDN
  real,save,allocatable,dimension(:)::PSENDS
  real,save,allocatable,dimension(:)::PRECVN
  real,save,allocatable,dimension(:)::PRECVS

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(PSENDW) )then
    allocate(PSENDW(max_width_y*2))
    allocate(PSENDE(max_width_y*2))
    allocate(PRECVE(max_width_y*2))
    allocate(PRECVW(max_width_y*2))
    
    allocate(PSENDN(max_width_x*2))
    allocate(PSENDS(max_width_x*2))
    allocate(PRECVN(max_width_x*2))
    allocate(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif
    
  if( MPI_DEBUG_FLAG )then
    call WriteBreak(mpi_comm_unit)
    write(mpi_comm_unit, '(a)') var_name
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = pardata(L)
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(PSENDW, length_arg, mpi_real8, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVE, length_arg, mpi_real8, nbr_east, nbr_east, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      pardata(L) = PRECVE(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDE, length_arg, mpi_real8, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVW, length_arg, mpi_real8, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = PRECVW(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDN, length_arg, mpi_real8, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVS, length_arg, mpi_real8, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = PRECVS(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDS, length_arg, mpi_real8, nbr_south, process_id,  DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVN, length_arg, mpi_real8, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = PRECVN(II)
    enddo
  endif

  END SUBROUTINE Communicate_1D_Real8
  !-----------------------------------------------------------------------
  ! @details Communicates lcm boundary values between processors
  ! @author Zander Mausolff
  ! @date 
  ! @parameters pardata(lcm)

  SUBROUTINE Communicate_1D_Integer(pardata, var_name)

  use Variables_MPI
  use MPI

  implicit none

  !***Read in variables
  integer, Intent(inout), dimension(lcm)::  pardata
  Character(len = *), Intent(in) :: var_name
  !***Local variables

  integer :: I, J, II, L
  integer :: length_arg
  integer(4) :: ierr, status_message(MPI_Status_Size)

  integer,save,allocatable,dimension(:)::PSENDW
  integer,save,allocatable,dimension(:)::PSENDE
  integer,save,allocatable,dimension(:)::PRECVE
  integer,save,allocatable,dimension(:)::PRECVW
  integer,save,allocatable,dimension(:)::PSENDN
  integer,save,allocatable,dimension(:)::PSENDS
  integer,save,allocatable,dimension(:)::PRECVN
  integer,save,allocatable,dimension(:)::PRECVS

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(PSENDW) )then
    allocate(PSENDW(max_width_y*2))
    allocate(PSENDE(max_width_y*2))
    allocate(PRECVE(max_width_y*2))
    allocate(PRECVW(max_width_y*2))
    allocate(PSENDN(max_width_x*2))
    allocate(PSENDS(max_width_x*2))
    allocate(PRECVN(max_width_x*2))
    allocate(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  if( MPI_DEBUG_FLAG )then
    call WriteBreak(mpi_comm_unit)
    write(mpi_comm_unit, '(a)') var_name
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = pardata(L)
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(PSENDW, length_arg, mpi_integer, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVE, length_arg, mpi_integer, nbr_east, nbr_east, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      pardata(L) = PRECVE(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDE, length_arg, mpi_integer, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVW, length_arg, mpi_integer, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = PRECVW(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDN, length_arg, mpi_integer, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVS, length_arg, mpi_integer, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = PRECVS(II)
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = pardata(L)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDS, length_arg, mpi_integer, nbr_south, process_id,  DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVN, length_arg, mpi_integer, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = PRECVN(II)
    enddo
  endif

  END SUBROUTINE Communicate_1D_Integer
  
  !-----------------------------------------------------------------------
  ! @details Communicates lcm boundary values between processors
  ! @author Zander Mausolff/Paul M. Craig
  ! @date 
  ! @parameters pardata(lcm)

  SUBROUTINE Communicate_1D_Logical(pardata)

  use Variables_MPI
  use MPI

  implicit none

  !***Read in variables
  logical, Intent(inout), dimension(lcm)::  pardata

  !***Local variables
  integer :: I, J, II, L
  integer :: length_arg
  integer(4) :: ierr, status_message(MPI_Status_Size)

  integer,save,allocatable,dimension(:)::PSENDW
  integer,save,allocatable,dimension(:)::PSENDE
  integer,save,allocatable,dimension(:)::PRECVE
  integer,save,allocatable,dimension(:)::PRECVW
  integer,save,allocatable,dimension(:)::PSENDN
  integer,save,allocatable,dimension(:)::PSENDS
  integer,save,allocatable,dimension(:)::PRECVN
  integer,save,allocatable,dimension(:)::PRECVS

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(PSENDW) )then
    allocate(PSENDW(max_width_y*2))
    allocate(PSENDE(max_width_y*2))
    allocate(PRECVE(max_width_y*2))
    allocate(PRECVW(max_width_y*2))
    allocate(PSENDN(max_width_x*2))
    allocate(PSENDS(max_width_x*2))
    allocate(PRECVN(max_width_x*2))
    allocate(PRECVS(max_width_x*2))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  if( MPI_DEBUG_FLAG )then
    call WriteBreak(mpi_comm_unit)
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = 0
      if( pardata(L) ) PSENDW(II) = 1
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(PSENDW, length_arg, mpi_integer, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVE, length_arg, mpi_integer, nbr_east, nbr_east, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      
      pardata(L) = .false.
      if( PRECVE(II) == 1 ) pardata(L) = .true. 
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = 0
      if( pardata(L) ) PSENDE(II) = 1
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDE, length_arg, mpi_integer, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVW, length_arg, mpi_integer, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      pardata(L) = .false.
      if( PRECVW(II) == 1 ) pardata(L) = .true. 
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = 0
      if( pardata(L) ) PSENDN(II) = 1
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDN, length_arg, mpi_integer, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVS, length_arg, mpi_integer, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      pardata(L) = .false.
      if( PRECVS(II) == 1 ) pardata(L) = .true. 
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = 0
      if( pardata(L) ) PSENDS(II) = 1
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(PSENDS, length_arg, mpi_integer, nbr_south, process_id,  DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
    enddo
    length_arg = II
    call MPI_RECV(PRECVN, length_arg, mpi_integer, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      pardata(L) = .false.
      if( PRECVN(II) == 1 ) pardata(L) = .true. 
    enddo
  endif

  END SUBROUTINE Communicate_1D_Logical

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

  implicit none

  !***passed in variables
  !      dimension  partem(lcm,kcm)
  real(4),intent(inout) ::partem(lcm, kcm)

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
  integer(4) :: ierr, status_message(MPI_Status_Size)
    
  Real,save,allocatable,dimension(:)::dsendw
  Real,save,allocatable,dimension(:)::dsende
  Real,save,allocatable,dimension(:)::drecve
  Real,save,allocatable,dimension(:)::drecvw
  Real,save,allocatable,dimension(:)::dsendn
  Real,save,allocatable,dimension(:)::dsends
  Real,save,allocatable,dimension(:)::drecvn
  Real,save,allocatable,dimension(:)::drecvs

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*4
  north_south_size = max_width_x*kcm*4
  
  if(.not.allocated(dsendw) )then
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
  endif
    
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsendw, length_arg, mpi_real, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVE, length_arg, mpi_real, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    ii = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecve(ii)
      enddo
    enddo
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsende, length_arg, mpi_real, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVW, length_arg, mpi_real, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVW(II)
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsendn, length_arg, mpi_real, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(drecvs, length_arg, mpi_real, nbr_south, nbr_south, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    ii = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do K = KSZ(L),KC
        ii = II + 1
        partem(l,k) = drecvs(ii) 
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    ii = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do K = KSZ(L),KC
        ii = II + 1
        dsends(ii) = partem(l,k)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsends, length_arg, mpi_real, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north  /=  -1 )then
    ii = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    call MPI_RECV(drecvn, length_arg, mpi_real, nbr_north, nbr_north, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    ii = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecvn(ii)
      enddo
    enddo
  endif


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

  implicit none

  !***passed in variables
  !      dimension  partem(lcm,kcm)
  real(8),intent(inout) ::partem(lcm, kcm)

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
  integer(4) :: ierr, status_message(MPI_Status_Size)
    
  Real,save,allocatable,dimension(:)::dsendw
  Real,save,allocatable,dimension(:)::dsende
  Real,save,allocatable,dimension(:)::drecve
  Real,save,allocatable,dimension(:)::drecvw
  Real,save,allocatable,dimension(:)::dsendn
  Real,save,allocatable,dimension(:)::dsends
  Real,save,allocatable,dimension(:)::drecvn
  Real,save,allocatable,dimension(:)::drecvs

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*4
  north_south_size = max_width_x*kcm*4
  
  if(.not.allocated(dsendw) )then
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

  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsendw, length_arg, mpi_real8, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVE, length_arg, mpi_real8, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    ii = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecve(ii)
      enddo
    enddo
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsende, length_arg, mpi_real8, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVW, length_arg, mpi_real8, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVW(II)
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsendn, length_arg, mpi_real8, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(drecvs, length_arg, mpi_real8, nbr_south, nbr_south, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    ii = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do K = KSZ(L),KC
        ii = II + 1
        partem(l,k) = drecvs(ii) 
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    ii = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do K = KSZ(L),KC
        ii = II + 1
        dsends(ii) = partem(l,k)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsends, length_arg, mpi_real8, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north  /=  -1 )then
    ii = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    call MPI_RECV(drecvn, length_arg, mpi_real8, nbr_north, nbr_north, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    ii = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do K = KSZ(L),KC
        ii = ii + 1
        partem(l,k) = drecvn(ii)
      enddo
    enddo
  endif

  End Subroutine Communicate_3d_Real8

  !-----------------------------------------------------------------------
  ! @details This routine exchanges data for values for arrays with
  !!         the indexing: (L, K) i.e. ==> (cell index, layer index)
  ! @author zander mausolff / Paul M. Craig
  ! @date 
  !---------------------------------------------------------------------------!

  subroutine Communicate_3d_Integer(partem)

  use variables_mpi
  use mpi

  implicit none

  !***passed in variables
  !      dimension  partem(lcm,kcm)
  integer,intent(inout) ::partem(lcm, kcm)

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer(4) :: ierr, status_message(MPI_Status_Size)

  integer,save,allocatable,dimension(:)::dsendw
  integer,save,allocatable,dimension(:)::dsende
  integer,save,allocatable,dimension(:)::drecve
  integer,save,allocatable,dimension(:)::drecvw
  integer,save,allocatable,dimension(:)::dsendn
  integer,save,allocatable,dimension(:)::dsends
  integer,save,allocatable,dimension(:)::drecvn
  integer,save,allocatable,dimension(:)::drecvs

  if( num_Processors == 1 ) return

  if(.not.allocated(dsendw) )then
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
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    ii = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(dsendw, length_arg, MPI_Integer, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II 

    ierr = 0
    call MPI_RECV(drecve, length_arg, MPI_Integer, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVE(II)
      enddo
    enddo
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsende, length_arg, MPI_Integer, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II   !  (jc-4)*2*KC
    
    ierr = 0
    call MPI_RECV(DRECVW, length_arg, MPI_Integer, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVW(II)
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    ii = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = PARTEM(L,K)
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsendn, length_arg, MPI_Integer, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(drecvs, length_arg, MPI_Integer, nbr_south, nbr_south, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVS(II) 
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDS(II) = PARTEM(L,K)
      enddo
    enddo
    LENGTH_ARG = II
    
    ierr = 0
    call MPI_SEND(dsends, length_arg, MPI_Integer, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    LENGTH_ARG = II
    
    ierr = 0
    call MPI_RECV(drecvn, length_arg, MPI_Integer, nbr_north, nbr_north, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = DRECVN(II)
      enddo
    enddo
  endif

  end subroutine Communicate_3d_Integer

  !-----------------------------------------------------------------------
  ! @details This routine exchanges data for values for arrays with
  !!         the indexing: (L, K) i.e. ==> (cell index, layer index)
  ! @author zander mausolff / Paul M. Craig
  ! @date 
  !---------------------------------------------------------------------------!

  subroutine Communicate_3d_Logical(partem)

  use variables_mpi
  use mpi

  implicit none

  !***passed in variables
  !      dimension  partem(lcm,kcm)
  logical,intent(inout) ::partem(lcm, kcm)

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer(4) :: ierr, status_message(MPI_Status_Size)

  integer,save,allocatable,dimension(:)::dsendw
  integer,save,allocatable,dimension(:)::dsende
  integer,save,allocatable,dimension(:)::drecve
  integer,save,allocatable,dimension(:)::drecvw
  integer,save,allocatable,dimension(:)::dsendn
  integer,save,allocatable,dimension(:)::dsends
  integer,save,allocatable,dimension(:)::drecvn
  integer,save,allocatable,dimension(:)::drecvs

  if( num_Processors == 1 ) return

  if(.not.allocated(dsendw) )then
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
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    ii = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDW(II) = 0
        if( PARTEM(L,K) ) DSENDW(II) = 1
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(dsendw, length_arg, MPI_Integer, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II 

    ierr = 0
    call MPI_RECV(drecve, length_arg, MPI_Integer, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = .FALSE.
        if( DRECVE(II) == 1 ) PARTEM(L,K) = .TRUE.
      enddo
    enddo
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDE(II) = 0
        if( PARTEM(L,K) ) DSENDE(II) = 1
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsende, length_arg, MPI_Integer, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II   !  (jc-4)*2*KC
    
    ierr = 0
    call MPI_RECV(DRECVW, length_arg, MPI_Integer, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = .FALSE.
        if( DRECVW(II) == 1 ) PARTEM(L,K) = .TRUE.
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    ii = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDN(II) = 0
        if( PARTEM(L,K) ) DSENDN(II) = 1
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(dsendn, length_arg, MPI_Integer, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(drecvs, length_arg, MPI_Integer, nbr_south, nbr_south, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = .FALSE.
        if( DRECVS(II) == 1 ) PARTEM(L,K) = .TRUE.
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do K = KSZ(L),KC
        II = II + 1
        DSENDS(II) = 0
        if( PARTEM(L,K) ) DSENDS(II) = 1
      enddo
    enddo
    LENGTH_ARG = II
    
    ierr = 0
    call MPI_SEND(dsends, length_arg, MPI_Integer, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    LENGTH_ARG = II
    
    ierr = 0
    call MPI_RECV(drecvn, length_arg, MPI_Integer, nbr_north, nbr_north, DSIcomm, status_message, ierr)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do K = KSZ(L),KC
        II = II + 1
        PARTEM(L,K) = .FALSE.
        if( DRECVN(II) == 1 ) PARTEM(L,K) = .TRUE.
      enddo
    enddo
  endif

  end subroutine Communicate_3d_Logical


  SUBROUTINE Communicate_4D_Real(partem, ndim)

  use GLOBAL
  use MPI
  use Variables_MPI

  implicit none

  !***Passed in variables
  integer, intent(in) :: ndim                         !< Dimension of the third array index
  real(4),intent(inout)  :: PARTEM(lcm, kcm, ndim)

  !***Local variables
  integer :: i, j, k, II, L, NW
  integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  integer(4) :: ierr, status_message(MPI_Status_Size)
  integer,save :: MAXLEN = -1
  
  real,save,allocatable,dimension(:)::DSENDW_4D
  real,save,allocatable,dimension(:)::DSENDE_4D
  real,save,allocatable,dimension(:)::DRECVE_4D
  real,save,allocatable,dimension(:)::DRECVW_4D
  real,save,allocatable,dimension(:)::DSENDN_4D
  real,save,allocatable,dimension(:)::DSENDS_4D
  real,save,allocatable,dimension(:)::DRECVN_4D
  real,save,allocatable,dimension(:)::DRECVS_4D

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(DSENDW_4D) .or. ndim > MAXLEN )then
    NWE = JC*KC*2*ndim
    if( MAXLEN > -1 )then
      deallocate(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      deallocate(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    endif
    allocate(DSENDW_4D(NWE))
    allocate(DSENDE_4D(NWE))
    allocate(DRECVE_4D(NWE))
    allocate(DRECVW_4D(NWE))
    NNS = IC*KC*2*ndim
    allocate(DSENDN_4D(NNS))
    allocate(DSENDS_4D(NNS))
    allocate(DRECVN_4D(NNS))
    allocate(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = ndim
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(DSENDW_4D, length_arg, mpi_real, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVE_4D, length_arg, mpi_real, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(DSENDE_4D, length_arg, mpi_real, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVW_4D, length_arg, mpi_real, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        enddo
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then 
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(DSENDN_4D, length_arg, mpi_real, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVS_4D, length_arg, mpi_real, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        enddo
      enddo
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(DSENDS_4D, length_arg, mpi_real, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVN_4D, length_arg, mpi_real, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        enddo
      enddo
    enddo
  endif

  END SUBROUTINE Communicate_4D_Real

  
  
  SUBROUTINE Communicate_4D_Real8(partem, ndim)

  use GLOBAL
  use MPI
  use Variables_MPI

  implicit none

  !***Passed in variables
  real(8), intent(inout) :: PARTEM(lcm, kcm, ndim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS = 23)
  integer, intent(in)    :: ndim !< Dimension of the third array index

  !***Local variables
  integer :: i, j, k, II, L, NW
  integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  integer(4) :: ierr, status_message(MPI_Status_Size)
  integer,save :: MAXLEN = -1
  
  real,save,allocatable,dimension(:)::DSENDW_4D
  real,save,allocatable,dimension(:)::DSENDE_4D
  real,save,allocatable,dimension(:)::DRECVE_4D
  real,save,allocatable,dimension(:)::DRECVW_4D
  real,save,allocatable,dimension(:)::DSENDN_4D
  real,save,allocatable,dimension(:)::DSENDS_4D
  real,save,allocatable,dimension(:)::DRECVN_4D
  real,save,allocatable,dimension(:)::DRECVS_4D

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(DSENDW_4D) .or. ndim > MAXLEN )then
    NWE = JC*KC*2*ndim
    if( MAXLEN > -1 )then
      deallocate(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      deallocate(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    endif
    allocate(DSENDW_4D(NWE))
    allocate(DSENDE_4D(NWE))
    allocate(DRECVE_4D(NWE))
    allocate(DRECVW_4D(NWE))
    NNS = IC*KC*2*ndim
    allocate(DSENDN_4D(NNS))
    allocate(DSENDS_4D(NNS))
    allocate(DRECVN_4D(NNS))
    allocate(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = ndim
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(DSENDW_4D, length_arg, mpi_real8, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVE_4D, length_arg, mpi_real8, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(DSENDE_4D, length_arg, mpi_real8, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVW_4D, length_arg, mpi_real8, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        enddo
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then 
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_SEND(DSENDN_4D, length_arg, mpi_real8, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVS_4D, length_arg, mpi_real8, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        enddo
      enddo
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(DSENDS_4D, length_arg, mpi_real8, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVN_4D, length_arg, mpi_real8, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        enddo
      enddo
    enddo
  endif

  END SUBROUTINE Communicate_4D_Real8

  ! @details geneic communication routine

  SUBROUTINE Communicate_4D_Integer(partem, ndim)

  use GLOBAL
  use MPI
  use Variables_MPI

  implicit none

  !***Passed in variables
  integer,intent(inout)  :: PARTEM(lcm, kcm, ndim)   ! (NCELLS, NLAYERS, NUM_WQ_VARS = 23)
  integer, intent(in)    :: ndim !< Dimension of the third array index

  !***Local variables
  integer :: i, j, k, II, L, NW, NWE, NNS
  integer :: length_arg  !< LENGTH OF VECTOR TO BE COMMUNICATED
  integer,save :: MAXLEN = -1
  integer(4) :: ierr, status_message(MPI_Status_Size)

  integer,save,allocatable,dimension(:)::DSENDW_4D
  integer,save,allocatable,dimension(:)::DSENDE_4D
  integer,save,allocatable,dimension(:)::DRECVE_4D
  integer,save,allocatable,dimension(:)::DRECVW_4D
  integer,save,allocatable,dimension(:)::DSENDN_4D
  integer,save,allocatable,dimension(:)::DSENDS_4D
  integer,save,allocatable,dimension(:)::DRECVN_4D
  integer,save,allocatable,dimension(:)::DRECVS_4D

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(DSENDW_4D) .or. ndim > MAXLEN )then
    NWE = JC*KC*2*ndim
    if( MAXLEN > -1 )then
      deallocate(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      deallocate(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    endif
    allocate(DSENDW_4D(NWE))
    allocate(DSENDE_4D(NWE))
    allocate(DRECVE_4D(NWE))
    allocate(DRECVW_4D(NWE))
    NNS = IC*KC*2*ndim
    allocate(DSENDN_4D(NNS))
    allocate(DSENDS_4D(NNS))
    allocate(DRECVN_4D(NNS))
    allocate(DRECVS_4D(NNS))

    DSENDW_4D = 0
    DSENDE_4D = 0
    DRECVE_4D = 0
    DRECVW_4D = 0
    DSENDN_4D = 0
    DSENDS_4D = 0
    DRECVN_4D = 0
    DRECVS_4D = 0
    MAXLEN = ndim
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II     ! (jc-4)*2*KC*ndim    ! *** Length of vector for East-West Communication
    
    ierr = 0
    call MPI_SEND(DSENDW_4D, length_arg, MPI_Integer, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II     ! (jc-4)*2*KC*ndim    ! *** Length of vector for East-West Communication
    ierr = 0
    call MPI_RECV(DRECVE_4D, length_arg, MPI_Integer, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(DSENDE_4D, length_arg, MPI_Integer, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    ierr = 0
    call MPI_RECV(DRECVW_4D, length_arg, MPI_Integer, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        enddo
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then 
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II     ! IC*2*KC*ndim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    call MPI_SEND(DSENDN_4D, length_arg, MPI_Integer, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVS_4D, length_arg, MPI_Integer, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II) ! unpack data from north
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(DSENDS_4D, length_arg, MPI_Integer, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVN_4D, length_arg, MPI_Integer, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do NW = 1,ndim
        do K = KSZ(L),KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        enddo
      enddo
    enddo
  endif

  END SUBROUTINE Communicate_4D_Integer

  SUBROUTINE Communicate_4D0_Real(partem, ndim, istart)

  use GLOBAL
  use MPI
  use Variables_MPI

  implicit none

  !***Passed in variables
  integer, intent(in) :: ndim                      !< Dimension of the third array index
  integer, intent(in) :: istart                    !< Starting index of KCM index
  real(4),intent(inout)  :: PARTEM(lcm, istart:kcm, ndim)

  !***Local variables
  integer :: i, j, k, II, L, NW
  integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  integer(4) :: ierr, status_message(MPI_Status_Size)
  integer,save :: MAXLEN = -1
  
  real,save,allocatable,dimension(:)::DSENDW_4D
  real,save,allocatable,dimension(:)::DSENDE_4D
  real,save,allocatable,dimension(:)::DRECVE_4D
  real,save,allocatable,dimension(:)::DRECVW_4D
  real,save,allocatable,dimension(:)::DSENDN_4D
  real,save,allocatable,dimension(:)::DSENDS_4D
  real,save,allocatable,dimension(:)::DRECVN_4D
  real,save,allocatable,dimension(:)::DRECVS_4D

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(DSENDW_4D) .or. ndim > MAXLEN )then
    NWE = JC*(KC+1)*2*ndim
    if( MAXLEN > -1 )then
      deallocate(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      deallocate(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    endif
    allocate(DSENDW_4D(NWE))
    allocate(DSENDE_4D(NWE))
    allocate(DRECVE_4D(NWE))
    allocate(DRECVW_4D(NWE))
    NNS = IC*(KC+1)*2*ndim
    allocate(DSENDN_4D(NNS))
    allocate(DSENDS_4D(NNS))
    allocate(DRECVN_4D(NNS))
    allocate(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = ndim
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II    ! (jc-4)*2*(KC+1)*ndim

    ierr = 0
    call MPI_SEND(DSENDW_4D, length_arg, mpi_real, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II     !  (jc-4)*2*(KC+1)*ndim    ! *** Length of vector for East-West Communication
    
    ierr = 0
    call MPI_RECV(DRECVE_4D, length_arg, mpi_real, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(DSENDE_4D, length_arg, mpi_real, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVW_4D, length_arg, mpi_real, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        enddo
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then 
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II    ! IC*2*(KC+1)*ndim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    call MPI_SEND(DSENDN_4D, length_arg, mpi_real, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II    ! IC*2*(KC+1)*ndim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    call MPI_RECV(DRECVS_4D, length_arg, mpi_real, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        enddo
      enddo
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(DSENDS_4D, length_arg, mpi_real, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVN_4D, length_arg, mpi_real, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        enddo
      enddo
    enddo
  endif

  END SUBROUTINE Communicate_4D0_Real

  SUBROUTINE Communicate_4D0_Real8(partem, ndim, istart)

  use GLOBAL
  use MPI
  use Variables_MPI

  implicit none

  !***Passed in variables
  integer, intent(in) :: ndim                      !< Dimension of the third array index
  integer, intent(in) :: istart                    !< Starting index of KCM index
  real(8),intent(inout)  :: PARTEM(lcm, istart:kcm, ndim)

  !***Local variables
  integer :: i, j, k, II, L, NW
  integer :: length_arg, NNS, NWE  !< LENGTH OF VECTOR TO BE Communicated
  integer(4) :: ierr, status_message(MPI_Status_Size)
  integer,save :: MAXLEN = -1
  
  real,save,allocatable,dimension(:)::DSENDW_4D
  real,save,allocatable,dimension(:)::DSENDE_4D
  real,save,allocatable,dimension(:)::DRECVE_4D
  real,save,allocatable,dimension(:)::DRECVW_4D
  real,save,allocatable,dimension(:)::DSENDN_4D
  real,save,allocatable,dimension(:)::DSENDS_4D
  real,save,allocatable,dimension(:)::DRECVN_4D
  real,save,allocatable,dimension(:)::DRECVS_4D

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(DSENDW_4D) .or. ndim > MAXLEN )then
    NWE = JC*(KC+1)*2*ndim
    if( MAXLEN > -1 )then
      deallocate(DSENDW_4D,DSENDE_4D,DRECVE_4D,DRECVW_4D)
      deallocate(DSENDN_4D,DSENDS_4D,DRECVN_4D,DRECVS_4D)
    endif
    allocate(DSENDW_4D(NWE))
    allocate(DSENDE_4D(NWE))
    allocate(DRECVE_4D(NWE))
    allocate(DRECVW_4D(NWE))
    NNS = IC*(KC+1)*2*ndim
    allocate(DSENDN_4D(NNS))
    allocate(DSENDS_4D(NNS))
    allocate(DRECVN_4D(NNS))
    allocate(DRECVS_4D(NNS))

    DSENDW_4D = 0.0
    DSENDE_4D = 0.0
    DRECVE_4D = 0.0
    DRECVW_4D = 0.0
    DSENDN_4D = 0.0
    DSENDS_4D = 0.0
    DRECVN_4D = 0.0
    DRECVS_4D = 0.0
    MAXLEN = ndim
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDW_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II    ! (jc-4)*2*(KC+1)*ndim

    ierr = 0
    call MPI_SEND(DSENDW_4D, length_arg, mpi_real8, nbr_west, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II     !  (jc-4)*2*(KC+1)*ndim    ! *** Length of vector for East-West Communication
    
    ierr = 0
    call MPI_RECV(DRECVE_4D, length_arg, mpi_real8, nbr_east, nbr_east, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVE_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDE_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    ierr = 0
    call MPI_SEND(DSENDE_4D, length_arg, mpi_real8, nbr_east, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVW_4D, length_arg, mpi_real8, nbr_west, nbr_west, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVW_4D(II)
        enddo
      enddo
    enddo
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then 
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDN_4D(II) = PARTEM(L,K,NW)
        enddo
      enddo
    enddo
    length_arg = II    ! IC*2*(KC+1)*ndim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    call MPI_SEND(DSENDN_4D, length_arg, mpi_real8, nbr_north, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II    ! IC*2*(KC+1)*ndim    ! *** Length of vector for North-South Communication
    
    ierr = 0
    call MPI_RECV(DRECVS_4D, length_arg, mpi_real8, nbr_south, nbr_south, DSIcomm, status_message, ierr)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVS_4D(II)
        enddo
      enddo
    enddo
  endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          DSENDS_4D(II) = PARTEM(L,K,NW) 
        enddo
      enddo
    enddo
    length_arg = II

    ierr = 0
    call MPI_SEND(DSENDS_4D, length_arg, mpi_real8, nbr_south, process_id, DSIcomm, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC+1)*ndim
    enddo
    length_arg = II
    
    ierr = 0
    call MPI_RECV(DRECVN_4D, length_arg, mpi_real8, nbr_north, nbr_north, DSIcomm, status_message, ierr)
    
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do NW = 1,ndim
        do K = 0,KC
          II = II + 1
          PARTEM(L,K,NW) = DRECVN_4D(II)
        enddo
      enddo
    enddo
  endif

  END SUBROUTINE Communicate_4D0_Real8

End module Communicate_Ghost_Routines
