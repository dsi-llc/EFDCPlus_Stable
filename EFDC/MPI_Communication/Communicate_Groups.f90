! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Communicates pressure boundary values amongs processors
! @author Zander Mausolff
! @date 9/4/2019
! @parameters pardata(LCM)

! ********************************************************************************************
! *** Sending values should always be from the active cells
! *** Receiving values should always be placed into the ghost cells
! ********************************************************************************************
  
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! @details Sets up the global active cell list for each direction
! @author Paul Craig
! @date 2020-03-23

  SUBROUTINE Communicate_Initialize()

  use Variables_MPI
  use MPI
  use Mod_DSI_SendRecv

  implicit none

  !***Local variables
  integer :: I, J, II, L

  !                 4
  !             nbr_north
  !          |-------------|
  !          |             |
  !          |             |
  ! nbr_west | process_id  |  nbr_east
  !     1    |             |     2
  !          |             |
  !          |-------------|
  !             nbr_south
  !                 3
  
  II = MAX(max_width_y,max_width_x)*2
  
  ! *** Comm_Cells(d1,d2,d3) - Active cell list
  ! ***                   d1 - Number of active cells
  ! ***                   d2 - Active (1) or Ghost cells (2)
  ! ***                   d3 - Face index, 1-West, 2-East, 3-South, 4-North
  
  allocate(Comm_Cells(ii,2,4))
  allocate(nComm_Cells(2,4))
  Comm_Cells = 0
  nComm_Cells = 0
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather list of west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = 3,4
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,1,1) = L
        endif
      enddo
    enddo
    nComm_Cells(1,1) = II
  endif

  ! *** Receive and populate list of East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = IC-1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,2,2) = L
        endif
      enddo
    enddo
    nComm_Cells(2,2) = II
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = IC-3,IC-2
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,1,2) = L
        endif
      enddo
    enddo
    nComm_Cells(1,2) = II
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = 1,2
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,2,1) = L
        endif
      enddo
    enddo
    nComm_Cells(2,1) = II
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do J = JC-3,JC-2
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,1,4) = L
        endif
      enddo
    enddo
    nComm_Cells(1,4) = II
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do J = 1,2
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,2,3) = L
        endif
      enddo
    enddo
    nComm_Cells(2,3) = II
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do J = 3,4
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,1,3) = L
        endif
      enddo
    enddo
    nComm_Cells(1,3) = II
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do J = JC-1,JC
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + 1
          Comm_Cells(II,2,4) = L
        endif
      enddo
    enddo
    nComm_Cells(2,4) = II
  endif

END SUBROUTINE Communicate_Initialize
  
  
Subroutine Communicate_PUV1()

  use GLOBAL
  use Variables_MPI
  use MPI
  use Mod_DSI_SendRecv

  implicit none

  !***Local variables
  integer :: I, J, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values

  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*2*4     ! *** 2 Columns and 4 Variables
  north_south_size = max_width_x*2*4     ! *** 2 Rows    and 4 Variables

  if(.not.ALLOCATED(PSENDW) )then
    allocate(PSENDW(east_west_size))
    allocate(PSENDE(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))
    
    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = SUB(L)
      II = II + 1
      PSENDW(II) = SVB(L)
      II = II + 1
      PSENDW(II) = SBX(L)
      II = II + 1
      PSENDW(II) = SBY(L)
    enddo
    length_arg = II

    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 4
    enddo
    length_arg = II
    call DSI_RECV(PRECVE, length_arg, nbr_east)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      SUB(L) = PRECVE(II)
      II = II + 1
      SVB(L) = PRECVE(II)
      II = II + 1
      SBX(L) = PRECVE(II)
      II = II + 1
      SBY(L) = PRECVE(II)
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
      PSENDE(II) = SUB(L)
      II = II + 1
      PSENDE(II) = SVB(L)
      II = II + 1
      PSENDE(II) = SBX(L)
      II = II + 1
      PSENDE(II) = SBY(L)
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 4
    enddo
    length_arg = II
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      SUB(L) = PRECVW(II)
      II = II + 1
      SVB(L) = PRECVW(II)
      II = II + 1
      SBX(L) = PRECVW(II)
      II = II + 1
      SBY(L) = PRECVW(II)
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
      PSENDN(II) = SUB(L)
      II = II + 1
      PSENDN(II) = SVB(L)
      II = II + 1
      PSENDN(II) = SBX(L)
      II = II + 1
      PSENDN(II) = SBY(L)
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 4
    enddo
    length_arg = II
    call DSI_RECV(PRECVS, length_arg, nbr_south)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      SUB(L) = PRECVS(II)
      II = II + 1
      SVB(L) = PRECVS(II)
      II = II + 1
      SBX(L) = PRECVS(II)
      II = II + 1
      SBY(L) = PRECVS(II)
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
      PSENDS(II) = SUB(L)
      II = II + 1
      PSENDS(II) = SVB(L)
      II = II + 1
      PSENDS(II) = SBX(L)
      II = II + 1
      PSENDS(II) = SBY(L)
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 4
    enddo
    length_arg = II
    call DSI_RECV(PRECVN, length_arg, nbr_north)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      SUB(L) = PRECVN(II)
      II = II + 1
      SVB(L) = PRECVN(II)
      II = II + 1
      SBX(L) = PRECVN(II)
      II = II + 1
      SBY(L) = PRECVN(II)
    enddo
  endif

END SUBROUTINE Communicate_PUV1

! @details Communicates pressure boundary values amongs processors
! @author Zander Mausolff
! @date 9/4/2019
! @parameters pardata(LCM)

! ********************************************************************************************
! *** Sending values should always be from the active cells
! *** Receiving values should always be placed into the ghost cells
! ********************************************************************************************
  
Subroutine Communicate_PUV3()

  use GLOBAL
  use Variables_MPI
  use MPI
  use Mod_DSI_SendRecv

  implicit none

  !***Read in variables
  !Real(4), Intent(inout), dimension(LCM)::  pardata
  !Character(len = *), Intent(in) :: var_name
  
  !***Local variables

  integer :: I, J, K, II, L
  integer :: length_arg, maxdim

  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  if(.not.ALLOCATED(PSENDW) )then
    maxdim = max_width_y*2*6 + max_width_y*2*2*kcm
    allocate(PSENDW(maxdim))
    allocate(PSENDE(maxdim))
    allocate(PRECVE(maxdim))
    allocate(PRECVW(maxdim))
    
    maxdim = max_width_x*2*6 + max_width_x*2*2*kcm
    allocate(PSENDN(maxdim))
    allocate(PSENDS(maxdim))
    allocate(PRECVN(maxdim))
    allocate(PRECVS(maxdim))
    
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = UHDYE(L)
      II = II + 1
      PSENDW(II) = VHDXE(L)
      II = II + 1
      PSENDW(II) = UHE(L)
      II = II + 1
      PSENDW(II) = VHE(L)
      II = II + 1
      PSENDW(II) = HU(L)
      II = II + 1
      PSENDW(II) = HV(L)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(1,1)
        L = Comm_Cells(I,1,1)
      
        do K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = SUB3D(L,K)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = SVB3D(L,K)
        enddo
      enddo  
    endif
    length_arg = II

    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 6
      if( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    enddo
    length_arg = II
    call DSI_RECV(PRECVE, length_arg, nbr_east)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      UHDYE(L) = PRECVE(II)
      II = II + 1
      VHDXE(L) = PRECVE(II)
      II = II + 1
      UHE(L) = PRECVE(II)
      II = II + 1
      VHE(L) = PRECVE(II)
      II = II + 1
      HU(L) = PRECVE(II)
      II = II + 1
      HV(L) = PRECVE(II)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(2,2)
        L = Comm_Cells(I,2,2)
      
        do K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVE(II)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVE(II)
        enddo
      enddo
    endif
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = UHDYE(L)
      II = II + 1
      PSENDE(II) = VHDXE(L)
      II = II + 1
      PSENDE(II) = UHE(L)
      II = II + 1
      PSENDE(II) = VHE(L)
      II = II + 1
      PSENDE(II) = HU(L)
      II = II + 1
      PSENDE(II) = HV(L)
    enddo

    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(1,2)
        L = Comm_Cells(I,1,2)
      
        do K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = SUB3D(L,K)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = SVB3D(L,K)
        enddo
      enddo
    endif
    length_arg = II
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 6
      if( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    enddo
    length_arg = II
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      UHDYE(L) = PRECVW(II)
      II = II + 1
      VHDXE(L) = PRECVW(II)
      II = II + 1
      UHE(L) = PRECVW(II)
      II = II + 1
      VHE(L) = PRECVW(II)
      II = II + 1
      HU(L) = PRECVW(II)
      II = II + 1
      HV(L) = PRECVW(II)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(2,1)
        L = Comm_Cells(I,2,1)
      
        do K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVW(II)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVW(II)
        enddo
      enddo
    endif
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      II = II + 1
      PSENDN(II) = UHDYE(L)
      II = II + 1
      PSENDN(II) = VHDXE(L)
      II = II + 1
      PSENDN(II) = UHE(L)
      II = II + 1
      PSENDN(II) = VHE(L)
      II = II + 1
      PSENDN(II) = HU(L)
      II = II + 1
      PSENDN(II) = HV(L)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(1,4)
        L = Comm_Cells(I,1,4)
      
        do K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = SUB3D(L,K)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = SVB3D(L,K)
        enddo
      enddo
    endif
    length_arg = II
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 6
      if( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    enddo
    length_arg = II
    call DSI_RECV(PRECVS, length_arg, nbr_south)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      UHDYE(L) = PRECVS(II)
      II = II + 1
      VHDXE(L) = PRECVS(II)
      II = II + 1
      UHE(L) = PRECVS(II)
      II = II + 1
      VHE(L) = PRECVS(II)
      II = II + 1
      HU(L) = PRECVS(II)
      II = II + 1
      HV(L) = PRECVS(II)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(2,3)
        L = Comm_Cells(I,2,3)
      
        do K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVS(II)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVS(II)
        enddo
      enddo
    endif
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      II = II + 1
      PSENDS(II) = UHDYE(L)
      II = II + 1
      PSENDS(II) = VHDXE(L)
      II = II + 1
      PSENDS(II) = UHE(L)
      II = II + 1
      PSENDS(II) = VHE(L)
      II = II + 1
      PSENDS(II) = HU(L)
      II = II + 1
      PSENDS(II) = HV(L)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(1,3)
        L = Comm_Cells(I,1,3)
      
        do K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = SUB3D(L,K)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = SVB3D(L,K)
        enddo
      enddo
    endif
    length_arg = II
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 6
      if( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    enddo
    length_arg = II
    call DSI_RECV(PRECVN, length_arg, nbr_north)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      UHDYE(L) = PRECVN(II)
      II = II + 1
      VHDXE(L) = PRECVN(II)
      II = II + 1
      UHE(L) = PRECVN(II)
      II = II + 1
      VHE(L) = PRECVN(II)
      II = II + 1
      HU(L) = PRECVN(II)
      II = II + 1
      HV(L) = PRECVN(II)
    enddo
    if( ISDRY > 0 )then
      do I = 1,nComm_Cells(2,4)
        L = Comm_Cells(I,2,4)
      
        do K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVN(II)
        enddo
        do K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVN(II)
        enddo
      enddo
    endif
  endif

END SUBROUTINE Communicate_PUV3

!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_UVW1()

  use GLOBAL
  use variables_mpi
  use mpi
  use Mod_DSI_SendRecv

  implicit none

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*2*2   ! *** 2 Columns and 2 Variables
  north_south_size = max_width_x*kcm*2*2   ! *** 2 Rows    and 2 Variables
  
  if(.not.allocated(PSENDW) )then
    allocate(PSENDW(east_west_size))
    allocate(PSENDE(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))

    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))

    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = UHDY(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = VHDX(L,K)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*2
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVE(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVE(II)
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
        PSENDE(II) = UHDY(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = VHDX(L,K)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*2
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVW(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVW(II)
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
        PSENDN(II) = UHDY(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = VHDX(L,K)
      enddo
    enddo
    length_arg = II  
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*2
    enddo
    length_arg = II
    
    call DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVS(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVS(II)
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
        PSENDS(II) = UHDY(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = VHDX(L,K)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*2
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVN(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVN(II)
      enddo
    enddo
  endif

end subroutine Communicate_UVW1

Subroutine Communicate_1D2(Var1, Var2)

  use GLOBAL
  use Variables_MPI
  use MPI
  use Mod_DSI_SendRecv

  implicit none
  Real, Intent(inout) :: Var1(LCM), Var2(LCM)
  
  !***Local variables
  integer :: I, J, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values

  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*2*2     ! *** 2 Columns and 2 Variables
  north_south_size = max_width_x*2*2     ! *** 2 Rows    and 2 Variables

  if(.not.ALLOCATED(PSENDW) )then
    allocate(PSENDW(east_west_size))
    allocate(PSENDE(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))
    
    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = Var1(L)
      II = II + 1
      PSENDW(II) = Var2(L)
    enddo
    length_arg = II

    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 2
    enddo
    length_arg = II
    call DSI_RECV(PRECVE, length_arg, nbr_east)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      Var1(L) = PRECVE(II)
      II = II + 1
      Var2(L) = PRECVE(II)
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
      PSENDE(II) = Var1(L)
      II = II + 1
      PSENDE(II) = Var2(L)
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 2
    enddo
    length_arg = II
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      Var1(L) = PRECVW(II)
      II = II + 1
      Var2(L) = PRECVW(II)
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
      PSENDN(II) = Var1(L)
      II = II + 1
      PSENDN(II) = Var2(L)
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 2
    enddo
    length_arg = II
    call DSI_RECV(PRECVS, length_arg, nbr_south)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      Var1(L) = PRECVS(II)
      II = II + 1
      Var2(L) = PRECVS(II)
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
      PSENDS(II) = Var1(L)
      II = II + 1
      PSENDS(II) = Var2(L)
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 2
    enddo
    length_arg = II
    call DSI_RECV(PRECVN, length_arg, nbr_north)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      Var1(L) = PRECVN(II)
      II = II + 1
      Var2(L) = PRECVN(II)
    enddo
  endif

END SUBROUTINE Communicate_1D2
  
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_UVW3()

  use GLOBAL
  use variables_mpi
  use mpi
  use Mod_DSI_SendRecv

  implicit none

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  Real,save,allocatable,dimension(:) :: PSENDE
  Real,save,allocatable,dimension(:) :: PSENDW
  Real,save,allocatable,dimension(:) :: PRECVE
  Real,save,allocatable,dimension(:) :: PRECVW
  
  Real,save,allocatable,dimension(:) :: PSENDN
  Real,save,allocatable,dimension(:) :: PSENDS
  Real,save,allocatable,dimension(:) :: PRECVN
  Real,save,allocatable,dimension(:) :: PRECVS
  
  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*2*4   ! *** 2 Columns and 4 Variables
  north_south_size = max_width_x*kcm*2*4   ! *** 2 Rows    and 4 Variables
  east_west_size   = east_west_size   + max_width_y*(kcm+1)*2*2   ! *** 2 Columns and 2 Variables
  north_south_size = north_south_size + max_width_x*(kcm+1)*2*2   ! *** 2 Rows    and 2 Variables
  
  if(.not.allocated(PSENDW) )then
    allocate(PSENDE(east_west_size))
    allocate(PSENDW(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))

    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))

    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = U(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = V(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = W(L,K)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVE(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVE(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVE(II)
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
        PSENDE(II) = U(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = V(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = W(L,K)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVW(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVW(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVW(II)
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
        PSENDN(II) = U(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = V(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = W(L,K)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVS(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVS(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVS(II)
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
        PSENDS(II) = U(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = V(L,K)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = W(L,K)
      enddo
    enddo
    length_arg = II  
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVN(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVN(II)
      enddo
      do K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVN(II)
      enddo
    enddo
  endif

end subroutine Communicate_UVW3

!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_CON1()

  use GLOBAL
  use variables_mpi
  use mpi
  use Mod_DSI_SendRecv

  implicit none

  !***local variables
  integer :: i, j, k, II, L, IW, NCLASS, IBL, IC1, IC2, NVAR, NS, NT
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  ! *** Water Column
  east_west_size   = max_width_y*kcm*2*NACTIVEWC                   ! *** 2 Columns and NACTIVEWC Variables
  north_south_size = max_width_x*kcm*2*NACTIVEWC                   ! *** 2 Rows    and NACTIVEWC Variables
  NVAR = KC*NACTIVEWC
  
  ! *** Morphology
  !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
  !  NVAR = NVAR + 1
  !  east_west_size   = east_west_size   + max_width_y*2            ! *** 2 Columns
  !  north_south_size = north_south_size + max_width_x*2            ! *** 2 Rows   
  !ENDIF
  
  if(.not.allocated(PSENDS) )then
    allocate(PSENDE(east_west_size))
    allocate(PSENDW(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))

    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))

    PSENDE = 0.0
    PSENDW = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = WCV(IW).VAL0(L,K)
        enddo
      enddo

    enddo

    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = 3,JC-2
    !    do I = 3,4
    !      L  = LIJ(i,j)
    !      if( L > 0 )then
    !        II = II + 1
    !        PSENDS(II) = BELV(L)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
    length_arg = II 
    
    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    enddo
    length_arg = II 
    call DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVE(II)
        enddo
      enddo
    enddo
    
    !! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = 3,JC-2
    !    do I = IC-1,IC
    !      L  = LIJ(i,j)
    !      if( L > 0 )then
    !        II = II + 1
    !        BELV(L) = PRECVE(II)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = WCV(IW).VAL0(L,K)
        enddo
      enddo
    enddo

    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = 3,JC-2
    !    do I = IC-3,IC-2
    !      L = LIJ(I,J)
    !      if( L > 0 )then
    !        II = II + 1
    !        PSENDE(II) = BELV(L)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
    length_arg = II 
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    enddo
    length_arg = II 
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVW(II)
        enddo
      enddo
    enddo
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = 3,JC-2
    !    do I = 1,2
    !      L = LIJ(I,J)
    !      if( L > 0 )then
    !        II = II + 1
    !        BELV(L) = PRECVW(II)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  if( nbr_north /= -1 )then
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = WCV(IW).VAL0(L,K)
        enddo
      enddo
    enddo
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = JC-3,JC-2
    !    do I = 1,IC
    !      L = LIJ(I,J)
    !      if( L > 0 )then
    !        II = II + 1
    !        PSENDN(II) = BELV(L)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
    length_arg = II 
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    enddo
    length_arg = II 
    call DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVS(II)
        enddo
      enddo
    enddo
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = 1,2
    !    do I = 1,IC
    !      L = LIJ(I,J)
    !      if( L > 0 )then
    !        II = II + 1
    !        BELV(L) = PRECVS(II)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
  endif

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  if( nbr_south /= -1 )then
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = WCV(IW).VAL0(L,K)
        enddo
      enddo
    enddo
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = 3,4
    !    do I = 1,IC
    !      L = LIJ(I,J)
    !      if( L > 0 )then
    !        II = II + 1
    !        PSENDS(II) = BELV(L)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
    length_arg = II 
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    enddo
    length_arg = II 
    call DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVN(II)
        enddo
      enddo
    enddo
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .or. ISTRAN(7) > 0) .and. IMORPH > 0 )then
    !  do J = JC-1,JC
    !    do I = 1,IC
    !      L = LIJ(I,J)
    !      if( L > 0 )then
    !        II = II + 1
    !        BELV(L) = PRECVN(II)
    !      endif
    !    enddo
    !  enddo
    !ENDIF
  endif

end subroutine Communicate_CON1
  
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_CON2()

  use GLOBAL
  use variables_mpi
  use mpi
  use Mod_DSI_SendRecv

  implicit none

  !***local variables
  integer :: i, j, k, II, L, IW
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*2*2*NACTIVEWC   ! *** 2 Columns and 2 Variables
  north_south_size = max_width_x*kcm*2*2*NACTIVEWC   ! *** 2 Rows    and 2 Variables
  east_west_size   = east_west_size   + max_width_y*(kcm+1)*2*1*NACTIVEWC   ! *** 2 Columns and 1 Variable
  north_south_size = north_south_size + max_width_x*(kcm+1)*2*1*NACTIVEWC   ! *** 2 Rows    and 1 Variable
  
  if(.not.allocated(PSENDS) )then
    allocate(PSENDE(east_west_size))
    allocate(PSENDW(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))

    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))

    PSENDE = 0.0
    PSENDW = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = FUHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = FVHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = FWUU(L,k,IW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    enddo
    length_arg = II 
    call DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVE(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVE(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVE(II)
        enddo
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
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = FUHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = FVHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = FWUU(L,k,IW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    enddo
    length_arg = II 
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVW(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVW(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVW(II)
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
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = FUHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = FVHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = FWUU(L,k,IW)
        enddo
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    enddo
    length_arg = II 
    call DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVS(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVS(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVS(II)
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
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = FUHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = FVHUD(L,k,IW)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = FWUU(L,k,IW)
        enddo
      enddo
    enddo
    length_arg = II  
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    enddo
    length_arg = II 
    call DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVN(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVN(II)
        enddo
      enddo
      do IW = 1,NACTIVEWC
        do K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVN(II)
        enddo
      enddo
    enddo
  endif

end subroutine Communicate_CON2
  
subroutine Communicate_BEDLOAD(I1, I2)

  use GLOBAL
  use variables_mpi
  use mpi
  use Mod_DSI_SendRecv

  implicit none
  integer, Intent(in) :: I1, I2
  
  !***local variables
  integer :: i, j, k, II, L, NS, nSize
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  real,save,allocatable,dimension(:) :: PSENDE
  real,save,allocatable,dimension(:) :: PSENDW
  real,save,allocatable,dimension(:) :: PSENDN
  real,save,allocatable,dimension(:) :: PSENDS

  real,save,allocatable,dimension(:) :: PRECVE
  real,save,allocatable,dimension(:) :: PRECVW
  real,save,allocatable,dimension(:) :: PRECVN
  real,save,allocatable,dimension(:) :: PRECVS

  if( num_Processors == 1 ) return

  nSize = I2 - I1 + 1
  east_west_size   = max_width_y*2*2*nSize   ! *** 2 Columns and 2 Variables and nSize sediment classes
  north_south_size = max_width_x*2*2*nSize   ! *** 2 Rows    and 2 Variables and nSize sediment classes
  
  if(.not.allocated(PSENDW) )then
    allocate(PSENDW(east_west_size))
    allocate(PSENDE(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))

    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))

    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      do NS = I1,I2
        II = II + 1
        PSENDW(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDW(II) = QSBDLDY(L,NS)
      enddo
    enddo
    length_arg = II
    
    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 2*nSize
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      do NS = I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVE(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVE(II)
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
      
      do NS = I1,I2
        II = II + 1
        PSENDE(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDE(II) = QSBDLDY(L,NS)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 2*nSize
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      do NS = I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVW(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVW(II)
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
      
      do NS = I1,I2
        II = II + 1
        PSENDN(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDN(II) = QSBDLDY(L,NS)
      enddo
    enddo
    length_arg = II  
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 2*nSize
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      do NS = I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVS(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVS(II)
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
      
      do NS = I1,I2
        II = II + 1
        PSENDS(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDS(II) = QSBDLDY(L,NS)
      enddo
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 2*nSize
    enddo
    length_arg = II 
    
    call DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      do NS = I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVN(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVN(II)
      enddo
    enddo
  endif

end subroutine Communicate_BEDLOAD

subroutine Communicate_QQ()

  use GLOBAL
  use variables_mpi
  use mpi
  use Mod_DSI_SendRecv

  implicit none

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  Real,save,allocatable,dimension(:) :: PSENDE
  Real,save,allocatable,dimension(:) :: PSENDW
  Real,save,allocatable,dimension(:) :: PRECVE
  Real,save,allocatable,dimension(:) :: PRECVW
  
  Real,save,allocatable,dimension(:) :: PSENDN
  Real,save,allocatable,dimension(:) :: PSENDS
  Real,save,allocatable,dimension(:) :: PRECVN
  Real,save,allocatable,dimension(:) :: PRECVS
  
  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*2*(kcm*3 + 4)   ! *** 2 Columns and 3 2D Variables and 4 1D Variables
  north_south_size = max_width_x*2*(kcm*3 + 4)   ! *** 2 Rows    and 3 2D Variables and 4 1D Variables
  
  if( .not. allocated(PSENDW) )then
    allocate(PSENDE(east_west_size))
    allocate(PSENDW(east_west_size))
    allocate(PRECVE(east_west_size))
    allocate(PRECVW(east_west_size))

    allocate(PSENDN(north_south_size))
    allocate(PSENDS(north_south_size))
    allocate(PRECVN(north_south_size))
    allocate(PRECVS(north_south_size))

    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  endif

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = QQ(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = QQ(L,K)
      enddo
          
      II = II + 1
      PSENDW(II) = QQL(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = QQL(L,K)
      enddo
          
      II = II + 1
      PSENDW(II) = DML(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = DML(L,K)
      enddo
      
      II = II + 1
      PSENDW(II) = TBX(L)
      II = II + 1
      PSENDW(II) = TBY(L)
      II = II + 1
      PSENDW(II) = UV(L)
      II = II + 1
      PSENDW(II) = VU(L)
      
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDW, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      QQ(L,0) = PRECVE(II)
      do K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVE(II)
      enddo

      II = II + 1
      QQL(L,0) = PRECVE(II)
      do K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVE(II)
      enddo
          
      II = II + 1
      DML(L,0) = PRECVE(II)
      do K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVE(II)
      enddo
      
      II = II + 1
      TBX(L) = PRECVE(II)
      II = II + 1
      TBY(L) = PRECVE(II)
      II = II + 1
      UV(L) = PRECVE(II)
      II = II + 1
      VU(L) = PRECVE(II)
      
    enddo
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = QQ(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = QQ(L,K)
      enddo

      II = II + 1
      PSENDE(II) = QQL(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = QQL(L,K)
      enddo
          
      II = II + 1
      PSENDE(II) = DML(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = DML(L,K)
      enddo

      II = II + 1
      PSENDE(II) = TBX(L)
      II = II + 1
      PSENDE(II) = TBY(L)
      II = II + 1
      PSENDE(II) = UV(L)
      II = II + 1
      PSENDE(II) = VU(L)
      
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDE, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      QQ(L,0) = PRECVW(II)
      do K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVW(II)
      enddo

      II = II + 1
      QQL(L,0) = PRECVW(II)
      do K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVW(II)
      enddo

      II = II + 1
      DML(L,0) = PRECVW(II)
      do K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVW(II)
      enddo
      
      II = II + 1
      TBX(L) = PRECVW(II)
      II = II + 1
      TBY(L) = PRECVW(II)
      II = II + 1
      UV(L) = PRECVW(II)
      II = II + 1
      VU(L) = PRECVW(II)
      
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
      PSENDN(II) = QQ(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = QQ(L,K)
      enddo
          
      II = II + 1
      PSENDN(II) = QQL(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = QQL(L,K)
      enddo
          
      II = II + 1
      PSENDN(II) = DML(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = DML(L,K)
      enddo
      
      II = II + 1
      PSENDN(II) = TBX(L)
      II = II + 1
      PSENDN(II) = TBY(L)
      II = II + 1
      PSENDN(II) = UV(L)
      II = II + 1
      PSENDN(II) = VU(L)
      
    enddo
    length_arg = II 
    
    call DSI_SEND(PSENDN, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      QQ(L,0) = PRECVS(II)
      do K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVS(II)
      enddo
          
      II = II + 1
      QQL(L,0) = PRECVS(II)
      do K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVS(II)
      enddo

      II = II + 1
      DML(L,0) = PRECVS(II)
      do K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVS(II)
      enddo
      
      II = II + 1
      TBX(L) = PRECVS(II)
      II = II + 1
      TBY(L) = PRECVS(II)
      II = II + 1
      UV(L) = PRECVS(II)
      II = II + 1
      VU(L) = PRECVS(II)
      
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
      PSENDS(II) = QQ(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = QQ(L,K)
      enddo
          
      II = II + 1
      PSENDS(II) = QQL(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = QQL(L,K)
      enddo
          
      II = II + 1
      PSENDS(II) = DML(L,0)
      do K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = DML(L,K)
      enddo
      
      II = II + 1
      PSENDS(II) = TBX(L)
      II = II + 1
      PSENDS(II) = TBY(L)
      II = II + 1
      PSENDS(II) = UV(L)
      II = II + 1
      PSENDS(II) = VU(L)      
      
    enddo
    length_arg = II  
    
    call DSI_SEND(PSENDS, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    enddo
    length_arg = II 
    call DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      QQ(L,0) = PRECVN(II)
      do K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVN(II)
      enddo

      II = II + 1
      QQL(L,0) = PRECVN(II)
      do K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVN(II)
      enddo

      II = II + 1
      DML(L,0) = PRECVN(II)
      do K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVN(II)
      enddo
      
      II = II + 1
      TBX(L) = PRECVN(II)
      II = II + 1
      TBY(L) = PRECVN(II)
      II = II + 1
      UV(L) = PRECVN(II)
      II = II + 1
      VU(L) = PRECVN(II)
      
    enddo
  endif

end subroutine Communicate_QQ
  
  
