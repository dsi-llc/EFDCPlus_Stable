! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
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

  Use Variables_MPI
  Use MPI
  Use Mod_DSI_SendRecv

  Implicit None

  !***Local variables
  Integer :: I, J, II, L

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
  
  ALLOCATE(Comm_Cells(ii,2,4))
  ALLOCATE(nComm_Cells(2,4))
  Comm_Cells = 0
  nComm_Cells = 0
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather list of west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = 3,4
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,1,1) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(1,1) = II
  ENDIF

  ! *** Receive and populate list of East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = IC-1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,2,2) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(2,2) = II
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = IC-3,IC-2
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,1,2) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(1,2) = II
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = 1,2
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,2,1) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(2,1) = II
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO J = JC-3,JC-2
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,1,4) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(1,4) = II
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO J = 1,2
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,2,3) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(2,3) = II
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO J = 3,4
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,1,3) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(1,3) = II
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO J = JC-1,JC
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + 1
          Comm_Cells(II,2,4) = L
        ENDIF
      ENDDO
    ENDDO
    nComm_Cells(2,4) = II
  ENDIF

END SUBROUTINE Communicate_Initialize
  
  
Subroutine Communicate_PUV1()

  Use GLOBAL
  Use Variables_MPI
  Use MPI
  Use Mod_DSI_SendRecv

  Implicit None

  !***Local variables
  Integer :: I, J, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*2*4     ! *** 2 Columns and 4 Variables
  north_south_size = max_width_x*2*4     ! *** 2 Rows    and 4 Variables

  IF(.NOT.ALLOCATED(PSENDW))THEN
    ALLOCATE(PSENDW(east_west_size))
    ALLOCATE(PSENDE(east_west_size))
    ALLOCATE(PRECVE(east_west_size))
    ALLOCATE(PRECVW(east_west_size))
    
    ALLOCATE(PSENDN(north_south_size))
    ALLOCATE(PSENDS(north_south_size))
    ALLOCATE(PRECVN(north_south_size))
    ALLOCATE(PRECVS(north_south_size))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = SUB(L)
      II = II + 1
      PSENDW(II) = SVB(L)
      II = II + 1
      PSENDW(II) = SBX(L)
      II = II + 1
      PSENDW(II) = SBY(L)
    ENDDO
    length_arg = II

    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 4
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      SUB(L) = PRECVE(II)
      II = II + 1
      SVB(L) = PRECVE(II)
      II = II + 1
      SBX(L) = PRECVE(II)
      II = II + 1
      SBY(L) = PRECVE(II)
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
      PSENDE(II) = SUB(L)
      II = II + 1
      PSENDE(II) = SVB(L)
      II = II + 1
      PSENDE(II) = SBX(L)
      II = II + 1
      PSENDE(II) = SBY(L)
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 4
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      SUB(L) = PRECVW(II)
      II = II + 1
      SVB(L) = PRECVW(II)
      II = II + 1
      SBX(L) = PRECVW(II)
      II = II + 1
      SBY(L) = PRECVW(II)
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
      PSENDN(II) = SUB(L)
      II = II + 1
      PSENDN(II) = SVB(L)
      II = II + 1
      PSENDN(II) = SBX(L)
      II = II + 1
      PSENDN(II) = SBY(L)
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 4
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      SUB(L) = PRECVS(II)
      II = II + 1
      SVB(L) = PRECVS(II)
      II = II + 1
      SBX(L) = PRECVS(II)
      II = II + 1
      SBY(L) = PRECVS(II)
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
      PSENDS(II) = SUB(L)
      II = II + 1
      PSENDS(II) = SVB(L)
      II = II + 1
      PSENDS(II) = SBX(L)
      II = II + 1
      PSENDS(II) = SBY(L)
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 4
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      SUB(L) = PRECVN(II)
      II = II + 1
      SVB(L) = PRECVN(II)
      II = II + 1
      SBX(L) = PRECVN(II)
      II = II + 1
      SBY(L) = PRECVN(II)
    ENDDO
  ENDIF

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

  Use GLOBAL
  Use Variables_MPI
  Use MPI
  Use Mod_DSI_SendRecv

  Implicit None

  !***Read in variables
  !Real(4), Intent(inout), DIMENSION(LCM)::  pardata
  !Character(len=*), Intent(in) :: var_name
  
  !***Local variables

  Integer :: I, J, K, II, L
  Integer :: length_arg, maxdim

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  IF(.NOT.ALLOCATED(PSENDW))THEN
    maxdim = max_width_y*2*6 + max_width_y*2*2*kcm
    ALLOCATE(PSENDW(maxdim))
    ALLOCATE(PSENDE(maxdim))
    ALLOCATE(PRECVE(maxdim))
    ALLOCATE(PRECVW(maxdim))
    
    maxdim = max_width_x*2*6 + max_width_x*2*2*kcm
    ALLOCATE(PSENDN(maxdim))
    ALLOCATE(PSENDS(maxdim))
    ALLOCATE(PRECVN(maxdim))
    ALLOCATE(PRECVS(maxdim))
    
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(1,1)
        L = Comm_Cells(I,1,1)
      
        DO K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = SUB3D(L,K)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = SVB3D(L,K)
        ENDDO
      ENDDO  
    ENDIF
    length_arg = II

    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 6
      IF( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(2,2)
        L = Comm_Cells(I,2,2)
      
        DO K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVE(II)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVE(II)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
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
    ENDDO

    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(1,2)
        L = Comm_Cells(I,1,2)
      
        DO K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = SUB3D(L,K)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = SVB3D(L,K)
        ENDDO
      ENDDO
    ENDIF
    length_arg = II
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 6
      IF( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    II = 0
    DO I = 1,nComm_Cells(2,1)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(2,1)
        L = Comm_Cells(I,2,1)
      
        DO K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVW(II)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVW(II)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,4)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(1,4)
        L = Comm_Cells(I,1,4)
      
        DO K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = SUB3D(L,K)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = SVB3D(L,K)
        ENDDO
      ENDDO
    ENDIF
    length_arg = II
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 6
      IF( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(2,3)
        L = Comm_Cells(I,2,3)
      
        DO K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVS(II)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVS(II)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,3)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(1,3)
        L = Comm_Cells(I,1,3)
      
        DO K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = SUB3D(L,K)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = SVB3D(L,K)
        ENDDO
      ENDDO
    ENDIF
    length_arg = II
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 6
      IF( ISDRY > 0 ) II = II + 2*(KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
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
    ENDDO
    IF( ISDRY > 0 )THEN
      DO I = 1,nComm_Cells(2,4)
        L = Comm_Cells(I,2,4)
      
        DO K = KSZ(L),KC
          II = II + 1
          SUB3D(L,K) = PRECVN(II)
        ENDDO
        DO K = KSZ(L),KC
          II = II + 1
          SVB3D(L,K) = PRECVN(II)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

END SUBROUTINE Communicate_PUV3

!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_UVW1()

  use GLOBAL
  use variables_mpi
  use mpi
  Use Mod_DSI_SendRecv

  Implicit none

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*2*2   ! *** 2 Columns and 2 Variables
  north_south_size = max_width_x*kcm*2*2   ! *** 2 Rows    and 2 Variables
  
  IF(.not.allocated(PSENDW))THEN
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
        PSENDW(II) = UHDY(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = VHDX(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*2
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVE(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVE(II)
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
        PSENDE(II) = UHDY(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = VHDX(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*2
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVW(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVW(II)
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
        PSENDN(II) = UHDY(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = VHDX(L,K)
      ENDDO
    ENDDO
    length_arg = II  
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*2
    ENDDO
    length_arg = II
    
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVS(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVS(II)
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
        PSENDS(II) = UHDY(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = VHDX(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*2
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO K = KSZ(L),KC
        II = II + 1
        UHDY(L,K) = PRECVN(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        VHDX(L,K) = PRECVN(II)
      ENDDO
    ENDDO
  ENDIF

end subroutine Communicate_UVW1

Subroutine Communicate_1D2(Var1, Var2)

  Use GLOBAL
  Use Variables_MPI
  Use MPI
  Use Mod_DSI_SendRecv

  Implicit None
  Real, Intent(inout) :: Var1(LCM), Var2(LCM)
  
  !***Local variables
  Integer :: I, J, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*2*2     ! *** 2 Columns and 2 Variables
  north_south_size = max_width_x*2*2     ! *** 2 Rows    and 2 Variables

  IF(.NOT.ALLOCATED(PSENDW))THEN
    ALLOCATE(PSENDW(east_west_size))
    ALLOCATE(PSENDE(east_west_size))
    ALLOCATE(PRECVE(east_west_size))
    ALLOCATE(PRECVW(east_west_size))
    
    ALLOCATE(PSENDN(north_south_size))
    ALLOCATE(PSENDS(north_south_size))
    ALLOCATE(PRECVN(north_south_size))
    ALLOCATE(PRECVS(north_south_size))
    PSENDW = 0.0
    PSENDE = 0.0
    PRECVE = 0.0
    PRECVW = 0.0
    PSENDN = 0.0
    PSENDS = 0.0
    PRECVN = 0.0
    PRECVS = 0.0
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = Var1(L)
      II = II + 1
      PSENDW(II) = Var2(L)
    ENDDO
    length_arg = II

    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 2
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      Var1(L) = PRECVE(II)
      II = II + 1
      Var2(L) = PRECVE(II)
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
      PSENDE(II) = Var1(L)
      II = II + 1
      PSENDE(II) = Var2(L)
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 2
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      Var1(L) = PRECVW(II)
      II = II + 1
      Var2(L) = PRECVW(II)
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
      PSENDN(II) = Var1(L)
      II = II + 1
      PSENDN(II) = Var2(L)
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 2
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      Var1(L) = PRECVS(II)
      II = II + 1
      Var2(L) = PRECVS(II)
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
      PSENDS(II) = Var1(L)
      II = II + 1
      PSENDS(II) = Var2(L)
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 2
    ENDDO
    length_arg = II
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      Var1(L) = PRECVN(II)
      II = II + 1
      Var2(L) = PRECVN(II)
    ENDDO
  ENDIF

END SUBROUTINE Communicate_1D2
  
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_UVW3()

  use GLOBAL
  use variables_mpi
  use mpi
  Use Mod_DSI_SendRecv

  Implicit none

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
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
  
  IF(.not.allocated(PSENDW))THEN
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
        PSENDW(II) = U(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = V(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = W(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVE(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVE(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVE(II)
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
        PSENDE(II) = U(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = V(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = W(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVW(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVW(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVW(II)
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
        PSENDN(II) = U(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = V(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = W(L,K)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVS(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVS(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVS(II)
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
        PSENDS(II) = U(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = V(L,K)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = W(L,K)
      ENDDO
    ENDDO
    length_arg = II  
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO K = KSZ(L),KC
        II = II + 1
        U(L,K) = PRECVN(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        V(L,K) = PRECVN(II)
      ENDDO
      DO K = KSZ(L),KC
        II = II + 1
        W(L,K) = PRECVN(II)
      ENDDO
    ENDDO
  ENDIF

end subroutine Communicate_UVW3

!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_CON1()

  use GLOBAL
  use variables_mpi
  use mpi
  Use Mod_DSI_SendRecv

  Implicit none

  !***local variables
  Integer :: i, j, k, II, L, IW, NCLASS, IBL, IC1, IC2, NVAR, NS, NT
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  ! *** Water Column
  east_west_size   = max_width_y*kcm*2*NACTIVEWC                   ! *** 2 Columns and NACTIVEWC Variables
  north_south_size = max_width_x*kcm*2*NACTIVEWC                   ! *** 2 Rows    and NACTIVEWC Variables
  NVAR = KC*NACTIVEWC
  
  ! *** Morphology
  !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
  !  NVAR = NVAR + 1
  !  east_west_size   = east_west_size   + max_width_y*2            ! *** 2 Columns
  !  north_south_size = north_south_size + max_width_x*2            ! *** 2 Rows   
  !ENDIF
  
  IF(.not.allocated(PSENDS))THEN
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
  ENDIF

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = WCV(IW).VAL0(L,K)
        ENDDO
      ENDDO

    ENDDO

    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = 3,JC-2
    !    DO I = 3,4
    !      L  = LIJ(i,j)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        PSENDS(II) = BELV(L)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
    length_arg = II 
    
    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVE(II)
        ENDDO
      ENDDO
    ENDDO
    
    !! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = 3,JC-2
    !    DO I = IC-1,IC
    !      L  = LIJ(i,j)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        BELV(L) = PRECVE(II)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = WCV(IW).VAL0(L,K)
        ENDDO
      ENDDO
    ENDDO

    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = 3,JC-2
    !    DO I = IC-3,IC-2
    !      L = LIJ(I,J)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        PSENDE(II) = BELV(L)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
    length_arg = II 
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVW(II)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = 3,JC-2
    !    DO I = 1,2
    !      L = LIJ(I,J)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        BELV(L) = PRECVW(II)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( nbr_north /= -1 )THEN
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(1,4)
      L = Comm_Cells(I,1,4)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = WCV(IW).VAL0(L,K)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = JC-3,JC-2
    !    DO I = 1,IC
    !      L = LIJ(I,J)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        PSENDN(II) = BELV(L)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
    length_arg = II 
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVS(II)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = 1,2
    !    DO I = 1,IC
    !      L = LIJ(I,J)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        BELV(L) = PRECVS(II)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( nbr_south /= -1 )THEN
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(1,3)
      L = Comm_Cells(I,1,3)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = WCV(IW).VAL0(L,K)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = 3,4
    !    DO I = 1,IC
    !      L = LIJ(I,J)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        PSENDS(II) = BELV(L)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
    length_arg = II 
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      NVAR = (KC-KSZ(L)+1)*NACTIVEWC
      II = II + NVAR
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    ! *** Water Column
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          WCV(IW).VAL0(L,K) = PRECVN(II)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** Morphology
    !IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. IMORPH > 0 )THEN
    !  DO J = JC-1,JC
    !    DO I = 1,IC
    !      L = LIJ(I,J)
    !      IF( L > 0 )THEN
    !        II = II + 1
    !        BELV(L) = PRECVN(II)
    !      ENDIF
    !    ENDDO
    !  ENDDO
    !ENDIF
  ENDIF

end subroutine Communicate_CON1
  
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!

subroutine Communicate_CON2()

  use GLOBAL
  use variables_mpi
  use mpi
  Use Mod_DSI_SendRecv

  Implicit none

  !***local variables
  Integer :: i, j, k, II, L, IW
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*kcm*2*2*NACTIVEWC   ! *** 2 Columns and 2 Variables
  north_south_size = max_width_x*kcm*2*2*NACTIVEWC   ! *** 2 Rows    and 2 Variables
  east_west_size   = east_west_size   + max_width_y*(kcm+1)*2*1*NACTIVEWC   ! *** 2 Columns and 1 Variable
  north_south_size = north_south_size + max_width_x*(kcm+1)*2*1*NACTIVEWC   ! *** 2 Rows    and 1 Variable
  
  IF(.not.allocated(PSENDS))THEN
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
  ENDIF

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = FUHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = FVHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDW(II) = FWUU(L,k,IW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVE(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVE(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVE(II)
        ENDDO
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
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = FUHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = FVHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDE(II) = FWUU(L,k,IW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVW(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVW(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVW(II)
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
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = FUHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = FVHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDN(II) = FWUU(L,k,IW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVS(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVS(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVS(II)
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
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = FUHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = FVHUD(L,k,IW)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          PSENDS(II) = FWUU(L,k,IW)
        ENDDO
      ENDDO
    ENDDO
    length_arg = II  
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 3*(KC-KSZ(L)+1)*NACTIVEWC
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FUHUD(L,k,IW) = PRECVN(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FVHUD(L,k,IW) = PRECVN(II)
        ENDDO
      ENDDO
      DO IW=1,NACTIVEWC
        DO K = KSZ(L),KC
          II = II + 1
          FWUU(L,k,IW) = PRECVN(II)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

end subroutine Communicate_CON2
  
subroutine Communicate_BEDLOAD(I1, I2)

  use GLOBAL
  use variables_mpi
  use mpi
  Use Mod_DSI_SendRecv

  Implicit none
  Integer, Intent(in) :: I1, I2
  
  !***local variables
  Integer :: i, j, k, II, L, NS, nSize
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSENDS

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PRECVS

  if( num_Processors == 1 ) return

  nSize = I2 - I1 + 1
  east_west_size   = max_width_y*2*2*nSize   ! *** 2 Columns and 2 Variables and nSize sediment classes
  north_south_size = max_width_x*2*2*nSize   ! *** 2 Rows    and 2 Variables and nSize sediment classes
  
  IF(.not.allocated(PSENDW))THEN
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
  ENDIF

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      DO NS=I1,I2
        II = II + 1
        PSENDW(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDW(II) = QSBDLDY(L,NS)
      ENDDO
    ENDDO
    length_arg = II
    
    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 2*nSize
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      DO NS=I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVE(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVE(II)
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
      
      DO NS=I1,I2
        II = II + 1
        PSENDE(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDE(II) = QSBDLDY(L,NS)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 2*nSize
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      DO NS=I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVW(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVW(II)
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
      
      DO NS=I1,I2
        II = II + 1
        PSENDN(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDN(II) = QSBDLDY(L,NS)
      ENDDO
    ENDDO
    length_arg = II  
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 2*nSize
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      DO NS=I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVS(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVS(II)
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
      
      DO NS=I1,I2
        II = II + 1
        PSENDS(II) = QSBDLDX(L,NS)
        II = II + 1
        PSENDS(II) = QSBDLDY(L,NS)
      ENDDO
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 2*nSize
    ENDDO
    length_arg = II 
    
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      DO NS=I1,I2
        II = II + 1
        QSBDLDX(L,NS) = PRECVN(II)
        II = II + 1
        QSBDLDY(L,NS) = PRECVN(II)
      ENDDO
    ENDDO
  ENDIF

end subroutine Communicate_BEDLOAD

subroutine Communicate_QQ()

  use GLOBAL
  use variables_mpi
  use mpi
  Use Mod_DSI_SendRecv

  Implicit none

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containing ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containing ghost cell values
    
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
  
  IF( .not. allocated(PSENDW) )THEN
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
  ENDIF

  !*** nbr_west, nbr_east, nbr_south, nbr_north set in Topology setup routines

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,1)
      L = Comm_Cells(I,1,1)
      
      II = II + 1
      PSENDW(II) = QQ(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = QQ(L,K)
      ENDDO
          
      II = II + 1
      PSENDW(II) = QQL(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = QQL(L,K)
      ENDDO
          
      II = II + 1
      PSENDW(II) = DML(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDW(II) = DML(L,K)
      ENDDO
      
      II = II + 1
      PSENDW(II) = TBX(L)
      II = II + 1
      PSENDW(II) = TBY(L)
      II = II + 1
      PSENDW(II) = UV(L)
      II = II + 1
      PSENDW(II) = VU(L)
      
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDW, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVE, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + 1
      QQ(L,0) = PRECVE(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVE(II)
      ENDDO

      II = II + 1
      QQL(L,0) = PRECVE(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVE(II)
      ENDDO
          
      II = II + 1
      DML(L,0) = PRECVE(II)
      DO K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVE(II)
      ENDDO
      
      II = II + 1
      TBX(L) = PRECVE(II)
      II = II + 1
      TBY(L) = PRECVE(II)
      II = II + 1
      UV(L) = PRECVE(II)
      II = II + 1
      VU(L) = PRECVE(II)
      
    ENDDO
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(1,2)
      L = Comm_Cells(I,1,2)
      
      II = II + 1
      PSENDE(II) = QQ(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = QQ(L,K)
      ENDDO

      II = II + 1
      PSENDE(II) = QQL(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = QQL(L,K)
      ENDDO
          
      II = II + 1
      PSENDE(II) = DML(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDE(II) = DML(L,K)
      ENDDO

      II = II + 1
      PSENDE(II) = TBX(L)
      II = II + 1
      PSENDE(II) = TBY(L)
      II = II + 1
      PSENDE(II) = UV(L)
      II = II + 1
      PSENDE(II) = VU(L)
      
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDE, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + 1
      QQ(L,0) = PRECVW(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVW(II)
      ENDDO

      II = II + 1
      QQL(L,0) = PRECVW(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVW(II)
      ENDDO

      II = II + 1
      DML(L,0) = PRECVW(II)
      DO K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVW(II)
      ENDDO
      
      II = II + 1
      TBX(L) = PRECVW(II)
      II = II + 1
      TBY(L) = PRECVW(II)
      II = II + 1
      UV(L) = PRECVW(II)
      II = II + 1
      VU(L) = PRECVW(II)
      
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
      PSENDN(II) = QQ(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = QQ(L,K)
      ENDDO
          
      II = II + 1
      PSENDN(II) = QQL(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = QQL(L,K)
      ENDDO
          
      II = II + 1
      PSENDN(II) = DML(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDN(II) = DML(L,K)
      ENDDO
      
      II = II + 1
      PSENDN(II) = TBX(L)
      II = II + 1
      PSENDN(II) = TBY(L)
      II = II + 1
      PSENDN(II) = UV(L)
      II = II + 1
      PSENDN(II) = VU(L)
      
    ENDDO
    length_arg = II 
    
    CALL DSI_SEND(PSENDN, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVS, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + 1
      QQ(L,0) = PRECVS(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVS(II)
      ENDDO
          
      II = II + 1
      QQL(L,0) = PRECVS(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVS(II)
      ENDDO

      II = II + 1
      DML(L,0) = PRECVS(II)
      DO K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVS(II)
      ENDDO
      
      II = II + 1
      TBX(L) = PRECVS(II)
      II = II + 1
      TBY(L) = PRECVS(II)
      II = II + 1
      UV(L) = PRECVS(II)
      II = II + 1
      VU(L) = PRECVS(II)
      
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
      PSENDS(II) = QQ(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = QQ(L,K)
      ENDDO
          
      II = II + 1
      PSENDS(II) = QQL(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = QQL(L,K)
      ENDDO
          
      II = II + 1
      PSENDS(II) = DML(L,0)
      DO K = KSZ(L),KC
        II = II + 1
        PSENDS(II) = DML(L,K)
      ENDDO
      
      II = II + 1
      PSENDS(II) = TBX(L)
      II = II + 1
      PSENDS(II) = TBY(L)
      II = II + 1
      PSENDS(II) = UV(L)
      II = II + 1
      PSENDS(II) = VU(L)      
      
    ENDDO
    length_arg = II  
    
    CALL DSI_SEND(PSENDS, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 7 + (KC-KSZ(L)+1)*3
    ENDDO
    length_arg = II 
    CALL DSI_RECV(PRECVN, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + 1
      QQ(L,0) = PRECVN(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQ(L,K) = PRECVN(II)
      ENDDO

      II = II + 1
      QQL(L,0) = PRECVN(II)
      DO K = KSZ(L),KC
        II = II + 1
        QQL(L,K) = PRECVN(II)
      ENDDO

      II = II + 1
      DML(L,0) = PRECVN(II)
      DO K = KSZ(L),KC
        II = II + 1
        DML(L,K) = PRECVN(II)
      ENDDO
      
      II = II + 1
      TBX(L) = PRECVN(II)
      II = II + 1
      TBY(L) = PRECVN(II)
      II = II + 1
      UV(L) = PRECVN(II)
      II = II + 1
      VU(L) = PRECVN(II)
      
    ENDDO
  ENDIF

end subroutine Communicate_QQ
  
  
