! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
!> @details This routine exchanges data for values for arrays with
!!         the indexing: (L, K) i.e. ==> (cell index, layer index)
!> @author zander mausolff - adapted from o'donncha's
!> @date 8/29/2019
!---------------------------------------------------------------------------!

  subroutine communicate_3d_0(partem)

  Use MPI
  Use Variables_mpi
  Use Mod_DSI_SendRecv
  
  Implicit none

  ! *** Passed in variables
  Real,intent(inout) ::partem(LCM, 0:kcm)

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
    
  Real,save,allocatable,dimension(:)::dsendw
  Real,save,allocatable,dimension(:)::dsende
  Real,save,allocatable,dimension(:)::drecve
  Real,save,allocatable,dimension(:)::drecvw
  Real,save,allocatable,dimension(:)::dsendn
  Real,save,allocatable,dimension(:)::dsends
  Real,save,allocatable,dimension(:)::drecvn
  Real,save,allocatable,dimension(:)::drecvs

  !east_west_size   = (jc-4)*(KC+1)*2 !*** Original
  east_west_size   = max_width_y*(kcm+1)*4  !*** New
  north_south_size = max_width_x*(kcm+1)*4
  
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
    DO J = 3,JC-2
      DO I = 3,4
        L  = LIJ(i,j)
        IF( L > 0 )THEN
          DO K = 0,KC
            II = II + 1
            dsendw(II) = partem(L,k)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    length_arg = II   !(jc-4)*2*(KC+1)
    
    CALL DSI_SEND(dsendw, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = IC-1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + KC + 1
        ENDIF
      ENDDO
    ENDDO
    length_arg = II 
    CALL DSI_RECV(drecve, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    DO J = 3,JC-2
      DO I = IC-1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          DO K = 0,KC
            II = II + 1
            partem(L,k) = drecve(II)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( nbr_east /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = IC-3,IC-2
        L = LIJ(I,J)
        IF( L > 0 )THEN
          DO K = 0,KC
            II = II + 1
            dsende(II) = partem(L,k)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    length_arg = II   !(jc-4)*2*(KC+1)
    
    CALL DSI_SEND(dsende, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO J = 3,JC-2
      DO I = 1,2
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + KC + 1
        ENDIF
      ENDDO
    ENDDO
    length_arg = II  ! (jc-4)*2*(KC+1)
    CALL DSI_RECV(DRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    DO J = 3,JC-2
      DO I = 1,2
        L = LIJ(I,J)
        IF( L > 0 )THEN
          DO K = 0,KC
            II = II + 1
            PARTEM(L,K) = DRECVW(II)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
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
          DO K = 0,KC
            II = II + 1
            dsendn(II) = partem(L,k)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    length_arg = II  ! ic*2*(KC+1)
    
    CALL DSI_SEND(dsendn, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south  /=  -1 )THEN
    II = 0
    DO J = 1,2
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + KC + 1
        ENDIF
      ENDDO
    ENDDO
    length_arg = II  !ic*2*(KC+1)
    CALL DSI_RECV(drecvs, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    DO J = 1,2
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          DO K = 0,KC
            II = II + 1
            partem(L,k) = drecvs(II) 
          ENDDO
        ENDIF
      ENDDO
    ENDDO
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
          DO K = 0,KC
            II = II + 1
            dsends(II) = partem(L,k)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    length_arg = II  ! ic*2*(KC+1)
    
    CALL DSI_SEND(dsends, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north  /=  -1 )THEN
    II = 0
    DO J = JC-1,JC
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          II = II + KC + 1
        ENDIF
      ENDDO
    ENDDO
    length_arg = II  ! ic*2*(KC+1)
    CALL DSI_RECV(drecvn, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    DO J = JC-1,JC
      DO I = 1,IC
        L = LIJ(I,J)
        IF( L > 0 )THEN
          DO K = 0,KC
            II = II + 1
            partem(L,k) = drecvn(II)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  end subroutine communicate_3d_0

    !-----------------------------------------------------------------------
  !---------------------------------------------------------------------------!
  !                     EFDC+ Developed by DSI, LLC.
  !---------------------------------------------------------------------------!
  ! @details This routine exchanges data for values for arrays with
  !!         the indexing: (L, K) i.e. ==> (cell index, layer index)
  !! Handles the case for arrays with a zero index starting
  ! @author zander mausolff - adapted from o'donncha's
  ! @date 3/31/2020
  !---------------------------------------------------------------------------!

  subroutine Communicate_3d_Real_Zero_LCM(partem)

  use mpi
  use variables_mpi
  Use Mod_DSI_SendRecv

  Implicit none

  Real,intent(inout) ::partem(0:LCM, kcm)

  !***local variables
  Integer :: i, j, k, II, L
  Integer :: length_arg
  Integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  Integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
    
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
    
    CALL DSI_SEND(dsendw, length_arg, nbr_west)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( nbr_east /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    CALL DSI_RECV(DRECVE, length_arg, nbr_east)
    
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
    
    CALL DSI_SEND(dsende, length_arg, nbr_east)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( nbr_west /= -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    CALL DSI_RECV(DRECVW, length_arg, nbr_west)
    
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
    
    CALL DSI_SEND(dsendn, length_arg, nbr_north)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( nbr_south  /=  -1 )THEN
    II = 0
    DO I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    
    CALL DSI_RECV(drecvs, length_arg, nbr_south)

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
    
    CALL DSI_SEND(dsends, length_arg, nbr_south)
  ENDIF

  ! *** Receive and populate North ghost cells if North is active
  IF( nbr_north  /=  -1 )THEN
    ii = 0
    DO I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)
    ENDDO
    length_arg = II
    CALL DSI_RECV(drecvn, length_arg, nbr_north)

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


  End Subroutine Communicate_3d_Real_Zero_LCM
