! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
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

  subroutine communicate_ghost_3d0(partem)

  use MPI
  use Variables_mpi
  use Mod_DSI_SendRecv
  
  implicit none

  ! *** Passed in variables
  Real,intent(inout) ::partem(LCM, 0:kcm)

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
    
  Real,save,allocatable,dimension(:)::dsendw
  Real,save,allocatable,dimension(:)::dsende
  Real,save,allocatable,dimension(:)::drecve
  Real,save,allocatable,dimension(:)::drecvw
  Real,save,allocatable,dimension(:)::dsendn
  Real,save,allocatable,dimension(:)::dsends
  Real,save,allocatable,dimension(:)::drecvn
  Real,save,allocatable,dimension(:)::drecvs

  if( num_Processors == 1 ) return

  east_west_size   = max_width_y*(kcm+1)*4  
  north_south_size = max_width_x*(kcm+1)*4
  
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
    do J = 3,JC-2
      do I = 3,4
        L  = LIJ(i,j)
        if( L > 0 )then
          do K = 0,KC
            II = II + 1
            dsendw(II) = partem(L,k)
          enddo
        endif
      enddo
    enddo
    length_arg = II   !(jc-4)*2*(KC+1)
    
    call DSI_SEND(dsendw, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = IC-1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + KC + 1
        endif
      enddo
    enddo
    length_arg = II 
    call DSI_RECV(drecve, length_arg, nbr_east)
    
    ! *** Populate ghost cells
    II = 0
    do J = 3,JC-2
      do I = IC-1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          do K = 0,KC
            II = II + 1
            partem(L,k) = drecve(II)
          enddo
        endif
      enddo
    enddo
  Endif
  
  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  if( nbr_east /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = IC-3,IC-2
        L = LIJ(I,J)
        if( L > 0 )then
          do K = 0,KC
            II = II + 1
            dsende(II) = partem(L,k)
          enddo
        endif
      enddo
    enddo
    length_arg = II   !(jc-4)*2*(KC+1)
    
    call DSI_SEND(dsende, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do J = 3,JC-2
      do I = 1,2
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + KC + 1
        endif
      enddo
    enddo
    length_arg = II  ! (jc-4)*2*(KC+1)
    call DSI_RECV(DRECVW, length_arg, nbr_west)
    
    ! *** Populate ghost cells
    II = 0
    do J = 3,JC-2
      do I = 1,2
        L = LIJ(I,J)
        if( L > 0 )then
          do K = 0,KC
            II = II + 1
            PARTEM(L,K) = DRECVW(II)
          enddo
        endif
      enddo
    enddo
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
          do K = 0,KC
            II = II + 1
            dsendn(II) = partem(L,k)
          enddo
        endif
      enddo
    enddo
    length_arg = II  ! ic*2*(KC+1)
    
    call DSI_SEND(dsendn, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south  /=  -1 )then
    II = 0
    do J = 1,2
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + KC + 1
        endif
      enddo
    enddo
    length_arg = II  !ic*2*(KC+1)
    call DSI_RECV(drecvs, length_arg, nbr_south)

    ! *** Populate ghost cells
    II = 0
    do J = 1,2
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          do K = 0,KC
            II = II + 1
            partem(L,k) = drecvs(II) 
          enddo
        endif
      enddo
    enddo
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
          do K = 0,KC
            II = II + 1
            dsends(II) = partem(L,k)
          enddo
        endif
      enddo
    enddo
    length_arg = II  ! ic*2*(KC+1)
    
    call DSI_SEND(dsends, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north  /=  -1 )then
    II = 0
    do J = JC-1,JC
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          II = II + KC + 1
        endif
      enddo
    enddo
    length_arg = II  ! ic*2*(KC+1)
    call DSI_RECV(drecvn, length_arg, nbr_north)

    ! *** Populate ghost cells
    II = 0
    do J = JC-1,JC
      do I = 1,IC
        L = LIJ(I,J)
        if( L > 0 )then
          do K = 0,KC
            II = II + 1
            partem(L,k) = drecvn(II)
          enddo
        endif
      enddo
    enddo
  endif

  end subroutine communicate_ghost_3d0

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

  subroutine Communicate_Ghost_LCM0(partem)

  use mpi
  use variables_mpi
  use Mod_DSI_SendRecv

  implicit none

  Real,intent(inout) ::partem(0:LCM, kcm)

  !***local variables
  integer :: i, j, k, II, L
  integer :: length_arg
  integer :: east_west_size   !< Size of the collapsed 1D array containg ghost cell values
  integer :: north_south_size !< Size of the collapsed 1D array containg ghost cell values
    
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
    
    call DSI_SEND(dsendw, length_arg, nbr_west)
  endif

  ! *** Receive and populate East ghost cells if East is active
  if( nbr_east /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,2)
      L = Comm_Cells(I,2,2)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    call DSI_RECV(DRECVE, length_arg, nbr_east)
    
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
    
    call DSI_SEND(dsende, length_arg, nbr_east)
  endif

  ! *** Receive and populate West ghost cells if West is active
  if( nbr_west /= -1 )then
    II = 0
    do I = 1,nComm_Cells(2,1)
      L = Comm_Cells(I,2,1)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    call DSI_RECV(DRECVW, length_arg, nbr_west)
    
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
    
    call DSI_SEND(dsendn, length_arg, nbr_north)
  endif

  ! *** Receive and populate South ghost cells if South is active
  if( nbr_south  /=  -1 )then
    II = 0
    do I = 1,nComm_Cells(2,3)
      L = Comm_Cells(I,2,3)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    
    call DSI_RECV(drecvs, length_arg, nbr_south)

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
    
    call DSI_SEND(dsends, length_arg, nbr_south)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( nbr_north  /=  -1 )then
    ii = 0
    do I = 1,nComm_Cells(2,4)
      L = Comm_Cells(I,2,4)
      
      II = II + (KC-KSZ(L)+1)
    enddo
    length_arg = II
    call DSI_RECV(drecvn, length_arg, nbr_north)

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


  End Subroutine Communicate_Ghost_LCM0
