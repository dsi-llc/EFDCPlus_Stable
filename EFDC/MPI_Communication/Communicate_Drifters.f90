! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Contains send/recv and packagin/unpackaging routines for communicating drifters
!! amongs processes

Module Mod_Communicate_Drifters

! ***
Use MPI
Use GLOBAL
Use Variables_MPI_Drifter
Use Variables_MPI

Implicit None

! *** Local variables

Save

Public :: Communicate_Drifters

! *** Set up the geneic interface for commuicating ghost cells
! *** This will select the the proper array size depending on what is passed in
Interface Communicate_Drifters

  Module Procedure Communicate_Drifters_RK8, &
                   Communicate_Drifters_Integer

End interface

! *** Technical Note:  The STATUS variable returned by MPI_SEND and MPI_RECV is the number of BYTES sent or received.

!***Next are each of the subroutines for the procedure
Contains

!---------------------------------------------------------------------------!
!> @details Communicates drifters that a in ghost cells
!
!> @author Zander Mausolff
!> @date 2/19/2020
!
!> @parm[inout] data_to_comm - passed in drifter datat that is to be communicated
!---------------------------------------------------------------------------!
Subroutine Communicate_Drifters_RK8(data_to_comm)

  Implicit None

  ! *** Dummy variables
  Real(RKD), intent(inout) :: data_to_comm(NPD)

  ! *** Local variables
  Real(RKD), Allocatable :: send_west_data (:)
  Real(RKD), Allocatable :: send_east_data (:)
  Real(RKD), Allocatable :: send_north_data(:)
  Real(RKD), Allocatable :: send_south_data(:)

  Real(RKD), Allocatable :: recv_west_data (:)
  Real(RKD), Allocatable :: recv_east_data (:)
  Real(RKD), Allocatable :: recv_north_data(:)
  Real(RKD), Allocatable :: recv_south_data(:)

  Integer :: ierr, status_message(MPI_Status_Size)

  Allocate(send_west_data (global_max_drifters_to_comm))
  Allocate(send_east_data (global_max_drifters_to_comm))
  Allocate(send_north_data(global_max_drifters_to_comm))
  Allocate(send_south_data(global_max_drifters_to_comm))

  Allocate(recv_west_data (global_max_drifters_to_comm))
  Allocate(recv_east_data (global_max_drifters_to_comm))
  Allocate(recv_north_data(global_max_drifters_to_comm))
  Allocate(recv_south_data(global_max_drifters_to_comm))

  send_west_data  = 0.0
  send_east_data  = 0.0
  send_north_data = 0.0
  send_south_data = 0.0

  recv_west_data  = 0.0
  recv_east_data  = 0.0
  recv_north_data = 0.0
  recv_south_data = 0.0

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** West domain active: Gather west active cells to send to the west
  IF( send_west > 0 )THEN
    Call Package_Drifters_in_Ghost_RK8(data_to_comm, drifter_ids_send_west, send_west_data)

    Call MPI_Send(send_west_data, global_max_drifters_to_comm, mpi_real8, nbr_west, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( recv_east > 0 )then
    Call MPI_RECV(recv_east_data, global_max_drifters_to_comm, mpi_real8, nbr_east, nbr_east, comm_2d, status_message, ierr)
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( send_east > 0 )THEN
    Call Package_Drifters_in_Ghost_RK8(data_to_comm, drifter_ids_send_east, send_east_data)

    Call MPI_Send(send_east_data, global_max_drifters_to_comm, mpi_real8, nbr_east, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( recv_west > 0 )THEN
    CALL MPI_Recv(recv_west_data, global_max_drifters_to_comm, mpi_real8, nbr_west, nbr_west, comm_2d, status_message, ierr)
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( send_north > 0 )THEN
    Call Package_Drifters_in_Ghost_RK8(data_to_comm, drifter_ids_send_north, send_north_data)

    Call MPI_Send(send_north_data, global_max_drifters_to_comm, mpi_real8, nbr_north, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( recv_south > 0  )THEN
    Call MPI_Recv(recv_south_data, global_max_drifters_to_comm, mpi_real8, nbr_south, nbr_south, comm_2d, status_message, ierr)
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( send_south > 0 )THEN
    Call Package_Drifters_in_Ghost_RK8(data_to_comm, drifter_ids_send_south, send_south_data)
    Call MPI_Send(send_south_data, global_max_drifters_to_comm, mpi_real8, nbr_south, process_id, comm_2d, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( recv_north > 0 )then
    Call MPI_Recv(recv_north_data, global_max_drifters_to_comm, mpi_real8, nbr_north, nbr_north, comm_2d, status_message, ierr)
  endif

  ! *** Unpack the received data and put into the correct global drifter
  if( recv_east > 0 )then
    Call Unpack_Recieved_Data_RK8(drifter_ids_recv_east, recv_east, recv_east_data, size(data_to_comm), data_to_comm, 0)
  end if

  if( recv_west > 0 )then
    Call Unpack_Recieved_Data_RK8(drifter_ids_recv_west, recv_west, recv_west_data, size(data_to_comm), data_to_comm, 0)
  end if

  if( recv_south > 0 )then
    Call Unpack_Recieved_Data_RK8(drifter_ids_recv_south, recv_south, recv_south_data, size(data_to_comm), data_to_comm, 0)
  end if

  if( recv_north > 0 )then
    Call Unpack_Recieved_Data_RK8(drifter_ids_recv_north, recv_north, recv_north_data, size(data_to_comm), data_to_comm, 0)
  end if

End subroutine Communicate_Drifters_RK8

!---------------------------------------------------------------------------!
!> @details Communicates integer drifter parameters and handles the case
!! of communicating the drifter map ids
!
!> @author Zander Mausolff
!
!> @parm[inout] data_to_comm - passed in drifter datat that is to be communicated
!---------------------------------------------------------------------------!
Subroutine Communicate_Drifters_Integer(data_to_comm, iLLA)

  Implicit None

  ! *** Dummy variables
  Integer, Intent(inout) :: data_to_comm(NPD)
  Integer, Intent(in), Optional :: iLLA
  
  ! *** Local variables
  Integer, Allocatable :: send_west_data (:)
  Integer, Allocatable :: send_east_data (:)
  Integer, Allocatable :: send_north_data(:)
  Integer, Allocatable :: send_south_data(:)

  Integer, Allocatable :: recv_west_data (:)
  Integer, Allocatable :: recv_east_data (:)
  Integer, Allocatable :: recv_north_data(:)
  Integer, Allocatable :: recv_south_data(:)

  Integer :: ierr, iMap, status_message(MPI_Status_Size)

  Allocate(send_west_data (global_max_drifters_to_comm))
  Allocate(send_east_data (global_max_drifters_to_comm))
  Allocate(send_north_data(global_max_drifters_to_comm))
  Allocate(send_south_data(global_max_drifters_to_comm))

  Allocate(recv_west_data (global_max_drifters_to_comm))
  Allocate(recv_east_data (global_max_drifters_to_comm))
  Allocate(recv_north_data(global_max_drifters_to_comm))
  Allocate(recv_south_data(global_max_drifters_to_comm))

  send_west_data  = 0
  send_east_data  = 0
  send_north_data = 0
  send_south_data = 0

  recv_west_data  = 0
  recv_east_data  = 0
  recv_north_data = 0
  recv_south_data = 0
  
  ! *** Check LLA flag
  if( Present(iLLA) )then
    iMap = 1
  else
    iMap = 0
  endif
  
  ! *** West domain active: Gather west active cells to send to the west
  IF( send_west > 0 )THEN
    ! *** Package up ghost drifter data into send array
    Call Package_Mapping(data_to_comm, drifter_ids_send_west, send_west_data, iMap)

    ! *** Send global mapping info for this drifter
    Call MPI_Send(send_west_data, global_max_drifters_to_comm, MPI_Integer, nbr_west, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate East ghost cells if East is active
  IF( recv_east > 0 )then
    ! *** Recieve global mapping info for this drifter
    Call MPI_Recv(recv_east_data, global_max_drifters_to_comm, MPI_Integer, nbr_east, nbr_east, comm_2d, status_message, ierr)

    if(MPI_DEBUG_FLAG == .TRUE. )then
      Call WriteBreak(mpi_log_unit)
      write(mpi_log_unit,'(a,30i5)') 'INT recv_east_data: ',recv_east_data
    end if
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** East domain active: Gather East active cells to send to the East
  IF( send_east > 0 )THEN
    ! *** Package up ghost drifter data into send array
    Call Package_Mapping(data_to_comm, drifter_ids_send_east, send_east_data, iMap)

    Call MPI_Send(send_east_data, global_max_drifters_to_comm, MPI_Integer, nbr_east, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate West ghost cells if West is active
  IF( recv_west > 0 )THEN
    Call MPI_Recv(recv_west_data, global_max_drifters_to_comm, MPI_Integer, nbr_west, nbr_west, comm_2d, status_message, ierr)

    if(MPI_DEBUG_FLAG == .TRUE. )then
      Call WriteBreak(mpi_log_unit)
      write(mpi_log_unit,'(a,30i5)') 'INT recv_west_data: ',recv_west_data
    end if
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** North domain active: Gather North active cells to send to the North
  IF( send_north > 0 )THEN
    ! *** Package up ghost drifter data into send array
    Call Package_Mapping(data_to_comm, drifter_ids_send_north, send_north_data, iMap)

    Call MPI_Send(send_north_data, global_max_drifters_to_comm, MPI_Integer, nbr_north, process_id, comm_2d, ierr)
  ENDIF

  ! *** Receive and populate South ghost cells if South is active
  IF( recv_south > 0  )THEN
    Call MPI_Recv(recv_south_data, global_max_drifters_to_comm, MPI_Integer, nbr_south, nbr_south, comm_2d, status_message, ierr)

    if(MPI_DEBUG_FLAG == .TRUE. )then
      Call WriteBreak(mpi_log_unit)
      write(mpi_log_unit,'(a,30i5)') 'INT recv_south_data: ',recv_south_data
    end if
  ENDIF

  ! ********************************************************************************************
  ! ********************************************************************************************
  ! *** South domain active: Gather South active cells to send to the South
  IF( send_south > 0 )THEN
    ! *** Package up ghost drifter data into send array
    Call Package_Mapping(data_to_comm, drifter_ids_send_south, send_south_data, iMap)

    Call MPI_Send(send_south_data, global_max_drifters_to_comm, MPI_Integer, nbr_south, process_id, comm_2d, ierr)
  endif

  ! *** Receive and populate North ghost cells if North is active
  if( recv_north > 0 )then
    Call MPI_Recv(recv_north_data, global_max_drifters_to_comm, MPI_Integer, nbr_north, nbr_north, comm_2d, status_message, ierr)

    if(MPI_DEBUG_FLAG == .TRUE. )then
      Call WriteBreak(mpi_log_unit)
      write(mpi_log_unit,'(a,30i5)') 'INT recv_north_data: ',recv_north_data
    end if
  endif

  ! *** Unpack the received data and put into the correct global drifter
  if( recv_east > 0 )then
    Call Unpack_Recieved_Data_Integer(drifter_ids_recv_east, recv_east, recv_east_data, size(data_to_comm), data_to_comm, iMap)
  end if

  if( recv_west > 0 )then
    Call Unpack_Recieved_Data_Integer(drifter_ids_recv_west, recv_west, recv_west_data, size(data_to_comm), data_to_comm, iMap)
  end if

  if( recv_south > 0 )then
    Call Unpack_Recieved_Data_Integer(drifter_ids_recv_south, recv_south, recv_south_data, size(data_to_comm), data_to_comm, iMap)
  end if

  if( recv_north > 0 )then
    Call Unpack_Recieved_Data_Integer(drifter_ids_recv_north, recv_north, recv_north_data, size(data_to_comm), data_to_comm, iMap)
  end if

End Subroutine Communicate_Drifters_Integer

!---------------------------------------------------------------------------!
!> @details Packages data so only the drifter ids that need to be communicated
!! are sent
!
!> @author Zander Mausolff
!
!> @parm[in] data_to_send - all passed in drifter info
!> @parm[in] drifter_mapping - drifter IDs that should be communicated
!> @parm[inout] packaged_data - drifter info for only the ghost cells
!---------------------------------------------------------------------------!
Subroutine Package_Drifters_in_Ghost_RK8(data_to_send, drifter_mapping, packaged_data)

  Implicit None

  ! *** Passed variables
  Real(RKD), intent(in)    :: data_to_send(NPD)
  Integer, intent(in)      :: drifter_mapping(NPD)
  Real(RKD), intent(inout) :: packaged_data(global_max_drifters_to_comm)

  ! *** Local variables
  Integer :: ii
  Integer :: local_drifter_id

  Do ii = 1, global_max_drifters_to_comm
    ! *** Get the local drifter ID that needs to be communicated
    local_drifter_id = drifter_mapping(ii)
    if(local_drifter_id > 0 )then
      ! *** Get the drifters in ghost cells and place them in array that will be sent off
      packaged_data(ii) = data_to_send(local_drifter_id)

      ! *** Flags to indicate drifter has left the current domain
      LLA_Process(local_drifter_id) = -1
      LLA(local_drifter_id) = 1
    end if
  End do

  if(MPI_DEBUG_FLAG == .TRUE. )then
    Call WriteBreak(mpi_log_unit)
    write(mpi_log_unit,'(a,30f15.4)') 'REAL sending: ',packaged_data
  end if

End Subroutine Package_Drifters_in_Ghost_RK8

!---------------------------------------------------------------------------!
!> @details Packages data containing drifter mapping and other integer arrays
!
!> @author Zander Mausolff
!
!> @parm[in] data_to_send - array containing all local drifters values that should be sent
!> @parm[in] drifter_mapping -  mapping that contains the local ghost cell drifter ids
!> @parm[inout] packaged_data - packaged data containing drifters in the ghost cells
!---------------------------------------------------------------------------!
Subroutine Package_Mapping(data_to_send, ghost_drifter_ids, packaged_data, iMap)

  Implicit None

  ! *** Passed variables
  Integer, intent(inout)   :: data_to_send(NPD)
  Integer, intent(in)      :: iMap, ghost_drifter_ids(global_max_drifters_to_comm)
  Integer, intent(inout)   :: packaged_data(global_max_drifters_to_comm)

  ! *** Local variables
  Integer :: ii
  Integer :: local_drifter_id

  Do ii = 1, global_max_drifters_to_comm
    ! *** Get the local drifter ID that needs to be communicated
    local_drifter_id = ghost_drifter_ids(ii)

    if( local_drifter_id > 0 )then
      ! *** Get the drifters in ghost cells and place them in array that will be sent off
      packaged_data(ii) = data_to_send(local_drifter_id)

      ! *** Flags to indicate drifter has left the current domain
      LLA_Process(local_drifter_id) = -1
      if( iMap > 0 )then
        LLA(local_drifter_id) = 1   
      endif
    end if

  End do

  if( MPI_DEBUG_FLAG == .TRUE. )then
    Call WriteBreak(mpi_log_unit)
    write(mpi_log_unit,'(a,30i5)') 'All data :       ',data_to_send
    write(mpi_log_unit,'(a,30i5)') 'INTEGER sending: ',packaged_data
  end if

End Subroutine Package_Mapping

!---------------------------------------------------------------------------!
!> @details Adds recieved data to the end of currently valid drifter IDs
!---------------------------------------------------------------------------!
Subroutine Unpack_Recieved_Data_Integer(mapping_id, size_r, recvd_data, size, local_drifter_data, iMap)

  Implicit None

  ! *** Dummy variables
  Integer, Intent(in)    :: iMap, size_r
  Integer, Intent(in)    :: mapping_id(size_r)

  Integer, Intent(in)    :: recvd_data(size_r)
  Integer, Intent(in)    :: size
  Integer, Intent(inout) :: local_drifter_data(size)

  ! *** Local variables
  Integer :: ii, NP
  Integer :: drifter_val

  Do ii = 1, size_r
    ! *** Add recvd drifter values to the end of the active list of drifters
    ! *** Get the global ID for the new drifter recieved
    drifter_val = recvd_data(ii)
    NP = mapping_id(ii)

    local_drifter_data(NP) = drifter_val
    
    ! *** Update variables to indicate drifter is active
    LLA_Process(NP) = process_id
    JSPD(NP) = 0
    if( iMap > 0 )then
      IF( drifter_val > 1 ) LLA(NP) = Map2Local(drifter_val).LL   ! ** Map to Local
    endif
    
  End do

End Subroutine Unpack_Recieved_Data_Integer

!---------------------------------------------------------------------------!
!> @details Adds recieved data to the end of currently valid drifter IDs
!
!> @author Zander Mausolff
!
!---------------------------------------------------------------------------!
Subroutine Unpack_Recieved_Data_RK8(mapping_id, size_r, recvd_data, size, local_drifter_data, iMap)

  Implicit None

  ! *** Dummy variables
  Integer, Intent(in) :: iMap, size_r
  Integer, Intent(in) :: mapping_id(size_r)

  Real(RKD), Intent(in)    :: recvd_data(size_r)
  Integer, Intent(in)      :: size
  Real(RKD), Intent(inout) :: local_drifter_data(size)

  ! *** Local variables
  Integer :: ii, NP
  Real(RKD) :: drifter_val

  Do ii = 1, size_r
    ! *** Add recvd drifter values to the end of the active list of drifters
    ! *** Get the global ID for the new drifter recieved
    drifter_val = recvd_data(ii)
    NP = mapping_id(ii)

    local_drifter_data(NP) = drifter_val

    LLA_Process(NP) = process_id
    JSPD(NP) = 0
    ! *** iMap is not used for real variables
  End do

End Subroutine Unpack_Recieved_Data_RK8

End Module Mod_Communicate_Drifters
  