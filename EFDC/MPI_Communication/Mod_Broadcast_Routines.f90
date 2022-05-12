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
! @ details Contains routines to broadcast scalars and arrays from the master 
!! to all other processes.  Does so with a generic interface. Assumes the 
!! broadcasts use the new communicator instantiated with by the MPI Topology
!! Adapted from the link below 
!! https://oceans11.lanl.gov/COSIMdownloads/pop/mpi/broadcast.F90
! @author Zander Mausolff
! @date 5/28/2019
!---------------------------------------------------------------------------!

#ifdef _MPI
Module Broadcast_Routines

    Implicit None

    Save

    Public :: Broadcast_Scalar, Broadcast_Array

!***Handles the broadcasting of scalars with different types and precision
    Interface Broadcast_Scalar

     Module Procedure broadcast_scalar_int,    &
                      broadcast_scalar_int8,   &
                      broadcast_scalar_real,   &
                      broadcast_scalar_real8,  &
                      broadcast_scalar_log,    &
                      broadcast_scalar_char
    End Interface

!***Handles the broadcasting of array dimensions (up to 4D)
!   Handles different precision of arrays and integers (single & double)
Interface Broadcast_Array

     Module Procedure broadcast_array_int_1d,     &
                      broadcast_array_int_2d,     &
                      broadcast_array_int_3d,     &
                      broadcast_array_real_1d,    &
                      broadcast_array_real_2d,    &
                      broadcast_array_real_3d,    &
                      broadcast_array_real_4d,    &
                      broadcast_array_int8_3d,    &
                      broadcast_array_real8_1d,   &
                      broadcast_array_real8_2d,   &
                      broadcast_array_real8_3d,   &
                      broadcast_array_real8_4d,   &
                      broadcast_array_char_1d,    &
                      broadcast_array_log_1d,     &
                      broadcast_array_log_2d,     &
                      broadcast_array_log_3d
End Interface

Contains

!---------------------------------------------------------------------------!
! @details Broadcasts a scalar real variable from one processor (root_pe)
!!  to all other processors. This is a specific instance of the generic
!!  broadcast\_scalar Interface.
Subroutine broadcast_scalar_real(scalar, root_pe)

   USE MPI 
   USE Variables_MPI
   
!***Dummy variables
   real(4), intent(in) :: scalar    !< scalar to be broadcast
   integer, intent(in) :: root_pe   !< processor number to broadcast from

!***local variables
   integer  :: ierr  !< local MPI error flag

   call MPI_BCAST(scalar, 1, mpi_real4, root_pe, comm_2d, ierr)
   
   call MPI_BARRIER(comm_2d, ierr)

 End Subroutine broadcast_scalar_real
!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
! @details Broadcasts a scalar real variable from one processor (root_pe)
!!         to all other processors. Handles double precision reals
! 
 Subroutine broadcast_scalar_real8(scalar, root_pe)

    Use MPI ! MPI Fortran include file
    Use Variables_MPI
   
!***Dummy variables
    real(8), intent(in) :: scalar    !< scalar to be broadcast
    integer, intent(in) :: root_pe   !< processor number to broadcast from

!***Local
    Integer  :: ierr  !< local MPI error flag
    
    ! *** Must specify MPI_Real8 
    Call MPI_BCAST(scalar, 1, MPI_REAL8, root_pe, comm_2d, ierr)
    Call MPI_BARRIER(comm_2d, ierr)

 End Subroutine broadcast_scalar_real8
!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
! @details Broadcasts a integer  real variable from one processor (root_pe)
!!         to all other processors. Handles single precision Integer
Subroutine broadcast_scalar_int(scalar, root_pe)

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
!***Dummy variables
   Integer, intent(in)       :: root_pe !< processor number to broadcast from
   Integer(4), intent(inout) :: scalar  !< scalar to be broadcast

!***Local variables
   Integer  :: ierr  !< local MPI error flag

   Call MPI_BCAST(scalar, 1, MPI_Integer4, root_pe, comm_2d, ierr)
   Call MPI_BARRIER(comm_2d, ierr)

 End Subroutine broadcast_scalar_int

!---------------------------------------------------------------------------!
! @details Broadcasts a scalar integer variable from one processor (root_pe)
!!  to all other processors. This is a specific instance of the generic
!!  broadcast\_scalar Interface.
Subroutine broadcast_scalar_int8(scalar, root_pe)

    Use MPI  ! MPI Fortran include file
    Use Variables_MPI
    
!***Dummy variables
    Integer, intent(in)       :: root_pe !< processor number to broadcast from
    Integer(8), intent(inout) :: scalar  !< scalar to be broadcast

!***Local variables
    integer  :: ierr  ! local MPI error flag

   call MPI_BCAST(scalar, 1, MPI_Integer8, root_pe, comm_2d,ierr)
   call MPI_BARRIER(comm_2d, ierr)

End Subroutine broadcast_scalar_int8
!---------------------------------------------------------------------------!
! @details Broadcasts a scalar logical variable from one processor (root_pe)
!!  to all other processors. This is a specific instance of the generic
!!  broadcast\_scalar Interface.
Subroutine broadcast_scalar_log(scalar, root_pe)

    Use MPI  ! MPI Fortran include file
    Use Variables_MPI
   
!***Dummy variables
    Integer, intent(in)    :: root_pe !< processor number to broadcast from
    Logical, intent(inout) :: scalar  !< scalar to be broadcast

!***Local Variables
    Integer :: itmp !< local temporary
    Integer :: ierr !< MPI error flag

!***Maps the logical to 0/1 and then broadcasts that back
   IF( scalar )THEN
     itmp = 1
   else
     itmp = 0
   Endif

   call MPI_BCAST(itmp, 1, MPI_Int, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!***Remaps the integer value recieved back to a logical value
   IF( itmp == 1 )THEN
     scalar = .true.
   else
     scalar = .false.
   Endif

 End Subroutine broadcast_scalar_log
!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
! @details Broadcasts a scalar character variable from one processor (root_pe)
!!  to all other processors. 
Subroutine broadcast_scalar_char(scalar, root_pe)

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
!***Dummy variables
    Integer , intent(in)         :: root_pe !< processor number to broadcast from
    Character (*), intent(inout) :: scalar  !< scalar to be broadcast

!***Local variables
   Integer :: clength !< length of character
   Integer :: ierr    !< MPI error flag

   clength = len(scalar)

   call MPI_BCAST(scalar, clength, MPI_CHARACTER, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

 End Subroutine broadcast_scalar_char
!--------------------------------------------------------------------

! @details Broadcasts a real vector from one processor (root_pe)
!!  to all other processors. This handles single precision Real
Subroutine broadcast_array_real_1d(array, root_pe)

   Use MPI  ! MPI Fortran include file
   Use GLOBAL
   Use Variables_MPI
   
!***Dummy variables
   Integer , intent(in)           :: root_pe !< processor number to broadcast from
   Real(4), dimension(:), intent(in) :: array   !< array to be broadcast

!***Local variables
   Integer :: nelements !< size of array to be broadcast
   Integer :: ierr      !< local MPI error flag
    
   nelements = size(array) ! get the size of the array for broadcasting
   
   Call MPI_BCAST(array, nelements, mpi_real4, root_pe, comm_2d, ierr)
   Call MPI_BARRIER(comm_2d, ierr)

End Subroutine broadcast_array_real_1d
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! @details Broadcasts a real vector from one processor (root_pe)
!!  to all other processors. This is a specific instance of the generic
!!!  broadcast\_array Interface.
Subroutine broadcast_array_real8_1d(array, root_pe)

   Use MPI  ! MPI Fortran include file
   Use GLOBAL
   Use Variables_MPI
   
!***Dummy variables
   Integer , intent(in)           :: root_pe !< processor number to broadcast from
   Real(8), dimension(:), intent(in) :: array   !< array to be broadcast

!***Local variables
   Integer :: nelements !< size of array to be broadcast
   Integer :: ierr      !< local MPI error flag
    
   nelements = size(array)
   call MPI_BCAST(array, nelements, MPI_REAL8, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

End Subroutine broadcast_array_real8_1d
!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
! @details Broadcasts an integer vector from one processor (root_pe)
!!  to all other processors. Handles 1D single preciscion integer array 
Subroutine broadcast_array_int_1d(array, root_pe)

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
    
   Integer, intent(in)                  :: root_pe !< processor number to broadcast from
   Integer, dimension(:), intent(inout) :: array   !< array to be broadcast

!***Local variables
   Integer :: nelements !< size of array to be broadcast
   Integer :: ierr      !< local MPI error flag
  
   nelements = size(array)
   call MPI_BCAST(array, nelements, MPI_Int, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)
       
 End Subroutine broadcast_array_int_1d

!---------------------------------------------------------------------------!
! @details  Broadcasts a character vector from one processor (root_pe)
!!  to all other processors. This is a specific instance of the generic
!!  broadcast\_array Interface.
Subroutine broadcast_array_char_1d(array,  root_pe)

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
!***Dummy variables
   Integer , intent(in)                   :: root_pe !< processor number to broadcast from
   Character, dimension(:), intent(inout) :: array   !< array to be broadcast
!***Local variables
   Integer :: nelements !< size of array to be broadcast
   Integer :: ierr      !< local MPI error flag

   nelements = size(array)

   Call MPI_BCAST(array, nelements, MPI_CHARACTER, root_pe, comm_2d, ierr)
   Call MPI_BARRIER(comm_2d, ierr)

 End Subroutine broadcast_array_char_1d
!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
Subroutine broadcast_array_log_1d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a logical vector from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:
   Integer , intent(in) :: root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:
   Logical , dimension(:), intent(inout) :: array                ! array to be broadcast
   Integer , dimension(:), allocatable :: array_int            ! temporary array for MPI bcast
   Integer  :: nelements         ! size of array to be broadcast
   Integer ::   ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)
   allocate(array_int(nelements))

   where (array)
     array_int = 1
   elsewhere
     array_int = 0
   End where

   call MPI_BCAST(array_int, nelements, MPI_Int, root_pe, &
                  comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

   where (array_int == 1)
     array = .true.
   elsewhere
     array = .false.
   End where

   deallocate(array_int)



 End Subroutine broadcast_array_log_1d
!---------------------------------------------------------------------------!
!---------------------------------------------------------------------------!
 Subroutine broadcast_array_real_2d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a real 2d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use mpi  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real(4) , dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)
   
   call MPI_BCAST(array, nelements, mpi_real4, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_real_2d

 !---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_real_2d
! !INTERFACE:

 Subroutine broadcast_array_real8_2d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a real 2d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real(8) , dimension(:,:), intent(inout) :: array                ! array to be broadcast
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)
   
   call MPI_BCAST(array, nelements, MPI_REAL8, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_real8_2d
!---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_int_2d
! !INTERFACE:

 Subroutine broadcast_array_int_2d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a 2d integer array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use mpi  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer(4) , dimension(:,:), intent(inout) :: array              ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   Integer :: nelements  ! size of array to be broadcast
   Integer ::   ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_Integer4, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_int_2d

!---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_log_2d
! !INTERFACE:

 Subroutine broadcast_array_log_2d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a logical 2d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical , dimension(:,:), intent(inout) :: &
     array                ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer , dimension(:,:), allocatable :: &
     array_int            ! temporary array for MPI bcast

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)
   allocate(array_int(size(array,dim=1),size(array,dim=2)))

   where (array)
     array_int = 1
   elsewhere
     array_int = 0
   End where

   call MPI_BCAST(array_int, nelements, MPI_Int, root_pe, &
                  comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

   where (array_int == 1)
     array = .true.
   elsewhere
     array = .false.
   End where

   deallocate(array_int)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_log_2d

!---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_real_3d
! !INTERFACE:

 Subroutine broadcast_array_real_3d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a real 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use mpi  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real(4) , dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpi_real4, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_real_3d
!---------------------------------------------------------------------------!
 !---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_real_3d
! !INTERFACE:

 Subroutine broadcast_array_real8_3d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a real 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real(8) , dimension(:,:,:), intent(inout) :: array  ! array to be broadcast
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: nelements ! size of array to be broadcast
   integer  :: ierr      ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_REAL8, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_real8_3d
!---------------------------------------------------------------------------!
 
!---------------------------------------------------------------------------!

  Subroutine broadcast_array_real_4d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a real 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real(4) , dimension(:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpi_real4, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_real_4d

!---------------------------------------------------------------------------!

  Subroutine broadcast_array_real8_4d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a real 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use mpi  ! MPI Fortran include file
 Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   real(8) , dimension(:,:,:,:), intent(inout) :: &
     array                ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, mpi_real8, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_real8_4d
 
!---------------------------------------------------------------------------! 

!---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_int_3d
! !INTERFACE:

 Subroutine broadcast_array_int_3d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts an integer 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer , dimension(:,:,:), intent(inout) :: array              ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   Integer  :: nelements  ! size of array to be broadcast
   Integer  :: ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_Int, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_int_3d

!---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_int_3d
! !INTERFACE:

 Subroutine broadcast_array_int8_3d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts an integer 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!  This one handles special case where integer arrays have to be *8 precision
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use MPI  ! MPI Fortran include file
   Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   integer(8) , dimension(:,:,:), intent(inout) ::array              ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)

   call MPI_BCAST(array, nelements, MPI_Integer8, root_pe, comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

!---------------------------------------------------------------------------!
!EOC

 End Subroutine broadcast_array_int8_3d


!---------------------------------------------------------------------------!
!BOP
! !IROUTINE: broadcast_array_log_3d
! !INTERFACE:

 Subroutine broadcast_array_log_3d(array, root_pe)

! !DESCRIPTION:
!  Broadcasts a logical 3d array from one processor (root_pe)
!  to all other processors. This is a specific instance of the generic
!  broadcast\_array Interface.
!
! !REVISION HISTORY:
!  same as Module

! !INCLUDES:

   Use mpi  ! MPI Fortran include file
 Use Variables_MPI
   
! !INPUT PARAMETERS:

   integer , intent(in) :: &
     root_pe              ! processor number to broadcast from

! !INPUT/OUTPUT PARAMETERS:

   logical , dimension(:,:,:), intent(inout) :: &
     array                ! array to be broadcast

!EOP
!BOC
!---------------------------------------------------------------------------!
!
!  local variables
!
!---------------------------------------------------------------------------!

   integer , dimension(:,:,:), allocatable :: &
     array_int            ! temporary array for MPI bcast

   integer  :: &
     nelements,          &! size of array to be broadcast
     ierr                 ! local MPI error flag

!---------------------------------------------------------------------------!

   nelements = size(array)
   allocate(array_int(size(array,dim=1), &
                      size(array,dim=2), &
                      size(array,dim=3)))

   where (array)
     array_int = 1
   elsewhere
     array_int = 0
   End where

   call MPI_BCAST(array_int, nelements, MPI_Int, root_pe, &
                  comm_2d, ierr)
   call MPI_BARRIER(comm_2d, ierr)

   where (array_int == 1)
     array = .true.
   elsewhere
     array = .false.
   End where

   deallocate(array_int)

!---------------------------------------------------------------------------!

End Subroutine broadcast_array_log_3d

!---------------------------------------------------------------------------!
End Module Broadcast_Routines
#endif
