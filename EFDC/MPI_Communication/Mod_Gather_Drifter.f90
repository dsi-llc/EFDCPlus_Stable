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
! Module: Mod_Gather_Drifter
!
!> @details contains routines for gatherin up all drifter information
!
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 2/18/2020
!---------------------------------------------------------------------------!
    
Module Mod_Gather_Drifter
    
    use GLOBAL
    use Variables_MPI
    use Variables_MPI_Drifter
    use MPI
    
    implicit none
    
    save
    
    Public Gather_Drifter
    
    interface Gather_Drifter
        
        Module Procedure Gather_Int_Drifter, &
                         Gather_Real_Drifter_RK8
                         
    end interface
    
    Contains
!---------------------------------------------------------------------------!  
!                     EFDC+ Developed by DSI, LLC. 
!---------------------------------------------------------------------------!  
! Module: Mod_Gather_Drifter
!
!> @details contains routines for gatherin up all drifter information
!
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 2/18/2020
!---------------------------------------------------------------------------!
Subroutine Gather_Int_Drifter(size_local_1d, Soln_Local_1D, size_global_1d, Soln_Global_1D)
    
    implicit none
        
    !***Passed in
    integer, Intent(in)    :: size_local_1d 
    integer, Intent(in)    :: Soln_Local_1D(size_local_1d)
    integer, Intent(in)    :: size_global_1d
    integer, Intent(inout) :: Soln_Global_1D(size_global_1d)
    
    !***Local variables
    integer :: ierr
    integer :: i
    integer, Allocatable :: displacements_index(:)
    integer, Allocatable :: recv_counts(:)
    
    If(.not.allocated(recv_counts) )then
        allocate(recv_counts(active_domains))
        allocate(displacements_index(active_domains))
    endif
    ! *** We only want to send the number of active cells in the 1D array
    
    send_size_1d = size_local_1d
    
    recv_counts(:) = 0
    
    !***make sure each process has the same recv counts array for the 1d case
    call MPI_Allgather(send_size_1d,   1, MPI_Int, recv_counts, 1, MPI_Int, DSIcomm, ierr)
    
    ! *** setup the array that keeps track of where data from all processes will be places
    displacements_index(1) = 0 ! 
    
    do i = 2, active_domains
        displacements_index(i) = displacements_index(i-1) + recv_counts(i-1)
    enddo 
    
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_1d, MPI_Integer, &                         ! *** Send buff
                     Soln_Global_1D, recv_counts,  displacements_index, MPI_Integer,   &  ! *** Recv buff
                     master_id, DSIcomm, ierr)

    if(MPI_DEBUG_FLAG )then
        write(mpi_log_unit,'(a,20i5)') 'active_domains: ',active_domains
        write(mpi_log_unit,'(a,20i5)') 'send_size_1d: ', send_size_1d
        write(mpi_log_unit,'(a,20i5)') 'recv_counts: ', recv_counts
        write(mpi_log_unit,'(a,20i5)') 'displacements_index: ', displacements_index
    endif

                     
End Subroutine Gather_Int_Drifter

!---------------------------------------------------------------------------!  
!                     EFDC+ Developed by DSI, LLC. 
!---------------------------------------------------------------------------!  
! Module: Mod_Gather_Drifter
!
!> @details contains routines for gatherin up all drifter information
!
!---------------------------------------------------------------------------!
!> @author Zander Mausolff
!> @date 2/18/2020
!---------------------------------------------------------------------------!
Subroutine Gather_Real_Drifter_RK8(size_local_1d,  Soln_Local_1D, size_global_1d, Soln_Global_1D)

    implicit none
    
    !***Read in
    integer,   Intent(in)    :: size_local_1d 
    real(RKD), Intent(in)    :: Soln_Local_1D(size_local_1d)
    integer,   Intent(in)    :: size_global_1d
    real(RKD), Intent(inout) :: Soln_Global_1D(size_global_1d) 
    
    !***Local variables
    integer :: ierr
    integer :: i
    integer, Allocatable :: displacements_index(:)
    integer, Allocatable :: recv_counts(:)
    
    If(.not.allocated(recv_counts) )then
        allocate(recv_counts(active_domains))
        allocate(displacements_index(active_domains))
    endif
    
    !***we only want to send the number of drifters in
    send_size_1d = size_local_1d
    
    recv_counts(:) = 0
    displacements_index(:) = 0
    
    !***make sure each process has the same recv counts array for the 1d case
    call MPI_Allgather(send_size_1d, 1, MPI_Int, recv_counts, 1, MPI_Int, DSIcomm, ierr)

    ! *** setup the array that keeps track of where data from all processes will be places
    displacements_index(1) = 0 ! 
    
    do i = 2, active_domains
        displacements_index(i) = displacements_index(i-1) + recv_counts(i-1)
    enddo 
    
    if( MPI_DEBUG_FLAG )then
        write(mpi_log_unit,'(a,3010i5)') 'displacements_index: ',displacements_index
    endif
    
    !***Gather solution onto global 1D array
    call MPI_Gatherv(Soln_Local_1D,  send_size_1d, mpi_real8, &                         ! *** Send buff
                     Soln_Global_1D, recv_counts,  displacements_index, mpi_real8,   &  ! *** Recv buff
                     master_id, DSIcomm, ierr)
    
End Subroutine Gather_Real_Drifter_RK8
    
End Module Mod_Gather_Drifter