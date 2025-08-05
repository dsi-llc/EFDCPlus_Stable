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
  ! @details Generate cell mappings for the sub domain
  ! @author Zander Mausolff - adopted from O'Donncha's implementation
  ! @date 8/20/2019
  !---------------------------------------------------------------------------!

  Subroutine Child_Grid

  use GLOBAL
  use Variables_MPI

  use MPI

  implicit none

  ! *** Read in variables

  ! *** Local variables
  integer :: i, ii, j, jj, ierr

  integer :: ghost_cells_i_start
  integer :: ghost_cells_j_start
  integer :: ghost_cells_i_end
  integer :: ghost_cells_j_end
  integer :: total_ghost_i
  integer :: total_ghost_j

  ! *** Allocate arrays that contain the number of ghost cells in each partition

  !  /\             Assumes processors         |----|----|----|
  !  |     N        are aligned like this      | P6 | P7 | P8 |
  !  j   W   E      when domain is decomposed: | P3 | P4 | P5 |
  !  |     S        nx = 3 , ny = 3                | P0 | P1 | P2 |
  !  \/  <--i-->                               |----|----|----|

  ! *** Map local I to global I
  ii = 0
  Do i = IB_Decomp(process_id),IE_Decomp(process_id)
    ii = ii + 1
    IL2IG(ii) = i
  enddo

  ! *** Map local J to global J
  jj = 0
  Do j = JB_Decomp(process_id),JE_Decomp(process_id)
    jj = jj + 1
    JL2JG(jj) = j
  enddo

  ! *** Write out to MPI debug file
  ! *** Write out IL2IG, JL2JG

  if( MPI_Write_Flag )then
    ! *** Create a mapping of [I,J] starting coordinates for use in global reconciliation
    write(mpi_log_unit, '(A)') 'IL2IG: '
    do i = 1, IC
      write(mpi_log_unit, '(I5)', ADVANCE = 'NO') IL2IG(i)
    enddo
    write(mpi_log_unit, '(A)') ' '

    write(mpi_log_unit, '(A,I5)') 'JL2JG '
    do j = 1, jc
      write(mpi_log_unit, '(I5)', ADVANCE = 'NO') JL2JG(j)
    enddo
    call MPI_Barrier(DSIcomm, ierr)
  endif

  End subroutine Child_Grid