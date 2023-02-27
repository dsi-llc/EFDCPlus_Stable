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
  ! @details Generate cell mappings for the sub domain
  ! @author Zander Mausolff - adopted from O'Donncha's implementation
  ! @date 8/20/2019
  !---------------------------------------------------------------------------!

  Subroutine Child_Grid

  Use Global
  Use Variables_MPI

#ifdef _MPI
  Use MPI
#endif
  Implicit None

  ! *** Read in variables

  ! *** Local variables
  Integer :: i, ii, j, jj, ierr

  Integer :: ghost_cells_i_start
  Integer :: ghost_cells_j_start
  Integer :: ghost_cells_i_end
  Integer :: ghost_cells_j_end
  Integer :: total_ghost_i
  Integer :: total_ghost_j

  ! *** Allocate arrays that contain the number of ghost cells in each partition

  !  /\             Assumes processors         |----|----|----|
  !  |     N        are aligned like this      | P6 | P7 | P8 |
  !  j   W   E      when domain is decomposed: | P3 | P4 | P5 |
  !  |     S        nx=3 , ny=3                | P0 | P1 | P2 |
  !  \/  <--i-->                               |----|----|----|

  ! *** Map local I to global I
  ii = 0
  Do i = IB_Decomp(process_id),IE_Decomp(process_id)
    ii = ii + 1
    IL2IG(ii) = i
  End do

  ! *** Map local J to global J
  jj = 0
  Do j = JB_Decomp(process_id),JE_Decomp(process_id)
    jj = jj + 1
    JL2JG(jj) = j
  End do

  ! *** Write out to MPI debug file
  ! *** Write out IL2IG, JL2JG

  if( MPI_Write_Flag == .TRUE. )then
    ! *** Create a mapping of [I,J] starting coordinates for use in global reconciliation
    write(mpi_log_unit, '(A)') 'IL2IG: '
    do i = 1, IC
      write(mpi_log_unit, '(I5)', ADVANCE='NO') IL2IG(i)
    end do
    write(mpi_log_unit, '(A)') ' '

    write(mpi_log_unit, '(A,I5)') 'JL2JG '
    do j = 1, jc
      write(mpi_log_unit, '(I5)', ADVANCE='NO') JL2JG(j)
    end do
    Call MPI_Barrier(comm_2d, ierr)
  End if

  End subroutine Child_Grid