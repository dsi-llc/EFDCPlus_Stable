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
  ! @details Lets each process get array that helps map
  !! from global 1 to max_i to 1 to max_local_i
  ! @author Zander Mausolff - adopted from O'Donncha's code.
  ! @date 8/13/2019
  !---------------------------------------------------------------------------!

  Subroutine Parent_Grid

  use GLOBAL
  use variables_mpi
  use Mod_Write_Cell_Map, only: Write_Xloc_Yloc

  implicit none

  ! *** Read in variables

  ! *** Local variables
  integer :: ii, jj, i ,j, ierr

  ic_pos = 0 
  jc_pos = 0

  ! *** Gets starting point for process in 'i' space.  ic_decomp and jc_decomp are excluding ghost cells
  Do ii = 2, x_id
    ic_pos = ic_pos + ic_decomp(ii-1)         
  enddo

  ! *** Gets starting point for process in 'j' space
  Do jj = 2, y_id
    jc_pos = jc_pos + jc_decomp(jj-1)
  enddo

  ! *** Map global I to local I
  ii = 0
  Do i = IB_Decomp(process_id),IE_Decomp(process_id)
    ii = ii + 1
    IG2IL(i) = ii
  enddo

  ! *** Map global J to local J
  jj = 0
  Do j = JB_Decomp(process_id),JE_Decomp(process_id)
    jj = jj + 1
    JG2JL(j) = jj
  enddo

  if( MPI_Write_Flag )then
    call WriteBreak(mpi_log_unit)
    call Write_Xloc_Yloc
    write(mpi_log_unit, '(a,I4)') 'X_id    = ', x_id
    write(mpi_log_unit, '(a,I4)') 'Y_id    = ', y_id
    write(mpi_log_unit, '(a,I4)') 'IC pos  = ', ic_pos
    write(mpi_log_unit, '(a,I4)') 'JC pos  = ', jc_pos
    write(mpi_log_unit, '(a,I4)') 'IC      = ', ic
    write(mpi_log_unit, '(a,I4)') 'JC      = ', jc
    call WriteBreak(mpi_log_unit)
  endif

  End subroutine Parent_Grid
