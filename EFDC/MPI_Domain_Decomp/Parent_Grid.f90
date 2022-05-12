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
  ! @details Lets each process get array that helps map
  !! from global 1 to max_i to 1 to max_local_i
  ! @author Zander Mausolff - adopted from O'Donncha's code.
  ! @date 8/13/2019
  !---------------------------------------------------------------------------!

  Subroutine Parent_Grid

  Use global
  Use variables_mpi
  Use Mod_Write_Cell_Map, only: Write_Xloc_Yloc

  Implicit none

  ! *** Read in variables

  ! *** Local variables
  Integer :: ii, jj, i ,j, ierr

  ic_pos = 0 
  jc_pos = 0

  ! *** Gets starting point for process in 'i' space.  ic_decomp and jc_decomp are excluding ghost cells
  Do ii = 2, x_id
    ic_pos = ic_pos + ic_decomp(ii-1)         
  End do

  ! *** Gets starting point for process in 'j' space
  Do jj = 2, y_id
    jc_pos = jc_pos + jc_decomp(jj-1)
  End do

  ! *** Map global I to local I
  ii = 0
  Do i = IB_Decomp(process_id),IE_Decomp(process_id)
    ii = ii + 1
    IG2IL(i) = ii
  End do

  ! *** Map global J to local J
  jj = 0
  Do j = JB_Decomp(process_id),JE_Decomp(process_id)
    jj = jj + 1
    JG2JL(j) = jj
  End do

  if( MPI_Write_Flag == .TRUE. )then
    Call WriteBreak(mpi_log_unit)
    Call Write_Xloc_Yloc
    write(mpi_log_unit, '(a,I4)') 'X_id   =', x_id
    write(mpi_log_unit, '(a,I4)') 'Y_id   =', y_id
    write(mpi_log_unit, '(a,I4)') 'IC pos =', ic_pos
    write(mpi_log_unit, '(a,I4)') 'JC pos =', jc_pos
    write(mpi_log_unit, '(a,I4)') 'IC     =', ic
    write(mpi_log_unit, '(a,I4)') 'JC     =', jc
    Call WriteBreak(mpi_log_unit)
  endif

  End subroutine Parent_Grid
