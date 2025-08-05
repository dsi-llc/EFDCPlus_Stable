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
!> @details Collects timing informations from all processes and reports the
!! min,max, and ratio max/min.  Writes all to the TIME.LOG file
!> @author Zander Mausolff
!> @date 2/11/2020
!---------------------------------------------------------------------------!

Subroutine Report_Max_Min_Timing(routine_time, time_id)

  use GLOBAL
  use Variables_MPI
  use MPI

  implicit none

  ! *** Passed in Variables
  real(RKD), intent(in)    :: routine_time
  character*8, intent(in)  :: time_id

  ! *** Local Variables
  real(RKD), Allocatable, Dimension(:) :: timing_frm_all_procs
  real(RKD) :: min_time !< Minimum time spent in a given routine
  real(RKD) :: max_time !< Maximum time spent in a given routine
  real(RKD) :: ratio_max_min !< Ratio of the max to the min to give idea to load balancing issue
  integer   :: process_min
  integer   :: process_max
  integer   :: time_log_unit, i, ierr
  logical, Save :: first_time = .TRUE.

  time_log_unit = 9 ! this is the unit for the TIME.LOG file
  ! ***  Write to the TIME.LOG file only on master

  if(process_id == master_id )then
    open(time_log_unit,FILE = OUTDIR//'TIME.LOG',POSITION = 'APPEND')
  endif

  if( first_time )then
    if(process_id == master_id )then
      call WriteBreak(time_log_unit)
      write(time_log_unit,'(A,/)') 'Reporting Max/Min from each process for load balancing analysis'
    endif

    first_time = .FALSE.

  endif

  allocate(timing_frm_all_procs(num_Processors))
  timing_frm_all_procs = 0.0
  
  ! *** Get the time each process has spent in a respective routine
  ! *** This is collected in an array so we can identify which process was the min/max
  call MPI_AllGather(routine_time, 1, MPI_Real8, timing_frm_all_procs, 1, MPI_Real8, comm_2d, ierr)
  
  ! ***  Write to the TIME.LOG file only on master
  if(process_id == master_id )then
    ! *** Get min/max
    min_time = minval(timing_frm_all_procs, 1)
    max_time = maxval(timing_frm_all_procs, 1)

    ! *** Find out which process had that min/max time
    process_min = findloc(timing_frm_all_procs, min_time, 1)
    process_max = findloc(timing_frm_all_procs, max_time, 1)


    ! *** Get the ratio between the max/min
    if(max_time > 1E-8 .and. min_time > 1E-8 )then
      ratio_max_min = max_time/min_time
    else
      ratio_max_min = 0.0
    endif

    do i = 1, num_Processors
      write(time_log_unit, '(a,a,a,f8.3,a,i5) ') 'Time for ', time_id, ' = ',timing_frm_all_procs(i), ' in process ', i
    enddo

    write(time_log_unit, '(a,a,a,f8.3,a,i5)') 'Min value for ',time_id,' = ', min_time, ' in process ', process_min
    write(time_log_unit, '(a,a,a,f8.3,a,i5)') 'Max value for ',time_id,' = ', max_time, ' in process ', process_max
    write(time_log_unit, '(a,f8.3)')'Ratio Max/Min = ', ratio_max_min
    write(time_log_unit, '(a,f8.3,a)')'Percent of max time to elapsed time = ', (max_time/TIME_END)*100., '%'
    call WriteBreak(time_log_unit)

  endif

  close(time_log_unit)

  deallocate(timing_frm_all_procs)

End subroutine Report_Max_Min_Timing
