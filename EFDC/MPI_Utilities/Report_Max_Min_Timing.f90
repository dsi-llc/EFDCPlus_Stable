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
!> @details Collects timing informations from all processes and reports the
!! min,max, and ratio max/min.  Writes all to the TIME.LOG file
!> @author Zander Mausolff
!> @date 2/11/2020
!---------------------------------------------------------------------------!

Subroutine Report_Max_Min_Timing(routine_time, time_id)

  Use GLOBAL
  Use Variables_MPI
  Use MPI

  Implicit None

  ! *** Passed in Variables
  Real(RKD), intent(in)    :: routine_time
  CHARACTER*8, intent(in)  :: time_id

  ! *** Local Variables
  Real(RKD), Allocatable, Dimension(:) :: timing_frm_all_procs
  Real(RKD) :: min_time !< Minimum time spent in a given routine
  Real(RKD) :: max_time !< Maximum time spent in a given routine
  Real(RKD) :: ratio_max_min !< Ratio of the max to the min to give idea to load balancing issue
  Integer   :: process_min
  Integer   :: process_max
  Integer   :: time_log_unit, i, ierr
  Logical, Save :: first_time = .TRUE.

  time_log_unit = 9 ! this is the unit for the TIME.LOG file
  ! ***  Write to the TIME.LOG file only on master

  if(process_id == master_id )then
    OPEN(time_log_unit,FILE=OUTDIR//'TIME.LOG',POSITION='APPEND')
  endif

  if( first_time == .TRUE.  )then
    if(process_id == master_id )then
      Call WriteBreak(time_log_unit)
      WRITE(time_log_unit,'(A,/)') 'Reporting Max/Min from each process for load balancing analysis'
    endif

    first_time = .FALSE.

  endif

  Allocate(timing_frm_all_procs(num_Processors))
  timing_frm_all_procs = 0.0
  
  ! *** Get the time each process has spent in a respective routine
  ! *** This is collected in an array so we can identify which process was the min/max
  Call MPI_AllGather(routine_time, 1, MPI_Real8, timing_frm_all_procs, 1, MPI_Real8, comm_2d, ierr)
  
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
      write(time_log_unit, '(a,a,a,f6.3,a,i5) ') 'Time for ', time_id, ' = ',timing_frm_all_procs(i), ' in process ', i
    end do

    write(time_log_unit, '(a,a,a,f6.3,a,i5)') 'Min value for ',time_id,'=', min_time, ' in process ', process_min
    write(time_log_unit, '(a,a,a,f6.3,a,i5)') 'Max value for ',time_id,'=', max_time, ' in process ', process_max
    write(time_log_unit, '(a,f6.3)')'Ratio Max/Min = ', ratio_max_min
    write(time_log_unit, '(a,f6.3,a)')'Percent of max time to elapsed time = ',  (max_time/TIME_END)*100, '%'
    Call WriteBreak(time_log_unit)

  end if

  close(time_log_unit)

  Deallocate(timing_frm_all_procs)

End subroutine Report_Max_Min_Timing
