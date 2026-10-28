! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!*****************************************************************************80
!
!! Task_Division divides tasks among processors.
!
!  Discussion:
!
!    This routine assigns each of T tasks to P processors, assuming that 
!    the assignment is to be beforehand.
!
!    In that case, we just want to make sure that we assign each task
!    to a processor, that we assign about the same number of tasks
!    to each processor, and that we assign each processor a contiguous
!    range of tasks, say tasks I_LO to I_HI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TASK_NUMBER, the number of tasks.
!
!    Input, integer PROC_FIRST, PROC_LAST, the first and last processors.
!
#ifdef ENABLE_MPI
subroutine Task_Division ( task_number, proc_first, proc_last )
 
    Use GLOBAL
    Use Variables_MPI
    Use Broadcast_Routines  
    Use MPI
    
    implicit none

    integer  i_hi
    integer  i_lo
    integer  i4_div_rounded
    integer  cur_p
    integer  proc
    integer  proc_first
    integer  proc_last

    integer  proc_remain
    integer  task_number
    integer  task_proc
    integer  task_remain

    integer :: total_num_constituents
    integer :: i
    Integer :: writeoutid
    Character(len=3) :: outfilename
    
    cur_p = proc_last + 1 - proc_first
    
    if(process_id == master_id) then 
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Task_Division'
        write ( *, '(a)' ) '  Divide T tasks among P processors.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Number of tasks T = ', task_number
        write ( *, '(a,i8)' ) '  Number of processors P = ', cur_p
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  P_FIRST = ', proc_first
        write ( *, '(a,i8)' ) '  P_LAST =  ', proc_last
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '             Number of   First      Last'
        write ( *, '(a)' ) ' Processor     Tasks     Task       Task'
        write ( *, '(a)' ) ' '
    End if
     
    writeoutid = process_id
    Write(outfilename, '(I3)') process_id
    Open(unit_num_MPI+writeoutid, File = 'MPI_info_proc_'//trim(outfilename)//'.txt', Status='Unknown')
    Close(unit_num_MPI+writeoutid, Status='Delete')
    Open(unit_num_MPI+writeoutid, File = 'MPI_info_proc_'//trim(outfilename)//'.txt', Position='Append')
     
    write (unit_num_MPI+writeoutid , '(a)' ) 'Task_Division'
    write (unit_num_MPI+writeoutid , '(a)' ) '  Divide T tasks among P processors.'
    write (unit_num_MPI+writeoutid , '(a)' ) ' '
    write (unit_num_MPI+writeoutid , '(a,i8)' ) '  Number of tasks T = ', task_number
    write (unit_num_MPI+writeoutid , '(a,i8)' ) '  Number of processors P = ', cur_p
    write (unit_num_MPI+writeoutid , '(a)' ) ' '
    write (unit_num_MPI+writeoutid , '(a,i8)' ) '  P_FIRST = ', proc_first
    write (unit_num_MPI+writeoutid , '(a,i8)' ) '  P_LAST =  ', proc_last
    write (unit_num_MPI+writeoutid , '(a)' ) ' '
    write (unit_num_MPI+writeoutid , '(a)' ) '             Number of   First      Last'
    write (unit_num_MPI+writeoutid , '(a)' ) ' Processor     Tasks     Task       Task'
    write (unit_num_MPI+writeoutid , '(a)' ) ' '
    
    i_hi = 0

    task_remain = task_number
    proc_remain = cur_p
    
    do proc = proc_first, proc_last
        task_proc = i4_div_rounded ( task_remain, proc_remain )

        proc_remain = proc_remain - 1
        task_remain = task_remain - task_proc
        i_lo = i_hi + 1
        i_hi = i_hi + task_proc
        ! Arrays to keep track of lower and upper bound for contituent transport per process
        lower_task_bound(proc+1)    = i_lo 
        upper_task_bound(proc+1)    = i_hi
       
        if(process_id == master_id) then
            write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) proc, task_proc, i_lo, i_hi
        end if
        num_task_per_process(proc+1) = task_proc 
    
       write (unit_num_MPI+writeoutid, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) proc, task_proc, i_lo, i_hi
    end do
    
    !***Global Variable
    total_num_constituents = upper_task_bound(num_proc) - lower_task_bound(0)
     
    !***Setup send counts and displacements_L_index array 
    Do i = 1, num_proc 
        num_elem_to_each_array(i) = LCM*KCM*num_task_per_process(i)
        If(i == 1) then ! Starting displacement should be zero
            displacements_L_index(i) = 0 
        else
            displacements_L_index(i) = displacements_L_index(i-1) + num_task_per_process(i-1) 
        End if
    End do
 
    write(unit_num_MPI+writeoutid, '(a)') 'displacements_L_index:'
    write(unit_num_MPI+writeoutid, '(i8)' ) (displacements_L_index(i), i=1, num_proc) 
 
    Do i = 1, num_proc
          recv_counts_array(i) = LCM*KCM*num_task_per_process(i)
    End do 
 
    write(unit_num_MPI+writeoutid, '(a)') 'Recieve counts array:'
    write(unit_num_MPI+writeoutid, '(i8)' ) (recv_counts_array(i), i=1, num_proc) 

    return

end subroutine Task_Division

function i4_div_rounded ( a, b )

!*  ****************************************************************************80
!
!!   I4_DIV_ROUNDED computes the rounded result of I4 division.
!
!    Discussion:
!
!      This routine computes C = A / B, where A, B and C are integers
!      and C is the closest integer value to the exact real result.
!
!    Licensing:
!
!      This code is distributed under the GNU LGPL license.
!
!    Modified:
!
!      23 October 2011
!
!    Author:
!
!      John Burkardt
!
!    Parameters:
!
!      Input, integer ( kind = 4 ) A, the number to be divided.
!
!      Input, integer ( kind = 4 ) B, the divisor, or the number of parts.
!
!      Output, integer ( kind = 4 ) I4_DIV_ROUNDED, the rounded result
!      of the division.
!
    implicit none

    integer  a
    integer  a_abs
    integer  b
    integer  b_abs
    integer  i4_div_rounded
    integer , parameter :: i4_huge = 2147483647
    integer  value

    if ( a == 0 .and. b == 0 ) then

      value = i4_huge
 
    else if ( a == 0 ) then

      value = 0

    else if ( b == 0 ) then

      if ( a < 0 ) then
        value = - i4_huge
      else
        value = + i4_huge
      end if

    else

      a_abs = abs ( a )
      b_abs = abs ( b )

      value = a_abs / b_abs
!
!    Round the value.
!
      if ( ( 2 * value + 1 ) * b_abs < 2 * a_abs ) then
        value = value + 1
      end if
!
!    Set the sign.
!
      if ( ( a < 0 .and. 0 < b ) .or. ( 0 < a .and. b < 0 ) ) then
        value = - value
      end if

    end if

    i4_div_rounded = value

    return
end
#endif
