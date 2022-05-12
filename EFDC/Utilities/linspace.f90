! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!< @details Generates evenly spaced numbers from `lower` to `upper` (inclusive)
!< @author Luis Bastidas 
!< @param[in] arr_size
!< @param[in] lower
!< @param[in] upper
!< @param[out] evenly_space
!---------------------------------------------------------------------------!
subroutine linspace(arr_size, lower, upper, evenly_space) 

    Use Global, only : RKD
    
    implicit none

    ! *** Dummy variables
    integer, intent(in) :: arr_size            !< size of array
    real(kind = RKD), intent(in) :: lower !< lower bound of the array
    real(kind = RKD), intent(in) :: upper !< upper bound of the array
    real(kind = RKD), intent(out) :: evenly_space(arr_size) !< returned array with evenly spaced values

    ! *** Local variables
    real(kind = RKD) :: range
    integer :: i

    range = upper - lower

    ! *** 
    IF( arr_size == 0 )then 
        return
    end if
    ! *** Handle array of size 1
    IF( arr_size == 1 )then
        evenly_space(1) = lower
        return
    end if

    ! *** generate the array of evenly spaced values
    do i = 1, arr_size
        evenly_space(i) = lower + range * (i - 1) / (arr_size - 1)
    end do

    return

end subroutine linspace