! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------

Subroutine WriteBreak(unit_num)

    Implicit none
    
    ! *** Passed in variables
    Integer, intent(in) :: unit_num
    ! *** Local Variables
    write(unit_num, '(a)') ' '
    write (unit_num, '(a)' ) '|------------------------------------------------------------------------------|'  
    write(unit_num, '(a)') ' '
    
End subroutine WriteBreak
    
    
Subroutine WriteInteger(int_to_write_out, unit_num,variable_desc)

    Implicit None
    
    ! *** Passed in variables
    Integer, intent(in) :: int_to_write_out
    Integer, intent(in) :: unit_num
    Character(20), intent(in) :: variable_desc
    
    
    write(unit_num, '(A20,I5)') variable_desc, int_to_write_out

End subroutine WriteInteger 
