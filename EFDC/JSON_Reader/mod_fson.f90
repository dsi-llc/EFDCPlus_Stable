! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit 
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or 
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.


! FSON MODULE 
!
! File:   fson.f95
! Author: Joseph A. Levin
!
! Created on March 6, 2012, 7:48 PM
!

module fson
    use mod_fson_value, fson_print => fson_value_print, fson_destroy => fson_value_destroy
    use mod_fson_string
    use mod_fson_path, fson_get => fson_path_get

    implicit none

    private

    public :: fson_parse, fson_value, fson_get, fson_print, fson_destroy

    ! FILE IOSTAT CODES
    integer, parameter :: end_of_file = -1
    integer, parameter :: end_of_record = -2

    ! PARSING STATES
    integer, parameter :: STATE_LOOKING_FOR_VALUE = 1
    integer, parameter :: STATE_IN_OBJECT = 2
    integer, parameter :: STATE_IN_PAIR_NAME = 3
    integer, parameter :: STATE_IN_PAIR_VALUE = 4

    ! POP/PUSH CHARACTER
    integer :: pushed_index = 0
    character (len = 10) :: pushed_char

contains

    !
    ! FSON PARSE
    !
    function fson_parse(file, unit, str) result(p)
        type(fson_value), pointer :: p
        integer, optional, intent(inout) :: unit
        character(len = *), optional, intent(in) :: file
        character(len = *),target, optional, intent(in) :: str
        character(len = :),target, allocatable :: strBuffer

        logical :: unit_available
        integer :: u
        
        ! init the pointer to null
        nullify(p)

        ! select the file unit to use
        IF( present(unit) .and. present(file) )then
            u = unit
        elseif (present(file) )then
            ! find the first available unit
            unit_available = .false.
            u = 20

            do while (.not.unit_available)
                inquire(unit = u, exist = unit_available)
                u = u + 1
            end do
        elseif (present(str) )then
            strBuffer = str
            u = 0
        else 
            print *, "ERROR: Need a file or a string"
            call Stopp('.')    ! exit(1)
        end if

        ! open the file
        IF( present(file) )then
            open (unit = u, file = file, status = "old", action = "read", form = "formatted", position = "rewind")
        end if
        
        ! 2021-06-14, NTL: fixed unallocated strBuffer
        if (.not.present(str) )then
            allocate(character(len=1024) :: strBuffer)
        end if
        
        ! create the value and associate the pointer        
        p => fson_value_create()

        ! parse as a value
        call parse_value(unit=u, str=strBuffer, value=p)

        ! close the file
        if( .not. present(unit) )then
            close (u)
        end if

        !if(allocated(strBuffer)) deallocate(strBuffer)

    end function fson_parse

    !
    ! PARSE_VALUE
    !
    recursive subroutine parse_value(unit, str, value)
        integer, intent(inout) :: unit
        character(len = *),target,intent(inout) :: str
        type(fson_value),intent(inout), pointer :: value
        
        logical :: eof
        character :: c

        ! pop the next non whitespace character off the file
        c = pop_char(unit, str, eof = eof, skip_ws = .true.)

        IF( eof )then
            return
        else
            select case (c)
            case ("{")
                ! start object                
                value.value_type = TYPE_OBJECT
                call parse_object(unit, str, value)
            case ("[")
                ! start array
                value.value_type = TYPE_ARRAY
                call parse_array(unit, str, value)
            case ("]")
                ! end an empty array
               call push_char(c)
               nullify(value)
            case ('"')
                ! string                                      
                value.value_type = TYPE_STRING
                value.value_string => parse_string(unit, str)
            case ("t")
                !true
                value.value_type = TYPE_LOGICAL
                call parse_for_chars(unit, str, "rue")
                value.value_logical = .true.
            case ("f")
                !false
                value.value_type = TYPE_LOGICAL
                value.value_logical = .false.
                call parse_for_chars(unit, str, "alse")
            case ("n")
                value.value_type = TYPE_NULL
                call parse_for_chars(unit, str, "ull")
            case("-", "0": "9")
                call push_char(c)
                call parse_number(unit, str, value)
            case default
                print *, "ERROR: Unexpected character while parsing value. '", c, "' ASCII=", iachar(c)
                call Stopp('.')    ! exit(1)
            end select
        end if

    end subroutine parse_value

    !
    ! PARSE OBJECT
    !    
    recursive subroutine parse_object(unit, str, parent)
        integer, intent(inout) :: unit
        character(len=*),target, intent(inout) :: str
        type(fson_value), pointer :: parent, pair


        logical :: eof
        character :: c

        ! pair name
        c = pop_char(unit, str, eof = eof, skip_ws = .true.)
        IF( eof )then
            print *, "ERROR: Unexpected end of file while parsing start of object."
            call Stopp('.')    ! exit(1)
        else IF( "}" == c )then
            ! end of an empty object
            return
        else IF( '"' == c )then
            pair => fson_value_create()
            pair.name => parse_string(unit, str)
        else
            print *, "ERROR: Expecting string: '", c, "'"
            call Stopp('.')    ! exit(1)
        end if

        ! pair value
        c = pop_char(unit, str, eof = eof, skip_ws = .true.)
        IF( eof )then
            print *, "ERROR: Unexpected end of file while parsing object member. 1"
            call Stopp('.')    ! exit(1)
        else IF( ":" == c )then
            ! parse the value                       
            call parse_value(unit, str, pair)
            call fson_value_add(parent, pair)
        else
            print *, "ERROR: Expecting : and then a value. ", c
            call Stopp('.')    ! exit(1)
        end if

        ! another possible pair
        c = pop_char(unit, str, eof = eof, skip_ws = .true.)
        IF( eof )then
            return
        else IF( "," == c )then
            ! read the next member            
            call parse_object(unit = unit, str=str, parent = parent)
        else IF( "}" == c )then
            return
        else
            print *, "ERROR: Expecting end of object.", c
            call Stopp('.')    ! exit(1)
        end if

    end subroutine parse_object

    !
    ! PARSE ARRAY
    !    
    recursive subroutine parse_array(unit, str, array)

      implicit none
      integer, intent(inout) :: unit
      character(len=*),target,  intent(inout) :: str
      type(fson_value), pointer :: array, element

      logical :: eof, finished
      character :: c

      finished = .false.
      do while (.not. finished)

         ! try to parse an element value
         element => fson_value_create()
         call parse_value(unit, str, element)

         ! parse value will disassociate an empty array value
         IF( associated(element) )then
            call fson_value_add(array, element)
         end if

         ! pop the next character
         c = pop_char(unit, str, eof = eof, skip_ws = .true.)

         IF( eof )then
            finished = .true.
         else IF( "]" == c )then
            ! end of array
            finished = .true.
         end if

      end do

    end subroutine parse_array

    !
    ! PARSE STRING
    !
    function parse_string(unit, str) result(string)
        integer, intent(inout) :: unit
        character(*), intent(inout) :: str
        type(fson_string), pointer :: string

        logical :: eof, escape
        character :: c

        string => fson_string_create()
        escape = .false.

        do
            c = pop_char(unit, str, eof = eof, skip_ws = .false.)
            IF( eof )then
               print *, "Expecting end of string"
               call Stopp('.')    ! exit(1)
            else IF( escape )then
              call fson_string_append(string,c)
              escape = .false.
            else
               IF( c == '\' )then
                  escape = .true.
               else IF( c == '"' )then
                  exit
               else
                  call fson_string_append(string,c)
               end if
            end if
        end do
    end function parse_string

    !
    ! PARSE FOR CHARACTERS
    !
    subroutine parse_for_chars(unit, str, chars)
        integer, intent(in) :: unit
        character(*), intent(inout) :: str
        character(len = *), intent(in) :: chars
        integer :: i, length
        logical :: eof
        character :: c

        length = len_trim(chars)

        do i = 1, length
            c = pop_char(unit, str, eof = eof, skip_ws = .true.)
            IF( eof )then
                print *, "ERROR: Unexpected end of file while parsing array."
                call Stopp('.')    ! exit(1)
            else IF( c .ne. chars(i:i) )then
                print *, "ERROR: Unexpected character.'", c,"'", chars(i:i)
                call Stopp('.')    ! exit(1)
            end if
        end do

    end subroutine parse_for_chars

    !
    ! PARSE NUMBER
    !
    subroutine parse_number(unit, str, value)
        integer, intent(inout) :: unit
        character(*), intent(inout) :: str
        type(fson_value), pointer :: value
        logical :: eof, negative, decimal, scientific
        character :: c
        integer :: exp, digit_count
        integer(kind=8) :: integral
        double precision :: frac


        ! first character is either - or a digit        
        c = pop_char(unit, str, eof = eof, skip_ws = .true.)
        IF( eof )then
            print *, "ERROR: Unexpected end of file while parsing number."
            call Stopp('.')    ! exit(1)
        else IF( "-" == c )then
            negative = .true.
        else
            negative = .false.
            call push_char(c)
        end if

        ! parse the integral
        integral = parse_integer(unit, str)

        decimal = .false.
        scientific = .false.
        exp = 0
        frac = 0.0d0

        do
            ! first character is either - or a digit        
            c = pop_char(unit, str, eof = eof, skip_ws = .true.)
            IF( eof )then
                print *, "ERROR: Unexpected end of file while parsing number."
                call Stopp('.')    ! exit(1)
            else
                select case (c)
                case (".")
                    ! this is already fractional number
                    IF( decimal )then
                        ! already found a decimal place
                        print *, "ERROR: Unexpected second decimal place while parsing number."
                        call Stopp('.')    ! exit(1)
                    end if
                    decimal = .true.
                    frac = parse_integer(unit, str, digit_count, allow_truncate = .true.)
                    frac = frac / (10.0d0 ** digit_count)
                case ("e", "E")
                    ! this is already an exponent number
                    IF( scientific )then
                        ! already found a e place
                        print *, "ERROR: Unexpected second exponent while parsing number."
                        call Stopp('.')    ! exit(1)
                    end if
                    scientific = .true.
                    decimal = .true.
                    ! this number has an exponent
                    exp = int(parse_integer(unit, str), kind = 4)
                case default
                    IF( decimal )then
                        ! add the integral
                        frac = frac + integral
                        IF( scientific )then
                            ! apply exponent
                            frac = frac * (10.0d0 ** exp)
                        end if
                        ! apply negative
                        IF( negative )then
                            frac = -frac
                        end if
                        value.value_type = TYPE_REAL
                        value.value_real = real(frac)
                        value.value_double = frac
                    else
                       IF( negative )then
                          ! apply negative
                          integral = -integral
                       end if
                       value.value_type = TYPE_INTEGER
                       value.value_integer = int(integral, kind = 4)
                       value.value_long_integer = integral
                    end if
                    call push_char(c)
                    exit
                end select
            end if
        end do

    end subroutine

    !
    ! PARSE INTEGER    
    !
    integer(kind=8) function parse_integer(unit, str, digit_count, allow_truncate) &
         result(integral)
        integer, intent(in) :: unit
        character(*), intent(inout) :: str
        integer, optional, intent(out) :: digit_count
        logical, optional, intent(in) :: allow_truncate
        logical :: eof, found_sign, found_digit
        character :: c
        integer :: tmp, icount, isign
        logical :: do_truncate, truncating
        integer, parameter :: max_integer_length = 18

        IF( present(allow_truncate) )then
           do_truncate = allow_truncate
        else
           do_truncate = .false.
        end if

        icount = 0
        integral = 0
        isign = 1
        found_sign = .false.
        found_digit = .false.
        truncating = .false.
        do
            c = pop_char(unit, str, eof = eof, skip_ws = .true.)
            IF( eof )then
                print *, "ERROR: Unexpected end of file while parsing digit."
                call Stopp('.')    ! exit(1)
            else
                select case(c)
                case ("+")
                    IF( found_sign.or.found_digit )then
                        print *, "ERROR: Misformatted number."
                        call Stopp('.')    ! exit(1)
                    end if
                    found_sign = .true.
                case ("-")
                    IF( found_sign.or.found_digit )then
                        print *, "ERROR: Misformatted number."
                        call Stopp('.')    ! exit(1)
                    end if
                    found_sign = .true.
                    isign = -1
                case ("0":"9")
                    found_sign = .true.
                    IF( (icount > max_integer_length) .and. (.not. truncating) )then
                       IF( do_truncate )then
                          truncating = .true.
                       else
                          print *, "ERROR: Too many digits for an integer."
                          call Stopp('.')    ! exit(1)
                       end if
                    end if
                    ! digit        
                    read (c, '(i1)') tmp
                    ! shift
                    IF( .not. truncating )then
                       IF( icount > 0 )then
                          integral = integral * 10
                       end if
                       ! add
                       integral = integral + tmp
                       ! increase the icount
                       icount = icount + 1
                    end if

                case default
                    IF( present(digit_count) )then
                        digit_count = icount
                    end if
                    call push_char(c)
                    integral = isign * integral
                    return
                end select
            end if
        end do

    end function parse_integer

    !
    ! POP CHAR
    !
    recursive character function pop_char(unit, str, eof, skip_ws) result(popped)
        integer, intent(in) :: unit
        character(*), intent(inout) :: str
        logical, intent(out) :: eof
        logical, intent(in), optional :: skip_ws

        integer :: ios
        character :: c
        logical :: ignore

        eof = .false.
        IF( .not.present(skip_ws) )then
            ignore = .false.
        else
            ignore = skip_ws
        end if

        do
            IF( pushed_index > 0 )then
                ! there is a character pushed back on, most likely from the number parsing                
                c = pushed_char(pushed_index:pushed_index)
                pushed_index = pushed_index - 1
                ios = 0
            else
                IF( unit .gt. 0 )then
                    read (unit = unit, fmt = "(a)", advance = "no", iostat = ios, EOR = 777, end = 999, err = 999) c
                else
                    read (unit = str, fmt = "(a)", iostat = ios) c
                    str = str(2:)
                endif
            end if
            IF( ios == end_of_record )then
777             cycle            
            else IF( ios == end_of_file )then
999             eof = .true.
                exit
            else IF( iachar(c) <= 31 )then
                ! non printing ascii characters
                cycle
            else IF( ignore .and. c == " " )then
                cycle
            else
                popped = c
                exit
            end if
        end do

    end function pop_char

    !
    ! PUSH CHAR
    !
    subroutine push_char(c)
        character, intent(inout) :: c
        pushed_index = pushed_index + 1
        pushed_char(pushed_index:pushed_index) = c

    end subroutine push_char

end module fson
