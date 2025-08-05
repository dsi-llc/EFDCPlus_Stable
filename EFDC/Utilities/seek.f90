! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SEEK(TAG, ISKIP)
  
  use GLOBAL,only: IKV,CARDNO
  use Variables_MPI
  
  implicit none
  
  integer,intent(IN)      :: ISKIP
  character,intent(INOUT) :: TAG*(*)

  integer :: I, J, K, L, M
  character*120 TEXT
  logical(4) :: OPN, ECHO
  
  ECHO = .TRUE.
  if( ISKIP > 0 ) ECHO = .FALSE.
  
  CARDNO = TAG
  INQUIRE(UNIT = mpi_efdc_out_unit,OPENED = OPN) 
  
  L = LEN(TAG)
  do I = 1,L
    J = ICHAR(TAG(I:I))
    if( 97 <= J .and. J <= 122 )then
      TAG(I:I) = CHAR(J-32)
    endif
  enddo  
  if( OPN ) WRITE(mpi_efdc_out_unit,'(A,A)')'SEEKING GROUP: ',TAG
  
  do K = 1,2
10  read(1,'(A)',END = 20)TEXT
    M = MAX(1,LEN_TRIM(TEXT))
    if( OPN .and. ECHO ) WRITE(mpi_efdc_out_unit,'(A)') TEXT(1:M)
    do while( M > L .and. TEXT(1:1) == '' )
      TEXT(1:M-1) = TEXT(2:M)
      TEXT(M:M) = ' '
      M = M-1
    enddo
    if( M < L )GO TO 10
    do I = 1,M
      J = ICHAR(TEXT(I:I))
      if( 97 <= J .and. J <= 122 )then
        TEXT(I:I) = CHAR(J-32)
      endif
    enddo
    if( TEXT(1:L) /= TAG )     GO TO 10
    if( TEXT(L+1:L+1) /= ' ' ) GO TO 10
  enddo

  return
 
  20 WRITE(*,'(A,A,A)') 'GROUP: ',TAG,' NOT FOUND BEFORE END OF FILE'
  write(mpi_error_file,'(A,A,A)') 'GROUP: ',TAG,' NOT FOUND BEFORE END OF FILE'
  
  PAUSE
  
  call STOPP('.')

END

