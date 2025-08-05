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
! @details This routine reads cell.inp and sets up IJCT
!! Not sure if this is needed or even used in O'Donncha's code...
!! Does setup LCM as LCM = LCM+4 --> this is utilized elsewhere.
! @author Zander Mausolff
! @date 9/2/2019
!---------------------------------------------------------------------------!
Subroutine Scan_Cell

  use GLOBAL
  use INFOMOD
  use Variables_MPI

  use Broadcast_Routines
  use MPI

  implicit none

  ! *** Read in variables

  ! *** Local variables
  integer ::  L, ierr
  integer, Allocatable, Dimension(:) :: IB
  integer, Allocatable, Dimension(:) :: IE
  integer, Allocatable, Dimension(:) :: JB
  integer, Allocatable, Dimension(:) :: JE
  integer, Allocatable, Dimension(:,:) :: IJCT_Read_In

  integer :: IS, IACROSS, IP, IT, IFIRST, J, I, JDUMY, ILAST, ISO
  integer :: NJ, II, NI, JBEG, JEND, IBEG, IEND, jj, ICUR, JCUR
  integer :: JF, JT, JLAST, JACROSS

  character :: STRC*650
  
  !---------------------------------------------------------------------------!

  allocate(IB(n_x_partitions))
  allocate(IE(n_x_partitions))
  allocate(JB(n_y_partitions))
  allocate(JE(n_y_partitions))
  allocate(IB_Decomp(0:active_domains))
  allocate(IE_Decomp(0:active_domains))
  allocate(JB_Decomp(0:active_domains))
  allocate(JE_Decomp(0:active_domains))

  allocate(IJCT_Read_In(IC,JC))
  
  MPI_Write_Flag = DEBUG           ! *** Write out MPI details
  
  ! *** Zero it out
  IB_Decomp = 0
  IE_Decomp = 0
  JB_Decomp = 0
  JE_Decomp = 0
  IJCT_Read_In(:,:) = 0

  ! *** Read the cell.inp file
  if( process_id == master_id )then
    ! *** read in CELL.INP data and allocate LCM based on maximum
    open(1,FILE = 'cell.inp',STATUS = 'UNKNOWN')

    STRC = READSTR(1)
    read(STRC,*) JDUMY

    if( JDUMY /= JC )then                                                                                    
      ! ***   READ OLD FILE FORMAT                                                                                                                                                                               
      JACROSS = JC
      if( JC > 640 ) JACROSS = 640
      do JT = 1,JC,JACROSS
        JF = JT
        JLAST = JT + JACROSS - 1
        if( JLAST > JC ) JLAST = JC
        WRITE (777,8) JF,JLAST
        do I = 1,IC
          read(1,6,IOSTAT = ISO) (IJCT_Read_In(I,J),J = JF,JLAST)
          if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
          WRITE (777,66) (IJCT_Read_In(I,J),J = JF,JLAST)
        enddo
      enddo

    else

      ! *** "New" format of the CELL.INP file
      if( IC > 640 )then
        IACROSS = 640
        do IT = 1,IC,IACROSS
          IFIRST = IT
          ILAST = IT + IACROSS - 1
          if( ILAST > IC ) ILAST = IC
          do J = JC,1,-1
            read(1,66,IOSTAT = ISO) JDUMY,(IJCT_Read_In(I,J),I = IFIRST,ILAST)
            if(ISO.GT.0 )then
              write(6,*)'  READ ERROR FOR FILE CELL.INP '
              STOP
            endif
          enddo
        enddo

      else !***IC < 640
      
        IFIRST = 1
        ILAST  = IC
        do J = JC,1,-1
          read(1,66,IOSTAT = ISO) JDUMY, (IJCT_Read_In(I,J),I = IFIRST,ILAST)
          if(ISO.GT.0 )then
            write(6,*) '  READ ERROR FOR FILE CELL.INP '
            STOP
          endif
        enddo
      endif
    endif
 66 FORMAT (I4,1X,640I1)
    
    close(1) !***Close the cell.inp file
  endif

  ! *** Broadcast the IJCT values
  call MPI_BCAST(IJCT_Read_In, size(IJCT_Read_In), MPI_Integer, master_id, comm_2d, ierr)
  call MPI_BARRIER(comm_2d, ierr)

  IB(1) = 1
  IE(1) = ic_decomp(1) 
  do NI = 2,n_x_partitions
    IB(NI) = IE(NI-1) + 1                      ! *** Beginning global I for current domain, without ghost cells
    IE(NI) = IB(NI)   + ic_decomp(NI) - 1      ! *** Ending    global I for current domain, without ghost cells
  enddo

  JB(1) = 1
  JE(1) = jc_decomp(1)
  do NJ = 2,n_y_partitions
    JB(NJ) = JE(NJ-1) + 1                      ! *** Beginning global J for current domain, without ghost cells
    JE(NJ) = JB(NJ)   + jc_decomp(NJ) - 1      ! *** Ending    global J for current domain, without ghost cells
  enddo

  write(mpi_log_unit, '(A)') 'Number of Active Cells in Each Partition (including ghost cells)' 
  write(mpi_log_unit, '(6a5,a10)') 'X','Y','IS','IE','JS','JE','LCM'

  LCM = 0
  ii  = 0
  IP = -1
  ! *** Get the starting and ending I an J's for each sub-domain, including ghost cells
  do NJ = 1,n_y_partitions
    do NI = 1,n_x_partitions
      L = 0
      IBEG = -1
      IEND = -1
      JBEG = -1
      JEND = -1
      if( process_map(NI,NJ) /= -1 )then
        IP = IP + 1
        
        ! *** Only need the J component once
        JBEG = JB(NJ)
        JEND = JE(NJ)
        if( process_map(NI,NJ-1) /= -1 )  JBEG = JBEG - n_ghost_rows   ! *** South Edge
        if( process_map(NI,NJ+1) /= -1 )  JEND = JEND + n_ghost_rows   ! *** North Edge
        JB_Decomp(IP) = JBEG
        JE_Decomp(IP) = JEND
    
        IBEG = IB(NI)
        IEND = IE(NI)
        if( process_map(NI-1,NJ) /= -1 )  IBEG = IBEG - n_ghost_rows   ! *** West Edge
        if( process_map(NI+1,NJ) /= -1 )  IEND = IEND + n_ghost_rows   ! *** East Edge
        IB_Decomp(IP) = IBEG
        IE_Decomp(IP) = IEND
      
        L = 1
        do J = JBEG,JEND
          do I = IBEG,IEND
            if( IJCT_Read_In(i,j) > 0 .and. IJCT_Read_In(i,j) < 9 )then
              L = L + 1
            endif
          enddo
        enddo
      endif
      
      write(mpi_log_unit, '(6I5,I10)') NI, NJ, IBEG, IEND, JBEG, JEND, L
      
      LCM = max(LCM,L)
    enddo
  enddo

  LCM = LCM + 4

  call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)')    'Writing out global grid parameters in SCANCELL Routine'
  write(mpi_log_unit, '(a,I6)') 'Max local LCM = ', LCM
  write(mpi_log_unit, '(a,I6)') 'Global IC     = ', IC
  write(mpi_log_unit, '(a,I6)') 'Global JC     = ', JC
  
  call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)')    'Writing out Local grid parameters for the respective subdomain'
  write(mpi_log_unit, '(a,I6)') 'Starting I cell =  ', IB_Decomp(process_id)
  write(mpi_log_unit, '(a,I6)') 'Ending I cell   =  ', IE_Decomp(process_id)
  write(mpi_log_unit, '(a,I6)') 'Starting J cell =  ', JB_Decomp(process_id)
  write(mpi_log_unit, '(a,I6)') 'Ending J cell   =  ', JE_Decomp(process_id)
  write(mpi_log_unit, '(a)')    'End of Scan_Cell routine'
  call WriteBreak(mpi_log_unit)

  deallocate(IJCT_Read_In)
  deallocate(IB,IE)
  deallocate(JB,JE)

5 FORMAT(10I5)
6 FORMAT(A10,40I5)
8 FORMAT ('   CELL TYPE ARRAY,J = ',I5,2X,'TO J = ',I5,//)                                                              

End Subroutine Scan_Cell
