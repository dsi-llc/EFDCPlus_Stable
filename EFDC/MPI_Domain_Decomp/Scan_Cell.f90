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
! @details This routine reads cell.inp and sets up IJCT
!! Not sure if this is needed or even used in O'Donncha's code...
!! Does setup LCM as LCM=LCM+4 --> this is utilized elsewhere.
! @author Zander Mausolff
! @date 9/2/2019
!---------------------------------------------------------------------------!
Subroutine Scan_Cell

  Use GLOBAL
  Use INFOMOD
  Use Variables_MPI

#ifdef _MPI
  Use Broadcast_Routines
  USE MPI
#endif
  Implicit None

  ! *** Read in variables

  ! *** Local variables
  Integer ::  L, ierr
  Integer, Allocatable, Dimension(:) :: IB
  Integer, Allocatable, Dimension(:) :: IE
  Integer, Allocatable, Dimension(:) :: JB
  Integer, Allocatable, Dimension(:) :: JE
  Integer, Allocatable, Dimension(:,:) :: IJCT_Read_In

  INTEGER :: IS, IACROSS, IP, IT, IFIRST, J, I, JDUMY, ILAST, ISO
  INTEGER :: NJ, II, NI, JBEG, JEND, IBEG, IEND, jj, ICUR, JCUR
  INTEGER :: JF, JT, JLAST, JACROSS

  CHARACTER :: STRC*650
  
  !---------------------------------------------------------------------------!

  ALLOCATE(IB(n_x_partitions))
  ALLOCATE(IE(n_x_partitions))
  ALLOCATE(JB(n_y_partitions))
  ALLOCATE(JE(n_y_partitions))
  ALLOCATE(IB_Decomp(0:active_domains))
  ALLOCATE(IE_Decomp(0:active_domains))
  ALLOCATE(JB_Decomp(0:active_domains))
  ALLOCATE(JE_Decomp(0:active_domains))

  ALLOCATE(IJCT_Read_In(IC,JC))

  ! *** Zero it out
  IB_Decomp = 0
  IE_Decomp = 0
  JB_Decomp = 0
  JE_Decomp = 0
  IJCT_Read_In(:,:) = 0

  ! *** Read the cell.inp file
  if( process_id == master_id )THEN
    ! *** read in CELL.INP data and allocate LCM based on maximum
    OPEN(1,FILE='cell.inp',STATUS='UNKNOWN')

    STRC = READSTR(1)
    READ(STRC,*) JDUMY

    IF( JDUMY /= JC )THEN                                                                                    
      ! ***   READ OLD FILE FORMAT                                                                                                                                                                               
      JACROSS = JC
      IF( JC > 640 ) JACROSS = 640
      DO JT=1,JC,JACROSS
        JF = JT
        JLAST = JT + JACROSS - 1
        IF( JLAST > JC ) JLAST = JC
        WRITE (777,8) JF,JLAST
        DO I=1,IC
          READ(1,6,IOSTAT=ISO) (IJCT_Read_In(I,J),J=JF,JLAST)
          IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE CELL.INP')
          WRITE (777,66) (IJCT_Read_In(I,J),J=JF,JLAST)
        ENDDO
      ENDDO

    ELSE

      ! *** "New" format of the CELL.INP file
      IF( IC > 640 )THEN
        IACROSS = 640
        DO IT=1,IC,IACROSS
          IFIRST = IT
          ILAST = IT + IACROSS - 1
          IF( ILAST > IC ) ILAST = IC
          DO J=JC,1,-1
            READ(1,66,IOSTAT=ISO) JDUMY,(IJCT_Read_In(I,J),I=IFIRST,ILAST)
            IF(ISO.GT.0 )THEN
              WRITE(6,*)'  READ ERROR FOR FILE CELL.INP '
              STOP
            END IF
          ENDDO
        ENDDO

      ELSE !***IC < 640
      
        IFIRST = 1
        ILAST  = IC
        DO J=JC,1,-1
          READ(1,66,IOSTAT=ISO) JDUMY, (IJCT_Read_In(I,J),I=IFIRST,ILAST)
          IF(ISO.GT.0 )THEN
            WRITE(6,*) '  READ ERROR FOR FILE CELL.INP '
            STOP
          END IF
        ENDDO
      ENDIF
    ENDIF
 66 FORMAT (I4,1X,640I1)
    
    CLOSE(1) !***Close the cell.inp file
  end if

  ! *** Broadcast the IJCT values
  Call MPI_BCAST(IJCT_Read_In, size(IJCT_Read_In), MPI_Integer, master_id, comm_2d, ierr)
  Call MPI_BARRIER(comm_2d, ierr)

  IB(1) = 1
  IE(1) = ic_decomp(1) 
  DO NI = 2,n_x_partitions
    IB(NI) = IE(NI-1) + 1                      ! *** Beginning global I for current domain, without ghost cells
    IE(NI) = IB(NI)   + ic_decomp(NI) - 1      ! *** Ending    global I for current domain, without ghost cells
  END DO

  JB(1) = 1
  JE(1) = jc_decomp(1)
  DO NJ = 2,n_y_partitions
    JB(NJ) = JE(NJ-1) + 1                      ! *** Beginning global J for current domain, without ghost cells
    JE(NJ) = JB(NJ)   + jc_decomp(NJ) - 1      ! *** Ending    global J for current domain, without ghost cells
  END DO

  write(mpi_log_unit, '(A)') 'Number of Active Cells in Each Partition (including ghost cells)' 
  write(mpi_log_unit, '(6a5,a10)') 'X','Y','IS','IE','JS','JE','LCM'

  LCM = 0
  ii  = 0
  IP = -1
  ! *** Get the starting and ending I an J's for each sub-domain, including ghost cells
  DO NJ = 1,n_y_partitions
    DO NI = 1,n_x_partitions
      L = 0
      IBEG = -1
      IEND = -1
      JBEG = -1
      JEND = -1
      IF( process_map(NI,NJ) /= -1 )THEN
        IP = IP + 1
        
        ! *** Only need the J component once
        JBEG = JB(NJ)
        JEND = JE(NJ)
        IF( process_map(NI,NJ-1) /= -1 )  JBEG = JBEG - n_ghost_rows   ! *** South Edge
        IF( process_map(NI,NJ+1) /= -1 )  JEND = JEND + n_ghost_rows   ! *** North Edge
        JB_Decomp(IP) = JBEG
        JE_Decomp(IP) = JEND
    
        IBEG = IB(NI)
        IEND = IE(NI)
        IF( process_map(NI-1,NJ) /= -1 )  IBEG = IBEG - n_ghost_rows   ! *** West Edge
        IF( process_map(NI+1,NJ) /= -1 )  IEND = IEND + n_ghost_rows   ! *** East Edge
        IB_Decomp(IP) = IBEG
        IE_Decomp(IP) = IEND
      
        L = 1
        DO J = JBEG,JEND
          DO I = IBEG,IEND
            IF( IJCT_Read_In(i,j) > 0 .AND. IJCT_Read_In(i,j) < 9 )THEN
              L = L + 1
            End if
          END DO
        END DO
      ENDIF
      
      write(mpi_log_unit, '(6I5,I10)') NI, NJ, IBEG, IEND, JBEG, JEND, L
      
      LCM = max(LCM,L)
    END DO
  END DO

  LCM = LCM + 4

  Call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)')    'Writing out global grid parameters in SCANCELL Routine'
  write(mpi_log_unit, '(a,I6)') 'Max local LCM = ', LCM
  write(mpi_log_unit, '(a,I6)') 'Global IC     = ', IC
  write(mpi_log_unit, '(a,I6)') 'Global JC     = ', JC
  
  Call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)')    'Writing out Local grid parameters for the respective subdomain'
  write(mpi_log_unit, '(a,I6)') 'Starting I cell =  ', IB_Decomp(process_id)
  write(mpi_log_unit, '(a,I6)') 'Ending I cell   =  ', IE_Decomp(process_id)
  write(mpi_log_unit, '(a,I6)') 'Starting J cell =  ', JB_Decomp(process_id)
  write(mpi_log_unit, '(a,I6)') 'Ending J cell   =  ', JE_Decomp(process_id)
  write(mpi_log_unit, '(a)')    'End of Scan_Cell routine'
  Call WriteBreak(mpi_log_unit)

  DEALLOCATE(IJCT_Read_In)
  DEALLOCATE(IB,IE)
  DEALLOCATE(JB,JE)

5 FORMAT(10I5)
6 FORMAT(A10,40I5)
8 FORMAT ('   CELL TYPE ARRAY,J=',I5,2X,'TO J=',I5,//)                                                              

End Subroutine Scan_Cell
