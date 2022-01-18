! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Module contains many subroutines for writing out specific
!! arrays, modified for MPI, during the calculation process
! @date 9/12/2019
! @author Zander Mausolff
    
Module Mod_Write_Cell_Map
    
    Use Variables_MPI
    Use GLOBAL
    
    Implicit None
    
    ! *** Local variables
    Integer,Private :: I, J, L
    
    Contains

!**********************************************************************
    Subroutine Write_LIJ(write_out)
    
    Implicit none
    
    ! *** Read in
    Logical, intent(in) :: write_out
    
    ! *** Local
    
    if(write_out == .TRUE. )THEN
    
        Call WriteBreak(mpi_log_unit)

        write(mpi_log_unit, '(A)') ' ' 
        write(mpi_log_unit, '(A)') 'IJCT Global: '
        ! *** Write out XLOC-YLOC - local to a process
        do j = jc_global,1,-1
            write(mpi_log_unit, 900) (IJCT_Global(i,j), i=1, ic_global)
        end do
        write(mpi_log_unit, '(A)') ' '
        Call WriteBreak(mpi_log_unit)
        
        write(mpi_log_unit, '(A)') ' ' 
        write(mpi_log_unit, '(A)') 'IJCT Local: '
        ! *** Write out IJCT - local to a process
        do j = jc,1,-1 
            write(mpi_log_unit, 900) (IJCT(i,j), i=1, ic)
        end do
        write(mpi_log_unit, '(A)') ' '
        !**Format for the 2d matrix, max width of 100   
        900 format (100(I3)) 
        
        ! *** Write out to MPI debug
        write(mpi_log_unit, '(a)') ' ' 
        write(mpi_log_unit, '(a)') 'In CellMap Routine'
        write(mpi_log_unit, '(a,i5,i5)') 'IGRID=',IL2IG(1), IL2IG(IC)
        write(mpi_log_unit, '(a,i5)') 'LA Global = ', LA_Global
        ! *** End new section for MPI domain decomp
        write(mpi_log_unit, '(a,i5)') 'LA  = ', LA
        write(mpi_log_unit, '(a,i5)') 'LC  = ', LC 
        write(mpi_log_unit, '(a,i5)') 'IC  = ', IC
        write(mpi_log_unit, '(a,i5)') 'JC  = ', JC 
        
        Call WriteBreak(mpi_log_unit)
        write(mpi_log_unit, '(a)') 'Global L indexing: '
        do j = JC_Global,1,-1
            write(mpi_log_unit, '(100I5)',advance="no") (lij_global(i,j), i= 1, ic_global)
            write(mpi_log_unit, ' ' )
        end do
    
        Call WriteBreak(mpi_log_unit)
        write(mpi_log_unit, '(a)') 'Local L indexing: '
        do j = JC, 1, -1
            write(mpi_log_unit, '(100I5)',advance="no") (lij(i,j), i=1, IC)
            write(mpi_log_unit, ' ' )
        end do
        Call WriteBreak(mpi_log_unit)
        
    end if
    
    End subroutine Write_LIJ
!**********************************************************************
! @details Writes out LIJ global, LIJ local, and boundary mapping arrays
    Subroutine Write_Cell_Indexing(write_out)
    
    Implicit None
    
    ! *** Read in variable
    Logical, intent(in) :: write_out !*** T/F T- write out, F- do not write to a file 
    
    if(write_out == .TRUE. )THEN
        ! *** ****************************
        ! *** Write EFDC CELL MAP LOGGING to each processes file
        WRITE(mpi_log_unit,'(7A6)') 'I','J','L','LWC','LEC','LSC','LNC'
        DO L=2,LA
            WRITE(mpi_log_unit,'(7I6)') IL(L),JL(L),L,LWC(L),LEC(L),LSC(L),LNC(L)
        ENDDO
        
    End if
    
    End Subroutine Write_Cell_Indexing
!**********************************************************************
! @details writes out xloc and yloc depending on if MPI_Write_Flag is turned on or off
!
    Subroutine Write_Xloc_Yloc
    
        Implicit none
        
        if(MPI_Write_Flag == .TRUE. )THEN
        
            Call WriteBreak(mpi_log_unit)
            
            write(mpi_log_unit, '(a)') 'IG2IL(i):'
        
            do i =1, ic_global
                write(mpi_log_unit, '(2I5)') I, IG2IL(I)
            end do
            
            write(mpi_log_unit, '(a)') ' '
            
            write(mpi_log_unit, '(a)') 'JG2JL(j)'
            do j = 1, jc_global
                write(mpi_log_unit, '(2I5)') J, JG2JL(J)
            end do
            
            Call WriteBreak(mpi_log_unit)
            
            ! *** Write out XLOC-YLOC as 2D matrix
            !do j = jc_global,1,-1 
            !    do i = 1, ic_global
            !        write(mpi_log_unit, '(A)', ADVANCE='NO') '('
            !        write(mpi_log_unit, 900,ADVANCE='NO'   )  IG2IL(i), JG2JL(j)           
            !        write(mpi_log_unit, '(A)', ADVANCE='NO') ')'
            !    end do
            !    write(mpi_log_unit, '(A)') ' ' 
            !end do
    
        End if
    
    End subroutine Write_Xloc_Yloc
    
End module Mod_Write_Cell_Map
