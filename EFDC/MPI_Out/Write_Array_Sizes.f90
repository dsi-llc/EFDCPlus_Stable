! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details writes arrays sizes to file, previously in Varalloc.f90
! @author Zander Mausolff
! @date 7/31/2019
    
Subroutine Write_Array_Sizes

    USE GLOBAL
    Use Variables_WQ
    
    USE Variables_MPI

    Implicit None
   
    Integer      :: Array_Sizes_Unit
    Character(3) :: process_id_txt
    
    WRITE(process_id_txt, "(I3.3)") process_id 

    Array_Sizes_Unit = 888 + process_id
    
    ! *** REPORT ARRAY SIZES
    OPEN( Array_Sizes_Unit, FILE=OUTDIR//'ARRAYSIZES_'//process_id_txt//'.LOG',STATUS='UNKNOWN')
    WRITE(Array_Sizes_Unit,'(A,//)') RUNTITLE
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM L ARRAY (LCM): ',LCM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM K ARRAY (KCM): ',KCM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM I ARRAY (ICM): ',ICM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM J ARRAY (JCM): ',JCM
    WRITE(Array_Sizes_Unit,'(A,I10)')
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF COHESIVES (NSCM): ',NSCM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF NON-COHESIVES (NSNM): ',NSNM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF TOXICS (NTXM): ',NTXM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF WATER QUALITY CONSTITUENTS (NWQV): ',NWQV
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM NUMBER OF CONSTITUENTS (NSTVM): ',NSTVM
    WRITE(Array_Sizes_Unit,'(A,I10)')
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM WITHDRAWAL/RETURN SERIES (NQWRM): ',NQWRM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM MTIDE CONSITIUENTS (MTM): ',MTM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CELLS: SOUTH (NPBS): ',NPBS
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CELLS: WEST  (NPBW): ',NPBW
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CELLS: EAST  (NPBE): ',NPBE
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CELLS: NORTH (NPBN): ',NPBN
    WRITE(Array_Sizes_Unit,'(A,I10)')
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF 3D FIELD OUTPUT UNSTRETCHED LAYERS (KPCM): ',KPCM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF PRESSURE CONTROLLED WITHDRAWAL/RETURN PAIRS (NQCTTM): ',NQCTTM
    WRITE(Array_Sizes_Unit,'(A,I10)')
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF GROUNDWATER SERIES (NGWSERM): ',NGWSERM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM NUMBER OF POINTS IN GROUNDWATER SERIES (NDGWSER): ',NDGWSER
    WRITE(Array_Sizes_Unit,'(A,I10)')
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER (NQCTTM): ',NQCTTM
    WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER (KSM): ',KSM

    ! DELME - PLEASE FINISH  
    CLOSE(Array_Sizes_Unit)

End subroutine Write_Array_Sizes
