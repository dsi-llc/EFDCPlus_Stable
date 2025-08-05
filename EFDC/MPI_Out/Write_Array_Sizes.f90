! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details writes arrays sizes to file, previously in Varalloc.f90
! @author Zander Mausolff
! @date 7/31/2019
    
Subroutine Write_Array_Sizes

    use GLOBAL
    use Variables_WQ
    
    use Variables_MPI

    implicit none
   
    integer      :: Array_Sizes_Unit
    Character(3) :: process_id_txt, STR*200
    
    Array_Sizes_Unit = 888 + process_id
    write(process_id_txt, "(I3.3)") process_id 
    open(Array_Sizes_Unit, FILE = OUTDIR//'log_arraysize_proc_'//process_id_txt//'.log',STATUS = 'UNKNOWN')
    
    write(STR, '("*** EFDC+ Filename: ",A,",  Version: ",A10," ***")') TRIM(EFDC_EXE), EFDC_VER
    write(Array_Sizes_Unit, '(A,/,A,/,A,/)') REPEAT('*',LEN_TRIM(STR)), TRIM(STR), REPEAT('*',LEN_TRIM(STR))
    
    ! *** Report array sizes
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,//)') RUNTITLE
    
    ! *** Grid
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM I CELL DIRECTION (ICM): ',ICM
    write(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM J CELL DIRECTION (JCM): ',JCM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF LAYERS (KCM): ',KCM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF ACTIVE CELLS (LCM): ',LCM
    
    ! *** Boundary Conditions
    
    ! *** Flow
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF FLOW BOUNDARY CELLS (NQSIJ): ', NQSIJ
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF FLOW TIME SERIES (NQSERM): ', NQSERM
    
    ! *** Jet/Plume
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF JET-PLUME SOURCES (NQJPIJ): ', NQJPIJ
    
    ! *** Hydraulic Structure
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF EQUATION DEFINED HYDRAULIC STRUCTURE GROUPS (NHYDST): ', NHYDST
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF PRESSURE CONTROLLED UPSTREAM/DOWNSTEAM SERIES (NQCTTM): ',NQCTTM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF TABLE DATA POINTS (NDQCLT): ', NDQCLT
    
    ! *** Withdrawal/Return
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'MAXIMUM WITHDRAWAL/RETURN SERIES (NQWRM): ',NQWRM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF WITHDRAWAL/RETURN SERIES (NQWRSRM): ', NQWRSRM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF POINTS IN A WITHDRAWAL/RETURN SERIES (NDQWRSR): ', NDQWRSR
    
    ! *** Open Boundary
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY ELEVATION/PRESSURE CELLS ON EAST  BOUNDARY (NPBE): ', NPBE
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY ELEVATION/PRESSURE CELLS ON NORTH BOUNDARY (NPBN): ', NPBN
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY ELEVATION/PRESSURE CELLS ON SOUTH BOUNDARY (NPBS): ', NPBS
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY ELEVATION/PRESSURE CELLS ON WEST  BOUNDARY (NPBW): ', NPBW
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CONCENTRATION CELLS ON EAST  BOUNDARY (NCBE): ', NCBE
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CONCENTRATION CELLS ON NORTH BOUNDARY (NCBN): ', NCBN
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CONCENTRATION CELLS ON SOUTH BOUNDARY (NCBS): ', NCBS
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY CONCENTRATION CELLS ON WEST  BOUNDARY (NCBW): ', NCBW
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF OPEN BOUNDARY SURFACE ELEVATION TIME SERIES (NPSER): ', NPSER
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF DATA POINTS IN ELEVATION/PRESSURE TIME SERIES (NDPSER): ', NDPSER
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF HARMONIC TIDAL CONSITIUENTS (MTM): ',MTM

    ! *** Atmospheric/Wind/Heat
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF ATMOSPHERIC CONDITION TIME SERIES (NASER): ', NASER
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF WIND TIME SERIES (NWSER): ', NWSER
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF ICE INPUT SERIES (NISER): ', NISER
    
    ! *** Groundwater
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF GROUNDWATER SERIES (NGWSERM): ',NGWSERM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF POINTS IN GROUNDWATER SERIES (NDGWSER): ',NDGWSER

    ! *** Constituents
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF ALL CONCENTRATION TIME SERIES (NCSERM): ', NCSERM

    ! *** Dye
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF DYE CLASSES (NDYM): ',NDYM

    ! *** Sediments  
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF SEDIMENT BED LAYERS (KB): ', KB
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF COHESIVES OR ALL SEDIMENTS FOR SEDZLJ (NSCM): ',NSCM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF NON-COHESIVES (NSNM): ',NSNM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF COHESIVES + NON-COHESIVE S(NSTM): ', NSTM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF REDOPSITED GRAINSIZE CLASSES DEFINED (NSICM): ', NSICM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF BANK EROSION CELLS(NBEPAIRM): ', NBEPAIRM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF BED EROSION TIME SERIES(NBESERM): ', NBESERM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF DATA POINTS FOR BED EROSION TIME SERIES(NDBESER): ', NDBESER
    
    ! *** Toxic
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF TOXICS (NTXM): ',NTXM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF PARTICLE MIXING TIME SERIES(NPMXPTSM): ', NPMXPTSM
    
    ! *** Water Quality
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF WATER QUALITY CONSTITUENTS (NWQV): ',NWQV
    write(Array_Sizes_Unit,'(A,I10)') 'TOTAL NUMBER OF ALL WC CONSTITUENTS (NSTVM): ',NSTVM

    ! *** Vegetation
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF VEGETATION CLASSES (NVEGTPM): ', NVEGTPM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF VEGETATION TIME SERIES (NVEGSERM): ', NVEGSERM

    ! *** MHK Devices
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF CELLS WITH MHK DEVICES (TCOUNT): ', TCOUNT
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF DIFFERENT TYPES OF MHK DEVICES (MHKTYP): ', MHKTYP
    
    ! *** Output
    ! *** WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF TIME SERIES OUTPUT LOCATIONS (NSMTS): ', NSMTS
    ! *** WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF HORIZONTAL LOCATIONS TO WRITE TIME SERIES OF SURF ELEV,VELOCITY, AND CONCENTRATION VARIABLES (MLTMSRM): ', MLTMSRM
    
    ! *** MPI
    write(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF BOUNDARY CONDITION CELLS (NMAXBC): ', NMAXBC
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF CONCENTRATION BOUNDARY CONDITIONS CELLS ON EAST OPEN BOUNDARIES (NBBEM): ', NBBEM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF CONCENTRATION BOUNDARY CONDITIONS CELLS ON NORTH OPEN BOUNDARIES (NBBNM): ', NBBNM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF CONCENTRATION BOUNDARY CONDITIONS CELLS ON SOUTH OPEN BOUNDARIES (NBBSM): ', NBBSM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF CONCENTRATION BOUNDARY CONDITIONS CELLS ON WEST OPEN BOUNDARIES (NBBWM): ', NBBWM
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF SUB GRID CHANNEL (NCHANM): ', NCHANM
    
    ! *** Misc.
    write(Array_Sizes_Unit,'(A,I10)')
    ! *** WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF TIME SERIES START-STOP SCENARIOS (NTSSTSP): ', NTSSTSP
    ! *** WRITE(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF STOP-START PAIRS FOR SCENARIO ISSS (MTSSTSP): ', MTSSTSP
    ! *** WRITE(Array_Sizes_Unit,'(A,I10)')
    write(Array_Sizes_Unit,'(A,I10)') 'NUMBER OF CELLS WHOSE U AND/OR V FACES ARE PARTIALLY BLOCKED (NBLOCKED): ', NBLOCKED
    
    close(Array_Sizes_Unit)

End subroutine Write_Array_Sizes
