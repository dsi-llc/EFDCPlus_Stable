! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details
! @author
! @date
    
Subroutine WriteCellMapLogging

    use GLOBAL
    use Variables_MPI
    
    implicit none
    
    ! *** Local variables 
    integer :: L
    real :: GRADW,GRADE,GRADS,GRADN

    ! **********************************************************************************
    ! *** EFDC CELL MAP LOGGING
    write(mpi_efdc_out_unit,'(3A6,4A10,4A14)') 'I','J','L','BELV','HP','DX','DY','GRADW','GRADE','GRADS','GRADN'
    do L = 2,LA
        GRADW = -999.
        GRADE = -999.
        GRADS = -999.
        GRADN = -999.
        if( SUBO(L) > 0.      ) GRADW = ( BELV(LWC(L))+HP(LWC(L) ) - ( BELV(L)+HP(L) )           )/DXU(L)
        if( SUBO(LEC(L)) > 0. ) GRADE = ( BELV(L)+HP(L)            - ( BELV(LEC(L))+HP(LEC(L)) ) )/DXU(LEC(L))
        if( SVBO(L) > 0.      ) GRADS = ( BELV(LSC(L))+HP(LSC(L) ) - ( BELV(L)+HP(L) )           )/DYV(L)
        if( SVBO(LNC(L)) > 0. ) GRADN = ( BELV(L)+HP(L)            - ( BELV(LNC(L))+HP(LNC(L)) ) )/DYV(LNC(L))
        write(mpi_efdc_out_unit,'(3I6,4F10.3,4F14.5)') IL(L),JL(L),L,BELV(L),HP(L),DXP(L),DYP(L),GRADW,GRADE,GRADS,GRADN
    enddo
    
End subroutine WriteCellMapLogging
