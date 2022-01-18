! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details
! @author
! @date
    
Subroutine WriteCellMapLogging

    Use GLOBAL
    Use Variables_MPI
    
    Implicit None
    
    ! *** Local variables 
    Integer :: L
    Real :: GRADW,GRADE,GRADS,GRADN

    ! **********************************************************************************
    ! *** EFDC CELL MAP LOGGING
    WRITE(8,'(3A6,4A10,4A14)') 'I','J','L','BELV','HP','DX','DY','GRADW','GRADE','GRADS','GRADN'
    DO L=2,LA
        GRADW = -999.
        GRADE = -999.
        GRADS = -999.
        GRADN = -999.
        IF( SUBO(L) > 0.      ) GRADW = ( BELV(LWC(L))+HP(LWC(L) ) - ( BELV(L)+HP(L) )           )/DXU(L)
        IF( SUBO(LEC(L)) > 0. ) GRADE = ( BELV(L)+HP(L)            - ( BELV(LEC(L))+HP(LEC(L)) ) )/DXU(LEC(L))
        IF( SVBO(L) > 0.      ) GRADS = ( BELV(LSC(L))+HP(LSC(L) ) - ( BELV(L)+HP(L) )           )/DYV(L)
        IF( SVBO(LNC(L)) > 0. ) GRADN = ( BELV(L)+HP(L)            - ( BELV(LNC(L))+HP(LNC(L)) ) )/DYV(LNC(L))
        WRITE(8,'(3I6,4F10.3,4F14.5)') IL(L),JL(L),L,BELV(L),HP(L),DXP(L),DYP(L),GRADW,GRADE,GRADS,GRADN
    ENDDO
    
End subroutine WriteCellMapLogging
