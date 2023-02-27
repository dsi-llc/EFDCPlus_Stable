! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Defines mapping for 3d graphics grid.  Does a read of gcellmp.in
!! for some specifics about the mapping.  Not clear how often this feature is even used
! @date 9/12/2019
! @author Zander Mausolff
! Not sure this routine is even needed...
    
Subroutine Write_Mapping_3Dgraphics
    
    USE INFOMOD,ONLY:READSTR
    Use GLOBAL
    Use Variables_MPI
#ifdef _MPI
    Use Broadcast_Routines
#endif

    Implicit none
    
    ! *** Local variables
    Integer :: NW
    CHARACTER(200) :: STR
    Integer :: IGTMP,JGTMP,ICOMP,JCOMP,LTMP

    
! **  DEFINE MAPPING TO 3D GRAPHICS GRID
   !IF( ISCLO == 0 .OR. NWGG == 0 )THEN
   !    IG=IC
   !    JG=JC
   !ELSE
   !    ! *** ****************************
   !    if( process_id == master_id )THEN
   !        OPEN(1,FILE='gcellmp.inp',STATUS='UNKNOWN')
   !        STR = READSTR(1)
   !        READ(1,*)IG,JG
   !    endif
   !
   !    Call Broadcast_Scalar(IG, master_id)
   !    Call Broadcast_Scalar(JG, master_id)
   !
   !    DO NW=1,NWGG
   !        READ(1,*)IGTMP,JGTMP,ICOMP,JCOMP
   !        LTMP=LIJ(ICOMP,JCOMP)
   !        IWGG(NW)=IGTMP
   !        JWGG(NW)=JGTMP
   !        LWGG(NW)=LTMP
   !    ENDDO
   !    CLOSE(1)
   !ENDIF
    
End subroutine Write_Mapping_3Dgraphics