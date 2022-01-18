! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
subroutine Map_Jet_Plume
    
    use Variables_MPI
    use Variables_MPI_Mapping
    use GLOBAL
    
    implicit none
    
    ! *** Local variables
    integer :: NQJPIJ_Global
    integer :: ll, ii, jj, local
    
    ! *** save the global value
    NQJPIJ_Global = NQJPIJ
    
    ! *** Initialize to zero since now we will be determining each processes local value
    
    NQJPIJ = 0
    local = 0
    do ll = 1, NQJPIJ_Global
        ! ***get local I/J
        ii = IG2IL(IQJP_GL(ll)) 
        jj = JG2JL(JQJP_GL(ll)) 
        
        ! *** map to local values, including ghost cells
        if(ii > 0 .and. ii <= IC )then
            if(jj > 0 .and. jj <= JC )then
                ! *** 
                NQJPIJ = NQJPIJ + 1
                local = local + 1
                ! *** map to local I/J
                IQJP(local) = ii
                JQJP(local) = jj
                
                KQJP(local)  = KQJP_GL(ll)
                NPORTJP(local) = NPORTJP_GL(ll)
                XJETL(local) = XJETL_GL(ll)
                YJETL(local) = YJETL_GL(ll)
                ZJET(local)  = ZJET_GL(ll)
                PHJET(local) = PHJET_GL(ll)
                THJET(local) = THJET_GL(ll)
                DJET(local)  = DJET_GL(ll)
                CFRD(local)  = CFRD_GL(ll)
                DJPER(local) = DJPER_GL(ll)
 
            end if
        end if
    end do
    
    
end subroutine Map_Jet_Plume