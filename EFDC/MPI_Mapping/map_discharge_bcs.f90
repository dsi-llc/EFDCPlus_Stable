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
! @details Maps flow type boundary condition values to local ones
! @author Zander Mausolff, adapted from O'Donncha's code
! @date 9/4/2019
!---------------------------------------------------------------------------!

Subroutine Map_River

  use GLOBAL
  use HYDSTRUCMOD

  use MPI
  use Variables_MPI
  use Variables_MPI_Mapping
  use MPI_All_Reduce

  implicit none

  ! *** Local
  real(8) :: TTDS, TWAIT
  integer :: II, LL, MS, K, III, IERR, JJJ, L, M
  integer :: MMAX, MMIN
  integer :: NQSIJ_GL
  integer :: LLSave(NQSIJM)

  NQSIJ_GL = NQSIJ
  LLSave = 0

  ! *** Determine the local NQSIJ value
  NGRPID = 0
  NQSIJ = 0
  II = 0
  do LL = 1, NQSIJ_GL
    BCPS_GL(LL).L = LIJ_Global(BCPS_GL(LL).I,BCPS_GL(LL).J)
    if( BCPS_GL(LL).L < 2 )then
      write(6,*) 'ERROR! FLOW CELL IS NOT IN GLOBAL DOMAIN FOR FLOW BOUNDARY: ',LL
      call STOPP('', 1)
    endif

    III = IG2IL(BCPS_GL(LL).I)
    JJJ = JG2JL(BCPS_GL(LL).J)
    if( III > 0 .and. III <= IC )then     ! *** Allow ghost cells containing inflow
      if( JJJ > 0 .and. JJJ <= JC )then   ! *** Allow ghost cells containing inflow
        NQSIJ   = NQSIJ + 1

        II = II + 1
        LLSave(II) = LL
        BCPS(II).I    = III
        BCPS(II).J    = JJJ

        BCPS(II).NQSMUL    = BCPS_GL(LL).NQSMUL
        BCPS(II).NQSMF     = BCPS_GL(LL).NQSMF
        BCPS(II).NQSERQ    = BCPS_GL(LL).NQSERQ
        BCPS(II).NCSERQ(:) = BCPS_GL(LL).NCSERQ(:)
        BCPS(II).QWIDTH    = BCPS_GL(LL).QWIDTH
        BCPS(II).QFACTOR   = BCPS_GL(LL).QFACTOR
        BCPS(II).GRPID     = BCPS_GL(LL).GRPID
        
        BCPS(II).QSSE    = BCPS_GL(LL).QSSE
        do K = 1,KC
          QSS(K,II) = BCPS(II).QSSE*DZCK(K)
        enddo

        ! *** Constituent Constants
        MMAX = 3 + NDYM + NTOX + NSED + NSND
        do MS = 1,MMAX
          do K  = 1,KC
            CQS(K,II,MS) = BCPS_GL(LL).CQSE(MS)
          enddo
        enddo
        
        ! *** GET NUMBER OF GROUPS.  GROUP NUMBERS ARE NOT REQUIRED TO BE CONTINGUOUS.
        if( BCPS(II).GRPID > NGRPID ) NGRPID = BCPS(II).GRPID
      endif
    endif
  enddo
  call DSI_All_Reduce(NGRPID, M, MPI_Max, TTDS, 1, TWAIT)
  NGRPID = M
  
  !! *** Constituent Series
  !DO L = 1,NQSIJ
  !  do N = 1,NTOX
  !    M = MSVTOX(N)
  !    BCPS(L).NCSERQ(M) = NTOXSRQ(L)
  !  enddo
  !  do N = 1,NSED
  !    M = MSVSED(N)
  !    BCPS(L).NCSERQ(M) = NSEDSRQ(L)
  !  enddo
  !  do N = 1,NSND
  !    M = MSVSND(N)
  !    BCPS(L).NCSERQ(M) = NSNDSRQ(L)
  !  enddo
  !END DO

  call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' FLOW BOUNDARIES '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3','N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3'
  do II = 1,NQSIJ
    LL = LLSave(II)
    write(mpi_mapping_unit,'(2(I5,2I8,F8.1,F8.4,F8.1,I8,1X))') LL, BCPS_GL(LL).I, BCPS_GL(LL).J, BCPS_GL(LL).QSSE, BCPS_GL(LL).QFACTOR, BCPS_GL(LL).CQSE(3), BCPS_GL(LL).NCSERQ(3),  &
                                                               II, BCPS(II).I,    BCPS(II).J,    BCPS(II).QSSE,    BCPS(II).QFACTOR,    CQS(KC,II,3),        BCPS(II).NCSERQ(3)
  enddo

  call WriteBreak(mpi_mapping_unit)

  return

End Subroutine Map_River

!---------------------------------------------------------------------------!
! @details Maps hydraulic structure boundary condition values to local ones
! @date 9/4/2019
!---------------------------------------------------------------------------!

Subroutine Map_Hydraulic_Structures

  use GLOBAL
  use HYDSTRUCMOD
  use Variables_MPI
  use Variables_MPI_Mapping

  implicit none

  ! *** Local
  integer :: NCTL, LL, MS, K, III, JJJ, L, M
  integer :: MMAX, MMIN
  integer :: NQCTL_GL
  integer :: LLSave(NQCTL)

  NQCTL_GL = NQCTL
  
  allocate(HYD_STR(0:NQCTLM))
  allocate(HSCTL(0:NQCTLM))
  do NCTL = 0,NQCTL
    HYD_STR(NCTL).IQCTLU = 0
    HYD_STR(NCTL).JQCTLU = 0
    HYD_STR(NCTL).IQCTLD = 0
    HYD_STR(NCTL).JQCTLD = 0
    HYD_STR(NCTL).NQCTYP = 0
    HYD_STR(NCTL).CURHEI = 0.0
    HYD_STR(NCTL).CURWID = 0.0
    HYD_STR(NCTL).CURSIL = 0.0
    
    HSCTL(NCTL).IREFUP = 0
    HSCTL(NCTL).JREFUP = 0
    HSCTL(NCTL).IREFDN = 0
    HSCTL(NCTL).JREFDN = 0
    HSCTL(NCTL).LUR    = 0
    HSCTL(NCTL).LDR    = 0
  enddo
  LLSave = 0

  ! *** Determine the local NQCTL value
  NQCTL = 0
  NCTL = 0
  do LL = 1,NQCTL_GL
    III = HYD_STR_GL(LL).IQCTLU
    JJJ = HYD_STR_GL(LL).JQCTLU
    if( LIJ_Global(III,JJJ) < 2 )then
      write(6,*) 'ERROR! UPSTREAM CELL IS NOT VALID FOR HYDRAULIC STRUCTURE ',LL
      call STOPP('', 1)
    endif

    III = IG2IL(HYD_STR_GL(LL).IQCTLU)
    JJJ = JG2JL(HYD_STR_GL(LL).JQCTLU)

    if( III > 0 .and. III <= IC )then     ! *** Allow ghost cells containing inflow
      if( JJJ > 0 .and. JJJ <= JC )then   ! *** Allow ghost cells containing inflow
        NQCTL = NQCTL + 1
        NCTL = NCTL + 1
        
        LLSave(NCTL) = LL
        HYD_STR(NCTL) = HYD_STR_GL(LL)
        HSCTL(NCTL)   = HSCTL_GL(LL)

        HYD_STR(NCTL).IQCTLU = III
        HYD_STR(NCTL).JQCTLU = JJJ
        HSCTL(NCTL).IREFUP   = IG2IL(HSCTL_GL(LL).IREFUP)   ! *** Assume reference cell is in the save subdomain
        HSCTL(NCTL).JREFUP   = JG2JL(HSCTL_GL(LL).JREFUP)
        if( HSCTL(NCTL).IREFUP > 0 .and. HSCTL(NCTL).IREFUP > 0 ) HSCTL(NCTL).LUR = LIJ(HSCTL(NCTL).IREFUP,HSCTL(NCTL).JREFUP)
        
        ! *** CONVERT DOWNSTREAM
        if( HYD_STR_GL(LL).IQCTLD > 0 .and. HYD_STR_GL(LL).IQCTLD > 0 )then
          III = IG2IL(HYD_STR_GL(LL).IQCTLD)
          JJJ = JG2JL(HYD_STR_GL(LL).JQCTLD)
          if( III == 0 .or. JJJ == 0 )then
            write(6,*) 'ERROR! DOWNSTREAM CELL IS NOT IN DOMAIN FOR HYDRAULIC STRUCTURE ',LL
            call STOPP('', 1)
          endif
          if( LIJ(III,JJJ) < 2 )then
            write(6,*) 'ERROR! DOWNSTREAM CELL IS NOT VALID FOR HYDRAULIC STRUCTURE ',LL
            call STOPP('', 1)
          endif
          
          HYD_STR(NCTL).IQCTLD = III
          HYD_STR(NCTL).JQCTLD = JJJ
          HSCTL(NCTL).IREFDN   = IG2IL(HSCTL_GL(LL).IREFDN)   ! *** Assume reference cell is in the save subdomain
          HSCTL(NCTL).JREFDN   = JG2JL(HSCTL_GL(LL).JREFDN)
          HSCTL(NCTL).LUR = LIJ(HSCTL(NCTL).IREFDN,HSCTL(NCTL).JREFDN)
        endif
      endif
    endif
  enddo

  call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' HYDR STRUCTURES '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IU','JU','ID','JD','NQCTYP','TABLEID','N','IU','JU','ID','JD','NQCTYP','TABLEID'
  do NCTL = 1,NQCTL
    LL = LLSave(NCTL)
    write(mpi_mapping_unit,'(2(I5,6I8,1X))')                                              &
          LL,   HYD_STR_GL(LL).IQCTLU, HYD_STR_GL(LL).JQCTLU, HYD_STR_GL(LL).IQCTLD, HYD_STR_GL(LL).JQCTLD, HYD_STR_GL(LL).NQCTYP, HYD_STR_GL(LL).NQCTLQ,  &
          NCTL, HYD_STR(NCTL).IQCTLU,  HYD_STR(NCTL).JQCTLU,  HYD_STR(NCTL).IQCTLD,  HYD_STR(NCTL).JQCTLD,  HYD_STR(NCTL).NQCTYP,  HYD_STR(NCTL).NQCTLQ
  enddo

  call WriteBreak(mpi_mapping_unit)

  return

End Subroutine Map_Hydraulic_Structures

!---------------------------------------------------------------------------!
! @details Maps withdrawal-return boundary condition values to local ones
! @date 9/4/2019
!---------------------------------------------------------------------------!

Subroutine Map_Withdrawal_Return

  use GLOBAL
  use HYDSTRUCMOD

  use Variables_MPI
  use Variables_MPI_Mapping

  implicit none

  ! *** Local
  integer :: NWR, LL, MS, K, III, JJJ, L, M
  integer :: MMAX, MMIN
  integer :: NQWR_GL
  integer :: LLSave(NQWR)
  Real, Allocatable :: CQWR_GL(:,:)
  
  allocate(CQWR_GL(NQWRM,NSTVM))
  CQWR_GL = 0.0
  
  NQWR_GL = NQWR
  
  allocate(WITH_RET(NQWR))
  allocate(WITH_RET_CTL(NQWR))
  do NWR = 1,NQWR
    WITH_RET(NWR).IQWRU = 0
    WITH_RET(NWR).JQWRU = 0
    WITH_RET(NWR).KQWRU = 0
    WITH_RET(NWR).IQWRD = 0
    WITH_RET(NWR).JQWRD = 0
    WITH_RET(NWR).KQWRD = 0
    
    WITH_RET_CTL(NWR).IREFUP = 0
    WITH_RET_CTL(NWR).JREFUP = 0
    WITH_RET_CTL(NWR).IREFDN = 0
    WITH_RET_CTL(NWR).JREFDN = 0
    
    do MS = 1,NSTVM
      CQWR_GL(NWR,MS) = CQWR(NWR,MS)
      CQWR(NWR,MS) = 0.
    enddo
  enddo

  LLSave = 0

  ! *** Determine the local NQWR value
  NQWR = 0
  NWR = 0
  do LL = 1,NQWR_GL
    III = WITH_RET_GL(LL).IQWRU
    JJJ = WITH_RET_GL(LL).JQWRU
    if( LIJ_Global(III,JJJ) < 2 )then
      write(6,*) 'ERROR! UPSTREAM CELL IS NOT VALID FOR WITHDRAWAL-RETURN BOUNDARY ',LL
      call STOPP('', 1)
    endif
    III = IG2IL(WITH_RET_GL(LL).IQWRU)
    JJJ = JG2JL(WITH_RET_GL(LL).JQWRU)
    
    if( III > 0 .and. III <= IC )then     ! *** Allow ghost cells containing inflow
      if( JJJ > 0 .and. JJJ <= JC )then   ! *** Allow ghost cells containing inflow
        NQWR = NQWR + 1
        NWR = NWR + 1
        
        LLSave(NWR) = LL
        WITH_RET(NWR) = WITH_RET_GL(LL)
        WITH_RET_CTL(NWR)    = WITH_RET_CTL_GL(LL)

        WITH_RET(NWR).IQWRU = III
        WITH_RET(NWR).JQWRU = JJJ
        WITH_RET_CTL(NWR).IREFUP   = IG2IL(WITH_RET_CTL_GL(LL).IREFUP)       ! *** Assume reference cell is in the same subdomain
        WITH_RET_CTL(NWR).JREFUP   = JG2JL(WITH_RET_CTL_GL(LL).JREFUP)

        ! *** CONVERT DOWNSTREAM
        if( WITH_RET_GL(LL).IQWRD > 0 .and. WITH_RET_GL(LL).IQWRD > 0 )then
          III = IG2IL(WITH_RET_GL(LL).IQWRD)
          JJJ = JG2JL(WITH_RET_GL(LL).JQWRD)
          if( III == 0 .or. JJJ == 0 )then
            write(6,*) 'ERROR!  FLOW RETURN CELLS NOT IN SAME DOMAIN, CHECK CARD 33 and DECOMP.INP ',LL
            call STOPP('', 1)
          endif
          WITH_RET(NWR).IQWRD = III
          WITH_RET(NWR).JQWRD = JJJ
          WITH_RET_CTL(NWR).IREFDN   = IG2IL(WITH_RET_CTL_GL(LL).IREFDN)     ! *** Assume reference cell is in the same subdomain
          WITH_RET_CTL(NWR).JREFDN   = JG2JL(WITH_RET_CTL_GL(LL).JREFDN)
        endif
        
        ! *** Convert constant rise/fall concentrations
        do MS = 1,NSTVM
          CQWR(NWR,MS) = CQWR_GL(LL,MS)
        enddo
        
      endif
    endif
  enddo

  call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' WR/RET BOUNDRY '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IU','JU','KU','ID','JD','TABLEID','N','IU','JU','KU','ID','JD','TABLEID'
  do NWR = 1,NQWR
    LL = LLSave(NWR)
    write(mpi_mapping_unit,'(2(I5,6I8,1X))')                                              &
          LL,  WITH_RET_GL(LL).IQWRU, WITH_RET_GL(LL).JQWRU, WITH_RET_GL(LL).KQWRU, WITH_RET_GL(LL).IQWRD, WITH_RET_GL(LL).JQWRD, WITH_RET_GL(LL).NQWRSERQ,  &
          NWR, WITH_RET(NWR).IQWRU,   WITH_RET(NWR).JQWRU,   WITH_RET(NWR).KQWRU,   WITH_RET(NWR).IQWRD,   WITH_RET(NWR).JQWRD,   WITH_RET(NWR).NQWRSERQ
  enddo

  call WriteBreak(mpi_mapping_unit)

  return

End Subroutine Map_Withdrawal_Return

!---------------------------------------------------------------------------!
! @details Maps jet-plume boundary condition values to local ones
!---------------------------------------------------------------------------!

Subroutine Map_Jet_Plume

  use GLOBAL

  use Variables_MPI
  use Variables_MPI_Mapping

  implicit none

  ! *** Local
  ! *** LOCAL VARIABLES
  integer :: NQJPIJ_Global
  integer :: III, JJJ, LL, IU, JU, LU, LOCAL, NJP
  integer :: LLSave(NQJPIJ)

  ! *** Save the global value
  NQJPIJ_Global = NQJPIJ
  allocate(JET_PLM(NQJPIJ))
  !DO NJP = 1,NJPSM
  !  call AllocateDSI( JET_PLM_GL(NJP).NCSERJP, -NSTVM2,    0) 
  !  call AllocateDSI( JET_PLM_GL(NJP).CWRCJP,  -NSTVM2,  0.0)
  !  call AllocateDSI( JET_PLM_GL(NJP).CQCJP,    KCM, -NSTVM2, 0.0)
  !ENDIF
  
  ! *** Initialize to zero since now we will be determining each processes local value
  NQJPIJ = 0
  LOCAL = 0
  LLSave = 0
  do LL = 1, NQJPIJ_Global
    ! *** Map to global I/J
    if( LIJ_Global(JET_PLM_GL(LL).IQJP,JET_PLM_GL(LL).JQJP) < 2 )then
      write(6,*) 'ERROR! DIFFUSER CELL IS NOT VALID FOR JET-PLUME BOUNDARY ',LL
      call STOPP('', 1)
    endif
    III = IG2IL(JET_PLM_GL(LL).IQJP)
    JJJ = JG2JL(JET_PLM_GL(LL).JQJP)

    ! *** Map to local values, including ghost cells
    if( III > 0 .and. III <= IC )then
      if( JJJ > 0 .and. JJJ <= JC )then
        ! *** Found a jet-plume cell
        NQJPIJ = NQJPIJ + 1
        LLSave(NQJPIJ) = LL
        
        JET_PLM(NQJPIJ) = JET_PLM_GL(LL)    ! *** Copy the entire structure
        
        ! *** Map to local I/J
        JET_PLM(NQJPIJ).IQJP = III
        JET_PLM(NQJPIJ).JQJP = JJJ

        if( JET_PLM(NQJPIJ).ICALJP == 2 )then
          if( LIJ_Global(JET_PLM_GL(LL).IUPCJP,JET_PLM_GL(LL).JUPCJP) < 2 )then
            write(6,*) 'ERROR! DIFFUSER WITHDRAWAL CELL IS NOT VALID FOR JET-PLUME BOUNDARY ',LL
            call STOPP('', 1)
          endif
          III = IG2IL(JET_PLM_GL(LL).IUPCJP)
          JJJ = JG2JL(JET_PLM_GL(LL).JUPCJP)
          if( LIJ(III,JJJ) < 2 )then
            write(6,*) 'ERROR! DIFFUSER WITHDRAWAL CELL IS NOT INSIDE DIFFUSER MPI DOMAIN ',LL
            call STOPP('', 1)
          endif
          
          ! *** Map to local I/J
          JET_PLM(NQJPIJ).IUPCJP = III
          JET_PLM(NQJPIJ).JUPCJP = JJJ
        endif
      endif
    endif
  enddo

  call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' JET-PLUME BNDY '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','ID','JD','KD','IU','JU','ICALJP','N','ID','JD','KD','IU','JU','ICALJP'
  do NJP = 1,NQJPIJ
    LL = LLSave(NJP)
    write(mpi_mapping_unit,'(2(I5,6I8,1X))')                                              &
          LL,  JET_PLM_GL(LL).IQJP, JET_PLM_GL(LL).JQJP, JET_PLM_GL(LL).KQJP, JET_PLM_GL(LL).IUPCJP, JET_PLM_GL(LL).JUPCJP, JET_PLM_GL(LL).ICALJP,  &
          NJP, JET_PLM(NJP).IQJP,   JET_PLM(NJP).JQJP,   JET_PLM(NJP).KQJP,   JET_PLM(NJP).IUPCJP,   JET_PLM(NJP).JUPCJP,   JET_PLM(NJP).ICALJP
  enddo

  call WriteBreak(mpi_mapping_unit)

  return

End Subroutine Map_Jet_Plume

!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! @details Maps water quality inflow boundary condition values to local ones
! @date 2020-03-13
!---------------------------------------------------------------------------!

Subroutine Map_WQ_PointSource

  use GLOBAL
  use Variables_WQ   !, only:IWQPSC, IWQPSV
  use Variables_MPI
  use Variables_MPI_Mapping

  implicit none

  ! *** Local
  integer :: II, LL, MS, K, III, JJJ, L, M, NT, NW
  integer :: MMAX, MMIN
  integer :: NWQPS_GL
  integer :: LPSL_GL (NWQPSM)
  integer :: ICPSL_GL(NWQPSM)
  integer :: JCPSL_GL(NWQPSM)
  integer :: KCPSL_GL(NWQPSM)
  integer :: MVPSL_GL(NWQPSM)
  real    :: XPSL_GL(0:NWQPSM,NWQVM)
  integer :: LLSave(NWQPSM)

  NWQPS_GL = NWQPS
  KCPSL_GL = KCPSL
  ICPSL_GL = ICPSL
  JCPSL_GL = JCPSL
  MVPSL_GL = MVPSL
  XPSL_GL  = WQWPSLC
  LLSave = 0
  
  ! *** Determine the local NWQPS value
  NWQPS = 0
  KCPSL = 0
  MVPSL = 0
  WQWPSLC = 0.
  II = 0
  do LL = 1, NWQPS_GL
    LPSL_GL(LL) = LIJ_Global(ICPSL_GL(LL),JCPSL_GL(LL))
    
    III = Map2Local(LPSL_GL(LL)).IL           ! *** Assumes WQ point source cells use the same order and number
    JJJ = Map2Local(LPSL_GL(LL)).JL           ! *** Assumes WQ point source cells use the same order and number
    if( III > 0 .and. III <= IC )then         ! *** Allow ghost cells containing inflow
      if( JJJ > 0 .and. JJJ <= JC )then       ! *** Allow ghost cells containing inflow
        NWQPS   = NWQPS + 1

        II = II + 1
        LLSave(II) = LL

        ICPSL(II) = III
        JCPSL(II) = JJJ
        KCPSL(II) = KCPSL_GL(LL)
        MVPSL(II) = MVPSL_GL(LL)
      
        ! *** ASSIGN GLOBAL CONCENTRATION TIME SERIES INDEX
        if( IWQPSL == 2 )then
          ! *** CONSTAND CONCENTRATIONS
          BCPS(II).NCSERQ(8) = NSERWQ(LL)          ! *** All WQ variables use same time series
          do NW = 1,NWQV
            if( ISTRWQ(NW) > 0 )then
              NT = MSVWQV(NW)
              do K  = 1,KC
                CQS(K,II,NT) = XPSL_GL(LL,NW)
              enddo
            endif
          enddo
        else
          ! *** Constant mass fluxes
          ! *** NW 1 TO 19 are already in mass (g/day) from WQ3DCONTROL [CONC (mg/l) and Q (m3/s)]
          ! *** NW = 20 already in moles
          WQWPSLC(II,1:NWQV) = XPSL_GL(LL,1:NWQV)
          
          ! *** Convert from mpn/l to mpn/day
          WQWPSLC(II,IFCB) = XPSL_GL(LL,IFCB) * 1000.
      
          ! *** Assign loading by layer
          L = LIJ(III,JJJ)
          if( IWQPSL == 0 ) FORALL(K = 1:KC) WQWPSL(L,K,1:NWQV) = WQWPSL(L,K,1:NWQV) + WQWPSLC(LL,1:NWQV)/DZI

          ! *** Uniform load in horizontal cell stack over all layers
          do K = 1,KC
            IWQPSC(L,K) = II
            IWQPSV(L,K) = MVPSL(II)
          enddo
        endif

      endif
    endif
  enddo

  call WriteBreak(mpi_mapping_unit)
  
  write(mpi_mapping_unit,'( 46("*"),A17,46("*") )' )  ' WQ POINT SRCS  '
  write(mpi_mapping_unit,'(2( " ****",2("********"),a8,3("********"),"|") )' ) 'GLOBAL ', 'LOCAL '
  write(mpi_mapping_unit,'(2(a5,6a8,1x))') 'N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3','N','IQS','JQS','QSSE','QFACTOR','CQS3','NCSR3'
  do II = 1,NWQPS
    LL = LLSave(II)
    write(mpi_mapping_unit,'(2(I5,2I8,F8.1,F8.4,F8.1,I8,1X))') LL, BCPS_GL(LL).I, BCPS_GL(LL).J, BCPS_GL(LL).QSSE, BCPS_GL(LL).QFACTOR, BCPS_GL(LL).CQSE(3), BCPS_GL(LL).NCSERQ(3),  &
                                                               II, BCPS(II).I,    BCPS(II).J,    BCPS(II).QSSE,    BCPS(II).QFACTOR,    CQS(KC,II,3),        BCPS(II).NCSERQ(3)
  enddo

  call WriteBreak(mpi_mapping_unit)

  return

End Subroutine Map_WQ_PointSource

  
