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
!---------------------------------------------------------------------------!
! @details Maps east/west & north/south connectors to local subdomain
!          Each connection must be within the same subdomain
! @date 2020-04-19
!---------------------------------------------------------------------------!

Subroutine Map_Connectors

  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Mapping
  use Variables_MPI_Write_Out
  use Broadcast_Routines
  use INFOMOD,only:SKIPCOM,READSTR

  implicit none
  
  ! *** Local
  integer :: ISO, LL, MS, K, L, M, NC, NP
  integer :: MMAX, MMIN
  integer :: NQWR_GL
  integer, Allocatable,Dimension(:)   :: LLSave
  integer, Allocatable,Dimension(:,:) :: III, JJJ

  integer :: NPEWBP_Global !< Global number of cells with East/West connections
  integer :: NPNSBP_Global !< Global Number of cells with North/South connections
  
  Character :: STR*200
  
  ! *******************************************************************************************************************
  ! *** IF ISCONNECT GE 2, READ IN EAST-WEST BOUNDARY CELLS FROM
  ! *** FILE MAPPGEW.INP TO SPECIFY EAST-WEST DIRECTION CELL CONNECTIONS
  if( ISCONNECT >= 2 )then
    allocate(I2D_Global(LCM_Global,4))
    I2D_Global = 0
    
    if( process_id == master_id )then

      write(*,'(A)')'READING MAPPGEW.INP'
      open(1,FILE = 'mappgew.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      
      read(1,*,IOSTAT = ISO) NPEWBP_Global
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGEW.INP')
      
      do NP = 1,NPEWBP_Global
        !READ(1,*,IOSTAT = ISO) IWPEW_Global(NP), JWPEW_Global(NP), IEPEW_Global(NP), JEPEW_Global(NP)
        read(1,*,IOSTAT = ISO) (I2D_Global(NP,K),K = 1,4)
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGEW.INP')
      enddo
      close(1)
    endif

    ! *** Broadcast the globally read in values.
    call Broadcast_Scalar(NPEWBP_Global, master_id)
    call Broadcast_Array(I2D_Global,     master_id)

    ! *** Now map global connectors to local domain
    allocate(III(NPEWBP_Global,2), JJJ(NPEWBP_Global,2), LLSave(NPEWBP_Global))
    III = 0
    JJJ = 0
    NC = 0
    LLSave = 0
    NPEWBP = 0
    
    ! *** Loop over all of the connectors
    do NP = 1,NPEWBP_Global
      M = 0                                         ! *** Valid cell count.  M must = 2 or 0 to be valid
      III(NP,1) = IG2IL(I2D_Global(NP,1))           ! *** IWPEW
      JJJ(NP,1) = JG2JL(I2D_Global(NP,2))           ! *** JWPEW
      if( III(NP,1) > 0 .and. III(NP,1) <= IC )then
        if( JJJ(NP,1) > 0 .and. JJJ(NP,1) <= JC )then
          M = M + 1
        endif
      endif
      
      III(NP,2) = IG2IL(I2D_Global(NP,3))           ! *** IEPEW
      JJJ(NP,2) = JG2JL(I2D_Global(NP,4))           ! *** JEPEW
      if( III(NP,2) > 0 .and. III(NP,2) <= IC )then
        if( JJJ(NP,2) > 0 .and. JJJ(NP,2) <= JC )then
          M = M + 1
        endif
      endif

      if( M == 1 )then
        write(6,*) 'ERROR!  E/W CELL CONNECTORS MUST BE IN THE SAME DOMAIN', NP
        call STOPP('', 1)
      elseif( M == 2 )then
        ! *** VALID CONNECTOR
        NPEWBP = NPEWBP + 1
        
        NC = NC + 1
        LLSave(NC) = NP
      endif
    enddo

    ! *** Now populate the local arrays
    if( NPEWBP > 0 )then
      allocate(IEPEW(NPEWBP))
      allocate(IWPEW(NPEWBP))
      allocate(JEPEW(NPEWBP))
      allocate(JWPEW(NPEWBP))
    
      do NC = 1,NPEWBP
        NP = LLSAVE(NC)
        IWPEW(NC) = III(NP,1)
        JWPEW(NC) = JJJ(NP,1)
        IEPEW(NC) = III(NP,2)
        JEPEW(NC) = JJJ(NP,2)
      enddo
    endif
    
    ! *** Generate the mapping report
    call WriteBreak(mpi_mapping_unit)
  
    write(mpi_mapping_unit,'( 30("*"),A24,30("*") )' )  ' EAST-WEST CONNECTORS '
    write(mpi_mapping_unit,'(2( " ****",1("********"),a8,2("********")," |") )' ) 'GLOBAL ', 'LOCAL '
    write(mpi_mapping_unit,'(2(a5,4a8,2x))') 'N','IW','JW','IE','JE','N','IW','JW','IE','JE'
    do NP = 1,NPEWBP
      LL = LLSave(NP)
      write(mpi_mapping_unit,'(2(I5,4I8,2X))') LL, (I2D_Global(LL,M),M = 1,4),  NC, IWPEW(NP), JWPEW(NP), IEPEW(NP), JEPEW(NP)
    enddo

    call WriteBreak(mpi_mapping_unit)

    deallocate(I2D_Global)
    deallocate(III, JJJ, LLSave)
    
  endif
  
  ! *******************************************************************************************************************
  ! *** IF ISCONNECT GE 1, READ IN NORTH-SOUTH BOUNDARY CELLS FROM
  ! *** FILE MAPPGNS.INP TO SPECIFY NORTH-SOUTH DIRECTION CELL CONNECTIONS
  if( ISCONNECT == 1 .or. ISCONNECT == 3 )then
    allocate(I2D_Global(LCM_Global,4))
    I2D_Global = 0
    
    if( process_id == master_id )then

      write(*,'(A)')'READING MAPPGNS.INP'
      open(1,FILE = 'mappgns.inp',STATUS = 'UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR = READSTR(1)
      
      read(1,*,IOSTAT = ISO) NPNSBP_Global
      if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGNS.INP')
      
      do NP = 1,NPNSBP_Global
        !READ(1,*,IOSTAT = ISO) ISPNS_Global(NP), JSPNS_Global(NP), INPNS_Global(NP), JNPNS_Global(NP)
        read(1,*,IOSTAT = ISO) (I2D_Global(NP,K),K = 1,4)
        if( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGNS.INP')
      enddo
      close(1)
    endif

    ! *** Broadcast the globally read in values.
    call Broadcast_Scalar(NPNSBP_Global, master_id)
    call Broadcast_Array(I2D_Global,     master_id)

    ! *** Now map global connectors to local domain
    allocate(III(NPNSBP_Global,2), JJJ(NPNSBP_Global,2), LLSave(NPNSBP_Global))
    III = 0
    JJJ = 0
    NC = 0
    LLSave = 0
    NPNSBP = 0
    
    ! *** Loop over all of the connectors
    do NP = 1,NPNSBP_Global
      M = 0                                         ! *** Valid cell count.  M must = 2 or 0 to be valid
      III(NP,1) = IG2IL(I2D_Global(NP,1))           ! *** ISPNS
      JJJ(NP,1) = JG2JL(I2D_Global(NP,2))           ! *** JSPNS
      if( III(NP,1) > 0 .and. III(NP,1) <= IC )then
        if( JJJ(NP,1) > 0 .and. JJJ(NP,1) <= JC )then
          M = M + 1
        endif
      endif
      
      III(NP,2) = IG2IL(I2D_Global(NP,3))           ! *** INPNS
      JJJ(NP,2) = JG2JL(I2D_Global(NP,4))           ! *** JNPNS
      if( III(NP,2) > 0 .and. III(NP,2) <= IC )then
        if( JJJ(NP,2) > 0 .and. JJJ(NP,2) <= JC )then
          M = M + 1
        endif
      endif

      if( M == 1 )then
        write(6,*) 'ERROR!  E/W CELL CONNECTORS MUST BE IN THE SAME DOMAIN', NP
        call STOPP('', 1)
      elseif( M == 2 )then
        ! *** VALID CONNECTOR
        NPNSBP = NPNSBP + 1
        
        NC = NC + 1
        LLSave(NC) = NP
      endif
    enddo

    ! *** Now populate the local arrays
    if( NPNSBP > 0 )then
      allocate(INPNS(NPNSBP))
      allocate(ISPNS(NPNSBP))
      allocate(JNPNS(NPNSBP))
      allocate(JSPNS(NPNSBP))
    
      do NC = 1,NPNSBP
        NP = LLSAVE(NC)
        ISPNS(NC) = III(NP,1)
        JSPNS(NC) = JJJ(NP,1)
        INPNS(NC) = III(NP,2)
        JNPNS(NC) = JJJ(NP,2)
      enddo
    endif
    
    call WriteBreak(mpi_mapping_unit)
  
    write(mpi_mapping_unit,'( 30("*"),A24,30("*") )' )  ' NORTH-SOUTH CONNECTORS '
    write(mpi_mapping_unit,'(2( " ****",1("********"),a8,2("********")," |") )' ) 'GLOBAL ', 'LOCAL '
    write(mpi_mapping_unit,'(2(a5,4a8,2x))') 'N','IW','JW','IE','JE','N','IW','JW','IE','JE'
    do NP = 1,NPNSBP
      LL = LLSave(NP)
      write(mpi_mapping_unit,'(2(I5,4I8,2X))') LL, (I2D_Global(LL,M),M = 1,4),  NC, ISPNS(NP), JSPNS(NP), INPNS(NP), JNPNS(NP)
    enddo

    call WriteBreak(mpi_mapping_unit)

    deallocate(I2D_Global)
    deallocate(III, JJJ, LLSave)
  endif

  return

End Subroutine Map_Connectors

