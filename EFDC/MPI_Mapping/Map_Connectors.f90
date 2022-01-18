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
!---------------------------------------------------------------------------!
! @details Maps east/west & north/south connectors to local subdomain
!          Each connection must be within the same subdomain
! @date 2020-04-19
!---------------------------------------------------------------------------!

Subroutine Map_Connectors

  USE GLOBAL
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  Use Broadcast_Routines
  USE INFOMOD,ONLY:SKIPCOM,READSTR

  Implicit None
  
  ! *** Local
  Integer :: ISO, LL, MS, K, L, M, NC, NP
  Integer :: MMAX, MMIN
  Integer :: NQWR_GL
  Integer, Allocatable,Dimension(:)   :: LLSave
  Integer, Allocatable,Dimension(:,:) :: III, JJJ

  Integer :: NPEWBP_Global !< Global number of cells with East/West connections
  Integer :: NPNSBP_Global !< Global Number of cells with North/South connections
  
  Character :: STR*200
  
  ! *******************************************************************************************************************
  ! *** IF ISCONNECT GE 2, READ IN EAST-WEST BOUNDARY CELLS FROM
  ! *** FILE MAPPGEW.INP TO SPECIFY EAST-WEST DIRECTION CELL CONNECTIONS
  IF( ISCONNECT >= 2 )THEN
    Allocate(I2D_Global(LCM_Global,4))
    I2D_Global = 0
    
    if( process_id == master_id )THEN

      WRITE(*,'(A)')'READING MAPPGEW.INP'
      OPEN(1,FILE='mappgew.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      
      READ(1,*,IOSTAT=ISO) NPEWBP_Global
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGEW.INP')
      
      DO NP=1,NPEWBP_Global
        !READ(1,*,IOSTAT=ISO) IWPEW_Global(NP), JWPEW_Global(NP), IEPEW_Global(NP), JEPEW_Global(NP)
        READ(1,*,IOSTAT=ISO) (I2D_Global(NP,K),K=1,4)
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGEW.INP')
      ENDDO
      CLOSE(1)
    endif

    ! *** Broadcast the globally read in values.
    Call Broadcast_Scalar(NPEWBP_Global, master_id)
    Call Broadcast_Array(I2D_Global,     master_id)

    ! *** Now map global connectors to local domain
    Allocate(III(NPEWBP_Global,2), JJJ(NPEWBP_Global,2), LLSave(NPEWBP_Global))
    III = 0
    JJJ = 0
    NC = 0
    LLSave = 0
    NPEWBP = 0
    
    ! *** Loop over all of the connectors
    DO NP = 1,NPEWBP_Global
      M = 0                                         ! *** Valid cell count.  M must = 2 or 0 to be valid
      III(NP,1) = IG2IL(I2D_Global(NP,1))           ! *** IWPEW
      JJJ(NP,1) = JG2JL(I2D_Global(NP,2))           ! *** JWPEW
      IF( III(NP,1) > 0 .AND. III(NP,1) <= IC )THEN
        IF( JJJ(NP,1) > 0 .AND. JJJ(NP,1) <= JC )THEN
          M = M + 1
        ENDIF
      ENDIF
      
      III(NP,2) = IG2IL(I2D_Global(NP,3))           ! *** IEPEW
      JJJ(NP,2) = JG2JL(I2D_Global(NP,4))           ! *** JEPEW
      IF( III(NP,2) > 0 .AND. III(NP,2) <= IC )THEN
        IF( JJJ(NP,2) > 0 .AND. JJJ(NP,2) <= JC )THEN
          M = M + 1
        ENDIF
      ENDIF

      IF( M == 1 )THEN
        WRITE(6,*) 'ERROR!  E/W CELL CONNECTORS MUST BE IN THE SAME DOMAIN', NP
        CALL STOPP('.')
      ELSEIF( M == 2 )THEN
        ! *** VALID CONNECTOR
        NPEWBP = NPEWBP + 1
        
        NC = NC + 1
        LLSave(NC) = NP
      ENDIF
    ENDDO

    ! *** Now populate the local arrays
    IF( NPEWBP > 0 )THEN
      ALLOCATE(IEPEW(NPEWBP))
      ALLOCATE(IWPEW(NPEWBP))
      ALLOCATE(JEPEW(NPEWBP))
      ALLOCATE(JWPEW(NPEWBP))
    
      DO NC = 1,NPEWBP
        NP = LLSAVE(NC)
        IWPEW(NC) = III(NP,1)
        JWPEW(NC) = JJJ(NP,1)
        IEPEW(NC) = III(NP,2)
        JEPEW(NC) = JJJ(NP,2)
      ENDDO
    ENDIF
    
    ! *** Generate the mapping report
    Call WriteBreak(mpi_mapping_unit)
  
    write(mpi_mapping_unit,'( 30("*"),A24,30("*") )' )  ' EAST-WEST CONNECTORS '
    write(mpi_mapping_unit,'(2( " ****",1("********"),a8,2("********")," |") )' ) 'GLOBAL ', 'LOCAL '
    write(mpi_mapping_unit,'(2(a5,4a8,2x))') 'N','IW','JW','IE','JE','N','IW','JW','IE','JE'
    DO NP = 1,NPEWBP
      LL = LLSave(NP)
      write(mpi_mapping_unit,'(2(I5,4I8,2X))') LL, (I2D_Global(LL,M),M=1,4),  NC, IWPEW(NP), JWPEW(NP), IEPEW(NP), JEPEW(NP)
    ENDDO

    Call WriteBreak(mpi_mapping_unit)

    DEALLOCATE(I2D_Global)
    DEALLOCATE(III, JJJ, LLSave)
    
  ENDIF
  
  ! *******************************************************************************************************************
  ! *** IF ISCONNECT GE 1, READ IN NORTH-SOUTH BOUNDARY CELLS FROM
  ! *** FILE MAPPGNS.INP TO SPECIFY NORTH-SOUTH DIRECTION CELL CONNECTIONS
  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    Allocate(I2D_Global(LCM_Global,4))
    I2D_Global = 0
    
    if( process_id == master_id )THEN

      WRITE(*,'(A)')'READING MAPPGNS.INP'
      OPEN(1,FILE='mappgns.inp',STATUS='UNKNOWN')

      ! *** SKIP OVER TITLE AND AND HEADER LINES
      STR=READSTR(1)
      
      READ(1,*,IOSTAT=ISO) NPNSBP_Global
      IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGEW.INP')
      
      DO NP=1,NPNSBP_Global
        !READ(1,*,IOSTAT=ISO) ISPNS_Global(NP), JSPNS_Global(NP), INPNS_Global(NP), JNPNS_Global(NP)
        READ(1,*,IOSTAT=ISO) (I2D_Global(NP,K),K=1,4)
        IF( ISO > 0 ) CALL STOPP('READ ERROR FOR FILE MAPPGEW.INP')
      ENDDO
      CLOSE(1)
    endif

    ! *** Broadcast the globally read in values.
    Call Broadcast_Scalar(NPNSBP_Global, master_id)
    Call Broadcast_Array(I2D_Global,     master_id)

    ! *** Now map global connectors to local domain
    Allocate(III(NPNSBP_Global,2), JJJ(NPNSBP_Global,2), LLSave(NPNSBP_Global))
    III = 0
    JJJ = 0
    NC = 0
    LLSave = 0
    NPNSBP = 0
    
    ! *** Loop over all of the connectors
    DO NP = 1,NPNSBP_Global
      M = 0                                         ! *** Valid cell count.  M must = 2 or 0 to be valid
      III(NP,1) = IG2IL(I2D_Global(NP,1))           ! *** ISPNS
      JJJ(NP,1) = JG2JL(I2D_Global(NP,2))           ! *** JSPNS
      IF( III(NP,1) > 0 .AND. III(NP,1) <= IC )THEN
        IF( JJJ(NP,1) > 0 .AND. JJJ(NP,1) <= JC )THEN
          M = M + 1
        ENDIF
      ENDIF
      
      III(NP,2) = IG2IL(I2D_Global(NP,3))           ! *** INPNS
      JJJ(NP,2) = JG2JL(I2D_Global(NP,4))           ! *** JNPNS
      IF( III(NP,2) > 0 .AND. III(NP,2) <= IC )THEN
        IF( JJJ(NP,2) > 0 .AND. JJJ(NP,2) <= JC )THEN
          M = M + 1
        ENDIF
      ENDIF

      IF( M == 1 )THEN
        WRITE(6,*) 'ERROR!  E/W CELL CONNECTORS MUST BE IN THE SAME DOMAIN', NP
        CALL STOPP('.')
      ELSEIF( M == 2 )THEN
        ! *** VALID CONNECTOR
        NPNSBP = NPNSBP + 1
        
        NC = NC + 1
        LLSave(NC) = NP
      ENDIF
    ENDDO

    ! *** Now populate the local arrays
    IF( NPNSBP > 0 )THEN
      ALLOCATE(INPNS(NPNSBP))
      ALLOCATE(ISPNS(NPNSBP))
      ALLOCATE(JNPNS(NPNSBP))
      ALLOCATE(JSPNS(NPNSBP))
    
      DO NC = 1,NPNSBP
        NP = LLSAVE(NC)
        ISPNS(NC) = III(NP,1)
        JSPNS(NC) = JJJ(NP,1)
        INPNS(NC) = III(NP,2)
        JNPNS(NC) = JJJ(NP,2)
      ENDDO
    ENDIF
    
    Call WriteBreak(mpi_mapping_unit)
  
    write(mpi_mapping_unit,'( 30("*"),A24,30("*") )' )  ' NORTH-SOUTH CONNECTORS '
    write(mpi_mapping_unit,'(2( " ****",1("********"),a8,2("********")," |") )' ) 'GLOBAL ', 'LOCAL '
    write(mpi_mapping_unit,'(2(a5,4a8,2x))') 'N','IW','JW','IE','JE','N','IW','JW','IE','JE'
    DO NP = 1,NPNSBP
      LL = LLSave(NP)
      write(mpi_mapping_unit,'(2(I5,4I8,2X))') LL, (I2D_Global(LL,M),M=1,4),  NC, ISPNS(NP), JSPNS(NP), INPNS(NP), JNPNS(NP)
    ENDDO

    Call WriteBreak(mpi_mapping_unit)

    DEALLOCATE(I2D_Global)
    DEALLOCATE(III, JJJ, LLSave)
  ENDIF

  RETURN

End Subroutine Map_Connectors

