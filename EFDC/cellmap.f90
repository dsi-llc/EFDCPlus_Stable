! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Generates global and decomposed domain cell mappings
! @date Modified to support MPI domain decomp 8/20/2019
! @author Zander Mausolff - MPI modifications

SUBROUTINE CELLMAP

  USE GLOBAL
  Use Allocate_Initialize      
  
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Broadcast_Routines
  Use Mod_Write_Cell_Map      ! Contains subroutines that write out to MPI status file
  Use Variables_MPI_Write_Out

  IMPLICIT NONE

  ! *** Local variables
  INTEGER :: I, II, J, NPN, L, LG, LE, LW, LS, LN, LF, LL, ND, NW, IGTMP, JGTMP, ICOMP, JCOMP, LTMP
  Integer :: ijctlt_wasp_id             !< Modified criteria for cell type for linkage with WASP

  ! *** New for MPI
  Integer :: iq
  Logical :: write_out
  Integer :: offset_east
  Integer :: offset_west
  Integer :: i_tmp, j_tmp, l_tmp
  
  Call AllocateDSI(IL_Global, LC_GLOBAL + 2, 0)
  Call AllocateDSI(JL_Global, LC_GLOBAL + 2, 0)
  Call AllocateDSI(IL_GL, LC_GLOBAL + 2, 0)
  Call AllocateDSI(JL_GL, LC_GLOBAL + 2, 0)

  !write(mpi_log_unit, '(a)') 'Starting cellmap'

  ! *** Decides if we are going to write to MPI status file
  write_out = .TRUE.

  ! *** Note, Removed child grid call from this subroutine
  IJCT(:,:) = 0
  IJCTLT(:,:) = 0
  DO J = 1,JC
    DO I = 1,IC
      if( IL2IG(i) <= 0  )THEN
        IJCT(I,J)   = 0
        IJCTLT(I,J) = 0
      elseif( JL2JG(j) <= 0  )THEN
        IJCT(I,J)   = 0
        IJCTLT(I,J) = 0
      else
        IJCT(I,J)   = IJCT_GLOBAL(IL2IG(I),  JL2JG(J))
        IJCTLT(I,J) = IJCTLT_GLOBAL(IL2IG(I),JL2JG(J))
      endif
    END DO
  END DO

  LIJ(:,:) = 0
  LIJ_GLOBAL(:,:) = 0
  
  ! *** Set the global LA value
  ! *** This becomes necessary for the reading in of Corners.inp

  L = 1
  DO J = 1,JC_global
    DO I = 1,IC_global
      IF( IJCT_GLOBAL(I,J) > 0 .AND. IJCT_GLOBAL(I,J) < 9 )THEN
        L = L + 1
        LIJ_Global(I,J) = L
        IL_GL(L)   = I
        JL_GL(L)   = J
      ENDIF
    ENDDO
  ENDDO
  LA_Global = L

  ! *** SET 1D CELL INDEX SEQUENCE AND MAPPINGS
  ! *** When using domain decomposition this routine gives the local LA value
  L = 1
  DO J=1,JC
    DO I=1,IC
      IF( IJCT(I,J) > 0 .AND. IJCT(I,J) < 9 )THEN
        L = L + 1
        IL_Global(L) = IL2IG(i)
        JL_Global(L) = JL2JG(j)
        IL(L)    = I
        JL(L)    = J

        LCT(L)   = IJCT(I,J)
        LIJ(I,J) = L

        ! *** Define mapping from local L to Global grid values
        Map2Global(L).PR = process_id
        Map2Global(L).IL = I
        Map2Global(L).JL = J
        
        Map2Global(L).LL = L
        Map2Global(L).IG = IL_Global(L)
        Map2Global(L).JG = JL_Global(L)
        Map2Global(L).LG = LIJ_GLOBAL( IL_Global(L), JL_Global(L) )
        
        LG = Map2Global(L).LG
        Map2Local(LG).PR = process_id
        Map2Local(LG).IL = I
        Map2Local(LG).JL = J
        Map2Local(LG).LL = L
        Map2Local(LG).IG = IL_Global(L)
        Map2Local(LG).JG = JL_Global(L)
        Map2Local(LG).LG = LIJ_GLOBAL( IL_Global(L), JL_Global(L) )

      ENDIF
    ENDDO
  ENDDO
  
  LA = L      ! *** Local LA
  LC = LA + 1
  write(mpi_log_unit,'(a,i10)') 'Local LA in scancell ',LA
  
  ! *** Create a local list of the active cells excluding the ghost cells
  Call Create_List_No_Ghost_Cells
  
  ! *** Create mapping of local i/j to global i/j to help with the mapping of the connectors
  ! *** Could remove?  Not used anymore.  Has a bug in it in a few test cases
  Call Setup_Local_to_Global
  
  ! *** Allocate the first and last values ofthe IL/JL arrays
  IL(1)= 0
  IL(LC)=0
  JL(1)= 0
  JL(LC)=0

  ! *** Wind Fetch Mask Assingments
  UMASK(0:1) = 1
  UMASK(LC)  = 1
  VMASK(0:1) = 1
  VMASK(LC)  = 1

  ! ****************************************************************************************
  ! *** Added MPI Domain Decomp
  IF( ISCONNECT > 0 )THEN
    Call Map_Connectors   ! *** Maps east/west and north/south connectors 
  ENDIF

  ! *** SET Map2Global TO LC FOR UNINITIALIZED CELLS
  DO J=1,JC
    DO I=1,IC
      L = LIJ(I,J)
      IF( L > 0 )THEN
        IF( Map2Global(L).LG == 0 ) Map2Global(L).LG = LA_Global + 1
      ENDIF
    ENDDO
  ENDDO
  Map2Global(1).LG  = LA_Global + 1
  Map2Global(LC).LG = LA_Global + 1

  ! *** Zeros out IL2IG and JL2JG if values are negative or beyond the domain size.
  ! This basically fixes the fact that sometimes arrays are indexed out of bounds because of the way
  ! IL2IG and JL2JG are setup:

  ! *** O'donncha's comments are in the next two lines
  ! *** Extract coordinates of computational footprint only (i.e. excluding ghost zones)
  ! *** for in particular conjugate gradient computation
  LL = 0
  DO I = 3, IC-2
    DO J = 3,JC-2
      IF(  IJCT(I,J) > 0 .AND. IJCT(I,J) < 9 )THEN
        LL = LL + 1
      END IF
    END DO
  END DO

  write(mpi_log_unit, '(a,I8)') 'Local Active Cells Excluding Ghost Cells ', LL
  write(mpi_log_unit, '(a,I8)') 'LA Local ', LA
  write(mpi_log_unit, '(a,I8)') 'IC Local ', IC
  write(mpi_log_unit, '(a,I8)') 'JC Local ', JC
  Call WriteBreak(mpi_log_unit)

  ! *** Write out to MPI status file
  IF( MPI_WRITE_FLAG == .TRUE. )THEN
    Call Write_LIJ(write_out)
  ENDIF
  
  ! *** Wasp linkage
  IF( ISWASP <= 5 )THEN
    ijctlt_wasp_id = 9 ! iwasp <=5
  End if

  IF( ISWASP >= 6 )THEN
    ijctlt_wasp_id = 7 ! iwasp >= 6
  End if

  ! *** Rewrote to remove duplication
  IF( ISWASP <= 5 .OR. ISWASP >= 6 )THEN
    L = 1
    DO J=1,JC
      DO I=1,IC
        IF( IJCTLT(I,J) > 0 .AND. IJCTLT(I,J) < ijctlt_wasp_id )THEN
          L = L + 1
          ILLT(L)=I
          JLLT(L)=J
          LCTLT(L)=IJCTLT(I,J)
          LIJLT(I,J)=L
        ENDIF
      ENDDO
    ENDDO
    LALT = L
    LCLT = L + 1
  ENDIF
  ! *** End wasp linkage

  ! *** SET NORTH AND SOUTH CELL IDENTIFIER ARRAYS
  ! Visualize the arrays being set in the 4 lines below 
  ! ---------------------
  ! | LNWC | LNC | LNEC |
  ! | LWC  | i,j | LEC  |
  ! | LSWC | LSC | LSEC |
  ! ---------------------

  LEC(1)   = LC
  LWC(1)   = LC
  LNC(1)   = LC
  LSC(1)   = LC
  LNEC(1)  = LC
  LNWC(1)  = LC
  LSEC(1)  = LC
  LSWC(1)  = LC

  LEC(LC)  = 1
  LWC(LC)  = 1
  LNC(LC)  = 1
  LSC(LC)  = 1
  LNEC(LC) = 1
  LNWC(LC) = 1
  LSEC(LC) = 1
  LSWC(LC) = 1

  ! *** BUILD CONNECTIONS IN THE 8 DIRECTIONS SURROUNDING CELL L
  DO L=2,LA
    I = IL(L) 
    J = JL(L)

    IF( IJCT(I+1,J) == 0 .OR. IJCT(I+1,J) == 9 )THEN
      LEC(L) = LC
    ELSE
      LEC(L) = MAX(LIJ(I+1,J), 1)
    ENDIF
    
    IF( IJCT(I-1,J) == 0 .OR. IJCT(I-1,J) == 9 )THEN
      LWC(L) = LC
    ELSE
      LWC(L) = MAX(LIJ(I-1,J), 1)
    ENDIF

    IF( IJCT(I,J+1) == 0 .OR. IJCT(I,J+1) == 9 )THEN
      LNC(L) = LC
    ELSE
      LNC(L) = MAX(LIJ(I,J+1), 1)
    ENDIF

    IF( IJCT(I,J-1) == 0 .OR. IJCT(I,J-1) > 8 )THEN
      LSC(L) = LC
    ELSE
      LSC(L) = MAX(LIJ(I,J-1), 1)
    ENDIF

    IF( IJCT(I+1,J+1) == 0 .OR. IJCT(I+1,J+1) == 9 )THEN
      LNEC(L) = LC
    ELSE
      LNEC(L) = MAX(LIJ(I+1,J+1),1)
    ENDIF

    IF( IJCT(I-1,J+1) == 0 .OR. IJCT(I-1,J+1) == 9 )THEN
      LNWC(L) = LC
    ELSE
      LNWC(L) = MAX(LIJ(I-1,J+1),1)
    ENDIF
    
    IF( IJCT(I+1,J-1) == 0 .OR. IJCT(I+1,J-1) == 9 )THEN
      LSEC(L) = LC
    ELSE
      LSEC(L) = MAX(LIJ(I+1,J-1),1)
    ENDIF

    IF( IJCT(I-1,J-1) == 0 .OR. IJCT(I-1,J-1) == 9 )THEN
      LSWC(L) = LC
    ELSE
      LSWC(L) = MAX(LIJ(I-1,J-1),1)
    ENDIF

  ENDDO
  
  ! **  MODIFY EAST-WEST CELL MAPPING FOR PERIOD GRID IN E-W DIRECTION
  IF( ISCONNECT >= 2 )THEN !==> both E-W and N-S connections, case where MAPPGEW.inp is being read
    offset_east = 1
    offset_west = 1
    DO NPN=1,NPEWBP !***NPEWBP is now local to a subdomain
        
        L = LIJ(IWPEW(NPN),JWPEW(NPN))
        LWC(L) = MAX(LIJ(IEPEW(NPN),JEPEW(NPN)), 1)
        
        IF( IJCT(IEPEW(NPN),JEPEW(NPN)-1) == 9 )THEN
            LSWC(L) = LC
        ELSE
            LSWC(L) = MAX(LIJ(IEPEW(NPN),JEPEW(NPN)-1), 1)
        ENDIF
        IF( IJCT(IEPEW(NPN),JEPEW(NPN)+1) == 9 )THEN
            LNWC(L) = LC
        ELSE
            LNWC(L) = MAX(LIJ(IEPEW(NPN),JEPEW(NPN)+1), 1)
        ENDIF

        L = LIJ(IEPEW(NPN),JEPEW(NPN))
        LEC(L) = MAX(LIJ(IWPEW(NPN),JWPEW(NPN)), 1)
        
        IF( IJCT(IWPEW(NPN),JWPEW(NPN)-1) == 9 )THEN
            LSEC(L) = LC
        ELSE
            LSEC(L) = MAX(LIJ(IWPEW(NPN),JWPEW(NPN)-1), 1)
        ENDIF
        IF( IJCT(IWPEW(NPN),JWPEW(NPN)+1) == 9 )THEN
            LNEC(L) = LC
        ELSE
            LNEC(L) = MAX(LIJ(IWPEW(NPN),JWPEW(NPN)+1), 1)
        ENDIF

    ENDDO !***End do over # of connectors
    
  ENDIF ! *** ISCONNECT > 2

  ! **  MODIFY NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION
  ! Case where we are reading MAPPGNS.inp
  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    DO NPN=1,NPNSBP

      LS = LIJ(ISPNS(NPN),JSPNS(NPN))
      LSC(LS) = MAX(LIJ(INPNS(NPN),JNPNS(NPN)), 1)

      IF( IJCT(INPNS(NPN)+1,JNPNS(NPN)) == 9 )THEN
        LSEC(LS) = LC
      ELSE
        LSEC(LS) = MAX(LIJ(INPNS(NPN)+1,JNPNS(NPN)), 1)
      ENDIF
      IF( IJCT(INPNS(NPN)-1,JNPNS(NPN)) == 9 )THEN
        LSWC(LS) = LC
      ELSE
        LSWC(LS) = MAX(LIJ(INPNS(NPN)-1,JNPNS(NPN)), 1)
      ENDIF

      LN = LIJ(INPNS(NPN),JNPNS(NPN))
      LNC(LN) = MAX(LIJ(ISPNS(NPN),JSPNS(NPN)), 1)

      IF( IJCT(ISPNS(NPN)+1,JSPNS(NPN)) == 9 )THEN
        LNEC(LN) = LC
      ELSE
        LNEC(LN) = MAX(LIJ(ISPNS(NPN)+1,JSPNS(NPN)), 1)
      ENDIF
      IF( IJCT(ISPNS(NPN)-1,JSPNS(NPN)) == 9 )THEN
        LNWC(LN) = LC
      ELSE
        LNWC(LN) = MAX(LIJ(ISPNS(NPN)-1,JSPNS(NPN)), 1)
      ENDIF
    ENDDO
  ENDIF ! isconnect = 1 or = 3

  ! **  SET LT NORTH AND SOUTH CELL IDENTIFIER ARRAYS
  LNCLT(1)=LCLT
  LSCLT(1)=LCLT
  LECLT(1)=LCLT
  LWCLT(1)=LCLT

  LECLT(LC)=1
  LWCLT(LC)=1
  LNCLT(LC)=1
  LSCLT(LC)=1

  DO L=2,LALT
    I=ILLT(L)
    J=JLLT(L)
    ! ** N & S
    IF( IJCTLT(I,J+1) == 9 )THEN
      LNCLT(L)=LCLT
    ELSE
      LNCLT(L)=LIJLT(I,J+1)
    ENDIF
    IF( IJCTLT(I,J-1) == 9 )THEN
      LSCLT(L)=LCLT
    ELSE
      LSCLT(L)=LIJLT(I,J-1)
    ENDIF

    ! ** W & E
    IF( IJCTLT(I+1,J) == 9 )THEN
      LECLT(L)=LCLT
    ELSE
      LECLT(L)=LIJLT(I+1,J)
    ENDIF
    IF( IJCTLT(I-1,J) == 9 )THEN
      LWCLT(L)=LCLT
    ELSE
      LWCLT(L)=LIJLT(I-1,J)
    ENDIF

  ENDDO

  ! **  MODIFY LT NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION
  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    DO NPN=1, NPNSBP
      ! SET NORTH CELL SOUTH OF SOUTH CELL
      LS=LIJLT(ISPNS(NPN),JSPNS(NPN))
      LSCLT(LS)=LIJLT(INPNS(NPN),JNPNS(NPN))

      !     SET SOUTH CELL NORTH OF NORTH CELL
      LN=LIJLT(INPNS(NPN),JNPNS(NPN))
      LNCLT(LN)=LIJLT(ISPNS(NPN),JSPNS(NPN))
    END DO
  END IF

  ! **  MODIFY 6 CELLS FOR W-E CONNECTIONS (IJCTLT)
  IF( ISCONNECT >= 2 )THEN
    DO NPN=1,NPEWBP
      LE = LIJLT(IEPEW(NPN),JEPEW(NPN))
      LECLT(LE) = LIJLT(IWPEW(NPN),JWPEW(NPN))

      LW = LIJLT(IWPEW(NPN),JWPEW(NPN))
      LWCLT(LW) = LIJLT(IEPEW(NPN),JEPEW(NPN))
    END DO
  ENDIF

  ! *** ****************************
  if( process_id == master_id )THEN
    if(write_out == .TRUE. .AND. MPI_Write_Flag == .TRUE. )THEN
      WRITE(7,1616)LALT,LCLT
      WRITE(8,1616)LALT,LCLT

      WRITE(7,601)LA
      WRITE(8,601)LA

601   FORMAT('  LA=',I10,//)
1616  FORMAT(2I10)

      WRITE(8,'(7A6)') 'I', 'J', 'L', 'LWC', 'LEC', 'LSC', 'LNC'
      DO L=2,LA
        WRITE(8,'(7I6)') IL(L), JL(L), L, LWC(L), LEC(L), LSC(L), LNC(L)
      ENDDO
    endif
  Endif

  if( active_domains == 1 )then
    ! *** Connections
    do L = 2,LA
      LWC_Global(L) = LWC(L)
      LEC_Global(L) = LEC(L)
      LSC_Global(L) = LSC(L)
      LNC_Global(L) = LNC(L)
    enddo
  endif
  
  ! *** Initialize ghost cell/active cell communication lists
  Call Communicate_Initialize

  ! *** Write out to MPI status file
  if( MPI_Write_Flag == .TRUE. )THEN
    Call Write_Cell_Indexing(write_out)
  endif

  RETURN

END Subroutine Cellmap
