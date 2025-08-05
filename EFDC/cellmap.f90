! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
! @details Generates global and decomposed domain cell mappings
! @date Modified to support MPI domain decomp 8/20/2019
! @author Zander Mausolff - MPI modifications

SUBROUTINE CELLMAP

  use GLOBAL
  use Allocate_Initialize      
  
  use Variables_MPI
  use Variables_MPI_Mapping
  use Broadcast_Routines
  use Mod_Write_Cell_Map      ! Contains subroutines that write out to MPI status file
  use Variables_MPI_Write_Out

  implicit none

  ! *** Local variables
  integer :: I, II, J, NPN, L, LG, LE, LW, LS, LN, LF, LL, ND, NW, IGTMP, JGTMP, ICOMP, JCOMP, LTMP
  integer :: ijctlt_wasp_id             !< Modified criteria for cell type for linkage with WASP

  ! *** New for MPI
  integer :: iq
  logical :: write_out
  integer :: offset_east
  integer :: offset_west
  integer :: i_tmp, j_tmp, l_tmp
  
  call AllocateDSI(IL_Global, LC_GLOBAL + 2, 0)
  call AllocateDSI(JL_Global, LC_GLOBAL + 2, 0)
  call AllocateDSI(IL_GL, LC_GLOBAL + 2, 0)
  call AllocateDSI(JL_GL, LC_GLOBAL + 2, 0)

  !write(mpi_log_unit, '(a)') 'Starting cellmap'

  ! *** Decides if we are going to write to MPI status file
  write_out = .TRUE.

  ! *** Note, Removed child grid call from this subroutine
  IJCT(:,:) = 0
  IJCTLT(:,:) = 0
  do J = 1,JC
    do I = 1,IC
      if( IL2IG(i) <= 0  )then
        IJCT(I,J)   = 0
        IJCTLT(I,J) = 0
      elseif( JL2JG(j) <= 0  )then
        IJCT(I,J)   = 0
        IJCTLT(I,J) = 0
      else
        IJCT(I,J)   = IJCT_GLOBAL(IL2IG(I),  JL2JG(J))
        IJCTLT(I,J) = IJCTLT_GLOBAL(IL2IG(I),JL2JG(J))
      endif
    enddo
  enddo

  LIJ(:,:) = 0
  LIJ_GLOBAL(:,:) = 0
  
  ! *** Set the global LA and the global I and J mapping
  LG = 1
  do J = 1,JC_global
    do I = 1,IC_global
      if( IJCT_GLOBAL(I,J) > 0 .and. IJCT_GLOBAL(I,J) < 9 )then
        LG = LG + 1
        LIJ_Global(I,J) = LG
        IL_GL(LG)   = I
        JL_GL(LG)   = J
        Map2Local(LG).IG = I
        Map2Local(LG).JG = J
      endif
    enddo
  enddo
  LA_Global = LG

  ! *** SET 1D CELL INDEX SEQUENCE AND MAPPINGS
  ! *** When using domain decomposition this routine gives the local LA value
  L = 1
  do J = 1,JC
    do I = 1,IC
      if( IJCT(I,J) > 0 .and. IJCT(I,J) < 9 )then
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

      endif
    enddo
  enddo
  
  LA = L      ! *** Local LA
  LC = LA + 1
  write(mpi_log_unit,'(a,i10)') 'Local LA in scancell ',LA
  
  ! *** Create a local list of the active cells excluding the ghost cells
  call Create_List_No_Ghost_Cells
  
  ! *** Create mapping of local i/j to global i/j to help with the mapping of the connectors
  ! *** Could remove?  Not used anymore.  Has a bug in it in a few test cases
  call Setup_Local_to_Global
  
  ! *** Allocate the first and last values ofthe IL/JL arrays
  IL(1)= 0
  IL(LC) = 0
  JL(1)= 0
  JL(LC) = 0

  ! *** Wind Fetch Mask Assingments
  UMASK(0:1) = 1
  UMASK(LC)  = 1
  VMASK(0:1) = 1
  VMASK(LC)  = 1

  ! ****************************************************************************************
  ! *** Added MPI Domain Decomp
  if( ISCONNECT > 0 )then
    call Map_Connectors   ! *** Maps east/west and north/south connectors 
  endif

  ! *** SET Map2Global TO LC FOR UNINITIALIZED CELLS
  do J = 1,JC
    do I = 1,IC
      L = LIJ(I,J)
      if( L > 0 )then
        if( Map2Global(L).LG == 0 ) Map2Global(L).LG = LA_Global + 1
      endif
    enddo
  enddo
  Map2Global(1).LG  = LA_Global + 1
  Map2Global(LC).LG = LA_Global + 1

  ! *** Zeros out IL2IG and JL2JG if values are negative or beyond the domain size.
  ! This basically fixes the fact that sometimes arrays are indexed out of bounds because of the way
  ! IL2IG and JL2JG are setup:

  ! *** O'donncha's comments are in the next two lines
  ! *** Extract coordinates of computational footprint only (i.e. excluding ghost zones)
  ! *** for in particular conjugate gradient computation
  LL = 0
  do I = 3, IC-2
    do J = 3,JC-2
      if( IJCT(I,J) > 0 .and. IJCT(I,J) < 9 )then
        LL = LL + 1
      endif
    enddo
  enddo

  write(mpi_log_unit, '(a,I8)') 'Local Active Cells Excluding Ghost Cells ', LL
  write(mpi_log_unit, '(a,I8)') 'LA Local ', LA
  write(mpi_log_unit, '(a,I8)') 'IC Local ', IC
  write(mpi_log_unit, '(a,I8)') 'JC Local ', JC
  call WriteBreak(mpi_log_unit)

  ! *** Write out to MPI status file
  if( MPI_Write_Flag )then
    call Write_LIJ(write_out)
  endif
  
  ! *** Wasp linkage
  if( ISWASP <= 5 )then
    ijctlt_wasp_id = 9 ! iwasp <=5
  endif

  if( ISWASP >= 6 )then
    ijctlt_wasp_id = 7 ! iwasp >= 6
  endif

  ! *** Rewrote to remove duplication
  if( ISWASP <= 5 .or. ISWASP >= 6 )then
    L = 1
    do J = 1,JC
      do I = 1,IC
        if( IJCTLT(I,J) > 0 .and. IJCTLT(I,J) < ijctlt_wasp_id )then
          L = L + 1
          ILLT(L) = I
          JLLT(L) = J
          LCTLT(L) = IJCTLT(I,J)
          LIJLT(I,J) = L
        endif
      enddo
    enddo
    LALT = L
    LCLT = L + 1
  endif
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
  do L = 2,LA
    I = IL(L) 
    J = JL(L)

    if( IJCT(I+1,J) == 0 .or. IJCT(I+1,J) == 9 )then
      LEC(L) = LC
    else
      LEC(L) = MAX(LIJ(I+1,J), 1)
    endif
    
    if( IJCT(I-1,J) == 0 .or. IJCT(I-1,J) == 9 )then
      LWC(L) = LC
    else
      LWC(L) = MAX(LIJ(I-1,J), 1)
    endif

    if( IJCT(I,J+1) == 0 .or. IJCT(I,J+1) == 9 )then
      LNC(L) = LC
    else
      LNC(L) = MAX(LIJ(I,J+1), 1)
    endif

    if( IJCT(I,J-1) == 0 .or. IJCT(I,J-1) > 8 )then
      LSC(L) = LC
    else
      LSC(L) = MAX(LIJ(I,J-1), 1)
    endif

    if( IJCT(I+1,J+1) == 0 .or. IJCT(I+1,J+1) == 9 )then
      LNEC(L) = LC
    else
      LNEC(L) = MAX(LIJ(I+1,J+1),1)
    endif

    if( IJCT(I-1,J+1) == 0 .or. IJCT(I-1,J+1) == 9 )then
      LNWC(L) = LC
    else
      LNWC(L) = MAX(LIJ(I-1,J+1),1)
    endif
    
    if( IJCT(I+1,J-1) == 0 .or. IJCT(I+1,J-1) == 9 )then
      LSEC(L) = LC
    else
      LSEC(L) = MAX(LIJ(I+1,J-1),1)
    endif

    if( IJCT(I-1,J-1) == 0 .or. IJCT(I-1,J-1) == 9 )then
      LSWC(L) = LC
    else
      LSWC(L) = MAX(LIJ(I-1,J-1),1)
    endif

  enddo
  
  ! ***  MODIFY EAST-WEST CELL MAPPING FOR PERIOD GRID IN E-W DIRECTION
  if( ISCONNECT >= 2 )then !==> both E-W and N-S connections, case where MAPPGEW.inp is being read
    offset_east = 1
    offset_west = 1
    do NPN = 1,NPEWBP !***NPEWBP is now local to a subdomain
        
        L = LIJ(IWPEW(NPN),JWPEW(NPN))
        LWC(L) = MAX(LIJ(IEPEW(NPN),JEPEW(NPN)), 1)
        
        if( IJCT(IEPEW(NPN),JEPEW(NPN)-1) == 9 )then
            LSWC(L) = LC
        else
            LSWC(L) = MAX(LIJ(IEPEW(NPN),JEPEW(NPN)-1), 1)
        endif
        if( IJCT(IEPEW(NPN),JEPEW(NPN)+1) == 9 )then
            LNWC(L) = LC
        else
            LNWC(L) = MAX(LIJ(IEPEW(NPN),JEPEW(NPN)+1), 1)
        endif

        L = LIJ(IEPEW(NPN),JEPEW(NPN))
        LEC(L) = MAX(LIJ(IWPEW(NPN),JWPEW(NPN)), 1)
        
        if( IJCT(IWPEW(NPN),JWPEW(NPN)-1) == 9 )then
            LSEC(L) = LC
        else
            LSEC(L) = MAX(LIJ(IWPEW(NPN),JWPEW(NPN)-1), 1)
        endif
        if( IJCT(IWPEW(NPN),JWPEW(NPN)+1) == 9 )then
            LNEC(L) = LC
        else
            LNEC(L) = MAX(LIJ(IWPEW(NPN),JWPEW(NPN)+1), 1)
        endif

    enddo !***End do over # of connectors
    
  endif ! *** ISCONNECT > 2

  ! ***  MODIFY NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION
  ! Case where we are reading MAPPGNS.inp
  if( ISCONNECT == 1 .or. ISCONNECT == 3 )then
    do NPN = 1,NPNSBP

      LS = LIJ(ISPNS(NPN),JSPNS(NPN))
      LSC(LS) = MAX(LIJ(INPNS(NPN),JNPNS(NPN)), 1)

      if( IJCT(INPNS(NPN)+1,JNPNS(NPN)) == 9 )then
        LSEC(LS) = LC
      else
        LSEC(LS) = MAX(LIJ(INPNS(NPN)+1,JNPNS(NPN)), 1)
      endif
      if( IJCT(INPNS(NPN)-1,JNPNS(NPN)) == 9 )then
        LSWC(LS) = LC
      else
        LSWC(LS) = MAX(LIJ(INPNS(NPN)-1,JNPNS(NPN)), 1)
      endif

      LN = LIJ(INPNS(NPN),JNPNS(NPN))
      LNC(LN) = MAX(LIJ(ISPNS(NPN),JSPNS(NPN)), 1)

      if( IJCT(ISPNS(NPN)+1,JSPNS(NPN)) == 9 )then
        LNEC(LN) = LC
      else
        LNEC(LN) = MAX(LIJ(ISPNS(NPN)+1,JSPNS(NPN)), 1)
      endif
      if( IJCT(ISPNS(NPN)-1,JSPNS(NPN)) == 9 )then
        LNWC(LN) = LC
      else
        LNWC(LN) = MAX(LIJ(ISPNS(NPN)-1,JSPNS(NPN)), 1)
      endif
    enddo
  endif ! isconnect = 1 or = 3

  ! *** SET WASP/LT NORTH AND SOUTH CELL IDENTIFIER ARRAYS
  LNCLT(1) = LCLT
  LSCLT(1) = LCLT
  LECLT(1) = LCLT
  LWCLT(1) = LCLT

  LECLT(LC) = 1
  LWCLT(LC) = 1
  LNCLT(LC) = 1
  LSCLT(LC) = 1

  DO L = 2,LALT
    I = ILLT(L)
    J = JLLT(L)
    ! *** N & S
    IF( IJCTLT(I,J+1) == 9 )THEN
      LNCLT(L) = LCLT
    ELSE
      LNCLT(L) = LIJLT(I,J+1)
    ENDIF
    IF( IJCTLT(I,J-1) == 9 )THEN
      LSCLT(L) = LCLT
    ELSE
      LSCLT(L) = LIJLT(I,J-1)
    ENDIF

    ! *** W & E
    IF( IJCTLT(I+1,J) == 9 )THEN
      LECLT(L) = LCLT
    ELSE
      LECLT(L) = LIJLT(I+1,J)
    ENDIF
    IF( IJCTLT(I-1,J) == 9 )THEN
      LWCLT(L) = LCLT
    ELSE
      LWCLT(L) = LIJLT(I-1,J)
    ENDIF

  ENDDO

  ! *** MODIFY WASP/LT NORTH-SOUTH CELL MAPPING FOR PERIOD GRID IN N-S DIRECTION
  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    DO NPN = 1, NPNSBP
      ! *** SET NORTH CELL SOUTH OF SOUTH CELL
      LS = LIJLT(ISPNS(NPN),JSPNS(NPN))
      LSCLT(LS) = LIJLT(INPNS(NPN),JNPNS(NPN))

      ! *** SET SOUTH CELL NORTH OF NORTH CELL
      LN = LIJLT(INPNS(NPN),JNPNS(NPN))
      LNCLT(LN) = LIJLT(ISPNS(NPN),JSPNS(NPN))
    END DO
  END IF

  ! *** MODIFY 6 CELLS FOR W-E CONNECTIONS (IJCTLT)
  IF( ISCONNECT >= 2 )THEN
    DO NPN = 1,NPEWBP
      LE = LIJLT(IEPEW(NPN),JEPEW(NPN))
      LECLT(LE) = LIJLT(IWPEW(NPN),JWPEW(NPN))

      LW = LIJLT(IWPEW(NPN),JWPEW(NPN))
      LWCLT(LW) = LIJLT(IEPEW(NPN),JEPEW(NPN))
    END DO
  ENDIF

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
  call Communicate_Initialize

  ! *** Write out to MPI status file
  if( DEBUG .or. MPI_Write_Flag )then
    call Write_Cell_Indexing(.true.)
  endif

  return

END Subroutine Cellmap
