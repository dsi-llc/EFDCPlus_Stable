! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  MODULE EFDCOUT

  USE IFPORT
  USE GLOBAL
  USE DRIFTER, ONLY:DRIFTER_OUT
  USE HYDSTRUCMOD
  USE SHELLFISHMOD
  USE FIELDS
  Use Variables_WQ
  Use WQ_RPEM_MODULE
  
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out

  IMPLICIT NONE

  ! *** DO NOT CHANGE THESE PRECISIONS AS THEY ARE LINKED TO EFDC_EXPLORER
  REAL(RKD),PRIVATE,SAVE :: EETIME, SUMTIME

  INTEGER(IK4),PRIVATE,SAVE :: NSXD,NBCCELLS,NCELLLIST,INITBCOUT=0
  INTEGER(IK4),PRIVATE,SAVE :: IWQ(40)
  INTEGER(IK4),PRIVATE,SAVE :: NSEDSTEPS, NBEDSTEPS, NWQVAR, CELL3D
  INTEGER(IK4),PRIVATE,SAVE,ALLOCATABLE,DIMENSION(:) :: BCCELLS

  REAL,PRIVATE,SAVE,ALLOCATABLE,DIMENSION(:,:) :: BCQSUM  ! @todo define!!!
  REAL,PRIVATE,SAVE,ALLOCATABLE,DIMENSION(:,:) :: BCUHDY2 ! @todo define!!!
  REAL,PRIVATE,SAVE,ALLOCATABLE,DIMENSION(:,:) :: BCVHDX2 ! @todo define!!!
  REAL,PRIVATE,SAVE,ALLOCATABLE,DIMENSION(:)   :: BCQSUME ! @todo define!!!

  CHARACTER(30) :: FILENAME

  Character(1) :: EE_PARALLEL_ID !< This is used to write out an integer based on the processor that is writing out

  INTEGER :: EE_UNIT = 95
  INTEGER :: IERRio, NERR, IORIGIN
  CONTAINS

  SUBROUTINE EE_LINKAGE(JSEXPLORER)
  !------------------------------------------------------------------------
  ! **  SUBROUTINE EE_LINKAGE (OLD EEXPOUT.FOR) WRITES BINARY OUTPUT FILES:
  ! **    EE_HYD    - WATER DEPTH AND VELOCITY
  ! **    EE_WC     - WATER COLUMN AND TOP LAYER OF SEDIMENTS
  ! **    EE_BC     - EFDC COMPUTED BOUNDARY FLOWS
  ! **    EE_BED    - SEDIMENT BED LAYER INFORMATION
  ! **    EE_WQ     - WATER QUALITY INFORMATION FOR THE WATER COLUMN
  ! **    EE_SD     - SEDIMENT DIAGENSIS INFORMATION
  ! **    EE_RPEM   - ROOTED PLANT AND EPIPHYTE MODEL
  ! **    EE_ARRAYS - GENERAL/USER DEFINED ARRAY DUMP. LINKED TO
  ! **                EFDC_EXPLORER FOR DISPLAY
  ! **    EE_SEDZLJ - SEDIMENT BED DATA FOR SEDZLJ SUB-MODEL
  !------------------------------------------------------------------------
  INTEGER(IK4),INTENT(IN) :: JSEXPLORER
  INTEGER :: I, IW, J, NS
  
  IF( ISDYNSTP == 0 )THEN
    DELT=DT
  ELSE
    DELT=DTDYN
  ENDIF

  IF( JSEXPLORER == 1 )THEN
    CALL HEADEROUT
  ELSEIF( JSEXPLORER == -1 )THEN
    ! *** FORCE ALL OUTPUT
    NSEDSTEPS=32000
    NBEDSTEPS=32000
  ENDIF

  ! *** SET TIME EE LINKAGE FILES
  EETIME = TIMESEC
  EETIME = EETIME/86400._8

  CALL WSOUT
  CALL VELOUT
  
  ! *** disable writing ou to EE_BC.OUT with MPI
  if(num_Processors > 1 )then
      
  else
    CALL BCOUT
  end if

  IF( ISSPH(8) >= 1 ) CALL WCOUT

  IF( LSEDZLJ )THEN
    CALL SEDZLJOUT
  ELSE
    IF( ISBEXP >= 1 .AND. KB > 1 )THEN
      IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) CALL BEDOUT
    ENDIF
  ENDIF

  IF( ISTRAN(8) > 0 )THEN
    CALL WQOUT
    IF( IWQBEN > 0 .AND. ISSDBIN /= 0 ) CALL SDOUT(JSEXPLORER)
    IF( ISRPEM > 0) CALL RPEMOUT(JSEXPLORER)
    IF( ISFFARM > 0) CALL SHELLFISHOUT()
  ENDIF

  IF( ISINWV == 2 .AND. JSEXPLORER /= 1 ) CALL ARRAYSOUT
  IF( ISPD > 0 .AND. JSEXPLORER == -1 ) CALL DRIFTER_OUT(.TRUE.)

  
  END SUBROUTINE

  SUBROUTINE HEADEROUT
  INTEGER(IK4) :: NACTIVE,VER,HSIZE,BSIZE,I,NS,MW,NW,L,K,ITYPE,ITIMEVAR,LL,N2D,N3D
  CHARACTER(8) :: ARRAYNAME

  if( iswqlvl == 0 .and. NFIXED > 0 )THEN
    NWQVAR = NWQVAR + 1
  ELSE
    NWQVAR = NWQV              ! *** Number of water quality variables no longer changes 
  ENDIF
  NACTIVE = LA_Global - 1

  ! *** NUMBER OF 3D CELLS
  CELL3D = 0
  DO L=2, LA_Global !***Modified for Global
    CELL3D = CELL3D + (KC - KSZ_Global(L) + 1)
  ENDDO

  ! *** WATER DEPTHS
  FILENAME=OUTDIR//'EE_WS.OUT'
  VER   = 8400
  HSIZE = 6*4
  BSIZE = (2 + 3 + NACTIVE)*4
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
  CLOSE(EE_UNIT,STATUS='DELETE')
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
  WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
  WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(NACTIVE,4)
  CLOSE(EE_UNIT,STATUS='KEEP')

  ! *** VELOCITY
  FILENAME=OUTDIR//'EE_VEL.OUT'
  VER   = 8400
  HSIZE = (9 + 4*LCM_Global)*4 !***Modified for MPI
  BSIZE = (2 + 2 + 3*CELL3D)*4
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
  CLOSE(EE_UNIT,STATUS='DELETE')
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
  WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
  WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(NACTIVE,4)
  WRITE(EE_UNIT) REAL(RSSBCE_Global,4)
  WRITE(EE_UNIT) REAL(RSSBCW_Global,4)
  WRITE(EE_UNIT) REAL(RSSBCS_Global,4)
  WRITE(EE_UNIT) REAL(RSSBCN_Global,4)
  CLOSE(EE_UNIT,STATUS='KEEP')

  IF( ISSPH(8) >= 1 )THEN
    ! *** WATER COLUMN AND TOP LAYER OF SEDIMENT
    FILENAME=OUTDIR//'EE_WC.OUT'
    NSXD = NSED+NSND
    VER = 11300
    HSIZE = (29 + NSXD + 1)*4      ! *** +1 is for NSED2
    BSIZE = BLOCKWC(CELL3D)

    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
    CLOSE(EE_UNIT,STATUS='DELETE')
    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
    WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
    WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
    WRITE(EE_UNIT) (INT(ISTRAN(I),4),I=1,7)
    WRITE(EE_UNIT) INT(NDYE,4),INT(NSED,4),INT(NSND,4),INT(NTOX,4),INT(NSED2,4)
    WRITE(EE_UNIT) INT(ISWAVE,4),INT(ISBEDSTR,4),LOGICAL(LSEDZLJ,4),INT(ICALC_BL,4),REAL(TEMBO,4)
    WRITE(EE_UNIT) INT(IEVAP,4),INT(ISGWIE,4),INT(ISICE,4)
    
    DO NS=1,NSXD
      WRITE(EE_UNIT) REAL(SEDDIA(NS),4)
    ENDDO
    CLOSE(EE_UNIT,STATUS='KEEP')
  ENDIF

  ! *** SEDFLUME MODEL RESULTS
  IF( LSEDZLJ )THEN
    ! *** SEDIMENT BED
    FILENAME=OUTDIR//'EE_SEDZLJ.OUT'
    VER=8401
    HSIZE = 22*4
    BSIZE = BLOCKSEDZLJ(CELL3D)

    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
    CLOSE(EE_UNIT,STATUS='DELETE')
    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
    WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
    WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
    WRITE(EE_UNIT) (INT(ISTRAN(I),4),I=1,7)
    WRITE(EE_UNIT) INT(NSCM,4),INT(ITBM,4),INT(NSICM,4),INT(NTOX,4),INT(ICALC_BL,4)
    CLOSE(EE_UNIT,STATUS='KEEP')
    NBEDSTEPS = ISBEXP - 1            ! *** The very first bed snapshot will be saved
  ELSE
    ! *** SEDIMENT BED LAYERS FOR ORIGINAL SEDIMENT TRANSPORT APPROACH
    IF( ISBEXP >= 1 .AND. KB > 1 )THEN
      IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
        FILENAME=OUTDIR//'EE_BED.OUT'
        VER=8400
        HSIZE = (18 + NSXD)*4
        BSIZE = BLOCKBED(CELL3D)

        OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
        CLOSE(EE_UNIT,STATUS='DELETE')
        OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
        WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
        WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
        WRITE(EE_UNIT) (INT(ISTRAN(I),4),I=1,7)
        WRITE(EE_UNIT) INT(NSED,4),INT(NSND,4),INT(NTOX,4)

        DO NS=1,NSXD
          WRITE(EE_UNIT)REAL(SEDDIA(NS),4)
        ENDDO
        CLOSE(EE_UNIT,STATUS='KEEP')
        NBEDSTEPS = ISBEXP - 1            ! *** The very bed first snapshot will be saved, regardless of the skip count
      ENDIF
    ENDIF
  ENDIF

  ! *** WATER QUALITY MODEL (HEM3D) RESULTS
  IF( ISTRAN(8) > 0 )THEN
    FILENAME=OUTDIR//'EE_WQ.OUT'

    IWQ = 0
    BSIZE = 0

    DO NW = 1,NWQVAR
      IWQ(NW) = ISKINETICS(NW)
      IF( IWQ(NW) > 0) BSIZE = BSIZE + 1
    ENDDO
    
    VER = 10300                 !< version 10.3
    HSIZE = (13 + NWQVAR)*4
    BSIZE = (2 + BSIZE*CELL3D)*4
    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
    CLOSE(EE_UNIT,STATUS='DELETE')
    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
    WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4),INT(IGRIDV,4),INT(CELL3D,4)
    WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
    WRITE(EE_UNIT) INT(NALGAE,4),INT(NZOOPL,4),INT(NWQVAR,4)
    WRITE(EE_UNIT) (INT(IWQ(NW),4),NW=1,NWQVAR)
    CLOSE(EE_UNIT,STATUS='KEEP')

    ! *** SAVE SEDIMENT DIAGENESIS RESULTS
    IF( ISSDBIN /= 0 )THEN
      FILENAME=OUTDIR//'EE_SD.OUT'
      VER=8400
      HSIZE = 8*4
      BSIZE = (2 + (6*3 + 21)*NACTIVE)*4
      OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
      CLOSE(EE_UNIT,STATUS='DELETE')
      OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
      WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
      CLOSE(EE_UNIT,STATUS='KEEP')
      NSEDSTEPS=-1
    ENDIF

    IF( ISRPEM > 0 )THEN
      FILENAME=OUTDIR//'EE_RPEM.OUT'
      VER = 8400
      HSIZE = (10 + NACTIVE)*4
      BSIZE = 0
      DO L=2,LA_Global
        IF( LMASKRPEM_Global(L) .OR. INITRPEM == 0 )THEN
          BSIZE = BSIZE + 1
        ENDIF
      ENDDO
      BSIZE = (2 + 4*BSIZE)*4

      OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
      CLOSE(EE_UNIT,STATUS='DELETE')
      OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
      WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(KB,4),INT(NACTIVE,4)
      WRITE(EE_UNIT) INT(NRPEM,4),INT(INITRPEM,4)
      WRITE(EE_UNIT) (LOGICAL(LMASKRPEM_Global(L),4),L=2,LA_Global)
      CLOSE(EE_UNIT,STATUS='KEEP')
    ENDIF

    IF( ISFFARM > 0 )THEN
      FILENAME=OUTDIR//'EE_SHF.OUT'
      VER = 10000
      HSIZE = (10 + NACTIVE)*4
      N2D = 0
      N3D = 0
      DO L=2,LA_Global
        IF( FARMCELL(L) > 0 )THEN
          N2D = N2D + 1
          N3D = N3D + (KC - KSZ_Global(L) + 1)
        ENDIF
      ENDDO
      BSIZE = (2 + NSF*(12*N2D + 11*N3D))*4

      OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
      CLOSE(EE_UNIT,STATUS='DELETE')
      OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      WRITE(EE_UNIT) INT(VER,4),INT(HSIZE,4),INT(BSIZE,4)
      WRITE(EE_UNIT) INT(IC_Global,4),INT(JC_Global,4),INT(KC,4),INT(NACTIVE,4)
      WRITE(EE_UNIT) INT(NSF,4),INT(NSFCELLS,4),INT(NSFHARV,4)
      WRITE(EE_UNIT) (INT(FARMCELL(L),4),L=2,LA_Global)
      CLOSE(EE_UNIT,STATUS='KEEP')
    ENDIF

  ENDIF

  ! *** BC OUTPUT: QSUM
  !CALL BCOUT_INITIALIZE !***Initializes some arrays that hold some info later on.. not clearly explained

  ! *** Set file name and write out to #output directory
  FILENAME=OUTDIR//'EE_BC.OUT'
  ! *** Setup binary file to write to
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
  CLOSE(EE_UNIT,STATUS='DELETE')
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')

  ! *** ONLY OUTPUT CELLS THAT HAVE BC'S
  IF( KC > 1 )THEN  ! *** If more than one layer
    NBCCELLS = NBCS + NBCSOP  ! Number of boundary condition cells?
    NBCCELLS = NBCCELLS + LA_Global                                      ! *** EVAP/RAINFALL/ICE
    IF( NGWSER > 0 .OR. ISGWIT /= 0 ) NBCCELLS = NBCCELLS + LA_Global    ! *** GROUNDWATER

    ALLOCATE(BCCELLS(NBCCELLS)) !*** @todo what is this

    ! *** BUILD CELL LIST FOR SELECTIVE QSUM
    NCELLLIST = 0
    ! *** BC CELLS
    DO LL=1, NBCS
      NCELLLIST = NCELLLIST+1
      BCCELLS(NCELLLIST) = LBCS(LL)
    ENDDO

    ! *** OPEN BC CELLS
    DO LL=1, NBCSOP
      NCELLLIST = NCELLLIST+1
      BCCELLS(NCELLLIST) = LOBCS(LL)
    ENDDO
  ELSE
    ! *** SINGLE LAYER
    NBCCELLS = LA_Global
    NCELLLIST = 0
    ALLOCATE(BCCELLS(1))
  ENDIF
  ! *** Setup values for writing to header

  VER   = 8400 !*** version number, @todo should this be updated??
  HSIZE = (20 + NCELLLIST + NBCS + NPBS + NPBW + NPBE + NPBN)*4 !***Header size, why is there a +20?
  BSIZE = BLOCKBC(CELL3D) !***Block size? what is cell3d supposed to do?

  WRITE(EE_UNIT) INT(VER,4), INT(HSIZE,4), INT(BSIZE,4)
  WRITE(EE_UNIT) INT(IC_Global,4), INT(JC_Global,4), INT(KC,4), INT(NACTIVE,4) !***Modified for global
  WRITE(EE_UNIT) INT(NBCS,4), INT(NBCCELLS,4), INT(NCELLLIST,4)
  WRITE(EE_UNIT) (INT(BCCELLS(L),4),L=1,NCELLLIST)
  WRITE(EE_UNIT) INT(NPBS,4),  INT(NPBW,4), INT(NPBE,4),    INT(NPBN,4)
  WRITE(EE_UNIT) INT(NQCTL,4), INT(NQWR,4), INT(NQCTLSER,4),INT(NQCTRULES,4)
  WRITE(EE_UNIT) INT(NGWSER,4),INT(ISGWIT,4)
  WRITE(EE_UNIT) (INT(LBCS(L),4),L=1,NBCS) 
  WRITE(EE_UNIT) (INT(LPBS(L),4),L=1,NPBS)
  WRITE(EE_UNIT) (INT(LPBW(L),4),L=1,NPBW)
  WRITE(EE_UNIT) (INT(LPBE(L),4),L=1,NPBE)
  WRITE(EE_UNIT) (INT(LPBN(L),4),L=1,NPBN)
  CLOSE(EE_UNIT,STATUS='KEEP')

  ! *** End BC OUT

  IF( ISINWV == 2 )THEN
    FILENAME=OUTDIR//'EE_ARRAYS.OUT'
    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
    CLOSE(EE_UNIT,STATUS='DELETE')
    OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')

    ! *** EE LINKAGE VERSION
    VER=7300
    WRITE(EE_UNIT) INT(VER,4)

    ! *** NUMBER OF TIME VARYING ARRAYS
    NS = 4 + 4
    IF( ISHDMF >= 1 ) NS = NS + 4  ! NOTE: THIS IS ONLY NEEDED AS LONG AS FMDUX AND FMDUY ARE BEING USED
    WRITE(EE_UNIT) INT(NS,4)

    ! FLAGS: ARRAY TYPE, TIME VARIABLE
    ! ARRAY TYPE:    0 = L            DIM'D
    !                1 = L,KC         DIM'D
    !                2 = L,0:KC       DIM'D
    !                3 = L,KB         DIM'D
    !                4 = L,KC,NCLASS  DIM'D
    ! TIME VARIABLE: 0 = NOT CHANGING
    !                1 = TIME VARYING

    ! *****************************************************************
    ! *** THE FOLLOWING ARRAYS ARE TIME INVARIANT SO ONLY WRITTEN ONCE
    ! *****************************************************************
    ITIMEVAR = 0
    ITYPE = 0

    WRITE(EE_UNIT) INT(ITYPE,4),INT(ITIMEVAR,4)
    ARRAYNAME='SBX'
    WRITE(EE_UNIT) ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT) REAL(SBXO(L),4)
    ENDDO

    WRITE(EE_UNIT) INT(ITYPE,4),INT(ITIMEVAR,4)
    ARRAYNAME='SBY'
    WRITE(EE_UNIT) ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT) REAL(SBYO(L),4)
    ENDDO

    ITYPE = 1
    WRITE(EE_UNIT) INT(ITYPE,4),INT(ITIMEVAR,4)
    ARRAYNAME='DZC'
    WRITE(EE_UNIT) ARRAYNAME
    DO K=1,KC
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(DZC(L,K),4)
      ENDDO
    ENDDO

    CLOSE(EE_UNIT,STATUS='KEEP')

  ENDIF
  END SUBROUTINE

  SUBROUTINE WSOUT

  IMPLICIT NONE

  ! ** OUTPUT WATER DEPTH
  INTEGER(IK4) :: VER, HSIZE, BSIZE, I, J, LG, LL, L2
  INTEGER(IK4) :: L, ITMP, LTMP, NN, NS, ISTAT
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME

  FILENAME=OUTDIR//'EE_WS.OUT'
  EE_UNIT = 96                        ! *** EE_WS.OUT

  ! *** If we are restarting and continuing a run
  IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_WS == 0 )THEN
    RSTFIRST_WS = 1
    WRITE(*,'(A)')'READING TO STARTING TIME FOR WS'
    FSIZE = FILESIZE(FILENAME)

    OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
    READ(EE_UNIT) VER,HSIZE,BSIZE
    IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
    OFFSET = HSIZE + 4

    NS = 0
    DO WHILE(OFFSET < FSIZE)
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      IF( EOF(EE_UNIT) )THEN
        WRITE(EE_UNIT) INT(N+NRESTART,4), EETIME,REAL(DELT,4), INT(LMINSTEP,4)
        EXIT
      ENDIF
      IF( ISTAT /= 0 ) EXIT

      READ(EE_UNIT) PTIME

      IF( DEBUG )THEN
        NS=NS+1
        IF( NS == 8 )THEN
          WRITE(*,'(F10.3)')PTIME
          NS=0
        ELSE
          WRITE(*,'(F10.3,\)')PTIME
        ENDIF
      ENDIF
      IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

      OFFSET = OFFSET + BSIZE
    ENDDO

    IF( DEBUG ) WRITE(*,'(" ")')
    WRITE(*,'(A)')'FINSIHED READING WS'

    ISTAT = FSEEK(EE_UNIT,8,1)

  ELSE !***Normal write out
    ! *** Writing out some header information for the current time step
    NERR = 0
    IORIGIN = 1
    100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
    WRITE(EE_UNIT) INT(N+NRESTART,4), EETIME, REAL(DELT,4),INT(LMINSTEP,4)
  ENDIF

  IF( ISRESTI == 0 .OR. ICONTINUE == 0) NRESTART=0

  ! *** Check for anti-reflection type BC's and replace the interior cell depths with the boundary values
  DO LL=1,NPBW_GL
    IF( ISPBW_GL(LL) == 4 .OR. ISPBW_GL(LL) == 5 )THEN
      I = IG2IL(IPBW_GL(LL))
      J = JG2JL(JPBW_GL(LL))
      IF(  I > 0 .AND. I <= IC_GLOBAL )THEN
        IF(  J > 0 .AND. J <= JC_GLOBAL )THEN
          LG = LIJ_GLOBAL(I,J) 
          L2 = LEC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        END IF
      END IF
    END IF
  END DO
  DO LL=1,NPBE_GL
    IF( ISPBE_GL(LL) == 4 .OR. ISPBE_GL(LL) == 5 )THEN
      I = IG2IL(IPBE_GL(LL))
      J = JG2JL(JPBE_GL(LL))
      IF(  I > 0 .AND. I <= IC_GLOBAL )THEN
        IF(  J > 0 .AND. J <= JC_GLOBAL )THEN
          LG = LIJ_GLOBAL(I,J) 
          L2 = LWC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        END IF
      END IF
    END IF
  END DO
  DO LL=1,NPBS_GL
    IF( ISPBS_GL(LL) == 4 .OR. ISPBS_GL(LL) == 5 )THEN
      I = IG2IL(IPBS_GL(LL))
      J = JG2JL(JPBS_GL(LL))
      IF(  I > 0 .AND. I <= IC_GLOBAL )THEN
        IF(  J > 0 .AND. J <= JC_GLOBAL )THEN
          LG = LIJ_GLOBAL(I,J) 
          L2 = LNC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        END IF
      END IF
    END IF
  END DO
  DO LL=1,NPBN_GL
    IF( ISPBN_GL(LL) == 4 .OR. ISPBN_GL(LL) == 5 )THEN
      I = IG2IL(IPBN_GL(LL))
      J = JG2JL(JPBN_GL(LL))
      IF(  I > 0 .AND. I <= IC_GLOBAL )THEN
        IF(  J > 0 .AND. J <= JC_GLOBAL )THEN
          LG = LIJ_GLOBAL(I,J) 
          L2 = LSC_GLOBAL(LG)
          HP_GLOBAL(L2) = (HP_GLOBAL(LG) + BELV_GLOBAL(LG)) - BELV_GLOBAL(L2)
        END IF
      END IF
    END IF
  END DO
 
  WRITE(EE_UNIT) (REAL(HP_Global(L),4), L=2,LA_Global)

  FLUSH(EE_UNIT)
  CLOSE(EE_UNIT,STATUS='KEEP')
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_WS.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF

END SUBROUTINE

SUBROUTINE VELOUT

  ! ** WATER DEPTH AND VELOCITY OUTPUT
  INTEGER(IK4) :: VER,HSIZE,BSIZE
  INTEGER(IK4) :: L,K,ITMP,NN,NS,ISTAT
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME

  FILENAME=OUTDIR//'EE_VEL.OUT'
  EE_UNIT = 97                        ! *** EE_VEL.OUT
  IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_VEL == 0 )THEN
    RSTFIRST_VEL=1
    WRITE(*,'(A)')'READING TO STARTING TIME FOR VEL'
    FSIZE = FILESIZE(FILENAME)
    OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
    READ(EE_UNIT) VER,HSIZE,BSIZE
    IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
    OFFSET = HSIZE + 4

    NS = 0
    DO WHILE(OFFSET < FSIZE)
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      IF( EOF(EE_UNIT) )THEN
        WRITE(EE_UNIT) INT(N+NRESTART,4),EETIME,REAL(DELT,4)
        EXIT
      ENDIF
      IF( ISTAT /= 0 ) EXIT

      READ(EE_UNIT) PTIME

      IF( DEBUG )THEN
        NS=NS+1
        IF( NS == 8 )THEN
          WRITE(*,'(F10.3)')PTIME
          NS=0
        ELSE
          WRITE(*,'(F10.3,\)')PTIME
        ENDIF
      ENDIF
      IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

      OFFSET = OFFSET + BSIZE
    ENDDO
    IF( DEBUG ) WRITE(*,'(" ")')
    WRITE(*,'(A)')'FINISHED READING VEL'

    ISTAT = FSEEK(EE_UNIT,4,1)

  ELSE
    NERR = 0
    IORIGIN = 2
    100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
    WRITE(EE_UNIT) INT(N+NRESTART,4), EETIME, REAL(DELT,4)
  ENDIF

  IF( ISRESTI == 0 .OR. ICONTINUE == 0) NRESTART=0
  
  WRITE(EE_UNIT) ((REAL(U_Global(L,K),4), K=KSZ_Global(L), KC), L=2, LA_Global)
  WRITE(EE_UNIT) ((REAL(V_Global(L,K),4), K=KSZ_Global(L), KC), L=2, LA_Global)
  WRITE(EE_UNIT) ((REAL(W_Global(L,K),4), K=KSZ_Global(L), KC), L=2, LA_Global)

  FLUSH(EE_UNIT)
  CLOSE(EE_UNIT,STATUS='KEEP')
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_VEL.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF
  
END SUBROUTINE

SUBROUTINE WCOUT

  ! ** WATER COLUMN OF CONSTITUENTS OUTPUT
  INTEGER(IK4) :: VER, I, ITMP, NSED4, NSND4, NTOX4
  INTEGER(IK4) :: NS, HSIZE, BSIZE, ISTAT
  INTEGER(IK4) :: L, K, NT, NX, MD
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP, SHEAR
  REAL(RKD)    :: PTIME
  LOGICAL      :: LTMP

  FILENAME=OUTDIR//'EE_WC.OUT'
  EE_UNIT = 98                        ! *** EE_WC.OUT
  IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_WC == 0 )THEN
    RSTFIRST_WC=1
    WRITE(*,'(A)')'READING TO STARTING TIME FOR WC'
    FSIZE = FILESIZE(FILENAME)
    OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
    READ(EE_UNIT) VER,HSIZE,BSIZE
    IF( VER /= 8400 .AND. VER /= 8500 .AND. VER /=11300) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
    OFFSET = HSIZE

    NS = 0
    DO WHILE(OFFSET < FSIZE)
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      IF( EOF(EE_UNIT) )THEN
        WRITE(EE_UNIT) EETIME
        EXIT
      ENDIF
      IF( ISTAT /= 0 ) EXIT

      READ(EE_UNIT) PTIME

      IF( DEBUG )THEN
        NS=NS+1
        IF( NS == 8 )THEN
          WRITE(*,'(F10.3)')PTIME
          NS=0
        ELSE
          WRITE(*,'(F10.3,\)')PTIME
        ENDIF
      ENDIF
      IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

      OFFSET = OFFSET + BSIZE
    ENDDO
    IF( DEBUG ) WRITE(*,'(" ")')
    WRITE(*,'(A)')'FINISHED READING WC'

  ELSE
    NERR = 0
    IORIGIN = 3
    !PRINT *, 'OPEN: NITER = ', NITER
    100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
    WRITE(EE_UNIT) EETIME
  ENDIF
  
  ! *** WRITE THE TOP LAYER INDEX
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    WRITE(EE_UNIT) (INT(KBT_Global(L),4), L=2,LA_Global)
  ENDIF

  ! *** TOTAL BED SHEAR STRESS
  DO L=2,LA_Global
    WRITE(EE_UNIT) REAL(SHEAR_Global(L),4)
  ENDDO

  IF( ISBEDSTR == 1 .AND. .NOT. LSEDZLJ )THEN
    ! *** Stress from non-cohesive components
    DO L=2,LA_Global
      WRITE(EE_UNIT) REAL(SHEAR_Global2(L),4)
    ENDDO
  ENDIF

  IF( ISWAVE >= 1 )THEN
    IF( LSEDZLJ )THEN
      QQWV3_Global = QQWV3_Global* 0.1 / (WATERDENS*1000.)   ! *** Convert from dyne/cm2 to Pa and density normalize
    ENDIF
    WRITE(EE_UNIT) (REAL(QQWV3_Global(L),4), L=2,LA_Global)  ! *** Bed Shear due to Waves Only
    ! *** Shear due to Current Only
    IF( LSEDZLJ )THEN
      DO L=2,LA_Global
        SHEAR = SHEAR_Global(L) - QQWV3_Global(L)
        !SHEAR = MAX(SHEAR,0.0)
        WRITE(EE_UNIT) SHEAR             ! *** Bed Shear due to Current Only
      ENDDO
    ELSE
      DO L=2,LA_Global
        SHEAR = ( RSSBCE_Global(L)*TBX_Global(LEC_Global(L)) + RSSBCW_Global(L)*TBX_Global(L) )**2  + &
                ( RSSBCN_Global(L)*TBY_Global(LNC_Global(L)) + RSSBCS_Global(L)*TBY_Global(L) )**2
        SHEAR = 0.5*SQRT(SHEAR)
        WRITE(EE_UNIT) SHEAR             ! *** Bed Shear due to Current Only
      ENDDO
    ENDIF
    IF( ISWAVE >= 3 )THEN
      WRITE(EE_UNIT) (REAL(WV_HEIGHT_Global(L),4), L=2,LA_Global)
      WRITE(EE_UNIT) (REAL(WV_FREQ_Global(L),4), L=2,LA_Global)
      WRITE(EE_UNIT) (REAL(WV_DIR_Global(L),4), L=2,LA_Global)
      IF( ISWAVE == 4 )THEN
        WRITE(EE_UNIT) (REAL(WV_DISSIPA_Global(L),4), L=2,LA_Global)      ! *** DISSIPATION
        WRITE(EE_UNIT) (REAL(WVHUU_Global(L,KC),4), L=2,LA_Global)        ! *** SXX (M3/S2)
        WRITE(EE_UNIT) (REAL(WVHVV_Global(L,KC),4), L=2,LA_Global)        ! *** SYY (M3/S2)
        WRITE(EE_UNIT) (REAL(WVHUV_Global(L,KC),4), L=2,LA_Global)        ! *** SXY (M3/S2)
      ENDIF
    ENDIF
  ENDIF

  IF( ISTRAN(1) >= 1 ) WRITE(EE_UNIT) ((REAL(SAL_Global(L,K),4), K=KSZ_Global(L),KC), L=2,LA_Global)

  IF( ISTRAN(2) >= 1 )THEN
    WRITE(EE_UNIT) ((REAL(TEM_Global(L,K),4), K=KSZ_Global(L),KC), L=2,LA_Global)
    IF( TEMBO > 0.) WRITE(EE_UNIT) (REAL(TEMB_Global(L),4), L=2,LA_Global)
    IF( IEVAP > 1 )THEN
      WRITE(EE_UNIT) (REAL(EVAPT_Global(L),4), L=2,LA_Global)
      WRITE(EE_UNIT) (REAL(RAINT_Global(L),4), L=2,LA_Global)
    ENDIF
  ENDIF

  IF( ISTRAN(3) >= 1 ) WRITE(EE_UNIT) (((REAL(DYE_Global(L,K,MD),4), K=KSZ_Global(L),KC), L=2,LA_Global), MD=1,NDYE)

  IF( ISTRAN(4) >= 1 ) WRITE(EE_UNIT) ((REAL(SFL_Global(L,K),4), K=KSZ_Global(L),KC), L=2,LA_Global)

  IF( ISTRAN(5) >= 1 )THEN
    WRITE(EE_UNIT) ((REAL(TOXB_Global(L,KBT_Global(L),NT),4), L=2,LA_Global), NT=1,NTOX)
    
    WRITE(EE_UNIT) (((REAL(TOX_Global(L,K,NT),4), K=KSZ_Global(L),KC), L=2,LA_Global), NT=1,NTOX)
  ENDIF

  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    WRITE(EE_UNIT) (REAL(BELV_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(HBED_Global(L,KBT_Global(L)),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(BDENBED_Global(L,KBT_Global(L)),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(PORBED_Global(L,KBT_Global(L)),4), L=2,LA_Global)
    IF( ISTRAN(6) >= 1 ) WRITE(EE_UNIT) ((REAL(SEDB_Global(L,KBT_Global(L),NS),4), L=2,LA_Global), NS=1,NSED)
    IF( ISTRAN(7) >= 1 ) WRITE(EE_UNIT) ((REAL(SNDB_Global(L,KBT_Global(L),NX),4), L=2,LA_Global), NX=1,NSND)
    WRITE(EE_UNIT) ((REAL(0.0,4), L=2,LA_Global), NS=1,NSED+NSND)
    IF( ISTRAN(6) >= 1 ) WRITE(EE_UNIT) (((REAL(SED_Global(L,K,NS),4), K=KSZ_Global(L),KC), L=2,LA_Global), NS=1,NSED2)
    IF( ISTRAN(7) >= 1 )THEN
      WRITE(EE_UNIT) (((REAL(SND_GLobal(L,K,NX),4), K=KSZ_Global(L),KC), L=2,LA_Global), NX=1,NSND)
      IF( ICALC_BL > 0 .AND. NSND > 0 )THEN
        ! *** CALCULATE EQUIVALENT CONCENTRATIONS
        WRITE(EE_UNIT) ((REAL(QSBDLDX_Global(L,NX),4), L=2,LA_Global), NX=1,NSND)
        WRITE(EE_UNIT) ((REAL(QSBDLDY_Global(L,NX),4), L=2,LA_Global), NX=1,NSND)
      ENDIF
    ENDIF
  ELSEIF (BATHY.IFLAG > 0 )THEN
      WRITE(EE_UNIT) (REAL(BELV_Global(L),4), L=2,LA_Global)
  ENDIF
  IF( ISGWIE > 0 )THEN
    WRITE(EE_UNIT) (REAL(EVAPSW_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(EVAPGW_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(QGW_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(AGWELV_Global(L),4), L=2,LA_Global)
  ENDIF

  IF( ISTRAN(2) > 0 .AND. ISICE  >= 3 )THEN
    WRITE(EE_UNIT) (REAL(ICETHICK_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(ICETEMP_Global(L),4), L=2,LA_Global)
  ENDIF

  FLUSH(EE_UNIT)
  CLOSE(EE_UNIT,STATUS='KEEP')
  !PRINT *, 'CLOSE: NITER = ', NITER
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_WC.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF

END SUBROUTINE

INTEGER FUNCTION BLOCKWC(CELL3D)

  INTEGER:: DSIZE, CELL2D, CELL3D

  CELL2D = LA_Global - 1

  DSIZE = 2 !**EETIME
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    DSIZE = DSIZE + CELL2D                      ! TOP LAYER
  ENDIF

  ! *** READ THE WATER COLUMN AND TOP LAYER OF SEDIMENT DATA, IF NEEDED
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    IF( LSEDZLJ )THEN
      DSIZE = DSIZE + CELL2D                  ! SHEAR
    ELSEIF( ISBEDSTR >= 1 )THEN
      DSIZE = DSIZE + CELL2D                  ! TAUBSED
      IF( ISBEDSTR == 1 )THEN
        DSIZE = DSIZE + CELL2D                ! TAUBSND
      ENDIF
    ELSE
      ! *** TOTAL BED SHEAR STRESS
      DSIZE = DSIZE + CELL2D                  ! TAUB
    ENDIF
  ELSE
    ! *** TOTAL BED SHEAR STRESS
    DSIZE = DSIZE + CELL2D                    ! SHEAR
  ENDIF
  IF( ISWAVE >= 1 )THEN
    DSIZE = DSIZE + 2*CELL2D                  ! QQWV3, SHEAR
    IF( ISWAVE >= 3 )THEN
      DSIZE = DSIZE + 3*CELL2D                ! WV(L).HEIGHT, WV.FREQ, WACCWE
      IF( ISWAVE == 4 )THEN
        DSIZE = DSIZE + 4*CELL2D              ! WV.DISSIPA, WVHUU, WVHVV, WVHUV
      ENDIF
    ENDIF
  ENDIF
  IF( ISTRAN(1) >= 1 ) DSIZE = DSIZE + CELL3D   ! SAL_Global
  IF( ISTRAN(2) >= 1 )THEN
    DSIZE = DSIZE + CELL3D                      ! TEM_Global
    IF( TEMBO > 0.) DSIZE = DSIZE + CELL2D      ! TEMB
    IF( IEVAP > 1 )THEN                         
      DSIZE = DSIZE + 2 * CELL2D                ! EVAPT, RAINT
    ENDIF
  ENDIF
  IF( ISTRAN(3) >= 1 ) DSIZE = DSIZE + NDYE*CELL3D   ! DYE
  IF( ISTRAN(4) >= 1 ) DSIZE = DSIZE + CELL3D        ! SFL
  IF( ISTRAN(5) >= 1 )THEN
    DSIZE = DSIZE + NTOX*CELL2D               ! TOXB
    DSIZE = DSIZE + NTOX*CELL3D               ! TOX
  ENDIF
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    DSIZE = DSIZE + 4*CELL2D                  ! BELV, HBED, BDENBED, PORBED
    IF( ISTRAN(6) == 1 )THEN
      DSIZE = DSIZE + 2*NSED*CELL2D           ! SEDB, VFRBED
      DSIZE = DSIZE + NSED2*CELL3D            ! SED
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DSIZE = DSIZE + 2*NSND*CELL2D           ! SNDB, VFRBED
      DSIZE = DSIZE + NSND*CELL3D             ! SND
      IF( ICALC_BL > 0 .AND. NSND > 0 )THEN
        DSIZE = DSIZE + 2*NSND*CELL2D         ! QSBDLDX, QSBDLDY
      ENDIF
    ENDIF
  ELSEIF (BATHY.IFLAG > 0 )THEN
      DSIZE = DSIZE + CELL2D                  ! BELV
  ENDIF
  ! ** NEW BLOCK
  IF( ISGWIE > 0 )THEN
    DSIZE = DSIZE + 4*CELL2D                  ! EVAPSW, EVAPGW, QGW, AGWELV
  ENDIF
  IF(ISTRAN(2) > 0 .AND. ISICE  >= 3)THEN
    DSIZE = DSIZE + 2*CELL2D                  ! ICETHICK, ICETEMP
  ENDIF
  BLOCKWC = 4*DSIZE                             ! SINGLE PRECISION
  END FUNCTION

  SUBROUTINE SEDZLJOUT

  INTEGER(IK4) :: VER, HSIZE, BSIZE, ISTAT
  INTEGER(IK4) :: L, NS, NT, K, ITMP
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP, SURFACE
  REAL(RKD)    :: PTIME

  NBEDSTEPS = NBEDSTEPS + 1
  IF( NBEDSTEPS >= ISBEXP )THEN
    FILENAME=OUTDIR//'EE_SEDZLJ.OUT'
    EE_UNIT = 99                        ! *** EE_SEDZLJ.OUT
    IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_SEDZLJ == 0 )THEN
      RSTFIRST_SEDZLJ=1
      WRITE(*,'(A)')'READING TO STARTING TIME FOR SEDZLJ'
      FSIZE = FILESIZE(FILENAME)
      OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
      READ(EE_UNIT) VER,HSIZE,BSIZE
      IF(VER /= 8401) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
      OFFSET = HSIZE

      NS = 0
      DO WHILE(OFFSET < FSIZE)
        ISTAT = FSEEK(EE_UNIT,OFFSET,0)

        IF( EOF(EE_UNIT) )THEN
          WRITE(EE_UNIT) EETIME
          EXIT
        ENDIF
        IF( ISTAT /= 0 ) EXIT

        READ(EE_UNIT) PTIME

        IF( DEBUG )THEN
          NS=NS+1
          IF( NS == 8 )THEN
            WRITE(*,'(F10.3)')PTIME
            NS=0
          ELSE
            WRITE(*,'(F10.3,\)')PTIME
          ENDIF
        ENDIF
        IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

        OFFSET = OFFSET + BSIZE
      ENDDO
      IF( DEBUG ) WRITE(*,'(" ")')
      WRITE(*,'(A)')'FINISHED READING SEDZLJ'

    ELSE
      NERR = 0
      IORIGIN = 5
      100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
      WRITE(EE_UNIT) EETIME
    ENDIF
    IF( DTSED == 0. ) DTSED = DT

    ! *** This is used to write out the active layers instead of using the LAYERACTIVE array. 
    DO L=2,LA_Global
      DO K=1,KB
        IF( TSED_GLOBAL(K,L) > 1E-8 )THEN
          WRITE(EE_UNIT) INT(1,4)                                                    ! *** LAYERACTIVE(KB,LCM) - This is = 1 when a bed layer (KB index) exists with mass
        ELSE
          WRITE(EE_UNIT) INT(0,4)
        ENDIF
      ENDDO
    ENDDO
  
    ! *** REAL*4  - Global arrays use reversed indicies for MPI mapping.  Writes the same order as OMP original
    WRITE(EE_UNIT) (REAL(TAU_Global(L),4), L=2,LA_Global)                                  ! *** TAU(LCM)      - Shear Stress in dynes/cm^2
    WRITE(EE_UNIT) ((REAL(TSED_Global(K,L),4), K=1,KB), L=2,LA_Global)                     ! *** TSED(KB,LCM)  - This is the mass in g/cm^2 in each layer
    WRITE(EE_UNIT) ((REAL(BULKDENS_Global(K,L),4), K=1,KB), L=2,LA_Global)                 ! *** BULKDENS(KB,LCM) - Dry Bulk density of each layer (g/cm^3)
    WRITE(EE_UNIT) (((REAL(PERSED_Global(NS,K,L),4), K=1,KB), L=2,LA_Global), NS=1,NSCM)   ! *** PERSED(NSCM,KB,LCM) - This is the mass percentage of each size class in a layer
    WRITE(EE_UNIT) (REAL(D50AVG_Global(L),4), L=2,LA_Global)                               ! *** D50AVG(LCM)   - Average particle size of bed surface (microns)
                                                                                      
    WRITE(EE_UNIT) (REAL(ETOTO_Global(L),4), L=2,LA_Global)                                ! *** ETOTO(LCM)    - Total erosion rate in the cell g/cm^2/s
    WRITE(EE_UNIT) (REAL(DEPO_Global(L),4), L=2,LA_Global)                                 ! *** DEPO(LCM)     - Total deposition rate in the cell g/cm^2/s
                                                                                      
    IF( ICALC_BL > 0 )THEN                                                              
      WRITE(EE_UNIT) ((REAL(CBL_Global(L,NS)*10000.,4), L=2,LA_Global), NS=1,NSCM)         ! *** Bedload mass inb g/m^2 for each size class
      WRITE(EE_UNIT) ((REAL(QSBDLDX_Global(L,NS),4), L=2,LA_Global), NS=1,NSCM)            ! *** Bedload flux in X direction (g/s)
      WRITE(EE_UNIT) ((REAL(QSBDLDY_Global(L,NS),4), L=2,LA_Global), NS=1,NSCM)            ! *** Bedload flux in Y direction (g/s)
    ENDIF

    IF( ISTRAN(5) > 0 )THEN
      DO NT=1,NTOX
        DO L=2,LA_Global
          DO K=1,KB
            IF( K < 3 .AND. TSED_Global(K,L) > 0. .AND. HBED_Global(L,KBT_Global(L)) > 0.)THEN  ! *** avoid division by HBED = 0 --- DKT
              TMP = 0.01*TSED_Global(K,L)/BULKDENS_Global(K,L)             ! *** HBED(L,K)
              TMP = TMP/HBED_Global(L,KBT_Global(L))
              TMP = TMP*TOXB_Global(L,KBT_Global(L),NT)
              WRITE(EE_UNIT) TMP
            ELSEIF( TSED_Global(K,L) > 0. )THEN
              WRITE(EE_UNIT) REAL(TOXB_Global(L,K,NT),4)
            ELSE
              WRITE(EE_UNIT) 0.0_4
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      IF( ICALC_BL > 0 )THEN
        WRITE(EE_UNIT) ((REAL(CBLTOX_Global(L,NT),4), L=2,LA_Global), NT=1,NTOX)    ! *** Bedload toxic concentration (mg/m^2)
      ENDIF
    ENDIF

    FLUSH(EE_UNIT)
    CLOSE(EE_UNIT,STATUS='KEEP')
    
    NBEDSTEPS = 0
  ENDIF
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_SEDZLJ.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF
  
END SUBROUTINE

INTEGER FUNCTION BLOCKSEDZLJ(CELL3D)

  INTEGER:: DSIZE, CELL2D, CELL3D

  CELL2D = LA_Global - 1
  DSIZE = 2                             ! *** EETIME

  ! *** INTEGER*4
  DSIZE = DSIZE + KB*CELL2D             ! *** LAYERACTIVE

  ! *** REAL*4
  DSIZE = DSIZE + 4*CELL2D              ! *** TAU,D50AVG,ETOTO,DEP
  DSIZE = DSIZE + 2*KB*CELL2D           ! *** TSED,BULKDENS
  DSIZE = DSIZE + KB*NSCM*CELL2D        ! *** PERSED
  IF( ICALC_BL > 0 )THEN
    DSIZE = DSIZE + 3*NSCM*CELL2D       ! *** CBL,QSBDLDX,QSBDLDY
  ENDIF
  IF( ISTRAN(5) >= 1 )THEN
    DSIZE = DSIZE + NTOX*KB*CELL2D      ! *** TOXB
    IF(  ICALC_BL > 0 )THEN
      DSIZE = DSIZE + NTOX*CELL2D       ! *** CBLTOX
    ENDIF
  ENDIF
  BLOCKSEDZLJ = 4*DSIZE                 ! *** CONVERT TO BYTES
  END FUNCTION

  SUBROUTINE BEDOUT

  INTEGER(IK4) :: VER, HSIZE, BSIZE, ISTAT
  INTEGER(IK4) :: I, NS, L, K, NX, NT, ITMP
  INTEGER(IK8) :: FSIZE,  OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME

  ! *** SEDIMENT BED LAYERS
  NBEDSTEPS = NBEDSTEPS + 1
  IF( NBEDSTEPS >= ISBEXP )THEN
    FILENAME=OUTDIR//'EE_BED.OUT'
    IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_BED == 0 )THEN
      RSTFIRST_BED=1
      WRITE(*,'(A)')'READING TO STARTING TIME FOR BED'
      FSIZE = FILESIZE(FILENAME)
      OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
      READ(EE_UNIT) VER,HSIZE,BSIZE
      IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
      OFFSET = HSIZE

      NS = 0
      DO WHILE(OFFSET < FSIZE)
        ISTAT = FSEEK(EE_UNIT,OFFSET,0)

        IF( EOF(EE_UNIT) )THEN
          WRITE(EE_UNIT) EETIME
          EXIT
        ENDIF
        IF( ISTAT /= 0 ) EXIT

        READ(EE_UNIT) PTIME

        IF( DEBUG )THEN
          NS=NS+1
          IF( NS == 8 )THEN
            WRITE(*,'(F10.3)')PTIME
            NS=0
          ELSE
            WRITE(*,'(F10.3,\)')PTIME
          ENDIF
        ENDIF
        IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

        OFFSET = OFFSET + BSIZE
      ENDDO
      IF( DEBUG ) WRITE(*,'(" ")')
      WRITE(*,'(A)')'FINISHED READING BED'
    ELSE
      NERR = 0
      IORIGIN = 4
      100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
      WRITE(EE_UNIT) EETIME
    ENDIF

    WRITE(EE_UNIT) (INT(KBT_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) ((REAL(HBED_Global(L,K),4), K=1,KB), L=2,LA_Global)
    WRITE(EE_UNIT) ((REAL(BDENBED_Global(L,K),4), K=1,KB), L=2,LA_Global)
    WRITE(EE_UNIT) ((REAL(PORBED_Global(L,K),4), K=1,KB), L=2,LA_Global)

    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED
        WRITE(EE_UNIT) ((REAL(SEDB_Global(L,K,NS),4), K=1,KB), L=2,LA_Global)
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NSED+NX
        WRITE(EE_UNIT) ((REAL(SNDB_Global(L,K,NX),4), K=1,KB), L=2,LA_Global)
      ENDDO
    ENDIF
    IF( ISTRAN(5) >= 1 )THEN
      DO NT=1,NTOX
        WRITE(EE_UNIT) ((REAL(TOXB_Global(L,K,NT),4), K=1,KB), L=2,LA_Global)
      ENDDO
    ENDIF
    FLUSH(EE_UNIT)
    CLOSE(EE_UNIT,STATUS='KEEP')
    
    NBEDSTEPS = 0
  ENDIF
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_BED.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF

END SUBROUTINE

INTEGER FUNCTION BLOCKBED(CELL3D)
  
  INTEGER:: DSIZE, CELL2D, CELL3D
  
  CELL2D = LA_Global - 1
  DSIZE = 2 !**EETIME
  DSIZE = DSIZE + CELL2D
  DSIZE = DSIZE + 3*KB*CELL2D
  IF( ISTRAN(6) >= 1 )THEN
    DSIZE = DSIZE + NSED*KB*CELL2D
  ENDIF
  IF( ISTRAN(7) >= 1 )THEN
    DSIZE = DSIZE + NSND*KB*CELL2D
  ENDIF
  IF( ISTRAN(5) >= 1 )THEN
    DSIZE = DSIZE + NTOX*KB*CELL2D
  ENDIF
  BLOCKBED = 4*DSIZE
  
  END FUNCTION

  SUBROUTINE SDOUT(JSEXPLORER)

  INTEGER(IK4) :: VER,HSIZE,BSIZE,ISTAT
  INTEGER(IK4) :: L,K,JSEXPLORER,ITMP,NS
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME

  NSEDSTEPS=NSEDSTEPS+1
  IF( NSEDSTEPS >= ABS(ISSDBIN) .OR. JSEXPLORER == 1 )THEN

    FILENAME=OUTDIR//'EE_SD.OUT'
    EE_UNIT = 101                      ! *** EE_SD.OUT
    IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_SD == 0 )THEN
      RSTFIRST_SD=1
      WRITE(*,'(A)')'READING TO STARTING TIME FOR SD'
      FSIZE = FILESIZE(FILENAME)
      OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
      READ(EE_UNIT) VER,HSIZE,BSIZE
      IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
      OFFSET = HSIZE

      NS = 0
      DO WHILE(OFFSET < FSIZE)
        ISTAT = FSEEK(EE_UNIT,OFFSET,0)

        IF( EOF(EE_UNIT) )THEN
          WRITE(EE_UNIT) EETIME
          EXIT
        ENDIF
        IF( ISTAT /= 0 ) EXIT

        READ(EE_UNIT) PTIME

        IF( DEBUG )THEN
          NS=NS+1
          IF( NS == 8 )THEN
            WRITE(*,'(F10.3)')PTIME
            NS=0
          ELSE
            WRITE(*,'(F10.3,\)')PTIME
          ENDIF
        ENDIF
        IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

        OFFSET = OFFSET + BSIZE
      ENDDO
      IF( DEBUG ) WRITE(*,'(" ")')
      WRITE(*,'(A)')'FINISHED READING SD'

    ELSE
      NERR = 0
      IORIGIN = 7
      100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
      WRITE(EE_UNIT) EETIME
    ENDIF

    WRITE(EE_UNIT) ((REAL(SMPON_Global(L,K),4), L=2,LA_Global), K=1,3)
    WRITE(EE_UNIT) ((REAL(SMPOP_Global(L,K),4), L=2,LA_Global), K=1,3)
    WRITE(EE_UNIT) ((REAL(SMPOC_Global(L,K),4), L=2,LA_Global), K=1,3)
    WRITE(EE_UNIT) ((REAL(SMDFN_Global(L,K),4), L=2,LA_Global), K=1,3)
    WRITE(EE_UNIT) ((REAL(SMDFP_Global(L,K),4), L=2,LA_Global), K=1,3)
    WRITE(EE_UNIT) ((REAL(SMDFC_Global(L,K),4), L=2,LA_Global), K=1,3)
    WRITE(EE_UNIT) (REAL(SM1NH4_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM2NH4_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM1NO3_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM2NO3_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM1PO4_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM2PO4_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM1H2S_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM2H2S_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM1SI_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SM2SI_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SMPSI_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SMBST_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SMT_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SMCSOD_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(SMNSOD_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(WQBFNH4_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(WQBFNO3_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(WQBFO2_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(WQBFCOD_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(WQBFPO4D_Global(L),4), L=2,LA_Global)
    WRITE(EE_UNIT) (REAL(WQBFSAD_Global(L),4), L=2,LA_Global)

    FLUSH(EE_UNIT)
    CLOSE(EE_UNIT,STATUS='KEEP')
    NSEDSTEPS=0
  ENDIF
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_SD.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF
  
  END SUBROUTINE

  SUBROUTINE RPEMOUT(JSEXPLORER)

  INTEGER(IK4) :: VER,HSIZE,BSIZE,ISTAT
  INTEGER(IK4) :: L,JSEXPLORER,NS,ITMP
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME
  LOGICAL(4)   :: LMASK

  ! *** RPEM
  IF( ISRPEM > 0 )THEN
    ! *** IF JSEXPLORER=1 THEN WRITE THE ARRAYS (I.E. IC'S)
    NRPEMSTEPS=NRPEMSTEPS+1

    IF( NRPEMSTEPS >= NRPEMEE .OR. JSEXPLORER == 1 )THEN
      FILENAME=OUTDIR//'EE_RPEM.OUT'
      EE_UNIT = 102                      ! *** EE_RPEM.OUT
      IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_RPEM == 0 )THEN
        RSTFIRST_RPEM=1
        WRITE(*,'(A)')'READING TO STARTING TIME FOR RPEM'
        FSIZE = FILESIZE(FILENAME)
        OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
        READ(EE_UNIT) VER,HSIZE,BSIZE
        IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
        OFFSET = HSIZE

        NS = 0
        DO WHILE(OFFSET < FSIZE)
          ISTAT = FSEEK(EE_UNIT,OFFSET,0)

          IF( EOF(EE_UNIT) )THEN
            WRITE(EE_UNIT) EETIME
            EXIT
          ENDIF
          IF( ISTAT /= 0 ) EXIT

          READ(EE_UNIT) PTIME

          IF( DEBUG )THEN
            NS=NS+1
            IF( NS == 8 )THEN
              WRITE(*,'(F10.3)')PTIME
              NS=0
            ELSE
              WRITE(*,'(F10.3,\)')PTIME
            ENDIF
          ENDIF
          IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

          OFFSET = OFFSET + BSIZE
        ENDDO
        IF( DEBUG ) WRITE(*,'(" ")')
        WRITE(*,'(A)')'FINISHED READING RPEM'

      ELSE
        NERR = 0
        IORIGIN = 9
        100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
        WRITE(EE_UNIT) EETIME
      ENDIF

      DO L=2,LA_Global
        IF( LMASKRPEM_Global(L) .OR. INITRPEM == 0 )THEN
          WRITE(EE_UNIT) REAL(WQRPS_Global(L),4)
        ENDIF
      ENDDO
      DO L=2,LA_Global
        IF( LMASKRPEM_Global(L) .OR. INITRPEM == 0 )THEN
          WRITE(EE_UNIT) REAL(WQRPR_Global(L),4)
        ENDIF
      ENDDO
      DO L=2,LA_Global
        IF( LMASKRPEM_Global(L) .OR. INITRPEM == 0 )THEN
          WRITE(EE_UNIT) REAL(WQRPE_Global(L),4)
        ENDIF
      ENDDO
      DO L=2,LA_Global
        IF( LMASKRPEM_Global(L) .OR. INITRPEM == 0 )THEN
          WRITE(EE_UNIT) REAL(WQRPD_Global(L),4)
        ENDIF
      ENDDO

      FLUSH(EE_UNIT)
      CLOSE(EE_UNIT,STATUS='KEEP')
      NRPEMSTEPS = 0
    ENDIF
  ENDIF
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_RPEM.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF
  
END SUBROUTINE

SUBROUTINE WQOUT

  INTEGER(IK4) :: VER,HSIZE,BSIZE,ISTAT
  INTEGER(IK4) :: L,K,MW,NW,ITMP,NS
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: WQ
  REAL(RKD)    :: PTIME

  ! *** ISWQLVL = 0
  ! ***  1) CHC - cyanobacteria
  ! ***  2) CHD - diatom algae
  ! ***  3) CHG - green algae
  ! ***  4) ROC - refractory particulate organic carbon
  ! ***  5) LOC - labile particulate organic carbon
  ! ***  6) DOC - dissolved organic carbon
  ! ***  7) ROP - refractory particulate organic phosphorus
  ! ***  8) LOP - labile particulate organic phosphorus
  ! ***  9) DOP - dissolved organic phosphorus
  ! *** 10) P4D - total phosphate
  ! *** 11) RON - refractory particulate organic nitrogen 22) macroalgae
  ! *** 12) LON - labile particulate organic nitrogen
  ! *** 13) DON - dissolved organic nitrogen
  ! *** 14) NHX - ammonia nitrogen
  ! *** 15) NOX - nitrate nitrogen
  ! *** 16) SUU - particulate biogenic silica
  ! *** 17) SAA - dissolved available silica
  ! *** 18) COD - chemical oxygen demand
  ! *** 19) DOX - dissolved oxygen
  ! *** 20) TAM - total active metal
  ! *** 21) FCB - fecal coliform bacteria

  ! *** ISWQLVL = 1
  ! ***  1) ROC - Refractory particulate organic carbon
  ! ***  2) LOC - Labile particulate organic carbon
  ! ***  3) DOC - Dissolved organic carbon
  ! ***  4) ROP - Refractory particulate organic phosphorus
  ! ***  5) LOP - Labile particulate organic phosphorus
  ! ***  6) DOP - Dissolved organic phosphorus
  ! ***  7) P4D - Total phosphate
  ! ***  8) RON - Refractory particulate organic nitrogen
  ! ***  9) LON - Labile particulate organic nitrogen
  ! *** 10) DON - Dissolved organic nitrogen
  ! *** 11) NHX - Ammonia nitrogen
  ! *** 12) NOX - Nitrate nitrogen
  ! *** 13) SUU - Particulate biogenic silica
  ! *** 14) SAA - Dissolved available silica
  ! *** 15) COD - Chemical oxygen demand
  ! *** 16) DOX - Dissolved oxygen
  ! *** 17) TAM - Total active metal
  ! *** 18) FCB - Fecal coliform bacteria
  ! *** 19) CO2 - Dissolved carbon dioxide
  ! *** 20) BIO - Any number of biota classes (NALGAE), phytoplankton, macrophytes and zooplankton
  ! ***           Can include green, blue-green, cyanobacteria, etc.  

  IF( ISRESTI == 1 )THEN
    IWQ = 0
    DO MW=1,NWQV
      IWQ(MW) = ISKINETICS(MW)
    ENDDO
  ENDIF

  FILENAME=OUTDIR//'EE_WQ.OUT'
  EE_UNIT = 100                      ! *** EE_WQ.OUT
  IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_WQ == 0 )THEN
    RSTFIRST_WQ=1
    WRITE(*,'(A)')'READING TO STARTING TIME FOR WQ'
    FSIZE = FILESIZE(FILENAME)
    OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
    READ(EE_UNIT) VER,HSIZE,BSIZE
    IF( VER /= 8400 .AND. VER /=10300 ) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
    OFFSET = HSIZE

    NS = 0
    DO WHILE(OFFSET < FSIZE)
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      IF( EOF(EE_UNIT) )THEN
        WRITE(EE_UNIT) EETIME
        EXIT
      ENDIF
      IF( ISTAT /= 0 ) EXIT

      READ(EE_UNIT) PTIME

      IF( DEBUG )THEN
        NS=NS+1
        IF( NS == 8 )THEN
          WRITE(*,'(F10.3)')PTIME
          NS=0
        ELSE
          WRITE(*,'(F10.3,\)')PTIME
        ENDIF
      ENDIF
      IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

      OFFSET = OFFSET + BSIZE
    ENDDO
    IF( DEBUG ) WRITE(*,'(" ")')
    WRITE(*,'(A)')'FINISHED READING WQ'

  ELSE
    NERR = 0
    IORIGIN = 6
    100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
    WRITE(EE_UNIT) EETIME
  ENDIF

  IF( ISWQLVL == 0 )THEN
    ! *** Starts from ROC (4 - 22) as algae are moved to the end
    DO NW = 4,22 
      IF( IWQ(NW) > 0 )THEN
        DO L = 2,LA_Global
          DO K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            WRITE(EE_UNIT) WQ
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    
    ! *** Algae (1,2,3) are now moved to the end
    DO NW = 1,3
      IF( IWQ(NW) > 0 )THEN
        DO L = 2,LA_Global
          DO K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            WRITE(EE_UNIT) WQ
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    
    ! *** 23) macroalgae
    IF( nfixed > 0 )THEN
      NW = 23 
      IF( IWQ(NW) > 0 )THEN
        DO L = 2,LA_Global
          DO K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            WRITE(EE_UNIT) WQ
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    ! *** Zooplankton
    IF( IWQZPL > 0 )THEN
      DO NW = 1,NZOOPL
        DO L = 2,LA_Global
          DO K = KSZ_Global(L),KC
            IF( NFIXED > 0 )THEN
              WQ = WQV_Global(L,K,22+NW)
            ELSE
              WQ = WQV_Global(L,K,21+NW)
            ENDIF
            WRITE(EE_UNIT) WQ
          ENDDO
        ENDDO
      ENDDO
    ENDIF  
    
  ELSE
    ! *** New DSI standard WQ
    DO NW = 1,NWQV
      IF( IWQ(NW) > 0 )THEN
        DO L = 2,LA_Global
          DO K = KSZ_Global(L),KC
            WQ = WQV_Global(L,K,NW)
            WRITE(EE_UNIT) WQ
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  FLUSH(EE_UNIT)
  CLOSE(EE_UNIT,STATUS='KEEP')

  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_WQ.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF
  
END SUBROUTINE

SUBROUTINE BCOUT
  ! *** OUTPUT QSUM AND OPEN BC FLOWS
  INTEGER(IK4) :: VER,HSIZE,BSIZE,ISTAT
  INTEGER(IK4) :: L,K,IBC,LL,LQ,LE,LN,ITMP,NS
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME

  FILENAME=OUTDIR//'EE_BC.OUT'
  EE_UNIT = 104                      ! *** EE_BC.OUT

  IF( ISRESTI /= 0 .AND. ICONTINUE == 1 .AND. RSTFIRST_BC == 0 )THEN
    RSTFIRST_BC=1
    WRITE(*,'(A)')'READING TO STARTING TIME FOR BC'
    FSIZE = FILESIZE(FILENAME)
    OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)
    READ(EE_UNIT) VER,HSIZE,BSIZE
    IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'
    OFFSET = HSIZE

    NS = 0
    DO WHILE(OFFSET < FSIZE)
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      IF( EOF(EE_UNIT) )THEN
        WRITE(EE_UNIT) EETIME
        EXIT
      ENDIF
      IF( ISTAT /= 0 ) EXIT

      READ(EE_UNIT) PTIME

      IF( DEBUG )THEN
        NS=NS+1
        IF( NS == 8 )THEN
          WRITE(*,'(F10.3)')PTIME
          NS=0
        ELSE
          WRITE(*,'(F10.3,\)')PTIME
        ENDIF
      ENDIF
      IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

      OFFSET = OFFSET + BSIZE
    ENDDO
    IF( DEBUG ) WRITE(*,'(" ")')
    WRITE(*,'(A)')'FINSIHED READING BC'

  ELSE
    NERR = 0
    IORIGIN = 0
    100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
    WRITE(EE_UNIT) EETIME
  ENDIF

  IF( ISRESTI == 0 .OR. ICONTINUE == 0) NRESTART=0

  ! *** GET AVERAGE FLOWS FOR SNAPSHOT
  ! @todo put bcout_average routine back in and make it worth with MPI
  !CALL BCOUT_AVERAGE

  ! *** OUTPUT SELECTIVE QSUM
  IF( KC > 1 )THEN

    WRITE(EE_UNIT) (REAL(QSUM_Global(L,KC),4), L=2,LA_Global)

    DO IBC=1, NBCS
      L=LBCS(IBC) !
      WRITE(EE_UNIT) (REAL(QSUM_Global(L,K),4),K=1,KC)
    ENDDO
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      WRITE(EE_UNIT) (REAL(QSUM_Global(L, KSZ_Global(L)),4), L=2,LA_Global)
    ENDIF
  ELSE
    ! *** SINGLE LAYER
    DO L=2,LA_Global
      WRITE(EE_UNIT) REAL(QSUM_Global(L,1),4)
    ENDDO
  ENDIF

  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  DO LL=1, NPBS
    LQ=LPBS(LL)
    LN=LNC_Global(LQ)
    WRITE(EE_UNIT) (REAL(VHDX2_Global(LN,K),4),K=1,KC)
  ENDDO
  DO LL=1, NPBW
    LQ=LPBW(LL)
    LE=LEC_Global(LQ)
    WRITE(EE_UNIT) (REAL(UHDY2_Global(LE,K),4),K=1,KC)
  ENDDO
  DO LL=1, NPBE
    LQ=LPBE(LL)
    WRITE(EE_UNIT) (REAL(UHDY2_Global(LQ,K),4),K=1,KC)
  ENDDO
  DO LL=1, NPBN
    LQ=LPBN(LL)
    WRITE(EE_UNIT) (REAL(VHDX2_Global(LQ,K),4),K=1,KC)
  ENDDO

  ! @todo verify this will work for MPI
  ! *** SAVE FOR INFORMATION ONLY.   QCTLT AND HSCTL.CUR.FLOW ALREADY INCLUDED IN QSUM
  IF( NQCTL > 0 )THEN

    WRITE(EE_UNIT) (((REAL(QCTLT(K,L,NS),4), K=1,KC), NS=1,2), L=1,NQCTL) !< do not need a global one since this was not decomposed

    IF(NQCTLSER > 0 .OR. NQCTRULES > 0 )THEN
      WRITE(EE_UNIT) (INT(HSCTL(L).CUR.STATE), L=1,NQCTL)
      WRITE(EE_UNIT) (INT(1), L=1,NQCTL)                     ! HSCTL(L).CUR.UNITS not used
      WRITE(EE_UNIT) (INT(HSCTL(L).CUR.ID), L=1,NQCTL)
      WRITE(EE_UNIT) (REAL(HSCTL(L).CUR.FLOW,4), L=1,NQCTL)
      WRITE(EE_UNIT) (REAL(HSCTL(L).CUR.HEIGHT,4), L=1,NQCTL)
      WRITE(EE_UNIT) (REAL(HSCTL(L).CUR.WIDTH,4), L=1,NQCTL)
      WRITE(EE_UNIT) (REAL(HSCTL(L).CUR.SILL,4), L=1,NQCTL)
    ENDIF
  ENDIF

  ! @todo verify this will work for MPI
  ! *** SAVE FOR INFORMATION ONLY.   WRCTL ALREADY INCLUDED IN QSUM
  IF( NQWR > 0 )THEN
    WRITE(EE_UNIT) (INT(WRCTL(L).CUR.STATE),   L=1,NQWR)
    WRITE(EE_UNIT) (REAL(WRCTL(L).CUR.FLOW,4), L=1,NQWR)
  ENDIF

  CLOSE(EE_UNIT,STATUS='KEEP')

  !SUMTIME = 0.
  !BCQSUM = 0.
  !BCQSUME = 0.
  !BCUHDY2 = 0.
  !BCVHDX2 = 0
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_BC.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF

  END SUBROUTINE

  SUBROUTINE BCOUT_ACCUMULATE(DELT)
  ! *** ACCUMULATE BOUNDARY INFLOWS/OUTFLOWS
  REAL,INTENT(IN) :: DELT
  INTEGER(IK4) :: L,K,LL,LQ,LE,LN,IBC

  IF( INITBCOUT == 0 ) CALL BCOUT_INITIALIZE

  SUMTIME = SUMTIME + DELT

  ! *** SUM SELECTIVE QSUM
  IF( KC > 1 )THEN
    DO L=2,LA_Global
      BCQSUM(L,KC) = BCQSUM(L,KC) + QSUM(L,KC)*DELT
    ENDDO
    DO IBC=1,NBCS
      L=LBCS(IBC)
      DO K=1,KS
        BCQSUM(L,K) = BCQSUM(L,K) + QSUM(L,K)*DELT
      ENDDO
    ENDDO
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=2,LA_Global
        K = KSZ_Global(L)
        BCQSUM(L,K) = BCQSUM(L,K) + QSUM(L,K)*DELT
      ENDDO
    ENDIF
  ELSE
    ! *** SINGLE LAYER
    DO L=2,LA_Global
      BCQSUME(L) = BCQSUME(L) + QSUME(L)*DELT
    ENDDO
  ENDIF

  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  DO LL=1,NPBS
    LQ=LPBS(LL)
    LN=LNC_Global(LQ)
    DO K=1,KC
      BCVHDX2(LN,K) = BCVHDX2(LN,K) + VHDX2(LN,K)*DELT
    ENDDO
  ENDDO
  DO LL=1,NPBW
    LQ=LPBW(LL)
    LE=LEC_Global(LQ)
    DO K=1,KC
      BCUHDY2(LE,K) = BCUHDY2(LE,K) + UHDY2(LE,K)*DELT
    ENDDO
  ENDDO
  DO LL=1,NPBE
    LQ=LPBE(LL)
    DO K=1,KC
      BCUHDY2(LQ,K) = BCUHDY2(LQ,K) + UHDY2(LQ,K)*DELT
    ENDDO
  ENDDO
  DO LL=1,NPBN
    LQ=LPBN(LL)
    DO K=1,KC
      BCVHDX2(LQ,K) = BCVHDX2(LQ,K) + VHDX2(LQ,K)*DELT
    ENDDO
  ENDDO

  RETURN

  END SUBROUTINE

  SUBROUTINE BCOUT_AVERAGE
  ! *** COMPUTE AVEARGE BOUNDARY INFLOWS/OUTFLOWS

  INTEGER(IK4) :: L,K,LL,LQ,LE,LN,IBC

  IF( SUMTIME <= 0. ) RETURN

  ! *** SUM SELECTIVE QSUM
  IF( KC > 1 )THEN
    DO L=2,LA_Global
      BCQSUM(L,KC) = BCQSUM(L,KC)/SUMTIME
    ENDDO
    DO IBC=1,NBCS
      L=LBCS(IBC)
      DO K=1,KS
        BCQSUM(L,K) = BCQSUM(L,K)/SUMTIME
      ENDDO
    ENDDO
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=2,LA_Global
        K = KSZ_Global(L)
        BCQSUM(L,K) = BCQSUM(L,K)/SUMTIME
      ENDDO
    ENDIF
  ELSE
    ! *** SINGLE LAYER
    DO L=2,LA_Global
      BCQSUME(L) = BCQSUME(L)/SUMTIME
    ENDDO
  ENDIF

  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  DO LL=1,NPBS
    LQ=LPBS(LL)
    LN=LNC_Global(LQ)
    DO K=1,KC
      BCVHDX2(LN,K) = BCVHDX2(LN,K)/SUMTIME
    ENDDO
  ENDDO
  DO LL=1,NPBW
    LQ=LPBW(LL)
    LE=LEC_Global(LQ)
    DO K=1,KC
      BCUHDY2(LE,K) = BCUHDY2(LE,K)/SUMTIME
    ENDDO
  ENDDO
  DO LL=1,NPBE
    LQ=LPBE(LL)
    DO K=1,KC
      BCUHDY2(LQ,K) = BCUHDY2(LQ,K)/SUMTIME
    ENDDO
  ENDDO
  DO LL=1,NPBN
    LQ=LPBN(LL)
    DO K=1,KC
      BCVHDX2(LQ,K) = BCVHDX2(LQ,K)/SUMTIME
    ENDDO
  ENDDO

  RETURN

  END SUBROUTINE

  SUBROUTINE BCOUT_INITIALIZE

  ! *** SETUP ACCUMULATION ARRAYS - @todo specify what the hell each is accumulating so a human can read...

  ALLOCATE(BCQSUM(LCM_Global,  KCM))  ! *** Modified so LCM is global
  ALLOCATE(BCUHDY2(LCM_Global, KCM))  ! *** Modified so LCM is global
  ALLOCATE(BCVHDX2(LCM_Global, KCM))  ! *** Modified so LCM is global
  ALLOCATE(BCQSUME(LCM_Global))       ! *** Modified so LCM is global

  SUMTIME = 0.
  BCQSUM = 0.
  BCUHDY2 = 0.
  BCVHDX2 = 0.
  BCQSUME = 0.

  INITBCOUT = 1

  END SUBROUTINE

  SUBROUTINE BCOUT_ZERO
  ! *** COMPUTE AVEARGE BOUNDARY INFLOWS/OUTFLOWS

  INTEGER(IK4) :: L,K,LL,LQ,LE,LN,IBC

  SUMTIME = 0.

  ! *** SUM SELECTIVE QSUM
  IF( KC > 1 )THEN
    DO L=2,LA_Global
      BCQSUM(L,KC) = 0.
    ENDDO
    DO IBC=1,NBCS
      L=LBCS(IBC)
      DO K=1,KS
        BCQSUM(L,K) = 0.
      ENDDO
    ENDDO
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DO L=2,LA_Global
        K = KSZ_Global(L)
        BCQSUM(L,K) = 0.
      ENDDO
    ENDIF
  ELSE
    ! *** SINGLE LAYER
    BCQSUME = 0.
  ENDIF

  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  DO LL=1,NPBS
    LQ=LPBS(LL)
    LN=LNC_Global(LQ)
    DO K=1,KC
      BCVHDX2(LN,K) = 0.
    ENDDO
  ENDDO
  DO LL=1,NPBW
    LQ=LPBW(LL)
    LE=LEC_Global(LQ)
    DO K=1,KC
      BCUHDY2(LE,K) = 0.
    ENDDO
  ENDDO
  DO LL=1,NPBE
    LQ=LPBE(LL)
    DO K=1,KC
      BCUHDY2(LQ,K) = 0.
    ENDDO
  ENDDO
  DO LL=1,NPBN
    LQ=LPBN(LL)
    DO K=1,KC
      BCVHDX2(LQ,K) = 0.
    ENDDO
  ENDDO

END SUBROUTINE

INTEGER FUNCTION BLOCKBC(CELL3D)
  INTEGER :: DSIZE,CELL2D,CELL3D
  CELL2D = LA_Global - 1
  DSIZE = 2                 ! *** EETIME
  IF( KC > 1 )THEN
    DSIZE = DSIZE + CELL2D
    DSIZE = DSIZE + KC*NBCS
    IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
      DSIZE = DSIZE + CELL2D
    ENDIF
  ELSE
    DSIZE = DSIZE + CELL2D  ! *** SINGLE LAYER
  ENDIF

  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
  DSIZE = DSIZE + KC*NPBS
  DSIZE = DSIZE + KC*NPBW
  DSIZE = DSIZE + KC*NPBE
  DSIZE = DSIZE + KC*NPBN
  IF( NQCTL > 0 )THEN
    DSIZE = DSIZE + NQCTL*2*KC
    IF(NQCTLSER > 0 .OR. NQCTRULES > 0 )THEN
      DSIZE = DSIZE + 7*NQCTL
    ENDIF
  ENDIF
  IF( NQWR > 0 )THEN
    DSIZE = DSIZE + 2*NQWR
  ENDIF
  BLOCKBC = DSIZE*4
END FUNCTION

SUBROUTINE SHELLFISHOUT()

  ! ** SHELLFISH OUTPUT
  INTEGER(IK4) :: VER,HSIZE,BSIZE
  INTEGER(IK4) :: I,L,K,ITMP,NN,NS,ISTAT
  INTEGER(IK8) :: FSIZE, OFFSET
  REAL(RK4)    :: TMP
  REAL(RKD)    :: PTIME

  FILENAME=OUTDIR//'EE_SHF.OUT'
  EE_UNIT = 103                      ! *** EE_SHF.OUT
  IF( ISRESTI /= 0 .AND. ICONTINUE == 1 )THEN
    WRITE(*,'(A)')'READING TO STARTING TIME FOR SHELLFISH'
    FSIZE = FILESIZE(FILENAME)
    OPEN(EE_UNIT,FILE=FILENAME,ACTION='READWRITE',STATUS='OLD',FORM='BINARY',SHARED)     
    READ(EE_UNIT) VER,HSIZE,BSIZE
    IF(VER /= 8400) WRITE(*,*)'FILE IS CORRUPTED OR VERSION IS INVALID!'    
    OFFSET = HSIZE
    
    NS = 0
    DO WHILE(OFFSET < FSIZE)
      ISTAT = FSEEK(EE_UNIT,OFFSET,0)

      IF( EOF(EE_UNIT) )THEN
        WRITE(EE_UNIT) EETIME
        EXIT
      ENDIF
      IF( ISTAT /= 0 ) EXIT
        
      READ(EE_UNIT) PTIME
      
      IF( DEBUG )THEN
        NS=NS+1
        IF( NS == 8 )THEN
          WRITE(*,'(F10.3)')PTIME
          NS=0
        ELSE
          WRITE(*,'(F10.3,\)')PTIME
        ENDIF
      ENDIF
      IF( ABS(PTIME-TIMEDAY) <= 1E-4 .OR. PTIME > TIMEDAY ) EXIT

      OFFSET = OFFSET + BSIZE
    ENDDO
    IF( DEBUG ) WRITE(*,'(" ")')
    WRITE(*,'(A)')'FINISHED READING SHELLFISH'
  ELSE
    NERR = 0
    IORIGIN = 8
    100 OPEN(EE_UNIT,FILE=FILENAME,STATUS='OLD',POSITION='APPEND',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)
    WRITE(EE_UNIT) EETIME
  ENDIF
  
  IF( ISRESTI == 0 .OR. ICONTINUE == 0) NRESTART=0 
  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).CELL_C(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI_C(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).SHLEN(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI_D(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).INDI_H(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
    DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).HARV_C(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).NP(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
    DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).PR(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).CRB(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO
    DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) (REAL(SF(I).SPAWN(L,K),4), K=KSZ(L),KC)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).WQAQ(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).WQEA(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).FR(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).BMG(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_RESPI(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_GRAZI(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_DEATH(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_FECAL(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_URINE(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_RPOC(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_LPOC(L),4)
    ENDDO
  ENDDO  
  DO I=1,NSF
    DO L=2,LA
      IF(FARMCELL(L) > 0) WRITE(EE_UNIT) REAL(SF(I).B_DOC(L),4)
    ENDDO
  ENDDO

  FLUSH(EE_UNIT)
  CLOSE(EE_UNIT,STATUS='KEEP') 
  
  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_SHF.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF
  
END SUBROUTINE

SUBROUTINE GRIDOUT
  USE INFOMOD,ONLY:READSTR
  IMPLICIT NONE
  INTEGER(IK4) :: NACTIVE,VER,HSIZE,BSIZE,K,L,CELL3D

  NACTIVE=LA_Global-1

  FILENAME=OUTDIR//'EE_GRD.OUT'
  VER   = 8400
  HSIZE = 8*4
  BSIZE = 0
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN')
  CLOSE(EE_UNIT,STATUS='DELETE')
  OPEN(EE_UNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY',SHARED)
  WRITE(EE_UNIT) VER,HSIZE,BSIZE
  WRITE(EE_UNIT) INT(IC,4),INT(JC,4),INT(KC,4),NACTIVE,IGRIDV

  !XCR(1:4,JMN:JMX,IMN:IMX),YCR(1:4,JMN:JMX,IMN:IMX)
  IF( IGRIDV > 0 )THEN ! Sigma-Zed level
    WRITE(EE_UNIT) (INT(KSZ_Global(L),4), L = 2,LA_Global)
    WRITE(EE_UNIT) ((REAL(DZC(L,K),4), K=1,KC), L = 2,LA_Global)
  ELSE  ! ** Standard sigma stretched level
    WRITE(EE_UNIT) (REAL(DZCK(K),4), K=1,KC)
  ENDIF
  CLOSE(EE_UNIT,STATUS='KEEP')

END SUBROUTINE

SUBROUTINE ARRAYSOUT

  ! *** TIME VARIABLE USER SPECIFIED ARRAYS IN THIS SECTION

  ! FLAGS: ARRAY TYPE, TIME VARIABLE
  ! ARRAY TYPE:    0 = L            DIM'D
  !                1 = L,KC         DIM'D
  !                2 = L,0:KC       DIM'D
  !                3 = L,KB         DIM'D
  !                4 = L,KC,NCLASS  DIM'D
  ! TIME VARIABLE: 0 = NOT CHANGING
  !                1 = TIME VARYING

  INTEGER(IK4) :: K,L,ITYPE,ITIMEVAR
  REAL(RK4)    :: ZERO
  REAL(RK4)    :: TMPVAL
  CHARACTER*8 ARRAYNAME

  IF( ISINWV == 2 )THEN
    ZERO=0.0

    FILENAME=OUTDIR//'EE_ARRAYS.OUT'
    EE_UNIT = 105                      ! *** EE_ARRAYS.OUT

    NERR = 0
    IORIGIN = 10
    100 OPEN(EE_UNIT,FILE=FILENAME,POSITION='APPEND',STATUS='UNKNOWN',FORM='BINARY',SHARED,ERR=999,IOSTAT=IERRio)

    ! *** TIME VARIABLE FLAG
    ITIMEVAR = 1
    
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='TBX'
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      TMPVAL=TBX(L)
      WRITE(EE_UNIT)TMPVAL
    ENDDO
    
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='TBY'
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      TMPVAL=TBY(L)
      WRITE(EE_UNIT)TMPVAL
    ENDDO
        
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='FUHDYE' 
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT)REAL(FUHDYE(L),4)
    ENDDO

    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='FVHDXE' 
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT)REAL(FVHDXE(L),4)
    ENDDO
    
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='P'
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT)REAL(PMCTESTX(1,L),4)
    ENDDO
    
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='PCG'
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT)REAL(PMCTESTX(2,L),4)
    ENDDO
    
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='RCG'
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT)REAL(PMCTESTX(3,L),4)
    ENDDO
    
    ITYPE = 0
    WRITE(EE_UNIT)ITYPE,ITIMEVAR
    ARRAYNAME='TMPCG'
    WRITE(EE_UNIT)ARRAYNAME
    DO L=2,LA_Global
      WRITE(EE_UNIT)REAL(PMCTESTX(4,L),4)
    ENDDO
    
    IF( ISHDMF >= 1 )THEN
      ITYPE = 1
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='FMDUX'
      WRITE(EE_UNIT)ARRAYNAME
      DO K=1,KC
        DO L=2,LA_Global
          TMPVAL=FMDUX(L,K)
          WRITE(EE_UNIT)TMPVAL
        ENDDO
      ENDDO

      ITYPE = 1
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='FMDUY'
      WRITE(EE_UNIT)ARRAYNAME
      DO K=1,KC
        DO L=2,LA_Global
          TMPVAL=FMDUY(L,K)
          WRITE(EE_UNIT)TMPVAL
        ENDDO
      ENDDO

      ITYPE = 1
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='FMDVX'
      WRITE(EE_UNIT)ARRAYNAME
      DO K=1,KC
        DO L=2,LA_Global
          TMPVAL=FMDVX(L,K)
          WRITE(EE_UNIT)TMPVAL
        ENDDO
      ENDDO

      ITYPE = 1
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='FMDVY'
      WRITE(EE_UNIT)ARRAYNAME
      DO K=1,KC
        DO L=2,LA_Global
          TMPVAL=FMDVY(L,K)
          WRITE(EE_UNIT)TMPVAL
        ENDDO
      ENDDO
    ENDIF

   
    
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='HU'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(HU(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='STBY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(STBY(L),4)
    !ENDDO



    if( .false. )then
      ! *** SURFACE HEAT FLUX ARRAYS
      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='HBLW'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(PMCTESTX(1,L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='HBCD'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(PMCTESTX(2,L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='HBCV'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(PMCTESTX(3,L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='HBEV'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(PMCTESTX(4,L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='HBNT'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(PMCTESTX(5,L),4)
      ENDDO
      ! *** END OF SURFACE HEAT FLUX ARRAYS

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='WINDST'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(WINDST(L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='WINDCD10'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(WINDCD10(L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='TSX'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(TSX(L),4)
      ENDDO

      ITYPE = 0
      WRITE(EE_UNIT)ITYPE,ITIMEVAR
      ARRAYNAME='TSY'
      WRITE(EE_UNIT)ARRAYNAME
      DO L=2,LA_Global
        WRITE(EE_UNIT)REAL(TSY(L),4)
      ENDDO
    endif

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='U'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(U(L,K),4)
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='V'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(V(L,K),4)
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='TBX'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(TBX(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='TBY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(TBY(L),4)
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='B'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(B(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='U1'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(U1(L,K),4)
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='V1'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(V1(L,K),4)
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='V1U'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(V1U(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='U1V'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(U1V(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='H1U'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(H1U(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='H1V'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(H1V(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='STBX'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(STBX(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='STBY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(STBY(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='LMASKDRY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  IF( LMASKDRY(L) )THEN
    !    WRITE(EE_UNIT)REAL(1.0,4)
    !  ELSE
    !    WRITE(EE_UNIT)REAL(0.0,4)
    !  ENDIF
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='TBX'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(TBX(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='TBY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(TBY(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='QSBDLDX'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(QSBDLDX(L,2),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='QSBDLDY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(QSBDLDY(L,2),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='UBL'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(UBL(L,2),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='VBL'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(VBL(L,2),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='BLDELTA'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL( ( QSBDLDX(L,2)-QSBDLDX(LEC_Global(L),2) + QSBDLDY(L,2)-QSBDLDY(LNC(L),2) ),4)
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FXMHK'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    IF( K < KSZ_Global(L) )THEN
    !      WRITE(EE_UNIT)REAL(-999,4)
    !    ELSE
    !      WRITE(EE_UNIT)REAL(FXMHK(L,k),4)
    !    ENDIF
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FYMHK'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    IF( K < KSZ_Global(L) )THEN
    !      WRITE(EE_UNIT)REAL(-999.,4)
    !    ELSE
    !      WRITE(EE_UNIT)REAL(FYMHK(L,K),4)
    !    ENDIF
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FXMHKE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(FXMHKE(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FYMHKE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(FYMHKE(L),4)
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FXSUP'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    IF( K < KSZ_Global(L) )THEN
    !      WRITE(EE_UNIT)REAL(-999,4)
    !    ELSE
    !      WRITE(EE_UNIT)REAL(FXSUP(L,k),4)
    !    ENDIF
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FYSUP'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    IF( K < KSZ_Global(L) )THEN
    !      WRITE(EE_UNIT)REAL(-999.,4)
    !    ELSE
    !      WRITE(EE_UNIT)REAL(FYSUP(L,K),4)
    !    ENDIF
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FXSUPE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(FXSUPE(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FYSUPE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(FYSUPE(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='ET'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(PMCTESTX(1,L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='CSHE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(PMCTESTX(2,L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='ICETEMP'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(PMCTESTX(3,L),4)
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FRAZILICE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(FRAZILICE(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='AB'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    TMPVAL=AB(L,K)*HP(L)
    !    WRITE(EE_UNIT)TMPVAL
    !  ENDDO
    !ENDDO

    !ITYPE = 2
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='W'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=0,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(W(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 2
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='W2'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=0,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(W2(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='DML'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(DML(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='UUU'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(UUU(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='VVV'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(VVV(L,K),4)
    !  ENDDO
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='UV'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(UV(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='VU'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(VU(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='RCX'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(RCX(L),4)
    !ENDDO

    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='RCY'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(RCY(L),4)
    !ENDDO

    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FXVEG'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(FXVEG(L,K),4)
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 1
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FYVEG'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO K=1,KC
    !  DO L=2,LA_Global
    !    WRITE(EE_UNIT)REAL(FYVEG(L,K),4)
    !  ENDDO
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FXVEGE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(FXVEGE(L),4)
    !ENDDO
    !
    !ITYPE = 0
    !WRITE(EE_UNIT)ITYPE,ITIMEVAR
    !ARRAYNAME='FYVEGE'
    !WRITE(EE_UNIT)ARRAYNAME
    !DO L=2,LA_Global
    !  WRITE(EE_UNIT)REAL(FYVEGE(L),4)
    !ENDDO

    ! FLUSH(EE_UNIT)
    CLOSE(EE_UNIT,STATUS='KEEP')

  ENDIF

  RETURN
  
999 IF( IERRio == 30 .AND. NERR < 3 )THEN
      ! *** Unit EE_UNIT busy.  Pause for 1 second
      NERR = NERR + 1
      PRINT *, 'File in use (EE_ARRAYS.OUT) error.  Try # ',nerr
      CALL SLEEP(5)
      IF( NERR == 3 )THEN
        CLOSE(EE_UNIT,STATUS='KEEP')
      ENDIF
      GOTO 100
    ELSE
      CALL FILE_ERROR
      CALL STOPP('.')
    ENDIF

END SUBROUTINE

INTEGER(IK8) FUNCTION FILESIZE(FNAME)

  USE IFPORT

  CHARACTER*(*), INTENT(IN) :: FNAME
  INTEGER(IK4) :: ISTAT
  INTEGER(IK8) :: FINFO(12)

  ISTAT = STAT(FNAME,FINFO)     !  ISTAT = FSTAT(IUNIT,FINFO)
  IF( ISTAT == 0 )THEN
    FILESIZE = FINFO(8)
  ELSE
    FILESIZE = 0
  ENDIF

  END FUNCTION
  
  SUBROUTINE FILE_ERROR
    
    SELECT CASE ( IORIGIN )
    CASE (0)
      WRITE(*,'(a)') ' Error tying to open EE_BC.OUT'

    CASE (1)
      WRITE(*,'(a)') ' Error tying to open EE_WS.OUT'

    CASE (2)
      WRITE(*,'(a)') ' Error tying to open EE_VEL.OUT'

    CASE (3)
      WRITE(*,'(a)') ' Error tying to open EE_WC.OUT'

    CASE (4)
      WRITE(*,'(a)') ' Error tying to open EE_BED.OUT'

    CASE (5)
      WRITE(*,'(a)') ' Error tying to open EE_SEDZLJ.OUT'

    CASE (6)
      WRITE(*,'(a)') ' Error tying to open EE_WQ.OUT'

    CASE (7)
      WRITE(*,'(a)') ' Error tying to open EE_SD.OUT'

    CASE (8)
      WRITE(*,'(a)') ' Error tying to open EE_SHF.OUT'

    CASE (9)
      WRITE(*,'(a)') ' Error tying to open EE_RPEM.OUT'

    CASE (10)
      WRITE(*,'(a)') ' Error tying to open EE_ARRAYS.OUT'

    CASE DEFAULT
      WRITE(*,'(a)') ' Error tying to open EE linkage files'
    
    END SELECT  
    
  END SUBROUTINE FILE_ERROR

  END MODULE
