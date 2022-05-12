! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!< @author Paul Craig
!< @date 2015-06 - Implemented SIGMA-Z (SGZ) IN EE7.3
!< @date 2022-01 - Added mass erosion class for propwash
!< @details Scans the input files and calls routines to allocate and zero out arrays
!< @details   does some other setup as well.

SUBROUTINE VARINIT

  USE GLOBAL
  Use Allocate_Initialize
  USE SCANINPMOD
  USE Variables_MPI
  USE MPI
  USE Broadcast_Routines
  
  IMPLICIT NONE

  INTEGER :: NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7
  INTEGER :: MW,NS,M,MD,L
  CHARACTER*3 :: NSTR

  Integer :: ierr ! MPI error return value

  KPCM     = 1
  MDVSM    = 1
  MTVSM    = 1
  NDDAM    = 1
  NDQCLT   = 1
  NDQCLT2  = 1
  NDVEGSER = 1
  NGLM     = 1
  NJPSM    = 1
  NJUNXM   = 1
  NJUNYM   = 1
  NLDAM    = 1
  NPDM     = 1
  NSMGM    = 1
  NSMTSM   = 1
  NSMZM    = 1
  NTSM     = 1
  NTSSMVM  = 1
  NVEGSERM = 1
  NVEGTPM  = 100
  NWGGM    = 1
  NWQPSM   = 1
  NWQPSRM  = 1
  NWQZM    = 1
  NXYSDATM = 1
  NWNDMAP  = 1
  NATMMAP  = 1
  NICEMAP  = 1

  Call AllocateDSI( NCSER,  8,  0)
  
#ifdef _MPI
  ! *** Scanning DECOMP file to setup variables for allocation
  !Call Scan_Decomp
  Call Scan_JSON_Decomp
  ! *** After we have the partitioning from the DECOMP file setup the MPI topology
  Call Setup_MPI_Topology
#endif

  ! *** SCAN MAIN CONTROL FILE - only on the master
  CALL SCANEFDC(NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7)

  ! ***  Read domain decomposition settings
  Call Read_JSON_Decomp
  
  ! *** Scans the cell.inp and sets up some variables used for allocations later on
  Call Scan_Cell

#ifdef _MPI 
  ! *** Send out the global values originally read in
  Call Broadcast_Scalar(ICM_Global, master_id)
  Call Broadcast_Scalar(JCM_Global, master_id)
  Call Broadcast_Scalar(LCM_Global, master_id)
  Call Broadcast_Scalar(LC_Global,  master_id)

  ! *** Set some variables for array allocation
  max_width_x = max_width_x + 4  ! Adding the plus 4 because of the 2 ghost rows on each side of domain
  max_width_y = max_width_y + 4

  global_max_width_x = icm + 4   ! So this is IC + 1 + 4
  global_max_width_y = jcm + 4   ! So this is IC + 1 + 4

  ! *** Remapping ICM and JCM to local values
  icm = max_width_x
  jcm = max_width_y

  Call WriteBreak(mpi_log_unit)
  write(mpi_log_unit, '(a)')    'Writing from VARINIT afer Scan_Cell routine '
  write(mpi_log_unit, '(a,i10)') 'Local ICM     = ', ICM
  write(mpi_log_unit, '(a,i10)') 'Local JCM     = ', JCM
  write(mpi_log_unit, '(a,i10)') 'Max local LCM = ', LCM
  write(mpi_log_unit, '(a,i10)') 'LCM Global    = ', LCM_Global
  write(mpi_log_unit, '(a,i10)') 'Global max X  = ', global_max_width_x
  write(mpi_log_unit, '(a,i10)') 'Global max Y  = ', global_max_width_y
  write(mpi_log_unit, '(a,i10)') 'Max width X   = ', max_width_x
  write(mpi_log_unit, '(a,i10)') 'Max width Y   = ', max_width_y
  write(mpi_log_unit,'(a)')     'Max X/Y widths should be equal to ICM/JCM at this point '
  Call WriteBreak(mpi_log_unit)

#else
  print *,'*****************Not performing MPI Calculation Routine*****************'
  global_max_width_x = icm
  global_max_width_y = jcm
  max_width_x = global_max_width_x
  max_width_y = global_max_width_y
  n_x_partitions = 1
  n_y_partitions = 1
  active_domains = 1
#endif

  ! *** Scan other files for allocation of variables
  IF( ISCHAN > 0  ) CALL SCANMODC
  IF( ISGWIT == 2 ) CALL SCANGWSR
  IF( NASER  > 0  ) CALL SCANASER
  IF( NCSER1 > 0 .AND. ISTRAN(1) >= 1 ) CALL SCANSSER
  IF( (NCSER2 > 0 .AND. ISTRAN(2) >= 1) .OR. (ISTRAN(5) > 0 .AND. ITOXTEMP > 1 ) ) CALL SCANTSER
  IF( NCSER3 > 0  .AND. ISTRAN(3) >= 1 ) CALL SCANDSER
  IF( NCSER4 > 0  .AND. ISTRAN(4) >= 1 ) CALL SCANSFSR
  IF( NQSER  > 0  ) CALL SCANQSER
  IF( NPSER  > 0  ) CALL SCANPSER
  IF( NWSER  > 0  ) CALL SCANWSER
  IF( NQCTLT >= 1 ) CALL SCANQCTL
  IF( NQWRSR > 0  ) CALL SCANQWSER
  IF( LSEDZLJ ) CALL SCANSEDZLJ
  IF( NCSER5 > 0 .OR. NCSER6 > 0 .OR. NCSER7 > 0 ) CALL SCNTXSED
  IF( ISPROPWASH > 0 ) CALL SCANPROPWASH

  Call MPI_Barrier(comm_2d, ierr)

  ! *** These need to be broadcast to allocate the following variables correctly
#ifdef _MPI
  Call Broadcast_Scalar(NDYM, master_id)
  Call Broadcast_Scalar(NSCM, master_id)
  Call Broadcast_Scalar(NSNM, master_id)
  Call Broadcast_Scalar(NTXM, master_id)
#endif

  NSND2 = NSND     ! *** Propwash does not include non-cohesives for now.
  NSNM2 = NSNM
  NSED2 = NSED
  NSCM2 = NSCM

  ! *** Time series pointers for different constituents
  Call AllocateDSI( MSVDYE,  NDYM,   0)
  Call AllocateDSI( MSVSED,  NSCM,   0)
  Call AllocateDSI( MSVSND,  NSNM,   0)
  Call AllocateDSI( MSVTOX,  NTXM,   0)

  M = 2             ! *** M=1 (SALINITY), M=2 (TEMPERATURE)
  DO NS=1,NDYM      ! *** MUST USE NDYM INSTEAD OF NDYE TO ENSURE THE BASE OF 4 FOR LEGACY MASS BALANCE ROUTINES
    M = M + 1
    MSVDYE(NS) = M
  ENDDO
  M = M + 1         ! *** SHELL FISH (SFL)
  DO NS=1,NTOX
    M = M + 1
    MSVTOX(NS) = M
  ENDDO
  DO NS=1,NSED
    M = M + 1
    MSVSED(NS) = M
  ENDDO
  DO NS=1,NSND
    M = M + 1
    MSVSND(NS) = M
  ENDDO

  IF( ISTRAN(8) > 0 )THEN
    CALL SCANWQ
    
    Call AllocateDSI( MSVWQV,  NWQV,  0)
    DO MW = 1,NWQV
      IF( ISKINETICS(MW) > 0 )THEN
        M = M + 1
        MSVWQV(MW) = M        ! *** 3 + NDYM + NTOX + NSED + NSND + NW
        IF( MW == IDOX ) MSVDOX = MSVWQV(MW)
      ENDIF
    ENDDO
    !ENDIF
  ENDIF

  KSM=KCM
  IGM=ICM+1
  JGM=JCM+1
  MGM=2*MTM
  NSTM   = MAX(3, NSCM  + NSNM  + NTXM)                      ! *** Maximum number of SedTox variables
  NSTVM  = MAX(7, 3 + NDYM + NSCM  + NSNM  + NTXM + NWQV)    ! *** Maximum number of constituents

  NQINFLM = MAX(1,NQSIJ+NQCTL+NQWR+2*MDCHH)

  ! *** CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENTS
  ISWQFLUX = 0
  IF( ISTRAN(4) >= 1 ) ISWQFLUX = 1
  IF( ISTRAN(8) >= 1 ) ISWQFLUX = 1
  IF( ISWASP >= 1 )    ISWQFLUX = 1
  IF( ISICM >= 1  )    ISWQFLUX = 1
  IF( ISSSMMT > 0 )    ISWQFLUX = 1
  IF( ISLSHA == 1)     ISWQFLUX = 1

  ! *** Allocate the main EFDC+ arrays
  CALL VARALLOC                                 ! *** Each process needs to allocate variables
  
  IF( LSEDZLJ ) CALL VARZEROSNL

  ! *** MPI subdomain mapping
  Map2Global(:).PR = -999
  Map2Global(:).IL = -999
  Map2Global(:).JL = -999
  Map2Global(:).LL = -999
  Map2Global(:).IG = -999
  Map2Global(:).JG = -999
  Map2Global(:).LG = -999

  Map2Local(:).PR = -999
  Map2Local(:).IL = -999
  Map2Local(:).JL = -999
  Map2Local(:).LL = -999
  Map2Local(:).IG = -999
  Map2Local(:).JG = -999
  Map2Local(:).LG = -999

  LWVMASK=.FALSE.
  
  !*** Vegetation parameters from Katul et al., 2003
  BETASUP_P = 1.
  BETASUP_D = 5.1
  BETAVEG_P = 1.
  BETAVEG_D = 5.1
  CE4SUP = 0.9
  CE4VEG = 0.9
  
  IF( ISICE > 0 )THEN
    MITLAST    = 0
    FRAZILICE  = 0.0
    FRAZILICE1 = 0.0
    ICECOVER   = 0.0
    ICERATE    = 0.0
    ICETHICK   = 0.0
    ICETHICK1  = 0.0
    ICEVOL     = 0.0
    ICETEMP    = 0.
    TCISER     = 0.
    TAISER     = 0.
    RICECOVT   = 0.
    RICETHKT   = 0.
    IF( NISER > 1) RICEWHT = 1.0
  ENDIF
  ICECELL  = .FALSE.
  RHOW = 1000.
  SHLIM = 44.         ! *** MAXIMUM ANGULAR WAVE NUMBER*DEPTH

  ! *** FOR WATER COLUMN POINTERS
  Call AllocateDSI( IACTIVEWC1,  100,  0)
  Call AllocateDSI( IACTIVEWC2,  100,  0)

  NACTIVEWC = 0
  IF( ISTRAN(1) > 0 )THEN
    NACTIVEWC = NACTIVEWC+1
    IACTIVEWC1(NACTIVEWC) = 1
    IACTIVEWC2(NACTIVEWC) = 1
    WCV(NACTIVEWC).ID = 'SAL'
    WCV(NACTIVEWC).VAL0 => SAL
    WCV(NACTIVEWC).VAL1 => SAL1
    WCV(NACTIVEWC).WCLIMIT = 0.0
  ENDIF
  IF( ISTRAN(2) > 0 )THEN
    NACTIVEWC = NACTIVEWC+1
    IACTIVEWC1(NACTIVEWC) = 2
    IACTIVEWC2(NACTIVEWC) = 2
    WCV(NACTIVEWC).ID = 'TEM'
    WCV(NACTIVEWC).VAL0 => TEM
    WCV(NACTIVEWC).VAL1 => TEM1
    WCV(NACTIVEWC).WCLIMIT = 0.0
  ENDIF
  IF( ISTRAN(3) > 0 )THEN
    DO MD = 1,NDYE
      NACTIVEWC = NACTIVEWC+1
      IACTIVEWC1(NACTIVEWC) = 3
      IACTIVEWC2(NACTIVEWC) = MSVDYE(MD)
      WRITE(NSTR,'(I3.3)') MD
      WCV(NACTIVEWC).ID = 'DYE'//NSTR
      WCV(NACTIVEWC).VAL0 => DYE(:,:,MD)
      WCV(NACTIVEWC).VAL1 => DYE1(:,:,MD)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    ENDDO
  ENDIF
  IF( ISTRAN(4) > 0 )THEN
    NACTIVEWC = NACTIVEWC+1
    IACTIVEWC1(NACTIVEWC) = 4
    IACTIVEWC2(NACTIVEWC) = 3 + NDYM
    WCV(NACTIVEWC).ID = 'SFL'
    WCV(NACTIVEWC).VAL0 => SFL
    WCV(NACTIVEWC).VAL1 => SFL2
    WCV(NACTIVEWC).WCLIMIT = 0.0
  ENDIF
  IF( ISTRAN(5) > 0 )THEN
    DO NS = 1,NTOX
      NACTIVEWC = NACTIVEWC+1
      IACTIVEWC1(NACTIVEWC) = 5
      IACTIVEWC2(NACTIVEWC) = MSVTOX(NS)
      WRITE(NSTR,'(I3.3)') NS
      WCV(NACTIVEWC).ID = 'TOX'//NSTR
      WCV(NACTIVEWC).VAL0 => TOX(:,:,NS)
      WCV(NACTIVEWC).VAL1 => TOX1(:,:,NS)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    ENDDO
  ENDIF
  IF( ISTRAN(6) > 0 )THEN
    DO NS = 1,NSED
      NACTIVEWC = NACTIVEWC+1
      IACTIVEWC1(NACTIVEWC) = 6
      IACTIVEWC2(NACTIVEWC) = MSVSED(NS)
      WRITE(NSTR,'(I3.3)') NS
      WCV(NACTIVEWC).ID = 'SED'//NSTR
      WCV(NACTIVEWC).VAL0 => SED(:,:,NS)
      WCV(NACTIVEWC).VAL1 => SED1(:,:,NS)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    ENDDO
  ENDIF
  IF( ISTRAN(7) > 0 )THEN
    DO NS = 1,NSND
      NACTIVEWC = NACTIVEWC+1
      IACTIVEWC1(NACTIVEWC) = 7
      IACTIVEWC2(NACTIVEWC) = MSVSND(NS)
      WRITE(NSTR,'(I3.3)') NS
      WCV(NACTIVEWC).ID = 'SND'//NSTR
      WCV(NACTIVEWC).VAL0 => SND(:,:,NS)
      WCV(NACTIVEWC).VAL1 => SND1(:,:,NS)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    ENDDO
  ENDIF
  IF( ISTRAN(8) > 0 )THEN
    NCSER(8) = NCSERM
    DO MW = 1,NWQV
      IF( ISTRWQ(MW) == 1 )THEN
        NACTIVEWC = NACTIVEWC+1
        IACTIVEWC1(NACTIVEWC) = 8
        IACTIVEWC2(NACTIVEWC) = MSVWQV(MW)
        WRITE(NSTR,'(I3.3)') MW
        WCV(NACTIVEWC).ID = 'WQV'//NSTR
        WCV(NACTIVEWC).VAL0 => WQV(:,:,MW)
        WCV(NACTIVEWC).VAL1 => WQV(:,:,MW)
        WCV(NACTIVEWC).WCLIMIT = 0.0
      ENDIF
    ENDDO
  ENDIF

  ! *** MPI DECLARATIONS
  Call AllocateDSI( CON2,  LCM, KCM, NACTIVEWC, 0.0)
  Call AllocateDSI( FUHUD, LCM, KCM, NACTIVEWC, 0.0)
  Call AllocateDSI( FVHUD, LCM, KCM, NACTIVEWC, 0.0)

  Call AllocateDSI( FWUU,  LCM, -KCM, NACTIVEWC, 0.0)

  ! *** Water column fluxes for output
  Call AllocateDSI( WC_UP, 100, KCM, NACTIVEWC, 0.0 )
  Call AllocateDSI( WC_DN, 100, KCM, NACTIVEWC, 0.0 )
  Call AllocateDSI( WC_QU, 100, KCM, 0.0 )
  Call AllocateDSI( WC_QD, 100, KCM, 0.0 )
  
  ! *** MPI Timing variables
  DSITIME = '        '
  DSITIME(1)  = 'CONGRAD'
  DSITIME(2)  = 'CALPUV'
  DSITIME(3)  = 'CALUVW'
  DSITIME(4)  = 'CALTRN1'
  DSITIME(5)  = 'CALTRN2'
  DSITIME(6)  = 'CALCONC'
  DSITIME(7)  = 'CALQQ'
  DSITIME(8)  = 'BEDLOAD'
  DSITIME(9)  = 'HDMT1'
  DSITIME(10) = 'TOTGHOST'
  DSITIME(11) = 'TOT_AR'
  DSITIME(12) = 'ARCONGRD'
  DSITIME(13) = 'ARCALSTP'
  DSITIME(14) = ''
  DSITIME(15) = ''
  DSITIME(16) = ''
   
END

