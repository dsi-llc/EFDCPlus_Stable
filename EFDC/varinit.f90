! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!< @author Paul Craig
!< @date 2015-06 - Implemented SIGMA-Z (SGZ) IN EE7.3
!< @date 2022-01 - Added mass erosion class for propwash
!< @details Scans the input files and calls routines to allocate and zero out arrays
!< @details   does some other setup as well.

SUBROUTINE VARINIT

  use GLOBAL
  use Allocate_Initialize
  use SCANINPMOD
  use Variables_MPI
  use MPI
  use Broadcast_Routines
  
  implicit none

  integer :: NCSER1, NCSER2, NCSER3, NCSER4, NCSER5, NCSER6, NCSER7
  integer :: MW, NS, NT, M, MD, L, NGOTM
  character*3 :: NSTR

  integer :: ierr ! MPI error return value

  NDQCLT   = 1
  NDVEGSER = 1
  NJPSM    = 1
  NPDM     = 1
  NSMGM    = 1
  NSMZM    = 1
  NVEGSERM = 1
  NVEGTPM  = 100
  NWQPSM   = 1
  NWQPSRM  = 1
  NWQZM    = 1
  NWNDMAP  = 1
  NATMMAP  = 1
  NICEMAP  = 1
  
  call AllocateDSI( NCSER,    8, 0)
  
  ! *** Scanning DECOMP file to setup variables for allocation
  call Scan_JSON_Decomp

  ! *** After we have the partitioning from the DECOMP file setup the MPI topology
  call Setup_MPI_Topology

  ! *** SCAN MAIN CONTROL FILE - only on the master
  call SCANEFDC(NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7)
  
  ! *** Initialize optimal thread counts
  NOPTIMAL = 1
  LDMOPT   = INT(FLOAT(LCM-1)/FLOAT(NOPTIMAL)) + 1
  call AllocateDSI( TIMES,  -50, NTHREADS, 864000.)

  if( ISGOTM > 0 )then
    NGOTM = 2
  else
    NGOTM = 0
  endif
  
  ! ***  Read domain decomposition settings
  call Read_JSON_Decomp
  
  ! *** Scans the cell.inp and sets up some variables used for allocations later on
  call Scan_Cell


  ! *** Initialize
  !global_max_width_x = icm
  !global_max_width_y = jcm
  !max_width_x = global_max_width_x
  !max_width_y = global_max_width_y
  !n_x_partitions = 1
  !n_y_partitions = 1
  !active_domains = 1

  ! *** Send out the global values originally read in
  call Broadcast_Scalar(ICM_Global, master_id)
  call Broadcast_Scalar(JCM_Global, master_id)
  call Broadcast_Scalar(LCM_Global, master_id)
  call Broadcast_Scalar(LC_Global,  master_id)

  ! *** Set some variables for array allocation
  max_width_x = max_width_x + 4  ! Adding the plus 4 because of the 2 ghost rows on each side of domain
  max_width_y = max_width_y + 4

  global_max_width_x = icm + 4   ! So this is IC + 1 + 4
  global_max_width_y = jcm + 4   ! So this is IC + 1 + 4

  ! *** Remapping ICM and JCM to local values
  icm = max_width_x
  jcm = max_width_y

  call WriteBreak(mpi_log_unit)
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
  call WriteBreak(mpi_log_unit)

  ! *** Scan other files for allocation of variables
  if( ISCHAN > 0  ) CALL SCANMODC
  if( ISGWIT == 2 ) CALL SCANGWSR
  if( NASER  > 0  ) CALL SCANASER
  if( NCSER1 > 0 .and. ISTRAN(1) >= 1 ) CALL SCANSSER
  if( (NCSER2 > 0 .and. ISTRAN(2) >= 1) .or. ((ISTRAN(3) > 0 .or. ISTRAN(5) > 0) .and. IVOLTEMP > 1 ) ) CALL SCANTSER
  if( NCSER3 > 0  .and. ISTRAN(3) >= 1 ) CALL SCANDSER
  if( NCSER4 > 0  .and. ISTRAN(4) >= 1 ) CALL SCANSFSR
  if( NQSER  > 0  ) CALL SCANQSER
  if( NPSER  > 0  ) CALL SCANPSER
  if( NWSER  > 0  ) CALL SCANWSER
  if( NQCTLT >= 1 ) CALL SCANQCTL
  if( NQWRSR > 0  ) CALL SCANWRSER
  if( LSEDZLJ ) CALL SCANSEDZLJ
  if( NCSER5 > 0 .or. NCSER6 > 0 .or. NCSER7 > 0 ) CALL SCNTXSED
  if( ISPROPWASH > 0 ) CALL SCANPROPWASH

  call MPI_Barrier(DSIcomm, ierr)

  ! *** These need to be broadcast to allocate the following variables correctly
  call Broadcast_Scalar(NDYM, master_id)
  call Broadcast_Scalar(NSCM, master_id)
  call Broadcast_Scalar(NSNM, master_id)
  call Broadcast_Scalar(NTXM, master_id)

  NSND2 = NSND     ! *** Propwash does not include non-cohesives for now.
  NSNM2 = NSNM
  NSED2 = NSED
  NSEDS2 = NSEDS
  if( ISPROPWASH > 0 .and. ISTRAN(6) > 0 .and. fraction_fast > 0 .and. fast_multiplier > 0. )then
    ! *** Add WC arrays for propwash eroded WC sediments
    NSED2 = NSED*2          ! *** Maximum array sizes.  May not need all of these but will be determined later
    NSEDS2 = NSEDS*2        ! *** Maximum array sizes.  May not need all of these but will be determined later
  endif

  ! *** Time series pointers for different constituents
  call AllocateDSI( MSVDYE,  NDYM,   0)
  call AllocateDSI( MSVSED,  NSEDS2,  0)
  call AllocateDSI( MSVSND,  NSNM2,  0)
  call AllocateDSI( MSVTOX,  NTXM,   0)

  M = 2             ! *** M = 1 (SALINITY), M = 2 (TEMPERATURE)
  do MD = 1,NDYM      ! *** MUST use NDYM INSTEAD OF NDYE TO ENSURE THE BASE OF 4 FOR LEGACY MASS BALANCE ROUTINES
    M = M + 1
    MSVDYE(MD) = M
  enddo
  M = M + 1         ! *** SHELL FISH (SFL)
  do NT = 1,NTOX
    M = M + 1
    MSVTOX(NT) = M
  enddo
  do NS = 1,NSED
    M = M + 1
    MSVSED(NS) = M
  enddo
  do NS = 1,NSND
    M = M + 1
    MSVSND(NS) = M
  enddo
  
  if( ISPROPWASH > 0 .and. ISTRAN(6) > 0 .and. NSED2 > NSED )then
    do NS = NSEDS+1, NSEDS2
      M = M + 1
      MSVSED(NS) = M
    enddo
  endif
  
  if( ISTRAN(8) > 0 )then
    call SCANWQ
    
    call AllocateDSI( MSVWQV,  NWQV,  0)
    do MW = 1,NWQV
      if( ISKINETICS(MW) > 0 )then
        M = M + 1
        MSVWQV(MW) = M        ! *** 3 + NDYM + NTOX + NSED + NSND + NW
        if( MW == IDOX ) MSVDOX = MSVWQV(MW)
      endif
    enddo
    !ENDIF
  endif

  KSM = KCM
  NSTM   = max(3, NSCM  + NSNM  + NTXM)                      ! *** Maximum number of SedTox variables
  NSTM2  = max(3, NSEDS2 + NSNM2 + NTXM)                      ! *** Maximum number of SedTox variables with fast settling classes
  NSTVM  = max(7, 3 + NDYM + NSCM  + NSNM  + NTXM + NWQV)    ! *** Maximum number of constituents
  NSTVM2 = max(7, 3 + NDYM + NSEDS2 + NSNM2 + NTXM + NWQV)    ! *** Maximum number of constituents with fast settling classes

  NQINFLM = max(1,NQSIJ+NQCTL+NQWR+2*MDCHH)

  ! *** CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENTS
  ISWQFLUX = 0
  if( ISTRAN(4) >= 1 ) ISWQFLUX = 1
  if( ISTRAN(8) >= 1 ) ISWQFLUX = 1
  if( ISWASP >= 1 )    ISWQFLUX = 1
  if( ISICM >= 1  )    ISWQFLUX = 1
  if( ISSSMMT > 0 )    ISWQFLUX = 1

  ! *** Allocate the main EFDC+ arrays
  call VARALLOC                                 ! *** Each process needs to allocate variables
  
  if( LSEDZLJ ) CALL VARZEROSNL

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

  LWVMASK = .FALSE.
  
  !*** Vegetation parameters from Katul et al., 2003
  BETASUP_P = 1.
  BETASUP_D = 5.1
  BETAVEG_P = 1.
  BETAVEG_D = 5.1
  CE4SUP = 0.9
  CE4VEG = 0.9
  
  if( ISICE > 0 )then
    MITLAST    = 0
    FRAZILICE  = 0.0
    FRAZILICE1 = 0.0
    ICECOVER   = 0.0
    ICERATE    = 0.0
    ICETHICK   = 0.0
    ICETHICK1  = 0.0
    ICEVOL     = 0.0
    ICETEMP    = 0.
    RICECOVT   = 0.
    RICETHKT   = 0.
    if( NISER > 1) RICEWHT = 1.0
  endif
  ICECELL  = .FALSE.
  RHOW = 1000.
  SHLIM = 44.         ! *** MAXIMUM ANGULAR WAVE NUMBER*DEPTH
  SWRATNF = 0.1       ! *** Lake Tahoe extinction coefficient. Prevent possible division by zero for models without solar radiation.

  ! *** FOR WATER COLUMN pointerS
  call AllocateDSI( IACTIVEWC1,  100,  0)
  call AllocateDSI( IACTIVEWC2,  100,  0)

  NACTIVEWC = 0
  if( ISTRAN(1) > 0 )then
    NACTIVEWC = NACTIVEWC + 1
    IACTIVEWC1(NACTIVEWC) = 1
    IACTIVEWC2(NACTIVEWC) = 1
    WCV(NACTIVEWC).ID = 'SAL'
    WCV(NACTIVEWC).VAL0 => SAL
    WCV(NACTIVEWC).VAL1 => SAL1
    WCV(NACTIVEWC).WCLIMIT = 0.0
  endif
  if( ISTRAN(2) > 0 )then
    NACTIVEWC = NACTIVEWC + 1
    IACTIVEWC1(NACTIVEWC) = 2
    IACTIVEWC2(NACTIVEWC) = 2
    WCV(NACTIVEWC).ID = 'TEM'
    WCV(NACTIVEWC).VAL0 => TEM
    WCV(NACTIVEWC).VAL1 => TEM1
    WCV(NACTIVEWC).WCLIMIT = 0.0
  endif
  if( ISTRAN(3) > 0 )then
    do MD = 1,NDYE
      NACTIVEWC = NACTIVEWC + 1
      IACTIVEWC1(NACTIVEWC) = 3
      IACTIVEWC2(NACTIVEWC) = MSVDYE(MD)
      write(NSTR,'(I3.3)') MD
      WCV(NACTIVEWC).ID = 'DYE'//NSTR
      WCV(NACTIVEWC).VAL0 => DYE(:,:,MD)
      WCV(NACTIVEWC).VAL1 => DYE1(:,:,MD)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    enddo
  endif
  if( ISTRAN(4) > 0 )then
    NACTIVEWC = NACTIVEWC + 1
    IACTIVEWC1(NACTIVEWC) = 4
    IACTIVEWC2(NACTIVEWC) = 3 + NDYM
    WCV(NACTIVEWC).ID = 'SFL'
    WCV(NACTIVEWC).VAL0 => SFL
    WCV(NACTIVEWC).VAL1 => SFL2
    WCV(NACTIVEWC).WCLIMIT = 0.0
  endif
  if( ISTRAN(5) > 0 )then
    do NT = 1,NTOX
      NACTIVEWC = NACTIVEWC + 1
      IACTIVEWC1(NACTIVEWC) = 5
      IACTIVEWC2(NACTIVEWC) = MSVTOX(NT)
      write(NSTR,'(I3.3)') NT
      WCV(NACTIVEWC).ID = 'TOX'//NSTR
      WCV(NACTIVEWC).VAL0 => TOX(:,:,NT)
      WCV(NACTIVEWC).VAL1 => TOX1(:,:,NT)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    enddo
  endif
  if( ISTRAN(6) > 0 )then
    do NS = 1,NSED
      NACTIVEWC = NACTIVEWC + 1
      IACTIVEWC1(NACTIVEWC) = 6
      IACTIVEWC2(NACTIVEWC) = MSVSED(NS)
      write(NSTR,'(I3.3)') NS
      WCV(NACTIVEWC).ID = 'SED'//NSTR
      WCV(NACTIVEWC).VAL0 => SED(:,:,NS)
      WCV(NACTIVEWC).VAL1 => SED1(:,:,NS)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    enddo
  endif
  if( ISTRAN(7) > 0 )then
    do NS = 1,NSND
      NACTIVEWC = NACTIVEWC + 1
      IACTIVEWC1(NACTIVEWC) = 7
      IACTIVEWC2(NACTIVEWC) = MSVSND(NS)
      write(NSTR,'(I3.3)') NS
      WCV(NACTIVEWC).ID = 'SND'//NSTR
      WCV(NACTIVEWC).VAL0 => SND(:,:,NS)
      WCV(NACTIVEWC).VAL1 => SND1(:,:,NS)
      WCV(NACTIVEWC).WCLIMIT = 0.0
    enddo
  endif
  if( ISTRAN(8) > 0 )then
    NCSER(8) = NCSERM
    do MW = 1,NWQV
      if( ISTRWQ(MW) == 1 )then
        NACTIVEWC = NACTIVEWC + 1
        IACTIVEWC1(NACTIVEWC) = 8
        IACTIVEWC2(NACTIVEWC) = MSVWQV(MW)
        write(NSTR,'(I3.3)') MW
        WCV(NACTIVEWC).ID = 'WQV'//NSTR
        WCV(NACTIVEWC).VAL0 => WQV(:,:,MW)
        WCV(NACTIVEWC).VAL1 => WQV(:,:,MW)
        WCV(NACTIVEWC).WCLIMIT = 0.0
      endif
    enddo
  endif
  
  ! *** Fast settling classes are at the end of the active constituent list to allow for easy triming later
  if( ISPROPWASH > 0 .and. ISTRAN(6) > 0 .and. NSED2 > NSED )then
    do NS = NSED+1, NSED2
      NACTIVEWC = NACTIVEWC + 1
      IACTIVEWC1(NACTIVEWC) = 6
      IACTIVEWC2(NACTIVEWC) = MSVSED(NS)
      write(NSTR,'(I3.3)') NS
      WCV(NACTIVEWC).ID = 'SED'//NSTR
      WCV(NACTIVEWC).VAL0 => SED(:,:,NS)
      WCV(NACTIVEWC).VAL1 => SED1(:,:,NS)
      if( PROCESS_ID == 0 ) PRINT *,'WC LIMIT SET FOR ', WCV(NACTIVEWC).ID
      WCV(NACTIVEWC).WCLIMIT = 1.0E-6
    enddo
  endif
  
  ! *** MPI DECLARATIONS
  call AllocateDSI( CON2,  LCM,  KCM, NACTIVEWC,         0.0)
  call AllocateDSI( FUHUD, LCM,  KCM, NACTIVEWC + NGOTM, 0.0)
  call AllocateDSI( FVHUD, LCM,  KCM, NACTIVEWC + NGOTM, 0.0)
  call AllocateDSI( FWUU,  LCM, -KCM, NACTIVEWC + NGOTM, 0.0)

  ! *** Water column fluxes for output
  call AllocateDSI( WC_UP, 100, KCM, NACTIVEWC, NSEDS2+2, 0.0 )
  call AllocateDSI( WC_DN, 100, KCM, NACTIVEWC, NSEDS2+2, 0.0 )
  call AllocateDSI( WC_QU, 100, KCM, 0.0 )
  call AllocateDSI( WC_QD, 100, KCM, 0.0 )
  
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

