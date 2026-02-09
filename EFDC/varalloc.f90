! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE VARALLOC

  !----------------------------------------------------------------------!
  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3

  use GLOBAL
  use Allocate_Initialize
  use HYDSTRUCMOD
  use Variables_WQ
  
  use MPI
  use Variables_MPI
  use Broadcast_Routines
#ifndef GNU  
  USE IFPORT
#endif

  implicit none

  integer :: L, ierr, NJP, NMAX
  integer :: NS, lower_bound, upper_bound
  logical Test
  
  ! *** Only write to LOG file from the master process
  if( process_id == master_id  )then
    write(*, '(A)')'ALLOCATING ARRAYS'
  endif

  if( LSEDZLJ )then
    NMAX = NSEDS
  else
    NMAX = NSCM
  endif

  ! *** Write arrays sizes to file
  call Write_Array_Sizes
  call MPI_Barrier(DSIcomm, ierr)

  ! *** Calculates max size of some arrays
  NMAXBC = MDCHH*2 + NQSIJ + NQCTL*2 + NQWR*2 + NPBS + NPBW + NPBE + NPBN + NQJPIJ*2

  ! *** Need to broadcast NMAXBC size of arrays to each process so they can allocate arrays of the correct size
  call broadcast_scalar_int(NMAXBC, master_id)

  allocate(Map2Global(LCM))
  allocate(Map2Local(LCM_Global))

  ! *** Character arrays
  allocate(CLTMSR(MLTMSRM))
  allocate(SYMBOL(MTM))
  CLTMSR = ' '
  SYMBOL = ' '
  
  ! *** Special Range arrays
  allocate(LIJ(-1:ICM,    -1:JCM))
  allocate(LIJLT(-1:ICM,  -1:JCM))
  allocate(IJCT(-1:ICM,   -1:JCM))
  allocate(IJCTLT(-1:ICM, -1:JCM))
  LIJ = 0
  LIJLT = 0
  IJCT = 0
  IJCTLT = 0
  
  ! *** Type Structures
  if( ISTRAN(3) > 0 .and. NDYE > 0 )then
    allocate(DYES(NDYE))
  endif
  allocate(WV(LCM))

  allocate(BCPS(NQSIJM))             ! *** Flow BC
  allocate(BCPS_GL(NQSIJM))          ! *** Flow BC
  do NS = 1, NQSIJM
    allocate(BCPS(NS).NCSERQ(NSTVM2))
    allocate(BCPS_GL(NS).NCSERQ(NSTVM2))
    BCPS(NS).NCSERQ    = 0
    BCPS_GL(NS).NCSERQ = 0

    allocate(BCPS(NS).CQSE(NSTVM2))
    allocate(BCPS_GL(NS).CQSE(NSTVM2))
    BCPS(NS).CQSE    = 0.0
    BCPS_GL(NS).QSSE = 0.0
  enddo
  
  allocate(HYD_STR_GL(NQCTLM))       ! *** Hydraulic structure configuration
  allocate(HSCTL_GL(NQCTLM))         ! *** Hydraulic structure control
                                    
  allocate(WITH_RET_GL(NQWRM))       ! *** Withdrawal-Return configuration
  allocate(WITH_RET_CTL_GL(NQWRM))   ! *** Withdrawal-Return control
  
  allocate(JET_PLM_GL(NJPSM))        ! *** Jet-Plume configuration
  do NJP = 1,NJPSM
    call AllocateDSI( JET_PLM_GL(NJP).NCSERJP, -NSTVM2,    0) 
    call AllocateDSI( JET_PLM_GL(NJP).CWRCJP,  -NSTVM2,  0.0)
    call AllocateDSI( JET_PLM_GL(NJP).CQCJP,    KCM, -NSTVM2, 0.0)
  enddo
  
  ! *** Each process allocates their own variable
  call AllocateDSI( AB,         LCM,        KCM,     0.0)
  call AllocateDSI( AGWELV,     LCM,        0.0)                ! *** GROUNDWATER INTERACTION WITH CELL WATER ELEVATION (M) IF ISGWIE >= 1
  call AllocateDSI( AGWELV1,    LCM,        0.0)     
  call AllocateDSI( AGWELV2,    LCM,        0.0)     
  call AllocateDSI( AH,         LCM,        KCM,     0.0)
  call AllocateDSI( AHC,        LCM,        KCM,     0.0)
  call AllocateDSI( AHDXY,      LCM,        0.0)     
  call AllocateDSI( AHOXY,      LCM,        0.0)     
  call AllocateDSI( APCG,       LCM,        0.0)     
  call AllocateDSI( AQ,         LCM,        KCM,     0.0)
  call AllocateDSI( AQCTL,      NQCTTM,     0.0)     
  call AllocateDSI( ATMP,       LCM,        0.0)     
  call AllocateDSI( AV,         LCM,        KCM,     0.0)
  call AllocateDSI( AVOXY,      LCM,        0.0)     
  call AllocateDSI( AVBXY,      LCM,        0.0)     
  call AllocateDSI( AVUI,       LCM,        KCM,     0.0)
  call AllocateDSI( AVVI,       LCM,        KCM,     0.0)
  call AllocateDSI( B,          LCM,        KCM,     0.0)
  call AllocateDSI( B1,         LCM,        KCM,     0.0)
  call AllocateDSI( BELAGW,     LCM,        0.0)     
  call AllocateDSI( BELV,       LCM,        0.0)     
  call AllocateDSI( BELV0,      LCM,        0.0)     
  call AllocateDSI( BELV1,      LCM,        0.0)     
  
  call AllocateDSI( CAC,        LCM,        KCM,     0.0)
  call AllocateDSI( CE,         LCM,        0.0)
  call AllocateDSI( CN,         LCM,        0.0)             
  call AllocateDSI( CS,         LCM,        0.0)
  call AllocateDSI( CW,         LCM,        0.0)
  call AllocateDSI( CBE,        NBBEM,      2,       NSTVM2,0.0)
  call AllocateDSI( CBN,        NBBNM,      2,       NSTVM2,0.0)
  call AllocateDSI( CBS,        NBBSM,      2,       NSTVM2,0.0)
  call AllocateDSI( CBW,        NBBWM,      2,       NSTVM2,0.0)
  call AllocateDSI( CC,         LCM,        0.0)
  call AllocateDSI( CCC,        LCM,        0.0)
  call AllocateDSI( CCCCHH,     NCHANM,     0.0)
  call AllocateDSI( CCCCHU,     NCHANM,     0.0)
  call AllocateDSI( CCCCHV,     NCHANM,     0.0)
  call AllocateDSI( CCCI,       LCM,        0.0)
  call AllocateDSI( CCE,        LCM,        0.0)
  call AllocateDSI( CCN,        LCM,        0.0)
  call AllocateDSI( CCNHTT,     LCM,        0.0)                           ! *** WIND FUNCTION FOR CONDUCTIVE HEAT TRANSFER
  call AllocateDSI( CCS,        LCM,        0.0)
  call AllocateDSI( CCW,        LCM,        0.0)
  call AllocateDSI( CDZKK,      LCM,        KCM, 0.0)
  call AllocateDSI( CDZKKP,     LCM,        KCM, 0.0)
  call AllocateDSI( CDZKMK,     LCM,        KCM, 0.0)
  call AllocateDSI( CHANFRIC,   NCHANM,     0.0)
  call AllocateDSI( CHANLEN,    NCHANM,     0.0)
  call AllocateDSI( CLEVAP,     LCM,        0.0)                        ! *** WIND FUNCTION FOR EVAPORATION CALCULATIONS
  call AllocateDSI( CLOUDT,     LCM, 0.0  )                  
  call AllocateDSI( CONGW,      LCM,        NSTVM2,    0.0)  
  call AllocateDSI( CONT,       LCM,        KCM,       0.0)  
  call AllocateDSI( CQWR,       NQWRM,      NSTVM2,    0.0)
  call AllocateDSI( CSERT,      KCM,       -NCSERM, -NSTVM2, 0.0)
  call AllocateDSI( CTAUC,      LCM,        0.0)
  call AllocateDSI( CTURBB1,    LCM,        KCM,     0.0)
  call AllocateDSI( CTURBB2,    LCM,        KCM,     0.0)
  call AllocateDSI( CU1,        LCM,        KCM,     0.0)
  call AllocateDSI( CU2,        LCM,        KCM,     0.0)
  call AllocateDSI( CUE,        LCM,        0.0)
  call AllocateDSI( CUN,        LCM,        0.0)
  call AllocateDSI( CUU,        LCM,        0.0)
  call AllocateDSI( CVE,        LCM,        0.0)
  call AllocateDSI( CVN,        LCM,        0.0)
  call AllocateDSI( CVV,        LCM,        0.0)
  call AllocateDSI( DLAT,       LCM,        0.0)
  call AllocateDSI( DLON,       LCM,        0.0)
  call AllocateDSI( DML,        LCM,       -KCM,      0.0)
  call AllocateDSI( DU,         LCM,        KCM,      0.0)
  call AllocateDSI( DV,         LCM,        KCM,      0.0)
  call AllocateDSI( DXDJ,       LCM,        0.0)
  call AllocateDSI( DXIU,       LCM,        0.0)
  call AllocateDSI( DXIV,       LCM,        0.0)
  call AllocateDSI( DXP,        LCM,        0.0)
  call AllocateDSI( DXU,        LCM,        0.0)
  call AllocateDSI( DXV,        LCM,        0.0)
  call AllocateDSI( DXYIP,      LCM,        0.0)
  call AllocateDSI( DXYIU,      LCM,        0.0)
  call AllocateDSI( DXYIV,      LCM,        0.0)
  call AllocateDSI( DXYP,       LCM,        0.0)
  call AllocateDSI( DXYU,       LCM,        0.0)
  call AllocateDSI( DXYV,       LCM,        0.0)
  call AllocateDSI( DYDI,       LCM,        0.0)
  call AllocateDSI( DYE,        LCM,        KCM,      NDYM, 0.0)
  call AllocateDSI( DYE1,       LCM,        KCM,      NDYM, 0.0)
  call AllocateDSI( DYEINIT,    LCM,        KCM,      NDYM, 0.0)
  call AllocateDSI( DYIU,       LCM,        0.0)     
  call AllocateDSI( DYIV,       LCM,        0.0)     
  call AllocateDSI( DYP,        LCM,        0.0)     
  call AllocateDSI( DYU,        LCM,        0.0)     
  call AllocateDSI( DYV,        LCM,        0.0)     
  call AllocateDSI( DZIC,       LCM,       -KCM,      0.0)
  call AllocateDSI( DZIG,       LCM,       -KCM,      0.0)
  call AllocateDSI( EVAPGW,     LCM,        0.0)                     ! *** TRANSPIRATION RATE OF WATER PER CELL (M3/S)
  call AllocateDSI( EVAPSW,     LCM,        0.0)                     ! *** EVAPORATION RATE OF WATER PER CELL (M3/S)
  call AllocateDSI( EVAPT,      LCM,        0.0)                     ! *** EVAPORATION RATE OF WATER PER CELL (M/S)
  call AllocateDSI( EVACOARE,   LCM,        0.0)                     ! *** EVAPORATION RATE OF WATER PER CELL (M/S)
  call AllocateDSI( FACBEDL,    LCM,        0.0)      
  call AllocateDSI( FACSUSL,    LCM,        0.0)      
  call AllocateDSI( FBBX,       LCM,        KCM,      0.0)
  call AllocateDSI( FBBY,       LCM,        KCM,      0.0)
  if( ISBODYF > 0 )then                             
    call AllocateDSI( FBODYFX,   LCM,       KCM,      0.0)
    call AllocateDSI( FBODYFY,   LCM,       KCM,      0.0)
  endif                                             
  call AllocateDSI( FCAX,       LCM,        KCM,      0.0)
  call AllocateDSI( FCAXE,      LCM,        0.0)      
  call AllocateDSI( FCAY,       LCM,        KCM,      0.0)
  call AllocateDSI( FCAY1,      LCM,        KCM,      0.0)
  call AllocateDSI( FCAY1E,     LCM,        0.0)      
  call AllocateDSI( FCAYE,      LCM,        0.0)      
  call AllocateDSI( FCORC,      LCM,        0.0)      
  call AllocateDSI( FP,         LCM,        0.0)      
  call AllocateDSI( FP1,        LCM,        0.0)      
  call AllocateDSI( FPGXE,      LCM,        0.0)      
  call AllocateDSI( FPGYE,      LCM,        0.0)      
  call AllocateDSI( FPROX,      LCM,       -KCM,      0.0)
  call AllocateDSI( FPTMP,      LCM,        0.0)      
  call AllocateDSI( FSCORTBCV,  LCM,        0.0)      
  call AllocateDSI( FUHDYE,     LCM,        0.0)      
  call AllocateDSI( FUHU,       LCM,        KCM,      0.0)
  call AllocateDSI( FUHV,       LCM,        KCM,      0.0)
  call AllocateDSI( FVHDXE,     LCM,        0.0)      
  call AllocateDSI( FVHU,       LCM,        KCM,      0.0)
  call AllocateDSI( FVHV,       LCM,        KCM,      0.0)
  call AllocateDSI( FWQQ,       LCM,        KCM,      0.0)
  call AllocateDSI( FWQQL,      LCM,        KCM,      0.0)
  call AllocateDSI( FWU,        LCM,       -KCM,      0.0)
  call AllocateDSI( FWV,        LCM,       -KCM,      0.0)
  call AllocateDSI( FX,         LCM,        KCM,      0.0)
  call AllocateDSI( FX1,        LCM,        KCM,      0.0)
  call AllocateDSI( FXE,        LCM,        0.0)      
  call AllocateDSI( FY,         LCM,        KCM,      0.0)
  call AllocateDSI( FY1,        LCM,        KCM,      0.0)
  call AllocateDSI( FYE,        LCM,        0.0)      
  call AllocateDSI( GWCSER,     NDGWSER,    NGWSERM,  NSTVM2, 0.0)
  call AllocateDSI( GWCSERT,   -NGWSERM,   -NSTVM2,   0.0)
  call AllocateDSI( GWFAC,      LCM,        0.0)      
  call AllocateDSI( GWSER,      NDGWSER,    NGWSERM,  0.0)
  call AllocateDSI( GWSERT,    -NGWSERM,    0.0)
  call AllocateDSI( H1P,        LCM,        0.0)
  call AllocateDSI( H1U,        LCM,        0.0)
  call AllocateDSI( H1UI,       LCM,        0.0)
  call AllocateDSI( H1V,        LCM,        0.0)
  call AllocateDSI( H1VI,       LCM,        0.0)
  call AllocateDSI( H2P,        LCM,        0.0)
  call AllocateDSI( H2WQ,       LCM,        0.0)
  call AllocateDSI( HBEDA,      LCM,        0.0)
  call AllocateDSI( HBEDA1,     LCM,        0.0)
  call AllocateDSI( HCTLDA,     NQCTTM,     0.0)
  call AllocateDSI( HCTLDM,     NQCTTM,     0.0)
  call AllocateDSI( HCTLUA,     NQCTTM,     0.0)
  call AllocateDSI( HCTLUM,     NQCTTM,     0.0)
  call AllocateDSI( HDFUFX,     LCM,        0.0)
  call AllocateDSI( HDFUFY,     LCM,        0.0)
  call AllocateDSI( HDFUF,      LCM,        0.0)
  call AllocateDSI( HDIFCTD,    NDQCLT,     NQCTTM,   0.0)
  call AllocateDSI( HDIFCTL,    NDQCLT,     NQCTTM,   0.0)
  call AllocateDSI( HGDH,       LCM,        0.0)
  call AllocateDSI( HMP,        LCM,        0.0)
  call AllocateDSI( HMPW,       LCM,        0.0)
  call AllocateDSI( HMU,        LCM,        0.0)
  call AllocateDSI( HMUW,       LCM,        0.0)
  call AllocateDSI( HMV,        LCM,        0.0)
  call AllocateDSI( HMVW,       LCM,        0.0)
  call AllocateDSI( HP,         LCM,        0.0)
  call AllocateDSI( HPI,        LCM,        0.0)
  call AllocateDSI( HPTMP,      LCM,        0.0)
  call AllocateDSI( HRU,        LCM,        0.0)
  call AllocateDSI( HRUO,       LCM,        0.0)            ! *** U INTERFACE METRIC (DY AT THE U FACE) / (DX AT THE U FACE) [DYU(L)*DXIU(L)] (DIMENSIONLESS)
  call AllocateDSI( HRV,        LCM,        0.0)
  call AllocateDSI( HRVO,       LCM,        0.0)            ! *** V INTERFACE METRIC (DX AT THE V FACE) / (DY AT THE V FACE) [DXV(L)*DYIV(L)] (DIMENSIONLESS)
  call AllocateDSI( HTMP,       LCM,        0.0)
  call AllocateDSI( HU,         LCM,        0.0)
  call AllocateDSI( HUI,        LCM,        0.0)
  call AllocateDSI( HUTMP,      LCM,        0.0)
  call AllocateDSI( HV,         LCM,        0.0)
  call AllocateDSI( HVI,        LCM,        0.0)
  call AllocateDSI( HVTMP,      LCM,        0.0)
  call AllocateDSI( HWQ,        LCM,        0.0)
  call AllocateDSI( HWQI,       LCM,        0.0)
                                         
  call AllocateDSI( ICBE,       NBBEM,        0)      
  call AllocateDSI( ICBN,       NBBNM,        0)      
  call AllocateDSI( ICBS,       NBBSM,        0)      
  call AllocateDSI( ICBW,       NBBWM,        0)      
  call AllocateDSI( ICFLMP,     LCM,          0)      
  call AllocateDSI( IGWSER,    -NGWSERM,      0)
  call AllocateDSI( IL,         LCM,          0)
  call AllocateDSI( ILLT,       LCM,          0)
  call AllocateDSI( ILTMSR,     MLTMSRM,      0)
  !ALLOCATE(IMASKDRY(LCM, 0. 0 )                        ! *** DEPRECATED - NO LONGER USED 2015-02
  call AllocateDSI( IMDCHH,     NCHANM,       0)
  call AllocateDSI( IMDCHU,     NCHANM,       0)
  call AllocateDSI( IMDCHV,     NCHANM,       0)
  call AllocateDSI( INTPSER,    NPSERM,       0)
  
  call AllocateDSI( IPBE,       NPBEM,        0)
  call AllocateDSI( IPBN,       NPBNM,        0)
  call AllocateDSI( IPBS,       NPBSM,        0)
  call AllocateDSI( IPBW,       NPBWM,        0)
  
  call AllocateDSI( ISCDRY,     LCM,          0)
  
  call AllocateDSI( ISPBE,      NPBEM,        0)
  call AllocateDSI( ISPBN,      NPBNM,        0)
  call AllocateDSI( ISPBS,      NPBSM,        0)
  call AllocateDSI( ISPBW,      NPBWM,        0)
  
  call AllocateDSI( ISPRE,      NPBEM,        0)
  call AllocateDSI( ISPRN,      NPBNM,        0)
  call AllocateDSI( ISPRS,      NPBSM,        0)
  call AllocateDSI( ISPRW,      NPBWM,        0)
  
  call AllocateDSI( ISSBCP,     LCM,          0)
  call AllocateDSI( ISUDPC,     NCHANM,       0)
  
  call AllocateDSI( JCBE,       NBBEM,        0)
  call AllocateDSI( JCBN,       NBBNM,        0)
  call AllocateDSI( JCBS,       NBBSM,        0)
  call AllocateDSI( JCBW,       NBBWM,        0)
  
  call AllocateDSI( JL,         LCM,          0)
  call AllocateDSI( JLLT,       LCM,          0)
  call AllocateDSI( JLTMSR,     MLTMSRM,      0)
  call AllocateDSI( JMDCHH,     NCHANM,       0)
  call AllocateDSI( JMDCHU,     NCHANM,       0)
  call AllocateDSI( JMDCHV,     NCHANM,       0)
  
  call AllocateDSI( JPBE,       NPBEM,        0)
  call AllocateDSI( JPBN,       NPBNM,        0)
  call AllocateDSI( JPBS,       NPBSM,        0)
  call AllocateDSI( JPBW,       NPBWM,        0)
  
  call AllocateDSI( KBT,        LCM,          0)
  call AllocateDSI( KEFFJP,     NJPSM,        0)
  call AllocateDSI( KUPW,       LCM,        KCM,      0)
  
  call AllocateDSI( LADJ,       9,          LCM,      0)
  call AllocateDSI( LBERC,      NMAXBC,       0)
  call AllocateDSI( LBNRC,      NMAXBC,       0)
  call AllocateDSI( LBSRC,      NMAXBC,       0)
  call AllocateDSI( LBWRC,      NMAXBC,       0)
  
  call AllocateDSI( LCBE,       NBBEM,        0)
  call AllocateDSI( LCBN,       NBBNM,        0)
  call AllocateDSI( LCBS,       NBBSM,        0)
  call AllocateDSI( LCBW,       NBBWM,        0)
  
  call AllocateDSI( LSC,        LCM,          1)
  call AllocateDSI( LSEC,       LCM,          1)
  call AllocateDSI( LSWC,       LCM,          1)
  call AllocateDSI( LEC,        LCM,          1)
  call AllocateDSI( LNC,        LCM,          1)
  call AllocateDSI( LNEC,       LCM,          1)
  call AllocateDSI( LNWC,       LCM,          1)
  call AllocateDSI( LWC,        LCM,          1)
  
  call AllocateDSI( LCT,        LCM,          0)
  call AllocateDSI( LMASKDRY,   LCM,    .false.)
  call AllocateDSI( LWVMASK,    LCM,    .false.)
  call AllocateDSI( LMDCHH,     NCHANM,       0)
  call AllocateDSI( LMDCHU,     NCHANM,       0)
  call AllocateDSI( LMDCHV,     NCHANM,       0)
  call AllocateDSI( LOPENBCDRY, LCM,    .false.)
  call AllocateDSI( LOBCS,      LCM,          0)
  call AllocateDSI( LOBCS2,     LCM,          0)
  call AllocateDSI( LBCS,       NMAXBC,       0)
  
  call AllocateDSI( LPBE,       NPBEM,        0)
  call AllocateDSI( LPBN,       NPBNM,        0)
  call AllocateDSI( LPBS,       NPBSM,        0)
  call AllocateDSI( LPBW,       NPBWM,        0)
  
  call AllocateDSI( LQSPATH,    NQSIJM,     100,      0)
  call AllocateDSI( LQSSAVED,    NQSIJM,       0)
  call AllocateDSI( LQSMOVED,    NQSIJM,       0)
  call AllocateDSI( LSBLBCD,    LCM,          0)
  call AllocateDSI( LSBLBCU,    LCM,          0)
  call AllocateDSI( MFDCHZ,     LCM,          0)
  call AllocateDSI( LUPU,       LCM,       KCM,       0)
  call AllocateDSI( LUPV,       LCM,       KCM,       0)
  
  call AllocateDSI( LCTLT,      LCM,          0)
  call AllocateDSI( LECLT,      LCM,          0)
  call AllocateDSI( LWCLT,      LCM,          0)
  call AllocateDSI( LNCLT,      LCM,          0)
  call AllocateDSI( LSCLT,      LCM,          0)
  
  call AllocateDSI( LWVCELL,    LCM,          0)
  call AllocateDSI( MCSER,      NCSERM,    NSTVM2,    0)
  call AllocateDSI( MTSCLAST,   NCSERM,    NSTVM2,    0)
  call AllocateDSI( MDCHTYP,    NCHANM,               0)
  call AllocateDSI( MGWSER,     NGWSERM,              0)
  call AllocateDSI( MTSGWLAST,  NGWSERM,              0)
  call AllocateDSI( MTSPLAST,   NPSERM,               0)
  call AllocateDSI( MQCTL,      NQCTTM,               0)
  call AllocateDSI( MTSQLAST,   NQSERM,               0)
  call AllocateDSI( MTMSRA,     MLTMSRM,              0)
  call AllocateDSI( MTMSRC,     MLTMSRM,              0)
  call AllocateDSI( MTMSRP,     MLTMSRM,              0)
  call AllocateDSI( MTMSRQ,     MLTMSRM,              0)
  call AllocateDSI( MTMSRQE,    MLTMSRM,              0)
  call AllocateDSI( MTMSRU,     MLTMSRM,              0)
  call AllocateDSI( MTMSRUE,    MLTMSRM,              0)
  call AllocateDSI( MTMSRUT,    MLTMSRM,              0)
  call AllocateDSI( MTSCUR,     NTSSTSPM,             0)
  call AllocateDSI( MTSSTSP,    NTSSTSPM,             0)
  call AllocateDSI( MVEGL,      LCM,                  0)
  call AllocateDSI( MVEGSER,    NVEGSERM,             0)
  call AllocateDSI( MVEGTLAST,  NVEGSERM,             0)
  call AllocateDSI( NATDRY,     LCM,                  0)
  call AllocateDSI( NCSERE,     NBBEM,     NSTVM2,   0)
  call AllocateDSI( NCSERN,     NBBNM,     NSTVM2,   0)
  call AllocateDSI( NCSERS,     NBBSM,     NSTVM2,   0)
  call AllocateDSI( NCSERW,     NBBWM,     NSTVM2,   0)
  call AllocateDSI( NSERWQ,     LCM,                  0)
  call AllocateDSI( NGWSL,      LCM,                  0)
  call AllocateDSI( NLOE,       NBBEM,     KCM,       NSTVM2, 0)
  call AllocateDSI( NLON,       NBBNM,     KCM,       NSTVM2, 0)
  call AllocateDSI( NLOS,       NBBSM,     KCM,       NSTVM2, 0)
  call AllocateDSI( NLOW,       NBBWM,     KCM,       NSTVM2, 0)
                                             
  call AllocateDSI( NPSERE,     NPBEM,        0)
  call AllocateDSI( NPSERN,     NPBNM,        0)
  call AllocateDSI( NPSERS,     NPBSM,        0)
  call AllocateDSI( NPSERW,     NPBWM,        0)
                                             
  call AllocateDSI( NPSERE1,    NPBEM,        0)
  call AllocateDSI( NPSERN1,    NPBNM,        0)
  call AllocateDSI( NPSERS1,    NPBSM,        0)
  call AllocateDSI( NPSERW1,    NPBWM,        0)
                                             
  call AllocateDSI( NTSCRE,     NBBEM,        0)
  call AllocateDSI( NTSCRN,     NBBNM,        0)
  call AllocateDSI( NTSCRS,     NBBSM,        0)
  call AllocateDSI( NTSCRW,     NBBWM,        0)
  call AllocateDSI( NTSSSS,     MLTMSRM,      0)
  call AllocateDSI( NWET,       LCM,          0)
  call AllocateDSI( OLDMASK,    LCM,    .false.)
  call AllocateDSI( P,          LCM,        0.0)
  call AllocateDSI( P1,         LCM,        0.0)
  call AllocateDSI( PATMT,      LCM,        0.0)
  call AllocateDSI( PCBE,       NPBEM,     MTM,       0.0)
  call AllocateDSI( PCBN,       NPBNM,     MTM,       0.0)
  call AllocateDSI( PCBS,       NPBSM,     MTM,       0.0)
  call AllocateDSI( PCBW,       NPBWM,     MTM,       0.0)
  call AllocateDSI( PSERZDF,   -NPSERM,    0.0)       
  call AllocateDSI( PSERZDS,   -NPSERM,    0.0)       
  call AllocateDSI( PFAM,       NPFORM,    MTM,       0.0)
  call AllocateDSI( PFPH,       NPFORM,    MTM,       0.0)
  call AllocateDSI( PMDCH,      NCHANM,    0.0)       
  call AllocateDSI( PNHYDS,     LCM,       KCM,       0.0)
  call AllocateDSI( PPH,        LCM,       0.0)       
  call AllocateDSI( PSBE,       NPBEM,     MTM,       0.0)
  call AllocateDSI( PSBN,       NPBNM,     MTM,       0.0)
  call AllocateDSI( PSBS,       NPBSM,     MTM,       0.0)
  call AllocateDSI( PSBW,       NPBWM,     MTM,       0.0)
  call AllocateDSI( PSERAVG,    NPSERM,    2,         0.0)
  call AllocateDSI( PSERT,     -NPSERM,    0.0)
  call AllocateDSI( PSERST,    -NPSERM,    0.0)
  call AllocateDSI( PSHADE,     LCM,       0.0)
  call AllocateDSI( QCELLCTR,   LCM,       0.0)
  call AllocateDSI( QCHANU,     NCHANM,    0.0)
  call AllocateDSI( QCHANUN,    NCHANM,    0.0)
  call AllocateDSI( QCHANV,     NCHANM,    0.0)
  call AllocateDSI( QCHANVN,    NCHANM,    0.0)
  call AllocateDSI( QCHNULP,    NCHANM,    0.0)
  call AllocateDSI( QCHNVLP,    NCHANM,    0.0)
  call AllocateDSI( QCTL,       NDQCLT,    NDQCLT,    KCM,  NQCTTM,  0.0)
  call AllocateDSI( QCTLST,     KCM,       NQCTLM,    0.0)
  call AllocateDSI( QCTLSTO,    KCM,       NQCTLM,    0.0)
  call AllocateDSI( QCTLT,      KCM,      -NQCTLM,    2,     0.0)
  call AllocateDSI( QCTLTLP,    KCM,       NQCTLM,    0.0)
  call AllocateDSI( QCTLTO,     KCM,      -NQCTLM,    0.0)
  call AllocateDSI( QDNEG,      LCM,       0.0)
  call AllocateDSI( QDWASTE,    LCM,       0.0)
  call AllocateDSI( QGW,        LCM,       0.0)                              ! *** GROUNDWATER FLUX TERM (M^3/S), 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT)
  call AllocateDSI( QJPENT,     KCM,       NJPSM,     0.0)
  call AllocateDSI( QJPENTT,    NJPSM,     0.0)       
  call AllocateDSI( QMORPH,     LCM,       0.0)       
  call AllocateDSI( QQ,         LCM,      -KCM,       0.0)
  call AllocateDSI( QQ1,        LCM,      -KCM,       0.0)
  call AllocateDSI( QQ2,        LCM,      -KCM,       0.0)
  call AllocateDSI( QQSQR,      LCM,      -KCM,       0.0)
  call AllocateDSI( QQL,        LCM,      -KCM,       0.0)
  call AllocateDSI( QQL1,       LCM,      -KCM,       0.0)
  call AllocateDSI( QQL2,       LCM,      -KCM,       0.0)
  if( ISGOTM > 0 )then
    call AllocateDSI( TKE3D,      LCM,      -KCM,       0.0)
    call AllocateDSI( TKE3D1,     LCM,      -KCM,       0.0)
    call AllocateDSI( EPS3D,      LCM,      -KCM,       0.0)
    call AllocateDSI( EPS3D1,     LCM,      -KCM,       0.0)
    call AllocateDSI( GL3D,       LCM,      -KCM,       0.0)
  endif
  call AllocateDSI( QQWC,       LCM,       0.0)
  call AllocateDSI( QQWCR,      LCM,       0.0)
  call AllocateDSI( QQWV1,      LCM,       0.0)
  call AllocateDSI( QQWV2,      LCM,       0.0)       
  call AllocateDSI( QQWV3,      LCM,       0.0)       
  call AllocateDSI( QRAIN,      LCM,       0.0)       
  call AllocateDSI( QSERT,      KCM,      -NQSERM,    0.0)
  call AllocateDSI( QSERCELL,   KCM,      -NQSIJM,    0.0)
  call AllocateDSI( QSRTLPN,    KCM,       NQSERM,    0.0)
  call AllocateDSI( QSRTLPP,    KCM,       NQSERM,    0.0)
  call AllocateDSI( QSSDPA,     LCM,       0.0)       
  call AllocateDSI( QSUM,       LCM,       KCM,       0.0)
  call AllocateDSI( QSUME,      LCM,       0.0)
  call AllocateDSI( QSUM1E,     LCM,       0.0)
  call AllocateDSI( CQS,        KCM,       NQSIJM,   NSTVM2, 0.0)
  call AllocateDSI( QSS,        KCM,       NQSIJM,     0.0)
  call AllocateDSI( QWATPA,     LCM,       0.0)
  call AllocateDSI( QWBDTOP,    LCM,       0.0)
  call AllocateDSI( RADBOT,     LCM,       KCM,        0.0)              ! *** SOLAR RADIATION AT THE BOTTOM OF THE LAYER  (W/M2)
  call AllocateDSI( RADNET,     LCM,       KCM,        0.0)              ! *** NET SOLAR RADIATION AT THE MIDPOINT OF THE LAYER (W/M2)
  call AllocateDSI( RADTOP,     LCM,      -KCM,        0.0)              ! *** SOLAR RADIATION AT THE TOP OF THE LAYER (W/M2)
  call AllocateDSI( RADKE,      LCM,       KCM,        0.0)              ! *** EXTINCTION COEFFICIENT AT THE MIDPOINT OF THE LAYER
  call AllocateDSI( RAINT,      LCM,       0.0)                          ! *** CURRENT RAINFALL PER CELL (M/S)
  call AllocateDSI( RBPSBL,     LCM,       0.0)        
  call AllocateDSI( RCX,        LCM,       0.0)        
  call AllocateDSI( RCY,        LCM,       0.0)        
  call AllocateDSI( RHAT,       LCM,       0.0)        
  call AllocateDSI( ROUSE,      LCM,       0.0)        
  call AllocateDSI( RSSBCE,     LCM,       0.0)        
  call AllocateDSI( RSSBCN,     LCM,       0.0)        
  call AllocateDSI( RSSBCS,     LCM,       0.0)        
  call AllocateDSI( RSSBCW,     LCM,       0.0)        
  call AllocateDSI( SAAX,       LCM,       0.0)        
  call AllocateDSI( SAAY,       LCM,       0.0)        
  call AllocateDSI( SAL,        LCM,       KCM,        0.0)
  call AllocateDSI( SAL1,       LCM,       KCM,        0.0)
  call AllocateDSI( SALINIT,    LCM,       KCM,        0.0)
  call AllocateDSI( SBX,        LCM,       0.0)        
  call AllocateDSI( SBXO,       LCM,       0.0)        
  call AllocateDSI( SBY,        LCM,       0.0)        
  call AllocateDSI( SBYO,       LCM,       0.0)        
  call AllocateDSI( SCAX,       LCM,       0.0)        
  call AllocateDSI( SCAY,       LCM,       0.0)        
  call AllocateDSI( SCB,        LCM,       0.0)        
  call AllocateDSI( SDX,        LCM,       0.0)        
  call AllocateDSI( SDY,        LCM,       0.0)        
                                                      
  call AllocateDSI( SFL,        LCM,       KCM,        0.0)
  call AllocateDSI( SFL2,       LCM,       KCM,        0.0)
  call AllocateDSI( SFLINIT,    LCM,       KCM,        0.0)
  call AllocateDSI( SFLSBOT,    LCM,       0.0)
                                         
  call AllocateDSI( SIGPHIA,    LCM,       0.0)
  call AllocateDSI( SNAPSHOTS,  5000,      0.0)  ! *** PMC - Hardwired for now
                                         
  call AllocateDSI( SOLSWRT,   LCM,        0.0)
  call AllocateDSI( SPB,       LCM,        0.0)

  call AllocateDSI( SSSS,      MTM,        0.0)
  call AllocateDSI( STBX,      LCM,        0.0)
  call AllocateDSI( STBXO,     LCM,        0.0)
  call AllocateDSI( STBY,      LCM,        0.0)
  call AllocateDSI( STBYO,     LCM,        0.0)
  call AllocateDSI( STCAP,     LCM,        0.0)
  call AllocateDSI( STCUV,     LCM,        0.0)
  call AllocateDSI( SUB,       LCM,        0.0)
  call AllocateDSI( SUBO,      LCM,        0.0)
  call AllocateDSI( SUBD,      LCM,        1.0)
  call AllocateDSI( UMASK,    -LCM,          0)
  call AllocateDSI( VMASK,    -LCM,          0)
  call AllocateDSI( SVB,       LCM,        0.0)
  call AllocateDSI( SVBO,      LCM,        0.0)
  call AllocateDSI( SVBD,      LCM,        1.0)
  call AllocateDSI( SVPAT,     LCM,        0.0)
  call AllocateDSI( SVPW,      LCM,        0.0)
  call AllocateDSI( SWB,       LCM,        0.0)
  call AllocateDSI( TAGWSER,   NGWSERM,    0.0)
  call AllocateDSI( TATMT,     LCM,        0.0)
  call AllocateDSI( TAUB,      LCM,        0.0)                   ! *** Density normalized bed shear stress (m2/s2)
  call AllocateDSI( TAVEGSER,  NVEGSERM,   0.0)
  call AllocateDSI( TBX,       LCM,        0.0)
  call AllocateDSI( TBX1,      LCM,        0.0)
  call AllocateDSI( TBY,       LCM,        0.0)
  call AllocateDSI( TBY1,      LCM,        0.0)
  call AllocateDSI( TCCSER,    NCSERM,     NSTVM,     0.0)
  call AllocateDSI( TCGWSER,   NGWSERM,    0.0)
  call AllocateDSI( TCP,       MTM,        0.0)
  call AllocateDSI( TCVEGSER,  NVEGSERM,   0.0)
  call AllocateDSI( TDEWT,     LCM,        0.0)
  call AllocateDSI( TEM,       LCM,        KCM,       0.0)
  call AllocateDSI( TEM1,      LCM,        KCM,       0.0)
  call AllocateDSI( TEMB,      LCM,        0.0)
  call AllocateDSI( TEMB1,     LCM,        0.0)
  call AllocateDSI( TEMINIT,   LCM,        KCM,       0.0)
  call AllocateDSI( TGWSER,    NDGWSER,    NGWSERM,   0.0)
                             
  call AllocateDSI( TSSTOP,    MTSSTSPM,   NTSSTSPM,  0.0)
  call AllocateDSI( TSSTRT,    MTSSTSPM,   NTSSTSPM,  0.0)
  call AllocateDSI( TSX,       LCM,        0.0)
  call AllocateDSI( TSX1,      LCM,        0.0)
  call AllocateDSI( TSY,       LCM,        0.0)
  call AllocateDSI( TSY1,      LCM,        0.0)
  call AllocateDSI( CDCOARE,   LCM,        0.0)
  call AllocateDSI( HS_OUT,    LCM,        0.0)
  call AllocateDSI( HL_OUT,    LCM,        0.0)
  call AllocateDSI( HW_OUT,    LCM,        0.0)
  call AllocateDSI( TVAR1E,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR1N,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR1S,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR1W,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR2C,    LCM,       -KCM,       0.0)
  call AllocateDSI( TVAR2E,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR2N,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR2S,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR2W,    LCM,        KCM,       0.0)
  call AllocateDSI( TVAR3C,    LCM,        0.0)
  call AllocateDSI( TVAR3E,    LCM,        0.0)
  call AllocateDSI( TVAR3N,    LCM,        0.0)
  call AllocateDSI( TVAR3S,    LCM,        0.0)
  call AllocateDSI( TVAR3W,    LCM,        0.0)
  call AllocateDSI( TVEGSER,   NDVEGSER,   NVEGSERM,  0.0)
  call AllocateDSI( TWATER,    LCM,        0.0)                      ! *** For outputing water temperature to tecplot.  Otherwise can use for temporary debugging
  call AllocateDSI( U,         LCM,        KCM,       0.0)
  call AllocateDSI( U1,        LCM,        KCM,       0.0)
  call AllocateDSI( U1V,       LCM,        0.0)      
  call AllocateDSI( U2,        LCM,        KCM,       0.0)
  call AllocateDSI( UCTR,      LCM,        KCM,       0.0)
  call AllocateDSI( UCTR1,     LCM,        KCM,       0.0)
  call AllocateDSI( UCTR2,     LCM,        KCM,       0.0)
  call AllocateDSI( UCELLCTR,  LCM,        0.0)      
  call AllocateDSI( UHDY1E,    LCM,        0.0)      
  call AllocateDSI( UHDY2E,    LCM,        0.0)      
  call AllocateDSI( UHDYE,     LCM,        0.0)      
  call AllocateDSI( UHDYEK,    LCM,        KCM,       0.0)
  call AllocateDSI( UHDY1EK,   LCM,        KCM,       0.0)
  call AllocateDSI( UHDY2EK,   LCM,        KCM,       0.0)
  call AllocateDSI( UHE,       LCM,        0.0)      
  call AllocateDSI( USTAR,     LCM,        0.0)      
  call AllocateDSI( USTARSED,  LCM,        0.0)      
  call AllocateDSI( USTARSND,  LCM,        0.0)      
  call AllocateDSI( UUU,       LCM,        KCM,       0.0)
  call AllocateDSI( UV,        LCM,        0.0)      
  call AllocateDSI( UWVSQ,     LCM,        0.0)      
  call AllocateDSI( V,         LCM,        KCM,       0.0)
  call AllocateDSI( V1,        LCM,        KCM,       0.0)
  call AllocateDSI( V1U,       LCM,        0.0)      
  call AllocateDSI( V2,        LCM,         KCM,       0.0)
  call AllocateDSI( VCTR,      LCM,         KCM,       0.0)
  call AllocateDSI( VCTR1,     LCM,         KCM,       0.0)
  call AllocateDSI( VCTR2,     LCM,         KCM,       0.0)
  call AllocateDSI( VCELLCTR, LCM,         0.0)      
  call AllocateDSI( VDWASTE,  LCM,         0.0)      
  call AllocateDSI( VEGK,     LCM,         0.0)      
  call AllocateDSI( VEGSERB,  NDVEGSER,    NVEGSERM,  0.0)
  call AllocateDSI( VEGSERBT, NVEGSERM,    0.0)
  call AllocateDSI( VEGSERH,  NDVEGSER,    NVEGSERM,  0.0)
  call AllocateDSI( VEGSERHT, NVEGSERM,    0.0)       
  call AllocateDSI( VEGSERR,  NDVEGSER,    NVEGSERM,  0.0)
  call AllocateDSI( VEGSERRT, NVEGSERM,    0.0)
  call AllocateDSI( VHDX1E,   LCM,         0.0)
  call AllocateDSI( VHDX2E,   LCM,         0.0)
  call AllocateDSI( VHDXE,    LCM,         0.0)
  call AllocateDSI( VHDXEK,   LCM,         KCM,       0.0)
  call AllocateDSI( VHDX1EK,  LCM,         KCM,       0.0)
  call AllocateDSI( VHDX2EK,  LCM,         KCM,       0.0)
  call AllocateDSI( VHE,      LCM,         0.0)

  call AllocateDSI( VPAT,     LCM,         0.0)
  call AllocateDSI( VU,       LCM,         0.0)
  call AllocateDSI( VVV,      LCM,         KCM,       0.0)
  call AllocateDSI( W,        LCM,        -KCM,       0.0)
  call AllocateDSI( W1,       LCM,        -KCM,       0.0)
  call AllocateDSI( W2,       LCM,        -KCM,       0.0)
  call AllocateDSI( WCOREST,  LCM,         0.0)
  call AllocateDSI( WCORWST,  LCM,         0.0)
  call AllocateDSI( WCORNTH,  LCM,         0.0)
  call AllocateDSI( WCORSTH,  LCM,         0.0)
  call AllocateDSI( WINDCD10, LCM,         0.0)              ! *** WIND DRAG COEFFICIENT FOR EACH CELL (NO DIM)
  call AllocateDSI( WINDST,   LCM,         0.0)              ! *** WIND SPEED AT 2 METERS FOR EACH CELL (M/S)
  call AllocateDSI( WINDSTKA, LCM,         0.0)              ! *** WIND SHELTERING COEFFICIENT
  call AllocateDSI( WINDSTKA_SAVE,         LCM,       0.0)
  call AllocateDSI( WKQ,      KCM,         0.0)       
  call AllocateDSI( WNDVELE,  LCM,         0.0)       
  call AllocateDSI( WNDVELN,  LCM,         0.0)       
  call AllocateDSI( WTCI,     KCM,         2,         0.0)
  call AllocateDSI( WVDTKEM,  LCM,         KCM,       0.0)       
  call AllocateDSI( WVDTKEP,  LCM,         KCM,       0.0)       
  call AllocateDSI( WVENEP,   LCM,         0.0)       
  call AllocateDSI( WVHUU,    LCM,         KCM,       0.0)
  call AllocateDSI( WVHUV,    LCM,         KCM,       0.0)
  call AllocateDSI( WVHVV,    LCM,         KCM,       0.0)
  call AllocateDSI( WVKHC,    LCM,         0.0)       
  call AllocateDSI( WVKHU,    LCM,         0.0)       
  call AllocateDSI( WVKHV,    LCM,         0.0)       
  call AllocateDSI( WVPP,     LCM,         KCM,       0.0)
  call AllocateDSI( WVPT,     LCM,        -KCM,       0.0)
  call AllocateDSI( WVPU,     LCM,         KCM,       0.0)
  call AllocateDSI( WVPV,     LCM,         KCM,       0.0)
  call AllocateDSI( WVTMP1,   LCM,         0.0)
  call AllocateDSI( WVTMP2,   LCM,         0.0)
  call AllocateDSI( WVTMP3,   LCM,         0.0)
  call AllocateDSI( WVTMP4,   LCM,         0.0)

  call AllocateDSI( WWW,      LCM,        -KCM,       0.0)
                                                     
  call AllocateDSI( Z,        LCM,         -KCM,      0.0)
  call AllocateDSI( ZBR,      LCM,         0.0)
  call AllocateDSI( ZBRE,     LCM,         0.0)
  call AllocateDSI( ZSRE,     LCM,         0.0)
  call AllocateDSI( ZEQ,      LCM,         0.0)
  call AllocateDSI( ZEQD,     LCM,         0.0)
  call AllocateDSI( ZEQDI,    LCM,         0.0)
  call AllocateDSI( ZEQI,     LCM,         0.0)
  call AllocateDSI( ZZ,       LCM,        -KCM,       0.0)
  call AllocateDSI( ZZC,     -KCM,         LCM,       0.0)

  ! *** NEXT SECTION ALLOCATES USAGE DEPENDENT VARIABLES
  !IF( ISDRY > 0 ) allocate(TIMEDRY(LCM, 0.0)

  if( ISINWV == 1 )then
    call AllocateDSI( CFLCAC,   LCM,       KCM,       0.0)
    call AllocateDSI( CFLUUU,   LCM,       KCM,       0.0)
    call AllocateDSI( CFLVVV,   LCM,       KCM,       0.0)
    call AllocateDSI( CFLWWW,   LCM,      -KCM,       0.0)
  endif

  if( ISVHEAT > 0 )then
    call AllocateDSI( LSVHTWINDE, LCM, .false.)
    call AllocateDSI( LSVHTWINDC, LCM, .false.)
    call AllocateDSI( SVREVC,     LCM,     0.0)
    call AllocateDSI( SVRCHC,     LCM,     0.0)
  endif
  call AllocateDSI( SVKEBACK, LCM, 0.0)

  if( ISHDMF >= 1 )then
    call AllocateDSI( FMDUX,    LCM,  KCM,  0.0)
    call AllocateDSI( FMDUY,    LCM,  KCM,  0.0)
    call AllocateDSI( FMDVX,    LCM,  KCM,  0.0)
    call AllocateDSI( FMDVY,    LCM,  KCM,  0.0)
    if( IS2TIM == 0 ) Call AllocateDSI( H1C, LCM, 0.0)
  endif

  if( ISVEG > 0 )then
    call AllocateDSI( FXVEG,     LCM,      KCM,   0.0)
    call AllocateDSI( FXVEGE,    LCM,      0.0)   
    call AllocateDSI( FYVEG,     LCM,      KCM,   0.0)
    call AllocateDSI( FYVEGE,    LCM,      0.0)
    call AllocateDSI( LVEG,      LCM,  .false.)
    call AllocateDSI( NVEGSERV, -NVEGTPM,    0)
    call AllocateDSI( ALPVEG,   -NVEGTPM,  0.0)
    call AllocateDSI( BDLPSQ,   -NVEGTPM,  0.0)
    call AllocateDSI( BPVEG,    -NVEGTPM,  0.0)
    call AllocateDSI( HPVEG,    -NVEGTPM,  0.0)
    call AllocateDSI( PVEGZ,    -NVEGTPM,  0.0)
    call AllocateDSI( RDLPSQ,   -NVEGTPM,  0.0)
    call AllocateDSI( SCVEG,    -NVEGTPM,  0.0)
  endif

  if( ISWAVE > 0 )then
    call AllocateDSI( FXWAVE,   LCM,  KCM, 0.0)
    call AllocateDSI( FYWAVE,   LCM,  KCM, 0.0)
    do L = 1, LCM
      call AllocateDSI( WV(L).DISSIPA, KCM, 0.0)
    enddo
  endif

  ! *** FOODCHAIN MODELLING OPTIONS
  if( ISTPOCB == 4 )then
    call AllocateDSI( PFPOCB,   LCM, KBM, 0.0)
    call AllocateDSI( FPOCB,    LCM, KBM, 0.0)
  else
    call AllocateDSI( PFPOCB,   1, 1, 0.0)
    call AllocateDSI( FPOCB,    1, 1, 0.0)
  endif

  ! *** TOXIC TRANSPORT VARIABLES
  call AllocateDSI( ISDIFBW,         NTXM,   0)
  call AllocateDSI( ISPMXZ,          NTXM,   0)
  call AllocateDSI( ISTOC,           NTXM,   0)
  call AllocateDSI( ITOXBU,          NTXM,   0)
  call AllocateDSI( ITOXKIN,         10,     NTXM,   0)
  call AllocateDSI( ITOXWU,          NTXM,   0)
  call AllocateDSI( ITXBDUT,         NTXM,   0)
  call AllocateDSI( ITXINT,          NTXM,   0)
  call AllocateDSI( ITXPARW,         NSTM2,  NTXM,   0)
  call AllocateDSI( ITXPARWC,        2,      NTXM,   0)
                                            
  call AllocateDSI( CONPARW,         NSTM2,  NTXM,  0.0)
  call AllocateDSI( CONPARWC,        2,      NTXM,  0.0)
  call AllocateDSI( DIFTOX,          NTXM,   0.0)
  call AllocateDSI( DIFTOXS,         NTXM,   0.0)
  call AllocateDSI( DIFTOXBW,        LCM,    NTXM,  0.0)
  call AllocateDSI( DPDIFTOX,        NTXM,   0.0)   
  call AllocateDSI( FPOCBST,         NSTM2,  NTXM,  0.0)
  call AllocateDSI( FPOCWST,         NSTM2,  NTXM,  0.0)
  call AllocateDSI( PDIFTOX,         NTXM,   0.0)
  call AllocateDSI( RKTOXP,          NTXM,   0.0)
  call AllocateDSI( SKTOXP,          NTXM,   0.0)
  call AllocateDSI( TOXINTB,         NTXM,   0.0)
  call AllocateDSI( TOXINTW,         NTXM,   0.0)
  
  call AllocateDSI( TOXPARB,         LCM,    NSTM2,  NTXM,   0.0)
  call AllocateDSI( TOXPARBC,        2,      NTXM,    0.0)
  call AllocateDSI( TOXPARW,         LCM,    NSTM2,  NTXM,   0.0)
  call AllocateDSI( TOXPARWC,        2,      NTXM,    0.0)
  call AllocateDSI( TOX,             LCM,    KCM,     NTXM,  0.0)
  call AllocateDSI( TOX1,            LCM,    KCM,     NTXM,  0.0)
  allocate(TOXS(NTXM))
  
  if( ISTRAN(5) >= 1 .or. NTOX > 0 )then
    call AllocateDSI( NSP2,          NTXM,   0)
    call AllocateDSI( STPOCW,        LCM,    KCM,     0.0)
    call AllocateDSI( TADFLUX,       LCM,    NTXM,    0.0)
    call AllocateDSI( TOXB,          LCM,    KBM,     NTXM,  0.0)
    call AllocateDSI( TOXB1,         LCM,    KBM,     NTXM,  0.0)
    call AllocateDSI( TOXBINIT,      LCM,    KBM, NTXM, 0.0)
    call AllocateDSI( TOXCDFB,       LCM,    KBM, NTXM, 0.0)
    call AllocateDSI( TOXCDFW,       LCM,    KCM, NTXM, 0.0)
    call AllocateDSI( TOXF,          LCM,   -KCM,     NTXM,  0.0)    ! *** TOXIC CONTAMINANT SETTLING AND BED EXCHANGE FLUX  (G/M2/S)
    call AllocateDSI( TOXFB,         LCM,    NTXM,    0.0)           ! *** TOXIC CONTAMINANT EROSIONAL FLUX AT THE BED/WATER INTERFACE  (G/M2/S)
    call AllocateDSI( TOXFBL,        LCM,    NTXM,    0.0)           ! *** TOXIC CONTAMINANT FLUX AT THE BED/WATER INTERFACE DUE TO BEDLOAD   (G/M2/S)
    call AllocateDSI( TOXFBEBKB,     LCM,    NTXM,    0.0)
    call AllocateDSI( TOXFBECHB,     LCM,    NTXM,    0.0)
    call AllocateDSI( TOXFBECHW,     LCM,    NTXM,    0.0)
    call AllocateDSI( TOXFDFB,       LCM,    KBM,     NTXM, 0.0)
    call AllocateDSI( TOXFDFW,       LCM,    KCM,     NTXM, 0.0)
    call AllocateDSI( TOXINIT,       LCM,    KCM,     NTXM, 0.0)
    call AllocateDSI( TOXPFB,        LCM,    KBM,     NSTM2+2, NTXM,  0.0)
    call AllocateDSI( TOXPFTB,       LCM,    KBM,     NTXM,    0.0)  
    call AllocateDSI( TOXPFTW,       LCM,    KCM,     NTXM,    0.0)  
    call AllocateDSI( TOXPFW,        LCM,    KCM,     NSTM2+2, NTXM,  0.0)
    call AllocateDSI( TOXTMP,        LCM,  -(KBM+1),  NTXM, 0.0)
  endif   ! *** END OF TOXIC VARIABLE DECLARATIONS
  
  ! *** SED VARIABLES
  call AllocateDSI( SED,      LCM,  KCM,  NSEDS2,  0.0)
  call AllocateDSI( SED1,     LCM,  KCM,  NSEDS2,  0.0)
  call AllocateDSI( SDF,      LCM,  KCM,   NSEDS,  0.0)
  if( ISTRAN(6) >= 1 .or. NSED > 0 )then
    call AllocateDSI( SED3DMAX,      NSEDS2,  0.0)
    call AllocateDSI( SED3DMIN,      NSEDS2,  0.0)
    call AllocateDSI( SEDBA,         LCM,     NMAX,   0.0)
    call AllocateDSI( SEDBINIT,      LCM,     KBM,    NMAX, 0.0)
    call AllocateDSI( SEDF,          LCM,    -KCM,   NSEDS2, 0.0)      ! *** VERTICAL SED FLUX ACROSS CELL LAYER FACES (G/M2/S)
    call AllocateDSI( SEDINIT,       LCM,     KCM,   NSEDS2, 0.0)
    call AllocateDSI( SED2,          LCM,     KCM,   NSEDS2, 0.0)
  endif   ! *** END OF SED VARIABLE DECLARATIONS

  ! *** SND VARIABLES
  call AllocateDSI( SND,      LCM, KCM, NSNM, 0.0)
  call AllocateDSI( SND1,     LCM, KCM, NSNM, 0.0)
  if( NSND > 0 )then
    call AllocateDSI( SBDLDA,        NSNM,    0.0)
    call AllocateDSI( SBDLDB,        NSNM,    0.0)
    call AllocateDSI( SBDLDG1,       NSNM,    0.0)
    call AllocateDSI( SBDLDG2,       NSNM,    0.0)
    call AllocateDSI( SBDLDG3,       NSNM,    0.0)
    call AllocateDSI( SBDLDG4,       NSNM,    0.0)
    call AllocateDSI( SBDLDP,        NSNM,    0.0)
  endif
  if( ISTRAN(7) >= 1 )then
    call AllocateDSI( SNDBA,         LCM,     NSNM,   0.0)
    call AllocateDSI( SNDBINIT,      LCM,     KBM,    NSNM,  0.0)
    call AllocateDSI( SNDEQ,         LCM,     0.0)          
    call AllocateDSI( SNDEQSAV,      LCM,     NSNM,   0.0)  
    call AllocateDSI( SNDEQB,        LCM,     0.0)          
    call AllocateDSI( SNDF,          LCM,    -KCM,    NSNM,  0.0)      ! *** VERTICAL SND FLUX ACROSS CELL LAYER FACES (G/M2/S)
    call AllocateDSI( SNDFBL,        LCM,     NSNM,   0.0)             ! *** BED/WATER INTERFACE SND FLUX DUE TO BEDLOAD (G/M2/S)
    call AllocateDSI( SNDINIT,       LCM,     KCM,    NSNM,  0.0)
    call AllocateDSI( SNDS,          LCM,     KCM,    NSNM,  0.0)
  endif   ! *** END OF SND VARIABLE DECLARATIONS

  ! *** EITHER SED OR SND
  call AllocateDSI( IBLTAUC,       NSEDS2,    0)
  call AllocateDSI( ISEDWU,        NSEDS2,    0)
  call AllocateDSI( IROUSE,        NSEDS2,    0)
  call AllocateDSI( ISBDLD,        NSEDS2,    0)
  call AllocateDSI( ISEDSCOR,      NSEDS2,    0)
  call AllocateDSI( ISLTAUC,       NSEDS2,    0)
  call AllocateDSI( ISNDEQ,        NSEDS2,    0)
  call AllocateDSI( ISNDM1,        NSEDS2,    0)
  call AllocateDSI( ISNDM2,        NSEDS2,    0)
  call AllocateDSI( IWRSPB,        NSEDS2,    0)
  call AllocateDSI( ISPROBDEP,     NSEDS2,    0)
                                         
  call AllocateDSI( COSEDHID,      NSTM2,  0.0)
  call AllocateDSI( RSNDM,         NSTM2,  0.0)
  call AllocateDSI( SDEN,          NSTM2,  0.0)               ! ***  SPECIFIC VOLUME OF EACH SEDIMENT CLASS (M3/G)
  call AllocateDSI( SEDN,          NSTM2,  0.0)
  call AllocateDSI( SEDO,          NSTM2,  0.0)
  call AllocateDSI( SEDBO,         NSTM2,  0.0)
  call AllocateDSI( SEDDIA,        NSTM2,  0.0)
  call AllocateDSI( SEXP,          NSTM2,  0.0)
  call AllocateDSI( SSG,           NSTM2,  0.0)
  call AllocateDSI( TAUD,          NSTM2,  0.0)
  call AllocateDSI( TAUN,          NSTM2,  0.0)
  call AllocateDSI( TAUR,          NSTM2,  0.0)
  call AllocateDSI( TCSHIELDS,     NSTM2,  0.0)
  call AllocateDSI( TEXP,          NSTM2,  0.0)
  call AllocateDSI( VDRDEPO,       NSTM2,  0.0)
  call AllocateDSI( WRSPO,         NSTM2,  0.0)
  call AllocateDSI( WSEDO,         NSTM2,  0.0)
  call AllocateDSI( VDRRSPO,       NSTM2,  0.0)
  
  if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
    call AllocateDSI( ACOEF,       LCM,      -KBM,     0.0)
    call AllocateDSI( ALOW,        LCM,      -(KBM+1), NTXM, 0.0)
    call AllocateDSI( BDENBED,     LCM,       KBM,     0.0)
    call AllocateDSI( BDENBED1,    LCM,       KBM,     0.0)
    call AllocateDSI( BDENBEDA,    LCM,       0.0)     
    call AllocateDSI( BDENBEDA1,   LCM,       0.0)     
    call AllocateDSI( BEDBINIT,    LCM,       KBM,     0.0)
    call AllocateDSI( BEDDINIT,    LCM,       KBM,     0.0)
    call AllocateDSI( BEDLINIT,    LCM,       KBM,     0.0)
    call AllocateDSI( CBEDTOTAL,   LCM,       0.0)     
    call AllocateDSI( COEFK,       LCM,       KBM,     0.0)
    call AllocateDSI( COEFSK,      LCM,       KBM,     0.0)
    call AllocateDSI( CQBEDLOADX,  LCM,       NSNM,    0.0)
    call AllocateDSI( CQBEDLOADY,  LCM,       NSNM,    0.0)
    call AllocateDSI( CSHIELDS50,  LCM,       0.0)     
    call AllocateDSI( CUPP,        LCM,      -(KBM+1), NTXM, 0.0)
    call AllocateDSI( DSTRSE,      LCM,       KBM,     0.0)
    call AllocateDSI( DZBTR,       LCM,       KBM,     0.0)
    call AllocateDSI( DZBTR1,      LCM,       KBM,     0.0)
    call AllocateDSI( FRACCOH,     LCM,       KBM,     0.0)
    call AllocateDSI( FRACNON,     LCM,       KBM,     0.0)
    call AllocateDSI( HBED,        LCM,       KBM,     0.0)
    call AllocateDSI( HBED1,       LCM,       KBM,     0.0)
    call AllocateDSI( HYDCN,       LCM,       KBM,     0.0)
                                             
    call AllocateDSI( ISEDBU,      NSEDS2,     0)
    call AllocateDSI( ISNDBU,      NSNM,      0)
    call AllocateDSI( ISDBLDIR,    LCM,       0)
                                             
    call AllocateDSI( PEXP,        LCM,       NSNM,             0.0)
    call AllocateDSI( PHID,        LCM,       NSNM,             0.0)
    call AllocateDSI( QSBDLDIN,    LCM,       NSNM,             0.0)                ! *** BEDLOAD FLUX IN  DUE TO FIXING BOUNDARY CONDITION CELLS TO ZERO BEDLOAD DELTA  (G/S)
    call AllocateDSI( QSBDLDOT,    LCM,       max(NSNM, NSEDS2), 0.0)    ! *** BEDLOAD FLUX OUT DUE TO FIXING BOUNDARY CONDITION CELLS TO ZERO BEDLOAD DELTA  (G/S)
    call AllocateDSI( QSBDLDX,     LCM,       max(NSNM, NSEDS2), 0.0)    ! *** U FACE SND FLUX DUE TO BEDLOAD (G/S)
    call AllocateDSI( QSBDLDY,     LCM,       max(NSNM, NSEDS2), 0.0)    ! *** V FACE SND FLUX DUE TO BEDLOAD (G/S)
    call AllocateDSI( QSBDLDP,     LCM,       0.0)                      ! *** CELL CENTER BED LOAD TRANSPORT RATE (G/S)
    call AllocateDSI( QSBDTOP,     LCM,       0.0)                      ! *** VOLUME OF SEDIMENT BED EXCHANGE (M/S)
    call AllocateDSI( QSEDBED,     LCM,       KBM,              NSTM2,  0.0)
    call AllocateDSI( QSEDBED1,    LCM,       KBM,              NSTM2,  0.0)
    call AllocateDSI( QSEDBEDA,    LCM,       NSTM2,            0.0)
    call AllocateDSI( QSEDBEDA1,   LCM,       NSTM2,            0.0)
                                             
    call AllocateDSI( SEDDIA50,    LCM,       KBM, 0.0)
    call AllocateDSI( SEDDIA90,    LCM,       KBM, 0.0)
                                             
    call AllocateDSI( SEDB,        LCM,       KBM,        NSEDS2, 0.0)
    call AllocateDSI( SEDB1,       LCM,       KBM,        NSEDS2, 0.0)
    call AllocateDSI( SEDBT,       LCM,       KBM,        0.0)
    call AllocateDSI( SEDFPA,      LCM,       NSEDS2,      0.0)
    call AllocateDSI( SEDT,        LCM,       KCM,        0.0)
                                                         
    call AllocateDSI( SNDB,        LCM,       KBM,        NSNM, 0.0)
    call AllocateDSI( SNDB1,       LCM,       KBM,        NSNM, 0.0)
    call AllocateDSI( SNDBT,       LCM,       KBM,        0.0)
    call AllocateDSI( SNDFPA,      LCM,       NSNM,       0.0)
    call AllocateDSI( SNDT,        LCM,       KCM,        0.0)
                                             
    call AllocateDSI( TAUBSED,     LCM,       0.0)                ! *** Grain separated density normalized bed shear stress (m2/s2)
    call AllocateDSI( TAUBSND,     LCM,       0.0)                ! *** Grain separated density normalized bed shear stress (m2/s2)
    call AllocateDSI( TAUCRCOH,    LCM,       KBM,     0.0)
    call AllocateDSI( TAURB,       LCM,       KBM,     0.0)
    call AllocateDSI( TAUDS,       LCM,       0.0)     
    call AllocateDSI( TAUNS,       LCM,       KBM,     0.0)
    call AllocateDSI( TAURS,       LCM,       KBM,     0.0)
    call AllocateDSI( TEXPS,       LCM,       KBM,     0.0)
                                                      
    call AllocateDSI( VDRBED,      LCM,       KBM,     0.0)
    call AllocateDSI( VDRBED1,     LCM,       KBM,     0.0)
    call AllocateDSI( VDRBED0,     LCM,       KBM,     0.0)
    call AllocateDSI( VDRBEDA,     LCM,       0.0)     
    call AllocateDSI( VDRBEDSED,   LCM,       KBM,     0.0)
    call AllocateDSI( VDRBEDSND,   LCM,       KBM,     0.0)
    call AllocateDSI( VFRBED,      LCM,       KBM,     NSTM2,  0.0)
    call AllocateDSI( VFRBED1,     LCM,       KBM,     NSTM2,  0.0)
    call AllocateDSI( WSETA,       LCM,      -KCM,     NSTM2,  0.0)
    call AllocateDSI( WRSPB,       LCM,       KBM,     0.0)
    call AllocateDSI( WRSPBA,      LCM,       0.0)     
    call AllocateDSI( WRSPS,       LCM,       KBM,     0.0)
    call AllocateDSI( ZBEDC,       LCM,       KBM,     0.0)
    call AllocateDSI( ZBEDG,       LCM,      -KBM,     0.0)
    call AllocateDSI( ZBEDGT,      LCM,       0.0)     
    if( ISBEDSTR == 3 )then                           
      call AllocateDSI( ZBRSED,    LCM,       0.0)     
    else                                              
      call AllocateDSI( ZBRSED,    1,         0.0)     
    endif                                             
    call AllocateDSI( ZELBED,      LCM,       KBM,     0.0)
    call AllocateDSI( ZELBED1,     LCM,       KBM,     0.0)
    call AllocateDSI( ZELBEDA,     LCM,       0.0)
    call AllocateDSI( ZELBEDA1,    LCM,       0.0)
                                             
    call AllocateDSI( PARTMIXZ,    LCM,       KBM,     0.0)
    call AllocateDSI( PMXDEPTH,    NPMXPTSM,  NPMXZM,  0.0)
    call AllocateDSI( PMXCOEF,     NPMXPTSM,  NPMXZM,  0.0)
    call AllocateDSI( PORBED,      LCM,       KBM,     0.0)
    call AllocateDSI( PORBED1,     LCM,       KBM,     0.0)
    call AllocateDSI( QCOEF,       LCM,      -KBM,     0.0)
    call AllocateDSI( QWTRBED,     LCM,      -KBM,     0.0)                ! *** VOLUME OF WATER EXCHANGE BETWEEN SEDIMENT LAYERS (M/S)
    call AllocateDSI( RRHS,        LCM,      -(KBM+1), NTXM,  0.0)
    call AllocateDSI( SDENAVG,     LCM,       KBM,     0.0)
    call AllocateDSI( SGSM1,       LCM,       KBM,     0.0)
    call AllocateDSI( SIGPHI,      LCM,       KBM,     0.0)
    call AllocateDSI( STDOCB,      LCM,       KBM,     0.0)
    call AllocateDSI( STDOCW,      LCM,       KCM,     0.0)
    call AllocateDSI( STFPOCB,     LCM,       KBM,     NSTM2,  0.0)
    call AllocateDSI( STFPOCW,     LCM,       KCM,     NSTM2,  0.0)
    call AllocateDSI( STPOCB,      LCM,       KBM,     0.0)
    call AllocateDSI( STRSE,       LCM,       KBM,     0.0)
                                  
    ! *** BANK EROSION            
    call AllocateDSI( IBANKBE,     NBEPAIRM,  0)
    call AllocateDSI( JBANKBE,     NBEPAIRM,  0)
    call AllocateDSI( ICHANBE,     NBEPAIRM,  0)
    call AllocateDSI( JCHANBE,     NBEPAIRM,  0)
    call AllocateDSI( NBESERN,     NBEPAIRM,  0)
    call AllocateDSI( MBESER,      NBESERM,   0)
    call AllocateDSI( MBETLAST,    NBESERM,   0)
    call AllocateDSI( LPMXZ,       LCM,       0)
  endif
  allocate( SEDS(NSEDS2) )

  ! *** CALMMT/RESIDUAL/MODEL LINKAGES VARIABLES
  if( ISSSMMT > 0 .or. ISWASP > 0 )then
    call AllocateDSI( EVPGLPF,     LCM,       0.0)                     ! *** MMT - AVERAGE TRANSPIRATION RATE (M3/S)
    call AllocateDSI( EVPSLPF,     LCM,       0.0)                     ! *** MMT - AVERAGE EVAPORATION RATE (M3/S)
    call AllocateDSI( GWLPF,       LCM,       0.0)                     ! *** MMT - AVERAGE GROUNDWATER ELEVATION BELOW CELL (M)
    call AllocateDSI( HLPF,        LCM,       0.0)                     ! *** MMT - AVERAGE DEPTH
    call AllocateDSI( QSUMELPF,    LCM,       0.0)
    call AllocateDSI( QSUMLPF,     LCM,       KCM,  0.0)
    call AllocateDSI( RAINLPF,     LCM,       0.0)                     ! *** MMT - ACCUMULATED RAINFALL (M3/S)
    call AllocateDSI( RINFLPF,     LCM,       0.0)                     ! *** ISGWIE > - MMT -ACCUMULATED INFILTRATION RATE (M3/S)
    call AllocateDSI( SALLPF,      LCM,       KCM,  0.0)
    call AllocateDSI( SFLLPF,      LCM,       KCM,  0.0)
    call AllocateDSI( TEMLPF,      LCM,       KCM,  0.0)
    call AllocateDSI( UELPF,       LCM,       0.0)                     ! *** MMT - AVERAGE DEPTH AVERAGED U VELOCITY
    call AllocateDSI( VELPF,       LCM,       0.0)                     ! *** MMT - AVERAGE DEPTH AVERAGED V VELOCITY
    call AllocateDSI( WLPF,        LCM,      -KCM,  0.0)
    call AllocateDSI( WTLPF,       LCM,       KCM,  0.0)

    if( ISTRAN(5) > 0 )then
      call AllocateDSI( TOXBLPF,   LCM,       KBM,  NTXM,     0.0)
      call AllocateDSI( TOXLPF,    LCM,       KCM,  NTXM,     0.0)
      call AllocateDSI( TXPFLPF,   LCM,       KCM,  NSTM2+2,  NTXM,  0.0)  ! ***
    endif
    if( ISTRAN(6) > 0 )then
      call AllocateDSI( SEDBLPF,   LCM,       KBM,  NSEDS2,    0.0)
      call AllocateDSI( SEDBTLPF,  LCM,       KBM,  0.0)      
      call AllocateDSI( SEDLPF,    LCM,       KCM,  NSEDS2,    0.0)
      call AllocateDSI( SEDTLPF,   LCM,       KCM,  0.0)
    endif
    if( ISTRAN(7) > 0 )then
      call AllocateDSI( SNDBLPF,   LCM,       KBM,  NSNM,     0.0)
      call AllocateDSI( SNDBTLPF,  LCM,       KBM,  0.0)      
      call AllocateDSI( SNDLPF,    LCM,       KCM,  NSNM,     0.0)
      call AllocateDSI( SNDTLPF,   LCM,       KCM,  0.0)
    endif
  endif

  ! *** WATER QUALITY VARIABLES
  call AllocateDSI( WQKEB,    1, 0.0)
  if( ISTRAN(8) >= 1 )then
    call WQ_Allocate
    if( IWQBEN > 0 ) Call SD_Allocate
    if( IWQZPL > 0 ) Call Zoo_Allocate
  endif   ! *** END OF WQ VARIABLE DECLARATIONS

  if( ISWQFLUX == 1 )then
    call AllocateDSI( ABEFF,       LCM,       KCM,     0.0)
    call AllocateDSI( ABLPF,       LCM,       KCM,     0.0)
    call AllocateDSI( ACCWFLD,     LCM,       2,       0.0)
    call AllocateDSI( AHULPF,      LCM,       KCM,     0.0)
    call AllocateDSI( AHVLPF,      LCM,       KCM,     0.0)
    call AllocateDSI( DYELPF,      LCM,       KCM,     NDYM,  0.0)
                                                     
    call AllocateDSI( UHLPF,       LCM,       KCM,     0.0)
    call AllocateDSI( UIRT,        LCM,       KCM,     0.0)
    call AllocateDSI( ULPF,        LCM,       KCM,     0.0)
    call AllocateDSI( UTLPF,       LCM,       KCM,     0.0)
    call AllocateDSI( UVPT,        LCM,       KCM,     0.0)
                                                     
    call AllocateDSI( VHLPF,       LCM,       KCM,     0.0)
    call AllocateDSI( VIRT,        LCM,       KCM,     0.0)
    call AllocateDSI( VLPF,        LCM,       KCM,     0.0)
    call AllocateDSI( VTLPF,       LCM,       KCM,     0.0)
    call AllocateDSI( VVPT,        LCM,       KCM,     0.0)
                                                     
    call AllocateDSI( WIRT,        LCM,       KCM,     0.0)
  endif                           
                                  
  if( ISWASP > 0 )then            
    call AllocateDSI( KCEFDC,      LCM*KCM,    0)
    call AllocateDSI( LCEFDC,      LCM*KCM,    0)
    call AllocateDSI( KFEFDC,      LCM*KCM,    0)
    call AllocateDSI( LFEFDC,      KCM*LCM,    0)
    call AllocateDSI( NCHNC,       KCM*LCM,    0)
    call AllocateDSI( LCHNC,     2*KCM*LCM,   10,       0)
  else                          
    call AllocateDSI( KCEFDC,      1,         0)
    call AllocateDSI( LCEFDC,      1,         0)
    call AllocateDSI( KFEFDC,      1,         0)
    call AllocateDSI( LFEFDC,      1,         0)
    call AllocateDSI( NCHNC,      -1,         0)
    call AllocateDSI( LCHNC,      -1,         1,        0)
  endif
  
  ! Begin SEDZLJ variables
  call AllocateDSI( FWVTP,          LCM,    0.0)
  call AllocateDSI( FWDIR,          LCM,    16,   0.0)
  call AllocateDSI( LWDIR,          LCM,    16,  -(ICM+JCM), 0)
  if( LSEDZLJ )then
    call AllocateDSI( ALPHA_PX,     LCM,    0.0)
    call AllocateDSI( ALPHA_PY,     LCM,    0.0)
    call AllocateDSI( ALPHA_RX,     LCM,    NSEDS2,   0.0)
    call AllocateDSI( ALPHA_RY,     LCM,    NSEDS2,   0.0)
    call AllocateDSI( BLFLAG,       LCM,    NSEDS2,     0)
    call AllocateDSI( BULKDENS,     KB,     LCM,      0.0)
    call AllocateDSI( D50,          NSEDS2, 0.0)
    call AllocateDSI( D50AVG,       LCM,    0.0)
    call AllocateDSI( DEP_SED_FLX,  LCM,    NSEDS2,   0.0)
    call AllocateDSI( DISTAR,       NSEDS2, 0.0)
    call AllocateDSI( DZBL,         LCM,    NSEDS2,   0.0)
    !ALLOCATE(DZBL_LAST(LCM, NS    CM2, 0.0)
    call AllocateDSI( DWS,          NSEDS2,  0.0)
    call AllocateDSI( DWSIN,        NSEDS2,  0.0)
    call AllocateDSI( ESUS,         LCM,    0.0)
    call AllocateDSI( ERO_SED_FLX,  LCM,    NSEDS2,   0.0)
    call AllocateDSI( HPCM,         LCM,    0.0)
    call AllocateDSI( KPART,        NTXM,   0.0)
    call AllocateDSI( LAYERACTIVE,  KB,     LCM,      0)
    call AllocateDSI( NCORENO,      IC,     JC,       0)
    call AllocateDSI( PERSED,       NSEDS,  KB, LCM,  0.0)
    call AllocateDSI( PSUS,         LCM,    NSEDS,    0.0)
    call AllocateDSI( QBLFLUX,      LCM,    NSEDS,    0.0)
    call AllocateDSI( SCND,         NSICM,  0.0)      
    call AllocateDSI( SH_SCALE,     LCM,    0.0)      
    call AllocateDSI( SSGI,         NSEDS2, 0.0)      
    call AllocateDSI( STWVHT,       LCM,    200,      0.0)
    call AllocateDSI( STWVTP,       LCM,    200,      0.0)
    call AllocateDSI( STWVDR,       LCM,    200,      0.0)
    call AllocateDSI( TAU,          LCM,    0.0)      
    call AllocateDSI( TAUCOR,       KB+1,   LCM,      0.0)
    call AllocateDSI( TAUCRITE,     NSICM,  0.0)
    call AllocateDSI( TCRE,         NSEDS,  0.0)
    call AllocateDSI( TCRSUS,       NSEDS,  0.0)
    call AllocateDSI( TRANS,        LCM,    NSEDS2,   0.0)
    call AllocateDSI( TSED,         KB,     LCM,      0.0)
    call AllocateDSI( TSED0,        KB,     LCM,      0.0)
    call AllocateDSI( TSET0T,       LCM,    0.0)      
    call AllocateDSI( TSEDT,        LCM,    0.0)      
    call AllocateDSI( USW,          LCM,    NSEDS,    0.0)     
    if( ICALC_BL > 0 )then                           
      call AllocateDSI( BLVEL,      LCM,    NSEDS,    0.0)
      call AllocateDSI( CBL,        LCM,    NSEDS,    0.0)
      call AllocateDSI( DBL,        LCM,    NSEDS,    0.0)
      call AllocateDSI( EBL,        LCM,    NSEDS,    0.0)
      call AllocateDSI( UBL,        LCM,    NSEDS,    0.0)
      call AllocateDSI( VBL,        LCM,    NSEDS,    0.0)
      if( ISTRAN(5) > 0 )then                        
        call AllocateDSI( CBLTOX,   LCM,    NTXM,     0.0)
        call AllocateDSI( CBLTXCON, LCM,    NSEDS,    NTXM,   0.0)
      endif
    endif
  endif
  ! End SEDZLJ variables

  ! *** OMP DECLARATIONS
  call AllocateDSI( CD,         LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( FQC,        LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( FUHVD,      LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( FVHVD,      LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( FWVV,       LCM,      -KCM,   NTHREADS,   0.0)
  call AllocateDSI( DUU,        LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( DVV,        LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( POS,        LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( UUUU,       LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( VVVV,       LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( WWWW,       LCM,      -KCM,   NTHREADS,   0.0)
  call AllocateDSI( CONQ,       LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( CMAX,       LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( CMIN,       LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( CONTD,      LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( CONTMN,     LCM,       KCM,   NTHREADS,   0.0)
  call AllocateDSI( CONTMX,     LCM,       KCM,   NTHREADS,   0.0)

  ! Begin MHK variables SCJ
  call AllocateDSI( IJLTURB,    TCOUNT,    3,     0)
  call AllocateDSI( CTMHK,      MHKTYP,    0.0)
  call AllocateDSI( CDSUP,      MHKTYP,    0.0)
  call AllocateDSI( DENMHK,     MHKTYP,    0.0)            ! *** DENSITY OF MHK DEVICES (#/CELL)
  call AllocateDSI( ESUP,       TCOUNT,    LCM,   0.0)     ! *** zero the allocatable arrays for MHK energy generation
  call AllocateDSI( EMHK,       TCOUNT,    LCM,   0.0)     ! *** zero the allocatable arrays for support energy dissipatio
  call AllocateDSI( FXMHK,      LCM,       KCM,   0.0)     ! *** X"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK DEVICE
  call AllocateDSI( FXMHKE,     LCM,       0.0)            ! *** COLUMN X"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK DEVICE
  call AllocateDSI( FXSUP,      LCM,       KCM,   0.0)     ! *** Y"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK DEVICE
  call AllocateDSI( FXSUPE,     LCM,       0.0)            ! *** COLUMN Y"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK DEVICE
  call AllocateDSI( FYMHK,      LCM,       KCM,   0.0)     ! *** X"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK SUPPORT
  call AllocateDSI( FYMHKE,     LCM,       0.0)            ! *** COLUMN X"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK SUPPORT
  call AllocateDSI( FYSUP,      LCM,       KCM,   0.0)     ! *** Y"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK SUPPORT
  call AllocateDSI( FYSUPE,     LCM,       0.0)            ! *** COLUMN Y"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK DEVICE
  call AllocateDSI( PMHK,       LCM,       KCM,   0.0)     ! *** array that accumulates layer-wise MHK power in x direction
  call AllocateDSI( PSUP,       LCM,       KCM,   0.0)     ! *** array that accumulates layer-wise MHK support power in x direction 
  call AllocateDSI( VMAXCUT,    MHKTYP,    0.0)   
  call AllocateDSI( VMINCUT,    MHKTYP,    0.0)   
  call AllocateDSI( WIDTHMHK,   MHKTYP,    0.0)   
  call AllocateDSI( WIDTHSUP,   MHKTYP,    0.0)   
  call AllocateDSI( ZMAXMHK,    MHKTYP,    LCM,   0.0)
  call AllocateDSI( ZMINMHK,    MHKTYP,    LCM,   0.0)
  call AllocateDSI( ZMAXSUP,    MHKTYP,    LCM,   0.0)
  call AllocateDSI( ZMINSUP,    MHKTYP,    LCM,   0.0)

  ! *** ALLOCATE MEMORY FOR VARIABLE TO STORE CONCENTRATIONS AT OPEN BOUNDARIES
  call AllocateDSI( FBESER,         NBEPAIRM,  0.0)
  call AllocateDSI( BESERT,         NBESERM,   0.0)
  call AllocateDSI( FWCBESERT,      NBESERM,   0.0)
  call AllocateDSI( TCBESER,        NBESERM,   0.0)
  call AllocateDSI( TABESER,        NBESERM,   0.0)
  call AllocateDSI( TBESER,         NDBESER,   NBESERM,  0.0)
  call AllocateDSI( BESER,          NDBESER,   NBESERM,  0.0)
  call AllocateDSI( FWCBESER,       NDBESER,   NBESERM,  0.0)
  call AllocateDSI( QSBDTOPBEBKB,   LCM,       0.0)
  call AllocateDSI( QSBDTOPBECHB,   LCM,       0.0)
  call AllocateDSI( QSBDTOPBECHW,   LCM,       0.0)
  call AllocateDSI( QWBDTOPBEBKB,   LCM,       0.0)
  call AllocateDSI( QWBDTOPBECHB,   LCM,       0.0)
  call AllocateDSI( QWBDTOPBECHW,   LCM,       0.0)
  call AllocateDSI( SEDFBEBKB,      LCM,       NSEDS2,    0.0)
  call AllocateDSI( SEDFBECHB,      LCM,       NSEDS2,    0.0)
  call AllocateDSI( SEDFBECHW,      LCM,       NSEDS2,    0.0)
  call AllocateDSI( SNDFBEBKB,      LCM,       NSNM,     0.0)
  call AllocateDSI( SNDFBECHB,      LCM,       NSNM,     0.0)
  call AllocateDSI( SNDFBECHW,      LCM,       NSNM,     0.0)

  ! *** VARIABLES FOR QCTL NQCTYP = 3 & 4
  call AllocateDSI( LOWCHORDU,      NQCTL,     0.0)
  call AllocateDSI( LOWCHORDV,      NQCTL,     0.0)
  call AllocateDSI( NLOWCHORD,      NQCTL,     0)
  call AllocateDSI( SAVESUB,        2,         NQCTL,    0.0)
  call AllocateDSI( SAVESVB,        2,         NQCTL,    0.0)

  ! *** HYDRAULIC STRUCTURE EQUATIONS
  NS = max(NHYDST,1)
  call AllocateDSI( HS_LENGTH,    -NS,         0.0)
  call AllocateDSI( HS_XSTYPE,    -NS,           0)
  call AllocateDSI( HS_REVERSE,   -NS,           0)
  call AllocateDSI( HS_WIDTH,     -NS,         0.0)
  call AllocateDSI( HS_HEIGHT,    -NS,         0.0)
  call AllocateDSI( HS_MANN,      -NS,         0.0)
  call AllocateDSI( HS_ANGLE,     -NS,         0.0)
  call AllocateDSI( HS_USELEV,    -NS,         0.0)
  call AllocateDSI( HS_DSELEV,    -NS,         0.0)
  call AllocateDSI( HS_COEFF,     -NS,         4,   0.0)

  
  ! *** Macrophyte/vegetation
  call AllocateDSI( MAC_CELL, LCM, .FALSE.)
  
  ! *** SIGMA-Z HSTRUCT
  call AllocateDSI( KSZ,        LCM,     1)
  call AllocateDSI( KSZU,       LCM,     1)
  call AllocateDSI( KSZV,       LCM,     1)
  call AllocateDSI( BI1W,       LCM,   0.0)
  call AllocateDSI( BI2W,       LCM,   0.0)
  call AllocateDSI( BEW,        LCM,   0.0)
  call AllocateDSI( BI1S,       LCM,   0.0)
  call AllocateDSI( BI2S,       LCM,   0.0)
  call AllocateDSI( BES,        LCM,   0.0)
  if( IGRIDV > 0 .or. ISBLOCKED > 0 )then
    call AllocateDSI( BELVE,      LCM,   0.0)
    call AllocateDSI( BELVW,      LCM,   0.0)
    call AllocateDSI( BELVN,      LCM,   0.0)
    call AllocateDSI( BELVS,      LCM,   0.0)
    call AllocateDSI( BE,         LCM,   KCM,   0.0)
    call AllocateDSI( BW,         LCM,   KCM,   0.0)
    call AllocateDSI( BN,         LCM,   KCM,   0.0)
    call AllocateDSI( BS,         LCM,   KCM,   0.0)
    call AllocateDSI( BEE,        LCM,   0.0)
    call AllocateDSI( BEN,        LCM,   0.0)
    call AllocateDSI( BI1E,       LCM,   0.0)
    call AllocateDSI( BI1N,       LCM,   0.0)
    call AllocateDSI( BI2E,       LCM,   0.0)
    call AllocateDSI( BI2N,       LCM,   0.0)
    call AllocateDSI( HPE,        LCM,   0.0)
    call AllocateDSI( HPW,        LCM,   0.0)
    call AllocateDSI( HPN,        LCM,   0.0)
    call AllocateDSI( HPS,        LCM,   0.0)
    call AllocateDSI( KSZE,       LCM,     1)
    call AllocateDSI( KSZW,       LCM,     1)
    call AllocateDSI( KSZN,       LCM,     1)
    call AllocateDSI( KSZS,       LCM,     1)
    call AllocateDSI( SGZKE,      KCM,   LCM,   0.0)
    call AllocateDSI( SGZKW,      KCM,   LCM,   0.0)
    call AllocateDSI( SGZKN,      KCM,   LCM,   0.0)
    call AllocateDSI( SGZKS,      KCM,   LCM,   0.0)
    call AllocateDSI( SGZE,       LCM,   KCM,   0.0)
    call AllocateDSI( SGZW,       LCM,   KCM,   0.0)
    call AllocateDSI( SGZN,       LCM,   KCM,   0.0)
    call AllocateDSI( SGZS,       LCM,   KCM,   0.0)
    call AllocateDSI( ZE,         LCM,  -KCM,   0.0)
    call AllocateDSI( ZW,         LCM,  -KCM,   0.0)
    call AllocateDSI( ZN,         LCM,  -KCM,   0.0)
    call AllocateDSI( ZS,         LCM,  -KCM,   0.0)
    call AllocateDSI( ZZE,       -KCM,   LCM,   0.0)
    call AllocateDSI( ZZW,       -KCM,   LCM,   0.0)
    call AllocateDSI( ZZN,       -KCM,   LCM,   0.0)
    call AllocateDSI( ZZS,       -KCM,   LCM,   0.0)
  endif
  call AllocateDSI( DZIGSD4U, LCM,   KCM,   0.0)
  call AllocateDSI( DZIGSD4V, LCM,   KCM,   0.0)
                                          
  call AllocateDSI( DZC,      LCM,   KCM,   0.0)
  call AllocateDSI( DZCK,     KCM,   0.0)  
  call AllocateDSI( DZCU,     LCM,   KCM,   0.0)      ! *** Not used
  call AllocateDSI( DZCV,     LCM,   KCM,   0.0)      ! *** Not used
  call AllocateDSI( DZG,      LCM,   KCM,   0.0)
  call AllocateDSI( DZGU,     LCM,   KCM,   0.0)
  call AllocateDSI( DZGV,     LCM,   KCM,   0.0)
  call AllocateDSI( LKSZ,     LCM,   KCM,   .false.)
                                          
  call AllocateDSI( CDZDU,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZDV,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZFU,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZFV,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZLU,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZLV,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZMU,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZMV,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZRU,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZRV,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZUU,    LCM,   KCM,   0.0)
  call AllocateDSI( CDZUV,    LCM,   KCM,   0.0)
                                          
  call AllocateDSI( SUB3D,   -LCM,   KCM,   0.0)
  call AllocateDSI( SUB3DO,  -LCM,   KCM,   0.0)
  call AllocateDSI( SVB3D,   -LCM,   KCM,   0.0)
  call AllocateDSI( SVB3DO,  -LCM,   KCM,   0.0)
                                          
  call AllocateDSI( SGZU,     LCM,   KCM,   0.0)
  call AllocateDSI( SGZV,     LCM,   KCM,   0.0)
  call AllocateDSI( FSGZU,    LCM,   KCM,   0.0)
  call AllocateDSI( FSGZV,    LCM,   KCM,   0.0)
  call AllocateDSI( LSGZU,    LCM,   KCM,   .false.)
  call AllocateDSI( LSGZV,    LCM,   KCM,   .false.)
                                          
  call AllocateDSI( HPK,      LCM,   KCM,   0.0)
  call AllocateDSI( HPKI,     LCM,   KCM,   0.0)
  call AllocateDSI( H1PK,     LCM,   KCM,   0.0)
  call AllocateDSI( H2PK,     LCM,   KCM,   0.0)
  call AllocateDSI( UHDY,     LCM,   KCM,   0.0)
  call AllocateDSI( VHDX,     LCM,   KCM,   0.0)
  call AllocateDSI( UHDY1,    LCM,   KCM,   0.0)
  call AllocateDSI( VHDX1,    LCM,   KCM,   0.0)
  call AllocateDSI( UHDY2,    LCM,   KCM,   0.0)
  call AllocateDSI( VHDX2,    LCM,   KCM,   0.0)
  call AllocateDSI( UHDYF,    LCM,   KCM,   0.0)
  call AllocateDSI( VHDXF,    LCM,   KCM,   0.0)
  call AllocateDSI( UHDYF1,   LCM,   KCM,   0.0)
  call AllocateDSI( VHDXF1,   LCM,   KCM,   0.0)
  call AllocateDSI( UHDYF2,   LCM,   KCM,   0.0)
  call AllocateDSI( VHDXF2,   LCM,   KCM,   0.0)

  ! *** WET/DRY BYPASS VARIABLES
  call AllocateDSI( LWET,    -LCM,   0)
  call AllocateDSI( LDRY,    -LCM,   0)

  ! *** ATMOS
  if( NASER >= 1 )then
    call AllocateDSI( IRELH,     NASER,  0)
    call AllocateDSI( MTSALAST,  NASER,  0)
  endif

  ! *** WIND
  if( NWSER >= 1 )then
    call AllocateDSI( MTSWLAST,  NWSER,  0)
    call AllocateDSI( WINDH,     NWSER,  0.0)
  endif

  ! *** ICE
  if( ISICE > 0 )then
    call AllocateDSI( MITLAST,     NISER,   0)           ! *** ISICE = 1  INPUT SERIES LAST INDEX FOR SERIES
    call AllocateDSI( FRAZILICE,   LCM,     KCM,   0.0)  ! *** ISICE = 4 FRAZIL ICE IN EACH CELL AND LAYER (M)
    call AllocateDSI( FRAZILICE1,  LCM,     KCM,   0.0)  ! *** ISICE = 4 FRAZIL ICE IN EACH CELL AND LAYER FOR PREVIOUS TIME STEP (M)
    call AllocateDSI( ICECOVER,    LCM,     0.0)         ! *** ISICE>0 FRACTION OF ICE COVER FOR EACH CELL (M)
    call AllocateDSI( ICERATE,     LCM,     0.0)         ! *** ISICE>2 WATER FLOW RATE TO/FROM WATER COLUMN (M3/S)
    call AllocateDSI( ICETHICK,    LCM,     0.0)         ! *** ISICE>0 ICE THICKNESS FOR EACH CELL (M)
    call AllocateDSI( ICETHICK1,   LCM,     0.0)         ! *** ISICE>0 ICE THICKNESS FOR EACH CELL FOR PREVIOUS TIME STEP (M)
    call AllocateDSI( ICETEMP,     LCM,     0.0)         ! *** ISICE>0 ICE TEMPERATURE OF EACH CELL (M)
    call AllocateDSI( ICEVOL,      LCM,     0.0)         ! *** ISICE>0 INCREMENTAL WATER VOLUMES FOR ICE FORMATION/MELT (M3)
    call AllocateDSI( RICECOVT,    NISER,   0.0)         ! *** ISICE = 1  INPUT ICE COVER AT CURRENT TIME
    call AllocateDSI( RICETHKT,    NISER,   0.0)         ! *** ISICE = 1  INPUT ICE THICKNESS AT CURRENT TIME
  endif
  call AllocateDSI( ICECELL,   LCM,  .false.)            ! *** ISICE>0  PERFORM ICE COMPUTATIONS
  call AllocateDSI( RHOW,      LCM,  KCM,   0.0)         ! *** DENSITY OF WATER KG/M3, FROM CALBUOY

  ! *** Withdrawal/Return
  call AllocateDSI( CQWRSERT,  -NQWRSRM,    NSTVM2,    0.0)
  call AllocateDSI( MTSWRLAST,  NQWRSRM,              0)
  call AllocateDSI( QWRSERT,   -NQWRSRM,   0.0)
  call AllocateDSI( QWRSERTLP, -NQWRSRM,   0.0)
  
  ! *** BLOCKED LAYER FACE OPTION
  if( NBLOCKED > 0 )then
    call AllocateDSI( KBBU,        NBLOCKED,     0)
    call AllocateDSI( KBBV,        NBLOCKED,     0)
    call AllocateDSI( KTBU,        NBLOCKED,     0)
    call AllocateDSI( KTBV,        NBLOCKED,     0)
    call AllocateDSI( LBLOCKED,    NBLOCKED,     0)
    call AllocateDSI( BLANCHORU,   NBLOCKED,   0.0)
    call AllocateDSI( BLANCHORV,   NBLOCKED,   0.0)
    call AllocateDSI( BLDRAFTUO,   NBLOCKED,   0.0)
    call AllocateDSI( BLDRAFTVO,   NBLOCKED,   0.0)
    call AllocateDSI( BLSILLU,     NBLOCKED,   0.0)
    call AllocateDSI( BLSILLV,     NBLOCKED,   0.0)
    call AllocateDSI( BLDRAFTU,    NBLOCKED,   0.0)
    call AllocateDSI( BLDRAFTV,    NBLOCKED,   0.0)
  endif

  call AllocateDSI( CLOE,       NBBEM,      KCM,     NSTVM2, 0.0)
  call AllocateDSI( CLON,       NBBNM,      KCM,     NSTVM2, 0.0)
  call AllocateDSI( CLOS,       NBBSM,      KCM,     NSTVM2, 0.0)
  call AllocateDSI( CLOW,       NBBWM,      KCM,     NSTVM2, 0.0)
  
  call AllocateDSI( TPCOORDE,      NPBEM,  0.0)
  call AllocateDSI( TPCOORDN,      NPBNM,  0.0)
  call AllocateDSI( TPCOORDS,      NPBSM,  0.0)
  call AllocateDSI( TPCOORDW,      NPBWM,  0.0)

  ! *** Temporary use only for debugging and other testing
  call AllocateDSI( DOMAIN_LIST,  LCM_Global,    2,   0)     !< Debugging cell list of each sub-domain
  call AllocateDSI( ITEST,       -LCM,           0)          !< Created and allocated for debug assignments - Not used for routine runs
  call AllocateDSI( TEST1D,      -LCM,         0.0)          !< Created and allocated for debug assignments - Not used for routine runs
  call AllocateDSI( TEST1,        -10,         LCM, 0.0)     !< Created and allocated for debug assignments - Not used for routine runs
  call AllocateDSI( TEST2,        lcm,         kcm, 0.0)     !< Created and allocated for debug assignments - Not used for routine runs
  
  
END SUBROUTINE VARALLOC

