! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE SCANINPMOD
  
  use GLOBAL
  use INFOMOD, only:SKIPCOM,READSTR
  use Allocate_Initialize
  use Variables_WQ
  use HYDSTRUCMOD
  use Variables_Propwash

  use Variables_MPI
  use Variables_MPI_Write_Out
  use Broadcast_Routines
  
  implicit none

  integer :: NSER(8)
  character(*) :: STR*200 ,CFILE*15
  character( * ), PRIVATE, parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  character( * ), PRIVATE, parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  contains

SUBROUTINE SCANEFDC(NCSER1, NCSER2, NCSER3, NCSER4, NCSER5, NCSER6, NCSER7)

  integer, intent(OUT) :: NCSER1, NCSER2, NCSER3, NCSER4, NCSER5, NCSER6, NCSER7

  !INTEGER :: K, NT, I, M, M1, NN
  !INTEGER :: NDUM, ISTOCNT, ISO, ITIDASM, NPFOR, NDATAPTS, ISTYP
  integer :: K, NT, I, M, M1, NN, ITYPE
  integer :: LDUM, NDUM, ISTOCNT, ISO, ITIDASM, NPFOR, NDATAPTS
  integer :: ITMPPMX, ITSSS, NS, IS, IITMP, IJTMP, NPMXZ, NPMXPTS
  integer :: ISBEDTEMI  ! *** DEPRECATED VARIABLE
  real    :: DIFTOXNT, DIFTOXSNT
  real    :: PDIFTOXNT, DPDIFTOXNT, R, TC1, TAV1, TMP
  real    :: R1TMP, R2TMP, R3TMP, R4TMP, R5TMP, PMIXSF

  ! *** ********************************************
  if( process_id == master_id )then
    call OPENFILE('EFDC.INP')

    call SEEK('C1',1)
    read(1,'(A100)') RUNTITLE

    call SEEK('C1A',1)
    read(1,*,err = 10) IS2TIM,IGRIDH,IGRIDV,KMINV,SGZHPDELTA
  endif

  call Broadcast_Scalar(IS2TIM,     master_id)
  call Broadcast_Scalar(IGRIDH,     master_id)
  call Broadcast_Scalar(IGRIDV,     master_id)
  call Broadcast_Scalar(KMINV,      master_id)
  call Broadcast_Scalar(SGZHPDELTA, master_id)

  if( process_id == master_id )then
    call SEEK('C2',1)
    read(1,*,err = 10) ISRESTI, ISRESTO, ISRESTR, ISGREGOR, ISLOG, ISDIVEX, ISNEGH, ISMMC, ldum, ICONTINUE, ISHOW
  endif

  call Broadcast_Scalar(ISRESTI,  master_id)
  call Broadcast_Scalar(ISRESTO,  master_id)
  call Broadcast_Scalar(ISRESTR,  master_id)
  call Broadcast_Scalar(ISGREGOR, master_id)
  call Broadcast_Scalar(ISLOG,    master_id)
  call Broadcast_Scalar(ISDIVEX,  master_id)
  call Broadcast_Scalar(ISNEGH,   master_id)
  call Broadcast_Scalar(ISMMC,    master_id)
  call Broadcast_Scalar(ICONTINUE,master_id)
  call Broadcast_Scalar(ISHOW, master_id)

  if( process_id == master_id )then
    call SEEK('C4',1)
    read(1,*,err = 10) ISLTMT,ISSSMMT,RESSTEP
  endif

  call Broadcast_Scalar(ISLTMT, master_id)
  call Broadcast_Scalar(ISSSMMT,master_id)
  call Broadcast_Scalar(RESSTEP,master_id)

  if( ISLTMT > 0 )then
    call STOPP('ISLTMT LONG TERM SIMULATION DISABLED!')
  endif

  if( process_id == master_id )then
    call SEEK('C5',1)
    read(1,*,err = 10) ISCDMA, ISHDMF, ldum, ISWASP, ISDRY, ICALTB, ISRLID, ISVEG, ISVEGL, ISITB, IHMDSUB, IINTPG
  endif

  call Broadcast_Scalar(ISCDMA,  master_id)
  call Broadcast_Scalar(ISHDMF,  master_id)
  call Broadcast_Scalar(ISWASP,  master_id)
  call Broadcast_Scalar(ISDRY,   master_id)
  call Broadcast_Scalar(ICALTB,  master_id)
  call Broadcast_Scalar(ISRLID,  master_id)
  call Broadcast_Scalar(ISVEG,   master_id)
  call Broadcast_Scalar(ISVEGL,  master_id)
  call Broadcast_Scalar(ISITB,   master_id)
  call Broadcast_Scalar(IHMDSUB, master_id)
  call Broadcast_Scalar(IINTPG,  master_id)

  if( ISCDMA == 10 )then
    call STOPP('EFDC 1D CHANNEL MODE IS DISABLED IN THIS VERSION OF EFDC')
  endif
  if( ISCDMA > 2 )then
    call STOPP('EXPERIMENTAL MOMENTUM OPTIONs ARE ARE DISABLED IN THIS VERSION OF EFDC')
  endif

  if( process_id == master_id )then
    call SEEK('C6',1)
    do NN = 0,8
      read(1,*,err = 10) ISTRAN(NN), ISTOPT(NN), ldum, ISADAC(NN), ISFCT(NN), ldum, ldum, ldum, ISCI(NN), ISCO(NN)
    enddo
  endif

  call Broadcast_Array(ISTRAN, master_id)
  call Broadcast_Array(ISTOPT, master_id)
  call Broadcast_Array(ISADAC, master_id)
  call Broadcast_Array(ISFCT , master_id)
  call Broadcast_Array(ISCI  , master_id)
  call Broadcast_Array(ISCO  , master_id)

  if( process_id == master_id )then
    call SEEK('C7',1)
    read(1,*,err = 10) NTC, NTSPTC, NLTC, NTTC, ldum, NTSTBC, ldum, NTCVB, NTSMMT, ldum, NDRYSTP
  endif

  call Broadcast_Scalar(NTC,     master_id)
  call Broadcast_Scalar(NTSPTC,  master_id)
  call Broadcast_Scalar(NLTC,    master_id)
  call Broadcast_Scalar(NTTC,    master_id)
  call Broadcast_Scalar(NTSTBC,  master_id)
  call Broadcast_Scalar(NTCVB,   master_id)
  call Broadcast_Scalar(NTSMMT,  master_id)
  call Broadcast_Scalar(NDRYSTP, master_id)

  if( process_id == master_id )then
    call SEEK('C9',1)
    read(1,*,err = 10,end = 30) IC, JC, LC, LVC, ldum, NDM, LDM, ISMASK, NBLOCKED, ISCONNECT, ldum, ldum, tmp, tmp
  endif

  call Broadcast_Scalar(IC,       master_id)
  call Broadcast_Scalar(JC,       master_id)
  call Broadcast_Scalar(LC,       master_id)
  call Broadcast_Scalar(LVC,      master_id)
  call Broadcast_Scalar(NDM,      master_id)
  call Broadcast_Scalar(LDM,      master_id)
  call Broadcast_Scalar(ISMASK,   master_id)
  call Broadcast_Scalar(NBLOCKED, master_id)
  call Broadcast_Scalar(ISCONNECT,master_id)

  ! *** TO ENSURE CONSISTENCY OF DECLARATIONS BETWEEN SUB-DOMAINS
  ISBLOCKED = 0
  if( NBLOCKED > 0 ) ISBLOCKED = 1
  
  if( IC >= 3 )then
    ICM = IC + 1
    ICM_Global = ICM
  else
    call STOPP('IC MUST BE AT LEAST 3')
  endif
  if( JC >= 3 )then
    JCM = JC + 1
    JCM_Global = JCM
  else
    call STOPP('IJ MUST BE AT LEAST 3')
  endif
  if( LC >= 3 )then
    LCM = LC + 1
    LCM_Global = LCM !***Save the 'true' global LCM value based on what was read in from the input
    LC_Global = LC
  else
    call STOPP('LC MUST BE AT LEAST 3')
  endif

  if( process_id == master_id )then
    call SEEK('C9A',1)
    read(1,*,IOSTAT = ISO) KC, ldum, ldum, tmp, tmp, ldum
  endif

  call Broadcast_Scalar(KC,      master_id)

  KC = ABS(KC)
  if( KC >= 1 )then
    KCM = KC+1
  else
    call STOPP('KC MUST BE AT LEAST 1')
  endif

  if( process_id == master_id )then
    call SEEK('C12A',1)
    read(1,*,IOSTAT = ISO) ISTOPT(0), ISSQL, ISAVBMX, ISFAVB, ISINWV, ISLLIM, IFPROX, ISVTURB, BC_EDGEFACTOR
  endif

  call Broadcast_Scalar(ISTOPT(0),     master_id)
  call Broadcast_Scalar(ISSQL,         master_id)
  call Broadcast_Scalar(ISAVBMX,       master_id)
  call Broadcast_Scalar(ISFAVB,        master_id)
  call Broadcast_Scalar(ISINWV,        master_id)
  call Broadcast_Scalar(ISLLIM,        master_id)
  call Broadcast_Scalar(IFPROX,        master_id)
  call Broadcast_Scalar(ISVTURB,       master_id)
  call Broadcast_Scalar(BC_EDGEFACTOR, master_id)
  
  if( process_id == master_id )then
    call SEEK('C12B',1)
    read(1,*,IOSTAT = ISO) ISGOTM, IFRICTION, ICALNN, ICALSS, CHARNOCK
  endif
  call Broadcast_Scalar(ISGOTM,        master_id)
  call Broadcast_Scalar(IFRICTION,     master_id)
  call Broadcast_Scalar(ICALNN,        master_id)
  call Broadcast_Scalar(ICALSS,        master_id)
  call Broadcast_Scalar(CHARNOCK,      master_id)
  
  if( process_id == master_id )then
    call SEEK('C14',1)
    read(1,*,err = 10,end = 30) MTIDE, NWSER, NASER, ISGWIT, ISCHAN, ISWAVE, ITIDASM, ISPERC, ISBODYF, ISPNHYDS, ISPROPWASH
  endif

  call Broadcast_Scalar(MTIDE,      master_id)
  call Broadcast_Scalar(NWSER,      master_id)
  call Broadcast_Scalar(NASER,      master_id)
  call Broadcast_Scalar(ISGWIT,     master_id)
  call Broadcast_Scalar(ISCHAN,     master_id)
  call Broadcast_Scalar(ISWAVE,     master_id)
  call Broadcast_Scalar(ITIDASM,    master_id)
  call Broadcast_Scalar(ISPERC,     master_id)
  call Broadcast_Scalar(ISBODYF,    master_id)
  call Broadcast_Scalar(ISPNHYDS,   master_id)
  call Broadcast_Scalar(ISPROPWASH, master_id)

  MTM = max(1,MTIDE)+1
  NWSERM = max(1,NWSER)
  NASERM = max(1,NASER)
  NGWSERM = 1
  NDASER = 1
  NDGWSER = 1

  if( process_id == master_id )then
    call SEEK('C16',1)
    read(1,*,err = 10,end = 30) NPBS, NPBW, NPBE, NPBN, NPFOR, NPFORT, NPSER, PDGINIT
  endif

  call Broadcast_Scalar(NPBS,   master_id)
  call Broadcast_Scalar(NPBW,   master_id)
  call Broadcast_Scalar(NPBE,   master_id)
  call Broadcast_Scalar(NPBN,   master_id)
  call Broadcast_Scalar(NPFOR,  master_id)
  call Broadcast_Scalar(NPFORT, master_id)
  call Broadcast_Scalar(NPSER,  master_id)
  call Broadcast_Scalar(PDGINIT,master_id)

  NPBSM = max(1,NPBS)
  NPBWM = max(1,NPBW)
  NPBEM = max(1,NPBE)
  NPBNM = max(1,NPBN)
  NPSERM = max(1,NPSER)
  NPFORM = max(1,NPFOR,NPSER)
  NDPSER = 1

  ! *** Allocate structures
  allocate(TSPS(NPSERM))

  if( process_id == master_id )then
    call SEEK('C22',1)
    read(1,*,err = 10,end = 30) NDYE, NTOX, NSED, NSND, NCSER1, NCSER2, NCSER3, NCSER4, NCSER5, NCSER6, NCSER7, ldum

    NDYE = max(1,NDYE)
    NDYM = max(1,NDYE)
    NTXM = max(1,NTOX)
    NSCM = max(1,NSED)
    NSNM = max(1,NSND)
    NSEDS = max(1,NSED + NSND)
    
    NCSERM = max(1,NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7)

  endif   ! *** End on master process

  NDCSER = 1

  call Broadcast_Scalar(JC, master_id)
  call Broadcast_Scalar(IC, master_id)
  call Broadcast_Scalar(JCM, master_id)
  call Broadcast_Scalar(ICM, master_id)

  call Broadcast_Scalar(NDYE, master_id)
  call Broadcast_Scalar(NTOX, master_id)
  call Broadcast_Scalar(NSED, master_id)
  call Broadcast_Scalar(NSND, master_id)
  call Broadcast_Scalar(NDYE, master_id)
  call Broadcast_Scalar(NDYM, master_id)
  call Broadcast_Scalar(NTXM, master_id)
  call Broadcast_Scalar(NSCM, master_id)
  call Broadcast_Scalar(NSNM, master_id)
  call Broadcast_Scalar(NSEDS,  master_id)
  call Broadcast_Scalar(NCSERM, master_id)
  call Broadcast_Scalar(NCSER1, master_id)
  call Broadcast_Scalar(NCSER2, master_id)
  call Broadcast_Scalar(NCSER3, master_id)
  call Broadcast_Scalar(NCSER4, master_id)
  call Broadcast_Scalar(NCSER5, master_id)
  call Broadcast_Scalar(NCSER6, master_id)
  call Broadcast_Scalar(NCSER7, master_id)

  ! *** Allocate structures
  allocate(TSSAL(NCSER1))
  allocate(TSTEM(NCSER2))
  allocate(TSDYE(NCSER3,NDYE))
  allocate(TSSFL(NCSER4))
  allocate(TSTOX(NCSER5,NTOX))
  allocate(TSSED(NCSER6,NSED))
  allocate(TSSND(NCSER7,NSND))

  NSER(1) = NCSER1
  NSER(2) = NCSER2
  NSER(3) = NCSER3
  NSER(4) = NCSER4
  NSER(5) = NCSER5
  NSER(6) = NCSER6
  NSER(7) = NCSER7

  allocate(TOXDEP(NTXM))

  if( process_id == master_id )then
    call SEEK('C23',1)
    read(1,*,err = 10,end = 30) NQSIJ, NQJPIJ, NQSER, NQCTL, NQCTLT, NHYDST, NQWR, NQWRSR, ISDIQ, NQCTLSER, NQCTRULES
  endif

  call Broadcast_Scalar(NQSIJ,    master_id)
  call Broadcast_Scalar(NQJPIJ,   master_id)
  call Broadcast_Scalar(NQSER,    master_id)
  call Broadcast_Scalar(NQCTL,    master_id)
  call Broadcast_Scalar(NQCTLT,   master_id)
  call Broadcast_Scalar(NHYDST,   master_id)
  call Broadcast_Scalar(NQWR,     master_id)
  call Broadcast_Scalar(NQWRSR,   master_id)
  call Broadcast_Scalar(ISDIQ,    master_id)
  call Broadcast_Scalar(NQCTLSER, master_id)
  call Broadcast_Scalar(NQCTRULES, master_id)

  NQSIJM = max(1,NQSIJ)
  NJPSM = max(1,NQJPIJ)
  NQSERM = max(1,NQSER)
  NQCTLM = max(1,NQCTL)
  NQCTTM = max(1,NQCTLT)
  NQWRM = max(1,NQWR)
  NQWRSRM = max(1,NQWRSR)
  NDQSER = 1   ! *** Flow              : Maximum number of  points in a series
  NDQWRSR = 1  ! *** Withdrawal/Return : Maximum number of  points in a series

  ! *** SET KB AND CHECK FOR SEDZLJ USAGE
  NSEDFLUME = 0
  NWARNING = 0
  LSEDZLJ = .FALSE.
  if( NSED > 0 .or. NSND > 0 )then
    if( process_id == master_id )then
      call SEEK('C36',1)
      read(1,*,err = 10,end = 30) ISEDINT, ISEDBINT, NSEDFLUME, ISMUD, ISBEDMAP, ISEDVW, ISNDVW, KB, ISDTXBUG
    endif

    call Broadcast_Scalar(ISEDINT,   master_id)
    call Broadcast_Scalar(ISEDBINT,  master_id)
    call Broadcast_Scalar(NSEDFLUME, master_id)
    call Broadcast_Scalar(ISMUD,     master_id)
    call Broadcast_Scalar(ISBEDMAP,  master_id)
    call Broadcast_Scalar(ISEDVW,    master_id)
    call Broadcast_Scalar(ISNDVW,    master_id)
    call Broadcast_Scalar(KB,        master_id)
    call Broadcast_Scalar(ISDTXBUG,  master_id)

    if( KB >= 1 )then
      KBM = KB+1
    else
      call STOPP('KB MUST BE AT LEAST 1')
    endif
    if( ISTRAN(6) > 0 )then
      if( NSEDFLUME == 98 .or. NSEDFLUME == 99 )then
        LSEDZLJ = .TRUE.
        NSEDFLUME = 1
      elseif( NSEDFLUME > 0 )then
        LSEDZLJ = .TRUE.
      endif
    endif

    ! *** LOOK FOR LEGACY APPROACH TO SET SEDZLJ, I.E. IWRSP(1) = 98
    if( .not. LSEDZLJ .and. ISTRAN(6) > 0 )then
      ISDTXBUG = 0
      !C40*  READ COHESIVE SEDIMENT parameter SET 2 REPEAT DATA LINE NSED TIMES
      if( process_id == master_id )then
        call SEEK('C40',1)
        read(1,*,err = 10,end = 30)ISDTXBUG
      endif

      call Broadcast_Scalar(ISDTXBUG, master_id)

      if( ISDTXBUG == 98 )then
        LSEDZLJ = .TRUE.
        NSEDFLUME = 1
      endif
    endif
    
    if( .not. LSEDZLJ .and. ISTRAN(7) > 0 .and. NSND > 0 )then
      if( process_id == master_id )then
        call SEEK('C42A',1)
        do NS = 1,1
          read(1,*,IOSTAT = ISO) ICALC_BL, R, R, R, R, R, R, R, R, R
        enddo
      endif
      
      call Broadcast_Scalar(ICALC_BL , master_id)
    endif
    
    !IF( .not. LSEDZLJ )then
    !  ! *** SET NSCM TO THE MAXIMUM NUMBER OF SEDIMENT CLASSES
    !  NSCM = max(NSCM, NSED+NSND)   delme
    !ENDIF
    
  else
    KBM = 1
  endif
  
  ITMPPMX = 0
  if( NTOX > 0 )then
  
    !C43C*  READ TOXIC TIMING AND VOLATILIZATION SWITCHES
    if(process_id == master_id )then
        call SEEK('C43C',1)
        read(1,*,IOSTAT = ISO) TOXSTEPW, TOXSTEPB
    endif
    
    call Broadcast_Scalar(TOXSTEPW, master_id)
    call Broadcast_Scalar(TOXSTEPB, master_id)
    
    if( process_id == master_id )then
      call SEEK('C44',1)
    endif
    ! *** NEED TO READ EVEN IF SEDZLJ IS BEING USED
    do NT = 1,NTOX
      if( process_id == master_id )then
        read(1,*,err = 10) NDUM, ISTOCNT, DIFTOXNT, DIFTOXSNT, PDIFTOXNT, DPDIFTOXNT
      endif

      call Broadcast_Scalar(ISTOCNT,   master_id)
      call Broadcast_Scalar(DIFTOXNT,  master_id)
      call Broadcast_Scalar(DIFTOXSNT, master_id)
      call Broadcast_Scalar(PDIFTOXNT, master_id)
      call Broadcast_Scalar(DPDIFTOXNT,master_id)

      if( PDIFTOXNT < 0. )ITMPPMX = 1
    enddo

    if( process_id == master_id )then
      call SEEK('C45A',1)
      read(1,*,err = 10) ISTDOCW, ISTPOCW, ISTDOCB, ISTPOCB, STDOCWC, STPOCWC, STDOCBC, STPOCBC
    endif

    call Broadcast_Scalar(ISTDOCW, master_id)
    call Broadcast_Scalar(ISTPOCW, master_id)
    call Broadcast_Scalar(ISTDOCB, master_id)
    call Broadcast_Scalar(ISTPOCB, master_id)
    call Broadcast_Scalar(STDOCWC, master_id)
    call Broadcast_Scalar(STPOCWC, master_id)
    call Broadcast_Scalar(STDOCBC, master_id)
    call Broadcast_Scalar(STPOCBC, master_id)

    if( process_id == master_id )then
      call SEEK('C45E',1)
    endif

    do NT = 1,NTOX
      if( process_id == master_id )then
        read(1,*,err = 10) NDUM, TOXDEP(NT).ITXDRY, TOXDEP(NT).TXDRY, TOXDEP(NT).ITXWET, TOXDEP(NT).TXWET
      endif

      call Broadcast_Scalar(TOXDEP(NT).ITXDRY, master_id)
      call Broadcast_Scalar(TOXDEP(NT).ITXWET, master_id)
    enddo
  endif

  if( process_id == master_id )then
    call SEEK('C46',1)
    read(1,*,ERR = 10,END = 30) BSC, IBSC, TEMO, tmp, tmp, NCBS, NCBW, NCBE, NCBN, IVOLTEMP, VOL_VEL_MAX, VOL_DEP_MIN
  endif

  call Broadcast_Scalar(BSC,      master_id)
  call Broadcast_Scalar(IBSC,     master_id)
  call Broadcast_Scalar(TEMO,     master_id)
  call Broadcast_Scalar(NCBS,     master_id)
  call Broadcast_Scalar(NCBW,     master_id)
  call Broadcast_Scalar(NCBE,     master_id)
  call Broadcast_Scalar(NCBN,     master_id)
  call Broadcast_Scalar(IVOLTEMP, master_id)

  NBBSM = max(1,NCBS)
  NBBWM = max(1,NCBW)
  NBBEM = max(1,NCBE)
  NBBNM = max(1,NCBN)

  if( process_id == master_id )then
    call SEEK('C46A',1)
    read(1,*,err = 10,end = 30) ISICE, NISER, TEMPICE, CDICE, ICETHMX, RICETHK0
  endif

  call Broadcast_Scalar(ISICE,    master_id)
  call Broadcast_Scalar(NISER,    master_id)
  call Broadcast_Scalar(TEMPICE,  master_id)
  call Broadcast_Scalar(CDICE,    master_id)
  call Broadcast_Scalar(ICETHMX,  master_id)
  call Broadcast_Scalar(RICETHK0, master_id)

  if( NASER > 0  )then
    if( process_id == master_id )then
      call SEEK('C46D',1)
      read(1,*,IOSTAT = ISO) IASWRAD, REVC, RCHC, ISVHEAT, SWRATNF, SWRATNS, FSWRATF, TEMTHKO, TEMBO, HTBED1, HTBED2
    endif

    call Broadcast_Scalar(IASWRAD, master_id)
    call Broadcast_Scalar(REVC,    master_id)
    call Broadcast_Scalar(RCHC,    master_id)
    call Broadcast_Scalar(ISVHEAT, master_id)
    call Broadcast_Scalar(SWRATNF, master_id)
    call Broadcast_Scalar(SWRATNS, master_id)
    call Broadcast_Scalar(FSWRATF, master_id)
    call Broadcast_Scalar(TEMTHKO, master_id)
    call Broadcast_Scalar(TEMBO,   master_id)
    call Broadcast_Scalar(HTBED1,  master_id)
    call Broadcast_Scalar(HTBED2,  master_id)

  endif

  if( process_id == master_id )then
    call SEEK('C67',1)
    read(1,*,err = 10) ISPD
  endif

  call Broadcast_Scalar(ISPD,   master_id)

  NPDM = max(1,NPD)

  if( process_id == master_id )then
    call SEEK('C84',1)
    read(1,*,err = 10,end = 30)ISTMSR,MLTMSR,ldum,ldum,NWTMSR,NTSSTSP,TCTMSR
  endif

  call Broadcast_Scalar(ISTMSR, master_id)
  call Broadcast_Scalar(MLTMSR, master_id)
  call Broadcast_Scalar(NWTMSR, master_id)
  call Broadcast_Scalar(NTSSTSP,master_id)
  call Broadcast_Scalar(TCTMSR, master_id)

  MLTMSRM = max(1,MLTMSR)
  NTSSTSPM = max(1,NTSSTSP)
  MTSSTSPM = 1
  if( process_id == master_id )then
    if( NTSSTSP > 0 )then
      call SEEK('C85',1)
      do ITSSS = 1,NTSSTSP
        read(1,*,err = 10,end = 30)I,M
        MTSSTSPM = max(MTSSTSPM,M)
      enddo
    endif

  endif

  call Broadcast_Scalar(MTSSTSPM, master_id)

  ! ***    
  if( process_id == master_id )then
    call SEEK('C91',1)
    read(1,*,ERR = 10,END = 30) NCDFOUT, DEFLEV, ROTA, BLK, UTMZ, HREST, BASEDATE, BASETIME, PROJ
    ! *** close efdc.inp
    close(1)
  endif
  
  call Broadcast_Scalar(NCDFOUT  , master_id)
  call Broadcast_Scalar(DEFLEV   , master_id)
  call Broadcast_Scalar(ROTA     , master_id)
  call Broadcast_Scalar(BLK      , master_id)
  call Broadcast_Scalar(UTMZ     , master_id)
  call Broadcast_Scalar(HREST    , master_id)
  call Broadcast_Scalar(BASEDATE , master_id)
  call Broadcast_Scalar(BASETIME , master_id)
  call Broadcast_Scalar(PROJ     , master_id)
  
  ! ***********************************************************************************************
  ! *** END OF EFDC.INP SCANNING
  ! ***********************************************************************************************

  if( ISVEG >= 1 )then

    if( process_id == master_id )then
      call OPENFILE('VEGE.INP')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*,err = 10,end = 30) MVEGTYP, MVEGOW, NVEGSER, UVEGSCL
    endif

    call Broadcast_Scalar(MVEGTYP,master_id)
    call Broadcast_Scalar(MVEGOW, master_id)
    call Broadcast_Scalar(NVEGSER,master_id)
    call Broadcast_Scalar(UVEGSCL,master_id)

    NVEGTPM = max(NVEGTPM,MVEGTYP)
    NVEGSERM = max(NVEGSERM,NVEGSER)

    if( process_id == master_id )then
      close(1)
    endif

    if( NVEGSER >= 1 )then
      if( process_id == master_id )then
        call OPENFILE('VEGSER.INP')
        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      endif

      do NS = 1,NVEGSER
        if( process_id == master_id )then
          read(1,*,err = 10,end = 30) M1,TC1,TAV1
        endif

        call Broadcast_Scalar(M1, master_id)

        NDVEGSER = max(NDVEGSER,M1)

        if( process_id == master_id )then
          do M = 1,M1
            read(1,*)
          enddo
        endif
      enddo
      if( process_id == master_id )then
        close(1)
      endif

    endif

    ! *** DETERMINE IF MHK IS USED
    LMHK = .FALSE.
    if( process_id == master_id )then
      call OPENFILE('DXDY.INP')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
    endif

    write(mpi_log_unit,'(a,i6)') 'LVC: ', LVC

    call Broadcast_Scalar(LVC, master_id)

    TCOUNT = 0
    do I = 1,LVC
      if( process_id == master_id )then
        read(1,*)IITMP,IITMP,R1TMP,R2TMP,R3TMP,R4TMP,R5TMP,IJTMP
      endif

      call Broadcast_Scalar(IJTMP, master_id)

      if( IJTMP > 90 )then
        LMHK = .TRUE.
        TCOUNT = TCOUNT+1
      endif
    enddo

    if( process_id == master_id )then
      close(1)
    endif

    if( LMHK )then
      if( process_id == master_id )then
        call OPENFILE('MHK.INP')
        call SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES
        read(1,*,err = 10,end = 30)MHKTYP
        close(1)
      endif

      call Broadcast_Scalar(MHKTYP, master_id)

    endif
  endif

  ! *** BANK EROSION
  if( ISBKERO >= 1 )then
    if( process_id == master_id )then
      call OPENFILE('BEMAP.INP')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*,err = 10,end = 30)NBEPAIR,NBESER
    endif

    call Broadcast_Scalar(NBEPAIR, master_id)
    call Broadcast_Scalar(NBESER,  master_id)

    NBEPAIRM = NBEPAIR
    NBESERM = NBESER

    if( process_id == master_id )then
      close(1)
    endif

    NDBESER = 1
    if( NBESER > 0 )then
      if( process_id == master_id )then
        call OPENFILE('BESER.INP')
      endif

      do NS = 1,NBESER
        if( process_id == master_id )then
100       read(1,*,err = 100,end = 30)M,R,R,R,R
        endif

        call Broadcast_Scalar(M, master_id)

        NDBESER = max(NDBESER,M)
        do I = 1,M
          if( process_id == master_id )then
            read(1,*,err = 10,end = 30)R,R,R
          endif
        enddo
      enddo

      if( process_id == master_id )then
        close(1)
      endif

    endif
  endif

  ! *** TOXIC SIMULATION FILES
  NPMXPTSM = 1
  NPMXZM = 1
  if( ISTRAN(5) > 0 .and. NTOX > 0 )then
    ! *** SEDIMENT MIXING
    if( ITMPPMX == 1 )then
      call OPENFILE('PARTMIX.INP')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*)NPMXZ,NPMXPTS,PMIXSF
      NPMXPTSM = NPMXPTS
      NPMXZM = NPMXZ
      close(1)
    endif

    ! *** TIME SERIES DRY DEPOSITION BEING USED
    if( SUM(TOXDEP(:).ITXDRY) > 0 )then
      call OPENFILE('TXDRY.INP')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*) ITYPE, NDATAPTS
      TXDRYSER(1).NREC = NDATAPTS

      allocate(TXDRYSER(1).TIM(NDATAPTS), TXDRYSER(1).VAL(NDATAPTS,NTXM))
      close(1)
    endif

    ! *** TIME SERIES WET DEPOSITION BEING USED
    if( SUM(TOXDEP(:).ITXWET) > 0 )then
      call OPENFILE('TXWET.INP')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*) ITYPE, NDATAPTS
      TXWETSER(1).NREC = NDATAPTS

      allocate(TXWETSER(1).TIM(NDATAPTS), TXWETSER(1).VAL(NDATAPTS,NTXM))
      close(1)
    endif
  endif

  if( ISCONNECT >= 2 )then
    call OPENFILE('MAPPGEW.INP')

    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR = READSTR(1)
    read(1,*,IOSTAT = ISO) NPEWBP
    close(1)
  endif

  if( ISCONNECT == 1 .or. ISCONNECT == 3 )then
    call OPENFILE('MAPPGNS.INP')

    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR = READSTR(1)
    read(1,*,IOSTAT = ISO) NPNSBP ! global value for now
    close(1)
  endif

  MDCHH = 0
  if( ISCHAN > 0 )then
    if( process_id == master_id )then
      call OPENFILE('MODCHAN.INP')

      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      read(1,*) MDCHH,MDCHHD,MDCHHD2
      close(1)
    endif
    call Broadcast_Scalar(MDCHH, master_id)
  endif

  !end if !***End on master process

  if( NWSER > 1 )then
    ! *** ****************************
    if( process_id == master_id )then
      call OPENFILE('WNDMAP.INP')
      STR = READSTR(1)
      read(1,*,IOSTAT = ISO) NWNDMAP
      if( ISO /= 0) CALL STOPP(' *** WNDMAP.INP: READING ERROR!')
      close(1)
    endif   ! *** End on master process

    call Broadcast_Scalar(NWSER, master_id)
    call Broadcast_Scalar(NWNDMAP, master_id)

    allocate(TWNDMAPBEG(NWNDMAP),TWNDMAPEND(NWNDMAP))
    allocate(WNDWHT(NWSER,LCM,NWNDMAP))
    allocate(WNDWHT_TEMP(NWSER))

  endif

  call Broadcast_Scalar(NASER, master_id)

  if( NASER > 1 )then
    ! *** ****************************
    if( process_id == master_id )then
      call OPENFILE('ATMMAP.INP')
      STR = READSTR(1)
      read(1,*,IOSTAT = ISO) NATMMAP
      if( ISO /= 0) CALL STOPP(' *** ATMMAP.INP: READING ERROR!')
      close(1)
    endif   ! *** End on master process

    call Broadcast_Scalar(NATMMAP, master_id)

    allocate(TATMMAPBEG(NATMMAP),TATMMAPEND(NATMMAP))
    allocate(ATMWHT(NASER,LCM,NATMMAP))
    allocate(ATMWHT_TEMP(NASER))     ! Temporary storage for remapping in input.f90 because of domain decomposition
  endif

  if( ISICE == 1 .and. NISER > 1 )then
    ! *** ****************************
    if( process_id == master_id )then
      call OPENFILE('ICEMAP.INP')

      STR = READSTR(1)
      read(1,*,IOSTAT = ISO) NICEMAP
      if( ISO /= 0) CALL STOPP(' *** ICEMAP.INP: READING ERROR!')
      close(1)
    endif   ! *** End on master process

    call Broadcast_Scalar(NICEMAP, master_id)
    call Broadcast_Scalar(NISER,   master_id)

    allocate(TICEMAPBEG(NICEMAP),TICEMAPEND(NICEMAP))
    allocate(RICEWHT(NICEMAP,LCM,NISER))
    allocate(RICEWHT_Global(NICEMAP,LCM_Global,NISER))
  endif

  return
  ! *** ERROR MESSAGES
10 WRITE(6,'(A)') '  READ ERROR '//CARDNO//' IN INPUT FILE: '//TRIM(CFILE)
  write(mpi_error_unit,'(A)') '  READ ERROR '//CARDNO//' IN INPUT FILE: '//TRIM(CFILE)
  call STOPP('.')

30 WRITE(6,'(A)') '  READ ERROR '//CARDNO//': UNEXPECTED END OF INPUT FILE: '//TRIM(CFILE)
  write(mpi_error_unit,'(A)') '  READ ERROR '//CARDNO//': UNEXPECTED END OF INPUT FILE: '//TRIM(CFILE)
  call STOPP('.')

END SUBROUTINE

SUBROUTINE SCANMODC
  integer :: M,I
  real    :: R

  call OPENFILE('MODCHAN.INP')

10 READ(1,*,err = 10,end = 40)M,I,I
  NCHANM = max(1,M)
  read(1,*,err = 20,end = 40)I,I,R
  close(1)
  return
20 WRITE(6,'(A)') '  READ ERROR IN INPUT FILE: MODCHAN.INP'
  write(mpi_error_unit,'(A)') '  READ ERROR IN INPUT FILE: MODCHAN.INP'
  call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: MODCHAN.INP'
  write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: MODCHAN.INP'
  call STOPP('.')

END SUBROUTINE

SUBROUTINE SCANGWSR
  integer :: I,J,M,NS
  real    :: R,T,F

  if( process_id == master_id )then
    call OPENFILE('GWSER.INP')

  10 READ(1,*,err = 10,end = 40) NGWSER
    NGWSERM = max(1,NGWSER)
    do NS = 1,NGWSER
      read(1,*,err = 20,end = 40)M,R,R,R,R,I
      NDGWSER = max(NDGWSER,M)
      do I = 1,M
        read(1,*,err = 20,end = 40)T,F,(R,J = 1,3+NDYE+NSED+NSND+NTOX)   ! DELME - WQ
      enddo
    enddo
    close(1)
  endif
  
  call Broadcast_Scalar(NGWSERM, master_id)
  call Broadcast_Scalar(NGWSER, master_id)
  call Broadcast_Scalar(NDGWSER, master_id)
  
  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: GWSER.INP IN SERIES:',NS,', POINT:',I
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: GWSER.INP IN SERIES:',NS,', POINT:',I
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: GWSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: GWSER.INP'
   call STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANASER
  use INFOMOD,ONLY:SKIPCOM
  integer :: M,I,NA,NS
  real    :: R
  character*120 LIN,STR*200

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('ASER.INP')
  endif   ! *** End on master process

  call Broadcast_Scalar(NASER, master_id)

  ! *** Allocate structures
  allocate(TSATM(NASER))

  ! *** ****************************
  if( process_id == master_id )then
    call SKIPCOM(1,'*')
    read(1,'(A)') STR
    STR = READSTR(1)
  endif

  do NA = 1,NASER
    ! *** ****************************
    if( process_id == master_id )then
      read(1,*,end = 40)M,R,R,I,R,R,R,R
    endif   ! *** End on master process

    call Broadcast_Scalar(M, master_id)

    NDASER = max(NDASER,M)
    allocate(TSATM(NA).TIM(M))
    allocate(TSATM(NA).VAL(M,7))
    
    TSATM(NA).NREC= M
    TSATM(NA).TIM = 0.0
    TSATM(NA).VAL = 0.0

    ! *** ****************************
    if( process_id == master_id )then
      do I = 1,M
        read(1,*,err = 20,end = 40)R,R,R,R,R,R,R,R
      enddo
    endif   ! *** End on master process

  enddo
  
  ! *** ****************************
  if( process_id == master_id )then
    close(1)
  endif   ! *** End on master process

  ! *** ****************************
  if( process_id == master_id )then
    if( ISTRAN(8) > 0 )then
      if( IWQSUN == 1 )then
        call OPENFILE('SUNDAY.INP')

        M = 0
        call SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES
        !DO I = 1,7
        !  read(1,'(A120)',err = 30,end = 40)LIN
        !END DO
        read(1,*,err = 30,end = 40)M,R,R,R,R
        close(1)
        NDASER = max(NDASER,M)
      endif
    endif
  endif   ! *** End on master process

  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: ASER.INP IN SERIES:',NS,', POINT:',I
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: ASER.INP IN SERIES:',NS,', POINT:',I
   call STOPP('.')
  
30 WRITE(6,'(A)') '  READ ERROR IN INPUT FILE'
   write(mpi_error_unit,'(A)') '  READ ERROR IN INPUT FILE'
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE'
   call STOPP('.')

END SUBROUTINE

SUBROUTINE SCANSSER
  integer :: NS,I,J,M,K
  real   :: R

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('SSER.INP')
  endif   ! *** End on master process

  do NS = 1,NSER(1)
    ! *** ****************************
    if( process_id == master_id )then
10    read(1,*,err = 10,end = 40)I,M,R,R,R,R
    endif   ! *** End on master process

    call Broadcast_Scalar(M, master_id)

    NDCSER = max(NDCSER,M)

    allocate(TSSAL(NS).VAL(M,KCM),TSSAL(NS).TIM(M))
    TSSAL(NS).NREC = M
    TSSAL(NS).VAL(:,:) = 0
    TSSAL(NS).TIM(:) = 0

    ! *** ****************************
    if( process_id == master_id )then
      if( I == 1 )then
        read(1,*,err = 20,end = 40)(R,K = 1,KC)
        do J = 1,M
          read(1,*,err = 20,end = 40)R,R
        enddo
      else
        do J = 1,M
          read(1,*,err = 20,end = 40)R,(R,K = 1,KC)
        enddo
      endif
    endif   ! *** End on master process

  enddo
  ! *** ****************************
  if( process_id == master_id )then
    close(1)
  endif   ! *** End on master process

  return  
   ! *** ****************************
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SSER.INP IN SERIES:',NS,', POINT:',J
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SSER.INP IN SERIES:',NS,', POINT:',J
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: SSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: SSER.INP'
   call STOPP('.')

END SUBROUTINE

SUBROUTINE SCANTSER
  integer   :: NS,I,J,K,M,NN
  real      :: R
  character(200) :: STR

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('TSER.INP')
    STR = READSTR(1)
  endif

  call Broadcast_Scalar(NSER(2), master_id)

  do NS = 1,NSER(2)
    ! *** ****************************
    if( process_id == master_id )then
      read(1,*,err = 991) I, M, R, R, R, R
    endif

    call Broadcast_Scalar(M, master_id)

    NDCSER = max(NDCSER,M)

    allocate(TSTEM(NS).VAL(M,KCM),TSTEM(NS).TIM(M))
    TSTEM(NS).NREC = M
    TSTEM(NS).VAL(:,:) = 0
    TSTEM(NS).TIM(:) = 0
    ! *** ****************************
    if( process_id == master_id )then
      if( I == 1 )then
        read(1,*,err = 991) (R,K = 1,KC)
        do J = 1,M
          read(1,*,err = 991) R, R
        enddo
      else
        do J = 1,M
          read(1,*,err = 991) R, (R,K = 1,KC)
        enddo
      endif
    endif   ! *** End on master process

  enddo
  
  ! *** ****************************
  if( process_id == master_id )then
    close(1)
  endif 

  if( ISICE == 1 .and. NISER >= 1 )then

    ! *** ****************************
    if( process_id == master_id )then
      call OPENFILE('ISER.INP')
      STR = READSTR(1)
    endif   ! *** End on master process

    allocate(TSICE(NISER))

    do NS = 1,NISER
      ! *** ****************************
      if( process_id == master_id )then
        read(1,*,err = 992) M,R,R,R
      endif
      
      call Broadcast_Scalar(M, master_id)

      allocate(TSICE(NS).VAL(M,2),TSICE(NS).TIM(M))  ! *** VAL(M,1): ICE THICKNESS; VAL(M,2): ICE COVER
      TSICE(NS).NREC= M
      TSICE(NS).VAL = 0
      TSICE(NS).TIM = 0
      ! *** ****************************
      if( process_id == master_id )then
        do J = 1,M
          read(1,*,err = 992)
        enddo
      endif  ! *** End on master process

    enddo

    ! *** ****************************
    if( process_id == master_id )then
      close(1)
    endif   ! *** End on master process

  elseif( ISICE == 2 )then
    ! *** ****************************
    if( process_id == master_id )then
      call OPENFILE('ISTAT.INP')
      STR = READSTR(1)
      read(1,*,err = 993) M,R,R   !MISER(NN),TCISER(NN),TAISER(NN)
    endif   ! *** End on master process

    call Broadcast_Scalar(M, master_id)
    NN = 1

    allocate(TSICE(NN))
    allocate(TSICE(NN).VAL(M,2),TSICE(NN).TIM(M))  ! *** VAL(M,1): ICE THICKNESS; VAL(M,2): ICE COVER
    TSICE(NN).NREC= M
    TSICE(NN).VAL = 0
    TSICE(NN).TIM = 0

    ! *** ****************************
    if( process_id == master_id )then
      close(1)
    endif   ! *** End on master process

  endif

  return
  
991 call STOPP('TSER.INP: READING ERROR')
992 call STOPP('ISER.INP: READING ERROR')
993 call STOPP('ISTAT.INP: READING ERROR')
END SUBROUTINE

SUBROUTINE SCANDSER
  integer :: NS,I,K,M,MD,ISTYP
  real   :: R
  character*120 SKIP
  character(*) :: STR*200

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('DSER.INP')
    ! *** SKIP HEADER
    STR = READSTR(1)
  endif
  

  do NS = 1,NSER(3)
    ! *** ****************************
    if( process_id == master_id )then
      read(1,*,err = 20,end = 40) ISTYP, M
    endif
    
    call Broadcast_Scalar(M, master_id)

    NDCSER = max(NDCSER,M)

    ! *** Allocate structure
    do MD = 1,NDYE
      allocate(TSDYE(NS,MD).VAL(M,KCM),TSDYE(NS,MD).TIM(M))
      TSDYE(NS,MD).NREC = M
      TSDYE(NS,MD).VAL(:,:) = 0.
      TSDYE(NS,MD).TIM(:) = 0.
    enddo

    if( process_id == master_id )then
      if( ISTYP == 1 )then
        read(1,*,err = 20,end = 40) (R, K = 1,KC)
      endif

      ! *** SKIP TO THE NEXT SERIES
      do I = 1,M
        do MD = 1,NDYE
          read(1,'(A)') SKIP
        enddo
      enddo
    endif
  enddo
  close(1)
  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: DSER.INP IN SERIES:',NS,', POINT:',I
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: DSER.INP IN SERIES:',NS,', POINT:',I
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: DSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: DSER.INP'
   call STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANSFSR
  integer :: NS,I,J,K,M
  real   :: R

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('SFSER.INP')
  endif

  do NS = 1,NSER(4)
    ! *** ****************************
    if( process_id == master_id )then
10    read(1,*,err = 10,end = 40)I,M,R,R,R,R
    endif
    call Broadcast_Scalar(M, master_id)

    NDCSER = max(NDCSER,M)

    allocate(TSSFL(NS).VAL(M,KCM),TSSFL(NS).TIM(M))
    TSSFL(NS).NREC = M
    TSSFL(NS).VAL(:,:) = 0
    TSSFL(NS).TIM(:) = 0
    ! *** ****************************
    if( process_id == master_id )then
      if( I == 1 )then
        read(1,*,err = 20,end = 40)(R,K = 1,KC)
        do J = 1,M
          read(1,*,err = 20,end = 40)R,R
        enddo
      else
        do J = 1,M
          read(1,*,err = 20,end = 40)R,(R,K = 1,KC)
        enddo
      endif
    endif   ! *** End on master process

  enddo
  ! *** ****************************
  if( process_id == master_id )then
    close(1)
  endif   ! *** End on master process

  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SFSER.INP IN SERIES:',NS,', POINT:',J
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SFSER.INP IN SERIES:',NS,', POINT:',J
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: SFSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: SFSER.INP'
   call STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANQSER
  integer :: NS,I,J,M,K
  real    :: R,T,Q

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('QSER.INP')
  endif

  allocate(TSFL(NQSER))

  do NS = 1,NQSER
    ! *** ****************************
    if( process_id == master_id )then
10    read(1,*,err = 10,end = 40)I,M,R,R,R,R,J
    endif

    call Broadcast_Scalar(M, master_id)

    ! *** Allocate structure
    allocate( TSFL(NS).VAL(M,KCM) )
    allocate( TSFL(NS).TIM(M) )
    TSFL(NS).NREC = M
    TSFL(NS).VAL(:,:) = 0.0
    TSFL(NS).TIM(:) = 0.0
    ! *******************************
    
    if( process_id == master_id )then
      if( I == 1 )then
        read(1,*,err = 20,end = 40)(R,K = 1,KC)
        do J = 1,M
          read(1,*,err = 20,end = 40)T,Q
        enddo
      else
        do J = 1,M
          read(1,*,err = 20,end = 40)T,(Q,K = 1,KC)
        enddo
      endif
    endif   ! *** End on master process

  enddo
  ! *** ****************************
  if( process_id == master_id )then
    close(1)
  endif   ! *** End on master process

  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QSER.INP IN SERIES:',NS,', POINT:',J
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QSER.INP IN SERIES:',NS,', POINT:',J
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: QSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: QSER.INP'
   call STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANWRSER
  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  integer :: NTMP,I,J,M,NV,NS
  real(RKD)   :: R,T,Q
  type(fson_value), Pointer :: json_data
  
  ! *** ****************************
    ! *** Allocate time series of Withdrawal/Return structure type
  allocate(TSWR(NQWRSR))
  
  if( process_id == master_id )then

    ! *** SETTING NUMBER OF CONSTITUENTS (NTMP)
    NTMP = 3 + NDYM + NSED + NSND + NTOX
    
    ! *** Handle Water Quality variables, if needed
    if( ISTRAN(8) > 0 )then
      json_data => fson_parse("wq_3dwc.jnp")
    
      call fson_get(json_data, "number_of_algae_groups",       NALGAE)
      call fson_get(json_data, "number_of_zooplankton_groups", NZOOPL)
      NTMP = NTMP + 19 + NALGAE + NZOOPL
    endif

    call OPENFILE('QWRS.INP')
    
    do NS = 1,NQWRSR
10    read(1,*,err = 10,end = 40) I ,M, R, R, R, R
      NDQWRSR = max(NDQWRSR,M,1)
    
      if( I == 0 )then
        ! *** Flow Only
        do J = 1,M
          read(1,*,err = 20,end = 40)T,Q
        enddo
      else
        ! *** Flow with Rise/Fall
        do J = 1,M
          read(1,*,err = 20,end = 40)T,Q,(R,NV = 1,NTMP)
        enddo
      endif
    enddo

    close(1)
  endif !*** Writing on master

  call Broadcast_Scalar(NTMP, master_id)
  call Broadcast_Scalar(NDQWRSR, master_id)
        
  ! *** Allocate structure
  do NS = 1,NQWRSR
    allocate(TSWR(NS).TIM(NDQWRSR))
    allocate(TSWR(NS).VAL(NDQWRSR,0:NTMP))
    TSWR(NS).NREC = NDQWRSR
    TSWR(NS).VAL(:,:) = 0.0
    TSWR(NS).TIM(:) = 0.0
  enddo
  
  return

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QWRS.INP IN SERIES:',NS,', POINT:',J
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QWRS.INP IN SERIES:',NS,', POINT:',J

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: QWRS.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: QWRS.INP'
   call STOPP('.')

END SUBROUTINE

SUBROUTINE SCANPSER
  integer :: NS,M,I,I1
  real   :: R,T,E

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('PSER.INP')
  endif

  do NS = 1,NPSER
    if( process_id == master_id )then
  10  read(1,*,err = 10,end = 40) I1, M, R, R, R, R
      NDPSER = max(NDPSER,M)
    endif

    call Broadcast_Scalar(M, master_id)

    ! *** Allocate structure
    allocate( TSPS(NS).VAL(M,2) )
    allocate( TSPS(NS).TIM(M) )
    TSPS(NS).NREC = M
    TSPS(NS).VAL(:,:) = 0.0
    TSPS(NS).TIM(:) = 0.0
    ! *******************************

    if( process_id == master_id )then
      if( I1 == 1)READ(1,*) R, R, R
      do I = 1,M
        if( I1 == 0 )then
          read(1,*,err = 20,end = 40) T, E
        elseif( I1 == 1 )then
          read(1,*,err = 20,end = 40) T, E, R
        endif
      enddo
    endif   ! *** End on master process

  enddo

  if( process_id == master_id )then
    close(1)
  endif   ! *** End on master process

  call Broadcast_Scalar(NDPSER, master_id)

  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: PSER.INP IN SERIES:',NS,', POINT:',I
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: PSER.INP IN SERIES:',NS,', POINT:',I
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: PSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: PSER.INP'
   call STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANWSER
  integer :: I,M,NS
  real    :: R

  ! *** ****************************

  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('WSER.INP')
  endif   ! *** End on master process

  allocate(TSWND(NWSER))

  do NS = 1,NWSER
    ! *** ****************************
    if( process_id == master_id )then
10    read(1,*,err = 10,end = 40)M,R,R,R,I
    endif

    call Broadcast_Scalar(M, master_id)

    allocate(TSWND(NS).TIM(M))
    allocate(TSWND(NS).VAL(M,2))
    TSWND(NS).NREC= M
    TSWND(NS).VAL = 0.0
    TSWND(NS).TIM = 0.0

    ! *** ****************************
    if( process_id == master_id )then
      do I = 1,M
        read(1,*,err = 20,end = 40)R,R,R
      enddo
    endif   ! *** End on master process

  enddo

  ! *** ****************************
  if( process_id == master_id )then
    close(1)
  endif   ! *** End on master process

  return
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: WSER.INP IN SERIES:',NS,', POINT:',I
   write(mpi_error_unit,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: WSER.INP IN SERIES:',NS,', POINT:',I
   call STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: WSER.INP'
   write(mpi_error_unit,'(A)') '  UNEXPECTED END OF INPUT FILE: WSER.INP'
   call STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANQCTL
  character*80 STR*200
  character*80 :: SKIP

  integer :: M, M1, M2, MP, IS, NS, ISO, ISTYP
  real    :: HUA, HUM, HDA, HDM, R, A, A1

  if( process_id == master_id )then
  
    call OPENFILE('QCTL.INP')

    ! *** FIND THE MAXIMUM NUMBER OF TABLE DATA POINTS
    NDQCLT = 0
    STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
    do NS = 1,NQCTLT
      read(1,*,IOSTAT = ISO) ISTYP, MP, HUA, HUM, HDA, HDM, R, A, A1
      NDQCLT = max(NDQCLT,MP)
      if( ISO > 0 )GOTO 20

      if( ISTYP == 0 )then
        do M = 1,MP
          read(1,'(A)')SKIP
        enddo
      elseif( ISTYP == 1 )then
        read(1,'(A)')SKIP
        do M = 1,MP
          read(1,'(A)')SKIP
        enddo
      elseif( ISTYP == 2 )then
        do M1 = 1,MP
          do M2 = 1,MP
            read(1,'(A)')SKIP
          enddo
        enddo
      elseif( ISTYP == 3 )then
        read(1,'(A)')SKIP
        do M1 = 1,MP
          do M2 = 1,MP
            read(1,'(A)')SKIP
          enddo
        enddo
      endif
    enddo

    close(1)
  endif
  call Broadcast_Scalar(NDQCLT, master_id)
  
  return
20 WRITE(6,'(A)') ' READ ERROR IN FILE: '//TRIM(CFILE)
   write(mpi_error_unit,'(A)') ' READ ERROR IN FILE: '//TRIM(CFILE)
   call STOPP('.')

END SUBROUTINE

SUBROUTINE SCANWQ

  use SHELLFISHMOD
  use WQ_RPEM_MODULE
  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  use INFOMOD, only:SKIPCOM 
  
  Character(len = 2)   :: SNUM 
  Character(len = 15)  :: FNWQSR(40)
  Character(len = 15)  :: FNWQSRALGAE(40)
  Character(len = 100) :: ERRSTR
  Character(len = 120) :: LINE
  
  integer :: NW, NPS, I, J, K, IDRAG, ITMP, M, NS, ISTYP, NAL, NDWQCSR
  real    :: XPSQ, TM, TA, RMULADJ, ADDADJ, T1, T2, TSMTSB, TSMTSE, SMTSDT
  logical :: BFLAG

  type(fson_value), Pointer :: json_data, item, pointsource, algaegroups

  ! *** ****************************
  if( process_id == master_id )then 
    write(*,'(A)') 'SCANNING EUTROPHICATION CONTROL FILE'
    
    json_data => fson_parse("wq_3dwc.jnp")
    
    call fson_get(json_data, "kinetics_option",                                  ISWQLVL)
    call fson_get(json_data, "number_of_kinetic_zones",                          NWQZ)
    call fson_get(json_data, "temperature_lookup_table_size",                    NWQTD)
    call fson_get(json_data, "number_of_time_series_output_locations",           NWQTS)
    call fson_get(json_data, "number_of_time_series_output_variables",           NTSWQV)
    call fson_get(json_data, "number_of_sediment_zones",                         NSMZ)
    call fson_get(json_data, "max_number_of_time_series_output_locations",       NSMTS)

    call fson_get(json_data, "point_source_load_option",                         IWQPSL)
    call fson_get(json_data, "mass_loading_point_sources.number_of_point_sources",NWQPS) 
    call fson_get(json_data, "mass_loading_point_sources.number_of_time_series", NPSTMSR) 
    
    call fson_get(json_data, "number_of_algae_groups",                           NALGAE)
                                                                                 
    call fson_get(json_data, "zooplankton_activate",                             IWQZPL)
    call fson_get(json_data, "number_of_zooplankton_groups",                     NZOOPL)
                                                                                 
    call fson_get(json_data, "shellfish_farm_activate",                          ISFFARM)
    call fson_get(json_data, "number_of_shellfish_species",                      NSF)
    
    ! *** Update number of WQ components
    NWQV = 19 + NALGAE + NZOOPL
    
    ! *** Set constituent index based on kinetic option
    if( ISWQLVL == 0 )then
      ICHC = 1
      ICHD = 2
      ICHG = 3
      IROC = 4
      ILOC = 5
      IDOC = 6
      IROP = 7
      ILOP = 8
      IDOP = 9
      IP4D = 10
      IRON = 11
      ILON = 12
      IDON = 13
      INHX = 14
      INOX = 15
      ISUU = 16
      ISAA = 17
      ICOD = 18
      IDOX = 19
      ITAM = 20
      IFCB = 21
      ICO2 = 22
    else
      IROC = 1
      ILOC = 2
      IDOC = 3
      IROP = 4
      ILOP = 5
      IDOP = 6
      IP4D = 7
      IRON = 8
      ILON = 9
      IDON = 10
      INHX = 11
      INOX = 12
      ISUU = 13
      ISAA = 14
      ICOD = 15
      IDOX = 16
      ITAM = 17
      IFCB = 18
      ICO2 = 19
    endif
    
    ! *** Nutrient components name
    WQCONSTIT(IROC) = 'ROC'
    WQCONSTIT(ILOC) = 'LOC'
    WQCONSTIT(IDOC) = 'DOC'
    WQCONSTIT(IROP) = 'ROP'
    WQCONSTIT(ILOP) = 'LOP'
    WQCONSTIT(IDOP) = 'DOP'
    WQCONSTIT(IP4D) = 'P4D'
    WQCONSTIT(IRON) = 'RON'
    WQCONSTIT(ILON) = 'LON'
    WQCONSTIT(IDON) = 'DON'
    WQCONSTIT(INHX) = 'NHX'
    WQCONSTIT(INOX) = 'NOX'
    WQCONSTIT(ISUU) = 'SUU'
    WQCONSTIT(ISAA) = 'SAA'
    WQCONSTIT(ICOD) = 'COD'
    WQCONSTIT(IDOX) = 'DOX'
    WQCONSTIT(ITAM) = 'TAM'
    WQCONSTIT(IFCB) = 'FCB'
    WQCONSTIT(ICO2) = 'CO2' 
    
    if( ISWQLVL == 0 )then
      ! *** Algae names if WQSKE0
      WQCONSTIT(ICHC) = 'CHC'
      WQCONSTIT(ICHD) = 'CHD'
      WQCONSTIT(ICHG) = 'CHG'
    else
    ! *** Algae names if WQSKE1
      do NW = 1,NALGAE
        if( NW < 10 )then
          write(SNUM,'(I1)') NW
          WQCONSTIT(19+NW) = 'ALG'//SNUM
        else
          write(SNUM,'(I2)') NW
          WQCONSTIT(19+NW) = 'ALG'//SNUM  
        endif
      enddo
    endif
    ! *** Zooplankton name
    if( IWQZPL > 0 )then
      do NW = 1,NZOOPL
        if( NW < 10 )then
          write(SNUM,'(I1)') NW
          WQCONSTIT(19+NALGAE+NW) = 'ZOO'//SNUM
        else
          write(SNUM,'(I2)') NW
          WQCONSTIT(19+NALGAE+NW) = 'ZOO'//SNUM 
        endif
      enddo
    endif
    
    call fson_get(json_data, "active_constituents.ROC",                          ISKINETICS(IROC))
    call fson_get(json_data, "active_constituents.LOC",                          ISKINETICS(ILOC))
    call fson_get(json_data, "active_constituents.DOC",                          ISKINETICS(IDOC))
    call fson_get(json_data, "active_constituents.ROP",                          ISKINETICS(IROP))
    call fson_get(json_data, "active_constituents.LOP",                          ISKINETICS(ILOP))
    call fson_get(json_data, "active_constituents.DOP",                          ISKINETICS(IDOP))
    call fson_get(json_data, "active_constituents.P4D",                          ISKINETICS(IP4D))
    call fson_get(json_data, "active_constituents.RON",                          ISKINETICS(IRON))
    call fson_get(json_data, "active_constituents.LON",                          ISKINETICS(ILON))
    call fson_get(json_data, "active_constituents.DON",                          ISKINETICS(IDON))
    call fson_get(json_data, "active_constituents.NHX",                          ISKINETICS(INHX))
    call fson_get(json_data, "active_constituents.NOX",                          ISKINETICS(INOX))
    call fson_get(json_data, "active_constituents.SUU",                          ISKINETICS(ISUU))
    call fson_get(json_data, "active_constituents.SAA",                          ISKINETICS(ISAA))
    call fson_get(json_data, "active_constituents.COD",                          ISKINETICS(ICOD))
    call fson_get(json_data, "active_constituents.DOX",                          ISKINETICS(IDOX))
    call fson_get(json_data, "active_constituents.TAM",                          ISKINETICS(ITAM))
    call fson_get(json_data, "active_constituents.FCB",                          ISKINETICS(IFCB))
    call fson_get(json_data, "active_constituents.CO2",                          ISKINETICS(ICO2))
    If( ISWQLVL == 0  )then                                                      
      call fson_get(json_data, "active_constituents.CHC",                        ISKINETICS(ICHC))
      call fson_get(json_data, "active_constituents.CHD",                        ISKINETICS(ICHD))
      call fson_get(json_data, "active_constituents.CHG",                        ISKINETICS(ICHG))
    endif
    
    call fson_get(json_data, "sediment_diagenesis.benthic_flux_option",          IWQBEN)
    call fson_get(json_data, "sediment_diagenesis.number_of_reactive_classes",   NSMG)
    call fson_get(json_data, "rpem_activate",                                    ISRPEM)

    do NW = 1,NWQV
      call fson_get(json_data,'number_of_time_series.'//WQCONSTIT(NW),           NWQCSR(NW))
    Enddo
    
    if( IWQPSL == 2 )then
      ! *** Mass loading point source options
      call fson_get(json_data, "mass_loading_point_sources.constant_point_sources",   pointsource)
      do NPS = 1, fson_value_count(pointsource)
        item => fson_value_get(pointsource, NPS)
        call fson_get(item, "I", I)
        call fson_get(item, "J", J)
        call fson_get(item, "K", K)
        call fson_get(item, "NSR", ITMP)
        call fson_get(item, "PSQ", XPSQ)
        NCSERM = max(1,NCSERM,ITMP)
      enddo
    endif
  endif   ! *** End on master process
  
  call Broadcast_Scalar(ISWQLVL, master_id)
  call Broadcast_Scalar(NWQZ,    master_id)
  call Broadcast_Scalar(NWQTD,   master_id)
  call Broadcast_Scalar(NWQTS,   master_id)
  call Broadcast_Scalar(NTSWQV,  master_id)
  call Broadcast_Scalar(NSMZ,    master_id)
  call Broadcast_Scalar(NSMTS,   master_id)

  call Broadcast_Scalar(IWQPSL,  master_id)  
  call Broadcast_Scalar(NWQPS,   master_id)
  call Broadcast_Scalar(NPSTMSR, master_id)
  
  call Broadcast_Scalar(NALGAE,  master_id)
  call Broadcast_Scalar(IWQZPL,  master_id)
  call Broadcast_Scalar(NZOOPL,  master_id)
  
  call Broadcast_Scalar(ISFFARM, master_id)
  call Broadcast_Scalar(NSF,     master_id)
  call Broadcast_Scalar(NWQV,    master_id)
  
  call Broadcast_Scalar(IROC,    master_id)
  call Broadcast_Scalar(ILOC,    master_id)
  call Broadcast_Scalar(IDOC,    master_id)
  call Broadcast_Scalar(IROP,    master_id)
  call Broadcast_Scalar(ILOP,    master_id)
  call Broadcast_Scalar(IDOP,    master_id)
  call Broadcast_Scalar(IP4D,    master_id)
  call Broadcast_Scalar(IRON,    master_id)
  call Broadcast_Scalar(ILON,    master_id)
  call Broadcast_Scalar(IDON,    master_id)
  call Broadcast_Scalar(INHX,    master_id)
  call Broadcast_Scalar(INOX,    master_id)
  call Broadcast_Scalar(ISUU,    master_id)
  call Broadcast_Scalar(ISAA,    master_id)
  call Broadcast_Scalar(ICOD,    master_id)
  call Broadcast_Scalar(IDOX,    master_id)
  call Broadcast_Scalar(ITAM,    master_id)
  call Broadcast_Scalar(IFCB,    master_id)
  call Broadcast_Scalar(ICO2,    master_id)
  
  call Broadcast_Scalar(IWQBEN,  master_id)
  call Broadcast_Scalar(NSMG,    master_id)
  call Broadcast_Scalar(ISRPEM,  master_id)
  
  call Broadcast_Scalar(IWQTS,   master_id)
  call Broadcast_Scalar(WQTSDT,  master_id)

  call Broadcast_Array(ISKINETICS, master_id)
  call Broadcast_Array(NWQCSR,     master_id)

  NWQZM   = max(1,NWQZ)
  NWQPSM  = max(1,NWQPS)
  NWQTDM  = max(1,NWQTD) 
  NWQTSM  = max(1,NWQTS)
  NSMGM   = max(1,NSMG)
  NSMZM   = max(1,NSMZ)
  NSMTSM  = max(1,NSMTS)
  NALGAEM = max(1,NALGAE)
  ! NTSWQVM = max(1,NTSWQV)
  NWQVM = NWQV
    
  !If (ISFFARM > 0) NWQV = NWQV + NSF
    
  NWQPSRM = max(1,NPSTMSR)
  do NW = 1,NWQV
    NWQCSRM = max(NWQCSRM,NWQCSR(NW))
  enddo
  NCSERM = max(NCSERM,NWQCSRM)

  ISMOB = 0
  MACDRAG = 0
  IVARSETL = 0
  if( process_id == master_id .and. NALGAE > 0 )then  ! *** Only read wq_biota when model simulates algae - DKT
    ! *** SCAN WQ_BIOTA.JNP FILE
    write(*,'(A)') 'SCANNING EUTROPHICATION BIOTA CONFIGURATION'

    json_data => fson_parse("wq_biota.jnp")
    call fson_get(json_data, "groups", algaegroups)
    do NAL = 1, fson_value_count(algaegroups)
      item => fson_value_get(algaegroups, NAL)
      call fson_get(item, "mobility_flag", ISMOB(NAL)) 

      call fson_get(item, "use_hydro_feedback", IDRAG)
      if( IDRAG > 0 .and. ISVEG == 0 )then
        MACDRAG = 1
        ISVEG = 2                                   ! *** Gets overwritten in INPUT but used to flag allocation of vegetation arrays
      endif
      
      call fson_get(item, "variable_settling", NW) 
      if( NW == 3 )then
        IVARSETL = 3
      endif
    enddo
  endif   ! *** End on master process
    
  call Broadcast_Array( ISMOB,     master_id)
  call Broadcast_Scalar( IVARSETL, master_id)
  
  NFIXED = 0
  
  do NW = 1, 19
    ISTRWQ(NW) = ISKINETICS(NW)       ! *** Transport flag for nutrients
  enddo
  
  if( NALGAE > 0 )then
    do NAL = 1,NALGAE
      ISKINETICS(NAL+19) = 1            ! *** Kinetics Flag
      ISTRWQ(NAL+19) = 1                ! *** Transport Flag
      if( ISMOB(NAL) == 0 )then
        NFIXED = NFIXED + 1
        ISTRWQ(NAL+19) = 0              ! *** Do not transport this constituent
      endif
    enddo
  endif
  
  if( NZOOPL > 0 )then
    do NW = 1,NZOOPL
      ISKINETICS(NW+NALGAE+19) = 1      ! *** Kinetics Flag
      ISTRWQ(NW+NALGAE+19) = 1          ! *** Transport Flag
    enddo
  endif
  
  call Broadcast_Scalar(NFIXED,  master_id) 

  if( process_id == master_id )then 
    ! *** SCAN THE TIME SERIES
    if( NPSTMSR >= 1 .and. IWQPSL /= 2 )then
      call OPENFILE('WQPSL.INP')
      call SKIPCOM(1,'*')
      ERRSTR = 'READING WPQSL.INP'
      do NS = 1,NPSTMSR
        read(1,*,err = 999,end = 20) M, TM, TA, RMULADJ, ADDADJ
        NDWQPSR = max(NDWQPSR,M)
        do J = 1,M
          do K = 1,3
            read(1,'(A120)')LINE
          enddo
        enddo
      enddo
      20 close(1)
    endif
  endif   ! *** End on master process

  call Broadcast_Scalar(NDWQPSR,  master_id) 
    
  ! *** SCAN THE OPEN BC TIME SERIES
  ! ***  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
  ! ***  TIME SERIES FROM THE FILES WQCSRNN.INP
  ! *** SCAN FOR NUMBER OF SERIES   
  !IF (ISWQLVL == 0 )then
  !  do NW = 1,40
  !    write(SNUM,'(I2.2)')NW
  !    FNWQSR(NW)  = 'wqcsr'//SNUM//'.inp' 
  !  enddo
  !  write(6,'(A)')'SCANNING INPUT FILE: WQCSRXX.INP'
  !ENDIF
    
  if( process_id == master_id )then 
    if( ISWQLVL == 1  )then  
      ! *** Nutrient
      do NW = 1,19
        write(SNUM,'(I2.2)') NW
        FNWQSR(NW) = 'wqcsr'//SNUM//'.inp'
      enddo
      ! *** Algae
      do NW = 1,NALGAE
        write(SNUM,'(I2.2)')NW
        FNWQSR(NW+19) = 'wqalgsr'//SNUM//'.inp'
      enddo
      ! *** Zooplankton
      do NW = 1,NZOOPL
         write(SNUM,'(I2.2)')NW
         FNWQSR(NW+19+NALGAE) = 'wqzoosr'//SNUM//'.inp'
      enddo
      
      write(6,'(A)')'SCANNING WQ CONCENTRATION TIME SERIES INPUT FILES'
    endif
        
    NSER(8) = 0
    do NW = 1,NWQV
      INQUIRE(FILE = FNWQSR(NW),EXIST = BFLAG)
      if( BFLAG )then
        open(1,FILE = FNWQSR(NW))
        
        ERRSTR = 'READING '//TRIM(FNWQSR(NW))
        call SKIPCOM(1,'*')
        do NS = 1,1000
          read(1,*,err = 999,end = 40)ISTYP,M,T1,T2,RMULADJ,ADDADJ
          if( ISTYP == 1 ) M = M + 1   ! *** SKIP THE LAYER SPLITS
          do J = 1,M
            read(1,'(A120)')LINE
          enddo
        enddo
        40 close(1)
        NSER(8) = max(NSER(8),NS-1)
      endif
    enddo
  endif   ! *** End on master process
    
  call Broadcast_Scalar(NSER(8),  master_id) 
  
  allocate(TSWQ(NSER(8),NWQV))

  ! *** SCAN FOR NUMBER OF POINTS IN EACH SERIES
  if( process_id == master_id )then
    do NW = 1,NWQV
      ! *** ****************************
      INQUIRE(FILE = FNWQSR(NW),EXIST = BFLAG)
      if( BFLAG )then
        ! *** ****************************
        call OPENFILE(FNWQSR(NW))
        call SKIPCOM(1,'*')

        do NS = 1,NSER(8)
          ! *** ****************************
          read(1,*,end = 45) ISTYP, M, T1, T2, RMULADJ, ADDADJ
          TSWQ(NS,NW).NREC = M
          
          if( ISTYP == 1 ) M = M+1   ! *** SKIP THE LAYER SPLITS
          ! *** ****************************
          do J = 1,M
            read(1,'(A120)')LINE
          enddo
        enddo
      endif
45    close(1)
    enddo
  endif ! *** End on master process
  
  NDWQCSR = 1
  do NW = 1,NWQV
    do NS = 1,NSER(8)
      call Broadcast_Scalar(TSWQ(NS,NW).NREC, master_id)
      M = TSWQ(NS,NW).NREC
      NDWQCSR = max(NDWQCSR,M)

      allocate(TSWQ(NS,NW).VAL(M,KCM), TSWQ(NS,NW).TIM(M))
      TSWQ(NS,NW).VAL = 0.
      TSWQ(NS,NW).TIM = 0.
    enddo
  enddo
  
  NDCSER = max(NDCSER,NDWQCSR)
  
  ! *** ****************************
  ! *** Sediment Diagenesis
  if( IWQBEN == 1 )then
    if( process_id == master_id )then
      write(*,'(A)') 'SCANNING SEDIMENT DIAGENESIS CONFIGURATION'
      json_data => fson_parse("wq_3dsd.jnp")

      call fson_get(json_data, "initial_condition_option", ISMICI)
      call fson_get(json_data, "number_of_spatial_zones", ISMZ)
      call fson_get(json_data, "write_restart_option", ISMRST)
      call fson_get(json_data, "benthic_stress.activate_hysteresis_in_benthic_mixing", ISMHYST)
      call fson_get(json_data, "activate_diagnostic_output", ISMZB)
      
      call fson_get(json_data, "time_series_output.number_of_locations", ISMTS)
    endif  ! *** End on master process

    call Broadcast_Scalar(ISMICI,  master_id)
    call Broadcast_Scalar(ISMZ,    master_id)
    call Broadcast_Scalar(ISMRST,  master_id)
    call Broadcast_Scalar(ISMHYST, master_id)
    call Broadcast_Scalar(ISMZB,   master_id)
    call Broadcast_Scalar(ISMTS,   master_id)

    NSMTSM = max(ISMTS,NSMTS)

  endif

  return
999 PRINT*,ERRSTR
  STOP

    
END SUBROUTINE SCANWQ


SUBROUTINE SCANSEDZLJ
  !  REVISION DATE :  May 24, 2006
  !  Craig Jones and Scott James
  ! *** ************************************************************
  integer :: IDUMMY,ERROR
  ! *** ****************************
  if( process_id == master_id )then
    call OPENFILE('BED.SDF')
    read(1,*,IOSTAT = ERROR) !SKIP THIS LINE
    if( ERROR == 1 )then
      write(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      write(mpi_error_unit,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    endif
    read(1,*,IOSTAT = ERROR) IDUMMY, KB, ICALC_BL, SEDSTEP, SEDSTART, IHTSTRT, IMORPH, ISWNWAVE, MAXDEPLIMIT, HPMIN
    if( ERROR == 1 )then
      write(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      write(mpi_error_unit,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    endif
    read(1,*,IOSTAT = ERROR) !SKIP THIS LINE
    if( ERROR == 1 )then
      write(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      write(mpi_error_unit,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    endif
    read(1,*,IOSTAT = ERROR) ITBM,NSICM
    if( ERROR == 1 )then
      write(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      write(mpi_error_unit,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    endif
    close(1)
  endif   ! *** End on master process

  call Broadcast_Scalar(KB,       master_id)
  call Broadcast_Scalar(ICALC_BL, master_id)
  call Broadcast_Scalar(ITBM,     master_id)
  call Broadcast_Scalar(NSICM,    master_id)

END SUBROUTINE SCANSEDZLJ

SUBROUTINE SCNTXSED

  integer :: NCSERNC,NLOOP,IS,NS,ISTYP,NDATAPTS,NN,NT,M
  character*120 SKIP
  character*10 INFILE
  character*80 STR*200


  ! *** NOW FIND MAX FOR TOXICS AND SEDIMENTS
  do NN = 1,3
    NCSERNC = 0
    ! *** ****************************
    if( process_id == master_id )then
      if( NN == 1 )then
        if( NTOX > 0 .and. NSER(5) > 0 .and. ISTRAN(5) > 0 )then
          call OPENFILE('TXSER.INP')
          NLOOP = NTOX
          NCSERNC = NSER(5)
        endif
      elseif( NN == 2 )then
        if( NSED > 0 .and. NSER(6) > 0 .and. ISTRAN(6) > 0 )then
          call OPENFILE('SDSER.INP')
          NLOOP = NSED
          NCSERNC = NSER(6)
        endif
      elseif( NN == 3 )then
        if( NSND > 0 .and. NSER(7) > 0 .and. ISTRAN(7) > 0 )then
          call OPENFILE('SNSER.INP')
          NLOOP = NSND
          NCSERNC = NSER(7)
        endif
      endif
    endif

    call Broadcast_Scalar(NLOOP, master_id)
    call Broadcast_Scalar(NCSERNC, master_id)

    if( NCSERNC > 0 )then
      ! *** ****************************
      if( process_id == master_id )then
        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      endif   ! *** End on master process

      ! *** LOOP OVER THE NUMBER OF SERIES
      do NS = 1,NSER(NN+4)  !NCSERNC
        ! *** ****************************
        if( process_id == master_id )then
          read(1,*,err = 20,end = 40)ISTYP,NDATAPTS   !,X1,X2,X3,X4
        endif   ! *** End on master process

        call Broadcast_Scalar(ISTYP, master_id)
        call Broadcast_Scalar(NDATAPTS, master_id)

        if( NN == 1 .and. NTOX > 0 .and. NSER(5) > 0 )then
          do NT = 1,NTOX
            allocate(TSTOX(NS,NT).VAL(NDATAPTS,KCM),TSTOX(NS,NT).TIM(NDATAPTS))
            TSTOX(NS,NT).NREC = NDATAPTS
            TSTOX(NS,NT).VAL = 0
            TSTOX(NS,NT).TIM = 0
          enddo

        elseif( NN == 2 .and. NSED > 0 .and. NSER(6) > 0 )then
          do NT = 1,NSED
            allocate(TSSED(NS,NT).VAL(NDATAPTS,KCM),TSSED(NS,NT).TIM(NDATAPTS))
            TSSED(NS,NT).NREC = NDATAPTS
            TSSED(NS,NT).VAL = 0
            TSSED(NS,NT).TIM = 0
          enddo
        elseif( NN == 3 .and. NSND > 0 .and. NSER(7) > 0 )then
          do NT = 1,NSND
            allocate(TSSND(NS,NT).VAL(NDATAPTS,KCM),TSSND(NS,NT).TIM(NDATAPTS))
            TSSND(NS,NT).NREC = NDATAPTS
            TSSND(NS,NT).VAL = 0
            TSSND(NS,NT).TIM = 0
          enddo
        endif
        ! *** ****************************
        if( process_id == master_id )then
          ! *** SKIP THE CONVERSIONS
          if( NN /= 1 )then
            do NT = 2,NLOOP
              read(1,'(A)') SKIP
            enddo
          endif

          NDCSER = max(NDCSER,NDATAPTS)
          if( ISTYP == 1 )then
            ! *** SKIP THE SPLITS
            read(1,'(A)') SKIP
          endif
          ! *** SKIP TO THE NEXT SERIES
          do M = 1,NDATAPTS
            do NT = 1,NLOOP
              read(1,'(A)') SKIP
            enddo
          enddo
        endif   ! *** End on master process

      enddo
      ! *** ****************************
      if( process_id == master_id )then
        close(1)
      endif   ! *** End on master process

    endif
  enddo
  return
   ! *** ****************************
20 WRITE(6,'(A)') '*** READ ERROR IN FILE: '//CFILE
   write(mpi_error_unit,'(A)')    '*** READ ERROR IN FILE: '//CFILE
   call STOPP('.')

40 WRITE(6,'(A)') '*** UNEXPECTED END OF FILE: '//CFILE
   write(mpi_error_unit,'(A)')    '*** UNEXPECTED END OF FILE: '//CFILE
   call STOPP('.')

END SUBROUTINE SCNTXSED

SUBROUTINE SCANPROPWASH

  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get

  type(fson_value), Pointer :: json_data
  
  ! *** ****************************
  if( process_id == master_id )then
    ! *** Scan the propwash config file
    write(*,'(A)') 'SCANNING PROPWASH CONFIGURATION'

    ! *** Open the propwash ship data file
    json_data => fson_parse("propwash_config.jnp")
    
    ! *** Fast Settling Classes
    call fson_get(json_data, "parms.fraction_fast",   fraction_fast)
    call fson_get(json_data, "parms.fast_multiplier", fast_multiplier)
  endif

  call Broadcast_Scalar(fraction_fast,   master_id)
  call Broadcast_Scalar(fast_multiplier, master_id)

END SUBROUTINE SCANPROPWASH

SUBROUTINE OPENFILE(INFILE)

  character(*), intent(IN) :: INFILE
  CFILE = StrUpCase(INFILE)

  write(6,'(A)')'SCANNING INPUT FILE: '//CFILE

  CFILE = StrLowCase(INFILE)
  open(1,FILE = CFILE,STATUS = 'OLD',err = 200)

  return
200 WRITE(6,'(A)') '  FILE DOES NOT EXIST:  '//TRIM(CFILE)
  write(mpi_error_unit,'(A)') '  FILE DOES NOT EXIST:  '//TRIM(CFILE)
  call STOPP('.')

END SUBROUTINE

FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )
  ! FROM "String_Utility" BY Paul van Delst
  ! -- Argument and result
  character( * ), intent( IN )     :: Input_String
  character( LEN( Input_String ) ) :: Output_String
  ! -- Local variables
  integer :: i, n

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  do i = 1, LEN( Output_String )
    ! -- Find location of letter in lower case constant string
    n = INDEX( LOWER_CASE, Output_String( i:i ) )
    ! -- If current substring is a lower case letter, make it upper case
    if( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
  enddo
END FUNCTION StrUpCase

FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String )
  ! FROM "String_Utility" BY Paul van Delst
  ! -- Argument and result
  character( * ), intent( IN )     :: Input_String
  character( LEN( Input_String ) ) :: Output_String
  ! -- Local variables
  integer :: i, n

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  do i = 1, LEN( Output_String )
    ! -- Find location of letter in upper case constant string
    n = INDEX( UPPER_CASE, Output_String( i:i ) )
    ! -- If current substring is an upper case letter, make it lower case
    if( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
  enddo
END FUNCTION StrLowCase

END MODULE
