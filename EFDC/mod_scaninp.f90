! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE SCANINPMOD
  
  USE GLOBAL
  USE INFOMOD, ONLY:SKIPCOM,READSTR
  Use Allocate_Initialize
  Use Variables_WQ
  USE HYDSTRUCMOD
  Use Variables_Propwash

  Use Variables_MPI
  Use Variables_MPI_Write_Out
  Use Broadcast_Routines
  
  IMPLICIT NONE

  INTEGER :: NSER(8)
  CHARACTER(*) :: STR*200 ,CFILE*15
  CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  CONTAINS

SUBROUTINE SCANEFDC(NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7)

  INTEGER, INTENT(OUT) :: NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7

  !INTEGER :: K, NT, I, M, M1, NN
  !INTEGER :: NDUM, ISTOCNT, ISO, ITIDASM, NPFOR, NDATAPTS, ISTYP
  INTEGER :: K, NT, I, M, M1, NN, ITYPE
  INTEGER :: NDUM, ISTOCNT, ISO, ITIDASM, NPFOR, NDATAPTS
  INTEGER :: ITMPPMX, ITSSS, NS, IS, IITMP, IJTMP, NPMXZ, NPMXPTS
  INTEGER :: ISBEDTEMI  ! *** DEPRECATED VARIABLE
  REAL    :: DIFTOXNT, DIFTOXSNT
  REAL    :: PDIFTOXNT, DPDIFTOXNT, R, TC1, TAV1
  REAL    :: R1TMP, R2TMP, R3TMP, R4TMP, R5TMP, PMIXSF

  ! *** ********************************************
  if( process_id == master_id )THEN
    CALL OPENFILE('EFDC.INP')

    CALL SEEK('C1')
    READ(1,'(A100)') RUNTITLE

    CALL SEEK('C1A')
    READ(1,*,ERR=10) IS2TIM,IGRIDH,IGRIDV,KMINV,SGZHPDELTA
  end if

  Call Broadcast_Scalar(IS2TIM,     master_id)
  Call Broadcast_Scalar(IGRIDH,     master_id)
  Call Broadcast_Scalar(IGRIDV,     master_id)
  Call Broadcast_Scalar(KMINV,      master_id)
  Call Broadcast_Scalar(SGZHPDELTA, master_id)

  if( process_id == master_id )THEN
    CALL SEEK('C2')
    READ(1,*,ERR=10) ISRESTI,ISRESTO,ISRESTR,ISGREGOR,ISLOG,ISDIVEX,ISNEGH,ISMMC,ISBAL,ICONTINUE,ISHOW
  end if

  Call Broadcast_Scalar(ISRESTI,  master_id)
  Call Broadcast_Scalar(ISRESTO,  master_id)
  Call Broadcast_Scalar(ISRESTR,  master_id)
  Call Broadcast_Scalar(ISGREGOR, master_id)
  Call Broadcast_Scalar(ISLOG,    master_id)
  Call Broadcast_Scalar(ISDIVEX,  master_id)
  Call Broadcast_Scalar(ISNEGH,   master_id)
  Call Broadcast_Scalar(ISMMC,    master_id)
  Call Broadcast_Scalar(ISBAL,    master_id)
  Call Broadcast_Scalar(ICONTINUE,master_id)
  Call Broadcast_Scalar(ISHOW, master_id)

  if( process_id == master_id )THEN
    CALL SEEK('C4')
    READ(1,*,ERR=10) ISLTMT,ISSSMMT,RESSTEP
  end if

  Call Broadcast_Scalar(ISLTMT, master_id)
  Call Broadcast_Scalar(ISSSMMT,master_id)
  Call Broadcast_Scalar(RESSTEP,master_id)

  IF( ISLTMT > 0 )THEN
    CALL STOPP('ISLTMT LONG TERM SIMULATION DISABLED!')
  ENDIF

  if( process_id == master_id )THEN
    CALL SEEK('C5')
    READ(1,*,ERR=10) ISCDMA,ISHDMF,ISDISP,ISWASP,ISDRY,ICALTB,ISRLID,ISVEG,ISVEGL,ISITB,IHMDSUB,IINTPG
  end if

  Call Broadcast_Scalar(ISCDMA,  master_id)
  Call Broadcast_Scalar(ISHDMF,  master_id)
  Call Broadcast_Scalar(ISDISP,  master_id)
  Call Broadcast_Scalar(ISWASP,  master_id)
  Call Broadcast_Scalar(ISDRY,   master_id)
  Call Broadcast_Scalar(ICALTB,  master_id)
  Call Broadcast_Scalar(ISRLID,  master_id)
  Call Broadcast_Scalar(ISVEG,   master_id)
  Call Broadcast_Scalar(ISVEGL,  master_id)
  Call Broadcast_Scalar(ISITB,   master_id)
  Call Broadcast_Scalar(IHMDSUB, master_id)
  Call Broadcast_Scalar(IINTPG,  master_id)

  IF( ISCDMA == 10 )THEN
    CALL STOPP('EFDC 1D CHANNEL MODE IS DISABLED IN THIS VERSION OF EFDC')
  ENDIF
  IF( ISCDMA > 2 )THEN
    CALL STOPP('EXPERIMENTAL MOMENTUM OPTIONs ARE ARE DISABLED IN THIS VERSION OF EFDC')
  ENDIF

  if( process_id == master_id )THEN
    CALL SEEK('C6')
    DO NN=0,8
      READ(1,*,ERR=10) ISTRAN(NN),ISTOPT(NN),ISCDCA(NN),ISADAC(NN),ISFCT(NN),ISPLIT(NN),ISADAH(NN),ISADAV(NN),ISCI(NN),ISCO(NN)
    ENDDO
  end if

  Call Broadcast_Array(ISTRAN, master_id)
  Call Broadcast_Array(ISTOPT, master_id)
  Call Broadcast_Array(ISCDCA, master_id)
  Call Broadcast_Array(ISADAC, master_id)
  Call Broadcast_Array(ISFCT , master_id)
  Call Broadcast_Array(ISPLIT, master_id)
  Call Broadcast_Array(ISADAH, master_id)
  Call Broadcast_Array(ISADAV, master_id)
  Call Broadcast_Array(ISCI  , master_id)
  Call Broadcast_Array(ISCO  , master_id)

  if( process_id == master_id )THEN
    CALL SEEK('C7')
    READ(1,*,ERR=10) NTC,NTSPTC,NLTC,NTTC,NTCPP,NTSTBC,NTCNB,NTCVB,NTSMMT,NFLTMT,NDRYSTP
  end if

  Call Broadcast_Scalar(NTC,     master_id)
  Call Broadcast_Scalar(NTSPTC,  master_id)
  Call Broadcast_Scalar(NLTC,    master_id)
  Call Broadcast_Scalar(NTTC,    master_id)
  Call Broadcast_Scalar(NTCPP,   master_id)
  Call Broadcast_Scalar(NTSTBC,  master_id)
  Call Broadcast_Scalar(NTCNB,   master_id)
  Call Broadcast_Scalar(NTCVB,   master_id)
  Call Broadcast_Scalar(NTSMMT,  master_id)
  Call Broadcast_Scalar(NFLTMT,  master_id)
  Call Broadcast_Scalar(NDRYSTP, master_id)

  NDDAM=NTC

  if( process_id == master_id )THEN
    CALL SEEK('C9')
    READ(1,*,ERR=10,END=30)IC,JC,LC,LVC,ISCLO,NDM,LDM,ISMASK,NBLOCKED,ISCONNECT,NSHMAX,NSBMAX,WSMH,WSMB
  end if

  Call Broadcast_Scalar(IC,       master_id)
  Call Broadcast_Scalar(JC,       master_id)
  Call Broadcast_Scalar(LC,       master_id)
  Call Broadcast_Scalar(LVC,      master_id)
  Call Broadcast_Scalar(ISCLO,    master_id)
  Call Broadcast_Scalar(NDM,      master_id)
  Call Broadcast_Scalar(LDM,      master_id)
  Call Broadcast_Scalar(ISMASK,   master_id)
  Call Broadcast_Scalar(NBLOCKED, master_id)
  Call Broadcast_Scalar(ISCONNECT,master_id)
  Call Broadcast_Scalar(NSHMAX, master_id)
  Call Broadcast_Scalar(NSBMAX, master_id)
  Call Broadcast_Scalar(WSMH,   master_id)
  Call Broadcast_Scalar(WSMB,   master_id)

  ! *** TO ENSURE CONSISTENCY OF DECLARATIONS BETWEEN SUB-DOMAINS
  ISBLOCKED = 0
  IF( NBLOCKED > 0 ) ISBLOCKED = 1
  
  IF( IC >= 3 )THEN
    ICM = IC + 1
    ICM_Global = ICM
  ELSE
    CALL STOPP('IC MUST BE AT LEAST 3')
  ENDIF
  IF( JC >= 3 )THEN
    JCM = JC + 1
    JCM_Global = JCM
  ELSE
    CALL STOPP('IJ MUST BE AT LEAST 3')
  ENDIF
  IF( LC >= 3 )THEN
    LCM = LC + 1
    LCM_Global = LCM !***Save the 'true' global LCM value based on what was read in from the input
    LC_Global = LC
  ELSE
    CALL STOPP('LC MUST BE AT LEAST 3')
  ENDIF

  if( process_id == master_id )THEN
    CALL SEEK('C9A')
    READ(1,*,IOSTAT=ISO)KC,KSIG,ISETGVC,SELVREF,BELVREF,ISGVCCK
  end if

  Call Broadcast_Scalar(KC,      master_id)
  Call Broadcast_Scalar(KSIG,    master_id)
  Call Broadcast_Scalar(ISETGVC, master_id)
  Call Broadcast_Scalar(SELVREF, master_id)
  Call Broadcast_Scalar(BELVREF, master_id)
  Call Broadcast_Scalar(ISGVCCK, master_id)

  KC=ABS(KC)
  IF( KC >= 1 )THEN
    KCM=KC+1
  ELSE
    CALL STOPP('KC MUST BE AT LEAST 1')
  ENDIF

  if( process_id == master_id )THEN
    CALL SEEK('C12A')
    READ(1,*,IOSTAT=ISO)ISTOPT(0),ISSQL,ISAVBMX,ISFAVB,ISINWV,ISLLIM,IFPROX,ISVTURB,BC_EDGEFACTOR
  end if

  Call Broadcast_Scalar(ISTOPT(0),     master_id)
  Call Broadcast_Scalar(ISSQL,         master_id)
  Call Broadcast_Scalar(ISAVBMX,       master_id)
  Call Broadcast_Scalar(ISFAVB,        master_id)
  Call Broadcast_Scalar(ISINWV,        master_id)
  Call Broadcast_Scalar(ISLLIM,        master_id)
  Call Broadcast_Scalar(IFPROX,        master_id)
  Call Broadcast_Scalar(ISVTURB,       master_id)
  Call Broadcast_Scalar(BC_EDGEFACTOR, master_id)

  if( process_id == master_id )THEN
    CALL SEEK('C14')
    READ(1,*,ERR=10,END=30) MTIDE, NWSER, NASER, ISGWIT, ISCHAN, ISWAVE, ITIDASM, ISPERC, ISBODYF, ISPNHYDS, ISPROPWASH
  end if

  Call Broadcast_Scalar(MTIDE,      master_id)
  Call Broadcast_Scalar(NWSER,      master_id)
  Call Broadcast_Scalar(NASER,      master_id)
  Call Broadcast_Scalar(ISGWIT,     master_id)
  Call Broadcast_Scalar(ISCHAN,     master_id)
  Call Broadcast_Scalar(ISWAVE,     master_id)
  Call Broadcast_Scalar(ITIDASM,    master_id)
  Call Broadcast_Scalar(ISPERC,     master_id)
  Call Broadcast_Scalar(ISBODYF,    master_id)
  Call Broadcast_Scalar(ISPNHYDS,   master_id)
  Call Broadcast_Scalar(ISPROPWASH, master_id)

  MTM=MAX(1,MTIDE)+1
  NWSERM=MAX(1,NWSER)
  NASERM=MAX(1,NASER)
  NGWSERM=1
  NDASER=1
  NDGWSER=1

  if( process_id == master_id )THEN
    CALL SEEK('C16')
    READ(1,*,ERR=10,END=30)NPBS,NPBW,NPBE,NPBN,NPFOR,NPFORT,NPSER,PDGINIT
  end if

  Call Broadcast_Scalar(NPBS,   master_id)
  Call Broadcast_Scalar(NPBW,   master_id)
  Call Broadcast_Scalar(NPBE,   master_id)
  Call Broadcast_Scalar(NPBN,   master_id)
  Call Broadcast_Scalar(NPFOR,  master_id)
  Call Broadcast_Scalar(NPFORT, master_id)
  Call Broadcast_Scalar(NPSER,  master_id)
  Call Broadcast_Scalar(PDGINIT,master_id)

  NPBSM=MAX(1,NPBS)
  NPBWM=MAX(1,NPBW)
  NPBEM=MAX(1,NPBE)
  NPBNM=MAX(1,NPBN)
  NPSERM=MAX(1,NPSER)
  NPFORM=MAX(1,NPFOR,NPSER)
  NDPSER=1

  if( process_id == master_id )THEN
    CALL SEEK('C22')
    READ(1,*,ERR=10,END=30)NDYE,NTOX,NSED,NSND,NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7,ISSBAL

    NDYE = MAX(1,NDYE)
    NDYM = MAX(1,NDYE)
    NTXM = MAX(1,NTOX)
    NSCM = MAX(1,NSED)
    NSNM = MAX(1,NSND)
    
    NCSERM = MAX(1,NCSER1,NCSER2,NCSER3,NCSER4,NCSER5,NCSER6,NCSER7)

  end if!***End on master process

  NDCSER=1

  Call Broadcast_Scalar(JC, master_id)
  Call Broadcast_Scalar(IC, master_id)
  Call Broadcast_Scalar(JCM, master_id)
  Call Broadcast_Scalar(ICM, master_id)

  Call Broadcast_Scalar(NDYE, master_id)
  Call Broadcast_Scalar(NTOX, master_id)
  Call Broadcast_Scalar(NSED, master_id)
  Call Broadcast_Scalar(NSND, master_id)
  Call Broadcast_Scalar(NDYE, master_id)
  Call Broadcast_Scalar(NDYM, master_id)
  Call Broadcast_Scalar(NTXM, master_id)
  Call Broadcast_Scalar(NSCM, master_id)
  Call Broadcast_Scalar(NSNM, master_id)
  Call Broadcast_Scalar(NCSERM, master_id)
  Call Broadcast_Scalar(NCSER1, master_id)
  Call Broadcast_Scalar(NCSER2, master_id)
  Call Broadcast_Scalar(NCSER3, master_id)
  Call Broadcast_Scalar(NCSER4, master_id)
  Call Broadcast_Scalar(NCSER5, master_id)
  Call Broadcast_Scalar(NCSER6, master_id)
  Call Broadcast_Scalar(NCSER7, master_id)

  ALLOCATE(TSSAL(NCSER1),TSTEM(NCSER2),TSDYE(NCSER3,NDYE),TSSFL(NCSER4))
  ALLOCATE(TSTOX(NCSER5,NTOX),TSSED(NCSER6,NSED),TSSND(NCSER7,NSND))

  NSER(1) = NCSER1
  NSER(2) = NCSER2
  NSER(3) = NCSER3
  NSER(4) = NCSER4
  NSER(5) = NCSER5
  NSER(6) = NCSER6
  NSER(7) = NCSER7

  ALLOCATE(TOXDEP(NTXM))

  if( process_id == master_id )THEN
    CALL SEEK('C23')
    READ(1,*,ERR=10,END=30) NQSIJ, NQJPIJ, NQSER, NQCTL, NQCTLT, NHYDST, NQWR, NQWRSR, ISDIQ, NQCTLSER, NQCTRULES
  end if

  Call Broadcast_Scalar(NQSIJ,    master_id)
  Call Broadcast_Scalar(NQJPIJ,   master_id)
  Call Broadcast_Scalar(NQSER,    master_id)
  Call Broadcast_Scalar(NQCTL,    master_id)
  Call Broadcast_Scalar(NQCTLT,   master_id)
  Call Broadcast_Scalar(NHYDST,   master_id)
  Call Broadcast_Scalar(NQWR,     master_id)
  Call Broadcast_Scalar(NQWRSR,   master_id)
  Call Broadcast_Scalar(ISDIQ,    master_id)
  Call Broadcast_Scalar(NQCTLSER, master_id)
  Call Broadcast_Scalar(NQCTRULES, master_id)

  NQSIJM=MAX(1,NQSIJ)
  NQJPM=MAX(1,NQJPIJ)
  NJPSM=NQJPM
  NQSERM=MAX(1,NQSER)
  NQCTLM=MAX(1,NQCTL)
  NQCTTM=MAX(1,NQCTLT)
  NQWRM=MAX(1,NQWR)
  NQWRSRM=MAX(1,NQWRSR)
  NDQSER=1   ! *** Flow              : Maximum number of  points in a series
  NDQWRSR=1  ! *** Withdrawal/Return : Maximum number of  points in a series

  ! *** SET KB AND CHECK FOR SEDZLJ USEAGE
  NSEDFLUME = 0
  NWARNING = 0
  LSEDZLJ = .FALSE.
  IF( NSED > 0 .OR. NSND > 0 )THEN
    if( process_id == master_id )THEN
      CALL SEEK('C36')
      READ(1,*,ERR=10,END=30)ISEDINT,ISEDBINT,NSEDFLUME,ISMUD,ISBEDMAP,ISEDVW,ISNDVW,KB,ISDTXBUG
    end if

    Call Broadcast_Scalar(ISEDINT,   master_id)
    Call Broadcast_Scalar(ISEDBINT,  master_id)
    Call Broadcast_Scalar(NSEDFLUME, master_id)
    Call Broadcast_Scalar(ISMUD,     master_id)
    Call Broadcast_Scalar(ISBEDMAP,  master_id)
    Call Broadcast_Scalar(ISEDVW,    master_id)
    Call Broadcast_Scalar(ISNDVW,    master_id)
    Call Broadcast_Scalar(KB,        master_id)
    Call Broadcast_Scalar(ISDTXBUG,  master_id)

    IF( KB >= 1 )THEN
      KBM=KB+1
    ELSE
      CALL STOPP('KB MUST BE AT LEAST 1')
    ENDIF
    IF( ISTRAN(6) > 0 )THEN
      IF( NSEDFLUME == 98 .OR. NSEDFLUME == 99 )THEN
        LSEDZLJ = .TRUE.
        NSEDFLUME = 1
      ELSEIF( NSEDFLUME > 0 )THEN
        LSEDZLJ = .TRUE.
      ENDIF
    ENDIF

    ! *** LOOK FOR LEGACY APPROACH TO SET SEDZLJ, I.E. IWRSP(1)=98
    IF( .NOT. LSEDZLJ .AND. ISTRAN(6) > 0 )THEN
      ISDTXBUG=0
      !C40*  READ COHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSED TIMES
      if( process_id == master_id )THEN
        CALL SEEK('C40')
        READ(1,*,ERR=10,END=30)ISDTXBUG
      end if

      Call Broadcast_Scalar(ISDTXBUG, master_id)

      IF( ISDTXBUG == 98 )THEN
        LSEDZLJ = .TRUE.
        NSEDFLUME=1
      ENDIF
    ENDIF
    
    IF( .NOT. LSEDZLJ .AND. ISTRAN(7) > 0 .AND. NSND > 0 )THEN
      if( process_id == master_id )THEN
        CALL SEEK('C42A')
        DO NS=1,1
          READ(1,*,IOSTAT=ISO) ICALC_BL, R, R, R, R, R, R, R, R, R
        ENDDO
      endif
      
      Call Broadcast_Scalar(ICALC_BL , master_id)
    ENDIF
    
    IF( .NOT. LSEDZLJ )THEN
      ! *** SET NSCM TO THE MAXIMUM NUMBER OF SEDIMENT CLASSES
      NSCM = MAX(NSCM, NSED+NSND)
    ENDIF
  ELSE
    KBM=1
  ENDIF
  
  ITMPPMX=0
  IF( NTOX > 0 )THEN
  
    !C43C*  READ TOXIC TIMING AND VOLATILIZATION SWITCHES
    if(process_id == master_id)THEN
        CALL SEEK('C43C')
        READ(1,*,IOSTAT=ISO) TOXSTEPW, TOXSTEPB, TOX_VEL_MAX, TOX_DEP_MAX, ITOXTEMP, TOXTEMP
    endif
    
    Call Broadcast_Scalar(TOXSTEPW, master_id)
    Call Broadcast_Scalar(TOXSTEPB, master_id)
    Call Broadcast_Scalar(TOX_VEL_MAX, master_id)
    Call Broadcast_Scalar(TOX_DEP_MAX, master_id)
    Call Broadcast_Scalar(ITOXTEMP, master_id)
    Call Broadcast_Scalar(TOXTEMP, master_id)
    
    if( process_id == master_id )THEN
      CALL SEEK('C44')
    endif
    ! *** NEED TO READ EVEN IF SEDZLJ IS BEING USED
    DO NT=1,NTOX
      if( process_id == master_id )THEN
        READ(1,*,ERR=10)NDUM,ISTOCNT,DIFTOXNT,DIFTOXSNT,PDIFTOXNT,DPDIFTOXNT
      end if

      Call Broadcast_Scalar(ISTOCNT,   master_id)
      Call Broadcast_Scalar(DIFTOXNT,  master_id)
      Call Broadcast_Scalar(DIFTOXSNT, master_id)
      Call Broadcast_Scalar(PDIFTOXNT, master_id)
      Call Broadcast_Scalar(DPDIFTOXNT,master_id)

      IF( PDIFTOXNT < 0. )ITMPPMX=1
    ENDDO

    if( process_id == master_id )THEN
      CALL SEEK('C45A')
      READ(1,*,ERR=10)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
    end if

    Call Broadcast_Scalar(ISTDOCW, master_id)
    Call Broadcast_Scalar(ISTPOCW, master_id)
    Call Broadcast_Scalar(ISTDOCB, master_id)
    Call Broadcast_Scalar(ISTPOCB, master_id)
    Call Broadcast_Scalar(STDOCWC, master_id)
    Call Broadcast_Scalar(STPOCWC, master_id)
    Call Broadcast_Scalar(STDOCBC, master_id)
    Call Broadcast_Scalar(STPOCBC, master_id)

    if( process_id == master_id )THEN
      CALL SEEK('C45E')
    end if

    DO NT=1,NTOX
      if( process_id == master_id )THEN
        READ(1,*,ERR=10) NDUM, TOXDEP(NT).ITXDRY, TOXDEP(NT).TXDRY, TOXDEP(NT).ITXWET, TOXDEP(NT).TXWET
      end if

      Call Broadcast_Scalar(TOXDEP(NT).ITXDRY, master_id)
      Call Broadcast_Scalar(TOXDEP(NT).ITXWET, master_id)
    ENDDO
  ENDIF

  if( process_id == master_id )THEN
    CALL SEEK('C46')
    READ(1,*,ERR=10,END=30) BSC,TEMO,HEQT,RKDYE,NCBS,NCBW,NCBE,NCBN
  end if

  Call Broadcast_Scalar(BSC,       master_id)
  Call Broadcast_Scalar(TEMO,      master_id)
  Call Broadcast_Scalar(HEQT,      master_id)
  Call Broadcast_Scalar(RKDYE,     master_id)
  Call Broadcast_Scalar(NCBS,      master_id)
  Call Broadcast_Scalar(NCBW,      master_id)
  Call Broadcast_Scalar(NCBE,      master_id)
  Call Broadcast_Scalar(NCBN,      master_id)

  NBBSM=MAX(1,NCBS)
  NBBWM=MAX(1,NCBW)
  NBBEM=MAX(1,NCBE)
  NBBNM=MAX(1,NCBN)

  if( process_id == master_id )THEN
    CALL SEEK('C46A')
    READ(1,*,ERR=10,END=30)ISICE,NISER,TEMPICE,CDICE,ICETHMX,RICETHK0
  end if

  call Broadcast_Scalar(ISICE,    master_id)
  call Broadcast_Scalar(NISER,    master_id)
  call Broadcast_Scalar(TEMPICE,  master_id)
  call Broadcast_Scalar(CDICE,    master_id)
  call Broadcast_Scalar(ICETHMX,  master_id)
  call Broadcast_Scalar(RICETHK0, master_id)

  IF(  NASER > 0  )THEN
    if( process_id == master_id )THEN
      CALL SEEK('C46D')
      READ(1,*,IOSTAT=ISO) IASWRAD,REVC,RCHC,ISVHEAT,SWRATNF,SWRATNS,FSWRATF,TEMTHKO,TEMBO,HTBED1,HTBED2
    end if

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

  ENDIF

  if( process_id == master_id )THEN
    CALL SEEK('C66A')
    READ(1,*,ERR=10) NLCDA,TSCDA,(ISCDA(K),K=1,7)
  end if

  Call Broadcast_Scalar(NLCDA, master_id)
  Call Broadcast_Scalar(TSCDA, master_id)
  Call Broadcast_Array (ISCDA, master_id)
  NLDAM=NLCDA

  if( process_id == master_id )THEN
    CALL SEEK('C67')
    READ(1,*,ERR=10) ISPD
  end if

  Call Broadcast_Scalar(ISPD,   master_id)

  NPDM=MAX(1,NPD)

  if( process_id == master_id )THEN
    CALL SEEK('C80')
    READ(1,*,ERR=10,END=30)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
  end if

  Call Broadcast_Scalar(IS3DO,  master_id)
  Call Broadcast_Scalar(ISR3DO, master_id)
  Call Broadcast_Scalar(NP3DO,  master_id)
  Call Broadcast_Scalar(KPC,    master_id)
  Call Broadcast_Scalar(NWGG,   master_id)
  Call Broadcast_Scalar(I3DMIN, master_id)
  Call Broadcast_Scalar(I3DMAX, master_id)
  Call Broadcast_Scalar(J3DMIN, master_id)
  Call Broadcast_Scalar(J3DMAX, master_id)
  Call Broadcast_Scalar(I3DRW,  master_id)
  Call Broadcast_Scalar(SELVMAX,master_id)
  Call Broadcast_Scalar(BELVMIN,master_id)

  KPCM=MAX(1,KPC)

  if( process_id == master_id )THEN
    CALL SEEK('C82')
    READ(1,*,ERR=10,END=30)ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
  end if

  Call Broadcast_Scalar(ISLSHA, master_id)
  Call Broadcast_Scalar(MLLSHA, master_id)
  Call Broadcast_Scalar(NTCLSHA,master_id)
  Call Broadcast_Scalar(ISLSTR, master_id)
  Call Broadcast_Scalar(ISHTA,  master_id)

  MLM=MAX(1,MLLSHA)

  if( process_id == master_id )THEN
    CALL SEEK('C84')
    READ(1,*,ERR=10,END=30)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR
  end if

  Call Broadcast_Scalar(ISTMSR, master_id)
  Call Broadcast_Scalar(MLTMSR, master_id)
  Call Broadcast_Scalar(NBTMSR, master_id)
  Call Broadcast_Scalar(NSTMSR, master_id)
  Call Broadcast_Scalar(NWTMSR, master_id)
  Call Broadcast_Scalar(NTSSTSP,master_id)
  Call Broadcast_Scalar(TCTMSR, master_id)

  MLTMSRM=MAX(1,MLTMSR)
  NTSSTSPM=MAX(1,NTSSTSP)
  MTSSTSPM=1
  if( process_id == master_id )THEN
    IF( NTSSTSP > 0 )THEN
      CALL SEEK('C85')
      DO ITSSS=1,NTSSTSP
        READ(1,*,ERR=10,END=30)I,M
        MTSSTSPM=MAX(MTSSTSPM,M)
      ENDDO
    ENDIF

  end if

  Call Broadcast_Scalar(MTSSTSPM, master_id)

  ! ***    
  if( process_id == master_id )THEN
    CALL SEEK('C91')
    !                                       NOT      NOT
    READ(1,*,ERR=10,END=30) NCDFOUT,DEFLEV,ROTA,BLK,UTMZ,HREST,BASEDATE,BASETIME,PROJ
    ! *** close efdc.inp
    CLOSE(1)
  end if
  
  Call Broadcast_Scalar(NCDFOUT, master_id)
  Call Broadcast_Scalar(DEFLEV,master_id)
  Call Broadcast_Scalar(ROTA,master_id)
  Call Broadcast_Scalar(BLK,master_id)
  Call Broadcast_Scalar(UTMZ,master_id)
  Call Broadcast_Scalar(HREST,master_id)
  Call Broadcast_Scalar(BASEDATE,master_id)
  Call Broadcast_Scalar(BASETIME,master_id)
  Call Broadcast_Scalar(PROJ,master_id)
  ! *** END OF EFDC.INP SCANNING


  IF( ISVEG >= 1 )THEN

    if( process_id == master_id )THEN
      CALL OPENFILE('VEGE.INP')
      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*,ERR=10,END=30) MVEGTYP, MVEGOW, NVEGSER, UVEGSCL
    endif

    Call Broadcast_Scalar(MVEGTYP,master_id)
    Call Broadcast_Scalar(MVEGOW, master_id)
    Call Broadcast_Scalar(NVEGSER,master_id)
    Call Broadcast_Scalar(UVEGSCL,master_id)

    NVEGTPM = MAX(NVEGTPM,MVEGTYP)
    NVEGSERM = MAX(NVEGSERM,NVEGSER)

    if( process_id == master_id )THEN
      CLOSE(1)
    endif

    IF( NVEGSER >= 1 )THEN
      if( process_id == master_id )THEN
        CALL OPENFILE('VEGSER.INP')
        STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      endif

      DO NS=1,NVEGSER
        if( process_id == master_id )THEN
          READ(1,*,ERR=10,END=30) M1,TC1,TAV1
        endif

        Call Broadcast_Scalar(M1, master_id)

        NDVEGSER=MAX(NDVEGSER,M1)

        if( process_id == master_id )THEN
          DO M=1,M1
            READ(1,*)
          ENDDO
        endif
      ENDDO
      if( process_id == master_id )THEN
        CLOSE(1)
      endif

    ENDIF

    ! *** DETERMINE IF MHK IS USED
    LMHK=.FALSE.
    if( process_id == master_id )THEN
      CALL OPENFILE('DXDY.INP')
      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
    endif

    write(mpi_log_unit,'(a,i6)') 'LVC: ', LVC

    Call Broadcast_Scalar(LVC, master_id)

    TCOUNT=0
    DO I=1,LVC
      if( process_id == master_id )THEN
        READ(1,*)IITMP,IITMP,R1TMP,R2TMP,R3TMP,R4TMP,R5TMP,IJTMP
      endif

      Call Broadcast_Scalar(IJTMP, master_id)

      IF( IJTMP > 90 )THEN
        LMHK=.TRUE.
        TCOUNT=TCOUNT+1
      ENDIF
    ENDDO

    if( process_id == master_id )THEN
      CLOSE(1)
    endif

    IF( LMHK )THEN
      if( process_id == master_id )THEN
        CALL OPENFILE('MHK.INP')
        CALL SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES
        READ(1,*,ERR=10,END=30)MHKTYP
        CLOSE(1)
      endif

      Call Broadcast_Scalar(MHKTYP, master_id)

    ENDIF
  ENDIF

  ! *** BANK EROSION
  IF( ISBKERO >= 1 )THEN
    if( process_id == master_id )THEN
      CALL OPENFILE('BEMAP.INP')
      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*,ERR=10,END=30)NBEPAIR,NBESER
    endif

    Call Broadcast_Scalar(NBEPAIR, master_id)
    Call Broadcast_Scalar(NBESER,  master_id)

    NBEPAIRM=NBEPAIR
    NBESERM=NBESER

    if( process_id == master_id )THEN
      CLOSE(1)
    end if

    NDBESER=1
    IF( NBESER > 0 )THEN
      if( process_id == master_id )THEN
        CALL OPENFILE('BESER.INP')
      end if

      DO NS=1,NBESER
        if( process_id == master_id )THEN
100       READ(1,*,ERR=100,END=30)M,R,R,R,R
        endif

        Call Broadcast_Scalar(M, master_id)

        NDBESER=MAX(NDBESER,M)
        DO I=1,M
          if( process_id == master_id )THEN
            READ(1,*,ERR=10,END=30)R,R,R
          endif
        ENDDO
      ENDDO

      if( process_id == master_id )THEN
        CLOSE(1)
      endif

    ENDIF
  ENDIF

  ! *** TOXIC SIMULATION FILES
  NPMXPTSM=1
  NPMXZM=1
  IF( ISTRAN(5) > 0 .AND. NTOX > 0 )THEN
    ! *** SEDIMENT MIXING
    IF( ITMPPMX == 1 )THEN
      CALL OPENFILE('PARTMIX.INP')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*)NPMXZ,NPMXPTS,PMIXSF
      NPMXPTSM=NPMXPTS
      NPMXZM=NPMXZ
      CLOSE(1)
    ENDIF

    ! *** TIME SERIES DRY DEPOSITION BEING USED
    IF( SUM(TOXDEP(:).ITXDRY) > 0 )THEN
      CALL OPENFILE('TXDRY.INP')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*) ITYPE, NDATAPTS
      TXDRYSER(1).NREC = NDATAPTS

      ALLOCATE(TXDRYSER(1).TIM(NDATAPTS), TXDRYSER(1).VAL(NDATAPTS,NTXM))
      CLOSE(1)
    ENDIF

    ! *** TIME SERIES WET DEPOSITION BEING USED
    IF( SUM(TOXDEP(:).ITXWET) > 0 )THEN
      CALL OPENFILE('TXWET.INP')

      STR = READSTR(1)             ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*) ITYPE, NDATAPTS
      TXWETSER(1).NREC = NDATAPTS

      ALLOCATE(TXWETSER(1).TIM(NDATAPTS), TXWETSER(1).VAL(NDATAPTS,NTXM))
      CLOSE(1)
    ENDIF
  ENDIF

  IF( ISCONNECT >= 2 )THEN
    CALL OPENFILE('MAPPGEW.INP')

    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR=READSTR(1)
    READ(1,*,IOSTAT=ISO) NPEWBP
    CLOSE(1)
  ENDIF

  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    CALL OPENFILE('MAPPGNS.INP')

    ! *** SKIP OVER TITLE AND AND HEADER LINES
    STR=READSTR(1)
    READ(1,*,IOSTAT=ISO) NPNSBP ! global value for now
    CLOSE(1)
  ENDIF

  MDCHH = 0
  IF( ISCHAN > 0 )THEN
    if( process_id == master_id )THEN
      CALL OPENFILE('MODCHAN.INP')

      STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      READ(1,*) MDCHH,MDCHHD,MDCHHD2
      CLOSE(1)
    ENDIF
    Call Broadcast_Scalar(MDCHH, master_id)
  ENDIF

  !end if !***End on master process

  IF( NWSER > 1 )THEN
    ! *** ****************************
    if( process_id == master_id )THEN
      CALL OPENFILE('WNDMAP.INP')
      STR=READSTR(1)
      READ(1,*,IOSTAT=ISO) NWNDMAP
      IF( ISO /= 0) CALL STOPP(' *** WNDMAP.INP: READING ERROR!')
      CLOSE(1)
    end if!***End on master process

    Call Broadcast_Scalar(NWSER, master_id)
    Call Broadcast_Scalar(NWNDMAP, master_id)

    ALLOCATE(TWNDMAPBEG(NWNDMAP),TWNDMAPEND(NWNDMAP))
    ALLOCATE(WNDWHT(NWSER,LCM,NWNDMAP))
    ALLOCATE(WNDWHT_TEMP(NWSER))

  ENDIF

  Call Broadcast_Scalar(NASER, master_id)

  IF( NASER > 1 )THEN
    ! *** ****************************
    if( process_id == master_id )THEN
      CALL OPENFILE('ATMMAP.INP')
      STR=READSTR(1)
      READ(1,*,IOSTAT=ISO) NATMMAP
      IF( ISO /= 0) CALL STOPP(' *** ATMMAP.INP: READING ERROR!')
      CLOSE(1)
    end if!***End on master process

    Call Broadcast_Scalar(NATMMAP, master_id)

    ALLOCATE(TATMMAPBEG(NATMMAP),TATMMAPEND(NATMMAP))
    ALLOCATE(ATMWHT(NASER,LCM,NATMMAP))
    ALLOCATE(ATMWHT_TEMP(NASER))     ! Temporary storage for remapping in input.f90 because of domain decomposition
  ENDIF

  IF( ISICE == 1 .AND. NISER > 1 )THEN
    ! *** ****************************
    if( process_id == master_id )THEN
      CALL OPENFILE('ICEMAP.INP')

      STR=READSTR(1)
      READ(1,*,IOSTAT=ISO) NICEMAP
      IF( ISO /= 0) CALL STOPP(' *** ICEMAP.INP: READING ERROR!')
      CLOSE(1)
    end if!***End on master process

    Call Broadcast_Scalar(NICEMAP, master_id)
    Call Broadcast_Scalar(NISER,   master_id)

    ALLOCATE(TICEMAPBEG(NICEMAP),TICEMAPEND(NICEMAP))
    ALLOCATE(RICEWHT(NICEMAP,LCM,NISER))
    ALLOCATE(RICEWHT_Global(NICEMAP,LCM_Global,NISER))
  ENDIF

  RETURN

  ! *** ERROR MESSAGES
10 WRITE(6,'(A)') '  READ ERROR '//CARDNO//' IN INPUT FILE: '//TRIM(CFILE)
  WRITE(8,'(A)') '  READ ERROR '//CARDNO//' IN INPUT FILE: '//TRIM(CFILE)
  CALL STOPP('.')

30 WRITE(6,'(A)') '  READ ERROR '//CARDNO//': UNEXPECTED END OF INPUT FILE: '//TRIM(CFILE)
  WRITE(8,'(A)') '  READ ERROR '//CARDNO//': UNEXPECTED END OF INPUT FILE: '//TRIM(CFILE)
  CALL STOPP('.')

END SUBROUTINE

SUBROUTINE SCANMODC
  INTEGER :: M,I
  REAL    :: R

  CALL OPENFILE('MODCHAN.INP')

10 READ(1,*,ERR=10,END=40)M,I,I
  NCHANM=MAX(1,M)
  READ(1,*,ERR=20,END=40)I,I,R
  CLOSE(1)
  RETURN

20 WRITE(6,'(A)') '  READ ERROR IN INPUT FILE: MODCHAN.INP'
  WRITE(8,'(A)') '  READ ERROR IN INPUT FILE: MODCHAN.INP'
  CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: MODCHAN.INP'
  WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: MODCHAN.INP'
  CALL STOPP('.')

END SUBROUTINE

SUBROUTINE SCANGWSR
  INTEGER :: I,J,M,NS
  REAL    :: R,T,F

  if( process_id == master_id )THEN
    CALL OPENFILE('GWSER.INP')

  10 READ(1,*,ERR=10,END=40) NGWSER
    NGWSERM = MAX(1,NGWSER)
    DO NS=1,NGWSER
      READ(1,*,ERR=20,END=40)M,R,R,R,R,I
      NDGWSER = MAX(NDGWSER,M)
      DO I=1,M
        READ(1,*,ERR=20,END=40)T,F,(R,J=1,3+NDYE+NSED+NSND+NTOX)   ! DELME - WQ
      ENDDO
    ENDDO
    CLOSE(1)
  endif
  
  Call Broadcast_Scalar(NGWSERM, master_id)
  Call Broadcast_Scalar(NGWSER, master_id)
  Call Broadcast_Scalar(NDGWSER, master_id)
  
  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: GWSER.INP IN SERIES:',NS,', POINT:',I
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: GWSER.INP IN SERIES:',NS,', POINT:',I
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: GWSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: GWSER.INP'
   CALL STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANASER
  USE INFOMOD,ONLY:SKIPCOM
  INTEGER :: M,I,NS
  REAL    :: R
  CHARACTER*120 LIN,STR*200

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('ASER.INP')
  end if!***End on master process

  Call Broadcast_Scalar(NASER, master_id)

  ALLOCATE(TSATM(NASER))

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL SKIPCOM(1,'*')
    READ(1,'(A)') STR
    STR=READSTR(1)
  endif

  DO NS=1,NASER

    ! *** ****************************
    if( process_id == master_id)then
      READ(1,*,END=40)M,R,R,I,R,R,R,R
    end if!***End on master process

    Call Broadcast_Scalar(M, master_id)

    NDASER=MAX(NDASER,M)
    ALLOCATE(TSATM(NS).TIM(M),TSATM(NS).VAL(M,7))
    TSATM(NS).NREC= M
    TSATM(NS).TIM = 0
    TSATM(NS).VAL = 0

    ! *** ****************************
    if( process_id == master_id )THEN
      DO I=1,M
        READ(1,*,ERR=20,END=40)R,R,R,R,R,R,R,R
      ENDDO
    endif!***End on master process

  ENDDO
  ! *** ****************************
  if( process_id == master_id )THEN
    CLOSE(1)
  End if!***End on master process

  ! *** ****************************
  if( process_id == master_id )THEN
    IF( ISTRAN(8) > 0 )THEN
      IF( IWQSUN == 1 )THEN
        CALL OPENFILE('SUNDAY.INP')

        M=0
        CALL SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES
        !DO I = 1,7
        !  READ(1,'(A120)',ERR=30,END=40)LIN
        !END DO
        READ(1,*,ERR=30,END=40)M,R,R,R,R
        CLOSE(1)
        NDASER=MAX(NDASER,M)
      ENDIF
    ENDIF
  end if!***End on master process

  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: ASER.INP IN SERIES:',NS,', POINT:',I
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: ASER.INP IN SERIES:',NS,', POINT:',I
   CALL STOPP('.')
  
30 WRITE(6,'(A)') '  READ ERROR IN INPUT FILE'
   WRITE(8,'(A)') '  READ ERROR IN INPUT FILE'
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE'
   CALL STOPP('.')

END SUBROUTINE

SUBROUTINE SCANSSER
  INTEGER :: NS,I,J,M,K
  REAL   :: R

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('SSER.INP')
  endif!***End on master process

  DO NS=1,NSER(1)
    ! *** ****************************
    if( process_id == master_id )THEN
10    READ(1,*,ERR=10,END=40)I,M,R,R,R,R
    endif!***End on master process

    Call Broadcast_Scalar(M, master_id)

    NDCSER=MAX(NDCSER,M)

    ALLOCATE(TSSAL(NS).VAL(M,KCM),TSSAL(NS).TIM(M))
    TSSAL(NS).NREC=M
    TSSAL(NS).VAL(:,:)=0
    TSSAL(NS).TIM(:)=0

    ! *** ****************************
    if( process_id == master_id )THEN
      IF( I == 1 )THEN
        READ(1,*,ERR=20,END=40)(R,K=1,KC)
        DO J=1,M
          READ(1,*,ERR=20,END=40)R,R
        ENDDO
      ELSE
        DO J=1,M
          READ(1,*,ERR=20,END=40)R,(R,K=1,KC)
        ENDDO
      ENDIF
    endif!***End on master process

  ENDDO
  ! *** ****************************
  if( process_id == master_id )THEN
    CLOSE(1)
  endif!***End on master process

  RETURN
  
   ! *** ****************************
20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SSER.INP IN SERIES:',NS,', POINT:',J
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SSER.INP IN SERIES:',NS,', POINT:',J
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: SSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: SSER.INP'
   CALL STOPP('.')

END SUBROUTINE

SUBROUTINE SCANTSER
  INTEGER   :: NS,I,J,K,M,NN
  REAL      :: R
  CHARACTER(200) :: STR

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('TSER.INP')
    STR=READSTR(1)
  endif

  Call Broadcast_Scalar(NSER(2), master_id)

  DO NS=1,NSER(2)
    ! *** ****************************
    if( process_id == master_id )THEN
      READ(1,*,ERR=991) I,M,R,R,R,R
    endif!***End on master process

    Call Broadcast_Scalar(M, master_id)

    NDCSER=MAX(NDCSER,M)

    ALLOCATE(TSTEM(NS).VAL(M,KCM),TSTEM(NS).TIM(M))
    TSTEM(NS).NREC=M
    TSTEM(NS).VAL(:,:)=0
    TSTEM(NS).TIM(:)=0
    ! *** ****************************
    if( process_id == master_id )THEN
      IF( I == 1 )THEN
        READ(1,*,ERR=991) (R,K=1,KC)
        DO J=1,M
          READ(1,*,ERR=991) R,R
        ENDDO
      ELSE
        DO J=1,M
          READ(1,*,ERR=991) R,(R,K=1,KC)
        ENDDO
      ENDIF
    endif!***End on master process

  ENDDO
  ! *** ****************************
  if( process_id == master_id )THEN
    CLOSE(1)
  endif!***End on master process

  IF( ISICE == 1 .AND. NISER >= 1 )THEN

    ! *** ****************************
    if( process_id == master_id )THEN
      CALL OPENFILE('ISER.INP')
      STR=READSTR(1)
    endif!***End on master process

    ALLOCATE(TSICE(NISER))


    DO NS=1,NISER
      ! *** ****************************
      if( process_id == master_id )THEN
        READ(1,*,ERR=992) M,R,R,R
      endif

      Call Broadcast_Scalar(M, master_id)

      ALLOCATE(TSICE(NS).VAL(M,2),TSICE(NS).TIM(M))  ! ** VAL(M,1): ICE THICKNESS; VAL(M,2): ICE COVER
      TSICE(NS).NREC= M
      TSICE(NS).VAL = 0
      TSICE(NS).TIM = 0
      ! *** ****************************
      if( process_id == master_id )THEN
        DO J=1,M
          READ(1,*,ERR=992)
        ENDDO
      endif!***End on master process

    ENDDO

    ! *** ****************************
    if( process_id == master_id )THEN
      CLOSE(1)
    endif!***End on master process

  ELSEIF( ISICE == 2 )THEN
    ! *** ****************************
    if( process_id == master_id )THEN
      CALL OPENFILE('ISTAT.INP')
      STR=READSTR(1)
      READ(1,*,ERR=993) M,R,R   !MISER(NN),TCISER(NN),TAISER(NN)
    endif!***End on master process

    Call Broadcast_Scalar(M, master_id)
    NN=1

    ALLOCATE(TSICE(NN))
    ALLOCATE(TSICE(NN).VAL(M,2),TSICE(NN).TIM(M))  ! ** VAL(M,1): ICE THICKNESS; VAL(M,2): ICE COVER
    TSICE(NN).NREC= M
    TSICE(NN).VAL = 0
    TSICE(NN).TIM = 0

    ! *** ****************************
    if( process_id == master_id )THEN
      CLOSE(1)
    endif!***End on master process

  ENDIF

  RETURN
991 CALL STOPP('TSER.INP: READING ERROR')
992 CALL STOPP('ISER.INP: READING ERROR')
993 CALL STOPP('ISTAT.INP: READING ERROR')
END SUBROUTINE

SUBROUTINE SCANDSER
  INTEGER :: NS,I,K,M,MD,ISTYP
  REAL   :: R
  CHARACTER*120 SKIP
  CHARACTER(*) :: STR*200

  CALL OPENFILE('DSER.INP')

  ! *** SKIP HEADER
  STR=READSTR(1)

  DO NS=1,NSER(3)
    READ(1,*,ERR=20,END=40)ISTYP,M
    NDCSER=MAX(NDCSER,M)

    DO MD=1,NDYE
      ALLOCATE(TSDYE(NS,MD).VAL(M,KCM),TSDYE(NS,MD).TIM(M))
      TSDYE(NS,MD).NREC=M
      TSDYE(NS,MD).VAL(:,:)=0.
      TSDYE(NS,MD).TIM(:)=0.
    ENDDO

    IF( ISTYP == 1 )THEN
      READ(1,*,ERR=20,END=40)(R,K=1,KC)
    ENDIF

    ! *** SKIP TO THE NEXT SERIES
    DO I=1,M
      DO MD=1,NDYE
        READ(1,'(A)') SKIP
      ENDDO
    ENDDO
  ENDDO
  CLOSE(1)
  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: DSER.INP IN SERIES:',NS,', POINT:',I
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: DSER.INP IN SERIES:',NS,', POINT:',I
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: DSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: DSER.INP'
   CALL STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANSFSR
  INTEGER :: NS,I,J,K,M
  REAL   :: R

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('SFSER.INP')
  endif

  DO NS=1,NSER(4)
    ! *** ****************************
    if( process_id == master_id )THEN
10    READ(1,*,ERR=10,END=40)I,M,R,R,R,R
    endif
    Call Broadcast_Scalar(M, master_id)

    NDCSER=MAX(NDCSER,M)

    ALLOCATE(TSSFL(NS).VAL(M,KCM),TSSFL(NS).TIM(M))
    TSSFL(NS).NREC=M
    TSSFL(NS).VAL(:,:)=0
    TSSFL(NS).TIM(:)=0
    ! *** ****************************
    if( process_id == master_id )THEN
      IF( I == 1 )THEN
        READ(1,*,ERR=20,END=40)(R,K=1,KC)
        DO J=1,M
          READ(1,*,ERR=20,END=40)R,R
        ENDDO
      ELSE
        DO J=1,M
          READ(1,*,ERR=20,END=40)R,(R,K=1,KC)
        ENDDO
      ENDIF
    endif!***End on master process

  ENDDO
  ! *** ****************************
  if( process_id == master_id )THEN
    CLOSE(1)
  endif!***End on master process

  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SFSER.INP IN SERIES:',NS,', POINT:',J
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: SFSER.INP IN SERIES:',NS,', POINT:',J
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: SFSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: SFSER.INP'
   CALL STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANQSER
  INTEGER :: NS,I,J,M,K
  REAL    :: R,T,Q

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('QSER.INP')
  endif

  ALLOCATE(QSER(NQSER))

  DO NS=1,NQSER
    ! *** ****************************
    if( process_id == master_id )THEN
10    READ(1,*,ERR=10,END=40)I,M,R,R,R,R,J
    endif

    Call Broadcast_Scalar(M, master_id)

    ALLOCATE(QSER(NS).VAL(M,KCM),QSER(NS).TIM(M))
    QSER(NS).NREC = M
    QSER(NS).VAL(:,:) = 0.0
    QSER(NS).TIM(:) = 0.0
    ! *** ****************************
    if( process_id == master_id )THEN
      !NDQSER=MAX(NDQSER,M)
      IF( I == 1 )THEN
        READ(1,*,ERR=20,END=40)(R,K=1,KC)
        DO J=1,M
          READ(1,*,ERR=20,END=40)T,Q
        ENDDO
      ELSE
        DO J=1,M
          READ(1,*,ERR=20,END=40)T,(Q,K=1,KC)
        ENDDO
      ENDIF
    endif!***End on master process

  ENDDO
  ! *** ****************************
  if( process_id == master_id )THEN
    CLOSE(1)
  endif!***End on master process

  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QSER.INP IN SERIES:',NS,', POINT:',J
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QSER.INP IN SERIES:',NS,', POINT:',J
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: QSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: QSER.INP'
   CALL STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANQWSER
  Use fson
  Use mod_fson_value, Only: fson_value_count, fson_value_get
  INTEGER :: NTMP,I,J,M,NV,NS
  REAL(RKD)   :: R,T,Q
  Type(fson_value), Pointer :: json_data
  
  ! *** ****************************
  if( process_id == master_id )THEN

    ! *** SETTING NUMBER OF CONSTITUENTS (NTMP)
    NTMP = 3 + NDYM + NSED + NSND + NTOX
    ! *** Handle Water Quality variables, if needed
    IF( ISTRAN(8) > 0 )THEN
      !CALL OPENFILE('WQ_3DWC.INP')
      !
      !CALL SEEK('C02')
      !READ(1,*) ISWQLVL,NWQV,NWQZ,NWQPS,NWQTD,NWQTS,NTSWQV,NSMG,NSMZ,NTSSMV,NSMTS,WQKINUPT
      !CLOSE(1)
      !
      !NTMP = NTMP + NWQV   ! *** NTMP INCLUDES WQ
      json_data => fson_parse("wq_3dwc.jnp")
    
      Call fson_get(json_data, "number_of_algae_groups",                          NALGAE)
      Call fson_get(json_data, "number_of_zooplankton_groups",                    NZOOPL)
      NTMP = NTMP + 19 + NALGAE + NZOOPL
    ENDIF

    CALL OPENFILE('QWRS.INP')

    DO NS=1,NQWRSR
10    READ(1,*,ERR=10,END=40)I,M,R,R,R,R
      NDQWRSR = MAX(NDQWRSR,M,1)
      IF( I == 0 )THEN
        ! *** Flow Only
        DO J=1,M
          READ(1,*,ERR=20,END=40)T,Q
        ENDDO
      ELSE
        ! *** Flow with Rise/Fall
        DO J=1,M
          READ(1,*,ERR=20,END=40)T,Q,(R,NV=1,NTMP)
        ENDDO
      ENDIF
    ENDDO

    CLOSE(1)
  End if !*** Writing on master

  Call Broadcast_Scalar(NDQWRSR, master_id)
  
  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QWRS.INP IN SERIES:',NS,', POINT:',J
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: QWRS.INP IN SERIES:',NS,', POINT:',J

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: QWRS.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: QWRS.INP'
   CALL STOPP('.')

END SUBROUTINE

SUBROUTINE SCANPSER
  INTEGER :: NS,M,I,I1
  REAL   :: R,T,E

  ! *** ****************************
  if( process_id == master_id )THEN

    CALL OPENFILE('PSER.INP')
    DO NS=1,NPSER
10    READ(1,*,ERR=10,END=40)I1,M,R,R,R,R
      NDPSER=MAX(NDPSER,M)
      IF( I1 == 1)READ(1,*)R,R,R
      DO I=1,M
        IF( I1 == 0 )THEN
          READ(1,*,ERR=20,END=40)T,E
        ELSEIF( I1 == 1 )THEN
          READ(1,*,ERR=20,END=40)T,E,R
        ENDIF
      ENDDO
    ENDDO
    CLOSE(1)
  end if!***End on master process

  Call Broadcast_Scalar(NDPSER, master_id)

  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: PSER.INP IN SERIES:',NS,', POINT:',I
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: PSER.INP IN SERIES:',NS,', POINT:',I
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: PSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: PSER.INP'
   CALL STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANWSER
  INTEGER :: I,M,NS
  REAL    :: R

  ! *** ****************************

  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('WSER.INP')
  endif!***End on master process

  ALLOCATE(TSWND(NWSER))

  DO NS=1,NWSER
    ! *** ****************************
    if( process_id == master_id )THEN
10    READ(1,*,ERR=10,END=40)M,R,R,R,I
    endif

    Call Broadcast_Scalar(M, master_id)

    ALLOCATE(TSWND(NS).TIM(M),TSWND(NS).VAL(M,2))
    TSWND(NS).NREC= M
    TSWND(NS).VAL = 0
    TSWND(NS).TIM = 0

    ! *** ****************************
    if( process_id == master_id )THEN
      DO I=1,M
        READ(1,*,ERR=20,END=40)R,R,R
      ENDDO
    endif!***End on master process

  ENDDO

  ! *** ****************************
  if( process_id == master_id )THEN
    CLOSE(1)
  endif!***End on master process

  RETURN

20 WRITE(6,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: WSER.INP IN SERIES:',NS,', POINT:',I
   WRITE(8,'(A,I5,A,I10)') '  READ ERROR IN INPUT FILE: WSER.INP IN SERIES:',NS,', POINT:',I
   CALL STOPP('.')

40 WRITE(6,'(A)') '  UNEXPECTED END OF INPUT FILE: WSER.INP'
   WRITE(8,'(A)') '  UNEXPECTED END OF INPUT FILE: WSER.INP'
   CALL STOPP('.')
   
END SUBROUTINE

SUBROUTINE SCANQCTL
  CHARACTER*80 STR*200
  CHARACTER*80 :: SKIP

  INTEGER :: M, M1, M2, MP, IS, NS, ISO, ISTYP
  REAL    :: HUA, HUM, HDA, HDM, R, A, A1

  if( process_id == master_id )then
  
    CALL OPENFILE('QCTL.INP')

    ! *** FIND THE MAXIMUM NUMBER OF TABLE DATA POINTS
    NDQCLT=0
    STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
    DO NS=1,NQCTLT
      READ(1,*,IOSTAT=ISO) ISTYP, MP, HUA, HUM, HDA, HDM, R, A, A1
      NDQCLT = MAX(NDQCLT,MP)
      IF( ISO > 0 )GOTO 20

      IF( ISTYP == 0 )THEN
        DO M=1,MP
          READ(1,'(A)')SKIP
        ENDDO
      ELSEIF( ISTYP == 1 )THEN
        READ(1,'(A)')SKIP
        DO M=1,MP
          READ(1,'(A)')SKIP
        ENDDO
      ELSEIF( ISTYP == 2 )THEN
        DO M1=1,MP
          DO M2=1,MP
            READ(1,'(A)')SKIP
          ENDDO
        ENDDO
      ELSEIF( ISTYP == 3 )THEN
        READ(1,'(A)')SKIP
        DO M1=1,MP
          DO M2=1,MP
            READ(1,'(A)')SKIP
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    CLOSE(1)
  endif
  Call Broadcast_Scalar(NDQCLT, master_id)
  
  NDQCLT2 = NDQCLT

  RETURN

20 WRITE(6,'(A)') ' READ ERROR IN FILE: '//TRIM(CFILE)
   WRITE(8,'(A)') ' READ ERROR IN FILE: '//TRIM(CFILE)
   CALL STOPP('.')

END SUBROUTINE

SUBROUTINE SCANWQ

  Use SHELLFISHMOD
  Use WQ_RPEM_MODULE
  Use fson
  Use mod_fson_value, Only: fson_value_count, fson_value_get
  Use INFOMOD, Only:SKIPCOM 
  
  Character(len=2)   :: SNUM 
  Character(len=15)  :: FNWQSR(40)
  Character(len=15)  :: FNWQSRALGAE(40)
  Character(len=100) :: ERRSTR
  Character(len=120) :: LINE
  
  Integer :: NW, NPS, I, J, K, IDRAG, ITMP, M, NS, ISTYP, NAL, NDWQCSR
  Real    :: XPSQ, TM, TA, RMULADJ, ADDADJ, T1, T2, TSMTSB, TSMTSE, SMTSDT
  LOGICAL :: BFLAG

  Type(fson_value), Pointer :: json_data, item, pointsource, algaegroups

  ! *** ****************************
  if( process_id == master_id )THEN 
    write(*,'(A)') 'SCANNING EUTROPHICATION CONTROL FILE'
    
    json_data => fson_parse("wq_3dwc.jnp")
    
    Call fson_get(json_data, "kinetics_option",                                  ISWQLVL)
    Call fson_get(json_data, "number_of_kinetic_zones",                          NWQZ)
    Call fson_get(json_data, "temperature_lookup_table_size",                    NWQTD)
    Call fson_get(json_data, "number_of_time_series_output_locations",           NWQTS)
    Call fson_get(json_data, "number_of_time_series_output_variables",           NTSWQV)
    Call fson_get(json_data, "number_of_sediment_zones",                         NSMZ)
    Call fson_get(json_data, "number_of_sediment_time_series_output_variables",  NTSSMV)
    Call fson_get(json_data, "max_number_of_time_series_output_locations",       NSMTS)

    Call fson_get(json_data, "point_source_load_option",                         IWQPSL)
    Call fson_get(json_data, "mass_loading_point_sources.number_of_point_sources",NWQPS) 
    Call fson_get(json_data, "mass_loading_point_sources.number_of_time_series", NPSTMSR) 
    
    Call fson_get(json_data, "number_of_algae_groups",                           NALGAE)
                                                                                 
    Call fson_get(json_data, "zooplankton_activate",                             IWQZPL)
    Call fson_get(json_data, "number_of_zooplankton_groups",                     NZOOPL)
                                                                                 
    Call fson_get(json_data, "shellfish_farm_activate",                          ISFFARM)
    Call fson_get(json_data, "number_of_shellfish_species",                      NSF)
    
    ! *** Update number of WQ components
    NWQV = 19 + NALGAE + NZOOPL
    
    ! *** Set constituent index based on kinetic option
    IF( ISWQLVL == 0 )THEN
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
    ELSE
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
    ENDIF
    
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
    
    IF( ISWQLVL == 0 )THEN
      ! *** Algae names if WQSKE0
      WQCONSTIT(ICHC) = 'CHC'
      WQCONSTIT(ICHD) = 'CHD'
      WQCONSTIT(ICHG) = 'CHG'
    ELSE
    ! *** Algae names if WQSKE1
      DO NW = 1,NALGAE
        WRITE(SNUM,'(I1)') NW
        WQCONSTIT(19+NW) = 'ALG'//SNUM
      ENDDO
    ENDIF
    ! *** Zooplankton name
    IF( IWQZPL > 0 )THEN
      DO NW = 1,NZOOPL
        WRITE(SNUM,'(I1)') NW
        WQCONSTIT(19+NALGAE+NW) = 'ZOO'//SNUM
      ENDDO
    ENDIF
    
    Call fson_get(json_data, "active_constituents.ROC",                          ISKINETICS(IROC))
    Call fson_get(json_data, "active_constituents.LOC",                          ISKINETICS(ILOC))
    Call fson_get(json_data, "active_constituents.DOC",                          ISKINETICS(IDOC))
    Call fson_get(json_data, "active_constituents.ROP",                          ISKINETICS(IROP))
    Call fson_get(json_data, "active_constituents.LOP",                          ISKINETICS(ILOP))
    Call fson_get(json_data, "active_constituents.DOP",                          ISKINETICS(IDOP))
    Call fson_get(json_data, "active_constituents.P4D",                          ISKINETICS(IP4D))
    Call fson_get(json_data, "active_constituents.RON",                          ISKINETICS(IRON))
    Call fson_get(json_data, "active_constituents.LON",                          ISKINETICS(ILON))
    Call fson_get(json_data, "active_constituents.DON",                          ISKINETICS(IDON))
    Call fson_get(json_data, "active_constituents.NHX",                          ISKINETICS(INHX))
    Call fson_get(json_data, "active_constituents.NOX",                          ISKINETICS(INOX))
    Call fson_get(json_data, "active_constituents.SUU",                          ISKINETICS(ISUU))
    Call fson_get(json_data, "active_constituents.SAA",                          ISKINETICS(ISAA))
    Call fson_get(json_data, "active_constituents.COD",                          ISKINETICS(ICOD))
    Call fson_get(json_data, "active_constituents.DOX",                          ISKINETICS(IDOX))
    Call fson_get(json_data, "active_constituents.TAM",                          ISKINETICS(ITAM))
    Call fson_get(json_data, "active_constituents.FCB",                          ISKINETICS(IFCB))
    Call fson_get(json_data, "active_constituents.CO2",                          ISKINETICS(ICO2))
    If( ISWQLVL == 0  )then                                                      
      Call fson_get(json_data, "active_constituents.CHC",                        ISKINETICS(ICHC))
      Call fson_get(json_data, "active_constituents.CHD",                        ISKINETICS(ICHD))
      Call fson_get(json_data, "active_constituents.CHG",                        ISKINETICS(ICHG))
    End if
    
    Call fson_get(json_data, "sediment_diagenesis.benthic_flux_option",          IWQBEN)
    Call fson_get(json_data, "sediment_diagenesis.number_of_reactive_classes",   NSMG)
    Call fson_get(json_data, "rpem_activate",                                    ISRPEM)

    Do NW = 1,NWQV
      Call fson_get(json_data,'number_of_time_series.'//WQCONSTIT(NW),           NWQCSR(NW))
    Enddo
    
    IF( IWQPSL == 2 )THEN
      ! *** Mass loading point source options
      Call fson_get(json_data, "mass_loading_point_sources.constant_point_sources",   pointsource)
      DO NPS = 1, fson_value_count(pointsource)
        item => fson_value_get(pointsource, NPS)
        Call fson_get(item, "I", I)
        Call fson_get(item, "J", J)
        Call fson_get(item, "K", K)
        Call fson_get(item, "NSR", ITMP)
        Call fson_get(item, "PSQ", XPSQ)
        NCSERM = MAX(1,NCSERM,ITMP)
      ENDDO
    ENDIF
  endif   ! *** End on master process
  
  Call Broadcast_Scalar(ISWQLVL, master_id)
  Call Broadcast_Scalar(NWQZ,    master_id)
  Call Broadcast_Scalar(NWQTD,   master_id)
  Call Broadcast_Scalar(NWQTS,   master_id)
  Call Broadcast_Scalar(NTSWQV,  master_id)
  Call Broadcast_Scalar(NSMZ,    master_id)
  Call Broadcast_Scalar(NTSSMV,  master_id)
  Call Broadcast_Scalar(NSMTS,   master_id)

  Call Broadcast_Scalar(IWQPSL,  master_id)  
  Call Broadcast_Scalar(NWQPS,   master_id)
  Call Broadcast_Scalar(NPSTMSR, master_id)
  
  Call Broadcast_Scalar(NALGAE,  master_id)
  Call Broadcast_Scalar(IWQZPL,  master_id)
  Call Broadcast_Scalar(NZOOPL,  master_id)
  
  Call Broadcast_Scalar(ISFFARM, master_id)
  Call Broadcast_Scalar(NSF,     master_id)
  Call Broadcast_Scalar(NWQV,    master_id)
  
  Call Broadcast_Scalar(IROC,    master_id)
  Call Broadcast_Scalar(ILOC,    master_id)
  Call Broadcast_Scalar(IDOC,    master_id)
  Call Broadcast_Scalar(IROP,    master_id)
  Call Broadcast_Scalar(ILOP,    master_id)
  Call Broadcast_Scalar(IDOP,    master_id)
  Call Broadcast_Scalar(IP4D,    master_id)
  Call Broadcast_Scalar(IRON,    master_id)
  Call Broadcast_Scalar(ILON,    master_id)
  Call Broadcast_Scalar(IDON,    master_id)
  Call Broadcast_Scalar(INHX,    master_id)
  Call Broadcast_Scalar(INOX,    master_id)
  Call Broadcast_Scalar(ISUU,    master_id)
  Call Broadcast_Scalar(ISAA,    master_id)
  Call Broadcast_Scalar(ICOD,    master_id)
  Call Broadcast_Scalar(IDOX,    master_id)
  Call Broadcast_Scalar(ITAM,    master_id)
  Call Broadcast_Scalar(IFCB,    master_id)
  Call Broadcast_Scalar(ICO2,    master_id)
  
  Call Broadcast_Scalar(IWQBEN,  master_id)
  Call Broadcast_Scalar(NSMG,    master_id)
  Call Broadcast_Scalar(ISRPEM,  master_id)
  
  Call Broadcast_Scalar(IWQTS,   master_id)
  Call Broadcast_Scalar(WQTSDT,  master_id)
  Call Broadcast_Scalar(ISCOMP,  master_id)

  Call Broadcast_Array(ISKINETICS, master_id)
  Call Broadcast_Array(NWQCSR,     master_id)

  NWQZM   = MAX(1,NWQZ)
  NWQPSM  = MAX(1,NWQPS)
  NWQTDM  = MAX(1,NWQTD) 
  NWQTSM  = MAX(1,NWQTS)
  NSMGM   = MAX(1,NSMG)
  NSMZM   = MAX(1,NSMZ)
  NTSSMVM = MAX(1,NTSSMV)
  NSMTSM  = MAX(1,NSMTS)
  NALGAEM = MAX(1,NALGAE)
  ! NTSWQVM = MAX(1,NTSWQV)
        
  ! *** Update kinetics bypass flag for alage + zooplankton
  IF( NALGAE > 0 .OR. NZOOPL > 0 )then
    Do NW = 20,NWQV
      ISKINETICS(NW) = 1
    Enddo
  End if
  NWQVM = NWQV
    
  !If (ISFFARM > 0) NWQV = NWQV + NSF
    
  NWQPSRM = MAX(1,NPSTMSR)
  DO NW = 1,NWQV
    NWQCSRM = MAX(NWQCSRM,NWQCSR(NW))
  ENDDO
  NCSERM = MAX(NCSERM,NWQCSRM)

  ISMOB = 0
  MACDRAG = 0
  if( process_id == master_id .and. NALGAE > 0 )THEN  ! *** Only read wq_biota when model simulates algae - DKT
    ! *** SCAN WQ_BIOTA.JNP FILE
    write(*,'(A)') 'SCANNING EUTROPHICATION BIOTA CONFIGURATION'

    json_data => fson_parse("wq_biota.jnp")
    CALL fson_get(json_data, "groups", algaegroups)
    DO NAL = 1, fson_value_count(algaegroups)
      item => fson_value_get(algaegroups, NAL)
      CALL fson_get(item, "mobility_flag", ISMOB(NAL)) 

      Call fson_get(item, "use_hydro_feedback", IDRAG)
      IF( IDRAG > 0 .AND. ISVEG == 0 )THEN
        MACDRAG = 1
        ISVEG = 2                                   ! *** Gets overwritten in INPUT but used to flag allocation of vegetation arrays
      ENDIF
    ENDDO
  endif   ! *** End on master process
    
  Call Broadcast_Array( ISMOB,    master_id)
  NFIXED = 0
  
  DO NW = 1, 19
    ISTRWQ(NW) = ISKINETICS(NW)       ! *** Transport flag for nutrients
  ENDDO
  
  DO NAL = 1,NALGAE
    ISKINETICS(NAL+19) = 1            ! *** Kinetics Flag
    ISTRWQ(NAL+19) = 1                ! *** Transport Flag
    IF( ISMOB(NAL) == 0 )THEN 
      NFIXED = NFIXED + 1
      ISTRWQ(NAL+19) = 0              ! *** Do not transport this constituent
    ENDIF
  ENDDO
  Call Broadcast_Scalar(NFIXED,  master_id) 

  if( process_id == master_id )THEN 
    ! *** SCAN THE TIME SERIES
    IF( NPSTMSR >= 1 .AND. IWQPSL /= 2 )THEN
      CALL OPENFILE('WQPSL.INP')
      CALL SKIPCOM(1,'*')
      ERRSTR = 'READING WPQSL.INP'
      DO NS = 1,NPSTMSR
        READ(1,*,ERR=999,END=20) M, TM, TA, RMULADJ, ADDADJ
        NDWQPSR = MAX(NDWQPSR,M)
        DO J = 1,M
          DO K = 1,3
            READ(1,'(A120)')LINE
          ENDDO
        ENDDO
      ENDDO
      20 CLOSE(1)
    ENDIF
  endif   ! *** End on master process

  Call Broadcast_Scalar(NDWQPSR,  master_id) 
    
  ! *** SCAN THE OPEN BC TIME SERIES
  ! **  READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE WQ CONCENTRATION
  ! **  TIME SERIES FROM THE FILES WQCSRNN.INP
  ! *** SCAN FOR NUMBER OF SERIES   
  !IF (ISWQLVL == 0 )THEN
  !  DO NW = 1,40
  !    WRITE(SNUM,'(I2.2)')NW
  !    FNWQSR(NW) ='wqcsr'//SNUM//'.inp' 
  !  ENDDO
  !  WRITE(6,'(A)')'SCANNING INPUT FILE: WQCSRXX.INP'
  !ENDIF
    
  if( process_id == master_id )THEN 
    IF( ISWQLVL == 1  )THEN  
      ! *** Nutrient
      DO NW = 1,19
        WRITE(SNUM,'(I2.2)') NW
        FNWQSR(NW) = 'wqcsr'//SNUM//'.inp'
      ENDDO
      ! *** Algae
      DO NW = 1,NALGAE
        WRITE(SNUM,'(I2.2)')NW
        FNWQSR(NW+19)='wqalgsr'//SNUM//'.inp'
      ENDDO
      ! *** Zooplankton
      DO NW = 1,NZOOPL
         WRITE(SNUM,'(I2.2)')NW
         FNWQSR(NW+19+NALGAE)='wqzoosr'//SNUM//'.inp'
      ENDDO
      
      WRITE(6,'(A)')'SCANNING WQ CONCENTRATION TIME SERIES INPUT FILES'
    ENDIF
        
    NSER(8) = 0
    DO NW = 1,NWQV
      INQUIRE(FILE=FNWQSR(NW),EXIST=BFLAG)
      IF( BFLAG )THEN
        OPEN(1,FILE=FNWQSR(NW))
        
        ERRSTR = 'READING '//TRIM(FNWQSR(NW))
        CALL SKIPCOM(1,'*')
        DO NS = 1,1000
          READ(1,*,ERR=999,END=40)ISTYP,M,T1,T2,RMULADJ,ADDADJ
          IF( ISTYP == 1 ) M = M + 1   ! *** SKIP THE LAYER SPLITS
          DO J = 1,M
            READ(1,'(A120)')LINE
          ENDDO
        ENDDO
        40 CLOSE(1)
        NSER(8) = MAX(NSER(8),NS-1)
      ENDIF
    ENDDO
  endif   ! *** End on master process
    
  Call Broadcast_Scalar(NSER(8),  master_id) 
  
  Allocate(TSWQ(NSER(8),NWQV))

  ! *** SCAN FOR NUMBER OF POINTS IN EACH SERIES
  if( process_id == master_id )THEN
    DO NW=1,NWQV
      ! *** ****************************
      INQUIRE(FILE=FNWQSR(NW),EXIST=BFLAG)
      IF( BFLAG )THEN
        ! *** ****************************
        CALL OPENFILE(FNWQSR(NW))
        CALL SKIPCOM(1,'*')

        DO NS=1,NSER(8)
          ! *** ****************************
          READ(1,*,END=45) ISTYP, M, T1, T2, RMULADJ, ADDADJ
          TSWQ(NS,NW).NREC = M
          
          IF( ISTYP == 1 ) M = M+1   ! *** SKIP THE LAYER SPLITS
          ! *** ****************************
          DO J=1,M
            READ(1,'(A120)')LINE
          ENDDO
        ENDDO
      ENDIF
45    CLOSE(1)
    ENDDO
  endif ! *** End on master process
  
  NDWQCSR = 1
  DO NW=1,NWQV
    DO NS=1,NSER(8)
      Call Broadcast_Scalar(TSWQ(NS,NW).NREC, master_id)
      M = TSWQ(NS,NW).NREC
      NDWQCSR = MAX(NDWQCSR,M)

      ALLOCATE(TSWQ(NS,NW).VAL(M,KCM), TSWQ(NS,NW).TIM(M))
      TSWQ(NS,NW).VAL = 0.
      TSWQ(NS,NW).TIM = 0.
    ENDDO
  ENDDO
  
  NDCSER = MAX(NDCSER,NDWQCSR)
  
  ! *** ****************************
  ! *** Sediment Diagenesis
  IF( IWQBEN == 1 )THEN
    if( process_id == master_id )THEN
      write(*,'(A)') 'SCANNING SEDIMENT DIAGENESIS CONFIGURATION'
      json_data => fson_parse("wq_3dsd.jnp")

      Call fson_get(json_data, "initial_condition_option", ISMICI)
      Call fson_get(json_data, "number_of_spatial_zones", ISMZ)
      Call fson_get(json_data, "write_restart_option", ISMRST)
      Call fson_get(json_data, "benthic_stress.activate_hysteresis_in_benthic_mixing", ISMHYST)
      Call fson_get(json_data, "activate_diagnostic_output", ISMZB)
      
      Call fson_get(json_data, "time_series_output.number_of_locations", ISMTS)
    endif  ! *** End on master process

    Call Broadcast_Scalar(ISMICI,  master_id)
    Call Broadcast_Scalar(ISMZ,    master_id)
    Call Broadcast_Scalar(ISMRST,  master_id)
    Call Broadcast_Scalar(ISMHYST, master_id)
    Call Broadcast_Scalar(ISMZB,   master_id)
    Call Broadcast_Scalar(ISMTS,   master_id)

    NSMTSM = MAX(ISMTS,NSMTS)

  ENDIF

  RETURN

999 PRINT*,ERRSTR
  STOP

    
END SUBROUTINE SCANWQ


SUBROUTINE SCANSEDZLJ
  !  REVISION DATE :  May 24, 2006
  !  Craig Jones and Scott James
  ! *** ************************************************************
  INTEGER :: IDUMMY,ERROR
  ! *** ****************************
  if( process_id == master_id )THEN
    CALL OPENFILE('BED.SDF')
    READ(1,*,IOSTAT=ERROR) !SKIP THIS LINE
    IF( ERROR == 1 )THEN
      WRITE(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      WRITE(8,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    ENDIF
    READ(1,*,IOSTAT=ERROR) IDUMMY, KB, ICALC_BL, SEDSTEP, SEDSTART, IHTSTRT, IMORPH, ISWNWAVE, MAXDEPLIMIT, HPMIN
    IF( ERROR == 1 )THEN
      WRITE(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      WRITE(8,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    ENDIF
    READ(1,*,IOSTAT=ERROR) !SKIP THIS LINE
    IF( ERROR == 1 )THEN
      WRITE(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      WRITE(8,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    ENDIF
    READ(1,*,IOSTAT=ERROR) ITBM,NSICM
    IF( ERROR == 1 )THEN
      WRITE(*,'("READ ERROR IN SEDZLJ INPUT FILE")')
      WRITE(8,'("READ ERROR IN SEDZLJ INPUT FILE")')
      STOP
    ENDIF
    CLOSE(1)
  endif!***End on master process

  Call Broadcast_Scalar(KB,       master_id)
  Call Broadcast_Scalar(ICALC_BL, master_id)
  Call Broadcast_Scalar(ITBM,     master_id)
  Call Broadcast_Scalar(NSICM,    master_id)

END SUBROUTINE SCANSEDZLJ

SUBROUTINE SCNTXSED

  INTEGER :: NCSERNC,NLOOP,IS,NS,ISTYP,NDATAPTS,NN,NT,M
  CHARACTER*120 SKIP
  CHARACTER*10 INFILE
  CHARACTER*80 STR*200


  ! *** NOW FIND MAX FOR TOXICS AND SEDIMENTS
  DO NN=1,3
    NCSERNC = 0
    ! *** ****************************
    if( process_id == master_id )THEN
      IF( NN == 1 )THEN
        IF( NTOX > 0 .AND. NSER(5) > 0 .AND. ISTRAN(5) > 0 )THEN
          CALL OPENFILE('TXSER.INP')
          NLOOP = NTOX
          NCSERNC = NSER(5)
        ENDIF
      ELSEIF( NN == 2 )THEN
        IF( NSED > 0 .AND. NSER(6) > 0 .AND. ISTRAN(6) > 0 )THEN
          CALL OPENFILE('SDSER.INP')
          NLOOP=NSED
          NCSERNC=NSER(6)
        ENDIF
      ELSEIF( NN == 3 )THEN
        IF( NSND > 0 .AND. NSER(7) > 0 .AND. ISTRAN(7) > 0 )THEN
          CALL OPENFILE('SNSER.INP')
          NLOOP=NSND
          NCSERNC=NSER(7)
        ENDIF
      ENDIF
    endif

    Call Broadcast_Scalar(NLOOP, master_id)
    Call Broadcast_Scalar(NCSERNC, master_id)

    IF( NCSERNC > 0 )THEN
      ! *** ****************************
      if( process_id == master_id )THEN
        STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      endif!***End on master process

      ! *** LOOP OVER THE NUMBER OF SERIES
      DO NS=1,NSER(NN+4)  !NCSERNC
        ! *** ****************************
        if( process_id == master_id )THEN
          READ(1,*,ERR=20,END=40)ISTYP,NDATAPTS   !,X1,X2,X3,X4
        endif!***End on master process

        Call Broadcast_Scalar(ISTYP, master_id)
        Call Broadcast_Scalar(NDATAPTS, master_id)

        IF( NN == 1 .AND. NTOX > 0 .AND. NSER(5) > 0 )THEN
          DO NT=1,NTOX
            ALLOCATE(TSTOX(NS,NT).VAL(NDATAPTS,KCM),TSTOX(NS,NT).TIM(NDATAPTS))
            TSTOX(NS,NT).NREC=NDATAPTS
            TSTOX(NS,NT).VAL=0
            TSTOX(NS,NT).TIM=0
          ENDDO

        ELSEIF( NN == 2 .AND. NSED > 0 .AND. NSER(6) > 0 )THEN
          DO NT=1,NSED
            ALLOCATE(TSSED(NS,NT).VAL(NDATAPTS,KCM),TSSED(NS,NT).TIM(NDATAPTS))
            TSSED(NS,NT).NREC=NDATAPTS
            TSSED(NS,NT).VAL=0
            TSSED(NS,NT).TIM=0
          ENDDO
        ELSEIF( NN == 3 .AND. NSND > 0 .AND. NSER(7) > 0 )THEN
          DO NT=1,NSND
            ALLOCATE(TSSND(NS,NT).VAL(NDATAPTS,KCM),TSSND(NS,NT).TIM(NDATAPTS))
            TSSND(NS,NT).NREC=NDATAPTS
            TSSND(NS,NT).VAL=0
            TSSND(NS,NT).TIM=0
          ENDDO
        ENDIF
        ! *** ****************************
        if( process_id == master_id )THEN
          ! *** SKIP THE CONVERSIONS
          IF( NN /= 1 )THEN
            DO NT=2,NLOOP
              READ(1,'(A)') SKIP
            ENDDO
          ENDIF

          NDCSER=MAX(NDCSER,NDATAPTS)
          IF( ISTYP == 1 )THEN
            ! *** SKIP THE SPLITS
            READ(1,'(A)') SKIP
          ENDIF
          ! *** SKIP TO THE NEXT SERIES
          DO M=1,NDATAPTS
            DO NT=1,NLOOP
              READ(1,'(A)') SKIP
            ENDDO
          ENDDO
        endif!***End on master process

      ENDDO
      ! *** ****************************
      if( process_id == master_id )THEN
        CLOSE(1)
      endif!***End on master process

    ENDIF
  ENDDO
  RETURN

   ! *** ****************************
20 WRITE(6,'(A)') '*** READ ERROR IN FILE: '//CFILE
   WRITE(8,'(A)')    '*** READ ERROR IN FILE: '//CFILE
   CALL STOPP('.')

40 WRITE(6,'(A)') '*** UNEXPECTED END OF FILE: '//CFILE
   WRITE(8,'(A)')    '*** UNEXPECTED END OF FILE: '//CFILE
   CALL STOPP('.')

END SUBROUTINE SCNTXSED

SUBROUTINE SCANPROPWASH

  Use fson
  Use mod_fson_value, Only: fson_value_count, fson_value_get

  Type(fson_value), Pointer :: json_data
  
  ! *** ****************************
  IF( process_id == master_id )THEN
    ! *** Scan the propwash config file
    write(*,'(A)') 'SCANNING PROPWASH CONFIGURATION'

    ! *** Open the propwash ship data file
    json_data => fson_parse("propwash_config.jnp")
    
    ! *** Fast Settling Classes
    call fson_get(json_data, "parms.fraction_fast",   fraction_fast)
    call fson_get(json_data, "parms.fast_multiplier", fast_multiplier)
  ENDIF

  Call Broadcast_Scalar(fraction_fast,   master_id)
  Call Broadcast_Scalar(fast_multiplier, master_id)

END SUBROUTINE SCANPROPWASH

SUBROUTINE OPENFILE(INFILE)

  CHARACTER(*), INTENT(IN) :: INFILE
  CFILE = StrUpCase(INFILE)

  WRITE(6,'(A)')'SCANNING INPUT FILE: '//CFILE

  CFILE = StrLowCase(INFILE)
  OPEN(1,FILE=CFILE,STATUS='OLD',ERR=200)

  RETURN

200 WRITE(6,'(A)') '  FILE DOES NOT EXIST:  '//TRIM(CFILE)
  WRITE(8,'(A)') '  FILE DOES NOT EXIST:  '//TRIM(CFILE)
  CALL STOPP('.')

END SUBROUTINE

FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )
  ! FROM "String_Utility" BY Paul van Delst
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  ! -- Local variables
  INTEGER :: i, n

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
    ! -- Find location of letter in lower case constant string
    n = INDEX( LOWER_CASE, Output_String( i:i ) )
    ! -- If current substring is a lower case letter, make it upper case
    IF(  n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
  END DO
END FUNCTION StrUpCase

FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String )
  ! FROM "String_Utility" BY Paul van Delst
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  ! -- Local variables
  INTEGER :: i, n

  ! -- Copy input string
  Output_String = Input_String
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
    ! -- Find location of letter in upper case constant string
    n = INDEX( UPPER_CASE, Output_String( i:i ) )
    ! -- If current substring is an upper case letter, make it lower case
    IF(  n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
  END DO
END FUNCTION StrLowCase

END MODULE
