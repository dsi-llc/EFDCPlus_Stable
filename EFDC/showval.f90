! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE SHOWVAL

  ! *** REWRITTEN BY PAUL M. CRAIG  ON DEC 2006
  ! ***
  ! *** 2010_06 CHANGED THE NSHTYPE TO CORRESPOND TO parameter LIST
       
  use GLOBAL
  use Variables_WQ
  
  use Variables_MPI
  use INFOMOD,only:SKIPCOM,READSTR

  implicit none

  integer :: ISREAD,NSKIP,JSHPRT,INFODT,L,LN,IZSURF,IVELEKC,IVELNKC,IAVKS,IABKS
  integer :: IVELEKB,IVELNKB,IAVKB,IABKB,ICKC,ICKB,NINFO,IQLO(0:5),IQHI(0:5)
  integer :: I, K, NWQ, LQ, LL, NOPEN
  integer, STATIC :: IPARAM,ISUB

  real :: TIME,ZSURF,UTMP,VTMP,VELEKC,VELNKC,VELEKB,VELNKB,AVKS,AVKB,ABKS,ABKB,CKC,CKB
  real :: T1,T2,TSPEED,ETA,QOUT,QIN,QOPEN,ETIME

  character UNITS*3, PARM*4, CSUB*1
  character*80 STR*200
   
  save      INFODT, JSHPRT, UNITS, PARM, NINFO

  real(RKD),external :: DSTIME 

  DATA ISREAD/0/
  DATA UNITS/'PPM'/
  
  integer :: ICSHOW_DUM
  integer :: JCSHOW_DUM
  integer :: i_tmp, j_tmp
  real    :: tmp
  logical :: show_val_inside
  
  if( ISDYNSTP == 0 )then  
    DELT = DT  
  else  
    DELT = DTDYN  
  endif  

  if( ISREAD == 0 )then
    ISREAD = 1
    open(1,FILE = 'show.inp',STATUS = 'OLD')
    STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
    read(1,*)NSHTYPE,NSHOWR,ICSHOW_DUM,JCSHOW_DUM,ISHPRT
    
    ICSHOW = IG2IL(ICSHOW_DUM) !*** Remapping for domain decomposition
    JCSHOW = JG2JL(JCSHOW_DUM) !*** Remapping for domain decomposition
    
    read(1,*)ZSSMIN,ZSSMAX,tmp
    close(1)
    NSHOWR = 20
    NSHOWC = NSHOWR
    if( ISHPRT < 1)ISHPRT = 1
    JSHPRT = ISHPRT
    NINFO = -1
    
    ! *** SET THE DISPLAY parameter
    if( NSHTYPE < 10 )then
      IPARAM = NSHTYPE
      ISUB = 0
      CSUB = '0'
    else
      IPARAM = NSHTYPE/100
      ISUB= MOD(INT(NSHTYPE,4),100)
      if( ISUB < 10)WRITE(CSUB,'(I1)')ISUB
    endif

    UNITS = 'PPM'
    if( IPARAM == 1 )then
      ! *** SALINITY
      UNITS = 'PPT'
      PARM = 'SAL'
    elseif( IPARAM == 2 )then
      ! *** TEMPERATURE
      UNITS = 'D:C'
      PARM = 'TEM'
    elseif( IPARAM == 3 )then
      ! *** DYE
      PARM = 'DYE'
    elseif( IPARAM == 5 )then
      ! *** TOXICS
      UNITS = 'PPB'
      if( ISUB > 0 .and. ISUB <= NTOX )then
        PARM = 'TX' // CSUB
      else
        PARM = 'TOX'
      endif
    elseif( IPARAM == 6 )then
      ! *** COHESIVES
      if( ISUB > 0 .and. ISUB <= NSED )then
        PARM = 'SD' // CSUB
      else
        PARM = 'SED'
      endif
    elseif( IPARAM == 7 )then
      ! *** NON-COHESIVES
      if( ISUB > 0 .and. ISUB <= NSND )then
        PARM = 'SN' // CSUB
      else
        PARM = 'SND'
      endif

    elseif( IPARAM == 0 )then
      ! *** TSS
      PARM = 'TSS'

    elseif( IPARAM == 8 .and. ISTRAN(8) > 0 )then
      ! *** WATER QUALITY
      if( ISUB < 1 .or. ISUB > NWQV ) ISUB = IDOX  ! *** DEFAULT TO D.O.
      PARM = WQCONSTIT(ISUB)

    endif

  endif

  NITERAT = NITERAT + 1
  
  ! *** SETUP THE VARIABLES
  SHOW_VAL_INSIDE = .FALSE.
  I_TMP = IG2IL(ICSHOW)
  J_TMP = JG2JL(JCSHOW)
  
  ! *** CHECK IF THE CELL TO SHOW IS IN THE CURRENT DOMAIN
  if( I_TMP > 0 .and. I_TMP < IC )then
    if( J_TMP > 0 .and. J_TMP < JC )then
      SHOW_VAL_INSIDE = .TRUE.
      L = LIJ(I_TMP,J_TMP)
    endif
  endif
  
  ! *** Set L to a value that definitely will be inside so that we don't index arrays outside of domain
  if(.NOT. SHOW_VAL_INSIDE)then  ! *** L is outside of domain
    L = 2
  endif
  
  TIME = TIMESEC/86400.

  ! *** DISPLAY THE HEADER
  if( NSHOWC >= NSHOWR )then
    NSHOWC = 0
    NINFO = NINFO + 1
    
    if( ISDYNSTP > 0 )then
      INFODT = INFODT+1
      if( INFODT > 10 .or. NITERAT == 1 )then
        if( DTSSDHDT > 0. )then
          write(6,9000) DTL1MN,L1LOC,DTL2MN,L2LOC,DTL3MN,L3LOC,DTL4MN,L4LOC
        else
          write(6,9000) DTL1MN,L1LOC,DTL2MN,L2LOC,DTL3MN,L3LOC
        endif
        write(6,*)' '
        INFODT = 0
      endif
    endif
    
    ! *** MODEL DETAILS
    if( NINFO == 0 .or. NINFO == 10 )then
      NINFO = 0
      
      QOUT = 0.
      QIN  = 0.
      do LQ = 2,LA
        if( QSUME(LQ) < 0. )then
          QOUT = QOUT - QSUME(LQ)
        else
          QIN  = QIN  + QSUME(LQ)
        endif
      enddo

      IQLO(1) = QOUT
      IQHI(1) = QIN
      
      ! ***  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES
      QOUT = 0.
      QIN  = 0.
      do K = 1,KC
        do LL = 1,NPBS
          LQ = LPBS(LL)
          LN = LNC(LQ)
          QOPEN = VHDX2(LN,K)
          if( QOPEN < 0. )then
            QOUT = QOUT + ABS(QOPEN)
          else
            QIN  = QIN  + ABS(QOPEN)
          endif
        enddo
      enddo
      IQLO(2) = QOUT
      IQHI(2) = QIN
  
      QOUT = 0.
      QIN  = 0.
      do K = 1,KC
        do LL = 1,NPBW
          LQ = LPBW(LL)
          QOPEN = UHDY2(LQ+1,K)
          if( QOPEN < 0. )then
            QOUT = QOUT + ABS(QOPEN)
          else
            QIN  = QIN  + ABS(QOPEN)
          endif
        enddo
      enddo
      IQLO(3) = QOUT
      IQHI(3) = QIN
      
      QOUT = 0.
      QIN  = 0.
      do K = 1,KC
        do LL = 1,NPBE
          LQ = LPBE(LL)
          QOPEN = UHDY2(LQ,K)
          if( QOPEN > 0. )then
            QOUT = QOUT + ABS(QOPEN)
          else
            QIN  = QIN  + ABS(QOPEN)
          endif
        enddo
      enddo
      IQLO(4) = QOUT
      IQHI(4) = QIN
  
      QOUT = 0.
      QIN  = 0.
      do K = 1,KC
        do LL = 1,NPBN
          LQ = LPBN(LL)
          QOPEN = VHDX2(LQ,K)
          if( QOPEN > 0. )then
            QOUT = QOUT + ABS(QOPEN)
          else
            QIN  = QIN  + ABS(QOPEN)
          endif
        enddo
      enddo
      IQLO(5) = QOUT
      IQHI(5) = QIN
      IQLO(0) = SUM(IQLO(1:5))
      IQHI(0) = SUM(IQHI(1:5))
      
      if( NWSER > 0 )then
        UTMP = ATAN2(WNDVELE(L),WNDVELN(L))/PI*180.
        if( UTMP < 0. )UTMP = UTMP+360.
      else
        UTMP = 0.
      endif
      if( SVPAT(L) > 0. )then
        VTMP = VPAT(L)/SVPAT(L)*100.
      else
        VTMP = 0.
      endif
      if( ISICE > 0 )then
        CKC = ICETHICK(L)*100.
      else
        CKC = 0.
      endif
      
      write(*,'(/)')
      write(*,'(A)')'    TIME WSPEED DIRTO  TAIR  RELH  SOLR  RAIN  EVAP ICETH    INFLOW   OUTFLOW'
      write(*,'(A)')'    DAYS    M/S   DEG     C   PER  W/M2  MM/D  MM/D    CM      M3/S      M3/S'
      write(*,'(A)')'-------------------------------------------------------------------------------'
      write(*,9300) TIME, INT(WINDST(L),4), INT(UTMP,4), INT(TATMT(L)+.5,4),INT(VTMP,4), INT(SOLSWRT(L),4),  & 
                          INT(RAINT(L)*86400000.,4), INT(EVAPT(L)*86400000.,4), INT(CKC,4),IQHI(0),IQLO(0)
      NOPEN = NPBS+NPBW+NPBE+NPBN
      if( NOPEN > 0 )then
        write(*,9400)'OPEN: NORTH',IQHI(5),IQLO(5)
        write(*,9400)'OPEN: EAST ',IQHI(4),IQLO(4)
        write(*,9400)'OPEN: WEST ',IQHI(3),IQLO(3)
        write(*,9400)'OPEN: SOUTH',IQHI(2),IQLO(2)
        write(*,9401)'OTHER BCs  ',IQHI(1),IQLO(1)
      else
        write(*,9402)
      endif
    endif
    
    ! *** ESTIMATE TIME TO COMPLETION
    if( N > 1 )then
      TCGRS = DSTIME(1)
      T1 = TBEGIN*TCON
      T2 = (TBEGIN*TCON + TIDALP*NTC)
      TSPEED = TCGRS/(TIMESEC - T1)
      ETA = (T2 - TIMESEC)*TSPEED/3600.
      ETA = MIN(ETA,99999.989)
      T1 = TCGRS/3600.
      write(*,'('' ** ELAPSED TIME: '',F8.2,'' (hr)   ESTIMATED TIME TO COMPLETION:'',F10.2,'' (hr)'')') T1,ETA
    endif
        
    write(*,'(A)')'--------------------------------------------------------------------------------'
    if( ISDYNSTP > 0 )then
      write(*,'(A)')'    TIME     TIME   ELEV VELE VELN  '//PARM//'   AV    AB  VELE  VELN  '//PARM//'   AV      '
      write(*,'(A)')'     IN      STEP   SURF SURF SURF  SUR  SURF  SURF  BOTT  BOTT  BOTT  BOTT     '
      write(*,'(A)')'    DAYS      SEC     CM CM/S CM/S  '//UNITS//'  CM/S  CM/S  CM/S  CM/S  '//UNITS//'  CM/S  LMIN'
    else
      write(*,'(A)')'    TIME     TIME   ELEV VELE VELN  '//PARM//'   AV    AB  VELE  VELN  '//PARM//'   AV    AB '
      write(*,'(A)')'     IN      STEP   SURF SURF SURF  SUR  SURF  SURF  BOTT  BOTT  BOT  BOTT  BOTT'
      write(*,'(A)')'    DAYS      SEC     CM CM/S CM/S  '//UNITS//'  CM/S  CM/S  CM/S  CM/S  '//UNITS//'  CM/S  CM/S'
    endif
    write(*,'(A)')'--------------------------------------------------------------------------------'
    ETIME = DT*N/86400.
    if( LMHK .and. ETIME > 0.01 )WRITE(*,'("SUPPORT ENERGY LOSS",F10.4," MW-hr")')SUM(ESUP(:,:))
    if( LMHK .and. ETIME > 0.01 )WRITE(*,'("MHK ENERGY OUTPUT  ",F10.4," MW-hr")')SUM(EMHK(:,:))
    if( LMHK .and. ETIME > 0.01 )WRITE(*,'("MHK POWER OUTPUT   ",F10.4," kW")')SUM(PMHK(:,:))*1E-3
    
  endif

  ! *** INCREMENT THE SCREEN COUNTER
  JSHPRT = JSHPRT + 1
  if( JSHPRT < ISHPRT ) return
  
  ! *** INCREMENT THE SCREEN COUNTER
  JSHPRT = 1
  NSHOWC = NSHOWC+1

  LN = LNC(L)
  ZSURF = (HP(L)+BELV(L))*100.
  UTMP = 0.5*STCUV(L)*(U(LEC(L),KC)+U(L,KC))*100.
  VTMP = 0.5*STCUV(L)*(V(LN,KC)+V(L,KC))*100.
  VELEKC = CUE(L)*UTMP+CVE(L)*VTMP

  
  VELNKC = CUN(L)*UTMP+CVN(L)*VTMP
  UTMP = 0.5*STCUV(L)*(U(LEC(L),KSZ(L))+U(L,KSZ(L)))*100.
  VTMP = 0.5*STCUV(L)*(V(LN,KSZ(L))+V(L,KSZ(L)))*100.
  VELEKB = CUE(L)*UTMP+CVE(L)*VTMP
  VELNKB = CUN(L)*UTMP+CVN(L)*VTMP
  AVKS = MIN(AV(L,KS)*10000.*HP(L),99999.)
  AVKB = MIN(AV(L,KSZ(L))*10000.*HP(L),99999.)
  ABKS = MIN(AB(L,KS)*10000.*HP(L),99999.)
  ABKB = MIN(AB(L,KSZ(L))*10000.*HP(L),99999.)

  IZSURF = NINT(ZSURF)
  IVELEKC = NINT(VELEKC)
  IVELNKC = NINT(VELNKC)
  IAVKS = NINT(AVKS)
  IABKS = NINT(ABKS)
  IVELEKB = NINT(VELEKB)
  IVELNKB = NINT(VELNKB)
  IAVKB = NINT(AVKB)
  IABKB = NINT(ABKB)
  
  ! *** CONTROL SIZE TO PREVENT FORMAT ERRORS
  IAVKS = MIN(IAVKS,99999)
  IABKS = MIN(IABKS,99999)
  IAVKB = MIN(IAVKB,99999)
  IABKB = MIN(IABKB,99999)

  ! *** CONSTITUENTS
  CKC = 0.
  CKB = 0.
  if( IPARAM == 1 .and. ISTRAN(1) > 0  )then
    CKC = SAL(L,KC)
    CKB = SAL(L,KSZ(L))

  elseif( IPARAM == 2 .and. ISTRAN(2) > 0  )then
    CKC = TEM(L,KC)
    CKB = TEM(L,KSZ(L))

  elseif( IPARAM == 3 .and. ISTRAN(3) > 0  )then
    CKC = DYE(L,KC,1)
    CKB = DYE(L,KSZ(L),1)

  elseif( IPARAM == 5 .and. ISTRAN(5) > 0 )then
    if( ISUB >0 .and. ISUB <= NTOX )then
      CKC = TOX(L,KC,ISUB)
      CKB = TOX(L,KSZ(L),ISUB)
    else
      CKC = 0.
      do I = 1,NTOX
        CKC = CKC+TOX(L,KC,I)
      enddo
      CKB = 0.
      do I = 1,NTOX
        CKB = CKB+TOX(L,KSZ(L),I)
      enddo
    endif

  elseif( IPARAM == 6 .and. ISTRAN(6) > 0 )then
    if( ISUB >0 .and. ISUB <= NSED )then
      CKC = SED(L,KC,ISUB)
      CKB = SED(L,KSZ(L),ISUB)
    else
      CKC = SEDT(L,KC)
      CKB = SEDT(L,KSZ(L))
    endif

  elseif( IPARAM == 7 .and. ISTRAN(7) > 0 )then
    if( ISUB >0 .and. ISUB <= NSND )then
      CKC = SND(L,KC,ISUB)
      CKB = SND(L,KSZ(L),ISUB)
    else
      CKC = SNDT(L,KC)
      CKB = SNDT(L,KSZ(L))
    endif

  elseif( IPARAM == 0 .and. (ISTRAN(6) > 0 .or. ISTRAN(7) > 0 ) )then
    CKC = (SEDT(L,KC)+SNDT(L,KC))
    CKB = (SEDT(L,KSZ(L))+SNDT(L,KSZ(L)))

  elseif( IPARAM == 8 .and. ISTRAN(8) > 0 )then
    CKC = WQV(L,KC,ISUB)
    CKB = WQV(L,KSZ(L),ISUB)

  endif

  ICKC = MIN(NINT(CKC),99999)
  ICKB = MIN(NINT(CKB),99999)
  if( ISDYNSTP > 0 )then
    write(*,9100)TIME,DELT,IZSURF,IVELEKC,IVELNKC,ICKC,IAVKS,IABKS,IVELEKB,IVELNKB,ICKB,IAVKB,LMINSTEP
  else
    write(*,9100)TIME,DELT,IZSURF,IVELEKC,IVELNKC,ICKC,IAVKS,IABKS,IVELEKB,IVELNKB,ICKB,IAVKB,IABKB
  endif

  return

  9000 FORMAT(/' AUTOSTEPPING SUMMARY (WITH SAFETY FACTOR):',/       &
               '   METHOD1:    MOMENTUM CHECK (DT,L): ',F10.4,I5,/   &
               '   METHOD2:   ADVECTION CHECK (DT,L): ',F10.4,I5,/   &
               '   METHOD3:   BTM FRICT CHECK (DT,L): ',F10.4,I5,/,: &
               '   METHOD4: LIMIT DH/DT CHECK (DT,L): ',F10.4,I5,/)
  9100 FORMAT(F9.3,F8.3,I7,3I5,4I6,I5,2I6)
  9300 FORMAT(F9.3,8I6,2I10,/)
  9400 FORMAT(46X,A11,2I10)
  9401 FORMAT(46X,A11,2I10,//)
  9402 FORMAT(//)

END
