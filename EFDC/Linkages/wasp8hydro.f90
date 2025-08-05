! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
  ! ******************************************************************************************
  ! *** Writes the hyd file for WINASP/WASP-HYDRO water quality model.

  Subroutine WASP8HYDRO

  ! *** ----- Previous subroutine name -- SUBROUTINE WASP7EPA
  ! *** ------------------------------------------------------------------------------
  ! ***  PURPOSE:
  ! ***
  ! ***    Subroutine WASPHYDRO writes the hyd file for WINASP water quality model.
  ! ***
  ! ***  MODIFICATION HISTORY:
  ! ***
  ! ***    Date       Author         Comments
  ! ***    ---------- -------------- -------------------------------------------------
  ! ***    2013-05-08 Paul M. Craig  Converted to f90 and cleaned code
  ! ***     5/23/2006 Hugo Rodriguez Add a limit maximum dispersion specific to certain cells
  ! ***     4/11/2006 Hugo Rodriguez Add qfactor
  ! ***    10/27/2005 Hugo Rodriguez Calculate the initial date of wasp instead of reading it
  ! ***    06/24/2005 Hugo Rodriguez Limit Ab to Abwmax
  ! ***    06/15/2005 Hugo Rodriguez Zero out AH coefficient transferred to wasp IF
  ! ***                              ISHDMF<2
  ! ***    05/04/2005 Hugo Rodriguez Limit maximum dispersion to cell volume/time step
  ! ***    05/02/2005 Hugo Rodriguez Eliminate dispersion for boundary flows.
  ! ***    30/10/2003 Hugo Rodriguez This version uses a dll to transfer the hydrodynamic
  ! ***               Tim Wool       data, including salinity and temperature, from EFDC
  ! ***                              to WASP version 6.1.
  ! *** ------------------------------------------------------------------------------
  ! ***
  use GLOBAL
#ifndef GNU  
  USE IFPORT
#endif
  use INFOMOD,only:SKIPCOM,READSTR

#ifdef WASPOUT
  implicit none

  integer(IK4) :: LT, L, IYEAR, IMON, IDAY, I, IS, LCLTM2, LTYPE, KWASP, LBELOW, LAUX, K, IONE,  NCTL
  integer(IK4) :: LWSPTMP, IZERO, IM1, IERROR, NJUN, IOS, KK, NCHNH, NCHNV, NCHN, NQ, NS, KMUL, LDTM, LUTM
  integer(IK4) :: LSLT, LE, LN, LW, LS, KMUL1, KMUL2, IPTMP, JMTMP, JPTMP, LCELTMP

  integer(IK4), save :: IBEGIN,IDAYS,nsg,nf

  real(RKD)    :: SVPT,SCALR,WSS1,WSS2,WSS3,VOLUME,DXYSUM,VELX,VELY,VELZ,VELMAG,ADDLW,ADDLS,FLOWX
  real(RKD)    :: UDDXTMP,ADDL,ADDL1,TMPVAL,FLOWY,VDDYTMP,FLOWZ,WDDZTMP,VOLUM,DEPTH,dt0, dtwasp

  !INCLUDE 'EFDC_WASPHYDRO.CMN'
  character*80 STR*200
  character(80)    :: errstring
  logical(4) RES

  ! *** DECLARE PERSISTENT VARIABLES
  real(RK4),    save    :: AD,ADCOEFF,ABWMAX
  real(RK4),    save,allocatable,dimension(:) :: ABWMX
  real(RK4),    save,allocatable,dimension(:) :: crnu        !  crnu(nf),brintt(nf),flow(nf)
  real(RK4),    save,allocatable,dimension(:) :: brintt
  real(RK4),    save,allocatable,dimension(:) :: flow
  real(RK4),    save,allocatable,dimension(:) :: SegVolume   !  SegVolume(nsg),SegDepth(nsg),SegVel(nsg),SegSalt(nsg),SegTemp(nsg)
  real(RK4),    save,allocatable,dimension(:) :: SegDepth
  real(RK4),    save,allocatable,dimension(:) :: SegVel
  real(RK4),    save,allocatable,dimension(:) :: SegSalt
  real(RK4),    save,allocatable,dimension(:) :: SegTemp

  character*6,save,allocatable,dimension(:) :: SSN          !  SSN(lcm)

  !parameter(nf = 31000,nsg = 12000)
  dimension LAUX(ICM,JCM,KCM)

  integer(IK4) LU,LD,IU,ID,JU,JD,NAUX
  integer(IK4) istartyear,istartmonth,istartday,istarthour,istartminute,istartsecond,Ihl_debug,Ihl_mode,inumsegconsts,j
  !INTEGER FLAGWASPBC(NQSERM,KCM)

  character*6  sn1
  character*23 segname
  character*23 SN2
  character*80 DESCRIPTION,MODELERNAME
  character*3 Itext, Jtext, Ktext
  character*80 RECHNR

  real(RKD)  AUX,thour,tmin,tsec
  real(RK4)  rinterval
  real(RK4)  vol1,vol2

  ! ***
  ! ***
  ! *** **********************************************************************! ***
  ! ***
  ! ***  READ CONTROL DATA FILE EFDC.WSP
  ! ***
  ! *** ----------------------------------------------------------------------!
  ! ***

  ! *** ALLOCATE LOCAL ARRAYS
  if( .not. allocated(ABWMX) )then
    nsg = lcm*kcm

    allocate(ABWMX(nsg))
    allocate(SSN(LCM))
    allocate(SegVolume(nsg))
    allocate(SegDepth(nsg))
    allocate(SegVel(nsg))
    allocate(SegSalt(nsg))
    allocate(SegTemp(nsg))
  endif

  SVPT = 1.
  if( RESSTEP < TIDALP ) SVPT = 0.
  ! *** Showing fixed time step 10.s in WASP interface when dynamic time step is used
  dt0 = DELT
  if( ISDYNSTP == 0 )then
    dtwasp = DELT
  else
    dtwasp = 10.
  endif
  ! ***
  if( JSWASP == 1 )then

    ! ***       for jswasp = 1 only first entry
    write(*,'(A)')'READING EFDC.WSP'
    open(1,FILE = 'EFDC.WSP',STATUS = 'UNKNOWN')
    write(6,*)'EFDC.WSP opened'

    ! *** C1**  READ CELL VOLUME parameterS (PMC - None of these paramters are being used.  Keep for legacy info)
    read(1,1)
    read(1,1)
    read(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
    write(*,*)'<EFDC.WSP1> ',IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP

    ! *** C2**  READ DIFFUSION parameterS (PMC - The 1st 4 of these paramters are not being used.  Keep for legacy info)
    read(1,1)
    read(1,1)
    read(1,*) NRFLD,SCALR,CONVR,ISNKH,ADCOEFF,ABWMAX
    write(*,*)'<EFDC.WSP2> ',NRFLD,SCALR,CONVR,ISNKH,ADCOEFF,ABWMAX
    do LT = 2,LALT
      l = lij(illt(lt),jllt(lt))
      ABWMX(l) = ABWMAX
    enddo

    ! *** C3**  READ ADVECTION parameterS (PMC - Only HYDFIL and IDAYS are being used.  Keep for legacy info)
    read(1,1)
    read(1,1)
    read(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS
    write(*,*)'<EFDC.WSP3> ', IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS

    ! *** C4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)  (PMC - Only DEPSED is being used.  Keep for legacy info)
    read(1,1)
    read(1,1)
    read(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3
    write(*,*)'<EFDC.WSP4> ',DEPSED,TDINTS,SEDIFF,WSS1,WSS2,WSS3

    ! *** C5**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
    read(1,1)
    read(1,1)
    read(1,*) iyear,imon,iday
    write(*,*)'<EFDC.WSP5> ',iyear,imon,iday

    do i = 1,5
      read(1,*,err = 11) RECHNR
      write(*,*)'<EFDC.WSP6> ',RECHNR
    enddo
11  continue

    do LT = 2,LALT
      read(1,*,err = 12)i,j,ABWMAX
      l = lij(i,j)
      ABWMX(l) = ABWMAX
      write(*,*)'<EFDC.WSP7> ',i,j,ABWMAX
    enddo
12  continue

    close(1)
    write(*,*)'EFDC.WSP read succesfully, now WRITE ABMax'

    ! *** CREATE THE WASP FOLDER
#ifdef GNU
    RES = SYSTEM( 'mkdir -p ./' // trim(OUTDIR)//'wasp')
#else       
    RES = MAKEDIRQQ(OUTDIR//'wasp')
#endif 

    open(1,FILE = OUTDIR//'wasp/ABmax.txt',STATUS = 'UNKNOWN')
    write(1,*)'    I    J     ABmax'
    do LT = 2,LALT
      l = lij(illt(lt),jllt(lt))
      write(1,21)illt(lt),jllt(lt),ABWMX(l)
    enddo
    close(1)
21  format(2I5,f10.6)
    write(6,*)'EFDC.WSP read succesfully and ABmax.txt written'
1   FORMAT (80X)

    ! *** read qser file to check for flows only to some layers  (PMC - Not used for anything.  Keep for Legacy)
    if( NQSER >= 1 )then
      open(1,FILE = 'QSER.INP',STATUS = 'UNKNOWN')
      STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES
    endif

    close(1)
    goto 862
860 WRITE(6,861)
861 FORMAT('  READ ERROR FOR FILE QSER.INP ')
862 continue

    ! *** **********************************************************************!

    ! ***  DEFINE EFDC-WASP CORRESPONDENCE AND INITIALIZE FILES

    ! *** **********************************************************************!
    ! ***  begin code inserted by Andy Stoddard 2-26-2008
    ! *** ------ Feb-26-2008 A Stoddard add code to WRITE WASP8SEG_EFDCIJK
    ! *** ------ file listing of EFDC L,IJK cells and Wasp segment/layer data
    ! *** ----------------------------------------------------------------------!
    ! *** ----------------------------------------------------------------------!
    ! *** ---- 7/1/2005 A Stoddard added code to WRITE WASP8SEG_EFDCIJK.DAT file
    ! *** ---- Segment Name file is linked for import to WASP6 & WASP7 & WASP8
    open(970,FILE = OUTDIR//'wasp\WASP8SEG_EFDCIJK.DAT',STATUS = 'UNKNOWN')  ! 7/1/2005 AS
    close(970,STATUS = 'DELETE')                                      ! 7/1/2005 AS

    ! ***    open(90,FILE = 'WASPP.OUT',STATUS = 'UNKNOWN')
    ! ***    open(93,FILE = 'WASPC.OUT',STATUS = 'UNKNOWN')
    ! ***    close(90,STATUS = 'DELETE')
    ! ***    close(93,STATUS = 'DELETE')
    ! ***    open(90,FILE = 'WASPP.OUT',STATUS = 'UNKNOWN')
    ! ***    open(93,FILE = 'WASPC.OUT',STATUS = 'UNKNOWN')
    open(970,FILE = OUTDIR//'wasp\WASP8SEG_EFDCIJK.DAT',STATUS = 'UNKNOWN')  ! 7/1/2005 AS

    write(970,9702)              ! WRITE header records

    LCLTM2 = LCLT-2
    LWASP = 0
    if( KC > 1 )then
      LTYPE = 1
      KWASP = 1
      do LT = 2,LALT
        LWASP = LWASP + 1
        LBELOW = LWASP + LCLTM2
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        DMULT = HLPF(L)*DZC(L,KC)
        VOLUME = DXYP(L)*HLPF(L)*DZC(L,KC)
        if( RESSTEP < TIDALP )then
          DMULT = HP(L)*DZC(L,KC)
          VOLUME = DXYP(L)*HP(L)*DZC(L,KC)
        endif
        LAUX(I,J,KC) = LWASP

        write(970,9701)L,I,J,KC,LWASP,KWASP         ! ,I,J,L,KC               ! 7/1/2005 AS
      enddo
      LTYPE = 2
      do K = KS,2,-1
        KWASP = KC-K + 1
        do LT = 2,LALT
          LWASP = LWASP + 1
          LBELOW = LWASP + LCLTM2
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          DMULT = HLPF(L)*DZC(L,K)
          VOLUME = DXYP(L)*HLPF(L)*DZC(L,K)
          if( RESSTEP < TIDALP )then
            DMULT = HP(L)*DZC(L,K)
            VOLUME = DXYP(L)*HP(L)*DZC(L,K)
          endif
          LAUX(I,J,K) = LWASP
          write(970,9701)L,I,J,K,LWASP,KWASP         ! ,I,J,L,K               ! 7/1/2005 AS
        enddo
      enddo
    endif
    LTYPE = 2
    if( KC == 1 ) LTYPE = 1
    KWASP = KC
    do LT = 2,LALT
      LWASP = LWASP + 1
      LBELOW = LWASP + LCLTM2
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      DMULT = HLPF(L)*DZC(L,KC)
      VOLUME = DXYP(L)*HLPF(L)*DZC(L,KC)
      if( RESSTEP < TIDALP )then
        DMULT = HP(L)*DZC(L,KC)
        VOLUME = DXYP(L)*HP(L)*DZC(L,KC)
      endif
      IONE = 1
      LAUX(I,J,1) = LWASP
      write(970,9701)L,I,J,IONE,LWASP,KWASP         ! ,I,J,L,IONE               ! 7/1/2005 AS
    enddo
    LTYPE = 3
    KWASP = KC + 1
    DXYSUM = 0.
    LWSPTMP = LWASP + 1
    do LT = 2,LALT
      LWSPTMP = LWSPTMP + 1
    enddo
    ! ***  The following the lower benthic layer.  All upper benthic layer segments
    ! ***  have this layer immediately below them:
    do LT = 2,LALT
      LWASP = LWASP + 1
      LBELOW = LWSPTMP
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      DXYSUM = DXYSUM + DXYP(L)
      VOLUME = DXYP(L)*DEPSED
      IZERO = 0
    enddo
    ! ***  Next do the lower benthic layer:
    LTYPE = 4
    KWASP = KC + 2
    LWASP = LWASP + 1
    LBELOW = 0
    DMULT = DEPSED
    VOLUME = DXYSUM*DEPSED
    IM1 = -1

    close(970)                           ! 7/1/2005

    ! ***  1001 FORMAT(6I5,2F10.4)
    ! ***  1031 FORMAT(2I5,F10.4,10X,A50)
    ! ***  1032 FORMAT(2F10.4)
    ! ***  1033 FORMAT(3I10,F10.1,4F10.3,'   !',4i5)


    ! ***  9701 FORMAT('I = ',I3,' J = ',I3,' L = ',I4 ,' K = ',I2,
    ! *** +       ' W',I4,' WK = ',I2)

9702 format('######################################################',/,  &
      '# File    = WASP8SEG_EFDCIJK.DAT',/,  &
      '# Field#1 = EFDC Cell L-Number',/,  &
      '# Field#2 = EFDC Cell I-Index',/,  &
      '# Field#3 = EFDC Cell J-Index',/,  &
      '# Field#4 = EFDC Cell K-Layer  [1 = Bottom;KC = Surface]',/,  &
      '# Field#5 = WASP Segment Number',/,  &
      '# Field#6 = WASP Segment Layer [KC = Bottom;1 = Surface]',/,  &
      '######################################################')

9701 FORMAT(I4,',', &            ! EFDC L-coordinate
      I4,',',     &        ! EFDC I-coordinate
      I5,',',     &        ! EFDC J-cell index
      I2,',',     &        ! EFDC K-layer number   1 = bottom;  KC = surface layer
      I5,',',     &        ! WASP Segment Number
      I2)                  ! WASP Layer Number     KC = bottom; 1= surface layer

    ! ***         WASPSEG.OUT data structure
    ! ***         Wasp8 wired for linkage of A30 character string file
    ! ***         W Wasp Segment Number
    ! ***         WK wasp Layer number 1 = sfc   K = KC = bottom
    ! ***         I,J = EFDC I, J index
    ! ***         L = EFDC grid cell number
    ! ***         K = EFDC layer K = KC sfc    K = 1 bottom

    ! ***  END code inserted by Andy Stoddard 2-26-2008
    ! *** **********************************************************************!

    ! *** ----------------------------------------------------------------------!
    ! ***
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! ***  hlopen parameters
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    Ihl_handle = 0
    Ihl_debug = 1
    Ihl_mode = 1  !Ihl_mode = 0 to READ from dll hyd file,  = 1 to WRITE to dll hyd file
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the debug flag 0 = No debug 1 = Debug (LOG.OUT)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetdebug(Ihl_debug)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Open the file
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlopen(HYDFIL,Ihl_mode,Ihl_handle,ierror)
    if( ierror > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000)ierror,errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the language to FORTRAN
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetlanguage(Ihl_handle,1,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Store a description string
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    DESCRIPTION = '   '//char(0)
    call hladddescription(Ihl_handle,0,DESCRIPTION,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Store the modeler name
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    MODELERNAME = 'Created by:  '//char(0)
    call hladddescription(Ihl_handle,1,MODELERNAME,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the creator
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetcreator(Ihl_handle, 1, ierror)   ! *** 1 for EFDC code
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the seed moment (start date of the hyd file)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Calculate istartmonth, istartday,istartyear
    j = 1461*(iyear + 4800 + int((imon-14)/12))/4 + 367*(imon-2-12*int((imon-14)/12))/12-(3*((iyear + 4800 + int((imon-14)/12) + 100)/100))/4 + iday-32075 + tbegin + IDAYS
    istartyear = 100*(int(4*(j + 68569)/146097)-49) + int(4000*(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4) + 1)/   &
      1461001) + int(int(80*int(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4)-1461*int(4000*(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4) + 1)/1461001)/4 + 31)/2447)/11)
    istartmonth = int(80*int(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4)-1461*int(4000*(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4) + 1)/1461001)/4 + 31)/2447) + 2-12*int(int(80*    &
      int(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4)-1461*int(4000*(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4) + 1)/1461001)/4 + 31)/2447)/11)
    istartday = int(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4)-1461*int(4000*(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4) + 1)/1461001)/4 + 31)-2447*int(80*int(  &
      int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4)-1461*int(4000*(int(j + 68569-(146097*int(4*int(j + 68569)/146097) + 3)/4) + 1)/1461001)/4 + 31)/2447)/80
    istarthour = 0
    istartminute = 0
    istartsecond = 0

    ! *** Adjust hour, minutes and seconds for the start day IF necessary
    IBEGIN = IDAYS*NTSPTC

    ! *** Update the number of WASP time step per hydrodynamic file read
    if( ISDYNSTP == 0 )then
      NMMT = RESSTEP/dtwasp
    else
      NMMT = NTC*TIDALP/dtwasp/(NWASPOUT - 1)
    endif

    AUX = FLOAT(IBEGIN)/FLOAT(NMMT)
    NAUX = INT(AUX)
    if( AUX > NAUX )then
      AUX = (NAUX + 1.-AUX)
      tsec = AUX*NMMT*TCON/NTSPTC
      istartsecond = INT(tsec)
      if( tsec > 60 )then
        TMIN = tsec/60.0
        istartminute = INT(TMIN)
        istartsecond = INT((TMIN-istartminute)*60.)
        if( IStartminute > 60 )then
          THOUR = FLOAT(istartminute)/60.0
          istarthour = INT(THOUR)
          istartminute = INT((THOUR-istarthour)*60.)
        endif
      endif
    endif

    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetseedmoment(Ihl_handle,istartmonth,istartday,istartyear,istarthour,istartminute,istartsecond,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 

    call hlsetnumlayers(Ihl_handle,INT(kc,4),ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the number of segments
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    NJUN = KC*(LCLT-2)
    if( njun > nsg )then
      write(6,500)NJUN,NSG
      STOP
    endif
500 FORMAT('THE NUMBER OF WASP SEGMENTS IN YOUR APPLICATION',I6,1x,'IS GREATER THAN THE ARRAY dimension:',I7)

    call hlsetnumsegments(Ihl_handle,njun,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! ***           Set WASP Segment Names (defaulft as I,J,K)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    do LT = 2,LALT
      l = lij(illt(lt),jllt(lt))
      SSN(l) = '      '//char(0)
    enddo
    open(94,FILE = 'segname.inp',iostat = ios,STATUS = 'old')
    if( ios == 0 )then
      STR = READSTR(94)  ! *** SKIP OVER TITLE AND AND HEADER LINES
      do kk = 1,la
        read(94,*,err = 111)I,J,SN1
        L = LIJ(I,J)
        SSN(l) = sn1
      enddo
    endif
111 close(94)
    I = 0
    do K = KC,1,-1
      do LT = 2,LALT
        I = I + 1
        LAUX(ILLT(LT),JLLT(LT),K) = I
        l = lij(illt(lt),jllt(lt))
        if( il(l) >= 100 )then
          write(itext,"(I3)")IL(L)
        elseif( il(l) >= 10 )then
          write(itext,"(I2)")IL(L)
        else
          write(itext,"(I1)")IL(L)
        endif
        if( jl(l) >= 100 )then
          write(jtext,"(I3)")jL(L)
        elseif( jl(l) >= 10 )then
          write(jtext,"(I2)")jL(L)
        else
          write(jtext,"(I1)")jL(L)
        endif
        if( k >= 100 )then
          write(ktext,"(I3)")k
        elseif( k >= 10 )then
          write(ktext,"(I2)")k
        else
          write(ktext,"(I1)")k
        endif
        sn2 = ' I = '//itext//' J = '//jtext//' K = '//ktext
        segname = sn2//char(0)
        call hlsetsegname(ihl_handle,i,segname,ierror)
      enddo
    enddo
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the number of flow paths
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 

    ! *** Initializing the horizontal U and V flow Connections
    NCHNH = 0
    NCHNV = 0
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      NCHNH = NCHNH + INT(SUBO(L))                   ! *** U face
      if( IJCTLT(I + 1,J) == 8 )then
        if( SUBO(LEC(L)) == 1. ) NCHNH = NCHNH + 1   ! *** Add open BC cell to the East
      endif
      NCHNH = NCHNH + INT(SVBO(L))                   ! *** V face
      if( IJCTLT(I,J + 1) == 8 )then
        if( SVBO(LNC(L)) == 1.) NCHNH = NCHNH + 1    ! *** Add open BC cell to the North
      endif
      NCHNV = NCHNV + INT(SWB(L))
    enddo
    NCHN = KC*NCHNH + (KC-1)*NCHNV                   ! *** delme - does not work with SGZ grids

    ! *** Remove flow BC's
    NQ = NQSIJ
    do L = 1,NQSIJ
      if( LIJLT(BCPS(L).I,BCPS(L).J) == 0 ) NQ = NQ-1
    enddo
    NCHN = NCHN + KC*NQ

    ! *** Remove hydraulic structures
    NQ = NQCTL
    do NCTL = 1,NQCTL
      if( LIJLT(HYD_STR(NCTL).IQCTLU,HYD_STR(NCTL).JQCTLU) == 0 )then
        if( LIJLT(HYD_STR(NCTL).IQCTLD,HYD_STR(NCTL).JQCTLD) == 0 ) NQ = NQ-1
      endif
    enddo
    NCHN = NCHN + KC*NQ

    if( .not. allocated(crnu) )then
      allocate(crnu(NCHN))
      allocate(brintt(NCHN))
      allocate(flow(NCHN))
    endif
    !  DELME
    !  if( NCHN > NF )then
    !    write(6,600)NCHN,NF
    !    STOP
    !  endif
    !600 FORMAT('THE NUMBER OF WASP FLOWS IN YOUR APPLICATION',I6,1X,'IS GREATER THAN THE ARRAY dimension:',I7)

    call hlsetnumflowpaths(Ihl_handle, NCHN, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif

    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the flow path and direction
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    LCLTM2 = LCLT-2
    LWASP = 0
    do K = KC,1,-1
      KMUL = KC-K
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        if( SUBO(L) == 1. )then
          LWASP = LWASP + 1
          LDTM = LT-1 + KMUL*LCLTM2
          LUTM = LDTM-1
          if( IJCTLT(I-1,J) == 8 ) LUTM = 0 ! *** West open boundary cells
          call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
          if( ierror  > 0 )then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            STOP
          endif
        endif

        if( IJCTLT(I + 1,J) == 8 )then ! *** East open boundary cells
          if( SUBO(LEC(L)) == 1. )then
            LWASP = LWASP + 1
            LDTM = 0
            LUTM = LT-1 + KMUL*LCLTM2
            call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
            if( ierror  > 0 )then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              STOP
            endif
          endif
        endif
      enddo

      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        if( SVBO(L) == 1. )then
          LWASP = LWASP + 1
          LSLT = LSCLT(LT)
          LDTM = LT-1 + KMUL*LCLTM2
          LUTM = LSLT-1 + KMUL*LCLTM2
          if( IJCTLT(I,J-1) == 8 ) LUTM = 0 ! *** South open boundary cells
          call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,2,ierror)
          if( ierror  > 0 )then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            STOP
          endif
        endif
        if( IJCTLT(I,J + 1) == 8 )then ! *** North open boundary cells
          LN = LNC(L)
          if( SVBO(LN) == 1. )then
            LWASP = LWASP + 1
            LSLT = LSCLT(LT)
            LDTM = 0
            LUTM = LT-1 + KMUL*LCLTM2
            call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,2,ierror)
            if( ierror  > 0 )then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              STOP
            endif
          endif
        endif
      enddo
    enddo

    do K = KC,1,-1
      do LT = 1,NQSIJ
        I = BCPS(LT).I
        J = BCPS(LT).J
        if( LIJLT(I,J) == 0 ) GOTO 100
        NS = BCPS(Lt).NQSERQ
        !IF(flagwaspbC(ns,k) == 1 ) GOTO 100
        LWASP = LWASP + 1
        LDTM = LAUX(I,J,K)
        LUTM = 0
        ! *** According to WASP document, LDTM and LUTM have to be a number 1 to number of segments - DKT
        call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif
100   enddo
    enddo

    do K = KC,1,-1
      do NCTL = 1,NQCTL
        I = HYD_STR(NCTL).IQCTLU
        J = HYD_STR(NCTL).JQCTLU
        LUTM = LAUX(I,J,K)
        I = HYD_STR(NCTL).IQCTLD
        J = HYD_STR(NCTL).JQCTLD
        LDTM = LAUX(I,J,K)
        if( LUTM == 0 .and. LDTM == 0 ) GOTO 200
        LWASP = LWASP + 1
        call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif
200   enddo
    enddo

    if( KC > 1 )then
      do K = KS,1,-1
        KMUL1 = KS-K
        KMUL2 = KMUL1 + 1
        do LT = 2,LALT
          I = ILLT(LT)
          J = JLLT(LT)
          L = LIJ(I,J)
          if( SWB(L) == 1. )then
            LWASP = LWASP + 1
            LUTM = LT-1 + KMUL1*LCLTM2
            LDTM = LT-1 + KMUL2*LCLTM2
            call hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,3,ierror)
            if( ierror  > 0 )then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              STOP
            endif
          endif
        enddo
      enddo
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the number of segment constituents
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 

    inumsegconsts = 3  !volume,depth,velocity
    ! ***   if( ISTRAN(2) >= 1 ) inumsegconsts = inumsegconsts + 1   !temperature modeled in EFDC and transfered
    ! ***   if( ISTRAN(1) >= 1 ) inumsegconsts = inumsegconsts + 1   !salinity modeled in EFDC and transfered
    if( ISTRAN(2) >= 1 ) inumsegconsts = 4   !temperature modeled in EFDC and transfered
    if( ISTRAN(1) >= 1 ) inumsegconsts = 5   !salinity modeled in EFDC and transfered

    call hlsetnumsegconsts(Ihl_handle, inumsegconsts, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set the number of flow path constituents
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! ***  when we add sed transport we need to add more
    call hlsetnumfpconsts(Ihl_handle, 3, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Now we will set all the segment constituent types
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetsegconsttype(Ihl_handle, 1, 0, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetsegconsttype(Ihl_handle, 2, 1, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetsegconsttype(Ihl_handle, 3, 2, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! ***  next only IF temperature is transfered
    if( ISTRAN(2) >= 1 )then
      call hlsetsegconsttype(Ihl_handle, 4, 3, ierror)
      if( ierror  > 0 )then
        call hlgetlasterror(errstring)
        write(6,6000) ierror, errstring
        STOP
      endif
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! ***  next only IF salinity is transfered
    if( ISTRAN(1) >= 1 )then
      call hlsetsegconsttype(Ihl_handle, 5, 4, ierror)
      if( ierror  > 0 )then
        call hlgetlasterror(errstring)
        write(6,6000) ierror, errstring
        STOP
      endif
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! *** Set all the flow constituent types
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetfpconsttype(Ihl_handle, 1, 0, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    call hlsetfpconsttype(Ihl_handle, 1, 1, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    call hlsetfpconsttype(Ihl_handle, 1, 2, ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsetvartimestep(Ihl_handle,0,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    call hlsethydtimestep(Ihl_handle,REAL(dtwasp,4),ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    rinterval = RESSTEP/86400.
    call hlsetupdateint(Ihl_handle,rinterval,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 

    call hlsethydtowaspratio(Ihl_handle,INT(NMMT,4),ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
    ! ***
    ! *** INITIAL CONDITIONS WHEN IDAYS = 0
    ! ***
    if( IDAYS == 0 )then
      if( ISRESTI == 0 )then

        ! *** INITIAL CONDITION FOR A COLD START
        LWASP = 0
        do K = KC,1,-1
          do LT = 2,LALT
            I = ILLT(LT)
            J = JLLT(LT)
            L = LIJ(I,J)
            LWASP = LWASP + 1
            SegVel(LWASP) = 0.0
            if( K>= KSZ(L) )then
              SegDepth(LWASP) = HP(L)*DZC(L,K)
              SegVolume(LWASP) = SegDepth(LWASP)*DXYP(L)
              SegSalt(LWASP) = SAL(L,K)
              SegTemp(LWASP) = TEM(L,K)
              if( RESSTEP < TIDALP )then
                SegDepth(LWASP) = HP(L)*DZC(L,K)
                SegVolume(LWASP) = SegDepth(LWASP)*DXYP(L)
                SegSalt(LWASP) = SALINIT(L,K)
                SegTemp(LWASP) = TEMINIT(L,K)
              endif
            else
              SegDepth(LWASP) = 0.
              SegVolume(LWASP) = 0.
              SegSalt(LWASP) = 0.
              SegTemp(LWASP) = 0.
            endif
          enddo
        enddo
        do I = 1,NCHN
          FLOW(I) = 0.0
          CRNU(I) = 0.0
          BRINTT(I) = 0.0
        enddo

        call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        ! ***     next only IF temperature is transfered
        if( ISTRAN(2) >= 1 )then
          call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
          if( ierror  > 0 )then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            STOP
          endif
        endif
        ! ***     next only IF salinity is transfered
        if( ISTRAN(1) >= 1 )then
          call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
          if( ierror  > 0 )then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            STOP
          endif
        endif
        ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
        call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlmomentcomplete(Ihl_Handle,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif
        ! *** END OF COLD START

      else

        ! *** INITIAL CONDITIONS FROM A RESTART FILE
        LWASP = 0
        ! *** Segment Properties from the RESTART file
        do K = KC,1,-1
          do LT = 2,LALT
            I = ILLT(LT)
            J = JLLT(LT)
            L = LIJ(I,J)
            LWASP = LWASP + 1
            VELX = 0.5*(U(L,K) + U(LEC(L),K))
            VELY = 0.5*(V(L,K) + V(LNC(L),K))
            VELZ = 0.5*(W(L,K-1) + W(L,K))
            VELMAG = SQRT(VELX*VELX + VELY*VELY + VELZ*VELZ)
            SegVel(LWASP) = VELMAG
            SegDepth(LWASP) = HP(L)*DZC(L,K)
            SegVolume(LWASP) = SegDepth(LWASP)*DXYP(L)
            SegSalt(LWASP) = SAL(L,K)
            SegTemp(LWASP) = TEM(L,K)
          enddo
        enddo
        
        ! ***  Advection and dispersion in the X-direction:
        LWASP = 0
        do K = KC,1,-1
          do LT = 2,LALT
            I = ILLT(LT)
            J = JLLT(LT)
            L = LIJ(I,J)
            ADDLW = 0.0

            if( SUB(L) == 1. )then
              LW = LWC(L)
              ADDLW = DYU(L)*0.5*(AH(L,K) + AH(LW,K))*DZC(L,K)*0.5*(HP(L) + HP(LW))*DXIU (L)
              vol1 = DXYP(L)* HP(L)* DZC(L,K)
              vol2 = DXYP(LW)*HP(LW)*DZC(L,K)
              if( vol1 < vol2) vol2 = vol1
              AD = vol2/dt0*ADCOEFF
              if( addlw > AD )then
                addlw = AD
              endif
              if( IShdmf < 2 )then
                addlw = 0.
              endif
            endif
            
            if( SUBO(L) == 1. )then
              LWASP = LWASP + 1
              ! ***           if( IJCTLT(I-1,J) == 8 ) addlw = 0.0
              FLOW(LWASP) = UHDY(L,K)
              CRNU(LWASP) = 2.*UHDY(L,K)*DYIU(L)*DXIU(L)/(HP(L) + HP(LWC(L)))
              BRINTT(LWASP) = ADDLW
            endif
            if( IJCTLT(I + 1,J) == 8 )then
              LE = LEC(L)
              if( SUBO(LE) == 1. )then
                LWASP = LWASP + 1
                FLOW(LWASP) = UHDY(LE,K)
                CRNU(LWASP) = 2.*UHDY(LE,K)*DYIU(LE)*DXIU(LE)/(HP(L) + HP(LE))
                BRINTT(LWASP) = 0.0
              endif
            endif
          enddo

          ! ***  Advection and dispersion in the Y-direction:
          do LT = 2,LALT
            I = ILLT(LT)
            J = JLLT(LT)
            L = LIJ(I,J)
            ADDLS = 0.0
            LS = LSC(L)
            if( SVB(L) == 1. )then
              ADDLS = DXV(L)*0.5*(AH(L,K) + AH(LS,K))*DZC(L,K)*0.5*(HP(L) +HP(LS))*DYIV (L)
              vol1 = DXYP(L)*HP(L)*DZC(L,K)
              vol2 = DXYP(Ls)*HP(Ls)*DZC(L,K)
              if( vol1 < vol2) vol2 = vol1
              AD = vol2/dt0*ADCOEFF
              if( addls > AD )then
                addls = AD
              endif
              if( IShdmf < 2 )then
                addls = 0.
              endif
            endif
            if( SVBO(L) == 1. )then
              LWASP = LWASP + 1
              ! ***            if( IJCTLT(I,J-1) == 8 ) addls = 0.0
              FLOW(LWASP) = VHDX(L,K)
              CRNU(LWASP) = 2.*VHDX(L,K)*DYIV(L)*DXIV(L)/(HP(L) + HP(LS))
              BRINTT(LWASP) = ADDLS
            endif
            if( IJCTLT(I,J + 1) == 8 )then
              LN = LNC(L)
              if( SVBO(LN) == 1. )then
                LWASP = LWASP + 1
                FLOW(LWASP) = VHDX(LN,K)
                CRNU(LWASP) = 2.*VHDX(LN,K)*DYIV(LN)*DXIV(LN)/(HP(L) + HP(LN))
                BRINTT(LWASP) = addls
              endif
            endif
          enddo
        enddo

        ! ***  Advection and dispersion in input flows
        do K = KC,1,-1
          do LT = 1,NQSIJ
            I = BCPS(LT).I
            J = BCPS(LT).J
            if( LIJLT(I,J) == 0 ) GOTO 310
            NS = BCPS(LT).NQSERQ
            L = BCPS(LT).L
            !IF(FLAGWASPBC(NS,K) == 1 ) GOTO 310
            LWASP = LWASP + 1
            FLOW(LWASP) = BCPS(LT).RQSMUL*(QSS(K,LT) + BCPS(LT).QFACTOR*QSERT(K,NS))
            CRNU(LWASP) = FLOW(LWASP)/DXP(L)/DYP(L)/(HPK(L,K))
            ! ***          BRINTT(LWASP) = DYP(L)*AH(L,K)*HPK(L,K)/DXP(L)
            BRINTT(LWASP) = 0.0
310       enddo
        enddo
        
        ! ***   ADVECTION AND DISPERSION IN STRUCTURE FLOWS
        do K = KC,1,-1
          do NCTL = 1,NQCTL
            IU = HYD_STR(NCTL).IQCTLU
            JU = HYD_STR(NCTL).JQCTLU
            LU = LIJ(IU,JU)
            ID = HYD_STR(NCTL).IQCTLD
            JD = HYD_STR(NCTL).JQCTLD
            LD = LIJ(ID,JD)
            if( LU == 0 .and. LD == 0 ) GOTO 410
            FLOWX = HYD_STR(NCTL).RQCMUL*QCTLT(K,NCTL,1)
            UDDXTMP = FLOWX/DXP(LU)/DYP(LU)/(HPK(LU,K))
            if( IU == ID )then
              ADDLS = DXV(LU)*AH(LU,K)*DZC(LU,K)*0.5*(HP(LU) + HP(LD))*DYIV(LU)    ! *** DELME - PMC - The LU and LD don't seem correct
            else                                                                   ! *** DELME - PMC - The LU and LD don't seem correct
              ADDLS = DYU(LU)*AH(LU,K)*DZC(LU,K)*0.5*(HP(LU) + HP(LD))*DXIU(LU)    ! *** DELME - PMC - The LU and LD don't seem correct
            endif
            if( ISHDMF < 2 )then
              ADDLS = 0.
            endif
            LWASP = LWASP + 1
            FLOW(LWASP) = FLOWX
            CRNU(LWASP) = UDDXTMP
            BRINTT(LWASP) = ADDLS
410       enddo
        enddo

        ! ***  Advection and dispersion in the Z-direction:
        if( KC > 1 )then
          do K = KS,1,-1
            do LT = 2,LALT
              I = ILLT(LT)
              J = JLLT(LT)
              L = LIJ(I,J)
              addl = 0.0
              ADDL1 = ab(l,k)*hp(l)  !hnr Ev = ab*H
              if( addl1 > ABWMX(l) )then
                addl1 = ABWMX(l)
              endif
              if( SPB(L) == 1. )then
                ADDL = DXYP(L)*addl1/hp(l)*DZIG(L,K)
                vol1 = DXYP(L)*HP(L)*DZC(L,K)
                vol2 = DXYP(L)*HP(L)*DZC(L,K + 1)
                if( vol1 < vol2) vol2 = vol1
                AD = ADCOEFF*vol2/dt0
                if( addl > AD )then
                  addl = AD
                endif
              endif
              if( SWB(L) == 1 )then
                LWASP = LWASP + 1
                FLOW(LWASP) = -DXYP(L)*W(L,K)
                CRNU(LWASP) = W(L,K)*DZIG(L,K)/HP(L)
                BRINTT(LWASP) = ADDL
              endif
            enddo
          enddo
        endif

        call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        ! ***  next only IF temperature is transfered
        if( ISTRAN(2) >= 1 )then
          call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
          if( ierror  > 0 )then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            STOP
          endif
        endif

        ! ***  next only IF salinity is transfered
        if( ISTRAN(1) >= 1 )then
          call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
          if( ierror  > 0 )then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            STOP
          endif
        endif

        ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
        call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif

        call hlmomentcomplete(Ihl_Handle,ierror)
        if( ierror  > 0 )then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          STOP
        endif
        ! *** END OF HOTSTART

      endif  ! *** END OF ISRESTI == 0

    endif    ! ***  END INITIAL CONDITIONS WHEN IDAYS = 0

    ! *** FINISH INITIALIZATION OF FILES (JSWASP = 1)
    GOTO 3000

  endif   ! *** END OF JSWASP = 1

  ! ***************************************************************************
  IBEGIN = IDAYS*NTSPTC
  if( N < IBEGIN) GOTO 3000

  ! *** ----------------------------------------------------------------------!
  ! ***
  ! ***  WRITE TIME STEP DATA
  
  ! ***  Advection and dispersion in the X-direction:
  LWASP = 0
  do K = KC,1,-1
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      ADDLW = 0.0
      if( SUB(L) == 1. )then
        LW = LWC(L)
        ADDLW = DYU(L)*AHULPF(L,K)*DZC(L,K)*0.5*(HLPF(L) +HLPF(LW))*DXIU (L)
        vol1 = DXYP(L)*HLPF(L)*DZC(L,K)
        vol2 = DXYP(LW)*HLPF(LW)*DZC(L,K)
        if( vol1 < vol2) vol2 = vol1
        AD = ADCOEFF*vol2/dt0
        if( addlw > AD )then
          addlw = AD
        endif
        if( IShdmf < 2 )then
          addlw = 0.
        endif
      endif
      if( SUBO(L) == 1. )then
        TMPVAL = UHLPF(L,K) + SVPT*UVPT(L,K)
        !FLOWX = DYU(L)*TMPVAL*DZC(L,K)
        FLOWX = DYU(L)*TMPVAL          ! *** DKT removed DZC(L,K)
        UDDXTMP = 2.*TMPVAL*DXIU(L)/(HLPF(L) + HLPF(LWC(L)))
        LWASP = LWASP + 1
        ! ***        if( IJCTLT(I-1,J) == 8 ) addlw = 0.0
        FLOW(LWASP) = FLOWX
        CRNU(LWASP) = UDDXTMP
        BRINTT(LWASP) = ADDLW
      endif
      
      ! *** East open boundary cells
      if( IJCTLT(I + 1,J) == 8 )then
        LE = LEC(L)
        if( SUBO(LE) == 1. )then
          TMPVAL = UHLPF(LE,K) + SVPT*UVPT(LE,K)
          !FLOWX = DYU(LE)*TMPVAL*DZC(L,K)
          FLOWX = DYU(LE)*TMPVAL          ! *** DKT removed DZC(L,K)
          UDDXTMP = 2.*TMPVAL*DXIU(LE)/(HLPF(LE) + HLPF(L))
          IPTMP = I + 1
          ! ***          addlw = 0.0
          LWASP = LWASP + 1
          FLOW(LWASP) = FLOWX
          CRNU(LWASP) = UDDXTMP
          BRINTT(LWASP) = ADDLW
        endif
      endif
    enddo

    ! ***  Advection and dispersion in the Y-direction:
    do LT = 2,LALT
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      ADDLS = 0.0
      if( SVB(L) == 1. )then
        LS = LSC(L)
        ADDLS = DXV(L)*AHVLPF(L,K)*DZC(L,K)*0.5*(HLPF(L) +HLPF(LS))*DYIV (L)
        vol1 = DXYP(L)*HLPF(L)*DZC(L,K)
        vol2 = DXYP(LS)*HLPF(LS)*DZC(L,K)
        if( vol1 < vol2) vol2 = vol1
        AD = ADCOEFF*vol2/dt0
        if( addls > AD )then
          addls = AD
        endif
        if( IShdmf < 2 )then
          addls = 0.
        endif
      endif
      if( SVBO(L) == 1. )then
        TMPVAL = VHLPF(L,K) + SVPT*VVPT(L,K)
        !FLOWY = DXV(L)*TMPVAL*DZC(L,K)
        FLOWY = DXV(L)*TMPVAL          ! *** DKT removed DZC(L,K)
        ! ***        FLOWY = DXV(L)*vlpf(l,k)*hlpf(l)*DZC(L,K)
        VDDYTMP = 2.*TMPVAL*DYIV(L)/(HLPF(L) + HLPF(LSC(L)))
        JMTMP = J-1
        ! ***        if( IJCTLT(I,J-1) == 8 ) addls = 0.0
        LWASP = LWASP + 1
        FLOW(LWASP) = FLOWY
        CRNU(LWASP) = VDDYTMP
        BRINTT(LWASP) = ADDLS
      endif
      
      ! *** North open boundary cells
      if( IJCTLT(I,J + 1) == 8 )then
        LN = LNC(L)
        if( SVBO(LN) == 1. )then
          TMPVAL = VHLPF(LN,K) + SVPT*VVPT(LN,K)
          !FLOWY = DXV(LN)*TMPVAL*DZC(L,K)
          FLOWY = DXV(LN)*TMPVAL          ! *** DKT removed DZC(L,K)
          VDDYTMP = 2.*TMPVAL*DYIV(LN)/(HLPF(LN) + HLPF(L))
          JPTMP = J + 1
          ! ***          addls = 0.0
          LWASP = LWASP + 1
          FLOW(LWASP) = FLOWY
          CRNU(LWASP) = VDDYTMP
          BRINTT(LWASP) = ADDLS
        endif
      endif
    enddo    ! *** ENDI OF LT = 2,LALT LOOP
  enddo      ! *** END OF KC LOOP

  ! ***  Advection and dispersion in input flows
  do K = KC,1,-1
    do LT = 1,NQSIJ
      I = BCPS(LT).I
      J = BCPS(LT).J
      if( LIJLT(I,J) == 0 ) GOTO 300
      NS = BCPS(LT).NQSERQ
      L = BCPS(LT).L
      !IF(FLAGWASPBC(NS,K) == 1 ) GOTO 300
      if(K >= KSZ(L) )then
        FLOWX = BCPS(LT).RQSMUL*(QSS(K,LT) + BCPS(LT).QFACTOR*QSERT(K,NS))
        UDDXTMP = FLOWX/DXP(L)/DYP(L)/(HPK(L,K))
      else
        FLOWX = 0.
        UDDXTMP = 0.
      endif
      ! ***      ADDLW = DYP(L)*AHULPF(L,K)*DZC(L,K)*HLPF(L)/DXP(L)  !HNR
      ADDLW = 0.0
      LWASP = LWASP + 1
      FLOW(LWASP) = FLOWX
      CRNU(LWASP) = UDDXTMP
      BRINTT(LWASP) = ADDLW
300 END DO
  enddo

  ! ***   Advection and dispersion in structure flows
  do K = KC,1,-1
    do NCTL = 1,NQCTL
      IU = HYD_STR(NCTL).IQCTLU
      JU = HYD_STR(NCTL).JQCTLU
      LU = LIJ(IU,JU)
      ID = HYD_STR(NCTL).IQCTLD
      JD = HYD_STR(NCTL).JQCTLD
      LD = LIJ(ID,JD)
      if( LU == 0 .and. LD == 0 ) GOTO 400
      FLOWX = HYD_STR(NCTL).RQCMUL*QCTLT(K,NCTL,1)
      UDDXTMP = FLOWX/DXP(LU)/DYP(LU)*(HPKI(LU,K))
      if( IU == ID )then           ! delme - IU and ID?  Shouldn't these be LU and LD?
        ADDLS = DXV(LU)*AHVLPF(LU,K)*DZC(LU,K)*0.5*(HLPF(LU) + HLPF(LD))*DYIV(LU)    ! *** DELME - PMC - The LU and LD don't seem correct
      else
        ADDLS = DYU(LU)*AHULPF(LU,K)*DZC(LU,K)*0.5*(HLPF(LU) + HLPF(LD))*DXIU(LU)    ! *** DELME - PMC - The LU and LD don't seem correct
      endif
      if( ISHDMF < 2 )then
        ADDLS = 0.
      endif
      LWASP = LWASP + 1
      FLOW(LWASP) = FLOWX
      CRNU(LWASP) = UDDXTMP
      BRINTT(LWASP) = ADDLS
400 END DO
  enddo

  ! ***  Advection and dispersion in the Z-direction:
  if( KC > 1 )then
    do K = KS,1,-1
      do LT = 2,LALT
        I = ILLT(LT)
        J = JLLT(LT)
        L = LIJ(I,J)
        ADDL = 0.0
        ADDL1 = ABLPF(L,K)*HLPF(L)  !HNR  EV = AB*HP
        if( ADDL1 > ABWMX(L) )then
          ADDL1 = ABWMX(L)   !HNR
        endif
        if( SPB(L) == 1. )then
          ADDL = DXYP(L)*ADDL1/HLPF(L)*DZIG(L,K)
          VOL1 = DXYP(L)*HLPF(L)*DZC(L,K)
          VOL2 = DXYP(L)*HLPF(L)*DZC(L,K + 1)
          if( VOL1 < VOL2 ) VOL2 = VOL1
          AD = ADCOEFF*VOL2/DT0
          if( ADDL > AD )then
            ADDL = AD
          endif
        endif
        if( SWB(L) == 1 )then
          TMPVAL = WLPF(L,K) + SVPT*WVPT(L,K)
          FLOWZ = -DXYP(L)*TMPVAL
          WDDZTMP = TMPVAL*DZIG(L,K)/HLPF(L)
          LWASP = LWASP + 1
          FLOW(LWASP) = FLOWZ
          CRNU(LWASP) = WDDZTMP
          BRINTT(LWASP) = ADDL
        endif
      enddo
    enddo
  endif

  ! ***  Segment Properties:
  LCELTMP = 0
  do K = KC,1,-1
    do LT = 2,LALT
      LCELTMP = LCELTMP + 1
      I = ILLT(LT)
      J = JLLT(LT)
      L = LIJ(I,J)
      LN = LNC(L)
      LW = LWC(L)
      if( K >= KSZ(L) )then
        VOLUM = DXYP(L)*HLPF(L)*DZC(L,K)
        if( RESSTEP < TIDALP ) VOLUM = DXYP(L)*HP(L)*DZC(L,K)
        DEPTH = HLPF(L)*DZC(L,K)
        VELX = 0.5*(UHLPF(L,K) + SVPT*UVPT(L,K) + UHLPF(LW,K) + SVPT*UVPT(LW,K))/DEPTH
        VELY = 0.5*(VHLPF(L,K) + SVPT*VVPT(L,K) + VHLPF(LN,K) + SVPT*VVPT(LN,K) )/DEPTH
        VELZ = 0.5*(WLPF(L,K-1) + SVPT*WVPT(L,K-1) + WLPF(L,K) + SVPT*WVPT(L,K) )
        VELMAG = SQRT(VELX*VELX + VELY*VELY + VELZ*VELZ)
      else
        VOLUM = 0.
        DEPTH = 0.
        VELMAG  = 0.
      endif
      SegSalt(LCELTMP) = SALLPF(L,K)
      SegTemp(LCELTMP) = TEMLPF(L,K)
      SegVolume(LCELTMP) = VOLUM
      SegDepth(LCELTMP) = DEPTH
      SegVel(LCELTMP) = VELMAG
    enddo
  enddo

  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

  ! ***  next only IF temperature is transfered
  if( ISTRAN(2) >= 1 )then
    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
  endif

  ! ***  next only IF salinity is transfered
  if( ISTRAN(1) >= 1 )then
    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
    if( ierror  > 0 )then
      call hlgetlasterror(errstring)
      write(6,6000) ierror, errstring
      STOP
    endif
  endif

  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  = 
  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

  call hlmomentcomplete(Ihl_Handle,ierror)
  if( ierror  > 0 )then
    call hlgetlasterror(errstring)
    write(6,6000) ierror, errstring
    STOP
  endif

3000 JSWASP = 0

6000 format('Error ',I10, ' : ', A)
9966 FORMAT(f12.4,I6,f12.4)
9977 FORMAT(3I6,6f12.4)

  ! *** CLOSE LINKAGE FILE
  if( N >= NTS )then
    call hlclose(IHL_HANDLE,ierror)
    IHL_HANDLE = 0
  endif
#endif

  ! *** ----------------------------------------------------------------------!
  ! *** END SUBROUTINE WASPHYDROLINK
  ! *** ----------------------------------------------------------------------!
  return

  END
