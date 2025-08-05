! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE FOODCHAIN(IFINISH)

  ! ***  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
  !
  ! ***  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001

  !----------------------------------------------------------------------C
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY                 DATE APPROVED    BY

  !----------------------------------------------------------------------C
  !
  ! ***  SUBROUTINES OUTPUT SPACE AND TIME AVERAGE TOXICS CONCENTRATIONS
  ! ***  FOR FOOD CHAIN MODEL

  ! *** *******************************************************************C

  use GLOBAL
  implicit none

  integer :: IFINISH,JSFDCH,M,NT,L,K,NS,NX,KSTOP,KTMP
  real    :: TIMFDCH,HBEDTMP,HBSTOP,TMPVAL,PORHINV,FDCHVAL
  
  integer,allocatable,dimension(:) :: KBFC

  real,allocatable,dimension(:,:) :: VALPOCW
  real,allocatable,dimension(:,:) :: TMPVOLW

  real,allocatable,dimension(:,:) :: WTBED
  real,allocatable,dimension(:,:) :: VALPOCB
  real,allocatable,dimension(:,:) :: TMPVOLB

  real,allocatable,dimension(:,:) :: TMPTXWF
  real,allocatable,dimension(:,:) :: TMPTXWC
  real,allocatable,dimension(:,:) :: TMPTXWP
  real,allocatable,dimension(:,:) :: TMPTXBF
  real,allocatable,dimension(:,:) :: TMPTXBC
  real,allocatable,dimension(:,:) :: TMPTXBP
  real,allocatable,dimension(:,:) :: TMPTXBPD

  real,allocatable,dimension(:,:) :: VALBCONC

  real,allocatable,dimension(:) :: TMPDOCW
  real,allocatable,dimension(:) :: TMPPOCW
  real,allocatable,dimension(:) :: TMPDOCB
  real,allocatable,dimension(:) :: TMPPOCB
  real,allocatable,dimension(:) :: VOLFCW
  real,allocatable,dimension(:) :: VOLFCB

  logical,allocatable,dimension(:) :: LMASKFC

  real,allocatable,dimension(:) :: FDCHDOCW
  real,allocatable,dimension(:) :: FDCHPOCW
  real,allocatable,dimension(:) :: FDCHDOCB
  real,allocatable,dimension(:) :: FDCHPOCB

  real,allocatable,dimension(:,:) :: FDCHTXWF
  real,allocatable,dimension(:,:) :: FDCHTXWC
  real,allocatable,dimension(:,:) :: FDCHTXWP
  real,allocatable,dimension(:,:) :: FDCHTXBF
  real,allocatable,dimension(:,:) :: FDCHTXBC
  real,allocatable,dimension(:,:) :: FDCHTXBP
  real,allocatable,dimension(:,:) :: FDCHTXBD

  if( .not. allocated(KBFC) )then
    allocate(KBFC(LCM))
    allocate(VALPOCW(LCM,KCM))
    allocate(TMPVOLW(LCM,KCM))

    allocate(WTBED(LCM,KBM))
    allocate(VALPOCB(LCM,KBM))
    allocate(VALBCONC(LCM,KBM))

    allocate(TMPVOLB(LCM,KBM))
    allocate(TMPTXWF(NFDCHZ,NTXM))
    allocate(TMPTXWC(NFDCHZ,NTXM))
    allocate(TMPTXWP(NFDCHZ,NTXM))
    allocate(TMPTXBF(NFDCHZ,NTXM))
    allocate(TMPTXBC(NFDCHZ,NTXM))
    allocate(TMPTXBP(NFDCHZ,NTXM))
    allocate(TMPTXBPD(NFDCHZ,NTXM))

    allocate(TMPDOCW(NFDCHZ))
    allocate(TMPPOCW(NFDCHZ))
    allocate(TMPDOCB(NFDCHZ))
    allocate(TMPPOCB(NFDCHZ))
    allocate(VOLFCW(NFDCHZ))
    allocate(VOLFCB(NFDCHZ))

    allocate(FDCHDOCW(NFDCHZ))
    allocate(FDCHPOCW(NFDCHZ))
    allocate(FDCHDOCB(NFDCHZ))
    allocate(FDCHPOCB(NFDCHZ))

    allocate(FDCHTXWF(NFDCHZ,NTXM))
    allocate(FDCHTXWC(NFDCHZ,NTXM))
    allocate(FDCHTXWP(NFDCHZ,NTXM))
    allocate(FDCHTXBF(NFDCHZ,NTXM))
    allocate(FDCHTXBC(NFDCHZ,NTXM))
    allocate(FDCHTXBP(NFDCHZ,NTXM))
    allocate(FDCHTXBD(NFDCHZ,NTXM))

    ! *** ALLOCATE LOCAL ARRAYS
    KBFC = 0
    VALPOCW = 0.
    TMPVOLW = 0.

    WTBED = 0.
    VALPOCB = 0.
    VALBCONC = 0.

    TMPVOLB = 0.
    TMPTXWF = 0.
    TMPTXWC = 0.
    TMPTXWP = 0.
    TMPTXBF = 0.
    TMPTXBC = 0.
    TMPTXBP = 0.
    TMPTXBPD = 0.

    TMPDOCW = 0.
    TMPPOCW = 0.
    TMPDOCB = 0.
    TMPPOCB = 0.
    VOLFCW = 0.
    VOLFCB = 0.

    FDCHDOCW = 0.
    FDCHPOCW = 0.
    FDCHDOCB = 0.
    FDCHPOCB = 0.

    FDCHTXWF = 0.
    FDCHTXWC = 0.
    FDCHTXWP = 0.
    FDCHTXBF = 0.
    FDCHTXBC = 0.
    FDCHTXBP = 0.
    FDCHTXBD = 0.
  endif

  ! *** *******************************************************************C

  if( IFINISH == 1 ) GO TO 2000
  if( JSFDCH == 0 )GO TO 1000

  !      write(mpi_efdc_out_unit,*)' FIRST ENTRY TO FOODCHAIN.FOR '

  if( DEBUG )then
    open(1,FILE = OUTDIR//'FOODCHAIN.OUT')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'FOODCHAIN.OUT')
    write(1,121)
    write(1,122)
    write(1,123)
    close(1)
  endif

  !     JSFDCH = 0

  do M = 1,NFDCHZ
    FDCHDOCW(M) = 0.
    FDCHPOCW(M) = 0.
    FDCHDOCB(M) = 0.
    FDCHPOCB(M) = 0.
  enddo

  do NT = 1,NTOX
    do M = 1,NFDCHZ
      FDCHTXWF(M,NT) = 0.
      FDCHTXWC(M,NT) = 0.
      FDCHTXWP(M,NT) = 0.
      FDCHTXBF(M,NT) = 0.
      FDCHTXBC(M,NT) = 0.
      FDCHTXBP(M,NT) = 0.
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      FDCHTXBD(M,NT) = 0.
    !####################################################################################
    enddo
  enddo

  TIMFDCH = 0.0

  ! *** *******************************************************************C

    1000 continue

  TIMFDCH = TIMFDCH + DTSED

  ! ***  INITIALIZE VOLUMES AND VOLUME AVERAGES

  do M = 1,NFDCHZ
    VOLFCW(M) = 0.
    VOLFCB(M) = 0.
  enddo

  do M = 1,NFDCHZ
    TMPDOCW(M) = 0.
    TMPPOCW(M) = 0.
    TMPDOCB(M) = 0.
    TMPPOCB(M) = 0.
  enddo

  do NT = 1,NTOX
    do M = 1,NFDCHZ
      TMPTXWF(M,NT) = 0.
      TMPTXWC(M,NT) = 0.
      TMPTXWP(M,NT) = 0.
      TMPTXBF(M,NT) = 0.
      TMPTXBC(M,NT) = 0.
      TMPTXBP(M,NT) = 0.
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      TMPTXBPD(M,NT) = 0.
    !####################################################################################
    enddo
  enddo

  ! ***  INITIALIZE MASK
  do L = 2,LA
    LMASKFC(L) = .FALSE.
  enddo
  !
  do L = 2,LA
    if( LMASKDRY(L) )then
      if( MFDCHZ(L) > 0 )LMASKFC(L) = .TRUE.
    endif
  enddo

  !----------------------------------------------------------------------C
  !
  ! ***  VOLUME WEIGHTED AVERAGE OVER WATER COLUMN ZONES
  !
  !     STDOCW(L,K) HAS UNITS: MG/L OR GM/M**3
  !     STPOCW(L,K) AND VALPOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
  if( ISTPOCW <= 1 )then
    do K = 1,KC
      do L = 2,LA
        if( LMASKFC(L) ) VALPOCW(L,K) = STPOCW(L,K)
      enddo
    enddo
  endif

  if( ISTPOCW >= 2 )then
    do K = 1,KC
      do L = 2,LA
        if( LMASKFC(L) ) VALPOCW(L,K) = 0.
      enddo
    enddo
    do NS = 1,NSED
      do K = 1,KC
        do L = 2,LA
          if( LMASKFC(L) ) VALPOCW(L,K) = VALPOCW(L,K) + SED(L,K,NS)*STFPOCW(L,K,NS)
        enddo
      enddo
    enddo
    do NX = 1,NSND
      NS = NX+NSED
      do K = 1,KC
        do L = 2,LA
          if( LMASKFC(L)) VALPOCW(L,K) = VALPOCW(L,K) + SND(L,K,NX)*STFPOCW(L,K,NS)
        enddo
      enddo
    enddo
  endif

  !     changed to areal weighting from voloume weighting
  do K = 1,KC
    do L = 2,LA
      if( LMASKFC(L)) TMPVOLW(L,K) = DXYP(L)*DZC(L,K)
    enddo
  enddo

  do K = 1,KC
    do L = 2,LA
      if( LMASKFC(L) )then
        M = MFDCHZ(L)
        VOLFCW(M) = VOLFCW(M)+TMPVOLW(L,K)
      endif
    enddo
  enddo

  !     TMPTXWF(M,NT)/TMPVOLW HAS UNITS: UG/L OR MG/M**3
  !     TMPTXWC(M,NT)/TMPVOLW HAS UNITS: UG/L OR MG/M**3
  !     TMPTXWP(M,NT)/TMPVOLW HAS UNITS: UG/MG
  !     TMPDOCW(M,NT)/TMPVOLW AND STDOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
  !     TMPPOCW(M,NT)/TMPVOLW AND VALPOCW(L,K) HAVE UNITS: MG/L OR GM/M**3

  do K = 1,KC
  do L = 2,LA
    if( LMASKFC(L) )then
    M = MFDCHZ(L)
    TMPDOCW(M) = TMPDOCW(M)+TMPVOLW(L,K)*STDOCW(L,K)
    TMPPOCW(M) = TMPPOCW(M)+TMPVOLW(L,K)*VALPOCW(L,K)
    endif
  enddo
  enddo

  do NT = 1,NTOX
  do K = 1,KC
  do L = 2,LA
    if( LMASKFC(L) )then
      M = MFDCHZ(L)
      TMPTXWF(M,NT) = TMPTXWF(M,NT) + TMPVOLW(L,K)*TOXFDFW(L,K,NT)*TOX(L,K,NT)
      TMPTXWC(M,NT) = TMPTXWC(M,NT) + TMPVOLW(L,K)*TOXCDFW(L,K,NT)*TOX(L,K,NT)
      if( VALPOCW(L,K) > 0.) TMPTXWP(M,NT) = TMPTXWP(M,NT) + TMPVOLW(L,K)*TOXPFTW(L,K,NT)*TOX(L,K,NT) / VALPOCW(L,K)
    endif
  enddo
  enddo
  enddo

  do M = 1,NFDCHZ
    if( VOLFCW(M) > 0.0 )then
      TMPDOCW(M) = TMPDOCW(M)/VOLFCW(M)
      TMPPOCW(M) = TMPPOCW(M)/VOLFCW(M)
    endif
  enddo


  do NT = 1,NTOX
    do M = 1,NFDCHZ
      if( VOLFCW(M) > 0.0 )then
        TMPTXWF(M,NT) = TMPTXWF(M,NT)/VOLFCW(M)
        TMPTXWC(M,NT) = TMPTXWC(M,NT)/VOLFCW(M)
        TMPTXWP(M,NT) = TMPTXWP(M,NT)/VOLFCW(M)
      endif
    enddo
  enddo

  !     CONVERT PARTICULATE FROM UG/MG TO UG/GM

  do NT = 1,NTOX
    do M = 1,NFDCHZ
      TMPTXWP(M,NT) = 1000.*TMPTXWP(M,NT)
    enddo
  enddo

  !----------------------------------------------------------------------C
  !
  ! ***  VOLUME WEIGHTED AVERAGE OVER BED ZONES
  !
  !     STDOCB(L,K) HAS UNITS: MG/L OR GM/M**3 (MASS PER VOLUME OF PORE WATER)
  !     STPOCB(L,K) AND VALPOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER TOTAL VOLUME)

  if( ISTPOCB <= 1 )then
    do K = 1,KB
      do L = 2,LA
        if( LMASKFC(L)) VALPOCB(L,K) = STPOCB(L,K)
      enddo
    enddo
  endif

  !      if( ISTPOCB >= 2 )then   ! PMC  FPOCB IS NOT INIITALIZED UNTIL ISTPOCB>3
  if( ISTPOCB >= 4 )then
    do K = 1,KB
      do L = 2,LA
        VALPOCB(L,K) = 0.
      enddo
    enddo
    do NS = 1,NSED
      do K = 1,KB
        do L = 2,LA
          if( LMASKFC(L) )then
            if( K <= KBT(L)) VALPOCB(L,K) = VALPOCB(L,K)+SEDB(L,K,NS)*FPOCB(L,K)/HBED(L,K)
            !####################################################################################
            ! RM 05/14/04
            ! Change to average using data-based foc rather than partitioning foc
            !     &                    +SEDB(L,K,NS)*STFPOCB(L,K,NS)/HBED(L,K)
            !####################################################################################
          endif
        enddo
      enddo
    enddo
    do NX = 1,NSND
      NS = NX+NSED
      do K = 1,KB
        do L = 2,LA
          if( LMASKFC(L) )then
            if( K <= KBT(L)) VALPOCB(L,K) = VALPOCB(L,K)+SNDB(L,K,NX)*FPOCB(L,K)/HBED(L,K)
            !####################################################################################
            ! RM 05/14/04
            ! Change to average using data-based foc rather than partitioning foc
            !     &                    +SNDB(L,K,NX)*STFPOCB(L,K,NS)/HBED(L,K)
            !####################################################################################
          endif
        enddo
      enddo
    enddo
  endif

  !####################################################################################
  ! RM 05/14/04
  ! Change to average dry weight PCB
  do K = 1,KB
    do L = 2,LA
        VALBCONC(L,K) = 0.
    enddo
  enddo
  do NS = 1,NSED
    do K = 1,KB
      do L = 2,LA
        if( LMASKFC(L) )then
          if( K <= KBT(L)) VALBCONC(L,K) = VALBCONC(L,K) + SEDB(L,K,NS)
        endif
      enddo
    enddo
  enddo
  do NX = 1,NSND
    NS = NX+NSED
    do K = 1,KB
      do L = 2,LA
        if( LMASKFC(L) )then
          if( K <= KBT(L)) VALBCONC(L,K) = VALBCONC(L,K) + SNDB(L,K,NX)
        endif
      enddo
    enddo
  enddo
  !####################################################################################

  do K = 1,KB
    do L = 2,LA
      WTBED(L,K) = 0.0
    enddo
  enddo

  do L = 2,LA
    if( LMASKFC(L) )then
      KBFC(L) = KBT(L)
      HBEDTMP = 0.0
      KSTOP = 0
      do K = KBT(L),1,-1
        HBEDTMP = HBEDTMP+HBED(L,K)
        if( HBEDTMP > HBFDCH .and. KSTOP == 0 )then
          KBFC(L) = K
          KSTOP = 1
          HBSTOP = HBED(L,K)-HBEDTMP+HBFDCH
          WTBED(L,K) = HBSTOP/HBFDCH
        endif
      enddo
      KTMP = KBFC(L)+1
      do K = KTMP,KBT(L)
        !####################################################################################
        ! RM 05/14/04
        ! Weightages greater than 1 could occur with this method of depth-weighting.
        ! When the thickness of the top layer is greater than HBFDCH (0.1524
        ! meters), the weightage assigned to this layer could become greater
        ! than 1. Need to confirm this with JH.
        !####################################################################################
        WTBED(L,K) = HBED(L,K)/HBFDCH
      enddo
    else
      KBFC(L) = 0
    endif
  enddo
  !
  if( JSFDCH == 1 .and. DEBUG )then
    open(1,FILE = OUTDIR//'FOODCHAIN.DIA')
    close(1,STATUS = 'DELETE')
    open(1,FILE = OUTDIR//'FOODCHAIN.DIA')
    do L = 2,LA
      if( LMASKFC(L) )then
        write(1,111)IL(L),JL(L),KBFC(L),KBT(L),(WTBED(L,K),K = KBFC(L),KBT(L))
        write(1,112)(HBED(L,K),K = KBFC(L),KBT(L))
      endif
    enddo
    close(1)
  endif
  
  do K = 1,KB
  do L = 2,LA
    if( LMASKFC(L) )then
      if( K >= KBFC(L) .and. K <= KBT(L) )then
        TMPVOLB(L,K) = DXYP(L)*WTBED(L,K)*HBED(L,K)
      endif
    endif
  enddo
  enddo

  do K = 1,KB
  do L = 2,LA
    if( LMASKFC(L) )then
      M = MFDCHZ(L)
      if( K >= KBFC(L) .and. K <= KBT(L) )then
        VOLFCB(M) = VOLFCB(M)+TMPVOLB(L,K)
      endif
    endif
  enddo
  enddo

  !     TMPTXBF(M,NT)/TMPVOLB HAS UNITS: UG/L OR MG/M**3  (MASS PER VOLUME PORE WATER)
  !     TMPTXBC(M,NT)/TMPVOLB HAS UNITS: UG/L OR MG/M**3  (MASS PER VOLUME PORE WATER)
  !     TMPTXWP(M,NT)/TMPVOLB HAS UNITS: UG/MG
  !     TMPDOCW(M,NT)/TMPVOLB AND STDOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER VOLUME PORE WATER)
  !     TMPPOCW(M,NT)/TMPVOLB AND VALPOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER TOTAL VOLUME)

  do K = KB,1,-1
  do L = 2,LA
    if( LMASKFC(L) )then
      M = MFDCHZ(L)
      if( K >= KBFC(L) .and. K <= KBT(L) )then
        TMPDOCB(M) = TMPDOCB(M)+TMPVOLB(L,K)*STDOCB(L,K)
        TMPPOCB(M) = TMPPOCB(M)+TMPVOLB(L,K)*VALPOCB(L,K)
      endif
    endif
  enddo
  enddo

  do NT = 1,NTOX
  do K = KB,1,-1
  do L = 2,LA
    if( LMASKFC(L) )then
      M = MFDCHZ(L)
      if( K >= KBFC(L) .and. K <= KBT(L) )then
        TMPVAL = HBED(L,K)*VALPOCB(L,K)
        PORHINV = 1.0/(HBED(L,K)*PORBED(L,K))
        TMPTXBF(M,NT) = TMPTXBF(M,NT)+TMPVOLB(L,K)*PORHINV*TOXFDFB(L,K,NT)*TOXB(L,K,NT)
        TMPTXBC(M,NT) = TMPTXBC(M,NT)+TMPVOLB(L,K)*PORHINV*TOXCDFB(L,K,NT)*TOXB(L,K,NT)
        
        if( TMPVAL > 0.) TMPTXBP(M,NT) = TMPTXBP(M,NT)+TMPVOLB(L,K)*TOXPFTB(L,K,NT)*TOXB(L,K,NT)/TMPVAL
        !####################################################################################
        ! RM 05/14/04
        ! Change to average dry weight PCBs
        TMPTXBPD(M,NT) = TMPTXBPD(M,NT) + TMPVOLB(L,K)*TOXPFTB(L,K,NT)*TOXB(L,K,NT)/VALBCONC(L,K)
        !####################################################################################
      endif
    endif
  enddo
  enddo
  enddo
  
  do M = 1,NFDCHZ
    if( VOLFCB(M) > 0.0 )then
      TMPDOCB(M) = TMPDOCB(M)/VOLFCB(M)
      TMPPOCB(M) = TMPPOCB(M)/VOLFCB(M)
    endif
  enddo
  !
  do NT = 1,NTOX
  do M = 1,NFDCHZ
    if( VOLFCB(M) > 0.0 )then
      TMPTXBF(M,NT) = TMPTXBF(M,NT)/VOLFCB(M)
      TMPTXBC(M,NT) = TMPTXBC(M,NT)/VOLFCB(M)
      TMPTXBP(M,NT) = TMPTXBP(M,NT)/VOLFCB(M)
      !####################################################################################
      ! RM 05/14/04
      ! Change to average dry weight PCBs
          TMPTXBPD(M,NT) = TMPTXBPD(M,NT)/VOLFCB(M)
      !####################################################################################
    endif
  enddo
  enddo

  !     CONVERT PARTICULATE FROM UG/MG TO UG/GM
  do NT = 1,NTOX
    do M = 1,NFDCHZ
      TMPTXBP(M,NT) = 1000.*TMPTXBP(M,NT)
    enddo
  enddo

  !----------------------------------------------------------------------C
  !
  ! ***  ACCUMULATE THE TIME AVERAGE
  do M = 1,NFDCHZ
    FDCHDOCW(M) = FDCHDOCW(M)+DTSED*TMPDOCW(M)
    FDCHPOCW(M) = FDCHPOCW(M)+DTSED*TMPPOCW(M)
    FDCHDOCB(M) = FDCHDOCB(M)+DTSED*TMPDOCB(M)
    FDCHPOCB(M) = FDCHPOCB(M)+DTSED*TMPPOCB(M)
  enddo
  !
  do NT = 1,NTOX
    do M = 1,NFDCHZ
      FDCHTXWF(M,NT) = FDCHTXWF(M,NT)+DTSED*TMPTXWF(M,NT)
      FDCHTXWC(M,NT) = FDCHTXWC(M,NT)+DTSED*TMPTXWC(M,NT)
      FDCHTXWP(M,NT) = FDCHTXWP(M,NT)+DTSED*TMPTXWP(M,NT)
      FDCHTXBF(M,NT) = FDCHTXBF(M,NT)+DTSED*TMPTXBF(M,NT)
      FDCHTXBC(M,NT) = FDCHTXBC(M,NT)+DTSED*TMPTXBC(M,NT)
      FDCHTXBP(M,NT) = FDCHTXBP(M,NT)+DTSED*TMPTXBP(M,NT)
      !####################################################################################
      ! RM 05/14/04
      ! Change to average dry weight PCBs
      FDCHTXBD(M,NT) = FDCHTXBD(M,NT)+DTSED*TMPTXBPD(M,NT)
      !####################################################################################
      enddo
  enddo

  JSFDCH = 0

  if( TIMFDCH < TFCAVG) return

  ! *** *******************************************************************C
  !
  ! ***  COMPLETE AVERAGING AND OUTPUT RESULTS
  !
    2000 continue

  FDCHVAL = 1./TIMFDCH
  do M = 1,NFDCHZ
    FDCHDOCW(M) = FDCHVAL*FDCHDOCW(M)
    FDCHPOCW(M) = FDCHVAL*FDCHPOCW(M)
    FDCHDOCB(M) = FDCHVAL*FDCHDOCB(M)
    FDCHPOCB(M) = FDCHVAL*FDCHPOCB(M)
  enddo
  !
  do NT = 1,NTOX
    do M = 1,NFDCHZ
      FDCHTXWF(M,NT) = FDCHVAL*FDCHTXWF(M,NT)
      FDCHTXWC(M,NT) = FDCHVAL*FDCHTXWC(M,NT)
      FDCHTXWP(M,NT) = FDCHVAL*FDCHTXWP(M,NT)
      FDCHTXBF(M,NT) = FDCHVAL*FDCHTXBF(M,NT)
      FDCHTXBC(M,NT) = FDCHVAL*FDCHTXBC(M,NT)
      FDCHTXBP(M,NT) = FDCHVAL*FDCHTXBP(M,NT)
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      FDCHTXBD(M,NT) = FDCHVAL*FDCHTXBD(M,NT)
    !####################################################################################
    enddo
  enddo

  if( DEBUG )then
    open(1,FILE = OUTDIR//'FOODCHAIN.OUT',POSITION = 'APPEND')
    write(1,101)TIMEDAY,NTOX,NFDCHZ,TIMFDCH
    do NT = 1,NTOX
      do M = 1,NFDCHZ
        write(1,102)NT,M, &
                      FDCHTXWF(M,NT),FDCHTXWC(M,NT),FDCHTXWP(M,NT), &
                      FDCHDOCW(M),FDCHPOCW(M),FDCHTXBF(M,NT), &
                      FDCHTXBC(M,NT),FDCHTXBP(M,NT),FDCHDOCB(M), &
                      FDCHPOCB(M),FDCHTXBD(M,NT)
      enddo
    enddo
    close(1)
  endif

  ! *** *******************************************************************C
  !
  ! ***  INITIALIZE FOR NEXT AVERAGING PERIOD
  do M = 1,NFDCHZ
    FDCHDOCW(M) = 0.
    FDCHPOCW(M) = 0.
    FDCHDOCB(M) = 0.
    FDCHPOCB(M) = 0.
  enddo
  !
  do NT = 1,NTOX
    do M = 1,NFDCHZ
      FDCHTXWF(M,NT) = 0.
      FDCHTXWC(M,NT) = 0.
      FDCHTXWP(M,NT) = 0.
      FDCHTXBF(M,NT) = 0.
      FDCHTXBC(M,NT) = 0.
      FDCHTXBP(M,NT) = 0.
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      FDCHTXBD(M,NT) = 0.
    !####################################################################################
    enddo
  enddo

  TIMFDCH = 0.0

  ! *** *******************************************************************C
    111 FORMAT(4I5,10F10.4)
    112 FORMAT(20X,10F10.4)
    101 FORMAT(F12.4,2I7,F12.3)
    102 FORMAT(1X,2I6,10E13.5)
    103 FORMAT('              TXWF         TXWC         TXWP', &
          '         DOCW         POCW         TXBF         TXBC', &
          '         TXBP (roc)   DOCB         POCB        TXBPD (r)')
    121 FORMAT('DATA: OUTPUT TIME (DAYS), NTOX, NZONES, ', &
          'AERAGING PERIOD (SECS)')
    122 FORMAT('DATA: NT    NZ   TXWF         TXWC         TXWP', &
          '         DOCW         POCW         TXBF         TXBC', &
          '         TXBP (roc)   DOCB         POCB        TXBPD (r)')
    123 FORMAT('DATA:            UG/L         UG/L         UG/GM', &
          '        MG/L         MG/L         UG/L         UG/L', &
          '         UG/GM OC     MG/L         MG/L        UG/GM Dry')
  !
  ! *** *******************************************************************C
  !
  return

END
