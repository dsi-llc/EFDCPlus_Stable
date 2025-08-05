! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE AINIT

  ! CHANGE RECORD
  !
  !  ALL ZEROING OF ARRAYS MOVED TO VARZEROINT AND VARZEROREAL
  
  use GLOBAL
  use Variables_WQ
  
  use Variables_MPI
  use turbulence, only: eps_min,k_min
  
  implicit none  
                                                                                                        
  integer :: L,I,J,LS,NT,LCHNV,IVAL,NS,K,NMD,LHOST,LCHNU,NV,NX,LW                                                       
  integer :: NTMPC,NTMPN                                                                                                   

  ! ***  INITIALIZE ARRAYS                                                                                                 
  ZBR(1) = ZBRADJ
  ZBRE(1) = ZBRADJ
  HMP(1) = HMIN
  HMU(1) = HMIN
  HMV(1) = HMIN
  HWQ(1) = HMIN
  H2WQ(1) = HMIN
  DXP(1) = DX
  DYP(1) = DY
  DXU(1) = DX
  DYU(1) = DY
  DXV(1) = DX
  DYV(1) = DY
  DXYP(1) = DX*DY
  MVEGL(1) = 1
  BELV(1) = BELV(2)
  
  ZBR(LC) = ZBRADJ
  ZBRE(LC) = ZBRADJ
  HMP(LC) = HMIN
  HMU(LC) = HMIN
  HMV(LC) = HMIN
  HWQ(LC) = HMIN
  H2WQ(LC) = HMIN
  DXP(LC) = DX
  DYP(LC) = DY
  DXU(LC) = DX
  DYU(LC) = DY
  DXV(LC) = DX
  DYV(LC) = DY
  DXYP(LC) = DX*DY
  MVEGL(LC) = 1
  BELV(LC) = BELV(LA)
  BELV0 = BELV
  
  if( ISGWIE == 0 ) DAGWZ = 0.
  do L = 2,LA
    I = IL(L)
    J = JL(L)
    KBT(L) = 1
    BELAGW(L) = BELV(L) - DAGWZ
    ZBRE(L) = ZBR(L)
    CUE(L) = 1.
    CVE(L) = 0.
    CUN(L) = 0.
    CVN(L) = 1.
  enddo
  KBT(1) = 1
  KBT(LC) = 1
  
  do L = 2,LA
    LW = LWC(L)
    LS = LSC(L)
    DXU(L) = 0.5*(DXP(L)+DXP(LW))
    DYU(L) = 0.5*(DYP(L)+DYP(LW))
    DXV(L) = 0.5*(DXP(L)+DXP(LS))
    DYV(L) = 0.5*(DYP(L)+DYP(LS))
  enddo

  do L = 2,LA
    LW = LWC(L)
    LS = LSC(L)
    HMU(L) = 0.5*(DXP(L)*DYP(L)*HMP(L)+DXP(LW)*DYP(LW)*HMP(LW))/(DXU(L)*DYU(L))
    HMV(L) = 0.5*(DXP(L)*DYP(L)*HMP(L)+DXP(LS)*DYP(LS)*HMP(LS))/(DXV(L)*DYV(L))
  enddo

  HMU(1) = HMU(2)    ! *** PMC
  HMV(1) = HMV(2)    ! *** PMC
  HMU(LC) = HMU(LA)  ! *** PMC
  HMV(LC) = HMV(LA)  ! *** PMC
  do L = 1,LC
    CC(L) = 1.
    CCC(L) = 1.
    P(L) = G*(HMP(L)+BELV(L))
    P1(L) = G*(HMP(L)+BELV(L))
    HP(L) = HMP(L)
    HU(L) = HMU(L)
    HV(L) = HMV(L)
    HPI(L) = 1./HP(L)
    HUI(L) = 1./HU(L)
    HVI(L) = 1./HV(L)
    HWQ(L) = HMP(L)
    H1P(L) = HMP(L)
    H1U(L) = HMU(L)
    H1V(L) = HMV(L)
    H1UI(L) = 1./H1U(L)
    H1VI(L) = 1./H1V(L)
    H2WQ(L) = HMP(L)
    SCB(L) = 1.
    SPB(L) = 1.
    SUB(L) = 1.
    SVB(L) = 1.
    SWB(L) = 1.
    STCUV(L) = 1.
    STCAP(L) = 1.
    STBX(L) = 1.
    STBY(L) = 1.
    SAAX(L) = 1.
    SAAY(L) = 1.
    SCAX(L) = 1.
    SCAY(L) = 1.
    SBX(L) = 1.
    SBY(L) = 1.
    SDX(L) = 1.
    SDY(L) = 1.
    LMASKDRY(L) = .TRUE.
  enddo
  
  LOPENBCDRY = .FALSE.
  RADKE = SWRATNF
  
  ! *** OPEN WATER DEFAULT SETTINGS
  if( ISVEG > 0 )then
    NV = 0
    PVEGZ(NV)  = 1.
  endif
  
  !      if( IS1DCHAN == 1 )then
  !        do L = 1,LC
  !          FADYP(L) = 1.
  !          FADYP1(L) = 1.
  !          FADYP2(L) = 1.
  !          WPDYP(L) = 1.
  !          WPDYP1(L) = 1.
  !          FADXP(L) = 1.
  !          FADXP1(L) = 1.
  !          FADXP2(L) = 1.
  !          WPDXP(L) = 1.
  !          WPDXP1(L) = 1.
  !          FADYU(L) = 1.
  !          FADYU1(L) = 1.
  !          WPDYU(L) = 1.
  !          WPDYU1(L) = 1.
  !          FADXV(L) = 1.
  !          FADXV1(L) = 1.
  !          WPDXV(L) = 1.
  !          WPDXV1(L) = 1.
  !          DADH(L) = 1.
  !          DADH1(L) = 1.
  !          SRFXP(L) = 0.
  !          SRFYP(L) = 0.  C
  !          SRFXP1(L) = 0.
  !          SRFYP1(L) = 0.
  !          SRFXV(L) = 0.
  !          SRFYU(L) = 0.
  !          SRFXV1(L) = 0.
  !          SRFYU1(L) = 0.
  !        enddo
  !      endif
  do K = 1,KS
    do L = 1,LC
      AV(L,K) = AVO
      AVVI(L,K) = 1./AVO
      AVUI(L,K) = 1./AVO
      AB(L,K) = ABO
      QQL(L,K) = QQLMIN
      QQL1(L,K) = QQLMIN
      QQL2(L,K) = QQLMIN
      DML(L,K) = DMLMIN
      ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
    enddo
  enddo
  do K = 1,KC
    do L = 1,LC
      AH(L,K) = AHO
      AHC(L,K) = AHO
      AQ(L,K) = AVO
      ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      CTURBB1(L,K) = CTURB
      CTURBB2(L,K) = CTURB2B
      ! *** TEMPERATURE INITIATION
      TEM(L,K) = TEMO
      TEM1(L,K) = TEMO
    enddo
  enddo
  if( ISWQFLUX == 1 )then
    do K = 1,KC
      do L = 1,LC
        AHULPF(L,K) = AHO
        AHVLPF(L,K) = AHO
      enddo
    enddo
  endif
  
  NTMPC = max(NSED,1)
  do NS = 1,NTMPC
    do K = 1,KC
      do L = 1,LC
        SED(L,K,NS) = SEDO(NS)
        SED1(L,K,NS) = SEDO(NS)
        ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      enddo
    enddo
  enddo

  NTMPN = max(NSND,1)
  do NX = 1,NTMPN
    NS = NX+NTMPC
    do K = 1,KC
      do L = 1,LC
        SND(L,K,NX) = SEDO(NS)
        SND1(L,K,NX) = SEDO(NS)
        ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      enddo
    enddo
  enddo
  do NT = 1,NTOX
    do K = 1,KC
      do L = 1,LC
        TOX(L,K,NT) = TOXINTW(NT)
        TOX1(L,K,NT) = TOXINTW(NT)
        ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      enddo
    enddo
  enddo

  do K = 0,KC
    do L = 1,LC
      ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      QQ(L,K) = QQMIN
      QQ1(L,K) = QQMIN
      QQ2(L,K) = QQMIN
      QQSQR(L,K) = SQRT(QQMIN)
      if( ISGOTM > 0 )then
        TKE3D(L,K)  = k_min
        TKE3D1(L,K) = k_min
        EPS3D(L,K)  = eps_min
        EPS3D1(L,K) = eps_min
      endif
    enddo
  enddo

  if( MDCHH >= 1 )then
    do NMD = 1,MDCHH
      LHOST = LMDCHH(NMD)
      LCHNU = LMDCHU(NMD)
      LCHNV = LMDCHV(NMD)
  !
  !         SET HOST DRYING DEPTH
  !
      if( PMDCH(NMD) < 0.0) PMDCH(NMD) = HWET
  !
  !         X-DIRECTION CHANNEL
  !
      if( MDCHTYP(NMD) == 1 )then
        if( CHANLEN(NMD) < 0.0 )then
          CHANLEN(NMD) = 0.25*DYP(LHOST)
        else
          CHANLEN(NMD) = CHANLEN(NMD)-0.5*DYP(LCHNU)
        endif
      endif
  !
  !         Y-DIRECTION CHANNEL
  !
      if( MDCHTYP(NMD) == 2 )then
        if( CHANLEN(NMD) < 0.0 )then
          CHANLEN(NMD) = 0.25*DXP(LHOST)
        else
          CHANLEN(NMD) = CHANLEN(NMD)-0.5*DXP(LCHNV)
        endif
      endif
    enddo
  endif

  ! *** INITIALIZE ORGANIC CARBON VARIABLES OF SEDIMENT-TOXICS
  IVAL = 0
  do NT = 1,NTOX
    if( ISTOC(NT) > 0 ) IVAL = 1
  enddo

  if( IVAL == 0 .and. ISTRAN(5) > 0 )then
    ! *** Kd APPROACH ONLY USED IN MODEL
    do NS = 1,NSED+NSND
      do K = 1,KB
        do L = 1,LC
          STFPOCB(L,K,NS) = 1.0
        enddo
      enddo
    enddo
    
    do NS = 1,NSED+NSND
      do K = 1,KC
        do L = 1,LC
          STFPOCW(L,K,NS) = 1.0
        enddo
      enddo
    enddo
  endif
  
  if( IVAL == 1 .and. ISTRAN(5) > 0 )then
    ! *** IF ANY fPOC, POC AND DOC OPTIONS USED INITIALIZE HERE
    
    ! *** SEDIMENT BED
    if( ISTDOCB == 0 )then  !   ISTDOCB == 1 STDOCB IS INITIALIZED FROM DOCB.INP
      do K = 1,KB
        do L = 1,LC
          STDOCB(L,K) = STDOCBC
        enddo
      enddo
    endif
    
    if( ISTPOCB == 0 )then  !   ISTPOCB == 1 STPOCB IS INITIALIZED FROM POCB.INP
      do K = 1,KB
        do L = 1,LC
          STPOCB(L,K) = STPOCBC
        enddo
      enddo
    endif
    
    if( ISTPOCB /= 3 )then
      ! *** SPATIALLY CONSTANT FPOC.  ISTPOCB == 3 STFPOCB IS INITIALIZED FROM FPOCB.INP
      do NS = 1,NSED+NSND
        do K = 1,KB
          do L = 1,LC
            STFPOCB(L,K,NS) = FPOCBST(NS,1)
          enddo
        enddo
      enddo
    endif
    
    ! *** WATER COLUMN
    if( ISTDOCW == 0 )then
      do K = 1,KC
        do L = 1,LC
          STDOCW(L,K) = STDOCWC
        enddo
      enddo
    endif
    
    if( ISTPOCW == 0 )then
      if( STPOCWC <= 0.0 ) STPOCWC = 1E-12
      do K = 1,KC
        
        do L = 1,LC
          STPOCW(L,K) = STPOCWC
        enddo
      enddo
    endif
    
    if( ISTPOCW /= 3 )then
      ! *** SPATIALLY CONSTANT FPOC.  ISTPOCW == 3 STFPOCB IS INITIALIZED FROM FPOCW.INP
      do NS = 1,NSED+NSND
        do K = 1,KC
          do L = 1,LC
            STFPOCW(L,K,NS) = FPOCWST(NS,1)
          enddo
        enddo
      enddo
    endif
  endif

  ! *** NEW VARIABLES FOR QCTL NQCTYP = 3 & 4 
  LOWCHORDU = -9999.
  LOWCHORDV = -9999.
  NLOWCHORD = 0

  return
  
END SUBROUTINE AINIT

