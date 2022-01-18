! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE VARZEROReal
  
  ! *** THIS SUBROUTINE ZERO'S ALL OF THE ARRAYS AFTER ALLOCATION

  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 

  USE GLOBAL  
  Use Variables_WQ
  
  USE Variables_MPI
  
  IMPLICIT NONE
  
  INTEGER L,K

  ! *** REAL ARRAYS
  if( process_id == master_id )THEN
    WRITE(*,'(A)')'INITIALIZING REAL ARRAYS'  
  end if
  
  AAU = 0.0
  AAV = 0.0
  AB = 0.0
  ACCWX = 0.0
  AGWELV = 0.0
  AGWELV1 = 0.0
  AGWELV2 = 0.0
  AH = 0.0
  AHDXY = 0.0
  AHOXY = 0.0
  AHC = 0.0
  AHU = 0.0     ! *** NOT USED
  AHV = 0.0     ! *** NOT USED
  AP = 0.0
  APCG = 0.0
  APT = 0.0
  AQ = 0.0
  AQCTL = 0.0
  ASURFEL = 0.0
  ATMP = 0.0
  AV = 0.0
  AVOXY = 0.
  AVBXY = 0.
  AVUI = 0.0
  AVVI = 0.0
  B = 0.0
  B1 = 0.0
  BBU = 0.0
  BBV = 0.0
  BCLSHA = 0.0
  BELAGW = 0.0
  BELSURF = 0.0
  BELV = 0.0
  BELV0 = 0.0
  BELV1 = 0.0
  BLSHA = 0.0
  BSLSHA = 0.0
  BTLSHA = 0.0
  BTMP = 0.0
  CAC = 0.0
  CBE = 0.0
  CBN = 0.0
  CBS = 0.0
  CBW = 0.0
  CC = 0.0
  CCC = 0.0
  CCCC = 0.0
  CCCCHH = 0.0
  CCCCHU = 0.0
  CCCCHV = 0.0
  CCCI = 0.0
  CCCOS = 0.0
  CCCOS1 = 0.0
  CCE = 0.0
  CCI = 0.0
  CCLSHA = 0.0
  CCN = 0.0
  CCNHTT = 0.0
  CCS = 0.0
  CCW = 0.0
  CDZKK = 0.0
  CDZKKP = 0.0
  CDZKMK = 0.0
  CE = 0.0
  CFRD = 0.0
  CHANFRIC = 0.0
  CHANLEN = 0.0
  CLEVAP = 0.0
  CLOE = 0.0
  CLON = 0.0
  CLOS = 0.0
  CLOUDT = 0.0
  CLOW = 0.0
  CLSHA = 0.0
  CMB = 0.0
  CN = 0.0
  CONGW = 0.0
  CONT = 0.0
  CPFAM0 = 0.0
  CPFAM1 = 0.0
  CPFAM2 = 0.0
  CQCJP = 0.0
  CQS = 0.0
  CQSE = 0.0
  CQWR = 0.0
  CQWRSER = 0.0
  CQWRSERT = 0.0
  CS = 0.0
  CSERT = 0.0
  CSLSHA = 0.0
  CTAUC = 0.0
  CTLSHA = 0.0
  CTURBB1 = 0.0
  CTURBB2 = 0.0
  CU1 = 0.0
  CU2 = 0.0
  CUE = 0.0
  CUN = 0.0
  CUU = 0.0
  CVE = 0.0
  CVN = 0.0
  CVV = 0.0
  CW = 0.0
  CWRCJP = 0.0
  DJET = 0.0
  DJPER = 0.0
  DLAT = 0.0
  DLON = 0.0
  DML = 0.0
  DTAUC = 0.0
  DU = 0.0
  DV = 0.0
  DXDJ = 0.0
  DXFP = 0.0
  DXIU = 0.0
  DXIV = 0.0
  DXP = 0.0
  DXU = 0.0
  DXV = 0.0
  DXXTCA = 0.0
  DXYIP = 0.0
  DXYIU = 0.0
  DXYIV = 0.0
  DXYP = 0.0
  DXYTCA = 0.0
  DXYU = 0.0
  DXYV = 0.0
  DYDI = 0.0
  DYE = 0.0
  DYE1 = 0.0
  DYEINIT = 0.0
  DYFP = 0.0
  DYIU = 0.0
  DYIV = 0.0
  DYP = 0.0
  DYU = 0.0
  DYV = 0.0
  DZIC = 0.0
  DZIG = 0.0
  !DZIGSD4 = 0.0
  EVAPGW = 0.0
  EVAPSW = 0.0
  EVAPT = 0.0
  FACBEDL = 0.0
  FACSUSL = 0.0
  FBBX = 0.0
  FBBY = 0.0
  FCAX = 0.0
  FCAXE = 0.0
  FCAY = 0.0
  FCAY1 = 0.0
  FCAY1E = 0.0
  FCAYE = 0.0
  FCORC = 0.0
  FP = 0.0
  FP1 = 0.0
  FPGXE = 0.0
  FPGYE = 0.0
  PMCTESTX = 0.0
  !PMCTESTY = 0.0
  FPROX = 0.0
  FPTMP = 0.0
  FSCORTBCV = 0.0
  FUDISP = 0.0
  FUHDYE = 0.0
  FUHU = 0.0
  FUHV = 0.0
  FVHDXE = 0.0
  FVHU = 0.0
  FVHV = 0.0
  FWQQ = 0.0
  FWQQL = 0.0
  FWU = 0.0
  FWV = 0.0
  FX = 0.0
  FXE = 0.0
  FY = 0.0
  FYE = 0.0
  FZU = 0.0
  FZV = 0.0
  GAMB = 0.0
  GLSHA = 0.0
  GWCSER = 0.0
  GWCSERT = 0.0
  GWFAC = 0.0
  GWSER = 0.0
  GWSERT = 0.0
  H1P = 0.0
  H1U = 0.0
  H1UI = 0.0
  H1V = 0.0
  H1VI = 0.0
  H2P = 0.0
  HBEDA = 0.0
  HBEDA1 = 0.0
  HCTLDA = 0.0
  HCTLDM = 0.0
  HCTLUA = 0.0
  HCTLUM = 0.0
  HDFUFX = 0.0
  HDFUFY = 0.0
  HDFUF = 0.0
  HDIFCTD = 0.0
  HDIFCTL = 0.0
  HGDH = 0.0
  HMP = 0.0
  HMPW = 0.0
  HMU = 0.0
  HMUW = 0.0
  HMV = 0.0
  HMVW = 0.0
  HP = 0.0
  HPI = 0.0
  HPTMP = 0.0
  HRU = 0.0
  HRUO = 0.0
  HRV = 0.0
  HRVO = 0.0
  HTMP = 0.0
  HU = 0.0
  HUDRY = 0.0
  HUI = 0.0
  HUTMP = 0.0
  HUWET = 0.0
  HV = 0.0
  HVDRY = 0.0
  HVI = 0.0
  HVTMP = 0.0
  HVWET = 0.0
  HWQ = 0.0
  HWQI = 0.0

  P = 0.0
  P1 = 0.0
  P1DT1 = 0.0
  PATMT = 0.0
  PCBE = 0.0
  PCBN = 0.0
  PCBS = 0.0
  PCBW = 0.0
  PCLSHA = 0.0
  PDIFTOX = 0.0
  PSERZDF = 0.0
  PSERZDS = 0.0
  PFAM = 0.0
  PFAM1 = 0.0
  PFAM2 = 0.0
  PFPH = 0.0
  PFPH1 = 0.0
  PFPH2 = 0.0
  PHASEE = 0.0
  PHASEU = 0.0
  PHASEV = 0.0
  PHJET = 0.0
  PLSHA = 0.0
  PMDCH = 0.0
  PNHYDS = 0.0
  PPH = 0.0
  PSBE = 0.0
  PSBN = 0.0
  PSBS = 0.0
  PSBW = 0.0
  PSER = 0.0
  PSERAVG = 0.0
  PSERS = 0.0
  PSERT = 0.0
  PSERST = 0.0
  PSHADE = 0.0
  PSLSHA = 0.0
  PTLSHA = 0.0
  PTMP = 0.0
  QCELLCTR = 0.0
  QCHANU = 0.0
  QCHANUN = 0.0
  QCHANV = 0.0
  QCHANVN = 0.0
  QCHNULP = 0.0
  QCHNVLP = 0.0
  QCTL = 0.0
  QCTLST = 0.0
  QCTLSTO = 0.0
  QCTLT = 0.0
  QCTLTLP = 0.0
  QCTLTO = 0.0
  QDWASTE = 0.0
  QFACTOR = 0.0
  QGW = 0.0
  QJPENT = 0.0
  QJPENTT = 0.0
  QMORPH = 0.0
  QQ = 0.0
  QQ1 = 0.0
  QQ2 = 0.0
  QQSQR = 0.0
  QQCJP = 0.0
  QQL = 0.0
  QQL1 = 0.0
  QQL2 = 0.0
  QQWC = 0.0
  QQWCR = 0.0
  QQWV1 = 0.0
  QQWV2 = 0.0
  QQWV3 = 0.0
  QRAIN = 0.0
  QSERT = 0.0
  QSERCELL = 0.0
  QSRTLPN = 0.0
  QSRTLPP = 0.0
  QSS = 0.0
  !QSSE = 0.0 !< Removed and placed into domain decomp because this is now local
  QSSDPA = 0.0
  QSUM = 0.0
  QSUME = 0.0
  QSUM1E = 0.0
  QWATPA = 0.0
  QWBDTOP = 0.0
  QWIDTH = 0.0
  QWRCJP = 0.0
  QWRSER = 0.0
  QWRSERT = 0.0
  QWRSERTLP = 0.0
  RADBOT = 0.0
  RADNET = 0.0
  RADTOP = 0.0
  RADKE  = 0.0
  RAINT = 0.0
  RBPSBL = 0.0
  RCX = 0.0
  RCY = 0.0
  RHAT = 0.0
  RHS = 0.0
  RMAJ = 0.0
  RMIN = 0.0
  ROUSE = 0.0
  RQSMUL = 0.0
  RSOL = 0.0
  RSSBCE = 0.0
  RSSBCN = 0.0
  RSSBCS = 0.0
  RSSBCW = 0.0
  SAAX = 0.0
  SAAY = 0.0
  SAL = 0.0
  SAL1 = 0.0
  SALINIT = 0.0
  SBX = 0.0
  SBXO = 0.0
  SBY = 0.0
  SBYO = 0.0
  SCAX = 0.0
  SCAY = 0.0
  SCB = 0.0
  SCLSHA = 0.0
  SDX = 0.0
  SDY = 0.0
  
  SFL = 0.0
  SFL2 = 0.0
  SFLINIT = 0.0
  SFLSBOT = 0.0
  
  SIGPHIA = 0.0
  SNAPSHOTS = 0.0
  
  SOLSWRT = 0.0
  SPB = 0.0
  SPBE1 = 0.0
  SPBE2 = 0.0
  SPBN1 = 0.0
  SPBN2 = 0.0
  SPBS1 = 0.0
  SPBS2 = 0.0
  SPBW1 = 0.0
  SPBW2 = 0.0
  SPFAM0 = 0.0
  SPFAM1 = 0.0
  SPFAM2 = 0.0
  SSLSHA = 0.0
  SSSIN = 0.0
  SSSIN1 = 0.0
  SSSS = 0.0
  STBX = 0.0
  STBXO = 0.0
  STBY = 0.0
  STBYO = 0.0
  STCAP = 0.0
  STCUV = 0.0
  STLSHA = 0.0
  SUB = 0.0
  SUBO = 0.0
  SVB = 0.0
  SVBO = 0.0
  SVPAT = 0.0
  SVPW = 0.0
  SWB = 0.0
  TACSER = 0.0
  TAGWSER = 0.0
  TAPSER = 0.0
  TAQSER = 0.0
  TAQWRSR = 0.0
  TATMT = 0.0
  TAUB = 0.0
  TAVEGSER = 0.0
  TBX = 0.0
  TBX1 = 0.0
  TBY = 0.0
  TBY1 = 0.0
  TCCSER = 0.0
  TCGWSER = 0.0
  TCP = 0.0
  TCPSER = 0.0
  TCQSER = 0.0
  TCQWRSR = 0.0
  TCVEGSER = 0.0
  TDEWT = 0.0
  TEM = 0.0
  TEM1 = 0.0
  TEMB = 0.0
  TEMB1 = 0.0
  TEMINIT = 0.0
  TGWSER = 0.0
  THJET = 0.0
  TMP3D = 0.0
  
  TPCOORDE = 0.0
  TPCOORDN = 0.0
  TPCOORDS = 0.0
  TPCOORDW = 0.0
  TPSER = 0.0
  TQWRSER = 0.0
  TSSTOP = 0.0
  TSSTRT = 0.0
  TSX = 0.0
  TSX1 = 0.0
  TSY = 0.0
  TSY1 = 0.0
  TVAR1E = 0.0
  TVAR1N = 0.0
  TVAR1S = 0.0
  TVAR1W = 0.0
  TVAR2C = 0.0
  TVAR2E = 0.0
  TVAR2N = 0.0
  TVAR2S = 0.0
  TVAR2W = 0.0
  TVAR3C = 0.0
  TVAR3E = 0.0
  TVAR3N = 0.0
  TVAR3S = 0.0
  TVAR3W = 0.0
  TVEGSER = 0.0
  TWATER = 0.0
  U = 0.0
  U0 = 0.0
  U1 = 0.0
  U1V = 0.0
  U2 = 0.0
  UCELLCTR = 0.0
  UCLSHA = 0.0
  UCOS = 0.0
  UE0 = 0.0
  UE1DT1 = 0.0
  UECLSHA = 0.0
  UELSHA = 0.0
  UESLSHA = 0.0
  UETLSHA = 0.0
  UHDY1E = 0.0
  UHDY2E = 0.0
  UHDYE = 0.0
  UHDYEK = 0.
  UHDY1EK = 0.
  UHDY2EK = 0.
  UHE = 0.0
  USIN = 0.0
  USLSHA = 0.0
  USTAR = 0.0
  USTARSED = 0.0
  USTARSND = 0.0
  UUU = 0.0
  UV = 0.0
  UWVSQ = 0.0
  V = 0.0
  V0 = 0.0
  V1 = 0.0
  V1U = 0.0
  V2 = 0.0
  VCELLCTR = 0.0
  VCLSHA = 0.0
  VCOS = 0.0
  VDWASTE = 0.0
  VE0 = 0.0
  VE1DT1 = 0.0
  VECLSHA = 0.0
  VEGK = 0.0
  VEGSERB = 0.0
  VEGSERBT = 0.0
  VEGSERH = 0.0
  VEGSERHT = 0.0
  VEGSERR = 0.0
  VEGSERRT = 0.0
  VELSHA = 0.0
  VESLSHA = 0.0
  VETLSHA = 0.0
  VHDX1E = 0.0
  VHDX2E = 0.0
  VHDXE = 0.0
  VHDXEK = 0.
  VHDX1EK = 0.
  VHDX2EK = 0.
  VHE = 0.0
  VOLPERC = 0.0
  VOLSEL = 0.0
  VPAT = 0.0
  VSIN = 0.0
  VSLSHA = 0.0
  VU = 0.0
  VVLSHA = 0.0
  VVV = 0.0
  W = 0.0
  W1 = 0.0
  W2 = 0.0
  WACCWE = 0.0
  WC = 0.0
  WC2 = 0.0
  WCOREST = 0.0
  WCORWST = 0.0
  WCORNTH = 0.0
  WCORSTH = 0.0
  WINDCD10 = 0.0
  WINDST = 0.0
  WINDSTKA = 0.0
  WINDSTKA_SAVE = 0.0
  WKQ = 0.0
  WLSHA = 0.0
  WNDVELE = 0.0
  WNDVELN = 0.0
  
  WTCI = 0.0
  WVDTKEM = 0.0
  WVDTKEP = 0.0
  WVENEP = 0.0
  WVHUU = 0.0
  WVHUV = 0.0
  WVHVV = 0.0
  WVKHC = 0.0
  WVKHU = 0.0
  WVKHV = 0.0
  WVPP = 0.0
  WVPT = 0.0
  WVPU = 0.0
  WVPV = 0.0
  WVTMP1 = 0.0
  WVTMP2 = 0.0
  WVTMP3 = 0.0
  WVTMP4 = 0.0
  WWW = 0.0
  XJETL = 0.0
  YJETL = 0.0
  Z = 0.0
  ZBR = 0.0
  ZBRE = 0.0

  ZEQ = 0.0
  ZEQD = 0.0
  ZEQDI = 0.0
  ZEQI = 0.0
  ZETATOP = 0.0
  ZJET = 0.0
  ZP = 0.0
  ZZ = 0.0
  ZZC = 0.0
  ZZP = 0.0

  IF( ISINWV == 1 )THEN
    CFLCAC = 0.0
    CFLUUU = 0.0
    CFLVVV = 0.0
    CFLWWW = 0.0
  ENDIF
  IF( ISDISP > 0 )THEN
    BDISP = 0.0
    CUDISPT = 0.0
    CVDISPT = 0.0
    DYXTCA = 0.0
    DYYTCA = 0.0
    FVDISP = 0.0
  ENDIF
  IF( ISBODYF > 0 )THEN
    FBODYFX = 0.0
    FBODYFY = 0.0
  ENDIF
  IF( ISHDMF >= 1 )THEN
    FMDUX = 0.0
    FMDUY = 0.0
    FMDVX = 0.0
    FMDVY = 0.0
    IF( IS2TIM == 0 )H1C = 0.0
  ENDIF    
  IF( ISVEG > 0 )THEN
    FXVEG = 0.0
    FXVEGE = 0.0
    FYVEG = 0.0
    FYVEGE = 0.0
    ALPVEG = 0.0
    BDLPSQ = 0.0
    !BETVEG = 0.0
    BPVEG = 0.0
    !GAMVEG = 0.0
    HPVEG = 0.0
    !PVEGX = 0.0
    !PVEGY = 0.0
    PVEGZ = 0.0
    RDLPSQ = 0.0
    SCVEG = 0.0
  ENDIF
  IF( ISWAVE > 0 )THEN
    FXWAVE = 0.0
    FYWAVE = 0.0
  ENDIF
  IF(  ISWAVE > 0  )THEN
    DO L=1,LCM
      WV(L).TWX = 0
      WV(L).TWY = 0
      WV(L).K = 0
      WV(L).UDEL = 0
      WV(L).LENGTH = 0
      WV(L).HEIGHT = 0
      WV(L).HEISIG = 0
      WV(L).DIR = 0
      WV(L).FREQ = 0
      WV(L).PERIOD = 0
      DO K=1,KCM
        WV(L).DISSIPA(K) = 0
      ENDDO
    ENDDO
  ENDIF
  
  IF( ISVHEAT > 0 )THEN
    LSVHTWINDE = .FALSE.
    LSVHTWINDC = .FALSE.
    SVREVC = 0.
    SVRCHC = 0.
  ENDIF
  SVKEBACK = 0.
  
  ! *** FOODCHAIN MODELLING OPTIONS
  IF( ISTPOCB == 4 )THEN
    PFPOCB = 0.0
    FPOCB = 0.0
  ELSE
    PFPOCB = 0.0
    FPOCB = 0.0
  ENDIF

  ! ** TOXIC TRANSPORT VARIABLES
  CONPARB = 0.0
  CONPARBC = 0.0
  CONPARW = 0.0
  CONPARWC = 0.0
  DIFTOX = 0.0
  DIFTOXS = 0.0
  DPDIFTOX = 0.0
  FPOCBST = 0.0
  FPOCWST = 0.0
  PDIFTOX = 0.0
  RKTOXP = 0.0
  SKTOXP = 0.0
  TOXPARB = 0.0
  TOXPARBC = 0.0
  TOXPARW = 0.0
  TOXPARWC = 0.0
  TOX = 0.0
  TOX_BIO_KB = 0.0
  TOX_BIO_KW = 0.0
  TOX_BIO_MXD = 0.0
  TOX_BIO_Q10B = 0.0
  TOX_BIO_Q10W = 0.0
  TOX_BIO_TB = 0.0
  TOX_BIO_TW = 0.0
  TOX_BLK_KB = 0.0
  TOX_BLK_KW = 0.0
  TOX_BLK_MXD = 0.0
  TOX_ADJ = 0.0
  TOX_ATM = 0.0
  TOX_HE = 0.0
  TOX_MW = 0.0
  TOX_KV_TCOEFF = 0.0
  TOX1 = 0.0

  IF( ISTRAN(5) >= 1 )THEN
    ISDIFBW = 0.0
    ISTOC = 0.0
    ISTOXR = 0.0
    ITXBDUT = 0.0
    ITXINT = 0.0
    ITXPARB = 0.0
    ITXPARBC = 0.0
    ITXPARW = 0.0
    ITXPARWC = 0.0
    NSP2 = 0.0
    RTOXERE2T = 0.0
    RTOXERE2TB = 0.0
    RTOXERE2TW = 0.0
    RTOXERO2T = 0.0
    RTOXERO2TB = 0.0
    RTOXERO2TW = 0.0
    SKTOXP = 0.0
    STPOCW = 0.0
    TADFLUX = 0.0
    TADFLUX2T = 0.0
    TOXATM = 0.0
    TOXB = 0.0
    TOXB1 = 0.0
    TOXBEG2T = 0.0
    TOXBEG2TB = 0.0
    TOXBEG2TW = 0.0
    TOXBINIT = 0.0
    TOXBLB = 0.0
    TOXBLB2T = 0.0
    TOXBMO2T = 0.0
    TOXBMO2TB = 0.0
    TOXBMO2TW = 0.0
    TOXBS = 0.0
    TOXCDFB = 0.0
    TOXCDFW = 0.0
    TOXEND2T = 0.0
    TOXEND2TB = 0.0
    TOXEND2TW = 0.0
    TOXERE2T = 0.0
    TOXERE2TB = 0.0
    TOXERE2TW = 0.0
    TOXERO2T = 0.0
    TOXERO2TB = 0.0
    TOXERO2TW = 0.0
    TOXERR2T = 0.0
    TOXERR2TB = 0.0
    TOXERR2TW = 0.0
    TOXF = 0.0
    TOXFB = 0.0
    TOXFBL = 0.0
    TOXFBL2T = 0.0
    TOXFBLT = 0.0
    TOXFBEBKB = 0.0
    TOXFBECHB = 0.0
    TOXFBECHW = 0.0
    TOXFDFB = 0.0
    TOXFDFW = 0.0
    TOXFLUXB2T = 0.0
    TOXFLUXW2T = 0.0
    TOXINIT = 0.0
    TOXINTB = 0.0
    TOXINTW = 0.0
    TOXOUT2T = 0.0
    TOXOUT2TB = 0.0
    TOXOUT2TW = 0.0
    TOXPFB = 0.0
    TOXPFTB = 0.0
    TOXPFTW = 0.0
    TOXPFW = 0.0
    TOXS = 0.0
    TOXTMP = 0.0
    TOXWBALN = 0.0
    TOXWBALO = 0.0
  ENDIF   ! *** END OF TOXIC VARIABLE DECLARATIONS
  
  ! ** SED VARIABLES
  SED = 0.0
  SED1 = 0.0
  IF( ISTRAN(6) >= 1 )THEN
    SED3DMAX = 0.0
    SED3DMIN = 0.0
    SEDB = 0.0
    SEDB1 = 0.0
    SEDBA = 0.0
    SEDBEG2T = 0.0
    SEDBEG2TB = 0.0
    SEDBEG2TW = 0.0
    SEDBINIT = 0.0
    SEDBMO2T = 0.0
    SEDBMO2TB = 0.0
    SEDBMO2TW = 0.0
    SEDEND2T = 0.0
    SEDEND2TB = 0.0
    SEDEND2TW = 0.0
    SEDERE2T = 0.0
    SEDERE2TB = 0.0
    SEDERE2TW = 0.0
    SEDERO2T = 0.0
    SEDERO2TB = 0.0
    SEDERO2TW = 0.0
    SEDERR2T = 0.0
    SEDERR2TB = 0.0
    SEDERR2TW = 0.0
    SEDF = 0.0
    SEDFLUX2T = 0.0
    SEDINIT = 0.0
    SEDOUT2T = 0.0
    SEDOUT2TB = 0.0
    SEDOUT2TW = 0.0
    SEDS = 0.0
  ENDIF   ! *** END OF SED VARIABLE DECLARATIONS
  
  ! ** SND VARIABLES
  SND = 0.0
  SND1 = 0.0
  IF( NSND > 0 )THEN
    SBDLDA = 0.0
    SBDLDB = 0.0
    SBDLDG1 = 0.0
    SBDLDG2 = 0.0
    SBDLDG3 = 0.0
    SBDLDG4 = 0.0
    SBDLDP = 0.0
  ENDIF
  IF( ISTRAN(7) >= 1 )THEN
    SBLOUT2T = 0.0
    SBLOUT2TB = 0.0
    SBLOUT2TW = 0.0
    SND3DMAX = 0.0
    SND3DMIN = 0.0
    SNDBA = 0.0
    SNDBA1 = 0.0
    SNDBAT = 0.0
    SNDBEG2T = 0.0
    SNDBEG2TB = 0.0
    SNDBEG2TW = 0.0
    SNDBINIT = 0.0
    SNDBMO2T = 0.0
    SNDBMO2TB = 0.0
    SNDBMO2TW = 0.0
    SNDBS = 0.0
    SNDEND2T = 0.0
    SNDEND2TB = 0.0
    SNDEND2TW = 0.0
    SNDEQ = 0.0
    SNDEQSAV = 0.0
    SNDEQB = 0.0
    SNDERE2T = 0.0
    SNDERE2TB = 0.0
    SNDERE2TW = 0.0
    SNDERO2T = 0.0
    SNDERO2TB = 0.0
    SNDERO2TW = 0.0
    SNDERR2T = 0.0
    SNDERR2TB = 0.0
    SNDERR2TW = 0.0
    SNDF = 0.0
    SNDFBL = 0.0
    SNDFBL2T = 0.0
    SNDFLUX2T = 0.0
    SNDINIT = 0.0
    SNDOUT2T = 0.0
    SNDOUT2TB = 0.0
    SNDOUT2TW = 0.0
    SNDS = 0.0
  ENDIF   ! *** END OF SND VARIABLE DECLARATIONS
  
  ! *** EITHER SED OR SND
  COSEDHID = 0.0
  RSNDM = 0.0
  SDBLV = 0.0
  SDEN = 0.0
  SEDN = 0.0
  SEDO = 0.0
  SEDBO = 0.0
  SEDDIA = 0.0
  SEXP = 0.0
  SSG = 0.0
  TAUD = 0.0
  TAUN = 0.0
  TAUR = 0.0
  TCSHIELDS = 0.0
  TEXP = 0.0
  VDRDEPO = 0.0
  VDRRSPO = 0.0
  WRSPO = 0.0
  WSEDO = 0.0
  WS = 0.0
  WS2 = 0.0

  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 .OR. ISBAL > 0 )THEN
    ACOEF = 0.0
    ALOW = 0.0
    BDENBED = 0.0
    BDENBED1 = 0.0
    BDENBEDA = 0.0
    BDENBEDA1 = 0.0
    BEDBINIT = 0.0
    BEDBKDSV = 0.0
    BEDDINIT = 0.0
    BEDLINIT = 0.0
    BEDPORSV = 0.0
    BEDTHKSV = 0.0
    BEDVDRSV = 0.0
    BMNN = 0.0
    CBEDTOTAL = 0.0
    COEFK = 0.0
    COEFSK = 0.0
    CQBEDLOADX = 0.0
    CQBEDLOADY = 0.0
    CSHIELDS50 = 0.0
    CUPP = 0.0
    DSTRSE = 0.0
    DZBTR = 0.0
    DZBTR1 = 0.0
    FRACCOH = 0.0
    FRACNON = 0.0
    GAMTMP = 0.0
    HBED = 0.0
    HBED1 = 0.0
    HYDCN = 0.0

    PEXP = 0.0
    PHID = 0.0
    QSBDLDIN = 0.0
    QSBDLDOT = 0.0
    QSBDLDX = 0.0
    QSBDLDY = 0.0
    QSBDLDP = 0.0
    QSBDTOP = 0.0
    QSEDBED = 0.0
    QSEDBED1 = 0.0
    QSEDBEDA = 0.0
    QSEDBEDA1 = 0.0

    SEDDIA50 = 0.0
    SEDDIA90 = 0.0

    SEDBT = 0.0
    SEDFPA = 0.0
    SEDT = 0.0
    
    SNDB = 0.0
    SNDB1 = 0.0
    SNDBT = 0.0
    SNDFPA = 0.0
    SNDT = 0.0

    TAUBSED = 0.0
    TAUBSND = 0.0
    TAUCRCOH = 0.0
    TAURB = 0.0
    TAUDS = 0.0
    TAUNS = 0.0
    TAURS = 0.0
    TEXPS = 0.0

    VDRBED = 0.0
    VDRBED1 = 0.0
    VDRBED2 = 0.0
    VDRBEDA = 0.0
    VDRBEDSED = 0.0
    VDRBEDSND = 0.0
    VDRHBEDA1 = 0.0
    VFRBED = 0.0
    VFRBED1 = 0.0
    VOLBW2 = 0.0
    VOLBW3 = 0.0
    
    WSETA = 0.0
    WRSPB = 0.0
    WRSPBA = 0.0
    WRSPS = 0.0
    WRSPSA = 0.0

    ZBEDC = 0.0
    ZBEDG = 0.0
    ZBEDGT = 0.0
    IF( ISBEDSTR == 3 )THEN
      ZBRSED = 0.0
    ELSE
      ZBRSED = 0.0
    ENDIF
    ZELBED = 0.0
    ZELBED1 = 0.0
    ZELBEDA = 0.0
    ZELBEDA1 = 0.0

    PARTMIXZ = 0.0
    PMXDEPTH = 0.0
    PMXCOEF = 0.0
    PORBED = 0.0
    PORBED1 = 0.0
    QCOEF = 0.0
    QWTRBED = 0.0
    RRHS = 0.0
    RSEDERE2T = 0.0
    RSEDERE2TB = 0.0
    RSEDERE2TW = 0.0
    RSEDERO2T = 0.0
    RSEDERO2TB = 0.0
    RSEDERO2TW = 0.0
    RSNDERE2T = 0.0
    RSNDERE2TB = 0.0
    RSNDERE2TW = 0.0
    RSNDERO2T = 0.0
    RSNDERO2TB = 0.0
    RSNDERO2TW = 0.0
    SDENAVG = 0.0
    SGSM1 = 0.0
    SIGPHI = 0.0
    STDOCB = 0.0
    STDOCW = 0.0
    STFPOCB = 0.0
    STFPOCW = 0.0
    STPOCB = 0.0
    STRSE = 0.0

    ! *** BANK EROSION
    IBANKBE = 0.0
    JBANKBE = 0.0
    ICHANBE = 0.0
    JCHANBE = 0.0
    NBESERN = 0.0
    MBESER = 0.0
    MBETLAST = 0.0
    LPMXZ = 0.0
  ENDIF
  
  ! *** CALMMT/RESIDUAL/MODEL LINKAGES VARIABLES
  IF( ISSSMMT > 0 .OR. ISWASP > 0 )THEN
    EVPGLPF = 0.0
    EVPSLPF = 0.0
    GWLPF = 0.0
    HLPF = 0.0
    QSUMELPF = 0.0
    QSUMLPF = 0.0
    RAINLPF = 0.0
    RINFLPF = 0.0
    SALLPF = 0.0
    SFLLPF = 0.0
    TEMLPF = 0.0
    UELPF = 0.0
    VELPF = 0.0
    WLPF = 0.0
    WTLPF = 0.0
    
    IF( ISTRAN(5) > 0 )THEN
      TOXBLPF = 0.0
      TOXLPF = 0.0
      TXPFLPF = 0.0
    ENDIF
    IF( ISTRAN(6) > 0 )THEN
      SEDBLPF = 0.0
      SEDBTLPF = 0.0
      SEDLPF = 0.0
      SEDTLPF = 0.0
    ENDIF
    IF( ISTRAN(7) > 0 )THEN
      SNDBLPF = 0.0
      SNDBTLPF = 0.0
      SNDLPF = 0.0
      SNDTLPF = 0.0
    ENDIF
  ENDIF
  
  
  ! ** WATER QUALITY VARIABLES
  WQKEB = 0.0
  
  IF( ISWQFLUX == 1 )THEN
    ABEFF = 0.0
    ABLPF = 0.0
    ACCWFLD = 0.0
    AHULPF = 0.0
    AHVLPF = 0.0
    DYELPF = 0.0

    UHLPF = 0.0
    UIRT = 0.0
    ULPF = 0.0
    ULSHA = 0.0
    UTLPF = 0.0
    UTLSHA = 0.0
    UVPT = 0.0
    U1DT1 = 0.0

    VHLPF = 0.0
    VIRT = 0.0
    VLPF = 0.0
    VLSHA = 0.0
    VTLPF = 0.0
    VTLSHA = 0.0
    VVPT = 0.0
    V1DT1 = 0.0

    WIRT = 0.0
  ENDIF
  
  IF( ISWASP > 0 )THEN
    LCEFDC = 0.0
    LCHNC = 0.0
    NCHNC = 0.0
  ENDIF
  
  ! Begin SEDZLJ variables not zeroed in varzerosnl
  FWVTP = 0.0
  FWDIR = 0.0
  ! End SEDZLJ variables

  ! *** OMP DECLARATIONS
  CD = 0.0
  FQC = 0.0
  FUHVD = 0.0
  FVHVD = 0.0
  FWVV = 0.0
  FQCPAD = 0.0
  QSUMNAD = 0.0
  QSUMPAD = 0.0
  DUU = 0.0
  DVV = 0.0
  POS = 0.0
  UUUU = 0.0
  VVVV = 0.0
  WWWW = 0.0
  CONQ = 0.0
  CMAX = 0.0
  CMIN = 0.0
  CONTD = 0.0
  CONTMN = 0.0
  CONTMX = 0.0
    
  ! Begin MHK variables SCJ
  IJLTURB = 0.0
  CTMHK = 0.0
  CDSUP = 0.0
  DENMHK = 0.0
  ESUP = 0.0
  EMHK = 0.0
  FXMHK = 0.0
  FXMHKE = 0.0
  FXSUP = 0.0
  FXSUPE = 0.0
  FYMHK = 0.0
  FYMHKE = 0.0
  FYSUP = 0.0
  FYSUPE = 0.0
  PMHK = 0.0
  PSUP = 0.0
  VMAXCUT = 0.0
  VMINCUT = 0.0

  ! *** ALLOCATE MEMORY FOR VARIABLE TO STORE CONCENTRATIONS AT OPEN BOUNDARIES
  
  FBESER = 0.0
  BESERT = 0.0
  FWCBESERT = 0.0
  TCBESER = 0.0
  TABESER = 0.0
  TBESER = 0.0
  BESER = 0.0
  FWCBESER = 0.0
  QSBDTOPBEBKB = 0.0
  QSBDTOPBECHB = 0.0
  QSBDTOPBECHW = 0.0
  QWBDTOPBEBKB = 0.0
  QWBDTOPBECHB = 0.0
  QWBDTOPBECHW = 0.0
  SEDFBEBKB = 0.0
  SEDFBECHB = 0.0
  SEDFBECHW = 0.0
  SNDFBEBKB = 0.0
  SNDFBECHB = 0.0
  SNDFBECHW = 0.0

  ! *** NEW VARIABLES FOR QCTL NQCTYP=3 & 4 
  LOWCHORDU = 0.0
  LOWCHORDV = 0.0
  NLOWCHORD = 0.0
  SAVESUB = 0.0
  SAVESVB = 0.0
  
  ! *** SIGMA-Z SGZ
  BI1W = 0.0
  BI2W = 0.0
  BEW = 0.0
  BI1S = 0.0
  BI2S = 0.0
  BES = 0.0
  SGZU = 0.0
  SGZV = 0.0
  IF( IGRIDV > 0 .OR. NBLOCKED > 0 )THEN
    BELVE = 0.0
    BELVW = 0.0
    BELVN = 0.0
    BELVS = 0.0
    BE = 0.0
    BW = 0.0
    BN = 0.0
    BS = 0.0
    BEE = 0.0
    BEN = 0.0
    BI1E = 0.0
    BI1N = 0.0
    BI2E = 0.0
    BI2N = 0.0
    HPE = 0.0
    HPW = 0.0
    HPN = 0.0
    HPS = 0.0
    SGZKE = 0.0
    SGZKW = 0.0
    SGZKN = 0.0
    SGZKS = 0.0
    SGZE = 0.0
    SGZW = 0.0
    SGZN = 0.0
    SGZS = 0.0
    ZE = 0.0
    ZW = 0.0
    ZN = 0.0
    ZS = 0.0
    ZZE = 0.0
    ZZW = 0.0
    ZZN = 0.0
    ZZS = 0.0
  ENDIF
  DZIGSD4U = 0.0
  DZIGSD4V = 0.0

  DZC = 0.0
  DZCK = 0.0
  DZG = 0.0
  DZGU = 0.0
  DZGV = 0.0

  CDZDU = 0.0
  CDZDV = 0.0
  CDZFU = 0.0
  CDZFV = 0.0
  CDZLU = 0.0
  CDZLV = 0.0
  CDZMU = 0.0
  CDZMV = 0.0
  CDZRU = 0.0
  CDZRV = 0.0
  CDZUU = 0.0
  CDZUV = 0.0

  SUB3D = 0.0
  SUB3DO = 0.0
  SVB3D = 0.0
  SVB3DO = 0.0
  
  FSGZU = 0.0
  FSGZV = 0.0

  HPK  = 0.0
  HPKI = 0.0
  H1PK = 0.0
  H2PK = 0.0
  H2WQ = 0.0
  UHDY = 0.0
  VHDX = 0.0
  UHDY1 = 0.0
  VHDX1 = 0.0
  UHDY2 = 0.0
  VHDX2 = 0.0
  UHDYF = 0.0
  VHDXF = 0.0
  UHDYF1 = 0.0
  VHDXF1 = 0.0
  UHDYF2 = 0.0
  VHDXF2 = 0.0

  ! *** HYDRAULIC STRUCTURE EQUATIONS
  IF( NHYDST > 0 )THEN
    HS_LENGTH=0.
    HS_XSTYPE=0
    HS_WIDTH=0.
    HS_HEIGHT=0.
    HS_MANN=0.
    HS_ANGLE=0.
  ENDIF    
  

  ! *** WET/DRY BYPASS VARIABLES
  LWET = 0.0
  LDRY = 0.0
  
  ! *** ATMOS
  IF( NASER > 0 )THEN
    TCASER = 0.0
    TAASER=0.0
    IF( NASER > 1  )THEN
        ATMWHT=0.0
        ATMWHT_TEMP=0.0
    ENDIF
    
  ENDIF
  
  ! *** WIND
  IF( NWSER > 0 )THEN
    TCWSER = 0.0
    ISWDINT = 0.0
    TAWSER=0.0
    IF( NWSER > 1 ) WNDWHT=0.0
  ENDIF
    
  ! *** ICE
  IF( ISICE > 0 )THEN
    MITLAST = 0.0
    FRAZILICE = 0.0
    FRAZILICE1 = 0.0
    ICECOVER = 0.0
    ICERATE = 0.0
    ICETHICK = 0.0
    ICETHICK1 = 0.0
    ICETEMP = 0.0
    ICEVOL = 0.0
    TCISER = 0.0
    TAISER = 0.0
    RICECOVT = 0.0
    RICETHKT = 0.0
    IF( ISICE == 1 .AND. NISER > 1 ) RICEWHT = 0.0
  ENDIF
  ICECELL = 0.0
  RHOW = 0.0

  IF( ISLSHA > 0 .OR. ISHTA >0 )THEN
    AMCP=0.0
    AMCU=0.0
    AMCUE=0.0
    AMCV=0.0
    AMCVE=0.0
    AMPU=0.0
    AMPV=0.0
    AMSP=0.0
    AMSU=0.0
    AMSUE=0.0
    AMSV=0.0
    AMSVE=0.0
  ENDIF
  
  IF( NBLOCKED > 0 )THEN
    BLANCHORU = 0.0
    BLANCHORV = 0.0
    BLDRAFTUO = 0.0
    BLDRAFTVO = 0.0
    BLSILLU = 0.0
    BLSILLV = 0.0
    BLDRAFTU = 0.0
    BLDRAFTV = 0.0
  ENDIF

END SUBROUTINE VARZEROReal
