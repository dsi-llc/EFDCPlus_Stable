! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE HIFREQOUT
  !PURPOSE:
  ! WRITE TIME SERIES IN HIGH FREQUENCY
  ! FOR A SET OF SPECIFIC CELLS
  ! BY DANG HUU CHUNG, JULY 2012
  use GLOBAL
  use Variables_MPI
  use Variables_MPI_Write_Out
  use Variables_WQ
  use WQ_RPEM_MODULE
  
  implicit none
  integer(IK4),PRIVATE,parameter :: UHFR = 400
  contains

  SUBROUTINE HFREHYOUT(NCAL,NS)
    integer(4  ),intent(IN) :: NCAL,NS
    integer :: L,LN,K,NP,HSIZE,LG,LE
    real   :: UTMPS,VTMPS
    real(8)   :: EETIME
    character(50) :: FNAME,NSTR*3
    logical(4) :: STAT
    integer(IK4) :: VER

    write(NSTR,'(I3.3)') NS
    FNAME  = OUTDIR//'HFRE_HYD_'//TRIM(NSTR)//'.OUT'
    if( NCAL ==  1 )then
      INQUIRE(UNIT = UHFR,OPENED = STAT)
      if( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      open(UHFR,FILE = FNAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      close(UHFR,STATUS = 'DELETE')
      open(UHFR,FILE = FNAME,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      !------ Header contains Support Information for Output Extraction ------
      VER = 10300
      HSIZE = 16+24*NPNT(NS)             ! Header Size = 2*4+1*4+1*4+NPNT(NSS)*(2*4+2*8)
      write(UHFR) VER,INT(HSIZE,4)     ! File Format Version and Header Size
      write(UHFR) KC
      write(UHFR) INT(NPNT(NS),4)      ! Number of stations
      do NP = 1,NPNT(NS) 
        write(UHFR) INT(HFREGRP(NS)%ICEL(NP),4),INT(HFREGRP(NS)%JCEL(NP),4),HFREGRP(NS)%XCEL(NP),HFREGRP(NS)%YCEL(NP)
      enddo
      !------ Header End ------
      close(UHFR)
      if( .not. allocated(UK) )then
        allocate(UK(KC),VK(KC))
      endif
      return
    else
      open(UHFR,FILE = FNAME,POSITION = 'APPEND',STATUS = 'OLD',FORM = FMT_BINARY)
      close(UHFR)
      open(UHFR,FILE = FNAME,POSITION = 'APPEND',STATUS = 'OLD',FORM = FMT_BINARY)
      EETIME = TIMEDAY
    endif

    WRITE (UHFR) EETIME 
    do NP = 1,NPNT(NS)
      LG = LIJ_Global(HFREGRP(NS)%ICEL(NP),HFREGRP(NS)%JCEL(NP))
      LN = LNC_Global(LG)
      LE = LEC_Global(LG)
      do K = KSZ_Global(LG),KC
        UTMPS = 0.5*(RSSBCE_Global(LG)*U_Global(LE,K) + RSSBCW_Global(LG)*U_Global(LG,K))  
        VTMPS = 0.5*(RSSBCN_Global(LG)*V_Global(LN,K) + RSSBCS_Global(LG)*V_Global(LG,K))  
        UK(K) = CUE_Global(LG)*UTMPS + CVE_Global(LG)*VTMPS  
        VK(K) = CUN_Global(LG)*UTMPS + CVN_Global(LG)*VTMPS  
      enddo
      write(UHFR) REAL(BELV_Global(LG),4),REAL(HP_Global(LG),4),(REAL(UK(K),4),REAL(VK(K),4),REAL(W_Global(LG,K),4),K = KSZ_Global(LG),KC)
    enddo
    FLUSH(UHFR)
    close(UHFR)
  END SUBROUTINE

  SUBROUTINE HFREWCOUT(NCAL,NSS)
    !FOR WATER COLUMN
    integer(IK4),intent(IN) :: NCAL,NSS
    integer :: NP,NACTIVE,NS,HSIZE
    integer :: I,L,K,NT,NX,MD,LG
    real(8) :: EETIME
    character(50) :: FNAME1,NSTR*3
    logical(4) :: STAT
    integer(IK4) :: VER

    if( ALL(ISTRAN(1:7) == 0)) return

    if( ISDYNSTP == 0 )then  
      DELT = DT  
    else  
      DELT = DTDYN  
    endif  
    EETIME = TIMESEC
    EETIME = EETIME/86400._8
    NACTIVE = LA-1
    write(NSTR,'(I3.3)') NSS
    FNAME1  = OUTDIR//'HFRE_WC_'//TRIM(NSTR)//'.OUT'

    if( NCAL == 1 )then
      INQUIRE(UNIT = UHFR,OPENED = STAT)
      if( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      open(UHFR,FILE = FNAME1,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      close(UHFR,STATUS = 'DELETE')
      open(UHFR,FILE = FNAME1,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      !------ Header contains Support Information for Output Extraction ------
      VER = 10300
      HSIZE = 56+24*NPNT(NSS)             ! Header Size = 2*4+7*4+4*4+1*4+NPNT(NSS)*(2*4+2*8)
      write(UHFR) VER,INT(HSIZE,4)      ! File Format Version and Header Size
      write(UHFR) (ISTRAN(I),I = 1,7)
      write(UHFR) KC,NSED,NSND,NTOX
      write(UHFR) INT(NPNT(NSS),4)             ! Number of stations
      do NP = 1,NPNT(NSS) 
        write(UHFR) INT(HFREGRP(NSS)%ICEL(NP),4),INT(HFREGRP(NSS)%JCEL(NP),4),HFREGRP(NSS)%XCEL(NP),HFREGRP(NSS)%YCEL(NP)
      enddo
      !------ Header End ------
      close(UHFR)
      return
    endif

    open(UHFR,FILE = FNAME1,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY)
    close(UHFR)
    open(UHFR,FILE = FNAME1,STATUS = 'OLD',POSITION = 'APPEND',FORM = FMT_BINARY)
    write(UHFR) EETIME,NACTIVE
    if( ISTRAN(6) == 1 .or. ISTRAN(7) >= 1 )then
      do NP = 1,NPNT(NSS) 
        LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        write(UHFR)KBT_Global(LG)
      enddo
    endif
    
    do NP = 1,NPNT(NSS) 
      LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      if( ISTRAN(1) == 1 ) WRITE(UHFR) (REAL(SAL_Global(LG,K),4),K = KSZ_Global(LG),KC)
      if( ISTRAN(2) == 1 ) WRITE(UHFR) (REAL(TEM_Global(LG,K),4),K = KSZ_Global(LG),KC)
      if( ISTRAN(3) == 1 ) WRITE(UHFR) ((REAL(DYE_Global(LG,K,MD),4),K = KSZ_Global(LG),KC),MD = 1,NDYE)
      if( ISTRAN(4) == 1 ) WRITE(UHFR) (REAL(SFL_Global(LG,K),4),K = KSZ_Global(LG),KC)
      if( ISTRAN(5) == 1 ) WRITE(UHFR) ((REAL(TOX_Global(LG,K,NT),4),K = KSZ_Global(LG),KC),NT = 1,NTOX)     
      if( ISTRAN(6) == 1 ) WRITE(UHFR) ((REAL(SED_Global(LG,K,NS),4),K = KSZ_Global(LG),KC),NS = 1,NSED)
      if( ISTRAN(7) == 1 ) WRITE(UHFR) ((REAL(SND_Global(LG,K,NX),4),K = KSZ_Global(LG),KC),NX = 1,NSND)
    enddo
    FLUSH(UHFR)
    close(UHFR)
  END SUBROUTINE

  SUBROUTINE HFREWQOUT(NCAL,NSS)
    !FOR WATER QUALITY
    integer(IK4),intent(IN) :: NCAL,NSS
    integer(IK4) :: VER,NWQVAR,IWQ(40)
    integer :: HSIZE,NP,NACTIVE
    integer :: NW,MW
    integer :: L,K,LG
    
    real(RK4)    :: WQ
    
    real(8) :: EETIME

    character(50) :: FNAME4,NSTR*3
    
    logical(4) :: STAT
    
    save IWQ
    save NWQVAR

    if( ISTRAN(8) == 0 ) return
  
    if( ISDYNSTP == 0 )then  
      DELT = DT  
    else  
      DELT = DTDYN  
    endif  
    EETIME = TIMESEC

    EETIME = EETIME/86400._8
    NACTIVE = LA-1
    write(NSTR,'(I3.3)') NSS
    FNAME4  = OUTDIR//'HFRE_WQ_'//TRIM(NSTR)//'.OUT'

    if( NCAL == 1 )then     
      INQUIRE(UNIT = UHFR,OPENED = STAT)
      if( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      open(UHFR,FILE = FNAME4,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      close(UHFR,STATUS = 'DELETE')
      open(UHFR,FILE = FNAME4,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      NWQVAR = NWQV
 
      !------ Header contains Support Information for Output Extraction ------
      VER = 10300
      HSIZE = 20 + 24*NPNT(NSS) + 4*NWQVAR    ! Header Size = 2*4+2*4+1*4 + NPNT(NSS)*(2*4+2*8) + NWQVAR*4
      write(UHFR) VER,INT(HSIZE,4)          ! File Format Version and Header Size
      write(UHFR) KC,NWQVAR
      write(UHFR) INT(NPNT(NSS),4)          ! Number of stations
      do NP = 1,NPNT(NSS) 
        LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        write(UHFR) INT(HFREGRP(NSS)%ICEL(NP),4),INT(HFREGRP(NSS)%JCEL(NP),4),HFREGRP(NSS)%XCEL(NP),HFREGRP(NSS)%YCEL(NP)
      enddo
      IWQ = 0
      do MW = 1,NWQVAR
        IWQ(MW) = ISKINETICS(MW)
      enddo
      
      write(UHFR) (IWQ(NW),NW = 1,NWQVAR)
      !------ Header End ------
      close(UHFR)  
      return
    endif
  
    open(UHFR,FILE = FNAME4,STATUS = 'UNKNOWN',POSITION = 'APPEND',FORM = FMT_BINARY)
    close(UHFR)
    open(UHFR,FILE = FNAME4,STATUS = 'UNKNOWN',POSITION = 'APPEND',FORM = FMT_BINARY)
    write(UHFR) EETIME
    do NP = 1,NPNT(NSS)
      LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      do K = KSZ_Global(LG),KC
        do NW = 1,NWQVAR
          if( IWQ(NW) > 0 )then
            WQ = WQV_Global(LG,K,NW)
            write(UHFR) WQ
          endif
        enddo
      enddo
    enddo
    FLUSH(UHFR)
    close(UHFR)   
  END SUBROUTINE

  SUBROUTINE HFRERPEMOUT(NCAL,NSS)
    integer(4),intent(IN) :: NCAL,NSS
    integer :: NP, L, LG, HSIZE
    integer :: NACTIVE
    real(8)     :: EETIME
    character(50) :: FNAME6,NSTR*3
    logical(4) :: STAT
    integer(IK4) :: VER

    if( ISTRAN(8) == 0 .or. ISRPEM == 0 ) return
  
    if( ISDYNSTP == 0 )then  
      DELT = DT  
    else  
      DELT = DTDYN  
    endif  
    EETIME = TIMESEC
    EETIME = EETIME/86400._8
    NACTIVE = LA-1

    write(NSTR,'(I3.3)') NSS
    FNAME6  = OUTDIR//'HFRE_RPEM_'//TRIM(NSTR)//'.OUT'

    if( NCAL == 1 )then      
      INQUIRE(UNIT = UHFR,OPENED = STAT)
      if( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      open(UHFR,FILE = FNAME6,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      close(UHFR,STATUS = 'DELETE')
      open(UHFR,FILE = FNAME6,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY)
      !------ Header contains Support Information for Output Extraction ------
      VER = 10300
      HSIZE = 20 + 28*NPNT(NSS)           ! Header Size = 2*4+2*4+1*4 + NPNT(NSS)*(2*4+2*8+4)
      write(UHFR) VER,INT(HSIZE,4)      ! File Format Version and Header Size
      write(UHFR) NRPEM,INITRPEM
      write(UHFR) INT(NPNT(NSS),4)      ! Number of stations
      do NP = 1,NPNT(NSS) 
        LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        write(UHFR) INT(HFREGRP(NSS)%ICEL(NP),4), INT(HFREGRP(NSS)%JCEL(NP),4), HFREGRP(NSS)%XCEL(NP), HFREGRP(NSS)%YCEL(NP), LMASKRPEM_Global(LG)
      enddo
      !------ Header End ------
      close(UHFR)
      return
    endif

    open(UHFR,FILE = FNAME6,STATUS = 'UNKNOWN',POSITION = 'APPEND',FORM = FMT_BINARY)
    close(UHFR)
    open(UHFR,FILE = FNAME6,STATUS = 'UNKNOWN',POSITION = 'APPEND',FORM = FMT_BINARY)
    write(UHFR) EETIME
    do NP = 1,NPNT(NSS)
      LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      if( LMASKRPEM_Global(LG) .or. INITRPEM == 0 )then                                                                                               
        write(UHFR) REAL(WQRPS_Global(LG),4), REAL(WQRPR_Global(LG),4), REAL(WQRPE_Global(LG),4), REAL(WQRPD_Global(LG),4)
      endif
    enddo       
    FLUSH(UHFR)
    close(UHFR)
  END SUBROUTINE
  
  
  
  SUBROUTINE Gather_High_Frequency
  !---------------------------------------------------------------------------!
  ! *** Gathers values for high frequency output
  !---------------------------------------------------------------------------!
  use Mod_Assign_Loc_Glob_For_Write
  use Mod_Map_Write_EE_Binary

  implicit none

  ! *** Local variables
  integer :: j
   
  ! *** Clear out the pointer array
  j = 0

  ! *** Water Elevation Section
  j = j + 1
  call Assign_Loc_Glob_For_Write(j, size(HP,1), HP, size(HP_Global,1), HP_Global)

  j = j + 1
  call Assign_Loc_Glob_For_Write(j, size(BELV,1), BELV, size(BELV_Global,1), BELV_Global)
  
  ! *** Velocity Section  
  j = j + 1
  call Assign_Loc_Glob_For_Write(j, size(U,1), size(U,2), U, size(U_Global,1), size(U_Global,2), U_Global)

  j = j + 1
  call Assign_Loc_Glob_For_Write(j, size(V,1), size(V,2), V, size(V_Global,1), size(V_Global,2), V_Global)

  j = j + 1
  call Assign_Loc_Glob_For_Write(j, size(W,1), size(W,2), W, size(W_Global,1), size(W_Global,2), W_Global)
  
  ! *** Water Column Section
  if( ISTRAN(1) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(SAL,1), size(SAL,2), SAL, &
      size(SAL_Global,1), size(SAL_Global,2), SAL_Global)
  endif

  if( ISTRAN(2) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(TEM,1), size(TEM,2), TEM, &
      size(TEM_Global,1), size(TEM_Global,2), TEM_Global)
  endif

  if( ISTRAN(3) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(DYE,1), size(DYE,2), size(DYE,3), DYE, &
      size(DYE_Global,1), size(DYE_Global,2), size(DYE_Global,3), DYE_Global)
  endif

  if( ISTRAN(4) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(SFL,1), size(SFL,2), SFL, &
      size(SFL_Global,1), size(SFL_Global,2), SFL_Global)
  endif

  if( ISTRAN(5) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(TOX,1), size(TOX,2), size(TOX,3), TOX,&
      size(TOX_Global,1), size(TOX_Global,2), size(TOX_Global,3), TOX_Global)
  endif

  if( ISTRAN(6) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(SED,1), size(SED,2), size(SED,3), SED,&
      size(SED_Global,1), size(SED_Global,2), size(SED_Global,3), SED_Global)
  endif

  if( ISTRAN(7) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(SND,1), size(SND,2), size(SND,3), SND,&
      size(SND_Global,1), size(SND_Global,2), size(SND_Global,3), SND_Global)
  endif
    
  if( ISTRAN(8) >= 1 )then
    j = j + 1
    call Assign_Loc_Glob_For_Write(j, size(WQV,1),        size(WQV,2),        size(WQV,3),        WQV, &
                                      size(WQV_Global,1), size(WQV_Global,2), size(WQV_Global,3), WQV_Global)
  endif

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  call Handle_Calls_MapGatherSort(num_arrays_to_write_out)
  
  END SUBROUTINE Gather_High_Frequency

END MODULE
