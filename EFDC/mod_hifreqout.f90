! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE HIFREQOUT
  !PURPOSE:
  ! WRITE TIME SERIES IN HIGH FREQUENCY
  ! FOR A SET OF SPECIFIC CELLS
  ! BY DANG HUU CHUNG, JULY 2012
  USE GLOBAL
  Use Variables_MPI
  Use Variables_MPI_Write_Out
  Use Variables_WQ
  Use WQ_RPEM_MODULE
  
  IMPLICIT NONE
  INTEGER(IK4),PRIVATE,PARAMETER :: UHFR=400
  CONTAINS

  SUBROUTINE HFREHYOUT(NCAL,NS)
    INTEGER(4  ),INTENT(IN) :: NCAL,NS
    INTEGER :: L,LN,K,NP,HSIZE,LG,LE
    REAL   :: UTMPS,VTMPS
    REAL(8)   :: EETIME
    CHARACTER(50) :: FNAME,NSTR*3
    LOGICAL(4) :: STAT
    INTEGER(IK4) :: VER

    WRITE(NSTR,'(I3.3)') NS
    FNAME =OUTDIR//'HFRE_HYD_'//TRIM(NSTR)//'.OUT'
    IF( NCAL ==  1 )THEN
      INQUIRE(UNIT=UHFR,OPENED=STAT)
      IF( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      OPEN(UHFR,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      CLOSE(UHFR,STATUS='DELETE')
      OPEN(UHFR,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      !------ Header contains Support Information for Output Extraction ------
      VER=10300
      HSIZE=16+24*NPNT(NS)             ! Header Size = 2*4+1*4+1*4+NPNT(NSS)*(2*4+2*8)
      WRITE(UHFR) VER,INT(HSIZE,4)     ! File Format Version and Header Size
      WRITE(UHFR) KC
      WRITE(UHFR) INT(NPNT(NS),4)      ! Number of stations
      DO NP=1,NPNT(NS) 
        WRITE(UHFR) INT(HFREGRP(NS)%ICEL(NP),4),INT(HFREGRP(NS)%JCEL(NP),4),HFREGRP(NS)%XCEL(NP),HFREGRP(NS)%YCEL(NP)
      ENDDO
      !------ Header End ------
      CLOSE(UHFR)
      IF(  .NOT. ALLOCATED(UK) )THEN
        ALLOCATE(UK(KC),VK(KC))
      ENDIF
      RETURN
    ELSE
      OPEN(UHFR,FILE=FNAME,POSITION='APPEND',STATUS='OLD',FORM='BINARY')
      CLOSE(UHFR)
      OPEN(UHFR,FILE=FNAME,POSITION='APPEND',STATUS='OLD',FORM='BINARY')
      EETIME = TIMEDAY
    ENDIF

    WRITE (UHFR) EETIME 
    DO NP=1,NPNT(NS)
      LG = LIJ_Global(HFREGRP(NS)%ICEL(NP),HFREGRP(NS)%JCEL(NP))
      LN = LNC_Global(LG)
      LE = LEC_Global(LG)
      DO K = KSZ_Global(LG),KC
        UTMPS = 0.5*(RSSBCE_Global(LG)*U_Global(LE,K) + RSSBCW_Global(LG)*U_Global(LG,K))  
        VTMPS = 0.5*(RSSBCN_Global(LG)*V_Global(LN,K) + RSSBCS_Global(LG)*V_Global(LG,K))  
        UK(K) = CUE_Global(LG)*UTMPS + CVE_Global(LG)*VTMPS  
        VK(K) = CUN_Global(LG)*UTMPS + CVN_Global(LG)*VTMPS  
      ENDDO
      WRITE(UHFR) REAL(BELV_Global(LG),4),REAL(HP_Global(LG),4),(REAL(UK(K),4),REAL(VK(K),4),REAL(W_Global(LG,K),4),K=KSZ_Global(LG),KC)
    ENDDO
    FLUSH(UHFR)
    CLOSE(UHFR)
  END SUBROUTINE

  SUBROUTINE HFREWCOUT(NCAL,NSS)
    !FOR WATER COLUMN
    INTEGER(IK4),INTENT(IN) :: NCAL,NSS
    INTEGER :: NP,NACTIVE,NS,HSIZE
    INTEGER :: I,L,K,NT,NX,MD,LG
    REAL(8) :: EETIME
    CHARACTER(50) :: FNAME1,NSTR*3
    LOGICAL(4) :: STAT
    INTEGER(IK4) :: VER

    IF( ALL(ISTRAN(1:7) == 0)) RETURN

    IF( ISDYNSTP == 0 )THEN  
      DELT=DT  
    ELSE  
      DELT=DTDYN  
    ENDIF  
    EETIME = TIMESEC
    EETIME=EETIME/86400._8
    NACTIVE=LA-1
    WRITE(NSTR,'(I3.3)') NSS
    FNAME1 =OUTDIR//'HFRE_WC_'//TRIM(NSTR)//'.OUT'

    IF( NCAL == 1 )THEN
      INQUIRE(UNIT=UHFR,OPENED=STAT)
      IF( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      OPEN(UHFR,FILE=FNAME1,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      CLOSE(UHFR,STATUS='DELETE')
      OPEN(UHFR,FILE=FNAME1,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      !------ Header contains Support Information for Output Extraction ------
      VER=10300
      HSIZE=56+24*NPNT(NSS)             ! Header Size = 2*4+7*4+4*4+1*4+NPNT(NSS)*(2*4+2*8)
      WRITE(UHFR) VER,INT(HSIZE,4)      ! File Format Version and Header Size
      WRITE(UHFR) (ISTRAN(I),I=1,7)
      WRITE(UHFR) KC,NSED,NSND,NTOX
      WRITE(UHFR) INT(NPNT(NSS),4)             ! Number of stations
      DO NP=1,NPNT(NSS) 
        WRITE(UHFR) INT(HFREGRP(NSS)%ICEL(NP),4),INT(HFREGRP(NSS)%JCEL(NP),4),HFREGRP(NSS)%XCEL(NP),HFREGRP(NSS)%YCEL(NP)
      ENDDO
      !------ Header End ------
      CLOSE(UHFR)
      RETURN
    ENDIF

    OPEN(UHFR,FILE=FNAME1,STATUS='OLD',POSITION='APPEND',FORM='BINARY')
    CLOSE(UHFR)
    OPEN(UHFR,FILE=FNAME1,STATUS='OLD',POSITION='APPEND',FORM='BINARY')
    WRITE(UHFR) EETIME,NACTIVE
    IF( ISTRAN(6) == 1 .OR. ISTRAN(7) >= 1 )THEN
      DO NP=1,NPNT(NSS) 
        LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        WRITE(UHFR)KBT_Global(LG)
      ENDDO
    ENDIF
    
    DO NP=1,NPNT(NSS) 
      LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      IF( ISTRAN(1) == 1 ) WRITE(UHFR) (REAL(SAL_Global(LG,K),4),K=KSZ_Global(LG),KC)
      IF( ISTRAN(2) == 1 ) WRITE(UHFR) (REAL(TEM_Global(LG,K),4),K=KSZ_Global(LG),KC)
      IF( ISTRAN(3) == 1 ) WRITE(UHFR) ((REAL(DYE_Global(LG,K,MD),4),K=KSZ_Global(LG),KC),MD=1,NDYE)
      IF( ISTRAN(4) == 1 ) WRITE(UHFR) (REAL(SFL_Global(LG,K),4),K=KSZ_Global(LG),KC)
      IF( ISTRAN(5) == 1 ) WRITE(UHFR) ((REAL(TOX_Global(LG,K,NT),4),K=KSZ_Global(LG),KC),NT=1,NTOX)     
      IF( ISTRAN(6) == 1 ) WRITE(UHFR) ((REAL(SED_Global(LG,K,NS),4),K=KSZ_Global(LG),KC),NS=1,NSED)
      IF( ISTRAN(7) == 1 ) WRITE(UHFR) ((REAL(SND_Global(LG,K,NX),4),K=KSZ_Global(LG),KC),NX=1,NSND)
    ENDDO
    FLUSH(UHFR)
    CLOSE(UHFR)
  END SUBROUTINE

  SUBROUTINE HFREWQOUT(NCAL,NSS)
    !FOR WATER QUALITY
    INTEGER(IK4),INTENT(IN) :: NCAL,NSS
    INTEGER(IK4) :: VER,NWQVAR,IWQ(40)
    INTEGER :: HSIZE,NP,NACTIVE
    INTEGER :: NW,MW
    INTEGER :: L,K,LG
    
    REAL(RK4)    :: WQ
    
    REAL(8) :: EETIME

    CHARACTER(50) :: FNAME4,NSTR*3
    
    LOGICAL(4) :: STAT
    
    SAVE IWQ
    SAVE NWQVAR

    IF( ISTRAN(8) == 0 ) RETURN
  
    IF( ISDYNSTP == 0 )THEN  
      DELT=DT  
    ELSE  
      DELT=DTDYN  
    ENDIF  
    EETIME=TIMESEC

    EETIME=EETIME/86400._8
    NACTIVE=LA-1
    WRITE(NSTR,'(I3.3)') NSS
    FNAME4 =OUTDIR//'HFRE_WQ_'//TRIM(NSTR)//'.OUT'

    IF( NCAL == 1 )THEN     
      INQUIRE(UNIT=UHFR,OPENED=STAT)
      IF( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      OPEN(UHFR,FILE=FNAME4,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      CLOSE(UHFR,STATUS='DELETE')
      OPEN(UHFR,FILE=FNAME4,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      NWQVAR = NWQV
 
      !------ Header contains Support Information for Output Extraction ------
      VER=10300
      HSIZE=20 + 24*NPNT(NSS) + 4*NWQVAR    ! Header Size = 2*4+2*4+1*4 + NPNT(NSS)*(2*4+2*8) + NWQVAR*4
      WRITE(UHFR) VER,INT(HSIZE,4)          ! File Format Version and Header Size
      WRITE(UHFR) KC,NWQVAR
      WRITE(UHFR) INT(NPNT(NSS),4)          ! Number of stations
      DO NP=1,NPNT(NSS) 
        LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        WRITE(UHFR) INT(HFREGRP(NSS)%ICEL(NP),4),INT(HFREGRP(NSS)%JCEL(NP),4),HFREGRP(NSS)%XCEL(NP),HFREGRP(NSS)%YCEL(NP)
      ENDDO
      IWQ = 0
      DO MW = 1,NWQVAR
        IWQ(MW) = ISKINETICS(MW)
      ENDDO
      
      WRITE(UHFR) (IWQ(NW),NW=1,NWQVAR)
      !------ Header End ------
      CLOSE(UHFR)  
      RETURN
    ENDIF
  
    OPEN(UHFR,FILE=FNAME4,STATUS='UNKNOWN',POSITION='APPEND',FORM='BINARY')
    CLOSE(UHFR)
    OPEN(UHFR,FILE=FNAME4,STATUS='UNKNOWN',POSITION='APPEND',FORM='BINARY')
    WRITE(UHFR) EETIME
    DO NP=1,NPNT(NSS)
      LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      DO K = KSZ_Global(LG),KC
        DO NW = 1,NWQVAR
          IF( IWQ(NW) > 0 )THEN
            WQ = WQV_Global(LG,K,NW)
            WRITE(UHFR) WQ
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    FLUSH(UHFR)
    CLOSE(UHFR)   
  END SUBROUTINE

  SUBROUTINE HFRERPEMOUT(NCAL,NSS)
    INTEGER(4),INTENT(IN) :: NCAL,NSS
    INTEGER :: NP, L, LG, HSIZE
    INTEGER :: NACTIVE
    REAL(8)     :: EETIME
    CHARACTER(50) :: FNAME6,NSTR*3
    LOGICAL(4) :: STAT
    INTEGER(IK4) :: VER

    IF( ISTRAN(8) == 0 .OR. ISRPEM == 0 ) RETURN
  
    IF( ISDYNSTP == 0 )THEN  
      DELT=DT  
    ELSE  
      DELT=DTDYN  
    ENDIF  
    EETIME=TIMESEC
    EETIME=EETIME/86400._8
    NACTIVE=LA-1

    WRITE(NSTR,'(I3.3)') NSS
    FNAME6 =OUTDIR//'HFRE_RPEM_'//TRIM(NSTR)//'.OUT'

    IF( NCAL == 1 )THEN      
      INQUIRE(UNIT=UHFR,OPENED=STAT)
      IF( STAT) STOP 'OPEN UNIT ERROR: UHFR IS EXISTING'
      OPEN(UHFR,FILE=FNAME6,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      CLOSE(UHFR,STATUS='DELETE')
      OPEN(UHFR,FILE=FNAME6,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')
      !------ Header contains Support Information for Output Extraction ------
      VER=10300
      HSIZE=20 + 28*NPNT(NSS)           ! Header Size = 2*4+2*4+1*4 + NPNT(NSS)*(2*4+2*8+4)
      WRITE(UHFR) VER,INT(HSIZE,4)      ! File Format Version and Header Size
      WRITE(UHFR) NRPEM,INITRPEM
      WRITE(UHFR) INT(NPNT(NSS),4)      ! Number of stations
      DO NP=1,NPNT(NSS) 
        LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        WRITE(UHFR) INT(HFREGRP(NSS)%ICEL(NP),4), INT(HFREGRP(NSS)%JCEL(NP),4), HFREGRP(NSS)%XCEL(NP), HFREGRP(NSS)%YCEL(NP), LMASKRPEM_Global(LG)
      ENDDO
      !------ Header End ------
      CLOSE(UHFR)
      RETURN
    ENDIF

    OPEN(UHFR,FILE=FNAME6,STATUS='UNKNOWN',POSITION='APPEND',FORM='BINARY')
    CLOSE(UHFR)
    OPEN(UHFR,FILE=FNAME6,STATUS='UNKNOWN',POSITION='APPEND',FORM='BINARY')
    WRITE(UHFR) EETIME
    DO NP=1,NPNT(NSS)
      LG = LIJ_Global(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      IF( LMASKRPEM_Global(LG) .OR. INITRPEM == 0 )THEN                                                                                               
        WRITE(UHFR) REAL(WQRPS_Global(LG),4), REAL(WQRPR_Global(LG),4), REAL(WQRPE_Global(LG),4), REAL(WQRPD_Global(LG),4)
      ENDIF
    ENDDO       
    FLUSH(UHFR)
    CLOSE(UHFR)
  END SUBROUTINE
  
  
  
  SUBROUTINE Gather_High_Frequency
  !---------------------------------------------------------------------------!
  ! *** Gathers values for high frequency output
  !---------------------------------------------------------------------------!
  Use Mod_Assign_Loc_Glob_For_Write
  Use Mod_Map_Write_EE_Binary

  Implicit None

  ! *** Local variables
  Integer :: j
   
  ! *** Clear out the pointer array
  j = 0

  ! *** Water Elevation Section
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(HP,1), HP, size(HP_Global,1), HP_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(BELV,1), BELV, size(BELV_Global,1), BELV_Global)
  
  ! *** Velocity Section  
  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(U,1), size(U,2), U, size(U_Global,1), size(U_Global,2), U_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(V,1), size(V,2), V, size(V_Global,1), size(V_Global,2), V_Global)

  j = j + 1
  Call Assign_Loc_Glob_For_Write(j, size(W,1), size(W,2), W, size(W_Global,1), size(W_Global,2), W_Global)
  
  ! *** Water Column Section
  if( ISTRAN(1) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SAL,1), size(SAL,2), SAL, &
      size(SAL_Global,1), size(SAL_Global,2), SAL_Global)
  endif

  if( ISTRAN(2) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TEM,1), size(TEM,2), TEM, &
      size(TEM_Global,1), size(TEM_Global,2), TEM_Global)
  endif

  if( ISTRAN(3) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(DYE,1), size(DYE,2), size(DYE,3), DYE, &
      size(DYE_Global,1), size(DYE_Global,2), size(DYE_Global,3), DYE_Global)
  endif

  if( ISTRAN(4) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SFL,1), size(SFL,2), SFL, &
      size(SFL_Global,1), size(SFL_Global,2), SFL_Global)
  endif

  if( ISTRAN(5) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(TOX,1), size(TOX,2), size(TOX,3), TOX,&
      size(TOX_Global,1), size(TOX_Global,2), size(TOX_Global,3), TOX_Global)
  endif

  if( ISTRAN(6) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SED,1), size(SED,2), size(SED,3), SED,&
      size(SED_Global,1), size(SED_Global,2), size(SED_Global,3), SED_Global)
  endif

  if( ISTRAN(7) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(SND,1), size(SND,2), size(SND,3), SND,&
      size(SND_Global,1), size(SND_Global,2), size(SND_Global,3), SND_Global)
  endif
    
  if( ISTRAN(8) >= 1 )then
    j = j + 1
    Call Assign_Loc_Glob_For_Write(j, size(WQV,1),        size(WQV,2),        size(WQV,3),        WQV, &
                                      size(WQV_Global,1), size(WQV_Global,2), size(WQV_Global,3), WQV_Global)
  endif

  num_arrays_to_write_out = j

  ! *** Call routine that maps, gathers, and sorts to produce the final Global value
  Call Handle_Calls_MapGatherSort(num_arrays_to_write_out)
  
  END SUBROUTINE Gather_High_Frequency

END MODULE
