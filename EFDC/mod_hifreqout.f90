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
      L = LIJ(HFREGRP(NS)%ICEL(NP),HFREGRP(NS)%JCEL(NP))
      LG = Map2Global(L).LG
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
        L = LIJ(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        LG = Map2Global(L).LG
        WRITE(UHFR)KBT_Global(LG)
      ENDDO
    ENDIF
    
    DO NP=1,NPNT(NSS) 
      L = LIJ(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      LG = Map2Global(L).LG
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
        L = LIJ(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
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
      L = LIJ(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      LG = Map2Global(L).LG
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
        L = LIJ(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
        LG = Map2Global(L).LG
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
      L = LIJ(HFREGRP(NSS)%ICEL(NP),HFREGRP(NSS)%JCEL(NP))
      LG = Map2Global(L).LG
      IF( LMASKRPEM_Global(LG) .OR. INITRPEM == 0 )THEN                                                                                               
        WRITE(UHFR) REAL(WQRPS_Global(LG),4), REAL(WQRPR_Global(LG),4), REAL(WQRPE_Global(LG),4), REAL(WQRPD_Global(LG),4)
      ENDIF
    ENDDO       
    FLUSH(UHFR)
    CLOSE(UHFR)
  END SUBROUTINE

END MODULE
