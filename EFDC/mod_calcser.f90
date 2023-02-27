! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE CALCSERMOD
! ***PURPOSE:
! ***INTERPOLATION FOR TIME SERIES
! ***
! ***DATE: JAN 2014

USE GLOBAL 
USE INFOMOD
Use Variables_WQ
USE Variables_MPI

IMPLICIT NONE

CONTAINS
  
  SUBROUTINE CALCSER
    ! CHANGE RECORD
    ! ***SUBROUTINE CALPSER UPDATES TIME VARIABLE SALINITY, TEMPERATURE
    ! ***DYE, SEDIMENT, AND SHELL FISH LARVAE
    ! ***BOUNDARY CONDITIONS AND INFLOW CONCENTRATIONS

    INTEGER :: NS, K, NT, NTT, M1, M2, NW, NC, MD, MS                                                                          
    REAL :: TIME, TDIFF, WTM1, WTM2
    
    ! *** INITIALIZE NULL SERIES CONCENTRATIONS                                                                             
    NTT = 3 + NDYM + NTOX + NSED + NSND
    DO NT=1,NTT
      CQWRSERT(0,NT)=0.
      DO K=1,KC
        CSERT(K,0,NT)=0.
      ENDDO
    ENDDO

    ! *** CONCENTRATION SERIES INTERPOLATION FOR SAL,TEM,DYE,SFL
    DO NC=1,4
      IF( ISTRAN(NC) == 0 ) CYCLE
      DO NS=1,NCSER(NC)
        TIME = TIMESEC/TCCSER(NS,NC)
      
        IF( NC == 1 )THEN
          CALL LIN_INTER_COEF(TIME,TSSAL(NS).TIM,NS,NC,M1,M2,WTM1,WTM2)
          DO K=1,KC
            CSERT(K,NS,NC)=WTM1*TSSAL(NS).VAL(M1,K)+WTM2*TSSAL(NS).VAL(M2,K)
          ENDDO       
        
        ELSEIF( NC == 2 )THEN
          CALL LIN_INTER_COEF(TIME,TSTEM(NS).TIM,NS,NC,M1,M2,WTM1,WTM2)
          DO K=1,KC
            CSERT(K,NS,NC)=WTM1*TSTEM(NS).VAL(M1,K)+WTM2*TSTEM(NS).VAL(M2,K)
          ENDDO      
        
        ELSEIF( NC == 4 )THEN
          NW = 3 + NDYM
          CALL LIN_INTER_COEF(TIME,TSSFL(NS).TIM,NS,NW,M1,M2,WTM1,WTM2)
          DO K=1,KC
            CSERT(K,NS,NW)=WTM1*TSSFL(NS).VAL(M1,K)+WTM2*TSSFL(NS).VAL(M2,K)
          ENDDO   
        
        ENDIF       
      ENDDO
    ENDDO
  
    ! *** CONCENTRATION SERIES INTERPOLATION FOR DYE
    IF( ISTRAN(3) >= 1 )THEN
      NC = 3
      DO MD=1,NDYE
        MS = MSVDYE(MD)
        DO NS=1,NCSER(NC)
          TIME = TIMESEC/TCCSER(NS,NC)

          CALL LIN_INTER_COEF(TIME,TSDYE(NS,MD).TIM,NS,NC,M1,M2,WTM1,WTM2)
          DO K=1,KC
            CSERT(K,NS,MS) = WTM1*TSDYE(NS,MD).VAL(M1,K) + WTM2*TSDYE(NS,MD).VAL(M2,K)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
        
    ! *** CONCENTRATION SERIES INTERPOLATION FOR TOX
    IF( ISTRAN(5) >= 1 )THEN
      NC = 5
      DO NT=1,NTOX
        MS = MSVTOX(NT)
        DO NS=1,NCSER(NC)
          TIME = TIMESEC/TCCSER(NS,NC)
        
          CALL LIN_INTER_COEF(TIME,TSTOX(NS,NT).TIM,NS,NC,M1,M2,WTM1,WTM2)       
          DO K=1,KC
            CSERT(K,NS,MS) = WTM1*TSTOX(NS,NT).VAL(M1,K) + WTM2*TSTOX(NS,NT).VAL(M2,K)
          ENDDO
        ENDDO
      
        ! *** ATMOSPHERIC DEPOSITION SERIES
        TIME = REAL(TIMEDAY,4)
        IF( TOXDEP(NT).ITXDRY > 0 )THEN
          CALL LIN_INTER_COEF_DEP(TIME,TXDRYSER(1).TIM,0,M1,M2,WTM1,WTM2)
          TOXDEP(NT).TXDRYCUR = WTM1*TXDRYSER(1).VAL(M1,NT) + WTM2*TXDRYSER(1).VAL(M2,NT)
        ENDIF
        IF( TOXDEP(NT).ITXWET > 0 )THEN
          CALL LIN_INTER_COEF_DEP(TIME,TXWETSER(1).TIM,1,M1,M2,WTM1,WTM2)
          TOXDEP(NT).TXWETCUR = WTM1*TXWETSER(1).VAL(M1,NT) + WTM2*TXWETSER(1).VAL(M2,NT)
        ENDIF
        
      ENDDO
    ENDIF
  
    ! *** CONCENTRATION SERIES INTERPOLATION FOR SED
    IF( ISTRAN(6) >= 1 )THEN
      NC = 6
      DO NT=1,NSED
        MS = MSVSED(NT)
        DO NS=1,NCSER(NC)
          TIME = TIMESEC/TCCSER(NS,NC)
        
          CALL LIN_INTER_COEF(TIME,TSSED(NS,NT).TIM,NS,NC,M1,M2,WTM1,WTM2)       
          DO K=1,KC
            CSERT(K,NS,MS) = WTM1*TSSED(NS,NT).VAL(M1,K) + WTM2*TSSED(NS,NT).VAL(M2,K)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  
    ! *** CONCENTRATION SERIES INTERPOLATION FOR SND
    IF( ISTRAN(7) >= 1 )THEN
      NC = 7
      DO NT=1,NSND
        MS = MSVSND(NT)
        DO NS=1,NCSER(NC)
          TIME = TIMESEC/TCCSER(NS,NC)
                
          CALL LIN_INTER_COEF(TIME,TSSND(NS,NT).TIM,NS,NC,M1,M2,WTM1,WTM2)       
          DO K=1,KC
            CSERT(K,NS,MS) = WTM1*TSSND(NS,NT).VAL(M1,K) + WTM2*TSSND(NS,NT).VAL(M2,K)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! *** CONCENTRATION SERIES INTERPOLATION FOR WATER QUALITY
    IF( ISTRAN(8) >= 1 )THEN   
      NC = 8
      DO NW=1,NWQV
        IF( ISTRWQ(NW) > 0 )THEN
          MS = MSVWQV(NW)
          DO NS=1,NCSER(NC)
            TIME = TIMESEC/TCCSER(NS,NC)
        
            CALL LIN_INTER_COEF(TIME, TSWQ(NS,NW).TIM, NS, NC, M1, M2, WTM1, WTM2)        
            DO K=1,KC
              CSERT(K,NS,MS) = WTM1*TSWQ(NS,NW).VAL(M1,K) + WTM2*TSWQ(NS,NW).VAL(M2,K)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  
    ! *** WRITE DIAGNOSTIC FILE FOR CSER INTERPOLATION
    IF( ISDIQ >= 1 .AND. N == 1 .AND. DEBUG )THEN
      OPEN(1,FILE=OUTDIR//'CDIAG.OUT',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'CDIAG.OUT',STATUS='UNKNOWN')
      DO NC=1,NTT
        WRITE(1,1001)NC
        DO NS=1,NCSER(NC)
          WRITE(1,1002)NS,(CSERT(K,NS,NC),K=1,KC)
        ENDDO
      ENDDO
      CLOSE(1)
    ENDIF
    1001 FORMAT(/' TRANSPORT VARIABLE ID =',I5/)
    1002 FORMAT(I5,2X,12E12.4)
 
    ! *** SHELL FISH LARVAE BEHAVIOR TIME SERIES INTERPOLATION
    IF( ISTRAN(4) >= 1 )THEN
      IF( ISTL == 2 )THEN
        IF( ISDYNSTP == 0 )THEN
          TIME=DT*(FLOAT(N)-0.5)/TCSFSER +TBEGIN*(TCON/TCSFSER)
        ELSE
          TIME=TIMESEC/TCSFSER
        ENDIF
      ELSE
        IF( ISDYNSTP == 0 )THEN
          TIME=DT*FLOAT(N-1)/TCSFSER +TBEGIN*(TCON/TCSFSER)
        ELSE
          TIME=TIMESEC/TCSFSER
        ENDIF
      ENDIF
      
      M1=MSFTLST
      DO WHILE(1)
        M2=M1+1
        IF( TIME > TSFSER(M2) )THEN
          M1=M2
        ELSE
          MSFTLST=M1
          EXIT
        ENDIF
      ENDDO
    
      TDIFF=TSFSER(M2)-TSFSER(M1)
      WTM1=(TSFSER(M2)-TIME)/TDIFF
      WTM2=(TIME-TSFSER(M1))/TDIFF
      RKDSFLT=WTM1*RKDSFL(M1)+WTM2*RKDSFL(M2)
      WSFLSTT=WTM1*WSFLST(M1)+WTM2*WSFLST(M2)
      WSFLSMT=WTM1*WSFLSM(M1)+WTM2*WSFLSM(M2)
      DSFLMNT=WTM1*DSFLMN(M1)+WTM2*DSFLMN(M2)
      DSFLMXT=WTM1*DSFLMX(M1)+WTM2*DSFLMX(M2)
      SFNTBET=WTM1*SFNTBE(M1)+WTM2*SFNTBE(M2)
      SFATBTT=WTM1*SFATBT(M1)+WTM2*SFATBT(M2)
    ENDIF

    ! *** CONCENTRATION SERIES INTERPOLATION FOR BANK EROSION
    IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
      IF( ISBKERO == 1 )THEN
        DO NS=1,NBESER
          IF( ISTL == 2 )THEN
            IF( ISDYNSTP == 0 )THEN
              TIME=DT*(FLOAT(N)-0.5)/TCBESER(NS) +TBEGIN*(TCON/TCBESER(NS))
            ELSE
              TIME=TIMESEC/TCBESER(NS)
            ENDIF
          ELSE
            IF( ISDYNSTP == 0 )THEN
              TIME=DT*FLOAT(N-1)/TCBESER(NS) +TBEGIN*(TCON/TCBESER(NS))
            ELSE
              TIME=TIMESEC/TCBESER(NS)
            ENDIF
          ENDIF

          M1=MBETLAST(NS)
          DO WHILE(1)
            M2=M1+1
            IF( TIME > TBESER(M2,NS) )THEN
              M1=M2
            ELSE
              MBETLAST(NS)=M1
              EXIT
            ENDIF
          ENDDO
        
          TDIFF=TBESER(M2,NS)-TBESER(M1,NS)
          WTM1=(TBESER(M2,NS)-TIME)/TDIFF
          WTM2=(TIME-TBESER(M1,NS))/TDIFF
          BESERT(NS)=WTM1*BESER(M1,NS)+WTM2*BESER(M2,NS)
          FWCBESERT(NS)=WTM1*FWCBESER(M1,NS)+WTM2*FWCBESER(M2,NS)

        ENDDO
      ENDIF
    ENDIF
  END SUBROUTINE

  SUBROUTINE LIN_INTER_COEF(TIME,TDAT,NS,NC,M1,M2,WTM1,WTM2)
    INTEGER,INTENT(IN ) :: NS,NC
    REAL,   INTENT(IN ) :: TIME,TDAT(:)
    INTEGER,INTENT(OUT) :: M1,M2
    REAL,   INTENT(OUT) :: WTM1,WTM2
    REAL   :: TDIFF

    M2 = MTSCLAST(NS,NC) 
    DO WHILE (TIME > TDAT(M2))   
      M2 = M2+1
      IF( M2 > MCSER(NS,NC) )THEN
        M2 = MCSER(NS,NC)
        EXIT
      ENDIF    
    END DO
    MTSCLAST(NS,NC) = M2  
    M1 = M2-1
    TDIFF= TDAT(M2)-TDAT(M1)     
    WTM1 = (TDAT(M2)-TIME)/TDIFF 
    WTM2 = (TIME-TDAT(M1))/TDIFF 
    
  END SUBROUTINE
  
  SUBROUTINE LIN_INTER_COEF_DEP(TIME,TDAT,ITYPE,M1,M2,WTM1,WTM2)
    INTEGER,INTENT(IN)  :: ITYPE
    REAL,   INTENT(IN)  :: TIME, TDAT(:)
    INTEGER,INTENT(OUT) :: M1, M2
    REAL,   INTENT(OUT) :: WTM1, WTM2
    INTEGER :: NPTS
    REAL    :: TDIFF

    IF( ITYPE == 0 )THEN
      M2   = TOXDEP(1).ITDRY    ! *** DRY DEPOSITION SERIES
      NPTS = TXDRYSER(1).NREC   ! *** NUMBER OF DATA POINTS IN SERIES
    ELSE
      M2   = TOXDEP(1).ITWET    ! *** WET DEPOSITION SERIES
      NPTS = TXWETSER(1).NREC   ! *** NUMBER OF DATA POINTS IN SERIES
    ENDIF

    DO WHILE (TIME > TDAT(M2))
      M2 = M2 + 1
      IF( M2 > NPTS )THEN
        M2 = NPTS
        EXIT
      ENDIF    
    END DO
    IF( ITYPE == 0 )THEN
      TOXDEP(1).ITDRY = M2      ! *** DRY DEPOSITION SERIES
    ELSE
      TOXDEP(1).ITWET = M2      ! *** WET DEPOSITION SERIES
    ENDIF

    M1 = M2-1
    TDIFF= TDAT(M2)-TDAT(M1)     
    WTM1 = (TDAT(M2)-TIME)/TDIFF 
    WTM2 = (TIME-TDAT(M1))/TDIFF 
    
  END SUBROUTINE
  
END MODULE
