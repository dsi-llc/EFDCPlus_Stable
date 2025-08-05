! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE CALCSERMOD
! ***PURPOSE:
! ***INTERPOLATION FOR TIME SERIES
! ***
! ***DATE: JAN 2014

use GLOBAL 
USE INFOMOD
Use Variables_WQ
USE Variables_MPI

implicit none

contains
  
  SUBROUTINE CALCSER
    ! CHANGE RECORD
    ! ***SUBROUTINE CALPSER UPDATES TIME VARIABLE SALINITY, TEMPERATURE
    ! ***DYE, SEDIMENT, AND SHELL FISH LARVAE
    ! ***BOUNDARY CONDITIONS AND INFLOW CONCENTRATIONS

    integer :: NS, K, NT, NTT, M1, M2, NW, NC, MD, MS                                                                          
    real :: TIME, TDIFF, WTM1, WTM2
    real(RKD) :: TIMESEC2      !< Current simulation time in seconds, adjusted for time stepping option

    ! *** Adjust current time based on time stepping option
    if( ISDYNSTP == 0 )then
      if( ISTL == 2 )then
        TIMESEC2 = TIMESEC - DT/2.    ! *** 3TL when ISTL = 2 and 2TL without dynamic timestepping
      else
        TIMESEC2 = TIMESEC - DT       ! *** 3TL when ISTL = 3
      endif
    else
      TIMESEC2 = TIMESEC              ! *** 2TL dynamic time stepping
    endif 
    
    ! *** INITIALIZE NULL SERIES CONCENTRATIONS                                                                             
    NTT = 3 + NDYM + NTOX + NSED + NSND
    do NT = 1,NTT
      CQWRSERT(0,NT) = 0.
      do K = 1,KC
        CSERT(K,0,NT) = 0.
      enddo
    enddo

    ! *** CONCENTRATION SERIES INTERPOLATION FOR SAL,TEM,DYE,SFL
    do NC = 1,4
      if( ISTRAN(NC) == 0 ) CYCLE
      do NS = 1,NCSER(NC)
        TIME = TIMESEC2/TCCSER(NS,NC)
      
        if( NC == 1 )then
          call LIN_INTER_COEF(TIME,TSSAL(NS).TIM,NS,NC,M1,M2,WTM1,WTM2)
          do K = 1,KC
            CSERT(K,NS,NC) = WTM1*TSSAL(NS).VAL(M1,K)+WTM2*TSSAL(NS).VAL(M2,K)
          enddo       
        
        elseif( NC == 2 )then
          call LIN_INTER_COEF(TIME,TSTEM(NS).TIM,NS,NC,M1,M2,WTM1,WTM2)
          do K = 1,KC
            CSERT(K,NS,NC) = WTM1*TSTEM(NS).VAL(M1,K)+WTM2*TSTEM(NS).VAL(M2,K)
          enddo      
        
        elseif( NC == 4 )then
          NW = 3 + NDYM
          call LIN_INTER_COEF(TIME,TSSFL(NS).TIM,NS,NW,M1,M2,WTM1,WTM2)
          do K = 1,KC
            CSERT(K,NS,NW) = WTM1*TSSFL(NS).VAL(M1,K)+WTM2*TSSFL(NS).VAL(M2,K)
          enddo   
        
        endif       
      enddo
    enddo
  
    ! *** CONCENTRATION SERIES INTERPOLATION FOR DYE
    if( ISTRAN(3) >= 1 )then
      NC = 3
      do MD = 1,NDYE
        MS = MSVDYE(MD)
        do NS = 1,NCSER(NC)
          TIME = TIMESEC2/TCCSER(NS,NC)

          call LIN_INTER_COEF(TIME,TSDYE(NS,MD).TIM,NS,NC,M1,M2,WTM1,WTM2)
          do K = 1,KC
            CSERT(K,NS,MS) = WTM1*TSDYE(NS,MD).VAL(M1,K) + WTM2*TSDYE(NS,MD).VAL(M2,K)
          enddo
        enddo
      enddo
    endif
        
    ! *** CONCENTRATION SERIES INTERPOLATION FOR TOX
    if( ISTRAN(5) >= 1 )then
      NC = 5
      do NT = 1,NTOX
        MS = MSVTOX(NT)
        do NS = 1,NCSER(NC)
          TIME = TIMESEC2/TCCSER(NS,NC)
        
          call LIN_INTER_COEF(TIME,TSTOX(NS,NT).TIM,NS,NC,M1,M2,WTM1,WTM2)       
          do K = 1,KC
            CSERT(K,NS,MS) = WTM1*TSTOX(NS,NT).VAL(M1,K) + WTM2*TSTOX(NS,NT).VAL(M2,K)
          enddo
        enddo
      
        ! *** ATMOSPHERIC DEPOSITION SERIES
        TIME = REAL(TIMEDAY,4)
        if( TOXDEP(NT).ITXDRY > 0 )then
          call LIN_INTER_COEF_DEP(TIME,TXDRYSER(1).TIM,0,M1,M2,WTM1,WTM2)
          TOXDEP(NT).TXDRYCUR = WTM1*TXDRYSER(1).VAL(M1,NT) + WTM2*TXDRYSER(1).VAL(M2,NT)
        endif
        if( TOXDEP(NT).ITXWET > 0 )then
          call LIN_INTER_COEF_DEP(TIME,TXWETSER(1).TIM,1,M1,M2,WTM1,WTM2)
          TOXDEP(NT).TXWETCUR = WTM1*TXWETSER(1).VAL(M1,NT) + WTM2*TXWETSER(1).VAL(M2,NT)
        endif
        
      enddo
    endif
  
    ! *** CONCENTRATION SERIES INTERPOLATION FOR SED
    if( ISTRAN(6) >= 1 )then
      NC = 6
      do NT = 1,NSED
        MS = MSVSED(NT)
        do NS = 1,NCSER(NC)
          TIME = TIMESEC2/TCCSER(NS,NC)
        
          call LIN_INTER_COEF(TIME,TSSED(NS,NT).TIM,NS,NC,M1,M2,WTM1,WTM2)       
          do K = 1,KC
            CSERT(K,NS,MS) = WTM1*TSSED(NS,NT).VAL(M1,K) + WTM2*TSSED(NS,NT).VAL(M2,K)
          enddo
        enddo
      enddo
    endif
  
    ! *** CONCENTRATION SERIES INTERPOLATION FOR SND
    if( ISTRAN(7) >= 1 )then
      NC = 7
      do NT = 1,NSND
        MS = MSVSND(NT)
        do NS = 1,NCSER(NC)
          TIME = TIMESEC2/TCCSER(NS,NC)
                
          call LIN_INTER_COEF(TIME,TSSND(NS,NT).TIM,NS,NC,M1,M2,WTM1,WTM2)       
          do K = 1,KC
            CSERT(K,NS,MS) = WTM1*TSSND(NS,NT).VAL(M1,K) + WTM2*TSSND(NS,NT).VAL(M2,K)
          enddo
        enddo
      enddo
    endif

    ! *** CONCENTRATION SERIES INTERPOLATION FOR WATER QUALITY
    if( ISTRAN(8) >= 1 )then   
      NC = 8
      do NW = 1,NWQV
        if( ISTRWQ(NW) > 0 )then
          MS = MSVWQV(NW)
          do NS = 1,NCSER(NC)
            TIME = TIMESEC2/TCCSER(NS,NC)
        
            call LIN_INTER_COEF(TIME, TSWQ(NS,NW).TIM, NS, NC, M1, M2, WTM1, WTM2)        
            do K = 1,KC
              CSERT(K,NS,MS) = WTM1*TSWQ(NS,NW).VAL(M1,K) + WTM2*TSWQ(NS,NW).VAL(M2,K)
            enddo
          enddo
        endif
      enddo
    endif
  
    ! *** WRITE DIAGNOSTIC FILE FOR CSER INTERPOLATION
    if( ISDIQ >= 1 .and. N == 1 .and. DEBUG )then
      open(1,FILE = OUTDIR//'CDIAG.OUT',STATUS = 'UNKNOWN')
      close(1,STATUS = 'DELETE')
      open(1,FILE = OUTDIR//'CDIAG.OUT',STATUS = 'UNKNOWN')
      do NC = 1,NTT
        write(1,1001)NC
        do NS = 1,NCSER(NC)
          write(1,1002)NS,(CSERT(K,NS,NC),K = 1,KC)
        enddo
      enddo
      close(1)
    endif
    1001 FORMAT(/' TRANSPORT VARIABLE ID  = ',I5/)
    1002 FORMAT(I5,2X,12E12.4)
 
    ! *** SHELL FISH LARVAE BEHAVIOR TIME SERIES INTERPOLATION
    if( ISTRAN(4) >= 1 )then
      TIME = TIMESEC2/TCSFSER
      
      M1 = MSFTLST
      do while (.true.)
        M2 = M1+1
        if( TIME > TSFSER(M2) )then
          M1 = M2
        else
          MSFTLST = M1
          EXIT
        endif
      enddo
    
      TDIFF = TSFSER(M2)-TSFSER(M1)
      WTM1 = (TSFSER(M2)-TIME)/TDIFF
      WTM2 = (TIME-TSFSER(M1))/TDIFF
      RKDSFLT = WTM1*RKDSFL(M1)+WTM2*RKDSFL(M2)
      WSFLSTT = WTM1*WSFLST(M1)+WTM2*WSFLST(M2)
      WSFLSMT = WTM1*WSFLSM(M1)+WTM2*WSFLSM(M2)
      DSFLMNT = WTM1*DSFLMN(M1)+WTM2*DSFLMN(M2)
      DSFLMXT = WTM1*DSFLMX(M1)+WTM2*DSFLMX(M2)
      SFNTBET = WTM1*SFNTBE(M1)+WTM2*SFNTBE(M2)
      SFATBTT = WTM1*SFATBT(M1)+WTM2*SFATBT(M2)
    endif

    ! *** CONCENTRATION SERIES INTERPOLATION FOR BANK EROSION
    if( ISTRAN(6) >= 1 .or. ISTRAN(7) >= 1 )then
      if( ISBKERO == 1 )then
        do NS = 1,NBESER
          TIME = TIMESEC2/TCBESER(NS)

          M1 = MBETLAST(NS)
          do while (.true.)
            M2 = M1+1
            if( TIME > TBESER(M2,NS) )then
              M1 = M2
            else
              MBETLAST(NS) = M1
              EXIT
            endif
          enddo
        
          TDIFF = TBESER(M2,NS)-TBESER(M1,NS)
          WTM1 = (TBESER(M2,NS)-TIME)/TDIFF
          WTM2 = (TIME-TBESER(M1,NS))/TDIFF
          BESERT(NS) = WTM1*BESER(M1,NS)+WTM2*BESER(M2,NS)
          FWCBESERT(NS) = WTM1*FWCBESER(M1,NS)+WTM2*FWCBESER(M2,NS)

        enddo
      endif
    endif
  END SUBROUTINE

  SUBROUTINE LIN_INTER_COEF( TIME, TDAT, NS, NC, M1, M2, WTM1, WTM2)
    integer,intent(IN ) :: NS, NC
    real,   intent(IN ) :: TIME, TDAT(:)
    integer,intent(OUT) :: M1, M2
    real,   intent(OUT) :: WTM1, WTM2
    real   :: TDIFF

    M2 = MTSCLAST(NS,NC) 
    do while (TIME > TDAT(M2))   
      M2 = M2+1
      if( M2 > MCSER(NS,NC) )then
        M2 = MCSER(NS,NC)
        EXIT
      endif    
    enddo
    MTSCLAST(NS,NC) = M2  
    M1 = M2-1
    TDIFF= TDAT(M2)-TDAT(M1)     
    WTM1 = (TDAT(M2)-TIME)/TDIFF 
    WTM2 = (TIME-TDAT(M1))/TDIFF 
    
  END SUBROUTINE
  
  SUBROUTINE LIN_INTER_COEF_DEP(TIME,TDAT,ITYPE,M1,M2,WTM1,WTM2)
    integer,intent(IN)  :: ITYPE
    real,   intent(IN)  :: TIME, TDAT(:)
    integer,intent(OUT) :: M1, M2
    real,   intent(OUT) :: WTM1, WTM2
    integer :: NPTS
    real    :: TDIFF

    if( ITYPE == 0 )then
      M2   = TOXDEP(1).ITDRY    ! *** DRY DEPOSITION SERIES
      NPTS = TXDRYSER(1).NREC   ! *** NUMBER OF DATA POINTS IN SERIES
    else
      M2   = TOXDEP(1).ITWET    ! *** WET DEPOSITION SERIES
      NPTS = TXWETSER(1).NREC   ! *** NUMBER OF DATA POINTS IN SERIES
    endif

    do while (TIME > TDAT(M2))
      M2 = M2 + 1
      if( M2 > NPTS )then
        M2 = NPTS
        EXIT
      endif    
    enddo
    if( ITYPE == 0 )then
      TOXDEP(1).ITDRY = M2      ! *** DRY DEPOSITION SERIES
    else
      TOXDEP(1).ITWET = M2      ! *** WET DEPOSITION SERIES
    endif

    M1 = M2-1
    TDIFF= TDAT(M2)-TDAT(M1)     
    WTM1 = (TDAT(M2)-TIME)/TDIFF 
    WTM2 = (TIME-TDAT(M1))/TDIFF 
    
  END SUBROUTINE
  
END MODULE
