! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE WQ_RPEM_MODULE

! *** EFDC_DSI  ROOTED PLANT AND EPIPHYTE MODEL (RPEM)
! *** MODULE: WQ_RPEM_MODULE

! CHANGE RECORD 
! DATE MODIFIED     BY               DESCRIPTION        
!-- ------------------------------------------------------------------
! 2012-02-25        Paul M. Craig    Fixed Temperature Dependency Table
! 2011-11-29        PAUL M. CRAIG &  OMP'd MODEL AND ADDED COMPUTATIONAL BYPASSES
!                   DANG CHUNG       RESTRUCTURED TO F90 AND USE OF MODULES 

  Use Variables_WQ
  USE WQ_DIAGENESIS
  
IMPLICIT NONE

  ! *** INTEGER ARRAYS                                                                                                     
  INTEGER, ALLOCATABLE :: IRPEMTS(:),JRPEMTS(:),LRPEMTS(:)

  ! *** INTEGER SCALARS                                                                                                    
  INTEGER   :: ISRPEM,INITRPEM,IRPEMWC,IRPEME,IJRPRS
  INTEGER   :: ISTOXRPE,JSRPEM,ISRPEMSPAC,ISRPEMSPFR
  INTEGER   :: ISRPEMTIME,ISRPEMTIFR,ISRPEMTILC,NSRPEMSPFR
  INTEGER   :: NCRPEMRST,IRPEMSD
  INTEGER   :: NRPEM, NRPEMEE
  INTEGER   :: NRPEMSTEPS

  ! *** REAL SCALARS                                                                                                       
  REAL       :: HRPEMIC,RPSO,RPRO,RPEO,RPDO,RPSOC,RPEOC,RPSNC
  REAL       :: RPRNC,RPENC,RPSPC,RPRPC,RPEPC,PMRPS,FPRPR,RMRPS,RLRPS,FRPSD
  REAL       :: RKHNPRPS,RMRPR,RLRPR,RLRPD,RJRPRSC,RKRPORS,ROSR,RKRPRS,RISSS
  REAL       :: PMRPE,RMRPE,RLRPE,RKHNPRPE,RKHNRPS,RKHNRPR,RKHNRPE,RKHPRPS
  REAL       :: RKHPRPR,RKHPRPE,TP1RPS,TP2RPS,RKTP1RPS,RKTP2RPS,TP1RPE
  REAL       :: RKTP1RPE,RKTP2RPE,HRPS,HOPT,RKERPE,CCHLRPE,RISSOM,RISSOEM
  REAL       :: STOXS,STOXE,TR1RPS,TR2RPS,RKTR1RPS,RKTR2RPS,TR1RPR,TR2RPR
  REAL       :: RKTR1RPR,RKTR2RPR,TR1RPE,TR2RPE,RKTR1RPE,RKTR2RPE,FCRRPS
  REAL       :: FCLRPS,FCDRPS,FCRLRPS,FCLLRPS,FCDLRPS,FCRRPR,FCLRPR,FCDRPR
  REAL       :: FCRLRPR,FCLLRPR,FCDLRPR,FCRRPE,FCLRPE,FCDRPE,FCRLRPE
  REAL       :: FCLLRPE,FCDLRPE,FCRLRPD,FCLLRPD,FCDLRPD,FPRRPS,FPLRPS
  REAL       :: FPDRPS,FPIRPS,FPRLRPS,FPLLRPS,FPDLRPS,FPILRPS,FPRRPR,FPLRPR
  REAL       :: FPDRPR,FPIRPR,FPRLRPR,FPLLRPR,FPDLRPR,FPILRPR,FPRRPE,FPLRPE
  REAL       :: FPDRPE,FPIRPE,FPRLRPE,FPLLRPE,FPDLRPE,FPILRPE,FPRLRPD
  REAL       :: FPLLRPD,FPDLRPD,FPILRPD,FNRRPS,FNLRPS,FNDRPS,FNIRPS,FNRLRPS
  REAL       :: FNLLRPS,FNDLRPS,FNILRPS,FNRRPR,FNLRPR,FNDRPR,FNIRPR,FNRLRPR
  REAL       :: FNLLRPR,FNDLRPR,FNILRPR,FNRRPE,FNLRPE,FNDRPE,FNIRPE,FNRLRPE
  REAL       :: FNLLRPE,FNDLRPE,FNILRPE,FNRLRPD,FNLLRPD,FNDLRPD,FNILRPD
  REAL       :: FRPRRPG1,FRPRRPG2,FRPRRPG3
  REAL       :: RKSH,RKHI

  REAL,TARGET, ALLOCATABLE :: WQRPS(:),     WQRPR(:),     WQRPE(:),     WQRPD(:)
  REAL, ALLOCATABLE :: PRPS(:),      RRPS(:),      RRPR(:),      PRPE(:),   RRPE(:)
  REAL, ALLOCATABLE :: WQRPSR(:),    WQRPSL(:),    WQRPER(:),    WQRPEL(:), WQRPDL(:)
  REAL, ALLOCATABLE :: WQRPSRP(:),   WQRPSLP(:),   WQRPERP(:),   WQRPELP(:)
  REAL, ALLOCATABLE :: WQRPDLP(:),   WQRPSRN(:),   WQRPSLN(:),   WQRPERN(:)
  REAL, ALLOCATABLE :: WQRPELN(:),   WQRPDLN(:),   WQRPRR(:),    WQRPRL(:)
  REAL, ALLOCATABLE :: WQRPRRP(:),   WQRPRLP(:),   WQRPRRN(:),   WQRPRLN(:)
  REAL, ALLOCATABLE :: FRPSPW(:),    FRPSNW(:),    PNRPS(:),     PNRPE(:),  RJRPRS(:)   
  REAL, ALLOCATABLE :: XLIMTPRPS(:), XLIMTPRPE(:), XLIMTRRPS(:), XLIMTRRPE(:)
  REAL, ALLOCATABLE :: XLIMTRRPR(:), XLIMNRPS(:),  XLIMNRPE(:),  XLIMLRPS(:)
  REAL, ALLOCATABLE :: XLIMLRPE(:),  RISS(:),      RISSO(:),     RISSOE(:)
  REAL, ALLOCATABLE :: RPEMTPRPS(:), RPEMTPRPE(:), RPEMTRRPS(:)
  REAL, ALLOCATABLE :: RPEMTRRPE(:), RPEMTRRPR(:)

  CONTAINS

!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine CAL_RPEM
!
!> @details  SIMULATES ROOTED PLANTS (SHOOTS AND ROOTS) 
!>           EPIPHYTES GROWING ON SHOOTS, AND SHOOT ORGANIC DETRITUS
!---------------------------------------------------------------------------!
!  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION
!---------------------------------------------------------------------------!
SUBROUTINE CAL_RPEM
                                            
  INTEGER :: L,IWQTRPEM,K,LL,LF,LP,ND,IOBC,LN                         
  REAL    :: RATION,TOP,TMPNIT,RATIOP,TMPPO4,WQAVGIO,TMP1,RATIOHP,HDRY2    
  REAL    :: RKESSAVG,ALPHATOP,ALPHABOT,TMPEXP,SOURSINK,FACIMP,DTDHWQ 
  REAL    :: TMPWAT,TMPBED,TOP1,TOP2,BOT1,BOT2,BOT3,TMPNH4S,TMPNO3S
  REAL    :: TMPNH4E,TMPNO3E,WQKESS                                                            
  REAL    :: RKESSTOP(LCM),RKESSBOT(LCM)
  REAL    :: WQBCV(NBCSOP,NWQV)
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)   :: NLRPEM   ! *** NUMBER OF RPEM ACTIVE CELLS FOR EACH LAYER BY DOMAIN
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:) :: LLRPEM   ! *** L INDEX FOR THE RPEM CELLS
  INTEGER,SAVE :: ICOUNT
  
  ! *** INPUT AND INITIALIZATION                                                                                           
  IF( JSRPEM == 1 )THEN
    NSRPEMSPFR = 0                                                                                                         
    NCRPEMRST = 0                                                                                                          
    
    ALLOCATE(NLRPEM(NDM))
    ALLOCATE(LLRPEM(LCM,NDM))
    NLRPEM = 0
    LLRPEM = 0
    ICOUNT = 0
  ENDIF
  ICOUNT = ICOUNT+1
  
  ! *** RPEM ONLY INTERACTS WITH THE BOTTOM LAYER K = KSZ(L)
  
  ! *** SAVE VALUES AT OPEN BOUNDARIES
  IF( IRPEMWC == 0 )THEN
    DO IOBC = 1,NBCSOP  
      L = LOBCS(IOBC)  
      WQBCV(IOBC,1:NWQV) = WQV(L,KSZ(L),1:NWQV)
    ENDDO  
  ENDIF
  HDRY2 = 2.*HDRY
  
  !$OMP PARALLEL DEFAULT(SHARED)
  
  ! *** OBTAIN THE ACTIVE CELL LIST
  IF( (ISDRY > 0 .AND. LADRY > 0) .OR. JSRPEM == 1 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LN,LP,L)
    DO ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = MIN(LF+LDMWET-1,LAWET)
      LN = 0
      DO LP = LF,LL
        L = LWET(LP)
        IF( LMASKRPEM(L) )THEN
          LN = LN+1
          LLRPEM(LN,ND) = L
          IF( ISICE > 2 )THEN
            IF( ICECELL(L) .AND. HP(L) < 3.*HDRY )THEN
              ! *** DEACTIVATE RPEM KINETICS IF THE CELL HAS ICE COVER AND VERY SHALLOW DEPTHS
              LN = LN-1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      NLRPEM(ND) = LN
    ENDDO
    !$OMP END DO
      
    ! *** Report only the fist instance
    IF( JSRPEM == 1 )THEN
      !$OMP SINGLE
      DO ND = 1,NDM  
        if( process_id == master_id ) WRITE(6,*) 'RPEM DOMAIN LIST (Thread,Count):', ND, NLRPEM(ND)  
        OPEN(mpi_log_unit,FILE = OUTDIR//mpi_log_file,POSITION = 'APPEND')
        WRITE(mpi_log_unit,*) 'RPEM DOMAIN LIST (Thread,Count):', ND, NLRPEM(ND) 
        CLOSE(mpi_log_unit)
      ENDDO
      !$OMP END SINGLE
    ENDIF
  ENDIF
  
  !$OMP DO PRIVATE(ND,LP,L,IWQTRPEM,K) &
  !$OMP PRIVATE(RATION,TOP,TMPNIT,RATIOP,TMPPO4,WQAVGIO,TMP1)  &
  !$OMP PRIVATE(RKESSAVG,ALPHATOP,ALPHABOT,TMPEXP,SOURSINK,FACIMP,DTDHWQ) & 
  !$OMP PRIVATE(TMPWAT,TMPBED,TOP1,TOP2,BOT1,BOT2,BOT3,TMPNH4S,TMPNO3S)   &
  !$OMP PRIVATE(TMPNH4E,TMPNO3E,WQKESS,RKESSTOP,RKESSBOT,RATIOHP)
  DO ND = 1,NDM  
    ! **********************************************************************C                                                
    ! *** SET TEMPERATURE DEPENDENCY FOR GROWTH AND RESPIRATION 
    ! *** F3(T), EQ.(13), LOOK-UP TABLE
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      IWQTRPEM = NINT((TWQ(L)-WQTDMIN)/WQTDINC)+1
      IF( IWQTRPEM < 1) IWQTRPEM = 1
      IF( IWQTRPEM > NWQTD) IWQTRPEM = NWQTD
      XLIMTPRPS(L) = RPEMTPRPS(IWQTRPEM)  ! Shoot Production (Growth)                           
      XLIMTPRPE(L) = RPEMTPRPE(IWQTRPEM)  ! Epiphyte Production (Growth)
      XLIMTRRPS(L) = RPEMTRRPS(IWQTRPEM)  ! Shoot Respiration
      XLIMTRRPE(L) = RPEMTRRPE(IWQTRPEM)  ! Epiphyte Respiration
      XLIMTRRPR(L) = RPEMTRRPR(IWQTRPEM)  ! Root Respiration
    ENDDO                                                                                                                  
   
    ! **********************************************************************C                                                
    ! *** SET NUTRIENT LIMITATION FOR PLANT SHOOT GROWTH                                                                     
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      K = KSZ(L)
      RATION = RKHNRPS/RKHNRPR                                       ! Ratio of N half saturation constants for water/bed  
      TOP = WQV(L,K,INHX) + WQV(L,K,INOX) + RATION*(SM2NH4(L)+SM2NO3(L))                                                           
      TMPNIT = TOP/(TOP+RKHNRPS)
      RATIOP = RKHPRPS/RKHPRPR
      TOP = WQPO4D(L,K)+RATIOP*SM2PO4(L)                                                                                   
      TMPPO4 = TOP/(TOP+RKHPRPS)
      XLIMNRPS(L) = MIN(TMPNIT,TMPPO4)                                ! F1(N), EQ. (6), Minimum of N or P limit
    ENDDO                                                                                                                  
   
    ! **********************************************************************C                                                
    ! *** SET NUTRIENT LIMITATION FOR EPIPHYTE GROWTH                                                                        
    IF( IRPEME > 0 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TOP = WQV(L,K,INHX) + WQV(L,K,INOX)                                                                                        
        TMPNIT = TOP/(TOP+RKHNRPE+1.E-18)                                                                                           
        TMPPO4 = WQPO4D(L,K)/(WQPO4D(L,K)+RKHPRPE+1.E-18)
        XLIMNRPE(L) = MIN(TMPNIT,TMPPO4)                                ! *** Minimum of N or P limits.  EQ. (20)
      ENDDO                                                                                                                  
    ENDIF
     
    ! **********************************************************************C                                                
    ! *** SET LIGHT LIMITATIONS                                                                                              
   
    ! *** NOTE THAT WQI0,WQI1,AND WQI2 ARE PHOTOSYNTHETIC SSW RADIATION                                                      
    !    FROM CURRENT AND PREVIOUS TWO TIME INCREMENTS IN LANGLEY/DAY                                                       
    !    AND ARE PROVIDED BY THE WATER QUALITY MODEL                                                                        
    WQAVGIO = WQCIA*WQI0+WQCIB*WQI1+WQCIC*WQI2                          ! EQ. (12)
   
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      ! *** COMPUTE TOTAL EXTINCTION COEFFICIENT

      ! *** Bottom PLUS 1 Layer
      IF( KSZ(L) /= KC )THEN
        RKESSTOP(L) = RADKE(L,KSZ(L)+1)   ! EQ. (24)
      ELSE
        RKESSTOP(L) = 0.                    ! EQ. (24)
      ENDIF
              
      ! *** Bottom Layer
      WQKESS = RADKE(L,KSZ(L))
      RKESSBOT(L) = WQKESS
    ENDDO                                                                                                                  
   
    ! *** LIGHT LIMITATION SECTION
    IF( IRPEME > 0 )THEN
      ! *** LIGHT LIMITATION FOR SHOOTS                                                                                        
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RKESSBOT(L) = RKESSBOT(L) + RKERPE*Z(L,KSZ(L))*HWQ(L)*WQRPE(L)/CCHLRPE   ! EQ. (10), INCLUDE EPIPHYTES SHADING
      ENDDO                                                                                                                  
    ENDIF
   
    IF( IWQSUN == 2 )THEN
      ! *** ASER TIME INTERVAL VARIABLE LIGHT
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSO(L) = RISSOM                                                ! Maximum value for optimum SSW for rooted plant growth (Langleys/day)
        RKESSAVG = 0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        WQAVGIO = PARADJ*2.065*RADTOP(L,KC)                              ! RADTOP INCLUDES SHADING AND ICE COVER
        ALPHATOP = -(WQAVGIO/RISSO(L)) *EXP(-RKESSTOP(L)*TMP1*HWQ(L))    ! EQ. (9), ASSUME HRPS IS WITHIN THE BOTTOM LA
        ALPHABOT = -(WQAVGIO/RISSO(L)) *EXP(-RKESSBOT(L)*HWQ(L))         ! EQ. (8)
        TMPEXP = EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPS(L) = 2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP         ! F2(I), EQ. (7)
        XLIMLRPS(L) = MIN(XLIMLRPS(L),1.0)
      ENDDO                                                                                                                  
    ELSE
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSO(L) = RISSOM                                                ! Maximum value for optimum SSW for rooted plant growth (Langleys/day)
        RKESSAVG = 0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        ALPHATOP = -(WQAVGIO/RISSO(L)) *EXP(-RKESSTOP(L)*TMP1*HWQ(L))    ! EQ. (9), ASSUME HRPS IS WITHIN THE BOTTOM LA
        ALPHABOT = -(WQAVGIO/RISSO(L)) *EXP(-RKESSBOT(L)*HWQ(L))         ! EQ. (8)
        TMPEXP = EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPS(L) = 2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP         ! F2(I), EQ. (7), BUT PRACTICALLY, IT DOES NOT
        XLIMLRPS(L) = MIN(XLIMLRPS(L),1.0)
      ENDDO                                                                                                                  
    ENDIF
    
    IF( IRPEME > 0 )THEN
      ! *** LIGHT LIMITATION FOR EPIPHYTES                                                                                     
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSOE(L) = RISSOEM                                                                                                  
        RKESSAVG = 0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        ALPHATOP = -(WQAVGIO/RISSOE(L))*EXP(-RKESSTOP(L)*TMP1*HWQ(L))   ! EQ. (21), ASSUME HRPS IS WITHIN THE BOTTOM L
        ALPHABOT = -(WQAVGIO/RISSOE(L))*EXP(-RKESSBOT(L)*HWQ(L))        ! EQ. (22)
        TMPEXP = EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPE(L) = 2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP        ! EQ. (21)
        XLIMLRPE(L) = MIN(XLIMLRPE(L),1.0)
      ENDDO                                                                                                                  
    ENDIF
  
    IF( IJRPRS == 2 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        TMP1 = 1.-Z(L,KSZ(L))
        RISS(L) = WQAVGIO*EXP(-RKESSTOP(L)*TMP1*HWQ(L))                 ! EQ. (11), PART OF IT, USED IF IJRPRS = 2 IN CAL_RPEM.INP C5
                                                                        ! HOPT IS NOT USED HERE!, CAL_RPEM.INP C9
      ENDDO                                                                                                                  
    ENDIF
     
    ! **********************************************************************C                                                
    ! *** UPDATE GROWTH AND RESPIRATION RATES                                                                                
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      
      ! *** Limitors     TEM          N/P         Light         Self Shading
      PRPS(L) = PMRPS*XLIMTPRPS(L)*XLIMNRPS(L)*XLIMLRPS(L)*EXP(-RKSH*WQRPS(L))  ! EQ. (5)   PRPS - Nutrient, light and temperature limited growth rate (/day)
                                                                                !                  Include self-shading of shoots, JI, 7/4/04
      RRPS(L) = RMRPS*XLIMTRRPS(L)                                              ! EQ. (15)  RRPS - Temperature limited shoot respiration (/day)
    ENDDO                                                                                                                  
   
    ! *** Limit RPEM with depth
    IF( ISDRY > 0 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RATIOHP = MIN( (MAX(HP(L)-HDRY,0.))/HDRY2, 1.0 )
        PRPS(L) = PRPS(L)*RATIOHP
        RRPS(L) = RRPS(L)*RATIOHP
      ENDDO    
    ENDIF
    
    ! *** ROOT RESPIRATION                                                                                                   
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      RRPR(L) = RMRPR*XLIMTRRPR(L)                                    ! EQ. (18)
    ENDDO                                                                                                                  
   
    ! *** EPIPHYTE GROWTH AND RESPIRATION                                                                                    
    IF( IRPEME > 0 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        PRPE(L) = PMRPE*XLIMTPRPE(L)*XLIMNRPE(L)*XLIMLRPE(L)            ! EQ. (19), F3T*F1N*F2I
        RRPE(L) = RMRPE*XLIMTRRPE(L)                                                                                         
      ENDDO                                                                                                                  
    ENDIF
     
    ! **********************************************************************C                                                
   
    ! *** UPDATE ROOT TO SHOOT FLUX                                                                                          
    IF( IJRPRS == 0 .AND. N < 5 )THEN                      ! *** IJRPRS = 0, CONSTANT TRANSPORT RATE
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RJRPRSC
      ENDDO                                                                                                                  
    ENDIF                                                                                                                  
   
    IF( IJRPRS == 1 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RKRPORS*(ROSR*WQRPR(L)-WQRPS(L))                                                                         
      ENDDO                                                                                                                  
    ENDIF                                                                                                                  
   
    IF( IJRPRS == 2 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RKRPRS*RISS(L)/(RISS(L)+RISSS+1E-18)                                                                           
      ENDDO                                                                                                                  
    ENDIF                                                                                                                  
   
    ! **********************************************************************C                                                
    ! ***  UPDATE SHOOT, ROOT, EPIPHYTE AND DETRITUS STATE VARIABLES                                                         
   
    ! *** UPDATE SHOOTS BIOMASS                                                                                                      
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      ! ***            Growth     Respiration   Other
      SOURSINK = (1.-FPRPR)*PRPS(L) - RRPS(L) - RLRPS                     ! EQ. (1)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPS(L) = FACIMP*(WQRPS(L) + DTWQ*RJRPRS(L))
      WQRPS(L) = MAX(WQRPS(L),0.2)   ! KEEP THE "SEED", SINCE IF RPS! = 0, IT WILL NEVER GROW AGAIN, ACCORDING TO EQ.
    ENDDO                                                                                                                  
   
    ! *** UPDATE ROOTS BIOMASS                                                                                                      
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      SOURSINK = -RRPR(L)-RLRPR                                       ! EQ. (2)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPR(L) = FACIMP*(WQRPR(L)+DTWQ*FPRPR*PRPS(L)*WQRPS(L)-DTWQ*RJRPRS(L))
    ENDDO                                                                                                                  

    !------------------------------                                                                                         
    !GO TO 904  ! SKIP EPHIPHYTES AND DETRITUS
    
    ! *** UPDATE EPIPHYTES                                                                                                   
    IF( IRPEME > 0 )THEN
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        SOURSINK = PRPE(L)-RRPE(L)-RLRPE                              ! EQ. (3)
        FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
        WQRPE(L) = FACIMP*WQRPE(L)
      ENDDO                                                                                                                  
    ENDIF
     
    ! *** UPDATE SHOOT DETRITUS IN WATER COLUMN                                                                              
    SOURSINK = -RLRPD                                                    ! EQ. (4)  RLRPD-Loss rate for plant detritus at bottom of water column (/day)
    DO LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPD(L) = FACIMP*(WQRPD(L)+DTWQ*FRPSD*RLRPS*WQRPS(L))
    ENDDO                                                                                                                  
    904   CONTINUE                                                                                                          
   
    ! **********************************************************************C                                                
    ! *** CALCULATE WATER COLUMN SHOOT, EPIPHYTE AND DETRITUS                                                                
    ! *** RESPIRATION AND NON-RESPIRATION LOSSES TO WATER QUALITY ORGANIC                                                    
    ! *** MATTER STATE VARIABLES AND TO PHOSPHATE AND AMMONIA                                                                
   
    !     GO TO 901   ! SKIP COUPLING WITH NUTRIENTS AND DO IN WATER COLUMN                                    
    IF( IRPEMWC == 0 )THEN

      ! *** UPDATE SHOOT, ROOT AND DETRITUS
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        
        ! *** LOSSES TO ORGANIC CARBON
        ! *** RLRPD - Carbon respiration loss rate of detritus (1/day)
        ! *** RLRPS - Carbon non-respiration loss rate of detritus (1/day)
        ! *** FRPSD - Non-respiration
        WQRPSR(L) = RRPS(L)*WQRPS(L)                                    ! *** WQRPSR - Shoot carbon biomass loss due to respiration
        WQRPSL(L) = (1.-FRPSD)*RLRPS*WQRPS(L)                           ! *** WQRPSL - Shoot carbon biomass loss due to non-respiration processes
        WQRPDL(L) = RLRPD*WQRPD(L)                                      ! *** WQRPDL - Detritus carbon biomass loss rate

        ! *** LOSSES TO ORGANIC PHOSPHOROUS 
        ! *** RPSPC - Phosphorus to Carbon ratio for shoots
        WQRPSRP(L) = RPSPC*WQRPSR(L)                                    ! *** WQRPSRP - Shoot phosphorus biomass loss due to respiration
        WQRPSLP(L) = RPSPC*WQRPSL(L)                                    ! *** WQRPSLP - Shoot phosphorus biomass loss due to non-respiration processes                                      
        WQRPDLP(L) = RPSPC*WQRPDL(L)                                    ! *** WQRPDLP - Detritus phosphorus biomass loss rate                                        

        ! *** LOSSES TO ORGANIC NITROGEN                                                                                         
        ! *** RPSNC - Nitrogen to Carbon ratio for shoots
        WQRPSRN(L) = RPSNC*WQRPSR(L)                                    ! *** WQRPSRN - Shoot nitrogen biomass loss due to respiration
        WQRPSLN(L) = RPSNC*WQRPSL(L)                                    ! *** WQRPSLN - Shoot nitrogen biomass loss due to non-respiration processes                                        
      ENDDO
                                                                                      
      ! *** UPDATE EPIPHYTES
      IF( IRPEME > 0 )THEN
        DO LP = 1,NLRPEM(ND)
          L = LLRPEM(LP,ND)
          
          ! *** LOSSES TO ORGANIC CARBON                                                                                           
          WQRPER(L) = RRPE(L)*WQRPE(L)                                                                                         
          WQRPEL(L) = RLRPE*WQRPE(L)                                                                                           

          ! *** LOSSES TO ORGANIC PHOSPHOROUS                                                                                      
          WQRPERP(L) = RPEPC*WQRPER(L)                                                                                         
          WQRPELP(L) = RPEPC*WQRPEL(L)                                                                                         

          ! *** LOSSES TO ORGANIC NITROGEN                                                                                         
          WQRPERN(L) = RPENC*WQRPER(L)                                                                                         
          WQRPELN(L) = RPENC*WQRPEL(L)                                                                                         
        ENDDO                                                                                                                  
      ENDIF
     
      ! *** UPDATE WATER COLUMN ORGANIC CARBON, PHOSPHOROUS AND NITROGEN                                                       
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))                                                                                        
     
        ! *** WATER COLUMN REFRACTORY PARTICULATE ORGANIC CARBON                                                                 
        WQV(L,K,IROC) = WQV(L,K,IROC) + DTDHWQ*( FCRRPS*WQRPSR(L) + FCRLRPS*WQRPSL(L)                &      ! *** RPOC  (FCRRPS, FCRLRPS, FCRRPE are carbon fractions loss to RPOC)
                                      + FCRRPE*WQRPER(L) + FCRLRPE*WQRPEL(L) + FCRLRPD*WQRPDL(L) )
     
        ! *** WATER COLUMN LABILE PARTICULATE ORGANIC CARBON                                                                     
        WQV(L,K,ILOC) = WQV(L,K,ILOC) + DTDHWQ*( FCLRPS*WQRPSR(L) + FCLLRPS*WQRPSL(L)                &      ! *** LPOC  (FCLRPS, FCLLRPS, FCLRPE are carbon fractions loss to LPOC)
                                      + FCLRPE*WQRPER(L) + FCLLRPE*WQRPEL(L) + FCLLRPD*WQRPDL(L) )
     
        ! *** WATER COLUMN DISSOLVED ORGANIC CARBON                                                                              
        WQV(L,K,IDOC) = WQV(L,K,IDOC) + DTDHWQ*( FCDRPS*WQRPSR(L) + FCDLRPS*WQRPSL(L)                   &   ! *** DOC  (FCDRPS, FCDLRPS, FCDRPE are carbon fractions loss to DOC)
                                      + FCDRPE*WQRPER(L) + FCDLRPE*WQRPEL(L) + FCDLRPD*WQRPDL(L) )
     
        ! *** WATER COLUMN REFRACTORY PARTICULATE ORGANIC PHOSPHOROUS                                                            
        WQV(L,K,IROP) = WQV(L,K,IROP) + DTDHWQ*( FPRRPS*WQRPSRP(L)+FPRLRPS*WQRPSLP(L)                   &   ! *** RPOP  (FPRRPS, FPRLRPS, FPRRPE are carbon fractions loss to RPOP)
                                      + FPRRPE*WQRPERP(L) + FPRLRPE*WQRPELP(L) + FPRLRPD*WQRPDLP(L) )
     
      ! *** WATER COLUMN LABILE PARTICULATE ORGANIC PHOSPHOROUS                                                                
        WQV(L,K,ILOP) = WQV(L,K,ILOP) + DTDHWQ*( FPLRPS*WQRPSRP(L) + FPLLRPS*WQRPSLP(L)                 &   ! *** LPOP  (FPLRPS, FPLLRPS, FPLRPE are carbon fractions loss to LPOP)
                                      + FPLRPE*WQRPERP(L) + FPLLRPE*WQRPELP(L) + FPLLRPD*WQRPDLP(L) )
     
        ! *** WATER COLUMN DISSOLVED ORGANIC PHOSPHOROUS                                                                         
        WQV(L,K,IDOP) = WQV(L,K,IDOP) + DTDHWQ*( FPDRPS*WQRPSRP(L) + FPDLRPS*WQRPSLP(L)                 &   ! *** DOP   (FPDRPS, FPDLRPS, FPDRPE are carbon fractions loss to DOP)
                                      + FPDRPE*WQRPERP(L) + FPDLRPE*WQRPELP(L) + FPDLRPD*WQRPDLP(L) )
     
        ! *** WATER COLUMN TOTAL PHOSPHATE                                                                                       
        WQV(L,K,IP4D) = WQV(L,K,IP4D) + DTDHWQ*( FPIRPS*WQRPSRP(L) + FPILRPS*WQRPSLP(L)                 &   ! *** PO4T  (FPIRPS, FPILRPS, FPIRPE are carbon fractions loss to PO4T)  (INCOMPLETE)
                                      + FPIRPE*WQRPERP(L) + FPILRPE*WQRPELP(L) + FPILRPD*WQRPDLP(L) )
     
        ! *** WATER COLUMN REFRACTORY PARTICULATE ORGANIC NITROGEN                                                               
        WQV(L,K,IRON) = WQV(L,K,IRON) + DTDHWQ*( FNRRPS*WQRPSRN(L) + FNRLRPS*WQRPSLN(L)                 &   ! *** RPON   (FNRRPS, FNRLRPS, FNRRPE are carbon fractions loss to RPON)
                                      + FNRRPE*WQRPERN(L) + FNRLRPE*WQRPELN(L) + FNRLRPD*WQRPDLN(L) )
     
        ! *** WATER COLUMN LABILE PARTICULATE ORGANIC NITROGEN                                                                   
        WQV(L,K,ILON) = WQV(L,K,ILON) + DTDHWQ*( FNLRPS*WQRPSRN(L) + FNLLRPS*WQRPSLN(L)                 &   ! *** LPON   (FNLRPS, FNLLRPS, FNLRPE are carbon fractions loss to LPON)
                                      + FNLRPE*WQRPERN(L) + FNLLRPE*WQRPELN(L) + FNLLRPD*WQRPDLN(L) )
     
        ! *** WATER COLUMN DISSOLVED ORGANIC NITROGEN                                                                            
        WQV(L,K,IDON) = WQV(L,K,IDON) + DTDHWQ*( FNDRPS*WQRPSRN(L) + FNDLRPS*WQRPSLN(L)                 &   ! *** DON    (FNDRPS, FNDLRPS, FNDRPE are carbon fractions loss to DON)
                                      + FNDRPE*WQRPERN(L) + FNDLRPE*WQRPELN(L) + FNDLRPD*WQRPDLN(L) )
     
        ! *** WATER COLUMN AMMONIA NITROGEN                                                                                      
        WQV(L,K,INHX) = WQV(L,K,INHX) + DTDHWQ*( FNIRPS*WQRPSRN(L) + FNILRPS*WQRPSLN(L)                 &   ! *** NH4, EQ. INCOMPLETE
                                      + FNIRPE*WQRPERN(L) + FNILRPE*WQRPELN(L) + FNILRPD*WQRPDLN(L) )
     
      ENDDO                                                                                                                  
     
      ! **********************************************************************C                                                
      ! *** CALCULATE WATER COLUMN SHOOT AND EPIPHYTE UPTAKE OF                                                                
      ! *** DISSOLVED PHOSPHATE, AMMONIA, AND NO3                                                                              
     
      ! *** CALCULATE THE FRACTION OF DISSOLVED PHOSPHATE UPTAKE FROM WATER                                                    
      ! *** COLUMN (FRPSPW) BY PLANT SHOOTS                                                                                    
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TMPWAT = RKHPRPR*WQV(L,K,10)                                    ! ??, BUT (38) MEANS PO4DW, DISSOLVED ONLY, NO
        TMPBED = RKHPRPS*SM2PO4(L)                                                                                           
        FRPSPW(L) = TMPWAT/(TMPWAT+TMPBED+1E-18)                        ! EQ. (38)
      ENDDO                                                                                                                  
   
      ! *** UPDATE TOTAL PHOSPHATE IN WATER COLUMN                                                                             
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))                                                                                           
        WQV(L,K,IP4D) = WQV(L,K,IP4D) - DTDHWQ*( RPEPC*PRPE(L)*WQRPE(L) + RPSPC*FRPSPW(L)*PRPS(L)*WQRPS(L) )  ! PO4T, EQ. COMPLETED
        WQV(L,K,IP4D) = MAX(WQV(L,K,IP4D),0.)
      ENDDO                                                                                                                  
     
      ! *** CALCULATE THE FRACTION OF AMMONIA AND NO3 UPTAKE FROM WATER                                                        
      ! *** COLUMN (FRPSNW) BY PLANT SHOOTS                                                                                    
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TMPWAT = RKHNRPR*(WQV(L,K,INHX) + WQV(L,K,INOX))
        TMPBED = RKHNRPS*(SM2NH4(L)+SM2NO3(L))                                                                               
        FRPSNW(L) = TMPWAT/(TMPWAT+TMPBED+1E-18)                              ! EQ. (45)
      ENDDO                                                                                                                  
     
      ! *** CALCULATE THE AMMONIA PREFERENCE FOR SHOOTS                                                                        
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TOP1 = (WQV(L,K,INHX) + SM2NH4(L))*(WQV(L,K,INOX) + SM2NO3(L))                ! WHY USE WATER COLUMN AND BED, BOTH OF THEM?
        TOP2 = RKHNPRPS*(WQV(L,K,INHX) + SM2NH4(L))                                                                              
        BOT1 = RKHNPRPS + (WQV(L,K,INHX) + SM2NH4(L))                                                                              
        BOT2 = RKHNPRPS + (WQV(L,K,INOX) + SM2NO3(L))                                                                              
        BOT3 = WQV(L,K,INHX) + SM2NH4(L) + (WQV(L,K,INOX) + SM2NO3(L))     
        BOT1 = MAX(BOT1,1E-12)                                                            
        BOT2 = MAX(BOT2,1E-12)                                                            
        BOT3 = MAX(BOT3,1E-12)                                                            
        PNRPS(L) = TOP1/(BOT1*BOT2) + TOP2/(BOT2*BOT3)                           ! EQ. (44A)
      ENDDO                                                                                                                  
   
      ! *** CALCULATE THE AMMONIA PREFERENCE FOR EPIPHYTES                                                                     
      IF( IRPEME > 0 )THEN
        DO LP = 1,NLRPEM(ND)
          L = LLRPEM(LP,ND)
          K = KSZ(L)
          ! ***    NH3        NO3                                                                                                 
          TOP1 = WQV(L,K,INHX)*WQV(L,K,INOX)                                                                                       
          TOP2 = RKHNPRPE*WQV(L,K,INHX)                                                                                          
          BOT1 = RKHNPRPE+WQV(L,K,INHX)                                                                                          
          BOT2 = RKHNPRPE+WQV(L,K,INOX)                                                                                          
          BOT3 = WQV(L,K,INHX)+WQV(L,K,INOX)     
          BOT1 = MAX(BOT1,1E-12)                                                               
          BOT2 = MAX(BOT2,1E-12)                                                               
          BOT3 = MAX(BOT3,1E-12)                                                               
          PNRPE(L) = TOP1/(BOT1*BOT2)+TOP2/(BOT2*BOT3)                    ! EQ. (44B)                                                                    
        ENDDO                                                                                                                  
      ENDIF
     
      ! *** UPDATE WATER COLUMN AMMONIA AND NO3 NITROGEN                                                                       
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))                                                                                           
        TMPNH4S = PNRPS(L)                                                                                                   
        TMPNO3S = 1.-PNRPS(L)                                                                                                
        TMPNH4E = PNRPE(L)                                                                                                   
        TMPNO3E = 1.-PNRPE(L)                                                                                                
        WQV(L,K,INHX) = WQV(L,K,INHX) - DTDHWQ*(RPENC*TMPNH4E*PRPE(L)*WQRPE(L) + RPSNC*TMPNH4S*FRPSNW(L)*PRPS(L)*WQRPS(L))  ! NH4, EQ. COMPLETED
        WQV(L,K,INOX) = WQV(L,K,INOX) - DTDHWQ*(RPENC*TMPNO3E*PRPE(L)*WQRPE(L) + RPSNC*TMPNO3S*FRPSNW(L)*PRPS(L)*WQRPS(L))  ! NO3, EQ. COMPLETED
        WQV(L,K,INHX) = MAX(WQV(L,K,INHX),0.)
        WQV(L,K,INOX) = MAX(WQV(L,K,INOX),0.)
      ENDDO                                                                                                                  
     
      ! **********************************************************************C                                                
      ! *** CALCULATE WATER COLUMN SOURCES OF DISSOLVED OXYGEN DUE                                                             
      ! *** TO SHOOT AND EPIPHYTE GROWTH                                                                                       
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))       
        ! *** RPSOC - Plant shoot oxygen to carbon ratio (gO2/gC)
        ! *** PRPS  - Net shoot production rate (/day)
        ! *** WQRPS - Shoot mass (gC/m2)
        ! *** RPEOC - Epiphyte oxygen to carbon ratio (gO2/gC)
        ! *** PRPS  - Net shoot production rate (/day)
        ! *** WQRPS - Shoot mass (gC/m2)
        ! *** WQV(L,K,19) - Dissolved  Oxygen (g/m3)
        WQV(L,K,IDOX) = WQV(L,K,IDOX) + DTDHWQ*( RPSOC*WQRPS(L)*(PRPS(L)-RRPS(L)) + RPEOC*WQRPE(L)*(PRPE(L)-RRPE(L)) )  ! *** D.O.
        WQV(L,K,IDOX) = MAX(WQV(L,K,IDOX),0.0)
      ENDDO                                                                                                                  
    ENDIF    ! *** END OF BLOCK FOR SKIPPING WATER COLUMN LINKAGE
   
    !901   CONTINUE                                                                                                          
   
    ! **********************************************************************C                                                
    ! *** CALCULATE SEDIMENT BED ROOT RESPIRATION AND NON-RESPIRATION                                                        
    ! *** LOSSES TO SEDIMENT DIAGENESIS ORGANIC MATTER STATE VARIABLES                                                       
    ! *** AND TO PHOSPHATE AND AMMONIA                                                                                       
   
    !     GO TO 902   ! SKIP COUPLING WITH THE SEDIMENT DIAGENESIS MODEL.  FULL DIAGENESIS MUST BE ON TO COUPLE
    
    IF( IRPEMSD == 0 .AND. IWQBEN == 1 )THEN   
      ! *** LOSSES TO ORGANIC CARBON                                                                                           
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRR(L) = RRPR(L)*WQRPR(L)                                    ! FROM EQ. (32B)
        WQRPRL(L) = RLRPR*WQRPR(L)                                                                                           
      ENDDO                                                                                                                  
     
      ! *** LOSSES TO ORGANIC PHOSPHOROUS                                                                                      
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRRP(L) = RPRPC*WQRPRR(L)                                    !  ! FROM EQ. (34B)
        WQRPRLP(L) = RPRPC*WQRPRL(L)                                                                                         
      ENDDO                                                                                                                  
     
      ! *** LOSSES TO ORGANIC NITROGEN                                                                                         
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRRN(L) = RPRNC*WQRPRR(L)                                    ! FROM EQ. (39B)
        WQRPRLN(L) = RPRNC*WQRPRL(L)                                                                                         
      ENDDO                                                                                                                  
     
      ! *** UPDATE SEDIMENT BED ORGANIC CARBON, PHOSPHOROUS AND NITROGEN                                                       
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        DTDHWQ = DTWQ/SMHSED(1)
     
        ! *** SEDIMENT BED REFRACTORY PARTICULATE ORGANIC CARBON FRPRRPG1  FRPRRPG2   FRPRRPG3                                    
     
        SMPOC(L,1) = SMPOC(L,1) + DTDHWQ*FRPRRPG1*(FCRRPR*WQRPRR(L)+FCRLRPR*WQRPRL(L) )  ! EQ. (32B), REFRACTORY
        SMPOC(L,2) = SMPOC(L,2) + DTDHWQ*FRPRRPG2*(FCRRPR*WQRPRR(L)+FCRLRPR*WQRPRL(L) )
        SMPOC(L,3) = SMPOC(L,3) + DTDHWQ*FRPRRPG3*(FCRRPR*WQRPRR(L)+FCRLRPR*WQRPRL(L) )
     
        ! *** SEDIMENT BED LABILE PARTICULATE ORGANIC CARBON                                                                     
        SMPOC(L,1) = SMPOC(L,1) + DTDHWQ*( FCLRPR*WQRPRR(L)+FCLLRPR*WQRPRL(L) )           ! EQ. (32C), LABILE
     
        ! *** SEDIMENT BED DISSOLVED ORGANIC CARBON                                                                              
        SMPOC(L,1) = SMPOC(L,1)+DTDHWQ*(FCDRPR*WQRPRR(L)+FCDLRPR*WQRPRL(L) )            ! EQ. (33B), DISSOLVED IN WATER COLUMN IS AS
     
        ! *** SEDIMENT BED REFRACTORY PARTICULATE ORGANIC PHOSPHOROUS                                                            
        SMPOP(L,1) = SMPOP(L,1) + DTDHWQ*FRPRRPG1*(FPRRPR*WQRPRRP(L)+FPRLRPR*WQRPRLP(L) )
        SMPOP(L,2) = SMPOP(L,2) + DTDHWQ*FRPRRPG2*(FPRRPR*WQRPRRP(L)+FPRLRPR*WQRPRLP(L) )
        SMPOP(L,3) = SMPOP(L,3) + DTDHWQ*FRPRRPG3*(FPRRPR*WQRPRRP(L)+FPRLRPR*WQRPRLP(L) )
     
        ! *** SEDIMENT BED LABILE PARTICULATE ORGANIC PHOSPHOROUS                                                                
        SMPOP(L,1) = SMPOP(L,1) + DTDHWQ*(FPLRPR*WQRPRRP(L)+FPLLRPR*WQRPRLP(L) )
     
        ! *** SEDIMENT BED DISSOLVED ORGANIC PHOSPHOROUS                                                                         
        SMPOP(L,1) = SMPOP(L,1) + DTDHWQ*(FPDRPR*WQRPRRP(L)+FPDLRPR*WQRPRLP(L) )
     
        ! *** SEDIMENT BED TOTAL PHOSPHATE                                                                                       
        SM2PO4(L) = SM2PO4(L) + DTDHWQ*(FPIRPR*WQRPRRP(L)+FPILRPR*WQRPRLP(L) )             ! EQ. (37B), INCOMPLETE
     
        ! *** SEDIMENT BED REFRACTORY PARTICULATE ORGANIC NITROGEN                                                               
        SMPON(L,1) = SMPON(L,1) + DTDHWQ*FRPRRPG1*(FNRRPR*WQRPRRN(L)+FNRLRPR*WQRPRLN(L) )
        SMPON(L,2) = SMPON(L,2) + DTDHWQ*FRPRRPG2*(FNRRPR*WQRPRRN(L)+FNRLRPR*WQRPRLN(L) )
        SMPON(L,3) = SMPON(L,3) + DTDHWQ*FRPRRPG3*(FNRRPR*WQRPRRN(L)+FNRLRPR*WQRPRLN(L) )
     
        ! *** SEDIMENT BED LABILE PARTICULATE ORGANIC NITROGEN                                                                   
        SMPON(L,1) = SMPON(L,1)+DTDHWQ*(FNLRPR*WQRPRRN(L)+FNLLRPR*WQRPRLN(L) )
     
        ! *** SEDIMENT BED DISSOLVED ORGANIC NITROGEN                                                                            
        SMPON(L,1) = SMPON(L,1)+DTDHWQ*(FNDRPR*WQRPRRN(L)+FNDLRPR*WQRPRLN(L) )
     
        ! *** SEDIMENT BED AMMONIA NITROGEN                                                                                      
        SM2NH4(L) = SM2NH4(L) + DTDHWQ*(FNIRPR*WQRPRRN(L)+FNILRPR*WQRPRLN(L) )              ! EQ. (42B), INCOMPLETE
      ENDDO                                                                                                                  
     
      ! **********************************************************************C                                                
      ! *** CALCULATE SEDIMENT BED ROOT UPTAKE OF DISSOLVED PHOSPHATE,                                                         
      ! *** AMMONIA AND NO3                                                                                                    
     
      ! *** UPDATE TOTAL PHOSPHATE IN SEDIMENT BED                                                                             
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))                                   
        SM2PO4(L) = SM2PO4(L)-DTDHWQ*RPRPC*(1.-FRPSPW(L))*PRPS(L)*WQRPS(L)    ! EQ. (37B), COMPLETED
        SM2PO4(L) = MAX(SM2PO4(L), 1.E-12)
      ENDDO                                                                                                                  
     
      ! *** UPDATE SEDIMENT BED AMMONIA AND NO3 NITROGEN                                                                       
      DO LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ    = DTWQ/(DZC(L,K)*HWQ(L))
        TMPNH4S   = PNRPS(L)
        TMPNO3S   = 1. - PNRPS(L)                                                                                                
        TMPBED    = 1. - FRPSNW(L)                                                                                                
        SM2NH4(L) = SM2NH4(L) - DTDHWQ*RPSNC*TMPNH4S*TMPBED*PRPS(L)*WQRPS(L)   ! EQ. (42B), COMPLETED
        SM2NO3(L) = SM2NO3(L) - DTDHWQ*RPSNC*TMPNO3S*TMPBED*PRPS(L)*WQRPS(L)   ! EQ. (43B)
        SM2NH4(L) = MAX(SM2NH4(L), 1.E-12)
        SM2NO3(L) = MAX(SM2NO3(L), 1.E-12)
      ENDDO                                                                                                                  
    ENDIF   ! *** END OF SEDIMENT BED LINKAGE BYPASS
  
    902   CONTINUE                                                                                                          

  ENDDO    ! *** END OF DOMAIN LOOP 
  !$OMP END DO
   
  !$OMP END PARALLEL

  ! **********************************************************************C 
  ! *** RESTORE OPEN BOUNDARY CONCENTRATIONS                                               
  IF( IRPEMWC == 0 )THEN
    DO IOBC = 1,NBCSOP  
      L = LOBCS(IOBC)
      K = KSZ(L)
      WQV(L,K,1:NWQV) = WQBCV(IOBC,1:NWQV)
    ENDDO  
  ENDIF
  
   
  ! **********************************************************************C                                                
  JSRPEM = 0                                                                                                             

  RETURN

END SUBROUTINE                                                                                                                    
                                                                                                                      
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine RPEMINP_JNP
!
!> @details  READING ALL INPUT FOR ROOTED PLANTS AND                                                        
!>           EPIPHYTES GROWING ON ROOTED PLANTS AND INITIALIZES VARIABLES                                                      
!>           INCLUDING TEMPERATURE GROWTH RELATIONSHIPS
!---------------------------------------------------------------------------!
! THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION  
! LAST MODIFIED BY JOHN HAMRICK ON 23 JUNE 2004
!---------------------------------------------------------------------------!                                            

SUBROUTINE RPEMINP_JNP

  Use fson
  Use fson_value_m, Only: fson_value_count, fson_value_get
  Use INFOMOD,ONLY:SKIPCOM,READSTR

  Character(len = 79), allocatable :: TITLE(:)
  Character*80 STR*200
  Type(fson_value), Pointer :: json_data, item
  
  Integer :: L,LDATA,LL,I,J,NT,LG
  Real    :: SMNH4,SMNO3,SMPO4,TP2RPE,RPST,RPRT,RPET,RPDT,TIMETMP,WTEMP 
  
  Real,allocatable :: FRPRRPG_TEMP(:)
  
  if( process_id == master_id )then
    json_data => fson_parse("wq_rpem.jnp")
    Write(*,'(A)') ' WQ: WQRPEM.JNP - RPEM CONTROL FILE'
         
    !Call fson_get(json_data, "title", TITLE)
    Call fson_get(json_data, "initial_condition_option", INITRPEM)
    Call fson_get(json_data, "water_column_nutrient_linkage", IRPEMWC)
    Call fson_get(json_data, "sediment_diagenesis_nutrient_linkage", IRPEMSD)
    Call fson_get(json_data, "initial_carbon_biomass.shoot", RPSO)
    Call fson_get(json_data, "initial_carbon_biomass.root", RPRO)
    Call fson_get(json_data, "initial_carbon_biomass.epiphyte", RPEO)
    Call fson_get(json_data, "initial_carbon_biomass.detritus", RPDO)
    
    Call fson_get(json_data, "const_bed_porewater_concentration.ammonia", SMNH4)
    Call fson_get(json_data, "const_bed_porewater_concentration.nitrate", SMNO3)
    Call fson_get(json_data, "const_bed_porewater_concentration.phosphate", SMPO4)
    
    Call fson_get(json_data, "rooted_plant.average_shoot_height_above_bed", HRPS)
    Call fson_get(json_data, "rooted_plant.stoichiometric_ratios.root.nitrogen_to_carbon", RPRNC)
    Call fson_get(json_data, "rooted_plant.stoichiometric_ratios.root.phosphorus_to_carbon", RPRPC)
    Call fson_get(json_data, "rooted_plant.stoichiometric_ratios.shoot.oxygen_to_carbon", RPSOC)
    Call fson_get(json_data, "rooted_plant.stoichiometric_ratios.shoot.nitrogen_to_carbon", RPSNC)
    Call fson_get(json_data, "rooted_plant.stoichiometric_ratios.shoot.phosphorus_to_carbon", RPSPC)
    
    
    Call fson_get(json_data, "rooted_plant.respiration.root.max_respiration_rate", RMRPR)
    Call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.suboptimal_temperature", TR1RPR)
    Call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.superoptimal_temperature", TR2RPR)
    Call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.suboptimal_temperature_effect", rKTR1RPR)
    Call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.superoptimal_temperature_effect", rKTR2RPR)
    
    
    Call fson_get(json_data, "rooted_plant.respiration.shoot.max_respiration_rate", RMRPS)
    Call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.suboptimal_temperature", TR1RPS)
    Call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.superoptimal_temperature", TR2RPS)
    Call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.suboptimal_temperature_effect", rKTR1RPS)
    Call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.superoptimal_temperature_effect", rKTR2RPS)
    
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_carbon.RPOC", FCRrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_carbon.LPOC", FCLrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_carbon.DOC", FCDrpr)
    
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.RPOP", FPRrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.LPOP", FPLrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.DOP", FPDrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.PO4", FPIrpr)
    
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.RPON", FNRrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.LPON", FNLrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.DON", FNDrpr)
    Call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.NH4", FNIrpr)
    
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_carbon.RPOC", FCRrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_carbon.LPOC", FCLrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_carbon.DOC", FCDrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.RPOP", FPRrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.LPOP", FPLrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.DOP", FPDrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.PO4", FPIrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.RPON", FNRrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.LPON", FNLrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.DON", FNDrps)
    Call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.NH4", FNIrps)
    
    Call fson_get(json_data, "rooted_plant.growth.max_rooted_depth", HRPEMIC)
    Call fson_get(json_data, "rooted_plant.growth.optimum_water_depth_for_growth", HOPT)
    Call fson_get(json_data, "rooted_plant.growth.max_SSW_for_optimum_growth", rISSOM)
    Call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.transport_option", iJRPRS)
    Call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.const_positive_carbon_transport", rJRPRSC)
    Call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.transfer_rate_ratio", rKRPORS)
    Call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.observed_ratio", ROSR)
    Call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.light_dependent_transfer_rate", rKRPRS)
    Call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.half_sat_solar_ratio_at_surface", rISSS)
    
    Call fson_get(json_data, "rooted_plant.growth.root.nutrient_limits.half_sat_const_for_nitrogen_uptake_from_bed", rKHNRPR)
    Call fson_get(json_data, "rooted_plant.growth.root.nutrient_limits.half_sat_const_for_phosphorus_uptake_from_bed", rKHPRPR)
    
    Call fson_get(json_data, "rooted_plant.growth.shoot.max_growth_rate", PMRPS)
    Call fson_get(json_data, "rooted_plant.growth.shoot.fraction_of_production_transferred_to_roots", FPRPR)
    Call fson_get(json_data, "rooted_plant.growth.shoot.half_sat_nitrogen_preference", rKHNPRPS)
    Call fson_get(json_data, "rooted_plant.growth.shoot.self_shading", rKSH)
    Call fson_get(json_data, "rooted_plant.growth.shoot.half_sat_const_for_irradiance", rKHI)
    Call fson_get(json_data, "rooted_plant.growth.shoot.nutrient_limits.half_sat_const_for_nitrogen_uptake_from_water", rKHNRPS)
    Call fson_get(json_data, "rooted_plant.growth.shoot.nutrient_limits.half_sat_const_for_phosphorus_uptake_from_water", rKHPRPS)
    
    
    Call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.suboptimal_temperature", TP1RPS)
    Call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.superoptimal_temperature", TP2RPS)
    Call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.suboptimal_temperature_effect", rKTP1RPS)
    Call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.superoptimal_temperature_effect", rKTP2RPS)
    
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.loss_rate", rLRPR)
    
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_carbon.RPOC", FCRLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_carbon.LPOC", FCLLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_carbon.DOC", FCDLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.RPOP", FPRLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.LPOP", FPLLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.DOP", FPDLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.PO4", FPILrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.RPON", FNRLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.LPON", FNLLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.DON", FNDLrpr)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.NH4", FNILrpr)
    
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.loss_rate", rLRPS)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_to_detritus", FRPSD)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_carbon.RPOC", FCRLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_carbon.LPOC", FCLLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_carbon.DOC", FCDLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.RPOP", FPRLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.LPOP", FPLLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.DOP", FPDLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.PO4", FPILrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.RPON", FNRLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.LPON", FNLLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.DON", FNDLrps)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.NH4", FNILrps)
    
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.loss_rate_at_bottom", rLRPD)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_carbon.RPOC", FCRLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_carbon.LPOC", FCLLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_carbon.DOC", FCDLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.RPOP", FPRLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.LPOP", FPLLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.DOP", FPDLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.PO4", FPILrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.RPON", FNRLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.LPON", FNLLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.DON", FNDLrpd)
    Call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.NH4", FNILrpd)
    
    Call fson_get(json_data, "rooted_plant.transfer_of_root_generated_RPOM_to_sediment", FRPRRPG_TEMP)
    FRPRRPG1 = FRPRRPG_TEMP(1)
    FRPRRPG2 = FRPRRPG_TEMP(2)
    FRPRRPG3 = FRPRRPG_TEMP(3)
    
    
    Call fson_get(json_data, "epiphyte.stoichiometric_ratios.oxygen_to_carbon", RPEOC)
    Call fson_get(json_data, "epiphyte.stoichiometric_ratios.nitrogen_to_carbon", RPENC)
    Call fson_get(json_data, "epiphyte.stoichiometric_ratios.phosphorus_to_carbon", RPEPC)
    Call fson_get(json_data, "epiphyte.stoichiometric_ratios.carbon_to_chlorophyll_ratio", CChlRPE)
    
    Call fson_get(json_data, "epiphyte.respiration.max_respiration_rate", RMRPE)
    Call fson_get(json_data, "epiphyte.respiration.temperature_effects.suboptimal_temperature", TR1RPE)
    Call fson_get(json_data, "epiphyte.respiration.temperature_effects.superoptimal_temperature", TR2RPE)
    Call fson_get(json_data, "epiphyte.respiration.temperature_effects.suboptimal_temperature_effect", rKTR1RPE)
    Call fson_get(json_data, "epiphyte.respiration.temperature_effects.superoptimal_temperature_effect", rKTR2RPE)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_carbon.RPOC", FCRrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_carbon.LPOC", FCLrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_carbon.DOC", FCDrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.RPOP", FPRrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.LPOP", FPLrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.DOP", FPDrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.PO4", FPIrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.RPON", FNRrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.LPON", FNLrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.DON", FNDrpe)
    Call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.NH4", FNIrpe)
    
    Call fson_get(json_data, "epiphyte.growth.growing_option", IRPEME)
    Call fson_get(json_data, "epiphyte.growth.max_growth_rate", PMRPE)
    Call fson_get(json_data, "epiphyte.growth.half_sat_nitrogen_preference", rKHNPRPE)
    Call fson_get(json_data, "epiphyte.growth.max_SSW_for_optimum_growth", rISSOEM)
    Call fson_get(json_data, "epiphyte.growth.light_extinction_coef_for_chlorophyll", rKeRPE)
    Call fson_get(json_data, "epiphyte.growth.nutrient_limits.half_sat_const_for_nitrogen_uptake", rKHNRPE)
    Call fson_get(json_data, "epiphyte.growth.nutrient_limits.half_sat_const_for_phosphorus_uptake", rKHPRPE)
    
    Call fson_get(json_data, "epiphyte.growth.temperature_effects.suboptimal_temperature", TP1RPE)
    Call fson_get(json_data, "epiphyte.growth.temperature_effects.superoptimal_temperature", TP2RPE)
    Call fson_get(json_data, "epiphyte.growth.temperature_effects.suboptimal_temperature_effect", rKTP1RPE)
    Call fson_get(json_data, "epiphyte.growth.temperature_effects.superoptimal_temperature_effect", rKTP2RPE)
    
    Call fson_get(json_data, "epiphyte.non_respiration_loss.loss_rate", rLRPE)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_carbon.RPOC", FCRLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_carbon.LPOC", FCLLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_carbon.DOC", FCDLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.RPOP", FPRLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.LPOP", FPLLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.DOP", FPDLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.PO4", FPILrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.RPON", FNRLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.LPON", FNLLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.DON", FNDLrpe)
    Call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.NH4", FNILrpe)
    
    Call fson_get(json_data, "salinity_toxicity.salinity_effect_option", iSTOXRPE)
    Call fson_get(json_data, "salinity_toxicity.STOXS", STOXS)
    Call fson_get(json_data, "salinity_toxicity.STOXE", STOXE)
    
    Call fson_get(json_data, "output.number_of_EE_WC_writes", NRPEMEE)
    Call fson_get(json_data, "output.spatial_output.output_option", ISRPEMSPAC)
    Call fson_get(json_data, "output.spatial_output.output_frequency", ISRPEMSPFR)
    
    Call fson_get(json_data, "output.time_series_output.output_option", ISRPEMTIME)
    Call fson_get(json_data, "output.time_series_output.output_frequency", ISRPEMTIFR)
    Call fson_get(json_data, "output.time_series_output.number_of_locations", ISRPEMTILC)
  
  endif   ! *** End of master_id
  
  ! *** Scalar Variables
  Call Broadcast_Scalar(INITRPEM, master_id)
  Call Broadcast_Scalar(IRPEMWC,  master_id)
  Call Broadcast_Scalar(IRPEMSD,  master_id)
  Call Broadcast_Scalar(HRPEMIC,  master_id)
  Call Broadcast_Scalar(RPSO,     master_id)
  Call Broadcast_Scalar(RPRO,     master_id)
  Call Broadcast_Scalar(RPEO,     master_id)
  Call Broadcast_Scalar(RPDO,     master_id)
  Call Broadcast_Scalar(SMNH4,    master_id)
  Call Broadcast_Scalar(SMNO3,    master_id)
  Call Broadcast_Scalar(SMPO4,    master_id)
  
  Call Broadcast_Scalar(RPSOC,    master_id)
  Call Broadcast_Scalar(RPEOC,    master_id)
  Call Broadcast_Scalar(RPSNC,    master_id)
  Call Broadcast_Scalar(RPRNC,    master_id)
  Call Broadcast_Scalar(RPENC,    master_id)
  Call Broadcast_Scalar(RPSPC,    master_id)
  Call Broadcast_Scalar(RPRPC,    master_id)
  Call Broadcast_Scalar(RPEPC,    master_id)
  
  Call Broadcast_Scalar(PMRPS,    master_id)
  Call Broadcast_Scalar(FPRPR,    master_id)
  Call Broadcast_Scalar(RMRPS,    master_id)
  Call Broadcast_Scalar(rLRPS,    master_id)
  Call Broadcast_Scalar(FRPSD,    master_id)
  Call Broadcast_Scalar(rKHNPRPS, master_id)
  Call Broadcast_Scalar(rKSH,     master_id)
  Call Broadcast_Scalar(rKHI,     master_id)
  
  Call Broadcast_Scalar(RMRPR,    master_id)
  Call Broadcast_Scalar(rLRPR,    master_id)
  Call Broadcast_Scalar(rLRPD,    master_id)
  
  Call Broadcast_Scalar(iJRPRS,   master_id)
  Call Broadcast_Scalar(rJRPRSC,  master_id)
  Call Broadcast_Scalar(rKRPORS,  master_id)
  Call Broadcast_Scalar(ROSR,     master_id)
  Call Broadcast_Scalar(rKRPRS,   master_id)
  Call Broadcast_Scalar(rISSS,    master_id)
  
  Call Broadcast_Scalar(IRPEME,   master_id)
  Call Broadcast_Scalar(PMRPE,    master_id)
  Call Broadcast_Scalar(RMRPE,    master_id)
  Call Broadcast_Scalar(rLRPE,    master_id)
  Call Broadcast_Scalar(rKHNPRPE, master_id)
  
  Call Broadcast_Scalar(rKHNRPS,  master_id)
  Call Broadcast_Scalar(rKHNRPR,  master_id)
  Call Broadcast_Scalar(rKHNRPE,  master_id)
  Call Broadcast_Scalar(rKHPRPS,  master_id)
  Call Broadcast_Scalar(rKHPRPR,  master_id)
  Call Broadcast_Scalar(rKHPRPE,  master_id)
  
  Call Broadcast_Scalar(TP1RPS,   master_id)
  Call Broadcast_Scalar(TP2RPS,   master_id)
  Call Broadcast_Scalar(rKTP1RPS, master_id)
  Call Broadcast_Scalar(rKTP2RPS, master_id)
  Call Broadcast_Scalar(TP1RPE,   master_id)
  Call Broadcast_Scalar(TP2RPE,   master_id)
  Call Broadcast_Scalar(rKTP1RPE, master_id)
  Call Broadcast_Scalar(rKTP2RPE, master_id)
  
  Call Broadcast_Scalar(HRPS,     master_id)
  Call Broadcast_Scalar(HOPT,     master_id)
  Call Broadcast_Scalar(rKeRPE,   master_id)
  Call Broadcast_Scalar(CChlRPE,  master_id)
  Call Broadcast_Scalar(rISSOM,   master_id)
  Call Broadcast_Scalar(rISSOEM,  master_id)
  
  Call Broadcast_Scalar(iSTOXRPE, master_id)
  Call Broadcast_Scalar(STOXS,    master_id)
  Call Broadcast_Scalar(STOXE,    master_id)
  
  Call Broadcast_Scalar(TR1RPS,   master_id)
  Call Broadcast_Scalar(TR2RPS,   master_id)
  Call Broadcast_Scalar(rKTR1RPS, master_id)
  Call Broadcast_Scalar(rKTR2RPS, master_id)
  Call Broadcast_Scalar(TR1RPR,   master_id)
  Call Broadcast_Scalar(TR2RPR,   master_id)
  Call Broadcast_Scalar(rKTR1RPR, master_id)
  Call Broadcast_Scalar(rKTR2RPR, master_id)
  
  Call Broadcast_Scalar(TR1RPE,   master_id)
  Call Broadcast_Scalar(TR2RPE,  master_id)
  Call Broadcast_Scalar(rKTR1RPE,master_id)
  Call Broadcast_Scalar(rKTR2RPE,master_id)
  
  Call Broadcast_Scalar(FCRrps,   master_id)
  Call Broadcast_Scalar(FCLrps,   master_id)
  Call Broadcast_Scalar(FCDrps,   master_id)
  Call Broadcast_Scalar(FCRLrps,  master_id)
  Call Broadcast_Scalar(FCLLrps,  master_id)
  Call Broadcast_Scalar(FCDLrps,  master_id)
  
  Call Broadcast_Scalar(FCRrpr,   master_id)
  Call Broadcast_Scalar(FCLrpr,   master_id)
  Call Broadcast_Scalar(FCDrpr,   master_id)
  Call Broadcast_Scalar(FCRLrpr,  master_id)
  Call Broadcast_Scalar(FCLLrpr,  master_id)
  Call Broadcast_Scalar(FCDLrpr,  master_id)
  
  Call Broadcast_Scalar(FCRrpe,   master_id)
  Call Broadcast_Scalar(FCLrpe,   master_id)
  Call Broadcast_Scalar(FCDrpe,   master_id)
  Call Broadcast_Scalar(FCRLrpe,  master_id)
  Call Broadcast_Scalar(FCLLrpe,  master_id)
  Call Broadcast_Scalar(FCDLrpe,  master_id)
  
  Call Broadcast_Scalar(FCRLrpd,  master_id)
  Call Broadcast_Scalar(FCLLrpd,  master_id)
  Call Broadcast_Scalar(FCDLrpd,  master_id)
  
  Call Broadcast_Scalar(FPRrps,   master_id)
  Call Broadcast_Scalar(FPLrps,   master_id)
  Call Broadcast_Scalar(FPDrps,   master_id)
  Call Broadcast_Scalar(FPIrps,   master_id)
  Call Broadcast_Scalar(FPRLrps,  master_id)
  Call Broadcast_Scalar(FPLLrps,  master_id)
  Call Broadcast_Scalar(FPDLrps,  master_id)
  Call Broadcast_Scalar(FPILrps,  master_id)
  
  Call Broadcast_Scalar(FPRrpr,   master_id)
  Call Broadcast_Scalar(FPLrpr,   master_id)
  Call Broadcast_Scalar(FPDrpr,   master_id)
  Call Broadcast_Scalar(FPIrpr,   master_id)
  Call Broadcast_Scalar(FPRLrpr,  master_id)
  Call Broadcast_Scalar(FPLLrpr,  master_id)
  Call Broadcast_Scalar(FPDLrpr,  master_id)
  Call Broadcast_Scalar(FPILrpr,  master_id)
  
  Call Broadcast_Scalar(FPRrpe,   master_id)
  Call Broadcast_Scalar(FPLrpe,   master_id)
  Call Broadcast_Scalar(FPDrpe,   master_id)
  Call Broadcast_Scalar(FPIrpe,   master_id)
  Call Broadcast_Scalar(FPRLrpe,  master_id)
  Call Broadcast_Scalar(FPLLrpe,  master_id)
  Call Broadcast_Scalar(FPDLrpe,  master_id)
  Call Broadcast_Scalar(FPILrpe,  master_id)
  
  Call Broadcast_Scalar(FPRLrpd,  master_id)
  Call Broadcast_Scalar(FPLLrpd,  master_id)
  Call Broadcast_Scalar(FPDLrpd,  master_id)
  Call Broadcast_Scalar(FPILrpd,  master_id)
  
  Call Broadcast_Scalar(FNRrps,   master_id)
  Call Broadcast_Scalar(FNLrps,   master_id)
  Call Broadcast_Scalar(FNDrps,   master_id)
  Call Broadcast_Scalar(FNIrps,   master_id)
  Call Broadcast_Scalar(FNRLrps,  master_id)
  Call Broadcast_Scalar(FNLLrps,  master_id)
  Call Broadcast_Scalar(FNDLrps,  master_id)
  Call Broadcast_Scalar(FNILrps,  master_id)
  
  Call Broadcast_Scalar(FNRrpr,   master_id)
  Call Broadcast_Scalar(FNLrpr,   master_id)
  Call Broadcast_Scalar(FNDrpr,   master_id)
  Call Broadcast_Scalar(FNIrpr,   master_id)
  Call Broadcast_Scalar(FNRLrpr,  master_id)
  Call Broadcast_Scalar(FNLLrpr,  master_id)
  Call Broadcast_Scalar(FNDLrpr,  master_id)
  Call Broadcast_Scalar(FNILrpr,  master_id)
  
  Call Broadcast_Scalar(FNRrpe,   master_id)
  Call Broadcast_Scalar(FNLrpe,   master_id)
  Call Broadcast_Scalar(FNDrpe,   master_id)
  Call Broadcast_Scalar(FNIrpe,   master_id)
  Call Broadcast_Scalar(FNRLrpe,  master_id)
  Call Broadcast_Scalar(FNLLrpe,  master_id)
  Call Broadcast_Scalar(FNDLrpe,  master_id)
  Call Broadcast_Scalar(FNILrpe,  master_id)
  
  Call Broadcast_Scalar(FNRLrpd,  master_id)
  Call Broadcast_Scalar(FNLLrpd,  master_id)
  Call Broadcast_Scalar(FNDLrpd,  master_id)
  Call Broadcast_Scalar(FNILrpd,  master_id)
  
  Call Broadcast_Scalar(FRPRRPG1, master_id)
  Call Broadcast_Scalar(FRPRRPG2, master_id)
  Call Broadcast_Scalar(FRPRRPG3, master_id)
  
  Call Broadcast_Scalar(ISRPEMSPAC, master_id)
  Call Broadcast_Scalar(ISRPEMSPFR, master_id)
  Call Broadcast_Scalar(ISRPEMTIME, master_id)
  Call Broadcast_Scalar(ISRPEMTIFR, master_id)
  Call Broadcast_Scalar(ISRPEMTILC, master_id)
  Call Broadcast_Scalar(NRPEMEE,    master_id)
  
  ! *** *******************************************************************C                                                
  ! *** SET INITIAL CONDITIONS                                                                          

  DO L = 1,LC
    WQRPS(L) = 0.0   ! *** Initial condition for shoot mass    ( gC/m^2 )
    WQRPR(L) = 0.0   ! *** Initial condition for root mass     ( gC/m^2 )                                                                                                           
    WQRPE(L) = 0.0   ! *** Initial condition for epiphyte mass ( gC/m^2 )                                                                                                         
    WQRPD(L) = 0.0   ! *** Initial condition for detritus mass ( gC/m^2 )                                                                                                         
    LMASKRPEM(L) = .FALSE.                                                                                                 
  ENDDO                                                                                                                  
  LMASKRPEM_Global(:) = .false.
  
  IF( INITRPEM > 0 )THEN
    ! *** SET SPATIALLY VARIABLE INITIAL CONDITIONS                                                                          
    if( process_id == master_id )then
      IF( INITRPEM == 1 )THEN
        WRITE(*,'(A)')' WQ: WQRPEMSIC.INP'
        OPEN(1,FILE = 'wqrpemsic.inp')
        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      
        READ(1,*) LDATA
        
        DO LL = 1,LDATA                                                                                                        
          READ(1,*) I, J, RPST, RPRT, RPET, RPDT                                                                                   
          LG = LIJ_Global(I,J)                                                                                                         
          WQRPS_Global(LG) = RPST     ! *** WQRPS
          WQRPR_Global(LG) = RPRT     ! *** WQRPR                                                                                                      
          WQRPE_Global(LG) = RPET     ! *** WQRPE                                                                                                      
          WQRPD_Global(LG) = RPDT     ! *** WQRPD
          ! *** Only cell defined in wqrpemsic.inp file can have mask - DKT
          IF( HP_Global(LG) <= 10. ) LMASKRPEM_Global(LG) = .TRUE.   ! *** Only include the cell if not too deep
        ENDDO
        CLOSE(1)
        
      ELSEIF( INITRPEM == 2 )THEN
        ! *** READ RESTART CONDITIONS        
        WRITE(*,'(A)')' WQ: WQRPEMRST.INP'
        OPEN(1,FILE = 'wqrpemrst.inp')
        
        CALL SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES 
        READ(1,*)TIMETMP

    100 CONTINUE                                                                                                        
        READ(1,*,END = 200) LG, RPST, RPRT, RPET, RPDT                                                                             

        WQRPS_Global(LG) = RPST     ! *** WQRPS
        WQRPR_Global(LG) = RPRT     ! *** WQRPR                                                                                                      
        WQRPE_Global(LG) = RPET     ! *** WQRPE                                                                                                      
        WQRPD_Global(LG) = RPDT     ! *** WQRPD     
        LMASKRPEM_Global(LG) = .TRUE.   ! *** Only include the cell if not too deep
        GOTO 100     
        
    200 CONTINUE                                                                                                        
        CLOSE(1)
      ENDIF
    endif   ! *** End of master_id

    Call Broadcast_Array(WQRPS_Global, master_id)
    Call Broadcast_Array(WQRPR_Global, master_id)
    Call Broadcast_Array(WQRPE_Global, master_id)
    Call Broadcast_Array(WQRPD_Global, master_id)

    ! *** Map to Local Domain
    DO LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      IF( L > 1 )THEN
        WQRPS(L) = WQRPS_Global(LG)
        WQRPR(L) = WQRPR_Global(LG)
        WQRPE(L) = WQRPE_Global(LG)
        WQRPD(L) = WQRPD_Global(LG)
        LMASKRPEM(L) = LMASKRPEM_Global(LG)
      ENDIF
    ENDDO
  
  ELSE
    ! *** *******************************************************************C                                                
    ! *** SET SPATIALLY CONSTANT INITIAL CONDITIONS                                                                          
    DO L = 1,LC
      IF( HP(L) <= HRPEMIC )THEN  ! uniform values at depth <= HRPEMIC                                                    
        WQRPS(L) = RPSO   ! *** Initial condition for shoot mass    ( gC/m^2 )
        WQRPR(L) = RPRO   ! *** Initial condition for root mass     ( gC/m^2 )                                                 
        WQRPE(L) = RPEO   ! *** Initial condition for epiphyte mass ( gC/m^2 )
        WQRPD(L) = RPDO   ! *** Initial condition for detritus mass ( gC/m^2 )
        LMASKRPEM(L) = .TRUE.                                                                                              
      ENDIF                                                                                                              
    ENDDO
    
  ENDIF                                                                                                                  

  ! *** Count active RPEM cells
  NRPEM = 0
  DO L = 1,LC
    IF( LMASKRPEM(L) ) NRPEM = NRPEM + 1
  ENDDO

  
  ! *** ***********************************************************************************                                               
  ! *** GENERATE TEMPERATURE DEPENDENCY FOR GROWTH AND RESPIRATION OVER WQTDMIN TO WQTDMAX                                                        
  DO NT = 1,NWQTD  
    WTEMP = WQTDTEMP(NT)
  
    ! *** Shoot Production (Growth)
    RPEMTPrps(NT) = 1.
    IF( WTEMP < TP1RPS )THEN
      RPEMTPrps(NT) = EXP(-rKTP1RPS*(WTEMP-TP1RPS)*(WTEMP-TP1RPS) )
    ENDIF
    IF( WTEMP > TP2RPS )THEN
      RPEMTPrps(NT) = EXP(-rKTP2RPS*(WTEMP-TP2RPS)*(WTEMP-TP2RPS) )
    ENDIF
  
    ! *** Epiphyte Production (Growth)
    RPEMTPrpe(NT) = 1.
    IF( WTEMP < TP1RPE )THEN
      RPEMTPrpe(NT) = EXP(-rKTP1RPE*(WTEMP-TP1RPE)*(WTEMP-TP1RPE) )
    ENDIF
    IF( WTEMP > TP2RPE )THEN
      RPEMTPrpe(NT) = EXP(-rKTP2RPE*(WTEMP-TP2RPE)*(WTEMP-TP2RPE) )
    ENDIF
  
    ! *** Shoot Respiration
    RPEMTRrps(NT) = 1.
    IF( WTEMP < TR1RPS )THEN
      RPEMTRrps(NT) = EXP(-rKTR1RPS*(WTEMP-TR1RPS)*(WTEMP-TR1RPS) )
    ENDIF
    IF( WTEMP > TR2RPS )THEN
      RPEMTRrps(NT) = EXP(-rKTR2RPS*(WTEMP-TR2RPS)*(WTEMP-TR2RPS) )
    ENDIF
  
    ! *** Epiphyte Respiration
    RPEMTRrpe(NT) = 1.
    IF( WTEMP < TR1RPE )THEN
      RPEMTRrpe(NT) = EXP(-rKTR1RPE*(WTEMP-TR1RPE)*(WTEMP-TR1RPE) )
    ENDIF
    IF( WTEMP > TR2RPE )THEN
      RPEMTRrpe(NT) = EXP(-rKTR2RPE*(WTEMP-TR2RPE)*(WTEMP-TR2RPE) )
    ENDIF
  
    ! *** Root Respiration
    RPEMTRrpr(NT) = 1.
    IF( WTEMP < TR1RPR )THEN
      RPEMTRrpr(NT) = EXP(-rKTR1RPR*(WTEMP-TR1RPR)*(WTEMP-TR1RPR) )
    ENDIF
    IF( WTEMP > TR2RPR )THEN
      RPEMTRrpr(NT) = EXP(-rKTR2RPR*(WTEMP-TR2RPR)*(WTEMP-TR2RPR) )
    ENDIF
  
  ENDDO                                                                                                                  

  ! *** Write out the temperature dependencies for RPEM
  WRITE(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for RPEM Temperature Dependency",                               &
                                            "WQTDMIN = ", WQTDMIN, "WQTDMAX = ", WQTDMAX, "WQTDINC = ", WQTDINC, "NWQTD = ", NWQTD
  WRITE(2,'(/,A5,A10,A15,20A10)') "IT", "TEMP", "RPEMTPrps", "RPEMTPrpe", "RPEMTRrps", "RPEMTRrpe", "RPEMTRrpr"
                                            
  DO NT = 1,NWQTD
    WRITE(2,'(I5,F10.3,F15.3,20F10.5)') NT, WQTDTEMP(NT), RPEMTPrps(NT), RPEMTPrpe(NT), RPEMTRrps(NT), RPEMTRrpe(NT), RPEMTRrpr(NT)
  ENDDO
  
  ! *** ASSIGN BED POREWATER CONCENTRATIONS IF SEDIMENT DIAGENESIS NOT SIMULATED
  IF( IWQBEN /= 1 )THEN
    DO L = 2,LA
      SM2NH4(L) = SMNH4
      SM2NO3(L) = SMNO3
      SM2PO4(L) = SMPO4
    ENDDO
  ENDIF
  
  ! *** *******************************************************************C                                                
  RETURN
  
END SUBROUTINE RPEMINP_JNP
!---------------------------------------------------------------------------!
! EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
! Subroutine: Subroutine INIT_RPEMVARS
!
!> @details  Allocating and initializing the arrays 
!---------------------------------------------------------------------------!
!
!---------------------------------------------------------------------------!
SUBROUTINE INIT_RPEMVARS

  ! *** ALLOCATE
  ALLOCATE(WQRPS(LCM),WQRPR(LCM),WQRPE(LCM),WQRPD(LCM), &
     PRPS(LCM),RRPS(LCM),RRPR(LCM),PRPE(LCM),RRPE(LCM), &
     WQRPSR(LCM),WQRPSL(LCM),WQRPER(LCM),WQRPEL(LCM),WQRPDL(LCM), &
     WQRPSRP(LCM),WQRPSLP(LCM),WQRPERP(LCM),WQRPELP(LCM), &
     WQRPDLP(LCM),WQRPSRN(LCM),WQRPSLN(LCM),WQRPERN(LCM), &
     WQRPELN(LCM),WQRPDLN(LCM),WQRPRR(LCM),WQRPRL(LCM), &
     WQRPRRP(LCM),WQRPRLP(LCM),WQRPRRN(LCM),WQRPRLN(LCM), &
     FRPSPW(LCM),FRPSNW(LCM),PNRPS(LCM),PNRPE(LCM),RJRPRS(LCM), &
     XLIMTPRPS(LCM),XLIMTPRPE(LCM),XLIMTRRPS(LCM),XLIMTRRPE(LCM), &
     XLIMTRRPR(LCM),XLIMNRPS(LCM),XLIMNRPE(LCM),XLIMLRPS(LCM), &
     XLIMLRPE(LCM),RISS(LCM),RISSO(LCM),RISSOE(LCM), &
     RPEMTPRPS(NWQTDM),RPEMTPRPE(NWQTDM),RPEMTRRPS(NWQTDM), &
     RPEMTRRPE(NWQTDM),RPEMTRRPR(NWQTDM)  )
     
  ALLOCATE(IRPEMTS(LCM),JRPEMTS(LCM),LRPEMTS(LCM))
  ALLOCATE(LMASKRPEM(LCM))
  ! *** Allocating fluxes if sediment diagenesis is not simulated
  IF(  .NOT. ALLOCATED(SM2NH4)  )THEN
      ALLOCATE(SM2NH4(LCM))
      ALLOCATE(SM2NO3(LCM))
      ALLOCATE(SM2PO4(LCM))
      SM2NH4 = 0.0
      SM2NO3 = 0.0
      SM2PO4 = 0.0 
  ENDIF
  
  WQRPS =0
  WQRPR =0
  WQRPE =0
  WQRPD =0
  PRPS =0
  RRPS =0
  RRPR =0
  PRPE =0
  RRPE =0
  WQRPSR =0
  WQRPSL =0
  WQRPER =0
  WQRPEL =0
  WQRPDL =0
  WQRPSRP =0
  WQRPSLP =0
  WQRPERP =0
  WQRPELP =0
  WQRPDLP =0
  WQRPSRN =0
  WQRPSLN =0
  WQRPERN =0
  WQRPELN =0
  WQRPDLN =0
  WQRPRR =0
  WQRPRL =0
  WQRPRRP =0
  WQRPRLP =0
  WQRPRRN =0
  WQRPRLN =0
  FRPSPW =0
  FRPSNW =0
  PNRPS =0
  PNRPE =0
  RJRPRS =0
  XLIMTPRPS =0
  XLIMTPRPE =0
  XLIMTRRPS =0
  XLIMTRRPE =0
  XLIMTRRPR =0
  XLIMNRPS =0
  XLIMNRPE =0
  XLIMLRPS =0
  XLIMLRPE =0
  RISS =0
  RISSO =0
  RISSOE =0
  RPEMTPRPS = 0
  RPEMTPRPE = 0
  RPEMTRRPS = 0
  RPEMTRRPE = 0
  RPEMTRRPR = 0

  JSRPEM = 1
  
END SUBROUTINE

END MODULE
