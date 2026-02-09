! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
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

  use Variables_WQ
  use WQ_DIAGENESIS
  
implicit none

  ! *** INTEGER ARRAYS                                                                                                     
  integer, allocatable :: IRPEMTS(:),JRPEMTS(:),LRPEMTS(:)

  ! *** INTEGER SCALARS                                                                                                    
  integer   :: ISRPEM,INITRPEM,IRPEMWC,IRPEME,IJRPRS
  integer   :: ISTOXRPE,JSRPEM,ISRPEMSPAC,ISRPEMSPFR
  integer   :: ISRPEMTIME,ISRPEMTIFR,ISRPEMTILC,NSRPEMSPFR
  integer   :: NCRPEMRST,IRPEMSD
  integer   :: NRPEM, NRPEMEE
  integer   :: NRPEMSTEPS

  ! *** REAL SCALARS                                                                                                       
  real       :: HRPEMIC, RPSO, RPRO, RPEO, RPDO, RPSOC, RPEOC, RPSNC
  real       :: RPRNC, RPENC, RPSPC, RPRPC, RPEPC, PMRPS, FPRPR, RMRPS, RLRPS, FRPSD
  real       :: RKHNPRPS, RMRPR, RLRPR, RLRPD, RJRPRSC, RKRPORS, ROSR, RKRPRS, RISSS
  real       :: PMRPE, RMRPE, RLRPE, RKHNPRPE, RKHNRPS, RKHNRPR, RKHNRPE, RKHPRPS
  real       :: RKHPRPR, RKHPRPE, TP1RPS, TP2RPS, RKTP1RPS, RKTP2RPS, TP1RPE
  real       :: RKTP1RPE, RKTP2RPE, HRPS, HOPT, RKERPE, CCHLRPE, RISSOM, RISSOEM
  real       :: STOXS, STOXE, TR1RPS, TR2RPS, RKTR1RPS, RKTR2RPS, TR1RPR, TR2RPR
  real       :: RKTR1RPR, RKTR2RPR, TR1RPE, TR2RPE, RKTR1RPE, RKTR2RPE, FCRRPS
  real       :: FCLRPS, FCDRPS, FCRLRPS, FCLLRPS, FCDLRPS, FCRRPR, FCLRPR, FCDRPR
  real       :: FCRLRPR, FCLLRPR, FCDLRPR, FCRRPE, FCLRPE, FCDRPE, FCRLRPE
  real       :: FCLLRPE, FCDLRPE, FCRLRPD, FCLLRPD, FCDLRPD, FPRRPS, FPLRPS
  real       :: FPDRPS, FPIRPS, FPRLRPS, FPLLRPS, FPDLRPS, FPILRPS, FPRRPR, FPLRPR
  real       :: FPDRPR, FPIRPR, FPRLRPR, FPLLRPR, FPDLRPR, FPILRPR, FPRRPE, FPLRPE
  real       :: FPDRPE, FPIRPE, FPRLRPE, FPLLRPE, FPDLRPE, FPILRPE, FPRLRPD
  real       :: FPLLRPD, FPDLRPD, FPILRPD, FNRRPS, FNLRPS, FNDRPS, FNIRPS, FNRLRPS
  real       :: FNLLRPS, FNDLRPS, FNILRPS, FNRRPR, FNLRPR, FNDRPR, FNIRPR, FNRLRPR
  real       :: FNLLRPR, FNDLRPR, FNILRPR, FNRRPE, FNLRPE, FNDRPE, FNIRPE, FNRLRPE
  real       :: FNLLRPE, FNDLRPE, FNILRPE, FNRLRPD, FNLLRPD, FNDLRPD, FNILRPD
  real       :: FRPRRPG1, FRPRRPG2, FRPRRPG3
  real       :: RKSH, RKHI

  real,target, allocatable :: WQRPS(:),     WQRPR(:),     WQRPE(:),     WQRPD(:)
  real, allocatable :: PRPS(:),      RRPS(:),      RRPR(:),      PRPE(:),   RRPE(:)
  real, allocatable :: WQRPSR(:),    WQRPSL(:),    WQRPER(:),    WQRPEL(:), WQRPDL(:)
  real, allocatable :: WQRPSRP(:),   WQRPSLP(:),   WQRPERP(:),   WQRPELP(:)
  real, allocatable :: WQRPDLP(:),   WQRPSRN(:),   WQRPSLN(:),   WQRPERN(:)
  real, allocatable :: WQRPELN(:),   WQRPDLN(:),   WQRPRR(:),    WQRPRL(:)
  real, allocatable :: WQRPRRP(:),   WQRPRLP(:),   WQRPRRN(:),   WQRPRLN(:)
  real, allocatable :: FRPSPW(:),    FRPSNW(:),    PNRPS(:),     PNRPE(:),  RJRPRS(:)   
  real, allocatable :: XLIMTPRPS(:), XLIMTPRPE(:), XLIMTRRPS(:), XLIMTRRPE(:)
  real, allocatable :: XLIMTRRPR(:), XLIMNRPS(:),  XLIMNRPE(:),  XLIMLRPS(:)
  real, allocatable :: XLIMLRPE(:),  RISS(:),      RISSO(:),     RISSOE(:)
  real, allocatable :: RPEMTPRPS(:), RPEMTPRPE(:), RPEMTRRPS(:)
  real, allocatable :: RPEMTRRPE(:), RPEMTRRPR(:)

  ! *** LOGICAL                                                                                                            
  logical(4), allocatable :: LMASKRPEM(:)
  
  contains

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
                                            
  integer :: L,IWQTRPEM,K,LL,LF,LP,ND,IOBC,LN                         
  real    :: RATION,TOP,TMPNIT,RATIOP,TMPPO4,WQAVGIO,TMP1,RATIOHP,HDRY2    
  real    :: RKESSAVG,ALPHATOP,ALPHABOT,TMPEXP,SOURSINK,FACIMP,DTDHWQ 
  real    :: TMPWAT,TMPBED,TOP1,TOP2,BOT1,BOT2,BOT3,TMPNH4S,TMPNO3S
  real    :: TMPNH4E,TMPNO3E,WQKESS                                                            
  real    :: RKESSTOP(LCM),RKESSBOT(LCM)
  real    :: WQBCV(NBCSOP,NWQV)
  
  integer,save,allocatable,dimension(:)   :: NLRPEM   ! *** NUMBER OF RPEM ACTIVE CELLS FOR EACH LAYER BY DOMAIN
  integer,save,allocatable,dimension(:,:) :: LLRPEM   ! *** L INDEX FOR THE RPEM CELLS
  integer,save :: ICOUNT
  
  ! *** INPUT AND INITIALIZATION                                                                                           
  if( JSRPEM == 1 )then
    NSRPEMSPFR = 0                                                                                                         
    NCRPEMRST = 0                                                                                                          
    
    allocate(NLRPEM(NDM))
    allocate(LLRPEM(LCM,NDM))
    NLRPEM = 0
    LLRPEM = 0
    ICOUNT = 0
  endif
  ICOUNT = ICOUNT+1
  
  ! *** RPEM ONLY INTERACTS WITH THE BOTTOM LAYER K = KSZ(L)
  
  ! *** save VALUES AT OPEN BOUNDARIES
  if( IRPEMWC == 0 )then
    do IOBC = 1,NBCSOP  
      L = LOBCS(IOBC)  
      WQBCV(IOBC,1:NWQV) = WQV(L,KSZ(L),1:NWQV)
    enddo  
  endif
  HDRY2 = 2.*HDRY
  
  ! *** OBTAIN THE ACTIVE CELL LIST
  if( (ISDRY > 0 .and. LADRY > 0) .or. JSRPEM == 1 )then
    do ND = 1,NDM  
      LF = (ND-1)*LDMWET+1  
      LL = min(LF+LDMWET-1,LAWET)
      LN = 0
      do LP = LF,LL
        L = LWET(LP)
        if( LMASKRPEM(L) )then
          LN = LN + 1
          LLRPEM(LN,ND) = L
          if( ISICE > 2 )then
            if( ICECELL(L) .and. HP(L) < 3.*HDRY )then
              ! *** DEACTIVATE RPEM KINETICS IF THE CELL HAS ICE COVER AND VERY SHALLOW DEPTHS
              LN = LN - 1
            endif
          endif
        endif
      enddo
      NLRPEM(ND) = LN
    enddo
      
    ! *** Report only the fist instance
    if( JSRPEM == 1 )then
      do ND = 1,NDM  
        if( process_id == master_id ) WRITE(6,*) 'RPEM DOMAIN LIST (Thread,Count):', ND, NLRPEM(ND)  
        open(mpi_efdc_out_unit,FILE = OUTDIR//mpi_efdc_out_file,POSITION = 'APPEND')
        write(mpi_efdc_out_unit,*) 'RPEM DOMAIN LIST (Thread,Count):', ND, NLRPEM(ND) 
        close(mpi_efdc_out_unit)
      enddo
    endif
  endif

  ! *** Main loop over all the active RPEM cells
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LP, L, IWQTRPEM, K)              &
  !$OMP PRIVATE(RATION, TOP, TMPNIT, RATIOP, TMPPO4, WQAVGIO, TMP1)              &
  !$OMP PRIVATE(RKESSAVG, ALPHATOP, ALPHABOT, TMPEXP, SOURSINK, FACIMP, DTDHWQ)  & 
  !$OMP PRIVATE(TMPWAT, TMPBED, TOP1, TOP2, BOT1, BOT2, BOT3, TMPNH4S, TMPNO3S)  &
  !$OMP PRIVATE(TMPNH4E, TMPNO3E, WQKESS, RKESSTOP, RKESSBOT, RATIOHP)
  do ND = 1,NDM  
    ! **********************************************************************C                                                
    ! *** SET TEMPERATURE DEPENDENCY FOR GROWTH AND RESPIRATION 
    ! *** F3(T), EQ.(13), LOOK-UP TABLE
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      IWQTRPEM = NINT((TWQ(L)-WQTDMIN)/WQTDINC)+1
      if( IWQTRPEM < 1) IWQTRPEM = 1
      if( IWQTRPEM > NWQTD) IWQTRPEM = NWQTD
      XLIMTPRPS(L) = RPEMTPRPS(IWQTRPEM)  ! *** Shoot Production (Growth)                           
      XLIMTPRPE(L) = RPEMTPRPE(IWQTRPEM)  ! *** Epiphyte Production (Growth)
      XLIMTRRPS(L) = RPEMTRRPS(IWQTRPEM)  ! *** Shoot Respiration
      XLIMTRRPE(L) = RPEMTRRPE(IWQTRPEM)  ! *** Epiphyte Respiration
      XLIMTRRPR(L) = RPEMTRRPR(IWQTRPEM)  ! *** Root Respiration
    enddo                                                                                                                  
   
    ! **********************************************************************C                                                
    ! *** SET NUTRIENT LIMITATION FOR PLANT SHOOT GROWTH                                                                     
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      K = KSZ(L)
      RATION = RKHNRPS/RKHNRPR                                       ! Ratio of N half saturation constants for water/bed  
      TOP = WQV(L,K,INHX) + WQV(L,K,INOX) + RATION*(SM2NH4(L)+SM2NO3(L))                                                           
      TMPNIT = TOP/(TOP+RKHNRPS)
      RATIOP = RKHPRPS/RKHPRPR
      TOP = WQPO4D(L,K)+RATIOP*SM2PO4(L)                                                                                   
      TMPPO4 = TOP/(TOP+RKHPRPS)
      XLIMNRPS(L) = min(TMPNIT,TMPPO4)                                ! F1(N), EQ. (6), Minimum of N or P limit
    enddo                                                                                                                  
   
    ! **********************************************************************C                                                
    ! *** SET NUTRIENT LIMITATION FOR EPIPHYTE GROWTH                                                                        
    if( IRPEME > 0 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TOP = WQV(L,K,INHX) + WQV(L,K,INOX)                                                                                        
        TMPNIT = TOP/(TOP+RKHNRPE+1.E-18)                                                                                           
        TMPPO4 = WQPO4D(L,K)/(WQPO4D(L,K)+RKHPRPE+1.E-18)
        XLIMNRPE(L) = min(TMPNIT,TMPPO4)                                ! *** Minimum of N or P limits.  EQ. (20)
      enddo                                                                                                                  
    endif
     
    ! **********************************************************************C                                                
    ! *** SET LIGHT LIMITATIONS                                                                                              
   
    ! *** NOTE THAT WQI0,WQI1,AND WQI2 ARE PHOTOSYNTHETIC SSW RADIATION                                                      
    !    FROM CURRENT AND PREVIOUS TWO TIME INCREMENTS IN LANGLEYS/DAY                                                       
    !    AND ARE PROVIDED BY THE WATER QUALITY MODEL                                                                        
    WQAVGIO = WQCIA*WQI0+WQCIB*WQI1+WQCIC*WQI2                          ! EQ. (12)
   
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      ! *** COMPUTE TOTAL EXTINCTION COEFFICIENT

      ! *** Bottom PLUS 1 Layer
      if( KSZ(L) /= KC )then
        RKESSTOP(L) = RADKE(L,KSZ(L)+1)   ! EQ. (24)
      else
        RKESSTOP(L) = 0.                    ! EQ. (24)
      endif
              
      ! *** Bottom Layer
      WQKESS = RADKE(L,KSZ(L))
      RKESSBOT(L) = WQKESS
    enddo                                                                                                                  
   
    ! *** LIGHT LIMITATION SECTION
    if( IRPEME > 0 )then
      ! *** LIGHT LIMITATION FOR SHOOTS                                                                                        
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RKESSBOT(L) = RKESSBOT(L) + RKERPE*Z(L,KSZ(L))*HWQ(L)*WQRPE(L)/CCHLRPE   ! EQ. (10), INCLUDE EPIPHYTES SHADING
      enddo                                                                                                                  
    endif
   
    if( IWQSUN == 2 )then
      ! *** ASER TIME INTERVAL VARIABLE LIGHT
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSO(L) = RISSOM                                                ! Maximum value for optimum SSW for rooted plant growth (Langleys/day)
        RKESSAVG = 0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        WQAVGIO = PARADJ*2.065*RADTOP(L,KC)                              ! RADTOP INCLUDES SHADING AND ICE COVER
        ALPHATOP = -(WQAVGIO/RISSO(L)) *EXP(-RKESSTOP(L)*TMP1*HWQ(L))    ! EQ. (9), ASSUME HRPS IS WITHIN THE BOTTOM LA
        ALPHABOT = -(WQAVGIO/RISSO(L)) *EXP(-RKESSBOT(L)*HWQ(L))         ! EQ. (8)
        TMPEXP = EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPS(L) = 2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP         ! F2(I), EQ. (7)
        XLIMLRPS(L) = min(XLIMLRPS(L),1.0)
      enddo                                                                                                                  
    else
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSO(L) = RISSOM                                                ! Maximum value for optimum SSW for rooted plant growth (Langleys/day)
        RKESSAVG = 0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        ALPHATOP = -(WQAVGIO/RISSO(L)) *EXP(-RKESSTOP(L)*TMP1*HWQ(L))    ! EQ. (9), ASSUME HRPS IS WITHIN THE BOTTOM LA
        ALPHABOT = -(WQAVGIO/RISSO(L)) *EXP(-RKESSBOT(L)*HWQ(L))         ! EQ. (8)
        TMPEXP = EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPS(L) = 2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP         ! F2(I), EQ. (7), BUT PRACTICALLY, IT DOES NOT
        XLIMLRPS(L) = min(XLIMLRPS(L),1.0)
      enddo                                                                                                                  
    endif
    
    if( IRPEME > 0 )then
      ! *** LIGHT LIMITATION FOR EPIPHYTES                                                                                     
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSOE(L) = RISSOEM                                                                                                  
        RKESSAVG = 0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        ALPHATOP = -(WQAVGIO/RISSOE(L))*EXP(-RKESSTOP(L)*TMP1*HWQ(L))   ! EQ. (21), ASSUME HRPS IS WITHIN THE BOTTOM L
        ALPHABOT = -(WQAVGIO/RISSOE(L))*EXP(-RKESSBOT(L)*HWQ(L))        ! EQ. (22)
        TMPEXP = EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPE(L) = 2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP        ! EQ. (21)
        XLIMLRPE(L) = min(XLIMLRPE(L),1.0)
      enddo                                                                                                                  
    endif
  
    if( IJRPRS == 2 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        TMP1 = 1.-Z(L,KSZ(L))
        RISS(L) = WQAVGIO*EXP(-RKESSTOP(L)*TMP1*HWQ(L))                 ! EQ. (11), PART OF IT, USED IF IJRPRS = 2 IN CAL_RPEM.INP C5
                                                                        ! HOPT IS NOT USED HERE!, CAL_RPEM.INP C9
      enddo                                                                                                                  
    endif
     
    ! **********************************************************************C                                                
    ! *** UPDATE GROWTH AND RESPIRATION RATES                                                                                
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      
      ! *** Limitors     TEM          N/P         Light         Self Shading
      PRPS(L) = PMRPS*XLIMTPRPS(L)*XLIMNRPS(L)*XLIMLRPS(L)*EXP(-RKSH*WQRPS(L))  ! EQ. (5)   PRPS - Nutrient, light and temperature limited growth rate (/day)
                                                                                !                  Include self-shading of shoots, JI, 7/4/04
      RRPS(L) = RMRPS*XLIMTRRPS(L)                                              ! EQ. (15)  RRPS - Temperature limited shoot respiration (/day)
    enddo                                                                                                                  
   
    ! *** Limit RPEM with depth
    if( ISDRY > 0 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RATIOHP = min( (MAX(HP(L)-HDRY,0.))/HDRY2, 1.0 )
        PRPS(L) = PRPS(L)*RATIOHP
        RRPS(L) = RRPS(L)*RATIOHP
      enddo    
    endif
    
    ! *** ROOT RESPIRATION                                                                                                   
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      RRPR(L) = RMRPR*XLIMTRRPR(L)                                    ! EQ. (18)
    enddo                                                                                                                  
   
    ! *** EPIPHYTE GROWTH AND RESPIRATION                                                                                    
    if( IRPEME > 0 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        PRPE(L) = PMRPE*XLIMTPRPE(L)*XLIMNRPE(L)*XLIMLRPE(L)            ! EQ. (19), F3T*F1N*F2I
        RRPE(L) = RMRPE*XLIMTRRPE(L)                                                                                         
      enddo                                                                                                                  
    endif
     
    ! **********************************************************************C                                                
   
    ! *** UPDATE ROOT TO SHOOT FLUX                                                                                          
    if( IJRPRS == 0 .and. N < 5 )then                      ! *** IJRPRS = 0, CONSTANT TRANSPORT RATE
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RJRPRSC
      enddo                                                                                                                  
    endif                                                                                                                  
   
    if( IJRPRS == 1 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RKRPORS*(ROSR*WQRPR(L)-WQRPS(L))                                                                         
      enddo                                                                                                                  
    endif                                                                                                                  
   
    if( IJRPRS == 2 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RKRPRS*RISS(L)/(RISS(L)+RISSS+1E-18)                                                                           
      enddo                                                                                                                  
    endif                                                                                                                  
   
    ! **********************************************************************C                                                
    ! ***  UPDATE SHOOT, ROOT, EPIPHYTE AND DETRITUS STATE VARIABLES                                                         
   
    ! *** UPDATE SHOOTS BIOMASS                                                                                                      
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      ! ***            Growth     Respiration   Other
      SOURSINK = (1.-FPRPR)*PRPS(L) - RRPS(L) - RLRPS                     ! EQ. (1)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPS(L) = FACIMP*(WQRPS(L) + DTWQ*RJRPRS(L))
      WQRPS(L) = max(WQRPS(L),0.2)   ! KEEP THE "SEED", SINCE IF RPS! = 0, IT WILL NEVER GROW AGAIN, ACCORDING TO EQ.
    enddo                                                                                                                  
   
    ! *** UPDATE ROOTS BIOMASS                                                                                                      
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      SOURSINK = -RRPR(L)-RLRPR                                       ! EQ. (2)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPR(L) = FACIMP*(WQRPR(L)+DTWQ*FPRPR*PRPS(L)*WQRPS(L)-DTWQ*RJRPRS(L))
    enddo                                                                                                                  

    !------------------------------                                                                                         
    !GO TO 904  ! SKIP EPHIPHYTES AND DETRITUS
    
    ! *** UPDATE EPIPHYTES                                                                                                   
    if( IRPEME > 0 )then
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        SOURSINK = PRPE(L)-RRPE(L)-RLRPE                              ! EQ. (3)
        FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
        WQRPE(L) = FACIMP*WQRPE(L)
      enddo                                                                                                                  
    endif
     
    ! *** UPDATE SHOOT DETRITUS IN WATER COLUMN                                                                              
    SOURSINK = -RLRPD                                                    ! EQ. (4)  RLRPD-Loss rate for plant detritus at bottom of water column (/day)
    do LP = 1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPD(L) = FACIMP*(WQRPD(L)+DTWQ*FRPSD*RLRPS*WQRPS(L))
    enddo                                                                                                                  
    904   continue                                                                                                          
   
    ! **********************************************************************C                                                
    ! *** CALCULATE WATER COLUMN SHOOT, EPIPHYTE AND DETRITUS                                                                
    ! *** RESPIRATION AND NON-RESPIRATION LOSSES TO WATER QUALITY ORGANIC                                                    
    ! *** MATTER STATE VARIABLES AND TO PHOSPHATE AND AMMONIA                                                                
   
    !     GO TO 901   ! SKIP COUPLING WITH NUTRIENTS AND DO IN WATER COLUMN                                    
    if( IRPEMWC == 0 )then

      ! *** UPDATE SHOOT, ROOT AND DETRITUS
      do LP = 1,NLRPEM(ND)
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
      enddo
                                                                                      
      ! *** UPDATE EPIPHYTES
      if( IRPEME > 0 )then
        do LP = 1,NLRPEM(ND)
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
        enddo                                                                                                                  
      endif
     
      ! *** UPDATE WATER COLUMN ORGANIC CARBON, PHOSPHOROUS AND NITROGEN                                                       
      do LP = 1,NLRPEM(ND)
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
     
      enddo                                                                                                                  
     
      ! **********************************************************************C                                                
      ! *** CALCULATE WATER COLUMN SHOOT AND EPIPHYTE UPTAKE OF                                                                
      ! *** DISSOLVED PHOSPHATE, AMMONIA, AND NO3                                                                              
     
      ! *** CALCULATE THE FRACTION OF DISSOLVED PHOSPHATE UPTAKE FROM WATER                                                    
      ! *** COLUMN (FRPSPW) BY PLANT SHOOTS                                                                                    
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TMPWAT = RKHPRPR*WQV(L,K,10)                                    ! ??, BUT (38) MEANS PO4DW, DISSOLVED ONLY, NO
        TMPBED = RKHPRPS*SM2PO4(L)                                                                                           
        FRPSPW(L) = TMPWAT/(TMPWAT+TMPBED+1E-18)                        ! EQ. (38)
      enddo                                                                                                                  
   
      ! *** UPDATE TOTAL PHOSPHATE IN WATER COLUMN                                                                             
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))                                                                                           
        WQV(L,K,IP4D) = WQV(L,K,IP4D) - DTDHWQ*( RPEPC*PRPE(L)*WQRPE(L) + RPSPC*FRPSPW(L)*PRPS(L)*WQRPS(L) )  ! PO4T, EQ. COMPLETED
        WQV(L,K,IP4D) = max(WQV(L,K,IP4D),0.)
      enddo                                                                                                                  
     
      ! *** CALCULATE THE FRACTION OF AMMONIA AND NO3 UPTAKE FROM WATER                                                        
      ! *** COLUMN (FRPSNW) BY PLANT SHOOTS                                                                                    
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TMPWAT = RKHNRPR*(WQV(L,K,INHX) + WQV(L,K,INOX))
        TMPBED = RKHNRPS*(SM2NH4(L)+SM2NO3(L))                                                                               
        FRPSNW(L) = TMPWAT/(TMPWAT+TMPBED+1E-18)                              ! EQ. (45)
      enddo                                                                                                                  
     
      ! *** CALCULATE THE AMMONIA PREFERENCE FOR SHOOTS                                                                        
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TOP1 = (WQV(L,K,INHX) + SM2NH4(L))*(WQV(L,K,INOX) + SM2NO3(L))                ! WHY USE WATER COLUMN AND BED, BOTH OF THEM?
        TOP2 = RKHNPRPS*(WQV(L,K,INHX) + SM2NH4(L))                                                                              
        BOT1 = RKHNPRPS + (WQV(L,K,INHX) + SM2NH4(L))                                                                              
        BOT2 = RKHNPRPS + (WQV(L,K,INOX) + SM2NO3(L))                                                                              
        BOT3 = WQV(L,K,INHX) + SM2NH4(L) + (WQV(L,K,INOX) + SM2NO3(L))     
        BOT1 = max(BOT1,1E-12)                                                            
        BOT2 = max(BOT2,1E-12)                                                            
        BOT3 = max(BOT3,1E-12)                                                            
        PNRPS(L) = TOP1/(BOT1*BOT2) + TOP2/(BOT2*BOT3)                           ! EQ. (44A)
      enddo                                                                                                                  
   
      ! *** CALCULATE THE AMMONIA PREFERENCE FOR EPIPHYTES                                                                     
      if( IRPEME > 0 )then
        do LP = 1,NLRPEM(ND)
          L = LLRPEM(LP,ND)
          K = KSZ(L)
          ! ***    NH3        NO3                                                                                                 
          TOP1 = WQV(L,K,INHX)*WQV(L,K,INOX)                                                                                       
          TOP2 = RKHNPRPE*WQV(L,K,INHX)                                                                                          
          BOT1 = RKHNPRPE+WQV(L,K,INHX)                                                                                          
          BOT2 = RKHNPRPE+WQV(L,K,INOX)                                                                                          
          BOT3 = WQV(L,K,INHX)+WQV(L,K,INOX)     
          BOT1 = max(BOT1,1E-12)                                                               
          BOT2 = max(BOT2,1E-12)                                                               
          BOT3 = max(BOT3,1E-12)                                                               
          PNRPE(L) = TOP1/(BOT1*BOT2)+TOP2/(BOT2*BOT3)                    ! EQ. (44B)                                                                    
        enddo                                                                                                                  
      endif
     
      ! *** UPDATE WATER COLUMN AMMONIA AND NO3 NITROGEN                                                                       
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))                                                                                           
        TMPNH4S = PNRPS(L)                                                                                                   
        TMPNO3S = 1.-PNRPS(L)                                                                                                
        TMPNH4E = PNRPE(L)                                                                                                   
        TMPNO3E = 1.-PNRPE(L)                                                                                                
        WQV(L,K,INHX) = WQV(L,K,INHX) - DTDHWQ*(RPENC*TMPNH4E*PRPE(L)*WQRPE(L) + RPSNC*TMPNH4S*FRPSNW(L)*PRPS(L)*WQRPS(L))  ! NH4, EQ. COMPLETED
        WQV(L,K,INOX) = WQV(L,K,INOX) - DTDHWQ*(RPENC*TMPNO3E*PRPE(L)*WQRPE(L) + RPSNC*TMPNO3S*FRPSNW(L)*PRPS(L)*WQRPS(L))  ! NO3, EQ. COMPLETED
        WQV(L,K,INHX) = max(WQV(L,K,INHX),0.)
        WQV(L,K,INOX) = max(WQV(L,K,INOX),0.)
      enddo                                                                                                                  
     
      ! **********************************************************************C                                                
      ! *** CALCULATE WATER COLUMN SOURCES OF DISSOLVED OXYGEN DUE                                                             
      ! *** TO SHOOT AND EPIPHYTE GROWTH                                                                                       
      do LP = 1,NLRPEM(ND)
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
        WQV(L,K,IDOX) = max(WQV(L,K,IDOX),0.0)
      enddo                                                                                                                  
    endif    ! *** END OF BLOCK FOR SKIPPING WATER COLUMN LINKAGE
   
    !901   continue                                                                                                          
   
    ! **********************************************************************C                                                
    ! *** CALCULATE SEDIMENT BED ROOT RESPIRATION AND NON-RESPIRATION                                                        
    ! *** LOSSES TO SEDIMENT DIAGENESIS ORGANIC MATTER STATE VARIABLES                                                       
    ! *** AND TO PHOSPHATE AND AMMONIA                                                                                       
   
    !     GO TO 902   ! SKIP COUPLING WITH THE SEDIMENT DIAGENESIS MODEL.  FULL DIAGENESIS MUST BE ON TO COUPLE
    
    if( IRPEMSD == 0 .and. IWQBEN == 1 )then   
      ! *** LOSSES TO ORGANIC CARBON                                                                                           
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRR(L) = RRPR(L)*WQRPR(L)                                    ! FROM EQ. (32B)
        WQRPRL(L) = RLRPR*WQRPR(L)                                                                                           
      enddo                                                                                                                  
     
      ! *** LOSSES TO ORGANIC PHOSPHOROUS                                                                                      
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRRP(L) = RPRPC*WQRPRR(L)                                    !  ! FROM EQ. (34B)
        WQRPRLP(L) = RPRPC*WQRPRL(L)                                                                                         
      enddo                                                                                                                  
     
      ! *** LOSSES TO ORGANIC NITROGEN                                                                                         
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRRN(L) = RPRNC*WQRPRR(L)                                    ! FROM EQ. (39B)
        WQRPRLN(L) = RPRNC*WQRPRL(L)                                                                                         
      enddo                                                                                                                  
     
      ! *** UPDATE SEDIMENT BED ORGANIC CARBON, PHOSPHOROUS AND NITROGEN                                                       
      do LP = 1,NLRPEM(ND)
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
      enddo                                                                                                                  
     
      ! **********************************************************************C                                                
      ! *** CALCULATE SEDIMENT BED ROOT UPTAKE OF DISSOLVED PHOSPHATE,                                                         
      ! *** AMMONIA AND NO3                                                                                                    
     
      ! *** UPDATE TOTAL PHOSPHATE IN SEDIMENT BED                                                                             
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))
        SM2PO4(L) = SM2PO4(L)-DTDHWQ*RPRPC*(1.-FRPSPW(L))*PRPS(L)*WQRPS(L)    ! EQ. (37B), COMPLETED
        SM2PO4(L) = max(SM2PO4(L), 1.E-12)
      enddo                                                                                                                  
     
      ! *** UPDATE SEDIMENT BED AMMONIA AND NO3 NITROGEN                                                                       
      do LP = 1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ    = DTWQ/(DZC(L,K)*HWQ(L))
        TMPNH4S   = PNRPS(L)
        TMPNO3S   = 1. - PNRPS(L)                                                                                                
        TMPBED    = 1. - FRPSNW(L)                                                                                                
        SM2NH4(L) = SM2NH4(L) - DTDHWQ*RPSNC*TMPNH4S*TMPBED*PRPS(L)*WQRPS(L)   ! EQ. (42B), COMPLETED
        SM2NO3(L) = SM2NO3(L) - DTDHWQ*RPSNC*TMPNO3S*TMPBED*PRPS(L)*WQRPS(L)   ! EQ. (43B)
        SM2NH4(L) = max(SM2NH4(L), 1.E-12)
        SM2NO3(L) = max(SM2NO3(L), 1.E-12)
      enddo                                                                                                                  
    endif   ! *** END OF SEDIMENT BED LINKAGE BYPASS
  
    902   continue                                                                                                          

  enddo    ! *** END OF DOMAIN LOOP 
  !$OMP END PARALLEL DO
   
  ! **********************************************************************C 
  ! *** RESTORE OPEN BOUNDARY CONCENTRATIONS                                               
  if( IRPEMWC == 0 )then
    do IOBC = 1,NBCSOP  
      L = LOBCS(IOBC)
      K = KSZ(L)
      WQV(L,K,1:NWQV) = WQBCV(IOBC,1:NWQV)
    enddo  
  endif
  
   
  ! **********************************************************************C                                                
  JSRPEM = 0                                                                                                             

  return

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

  use fson
  use mod_fson_value, only: fson_value_count, fson_value_get
  use INFOMOD,only:SKIPCOM,READSTR

  Character(len = 79), allocatable :: TITLE(:)
  Character*80 STR*200
  type(fson_value), Pointer :: json_data, item
  
  integer :: L,LDATA,LL,I,J,NT,LG
  real    :: SMNH4,SMNO3,SMPO4,TP2RPE,RPST,RPRT,RPET,RPDT,TIMETMP,WTEMP 
  
  Real,allocatable :: FRPRRPG_TEMP(:)
  
  if( process_id == master_id )then
    json_data => fson_parse("wq_rpem.jnp")
    Write(*,'(A)') ' WQ: WQRPEM.JNP - RPEM CONTROL FILE'
         
    !Call fson_get(json_data, "title", TITLE)
    call fson_get(json_data, "initial_condition_option", INITRPEM)
    call fson_get(json_data, "water_column_nutrient_linkage", IRPEMWC)
    call fson_get(json_data, "sediment_diagenesis_nutrient_linkage", IRPEMSD)
    call fson_get(json_data, "initial_carbon_biomass.shoot", RPSO)
    call fson_get(json_data, "initial_carbon_biomass.root", RPRO)
    call fson_get(json_data, "initial_carbon_biomass.epiphyte", RPEO)
    call fson_get(json_data, "initial_carbon_biomass.detritus", RPDO)
    
    call fson_get(json_data, "const_bed_porewater_concentration.ammonia", SMNH4)
    call fson_get(json_data, "const_bed_porewater_concentration.nitrate", SMNO3)
    call fson_get(json_data, "const_bed_porewater_concentration.phosphate", SMPO4)
    
    call fson_get(json_data, "rooted_plant.average_shoot_height_above_bed", HRPS)
    call fson_get(json_data, "rooted_plant.stoichiometric_ratios.root.nitrogen_to_carbon", RPRNC)
    call fson_get(json_data, "rooted_plant.stoichiometric_ratios.root.phosphorus_to_carbon", RPRPC)
    call fson_get(json_data, "rooted_plant.stoichiometric_ratios.shoot.oxygen_to_carbon", RPSOC)
    call fson_get(json_data, "rooted_plant.stoichiometric_ratios.shoot.nitrogen_to_carbon", RPSNC)
    call fson_get(json_data, "rooted_plant.stoichiometric_ratios.shoot.phosphorus_to_carbon", RPSPC)
    
    
    call fson_get(json_data, "rooted_plant.respiration.root.max_respiration_rate", RMRPR)
    call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.suboptimal_temperature", TR1RPR)
    call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.superoptimal_temperature", TR2RPR)
    call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.suboptimal_temperature_effect", rKTR1RPR)
    call fson_get(json_data, "rooted_plant.respiration.root.temperature_effects.superoptimal_temperature_effect", rKTR2RPR)
    
    
    call fson_get(json_data, "rooted_plant.respiration.shoot.max_respiration_rate", RMRPS)
    call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.suboptimal_temperature", TR1RPS)
    call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.superoptimal_temperature", TR2RPS)
    call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.suboptimal_temperature_effect", rKTR1RPS)
    call fson_get(json_data, "rooted_plant.respiration.shoot.temperature_effects.superoptimal_temperature_effect", rKTR2RPS)
    
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_carbon.RPOC", FCRrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_carbon.LPOC", FCLrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_carbon.DOC", FCDrpr)
    
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.RPOP", FPRrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.LPOP", FPLrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.DOP", FPDrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_phosphorus.PO4", FPIrpr)
    
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.RPON", FNRrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.LPON", FNLrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.DON", FNDrpr)
    call fson_get(json_data, "rooted_plant.excretion.root.fraction_of_nitrogen.NH4", FNIrpr)
    
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_carbon.RPOC", FCRrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_carbon.LPOC", FCLrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_carbon.DOC", FCDrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.RPOP", FPRrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.LPOP", FPLrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.DOP", FPDrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_phosphorus.PO4", FPIrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.RPON", FNRrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.LPON", FNLrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.DON", FNDrps)
    call fson_get(json_data, "rooted_plant.excretion.shoot.fraction_of_nitrogen.NH4", FNIrps)
    
    call fson_get(json_data, "rooted_plant.growth.max_rooted_depth", HRPEMIC)
    call fson_get(json_data, "rooted_plant.growth.optimum_water_depth_for_growth", HOPT)
    call fson_get(json_data, "rooted_plant.growth.max_SSW_for_optimum_growth", rISSOM)
    call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.transport_option", iJRPRS)
    call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.const_positive_carbon_transport", rJRPRSC)
    call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.transfer_rate_ratio", rKRPORS)
    call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.observed_ratio", ROSR)
    call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.light_dependent_transfer_rate", rKRPRS)
    call fson_get(json_data, "rooted_plant.growth.root_to_shoot_carbon_transport.half_sat_solar_ratio_at_surface", rISSS)
    
    call fson_get(json_data, "rooted_plant.growth.root.nutrient_limits.half_sat_const_for_nitrogen_uptake_from_bed", rKHNRPR)
    call fson_get(json_data, "rooted_plant.growth.root.nutrient_limits.half_sat_const_for_phosphorus_uptake_from_bed", rKHPRPR)
    
    call fson_get(json_data, "rooted_plant.growth.shoot.max_growth_rate", PMRPS)
    call fson_get(json_data, "rooted_plant.growth.shoot.fraction_of_production_transferred_to_roots", FPRPR)
    call fson_get(json_data, "rooted_plant.growth.shoot.half_sat_nitrogen_preference", rKHNPRPS)
    call fson_get(json_data, "rooted_plant.growth.shoot.self_shading", rKSH)
    call fson_get(json_data, "rooted_plant.growth.shoot.half_sat_const_for_irradiance", rKHI)
    call fson_get(json_data, "rooted_plant.growth.shoot.nutrient_limits.half_sat_const_for_nitrogen_uptake_from_water", rKHNRPS)
    call fson_get(json_data, "rooted_plant.growth.shoot.nutrient_limits.half_sat_const_for_phosphorus_uptake_from_water", rKHPRPS)
    
    
    call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.suboptimal_temperature", TP1RPS)
    call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.superoptimal_temperature", TP2RPS)
    call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.suboptimal_temperature_effect", rKTP1RPS)
    call fson_get(json_data, "rooted_plant.growth.shoot.temperature_effects.superoptimal_temperature_effect", rKTP2RPS)
    
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.loss_rate", rLRPR)
    
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_carbon.RPOC", FCRLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_carbon.LPOC", FCLLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_carbon.DOC", FCDLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.RPOP", FPRLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.LPOP", FPLLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.DOP", FPDLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_phosphorus.PO4", FPILrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.RPON", FNRLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.LPON", FNLLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.DON", FNDLrpr)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.root.fraction_of_nitrogen.NH4", FNILrpr)
    
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.loss_rate", rLRPS)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_to_detritus", FRPSD)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_carbon.RPOC", FCRLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_carbon.LPOC", FCLLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_carbon.DOC", FCDLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.RPOP", FPRLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.LPOP", FPLLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.DOP", FPDLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_phosphorus.PO4", FPILrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.RPON", FNRLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.LPON", FNLLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.DON", FNDLrps)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.shoot.fraction_of_nitrogen.NH4", FNILrps)
    
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.loss_rate_at_bottom", rLRPD)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_carbon.RPOC", FCRLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_carbon.LPOC", FCLLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_carbon.DOC", FCDLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.RPOP", FPRLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.LPOP", FPLLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.DOP", FPDLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_phosphorus.PO4", FPILrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.RPON", FNRLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.LPON", FNLLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.DON", FNDLrpd)
    call fson_get(json_data, "rooted_plant.non_respiration_loss.detritus.fraction_of_nitrogen.NH4", FNILrpd)
    
    call fson_get(json_data, "rooted_plant.transfer_of_root_generated_RPOM_to_sediment", FRPRRPG_TEMP)
    FRPRRPG1 = FRPRRPG_TEMP(1)
    FRPRRPG2 = FRPRRPG_TEMP(2)
    FRPRRPG3 = FRPRRPG_TEMP(3)
    
    
    call fson_get(json_data, "epiphyte.stoichiometric_ratios.oxygen_to_carbon", RPEOC)
    call fson_get(json_data, "epiphyte.stoichiometric_ratios.nitrogen_to_carbon", RPENC)
    call fson_get(json_data, "epiphyte.stoichiometric_ratios.phosphorus_to_carbon", RPEPC)
    call fson_get(json_data, "epiphyte.stoichiometric_ratios.carbon_to_chlorophyll_ratio", CChlRPE)
    
    call fson_get(json_data, "epiphyte.respiration.max_respiration_rate", RMRPE)
    call fson_get(json_data, "epiphyte.respiration.temperature_effects.suboptimal_temperature", TR1RPE)
    call fson_get(json_data, "epiphyte.respiration.temperature_effects.superoptimal_temperature", TR2RPE)
    call fson_get(json_data, "epiphyte.respiration.temperature_effects.suboptimal_temperature_effect", rKTR1RPE)
    call fson_get(json_data, "epiphyte.respiration.temperature_effects.superoptimal_temperature_effect", rKTR2RPE)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_carbon.RPOC", FCRrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_carbon.LPOC", FCLrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_carbon.DOC", FCDrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.RPOP", FPRrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.LPOP", FPLrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.DOP", FPDrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_phosphorus.PO4", FPIrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.RPON", FNRrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.LPON", FNLrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.DON", FNDrpe)
    call fson_get(json_data, "epiphyte.excretion.fraction_of_nitrogen.NH4", FNIrpe)
    
    call fson_get(json_data, "epiphyte.growth.growing_option", IRPEME)
    call fson_get(json_data, "epiphyte.growth.max_growth_rate", PMRPE)
    call fson_get(json_data, "epiphyte.growth.half_sat_nitrogen_preference", rKHNPRPE)
    call fson_get(json_data, "epiphyte.growth.max_SSW_for_optimum_growth", rISSOEM)
    call fson_get(json_data, "epiphyte.growth.light_extinction_coef_for_chlorophyll", rKeRPE)
    call fson_get(json_data, "epiphyte.growth.nutrient_limits.half_sat_const_for_nitrogen_uptake", rKHNRPE)
    call fson_get(json_data, "epiphyte.growth.nutrient_limits.half_sat_const_for_phosphorus_uptake", rKHPRPE)
    
    call fson_get(json_data, "epiphyte.growth.temperature_effects.suboptimal_temperature", TP1RPE)
    call fson_get(json_data, "epiphyte.growth.temperature_effects.superoptimal_temperature", TP2RPE)
    call fson_get(json_data, "epiphyte.growth.temperature_effects.suboptimal_temperature_effect", rKTP1RPE)
    call fson_get(json_data, "epiphyte.growth.temperature_effects.superoptimal_temperature_effect", rKTP2RPE)
    
    call fson_get(json_data, "epiphyte.non_respiration_loss.loss_rate", rLRPE)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_carbon.RPOC", FCRLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_carbon.LPOC", FCLLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_carbon.DOC", FCDLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.RPOP", FPRLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.LPOP", FPLLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.DOP", FPDLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_phosphorus.PO4", FPILrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.RPON", FNRLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.LPON", FNLLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.DON", FNDLrpe)
    call fson_get(json_data, "epiphyte.non_respiration_loss.fraction_of_nitrogen.NH4", FNILrpe)
    
    call fson_get(json_data, "salinity_toxicity.salinity_effect_option", iSTOXRPE)
    call fson_get(json_data, "salinity_toxicity.STOXS", STOXS)
    call fson_get(json_data, "salinity_toxicity.STOXE", STOXE)
    
    call fson_get(json_data, "output.number_of_EE_WC_writes", NRPEMEE)
    call fson_get(json_data, "output.spatial_output.output_option", ISRPEMSPAC)
    call fson_get(json_data, "output.spatial_output.output_frequency", ISRPEMSPFR)
    
    call fson_get(json_data, "output.time_series_output.output_option", ISRPEMTIME)
    call fson_get(json_data, "output.time_series_output.output_frequency", ISRPEMTIFR)
    call fson_get(json_data, "output.time_series_output.number_of_locations", ISRPEMTILC)
  
  endif   ! *** End of master_id
  
  ! *** Scalar Variables
  call Broadcast_Scalar(INITRPEM, master_id)
  call Broadcast_Scalar(IRPEMWC,  master_id)
  call Broadcast_Scalar(IRPEMSD,  master_id)
  call Broadcast_Scalar(HRPEMIC,  master_id)
  call Broadcast_Scalar(RPSO,     master_id)
  call Broadcast_Scalar(RPRO,     master_id)
  call Broadcast_Scalar(RPEO,     master_id)
  call Broadcast_Scalar(RPDO,     master_id)
  call Broadcast_Scalar(SMNH4,    master_id)
  call Broadcast_Scalar(SMNO3,    master_id)
  call Broadcast_Scalar(SMPO4,    master_id)
  
  call Broadcast_Scalar(RPSOC,    master_id)
  call Broadcast_Scalar(RPEOC,    master_id)
  call Broadcast_Scalar(RPSNC,    master_id)
  call Broadcast_Scalar(RPRNC,    master_id)
  call Broadcast_Scalar(RPENC,    master_id)
  call Broadcast_Scalar(RPSPC,    master_id)
  call Broadcast_Scalar(RPRPC,    master_id)
  call Broadcast_Scalar(RPEPC,    master_id)
  
  call Broadcast_Scalar(PMRPS,    master_id)
  call Broadcast_Scalar(FPRPR,    master_id)
  call Broadcast_Scalar(RMRPS,    master_id)
  call Broadcast_Scalar(rLRPS,    master_id)
  call Broadcast_Scalar(FRPSD,    master_id)
  call Broadcast_Scalar(rKHNPRPS, master_id)
  call Broadcast_Scalar(rKSH,     master_id)
  call Broadcast_Scalar(rKHI,     master_id)
  
  call Broadcast_Scalar(RMRPR,    master_id)
  call Broadcast_Scalar(rLRPR,    master_id)
  call Broadcast_Scalar(rLRPD,    master_id)
  
  call Broadcast_Scalar(iJRPRS,   master_id)
  call Broadcast_Scalar(rJRPRSC,  master_id)
  call Broadcast_Scalar(rKRPORS,  master_id)
  call Broadcast_Scalar(ROSR,     master_id)
  call Broadcast_Scalar(rKRPRS,   master_id)
  call Broadcast_Scalar(rISSS,    master_id)
  
  call Broadcast_Scalar(IRPEME,   master_id)
  call Broadcast_Scalar(PMRPE,    master_id)
  call Broadcast_Scalar(RMRPE,    master_id)
  call Broadcast_Scalar(rLRPE,    master_id)
  call Broadcast_Scalar(rKHNPRPE, master_id)
  
  call Broadcast_Scalar(rKHNRPS,  master_id)
  call Broadcast_Scalar(rKHNRPR,  master_id)
  call Broadcast_Scalar(rKHNRPE,  master_id)
  call Broadcast_Scalar(rKHPRPS,  master_id)
  call Broadcast_Scalar(rKHPRPR,  master_id)
  call Broadcast_Scalar(rKHPRPE,  master_id)
  
  call Broadcast_Scalar(TP1RPS,   master_id)
  call Broadcast_Scalar(TP2RPS,   master_id)
  call Broadcast_Scalar(rKTP1RPS, master_id)
  call Broadcast_Scalar(rKTP2RPS, master_id)
  call Broadcast_Scalar(TP1RPE,   master_id)
  call Broadcast_Scalar(TP2RPE,   master_id)
  call Broadcast_Scalar(rKTP1RPE, master_id)
  call Broadcast_Scalar(rKTP2RPE, master_id)
  
  call Broadcast_Scalar(HRPS,     master_id)
  call Broadcast_Scalar(HOPT,     master_id)
  call Broadcast_Scalar(rKeRPE,   master_id)
  call Broadcast_Scalar(CChlRPE,  master_id)
  call Broadcast_Scalar(rISSOM,   master_id)
  call Broadcast_Scalar(rISSOEM,  master_id)
  
  call Broadcast_Scalar(iSTOXRPE, master_id)
  call Broadcast_Scalar(STOXS,    master_id)
  call Broadcast_Scalar(STOXE,    master_id)
  
  call Broadcast_Scalar(TR1RPS,   master_id)
  call Broadcast_Scalar(TR2RPS,   master_id)
  call Broadcast_Scalar(rKTR1RPS, master_id)
  call Broadcast_Scalar(rKTR2RPS, master_id)
  call Broadcast_Scalar(TR1RPR,   master_id)
  call Broadcast_Scalar(TR2RPR,   master_id)
  call Broadcast_Scalar(rKTR1RPR, master_id)
  call Broadcast_Scalar(rKTR2RPR, master_id)
  
  call Broadcast_Scalar(TR1RPE,   master_id)
  call Broadcast_Scalar(TR2RPE,  master_id)
  call Broadcast_Scalar(rKTR1RPE,master_id)
  call Broadcast_Scalar(rKTR2RPE,master_id)
  
  call Broadcast_Scalar(FCRrps,   master_id)
  call Broadcast_Scalar(FCLrps,   master_id)
  call Broadcast_Scalar(FCDrps,   master_id)
  call Broadcast_Scalar(FCRLrps,  master_id)
  call Broadcast_Scalar(FCLLrps,  master_id)
  call Broadcast_Scalar(FCDLrps,  master_id)
  
  call Broadcast_Scalar(FCRrpr,   master_id)
  call Broadcast_Scalar(FCLrpr,   master_id)
  call Broadcast_Scalar(FCDrpr,   master_id)
  call Broadcast_Scalar(FCRLrpr,  master_id)
  call Broadcast_Scalar(FCLLrpr,  master_id)
  call Broadcast_Scalar(FCDLrpr,  master_id)
  
  call Broadcast_Scalar(FCRrpe,   master_id)
  call Broadcast_Scalar(FCLrpe,   master_id)
  call Broadcast_Scalar(FCDrpe,   master_id)
  call Broadcast_Scalar(FCRLrpe,  master_id)
  call Broadcast_Scalar(FCLLrpe,  master_id)
  call Broadcast_Scalar(FCDLrpe,  master_id)
  
  call Broadcast_Scalar(FCRLrpd,  master_id)
  call Broadcast_Scalar(FCLLrpd,  master_id)
  call Broadcast_Scalar(FCDLrpd,  master_id)
  
  call Broadcast_Scalar(FPRrps,   master_id)
  call Broadcast_Scalar(FPLrps,   master_id)
  call Broadcast_Scalar(FPDrps,   master_id)
  call Broadcast_Scalar(FPIrps,   master_id)
  call Broadcast_Scalar(FPRLrps,  master_id)
  call Broadcast_Scalar(FPLLrps,  master_id)
  call Broadcast_Scalar(FPDLrps,  master_id)
  call Broadcast_Scalar(FPILrps,  master_id)
  
  call Broadcast_Scalar(FPRrpr,   master_id)
  call Broadcast_Scalar(FPLrpr,   master_id)
  call Broadcast_Scalar(FPDrpr,   master_id)
  call Broadcast_Scalar(FPIrpr,   master_id)
  call Broadcast_Scalar(FPRLrpr,  master_id)
  call Broadcast_Scalar(FPLLrpr,  master_id)
  call Broadcast_Scalar(FPDLrpr,  master_id)
  call Broadcast_Scalar(FPILrpr,  master_id)
  
  call Broadcast_Scalar(FPRrpe,   master_id)
  call Broadcast_Scalar(FPLrpe,   master_id)
  call Broadcast_Scalar(FPDrpe,   master_id)
  call Broadcast_Scalar(FPIrpe,   master_id)
  call Broadcast_Scalar(FPRLrpe,  master_id)
  call Broadcast_Scalar(FPLLrpe,  master_id)
  call Broadcast_Scalar(FPDLrpe,  master_id)
  call Broadcast_Scalar(FPILrpe,  master_id)
  
  call Broadcast_Scalar(FPRLrpd,  master_id)
  call Broadcast_Scalar(FPLLrpd,  master_id)
  call Broadcast_Scalar(FPDLrpd,  master_id)
  call Broadcast_Scalar(FPILrpd,  master_id)
  
  call Broadcast_Scalar(FNRrps,   master_id)
  call Broadcast_Scalar(FNLrps,   master_id)
  call Broadcast_Scalar(FNDrps,   master_id)
  call Broadcast_Scalar(FNIrps,   master_id)
  call Broadcast_Scalar(FNRLrps,  master_id)
  call Broadcast_Scalar(FNLLrps,  master_id)
  call Broadcast_Scalar(FNDLrps,  master_id)
  call Broadcast_Scalar(FNILrps,  master_id)
  
  call Broadcast_Scalar(FNRrpr,   master_id)
  call Broadcast_Scalar(FNLrpr,   master_id)
  call Broadcast_Scalar(FNDrpr,   master_id)
  call Broadcast_Scalar(FNIrpr,   master_id)
  call Broadcast_Scalar(FNRLrpr,  master_id)
  call Broadcast_Scalar(FNLLrpr,  master_id)
  call Broadcast_Scalar(FNDLrpr,  master_id)
  call Broadcast_Scalar(FNILrpr,  master_id)
  
  call Broadcast_Scalar(FNRrpe,   master_id)
  call Broadcast_Scalar(FNLrpe,   master_id)
  call Broadcast_Scalar(FNDrpe,   master_id)
  call Broadcast_Scalar(FNIrpe,   master_id)
  call Broadcast_Scalar(FNRLrpe,  master_id)
  call Broadcast_Scalar(FNLLrpe,  master_id)
  call Broadcast_Scalar(FNDLrpe,  master_id)
  call Broadcast_Scalar(FNILrpe,  master_id)
  
  call Broadcast_Scalar(FNRLrpd,  master_id)
  call Broadcast_Scalar(FNLLrpd,  master_id)
  call Broadcast_Scalar(FNDLrpd,  master_id)
  call Broadcast_Scalar(FNILrpd,  master_id)
  
  call Broadcast_Scalar(FRPRRPG1, master_id)
  call Broadcast_Scalar(FRPRRPG2, master_id)
  call Broadcast_Scalar(FRPRRPG3, master_id)
  
  call Broadcast_Scalar(ISRPEMSPAC, master_id)
  call Broadcast_Scalar(ISRPEMSPFR, master_id)
  call Broadcast_Scalar(ISRPEMTIME, master_id)
  call Broadcast_Scalar(ISRPEMTIFR, master_id)
  call Broadcast_Scalar(ISRPEMTILC, master_id)
  call Broadcast_Scalar(NRPEMEE,    master_id)
  
  ! *** *******************************************************************C                                                
  ! *** SET INITIAL CONDITIONS                                                                          

  do L = 1,LC
    WQRPS(L) = 0.0   ! *** Initial condition for shoot mass    ( gC/m^2 )
    WQRPR(L) = 0.0   ! *** Initial condition for root mass     ( gC/m^2 )                                                                                                           
    WQRPE(L) = 0.0   ! *** Initial condition for epiphyte mass ( gC/m^2 )                                                                                                         
    WQRPD(L) = 0.0   ! *** Initial condition for detritus mass ( gC/m^2 )                                                                                                         
    LMASKRPEM(L) = .FALSE.                                                                                                 
  enddo                                                                                                                  
  LMASKRPEM_Global(:) = .false.
  
  if( INITRPEM > 0 )then
    ! *** SET SPATIALLY VARIABLE INITIAL CONDITIONS                                                                          
    if( process_id == master_id )then
      if( INITRPEM == 1 )then
        write(*,'(A)')' WQ: WQRPEMSIC.INP'
        open(1,FILE = 'wqrpemsic.inp')
        STR = READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
      
        read(1,*) LDATA
        
        do LL = 1,LDATA                                                                                                        
          read(1,*) I, J, RPST, RPRT, RPET, RPDT                                                                                   
          LG = LIJ_Global(I,J)                                                                                                         
          WQRPS_Global(LG) = RPST     ! *** WQRPS
          WQRPR_Global(LG) = RPRT     ! *** WQRPR                                                                                                      
          WQRPE_Global(LG) = RPET     ! *** WQRPE                                                                                                      
          WQRPD_Global(LG) = RPDT     ! *** WQRPD
          ! *** Only cell defined in wqrpemsic.inp file can have mask - DKT
          if( HP_Global(LG) <= 10. ) LMASKRPEM_Global(LG) = .TRUE.   ! *** Only include the cell if not too deep
        enddo
        close(1)
        
      elseif( INITRPEM == 2 )then
        ! *** READ RESTART CONDITIONS        
        write(*,'(A)')' WQ: WQRPEMRST.INP'
        open(1,FILE = 'wqrpemrst.inp')
        
        call SKIPCOM(1,'*')  ! *** SKIP OVER TITLE AND AND HEADER LINES 
        read(1,*)TIMETMP

    100 continue                                                                                                        
        read(1,*,END = 200) LG, RPST, RPRT, RPET, RPDT                                                                             

        WQRPS_Global(LG) = RPST     ! *** WQRPS
        WQRPR_Global(LG) = RPRT     ! *** WQRPR                                                                                                      
        WQRPE_Global(LG) = RPET     ! *** WQRPE                                                                                                      
        WQRPD_Global(LG) = RPDT     ! *** WQRPD     
        LMASKRPEM_Global(LG) = .TRUE.   ! *** Only include the cell if not too deep
        GOTO 100     
        
    200 continue                                                                                                        
        close(1)
      endif
    endif   ! *** End of master_id

    call Broadcast_Array(WQRPS_Global, master_id)
    call Broadcast_Array(WQRPR_Global, master_id)
    call Broadcast_Array(WQRPE_Global, master_id)
    call Broadcast_Array(WQRPD_Global, master_id)

    ! *** Map to Local Domain
    do LG = 2,LA_GLOBAL
      L = Map2Local(LG).LL
      if( L > 1 )then
        WQRPS(L) = WQRPS_Global(LG)
        WQRPR(L) = WQRPR_Global(LG)
        WQRPE(L) = WQRPE_Global(LG)
        WQRPD(L) = WQRPD_Global(LG)
        LMASKRPEM(L) = LMASKRPEM_Global(LG)
      endif
    enddo
  
  else
    ! *** *******************************************************************C                                                
    ! *** SET SPATIALLY CONSTANT INITIAL CONDITIONS                                                                          
    do L = 1,LC
      if( HP(L) <= HRPEMIC )then  ! uniform values at depth <= HRPEMIC                                                    
        WQRPS(L) = RPSO   ! *** Initial condition for shoot mass    ( gC/m^2 )
        WQRPR(L) = RPRO   ! *** Initial condition for root mass     ( gC/m^2 )                                                 
        WQRPE(L) = RPEO   ! *** Initial condition for epiphyte mass ( gC/m^2 )
        WQRPD(L) = RPDO   ! *** Initial condition for detritus mass ( gC/m^2 )
        LMASKRPEM(L) = .TRUE.                                                                                              
      endif                                                                                                              
    enddo
    
  endif                                                                                                                  

  ! *** Count active RPEM cells
  NRPEM = 0
  do L = 1,LC
    if( LMASKRPEM(L) ) NRPEM = NRPEM + 1
  enddo

  
  ! *** ***********************************************************************************                                               
  ! *** GENERATE TEMPERATURE DEPENDENCY FOR GROWTH AND RESPIRATION OVER WQTDMIN TO WQTDMAX                                                        
  do NT = 1,NWQTD  
    WTEMP = WQTDTEMP(NT)
  
    ! *** Shoot Production (Growth)
    RPEMTPrps(NT) = 1.
    if( WTEMP < TP1RPS )then
      RPEMTPrps(NT) = EXP(-rKTP1RPS*(WTEMP-TP1RPS)*(WTEMP-TP1RPS) )
    endif
    if( WTEMP > TP2RPS )then
      RPEMTPrps(NT) = EXP(-rKTP2RPS*(WTEMP-TP2RPS)*(WTEMP-TP2RPS) )
    endif
  
    ! *** Epiphyte Production (Growth)
    RPEMTPrpe(NT) = 1.
    if( WTEMP < TP1RPE )then
      RPEMTPrpe(NT) = EXP(-rKTP1RPE*(WTEMP-TP1RPE)*(WTEMP-TP1RPE) )
    endif
    if( WTEMP > TP2RPE )then
      RPEMTPrpe(NT) = EXP(-rKTP2RPE*(WTEMP-TP2RPE)*(WTEMP-TP2RPE) )
    endif
  
    ! *** Shoot Respiration
    RPEMTRrps(NT) = 1.
    if( WTEMP < TR1RPS )then
      RPEMTRrps(NT) = EXP(-rKTR1RPS*(WTEMP-TR1RPS)*(WTEMP-TR1RPS) )
    endif
    if( WTEMP > TR2RPS )then
      RPEMTRrps(NT) = EXP(-rKTR2RPS*(WTEMP-TR2RPS)*(WTEMP-TR2RPS) )
    endif
  
    ! *** Epiphyte Respiration
    RPEMTRrpe(NT) = 1.
    if( WTEMP < TR1RPE )then
      RPEMTRrpe(NT) = EXP(-rKTR1RPE*(WTEMP-TR1RPE)*(WTEMP-TR1RPE) )
    endif
    if( WTEMP > TR2RPE )then
      RPEMTRrpe(NT) = EXP(-rKTR2RPE*(WTEMP-TR2RPE)*(WTEMP-TR2RPE) )
    endif
  
    ! *** Root Respiration
    RPEMTRrpr(NT) = 1.
    if( WTEMP < TR1RPR )then
      RPEMTRrpr(NT) = EXP(-rKTR1RPR*(WTEMP-TR1RPR)*(WTEMP-TR1RPR) )
    endif
    if( WTEMP > TR2RPR )then
      RPEMTRrpr(NT) = EXP(-rKTR2RPR*(WTEMP-TR2RPR)*(WTEMP-TR2RPR) )
    endif
  
  enddo                                                                                                                  
  
  if( process_id == master_id )then
    ! *** Write out the temperature dependencies for RPEM
    write(2,'(//,A,//,3(A11,F10.5),A10,I5)') "Look-up Table for RPEM Temperature Dependency",                               &
                                              "WQTDMIN = ", WQTDMIN, "WQTDMAX = ", WQTDMAX, "WQTDINC = ", WQTDINC, "NWQTD = ", NWQTD
    write(2,'(/,A5,A10,A15,20A10)') "IT", "TEMP", "RPEMTPrps", "RPEMTPrpe", "RPEMTRrps", "RPEMTRrpe", "RPEMTRrpr"
                                              
    do NT = 1,NWQTD
      write(2,'(I5,F10.3,F15.3,20F10.5)') NT, WQTDTEMP(NT), RPEMTPrps(NT), RPEMTPrpe(NT), RPEMTRrps(NT), RPEMTRrpe(NT), RPEMTRrpr(NT)
    enddo
  endif
  
  ! *** ASSIGN BED POREWATER CONCENTRATIONS IF SEDIMENT DIAGENESIS NOT SIMULATED
  if( IWQBEN /= 1 )then
    do L = 2,LA
      SM2NH4(L) = SMNH4
      SM2NO3(L) = SMNO3
      SM2PO4(L) = SMPO4
    enddo
  endif
  
  ! *** *******************************************************************C                                                
  return
  
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
  allocate(WQRPS(LCM),WQRPR(LCM),WQRPE(LCM),WQRPD(LCM), &
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
     
  allocate(IRPEMTS(LCM),JRPEMTS(LCM),LRPEMTS(LCM))
  allocate(LMASKRPEM(LCM))
  ! *** Allocating fluxes if sediment diagenesis is not simulated
  if( .not. allocated(SM2NH4)  )then
      allocate(SM2NH4(LCM))
      allocate(SM2NO3(LCM))
      allocate(SM2PO4(LCM))
      SM2NH4 = 0.0
      SM2NO3 = 0.0
      SM2PO4 = 0.0 
  endif
  
  WQRPS  = 0
  WQRPR  = 0
  WQRPE  = 0
  WQRPD  = 0
  PRPS  = 0
  RRPS  = 0
  RRPR  = 0
  PRPE  = 0
  RRPE  = 0
  WQRPSR  = 0
  WQRPSL  = 0
  WQRPER  = 0
  WQRPEL  = 0
  WQRPDL  = 0
  WQRPSRP  = 0
  WQRPSLP  = 0
  WQRPERP  = 0
  WQRPELP  = 0
  WQRPDLP  = 0
  WQRPSRN  = 0
  WQRPSLN  = 0
  WQRPERN  = 0
  WQRPELN  = 0
  WQRPDLN  = 0
  WQRPRR  = 0
  WQRPRL  = 0
  WQRPRRP  = 0
  WQRPRLP  = 0
  WQRPRRN  = 0
  WQRPRLN  = 0
  FRPSPW  = 0
  FRPSNW  = 0
  PNRPS  = 0
  PNRPE  = 0
  RJRPRS  = 0
  XLIMTPRPS  = 0
  XLIMTPRPE  = 0
  XLIMTRRPS  = 0
  XLIMTRRPE  = 0
  XLIMTRRPR  = 0
  XLIMNRPS  = 0
  XLIMNRPE  = 0
  XLIMLRPS  = 0
  XLIMLRPE  = 0
  RISS  = 0
  RISSO  = 0
  RISSOE  = 0
  RPEMTPRPS = 0
  RPEMTPRPE = 0
  RPEMTRRPS = 0
  RPEMTRRPE = 0
  RPEMTRRPR = 0

  JSRPEM = 1
  
END SUBROUTINE

END MODULE
