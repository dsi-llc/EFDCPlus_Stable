! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
SUBROUTINE S_TECPLOT
  
  use GLOBAL
  
  implicit none
  
  integer :: COUNTER,M
  real :: UPVEL,DNVEL,ATVEL,ZMIDLAYER
  integer :: I,J,L,LE,K,ITEMPMSK,ILL,IDOWN
  integer :: ICOUNT(TCOUNT,3)
  real :: UTMPS,VTMPS,UTMP,VTMP
  real :: MAG,MCHANGE
  real,dimension(5) :: MB
  real,dimension(TCOUNT,2) :: ESTORAGE
  real,dimension(LC) :: AVGSED
  real,dimension(KC) :: CTEMP1
  real,dimension(IC-4) :: MAGREF1,MAGREF2
  real :: RHOH2O, HTCAP, HTCONT, HTEMP              !VB HEAT CONTENT VARIABLES
  real :: WQV2, WQV3, WQV19, WQV22, WQV16,WQV17, WQV6, TIMESTEP      !VB TEMPORARY VARIABLES FOR TIMESERIES O/P
  real :: WTEMP, WVOL, VOLTEMP, WTMP, WQV10, WQV14, WQV15
  real :: VELPDF !velocity PDF variable for ocean model
  real,save :: ELAST,TLAST
  integer,save :: nstep
  logical,save :: FIRSTTIME = .FALSE.
  integer :: LL1,LL2
  real,dimension(IC,JC) :: PUPSTREAM,PDNSTREAM,PUPKIN,PDNKIN,PUPPOT,PDNPOT,DELTA
  real,dimension(IC,JC,KC) :: MDOTU,MDOTD
  real,dimension(IC) :: PUP,PDN,KUP,KDN,PPUP,PPDN
  character(LEN = 66) :: timeline
  
  DOUBLE PRECISION,save,allocatable,dimension(:)     :: CBLTOT        !(LCM)
  DOUBLE PRECISION,allocatable,dimension(:)     :: THCK          !(LCM)
   DOUBLE PRECISION,allocatable,dimension(:,:)   :: UVEL          !(LCM,KC)
  DOUBLE PRECISION,allocatable,dimension(:,:)   :: VVEL          !(LCM,KC)
   allocate(CBLTOT(LC))
    allocate(THCK(LC))
    allocate(UVEL(LC,KC))
    allocate(VVEL(LC,KC))
  CBLTOT = 0.0  ! tecplot only - move
  THCK = 0.0        !(LCM)
    UVEL = 0.0        !(LCM,KC)
  VVEL = 0.0        !(LCM,KC)

  
  timeline = 'TEXT X = 9, Y = 84, T = "&(ZONENAME[0001])", F = TIMES, CS = FRAME, H = 3, ZN = '
  if( .not. FIRSTTIME )then
! This opens the Tecplot output file
    if( MAXVAL(MVEGL(2:LA))>90 )then !MHK devices exist
      open(UNIT = 222,FILE = OUTDIR//'powerout.dat')
      ITEMPMSK = 1
      do ILL = 2,LA
        if( MVEGL(ILL)>90 )then
          ICOUNT(ITEMPMSK,1) = ITEMPMSK
          ICOUNT(ITEMPMSK,2) = IL(ILL)
          ICOUNT(ITEMPMSK,3) = JL(ILL)
          ITEMPMSK = ITEMPMSK+1
        endif
      enddo
      write(222,'("TURBINE",100(I6,6X))')(ICOUNT(ILL,1),ILL = 1,TCOUNT)
      write(222,'(7X,100(3X,I3,3X,I3))') (ICOUNT(ILL,2),ICOUNT(ILL,3),ILL = 1,TCOUNT)
    endif
    OPEN (UNIT = 111,FILE = OUTDIR//'tecplot2d.dat')
    if( ISTRAN(8) == 1 )then !WQ data
      write(111,'(A36)')'TITLE = "EFDC 2D Tecplot Algae Data"'
      write(111,'(A77)')'VARIABLES = "X","Y","U","V","SPEED","CHG(g)","P4D","NHX","NOX","DO(g)","CO2(g)"' !,"HEAT(kJ)","VOL(M3)"'
    elseif( LSEDZLJ )then !SEDZLJ data
      write(111,*)'TITLE = "EFDC 2D Tecplot Sediment Data"'
      write(111,*)'VARIABLES = "I","J","X","Y","TAU","D50","CBL","SED","U","V","THICK1","SPEED"'
      OPEN (UNIT = 112,FILE = 'massbal.dat')
      write(112,*)'TITLE = "EFDC Mass Balance Data"'
      write(112,*)'VARIABLES = "Time","MB1","MB2","MB3","MB4","MB5","ERATE","D50"'
    else !Flow data
      write(111,'(A36)')'TITLE = "EFDC 2D Tecplot Flow Data"'
      write(111,'(A77)')'VARIABLES = "X","Y","U","V","SPEED","SHEAR","DEPTH"' !,"HEAT(kJ)","VOL(M3)"'
    endif
    if( OUTPUTFLAG == 4 )then   
      open(444,FILE = OUTDIR//'CALIBRATION.DAT')
      write(444,'("TIME        VEL-FED DEFECT  UPVEL   ATVEL   DNVEL")')
    elseif( OUTPUTFLAG == 3 )then
      open(864,FILE = OUTDIR//'TIDALREF.DAT')
      write(864,'("3 ROWS OF AVERAGE VELOCITY, TOP VELOCITY, AND DEPTH")')
    elseif( OUTPUTFLAG == 2 )then
      open(468,FILE = OUTDIR//'VELPDF.DAT')
      write(468,'("TIME AVERAGE VELOCITY ACROSS NARROWS I = 60?")')
      open(579,FILE = OUTDIR//'ZPROF.DAT')
      write(579,'("TIME VELOCITIES FOR EACH LAYER FOR PROFILES")')
    elseif( OUTPUTFLAG == 1 )then
      open(765,FILE = OUTDIR//'POWERDIF.DAT')
      write(765,'("TIME   POWER DIFFERENCES ACROSS DIFFERENT I-CELL COLUMNS")')
      open(567,FILE = OUTDIR//'POTENTIAL.DAT')
      write(567,'("TIME   POTENTIAL DIFFERENCES ACROSS DIFFERENT I-CELL COLUMNS")')
      open(678,FILE = OUTDIR//'KINETIC.DAT')
      write(678,'("TIME   KINETIC DIFFERENCES ACROSS DIFFERENT I-CELL COLUMNS")')
    endif
    FIRSTTIME = .TRUE.
  endif
  TIMESTEP = TIMEDAY
!  write(765,*)TIMESTEP,SUM(WQV(3,1:KC,3)*DZC(L,1:KC))
  ITEMPMSK = 1
  do ILL = 2,LA
    if( MVEGL(ILL)>90 )then
      ESTORAGE(ITEMPMSK,1) = SUM(ESUP(:,ILL))
      ESTORAGE(ITEMPMSK,2) = SUM(EMHK(:,ILL))
      ITEMPMSK = ITEMPMSK+1
    endif
  enddo  
  if( MAXVAL(MVEGL(2:LA))>90)WRITE(222,'(F7.3,100(F7.4,1X))')TIMESTEP,(ESTORAGE(ILL,1),ESTORAGE(ILL,2),ILL = 1,TCOUNT)
  nstep = nstep+1
  if( nstep>9999)PRINT*,'Tecplot timestamp is greater than 9999, increase field width or reduce writing frequency'                
  write(timeline(31:34),'(I4.4)')nstep
  write(111,'(A66,I4.4)')timeline,nstep
  write(111,*)'ZONE T = "',TIMESTEP,'" I= ' ,IC-4,' J= ' ,JC-4,' F = POINT'
  if( ISTRAN(8) == 1 )then

    !VB HEAT CONTENT VARIABLE
    RHOH2O = 1000     !KG/M3
    HTCAP = 4.187     !kJ/KG-K
    HTCONT = 0.0
    ! 2 Dimensional Output !VB REWRITTEN TO OUTPUT TIME SERIES PLOTS
    do J = 3,JC-2
      do I = 3,IC-2
        if( LIJ(I,J) > 0 )then
!VB       TWATER = SUM(TEM(LIJ(I,J),1:KC)*DZC(L,1:KC))
!        TALT = MAXVAL(TWQ(LIJ(I,J))
          L = LIJ(I,J)
        UTMPS = SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
        VTMPS = SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
        MAG = SQRT((UTMPS)**2+(VTMPS)**2)
        WTEMP = SUM((TEM(L,1:KC)+273.15)*DZC(L,1:KC)*DXYP(L)*HP(L))
        WTMP = WTMP+WTEMP
        HTEMP = RHOH2O*HTCAP*WTEMP
        HTCONT = HTCONT+HTEMP            !HEAT CONTENT IN KILOJOULES                    
        WQV2 = SUM(WQV(L,1:KC,2)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WQV3 = SUM(WQV(L,1:KC,3)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WQV6 = SUM(WQV(L,1:KC,6)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV10 = SUM(WQV(L,1:KC,10)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV14 = SUM(WQV(L,1:KC,14)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV15 = SUM(WQV(L,1:KC,15)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV16 = SUM(WQV(L,1:KC,16)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV17 = SUM(WQV(L,1:KC,17)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV19 = SUM(WQV(L,1:KC,19)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WQV22 = SUM(WQV(L,1:KC,22)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        write(111,'(5(1X,F11.2), 7(1X,F11.4))')DLON(L),DLAT(L),UTMPS,VTMPS,MAG,WQV3,WQV10,WQV14,WQV15,WQV19,WQV22
      endif
      enddo
    enddo
    do L = 2, LA
      VOLTEMP = SUM(DZC(L,1:KC)*DXYP(L)*HP(L))
      WVOL = WVOL + VOLTEMP
    enddo
  elseif( LSEDZLJ )then      
! OUTPUT FOR TECPLOT, SEDIMENT DATA
    MB(1) = 0.0 !Mass in bedload (kg)
    MB(2) = 0.0 !Mass in suspension (kg)
    !MB(3)     !Mass exiting in suspension (stored)
    !MB(4)     !Mass exiting as bedload (stored)
    MB(5) = 0.0 !Mass eroded from bed (kg)
    FORALL(L = 2:LA) !added to convert THCK from g/cm^2 to cm
      TSEDT(L) = SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
    ENDFORALL
    FORALL(L = 2:LA) THCK(L) = TSEDT(L) - TSET0T(L)
    do J = 3,JC-2
      do I = 3,IC-2
        L = LIJ(I,J)
        CBLTOT(L) = 10.0*SUM(CBL(L,1:NSEDS))*DXYP(L) !g/cm^3*cm*m^2*(0.001*100*100)
        CTEMP1(1:KC) = 0.001*SUM(SED(L,1:KC,1:NSEDS))
        MB(1) = MB(1) + CBLTOT(L)
        do K = 1,KC
          MB(2) = MB(2) + CTEMP1(K)*DZC(L,K)*DXYP(L)*HP(L) !g/m^3*m^2*m*(0.001)
          UVEL(L,K) = U(L,K)*CUE(L) + V(L,K)*CVE(L)
          VVEL(L,K) = U(L,K)*CUN(L) + V(L,K)*CVN(L)
          if( I == IC-3 )then
            !MB(3) = MB(3) - DT*VVEL(LIJ(I,J),K)*DYP(LIJ(I,J))*DZC(L,K)*HP(LIJ(I,J))*SUM(SED(LIJ(I,J),K,1:NSEDS))*0.001 !Dt*m/s*m*m*g/m^3*(0.001)
            MB(3) = MB(3) - 0.25*DT*VVEL(L,K)*DYP(L)*DZC(L,K)*(HP(L) + HP(LIJ(I+1,J)))*(SUM(SED(L,K,1:NSEDS)) + SUM(SED(LIJ(I+1,J),K,1:NSEDS)))*0.001 !Dt*m/s*m*m*g/m^3*(0.001)
          endif
        enddo
        if( I == IC-3 )then
          !MB(4) = MB(4) - 0.05*DT*DYP(L)*SUM((UBL(L,1:NSEDS)*CUN(L) + VBL(L,1:NSEDS)*CVN(L))*(CBL(L,1:NSEDS) + CBL(LIJ(I+1,J),1:NSEDS),1:NSEDS)) !Dt*m*cm/s*g/cm^3*cm*(0.001*100)
        endif
        MB(5) = MB(5) - 10.0*THCK(L)*DXYP(L) !g/cm^2*m^2*(0.001*100*100)
        AVGSED(L) = SUM(CTEMP1(1:KC)*DZC(L,1:KC))*DXYP(L)*HP(L)
        !WRITE(111,'(I4,1X,I4,1X,10(E13.4,1X))')I, J, DLON(LIJ(I,J)), DLAT(LIJ(I,J)), TAU(LIJ(I,J)), D50AVG(LIJ(I,J)), CBLTOT(LIJ(I,J)), AVGSED(LIJ(I,J)), SUM((U(LIJ(I,J),1:KC)*CUE(LIJ(I,J)) + V(LIJ(I,J),1:KC)*CVE(LIJ(I,J)))*DZC(L,1:KC)), SUM((U(LIJ(I,J),1:KC)*CUN(LIJ(I,J)) + V(LIJ(I,J),1:KC)*CVN(LIJ(I,J)))*DZC(L,1:KC)), SUM(TSED(1:KB,LIJ(I,J)))/1.4625, THCK(LIJ(I,J))/1.4625
        UTMPS = SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
        VTMPS = SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
        MAG = SQRT((UTMPS)**2+(VTMPS)**2)
        write(111,'(I4,1X,I4,1X,10(E13.4,1X))')I, J, DLON(L), DLAT(L), TAU(L), D50AVG(L), CBLTOT(L), AVGSED(L), SUM((U(L,1:KC)*CUE(L) + V(L,1:KC)*CVE(L))*DZC(L,1:KC)), SUM((U(L,1:KC)*CUN(L) + V(L,1:KC)*CVN(L))*DZC(L,1:KC)), THCK(L), MAG
      enddo
    enddo
    !PRINT*,TBEGIN + FLOAT(N-1)*DT,MB(1),MB(2),MB(3),MB(4),MB(5),SUM(TAU(2:LA))/FLOAT(LA-1),MAXVAL(TAU(2:LA))
    !WRITE(112,'(8(E13.4,1X))')TBEGIN + FLOAT(N-1)*DT,MB(1),MB(2),MB(3),MB(4),MB(5),SUM(TAU(2:LA))/FLOAT(LA-1),MAXVAL(TAU(2:LA))
    write(112,'(8(E13.4,1X))')TIMESTEP,MB(1),MB(2),MB(3),MB(4),MB(5),(MB(5)-(MB(1) + MB(2))-ELAST)/((TBEGIN + FLOAT(N-1)*DT)-TLAST),SUM(D50AVG(:))/(FLOAT(LA-1))
    ELAST = MB(5)-(MB(1) + MB(2)) !Erosion is saved (as mass eroded, minus the mass in the water column)
    TLAST = TBEGIN + FLOAT(N-1)*DT
    !print *,SNGL(CBL(LIJ(13,4),1:7))
    !print *,SNGL(SED(LIJ(13,4),KC,1:7))
    !print *, maxval(tau(2:LA)), minval(tau(2:LA))
  else
    do J = 3,JC-2
      do I = 3,IC-2
        if( LIJ(I,J) > 0 )then
          L = LIJ(I,J)
          UTMPS = SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
          VTMPS = SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
          MAG = SQRT((UTMPS)**2+(VTMPS)**2)
          write(111,'(5(1X,F11.2), 7(1X,F11.4))')DLON(L),DLAT(L),UTMPS,VTMPS,MAG,TAUB(L),HP(L)
        endif
      enddo
    enddo
  endif	
  if( OUTPUTFLAG == 1 )then !MHK look at power up- and down-stream of the MHK device
    PUPSTREAM = 0.0;PDNSTREAM = 0.0;PUPPOT = 0.0;PDNPOT = 0.0;PUPKIN = 0.0;PDNKIN = 0.0;PUP = 0.0;PDN = 0.0;KUP = 0.0;KDN = 0.0;PPUP = 0.0;PPDN = 0.0;MDOTU = 0.0;MDOTD = 0.0
    do I = 1,1
      do J = 3,JC-2
        LL1 = LIJ(5-I,J) !looking where whatever integer is in this expression and the next
        LL2 = LIJ(I+6,J)
        do K = 1,KC  
          MDOTU(I,J,K) = 1024.0*U(LL1,K)*DYU(LL1)*DZC(L,K)*HU(LL1)
          MDOTD(I,J,K) = 1024.0*U(LL2,K)*DYU(LL2)*DZC(L,K)*HU(LL2)
          PUPKIN(I,J) = PUPKIN(I,J)+0.5*MDOTU(I,J,K)*U(LL1,K)**2
          PDNKIN(I,J) = PDNKIN(I,J)+0.5*MDOTD(I,J,K)*U(LL2,K)**2
          PUPPOT(I,J) = PUPPOT(I,J)+MDOTU(I,J,K)*(FLOAT(K)-0.5)/DZI*G*HU(LL1)*DZC(L,K)
          PDNPOT(I,J) = PDNPOT(I,J)+MDOTD(I,J,K)*(FLOAT(K)-0.5)/DZI*G*HU(LL2)*DZC(L,K)
          continue
        enddo
        DELTA(I,J) = SUM(MDOTU(I,J,1:KC))-SUM(MDOTD(I,J,1:KC))
 !      EUPSTREAM(I,J) = 0.5*DXP(LL1)*HP(LL1)*1024.*SUM((U(LL1,1:KC)+U(LL1+1,1:KC))*DZC(L,1:KC)*(0.5*(0.5*(U(LL1,1:KC)+U(LL1+1,1:KC)))**2+9.8106*HP(LL1)))
  !     EDNSTREAM(I,J) = 0.5*DXP(LL2)*HP(LL2)*1024.*SUM((U(LL2,1:KC)+U(LL2+1,1:KC))*DZC(L,1:KC)*(0.5*(0.5*(U(LL2,1:KC)+U(LL2+1,1:KC)))**2+9.8106*HP(LL2)))
        PUPSTREAM(I,J) = PUPPOT(I,J)+PUPKIN(I,J)
        PDNSTREAM(I,J) = PDNPOT(I,J)+PDNKIN(I,J)
      enddo
      MCHANGE = SUM(MDOTU(I,3:JC-2,1:KC))-SUM(MDOTD(I,3:JC-2,1:KC))
      PUP(I) = SUM(PUPPOT(I,3:JC-2))
      PDN(I) = SUM(PDNPOT(I,3:JC-2))
      KUP(I) = SUM(PUPKIN(I,3:JC-2))
      KDN(I) = SUM(PDNKIN(I,3:JC-2))
      PPUP(I) = SUM(PUPSTREAM(I,3:JC-2))
      PPDN(I) = SUM(PDNSTREAM(I,3:JC-2))
      continue
    enddo
    write(765,'(E10.3,8(1X,F10.1))')TIMESTEP,(SUM(PUPSTREAM(I,3:JC-2))-SUM(PDNSTREAM(I,3:JC-2)),I = 2,8)
    write(567,'(E10.3,16(1X,F10.0))')TIMESTEP,(SUM(PUPPOT(I,3:JC-2)),I = 2,8),(SUM(PDNPOT(I,3:JC-2)),I = 2,8)
    write(678,'(E10.3,16(1X,F10.1))')TIMESTEP,(SUM(PUPKIN(I,3:JC-2)),I = 2,8),(SUM(PDNKIN(I,3:JC-2)),I = 2,8)
  endif
 !OUTPUTFLAG = 2 !Tidal reference model average and z-profile velocities
  if( OUTPUTFLAG == 2 )then
    I = 60
    VELPDF = 0.0
    do J = 3,JC-2
      L = LIJ(I,J)
      UTMPS = SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
      VTMPS = SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
      MAG = SQRT((UTMPS)**2+(VTMPS)**2)
      VELPDF = VELPDF+MAG
    enddo
    VELPDF = VELPDF/FLOAT(JC-4)
    write(468,*)TIMESTEP,VELPDF
    write(579,'(E10.3,10(1X,F6.3))')TIMESTEP,((U(L,K)*CUE(L)+V(L,K)*CVE(L)),K = 1,KC)
  endif
 if( OUTPUTFLAG == 3 )then !average and surface velocities for tidal reference model
   J = 20
   do I = 3,IC-2
     L = LIJ(I,J)
     LE = LIJ(I+1,J)
     UTMPS = SUM(0.5*((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))+0.5*((U(LE,1:KC)*CUE(LE)+V(LE,1:KC)*CVE(LE))*DZC(L,1:KC)))
     VTMPS = SUM(0.5*((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))+0.5*((U(LE,1:KC)*CUN(LE)+V(LE,1:KC)*CVN(LE))*DZC(L,1:KC)))
     MAGREF1(I-2) = VTMPS
     UTMP=     0.5*(U(L,KC) *CUE(L) +V(L,KC) *CVE(L))
     UTMP = UTMP+0.5*(U(LE,KC)*CUE(LE)+V(LE,KC)*CVE(LE))
     UTMPS = UTMP
     VTMP=     0.5*(U(L,KC) *CUN(L)+ V(L,KC)* CVN(L))
     VTMP = VTMP+0.5*(U(LE,KC)*CUN(LE)+V(LE,KC)*CVN(LE))
     VTMPS = VTMP
     MAGREF2(I-2) = VTMPS
   enddo
   write(864,'(84(F6.3,1X))')(MAGREF1(I-2),I = 3,IC-2)
   write(864,'(84(F6.3,1X))')(MAGREF2(I-2),I = 3,IC-2)
   write(864,'(84(F6.3,1X))')(HP(I-2),I = 3,IC-2)
 ENDIF
 if( OUTPUTFLAG == 4 )then !this interrogates the W-2-E straight flow channel for wake calibration
   UPVEL = 0.0;DNVEL = 0.0;ATVEL = 0.0;COUNTER = 0
   do L = 2,LA
     if( MVEGL(L)>90 )then
       M = MVEGL(L)-90
       do K = 1,KC
         ZMIDLAYER = HP(L)*(SUM(DZC(L,1:K))-0.5*DZC(L,K))+BELV(L) !midlayer height
         if( ZMIDLAYER >= ZMINMHK(M,L) .and. ZMIDLAYER <= ZMAXMHK(M,L) )then !MHK device exists in this layer
           COUNTER = COUNTER+1 
           IDOWN = NINT(20.*WIDTHMHK(M)/DXP(L))
           ATVEL = ATVEL+U(LEC(L),K)    !1 cell downstream
           UPVEL = UPVEL+U(L-IDOWN/4,K) !20 cells upstream
           DNVEL = DNVEL+U(L+IDOWN,K)   !80 cells downstream for a 4-m turbine
         endif
       enddo
     endif
   enddo
   ATVEL = ATVEL/FLOAT(COUNTER)
   UPVEL = UPVEL/FLOAT(COUNTER)
   DNVEL = DNVEL/FLOAT(COUNTER)
   write(444,'(e10.2,1x,5(f7.4,1x))')TIMESTEP,1.0-DNVEL/UPVEL,1.0-ATVEL/UPVEL,UPVEL,ATVEL,DNVEL !calculate velocity ratio, we are looking for 90% recovery at 20D downstream
 ENDIF
 RETURN
END SUBROUTINE S_TECPLOT
