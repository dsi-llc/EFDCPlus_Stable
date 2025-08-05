! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE GETSWANMOD
  ! *** PURPOSE:
  ! *** DIRECTLY READ SWAN OUTPUT CONTAINED IN
  ! *** FRM/GRP FILE OR TBL FILE CORRESPONDING TO LOC FILE
  ! *** SWN FILE: SET INRHOG 1 TO GET DISSIPATION IN W/M2
  ! *** IFWAVE  = 1   :C14A
  ! *** SWANGRP = 1/0 :C14A
  ! *** OTHER PARAM :C14B
  ! *** AUTHOR: DANG HUU CHUNG

  use GLOBAL 
  use INFOMOD
  use WAVELENGTH
  use XYIJCONV
  use WINDWAVE,only:SHLIM,WHMI
  Use Variables_MPI_Write_Out
  use Variables_MPI
  use Broadcast_Routines
  
  implicit none

  contains

  SUBROUTINE GETSWAN_GRP
  ! *** READ OUTPUT BY SWAN FOR THE WHOLE EFDC GRID
  ! *** AND THE TERMS IN SWN FILE DETERMINED BY:
  ! $************ OUTPUT REQUESTS *************************
  ! *** GROUP  'GRP1' SUBGrid 0 IC-1  0 JC-1    
  ! *** TABLE  'GRP1' HEADer 'EXAMPLE.grp' HS PDIR RTP WLEN DEPTH DISSIP
  !
  ! *** OUTPUT:
  ! *** GENERATE WAVE parameterS FOR EVERY CELL L = 2:LA OF EFDC

  character(200) :: STR, TERM
  integer :: NL, I, J, L, ISO, NW, IPOS, LG
  
  integer(4),save :: SWAN_M1, SWAN_N1, SWAN_M2, SWAN_N2
  integer(4),save :: NHS, NPK, NRT, NWL, NHP, NDI, NCOL
  real(RKD),save  :: VAL(20)
  
  real(RKD)   :: WDIR, WPRD, EDIS
  real(RKD)   :: WVDX, WVDY, WVCX, WVCY
  if( .not. allocated(LWVCELL_Global) )then
    allocate(LWVCELL_Global(LCM_Global))
    allocate(LWVMASK_Global(LCM_Global))
    LWVMASK_Global = .false.
  endif
  
  if( JSWAVE == 0 )then
    ! *** First call
    ! *** Identify M1,N1,M2,N2 from SWAN.EE - needs to use SCAN intrinsic function
    if( process_id == master_id )then 
      open(99,FILE = 'SWAN.EE',ACTION = 'READ')
      read(99,*) STR
      read(99,*) STR
      IPOS = SCAN(STR, "=", BACK = .TRUE.)
      read(STR(1+IPOS:),*) SWAN_M1
      read(99,*) STR
      IPOS = SCAN(STR, "=", BACK = .TRUE.)
      read(STR(1+IPOS:),*) SWAN_N1
      read(99,*) STR
      IPOS = SCAN(STR, "=", BACK = .TRUE.)
      read(STR(1+IPOS:),*) SWAN_M2
      read(99,*) STR
      IPOS = SCAN(STR, "=", BACK = .TRUE.)
      read(STR(1+IPOS:),*) SWAN_N2
    endif
    call Broadcast_Scalar(SWAN_M1, master_id)
    call Broadcast_Scalar(SWAN_N1, master_id)
    call Broadcast_Scalar(SWAN_M2, master_id)
    call Broadcast_Scalar(SWAN_N2, master_id)
    
    if( process_id == master_id )then
      write(*,'(A)')'READING SWAN_GRP.INP'
      open(UGRP,FILE = 'swan_grp.inp',ACTION = 'READ')  !FRMFILE
      do NL = 1,5
        read(UGRP,'(A)',IOSTAT = ISO) STR
        if( ISO > 0 ) STOP 'SWAN_GRP.INP: READING ERROR'  
      enddo
      STR  = ADJUSTL(STR)
      TERM = ADJUSTL(STR(2:))
      STR = READSTR(UGRP)
      NCOL = NUMCOL(STR)
      
      NHS = FINDSTR(TERM,'Hsi',NCOL)   ! *** INCIDENT OFFSHORE SIGNIFICANT WAVE HEIGHT (m)
      NPK = FINDSTR(TERM,'PkD',NCOL)   ! *** WAVE ANGLE (DEG)
      NRT = FINDSTR(TERM,'RTp',NCOL)   ! *** TIME TO PEAK (SEC)
      NWL = FINDSTR(TERM,'Wle',NCOL)   ! *** WAVE LENGTH (M)
      NHP = FINDSTR(TERM,'Dep',NCOL)   ! *** WATER DEPTH (M)
      NDI = FINDSTR(TERM,'Dis',NCOL)   ! *** WAVE DISSIPATION (M)
      
      ! *** READ BUFFER DATA 
      do NW = 1,IWVCOUNT-1
        do J = 1,JC_Global
          do I = 1,IC_Global
            read(UGRP,*,IOSTAT = ISO) VAL(1:NCOL)
            if( ISO > 0 ) STOP 'SWAN_GRP.INP: READING ERROR' 
          enddo
        enddo
      enddo
      
      NWVCELLS = 0
      LOOPJ: DO J = SWAN_N1+1, SWAN_N2+1
        do I = SWAN_M1+1, SWAN_M2+1 
          read(UGRP,*,IOSTAT = ISO) VAL(1:NCOL)   
          if( ISO > 0 ) STOP 'SWAN_GRP.INP: READING ERROR'  
      
          LG = LIJ_GLOBAL(I,J)    
          if( LG <= 1 .or. LG > LA_Global ) CYCLE
        
          NWVCELLS = NWVCELLS+1
          LWVCELL_Global(NWVCELLS) = LG   ! *** PMC - 2015-12 - This was pointing to L-1 but LWVCELL is used as an L array reference.  This has to be L not L-1.
          LWVMASK_Global(LG) = .TRUE.
          
          ! *** WAVE HEIGHT (M), THEN ZERO WAVE IMPACTS FOR VEGETATION
          WV_Global(LG).HEISIG = VAL(NHS)
          if( MVEG_Global(LG) /= MVEGOW  )then
            WV_Global(LG).HEISIG = 0.
          else
            if( (MVEG_Global(LWC_Global(LG)) /= MVEGOW) .or. &
                (MVEG_Global(LEC_Global(LG)) /= MVEGOW) .or. &
                (MVEG_Global(LSC_Global(LG)) /= MVEGOW) .or. &
                (MVEG_Global(LNC_Global(LG)) /= MVEGOW)      )  WV_Global(LG).HEISIG = 0.5*WV_Global(LG).HEISIG
          endif
          
          WDIR    = VAL(NPK)                                ! *** WAVE ANGLE = (TRUE EAST,WAVE) ANTI-CLOCKWISE [0,360] [D]
          WPRD    = VAL(NRT)                                ! *** WAVE PERIOD [S] 
          WV_Global(LG).LENGTH  = VAL(NWL)                          ! *** WAVE LENGTH [M] 
          EDIS    = VAL(NDI)                                ! *** WAVE DISSIPATION  [ W/m2 = Kg/s3]
          !WV(L).DISSIPA(KC) = EDIS/RHO                     ! *** [M3/S3]   
          
          if( WV_Global(LG).HEISIG >= WHMI )then                       
            WVDX= COS(WDIR*PI/180)
            WVDY= SIN(WDIR*PI/180)
            WVCX =  CVN_Global(LG)*WVDX - CVE_Global(LG)*WVDY               ! *** TO LOCAL CURVI-LINEAR SYSTEM
            WVCY = -CUN_Global(LG)*WVDX + CUE_Global(LG)*WVDY
            WV_Global(LG).DIR  = ATAN2(WVCY,WVCX)                   ! *** CELL-EAST,WAVE (COUNTER-CLOCKWISE [-pi,pi])
            WV_Global(LG).FREQ = 2.*PI/WPRD
            !WV_Global(LG).DISSIPA(KC) = 0.25*g*WV_Global(LG).HEISIG**2/WPRD ! *** Energy dissipation due to breaking wave  
          else
            WV_Global(LG).DIR  = 0.
            WV_Global(LG).LENGTH     = 0.
            WV_Global(LG).FREQ  = 1.
            WV_Global(LG).HEISIG   = 0.
            WV_Global(LG).DISSIPA(KC) = 0
          endif
        enddo
      enddo LOOPJ
      WV_Global(2:LA_Global).HEIGHT = WV_Global(2:LA_Global).HEISIG  ! *** ASSIGN INCIDENT SIGNIFICANT WAVE HEIGHT TO APPLIED WAVE HEIGHT
    endif
  endif
  ! *** Broadcast wave cells properties
  Call Broadcast_Array(LWVCELL_Global,   master_id)
  Call Broadcast_Array(LWVMASK_Global,   master_id)
  Call Broadcast_Array(WV_Global.HEIGHT, master_id)
  Call Broadcast_Array(WV_Global.DIR,    master_id)
  Call Broadcast_Array(WV_Global.LENGTH, master_id)
  Call Broadcast_Array(WV_Global.FREQ,   master_id)
  Call Broadcast_Array(WV_Global.HEISIG, master_id)
  
  ! *** Map wave cells properties to each process
  NWVCELLS = 0
  do LG = 2, LA_Global
    ! *** Get local L
    L = Map2Local(LG).LL
    ! *** Valid if greater than zero
    if(L > 0 )then
      ! *** set local copy of this
      NWVCELLS = NWVCELLS + 1

      LWVCELL(NWVCELLS) = L
      LWVMASK(L)   = LWVMASK_Global(LG)
      WV(L).DIR    = WV_Global(LG).DIR
      WV(L).LENGTH = WV_Global(LG).LENGTH
      WV(L).FREQ   = WV_Global(LG).FREQ  
      WV(L).HEIGHT = WV_Global(LG).HEIGHT
      WV(L).HEISIG = WV_Global(LG).HEISIG
      
      ! *** Computing Energy dissipation
      WV(L).DISSIPA(KC) = 0.25*g*WV(L).HEIGHT**2*WV(L).FREQ/(2.*PI)
    endif
  enddo
  
  if( IWVCOUNT == NWVTIM )then
    close(UGRP)
    write(*,'(A)') '** WAVE: END OF WAVE RECORD'
  endif 
  END SUBROUTINE

  SUBROUTINE GETSWAN_LOC
  ! *** READ OUTPUT BY SWAN BASED ON A LOC FILE
  ! *** AND THE TERMS IN SWN FILE DETERMINED BY:
  ! 
  ! *** POINTS 'loc' FILE   'EXAMPLE.loc'
  ! *** TABLE  'loc' HEAD   'EXAMPLE.tbl'  HS PDIR RTP WLEN DEPTH DISSIP 
  !
  ! *** OUTPUT:
  ! *** WAVE parameterS FOR THE CELLS OF EFDC
  
  character(200) :: STR,TERM
  integer :: N,L,ISO,NP,NW
  
  integer(4),save :: SWAN_M1, SWAN_N1, SWAN_M2, SWAN_N2
  integer(4),save :: NHS, NPK, NRT, NWL, NHP, NDI, NCOL
  real(RKD),save  :: VAL(20)
  
  real(RKD)   :: WDIR, WPRD, EDIS
  real(RKD)   :: WVDX, WVDY, WVCX, WVCY

  if( JSWAVE == 0 )then
    
    if( process_id == master_id )then
        !** FIRST CALL
        write(*,'(A)')'READING SWAN_LOC.INP'
        open(ULOC,FILE = 'swan_loc.inp',ACTION = 'READ')
        STR = READSTR(ULOC)  
        NP = 0
        do while(.TRUE.)
          read(ULOC,'(A)',end = 250) STR
          NP = NP+1
        enddo
        250 REWIND(ULOC) 
        NLOC = NP
    endif
    
    allocate(SWNLOC%ICEL(NLOC),SWNLOC%JCEL(NLOC))
    allocate(SWNLOC%XCEL(NLOC),SWNLOC%YCEL(NLOC))

    if( process_id == master_id )then
        STR = READSTR(ULOC)  
        do N = 1,NLOC
          read(ULOC,*) SWNLOC%XCEL(N),SWNLOC%YCEL(N)
        enddo
        close(ULOC)
        call XY2IJ(SWNLOC)  

        write(*,'(A)')'READING SWAN_TBL.INP'
        open(UTBL,FILE = 'swan_tbl.inp',ACTION = 'READ')
        do N = 1,5
          read(UTBL,'(A)',IOSTAT = ISO) STR
          if( ISO > 0 )then
            STOP 'SWAN_TBL.INP: READING ERROR'  
          elseif( ISO < 0 )then
            PRINT*,'END OF SWAN_TBL.INP'
          endif
        enddo
        STR  = ADJUSTL(STR)
        TERM = ADJUSTL(STR(2:))
        STR = READSTR(UTBL)   !NUMBER
        NCOL = NUMCOL(STR)
        
        !%       Hsig          PkDir         RTpeak        Wlen          Depth         Dissip        Genera        Redist        Radstr   
        !%       [m]           [degr]        [sec]         [m]           [m]           [m2/s]        [m2/s]        [m2/s]        [m2/s]   
        
        NHS = FINDSTR(TERM,'Hsi',NCOL)   ! *** INCIDENT OFFSHORE SIGNIFICANT WAVE HEIGHT (m)
        NPK = FINDSTR(TERM,'PkD',NCOL)   ! *** WAVE ANGLE (DEG)
        NRT = FINDSTR(TERM,'RTp',NCOL)   ! *** TIME TO PEAK (SEC)
        NWL = FINDSTR(TERM,'Wle',NCOL)   ! *** WAVE LENGTH (M)
        NHP = FINDSTR(TERM,'Dep',NCOL)   ! *** WATER DEPTH (M)
        NDI = FINDSTR(TERM,'Dis',NCOL)   ! *** WAVE DISSIPATION (M)
    
        ! *** READ BUFFER DATA 
        do NW = 1,IWVCOUNT-1
          do N = 1,NLOC
            read(UTBL,*,IOSTAT = ISO) VAL(1:NCOL)
            if( ISO > 0 ) STOP 'SWAN_TBL.INP: READING ERROR' 
          enddo
        enddo
    endif
    
    !Call Broadcast_Array(SWNLOC%ICEL, master_id)
    !Call Broadcast_Array(SWNLOC%JCEL, master_id)
    !Call Broadcast_Array(SWNLOC%XCEL, master_id)
    !Call Broadcast_Array(SWNLOC%YCEL, master_id)

  endif

  NWVCELLS = 0
  do N = 1,NLOC
    read(UTBL,*,IOSTAT = ISO) VAL(1:NCOL) 
    if( ISO > 0 ) STOP 'SWAN_TBL.INP: READING ERROR'  
  
    L= LIJ(SWNLOC%ICEL(N),SWNLOC%JCEL(N))
    if( L <= 1 .or. L>LA) CYCLE
 
    NWVCELLS = NWVCELLS+1
    LWVCELL(NWVCELLS) = L   ! *** PMC - 2015-12 - This was pointing to L-1 but LWVCELL is used as an L array reference.  This has to be L not L-1.
    LWVMASK(L) = .TRUE.
    WV(L).HEISIG = VAL(NHS)                          ! *** SIGNIFICANT WAVE HEIGHT (M)
    
    ! *** VEGETATION EFFECT
    if( ISVEG > 0 )then
      if( MVEGL(L) /= MVEGOW  )then
        WV(L).HEIGHT = 0.
      else
        if( (MVEGL(LWC(L)) /= MVEGOW) .or. &
            (MVEGL(LEC(L)) /= MVEGOW) .or. &
            (MVEGL(LSC(L)) /= MVEGOW) .or. &
            (MVEGL(LNC(L)) /= MVEGOW)      )  WV(L).HEIGHT = 0.5*WV(L).HEIGHT
      endif
    endif

    WDIR    = VAL(NPK)                          ![D] WAVE ANGLE = (TRUE EAST,WAVE) ANTI-CLOCKWISE [0,360]
    WPRD    = VAL(NRT)                          ![S] PERIOD
    WV(L).LENGTH  = VAL(NWL)                    ![M] W LENGTH
    EDIS    = VAL(NDI)                          ! W/m2 = Kg/s3: WAVE DISSIPATION
    !WV(L).DISSIPA(KC) = EDIS/RHO               ![M3/S3] 
    if( WV(L).HEISIG >= WHMI .and. WVPRD > 0 )then                    
      WVDX= COS(WDIR*PI/180)
      WVDY= SIN(WDIR*PI/180)
      WVCX =  CVN(L)*WVDX - CVE(L)*WVDY         ! TO LOCAL CURVI-LINEAR SYSTEM
      WVCY = -CUN(L)*WVDX + CUE(L)*WVDY
      WV(L).DIR= ATAN2(WVCY,WVCX)               ! CELL-EAST,WAVE (COUNTER-CLOCKWISE [-pi,pi])
      WV(L).FREQ = 2.*PI/WPRD
      WV(L).DISSIPA(KC) = 0.25*g*WV(L).HEISIG**2/WPRD    !Energy dissipation due to breaking wave 
    else
      WV(L).DIR  = 0.
      WV(L).LENGTH     = 0.
      WV(L).FREQ  = 1.
      WV(L).HEISIG   = 0.
      WV(L).DISSIPA(KC) = 0
    endif    
  enddo
  
  WV(2:LA).HEIGHT = WV(2:LA).HEISIG  ! *** ASSIGN INCIDENT SIGNIFICANT WAVE HEIGHT TO APPLIED WAVE HEIGHT
  
  if( IWVCOUNT == NWVTIM )then
    close(UTBL)
    write(*,'(A)') '** WAVE: END OF WAVE RECORD'
  endif 
  END SUBROUTINE

  SUBROUTINE STRPROC(STR,TERM,N)
  character(*),intent(IN) :: STR
  character(*),intent(OUT) :: TERM(:)
  integer(IK4),intent(OUT) :: N
  integer(IK4) :: STRL,I
  character(120) :: SSTR

  SSTR = ADJUSTL(STR)
  STRL = LEN_TRIM(SSTR)
  SSTR = SSTR(2:STRL)
  N = 0
  do while(.TRUE.)
    SSTR = ADJUSTL(SSTR)
    STRL = LEN_TRIM(SSTR)
    if( STRL > 0 )then
      N = N+1
      do I = 1,STRL
        if( ICHAR(SSTR(I:I)) == 32) EXIT
      enddo
      TERM(N) = SSTR(1:I-1)
      SSTR = SSTR(I:STRL)
    else
      EXIT
    endif
  enddo

  END SUBROUTINE

  SUBROUTINE GETWAVEDAY
  
   
  ! *** READ WAVETIME.INP 
  ! *** OUTPUT:
  ! *** WAVEDAY(1:NWVTIM)
  character(200) :: STR
  integer :: NP,IOS

  if( process_id == master_id )then
    write(*,'(A)')'READING WAVETIME.INP'
    open(1,FILE = 'wavetime.inp',ACTION = 'READ')
    STR = READSTR(1)  
    NP = 0
    do while(.TRUE.)
      read(1,'(A)',end = 100,IOSTAT = IOS) STR
      if( IOS > 0 ) STOP 'WAVETIME.INP: READING ERROR'  
      NP = NP + 1
    enddo
    100 REWIND(1) 
    NWVTIM = NP
  endif
  call Broadcast_Scalar(NWVTIM, master_id)
  
  allocate(WAVEDAY(NWVTIM+1))
  
  if( process_id == master_id )then 
    STR = READSTR(1)  
    do NP = 1,NWVTIM
      read(1,*) WAVEDAY(NP)
    enddo
    close(1)
  endif
  
 Call Broadcast_Array (WAVEDAY, master_id)
  
  END SUBROUTINE
    
END MODULE
