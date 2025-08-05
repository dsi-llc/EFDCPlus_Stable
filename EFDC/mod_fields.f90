! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2024 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE FIELDS
  use GLOBAL
  use INFOMOD
  use Variables_MPI
  use Variables_MPI_Write_Out
  use Variables_MPI_Mapping
  
  use Broadcast_Routines

  implicit none

  type TFIELD                                ! Time and space variable field data
    character :: FNAME*12                    ! Filename for debugging
    integer :: IFLAG                         ! File format (0: Not used, 1: ASCII format, 2: Binary format)
    integer :: ITYPE                         ! Input type (0: value array, 1: LC,IC,JC,values, ...)
    integer :: NT                            ! Number of time steps
    integer :: NC                            ! Number of classes/components
    integer :: NK                            ! Number of layers
    integer :: ITOPT                         ! Interpolation option (0 = No interpolation, 1 = Linear interpolation)
    integer :: IUPDATE                       ! Update flag (0 = Replace, 1 = Add, 2 = Min., 3 = Max.)
    integer :: IDIST                         ! Distribution option (0 = Not multiplied with area, 1 = Multiplied with area)
    integer :: YY,MM,DD                      ! Base date
    integer :: IUNIT                         ! File unit (used with binadry format)
    integer :: NI                            ! Time index for interpolation
    integer, allocatable :: CFLAG(:)         ! Flag for the cells using field data
    real(RK4) :: NODATA                      ! No data value (binary input is fixed for REAL*4)
    real(RK4) :: TSCALE                      ! Time conversion factor to seconds (should > 0, default 86400)
    real(RK4) :: TSHIFT                      ! Time shift (the same units as times, default 0)
    real(RK4) :: VSCALE                      ! Value conversion factor to model units (should > 0, default 1)
    real(RK4) :: VSHIFT                      ! Value shift (the same units as values, default 0)
    real(RKD),allocatable :: TIMES(:)        ! Times (binary input is fixed for REAL*8)
    real(RK4),allocatable :: VALUES(:,:,:,:) ! Values (times,components,cells,layers) (binary input is fixed for REAL*4)
  end type 

  type(TFIELD) :: BATHY                      ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
  type(TFIELD) :: ROUGH                      ! Roughness (e.g., seasonal roughness)
  type(TFIELD) :: VEGE                       ! Vegetation field
  type(TFIELD) :: GWSP                       ! Seepage/groundwater flow
  type(TFIELD) :: WIND                       ! Wind field (e.g., cyclone)
  type(TFIELD) :: PRESS                      ! Barometric pressure field (e.g., cyclone)
  type(TFIELD) :: RAIN                       ! Rainfall
  type(TFIELD) :: EVAP                       ! Evaporation
  type(TFIELD) :: SHELT                      ! Wind shelter field
  type(TFIELD) :: SHADE                      ! Sun shading field
  type(TFIELD) :: SNOW                       ! Snow (snow depth, snowfall)
  type(TFIELD) :: ICETHK                     ! Ice thickness
  type(TFIELD) :: SEDZLJER                   ! SedZLJ erosion rate (NC = 2 includes both multiplier and exponent)

  real(RK4), Allocatable, Dimension(:,:,:,:) :: VALUES_TEMP    !< Temporary variable to help with reading in and remapping to subdomains

contains

  SUBROUTINE INITFIELDS()
    BATHY.IFLAG = 0                        ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
    ROUGH.IFLAG = 0                        ! Roughness (e.g., seasonal roughness)
    VEGE.IFLAG = 0                         ! Vegetation field
    GWSP.IFLAG = 0                         ! Seepage/groundwater flow
    WIND.IFLAG = 0                         ! Wind field (e.g., cyclone)
    PRESS.IFLAG = 0                        ! Barometric pressure field (e.g., cyclone)
    RAIN.IFLAG = 0                         ! Rainfall
    EVAP.IFLAG = 0                         ! Evaporation
    SHELT.IFLAG = 0                        ! Wind shelter field
    SHADE.IFLAG = 0                        ! Sun shading field
    SNOW.IFLAG = 0                         ! Snow depth
    ICETHK.IFLAG = 0                       ! Ice thickness
    SEDZLJER.IFLAG = 0                     ! SedZLJ erosion rate multiplier & exponent
  END SUBROUTINE

  SUBROUTINE READFIELDS()

    call READFIELD(BATHY, "bathfld")         ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
    call READFIELD(ROUGH, "roughfld")        ! Roughness (e.g., seasonal roughness)
    call READFIELD(VEGE, "vegefld")          ! Vegetation field (this can have multiple components with NC>1)
    call READFIELD(GWSP, "gwspfld")          ! Groundwater/seepage flow
    call READFIELD(WIND, "windfld")          ! Wind field (e.g., cyclone)
    call READFIELD(PRESS, "presfld")         ! Barometric pressure field (e.g., cyclone)
    call READFIELD(RAIN, "rainfld")          ! Rainfall
    call READFIELD(EVAP, "evapfld")          ! Evaporation
    call READFIELD(SHELT, "shelfld")         ! Wind shelter field
    call READFIELD(SHADE, "shadfld")         ! Sun shading field
    call READFIELD(SNOW, "snowfld")          ! Snow depth (these two can be combined with NC = 2)
    call READFIELD(ICETHK, "icefld")         ! Ice thickness
    call READFIELD(SEDZLJER, "sedzljfld")    ! SedZLJ erosion rate multiplier (these two can be combined with NC = 2)
  END SUBROUTINE

  SUBROUTINE FREEFIELDS()
    call FREEFIELD(BATHY)                    ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
    call FREEFIELD(ROUGH)                    ! Roughness (e.g., seasonal roughness)
    call FREEFIELD(VEGE)                     ! Vegetation field
    call FREEFIELD(GWSP)                     ! Seepage/groundwater flow
    call FREEFIELD(WIND)                     ! Wind field (e.g., cyclone)
    call FREEFIELD(PRESS)                    ! Barometric pressure field (e.g., cyclone)
    call FREEFIELD(RAIN)                     ! Rainfall
    call FREEFIELD(EVAP)                     ! Evaporation
    call FREEFIELD(SHELT)                    ! Wind shelter field
    call FREEFIELD(SHADE)                    ! Sun shading field
    call FREEFIELD(SNOW)                     ! Snow depth
    call FREEFIELD(ICETHK)                   ! Ice thickness
    call FREEFIELD(SEDZLJER)                 ! SedZLJ erosion rate multiplier & exponent
  END SUBROUTINE

  SUBROUTINE READFIELD(DAT,FILENAME)

  type(TFIELD), intent(INOUT) :: DAT
  character(*), intent(IN) :: FILENAME
  character(127) :: LINE
  integer :: I,K,L,LCC,M,ISO,II,JJ,LL,NN,J
  INTEGER*4 :: SIGNATURE,L2
  real(RK4) :: VALUE
  real,allocatable :: VALUES(:,:)
  logical :: ISEXIST
  integer :: ierr
  
  if( DAT.IFLAG == 0 ) return
    
  allocate(DAT.TIMES(2))
  allocate(DAT.CFLAG(LA_Global))
  DAT.CFLAG = 0

  ! *** Only read/write to screen on the master
  if( process_id == master_id )then

    write(*,*) 'READING ',FILENAME,'...'
    if( DAT.IFLAG == 1 )then
      ! *** ASCII
      LINE = TRIM(FILENAME)//'.inp'
    
      INQUIRE(FILE = LINE,EXIST = ISEXIST)
      if( ISEXIST )then
        ! *** ASCII FIELD FILE
        open(NEWUNIT = DAT.IUNIT,FILE = LINE,STATUS = 'UNKNOWN')
      
        ! *** SKIP OVER HEADER COMMENTS
#ifdef GNU    
        do while(.true.)
          read(DAT.IUNIT,'(A)',end=100) LINE
          if(LINE(1:1)/='*' .AND. LINE(1:1)/='$' .AND. LINE(1:1)/='!' .AND. LINE(1:1)/='C') EXIT
        enddo
100     continue
#else        
        do while(.not. EOF(DAT.IUNIT))
          read(DAT.IUNIT,'(A)') LINE
          if(LINE(1:1)/='*' .and. LINE(1:1)/='$' .and. LINE(1:1)/='!' .and. LINE(1:1)/='C') EXIT
        enddo
#endif 
        read(LINE,*,IOSTAT = ISO) DAT.ITYPE,  DAT.NT,     DAT.NC, LCC, DAT.NK,     DAT.ITOPT, DAT.IUPDATE, DAT.IDIST, &
                                DAT.NODATA, DAT.TSCALE, DAT.TSHIFT,  DAT.VSCALE, DAT.VSHIFT !, DAT.YY, DAT.MM, DAT.DD

        ! *** QC
        if( LCC /= LA_Global - 1 )then
          write(*,*) 'WRONG CELL COUNT FOR ',FILENAME,'...'
          close(DAT.IUNIT)
          return
        endif
        if( DAT.TSCALE <= 0 )then
          write(*,*) 'INVALID FIELD TSCALE. BOTH TSCALE AND TSHIFT ARE RESET.'
          DAT.TSCALE = 86400.
          DAT.TSHIFT = 0.
        endif
        if( ABS(DAT.VSCALE) <= 1.0e-15 )then
          write(*,*) 'INVALID FIELD VSCALE. BOTH VSCALE AND VSHIFT ARE RESET.'
          DAT.VSCALE = 1.
          DAT.VSHIFT = 0.
        endif

        allocate(VALUES_TEMP(2,DAT.NC,LC_Global,DAT.NK))
        allocate(VALUES(DAT.NC,DAT.NK))
        VALUES_TEMP = DAT.NODATA

        ! *** READ THE FIRST TWO FIELD SNAPSHOTS
        do NN = 1,2
          read(DAT.IUNIT,*,IOSTAT = ISO) DAT.TIMES(NN),LCC
          DAT.TIMES(NN) = (DAT.TIMES(NN) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days
          
          if( NN > 1 )then
            if( DAT.TIMES(NN) < DAT.TIMES(NN-1) )then
              write(*,*) 'TIME AT ',NN-1,' = ',DAT.TIMES(NN-1),', AT ',NN,' = ',DAT.TIMES(NN)
            endif
          endif
          
          if( DAT.ITYPE > 0 )then
            ! *** L,I,J,VALUE FORMAT
            do LL = 1,LCC
              read(DAT.IUNIT,*,IOSTAT = ISO) L2,II,JJ,((VALUES(M,K), K = 1,DAT.NK), M = 1,DAT.NC)
              L = LIJ_Global(II,JJ)
              DAT.CFLAG(L) = 1
              do M = 1,DAT.NC
                do K = 1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES(M,K) + DAT.VSHIFT)*DAT.VSCALE
                enddo
              enddo
            enddo
          else
            ! *** READ ALL CELLS IN L ARRAY ORDER
            read(DAT.IUNIT,*,IOSTAT = ISO) (((VALUES_TEMP(NN,M,L,K),K = 1,DAT.NK),L = 2,LA_Global),M = 1,DAT.NC)
            do M = 1,DAT.NC
              do L = 2,LA_Global
                do K = 1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES_TEMP(NN,M,L,K) + DAT.VSHIFT)*DAT.VSCALE
                enddo
              enddo
            enddo
          endif
          
        enddo
        deallocate(VALUES)
      else
        write(*,*) 'FILE NOT FOUND: ',FILENAME,'.inp'
        DAT.IFLAG = 0
        return
      endif
    else
      ! *** BINARY
      LINE = TRIM(FILENAME)//'.fld'
      INQUIRE(FILE = LINE,EXIST = ISEXIST)
      if(ISEXIST )then
        open(NEWUNIT = DAT.IUNIT,FILE = LINE,STATUS = 'UNKNOWN',ACCESS = 'SEQUENTIAL',FORM = FMT_BINARY,ACTION = 'READ',SHARE  = 'DENYNONE')
        read(DAT.IUNIT) SIGNATURE
        read(DAT.IUNIT) DAT.ITYPE
        read(DAT.IUNIT) DAT.NT
        read(DAT.IUNIT) DAT.NC
        read(DAT.IUNIT) LCC
        read(DAT.IUNIT) DAT.NK
        read(DAT.IUNIT) DAT.ITOPT
        read(DAT.IUNIT) DAT.IUPDATE
        read(DAT.IUNIT) DAT.IDIST
        read(DAT.IUNIT) DAT.NODATA
        read(DAT.IUNIT) DAT.TSCALE
        read(DAT.IUNIT) DAT.TSHIFT
        read(DAT.IUNIT) DAT.VSCALE
        read(DAT.IUNIT) DAT.VSHIFT
        read(DAT.IUNIT) DAT.YY
        read(DAT.IUNIT) DAT.MM
        read(DAT.IUNIT) DAT.DD
        read(DAT.IUNIT) L2
        read(DAT.IUNIT) L2
        read(DAT.IUNIT) L2

        ! *** QC
        if( LCC /= LA_Global - 1 )then
          write(*,*) 'WRONG CELL COUNT FOR ',FILENAME,'...'
          close(DAT.IUNIT)
          return
        endif
        if( DAT.TSCALE <= 0 )then
          write(*,*) 'INVALID FIELD TSCALE. BOTH TSCALE AND TSHIFT ARE RESET.'
          DAT.TSCALE = 86400.
          DAT.TSHIFT = 0.
        endif
        if( ABS(DAT.VSCALE) <= 1.0e-15 )then
          write(*,*) 'INVALID FIELD VSCALE. BOTH VSCALE AND VSHIFT ARE RESET.'
          DAT.VSCALE = 1.
          
          DAT.VSHIFT = 0.
        endif

        allocate(VALUES_TEMP(2,DAT.NC,LC_global,DAT.NK))
        VALUES_TEMP(:,:,:,:) = DAT.NODATA

        ! *** READ THE FIRST TWO FIELD SNAPSHOTS
        do NN = 1,2
          read(DAT.IUNIT) DAT.TIMES(NN)
          read(DAT.IUNIT) LCC
          DAT.TIMES(NN) = (DAT.TIMES(NN) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days

          if( NN > 1 )then
            if( DAT.TIMES(NN) < DAT.TIMES(NN-1) )then
              write(*,*) 'TIME AT ',NN-1,' = ',DAT.TIMES(NN-1),', AT ',NN,' = ',DAT.TIMES(NN)
            endif
          endif
                    
          if( DAT.ITYPE > 0 )then
            ! *** READ ONLY SPECIFIC CELLS (FORMAT: L, I, J, VALUES)
            do LL = 1,LCC
              read(DAT.IUNIT) L2
              read(DAT.IUNIT) II
              read(DAT.IUNIT) JJ
              L = LIJ_Global(II,JJ)
              read(DAT.IUNIT) ((VALUES_TEMP(NN,M,L,K), K = 1,DAT.NK), M = 1,DAT.NC)
              do M = 1,DAT.NC
                do K = 1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES_TEMP(NN,M,L,K) + DAT.VSHIFT)*DAT.VSCALE
                enddo
              enddo
            enddo
          else
            ! *** READ ALL CELLS
            read(DAT.IUNIT) (((VALUES_TEMP(NN,M,L,K), K = 1,DAT.NK), L = 2,LA_Global), M = 1,DAT.NC)
            do M = 1,DAT.NC
              do L = 2,LA_Global
                do K = 1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES_TEMP(NN,M,L,K) + DAT.VSHIFT)*DAT.VSCALE
                enddo
              enddo
            enddo
          endif
        enddo
      else
        write(*,*) 'FILE NOT FOUND: ',FILENAME,'.fld'
        DAT.IFLAG = 0
        call STOPP('.')
      endif
    endif
  endif
  
  ! *** broadcast header info just read
  call Broadcast_Scalar(DAT.ITYPE,   master_id)
  call Broadcast_Scalar(DAT.NT,      master_id)
  call Broadcast_Scalar(DAT.NC,      master_id)
  call Broadcast_Scalar(DAT.NK,      master_id)
  call Broadcast_Scalar(DAT.ITOPT,   master_id)
  call Broadcast_Scalar(DAT.IUPDATE, master_id)
  call Broadcast_Scalar(DAT.IDIST,   master_id)
  call Broadcast_Scalar(DAT.NODATA,  master_id)
  call Broadcast_Scalar(DAT.TSCALE,  master_id)
  call Broadcast_Scalar(DAT.TSHIFT,  master_id)
  call Broadcast_Scalar(DAT.VSCALE,  master_id)
  call Broadcast_Scalar(DAT.VSHIFT,  master_id)
  call Broadcast_Array(DAT.TIMES,    master_id)
  call Broadcast_Array(DAT.CFLAG,    master_id)
  
  if( process_id /= master_id )then
    allocate(VALUES_TEMP(2,DAT.NC,LC_Global,DAT.NK))
    VALUES_TEMP = DAT.NODATA
  endif
  call Broadcast_Array(VALUES_TEMP, master_id)
  
  ! *** NOW ALL PROCESSES HAVE ALL THE INPUTS.  ASSIGN LOCAL VALUES
  allocate(DAT.VALUES(2,DAT.NC,LC,DAT.NK))
  DAT.VALUES = DAT.NODATA
  
  ! *** Map from global --> local values
  do NN = 1,2
    do M = 1,DAT.NC
      do L = 2,LA
        do K = 1,DAT.NK
          L2 = Map2Global(L).LG     ! *** Global L that corresponds to Local L
          DAT.VALUES(NN,M,L,K) = VALUES_TEMP(NN,M,L2,K) 
        enddo
      enddo
    enddo
  enddo  
  deallocate(VALUES_TEMP)
  
  DAT.NI = 1

END SUBROUTINE

  SUBROUTINE FREEFIELD(DAT)

    type(TFIELD), intent(INOUT) :: DAT
    if(DAT.IFLAG == 0) return
    deallocate(DAT.TIMES,DAT.VALUES)
    if(DAT.NI > 0) close(DAT.IUNIT)

  END SUBROUTINE

  SUBROUTINE UPDATEFIELD(DAT,TIME,M,VALUE)

    ! *** SETS THE ARRAY "VALUE" FOR THE CURRENT TIME STEP.
    ! *** READS FROM FILE IF NEEDED AND APPLIES AREAL CONVERSIONS
    ! *** DELME - THIS USES A LOT OF LOOPS THAT ARE IN THE OMP "SINGLE" REGION.  MAY NEED TO MOVE TO BETTER UTILIZE OMP

    ! *** DAT   - FIELD STRUCTURE
    ! *** TIME  - CURRENT TIME IN SECONDS
    ! *** M     - COMPONENT
    ! *** VALUE - FINAL VALUE OF ARRAY AFTER ALL FIELD OPERATIONS

  type(TFIELD), intent(INOUT) :: DAT
  real(RKD), intent(IN) :: TIME
  integer, intent(IN) :: M
  real, intent(INOUT) :: VALUE(LCM)
  real :: RES(LC_Global)
  real*4,allocatable :: VALUES(:,:)
  real(RKD) :: DELTF, FAC
  integer :: I,K,L,NCC,ISO,LCC,II,JJ,LL,L2

  if( DAT.IFLAG == 0 ) return
  if( TIME < DAT.TIMES(1) ) return
   
  RES(:) = DAT.NODATA
   
  ! *** ENSURE THE CURRENT TIME IS BETWEEN THE TWO FIELD SNAPSHOTS
  do while (TIME > DAT.TIMES(2) .and. DAT.NI < DAT.NT)
    ! *** CURRENT TIME IS OUT OF RANGE

    ! *** ADVANCE THE TIME
    DAT.TIMES(1) = DAT.TIMES(2)
    DAT.VALUES(1,:,:,:) = DAT.VALUES(2,:,:,:)

    allocate(VALUES_TEMP(2,DAT.NC,LC_Global,DAT.NK))
    VALUES_TEMP(:,:,:,:) = DAT.NODATA

    if( process_id == master_id )then
      ! *** READ IN THE NEXT TIME STEP
      if( DAT.IFLAG == 1 )then
        ! *** ASCII
        read(DAT.IUNIT,*,IOSTAT = ISO) DAT.TIMES(2),LCC
        DAT.TIMES(2) = (DAT.TIMES(2) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days

        allocate(VALUES(DAT.NC,DAT.NK))

        if( DAT.ITYPE > 0 )then
          ! *** READ ONLY SPECIFIC CELLS (FORMAT: L, I, J, VALUES)
          !READ(DAT.IUNIT,*,IOSTAT = ISO) L,II,JJ,(((DAT.VALUES(2,NCC,L,K), K = 1,DAT.NK), LL = 1,LCC), NCC = 1,DAT.NC)
          do LL = 1,LCC
            read(DAT.IUNIT,*,IOSTAT = ISO) L2,II,JJ,((VALUES(NCC,K), K = 1,DAT.NK), NCC = 1,DAT.NC)
            L = LIJ_Global(II,JJ)
            do NCC = 1,DAT.NC
              do K = 1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES(NCC,K) + DAT.VSHIFT)*DAT.VSCALE
              enddo
            enddo
          enddo
        else
          ! *** READ ALL CELLS
          read(DAT.IUNIT,*,IOSTAT = ISO) (((VALUES_TEMP(2,NCC,L,K), K = 1,DAT.NK), L = 2,LA_Global), NCC = 1,DAT.NC)
          do NCC = 1,DAT.NC
            do L = 2,LA_Global
              do K = 1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES_TEMP(2,NCC,L,K) + DAT.VSHIFT)*DAT.VSCALE
              enddo
            enddo
          enddo
        endif
        deallocate(VALUES)
      else
        ! *** BINARY
        read(DAT.IUNIT) DAT.TIMES(2)
        read(DAT.IUNIT) LCC
        DAT.TIMES(2) = (DAT.TIMES(2) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days

        if( DAT.ITYPE > 0 )then
          ! *** READ ONLY SPECIFIC CELLS (FORMAT: L, I, J, VALUES)
          !READ(DAT.IUNIT) L,II,JJ,(((DAT.VALUES(2,NCC,L,K), K = 1,DAT.NK), LL = 1,LCC), NCC = 1,DAT.NC)
          do LL = 1,LCC
            read(DAT.IUNIT) L2
            read(DAT.IUNIT) II
            read(DAT.IUNIT) JJ
            L = LIJ_Global(II,JJ)
            read(DAT.IUNIT) ((VALUES_TEMP(2,NCC,L,K), K = 1,DAT.NK), NCC = 1,DAT.NC)
            do NCC = 1,DAT.NC
              do K = 1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES_TEMP(2,NCC,L,K) + DAT.VSHIFT)*DAT.VSCALE
              enddo
            enddo
          enddo
        else
          ! *** READ ALL CELLS
          read(DAT.IUNIT) (((VALUES_TEMP(2,NCC,L,K), K = 1,DAT.NK), L = 2,LA_Global), NCC = 1,DAT.NC)
          do NCC = 1,DAT.NC
            do L = 2,LA_Global
              do K = 1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES_TEMP(2,NCC,L,K) + DAT.VSHIFT)*DAT.VSCALE
              enddo
            enddo
          enddo
        endif
      endif
    endif !***End on master
    
    ! *** broadcast header info just read
    call Broadcast_Scalar(DAT.TIMES(2), master_id)
    call Broadcast_Array(VALUES_TEMP,  master_id)

    ! *** NOW ALL PROCESSES HAVE ALL THE INPUTS.  ASSIGN LOCAL VALUES
    do NCC = 1,DAT.NC
      do L = 2,LA
        do K = 1,DAT.NK
          L2 = Map2Global(L).LG                              ! *** Global L that corresponds to Local L
          DAT.VALUES(2,NCC,L,K) = VALUES_TEMP(2,NCC,L2,K)
        enddo
      enddo
    enddo
    deallocate(VALUES_TEMP)
  
    DAT.NI = DAT.NI + 1
  enddo

    K = 1   ! *** TODO: LAYERED DATA

    if( DAT.ITOPT == 1  )then
      ! *** INTERPOLATE FIELD
      I = 1
      DELTF = DAT.TIMES(I+1) - DAT.TIMES(I)
      if( DELTF > 1.0E-9 )then
        FAC = (TIME - DAT.TIMES(I))/DELTF
      else
        FAC = 0.5
      endif
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(L)
      do L = 2,LA
        if(DAT.VALUES(I,M,L,K) == DAT.NODATA .or. DAT.VALUES(I+1,M,L,K) == DAT.NODATA )then
          RES(L) = DAT.NODATA
        else
          RES(L) = DAT.VALUES(I,M,L,K) + FAC*(DAT.VALUES(I+1,M,L,K)-DAT.VALUES(I,M,L,K))
        endif
      enddo
      !$OMP END PARALLEL DO
    else
      ! *** NO INTERPOLATION
      RES(:) = DAT.VALUES(1,M,:,K)
    endif

  ! *** ADJUST FOR AREA, IF NEEDED (DEFAULT IS FLUX: IDIST == 0)
  if( DAT.IDIST == 1 )then
    ! *** INPUT DATA IS A RATE
    do L = 2,LA
      if( RES(L) /= DAT.NODATA )then
        RES(L) = RES(L)*DXYP(L)
      endif
    enddo
  elseif( DAT.IDIST == 2 )then
    ! *** INPUT DATA IS A TOTAL VOLUME
    do L = 2,LA
      if( RES(L) /= DAT.NODATA )then
        RES(L) = RES(L)/DXYP(L)
      endif
    enddo
  endif

    ! UPDATE FIELD
    SELECT CASE (DAT.IUPDATE)
    CASE (1)      ! *** ADD
      do L = 1,LC
        if( RES(L) /= DAT.NODATA )then
          VALUE(L) = VALUE(L) + RES(L)
        endif
      enddo
    CASE (2)      ! *** MIN
      do L = 1,LC
        if( RES(L) /= DAT.NODATA .and. RES(L) < VALUE(L) )then
          VALUE(L) = RES(L)
        endif
      enddo
    CASE (3)      ! *** MAX
      do L = 1,LC
        if( RES(L) /= DAT.NODATA .and. RES(L) > VALUE(L) )then
          VALUE(L) = RES(L)
        endif
      enddo
    CASE DEFAULT  ! *** REPLACE
      do L = 1,LC
        if( RES(L) /= DAT.NODATA )then
          VALUE(L) = RES(L)
        endif
      enddo
    END SELECT
END SUBROUTINE

  SUBROUTINE UPDATETOPO(TIME,BELV,HP,HDRY)
    real,dimension(:),intent(INOUT) :: BELV,HP
    real(RKD), intent(IN) :: TIME
    real, intent(IN) :: HDRY
    real,allocatable,dimension(:) :: WSEL
    integer :: L

    ! *** Delme - This approach does not maintain mass balance.  use with care!
    allocate(WSEL(LA))
    do L = 2,LA
      WSEL(L) = HP(L) + BELV(L)
    enddo
    call UPDATEFIELD(BATHY,TIME,1,BELV)
    do L = 2,LA
      if(LMASKDRY(L) .and. BATHY.CFLAG(L) > 0 )then
        HP(L) = WSEL(L) - BELV(L)
        if( HP(L) < HDRY) HP(L) = HDRY
      endif
    enddo
    !L = 1438
    !OPEN(910,FILE = OUTDIR//'FIELD.LOG',POSITION = 'APPEND')
    !write(910,*) 'T = ',TIME,', B = ',BELV(L),', H = ',HP(L),', H = ',WSEL(L)
    !CLOSE(910)
    deallocate(WSEL)
END SUBROUTINE

END MODULE
