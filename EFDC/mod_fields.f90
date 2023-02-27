! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE FIELDS
  USE GLOBAL
  USE INFOMOD
  Use Variables_MPI
  Use Variables_MPI_Write_Out
  Use Variables_MPI_Mapping
  
#ifdef _MPI
  Use Broadcast_Routines
#endif

  IMPLICIT NONE

  TYPE TFIELD                                ! Time and space variable field data
    CHARACTER :: FNAME*12                    ! Filename for debugging
    INTEGER :: IFLAG                         ! File format (0: Not used, 1: ASCII format, 2: Binary format)
    INTEGER :: ITYPE                         ! Input type (0: value array, 1: LC,IC,JC,values, ...)
    INTEGER :: NT                            ! Number of time steps
    INTEGER :: NC                            ! Number of classes/components
    INTEGER :: NK                            ! Number of layers
    INTEGER :: ITOPT                         ! Interpolation option (0 = No interpolation, 1 = Linear interpolation)
    INTEGER :: IUPDATE                       ! Update flag (0 = Replace, 1 = Add, 2 = Min., 3 = Max.)
    INTEGER :: IDIST                         ! Distribution option (0 = Not multiplied with area, 1 = Multiplied with area)
    INTEGER :: YY,MM,DD                      ! Base date
    INTEGER :: IUNIT                         ! File unit (used with binadry format)
    INTEGER :: NI                            ! Time index for interpolation
    INTEGER, ALLOCATABLE :: CFLAG(:)         ! Flag for the cells using field data
    REAL(RK4) :: NODATA                      ! No data value (binary input is fixed for REAL*4)
    REAL(RK4) :: TSCALE                      ! Time conversion factor to seconds (should > 0, default 86400)
    REAL(RK4) :: TSHIFT                      ! Time shift (the same units as times, default 0)
    REAL(RK4) :: VSCALE                      ! Value conversion factor to model units (should > 0, default 1)
    REAL(RK4) :: VSHIFT                      ! Value shift (the same units as values, default 0)
    REAL(RKD),ALLOCATABLE :: TIMES(:)        ! Times (binary input is fixed for REAL*8)
    REAL(RK4),ALLOCATABLE :: VALUES(:,:,:,:) ! Values (times,components,cells,layers) (binary input is fixed for REAL*4)
  END TYPE

  TYPE(TFIELD) :: BATHY                      ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
  TYPE(TFIELD) :: ROUGH                      ! Roughness (e.g., seasonal roughness)
  TYPE(TFIELD) :: VEGE                       ! Vegetation field
  TYPE(TFIELD) :: GWSP                       ! Seepage/groundwater flow
  TYPE(TFIELD) :: WIND                       ! Wind field (e.g., cyclone)
  TYPE(TFIELD) :: PRESS                      ! Barometric pressure field (e.g., cyclone)
  TYPE(TFIELD) :: RAIN                       ! Rainfall
  TYPE(TFIELD) :: EVAP                       ! Evaporation
  TYPE(TFIELD) :: SHELT                      ! Wind shelter field
  TYPE(TFIELD) :: SHADE                      ! Sun shading field
  TYPE(TFIELD) :: SNOW                       ! Snow (snow depth, snowfall)
  TYPE(TFIELD) :: ICETHK                     ! Ice thickness
  TYPE(TFIELD) :: SEDZLJER                   ! SedZLJ erosion rate (NC=2 includes both multiplier and exponent)

  Real(RK4), Allocatable, Dimension(:,:,:,:) :: VALUES_TEMP    !< Temporary variable to help with reading in and remapping to subdomains

CONTAINS

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

    CALL READFIELD(BATHY, "bathfld")         ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
    CALL READFIELD(ROUGH, "roughfld")        ! Roughness (e.g., seasonal roughness)
    CALL READFIELD(VEGE, "vegefld")          ! Vegetation field (this can have multiple components with NC>1)
    CALL READFIELD(GWSP, "gwspfld")          ! Groundwater/seepage flow
    CALL READFIELD(WIND, "windfld")          ! Wind field (e.g., cyclone)
    CALL READFIELD(PRESS, "presfld")         ! Barometric pressure field (e.g., cyclone)
    CALL READFIELD(RAIN, "rainfld")          ! Rainfall
    CALL READFIELD(EVAP, "evapfld")          ! Evaporation
    CALL READFIELD(SHELT, "shelfld")         ! Wind shelter field
    CALL READFIELD(SHADE, "shadfld")         ! Sun shading field
    CALL READFIELD(SNOW, "snowfld")          ! Snow depth (these two can be combined with NC=2)
    CALL READFIELD(ICETHK, "icefld")         ! Ice thickness
    CALL READFIELD(SEDZLJER, "sedzljfld")    ! SedZLJ erosion rate multiplier (these two can be combined with NC=2)
  END SUBROUTINE

  SUBROUTINE FREEFIELDS()
    CALL FREEFIELD(BATHY)                    ! Bathymetric data (e.g., dredging/dumping, land erosion/reclamation)
    CALL FREEFIELD(ROUGH)                    ! Roughness (e.g., seasonal roughness)
    CALL FREEFIELD(VEGE)                     ! Vegetation field
    CALL FREEFIELD(GWSP)                     ! Seepage/groundwater flow
    CALL FREEFIELD(WIND)                     ! Wind field (e.g., cyclone)
    CALL FREEFIELD(PRESS)                    ! Barometric pressure field (e.g., cyclone)
    CALL FREEFIELD(RAIN)                     ! Rainfall
    CALL FREEFIELD(EVAP)                     ! Evaporation
    CALL FREEFIELD(SHELT)                    ! Wind shelter field
    CALL FREEFIELD(SHADE)                    ! Sun shading field
    CALL FREEFIELD(SNOW)                     ! Snow depth
    CALL FREEFIELD(ICETHK)                   ! Ice thickness
    CALL FREEFIELD(SEDZLJER)                 ! SedZLJ erosion rate multiplier & exponent
  END SUBROUTINE

  SUBROUTINE READFIELD(DAT,FILENAME)

  TYPE(TFIELD), INTENT(INOUT) :: DAT
  CHARACTER(*), INTENT(IN) :: FILENAME
  CHARACTER(127) :: LINE
  INTEGER :: I,K,L,LCC,M,ISO,II,JJ,LL,NN,J
  INTEGER*4 :: SIGNATURE,L2
  REAL(RK4) :: VALUE
  REAL,ALLOCATABLE :: VALUES(:,:)
  LOGICAL :: ISEXIST
  Integer :: ierr
  
  IF( DAT.IFLAG == 0 ) RETURN
    
  ALLOCATE(DAT.TIMES(2))
  ALLOCATE(DAT.CFLAG(LA_Global))
  DAT.CFLAG = 0

  ! *** Only read/write to screen on the master
  if( process_id == master_id )THEN

    WRITE(*,*) 'READING ',FILENAME,'...'
    IF( DAT.IFLAG == 1 )THEN
      ! *** ASCII
      LINE = TRIM(FILENAME)//'.inp'
    
      INQUIRE(FILE=LINE,EXIST=ISEXIST)
      IF( ISEXIST )THEN
        ! *** ASCII FIELD FILE
        OPEN(NEWUNIT=DAT.IUNIT,FILE=LINE,STATUS='UNKNOWN')
      
        ! *** SKIP OVER HEADER COMMENTS
        DO WHILE(.NOT. EOF(DAT.IUNIT))
          READ(DAT.IUNIT,'(A)') LINE
          IF(LINE(1:1)/='*' .AND. LINE(1:1)/='$' .AND. LINE(1:1)/='!' .AND. LINE(1:1)/='C') EXIT
        ENDDO
        READ(LINE,*,IOSTAT=ISO) DAT.ITYPE,  DAT.NT,     DAT.NC, LCC, DAT.NK,     DAT.ITOPT, DAT.IUPDATE, DAT.IDIST, &
                                DAT.NODATA, DAT.TSCALE, DAT.TSHIFT,  DAT.VSCALE, DAT.VSHIFT !, DAT.YY, DAT.MM, DAT.DD

        ! *** QC
        IF( LCC /= LA_Global - 1 )THEN
          WRITE(*,*) 'WRONG CELL COUNT FOR ',FILENAME,'...'
          CLOSE(DAT.IUNIT)
          RETURN
        ENDIF
        IF( DAT.TSCALE <= 0 )THEN
          WRITE(*,*) 'INVALID FIELD TSCALE. BOTH TSCALE AND TSHIFT ARE RESET.'
          DAT.TSCALE = 86400.
          DAT.TSHIFT = 0.
        ENDIF
        IF( ABS(DAT.VSCALE) <= 1.0e-15 )THEN
          WRITE(*,*) 'INVALID FIELD VSCALE. BOTH VSCALE AND VSHIFT ARE RESET.'
          DAT.VSCALE = 1.
          DAT.VSHIFT = 0.
        ENDIF

        ALLOCATE(VALUES_TEMP(2,DAT.NC,LC_Global,DAT.NK))
        ALLOCATE(VALUES(DAT.NC,DAT.NK))
        VALUES_TEMP = DAT.NODATA

        ! *** READ THE FIRST TWO FIELD SNAPSHOTS
        DO NN=1,2
          READ(DAT.IUNIT,*,IOSTAT=ISO) DAT.TIMES(NN),LCC
          DAT.TIMES(NN) = (DAT.TIMES(NN) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days
          
          IF( NN > 1 )THEN
            IF( DAT.TIMES(NN) < DAT.TIMES(NN-1) )THEN
              WRITE(*,*) 'TIME AT ',NN-1,'=',DAT.TIMES(NN-1),', AT ',NN,'=',DAT.TIMES(NN)
            ENDIF
          ENDIF
          
          IF( DAT.ITYPE > 0 )THEN
            ! *** L,I,J,VALUE FORMAT
            DO LL=1,LCC
              READ(DAT.IUNIT,*,IOSTAT=ISO) L2,II,JJ,((VALUES(M,K), K=1,DAT.NK), M=1,DAT.NC)
              L = LIJ_Global(II,JJ)
              DAT.CFLAG(L) = 1
              DO M=1,DAT.NC
                DO K=1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES(M,K) + DAT.VSHIFT)*DAT.VSCALE
                ENDDO
              ENDDO
            ENDDO
          ELSE
            ! *** READ ALL CELLS IN L ARRAY ORDER
            READ(DAT.IUNIT,*,IOSTAT=ISO) (((VALUES_TEMP(NN,M,L,K),K=1,DAT.NK),L=2,LA_Global),M=1,DAT.NC)
            DO M=1,DAT.NC
              DO L=2,LA_Global
                DO K=1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES_TEMP(NN,M,L,K) + DAT.VSHIFT)*DAT.VSCALE
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          
        ENDDO
        DEALLOCATE(VALUES)
      ELSE
        WRITE(*,*) 'FILE NOT FOUND: ',FILENAME,'.inp'
        DAT.IFLAG = 0
        RETURN
      ENDIF
    ELSE
      ! *** BINARY
      LINE = TRIM(FILENAME)//'.fld'
      INQUIRE(FILE=LINE,EXIST=ISEXIST)
      IF(ISEXIST )THEN
        OPEN(NEWUNIT=DAT.IUNIT,FILE=LINE,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY',ACTION='READ',SHARE ='DENYNONE')
        READ(DAT.IUNIT) SIGNATURE
        READ(DAT.IUNIT) DAT.ITYPE
        READ(DAT.IUNIT) DAT.NT
        READ(DAT.IUNIT) DAT.NC
        READ(DAT.IUNIT) LCC
        READ(DAT.IUNIT) DAT.NK
        READ(DAT.IUNIT) DAT.ITOPT
        READ(DAT.IUNIT) DAT.IUPDATE
        READ(DAT.IUNIT) DAT.IDIST
        READ(DAT.IUNIT) DAT.NODATA
        READ(DAT.IUNIT) DAT.TSCALE
        READ(DAT.IUNIT) DAT.TSHIFT
        READ(DAT.IUNIT) DAT.VSCALE
        READ(DAT.IUNIT) DAT.VSHIFT
        READ(DAT.IUNIT) DAT.YY
        READ(DAT.IUNIT) DAT.MM
        READ(DAT.IUNIT) DAT.DD
        READ(DAT.IUNIT) L2
        READ(DAT.IUNIT) L2
        READ(DAT.IUNIT) L2

        ! *** QC
        IF( LCC /= LA_Global - 1 )THEN
          WRITE(*,*) 'WRONG CELL COUNT FOR ',FILENAME,'...'
          CLOSE(DAT.IUNIT)
          RETURN
        ENDIF
        IF( DAT.TSCALE <= 0 )THEN
          WRITE(*,*) 'INVALID FIELD TSCALE. BOTH TSCALE AND TSHIFT ARE RESET.'
          DAT.TSCALE = 86400.
          DAT.TSHIFT = 0.
        ENDIF
        IF( ABS(DAT.VSCALE) <= 1.0e-15 )THEN
          WRITE(*,*) 'INVALID FIELD VSCALE. BOTH VSCALE AND VSHIFT ARE RESET.'
          DAT.VSCALE = 1.
          
          DAT.VSHIFT = 0.
        ENDIF

        ALLOCATE(VALUES_TEMP(2,DAT.NC,LC_global,DAT.NK))
        VALUES_TEMP(:,:,:,:) = DAT.NODATA

        ! *** READ THE FIRST TWO FIELD SNAPSHOTS
        DO NN=1,2
          READ(DAT.IUNIT) DAT.TIMES(NN)
          READ(DAT.IUNIT) LCC
          DAT.TIMES(NN) = (DAT.TIMES(NN) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days

          IF( NN > 1 )THEN
            IF( DAT.TIMES(NN) < DAT.TIMES(NN-1) )THEN
              WRITE(*,*) 'TIME AT ',NN-1,'=',DAT.TIMES(NN-1),', AT ',NN,'=',DAT.TIMES(NN)
            ENDIF
          ENDIF
                    
          IF( DAT.ITYPE > 0 )THEN
            ! *** READ ONLY SPECIFIC CELLS (FORMAT: L, I, J, VALUES)
            DO LL=1,LCC
              READ(DAT.IUNIT) L2
              READ(DAT.IUNIT) II
              READ(DAT.IUNIT) JJ
              L = LIJ_Global(II,JJ)
              READ(DAT.IUNIT) ((VALUES_TEMP(NN,M,L,K), K=1,DAT.NK), M=1,DAT.NC)
              DO M=1,DAT.NC
                DO K=1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES_TEMP(NN,M,L,K) + DAT.VSHIFT)*DAT.VSCALE
                ENDDO
              ENDDO
            ENDDO
          ELSE
            ! *** READ ALL CELLS
            READ(DAT.IUNIT) (((VALUES_TEMP(NN,M,L,K), K=1,DAT.NK), L=2,LA_Global), M=1,DAT.NC)
            DO M=1,DAT.NC
              DO L=2,LA_Global
                DO K=1,DAT.NK
                  VALUES_TEMP(NN,M,L,K) = (VALUES_TEMP(NN,M,L,K) + DAT.VSHIFT)*DAT.VSCALE
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE
        WRITE(*,*) 'FILE NOT FOUND: ',FILENAME,'.fld'
        DAT.IFLAG = 0
        CALL STOPP('.')
      ENDIF
    ENDIF
  ENDIF
  
  ! *** broadcast header info just read
  Call Broadcast_Scalar(DAT.ITYPE,   master_id)
  Call Broadcast_Scalar(DAT.NT,      master_id)
  Call Broadcast_Scalar(DAT.NC,      master_id)
  Call Broadcast_Scalar(DAT.NK,      master_id)
  Call Broadcast_Scalar(DAT.ITOPT,   master_id)
  Call Broadcast_Scalar(DAT.IUPDATE, master_id)
  Call Broadcast_Scalar(DAT.IDIST,   master_id)
  Call Broadcast_Scalar(DAT.NODATA,  master_id)
  Call Broadcast_Scalar(DAT.TSCALE,  master_id)
  Call Broadcast_Scalar(DAT.TSHIFT,  master_id)
  Call Broadcast_Scalar(DAT.VSCALE,  master_id)
  Call Broadcast_Scalar(DAT.VSHIFT,  master_id)
  Call Broadcast_Array(DAT.TIMES,    master_id)
  Call Broadcast_Array(DAT.CFLAG,    master_id)
  
  if( process_id /= master_id )THEN
    ALLOCATE(VALUES_TEMP(2,DAT.NC,LC_Global,DAT.NK))
    VALUES_TEMP = DAT.NODATA
  endif
  Call Broadcast_Array(VALUES_TEMP, master_id)
  
  ! *** NOW ALL PROCESSES HAVE ALL THE INPUTS.  ASSIGN LOCAL VALUES
  ALLOCATE(DAT.VALUES(2,DAT.NC,LC,DAT.NK))
  DAT.VALUES = DAT.NODATA
  
  ! *** Map from global --> local values
  DO NN=1,2
    DO M=1,DAT.NC
      DO L=2,LA
        DO K=1,DAT.NK
          L2 = Map2Global(L).LG     ! *** Global L that corresponds to Local L
          DAT.VALUES(NN,M,L,K) = VALUES_TEMP(NN,M,L2,K) 
        ENDDO
      ENDDO
    ENDDO
  ENDDO  
  DEALLOCATE(VALUES_TEMP)
  
  DAT.NI = 1

END SUBROUTINE

  SUBROUTINE FREEFIELD(DAT)

    TYPE(TFIELD), INTENT(INOUT) :: DAT
    IF(DAT.IFLAG == 0) RETURN
    DEALLOCATE(DAT.TIMES,DAT.VALUES)
    IF(DAT.NI > 0) CLOSE(DAT.IUNIT)

  END SUBROUTINE

  SUBROUTINE UPDATEFIELD(DAT,TIME,M,VALUE)

    ! *** SETS THE ARRAY "VALUE" FOR THE CURRENT TIME STEP.
    ! *** READS FROM FILE IF NEEDED AND APPLIES AREAL CONVERSIONS
    ! *** DELME - THIS USES A LOT OF LOOPS THAT ARE IN THE OMP "SINGLE" REGION.  MAY NEED TO MOVE TO BETTER UTILIZE OMP

    ! *** DAT   - FIELD STRUCTURE
    ! *** TIME  - CURRENT TIME IN SECONDS
    ! *** M     - COMPONENT
    ! *** VALUE - FINAL VALUE OF ARRAY AFTER ALL FIELD OPERATIONS

  TYPE(TFIELD), INTENT(INOUT) :: DAT
  REAL(RKD), INTENT(IN) :: TIME
  INTEGER, INTENT(IN) :: M
  REAL, INTENT(INOUT) :: VALUE(LCM)
  REAL :: RES(LC_Global)
  REAL*4,ALLOCATABLE :: VALUES(:,:)
  REAL(RKD) :: DELTF, FAC
  INTEGER :: I,K,L,NCC,ISO,LCC,II,JJ,LL,L2

  IF( DAT.IFLAG == 0 ) RETURN
  IF( TIME < DAT.TIMES(1) ) RETURN
   
  RES(:) = DAT.NODATA
   
  ! *** ENSURE THE CURRENT TIME IS BETWEEN THE TWO FIELD SNAPSHOTS
  DO WHILE (TIME > DAT.TIMES(2) .AND. DAT.NI < DAT.NT)
    ! *** CURRENT TIME IS OUT OF RANGE

    ! *** ADVANCE THE TIME
    DAT.TIMES(1) = DAT.TIMES(2)
    DAT.VALUES(1,:,:,:) = DAT.VALUES(2,:,:,:)

    ALLOCATE(VALUES_TEMP(2,DAT.NC,LC_Global,DAT.NK))
    VALUES_TEMP(:,:,:,:) = DAT.NODATA

    if( process_id == master_id )then
      ! *** READ IN THE NEXT TIME STEP
      IF( DAT.IFLAG == 1 )THEN
        ! *** ASCII
        READ(DAT.IUNIT,*,IOSTAT=ISO) DAT.TIMES(2),LCC
        DAT.TIMES(2) = (DAT.TIMES(2) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days

        ALLOCATE(VALUES(DAT.NC,DAT.NK))

        IF( DAT.ITYPE > 0 )THEN
          ! *** READ ONLY SPECIFIC CELLS (FORMAT: L, I, J, VALUES)
          !READ(DAT.IUNIT,*,IOSTAT=ISO) L,II,JJ,(((DAT.VALUES(2,NCC,L,K), K=1,DAT.NK), LL=1,LCC), NCC=1,DAT.NC)
          DO LL=1,LCC
            READ(DAT.IUNIT,*,IOSTAT=ISO) L2,II,JJ,((VALUES(NCC,K), K=1,DAT.NK), NCC=1,DAT.NC)
            L = LIJ_Global(II,JJ)
            DO NCC=1,DAT.NC
              DO K=1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES(NCC,K) + DAT.VSHIFT)*DAT.VSCALE
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ! *** READ ALL CELLS
          READ(DAT.IUNIT,*,IOSTAT=ISO) (((VALUES_TEMP(2,NCC,L,K), K=1,DAT.NK), L=2,LA_Global), NCC=1,DAT.NC)
          DO NCC=1,DAT.NC
            DO L=2,LA_Global
              DO K=1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES_TEMP(2,NCC,L,K) + DAT.VSHIFT)*DAT.VSCALE
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        DEALLOCATE(VALUES)
      ELSE
        ! *** BINARY
        READ(DAT.IUNIT) DAT.TIMES(2)
        READ(DAT.IUNIT) LCC
        DAT.TIMES(2) = (DAT.TIMES(2) + DAT.TSHIFT)*DAT.TSCALE/86400.    ! Time in Julian days

        IF( DAT.ITYPE > 0 )THEN
          ! *** READ ONLY SPECIFIC CELLS (FORMAT: L, I, J, VALUES)
          !READ(DAT.IUNIT) L,II,JJ,(((DAT.VALUES(2,NCC,L,K), K=1,DAT.NK), LL=1,LCC), NCC=1,DAT.NC)
          DO LL=1,LCC
            READ(DAT.IUNIT) L2
            READ(DAT.IUNIT) II
            READ(DAT.IUNIT) JJ
            L = LIJ_Global(II,JJ)
            READ(DAT.IUNIT) ((VALUES_TEMP(2,NCC,L,K), K=1,DAT.NK), NCC=1,DAT.NC)
            DO NCC=1,DAT.NC
              DO K=1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES_TEMP(2,NCC,L,K) + DAT.VSHIFT)*DAT.VSCALE
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ! *** READ ALL CELLS
          READ(DAT.IUNIT) (((VALUES_TEMP(2,NCC,L,K), K=1,DAT.NK), L=2,LA_Global), NCC=1,DAT.NC)
          DO NCC=1,DAT.NC
            DO L=2,LA_Global
              DO K=1,DAT.NK
                VALUES_TEMP(2,NCC,L,K) = (VALUES_TEMP(2,NCC,L,K) + DAT.VSHIFT)*DAT.VSCALE
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    endif !***End on master
    
    ! *** broadcast header info just read
    Call Broadcast_Scalar(DAT.TIMES(2), master_id)
    Call Broadcast_Array(VALUES_TEMP,  master_id)

    ! *** NOW ALL PROCESSES HAVE ALL THE INPUTS.  ASSIGN LOCAL VALUES
    DO NCC=1,DAT.NC
      DO L=2,LA
        DO K=1,DAT.NK
          L2 = Map2Global(L).LG                              ! *** Global L that corresponds to Local L
          DAT.VALUES(2,NCC,L,K) = VALUES_TEMP(2,NCC,L2,K)
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE(VALUES_TEMP)
  
    DAT.NI = DAT.NI + 1
  ENDDO

    K = 1   ! *** TODO: LAYERED DATA

    IF( DAT.ITOPT == 1  )THEN
      ! *** INTERPOLATE FIELD
      I = 1
      DELTF = DAT.TIMES(I+1) - DAT.TIMES(I)
      IF( DELTF > 1.0E-9 )THEN
        FAC = (TIME - DAT.TIMES(I))/DELTF
      ELSE
        FAC = 0.5
      ENDIF
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(L)
      DO L=2,LA
        IF(DAT.VALUES(I,M,L,K) == DAT.NODATA .OR. DAT.VALUES(I+1,M,L,K) == DAT.NODATA )THEN
          RES(L) = DAT.NODATA
        ELSE
          RES(L) = DAT.VALUES(I,M,L,K) + FAC*(DAT.VALUES(I+1,M,L,K)-DAT.VALUES(I,M,L,K))
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO
    ELSE
      ! *** NO INTERPOLATION
      RES(:) = DAT.VALUES(1,M,:,K)
    ENDIF

  ! *** ADJUST FOR AREA, IF NEEDED (DEFAULT IS FLUX: IDIST == 0)
  IF( DAT.IDIST == 1 )THEN
    ! *** INPUT DATA IS A RATE
    DO L=2,LA
      IF( RES(L) /= DAT.NODATA )THEN
        RES(L) = RES(L)*DXYP(L)
      ENDIF
    ENDDO
  ELSEIF( DAT.IDIST == 2 )THEN
    ! *** INPUT DATA IS A TOTAL VOLUME
    DO L=2,LA
      IF( RES(L) /= DAT.NODATA )THEN
        RES(L) = RES(L)/DXYP(L)
      ENDIF
    ENDDO
  ENDIF

    ! UPDATE FIELD
    SELECT CASE (DAT.IUPDATE)
    CASE (1)      ! *** ADD
      DO L=1,LC
        IF( RES(L) /= DAT.NODATA )THEN
          VALUE(L) = VALUE(L) + RES(L)
        ENDIF
      ENDDO
    CASE (2)      ! *** MIN
      DO L=1,LC
        IF( RES(L) /= DAT.NODATA .AND. RES(L) < VALUE(L) )THEN
          VALUE(L) = RES(L)
        ENDIF
      ENDDO
    CASE (3)      ! *** MAX
      DO L=1,LC
        IF( RES(L) /= DAT.NODATA .AND. RES(L) > VALUE(L) )THEN
          VALUE(L) = RES(L)
        ENDIF
      ENDDO
    CASE DEFAULT  ! *** REPLACE
      DO L=1,LC
        IF( RES(L) /= DAT.NODATA )THEN
          VALUE(L) = RES(L)
        ENDIF
      ENDDO
    END SELECT
END SUBROUTINE

  SUBROUTINE UPDATETOPO(TIME,BELV,HP,HDRY)
    REAL,DIMENSION(:),INTENT(INOUT) :: BELV,HP
    REAL(RKD), INTENT(IN) :: TIME
    REAL, INTENT(IN) :: HDRY
    REAL,ALLOCATABLE,DIMENSION(:) :: WSEL
    INTEGER :: L

    ! *** Delme - This approach does not maintain mass balance.  Use with care!
    ALLOCATE(WSEL(LA))
    DO L=2,LA
      WSEL(L) = HP(L) + BELV(L)
    ENDDO
    CALL UPDATEFIELD(BATHY,TIME,1,BELV)
    DO L=2,LA
      IF(LMASKDRY(L) .AND. BATHY.CFLAG(L) > 0 ) THEN
        HP(L) = WSEL(L) - BELV(L)
        IF( HP(L) < HDRY) HP(L) = HDRY
      ENDIF
    ENDDO
    !L = 1438
    !OPEN(910,FILE=OUTDIR//'FIELD.LOG',POSITION='APPEND')
    !write(910,*) 'T=',TIME,', B=',BELV(L),', H=',HP(L),', H=',WSEL(L)
    !CLOSE(910)
    DEALLOCATE(WSEL)
END SUBROUTINE

END MODULE
