! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
MODULE MOD_NETCDF

  USE GLOBAL
  USE JULIANMOD
#ifdef NCOUT 
  USE CONVERTWGS84
  USE NETCDF
#endif
  Use Variables_MPI
  Use Variables_MPI_Mapping
  Use Variables_MPI_Write_Out
  Use Variables_MPI_Drifter
  !USE DRIFTER,ONLY:XYZSCL,MOC
  Use Mod_Map_Write_NetCDF

  IMPLICIT NONE

  TYPE NC_VARIABLE
    INTEGER(2)     :: ID = 0
    CHARACTER(15)  :: NAME,UNITS,POSITIVE,GRID_MAPPING,COORDINATES
    INTEGER(2)     :: DIM_TYPE
    REAL(4)        :: FILLVALUE
    CHARACTER(100) :: STANDARD_NAME,LONG_NAME
  END TYPE

#ifdef NCOUT 
  INTEGER(4), PARAMETER :: CNRCNT=4, TIMECNT=NF90_UNLIMITED, NC_VAR_CNT = 57
#endif
  REAL(4),    PARAMETER :: MISSING_VALUE=-999.
  INTEGER(4) :: NC_ID, NTI, NC_VAR_CNTM
  INTEGER(4) :: TIME_DIM, COL_DIM, ROW_DIM, LYR_DIM, CNR_DIM, DYE_DIM, SED_DIM, SND_DIM, TOX_DIM, CWQ_DIM, LPT_DIM, KB_DIM, SXD_DIM
  INTEGER(4) :: COLCNT, ROWCNT, LYRCNT, NCDFTYPE
  INTEGER(4) :: NC_XC, NC_YC, NC_XB, NC_YB, NC_CRS, NC_SIG, NC_KC, NC_TIM
  INTEGER(4) :: NC_TRK_XC,NC_TRK_YC,NC_TRK_CRS,NC_TRK_SIG,NC_TRK_VOL

  CHARACTER(10)                  :: NC_DATESTAMP
  INTEGER(4), ALLOCATABLE        :: NC_IDX(:)
  TYPE(NC_VARIABLE), ALLOCATABLE :: NC_VAR(:)
  REAL(4),ALLOCATABLE            :: TRK_LON(:),TRK_LAT(:),TRK_LEV(:),TRK_VOL(:)
  DATA NCDFTYPE /0/

  CONTAINS

#ifdef NCOUT 

    SUBROUTINE NC_DEFINE_VARS(NC_ID,IDX)
  
      INTEGER :: NC_ID,IDX,ID,STAT
    
      IF( NCDFTYPE==0 )THEN
        NCDFTYPE = 1
        ALLOCATE(NC_VAR(NC_VAR_CNT),NC_IDX(NC_VAR_CNT))
        ALLOCATE(TRK_LON(NPD),TRK_LAT(NPD),TRK_LEV(NPD),TRK_VOL(NPD))

        NC_IDX = 0
        
        NC_VAR(1:40) = (/&
          NC_VARIABLE( 1,'Bottom','m',  'down',  'crs','lat lng',2,MISSING_VALUE,'','Bottom Bathymetry'),&
          NC_VARIABLE( 2,'WSEL',  'm',  'up',  'crs','lat lng',3,MISSING_VALUE,'','Water Surface Elevation'),&
          NC_VARIABLE( 3,'Vx',    'm/s','',    'crs','lat lng',4,MISSING_VALUE,'','Eastward Water Velocity'),&
          NC_VARIABLE( 4,'Vy',    'm/s','',    'crs','lat lng',4,MISSING_VALUE,'','Northward Water Velocity'),&
          NC_VARIABLE( 5,'Vz',    'm/s','',  'crs','lat lng',4,MISSING_VALUE,'','Upward Water Velocity'),&
          NC_VARIABLE( 6,'TauB',  'N/m2','',   'crs','lat lng',3,MISSING_VALUE,'','Bottom Shear Stress'),&
          NC_VARIABLE( 7,'Wx',    'm/s','',    'crs','lat lng',3,MISSING_VALUE,'','Eastward Wind Velocity'),&
          NC_VARIABLE( 8,'Wy',    'm/s','',    'crs','lat lng',3,MISSING_VALUE,'','Northward Wind Velocity'),&
          NC_VARIABLE( 9,'Hs',    'm','',     'crs','lat lng',3,MISSING_VALUE,'','Significant Wave Height'),&
          NC_VARIABLE(10,'Dp',    'degree','','crs','lat lng',3,MISSING_VALUE,'','Wave Direction'),&
          NC_VARIABLE(11,'Tp',    's','',     'crs','lat lng',3,MISSING_VALUE,'','Wave Period'),&      
          NC_VARIABLE(12,'Salt',  'ppt','',    'crs','lat lng',4,MISSING_VALUE,'','Salinity'),&
          NC_VARIABLE(13,'Temp',  'C','', 'crs','lat lng',4,MISSING_VALUE,'','Temperature'),&
          NC_VARIABLE(14,'Dye',   'mg/l','',    'crs','lat lng',10,MISSING_VALUE,'','Dye'),&
          NC_VARIABLE(15,'SFL',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Shellfish Larvae'),&
          NC_VARIABLE(16,'Toxic', 'mg/l','',    'crs','lat lng',7,MISSING_VALUE,'','Toxics'),&
          NC_VARIABLE(17,'SED',   'mg/l','',    'crs','lat lng',5,MISSING_VALUE,'','Cohesive Sediment'),&
          NC_VARIABLE(18,'SND',   'mg/l','',    'crs','lat lng',6,MISSING_VALUE,'','Non-cohesive Sediment'),&
          NC_VARIABLE(19,'CHC',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Cyanobacteria'),&
          NC_VARIABLE(20,'CHD',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Diatom Algae'),&
          NC_VARIABLE(21,'CHG',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Green Algae'),&
          NC_VARIABLE(22,'ROC',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Refractory Particulate Organic Carbon'),&
          NC_VARIABLE(23,'LOC',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Labile Particulate Organic Carbon'),&
          NC_VARIABLE(24,'DOC',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Organic Carbon'),&
          NC_VARIABLE(25,'ROP',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Refractory Particulate Organic Phosphorus'),&
          NC_VARIABLE(26,'LOP',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Labile Particulate Organic Phosphorus'),&
          NC_VARIABLE(27,'DOP',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Organic Phosphorus'),&
          NC_VARIABLE(28,'P4D',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Total Phosphate'),&
          NC_VARIABLE(29,'RON',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Refractory Particulate Organic Nitrogen'),&
          NC_VARIABLE(30,'LON',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Labile Particulate Organic Nitrogen'),&
          NC_VARIABLE(31,'DON',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Organic Nitrogen'),&
          NC_VARIABLE(32,'NHX',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Ammonia Nitrogen'),&
          NC_VARIABLE(33,'NOX',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Nitrate Nitrogen'),&
          NC_VARIABLE(34,'SUU',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Particulate Biogenic Silica'),&
          NC_VARIABLE(35,'SAA',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Available Silica'),&
          NC_VARIABLE(36,'COD',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Chemical Oxygen Demand'),&
          NC_VARIABLE(37,'DOX',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Oxygen'),&
          NC_VARIABLE(38,'TAM',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Total Active Metal'),&
          NC_VARIABLE(39,'FCB',   'mpn/100ml','','crs','lat lng',4,MISSING_VALUE,'','Fecal Coliform Bacteria'),&
          NC_VARIABLE(40,'MAC',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Macroalgae') /)
      
        IF(  ISPD > 0  )THEN
          IF( ANY(ISOILSPI == 1)) NC_VAR(41) = NC_VARIABLE(41,'MOC', 'kg','',  'crs','lat lng',3,MISSING_VALUE,'','Oil Mass')
        ENDIF
                  
        IF( LSEDZLJ )THEN
          STAT = 11
        ELSE
          STAT = 3
        ENDIF
        NC_VAR(42) = NC_VARIABLE(42,'bed_top',  'layer',  '',  'crs','lat lng',STAT,MISSING_VALUE,'','Sediment Bed Top Layer')
        NC_VAR(43) = NC_VARIABLE(43,'bed_thickness',  'm',  '',  'crs','lat lng',11,MISSING_VALUE,'','Sediment Bed Layer Thickness')
        NC_VAR(44) = NC_VARIABLE(44,'bed_wet_density',  'kg/m3',  '',  'crs','lat lng',11,MISSING_VALUE,'','Sediment Bed Wet Density')
        NC_VAR(45) = NC_VARIABLE(45,'bed_porosity',  '-',  '',  'crs','lat lng',11,MISSING_VALUE,'','Sediment Bed Porosity')
        IF( LSEDZLJ )THEN
          STAT = 11
        ELSE
          STAT = 12
        ENDIF
        NC_VAR(46) = NC_VARIABLE(46,'bed_mass',  'g/m2',  '',  'crs','lat lng',STAT,MISSING_VALUE,'','Sediment Bed Mass')        
        NC_VAR(47) = NC_VARIABLE(47,'bed_tau',  'dynes/cm^2',  '',  'crs','lat lng',3,MISSING_VALUE,'','Shear Stress')
        NC_VAR(48) = NC_VARIABLE(48,'bed_dry_density',  'g/cm^3',  '',  'crs','lat lng',11,MISSING_VALUE,'','Dry Bulk Density')
        NC_VAR(49) = NC_VARIABLE(49,'bed_mass_percentage',  '-',  '',  'crs','lat lng',12,MISSING_VALUE,'','Sediment Bed Mass Percentage')
        NC_VAR(50) = NC_VARIABLE(50,'bed_d50',  'microns',  '',  'crs','lat lng',3,MISSING_VALUE,'','Particle Size of Bed Surface')
        NC_VAR(51) = NC_VARIABLE(51,'bed_erosion_rate',  'g/cm^2/s',  '',  'crs','lat lng',3,MISSING_VALUE,'','Total Erosion Rate')
        NC_VAR(52) = NC_VARIABLE(52,'bed_deposition_rate',  'g/cm^2/s',  '',  'crs','lat lng',3,MISSING_VALUE,'','Total Deposition Rate')
        NC_VAR(53) = NC_VARIABLE(53,'bedload_mass',  'g/m^2',  '',  'crs','lat lng',13,MISSING_VALUE,'','Bedload Mass')
        NC_VAR(54) = NC_VARIABLE(54,'bedload_x_flux',  'g/s',  '',  'crs','lat lng',13,MISSING_VALUE,'','Bedload Flux in X-direction')
        NC_VAR(55) = NC_VARIABLE(55,'bedload_y_flux',  'g/s',  '',  'crs','lat lng',13,MISSING_VALUE,'','Bedload Flux in Y-direction')
        
        NC_VAR(56) = NC_VARIABLE(56,'bed_toxics',  'mg/m^2',  '',  'crs','lat lng',14,MISSING_VALUE,'','Sediment Bed Toxics')
        NC_VAR(57) = NC_VARIABLE(57,'bedload_toxics',  'mg/m^2',  '',  'crs','lat lng',15,MISSING_VALUE,'','Bedload Toxics')
      ENDIF
  
      IF( NC_VAR(IDX).ID > 0 )THEN
        SELECT CASE (NC_VAR(IDX).DIM_TYPE)
          CASE (1)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/TIME_DIM/),ID)
          CASE (2)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/ROW_DIM,COL_DIM/),ID)
          CASE (3)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (4)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/LYR_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (5)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/SED_DIM,LYR_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (6)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/SND_DIM,LYR_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (7)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/TOX_DIM,LYR_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (8)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/CWQ_DIM,LYR_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID) 
          CASE (9)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/LPT_DIM,TIME_DIM/),ID)      
          CASE (10)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/DYE_DIM,LYR_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (11)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/KB_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (12)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/SXD_DIM,KB_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (13)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/SXD_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (14)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/TOX_DIM,KB_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
          CASE (15)
            STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/TOX_DIM,ROW_DIM,COL_DIM,TIME_DIM/),ID)
        END SELECT           
            
      STAT=CHECK_ERR(STAT,'def_'//TRIM(NC_VAR(IDX).NAME))
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,ID,1,1,DEFLEV)
      IF( LEN_TRIM(NC_VAR(IDX).STANDARD_NAME) > 0 )THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'standard_name',TRIM(NC_VAR(IDX).STANDARD_NAME))
      ENDIF
      IF( LEN_TRIM(NC_VAR(IDX).LONG_NAME) > 0 )THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'long_name',TRIM(NC_VAR(IDX).LONG_NAME))
      ENDIF
      IF( LEN_TRIM(NC_VAR(IDX).UNITS) > 0 )THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'units',TRIM(NC_VAR(IDX).UNITS))
      ENDIF
      IF( LEN_TRIM(NC_VAR(IDX).POSITIVE) > 0 )THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'positive',TRIM(NC_VAR(IDX).POSITIVE))
        STAT=NF90_PUT_ATT(NC_ID,ID,'_coordinatezispositive',TRIM(NC_VAR(IDX).POSITIVE))
      ENDIF
      IF( LEN_TRIM(NC_VAR(IDX).GRID_MAPPING) > 0 )THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'grid_mapping',TRIM(NC_VAR(IDX).GRID_MAPPING))
      ENDIF
      IF( LEN_TRIM(NC_VAR(IDX).COORDINATES) > 0 )THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'coordinates',TRIM(NC_VAR(IDX).COORDINATES))
      ENDIF
      STAT=NF90_PUT_ATT(NC_ID,ID,'_fillvalue',NC_VAR(IDX).FILLVALUE)
      NC_VAR(IDX).ID = ID
      NC_IDX(IDX) = ID
    ENDIF
    
    END SUBROUTINE

    SUBROUTINE NC_SET_VARS(NC_ID)
      INTEGER :: NC_ID,K
  
      ! ** ALWAYS EXPORT: ZBOT,WSEL,VELX,VELY,VELZ
      DO K=1,5
        CALL NC_DEFINE_VARS(NC_ID,K)                    ! K=1:5 
      ENDDO
    
      ! ** OPTIONALLY EXPORT
      IF( ISNCDF(10) >= 1) CALL NC_DEFINE_VARS(NC_ID,6) ! 6:TAUB
    
      IF( ISNCDF(11) >= 1 )THEN
        CALL NC_DEFINE_VARS(NC_ID,7)                    ! 7:WINX
        CALL NC_DEFINE_VARS(NC_ID,8)                    ! 8:WINY
      ENDIF 

    IF( ISNCDF(12) >= 1 )THEN
      DO K=9,11
        CALL NC_DEFINE_VARS(NC_ID,K)                  ! WINDWAVE: 9:11 WHEI,WANG,WPER
      ENDDO
    ENDIF
    
    DO K=1,7
      IF( ISNCDF(K) >= 1) CALL NC_DEFINE_VARS(NC_ID,  11+K)  ! 12-18: SAL,TEM,...SND
    ENDDO
    
    DO K=1,NWQV
      IF(ISTRWQ(K) > 0 .AND. ISNCDF(8) >= 1) CALL NC_DEFINE_VARS(NC_ID, 18 + K)   ! 19-40: WQV
    ENDDO
    
    ! ** OIL LPT
    IF( ISNCDF(9) >= 1 )THEN  
      IF( ANY(ISOILSPI == 1)) CALL NC_DEFINE_VARS(NC_ID, 41)  
    ENDIF
    
    IF( LSEDZLJ )THEN
      CALL NC_DEFINE_VARS(NC_ID,42)
      DO K=46,52
        CALL NC_DEFINE_VARS(NC_ID,K)
      ENDDO
      IF( NCALC_BL > 0 )THEN
        DO K=53,55
          CALL NC_DEFINE_VARS(NC_ID,K)
        ENDDO
      ENDIF
        
      IF( ISTRAN(5) > 0 )THEN
        CALL NC_DEFINE_VARS(NC_ID,56)
        IF( NCALC_BL > 0 )THEN
          CALL NC_DEFINE_VARS(NC_ID,57)
        ENDIF
      ENDIF
    ELSEIF( ISBEXP >= 1 .AND. KB > 1 .AND. ( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) )THEN
      DO K=42,45
        CALL NC_DEFINE_VARS(NC_ID,K)
      ENDDO
      IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
        CALL NC_DEFINE_VARS(NC_ID,46)
      ENDIF
      IF( ISTRAN(5) >= 1 )THEN
        CALL NC_DEFINE_VARS(NC_ID,56)
      ENDIF
    ENDIF

  END SUBROUTINE

    INTEGER FUNCTION NC_CREATE_FILE(NC_ID, FILENAME)
      ! ** VARIABLE DECLARATIONS
      CHARACTER(*),INTENT(IN) :: FILENAME
      INTEGER,INTENT(OUT) :: NC_ID
      CHARACTER(20) :: BDATE,ZONESTR*3
      INTEGER :: STAT,NSXD
    
      ! ** ENTER DEFINE MODE
      STAT=CHECK_ERR(NF90_CREATE(FILENAME, NF90_CLOBBER+NF90_NETCDF4, NC_ID),'nc_create');!NF90_HDF5
      IF(STAT /= NF90_NOERR )THEN
          NC_CREATE_FILE = STAT
          RETURN
      ENDIF
      ! *** modified to handle global
      COLCNT = IMX_Global - IMN_Global +1 
      ROWCNT = JMX_Global - JMN_Global +1 
      LYRCNT = KC
      BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)
      IF( LSEDZLJ )THEN
        NSXD = NSCM
      ELSE
        NSXD = NSED + NSND    !NSCM
      ENDIF
          
      CALL UTMPARS

      ! ** DEFINE DIMENSIONS
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'time', TIMECNT, TIME_DIM),'time_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'imax', COLCNT, COL_DIM),'col_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'jmax', ROWCNT, ROW_DIM),'row_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'kmax', LYRCNT, LYR_DIM),'lyr_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'cnr', CNRCNT, CNR_DIM),'cnr_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'ndye', NDYE, DYE_DIM),'dye_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'ntox', NTOX, TOX_DIM),'tox_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nsed', NSED, SED_DIM),'sed_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nsnd', NSND, SND_DIM),'snd_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nsxd', NSXD, SXD_DIM),'sxd_dim');      
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nwqv', NWQV, CWQ_DIM),'cwq_dim');
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'npd', NPD, LPT_DIM),'lpt_dim');   
      STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'kb', KB, KB_DIM),'kb_dim');             
    
      ! ** DEFINE VARIABLES
      STAT=CHECK_ERR(NF90_DEF_VAR(NCID=NC_ID,NAME='crs',XTYPE=NF90_CHAR,VARID=NC_CRS),'def_crs') !DIMIDS=0,
      STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'grid_mapping_name','latitude_longitude')
      STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'longitude_of_prime_meridian',0.)
      STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'semi_major_axis',R_MJR)
      STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'inverse_flattening',1._8/FLA)

      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lng',NF90_FLOAT,(/ROW_DIM,COL_DIM/),NC_XC),'def_lng')
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_XC,1,1,DEFLEV)
      STAT=NF90_PUT_ATT(NC_ID,NC_XC,'standard_name','longitude')
      STAT=NF90_PUT_ATT(NC_ID,NC_XC,'long_name','longitude at grid cell center')
      STAT=NF90_PUT_ATT(NC_ID,NC_XC,'units','degrees_east')
      STAT=NF90_PUT_ATT(NC_ID,NC_XC,'bounds','lng_bnds')
      STAT=NF90_PUT_ATT(NC_ID,NC_XC,'projection','geographic')
      STAT=NF90_PUT_ATT(NC_ID,NC_XC,'_FillValue',MISSING_VALUE)

      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lat',NF90_FLOAT,(/ROW_DIM,COL_DIM/),NC_YC),'def_lat')
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_YC,1,1,DEFLEV)
      STAT=NF90_PUT_ATT(NC_ID,NC_YC,'standard_name','latitude')
      STAT=NF90_PUT_ATT(NC_ID,NC_YC,'long_name','latitude at grid cell center')
      STAT=NF90_PUT_ATT(NC_ID,NC_YC,'units','degrees_north')
      STAT=NF90_PUT_ATT(NC_ID,NC_YC,'bounds','lat_bnds')
      STAT=NF90_PUT_ATT(NC_ID,NC_YC,'projection','geographic')
      STAT=NF90_PUT_ATT(NC_ID,NC_YC,'_FillValue',MISSING_VALUE)

      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lng_bnds',NF90_FLOAT,(/CNR_DIM,ROW_DIM,COL_DIM/),NC_XB),'def_xb')
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_XB,1,1,DEFLEV)
      STAT=NF90_PUT_ATT(NC_ID,NC_XB,'_FillValue',MISSING_VALUE)

      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lat_bnds',NF90_FLOAT,(/CNR_DIM,ROW_DIM,COL_DIM/),NC_YB),'def_yb')
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_YB,1,1,DEFLEV)
      STAT=NF90_PUT_ATT(NC_ID,NC_YB,'_FillValue',MISSING_VALUE)

      IF(IGRIDV > 0 )THEN
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'layers',NF90_BYTE,(/ROW_DIM,COL_DIM/),NC_KC),'def_lat')
        STAT=NF90_PUT_ATT(NC_ID,NC_KC,'long_name','number of vertical layers')
      
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'sigma',NF90_FLOAT,(/LYR_DIM,ROW_DIM,COL_DIM/), NC_SIG),'def_sigma')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'standard_name','ocean_sigma_coordinate')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'long_name','sigma at layer midpoints')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'units','sigma_level')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'positive','up')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'formula_terms','sigma: sigma eta: WSEL depth: Bottom')
      ELSE
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'sigma',NF90_FLOAT,(/LYR_DIM/), NC_SIG),'def_sigma')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'standard_name','ocean_sigma_coordinate')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'long_name','sigma at layer midpoints')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'units','sigma_level')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'positive','up')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'formula_terms','sigma: sigma eta: WSEL depth: Bottom')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateZisPositive','up')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateTransformType','Vertical')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateAxisType','GeoZ')
        STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateAxes','sigma')
      ENDIF
    
      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'time',NF90_DOUBLE,(/TIME_DIM/), NC_TIM),'def_time')
      STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'standard_name','time')
      STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'long_name','time')
      STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'units','days since '//TRIM(BDATE))
      STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'calendar','julian')
      STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'axis','T')
    
      ! *** LPT
      IF( ISNCDF(9) >= 1 )THEN 
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_LAT',NF90_FLOAT,(/LPT_DIM,TIME_DIM/),NC_TRK_YC),'def_trk_lat')
        STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_TRK_YC,1,1,DEFLEV)
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'standard_name','')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'long_name','Drifter Latitude')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'units','degrees_north')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'valid_min','-90')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'valid_max','90')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'grid_mapping','crs')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'_FillValue',MISSING_VALUE)
    
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_LON',NF90_FLOAT,(/LPT_DIM,TIME_DIM/),NC_TRK_XC),'def_trk_lon')
        STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_TRK_XC,1,1,DEFLEV)
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'standard_name','')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'long_name','Drifter Longitude')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'units','degrees_east')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'valid_min','-180')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'valid_max','180')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'grid_mapping','crs')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'_FillValue',MISSING_VALUE)

        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_LEV',NF90_FLOAT,(/LPT_DIM,TIME_DIM/), NC_TRK_SIG),'def_trk_sigma')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'standard_name','')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'long_name','Drifter Elevation')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'units','m')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'grid_mapping','crs')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'coordinates','TRK_LAT TRK_LON')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateZisPositive','up')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateTransformType','Vertical')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateAxisType','GeoZ')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateAxes','lev')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_FillValue',MISSING_VALUE)

      IF( ANY(ISOILSPI == 1) )THEN
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_VOL',NF90_FLOAT,(/LPT_DIM,TIME_DIM/), NC_TRK_VOL),'def_trk_vol')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'standard_name','')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'long_name','Drifter Oil Volume')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'units','m3')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'grid_mapping','crs')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'coordinates','TRK_LAT TRK_LON')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'_FillValue',MISSING_VALUE)
      ENDIF
    ENDIF
    ! *** END LPT

      CALL NC_SET_VARS(NC_ID)

      ! ** ASSIGN GLOBAL ATTRIBUTES

      WRITE(ZONESTR,'(I3)') ABS(UTMZ)
      STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'Conventions','CF-1.4'),'att_cf')
      STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'Base_date',BDATE),'att_base')
      STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'Project',TRIM(PROJ)),'att_prj')
      IF( UTMZ > 0 )THEN
          STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'utm_zone','UTM Zone '//TRIM(ZONESTR)//' Northern Hemisphere'),'')
      ELSE
          STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'utm_zone','UTM Zone '//TRIM(ZONESTR)//' Southern Hemisphere'),'')
      ENDIF

      ! ** LEAVE DEFINE MODE
      STAT=CHECK_ERR(NF90_ENDDEF(NC_ID),'end_def');
      NC_CREATE_FILE = STAT
    END FUNCTION
  
    SUBROUTINE NC_CLOSE_FILE(NC_ID)
      INTEGER NC_ID,STATUS
      IF( NCDFOUT > 1 )THEN
          STATUS=CHECK_ERR(NF90_CLOSE(NC_ID),'close')
          IF(STATUS /= NF90_NOERR )THEN
              PRINT * ,'Cannot close nc file!'
              RETURN
          ENDIF
      ENDIF
    END SUBROUTINE

    INTEGER FUNCTION CHECK_ERR(STATUS,MSG)
      INTEGER :: STATUS
      CHARACTER(*) :: MSG
      CHECK_ERR = STATUS
      IF( STATUS .NE. NF90_NOERR )THEN
          PRINT *, MSG//': '//TRIM(NF90_STRERROR(STATUS))
          CALL STOPP('.')
      END IF
    END FUNCTION 

    INTEGER FUNCTION NC_WRITE_SCALAR1D(NC_ID,ID,SNAPSHOT,FIELD1D)
      ! ** WRITES 1D SCALAR FIELD TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT
      REAL(4),INTENT(IN) :: FIELD1D(:)
      INTEGER :: STAT,STR1D(2),CNT1D(2)

      STR1D = (/1, SNAPSHOT/)
      CNT1D = (/NPD,1/)
      STAT=NF90_PUT_VAR(NC_ID, ID, FIELD1D, STR1D, CNT1D)
      NC_WRITE_SCALAR1D=STAT
    END FUNCTION
  
    INTEGER FUNCTION NC_WRITE_SCALAR2D(NC_ID,ID,FIELD2D,SNAPSHOT)
      ! ** WRITES 2D SCALAR FIELD TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT
      REAL(4),INTENT(IN) :: FIELD2D(:,:)
      INTEGER :: STAT,STR2D(3),CNT2D(3)

      STR2D = (/1, 1, SNAPSHOT/)
      CNT2D = (/ROWCNT,COLCNT,1/)
      STAT=NF90_PUT_VAR(NC_ID, ID, FIELD2D, STR2D, CNT2D)
      NC_WRITE_SCALAR2D=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_SCALAR3D(NC_ID,ID,FIELD3D,SNAPSHOT,NK)
      ! ** WRITES 3D SCALAR FIELD TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT,NK
      REAL(4),INTENT(IN) :: FIELD3D(:,:,:)
      INTEGER :: STAT,STR3D(4),CNT3D(4)
      STR3D = (/1, 1, 1, SNAPSHOT/)
      CNT3D = (/NK,ROWCNT,COLCNT,1/)
      STAT=NF90_PUT_VAR(NC_ID, ID, FIELD3D, STR3D, CNT3D)
      NC_WRITE_SCALAR3D=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_CLASS3D(NC_ID,ID,FIELD3D,SNAPSHOT,CLSCNT,NK)
      ! ** WRITES 3D SCALAR FIELD WITH CLASSES TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT,CLSCNT,NK
      REAL(4),INTENT(IN) :: FIELD3D(:,:,:,:)
      INTEGER :: STAT,STR3D(5),CNT3D(5)
      STR3D = (/1, 1, 1, 1, SNAPSHOT/)
      CNT3D = (/CLSCNT,NK,ROWCNT,COLCNT,1/)
      STAT=NF90_PUT_VAR(NC_ID, ID, FIELD3D, STR3D, CNT3D)
      NC_WRITE_CLASS3D=STAT
    END FUNCTION

    SUBROUTINE READCORN
      ! ** USE GLOBAL ,ONLY:XCR,YCR,UCOR,IC,JC
      USE INFOMOD,ONLY:READSTR
      IMPLICIT NONE
      INTEGER(4)::I,J,M
      CHARACTER*80 :: STR*200
    
      OPEN(UCOR,FILE='corners.inp',ACTION='READ')
      STR = READSTR(UCOR)
      ALLOCATE(XCR(4,JC_global,IC_global),YCR(4,JC_global,IC_global))
      XCR = MISSING_VALUE
      YCR = MISSING_VALUE
      DO WHILE(1)
          READ(UCOR,*,END=100,ERR=998) I,J,(XCR(M,J,I),YCR(M,J,I),M=1,4)
      ENDDO
      100 CLOSE(UCOR)
    
      RETURN
      998 CLOSE(UCOR)
      CALL STOPP('CORNERS.INP READING ERROR!')
    END SUBROUTINE
    
    INTEGER FUNCTION NC_WRITE_GRID(NC_ID,GRD_X,GRD_Y)
      ! ** WRITE HORIZONTAL GRID TO NETCDF FILE
      REAL,INTENT(IN) :: GRD_X(:,:,:),GRD_Y(:,:,:)
      REAL(4) :: LON(4,ROWCNT,COLCNT), LAT(4,ROWCNT,COLCNT)
      REAL(4) :: XC(ROWCNT,COLCNT),YC(ROWCNT,COLCNT)
      REAL(8)   :: XM(1),YM(1),XLL(1),YLL(1)
      INTEGER,INTENT(IN) :: NC_ID
      INTEGER :: STAT,I,J,M

      LON = MISSING_VALUE
      LAT = MISSING_VALUE
      XC  = MISSING_VALUE
      YC  = MISSING_VALUE

      DO I=1,COLCNT
          DO J=1,ROWCNT
              IF( GRD_X(1,J,I) > 0 )THEN
                  XC(J,I)=0.
                  YC(J,I)=0.
                  DO M=1,CNRCNT
                      XM(1)=GRD_X(M,J,I)
                      YM(1)=GRD_Y(M,J,I)
                      CALL UTMR_WGS84(XM,YM,XLL,YLL)
                      LON(M,J,I) = XLL(1)
                      LAT(M,J,I) = YLL(1)
                      XC(J,I)=XC(J,I)+LON(M,J,I)
                      YC(J,I)=YC(J,I)+LAT(M,J,I)
                  ENDDO
                  XC(J,I)=XC(J,I)/CNRCNT
                  YC(J,I)=YC(J,I)/CNRCNT
              ENDIF
          ENDDO
      ENDDO

      STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_XB,LON),'xb_put') 
      STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_YB,LAT),'yb_put') 
      STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_XC,XC),'xc_put') 
      STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_YC,YC),'yc_put') 
      NC_WRITE_GRID=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_ZBOT(NC_ID)
      ! ** WRITE BOTTOM DEPTH AND SIGMA LEVEL TO A NETCDF FILE
      ! 2017-09-07, NTL: Updated NetCDF layers from top=1 to bottom = KC as the
      !                  ocean_sigma_coordinate ranges from 0 at surface to -1 at bottom 
      INTEGER NC_ID,STAT,I,J,K,L,M
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: BATHY
      REAL(4) :: SIGMA(LYRCNT), SUMDZC
      BYTE, ALLOCATABLE, DIMENSION(:,:) :: KCSGZ
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: DZSGZ

      IF(IGRIDV > 0 )THEN ! Sigma-Zed level
          ALLOCATE(KCSGZ(ROWCNT,COLCNT),DZSGZ(LYRCNT,ROWCNT,COLCNT)) 
          KCSGZ = KMINV
          DZSGZ = MISSING_VALUE
          DO L = 2,LA_Global ! *** Modified to handle global
              I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
              J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
              KCSGZ(J,I) = KC - KSZ_Global(L) + 1    ! Number of vertical layers
              SUMDZC = 0.0
              DO M=1,LYRCNT
                  K = LYRCNT - M + 1
                  DZSGZ(M,J,I) = -(SUMDZC + 0.5*DZC(L,K))
                  SUMDZC = SUMDZC + DZC(L,K)
              ENDDO
          ENDDO     
          STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_KC,KCSGZ),'put_kc')
          STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_SIG,DZSGZ),'put_dz')
          DEALLOCATE(KCSGZ,DZSGZ)
      ELSE  ! ** Standard sigma stretched level
          DO M=1,LYRCNT
              K = LYRCNT - M + 1
              SIGMA(M) = SUM(DZCK(1:K)) - 0.5*DZCK(K) - 1.0
          ENDDO
          STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_SIG,SIGMA),'put_sig')
      ENDIF

      ! ** Write bottom bathymetry
      ALLOCATE(BATHY(ROWCNT,COLCNT))
      BATHY = MISSING_VALUE
      DO L = 2,LA_Global ! *** Modified to handle global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          BATHY(J,I) = - BELV_Global(L)
      ENDDO     
      STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_IDX(1),BATHY),'put_zb')
      DEALLOCATE(BATHY)
      NC_WRITE_ZBOT=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_WSEL(NC_ID, SNAPSHOT)
      ! ** WRITE WATER SURFACE ELEVATION
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: WSEL
      INTEGER :: STAT, I,J,L
  
      ALLOCATE(WSEL(ROWCNT,COLCNT))
      WSEL = MISSING_VALUE
      DO L = 2,LA_Global ! *** Modified to handle global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              WSEL(J,I) = BELV_Global(L) + HP_Global(L)
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(2),WSEL,SNAPSHOT)
      DEALLOCATE(WSEL)
      NC_WRITE_WSEL=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_VEL(NC_ID, SNAPSHOT)
      ! ** WRITES VELOCITIES TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: VX, VY, VZ
      REAL(4) :: UTMP,VTMP
      INTEGER :: I,J,L,K,M,LN,STAT

      ALLOCATE(VX(LYRCNT,ROWCNT,COLCNT),VY(LYRCNT,ROWCNT,COLCNT),VZ(LYRCNT,ROWCNT,COLCNT))
      VX = MISSING_VALUE
      VY = MISSING_VALUE
      VZ = MISSING_VALUE
      DO L = 2,LA_GLobal
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_GLobal(L) > HDRY )THEN  
              IF( ROTA == 1 )THEN
                  LN = LNC_GLobal(L)  
                  DO M=1,KC
                      K = KC - M + 1
                      UTMP = 0.5*(RSSBCE_GLobal(L)*U_GLobal(L+1,K)+RSSBCW_GLobal(L)*U_GLobal(L,K))  
                      VTMP = 0.5*(RSSBCN_GLobal(L)*V_GLobal(LN ,K)+RSSBCS_GLobal(L)*V_GLobal(L,K))  
                      VX(M,J,I) = CUE_Global(L)*UTMP+CVE_Global(L)*VTMP  
                      VY(M,J,I) = CUN_Global(L)*UTMP+CVN_Global(L)*VTMP  
                      VZ(M,J,I) = W_Global(L,K)
                  ENDDO
              ELSE
                  DO M=1,KC
                      K = KC - M + 1
                      VX(M,J,I) = U_Global(L,K)
                      VY(M,J,I) = V_Global(L,K)
                      VZ(M,J,I) = W_Global(L,K)
                  ENDDO
              ENDIF
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR3D(NC_ID, NC_IDX(3), VX, SNAPSHOT, LYRCNT)
      STAT=NC_WRITE_SCALAR3D(NC_ID, NC_IDX(4), VY, SNAPSHOT, LYRCNT)
      STAT=NC_WRITE_SCALAR3D(NC_ID, NC_IDX(5), VZ, SNAPSHOT, LYRCNT)
      DEALLOCATE(VX, VY, VZ)
      NC_WRITE_VEL=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_SHEAR(NC_ID, SNAPSHOT)
    ! ** WRITES BOTTOM SHEAR STRESSES TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: SHEAR
    INTEGER :: I,J,L,STAT
  
    ALLOCATE(SHEAR(ROWCNT,COLCNT))
    SHEAR = MISSING_VALUE
    DO L = 2,LA_GLobal
      I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
      J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
      IF( HP_Global(L) > HDRY )THEN  
        IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
          IF( LSEDZLJ )THEN
            SHEAR(J,I) = TAU_GLobal(L)*0.1             ! DYN/CM2 TO N/M2
          ELSEIF( ISBEDSTR >= 1 )THEN
            SHEAR(J,I) = 1000.*TAUBSED_Global(L)       ! M2/S2 TO N/M2
            IF( ISBEDSTR == 1 )THEN
              SHEAR(J,I) = 1000.*TAUBSND_Global(L)
            ENDIF
          ELSE
            SHEAR(J,I) = 1000.*TAUB(L)
          ENDIF
        ELSE
          SHEAR(J,I) = 1000.*MAX(QQ_Global(L,0),QQMIN)/CTURB2
        ENDIF
      ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(6),SHEAR,SNAPSHOT)
    DEALLOCATE(SHEAR)
    NC_WRITE_SHEAR=STAT
  END FUNCTION

    INTEGER FUNCTION NC_WRITE_WIND(NC_ID, SNAPSHOT)
      ! ** WRITES WIND VELOCITIES TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: WX, WY
      INTEGER :: I,J,L,STAT

      ALLOCATE(WX(ROWCNT,COLCNT),WY(ROWCNT,COLCNT))
      WX = MISSING_VALUE
      WY = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          WX(J,I) = WNDVELE(L)  
          WY(J,I) = WNDVELN(L)
      ENDDO
      STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(7), WX, SNAPSHOT)
      STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(8), WY, SNAPSHOT)
      DEALLOCATE(WX, WY)
      NC_WRITE_WIND=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_WAVE(NC_ID, SNAPSHOT)
      ! ** WRITES WAVE PARAMETERS TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: HS,TP,DIR
      INTEGER :: STAT, I,J,L
  
      ALLOCATE(HS(ROWCNT,COLCNT),TP(ROWCNT,COLCNT),DIR(ROWCNT,COLCNT))
      HS = MISSING_VALUE
      TP = MISSING_VALUE
      DIR = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              HS(J,I) = WV_GLobal(L).HEIGHT
              TP(J,I) = WV_GLobal(L).FREQ
              DIR(J,I) = WV_GLobal(L).DIR
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(9),HS,SNAPSHOT)
      STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(10),DIR,SNAPSHOT)
      STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(11),TP,SNAPSHOT)

      DEALLOCATE(HS,TP,DIR)
      NC_WRITE_WAVE=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_CONS(NC_ID, SNAPSHOT, IDX, VALUES, NK, IK)
      ! ** WRITES SAL/TEM TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, NK, IK
      REAL,INTENT(IN) :: VALUES(:,:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: CONC
      INTEGER :: I,J,L,K,M,STAT

      ALLOCATE(CONC(NK,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              DO M=1,NK
                  IF( IK == 1 )THEN
                    K = NK - M + 1
                  ELSE
                    K = M
                  ENDIF
                  CONC(M,J,I) = VALUES(L,K)
              ENDDO
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT, NK)
      DEALLOCATE(CONC)
      NC_WRITE_CONS=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_SZLJ_CONS(NC_ID, SNAPSHOT, IDX, VALUES, NK)
      ! ** WRITES SAL/TEM TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, NK
      REAL(RKD),INTENT(IN) :: VALUES(:,:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: CONC
      INTEGER :: I,J,L,K,STAT

      ALLOCATE(CONC(NK,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              DO K=1,NK
                  CONC(K,J,I) = VALUES(K,L)
              ENDDO
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT, NK)
      DEALLOCATE(CONC)
      NC_WRITE_SZLJ_CONS=STAT
    END FUNCTION    
    INTEGER FUNCTION NC_WRITE_SZLJ_BEDTOP(NC_ID, SNAPSHOT, IDX, VALUES, NK)
      ! ** WRITES SAL/TEM TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, NK
      INTEGER,INTENT(IN) :: VALUES(:,:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: CONC
      INTEGER :: I,J,L,K,STAT

      ALLOCATE(CONC(NK,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              DO K=1,NK
                  CONC(K,J,I) = VALUES(K,L)
              ENDDO
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT, NK)
      DEALLOCATE(CONC)
      NC_WRITE_SZLJ_BEDTOP=STAT
    END FUNCTION

    INTEGER FUNCTION NC_WRITE_CONM(NC_ID, SNAPSHOT, IDX, VALUES, CLSCNT, NK, IK)
      ! ** WRITES DYE/TOX/SED/SND TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, CLSCNT, NK, IK
      REAL,INTENT(IN) :: VALUES(:,:,:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:,:) :: CONC
      INTEGER :: I,J,L,K,M,NS,STAT

      ALLOCATE(CONC(CLSCNT,NK,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              DO M=1,NK
                  IF( IK == 1 )THEN
                    K = NK - M + 1
                  ELSE
                    K = M
                  ENDIF
                  DO NS=1,CLSCNT
                      CONC(NS,M,J,I) = VALUES(L,K,NS)
                  ENDDO
              ENDDO
          ENDIF
      ENDDO
      STAT=NC_WRITE_CLASS3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT,CLSCNT,NK)
      DEALLOCATE(CONC)
      NC_WRITE_CONM=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_SZLJ_CONM(NC_ID, SNAPSHOT, IDX, VALUES, CLSCNT, NK)
      ! ** WRITES DYE/TOX/SED/SND TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, CLSCNT, NK
      REAL(RKD),INTENT(IN) :: VALUES(:,:,:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:,:) :: CONC
      INTEGER :: I,J,L,K,NS,STAT

      ALLOCATE(CONC(CLSCNT,NK,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
              DO K=1,NK
                  DO NS=1,CLSCNT
                      CONC(NS,K,J,I) = VALUES(NS,K,L)
                  ENDDO
              ENDDO
          ENDIF
      ENDDO
      STAT=NC_WRITE_CLASS3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT,CLSCNT,NK)
      DEALLOCATE(CONC)
      NC_WRITE_SZLJ_CONM=STAT
    END FUNCTION
      
    INTEGER FUNCTION NC_WRITE_LPT(NC_ID, SNAPSHOT) 
      ! ** WRITES LPT VARIABLES
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
      INTEGER :: NP,STAT
      REAL(8)   :: XM(1),YM(1),XLL(1),YLL(1)
  
      TRK_LON = MISSING_VALUE
      TRK_LAT = MISSING_VALUE
      TRK_LEV = MISSING_VALUE
      TRK_VOL = MISSING_VALUE
 
      DO NP=1,NPD
        IF( LLA_Global(NP) < 2 .OR. LLA_Global(NP) > LA_Global) CYCLE ! *** Modified for global
        XM(1)=XLA(NP)
        YM(1)=YLA(NP)
        CALL UTMR_WGS84(XM,YM,XLL,YLL)
        TRK_LON(NP) = XLL(1)
        TRK_LAT(NP) = YLL(1)
        TRK_LEV(NP) = ZLA(NP)
        IF( ANY(ISOILSPI == 1)) TRK_VOL(NP) = DVOL(NP)
      ENDDO
    
      STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_XC,SNAPSHOT,TRK_LON) 
      STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_YC,SNAPSHOT,TRK_LAT)
      STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_SIG,SNAPSHOT,TRK_LEV)
      IF( ANY(ISOILSPI == 1)) STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_VOL,SNAPSHOT,TRK_VOL)   
      NC_WRITE_LPT=STAT    
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_FIELD(NC_ID, SNAPSHOT, IDX, FIELD)
      ! ** WRITES FIELD TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT, IDX
      REAL,INTENT(IN) :: FIELD(:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: VALUES
      INTEGER :: I,J,L,STAT

      ALLOCATE(VALUES(ROWCNT,COLCNT))
      VALUES = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
            VALUES(J,I) = FIELD(L)  
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(IDX), VALUES, SNAPSHOT)
      DEALLOCATE(VALUES)
      NC_WRITE_FIELD=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_ORG_BEDTOP(NC_ID, SNAPSHOT, IDX, FIELD)
      ! ** WRITES FIELD TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT, IDX
      INTEGER,INTENT(IN) :: FIELD(:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: VALUES
      INTEGER :: I,J,L,STAT

      ALLOCATE(VALUES(ROWCNT,COLCNT))
      VALUES = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN
            VALUES(J,I) = KB - FIELD(L) + 1
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(IDX), VALUES, SNAPSHOT)
      DEALLOCATE(VALUES)
      NC_WRITE_ORG_BEDTOP=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_SZLJ_FIELD(NC_ID, SNAPSHOT, IDX, FIELD)
      ! ** WRITES FIELD TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT, IDX
      REAL(RKD),INTENT(IN) :: FIELD(:)
      REAL(4), ALLOCATABLE, DIMENSION(:,:) :: VALUES
      INTEGER :: I,J,L,STAT

      ALLOCATE(VALUES(ROWCNT,COLCNT))
      VALUES = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN  
            VALUES(J,I) = FIELD(L)  
          ENDIF
      ENDDO
      STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(IDX), VALUES, SNAPSHOT)
      DEALLOCATE(VALUES)
      NC_WRITE_SZLJ_FIELD=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_SZLJ_CLS(NC_ID, SNAPSHOT, IDX, VALUES, CLSCNT, FACTOR)
      ! ** WRITES BED TOX/SED/SND TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, CLSCNT
      REAL,INTENT(IN) :: VALUES(:,:),FACTOR
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: CONC
      INTEGER :: I,J,L,NS,STAT,STR2D(4),CNT2D(4)

      ALLOCATE(CONC(CLSCNT,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN
            DO NS=1,CLSCNT
              CONC(NS,J,I) = VALUES(L,NS)*FACTOR
            ENDDO
          ENDIF
      ENDDO
      
      STR2D = (/1, 1, 1, SNAPSHOT/)
      CNT2D = (/CLSCNT,ROWCNT,COLCNT,1/)
      STAT=NF90_PUT_VAR(NC_ID, NC_IDX(IDX), CONC, STR2D, CNT2D)
      
      DEALLOCATE(CONC)
      NC_WRITE_SZLJ_CLS=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_SZLJ_CLS_DBL(NC_ID, SNAPSHOT, IDX, VALUES, CLSCNT, FACTOR)
      ! ** WRITES BED TOX/SED/SND TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, CLSCNT
      REAL(RKD),INTENT(IN) :: VALUES(:,:),FACTOR
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: CONC
      INTEGER :: I,J,L,NS,STAT,STR2D(4),CNT2D(4)

      ALLOCATE(CONC(CLSCNT,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN
            DO NS=1,CLSCNT
              CONC(NS,J,I) = VALUES(L,NS)*FACTOR
            ENDDO
          ENDIF
      ENDDO
      
      STR2D = (/1, 1, 1, SNAPSHOT/)
      CNT2D = (/CLSCNT,ROWCNT,COLCNT,1/)
      STAT=NF90_PUT_VAR(NC_ID, NC_IDX(IDX), CONC, STR2D, CNT2D)
      
      DEALLOCATE(CONC)
      NC_WRITE_SZLJ_CLS_DBL=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_SZLJ_TOX(NC_ID, SNAPSHOT)
      ! ** WRITES BED TOX/SED/SND TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:,:) :: CONC
      REAL :: TMP
      INTEGER :: I,J,L,K,NT,STAT

      ALLOCATE(CONC(NTOX,KB,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO NT=1,NTOX
        DO L = 2,LA_Global
          I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
          J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
          IF( HP_Global(L) > HDRY )THEN
            DO K=1,KB
              IF( K < 3 .AND. TSED_Global(K,L) > 0. )THEN
                TMP = 0.01*TSED_Global(K,L)/BULKDENS_Global(K,L)             ! *** HBED(L,K)
                TMP = TMP/HBED_Global(L,KBT_Global(L))
                TMP = TMP*TOXB_Global(L,KBT_Global(L),NT)
                CONC(NT,K,J,I) = TMP
              ELSEIF( TSED(K,L) > 0. )THEN
                CONC(NT,K,J,I) = TOXB_Global(L,K,NT)
              ELSE
                CONC(NT,K,J,I) = 0.0_4
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      STAT=NC_WRITE_CLASS3D(NC_ID,NC_IDX(56),CONC,SNAPSHOT, NTOX, KB)
      DEALLOCATE(CONC)
      NC_WRITE_SZLJ_TOX=STAT
    END FUNCTION
    
    INTEGER FUNCTION NC_WRITE_BED_MASS(NC_ID, SNAPSHOT)
      ! ** WRITES BED TOX/SED/SND TO A NETCDF FILE
      INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT
      REAL(4), ALLOCATABLE, DIMENSION(:,:,:,:) :: CONC
      INTEGER :: I,J,L,K,M,NS,NSXD,STAT
      
      NSXD = 0
      IF( ISTRAN(6) >= 1 )THEN
        NSXD = NSXD + NSED
      ENDIF
      IF( ISTRAN(7) >= 1 )THEN
        NSXD = NSXD + NSND
      ENDIF

      ALLOCATE(CONC(NSXD,KB,ROWCNT,COLCNT))
      CONC = MISSING_VALUE
      DO L = 2,LA_Global
        I = sorted_loc_to_glob(L).global_i - IMN_Global + 1
        J = sorted_loc_to_glob(L).global_j - JMN_Global + 1
        IF( HP_Global(L) > HDRY )THEN
          DO K=1,KB
            M = KB - K + 1
            IF( ISTRAN(6) >= 1 )THEN
              DO NS=1,NSED
                CONC(NS,M,J,I) = SEDB_Global(L,K,NS)
              ENDDO
            ENDIF
            IF( ISTRAN(7) >= 1 )THEN
              DO NS=1,NSND
                CONC(NS+NSED,M,J,I) = SNDB_Global(L,K,NS)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      
      STAT=NC_WRITE_CLASS3D(NC_ID,NC_IDX(46),CONC,SNAPSHOT, NSXD, KB)
      DEALLOCATE(CONC)
      NC_WRITE_BED_MASS=STAT
    END FUNCTION
    
    SUBROUTINE NC_NEW_FILE(NC_ID,FILENAME)
      INTEGER(4) :: NC_ID
      CHARACTER(*),INTENT(IN) :: FILENAME
      INTEGER(4) :: STATUS

      STATUS = NC_CREATE_FILE(NC_ID,FILENAME) 
      STATUS = NC_WRITE_GRID(NC_ID,XCR(1:4,JMN_Global:JMX_Global,IMN_Global:IMX_Global),YCR(1:4,JMN_Global:JMX_Global,IMN_Global:IMX_Global))
      STATUS = NC_WRITE_ZBOT(NC_ID)     
      NTI = 0
    END SUBROUTINE 

    SUBROUTINE NC_WRITE_SNAPSHOT(NC_ID, TIMESEC)

      DOUBLE PRECISION,INTENT(IN) :: TIMESEC
      INTEGER(4),INTENT(IN) :: NC_ID
      INTEGER :: STATUS, M
      REAL(8) :: TI(1)
  
      NTI = NTI + 1
      TI = TIMEDAY
      WRITE(*,'(A22,F15.3)') 'NETCDF OUTPUT @ HOURS: ',TI(1)*24
      STATUS = NF90_PUT_VAR(NC_ID, NC_TIM, TI, (/NTI/), (/1/))
      STATUS = NC_WRITE_WSEL(NC_ID, NTI)
      STATUS = NC_WRITE_VEL(NC_ID, NTI)
      IF( ISNCDF(10) == 1) STATUS = NC_WRITE_SHEAR(NC_ID, NTI)
      IF( ISNCDF(11) == 1) STATUS = NC_WRITE_WIND(NC_ID, NTI)
      IF( ISNCDF(12) == 1) STATUS = NC_WRITE_WAVE(NC_ID, NTI)
	   
      IF( ISNCDF(1) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 12, SAL_Global, KC, 1)
      IF( ISNCDF(2) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 13, TEM_Global, KC, 1)
      IF( ISNCDF(3) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 14, DYE_Global, NDYE, KC, 1)
      IF( ISNCDF(4) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 15, SFL_Global, KC, 1)	  
      IF( ISNCDF(5) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 16, TOX_Global, NTOX, KC, 1)
      IF( ISNCDF(6) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 17, SED_Global, NSED, KC, 1)
      IF( ISNCDF(7) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 18, SND_Global, NSND, KC, 1)
    
      IF( ISNCDF(8) == 1 )THEN
        DO M=1,NWQV
          IF(ISTRWQ(M) > 0) STATUS = NC_WRITE_CONS(NC_ID, NTI, 18 + M, WQV(:,:,M), KC, 1)
        ENDDO
      ENDIF
  
      ! ** LPT
      IF(ISNCDF(9) == 1 )THEN     
        IF( ANY(ISOILSPI == 1)) STATUS = NC_WRITE_FIELD(NC_ID, NTI, 41, MOC)  
        STATUS = NC_WRITE_LPT(NC_ID, NTI) 
      ENDIF
      
      IF( LSEDZLJ )THEN
        ! *** SEDZLJ SEDIMENT MODEL
        STATUS = NC_WRITE_SZLJ_BEDTOP(NC_ID, NTI, 42, LAYERACTIVE, KB)
        STATUS = NC_WRITE_SZLJ_FIELD(NC_ID, NTI, 47, TAU_Global)
        STATUS = NC_WRITE_SZLJ_CONS(NC_ID, NTI, 46, TSED_Global, KB)
        STATUS = NC_WRITE_SZLJ_CONS(NC_ID, NTI, 48, BULKDENS_Global, KB)
        STATUS = NC_WRITE_SZLJ_CONM(NC_ID, NTI, 49, PERSED_Global, NSCM, KB)
        STATUS = NC_WRITE_SZLJ_FIELD(NC_ID, NTI, 50, D50AVG_Global)
        STATUS = NC_WRITE_SZLJ_FIELD(NC_ID, NTI, 51, ETOTO_Global)
        STATUS = NC_WRITE_SZLJ_FIELD(NC_ID, NTI, 52, DEPO_Global)
        IF( NCALC_BL > 0 )THEN
          STATUS = NC_WRITE_SZLJ_CLS_DBL(NC_ID, NTI, 53, CBL_Global, NSCM, 10000.0_8)
          STATUS = NC_WRITE_SZLJ_CLS(NC_ID, NTI, 54, QSBDLDX_Global, NSCM, 1.)
          STATUS = NC_WRITE_SZLJ_CLS(NC_ID, NTI, 55, QSBDLDY_Global, NSCM, 1.)
        ENDIF
        
        IF(ISNCDF(5) == 1 .AND. ISTRAN(5) > 0 )THEN
          STATUS = NC_WRITE_SZLJ_TOX(NC_ID, NTI)
          IF( NCALC_BL > 0 )THEN
            STATUS = NC_WRITE_SZLJ_CLS_DBL(NC_ID, NTI, 57, CBLTOX_Global, NTOX, 1.0_8)
          ENDIF
        ENDIF
      ELSEIF( ISBEXP >= 1 .AND. KB > 1 .AND. ( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) )THEN
        ! *** ORIGINAL SEDIMENT MODEL
          STATUS = NC_WRITE_ORG_BEDTOP(NC_ID, NTI, 42, KBT_Global)
          STATUS = NC_WRITE_CONS(NC_ID, NTI, 43, HBED_Global, KB, 1)
          STATUS = NC_WRITE_CONS(NC_ID, NTI, 44, BDENBED_Global, KB, 1)
          STATUS = NC_WRITE_CONS(NC_ID, NTI, 45, PORBED_Global, KB, 1)
          IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
            STATUS = NC_WRITE_BED_MASS(NC_ID, NTI)
          ENDIF
          IF(ISNCDF(5) == 1 .AND. ISTRAN(5) >= 1 )THEN
            STATUS = NC_WRITE_CONM(NC_ID, NTI, 56, TOXB_Global, NTOX, KB, 1)
          ENDIF
      ENDIF
            
    END SUBROUTINE 
  
    SUBROUTINE NETCDF_WRITE(NC_ID)
  
      INTEGER(4),INTENT(IN) :: NC_ID

      REAL(8) :: TIMEHOUR
    
      CHARACTER(10) :: SDATE,DATESTAMP
      CHARACTER(80) :: FILENAME
      INTEGER, SAVE :: IFRST = 1
    
      ! *** Call subroutine to gather and map variables in prep for writing to the NetCDF file(s)
      !Call Map_Write_NetCDF
 
      if(process_id == master_id )then
          ! *** Single netcdf file
          IF( ISSGLFIL == 1 )THEN
            ! *** only create the file at the first call TOGREGOR this subroutine
            IF( IFRST == 1 )THEN
              FILENAME = './/#output//DSI_'//TRIM(PROJ)//'.nc'
              CALL NC_NEW_FILE(NC_ID, FILENAME)
              IFRST = 0
            ENDIF    
        
          ! *** multiple netcdf files 
          ELSEIF (ISSGLFIL == 0 )THEN
            TIMEHOUR = TIMEDAY*24
            CALL TOGREGOR(SDATE,TIMEDAY)
            DATESTAMP = SDATE(1:4)//'_'//SDATE(5:6)//'_'//SDATE(7:8)
            IF((NC_DATESTAMP /= DATESTAMP) )THEN
              IF( IFRST == 0 )THEN
                CALL NC_WRITE_SNAPSHOT(NC_ID, TIMESEC)
                CALL NC_CLOSE_FILE(NC_ID)  
                IF( TIMEDAY >= TENDNCDF) RETURN
              ENDIF
              FILENAME = './/#output//DSI_'//TRIM(PROJ)//'_'//DATESTAMP//'.nc'
              CALL NC_NEW_FILE(NC_ID, FILENAME)
              IFRST = 0
            ENDIF
            NC_DATESTAMP = DATESTAMP
          ENDIF 
      
          CALL NC_WRITE_SNAPSHOT(NC_ID, TIMESEC)
          
      end if
          
      END SUBROUTINE 
#endif

END MODULE
