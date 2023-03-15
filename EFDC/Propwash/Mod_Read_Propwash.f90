! ----------------------------------------------------------------------
!   This file is a part of EFDC+
!   Website:  https://eemodelingsystem.com/
!   Repository: https://github.com/dsi-llc/EFDC_Plus.git
! ----------------------------------------------------------------------
! Copyright 2021-2022 DSI, LLC
! Distributed under the GNU GPLv2 License.
! ----------------------------------------------------------------------
!---------------------------------------------------------------------------!
!                     EFDC+ Developed by DSI, LLC.
!---------------------------------------------------------------------------!
!< @details Contains subroutines to read propwash files
!< @author  Zander Mausolff & Paul Craig
!---------------------------------------------------------------------------!
Module Mod_Read_Propwash

  Use GLOBAL, only : timeday
  Use Variables_Propwash
  Use Variables_Ship

  Use Mod_Active_Ship
  Use Mod_All_Tracks
  USE Variables_MPI
  
  implicit none

  Integer :: total_prop_ships
  Integer :: total_track_ships

  Contains

  !---------------------------------------------------------------------------!
  !< @details Broadcasts the propwash parameters to the other domains
  !---------------------------------------------------------------------------!
  Subroutine Broadcast_Propwash
    USE Variables_MPI
    Use Broadcast_Routines
  
    Integer :: i, j, k
    
    ! *** Universal parameters
    Call Broadcast_Scalar(total_ships,         master_id)
    Call Broadcast_Scalar(num_radial_elems,    master_id)
    Call Broadcast_Scalar(num_axial_elems,     master_id)
    Call Broadcast_Scalar(mesh_width,          master_id)
    Call Broadcast_Scalar(mesh_Length,         master_id)
    Call Broadcast_Scalar(efflux_zone_mult,    master_id)
    Call Broadcast_Scalar(flow_est_zone1_mult, master_id)
    Call Broadcast_Scalar(flow_est_zone2_mult, master_id)
    Call Broadcast_Scalar(efflux_mag_mult,     master_id)

    if( process_id /= master_id .and. total_ships > 0 )then
      ! *** Ships
      Allocate(all_new_ships(total_ships))
      
      ! *** Tracks
      allocate(all_read_tracks(total_ships))
    endif

    ! *** Ship and Ship Tracks

    do i = 1,total_ships
      Call Broadcast_Scalar(all_new_ships(i).mmsi               , master_id)
      Call Broadcast_Scalar(all_new_ships(i).name               , master_id)

      Call Broadcast_Scalar(all_new_ships(i).length             , master_id)
      Call Broadcast_Scalar(all_new_ships(i).beam               , master_id)
      Call Broadcast_Scalar(all_new_ships(i).max_draft          , master_id)
      Call Broadcast_Scalar(all_new_ships(i).min_draft          , master_id)
      Call Broadcast_Scalar(all_new_ships(i).max_power          , master_id)
      Call Broadcast_Scalar(all_new_ships(i).max_rps            , master_id)
      Call Broadcast_Scalar(all_new_ships(i).ais_to_stern       , master_id)
      Call Broadcast_Scalar(all_new_ships(i).dist_between_props , master_id)
      Call Broadcast_Scalar(all_new_ships(i).dist_from_stern    , master_id)
      Call Broadcast_Scalar(all_new_ships(i).prop_offset        , master_id)

      Call Broadcast_Scalar(all_new_ships(i).num_props          , master_id)
      Call Broadcast_Scalar(all_new_ships(i).num_blades         , master_id)
      Call Broadcast_Scalar(all_new_ships(i).ducted             , master_id)
      Call Broadcast_Scalar(all_new_ships(i).thrust_coeff       , master_id)
      Call Broadcast_Scalar(all_new_ships(i).prop_diam          , master_id)
      Call Broadcast_Scalar(all_new_ships(i).prop_hub_diam      , master_id)
      Call Broadcast_Scalar(all_new_ships(i).blade_area_ratio   , master_id)
      Call Broadcast_Scalar(all_new_ships(i).pitch_ratio        , master_id)

      Call Broadcast_Scalar(all_new_ships(i).freq_out           , master_id)
      
      Call Broadcast_Scalar(all_read_tracks(i).num_tracks       , master_id)
      
      if( process_id /= master_id .and. all_read_tracks(i).num_tracks  > 0 )then
        ! *** Setup array containing all ships to be simulated
        allocate(all_read_tracks(i).all_ship_tracks(all_read_tracks(i).num_tracks))
      endif
    
      do j = 1,all_read_tracks(i).num_tracks
        Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).num_positions, master_id)

        if( process_id /= master_id .and. all_read_tracks(i).all_ship_tracks(j).num_positions  > 0 )then
          ! *** Track positions
          allocate(all_read_tracks(i).all_ship_tracks(j).track_pos(all_read_tracks(i).all_ship_tracks(j).num_positions))
        endif
        
        do k = 1,all_read_tracks(i).all_ship_tracks(j).num_positions
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).time,     master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).x_pos,    master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).y_pos,    master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).z_pos,    master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).speed,    master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).course,   master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).heading,  master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).draft,    master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).rps,      master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).power,    master_id)
          Call Broadcast_Scalar(all_read_tracks(i).all_ship_tracks(j).track_pos(k).freq_out, master_id)
        enddo
      enddo
    enddo

  End Subroutine Broadcast_Propwash
  
  !---------------------------------------------------------------------------!
  !< @details Reads the JSON formatted propwash parameters
  !< @param[in] debug (optional) argument
  !---------------------------------------------------------------------------!
  subroutine Read_propwash_config(debug)

    use Variables_Propwash
    use fson
    use mod_fson_value, only: fson_value_count, fson_value_get

    implicit none

    ! *** Dummy variables
    logical, optional :: debug

    ! *** Local variables
    type(fson_value), pointer :: json_data, prop_item
    integer :: i, ns, nt, L

    write(*,'(A)') 'READING Propwash Parameters'

    ! *** Open the propwash ship data file
    json_data => fson_parse("propwash_config.jnp")

    ! *** mesh portion
    call fson_get(json_data, "mesh.num_axial_elems",      num_axial_elems)
    call fson_get(json_data, "mesh.num_radial_elems",     num_radial_elems)
    call fson_get(json_data, "mesh.mesh_length",          mesh_length)
    call fson_get(json_data, "mesh.mesh_width",           mesh_width)

    ! *** parameters that define propwash velocity field
    call fson_get(json_data, "parms.efflux_zone_mult",    efflux_zone_mult)
    call fson_get(json_data, "parms.flow_est_zone1_mult", flow_est_zone1_mult)
    call fson_get(json_data, "parms.flow_est_zone2_mult", flow_est_zone2_mult)

    ! *** Efflux velocity multiplier
    call fson_get(json_data, "parms.efflux_mag_mult",     efflux_mag_mult)

    ! *** Fast Settling Classes
    call fson_get(json_data, "parms.fraction_fast",       fraction_fast)
    call fson_get(json_data, "parms.fast_multiplier",     fast_multiplier)

    ! *** write out what was just read in for debugging purposes
    !if(present(debug) )then
    if( debug )then
      OPEN(prop_log_unit,FILE=OUTDIR//'propeller_mesh.log',STATUS='UNKNOWN')
      write(prop_log_unit, '(a, f10.2)')' mesh_length         : ', mesh_length
      write(prop_log_unit, '(a, f10.2)')' mesh_width          : ', mesh_width
      write(prop_log_unit, '(a, i7)' )  ' num_axial_elems     : ', num_axial_elems
      write(prop_log_unit, '(a, i7)' )  ' num_radial_elems    : ', num_radial_elems
      write(prop_log_unit, '(a, f10.2)')' efflux_zone_mult    : ', efflux_zone_mult
      write(prop_log_unit, '(a, f10.2)')' flow_est_zone1_mult : ', flow_est_zone1_mult
      write(prop_log_unit, '(a, f10.2)')' flow_est_zone2_mult : ', flow_est_zone2_mult
      write(prop_log_unit, '(a)' ) '-----------------------------------'
      write(prop_log_unit, '(a, f10.3)')' efflux_mag_mult     : ', efflux_mag_mult
      write(prop_log_unit, '(a)' ) '-----------------------------------'
      write(prop_log_unit, '(a, f10.3)')' fraction_fast       : ', fraction_fast
      write(prop_log_unit, '(a, f10.3)')' fast_multiplier     : ', fast_multiplier
      write(prop_log_unit, '(a)' ) '-----------------------------------'
      close(prop_log_unit)
    end if
    
    ! *** Build the lookup tables for fast class mass erosion
    IFASTCLASS = 0
    NSCM2 = NSCM
    IF( fraction_fast > 0 .and. fast_multiplier > 0. )THEN
      DO NS = 1,NSCM
        IF( WSEDO(NS) > 0.0 .AND. SEDDIA(ns) < 65./1e6 )THEN   ! *** "Fast" settling classes due to mass erosion of cohesives.  Ignore washload.
          NSCM2 = NSCM2 + 1
          IWC2BED(NSCM2) = NS                         ! *** Map WC class to Bed class
          IBED2WC(NS) = NSCM2                         ! *** Map bed class to WC class
          WSEDO(NSCM2) = WSEDO(NS)*fast_multiplier
          IF( LSEDZLJ ) DWS(NSCM2) = DWS(NS)*fast_multiplier
        ELSE
          IWC2BED(NS) = NS                            ! *** Map WC class to Bed class
          IBED2WC(NS) = NS                            ! *** Map bed class to WC class
        ENDIF
      ENDDO
      
      ! *** Trim the active constituent list from unneeded fast classes
      NACTIVEWC = NACTIVEWC - (NSED2 - NSCM2 - NSND)
      
      ! *** Chemical Fate and Transport Parameters
      IF( ISTRAN(5) > 0 .AND. NTOX > 0 )THEN
        ! *** Set the new fast class CFT variables based on the source sediment class
        ! *** ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT), STFPOCB, STFPOCW
        DO NT = 1,NTOX
          DO NS = NSCM+1,NSCM2
            ITXPARW(NS,NT) = ITXPARW(IWC2BED(NS),NT)
            CONPARW(NS,NT) = CONPARW(IWC2BED(NS),NT)
            TOXPARW(2:LA,NS,NT) = TOXPARW(2:LA,IWC2BED(NS),NT)
            TOXPARB(2:LA,NS,NT) = TOXPARB(2:LA,IWC2BED(NS),NT)
            
            ! *** SPATIALLY VARIABLE PARAMETERS
            DO L = 2,LA
              STFPOCW(L,:,NS) = STFPOCW(L,:,IWC2BED(NS))
              STFPOCB(L,:,NS) = STFPOCB(L,:,IWC2BED(NS))
            ENDDO
          ENDDO
          
          ! *** Update each solids class for each toxic: NSP2
          NSP2(NT) = NSCM2                              ! *** Kd  Approach ISTOC(NT)=0
          IF( ISTOC(NT) == 1 ) NSP2(NT) = NSP2(NT) + 2  ! *** DOC and POC (POC non-sediment related)  (3 Phase)
          IF( ISTOC(NT) == 2 ) NSP2(NT) = NSP2(NT) + 1  ! *** DOC and POC fractionally distributed    (3 Phase)
          ! *** ISTOC(NT)=0 and ISTOC(NT)=3                   POC fOC*SED/SND BASED ONLY              (2 Phase)
        ENDDO
        
      ENDIF
      
    ENDIF
    NSED2 = NSCM2 - NSND           ! *** Mass eroded classes are added at the end of the NSED + NSND (original) or NSCM (SEDZLJ)
    
  end subroutine Read_propwash_config

  !---------------------------------------------------------------------------!
  !< @details Reads the JSON formatted ship data
  !< @param[in] debug (optional) argument
  !---------------------------------------------------------------------------!
  Subroutine Read_Ship_Data(debug)

    use Variables_Propwash
    Use Allocate_Initialize      
    use fson
    use mod_fson_value, only: fson_value_count, fson_value_get
    Use Variables_MPI, only : mpi_log_unit

    implicit none

    ! *** Dummy variables
    logical, optional :: debug

    ! *** Local variables
    type(fson_value), pointer :: json_ship_data, prop_item, all_ship_array
    type(ship_type)           :: new_ship

    integer :: i, j

    write(*,'(A)') 'READING Ship Data'

    ! *** Open the propwash ship data file
    json_ship_data => fson_parse("propwash_ships.jnp")

    ! *** get all of the ship data
    call fson_get(json_ship_data, "all_ships", all_ship_array)

    ! *** get the total number of ships in the ship data file
    total_prop_ships = fson_value_count(all_ship_array)

    ! *** Set the size of the array of types that will contain all of the ship data
    if( total_prop_ships > 0 )then
      Allocate(all_new_ships(total_prop_ships))
    else
      Call Stopp("***WARNING*** no ships are present in the propwash_ships.jnp file!")
    end if

    ! *** Loop over all ships in the input file
    do i = 1, total_prop_ships

      ! *** Instantiate ship object
      new_ship = ship_type()

      ! *** get a single ship from the array of ships
      prop_item => fson_value_get(all_ship_array, i)

      ! *** extract the values for that particular ship
      call fson_get(prop_item, "mmsi",               new_ship.mmsi)
      call fson_get(prop_item, "name",               new_ship.name )
      call fson_get(prop_item, "length",             new_ship.length )
      call fson_get(prop_item, "min_draft",          new_ship.min_draft)
      call fson_get(prop_item, "max_draft",          new_ship.max_draft)
      call fson_get(prop_item, "max_power",          new_ship.max_power)
      call fson_get(prop_item, "max_rps",            new_ship.max_rps)
      call fson_get(prop_item, "ais_to_stern",       new_ship.ais_to_stern)
      call fson_get(prop_item, "dist_from_stern",    new_ship.dist_from_stern)
      call fson_get(prop_item, "dist_between_props", new_ship.dist_between_props)
      call fson_get(prop_item, "prop_offset",        new_ship.prop_offset)
      call fson_get(prop_item, "freq_out",           new_ship.freq_out)

      call fson_get(prop_item, "num_props",          new_ship.num_props)
      call fson_get(prop_item, "num_blades",         new_ship.num_blades)
      call fson_get(prop_item, "ducted",             new_ship.ducted)
      call fson_get(prop_item, "thrust_coeff",       new_ship.thrust_coeff)
      call fson_get(prop_item, "prop_diam",          new_ship.prop_diam)
      call fson_get(prop_item, "prop_hub_diam",      new_ship.prop_hub_diam)
      call fson_get(prop_item, "blade_area_ratio",   new_ship.blade_area_ratio)
      call fson_get(prop_item, "pitch_ratio",        new_ship.pitch_ratio)

      !Call AllocateDSI(new_ship.mesh_count, 0:LCM, 0)
      Allocate(new_ship.mesh_count(0:LCM))
      new_ship.mesh_count = 0

      ! *** Special Cases
      call fson_get(prop_item, "num_fixed_cells",    new_ship.num_fixed_cells)
      if( new_ship.num_fixed_cells > 0 )then
        Call AllocateDSI(new_ship.fixed_cells, new_ship.num_fixed_cells, 0)
        Call AllocateDSI(new_ship.fixed_frac,  new_ship.num_fixed_cells, 0.0)
        call fson_get(prop_item, "fixed_cells",      new_ship.fixed_cells)
        call fson_get(prop_item, "fixed_frac",       new_ship.fixed_frac)
      endif
      
      ! *** QC missing data
      if( new_ship.prop_diam == -999. )then
        print '(a,2x,a)','Propeller diameter missing for ship:', new_ship.name
        Call Stopp('.')
      endif
      if( new_ship.blade_area_ratio == -999. )then
        print '(a,2x,a)','Propeller blade area ratio missing for ship:', new_ship.name
        Call Stopp('.')
      endif
      if( new_ship.pitch_ratio == -999. )then
        print '(a,2x,a)','Propeller blade pitch ratio missing for ship:', new_ship.name
        Call Stopp('.')
      endif
      if( new_ship.num_blades == -999. )then
        print '(a,2x,a)','Number of propeller blades missing for ship:', new_ship.name
        Call Stopp('.')
      endif

      if( new_ship.length == -999. )          new_ship.length          = 30.                 ! *** 30 meters for typical tug
      if( new_ship.ais_to_stern == -999. )    new_ship.ais_to_stern    = new_ship.length/2.
      if( new_ship.dist_from_stern == -999. ) new_ship.dist_from_stern = 0.                  ! *** Assume propeller at stern
      if( new_ship.prop_offset == -999. )     new_ship.prop_offset     = 0.                  ! *** Assume bottom of prop at reported draft
      if( new_ship.min_draft == -999. )       new_ship.min_draft       = 0.
      if( new_ship.max_draft == -999. )       new_ship.max_draft       = 100.
      if( new_ship.num_props == -999. )       new_ship.num_props       = 1
      if( new_ship.prop_hub_diam == -999. )   new_ship.prop_hub_diam   = 0.1*new_ship.prop_diam
      if( new_ship.max_rps == -999. )         new_ship.max_rps         = 1                   ! *** Set to minimum to prevent / by 0
      if( new_ship.ducted == -999. )then
        if( new_ship.max_power > 1000 )then
          new_ship.ducted = 1
        else
          new_ship.ducted = 0
        endif
      endif
      if( new_ship.num_props > 1 .and. new_ship.dist_between_props == -999. ) new_ship.dist_between_props = 2.*new_ship.prop_diam

      ! *** Add the new ship info to the global list of all the ship data
      all_new_ships(i) = new_ship

    end do

    ! *** Write out what we just read for debugging purposes
    if( debug )then
      do i = 1, total_prop_ships
        call all_new_ships(i).write_out(mpi_log_unit)
      end do
    end if

    return

  End Subroutine Read_ship_Data

  !---------------------------------------------------------------------------!
  !< @details Reads the JSON formatted ship track data
  !< @param[in] debug (optional) argument
  !---------------------------------------------------------------------------!
  Subroutine Read_Ship_Tracks(debug)

    use Variables_Propwash
    use fson
    use mod_fson_value, only: fson_value_count, fson_value_get
    Use Variables_MPI, only : mpi_log_unit

    implicit none

    ! *** Dummy variables
    logical, optional :: debug                 !< will print some read in parameters to log file

    ! *** Local variables
    integer :: i, j, k, nValid
    integer :: total_track_ships               !< total number of ships to be simulated
    integer :: total_tracks                    !< total number of tracks for a respective ship
    integer :: total_pts_in_track              !< total number of points in a given track for a respective ship
    type(fson_value), pointer :: json_data, ship_array, track_array,track_pos_array, item_ship, item_track, item_track_pos
    type(position_cell), Allocatable, Dimension(:) :: read_track_position !< Array of all track positions for the ship
    type(position_cell) :: new_position
    type(all_tracks), allocatable, dimension(:) :: all_track_pos

    integer :: track_ship_id
    integer :: test, pointcount

    ! *** Read the input data file

    write(*,'(A)') 'READING Propwash Track Data'
    
    !write(*,'(a)') '*** Getting JSON data'  ! delme
    json_data => fson_parse('propwash_tracks.jnp')
    !write(*,'(a)') '*** Succesfully retrieved JSON data'  ! delme

    ! *** Get array of all ships tracks
    call fson_get(json_data, 'all_ship_tracks', ship_array)

    ! *** Get total number of ships
    total_track_ships = fson_value_count(ship_array)

    ! *** Check to make sure there is consistency between ship and track input files
    if( total_track_ships > total_prop_ships )then
      Call Stopp('***WARNING*** there is a mismatch in the number of ships to be simulated in the propwash input files')
    endif

    ! *** Set total_ships number that can be used by other subroutines
    total_ships = total_prop_ships

    ! *** Setup array containing all ships to be simulated
    allocate(all_read_tracks(total_ships))

    ! *** Loop over all of the ships
    do i = 1, total_track_ships

      ! *** Select a single ship
      item_ship => fson_value_get(ship_array, i)

      call fson_get(item_ship, 'mmsi', track_ship_id)

      ! *** Get single ships array of tracks
      call fson_get(item_ship, 'ship_tracks', track_array)

      ! *** Get total number tracks for that ship, each ship can have a differet number of tracks
      total_tracks = fson_value_count(track_array)

      Allocate(all_track_pos(total_tracks))
      ! *** Set the number of tracks per ship

      !write(900,'(a,i5,i11,i5)') '*** New Ship I, MMSI, Number of Tracks', i, track_ship_id, total_tracks  ! delme
      !write(*,'(a,i5,i11,i5)') '*** New Ship I, MMSI, Number of Tracks', i, track_ship_id, total_tracks  ! delme
      
      ! *** Loop over the number of tracks for this ship
      do j = 1, total_tracks

        ! *** Get the indvidual track
        item_track => fson_value_get(track_array, j)

        call fson_get(item_track, 'track', track_pos_array)

        ! *** Determine the number of points in a given track
        total_pts_in_track = fson_value_count(track_pos_array)
        
        !write(900,'(a,3i5)') '*** Track ', i, j, total_pts_in_track  ! delme
        !write(*,'(a,3i5)') '*** Track ', i, j, total_pts_in_track  ! delme

        ! *** Setup array that will hold all of the track positions for that given ship
        allocate(all_track_pos(j).track_pos(total_pts_in_track))

        ! *** save the number of points in this current track
        all_track_pos(j).num_positions = total_pts_in_track

        ! *** Loop over the number of positions in a track
        nValid = 0
        do k = 1, total_pts_in_track

          ! ** get single position
          item_track_pos => fson_value_get(track_pos_array, k)

          ! *** instantiate new track position
          new_position = position_cell()

          ! *** read the properties for that single track
          call fson_get(item_track_pos, 't', new_position.time)
          call fson_get(item_track_pos, 'x', new_position.x_pos)
          call fson_get(item_track_pos, 'y', new_position.y_pos)

          call fson_get(item_track_pos, 's', new_position.speed)
          call fson_get(item_track_pos, 'c', new_position.course)
          call fson_get(item_track_pos, 'h', new_position.heading)
          call fson_get(item_track_pos, 'd', new_position.draft)
          call fson_get(item_track_pos, 'p', new_position.power)
          call fson_get(item_track_pos, 'o', new_position.freq_out)

          !write(900,'(3i5,f12.5,2f10.1,8f10.3)') i, j, k,                         &
          !           new_position.time, new_position.x_pos, new_position.y_pos, new_position.speed, new_position.course, &  ! delme
          !           new_position.heading, new_position.draft, new_position.power, new_position.freq_out
          
          ! *** QC track points.  Other QC will be done when creating active_ship object
          if( new_position.time == -999. .or. new_position.x_pos == -999. .or. new_position.y_pos == -999 ) CYCLE
          if( new_position.power /= -999. .and. abs(new_position.power) > 1.0 )then
            if( new_position.power < -100. .or. new_position.power > 1. )then    ! delme - testing
              print *, 'Invalid applied power fraction for ship, Track, Pt:  ' , i, j, k, all_new_ships(i).name 
              Call Stopp('.')
            endif
          endif
          if( all_new_ships(i).max_power == -999. .and. new_position.power > 0.0 )then
            print *, 'Invalid applied HP for ship, Track, Pt:  ' , i, j, k, all_new_ships(i).name 
            Call Stopp('.')
          endif

          ! *** Append the array of all track positions
          nValid = nValid + 1
          all_track_pos(j).track_pos(nValid) = new_position

        end do ! over k positions
        all_track_pos(j).num_positions = nValid

      end do ! over j tracks

      all_read_tracks(i).all_ship_tracks = all_track_pos
      all_read_tracks(i).num_tracks      = total_tracks     ! *** total tracks for this ship

      ! *** deallocate track position array so it can change for next ship
      deallocate(all_track_pos)

    end do ! over i boats

    ! *** Write out track info if we are trying to debug
    !if(present(debug) )then
    if( debug )then
      do i = 1, total_track_ships
        write(mpi_log_unit, '(a, i9)') 'Ship #: ', i

        total_tracks = size(all_read_tracks(i).all_ship_tracks)

        do j = 1, total_tracks

          write(mpi_log_unit, '(a, i9)') 'Track #: ', j
          total_pts_in_track = size(all_read_tracks(i).all_ship_tracks(j).track_pos)

          do k = 1, total_pts_in_track
            ! *** Write out each track
            Call all_read_tracks(i).all_ship_tracks(j).track_pos(k).write_out(mpi_log_unit, j)
          end do

        end do
      end do
    end if

  End Subroutine Read_Ship_Tracks

End module Mod_Read_Propwash
!******************************************************************************************
!******************************************************************************************
