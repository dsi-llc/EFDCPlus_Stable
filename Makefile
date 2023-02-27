# Select the fortran compiler to use
FC=/opt/intel/oneapi/mpi/latest/bin/mpiifort


# Compiler options
#FFLAGS=-g -traceback -fpp -DENABLE_MPI -D_MPI -nofixed -qopenmp -DNCOUT
FFLAGS=-g -traceback -qopenmp -O3 -fpp -DENABLE_MPI -D_MPI -DLINUX -m64 -real-size 32 -fpe0 -fp-speculation=safe -multiple-processes=12
# Compilation options for debug mode
#FFLAGS=-g -traceback -qopenmp -debug all -check all -fpp -DENABLE_MPI -D_MPI -DDEBUGGING -fpic -real-size 64 -fpe0 -fp-speculation=safe -multiple-processes=8

# Setup the OpenMP compiler options
OMPFLAGS=#-qopenmp
FPPFLAGS= 

# Preprocessor

# Libs
LIBS=

PROGRAM = efdc
PROGEXT = .x

PROJHOME= ..

# Suffix rules
.SUFFIXES: .f90 .F90 .for .FOR .o

INCLUD   =   -L/usr/local/include -LWASP 
CDFDIR   =   -L/usr/local/lib
CDFLIBS  =   #-Bstatic -lnetcdf -lnetcdff -lhdf5_hl -lhdf5

# Include dependency list created by makedepf90
#include .depend 


# list directories containing other files
SRCS=*.f90 */*.f90 
	  
OBJS= EFDC/Propwash/Setup_Ships.o  EFDC/aaefdc.o EFDC/csedress.o EFDC/fsedmode.o EFDC/fdstrse.o EFDC/MPI_Out/FormatOutputMPI.o EFDC/csedvis.o EFDC/csndset.o \
      EFDC/fsbdld.o EFDC/csndeqc.o EFDC/fstrse.o EFDC/svbksb.o EFDC/setshld.o EFDC/csedresb.o EFDC/csedtaus.o EFDC/fprobdep.o \
      EFDC/setstvel.o EFDC/csndzeq.o EFDC/csedset.o EFDC/ludcmp.o EFDC/fhydcn.o EFDC/csedtaub.o EFDC/lubksb.o \
      EFDC/budget.o EFDC/welcome.o \
	  EFDC/subchan.o EFDC/rsalpltv.o EFDC/calvegser.o EFDC/caltox_kinetics.o \
	  EFDC/caldiff.o EFDC/bankerosed.o EFDC/Propwash/Propwash_Calc_Sequence.o EFDC/MPI_Out/Write_Array_Sizes.o EFDC/bal2t.o EFDC/MPI_Utilities/Report_Max_Min_Timing.o EFDC/MPI_Communication/Communicate_Groups.o \
	  EFDC/s_slope.o  EFDC/s_bedload.o EFDC/caltox.o \
	  EFDC/calsnd.o EFDC/caldisp2.o EFDC/MPI_Mapping/Setup_MapBackToGlobal.o EFDC/MPI_Domain_Decomp/Scan_Cell.o \
	  EFDC/s_sedic.o EFDC/rcahq.o EFDC/cbalod.o EFDC/calpuv2c.o EFDC/calbuoy.o EFDC/MPI_Domain_Decomp/Scan_JSON_Decomp.o EFDC/MPI_Mapping/Calculate_Local_Shear.o \
	  EFDC/wavebl.o EFDC/showval.o EFDC/out3d.o \
	  EFDC/calsed.o EFDC/calqq1.o EFDC/ainit.o \
	  EFDC/wasp4.o EFDC/varinit.o EFDC/rsurfplt.o EFDC/cellmap.o \
	  EFDC/calfqc.o EFDC/Propwash/Det_Adjacent_Cells.o EFDC/Eutrophication/wwqnc.o EFDC/calbed9.o EFDC/MPI_Utilities/Setup_MPI_Topology.o \
	  EFDC/MPI_Mapping/Map_Connectors.o \
	  EFDC/lsqharm.o EFDC/caltran_ad.o EFDC/caltoxb.o \
	  EFDC/calpser.o EFDC/calbed.o EFDC/wasp6.o \
	  EFDC/ceqicm.o EFDC/calhdmf3.o EFDC/calblay.o EFDC/MPI_Utilities/Initialize_MPI.o \
	  EFDC/MPI_Domain_Decomp/Child_Grid.o EFDC/MPI_Communication/Communicate_3D_0.o \
	  EFDC/varalloc.o EFDC/s_sedzlj.o EFDC/negdep.o EFDC/calbal.o \
	  EFDC/MPI_Mapping/Create_List_No_Ghost_Cells.o EFDC/MPI_Domain_Decomp/Read_JSON_Decomp.o \
	  EFDC/ssedtox.o EFDC/s_tecplot.o EFDC/setbcs.o EFDC/mhkpwr.o EFDC/caluvw.o \
	  EFDC/calqvs.o EFDC/calpnhs.o EFDC/MPI_Domain_Decomp/Congrad_MPI.o \
	  EFDC/tmsr.o EFDC/timelog.o EFDC/setfpocb.o EFDC/cellmask.o EFDC/cbalev.o \
	  EFDC/calhta.o EFDC/calebi.o EFDC/MPI_Mapping/Map_Jet_Plume.o \
	  EFDC/MPI_Domain_Decomp/Allocate_Domain_Decomp.o EFDC/wasp7.o EFDC/rvelplth.o EFDC/runcontrol.o \
	  EFDC/caltsxy.o EFDC/caltranice.o EFDC/calexp2t.o EFDC/bedload.o EFDC/MPI_Mapping/Map_River.o \
	  EFDC/MPI_Mapping/Map_OpenBC_Pressure.o EFDC/MPI_Domain_Decomp/Parent_Grid.o \
	  EFDC/wasp5.o EFDC/svdcmp.o EFDC/s_main.o EFDC/calstepd.o EFDC/calmmt.o EFDC/bedinit.o \
	  EFDC/Propwash/Calc_Prop_Erosion.o EFDC/MPI_Out/Write_Mapping_3Dgraphics.o EFDC/wasp7hydro.o EFDC/varzerosnl.o EFDC/s_massbalance.o \
	  EFDC/hdmt2t.o EFDC/hdmt.o EFDC/calqq2t.o EFDC/calhdmf.o EFDC/caldye.o EFDC/Utilities/seek.o \
	  EFDC/MPI_Mapping/ReMap_RSSBC.o EFDC/MPI_Mapping/Map_OpenBC_Conc.o \
	  EFDC/rsalplth.o EFDC/calexp.o EFDC/calconc.o EFDC/wavesxy.o \
	  EFDC/setopenbc.o EFDC/rvelpltv.o EFDC/congrad.o EFDC/calsft.o EFDC/calavb.o EFDC/Utilities/dstime.o \
	  EFDC/MPI_Utilities/Setup_MPI_Debug_File.o \
	  EFDC/input.o EFDC/MPI_Mapping/Setup_Local_to_Global.o \
	  EFDC/s_shear.o EFDC/rout3d.o EFDC/foodchain.o EFDC/caltran.o \
	  EFDC/wasp8hydro.o \
	  EFDC/partmix.o EFDC/congradc.o EFDC/caltbxy.o  EFDC/calpuv9c.o \
	  EFDC/caldisp3.o \
	  
OBJSLINK= Setup_Ships.o aaefdc.o csedress.o fsedmode.o fdstrse.o FormatOutputMPI.o csedvis.o csndset.o fsbdld.o \
       csndeqc.o fstrse.o svbksb.o setshld.o csedresb.o csedtaus.o fprobdep.o \
       setstvel.o csndzeq.o csedset.o ludcmp.o fhydcn.o csedtaub.o lubksb.o budget.o \
       welcome.o \
	  subchan.o rsalpltv.o calvegser.o caltox_kinetics.o \
	  caldiff.o bankerosed.o Propwash_Calc_Sequence.o Write_Array_Sizes.o bal2t.o Report_Max_Min_Timing.o Communicate_Groups.o \
	  s_slope.o s_bedload.o caltox.o \
	  calsnd.o caldisp2.o Setup_MapBackToGlobal.o Scan_Cell.o \
	  s_sedic.o rcahq.o cbalod.o calpuv2c.o calbuoy.o Scan_JSON_Decomp.o Calculate_Local_Shear.o \
	  wavebl.o showval.o out3d.o \
	  calsed.o calqq1.o ainit.o \
	  wasp4.o varinit.o rsurfplt.o cellmap.o \
	  calfqc.o Det_Adjacent_Cells.o wwqnc.o calbed9.o Setup_MPI_Topology.o \
	  Map_Connectors.o \
		lsqharm.o caltran_ad.o caltoxb.o \
	  calpser.o calbed.o wasp6.o \
	  ceqicm.o calhdmf3.o calblay.o Initialize_MPI.o \
	  Child_Grid.o Communicate_3D_0.o \
	  varalloc.o s_sedzlj.o negdep.o calbal.o \
	  Create_List_No_Ghost_Cells.o Read_JSON_Decomp.o \
	  ssedtox.o s_tecplot.o setbcs.o mhkpwr.o caluvw.o \
	  calqvs.o calpnhs.o Congrad_MPI.o \
	  tmsr.o timelog.o setfpocb.o cellmask.o cbalev.o \
	  calhta.o calebi.o Map_Jet_Plume.o \
	  Allocate_Domain_Decomp.o wasp7.o rvelplth.o runcontrol.o \
	  caltsxy.o caltranice.o calexp2t.o bedload.o Map_River.o \
	  Map_OpenBC_Pressure.o Parent_Grid.o \
	  wasp5.o svdcmp.o s_main.o calstepd.o calmmt.o bedinit.o \
	  Calc_Prop_Erosion.o Write_Mapping_3Dgraphics.o wasp7hydro.o varzerosnl.o s_massbalance.o \
	  hdmt2t.o hdmt.o calqq2t.o calhdmf.o caldye.o seek.o \
	  ReMap_RSSBC.o Map_OpenBC_Conc.o \
	  rsalplth.o calexp.o calconc.o wavesxy.o \
	  setopenbc.o rvelpltv.o congrad.o calsft.o calavb.o dstime.o \
	  Setup_MPI_Debug_File.o \
	  input.o Setup_Local_to_Global.o \
	  s_shear.o rout3d.o foodchain.o caltran.o wasp8hydro.o \
	  partmix.o congradc.o caltbxy.o calpuv9c.o \
	  caldisp3.o \
  	  	 	  	  	 	 	  	 							 
MODS    = EFDC/JSON_Reader/fson_string_m.o EFDC/mod_var_global.o EFDC/MPI_Utilities/Variables_MPI_Mapping.o EFDC/Utilities/convertwgs84.o EFDC/MPI_Utilities/Variables_MPI_MapGatherSort.o \
              EFDC/Propwash/Mod_Erosive_Flux.o EFDC/Utilities/mod_julian.o EFDC/Propwash/Variables_Propwash.o EFDC/Drifter/Variables_MPI_Drifter.o \
	      EFDC/Utilities/mod_info.o EFDC/MPI_Utilities/Variables_MPI_Write_Out.o EFDC/Utilities/mod_allocate.o EFDC/Eutrophication/mod_wq_vars.o \
              EFDC/wavelength.o EFDC/Propwash/Mod_Ship.o EFDC/JSON_Reader/fson_value_m.o EFDC/MPI_Utilities/Variables_MPI.o EFDC/MPI_Mapping/Mod_Map_Soln.o \
		  EFDC/MPI_Mapping/Mod_Assign_Loc_Glob_For_Write.o EFDC/MPI_Out/Mod_Write_Cell_Map_.o EFDC/MPI_Communication/Communicate_Ghost_Routines.o \
		  EFDC/MPI_Communication/Mod_SEND_RECV.o EFDC/mod_diffuser.o EFDC/mod_calcser.o EFDC/mod_hydstructure.o EFDC/Eutrophication/mod_macrohydro.o \
                  EFDC/MPI_Communication/Mod_Broadcast_Routines.o EFDC/JSON_Reader/fson_path_m.o EFDC/Utilities/xyijconv.o \
		  EFDC/MPI_Communication/Mod_Gather_Drifter.o EFDC/MPI_Utilities/Mod_MPI_Helper_Functions.o EFDC/MPI_Mapping/Mod_Sort_Global_Soln.o \
		  EFDC/MPI_Communication/Mod_Gather_Soln.o EFDC/MPI_Communication/mod_allreduce.o EFDC/JSON_Reader/fson.o EFDC/Eutrophication/mod_shellfish.o EFDC/mod_fields.o \
		  EFDC/Propwash/Mod_Position.o EFDC/MPI_Mapping/Mod_Map_Gather_Sort.o EFDC/mod_cyclone.o EFDC/Propwash/Mod_Mesh.o EFDC/Eutrophication/mod_zoopl.o \
		  EFDC/MPI_Mapping/Mod_Map_Global_to_Local.o  EFDC/Eutrophication/mod_diagen.o EFDC/Propwash/Mod_Position_Cell.o EFDC/Propwash/Mod_Position_Elev.o \
                  EFDC/Propwash/Mod_All_Tracks.o EFDC/Eutrophication/mod_rpem.o EFDC/Propwash/Mod_All_Ship_Tracks.o EFDC/Propwash/Mod_Active_Ship.o \
                  EFDC/mod_scaninp.o EFDC/hifreqout.o EFDC/MPI_Out/Mod_Map_Write_EE_Binary.o EFDC/Propwash/Variables_Ship.o \
		  EFDC/mod_restart.o EFDC/MPI_Communication/Communicate_Drifters.o EFDC/Drifter/mod_drifter.o EFDC/mod_windwave.o EFDC/Eutrophication/mod_wq.o \
                  EFDC/Propwash/Mod_Read_Propwash.o EFDC/mod_heat.o \
		  EFDC/MOD_efdcout.o EFDC/mod_getswan.o EFDC/MPI_Out/Mod_Map_Write_NetCDF.o EFDC/mod_netcdf.o \
		  EFDC/MPI_Utilities/Variables_MPI_Concentration.o \

MODSLINK  = fson_string_m.o mod_var_global.o Variables_MPI_Mapping.o convertwgs84.o Variables_MPI_MapGatherSort.o Mod_Erosive_Flux.o \
	      mod_julian.o Variables_Propwash.o Variables_MPI_Drifter.o mod_info.o Variables_MPI_Write_Out.o \
              mod_allocate.o mod_wq_vars.o wavelength.o Mod_Ship.o fson_value_m.o Variables_MPI.o Mod_Map_Soln.o Mod_Assign_Loc_Glob_For_Write.o \
		  Mod_Write_Cell_Map_.o Communicate_Ghost_Routines.o Mod_SEND_RECV.o mod_diffuser.o mod_calcser.o mod_hydstructure.o \
		  mod_macrohydro.o Mod_Broadcast_Routines.o fson_path_m.o xyijconv.o \
		  Mod_Gather_Drifter.o Communicate_Drifters.o Mod_MPI_Helper_Functions.o Mod_Sort_Global_Soln.o Mod_Gather_Soln.o mod_allreduce.o \
		  fson.o mod_shellfish.o mod_fields.o Mod_Position.o Mod_Map_Gather_Sort.o mod_cyclone.o \
                  Mod_Mesh.o mod_zoopl.o Mod_Map_Global_to_Local.o mod_diagen.o Mod_Position_Cell.o Mod_Position_Elev.o \
                  Mod_All_Tracks.o mod_rpem.o Mod_All_Ship_Tracks.o Mod_Active_Ship.o mod_scaninp.o hifreqout.o \
                  Mod_Map_Write_EE_Binary.o Variables_Ship.o \
		  mod_drifter.o MOD_efdcout.o Mod_Map_Write_NetCDF.o mod_restart.o mod_windwave.o mod_netcdf.o mod_wq.o Mod_Read_Propwash.o mod_heat.o \
		  mod_getswan.o Variables_MPI_Concentration.o \

MAIN    = 


EXEPROG = $(PROGRAM)$(PROGEXT)


###### FOR LINUX #######################################################

$(EXEPROG) : $(OBJS) $(MODS) $(MAIN) 
  
	$(FC) $(FFLAGS) $(CDFDIR) $(CDFLIBS) $(INCLUD) -o  $@  $(OBJSLINK) $(MODSLINK) $(MAIN)


# Recompile only source file of change #################################

.f90.o: 
	$(FC) $(FFLAGS) $(INCLUD) -c $<

.F90.o: 
	$(FC) $(FFLAGS) $(INCLUD) -c $<

.for.o: 
	$(FC) $(FFLAGS) $(INCLUD) -c $<

.FOR.o: 
	$(FC) $(FFLAGS) $(INCLUD) -c $<

### Special dependencies ###############################################

$(OBJS)              : $(MODS)
$(MAIN)              : $(OBJS)
 
clean:
	rm -f *.o *.mod

# Make dependency list
#depend .depend :
#	makedepf90 -W -DENABLE_MPI -o aaefdc $(SRCS) > .depend

