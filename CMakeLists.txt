cmake_minimum_required(VERSION 3.15)

# Project name and languages (C needed for HDF5/NetCDF detection)
project(EFDC LANGUAGES C Fortran)

option(USE_GFORTRAN "Use gfortran instead of Intel Fortran compiler" OFF)
option(ENABLE_PROFILING "Enable profiling support" OFF)
option(ENABLE_FAST_MATH "Enable fast math (may reduce accuracy slightly)" OFF)
option(USE_STATIC_NETCDF "Link NetCDF libraries statically" OFF)
option(FORCE_STATIC_SYSTEM_LIBS "Force static linking of system libraries" OFF)

# Allow custom executable name
set(CUSTOM_EXE_NAME "" CACHE STRING "Custom executable name (without extension). Default: efdc")

# Set build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Detect CPU features
# Detect SUSE Linux for special handling
if(EXISTS "/etc/SuSE-release" OR EXISTS "/etc/SUSE-brand")
    set(IS_SUSE_LINUX TRUE)
    message(STATUS "Detected SUSE Linux - applying compatibility fixes")
endif()

# Set compilers and optimization flags
if(USE_GFORTRAN)
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_Fortran_COMPILER "gfortran")
    
    # GFortran base flags
    set(BASE_FLAGS "-fopenmp -ffree-form -ffree-line-length-none -freal-4-real-8")
    set(PARALLEL_FLAGS "-DENABLE_MPI -D_MPI -DLINUX")
    
    # Build-type specific flags for GFortran
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(DEBUG_FLAGS "-g -O0 -fcheck=all -fbacktrace -Wall -Wextra -Wconversion -pedantic")
        set(DEBUG_FLAGS "${DEBUG_FLAGS} -ffpe-trap=invalid,zero,overflow -finit-real=snan")
        set(CMAKE_Fortran_FLAGS "${BASE_FLAGS} ${DEBUG_FLAGS} ${PARALLEL_FLAGS}")
        message(STATUS "GFortran Debug build configured with extensive checking")
    else()
        # Release flags - balanced speed/accuracy
        set(SAFETY_FLAGS "-fcheck=bounds -Wall -Wextra -Wconversion")
        set(PERFORMANCE_FLAGS "-O3 -funroll-loops -march=native -mtune=native -flto")
        
        if(ENABLE_FAST_MATH)
            set(PERFORMANCE_FLAGS "${PERFORMANCE_FLAGS} -ffast-math -fno-protect-parens")
        endif()
        
        if(ENABLE_PROFILING)
            set(PERFORMANCE_FLAGS "${PERFORMANCE_FLAGS} -pg")
        endif()
        
        set(CMAKE_Fortran_FLAGS "${BASE_FLAGS} ${SAFETY_FLAGS} ${PERFORMANCE_FLAGS} ${PARALLEL_FLAGS}")
    endif()

else()
    # Intel Compiler Setup - detect from environment if available
    if(DEFINED ENV{CC} AND "$ENV{CC}" MATCHES "icx")
        set(CMAKE_C_COMPILER "$ENV{CC}")
        message(STATUS "Using Intel C compiler from environment: $ENV{CC}")
    else()
        set(CMAKE_C_COMPILER "icx")
        message(STATUS "Using default Intel C compiler: icx")
    endif()
    
    if(DEFINED ENV{FC} AND "$ENV{FC}" MATCHES "ifx")
        set(CMAKE_Fortran_COMPILER "$ENV{FC}")
        message(STATUS "Using Intel Fortran compiler from environment: $ENV{FC}")
    elseif(DEFINED ENV{FC} AND "$ENV{FC}" MATCHES "mpiifx")
        set(CMAKE_Fortran_COMPILER "$ENV{FC}")
        message(STATUS "Using Intel MPI Fortran compiler from environment: $ENV{FC}")
    else()
        # Try to find mpiifx first, then fall back to ifx
        find_program(MPIIFX_COMPILER 
            NAMES mpiifx
            PATHS /opt/intel/oneapi/mpi/latest/bin $ENV{I_MPI_ROOT}/bin
            NO_DEFAULT_PATH
        )
        if(MPIIFX_COMPILER)
            set(CMAKE_Fortran_COMPILER "${MPIIFX_COMPILER}")
            message(STATUS "Found Intel MPI Fortran compiler: ${MPIIFX_COMPILER}")
        else()
            set(CMAKE_Fortran_COMPILER "ifx")
            message(STATUS "Using default Intel Fortran compiler: ifx")
        endif()
    endif()
    
    # Intel base flags
    set(BASE_FLAGS "-fiopenmp -pthread -real-size 64")
    set(PREPROCESS_FLAGS "-fpp -DLINUX -DENABLE_MPI -D_MPI")
    
    # Build-type specific flags for Intel
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(DEBUG_FLAGS "-g -O0 -debug all -traceback -check all")
        set(DEBUG_FLAGS "${DEBUG_FLAGS} -warn shape -stand f18 -fpe0 -init=snan,arrays -ftrapuv")
        set(CMAKE_Fortran_FLAGS "${BASE_FLAGS} ${PREPROCESS_FLAGS} -DDEBUGGING ${DEBUG_FLAGS}")
        message(STATUS "Intel Debug build configured with comprehensive debugging")
    else()
        # Release flags - MODIFIED for SUSE compatibility
        set(SAFETY_FLAGS "-traceback -fp-model=source -fp-speculation=safe -fpe0")
        set(PERFORMANCE_FLAGS "-O3 -xSSE4.2 -ipo -unroll -no-prec-div")
        set(PROCESS_FLAGS "-multiple-processes=12")
        
        # UPDATED: More conservative static linking for Intel on SUSE
        if(USE_STATIC_NETCDF)
            set(STATIC_FLAGS "-static-intel")
            
            # Add Intel-specific workarounds for SUSE
            if(IS_SUSE_LINUX)
                set(STATIC_FLAGS "${STATIC_FLAGS} -Wl,--allow-multiple-definition")
            endif()
        endif()
            
        if(ENABLE_PROFILING)
            set(PERFORMANCE_FLAGS "${PERFORMANCE_FLAGS} -profile-functions -profile-loops=all")
        endif()
        
        # Memory optimization for large simulations
        set(MEMORY_FLAGS "-heap-arrays 64 -assume byterecl -assume buffered_io")
        
        set(CMAKE_Fortran_FLAGS "${BASE_FLAGS} ${PREPROCESS_FLAGS} ${SAFETY_FLAGS} ${PERFORMANCE_FLAGS} ${PROCESS_FLAGS} ${MEMORY_FLAGS} ${STATIC_FLAGS}")
    endif()
endif()

# Define the executable name with build type suffix
if(CUSTOM_EXE_NAME)
    set(PROGRAM_NAME "${CUSTOM_EXE_NAME}")
else()
    set(PROGRAM_NAME "efdc")
endif()

# Add build type to executable name for debug builds
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(EXE_NAME "${PROGRAM_NAME}_debug")
else()
    set(EXE_NAME "${PROGRAM_NAME}")
endif()

# NetCDF paths integration with modular setup system
set(NETCDF_SEARCH_PATHS)

# First priority: Environment variables from setup_netcdf_env.sh script
if(DEFINED ENV{NETCDF_FORTRAN_ROOT})
    list(APPEND NETCDF_SEARCH_PATHS "$ENV{NETCDF_FORTRAN_ROOT}")
    message(STATUS "Using NetCDF-Fortran from environment: $ENV{NETCDF_FORTRAN_ROOT}")
endif()
if(DEFINED ENV{NETCDF_ROOT})
    list(APPEND NETCDF_SEARCH_PATHS "$ENV{NETCDF_ROOT}")
    message(STATUS "Using NetCDF-C from environment: $ENV{NETCDF_ROOT}")
endif()
if(DEFINED ENV{HDF5_ROOT})
    list(APPEND NETCDF_SEARCH_PATHS "$ENV{HDF5_ROOT}")
    message(STATUS "Using HDF5 from environment: $ENV{HDF5_ROOT}")
endif()

# Second priority: Default paths from setup.sh script
if(DEFINED ENV{HOME})
    list(APPEND NETCDF_SEARCH_PATHS 
        "$ENV{HOME}/software/netcdf-fortran"
        "$ENV{HOME}/software/netcdf-c" 
        "$ENV{HOME}/software/hdf5"
    )
    message(STATUS "Checking setup.sh default paths: $ENV{HOME}/software/")
endif()

# Third priority: Standard system paths
list(APPEND NETCDF_SEARCH_PATHS "/usr/local" "/opt/netcdf" "/usr")

# Check if we found setup.sh-built libraries
if(EXISTS "$ENV{HOME}/software/setup_netcdf_env.sh")
    message(STATUS "Found setup.sh environment script: $ENV{HOME}/software/setup_netcdf_env.sh")
    message(STATUS "Hint: Source this script before running cmake for best results")
else()
    message(STATUS "Note: setup.sh environment script not found. Run setup/setup.sh first for optimal configuration.")
endif()

# Set CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH ${NETCDF_SEARCH_PATHS})

# Force static library preference when USE_STATIC_NETCDF is ON
if(USE_STATIC_NETCDF)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
    message(STATUS "Static NetCDF linking enabled - preferring .a libraries")
endif()

# Find required packages
find_package(PkgConfig REQUIRED)

# Manual detection for NetCDF-C
find_path(NetCDF_C_INCLUDE_DIR netcdf.h
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

find_library(NetCDF_C_LIBRARY netcdf
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
)

if(NetCDF_C_INCLUDE_DIR AND NetCDF_C_LIBRARY)
    set(NetCDF_C_FOUND TRUE)
    set(NetCDF_C_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIR})
    set(NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY})
    get_filename_component(NetCDF_C_LIBRARY_DIR ${NetCDF_C_LIBRARY} DIRECTORY)
    set(NetCDF_C_LIBRARY_DIRS ${NetCDF_C_LIBRARY_DIR})
endif()

# Manual detection for NetCDF-Fortran
find_path(NetCDF_Fortran_INCLUDE_DIR netcdf.mod
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

find_library(NetCDF_Fortran_LIBRARY netcdff
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
)

if(NetCDF_Fortran_INCLUDE_DIR AND NetCDF_Fortran_LIBRARY)
    set(NetCDF_Fortran_FOUND TRUE)
    set(NetCDF_Fortran_INCLUDE_DIRS ${NetCDF_Fortran_INCLUDE_DIR})
    set(NetCDF_Fortran_LIBRARIES ${NetCDF_Fortran_LIBRARY})
    get_filename_component(NetCDF_Fortran_LIBRARY_DIR ${NetCDF_Fortran_LIBRARY} DIRECTORY)
    set(NetCDF_Fortran_LIBRARY_DIRS ${NetCDF_Fortran_LIBRARY_DIR})
endif()

# Manual HDF5 detection
find_path(HDF5_INCLUDE_DIR hdf5.h
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

find_library(HDF5_C_LIBRARY hdf5
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
)

find_library(HDF5_HL_LIBRARY hdf5_hl
    HINTS ${NETCDF_SEARCH_PATHS}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
)

if(USE_STATIC_NETCDF)
    find_library(HDF5_Fortran_LIBRARY hdf5_fortran
        HINTS ${NETCDF_SEARCH_PATHS}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH
    )
    
    find_library(HDF5_HL_Fortran_LIBRARY hdf5_hl_fortran
        HINTS ${NETCDF_SEARCH_PATHS}
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH
    )
endif()

# Set HDF5 variables
if(HDF5_INCLUDE_DIR AND HDF5_C_LIBRARY AND HDF5_HL_LIBRARY)
    set(HDF5_FOUND TRUE)
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    
    set(HDF5_LIBRARIES)
    if(USE_STATIC_NETCDF AND HDF5_Fortran_LIBRARY)
        list(APPEND HDF5_LIBRARIES ${HDF5_Fortran_LIBRARY})
    endif()
    if(USE_STATIC_NETCDF AND HDF5_HL_Fortran_LIBRARY)
        list(APPEND HDF5_LIBRARIES ${HDF5_HL_Fortran_LIBRARY})
    endif()
    list(APPEND HDF5_LIBRARIES ${HDF5_HL_LIBRARY} ${HDF5_C_LIBRARY})
    
    get_filename_component(HDF5_LIBRARY_DIR ${HDF5_C_LIBRARY} DIRECTORY)
    set(HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIR})
else()
    message(FATAL_ERROR "Could not find HDF5 libraries")
endif()

# Error if dependencies not found
if(NOT NetCDF_C_FOUND)
    message(FATAL_ERROR "NetCDF-C not found! Please build NetCDF first.")
endif()

if(NOT NetCDF_Fortran_FOUND)
    message(FATAL_ERROR "NetCDF-Fortran not found! Please build NetCDF-Fortran first.")
endif()

# Find libunwind for debug builds only
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_library(LIBUNWIND_LIBRARY unwind
        PATHS /usr/lib64 /usr/lib /lib64 /lib /usr/local/lib
        PATH_SUFFIXES x86_64-linux-gnu
        REQUIRED
    )
    message(STATUS "Debug mode: Found libunwind: ${LIBUNWIND_LIBRARY}")
endif()

# UPDATED: Find system libraries for static linking
if(USE_STATIC_NETCDF AND FORCE_STATIC_SYSTEM_LIBS)
    # Try to find static system libraries
    find_library(ZLIB_STATIC_LIBRARY 
        NAMES libz.a z
        PATHS /usr/lib64 /usr/lib /lib64 /lib
        PATH_SUFFIXES x86_64-linux-gnu
    )
    
    find_library(MATH_STATIC_LIBRARY
        NAMES libm.a m
        PATHS /usr/lib64 /usr/lib /lib64 /lib
        PATH_SUFFIXES x86_64-linux-gnu
    )
    
    if(NOT ZLIB_STATIC_LIBRARY)
        message(WARNING "Static zlib not found. Install zlib-devel-static or libz-dev")
        set(SYSTEM_LIBS z m)  # Fall back to dynamic
    else()
        set(SYSTEM_LIBS ${ZLIB_STATIC_LIBRARY} ${MATH_STATIC_LIBRARY})
    endif()
else()
    # Use dynamic system libraries
    set(SYSTEM_LIBS z m)
endif()

# Include directories
include_directories(${NetCDF_C_INCLUDE_DIRS})
include_directories(${NetCDF_Fortran_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(/usr/local/include/WASP)

# Add source files
file(GLOB_RECURSE SOURCES "*.f90" "*/**.f90")

# Filter out excluded files
list(FILTER SOURCES EXCLUDE REGEX ".*/test/.*")
list(FILTER SOURCES EXCLUDE REGEX ".*/examples/.*")

set(EXCLUDED_FILES
    "EFDC/GOTM_Turbulence/variances.f90"
    "EFDC/GOTM_Util/diff_center.f90"
    "EFDC/MPI_Communication/communicate_w.f90"
    "EFDC/MPI_Utilities/task_division.f90"
    "EFDC/GOTM_Turbulence/cmue_a.f90"
    "EFDC/GOTM_Turbulence/r_ratio.f90"
)
foreach(EXCLUDED_FILE ${EXCLUDED_FILES})
    list(FILTER SOURCES EXCLUDE REGEX ".*${EXCLUDED_FILE}$")
endforeach()

# Add precompiled objects/modules
set(PRECOMPILED_OBJECTS)
set(PRECOMPILED_MODULES)

# Add the executable target
add_executable(${EXE_NAME} ${SOURCES} ${PRECOMPILED_OBJECTS} ${PRECOMPILED_MODULES})

# UPDATED: Improved linking strategy for SUSE Linux
if(USE_STATIC_NETCDF)
    # Approach 1: Conservative static linking (recommended for SUSE)
    target_link_libraries(${EXE_NAME} PRIVATE 
        ${NetCDF_Fortran_LIBRARIES}
        ${NetCDF_C_LIBRARIES}
        ${HDF5_LIBRARIES}
        ${SYSTEM_LIBS}
    )
    
    # SUSE-specific linker options
    if(IS_SUSE_LINUX)
        target_link_options(${EXE_NAME} PRIVATE
            "LINKER:--allow-multiple-definition"
            "LINKER:--as-needed"
        )
    endif()
    
    # Set static linking properties (but not for debug builds to avoid issues)
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        set_target_properties(${EXE_NAME} PROPERTIES
            LINK_SEARCH_START_STATIC ON
            LINK_SEARCH_END_STATIC ON
        )
    endif()
    
    message(STATUS "Configured for static NetCDF linking with SUSE compatibility")
else()
    # Dynamic linking
    target_link_libraries(${EXE_NAME} PRIVATE 
        ${NetCDF_Fortran_LIBRARIES}
        ${NetCDF_C_LIBRARIES}
        ${HDF5_LIBRARIES}
    )
    
    set_target_properties(${EXE_NAME} PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
    
    message(STATUS "Configured for dynamic NetCDF linking")
endif()

# Add libunwind for debug builds
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_libraries(${EXE_NAME} PRIVATE ${LIBUNWIND_LIBRARY})
endif()

# Add library directories
target_link_directories(${EXE_NAME} PRIVATE
    ${NetCDF_Fortran_LIBRARY_DIRS}
    ${NetCDF_C_LIBRARY_DIRS}
    ${HDF5_LIBRARY_DIRS}
)

# Debug-specific target properties
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Disable optimizations that can interfere with debugging
    set_target_properties(${EXE_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION OFF
    )
    
    # Add debug information to the executable
    target_compile_options(${EXE_NAME} PRIVATE -g)
    target_link_options(${EXE_NAME} PRIVATE -g)
endif()

# Clean target
add_custom_target(clean-all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMAND ${CMAKE_COMMAND} -E remove -f ${EXE_NAME}
    COMMAND ${CMAKE_COMMAND} -E remove -f ${PROGRAM_NAME}.x ${PROGRAM_NAME}_debug.x
    COMMAND ${CMAKE_COMMAND} -E remove -f efdc.x efdc_debug.x
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove -f CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove -f cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove -f Makefile
    COMMENT "Cleaning all build files"
)

# Check linking target
add_custom_target(check-linking
    COMMAND ldd ${EXE_NAME} || echo "Statically linked binary"
    DEPENDS ${EXE_NAME}
    COMMENT "Checking library dependencies"
)

# Debug convenience targets
add_custom_target(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build .
    COMMENT "Building debug version"
)

add_custom_target(release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build .
    COMMENT "Building release version"
)

# Print configuration summary
message(STATUS "")
message(STATUS "EFDC Configuration Summary:")
message(STATUS "  System: ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "  Flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "  Executable: ${EXE_NAME}")
if(CUSTOM_EXE_NAME)
    message(STATUS "  Custom Name: ${CUSTOM_EXE_NAME} (default: efdc)")
endif()
message(STATUS "  Use GFortran: ${USE_GFORTRAN}")
message(STATUS "  Static NetCDF: ${USE_STATIC_NETCDF}")
message(STATUS "  Force Static System Libs: ${FORCE_STATIC_SYSTEM_LIBS}")
if(IS_SUSE_LINUX)
    message(STATUS "  SUSE Linux detected: Applied compatibility fixes")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "  Debug build: Extensive checking and debugging enabled")
    message(STATUS "  libunwind: ${LIBUNWIND_LIBRARY}")
endif()
if(USE_STATIC_NETCDF)
    message(STATUS "  NetCDF-Fortran: ${NetCDF_Fortran_LIBRARY}")
    message(STATUS "  NetCDF-C: ${NetCDF_C_LIBRARY}")
    message(STATUS "  HDF5: ${HDF5_C_LIBRARY}")
    message(STATUS "  System libs: ${SYSTEM_LIBS}")
endif()
message(STATUS "")