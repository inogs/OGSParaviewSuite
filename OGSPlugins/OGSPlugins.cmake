##############################################
###     OGSPLUGINS CMAKE FOR PARAVIEW      ###
##############################################

# Variables
set( VECTORIZATION ON)
set( OPENMP_PARALL ON)
set( COMPILER      "GCC") # Either GCC or INTEL
set( DEBUGGING     OFF)
set( DOUBLEARRAYS  ON) # Use double arrays instead of float arrays

set( MAIN_FOLDER "${CMAKE_CURRENT_LIST_DIR}")

# Required external packages
find_package(ospray REQUIRED) # ParaView 5.6.3 fix

# Optimized compiling flags
if(NOT "${COMPILER}" STREQUAL "INTEL")
	# Debugging flags
	if ( ${DEBUGGING} )
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O0 -g -rdynamic") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -rdynamic")
	else()
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O2")
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
	endif()
	# Vectorization flags
	if ( ${VECTORIZATION} )
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -ffast-math -march=native -ftree-vectorize") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -march=native -ftree-vectorize")
	endif()
	# OpenMP flag
	if ( ${OPENMP_PARALL} )
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fopenmp -DUSE_OMP") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DUSE_OMP")
	endif()
else()
	# Debugging flags
	if ( ${DEBUGGING} )
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O0 -g -traceback") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -traceback")
	else()
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O2") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
	endif()
	# Vectorization flags
	if ( ${VECTORIZATION} )
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -xHost -mtune=skylake") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -mtune=skylake")
	endif()
	# OpenMP flag
	if ( ${OPENMP_PARALL} )
		set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -qopenmp -DUSE_OMP") 
		set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -DUSE_OMP")
	endif()
endif()
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Add the flag to use double arrays instead of float arrays
if ( ${DOUBLEARRAYS} )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVTK_USE_DOUBLE")
endif()

# Now define a set of variables containing various code
# that might be linked when compiling
set( UTILS  "${MAIN_FOLDER}/_utils/utils.cmake")         # General Utilities
set( PROJ   "${MAIN_FOLDER}/_utils/proj/proj.cmake")     # Proj library for projections
set( PUGI   "${MAIN_FOLDER}/_utils/pugixml/pugi.cmake")  # PUGI library for XML parsing
set( CNPY   "${MAIN_FOLDER}/_utils/cnpy/cnpy.cmake")     # CNPY library for numpy NPZ writing
set( LAPACK "${MAIN_FOLDER}/_utils/lapack/lapack.cmake") # LAPACK library for matrix operations
#set( VTK    "${MAIN_FOLDER}/VTK/vtk.cmake")


#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
include_directories(${VTK_INCLUDE_DIRS})

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

if (PARAVIEW_USE_MPI)
  include(vtkMPI)
endif()

# Set a consistent MACOSX_RPATH default across all CMake versions.
# When CMake 2.8.12 is required, change this default to 1.
# When CMake 3.0.0 is required, remove this block (see CMP0042).
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

