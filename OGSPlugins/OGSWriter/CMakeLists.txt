# create a plugin that adds a reader to the ParaView GUI
# it is added in the file dialog when doing opens/saves.

cmake_minimum_required(VERSION 3.3)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})

IF(NOT ParaView_BINARY_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF()

# Set a consistent MACOSX_RPATH default across all CMake versions.
# When CMake 2.8.12 is required, change this default to 1.
# When CMake 3.0.0 is required, remove this block (see CMP0042).
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

if(PARAVIEW_USE_MPI)
  include(vtkMPI)
endif()

if(NOT APPLE)
# OpenMP does not work with Apple's clang so we deactivate it.
	set ( CMAKE_C_FLAGS   "-fopenmp")
	set ( CMAKE_CXX_FLAGS "-fopenmp")
endif()

set ( CMAKE_C_FLAGS   "-lz")
set ( CMAKE_CXX_FLAGS "-lz")

ADD_PARAVIEW_PLUGIN(OGSWriter "1.0"
  SERVER_MANAGER_XML OGSWriter.xml
  SERVER_MANAGER_SOURCES ../_utils/OGS.cpp ../_utils/netcdfio.cpp ../_utils/vtkFields.cpp ../_utils/vtkOperations.cpp cnpy.cpp vtkOGSNPZWriter.cxx vtkOGSNPZTableWriter.cxx vtkOGSFieldWriter.cxx vtkOGSNetCDFWriter.cxx
  REQUIRED_ON_SERVER
  )