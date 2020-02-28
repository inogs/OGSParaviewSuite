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

# UTILS detection
include(../_utils/utils.cmake)
# PROJ detection
include(../_utils/proj/proj.cmake)

if(NOT APPLE)
# OpenMP does not work with Apple's clang so we deactivate it.
  set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fopenmp")
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# Fix for ParaView-5.6.3 and OSPRAY
find_package(ospray REQUIRED)

ADD_PARAVIEW_PLUGIN(OGSVortexDetection "1.0"
  SERVER_MANAGER_XML OGSVortexDetection.xml
  SERVER_MANAGER_SOURCES ${VTKOPS} ${VTKFIELDS} vtkOGSVortexDetection.cxx
  REQUIRED_ON_SERVER
  )
