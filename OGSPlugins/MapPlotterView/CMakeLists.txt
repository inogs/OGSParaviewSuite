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


INCLUDE_DIRECTORIES(${PARAVIEW_GUI_INCLUDE_DIRS})
#INCLUDE (ParaViewQt5)
QT5_WRAP_CPP(MOC_SRCS MapPlotterView.h)

ADD_PQPROXY(VIEW_IFACE VIEW_IFACE_SRCS
    TYPE MapPlotterView
    XML_GROUP views
    XML_NAME MapPlotterView
    )

ADD_PARAVIEW_PLUGIN(MapPlotterView "1.0"
  SERVER_MANAGER_XML MapPlotterView.xml
  SERVER_MANAGER_SOURCES vtkSMMapPlotterViewProxy.cxx vtkMapPlotterView.cxx
  REQUIRED_ON_SERVER
  GUI_INTERFACES
    ${IFACES} ${VIEW_IFACE}
  GUI_SOURCES
    ${IFACE_SRCS} ${VIEW_IFACE_SRCS} ${UI_SRCS}
    ${MOC_SRCS}
    MapPlotterView.cxx
  )