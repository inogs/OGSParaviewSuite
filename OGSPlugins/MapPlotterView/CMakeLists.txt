# create a plugin that adds a reader to the ParaView GUI
# it is added in the file dialog when doing opens/saves.
cmake_minimum_required(VERSION 3.3)

#--------------------------------------------------
# Main CMAKE file detection
#--------------------------------------------------
include(./../OGSPlugins.cmake)

#--------------------------------------------------
# Libraries
#--------------------------------------------------


#--------------------------------------------------
# GUI Proxis
#--------------------------------------------------
include_directories(${PARAVIEW_GUI_INCLUDE_DIRS})
qt5_wrap_cpp(MOC_SRCS MapPlotterView.h)

add_pqproxy(VIEW_IFACE VIEW_IFACE_SRCS
    TYPE      MapPlotterView
    XML_GROUP views
    XML_NAME  MapPlotterView
)

#--------------------------------------------------
# ParaView plugin
#--------------------------------------------------
add_paraview_plugin(MapPlotterView "1.0"
  SERVER_MANAGER_XML     MapPlotterView.xml
  SERVER_MANAGER_SOURCES vtkSMMapPlotterViewProxy.cxx vtkMapPlotterView.cxx
  REQUIRED_ON_SERVER
  GUI_INTERFACES         ${IFACES} ${VIEW_IFACE}
  GUI_SOURCES
    ${IFACE_SRCS} ${VIEW_IFACE_SRCS} ${UI_SRCS}
    ${MOC_SRCS}
    MapPlotterView.cxx
)