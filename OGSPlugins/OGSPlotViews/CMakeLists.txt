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
qt5_wrap_cpp(MOC_SRCS OGSVerticalProfilePlot.h OGSHovmoellerPlot.h OGSSpaghettiPlot.h OGSMapPlot.h)

add_pqproxy(VIEW_IFACE VIEW_IFACE_SRCS
    TYPE      OGSVerticalProfilePlot
    XML_GROUP views
    XML_NAME  OGSVerticalProfilePlot
)

add_pqproxy(VIEW_IFACE VIEW_IFACE_SRCS
    TYPE      OGSHovmoellerPlot
    XML_GROUP views
    XML_NAME  OGSHovmoellerPlot
)

add_pqproxy(VIEW_IFACE VIEW_IFACE_SRCS
    TYPE      OGSSpaghettiPlot
    XML_GROUP views
    XML_NAME  OGSSpaghettiPlot
)

add_pqproxy(VIEW_IFACE VIEW_IFACE_SRCS
    TYPE      OGSMapPlot
    XML_GROUP views
    XML_NAME  OGSMapPlot
)

#--------------------------------------------------
# ParaView plugin
#--------------------------------------------------
add_paraview_plugin(OGSPlotViews "1.0"
  SERVER_MANAGER_XML     OGSPlotViews.xml
  SERVER_MANAGER_SOURCES vtkSMOGSVerticalProfilePlotProxy.cxx  vtkOGSVerticalProfilePlot.cxx vtkSMOGSHovmoellerPlotProxy.cxx vtkOGSHovmoellerPlot.cxx vtkSMOGSSpaghettiPlotProxy.cxx vtkOGSSpaghettiPlot.cxx vtkSMOGSMapPlotProxy.cxx vtkOGSMapPlot.cxx
  REQUIRED_ON_SERVER
  GUI_INTERFACES         ${IFACES} ${VIEW_IFACE}
  GUI_SOURCES
    ${IFACE_SRCS} ${VIEW_IFACE_SRCS} ${UI_SRCS}
    ${MOC_SRCS}
    OGSVerticalProfilePlot.cxx OGSHovmoellerPlot.cxx OGSSpaghettiPlot.cxx OGSMapPlot.cxx

)