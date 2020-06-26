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
include(${UTILS}) # UTILS detection

#--------------------------------------------------
# ParaView plugin
#--------------------------------------------------
add_paraview_plugin(OGSTimeStatistics "1.0"
  SERVER_MANAGER_XML     OGSTimeStatistics.xml
  SERVER_MANAGER_SOURCES ${VTKFIELDS} ${OGS} vtkOGSTimeCommons.cxx vtkOGSTimeAverage.cxx vtkOGSTimeAggregator.cxx vtkOGSClimatology.cxx
  REQUIRED_ON_SERVER
)