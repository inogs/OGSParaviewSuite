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
include(${PROJ})  # PROJ detection
include(${CNPY})  # CNPY detection

#--------------------------------------------------
# ParaView plugin
#--------------------------------------------------
add_paraview_plugin(OGSWriter "1.0"
  SERVER_MANAGER_XML     OGSWriter.xml
  SERVER_MANAGER_SOURCES ${OGS} ${VTKFIELDS} ${VTKOPS} ${CNPY} vtkOGSNPZWriter.cxx vtkOGSNPZTableWriter.cxx vtkOGSFieldWriter.cxx vtkOGSNetCDFWriter.cxx
  REQUIRED_ON_SERVER
)
