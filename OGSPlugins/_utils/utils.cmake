##############################################
###        UTILS CMAKE FOR PARAVIEW        ###
##############################################

# List the folder where all the headers are
set( UTILS_FOLDER     "${CMAKE_CURRENT_LIST_DIR}")
set( UTILS_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}")

# Add the folders in the include path
set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${UTILS_INC_FOLDER}") 
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${UTILS_INC_FOLDER}")

# Now define a set of variables containing various code
# that might be linked when compiling

# VTK field creation
set( VTKFIELDS
	${UTILS_FOLDER}/vtkFields.cpp
	)

# Operations with field arrays
set( FIELDOPS
	${UTILS_FOLDER}/fieldOperations.cpp
	)

# Operations with VTK fields
set( VTKOPS
	${UTILS_FOLDER}/vtkOperations.cpp
	)

# Field and VTK operations combined
set( OPS
	${UTILS_FOLDER}/fieldOperations.cpp
	${UTILS_FOLDER}/vtkOperations.cpp
	)

# NetCDF Input/Output
set( NETCDFIO 
	${UTILS_FOLDER}/netcdfio.cpp
	)

# OGS module
set( OGS 
	${UTILS_FOLDER}/netcdfio.cpp
	${UTILS_FOLDER}/OGS.cpp
	)
