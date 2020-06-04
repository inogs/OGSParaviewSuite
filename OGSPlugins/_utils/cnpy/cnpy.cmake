##############################################
###         CNPY CMAKE FOR PARAVIEW        ###
##############################################

# List the folder where all the headers are
set( CNPY_FOLDER "${CMAKE_CURRENT_LIST_DIR}")
set( CNPY_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}")

# Add the folders in the include path
set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${CNPY_INC_FOLDER}") 
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${CNPY_INC_FOLDER}")

# Add lz in the linking flags
set ( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} -lz") 
set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lz") 
set ( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -lz")

# Now define a set of variables containing various code
# that might be linked when compiling
set( CNPY
	${CNPY_FOLDER}/cnpy.cpp
	)
