##############################################
###        PUGIXML CMAKE FOR PARAVIEW      ###
##############################################

# List the folder where all the headers are
set( PUGI_FOLDER "${CMAKE_CURRENT_LIST_DIR}")
set( PUGI_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}")

# Add the folders in the include path
set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${PUGI_INC_FOLDER}") 
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PUGI_INC_FOLDER}")

# Now define a set of variables containing various code
# that might be linked when compiling

set( PUGIXML
	${PUGI_FOLDER}/pugixml.cpp
	)
