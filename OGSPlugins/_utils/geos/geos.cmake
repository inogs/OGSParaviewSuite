##############################################
###         GEOS CMAKE FOR PARAVIEW        ###
##############################################

set( GEOS_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}/include")
set( GEOS_LIB_FOLDER "${CMAKE_CURRENT_LIST_DIR}/lib")

set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${GEOS_INC_FOLDER}") 
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${GEOS_INC_FOLDER}")

if(NOT APPLE)
	set( GEOS_LINK_FLAGS "-Wl,-whole-archive ${GEOS_LIB_FOLDER}/libgeos.a -Wl,-no-whole-archive")
else()
	set( GEOS_LINK_FLAGS "-Wl,${GEOS_LIB_FOLDER}/libgeos.a")
endif()

set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${GEOS_LINK_FLAGS}") 
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${GEOS_LINK_FLAGS}") 
set( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${GEOS_LINK_FLAGS}") 
