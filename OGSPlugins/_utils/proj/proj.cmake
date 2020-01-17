##############################################
###         PROJ CMAKE FOR PARAVIEW        ###
##############################################

set( PROJ_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}/include")

# CMAKE will change according to UNIX or MACOSX
# system
if(NOT APPLE)
    # Hence we are in UNIX
    set( PROJ_LIB_FOLDER "${CMAKE_CURRENT_LIST_DIR}/unix")
else()
    # Hence we are in MACOSX
    set( PROJ_LIB_FOLDER "${CMAKE_CURRENT_LIST_DIR}/macos")
endif()

set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${PROJ_INC_FOLDER}") 
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PROJ_INC_FOLDER}")

set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-whole-archive ${PROJ_LIB_FOLDER}/libproj.a -Wl,-no-whole-archive") 

set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-whole-archive ${PROJ_LIB_FOLDER}/libproj.a -Wl,-no-whole-archive") 

set ( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -Wl,-whole-archive ${PROJ_LIB_FOLDER}/libproj.a -Wl,-no-whole-archive") 
