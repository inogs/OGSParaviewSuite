##############################################
###       LAPACK CMAKE FOR PARAVIEW        ###
##############################################

set( LAPACK_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}/include")
set( LAPACK_LIB_FOLDER "${CMAKE_CURRENT_LIST_DIR}/lib")

# Apple cannot link with gfortran so deactivate the use of LAPACK
if(NOT APPLE)
	set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${LAPACK_INC_FOLDER} -DUSE_LAPACK") 
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${LAPACK_INC_FOLDER} -DUSE_LAPACK")
	
	# Compiler flags for linux
	set( LAPACK_LINK_FLAGS "-Wl,-whole-archive ${LAPACK_LIB_FOLDER}/liblapacke.a ${LAPACK_LIB_FOLDER}/liblapack.a -Wl,-no-whole-archive ${LAPACK_LIB_FOLDER}/libblas.a")
	
	# Linking flags according to the compiler
	if(NOT "${COMPILER}" STREQUAL "INTEL")
		set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINK_FLAGS} -lm -lgfortran") 
		set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINK_FLAGS} -lm -lgfortran") 
		set( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${LAPACK_LINK_FLAGS} -lm -lgfortran")
	else()
		set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINK_FLAGS} -lm -lifcore") 
		set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINK_FLAGS} -lm -lifcore") 
		set( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${LAPACK_LINK_FLAGS} -lm -lifcore")
	endif()
endif()
