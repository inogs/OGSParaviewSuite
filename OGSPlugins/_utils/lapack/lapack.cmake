##############################################
###       LAPACK CMAKE FOR PARAVIEW        ###
##############################################

set( LAPACK_INC_FOLDER "${CMAKE_CURRENT_LIST_DIR}/include")

# CMAKE will change according to UNIX or MACOSX
# system
if(NOT APPLE)
	# Hence we are in UNIX
	set( LAPACK_LIB_FOLDER "${CMAKE_CURRENT_LIST_DIR}/unix")
else()
	# Hence we are in MACOSX
	set( LAPACK_LIB_FOLDER "${CMAKE_CURRENT_LIST_DIR}/macos")
endif()

set ( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -I${LAPACK_INC_FOLDER} -DLAPACK -lm -lgfortran") 
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${LAPACK_INC_FOLDER} -DLAPACK -lm -lgfortran")

SET(SRC_LAPACK
	${LAPACK_LIB_FOLDER}/liblapacke.a
	${LAPACK_LIB_FOLDER}/liblapack.a
	${LAPACK_LIB_FOLDER}/libblas.a
)