#!/bin/bash
#
# SCRIPT to deploy the PROJ libraries compiled
# as static libraries
#
# Arnau Miro, OGS (2020)

VERS=${1}
VERS_DATV=${2}
INSTALL_PREFIX=${3}
CCOMPILER=${4}
CFLAGS=${5}
CXXCOMPILER=${6}
CXXFLAGS=${7}

PROJ_DIR="proj-${VERS}"
PROJ_TAR="${PROJ_DIR}.tar.gz"
PROJ_DATG="proj-datumgrid-${VERS_DATV}.zip"

# Check if the Lapack libraries have been deployed
# if not, compile
if test -f "${INSTALL_PREFIX}/lib/libproj.a"; then
    echo "PROJ already deployed!"
    echo "Skipping build..."
else
    echo "PROJ not deployed!"
    echo "Starting build..."

    echo "Version ${VERS}"
    echo "C compiler '${CCOMPILER}' with flags '${CFLAGS}'"
    echo "CXX compiler '${CXXCOMPILER}' with flags '${CXXFLAGS}'"
    echo "Install path ${INSTALL_PREFIX}"

	wget -O ${PROJ_TAR} "https://download.osgeo.org/proj/${PROJ_TAR}"
	tar xvf $PROJ_TAR && cd $PROJ_DIR
	wget -O ${PROJ_DATG} "https://download.osgeo.org/proj/$PROJ_DATG"
	unzip $PROJ_DATG -d data/
	cd ../

	# Configure 
	mkdir build
	cd build
	cmake ../$PROJ_DIR \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
		-DPROJ_TESTS=OFF \
		-DBUILD_LIBPROJ_SHARED=OFF -DBUILD_LIBPROJ_STATIC=ON \
		-DCMAKE_C_COMPILER="${CCOMPILER}" -DCMAKE_C_FLAGS="${CFLAGS}" \
		-DCMAKE_CXX_COMPILER="${CXXCOMPILER}" -DCMAKE_CXX_FLAGS="${CXXFLAGS}" 

	# Build
	make -j $(getconf _NPROCESSORS_ONLN)
	make install

	# Cleanup
	cd ..
	rm -rf build $PROJ_TAR $PROJ_DIR
fi