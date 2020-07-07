#!/bin/bash
#
# SCRIPT to deploy the GEOS libraries compiled
# as static libraries
#
# Arnau Miro, OGS (2020)

VERS=${1}
INSTALL_PREFIX=${2}
CCOMPILER=${3}
CFLAGS=${4}
CXXCOMPILER=${5}
CXXFLAGS=${6}

GEOS_DIR="geos-${VERS}"
GEOS_TAR="${GEOS_DIR}.tar.bz2"

# Check if the Lapack libraries have been deployed
# if not, compile
if test -f "${INSTALL_PREFIX}/lib/libgeos.a"; then
    echo "GEOS already deployed!"
    echo "Skipping build..."
else
    echo "GEOS not deployed!"
    echo "Starting build..."

    echo "Version ${VERS}"
    echo "C compiler '${CCOMPILER}' with flags '${CFLAGS}'"
    echo "CXX compiler '${CXXCOMPILER}' with flags '${CXXFLAGS}'"
    echo "Install path ${INSTALL_PREFIX}"

	wget -O ${GEOS_TAR} "http://download.osgeo.org/geos/${GEOS_TAR}"
	tar xf $GEOS_TAR && cd $GEOS_DIR
	cd ../

	# Configure 
	mkdir build
	cd build
	cmake ../$GEOS_DIR \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
		-DCMAKE_INSTALL_LIBDIR=$INSTALL_PREFIX/lib \
		-DGEOS_ENABLE_TESTS=OFF \
		-DGEOS_BUILD_SHARED=ON -DGEOS_BUILD_STATIC=ON \
		-DCMAKE_C_COMPILER="${CCOMPILER}" -DCMAKE_C_FLAGS="${CFLAGS}" \
		-DCMAKE_CXX_COMPILER="${CXXCOMPILER}" -DCMAKE_CXX_FLAGS="${CXXFLAGS}" 

	# Build
	make -j $(getconf _NPROCESSORS_ONLN)
	make install

	# Cleanup
	cd ..
	rm -rf build $GEOS_TAR $GEOS_DIR
fi