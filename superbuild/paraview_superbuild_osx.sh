#!/bin/bash
#
# Automated script to configure and install redistributable 
# PARAVIEW binaries using the superbuild (for OSX).
#
# Issues:
#  - PYTHONPATH must be unset, otherwise python dependent libraries
#    may fail to compile.
#
#  - vortexfinder2 issue with CHAR_WIDTH; the solution is to change
#    CHAR_WIDTH by CHAR_SIZE, as described in:
#    https://github.com/fmtlib/fmt/commit/abbefd71666055daac9e14e78262620f9e845850
#
#  - vrpn issue with union wait; the solution is to change union wait to int, 
#    as described in: https://github.com/opensgct/sgct/issues/13
#    This has been solved using the release version 07.34
#
#  - paraviewsdk is incompatible with building qt5
#
#  - as of now, CUDA is deactivated since it does not support gcc >= 8
#
#  - netcdf must be OFF since the OGS plugins are built on top of VTK NetCDF
#
#  - qt5 (5.8, 5.9 and 5.10) will fail with the new SDK. Solution is to use qt vers 5.12
#
#  - vrpn cannot find pthread.h in OS Mojave. Solved using the newest version.
#
# Arnau Miro, OGS 2020

PV_VERS=${1}
OSX_SDK=${2}
QT5_VERS=${3}
PROJ_VERS=${4}
PROJ_DATV=${5}
GEOS_VERS=${6}
CCOMPILER=${7}
CFLAGS=${8}
CXXCOMPILER=${9}
CXXFLAGS=${10}
FCOMPILER=${11}
FFLAGS=${12}

NPROCS=$(getconf _NPROCESSORS_ONLN)
INSTALL_PREFIX="${PWD}/../paraview-${PV_VERS}"
APP_DIR="${INSTALL_PREFIX}/ParaView-${PV_VERS}.app/Contents/"
SUITEDIR="${PWD}"
SUPERBUILD_DIR="${PWD}/../paraview-superbuild"
BUILD_DIR="${PWD}/../paraview-build"

VERSIONSDIR="${SUITEDIR}/superbuild/versions"
MAPPLOTLIBDIR="${SUITEDIR}/superbuild/matplotlib_py27"
PLUGINDIR="${SUITEDIR}/OGSPlugins/"
BINDIR="${SUITEDIR}/bin"
BITSEADIR="${SUITEDIR}/bit.sea"

# Load modules - provide a basic building environment
if [ $(command -v module) ]; then
   module purge
   module load gcc openmpi cmake
fi

# Obtain superbuild
git clone --recursive https://gitlab.kitware.com/paraview/paraview-superbuild.git $SUPERBUILD_DIR
cd $SUPERBUILD_DIR
git fetch origin
git checkout "v$PV_VERS"
git submodule update
# FIX: versions
cp $VERSIONSDIR/versions_new.cmake versions.cmake
cp $VERSIONSDIR/versions_superbuild_py27.cmake superbuild/versions.cmake
# FIX: python
cp superbuild/projects/zlib.cmake superbuild/projects/apple/zlib.cmake
cp superbuild/projects/bzip2.cmake superbuild/projects/apple/bzip2.cmake
cp superbuild/projects/patches/zlib* superbuild/projects/apple/patches
cp superbuild/projects/patches/bzip2* superbuild/projects/apple/patches
cp $SUITEDIR/superbuild/projects_apple/python.cmake superbuild/projects/apple/python.cmake
cp $SUITEDIR/superbuild/projects_apple/python-ssl.patch superbuild/projects/apple/patches
# FIX: matplotlib
cp $MAPPLOTLIBDIR/matplotlib.cmake superbuild/projects
cp $SUITEDIR/superbuild/projects_apple/matplotlib.cmake superbuild/projects/apple/matplotlib.cmake
cp $MAPPLOTLIBDIR/matplotlib-kiwisolver.patch superbuild/projects/patches
# Return to main folder
cd $SUITEDIR

# Build
mkdir -p $BUILD_DIR && cd $BUILD_DIR
cmake $SUPERBUILD_DIR \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
   -DBUILD_TESTING=OFF \
   -Dparaview_SOURCE_SELECTION=$PV_VERS \
   -Dqt5_SOURCE_SELECTION=$QT5_VERS \
   -DCMAKE_BUILD_TYPE_paraview=Release \
   -DCMAKE_OSX_SDK=$OSX_SDK \
   -DENABLE_mpi=OFF \
   -DUSE_SYSTEM_mpi=OFF \
   -DENABLE_cuda=OFF \
   -DENABLE_boost=ON \
   -DENABLE_lapack=ON \
   -DENABLE_ffmpeg=ON \
   -DENABLE_bzip2=ON \
   -DENABLE_libxml2=ON \
   -DENABLE_png=ON \
   -DENABLE_qt5=ON \
   -DENABLE_mesa=ON \
   -DENABLE_osmesa=OFF \
   -DENABLE_tbb=ON \
   -DENABLE_embree=ON \
   -DENABLE_ospray=ON \
   -DENABLE_python=ON \
   -Dpython_USE_UNICODE=UCS2 \
   -DUSE_SYSTEM_python=OFF \
   -DENABLE_numpy=ON \
   -DENABLE_scipy=ON \
   -DENABLE_matplotlib=ON \
   -DENABLE_vtkm=OFF \
   -DENABLE_netcdf=OFF \
   -DENABLE_vrpn=ON \
   -DENABLE_nvidiaindex=OFF \
   -DENABLE_vortexfinder2=OFF \
   -DENABLE_paraview=ON \
   -DENABLE_paraviewsdk=OFF \
   -DCMAKE_C_COMPILER=${CCOMPILER} \
   -DCMAKE_CXX_COMPILER=${CXXCOMPILER} \
   -DCMAKE_Fortran_COMPILER=${FCOMPILER} \
   -Dparaview_PLUGINS_EXTERNAL="MapPlotterView;OGSAnnotateDateTime;OGSCompareVariables;OGSDensity;OGSDepthProfile;OGSDerivatives;OGSHovmoeller;OGSMixingLayerDepth;OGSPlotViews;OGSReader;OGSRossbyRadius;OGSSelectTools;OGSSpaghetti;OGSSpatialStatistics;OGSTimeStatistics;OGSUtils;OGSVariableAggregator;OGSVortexDetection;OGSVortexIdentification;OGSWriter" \
   -Dparaview_PLUGINS_AUTOLOAD="MapPlotterView;OGSAnnotateDateTime;OGSCompareVariables;OGSDensity;OGSDepthProfile;OGSDerivatives;OGSHovmoeller;OGSMixingLayerDepth;OGSPlotViews;OGSReader;OGSRossbyRadius;OGSSelectTools;OGSSpaghetti;OGSSpatialStatistics;OGSTimeStatistics;OGSUtils;OGSVariableAggregator;OGSVortexDetection;OGSVortexIdentification;OGSWriter" \
   -Dparaview_PLUGIN_MapPlotterView_PATH=$PLUGINDIR/MapPlotterView \
   -Dparaview_PLUGIN_OGSAnnotateDateTime_PATH=$PLUGINDIR/OGSAnnotateDateTime \
   -Dparaview_PLUGIN_OGSCompareVariables_PATH=$PLUGINDIR/OGSCompareVariables \
   -Dparaview_PLUGIN_OGSDensity_PATH=$PLUGINDIR/OGSDensity \
   -Dparaview_PLUGIN_OGSDepthProfile_PATH=$PLUGINDIR/OGSDepthProfile \
   -Dparaview_PLUGIN_OGSDerivatives_PATH=$PLUGINDIR/OGSDerivatives \
   -Dparaview_PLUGIN_OGSHovmoeller_PATH=$PLUGINDIR/OGSHovmoeller \
   -Dparaview_PLUGIN_OGSMixingLayerDepth_PATH=$PLUGINDIR/OGSMixingLayerDepth \
   -Dparaview_PLUGIN_OGSPlotViews_PATH=$PLUGINDIR/OGSPlotViews \
   -Dparaview_PLUGIN_OGSReader_PATH=$PLUGINDIR/OGSReader \
   -Dparaview_PLUGIN_OGSRossbyRadius_PATH=$PLUGINDIR/OGSRossbyRadius \
   -Dparaview_PLUGIN_OGSSelectTools_PATH=$PLUGINDIR/OGSSelectTools \
   -Dparaview_PLUGIN_OGSSpaghetti_PATH=$PLUGINDIR/OGSSpaghetti \
   -Dparaview_PLUGIN_OGSSpatialStatistics_PATH=$PLUGINDIR/OGSSpatialStatistics \
   -Dparaview_PLUGIN_OGSTimeStatistics_PATH=$PLUGINDIR/OGSTimeStatistics \
   -Dparaview_PLUGIN_OGSUtils_PATH=$PLUGINDIR/OGSUtils \
   -Dparaview_PLUGIN_OGSVariableAggregator_PATH=$PLUGINDIR/OGSVariableAggregator \
   -Dparaview_PLUGIN_OGSVortexDetection_PATH=$PLUGINDIR/OGSVortexDetection \
   -Dparaview_PLUGIN_OGSVortexIdentification_PATH=$PLUGINDIR/OGSVortexIdentification \
   -Dparaview_PLUGIN_OGSWriter_PATH=$PLUGINDIR/OGSWriter 

# Prompt user to check configuration
printf "Deploying ParaView ${PV_VERS} in ${INSTALL_PREFIX}.\n"
#read -s -p "Please check install configuration and press [enter]..."

# Make
make -j $NPROCS
# Do some fixing...
ln -s $BUILD_DIR/install/lib/python2.7/site-packages/matplotlib-*/matplotlib/ $BUILD_DIR/install/lib/python2.7/site-packages
# Install
make -j $NPROCS install
# FIX: matplotlib depends on backports and kiwisolver
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/backports.*/backports $APP_DIR/Python/
cp $BUILD_DIR/install/lib/python2.7/site-packages/kiwisolver-*/kiwisolver.so $APP_DIR/Python/

# Load environment
export PATH=$PATH:$APP_DIR/bin
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:$APP_DIR/Libraries:$APP_DIR/Libraries/paraview-5.6
#export C_INCLUDE_PATH=$C_INCLUDE_PATH:$APP_DIR/include
export PYTHONPATH=$PYTHONPATH:$APP_DIR/Python

# Deploy GEOS library
bash $PLUGINDIR/_utils/geos/install_geos.sh "${GEOS_VERS}" "${BUILD_DIR}/install" "${CCOMPILER}" "${CFLAGS}" "${CXXCOMPILER}" "${CXXFLAGS}"
cp -r $BUILD_DIR/install/lib/libgeos* $APP_DIR/Libraries/

# Deploy PROJ library
bash $PLUGINDIR/_utils/proj/install_proj.sh "${PROJ_VERS}" "${PROJ_DATV}" "SHARED" "${BUILD_DIR}/install" "${CCOMPILER}" "${CFLAGS}" "${CXXCOMPILER}" "${CXXFLAGS}"
cp -r $BUILD_DIR/install/lib/libproj* $APP_DIR/Libraries/

# Install netCDF4, configparser, cython
$BUILD_DIR/install/bin/pip install --upgrade pip
$BUILD_DIR/install/bin/pip install kiwisolver requests netcdf4 configparser cython cartopy pyepsg
#cp -r $BUILD_DIR/install/include/python* $INSTALL_PREFIX/include
#cp -r $BUILD_DIR/install/include/paraview-* $INSTALL_PREFIX/include
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/kiwisolver $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/cftime $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/netCDF4 $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/setuptools $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/setuptools.pth $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/configparser.* $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/backports $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/backports.functools_lru_cache-1.6.1-py2.7.egg/backports/* $APP_DIR/Python/backports/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/cython.* $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/Cython $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/lib-dynload $APP_DIR/Libraries/
cp -r $BUILD_DIR/install/lib/python2.7/lib-tk $APP_DIR/Libraries/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/cartopy $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/shape* $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/pyepsg.* $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/requests $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/urllib3 $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/chardet $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/certifi $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/idna $APP_DIR/Python/
cp -r $BUILD_DIR/install/lib/python2.7/site-packages/numpy-*/numpy/core/include $APP_DIR/Python/numpy/core

# Copy ffmpeg
cp $BUILD_DIR/install/bin/ffmpeg $APP_DIR/bin
cd ..
printf "Install done successfully!\n"

# Deploy img2video
printf "Deploying img2video... "
cp $BINDIR/img2video $APP_DIR/bin/
printf "OK\n"

# Deploy bit.sea inside the ParaView installation
printf "Deploying bit.sea... "
cp -r $BITSEADIR/commons    $APP_DIR/Python
cp -r $BITSEADIR/basins     $APP_DIR/Python
cp -r $BITSEADIR/MapPlotter $APP_DIR/Python
printf "OK\n"

# Deploy OGSMesh and OGS2Paraview inside the installation
printf "Deploying OGSMesh and OGS2Paraview... "
mv $SUITEDIR/libOGS.dylib $APP_DIR/Libraries
cp $SUITEDIR/superbuild/env-linux.sh $APP_DIR/env.sh
cp $PLUGINDIR/_utils/python/OGSmesh.py $APP_DIR/Python
cp $PLUGINDIR/_utils/python/OGSlonlat2m.py $APP_DIR/bin
cp $PLUGINDIR/_utils/python/OGS2Paraview.py $APP_DIR/bin
cp $PLUGINDIR/_utils/python/default.ini $APP_DIR/bin
printf "OK\n"

# Clean-up
cd $SUITEDIR
#rm -rf $BUILD_DIR $SUPERBUILD_DIR
tar cvzf "${INSTALL_PREFIX}/ParaView-${PV_VERS}.tar.gz" ${INSTALL_PREFIX}/ParaView-${PV_VERS}.app
