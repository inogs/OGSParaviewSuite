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
# Arnau Miro, OGS 2019

# Paraview version
PV_VERS=$1
INSTALL_PREFIX=$PWD/paraview-$PV_VERS
OSX_SDK=macosx10.14
SUITEDIR='OGSParaviewSuite'
NPROCS=$(getconf _NPROCESSORS_ONLN)

# Obtain superbuild
git clone --recursive https://gitlab.kitware.com/paraview/paraview-superbuild.git
cd paraview-superbuild
git fetch origin
git checkout "v$PV_VERS"
git submodule update

# Clone the OGS ParaView suite
git clone https://github.com/inogs/OGSParaviewSuite.git
cd $SUITEDIR
git submodule update --init --recursive
cd bit.sea
git pull origin master
cd ../..

# FIX: versions
cp $SUITEDIR/superbuild/versions/versions_new.cmake versions.cmake
cp $SUITEDIR/superbuild/versions/versions_superbuild_py27.cmake superbuild/versions.cmake
# FIX: patches
cp $SUITEDIR/superbuild/projects_apple/matplotlib.cmake superbuild/projects/apple/matplotlib.cmake

# Build
cd ../ && mkdir -p paraview-build && cd paraview-build
cmake ../paraview-superbuild/ \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
   -DBUILD_TESTING=OFF \
   -Dparaview_SOURCE_SELECTION=$PV_VERS \
   -Dqt5_SOURCE_SELECTION=5.12 \
   -DCMAKE_BUILD_TYPE_paraview=Release \
   -DCMAKE_OSX_SDK=$OSX_SDK \
   -DENABLE_mpi=ON \
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
   -DUSE_SYSTEM_python=OFF \
   -DENABLE_numpy=ON \
   -DENABLE_scipy=ON \
   -DENABLE_matplotlib=ON \
   -DENABLE_vtkm=OFF \
   -DENABLE_netcdf=OFF \
   -DENABLE_vrpn=ON \
   -DENABLE_vortexfinder2=OFF \
   -DENABLE_paraview=ON \
   -DENABLE_paraviewsdk=OFF \
   -Dparaview_PLUGINS_EXTERNAL="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSHovmoeller;OGSVortexIdentification;OGSPlotViews;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStatistics;OGSTimeStatistics;OGSVariableAggregator;OGSUtils;OGSWriter;OGSCompareVariables" \
   -Dparaview_PLUGINS_AUTOLOAD="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSHovmoeller;OGSVortexIdentification;OGSPlotViews;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStatistics;OGSTimeStatistics;OGSVariableAggregator;OGSUtils;OGSWriter;OGSCompareVariables" \
   -Dparaview_PLUGIN_OGSAnnotateDateTime_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSAnnotateDateTime \
   -Dparaview_PLUGIN_OGSDepthProfile_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSDepthProfile \
   -Dparaview_PLUGIN_OGSDerivatives_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSDerivatives \
   -Dparaview_PLUGIN_OGSHovmoeller_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSHovmoeller \
   -Dparaview_PLUGIN_OGSPlotViews_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSPlotViews \
   -Dparaview_PLUGIN_OGSVortexIdentification_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSVortexIndentification \
   -Dparaview_PLUGIN_OGSReader_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSReader \
   -Dparaview_PLUGIN_OGSSelectTools_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSelectTools \
   -Dparaview_PLUGIN_OGSSpaghetti_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSpaghetti \
   -Dparaview_PLUGIN_OGSSpatialStatistics_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSpatialStatistics \
   -Dparaview_PLUGIN_OGSTimeStatistics_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSTimeStatistics \
   -Dparaview_PLUGIN_OGSVariableAggregator_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSVariableAggregator \
   -Dparaview_PLUGIN_OGSUtils_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSUtils \
   -Dparaview_PLUGIN_OGSWriter_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSWriter \
   -Dparaview_PLUGIN_OGSCompareVariables_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSCompareVariables \
   -DCMAKE_C_COMPILER=gcc \
   -DCMAKE_CXX_COMPILER=g++ \
   -DCMAKE_Fortran_COMPILER=gfortran

# Prompt user to check configuration
printf "Deploying ParaView $PV_VERS in $INSTALL_PREFIX.\n"
read -s -p "Please check install configuration and press [enter]..."

# Make
make -j $NPROCS
# Do some fixing...
ln -s $PWD/install/lib/python2.7/site-packages/matplotlib-2.2.3-py2.7-macosx-10.14-intel.egg/matplotlib/ $PWD/install/lib/python2.7/site-packages
# Install
make -j $NPROCS install

# FIX: matplotlib depends on backports and kiwisolver
# backports is already compiled during installation
cp -r $PWD/install/lib/python2.7/site-packages/backports.functools_lru_cache-1.5-py2.7.egg/backports $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/

export PATH=$PATH:$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries
export PYTHONPATH=$PYTHONPATH:$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python

# Install netCDF4, configparser, cython
./install/bin/pip install kiwisolver netcdf4 configparser cython
cp -r install/lib/python2.7/site-packages/kiwisolver $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/cftime $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/netCDF4 $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/configparser.* $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/backports $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/backports.functools_lru_cache-1.6.1-py2.7.egg/backports/* $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/backports/
cp -r install/lib/python2.7/site-packages/cython.* $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/lib-dynload $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries/
cp -r install/lib/python2.7/lib-tk $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries/


# Install geos
GEOS_VERS=3.7.0
GEOS_DIR=geos-$GEOS_VERS
GEOS_TAR=$GEOS_DIR.tar.bz2

wget http://download.osgeo.org/geos/$GEOS_TAR
tar xf $GEOS_TAR && cd $GEOS_DIR

./configure --prefix=$INSTALL_PREFIX
make -j $NPROCS && make -j $NPROCS install
./configure --prefix=$PWD/../install
make -j $NPROCS && make -j $NPROCS install
cd ../ && rm -rf $GEOS_TAR

# Install proj libraries
PROJ_VERS=5.2.0
PROJ_DATV=1.8
PROJ_DIR=proj-$PROJ_VERS
PROJ_TAR=$PROJ_DIR.tar.gz
PROJ_DATG=proj-datumgrid-$PROJ_DATV.zip

wget https://download.osgeo.org/proj/$PROJ_TAR
tar xvf $PROJ_TAR && cd $PROJ_DIR

wget https://download.osgeo.org/proj/$PROJ_DATG
unzip $PROJ_DATG -d data/

./configure --prefix=$INSTALL_PREFIX
make -j $NPROCS && make -j $NPROCS install
./configure --prefix=$PWD/../install
make -j $NPROCS && make -j $NPROCS install
cd ../ && rm -rf $PROJ_TAR

# Install cartopy and pyepsg
#https://github.com/SciTools/cartopy/archive/v0.17.0.tar.gz
./install/bin/pip install cartopy pyepsg
cp -r install/lib/python2.7/site-packages/cartopy $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/shape* $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/pyepsg.* $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/requests $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/urllib3 $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/chardet $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/certifi $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r install/lib/python2.7/site-packages/idna $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/

# Copy ffmpeg
cp install/bin/ffmpeg $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
cd ..
printf "Install done successfully!\n"

# Deploy img2video
printf "Deploying img2video... "
cp paraview-superbuild/$SUITEDIR/bin/img2video $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
printf "OK\n"

# Deploy bit.sea inside the ParaView installation
printf "Deploying bit.sea... "
cp -r paraview-superbuild/$SUITEDIR/bit.sea/commons $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r paraview-superbuild/$SUITEDIR/bit.sea/basins  $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
printf "OK\n"

# Deploy OGStools inside the installation
printf "Deploying OGS tools... "
# libOGS
cd paraview-superbuild/$SUITEDIR/OGSPlugins/_utils
g++ -shared -o $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries/libOGS.dylib -fPIC -std=c++11 OGS.cpp -DOGS_NO_NETCDF
cd ../../../../
# Python utils
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGSmesh.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGSlonlat2m.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGS2Paraview.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/default.ini $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
printf "OK\n"