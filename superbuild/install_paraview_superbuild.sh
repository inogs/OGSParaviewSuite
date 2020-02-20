#!/bin/bash
#
# Automated script to configure and install redistributable 
# PARAVIEW binaries using the superbuild.
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
#  - mesa version interferes with QT5 > 5.10
#
#  - under these circumstances VTKm fails.
#
# Arnau Miro, OGS 2019

# Paraview version
PV_VERS=$1
INSTALL_PREFIX=$PWD/paraview-$PV_VERS
SUITEDIR='OGSParaviewSuite'
NPROCS=$(getconf _NPROCESSORS_ONLN)

# Load modules
# provide a basic building environment
module purge
module load gcc openmpi cmake

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

# Build
cd ../ && mkdir -p paraview-build && cd paraview-build
cmake ../paraview-superbuild/ \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
   -DBUILD_TESTING=OFF \
   -Dparaview_SOURCE_SELECTION=$PV_VERS \
   -Dqt5_SOURCE_SELECTION=5.10 \
   -DCMAKE_BUILD_TYPE_paraview=Release \
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
   -Dparaview_PLUGIN_OGSCompareVariables_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSCompareVariables 

# Prompt user to check configuration
printf "Deploying ParaView $PV_VERS in $INSTALL_PREFIX.\n"
read -s -p "Please check install configuration and press [enter]..."

# Make
make -j $NPROCS
# Do some fixing...
ln -s $PWD/install/lib/python2.7/site-packages/matplotlib-2.2.3-py2.7-linux-x86_64.egg/matplotlib/ $PWD/install/lib/python2.7/site-packages
# Install
make -j $NPROCS install
# FIX: matplotlib depends on backports and kiwisolver
cp -r install/lib/python2.7/site-packages/backports.functools_lru_cache-1.5-py2.7.egg/backports $INSTALL_PREFIX/lib/python2.7/site-packages/
cp install/lib/python2.7/site-packages/kiwisolver-1.1.0-py2.7-linux-x86_64.egg/kiwisolver.so $INSTALL_PREFIX/lib/python2.7/site-packages/

# Load environment
export PATH=$PATH:$INSTALL_PREFIX/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PREFIX/lib:$INSTALL_PREFIX/lib/paraview-5.6
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$INSTALL_PREFIX/include
export PYTHONPATH=$PYTHONPATH:$INSTALL_PREFIX/lib/python2.7/site-packages

# Install netCDF4, configparser, cython
./install/bin/pip install netcdf4 configparser cython
cp -r install/lib/python2.7/site-packages/cftime $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/netCDF4 $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/configparser.* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/backports $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/backports.functools_lru_cache-1.6.1-py2.7.egg/backports/* $INSTALL_PREFIX/lib/python2.7/site-packages/backports/
cp -r install/lib/python2.7/site-packages/cython.* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/lib-dynload $INSTALL_PREFIX/lib/python2.7/
cp -r install/lib/python2.7/lib-tk $INSTALL_PREFIX/lib/python2.7/

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
cp -r install/lib/python2.7/site-packages/cartopy $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/shape* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/pyepsg.* $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/requests $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/urllib3 $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/chardet $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/certifi $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/idna $INSTALL_PREFIX/lib/python2.7/site-packages/

# Copy ffmpeg
cp install/bin/ffmpeg $INSTALL_PREFIX/bin
cd ..
printf "Install done successfully!\n"

# Deploy img2video
printf "Deploying img2video... "
cp paraview-superbuild/$SUITEDIR/bin/img2video paraview-$PV_VERS/bin
printf "OK\n"

# Deploy bit.sea inside the ParaView installation
printf "Deploying bit.sea... "
cp -r paraview-superbuild/$SUITEDIR/bit.sea/commons paraview-$PV_VERS/lib/python*/site-packages
cp -r paraview-superbuild/$SUITEDIR/bit.sea/basins  paraview-$PV_VERS/lib/python*/site-packages
printf "OK\n"

# Compile OGSMesh
printf "Compiling OGSMesh... "
cd paraview-superbuild/$SUITEDIR/OGSPlugins/_utils
g++ -shared -Wl,-soname,libOGS -o libOGS.so -fPIC -fopenmp OGS.cpp netcdfio.cpp -I../../../../paraview-build/install/include/paraview-5.6
cd ../../../../
printf "OK\n"

# Deploy OGSMesh and OGS2Paraview inside the installation
printf "Deploying OGSMesh and OGS2Paraview... "
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/libOGS.so $INSTALL_PREFIX/lib
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGSmesh.py $INSTALL_PREFIX/lib/python*/site-packages
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGSlonlat2m.py $INSTALL_PREFIX/bin
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGS2Paraview.py $INSTALL_PREFIX/bin
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/default.ini $INSTALL_PREFIX/bin

cd ..
printf "OK\n"
