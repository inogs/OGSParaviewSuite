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
NPROCS=2

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
cd ../ && mkdir paraview-build && cd paraview-build
cmake ../paraview-superbuild/ \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
   -DBUILD_TESTING=OFF \
   -Dparaview_SOURCE_SELECTION=$PV_VERS \
   -Dqt5_SOURCE_SELECTION=5.12 \
   -DCMAKE_BUILD_TYPE_paraview=Release \
   -DCMAKE_OSX_SDK=$OSX_SDK \
   -DENABLE_mpi=OFF \
   -DENABLE_cuda=OFF \
   -DENABLE_boost=ON \
   -DENABLE_lapack=ON \
   -DENABLE_ffmpeg=ON \
   -DENABLE_bzip2=ON \
   -DENABLE_libxml2=ON \
   -DENABLE_png=ON \
   -DENABLE_qt5=ON \
   -DENABLE_tbb=ON \
   -DENABLE_embree=ON \
   -DENABLE_ospray=ON \
   -DENABLE_python=ON \
   -DENABLE_numpy=ON \
   -DENABLE_scipy=ON \
   -DENABLE_matplotlib=ON \
   -DENABLE_vtkm=ON \
   -DENABLE_paraview=ON \
   -DENABLE_netcdf=OFF \
   -DENABLE_vrpn=ON \
   -DENABLE_vortexfinder2=OFF \
   -Dparaview_PLUGINS_EXTERNAL="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSGradient;OGSHovmoeller;OGSLambda2Criterion;OGSOkuboWeiss;OGSOmegaCriterion;OGSPlotViews;OGSQCriterion;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStats;OGSSpatialStatsFromFile;OGSTimeStatistics;OGSVariableAggregator;OGSUtils" \
   -Dparaview_PLUGINS_AUTOLOAD="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSGradient;OGSHovmoeller;OGSLambda2Criterion;OGSOkuboWeiss;OGSOmegaCriterion;OGSPlotViews;OGSQCriterion;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStats;OGSSpatialStatsFromFile;OGSTimeStatistics;OGSVariableAggregator;OGSUtils" \
   -Dparaview_PLUGIN_OGSAnnotateDateTime_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSAnnotateDateTime \
   -Dparaview_PLUGIN_OGSDepthProfile_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSDepthProfile \
   -Dparaview_PLUGIN_OGSDerivatives_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSDerivatives \
   -Dparaview_PLUGIN_OGSGradient_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSGradient \
   -Dparaview_PLUGIN_OGSHovmoeller_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSHovmoeller \
   -Dparaview_PLUGIN_OGSLambda2Criterion_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSLambda2Criterion \
   -Dparaview_PLUGIN_OGSOkuboWeiss_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSOkuboWeiss \
   -Dparaview_PLUGIN_OGSOmegaCriterion_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSOmegaCriterion \
   -Dparaview_PLUGIN_OGSPlotViews_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSPlotViews \
   -Dparaview_PLUGIN_OGSQCriterion_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSQCriterion \
   -Dparaview_PLUGIN_OGSReader_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSReader \
   -Dparaview_PLUGIN_OGSSelectTools_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSelectTools \
   -Dparaview_PLUGIN_OGSSpaghetti_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSpaghetti \
   -Dparaview_PLUGIN_OGSSpatialStats_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSpatialStats \
   -Dparaview_PLUGIN_OGSSpatialStatsFromFile_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSSpatialStatsFromFile \
   -Dparaview_PLUGIN_OGSTimeStatistics_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSTimeStatistics \
   -Dparaview_PLUGIN_OGSVariableAggregator_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSVariableAggregator \
   -Dparaview_PLUGIN_OGSUtils_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSUtils \
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

# Dowload extra sources from google drive
fileid="13L4jjYfTa85cAluJdqZYa0MhU4j9pLP4"
filename="extra_src.tar.gz"
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
tar xzf extra_src.tar.gz

# FIX: matplotlib depends on backports and kiwisolver
# backports is already compiled during installation
cp -r $PWD/install/lib/python2.7/site-packages/backports.functools_lru_cache-1.5-py2.7.egg/backports/ $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
# Unpack kiwisolver wheel
unzip extra_src/kiwisolver-1.1.0.whl -d $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
# Unpack cftime wheel
unzip extra_src/cftime-1.0.3.4.whl -d $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
# Unpack netCDF4 wheel
unzip extra_src/netCDF4-1.4.2.whl -d $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
# Compile and install geos
tar xf extra_src/geos-3.7.0.tar.bz2 && cd geos-3.7.0
./configure --prefix=$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/ --libdir=$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries
make -j $NPROCS install
cd ../
# Compile and install basemap
tar xfz extra_src/basemap-1.1.0.tar.gz && cd basemap-1.1.0
export GEOS_DIR=$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/
ln -s $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/lib
ln -s $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/lib64
python2.7 setup.py build
cp build/lib.macosx-10.14-intel-2.7/_geoslib.so $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp -r build/lib.macosx-10.14-intel-2.7/mpl_toolkits/basemap $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/mpl_toolkits
cd ../
# Compile and install proj
tar xzf extra_src/proj-6.1.0.tar.gz && cd proj-6.1.0
./configure --prefix=$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/ --libdir=$INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries
unzip extra_src/proj-datumgrid-1.8.zip -d data/
make -j $NPROCS install
cd ../
# Compile and install pyproj
tar xzf extra_src/pyproj-1.9.6.tar.gz && cd pyproj-1.9.6
python2.7 setup.py build
cp -r build/lib.macosx-10.14-intel-2.7/pyproj $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cd ../
# Remove dist-info folders
rm -rf $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/*.dist-info
# Also copy ffmpeg
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
g++ -shared -o $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries/libOGS.dylib -fPIC -std=c++11 OGS.cpp netcdfio.cpp -L $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Libraries/ -l vtkNetCDF-pv5.6
cd ../../../../
# Python utils
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGSmesh.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGSlonlat2m.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/OGS2Paraview.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/Python/
cp paraview-superbuild/$SUITEDIR/OGSPlugins/_utils/python/default.ini $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
# Symbolic links
ln -s paraview-$PV_VERS/ParaView-$PV_VERS.app/Contents/Python/OGS2Paraview.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
chmod +x $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/OGS2Paraview.py
ln -s paraview-$PV_VERS/ParaView-$PV_VERS.app/Contents/Python/OGSlonlat2m.py $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/
chmod +x $INSTALL_PREFIX/ParaView-$PV_VERS.app/Contents/bin/OGSlonlat2m.py
printf "OK\n"