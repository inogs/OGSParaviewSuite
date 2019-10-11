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
NPROCS=24

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
   -DENABLE_qt5=OFF \
   -DENABLE_mesa=OFF \
   -DENABLE_osmesa=ON \
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
   -Dparaview_PLUGINS_EXTERNAL="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSGradient;OGSHovmoeller;OGSLambda2Criterion;OGSOkuboWeiss;OGSOmegaCriterion;OGSPlotViews;OGSQCriterion;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStats;OGSSpatialStatsFromFile;OGSTimeStatistics;OGSVariableAggregator;OGSUtils;OGSWriter;OGSCompareVariables" \
   -Dparaview_PLUGINS_AUTOLOAD="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSGradient;OGSHovmoeller;OGSLambda2Criterion;OGSOkuboWeiss;OGSOmegaCriterion;OGSPlotViews;OGSQCriterion;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStats;OGSSpatialStatsFromFile;OGSTimeStatistics;OGSVariableAggregator;OGSUtils;OGSWriter;OGSCompareVariables" \
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

# Dowload extra sources from google drive
fileid="13L4jjYfTa85cAluJdqZYa0MhU4j9pLP4"
filename="extra_src.tar.gz"
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
tar xzf extra_src.tar.gz

# Load environment
export PATH=$PATH:$INSTALL_PREFIX/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_PREFIX/lib:$INSTALL_PREFIX/lib/paraview-5.6
export PYTHONPATH=$PYTHONPATH:$INSTALL_PREFIX/lib/python2.7/site-packages

# Install netCDF4
./install/bin/pip install netcdf4
cp -r install/lib/python2.7/site-packages/cftime $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r install/lib/python2.7/site-packages/netCDF4 $INSTALL_PREFIX/lib/python2.7/site-packages/

# Compile and install geos
tar xf extra_src/geos-3.7.0.tar.bz2 && cd geos-3.7.0
./configure --prefix=$PWD/install
export GEOS_DIR=$PWD/install
make -j $NPROCS install
ln -s $PWD/install/lib $PWD/install/lib64
cd ../
# Compile and install basemap
tar xfz extra_src/basemap-1.1.0.tar.gz && cd basemap-1.1.0
python setup.py build
cp build/lib.linux-x86_64-2.7/_geoslib.so $INSTALL_PREFIX/lib/python2.7/site-packages/
cp -r build/lib.linux-x86_64-2.7/mpl_toolkits/basemap $INSTALL_PREFIX/lib/python2.7/site-packages/mpl_toolkits/
cd ../

# Compile and install proj
tar xzf extra_src/proj-6.1.0.tar.gz && cd proj-6.1.0
./configure --prefix=$INSTALL_PREFIX
unzip extra_src/proj-datumgrid-1.8.zip -d data/
make -j $NPROCS install
cd ../
# Compile and install pyproj
tar xzf extra_src/pyproj-1.9.6.tar.gz && cd pyproj-1.9.6
python setup.py build
cp -r build/lib.linux-x86_64-2.7/pyproj $INSTALL_PREFIX/lib/python2.7/site-packages/
cd ../

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
g++ -shared -Wl,-soname,libOGS -o libOGS.so -fPIC OGS.cpp netcdfio.cpp -I../../../../paraview-build/install/include/paraview-5.6
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
