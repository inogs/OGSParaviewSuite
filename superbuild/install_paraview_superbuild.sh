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
# Arnau Miro, OGS 2019

# Paraview version
PV_VERS=$1
INSTALL_PREFIX=../paraview-$PV_VERS
SUITEDIR='OGSParaviewSuite'
NPROCS=28

# Load modules
# provide a basic building environment
module purge
#module load cuda cmake
module load gcc openmpi cuda cmake

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
cd ../ && mkdir paraview-build && cd paraview-build
cmake ../paraview-superbuild/ \
   -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX \
   -DBUILD_TESTING=OFF \
   -Dparaview_SOURCE_SELECTION=$PV_VERS \
   -Dqt5_SOURCE_SELECTION=5.12 \
   -DCMAKE_BUILD_TYPE_paraview=Release \
   -DENABLE_mpi=ON \
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
   -DENABLE_paraviewsdk=OFF \
   -DENABLE_netcdf=OFF \
   -DENABLE_vrpn=ON \
   -DENABLE_vortexfinder2=OFF \
   -Dparaview_PLUGINS_EXTERNAL="OGSAnnotateDateTime;OGSDepthProfile;OGSDerivatives;OGSGradient;OGSHovmoeller;OGSLambda2Criterion;OGSOkuboWeiss;OGSOmegaCriterion;OGSPlotViews;OGSQCriterion;OGSReader;OGSSelectTools;OGSSpaghetti;OGSSpatialStats;OGSSpatialStatsFromFile;OGSTimeStatistics;OGSVariableAggregator" \
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
   -Dparaview_PLUGIN_OGSVariableAggregator_PATH=../paraview-superbuild/$SUITEDIR/OGSPlugins/OGSVariableAggregator 

# Prompt user to check configuration
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

cd ..
